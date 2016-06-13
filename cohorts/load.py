# Copyright (c) 2016. Mount Sinai School of Medicine
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.

from __future__ import print_function

from os import path, makedirs
from shutil import rmtree
import pandas as pd
import seaborn as sb
import json
import warnings

# pylint doesn't like this line
# pylint: disable=no-name-in-module
import six.moves.cPickle as pickle
from types import FunctionType

import vap  ## vcf-annotate-polyphen
from sqlalchemy import create_engine

import varcode
from varcode import VariantCollection, EffectCollection
from mhctools import NetMHCcons, EpitopeCollection
from topiary import predict_epitopes_from_variants, epitopes_to_dataframe
from topiary.sequence_helpers import contains_mutant_residues
from isovar.protein_sequence import variants_to_protein_sequences_dataframe
from pysam import AlignmentFile

from .survival import plot_kmf
from .plot import mann_whitney_plot, fishers_exact_plot, roc_curve_plot
from .collection import Collection
from .varcode_utils import filter_variants_with_metadata, filter_effects_with_metadata

class InvalidDataError(ValueError):
    pass

def require_id_str(id):
    if type(id) != str:
        raise ValueError("Expected ID string, but id = %s" % str(id))

class Sample(object):
    """
    Represents a single tumor or normal sample. It can point to DNA and/or
    RNA reads.

    Parameters
    __________
    is_tumor : bool
        Does this `Sample` represent a tumor sample?
    bam_path_dna : str
        Path to the DNA BAM file.
    bam_path_rna : str
        Path to the RNA BAM file.
    cufflinks_path : str
        Path to the Cufflinks output file.
    """
    def __init__(self,
                 is_tumor,
                 bam_path_dna=None,
                 bam_path_rna=None,
                 cufflinks_path=None):
        self.is_tumor = is_tumor
        self.bam_path_dna = bam_path_dna
        self.bam_path_rna = bam_path_rna
        self.cufflinks_path = cufflinks_path

class Patient(object):
    """
    Represents a patient, which contains zero or more of: normal `Sample`,
    tumor `Sample`.

    Parameters
    __________
    id : str
        ID of the patient.
    os : int
        Overall survival in days.
    pfs : int
        Progression-free survival in days.
    deceased : bool
        Is the patient deceased?
    progressed : bool
        Has the patient progressed?
    progressed_or_deceased : bool
        Has the patient either progressed or passed away?
    benefit : bool
        Has the patient seen a durable clinical benefit?
    snv_vcf_paths : list
        List of paths to SNV VCFs or this patient; multple VCFs get merged.
    indel_vcf_paths : list
        List of paths to indel VCFs for this patient; multple VCFs get merged.
    normal_sample : Sample
        This patient's normal `Sample`.
    tumor_sample: Sample
        This patient's tumor `Sample`.
    hla_alleles : list
        A list of this patient's HLA class I alleles.
    additional_data : dict
        A dictionary of additional data: name of datum mapping to value.
    """
    def __init__(self,
                 id,
                 os,
                 pfs,
                 deceased,
                 progressed=None,
                 progressed_or_deceased=None,
                 benefit=None,
                 snv_vcf_paths=[],
                 indel_vcf_paths=[],
                 normal_sample=None,
                 tumor_sample=None,
                 hla_alleles=None,
                 additional_data=None):
        require_id_str(id)
        self.id = id
        self.os = os
        self.pfs = pfs
        self.deceased = deceased
        self.progressed = progressed
        self.progressed_or_deceased = progressed_or_deceased
        self.benefit = benefit
        self.snv_vcf_paths = snv_vcf_paths
        self.indel_vcf_paths = indel_vcf_paths
        self.normal_sample = normal_sample
        self.tumor_sample = tumor_sample
        self.hla_alleles = hla_alleles
        self.additional_data = additional_data

        self.add_progressed_or_deceased()

    def add_progressed_or_deceased(self):
        assert self.progressed is not None or self.progressed_or_deceased is not None, (
            "Need at least one of progressed and progressed_or_deceased")
        if self.progressed_or_deceased is None:
            self.progressed_or_deceased = self.progressed or self.deceased

        # If we have both of these, ensure that they're in sync
        if self.progressed_or_deceased is not None and self.progressed is not None:
            assert self.progressed_or_deceased == self.progressed or self.deceased, (
                    "progressed_or_deceased should equal progressed || deceased")

    def __hash__(self):
        return hash(self.id)

    def __eq__(self, other):
        if self is other:
            return True
        return (
            self.__class__ == other.__class__ and
            self.id == other.id)

class DataFrameLoader(object):
    """
    Wraps a `DataFrame` with some information on how to join it.

    Parameters
    __________
    name : str
        The name of the dataframe, to easily reference it.
    load_dataframe : function
        A function that returns the `DataFrame` object.
    join_on : str
        The column of the `DataFrame` to join on (i.e. the patient
        ID column name).
    """
    def __init__(self,
                 name,
                 load_dataframe,
                 join_on="patient_id"):
        self.name = name
        self.load_dataframe = load_dataframe
        self.join_on = join_on

class Cohort(Collection):
    """
    Represents a cohort of `Patient`s.

    Parameters
    __________
    patients : List
        A list of `Patient`s for this cohort.
    cache_dir : str
        Path to store cached results, e.g. cached variant effects.
    cache_results : bool
        Whether or not to cache results.
    extra_df_loaders : List
        List of `DataFrameLoader`s to include as join options `Cohort`.
    join_with : str or List
        The name of one or more `DataFrameLoader`s to join with by default.
    join_how : str
        What type of default join to use for joining `DataFrameLoader`s.
    """
    def __init__(self,
                 patients,
                 cache_dir,
                 cache_results=True,
                 extra_df_loaders=[],
                 join_with=None,
                 join_how="inner",
                 check_provenance=False,
                 polyphen_dump_path=None):
        Collection.__init__(
            self,
            elements=patients)
        self.cache_dir = cache_dir
        self.cache_results = cache_results

        df_loaders = [
            DataFrameLoader("cufflinks", self.load_cufflinks)]
        df_loaders.extend(extra_df_loaders)
        self.df_loaders = df_loaders
        self.join_with = join_with
        self.join_how = join_how
        self.check_provenance = check_provenance
        self.polyphen_dump_path = polyphen_dump_path

        self.verify_id_uniqueness()
        self.verify_survival()

        self.cache_names = {"variant": "cached-variants",
                            "effect": "cached-effects",
                            "nonsynonymous_effect": "cached-nonsynonymous-effects",
                            "neoantigen": "cached-neoantigens",
                            "expressed_neoantigen": "cached-expressed-neoantigens",
                            "polyphen": "cached-polyphen-annotations",
                            "isovar": "cached-isovar-output"}

    def verify_id_uniqueness(self):
        patient_ids = set([patient.id for patient in self])
        if len(patient_ids) != len(self):
            raise ValueError("Non-unique patient IDs")

    def verify_survival(self):
        cohort_dataframe = self.as_dataframe()
        if not (cohort_dataframe["pfs"] <=
                cohort_dataframe["os"]).all():
            raise InvalidDataError("PFS should be <= OS, but PFS is larger than OS for some patients.")

        def func(row):
            if row["pfs"] < row["os"]:
                if not row["progressed_or_deceased"]:
                    raise InvalidDataError(
                        "A patient did not progress despite PFS being less than OS. "
                        "Full row: %s" % row)
        cohort_dataframe.apply(func, axis=1)

    def _as_dataframe_unmodified(self, join_with, join_how):
        # Use join_with if specified, otherwise fall back to what is defined in the class
        join_with = first_not_none_param([join_with, self.join_with], default=[])
        if type(join_with) == str:
            join_with = [join_with]

        # Convert strings to DataFrameLoader objects
        df_loaders = [df_loader for df_loader in self.df_loaders if df_loader.name in join_with]

        # Use join_how if specified, otherwise fall back to what is defined in the class
        join_how = first_not_none_param([join_how, self.join_how], default="inner")

        patient_rows = []
        for patient in self:
            row = {} if patient.additional_data is None else patient.additional_data.copy()
            row['patient_id'] = patient.id
            for clinical_col in ["benefit", "os", "pfs", "deceased",
                                "progressed", "progressed_or_deceased"]:
                row[clinical_col] = getattr(patient, clinical_col)

            patient_rows.append(row)
        df = pd.DataFrame.from_records(patient_rows)

        for df_loader in df_loaders:
            old_len_df = len(df)
            df = df.merge(
                df_loader.load_dataframe(),
                left_on="patient_id",
                right_on=df_loader.join_on,
                how=join_how)
            print("%s join with %s: %d to %d rows" % (
                join_how,
                df_loader.name,
                old_len_df,
                len(df)))
        return df

    def as_dataframe(self, on=None, col=None, join_with=None, join_how=None, **kwargs):
        """
        Return this Cohort as a DataFrame, and optionally include additional columns
        using `on`.

        on : function or list or dict or str, optional
            - A function that creates a new column for comparison, e.g. count.snv_count.
            - Or a list of column-generating functions.
            - Or a map of new column names to their column-generating functions.
            - A column name that gets returned with the original dataframe as (col, df).
        col : str, optional
            If `on` is a function generating a column, col is the name of that column.
            If `on` is a str specifying a column, col is the name of a copy of that column.
            If None, defaults to the name of the function.

        If `on` is a function or functions, kwargs is passed to those functions.
        Otherwise kwargs is ignored.

        Return : tuple or DataFrame
            <dataframe>
            or (<name of new column>, <dataframe with that new column>)
            or (<list of names of new columns>, <dataframe with those new columns>)
        """
        df = self._as_dataframe_unmodified(join_with=join_with, join_how=join_how)
        if on is None:
            return df

        if type(on) == str:
            if col is not None:
                df[col] = df[on]
                return (col, df)
            return (on, df)

        def apply_func(on, col, df):
            """
            Sometimes we have functions that, by necessity, have more parameters
            than just `row`. For a function like that, we assume that it has two
            parameters: `row` and `cohort`. We use it to construct another function
            that has only the `row` parameter so it can be sent to `DataFrame.apply`.
            """
            if col is None:
                # Use the function name, or "column" for lambdas, if no
                # name is provided for the newly created column.
                col = on.__name__ if not is_lambda(on) else "column"
            on_argcount = on.__code__.co_argcount
            if on_argcount == 1:
                func = lambda row: on(row, **kwargs)
            elif on_argcount == 2:
                func = lambda row: on(row, self, **kwargs)
            else:
                raise ValueError("Function %s has too many arguments" % on)
            df[col] = df.apply(func, axis=1)
            return (col, df)

        def is_lambda(func):
            return func.__name__ == (lambda: None).__name__

        if type(on) == FunctionType:
            return apply_func(on, col, df)
        if type(on) == dict:
            cols = []
            for key, value in on.iteritems():
                col, df = apply_func(on=value, col=key, df=df)
                cols.append(col)
            return (cols, df)
        if type(on) == list:
            cols = []
            for i, elem in enumerate(on):
                col = elem.__name__ if not is_lambda(elem) else "column_%d" % i
                col, df = apply_func(on=elem, col=col, df=df)
                cols.append(col)
            return (cols, df)

    def load_dataframe(self, df_loader_name):
        """
        Instead of joining a DataFrameJoiner with the Cohort in `as_dataframe`, sometimes
        we may want to just directly load a particular DataFrame.
        """
        # Get the DataFrameLoader object corresponding to this name.
        df_loaders = [df_loader for df_loader in self.df_loaders if df_loader.name == df_loader_name]

        if len(df_loaders) == 0:
            raise ValueError("No DataFrameLoader with name %s" % df_loader_name)
        if len(df_loaders) > 1:
            raise ValueError("Multiple DataFrameLoaders with name %s" % df_loader_name)

        return df_loaders[0].load_dataframe()

    def generate_provenance(self):
        module_names = ["cohorts", "pyensembl", "varcode", "mhctools", "topiary", "isovar", "scipy", "numpy", "pandas"]
        module_versions = [__import__(module_name).__version__ for module_name in module_names]
        return dict(zip(module_names, module_versions))

    def load_provenance(self, patient_cache_dir):
        with open(path.join(patient_cache_dir, "PROVENANCE"), "r") as f:
            return json.load(f)

    def save_provenance(self, patient_cache_dir, provenance):
        with open(path.join(patient_cache_dir, "PROVENANCE"), "w") as f:
            json.dump(provenance, f)

    def load_from_cache(self, cache_name, patient_id, file_name):
        if not self.cache_results:
            return None

        cache_dir = path.join(self.cache_dir, cache_name)
        patient_cache_dir = path.join(cache_dir, str(patient_id))
        cache_file = path.join(patient_cache_dir, file_name)

        if not path.exists(cache_file):
            return None

        if self.check_provenance:
            provenance = set(self.generate_provenance().items())
            provenance_previous = set(self.load_provenance(patient_cache_dir).items())

            def provenance_str(provenance):
                return ["%s==%s" % (key, value) for (key, value) in provenance]

            # Two-way diff: are any modules introduced, and are any modules lost?
            new_diff = provenance.difference(provenance_previous)
            old_diff = provenance_previous.difference(provenance)
            warn_str = ""
            if len(new_diff) > 0:
                warn_str += "In current environment but not cached for patient %s: %s" % (
                    patient_id, provenance_str(new_diff))
            if len(old_diff) > 0:
                warn_str += "In cached environment for patient %s but not current: %s" % (
                    patient_id, provenance_str(old_diff))

            if len(warn_str) > 0:
                warnings.warn(warn_str, Warning)

        if path.splitext(cache_file)[1] == ".csv":
            return pd.read_csv(cache_file, dtype={"patient_id": object})
        else:
            with open(cache_file, "rb") as f:
                return pickle.load(f)

    def save_to_cache(self, obj, cache_name, patient_id, file_name):
        if not self.cache_results:
            return

        cache_dir = path.join(self.cache_dir, cache_name)
        patient_cache_dir = path.join(cache_dir, str(patient_id))
        cache_file = path.join(patient_cache_dir, file_name)

        if not path.exists(patient_cache_dir):
            makedirs(patient_cache_dir)

        if type(obj) == pd.DataFrame:
            obj.to_csv(cache_file, index=False)
        else:
            with open(cache_file, "wb") as f:
                pickle.dump(obj, f)

        provenance = self.generate_provenance()
        self.save_provenance(patient_cache_dir, provenance)

    def iter_patients(self, patients):
        if patients is None:
            return self
        return iter(patients)

    def patient_from_id(self, id):
        for patient in self:
            if patient.id == id:
                return patient
        raise ValueError("No patient with ID %s found" % id)

    def load_variants(self, patients=None, variant_type="snv", merge_type="union", filter_fn=None):
        """Load a dictionary of patient_id to varcode.VariantCollection

        Parameters
        ----------
        patients : str, optional
            Filter to a subset of patients
        variant_type : {'snv', 'indel'}, optional
            Load variants of a specific type, default 'snv'
        merge_type : {'union', 'intersection'}, optional
            Use this method to merge multiple variant sets for a single patient, default 'union'
        filter_fn: function
            Takes a variant and it's metadata and returns a boolean. Only variants returning True are preserved.

        Returns
        -------
        merged_variants
            Dictionary of patient_id to VariantCollection
        """
        assert variant_type in ["snv", "indel"], "Unknown variant type: %s" % variant_type
        patient_variants = {}

        for patient in self.iter_patients(patients):
            variants = self._load_single_patient_variants(patient, variant_type, merge_type, filter_fn)
            patient_variants[patient.id] = variants
        return patient_variants

    def _load_single_patient_variants(self, patient, variant_type, merge_type, filter_fn=None):
        failed_io = False
        try:
            cached_file_name = "%s-%s-variants.pkl" % (variant_type, merge_type)
            cached = self.load_from_cache(self.cache_names["variant"], patient.id, cached_file_name)
            if cached is not None:
                return filter_variants_with_metadata(cached, filter_fn)
            vcf_paths = patient.snv_vcf_paths if variant_type == "snv" else patient.indel_vcf_paths
            variant_collections = [
                (vcf_path, varcode.load_vcf_fast(vcf_path)) 
                for vcf_path in vcf_paths
            ]
        except IOError: 
            failed_io = True

        if failed_io or len(variant_collections) == 0:
            print("Variants did not exist for patient %s" % patient.id)
            return VariantCollection([])

        if len(variant_collections) == 1:
            # There is nothing to merge
            vcf_path, variants = variant_collections[0]
            merged_variants = variants
        else:
            merged_variants = self._merge_variant_collections(dict(variant_collections), merge_type)

        self.save_to_cache(merged_variants, self.cache_names["variant"], patient.id, cached_file_name)

        return filter_variants_with_metadata(merged_variants, filter_fn)

    def _merge_variant_collections(self, vcf_to_variant_collections, merge_type):
        assert merge_type in ["union", "intersection"], "Unknown merge type: %s" % merge_type
        combined_variants = [set(vc.elements) for (path, vc) in vcf_to_variant_collections.items()]
        if merge_type == "union":
            merged_variants = VariantCollection(set.union(*combined_variants))
        elif merge_type == "intersection":
            merged_variants = VariantCollection(set.intersection(*combined_variants))

        for variant in merged_variants:
            metadata = [
                (vcf_path, vc.metadata[variant]) 
                for (vcf_path, vc) in vcf_to_variant_collections.items() 
                if variant in vc.metadata
            ]

            merged_variants.metadata[variant] = dict(metadata)

        return merged_variants

    def load_polyphen_annotations(self, as_dataframe=False):
        """Load a dataframe containing polyphen2 annotations for all variants

        Parameters
        ----------
        database_file : string, sqlite
            Path to the WHESS/Polyphen2 SQLite database.
            Can be downloaded and bunzip2'ed from http://bit.ly/208mlIU

        Returns
        -------
        annotations
            Dictionary of patient_id to a DataFrame that contains annotations
        """
        patient_annotations = {}
        for patient in self:
            annotations = self._load_single_patient_polyphen(patient)
            annotations["patient_id"] = patient.id
            patient_annotations[patient.id] = annotations
        if as_dataframe:
            return pd.concat(patient_annotations.values())
        return patient_annotations

    def _load_single_patient_polyphen(self, patient):
        cache_name = self.cache_names["polyphen"]
        cached_file_name = "polyphen-annotations.csv"
        cached = self.load_from_cache(cache_name, patient.id, cached_file_name)
        if cached is not None:
            return cached

        engine = create_engine("sqlite:///{}".format(self.polyphen_dump_path))
        conn = engine.connect()

        variants = self._load_single_patient_variants(patient,
                                                      variant_type="snv",
                                                      merge_type="union")
        df = pd.DataFrame(columns=["chrom", "pos", "ref", "alt",
                                   "annotation_found", "gene", "protein",
                                   "aa_change", "hvar_pred", "hvar_prob",
                                   "hdiv_pred", "hdiv_prob"])
        for variant in variants:
            chrom = "chr{}".format(getattr(variant, "contig", None))
            pos = getattr(variant, "start", None)
            ref = getattr(variant, "ref", None)
            alt = getattr(variant, "alt",  None)
            annotation = vap.annotate_variant(conn, chrom, pos, ref, alt)
            datum = {"chrom": chrom,
                     "pos": pos,
                     "ref": ref,
                     "alt": alt,
                     "annotation_found": annotation is not None}
            attributes = ["gene", "protein", "aa_change",
                          "hvar_pred", "hvar_prob",
                          "hdiv_pred", "hdiv_prob"]
            for attr in attributes:
                datum[attr] = getattr(annotation, attr, None)
            df = df.append(datum, ignore_index=True)
        df["pos"] = df["pos"].astype("int")
        df["annotation_found"] = df["annotation_found"].astype("bool")
        self.save_to_cache(df, cache_name, patient.id, cached_file_name)
        return df

    def load_effects(self, patients=None, only_nonsynonymous=False, variant_type="snv", merge_type="union", filter_fn=None):
        """Load a dictionary of patient_id to varcode.EffectCollection
        Parameters
        ----------
        patients : str, optional
            Filter to a subset of patients
        only_nonsynonymous : bool, optional
            If true, load only nonsynonymous effects, default False
        variant_type : {'snv', 'indel'}, optional
            Load variants of a specific type, default 'snv'
        merge_type : {'union', 'intersection'}, optional
            Use this method to merge multiple variant sets for a single patient, default 'union'
        filter_fn: function
            Takes an effect and it's variant's metadata and returns a boolean. Only effects returning True are preserved.

        Returns
        -------
        effects
             Dictionary of patient_id to varcode.EffectCollection
        """
        patient_effects = {}
        for patient in self.iter_patients(patients):
            effects = self._load_single_patient_effects(
                patient, only_nonsynonymous, variant_type, merge_type, filter_fn)
            patient_effects[patient.id] = effects
        return patient_effects

    def _load_single_patient_effects(self, patient, only_nonsynonymous,
                                     variant_type, merge_type, filter_fn):
        cached_file_name = "%s-%s-effects.pkl" % (variant_type, merge_type)
        variants = self._load_single_patient_variants(
            patient, variant_type, merge_type)
        if only_nonsynonymous:
            cached = self.load_from_cache(self.cache_names["nonsynonymous_effect"], patient.id, cached_file_name)
        else:
            cached = self.load_from_cache(self.cache_names["effect"], patient.id, cached_file_name)
        if cached is not None:
            return filter_effects_with_metadata(cached, variants.metadata, filter_fn)

        effects = variants.effects()
        nonsynonymous_effects = EffectCollection(
            effects.drop_silent_and_noncoding().top_priority_effect_per_variant().values())

        self.save_to_cache(effects, self.cache_names["effect"], patient.id, cached_file_name)
        self.save_to_cache(nonsynonymous_effects, self.cache_names["nonsynonymous_effect"], patient.id, cached_file_name)

        if only_nonsynonymous:
            return filter_effects_with_metadata(nonsynonymous_effects, variants.metadata, filter_fn)
        return filter_effects_with_metadata(effects, variants.metadata, filter_fn)

    def load_cufflinks(self, filter_ok=True):
        """
        Load a Cufflinks gene expression data for a cohort

        Parameters
        ----------
        filter_ok : bool, optional
            If true, filter Cufflinks data to row with FPKM_status == 'OK'

        Returns
        -------
        cufflinks_data : Pandas dataframe
            Pandas dataframe with Cufflinks data for all patients
            columns include patient_id, gene_id, gene_short_name, FPKM, FPKM_conf_lo, FPKM_conf_hi
        """
        return \
            pd.concat(
                [self._load_single_patient_cufflinks(patient, filter_ok) for patient in self],
                copy=False
        )

    def _load_single_patient_cufflinks(self, patient, filter_ok):
        """
        Load Cufflinks gene quantification given a patient

        Parameters
        ----------
        patient : Patient
        filter_ok : bool, optional
            If true, filter Cufflinks data to row with FPKM_status == 'OK'

        Returns
        -------
        cufflinks_data: Pandas dataframe
            Pandas dataframe of sample's Cufflinks data
            columns include patient_id, gene_id, gene_short_name, FPKM, FPKM_conf_lo, FPKM_conf_hi
        """
        data = pd.read_csv(patient.tumor_sample.cufflinks_path, sep="\t")
        data["patient_id"] = patient.id

        if filter_ok:
            # Filter to OK FPKM counts
            data = data[data["FPKM_status"] == "OK"]
        return data

    def load_neoantigens(self, patients=None, variant_type="snv", merge_type="union",
                         only_expressed=False, epitope_lengths=[8, 9, 10, 11],
                         ic50_cutoff=500, process_limit=10, max_file_records=None):
        dfs = []
        for patient in self.iter_patients(patients):
            df_epitopes = self._load_single_patient_neoantigens(
                patient=patient,
                variant_type=variant_type,
                merge_type=merge_type,
                only_expressed=only_expressed,
                epitope_lengths=epitope_lengths,
                ic50_cutoff=ic50_cutoff,
                process_limit=process_limit,
                max_file_records=max_file_records)
            dfs.append(df_epitopes)
        return pd.concat(dfs)

    def _load_single_patient_neoantigens(self, patient, variant_type,
                                         merge_type, only_expressed, epitope_lengths,
                                         ic50_cutoff, process_limit, max_file_records):
        cached_file_name = "%s-%s-neoantigens.csv" % (variant_type, merge_type)
        if only_expressed:
            cached = self.load_from_cache(self.cache_names["expressed_neoantigen"], patient.id, cached_file_name)
        else:
            cached = self.load_from_cache(self.cache_names["neoantigen"], patient.id, cached_file_name)
        if cached is not None:
            return cached

        variants = self._load_single_patient_variants(
            patient, variant_type, merge_type)
        mhc_model = NetMHCcons(
            alleles=patient.hla_alleles,
            epitope_lengths=epitope_lengths,
            max_file_records=max_file_records,
            process_limit=process_limit)
        if only_expressed:
            isovar_cached_file_name = "%s-%s-isovar.csv" % (variant_type, merge_type)
            df_isovar = self.load_from_cache(self.cache_names["isovar"], patient.id, isovar_cached_file_name)
            if df_isovar is None:
                import logging
                logging.disable(logging.INFO)
                if patient.tumor_sample is None:
                    raise ValueError("Patient %s has no tumor sample" % patient.id)
                if patient.tumor_sample.bam_path_rna is None:
                    raise ValueError("Patient %s has no tumor RNA BAM path" % patient.id)
                rna_bam_file = AlignmentFile(patient.tumor_sample.bam_path_rna)
                from isovar.default_parameters import (
                    MIN_TRANSCRIPT_PREFIX_LENGTH,
                    MAX_REFERENCE_TRANSCRIPT_MISMATCHES
                )

                # To ensure that e.g. 8-11mers overlap substitutions, we need at least this
                # sequence length: (max peptide length * 2) - 1
                # Example:
                # 123456789AB
                #           123456789AB
                # AAAAAAAAAAVAAAAAAAAAA
                protein_sequence_length = (max(epitope_lengths) * 2) - 1
                df_isovar = variants_to_protein_sequences_dataframe(
                    variants=variants,
                    samfile=rna_bam_file,
                    protein_sequence_length=protein_sequence_length,
                    min_reads_supporting_rna_sequence=3, # Per Alex R.'s suggestion
                    min_transcript_prefix_length=MIN_TRANSCRIPT_PREFIX_LENGTH,
                    max_transcript_mismatches=MAX_REFERENCE_TRANSCRIPT_MISMATCHES,
                    max_protein_sequences_per_variant=1, # Otherwise we might have too much neoepitope diversity
                    min_mapping_quality=0)
                self.save_to_cache(df_isovar, self.cache_names["isovar"], patient.id, isovar_cached_file_name)

            # Map from isovar rows to protein sequences
            isovar_rows_to_protein_sequences = dict([
                (frozenset(row.to_dict().items()), row["amino_acids"]) for (i, row) in df_isovar.iterrows()])

            # MHC binding prediction
            epitopes = mhc_model.predict(isovar_rows_to_protein_sequences)

            # Call `get_filtered_isovar_epitopes` in order to only include peptides that
            # overlap a variant; without this filter, when we use
            # protein_sequence_length above, some 8mers generated from a 21mer source will
            # not overlap a variant.
            df_epitopes = self.get_filtered_isovar_epitopes(
                epitopes, ic50_cutoff=ic50_cutoff).dataframe()
            df_epitopes["patient_id"] = patient.id

            self.save_to_cache(df_epitopes, self.cache_names["expressed_neoantigen"], patient.id, cached_file_name)
        else:
            epitopes = predict_epitopes_from_variants(
                variants=variants,
                mhc_model=mhc_model,
                ic50_cutoff=ic50_cutoff,
                # Only include peptides with a variant
                only_novel_epitopes=True)
            df_epitopes = epitopes_to_dataframe(epitopes)
            df_epitopes["patient_id"] = patient.id

            self.save_to_cache(df_epitopes, self.cache_names["neoantigen"], patient.id, cached_file_name)

        return df_epitopes

    def get_filtered_isovar_epitopes(self, epitopes, ic50_cutoff):
        """
        Mostly replicates topiary.build_epitope_collection_from_binding_predictions

        Note: topiary needs to do fancy stuff like subsequence_protein_offset + binding_prediction.offset
        in order to figure out whether a variant is in the peptide because it only has the variant's
        offset into the full protein; but isovar gives us the variant's offset into the protein subsequence
        (dictated by protein_sequence_length); so all we need to do is map that onto the smaller 8-11mer
        peptides generated by mhctools.
        """
        mutant_binding_predictions = []
        for binding_prediction in epitopes:
            peptide = binding_prediction.peptide
            peptide_offset = binding_prediction.offset
            isovar_row = dict(binding_prediction.source_sequence_key)
            is_mutant = contains_mutant_residues(
                peptide_start_in_protein=peptide_offset,
                peptide_length=len(peptide),
                mutation_start_in_protein=isovar_row["variant_aa_interval_start"],
                mutation_end_in_protein=isovar_row["variant_aa_interval_end"])
            if is_mutant and binding_prediction.value <= ic50_cutoff:
                mutant_binding_predictions.append(binding_prediction)
        return EpitopeCollection(mutant_binding_predictions)

    def clear_caches(self):
        for cache in self.cache_names.keys():
            self.clear_cache(cache)

    def clear_cache(self, cache):
        cache_path = path.join(self.cache_dir, self.cache_names[cache])
        if path.exists(cache_path):
            rmtree(cache_path)

    def cohort_columns(self):
        cohort_dataframe = self.as_dataframe()
        column_types = [cohort_dataframe[col].dtype for col in cohort_dataframe.columns]
        return dict(zip(list(cohort_dataframe.columns), column_types))

    def plot_roc_curve(self, on, bootstrap_samples=100, **kwargs):
        """Plot an ROC curve for benefit and a given variable

        Parameters
        ----------
        on : str or function
            See `cohort.load.as_dataframe`
        bootstrap_samples : int, optional
            Number of boostrap samples to use to compute the AUC

        Returns
        -------
        (mean_auc_score, plot): (float, matplotlib plot)
            Returns the average AUC for the given predictor over `bootstrap_samples`
            and the associated ROC curve
        """
        plot_col, df = self.as_dataframe(on, col=None, **kwargs)
        df = filter_not_null(df, "benefit")
        df = filter_not_null(df, plot_col)
        df.benefit = df.benefit.astype(bool)
        return roc_curve_plot(df, plot_col, 'benefit', bootstrap_samples)

    def plot_benefit(self, on, col=None, benefit_col="benefit",
                     mw_alternative="two-sided", **kwargs):
        """Plot a comparison of benefit/response in the cohort on a given variable
        """
        return self.plot_boolean(on=on,
                                 boolean_col=benefit_col,
                                 col=col,
                                 mw_alternative=mw_alternative,
                                 **kwargs)

    def plot_boolean(self, on, boolean_col, col=None,
                     mw_alternative="two-sided", **kwargs):
        """Plot a comparison of `boolean_col` in the cohort on a given variable via
        `on` or `col`.

        If the variable (through `on` or `col` is binary) this will compare
        odds-ratios and perform a Fisher's exact test.

        If the variable is numeric, this will compare the distributions through
        a Mann-Whitney test and plot the distributions with box-strip plot

        Parameters
        ----------
        on : str or function
            See `cohort.load.as_dataframe`
        col : str, optional
            If specified, store the result of `on`. See `cohort.load.as_dataframe`
        mw_alternative : str, optional
            Choose the sidedness of the mannwhitneyu test.

        Returns
        -------
        (Test statistic, p-value): (float, float)

        """
        plot_col, df = self.as_dataframe(on, col, **kwargs)
        df = filter_not_null(df, boolean_col)
        df = filter_not_null(df, plot_col)
        df[boolean_col] = df[boolean_col].astype(bool)
        if df[plot_col].dtype == "bool":
            results = fishers_exact_plot(
                data=df,
                condition1=boolean_col,
                condition2=plot_col)
        else:
            results = mann_whitney_plot(
                data=df,
                condition=boolean_col,
                distribution=plot_col,
                alternative=mw_alternative)
        return results

    def plot_survival(self, on, col=None, how="os",
                      threshold=None, **kwargs):
        """Plot a Kaplan Meier survival curve by splitting the cohort into two groups
        Parameters
        ----------
        on : str or function
            See `cohort.load.as_dataframe`
        col : str, optional
            If specified, store the result of `on`. See `cohort.load.as_dataframe`
        how : {'os', 'pfs'}, optional
            Whether to plot OS (overall survival) or PFS (progression free survival)
        threshold : int or 'median', optional
            Threshold of `col` on which to split the cohort
        """
        assert how in ["os", "pfs"], "Invalid choice of survival plot type %s" % how
        plot_col, df = self.as_dataframe(on, col, **kwargs)
        if df[plot_col].dtype == "bool":
            default_threshold = None
        else:
            default_threshold = "median"
        results = plot_kmf(
            df=df,
            condition_col=plot_col,
            xlabel="Overall Survival" if how == "os" else "Progression-Free Survival",
            censor_col="deceased" if how == "os" else "progressed_or_deceased",
            survival_col=how,
            threshold=threshold if threshold is not None else default_threshold)
        print(results)

    def plot_joint(self, on, on_two=None, **kwargs):
        """Plot a jointplot.

        Parameters
        ----------
        on : function or list or map of functions
            See `cohort.load.as_dataframe`
        on_two : function, optional
            Can specify the second function here rather than creating a list.
        """
        if on_two is not None:
            on = [on, on_two]
        plot_cols, df = self.as_dataframe(on, **kwargs)
        sb.jointplot(data=df, x=plot_cols[0], y=plot_cols[1])

def first_not_none_param(params, default):
    """
    Given a list of `params`, use the first param in the list that is
    not None. If all are None, fall back to `default`.
    """
    for param in params:
        if param is not None:
            return param
    return default

def filter_not_null(df, col):
    original_len = len(df)
    df = df[df[col].notnull()]
    updated_len = len(df)
    if updated_len < original_len:
        print("Missing %s for %d patients: from %d to %d" % (col, original_len - updated_len, original_len, updated_len))
    return df
