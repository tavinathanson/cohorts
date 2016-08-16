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
import pprint

# pylint doesn't like this line
# pylint: disable=no-name-in-module
import six.moves.cPickle as pickle
from types import FunctionType

import vap  ## vcf-annotate-polyphen
from sqlalchemy import create_engine

from pyensembl import cached_release

import varcode
from varcode import VariantCollection, EffectCollection, Variant
from mhctools import NetMHCcons, EpitopeCollection
from topiary import predict_epitopes_from_variants, epitopes_to_dataframe
from topiary.sequence_helpers import contains_mutant_residues
from isovar.protein_sequence import variants_to_protein_sequences_dataframe
from pysam import AlignmentFile

from .utils import strip_column_names as _strip_column_names
from .survival import plot_kmf
from .plot import mann_whitney_plot, fishers_exact_plot, roc_curve_plot
from .collection import Collection
from .varcode_utils import (filter_variants, filter_effects,
                            filter_neoantigens, filter_polyphen)
from .variant_filters import no_filter
from .styling import set_styling
from . import variant_filters

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
                 cufflinks_path=None,
                 kallisto_path=None):
        self.is_tumor = is_tumor
        self.bam_path_dna = bam_path_dna
        self.bam_path_rna = bam_path_rna
        self.cufflinks_path = cufflinks_path
        self.kallisto_path = kallisto_path

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
                 additional_data=None,
                 cohort=None):
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

        # TODO: This can be removed once all patient-specific functions are
        # removed from Cohort.
        self.cohort = cohort

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
    kallisto_ensembl_version : int
        Cached release version to use from pyensembl to annotate Kallisto data
    cache_results : bool
        Whether or not to cache results.
    extra_df_loaders : List
        List of `DataFrameLoader`s to include as join options `Cohort`.
    join_with : str or List
        The name of one or more `DataFrameLoader`s to join with by default.
    join_how : str
        What type of default join to use for joining `DataFrameLoader`s.
    filter_fn : Function
        Specify a default filter function for `load_variants`, `load_effects`,
        `functions.missense_snv_count`, etc.
    normalized_per_mb : bool
        Whether or not to normalize by number of loci.
    min_coverage_depth_normal : int
        When counting number of exonic loci, only count loci with at least this much normal depth.
    min_coverage_depth_tumor : int
        When counting number of exonic loci, only count loci with at least this much tumor depth.
    responder_pfs_equals_os : bool
        Ensure that the PFS values for responders (not progressed) are equal to
        OS values.
    check_provenance : bool
        Verify that the cached provenance is equal to the current environment.
    print_provenance : bool
        Print a summary of cache file provenance.
    polyphen_dump_path : str
        Path to a Polyphen database dump.
    pageant_coverage_path : str
        Path to Pageant CoverageDepth output.
    variant_type : {"snv", "indel"}, optional
        Load variants of a specific type, default "snv"
    merge_type : {"union", "intersection"}, optional
        Use this method to merge multiple variant sets for a single patient, default "union"
    """
    def __init__(self,
                 patients,
                 cache_dir,
                 kallisto_ensembl_version=None,
                 cache_results=True,
                 extra_df_loaders=[],
                 join_with=None,
                 join_how="inner",
                 filter_fn=None,
                 normalized_per_mb=False,
                 min_coverage_normal_depth=0,
                 min_coverage_tumor_depth=0,
                 responder_pfs_equals_os=False,
                 check_provenance=False,
                 print_provenance=True,
                 polyphen_dump_path=None,
                 pageant_coverage_path=None,
                 variant_type="snv",
                 merge_type="union"):
        Collection.__init__(
            self,
            elements=patients)
        # TODO: Patients shouldn't actually need to reference their Cohort; remove
        # this when patient-specific functions all live in Patient.
        for patient in patients:
            patient.cohort = self
        self.cache_dir = cache_dir
        self.cache_results = cache_results
        self.kallisto_ensembl_version = kallisto_ensembl_version

        df_loaders = [
            DataFrameLoader("kallisto", self.load_kallisto),
            DataFrameLoader("cufflinks", self.load_cufflinks),
            DataFrameLoader("ensembl_coverage", self.load_ensembl_coverage)]
        df_loaders.extend(extra_df_loaders)
        self.df_loaders = df_loaders
        self.join_with = join_with
        self.join_how = join_how
        self.filter_fn = filter_fn
        self.normalized_per_mb = normalized_per_mb
        self.min_coverage_normal_depth = min_coverage_normal_depth
        self.min_coverage_tumor_depth = min_coverage_tumor_depth
        self.responder_pfs_equals_os = responder_pfs_equals_os
        self.check_provenance = check_provenance
        self.polyphen_dump_path = polyphen_dump_path
        self.pageant_coverage_path = pageant_coverage_path
        self.variant_type = variant_type
        self.merge_type = merge_type

        assert variant_type in ["snv", "indel"], "Unknown variant type: %s" % variant_type

        self.verify_id_uniqueness()
        self.verify_survival()
        self.dataframe_hash = None

        self.cache_names = {"variant": "cached-variants",
                            "effect": "cached-effects",
                            "nonsynonymous_effect": "cached-nonsynonymous-effects",
                            "neoantigen": "cached-neoantigens",
                            "expressed_neoantigen": "cached-expressed-neoantigens",
                            "polyphen": "cached-polyphen-annotations",
                            "isovar": "cached-isovar-output"}
        if print_provenance:
            pprint.pprint(self.summarize_data_sources())

        set_styling()

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
        if self.responder_pfs_equals_os:
            cohort_dataframe.apply(func, axis=1)

    def _as_dataframe_unmodified(self, join_with = None, join_how = None):
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
            row["patient_id"] = patient.id
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
        self.dataframe_hash = hash(str(df.sort_values("patient_id")))
        return df

    def as_dataframe(self, on=None, col=None, join_with=None, join_how=None,
                     rename_cols=False, keep_paren_contents=True, **kwargs):
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

        Other parameters
        ----------------
        `rename_cols`: (bool)
            if True, then return columns using "stripped" column names
            ("stripped" means lower-case names without punctuation other than `_`)
            See `utils.strip_column_names` for more details
            defaults to False
        `keep_paren_contents`: (bool)
            if True, then contents of column names within parens are kept.
            if False, contents of column names within-parens are dropped.
            Defaults to True
        ----------

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
            than just `row`. We construct a function with just the `row` parameter
            so it can be sent to `DataFrame.apply`. We hackishly pass `cohort`
            (as `self`) along if the function accepts a `cohort` argument.
            """
            if col is None:
                # Use the function name, or "column" for lambdas, if no
                # name is provided for the newly created column.
                col = on.__name__ if not is_lambda(on) else "column"
            on_argnames = on.__code__.co_varnames
            if "cohort" not in on_argnames:
                func = lambda row: on(row=row, **kwargs)
            else:
                func = lambda row: on(row=row, cohort=self, **kwargs)
            df[col] = df.apply(func, axis=1)
            return (col, df)

        def is_lambda(func):
            return func.__name__ == (lambda: None).__name__

        if type(on) == FunctionType:
            return apply_func(on, col, df)

        # For multiple functions, don't allow kwargs since we won't know which functions
        # they apply to.
        if len(kwargs) > 0:
            raise ValueError("kwargs are not supported when collecting multiple functions "
                             "as we don't know which function they apply to.")

        if type(on) == dict:
            cols = []
            for key, value in on.iteritems():
                col, df = apply_func(on=value, col=key, df=df)
                cols.append(col)
        if type(on) == list:
            cols = []
            for i, elem in enumerate(on):
                col = elem.__name__ if not is_lambda(elem) else "column_%d" % i
                col, df = apply_func(on=elem, col=col, df=df)
                cols.append(col)

        if rename_cols:
            rename_dict = _strip_column_names(df.columns, keep_paren_contents=keep_paren_contents)
            df.rename(columns=rename_dict, inplace=True)
            cols = [rename_dict[col] for col in cols]
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
            num_discrepant = compare_provenance(
                this_provenance = self.generate_provenance(), 
                other_provenance = self.load_provenance(patient_cache_dir),
                left_outer_diff = "In current environment but not cached in %s for patient %s" % (cache_name, patient_id),
                right_outer_diff = "In cached %s for patient %s but not current" % (cache_name, patient_id)
                )

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

    def load_variants(self, patients=None, filter_fn=None):
        """Load a dictionary of patient_id to varcode.VariantCollection

        Parameters
        ----------
        patients : str, optional
            Filter to a subset of patients
        filter_fn : function
            Takes a FilterableVariant and returns a boolean. Only variants returning True are preserved.
            Overrides default self.filter_fn. `None` passes through to self.filter_fn.

        Returns
        -------
        merged_variants
            Dictionary of patient_id to VariantCollection
        """
        filter_fn = first_not_none_param([filter_fn, self.filter_fn], no_filter)
        patient_variants = {}

        for patient in self.iter_patients(patients):
            variants = self._load_single_patient_variants(patient, filter_fn)
            if variants is not None:
                patient_variants[patient.id] = variants
        return patient_variants

    def _load_single_patient_variants(self, patient, filter_fn):
        failed_io = False
        try:
            cached_file_name = "%s-%s-variants.pkl" % (self.variant_type, self.merge_type)
            cached = self.load_from_cache(self.cache_names["variant"], patient.id, cached_file_name)
            if cached is not None:
                return filter_variants(variant_collection=cached,
                                       patient=patient,
                                       filter_fn=filter_fn)
            vcf_paths = patient.snv_vcf_paths if self.variant_type == "snv" else patient.indel_vcf_paths
            variant_collections = [
                (vcf_path, varcode.load_vcf_fast(vcf_path))
                for vcf_path in vcf_paths
            ]
        except IOError:
            failed_io = True

        # Note that this is the number of variant collections and not the number of
        # variants. 0 variants will lead to 0 neoantigens, for example, but 0 variant
        # collections will lead to NaN variants and neoantigens.
        if failed_io or len(variant_collections) == 0:
            print("Variants did not exist for patient %s" % patient.id)
            return None

        if len(variant_collections) == 1:
            # There is nothing to merge
            vcf_path, variants = variant_collections[0]
            merged_variants = variants
        else:
            merged_variants = self._merge_variant_collections(dict(variant_collections), self.merge_type)

        self.save_to_cache(merged_variants, self.cache_names["variant"], patient.id, cached_file_name)
        return filter_variants(variant_collection=merged_variants,
                               patient=patient,
                               filter_fn=filter_fn)

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

    def load_polyphen_annotations(self, as_dataframe=False,
                                  filter_fn=None):
        """Load a dataframe containing polyphen2 annotations for all variants

        Parameters
        ----------
        database_file : string, sqlite
            Path to the WHESS/Polyphen2 SQLite database.
            Can be downloaded and bunzip2"ed from http://bit.ly/208mlIU
        filter_fn : function
            Takes a FilterablePolyphen and returns a boolean.
            Only annotations returning True are preserved.
            Overrides default self.filter_fn. `None` passes through to self.filter_fn.

        Returns
        -------
        annotations
            Dictionary of patient_id to a DataFrame that contains annotations
        """
        filter_fn = first_not_none_param([filter_fn, self.filter_fn], no_filter)
        patient_annotations = {}
        for patient in self:
            annotations = self._load_single_patient_polyphen(
                patient,
                filter_fn=filter_fn)
            if annotations is not None:
                annotations["patient_id"] = patient.id
                patient_annotations[patient.id] = annotations
        if as_dataframe:
            return pd.concat(patient_annotations.values())
        return patient_annotations

    def _load_single_patient_polyphen(self, patient, filter_fn):
        cache_name = self.cache_names["polyphen"]
        cached_file_name = "polyphen-annotations.csv"

        # Don't filter here, as these variants are used to generate the
        # PolyPhen cache; and cached items are never filtered.
        variants = self._load_single_patient_variants(patient,
                                                      filter_fn=None)
        if variants is None:
            return None

        cached = self.load_from_cache(cache_name, patient.id, cached_file_name)
        if cached is not None:
            return filter_polyphen(polyphen_df=cached,
                                   variant_collection=variants,
                                   patient=patient,
                                   filter_fn=filter_fn)

        engine = create_engine("sqlite:///{}".format(self.polyphen_dump_path))
        conn = engine.connect()

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
        return filter_polyphen(polyphen_df=df,
                               variant_collection=variants,
                               patient=patient,
                               filter_fn=filter_fn)

    def load_effects(self, patients=None, only_nonsynonymous=False,
                     filter_fn=None):
        """Load a dictionary of patient_id to varcode.EffectCollection

        Note that this only loads one effect per variant.

        Parameters
        ----------
        patients : str, optional
            Filter to a subset of patients
        only_nonsynonymous : bool, optional
            If true, load only nonsynonymous effects, default False
        filter_fn : function
            Takes a FilterableEffect and returns a boolean. Only effects returning True are preserved.
            Overrides default self.filter_fn. `None` passes through to self.filter_fn.

        Returns
        -------
        effects
             Dictionary of patient_id to varcode.EffectCollection
        """
        filter_fn = first_not_none_param([filter_fn, self.filter_fn], no_filter)
        patient_effects = {}
        for patient in self.iter_patients(patients):
            effects = self._load_single_patient_effects(
                patient, only_nonsynonymous, filter_fn)
            if effects is not None:
                patient_effects[patient.id] = effects
        return patient_effects

    def _load_single_patient_effects(self, patient, only_nonsynonymous, filter_fn):
        cached_file_name = "%s-%s-effects.pkl" % (self.variant_type, self.merge_type)

        # Don't filter here, as these variants are used to generate the
        # effects cache; and cached items are never filtered.
        variants = self._load_single_patient_variants(patient, filter_fn=None)
        if variants is None:
            return None

        if only_nonsynonymous:
            cached = self.load_from_cache(self.cache_names["nonsynonymous_effect"], patient.id, cached_file_name)
        else:
            cached = self.load_from_cache(self.cache_names["effect"], patient.id, cached_file_name)
        if cached is not None:
            return filter_effects(effect_collection=cached,
                                  variant_collection=variants,
                                  patient=patient,
                                  filter_fn=filter_fn)

        effects = variants.effects()

        # Always take the top priority effect per variant so we end up with a single
        # effect per variant.
        nonsynonymous_effects = EffectCollection(
            effects.drop_silent_and_noncoding().top_priority_effect_per_variant().values())
        effects = EffectCollection(effects.top_priority_effect_per_variant().values())

        self.save_to_cache(effects, self.cache_names["effect"], patient.id, cached_file_name)
        self.save_to_cache(nonsynonymous_effects, self.cache_names["nonsynonymous_effect"], patient.id, cached_file_name)
        return filter_effects(
            effect_collection=(
                nonsynonymous_effects if only_nonsynonymous else effects),
            variant_collection=variants,
            patient=patient,
            filter_fn=filter_fn)

    def load_kallisto(self):
        """
        Load Kallisto transcript quantification data for a cohort

        Parameters
        ----------

        Returns
        -------
        kallisto_data : Pandas dataframe
            Pandas dataframe with Kallisto data for all patients
            columns include patient_id, gene_name, est_counts
        """
        kallisto_data = pd.concat(
            [self._load_single_patient_kallisto(patient) for patient in self],
            copy=False
        )

        if self.kallisto_ensembl_version is None:
            raise ValueError("Required a kallisto_ensembl_version but none was specified")

        ensembl_release = cached_release(self.kallisto_ensembl_version)

        kallisto_data['gene_name'] = \
            kallisto_data['target_id'].map(lambda t: ensembl_release.gene_name_of_transcript_id(t))

        # sum counts across genes
        kallisto_data = \
            kallisto_data.groupby(['patient_id', 'gene_name'])[['est_counts']].sum().reset_index()

        return kallisto_data

    def _load_single_patient_kallisto(self, patient):
        """
        Load Kallisto gene quantification given a patient

        Parameters
        ----------
        patient : Patient

        Returns
        -------
        data: Pandas dataframe
            Pandas dataframe of sample's Kallisto data
            columns include patient_id, target_id, length, eff_length, est_counts, tpm
        """
        data = pd.read_csv(patient.tumor_sample.kallisto_path, sep="\t")
        data["patient_id"] = patient.id
        return data

    def load_cufflinks(self, filter_ok=True):
        """
        Load a Cufflinks gene expression data for a cohort

        Parameters
        ----------
        filter_ok : bool, optional
            If true, filter Cufflinks data to row with FPKM_status == "OK"

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
            If true, filter Cufflinks data to row with FPKM_status == "OK"

        Returns
        -------
        data: Pandas dataframe
            Pandas dataframe of sample's Cufflinks data
            columns include patient_id, gene_id, gene_short_name, FPKM, FPKM_conf_lo, FPKM_conf_hi
        """
        data = pd.read_csv(patient.tumor_sample.cufflinks_path, sep="\t")
        data["patient_id"] = patient.id

        if filter_ok:
            # Filter to OK FPKM counts
            data = data[data["FPKM_status"] == "OK"]
        return data

    def load_neoantigens(self, patients=None, only_expressed=False,
                         epitope_lengths=[8, 9, 10, 11], ic50_cutoff=500,
                         process_limit=10, max_file_records=None,
                         filter_fn=None):
        filter_fn = first_not_none_param([filter_fn, self.filter_fn], no_filter)

        dfs = {}
        for patient in self.iter_patients(patients):
            df_epitopes = self._load_single_patient_neoantigens(
                patient=patient,
                only_expressed=only_expressed,
                epitope_lengths=epitope_lengths,
                ic50_cutoff=ic50_cutoff,
                process_limit=process_limit,
                max_file_records=max_file_records,
                filter_fn=filter_fn)
            if df_epitopes is not None:
                dfs[patient.id] = df_epitopes
        return dfs

    def _load_single_patient_neoantigens(self, patient, only_expressed, epitope_lengths,
                                         ic50_cutoff, process_limit, max_file_records,
                                         filter_fn):
        cached_file_name = "%s-%s-neoantigens.csv" % (self.variant_type, self.merge_type)

        # Don't filter here, as these variants are used to generate the
        # neoantigen cache; and cached items are never filtered.
        variants = self._load_single_patient_variants(patient, filter_fn=None)
        if variants is None:
            return None

        if patient.hla_alleles is None:
            print("HLA alleles did not exist for patient %s" % patient.id)
            return None

        if only_expressed:
            cached = self.load_from_cache(self.cache_names["expressed_neoantigen"], patient.id, cached_file_name)
        else:
            cached = self.load_from_cache(self.cache_names["neoantigen"], patient.id, cached_file_name)
        if cached is not None:
            return filter_neoantigens(neoantigens_df=cached,
                                      variant_collection=variants,
                                      patient=patient,
                                      filter_fn=filter_fn)

        mhc_model = NetMHCcons(
            alleles=patient.hla_alleles,
            epitope_lengths=epitope_lengths,
            max_file_records=max_file_records,
            process_limit=process_limit)
        if only_expressed:
            df_isovar = self.load_single_patient_isovar(patient=patient,
                                                         variants=variants,
                                                         epitope_lengths=epitope_lengths)

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
            # Store chr/pos/ref/alt in the cached DataFrame so we can filter based on
            # the variant later.
            for variant_column in ["chr", "pos", "ref", "alt"]:
                # Be consistent with Topiary's output of "start" rather than "pos".
                # Isovar, on the other hand, outputs "pos".
                # See https://github.com/hammerlab/topiary/blob/5c12bab3d47bd86d396b079294aff141265f8b41/topiary/converters.py#L50
                df_column = "start" if variant_column == "pos" else variant_column
                df_epitopes[df_column] = df_epitopes.source_sequence_key.apply(
                    lambda key: dict(key)[variant_column])
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

        return filter_neoantigens(neoantigens_df=df_epitopes,
                                  variant_collection=variants,
                                  patient=patient,
                                  filter_fn=filter_fn)

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

    def load_single_patient_isovar(self, patient, variants, epitope_lengths):
        # TODO: different epitope lengths, and other parameters, should result in
        # different caches
        isovar_cached_file_name = "%s-%s-isovar.csv" % (self.variant_type, self.merge_type)
        df_isovar = self.load_from_cache(self.cache_names["isovar"], patient.id, isovar_cached_file_name)
        if df_isovar is not None:
            return df_isovar

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
            min_mapping_quality=1)
        self.save_to_cache(df_isovar, self.cache_names["isovar"], patient.id, isovar_cached_file_name)
        return df_isovar

    def load_ensembl_coverage(self):
        if self.pageant_coverage_path is None:
            raise ValueError("Need a Pageant CoverageDepth path to load ensembl coverage values")
        return variant_filters.load_ensembl_coverage(
            cohort=self,
            coverage_path=self.pageant_coverage_path,
            min_normal_depth=self.min_coverage_normal_depth,
            min_tumor_depth=self.min_coverage_tumor_depth)

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

    def plot_roc_curve(self, on, bootstrap_samples=100, col=None, ax=None, **kwargs):
        """Plot an ROC curve for benefit and a given variable

        Parameters
        ----------
        on : str or function
            See `cohort.load.as_dataframe`
        bootstrap_samples : int, optional
            Number of boostrap samples to use to compute the AUC
        ax : Axes, default None
            Axes to plot on

        Returns
        -------
        (mean_auc_score, plot): (float, matplotlib plot)
            Returns the average AUC for the given predictor over `bootstrap_samples`
            and the associated ROC curve
        """
        plot_col, df = self.as_dataframe(on, col=col, **kwargs)
        df = filter_not_null(df, "benefit")
        df = filter_not_null(df, plot_col)
        df.benefit = df.benefit.astype(bool)
        return roc_curve_plot(df, plot_col, "benefit", bootstrap_samples, ax=ax)

    def plot_benefit(self, on, col=None, benefit_col="benefit", label="Response", ax=None,
                     alternative="two-sided", **kwargs):
        """Plot a comparison of benefit/response in the cohort on a given variable
        """
        return self.plot_boolean(on=on,
                                 boolean_col=benefit_col,
                                 col=col,
                                 alternative=alternative,
                                 boolean_label=label,
                                 boolean_value_map={True: "Benefit", False: "No Benefit"},
                                 order=["No Benefit", "Benefit"],
                                 ax=ax,
                                 **kwargs)

    def plot_boolean(self,
                     on,
                     boolean_col,
                     plot_col=None,
                     boolean_label=None,
                     boolean_value_map={},
                     col=None,
                     order=None,
                     ax=None,
                     alternative="two-sided",
                     **kwargs):
        """Plot a comparison of `boolean_col` in the cohort on a given variable via
        `on` or `col`.

        If the variable (through `on` or `col`) is binary this will compare
        odds-ratios and perform a Fisher's exact test.

        If the variable is numeric, this will compare the distributions through
        a Mann-Whitney test and plot the distributions with box-strip plot

        Parameters
        ----------
        on : str or function or list or dict
            See `cohort.load.as_dataframe`
        plot_col : str, optional
            If on has many columns, this is the one whose values we are plotting.
            If on has a single column, this is unnecessary.
        boolean_col : str
            Column name of boolean column to plot or compare against.
        boolean_label : None, optional
            Label to give boolean column in the plot
        boolean_value_map : dict, optional
            Map of conversions for values in the boolean column, i.e. {True: 'High', False: 'Low'}
        col : str, optional
            If specified, store the result of `on`. See `cohort.load.as_dataframe`
        order : None, optional
            Order of the labels on the x-axis
        ax : None, optional
            Axes to plot on
        alternative : str, optional
            Choose the sidedness of the mannwhitneyu or Fisher's Exact test.

        Returns
        -------
        (Test statistic, p-value): (float, float)

        """
        cols, df = self.as_dataframe(on, col, **kwargs)
        if type(cols) == str:
            if plot_col is not None:
                raise ValueError("plot_col is specified when it isn't ndeeded.")
            plot_col = cols
        elif type(cols) == list:
            if plot_col is None:
                raise ValueError("plot_col must be specified when multiple `on`s are present.")

        df = filter_not_null(df, boolean_col)
        df = filter_not_null(df, plot_col)

        if boolean_label:
            df[boolean_label] = df[boolean_col]
            boolean_col = boolean_label

        condition_value = None
        if boolean_value_map:
            assert set(boolean_value_map.keys()) == set([True, False]), \
                "Improper mapping of boolean column provided"
            df[boolean_col] = df[boolean_col].map(lambda v: boolean_value_map[v])
            condition_value = boolean_value_map[True]

        if df[plot_col].dtype == "bool":
            results = fishers_exact_plot(
                data=df,
                condition1=boolean_col,
                condition2=plot_col,
                condition1_value=condition_value,
                alternative=alternative,
                order=order,
                ax=ax)
        else:
            results = mann_whitney_plot(
                data=df,
                condition=boolean_col,
                distribution=plot_col,
                condition_value=condition_value,
                alternative=alternative,
                order=order,
                ax=ax)
        return results

    def plot_survival(self, on, col=None, how="os", ax=None,
                      threshold=None, **kwargs):
        """Plot a Kaplan Meier survival curve by splitting the cohort into two groups
        Parameters
        ----------
        on : str or function
            See `cohort.load.as_dataframe`
        col : str, optional
            If specified, store the result of `on`. See `cohort.load.as_dataframe`
        how : {"os", "pfs"}, optional
            Whether to plot OS (overall survival) or PFS (progression free survival)
        threshold : int or "median", optional
            Threshold of `col` on which to split the cohort
        """
        assert how in ["os", "pfs"], "Invalid choice of survival plot type %s" % how
        plot_col, df = self.as_dataframe(on, col, **kwargs)
        df = filter_not_null(df, plot_col)
        if df[plot_col].dtype == "bool":
            default_threshold = None
        else:
            default_threshold = "median"
        results = plot_kmf(
            df=df,
            condition_col=plot_col,
            xlabel="Days",
            ylabel="Overall Survival (%)" if how == "os" else "Progression-Free Survival (%)",
            censor_col="deceased" if how == "os" else "progressed_or_deceased",
            survival_col=how,
            threshold=threshold if threshold is not None else default_threshold,
            ax=ax)
        return results

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
        for plot_col in plot_cols:
            df = filter_not_null(df, plot_col)
        p = sb.jointplot(data=df, x=plot_cols[0], y=plot_cols[1])
        return p

    def _list_patient_ids(self):
        """ Utility function to return a list of patient ids in the Cohort
        """
        results = []
        for patient in self:
            results.append(patient.id)
        return(results)

    def summarize_provenance_per_cache(self):
        """Utility function to summarize provenance files for cached items used by a Cohort,
        for each cache_dir that exists. Only existing cache_dirs are summarized.

        This is a summary of provenance files because the function checks to see whether all
        patients data have the same provenance within the cache dir. The function assumes
        that it will be desireable to have all patients data generated using the same
        environment, for each cache type.

        At the moment, most PROVENANCE files contain details about packages used to generat
        e the cached data file. However, this function is generic & so it summarizes the
        contents of those files irrespective of their contents.

        Returns
        ----------
        Dict containing summarized provenance for each existing cache_dir, after checking
        to see that provenance files are identical among all patients in the data frame for
        that cache_dir.

        If conflicting PROVENANCE files are discovered within a cache-dir:
         - a warning is generated, describing the conflict
         - and, a value of `None` is returned in the dictionary for that cache-dir

        See also
        -----------
        * `?cohorts.Cohort.summarize_provenance` which summarizes provenance files among
        cache_dirs.
        * `?cohorts.Cohort.summarize_dataframe` which hashes/summarizes contents of the data
        frame for this cohort.
        """
        provenance_summary = {}
        df = self.as_dataframe()
        for cache in self.cache_names:
            cache_name = self.cache_names[cache]
            cache_provenance = None
            num_discrepant = 0
            this_cache_dir = path.join(self.cache_dir, cache_name)
            if path.exists(this_cache_dir):
                for patient_id in self._list_patient_ids():
                    patient_cache_dir = path.join(this_cache_dir, patient_id)
                    try:
                        this_provenance = self.load_provenance(patient_cache_dir = patient_cache_dir)
                    except:
                        this_provenance = None
                    if this_provenance:
                        if not(cache_provenance):
                            cache_provenance = this_provenance
                        else:
                            num_discrepant += compare_provenance(this_provenance, cache_provenance)
                if num_discrepant == 0:
                    provenance_summary[cache_name] = cache_provenance
                else:
                    provenance_summary[cache_name] = None
        return(provenance_summary)

    def summarize_dataframe(self):
        """Summarize default dataframe for this cohort using a hash function.
        Useful for confirming the version of data used in various reports, e.g. ipynbs
        """
        if self.dataframe_hash:
            return(self.dataframe_hash)
        else:
            df = self._as_dataframe_unmodified()
            return(self.dataframe_hash)

    def summarize_provenance(self):
        """Utility function to summarize provenance files for cached items used by a Cohort.

        At the moment, most PROVENANCE files contain details about packages used to
        generate files. However, this function is generic & so it summarizes the contents
        of those files irrespective of their contents.

        Returns
        ----------
        Dict containing summary of provenance items, among all cache dirs used by the Cohort.

        IE if all provenances are identical across all cache dirs, then a single set of
        provenances is returned. Otherwise, if all provenances are not identical, the provenance
        items per cache_dir are returned.

        See also
        ----------
        `?cohorts.Cohort.summarize_provenance_per_cache` which is used to summarize provenance
        for each existing cache_dir.
        """
        provenance_per_cache = self.summarize_provenance_per_cache()
        summary_provenance = None
        num_discrepant = 0
        for cache in provenance_per_cache:
            if not(summary_provenance):
                ## pick arbitrary provenance & call this the "summary" (for now)
                summary_provenance = provenance_per_cache[cache]
                summary_provenance_name = cache
            ## for each cache, check equivalence with summary_provenance
            num_discrepant += compare_provenance(
                provenance_per_cache[cache],
                summary_provenance,
                left_outer_diff = "In %s but not in %s" % (cache, summary_provenance_name),
                right_outer_diff = "In %s but not in %s" % (summary_provenance_name, cache)
            )
        ## compare provenance across cached items
        if num_discrepant == 0:
            prov = summary_provenance ## report summary provenance if exists
        else:
            prov = provenance_per_cache ## otherwise, return provenance per cache
        return(prov)

    def summarize_data_sources(self):
        """Utility function to summarize data source status for this Cohort, useful for confirming
        the state of data used for an analysis

        Returns
        ----------
        Dictionary with summary of data sources

        Currently contains
        - dataframe_hash: hash of the dataframe (see `?cohorts.Cohort.summarize_dataframe`)
        - provenance_file_summary: summary of provenance file contents (see `?cohorts.Cohort.summarize_provenance`)
        """
        provenance_file_summary = self.summarize_provenance()
        dataframe_hash = self.summarize_dataframe()
        results = {
            "provenance_file_summary": provenance_file_summary,
            "dataframe_hash": dataframe_hash
        }
        return(results)

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

def _provenance_str(provenance):
    """Utility function used by compare_provenance to print diff
    """
    return ["%s==%s" % (key, value) for (key, value) in provenance]

def compare_provenance(
        this_provenance, other_provenance,
        left_outer_diff = "In current but not comparison",
        right_outer_diff = "In comparison but not current"):
    """Utility function to compare two abritrary provenance dicts
    returns number of discrepancies.

    Parameters
    ----------
    this_provenance: provenance dict (to be compared to "other_provenance")
    other_provenance: comparison provenance dict

    (optional)
    left_outer_diff: description/prefix used when printing items in this_provenance but not in other_provenance
    right_outer_diff: description/prefix used when printing items in other_provenance but not in this_provenance

    Returns
    -----------
    Number of discrepancies (0: None)
    """
    ## if either this or other items is null, return 0
    if (not this_provenance or not other_provenance):
        return 0
    
    this_items = set(this_provenance.items())
    other_items = set(other_provenance.items())

    # Two-way diff: are any modules introduced, and are any modules lost?
    new_diff = this_items.difference(other_items)
    old_diff = other_items.difference(this_items)
    warn_str = ""
    if len(new_diff) > 0:
        warn_str += "%s: %s" % (
            left_outer_diff,
            _provenance_str(new_diff))
    if len(old_diff) > 0:
        warn_str += "%s: %s" % (
            right_outer_diff,
            _provenance_str(old_diff))

    if len(warn_str) > 0:
        warnings.warn(warn_str, Warning)

    return(len(new_diff)+len(old_diff))

