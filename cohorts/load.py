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

from os import path, makedirs, listdir
from shutil import rmtree
import pandas as pd
import numpy as np

# pylint doesn't like this line
# pylint: disable=no-name-in-module
import six.moves.cPickle as pickle
from types import FunctionType
from collections import defaultdict

import varcode
from varcode import VariantCollection, EffectCollection
from mhctools import NetMHCcons, EpitopeCollection
from topiary import predict_epitopes_from_variants, epitopes_to_dataframe
from topiary.sequence_helpers import contains_mutant_residues
from isovar.protein_sequence import variants_to_protein_sequences_dataframe
from pysam import AlignmentFile

from .survival import plot_kmf
from .plot import mann_whitney_plot, fishers_exact_plot
from .collection import Collection

class InvalidDataError(ValueError):
    pass

def require_id_str(id):
    if type(id) != str:
        raise ValueError("Expected ID string, but id = %s" % str(id))

class Sample(object):
    """
    Represents a single tumor or normal sample. It can point to DNA and/or
    RNA reads.
    """
    def __init__(self,
                 id,
                 is_tumor,
                 bam_path_dna=None,
                 bam_path_rna=None,
                 cufflinks_path=None):
        require_id_str(id)
        self.id = id
        self.is_tumor = is_tumor
        self.bam_path_dna = bam_path_dna
        self.bam_path_rna = bam_path_rna
        self.cufflinks_path = cufflinks_path

    def __hash__(self):
        return hash(self.id)

    def __eq__(self, other):
        if self is other:
            return True
        return (
            self.__class__ == other.__class__ and
            self.id == other.id)

class PairedSample(object):
    """
    Represents a tumor-normal sample pair, as well as pointing to data derived
    from that pair (e.g. VCFs).
    """
    def __init__(self,
                 id,
                 snv_vcf_paths=[],
                 indel_vcf_paths=[],
                 normal_sample=None,
                 tumor_sample=None):
        require_id_str(id)
        self.id = id
        self.snv_vcf_paths = snv_vcf_paths
        self.indel_vcf_paths = indel_vcf_paths
        self.normal_sample = normal_sample
        self.tumor_sample = tumor_sample
        self.samples = [normal_sample, tumor_sample]

    def __hash__(self):
        return hash(self.id)

    def __eq__(self, other):
        if self is other:
            return True
        return (
            self.__class__ == other.__class__ and
            self.id == other.id)

class Patient(object):
    """
    Represents a patient, which contains one or more `PairedSample`s.
    """
    def __init__(self,
                 id,
                 os,
                 pfs,
                 deceased,
                 progressed=None,
                 progressed_or_deceased=None,
                 benefit=None,
                 paired_samples=[],
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
        self.paired_samples = paired_samples
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

class CohortDataFrame(object):
    """
    Wraps a DataFrame with some information on how to join it.
    """
    def __init__(self,
                 name,
                 load_dataframe,
                 group_by,
                 id_col):
        self.name = name
        self.load_dataframe = load_dataframe
        self.group_by = group_by
        self.id_col = id_col

class Cohort(Collection):
    """
    Represents a cohort of `Patient`s.
    """
    def __init__(self,
                 patients,
                 cache_dir,
                 cache_results=True,
                 extra_dataframes=[],
                 join_with=None,
                 join_how="outer"):
        Collection.__init__(
            self,
            elements=patients)
        self.cache_dir = cache_dir
        self.cache_results = cache_results

        joinable_dataframes = [
            CohortDataFrame("cufflinks", self.load_cufflinks, group_by="patient", id_col="sample_id")]
        joinable_dataframes.extend(extra_dataframes)
        self.joinable_dataframes = joinable_dataframes
        self.join_with = join_with
        self.join_how =  join_how

        self.verify_id_uniqueness()
        self.verify_survival()

        self.cache_names = {"variant": "cached-variants",
                            "effect": "cached-effects",
                            "nonsynonymous_effect": "cached-nonsynonymous-effects",
                            "neoantigen": "cached-neoantigens",
                            "expressed_neoantigen": "cached-expressed-neoantigens",
                            "isovar": "cached-isovar-output"}
        self.verify_cache(self.cache_names)

    def verify_id_uniqueness(self):
        patient_ids = set()
        paired_sample_ids = set()
        tumor_sample_ids = set()
        normal_sample_ids = set()
        for patient in self.elements:
            patient_ids.add(patient.id)
            for paired_sample in patient.paired_samples:
                paired_sample_ids.add(paired_sample.id)
                tumor_sample_ids.add(paired_sample.tumor_sample)
                normal_sample_ids.add(paired_sample.normal_sample)
        if len(patient_ids) != len(self.elements):
            raise ValueError("Non-unique patient IDs")
        if len(paired_sample_ids) != len(list(self.iter_paired_samples())):
            raise ValueError("Non-unique paired sample IDs")
        if len(tumor_sample_ids) != len(list(self.iter_tumor_samples())):
            raise ValueError("Non-unique tumor sample IDs")
        if len(normal_sample_ids) != len(list(self.iter_normal_samples())):
            raise ValueError("Non-unique normal sample IDs")

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

    def verify_cache(self, cache_names):
        bad_caches = []
        for cache_name in cache_names.values():
            cache_dir = path.join(self.cache_dir, cache_name)
            if path.exists(cache_dir):
                cache_subdirs = set(listdir(cache_dir))
                cache_int_subdirs = set([int(name) for name in cache_subdirs])
                if len(cache_subdirs) != len(cache_int_subdirs):
                    bad_caches.append(cache_name)

        if len(bad_caches) > 0:
            raise ValueError("Caches %s have duplicate int/str directories" %
                             str(bad_caches))

    def iter_paired_samples(self, include_parents=False):
        for patient in self.elements:
            for paired_sample in patient.paired_samples:
                if include_parents:
                    yield (patient, paired_sample)
                else:
                    yield paired_sample

    def iter_tumor_samples(self, include_parents=False):
        return self._iter_samples_attr(attr="tumor_sample",
                                       include_parents=include_parents)

    def iter_normal_samples(self, include_parents=False):
        return self._iter_samples_attr(attr="normal_sample",
                                       include_parents=include_parents)

    def _iter_samples_attr(self, attr, **kwargs):
        for patient in self.elements:
            for paired_sample in patient.paired_samples:
                sample = getattr(paired_sample, attr)
                if kwargs["include_parents"]:
                    yield (patient, paired_sample, sample)
                else:
                    yield sample

    def as_dataframe(self, group_by="patient",
                     join_with=None, join_how=None):
        verify_group_by(group_by)

        df = pd.DataFrame()
        join_dataframes = []

        if not join_with:
            current_join_with = self.join_with
        elif type(join_with) == list:
            current_join_with = join_with
        else:
            current_join_with = [join_with]

        current_join_how = self.join_how if join_how is None else join_how
        current_join_how = "outer" if current_join_how is None else current_join_how

        for dataframe in self.joinable_dataframes:
            if dataframe.name in current_join_with:
                if dataframe.group_by != group_by:
                    raise ValueError("Requested extra DataFrame %s, but wrong group_by %s" % (dataframe.name, dataframe.group_by))
                join_dataframes.append(dataframe)

        if group_by == "paired_sample":
            patient_samples = [(patient, sample) for (patient, sample)
                               in self.iter_paired_samples(include_parents=True)]
            df["patient_id"] = pd.Series([patient.id for (patient, sample)
                                          in patient_samples])
            df["sample_id"] = pd.Series([sample.id for (patient, sample)
                                         in patient_samples])
            patients = [patient for (patient, sample) in patient_samples]
            num_patients_with_samples = len(set(patients))
            num_patients = len(set(self.elements))
            num_patients_without_samples = len(set(self.elements).difference(set(patients)))
            if num_patients != num_patients_with_samples:
                print("Missing samples for %d patients: from %d to %d" % (
                    num_patients_without_samples, num_patients, num_patients_with_samples))
        else:
            patients = self.elements
            df["patient_id"] = pd.Series([patient.id for patient in patients])

        for clinical_col in ["benefit", "os", "pfs", "deceased",
                             "progressed", "progressed_or_deceased"]:
            df[clinical_col] = pd.Series([getattr(patient, clinical_col) for patient in patients])

        additional_data_all_patients = defaultdict(list)
        for patient in patients:
            if patient.additional_data is not None:
                for key, value in patient.additional_data.items():
                    additional_data_all_patients[key].append(value)

        if len(additional_data_all_patients) > 0:
            df = df.merge(pd.DataFrame(additional_data_all_patients), on="patient_id", how="left")

        for dataframe in join_dataframes:
            left_on = "sample_id" if group_by == "paired_sample" else "patient_id"
            old_len_df = len(df)
            df = df.merge(
                dataframe.load_dataframe(),
                left_on=left_on,
                right_on=dataframe.id_col,
                how=current_join_how)
            print("%s join with %s: %d to %d rows" % (
                current_join_how,
                dataframe.name,
                old_len_df,
                len(df)))
        return df

    def load_from_cache(self, cache_name, sample_id, file_name):
        if not self.cache_results:
            return None

        cache_dir = path.join(self.cache_dir, cache_name)
        sample_cache_dir = path.join(cache_dir, str(sample_id))
        cache_file = path.join(sample_cache_dir, file_name)

        if not path.exists(cache_file):
            return None

        if path.splitext(cache_file)[1] == ".csv":
            return pd.read_csv(cache_file, dtype={"sample_id": object,
                                                  "patient_id": object})
        else:
            with open(cache_file, "rb") as f:
                return pickle.load(f)

    def save_to_cache(self, obj, cache_name, sample_id, file_name):
        if not self.cache_results:
            return

        cache_dir = path.join(self.cache_dir, cache_name)
        sample_cache_dir = path.join(cache_dir, str(sample_id))
        cache_file = path.join(sample_cache_dir, file_name)

        if not path.exists(sample_cache_dir):
            makedirs(sample_cache_dir)

        if type(obj) == pd.DataFrame:
            obj.to_csv(cache_file, index=False)
        else:
            with open(cache_file, "wb") as f:
                pickle.dump(obj, f)

    def load_variants(self, variant_type="snv", merge_type="union"):
        assert variant_type in ["snv", "indel"], "Unknown variant type: %s" % variant_type
        sample_variants = {}

        for sample in self.iter_paired_samples():
            try:
                variants = self._load_single_sample_variants(
                    sample, variant_type, merge_type)
            except IOError:
                print("Variants did not exist for %s" % sample)
                continue

            sample_variants[sample.id] = variants
        return sample_variants

    def _load_single_sample_variants(self, sample, variant_type, merge_type):
        cached_file_name = "%s-%s-variants.pkl" % (variant_type, merge_type)
        cached = self.load_from_cache(self.cache_names["variant"], sample.id, cached_file_name)
        if cached is not None:
            return cached

        combined_variants = []
        vcf_paths = sample.snv_vcf_paths if variant_type == "snv" else sample.indel_vcf_paths
        for vcf_path in vcf_paths:
            variants = varcode.load_vcf_fast(vcf_path)
            combined_variants.append(set(variants.elements))

        if len(combined_variants) == 1:
            # There is nothing to merge
            merged_variants =  VariantCollection(combined_variants[0])
        else:
            assert merge_type in ["union", "intersection"], "Unknown merge type: %s" % merge_type
            if merge_type == "union":
                merged_variants = VariantCollection(set.union(*combined_variants))
            elif merge_type == "intersection":
                merged_variants = VariantCollection(set.intersection(*combined_variants))

        self.save_to_cache(merged_variants, self.cache_names["variant"], sample.id, cached_file_name)

        return merged_variants

    def load_effects(self, only_nonsynonymous=False, variant_type="snv", merge_type="union"):
        sample_effects = {}
        for sample in self.iter_paired_samples():
            effects = self._load_single_sample_effects(
                sample, only_nonsynonymous, variant_type, merge_type)
            sample_effects[sample.id] = effects
        return sample_effects

    def _load_single_sample_effects(self, sample, only_nonsynonymous,
                                    variant_type, merge_type):
        cached_file_name = "%s-%s-effects.pkl" % (variant_type, merge_type)
        if only_nonsynonymous:
            cached = self.load_from_cache(self.cache_names["nonsynonymous_effect"], sample.id, cached_file_name)
        else:
            cached = self.load_from_cache(self.cache_names["effect"], sample.id, cached_file_name)
        if cached is not None:
            return cached

        variants = self._load_single_sample_variants(
            sample, variant_type, merge_type)
        effects = variants.effects()
        nonsynonymous_effects = EffectCollection(
            effects.drop_silent_and_noncoding().top_priority_effect_per_variant().values())

        self.save_to_cache(effects, self.cache_names["effect"], sample.id, cached_file_name)
        self.save_to_cache(nonsynonymous_effects, self.cache_names["nonsynonymous_effect"], sample.id, cached_file_name)

        if only_nonsynonymous:
            return nonsynonymous_effects
        return effects

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
            Pandas dataframe with Cufflinks data for all samples
            columns include sample_id, gene_id, gene_short_name, FPKM, FPKM_conf_lo, FPKM_conf_hi
        """
        return \
            pd.concat(
                [self._load_single_sample_cufflinks(sample, filter_ok)
                 for sample in self.iter_tumor_samples()],
                copy=False
        )

    def _load_single_sample_cufflinks(self, sample, filter_ok):
        """
        Load Cufflinks gene quantification given a sample_id

        Parameters
        ----------
        sample : Sample
        filter_ok : bool, optional
            If true, filter Cufflinks data to row with FPKM_status == 'OK'

        Returns
        -------
        cufflinks_data: Pandas dataframe
            Pandas dataframe of sample's Cufflinks data
            columns include sample_id, gene_id, gene_short_name, FPKM, FPKM_conf_lo, FPKM_conf_hi
        """
        data = pd.read_csv(sample.cufflinks_path, sep="\t")
        data["sample_id"] = sample.id

        if filter_ok:
            # Filter to OK FPKM counts
            data = data[data["FPKM_status"] == "OK"]
        return data

    def load_neoantigens(self, variant_type="snv", merge_type="union",
                         only_expressed=False, epitope_lengths=[8, 9, 10, 11],
                         ic50_cutoff=500, process_limit=10, max_file_records=None):
        dfs = []
        for patient, sample in self.iter_paired_samples(include_parents=True):
            df_epitopes = self._load_single_sample_neoantigens(
                patient=patient,
                sample=sample,
                variant_type=variant_type,
                merge_type=merge_type,
                only_expressed=only_expressed,
                epitope_lengths=epitope_lengths,
                ic50_cutoff=ic50_cutoff,
                process_limit=process_limit,
                max_file_records=max_file_records)
            dfs.append(df_epitopes)
        return pd.concat(dfs)

    def _load_single_sample_neoantigens(self, patient, sample, variant_type,
                                        merge_type, only_expressed, epitope_lengths,
                                        ic50_cutoff, process_limit, max_file_records):
        cached_file_name = "%s-%s-neoantigens.csv" % (variant_type, merge_type)
        if only_expressed:
            cached = self.load_from_cache(self.cache_names["expressed_neoantigen"], sample.id, cached_file_name)
        else:
            cached = self.load_from_cache(self.cache_names["neoantigen"], sample.id, cached_file_name)
        if cached is not None:
            return cached

        variants = self._load_single_sample_variants(
            sample, variant_type, merge_type)
        mhc_model = NetMHCcons(
            alleles=patient.hla_alleles,
            epitope_lengths=epitope_lengths,
            max_file_records=max_file_records,
            process_limit=process_limit)
        if only_expressed:
            isovar_cached_file_name = "%s-%s-isovar.csv" % (variant_type, merge_type)
            df_isovar = self.load_from_cache(self.cache_names["isovar"], sample.id, isovar_cached_file_name)
            if df_isovar is None:
                import logging
                logging.disable(logging.INFO)
                rna_bam_file = AlignmentFile(sample.rna_bam.path)
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
                self.save_to_cache(df_isovar, self.cache_names["isovar"], sample.id, isovar_cached_file_name)

            # Map from isovar rows to protein sequences
            isovar_rows_to_protein_sequences = dict([
                (frozenset(row.to_dict().items()), row["amino_acids"]) for (i, row) in df_isovar.iterrows()])

            # MHC binding prediction
            epitopes = mhc_model.predict(isovar_rows_to_protein_sequences)

            # Only include peptides that overlap a variant; without this filter, when we use
            # protein_sequence_length above, some 8mers generated from a 21mer source will
            # not overlap a variant.
            df_epitopes = self.get_filtered_isovar_epitopes(
                epitopes, ic50_cutoff=ic50_cutoff).dataframe()
            df_epitopes["sample_id"] = sample.id

            self.save_to_cache(df_epitopes, self.cache_names["expressed_neoantigen"], sample.id, cached_file_name)
        else:
            epitopes = predict_epitopes_from_variants(
                variants=variants,
                mhc_model=mhc_model,
                ic50_cutoff=ic50_cutoff,
                # Only include peptides with a variant
                only_novel_epitopes=True)
            df_epitopes = epitopes_to_dataframe(epitopes)
            df_epitopes["sample_id"] = sample.id

            self.save_to_cache(df_epitopes, self.cache_names["neoantigen"], sample.id, cached_file_name)

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

    def clinical_func(self, row_func, col):
        df = self.as_dataframe()
        df[col] = df.apply(row_func, axis=1)
        return col, df

    def plot_init(self, on, col, col_equals, **kwargs):
        """
        `on` is either:
        - a function that creates a new column for comparison, e.g. count.snv_count
        - a function that takes a clinical dataframe row as input and returns a boolean or quantity based on that row (with col being the name of that new boolean or quantity)
        - a string representing an existing column
        - a string representing a new column name, created by comparing column `col` with the value `col_equals`
        """
        if type(on) == FunctionType:
            try:
                return on(self)
            except TypeError:
                col = col if col is not None else "untitled"
                return self.clinical_func(on, col)
        if type(on) == str:
            cohort_dataframe = self.as_dataframe(**kwargs)
            if on in cohort_dataframe.columns:
                return (on, cohort_dataframe)
            else:
                return col_func(self, on, col, col_equals)

    def plot_benefit(self, on, col=None, col_equals=None, **kwargs):
        """Plot a comparison of benefit/response in the cohort on a given variable

        If the variable (through `on` or `col` is binary) this will compare
        odds-ratios and perform a Fisher's exact test.

        If the variable is numeric, this will compare the distributions through
        a Mann-Whitney test and plot the distributions with box-strip plot

        Parameters
        ----------
        on : See `cohort.load.plot_init`
        col : str, optional
            If specified, store the result of `on` or compare the column specified with `col_equals`
            See `cohort.load.plot_init`
        col_equals : str, optional
            Fixed value of `col` on which to split the cohort

        Returns
        -------
        (Test statistic, p-value): (float, float)

        """
        plot_col, df = self.plot_init(on, col, col_equals, **kwargs)
        original_len = len(df)
        df = df[df.benefit.notnull()]
        updated_len = len(df)
        df.benefit = df.benefit.apply(bool)
        if updated_len < original_len:
            print("Missing benefit for %d samples: from %d to %d" % (original_len - updated_len, original_len, updated_len))
        if df[plot_col].dtype == "bool":
            results = fishers_exact_plot(
                data=df,
                condition1="benefit",
                condition2=plot_col)
        else:
            results = mann_whitney_plot(
                data=df,
                condition="benefit",
                distribution=plot_col)

        return results

    def plot_survival(self, on, col=None, col_equals=None, how="os",
                      threshold=None, **kwargs):
        """Plot a Kaplan Meier survival curve by splitting the cohort into two groups

        Parameters
        ----------
        on :  See `cohort.load.plot_init`
        col : str, optional
            If specified, store the result of `on` or compare the column specified with `col_equals`
            See `cohort.load.plot_init`
        col_equals : str, optional
            Fixed value of `col` on which to split the cohort
        how : {'os', 'pfs'}, optional 
            Whether to plot OS (overall survival) or PFS (progression free survival)
        threshold : int or 'median', optional
            Threshold of `col` on which to split the cohort
        """
        assert how in ["os", "pfs"], "Invalid choice of survival plot type %s" % how
        plot_col, df = self.plot_init(on, col, col_equals, **kwargs)
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

def col_func(cohort, on, col, col_equals):
    df = cohort.as_dataframe()
    df[on] = df[col] == col_equals
    return on, df

def verify_group_by(group_by):
    if group_by not in ["patient", "paired_sample"]:
        raise ValueError("Invalid group_by: %s" % group_by)
