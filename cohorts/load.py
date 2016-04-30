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
import numpy as np

# pylint doesn't like this line
# pylint: disable=no-name-in-module
import six.moves.cPickle as pickle
from types import FunctionType
from collections import defaultdict

import varcode
from varcode import VariantCollection, EffectCollection, Collection
from mhctools import NetMHCcons, EpitopeCollection
from topiary import predict_epitopes_from_variants, epitopes_to_dataframe
from topiary.sequence_helpers import contains_mutant_residues
from isovar.protein_sequence import variants_to_protein_sequences_dataframe
from pysam import AlignmentFile

from .survival import plot_kmf
from .plot import mann_whitney_plot, fishers_exact_plot

class InvalidDataError(ValueError):
    pass

def require_id_str(id):
    if type(id) != str:
        raise ValueError("Expected ID string, but id = %s" % str(id))

class DataFile(object):
    def __init__(self,
                 id,
                 path):
        self.id = id
        self.path = path

class Sample(object):
    """
    Tumor vs normal
    DNA vs. RNA

    one sample is tumor
    one sample is normal

    but VCFS are both tumor and normal..
     * though they can be for multiple samples
    """
    def __init__(self,
                 id,
                 is_tumor,
                 dna_bam=None,
                 rna_bam=None):
        require_id_str(id)
        self.id = id
        self.is_tumor = is_tumor
        self.dna_bam = dna_bam
        self.rna_bam = rna_bam

class PairedSample(object):
    def __init__(self,
                 id,
                 vcfs=None,
                 normal_sample=None,
                 tumor_sample=None):
        require_id_str(id)
        self.id = id
        self.vcfs = vcfs
        self.normal_sample = normal_sample
        self.tumor_sample = tumor_sample
        self.samples = [normal_sample, tumor_sample]

class Patient(object):
    def __init__(self,
                 id,
                 os,
                 pfs,
                 deceased,
                 progressed=None,
                 progressed_or_deceased=None,
                 benefit=None,
                 paired_samples=[],
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

class Cohort(Collection):
    """Represents a cohort of patients."""

    def __init__(self,
                 patients,
                 cache_dir,
                 cache_results=True):
        Collection.__init__(
            self,
            elements=patients,
            path=None,
            distinct=False,
            sort_key=None)

        self.verify_id_uniqueness()

        self.cache_dir = cache_dir
        self.cache_results = cache_results
        self.clinical_dataframe = self.as_dataframe()

        self.verify_survival()

        self.cache_names = {"variant": "cached-variants",
                            "effect": "cached-effects",
                            "nonsynonymous_effect": "cached-nonsynonymous-effects",
                            "neoantigen": "cached-neoantigens",
                            "expressed_neoantigen": "cached-expressed-neoantigens",
                            "isovar": "cached-isovar-output"}

    def verify_id_uniqueness(self):
        patient_ids = set()
        paired_sample_ids = set()
        sample_ids = set()
        for patient in self.elements:
            patient_ids.add(patient.id)
            for paired_sample in patient.paired_samples:
                paired_sample_ids.add(paired_sample.id)
                for sample in paired_sample.samples:
                    sample_ids.add(sample.id)
        if len(patient_ids) != len(self.elements):
            raise ValueError("Non-unique patient IDs")
        if len(paired_sample_ids) != len(list(self.iter_paired_samples())):
            raise ValueError("Non-unique paired sample IDs")
        if len(sample_ids) != len(list(self.iter_samples())):
            raise ValueError("Non-unique sample IDs")

    def iter_paired_samples(self):
        for patient in self.elements:
            for paired_sample in patient.paired_samples:
                yield paired_sample

    def iter_samples(self):
        for patient in self.elements:
            for paired_sample in patient.paired_samples:
                for sample in paierd_sample.samples:
                    yield sample

    def as_dataframe(self):
        df =  pd.DataFrame({"patient_id": [patient.id for patient in self.elements],
                            "benefit": [patient.benefit for patient in self.elements],
                            "os": [patient.os for patient in self.elements],
                            "pfs": [patient.pfs for patient in self.elements],
                            "deceased": [patient.deceased for patient in self.elements],
                            "progressed": [patient.progressed for patient in self.elements],
                            "progressed_or_deceased": [patient.progressed_or_deceased
                                                       for patient in self.elements]})

        additional_data_all_patients = defaultdict(list)
        for patient in self.elements:
            if patient.additional_data is not None:
                for key, value in patient.additional_data.items():
                    additional_data_all_patients[key].append(value)

        if len(additional_data_all_patients) > 0:
            return df.merge(pd.DataFrame(additional_data_all_patients), on="patient_id", how="left")
        return df

    def verify_survival(self):
        if not (self.clinical_dataframe["pfs"] <=
                self.clinical_dataframe["os"]).all():
            raise InvalidDataError("PFS should be <= OS, but PFS is larger than OS for some patients.")

        def func(row):
            if row["pfs"] < row["os"]:
                if not row["progressed_or_deceased"]:
                    raise InvalidDataError(
                        "A patient did not progress despite PFS being less than OS. "
                        "Full row: %s" % row)
        self.clinical_dataframe.apply(func, axis=1)

    def load_from_cache(self, cache_name, sample_id, file_name):
        if not self.cache_results:
            return None

        cache_dir = path.join(self.cache_dir, cache_name)
        sample_cache_dir = path.join(cache_dir, str(sample_id))
        cache_file = path.join(sample_cache_dir, file_name)

        if not path.exists(cache_file):
            return None

        if path.splitext(cache_file)[1] == ".csv":
            return pd.read_csv(cache_file)
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
        for vcf in sample.vcfs:
            variants = varcode.load_vcf_fast(vcf.path)
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

    def load_neoantigens(self, variant_type="snv", merge_type="union",
                         only_expressed=False, epitope_lengths=[8, 9, 10, 11],
                         ic50_cutoff=500, process_limit=10, max_file_records=None):
        assert self.hla_alleles is not None, "Cannot predict neoantigens without HLA alleles"

        dfs = []
        for sample in self.iter_paired_samples():
            df_epitopes = self._load_single_sample_neoantigens(
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

    def _load_single_sample_neoantigens(self, sample, variant_type, merge_type,
                                        only_expressed, epitope_lengths,
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
        hla_alleles = self.sample.patient.hla
        mhc_model = NetMHCcons(
            alleles=hla_alleles,
            epitope_lengths=epitope_lengths,
            max_file_records=max_file_records,
            process_limit=process_limit)
        if only_expressed:
            df_isovar = self.load_from_cache(self.cache_names["isovar"], sample.id, cached_file_name)
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
                self.save_to_cache(df_isovar, self.cache_names["isovar"], sample.id, cached_file_name)

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

    def clinical_columns(self):
        column_types = [self.clinical_dataframe[col].dtype for col in self.clinical_dataframe.columns]
        return dict(zip(list(self.clinical_dataframe.columns), column_types))

    def clinical_func(self, row_func, col):
        df = self.clinical_dataframe.copy()
        df[col] = df.apply(row_func, axis=1)
        return col, df

    def plot_init(self, on, col, col_equals):
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
            except:
                col = col if col is not None else "untitled"
                return self.clinical_func(on, col)
        if type(on) == str:
            if on in self.clinical_dataframe.columns:
                return (on, self.clinical_dataframe)
            else:
                return col_func(self, on, col, col_equals)

    def plot_benefit(self, on, col=None, col_equals=None):
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
        plot_col, df = self.plot_init(on, col, col_equals)
        original_len = len(df)
        df = df[df[self.benefit_col].notnull()]
        updated_len = len(df)
        df[self.benefit_col] = df[self.benefit_col].apply(bool)
        if updated_len < original_len:
            print("Missing benefit for %d samples: from %d to %d" % (original_len - updated_len, original_len, updated_len))
        if df[plot_col].dtype == "bool":
            results = fishers_exact_plot(
                data=df,
                condition1=self.benefit_col,
                condition2=plot_col)
        else:
            results = mann_whitney_plot(
                data=df,
                condition=self.benefit_col,
                distribution=plot_col)

        return results

    def plot_survival(self, on, col=None, col_equals=None, how="os", threshold=None):
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
        plot_col, df = self.plot_init(on, col, col_equals)
        if df[plot_col].dtype == "bool":
            default_threshold = None
        else:
            default_threshold = "median"
        results = plot_kmf(
            df=df,
            condition_col=plot_col,
            xlabel='Overall Survival' if how == "os" else 'Progression-Free Survival',
            censor_col=self.deceased_col if how == "os" else self.progressed_or_deceased_col,
            survival_col=self.os_col if how == "os" else self.pfs_col,
            threshold=threshold if threshold is not None else default_threshold)
        print(results)

def col_func(cohort, on, col, col_equals):
    df = cohort.clinical_dataframe.copy()
    df[on] = df[col] == col_equals
    return on, df
