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
import six.moves.cPickle as pickle
from types import FunctionType

import varcode
from varcode import VariantCollection, EffectCollection
from mhctools import NetMHCcons
from topiary import predict_epitopes_from_variants, epitopes_to_dataframe

from .survival import plot_kmf
from .plot import mann_whitney_plot, fishers_exact_plot

class Cohort(object):
    """Represents a cohort of patients."""

    def __init__(self,
                 data_dir,
                 cache_dir,
                 sample_ids,
                 normal_bam_ids,
                 tumor_bam_ids,
                 clinical_dataframe=None,
                 clinical_dataframe_id_col=None,
                 benefit_col=None,
                 os_col=None,
                 pfs_col=None,
                 dead_col=None,
                 progressed_col=None,
                 progressed_or_dead_col=None,
                 hla_alleles=None,
                 cache_results=True,
                 snv_file_format_funcs=None,
                 indel_file_format_funcs=None):
        self.data_dir = data_dir
        self.cache_dir = cache_dir
        self.normal_bam_ids = normal_bam_ids
        self.tumor_bam_ids = tumor_bam_ids
        self.clinical_dataframe = clinical_dataframe.copy()
        self.clinical_dataframe_id_col = clinical_dataframe_id_col
        self.benefit_col = benefit_col
        self.os_col = os_col
        self.pfs_col = pfs_col
        self.dead_col = dead_col
        self.progressed_col = progressed_col
        self.progressed_or_dead_col = progressed_or_dead_col
        self.hla_alleles = hla_alleles
        self.cache_results = cache_results
        self.snv_file_format_funcs = snv_file_format_funcs
        self.indel_file_format_funcs = indel_file_format_funcs
        self.sample_ids = sample_ids

        self.verify_survival()
        self.add_progressed_or_dead_col()

        variant_type_to_format_funcs = {}
        if self.snv_file_format_funcs is not None:
            variant_type_to_format_funcs["snv"] = self.snv_file_format_funcs
        if self.indel_file_format_funcs is not None:
            variant_type_to_format_funcs["indel"] = self.indel_file_format_funcs
        self.variant_type_to_format_funcs = variant_type_to_format_funcs

        self.variant_cache_name = "cached-variants"
        self.effect_cache_name = "cached-effects"
        self.nonsynonymous_effect_cache_name = "cached-nonsynonymous-effects"
        self.neoantigen_cache_name = "cached-neoantigens"

    def verify_survival(self):
        assert (self.clinical_dataframe[self.pfs_col] <=
                self.clinical_dataframe[self.os_col]).all(), (
                    "PFS should be <= OS, but PFS is larger than OS for some patients.")

        def func(row):
            if row[self.pfs_col] < row[self.os_col]:
                assert row[self.progressed_col], (
                    "A patient did not progress despite PFS being less than OS. "
                    "Full row: %s" % row)
        self.clinical_dataframe.apply(func, axis=1)

    def add_progressed_or_dead_col(self):
        if self.progressed_or_dead_col is None:
            self.progressed_or_dead_col = "progressed_or_dead"
            self.clinical_dataframe[self.progressed_or_dead_col] = (
                self.clinical_dataframe[self.progressed_col] | self.clinical_dataframe[self.dead_col])
        else:
            assert self.clinical_dataframe[self.progressed_or_dead_col].equals(
                self.clinical_dataframe[self.progressed_col] | self.clinical_dataframe[self.dead_col]), (
                    "progressed_or_dead_col should equal progressed_col || dead_col")

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

        for i, sample_id in enumerate(self.sample_ids):
            try:
                variants = self._load_single_sample_variants(
                    i, self.variant_type_to_format_funcs[variant_type], variant_type, merge_type) 
            except IOError:
                print("Variants did not exist for %s" % sample_id)
                continue

            sample_variants[sample_id] = variants
        return sample_variants

    def _load_single_sample_variants(self, sample_idx, file_format_funcs, variant_type, merge_type):
        sample_id = self.sample_ids[sample_idx]
        normal_bam_id = self.normal_bam_ids[sample_idx]
        tumor_bam_id = self.tumor_bam_ids[sample_idx]

        cached_file_name = "%s-%s-variants.pkl" % (variant_type, merge_type)
        cached = self.load_from_cache(self.variant_cache_name, sample_id, cached_file_name)
        if cached is not None:
            return cached

        combined_variants = []
        for file_format_func in file_format_funcs:
            file_name = file_format_func(
                sample_id, normal_bam_id, tumor_bam_id)
            variants = varcode.load_vcf_fast(path.join(self.data_dir, file_name))
            combined_variants.append(set(variants.elements))

        if len(combined_variants) == 1:
            assert merge_type == None, "Cannot specify a merge type when there is nothing to merge"
            merged_variants =  VariantCollection(combined_variants)
        else:
            assert merge_type in ["union", "intersection"], "Unknown merge type: %s" % merge_type
            if merge_type == "union":
                merged_variants = VariantCollection(set.union(*combined_variants))
            elif merge_type == "intersection":
                merged_variants = VariantCollection(set.intersection(*combined_variants))

        self.save_to_cache(merged_variants, self.variant_cache_name, sample_id, cached_file_name)

        return merged_variants

    def load_effects(self, only_nonsynonymous=False, variant_type="snv", merge_type="union"):
        sample_effects = {}
        for i, sample_id in enumerate(self.sample_ids):
            effects = self._load_single_sample_effects(
                i, self.variant_type_to_format_funcs, only_nonsynonymous, variant_type, merge_type)
            sample_effects[sample_id] = effects
        return sample_effects

    def _load_single_sample_effects(self, sample_idx, file_format_funcs,
                                    only_nonsynonymous, variant_type, merge_type):
        sample_id = self.sample_ids[sample_idx]

        cached_file_name = "%s-%s-effects.pkl" % (variant_type, merge_type)
        if only_nonsynonymous:
            cached = self.load_from_cache(self.nonsynonymous_effect_cache_name, sample_id, cached_file_name)
        else:
            cached = self.load_from_cache(self.effect_cache_name, sample_id, cached_file_name)
        if cached is not None:
            return cached

        variants = self._load_single_sample_variants(
            sample_idx, self.variant_type_to_format_funcs[variant_type], variant_type, merge_type)
        effects = variants.effects()
        nonsynonymous_effects = EffectCollection(
            effects.drop_silent_and_noncoding().top_priority_effect_per_variant().values())

        self.save_to_cache(effects, self.effect_cache_name, sample_id, cached_file_name)
        self.save_to_cache(nonsynonymous_effects, self.nonsynonymous_effect_cache_name, sample_id, cached_file_name)

        if only_nonsynonymous:
            return nonsynonymous_effects
        return effects

    def load_neoantigens(self, variant_type="snv", merge_type="union",
                         epitope_lengths=[8, 9, 10, 11], ic50_cutoff=500,
                         process_limit=10, max_file_records=None):
        assert self.hla_alleles is not None, "Cannot predict neoantigens without HLA alleles"

        dfs = []
        for i, sample_id in enumerate(self.sample_ids):
            df_epitopes = self._load_single_sample_neoantigens(
                i, variant_type, merge_type, epitope_lengths, ic50_cutoff,
                process_limit, max_file_records)
            dfs.append(df_epitopes)
        return pd.concat(dfs)

    def _load_single_sample_neoantigens(self, sample_idx, variant_type, merge_type,
                                        epitope_lengths, ic50_cutoff, process_limit, max_file_records):
        sample_id = self.sample_ids[sample_idx]

        cached_file_name = "%s-%s-neoantigens.csv" % (variant_type, merge_type)
        cached = self.load_from_cache(self.neoantigen_cache_name, sample_id, cached_file_name)
        if cached is not None:
            return cached

        variants = self._load_single_sample_variants(
            sample_idx, self.variant_type_to_format_funcs[variant_type], variant_type, merge_type)
        hla_alleles = self.hla_alleles[sample_idx]
        mhc_model = NetMHCcons(
            alleles=hla_alleles,
            epitope_lengths=epitope_lengths,
            max_file_records=max_file_records,
            process_limit=process_limit)
        epitopes = predict_epitopes_from_variants(
            variants=variants,
            mhc_model=mhc_model,
            ic50_cutoff=ic50_cutoff,
            # Only include peptides with a variant
            only_novel_epitopes=True)
        df_epitopes = epitopes_to_dataframe(epitopes)
        df_epitopes["sample_id"] = sample_id

        self.save_to_cache(df_epitopes, self.neoantigen_cache_name, sample_id,  cached_file_name)

        return df_epitopes

    def clear_caches(self):
        self.clear_variant_cache()
        self.clear_neoantigen_cache()

    def clear_cache(self, cache_name):
        cache_path = path.join(self.cache_dir, cache_name)
        if path.exists(cache_path):
            rmtree(cache_path)

    def clear_variant_cache(self):
        self.clear_cache(self.variant_cache_name)

    def clear_neoantigen_cache(self):
        self.clear_cache(self.neoantigen_cache_name)

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
            If `on` is not specified, the column name on which to split the cohort
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
            If `on` is not specified, the column name on which to split the cohort
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
            censor_col=self.dead_col if how == "os" else self.progressed_or_dead_col,
            survival_col=self.os_col if how == "os" else self.pfs_col,
            threshold=threshold if threshold is not None else default_threshold)
        print(results)
              
def col_func(cohort, on, col, col_equals):
    df = cohort.clinical_dataframe.copy()
    df[on] = df[col] == col_equals
    return on, df
