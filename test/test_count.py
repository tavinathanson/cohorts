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

from varcode import ExonicSpliceSite, Substitution
from cohorts.variant_filters import no_filter

from .data_generate import generate_vcfs
from .functions import *

from nose.tools import raises, eq_, ok_
from mock import MagicMock
from os import path
from shutil import rmtree
from collections import defaultdict
import pandas as pd

from .test_basic import make_simple_cohort

FILE_FORMAT_1 = "patient_format1_%s.vcf"
FILE_FORMAT_2 = "patient_format2_%s.vcf"
FILE_FORMAT_3 = "patient_format3_%s.vcf"

def make_cohort(file_formats, merge_type="union"):
    cohort = make_simple_cohort(merge_type=merge_type)
    patient_ids = [patient.id for patient in cohort]
    vcf_dir = generate_vcfs(id_to_mutation_count=dict(zip(patient_ids, [3, 3, 6])),
                            file_format=FILE_FORMAT_1,
                            template_name="vcf_template_1.vcf")
    _ = generate_vcfs(id_to_mutation_count=dict(zip(patient_ids, [4, 1, 5])),
                      file_format=FILE_FORMAT_2,
                      template_name="vcf_template_1.vcf")
    _ = generate_vcfs(id_to_mutation_count=dict(zip(patient_ids, [5, 2, 3])),
                      file_format=FILE_FORMAT_3,
                      template_name="vcf_template_2.vcf")
    for patient in cohort:
        vcf_paths = []
        for file_format in file_formats:
            vcf_filename = (file_format % patient.id)
            vcf_path = path.join(vcf_dir, vcf_filename)
            vcf_paths.append(vcf_path)
        patient.snv_vcf_paths = vcf_paths
    return vcf_dir, cohort

def test_snv_counts():
    """
    Generate VCFs per-sample, and confirm that the counting functions work as expected.
    """
    vcf_dir, cohort = None, None
    try:
        # Use all three VCF sources
        vcf_dir, cohort = make_cohort([FILE_FORMAT_1])

        # The SNV count should be exactly what we generated
        count_col, df = cohort.as_dataframe(snv_count)
        eq_(len(df), 3)
        eq_(list(df[count_col]), [3, 3, 6])

        count_col, df = cohort.as_dataframe(missense_snv_count)
        eq_(len(df), 3)
        eq_(list(df[count_col]), [2, 2, 4])
    finally:
        if vcf_dir is not None and path.exists(vcf_dir):
            rmtree(vcf_dir)
        if cohort is not None:
            cohort.clear_caches()

def test_merge_three():
    """
    Generate three VCFs per-sample and confirm that merging works as expected.
    """
    vcf_dir, cohort = None, None
    try:
        # Use all three VCF sources
        vcf_dir, cohort = make_cohort([FILE_FORMAT_1, FILE_FORMAT_2, FILE_FORMAT_3],
                                      merge_type="union")

        # [3, 3, 6] and [4, 1, 5] use the same template, resulting in a union of [4, 3, 6] unique variants
        # [5, 2, 3] uses a separate template, resulting in a union of [4, 3, 6] + [5, 2, 3] = [9, 5, 9] unique variants
        count_col, df = cohort.as_dataframe(snv_count)
        eq_(len(df), 3)
        eq_(list(df[count_col]), [9, 5, 9])

        # For intersection, variants need to appear in *all*, here. None of them do.
        vcf_dir, cohort = make_cohort([FILE_FORMAT_1, FILE_FORMAT_2, FILE_FORMAT_3],
                                      merge_type="intersection")
        count_col, df = cohort.as_dataframe(snv_count)
        eq_(list(df[count_col]), [0, 0, 0])
    finally:
        if vcf_dir is not None and path.exists(vcf_dir):
            rmtree(vcf_dir)
        if cohort is not None:
            cohort.clear_caches()

def test_merge_two():
    """
    Generate two VCFs per-sample and confirm that merging works as expected.
    """
    vcf_dir, cohort = None, None
    try:
        # Now, with only two VCF sources
        vcf_dir, cohort = make_cohort([FILE_FORMAT_1, FILE_FORMAT_2],
                                      merge_type="union")

        # [3, 3, 6] and [4, 1, 5] use the same template, resulting in a union of [4, 3, 6] unique variants
        count_col, df = cohort.as_dataframe(snv_count)
        eq_(len(df), 3)
        eq_(list(df[count_col]), [4, 3, 6])

        # For intersection, some variants do appear in both.
        vcf_dir, cohort = make_cohort([FILE_FORMAT_1, FILE_FORMAT_2],
                                      merge_type="intersection")
        count_col, df = cohort.as_dataframe(snv_count)
        eq_(len(df), 3)
        eq_(list(df[count_col]), [3, 1, 5])

        cohort_variants = cohort.load_variants(filter_fn=None)
        for (sample, variants) in cohort_variants.items():
            for variant in variants:
                metadata = variants.metadata[variant]
                eq_(len(metadata), 2) # Each variant has two metadata entries

    finally:
        if vcf_dir is not None and path.exists(vcf_dir):
            rmtree(vcf_dir)
        if cohort is not None:
            cohort.clear_caches()

def test_filter_variants():
    vcf_dir, cohort = None, None
    try:
        vcf_dir, cohort = make_cohort([FILE_FORMAT_1, FILE_FORMAT_2])

        def filter_g_variants(filterable_variant):
            return filterable_variant.variant.ref == 'G'
        g_variants = {'1': 2, '4': 1, '5': 3}

        cohort_variants = cohort.load_variants(filter_fn=filter_g_variants)

        for (sample, variants) in cohort_variants.items():
            eq_(len(variants), g_variants[sample])

    finally:
        if vcf_dir is not None and path.exists(vcf_dir):
            rmtree(vcf_dir)
        if cohort is not None:
            cohort.clear_caches()

def test_filter_effects():
    vcf_dir, cohort = None, None
    try:
        vcf_dir, cohort = make_cohort([FILE_FORMAT_1])

        def filter_substitution_effects(filterable_effect):
            return type(filterable_effect.effect) == Substitution
        missense_counts = {'1': 2, '4': 2, '5': 4}

        cohort_effects = cohort.load_effects(only_nonsynonymous=True, filter_fn=filter_substitution_effects)
        for (sample, effects) in cohort_effects.items():
            eq_(len(effects), missense_counts[sample])

        def filter_exonic_splice_site_effects(filterable_effect):
            return type(filterable_effect.effect) == ExonicSpliceSite
        splice_site_counts = {'1': 1, '4': 1, '5': 2}

        cohort_effects = cohort.load_effects(only_nonsynonymous=True, filter_fn=filter_exonic_splice_site_effects)
        for (sample, effects) in cohort_effects.items():
            eq_(len(effects), splice_site_counts[sample])

    finally:
        if vcf_dir is not None and path.exists(vcf_dir):
            rmtree(vcf_dir)
        if cohort is not None:
            cohort.clear_caches()

def test_multiple_effects():
    vcf_dir, cohort = None, None
    try:
        vcf_dir, cohort = make_cohort([FILE_FORMAT_1])
        effects = cohort.load_effects(only_nonsynonymous=False)
        for patient in cohort:
            variant_to_effect_set = defaultdict(set)
            for effect in effects[patient.id]:
                variant_to_effect_set[effect.variant].add(effect)
            for variant, effect_set in variant_to_effect_set.items():
                eq_(len(effect_set), 1,
                    "Variant %s should only have 1 effect but it has %s"
                    % (variant, len(effect_set)))
    finally:
        if vcf_dir is not None and path.exists(vcf_dir):
            rmtree(vcf_dir)
        if cohort is not None:
            cohort.clear_caches()

def test_cohort_default():
    """
    Checks that the Cohort default is used when intended, and not used when overridden.
    """
    vcf_dir, cohort = None, None
    try:
        vcf_dir, cohort = make_cohort([FILE_FORMAT_1])

        def default_filter_fn(filterable_variant):
            return filterable_variant.variant.ref == "A"

        # Mock out load_single_patient_isovar since we don't actually have that data,
        # yet it's needed for expressed_missense_snv_count.
        variants = cohort.load_variants()
        isovar_mocked = {}
        for patient in cohort:
            chr_list = []
            pos_list = []
            ref_list = []
            alt_list = []
            patient_variants = variants[patient.id]
            for variant in patient_variants:
                chr_list.append(variant.contig)
                pos_list.append(variant.start)
                ref_list.append(variant.ref)
                alt_list.append(variant.alt)
            df_isovar_mocked = pd.DataFrame({"chr": chr_list,
                                             "pos": pos_list,
                                             "ref": ref_list,
                                             "alt": alt_list})
            isovar_mocked[patient.id] = df_isovar_mocked

        def mock_function(patient, **kwargs):
            return isovar_mocked[patient.id]

        cohort.load_single_patient_isovar = MagicMock(side_effect=mock_function)

        # Especially need to include expressed_missense_snv_count, which works a little differently.
        for count_func in [snv_count, missense_snv_count, expressed_missense_snv_count]:
            # No Cohort default
            cohort.filter_fn = None
            count_col_no_arg, df_no_arg  = cohort.as_dataframe(count_func)
            count_col_default_arg, df_default_arg = cohort.as_dataframe(count_func, filter_fn=default_filter_fn)
            count_col_none_arg, df_none_arg = cohort.as_dataframe(count_func, filter_fn=None)
            count_col_no_filter_arg, df_no_filter_arg = cohort.as_dataframe(count_func, filter_fn=no_filter)
            ok_((df_no_arg[count_col_no_arg] == df_none_arg[count_col_none_arg]).all())
            ok_((df_default_arg[count_col_default_arg] != df_no_arg[count_col_no_arg]).any())
            ok_((df_no_filter_arg[count_col_no_filter_arg] == df_no_arg[count_col_no_arg]).all())

            # With a Cohort default
            cohort.filter_fn = default_filter_fn
            count_col_no_arg, df_no_arg  = cohort.as_dataframe(count_func)
            count_col_default_arg, df_default_arg = cohort.as_dataframe(count_func, filter_fn=default_filter_fn)
            count_col_none_arg, df_none_arg = cohort.as_dataframe(count_func, filter_fn=None)
            count_col_no_filter_arg, df_no_filter_arg = cohort.as_dataframe(count_func, filter_fn=no_filter)
            ok_((df_no_arg[count_col_no_arg] == df_none_arg[count_col_none_arg]).all())
            ok_((df_default_arg[count_col_default_arg] == df_no_arg[count_col_no_arg]).all())
            ok_((df_no_filter_arg[count_col_no_filter_arg] != df_no_arg[count_col_no_arg]).any())
    finally:
        if vcf_dir is not None and path.exists(vcf_dir):
            rmtree(vcf_dir)
        if cohort is not None:
            cohort.clear_caches()
