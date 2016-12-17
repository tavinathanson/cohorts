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

from varcode.effects.effect_classes import ExonicSpliceSite, Substitution
from varcode import Variant, VariantCollection
from cohorts.variant_filters import no_filter
from cohorts.functions import *

from .data_generate import generate_vcfs
from . import data_path

from nose.tools import eq_, ok_
from mock import MagicMock
from os import path
from shutil import rmtree
from collections import defaultdict
import pandas as pd

from .test_basic import make_simple_cohort

FILE_FORMAT_1 = "patient_format1_%s.vcf"
FILE_FORMAT_2 = "patient_format2_%s.vcf"
FILE_FORMAT_3 = "patient_format3_%s.vcf"

MAF_FILE = "test_maf.maf"

def make_cohort(file_formats, merge_type="union", **kwargs):
    cohort = make_simple_cohort(merge_type=merge_type, **kwargs)
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
        patient.variants = vcf_paths
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
        df = cohort.as_dataframe(snv_count)
        eq_(len(df), 3)
        eq_(list(df["snv_count"]), [3, 3, 6])

        df = cohort.as_dataframe(missense_snv_count)
        eq_(len(df), 3)
        eq_(list(df["missense_snv_count"]), [2, 2, 4])
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
        df = cohort.as_dataframe(snv_count)
        eq_(len(df), 3)
        eq_(list(df["snv_count"]), [9, 5, 9])

        # For intersection, variants need to appear in *all*, here. None of them do.
        vcf_dir, cohort = make_cohort([FILE_FORMAT_1, FILE_FORMAT_2, FILE_FORMAT_3],
                                      merge_type="intersection")
        df = cohort.as_dataframe(snv_count)
        eq_(list(df["snv_count"]), [0, 0, 0])
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
        df = cohort.as_dataframe(snv_count)
        eq_(len(df), 3)
        eq_(list(df["snv_count"]), [4, 3, 6])

        # For intersection, some variants do appear in both.
        vcf_dir, cohort = make_cohort([FILE_FORMAT_1, FILE_FORMAT_2],
                                      merge_type="intersection")
        df = cohort.as_dataframe(snv_count)
        eq_(len(df), 3)
        eq_(list(df["snv_count"]), [3, 1, 5])

        cohort_variants = cohort.load_variants(filter_fn=None)
        for (sample, variants) in cohort_variants.items():
            for variant in variants:
                sources = variants.sources
                eq_(len(sources), 2) # Each variant has two sources entries

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

def test_multiple_variant_forms():
    """
    Load VCF, MAF and VariantCollection together.
    """
    vcf_dir, cohort = None, None
    try:
        vcf_dir, cohort = make_cohort([FILE_FORMAT_1])
        patient = cohort[0]
        patient.variants.append(data_path(MAF_FILE))
        # Make sure listing the file twice has no effect.
        patient.variants.append(data_path(MAF_FILE))
        variant = Variant(start=1000000, ref="A", alt="T", contig=1, ensembl=75)
        patient.variants.append(VariantCollection([variant]))

        cohort_variants = cohort.load_variants(patients=[patient])

        # Make sure the VariantCollection was included.
        eq_(len(cohort_variants[patient.id].filter(lambda v: v.start == 1000000)), 1)

        # Make sure the VCF was included.
        eq_(len(cohort_variants[patient.id].filter(lambda v: v.start == 53513530)), 1)

        # Make sure the MAF was included.
        eq_(len(cohort_variants[patient.id].filter(lambda v: v.start == 1650797)), 1)

        # Make sure a non-existant variant is not included.
        eq_(len(cohort_variants[patient.id].filter(lambda v: v.start == 1650798)), 0)
    finally:
        if vcf_dir is not None and path.exists(vcf_dir):
            rmtree(vcf_dir)
        if cohort is not None:
            cohort.clear_caches()
