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

from . import data_path, generated_data_path, DATA_DIR
from .data_generate import generate_vcfs

from cohorts import Cohort
from cohorts.count import *

import pandas as pd
from nose.tools import raises, eq_
from os import path
from shutil import rmtree

from .test_basic import make_simple_cohort

FILE_FORMAT_1 = "patient_format1_%s.vcf"
FILE_FORMAT_2 = "patient_format2_%s.vcf"
FILE_FORMAT_3 = "patient_format3_%s.vcf"

def make_cohort(file_formats):
    cohort = make_simple_cohort()
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
        count_col, df = snv_count(cohort)
        eq_(len(df), 3)
        eq_(list(df[count_col]), [3, 3, 6])

        count_col, df = missense_snv_count(cohort)
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
        vcf_dir, cohort = make_cohort([FILE_FORMAT_1, FILE_FORMAT_2, FILE_FORMAT_3])

        # [3, 3, 6] and [4, 1, 5] use the same template, resulting in a union of [4, 3, 6] unique variants
        # [5, 2, 3] uses a separate template, resulting in a union of [4, 3, 6] + [5, 2, 3] = [9, 5, 9] unique variants
        count_col, df = snv_count(cohort)
        eq_(len(df), 3)
        eq_(list(df[count_col]), [9, 5, 9])

        # For intersection, variants need to appear in *all*, here. None of them do.
        count_col, df = snv_count(cohort, merge_type="intersection")
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
        vcf_dir, cohort = make_cohort([FILE_FORMAT_1, FILE_FORMAT_2])

        # [3, 3, 6] and [4, 1, 5] use the same template, resulting in a union of [4, 3, 6] unique variants
        count_col, df = snv_count(cohort)
        eq_(len(df), 3)
        eq_(list(df[count_col]), [4, 3, 6])

        # For intersection, some variants do appear in both.
        count_col, df = snv_count(cohort, merge_type="intersection")
        eq_(len(df), 3)
        eq_(list(df[count_col]), [3, 1, 5])
    finally:
        if vcf_dir is not None and path.exists(vcf_dir):
            rmtree(vcf_dir)
        if cohort is not None:
            cohort.clear_caches()
