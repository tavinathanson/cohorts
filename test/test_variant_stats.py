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
from nose.tools import ok_, raises
from shutil import rmtree
from os import path

from .data_generate import generate_vcfs
from .test_basic import make_simple_cohort
from cohorts.variant_stats import (strelka_somatic_variant_stats, 
                                   mutect_somatic_variant_stats,
                                   variant_stats_from_variant)

FILE_FORMAT_1 = "patient_format1_%s.mutect.vcf"
FILE_FORMAT_2 = "patient_format2_%s.strelka.vcf"
FILE_FORMAT_3 = "patient_format3_%s.vcf"

def make_cohort(file_formats):
    cohort = make_simple_cohort()
    patient_ids = [patient.id for patient in cohort]
    vcf_dir = generate_vcfs(id_to_mutation_count=dict(zip(patient_ids, [3, 3, 6])),
                            file_format=FILE_FORMAT_1,
                            template_name="vcf_template_1.vcf")
    _ = generate_vcfs(id_to_mutation_count=dict(zip(patient_ids, [5, 2, 3])),
                      file_format=FILE_FORMAT_2,
                      template_name="vcf_template_2.vcf")

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

def extraction(cohort, extractor):
    variants = cohort.load_variants()
    for (sample, sample_variants) in variants.items():
        ok_(len(sample_variants) > 0)
        for variant in sample_variants:
            somatic_stats = extractor(variant, sample_variants.metadata[variant])

            ok_(somatic_stats.tumor_stats.depth > 10)
            ok_(somatic_stats.normal_stats.depth > 10)

            ok_(somatic_stats.tumor_stats.variant_allele_frequency > 0.0)
            ok_(somatic_stats.tumor_stats.variant_allele_frequency < 0.5)

            ok_(somatic_stats.normal_stats.variant_allele_frequency == 0)

@raises(ValueError)
def test_extract_unsupported_stats():
    vcf_dir, cohort = None, None
    try:
        # Use Strelka and unsupported VCF format
        vcf_dir, cohort = make_cohort([FILE_FORMAT_2, FILE_FORMAT_3])
        extraction(cohort, extractor=variant_stats_from_variant)

    finally:
        if vcf_dir is not None and path.exists(vcf_dir):
            rmtree(vcf_dir)
        if cohort is not None:
            cohort.clear_caches()

def test_extract_strelka_stats():
    vcf_dir, cohort = None, None
    try:
        # Use Strelka VCF format
        vcf_dir, cohort = make_cohort([FILE_FORMAT_2])
        extraction(cohort, extractor=strelka_somatic_variant_stats)

    finally:
        if vcf_dir is not None and path.exists(vcf_dir):
            rmtree(vcf_dir)
        if cohort is not None:
            cohort.clear_caches()

def test_extract_mutect_stats():
    vcf_dir, cohort = None, None
    try:
        # Use Mutect VCF format
        vcf_dir, cohort = make_cohort([FILE_FORMAT_1])
        extraction(cohort, extractor=mutect_somatic_variant_stats)

    finally:
        if vcf_dir is not None and path.exists(vcf_dir):
            rmtree(vcf_dir)
        if cohort is not None:
            cohort.clear_caches()

def test_extract_strelka_mutect_stats():
    vcf_dir, cohort = None, None
    try:
        # Use Mutect and Strelka VCF format
        vcf_dir, cohort = make_cohort([FILE_FORMAT_1, FILE_FORMAT_2])
        extraction(cohort, extractor=variant_stats_from_variant)

    finally:
        if vcf_dir is not None and path.exists(vcf_dir):
            rmtree(vcf_dir)
        if cohort is not None:
            cohort.clear_caches()
