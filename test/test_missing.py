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

from . import data_path, generated_data_path
from .data_generate import generate_vcfs

from cohorts.cohort import filter_not_null
from cohorts.functions import *

import pandas as pd
import numpy as np
from nose.tools import raises, eq_, ok_
from os import path, makedirs
from shutil import rmtree

from .test_basic import make_simple_cohort

def make_missing_vcf_cohort(patient_ids_with_missing_paths, missing_paths):
    file_format = "patient_format_%s.vcf"
    cohort = make_simple_cohort()
    patient_ids = [patient.id for patient in cohort]
    vcf_dir = generate_vcfs(id_to_mutation_count=dict(zip(patient_ids, [3, 3, 6])),
                            file_format=file_format,
                            template_name="vcf_template_1.vcf")
    for patient in cohort:
        vcf_paths = []
        vcf_filename = (file_format % patient.id)
        vcf_path = path.join(vcf_dir, vcf_filename)
        vcf_paths.append(vcf_path)
        if patient.id in patient_ids_with_missing_paths:
            patient.variants = missing_paths
        else:
            patient.variants = vcf_paths
        patient.hla_alleles = ["HLA-A02:01"]
    return vcf_dir, cohort

def test_broken_vcf_path():
    """
    Generate VCFs per-sample, and confirm that the NaNs are returned when
    VCF paths are broken.
    """
    vcf_dir, cohort = None, None
    try:
        vcf_dir, cohort = make_missing_vcf_cohort(
            patient_ids_with_missing_paths=["1"],
            missing_paths=["nonexistant_path.vcf"])

        df = cohort.as_dataframe(snv_count)
        eq_(len(df), 3)
        ok_(np.isnan(list(df["snv_count"])[0]))

        df = cohort.as_dataframe(missense_snv_count)
        eq_(len(df), 3)
        ok_(np.isnan(list(df["missense_snv_count"])[0]))
    finally:
        if vcf_dir is not None and path.exists(vcf_dir):
            rmtree(vcf_dir)
        if cohort is not None:
            cohort.clear_caches()

def test_missing_vcf_path():
    """
    Generate VCFs per-sample, and confirm that the NaNs are returned when
    VCF paths are broken.
    """
    vcf_dir, cohort = None, None
    try:
        vcf_dir, cohort = make_missing_vcf_cohort(
            patient_ids_with_missing_paths=["1"],
            missing_paths=[])

        df = cohort.as_dataframe(snv_count)
        eq_(len(df), 3)
        ok_(np.isnan(list(df["snv_count"])[0]))

        df = cohort.as_dataframe(missense_snv_count)
        eq_(len(df), 3)
        ok_(np.isnan(list(df["missense_snv_count"])[0]))
    finally:
        if vcf_dir is not None and path.exists(vcf_dir):
            rmtree(vcf_dir)
        if cohort is not None:
            cohort.clear_caches()

def test_no_neoantigens():
    """
    Make sure that exactly 0 neoantigens is treated differently than no calculated
    neoantigens.
    """
    vcf_dir, cohort = None, None
    try:
        vcf_dir, cohort = make_missing_vcf_cohort(
            patient_ids_with_missing_paths=["1"],
            missing_paths=[])

        generate_empty_neoantigens(cohort,
                                   patient_ids_with_zero_neoantigens=["4", "5"],
                                   patient_ids_with_empty_neoantigens=["1"])

        df = cohort.as_dataframe(neoantigen_count)
        eq_(len(df), 3)

        # This is used internally by plotting functions.
        # Ensure it filters NaN but not 0.
        eq_(len(filter_not_null(df, "neoantigen_count")), 2)
    finally:
        if vcf_dir is not None and path.exists(vcf_dir):
            rmtree(vcf_dir)
        if cohort is not None:
            cohort.clear_caches()

def generate_empty_neoantigens(cohort,
                               patient_ids_with_zero_neoantigens,
                               patient_ids_with_empty_neoantigens):
    for patient in cohort:
        if patient.id in patient_ids_with_empty_neoantigens:
            continue
        elif patient.id in patient_ids_with_zero_neoantigens:
            neoantigen_path = generated_data_path(
                path.join("cache", "cached-neoantigens",
                          patient.id, "union-neoantigens.csv"))
            if not path.exists(path.dirname(neoantigen_path)):
                makedirs(path.dirname(neoantigen_path))
                with open(neoantigen_path, "w") as f:
                    df_neoantigens = pd.read_csv(data_path("empty-neoantigens.csv"))
                    # pylint: disable=no-member
                    # pylint gets confused by to_csv
                    df_neoantigens.to_csv(neoantigen_path)
        else:
            raise ValueError("Patient ID needs to be empty or zero")
