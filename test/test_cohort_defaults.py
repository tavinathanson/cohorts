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
from cohorts import Cohort
from cohorts.variant_filters import no_filter
from cohorts.functions import *

from .data_generate import generate_vcfs

from nose.tools import eq_, ok_
from mock import MagicMock
from os import path
from shutil import rmtree
from collections import defaultdict
import pandas as pd

from .test_count import make_cohort

FILE_FORMAT_1 = "patient_format1_%s.vcf"

def test_default_filter_fn():
    """
    Test that filter_fn falls back to the Cohort default but can be overridden.
    """
    vcf_dir, cohort = None, None
    try:
        def default_filter_fn(filterable_variant):
            return filterable_variant.variant.start != 53513530

        vcf_dir, cohort = make_cohort([FILE_FORMAT_1])

        df = cohort.as_dataframe(snv_count)
        eq_(len(df), 3)
        eq_(list(df["snv_count"]), [3, 3, 6])

        cohort.filter_fn = default_filter_fn
        df = cohort.as_dataframe(snv_count)
        eq_(len(df), 3)
        eq_(list(df["snv_count"]), [2, 2, 5])

        df = cohort.as_dataframe(snv_count, filter_fn=None)
        eq_(list(df["snv_count"]), [2, 2, 5])

        df = cohort.as_dataframe(snv_count, filter_fn=no_filter)
        eq_(list(df["snv_count"]), [3, 3, 6])

        def another_filter_fn(filterable_variant):
            return default_filter_fn(filterable_variant) and filterable_variant.variant.start != 49658590
        df = cohort.as_dataframe(snv_count, filter_fn=another_filter_fn)
        eq_(list(df["snv_count"]), [1, 1, 4])
    finally:
        if vcf_dir is not None and path.exists(vcf_dir):
            rmtree(vcf_dir)
        if cohort is not None:
            cohort.clear_caches()

def test_default_normalized_per_mb():
    """
    Test that normalized_per_mb falls back to the Cohort default but can be overridden.
    """
    vcf_dir, cohort = None, None
    try:
        vcf_dir, cohort = make_cohort([FILE_FORMAT_1])

        df = cohort.as_dataframe(snv_count)
        eq_(len(df), 3)
        eq_(list(df["snv_count"]), [3, 3, 6])

        def load_ensembl_coverage(self):
            return pd.DataFrame({"patient_id": ["1", "4", "5"],
                                 "Num Loci": [2000000, 5000000, 8000000],
                                 "MB": [2, 5, 8]})
        Cohort.load_ensembl_coverage = load_ensembl_coverage
        vcf_dir, cohort = make_cohort([FILE_FORMAT_1])
        cohort.normalized_per_mb = True

        df = cohort.as_dataframe(snv_count)
        eq_(len(df), 3)
        # 3 / 2, 3 / 5, 6 / 8
        eq_(list(df["snv_count"]), [1.5, 0.6, 0.75])

        df = cohort.as_dataframe(snv_count, normalized_per_mb=None)
        eq_(list(df["snv_count"]), [1.5, 0.6, 0.75])

        df = cohort.as_dataframe(snv_count, normalized_per_mb=False)
        eq_(list(df["snv_count"]), [3, 3, 6])
    finally:
        if vcf_dir is not None and path.exists(vcf_dir):
            rmtree(vcf_dir)
        if cohort is not None:
            cohort.clear_caches()

def test_default_filter_fn_multiple_functions():
    """
    Across multiple functions, test that filter_fn falls back to the Cohort default.
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
            count_col = count_func.__name__

            # No Cohort default
            cohort.filter_fn = None
            df_no_arg  = cohort.as_dataframe(count_func)
            df_default_arg = cohort.as_dataframe(count_func, filter_fn=default_filter_fn)
            df_none_arg = cohort.as_dataframe(count_func, filter_fn=None)
            df_no_filter_arg = cohort.as_dataframe(count_func, filter_fn=no_filter)
            ok_((df_no_arg[count_col] == df_none_arg[count_col]).all())
            ok_((df_default_arg[count_col] != df_no_arg[count_col]).any())
            ok_((df_no_filter_arg[count_col] == df_no_arg[count_col]).all())

            # With a Cohort default
            cohort.filter_fn = default_filter_fn
            df_no_arg  = cohort.as_dataframe(count_func)
            df_default_arg = cohort.as_dataframe(count_func, filter_fn=default_filter_fn)
            df_none_arg = cohort.as_dataframe(count_func, filter_fn=None)
            df_no_filter_arg = cohort.as_dataframe(count_func, filter_fn=no_filter)
            ok_((df_no_arg[count_col] == df_none_arg[count_col]).all())
            ok_((df_default_arg[count_col] == df_no_arg[count_col]).all())
            ok_((df_no_filter_arg[count_col] != df_no_arg[count_col]).any())
    finally:
        if vcf_dir is not None and path.exists(vcf_dir):
            rmtree(vcf_dir)
        if cohort is not None:
            cohort.clear_caches()
