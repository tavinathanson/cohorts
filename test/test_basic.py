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
from .functions import *

from cohorts import Cohort, Patient
from cohorts.load import InvalidDataError

import pandas as pd
from nose.tools import raises, eq_, ok_

def make_simple_clinical_dataframe(
        os_list=None,
        pfs_list=None,
        deceased_list=None,
        progressed_or_deceased_list=None):
    return pd.DataFrame({"id": ["1", "4", "5"],
                         "age": [15, 20, 25],
                         "OS": [100, 150, 120] if os_list is None else os_list,
                         "PFS": [50, 40, 120] if pfs_list is None else pfs_list,
                         "deceased": [True, False, False] if deceased_list is None else deceased_list,
                         "progressed_or_deceased": [True, True, False] if progressed_or_deceased_list is None else progressed_or_deceased_list})

def make_simple_cohort(merge_type="union", **kwargs):
    clinical_dataframe = make_simple_clinical_dataframe(**kwargs)
    patients = []
    for i, row in clinical_dataframe.iterrows():
        patient = Patient(id=row["id"],
                          os=row["OS"],
                          pfs=row["PFS"],
                          deceased=row["deceased"],
                          progressed_or_deceased=row["progressed_or_deceased"],
                          additional_data=row
                          )
        patients.append(patient)

    return Cohort(
        patients=patients,
        responder_pfs_equals_os=True,
        merge_type=merge_type,
        cache_dir=generated_data_path("cache"))

def test_pfs_equal_to_os():
    # Should not error
    make_simple_cohort(pfs_list=[100, 150, 120])

@raises(InvalidDataError)
def test_pfs_greater_than_os():
    make_simple_cohort(pfs_list=[120, 150, 120])

@raises(InvalidDataError)
def test_progressed_vs_pfs():
    make_simple_cohort(progressed_or_deceased_list=[True, False, False])

def test_simple_cohort():
    cohort = make_simple_cohort()
    eq_(len(cohort.as_dataframe()), 3)

    columns = set(cohort.as_dataframe().columns)
    ok_("id" in columns)
    ok_("patient_id" in columns)
    ok_("age" in columns)
    ok_("pfs" in columns)
    ok_("os" in columns)
