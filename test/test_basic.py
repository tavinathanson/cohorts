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
from data_generate import generate_vcfs

from cohorts import Cohort
from cohorts.load import InvalidDataError

import pandas as pd 
from nose.tools import raises, eq_

def make_simple_cohort(os_list=None,
                       pfs_list=None,
                       dead_list=None,
                       progressed_or_dead_list=None):
    return Cohort(
        data_dir=DATA_DIR,
        cache_dir=generated_data_path("cache"),
        sample_ids=[1, 4, 5],
        clinical_dataframe=pd.DataFrame({"id": [1, 4, 5],
                                         "OS": [100, 150, 120] if os_list is None else os_list,
                                         "PFS": [50, 40, 120] if pfs_list is None else pfs_list,
                                         "dead": [True, False, False] if dead_list is None else dead_list,
                                         "progressed_or_dead": [True, True, False] if progressed_or_dead_list is None else progressed_or_dead_list}),
        clinical_dataframe_id_col="id",
        os_col="OS",
        pfs_col="PFS",
        dead_col="dead",
        progressed_or_dead_col="progressed_or_dead")

def test_pfs_equal_to_os():
    # Should not error
    make_simple_cohort(pfs_list=[100, 150, 120])

@raises(InvalidDataError)
def test_pfs_greater_than_os():
    make_simple_cohort(pfs_list=[120, 150, 120])

@raises(InvalidDataError)
def test_progressed_vs_pfs():
    make_simple_cohort(progressed_or_dead_list=[True, False, False])

def test_simple_cohort():
    cohort = make_simple_cohort()
    eq_(len(cohort.clinical_dataframe), 3)
