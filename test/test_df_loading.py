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

from cohorts import Cohort, DataFrameLoader

import pandas as pd
from nose.tools import eq_

from .test_basic import make_simple_cohort

def test_df_loading():
    cohort = make_simple_cohort()
    df_hello = pd.DataFrame({"the_id": ["1", "5", "7"],
                             "hello_value": ["hello", "goodbye", "hello"]})
    def load_df():
        return df_hello
    df_loader = DataFrameLoader("hello", load_df, join_on="the_id")
    cohort.df_loaders = [df_loader]

    df = cohort.as_dataframe(join_with="hello")
    eq_(len(df), 2)
    # pylint: disable=no-member
    # pylint gets confused by as_dataframe's return type
    eq_(set(df.patient_id), set(["1", "5"]))
