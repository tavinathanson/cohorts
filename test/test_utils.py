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

import pandas as pd
from cohorts.utils import strip_column_names, _strip_column_name

from nose.tools import eq_, ok_
from .test_basic import make_simple_cohort

def test_strip_single_column_name():
    res = [_strip_column_name(col) for col in ['PD-L1', 'PD L1', 'PD L1_']]
    eq_(res, ['pd_l1', 'pd_l1', 'pd_l1'])


def test_strip_column_names():
    d = {'one': pd.Series([1., 2., 3.], index=['a', 'b', 'c']),
         'two': pd.Series([1., 2., 3., 4.], index=['a', 'b', 'c', 'd']),
         'PD L1 (val)': pd.Series([1., 2., 3., 4.], index=['a', 'b', 'c', 'd']),
         'PD L1 (>1)': pd.Series([0., 1., 1., 1.], index=['a', 'b', 'c', 'd']),
         }
    df = pd.DataFrame(d)

    # should not error & should rename columns
    df2 = df.rename(columns=strip_column_names(df.columns))
    ok_((df2.columns != df.columns).any())

    # should not rename columns -- should error
    df3 = df.rename(columns=strip_column_names(
        df.columns, keep_paren_contents=False))
    ok_((df3.columns == df.columns).all())

def test_as_dataframe():
    cohort = make_simple_cohort()
    df_hello = pd.DataFrame({'one': pd.Series([1., 2., 3.], index=['a', 'b', 'c']),
         'two': pd.Series([1., 2., 3., 4.], index=['a', 'b', 'c', 'd']),
         'PD L1 (val)': pd.Series([1., 2., 3., 4.], index=['a', 'b', 'c', 'd']),
         'PD L1 (>1)': pd.Series([0., 1., 1., 1.], index=['a', 'b', 'c', 'd']),
         'the_id': ['1', '5', '7'],
         })
    def load_df():
        return df_hello
    df_loader = DataFrameLoader("hello", load_df, join_on="the_id")
    cohort.df_loaders = [df_loader]

    # test that column names haven't changed
    df = cohort.as_dataframe(join_with="hello")
    # column names should match those in df_hello
    ok_((df.columns == df_hello.columns).all())

    # test behavior with rename_cols=True
    df2 = cohort.as_dataframe(rename_cols=True)
    eq_(df2.columns, ['one','two','pd_l1_val','pd_l1_gt_1','the_id'])

    # test behavior with keep_paren_contents=False
    df3 = cohort.as_dataframe(rename_cols=True, keep_paren_contents=False)
    # column names should match those in df_hello
    ok_((df3.columns == df_hello.columns).all())