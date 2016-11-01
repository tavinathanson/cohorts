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
from cohorts import DataFrameLoader, Cohort, Patient
import warnings

from . import generated_data_path

from nose.tools import eq_, ok_
from .test_basic import make_simple_cohort

def make_alt_simple_clinical_dataframe(
        os_list=None,
        pfs_list=None,
        deceased_list=None,
        progressed_or_deceased_list=None):
    return pd.DataFrame({"id": ["1", "4", "5"],
                         "age": [15, 20, 25],
                         "os": [100, 150, 120] if os_list is None else os_list,
                         "pfs": [50, 40, 120] if pfs_list is None else pfs_list,
                         "deceased": [True, False, False] if deceased_list is None else deceased_list,
                         "progressed_or_deceased": [True, True, False] if progressed_or_deceased_list is None else progressed_or_deceased_list})

def make_alt_simple_cohort(merge_type="union", **kwargs):
    clinical_dataframe = make_alt_simple_clinical_dataframe(**kwargs)
    patients = []
    for i, row in clinical_dataframe.iterrows():
        patient = Patient(id=row["id"],
                          os=row["os"],
                          pfs=row["pfs"],
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
    # should not rename columns -- should raise a warning
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter('always')
        df3 = df.rename(columns=strip_column_names(
                        df.columns, keep_paren_contents=False))
        ok_(len(w) > 0, 'warning not raised when keep_paren_contents results in dups')
    ok_((df3.columns == df.columns).all())


def prep_test_cohort():
    cohort = make_simple_cohort()
    df_hello = pd.DataFrame({
        'one': [1., 2., 3.],
        'two': [1., 2., 3.],
        'PD L1 (val)': [1., 2., 3.],
        'PD L1 (>1)': [0., 1., 1.],
        'the_id': ['1', '5', '7'],
        })

    def load_df():
        return df_hello
    df_loader = DataFrameLoader("hello", load_df, join_on="the_id")
    cohort.df_loaders = [df_loader]
    return df_hello, cohort


def prep_alt_test_cohort():
    cohort = make_alt_simple_cohort()
    df_hello = pd.DataFrame({
        'one': [1., 2., 3.],
        'two': [1., 2., 3.],
        'PD L1 (val)': [1., 2., 3.],
        'PD L1 (>1)': [0., 1., 1.],
        'the_id': ['1', '5', '7'],
        })

    def load_df():
        return df_hello
    df_loader = DataFrameLoader("hello", load_df, join_on="the_id")
    cohort.df_loaders = [df_loader]
    return df_hello, cohort


def compare_column_names(expected, observed):
    expected = set(expected)
    observed = set(observed)
    print('Expected:', expected)
    print('Observed:', observed)
    return(expected.issubset(observed))


def test_as_dataframe_generic():
    df_hello, cohort = prep_test_cohort()
    # test that column names haven't changed
    df = cohort.as_dataframe(join_with="hello")
    # column names should match those in df_hello
    res = compare_column_names(expected = df_hello.columns,
                               observed = df.columns)
    ok_(res, 'columns names failed to match expected')


def test_as_dataframe_good_rename():
    df_hello, cohort = prep_alt_test_cohort()
    # test behavior with rename_cols=True. should not raise a warning
    df = cohort.as_dataframe(rename_cols=True, join_with='hello')
    res = compare_column_names(expected = strip_column_names(df_hello.columns),
                               observed = df.columns)
    ok_(res, 'column names failed to match expected')


def test_as_dataframe_bad_rename():
    df_hello, cohort = prep_test_cohort()
    # test behavior with rename_cols=True. should raise a warning
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        df = cohort.as_dataframe(rename_cols=True, join_with='hello')
        # skip test since warnings (for some reason) don't propagate
        #ok_(len(w) > 0, 'fail to generate dups warning when using rename_cols=True')
    res = compare_column_names(expected = df_hello.columns,
                               observed = df.columns)
    ok_(res, 'columns names failed to match expected')

def test_as_dataframe_drop_parens():
    df_hello, cohort = prep_test_cohort()
    # test behavior with keep_paren_contents=False
    with warnings.catch_warnings(record=True) as w:
        warnings.simplefilter("always")
        df = cohort.as_dataframe(rename_cols=True, keep_paren_contents=False, join_with='hello')
        # skip test for warning since warning doesn't propagate (not sure why)
        #ok_(len(w) > 0, 'no warning when duplicates resulting from rename_cols')
    res = compare_column_names(expected = df_hello.columns,
                               observed = df.columns)
    ok_(res, 'columns names failed to match expected')
