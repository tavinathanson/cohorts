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
from cohorts.utils import clean_column_names, _clean_column_name

from nose.tools import raises, eq_, ok_

def test_clean_column_name():
    res = [_clean_column_name(col) for col in ['PD-L1','PD L1','PD L1_']]
    eq_(res, ['PD_L1','PD_L1','PD_L1'])

def test_column_names():
    d = {'one' : pd.Series([1., 2., 3.], index=['a', 'b', 'c']),
      'two' : pd.Series([1., 2., 3., 4.], index=['a', 'b', 'c', 'd']),
      'PD L1 (value)': pd.Series([1., 2., 3., 4.], index=['a', 'b', 'c', 'd']),
      'PD L1 (>1)': pd.Series([0., 1., 1., 0.], index=['a', 'b', 'c', 'd']),
      }
    d = pd.DataFrame(d)
    ## should not error & should rename columns
    d2 = d.rename(columns = clean_column_names(d.columns))
    ok_((d2.columns != d.columns).any())

    ## should not rename columns -- should error
    d3 = d.rename(columns = clean_column_names(d.columns, keep_parens = False))
    ok_((d3.columns == d.columns).all())


