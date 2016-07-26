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

import re
import warnings


def _strip_column_name(col_name, keep_paren_contents=True):
    """
    Utility script applying several regexs to a string.
    Intended to be used by `strip_column_names`.

    This function will:
        1. replace informative punctuation components with text
        2. (optionally) remove text within parentheses
        3. replace remaining punctuation/whitespace with _
        4. strip leading/trailing punctuation/whitespace

    Parameters
    ----------
    col_name (str): input character string
    keep_paren_contents (logical):
        controls behavior of within-paren elements of text
         - if True, (the default) all text within parens retained
         - if False, text within parens will be removed from the field name

    Returns
    --------
    modified string for new field name

    Examples
    --------
    > print([_strip_column_name(col) for col in ['PD-L1','PD L1','PD L1_']])
    """
    # start with input
    new_col_name = col_name
    # replace meaningful punctuation with text equivalents
    # surround each with whitespace to enforce consistent use of _
    punctuation_to_text = {
        '<=': 'le',
        '>=': 'ge',
        '=<': 'le',
        '=>': 'ge',
        '<': 'lt',
        '>': 'gt',
        '#': 'num'
    }
    for punctuation, punctuation_text in punctuation_to_text.items():
        new_col_name = new_col_name.replace(punctuation, punctuation_text)

    # remove contents within ()
    if not(keep_paren_contents):
        new_col_name = re.sub('\([^)]*\)', '', new_col_name)

    # replace remaining punctuation/whitespace with _
    punct_pattern = '[\W_]+'
    punct_replacement = '_'
    new_col_name = re.sub(punct_pattern, punct_replacement, new_col_name)

    # remove leading/trailing _ if it exists (if last char was punctuation)
    new_col_name = new_col_name.strip("_")

    # TODO: check for empty string

    # return lower-case version of column name
    return new_col_name.lower()


def strip_column_names(cols, keep_paren_contents=True):
    """
    Utility script for renaming pandas columns to patsy-friendly names.

    Revised names have been:
        - stripped of all punctuation and whitespace (converted to text or `_`)
        - converted to lower case

    Takes a list of column names, returns a dict mapping
    names to revised names.

    If there are any concerns with the conversion, this will
    print a warning & return original column names.

    Parameters
    ----------

    cols (list): list of strings containing column names
    keep_paren_contents (logical):
        controls behavior of within-paren elements of text
         - if True, (the default) all text within parens retained
         - if False, text within parens will be removed from the field name

    Returns
    -------

    dict mapping col_names -> new_col_names

    Example
    -------

    > df = {'one' : pd.Series([1., 2., 3.], index=['a', 'b', 'c']),
      'two' : pd.Series([1., 2., 3., 4.], index=['a', 'b', 'c', 'd']),
      'PD L1 (value)': pd.Series([1., 2., 3., 4.], index=['a', 'b', 'c', 'd']),
      'PD L1 (>1)': pd.Series([0., 1., 1., 0.], index=['a', 'b', 'c', 'd']),
      }
    > df = pd.DataFrame(df)
    > df = df.rename(columns = strip_column_names(df.columns))

    ## observe, by comparison
    > df2 = df.rename(columns = strip_column_names(df.columns,
        keep_paren_contents=False))
    """

    # strip/replace punctuation
    new_cols = [
        _strip_column_name(col, keep_paren_contents=keep_paren_contents)
        for col in cols]

    if len(new_cols) != len(set(new_cols)):
        warn_str = 'Warning: strip_column_names (if run) would introduce duplicate names.'
        warn_str += ' Reverting column names to the original.'

        warnings.warn(warn_str, Warning)
        print('Warning: strip_column_names would introduce duplicate names. Please fix & try again.')
        return dict(zip(cols, cols))

    return dict(zip(cols, new_cols))
