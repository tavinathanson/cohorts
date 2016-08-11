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

def first_not_none_param(params, default):
    """
    Given a list of `params`, use the first param in the list that is
    not None. If all are None, fall back to `default`.
    """
    for param in params:
        if param is not None:
            return param
    return default

def filter_not_null(df, col):
    original_len = len(df)
    df = df[df[col].notnull()]
    updated_len = len(df)
    if updated_len < original_len:
        print("Missing %s for %d patients: from %d to %d" % (col, original_len - updated_len, original_len, updated_len))
    return df

def _provenance_str(provenance):
    """Utility function used by compare_provenance to print diff
    """
    return ["%s==%s" % (key, value) for (key, value) in provenance]

def compare_provenance(
        this_provenance, other_provenance,
        left_outer_diff = "In current but not comparison",
        right_outer_diff = "In comparison but not current"):
    """Utility function to compare two abritrary provenance dicts
    returns number of discrepancies.

    Parameters
    ----------
    this_provenance: provenance dict (to be compared to "other_provenance")
    other_provenance: comparison provenance dict

    (optional)
    left_outer_diff: description/prefix used when printing items in this_provenance but not in other_provenance
    right_outer_diff: description/prefix used when printing items in other_provenance but not in this_provenance

    Returns
    -----------
    Number of discrepancies (0: None)
    """
    ## if either this or other items is null, return 0
    if (not this_provenance or not other_provenance):
        return 0

    this_items = set(this_provenance.items())
    other_items = set(other_provenance.items())

    # Two-way diff: are any modules introduced, and are any modules lost?
    new_diff = this_items.difference(other_items)
    old_diff = other_items.difference(this_items)
    warn_str = ""
    if len(new_diff) > 0:
        warn_str += "%s: %s" % (
            left_outer_diff,
            _provenance_str(new_diff))
    if len(old_diff) > 0:
        warn_str += "%s: %s" % (
            right_outer_diff,
            _provenance_str(old_diff))

    if len(warn_str) > 0:
        warnings.warn(warn_str, Warning)

    return(len(new_diff)+len(old_diff))

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
