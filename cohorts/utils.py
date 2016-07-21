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

def _clean_column_name(x, keep_parens = True):
    """
        Utility script applying several regexs to a string.  Intended to be used by `clean_column_names`.

        This function will: 
            1. replace informative punctuation components with text
            2. (optionally) remove text within parentheses
            3. replace remaining punctuation/whitespace with _
            4. remove trailing punctuation/whitespace

        Parameters 
        -----------
        str (str): input character string
        keep_parens (logical): control behavior of within-paren elements of text
             - if True, (the default) all text within parens retained
             - if False, text within parens will be removed from the field name

        Returns
        --------
        modified string for new field name 


        Examples
        --------- 
        > print([clean_column_name(col) for col in ['PD-L1','PD L1','PD L1_']])

    """
    ## start with input
    newx = x
    
    ## replace meaningful punctuation with text equivalents
    ## surround each with whitespace to enforce consistent use of _
    newx = re.sub('<=', ' le ', newx)
    newx = re.sub('=<', ' le ', newx)
    newx = re.sub('>=', ' ge ', newx)
    newx = re.sub('=>', ' ge ', newx)
    newx = re.sub('<', ' lt ', newx)
    newx = re.sub('>', ' gt ', newx)
    newx = re.sub('#', ' num ', newx)
    
    ## remove contents within ()
    if not(keep_parens):
        newx = re.sub('\([^)]*\)', '', newx)
    
    ## replace remaining punctuation/whitespace with _
    punct_pattern = '[\W_]+'
    punct_replacement = '_'
    newx = re.sub(punct_pattern, punct_replacement, newx)
    
    ## remove trailing _ if it exists (if last char was punctuation)
    newx = re.sub('_$', '', newx)
    return newx

def clean_column_names(cols, keep_parens = True):
    """
        Utility script for renaming pandas columns to patsy-friendly names.
        
        Revised names have:
            - punctuation and whitespace -> text or _
            - all text in lower case 

        Takes a list of column names, returns a dictionary mapping names to revised names.

        If there are any concerns with the conversion, this will print a warning & return original column names. 
        
        Parameters 
        --------------
            cols (list): list of strings containing column names
            keep_parens (logical): whether to retain text within parentheses (defalt: True)

        Returns
        --------
            dict mapping cols -> newcols

        Example 
        -------- 
        > d = {'one' : pd.Series([1., 2., 3.], index=['a', 'b', 'c']),
              'two' : pd.Series([1., 2., 3., 4.], index=['a', 'b', 'c', 'd']),
              'PD L1 (value)': pd.Series([1., 2., 3., 4.], index=['a', 'b', 'c', 'd']),
              'PD L1 (>1)': pd.Series([0., 1., 1., 0.], index=['a', 'b', 'c', 'd']),
              }
        > d = pd.DataFrame(d)
        > d = d.rename(columns = clean_column_names(d.columns))
        
        ## observe, by comparison
        > d2 = d.rename(columns = clean_column_names(d.columns, keep_parens = False))

    """
    ## clean/replace punctuation
    newcols = [_clean_column_name(col, keep_parens = keep_parens) for col in cols]
    
    ## confirm no duplicates in remaining columns
    if not(len(newcols) == len(cols)):
        print('Warning: length of newcols not the same as cols')
        return dict(zip(cols, cols))
    
    if not(len(newcols) == len(set(newcols))):
        print('Warning: cleanup introduces duplicate names. Please resolve & try again.')
        return dict(zip(cols, cols))
    
    ## modify all column names to lower case
    newcols_lc = [col.lower() for col in newcols]
    if not(len(newcols_lc) == len(set(newcols_lc))):
        print('Warning: changing case introduces duplicate names. Skipping this.')
    else:
        newcols = newcols_lc
    
    return dict(zip(cols, newcols))