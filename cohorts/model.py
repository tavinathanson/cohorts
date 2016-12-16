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
from types import FunctionType
import numpy as np
import pandas as pd
from sklearn.utils import resample
from sklearn.metrics import roc_auc_score
import lifelines as ll
import patsy

from .utils import get_logger

logger = get_logger(__name__)

def is_single_class(arr, col):
    if type(arr) == pd.DataFrame:
        arr_positive_indices = list(arr[arr[col] == 1][col])
        arr_negative_indices = list(arr[arr[col] == 0][col])
    else:
        arr_positive_indices = (arr == 1).nonzero()[0]
        arr_negative_indices = (arr == 0).nonzero()[0]
    return len(arr_positive_indices) == 0 or len(arr_negative_indices) == 0

def bootstrap_auc(df, col, pred_col, n_bootstrap=1000):
    """
    Calculate the boostrapped AUC for a given col trying to predict a pred_col.

    Parameters
    ----------
    df : pandas.DataFrame
    col : str
        column to retrieve the values from
    pred_col : str
        the column we're trying to predict
    n_boostrap : int
        the number of bootstrap samples

    Returns
    -------
    list : AUCs for each sampling
    """
    scores = np.zeros(n_bootstrap)
    old_len = len(df)
    df.dropna(subset=[col], inplace=True)
    new_len = len(df)
    if new_len < old_len:
        logger.info("Dropping NaN values in %s to go from %d to %d rows" % (col, old_len, new_len))
    preds = df[pred_col].astype(int)
    for i in range(n_bootstrap):
        sampled_counts, sampled_pred = resample(df[col], preds)
        if is_single_class(sampled_pred, col=pred_col):
            continue
        scores[i] = roc_auc_score(sampled_pred, sampled_counts)
    return scores

def cohort_bootstrap_auc(cohort, on, pred_col="is_benefit", n_bootstrap=1000, **kwargs):
    col, df = cohort.as_dataframe(on, return_cols=True, **kwargs)
    return bootstrap_auc(df=df,
                         col=col,
                         pred_col=pred_col,
                         n_bootstrap=n_bootstrap)

def cohort_mean_bootstrap_auc(cohort, on, pred_col="is_benefit", n_bootstrap=1000, **kwargs):
    return cohort_bootstrap_auc(cohort, on, pred_col, n_bootstrap, **kwargs).mean()

def coxph_model(formula, data, time_col, event_col, **kwargs):
    # pylint: disable=no-member
    # pylint gets confused by dmatrix
    sdata = patsy.dmatrix(
        formula,
        data = data,
        return_type = "dataframe").join(data[[time_col, event_col]])
    sdata = sdata.ix[:, sdata.columns != "Intercept"]
    if not(hasattr(kwargs, "penalizer")):
        kwargs["penalizer"] = 0.1
    if not(hasattr(kwargs, "normalize")):
        kwargs['normalize'] = False
    cf = ll.CoxPHFitter(**kwargs)
    cf.fit(sdata, time_col, event_col)
    cf.print_summary()
    return cf

def cohort_coxph(cohort, on, formula=None, how="pfs"):
    # If not specified, the formula is just the names of the functions
    cols, df = cohort.as_dataframe(on=on, return_cols=True)
    if formula is None:
        # Single col
        if type(cols) == str:
            formula = cols
        elif type(cols) == list:
            formula = " + ".join(cols)
    if how == "pfs":
        event_col = "is_progressed_or_deceased"
        time_col = "pfs"
    elif how == "os":
        event_col = "is_deceased"
        time_col = "os"
    return coxph_model(formula=formula,
                       data=df,
                       time_col=time_col,
                       event_col=event_col)
