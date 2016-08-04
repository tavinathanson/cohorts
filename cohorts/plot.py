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

from scipy.stats import mannwhitneyu, fisher_exact
import seaborn as sb
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve
from .model import bootstrap_auc

def stripboxplot(x, y, data, ax=None, **kwargs):
    """
    Overlay a stripplot on top of a boxplot.
    """
    ax = sb.boxplot(
        x=x,
        y=y,
        data=data,
        ax=ax, 
        fliersize=0
    )

    return sb.stripplot(
        x=x,
        y=y,
        data=data,
        ax=ax,
        jitter=kwargs.pop("jitter", 0.05),
        color=kwargs.pop("color", "0.3"),
        **kwargs
    )

def fishers_exact_plot(data, condition1, condition2, ax=None):
    """
    Perform a Fisher's exact test to compare to binary columns

    Parameters
    ----------
    data: Pandas dataframe
        Dataframe to retrieve information from

    condition1: str
        First binary column compare

    condition2: str
        Second binary column to compare

    ax : Axes, default None
        Axes to plot on
    """
    plot = sb.barplot(
        x=condition1,
        y=condition2,
        ax=ax,
        data=data
    )
    count_table = pd.crosstab(data[condition1], data[condition2])
    print(count_table)
    oddsratio, pvalue = fisher_exact(count_table)
    print("Fisher's Exact Test: OR: {}, p-value={}".format(oddsratio, pvalue))
    return (oddsratio, pvalue, plot)

def mann_whitney_plot(data, condition, distribution, ax=None,
                      condition_value=None, alternative="two-sided",
                      skip_plot=False):
    """
    Create a box plot comparing a condition and perform a
    Mann Whitney test to compare the distribution in condition A v B

    Parameters
    ----------
    data: Pandas dataframe
        Dataframe to retrieve information from

    condition: str
        Column to use as the splitting criteria

    distribution: str
        Column to use as the Y-axis or distribution in the test

    ax : Axes, default None
        Axes to plot on

    condition_value:
        If `condition` is not a binary column, split on =/!= to condition_value

    alternative:
        Specify the sidedness of the Mann-Whitney test: "two-sided", "less"
        or "greater"

    skip_plot:
        Calculate the test statistic and p-value, but don't plot.
    """
    plot = None
    if not skip_plot:
        plot = stripboxplot(
            x=condition,
            y=distribution,
            data=data,
            ax=ax
        )

    if condition_value:
        condition_mask = data[condition] == condition_value
    else:
        condition_mask = data[condition]
    U, pvalue = mannwhitneyu(
        data[condition_mask][distribution],
        data[~condition_mask][distribution],
        alternative=alternative
    )

    if alternative is None:
        raise ValueError("Must pick a mannwhitneyu alternative")
    if alternative == "two-sided":
        sided_str = alternative
    else:
        # alternative hypothesis: condition is 'less' or 'greater' than no-condition
        op_str = ">" if alternative == "greater" else "<"
        sided_str = "one-sided: %s %s not %s" % (condition, op_str, condition)
    print("Mann-Whitney test: U={}, p-value={} ({})".format(U, pvalue, sided_str))
    return (U, pvalue, plot)

def roc_curve_plot(data, value_column, outcome_column, bootstrap_samples=100, ax=None):
    """Create a ROC curve and compute the bootstrap AUC for the given variable and outcome

    Parameters
    ----------
    data : Pandas dataframe
        Dataframe to retrieve information from
    value_column : str
        Column to retrieve the values from
    outcome_column : str
        Column to use as the outcome variable
    bootstrap_samples : int, optional
        Number of bootstrap samples to use to compute the AUC
    ax : Axes, default None
        Axes to plot on

    Returns
    -------
    (mean_bootstrap_auc, roc_plot) : (float, matplotlib plot)
        Mean AUC for the given number of bootstrap samples and the plot
    """
    scores = bootstrap_auc(df=data,
                           col=value_column,
                           pred_col=outcome_column,
                           n_bootstrap=bootstrap_samples)
    mean_bootstrap_auc = scores.mean()
    print("{}, Bootstrap (samples = {}) AUC:{}, std={}".format(
        value_column, bootstrap_samples, mean_bootstrap_auc, scores.std()))

    outcome = data[outcome_column].astype(int)
    values = data[value_column]
    fpr, tpr, thresholds = roc_curve(outcome, values)

    if ax is None:
        ax = plt.gca()

    roc_plot = ax.plot(fpr, tpr, lw=1, label=value_column)

    ax.set_xlim([-0.05, 1.05])
    ax.set_ylim([-0.05, 1.05])
    ax.set_xlabel('False Positive Rate')
    ax.set_ylabel('True Positive Rate')
    ax.legend(loc=2, borderaxespad=0.)
    ax.set_title('{} ROC Curve (n={})'.format(value_column, len(values)))

    return (mean_bootstrap_auc, roc_plot)
