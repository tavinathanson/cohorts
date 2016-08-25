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
from collections import namedtuple

from scipy.stats import mannwhitneyu, fisher_exact
import seaborn as sb
import pandas as pd
import matplotlib.pyplot as plt
from sklearn.metrics import roc_curve
from .model import bootstrap_auc

def vertical_percent(plot, percent=0.1):
    """
    Using the size of the y axis, return a fraction of that size.
    """
    plot_bottom, plot_top = plot.get_ylim()
    return percent * (plot_top - plot_bottom)

def as_numeric(text):
    try:
        return float(text)
    except:
        return None

def hide_ticks(plot, min_tick_value=None, max_tick_value=None):
    """Hide tick values that are outside of [min_tick_value, max_tick_value]"""
    for tick, tick_value in zip(plot.get_yticklabels(), plot.get_yticks()):
        tick_label = as_numeric(tick_value)
        if tick_label:
            if (min_tick_value is not None and tick_label < min_tick_value or 
                 max_tick_value is not None and tick_label > max_tick_value):
                tick.set_visible(False)

def hide_negative_y_ticks(plot):
    hide_ticks(plot, min_tick_value=0)

def only_percentage_ticks(plot):
    """
    Only show ticks from 0.0 to 1.0.
    """
    hide_ticks(plot, min_tick_value=0, max_tick_value=1.0)

def add_significance_indicator(plot, col_a=0, col_b=1, significant=False):
    """
    Add a p-value significance indicator.
    """
    plot_bottom, plot_top = plot.get_ylim()
    # Give the plot a little room for the significance indicator
    line_height = vertical_percent(plot, 0.1)
    # Add some extra spacing below the indicator
    plot_top = plot_top + line_height
    # Add some extra spacing above the indicator
    plot.set_ylim(top=plot_top + line_height * 2)
    color = "black"
    line_top = plot_top + line_height
    plot.plot([col_a, col_a, col_b, col_b], [plot_top, line_top, line_top, plot_top], lw=1.5, color=color)
    indicator = "*" if significant else "ns"
    plot.text((col_a + col_b) * 0.5, line_top, indicator, ha="center", va="bottom", color=color)

def stripboxplot(x, y, data, ax=None, significant=None, **kwargs):
    """
    Overlay a stripplot on top of a boxplot.
    """
    ax = sb.boxplot(
        x=x,
        y=y,
        data=data,
        ax=ax,
        fliersize=0,
        **kwargs
    )

    plot = sb.stripplot(
        x=x,
        y=y,
        data=data,
        ax=ax,
        jitter=kwargs.pop("jitter", 0.05),
        color=kwargs.pop("color", "0.3"),
        **kwargs
    )

    if data[y].min() >= 0:
        hide_negative_y_ticks(plot)
    if significant is not None:
        add_significance_indicator(plot=plot, significant=significant)

    return plot

def sided_str_from_alternative(alternative, condition):
    if alternative is None:
        raise ValueError("Must pick an alternative")
    if alternative == "two-sided":
        return alternative
    # alternative hypothesis: condition is 'less' or 'greater' than no-condition
    op_str = ">" if alternative == "greater" else "<"
    return "one-sided: %s %s not %s" % (condition, op_str, condition)

def get_condition_mask(df, condition, condition_value):
    if condition_value:
        condition_mask = df[condition] == condition_value
    else:
        # This is necessary in the event that condition has a non-bool dtype,
        # such as object. This may happen if a function returns np.nan in
        # addition to True/False (later filtered down to just True/False).
        # ~condition_mask will behave incorrectly if dtype is not bool.
        condition_mask = df[condition].astype("bool")
    return condition_mask

class FishersExactResults(namedtuple("FishersExactResults", ["oddsratio", "p_value", "sided_str", "with_condition1_series", "without_condition1_series", "plot"])):
    def __str__(self):
        return "FishersExactResults(oddsratio=%s, p_value=%s, sided_str='%s')" % (
            self.oddsratio, self.p_value, self.sided_str)

    def __repr__(self):
        return self.__str__()

def fishers_exact_plot(data, condition1, condition2, ax=None,
                       condition1_value=None,
                       alternative="two-sided", **kwargs):
    """
    Perform a Fisher's exact test to compare to binary columns

    Parameters
    ----------
    data: Pandas dataframe
        Dataframe to retrieve information from

    condition1: str
        First binary column to compare (and used for test sidedness)

    condition2: str
        Second binary column to compare

    ax : Axes, default None
        Axes to plot on

    condition1_value:
        If `condition1` is not a binary column, split on =/!= to condition1_value

    alternative:
        Specify the sidedness of the test: "two-sided", "less"
        or "greater"
    """
    plot = sb.barplot(
        x=condition1,
        y=condition2,
        ax=ax,
        data=data,
        **kwargs
    )

    plot.set_ylabel("Percent %s" % condition2)
    condition1_mask = get_condition_mask(data, condition1, condition1_value)
    count_table = pd.crosstab(data[condition1], data[condition2])
    print(count_table)
    oddsratio, p_value = fisher_exact(count_table, alternative=alternative)
    add_significance_indicator(plot=plot, significant=p_value <= 0.05)
    only_percentage_ticks(plot)

    if alternative != "two-sided":
        raise ValueError("We need to better understand the one-sided Fisher's Exact test")
    sided_str = "two-sided"
    print("Fisher's Exact Test: OR: {}, p-value={} ({})".format(oddsratio, p_value, sided_str))
    return FishersExactResults(oddsratio=oddsratio,
                               p_value=p_value,
                               sided_str=sided_str,
                               with_condition1_series=data[condition1_mask][condition2],
                               without_condition1_series=data[~condition1_mask][condition2],
                               plot=plot)

class MannWhitneyResults(namedtuple("MannWhitneyResults", ["U", "p_value", "sided_str", "with_condition_series", "without_condition_series", "plot"])):
    def __str__(self):
        return "MannWhitneyResults(U=%s, p_value=%s, sided_str='%s')" % (
            self.U, self.p_value, self.sided_str)

    def __repr__(self):
        return self.__str__()

def mann_whitney_plot(data,
                      condition,
                      distribution,
                      ax=None,
                      condition_value=None,
                      alternative="two-sided",
                      skip_plot=False,
                      **kwargs):
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
    condition_mask = get_condition_mask(data, condition, condition_value)
    U, p_value = mannwhitneyu(
        data[condition_mask][distribution],
        data[~condition_mask][distribution],
        alternative=alternative
    )

    plot = None
    if not skip_plot:
        plot = stripboxplot(
            x=condition,
            y=distribution,
            data=data,
            ax=ax,
            significant=p_value <= 0.05,
            **kwargs
        )

    sided_str = sided_str_from_alternative(alternative, condition)
    print("Mann-Whitney test: U={}, p-value={} ({})".format(U, p_value, sided_str))
    return MannWhitneyResults(U=U,
                              p_value=p_value,
                              sided_str=sided_str,
                              with_condition_series=data[condition_mask][distribution],
                              without_condition_series=data[~condition_mask][distribution],
                              plot=plot)

class CorrelationResults(namedtuple("CorrelationResults", ["coeff", "p_value", "stat_func", "series_x", "series_y", "plot"])):
    def __str__(self):
        return "CorrelationResults(coeff=%s, p_value=%s, stat_func=%s)" % (
            self.coeff, self.p_value, self.stat_func.__name__)

    def __repr__(self):
        return self.__str__()

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
