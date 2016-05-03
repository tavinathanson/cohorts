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

def stripboxplot(x, y, data, **kwargs):
    """
    Overlay a stripplot on top of a boxplot.
    """
    ax = sb.boxplot(
        x=x,
        y=y,
        data=data
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

def fishers_exact_plot(data, condition1, condition2):
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
    """
    plot = sb.factorplot(
        x=condition1,
        y=condition2,
        kind='bar',
        data=data
    )
    count_table = pd.crosstab(data[condition1], data[condition2])
    print(count_table)
    oddsratio, pvalue = fisher_exact(count_table)
    print("Fisher's Exact Test: OR: {}, p-value={}".format(oddsratio, pvalue))
    return (oddsratio, pvalue, plot)

def mann_whitney_plot(data, condition, distribution,
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
            data=data
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

    print("Mann-Whitney test: U={}, p-value={}".format(U, pvalue))
    return (U, pvalue, plot)
