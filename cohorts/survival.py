# -*- coding: utf-8 -*-

# Copyright (c) 2017. Mount Sinai School of Medicine
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

from lifelines import KaplanMeierFitter, CoxPHFitter
from lifelines.statistics import logrank_test
import matplotlib.colors as colors
from matplotlib import pyplot as plt
import seaborn as sb
import patsy
from .rounding import float_str

def plot_kmf(df,
             condition_col,
             censor_col,
             survival_col,
             threshold=None,
             title=None,
             xlabel=None,
             ylabel=None,
             ax=None,
             with_condition_color="#B38600",
             no_condition_color="#A941AC",
             with_condition_label=None,
             no_condition_label=None,
             color_map=None,
             label_map=None,
             color_palette="Set1",
             ci_show=False,
             print_as_title=False):
    """
    Plot survival curves by splitting the dataset into two groups based on
    condition_col

    if threshold is defined, the groups are split based on being > or <
    condition_col

    if threshold == 'median', the threshold is set to the median of condition_col

    Parameters
    ----------
        df: dataframe
        condition_col: string, column which contains the condition to split on
        survival_col: string, column which contains the survival time
        censor_col: string,
        threshold: int or string, if int, condition_col is thresholded,
                                  if 'median', condition_col thresholded
                                  at its median
        title: Title for the plot, default None
        ax: an existing matplotlib ax, optional, default None
        with_condition_color: str, hex code color for the with-condition curve
        no_condition_color: str, hex code color for the no-condition curve
        with_condition_label: str, optional, label for TRUE condition case
        no_condition_label: str, optional, label for FALSE condition case
        color_map: dict, optional, mapping of hex-values to condition text
          in the form of {value_name: color_hex_code}.
          defaults to `sb.color_palette` using `default_color_palette` name,
          or *_condition_color options in case of boolean operators.
        label_map: dict, optional, mapping of labels to condition text.
          defaults to "condition_name = condition_value", or *_condition_label
          options in case of boolean operators.
        color_palette: str, optional, name of sb.color_palette to use
          if color_map not provided.
        print_as_title: bool, optional, whether or not to print text
          within the plot's title vs. stdout, default False
    """
    if colors.is_color_like(with_condition_color):
        with_condition_color = colors.to_hex(with_condition_color)
    if colors.is_color_like(no_condition_color):
        no_condition_color = colors.to_hex(no_condition_color)
    kmf = KaplanMeierFitter()
    if threshold is not None:
        is_median = threshold == "median"
        if is_median:
            threshold = df[condition_col].median()
        label_suffix = float_str(threshold)
        condition = df[condition_col] > threshold
        default_label_no_condition = "%s ≤ %s" % (condition_col, label_suffix)
        if is_median:
            label_suffix += " (median)"
        default_label_with_condition = "%s > %s" % (condition_col, label_suffix)
        with_condition_label = with_condition_label or default_label_with_condition
        no_condition_label = no_condition_label or default_label_no_condition
        if not label_map:
            label_map = {False: no_condition_label,
                         True: with_condition_label}
        if not color_map:
            color_map = {False: no_condition_color,
                         True: with_condition_color}
    elif df[condition_col].dtype == 'O' or df[condition_col].dtype.name == "category":
        condition = df[condition_col].astype("category")
        if not label_map:
            label_map = dict()
            [label_map.update({condition_value: '{} = {}'.format(condition_col,
                                                        condition_value)})
                     for condition_value in condition.unique()]
        if not color_map:
            rgb_values = sb.color_palette(color_palette, len(label_map.keys()))
            hex_values = [colors.to_hex(col) for col in rgb_values]
            color_map = dict(zip(label_map.keys(), hex_values))
    elif df[condition_col].dtype == 'bool':
        condition = df[condition_col]
        default_label_with_condition = "= {}".format(condition_col)
        default_label_no_condition = "¬ {}".format(condition_col)
        with_condition_label = with_condition_label or default_label_with_condition
        no_condition_label = no_condition_label or default_label_no_condition
        if not label_map:
            label_map = {False: no_condition_label,
                         True: with_condition_label}
        if not color_map:
            color_map = {False: no_condition_color,
                         True: with_condition_color}
    else:
        raise ValueError('Don\'t know how to plot data of type\
                         {}'.format(df[condition_col].dtype))

    grp_desc = list()
    grp_survival_data = dict()
    grp_event_data = dict()
    grp_names = list(condition.unique())
    for grp_name, grp_df in df.groupby(condition):
        grp_survival = grp_df[survival_col]
        grp_event = (grp_df[censor_col].astype(bool))
        grp_label = label_map[grp_name]
        grp_color = color_map[grp_name]
        kmf.fit(grp_survival, grp_event, label=grp_label)
        desc_str = "# {}: {}".format(grp_label, len(grp_survival))
        grp_desc.append(desc_str)
        grp_survival_data[grp_name] = grp_survival
        grp_event_data[grp_name] = grp_event
        if ax:
            ax = kmf.plot(ax=ax, show_censors=True, ci_show=ci_show, color=grp_color)
        else:
            ax = kmf.plot(show_censors=True, ci_show=ci_show, color=grp_color)

    # Set the y-axis to range 0 to 1
    ax.set_ylim(0, 1)
    y_tick_vals = ax.get_yticks()
    ax.set_yticklabels(["%d" % int(y_tick_val * 100) for y_tick_val in y_tick_vals])

    if title:
        ax.set_title(title)
    elif print_as_title:
        ax.set_title(' | '.join(grp_desc))
    else:
        [print(desc) for desc in grp_desc]

    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)

    if len(grp_names) == 2:
        results = logrank_test(grp_survival_data[grp_names[0]],
                               grp_survival_data[grp_names[1]],
                               event_observed_A=grp_event_data[grp_names[0]],
                               event_observed_B=grp_event_data[grp_names[1]])
    else:
        cf = CoxPHFitter()
        cox_df = patsy.dmatrix('+'.join([condition_col, survival_col,
                                         censor_col]),
                               df, return_type='dataframe')
        del cox_df['Intercept']
        results = cf.fit(cox_df, survival_col, event_col=censor_col)
        results.print_summary()
    results.survival_data_series = grp_survival_data
    results.event_data_series = grp_event_data
    return results


def logrank(df,
            condition_col,
            censor_col,
            survival_col,
            threshold=None):
    if threshold is not None:
        if threshold == 'median':
            threshold = df[condition_col].median()
        condition = df[condition_col] > threshold
    else:
        condition = df[condition_col]
    df_with_condition = df[condition]
    df_no_condition = df[~condition]
    survival_no_condition = df_no_condition[survival_col]
    survival_with_condition = df_with_condition[survival_col]
    event_no_condition = (df_no_condition[censor_col].astype(bool))
    event_with_condition = (df_with_condition[censor_col].astype(bool))
    return logrank_test(survival_no_condition,
                        survival_with_condition,
                        event_observed_A=event_no_condition,
                        event_observed_B=event_with_condition)
