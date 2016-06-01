import seaborn as sb
import matplotlib.pyplot as plt
import math

def _bar_plot(data,
             title=None,
             figsize=None,
             colormap=None,
             ax=None):
    plot = data.plot(title=title,
           kind='bar', 
           stacked=True, 
           figsize=figsize,
           ax=ax,
           colormap=colormap).axes.get_xaxis().set_visible(False)
    return plot

def _binned_bar_plot(data,
             sample_col,
             bin_by,
             title=None,
             figsize=None,
             colormap=None,
             ax=None):

    by_bin = data.groupby(
        [sample_col, bin_by]
    )[[sample_col]].count().unstack()

    by_bin.columns = by_bin.columns.get_level_values(1)
    return _bar_plot(by_bin, 
                    title, 
                    figsize,
                    colormap,
                    ax)

def _indicator_plot(data,
                   sample_col,
                   indicator_col,
                   colormap=None,
                   figsize=None,
                   ax=None):
    indicator_data = data.set_index([sample_col])[[indicator_col]].T
    indicator_plot = sb.heatmap(indicator_data,
               square=True,
               cbar=None,
               xticklabels=True,
               linewidths=1,
               cmap=colormap,
               ax=ax)
    plt.setp(indicator_plot.axes.get_yticklabels(), rotation=0)
    return indicator_plot

def landscape_plot(cohort,
                   effects_df,
                   sample_col,
                   width=10,
                   bar_height=4,
                   bin_columns=[],
                   indicator_columns=[],
                   value_columns=[],):
    """Build a figure with a plot from each of the columns specified

    Parameters
    ----------
    cohort : Cohort
    effects_df : str
        Description
    sample_col : str
        Column that represents a sample id
    width : int, optional
        Width of the full figure
    bar_height : int, optional
        Height of each bar plot
    bin_columns : list, optional
         Columns to build a binned/stacked bar plot from
    indicator_columns : list, optional
         Columns to build a indicator bar from 
    value_columns : list, optional
        Columns to build a bar plot from
    """
    cohort_size = len(cohort)
    min_square_size = float(width) / (.75 * cohort_size)
    num_bar_plots = len(bin_columns) + len(value_columns)
    height = len(indicator_columns) * min_square_size + num_bar_plots * bar_height 

    grid_rows = int(math.ceil(float(height) / min_square_size)) + num_bar_plots
    indicator_column_rows = len(indicator_columns)

    bar_rows = int((grid_rows - indicator_column_rows - num_bar_plots) / num_bar_plots)
    gridsize = (grid_rows, 1)

    plt.figure(0, figsize=(width - 1, height))

    current_row = 0
    for bin_by_col in bin_columns:
        ax = plt.subplot2grid(gridsize,
                       (current_row, 0),
                       colspan=1,
                       rowspan=bar_rows)
        _binned_bar_plot(effects_df,
                        sample_col,
                        bin_by_col,
                        ax=ax,)
        current_row += bar_rows + 1

    for on in value_columns:
        ax = plt.subplot2grid(gridsize,
                              (current_row, 0),
                              colspan=1,
                              rowspan=bar_rows)
        plot_col, df = cohort.as_dataframe(on)
        _bar_plot(df[plot_col],
                  ax=ax,
                  title=plot_col,)
        current_row += bar_rows + 1

    current_row -= 1
    for (idx, indicator_on) in enumerate(indicator_columns):
        ax = plt.subplot2grid(gridsize,
                              (current_row, 0))
        indicator_col, df = cohort.as_dataframe(indicator_on)
        ip = _indicator_plot(df,
                            sample_col,
                            indicator_col,
                            ax=ax,)
        current_row += 1
        if idx != len(indicator_columns) - 1:
            ip.axes.xaxis.set_visible(False)
