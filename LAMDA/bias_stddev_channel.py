import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pyGSI.gsi_stat import GSIstat
from emcpy.plots.plots import LinePlot, HorizontalLine
from emcpy.plots import CreatePlot


def _plot_bias_stddev_channel(df, config, outdir):
    """
    Plot the bias and standard deviation per channel
    for a single cycle and an average for multiple
    cycles.
    """
    # Create single cycle and average dataframe
    cycles = df.index.get_level_values(0).unique()
    n_days = len(cycles)/2
    scyc_df = df.loc[cycles[-1]]
    avg_df = df.groupby(level=[1, 2, 3, 4]).mean()

    # Get data
    channels = df.index.get_level_values(-1).unique().to_numpy()
    omf_bc = scyc_df['OmF_bc']
    omf_std = scyc_df['std']
    avg_bc = avg_df['OmF_bc']
    avg_std = avg_df['std']

    # Create bias correction object
    bc_line = LinePlot(channels, omf_bc)
    bc_line.color = 'tab:green'
    bc_line.label = 'OmF'
    bc_line.marker = 'o'
    bc_line.markersize = 4
    bc_line.linewidth = 3

    # Create standard deviation object
    std_dev = LinePlot(channels, omf_std)
    std_dev.color = 'tab:orange'
    std_dev.label = 'Std Dev'
    std_dev.marker = 'o'
    std_dev.markersize = 4
    std_dev.linewidth = 3

    # Create average bias correction object
    avg_bc_line = LinePlot(channels, avg_bc)
    avg_bc_line.color = 'tab:green'
    avg_bc_line.label = f'{n_days} Day OmF Average'
    avg_bc_line.marker = 'o'
    avg_bc_line.markersize = 4
    avg_bc_line.linewidth = 3
    avg_bc_line.linestyle = '--'

    # Create average standard deviation object
    avg_std_dev = LinePlot(channels, avg_std)
    avg_std_dev.color = 'tab:orange'
    avg_std_dev.label = f'{n_days} Day Std Average'
    avg_std_dev.marker = 'o'
    avg_std_dev.markersize = 4
    avg_std_dev.linewidth = 3
    avg_std_dev.linestyle = '--'

    # Create a line at 0
    zero_line = HorizontalLine(y=0)
    zero_line.linewidth = 1.25

    # Create plot and draw data
    myplt = CreatePlot(figsize=(10, 6))
    plt_list = [bc_line, avg_bc_line, std_dev, avg_std_dev, zero_line]
    myplt.draw_data(plt_list)

    # Add features
    myplt.set_ylim(-3, 3)
    myplt.set_xticks(channels)
    myplt.add_grid(linewidth=0.5, color='grey', linestyle='--')
    myplt.add_legend(loc='lower right', fontsize='large')
    myplt.add_xlabel('Channels')

    str_cycles = cycles.strftime('%Y%m%d%H').to_numpy()
    left_title = (f"{config['bias type']} (with Bias Correction) - Observed (K)"
                  f"\n{config['sensor']}_{config['satellite']}")
    right_title = f"{str_cycles[-1]}"
    myplt.add_title(left_title, loc='left', fontsize=14)
    myplt.add_title(right_title, loc='right', fontsize=12, fontweight='semibold')

    # Return matplotlib figure
    fig = myplt.return_figure()
    savefile = (f"{str_cycles[-1]}_{config['sensor']}_{config['satellite']}"
                f"_channel_{config['channel']}_{config['bias type']}"
                f"_bias_std_channel.png")
    fig.savefig('./' + savefile, bbox_inches='tight',
                pad_inches=0.1)
    plt.close('all')


def bias_stddev_channel(df, config, outdir):
    """
    Plots bias and standard deviation per channel.

    Args:
        df : (pandas dataframe) multidimensional pandas dataframe
             with several cycles of gsi stats data
        config : (dict) dictionary including informaton about the data
                 being plotted
        outdir : (str) path to output figures
    """
    # Select data by satellite and channel
    for idx_col in ['satellite']:
        indx = df.index.get_level_values(idx_col) == ''
        indx = np.ma.logical_or(indx, df.index.get_level_values(idx_col) == config[idx_col])
        df = df.iloc[indx]

    # Create omf and oma df
    indx = df.index.get_level_values(idx_col) == ''
    omf_indx = np.ma.logical_or(indx, df.index.get_level_values('it') == 1)
    omf_df = df.iloc[omf_indx]

    oma_indx = np.ma.logical_or(indx, df.index.get_level_values('it') == 3)
    oma_df = df.iloc[oma_indx]

    # Plot omf
    config['bias type'] = 'Ges'
    _plot_bias_stddev_channel(omf_df, config, outdir)

    # Plot oma
    config['bias type'] = 'Anl'
    _plot_bias_stddev_channel(oma_df, config, outdir)
