import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pyGSI.gsi_stat import GSIstat
from emcpy.plots.plots import LinePlot, HorizontalLine
from emcpy.plots import CreatePlot


def _plot_bias_rmse_timeseries(df, config, outdir):
    """
    Used indexed df to plot rmse and bias.
    """
    df = df.reset_index()

    # Grab data from dataframe
    cycles = df['date'].dt.strftime("%d %b\n%Y %Hz")
    omf_bc = df['OmF_bc']
    omf_wobc = df['OmF_wobc']
    rmse = df['rms']

    # Create bias correction object
    bc_line = LinePlot(cycles, omf_bc)
    bc_line.label = 'Bias w/ BC'
    bc_line.marker = 'o'
    bc_line.markersize = 3

    # Create object without bias correction
    wobc_line = LinePlot(cycles, omf_wobc)
    wobc_line.color = 'tab:green'
    wobc_line.label = 'Bias w/o BC'
    wobc_line.marker = 'o'
    wobc_line.markersize = 3

    # Create rmse
    rmse_line = LinePlot(cycles, rmse)
    rmse_line.color = 'tab:red'
    rmse_line.label = 'RMSE'
    rmse_line.marker = 'o'
    rmse_line.markersize = 3

    # Create a line at 0
    zero_line = HorizontalLine(y=0)

    # Create plot and draw data
    myplt = CreatePlot(figsize=(10, 6))
    plt_list = [bc_line, wobc_line, rmse_line, zero_line]
    myplt.draw_data(plt_list)

    # Add features
    myplt.set_ylim(-5, 5)
    myplt.add_grid(linewidth=0.5, color='grey', linestyle='--')
    myplt.add_legend(loc='lower right', fontsize='large')

    title = (f"{config['bias type']} RMSE and Bias Time Series\n{config['sensor']} "
             f"{config['satellite']} Channel {config['channel']} tm0{config['tm']}")
    myplt.add_title(title, fontsize=14)

    # Return matplotlib figure
    fig = myplt.return_figure()

    # Save figure
    save_cycles = df['date'].dt.strftime('%Y%m%d%H').to_numpy()
    savefile = (f"{save_cycles[0]}_{save_cycles[-1]}_{config['sensor']}"
                f"_{config['satellite']}_channel_{config['channel']}"
                f"_{config['bias type']}_tm0{config['tm']}_rmse_bias_timeseries.png")
    fig.savefig('./' + savefile, bbox_inches='tight',
                pad_inches=0.1)
    plt.close('all')


def bias_rmse_timeseries(df, config, outdir):
    """
    Plots a timeseries of bias and rmse.

    Args:
        df : (pandas dataframe) multidimensional pandas dataframe
             with several cycles of gsi stats data
        config : (dict) dictionary including informaton about the data
                 being plotted
        outdir : (str) path to output figures
    """
    # Select data by satellite and channel
    for idx_col in ['satellite', 'channel']:
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
    config['bias type'] = 'OmF'
    _plot_bias_rmse_timeseries(omf_df, config, outdir)

    # Plot oma
    config['bias type'] = 'OmA'
    _plot_bias_rmse_timeseries(oma_df, config, outdir)
