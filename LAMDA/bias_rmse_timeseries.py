import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pyGSI.gsi_stat import GSIstat
from emcpy.plots.plots import LinePlot, HorizontalLine
from emcpy.plots import CreatePlot


def _get_conventional_data(df, config):
    """
    Grabs appropriate data from the dataframe
    and creates EMCPy plot list, title and
    save file name for conventional data.
    """

    indx = df.index.get_level_values('stat') == ''
    bias_indx = np.ma.logical_or(indx, df.index.get_level_values('stat') == 'bias')
    rmse_indx = np.ma.logical_or(indx, df.index.get_level_values('stat') == 'rms')
    bias_df = df.iloc[bias_indx].reset_index()
    rmse_df = df.iloc[rmse_indx].reset_index()

    data_col = rmse_df.columns[-1]

    # Grab data from dataframe
    cycles = rmse_df['date'].dt.strftime("%d %b\n%Y %Hz")
    rmse = rmse_df[data_col]
    bias = bias_df[data_col]

    # Create bias object
    bias_line = LinePlot(cycles, bias)
    bias_line.label = 'Bias'
    bias_line.marker = 'o'
    bias_line.markersize = 3

    # Create rmse
    rmse_line = LinePlot(cycles, rmse)
    rmse_line.color = 'tab:green'
    rmse_line.label = 'RMSE'
    rmse_line.marker = 'o'
    rmse_line.markersize = 3

    # Create a line at 0
    zero_line = HorizontalLine(y=0)

    # Add objects to plot list
    plot_list = [bias_line, rmse_line, zero_line]

    # Create title
    title = (f"{config['experiment name']} {config['bias type']} "
             f"RMSE and Bias Time Series\nOb Type: {config['ob type']} "
             f"Type: {config['obsid']} Subtype: {config['subtype']} "
             f"tm0{config['tm']}")
    config['title'] = title

    # Create save file name
    save_cycles = bias_df['date'].dt.strftime('%Y%m%d%H').to_numpy()
    savefile = (f"{save_cycles[0]}_{save_cycles[-1]}_{config['experiment name']}_"
                f"{config['ob type']}_{config['obsid']}_{config['subtype']}_"
                f"{config['bias type']}_tm0{config['tm']}_rmse_bias_timeseries.png")
    config['save file'] = savefile

    return plot_list, config


def _get_radiance_data(df, config):
    """
    Grabs appropriate data from the dataframe
    and creates EMCPy plot list, title and
    save file name for radiance data.
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
    bc_line.linestyle = '-.'

    # Create object without bias correction
    wobc_line = LinePlot(cycles, omf_wobc)
    wobc_line.color = 'tab:green'
    wobc_line.label = 'Bias w/o BC'
    wobc_line.marker = 'o'
    wobc_line.markersize = 3
    wobc_line.linestyle = '--'

    # Create rmse
    rmse_line = LinePlot(cycles, rmse)
    rmse_line.color = 'tab:brown'
    rmse_line.label = 'RMSE'
    rmse_line.marker = 'o'
    rmse_line.markersize = 3

    # Create a line at 0
    zero_line = HorizontalLine(y=0)

    # Add objects to plot list
    plot_list = [bc_line, wobc_line, rmse_line, zero_line]

    # Create title
    title = (f"{config['experiment name']} {config['bias type']} "
             f"RMSE and Bias Time Series \n{config['sensor']} "
             f"{config['satellite']} Channel {config['channel']} "
             f"tm0{config['tm']}")
    config['title'] = title

    # Create save file name
    save_cycles = df['date'].dt.strftime('%Y%m%d%H').to_numpy()

    savefile = (f"{save_cycles[0]}_{save_cycles[-1]}_{config['experiment name']}_"
                f"{config['sensor']}_{config['satellite']}_channel_{config['channel']}"
                f"_{config['bias type']}_tm0{config['tm']}_rmse_bias_timeseries.png")
    config['save file'] = savefile

    return plot_list, config


def _plot_bias_rmse_timeseries(df, config, outdir):
    """
    Used indexed df to plot rmse and bias.
    """

    if config['data type'] == 'radiance':
        plot_list, config = _get_radiance_data(df, config)

    elif config['data type'] == 'conventional':
        plot_list, config = _get_conventional_data(df, config)

    # Create plot and draw data
    myplt = CreatePlot(figsize=(10, 6))
    myplt.draw_data(plot_list)

    # Add features
    myplt.set_ylim(-5, 5)
    myplt.add_grid(linewidth=0.5, color='grey', linestyle='--')
    myplt.add_legend(loc='lower right', fontsize='large')
    myplt.add_title(config['title'], fontsize=14)

    # Return matplotlib figure and save
    fig = myplt.return_figure()
    fig.savefig(outdir + config['save file'], bbox_inches='tight',
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

    if config['data type'] == 'radiance':

        # Select data by satellite and channel
        for idx_col in ['satellite', 'channel']:
            indx = df.index.get_level_values(idx_col) == ''
            indx = np.ma.logical_or(indx, df.index.get_level_values(idx_col) == config[idx_col])
            df = df.iloc[indx]

    elif config['data type'] == 'conventional':
        for idx_col in ['typ', 'use']:
            d = {
                'typ': config['obsid'],
                'styp': config['subtype'],
                'use': 'asm'
            }

            config_val = d[idx_col]
            indx = df.index.get_level_values(idx_col) == ''
            indx = np.ma.logical_or(indx, df.index.get_level_values(idx_col) == config_val)
            df = df.iloc[indx]

    # Create omf and oma df
    indx = df.index.get_level_values('it') == ''
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
