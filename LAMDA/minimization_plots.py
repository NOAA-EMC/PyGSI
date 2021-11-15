import numpy as np
import pandas as pd
from datetime import datetime
import matplotlib.pyplot as plt
from emcpy.plots.plots import LinePlot
from emcpy.plots import CreatePlot

__all__ = ['plot_minimization']


def _plot_cost_function(df, config, outdir):
    """
    Use data from dataframe to plot the cost function.
    """
    # Get cycle info
    cycle = df.index.get_level_values(0)[-1]
    cyclestr = datetime.strftime(cycle, '%Y%m%d%H')
    n_cycles = len(np.unique(df.index.get_level_values(0)))

    # Create approriate dataframes
    current_cycle_df = df.loc[cycle]
    avg_df = df.groupby(level=[1, 2]).mean()

    # Grab data
    j = current_cycle_df['J'].to_numpy()
    x = np.arange(len(j))

    avg_j = avg_df['J'].to_numpy()
    avg_x = np.arange(len(avg_j))

    # Create LinePlot objects
    cost_plot = LinePlot(x, j)
    cost_plot.linewidth = 2
    cost_plot.label = 'Cost'

    avg_cost_plot = LinePlot(avg_x, avg_j)
    avg_cost_plot.color = 'tab:red'
    avg_cost_plot.linewidth = 2
    avg_cost_plot.label = f'Last {n_cycles} cycles Average'

    # Create Plot
    myplot = CreatePlot()
    myplot.draw_data([cost_plot, avg_cost_plot])
    myplot.set_yscale('log')
    myplot.add_grid()
    myplot.set_xlim(0, len(j)-1)
    myplot.add_xlabel('Iterations')
    myplot.add_ylabel('log (J)')
    myplot.add_legend(loc='upper right',
                      fontsize='large')

    title = f"{config['experiment']} Cost - {config['tm']}"
    myplot.add_title(title, loc='left')
    myplot.add_title(cyclestr, loc='right',
                     fontweight='semibold')

    fig = myplot.return_figure()

    savefile = (f"{cyclestr}_{config['experiment']}_" +
                f"tm0{config['tm']}_cost_function.png")
    plt.savefig(outdir + savefile, bbox_inches='tight',
                pad_inches=0.1)
    plt.close('all')


def _plot_gnorm(df, config, outdir):
    """
    Use date from dataframe to plot gnorm.
    """
    # Get cycle info
    cycle = df.index.get_level_values(0)[-1]
    cyclestr = datetime.strftime(cycle, '%Y%m%d%H')
    n_cycles = len(np.unique(df.index.get_level_values(0)))

    # Create approriate dataframes
    current_cycle_df = df.loc[cycle]
    avg_df = df.groupby(level=[1, 2]).mean()

    # Grab data
    gJ = current_cycle_df['gJ'].to_numpy()
    gJ = np.log(gJ/gJ[0])
    x = np.arange(len(gJ))

    avg_gJ = avg_df['gJ'].to_numpy()
    avg_gJ = np.log(avg_gJ/avg_gJ[0])
    avg_x = np.arange(len(avg_gJ))

    # Create LinePlot objects
    gnorm = LinePlot(x, gJ)
    gnorm.linewidth = 2
    gnorm.label = 'gnorm'

    avg_gnorm = LinePlot(avg_x, avg_gJ)
    avg_gnorm.color = 'tab:red'
    avg_gnorm.linewidth = 2
    avg_gnorm.linestyle = '--'
    avg_gnorm.label = f'Last {n_cycles} cycles Average'

    # Create Plot
    myplot = CreatePlot()
    myplot.draw_data([gnorm, avg_gnorm])
    # myplot.set_yscale('log')
    myplot.add_grid()
    myplot.set_xlim(0, len(gJ)-1)
    myplot.add_xlabel('Iterations')
    myplot.add_ylabel('log (gnorm)')
    myplot.add_legend(loc='upper right',
                      fontsize='large')

    title = f"{config['experiment']} gnorm - {config['tm']}"
    myplot.add_title(title, loc='left')
    myplot.add_title(cyclestr, loc='right',
                     fontweight='semibold')

    fig = myplot.return_figure()

    savefile = (f"{cyclestr}_{config['experiment']}_" +
                f"tm0{config['tm']}_gnorm.png")
    plt.savefig(outdir + savefile, bbox_inches='tight',
                pad_inches=0.1)
    plt.close('all')


def plot_minimization(df, plotting_config, outdir):
    """
    Plot minimization plots including gnorm and cost function
    and save them to outdir.

    Args:
        df : (pandas dataframe) dataframe with appropriate information
             from GSI Stat file
        plotting_config : (dict) dictionary with information about
                          minimization period
        outdir : (str) path to output diagnostics
    """

    _plot_gnorm(df, plotting_config, outdir)
    _plot_cost_function(df, plotting_config, outdir)
