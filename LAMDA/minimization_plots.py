import numpy as np
from datetime import datetime
import matplotlib.pyplot as plt
from pyGSI.gsi_stat import GSIstat
from emcpy.plots.plots import LinePlot
from emcpy.plots import CreatePlot

__all__ = ['plot_minimization']


def _plot_cost_single_cycle(cost_df, cyclestr, plotdir):
    """
    Uses cost dataframe to plot a single cycle of cost.
    """
    j = cost_df['J'].to_numpy()
    iterations = np.arange(len(j))

    cost_plot = LinePlot(iterations, j)
    cost_plot.linewidth = 2
    cost_plot.label = 'cost'

    myplot = CreatePlot()
    myplot.draw_data([cost_plot])
    myplot.set_yscale('log')
    myplot.add_grid()
    myplot.set_xlim(0, len(j)-1)
    myplot.add_xlabel('Iterations')
    myplot.add_ylabel('log (J)')
    myplot.add_legend(loc='upper right',
                      fontsize='large')
    myplot.add_title('FV3LAM Cost', loc='left')
    myplot.add_title(cyclestr, loc='right',
                     fontweight='semibold')

    fig = myplot.return_fig()

    plt.savefig(plotdir + f"single_cycle_cost_plot_{cyclestr}.png",
                bbox_inches='tight', pad_inches=0.1)
    plt.close('all')


def _plot_gnorm_single_cycle(cost_df, cyclestr, plotdir):
    """
    Uses cost dataframe to plot a single cycle of gnorm.
    """
    gj = cost_df['gJ'].to_numpy()
    gj = np.log(gj/gj[0])
    x = np.arange(len(j))

    cost_plot = LinePlot(x, gj)
    cost_plot.linewidth = 2
    cost_plot.label = 'gnorm'

    myplot = CreatePlot()
    myplot.draw_data([cost_plot])
    myplot.add_grid()
    myplot.set_xlim(0, len(j)-1)
    myplot.add_xlabel('Iterations')
    myplot.add_ylabel('log (gnorm)')
    myplot.add_legend(loc='upper right',
                      fontsize='large')
    myplot.add_title('FV3LAM gnorm', loc='left')
    myplot.add_title(cyclestr, loc='right',
                     fontweight='semibold')

    fig = myplot.return_fig()

    plt.savefig(plotdir + f"single_cycle_gnorm_plot_{cyclestr}.png",
                bbox_inches='tight', pad_inches=0.1)
    plt.close('all')


def plot_minimization(inputfile, cyclestr, plotdir):
    """
    Create minimization plots from gsistat file.

    Args:
        inputfile : (str) path to GSI stat file
        cyclestr : (str) cycle from GSI stat file
        plotdir : (str) Path to output plot directory
    """
    cycle = datatime.strptime(cyclestr, '%Y%m%d%H')

    # Get gdas object and extract cost info
    gdas = GSIstat(inputfile, cyclestr)
    cost_df = gdas.extract('cost')

    # Plot single cycle cost
    _plot_cost_single_cycle(cost_df, cyclestr, plotdir)

    # Plot single cycle gnorm
    _plot_gnorm_single_cycle(cost_df, cyclestr, plotdir)
