#!/usr/bin/env python
import argparse
import datetime
import matplotlib
from matplotlib import pyplot as plt
from matplotlib import rcParams, ticker
from matplotlib import gridspec as gspec
import os
import numpy as np
from pyGSI.gsi_stat import GSIstat

it = 'it == 1'
plottype = 'mean'
obtypes = [120, 220]
obvars = ['t', 'uv', 'q']
levels = [1000, 900, 800, 600, 400, 300, 250, 200, 150, 100, 50, 0]


def gen_figure(datadict, datatypestr, stattype, labels, sdate, edate, save, plotdir):
    # Line/marker colors for experiments ('k' is the first)
    mc = ['k', 'r', 'g', 'b', 'm', 'c', 'y']

    # set figure params one time only.
    rcParams['figure.subplot.left'] = 0.1
    rcParams['figure.subplot.top'] = 0.85
    rcParams['legend.fontsize'] = 12
    rcParams['axes.grid'] = True

    fig1 = plt.figure(figsize=(10, 8))
    plt.subplots_adjust(hspace=0.3)
    gs = gspec.GridSpec(1, 3)

    for v, var in enumerate(obvars):
        xmin = 999
        xmax = 0
        ax = plt.subplot(gs[v])
        for e, expid in enumerate(labels):
            profile = datadict[expid][stattype][var][:-1]
            ax.plot(profile, levels[:-1], marker='o', color=mc[e],
                    mfc=mc[e], mec=mc[e], label=labels[e])
            if (var in ['q']):
                xmin_, xmax_ = np.min(profile[:-1]), np.max(profile[:-1])
            else:
                xmin_, xmax_ = np.min(profile), np.max(profile)
            if (xmin_ < xmin):
                xmin = xmin_
            if (xmax_ > xmax):
                xmax = xmax_
        if (v in [0]):
            plt.legend(loc=0, numpoints=1)
        if (v in [0]):
            plt.ylabel('pressure (hPa)')

        if (var == 'uv'):
            var_unit = 'm/s'
            var_name = 'Winds'
        elif (var == 't'):
            var_unit = 'K'
            var_name = 'Temperature'
        elif (var == 'q'):
            var_unit = '%'
            var_name = 'Relative Humidity'

        if (stattype == 'sum'):
            plt.xlabel('count')
        else:
            plt.xlabel('magnitude (%s)' % var_unit)

        plt.title(var_name, fontsize=14)
        plt.ylim(1020, 50)
        ax.set_yscale('log')
        ax.yaxis.set_major_locator(ticker.LogLocator(base=10.0,
                                   subs=np.arange(1, 10)))
        ax.yaxis.set_major_formatter(ticker.FormatStrFormatter("%g"))
        xmin = xmin - (xmax-xmin)*0.1
        xmax = xmax + (xmax-xmin)*0.1
        plt.xlim(xmin, xmax)

    sdatestr = sdate.strftime('%Y%m%d%H')
    edatestr = edate.strftime('%Y%m%d%H')
    plt.figtext(0.5, 0.93, '%s O-F (%s-%s)'
                % (datatypestr, sdatestr, edatestr),
                horizontalalignment='center', fontsize=18)

    if (save_figure):
        fname = 'gsistat_uvtq'
        plt.savefig(plotdir+'/%s_%s.pdf' % (fname, datatypestr))
        plt.savefig(plotdir+'/%s_%s.png' % (fname, datatypestr))
    else:
        plt.show()


def get_gsistat_list(startdate, enddate):
    statfiles = []
    cycles = []
    mydate = startdate
    while mydate <= enddate:
        cycle = mydate.strftime('%Y%m%d%H')
        fname = 'gsistat.gdas.' + cycle
        statfiles.append(fname)
        cycles.append(cycle)
        mydate = mydate + datetime.timedelta(hours=6)
    return statfiles, cycles


if __name__ == '__main__':
    # get command line arguments
    parser = argparse.ArgumentParser(description=('Plots a comparison of GSI ',
                                                  'statistics for ',
                                                  'experiments compared to a ',
                                                  'reference control run.'))
    parser.add_argument('-d', '--gsistats',
                        help='list of directories containing GSI stat files',
                        nargs='+', required=True)
    parser.add_argument('-l', '--label',
                        help='list of labels for experiment IDs',
                        nargs='+', required=False)
    parser.add_argument('-f', '--save_figure',
                        help='save figures as png and pdf',
                        action='store_true', required=False)
    parser.add_argument('-p', '--plotdir',
                        help='path to where to save figures',
                        default='./', required=False)
    parser.add_argument('-s', '--start_date', help='starting date',
                        type=str, metavar='YYYYMMDDHH', required=True)
    parser.add_argument('-e', '--end_date', help='ending date',
                        type=str, metavar='YYYYMMDDHH', required=True)
    args = parser.parse_args()

    save_figure = args.save_figure
    if (save_figure):
        matplotlib.use('Agg')

    sdate = datetime.datetime.strptime(args.start_date, '%Y%m%d%H')
    edate = datetime.datetime.strptime(args.end_date, '%Y%m%d%H')

    statfiles, cycles = get_gsistat_list(sdate, edate)

    if args.label:
        labels = args.label
    else:
        labels = [g.rstrip('/').split('/')[-1] for g in args.gsistats]

    # loop through all files and variables and grab statistics
    rmses = {}
    counts = {}
    biases = {}
    for exp, gsistats in zip(labels, args.gsistats):
        rmses[exp] = {}
        counts[exp] = {}
        biases[exp] = {}
        for gsistat, cycle in zip(statfiles, cycles):
            rmses[exp][cycle] = {}
            counts[exp][cycle] = {}
            biases[exp][cycle] = {}
            inputfile = os.path.join(gsistats, gsistat)
            try:
                gdas = GSIstat(inputfile, cycle)
            except FileNotFoundError:
                raise FileNotFoundError(
                      f'Unable to find {inputfile} for cycle {cycle}')
            # now loop through variables
            for var in obvars:
                stat = gdas.extract(var)  # t, uv, q, etc.
                stat = stat.query(it)  # ges (1) or anl (3) ?
                tmpstat = stat[stat.index.isin(obtypes, level='typ')]
                tmpstat = tmpstat[tmpstat.index.isin(['asm'], level='use')]
                rmses[exp][cycle][var] = tmpstat[tmpstat.index.isin(
                                         ['rms'],
                                         level='stat')]
                counts[exp][cycle][var] = tmpstat[tmpstat.index.isin(
                                          ['count'],
                                          level='stat')]
                biases[exp][cycle][var] = tmpstat[tmpstat.index.isin(
                                          ['bias'],
                                          level='stat')]

    # now aggregate stats
    for exp in labels:
        rmses[exp]['mean'] = {}
        rmses[exp]['aggr'] = {}
        biases[exp]['mean'] = {}
        biases[exp]['aggr'] = {}
        counts[exp]['sum'] = {}
        for var in obvars:
            rmse_var = np.empty([len(cycles), len(levels)])
            bias_var = np.empty([len(cycles), len(levels)])
            counts_var = np.empty([len(cycles), len(levels)], dtype=int)
            for i, cycle in enumerate(cycles):
                rmse_var[i, :] = rmses[exp][cycle][var].values[0]
                bias_var[i, :] = biases[exp][cycle][var].values[0]
                counts_var[i, :] = counts[exp][cycle][var].values[0]
            # Compute mean rms, bias
            rmses[exp]['mean'][var] = rmse_var.mean(axis=0)
            biases[exp]['mean'][var] = bias_var.mean(axis=0)
            # Compute aggregate rms, bias
            ar = np.asarray([])
            ab = np.asarray([])
            for j in range(np.ma.size(rmse_var, axis=1)):
                r = rmse_var[:, j].squeeze()
                b = bias_var[:, j].squeeze()
                c = counts_var[:, j].squeeze()
                if (np.sum(c) > 0):
                    ar = np.append(ar, np.sqrt(np.sum(np.multiply(c, r**2.))/np.sum(c)))
                    ab = np.append(ab, np.sum(np.multiply(c, b))/np.sum(c))
                else:
                    ar = np.append(ar, np.nan)
                    ab = np.append(ab, np.nan)
            rmses[exp]['aggr'][var] = ar
            biases[exp]['aggr'][var] = ab
            # Compute summed counts
            counts[exp]['sum'][var] = counts_var.sum(axis=0)

    # make figures
    gen_figure(rmses, 'RMSE', plottype, labels, sdate, edate, save_figure, args.plotdir)
    gen_figure(biases, 'Bias', plottype, labels, sdate, edate, save_figure, args.plotdir)
    gen_figure(counts, 'Count', 'sum', labels, sdate, edate,
               save_figure, args.plotdir)
