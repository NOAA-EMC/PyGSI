#!/usr/bin/env python3
# generate bar graphs of obs counts
# from GSI stat ASCII files
import numpy as np
import datetime as dt
import matplotlib.pyplot as plt
import argparse
from pyGSI.gsi_stat import GSIstat
import matplotlib
matplotlib.use('agg')

# rad_instruments = ['abi','ahi','airs','amsr2','amsre_hig','amsre_low',
#                   'amsre_mid', 'amsua','amsub','atms','avhrr','cris',
#                   'cris-fsr','gmi','goes_img','hirs3','hirs4','iasi',
#                   'mhs','saphir','seviri','sndrd1','sndrd2',
#                   'sndrd3','sndrd4','ssmi','ssmis']
rad_instruments = ['abi', 'ahi',
                   'amsua', 'atms', 'avhrr', 'cris-fsr',
                   'iasi', 'mhs', 'seviri', 'ssmis']
oz_instruments = ['omi', 'ompsnp', 'ompstc8', 'sbuv2']

conv_instruments = {
    "aircraft": [130, 131, 132, 133, 134, 135, 136, 137, 138, 139,
                 230, 231, 232, 233, 234, 235, 236, 237, 238, 239],
    "sondes": [120, 121, 122, 220, 221, 222],
    "rass": [126],
    "sfcship": [180, 183, 280, 282, 284],
    "sfc": [181, 187, 281, 287],
    "gps": [3, 4, 42, 43, 745, 825],
    "sst": [181, 182, 183, 189, 190, 191, 192, 193, 194, 195,
            196, 197, 198, 199, 200, 201, 202],
    "scatwind": [290],
    "windprof": range(227, 230),
    "vadwind": [224],
    "satwind": range(240, 261),
}


def plot_obscounts(inputfile, cyclestr, plotdir):
    cycle = dt.datetime.strptime(cyclestr, '%Y%m%d%H')
    gdas = GSIstat(inputfile, cyclestr)
    # empty lists to append
    obsname = []
    nassim = []
    nread = []
    nkeep = []
    # radiance obs
    radstat = gdas.extract('rad')
    it1 = radstat.query('it == 1')
    for rad in rad_instruments:
        radvar = it1.query('instrument == "'+rad+'"')
        obsname.append(rad)
        nassim.append(radvar['assim'].sum())
        nread.append(radvar['read'].sum())
        nkeep.append(radvar['keep'].sum())
    # ozone obs
    ozstat = gdas.extract('oz')
    it1 = ozstat.query('it == 1')
    for oz in oz_instruments:
        ozvar = it1.query('instrument == "'+oz+'"')
        obsname.append(oz)
        nassim.append(ozvar['assim'].sum())
        nread.append(ozvar['read'].sum())
        nkeep.append(ozvar['keep'].sum())
    # conventional obs
    tstat = gdas.extract('t')
    tstat = tstat.query('it == 1')
    uvstat = gdas.extract('uv')
    uvstat = uvstat.query('it == 1')
    qstat = gdas.extract('q')
    qstat = qstat.query('it == 1')
    gpsstat = gdas.extract('gps')
    gpsstat = gpsstat.query('it == 1')
    psstat = gdas.extract('ps')
    psstat = psstat.query('it == 1')
    for conv, codes in conv_instruments.items():
        tmpassim = []
        tmpread = []
        # temperature
        tmptstat = tstat[tstat.index.isin(codes, level='typ')]
        tmptstat = tmptstat[tmptstat.index.isin(['count'], level='stat')]
        tmpread.append(tmptstat['column'].sum())
        tmptstat = tmptstat[tmptstat.index.isin(['asm'], level='use')]
        tmpassim.append(tmptstat['column'].sum())
        # winds
        tmpuvstat = uvstat[uvstat.index.isin(codes, level='typ')]
        tmpuvstat = tmpuvstat[tmpuvstat.index.isin(['count'], level='stat')]
        tmpread.append(tmpuvstat['column'].sum())
        tmpuvstat = tmpuvstat[tmpuvstat.index.isin(['asm'], level='use')]
        tmpassim.append(tmpuvstat['column'].sum())
        # humidity
        tmpqstat = qstat[qstat.index.isin(codes, level='typ')]
        tmpqstat = tmpqstat[tmpqstat.index.isin(['count'], level='stat')]
        tmpread.append(tmpqstat['column'].sum())
        tmpqstat = tmpqstat[tmpqstat.index.isin(['asm'], level='use')]
        tmpassim.append(tmpqstat['column'].sum())
        # gps
        tmpgpsstat = gpsstat[gpsstat.index.isin(codes, level='typ')]
        tmpgpsstat = tmpgpsstat[tmpgpsstat.index.isin(['count'], level='stat')]
        tmpread.append(tmpgpsstat['column'].sum())
        tmpgpsstat = tmpgpsstat[tmpgpsstat.index.isin(['asm'], level='use')]
        tmpassim.append(tmpgpsstat['column'].sum())
        # ps
        tmppsstat = psstat[psstat.index.isin(codes, level='typ')]
        tmpread.append(tmppsstat['count'].sum())
        tmppsstat = tmppsstat[tmppsstat.index.isin(['asm'], level='use')]
        tmpassim.append(tmppsstat['count'].sum())
        # sum it up
        obsname.append(conv)
        nassim.append(np.sum(np.array(tmpassim)))
        nread.append(np.sum(np.array(tmpread)))
        nkeep.append(np.sum(np.array(tmpread)))

    # sort everything by number of obs read
    yx = sorted(zip(nread, obsname))
    x_sorted = [x for y, x in yx]
    nread_sorted = [y for y, x in yx]
    yx = sorted(zip(nread, nassim))
    nassim_sorted = [x for y, x in yx]
    x_pos = [i for i, _ in enumerate(x_sorted)]
    # read obs
    plt.figure(figsize=(8, 11))
    plt.barh(x_pos, nread_sorted, height=1, log=True,
             edgecolor='black', color='firebrick')
    plt.xlabel('# of Observations Read')
    plt.title('GDAS Input Observation Counts ' +
              cycle.strftime('%Y-%m-%d %H:%MZ Cycle'))
    plt.yticks(x_pos, x_sorted)
    for i, v in enumerate(nread_sorted):
        plt.text(v * 1.5, i-.15, str(int(v)), color='black')
    plt.savefig(plotdir+'/nobs_read_'+cycle.strftime('%Y%m%d%H')+'.png')
    # assimilated obs
    plt.figure(figsize=(8, 11))
    plt.barh(x_pos, nassim_sorted, height=1, log=True,
             edgecolor='black', color='gray')
    plt.xlabel('# of Observations Assimilated')
    plt.title('GDAS Assimilated Observation Counts ' +
              cycle.strftime('%Y-%m-%d %H:%MZ Cycle'))
    plt.yticks(x_pos, x_sorted)
    for i, v in enumerate(nassim_sorted):
        plt.text(v * 1.5, i-.15, str(int(v)), color='black')
    plt.savefig(plotdir+'/nobs_assim_'+cycle.strftime('%Y%m%d%H')+'.png')
    # percentage of read obs assimilated
    plt.figure(figsize=(8, 11))
    plt.barh(x_pos, (np.array(nassim_sorted)/np.array(nread_sorted))
             * 100., height=1, edgecolor='black', color='firebrick')
    plt.xlabel('Percentage of Observations Assimilated')
    plt.title('GDAS Observations Assimilated Percentage ' +
              cycle.strftime('%Y-%m-%d %H:%MZ Cycle'))
    plt.yticks(x_pos, x_sorted)

    cond = (np.array(nassim_sorted)/np.array(nread_sorted))*100.
    for i, v in enumerate(cond):
        if not np.isinf(v):
            plt.text(v + 1, i-.15, str(round(v, 2)), color='black')
    plt.savefig(plotdir+'/nobs_pct_'+cycle.strftime('%Y%m%d%H')+'.png')
    # read and assimilated
    plt.figure(figsize=(8, 11))
    plt.barh(x_pos, nassim_sorted, height=0.5,
             log=True, edgecolor='black', color='gray')
    plt.barh(np.array(x_pos)+0.5, nread_sorted, height=0.5,
             log=True, edgecolor='black', color='firebrick')
    plt.xlabel('# of Observations Read/Assimilated')
    plt.title('GDAS Input and Assimilated Observation Counts ' +
              cycle.strftime('%Y-%m-%d %H:%MZ Cycle'))
    plt.yticks(np.array(x_pos)+0.25, x_sorted)
    for i, v in enumerate(nread_sorted):
        plt.text(v * 1.5, i+0.33, str(int(v)), color='black')
    for i, v in enumerate(nassim_sorted):
        plt.text(v * 1.5, i-0.20, str(int(v)), color='black')
    plt.savefig(plotdir+'/nobs_read_assim_'+cycle.strftime('%Y%m%d%H')+'.png')


if __name__ == '__main__':
    # called from command line
    ap = argparse.ArgumentParser()
    ap.add_argument("-i", "--input",
                    help="Path to GSI stat file", required=True)
    ap.add_argument("-c", "--cycle",
                    help="YYYYMMDDHH analysis cycle", required=True)
    ap.add_argument("-p", "--plotdir",
                    help="Path to output plot dir", default='./')
    MyArgs = ap.parse_args()

    plot_obscounts(MyArgs.input, MyArgs.cycle, MyArgs.plotdir)
