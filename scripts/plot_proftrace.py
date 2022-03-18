import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from pyGSI.gsi_stat import GSIstat
import os

def get_datelist(startdate, enddate):
    import datetime
    datelist = []
    mydate = startdate
    while mydate <= enddate:
        date = datetime.datetime.strptime(mydate,'%Y%m%d%H')
        datelist.append(mydate)
        nextdate = date + datetime.timedelta(hours=6)
        mydate = datetime.datetime.strftime(nextdate,'%Y%m%d%H')
    return datelist

########################
#
# USER DEFINED SETTINGS
#
# Each member of the plot is provided as a separate dictionary to setdict with the following keys:
#    name: Name of member, used in the legend of the plot to distinguish lines or bars belonging to this set
#    var: Variable being extracted ('uv', 't', or 'q')
#    it: Iteration(s) (1, 2, or 3), can be an individual value or a list of values (e.g. [1,3]), or None to collect all values
#    use: Use-type ('asm', 'mon', or 'rej'), can be an individual value or a list of values (e.g. ['asm','rej']), or None to collect all values
#    typ: Observation type number, can be an individual value or a list of values, or None to collect all values
#    styp: Observation subtype string, can be an individual value or a list of values, or None to collect all values
#    statdir: String identifying full-path to directory containing gsistat files (ends in '/')
#
# You can define the list of cycles to extract data from, which will be aggregated for profiles and plotted individually for traces by:
#    1) Specifying a list of cycles as '%Y%m%d%H' formatted strings, e.g.:
#        cycle=['2020080100','2020080106','2020080112']
#    2) Use get_datelist() to automatically generate a list of cycles at 6-hourly steps between a beginning and ending date, e.g.:
#        cycle=get_datelist('2020080100','2020080112')
#
# Figures are saved to profs_filename and trace_filename, respectively
#
# tskip defines the number of time-periods in the trace that are skipped when adding tick-labels. Tick-labels are dates of the
# form (e.g.) 'Aug01-00z', so for long timeseries traces the labels will be written over each other unless some of them are
# skipped. I usually use some multiple of 4, so that the tick-labels always correspond to the same analysis-period.
#
########################
#
setdict=[]
setdict.append({ 'name' : 'RAOB-ges', 
                 'var' : 'uv', 
                 'it' : 1, 
                 'use' : 'asm', 
                 'typ' : 220, 
                 'styp' : None,
                 'statdir' : '<full-path-to-gsistat-files>/'
               })
setdict.append({ 'name' : 'RAOB-anl', 
                 'var' : 'uv', 
                 'it' : 3, 
                 'use' : 'asm', 
                 'typ' : 220, 
                 'styp' : None, 
                 'statdir' : '<full-path-to-gsistat-files>/'
               })
cycles=get_datelist('2020090800', '2020093018')
profs_filename = 'gsistat_profs_RAOBS.png'
trace_filename = 'gsistat_trace_RAOBS.png'
tskip = 8
#
########################

rmses = {}
counts = {}
biases = {}
levels={}
for sd in setdict:
    setname = sd['name']
    statdir = sd['statdir']
    rmses[setname] = {}
    counts[setname] = {}
    biases[setname] = {}
    for cycle in cycles:
        rmses[setname][cycle] = {}
        counts[setname][cycle] = {}
        biases[setname][cycle] = {}
        inputfile = os.path.join(statdir, 'gsistat.gdas.'+cycle)
        try:
            gdas = GSIstat(inputfile, cycle)
        except FileNotFoundError:
            raise FileNotFoundError(f'Unable to find {inputfile} for cycle {cycle}')
        # statistics are drawn for a single variable
        var=sd['var']
        # all other filters (it,use,typ,styp) could contain multiple entries in a list
        # if they are not a list, assert them as a list
        it=sd['it']
        if type(it) != list:
            it = [it]
        use=sd['use']
        if type(use) != list:
            use = [use]
        typ=sd['typ']
        if type(typ) != list:
            typ = [typ]
        styp=sd['styp']
        if type(styp) != list:
            styp = [styp]
        # extract variable from gdas
        stat = gdas.extract(var)  # t, uv, q, etc.
        #                date         it         obs   use    typ  styp    stat
        # pull entire dataframe
        s = stat.loc[(slice(None), slice(None), slice(None), slice(None), slice(None), slice(None), slice(None)),:]
        # Filters:     date         it           obs         use          typ          styp         stat
        # Add try/except to filters: If any filter breaks, return 0-values for rmse, bias, and count
        try:
            if (it[0] != None):
                s = s[s.index.isin(it, level='it')]
            if (use[0] != None):
                s = s[s.index.isin(use, level='use')]
            if (typ[0] != None):
                s = s[s.index.isin(typ, level='typ')]
            if (styp[0] != None):
                s = s[s.index.isin(styp, level='styp')]
        except:
            # Pull entire variable frame again
            s = stat.loc[(slice(None), slice(None), slice(None), slice(None), slice(None), slice(None), slice(None)),:]   
            # Filter to just the first observation type in the frame (index=[0:3] is count, rms, and bias of first entry)
            s = s.iloc[0:3, :]
            # Zero-out all elements
            s.loc[(slice(None), slice(None), slice(None), slice(None), slice(None), slice(None), slice(None)),s.columns]=0.
        # If successfully passed filters but has zero rows, return 0-values for rmse, bias, and count
        if s.shape[0] == 0:
            # Pull entire variable frame again
            s = stat.loc[(slice(None), slice(None), slice(None), slice(None), slice(None), slice(None), slice(None)),:]   
            # Filter to just the first observation type in the frame (index=[0:3] is count, rms, and bias of first entry)
            s = s.iloc[0:3, :]
            # Zero-out all elements
            s.loc[(slice(None), slice(None), slice(None), slice(None), slice(None), slice(None), slice(None)),s.columns]=0.
        # Levels are contained in all but the final value, which is a statistic for the full column
        levs=[]
        levkeys=list(stat.keys())
        for l in levkeys[:-1]:
            levs.append(float(l))
        levels[setname] = levs
        rmses[setname][cycle]['lev'] = s.loc[(slice(None), slice(None), slice(None), slice(None), slice(None), slice(None), 'rms'),s.columns[:-1]]
        counts[setname][cycle]['lev'] = s.loc[(slice(None), slice(None), slice(None), slice(None), slice(None), slice(None), 'count'),s.columns[:-1]]
        biases[setname][cycle]['lev'] = s.loc[(slice(None), slice(None), slice(None), slice(None), slice(None), slice(None), 'bias'),s.columns[:-1]]
        rmses[setname][cycle]['col'] = s.loc[(slice(None), slice(None), slice(None), slice(None), slice(None), slice(None), 'rms'),s.columns[-1]]
        counts[setname][cycle]['col'] = s.loc[(slice(None), slice(None), slice(None), slice(None), slice(None), slice(None), 'count'),s.columns[-1]]
        biases[setname][cycle]['col'] = s.loc[(slice(None), slice(None), slice(None), slice(None), slice(None), slice(None), 'bias'),s.columns[-1]]

def plot_stat_profiles(rmsList,biasList,countList,nameList,plevList,colMap=['tab10', 10]):
    ######################################################################################################
    # Generates 2-panel plot:
    #    Left: Profiles of rms and bias for each set
    #    Right: Horizontal bar-chart of ob-counts for each set
    #
    # INPUTS
    #    rmsList: list of numpy arrays (nlev,) of rms by pressure level for each set
    #    biasList: list of numpy arrays (nlev,) of bias by pressure level for each set
    #    countList: list of numpy arrays (nlev,) of ob-count by pressure level for each set
    #    nameList: list of names for each set (for figure legend)
    #    plevList: list of profile pressure levels (str or float)
    #    colMap: 2-element list containing colormap and colormap-range for panels (default: ['tab10',10]))
    # OUTPUTS
    #    plotFig: plot figure, as plt.fig()
    #
    # DEPENDENCIES
    # numpy
    # matplotlib.rc
    # matplotlib.pyplot
    # matplotlib.cm
    ######################################################################################################
    import numpy as np
    import matplotlib
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    #
    # Set font size and type
    #
    font = {'family' : 'DejaVu Sans',
            'weight' : 'bold',
            'size'   : 22}

    matplotlib.rc('font', **font)
    #
    # Identify number of profile sets (should be indentical for rmsList, biasList, countList)
    #
    n_profiles = len(rmsList)
    #
    # Define colormap: The default settings are to select a range of 10 on 'tab10', so that 10 pairs of 
    # rms/bias profiles can be produced for the left panel that correspond to 10 ob-counts on the right 
    # panel. The user can select a different colormap and range with colMapLeft and colMapRight options. 
    # If you want to sample the entire colorbar, set the range to the number of datasets plotted.
    #
    # scalarMapList is used to select colors for each profile/bar
    scalarMap=cm.ScalarMappable(cmap=colMap[0])
    scalarMapList=scalarMap.to_rgba(range(colMap[1]))
    #
    # Generate figure
    plt.figure(figsize=(18, 18))
    # Define offset values (needed to define y-axis limits on both plots, for consistency)
    y_offset=8.0*(np.arange(n_profiles)-np.floor(n_profiles/2))
    # LEFT PANEL: rms and bias scores
    plt.subplot(121)
    # For each set, plot rms and bias: plot both lines and circular markers
    legend_list=[]
    for i in range(n_profiles):
        rms=rmsList[i]
        bias=biasList[i]
        levs=plevList[i]
        n_levs=np.size(levs)
        # Define y-axis limits
        y_min = min(levs)+min(y_offset)-8.0
        y_max = max(levs)+max(y_offset)+8.0
        # rms profile
        legend_list.append(nameList[i]+' rms')
        prof_color=list(scalarMapList[i][0:3])
        plt.plot(rms, levs, color=prof_color, linewidth=3)
        plt.plot(rms, levs, 'o',color=prof_color, markersize=8, label='_nolegend_')
        # bias profile
        legend_list.append(nameList[i]+' bias')
        prof_color=list(scalarMapList[i][0:3])
        plt.plot(bias, levs, color=prof_color, linewidth=3, linestyle='dashdot')
        plt.plot(bias, levs, 'o', color=prof_color, markersize=8, label='_nolegend_')
    # Zero-line
    plt.plot(np.zeros((n_levs, )), levs, color='k', linewidth=1, linestyle='dashed', label='_nolegend_')
    # Set y-limits to entire range
    plt.ylim((y_min, y_max))
    plt.yticks(levs)
    # Reverse y-axis, if levs is in descending-order (often the case with pressure coordinate data)
    if (levs[1]<levs[0]):
        plt.gca().invert_yaxis()
    # Set legend
    plt.legend(legend_list, frameon=False, fontsize=10)
    # Set x-label
    plt.xlabel('RMS or Bias')
    # RIGHT PANEL: ob-counts
    plt.subplot(122)
    # For each set, plot ob-counts and generate legend list
    # Counts are asserted to be in thousands
    legend_list=[]
    for i in range(n_profiles):
        count=0.001*countList[i]
        bar_color=list(scalarMapList[i][0:3])
        plt.barh(levs+y_offset[i], count, height=8.0, color=bar_color)
    # Set y-limits to entire range
    plt.ylim((y_min, y_max))
    plt.yticks(levs)
    # Reverse y-axis, if levs is in descending-order (often the case with pressure coordinate data)
    if (levs[1]<levs[0]):
        plt.gca().invert_yaxis()
    # Set legend
    plt.legend(nameList, frameon=False, fontsize=10)
    # Set x-label
    plt.xlabel('Ob Count (Thousands)')
    # Turn off interactive-mode to suppress plotting figure
    plt.ioff()
    # Return
    return plt.gcf()

def plot_stat_traces(rmsList, biasList, countList, nameList, dateList, tskip=4, colMap=['tab10', 10]):
    ######################################################################################################
    # Generates 2-panel plot:
    #    Top: Trace of rms and bias for each set
    #    Bottom: Bar-chart of ob-counts for each set
    #
    # INPUTS
    #    rmsList: list of numpy arrays (nlev,) of rms by pressure level for each set
    #    biasList: list of numpy arrays (nlev,) of bias by pressure level for each set
    #    countList: list of numpy arrays (nlev,) of ob-count by pressure level for each set
    #    nameList: list of names for each set (for figure legend)
    #    dateList: list of dates (str in '%Y%m%d%H' format)
    #    tskip: number of ticks to skip when writing tick-labels (default: 4)
    #    colMap: 2-element list containing colormap and colormap-range for panels (default: ['tab10',10]))
    # OUTPUTS
    #    plotFig: plot figure, as plt.fig()
    #
    # DEPENDENCIES
    # numpy
    # datetime
    # matplotlib.rc
    # matplotlib.pyplot
    # matplotlib.cm
    ######################################################################################################
    import numpy as np
    from datetime import datetime
    import matplotlib
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    #
    # Set font size and type
    #
    font = {'family' : 'DejaVu Sans',
            'weight' : 'bold',
            'size'   : 22}

    matplotlib.rc('font', **font)
    #
    # Identify number of trace sets (should be indentical for rmsList, biasList, countList)
    #
    n_trace = len(rmsList)
    #
    # Define colormap: The default settings are to select a range of 10 on 'tab10', so that 10 pairs of 
    # rms/bias profiles can be produced for the left panel that correspond to 10 ob-counts on the right 
    # panel. The user can select a different colormap and range with colMapLeft and colMapRight options. 
    # If you want to sample the entire colorbar, set the range to the number of datasets plotted.
    #
    # scalarMapList is used to select colors for each profile/bar
    scalarMap=cm.ScalarMappable(cmap=colMap[0])
    scalarMapList=scalarMap.to_rgba(range(colMap[1]))
    #
    # Generate figure
    plt.figure(figsize=(18, 18))
    # Define offset values (needed to define x-axis limits on both plots, for consistency)
    offs=0.8/n_trace
    x_offset=offs*(np.arange(n_trace)-np.floor(n_trace/2))
    # TOP PANEL: rms and bias traces
    plt.subplot(211)
    # For each set, plot rms and bias: plot both lines and circular markers
    legend_list=[]
    for i in range(n_trace):
        rms=rmsList[i]
        bias=biasList[i]
        dates=dateList[i]
        # Convert dates from %Y%m%d%H format to %b%d:%HZ format for tick labels
        dstr=[]
        for d in dates:
            dt=datetime.strptime(d, '%Y%m%d%H')
            dstr.append(datetime.strftime(dt, '%b%d:%HZ'))
        n_dates=np.size(dates)
        # Define x-axis limits
        x_min = min(x_offset)-offs
        x_max = n_dates+max(x_offset)+offs
        x_rng = np.arange(1, n_dates+1.0E-05)
        # rms trace
        legend_list.append(nameList[i]+' rms')
        prof_color=list(scalarMapList[i][0:3])
        plt.plot(x_rng, rms, color=prof_color, linewidth=3)
        plt.plot(x_rng, rms, 'o', color=prof_color, markersize=8, label='_nolegend_')
        # bias profile
        legend_list.append(nameList[i]+' bias')
        prof_color=list(scalarMapList[i][0:3])
        plt.plot(x_rng, bias, color=prof_color, linewidth=3, linestyle='dashdot')
        plt.plot(x_rng, bias, 'o', color=prof_color, markersize=8, label='_nolegend_')
    # Zero-line
    plt.plot(x_rng, np.zeros((n_dates, )), color='k', linewidth=1, linestyle='dashed', label='_nolegend_')
    # Set x-limits to entire range
    plt.xlim((x_min, x_max))
    plt.xticks(ticks=x_rng[::tskip], labels=dstr[::tskip], fontsize=10)
    # Set legend
    plt.legend(legend_list, frameon=False, fontsize=10)
    # Set y-label
    plt.ylabel('RMS or Bias')
    # BOTTOM PANEL: ob-counts
    plt.subplot(212)
    # For this plot, bars need to be offset from each other so that they all cluster around the pressure
    # level. This is accomplished with an offset value added to each bar
    legend_list=[]
    for i in range(n_trace):
        count=0.001*countList[i]
        bar_color=list(scalarMapList[i][0:3])
        plt.bar(x_rng+x_offset[i], count, width=offs, color=bar_color)
    # Set x-limits to entire range
    plt.xlim((x_min, x_max))
    plt.xticks(ticks=x_rng[::tskip], labels=dstr[::tskip], fontsize=10)
    # Set legend
    plt.legend(nameList, frameon=False, fontsize=10)
    # Set y-label
    plt.ylabel('Ob Count (Thousands)')
    # Turn off interactive-mode to suppress plotting figure
    plt.ioff()
    # Return
    return plt.gcf()

#
# Generate plotting data by aggregating rmse, bias, and count data across times (for profiles), or
# aggregate full-column data (for traces)
#
rmse_profs=[]
rmse_trace=[]
bias_profs=[]
bias_trace=[]
nobs_profs=[]
nobs_trace=[]
levl_profs=[]
date_trace=[]
name_list=[]
for setname in rmses.keys():
    name_list.append(setname)
    plevs=levels[setname]
    levl_profs.append(np.asarray(plevs).squeeze())
    nlev=np.size(plevs)
    rmse=np.zeros((nlev, ))
    bias=np.zeros((nlev, ))
    nobs=np.zeros((nlev, ))
    dates=list(rmses[setname].keys())
    date_trace.append(dates)
    ndate=len(dates)
    r_trace=np.zeros((ndate, ))
    b_trace=np.zeros((ndate, ))
    n_trace=np.zeros((ndate, ))
    for i in range(ndate):
        date=dates[i]
        n_ele = np.size(rmses[setname][date]['col'].values)
        for j in range(n_ele):
            # Aggregate rms, bias, and count for full-column values at each time to define traces
            rcol=rmses[setname][date]['col'].values[j]
            bcol=biases[setname][date]['col'].values[j]
            ncol=counts[setname][date]['col'].values[j]
            ncol2 = n_trace[i] + ncol
            if ncol2 > 0:
                rcol2 = np.sqrt((n_trace[i] * r_trace[i] ** 2. + ncol * rcol ** 2.) / ncol2)
                bcol2 = (n_trace[i] * b_trace[i] + ncol * bcol) / ncol2
            else:
                rcol2 = 0.
                bcol2 = 0.
            r_trace[i] = rcol2
            b_trace[i] = bcol2
            n_trace[i] = ncol2
            # Aggregate rms, bias, and count at each level across all times to define profiles
            r=rmses[setname][date]['lev'].values[j]
            b=biases[setname][date]['lev'].values[j]
            c=counts[setname][date]['lev'].values[j]
            for k in range(nlev):
                c2 = nobs[k] + c[k]
                if c2 > 0:
                    r2 = np.sqrt((nobs[k] * rmse[k] ** 2. + c[k] * r[k] ** 2.) / c2)
                    b2 = (nobs[k] * bias[k] + c[k] * b[k]) / c2
                else:
                    r2 = 0.
                    b2 = 0.
                rmse[k] = r2
                bias[k] = b2
                nobs[k] = c2
    # Any score with 0 obs should be changed to NaN before entering into record
    rmse[nobs == 0] = np.nan
    bias[nobs == 0] = np.nan
    r_trace[n_trace == 0] = np.nan
    b_trace[n_trace == 0] = np.nan
    rmse_trace.append(r_trace)
    bias_trace.append(b_trace)
    nobs_trace.append(n_trace)
    rmse_profs.append(rmse)
    bias_profs.append(bias)
    nobs_profs.append(nobs)
#
# Generate profile plot
#
fig_prof=plot_stat_profiles(rmse_profs, bias_profs, nobs_profs, name_list, levl_profs)
plt.ioff()
fig_prof.savefig(profs_filename, bbox_inches='tight', facecolor='w')
#
# Generate trace plot
#
fig_trace=plot_stat_traces(rmse_trace, bias_trace, nobs_trace, name_list, date_trace, tskip=tskip)
plt.ioff()
fig_trace.savefig(trace_filename, bbox_inches='tight', facecolor='w')


