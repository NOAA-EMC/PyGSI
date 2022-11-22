##########################################################################
#
# YAML FILE CONTENTS
#
# yaml file contents come in three types, based on the key name:
#
# figures: Defines settings for profile and trace figures
#    pname: Name of profile figure-file
#    tname: Name of trace figure-file
#    tskip: Number of ticks that are skipped when plotting tick-labels on
#           trace-figure (should be large enough to keep labels from
#           stepping on each other)
#
# dates: Defines settings for time-period over which statistics are
#        collected and plotted
#    dt_beg: Beginning date of time-series (%Y%m%d%H format)
#    dt_end: Ending date of time-series (%Y%m%d%H format)
#    hrs: Number of hours between each time-period between dt_beg and
#         dt_end (usually 6 for collecting every cycle)
#
# <any other key name>: Each member of the plot is provided as a separate
# dictionary with the following keys:
#    name: Name of member, used in the legend of the plot to distinguish
#          lines or bars belonging to this set
#    var: Variable being extracted ('uv', 't', or 'q')
#    it: Iteration(s) (1, 2, or 3), can be an individual value or a list
#        of values (e.g. [1,3]), or None to collect all values
#    use: Use-type ('asm', 'mon', or 'rej'), can be an individual value
#         or a list of values (e.g. ['asm','rej']), or None to collect all
#         values
#    typ: Observation type number, can be an individual value or a list of
#         values, or None to collect all values
#    styp: Observation subtype string, can be an individual value or a
#          list of values, or None to collect all values
#    statdir: String identifying full-path to directory containing gsistat
#             files (ends in '/')
##########################################################################
# The list of cycles to extract data from, which will be aggregated for
# profiles and plotted individually for traces, is defined by:
#    1) Use of get_datelist() to automatically generate a list of cycles
#       at regular steps between a beginning and ending date, e.g.:
#
#       cycles=get_datelist(dt_beg,dt_end,hrs)
#
#    2) During search for gsistat files for each cycle, if a gsistat file
#       is missing, the date is removed from cycles
#    3) During collection of statistics, if there are no valid statistics
#       for a given cycle, NaN values are returned and plotted
#
# Figures are saved to pname and tname, respectively
#
# tskip defines the number of time-periods in the trace that are skipped
# when adding tick-labels. Tick-labels are dates of the form (e.g.)
# 'Aug01-00z', so for long timeseries traces the labels will be written
# over each other unless some of them are skipped. I usually use some
# multiple of 4, so that the tick-labels always correspond to the same
# analysis-period.
##########################################################################
def get_datelist(startdate, enddate, hrs):
    import datetime
    datelist = []
    mydate = startdate
    while mydate <= enddate:
        date = datetime.datetime.strptime(mydate, '%Y%m%d%H')
        datelist.append(mydate)
        nextdate = date + datetime.timedelta(hours=hrs)
        mydate = datetime.datetime.strftime(nextdate, '%Y%m%d%H')
    return datelist


def parse_yaml(yaml_file):
    import yaml
    # YAML entries come in three types:
    # key='dates' : provides beginning (dt_beg)/ending (dt_end) dates
    #               for get_datelist in %Y%m%d%H format, and
    #               timedelta (hrs)
    # key='figures' : provides figure filenames for profile figure
    #                 (pname) and trace figure (tname), and number of
    #                 times to skip on tick-labels (tskip)
    # any other key: provides filtering dictionaries for a dataset to
    #                be plotted
    with open(yaml_file, 'r') as stream:
        try:
            parsed_yaml = yaml.safe_load(stream)
        except yaml.YAMLError as YAMLError:
            parsed_yaml = None
            print(f'YAMLError: {YAMLError}')
        if parsed_yaml is not None:
            # Extract 'dates' data
            try:
                dt_beg = parsed_yaml['dates']['dt_beg']
                dt_end = parsed_yaml['dates']['dt_end']
                hrs = parsed_yaml['dates']['hrs']
            except KeyError as MissingDatesError:
                dt_beg = None
                dt_end = None
                hrs = None
                print(f'MissingDatesError: {MissingDatesError}')
            # Extract 'figures' data
            try:
                pname = parsed_yaml['figures']['pname']
                tname = parsed_yaml['figures']['tname']
                tskip = parsed_yaml['figures']['tskip']
            except KeyError as MissingFiguresError:
                pname = None
                tname = None
                tskip = None
                print(f'MissingFiguresError: {MissingFiguresError}')
            # Extract all other keys as filtering dictionaries, store
            # in setdict list
            setdict = []
            fdicts = {x: parsed_yaml[x] for x in parsed_yaml
                      if x not in {'dates', 'figures'}}
            fkeys = ['it', 'use', 'typ', 'styp']
            # Check filter keys for 'None' and change to None as
            # appropriate before storing
            for fd in fdicts.keys():
                fdict = fdicts[fd]
                for key in fkeys:
                    try:
                        fdict[key] = None if fdict[key] == 'None' else\
                                     fdict[key]
                    except KeyError as FilterKeyError:
                        print(f'FilterKeyError: {FilterKeyError}')
                setdict.append(fdict)
    return dt_beg, dt_end, hrs, pname, tname, tskip, setdict


def collect_statistics(setdict):
    ######################################################################
    #
    # Collect rmse, bias, and ob-count statistics for each data (sub)set
    # in setdict. Statistics are provided for each individual pressure-bin
    # as well as a full-column value.
    #
    # If a gsistat file is missing for a given cycle, it raises a
    # FileNotFoundError and the cycle is skipped, not appearing at all in
    # the time-series
    #
    # If a gsistat file exists for a given cycle but no data for a data
    # (sub)set is present, the rmse, bias, and ob-count are set to zero but
    # the cycle is not skipped. Statistics with zero ob-count are
    # reassigned to NaN later.
    #
    # INPUTS:
    #    setdict: List of dictionaries for each subset, including name,
    #             full-path to directory containing gsistat files,
    #             selected variable, and subsetting filters
    #
    # OUTPUTS:
    #    rmses: rms-difference between subset of obs and model
    #           (unadjusted background, adusted background, or
    #           analysis, depending on value of it) for each cycle
    #           where a gsistat file exists
    #    biases: average difference between subset of obs and model
    #            for each cycle where a gsistat file exists
    #    counts: number of subset obs  for each cycle where a gsistat
    #            file exists
    #    levels: pressure levels for each subset of obs
    #
    # DEPENDENCIES:
    #    os
    #    pyGSI.gsi_stat.GSIstat
    #    pandas
    ######################################################################
    import os
    from pyGSI.gsi_stat import GSIstat
    import pandas
    # Initialize outputs as empty dictionaries
    rmses = {}
    counts = {}
    biases = {}
    levels = {}
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
            # Test if inputfile exists, if not, set all values 
            # for rmses/counts/biases to None for this setname
            # and cycle for both 'lev' and 'col' values
            if os.path.exists(inputfile):
                gdas = GSIstat(inputfile, cycle)
            else:
                rmses[setname][cycle]['lev'] = None
                counts[setname][cycle]['lev'] = None
                biases[setname][cycle]['lev'] = None
                rmses[setname][cycle]['col'] = None
                counts[setname][cycle]['col'] = None
                biases[setname][cycle]['col'] = None
                # Skip rest of this for-loop entry
                continue
            # statistics are drawn for a single variable
            var = sd['var']
            # all other filters (it,use,typ,styp) could contain multiple
            # entries in a list
            # if they are not a list, assert them as a list
            it = sd['it']
            it = [it] if type(it) != list else it
            use = sd['use']
            use = [use] if type(use) != list else use
            typ = sd['typ']
            typ = [typ] if type(typ) != list else typ
            styp = sd['styp']
            styp = [styp] if type(styp) != list else styp
            # extract variable from gdas
            stat = gdas.extract(var)  # t, uv, q, etc.
            #             date         it           obs          use
            #             typ          styp         stat
            # pull entire dataframe
            s = stat.loc[(slice(None), slice(None), slice(None), slice(None),
                          slice(None), slice(None), slice(None)), :]
            # Filters:    date         it           obs          use
            #             typ          styp         stat
            # Add try/except to filters: If any filter breaks, return 0-values
            # for rmse, bias, and count
            try:
                if (it[0] is not None):
                    s = s[s.index.isin(it, level='it')]
                if (use[0] is not None):
                    s = s[s.index.isin(use, level='use')]
                if (typ[0] is not None):
                    s = s[s.index.isin(typ, level='typ')]
                if (styp[0] is not None):
                    s = s[s.index.isin(styp, level='styp')]
            except Exception as FilterError:
                print(f'FilterError {FilterError} for {cycle}')
                # Pull entire variable frame again
                s = stat.loc[(slice(None), slice(None), slice(None),
                              slice(None), slice(None), slice(None),
                              slice(None)), :]
                # Filter to just the first observation type in the frame
                # (index=[0:3] is count, rms, and bias of first entry)
                s = s.iloc[0:3, :]
                # Zero-out all elements
                s.loc[(slice(None), slice(None), slice(None),
                       slice(None), slice(None), slice(None),
                       slice(None)), s.columns] = 0.
            # If successfully passed filters but has zero rows, return 0-values
            # for rmse, bias, and count
            if s.shape[0] == 0:
                # Pull entire variable frame again
                s = stat.loc[(slice(None), slice(None), slice(None),
                              slice(None), slice(None), slice(None),
                              slice(None)), :]
                # Filter to just the first observation type in the frame
                # (index=[0:3] is count, rms, and bias of first entry)
                s = s.iloc[0:3, :]
                # Zero-out all elements
                s.loc[(slice(None), slice(None), slice(None),
                       slice(None), slice(None), slice(None),
                       slice(None)), s.columns] = 0.
            # Levels are contained in all but the final value, which is a
            # statistic for the full column
            levs = []
            levkeys = list(stat.keys())
            for lk in levkeys[:-1]:
                levs.append(float(lk))
            levels[setname] = levs
            rmses[setname][cycle]['lev'] = s.loc[(slice(None), slice(None),
                                                  slice(None), slice(None),
                                                  slice(None), slice(None),
                                                  'rms'), s.columns[:-1]]
            counts[setname][cycle]['lev'] = s.loc[(slice(None), slice(None),
                                                   slice(None), slice(None),
                                                   slice(None), slice(None),
                                                   'count'), s.columns[:-1]]
            biases[setname][cycle]['lev'] = s.loc[(slice(None), slice(None),
                                                   slice(None), slice(None),
                                                   slice(None), slice(None),
                                                   'bias'), s.columns[:-1]]
            rmses[setname][cycle]['col'] = s.loc[(slice(None), slice(None),
                                                  slice(None), slice(None),
                                                  slice(None), slice(None),
                                                  'rms'), s.columns[-1]]
            counts[setname][cycle]['col'] = s.loc[(slice(None), slice(None),
                                                   slice(None), slice(None),
                                                   slice(None), slice(None),
                                                   'count'), s.columns[-1]]
            biases[setname][cycle]['col'] = s.loc[(slice(None), slice(None),
                                                   slice(None), slice(None),
                                                   slice(None), slice(None),
                                                   'bias'), s.columns[-1]]
    # Return output dictionaries
    return rmses, biases, counts, levels


def aggregate_figure_data(rmses, biases, counts, levs):
    ######################################################################
    #
    # Generate plotting data by aggregating rmse, bias, and count data
    # across times (for profiles), or aggregate full-column data (for
    # traces). Returned lists of profile and trace data serve as inputs to
    # plotting functions.
    #
    # INPUTS:
    #    rmses: Dictionary of rmse data for each subset
    #    biases: Dictionary of bias data for each subset
    #    counts: Dictionary of ob-count data for each subset
    #    levs: Dictionary of pressure-bin level data for each subset
    #
    # OUTPUTS:
    #    rmse_profs: List of rmse profiles to plot
    #    rmse_trace: List of rmse time-series traces to plot
    #    bias_profs: List of bias profiles to plot
    #    bias_trace: List of bias time-series traces to plot
    #    nobs_profs: List of ob-count profiles to plot
    #    nobs_trace: List of ob-count time-series traces to plot
    #    levl_profs: List of pressure-bin levels for profiles
    #    date_trace: List of dates for time-series traces
    #    name_list: List of dataset names for plot legends
    #
    # DEPENDENCIES:
    #    numpy
    #    pandas
    ######################################################################
    import numpy as np
    import pandas
    # Initialize outputs as empty lists
    rmse_profs = []
    rmse_trace = []
    bias_profs = []
    bias_trace = []
    nobs_profs = []
    nobs_trace = []
    levl_profs = []
    date_trace = []
    name_list = []
    # Loop through each set
    for setname in rmses.keys():
        name_list.append(setname)
        plevs = levels[setname]
        levl_profs.append(np.asarray(plevs).squeeze())
        nlev = np.size(plevs)
        rmse = np.zeros((nlev, ))
        bias = np.zeros((nlev, ))
        nobs = np.zeros((nlev, ))
        dates = list(rmses[setname].keys())
        date_trace.append(dates)
        ndate = len(dates)
        r_trace = np.zeros((ndate, ))
        b_trace = np.zeros((ndate, ))
        n_trace = np.zeros((ndate, ))
        for i in range(ndate):
            date = dates[i]
            # If rmses[setname][date]['col] is NoneType, set n_ele=0
            if type(rmses[setname][date]['col'])==type(None):
                n_ele=0
            # If not NoneType, define n_ele by values in dataframe
            else:
                n_ele = np.size(rmses[setname][date]['col'].values)
            for j in range(n_ele):
                # Aggregate rms, bias, and count for full-column values at each
                # time to define traces
                rcol = rmses[setname][date]['col'].values[j]
                bcol = biases[setname][date]['col'].values[j]
                ncol = counts[setname][date]['col'].values[j]
                ncol2 = n_trace[i] + ncol
                if ncol2 > 0:
                    rcol2 = np.sqrt((n_trace[i] * r_trace[i] ** 2. + ncol *
                                     rcol ** 2.) / ncol2)
                    bcol2 = (n_trace[i] * b_trace[i] + ncol * bcol) / ncol2
                else:
                    rcol2 = 0.
                    bcol2 = 0.
                r_trace[i] = rcol2
                b_trace[i] = bcol2
                n_trace[i] = ncol2
                # Aggregate rms, bias, and count at each level across all times
                # to define profiles
                r = rmses[setname][date]['lev'].values[j]
                b = biases[setname][date]['lev'].values[j]
                c = counts[setname][date]['lev'].values[j]
                for k in range(nlev):
                    c2 = nobs[k] + c[k]
                    if c2 > 0:
                        r2 = np.sqrt((nobs[k] * rmse[k] ** 2. +
                                     c[k] * r[k] ** 2.) / c2)
                        b2 = (nobs[k] * bias[k] + c[k] * b[k]) / c2
                    else:
                        r2 = 0.
                        b2 = 0.
                    rmse[k] = r2
                    bias[k] = b2
                    nobs[k] = c2
        # Any score with 0 obs should be changed to NaN before entering into
        # record
        rmse[nobs == 0] = np.nan
        bias[nobs == 0] = np.nan
        r_trace[n_trace == 0] = np.nan
        b_trace[n_trace == 0] = np.nan
        # Append profile and trace data to output lists
        rmse_trace.append(r_trace)
        bias_trace.append(b_trace)
        nobs_trace.append(n_trace)
        rmse_profs.append(rmse)
        bias_profs.append(bias)
        nobs_profs.append(nobs)
    # Return output lists
    return (rmse_profs, rmse_trace, bias_profs, bias_trace, nobs_profs,
            nobs_trace, levl_profs, date_trace, name_list)


def plot_stat_profiles(rmsList, biasList, countList, nameList, plevList,
                       colMap=['tab10', 10]):
    ######################################################################
    # Generates 2-panel plot:
    #    Left: Profiles of rms and bias for each set
    #    Right: Horizontal bar-chart of ob-counts for each set
    #
    # INPUTS
    #    rmsList: list of numpy arrays (nlev,) of rms by pressure level
    #             for each set
    #    biasList: list of numpy arrays (nlev,) of bias by pressure level
    #              for each set
    #    countList: list of numpy arrays (nlev,) of ob-count by pressure
    #               level for each set
    #    nameList: list of names for each set (for figure legend)
    #    plevList: list of profile pressure levels (str or float)
    #    colMap: 2-element list containing colormap and colormap-range for
    #            panels (default: ['tab10',10]))
    # OUTPUTS
    #    plotFig: plot figure, as plt.fig()
    #
    # DEPENDENCIES
    # numpy
    # matplotlib.rc
    # matplotlib.pyplot
    # matplotlib.cm
    ######################################################################
    import numpy as np
    import matplotlib
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    #
    # Set font size and type
    #
    font = {'family': 'DejaVu Sans',
            'weight': 'bold',
            'size': 22}

    matplotlib.rc('font', **font)
    #
    # Identify number of profile sets (should be indentical for rmsList,
    # biasList, countList)
    #
    n_profiles = len(rmsList)
    #
    # Define colormap: The default settings are to select a range of 10 on
    # 'tab10', so that 10 pairs of rms/bias profiles can be produced for
    # the left panel that correspond to 10 ob-counts on the right panel.
    # The user can select a different colormap and range with colMapLeft
    # and colMapRight options.
    # If you want to sample the entire colorbar, set the range to the
    # number of datasets plotted.
    #
    # scalarMapList is used to select colors for each profile/bar
    scalarMap = cm.ScalarMappable(cmap=colMap[0])
    scalarMapList = scalarMap.to_rgba(range(colMap[1]))
    #
    # Generate figure
    plt.figure(figsize=(18, 18))
    # Define offset values (needed to define y-axis limits on both plots,
    # for consistency)
    y_offset = 8.0*(np.arange(n_profiles)-np.floor(n_profiles/2))
    # LEFT PANEL: rms and bias scores
    plt.subplot(121)
    # For each set, plot rms and bias: plot lines and circular markers
    legend_list = []
    for i in range(n_profiles):
        rms = rmsList[i]
        bias = biasList[i]
        levs = plevList[i]
        n_levs = np.size(levs)
        # Define y-axis limits
        y_min = min(levs)+min(y_offset)-8.0
        y_max = max(levs)+max(y_offset)+8.0
        # rms profile
        legend_list.append(nameList[i]+' rms')
        prof_color = list(scalarMapList[i][0:3])
        plt.plot(rms, levs, color=prof_color, linewidth=3)
        plt.plot(rms, levs, 'o', color=prof_color, markersize=8,
                 label='_nolegend_')
        # bias profile
        legend_list.append(nameList[i]+' bias')
        prof_color = list(scalarMapList[i][0:3])
        plt.plot(bias, levs, color=prof_color, linewidth=3,
                 linestyle='dashdot')
        plt.plot(bias, levs, 'o', color=prof_color, markersize=8,
                 label='_nolegend_')
    # Zero-line
    plt.plot(np.zeros((n_levs, )), levs, color='k', linewidth=1,
             linestyle='dashed', label='_nolegend_')
    # Set y-limits to entire range
    plt.ylim((y_min, y_max))
    plt.yticks(levs)
    # Reverse y-axis, if levs is in descending-order (often the case with
    # pressure coordinate data)
    if (levs[1] < levs[0]):
        plt.gca().invert_yaxis()
    # Set legend
    plt.legend(legend_list, frameon=False, fontsize=10)
    # Set x-label
    plt.xlabel('RMS or Bias')
    # RIGHT PANEL: ob-counts
    plt.subplot(122)
    # For each set, plot ob-counts and generate legend list
    # Counts are asserted to be in thousands
    legend_list = []
    for i in range(n_profiles):
        count = 0.001*countList[i]
        bar_color = list(scalarMapList[i][0:3])
        plt.barh(levs+y_offset[i], count, height=8.0, color=bar_color)
    # Set y-limits to entire range
    plt.ylim((y_min, y_max))
    plt.yticks(levs)
    # Reverse y-axis, if levs is in descending-order (often the case with
    # pressure coordinate data)
    if (levs[1] < levs[0]):
        plt.gca().invert_yaxis()
    # Set legend
    plt.legend(nameList, frameon=False, fontsize=10)
    # Set x-label
    plt.xlabel('Ob Count (Thousands)')
    # Turn off interactive-mode to suppress plotting figure
    plt.ioff()
    # Return
    return plt.gcf()


def plot_stat_traces(rmsList, biasList, countList, nameList, dateList,
                     tskip=4, colMap=['tab10', 10]):
    ######################################################################
    # Generates 2-panel plot:
    #    Top: Trace of rms and bias for each set
    #    Bottom: Bar-chart of ob-counts for each set
    #
    # INPUTS
    #    rmsList: list of numpy arrays (nlev,) of rms by pressure level
    #             for each set
    #    biasList: list of numpy arrays (nlev,) of bias by pressure level
    #              for each set
    #    countList: list of numpy arrays (nlev,) of ob-count by pressure
    #               level for each set
    #    nameList: list of names for each set (for figure legend)
    #    dateList: list of dates (str in '%Y%m%d%H' format)
    #    tskip: number of ticks to skip when writing tick-labels (default: 4)
    #    colMap: 2-element list containing colormap and colormap-range for
    #            panels (default: ['tab10',10]))
    # OUTPUTS
    #    plotFig: plot figure, as plt.fig()
    #
    # DEPENDENCIES
    # numpy
    # datetime
    # matplotlib.rc
    # matplotlib.pyplot
    # matplotlib.cm
    ######################################################################
    import numpy as np
    from datetime import datetime
    import matplotlib
    import matplotlib.pyplot as plt
    import matplotlib.cm as cm
    #
    # Set font size and type
    #
    font = {'family': 'DejaVu Sans',
            'weight': 'bold',
            'size': 22}

    matplotlib.rc('font', **font)
    #
    # Identify number of trace sets (should be indentical for rmsList,
    # biasList, countList)
    #
    n_trace = len(rmsList)
    #
    # Define colormap: The default settings are to select a range of 10 on
    # 'tab10', so that 10 pairs of rms/bias profiles can be produced for
    # the left panel that correspond to 10 ob-counts on the right panel.
    # The user can select a different colormap and range with colMapLeft
    # and colMapRight options. If you want to sample the entire colorbar,
    # set the range to the number of datasets plotted.
    #
    # scalarMapList is used to select colors for each profile/bar
    scalarMap = cm.ScalarMappable(cmap=colMap[0])
    scalarMapList = scalarMap.to_rgba(range(colMap[1]))
    #
    # Generate figure
    plt.figure(figsize=(18, 18))
    # Define offset values (needed to define x-axis limits on both plots,
    # for consistency)
    offs = 0.8/n_trace
    x_offset = offs*(np.arange(n_trace)-np.floor(n_trace/2))
    # TOP PANEL: rms and bias traces
    plt.subplot(211)
    # For each set, plot rms and bias: plot both lines and circular markers
    legend_list = []
    for i in range(n_trace):
        rms = rmsList[i]
        bias = biasList[i]
        dates = dateList[i]
        # Convert dates from %Y%m%d%H format to %b%d:%HZ tick-label format
        dstr = []
        for d in dates:
            dt = datetime.strptime(d, '%Y%m%d%H')
            dstr.append(datetime.strftime(dt, '%b%d:%HZ'))
        n_dates = np.size(dates)
        # Define x-axis limits
        x_min = min(x_offset)-offs
        x_max = n_dates+max(x_offset)+offs
        x_rng = np.arange(1, n_dates+1.0E-05)
        # rms trace
        legend_list.append(nameList[i]+' rms')
        prof_color = list(scalarMapList[i][0:3])
        plt.plot(x_rng, rms, color=prof_color, linewidth=3)
        plt.plot(x_rng, rms, 'o', color=prof_color, markersize=8,
                 label='_nolegend_')
        # bias profile
        legend_list.append(nameList[i]+' bias')
        prof_color = list(scalarMapList[i][0:3])
        plt.plot(x_rng, bias, color=prof_color, linewidth=3,
                 linestyle='dashdot')
        plt.plot(x_rng, bias, 'o', color=prof_color, markersize=8,
                 label='_nolegend_')
    # Zero-line
    plt.plot(x_rng, np.zeros((n_dates, )), color='k', linewidth=1,
             linestyle='dashed', label='_nolegend_')
    # Set x-limits to entire range
    plt.xlim((x_min, x_max))
    plt.xticks(ticks=x_rng[::tskip], labels=dstr[::tskip], fontsize=10)
    # Set legend
    plt.legend(legend_list, frameon=False, fontsize=10)
    # Set y-label
    plt.ylabel('RMS or Bias')
    # BOTTOM PANEL: ob-counts
    plt.subplot(212)
    # For this plot, bars need to be offset from each other so that they
    # all cluster around the pressure
    # level. This is accomplished with an offset value added to each bar
    legend_list = []
    for i in range(n_trace):
        count = 0.001*countList[i]
        bar_color = list(scalarMapList[i][0:3])
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


if __name__ == "__main__":
    ######################################################################
    import matplotlib.pyplot as plt
    #
    # Parse yaml file and define settings for plotting
    #
    (dt_beg, dt_end, hrs, pname, tname, tskip,
     setdict) = parse_yaml('proftrace_yaml.yaml')
    cycles = get_datelist(dt_beg, dt_end, hrs)
    profs_filename = pname
    trace_filename = tname
    #
    ######################################################################
    #
    # Collect rmse, bias, and ob-count statistics for each data (sub)set
    # in setdict
    #
    # If a gsistat file is missing for a given cycle, it raises a
    # FileNotFoundError and the cycle is skipped, not appearing at all in
    # the time-series
    #
    # If a gsistat file exists for a given cycle but no data for a data
    # (sub)set is present, the rmse, bias, and ob-count are set to zero
    # but the cycle is not skipped. Statistics with zero ob-count are
    # reassigned to NaN later.
    #
    rmses, biases, counts, levels = collect_statistics(setdict)
    ######################################################################
    #
    # Generate plotting data by aggregating rmse, bias, and count data
    # across times (for profiles), or aggregate full-column data (for
    # traces)
    #
    (rmse_profs, rmse_trace, bias_profs, bias_trace,
     nobs_profs, nobs_trace, levl_profs, date_trace,
     name_list) = aggregate_figure_data(rmses, biases, counts, levels)
    ######################################################################
    #
    # Generate profile plot
    #
    fig_prof = plot_stat_profiles(rmse_profs, bias_profs, nobs_profs,
                                  name_list, levl_profs)
    plt.ioff()
    fig_prof.savefig(profs_filename, bbox_inches='tight', facecolor='w')
    #
    # Generate trace plot
    #
    fig_trace = plot_stat_traces(rmse_trace, bias_trace, nobs_trace,
                                 name_list, date_trace, tskip=tskip)
    plt.ioff()
    fig_trace.savefig(trace_filename, bbox_inches='tight', facecolor='w')
    ######################################################################
