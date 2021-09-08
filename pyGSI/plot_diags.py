#!/usr/bin/env python
# Created by Kevin Dougherty
# October 2020
# Updated January 2021

from datetime import datetime
from netCDF4 import Dataset
import cartopy.feature as cfeature
import cartopy.crs as ccrs
import matplotlib.colors as mcolors
import matplotlib.pyplot as plt
import numpy as np
from textwrap import TextWrapper
import matplotlib
matplotlib.use('agg')

def _get_map_params(region):
        
    map_params = {'global': {'figdims': (15,12),
                             'cenlon': 0,
                             'extent': [-180, 180, -90, 90],
                             'textloc': [185, 55],
                             'titlesize': 14,
                             'textsize': 14,
                             'markersize': 30},
                  'CONUS':  {'figdims': (10,8),
                             'cenlon': 0,
                             'extent': [-130, -60, 10, 60],
                             'textloc': [-58, 53],
                             'titlesize': 12,
                             'textsize': 12,
                             'markersize': 10},
                  'North America': {'figdims': (10,8),
                             'cenlon': 0,
                             'extent': [-170, -50, 7.5, 75],
                             'textloc': [-48, 65],
                             'titlesize': 12,
                             'textsize': 12,
                             'markersize': 15},
                  'Alaska': {'figdims': (10,8),
                             'cenlon': 0,
                             'extent': [-170, -140, 50, 75],
                             'textloc': [-139, 70],
                             'titlesize': 12,
                             'textsize': 12,
                             'markersize': 10},
                  'Europe': {'figdims': (10,8),
                             'cenlon': 0,
                             'extent': [-12.5, 40, 30, 70],
                             'textloc': [41, 64],
                             'titlesize': 12,
                             'textsize': 12,
                             'markersize': 15},
                  'Africa': {'figdims': (10,8),
                             'cenlon': 0,
                             'extent': [-20, 55, -35, 40],
                             'textloc': [57, 30],
                             'titlesize': 12,
                             'textsize': 12,
                             'markersize': 15},
                  'Asia':   {'figdims': (10,8),
                             'cenlon': 0,
                             'extent': [30, 160, -5, 80],
                             'textloc': [163, 66],
                             'titlesize': 12,
                             'textsize': 12,
                             'markersize': 20},
                  'Australia': {'figdims': (10,8),
                             'cenlon': 0,
                             'extent': [110, 155, -10, -42.5],
                             'textloc': [156, -17],
                             'titlesize': 12,
                             'textsize': 12,
                             'markersize': 15},
                  'South America': {'figdims': (10,8),
                             'cenlon': 0,
                             'extent': [-100, -30, 15, -55],
                             'textloc': [-28, 3],
                             'titlesize': 12,
                             'textsize': 12,
                             'markersize': 15},
                  'Baltic': {'figdims': (10,8),
                             'cenlon': 0,
                             'extent': [5, 35, 50, 70],
                             'textloc': [35.5, 66],
                             'titlesize': 12,
                             'textsize': 12,
                             'markersize': 10},
                  'Japan':  {'figdims': (10,8),
                             'cenlon': 0,
                             'extent': [125, 155, 25, 55],
                             'textloc': [156, 50],
                             'titlesize': 12,
                             'textsize': 12,
                             'markersize': 10}
                 }    
    
    try:
        mapdict = map_params[region]
    except KeyError:
        mapdict = map_params['global']
        print('Region entered is not included in map params. Returning global view ...')
    
    return mapdict

def _get_obs_type(obs_id):
    """
    Grabs the full name of a specific observation.
    Input:
        obs_id: ID number(s) for conventional diagnostic
    Output:
        list of observation names. If ID in the list, it will return,
        the proper name. If no observation ID, returns 'All Observations'.
        If ID not in the list, returns list of string of the ID number.
    """

    obs_indicators = {
        120: "Rawinsonde",
        126: "RASS",
        130: "Aircraft: AIREP and PIREP",
        131: "Aircraft: AMDAR",
        132: "Flight-Level Reconnaissance and Profile Dropsonde",
        133: "Aircraft: MDCRS ACARS",
        134: "Aircraft: TAMDAR",
        135: "Aircraft: Canadian AMDAR",
        153: "GPS-Integrated Precipitable Water",
        180: "Surface Marine w/ Station Pressure (Ship, Buoy, C-MAN, Tide Guage)",
        181: "Surface Land w/ Station Pressure (Synoptic, METAR)",
        182: "Splash-Level Dropsonde Over Ocean",
        183: "Surface Marine or Land - Missing Station Pressure",
        187: "Surface Land - Missing Station Pressure",
        210: "Synthetic Tropical Cyclone",
        220: "Rawinsonde",
        221: "PIBAL",
        224: "NEXRAD Vertical Azimuth Display",
        228: "Wind Profiler: JMA",
        229: "Wind Profiler: PIBAL",
        230: "Aircraft: AIREP and PIREP",
        231: "Aircraft: AMDAR",
        232: "Flight-Level Reconnaissance and Profile Dropsonde",
        233: "Aircraft: MDCRS ACARS",
        234: "Aircraft: TAMDAR",
        235: "Aircraft: Canadian AMDAR",
        242: "JMA IR (Longwave) and Visible Cloud Drift Below 850mb (GMS, MTSAT, HIMAWARI)",
        243: "EUMETSAT IR (Longwave) and Visible Cloud Drift Below 850mb (METEOSAT)",
        244: "AVHRR/POES IR (Longwave) Cloud Drift",
        245: "NESDIS IR (Longwave) Cloud Drift (All Levels)",
        246: "NESDIS Imager Water Vapor (All Levels) - Cloud Top (GOES)",
        250: "JMA Imager Water Vapor (All Levels) - Cloud Top & Deep Layer (GMS, MTSAT, HIMAWARI)",
        251: "NESDIS Visible Cloud Drift (All Levels) (GOES)",
        252: "JMA IR (Longwave) and Visible Cloud Drift Above 850mb (GMS, MTSAT, HIMAWARI)",
        253: "EUMETSAT IR (Longwave) and Visible Cloud Drift Above 850mb (METEOSAT)",
        254: "EUMETSAT Imager Water Vapor (All Levels) - Cloud Top & Deep Layer (METEOSAT)",
        257: "MODIS/POES IR (Longwave) Cloud Drift (All Levels) (AQUA, TERRA)",
        258: "MODIS/POES Imager Water Vapor (All Levels) - Cloud Top (AQUA, TERRA)",
        259: "MODIS/POES Imager Water Vapor (All Levels) - Deep Layer (AQUA, TERRA)",
        280: "Surface Marine w/ Station Pressure (Ship, Buoy, C-MAN, Tide Guage)",
        281: "Surface Land w/ Station Pressure (Synoptic, METAR)",
        282: "ATLAS Buoy",
        284: "Surface Marine or Land - Missing Station Pressure",
        287: "Surface Land (METAR) - Missing Station Pressure",
        289: "SUPEROBED (1.0 Lat/Lon) Scatterometer Winds over Ocean",
        290: "Non-SUPEROBED Scatterometer Winds over Ocean"
    }

    descripts = list()
    for ids in obs_id:
        if (ids in obs_indicators.keys()) == True:
            descripts.append(obs_indicators[ids])

    if descripts:
        return descripts
    elif not obs_id:
        return ['All Observations']
    else:
        return [str(x) for x in obs_id]

def _conv_features_dict(stats):

    conv_features_dict = {'t': {'cmap': 'jet',
                               'extend': 'both',
                               'upperbound': np.round(stats['Mean']) + (np.round(stats['Std']*2)/2)*2,
                               'lowerbound': np.round(stats['Mean']) - (np.round(stats['Std']*2)/2)*2,
                               'bins': ((np.round(stats['Mean']) + (np.round(stats['Std']*2)/2)*2) -
                                        (np.round(stats['Mean']) - (np.round(stats['Std']*2)/2)*2))/20},
                         'q': {'cmap': 'GnBu',
                               'extend': 'max',
                               'upperbound': stats['Std']*5,
                               'lowerbound': 0,
                               'bins': (stats['Std']*5)/10},
                         'ps': {'cmap': 'jet',
                                'extend': 'both',
                                'upperbound': 1050,
                                'lowerbound': 750,
                                'bins': (1050 - 750)/10},
                         'sst': {'cmap': 'Spectral_r',
                                 'extend': 'both',
                                 'upperbound': np.round(stats['Mean']) + (np.round(stats['Std']*2)/2)*2,
                                 'lowerbound': np.round(stats['Mean']) - (np.round(stats['Std']*2)/2)*2,
                                 'bins': ((np.round(stats['Mean']) + (np.round(stats['Std']*2)/2)*2) -
                                          (np.round(stats['Mean']) - (np.round(stats['Std']*2)/2)*2))/20},
                         'pw': {'cmap': 'GnBu',
                                'extend': 'max',
                                'upperbound': stats['Std']*5,
                                'lowerbound': 0,
                                'bins': (stats['Std']*5)/10},
                         'tcp': {'cmap': 'Reds_r',
                                 'extend': 'both',
                                 'upperbound': 1050,
                                 'lowerbound': 950,
                                 'bins': (1050 - 950)/10},
                         'u': {'cmap': 'viridis',
                               'extend': 'both',
                               'upperbound': np.round(stats['Std'])*5,
                               'lowerbound': np.round(stats['Std'])*-5,
                               'bins': ((np.round(stats['Std'])*5) - (np.round(stats['Std'])*-5))/10},
                         'v': {'cmap': 'viridis',
                               'extend': 'both',
                               'upperbound': np.round(stats['Std'])*5,
                               'lowerbound': np.round(stats['Std'])*-5,
                               'bins': ((np.round(stats['Std'])*5) - (np.round(stats['Std'])*-5))/10},
                         'windspeed': {'cmap': 'Spectral_r',
                                       'extend': 'max',
                                       'upperbound': np.round(stats['Std'])*5,
                                       'lowerbound': 0,
                                       'bins': (stats['Std']*5)/10}
                         }

    return conv_features_dict

def _myround(x, base):
    """
    x    : number needing to be rounded
    base : the interval you would to round to
    """
    if isinstance(base, float):
        return base * round(float(x)/base)
    else:
        return int(base * round(float(x)/base))

def _get_xlabel(metadata):
    """
    Creates and returns x label for a plot.
    """

    conv_xlabels = {'t': "Temperature (K)",
                    'q': "Specific Humidity (kg/kg)",
                    'sst': "Sea Surface Temperature (K)",
                    'pw': "Precipitable Water (mm)",
                    'ps': "Pressure (hPa)",
                    'tcp': "Pressure (hPa)",
                    'u': "Windspeed (m/s)",
                    'v': "Windspeed (m/s)",
                    'windspeed': "Windspeed (m/s)"
                   }

    # Observation minus Forecast
    if metadata['Diag Type'] in ['omf', 'o-f']:
        xlabel = 'Observation - Forecast'

    # Observation minus Analysis
    elif metadata['Diag Type'] in ['oma', 'o-a']:
        xlabel = 'Observation - Analysis'

    # Conventional Data
    elif metadata['Diag File Type'] == 'conventional':
        xlabel = conv_xlabels[metadata['Variable']]

    # Land, Water, Snow, Ice, or Cloud Fraction Data
    elif metadata['Diag Type'].split('_')[-1] == 'fraction':
        xlabel = '%s %s' % (metadata['Diag Type'].split(
            '_')[0], metadata['Diag Type'].split('_')[-1])

    # Ozone Data
    elif metadata['Diag File Type'] == 'ozone':
        xlabel = "Dobson Units"

    # Radiance Data
    elif metadata['Diag File Type'] == 'radiance':
        xlabel = "Brightness Temperature (K)"

    else:
        xlabel = None


    return xlabel

def _get_stats_labels(metadata,stats):

    if metadata['Diag File Type'] == 'conventional' and metadata['Variable'] == 'q':
        roundnum = 6
    else:
        roundnum = 3

    t = ('n: %s\nstd: %s\nmean: %s\nmax: %s\nmin: %s' % (stats['N'],
                                                         np.round(
                                                             stats['Std'], roundnum),
                                                         np.round(
                                                             stats['Mean'], roundnum),
                                                         np.round(
                                                             stats['Max'], roundnum),
                                                         np.round(stats['Min'], roundnum)))

    return t

def _get_title(metadata):

    # Handles analysis use data
    anl_use = True if 'Anl Use' in metadata and metadata['Anl Use'] else False

    # Left Title label
    if anl_use:
        if metadata['Diag File Type'] == 'conventional':

            conv_title_dict = {'assimilated' : {'title': '{Obs Type}: {Variable} - {Diag Type} - Data Assimilated\n'.format(**metadata) + \
                                                '%s' % '\n'.join(metadata['ObsID Name'])},
                               'rejected'    : {'title': '{Obs Type}: {Variable} - {Diag Type} - Data Rejected\n'.format(**metadata) + \
                                                '%s' % '\n'.join(metadata['ObsID Name'])},
                               'monitored'   : {'title': '{Obs Type}: {Variable} - {Diag Type} - Data Monitored\n'.format(**metadata) + \
                                                '%s' % '\n'.join(metadata['ObsID Name'])}
                              }

            title = conv_title_dict[metadata['Anl Use Type']]['title']

        elif metadata['Diag File Type'] == 'ozone':

            ozone_title_dict = {'assimilated': {'title': '{Obs Type}: {Satellite} - {Diag Type} - Data Assimilated\nLayer: {Layer}'.format(
                                                **metadata)},
                                'rejected'   : {'title': '{Obs Type}: {Satellite} - {Diag Type} - Data Rejected\nLayer: {Layer}'.format(
                                                **metadata)},
                                'monitored'  : {'title': '{Obs Type}: {Satellite} - {Diag Type} - Data Monitored\nLayer: {Layer}'.format(
                                                **metadata)}
                                }

            title = ozone_title_dict[metadata['Anl Use Type']]['title']

    else:
        if metadata['Diag File Type'] == 'conventional':
            title = '{Obs Type}: {Variable} - {Diag Type}\n'.format(**metadata) + \
            '%s' % '\n'.join(metadata['ObsID Name'])

        elif metadata['Diag File Type'] == 'ozone':
            title = '{Obs Type}: {Satellite} - {Diag Type} \nLayer: {Layer}'.format(
            **metadata)

        else:
            title = '{Obs Type}: {Satellite} - {Diag Type}\n'.format(
            **metadata) + 'Channels: %s' % ' '.join(str(x) for x in metadata['Channels'])

    return title

def _get_save_file(metadata):

    # Save file label
    if metadata['Diag File Type'] == 'conventional':
        save_file = '{Date:%Y%m%d%H}_{Obs Type}_{Variable}_{Diag Type}_'.format(
            **metadata) + '%s' % '_'.join(str(x) for x in metadata['ObsID'])

    elif metadata['Diag File Type'] == 'ozone':
        layer = '_'.join(metadata['Layer'].split()) if metadata['Layer'] == 'column total' else f"{metadata['Layer']:.3f}"

        save_file = '{Date:%Y%m%d%H}_{Obs Type}_{Satellite}_{Diag Type}_'.format(
            **metadata) + '%s' % layer

    else:
        save_file = '{Date:%Y%m%d%H}_{Obs Type}_{Satellite}_{Diag Type}_'.format(
            **metadata) + 'channels_%s' % '_'.join(str(x) for x in metadata['Channels'])

    # Handles analysis use data
    anl_use = True if 'Anl Use' in metadata and metadata['Anl Use'] else False

    if anl_use:
        if metadata['Anl Use Type'] == 'assimilated':
            save_file = save_file + '_assimilated'
        elif metadata['Anl Use Type'] == 'rejected':
            save_file = save_file + '_rejected'
        else:
            save_file = save_file + '_monitored'

    return save_file


def _plot_labels(metadata, stats):

    # Stats label
    if not stats:
        t = None
    else:
        t = _get_stats_labels(metadata, stats)

    # Date label
    date_title = metadata['Date'].strftime("%d %b %Y %Hz")

    # X Label
    xlabel = _get_xlabel(metadata)

    # Left Title label
    left_title = _get_title(metadata)

    # Save file label
    save_file = _get_save_file(metadata)

    labels = {'stat text': t,
              'x label': xlabel,
              'left title': left_title,
              'date title': date_title,
              'save file': save_file
              }

    return labels

def _colorbar_features(metadata, stats):
    """
    Returns colormaps, the boundary norm (generates a colormap
    index based on discrete intervals), and how to properly
    extend the colorbar based on the type of data being used.
    """

    # Get cmap, bins, norm and extend for O-F and O-A
    if metadata['Diag Type'] in ['omf', 'o-f', 'omb', 'o-b', 'oma', 'o-a']:
        cmap = 'bwr'

        upperbound = (np.round(stats['Std']*2)/2)*5
        if upperbound == 0:
            upperbound = stats['Std']*5

        lowerbound = 0-upperbound
        bins = (upperbound - lowerbound)/10

        norm = mcolors.BoundaryNorm(boundaries=np.arange(
            lowerbound, upperbound+bins, bins), ncolors=256)

        extend = 'both'

    # Get cmap, bins, norm and extend for observations, hofx, and windspeed
    elif metadata['Diag Type'] in ['observation', 'hofx', 'windspeed']:
        if metadata['Diag File Type'] == 'conventional':

            conv_features_dict = _conv_features_dict(stats)

            cmap = conv_features_dict[metadata['Variable']]['cmap']
            extend = conv_features_dict[metadata['Variable']]['extend']
            upperbound = conv_features_dict[metadata['Variable']]['upperbound']
            lowerbound = conv_features_dict[metadata['Variable']]['lowerbound']
            bins = conv_features_dict[metadata['Variable']]['bins']
        else:
            cmap = 'viridis'
            extend = 'both'
            upperbound = np.round(stats['Mean']) + \
                (np.round(stats['Std']*2)/2)*2
            lowerbound = np.round(stats['Mean']) - \
                (np.round(stats['Std']*2)/2)*2
            bins = (upperbound - lowerbound)/20

        norm = mcolors.BoundaryNorm(boundaries=np.arange(
            lowerbound, upperbound+bins, bins), ncolors=256)

    else:
        cmap = 'bwr'

        upperbound = 1
        lowerbound = 0
        bins = 0.5

        norm = mcolors.BoundaryNorm(boundaries=np.arange(
            lowerbound, upperbound+bins, bins), ncolors=256)

        extend = 'neither'

    return cmap, norm, extend


def _calculate_stats(data):
    """
    Calculates n, mean, min, max,
    standard deviation, and RMSE.
    Returns dictionary with stats.
    """

    n = np.count_nonzero(~np.isnan(data))

    mean = np.nanmean(data)
    std = np.nanstd(data)
    mx = np.nanmax(data)
    mn = np.nanmin(data)

    rmse = np.sqrt(np.nanmean(np.square(data)))

    stats = {'N': n,
             'Min': mn,
             'Max': mx,
             'Mean': mean,
             'Std': std,
             'RMSE': rmse
             }

    return stats

def _create_histogram_plot(data, metadata, outdir='./'):

    # Create Figure
    fig = plt.figure(figsize=(8, 6))
    ax = fig.add_subplot(111)

    # Plots 'No Data' if no obs
    if len(data) == 0:
        stats = False
        labels = _plot_labels(metadata, stats)
        ax.text(0.5, 0.5, 'No Data', fontsize=32,
                alpha=0.6, ha='center')

    # Plots 'Single Observation' if one ob
    elif len(data) == 1:
        stats = _calculate_stats(data)
        labels = _plot_labels(metadata, stats)

        ax.text(0.75, .7, labels['stat text'],
                fontsize=14, transform=ax.transAxes)
        ax.text(0.5, 0.5, 'Single Observation',
                fontsize=32, alpha=0.6, ha='center')

    else:
        # Calculate stats
        stats = _calculate_stats(data)

        # Get binsize
        binsize = (stats['Max']-stats['Min'])/np.sqrt(stats['N'])

        if metadata['Diag File Type'] == 'conventional' and metadata['Variable'] == 'windspeed':
            bins = np.arange(
                0, stats['Mean']+(4*stats['Std']), binsize)
        elif metadata['Diag Type'] in ['omf', 'o-f', 'omb', 'o-b', 'oma', 'o-a']:
            bins = np.arange(
                0-(4*stats['Std']), 0+(4*stats['Std']), binsize)
        else:
            bins = np.arange(
                stats['Mean']-(4*stats['Std']), stats['Mean']+(4*stats['Std']), binsize)

        # Plots histogram data
        plt.hist(data, bins=bins)

        # Plots data mean with red line
        plt.axvline(stats['Mean'], color='r',
                    linestyle='solid', linewidth=1)

        # If omf or oma, plots a line a black line at 0
        if metadata['Diag Type'] in ['omf', 'o-f', 'omb', 'o-b', 'oma', 'o-a']:
            plt.axvline(0, color='k', linestyle='dashed', linewidth=1)

        # Get labels
        labels = _plot_labels(metadata, stats)

        # Plots text of stats
        ax.text(0.75, .7, labels['stat text'],
                fontsize=14, transform=ax.transAxes)

    # Plots x and y label
    plt.xlabel(labels['x label'])
    plt.ylabel('Count')

    # Makes sure plot title does not run off if long observation name
    wrapper = TextWrapper(
        width=70, break_long_words=False, replace_whitespace=False)
    left_title = '\n'.join(wrapper.wrap(labels['left title']))

    # Plots title and saves fig
    plt.title(left_title, loc='left', fontsize=14)
    plt.title(labels['date title'], loc='right',
              fontweight='semibold', fontsize=14)
    plt.savefig(outdir+'/%s_histogram.png' %
                labels['save file'], bbox_inches='tight', pad_inches=0.1)
    plt.close('all')

    return

def _no_data_spatial(metadata, outdir='./', region='global'):
    """
    Plots spatial map with 'No Data' printed across
    middle of plot when there is no data returned
    """

    stats = False
    labels = _plot_labels(metadata, stats)
    figdict = _get_map_params(region)

    fig = plt.figure(figsize=figdict['figdims'])
    ax = fig.add_subplot(
        1, 1, 1,
        projection=ccrs.PlateCarree(central_longitude=figdict['cenlon']))

    ax.add_feature(cfeature.GSHHSFeature(scale='auto'))
    ax.set_extent(figdict['extent'])
    ax.text(0, 0, 'No Data', fontsize=32, alpha=0.6, ha='center')
    plt.title(labels['left title'], loc='left', fontsize=14)
    plt.title(labels['date title'], loc='right',
              fontweight='semibold', fontsize=14)
    plt.savefig(outdir+'/%s_spatial.png' %
                labels['save file'], bbox_inches='tight', pad_inches=0.1)
    plt.close('all')


def _create_spatial_plot(data, metadata, lats, lons,
                         outdir='./', region='global'):

    if len(data) == 0:
        _no_data_spatial(metadata, outdir, region=region)
        return
    else:
        stats = _calculate_stats(data)

        if len(data) == 1:
            cmap = 'Blues'
            norm = None
            extend = 'neither'
        else:
            cmap, norm, extend = _colorbar_features(metadata, stats)

    figdict = _get_map_params(region)
    # Creates Figure
    fig = plt.figure(figsize=figdict['figdims'])
    ax = fig.add_subplot(
        1, 1, 1,
        projection=ccrs.PlateCarree(central_longitude=figdict['cenlon']))

    # Adds plot features (country borders) and set extent of map
    ax.add_feature(cfeature.GSHHSFeature(scale='auto'))
    ax.set_extent(figdict['extent'])

    #plot the data
    cs = plt.scatter(lons, lats, c=data, s=figdict['markersize'],
                     norm=norm, cmap=cmap,
                     transform=ccrs.PlateCarree())

    # Get and plot labels
    labels = _plot_labels(metadata, stats)

    ax.text(figdict['textloc'][0], figdict['textloc'][1], labels['stat text'], fontsize=figdict['textsize'])

    plt.title(labels['left title'], loc='left', fontsize=figdict['titlesize'])
    plt.title(labels['date title'], loc='right',
              fontweight='semibold', fontsize=figdict['titlesize'])

    # Plot colorbar
    cb = plt.colorbar(cs, orientation='horizontal',
                      shrink=0.5, pad=.04, extend=extend)
    cb.set_label(labels['x label'], fontsize=12)

    # Save figure
    plt.savefig(outdir+'/%s_spatial.png' %
                labels['save file'], bbox_inches='tight', pad_inches=0.1)
    plt.close('all')


def _binned_plot_features(binned_var, metadata, stats):
    """
    Returns colormaps, the boundary norm (generates a colormap
    index based on discrete intervals), and how to properly
    extend the colorbar based on the type of data being used.
    """

    features_dict = {'binned_nobs': {'cmap': 'plasma',
                                     'extend': 'max',
                                     'upperbound': np.round(stats['Max']),
                                     'lowerbound': 0,
                                     'bins': _myround(stats['Max'], 25)/10,
                                     'x label': '# of Observations',
                                     'title': 'Binned Number of Observations'},
                     'binned_mean': {'x label': 'Binned Average',
                                     'title': 'Binned Mean'},
                     'binned_max':  {'cmap': 'Reds',
                                     'extend': 'max',
                                     'upperbound': np.round(stats['Max']),
                                     'lowerbound': np.round(stats['Min']),
                                     'bins': _myround(stats['Std'], 2),
                                     'x label': 'Binned Max',
                                     'title': 'Binned Max'
                                    },
                     'binned_min':  {'cmap': 'Blues_r',
                                     'extend': 'min',
                                     'upperbound': np.round(stats['Max']),
                                     'lowerbound': np.round(stats['Min']),
                                     'bins': _myround(stats['Std'], 2),
                                     'x label': 'Binned Min',
                                     'title': 'Binned Min'
                                    },
                     'binned_std':  {'cmap': 'plasma',
                                     'extend': 'max',
                                     'upperbound': np.round(stats['Max']),
                                     'lowerbound': np.round(stats['Min']),
                                     'bins':  _myround(stats['Std'], 0.5),
                                     'x label': 'Binned Std. Dev.',
                                     'title': 'Binned Standard Deviation'
                                    },
                     'binned_rmse': {'cmap': 'plasma',
                                     'extend': 'max',
                                     'upperbound': np.round(stats['Max']),
                                     'lowerbound': np.round(stats['Min']),
                                     'bins': _myround(stats['Std'], 0.5),
                                     'x label': 'Binned RMSE',
                                     'title': 'Binned Root Mean Square Error'
                                    }
                    }

    # Get cmap, bins, norm and extend for O-F and O-A
    if metadata['Diag Type'] in ['omf', 'o-f', 'omb', 'o-b', 'oma', 'o-a'] and binned_var =='binned_mean':
        cmap = 'bwr'

        upperbound = (np.round(stats['Std']*2)/2)*5
        if upperbound == 0:
            upperbound = stats['Std']*5

        lowerbound = 0-upperbound
        bins = (upperbound - lowerbound)/10

        norm = mcolors.BoundaryNorm(boundaries=np.arange(
            lowerbound, upperbound+bins, bins), ncolors=256)

        extend = 'both'


    else:
        cmap = features_dict[binned_var]['cmap']
        extend = features_dict[binned_var]['extend']

        upperbound = features_dict[binned_var]['upperbound']
        lowerbound = features_dict[binned_var]['lowerbound']
        bins = features_dict[binned_var]['bins']

        norm = mcolors.BoundaryNorm(boundaries=np.arange(
                lowerbound, upperbound+bins, bins), ncolors=256)


    # Stats label
    if not stats:
        t = None
    else:
        t = _get_stats_labels(metadata, stats)

    # Date label
    date_title = metadata['Date'].strftime("%d %b %Y %Hz")

    # Left Title label
    left_title = _get_title(metadata)
    left_title = '%s\n' % features_dict[binned_var]['title'] + left_title

    # X label
    if binned_var == 'binned_nobs':
        xlabel = features_dict[binned_var]['x label']
    else:
        xlabel = _get_xlabel(metadata)
        xlabel = xlabel + ' (%s)' % features_dict[binned_var]['x label']

    # Save file
    save_file = _get_save_file(metadata)

    labels = {'stat text': t,
              'x label': xlabel,
              'left title': left_title,
              'date title': date_title,
              'save file': save_file
              }


    return labels, cmap, norm, extend


def plot_histogram(data, metadata, outdir='./'):
    if metadata['Diag File Type'] == 'conventional':
        metadata['ObsID Name'] = _get_obs_type(metadata['ObsID'])

    # Handles analysis use data
    anl_use = True if 'Anl Use' in metadata and metadata['Anl Use'] else False

    # Handles uv data
    if metadata['Diag File Type'] == 'conventional' and metadata['Variable'] == 'uv':

        if anl_use:
            for anl_type in data.keys():
                for variable in data[anl_type].keys():

                    # Add variables to metadata
                    metadata['Variable'] = variable
                    metadata['Anl Use Type'] = anl_type

                    _create_histogram_plot(data[anl_type][variable], metadata, outdir=outdir)
            
            #metadata was being saved as windspeed, need to revert back to 'uv' for other plots
            metadata['Variable'] = 'uv'

        else:
            for variable in data.keys():

                metadata['Variable'] = variable

                _create_histogram_plot(data[variable], metadata, outdir=outdir)
            
            #metadata was being saved as windspeed, need to revert back to 'uv' for other plots
            metadata['Variable'] = 'uv'


    else:

        if anl_use:
            for anl_type in data.keys():
                metadata['Anl Use Type'] = anl_type

                _create_histogram_plot(data[anl_type], metadata, outdir=outdir)

        else:
            _create_histogram_plot(data, metadata, outdir=outdir)


    return

def plot_spatial(data, metadata, lats, lons, outdir='./', region='global'):
    
    if metadata['Diag File Type'] == 'conventional':
        metadata['ObsID Name'] = _get_obs_type(metadata['ObsID'])

    # Handles analysis use data
    anl_use = True if 'Anl Use' in metadata and metadata['Anl Use'] else False
    
    # Handles uv data
    if metadata['Diag File Type'] == 'conventional' and metadata['Variable'] == 'uv':
        if anl_use:
            for anl_type in data.keys():
                for variable in data[anl_type].keys():

                    # Add variables to metadata
                    metadata['Variable'] = variable
                    metadata['Anl Use Type'] = anl_type

                    _create_spatial_plot(data[anl_type][variable], metadata,
                                         lats[anl_type], lons[anl_type],
                                         outdir=outdir, region=region)
            #metadata was being saved as windspeed, need to revert back to 'uv' for other plots
            metadata['Variable'] = 'uv'

        else:
            for variable in data.keys():
                metadata['Variable'] = variable
                _create_spatial_plot(data[variable], metadata,
                                     lats, lons,
                                     outdir=outdir, region=region)
            #metadata was being saved as windspeed, need to revert back to 'uv' for other plots
            metadata['Variable'] = 'uv'

    else:

        if anl_use:
            for anl_type in data.keys():

                metadata['Anl Use Type'] = anl_type

                _create_spatial_plot(data[anl_type], metadata,
                                     lats[anl_type], lons[anl_type],
                                     outdir=outdir, region=region)

        else:
            _create_spatial_plot(data, metadata, lats, lons,
                                 outdir=outdir, region=region)

    return


def plot_binned_spatial(data, metadata, binned_var=None, binsize='1x1', outdir='./'):

    if binned_var is None:
        print('Please select a binned variable type i.e. binned_nobs, binned_mean,\n',
              'binned_max, binned_min, binned_std, binned_rmse.')
        return

    if metadata['Diag File Type'] == 'conventional':
        metadata['ObsID Name'] = _get_obs_type(metadata['ObsID'])

    # Create lats and lons based on binsize
    lonlen = 360
    latlen = 180

    lon_lowerlim = 0
    lon_upperlim = 360

    lat_lowerlim = -90
    lat_upperlim = 90

    if binsize.split('x')[0] != binsize.split('x')[1]:
        print('ERROR: Binsize must be square i.e. 1x1, 2x2, 5x5 etc. Please use different binsize.')

    binsize = int(binsize.split('x')[0])

    if latlen % binsize == 0 and lonlen % binsize == 0:
        latbin = int(latlen/binsize)
        lonbin = int(lonlen/binsize)
        n_deg = binsize/2

        ll_lats = np.linspace(lat_lowerlim+(n_deg),
                              lat_upperlim-(n_deg),
                              latbin)

        ll_lons = np.linspace(lon_lowerlim+(n_deg),
                              lon_upperlim-(n_deg),
                              lonbin)

    xx, yy = np.meshgrid(ll_lons, ll_lats)

    stats = _calculate_stats(data[binned_var])

    labels, cmap, norm, extend = _binned_plot_features(binned_var, metadata, stats)

    fig = plt.figure(figsize=(15, 12))
    ax = fig.add_subplot(
        1, 1, 1, projection=ccrs.PlateCarree(central_longitude=0))

    ax.add_feature(cfeature.GSHHSFeature(scale='auto'))
    ax.set_extent([-180, 180, -90, 90])

    cs = plt.pcolormesh(xx, yy, data[binned_var], #cmap='bwr',
                        norm=norm, cmap=cmap,
                        transform=ccrs.PlateCarree())

    ax.text(185, 55, labels['stat text'], fontsize=14)

    cb = plt.colorbar(cs, orientation='horizontal',
                      shrink=0.5, pad=.04, extend=extend)
    cb.ax.tick_params(labelsize=12)
    cb.set_label(labels['x label'], fontsize=13)

    plt.title(labels['left title'], loc='left', fontsize=14)
    plt.title(labels['date title'], loc='right',
              fontweight='semibold', fontsize=14)

    plt.savefig(outdir+'/%s_binned_spatial.png' %
                labels['save file'], bbox_inches='tight', pad_inches=0.1)
    plt.close('all')
