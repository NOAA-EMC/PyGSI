from datetime import datetime
import numpy as np
from textwrap import TextWrapper
import matplotlib.pyplot as plt
from emcpy.plots.plots import Scatter, Histogram, VerticalLine
from emcpy.plots.map_plots import MapScatter
from emcpy.plots import CreateMap, CreatePlot, VariableSpecs
from emcpy.plots.map_tools import Domain, MapProjection
import matplotlib
matplotlib.use('agg')

__all__ = ['get_obs_type', 'get_bins', 'calculate_stats',
           'get_labels', 'varspecs_name']


def qc_flag_colors(qcflag):
    """
    Dictionary of colors to use for specific QC flag numbers.

    Args:
        qcflag : (int) qc flag number
    Returns:
        color : (str) color based on qc flag number
    """
    qc_flag_colors = {
        0: 'tab:green',
        1: 'black',
        2: 'tab:pink',
        3: 'tab:red',
        4: 'tab:orange',
        5: 'tab:brown',
        6: 'tab:olive',
        7: 'skyblue',
        8: 'tab:blue',
        9: 'maroon',
        10: 'tab:cyan',
        12: 'navy',
        50: 'indigo',
        51: 'tab:purple',
        52: 'magenta',
        53: 'cornflowerblue',
        54: 'pink',
        55: 'lightgrey',
        56: 'darkgrey',
        57: 'dimgrey'
    }

    return qc_flag_colors[qcflag]


def get_obs_type(obs_id):
    """
    Grabs the full name of a specific observation.

    Args:
        obs_id: (list of ints) ID number(s) for conventional
                diagnostic
    Returns:
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
        180: ("Surface Marine w/ Station Pressure (Ship, Buoy, C-MAN, "
              "Tide Guage)"),
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
        242: ("JMA IR (Longwave) and Visible Cloud Drift Below 850mb "
              "(GMS, MTSAT, HIMAWARI)"),
        243: ("EUMETSAT IR (Longwave) and Visible Cloud Drift Below "
              " 850mb (METEOSAT)"),
        244: "AVHRR/POES IR (Longwave) Cloud Drift",
        245: "NESDIS IR (Longwave) Cloud Drift (All Levels)",
        246: "NESDIS Imager Water Vapor (All Levels) - Cloud Top (GOES)",
        250: ("JMA Imager Water Vapor (All Levels) - Cloud Top & Deep "
              "Layer (GMS, MTSAT, HIMAWARI)"),
        251: "NESDIS Visible Cloud Drift (All Levels) (GOES)",
        252: ("JMA IR (Longwave) and Visible Cloud Drift Above 850mb "
              "(GMS, MTSAT, HIMAWARI)"),
        253: ("EUMETSAT IR (Longwave) and Visible Cloud Drift Above "
              "850mb (METEOSAT)"),
        254: ("EUMETSAT Imager Water Vapor (All Levels) - Cloud Top "
              " & Deep Layer (METEOSAT)"),
        257: ("MODIS/POES IR (Longwave) Cloud Drift (All Levels) "
              "(AQUA, TERRA)"),
        258: ("MODIS/POES Imager Water Vapor (All Levels) - Cloud Top "
              "(AQUA, TERRA)"),
        259: ("MODIS/POES Imager Water Vapor (All Levels) - Deep Layer "
              "(AQUA, TERRA)"),
        280: ("Surface Marine w/ Station Pressure (Ship, Buoy, C-MAN, "
              "Tide Guage)"),
        281: "Surface Land w/ Station Pressure (Synoptic, METAR)",
        282: "ATLAS Buoy",
        284: "Surface Marine or Land - Missing Station Pressure",
        287: "Surface Land (METAR) - Missing Station Pressure",
        289: "SUPEROBED (1.0 Lat/Lon) Scatterometer Winds over Ocean",
        290: "Non-SUPEROBED Scatterometer Winds over Ocean"
    }

    descripts = list()
    if obs_id is None:
        return ['All Observations']

    else:
        for ids in obs_id:
            if ids in obs_indicators.keys():
                descripts.append(obs_indicators[ids])

            else:
                descripts.append(str(ids))

        return descripts


def get_bins(eval_type, stats):
    """
    Calculates bins to use for histogram plot.

    Args:
        eval_type : (str) determines how to create bins
                    whether 'diff' or 'magnitude'
        stats : (dict) dictionary of stats from
                `calculate_stats()`
    Returns:
        bins : (array-like)
    """
    binsize = stats['Std']/5

    if eval_type == 'diff':
        bins = np.arange(0-(4*stats['Std']),
                         0+(4*stats['Std']),
                         binsize)

    else:
        bins = np.arange(stats['Mean']-(4*stats['Std']),
                         stats['Mean']+(4*stats['Std']),
                         binsize)

    return bins


def calculate_stats(data):
    """
    Calculates n, mean, min, max,
    standard deviation, and RMSE.
    Returns dictionary with stats.

    Args:
        data : (array-like) array of data to calculate
               stats
    Returns:
        stats : (dict) dictionary of stats
    """

    n = np.count_nonzero(~np.isnan(data))

    mean = np.nanmean(data)
    std = np.nanstd(data)
    mx = np.nanmax(data)
    mn = np.nanmin(data)

    rmse = np.sqrt(np.nanmean(np.square(data)))

    stats = {'Nobs': n,
             'Min': np.round(mn, 3),
             'Max': np.round(mx, 3),
             'Mean': np.round(mean, 3),
             'Std': np.round(std, 3),
             'RMSE': np.round(rmse, 3)
             }

    return stats


def get_labels(metadata):
    """
    Creates a dictionary of title, date title, and the save
    file name using information from the metadata.

    Args:
        metadata : (dict) metadata dictionary created by diags.py
    Returns:
        labels : (dict) dictionary of labels
    """

    var = metadata['Satellite'] if 'Satellite' in metadata \
        else metadata['Variable']

    # Get title and save file name
    if metadata['Anl Use Type'] is not None:
        title = (f"{metadata['Obs Type']}: {var} - {metadata['Diag Type']}"
                 f" - Data {metadata['Anl Use Type']}")

        save_file = (f"{metadata['Date']:%Y%m%d%H}_{metadata['Obs Type']}_"
                     f"{var}_{metadata['Diag Type']}_"
                     f"{metadata['Anl Use Type']}")

    else:
        title = f"{metadata['Obs Type']}: {var} - {metadata['Diag Type']}"
        save_file = (f"{metadata['Date']:%Y%m%d%H}_{metadata['Obs Type']}_"
                     f"{var}_{metadata['Diag Type']}_")

    # Adds on specific obsid, channel, or layer info to title/save file
    if metadata['Diag File Type'] == 'conventional':
        title = title + '\n%s' % '\n'.join(metadata['ObsID Name'])
        save_file = save_file + '%s' % '_'.join(
            str(x) for x in metadata['ObsID Name'])

    elif metadata['Diag File Type'] == 'radiance':
        if metadata['Channels'] == 'All Channels':
            title = title + '\nAll Channels'
            save_file = save_file + 'All_Channels'

        else:
            title = title + '\nChannels: %s' % ', '.join(
                str(x) for x in metadata['Channels'])
            save_file = save_file + 'channels_%s' % '_'.join(
                str(x) for x in metadata['Channels'])

    else:
        layer = metadata['Layer']
        title = title + f'\nLayer: {layer}'
        save_file = save_file + '%s' % layer

    # Get date label
    date_title = metadata['Date'].strftime("%d %b %Y %Hz")

    labels = {
        'title': title,
        'date title': date_title,
        'save file': save_file
    }

    return labels


def varspecs_name(variable):
    """
    Grabs the specific variable name to utlitize emcpy's VariableSpecs.

    Args:
        variable : (str) variable name
    Returns:
        spec_variable : (str) name of spec variable
    """
    vardict = {
            'temperature': ['air temperature', 'tmp', 'temp', 't'],
            'specific humidity': ['q', 'spfh'],
            'u': ['eastward wind', 'ugrd', 'zonal wind'],
            'v': ['northward wind', 'vgrd', 'meridional wind'],
            'wind speed': ['windspeed'],
            'brightness temperature': ['bt'],
            'integrated layer ozone in air': ['ozone', 'o3']
        }

    # makes lower case and replaces underscore with space
    spec_variable = variable.lower().replace('_', ' ')

    for key in vardict.keys():
        spec_variable = key if spec_variable in vardict[key] \
            else spec_variable

    return spec_variable
