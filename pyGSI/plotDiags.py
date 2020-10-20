#!/usr/bin/env python
## Created by Kevin Dougherty
## October 2020

import numpy as np
from textwrap import TextWrapper
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature

def get_obs_type(obs_id):
    
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
    
    
def get_xlabel(metadata):
    conv_xlabel = {'t'  : "Temperature (K)",
                   'q'  : "Specific Humidity (kg/kg)",
                   'sst': "Sea Surface Temperature (K)",
                   'pw' : "Precipitable Water (mm)",
                   'tcp': "Pressure (hPa)",
                   'u'  : "Windspeed (m/s)",
                   'v'  : "Windspeed (m/s)",
                   'windspeed': "Windspeed (m/s)"
                  }
    
    if metadata['Diag_type'] == 'conv':
        xlabel = conv_xlabel[metadata['Variable']]
    else:
        xlabel = "Brightness Temperature (K)"
        
    return xlabel

def plot_labels(metadata, stats):
    
    # Stats label
    if stats == None:
        t = None
    else:
        if metadata['Diag_type'] == 'conv' and metadata['Variable'] == 'q':
            roundnum = 6
        else:
            roundnum = 3
            
        t = ('n: %s\nstd: %s\nmean: %s\nmax: %s\nmin: %s' % (stats['N'],
                                                             np.round(stats['Std'],roundnum),
                                                             np.round(stats['Mean'],roundnum),
                                                             np.round(stats['Max'],roundnum),
                                                             np.round(stats['Min'],roundnum)))

    # Date label
    date_title = metadata['Date'].strftime("%d %b %Y %Hz")
    
    
    # X label
    if metadata['Data_type'] == 'O-F':
        xlabel = "Observation - Forecast"
        dataType = 'OmF'
    elif metadata['Data_type'] == 'O-A':
        xlabel = "Observation - Analysis"
        dataType = 'OmA'
    elif metadata['Data_type'] == 'H(x)':
        xlabel = get_xlabel(metadata)
        dataType = 'hofx'
    else:
        xlabel = get_xlabel(metadata)
        dataType = metadata['Data_type']
        
        
    # Save file lable
    if metadata['Diag_type'] == 'conv':
        save_file = '{Date:%Y%m%d%H}_{Diag_type}_{Variable}_'.format(**metadata) + '%s_' % dataType + '%s' % '_'.join(str(x) for x in metadata['ObsID'])
    else:
        save_file = '{Date:%Y%m%d%H}_{Diag_type}_{Satellite}_'.format(**metadata) + '%s_' % dataType + 'channels_%s' % '_'.join(str(x) for x in metadata['Channels'])
        
            
    # Left Title label    
    if metadata['Diag_type'] == 'conv':
        convTitle_dict = {'yes': {'left_title': '{Data_type}: {Variable} - Data Assimilated\n'.format(**metadata) + '%s' % '\n'.join(metadata['Obs_Type']),
                                  'save_file' :  save_file + '_assimilated'},
                          'no' : {'left_title': '{Data_type}: {Variable} - Data Monitored\n'.format(**metadata) + '%s' % '\n'.join(metadata['Obs_Type']),
                                  'save_file' :  save_file + '_monitored'},
                          'n/a': {'left_title': '{Data_type}: {Variable}\n'.format(**metadata) + '%s' % '\n'.join(metadata['Obs_Type']),
                                  'save_file' :  save_file}
                         }
        
        left_title = convTitle_dict[metadata['assimilated']]['left_title']
        save_file  = convTitle_dict[metadata['assimilated']]['save_file']

    else:
        left_title = '{Data_type}: {Diag_type} {Satellite}\n'.format(**metadata) + 'Channels: %s' % ' '.join(str(x) for x in metadata['Channels'])
            
             
    labels = {'statText'  : t,
              'xLabel'    : xlabel,
              'leftTitle' : left_title,
              'dateTitle' : date_title,
              'saveFile'  : save_file
             }
             
             
    return labels


def calculate_stats(data):
    
    n = np.count_nonzero(~np.isnan(data))
    
    mean = np.nanmean(data)
    std = np.nanstd(data)
    mx = np.nanmax(data)
    mn = np.nanmin(data)
    
    rmse = np.sqrt(np.nanmean(np.square(data)))
    
    stats = {'N'    : n,
             'Min'  : mn,
             'Max'  : mx,
             'Mean' : mean,
             'Std'  : std,
             'RMSE' : rmse
            }
    
    return stats

def plot_features(dtype, stats, metadata):
    
    convfeatures_dict = {'t'  : {'cmap': 'jet',
                                 'extend': 'both',
                                 'upperbound': np.round(stats['Mean']) + (np.round(stats['Std']*2)/2)*2,
                                 'lowerbound': np.round(stats['Mean']) - (np.round(stats['Std']*2)/2)*2,
                                 'bins': ((np.round(stats['Mean']) + (np.round(stats['Std']*2)/2)*2) -\
                                          (np.round(stats['Mean']) - (np.round(stats['Std']*2)/2)*2))/20},
                         'q'  : {'cmap': 'GnBu',
                                 'extend': 'max',
                                 'upperbound': stats['Std']*5,
                                 'lowerbound': 0,
                                 'bins': (stats['Std']*5)/10},
                         'ps' : {'cmap': 'jet',
                                 'extend': 'both',
                                 'upperbound': 1050,
                                 'lowerbound': 750,
                                 'bins': (1050 - 750)/10},
                         'sst': {'cmap': 'Spectral_r',
                                 'extend': 'both',
                                 'upperbound': np.round(stats['Mean']) + (np.round(stats['Std']*2)/2)*2,
                                 'lowerbound': np.round(stats['Mean']) - (np.round(stats['Std']*2)/2)*2,
                                 'bins': ((np.round(stats['Mean']) + (np.round(stats['Std']*2)/2)*2) -\
                                          (np.round(stats['Mean']) - (np.round(stats['Std']*2)/2)*2))/20},
                         'pw' : {'cmap': 'GnBu',
                                 'extend': 'max',
                                 'upperbound': stats['Std']*5,
                                 'lowerbound': 0,
                                 'bins': (stats['Std']*5)/10},
                         'tcp': {'cmap': 'Reds_r',
                                 'extend': 'both',
                                 'upperbound': 1050,
                                 'lowerbound': 950,
                                 'bins': (1050 - 950)/10},
                         'u'  : {'cmap': 'viridis',
                                 'extend': 'both',
                                 'upperbound': np.round(stats['Std'])*5,
                                 'lowerbound': np.round(stats['Std'])*-5,
                                 'bins': ((np.round(stats['Std'])*5) - (np.round(stats['Std'])*-5))/10},
                         'v'  : {'cmap': 'viridis',
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
    
    # Get cmap, bins, norm and extend for O-F and O-A
    if dtype in ['O-F', 'O-A']:
        cmap = 'bwr'
        
        upperbound = (np.round(stats['Std']*2)/2)*5
        if upperbound == 0:
            upperbound = stats['Std']*5
        
        lowerbound = 0-upperbound
        bins = (upperbound - lowerbound)/10

        norm = mcolors.BoundaryNorm(boundaries=np.arange(lowerbound, upperbound+bins, bins), ncolors=256)
        
        extend='both'
    
    # Get cmap, bins, norm and extend for O-F and O-A
    elif dtype in ['Observation', 'H(x)', 'windspeed']:
        if metadata['Diag_type'] == 'conv':
            cmap = convfeatures_dict[metadata['Variable']]['cmap']
            extend = convfeatures_dict[metadata['Variable']]['extend']
            upperbound = convfeatures_dict[metadata['Variable']]['upperbound']
            lowerbound = convfeatures_dict[metadata['Variable']]['lowerbound']
            bins = convfeatures_dict[metadata['Variable']]['bins']
        else:
            cmap = 'viridis'
            extend = 'both'
            upperbound = np.round(stats['Mean']) + (np.round(stats['Std']*2)/2)*2
            lowerbound = np.round(stats['Mean']) - (np.round(stats['Std']*2)/2)*2
            bins = (upperbound - lowerbound)/20
            
        
        norm = mcolors.BoundaryNorm(boundaries=np.arange(lowerbound, upperbound+bins, bins), ncolors=256)
        
    else:
        cmap = 'bwr'
        
        upperbound = 1
        lowerbound = 0
        bins = 1
        
        norm = mcolors.BoundaryNorm(boundaries=np.arange(lowerbound, upperbound+bins, bins), ncolors=256)
        
        extend = 'neither'
        
    
    return cmap, norm, extend

def no_data_spatial(metadata, outDir='./'):
    fig = plt.figure(figsize=(15,12))
    ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree(central_longitude=0))

    ax.add_feature(cfeature.GSHHSFeature(scale='auto'))
    ax.set_extent([-180, 180, -90, 90])

    stats = None
    labels = plot_labels(metadata, stats)

    ax.text(0,0, 'No Data', fontsize=32, alpha=0.6, ha='center')
    plt.title(labels['leftTitle'], loc='left', fontsize=14)
    plt.title(labels['dateTitle'], loc='right', fontweight='semibold', fontsize=14)
    plt.savefig(outDir+'/%s_spatial.png' % labels['saveFile'], bbox_inches='tight', pad_inches=0.1)
    plt.close('all')
    
    return


def plot_histogram(data, metadata, outDir='./'):
    if metadata['Diag_type'] == 'conv':
        metadata['Obs_Type'] = get_obs_type(metadata['ObsID'])
        
    if metadata['Diag_type'] == 'conv' and metadata['Variable'] == 'uv':
        for i in data.keys():
            fig = plt.figure(figsize=(8,6))
            ax = fig.add_subplot(111)
            
            if len(data[i]) <= 1:            
                if len(data[i]) == 0:
                    stats = None
                    labels = plot_labels(metadata, stats)
                    ax.text(0.5, 0.5, 'No Data', fontsize=32, alpha=0.6, ha='center')
                else:
                    stats = calculate_stats(data[i])
                    labels = plot_labels(metadata, stats)

                    ax.text(0.75,.7, labels['statText'], fontsize=14, transform=ax.transAxes)
                    ax.text(0.5, 0.5, 'Single Observation', fontsize=32, alpha=0.6, ha='center')
                
            else:       
                metadata['Variable'] = i

                stats = calculate_stats(data[i])

                binsize = (stats['Max']-stats['Min'])/np.sqrt(stats['N'])

                if i == 'windspeed':
                    bins = np.arange(0,stats['Mean']+(4*stats['Std']),binsize)
                else:
                    bins = np.arange(0-(4*stats['Std']),0+(4*stats['Std']),binsize)

                plt.hist(data[i], bins=bins)
                plt.axvline(stats['Mean'], color='r', linestyle='solid', linewidth=1)
                if metadata['Data_type'] == 'O-F' or metadata['Data_type'] == 'O-A':
                    plt.axvline(0, color='k', linestyle='dashed', linewidth=1)

                labels = plot_labels(metadata, stats)

                ax.text(0.75,.7, labels['statText'], fontsize=14, transform=ax.transAxes)

            plt.xlabel(labels['xLabel'])
            plt.ylabel('Count')

            wrapper = TextWrapper(width=70,break_long_words=False,replace_whitespace=False)
            leftTitle = '\n'.join(wrapper.wrap(labels['leftTitle']))

            plt.title(leftTitle, loc='left', fontsize=14)
            plt.title(labels['dateTitle'], loc='right', fontweight='semibold', fontsize=14)
            plt.savefig(outDir+'/%s_histogram.png' % labels['saveFile'], bbox_inches='tight', pad_inches=0.1)
            plt.close('all')
            
    else:
        fig = plt.figure(figsize=(8,6))
        ax = fig.add_subplot(111)
        
        if len(data) <= 1:
            if len(data) == 0:
                stats = None
                labels = plot_labels(metadata, stats)
                ax.text(0.5, 0.5, 'No Data', fontsize=32, alpha=0.6, ha='center')
            else:
                stats = calculate_stats(data)
                labels = plot_labels(metadata, stats)
                
                ax.text(0.75,.7, labels['statText'], fontsize=14, transform=ax.transAxes)
                ax.text(0.5, 0.5, 'Single Observation', fontsize=32, alpha=0.6, ha='center')
        
        else:
            stats = calculate_stats(data) 

            binsize = (stats['Max']-stats['Min'])/np.sqrt(stats['N'])
            if metadata['Data_type'] == 'O-F' or metadata['Data_type'] == 'O-A': 
                bins = np.arange(0-(4*stats['Std']),0+(4*stats['Std']),binsize)
            else:
                bins = np.arange(stats['Mean']-(4*stats['Std']),stats['Mean']+(4*stats['Std']),binsize)

            plt.hist(data, bins=bins)
            plt.axvline(stats['Mean'], color='r', linestyle='solid', linewidth=1)
            if metadata['Data_type'] == 'O-F' or metadata['Data_type'] == 'O-A':
                plt.axvline(0, color='k', linestyle='dashed', linewidth=1)

            labels = plot_labels(metadata, stats)

            ax.text(0.75,.7, labels['statText'], fontsize=14, transform=ax.transAxes)

        plt.xlabel(labels['xLabel'])
        plt.ylabel('Count')

        wrapper = TextWrapper(width=70,break_long_words=False,replace_whitespace=False)
        leftTitle = '\n'.join(wrapper.wrap(labels['leftTitle']))

        plt.title(leftTitle, loc='left', fontsize=14)
        plt.title(labels['dateTitle'], loc='right', fontweight='semibold', fontsize=14)
        plt.savefig(outDir+'/%s_histogram.png' % labels['saveFile'], bbox_inches='tight', pad_inches=0.1)
        plt.close('all')
    
    return


def plot_spatial(data, metadata, lats, lons, outDir='./'):
    
    if metadata['Diag_type'] == 'conv':
        metadata['Obs_Type'] = get_obs_type(metadata['ObsID'])
    
    if metadata['Diag_type'] == 'conv' and metadata['Variable'] == 'uv':
        for i in data.keys():
            if len(data[i]) == 0:
                no_data_spatial(metadata, outDir)
                
            else:
                
                metadata['Variable'] = i

                stats = calculate_stats(data[i])

                if len(data) == 1:
                    cmap = 'Blues'
                    norm = None
                    extend = 'neither'
                
                elif str(i) == 'windspeed':              
                    cmap, norm, extend = plot_features('windspeed', stats, metadata)
                else:
                    cmap, norm, extend = plot_features(metadata['Data_type'], stats, metadata)

                fig = plt.figure(figsize=(15,12))
                ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree(central_longitude=0))

                ax.add_feature(cfeature.GSHHSFeature(scale='auto'))
                ax.set_extent([-180, 180, -90, 90])

                cs = plt.scatter(lons, lats, c=data[i], s=30,
                                 norm=norm, cmap=cmap,
                                 transform=ccrs.PlateCarree())

                labels = plot_labels(metadata, stats)

                ax.text(185, 55, labels['statText'], fontsize=14)

                cb = plt.colorbar(cs, orientation='horizontal', shrink=0.5, pad=.04, extend=extend)
                cb.set_label(labels['xLabel'], fontsize=12)

                plt.title(labels['leftTitle'], loc='left', fontsize=14)
                plt.title(labels['dateTitle'], loc='right', fontweight='semibold', fontsize=14)
                plt.savefig(outDir+'/%s_spatial.png' % labels['saveFile'], bbox_inches='tight', pad_inches=0.1)
                plt.close('all')
    
    else:
        
        if len(data) == 0:
            no_data_spatial(metadata, outDir)
        
        else:
            stats = calculate_stats(data)
            
            if len(data) == 1:
                cmap = 'Blues'
                norm = None
                extend = 'neither'
            else:
                cmap, norm, extend = plot_features(metadata['Data_type'], stats, metadata)

            fig = plt.figure(figsize=(15,12))
            ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree(central_longitude=0))

            ax.add_feature(cfeature.GSHHSFeature(scale='auto'))
            ax.set_extent([-180, 180, -90, 90])

            cs = plt.scatter(lons, lats, c=data, s=30,
                             norm=norm, cmap=cmap,
                             transform=ccrs.PlateCarree())

            labels = plot_labels(metadata, stats)

            ax.text(185, 55, labels['statText'], fontsize=14)

            cb = plt.colorbar(cs, orientation='horizontal', shrink=0.5, pad=.04, extend=extend)
            cb.set_label(labels['xLabel'], fontsize=12)

            plt.title(labels['leftTitle'], loc='left', fontsize=14)
            plt.title(labels['dateTitle'], loc='right', fontweight='semibold', fontsize=14)  
            plt.savefig(outDir+'/%s_spatial.png' % labels['saveFile'], bbox_inches='tight', pad_inches=0.1)
            plt.close('all')
    
    return
