## Created by Kevin Dougherty
## October 2020

import numpy as np
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
    if metadata['Diag_type'] == 'conv' and metadata['Variable'] == 't':
        xlabel = "Temperature (K)"
    elif metadata['Diag_type'] == 'conv' and metadata['Variable'] == 'q':
        xlabel = "Specific Humidity (kg/kg)"
    elif metadata['Diag_type'] == 'conv' and metadata['Variable'] == 'sst':
        xlabel = "Sea Surface Temperature (K)"
    elif metadata['Diag_type'] == 'conv' and metadata['Variable'] == 'pw':
        xlabel = "Precipitable Water (mm)"
    elif metadata['Diag_type'] == 'conv' and metadata['Variable'] == 'tcp':
        xlabel = "Pressure (hPa)"
    else:
        xlabel = 'Radiance'

    return xlabel

def plot_labels(metadata, stats):
    
    if stats == None:
        t = None
    elif metadata['Diag_type'] == 'conv' and metadata['Variable'] == 'q':
        t = ('n: %s\nstd: %s\nmean: %s\nmax: %s\nmin: %s' % (stats['N'],np.round(stats['Std'],6),np.round(stats['Mean'],6), np.round(stats['Max'],6), np.round(stats['Min'],6)))
    else:
        t = ('n: %s\nstd: %s\nmean: %s\nmax: %s\nmin: %s' % (stats['N'],np.round(stats['Std'],3),np.round(stats['Mean'],3), np.round(stats['Max'],3), np.round(stats['Min'],3)))
             
    
    date_title = metadata['Date'].strftime("%d %b %Y %Hz")
    
    
    if metadata['File_type'] == 'ges':
        
        if metadata['Data_type'] == 'O-F':
            
            if metadata['Diag_type'] == 'conv' and metadata['Variable'] == 'windspeed':
                xlabel = "Wind Speed (m/s)"
                save_file  = '{Date:%Y%m%d%H}_{Diag_type}_{Variable}_O_minus_F'.format(**metadata) + '%s' % '_'.join(str(metadata['ObsID']))
            else:
                xlabel = 'Observation - Forecast'
                save_file  = '{Date:%Y%m%d%H}_{Diag_type}_{Variable}_O_minus_F'.format(**metadata) + '%s' % '_'.join(str(metadata['ObsID']))
                
        elif metadata['Data_type'] == 'H(x)':
            xlabel = get_xlabel(metadata)
            save_file  = '{Date:%Y%m%d%H}_{Diag_type}_{Variable}_H_of_x'.format(**metadata) + '%s' % '_'.join(metadata['ObsID'])
            
        else:
            xlabel = get_xlabel(metadata)
            save_file  = '{Date:%Y%m%d%H}_{Diag_type}_{Variable}_{Data_type}'.format(**metadata) + '%s' % '_'.join(metadata['ObsID'])
            
            
        if metadata['Diag_type'] == 'conv':
            left_title = '{Data_type}: {Variable}\n'.format(**metadata) + '%s' % '\n'.join(metadata['Obs_Type'])
        else:
            left_title = '{Data_type}: {Satellite} {Diag_type}'.format(**metadata)
             
    elif metadata['File_type'] == 'anl':
             
        if metadata['Data_type'] == 'O-A':
            
            if metadata['Diag_type'] == 'conv' and metadata['Variable'] == 'windspeed':
                xlabel = "Wind Speed (m/s)"
                save_file  = '{Date:%Y%m%d%H}_{Diag_type}_{Variable}_O_minus_A'.format(**metadata) + '%s' % '_'.join(metadata['ObsID'])
            else:
                xlabel = 'Observation - Forecast'
                save_file  = '{Date:%Y%m%d%H}_{Diag_type}_{Variable}_O_minus_A'.format(**metadata)
                
        elif metadata['Data_type'] == 'H(x)':
            xlabel = get_xlabel(metadata)
            save_file  = '{Date:%Y%m%d%H}_{Diag_type}_{Variable}_H_of_x'.format(**metadata)
            
        else:
            xlabel = get_xlabel(metadata)
            save_file  = '{Date:%Y%m%d%H}_{Diag_type}_{Variable}_{Data_type}'.format(**metadata)
            
            
        if metadata['Diag_type'] == 'conv':
            left_title = '{Data_type}: {Variable}\n'.format(**metadata) + '%s' % '\n'.join(metadata['Obs_Type'])
        else:
            left_title = '{Data_type}: {Satellite} {Diag_type}'.format(**metadata)
            
             
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
    
    if dtype == 'O-F' or dtype == 'O-A':
        cmap = 'bwr'
        
        upperbound = (np.round(stats['Std']*2)/2)*5
        if upperbound == 0:
            upperbound = stats['Std']*5
        
        lowerbound = 0-upperbound
        bins = (upperbound - lowerbound)/10

        norm = mcolors.BoundaryNorm(boundaries=np.arange(lowerbound, upperbound+bins, bins), ncolors=256)
        
        extend='both'
    
    elif dtype == 'Observation' or dtype == 'H(x)' or dtype == 'windspeed':
        if dtype == 'windspeed':
            cmap = 'Spectral_r'
            
            upperbound = np.round(stats['Std'])*5
            lowerbound = 0
            bins = (upperbound - lowerbound)/10
            
            extend = 'max'
            
        elif metadata['Diag_type'] == 'conv' and metadata['Variable'] == 'q':
            cmap = 'GnBu'  
            
            upperbound = stats['Std']*5
            lowerbound = 0
            bins = (upperbound - lowerbound)/10
            
            extend = 'max'
            
        elif metadata['Diag_type'] == 'conv' and metadata['Variable'] == 't':
            cmap = 'jet' 
            
            upperbound = np.round(stats['Mean']) + (np.round(stats['Std']*2)/2)*2
            lowerbound = np.round(stats['Mean']) - (np.round(stats['Std']*2)/2)*2
            bins = (upperbound - lowerbound)/20
            
            extend = 'both'
            
        elif metadata['Diag_type'] == 'conv' and metadata['Variable'] == 'sst':
            cmap = 'Spectral_r'
            
            upperbound = np.round(stats['Mean']) + (np.round(stats['Std']*2)/2)*2
            lowerbound = np.round(stats['Mean']) - (np.round(stats['Std']*2)/2)*2
            bins = (upperbound - lowerbound)/20
            
            extend = 'both'
            
        elif metadata['Diag_type'] == 'conv' and metadata['Variable'] == 'pw':
            cmap = 'GnBu'  
            
            upperbound = stats['Std']*5
            lowerbound = 0
            bins = (upperbound - lowerbound)/10
            
            extend = 'max'
        
        elif metadata['Diag_type'] == 'conv' and metadata['Variable'] == 'tcp':
            cmap = 'Reds_r'  
            
            upperbound = 1050
            lowerbound = 950
            bins = (upperbound - lowerbound)/10
            
            extend = 'both'
            
        elif metadata['Diag_type'] == 'conv' and metadata['Variable'] == 'ps':
            cmap = 'jet'  
            
            upperbound = 1050
            lowerbound = 750
            bins = (upperbound - lowerbound)/50
            
            extend = 'both'
            
        
        norm = mcolors.BoundaryNorm(boundaries=np.arange(lowerbound, upperbound+bins, bins), ncolors=256)
        
    else:
        cmap = 'bwr'
        
        upperbound = 1
        lowerbound = 0
        bins = 1
        
        norm = mcolors.BoundaryNorm(boundaries=np.arange(lowerbound, upperbound+bins, bins), ncolors=256)
        
        extend = 'neither'
        
    
    return cmap, norm, extend


def plot_spatial(data, metadata, lats, lons):
    
    if metadata['Diag_type'] == 'conv':
        metadata['Obs_Type'] = get_obs_type(metadata['ObsID'])
    
    if metadata['Diag_type'] == 'conv' and metadata['Variable'] == 'uv':
        for i in data.keys():
            if len(data[i]) == 0:
                
                fig = plt.figure(figsize=(15,12))
                ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree(central_longitude=0))

                ax.add_feature(cfeature.GSHHSFeature(scale='auto'))
                ax.set_extent([-180, 180, -90, 90])
                
                stats = None
                labels = plot_labels(metadata, stats)

                ax.text(0,0, 'No Data', fontsize=32, alpha=0.6, ha='center')
                plt.title(labels['leftTitle'], loc='left', fontsize=14)
                plt.title(labels['dateTitle'], loc='right', fontweight='semibold', fontsize=14)
                plt.savefig('%s_spatial.png' % labels['saveFile'], bbox_inches='tight', pad_inches=0.1)
                
            else:
                metadata['Variable'] = i

                stats = calculate_stats(data[i])

                if str(i) == 'windspeed':              
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
                plt.savefig('%s_spatial.png' % labels['saveFile'], bbox_inches='tight', pad_inches=0.1)
    
    else:
        
        if len(data) == 0:
            fig = plt.figure(figsize=(15,12))
            ax = fig.add_subplot(1, 1, 1, projection=ccrs.PlateCarree(central_longitude=0))

            ax.add_feature(cfeature.GSHHSFeature(scale='auto'))
            ax.set_extent([-180, 180, -90, 90])
            
            stats = None
            labels = plot_labels(metadata, stats)
            
            ax.text(0,0, 'No Data', fontsize=32, alpha=0.6, ha='center')
            plt.title(labels['leftTitle'], loc='left', fontsize=14)
            plt.title(labels['dateTitle'], loc='right', fontweight='semibold', fontsize=14)  
            plt.savefig('%s_spatial.png' % labels['saveFile'], bbox_inches='tight', pad_inches=0.1)
        
        else:
            stats = calculate_stats(data)
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
            plt.savefig('%s_spatial.png' % labels['saveFile'], bbox_inches='tight', pad_inches=0.1)
    
    return