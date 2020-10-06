import numpy as np
from datetime import datetime
import yaml
import sys
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import argparse

sys.path.append('../pyGSI')

from Diags import conventional, radiance


def plot_labels(metadata, stats):
    
    if metadata['Diag_type'] == 'conv' and metadata['Variable'] == 'q':
        t = ('n: %s\nstd: %s\nmean: %s\nmax: %s\nmin: %s' % (stats['N'],np.round(stats['Std'],6),np.round(stats['Mean'],6), np.round(stats['Max'],6), np.round(stats['Min'],6)))
    else:
        t = ('n: %s\nstd: %s\nmean: %s\nmax: %s\nmin: %s' % (stats['N'],np.round(stats['Std'],3),np.round(stats['Mean'],3), np.round(stats['Max'],3), np.round(stats['Min'],3)))
             
    
    if metadata['File_type'] == 'ges':
             
        xlabel     = " ".join(metadata['Data_type'])
        date_title = metadata['Date'].strftime("%d %b %Y %Hz")
             
        if metadata['Diag_type'] == 'conv':
            left_title = '{Data_type}: {Variable}'.format(**metadata)
            save_file  = '{Date:%Y%m%d%H}_{Diag_type}_{Variable}_O_minus_F'.format(**metadata)
        else:
            left_title = '{Data_type}: {Satellite}'.format(**metadata)
            save_file  = '{Date:%Y%m%d%H}_{Diag_type}_{Satellite}_O_minus_F'.format(**metadata)
             
    elif metadata['File_type'] == 'anl':
             
        xlabel     = " ".join(metadata['Data_type'])
        date_title = metadata['Date'].strftime("%d %b %Y %Hz")
             
        if metadata['Diag_type'] == 'conv':
            left_title = '{Data_type}: {Variable}'.format(**metadata)
            save_file  = '{Date:%Y%m%d%H}_{Diag_type}_{Variable}_O_minus_A'.format(**metadata)
        else:
            left_title = '{Data_type}: {Satellite}'.format(**metadata)
            save_file  = '{Date:%Y%m%d%H}_{Diag_type}_{Satellite}_O_minus_A'.format(**metadata)
             
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


def plot_histogram(data, metadata):
    
    stats = calculate_stats(data) 
    
    binsize = (stats['Max']-stats['Min'])/np.sqrt(stats['N'])
    bins = np.arange(0-(4*stats['Std']),0+(4*stats['Std']),binsize)
    
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111)
    plt.hist(data, bins=bins)
             
    labels = plot_labels(metadata, stats)
             
    ax.text(0.75,.7, labels['statText'], fontsize=14, transform=ax.transAxes)
             
    plt.xlabel(labels['xLabel'])
    plt.ylabel('Count')
    
    title_split = labels['leftTitle'].split('\n')
    plt.title(labels['leftTitle'], loc='left', fontsize=14)
    plt.title(labels['dateTitle'], loc='right', fontweight='semibold', fontsize=14)
    plt.savefig('%s_histogram.png' % labels['saveFile'], bbox_inches='tight', pad_inches=0.1)
    
    return


def plot_spatial(data, metadata, lats, lons):
    
    stats = calculate_stats(data)
    
    plt.figure(figsize=(15,12))
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
    ax.add_feature(cfeature.GSHHSFeature(scale='auto'))
    ax.set_extent([-180, 180, -90, 90])
    
    upperbound = (np.round(stats['Std']*2)/2)*5
    lowerbound = 0-upperbound
    bins = (upperbound - lowerbound)/10
    
    norm = mcolors.BoundaryNorm(boundaries=np.arange(lowerbound, upperbound+bins, bins), ncolors=256)
    
    cs = plt.scatter(lons, lats, c=data, s=30,
                norm=norm, cmap='bwr', #edgecolors='gray', linewidth=0.25,
                transform=ccrs.PlateCarree())
    
    labels = plot_labels(metadata, stats)
             
    ax.text(-175, -70, labels['statText'], fontsize=14, transform=ccrs.PlateCarree())
    
    cb = plt.colorbar(cs, shrink=0.5, pad=.04, extend='both')
    cb.set_label(labels['xLabel'], fontsize=12)
    
    plt.title(labels['leftTitle'], loc='left', fontsize=14)
    plt.title(labels['dateTitle'], loc='right', fontweight='semibold', fontsize=14)
    plt.savefig('%s_spatial.png' % labels['saveFile'], bbox_inches='tight', pad_inches=0.1)

    
    return


def main(parsed_yaml_file):
    
    for group in parsed_yaml_file['diagnostic']:
        for groupType in group.keys():
            
            if groupType == 'conventional input':
        
                nc_file   = group[groupType]['path'][0]
                obs_id    = group[groupType]['observation id']
                qc_flag   = group[groupType]['qc flag']
                DATA_TYPE = group[groupType]['data type'][0]

                diag = conventional(nc_file)

                data = diag.getData(DATA_TYPE, obs_id, qc_flag)

                lats, lons = diag.get_lat_lon(obs_id, qc_flag)
        
            elif groupType == 'radiance input':
                
                nc_file   = group[groupType]['path'][0]
                channel   = group[groupType]['channel']
                qc_flag   = group[groupType]['qc flag']
                DATA_TYPE = group[groupType]['data type'][0]

                diag = radiance(nc_file)

                data = diag.getData(DATA_TYPE, channel, qc_flag)

                lats, lons = diag.get_lat_lon(channel, qc_flag)
        
            else:
                print('File type not recognized. Please address in yaml file.')
                return
        
    
            metadata = diag.get_metadata()
            metadata['Data_type'] = DATA_TYPE

            plot_histogram(data, metadata)

            plot_spatial(data, metadata, lats, lons)
    
    
    return

#########################################################

parser = argparse.ArgumentParser(description='Plot GSI diagnostics histogram and spatial plots.')          
parser.add_argument('-f', '--file', type=str,                                                                                                        
                    help='path to .yaml file.', required=True)

args = parser.parse_args()

f = args.file

file = open(f)
try:
    parsed_yaml_file = yaml.load(file, Loader=yaml.FullLoader)
except AttributeError:
    parsed_yaml_file = yaml.load(file)
main(parsed_yaml_file)
