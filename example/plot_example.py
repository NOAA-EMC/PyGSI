import numpy as np
from datetime import datetime
import yaml
import sys
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import argparse

sys.path.append('/home/Kevin.Dougherty/GSI_plots')

from pyGSI.diags import conventional, satellite


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
            save_file  = '{Date:%Y%m%d%H}_{Diag_type}_{Variable}_O_minus_F.png'.format(**metadata)
        else:
            left_title = '{Data_type}: {Satellite}'.format(**metadata)
            save_file  = '{Date:%Y%m%d%H}_{Diag_type}_{Satellite}_O_minus_F.png'.format(**metadata)
             
    elif metadata['File_type'] == 'anl':
             
        xlabel     = " ".join(metadata['Data_type'])
        date_title = metadata['Date'].strftime("%d %b %Y %Hz")
             
        if metadata['Diag_type'] == 'conv':
            left_title = '{Data_type}: {Variable}'.format(**metadata)
            save_file  = '{Date:%Y%m%d%H}_{Diag_type}_{Variable}_O_minus_A.png'.format(**metadata)
        else:
            left_title = '{Data_type}: {Satellite}'.format(**metadata)
            save_file  = '{Date:%Y%m%d%H}_{Diag_type}_{Satellite}_O_minus_A.png'.format(**metadata)
             
    labels = {'statText'  : t,
              'xLabel'    : xlabel,
              'leftTitle' : left_title,
              'dateTitle' : date_title,
              'saveFile'  : save_file
             }
             
             
    return labels


def calculate_stats(data):
    
    n = len(data)
    
    mean = np.mean(data)
    var_list = [(x-mean)**2 for x in data]
    variance = np.sum(var_list)/(len(var_list)-1)
    std = np.sqrt(variance)
    mx = max(data)
    mn = min(data)
    
    rmse = np.sqrt(np.mean(np.square(data)))
    
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
    plt.savefig(labels['saveFile'], bbox_inches='tight', pad_inches=0.1)
    
    return


def main(parsed_yaml_file):
    if parsed_yaml_file['conventional input']['path']:
        
        nc_file   = parsed_yaml_file['conventional input']['path'][0]
        obs_id    = parsed_yaml_file['conventional input']['observation id']
        qc_flag   = parsed_yaml_file['conventional input']['qc flag']
        DATA_TYPE = parsed_yaml_file['conventional input']['data type'][0]
        
        diag = conventional(nc_file)
        
        data = diag.getData(DATA_TYPE, obs_id, qc_flag)
        
    elif parsed_yaml_file['satellite input']['path']:
        
        nc_file   = parsed_yaml_file['satellite input']['path']
        channel   = parsed_yaml_file['satellite input']['channel']
        qc_flag   = parsed_yaml_file['satellite input']['qc flag']
        DATA_TYPE = parsed_yaml_file['satellite input']['data type'][0]
        
        diag = satellite(nc_file)
        
        data = diag.getData(DATA_TYPE, channel, qc_flag)
        
    else:
        print('File type not recognized. Please address in yaml file.')
        
    
    metadata = diag.get_metadata()
    metadata['Data_type'] = DATA_TYPE
    print(metadata)
    
    plot_histogram(data, metadata)
        
    
    
    
    return

#########################################################

parser = argparse.ArgumentParser(description='Get basic statistics from GSI diagnostic file and save as .csv file.')          
parser.add_argument('-f', '--file', type=str,                                                                                                        
                    help='path to .yaml file.', required=True)

args = parser.parse_args()

f = args.file

file = open(f)
parsed_yaml_file = yaml.load(file, Loader=yaml.FullLoader)
main(parsed_yaml_file)