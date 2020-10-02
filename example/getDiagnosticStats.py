import numpy as np
from datetime import datetime
import yaml
import csv
import sys
import argparse

sys.path.append('/home/Kevin.Dougherty/GSI_plots')

from runDiagnostics import conventional, satellite


def calculate_stats(data):
    
    n = len(data)
    
    mean = np.mean(data)
    var_list = [(x-mean)**2 for x in data]
    variance = np.sum(var_list)/(len(var_list)-1)
    std = np.sqrt(variance)
    mx = max(data)
    mn = min(data)
    
    rmse = np.sqrt(np.mean(np.square(data)))
    
    return mn, mx, mean, std, rmse

def write_csv(nc_file, csv_dict):
    
    fields = ['Date', 'min', 'max', 'mean', 'std', 'rmse']
    
    # filename of csv file
    filename = csv_dict[0]['Date'] + '_' + nc_file.split('/')[-1].split('.')[0] +'.csv'
    
    # write to csv file
    with open(filename, 'w') as csvfile:
        # creating a csv dict writer object
        writer = csv.DictWriter(csvfile, fieldnames=fields)

        # writing headers (filed names)
        writer.writeheader()

        # writing data rows
        writer.writerows(csv_dict)
        
        return
        

def main(parsed_yaml_file):
    """
    Main function that reads diagnostic files, calculates min, max, mean,
    standard deviation, and RMSE, then outputs results to csv. CSV can be
    concatenated in shell script for multiple dates to one csv that can 
    be read in by pandas dataframe and plotted.
    
    CSV is saved as:
        
        date + diagnostic file name + .csv
        
    Need to add output file. Make default '/.'
    """
    
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
            
            elif groupType == 'radiance input':
                
                nc_file   = group[groupType]['path'][0]
                channel   = group[groupType]['channel']
                qc_flag   = group[groupType]['qc flag']
                DATA_TYPE = group[groupType]['data type'][0]

                diag = satellite(nc_file)

                data = diag.getData(DATA_TYPE, channel, qc_flag)

                lats, lons = diag.get_lat_lon(channel, qc_flag)
        
            else:
                print('File type not recognized. Please address in yaml file.')
                return
            
    
            metadata = diag.get_metadata()

            mn, mx, mean, std, rmse = calculate_stats(data)

            date = metadata['Date'].strftime("%Y%m%d%H")

            # dictionary of values needed in csv
            csv_dict = [{'Date' : date,
                        'min'  : mn,
                        'max'  : mx,
                        'mean' : mean,
                        'std'  : std,
                        'rmse' : rmse
                       }]

            write_csv(nc_file, csv_dict)
        
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
