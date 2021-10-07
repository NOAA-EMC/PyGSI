import argparse
import numpy as np
import yaml
from multiprocessing import Pool
import sys
from pyGSI.diags import Ozone
from pyGSI.plot_diags import plot_map, plot_histogram
from datetime import datetime

start_time = datetime.now()


def plotting(ozone_config):

    diagfile = ozone_config['ozone input']['path'][0]
    diag_type = ozone_config['ozone input']['data type'][0].lower()
    analysis_use = ozone_config['ozone input']['analysis use'][0]
    layer = ozone_config['ozone input']['layer'][0]
    plot_type = ozone_config['ozone input']['plot type']
    outdir = ozone_config['outdir']

    diag = Ozone(diagfile)

    df = diag.get_data(analysis_use=analysis_use)
    metadata = diag.metadata
    metadata['Diag Type'] = diag_type
    
    column = f'{diag_type}' if diag_type in ['observation'] \
        else f'{diag_type}_adjusted'

    if layer == 0:
        dict_key = 'column total'
    else:
        dict_key = list(df)[layer]
    metadata['Layer'] = dict_key
    
    if analysis_use:
        lats = {
            'assimilated': df[dict_key]['assimilated']['latitude'].to_numpy(),
            'monitored': df[dict_key]['monitored']['latitude'].to_numpy()
        }
        lons = {
            'assimilated': df[dict_key]['assimilated']['longitude'].to_numpy(),
            'monitored': df[dict_key]['monitored']['longitude'].to_numpy()
        }
        
        data = {
            'assimilated': df[dict_key]['assimilated'][column].to_numpy(),
            'monitored': df[dict_key]['monitored'][column].to_numpy()
        }
        
    else:
        lats = df[dict_key]['latitude'].to_numpy()
        lons = df[dict_key]['longitude'].to_numpy()
    
        data = df[dict_key][column].to_numpy()

    if np.isin('histogram', plot_type):
        plot_histogram(data, metadata, outdir)
    if np.isin('spatial', plot_type):
        plot_map(lats, lons, data, metadata, outdir)


###############################################

# Parse command line
ap = argparse.ArgumentParser()
ap.add_argument("-n", "--nprocs",
                help="Number of tasks/processors for multiprocessing")
ap.add_argument("-y", "--yaml",
                help="Path to yaml file with diag data")
ap.add_argument("-o", "--outdir",
                help="Out directory where files will be saved")

myargs = ap.parse_args()

if myargs.nprocs:
    nprocs = int(myargs.nprocs)
else:
    nprocs = 1

input_yaml = myargs.yaml
outdir = myargs.outdir

with open(input_yaml, 'r') as file:
    parsed_yaml_file = yaml.load(file, Loader=yaml.FullLoader)


work = (parsed_yaml_file['diagnostic'])

# Add outdir to yaml dict
for w in work:
    w['outdir'] = outdir

p = Pool(processes=nprocs)
p.map(plotting, work)

print(datetime.now() - start_time)
