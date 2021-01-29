#!/usr/bin/env python

import argparse
import numpy as np
import yaml
import sys
import itertools
from multiprocessing import Pool
from pyGSI.diags import Ozone
from pyGSI.netcdf_diags import write_netcdf
from pyGSI.spatial_bin import spatial_bin
from datetime import datetime

start_time = datetime.now()


def create_netcdf(ozone_config):
    
    diagfile = ozone_config['ozone input']['path'][0]
    diag_type = ozone_config['ozone input']['data type'][0].lower()
    analysis_use = ozone_config['ozone input']['analysis use'][0]
    layer = ozone_config['ozone input']['layer'][0]
    plot_type = ozone_config['ozone input']['plot type']
    outdir = ozone_config['outdir']
    
    diag = Ozone(diagfile)

    data = diag.get_data(diag_type, analysis_use=analysis_use)
        
    lats, lons = diag.get_lat_lon(analysis_use=analysis_use)
        
    if layer == 0:
        dict_key = 'column total'
    else:
        dict_key = list(data)[layer]

    metadata = diag.metadata
    metadata['Layer'] = dict_key
    
    if analysis_use:
        binned_data = spatial_bin(data[dict_key]['assimilated'], lats[dict_key]['assimilated'], lons[dict_key]['assimilated'], binsize='1x1')

        write_netcdf(data[dict_key]['assimilated'], binned_data, metadata, outdir)
    
    else:
        binned_data = spatial_bin(data[dict_key], lats[dict_key], lons[dict_key], binsize='1x1')

        write_netcdf(data[dict_key], binned_data, metadata, outdir)

    return

###############################################

# Parse command line
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

for w in parsed_yaml_file['diagnostic']:
    w['outdir'] = outdir

work = (parsed_yaml_file['diagnostic'])

p = Pool(processes=nprocs)
p.map(create_netcdf, work)

print(datetime.now() - start_time)