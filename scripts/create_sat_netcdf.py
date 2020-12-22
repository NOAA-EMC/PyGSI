#!/usr/bin/env python

import argparse
import numpy as np
import yaml
import sys
import itertools
from multiprocessing import Pool
from pyGSI.diags import radiance
from pyGSI.netcdf_diags import write_netcdf
from pyGSI.spatial_bin import spatial_bin
from datetime import datetime

start_time = datetime.now()


def create_netcdf(YAML):

    diagfile = YAML['radiance input']['path'][0]
    data_type = YAML['radiance input']['data type'][0]
    channel = YAML['radiance input']['channel']
    qcflag = YAML['radiance input']['qc flag']
    outdir = YAML['outDir']

    diag = radiance(diagfile)

    data = diag.get_data(data_type, channel=channel, qcflag=qcflag)
    lats, lons = diag.get_lat_lon(channel=channel, qcflag=qcflag)

    metadata = diag.get_metadata()
    metadata['Data_type'] = data_type
    metadata['Channel'] = channel
    metadata['outDir'] = outdir

    binned_data = spatial_bin(data, lats, lons, binsize='1x1')

    write_netcdf(data, binned_data, metadata)

    return


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

YAML = myargs.yaml
outdir = myargs.outdir

file = open(YAML)
parsed_yaml_file = yaml.load(file, Loader=yaml.FullLoader)

for w in parsed_yaml_file['diagnostic']:
    w['outDir'] = outdir

work = (parsed_yaml_file['diagnostic'])

p = Pool(processes=nprocs)
p.map(create_netcdf, work)

print(datetime.now() - start_time)
