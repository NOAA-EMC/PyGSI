#!/usr/bin/env python

import argparse
import numpy as np
import yaml
from multiprocessing import Pool
import sys
from pyGSI.diags import radiance
from pyGSI.plot_diags import plot_spatial, plot_histogram
from datetime import datetime

start_time = datetime.now()


def plotting(YAML):

    diagfile = YAML['radiance input']['path'][0]
    data_type = YAML['radiance input']['data type'][0]
    channel = YAML['radiance input']['channel']
    qcflag = YAML['radiance input']['qc flag']
    plot_type = YAML['radiance input']['plot type']
    outdir = YAML['outdir']

    diag = radiance(diagfile)

    data = diag.get_data(data_type, channel=channel, qcflag=qcflag)
    lats, lons = diag.get_lat_lon(channel=channel, qcflag=qcflag)

    metadata = diag.get_metadata()
    metadata['Data_type'] = data_type
    metadata['channels'] = channel

    if np.isin('histogram', plot_type):
        plot_histogram(data, metadata, outdir)
    if np.isin('spatial', plot_type):
        plot_spatial(data, metadata, lats, lons, outdir)


def work_log(YAML):
    plotting(YAML)


def pool_handler(nprocs):
    p = Pool(processes=nprocs)
    p.map(work_log, work)


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


work = (parsed_yaml_file['diagnostic'])

# Add outdir to yaml dict
for w in work:
    w['outdir'] = outdir

pool_handler(nprocs)

print(datetime.now() - start_time)
