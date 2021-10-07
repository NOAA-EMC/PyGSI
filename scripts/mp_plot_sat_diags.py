#!/usr/bin/env python

import argparse
import numpy as np
import yaml
from multiprocessing import Pool
import sys
from pyGSI.diags import Radiance
from pyGSI.plot_diags import plot_map, plot_histogram
from datetime import datetime

start_time = datetime.now()


def plotting(sat_config):

    diagfile = sat_config['radiance input']['path'][0]
    diag_type = sat_config['radiance input']['data type'][0].lower()
    channel = sat_config['radiance input']['channel']
    qcflag = sat_config['radiance input']['qc flag']
    analysis_use = sat_config['radiance input']['analysis use'][0]
    plot_type = sat_config['radiance input']['plot type']
    outdir = sat_config['outdir']

    diag = Radiance(diagfile)

    df = diag.get_data(channel=channel, qcflag=qcflag,
                       analysis_use=analysis_use)
    metadata = diag.metadata
    metadata['Diag Type'] = diag_type

    column = f'{diag_type}' if diag_type in ['observation'] \
        else f'{diag_type}_adjusted'

    if analysis_use:
        lats = {
            'assimilated': df['assimilated']['latitude'].to_numpy(),
            'rejected': df['rejected']['latitude'].to_numpy(),
            'monitored': df['monitored']['latitude'].to_numpy()
        }
        lons = {
            'assimilated': df['assimilated']['longitude'].to_numpy(),
            'rejected': df['rejected']['longitude'].to_numpy(),
            'monitored': df['monitored']['longitude'].to_numpy()
        }

        data = {
            'assimilated': df['assimilated'][column].to_numpy(),
            'rejected': df['rejected'][column].to_numpy(),
            'monitored': df['monitored'][column].to_numpy()
        }

        for key in data.keys():
            data[key][data[key] > 1e5] = np.nan

    else:
        lats = df['latitude'].to_numpy()
        lons = df['longitude'].to_numpy()

        data = df[column].to_numpy()

        data[data > 1e5] = np.nan

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
