#!/usr/bin/env python

import argparse
import numpy as np
import yaml
from multiprocessing import Pool
import sys
from pyGSI.diags import Conventional
from pyGSI.plot_diags import plot_map, plot_histogram
from datetime import datetime

start_time = datetime.now()


def _uv_data(df, diag_type, analysis_use, lats, lons,
             metadata, plot_type, outdir):

    # Creates column name fo u and v data
    u = f'u_{diag_type}' if diag_type in ['observation'] \
        else f'u_{diag_type}_adjusted'
    v = f'v_{diag_type}' if diag_type in ['observation'] \
        else f'v_{diag_type}_adjusted'

    if analysis_use:
        data = {
            'u': {
                'assimilated': df['assimilated'][u].to_numpy(),
                'rejected': df['rejected'][u].to_numpy(),
                'monitored': df['monitored'][u].to_numpy()
            },
            'v': {
                'assimilated': df['assimilated'][v].to_numpy(),
                'rejected': df['rejected'][v].to_numpy(),
                'monitored': df['monitored'][v].to_numpy()
            }
        }

    else:
        data = {
            'u': df[u].to_numpy(),
            'v': df[v].to_numpy(),
        }

    for key in data.keys():
        metadata['Variable'] = key
        if np.isin('histogram', plot_type):
            plot_histogram(data[key], metadata, outdir)
        if np.isin('spatial', plot_type):
            plot_map(lats, lons, data[key], metadata, outdir)


def plotting(conv_config):

    diagfile = conv_config['conventional input']['path'][0]
    diag_type = conv_config['conventional input']['data type'][0].lower()
    obsid = conv_config['conventional input']['observation id']
    obsid = None if obsid in [['None'], ['none']] else obsid
    analysis_use = conv_config['conventional input']['analysis use'][0]
    plot_type = conv_config['conventional input']['plot type']
    outdir = conv_config['outdir']

    diag = Conventional(diagfile)

    if analysis_use:

        df_dict = diag.get_data(obsid=obsid, analysis_use=analysis_use)
        metadata = diag.metadata
        metadata['Diag Type'] = diag_type

        # Get lat/lons
        lats = {
            'assimilated': df_dict['assimilated']['latitude'].to_numpy(),
            'rejected': df_dict['rejected']['latitude'].to_numpy(),
            'monitored': df_dict['monitored']['latitude'].to_numpy()
        }
        lons = {
            'assimilated': df_dict['assimilated']['longitude'].to_numpy(),
            'rejected': df_dict['rejected']['longitude'].to_numpy(),
            'monitored': df_dict['monitored']['longitude'].to_numpy()
        }

        if metadata['Obs Type'] == 'conv' and metadata['Variable'] == 'uv':
            _uv_data(df_dict, diag_type, analysis_use, lats,
                     lons, metadata, plot_type, outdir)

        else:
            column = f'{diag_type}' if diag_type in ['observation'] \
                else f'{diag_type}_adjusted'

            data = {
                'assimilated': df_dict['assimilated'][column].to_numpy(),
                'rejected': df_dict['rejected'][column].to_numpy(),
                'monitored': df_dict['monitored'][column].to_numpy()
            }

            if np.isin('histogram', plot_type):
                plot_histogram(data, metadata, outdir)
            if np.isin('spatial', plot_type):
                plot_map(lats, lons, data, metadata, outdir)

    else:

        df = diag.get_data(diag_type, obsid=obsid,
                           analysis_use=analysis_use)
        metadata = diag.metadata
        metadata['Diag Type'] = diag_type

        lats = df['latitude'].to_numpy()
        lons = df['longitude'].to_numpy()

        if metadata['Obs Type'] == 'conv' and metadata['Variable'] == 'uv':
            _uv_data(df_dict, diag_type, analysis_use, lats,
                     lons, metadata, plot_type, outdir)

        else:
            column = f'{diag_type}' if diag_type in ['observation'] \
                else f'{diag_type}_adjusted'
            data = df[column].to_numpy()

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

for w in work:
    w['outdir'] = outdir

p = Pool(processes=nprocs)
p.map(plotting, work)

print(datetime.now() - start_time)
