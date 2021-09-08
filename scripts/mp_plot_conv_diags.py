#!/usr/bin/env python

import argparse
import numpy as np
import yaml
from multiprocessing import Pool
import sys
from pyGSI.diags import Conventional
from pyGSI.plot_diags import plot_spatial, plot_histogram
from datetime import datetime

start_time = datetime.now()


def plotting(conv_config):

    diagfile = conv_config['conventional input']['path'][0]
    diag_type = conv_config['conventional input']['data type'][0].lower()
    obsid = conv_config['conventional input']['observation id']
    analysis_use = conv_config['conventional input']['analysis use'][0]
    plot_type = conv_config['conventional input']['plot type']
    outdir = conv_config['outdir']

    diag = Conventional(diagfile)

    if analysis_use:
        diag_components = diagfile.split('/')[-1].split('.')[0].split('_')
        if diag_components[1] == 'conv' and diag_components[2] == 'uv':
            u, v = diag.get_data(diag_type, obsid=obsid,
                                 analysis_use=analysis_use)
            
            data = {'assimilated': {'u': u['assimilated'],
                                    'v': v['assimilated'],
                                    'windspeed': np.sqrt(np.square(u['assimilated']) + np.square(v['assimilated']))
                                   },
                    'rejected':    {'u': u['rejected'],
                                    'v': v['rejected'],
                                    'windspeed': np.sqrt(np.square(u['rejected']) + np.square(v['rejected']))
                                   },
                    'monitored':   {'u': u['monitored'],
                                    'v': v['monitored'],
                                    'windspeed': np.sqrt(np.square(u['monitored']) + np.square(v['monitored']))
                                   }
                   }

        else:
            data = diag.get_data(diag_type, obsid=obsid,
                                 analysis_use=analysis_use)
            
            data = {'assimilated': data['assimilated'],
                    'rejected': data['rejected'],
                    'monitored': data['monitored']
                   }

        lats, lons = diag.get_lat_lon(obsid=obsid, analysis_use=analysis_use)
        
        metadata = diag.metadata

        if np.isin('histogram', plot_type):
            plot_histogram(data, metadata, outdir)
        if np.isin('spatial', plot_type):
            plot_spatial(data, metadata, lats, lons, outdir)

    else:

        diag_components = diagfile.split('/')[-1].split('.')[0].split('_')
        if diag_components[1] == 'conv' and diag_components[2] == 'uv':
            u, v = diag.get_data(diag_type, obsid=obsid,
                                 analysis_use=analysis_use)
            data = {'u': u,
                    'v': v,
                    'windspeed': np.sqrt(np.square(u) + np.square(v))
                    }
        else:
            data = diag.get_data(diag_type, obsid=obsid)

        lats, lons = diag.get_lat_lon(obsid=obsid)

        metadata = diag.metadata

        if np.isin('histogram', plot_type):
            plot_histogram(data, metadata, outdir)
        if np.isin('spatial', plot_type):
            plot_spatial(data, metadata, lats, lons, outdir)


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
