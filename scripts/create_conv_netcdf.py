#!/usr/bin/env python

import argparse
import numpy as np
import yaml
import sys
import itertools
from multiprocessing import Pool
from pyGSI.diags import Conventional
from pyGSI.netcdf_diags import write_netcdf
from pyGSI.spatial_bin import spatial_bin
from datetime import datetime

start_time = datetime.now()


def first_occurrence(worklist):

    firstlist = []
    repeatinglist = []
    obsid = None

    for w in worklist:
        if obsid == w['conventional input']['observation id'][0]:
            repeatinglist.append(w)

        else:
            firstlist.append(w)

        obsid = w['conventional input']['observation id'][0]

    return firstlist, repeatinglist


def create_netcdf(conv_config):

    diagfile = conv_config['conventional input']['path'][0]
    diag_type = conv_config['conventional input']['data type'][0].lower()
    obsid = conv_config['conventional input']['observation id']
    subtype = conv_config['conventional input']['observation subtype']
    analysis_use = conv_config['conventional input']['analysis use'][0]
    outdir = conv_config['outdir']

    diag = Conventional(diagfile)

    if analysis_use:
        diag_components = diagfile.split('/')[-1].split('.')[0].split('_')
        if diag_components[1] == 'conv' and diag_components[2] == 'uv':
            u, v = diag.get_data(diag_type, obsid=obsid,
                                 subtype=subtype, analysis_use=analysis_use)

            data = {'u': u['assimilated'],
                    'v': v['assimilated'],
                    'windspeed': np.sqrt(
                        np.square(u['assimilated']) +
                        np.square(v['assimilated']))
                    }

        else:
            data = diag.get_data(diag_type, obsid=obsid,
                                 subtype=subtype, analysis_use=analysis_use)

            data = data['assimilated']

        lats, lons = diag.get_lat_lon(
            obsid=obsid, subtype=subtype, analysis_use=analysis_use)
        lats = lats['assimilated']
        lons = lons['assimilated']

        pressure = diag.get_pressure(
            obsid=obsid, subtype=subtype, analysis_use=analysis_use)
        pressure = pressure['assimilated']

        metadata = diag.metadata

        metadata['Anl Use Type'] = 'assimilated'

        # Get binned data
        if diag_components[1] == 'conv' and diag_components[2] == 'uv':
            binned_data = spatial_bin(
                data, lats, lons, binsize='1x1', uv_data=True,
                pressure=pressure)
        else:
            binned_data = spatial_bin(
                data, lats, lons, binsize='1x1', pressure=pressure)

        write_netcdf(data, binned_data, metadata, outdir)

    else:
        diag_components = diagfile.split('/')[-1].split('.')[0].split('_')
        if diag_components[1] == 'conv' and diag_components[2] == 'uv':
            u, v = diag.get_data(diag_type, obsid=obsid,
                                 subtype=subtype, analysis_use=analysis_use)

            data = {'u': u,
                    'v': v,
                    'windspeed': np.sqrt(np.square(u) + np.square(v))
                    }
        else:
            data = diag.get_data(diag_type, obsid=obsid,
                                 subtype=subtype, analysis_use=analysis_use)

        lats, lons = diag.get_lat_lon(
            obsid=obsid, subtype=subtype, analysis_use=analysis_use)
        pressure = diag.get_pressure(
            obsid=obsid, subtype=subtype, analysis_use=analysis_use)

        metadata = diag.metadata

        # Get binned data
        if diag_components[1] == 'conv' and diag_components[2] == 'uv':
            binned_data = spatial_bin(
                data, lat, lon, binsize='1x1', uv_data=True, pressure=pressure)
        else:
            binned_data = spatial_bin(
                data, lat, lon, binsize='1x1', pressure=pressure)

        write_netcdf(data, binned_data, metadata, outdir)

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

input_yaml = myargs.yaml
outdir = myargs.outdir

with open(input_yaml, 'r') as file:
    parsed_yaml_file = yaml.load(file, Loader=yaml.FullLoader)

worklist = (parsed_yaml_file['diagnostic'])

for w in parsed_yaml_file['diagnostic']:
    w['outdir'] = outdir

condition = True

while condition:
    worklist, repeatinglist = first_occurrence(worklist)

    # run multiprocessing pool with worklist
    p = Pool(processes=nprocs)
    p.map(create_netcdf, worklist)

    if len(repeatinglist) == 0:
        condition = False

    else:
        worklist = repeatinglist
        condition = True

print(datetime.now() - start_time)
