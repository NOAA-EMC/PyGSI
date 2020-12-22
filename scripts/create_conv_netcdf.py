#!/usr/bin/env python

import argparse
import numpy as np
import yaml
import sys
import itertools
from multiprocessing import Pool
from pyGSI.diags import conventional
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


def create_netcdf(YAML):

    diagfile = YAML['conventional input']['path'][0]
    data_type = YAML['conventional input']['data type'][0]
    obsid = YAML['conventional input']['observation id']
    subtype = YAML['conventional input']['observation subtype']
    analysis_use = YAML['conventional input']['analysis use'][0]
    outdir = YAML['outDir']

    diag = conventional(diagfile)

    if analysis_use == True:
        diag_components = diagfile.split('/')[-1].split('.')[0].split('_')
        if diag_components[1] == 'conv' and diag_components[2] == 'uv':
            u, v = diag.get_data(data_type, obsid=obsid,
                                 subtype=subtype, analysis_use=analysis_use)

            data = {'u': u['assimilated'],
                    'v': v['assimilated'],
                    'windspeed': np.sqrt(np.square(u['assimilated']) + np.square(v['assimilated']))
                    }

        else:
            data = diag.get_data(data_type, obsid=obsid,
                                 subtype=subtype, analysis_use=analysis_use)

            data = data['assimilated']

        lats, lons = diag.get_lat_lon(
            obsid=obsid, subtype=subtype, analysis_use=analysis_use)
        pressure = diag.get_pressure(
            obsid=obsid, subtype=subtype, analysis_use=analysis_use)
        pressure = pressure['assimilated']

        metadata = diag.get_metadata()

        metadata['Data_type'] = data_type
        metadata['Obsid'] = obsid
        metadata['Subtype'] = subtype
        metadata['outDir'] = outdir

        metadata['assimilated'] = 'yes'
        lat = lats['assimilated']
        lon = lons['assimilated']

        # Get binned data
        if diag_components[1] == 'conv' and diag_components[2] == 'uv':
            binned_data = spatialBin(
                data, lat, lon, binsize='1x1', uvData=True, pressure=pressure)
        else:
            binned_data = spatialBin(
                data, lat, lon, binsize='1x1', pressure=pressure)

        write_netcdf(data, binned_data, metadata)

    else:
        diag_components = diagfile.split('/')[-1].split('.')[0].split('_')
        if diag_components[1] == 'conv' and diag_components[2] == 'uv':
            u, v = diag.get_data(data_type, obsid=obsid,
                                 subtype=subtype, analysis_use=analysis_use)

            data = {'u': u,
                    'v': v,
                    'windspeed': np.sqrt(np.square(u) + np.square(v))
                    }
        else:
            data = diag.get_data(data_type, obsid=obsid,
                                 subtype=subtype, analysis_use=analysis_use)

        lats, lons = diag.get_lat_lon(
            obsid=obsid, subtype=subtype, analysis_use=analysis_use)
        pressure = diag.get_pressure(
            obsid=obsid, subtype=subtype, analysis_use=analysis_use)

        metadata = diag.get_metadata()

        metadata['Data_type'] = data_type
        metadata['obsid'] = obsid
        metadata['subtype'] = subtype
        metadata['outDir'] = outdir

        # Get binned data
        print('Binning..')
        if diag_components[1] == 'conv' and diag_components[2] == 'uv':
            binned_data = spatial_bin(
                data, lat, lon, binsize='1x1', uvData=True, pressure=pressure)
        else:
            binned_data = spatial_bin(
                data, lat, lon, binsize='1x1', pressure=pressure)

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

worklist = (parsed_yaml_file['diagnostic'])

for w in parsed_yaml_file['diagnostic']:
    w['outDir'] = outdir

condition = True

while condition == True:
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
