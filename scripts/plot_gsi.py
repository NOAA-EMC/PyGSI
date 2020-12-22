import numpy as np
import sys
import argparse
import yaml
from pyGSI.diags import conventional, radiance
from pyGSI.plot_diags import plot_spatial, plot_histogram
from datetime import datetime


def plotting_conventional(YAML):

    diagfile = YAML['conventional input']['path'][0]
    data_type = YAML['conventional input']['data type'][0]
    obsid = YAML['conventional input']['observation id']
    analysis_use = YAML['conventional input']['analysis use'][0]
    plot_type = YAML['conventional input']['plot type']
    outdir = YAML['outDir']

    diag = conventional(diagfile)

    if analysis_use == True:
        diag_components = diagfile.split('/')[-1].split('.')[0].split('_')
        if diag_components[1] == 'conv' and diag_components[2] == 'uv':
            u, v = diag.get_data(data_type, obsid=obsid,
                                 analysis_use=analysis_use)

            assimilated_data = {'u': u['assimilated'],
                                'v': v['assimilated'],
                                'windspeed': np.sqrt(np.square(u['assimilated']) + np.square(v['assimilated']))
                                }

            monitored_data = {'u': u['monitored'],
                              'v': v['monitored'],
                              'windspeed': np.sqrt(np.square(u['monitored']) + np.square(v['monitored']))
                              }
        else:
            data = diag.get_data(data_type, obsid=obsid,
                                 analysis_use=analysis_use)

            assimilated_data = data['assimilated']
            monitored_data = data['monitored']

        lats, lons = diag.get_lat_lon(obsid=obsid, analysis_use=analysis_use)

        for i, data in enumerate([assimilated_data, monitored_data]):
            for plot in plot_type:
                metadata = diag.get_metadata()

                metadata['Data_type'] = data_type
                metadata['ObsID'] = obsid

                if i == 0:
                    metadata['assimilated'] = 'yes'
                    lat = lats['assimilated']
                    lon = lons['assimilated']
                else:
                    metadata['assimilated'] = 'no'
                    lat = lats['monitored']
                    lon = lons['monitored']

                if plot == 'histogram':
                    plot_histogram(data, metadata, outdir)
                if plot == 'spatial':
                    plot_spatial(data, metadata, lat, lon, outdir)

    else:

        diag_components = diagfile.split('/')[-1].split('.')[0].split('_')
        if diag_components[1] == 'conv' and diag_components[2] == 'uv':
            u, v = diag.get_data(data_type, obsid=obsid,
                                 analysis_use=analysis_use)
            data = {'u': u,
                    'v': v,
                    'windspeed': np.sqrt(np.square(u) + np.square(v))
                    }
        else:
            data = diag.get_data(data_type, obsid=obsid)

        lats, lons = diag.get_lat_lon(obsid=obsid)

        metadata = diag.get_metadata()

        metadata['Data_type'] = data_type
        metadata['ObsID'] = obsid
        metadata['assimilated'] = 'n/a'

        if np.isin('histogram', plot_type):
            plot_histogram(data, metadata, outdir)
        if np.isin('spatial', plot_type):
            plot_spatial(data, metadata, lats, lons, outdir)


def plotting_radiance(YAML):

    diagfile = YAML['radiance input']['path'][0]
    data_type = YAML['radiance input']['data type'][0]
    channel = YAML['radiance input']['channel']
    qcflag = YAML['radiance input']['qc flag']
    plot_type = YAML['radiance input']['plot type']
    outdir = YAML['outDir']

    diag = radiance(diagfile)

    data = diag.get_data(data_type, channel=channel, qcflag=qcflag)
    lats, lons = diag.get_lat_lon(channel=channel, qcflag=qcflag)

    metadata = diag.get_metadata()
    metadata['Data_type'] = data_type
    metadata['Channels'] = channel

    if np.isin('histogram', plot_type):
        plot_histogram(data, metadata, outdir)
    if np.isin('spatial', plot_type):
        plot_spatial(data, metadata, lats, lons, outdir)


###############################################

# Parse command line
ap = argparse.ArgumentParser()
ap.add_argument("-y", "--yaml",
                help="Path to yaml file with diag data",
                required=True)
ap.add_argument("-o", "--outdir",
                help="Out directory where files will be saved",
                default="./")

myargs = ap.parse_args()

YAML = myargs.yaml
outdir = myargs.outdir

if outdir == None:
    outdir = './'

file = open(YAML)
parsed_yaml_file = yaml.load(file, Loader=yaml.FullLoader)

work = (parsed_yaml_file['diagnostic'])
for w in work:
    w['outDir'] = outdir

diagType = next(iter(parsed_yaml_file['diagnostic'][0].items()))[0]
if diagType == 'conventional input':
    plotting_conventional(work[0])
elif diagType == 'radiance input':
    plotting_radiance(work[0])
else:
    print('YAML entry incorrect. Please correct.')
