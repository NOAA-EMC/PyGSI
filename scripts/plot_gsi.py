import numpy as np
import sys
import argparse
import yaml
from pyGSI.diags import Conventional, Radiance
from pyGSI.plot_diags import plot_spatial, plot_histogram
from datetime import datetime


def plotting_conventional(conv_config):

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
                    'monitored':   {'u': u['monitored'],
                                    'v': v['monitored'],
                                    'windspeed': np.sqrt(np.square(u['monitored']) + np.square(v['monitored']))
                                   }
                   }

        else:
            data = diag.get_data(diag_type, obsid=obsid,
                                 analysis_use=analysis_use)
            
            data = {'assimilated': data['assimilated'],
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
            
        return


def plotting_radiance(sat_config):

    diagfile = sat_config['radiance input']['path'][0]
    diag_type = sat_config['radiance input']['data type'][0].lower()
    channel = sat_config['radiance input']['channel']
    qcflag = sat_config['radiance input']['qc flag']
    plot_type = sat_config['radiance input']['plot type']
    outdir = sat_config['outdir']

    diag = Radiance(diagfile)

    data = diag.get_data(diag_type, channel=channel, qcflag=qcflag)
    lats, lons = diag.get_lat_lon(channel=channel, qcflag=qcflag)

    metadata = diag.metadata

    if np.isin('histogram', plot_type):
        plot_histogram(data, metadata, outdir)
    if np.isin('spatial', plot_type):
        plot_spatial(data, metadata, lats, lons, outdir)
        
    return


###############################################

# # Parse command line
ap = argparse.ArgumentParser()
ap.add_argument("-y", "--yaml",
                help="Path to yaml file with diag data",
                required=True)
ap.add_argument("-o", "--outdir",
                help="Out directory where files will be saved",
                default="./")

myargs = ap.parse_args()

input_yaml = myargs.yaml
outdir = myargs.outdir

if outdir == None:
    outdir = './'

with open(input_yaml, 'r') as file:
    parsed_yaml_file = yaml.load(file, Loader=yaml.FullLoader)

work = (parsed_yaml_file['diagnostic'])
for w in work:
    w['outdir'] = outdir

diagType = next(iter(parsed_yaml_file['diagnostic'][0].items()))[0]
if diagType == 'conventional input':
    plotting_conventional(work[0])
elif diagType == 'radiance input':
    plotting_radiance(work[0])
else:
    print('YAML entry incorrect. Please correct.')
