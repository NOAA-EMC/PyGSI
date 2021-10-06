import argparse
import numpy as np
import yaml
from multiprocessing import Pool
import sys
from pyGSI.diags import Ozone
from pyGSI.plot_diags import plot_spatial, plot_histogram
from datetime import datetime

start_time = datetime.now()


def plotting(ozone_config):

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

    if np.isin('histogram', plot_type):
        plot_histogram(data[dict_key], metadata, outdir)
    if np.isin('spatial', plot_type):
        plot_spatial(data[dict_key], metadata,
                     lats[dict_key], lons[dict_key], outdir)


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
