#!/usr/bin/env python

import argparse
import numpy as np
import yaml
from multiprocessing import Pool
import sys
from pyGSI.Diags import radiance
from plotting.plotDiags import plot_spatial, plot_histogram
from datetime import datetime

startTime = datetime.now()

def plotting(YAML):
    
    diagFile = YAML['radiance input']['path'][0]
    DataType = YAML['radiance input']['data type'][0]
    Channel  = YAML['radiance input']['channel']
    QCFlag   = YAML['radiance input']['qc flag'][0]
    plotType = YAML['radiance input']['plot type']
    
    plotType = ['histogram', 'spatial']

    diag = radiance(diagFile)
    
    data = diag.getData(DataType, channel=Channel, qcflag=QCFlag)
    lats, lons = diag.get_lat_lon(channel=Channel, qcflag=QCFlag)
    
    metadata = diag.get_metadata()
    metadata['Data_type'] = DataType  
    metadata['Channels'] = Channel
    
    if np.isin('histogram', plotType):
        plot_histogram(data, metadata)
    if np.isin('spatial', plotType):
        plot_spatial(data, metadata, lats, lons)
    

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

MyArgs = ap.parse_args()

if MyArgs.nprocs:
    nprocs = int(MyArgs.nprocs)
else:
    nprocs = 1

YAML = MyArgs.yaml
    
file = open(YAML)
parsed_yaml_file = yaml.load(file, Loader=yaml.FullLoader)

work = (parsed_yaml_file['diagnostic'])

pool_handler(nprocs)
    
print(datetime.now() - startTime)