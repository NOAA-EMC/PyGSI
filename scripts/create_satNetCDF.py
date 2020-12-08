#!/usr/bin/env python

import argparse
import numpy as np
import yaml
import sys
import itertools
from multiprocessing import Pool
from pyGSI.Diags import radiance
from pyGSI.netcdfDiags import writeNetCDF
from pyGSI.spatialBin import spatialBin
from datetime import datetime

startTime = datetime.now()

def createNetCDF(YAML):
    
    diagFile = YAML['radiance input']['path'][0]
    DataType = YAML['radiance input']['data type'][0]
    Channel  = YAML['radiance input']['channel']
    QCFlag   = YAML['radiance input']['qc flag']
    plotType = YAML['radiance input']['plot type']   
    outDir   = YAML['outDir']

    diag = radiance(diagFile)
    
    data = diag.getData(DataType, channel=Channel, qcflag=QCFlag)
    lats, lons = diag.get_lat_lon(channel=Channel, qcflag=QCFlag)
    
    metadata = diag.get_metadata()
    metadata['Data_type'] = DataType  
    metadata['Channel'] = Channel
    metadata['outDir'] = outDir
    
    binnedData = spatialBin(data, lats, lons, binsize='1x1')
    
    writeNetCDF(data, binnedData, metadata)
    
    return

    
###############################################
    
#Parse command line
ap = argparse.ArgumentParser()  
ap.add_argument("-n", "--nprocs",
                help="Number of tasks/processors for multiprocessing")  
ap.add_argument("-y", "--yaml",
                help="Path to yaml file with diag data")
ap.add_argument("-o", "--outdir",
                help="Out directory where files will be saved")

MyArgs = ap.parse_args()

if MyArgs.nprocs:
    nprocs = int(MyArgs.nprocs)
else:
    nprocs = 1

YAML = MyArgs.yaml
outDir = MyArgs.outdir
    
file = open(YAML)
parsed_yaml_file = yaml.load(file, Loader=yaml.FullLoader)

for w in parsed_yaml_file['diagnostic']:
    w['outDir'] = outDir

work = (parsed_yaml_file['diagnostic'])
    
p = Pool(processes=nprocs)
p.map(createNetCDF, work)

print(datetime.now() - startTime)