#!/usr/bin/env python

import argparse
import numpy as np
import yaml
import sys
import itertools
from multiprocessing import Pool
from pyGSI.Diags import conventional
from pyGSI.netcdfDiags import writeNetCDF
from pyGSI.spatialBin import spatialBin
from datetime import datetime

startTime = datetime.now()

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


def createNetCDF(YAML):
    
    diagFile = YAML['conventional input']['path'][0]
    DataType = YAML['conventional input']['data type'][0]
    ObsID    = YAML['conventional input']['observation id']
    Subtype  = YAML['conventional input']['observation subtype']
    AnlUse   = YAML['conventional input']['analysis use'][0]
    plotType = YAML['conventional input']['plot type']
    outDir   = YAML['outDir']

    diag = conventional(diagFile)

    if AnlUse == True:
        diagComponents = diagFile.split('/')[-1].split('.')[0].split('_')
        if diagComponents[1] == 'conv' and diagComponents[2] == 'uv':
            u,v = diag.getData(DataType, obsid=ObsID, subtype=Subtype, analysis_use=AnlUse)

            data = {'u': u['assimilated'],
                    'v': v['assimilated'],
                    'windspeed': np.sqrt(np.square(u['assimilated']) + np.square(v['assimilated']))
                   }

        else:
            data = diag.getData(DataType, obsid=ObsID, subtype=Subtype, analysis_use=AnlUse)

            data = data['assimilated']

        lats, lons = diag.get_lat_lon(obsid=ObsID, subtype=Subtype, analysis_use=AnlUse)
        pressure = diag.get_pressure(obsid=ObsID, subtype=Subtype, analysis_use=AnlUse)
        pressure = pressure['assimilated']

        metadata = diag.get_metadata()

        metadata['Data_type'] = DataType   
        metadata['ObsID'] = ObsID
        metadata['Subtype'] = Subtype
        metadata['outDir'] = outDir

        metadata['assimilated'] = 'yes'
        lat = lats['assimilated']
        lon = lons['assimilated']

        # Get binned data
        if diagComponents[1] == 'conv' and diagComponents[2] == 'uv':
            binnedData = spatialBin(data, lat, lon, binsize='1x1', uvData=True, pressure=pressure)
        else:
            binnedData = spatialBin(data, lat, lon, binsize='1x1', pressure=pressure)


        writeNetCDF(data, binnedData, metadata) 

    else:
        diagComponents = diagFile.split('/')[-1].split('.')[0].split('_')
        if diagComponents[1] == 'conv' and diagComponents[2] == 'uv':
            u,v = diag.getData(DataType, obsid=ObsID, subtype=Subtype, analysis_use=AnlUse)

            data = {'u': u,
                    'v': v,
                    'windspeed': np.sqrt(np.square(u) + np.square(v))
                   }
        else:
            data = diag.getData(DataType, obsid=ObsID, subtype=Subtype, analysis_use=AnlUse)

        lats, lons = diag.get_lat_lon(obsid=ObsID, subtype=Subtype, analysis_use=AnlUse)
        pressure = diag.get_pressure(obsid=ObsID, subtype=Subtype, analysis_use=AnlUse)

        metadata = diag.get_metadata()

        metadata['Data_type'] = DataType   
        metadata['ObsID'] = ObsID
        metadata['Subtype'] = Subtype
        metadata['outDir'] = outDir

        # Get binned data
        print('Binning..')
        if diagComponents[1] == 'conv' and diagComponents[2] == 'uv':
            binnedData = spatialBin(data, lat, lon, binsize='1x1', uvData=True, pressure=pressure)
        else:
            binnedData = spatialBin(data, lat, lon, binsize='1x1', pressure=pressure)

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

worklist = (parsed_yaml_file['diagnostic'])

for w in parsed_yaml_file['diagnostic']:
    w['outDir'] = outDir

condition = True

while condition == True:
    worklist, repeatinglist = first_occurrence(worklist)

    #run multiprocessing pool with worklist
    p = Pool(processes=nprocs)
    p.map(createNetCDF, worklist)
    
    if len(repeatinglist) == 0:
        condition = False
        
    else:
        worklist = repeatinglist
        condition = True
    
print(datetime.now() - startTime)