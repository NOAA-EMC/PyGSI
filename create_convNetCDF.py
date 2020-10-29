#!/usr/bin/env python

import argparse
import numpy as np
import yaml
import sys
from pyGSI.Diags import conventional
from pyGSI.netcdfDiags import writeNetCDF
from datetime import datetime

startTime = datetime.now()

def createNetCDF(YAML):
    for y in YAML:
        diagFile = y['conventional input']['path'][0]
        DataType = y['conventional input']['data type'][0]
        ObsID    = y['conventional input']['observation id']
        Subtype  = y['conventional input']['observation subtype']
        AnlUse   = y['conventional input']['analysis use'][0]
        plotType = y['conventional input']['plot type']
        outDir   = y['outDir']

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

            metadata = diag.get_metadata()

            metadata['Data_type'] = DataType   
            metadata['ObsID'] = ObsID
            metadata['Subtype'] = Subtype
            metadata['outDir'] = outDir

            metadata['assimilated'] = 'yes'
            lat = lats['assimilated']
            lon = lons['assimilated']


            writeNetCDF(data, metadata, lat, lon) 

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

            metadata = diag.get_metadata()

            metadata['Data_type'] = DataType   
            metadata['ObsID'] = ObsID
            metadata['Subtype'] = Subtype
            metadata['outDir'] = outDir

            writeNetCDF(data, metadata, lats, lons)
    
    return
        
        
###############################################
    
# Parse command line
ap = argparse.ArgumentParser()  
ap.add_argument("-y", "--yaml",
                help="Path to yaml file with diag data")
ap.add_argument("-o", "--outdir",
                help="Out directory where files will be saved")

MyArgs = ap.parse_args()

YAML = MyArgs.yaml
outDir = MyArgs.outdir
    
file = open(YAML)
parsed_yaml_file = yaml.load(file, Loader=yaml.FullLoader)

for w in parsed_yaml_file['diagnostic']:
    w['outDir'] = outDir

createNetCDF(parsed_yaml_file['diagnostic'])
    
print(datetime.now() - startTime)