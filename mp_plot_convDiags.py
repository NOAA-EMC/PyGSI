#!/usr/bin/env python

import argparse
import numpy as np
import yaml
from multiprocessing import Pool
import sys
from pyGSI.Diags import conventional
from plotting.plotDiags import plot_spatial, plot_histogram
from datetime import datetime

startTime = datetime.now()

def plotting(YAML):
    
    diagFile = YAML['conventional input']['path'][0]
    DataType = YAML['conventional input']['data type'][0]
    ObsID    = YAML['conventional input']['observation id']
    AnlUse   = YAML['conventional input']['analysis use'][0]
    plotType = YAML['conventional input']['plot type']

    
    diag = conventional(diagFile)
    
#     print(AnlUse)
    
    if AnlUse == True:
        diagComponents = diagFile.split('/')[-1].split('.')[0].split('_')
        if diagComponents[1] == 'conv' and diagComponents[2] == 'uv':
            u,v = diag.getData(DataType, obsid=ObsID, analysis_use=AnlUse)

            assimilated_data = {'u': u['assimilated'],
                                'v': v['assimilated'],
                                'windspeed': np.sqrt(np.square(u['assimilated']) + np.square(v['assimilated']))
                               }

            monitored_data = {'u': u['monitored'],
                              'v': v['monitored'],
                              'windspeed': np.sqrt(np.square(u['monitored']) + np.square(v['monitored']))
                             }
        else:
            data = diag.getData(DataType, obsid=ObsID, analysis_use=AnlUse)

            assimilated_data = data['assimilated']
            monitored_data = data['monitored']

        lats, lons = diag.get_lat_lon(obsid=ObsID, analysis_use=AnlUse)

        for i, data in enumerate([assimilated_data, monitored_data]):
            for plot in plotType:
                metadata = diag.get_metadata()

                metadata['Data_type'] = DataType   
                metadata['ObsID'] = ObsID

                if i == 0:
                    metadata['assimilated'] = 'yes'
                    lat = lats['assimilated']
                    lon = lons['assimilated']
                else:
                    metadata['assimilated'] = 'no'
                    lat = lats['monitored']
                    lon = lons['monitored']
                    
                if plot == 'histogram':
                    plot_histogram(data, metadata)
                if plot == 'spatial':
                    plot_spatial(data, metadata, lat, lon)
 

    else:
        
        diagComponents = diagFile.split('/')[-1].split('.')[0].split('_')
        if diagComponents[1] == 'conv' and diagComponents[2] == 'uv':
            u,v = diag.getData(DataType, obsid=ObsID, analysis_use=AnlUse)
            data = {'u': u,
                    'v': v,
                    'windspeed': np.sqrt(np.square(u) + np.square(v))
                   }
        else:
            data = diag.getData(DataType, obsid=ObsID)

        lats, lons = diag.get_lat_lon(obsid=ObsID)

        metadata = diag.get_metadata()

        metadata['Data_type'] = DataType   
        metadata['ObsID'] = ObsID
        metadata['assimilated'] = 'n/a'
            
        if np.isin('histogram', plotType):
            plot_histogram(data, metadata)
        if np.isin('spatial', plotType):
            plot_spatial(data, metadata, lat, lon)

    
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
