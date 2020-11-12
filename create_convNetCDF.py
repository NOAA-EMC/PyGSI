#!/usr/bin/env python

import argparse
import numpy as np
import yaml
import sys
import itertools
from multiprocessing import Pool
from pyGSI.Diags import conventional
from pyGSI.netcdfDiags import writeNetCDF
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


def spatialBin(data, lat, lon):
    # Creates new 1x1 global map
    ll_lats = np.linspace(-89.5,89.5,180)
    ll_lons = np.linspace(0.5,359.5,360)
    
    #Generate a list of tuples where each tuple is a combination of parameters.
    #The list will contain all possible combinations of parameters.
    paramlist = list(itertools.product(ll_lats, ll_lons))
    
    # Create 2D arrays for each statistic
    binned_nobs = np.zeros((180,360))
    binned_mean = np.zeros((180,360))
    binned_max = np.zeros((180,360))
    binned_min = np.zeros((180,360))
    binned_std = np.zeros((180,360))
    binned_rmse = np.zeros((180,360))
    
    # Loop through all combinations in paramlist to find data values within
    # a 1x1 lat lon map. Get appropriate indexes and calculate statistics and
    # append them to above 2D arrays
    for val in paramlist:
        # Get index of 1x1 grid lat and lon
        latidx = np.where(ll_lats == val[0])
        lonidx = np.where(ll_lons == val[1])
        # values of the 1x1 grid lat and lon
        binnedlons = val[1]
        binnedlats = val[0]

        n_deg=0.5

        # find instances where data is within 1x1 grid point of orginal data
        data_idx = np.where((lon >= binnedlons - n_deg) & (lon <= binnedlons + n_deg) & \
                            (lat >= binnedlats - n_deg) & (lat <= binnedlats + n_deg))

        latlon_idx = [latidx[0][0],lonidx[0][0]]

        # calculate stats if there is data at this grid point, else append np.nan
        if len(data_idx[0]) > 0:
            try:
                n = len(data[data_idx])
            except:
                print(data[data_idx])
            mean = np.mean(data[data_idx])
            mx = np.max(data[data_idx])
            mn = np.min(data[data_idx])
            std = np.std(data[data_idx])
            rmse = np.std(data[data_idx])

            binned_nobs[latlon_idx[0], latlon_idx[1]] = n
            binned_mean[latlon_idx[0], latlon_idx[1]] = mean
            binned_max[latlon_idx[0], latlon_idx[1]] = mx
            binned_min[latlon_idx[0], latlon_idx[1]] = mn
            binned_std[latlon_idx[0], latlon_idx[1]] = std
            binned_rmse[latlon_idx[0], latlon_idx[1]] = rmse

        else:
            binned_nobs[latlon_idx[0], latlon_idx[1]] = np.nan
            binned_mean[latlon_idx[0], latlon_idx[1]] = np.nan
            binned_max[latlon_idx[0], latlon_idx[1]] = np.nan
            binned_min[latlon_idx[0], latlon_idx[1]] = np.nan
            binned_std[latlon_idx[0], latlon_idx[1]] = np.nan
            binned_rmse[latlon_idx[0], latlon_idx[1]] = np.nan
            
    binnedData = {'binned_nobs': binned_nobs,
                  'binned_mean': binned_mean,
                  'binned_max': binned_max,
                  'binned_min': binned_min,
                  'binned_std': binned_std,
                  'binned_rmse': binned_rmse
                 }

    return binnedData


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

        metadata = diag.get_metadata()

        metadata['Data_type'] = DataType   
        metadata['ObsID'] = ObsID
        metadata['Subtype'] = Subtype
        metadata['outDir'] = outDir

        metadata['assimilated'] = 'yes'
        lat = lats['assimilated']
        lon = lons['assimilated']

        # Get binned data
        print(diagFile, ObsID, Subtype)
        binnedData = spatialBin(data, lat, lon)


        writeNetCDF(data, binnedData, metadata, lat, lon) 

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

        # Get binned data
        binnedData = spatialBin(data, lats, lons)

        writeNetCDF(data, binnedData, metadata, lats, lons)
    
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
#     print('Working list: ', len(worklist))
#     print("Repeating list: ", len(repeatinglist))
    #run multiprocessing pool with worklist
    p = Pool(processes=nprocs)
    p.map(createNetCDF, worklist)
    
    if len(repeatinglist) == 0:
        condition = False
        
    else:
        worklist = repeatinglist
        condition = True
    
print(datetime.now() - startTime)