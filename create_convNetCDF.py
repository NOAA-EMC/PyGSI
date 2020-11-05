#!/usr/bin/env python

import argparse
import numpy as np
import yaml
import sys
import itertools
from pyGSI.Diags import conventional
from pyGSI.netcdfDiags import writeNetCDF
from datetime import datetime

startTime = datetime.now()

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
            n = len(data[data_idx])
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
            
            # Get binned data
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