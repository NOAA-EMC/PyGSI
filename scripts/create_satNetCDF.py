#!/usr/bin/env python

import argparse
import numpy as np
import yaml
import sys
import itertools
from multiprocessing import Pool
from pyGSI.Diags import radiance
from pyGSI.netcdfDiags import writeNetCDF
from datetime import datetime

startTime = datetime.now()

def binData(data, lat, lon, binsize='1x1'):
    """
    The main function to spatially bin the data.
    Inputs:
        data      : data to binned
        lat       : original data lats
        lon       : original data lons
        binsize   : string of the size of binning, must be square.
                    Examples: '1x1', '2x2', '5x5' (Default = 1x1)
    Outputs:
        binned_data
    """
    
    # Create lats and lons based on binsize
    lonlen = 360
    latlen = 180

    lon_lowerlim = 0
    lon_upperlim = 360

    lat_lowerlim = -90
    lat_upperlim = 90

    if binsize.split('x')[0] != binsize.split('x')[1]:
        print('ERROR: Binsize must be square i.e. 1x1, 2x2, 5x5 etc. Please use different binsize.')
        return

    binsize = int(binsize.split('x')[0])

    if latlen % binsize == 0 and lonlen % binsize == 0:
        latbin = int(latlen/binsize)
        lonbin = int(lonlen/binsize)
        n_deg = binsize/2

        ll_lats = np.linspace(lat_lowerlim+(n_deg),
                              lat_upperlim-(n_deg),
                              latbin)

        ll_lons = np.linspace(lon_lowerlim+(n_deg),
                              lon_upperlim-(n_deg),
                              lonbin)

    else:
        print('ERROR: Binsize does not work for grid shape (180,360). Please use different binsize.')
        return

    paramlist = list(itertools.product(ll_lats, ll_lons))
    
    # Bin Data
    binned_data = np.full((latbin,lonbin), np.nan, dtype=object)


    for val in paramlist:
        # Get index of grid lat and lon
        latidx = np.where(ll_lats == val[0])
        lonidx = np.where(ll_lons == val[1])
        # values of the 1x1 grid lat and lon
        binnedlons = val[1]
        binnedlats = val[0]

        # find instances where data is within 1x1 grid point of orginal data
        data_idx = np.where((lon >= binnedlons - n_deg) & (lon <= binnedlons + n_deg) & \
                            (lat >= binnedlats - n_deg) & (lat <= binnedlats + n_deg))

        latlon_idx = [latidx[0][0],lonidx[0][0]]

        # calculate stats if there is data at this grid point
        if len(data_idx[0]) > 0:
            d = data[data_idx]
            binned_data[latlon_idx[0], latlon_idx[1]] = d


    return binned_data



def spatialBin(data, lat, lon, binsize='1x1'):

    binned_data = binData(data, lat, lon, binsize=binsize)
    
    rows = binned_data.shape[0]
    cols = binned_data.shape[1]

    binned_nobs = np.full((rows, cols), np.nan)
    binned_mean = np.full((rows, cols), np.nan)
    binned_max = np.full((rows, cols), np.nan)
    binned_min = np.full((rows, cols), np.nan)
    binned_std = np.full((rows, cols), np.nan)
    binned_rmse = np.full((rows, cols), np.nan)

    for x in range(0, rows):
        for y in range(0, cols):
            if np.isnan(binned_data[x,y]).any() == False:
                binned_nobs[x,y] = len(binned_data[x,y])
                binned_mean[x,y] = np.mean(binned_data[x,y])
                binned_max[x,y] = np.max(binned_data[x,y])
                binned_min[x,y] = np.min(binned_data[x,y])
                binned_std[x,y] = np.std(binned_data[x,y])
                binned_rmse[x,y] = np.sqrt(np.nanmean(np.square(binned_data[x,y])))

    binnedData = {'binned_nobs': binned_nobs,
                  'binned_mean': binned_mean,
                  'binned_max': binned_max,
                  'binned_min': binned_min,
                  'binned_std': binned_std,
                  'binned_rmse': binned_rmse
                 }

    return binnedData

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