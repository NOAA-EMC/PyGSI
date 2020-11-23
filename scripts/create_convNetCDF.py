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


def binData(data, lat, lon, pressure, binsize='1x1', uvData=False):
    """
    The main function to spatially bin the data.
    Inputs:
        data      : data to binned
        lat       : original data lats
        lon       : original data lons
        pressure  : original data pressure values
        paramlist : list of all combinations of lat and lon values
                    of the new grid data is being binned to
        binsize   : string of the size of binning, must be square.
                    Examples: '1x1', '2x2', '5x5' (Default = 1x1)
        uvData    : if using uvData, will be True (Defalt = False) 
    Outputs:
        binned_data, binned_pressure
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
    if uvData == True:
        binned_u_data = np.full((latbin,lonbin), np.nan, dtype=object)
        binned_v_data = np.full((latbin,lonbin), np.nan, dtype=object)
        binned_pressure = np.full((latbin,lonbin), np.nan, dtype=object)
        
        for val in paramlist:
            # Get index of 1x1 grid lat and lon
            latidx = np.where(ll_lats == val[0])
            lonidx = np.where(ll_lons == val[1])
            # values of the 1x1 grid lat and lon
            binnedlons = val[1]
            binnedlats = val[0]

            # find instances where data is within 1x1 grid point of orginal data
            data_idx = np.where((lon >= binnedlons - n_deg) & (lon <= binnedlons + n_deg) & \
                                (lat >= binnedlats - n_deg) & (lat <= binnedlats + n_deg))

            latlon_idx = [latidx[0][0],lonidx[0][0]]

            # calculate stats if there is data at this grid point, else append np.nan        
            if len(data_idx[0]) > 0:
                u = data['u'][data_idx]
                v = data['v'][data_idx]
                p = pressure[data_idx]

                binned_u_data[latlon_idx[0], latlon_idx[1]] = u
                binned_v_data[latlon_idx[0], latlon_idx[1]] = v
                binned_pressure[latlon_idx[0], latlon_idx[1]] = p
                
                
        return binned_u_data, binned_v_data, binned_pressure
        
    else:
        binned_data = np.full((latbin,lonbin), np.nan, dtype=object)
        binned_pressure = np.full((latbin,lonbin), np.nan, dtype=object)
    
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
                try:
                    p = pressure[data_idx]
                except:
                    print(data_idx)
                    print(pressure)

                binned_data[latlon_idx[0], latlon_idx[1]] = d
                binned_pressure[latlon_idx[0], latlon_idx[1]] = p
            

        return binned_data, binned_pressure


def spatialBin(data, lat, lon, pressure, binsize='1x1', uvData=False):
    
    if uvData == True:
        binned_u_data, binned_v_data, binned_pressure = binData(data, lat, lon,
                                                                pressure,
                                                                binsize=binsize,
                                                                uvData=True)
        
        rows = binned_u_data.shape[0]
        cols = binned_u_data.shape[1]

        binned_u_nobs = np.full((rows,cols, 9), np.nan)
        binned_u_mean = np.full((rows,cols, 9), np.nan)
        binned_u_max = np.full((rows,cols, 9), np.nan)
        binned_u_min = np.full((rows,cols, 9), np.nan)
        binned_u_std = np.full((rows,cols, 9), np.nan)
        binned_u_rmse = np.full((rows,cols, 9), np.nan)
        
        binned_v_nobs = np.full((rows,cols, 9), np.nan)
        binned_v_mean = np.full((rows,cols, 9), np.nan)
        binned_v_max = np.full((rows,cols, 9), np.nan)
        binned_v_min = np.full((rows,cols, 9), np.nan)
        binned_v_std = np.full((rows,cols, 9), np.nan)
        binned_v_rmse = np.full((rows,cols, 9), np.nan)

        pressure_list = [None, 0, 100, 250, 500, 700, 850, 925, 1000, 1100]

        for i, pressure in enumerate(pressure_list[:-1]):
            for x in range(0, rows):
                for y in range(0, cols):
                    if np.isnan(binned_u_data[x,y]).any() == False:
                        if i == 0:
                            binned_u_nobs[x,y,i] = len(binned_u_data[x,y])
                            binned_u_mean[x,y,i] = np.mean(binned_u_data[x,y])
                            binned_u_max[x,y,i] = np.max(binned_u_data[x,y])
                            binned_u_min[x,y,i] = np.min(binned_u_data[x,y])
                            binned_u_std[x,y,i] = np.std(binned_u_data[x,y])
                            binned_u_rmse[x,y,i] = np.sqrt(np.nanmean(np.square(binned_u_data[x,y])))
                            
                            binned_v_nobs[x,y,i] = len(binned_v_data[x,y])
                            binned_v_mean[x,y,i] = np.mean(binned_v_data[x,y])
                            binned_v_max[x,y,i] = np.max(binned_v_data[x,y])
                            binned_v_min[x,y,i] = np.min(binned_v_data[x,y])
                            binned_v_std[x,y,i] = np.std(binned_v_data[x,y])
                            binned_v_rmse[x,y,i] = np.sqrt(np.nanmean(np.square(binned_v_data[x,y])))
                        else:
                            pressure_idx = np.where((binned_pressure[x,y] > pressure_list[i]) & (binned_pressure[x,y] < pressure_list[i+1]))
                            if len(pressure_idx[0]) > 0:
                                binned_u_nobs[x,y,i] = len(binned_u_data[x,y][pressure_idx])
                                binned_u_mean[x,y,i] = np.mean(binned_u_data[x,y][pressure_idx])
                                binned_u_max[x,y,i] = np.max(binned_u_data[x,y][pressure_idx])
                                binned_u_min[x,y,i] = np.min(binned_u_data[x,y][pressure_idx])
                                binned_u_std[x,y,i] = np.std(binned_u_data[x,y][pressure_idx])
                                binned_u_rmse[x,y,i] = np.sqrt(np.nanmean(np.square(binned_u_data[x,y][pressure_idx])))
                                
                                binned_v_nobs[x,y,i] = len(binned_v_data[x,y][pressure_idx])
                                binned_v_mean[x,y,i] = np.mean(binned_v_data[x,y][pressure_idx])
                                binned_v_max[x,y,i] = np.max(binned_v_data[x,y][pressure_idx])
                                binned_v_min[x,y,i] = np.min(binned_v_data[x,y][pressure_idx])
                                binned_v_std[x,y,i] = np.std(binned_v_data[x,y][pressure_idx])
                                binned_v_rmse[x,y,i] = np.sqrt(np.nanmean(np.square(binned_v_data[x,y][pressure_idx])))
                                
        binnedData = {'u': {'binned_nobs': binned_u_nobs,
                            'binned_mean': binned_u_mean,
                            'binned_max':  binned_u_max,
                            'binned_min':  binned_u_min,
                            'binned_std':  binned_u_std,
                            'binned_rmse': binned_u_rmse
                           },
                      'v': {'binned_nobs': binned_v_nobs,
                            'binned_mean': binned_v_mean,
                            'binned_max':  binned_v_max,
                            'binned_min':  binned_v_min,
                            'binned_std':  binned_v_std,
                            'binned_rmse': binned_v_rmse
                           }
                     }
        
        return binnedData
    
    else:
        binned_data, binned_pressure = binData(data, lat, lon,
                                               pressure,
                                               binsize=binsize, 
                                               uvData=False)
        rows = binned_data.shape[0]
        cols = binned_data.shape[1]

        binned_nobs = np.full((rows, cols, 9), np.nan)
        binned_mean = np.full((rows, cols, 9), np.nan)
        binned_max = np.full((rows, cols, 9), np.nan)
        binned_min = np.full((rows, cols, 9), np.nan)
        binned_std = np.full((rows, cols, 9), np.nan)
        binned_rmse = np.full((rows, cols, 9), np.nan)

        pressure_list = [None, 0, 100, 250, 500, 700, 850, 925, 1000, 1100]

        for i, pressure in enumerate(pressure_list[:-1]):
            for x in range(0, rows):
                for y in range(0, cols):
                    if np.isnan(binned_data[x,y]).any() == False:
                        if i == 0:
                            binned_nobs[x,y,i] = len(binned_data[x,y])
                            binned_mean[x,y,i] = np.mean(binned_data[x,y])
                            binned_max[x,y,i] = np.max(binned_data[x,y])
                            binned_min[x,y,i] = np.min(binned_data[x,y])
                            binned_std[x,y,i] = np.std(binned_data[x,y])
                            binned_rmse[x,y,i] = np.sqrt(np.nanmean(np.square(binned_data[x,y])))
                        else:
                            pressure_idx = np.where((binned_pressure[x,y] > pressure_list[i]) & (binned_pressure[x,y] < pressure_list[i+1]))
                            if len(pressure_idx[0]) > 0:
                                binned_nobs[x,y,i] = len(binned_data[x,y][pressure_idx])
                                binned_mean[x,y,i] = np.mean(binned_data[x,y][pressure_idx])
                                binned_max[x,y,i] = np.max(binned_data[x,y][pressure_idx])
                                binned_min[x,y,i] = np.min(binned_data[x,y][pressure_idx])
                                binned_std[x,y,i] = np.std(binned_data[x,y][pressure_idx])
                                binned_rmse[x,y,i] = np.sqrt(np.nanmean(np.square(binned_data[x,y][pressure_idx])))
    
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
            binnedData = spatialBin(data, lat, lon, pressure, binsize='1x1', uvData=True)
        else:
            binnedData = spatialBin(data, lat, lon, pressure, binsize='1x1')


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
            binnedData = spatialBin(data, lat, lon, pressure, binsize='1x1', uvData=True)
        else:
            binnedData = spatialBin(data, lat, lon, pressure, binsize='1x1')

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