#!/usr/bin/env python
## Created by Kevin Dougherty
## October 2020

import numpy as np
import os.path
from datetime import datetime
from netCDF4 import Dataset

def calculate_stats(data):
    
    n = len(data)
    
    if n == 0:
    
        statdict = {'n'    : 0,
                    'min'  : np.nan,
                    'max'  : np.nan,
                    'mean' : np.nan,
                    'std'  : np.nan,
                    'q25'  : np.nan,
                    'q50'  : np.nan,
                    'q75'  : np.nan,
                    'rmse' : np.nan
                   }

        return statdict
    
    else:
    
        mean = np.nanmean(data)
        std = np.nanstd(data)
        mx = np.nanmax(data)
        mn = np.nanmin(data)
        q25 = np.nanquantile(data, 0.25)
        q50 = np.nanquantile(data, 0.50)
        q75 = np.nanquantile(data, 0.75)
        rmse = np.sqrt(np.nanmean(np.square(data)))

        statdict = {'n'    : n,
                    'min'  : mn,
                    'max'  : mx,
                    'mean' : mean,
                    'std'  : std,
                    'q25'  : q25,
                    'q50'  : q50,
                    'q75'  : q75,
                    'rmse' : rmse
                   }

        return statdict
    
def write_ncfile(statDict, ncfilename):
    
    ncfile = Dataset(ncfilename, 'w', format='NETCDF4')

    time_dim = ncfile.createDimension("time", None)
    lat_dim = ncfile.createDimension("lats", 180)
    lons_dim = ncfile.createDimension("lons", 360)

    ncfile.diag_type='Conventional'
    ncfile.variable='t'
    ncfile.observationID='120'
    # Latitude
    lat = ncfile.createVariable('lat', np.float32, ('lats',))
    lat.units='degrees'
    lat.long_name='latitude'
    # Longitude
    lon = ncfile.createVariable('lon', np.float32, ('lons',))
    lon.units='degrees'
    lon.long_name='longitude'
    
    subtype = ncfile.createVariable('subtype', np.int64, ('time'))#('subtype'))
    subtype.units = 'subtype'
    subtype.long_name = 'subtype of specified variable'

    date = ncfile.createVariable('date', np.int64, ('time'))#, 'subtype'))
    date.units = 'time'
    date.long_name = 'integer time represented as YYYYMMDDHH'

    epoch = ncfile.createVariable('epoch', np.int64, ('time'))#, 'subtype'))
    epoch.units = 'seconds since 01-01-1970 00:00:00'
    epoch.long_name = 'epoch time'

    obscount = ncfile.createVariable('obscount', np.int32, ('time'))#, 'subtype'))
    obscount.standard_name = 'observation count'

    mean = ncfile.createVariable('mean', np.float32, ('time'))#, 'subtype'))
    mean.standard_name = 'mean'

    mx = ncfile.createVariable('max', np.float32, ('time'))#, 'subtype'))
    mx.standard_name = 'maximum'

    mn = ncfile.createVariable('min', np.float32, ('time'))#, 'subtype'))
    mn.standard_name = 'minimum'

    std = ncfile.createVariable('std', np.float32, ('time'))#, 'subtype'))
    std.standard_name = 'standard deviation'

    rmse = ncfile.createVariable('rmse', np.float32, ('time'))#, 'subtype'))
    rmse.standard_name = 'root mean square error'

    q25 = ncfile.createVariable('q25', np.float32, ('time'))#, 'subtype'))
    q25.standard_name = '25th quantile'

    q50 = ncfile.createVariable('q50', np.float32, ('time'))#, 'subtype'))
    q50.standard_name = '50th quantile'

    q75 = ncfile.createVariable('q75', np.float32, ('time'))#, 'subtype'))
    q75.standard_name = '75th quantile'
    
    ######################################################################
    
    obscount_ll = ncfile.createVariable('obscount', np.int32, ('time', 'lon', 'lat'))
    obscount_ll.standard_name = 'observation count'

    mean_ll = ncfile.createVariable('mean', np.float32, ('time', 'lon', 'lat'))
    mean_ll.standard_name = 'mean'

    mx_ll = ncfile.createVariable('max', np.float32, ('time', 'lon', 'lat'))
    mx_ll.standard_name = 'maximum'

    mn_ll = ncfile.createVariable('min', np.float32, ('time', 'lon', 'lat'))
    mn_ll.standard_name = 'minimum'

    std_ll = ncfile.createVariable('std', np.float32, ('time', 'lon', 'lat'))
    std_ll.standard_name = 'standard deviation'

    rmse_ll = ncfile.createVariable('rmse', np.float32, ('time', 'lon', 'lat'))
    rmse_ll.standard_name = 'root mean square error'

    q25_ll = ncfile.createVariable('q25', np.float32, ('time', 'lon', 'lat'))
    q25_ll.standard_name = '25th quantile'

    q50_ll = ncfile.createVariable('q50', np.float32, ('time', 'lon', 'lat'))
    q50_ll.standard_name = '50th quantile'

    q75_ll = ncfile.createVariable('q75', np.float32, ('time', 'lon', 'lat'))
    q75_ll.standard_name = '75th quantile'
    
    ######################################################################

    lat[:] = statDict['lats']
    lon[:] = statDict['lons']
    
    subtype[:] = statDict['Subtype']
    date[:] = statDict['Date']
    epoch[:] = statDict['epoch']
    obscount[:] = statDict['n']
    mean[:] = statDict['mean']
    mx[:] = statDict['max']
    mn[:] = statDict['min']
    std[:] = statDict['std']
    rmse[:] = statDict['rmse']
    q25[:] = statDict['q25']
    q50[:] = statDict['q50']
    q75[:] = statDict['q75']
    
    ncfile.close()
    
    return

def append_ncfile(statDict, ncfilename):
    
    ncfile = Dataset(ncfilename, 'a', format='NETCDF4')
    
    subtype= ncfile.variables['subtype']
    date = ncfile.variables['date']
    epoch = ncfile.variables['epoch']
    obscount = ncfile.variables['obscount']
    mean = ncfile.variables['mean']
    mx = ncfile.variables['max']
    mn = ncfile.variables['min']
    std = ncfile.variables['std']
    rmse = ncfile.variables['rmse']
    q25 = ncfile.variables['q25']
    q50 = ncfile.variables['q50']
    q75 = ncfile.variables['q75']
    
    idx = len(date)
    
    lat[idx] = statDict['lats']
    lon[idx] = statDict['lons']
    
    subtype[idx] = statDict['Subtype']
    date[idx] = statDict['Date']
    epoch[idx] = statDict['epoch']
    obscount[idx] = statDict['n']
    mean[idx] = statDict['mean']
    mx[idx] = statDict['max']
    mn[idx] = statDict['min']
    std[idx] = statDict['std']
    rmse[idx] = statDict['rmse']
    q25[idx] = statDict['q25']
    q50[idx] = statDict['q50']
    q75[idx] = statDict['q75']
    
    ncfile.close()
    
    return


def get_filename(metadata):
    
    if metadata['Data_type'] == 'O-F':
        if metadata['Diag_type'] == 'conv':
            ncfilename = '{Diag_type}_{Variable}_{ObsID[0]}_OmF.nc'.format(**metadata)
        else:
            ncfilename = '{Diag_type}_{Satellite}_Ch{Channel[0]}_OmF.nc'.format(**metadata)
            
    elif metadata['Data_type'] == 'O-A':
        if metadata['Diag_type'] == 'conv':
            ncfilename = '{Diag_type}_{Variable}_{ObsID[0]}_OmA.nc'.format(**metadata)
        else:
            ncfilename = '{Diag_type}_{Satellite}_Ch{Channel[0]}_OmA.nc'.format(**metadata)
            
    elif metadata['Data_type'] == 'H(x)':
        if metadata['Diag_type'] == 'conv':
            ncfilename = '{Diag_type}_{Variable}_{ObsID[0]}_HofX.nc'.format(**metadata)
        else:
            ncfilename = '{Diag_type}_{Satellite}_Ch{Channel[0]}_HofX.nc'.format(**metadata)
            
    else:
        if metadata['Diag_type'] == 'conv':
            ncfilename = '{Diag_type}_{Variable}_{ObsID[0]}_{Diag_type}.nc'.format(**metadata)
        else:
            ncfilename = '{Diag_type}_{Satellite}_Ch{Channel[0]}_{Diag_type}.nc'.format(**metadata)
            
    return ncfilename


def writeNetCDF(data, metadata, lat, lon):
    
    if metadata['Diag_type'] == 'conv' and metadata['Variable'] == 'uv':
            
        for i in data.keys():
            metadata['Variable'] = i
            statDict = calculate_stats(data[i])

            statDict['Date'] = metadata['Date'].strftime("%Y%m%d%H")
            statDict['epoch'] = metadata['Date'].timestamp()
            statDict['Subtype'] = int(metadata['Subtype'][0])
            
            lats = np.linspace(-89.5,89.5,180)
            lons = np.linspace(0.5,359.5,360)
        
            statDict['lats'] = lats
            statDict['lons'] = lons

            # Get the filename based on metadata
            ncfilename = get_filename(metadata)

            ncpath = metadata['outDir']+'/'+ncfilename

            # See if the file exists. If it doesn't, create it
            if os.path.isfile(ncpath) != True:
                write_ncfile(statDict, ncpath)

            # else, append data to that file
            else:
                append_ncfile(statDict, ncpath)

    else:
    
        statDict = calculate_stats(data)

        statDict['Date'] = int(metadata['Date'].strftime("%Y%m%d%H"))
        statDict['epoch'] = metadata['Date'].timestamp()
        statDict['Subtype'] = int(metadata['Subtype'][0])
        
        lats = np.linspace(-89.5,89.5,180)
        lons = np.linspace(0.5,359.5,360)
        
        statDict['lats'] = lats
        statDict['lons'] = lons

        # Get the filename based on metadata
        ncfilename = get_filename(metadata)

        ncpath = metadata['outDir']+'/'+ncfilename

        # See if the file exists. If it doesn't, create it
        if os.path.isfile(ncpath) != True:
            write_ncfile(statDict, ncpath)

        # else, append data to that file
        else:
            append_ncfile(statDict, ncpath)
        
    return