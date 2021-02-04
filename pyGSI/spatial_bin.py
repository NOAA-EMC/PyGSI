#!/usr/bin/env python

import numpy as np
import itertools


def bin_data(data, lat, lon, binsize=1, uv_data=False, pressure=None):
    """
    The main function to spatially bin the data.
    Inputs:
        data      : data to binned
        lat       : original data lats
        lon       : original data lons
        binsize   : string of the size of binning, must be square.
                    Examples: '1x1', '2x2', '5x5' (Default = 1x1)
        uv_data    : if using uv_data, will be True (Default = False)
        pressure  : original data pressure values (Default = None)
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
    if uv_data == True:
        binned_u_data = np.full((latbin, lonbin), np.nan, dtype=object)
        binned_v_data = np.full((latbin, lonbin), np.nan, dtype=object)

        if pressure is not None:
            binned_pressure = np.full((latbin, lonbin), np.nan, dtype=object)

        for val in paramlist:
            # Get index of 1x1 grid lat and lon
            latidx = np.where(ll_lats == val[0])
            lonidx = np.where(ll_lons == val[1])
            # values of the 1x1 grid lat and lon
            binnedlons = val[1]
            binnedlats = val[0]

            # find instances where data is within 1x1 grid point of orginal data
            data_idx = np.where((lon >= binnedlons - n_deg) & (lon <= binnedlons + n_deg) &
                                (lat >= binnedlats - n_deg) & (lat <= binnedlats + n_deg))

            latlon_idx = [latidx[0][0], lonidx[0][0]]

            # calculate stats if there is data at this grid point, else append np.nan
            if len(data_idx[0]) > 0:
                u = data['u'][data_idx]
                v = data['v'][data_idx]

                binned_u_data[latlon_idx[0], latlon_idx[1]] = u
                binned_v_data[latlon_idx[0], latlon_idx[1]] = v

                if pressure is not None:
                    p = pressure[data_idx]
                    binned_pressure[latlon_idx[0], latlon_idx[1]] = p

        if pressure is not None:
            return binned_u_data, binned_v_data, binned_pressure

        else:
            return binned_u_data, binned_v_data

    else:
        binned_data = np.full((latbin, lonbin), np.nan, dtype=object)
        if pressure is not None:
            binned_pressure = np.full((latbin, lonbin), np.nan, dtype=object)

        for val in paramlist:
            # Get index of grid lat and lon
            latidx = np.where(ll_lats == val[0])
            lonidx = np.where(ll_lons == val[1])
            # values of the 1x1 grid lat and lon
            binnedlons = val[1]
            binnedlats = val[0]

            # find instances where data is within 1x1 grid point of orginal data
            data_idx = np.where((lon >= binnedlons - n_deg) & (lon <= binnedlons + n_deg) &
                                (lat >= binnedlats - n_deg) & (lat <= binnedlats + n_deg))

            latlon_idx = [latidx[0][0], lonidx[0][0]]

            # calculate stats if there is data at this grid point
            if len(data_idx[0]) > 0:
                d = data[data_idx]
                binned_data[latlon_idx[0], latlon_idx[1]] = d

                if pressure is not None:
                    p = pressure[data_idx]
                    binned_pressure[latlon_idx[0], latlon_idx[1]] = p

        if pressure is not None:
            return binned_data, binned_pressure

        else:
            return binned_data


def spatial_bin(data, metadata, lat, lon, binsize=1, pressure=None, pbins=None):
    """
    Function to spatially bin data. Has the option to bin by pressures as well.
    Inputs:
        data      : data to binned
        lat       : original data lats
        lon       : original data lons
        binsize   : string of the size of binning, must be square.
                    Examples: '1x1', '2x2', '5x5' (Default = 1x1)
        ** ONLY APPLICABLE FOR CONVENTIONAL DATA **
        uv_data    : if using uv_data, will be True (Default = False)
        pressure  : original data pressure values (Default = None)
        pbins     : list of pressures in ascending order of what pressure
                    levels data will be binned into. Must provide data for
                    pressure to be used. (Default = None)
    Outputs:
        binned_data : a dictionary of statistics calculated for the newly
                      binned data. If binning by pressure for conventional
                      data, will return 3D data with the first index being
                      all pressure levels, followed by n-1 amount of indexes
                      based on the number of pressure levels given in pbins
    """
    
    
    uv_data = True if 'Variable' in metadata and metadata['Variable'] == 'uv' else False 

    if pbins is None:
        pressure_list = [None, 0, 100, 250, 500, 700, 850, 925, 1000, 1100]
    else:
        pbins.insert(0, None)
        pressure_list = pbins
        
    if uv_data:
        if pressure.any() == None:
            binned_u_data, binned_v_data = bin_data(data, lat, lon,
                                                   binsize=binsize,
                                                   uv_data=uv_data,
                                                   pressure=pressure)

            rows = binned_u_data.shape[0]
            cols = binned_u_data.shape[1]

            binned_u_nobs = np.full((rows, cols), np.nan)
            binned_u_mean = np.full((rows, cols), np.nan)
            binned_u_max = np.full((rows, cols), np.nan)
            binned_u_min = np.full((rows, cols), np.nan)
            binned_u_std = np.full((rows, cols), np.nan)
            binned_u_rmse = np.full((rows, cols), np.nan)

            binned_v_nobs = np.full((rows, cols), np.nan)
            binned_v_mean = np.full((rows, cols), np.nan)
            binned_v_max = np.full((rows, cols), np.nan)
            binned_v_min = np.full((rows, cols), np.nan)
            binned_v_std = np.full((rows, cols), np.nan)
            binned_v_rmse = np.full((rows, cols), np.nan)

            for x in range(0, rows):
                for y in range(0, cols):
                    if np.isnan(binned_u_data[x, y]).any() == False:
                        binned_u_nobs[x, y] = len(binned_u_data[x, y])
                        binned_u_mean[x, y] = np.mean(binned_u_data[x, y])
                        binned_u_max[x, y] = np.max(binned_u_data[x, y])
                        binned_u_min[x, y] = np.min(binned_u_data[x, y])
                        binned_u_std[x, y] = np.std(binned_u_data[x, y])
                        binned_u_rmse[x, y] = np.sqrt(
                            np.nanmean(np.square(binned_u_data[x, y])))

                        binned_v_nobs[x, y] = len(binned_v_data[x, y])
                        binned_v_mean[x, y] = np.mean(binned_v_data[x, y])
                        binned_v_max[x, y] = np.max(binned_v_data[x, y])
                        binned_v_min[x, y] = np.min(binned_v_data[x, y])
                        binned_v_std[x, y] = np.std(binned_v_data[x, y])
                        binned_v_rmse[x, y] = np.sqrt(
                            np.nanmean(np.square(binned_v_data[x, y])))

        else:
            binned_u_data, binned_v_data, binned_pressure = bin_data(data, lat, lon,
                                                                    binsize=binsize,
                                                                    uv_data=uv_data,
                                                                    pressure=pressure)

            rows = binned_u_data.shape[0]
            cols = binned_u_data.shape[1]
            
            n_plevs = len(pressure_list)-1

            binned_u_nobs = np.full((rows, cols, n_plevs), np.nan)
            binned_u_mean = np.full((rows, cols, n_plevs), np.nan)
            binned_u_max = np.full((rows, cols, n_plevs), np.nan)
            binned_u_min = np.full((rows, cols, n_plevs), np.nan)
            binned_u_std = np.full((rows, cols, n_plevs), np.nan)
            binned_u_rmse = np.full((rows, cols, n_plevs), np.nan)

            binned_v_nobs = np.full((rows, cols, n_plevs), np.nan)
            binned_v_mean = np.full((rows, cols, n_plevs), np.nan)
            binned_v_max = np.full((rows, cols, n_plevs), np.nan)
            binned_v_min = np.full((rows, cols, n_plevs), np.nan)
            binned_v_std = np.full((rows, cols, n_plevs), np.nan)
            binned_v_rmse = np.full((rows, cols, n_plevs), np.nan)

            for i, pressure in enumerate(pressure_list[:-1]):
                for x in range(0, rows):
                    for y in range(0, cols):
                        if np.isnan(binned_u_data[x, y]).any() == False:
                            if i == 0:
                                binned_u_nobs[x, y, i] = len(
                                    binned_u_data[x, y])
                                binned_u_mean[x, y, i] = np.mean(
                                    binned_u_data[x, y])
                                binned_u_max[x, y, i] = np.max(
                                    binned_u_data[x, y])
                                binned_u_min[x, y, i] = np.min(
                                    binned_u_data[x, y])
                                binned_u_std[x, y, i] = np.std(
                                    binned_u_data[x, y])
                                binned_u_rmse[x, y, i] = np.sqrt(
                                    np.nanmean(np.square(binned_u_data[x, y])))

                                binned_v_nobs[x, y, i] = len(
                                    binned_v_data[x, y])
                                binned_v_mean[x, y, i] = np.mean(
                                    binned_v_data[x, y])
                                binned_v_max[x, y, i] = np.max(
                                    binned_v_data[x, y])
                                binned_v_min[x, y, i] = np.min(
                                    binned_v_data[x, y])
                                binned_v_std[x, y, i] = np.std(
                                    binned_v_data[x, y])
                                binned_v_rmse[x, y, i] = np.sqrt(
                                    np.nanmean(np.square(binned_v_data[x, y])))
                            else:
                                pressure_idx = np.where((binned_pressure[x, y] > pressure_list[i]) & (
                                    binned_pressure[x, y] < pressure_list[i+1]))
                                if len(pressure_idx[0]) > 0:
                                    binned_u_nobs[x, y, i] = len(
                                        binned_u_data[x, y][pressure_idx])
                                    binned_u_mean[x, y, i] = np.mean(
                                        binned_u_data[x, y][pressure_idx])
                                    binned_u_max[x, y, i] = np.max(
                                        binned_u_data[x, y][pressure_idx])
                                    binned_u_min[x, y, i] = np.min(
                                        binned_u_data[x, y][pressure_idx])
                                    binned_u_std[x, y, i] = np.std(
                                        binned_u_data[x, y][pressure_idx])
                                    binned_u_rmse[x, y, i] = np.sqrt(np.nanmean(
                                        np.square(binned_u_data[x, y][pressure_idx])))

                                    binned_v_nobs[x, y, i] = len(
                                        binned_v_data[x, y][pressure_idx])
                                    binned_v_mean[x, y, i] = np.mean(
                                        binned_v_data[x, y][pressure_idx])
                                    binned_v_max[x, y, i] = np.max(
                                        binned_v_data[x, y][pressure_idx])
                                    binned_v_min[x, y, i] = np.min(
                                        binned_v_data[x, y][pressure_idx])
                                    binned_v_std[x, y, i] = np.std(
                                        binned_v_data[x, y][pressure_idx])
                                    binned_v_rmse[x, y, i] = np.sqrt(np.nanmean(
                                        np.square(binned_v_data[x, y][pressure_idx])))

        binned_data = {'u': {'binned_nobs': binned_u_nobs,
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

        return binned_data

    else:
        if pressure is None:
            binned_data = bin_data(data, lat, lon, binsize=binsize)

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
                    if np.isnan(binned_data[x, y]).any() == False:
                        binned_nobs[x, y] = len(binned_data[x, y])
                        binned_mean[x, y] = np.mean(binned_data[x, y])
                        binned_max[x, y] = np.max(binned_data[x, y])
                        binned_min[x, y] = np.min(binned_data[x, y])
                        binned_std[x, y] = np.std(binned_data[x, y])
                        binned_rmse[x, y] = np.sqrt(
                            np.nanmean(np.square(binned_data[x, y])))

        else:
            binned_data, binned_pressure = bin_data(data, lat, lon,
                                                   binsize=binsize,
                                                   uv_data=False,
                                                   pressure=pressure)
            rows = binned_data.shape[0]
            cols = binned_data.shape[1]
            
            n_plevs = len(pressure_list)-1

            binned_nobs = np.full((rows, cols, n_plevs), np.nan)
            binned_mean = np.full((rows, cols, n_plevs), np.nan)
            binned_max = np.full((rows, cols, n_plevs), np.nan)
            binned_min = np.full((rows, cols, n_plevs), np.nan)
            binned_std = np.full((rows, cols, n_plevs), np.nan)
            binned_rmse = np.full((rows, cols, n_plevs), np.nan)

            for i, pressure in enumerate(pressure_list[:-1]):
                for x in range(0, rows):
                    for y in range(0, cols):
                        if np.isnan(binned_data[x, y]).any() == False:
                            if i == 0:
                                binned_nobs[x, y, i] = len(binned_data[x, y])
                                binned_mean[x, y, i] = np.mean(
                                    binned_data[x, y])
                                binned_max[x, y, i] = np.max(binned_data[x, y])
                                binned_min[x, y, i] = np.min(binned_data[x, y])
                                binned_std[x, y, i] = np.std(binned_data[x, y])
                                binned_rmse[x, y, i] = np.sqrt(
                                    np.nanmean(np.square(binned_data[x, y])))
                            else:
                                pressure_idx = np.where((binned_pressure[x, y] > pressure_list[i]) &
                                                        (binned_pressure[x, y] < pressure_list[i+1]))
                                if len(pressure_idx[0]) > 0:
                                    binned_nobs[x, y, i] = len(
                                        binned_data[x, y][pressure_idx])
                                    binned_mean[x, y, i] = np.mean(
                                        binned_data[x, y][pressure_idx])
                                    binned_max[x, y, i] = np.max(
                                        binned_data[x, y][pressure_idx])
                                    binned_min[x, y, i] = np.min(
                                        binned_data[x, y][pressure_idx])
                                    binned_std[x, y, i] = np.std(
                                        binned_data[x, y][pressure_idx])
                                    binned_rmse[x, y, i] = np.sqrt(np.nanmean(
                                        np.square(binned_data[x, y][pressure_idx])))

        binned_data = {'binned_nobs': binned_nobs,
                      'binned_mean': binned_mean,
                      'binned_max': binned_max,
                      'binned_min': binned_min,
                      'binned_std': binned_std,
                      'binned_rmse': binned_rmse
                      }

        return binned_data
