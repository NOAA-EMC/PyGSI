#!/usr/bin/env python
# Created by Kevin Dougherty
# October 2020

import numpy as np
import os.path
from datetime import datetime
from netCDF4 import Dataset


def _calculate_stats(data):

    n = len(data)

    if n == 0:

        stat_dict = {'n': 0,
                     'min': np.nan,
                     'max': np.nan,
                     'mean': np.nan,
                     'std': np.nan,
                     'q25': np.nan,
                     'q50': np.nan,
                     'q75': np.nan,
                     'rmse': np.nan
                     }

        return stat_dict

    else:

        mean = np.nanmean(data)
        std = np.nanstd(data)
        mx = np.nanmax(data)
        mn = np.nanmin(data)
        q25 = np.nanquantile(data, 0.25)
        q50 = np.nanquantile(data, 0.50)
        q75 = np.nanquantile(data, 0.75)
        rmse = np.sqrt(np.nanmean(np.square(data)))

        stat_dict = {'n': n,
                     'min': mn,
                     'max': mx,
                     'mean': mean,
                     'std': std,
                     'q25': q25,
                     'q50': q50,
                     'q75': q75,
                     'rmse': rmse
                     }

        return stat_dict


def _write_ncfile(stat_dict, binned_data, ncfilename, obs_type):

    ncfile = Dataset(ncfilename, 'w', format='NETCDF4')

    time_dim = ncfile.createDimension("time", None)
    lat_dim = ncfile.createDimension("lats", 180)
    lons_dim = ncfile.createDimension("lons", 360)

    date = ncfile.createVariable('date', np.int64, ('time'))
    date.units = 'time'
    date.long_name = 'integer time represented as YYYYMMDDHH'

    epoch = ncfile.createVariable('epoch', np.int64, ('time'))
    epoch.units = 'seconds since 01-01-1970 00:00:00'
    epoch.long_name = 'epoch time'

    obscount = ncfile.createVariable('obscount', np.int32, ('time'))
    obscount.standard_name = 'observation count'

    mean = ncfile.createVariable('mean', np.float32, ('time'))
    mean.standard_name = 'mean'

    mx = ncfile.createVariable('max', np.float32, ('time'))
    mx.standard_name = 'maximum'

    mn = ncfile.createVariable('min', np.float32, ('time'))
    mn.standard_name = 'minimum'

    std = ncfile.createVariable('std', np.float32, ('time'))
    std.standard_name = 'standard deviation'

    rmse = ncfile.createVariable('rmse', np.float32, ('time'))
    rmse.standard_name = 'root mean square error'

    q25 = ncfile.createVariable('q25', np.float32, ('time'))
    q25.standard_name = '25th quantile'

    q50 = ncfile.createVariable('q50', np.float32, ('time'))
    q50.standard_name = '50th quantile'

    q75 = ncfile.createVariable('q75', np.float32, ('time'))
    q75.standard_name = '75th quantile'

    date[:] = stat_dict['Date']
    epoch[:] = stat_dict['epoch']
    obscount[:] = stat_dict['n']
    mean[:] = stat_dict['mean']
    mx[:] = stat_dict['max']
    mn[:] = stat_dict['min']
    std[:] = stat_dict['std']
    rmse[:] = stat_dict['rmse']
    q25[:] = stat_dict['q25']
    q50[:] = stat_dict['q50']
    q75[:] = stat_dict['q75']

    ######################################################################

    if obs_type == 'conv':

        pressure_dim = ncfile.createDimension("pressure", 9)

        subtype = ncfile.createVariable('subtype', np.int64, ('time'))
        subtype.units = 'subtype'
        subtype.long_name = 'subtype of specified variable'

        binned_obscount = ncfile.createVariable(
            'binned_obscount', np.int32, ('time', 'lats', 'lons', 'pressure'))
        binned_obscount.standard_name = 'binned observation count'

        binned_mean = ncfile.createVariable(
            'binned_mean', np.float32, ('time', 'lats', 'lons', 'pressure'))
        binned_mean.standard_name = 'binned mean'

        binned_mx = ncfile.createVariable(
            'binned_max', np.float32, ('time', 'lats', 'lons', 'pressure'))
        binned_mx.standard_name = 'binned maximum'

        binned_mn = ncfile.createVariable(
            'binned_min', np.float32, ('time', 'lats', 'lons', 'pressure'))
        binned_mn.standard_name = 'binned minimum'

        binned_std = ncfile.createVariable(
            'binned_std', np.float32, ('time', 'lats', 'lons', 'pressure'))
        binned_std.standard_name = 'binned standard deviation'

        binned_rmse = ncfile.createVariable(
            'binned_rmse', np.float32, ('time', 'lats', 'lons', 'pressure'))
        binned_rmse.standard_name = 'binned root mean square error'

        subtype[:] = stat_dict['Subtype']

        binned_obscount[0, :, :, :] = binned_data['binned_nobs']
        binned_mean[0, :, :, :] = binned_data['binned_mean']
        binned_mx[0, :, :, :] = binned_data['binned_max']
        binned_mn[0, :, :, :] = binned_data['binned_min']
        binned_std[0, :, :, :] = binned_data['binned_std']
        binned_rmse[0, :, :, :] = binned_data['binned_rmse']

    ######################################################################

    else:

        binned_obscount = ncfile.createVariable(
            'binned_obscount', np.int32, ('time', 'lats', 'lons'))
        binned_obscount.standard_name = 'binned observation count'

        binned_mean = ncfile.createVariable(
            'binned_mean', np.float32, ('time', 'lats', 'lons'))
        binned_mean.standard_name = 'binned mean'

        binned_mx = ncfile.createVariable(
            'binned_max', np.float32, ('time', 'lats', 'lons'))
        binned_mx.standard_name = 'binned maximum'

        binned_mn = ncfile.createVariable(
            'binned_min', np.float32, ('time', 'lats', 'lons'))
        binned_mn.standard_name = 'binned minimum'

        binned_std = ncfile.createVariable(
            'binned_std', np.float32, ('time', 'lats', 'lons'))
        binned_std.standard_name = 'binned standard deviation'

        binned_rmse = ncfile.createVariable(
            'binned_rmse', np.float32, ('time', 'lats', 'lons'))
        binned_rmse.standard_name = 'binned root mean square error'

        binned_obscount[0, :, :] = binned_data['binned_nobs']
        binned_mean[0, :, :] = binned_data['binned_mean']
        binned_mx[0, :, :] = binned_data['binned_max']
        binned_mn[0, :, :] = binned_data['binned_min']
        binned_std[0, :, :] = binned_data['binned_std']
        binned_rmse[0, :, :] = binned_data['binned_rmse']

    ncfile.close()

    return


def _append_ncfile(stat_dict, binned_data, ncfilename, obs_type):

    ncfile = Dataset(ncfilename, 'a', format='NETCDF4')

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

    binned_obscount = ncfile.variables['binned_obscount']
    binned_mean = ncfile.variables['binned_mean']
    binned_max = ncfile.variables['binned_max']
    binned_min = ncfile.variables['binned_min']
    binned_std = ncfile.variables['binned_std']
    binned_rmse = ncfile.variables['binned_rmse']

    idx = len(date)

    date[idx] = stat_dict['Date']
    epoch[idx] = stat_dict['epoch']
    obscount[idx] = stat_dict['n']
    mean[idx] = stat_dict['mean']
    mx[idx] = stat_dict['max']
    mn[idx] = stat_dict['min']
    std[idx] = stat_dict['std']
    rmse[idx] = stat_dict['rmse']
    q25[idx] = stat_dict['q25']
    q50[idx] = stat_dict['q50']
    q75[idx] = stat_dict['q75']

    if obs_type == 'conv':

        subtype = ncfile.variables['subtype']
        subtype[idx] = stat_dict['Subtype']

        binned_obscount[idx, :, :, :] = binned_data['binned_nobs']
        binned_mean[idx, :, :, :] = binned_data['binned_mean']
        binned_max[idx, :, :, :] = binned_data['binned_max']
        binned_min[idx, :, :, :] = binned_data['binned_min']
        binned_std[idx, :, :, :] = binned_data['binned_std']
        binned_rmse[idx, :, :, :] = binned_data['binned_rmse']

    else:

        binned_obscount[idx, :, :] = binned_data['binned_nobs']
        binned_mean[idx, :, :] = binned_data['binned_mean']
        binned_max[idx, :, :] = binned_data['binned_max']
        binned_min[idx, :, :] = binned_data['binned_min']
        binned_std[idx, :, :] = binned_data['binned_std']
        binned_rmse[idx, :, :] = binned_data['binned_rmse']

    ncfile.close()

    return


def _write_ncfile_uv(uv_stat_dict, binned_data, ncfilename):

    ncfile = Dataset(ncfilename, 'w', format='NETCDF4')

    time_dim = ncfile.createDimension("time", None)
    lat_dim = ncfile.createDimension("lats", 180)
    lons_dim = ncfile.createDimension("lons", 360)
    pressure_dim = ncfile.createDimension("pressure", 9)

    subtype = ncfile.createVariable('subtype', np.int64, ('time'))
    subtype.units = 'subtype'
    subtype.long_name = 'subtype of specified variable'

    date = ncfile.createVariable('date', np.int64, ('time'))
    date.units = 'time'
    date.long_name = 'integer time represented as YYYYMMDDHH'

    epoch = ncfile.createVariable('epoch', np.int64, ('time'))
    epoch.units = 'seconds since 01-01-1970 00:00:00'
    epoch.long_name = 'epoch time'

    u_obscount = ncfile.createVariable('obscount_u', np.int32, ('time'))
    u_obscount.standard_name = 'observation count for u wind'

    u_mean = ncfile.createVariable('mean_u', np.float32, ('time'))
    u_mean.standard_name = 'mean for u wind'

    u_mx = ncfile.createVariable('max_u', np.float32, ('time'))
    u_mx.standard_name = 'maximum for u wind'

    u_mn = ncfile.createVariable('min_u', np.float32, ('time'))
    u_mn.standard_name = 'minimum for u wind'

    u_std = ncfile.createVariable('std_u', np.float32, ('time'))
    u_std.standard_name = 'standard deviation for u wind'

    u_rmse = ncfile.createVariable('rmse_u', np.float32, ('time'))
    u_rmse.standard_name = 'root mean square error for u wind'

    u_q25 = ncfile.createVariable('q25_u', np.float32, ('time'))
    u_q25.standard_name = '25th quantile for u wind'

    u_q50 = ncfile.createVariable('q50_u', np.float32, ('time'))
    u_q50.standard_name = '50th quantile for u wind'

    u_q75 = ncfile.createVariable('q75_u', np.float32, ('time'))
    u_q75.standard_name = '75th quantile for u wind'

    v_obscount = ncfile.createVariable('obscount_v', np.int32, ('time'))
    v_obscount.standard_name = 'observation count for v wind'

    v_mean = ncfile.createVariable('mean_v', np.float32, ('time'))
    v_mean.standard_name = 'mean for v wind'

    v_mx = ncfile.createVariable('max_v', np.float32, ('time'))
    v_mx.standard_name = 'maximum for v wind'

    v_mn = ncfile.createVariable('min_v', np.float32, ('time'))
    v_mn.standard_name = 'minimum for v wind'

    v_std = ncfile.createVariable('std_v', np.float32, ('time'))
    v_std.standard_name = 'standard deviation for v wind'

    v_rmse = ncfile.createVariable('rmse_v', np.float32, ('time'))
    v_rmse.standard_name = 'root mean square error for v wind'

    v_q25 = ncfile.createVariable('q25_v', np.float32, ('time'))
    v_q25.standard_name = '25th quantile for v wind'

    v_q50 = ncfile.createVariable('q50_v', np.float32, ('time'))
    v_q50.standard_name = '50th quantile for v wind'

    v_q75 = ncfile.createVariable('q75_v', np.float32, ('time'))
    v_q75.standard_name = '75th quantile for v wind'

    ######################################################################

    u_binned_obscount = ncfile.createVariable(
        'binned_obscount_u', np.int32, ('time', 'lats', 'lons', 'pressure'))
    u_binned_obscount.standard_name = 'binned observation count for u wind'

    u_binned_mean = ncfile.createVariable(
        'binned_mean_u', np.float32, ('time', 'lats', 'lons', 'pressure'))
    u_binned_mean.standard_name = 'binned mean for u wind'

    u_binned_mx = ncfile.createVariable(
        'binned_max_u', np.float32, ('time', 'lats', 'lons', 'pressure'))
    u_binned_mx.standard_name = 'binned maximum for u wind'

    u_binned_mn = ncfile.createVariable(
        'binned_min_u', np.float32, ('time', 'lats', 'lons', 'pressure'))
    u_binned_mn.standard_name = 'binned minimum for u wind'

    u_binned_std = ncfile.createVariable(
        'binned_std_u', np.float32, ('time', 'lats', 'lons', 'pressure'))
    u_binned_std.standard_name = 'binned standard deviation for u wind'

    u_binned_rmse = ncfile.createVariable(
        'binned_rmse_u', np.float32, ('time', 'lats', 'lons', 'pressure'))
    u_binned_rmse.standard_name = 'binned root mean square error for u wind'

    v_binned_obscount = ncfile.createVariable(
        'binned_obscount_v', np.int32, ('time', 'lats', 'lons', 'pressure'))
    v_binned_obscount.standard_name = 'binned observation count for v wind'

    v_binned_mean = ncfile.createVariable(
        'binned_mean_v', np.float32, ('time', 'lats', 'lons', 'pressure'))
    v_binned_mean.standard_name = 'binned mean for v wind'

    v_binned_mx = ncfile.createVariable(
        'binned_max_v', np.float32, ('time', 'lats', 'lons', 'pressure'))
    v_binned_mx.standard_name = 'binned maximum for u wind'

    v_binned_mn = ncfile.createVariable(
        'binned_min_v', np.float32, ('time', 'lats', 'lons', 'pressure'))
    v_binned_mn.standard_name = 'binned minimum for v wind'

    v_binned_std = ncfile.createVariable(
        'binned_std_v', np.float32, ('time', 'lats', 'lons', 'pressure'))
    v_binned_std.standard_name = 'binned standard deviation for v wind'

    v_binned_rmse = ncfile.createVariable(
        'binned_rmse_v', np.float32, ('time', 'lats', 'lons', 'pressure'))
    v_binned_rmse.standard_name = 'binned root mean square error for v wind'

    ######################################################################

    subtype[:] = uv_stat_dict['Subtype']
    date[:] = uv_stat_dict['Date']
    epoch[:] = uv_stat_dict['epoch']

    u_obscount[:] = uv_stat_dict['u']['n']
    u_mean[:] = uv_stat_dict['u']['mean']
    u_mx[:] = uv_stat_dict['u']['max']
    u_mn[:] = uv_stat_dict['u']['min']
    u_std[:] = uv_stat_dict['u']['std']
    u_rmse[:] = uv_stat_dict['u']['rmse']
    u_q25[:] = uv_stat_dict['u']['q25']
    u_q50[:] = uv_stat_dict['u']['q50']
    u_q75[:] = uv_stat_dict['u']['q75']

    v_obscount[:] = uv_stat_dict['v']['n']
    v_mean[:] = uv_stat_dict['v']['mean']
    v_mx[:] = uv_stat_dict['v']['max']
    v_mn[:] = uv_stat_dict['v']['min']
    v_std[:] = uv_stat_dict['v']['std']
    v_rmse[:] = uv_stat_dict['v']['rmse']
    v_q25[:] = uv_stat_dict['v']['q25']
    v_q50[:] = uv_stat_dict['v']['q50']
    v_q75[:] = uv_stat_dict['v']['q75']

    u_binned_obscount[0, :, :, :] = binned_data['u']['binned_nobs']
    u_binned_mean[0, :, :, :] = binned_data['u']['binned_mean']
    u_binned_mx[0, :, :, :] = binned_data['u']['binned_max']
    u_binned_mn[0, :, :, :] = binned_data['u']['binned_min']
    u_binned_std[0, :, :, :] = binned_data['u']['binned_std']
    u_binned_rmse[0, :, :, :] = binned_data['u']['binned_rmse']

    v_binned_obscount[0, :, :, :] = binned_data['v']['binned_nobs']
    v_binned_mean[0, :, :, :] = binned_data['v']['binned_mean']
    v_binned_mx[0, :, :, :] = binned_data['v']['binned_max']
    v_binned_mn[0, :, :, :] = binned_data['v']['binned_min']
    v_binned_std[0, :, :, :] = binned_data['v']['binned_std']
    v_binned_rmse[0, :, :, :] = binned_data['v']['binned_rmse']

    ncfile.close()

    return


def _append_ncfile_uv(uv_stat_dict, binned_data, ncfilename):

    ncfile = Dataset(ncfilename, 'a', format='NETCDF4')

    subtype = ncfile.variables['subtype']
    date = ncfile.variables['date']
    epoch = ncfile.variables['epoch']

    u_obscount = ncfile.variables['obscount_u']
    u_mean = ncfile.variables['mean_u']
    u_mx = ncfile.variables['max_u']
    u_mn = ncfile.variables['min_u']
    u_std = ncfile.variables['std_u']
    u_rmse = ncfile.variables['rmse_u']
    u_q25 = ncfile.variables['q25_u']
    u_q50 = ncfile.variables['q50_u']
    u_q75 = ncfile.variables['q75_u']

    v_obscount = ncfile.variables['obscount_v']
    v_mean = ncfile.variables['mean_v']
    v_mx = ncfile.variables['max_v']
    v_mn = ncfile.variables['min_v']
    v_std = ncfile.variables['std_v']
    v_rmse = ncfile.variables['rmse_v']
    v_q25 = ncfile.variables['q25_v']
    v_q50 = ncfile.variables['q50_v']
    v_q75 = ncfile.variables['q75_v']

    u_binned_obscount = ncfile.variables['binned_obscount_u']
    u_binned_mean = ncfile.variables['binned_mean_u']
    u_binned_max = ncfile.variables['binned_max_u']
    u_binned_min = ncfile.variables['binned_min_u']
    u_binned_std = ncfile.variables['binned_std_u']
    u_binned_rmse = ncfile.variables['binned_rmse_u']

    v_binned_obscount = ncfile.variables['binned_obscount_v']
    v_binned_mean = ncfile.variables['binned_mean_v']
    v_binned_max = ncfile.variables['binned_max_v']
    v_binned_min = ncfile.variables['binned_min_v']
    v_binned_std = ncfile.variables['binned_std_v']
    v_binned_rmse = ncfile.variables['binned_rmse_v']

    idx = len(date)

    subtype[idx] = uv_stat_dict['Subtype']
    date[idx] = uv_stat_dict['Date']
    epoch[idx] = uv_stat_dict['epoch']

    u_obscount[idx] = uv_stat_dict['u']['n']
    u_mean[idx] = uv_stat_dict['u']['mean']
    u_mx[idx] = uv_stat_dict['u']['max']
    u_mn[idx] = uv_stat_dict['u']['min']
    u_std[idx] = uv_stat_dict['u']['std']
    u_rmse[idx] = uv_stat_dict['u']['rmse']
    u_q25[idx] = uv_stat_dict['u']['q25']
    u_q50[idx] = uv_stat_dict['u']['q50']
    u_q75[idx] = uv_stat_dict['u']['q75']

    v_obscount[idx] = uv_stat_dict['v']['n']
    v_mean[idx] = uv_stat_dict['v']['mean']
    v_mx[idx] = uv_stat_dict['v']['max']
    v_mn[idx] = uv_stat_dict['v']['min']
    v_std[idx] = uv_stat_dict['v']['std']
    v_rmse[idx] = uv_stat_dict['v']['rmse']
    v_q25[idx] = uv_stat_dict['v']['q25']
    v_q50[idx] = uv_stat_dict['v']['q50']
    v_q75[idx] = uv_stat_dict['v']['q75']

    u_binned_obscount[idx, :, :, :] = binned_data['u']['binned_nobs']
    u_binned_mean[idx, :, :, :] = binned_data['u']['binned_mean']
    u_binned_max[idx, :, :, :] = binned_data['u']['binned_max']
    u_binned_min[idx, :, :, :] = binned_data['u']['binned_min']
    u_binned_std[idx, :, :, :] = binned_data['u']['binned_std']
    u_binned_rmse[idx, :, :, :] = binned_data['u']['binned_rmse']

    v_binned_obscount[idx, :, :, :] = binned_data['v']['binned_nobs']
    v_binned_mean[idx, :, :, :] = binned_data['v']['binned_mean']
    v_binned_max[idx, :, :, :] = binned_data['v']['binned_max']
    v_binned_min[idx, :, :, :] = binned_data['v']['binned_min']
    v_binned_std[idx, :, :, :] = binned_data['v']['binned_std']
    v_binned_rmse[idx, :, :, :] = binned_data['v']['binned_rmse']

    ncfile.close()

    return


def _get_filename(metadata):

    metadata['Diag Type'] = 'omf' if metadata['Diag Type'] in [
        'o-f', 'omb', 'o-b'] else metadata['Diag Type']
    metadata['Diag Type'] = 'oma' if metadata['Diag Type'] in [
        'o-a'] else metadata['Diag Type']

    if metadata['Diag File Type'] == 'conventional':
        ncfilename = '{Obs Type}_{Variable}_{ObsID[0]}_{Diag Type}.nc'.format(
            **metadata)

    elif metadata['Diag File Type'] == 'ozone':
        layer = '_'.join(metadata['Layer'].split(
        )) if metadata['Layer'] == 'column total' else \
            f"{metadata['Layer']:.3f}"

        ncfilename = '{Obs Type}_{Satellite}_{Diag Type}_'.format(
            **metadata) + '%s.nc' % layer
    else:
        ncfilename = ('{Obs Type}_{Satellite}_Ch{Channel[0]}_'
                      '{Diag Type}.nc').format(**metadata)

    return ncfilename


def write_netcdf(data, binned_data, metadata, outdir):

    if metadata['Diag File Type'] == 'conventional' \
            and metadata['Variable'] == 'uv':

        uv_stat_dict = {}

        u_stat_dict = _calculate_stats(data['u'])
        v_stat_dict = _calculate_stats(data['v'])

        uv_stat_dict['u'] = u_stat_dict
        uv_stat_dict['v'] = v_stat_dict
        uv_stat_dict['Date'] = metadata['Date'].strftime("%Y%m%d%H")
        uv_stat_dict['epoch'] = metadata['Date'].timestamp()
        uv_stat_dict['Subtype'] = int(metadata['Subtype'][0])

        # Get the filename based on metadata
        ncfilename = _get_filename(metadata)

        ncpath = outdir+'/'+ncfilename

        # See if the file exists. If it doesn't, create it
        if not os.path.isfile(ncpath):
            _write_ncfile_uv(uv_stat_dict, binned_data, ncpath)

        # else, append data to that file
        else:
            _append_ncfile_uv(uv_stat_dict, binned_data, ncpath)

    else:

        stat_dict = _calculate_stats(data)

        stat_dict['Date'] = int(metadata['Date'].strftime("%Y%m%d%H"))
        stat_dict['epoch'] = metadata['Date'].timestamp()

        if metadata['Diag File Type'] == 'conventional':
            stat_dict['Subtype'] = int(metadata['Subtype'][0])

        # Get the filename based on metadata
        ncfilename = _get_filename(metadata)

        ncpath = outdir+'/'+ncfilename

        # See if the file exists. If it doesn't, create it
        if not os.path.isfile(ncpath):
            _write_ncfile(stat_dict, binned_data, ncpath, metadata['Obs Type'])

        # else, append data to that file
        else:
            _append_ncfile(stat_dict, binned_data,
                           ncpath, metadata['Obs Type'])

    return
