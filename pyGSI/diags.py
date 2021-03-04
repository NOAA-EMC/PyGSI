import os
import xarray as xr
import pandas as pd
import numpy as np
from netCDF4 import Dataset
from datetime import datetime
from pathlib import Path


class GSIdiag:

    def __init__(self, path):
        """
        Initialize a GSI diagnostic object
        INPUT:
            path : path to GSI diagnostic object
        """

        self.path = path
        self.filename = os.path.splitext(Path(self.path).stem)[0]
        self.obs_type = self.filename.split('_')[1]
        self.variable = self.filename.split('_')[2]
        self.ftype = self.filename.split('_')[-1]

        _str_date = os.path.basename(self.path).split('.')[1]
        # Checks if '_ensmean' is included in file name
        self.date = datetime.strptime(_str_date.split('_')[0], '%Y%m%d%H')

        _var_key_name = 'Variable' if self.obs_type == 'conv' else 'Satellite'
        self.metadata = {'Obs Type': self.obs_type,
                         _var_key_name: self.variable,
                         'Date': self.date,
                         'File Type': self.ftype
                         }

    def __len__(self):
        return len(self.lats)

    def query_diag_type(self, diag_type, indexed_df):
        """
        Query the data type being requested and returns
        the appropriate indexed data
        """

        diag_type = diag_type.lower()

        # Determines the diagnostic type
        omf = True if diag_type in ['omf', 'o-f', 'omb', 'o-b'] else False
        oma = True if diag_type in ['oma', 'o-a'] else False
        obs = True if diag_type in ['observation'] else False
        hofx = True if diag_type in ['hofx'] else False
        water_fraction = True if diag_type in ['water_fraction'] else False
        land_fraction = True if diag_type in ['land_fraction'] else False
        cloud_fraction = True if diag_type in ['cloud_fraction'] else False
        snow_fraction = True if diag_type in ['snow_fraction'] else False
        ice_fraction = True if diag_type in ['ice_fraction'] else False

        # Determines the file type
        is_ges = True if self.ftype == 'ges' else False
        is_anl = True if self.ftype == 'anl' else False

        if omf and is_ges:
            if self.variable == 'uv':
                u = indexed_df['u_Obs_Minus_Forecast_adjusted'].to_numpy()
                v = indexed_df['v_Obs_Minus_Forecast_adjusted'].to_numpy()

                return u, v

            else:
                data = indexed_df['Obs_Minus_Forecast_adjusted'].to_numpy()
                return data

        if oma and is_anl:
            if self.variable == 'uv':
                u = indexed_df['u_Obs_Minus_Forecast_adjusted'].to_numpy()
                v = indexed_df['v_Obs_Minus_Forecast_adjusted'].to_numpy()

                return u, v

            else:
                data = indexed_df['Obs_Minus_Forecast_adjusted'].to_numpy()
                return data

        if obs:
            if self.variable == 'uv':
                u = indexed_df['u_Observation'].to_numpy()
                v = indexed_df['v_Observation'].to_numpy()

                return u, v

            else:
                data = indexed_df['Observation'].to_numpy()
                return data

        if hofx:
            if self.variable == 'uv':
                u = indexed_df['u_Observation'] - \
                    indexed_df['u_Obs_Minus_Forecast_adjusted']
                v = indexed_df['v_Observation'] - \
                    indexed_df['v_Obs_Minus_Forecast_adjusted']

                return u.to_numpy(), v.to_numpy()

            else:
                data = indexed_df['Observation'].to_numpy() - \
                    indexed_df['Obs_Minus_Forecast_adjusted'].to_numpy()
                return data

        if water_fraction:
            data = indexed_df['Water_Fraction'].to_numpy()
            return data

        if land_fraction:
            data = indexed_df['Land_Fraction'].to_numpy()
            return data

        if ice_fraction:
            data = indexed_df['Ice_Fraction'].to_numpy()
            return data

        if snow_fraction:
            data = indexed_df['Snow_Fraction'].to_numpy()
            return data

        if cloud_fraction:
            data = indexed_df['Cloud_Frac'].to_numpy()
            return data

        else:
            raise Exception(f'Unrecognizable diag type input: {diag_type}')
            return None


class Conventional(GSIdiag):

    def __init__(self, path):
        """
        Initialize a conventional GSI diagnostic object
        INPUT:
            path   : path to conventional GSI diagnostic object
        RESULT:
            self   : GSI diag conventional object containing the path to extract data
        """
        super().__init__(path)

        self._read_obs()
        self.metadata['Diag File Type'] = 'conventional'

    def __str__(self):
        return "Conventional GSI diagnostic object"

    def _read_obs(self):
        """
        Reads the data from the conventional diagnostic file during initialization.
        Stores the data as a pandas dataframe.
        """
        ds = xr.open_dataset(self.path)
        df = ds.to_dataframe()

        # Set nobs as index
        df = df.reset_index().set_index('nobs')

        # only found this on t
        if 'Bias_Correction_Terms_arr_dim' in df:
            df = df.loc[df['Bias_Correction_Terms_arr_dim'] == 0]

        self.data_df = df

    def get_data(self, diag_type, obsid=None, subtype=None, station_id=None, analysis_use=False, plvls=None):
        """
        Given parameters, get the data from a conventional diagnostic file
        INPUT:
            required:
                diag_type  : type of data to extract i.e. observation, omf, oma, hofx

            optional:    
                obsid        : observation measurement ID number; default=None
                subtype      : observation measurement ID subtype number, default=None
                station_id   : station id, default=None
                analysis_use : if True, will return two sets of data: assimilated
                               (analysis_use_flag=1), and monitored (analysis_use
                               _flag=-1); default = False
                plvls        : List of pressure levels i.e. [250,500,750,1000]. Will 
                               return a dictionary with subsetted pressure levels where
                               data is seperated within those levels:

                               dict = {250-500: <data>,
                                       500-750: <data>,
                                       750-1000: <data>}

        OUTPUT:
            data   : requested data
        """

        self.metadata['Diag Type'] = diag_type
        self.metadata['ObsID'] = obsid
        self.metadata['Subtype'] = subtype
        self.metadata['Station ID'] = station_id
        self.metadata['Anl Use'] = analysis_use

        indexed_df = self._get_idx_conv(obsid, subtype, station_id)

        if plvls is not None:
            binned_pressure = self._pressure_indexing(
                diag_type, indexed_df, plvls, analysis_use=analysis_use)

            return binned_pressure

        else:
            if analysis_use:
                assimilated_df = indexed_df.loc[df['Analysis_Use_Flag'] == 1]
                monitored_df = indexed_df.loc[df['Analysis_Use_Flag'] == -1]

                if self.variable == 'uv':
                    u_assimilated, v_assimilated = self.query_diag_type(
                        diag_type, assimilated_df)
                    u_monitored, v_monitored = self.query_diag_type(
                        diag_type, monitored_df)

                    u = {'assimilated': u_assimilated,
                         'monitored': u_monitored}
                    v = {'assimilated': v_assimilated,
                         'monitored': v_monitored}

                    return u, v

                else:
                    assimilated_data = self.query_diag_type(
                        diag_type, assimilated_df)
                    monitored_data = self.query_diag_type(
                        diag_type, monitored_df)

                    data = {'assimilated': assimilated_data,
                            'monitored': monitored_data
                            }

                    return data

            else:
                indexed_df = self._get_idx_conv(obsid, subtype, station_id)

                if self.variable == 'uv':
                    u, v = self.query_diag_type(diag_type, indexed_df)

                    return u, v

                else:
                    data = self.query_diag_type(diag_type, indexed_df)

                    return data

    def _get_idx_conv(self, obsid=None, subtype=None, station_id=None):
        """
        Given parameters, get the indices of the observation
        locations from a conventional diagnostic file
        INPUT:
            obsid   : observation measurement ID number
            subtype : subtype number (default: None)
            station_id : station id tag (default: None)
        OUTPUT:
            idx    : indices of the requested data in the file
        """

        indexed_df = self.data_df

        if obsid != None:
            indexed_df = indexed_df.loc[indexed_df['Observation_Type'] == obsid]
        if subtype != None:
            indexed_df = indexed_df.loc[indexed_df['Observation_Subtype'] == subtype]
        if station_id != None:
            # Turns byte string to string and removes white spaces
            indexed_df['Station_ID'] = indexed_df['Station_ID'].str.decode(
                'UTF-8').str.strip()
            indexed_df = indexed_df.loc[indexed_df['Station_ID'] == station_id]

        return indexed_df

    def _pressure_indexing(self, diag_type, df, plvls, analysis_use=False, latlon=False):
        """
        Indexes inputted dataframe into the pressure levels given into a dictionary.
        i.e. if plvls=[250,500,750,1000], returns dictionary:

        dict = {250-500: <data>,
                500-750: <data>,
                750-1000: <data>}
        """

        if latlon:
            pressure_lats = {}
            pressure_lons = {}

        else:
            binned_pressure = {}

        if analysis_use:
            assimilated_df = df.loc[df['Analysis_Use_Flag'] == 1]
            monitored_df = df.loc[df['Analysis_Use_Flag'] == -1]

            for i, pressure in enumerate(plvls[:-1]):
                # Find where pressure is >= to first pressure and < second pressure
                assimilated_pressure_df = assimilated_df.loc[(
                    assimilated_df['Pressure'] >= plvls[i]) & (assimilated_df['Pressure'] < plvls[i+1])]
                monitored_pressure_df = monitored_df.loc[(
                    monitored_df['Pressure'] >= plvls[i]) & (monitored_df['Pressure'] < plvls[i+1])]

                if latlon:
                    lats = {'assimilated': assimilated_pressure_df['Latitude'].to_numpy(),
                            'monitored': monitored_pressure_df['Latitude'].to_numpy()}
                    lons = {'assimilated': assimilated_pressure_df['Longitude'].to_numpy(),
                            'monitored': monitored_pressure_df['Longitude'].to_numpy()}

                    pressure_lats[f'{plvls[i]}-{plvls[i+1]}'] = lats
                    pressure_lons[f'{plvls[i]}-{plvls[i+1]}'] = lons

                elif self.variable == 'uv':
                    u_assimilated, v_assimilated = self.query_diag_type(
                        diag_type, assimilated_pressure_df)
                    u_monitored, v_monitored = self.query_diag_type(
                        diag_type, monitored_pressure_df)

                    data = {'u': {'assimilated': u_assimilated,
                                  'monitored': u_monitored},
                            'v': {'assimilated': v_assimilated,
                                  'monitored': v_monitored}
                            }

                    binned_pressure[f'{plvls[i]}-{plvls[i+1]}'] = data

                else:
                    assimilated_data = self.query_diag_type(
                        diag_type, assimilated_pressure_df)
                    monitored_data = self.query_diag_type(
                        diag_type, monitored_pressure_df)

                    data = {'assimilated': assimilated_data,
                            'monitored': monitored_data
                            }

                    binned_pressure[f'{plvls[i]}-{plvls[i+1]}'] = data

        else:

            for i, pressure in enumerate(plvls[:-1]):
                pressure_df = df.loc[(df['Pressure'] >= plvls[i]) & (
                    df['Pressure'] < plvls[i+1])]

                if latlon:
                    pressure_lats[f'{plvls[i]}-{plvls[i+1]}'] = pressure_df['Latitude'].to_numpy()
                    pressure_lons[f'{plvls[i]}-{plvls[i+1]}'] = pressure_df['Longitude'].to_numpy()

                elif self.variable == 'uv':
                    u, v = self.query_diag_type(diag_type, pressure_df)

                    data = {'u': u,
                            'v': v}

                    binned_pressure[f'{plvls[i]}-{plvls[i+1]}'] = data

                else:
                    data = self.query_diag_type(diag_type, pressure_df)

                    binned_pressure[f'{plvls[i]}-{plvls[i+1]}'] = data

        if latlon:
            return pressure_lats, pressure_lons

        else:
            return binned_pressure

    def get_lat_lon(self, obsid=None, subtype=None, station_id=None, analysis_use=False, plvls=None):
        """
        Gets lats and lons with desired indices
        """

        indexed_df = self._get_idx_conv(obsid, subtype, station_id)

        if plvls is not None:
            binned_lat, binned_lon = self._pressure_indexing(diag_type, indexed_df, plvls,
                                                             analysis_use=analysis_use,
                                                             latlon=True)

            return binned_lat, binned_lon

        else:
            if analysis_use:
                assimilated_df, monitored_df = self._get_idx_conv(
                    obsid, subtype, station_id, analysis_use)

                lats = {'assimilated': assimilated_pressure_df['Latitude'].to_numpy(),
                        'monitored': monitored_pressure_df['Latitude'].to_numpy()}
                lons = {'assimilated': assimilated_pressure_df['Longitude'].to_numpy(),
                        'monitored': monitored_pressure_df['Longitude'].to_numpy()}
                return lats, lons

            else:
                indexed_df = self._get_idx_conv(
                    obsid, subtype, station_id, analysis_use)

                return indexed_df['Latitude'].to_numpy(), indexed_df['Longitude'].to_numpy()

    def get_pressure(self, obsid=None, subtype=None, station_id=None, analysis_use=False):

        if analysis_use:
            assimilated_df, monitored_df = self._get_idx_conv(
                obsid, subtype, station_id, analysis_use)
            pressure = {'assimilated': assimilated_df['Pressure'].to_numpy(),
                        'monitored': monitored_df['Pressure'].to_numpy()}

            return pressure
        else:
            indexed_df = self._get_idx_conv(
                obsid, subtype, station_id, analysis_use)
            return indexed_df['Pressure'].to_numpy()


class Radiance(GSIdiag):

    def __init__(self, path):
        """
        Initialize a radiance GSI diagnostic object
        INPUT:
            path   : path to conventional GSI diagnostic object
        RESULT:
            self   : GSI diag radiance object containing the path to extract data
        """
        super().__init__(path)

        self._read_obs()
        self.metadata['Diag File Type'] = 'radiance'

    def __str__(self):
        return "Radiance GSI diagnostic object"

    def _read_obs(self):
        """
        Reads the data from the radiance diagnostic file during initialization.
        """

        with Dataset(self.path, mode='r') as f:
            self.channel_idx = f.variables['Channel_Index'][:]
            self.sensor_chan = f.variables['sensor_chan'][:]
            self.chaninfo_idx = f.variables['chaninfoidx'][:]
            self.lons = f.variables['Longitude'][:]
            self.lats = f.variables['Latitude'][:]
            self.time = f.variables['Obs_Time'][:]
            self.o = f.variables['Observation'][:]
            self.omf = f.variables['Obs_Minus_Forecast_adjusted'][:]
            self.qc_flag = f.variables['QC_Flag'][:]
            self.water_frac = f.variables['Water_Fraction'][:]
            self.land_frac = f.variables['Land_Fraction'][:]
            self.ice_frac = f.variables['Ice_Fraction'][:]
            self.snow_frac = f.variables['Snow_Fraction'][:]
            self.cloud_frac = f.variables['Cloud_Frac'][:]
            self.inv_ob_err = f.variables['Inverse_Observation_Error'][:]

    def get_data(self, diag_type, channel=None, qcflag=None,
                 separate_channels=False, separate_qc=False, errcheck=True):
        """
        Given parameters, get the data from a radiance 
        diagnostic file.
        INPUT:
            required:
                diag_type : type of data to extract i.e. observation, omf, oma, hofx, water_fraction,
                        land_fraction, ice_fraction, snow_fraction, cloud_fraction

            optional:  
                channel           : observation channel number
                qcflag            : qc flag (default: None) i.e. 0, 1
                separate_channels : if True, calls _get_data_special() and returns dictionary
                                    of separate data by specified channels
                separate_qc       : if True, calls _get_data_special() and returns dictionary
                                    of separate data by specified QC flags
                errcheck          : when true, and qc==0, will toss out obs where inverse
                                    obs error is zero (i.e. not assimilated in GSI)
        OUTPUT:
            data : requested data

        """

        self.metadata['Diag Type'] = diag_type
        self.metadata['QC Flag'] = qcflag
        self.metadata['Channels'] = channel

        if separate_channels or separate_qc:
            data = self._get_data_special(
                diag_type, channel, qcflag, separate_channels, separate_qc, errcheck=errcheck)
            return data

        else:
            idx = self._get_idx_sat(channel, qcflag, errcheck=errcheck)

            data = self.query_diag_type(diag_type, idx)

            data[data > 1e5] = np.nan

            return data

    def _get_idx_sat(self, channel=None, qcflag=None, errcheck=True):
        """
        Given parameters, get the indices of the observation
        locations from a radiance diagnostic file.
        """

        if not channel or channel == [None]:
            channel = None
        if not qcflag or qcflag == [None]:
            qcflag = None

        idx = self.channel_idx
        valid_idx = np.full_like(idx, True, dtype=bool)
        if qcflag != None:
            idxqc = self.qc_flag
            valid_idx_qc = np.isin(idxqc, qcflag)
            valid_idx = np.logical_and(valid_idx, valid_idx_qc)
        if channel != None:
            chidx = np.where(self.sensor_chan == channel)
            if len(chidx) > 0 and len(chidx[0]) > 0:
                channel = self.chaninfo_idx[chidx[0][0]]
            else:
                print('Channel specified not in sensor_chan, using relative index')
            valid_idx_ch = np.isin(self.channel_idx, channel)
            valid_idx = np.logical_and(valid_idx, valid_idx_ch)
        if errcheck and qcflag == 0:
            valid_idx_err = np.isin(self.inv_ob_err, 0, invert=True)
            valid_idx = np.logical_and(valid_idx, valid_idx_err)

        idx = np.where(valid_idx)
        return idx

    def _get_data_special(self, diag_type, channel, qcflag,
                          separate_channels, separate_qc, errcheck=True):
        """
        Creates a dictionary that separates channels and qc flags
        depending on the conditions of seperate_channels and
        separate_qc
        """
        data_dict = {}

        if separate_channels and not separate_qc:
            for c in channel:
                idx = self._get_idx_sat(c, qcflag, errcheck=errcheck)

                data = self.query_diag_type(diag_type, idx)

                data[data > 1e5] = np.nan

                data_dict['Channel %s' % c] = data

            return data_dict

        if not separate_channels and separate_qc:
            for qc in qcflag:
                idx = self._get_idx_sat(channel, qc, errcheck=errcheck)

                data = self.query_diag_type(diag_type, idx)

                data[data > 1e5] = np.nan

                data_dict['QC Flag %s' % qc] = data

            return data_dict

        if separate_channels and separate_qc:
            for c in channel:
                data_dict['Channel %s' % c] = {}
                for qc in qcflag:
                    idx = self._get_idx_sat(c, qc, errcheck=errcheck)

                    data = self.query_diag_type(diag_type, idx)

                    data[data > 1e5] = np.nan

                    data_dict['Channel %s' % c]['QC Flag %s' % qc] = data

            return data_dict

    def get_lat_lon(self, channel=None, qcflag=None, errcheck=True):
        """
        Gets lats and lons with desired indices
        """
        idx = self._get_idx_sat(channel, qcflag, errcheck=errcheck)
        return self.lats[idx], self.lons[idx]


class Ozone(GSIdiag):

    def __init__(self, path):
        """
        Initialize an ozone GSI diagnostic object
        INPUT:
            path   : path to ozone GSI diagnostic object
        RESULT:
            self   : GSI ozone diag object containing the path to extract data
        """
        super().__init__(path)

        self._read_obs()
        self.metadata['Diag File Type'] = 'ozone'

    def __str__(self):
        return "Ozone GSI diagnostic object"

    def _read_obs(self):
        """
        Reads the data from the ozone diagnostic file during initialization.
        """

        with Dataset(self.path, mode='r') as f:
            self.lons = f.variables['Longitude'][:]
            self.lats = f.variables['Latitude'][:]
            self.ref_pressure = f.variables['Reference_Pressure'][:]
            self.time = f.variables['Time'][:]
            self.anl_use = f.variables['Analysis_Use_Flag'][:]
            self.o = f.variables['Observation'][:]
            self.omf = f.variables['Obs_Minus_Forecast_adjusted'][:]
            self.inv_ob_err = f.variables['Inverse_Observation_Error'][:]

    def get_data(self, diag_type, analysis_use=False, errcheck=True):
        """
        Given parameters, get the data from an ozone diagnostic file
        INPUT:
            required:
                data_type  : type of data to extract i.e. observation, omf, oma, hofx

            optional:    
                analysis_use : if True, will return two sets of data: assimlated
                               (analysis_use_flag=1), and monitored (analysis_use
                               _flag=-1); default = False
                errcheck     : when True, will toss out obs where inverse
                               obs error is zero (i.e. not assimilated in GSI)

        OUTPUT:
            data   : requested data as a dictionary
        """

        self.metadata['Diag Type'] = diag_type
        self.metadata['Anl Use'] = analysis_use

        data_dict = {}

        if not analysis_use:

            for layer in self.ref_pressure:
                if layer != 0:
                    layer_idx = np.isin(self.ref_pressure, layer)
                    idx = self._get_idx_ozone(layer_idx, errcheck=errcheck)
                    data = self.query_diag_type(diag_type, idx)

                    data_dict[layer] = data
                else:
                    break

            column_total_idx = np.isin(self.ref_pressure, 0)
            idx = self._get_idx_ozone(column_total_idx, errcheck=errcheck)

            data = self.query_diag_type(diag_type, idx)
            data_dict['column total'] = data

        else:
            for layer in self.ref_pressure:
                if layer != 0.:
                    layer_idx = np.isin(self.ref_pressure, layer)

                    assimilated_idx, monitored_idx = self._get_idx_ozone(
                        layer_idx, analysis_use=analysis_use, errcheck=errcheck)

                    assimilated_data = self.query_diag_type(
                        diag_type, assimilated_idx)
                    monitored_data = self.query_diag_type(
                        diag_type, monitored_idx)

                    data_dict[layer] = {'assimilated': assimilated_data,
                                        'monitored': monitored_data
                                        }
                else:
                    break

            column_total_idx = np.isin(self.ref_pressure, 0)

            assimilated_idx, monitored_idx = self._get_idx_ozone(
                column_total_idx, analysis_use=analysis_use, errcheck=errcheck)

            assimilated_data = self.query_diag_type(diag_type, assimilated_idx)
            monitored_data = self.query_diag_type(diag_type, monitored_idx)

            data_dict['column total'] = {'assimilated': assimilated_data,
                                         'monitored': monitored_data
                                         }

        return data_dict

    def _get_idx_ozone(self, index, analysis_use=False, errcheck=True):
        """
        Returns the index of data that was assimilated
        and monitored. 
        INPUT:
            index    : array of booleans that is used to
                       to find where they match assimilated
                       or monitored data
            errcheck : when True, will toss out obs where inverse
                       obs error is zero (i.e. not assimilated in GSI)
        OUTPUT:
            assimilated_idx, monitored_idx : if analysis_use is True
            index : if analysis_use is False
        """

        if analysis_use:

            valid_assimilated_idx = np.isin(self.anl_use, 1)
            valid_monitored_idx = np.isin(self.anl_use, -1)

            if errcheck:
                valid_idx_err = np.isin(self.inv_ob_err, 0, invert=True)
                index = np.logical_and(index, valid_idx_err)

            valid_assimilated_idx = np.logical_and(
                valid_assimilated_idx, index)
            valid_monitored_idx = np.logical_and(
                valid_monitored_idx, index)

            assimilated_idx = np.where(valid_assimilated_idx)
            monitored_idx = np.where(valid_monitored_idx)

            return assimilated_idx, monitored_idx

        else:
            if errcheck:
                valid_idx_err = np.isin(self.inv_ob_err, 0, invert=True)
                index = np.logical_and(index, valid_idx_err)

            return index

    def get_lat_lon(self, analysis_use=False, errcheck=True):
        """
        Gets lats and lons with desired indices. Returns dictionary for
        each level in the ozone data.
        """

        lats_dict = {}
        lons_dict = {}

        for layer in self.ref_pressure:
            if layer != 0:
                layer_idx = np.isin(self.ref_pressure, layer)

                if analysis_use:
                    assimilated_idx, monitored_idx = self._get_idx_ozone(
                        layer_idx, analysis_use=analysis_use, errcheck=errcheck)

                    lats = {'assimilated': self.lats[assimilated_idx],
                            'monitored': self.lats[monitored_idx]
                            }
                    lons = {'assimilated': self.lons[assimilated_idx],
                            'monitored': self.lons[monitored_idx]
                            }

                else:
                    idx = self._get_idx_ozone(layer_idx, errcheck=errcheck)
                    lats = self.lats[idx]
                    lons = self.lons[idx]

                lats_dict[layer] = lats
                lons_dict[layer] = lons
            else:
                break

        column_total_idx = np.isin(self.ref_pressure, 0)

        if analysis_use:
            assimilated_idx, monitored_idx = self._get_idx_ozone(
                column_total_idx, analysis_use=analysis_use, errcheck=errcheck)

            lats = {'assimilated': self.lats[assimilated_idx],
                    'monitored': self.lats[monitored_idx]
                    }
            lons = {'assimilated': self.lons[assimilated_idx],
                    'monitored': self.lons[monitored_idx]
                    }
        else:

            idx = self._get_idx_ozone(column_total_idx, errcheck=errcheck)

            lats = self.lats[idx]
            lons = self.lons[idx]

        lats_dict['column total'] = lats
        lons_dict['column total'] = lons

        return lats_dict, lons_dict
