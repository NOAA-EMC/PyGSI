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

    def get_diag_type_data(self, diag_type, df, bias_correction=False):
        """
        Based on the inputted diag_type ('omf', 'observation', 'hofx', etc.), returns
        the subsetted data from the inputted dataframe.
        """
        
        diag_type = diag_type.lower()
        
        if diag_type in ['o-f', 'omb', 'o-b', 'o-a', 'oma']:
            diag_type='omf'
            
        bias = 'adjusted' if bias_correction else 'unadjusted'
        
        try:
            if self.variable == 'uv':
                try:
                    u = df[f'u_{diag_type}_{bias}']
                    v = df[f'v_{diag_type}_{bias}']
                except:
                    u = df[f'u_{diag_type}']
                    v = df[f'v_{diag_type}']
                    
                return u, v
            else:
                try:
                    data = df[f'{diag_type}_{bias}']
                except:
                    data = df[f'{diag_type}']
                
                return data
        except:
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
        df = ds.to_dataframe().reset_index()

        # only found this on t
        if 'Bias_Correction_Terms_arr_dim' in df:
            df = df.loc[df['Bias_Correction_Terms_arr_dim'] == 0]
        
        # Convert byte strings to strings
        byte_cols = ['Station_ID', 'Observation_Class']
        for col in byte_cols:
            df[col] = df[col].str.decode("utf-8").str.strip()

        # Creates multidimensional indexed dataframe
        indices = ['Station_ID', 'Observation_Class', 'Observation_Type',
                   'Observation_Subtype', 'Pressure', 'Analysis_Use_Flag']
        df.set_index(indices, inplace=True)

        # Rename columns
        df.columns= df.columns.str.lower()
        if self.variable == 'uv':
            df = df.rename(columns={'u_obs_minus_forecast_unadjusted': 'u_omf_unadjusted',
                                    'u_obs_minus_forecast_adjusted'  : 'u_omf_adjusted',
                                    'v_obs_minus_forecast_unadjusted': 'v_omf_unadjusted',
                                    'v_obs_minus_forecast_adjusted'  : 'v_omf_adjusted'
                                   })
            # Create hofx columns
            df['u_hofx_unadjusted'] = df['u_observation'] - df['u_omf_unadjusted']
            df['v_hofx_unadjusted'] = df['v_observation'] - df['v_omf_unadjusted']
            df['u_hofx_adjusted'] = df['u_observation'] - df['u_omf_adjusted']
            df['v_hofx_adjusted'] = df['v_observation'] - df['v_omf_adjusted']

        else:
            df = df.rename(columns={'obs_minus_forecast_unadjusted': 'omf_unadjusted',
                                    'obs_minus_forecast_adjusted'  : 'omf_adjusted',
                                   }) 
            # Create hofx columns
            df['hofx_unadjusted'] = df['observation'] - df['omf_unadjusted']
            df['hofx_adjusted'] = df['observation'] - df['omf_adjusted']
            
        self.data_df = df
        
        
    def get_data(self, diag_type, obsid=None, subtype=None, station_id=None,
                 analysis_use=False, bias_correction=False):
        """
        Given parameters, get the data from a conventional diagnostic file
        INPUT:
        required:
            diag_type  : type of data to extract i.e. observation, omf, oma, hofx

        optional:    
            obsid        : observation measurement ID number; default=None
            subtype      : observation measurement ID subtype number, default=None
            station_id   : station id, default=None
            analysis_use : if True, will return two sets of data: assimlated
                           (analysis_use_flag=1), and monitored (analysis_use
                           _flag=-1); default = False
            bias_correction : if True, will return '<diag_type>_adjusted' data
                              which includes bias correction; default = False
        OUTPUT:
            data   : requested data
        """
        
        # Store metadata 
        self.metadata['Diag Type'] = diag_type
        self.metadata['ObsID'] = obsid
        self.metadata['Subtype'] = subtype
        self.metadata['Station ID'] = station_id
        self.metadata['Anl Use'] = analysis_use
        
        if analysis_use:
            assimilated_df, monitored_df = self._select_conv(obsid, subtype, station_id, analysis_use)
            
            if self.variable == 'uv':
                u_assimilated, v_assimilated = self.get_diag_type_data(
                    diag_type, assimilated_df, bias_correction)
                u_monitored, v_monitored = self.get_diag_type_data(
                    diag_type, monitored_df, bias_correction)

                u = {'assimilated': u_assimilated,
                     'monitored': u_monitored}
                v = {'assimilated': v_assimilated,
                     'monitored': v_monitored}

                return u, v
            else:
                assimilated_data = self.get_diag_type_data(
                    diag_type, assimilated_df, bias_correction)

                monitored_data = self.get_diag_type_data(
                    diag_type, monitored_df, bias_correction)

                data = {'assimilated': assimilated_data,
                        'monitored': monitored_data}

                return data
        
        else:
            df = self._select_conv(obsid, subtype, station_id, analysis_use)
            
            if self.variable == 'uv':
                u, v = self.get_diag_type_data(diag_type, df, bias_correction)
                
                return u, v
            else:
                data = self.get_diag_type_data(diag_type, df, bias_correction)

                return data
    
    
    def _select_conv(self, obsids=None, subtypes=None, station_ids=None, analysis_use=False):
        """
        Slice a dataframe given input parameters obsid, subtype, station_id,
        and analysis_use
        """
        
        df = self.data_df
        
        if obsids is not None:
            indx = df.index.get_level_values('Observation_Type') == ''
            for obsid in obsids:
                indx = np.ma.logical_or(indx, df.index.get_level_values('Observation_Type') == obsid)
            df = df.iloc[indx]
        if subtypes is not None:
            indx = df.index.get_level_values('Observation_Subtype') == ''
            for subtype in subtypes:
                indx = np.ma.logical_or(indx, df.index.get_level_values('Observation_Subtype') == subtype)
            df = df.iloc[indx]
        if station_ids is not None:
            indx = df.index.get_level_values('Station_ID') == ''
            for station_id in station_ids:
                indx = np.ma.logical_or(indx, df.index.get_level_values('Station_ID') == station_id)
            df = df.iloc[indx]
        
        if analysis_use:
            indx = df.index.get_level_values('Analysis_Use_Flag') == ''
            assimilated_indx = np.ma.logical_or(indx, df.index.get_level_values('Analysis_Use_Flag') == 1)
            monitored_indx = np.ma.logical_or(indx, df.index.get_level_values('Analysis_Use_Flag') == -1)
            
            assimilated_df = df.iloc[assimilated_indx]
            monitored_df = df.iloc[monitored_indx]
            
            return assimilated_df, monitored_df
        
        else:
            return df
        
    
    def get_lat_lon(self, obsid=None, subtype=None, station_id=None, analysis_use=False):
        """
        Gets lats and lons with desired indices
        """
            
        if analysis_use:
            assimilated_df, monitored_df = self._select_conv(
                obsid, subtype, station_id, analysis_use)

            lats = {'assimilated': assimilated_df['latitude'],
                    'monitored': monitored_df['latitude']}
            lons = {'assimilated': assimilated_df['longitude'],
                    'monitored': monitored_df['longitude']}
            return lats, lons

        else:
            df = self._select_conv(obsid, subtype, station_id, analysis_use)

            return df['latitude'], df['longitude']
        
        
    def pressure_binning(self, data, p1=None, p2=None):
        """
        Returns a series of data that is indexed between two 
        pressures
        Inputs:
            data : series of data
            p1   : lower pressure (i.e. 250)
            p2   : higher pressure (i.e. 500)
        Outputs:
            data : data indexed between p1 and p2
        """
        indx = data.index.get_level_values('Pressure')
        indx = np.where((indx >= p1) & (indx < p2))
        data = data.iloc[indx]

        return data


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
