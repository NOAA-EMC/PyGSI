import os
import pandas as pd
import numpy as np
from netCDF4 import Dataset
from datetime import datetime
from pathlib import Path

_VALID_LVL_TYPES = ["pressure", "height"]
_VALID_CONV_DIAG_TYPES = ["omf", "oma", "observation", "hofx"]
_VALID_RADIANCE_DIAG_TYPES = ["omf", "oma", "observation", "hofx",
                              "water_fraction", "land_fraction",
                              "cloud_fraction", "snow_fraction",
                              "ice_fraction"]


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

    def _query_diag_type(self, df, diag_type, bias_correction):
        """
        Query the data type being requested and returns
        the appropriate indexed data
        """

        bias = 'adjusted' if bias_correction else 'unadjusted'

        if self.variable == 'uv':
            if diag_type in ['Observation']:
                u = df[f'u_{diag_type}']
                v = df[f'v_{diag_type}']
            else:
                u = df[f'u_{diag_type}_{bias}']
                v = df[f'v_{diag_type}_{bias}']

            return u.to_numpy(), v.to_numpy()

        else:
            if diag_type in ['Observation']:
                data = df[f'{diag_type}']
            else:
                data = df[f'{diag_type}_{bias}']

            return data.to_numpy()


class Conventional(GSIdiag):

    def __init__(self, path):
        """
        Initialize a conventional GSI diagnostic object

        Args:
            path : (str) path to conventional GSI diagnostic object
        Returns:
            self : GSI diag conventional object containing the path
                   to extract data
        """
        super().__init__(path)

        self._read_obs()
        self.metadata['Diag File Type'] = 'conventional'

    def __str__(self):
        return "Conventional GSI diagnostic object"

    def _read_obs(self):
        """
        Reads the data from the conventional diagnostic file during
        initialization into a multidimensional pandas dataframe.
        """
        df_dict = {}

        # Open netCDF file, store data into dictionary
        with Dataset(file, mode='r') as f:
            for var in f.variables:
                # Station_ID and Observatio_Class variables need
                # to be converted from byte string to string
                if var in ['Station_ID', 'Observation_Class']:
                    data = f.variables[var][:]
                    data = [i.tobytes(fill_value='/////', order='C')
                            for i in data]
                    data = np.array(
                        [''.join(i.decode('UTF-8', 'ignore').split())
                         for i in data])
                    df_dict[var] = data

                # Grab variables with only 'nobs' dimension
                elif len(f.variables[var].shape) == 1:
                    df_dict[var] = f.variables[var][:]

        # Create pandas dataframe from dict
        df = pd.DataFrame(df_dict)

        # Creates multidimensional indexed dataframe
        indices = ['Station_ID', 'Observation_Class', 'Observation_Type',
                   'Observation_Subtype', 'Pressure', 'Analysis_Use_Flag']
        df.set_index(indices, inplace=True)

        # Rename columns
        df.columns = df.columns.str.lower()
        if self.variable == 'uv':
            df = df.rename(columns={
                'u_obs_minus_forecast_unadjusted': 'u_omf_unadjusted',
                'u_obs_minus_forecast_adjusted': 'u_omf_adjusted',
                'v_obs_minus_forecast_unadjusted': 'v_omf_unadjusted',
                'v_obs_minus_forecast_adjusted': 'v_omf_adjusted'
                })
            # Create hofx columns
            df['u_hofx_unadjusted'] = df['u_observation'] - \
                df['u_omf_unadjusted']
            df['v_hofx_unadjusted'] = df['v_observation'] - \
                df['v_omf_unadjusted']
            df['u_hofx_adjusted'] = df['u_observation'] - \
                df['u_omf_adjusted']
            df['v_hofx_adjusted'] = df['v_observation'] - \
                df['v_omf_adjusted']

        else:
            df = df.rename(columns={
                'obs_minus_forecast_unadjusted': 'omf_unadjusted',
                'obs_minus_forecast_adjusted': 'omf_adjusted',
                })
            # Create hofx columns
            df['hofx_unadjusted'] = df['observation'] - df['omf_unadjusted']
            df['hofx_adjusted'] = df['observation'] - df['omf_adjusted']

        self.data_df = df

    def get_data(self, diag_type, obsid=None, subtype=None, station_id=None,
                 analysis_use=False, lvls=None, lvl_type='pressure',
                 bias_correction=True):
        """
        Given parameters, get the data from a conventional diagnostic file

        Args:
            diag_type : (str; Required) type of data to extract
                        i.e. observation, omf, oma, hofx
            obsid : (list of ints; optional; default=None) observation
                    measurement ID number; default=None
            subtype : (list of ints; optional; default=None) observation
                      measurement ID subtype number, default=None
            station_id : (list of str; optional; default=None)
                         station id, default=None
            analysis_use : (bool; defaul=False) if True, will return
                           three sets of data:
                           assimlated (analysis_use_flag=1, qc<7),
                           rejected (analysis_use_flag=-1, qc<8),
                           monitored (analysis_use_flag=-1, qc>7)
            lvls : (list type; default=None) List of pressure or height
                   levels i.e. [250,500,750,1000]. List must be arranged
                   low to high.  Will return a dictionary with subsetted
                   pressure or height levels where data is separated
                   within those levels:

                   dict = {250-500: <data>,
                           500-750: <data>,
                           750-1000: <data>}
            lvl_type : (str; default='pressure') lvls definition as
                       'pressure' or 'height'.
            bias_correction : (bool; default=True) If True, will return bias
                              corrected data.
        Returns:
            data   : requested indexed data
        """

        if diag_type not in _VALID_CONV_DIAG_TYPES:
            raise ValueError((f'{diag_type} wrong. Valid choices are: '
                              f'{" | ".join(_VALID_DIAG_TYPES)}'))

        self.metadata['Diag Type'] = diag_type
        self.metadata['ObsID'] = obsid
        self.metadata['Subtype'] = subtype
        self.metadata['Station ID'] = station_id
        self.metadata['Anl Use'] = analysis_use
        self.metadata['Levels'] = lvls
        self.metadata['Levels Type'] = lvl_type

        if analysis_use:
            assimilated_df, rejected_df, monitored_df = self._select_conv(
                obsid, subtype, station_id, analysis_use)

            if self.variable == 'uv':
                u_assimilated, v_assimilated = self.query_diag_type(
                    assimilated_df, diag_type, bias_correction)
                u_rejected, v_rejected = self.query_diag_type(
                    rejected_df, diag_type, bias_correction)
                u_monitored, v_monitored = self.query_diag_type(
                    monitored_df, diag_type, bias_correction)

                u = {'assimilated': u_assimilated,
                     'rejected': u_rejected,
                     'monitored': u_monitored}
                v = {'assimilated': v_assimilated,
                     'rejected': v_rejected,
                     'monitored': v_monitored}

                return u, v

            else:
                assimilated_data = self.query_diag_type(
                    assimilated_df, diag_type, bias_correction)
                rejected_data = self.query_diag_type(
                    rejected_df, diag_type, bias_correction)
                monitored_data = self.query_diag_type(
                    monitored_df, diag_type, bias_correction)

                data = {'assimilated': assimilated_data,
                        'rejected': rejected_data,
                        'monitored': monitored_data
                        }

                return data

        else:
            indexed_df = self._select_conv(obsid, subtype, station_id)

            if self.variable == 'uv':
                u, v = self._query_diag_type(
                    indexed_df, diag_type, bias_correction)

                return u, v

            else:
                data = self._query_diag_type(
                    indexed_df, diag_type, bias_correction)

                return data

    def _select_conv(self, obsid=None, subtype=None, station_id=None,
                     analysis_use=False):
        """
        Given parameters, multidimensional dataframe is indexed
        to only include selected locations from a conventional
        diagnostic file.

        Args:
            obsid : (list of ints; default=None) observation measurement
                    ID number
            subtype : (list of ints; default=None) subtype number
            station_id : (list of str; default=None) station id tag
            analysis_use : (bool; deafault=False) if True, will separate
                           into three indexed dataframes: assimilated,
                           rejected, monitored
        Returns:
            df : (pandas dataframe) indexed multidimentsional
                 dataframe from selected data
        """

        df = self.data_df

        # select data by obsid, subtype, and station ids
        if obsid is not None:
            idx_col = 'Observation_Type'
            indx = df.index.get_level_values(idx_col) == ''
            for obid in obsid:
                indx = np.ma.logical_or(
                    indx, df.index.get_level_values(idx_col) == obid)
            df = df.iloc[indx]
        if subtype is not None:
            idx_col = 'Observation_Subtype'
            indx = df.index.get_level_values(idx_col) == ''
            for stype in subtype:
                indx = np.ma.logical_or(
                    indx, df.index.get_level_values(idx_col) == stype)
            df = df.iloc[indx]
        if station_id is not None:
            idx_col = 'Station_ID'
            indx = df.index.get_level_values(idx_col) == ''
            for stn_id in station_id:
                indx = np.ma.logical_or(
                    indx, df.index.get_level_values(idx_col) == stn_id)
            df = df.iloc[indx]

        if analysis_use:
            # Separate into 3 dataframes; assimilated, rejected, and monitored
            indx = df.index.get_level_values('Analysis_Use_Flag') == ''

            assimilated_indx = np.ma.logical_or(
                indx, df.index.get_level_values('Analysis_Use_Flag') == 1)
            rejected_indx = np.ma.logical_or(
                indx, df.index.get_level_values('Analysis_Use_Flag') == -1)
            monitored_indx = np.ma.logical_or(
                indx, df.index.get_level_values('Analysis_Use_Flag') == -1)

            assimilated_df = df.iloc[assimilated_indx]
            rejected_df = df.iloc[rejected_indx]
            monitored_df = df.iloc[monitored_indx]

            # Find rejected and monitored based on Prep_QC_Mark
            try:
                assimilated_df = assimilated_df.loc[
                    assimilated_df['Prep_QC_Mark'] < 7]
                rejected_df = rejected_df.loc[
                    rejected_df['Prep_QC_Mark'] < 8]
                monitored_df = monitored_df.loc[
                    monitored_df['Prep_QC_Mark'] > 7]
            except KeyError:
                assimilated_df = assimilated_df.loc[
                    assimilated_df['Setup_QC_Mark'] < 7]
                rejected_df = rejected_df.loc[
                    rejected_df['Setup_QC_Mark'] < 8]
                monitored_df = monitored_df.loc[
                    monitored_df['Setup_QC_Mark'] > 7]

            return assimilated_df, rejected_df, monitored_df

        else:
            return df
