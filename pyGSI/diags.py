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
                              "ice_fraction", "vegetation_fraction"]
_VALID_OZONE_DIAG_TYPES = ["omf", "oma", "observation", "hofx"]


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
        return self.data_df.shape[0]

    def _query_diag_type(self, df, diag_type, bias_correction):
        """
        Query the data type being requested and returns
        the appropriate indexed data
        """

        bias = 'adjusted' if bias_correction else 'unadjusted'
        diag_type = 'omf' if diag_type in ['oma'] else diag_type

        if self.variable == 'uv':
            if diag_type in ['observation']:
                u = df[f'u_{diag_type}']
                v = df[f'v_{diag_type}']
            else:
                u = df[f'u_{diag_type}_{bias}']
                v = df[f'v_{diag_type}_{bias}']

            return u.to_numpy(), v.to_numpy()

        else:
            # uv is not a variable for radiance so the fraction
            # variables only need to be considered here
            if diag_type in ['observation', "water_fraction",
                             "land_fraction", "cloud_fraction",
                             "snow_fraction", "ice_fraction",
                             "vegetation_fraction"]:
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
        with Dataset(self.path, mode='r') as f:
            for var in f.variables:
                # Station_ID and Observation_Class variables need
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
                   'Observation_Subtype', 'Pressure', 'Height',
                   'Analysis_Use_Flag']
        df.set_index(indices, inplace=True)

        # Rename columns
        df.columns = df.columns.str.lower()
        if self.variable == 'uv':
            for wind_type in ['u', 'v']:
                for bias_type in ['unadjusted', 'adjusted']:
                    df = df.rename(columns={
                        f'{wind_type}_obs_minus_forecast_{bias_type}':
                        f'{wind_type}_omf_{bias_type}'
                        })
                    # Create hofx columns
                    df[f'{wind_type}_hofx_{bias_type}'] = \
                        df[f'{wind_type}_observation'] - \
                        df[f'{wind_type}_omf_{bias_type}']

        else:
            for bias_type in ['unadjusted', 'adjusted']:
                df = df.rename(columns={
                    f'obs_minus_forecast_{bias_type}': f'omf_{bias_type}',
                    })
                # Create hofx columns
                df[f'hofx_{bias_type}'] = df['observation'] - \
                    df[f'omf_{bias_type}']

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
                    type ID number; default=None
            subtype : (list of ints; optional; default=None) observation
                      measurement ID subtype number, default=None
            station_id : (list of str; optional; default=None)
                         station id, default=None
            analysis_use : (bool; default=False) if True, will return
                           three sets of data:
                           assimilated (analysis_use_flag=1, qc<7),
                           rejected (analysis_use_flag=-1, qc<8),
                           monitored (analysis_use_flag=-1, qc>7)
            lvls : (list type; default=None) List of pressure or height
                   levels i.e. [250,500,750,1000]. List must be arranged
                   low to high.
                   For pressure, will return a dictionary of data with data
                   greater than the low bound, and less than or equal to the
                   high bound.
                   For height, will return a dictionary of data with data
                   greater than or equal to the low bound, and less than the
                   high bound.
            lvl_type : (str; default='pressure') lvls definition as
                       'pressure' or 'height'.
            bias_correction : (bool; default=True) If True, will return bias
                              corrected data.
        Returns:
            data : requested indexed data
        """

        if diag_type not in _VALID_CONV_DIAG_TYPES:
            raise ValueError((f'{diag_type} is not a valid diag_type. '
                              'Valid choices are: '
                              f'{" | ".join(_VALID_CONV_DIAG_TYPES)}'))

        self.metadata['Diag Type'] = diag_type
        self.metadata['ObsID'] = obsid
        self.metadata['Subtype'] = subtype
        self.metadata['Station ID'] = station_id
        self.metadata['Anl Use'] = analysis_use
        self.metadata['Levels'] = lvls
        self.metadata['Levels Type'] = lvl_type

        # Selects proper levels
        if lvls is not None:
            # Check if level type is valid
            if lvl_type not in _VALID_LVL_TYPES:
                raise ValueError((f'{lvl_type} is not a valid lvl_type. '
                                  'Valid choices are: '
                                  f'{" | ".join(_VALID_LVL_TYPES)}'))

            data = self._get_lvl_data(
                diag_type, obsid, subtype, station_id,
                analysis_use, lvls, lvl_type, bias_correction)

            return data

        else:
            if analysis_use:
                assimilated_df, rejected_df, monitored_df = self._select_conv(
                    obsid, subtype, station_id, analysis_use)

                if self.variable == 'uv':
                    u_assimilated, v_assimilated = self._query_diag_type(
                        assimilated_df, diag_type, bias_correction)
                    u_rejected, v_rejected = self._query_diag_type(
                        rejected_df, diag_type, bias_correction)
                    u_monitored, v_monitored = self._query_diag_type(
                        monitored_df, diag_type, bias_correction)

                    u = {'assimilated': u_assimilated,
                         'rejected': u_rejected,
                         'monitored': u_monitored}
                    v = {'assimilated': v_assimilated,
                         'rejected': v_rejected,
                         'monitored': v_monitored}

                    return u, v

                else:
                    assimilated_data = self._query_diag_type(
                        assimilated_df, diag_type, bias_correction)
                    rejected_data = self._query_diag_type(
                        rejected_df, diag_type, bias_correction)
                    monitored_data = self._query_diag_type(
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
            obsid : (list of ints; default=None) observation
                    type ID number
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
                    assimilated_df['prep_qc_mark'] < 7]
                rejected_df = rejected_df.loc[
                    rejected_df['prep_qc_mark'] < 8]
                monitored_df = monitored_df.loc[
                    monitored_df['prep_qc_mark'] > 7]
            except KeyError:
                assimilated_df = assimilated_df.loc[
                    assimilated_df['setup_qc_mark'] < 7]
                rejected_df = rejected_df.loc[
                    rejected_df['setup_qc_mark'] < 8]
                monitored_df = monitored_df.loc[
                    monitored_df['setup_qc_mark'] > 7]

            return assimilated_df, rejected_df, monitored_df

        else:
            return df

    def _get_lvl_data(self, diag_type, obsid=None, subtype=None,
                      station_id=None, analysis_use=False, lvls=None,
                      lvl_type='pressure', bias_correction=True):
        """
        Given a list of levels, will create a dictionary of data that is
        selected between each level. Will return a dictionary with subsetted
        pressure or height levels where data is separated within those levels:

        dict = {250-500: <data>,
                500-750: <data>,
                750-1000: <data>}
        """
        binned_data = {}

        for i, low_bound in enumerate(lvls[:-1]):
            high_bound = lvls[i+1]

            if analysis_use:
                assimilated_df, rejected_df, monitored_df = self._select_conv(
                    obsid, subtype, station_id, analysis_use)

                assimilated_lvl_df = self._select_levels(
                    assimilated_df, low_bound, high_bound, lvl_type)
                rejected_lvl_df = self._select_levels(
                    rejected_df, low_bound, high_bound, lvl_type)
                monitored_lvl_df = self._select_levels(
                    monitored_df, low_bound, high_bound, lvl_type)

                if self.variable == 'uv':
                    u_assimilated, v_assimilated = self._query_diag_type(
                        assimilated_lvl_df, diag_type, bias_correction)
                    u_rejected, v_rejected = self._query_diag_type(
                        rejected_lvl_df, diag_type, bias_correction)
                    u_monitored, v_monitored = self._query_diag_type(
                        monitored_lvl_df, diag_type, bias_correction)

                    u = {'assimilated': u_assimilated,
                         'rejected': u_rejected,
                         'monitored': u_monitored}
                    v = {'assimilated': v_assimilated,
                         'rejected': v_rejected,
                         'monitored': v_monitored}

                    data = {'u': u,
                            'v': v}
                else:
                    assimilated_data = self._query_diag_type(
                        assimilated_lvl_df, diag_type, bias_correction)
                    rejected_data = self._query_diag_type(
                        rejected_lvl_df, diag_type, bias_correction)
                    monitored_data = self._query_diag_type(
                        monitored_lvl_df, diag_type, bias_correction)

                    data = {'assimilated': assimilated_data,
                            'rejected': rejected_data,
                            'monitored': monitored_data
                            }

                binned_data[f'{low_bound}-{high_bound}'] = data

            else:
                indexed_df = self._select_conv(obsid, subtype, station_id)

                lvl_df = self._select_levels(
                    indexed_df, low_bound, high_bound, lvl_type)

                if self.variable == 'uv':
                    u, v = self._query_diag_type(
                        lvl_df, diag_type, bias_correction)

                    data = {'u': u,
                            'v': v}
                else:
                    data = self._query_diag_type(
                        lvl_df, diag_type, bias_correction)

                binned_data[f'{low_bound}-{high_bound}'] = data

        return binned_data

    def _select_levels(self, df, low_bound, high_bound, lvl_type):
        """
        Selects data between two level bounds from given dataframe.
        """

        if lvl_type == 'pressure':
            # Grab data greater than low bound and less than or
            # equal to high_bound
            df = df.query(
                f'(Pressure > {low_bound}) and (Pressure <= {high_bound})')

        else:
            # Grab data greater than or equal to low bound and
            # less than high_bound
            df = df.query(
                f'(Height >= {low_bound}) and (Height < {high_bound})')

        return df

    def get_pressure(self, obsid=None, subtype=None, station_id=None,
                     analysis_use=False):
        """
        Grabs indexed pressure data.

        Args:
            obsid : (list of ints; optional; default=None) observation
                    type ID number; default=None
            subtype : (list of ints; optional; default=None) observation
                      measurement ID subtype number, default=None
            station_id : (list of str; optional; default=None)
                         station id, default=None
            analysis_use : (bool; default=False) if True, will return
                           three sets of data:
                           assimilated (analysis_use_flag=1, qc<7),
                           rejected (analysis_use_flag=-1, qc<8),
                           monitored (analysis_use_flag=-1, qc>7)
        Returns:
            height : indexed pressure values
        """
        if analysis_use:
            assimilated_df, rejected_df, monitored_df = self._select_conv(
                    obsid, subtype, station_id, analysis_use)

            pressure = {
                'assimilated':
                    assimilated_df.reset_index()['Pressure'].to_numpy(),
                'rejected': rejected_df.reset_index()['Pressure'].to_numpy(),
                'monitored': monitored_df.reset_index()['Pressure'].to_numpy()
            }

        else:
            indexed_df = self._select_conv(obsid, subtype, station_id)
            pressure = indexed_df.reset_index()['Pressure'].to_numpy()

        return pressure

    def get_height(self, obsid=None, subtype=None, station_id=None,
                   analysis_use=False):
        """
        Grabs indexed height data.

        Args:
            obsid : (list of ints; optional; default=None) observation
                    type ID number; default=None
            subtype : (list of ints; optional; default=None) observation
                      measurement ID subtype number, default=None
            station_id : (list of str; optional; default=None)
                         station id, default=None
            analysis_use : (bool; default=False) if True, will return
                           three sets of data:
                           assimilated (analysis_use_flag=1, qc<7),
                           rejected (analysis_use_flag=-1, qc<8),
                           monitored (analysis_use_flag=-1, qc>7)
        Returns:
            height : indexed height values
        """
        if analysis_use:
            assimilated_df, rejected_df, monitored_df = self._select_conv(
                    obsid, subtype, station_id, analysis_use)

            height = {
                'assimilated':
                    assimilated_df.reset_index()['Height'].to_numpy(),
                'rejected': rejected_df.reset_index()['Height'].to_numpy(),
                'monitored': monitored_df.reset_index()['Height'].to_numpy()
            }

        else:
            indexed_df = self._select_conv(obsid, subtype, station_id)
            height = indexed_df.reset_index()['Height'].to_numpy()

        return height

    def get_lat_lon(self, obsid=None, subtype=None, station_id=None,
                    analysis_use=False, lvls=None, lvl_type='pressure'):
        """
        Grabs indexed lats and lons from inputs.

        Args:
            obsid : (list of ints; optional; default=None) observation
                    type ID number; default=None
            subtype : (list of ints; optional; default=None) observation
                      measurement ID subtype number, default=None
            station_id : (list of str; optional; default=None)
                         station id, default=None
            analysis_use : (bool; default=False) if True, will return
                           three sets of data:
                           assimilated (analysis_use_flag=1, qc<7),
                           rejected (analysis_use_flag=-1, qc<8),
                           monitored (analysis_use_flag=-1, qc>7)
            lvls : (list type; default=None) List of pressure or height
                   levels i.e. [250,500,750,1000]. List must be arranged
                   low to high.
                   For pressure, will return a dictionary of data with data
                   greater than the low bound, and less than or equal to the
                   high bound.
                   For height, will return a dictionary of data with data
                   greater than or equal to the low bound, and less than the
                   high bound.
            lvl_type : (str; default='pressure') lvls definition as
                       'pressure' or 'height'.
        Returns:
            lat, lon : (array like) requested indexed latitude and longitude
        """
        # Selects proper levels
        if lvls is not None:
            # Check if level type is valid
            if lvl_type not in _VALID_LVL_TYPES:
                raise ValueError((f'{lvl_type} is not a valid lvl_type. '
                                  'Valid choices are: '
                                  f'{" | ".join(_VALID_LVL_TYPES)}'))

        if analysis_use:
            assimilated_df, rejected_df, monitored_df = self._select_conv(
                        obsid, subtype, station_id, analysis_use)

            # select by levels
            if lvls is not None:
                binned_lats = {}
                binned_lons = {}

                for i, low_bound in enumerate(lvls[:-1]):
                    high_bound = lvls[i+1]

                    assimilated_lvl_df = self._select_levels(
                        assimilated_df, low_bound, high_bound, lvl_type)
                    rejected_lvl_df = self._select_levels(
                        rejected_df, low_bound, high_bound, lvl_type)
                    monitored_lvl_df = self._select_levels(
                        monitored_df, low_bound, high_bound, lvl_type)

                    lats = {'assimilated':
                            assimilated_lvl_df['latitude'].to_numpy(),
                            'rejected':
                            rejected_lvl_df['latitude'].to_numpy(),
                            'monitored':
                            monitored_lvl_df['latitude'].to_numpy()
                            }
                    lons = {'assimilated':
                            assimilated_lvl_df['longitude'].to_numpy(),
                            'rejected':
                            rejected_lvl_df['longitude'].to_numpy(),
                            'monitored':
                            monitored_lvl_df['longitude'].to_numpy()
                            }

                    binned_lats[f'{low_bound}-{high_bound}'] = lats
                    binned_lons[f'{low_bound}-{high_bound}'] = lons

                return binned_lats, binned_lons

            else:
                lats = {'assimilated': assimilated_df['latitude'].to_numpy(),
                        'rejected': rejected_df['latitude'].to_numpy(),
                        'monitored': monitored_df['latitude'].to_numpy()
                        }
                lons = {'assimilated': assimilated_df['longitude'].to_numpy(),
                        'rejected': rejected_df['longitude'].to_numpy(),
                        'monitored': monitored_df['longitude'].to_numpy()
                        }

                return lats, lons

        else:
            indexed_df = self._select_conv(obsid, subtype, station_id)

            # select by levels
            if lvls is not None:
                binned_lats = {}
                binned_lons = {}

                for i, low_bound in enumerate(lvls[:-1]):
                    high_bound = lvls[i+1]

                    lvl_df = self._select_levels(
                        indexed_df, low_bound, high_bound, lvl_type)

                    lats = lvl_df['latitude'].to_numpy()
                    lons = lvl_df['longitude'].to_numpy()

                    binned_lats[f'{low_bound}-{high_bound}'] = lats
                    binned_lons[f'{low_bound}-{high_bound}'] = lons

                return binned_lats, binned_lons

            else:
                lats = indexed_df['latitude'].to_numpy()
                lons = indexed_df['longitude'].to_numpy()

                return lats, lons

	def return_df(self):
		"""Returns working dataframe."""
		return self.data_df


class Radiance(GSIdiag):

    def __init__(self, path):
        """
        Initialize a Radiance GSI diagnostic object

        Args:
            path : (str) path to radiance GSI diagnostic object
        Returns:
            self : GSI diag radiance object containing the path
                   to extract data
        """
        super().__init__(path)

        self._read_obs()
        self.metadata['Diag File Type'] = 'radiance'

    def __str__(self):
        return "Radiance GSI diagnostic object"

    def _read_obs(self):
        """
        Reads the data from the radiance diagnostic file during
        initialization into a multidimensional pandas dataframe.
        """
        df_dict = {}
        chan_info = {}

        with Dataset(self.path, mode='r') as f:

            # Grab dimensions to get lens
            nchans = f.dimensions['nchans']
            nobs = f.dimensions['nobs']

            for var in f.variables:
                if len(f.variables[var].shape) == 1:
                    # Add channel info to own dict
                    if len(f.variables[var][:]) == len(nchans):
                        chan_info[var] = f.variables[var][:]
                    elif len(f.variables[var][:]) == len(nobs):
                        df_dict[var] = f.variables[var][:]

        self.chan_info = chan_info

        # Sets correct channel number to indexed channel
        nchans = len(chan_info['chaninfoidx'])
        iters = int(len(df_dict['Channel_Index'])/nchans)

        for a in range(iters):
            df_dict['Channel_Index'][a*nchans:(a+1)*nchans] = \
                chan_info['sensor_chan']
        df_dict['Channel'] = df_dict['Channel_Index']

        # Create pandas dataframe from dict
        df = pd.DataFrame(df_dict)

        # Creates multidimensional indexed dataframe
        indices = ['Channel', 'QC_Flag']
        df.set_index(indices, inplace=True)

        # Rename columns
        df.columns = df.columns.str.lower()

        # Rename cloud_frac to cloud_fraction
        df = df.rename(columns={'cloud_frac': 'cloud_fraction'})

        for bias_type in ['unadjusted', 'adjusted']:
            df = df.rename(columns={
                f'obs_minus_forecast_{bias_type}': f'omf_{bias_type}',
                })
            # Create hofx columns
            df[f'hofx_{bias_type}'] = df['observation'] - \
                df[f'omf_{bias_type}']

        self.data_df = df

    def get_data(self, diag_type, channel=None, qcflag=None, use_flag=False,
                 separate_channels=False, separate_qc=False, errcheck=True,
                 bias_correction=True):
        """
        Given parameters, get the data from a radiance diagnostic file

        Args:
            diag_type : (str; Required) type of data to extract
                        i.e. observation, omf, oma, hofx
            channel : (list of ints; default=None) observation channel number
            qcflag : (list of ints; default=None) qc flag number
            separate_channels : (bool; default=False) if True, returns
                                dict of separate data by specified channels
            separate_qc : (bool; default=False) if True, returns dict of
                          separate data by specified qc flag
            use_flag : (bool; default=False) if True, will only return where
                       use_flag==1
            errcheck : (bool; default=True) when True and qcflag==0, will
                       toss out obs where inverse obs error is zero (i.e.
                       not assimilated in GSI)
            bias_correction : (bool; default=True) If True, will return bias
                              corrected data.
        Returns:
            data : requested indexed data
        """
        if diag_type not in _VALID_RADIANCE_DIAG_TYPES:
            raise ValueError((f'{diag_type} is not a valid diag_type. '
                              'Valid choices are: '
                              f'{" | ".join(_VALID_RADIANCE_DIAG_TYPES)}'))

        self.metadata['Variable'] = 'brightness_temperature'
        self.metadata['Diag Type'] = diag_type
        self.metadata['Anl Use'] = False

        # If no channels given, return all channels
        self.metadata['Channels'] = 'All Channels' if channel is None \
            else channel

        # If no qc flags given, return all qc flags
        self.metadata['QC Flag'] = 'All QC Flags' if qcflag is None \
            else qcflag

        if separate_channels or separate_qc:
            data = self._get_data_special(
                diag_type, channel, qcflag, use_flag, separate_channels,
                separate_qc, errcheck, bias_correction)
            return data

        else:
            indexed_df = self._select_radiance(
                channel, qcflag, use_flag, errcheck)
            data = self._query_diag_type(
                indexed_df, diag_type, bias_correction)

            data[data > 1e5] = np.nan

            return data

    def _select_radiance(self, channel=None, qcflag=None, use_flag=False,
                         errcheck=True):
        """
        Given parameters, get the indices of the observation
        locations from a radiance diagnostic file.
        """
        df = self.data_df

        if use_flag:
            use_flag_indx = np.where(self.chan_info['use_flag'] == 1)
            use_flag_channel = self.chan_info[
                'sensor_chan'][use_flag_indx].tolist()

            self.metadata['Channels'] = use_flag_channel

            idx_col = 'Channel'
            indx = df.index.get_level_values(idx_col) == ''

            for chan in use_flag_channel:
                indx = np.ma.logical_or(
                    indx, df.index.get_level_values(idx_col) == chan)

            df = df.iloc[indx]

        # index dataframe by channel
        if channel is not None:
            idx_col = 'Channel'
            indx = df.index.get_level_values(idx_col) == ''
            for chan in channel:
                indx = np.ma.logical_or(
                    indx, df.index.get_level_values(idx_col) == chan)

                # If channel not valid, raise TypeError
                if not any(indx):
                    VALIDCHANS = df.index.get_level_values(
                        'Channel').unique().to_numpy()
                    raise TypeError(f'Channel {chan} is not a valid channel. '
                                    'Valid channels include: '
                                    f'{", ".join(str(i) for i in VALIDCHANS)}')
            df = df.iloc[indx]

        # index dataframe by qcflag
        if qcflag is not None:
            idx_col = 'QC_Flag'
            indx = df.index.get_level_values(idx_col) == ''
            for qcf in qcflag:
                indx = np.ma.logical_or(
                    indx, df.index.get_level_values(idx_col) == qcf)

            df = df.iloc[indx]

            # remove obs where inverse obs error is zero
            if errcheck and 0 in qcflag:
                indx = df.index.get_level_values(idx_col) == ''

                # Grab index where inverse ob error is not zero
                err_indx = np.isin(
                    df['inverse_observation_error'], 0, invert=True)
                indx = np.ma.logical_or(indx, err_indx)

                df = df.iloc[indx]

        return df

    def _get_data_special(self, diag_type, channel, qcflag, use_flag,
                          separate_channels, separate_qc, errcheck,
                          bias_correction):
        """
        Creates a dictionary that separates channels and qc flags
        depending on the conditions of seperate_channels and
        separate_qc
        """
        data_dict = {}

        if separate_channels and not separate_qc:
            for c in channel:
                indexed_df = self._select_radiance(
                    [c], qcflag, errcheck=errcheck)
                data = self._query_diag_type(
                    indexed_df, diag_type, bias_correction)
                data[data > 1e5] = np.nan

                data_dict['Channel %s' % c] = data

        if not separate_channels and separate_qc:
            for qc in qcflag:
                indexed_df = self._select_radiance(
                    channel, [qc], errcheck=errcheck)
                data = self._query_diag_type(diag_type, idx)
                data[data > 1e5] = np.nan

                data_dict['QC Flag %s' % qc] = data

        if separate_channels and separate_qc:
            for c in channel:
                data_dict['Channel %s' % c] = {}
                for qc in qcflag:
                    indexed_df = self._select_radiance(
                        [c], [qc], errcheck=errcheck)
                    data = self._query_diag_type(
                        indexed_df, diag_type, bias_correction)
                    data[data > 1e5] = np.nan

                    data_dict['Channel %s' % c]['QC Flag %s' % qc] = data

        return data_dict

    def get_lat_lon(self, channel=None, qcflag=None, use_flag=False,
                    errcheck=True):
        """
        Gets lats and lons with desired indices.

        Args:
            channel : (list of ints; default=None) observation channel number
            qcflag : (list of ints; default=None) qc flag number
            use_flag : (bool; default=False) if True, will only return where
                       use_flag==1
            errcheck : (bool; default=True) when True and qcflag==0, will
                       toss out obs where inverse obs error is zero (i.e.
                       not assimilated in GSI)
        Returns:
            lat, lon : (array like) indexed latitude and longitude values
        """
        indexed_df = self._select_radiance(channel, qcflag, use_flag, errcheck)
        lats = indexed_df['latitude'].to_numpy()
        lons = indexed_df['longitude'].to_numpy()

        return lats, lons

	def return_df(self):
        """Returns working dataframe."""
        return self.data_df


class Ozone(GSIdiag):

    def __init__(self, path):
        """
        Initialize an Ozone GSI diagnostic object

        Args:
            path : (str) path to ozone GSI diagnostic object
        Returns:
            self : GSI diag ozone object containing the path
                   to extract data
        """
        super().__init__(path)

        self._read_obs()
        self.metadata['Diag File Type'] = 'ozone'

    def __str__(self):
        return "Ozone GSI diagnostic object"

    def _read_obs(self):
        """
        Reads the data from the ozone diagnostic file during
        initialization into a multidimensional pandas dataframe.
        """
        df_dict = {}

        with Dataset(self.path, mode='r') as f:
            for var in f.variables:
                df_dict[var] = f.variables[var][:]

        # Create pandas dataframe from dict
        df = pd.DataFrame(df_dict)
        indices = ['Reference_Pressure', 'Analysis_Use_Flag']
        df.set_index(indices, inplace=True)

        # Rename columns
        df.columns = df.columns.str.lower()
        for bias_type in ['unadjusted', 'adjusted']:
            df = df.rename(columns={
                f'obs_minus_forecast_{bias_type}': f'omf_{bias_type}',
                })
            # Create hofx columns
            df[f'hofx_{bias_type}'] = df['observation'] - \
                df[f'omf_{bias_type}']

        self.data_df = df

    def get_data(self, diag_type, analysis_use=False, errcheck=True,
                 bias_correction=True):
        """
        Given parameters, get the data from an ozone diagnostic file

        Args:
            diag_type : (str; Required) type of data to extract
                        i.e. observation, omf, oma, hofx
            analysis_use : (bool; default=False) if True, will return
                           two sets of data:
                           assimilated (analysis_use_flag=1),
                           monitored (analysis_use_flag=-1)
            errcheck : (bool; default=True) when True, will toss out
                       obs where inverse obs error is zero (i.e.
                       not assimilated in GSI)
        Returns:
            data_dict : (dict) requested indexed data
        """
        if diag_type not in _VALID_OZONE_DIAG_TYPES:
            raise ValueError((f'{diag_type} is not a valid diag_type. '
                              'Valid choices are: '
                              f'{" | ".join(_VALID_OZONE_DIAG_TYPES)}'))

        self.metadata['Variable'] = 'Ozone'
        self.metadata['Diag Type'] = diag_type
        self.metadata['Anl Use'] = analysis_use

        data_dict = {}

        pressures = self.data_df.index.get_level_values(
            'Reference_Pressure').unique().to_numpy()

        # Loop through all pressures. If pressure is 0, save in
        # data_dict as 'column total', else save as pressure level
        for p in pressures:
            if analysis_use:
                assimilated_df, monitored_df = self._select_ozone(
                    p, analysis_use, errcheck)

                assimilated_data = self._query_diag_type(
                    assimilated_df, diag_type, bias_correction)
                monitored_data = self._query_diag_type(
                    monitored_df, diag_type, bias_correction)

                if p == 0:
                    data_dict['column total'] = {
                        'assimilated': assimilated_data,
                        'monitored': monitored_data
                    }
                else:
                    data_dict[p] = {
                        'assimilated': assimilated_data,
                        'monitored': monitored_data
                    }

            else:
                indexed_df = self._select_ozone(
                    p, analysis_use, errcheck)

                data = self._query_diag_type(
                    indexed_df, diag_type, bias_correction)

                if p == 0:
                    data_dict['column total'] = data
                else:
                    data_dict[p] = data

        return data_dict

    def _select_ozone(self, pressure, analysis_use, errcheck):
        """
        Given parameters, multidimensional dataframe is indexed
        to only include selected locations from an ozone
        diagnostic file.

        Args:
            pressure : (float) pressure level of ozone data
            analysis_use : (bool) if True, will separate into two
                           indexed dataframes: assimilated, monitored
            errcheck : (bool) when True, will toss out obs where
                       inverse obs error is zero (i.e. not
                       assimilated in GSI)
        Returns:
            df : (pandas dataframe) indexed multidimentsional
                 dataframe from selected data
        """
        df = self.data_df

        idx_col = 'Reference_Pressure'
        indx = df.index.get_level_values(idx_col) == ''
        indx = np.ma.logical_or(
            indx, df.index.get_level_values(idx_col) == pressure)

        df = df[indx]

        if analysis_use:
            assimilated_indx = (
                df.index.get_level_values('Analysis_Use_Flag') == 1)
            monitored_indx = (
                df.index.get_level_values('Analysis_Use_Flag') == -1)

            assimilated_df = df.iloc[assimilated_indx]
            monitored_df = df.iloc[monitored_indx]

            if errcheck:
                assimilated_err_indx = np.isin(
                    assimilated_df['inverse_observation_error'],
                    0, invert=True)
                monitored_err_indx = np.isin(
                    monitored_df['inverse_observation_error'], 0, invert=True)

                assimilated_df = assimilated_df.iloc[assimilated_err_indx]
                monitored_df = monitored_df.iloc[monitored_err_indx]

            return assimilated_df, monitored_df

        else:
            if errcheck:
                err_indx = np.isin(
                    df['inverse_observation_error'], 0, invert=True)
                df = df.iloc[err_indx]

            return df

    def get_lat_lon(self, analysis_use=False, errcheck=True):
        """
        Gets lats and lons with desired indices. Returns dictionary for
        each level in the ozone data.
        """
        lats_dict = {}
        lons_dict = {}

        pressures = self.data_df.index.get_level_values(
            'Reference_Pressure').unique().to_numpy()

        # Loop through all pressures. If pressure is 0, save in
        # dicts as 'column total', else save as pressure level
        for p in pressures:
            if analysis_use:
                assimilated_df, monitored_df = self._select_ozone(
                    p, analysis_use, errcheck)

                if p == 0:
                    lats_dict['column total'] = {
                        'assimilated': assimilated_df['latitude'].to_numpy(),
                        'monitored': monitored_df['latitude'].to_numpy()
                        }
                    lons_dict['column total'] = {
                        'assimilated': assimilated_df['longitude'].to_numpy(),
                        'monitored': monitored_df['longitude'].to_numpy()
                        }
                else:
                    lats_dict[p] = {
                        'assimilated': assimilated_df['latitude'].to_numpy(),
                        'monitored': monitored_df['latitude'].to_numpy()
                    }
                    lons_dict[p] = {
                      'assimilated': assimilated_df['longitude'].to_numpy(),
                      'monitored': monitored_df['longitude'].to_numpy()
                    }

            else:
                indexed_df = self._select_ozone(
                    p, analysis_use, errcheck)

                if p == 0:
                    lats_dict['column total'] = \
                        indexed_df['latitude'].to_numpy()
                    lons_dict['column total'] = \
                        indexed_df['longitude'].to_numpy()
                else:
                    lats_dict[p] = indexed_df['latitude'].to_numpy()
                    lons_dict[p] = indexed_df['longitude'].to_numpy()

        return lats_dict, lons_dict

	def return_df(self):
        """Returns working dataframe."""
        return self.data_df

