import os
import pandas as pd
import numpy as np
from netCDF4 import Dataset
from datetime import datetime
from pathlib import Path

_VALID_LVL_TYPES = ["pressure", "height"]


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

        # Get index values
        self.obs_ids = df.index.get_level_values(
            'Observation_Type').unique().to_numpy()
        self.stn_ids = df.index.get_level_values(
            'Station_ID').unique().to_numpy()

        self.data_df = df

    def get_data(self, obsid=None, subtype=None, station_id=None,
                 analysis_use=False, lvls=None, lvl_type='pressure'):
        """
        Given parameters, get indexed dataframe from a conventional diagnostic
        file.

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
            indexed_df : requested indexed dataframe
        """

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

            binned_dfs = self._get_lvl_data(obsid, subtype, station_id,
                                            analysis_use, lvls, lvl_type)

            return binned_dfs

        else:
            if analysis_use:
                assimilated_df, rejected_df, monitored_df = self._select_conv(
                    obsid, subtype, station_id, analysis_use)

                indexed_df = {
                    'assimilated': assimilated_df,
                    'rejected': rejected_df,
                    'monitored': monitored_df
                }

            else:
                indexed_df = self._select_conv(obsid, subtype, station_id,
                                               analysis_use)

            return indexed_df

    def _select_conv(self, obsid, subtype, station_id, analysis_use):
        """
        Given parameters, multidimensional dataframe is indexed
        to only include selected locations from a conventional
        diagnostic file.

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

    def _get_lvl_data(self, obsid, subtype, station_id,
                      analysis_use, lvls, lvl_type):
        """
        Given a list of levels, will create a dictionary of data that is
        selected between each level. Will return a dictionary with subsetted
        pressure or height levels where data is separated within those levels:

        dict = {250-500: <data>,
                500-750: <data>,
                750-1000: <data>}
        """
        binned_dfs = {}

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

                lvl_df = {
                    'assimilated': assimilated_lvl_df,
                    'rejected': rejected_lvl_df,
                    'monitored': monitored_lvl_df
                }

                binned_dfs[f'{low_bound}-{high_bound}'] = lvl

            else:
                indexed_df = self._select_conv(obsid, subtype, station_id)

                lvl_df = self._select_levels(
                    indexed_df, low_bound, high_bound, lvl_type)

                binned_dfs[f'{low_bound}-{high_bound}'] = lvl_df

        return binned_dfs

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

    def list_obsids(self):
        """
        Prints all the unique observation ids in the diagnostic file.
        """
        obsids = sorted(self.obs_ids)

        print('Observation IDs in this diagnostic file include:\n'
              f'{", ".join([str(x) for x in obsids])}')

    def list_stationids(self):
        """
        Prints all the unique station ids in the diagnostic file.
        """
        stnids = sorted(self.stn_ids)

        print('Station IDs in this diagnostic file include:\n'
              f'{", ".join([str(x) for x in stnids])}')


class Radiance(GSIdiag):

    def __init__(self, path, var_names=None, read_jac=False):
        """
        Initialize a Radiance GSI diagnostic object

        Args:
            path      : (str) path to radiance GSI diagnostic object
            var_names : (string array, optional) Names of state
                        variables in Jacobians.
            read_jac  : (bool, optional) Option to read in Jacobians.
                        Default = False.
        Returns:
            self : GSI diag radiance object containing the path
                   to extract data
        """
        super().__init__(path)

        # The next two commands are for when the Jacobian is read in.
        # This is optional and defaults to false.  It also requires the
        # state variables to be defined as they are not present in the
        # diagnostic file.
        self.var_names = var_names or \
            ['sst', 'u', 'v', 'tv', 'q', 'oz', 'ql', 'qi']
        self.read_jac = read_jac
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
        jacobians_present = False
        jacstart = {}
        jacend = {}
        jac = {}

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
                elif len(f.variables[var].shape) == 2:
                    if self.read_jac:
                        if (var == 'Observation_Operator_Jacobian_stind'):
                            jacstart = f.variables[var][:, :]
                        if (var == 'Observation_Operator_Jacobian_endind'):
                            jacend = f.variables[var][:, :]
                        if (var == 'Observation_Operator_Jacobian_val'):
                            jac = f.variables[var][:, :]
                            jacobians_present = True

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

        # Get index values
        self.channels = df.index.get_level_values(
            'Channel').unique().to_numpy()
        self.qc_flags = df.index.get_level_values(
            'QC_Flag').unique().to_numpy()

        # Process Jacobians if present
        if jacobians_present:

            # The jacstart and jacend arrays from the netCDF file
            # actually refer to internal GSI indicies.  We can use them
            # to infer the length of each variable's section of the
            # Jacobian and set up indicies within the Jac array.

            # Find the length of each variable's Jacobian assuming all
            # channels/locs are the same (they are)
            len_jacs = np.array(jacend[0, :] - jacstart[0, :] + 1)
            # Make jac_indices dataframe with positions of each Jacoboian
            end_index = np.cumsum(len_jacs)-1
            start_index = np.roll(end_index, 1)+1
            start_index[0] = 0
            jac_indices = pd.DataFrame([start_index, end_index],
                                       columns=self.var_names)
            jacobians = {}
            for var in self.var_names:
                jac_range = jac_indices[var][1]-jac_indices[var][0]+1
                jacobians[var] = \
                    pd.DataFrame(jac[:, jac_indices[var][0]:
                                     jac_indices[var][1]+1],
                                 columns=np.arange(jac_range))
                jacobians[var].index = df_dict['Channel_Index']
                jacobians[var].insert(0, 'Latitude', df_dict['Latitude'])
                jacobians[var].insert(1, 'Longitude', df_dict['Longitude'])

            del jac

            self.jacobians = jacobians

        self.data_df = df

    def get_data(self, channel=None, qcflag=None, analysis_use=False,
                 use_flag=False, separate_channels=False, separate_qc=False,
                 errcheck=True):
        """
        Given parameters, get indexed dataframe from a radiance diagnostic file

        Args:
            channel : (list of ints; default=None) observation channel number
            qcflag : (list of ints; default=None) qc flag number
            analysis_use : (bool; default=False) if True, will return
                           three sets of data:
                           assimilated (QC_Flag=0, inv_observation_error!=0),
                           rejected (QC_Flag!=0),
                           monitored (use_flag!=1)
            separate_channels : (bool; default=False) if True, returns
                                dict of separate data by specified channels
            separate_qc : (bool; default=False) if True, returns dict of
                          separate data by specified qc flag
            use_flag : (bool; default=False) if True, will only return where
                       use_flag==1
            errcheck : (bool; default=True) when True and qcflag==0, will
                       toss out obs where inverse obs error is zero (i.e.
                       not assimilated in GSI)
        Returns:
            indexed_df : requested indexed dataframe
        """

        self.metadata['Variable'] = 'brightness_temperature'
        self.metadata['Anl Use'] = analysis_use

        # If no channels given, return all channels
        self.metadata['Channels'] = 'All Channels' if channel is None \
            else channel

        # If no qc flags given, return all qc flags
        self.metadata['QC Flag'] = 'All QC Flags' if qcflag is None \
            else qcflag

        if separate_channels or separate_qc:
            df_dict = self._get_data_special(
                channel, qcflag, use_flag, separate_channels,
                separate_qc, errcheck)

            return df_dict

        else:
            if analysis_use:
                assimilated_df, rejected_df, monitored_df = \
                    self._select_radiance(channel, qcflag,
                                          analysis_use, use_flag,
                                          errcheck)

                indexed_df = {
                    'assimilated': assimilated_df,
                    'rejected': rejected_df,
                    'monitored': monitored_df
                }

            else:
                indexed_df = self._select_radiance(
                    channel, qcflag, analysis_use, use_flag, errcheck)

            return indexed_df

    def _select_radiance(self, channel, qcflag, analysis_use,
                         use_flag, errcheck):
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

        if analysis_use:
            # Separate into 3 dataframes; assimilated, rejected, and monitored
            # First step, separate use_flag=1 and useflag!=1
            use_indx = np.where(self.chan_info['use_flag'] == 1)
            monitor_indx = np.where(self.chan_info['use_flag'] != 1)

            assimilated_channels = self.chan_info[
                'sensor_chan'][use_indx].tolist()
            monitored_channels = self.chan_info[
                'sensor_chan'][monitor_indx].tolist()

            # Create assimilated channels
            idx_col = 'Channel'
            indx = df.index.get_level_values(idx_col) == ''

            for chan in assimilated_channels:
                indx = np.ma.logical_or(
                    indx, df.index.get_level_values(idx_col) == chan)

            good_channels_df = df.iloc[indx]

            # reset index and get monitored dataframe
            indx = df.index.get_level_values(idx_col) == ''

            for chan in monitored_channels:
                indx = np.ma.logical_or(
                    indx, df.index.get_level_values(idx_col) == chan)

            monitored_df = df.iloc[indx]

            # Get assimilated data where QC Flag = 0 and inverse ob error != 0
            idx_col = 'QC_Flag'
            indx = good_channels_df.index.get_level_values(idx_col) == ''

            assimilated_indx = np.ma.logical_or(
                indx, good_channels_df.index.get_level_values('QC_Flag') == 0)
            err_indx = np.isin(
                good_channels_df['inverse_observation_error'], 0, invert=True)
            assimilated_indx = np.ma.logical_or(assimilated_indx, err_indx)

            assimilated_df = good_channels_df.iloc[assimilated_indx]

            # Get rejected data where QC != 1
            rejected_indx = np.ma.logical_or(
                indx, good_channels_df.index.get_level_values('QC_Flag') != 0)

            rejected_df = good_channels_df.iloc[rejected_indx]

            return assimilated_df, rejected_df, monitored_df

        else:
            return df

    def _get_data_special(self, channel, qcflag, use_flag,
                          separate_channels, separate_qc, errcheck):
        """
        Creates a dictionary that separates channels and qc flags
        depending on the conditions of seperate_channels and
        separate_qc
        """
        df_dict = {}

        if separate_channels and not separate_qc:
            for c in channel:
                indexed_df = self._select_radiance(
                    [c], qcflag, errcheck=errcheck)

                df_dict['Channel %s' % c] = indexed_df

        if not separate_channels and separate_qc:
            for qc in qcflag:
                indexed_df = self._select_radiance(
                    channel, [qc], errcheck=errcheck)

                df_dict['QC Flag %s' % qc] = indexed_df

        if separate_channels and separate_qc:
            for c in channel:
                df_dict['Channel %s' % c] = {}
                for qc in qcflag:
                    indexed_df = self._select_radiance(
                        [c], [qc], errcheck=errcheck)

                    df_dict['Channel %s' % c]['QC Flag %s' % qc] = indexed_df

        return df_dict

    def list_channels(self):
        """
        Prints all the unique channels in the diagnostic file.
        """
        chans = sorted(self.channels)

        print('Channel numbers in this diagnostic file include:\n'
              f'{", ".join([str(x) for x in chans])}')

    def list_qcflags(self):
        """
        Prints all the unique qcflags in the diagnostic file.
        """
        qcflags = sorted(self.qc_flags)

        print('QC Flags in this diagnostic file include:\n'
              f'{", ".join([str(x) for x in qcflags])}')


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

        # Get index values
        self.pressures = df.index.get_level_values(
            'Reference_Pressure').unique().to_numpy()

        self.data_df = df

    def get_data(self, analysis_use=False, errcheck=True):
        """
        Given parameters, get indexed dataframes from an ozone diagnostic file

        Args:
            analysis_use : (bool; default=False) if True, will return
                           two sets of data:
                           assimilated (analysis_use_flag=1),
                           monitored (analysis_use_flag=-1)
            errcheck : (bool; default=True) when True, will toss out
                       obs where inverse obs error is zero (i.e.
                       not assimilated in GSI)
        Returns:
            df_dict : (dict) requested indexed dataframes
        """

        self.metadata['Variable'] = 'Ozone'
        self.metadata['Anl Use'] = analysis_use

        df_dict = {}

        pressures = self.pressures

        # Loop through all pressures. If pressure is 0, save in
        # df_dict as 'column total', else save as pressure level
        for p in pressures:
            if analysis_use:
                assimilated_df, monitored_df = self._select_ozone(
                    p, analysis_use, errcheck)

                if p == 0:
                    df_dict['column total'] = {
                        'assimilated': assimilated_df,
                        'monitored': monitored_df
                    }
                else:
                    df_dict[p] = {
                        'assimilated': assimilated_df,
                        'monitored': monitored_df
                    }

            else:
                indexed_df = self._select_ozone(
                    p, analysis_use, errcheck)

                if p == 0:
                    df_dict['column total'] = indexed_df
                else:
                    df_dict[p] = indexed_df

        return df_dict

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

    def list_pressures(self):
        """
        Prints all the unique pressure levels in the diagnostic file.
        """
        plvls = sorted(self.pressures)
        plvls = ["Column Total" if x == 0 else str(x) for x in plvls]

        print('Pressure levels in this diagnostic file include:\n'
              f'{", ".join([x for x in plvls])}')
