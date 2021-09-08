# 2021-09-02 X.Su : modified _read_obs such as setupqc, vqc, and all error items
# 2021-09-02 X.Su : Added one data category: rejected to the functions: get_data, 
#                   _get_idx_conv, get_lat_lon, get_pressure,get_height for 
#                   conventional data. 
import os
import numpy as np
from netCDF4 import Dataset
from datetime import datetime
from pathlib import Path

_VALID_LVL_TYPES = ["pressure", "height"]
_VALID_DIAG_TYPES = ["omf", "o-f", "omb", "o-b", "oma", "o-a", "observation", "hofx"]

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

    def query_diag_type(self, diag_type, idx):
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
                u = self.u_omf[idx]
                v = self.v_omf[idx]

                return u, v

            else:
                data = self.omf[idx]
                return data

        if oma and is_anl:
            if self.variable == 'uv':
                u = self.u_omf[idx]
                v = self.v_omf[idx]

                return u, v

            else:
                data = self.omf[idx]
                return data

        if obs:
            if self.variable == 'uv':
                u = self.u_o[idx]
                v = self.v_o[idx]

                return u, v

            else:
                data = self.o[idx]
                return data

        if hofx:
            if self.variable == 'uv':
                u = self.u_o[idx] - self.u_omf[idx]
                v = self.v_o[idx] - self.v_omf[idx]

                return u, v

            else:
                hx = self.o - self.omf
                data = hx[idx]
                return data

        if water_fraction:
            data = self.water_frac[idx]
            return data

        if land_fraction:
            data = self.land_frac[idx]
            return data

        if ice_fraction:
            data = self.ice_frac[idx]
            return data

        if snow_fraction:
            data = self.snow_frac[idx]
            return data

        if cloud_fraction:
            data = self.cloud_frac[idx]
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
        """

        with Dataset(self.path, mode='r') as f:
            self.o_type = f.variables['Observation_Type'][:]
            try:
                self.o_stype = f.variables['Observation_Subtype'][:]
            except:
                print('Subtype is not defined')           

            self.lons = f.variables['Longitude'][:]
            self.lats = f.variables['Latitude'][:]
            self.press = f.variables['Pressure'][:]
            self.height = f.variables['Height'][:]
            self.time = f.variables['Time'][:]
            self.anl_use = f.variables['Analysis_Use_Flag'][:]
            self.stnid = f.variables['Station_ID'][:]
            try:
                self.prepqc = f.variables['Prep_QC_Mark'][:]
            except:
                self.prepqc = f.variables['Setup_QC_Mark'][:]
            self.setupqc = f.variables['Prep_Use_Flag'][:]
            self.vqc = f.variables['Nonlinear_QC_Rel_Wgt'][:]
            self.input_err = f.variables['Errinv_Input'][:]
            self.adjst_err = f.variables['Errinv_Adjust'][:]
            self.finl_err = f.variables['Errinv_Final'][:]

            try:
                self.stnelev = f.variables['Station_Elevation'][:]
            except:
                self.modelelev = f.variables['Model_Elevation'][:]

            if self.variable == 'uv':
                self.u_o = f.variables['u_Observation'][:]
                self.v_o = f.variables['v_Observation'][:]

                self.u_omf = f.variables['u_Obs_Minus_Forecast_adjusted'][:]
                self.v_omf = f.variables['v_Obs_Minus_Forecast_adjusted'][:]

            else:
                self.o = f.variables['Observation'][:]
                self.omf = f.variables['Obs_Minus_Forecast_adjusted'][:]



    def get_data(self, diag_type, obsid=None, subtype=None, station_id=None, analysis_use=False, lvls=None, lvl_type='pressure'):
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
                qcflag       : used:qcflag<7 and analysis_use_flag=1, rejected: analysis_use_flag=-1;monitored: qcflag>7
                vqcflag      : used:vqcflag>=1;rejected by vqc:  vqcflag<1
                lvls         : List of pressure or height levels i.e. [250,500,750,1000]. List 
			       must be arranged low to high.  Will return a dictionary with
                               subsetted pressure or height levels where data is seperated 
                               within those levels:

                               dict = {250-500: <data>,
                                       500-750: <data>,
                                       750-1000: <data>}
                lvl_type     : lvls definition as 'pressure' or 'height'.  Default is 'pressure'.

        OUTPUT:
            data   : requested data
        """
       
        if diag_type not in _VALID_DIAG_TYPES:
            raise ValueError(f'{diag_type} wrong. Valid choices are: {" | ".join(_VALID_DIAG_TYPES)}') 

        self.metadata['Diag Type'] = diag_type
        self.metadata['ObsID'] = obsid
        self.metadata['Subtype'] = subtype
        self.metadata['Station ID'] = station_id
        self.metadata['Anl Use'] = analysis_use

        if lvls is not None:
       
            if lvl_type not in _VALID_LVL_TYPES:
                raise ValueError('{lvl_type} wrong, use "pressure" or "height" for input lvl_type'.format(lvl_type=repr(lvl_type)))

            level_list = lvls
            binned_data = {}

            if analysis_use:
                assimilated_idx, rejected_idx, monitored_idx = self._get_idx_conv(
                    obsid, subtype, station_id, analysis_use)


                for i, low_bound in enumerate(level_list[:-1]):
                    if lvl_type == 'height':
                       hght_idx = np.where((self.height >= level_list[i]) & (
                           self.height < level_list[i+1]))
                       valid_assimilated_idx = np.isin(
                           assimilated_idx[0], hght_idx[0])
                       valid_rejected_idx = np.isin(
                           rejected_idx[0], hght_idx[0])
                       valid_monitored_idx = np.isin(
                           monitored_idx[0], hght_idx[0])
                    else:
                       pres_idx = np.where((self.press > level_list[i]) & (
                           self.press <= level_list[i+1]))
                       valid_assimilated_idx = np.isin(
                           assimilated_idx[0], pres_idx[0])
                       valid_rejected_idx = np.isin(
                           rejected_idx[0], pres_idx[0])
                       valid_monitored_idx = np.isin(
                           monitored_idx[0], pres_idx[0])

                    assimilated_pidx = np.where(valid_assimilated_idx)
                    rejected_pidx = np.where(valid_rejected_idx)
                    monitored_pidx = np.where(valid_monitored_idx)

                    if self.variable == 'uv':
                        u_assimilated, v_assimilated = self.query_diag_type(
                            diag_type, assimilated_pidx)
                        u_rejected, v_rejected = self.query_diag_type(
                            diag_type, rejected_pidx)
                        u_monitored, v_monitored = self.query_diag_type(
                            diag_type, monitored_pidx)

                        data = {'u': {'assimilated': u_assimilated,
                                      'rejected': u_rejected,
                                      'monitored': u_monitored},
                                'v': {'assimilated': v_assimilated,
                                      'rejected': v_rejected,
                                      'monitored': v_monitored}
                                }

                        binned_data['%s-%s' %
                                        (level_list[i], level_list[i+1])] = data

                    else:
                        assimilated_data = self.query_diag_type(
                            diag_type, assimilated_pidx)
                        rejected_data = self.query_diag_type(
                            diag_type, rejected_pidx)
                        monitored_data = self.query_diag_type(
                            diag_type, monitored_pidx)

                        data = {'assimilated': assimilated_data,
                                'rejected': rejected_data,
                                'monitored': monitored_data
                                }

                        binned_data['%s-%s' %
                                        (level_list[i], level_list[i+1])] = data

                return binned_data

            else:
                idx = self._get_idx_conv(
                    obsid, subtype, station_id, analysis_use)

                for i, low_bound in enumerate(level_list[:-1]):

                    if lvl_type == 'height':
                       hght_idx = np.where((self.height >= level_list[i]) & (
                           self.height < level_list[i+1]))
                       valid_idx = np.isin(idx[0], hght_idx[0])
                       pidx = np.where(valid_idx)

                    else:
                       pres_idx = np.where((self.press > level_list[i]) & (
                           self.press <= level_list[i+1]))
                       valid_idx = np.isin(idx[0], pres_idx[0])
                       pidx = np.where(valid_idx)


                    if self.variable == 'uv':
                        u, v = self.query_diag_type(diag_type, pidx)

                        data = {'u': u,
                                'v': v}

                        binned_data['%s-%s' %
                                        (level_list[i], level_list[i+1])] = data

                    else:
                        data = self.query_diag_type(diag_type, pidx)

                        binned_data['%s-%s' %
                                        (level_list[i], level_list[i+1])] = data

                return binned_data

        else:

            if analysis_use:
                assimilated_idx, rejected_idx, monitored_idx = self._get_idx_conv(
                    obsid, subtype, station_id, analysis_use)

                if self.variable == 'uv':
                    u_assimilated, v_assimilated = self.query_diag_type(
                        diag_type, assimilated_idx)
                    u_rejected, v_rejected = self.query_diag_type(
                        diag_type, rejected_idx)
                    u_monitored, v_monitored = self.query_diag_type(
                        diag_type, monitored_idx)

                    u = {'assimilated': u_assimilated,
                         'rejected': u_rejected,
                         'monitored': u_monitored}
                    v = {'assimilated': v_assimilated,
                         'rejected': v_rejected,
                         'monitored': v_monitored}

                    return u, v
                else:
                    assimilated_data = self.query_diag_type(
                        diag_type, assimilated_idx)
                    rejected_data = self.query_diag_type(
                        diag_type, rejected_idx)
                    monitored_data = self.query_diag_type(
                        diag_type, monitored_idx)

                    data = {'assimilated': assimilated_data,
                            'rejected': rejected_data,
                            'monitored': monitored_data
                            }

                    return data

            else:
                idx = self._get_idx_conv(
                    obsid, subtype, station_id, analysis_use)

                if self.variable == 'uv':
                    u, v = self.query_diag_type(diag_type, idx)

                    return u, v

                else:
                    data = self.query_diag_type(diag_type, idx)

                    return data

    def _get_idx_conv(self, obsid=None, subtype=None, station_id=None, analysis_use=False):
        """
        Given parameters, get the indices of the observation
        locations from a conventional diagnostic file
        INPUT:
            obsid   : observation measurement ID number
            qcflag  : qc flag (default: None) i.e. <7, >7
            vqcflag : vqc flag (default: None) i.e. <1, >=1
            subtype : subtype number (default: None)
            station_id : station id tag (default: None)
        OUTPUT:
            idx    : indices of the requested data in the file
        """

        if not obsid or obsid == [None]:
            obsid = None
        if not subtype or subtype == [None]:
            subtype = None
        if not station_id or station_id == [None]:
            station_id = None

        if not analysis_use:
            idxobs = self.o_type
            valid_idx = np.full_like(idxobs, True, dtype=bool)
            if obsid != None:
                valid_obs_idx = np.isin(idxobs, obsid)
                valid_idx = np.logical_and(valid_idx, valid_obs_idx)
            if subtype != None:
                subtypeidx = self.o_stype
                valid_subtype_idx = np.isin(subtypeidx, subtype)
                valid_idx = np.logical_and(valid_idx, valid_subtype_idx)
            if station_id != None:
                stnididx = self.stnid
                stnididx = [i.tobytes(fill_value='/////', order='C') for i in stnididx]
                stnididx = np.array([''.join(i.decode('UTF-8', 'ignore').split()) for i in stnididx])
                valid_stnid_idx = np.isin(stnididx, station_id)
                valid_idx = np.logical_and(valid_idx, valid_stnid_idx)

            idx = np.where(valid_idx)

            return idx

        else:
            valid_assimilated_idx = np.isin(self.anl_use, 1)
            valid_rejected_idx = np.isin(self.anl_use, -1)
            valid_monitored_idx = np.isin(self.anl_use, -1)

            valid_qc_idx = self.prepqc < 8.0
            valid_qc2_idx = self.prepqc > 7.0
            valid_rejected_idx = np.logical_and(valid_rejected_idx, valid_qc_idx)
            valid_monitored_idx = np.logical_and(valid_monitored_idx, valid_qc2_idx)

            if obsid != None:
                idxobs = self.o_type
                valid_obs_idx = np.isin(idxobs, obsid)

                valid_assimilated_idx = np.logical_and(
                    valid_assimilated_idx, valid_obs_idx)
                valid_rejected_idx = np.logical_and(
                    valid_rejected_idx, valid_obs_idx)
                valid_monitored_idx = np.logical_and(
                    valid_monitored_idx, valid_obs_idx)

            if subtype != None:
                subtypeidx = self.o_stype
                valid_subtype_idx = np.isin(subtypeidx, subtype)

                valid_assimilated_idx = np.logical_and(
                    valid_assimilated_idx, valid_subtype_idx)
                valid_rejected_idx = np.logical_and(
                    valid_rejected_idx, valid_subtype_idx)
                valid_monitored_idx = np.logical_and(
                    valid_monitored_idx, valid_subtype_idx)

            if station_id != None:
                stnididx = self.stnid
                stnididx = np.array(
                    [''.join(i.tostring().decode().split()) for i in stnididx])
                valid_stnid_idx = np.isin(stnididx, station_id)

                valid_assimilated_idx = np.logical_and(
                    valid_assimilated_idx, valid_stnid_idx)
                valid_rejected_idx = np.logical_and(
                    valid_rejected_idx, valid_stnid_idx)
                valid_monitored_idx = np.logical_and(
                    valid_monitored_idx, valid_stnid_idx)

            assimilated_idx = np.where(valid_assimilated_idx)
            rejected_idx = np.where(valid_rejected_idx)
            monitored_idx = np.where(valid_monitored_idx)

            return assimilated_idx, rejected_idx, monitored_idx

    def get_lat_lon(self, obsid=None, subtype=None, station_id=None, analysis_use=False, lvls=None, lvl_type='pressure'):
        """
        Gets lats and lons with desired indices
        """

        if lvls is not None:

            if lvl_type not in _VALID_LVL_TYPES:
                raise ValueError('{lvl_type} wrong, use "pressure" or "height" for input lvl_type'.format(lvl_type=repr(lvl_type)))

            level_list = lvls
            binned_lats = {}
            binned_lons = {}

            if analysis_use:
                assimilated_idx, rejected_idx, monitored_idx = self._get_idx_conv(
                    obsid, subtype, station_id, analysis_use)

                for i, low_bound in enumerate(level_list[:-1]):
                    if lvl_type == 'height':
                       hght_idx = np.where((self.height >= level_list[i]) & (
                           self.height < level_list[i+1]))
                       valid_assimilated_idx = np.isin(
                           assimilated_idx[0], hght_idx[0])
                       valid_rejected_idx = np.isin(
                           rejected_idx[0], hght_idx[0])
                       valid_monitored_idx = np.isin(
                           monitored_idx[0], hght_idx[0])
                    else:
                       pres_idx = np.where((self.press > level_list[i]) & (
                           self.press <= level_list[i+1]))
                       valid_assimilated_idx = np.isin(
                           assimilated_idx[0], pres_idx[0])
                       valid_rejected_idx = np.isin(
                           rejected_idx[0], pres_idx[0]) 
                       valid_monitored_idx = np.isin(
                           monitored_idx[0], pres_idx[0])

                    assimilated_pidx = np.where(valid_assimilated_idx)
                    rejected_pidx = np.where(valid_rejected_idx)
                    monitored_pidx = np.where(valid_monitored_idx)

                    lats = {'assimilated': self.lats[assimilated_pidx],
                            'rejected': self.lats[rejected_pidx],
                            'monitored': self.lats[monitored_pidx]}
                    lons = {'assimilated': self.lons[assimilated_pidx],
                            'rejected': self.lons[rejected_pidx],
                            'monitored': self.lons[monitored_pidx]}

                    binned_lats[low_bound] = lats
                    binned_lons[low_bound] = lons

                return binned_lats, binned_lons

            else:
                binned_lats = {}
                binned_lons = {}

                idx = self._get_idx_conv(
                    obsid, subtype, station_id, analysis_use)

                for i, low_bound in enumerate(level_list[:-1]):
                    if lvl_type == 'height':
                        hght_idx = np.where((self.height >= level_list[i]) & (
                            self.height < level_list[i+1]))
                        valid_idx = np.isin(idx[0], hght_idx[0])
                        pidx = np.where(valid_idx)
                    else:
                        pres_idx = np.where((self.press > level_list[i]) & (
                            self.press <= level_list[i+1]))
                        valid_idx = np.isin(idx[0], pres_idx[0])
                        pidx = np.where(valid_idx)

                    binned_lats[low_bound] = self.lats[pidx]
                    binned_lons[low_bound] = self.lons[pidx]

                return binned_lats, binned_lons

        else:
            if analysis_use:
                assimilated_idx, rejected_idx, monitored_idx = self._get_idx_conv(
                    obsid, subtype, station_id, analysis_use)
                lats = {'assimilated': self.lats[assimilated_idx],
                        'rejected': self.lats[rejected_idx],
                        'monitored': self.lats[monitored_idx]}
                lons = {'assimilated': self.lons[assimilated_idx],
                        'rejected': self.lons[rejected_idx],
                        'monitored': self.lons[monitored_idx]}
                return lats, lons
            else:
                idx = self._get_idx_conv(
                    obsid, subtype, station_id, analysis_use)
                return self.lats[idx], self.lons[idx]

    def get_pressure(self, obsid=None, subtype=None, station_id=None, analysis_use=False):
        if analysis_use:
            assimilated_idx, rejected_idx, monitored_idx = self._get_idx_conv(
                obsid, subtype, station_id, analysis_use)
            pressure = {'assimilated': self.press[assimilated_idx],
                        'rejected': self.press[rejected_idx],
                        'monitored': self.press[monitored_idx]}

            return pressure
        else:
            idx = self._get_idx_conv(obsid, subtype, station_id, analysis_use)
            return self.press[idx]

    def get_height(self, obsid=None, subtype=None, station_id=None, analysis_use=False):
        if analysis_use:
            assimilated_idx, rejected_idx, monitored_idx = self._get_idx_conv(
                obsid, subtype, station_id, analysis_use)
            height = {'assimilated': self.height[assimilated_idx],
                      'rejected': self.height[rejected_idx],
                      'monitored': self.height[monitored_idx]}

            return height
        else:
            idx = self._get_idx_conv(obsid, subtype, station_id, analysis_use)
            return self.height[idx]


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
