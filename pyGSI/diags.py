import numpy as np
from netCDF4 import Dataset
from datetime import datetime


class GSIdiag:

    def __init__(self, path):
        """
        Initialize a GSI diagnostic object
        INPUT:
            path : path to GSI diagnostic object
        """

        self.path = path
        self.variable = path.split('/')[-1].split('.')[0].split('_')[2]

    def __len__(self):
        return len(self.lats)

    def get_metadata(self):
        """
        Grabs metadata from the diagnostic filename
        """

        dtype = self.path.split('/')[-1].split('.')[0].split('_')[1]
        str_date = self.path.split('/')[-1].split('.')[1]
        date = datetime.strptime(str_date, '%Y%m%d%H')
        ftype = self.path.split('/')[-1].split('.')[0].split('_')[-1]

        if dtype == 'conv':
            metadata = {'Diag_type': dtype,
                        'Variable': self.variable,
                        'Date': date,
                        'File_type': ftype
                        }
        else:
            satellite = self.variable

            metadata = {'Diag_type': dtype,
                        'Satellite': satellite,
                        'Date': date,
                        'File_type': ftype
                        }

        return metadata

    def query_data_type(self, dtype, idx):
        """
        Query the data type being requested and returns
        the appropriate indexed data
        """
        if dtype == 'O-F' or dtype == 'o-f':
            if self.variable == 'uv':
                u = self.u_omf[idx]
                v = self.v_omf[idx]

                return u, v

            else:
                data = self.omf[idx]
                return data

        elif dtype == 'Observation' or dtype == 'observation':
            if self.variable == 'uv':
                u = self.u_o[idx]
                v = self.v_o[idx]

                return u, v

            else:
                data = self.o[idx]
                return data

        elif dtype == 'O-A' or dtype == 'o-a':
            check = self.check_o_minus_a()
            if check:
                data = self.omf[idx]
                return data
            else:
                print('File type does not support O-A')
                return None

        elif dtype == 'H(x)' or dtype == 'h(x)' or dtype == 'hofx':
            if self.variable == 'uv':
                u = self.u_o[idx] - self.u_omf[idx]
                v = self.v_o[idx] - self.v_omf[idx]
                
                return u,v
            
            else:
                hx = self.o - self.omf
                data = hx[idx]
                return data

        elif dtype == 'Water_Fraction' or dtype == 'water_fraction':
            data = self.water_frac[idx]
            return data

        elif dtype == 'Land_Fraction' or dtype == 'land_fraction':
            data = self.land_frac[idx]
            return data

        elif dtype == 'Ice_Fraction' or dtype == 'ice_fraction':
            data = self.ice_frac[idx]
            return data

        elif dtype == 'Snow_Fraction' or dtype == 'snow_fraction':
            data = self.snow_frac[idx]
            return data

        elif dtype == 'Cloud_Fraction' or dtype == 'cloud_fraction':
            data = self.cloud_frac[idx]
            return data

        else:
            raise Exception(f'Unrecognizable data type: {dtype}')
            return None

    def check_o_minus_a(self):
        """
        Checks if the diagnostic file is an analysis file.
        """

        ftype = self.path.split('/')[-1].split('.')[0].split('_')[-1]

        if ftype == 'anl':
            return True
        else:
            return False


class Conventional(GSIdiag):

    def __init__(self, path):
        """
        Initialize a conventional GSI diagnostic object
        INPUT:
            path   : path to conventional GSI diagnostic object
        RESULT:
            self   : GSIdiag object containing the path to extract data
        """
        super().__init__(path)

        self.read_conv_obs()
        self.variable

    def __str__(self):
        return "Conventional object"

    def read_conv_obs(self):
        """
        Reads the data from the conventional diagnostic file during initialization.
        """

        f = Dataset(self.path, mode='r')

        self.o_type = f.variables['Observation_Type'][:]
        self.o_stype = f.variables['Observation_Subtype'][:]
        self.lons = f.variables['Longitude'][:]
        self.lats = f.variables['Latitude'][:]
        self.press = f.variables['Pressure'][:]
        self.time = f.variables['Time'][:]
        self.anl_use = f.variables['Analysis_Use_Flag'][:]
        try:
            self.stnelev = f.variables['Station_Elevation'][:]
        except:
            self.modelelev = f.variables['Model_Elevation'][:]

        if self.path.split('/')[-1].split('.')[0].split('_')[2] == 'uv':
            self.u_o = f.variables['u_Observation'][:]
            self.v_o = f.variables['v_Observation'][:]

            self.u_omf = f.variables['u_Obs_Minus_Forecast_adjusted'][:]
            self.v_omf = f.variables['v_Obs_Minus_Forecast_adjusted'][:]

        else:
            self.o = f.variables['Observation'][:]
            self.omf = f.variables['Obs_Minus_Forecast_adjusted'][:]

        f.close()

    def get_data(self, dtype, obsid=None, subtype=None, analysis_use=False, plvls=False):
        """
        Given parameters, get the data from a conventional diagnostic file
        INPUT:
            required:
                dtype  : type of data to extract i.e. Observation, O-F, O-A, H(x)

            optional:    
                obsid        : observation measurement ID number; default=None
                subtype      : observation measurement ID subtype number, default=None
                analysis_use : if True, will return two sets of data: assimlated
                               (analysis_use_flag=1), and monitored (analysis_use
                               _flag=-1); default = False
                plvls        : if True, will return a dictionary of data subsetting 
                               into pressure levels

        OUTPUT:
            data   : requested data
        """

        if plvls == True:
            pressure_list = [0, 100, 250, 500, 700, 850, 925, 1000, 1100]
            binned_pressure = {}

            if analysis_use == True:
                assimilated_idx, monitored_idx = self.get_idx_conv(
                    obsid, subtype, analysis_use)

                for i, pressure in enumerate(pressure_list[:-1]):
                    pres_idx = np.where((self.press > pressure_list[i]) & (
                        self.press < pressure_list[i+1]))
                    valid_assimilated_idx = np.isin(
                        assimilated_idx[0], pres_idx[0])
                    valid_monitored_idx = np.isin(
                        monitored_idx[0], pres_idx[0])

                    assimilated_pidx = np.where(valid_assimilated_idx)
                    monitored_pidx = np.where(valid_monitored_idx)

                    if self.path.split('/')[-1].split('.')[0].split('_')[2] == 'uv':
                        u_assimilated, v_assimilated = self.query_data_type(
                            dtype, assimilated_pidx)
                        u_monitored, v_monitored = self.query_data_type(
                            dtype, monitored_pidx)

                        u = {'assimilated': u_assimilated,
                             'monitored': u_monitored}
                        v = {'assimilated': v_assimilated,
                             'monitored': v_monitored}

                        binned_pressure[pressure]['u'] = u
                        binned_pressure[pressure]['v'] = v

                    else:
                        assimilated_data = self.query_data_type(
                            dtype, assimilated_pidx)
                        monitored_data = self.query_data_type(
                            dtype, monitored_pidx)

                        data = {'assimilated': assimilated_data,
                                'monitored': monitored_data
                                }

                        binned_pressure[pressure] = data

                return binned_pressure

            else:
                idx = self.get_idx_conv(obsid, subtype, analysis_use)

                for i, pressure in enumerate(pressure_list[:-1]):
                    pres_idx = np.where((self.press > pressure_list[i]) & (
                        self.press < pressure_list[i+1]))
                    valid_idx = np.isin(idx[0], pres_idx[0])
                    pidx = np.where(valid_idx)

                    if self.path.split('/')[-1].split('.')[0].split('_')[2] == 'uv':
                        u, v = self.query_data_type(dtype, pidx)
                        binned_pressure[pressure]['u'] = u
                        binned_pressure[pressure]['v'] = v
                    else:
                        data = self.query_data_type(dtype, pidx)
                        binned_pressure[pressure] = data

                return binned_pressure

        else:

            if analysis_use == True:
                assimilated_idx, monitored_idx = self.get_idx_conv(
                    obsid, subtype, analysis_use)

                if self.path.split('/')[-1].split('.')[0].split('_')[2] == 'uv':
                    u_assimilated, v_assimilated = self.query_data_type(
                        dtype, assimilated_idx)
                    u_monitored, v_monitored = self.query_data_type(
                        dtype, monitored_idx)

                    u = {'assimilated': u_assimilated,
                         'monitored': u_monitored}
                    v = {'assimilated': v_assimilated,
                         'monitored': v_monitored}

                    return u, v
                else:
                    assimilated_data = self.query_data_type(
                        dtype, assimilated_idx)
                    monitored_data = self.query_data_type(dtype, monitored_idx)

                    data = {'assimilated': assimilated_data,
                            'monitored': monitored_data
                            }

                    return data

            else:
                idx = self.get_idx_conv(obsid, subtype, analysis_use)

                if self.path.split('/')[-1].split('.')[0].split('_')[2] == 'uv':
                    u, v = self.query_data_type(dtype, idx)

                    return u, v

                else:
                    data = self.query_data_type(dtype, idx)

                    return data

    def get_idx_conv(self, obsid=None, subtype=None, analysis_use=False):
        """
        Given parameters, get the indices of the observation
        locations from a conventional diagnostic file
        INPUT:
            obsid   : observation measurement ID number
            qcflag  : qc flag (default: None) i.e. 0, 1
            subtype : subtype number (default: None)
        OUTPUT:
            idx    : indices of the requested data in the file
        """

        if not obsid or obsid == [None]:
            obsid = None
        if not subtype or subtype == [None]:
            subtype = None

        if analysis_use == False:
            idxobs = self.o_type
            valid_idx = np.full_like(idxobs, True, dtype=bool)
            if obsid != None:
                valid_obs_idx = np.isin(idxobs, obsid)
                valid_idx = np.logical_and(valid_idx, valid_obs_idx)
            if subtype != None:
                subtypeidx = self.o_stype
                valid_idx_subtype = np.isin(subtypeidx, subtype)
                valid_idx = np.logical_and(valid_idx, valid_idx_subtype)

            idx = np.where(valid_idx)

            return idx

        else:
            valid_assimilated_idx = np.isin(self.anl_use, 1)
            valid_monitored_idx = np.isin(self.anl_use, -1)

            if obsid != None:
                idxobs = self.o_type
                valid_obs_idx = np.isin(idxobs, obsid)

                valid_assimilated_idx = np.logical_and(
                    valid_assimilated_idx, valid_obs_idx)
                valid_monitored_idx = np.logical_and(
                    valid_monitored_idx, valid_obs_idx)

            if subtype != None:
                subtypeidx = self.o_stype
                valid_subtype_idx = np.isin(subtypeidx, subtype)

                valid_assimilated_idx = np.logical_and(
                    valid_assimilated_idx, valid_subtype_idx)
                valid_monitored_idx = np.logical_and(
                    valid_monitored_idx, valid_subtype_idx)

            assimilated_idx = np.where(valid_assimilated_idx)
            monitored_idx = np.where(valid_monitored_idx)

            return assimilated_idx, monitored_idx

    def get_lat_lon(self, obsid=None, subtype=None, analysis_use=False, plvls=False):
        """
        Gets lats and lons with desired indices
        """

        if plvls == True:
            pressure_list = [0, 100, 250, 500, 700, 850, 925, 1000, 1100]
            pressure_lats = {}
            pressure_lons = {}

            if analysis_use == True:
                assimilated_idx, monitored_idx = self.get_idx_conv(
                    obsid, subtype, analysis_use)

                for i, pressure in enumerate(pressure_list[:-1]):
                    pres_idx = np.where((self.press > pressure_list[i]) & (
                        self.press < pressure_list[i+1]))
                    valid_assimilated_idx = np.isin(
                        assimilated_idx[0], pres_idx[0])
                    valid_monitored_idx = np.isin(
                        monitored_idx[0], pres_idx[0])

                    assimilated_pidx = np.where(valid_assimilated_idx)
                    monitored_pidx = np.where(valid_monitored_idx)

                    lats = {'assimilated': self.lats[assimilated_pidx],
                            'monitored': self.lats[monitored_pidx]}
                    lons = {'assimilated': self.lons[assimilated_pidx],
                            'monitored': self.lons[monitored_pidx]}

                    pressure_lats[pressure] = lats
                    pressure_lons[pressure] = lons

                return pressure_lats, pressure_lons

            else:
                idx = self.get_idx_conv(obsid, subtype, analysis_use)

                for i, pressure in enumerate(pressure_list[:-1]):
                    pres_idx = np.where((self.press > pressure_list[i]) & (
                        self.press < pressure_list[i+1]))
                    valid_idx = np.isin(idx[0], pres_idx[0])
                    pidx = np.where(valid_idx)

                    pressure_lats[pressure] = self.lats[pidx]
                    pressure_lons[pressure] = self.lons[pidx]

                return pressure_lats, pressure_lons

        else:
            if analysis_use == True:
                assimilated_idx, monitored_idx = self.get_idx_conv(
                    obsid, subtype, analysis_use)
                lats = {'assimilated': self.lats[assimilated_idx],
                        'monitored': self.lats[monitored_idx]}
                lons = {'assimilated': self.lons[assimilated_idx],
                        'monitored': self.lons[monitored_idx]}
                return lats, lons
            else:
                idx = self.get_idx_conv(obsid, subtype, analysis_use)
                return self.lats[idx], self.lons[idx]

    def get_pressure(self, obsid=None, subtype=None, analysis_use=False):
        if analysis_use == True:
            assimilated_idx, monitored_idx = self.get_idx_conv(
                obsid, subtype, analysis_use)
            pressure = {'assimilated': self.press[assimilated_idx],
                        'monitored': self.press[monitored_idx]}

            return pressure
        else:
            idx = self.get_idx_conv(obsid, subtype, analysis_use)
            return self.press[idx]

    def metadata(self):
        dic = self.get_metadata()
        return dic


class Radiance(GSIdiag):

    def __init__(self, path):
        """
        Initialize a radiance GSI diagnostic object
        INPUT:
            path   : path to conventional GSI diagnostic object
        RESULT:
            self   : GSIdiag object containing the path to extract data
        """
        super().__init__(path)

        self.read_radiance_obs()

    def __str__(self):
        return "radiance object"

    def read_radiance_obs(self):
        """
        Reads the data from the radiance diagnostic file during initialization.
        """

        f = Dataset(self.path, mode='r')

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

        f.close()

    def get_data(self, dtype, channel=None, qcflag=None,
                separate_channels=False, separate_qc=False, errcheck=True):
        """
        Given parameters, get the data from a radiance 
        diagnostic file.
        INPUT:
            required:
                dtype : type of data to extract i.e. observation, O-F, O-A, H(x), Water_Fraction,
                        Land_Fraction, Ice_Fraction, Snow_Fraction

            optional:  
                channel           : observation channel number
                qcflag            : qc flag (default: None) i.e. 0, 1
                separate_channels : if True, calls get_data_special() and returns dictionary
                                    of separate data by specified channels
                separate_qc       : if True, calls get_data_special() and returns dictionary
                                    of separate data by specified QC flags
                errcheck          : when true, and qc==0, will toss out obs where inverse
                                    obs error is zero (i.e. not assimilated in GSI)
        OUTPUT:
            data : requested data

        """
        if separate_channels == True or separate_qc == True:
            data = self.get_data_special(
                dtype, channel, qcflag, separate_channels, separate_qc, errcheck=errcheck)
            return data

        else:
            idx = self.get_idx_sat(channel, qcflag, errcheck=errcheck)

            data = self.query_data_type(dtype, idx)

            data[data > 1e5] = np.nan

            return data

    def get_idx_sat(self, channel=None, qcflag=None, errcheck=True):
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

    def get_data_special(self, dtype, channel, qcflag,
                        separate_channels, separate_qc, errcheck=True):
        """
        Creates a dictionary that separates channels and qc flags
        depending on the conditions of seperate_channels and
        separate_qc
        """
        data_dict = {}

        if separate_channels == True and separate_qc == False:
            for c in channel:
                idx = self.get_idx_sat(c, qcflag, errcheck=errcheck)

                data = self.query_data_type(dtype, idx)

                data[data > 1e5] = np.nan

                data_dict['Channel_%s' % c] = data

            return data_dict

        if separate_channels == False and separate_qc == True:
            for qc in qcflag:
                idx = self.get_idx_sat(channel, qc, errcheck=errcheck)

                data = self.query_data_type(dtype, idx)

                data[data > 1e5] = np.nan

                data_dict['QC_Flag_%s' % qc] = data

            return data_dict

        if separate_channels == True and separate_qc == True:
            for c in channel:
                data_dict['Channel_%s' % c] = {}
                for qc in qcflag:
                    idx = self.get_idx_sat(c, qc, errcheck=errcheck)

                    data = self.query_data_type(dtype, idx)

                    data[data > 1e5] = np.nan

                    data_dict['Channel_%s' % c]['QC_Flag_%s' % qc] = data

            return data_dict

    def get_lat_lon(self, channel=None, qcflag=None, errcheck=True):
        """
        Gets lats and lons with desired indices
        """
        idx = self.get_idx_sat(channel, qcflag, errcheck=errcheck)
        return self.lats[idx], self.lons[idx]
