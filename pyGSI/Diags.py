import numpy as np
from netCDF4 import Dataset
from datetime import datetime

class gsidiag:
    
    def __init__(self, path):
        """
        Initialize a GSI diagnostic object
        INPUT:
            path : path to GSI diagnostic object
        """

        self.path = path
        
    def __len__(self):
        return len(self.lats)
        
        
    def get_metadata(self):
        """
        Grabs metadata from the diagnostic filename **NEED TO ADDRESS STILL**
        """
        
        dtype = self.path.split('/')[-1].split('.')[0].split('_')[1]
        str_date = self.path.split('/')[-1].split('.')[1]
        date = datetime.strptime(str_date, '%Y%m%d%H')
        ftype = self.path.split('/')[-1].split('.')[0].split('_')[-1]
        
        if dtype == 'conv':
            variable = self.path.split('/')[-1].split('.')[0].split('_')[2]
        
            metadata = {'Diag_type' : dtype,
                        'Variable'  : variable,
                        'Date'      : date,
                        'File_type' : ftype
                       }
        else:
            satellite = self.path.split('/')[-1].split('.')[0].split('_')[2]
            
            metadata = {'Diag_type' : dtype,
                        'Satellite' : satellite,
                        'Date'      : date,
                        'File_type' : ftype
                       }
        
        return metadata
        
    def query_dataType(self, dtype, idx):
        """
        Query the data type being requested and returns
        the appropriate indexed data
        """
        if dtype == 'O-F':
            if self.path.split('/')[-1].split('.')[0].split('_')[2] == 'uv':
                u = self.u_omf[idx]
                v = self.v_omf[idx]
                
                return u,v
            
            else:
                data = self.omf[idx]
                return data
        
        elif dtype == 'Observation':
            if self.path.split('/')[-1].split('.')[0].split('_')[2] == 'uv':
                u = self.u_o[idx]
                v = self.v_o[idx]
                
                return u,v
            
            else:
                data = self.o[idx]
                return data
        
        elif dtype == 'O-A':
            # Not sure if I need to do this because 'ges' and 'anl' files are both f.variable['Obs_Minus_Forecast_adjusted'][:]
            # but I will leave it for now
            check = check_o_minus_a()
            if check:
                data = self.omf[idx]
                return data
            else:
                print('File type does not support O-A')
        
        elif dtype == 'H(x)':
            Hx = self.o - self.omf
            data = Hx[idx]
            return data
        
        elif dtype == 'Water_Fraction':
            data = self.water_frac[idx]
            return data
        
        elif dtype == 'Land_Fraction':
            data = self.land_frac[idx]
            return data
        
        elif dtype == 'Ice_Fraction':
            data = self.ice_frac[idx]
            return data
        
        elif dtype == 'Snow_Fraction':
            data = self.snow_frac[idx]
            return data
        
        else:
            raise Exception(f'Unrecognizable data type: {dtype}')
            return None
            
        
    def check_o_minus_a(self):
        """
        Checks if the diasnostic file is an analyses file.
        """
        
        ftype = self.path.split('/')[-1].split('.')[0].split('_')[-1]
        
        if ftype == 'anl':
            return True
        else:
            return False        
        
    
class conventional(gsidiag):
    
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
    
    def __str__(self):
        return "Conventional object"
        
        
    def read_conv_obs(self):
        """
        Reads the data from the conventional diagnostic file during initialization.
        """
        
        f = Dataset(self.path, mode='r')
        
        self.o_type = f.variables['Observation_Type'][:]
        self.lons = f.variables['Longitude'][:]
        self.lats = f.variables['Latitude'][:]
        self.press = f.variables['Pressure'][:]
        self.time = f.variables['Time'][:]
        self.anl_use = f.variables['Analysis_Use_Flag'][:]
#         self.prepqc = f.variables['Prep_QC_Mark'][:]
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

            
    def getData(self, dtype, obsid=None, analysis_use=False):
        """
        Given parameters, get the data from a conventional diagnostic file
        INPUT:
            required:
                dtype  : type of data to extract i.e. observation, O-F, O-A, H(x)
                
            optional:    
                obsid        : observation measurement ID number; default=None
                analysis_use : if True, will return two sets of data: assimlated
                               (analysis_use_flag=1), and monitored (analysis_use
                               _flag=-1); default = False
                
        OUTPUT:
            data   : requested data
        """
        if analysis_use == True:
            assimilated_idx, monitored_idx = self.get_idx_conv(obsid, analysis_use)
            if self.path.split('/')[-1].split('.')[0].split('_')[2] == 'uv':
                u_assimilated, v_assimilated = self.query_dataType(dtype, assimilated_idx)
                u_monitored, v_monitored = self.query_dataType(dtype, monitored_idx)
                
                u = {'assimilated': u_assimilated,
                     'monitored': u_monitored}
                v = {'assimilated': v_assimilated,
                     'monitored': v_monitored}
                
                return u, v
            else:
                assimilated_data = self.query_dataType(dtype, assimilated_idx)
                monitored_data = self.query_dataType(dtype, monitored_idx)
                
                data = {'assimilated': assimilated_data,
                        'monitored'  : monitored_data
                       }
                
                return data
            
        else:
            idx = self.get_idx_conv(obsid, analysis_use)

            if self.path.split('/')[-1].split('.')[0].split('_')[2] == 'uv':
                u, v = self.query_dataType(dtype, idx)

                return u, v

            else:
                data = self.query_dataType(dtype, idx)

                return data
    
    
    def get_idx_conv(self, obsid=None, analysis_use=False):
        """
        Given parameters, get the indices of the observation
        locations from a conventional diagnostic file
        INPUT:
            obsid  : observation measurement ID number
            qcflag : qc flag (default: None) i.e. 0, 1
        OUTPUT:
            idx    : indices of the requested data in the file
        """
        
        if not obsid or obsid == [None]:
            obsid=None
        
        idx = self.o_type
        
        if analysis_use == False:
            if obsid != None:
                valid_idxs = np.isin(idx, obsid)
                idx = np.where(valid_idxs)

            else:
                idx = np.where(idx)
        
            return idx
        
        else:
            if obsid != None:
                obs_idx = np.isin(idx, obsid)
                
                assimilated = np.isin(self.anl_use, 1)
                monitored = np.isin(self.anl_use, -1)
                
                valid_assimilated = np.logical_and(obs_idx, assimilated)
                valid_monitored = np.logical_and(obs_idx, monitored)
                
                assimilated_idx = np.where(valid_assimilated)
                monitored_idx = np.where(valid_monitored)
            
            else:
                assimilated = np.isin(self.anl_use, 1)
                monitored = np.isin(self.anl_use, -1)
                
                assimilated_idx = np.where(assimilated)
                monitored_idx = np.where(monitored)
                
            return assimilated_idx, monitored_idx
            
    
    def get_lat_lon(self, obsid=None, analysis_use=False):
        """
        Gets lats and lons with desired indices
        """
        if analysis_use == True:
            assimilated_idx, monitored_idx = self.get_idx_conv(obsid, analysis_use)
            lats = {'assimilated': self.lats[assimilated_idx],
                    'monitored': self.lats[monitored_idx]}
            lons = {'assimilated': self.lons[assimilated_idx],
                    'monitored': self.lons[monitored_idx]}
            return lats, lons
        else:
            idx = self.get_idx_conv(obsid, analysis_use)
            return self.lats[idx], self.lons[idx]
    
    def metadata(self):
        dic = self.get_metadata()
        return dic
    
    
class radiance(gsidiag):
    
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
        self.inv_ob_err = f.variables['Inverse_Observation_Error'][:]
        
        f.close()
        
        
    def getData(self, dtype, channel=None, qcflag=None,
                separate_channels=False, separate_qc=False):
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
                separate_channels : if True, calls getData_special() and returns dictionary
                                    of separate data by specified channels
                separate_qc       : if True, calls getData_special() and returns dictionary
                                    of separate data by specified QC flags
        OUTPUT:
            data : requested data
            
        """
        if separate_channels == True or separate_qc == True:
            data = self.getData_special(dtype, channel, qcflag, separate_channels, separate_qc)
            return data
        
        else:
            idx = self.get_idx_sat(channel, qcflag)

            data = self.query_dataType(dtype, idx)

            data[data > 1e5] = np.nan

            return data
    
    def get_idx_sat(self, channel=None, qcflag=None):
        """
        Given parameters, get the indices of the observation
        locations from a radiance diagnostic file.
        """
        
        if not channel or channel == [None]:
            channel=None
        if not qcflag or qcflag == [None]:
            qcflag=None
            
        
        if channel != None and qcflag != None:
            chidx = np.where(self.sensor_chan == channel)
            if len(chidx) > 0 and len(chidx[0]) > 0:
                channel = self.chaninfo_idx[chidx[0][0]]
            else:
                print('Channel specified not in sensor_chan, using relative index')
            idx = self.channel_idx
            chan_idxs = np.isin(self.channel_idx, channel)
    
            qc_idxs = np.isin(self.qc_flag, qcflag)

            valid_idxs = np.logical_and(chan_idxs, qc_idxs)
            idx = np.where(valid_idxs)
            
        elif channel != None and qcflag == None:
            chidx = np.where(self.sensor_chan == channel)
            if len(chidx) > 0 and len(chidx[0]) > 0:
                channel = self.chaninfo_idx[chidx[0][0]]
            else:
                print('Channel specified not in sensor_chan, using relative index')
            idx = self.channel_idx
            valid_idxs = np.isin(idx, channel)
            idx = np.where(valid_idxs)

        elif channel == None and qcflag != None:
            idx = self.qc_flag
            valid_idxs = np.isin(idx, qcflag)
            idx = np.where(valid_idxs)

        elif channel == None and qcflag == None:
            idx = np.where(idx)

        return idx
    
    def getData_special(self, dtype, channel, qcflag,
                        separate_channels, separate_qc):
        """
        Creates a dictionary that separates channels and qc flags
        depending on the conditions of seperate_channels and
        separate_qc
        """
        data_dict = {}
        
        if separate_channels == True and separate_qc == False:
            for c in channel:
                idx = self.get_idx_sat(c, qcflag)

                data = self.query_dataType(dtype, idx)

                data[data > 1e5] = np.nan

                data_dict['Channel_%s' % c] = data

            return data_dict
        
        if separate_channels == False and separate_qc == True:
            for qc in qcflag:
                idx = self.get_idx_sat(channel, qc)

                data = self.query_dataType(dtype, idx)

                data[data > 1e5] = np.nan

                data_dict['QC_Flag_%s' % qc] = data

            return data_dict
        
        if separate_channels == True and separate_qc == True:
            for c in channel:
                data_dict['Channel_%s' % c] = {}
                for qc in qcflag:
                    idx = self.get_idx_sat(c, qc)

                    data = self.query_dataType(dtype, idx)

                    data[data > 1e5] = np.nan

                    data_dict['Channel_%s' % c]['QC_Flag_%s' % qc] = data
                    
            return data_dict
                    
                    
    
    # How can we get this method to be used for both conv
    # and radiance
    def get_lat_lon(self, channel=None, qcflag=None):
        """
        Gets lats and lons with desired indices
        """
        idx = self.get_idx_sat(channel, qcflag)
        return self.lats[idx], self.lons[idx]
    
