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
            data = self.omf[idx]
            return data
        
        elif dtype == 'observation':
            data = self.observation[idx]
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
            Hx = get_H_x()
            data = Hx[idx]
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
        
    def get_H_x(self):
        """
        Calculates H(x): O - (O-F) = F = H(x)
        """
        
        Hx = self.o - self.omf
        return Hx        
        
    
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
        
#         self.get_metadata()
        self.read_conv_obs()
        
        
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
        self.stnelev = f.variables['Station_Elevation'][:]
        self.o = f.variables['Observation'][:]
        
        if self.path.split('/')[-1].split('.')[0].split('_')[2] == 'uv':
            self.u_o = f.variables['u_Observation'][:]
            self.v_o = f.variables['v_Observation'][:]
            
            self.u_omf = f.variables['u_Obs_Minus_Forecast_adjusted'][:]
            self.v_omf = f.variables['v_Obs_Minus_Forecast_adjusted'][:]
            
        else:
            self.o = f.variables['Observation'][:]
            self.omf = f.variables['Obs_Minus_Forecast_adjusted'][:]
        
        f.close()

            
    def getData(self, dtype, obsid=None, qcflag=None):
        """
        Given parameters, get the data from a conventional diagnostic file
        INPUT:
            dtype  : type of data to extract i.e. observation, O-F, O-A, H(x)
            obsid  : observation measurement ID number
            qcflag : qc flag (default: None) i.e. 0, 1
        OUTPUT:
            data   : requested data
        """
        
        idx = self.get_idx_conv(obsid, qcflag=None)
        
        data = self.query_dataType(dtype, idx)
        
        return data
    
    
    def get_idx_conv(self, obsid=None, qcflag=None):
        """
        Given parameters, get the indices of the observation
        locations from a conventional diagnostic file
        INPUT:
            obsid  : observation measurement ID number
            qcflag : qc flag (default: None) i.e. 0, 1
        OUTPUT:
            idx    : indices of the requested data in the file
        """
        
        if not obsid:
            obsid=None
        if not qcflag:
            qcflag=None
        
        idx = self.o_type
        
        if obsid != None:
            valid_idxs = np.isin(idx, obsid)
            idx = np.where(valid_idxs)
        if qcflag != None:
            valid_idxs = np.isin(idx, qcflag)
            idx = np.where(valid_idxs)
        
        return idx
    
    def get_lat_lon(self, obsid=None, qcflag=None):
        """
        Gets lats and lons with desired indices
        """
        idx = self.get_idx_conv(obsid, qcflag)
        return self.lats[idx], self.lons[idx]
    
    def metadata(self):
        dic = self.get_metadata()
        return dic
    
    
class satellite(gsidiag):
    
    def __init__(self, path):
        """
        Initialize a satellite GSI diagnostic object
        INPUT:
            path   : path to conventional GSI diagnostic object
        RESULT:
            self   : GSIdiag object containing the path to extract data
        """
        super().__init__(path)

        self.read_satellite_obs()
        
    def read_satellite_obs(self):
        """
        Reads the data from the satellite diagnostic file during initialization.
        """
        
        f = Dataset(self.path, mode='r')
        
        self.channel_idx = f.variables['Channel_Index'][:]
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
        
        f.close()
        
        
    def getData(self, dtype, channel=None, qcflag=None):
        """
        Given parameters, get the data from a satellite 
        diagnostic file.
        INPUT:
            dtype   : type of data to extract i.e. observation, O-F, O-A, H(x)
            channel : observation channel number
            qcflag  : qc flag (default: None) i.e. 0, 1
        OUTPUT:
            data    : requested data
            
        """
        
        idx = self.get_idx_sat(channel, qc_flag)
        
        data = self.query_dataType(dtype, idx)
        
        return data
    
    def get_idx_sat(self, channel=None, qcflag=None):
        """
        Given parameters, get the indices of the observation
        locations from a satellite diagnostic file.
        """
        
        if not channel:
            channel=None
        if not qcflag:
            qcflag=None
            
        idx = self.channel_idx
        
        if channel != None:
            valid_idxs = np.isin(idx, channel)
            idx = np.where(valid_idxs)
        if qcflag != None:
            valid_idxs = np.isin(self.qc_flag[idx], qcflag)
            idx = np.where(valid_idxs)
            
        return idx
    
    # How can we get this method to be used for both conv
    # and satellite
    def get_lat_lon(self, channel=None, qcflag=None):
        """
        Gets lats and lons with desired indices
        """
        idx = self.get_idx_sat(channel, qcflag)
        return self.lats[idx], self.lons[idx]
    