import netCDF4 as nc
import numpy as np
import datetime as dt


class Increment:

    def __init__(self, path):
        """
        Initialize a GSI increment file object
        INPUT:
            path : path to netCDF GSI increment file
        """

        self.path = path
        self.get_metadata()

    def get_metadata(self):
        """
        Grabs metadata from the netCDF file
        """

        inc = nc.Dataset(self.path, 'r')
        # newer versions of the GSI increment file
        # include metadata in the global attributes
        # if these don't exist, just return dummy vars
        try:
            anltimeint = inc.getncattr('analysis_time')
            anltime = dt.datetime.strptime(str(anltimeint), '%Y%m%d%H')
            IAUhr = int(inc.getncattr('IAU_hour_from_guess'))
            tdiff = 6-IAUhr
            if tdiff > 0:
                validtime = anltime + dt.timedelta(hours=tdiff)
            elif tdiff < 0:
                validtime = anltime - dt.timedelta(hours=abs(tdiff))
            else:
                validtime = anltime
        except:
            validtime = dt.datetime(2000, 1, 1)
            anltime = dt.datetime(2000, 1, 1)
            IAUhr = -9999
        inc.close()
        self.validtime = validtime
        self.analysistime = anltime
        self.IAUhr = IAUhr

    def get_latlonlevs(self):
        """
        Get lat/lon/levs values from increment file
        """
        inc = nc.Dataset(self.path, 'r')
        lat = inc.variables['lat'][:]
        lon = inc.variables['lon'][:]
        lon[lon > 180.] = lon[lon > 180.]-360.
        lon = np.roll(lon, int(len(lon)/2), axis=0)
        lons, lats = np.meshgrid(lon, lat)
        levs = inc.variables['lev'][:]
        inc.close()
        return lats, lons, levs

    def get_increment(self, varname, lev=None):
        """
        Get increment field in either 2 or 3D
        call it as so:
        for 2D
        T_inc = incobject.get_increment('T_inc', lev=64)
        for 3D
        T_inc3D = incobject.get_increment('T_inc')
        """
        inc = nc.Dataset(self.path, 'r')
        if lev:
            vardata = inc.variables[varname][lev, ...]
        else:
            vardata = inc.variables[varname][:]
        lon = inc.variables['lon'][:]
        vardata = np.roll(vardata, int(len(lon)/2), axis=1)
        return vardata
