import netCDF4 as nc
import numpy as np
import datetime as dt

class increment:

    def __init__(self, path):
        """
        Initialize a GSI increment file object
        INPUT:
            path : path to netCDF GSI increment file
        """

        self.path = path

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
            anltime = dt.datetime.strptime(str(anltimeint),'%Y%m%d%H')
            IAUhr = int(inc.getncattr('IAU_hour_from_guess'))
            tdiff = 6-IAUhr
            if tdiff > 0:
                validtime = anltime + dt.timedelta(hours=tdiff)
            elif tdiff < 0:
                validtime = anltime - dt.timedelta(hours=abs(tdiff))
            else:
                validtime = anltime
        except:
            validtime = dt.datetime(2000,1,1)
            anltime = dt.datetime(2000,1,1)
            IAUhr = -9999
        metadata = {
                    'valid time': validtime,
                    'analysis time': anltime,
                    'IAU forecast hour': IAUhr,
                   }
        return metadata
