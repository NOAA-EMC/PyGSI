import sys

sys.path.append('/home/Kevin.Dougherty/GSI_plots')

from runDiagnostics import conventional



path = '/scratch1/NCEPDEV/stmp2/Cory.R.Martin/ICs_tmp/v16para/gdas.20200531/00/'
nc_file = path+'diag_conv_t_ges.2020053100.nc4'

obs_id  = []
qc_flag = []

# DATA_TYPE can be: 'O-F', 'O-A', 'observation', 'H(x)'  (so far)
DATA_TYPE = 'O-F'

# Create Conventional Diagnostic data object
diag = conventional(nc_file)

# Get data from specified DATA_TYPE, include specified observation ids and qc flags (Don't know if there are any for conventional data)
data = diag.getData(DATA_TYPE, obs_id, qc_flag)

# Get lat and lon of data for spatial plots. Include specified observation ids and qc flags
lat, lon = diag.get_lat_lon()