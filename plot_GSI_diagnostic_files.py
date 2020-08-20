import numpy as np
from netCDF4 import Dataset
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
from scipy import spatial
import cartopy.crs as ccrs
import cartopy.feature as cfeature
import argparse

parser = argparse.ArgumentParser(description='Plot O-F histogram and spatial plot with GSI diagnostic files.')          
parser.add_argument('-f', '--file', type=str,                                                                                                        
                    help='path to GSI diagnostic file.', required=True)
parser.add_argument('-o', '--obs', type=int, default=None,                                                                                                        
                    help='the specific observation code.')

args = parser.parse_args()

nc_file = args.file
obs_id = args.obs


################################################################## 

def get_metadata(data_path):
    """
    Input: path to data
    Output:
        var: variable being analyzed
        date: date 'YYYYMMDDHH'
        hour: hour HH (zulu)
        file_type: either ges or anl
    """
    split = data_path.split('/')[-1].split('.')
    var = split[0].split('_')[2]
    date = split[1].split('_')[0]
    hour = date[-2:]
    file_type = split[0].split('_')[-1]
    
    return var, date, hour, file_type

def calculate_stats(o_f):
    mean = np.mean(o_f)
    var_list = [(x-mean)**2 for x in o_f]
    variance = np.sum(var_list)/(len(var_list)-1)
    std = np.sqrt(variance)
    
    return mean, std


def plot_histogram(o_f, bins, meta_data):
    # get count and calculate mean and standard deviation of o_f
    n = len(o_f) 
    mean, std = calculate_stats(o_f)
    
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111)
    plt.hist(o_f, bins=bins)
    if meta_data['Variable'] == 'q':
        t = ('n: %s\nstd: %s\nmean: %s' % (n,np.round(std,6),np.round(mean,6)))
        ax.text(0.75,.8, t, fontsize=14, transform=ax.transAxes)
    else:
        t = ('n: %s\nstd: %s\nmean: %s' % (n,np.round(std,3),np.round(mean,3)))
        ax.text(0.75,.8, t, fontsize=14, transform=ax.transAxes)
    plt.xlabel('O - F')
    plt.ylabel('Count')
    plt.title('%s%s_%s:%s,O-F all data on %s' % (meta_data['Variable'],meta_data['Obs_ID'],meta_data['Hour'],meta_data['File_type'],meta_data['Date']), fontsize=14)
    plt.savefig('%s_%s_%s_O_minus_F_histogram.png' % (meta_data['Date'],meta_data['Variable'],meta_data['Obs_ID']), bbox_inches='tight', pad_inches=0.1)
    
    return 0

def plot_spatial(o_f, bounds, meta_data, lons, lats):
    
    n = len(o_f) 
    mean, std = calculate_stats(o_f)
    
    plt.figure(figsize=(15,12))
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
    ax.coastlines()
    norm = mcolors.BoundaryNorm(boundaries=bounds, ncolors=256)
    
    cs = plt.scatter(lons, lats, c=o_f, s=30,
                norm=norm, cmap='bwr', #edgecolors='gray', linewidth=0.25,
                transform=ccrs.PlateCarree())
    
    t = ('n: %s\nstd: %s\nmean: %s' % (n,std,mean))
    ax.text(-175, -70, t, fontsize = 16, transform=ccrs.PlateCarree())

    cb = plt.colorbar(cs, shrink=0.5, pad=.04, extend='both')
    cb.set_label('O - F', fontsize=12)
    plt.title('%s%s_%s:%s,O-F all data on %s' % (meta_data['Variable'],meta_data['Obs_ID'],meta_data['Hour'],meta_data['File_type'],meta_data['Date']), fontsize=14)
    plt.savefig('%s_%s_%s_O_minus_F_spatial.png' % (meta_data['Date'],meta_data['Variable'],meta_data['Obs_ID']), bbox_inches='tight', pad_inches=0.1)
    
    return 0
    
    
##################################################################    

var, date, hour, file_type = get_metadata(nc_file)

if obs_id != None:
                
    meta_data = {"Variable": var,
                 "Date": date,
                 "Hour": hour,
                 "File_type": file_type,
                 "Obs_ID": obs_id
                }
else:
    meta_data = {"Variable": var,
                 "Date": date,
                 "Hour": hour,
                 "File_type": file_type,
                 "Obs_ID": "Total"
                }


## Read data
f = Dataset(nc_file, mode='r')

lons = f.variables['Longitude'][:]
lats = f.variables['Latitude'][:]
o_f  = f.variables['Obs_Minus_Forecast_adjusted'][:]
o_type = f.variables['Observation_Type'][:]
f.close()
    

## Find data with indicated observation type
if obs_id != None:    
    idx = np.where(o_type == obs_id)
    o_f = o_f[idx]
    
    lons = lons[idx]
    lats = lats[idx]

if o_f.size == 0:
    print("No observations for %s from Observation ID: %s" % (var, obs_id))

# Temperature
if var == 't':
    bins = np.arange(-10,10.1,0.1)
    plot_histogram(o_f, bins, meta_data)
    
    bounds = np.arange(-10,12.5,2.5)
    plot_spatial(o_f, bounds, meta_data, lons, lats)

# Specific Humidity
if var == 'q':
    bins = np.arange(-0.01,0.0101,0.0001)
    plot_histogram(o_f, bins, meta_data)
    
    bounds = np.arange(-0.01,0.011,0.0025)
    plot_spatial(o_f, bounds, meta_data, lons, lats)

# Pressure
if var == 'ps':
    bins = np.arange(-500,510,10)
    plot_histogram(o_f, bins, meta_data)
    
    bounds = np.arange(-500,510,100)
    plot_spatial(o_f, bounds, meta_data, lons, lats)

