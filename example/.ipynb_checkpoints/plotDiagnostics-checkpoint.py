import numpy as np
from netCDF4 import Dataset
import yaml
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import cartopy.crs as ccrs
import cartopy.feature as cfeature
from textwrap import wrap
import sys


sys.path.append('/home/Kevin.Dougherty/GSI_plots')

from runDiagnostics import conventional

file = open('test_YAML.yaml')
parsed_yaml_file = yaml.load(file, Loader=yaml.FullLoader)


def calculate_stats(diff):
    mean = np.mean(diff)
    var_list = [(x-mean)**2 for x in diff]
    variance = np.sum(var_list)/(len(var_list)-1)
    std = np.sqrt(variance)
    mx = max(diff)
    mn = min(diff)
    
    return mean, std, mx, mn   

def get_obs_type(obs_id):
    
    obs_indicators = {
        120: "Rawinsonde",
        126: "RASS",
        130: "Aircraft: AIREP and PIREP",
        131: "Aircraft: AMDAR",
        132: "Flight-Level Reconnaissance and Profile Dropsonde",
        133: "Aircraft: MDCRS ACARS",
        134: "Aircraft: TAMDAR",
        135: "Aircraft: Canadian AMDAR",
        153: "GPS-Integrated Precipitable Water",
        180: "Surface Marine w/ Station Pressure (Ship, Buoy, C-MAN, Tide Guage)",
        181: "Surface Land w/ Station Pressure (Synoptic, METAR)",
        182: "Splash-Level Dropsonde Over Ocean",
        183: "Surface Marine or Land - Missing Station Pressure",
        187: "Surface Land - Missing Station Pressure",
        210: "Synthetic Tropical Cyclone",
        220: "Rawinsonde",
        221: "PIBAL",
        224: "NEXRAD Vertical Azimuth Display",
        228: "Wind Profiler: JMA",
        229: "Wind Profiler: PIBAL",
        230: "Aircraft: AIREP and PIREP",
        231: "Aircraft: AMDAR",
        232: "Flight-Level Reconnaissance and Profile Dropsonde",
        233: "Aircraft: MDCRS ACARS",
        234: "Aircraft: TAMDAR",
        235: "Aircraft: Canadian AMDAR",
        242: "JMA IR (Longwave) and Visible Cloud Drift Below 850mb (GMS, MTSAT, HIMAWARI)",
        243: "EUMETSAT IR (Longwave) and Visible Cloud Drift Below 850mb (METEOSAT)",
        245: "NESDIS IR (Longwave) Cloud Drift (All Levels)",
        246: "NESDIS Imager Water Vapor (All Levels) - Cloud Top (GOES)",
        250: "JMA Imager Water Vapor (All Levels) - Cloud Top & Deep Layer (GMS, MTSAT, HIMAWARI)",
        251: "NESDIS Visible Cloud Drift (All Levels) (GOES)",
        252: "JMA IR (Longwave) and Visible Cloud Drift Above 850mb (GMS, MTSAT, HIMAWARI)",
        253: "EUMETSAT IR (Longwave) and Visible Cloud Drift Above 850mb (METEOSAT)",
        254: "EUMETSAT Imager Water Vapor (All Levels) - Cloud Top & Deep Layer (METEOSAT)",
        257: "MODI/POES IR (Longwave) Cloud Drift (All Levels) (AQUA, TERRA)",
        258: "MODIS/POES Imager Water Vapor (All Levels) - Cloud Top (AQUA, TERRA)",
        259: "MODIS/POES Imager Water Vapor (All Levels) - Deep Layer (AQUA, TERRA)",
        280: "Surface Marine w/ Station Pressure (Ship, Buoy, C-MAN, Tide Guage)",
        281: "Surface Land w/ Station Pressure (Synoptic, METAR)",
        282: "ATLAS Buoy",
        284: "Surface Marine or Land - Missing Station Pressure",
        287: "Surface Land (METAR) - Missing Station Pressure",
        289: "SUPEROBED (1.0 Lat/Lon) Scatterometer Winds over Ocean",
        290: "Non-SUPEROBED Scatterometer Winds over Ocean"
    }
    
    descripts = list()
    for ids in obs_id:
        if (ids in obs_indicators.keys()) == True:
            descripts.append(obs_indicators[ids])
            
    if descripts:
        return descripts
    elif not obs_id:
        return ['All Observations']
    else:
        return [str(x) for x in obs_id]


def plot_histogram(diag):
    # get count and calculate mean and standard deviation of diff
    n = len(diag.diff) 
    mean, std, mx, mn = calculate_stats(diag.diff)
    
    obs_descript = get_obs_type(diag.obs_id)
    diag.obs_descript = obs_descript
    
    fig = plt.figure(figsize=(8,6))
    ax = fig.add_subplot(111)
    plt.hist(diag.diff, bins=diag.variables[diag.variable]['binsize'])
    
    if diag.dtype == 'conv' and diag.variable == 'q':
            t = ('n: %s\nstd: %s\nmean: %s' % (n,np.round(std,6),np.round(mean,6)))
            ax.text(0.75,.8, t, fontsize=14, transform=ax.transAxes)
    else:
        t = ('n: %s\nstd: %s\nmean: %s\nmax: %s\nmin: %s' % (n,np.round(std,3),np.round(mean,3), np.round(mx,3), np.round(mn,3)))
        ax.text(0.75,.7, t, fontsize=14, transform=ax.transAxes)
    
    labels = diag.get_labels()
        
    plt.xlabel(labels[diag.ftype]['xlabel'])
    plt.ylabel('Count')
    
    title_split = labels[diag.ftype]['left_title'].split('\n')
    plt.title("%s\n%s" % (title_split[0], '\n'.join(wrap(title_split[-1], 40))), loc='left', fontsize=14)
    plt.title(labels[diag.ftype]['right_title'], loc='right', fontweight='semibold', fontsize=14)
    plt.savefig(labels[diag.ftype]['save_title'], bbox_inches='tight', pad_inches=0.1)
    
    return 0

def plot_spatial(diag):
    
    n = len(diag.diff) 
    mean, std, mx, mn = calculate_stats(diag.diff)
    
    plt.figure(figsize=(15,12))
    ax = plt.axes(projection=ccrs.PlateCarree(central_longitude=0))
    ax.add_feature(cfeature.GSHHSFeature(scale='auto'))
    ax.set_extent([-180, 180, -90, 90])
    norm = mcolors.BoundaryNorm(boundaries=diag.variables[diag.variable]['boundsize'], ncolors=256)
    
    cs = plt.scatter(diag.lons, diag.lats, c=diag.diff, s=30,
                norm=norm, cmap='bwr', #edgecolors='gray', linewidth=0.25,
                transform=ccrs.PlateCarree())
    
    if diag.dtype == 'conv' and diag.variable == 'q':
        t = ('n: %s\nstd: %s\nmean: %s' % (n,np.round(std,6),np.round(mean,6)))
        ax.text(-175, -70, t, fontsize=16, transform=ccrs.PlateCarree())
    else:
        t = ('n: %s\nstd: %s\nmean: %s\nmax: %s\nmin: %s' % (n,np.round(std,3),np.round(mean,3), np.round(mx,3), np.round(mn,3)))
        ax.text(-175, -70, t, fontsize=16, transform=ccrs.PlateCarree())


    labels = diag.get_labels()
    
    cb = plt.colorbar(cs, shrink=0.5, pad=.04, extend='both')
    cb.set_label(labels[diag.ftype]['xlabel'], fontsize=12)
    
    title_split = labels[diag.ftype]['left_title'].split('\n')
    plt.title("%s\n%s" % (title_split[0], '\n'.join(wrap(title_split[-1], 40))), loc='left', fontsize=14)
    plt.title(labels[diag.ftype]['right_title'], loc='right', fontweight='semibold', fontsize=14)
    plt.savefig(labels[diag.ftype]['save_title'], bbox_inches='tight', pad_inches=0.1)
    
    return 0
        
#####################################################
def main(parsed_yaml_file):
    if parsed_yaml_file['conventional input']['path']:
        
        nc_file = parsed_yaml_file['conventional input']['path']
        obs_id = parsed_yaml_file['conventional input']['observation id']
        
        diag = conventional(nc_file, obs_id)
        
        if np.isin('histogram',parsed_yaml_file['conventional input']['plot type']) == True:
            plot_histogram(diag)
            
        if np.isin('spatial',parsed_yaml_file['conventional input']['plot type']) == True:
            plot_spatial(diag)

main(parsed_yaml_file)