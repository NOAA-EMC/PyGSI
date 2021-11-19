import numpy as np
import yaml
import os
from pathlib import Path
import matplotlib.pyplot as plt
import LAMDA.plot_features as features
from LAMDA.no_data_plots import no_data_map
from emcpy.plots.map_plots import MapScatter
from emcpy.plots import CreateMap
from emcpy.plots.map_tools import Domain, MapProjection
from pyGSI.diags import Conventional, Radiance, Ozone
from emcpy.plots.plots import Histogram, LinePlot
from emcpy.plots.create_plots import CreatePlot
import matplotlib.mlab as mlab
import scipy.stats as stats

def create_bins(data,nbins=50):
   """
   Creates an array of the bins to be used.
   
   Args:
       data : (array) data being used in histogram
       nbins: (int;default=50) number of bins to use
   """
   _max = np.max(np.abs(data))
   bins = np.linspace(-_max,_max, nbins)

   return bins
   

def _create_hist_layer(df_ges,df_anl, domain,
                   metadata, outdir):
    """
    Create the ilayered histogram figure and plot data omf vs oma.
    """
    plot_objects = []

    lats = df_ges['latitude'].to_numpy()
    lons = df_anl['longitude'].to_numpy()
    omf = df_ges['omf_adjusted'].to_numpy()
    oma = df_anl['omf_adjusted'].to_numpy()

    omf_density = stats.gaussian_kde(omf)
    oma_density = stats.gaussian_kde(oma)

    omf_bins = create_bins(omf,50)
    oma_bins = create_bins(oma,50)


    # Create histogram objects
    hst_ges = Histogram(omf)
    hst_ges.color = 'tab:green'
    hst_ges.alpha = 0.7
    hst_ges.label = 'O-B'
    hst_ges.bins = 50
    hst_ges.density = True
    hst_ges.histtype = 'step'

    hst_anl = Histogram(oma)
    hst_anl.color = 'tab:purple'
    hst_anl.alpha = 0.7
    hst_anl.label = 'O-A'
    hst_anl.bins = 50
    hst_anl.density = True
    hst_anl.histtype = 'step'

    # Creat Line object using histogram bins and the density
    omf_line = LinePlot(omf_bins, omf_density(omf_bins))
    omf_line.color='tab:green'
    omf_line.label='O-B'
    oma_line = LinePlot(oma_bins, oma_density(oma_bins))
    oma_line.color='tab:purple'
    oma_line.label='O-A'

    # Create histogram plot and draw data
    myplt = CreatePlot()
    #plt_list = [hst_ges, hst_anl]
    plt_list = [omf_line, oma_line]
    myplt.draw_data(plt_list)

    # Add features
    labels = features.get_labels(metadata)
    # Titles
    labels = features.get_labels(metadata)
    myplt.add_title(labels['title'], loc='left', fontsize=12)
    myplt.add_title(labels['date title'], loc='right', fontsize=12,
                    fontweight='semibold')
    myplt.add_xlabel(xlabel='FG Departure (K)')
    myplt.add_ylabel(ylabel='Normalized Count')
    myplt.add_legend()
   
    # Return matplotlib figure
    fig = myplt.return_figure()
    str_domain = domain.replace(" ", "_")
    fig.savefig(outdir + f"{labels['save file']}_{str_domain}_layered_histogram.png",
                bbox_inches='tight', pad_inches=0.1)
    plt.close('all')
    return

def layer_histogram(config):
    """
    Create layered histogram plots with omf and oma.

    Args:
        config : (dict) configuration file that includes the
                 appropriate inputs based on file type (i.e.
                 conventional or radiance data)
    """

    # Get filename to determing what the file type is
    filename = os.path.splitext(Path(config['diag file']).stem)[0]
    filetype = filename.split('_')[1]

    anl_filename = config['diag file'].replace('ges', 'anl') #get the anl diag file
    print('anl_filename=',anl_filename)

    if filetype == 'conv':
        diag = Conventional(config['diag file'])

        df = diag.get_data(obsid=[config['observation id']],
                           subtype=[config['observation subtype']],
                           analysis_use=config['analysis use'])
        metadata = diag.metadata
        metadata['ObsID Name'] = features.get_obs_type(
            [config['observation id']])

    else:
        diag_ges = Radiance(config['diag file'])
        print('diag_ges=',diag_ges)
        diag_anl = Radiance(anl_filename)
        print('diag_anl=',diag_anl)

        df_ges = diag_ges.get_data(channel=[config['channel']], qcflag=[config['qc flag']],
                           analysis_use=config['analysis use']) 
        df_anl = diag_anl.get_data(channel=[config['channel']], qcflag=[config['qc flag']],
                           analysis_use=config['analysis use']) 
        print('df_anl=',df_anl)
        print('df_ges=',df_ges)
        metadata = diag_ges.metadata

        # Grab qc flags
        qc_unique = sorted(np.unique(np.abs(diag_ges.qc_flags)))

    metadata['Diag Type'] = 'Obs_Minus_Forecast_adjusted'

    # Handles analysis use data
    anl_use = metadata['Anl Use']

    if anl_use:
        for anl_type in df_ges.keys():
            metadata['Anl Use Type'] = anl_type

            if filetype == 'conv':
                # Need to grab qc_unique here for conv based on
                # analysis type (assimilated, rejected, monitored)
                qc_unique = sorted(np.unique(
                    np.abs(df[anl_type]['prep_qc_mark'])))

            _create_hist_layer(df_ges['assimilated'],df_anl['assimilated'],config['domain'],
                               metadata,config['outdir'])

    else:
        metadata['Anl Use Type'] = None

        if filetype == 'conv':
            # Need to grab qc_unique here for conv based on
            # analysis type (assimilated, rejected, monitored)
            qc_unique = sorted(np.unique(
                np.abs(df['prep_qc_mark'])))

            _create_hist_layer(df_ges['assimilated'],df_anl['assimilated'],config['domain'],
                               metadata,config['outdir'])
