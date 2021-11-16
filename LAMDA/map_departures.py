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


def _create_map_departures(df, qc_unique, domain, projection,
                           metadata, outdir):
    """
    Create the map figure and plot data.
    """
    plot_objects = []


def map_departures(config):
    """
    Create map of departures (O-F and O-A)

    Args:
        config : (dict) configuration file that includes the
                 appropriate inputs based on file type (i.e.
                 conventional or radiance data)
    """

    # Get filename to determing what the file type is
    filename = os.path.splitext(Path(config['diag file']).stem)[0]
    filetype = filename.split('_')[1]

    if filetype == 'conv':
        diag = Conventional(config['diag file'])

        df = diag.get_data(obsid=[config['observation id']],
                           subtype=[config['observation subtype']],
                           analysis_use=config['analysis use'])
        metadata = diag.metadata
        metadata['ObsID Name'] = features.get_obs_type(
            [config['observation id']])

    else:
        diag = Radiance(config['diag file'])

        df = diag.get_data(channel=[config['channel']],
                           qcflag=[config['qc flag']],
                           analysis_use=config['analysis use'])
        metadata = diag.metadata

    metadata['Diag Type'] = 'Departures'

    # Handles analysis use data
    anl_use = metadata['Anl Use']
