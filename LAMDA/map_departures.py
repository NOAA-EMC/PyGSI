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


def _create_map_departures(df, domain, projection,
                           metadata, outdir):
    """
    Create the map figure and plot data.
    """
    # this for now puts all O-F/O-A together, might want to split by
    # monitored/rejected/assim later
    # grab variables for plotting
    lats = df['latitude'].to_numpy()
    lons = df['longitude'].to_numpy()
    omf = df['omf_adjusted'].to_numpy()

    plot_objects = []

    # Create map object
    mymap = CreateMap(figsize=(12, 8),
                      domain=Domain(domain),
                      proj_obj=MapProjection(projection))
    # Add coastlines and states
    mymap.add_features(['coastlines', 'states'])

    # determine vmin/vmax for colorbar and must be symmetric
    if len(lats) > 0:
        omfscatter = MapScatter(latitude=lats,
                                longitude=lons,
                                data=omf)
        omfscatter.markersize = 2
        omfscatter.cmap = 'coolwarm'
        # determine min/max of values
        maxval = np.nanmax(omf)
        minval = np.nanmin(omf)
        vmax = max(abs(minval), maxval)
        omfscatter.vmin = vmax * -1
        omfscatter.vmax = vmax
        # add data to plot objects
        plot_objects.append(omfscatter)
        # labels and annotations
        labels = features.get_labels(metadata)
        # draw data
        mymap.draw_data(plot_objects)
        # add labels
        mymap.add_colorbar(label=metadata['Diag Type'],
                           label_fontsize=12, extend='neither')
        mymap.add_title(labels['title'], loc='left', fontsize=12)
        mymap.add_title(labels['date title'], loc='right', fontsize=12,
                        fontweight='semibold')
        mymap.add_xlabel("Longitude", fontsize=12)
        mymap.add_ylabel("Latitude", fontsize=12)
        # add stats to figure
        stats_dict = {
            'nobs': len(lats),
            'min': str(np.round(minval, 4)),
            'max': str(np.round(maxval, 4)),
        }
        mymap.add_stats_dict(stats_dict=stats_dict)
        # save figure
        fig = mymap.return_figure()
        str_domain = domain.replace(" ", "_")
        plt.savefig(outdir + f"{labels['save file']}_{str_domain}.png",
                    bbox_inches='tight', pad_inches=0.1)
        plt.close('all')

    else:
        # no data to plot
        fig = no_data_map(mymap, Domain(domain), metadata)


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
        metadata['ObsID Name'] = features.get_obs_type([config['observation id']])

    else:
        diag = Radiance(config['diag file'])
        df = diag.get_data(channel=[config['channel']],
                           analysis_use=config['analysis use'])
        metadata = diag.metadata

    metadata['Diag Type'] = config['data type']

    anl_use = config['analysis use']
    if anl_use:
        for anl_type in df.keys():
            metadata['Anl Use Type'] = anl_type
            _create_map_departures(df[anl_type], config['domain'], config['projection'],
                                   metadata, config['outdir'])
    else:
        metadata['Anl Use Type'] = None
        _create_map_departures(df, config['domain'], config['projection'],
                               metadata, config['outdir'])
