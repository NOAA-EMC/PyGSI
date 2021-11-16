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


def _create_map_qc(df, qc_unique, domain, projection,
                   metadata, outdir):
    """
    Create the map figure and plot data.
    """
    plot_objects = []
    # Loop through unique QC Flags to create plot objects
    for i, flag in enumerate(qc_unique):
        if metadata['Diag File Type'] == 'conventional':
            indx = df.index[df['prep_qc_mark'] == flag]
        else:
            indx = (df.index.get_level_values('QC_Flag') == flag)

        lats = df['latitude'][indx].to_numpy()
        lons = df['longitude'][indx].to_numpy()

        # If data is not empty, creates scatter object
        if len(lats) > 0:
            plotobj = MapScatter(latitude=lats,
                                 longitude=lons)
            plotobj.color = features.qc_flag_colors(flag)
            plotobj.label = flag
            plot_objects.append(plotobj)

    # Create map object
    mymap = CreateMap(figsize=(12, 8),
                      domain=Domain(domain),
                      proj_obj=MapProjection(projection))
    # Add coastlines
    mymap.add_features(['coastlines', 'states'])

    # If there is no data, create figure that displays 'No Data'
    if len(plot_objects) == 0:
        fig = no_data_map(mymap, Domain(domain), metadata)

    else:
        # Draw data
        mymap.draw_data(plot_objects)

        # Add legend
        ncol = 2 if len(plot_objects) > 4 else 1
        legend = mymap.add_legend(loc='lower left', ncol=ncol,
                                  title='QC Flags')

        # Titles
        labels = features.get_labels(metadata)
        mymap.add_title(labels['title'], loc='left', fontsize=12)
        mymap.add_title(labels['date title'], loc='right', fontsize=12,
                        fontweight='semibold')

        # X and Y labels
        mymap.add_xlabel("Longitude", fontsize=12)
        mymap.add_ylabel("Latitude", fontsize=12)

        # Return figure
        fig = mymap.return_figure()

        str_domain = domain.replace(" ", "_")
        plt.savefig(outdir + f"{labels['save file']}_{str_domain}.png",
                    bbox_inches='tight', pad_inches=0.1)
        plt.close('all')

    return


def map_qc_flags(config):
    """
    Create map plotting location of qcflags.

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

        # Grab qc flags
        qc_unique = sorted(np.unique(np.abs(diag.qc_flags)))

    metadata['Diag Type'] = 'QC_Flags'

    # Handles analysis use data
    anl_use = metadata['Anl Use']

    if anl_use:
        for anl_type in df.keys():
            metadata['Anl Use Type'] = anl_type

            if filetype == 'conv':
                # Need to grab qc_unique here for conv based on
                # analysis type (assimilated, rejected, monitored)
                qc_unique = sorted(np.unique(
                    np.abs(df[anl_type]['prep_qc_mark'])))

            _create_map_qc(df[anl_type], qc_unique,
                           config['domain'], config['projection'],
                           metadata, config['outdir'])

    else:
        metadata['Anl Use Type'] = None

        if filetype == 'conv':
            # Need to grab qc_unique here for conv based on
            # analysis type (assimilated, rejected, monitored)
            qc_unique = sorted(np.unique(
                np.abs(df['prep_qc_mark'])))

        _create_map_qc(df, qc_unique, config['domain'],
                       config['projection'], metadata,
                       config['outdir'])
