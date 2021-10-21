import numpy as np
import yaml
import os
from pathlib import Path
import matplotlib.pyplot as plt
import plot_features as features
from no_data_plots import no_data_map
from emcpy.plots.map_plots import MapScatter
from emcpy.plots import CreateMap
from pyGSI.diags import Conventional, Radiance, Ozone


def _create_map_qc(df, qc_unique, domain, projection,
                   metadata, plotdir):
    """
    Create the map figure and plot data.
    """
    plot_objects = []
    # Loop through unique QC Flags to create plot objects
    for i, flag in enumerate(qc_unique):
        if metadata['Diag File Type'] == 'conventional':
            indx = data_df.index[df['prep_qc_mark'] == flag]
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
    mymap.add_features(['coastlines'])

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

    plt.savefig(plotdir + f"{labels['save file']}_map.png",
                bbox_inches='tight', pad_inches=0.1)
    plt.close('all')

    return


def map_qc_flags(config_file):
    """
    Create map plotting location of qcflags.

    Args:
        config_file : (dict) configuration file that includes the
                      appropriate inputs based on file type (i.e.
                      conventional or radiance data)
        inputfile : (str) path to diagnostic file
        channel : (list of ints; default=None) channel number
                  to plot
        qcflag : (list of ints; default=None) qc flags to
                 plot
        analysis_use : (bool; default=False) if True, will return
                       three sets of data:
                       assimilated (QC_Flag=0, inv_observation_error!=0),
                       rejected (QC_Flag!=0),
                       monitored (use_flag!=1)
        domain : (str; default='conus') domain in which to plot data
        projection : (str; default='plcarr') projection of map to plot
                     data
        plotdir : (str; default='./') path to where figures should be saved
    """

    # Get filename to determing what the file type is
    filename = os.path.splitext(Path(inputfile).stem)[0]
    filetype = filename.split('_')[1]

    if filetype == 'conv':
        diag = Conventional(inputfile)

        df = diag.get_data(obsid=obsid, subtype=subtype, station_id=station_id,
                           analysis_use=analysis_use)

        qc_unique = sorted(np.unique(np.abs(df['prep_qc_mark'])))

    else:
        diag = Radiance(inputfile)

        df = diag.get_data(channel=channel, qcflag=qcflag,
                           analysis_use=analysis_use)

        # Grab qc flags
        qc_unique = sorted(np.unique(np.abs(diag.qc_flags)))

    metadata = diag.metadata
    metadata['Diag Type'] = 'QC Flags'

    # Handles analysis use data
    anl_use = metadata['Anl Use']

    if anl_use:
        for anl_type in data.keys():
            metadata['Anl Use Type'] = anl_type

            _create_map_qc(df[anl_type], qc_unique,
                           domain, projection, metadata,
                           plotdir)

    else:
        metadata['Anl Use Type'] = None
        _create_map_qc(df, qc_unique, domain,
                       projection, metadata, plotdir)
