import numpy as np
import matplotlib.pyplot as plt
import LAMDA.plot_features as features


def no_data_map(plotmap, domain, metadata):
    """
    Creates a plot with 'No Data' across the map if
    there is no data.

    Args:
        plotmap : (object) Map object created from emcpy
        domain : (object) Domain object created from emcpy
        metadata : (dict) metadata dictionary created from
                   PyGSI
    Returns:
        fig : (matplotlib figure) figure with 'No Data' displayed
    """
    # Titles
    labels = features.get_labels(metadata)
    plotmap.add_title(labels['title'], loc='left', fontsize=12)
    plotmap.add_title(label=labels['date title'],
                      loc='right', fontsize=12,
                      fontweight='semibold')

    # Get center of map location
    lon1 = domain.extent[0]+180
    lon2 = domain.extent[1]+180
    xloc = (lon2 - ((lon2-lon1)/2)) - 180

    lat1 = domain.extent[2]+90
    lat2 = domain.extent[3]+90
    yloc = (lat2 - (lat2-lat1)/2) - 90

    # Plot text
    plotmap.add_text(xloc, yloc, 'No Data', fontsize=32,
                     alpha=0.6, horizontalalignment='center')

    # Return figure
    fig = plotmap.return_figure()

    return fig
