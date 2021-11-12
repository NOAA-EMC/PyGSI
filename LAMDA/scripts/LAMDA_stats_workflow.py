import yaml
import glob
import os
import pandas as pd
import numpy as np
import itertools
from multiprocessing import Pool
from functools import partial
from datetime import datetime
from pyGSI.gsi_stat import GSIstat

# import the plotting scripts i.e:
# from LAMDA.obs_count import plot_obscount
# from LAMDA.minimization import plot_minimization


def concatenate_dfs(files, variable, cycles, data_type):
    """
    Reads in list of files and creates a concatenated
    dataframe of all cycles.

    Args:
        files : (list) list of files
        variable : (str) data obs type i.e. t, uv, q etc.
        cycles : (list) list of str cycles
        data_type : (str) the input data type i.e conventional
                    or radiance
    Returns:
        concatenated_df : concatenated dataframe over cycles
    """
    dfs = []
    for i, file in enumerate(files):
        gstat = GSIstat(file, cycles[i])

        if data_type == 'radiance':
            df = gstat.extract_instrument('rad', variable)
        else:
            df = gstat.extract(variable)

        dfs.append(df)

    concatenated_df = pd.concat(dfs)

    return concatenated_df


def plotting(config, data_dict, outdir, data_type, ob_type):
    """
    Main plotting function that gets all inputs for plotting scripts
    in order which includes fits_df, plotting_config, and outdir.

    Args:
        config : (list of tuples) the multiprocessing inputs
                 including the plot type to create and the data
                 inputs (either conventional (obsids, subtype, obuse)
                 or radiance (channel, obuse))
        data_dir : (dict) dictionary that includes fits data
        outdir : (str) path to where figures should be outputted
        data_type : (str) determines if data is conventional or
                    radiance
        ob_type : (str) the observation type
    """

    data_inputs = config[0]
    plot_type = config[-1]

    plotting_config = {}
    plotting_config['ob_type'] = ob_type

    if data_type == 'conventional':
        plotting_config['obsid'] = data_inputs[0]
        plotting_config['subtype'] = data_inputs[1]
        plotting_config['obuse'] = data_inputs[-1]

    elif data_type == 'radiance':
        # change ob_type to just the sensor to utilize
        # GSIstat extract_sensor()
        ob_type = ob_type.split('_')[0]
        plotting_config['channel'] = data_inputs[0]
        plotting_config['obuse'] = data_inputs[-1]

    # Loop through t-minus hours to grab proper data and
    # generate plots
    for tm in np.arange(7):
        fits_data = []
        cycles = []

        for cycle in data_dict.keys():
            fits_data.append(data_dict[cycle]['fits'][tm])
            cycles.append(cycle)

        # Concatenate all files into one dataframe
        fits_df = concatenate_dfs(fits_data, ob_type, cycles, data_type)

        plot_dict = {
            'obs count': plot_obscount,
        }

        plot_dict[plot_type](fits_df, plotting_config, outdir)


def create_minimization_plots(data_dict, outdir):
    """
    Since it is stand alone from other plotting scripts, this
    functions purpose is to generate the minimization plots.

    Args:
        data_dir : (dict) dictionary that includes fits2 data
        outdir : (str) path to where figures should be outputted
    """

    # Loop through t-minus hours to grab proper data and
    # generate plots
    for tm in np.arange(7):
        fits2_data = []
        cycles = []

        for cycle in data_dict.keys():
            fits2_data.append(data_dict[cycle]['fits2'][tm])
            cycles.append(cycle)

        # Concatenate all files into one dataframe
        fits2_df = concatenate_dfs(fits_data, 'cost', cycles,
                                   data_type='cost')

        # Create plot by calling plotting script
        minimization_plots(fits2_df, outdir)


def stats_workflow(config_yaml, nprocs, outdir):
    """
    Main function that reads input yaml and sorts
    GSI stat data correctly.

    Args:
        config_yaml : (yaml) input yaml that includes
                      path to stat files directory, variables
                      plot types, and out directory
        nprocs : (int) number of processors to use for
                 multiprocessing
        outdir : (str) path to output diagnostics
    """

    # Open config yaml files
    with open(config_yaml, 'r') as file:
        config_yaml = yaml.load(file, Loader=yaml.FullLoader)

    statdir = config_yaml['stat']['stat dir']
    data_type = config_yaml['stat']['data type']
    ob_type = config_yaml['stat']['ob type']
    plot_types = config_yaml['plot types']

    # Grabs all subdirectories of stats files (in date subdirectory)
    subdirs = sorted([f.path for f in os.scandir(statdir) if f.is_dir()])

    # Create dictionary of all stats files sorted by date
    data_dict = {}

    for subdir in subdirs:
        cycle = subdir.split('/')[-1]
        fits = glob.glob(subdir + '/*fits.*')
        fits2 = glob.glob(subdir + '/*fits2.*')

        data_dict[cycle] = {}
        data_dict[cycle]['fits'] = sorted(fits)
        data_dict[cycle]['fits2'] = sorted(fits2)

    # minimization plots are on their own so need to separate
    # remove from plot list and then call function to plot
    if 'minimization' in plot_types:
        plot_types.remove('minimization')
        create_minimization_plots(data_dict, outdir)

    # Create multiprocessing lists
    if data_type == 'conventional':
        conv_inputs = list(zip(config_yaml['stat']['observation id'],
                               config_yaml['stat']['observation subtype'],
                               config_yaml['stat']['obuse']))
        work_list = list(itertools.product(conv_inputs, plot_types))

    elif data_type == 'radiance':
        sat_inputs = list(zip(config_yaml['stat']['channels'],
                              config_yaml['stat']['obuse']))
        work_list = list(itertools.product(sat_inputs, plot_types))

    # Create multiprocessing Pool
    p = Pool(processes=nprocs)
    p.map(partial(plotting, data_dict=data_dict, outdir=outdir,
                  ob_type=ob_type, data_type=data_type), work_list)
