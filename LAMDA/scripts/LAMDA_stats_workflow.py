import glob
import os
import pandas as pd
import numpy as np
import itertools
from multiprocessing import Pool
from functools import partial
from datetime import datetime
from pyGSI.gsi_stat import GSIstat
# Add imports of specific plotting scripts


def concatenate_dfs(files, variable, cycles):
    """
    Reads in list of files and creates a concatenated
    dataframe of all cycles.
    """
    dfs = []
    for i, file in enumerate(files):
        gdas = GSIstat(file, cycles[i])
        df = gdas.extract(variable)
        dfs.append(df)

    concatenated_df = pd.concat(dfs)

    return concatenated_df


def plotting_func(config, data_dict, outdir):

    variable = config[0]
    plot_type = config[-1]

    for tm in np.arange(6):
        fits_data = []
        fits2_data = []
        cycles = []

        for cycle in data_dict.keys():
            fits_data.append(data_dict[cycle]['fits'][tm])
            fits2_data.append(data_dict[cycle]['fits2'][tm])
            cycles.append(cycle)

        if plot_type == 'minimization':
            # Minimization plots
            fits2_df = concatenate_dfs(fits2_data, variable, cycles)
            minimization_plots(fits2_df, outdir)

        else:
            # Other stats plots
            fits_df = concatenate_dfs(fits_data, var, cycles)

            plot_dict = {
                'obs count': obs_count,
                'rmse bias': rmse_bias,
                'error standard dev': error_std_dev
            }

            plot_dict[plot_type](fits_df, outdir)


def stats_workflow(config_yaml, nprocs):
    """
    Main function that reads input yaml and sorts
    GSI stat data correctly.

    Args:
        config_yaml : (yaml) input yaml that includes
                      path to stat files directory, variables
                      plot types, and out directory
        nprocs : (int) number of processors to use for
                 multiprocessing
    """

    inpath = config_yaml['inpath']
    outdir = config_yaml['outdir']
    variables = config_yaml['variables']
    plot_types = config_yaml['plot types']

    # Grabs all subdirectories of stats files (in date subdirectory)
    subdirs = sorted([f.path for f in os.scandir(inpath) if f.is_dir()])

    # Create dictionary of all stats files sorted by date
    data_dict = {}

    for subdir in subdirs:
        cycle = subdir.split('/')[-1]

        fits = glob.glob(subdir + '/*fits.*')
        fits2 = glob.glob(subdir + '/*fits2.*')

        data_dict[cycle] = {}
        data_dict[cycle]['fits'] = sorted(fits)
        data_dict[cycle]['fits2'] = sorted(fits2)

    # Create work list
    if 'minimization' in plot_types:
        plot_types.remove('minimization')
        work_list = list(itertools.product(variables, plot_types))
        work_list.append(('cost', 'minimization'))

    else:
        work_list = list(itertools.product(variables, plot_types))

    # Create multiprocessing Pool
    p = Pool(processes=nprocs)
    p.map(partial(plotting_func, data_dict=data_dict,
                  outdir=outdir), work)


if __name__ == '__main__':
    # Parse command line
    ap = argparse.ArgumentParser()
    ap.add_argument("-n", "--nprocs",
                    help="Number of tasks/processors for multiprocessing")
    ap.add_argument("-y", "--yaml_input",
                    help="Path to yaml file")

    myargs = ap.parse_args()
    nprocs = int(myargs.nprocs) if myargs.nprocs else 1

    stats_workflow(myargs.input_yaml, nprocs)
