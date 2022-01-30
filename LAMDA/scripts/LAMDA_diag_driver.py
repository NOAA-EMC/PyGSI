#!/usr/bin/env python3
import argparse
import yaml
from multiprocessing import Pool
from functools import partial

# import the plotting functions:
from LAMDA.minimization_plots import minimization_plots
from LAMDA.bias_rmse_timeseries import bias_rmse_timeseries
from LAMDA.bias_stddev_channel import bias_stddev_channel
from LAMDA.map_qc_flags import map_qc_flags
from LAMDA.map_departures import map_departures
from LAMDA.layer_histogram import layer_histogram


diag_dict_stats = {
    'minimization': minimization_plots,
    'bias_rmse_timeseries': bias_rmse_timeseries,
    'bias_stddev_channel': bias_stddev_channel,
}
diag_dict_ncdiags = {
    'map_qc_flags': map_qc_flags,
    'map_departures': map_departures,
    'layer_histogram': layer_histogram,
}


def create_diag(diag_config, exp_config, ref_config, out_config):
    """
    Main function to create diagnostic figure based on
    input configuration dictionaries
    
    Args:
        diag_config : (dict) dictionary to specify configuration
                      for each diagnostic figure listed in input YAML
        exp_config : (dict) dictionary containing configuration
                     for the desired experiment 
        ref_config : (dict or False) dictionary containing configuration
                     for the specified reference experiment
        out_config : (dict) dictionary containing configuration
                     for the output of diagnostic figures
    """
    print(diag_config)


def LAMDA_diag_driver(nprocs, config_yaml):
    """
    Top-level driver function for processing
    diagnostic figures for regional FV3 GSI experiments

    Args:
        config_yaml : (yaml) path to input yaml
                      that includes all configuration
                      for producing diagnostic plots
        nprocs : (int) number of processors to use for
                 multiprocessing
    """
    # Get config dict from YAML file
    with open(config_yaml, 'r') as file:
        config = yaml.load(file, Loader=yaml.FullLoader)
    
    # pull out config sections into separate dicts
    exp_config = config['experiment']
    if 'reference' in config:
        ref_config = config['reference']
    else:
        ref_config = False
    out_config = config['output']
    diag_list = config['diagnostics']
    
    # create multiprocessing pool
    p = Pool(processes=nprocs)
    # run multiprocessing pool and call generic function
    p.map(partial(create_diag, exp_config=exp_config,
                  ref_config=ref_config, out_config=out_config),
          diag_list)


if __name__ == '__main__':
    # Parse command line
    ap = argparse.ArgumentParser()
    ap.add_argument("-n", "--nprocs",
                    help="Number of tasks/processors for multiprocessing",
                    type=int, default=1)
    ap.add_argument("-y", "--yaml", required=True,
                    help="Path to yaml file with configuration")

    myargs = ap.parse_args()

    nprocs = myargs.nprocs
    config_yaml = myargs.yaml

    LAMDA_diag_driver(nprocs, config_yaml)
