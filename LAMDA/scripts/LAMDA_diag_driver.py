#!/usr/bin/env python3
import yaml

# import the plotting functions:
from LAMDA.minimization_plots import minimization_plots
from LAMDA.bias_rmse_timeseries import bias_rmse_timeseries
from LAMDA.bias_stddev_channel import bias_stddev_channel
from LAMDA.map_qc_flags import map_qc_flags
from LAMDA.map_departures import map_departures
from LAMDA.layer_histogram import layer_histogram

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