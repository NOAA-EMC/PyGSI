import yaml
import numpy as np
from multiprocessing import Pool
from functools import partial
import itertools
import argparse
from LAMDA.map_qc_flags import map_qc_flags
from LAMDA.map_departures import map_departures
from LAMDA.layer_histogram import layer_histogram


def create_mp_work_list(diag_inputs, plotting_config):
    """
    Create working list of every possible combination of
    diag inputs and plotting features.
    """
    work = []
    for d in diag_inputs:
        keys, values = zip(*d.items())
        plotkeys, plotvals = zip(*plotting_config.items())

        keys = keys+plotkeys
        values = values+plotvals

        work_dicts = [dict(zip(keys, v)) for v in itertools.product(*values)]
        work.extend(work_dicts)

    return work


def plotting_func(config, diag_file, data_type='omf',
                  diag_type=None, outdir='./'):
    """
    Takes information required to create plots using
    multiprocessing.
    """

    config['diag file'] = diag_file
    config['data type'] = data_type
    config['diag type'] = diag_type
    config['outdir'] = outdir

    plot_dict = {
        'qcflags': map_qc_flags,
        'layer_histogram': layer_histogram,
        'map departures': map_departures,
    }

    plot_dict[config['plot type']](config)


def workflow(data_config, plotting_config, nprocs, outdir='./'):
    """
    Main workflow function for LAMDA diagnostic files.

    Args:
        data_config : (yaml) yaml file containing info about diag
                      files
        plotting_config : (yaml) yaml file containing info about
                          plotting information i.e. plot type,
                          domain, projection etc.
        nprocs : (int) number of processors to use for
                 multiprocessing
        outdir : (str) path to output diagnostics
    """

    # Open config yaml files
    with open(data_config, 'r') as file:
        config_yaml = yaml.load(file, Loader=yaml.FullLoader)

    with open(plotting_config, 'r') as file:
        plotting_yaml = yaml.load(file, Loader=yaml.FullLoader)

    diag_file = config_yaml['diagnostic']['path']
    data_type = config_yaml['diagnostic']['data type']

    for dtype in ['conventional', 'radiance', 'ozone']:
        if dtype in config_yaml['diagnostic'].keys():
            diag_type = dtype

    work = create_mp_work_list(config_yaml['diagnostic'][diag_type],
                               plotting_yaml['plotting'])

    # Create multiprocessing Pool
    p = Pool(processes=nprocs)
    p.map(partial(plotting_func, diag_file=diag_file,
                  data_type=data_type, diag_type=diag_type,
                  outdir=outdir), work)


if __name__ == '__main__':
    # Parse command line
    ap = argparse.ArgumentParser()
    ap.add_argument("-n", "--nprocs",
                    help="Number of tasks/processors for multiprocessing",
                    type=int, default=1)
    ap.add_argument("-d", "--data_yaml", required=True,
                    help="Path to yaml file with diag data")
    ap.add_argument("-p", "--plotting_yaml", required=True,
                    help="Path to yaml file with plotting info")
    ap.add_argument("-o", "--outdir", default='./',
                    help="Out directory where files will be saved")

    myargs = ap.parse_args()

    nprocs = myargs.nprocs
    data_yaml = myargs.data_yaml
    plotting_yaml = myargs.plotting_yaml
    outdir = myargs.outdir

    workflow(data_yaml, plotting_yaml, nprocs, outdir)
