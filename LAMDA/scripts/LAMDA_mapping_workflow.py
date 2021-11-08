import yaml
import numpy as np
from multiprocessing import Pool
from functools import partial
import itertools
import argparse
import sys
sys.path.append('/scratch1/NCEPDEV/da/Kevin.Dougherty/PyGSI/')
from LAMDA.map_qc_flags import map_qc_flags


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


def plotting_func(config, diag_dir='./', data_type='omf',
                  diag_type=None, plot_dir='./'):
    """
    Takes information required to create plots using
    multiprocessing.
    """
    
    config['diag dir'] = diag_dir
    config['data type'] = data_type
    config['diag type'] = diag_type
    config['plot dir'] = plot_dir

    if diag_type == 'conventional':

        plot_dict = {
            'qcflags': map_qc_flags,
        }

        plot_dict[config['plot type']](config)

    elif diag_type == 'radiance':
        channel = config['channel']
        qc_flag = config['qc flag']
        analysis_use = config['analysis use']
        bias_correction = config['bias correction']

#         plot_dict = {
#             'qc flags': rad_map_qc_flags(
#                 diag_dir, channel, analysis_use,
#                 domain, projection, plot_dir)
#         }

#         plot_dict[plot_type]

    else:
        layer = config['layer']
        bias_correction = config['bias correction']

#         plot_dict = {}

#         plot_dict[plot_type]


def workflow(data_config, plotting_config, nprocs=1, plot_dir='./'):

    # Open config yaml files
    with open(data_config, 'r') as file:
        config_yaml = yaml.load(file, Loader=yaml.FullLoader)

    with open(plotting_config, 'r') as file:
        plotting_yaml = yaml.load(file, Loader=yaml.FullLoader)

    diag_dir = config_yaml['diagnostic']['path']
    data_type = config_yaml['diagnostic']['data type']

    for dtype in ['conventional', 'radiance', 'ozone']:
        if dtype in config_yaml['diagnostic'].keys():
            diag_type = dtype

    work = create_mp_work_list(config_yaml['diagnostic'][diag_type],
                               plotting_yaml['plotting'])
    
#     print(diag_dir)
#     plotting_func(work, diag_dir, data_type, diag_type, plot_dir)
    
    
    # Create multiprocessing Pool
    p = Pool(processes=nprocs)
    p.map(partial(plotting_func, diag_dir=diag_dir,
                  data_type=data_type, diag_type=diag_type,
                  plot_dir=plot_dir), work)


if __name__ == '__main__':
    # Parse command line
    ap = argparse.ArgumentParser()
    ap.add_argument("-n", "--nprocs",
                    help="Number of tasks/processors for multiprocessing")
    ap.add_argument("-d", "--data_yaml",
                    help="Path to yaml file with diag data")
    ap.add_argument("-p", "--plotting_yaml",
                    help="Path to yaml file with plotting info")
    ap.add_argument("-o", "--outdir",
                    help="Out directory where files will be saved")

    myargs = ap.parse_args()

    nprocs = int(myargs.nprocs) if myargs.nprocs else 1

    data_yaml = myargs.data_yaml
    plotting_yaml = myargs.plotting_yaml
    outdir = myargs.outdir

    workflow(data_yaml, plotting_yaml, nprocs, outdir)