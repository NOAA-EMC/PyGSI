import yaml
import numpy as np
from multiprocessing import Pool
from LAMDA import map_qc_flags


def create_mp_work_list(diag_inputs, plotting_yaml):
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

    domain = config['domain']
    projection = config['projection']
    plot_type = config['plot type']

    if diag_type == 'conventional':
        obsid = config['observation id']
        subtype = config['observation subtype']
        analysis_use = config['analysis use']
        bias_correction = config['bias correction']

        plot_dict = {
            'qc flags': conv_map_qc_flags(
                diag_dir, obsid, subtype, analysis_use,
                domain, projection, plot_dir)
        }

        plot_dict[plot_type]

    elif diag_type == 'radiance':
        channel = config['channel']
        qc_flag = config['qc flag']
        analysis_use = config['analysis use']
        bias_correction = config['bias correction']

        plot_dict = {
            'qc flags': rad_map_qc_flags(
                diag_dir, channel, analysis_use,
                domain, projection, plot_dir)
        }

        plot_dict[plot_type]

    else:
        layer = config['layer']
        bias_correction = config['bias correction']

        plot_dict = {}

        plot_dict[plot_type]


def workflow(data_config, plotting_config, nprocs=1, plot_dir):

    # Open config yaml files
    with open(data_config, 'r') as file:
        config_yaml = yaml.load(file, Loader=yaml.FullLoader)

    with open(plotting_config, 'r') as file:
        plotting_yaml = yaml.load(file, Loader=yaml.FullLoader)

    inputfile = config_yaml['diagnostic']['path']
    datatype = config_yaml['diagnostic']['data type']

    for dtype in ['conventional', 'radiance', 'ozone']:
        if dtype in config_yaml['diagnostic'].keys():
            diag_type = dtype

    work = create_mp_work_list(config_yaml[diag_type],
                               plotting_yaml)

    # Create mulitprocessing Pool
    p = Pool(processes=nprocs)
    results = p.map(partial(plotting_func, diag_dir=diag_dir,
                            data_type=data_type, diag_type=diag_type,
                            plot_dir=plot_dir), work)
