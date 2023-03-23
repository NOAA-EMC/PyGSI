import argparse
import numpy as np
import yaml
from multiprocessing import Pool
from functools import partial
import sys
from pyGSI.diags import Ozone
from pyGSI.plot_diags import plot_map, plot_histogram
from datetime import datetime

start_time = datetime.now()


def plotting(ozone_config, diag_file, data_type, plot_type, var_yaml, outdir):

    analysis_use = ozone_config['analysis use'][0]
    layer = ozone_config['layer'][0]
    bias_correction = ozone_config['bias correction'][0]

    diag = Ozone(diag_file)

    df = diag.get_data(analysis_use=analysis_use)
    metadata = diag.metadata
    metadata['Diag Type'] = data_type

    bias = 'adjusted' if bias_correction else 'unadjusted'
    column = f'{data_type}' if data_type in ['observation'] \
        else f'{data_type}_{bias}'

    if layer == 0:
        dict_key = 'column total'
    else:
        dict_key = list(df)[layer]
    metadata['Layer'] = dict_key

    if analysis_use:
        lats = {
            'assimilated': df[dict_key]['assimilated']['latitude'].to_numpy(),
            'monitored': df[dict_key]['monitored']['latitude'].to_numpy()
        }
        lons = {
            'assimilated': df[dict_key]['assimilated']['longitude'].to_numpy(),
            'monitored': df[dict_key]['monitored']['longitude'].to_numpy()
        }

        data = {
            'assimilated': df[dict_key]['assimilated'][column].to_numpy(),
            'monitored': df[dict_key]['monitored'][column].to_numpy()
        }

    else:
        lats = df[dict_key]['latitude'].to_numpy()
        lons = df[dict_key]['longitude'].to_numpy()

        data = df[dict_key][column].to_numpy()

    if np.isin('histogram', plot_type):
        plot_histogram(data, metadata, var_yaml, outdir)
    if np.isin('spatial', plot_type):
        plot_map(lats, lons, data, metadata, var_yaml, outdir)


if __name__ == '__main__':

    # Parse command line
    ap = argparse.ArgumentParser()
    ap.add_argument("-n", "--nprocs",
                    help="Number of tasks/processors for multiprocessing")
    ap.add_argument("-y", "--yaml",
                    help="Path to yaml file with diag data")
    ap.add_argument("-o", "--outdir",
                    help="Out directory where files will be saved")
    ap.add_argument("-v", "--varyaml", required=False,
                    help="Path to yaml file with specific variable information")

    myargs = ap.parse_args()

    if myargs.nprocs:
        nprocs = int(myargs.nprocs)
    else:
        nprocs = 1

    input_yaml = myargs.yaml
    outdir = myargs.outdir

    with open(input_yaml, 'r') as file:
        parsed_yaml_file = yaml.load(file, Loader=yaml.FullLoader)

    work = parsed_yaml_file['diagnostic']['ozone']
    data_type = parsed_yaml_file['diagnostic']['data type']
    data_path = parsed_yaml_file['diagnostic']['path']
    try:
        plot_type = parsed_yaml_file['diagnostic']['plot types']
    except KeyError:
        raise Exception("'plot types' key not included in input yaml. "
                        "Please add key 'plot types' to yaml and list "
                        "of the plot types you would like to create. "
                        "i.e. ['histogram', 'spatial']")

    variable_yaml = myargs.varyaml if myargs.varyaml else None

    p = Pool(processes=nprocs)
    p.map(partial(plotting, diag_file=data_path,
                  data_type=data_type,
                  plot_type=plot_type,
                  var_yaml=variable_yaml,
                  outdir=outdir), work)

    print(datetime.now() - start_time)
