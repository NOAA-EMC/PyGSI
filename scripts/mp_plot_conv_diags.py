import argparse
import numpy as np
import yaml
from multiprocessing import Pool
from functools import partial
from datetime import datetime
from pyGSI.diags import Conventional
from pyGSI.plot_diags import plot_map, plot_histogram

start_time = datetime.now()


def _uv_plotting(df, lats, lons, data_type, analysis_use, bias_correction,
                 metadata, plot_type, outdir):

    bias = 'adjusted' if bias_correction else 'unadjusted'
    # Creates column name fo u and v data
    u = f'u_{data_type}' if data_type in ['observation'] \
        else f'u_{data_type}_{bias}'
    v = f'v_{data_type}' if data_type in ['observation'] \
        else f'v_{data_type}_{bias}'

    if analysis_use:
        data = {
            'u': {
                'assimilated': df['assimilated'][u].to_numpy(),
                'rejected': df['rejected'][u].to_numpy(),
                'monitored': df['monitored'][u].to_numpy()
            },
            'v': {
                'assimilated': df['assimilated'][v].to_numpy(),
                'rejected': df['rejected'][v].to_numpy(),
                'monitored': df['monitored'][v].to_numpy()
            }
        }
    else:
        data = {
                'u': df[u].to_numpy(),
                'v': df[v].to_numpy(),
            }

    for key in data.keys():
        metadata['Variable'] = key
        if np.isin('histogram', plot_type):
            plot_histogram(data[key], metadata, outdir)
        if np.isin('spatial', plot_type):
            plot_map(lats, lons, data[key], metadata, outdir)


def plotting(conv_config, diag_file, data_type, plot_type, outdir):

    obsid = conv_config['observation id']
    obsubtype = conv_config['observation subtype']
    analysis_use = conv_config['analysis use'][0]
    bias_correction = conv_config['bias correction'][0]

    diag = Conventional(diag_file)

    df = diag.get_data(obsid=obsid, subtype=obsubtype,
                       analysis_use=analysis_use)
    metadata = diag.metadata
    metadata['Diag Type'] = data_type

    bias = 'adjusted' if bias_correction else 'unadjusted'
    column = f'{data_type}' if data_type in ['observation'] \
        else f'{data_type}_{bias}'

    if analysis_use:
        lats = {
            'assimilated': df['assimilated']['latitude'].to_numpy(),
            'rejected': df['rejected']['latitude'].to_numpy(),
            'monitored': df['monitored']['latitude'].to_numpy()
        }
        lons = {
            'assimilated': df['assimilated']['longitude'].to_numpy(),
            'rejected': df['rejected']['longitude'].to_numpy(),
            'monitored': df['monitored']['longitude'].to_numpy()
        }

        if metadata['Obs Type'] == 'conv' and metadata['Variable'] == 'uv':
            _uv_plotting(df, lats, lons, data_type, analysis_use,
                         metadata, plot_type, outdir)

        else:
            data = {
                'assimilated': df['assimilated'][column].to_numpy(),
                'rejected': df['rejected'][column].to_numpy(),
                'monitored': df['monitored'][column].to_numpy()
            }

            if np.isin('histogram', plot_type):
                plot_histogram(data, metadata, outdir)
            if np.isin('spatial', plot_type):
                plot_map(lats, lons, data, metadata, outdir)

    else:
        lats = df['latitude'].to_numpy()
        lons = df['longitude'].to_numpy()

        if metadata['Obs Type'] == 'conv' and metadata['Variable'] == 'uv':
            _uv_plotting(df, lats, lons, data_type, analysis_use,
                         bias_correction, metadata, plot_type, outdir)

        else:
            data = df[column].to_numpy()

            if np.isin('histogram', plot_type):
                plot_histogram(data, metadata, outdir)
            if np.isin('spatial', plot_type):
                plot_map(lats, lons, data, metadata, outdir)


if __name__ == '__main__':

    # Parse command line
    ap = argparse.ArgumentParser()
    ap.add_argument("-n", "--nprocs",
                    help="Number of tasks/processors for multiprocessing")
    ap.add_argument("-y", "--yaml",
                    help="Path to yaml file with diag data")
    ap.add_argument("-o", "--outdir",
                    help="Out directory where files will be saved")

    myargs = ap.parse_args()

    if myargs.nprocs:
        nprocs = int(myargs.nprocs)
    else:
        nprocs = 1

    input_yaml = myargs.yaml
    outdir = myargs.outdir

    with open(input_yaml, 'r') as file:
        parsed_yaml_file = yaml.load(file, Loader=yaml.FullLoader)

    work = parsed_yaml_file['diagnostic']['conventional']
    data_type = parsed_yaml_file['diagnostic']['data type']
    data_path = parsed_yaml_file['diagnostic']['path']
    try:
        plot_type = parsed_yaml_file['diagnostic']['plot types']
    except KeyError:
        raise Exception("'plot types' key not included in input yaml. "
                        "Please add key 'plot types' to yaml and list "
                        "of the plot types you would like to create. "
                        "i.e. ['histogram', 'spatial']")

    p = Pool(processes=nprocs)
    p.map(partial(plotting, diag_file=data_path,
                  data_type=data_type,
                  plot_type=plot_type,
                  outdir=outdir), work)

    print(datetime.now() - start_time)
