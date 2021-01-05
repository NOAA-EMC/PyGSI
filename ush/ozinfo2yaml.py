#!/usr/bin/env python
# ozinfo2yaml.py
# generate YAML from input ozinfo file
# for specified input GSI diag file(s)
import argparse
import yaml
import glob
import csv


def read_ozinfo(infofile):
    # read in ozinfo file
    sensor = []  # string used in the filename
    layer = []  # integer value of layer
    obuse = []  # use 1, monitor -1
    cv = open(infofile)
    rdcv = csv.reader(filter(lambda row: row[0] != '!', cv))
    for row in rdcv:
        try:
            rowsplit = row[0].split()
            sensor.append(rowsplit[0])
            layer.append(int(rowsplit[1]))
            obuse.append(int(rowsplit[2]))
        except IndexError:
            pass  # end of file
    cv.close()
    # loop through and set the last layer to 0, assuming total column
    # for sensors with multiple layers
    for i in range(len(layer)-1):
        if layer[i+1] <= layer[i]:
            layer[i] = 0  # indicates total column value
    return sensor, layer, obuse


def main(config):
    # call function to get lists from ozinfo file
    sensor, layer, obuse = read_ozinfo(config['ozinfo'])
    # get list of diagnostic files available
    diagpath = (f"{config['diagdir']}/diag_*_"
                f"{config['loop']}.{config['cycle']}.nc4")
    diagfiles = glob.glob(diagpath)

    # initialize YAML dictionary for output
    yamlout = {'diagnostic': []}
    if config['variable'] == 'obs':
        diagtype = 'observation'
    elif config['variable'] == 'hofx':
        diagtype = 'hofx'
    else:
        diagtype = 'O-A' if config['loop'] == 'anl' else 'O-F'
    figs = ['histogram', 'spatial']

    # loop through obtypes
    for isensor, iuse, ichan in zip(sensor, obuse, layer):
        # first get filename and verify it exists
        diagfile = (f"{config['diagdir']}/"
                    f"diag_{isensor}_{config['loop']}"
                    f".{config['cycle']}.nc4")
        if diagfile not in diagfiles:
            continue  # skip if diag file is missing
        if iuse != 1 and not config['monitor']:
            continue  # only process assimilated obs for now
        dictloop = {
                   'path': [diagfile],
                   'layer': [ichan],
                   'qc flag': [0],
                   'data type': [diagtype],
                   'plot type': figs,
                   }
        yamlout['diagnostic'].append({'ozone input': dictloop})

    # write out the YAML
    with open(config['yaml'], 'w') as file:
        yaml.dump(yamlout, file, default_flow_style=False)
    print('YAML written to '+config['yaml'])


parser = argparse.ArgumentParser(description=('Given an input ozinfo ',
                                              'GSI file and path to ',
                                              'GSI diags, generate an output ',
                                              'YAML file for use by PyGSI'))
parser.add_argument('-d', '--diagdir', type=str,
                    help='path to GSI netCDF diags', required=True)
parser.add_argument('-c', '--cycle', type=str,
                    help='cycle YYYYMMDDHH', required=True)
parser.add_argument('-i', '--ozinfo', type=str,
                    help='path to GSI ozinfo file', required=True)
parser.add_argument('-y', '--yaml', type=str,
                    help='path to output YAML file', required=True)
parser.add_argument('-l', '--loop', type=str,
                    help='ges|anl default ges', default='ges')
parser.add_argument('-v', '--variable', type=str,
                    help='read departures, obs, or H(x): omf | obs | hofx',
                    default='omf')
parser.add_argument('-m', '--monitor', action='store_true',
                    help='include monitored obs?')
args = parser.parse_args()

config = vars(args)
config['diagdir'] = config['diagdir'].rstrip('/')

main(config)
