#!/usr/bin/env python
# satinfo2yaml.py
# generate YAMLs from input satinfo file
# for specified input GSI diag file(s)
import argparse
import yaml
import glob
import csv
import os
import numpy as np


def read_satinfo(infofile):
    # read in satinfo file
    sensor = []  # string used in the filename
    channel = []  # integer value of channel
    obuse = []  # use 1, monitor -1
    cv = open(infofile)
    rdcv = csv.reader(filter(lambda row: row[0] != '!', cv))
    for row in rdcv:
        try:
            rowsplit = row[0].split()
            sensor.append(rowsplit[0])
            channel.append(int(rowsplit[1]))
            obuse.append(int(rowsplit[2]))
        except IndexError:
            pass  # end of file
    cv.close()
    return sensor, channel, obuse


def main(config):
    # call function to get lists from satinfo file
    sensor, channel, obuse = read_satinfo(config['satinfo'])
    # get list of diagnostic files available
    diagpath = os.path.join(config['diagdir'], 'diag_*')
    diagfiles = glob.glob(diagpath)
    # compute suffix for files
    tmpsuffix = os.path.basename(diagfiles[0]).split('.')
    suffix = '.'.join((tmpsuffix[-2][10:], tmpsuffix[-1]))

    if config['variable'] == 'obs':
        diagtype = 'observation'
    elif config['variable'] == 'hofx':
        diagtype = 'hofx'
    else:
        diagtype = 'oma' if config['loop'] == 'anl' else 'omf'

    unique_sensors = np.unique(sensor)

    for isensor in unique_sensors:
        yamlout = {'diagnostic': {}}

        sensorindx = np.where(np.array(sensor) == isensor)

        # first get filename and verify it exists
        diagfile = os.path.join(f"{config['diagdir']}",
                                (f"diag_{isensor}_{config['loop']}"
                                 f".{config['cycle']}{suffix}"))
        if diagfile not in diagfiles:
            continue  # skip if diag file is missing

        yamlout['diagnostic']['path'] = diagfile
        yamlout['diagnostic']['data type'] = diagtype
        yamlout['diagnostic']['radiance'] = []

        # loop through obtypes
        for iuse, ichan in zip(np.array(obuse)[sensorindx],
                               np.array(channel)[sensorindx]):

            if iuse != 1 and not config['monitor']:
                continue  # only process assimilated obs for now

            dictloop = {
                'channel': [int(ichan)],
                'qc flag': [0],
                'analysis use': [True]
            }

            yamlout['diagnostic']['radiance'].append(dictloop)

        filetype = os.path.basename(diagfile).split('.')[0]
        output_yamlfile = config['yaml'] + filetype + '.yaml'

        # write out the YAML
        with open(output_yamlfile, 'w') as file:
            yaml.dump(yamlout, file, default_flow_style=False)
        print('YAML written to ' + output_yamlfile)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=(
        'Given an input satinfo ',
        'GSI file and path to ',
        'GSI diags, generate an output ',
        'YAML file for use by PyGSI'))
    parser.add_argument('-d', '--diagdir', type=str,
                        help='path to GSI netCDF diags', required=True)
    parser.add_argument('-c', '--cycle', type=str,
                        help='cycle YYYYMMDDHH', required=True)
    parser.add_argument('-i', '--satinfo', type=str,
                        help='path to GSI satinfo file', required=True)
    parser.add_argument('-y', '--yaml', type=str,
                        help='path to output YAML file', required=True)
    parser.add_argument('-l', '--loop', type=str,
                        help='guess or analysis?',
                        choices=['ges', 'anl'], default='ges',
                        required=False)
    parser.add_argument('-v', '--variable', type=str,
                        help='read departures, obs, or H(x)',
                        choices=['omf', 'obs', 'hofx'],
                        default='omf', required=False)
    parser.add_argument('-m', '--monitor', action='store_true',
                        help='include monitored obs?', required=False)
    args = parser.parse_args()

    config = vars(args)

    main(config)
