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
    sensors, channel, obuse = read_satinfo(config['satinfo'])

    unique_sensors = np.unique(sensors)

    for sensor in unique_sensors:
        yamlout = {'stat': {}}

        sensorindx = np.where(np.array(sensors) == sensor)

        channel_list = np.array(channel)[sensorindx].tolist()
        obuse_list = np.array(obuse)[sensorindx].tolist()

        yamlout['stat']['stat dir'] = config['statdir']
        yamlout['stat']['sensor'] = sensor
        yamlout['stat']['data type'] = 'radiance'
        yamlout['stat']['channels'] = channel_list
        yamlout['stat']['obuse'] = obuse_list

        filetype = f'{sensor}_stats'
        output_yamlfile = config['yaml'] + filetype + '.yaml'

        # write out the YAML
        with open(output_yamlfile, 'w') as file:
            yaml.dump(yamlout, file, default_flow_style=False)
        print('YAML written to ' + output_yamlfile)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(
        description=('Given an input convinfo GSI file and path to ',
                     'GSI stat files, generate an output YAML file ',
                     'for LAMDA monitoring.'))
    parser.add_argument('-s', '--statdir', type=str,
                        help='path to GSI stat files',
                        required=True)
    parser.add_argument('-i', '--convinfo', type=str,
                        help='path to GSI convinfo file',
                        required=True)
    parser.add_argument('-y', '--yaml', type=str,
                        help='path to output YAML file',
                        required=True)
    args = parser.parse_args()

    config = vars(args)

    main(config)