import argparse
import yaml
import glob
import csv
import os
import numpy as np


def read_convinfo(infofile):
    # read in convinfo file
    obtype = []  # string used in the filename
    typeint = []  # integer value of observation type
    subtypeint = []  # integer value of observation subtype
    obuse = []  # use 1, monitor -1
    cv = open(infofile)
    rdcv = csv.reader(filter(lambda row: row[0] != '!', cv))
    for row in rdcv:
        try:
            rowsplit = row[0].split()
            obtype.append(rowsplit[0])
            typeint.append(int(rowsplit[1]))
            subtypeint.append(int(rowsplit[2]))
            obuse.append(int(rowsplit[3]))
        except IndexError:
            pass  # end of file
    cv.close()
    return obtype, typeint, subtypeint, obuse


def main(config):
    # call function to get lists from convinfo file
    obtype, typeint, subtypeint, obuse = read_convinfo(config['convinfo'])

    unique_obstypes = ['uv', 't', 'q', 'gps', 'ps']

    for ob in unique_obstypes:
        yamlout = {'stat': {}}

        obindx = np.where(np.array(obtype) == ob)

        obsid_list = np.array(typeint)[obindx].tolist()
        subtype_list = np.array(subtypeint)[obindx].tolist()
        obuse_list = np.array(obuse)[obindx].tolist()

        yamlout['stat']['stat dir'] = config['statdir']
        yamlout['stat']['ob type'] = ob
        yamlout['stat']['data type'] = 'conventional'
        yamlout['stat']['observation id'] = obsid_list
        yamlout['stat']['observation subtype'] = subtype_list
        yamlout['stat']['obuse'] = obuse_list

        filetype = f'conv_{ob}_stats'
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
