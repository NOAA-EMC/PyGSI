#!/usr/bin/env python
# convinfo2yaml.py
# generate YAML from input convinfo file
# for specified input GSI diag file(s)
import argparse
import yaml
import glob
import csv
import os


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
    # get list of conventional diagnostic files available
    diagpath = os.path.join(config['diagdir'], 'diag_conv*')
    diagfiles = glob.glob(diagpath)
    # compute suffix for files
    tmpsuffix = os.path.basename(diagfiles[0]).split('.')
    suffix = '.'.join((tmpsuffix[-2][10:], tmpsuffix[-1]))

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
    for iobtype, itype, isub, iuse in zip(obtype, typeint, subtypeint, obuse):
        # first get filename and verify it exists
        diagfile = os.path.join(f"{config['diagdir']}",
                                (f"diag_conv_{iobtype}_{config['loop']}"
                                 f".{config['cycle']}{suffix}"))
        if diagfile not in diagfiles:
            continue  # skip if diag file is missing
        if iuse != 1 and not config['monitor']:
            continue  # only process assimilated obs for now
        dictloop = {
                   'path': [diagfile],
                   'observation id': [itype],
                   'observation subtype': [isub],
                   'analysis use': [True],
                   'data type': [diagtype],
                   'plot type': figs,
                   }
        yamlout['diagnostic'].append({'conventional input': dictloop})

    # write out the YAML
    with open(config['yaml'], 'w') as file:
        yaml.dump(yamlout, file, default_flow_style=False)
    print('YAML written to '+config['yaml'])


parser = argparse.ArgumentParser(description=('Given an input convinfo ',
                                              'GSI file and path to ',
                                              'GSI diags, generate an ',
                                              'output YAML file',
                                              'for use by PyGSI'))
parser.add_argument('-d', '--diagdir', type=str,
                    help='path to GSI netCDF diags', required=True)
parser.add_argument('-c', '--cycle', type=str,
                    help='cycle YYYYMMDDHH', required=True)
parser.add_argument('-i', '--convinfo', type=str,
                    help='path to GSI convinfo file', required=True)
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
