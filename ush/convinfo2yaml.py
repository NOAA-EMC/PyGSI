#!/usr/bin/env python
# convinfo2yaml.py
# generate YAML from input convinfo file
# for specified input GSI diag file(s)
import argparse
import yaml
import glob
import csv


def read_convinfo(infofile):
    # read in convinfo file
    obtype = []  # string used in the filename
    typeint = []  # integer value of observation type
    subtypeint = []  # integer value of observation subtype
    obuse = []  # use 1, monitor -1
    cv = open(infofile)
    rdcv = csv.reader(filter(lambda row: row[0] != '!', cv))
    cv.close()
    for row in rdcv:
        try:
            rowsplit = row[0].split()
            obtype.append(int(rowsplit[0]))
            typeint.append(int(rowsplit[1]))
            subtypeint.append(int(rowsplit[2]))
            obuse.append(int(rowsplit[3]))
        except IndexError:
            pass  # end of file


def main(config):
    # call function to get lists from convinfo file
    obtype, typeint, subtypeint, obuse = read_convinfo(config['convinfo'])
    # get list of conventional diagnostic files available
    diagpath = (f"{config['diagdir']}/diag_conv_*_",
                f"{config['loop']}.{config['cycle']}.nc4")
    diagfiles = glob.glob(diagpath)

    # initialize YAML dictionary for output
    yamlout = {'diagnostic': []}
    if config['variable'] == 'obs':
        diagtype = 'observation'
    else if config['variable'] == 'hofx':
        diagtype = 'hofx'
    else:
        diagtype = 'O-A' if config['loop'] == 'anl' else 'O-F'
    figs = ['histogram', 'spatial']

    # loop through obtypes
    for iobtype, itype, isub, iuse in zip(obtype, otypeint, subtypeint, obuse):
        # first get filename and verify it exists
        diagfile = (f"{config['diagdir']}/",
                    f"diag_conv_{iobtype}_{config['loop']}",
                    f".{config['cycle']}.nc4")
        if diagfile not in diagfiles:
            continue  # skip if diag file is missing
        if iuse != 1:
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
                    help='ges|anl default ges', default='ges')
parser.add_argument('-v', '--variable', type=str,
                    help='read departures, obs, or H(x): omf | obs | hofx',
                    default='omf')
args = parser.parse_args()

config = {}
config['diagdir'] = args.diagdir.rstrip('/')
config['cycle'] = args.cycle
config['convinfo'] = args.convinfo
config['yaml'] = args.yaml
config['loop'] = args.loop
config['variable'] = args.variable

main(config)
