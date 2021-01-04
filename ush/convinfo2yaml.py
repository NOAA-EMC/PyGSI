#!/usr/bin/env python
# convinfo2yaml.py
# generate YAML from input convinfo file
# for specified input GSI diag file(s)
import argparse
import yaml
import glob
import csv


def main(config):
    # read in convinfo file
    obtype = []  # string used in the filename
    typeint = []  # integer value of observation type
    subtypeint = []  # integer value of observation subtype
    obuse = []  # use 1, monitor -1
    cv = open(config['convinfo'])
    rdcv = csv.reader(filter(lambda row: row[0] != '!', cv))
    for row in rdcv:
        try:
            rowsplit = row[0].split()
            obtype.append(rowsplit[0])
            typeint.append(rowsplit[1])
            subtypeint.append(rowsplit[2])
            obuse.append(rowsplit[3])
        except IndexError:
            pass  # end of file
    cv.close()

    # get list of conventional diagnostic files available
    diagpath = '%s/diag_conv_*_%s.%s.nc4' % (config['diagdir'],
                                             config['loop'], config['cycle'])
    diagfiles = glob.glob(diagpath)

    # initialize YAML dictionary for output
    yamlout = {'diagnostic': []}
    diagtype = 'O-A' if config['loop'] == 'anl' else 'O-F'
    figs = ['histogram', 'spatial']

    # loop through obtypes
    for i in range(len(obtype)):
        # first get filename and verify it exists
        diagfile = '%s/%s_%s_%s.%s.nc4' % ('diag_conv',
                                           config['diagdir'].rstrip('/'),
                                           obtype[i], config['loop'],
                                           config['cycle'])
        if diagfile not in diagfiles:
            continue  # skip if diag file is missing
        if int(obuse[i]) != 1:
            continue  # only process assimilated obs for now
        dictloop = {
                   'path': [diagfile],
                   'observation id': [int(typeint[i])],
                   'observation subtype': [int(subtypeint[i])],
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
args = parser.parse_args()

config = {}
config['diagdir'] = args.diagdir
config['cycle'] = args.cycle
config['convinfo'] = args.convinfo
config['yaml'] = args.yaml
config['loop'] = args.loop

main(config)
