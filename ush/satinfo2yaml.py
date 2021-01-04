#!/usr/bin/env python
# satinfo2yaml.py
# generate YAML from input satinfo file
# for specified input GSI diag file(s)
import argparse
import yaml
import glob
import csv


def read_satinfo(infofile):
    # read in satinfo file
    sensor = []  # string used in the filename
    channel = []  # integer value of channel
    obuse = []  # use 1, monitor -1
    cv = open(infofile)
    rdcv = csv.reader(filter(lambda row: row[0] != '!', cv))
    cv.close()
    for row in rdcv:
        try:
            rowsplit = row[0].split()
            sensor.append(int(rowsplit[0]))
            channel.append(int(rowsplit[1]))
            obuse.append(int(rowsplit[2]))
        except IndexError:
            pass  # end of file
    return sensor, channel, obuse


def main(config):
    # call function to get lists from satinfo file
    sensor, channel, obuse = read_satinfo(config['satinfo'])
    # get list of diagnostic files available
    diagpath = '%s/diag_*_%s.%s.nc4' % (config['diagdir'], config['loop'],
                                        config['cycle'])
    diagpath = (f"{config['diagdir']}/diag_*_",
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
    for isensor, iuse, ichan in zip(sensor, obuse, channel):
        # first get filename and verify it exists
        diagfile = (f"{config['diagdir']}/",
                    f"diag_{isensor}_{config['loop']}",
                    f".{config['cycle']}.nc4")
        if diagfile not in diagfiles:
            continue  # skip if diag file is missing
        if iuse != 1:
            continue  # only process assimilated obs for now
        dictloop = {
                   'path': [diagfile],
                   'channel': [ichan],
                   'qc flag': [0],
                   'data type': [diagtype],
                   'plot type': figs,
                   }
        yamlout['diagnostic'].append({'radiance input': dictloop})

    # write out the YAML
    with open(config['yaml'], 'w') as file:
        yaml.dump(yamlout, file, default_flow_style=False)
    print('YAML written to '+config['yaml'])


parser = argparse.ArgumentParser(description=('Given an input satinfo ',
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
                    help='ges|anl default ges', default='ges')
parser.add_argument('-v', '--variable', type=str,
                    help='read departures, obs, or H(x): omf | obs | hofx',
                    default='omf')
args = parser.parse_args()

config = {}
config['diagdir'] = args.diagdir.rstrip('/')
config['cycle'] = args.cycle
config['satinfo'] = args.satinfo
config['yaml'] = args.yaml
config['loop'] = args.loop
config['variable'] = args.variable

main(config)
