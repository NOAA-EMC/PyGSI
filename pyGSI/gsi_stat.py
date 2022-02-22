"""
utility functions for GSI to read ASCII stat files
from Rahul.Mahajan
"""

__all__ = ['GSIstat']

import sys
import numpy as np
import pandas as pd
import re
from datetime import datetime


class GSIstat(object):
    """
    Object containing the GSI statistics
    """

    def __init__(self, filename, adate):
        """
        Initialize the GSIstat object
        INPUT:
            filename = filename of the gsistat file
            adate = analysis date
        OUTPUT:
            GSIstat: object containing the contents of the filename
        """

        self.filename = filename
        self.analysis_date = datetime.strptime(adate, '%Y%m%d%H')

        with open(self.filename, 'r', encoding='utf-8', errors='ignore') as fh:
            self._lines = fh.readlines()

        # Initialize cache for fast parsing
        self._cache = {}

        return

    def extract(self, name):
        """
        From the gsistat file, extract information:
        INPUT:
            name = information seeked
            Valid options are:
                ps, oz, uv, t, q, gps, rad, cost
        OUTPUT:
            df = dataframe containing information
        """

        # If name has already been parsed,
        # just return it from cache
        if name in self._cache:
            df = self._cache[name]
            return df

        if name in ['ps']:
            df = self._get_ps()
        elif name in ['oz']:
            df = self._get_ozone()
        elif name in ['uv', 't', 'q', 'gps']:
            df = self._get_conv(name)
        elif name in ['rad']:
            df = self._get_radiance()
        elif name in ['cost']:
            df = self._get_cost()
        else:
            raise IOError(f'option {name} is not defined')

        # Drop the o-g from the indices list
        if 'o-g' in list(df.index.names):
            df.reset_index(level='o-g', drop=True, inplace=True)

        # Add datetime index
        df = self._add_datetime_index(df)

        # Cache it for faster access
        self._cache[name] = df

        return df

    def list_instruments(self):
        """
        From the gsistat file, list all available instruments
        """

        for rtype in ['rad', 'oz']:
            otype = self.extract(rtype)
            instruments = sorted(
                otype.index.get_level_values('instrument').unique())
            satellites = sorted(
                otype.index.get_level_values('satellite').unique())

            print(f'Available {rtype} instruments:')
            print(' '.join(instruments))
            print()
            print(f'Available {rtype} satellites:')
            print(' '.join(satellites))
            print()

        return

    def extract_instrument(self, obtype, instrument):
        """
        From the gsistat file, extract detailed information on an instrument:
        INPUT:
            obtype     = observation type to extract (rad or oz)
            instrument = instrument name [must be in the observation type]
            E.g.:
                amsua, mhs, iasi, hirs, etc
        OUTPUT:
            df = dataframe containing information
        """

        # If instrument has already been parsed,
        # just return it from cache
        if instrument in self._cache:
            df = self._cache[instrument]
            return df

        # Ensure obtype is already called,
        # if not call it and cache it
        if obtype in list(self._cache.keys()):
            otype = self._cache[obtype]
        else:
            otype = self.extract(obtype)
            self._cache[obtype] = otype

        instruments = sorted(
            otype.index.get_level_values('instrument').unique())
        satellites = sorted(otype.index.get_level_values('satellite').unique())

        if instrument not in instruments:
            print('Instrument %s not found!' % instrument)
            print('%s contains ...' % self.filename)
            print(', '.join(str(x) for x in instruments))
            return None

        tmp = []
        pattern = r"\s+\d+\s+\d+\s+%s_\S+\s+\d+\s+\d+\s+" % instrument
        for line in self._lines:
            if re.match(pattern, line):
                tst = line.strip().split()
                tst = tst[:2] + tst[2].split('_') + tst[3:]
                if obtype == 'oz':
                    tsttmp = tst
                    tst = tsttmp[0:7] + tsttmp[8:]
                tmp.append(tst)

        columns = ['it', 'channel', 'instrument', 'satellite',
                   'nassim', 'nrej', 'oberr', 'OmF_bc', 'OmF_wobc',
                   'col1', 'col2', 'col3']
        df = pd.DataFrame(data=tmp, columns=columns)
        df.drop(['col1', 'col2', 'col3'], inplace=True, axis=1)
        df[['channel', 'nassim', 'nrej']] = df[[
            'channel', 'nassim', 'nrej']].astype(np.int)
        df[['oberr', 'OmF_bc', 'OmF_wobc']] = df[[
            'oberr', 'OmF_bc', 'OmF_wobc']].astype(np.float)

        # Since iteration number is not readily available, make one
        lendf = len(df)
        nouter = np.int(lendf / len(df['it'].unique()))
        douter = np.int(lendf / nouter)
        it = np.zeros(lendf, dtype=int)
        for i in range(nouter):
            its = douter * i
            ite = douter * (i + 1)
            it[its:ite] = i + 1
        df['it'] = it

        df = df[['it', 'instrument', 'satellite', 'channel',
                 'nassim', 'nrej', 'oberr', 'OmF_bc', 'OmF_wobc']]
        df.set_index(['it', 'instrument', 'satellite',
                      'channel'], inplace=True)

        # Add datetime index
        df = self._add_datetime_index(df)

        # Cache it for faster access
        self._cache[instrument] = df

        return df

    def _add_datetime_index(self, df):
        """
        Add the datetime as the first index
        INPUT:
            df = dataframe without datetime index
        OUTPUT:
            df = dataframe with datetime as the 1st index
        """

        # If date is already present, return
        if 'date' in list(df.index.names):
            return df

        indices = ['date'] + list(df.index.names)
        df['date'] = self.analysis_date
        df.set_index('date', append=True, inplace=True)
        df = df.reorder_levels(indices)

        return df

    # Surface pressure Fit
    def _get_ps(self):
        """
        Search for surface pressure
        """

        header = None
        pattern = (r'o-g\s+it\s+obs\s+use\s+typ\s+styp\s+count\s+bias\s'
                   r'+rms\s+cpen\s+qcpen')
        for line in self._lines:
            if re.search(pattern, line):
                header = line.strip()
                break

        if header is None:
            print('No matching header for PS')
            return None

        tmp = []
        pattern = r' o-g\s+\d{2}\s+ps\s'
        for line in self._lines:
            if re.match(pattern, line):
                # don't add monitored or rejected data
                # if any(x in line for x in ['mon', 'rej']):
                #    continue
                tmp.append(line.strip().split())

        columns = header.split()
        df = pd.DataFrame(data=tmp, columns=columns)
        df[['it', 'typ', 'count']] = df[['it', 'typ', 'count']].astype(np.int)
        df[['bias', 'rms', 'cpen', 'qcpen']] = df[[
            'bias', 'rms', 'cpen', 'qcpen']].astype(np.float)
        df.set_index(columns[:5], inplace=True)

        return df

    # Conventional Observation Fits
    def _get_conv(self, name):
        """
        Search for uv, t, q, or gps
        """

        # Get pressure levels
        ptops = []
        pattern = r'o-g\s+ptop\s'
        for line in self._lines:
            if re.search(pattern, line):
                header = line.strip()
                ptops = np.array(header.split()[2:-1], dtype=np.float)
                break
        if ptops is []:
            print(f'No matching ptop for {name}')
            sys.exit(1)

        header = None
        pbots = []
        pattern = r'o-g\s+it\s+obs\s+use\s+typ\s+styp\s+pbot\s'
        for line in self._lines:
            if re.search(pattern, line):
                header = line.strip()
                header = re.sub('pbot', 'stat', header)
                header = re.sub('0.200E' + r'\+04', 'column', header)
                pbots = np.array(header.split()[7:-1], dtype=np.float)
                break
        if pbots is []:
            print(f'No matching pbot for {name}')
            sys.exit(1)
        if header is None:
            print('No matching header for PS')
            sys.exit(1)
        tmp = []
        pattern = r' o-g\s+(\d\d)\s+%s\s' % name
        for line in self._lines:
            if re.match(pattern, line):
                # don't add monitored or rejected data
                # if any(x in line for x in ['mon', 'rej']):
                #    continue
                # don't add cpen or qcpen either
                # careful here, cpen here also removes qcpen
                # hence the extra space before qcpen and cpen
                if any(x in line for x in [' qcpen', ' cpen']):
                    continue
                tmp.append(line.strip().split())

        columns = header.split()
        df = pd.DataFrame(data=tmp, columns=columns)
        df[['it', 'typ']] = df[['it', 'typ']].astype(np.int)
        df.set_index(columns[:7], inplace=True)
        df = df.astype(np.float)

        return df

    # Ozone Fits
    def _get_ozone(self):
        """
        Search for ozone summary statistics
        """

        # Get header
        header = None
        pattern = r'it\s+sat\s+inst\s+'
        for line in self._lines:
            if re.search(pattern, line):
                header = re.sub('#', ' ', line)
                header = 'o-g ' + header.strip()
                break
        if header is None:
            print('No matching header for ozone')
            sys.exit(1)

        tmp = []
        pattern = r'o-g\s+(\d\d)\s+oz\s'
        for line in self._lines:
            if re.match(pattern, line):
                # don't add monitored or rejected data
                # if any(x in line for x in ['mon', 'rej']):
                #    continue
                line = re.sub('oz', ' ', line)
                tmp.append(line.strip().split())

        columns = header.split()
        df = pd.DataFrame(data=tmp, columns=columns)
        df[['it', 'read', 'keep', 'assim']] = df[[
            'it', 'read', 'keep', 'assim']].astype(np.int)
        df[['penalty', 'cpen', 'qcpen', 'qcfail']] = df[[
            'penalty', 'cpen', 'qcpen', 'qcfail']].astype(np.float)
        df.set_index(columns[:4], inplace=True)
        df = df.swaplevel('sat', 'inst')
        df.index.rename(['satellite', 'instrument'], level=[
                        'sat', 'inst'], inplace=True)

        return df

    # Radiances
    def _get_radiance(self):
        """
        Search for radiance summary statistics
        """

        # Get header
        header = None
        pattern = r'it\s+satellite\s+instrument\s+'
        for line in self._lines:
            if re.search(pattern, line):
                header = re.sub('#', ' ', line)
                header = 'o-g ' + header.strip()
                break
        if header is None:
            print('No matching header for radiance')
            sys.exit(1)

        tmp = []
        pattern = r'o-g (\d\d) %3s' % 'rad'
        for line in self._lines:
            if re.match(pattern, line):
                # don't add monitored or rejected data
                # if any(x in line for x in ['mon', 'rej']):
                #    continue
                line = re.sub('rad', ' ', line)
                tmp.append(line.strip().split())

        columns = header.split()
        df = pd.DataFrame(data=tmp, columns=columns)
        df[['it', 'read', 'keep', 'assim']] = df[[
            'it', 'read', 'keep', 'assim']].astype(np.int)
        df[['penalty', 'qcpnlty', 'cpen', 'qccpen']] = df[[
            'penalty', 'qcpnlty', 'cpen', 'qccpen']].astype(np.float)
        df.set_index(columns[:4], inplace=True)
        df = df.swaplevel('satellite', 'instrument')

        return df

    # Minimization
    def _get_cost(self):
        """
        Search for minimization and cost function information
        """

        # tmp = []
        # pattern = 'costterms Jb,Jo,Jc,Jl'
        # for line in self._lines:
        #    if re.match(pattern, line):
        #        tmp.append(line.strip().split('=')[-1].split())

        # if tmp is not []:
        #    columns = ['Outer', 'Inner', 'Jb', 'Jo', 'Jc', 'Jl']
        #    df = pd.DataFrame(data=tmp, columns=columns)
        #    df[['Outer', 'Inner', ]] = df[['Outer', 'Inner']].astype(np.int)
        #    df.set_index(columns[:2], inplace=True)
        #    df = df.astype(np.float)
        #    df['J'] = df.sum(axis=1)

        tmp = []
        pattern = 'cost,grad,step,b,step'
        for line in self._lines:
            if re.match(pattern, line):
                tmp.append(line.strip().split('=')[-1].split()[:-3])

        columns = ['Outer', 'Inner', 'J', 'gJ']
        df = pd.DataFrame(data=tmp, columns=columns)
        df[['Outer', 'Inner', ]] = df[['Outer', 'Inner']].astype(np.int)
        df.set_index(columns[:2], inplace=True)
        df = df.astype(np.float)

        return df
