#!/usr/bin/env python3

from matplotlib import pyplot
from pyGSI import GSIStat


adate = '2020010100'
fname = 'gdas.t00z.gsistat'

gdas = GSIstat(fname, adate)

ps = gdas.extract('ps')
print(ps)
uv = gdas.extract('uv')
print(uv)
oz = gdas.extract('oz')
print(oz)
gps = gdas.extract('gps')
print(gps)
rad = gdas.extract('rad')
print(rad)
cost = gdas.extract('cost')
print(cost)
cost.plot(y='gJ')

gdas.list_instruments()

amsua = gdas.extract_instrument('rad', 'amsua')
print(amsua)

pyplot.savefig('cost.png')
