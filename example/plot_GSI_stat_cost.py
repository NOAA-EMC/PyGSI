#!/usr/bin/env python3

from matplotlib import pyplot
from pyGSI.GSIStat import GSIstat


adate = '2020102606'
fname = 'gdas.t06z.gsistat'

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
print(cost.index)
#cost.plot(y='gJ',x=cost.index)
pyplot.plot(cost.values[:,1])

gdas.list_instruments()

amsua = gdas.extract_instrument('rad', 'amsua')
print(amsua)

pyplot.savefig('cost.png')
