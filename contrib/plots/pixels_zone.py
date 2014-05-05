#!/usr/bin/python
import sys
import numpy as np
from pycasso import fitsQ3DataCube

fitsfile = sys.argv[1]

K = fitsQ3DataCube(fitsfile)

total = K.zoneArea_pix.size
onepix = np.where(K.zoneArea_pix == 1)[0].size
morethanten = np.where(K.zoneArea_pix > 10)[0].size

onepixperc = onepix / np.double(total)
morethantenperc = morethanten / np.double(total)
z = np.double(K.masterListData['z'])

arc2rad = 0.0000048481368111
c = 3.e5
H0 = 73.

umpixpc = 1. * arc2rad * (z * c / H0) * 1.e6

print '%s & %s & $%d$ & $%d$ & $%.2f$ & $%d$ & $%.2f$ & $%.2f$ \\\\' % (K.galaxyName, K.califaID, total, onepix, onepixperc, morethanten, morethantenperc, umpixpc)
