#!/usr/bin/python
import sys
import numpy as np
import PCAlifa as PCA
import matplotlib as mpl
from matplotlib import pyplot as plt

fitsfile = sys.argv[1]
scinot = 1e-16
output_dir = '.'

if len(sys.argv) > 2:
    zone = np.int(sys.argv[2])
else:
    zone = 0

P = PCA.PCAlifa(fitsFile = fitsfile, quantilQFlag = 0.95, lc = [3800, 6850])
P.setStarlightMaskFile('/home/lacerda/CALIFA/Mask.mC')
mask = P.maskLambdaConstrains & P.maskQFlag
maskNoEmLin = P.maskLambdaConstrains & P.maskQFlag & P.maskEmLines

K = P.K

# for i in range(66, K.N_zone):
#    zone = i
f_obs = K.f_obs[:, zone] / scinot
f_obs_mask = np.ma.masked_array(f_obs, mask = ~mask)
f_obs_maskEmLin = np.ma.masked_array(f_obs, mask = ~maskNoEmLin)

outformat = 'pdf'

mf = K.f_obs[:, zone].mean()

mpl.rcParams['font.size'] = 16.0

f, (ax1, ax2, ax3) = plt.subplots(3, sharex = True, sharey = True)

f.set_size_inches((15, 10))
f.set_dpi(100)
ax1.plot(P.K.l_obs, f_obs, 'k')
ax2.plot(P.K.l_obs, f_obs_mask, 'k')
ax3.plot(P.K.l_obs, f_obs_maskEmLin, 'k')
ax2.set_ylabel(r'$F_{obs}\ [10^{-16} erg/s/cm^2/\AA]$')
ax3.set_xlabel(r'$\lambda\ [\AA]$')
ax1.set_title(r'Espectro zona %d - NGC 2916' % zone)
plt.setp(ax2.get_yticklabels(), visible = False)
plt.setp(ax3.get_yticklabels(), visible = False)
ax3.xaxis.set_minor_locator(mpl.ticker.MaxNLocator(nbins = 35))
ax2.yaxis.labelpad = 60

for xmin in ax3.xaxis.get_minorticklocs():
    ax1.axvline(x = xmin, ls = ':', c = 'grey')
    ax2.axvline(x = xmin, ls = ':', c = 'grey')
    ax3.axvline(x = xmin, ls = ':', c = 'grey')

ax1.grid()
ax2.grid()
ax3.grid()
f.subplots_adjust(hspace = 0)
f.savefig('%s/%s-constant_inital_mask-%d.%s' % (output_dir, K.califaID, zone, outformat))
