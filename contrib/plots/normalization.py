#!/usr/bin/python
import sys
import numpy as np
import PCAlifa as PCA
import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec

output_fmt = 'pdf'

fitsfile = sys.argv[1]

if len(sys.argv) > 2:
    output_dir = sys.argv[2]
else:
    output_dir = '../figuras'

print('Output directory: %s' % output_dir)

P = PCA.PCAlifa(fitsFile = fitsfile, quantilQFlag = 0.95, lc = [3800, 6850])
P.setStarlightMaskFile('/home/lacerda/CALIFA/Mask.mC')

P.PCA_obs()
P.tomograms_obs()

K = P.K

nPCs = 5
hfig = nPCs * 5

f = plt.figure(figsize = (15, 5), dpi = 100)
gs = gridspec.GridSpec(1, 2, width_ratios = [4, 7])
ax1 = plt.subplot(gs[0])
ax2 = plt.subplot(gs[1])
ax2.set_title('Fluxo de normaliz. por zona\n', fontsize=16)
fobs_norm__yx = K.zoneToYX(K.fobs_norm, extensive = False)
im = ax1.imshow(fobs_norm__yx, origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'hot_r')
f.colorbar(ax = ax1, mappable = im, use_gridspec = True)
ax2.plot(K.fobs_norm / 1.e-16, 'k')
ax2.set_xlim([0, K.N_zone])
ax2.set_ylabel(r'$F_{\lambda 5635}\ [10^{-16} erg/s/cm^2/\AA]$')
ax2.set_xlabel(r'zona')
ax2.grid()
f.tight_layout()
f.savefig('%s/K0277-fobs_norm.%s' % (output_dir, output_fmt))

nPCs = 5
hfig = nPCs * 5

f = plt.figure(figsize = (15, hfig))
gs = gridspec.GridSpec(nPCs, 2, width_ratios = [5, 8])

for i in range(nPCs):
    tn = i + 1
    gsi = 2 * i

    if i == 0 or i == 2:
        tomo_pc = P.tomo_obs__kyx[i, :, :]
        pc = P.eigVec_obs__lk[:, i]
    else:
        tomo_pc = -1. * P.tomo_obs__kyx[i, :, :]
        pc = -1. * P.eigVec_obs__lk[:, i]

    eval = P.eigVal_obs__k[i]
    evals = P.eigVal_obs__k
    f_obs_mean = P.ms_obs__l
    l = P.l_obs

    ax1 = plt.subplot(gs[gsi])
    ax2 = plt.subplot(gs[gsi + 1])
    ax1.set_title(r'tomograma $%02i$' % tn)

    im = ax1.imshow(tomo_pc, origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'hot_r')
    f.colorbar(ax = ax1, mappable = im, use_gridspec = True)
    eval_perc = 100. * eval / evals.sum()

    if i == 0:
        ax2.set_title(r'PCA com $F_{obs}$. - var: $%.2f\ \%%$' % eval_perc)
    else:
        ax2.set_title(r'var: $%.2f\ \%%$' % eval_perc)

    ax2.plot(l, pc, 'k')
    ax2.set_ylabel(r'$PC %02i$' % tn)
    plt.setp(ax2.get_xticklabels(), rotation = 45)
    ax2.grid()
    ax3 = ax2.twinx()
    ax3.plot(l, f_obs_mean, color = '0.65')
    ax3.set_ylabel(r'Espectro medio')
    ax2.xaxis.set_minor_locator(mpl.ticker.MaxNLocator(nbins = 35))

    for xmin in ax2.xaxis.get_minorticklocs():
        ax2.axvline(x = xmin, ls = ':', c = 'grey')

f.tight_layout()
f.savefig('%s/K0277-tomo-obs.%s' % (output_dir, output_fmt))
