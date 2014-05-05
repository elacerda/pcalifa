#!/usr/bin/python
from os.path import expanduser
import numpy as np
import pcalifa as PCA
import matplotlib as mpl
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from parser_opt import *

args = parser_args()

print('Output directory: %s' % args.outputDir)

P = PCA.PCAlifa(fitsFile = args.fitsfile, quantilQFlag = 0.95, lc = args.lc)
P.setStarlightMaskFile('%s/CALIFA/Mask.mC' % expanduser('~'))

P.PCA_syn_norm()
P.tomograms()

K = P.K

# xxx

firstPC = args.firstTomo
nPCs = args.nTomo
hfig = nPCs * 5

f = plt.figure(figsize = (15, hfig))
gs = gridspec.GridSpec(nPCs, 2, width_ratios = [5, 8])

for i in range(nPCs):
    iTomo = firstPC + i
    tn = iTomo + 1
    gsi = 2 * iTomo

    if iTomo > 1:
        tomo_pc = P.tomo__kyx[iTomo, :, :]
        pc = P.eigVec__lk[:, iTomo]
    else:
        tomo_pc = -1. * P.tomo__kyx[iTomo, :, :]
        pc = -1. * P.eigVec__lk[:, iTomo]

    eval = P.eigVal__k[iTomo]
    evals = P.eigVal__k
    f_syn_mean = P.ms__l
    l = P.l_syn

    ax1 = plt.subplot(gs[gsi])
    ax2 = plt.subplot(gs[gsi + 1])
    ax1.set_title(r'tomograma $%02i$' % tn)

    im = ax1.imshow(tomo_pc, origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'hot_r')
    f.colorbar(ax = ax1, mappable = im, use_gridspec = True)
    eval_perc = 100. * eval / evals.sum()

    if i == 0:
        ax2.set_title(r'PCA com $f_{syn}$. - var: $%.2f\ \%%$' % eval_perc)
    else:
        ax2.set_title(r'var: $%.2f\ \%%$' % eval_perc)

    ax2.plot(l, pc, 'k')
    ax2.set_ylabel(r'$PC %02i$' % tn)
#    plt.setp(ax2.get_xticklabels(), rotation = 45)
    ax2.set_ylabel(r'$\lambda$ [\AA]')
    ax2.grid()
    ax3 = ax2.twinx()
    ax3.plot(l, f_syn_mean, color = '0.65')
    ax3.set_ylabel(r'Espectro medio')
    ax2.xaxis.set_minor_locator(mpl.ticker.MaxNLocator(nbins = 35))

    for xmin in ax2.xaxis.get_minorticklocs():
        ax2.axvline(x = xmin, ls = ':', c = 'grey')

f.tight_layout()
f.savefig('%s/%s-tomo-syn-norm.%s' % (args.outputDir, K.califaID, args.outputImgSuffix))
