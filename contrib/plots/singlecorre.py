#!/usr/bin/python
# -*- coding: utf-8 -*-
from os.path import expanduser
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
from scipy import stats as st
import pcalifa as PCA
from parser_opt import *

mpl.rcParams['font.size'] = 15

def plot_corr_axes(x, y, ax):
    rhoPearson, pvalPearson = st.pearsonr(x, y)
    rhoSpearman, pvalSpearman = st.spearmanr(x, y)
    txt = 's: %.2f' % rhoSpearman

    ax.scatter(x, y, marker = 'o', s = 0.1)

    textbox = dict(boxstyle = 'round', facecolor = 'wheat', alpha = 0.)

    ax.text(0.85, 0.10, txt,
            fontsize = 15,
            transform = ax.transAxes,
            verticalalignment = 'top',
            bbox = textbox)

    # plt.setp(ax.get_yticklabels(), visible = False)

args = parser_args()

print('Output directory: %s' % args.outputDir)

P = PCA.PCAlifa(fitsFile = args.fitsfile, quantilQFlag = 0.95, lc = args.lc)
P.setStarlightMaskFile('%s/CALIFA/Mask.mC' % expanduser('~'))

K = P.K

P.PCA_obs_norm()
P.tomograms()

prop = {
    'arr'   : [ K.at_flux__z, np.log10(K.aZ_flux__z / 0.019), K.A_V, K.v_0, K.v_d ],
    'label' : [ r'$\langle \log\ t \rangle_L\ [yr]$', r'$\log\ \langle Z \rangle_L\ [Z_\odot]$', r'$A_V\ [mag]$', r'$v_\star\ [km/s]$', r'$\sigma_\star\ [km/s]$', r'eigenvector' ],
}

iPC = 2
iProp = 3

f = plt.figure()
ax = plt.gca()
x = prop['arr'][iProp]
y = P.tomo__zk[:, iPC]
ax.set_ylim(y.min(), y.max())
ax.set_ylabel(r'PC%d' % (iPC + 1))
ax.set_xlabel(prop['label'][iProp])
plot_corr_axes(x, y, ax)

if prop['label'][iProp] == r'$\sigma_\star\ [km/s]$':
    ax.set_xlim(0, 230)

plt.tight_layout()

f.savefig('%s/%s-f_obs_norm_pc_%d_prop_%d.%s' % (args.outputDir, K.califaID, (iPC + 1), iProp, args.outputImgSuffix))
