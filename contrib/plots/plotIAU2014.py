#!/usr/bin/python
import numpy as np
import pcalifa
from pycasso import fitsQ3DataCube
import matplotlib as mpl
from matplotlib import pyplot as plt
from scipy import stats as st

mpl.rcParams['font.size']       = 24
mpl.rcParams['axes.labelsize']  = 24
mpl.rcParams['axes.titlesize']  = 26
mpl.rcParams['xtick.labelsize'] = 15
mpl.rcParams['ytick.labelsize'] = 15 
mpl.rcParams['font.family']     = 'serif'
mpl.rcParams['font.serif']      = 'Times New Roman'

gal = 'K0802'
#fUsed = 'f_obs_norm'
fUsed = 'f_syn_norm'

P = pcalifa.PCAlifa()
K = fitsQ3DataCube('/Users/lacerda/CALIFA/gal_fits/px1_q043.d14a/%s_synthesis_eBR_px1_q043.d14a512.ps03.k1.mE.CCM.Bgsd6e.fits' % gal)
P.readPCAlifaFits('/Users/lacerda/LOCAL/data/%s_synthesis_eBR_px1_q043.d14a512.ps03.k1.mE.CCM.Bgsd6e_%s_PCA_0.95_STMASK.fits' % (gal, fUsed))

prop = {
    'arr'   : [ K.at_flux__z, np.log10(K.aZ_flux__z / 0.019), K.A_V, K.v_0, K.v_d ],
    'arr__yx'   : [ K.at_flux__yx, np.log10(K.aZ_flux__yx / 0.019), K.A_V__yx, K.v_0__yx, K.v_d__yx ],
    'label' : [ r'$\langle \log\ t \rangle_L\ [yr]$', r'$\log\ \langle Z \rangle_L\ [Z_\odot]$', r'$A_V\ [mag]$', r'$v_\star\ [km/s]$', r'$\sigma_\star\ [km/s]$', r'eigenvector' ],
    'filePrefix' : [ 'at_flux',  'logaZ_flux' ,'A_V', 'v_0', 'v_d' ]
}

iTomo = 2
iProp = 2
#inv = 1
inv = -1

propLabel = prop['label'][iProp]
prop__z = prop['arr'][iProp]
prop__yx = prop['arr__yx'][iProp]
tomo__z = inv * P.tomo__zk[:,iTomo]
tomo__yx = inv * P.tomo__kyx[iTomo, :, :]
tomoLim = [-2, 4]
propLim = [0, 2.5] 

f, axArr = plt.subplots(1, 3)
f.set_size_inches(15, 5)
#tomograma
ax = axArr[0]
im = ax.imshow(tomo__yx, origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'hot', vmin = tomoLim[0], vmax = tomoLim[1])
#im = ax.imshow(tomo__yx, origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'hot')
f.colorbar(ax = ax, mappable = im)
ax.set_title(r'Tomogram %i' % (iTomo + 1))

ax = axArr[1]
im = ax.imshow(prop__yx, origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'hot', vmin = propLim[0], vmax = propLim[1])
#im = ax.imshow(prop__yx, origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'hot')
f.colorbar(ax = ax, mappable = im)
ax.set_title(prop['label'][iProp])

ax = axArr[2]
ax.scatter(prop__z, tomo__z, c = 'k', marker = 'o', s = 0.1)
rhoPearson, pvalPearson = st.pearsonr(prop__z, tomo__z)
rhoSpearman, pvalSpearman = st.spearmanr(prop__z, tomo__z)
txt = 'Rs: %.2f' % rhoSpearman
textbox = dict(boxstyle = 'round', facecolor = 'wheat', alpha = 0.)
ax.text(0.75, 0.95, txt,
        fontsize = 16,
        transform = ax.transAxes,
        verticalalignment = 'top',
        bbox = textbox)
ax.set_xlabel(propLabel)
ax.set_ylabel(r'PC%1i' % (iTomo + 1))
ax.set_xlim(propLim)
ax.set_ylim(tomoLim)

f.tight_layout()
plt.subplots_adjust(wspace = 0.20, top = 0.90, bottom = 0.15)
f.savefig('%s-%s-PC%d-%s.pdf' % (K.califaID, prop['filePrefix'][iProp], (iTomo + 1), fUsed))