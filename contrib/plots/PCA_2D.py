#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
import numpy as np
import matplotlib.pyplot as plt
from pycasso import fitsQ3DataCube
import scipy.linalg as splin

fitsfile = sys.argv[1]
# fitsfile = '/home/lacerda/CALIFA/K0277/K0277_synthesis_eBR_v20_q036.d13c512.ps03.k2.mC.CCM.Bgsd61.fits'

K = fitsQ3DataCube(fitsfile)

l1 = K.f_obs[K.l_obs == 4464][0] / K.fobs_norm
l2 = K.f_obs[K.l_obs == 4470][0] / K.fobs_norm

f_obs = np.copy(np.asarray([l1, l2]).T)
f_obs_mean = f_obs.mean(axis = 0)
f_obs_dev = f_obs - f_obs_mean

covMat = np.cov(f_obs_dev, rowvar = 0)

w, e = splin.eigh(covMat)
S = np.argsort(w)[::-1]
wS = w[S]
eS = e[:, S]

f = plt.figure()
plt.scatter(x = f_obs_dev.T[0] , y = f_obs_dev.T[1], marker = 'o', s = 0.4, c = 'lightgray')

vectscale = 0.25

# PC1
PC1_x0 = 0.
PC1_y0 = 0.
PC1_x1 = vectscale * (PC1_x0 + eS[0, 0])
PC1_y1 = vectscale * (PC1_y0 + eS[1, 0])
plt.arrow(PC1_x0, PC1_y0, PC1_x1, PC1_y1, head_width = 0.01, head_length = 0.02, fc = 'k', ec = 'k')
plt.text(PC1_x1, PC1_y1 + 0.03, 'PC1', fontsize = 12, backgroundcolor = 'w', verticalalignment = 'center', horizontalalignment = 'center')

# PC2
PC2_x0 = 0.
PC2_y0 = 0.
PC2_x1 = vectscale * (PC2_x0 + eS[0, 1])
PC2_y1 = vectscale * (PC2_y0 + eS[1, 1])
plt.arrow(PC2_x0, PC2_y0, PC2_x1, PC2_y1, head_width = 0.01, head_length = 0.02, fc = 'k', ec = 'k')
plt.text(PC2_x1, PC2_y1 + 0.04, 'PC2', fontsize = 12, backgroundcolor = 'w', verticalalignment = 'center', horizontalalignment = 'center')

plt.axis('scaled')
plt.xlim([-0.4, 0.4])
plt.ylim([-0.4, 0.4])
plt.xlabel(r'$\lambda_1$')
plt.ylabel(r'$\lambda_2$')
plt.tight_layout()
f.savefig('PCA2D.pdf')
