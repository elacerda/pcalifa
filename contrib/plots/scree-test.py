#!/usr/bin/python
# -*- coding: utf-8 -*-
import numpy as np
import matplotlib as mpl
from matplotlib import pyplot as plt
import PCAlifa as PCA
import matplotlib.gridspec as gridspec
from parser_opt import *

args = parser_args()

print('Output directory: %s' % args.outputdir)

P = PCA.PCAlifa(fitsFile = args.fitsfile, quantilQFlag = 0.95, lc = args.lc)
P.setStarlightMaskFile('/home/lacerda/CALIFA/Mask.mC')

P.PCA_obs()
P.PCA_obs_norm()
P.PCA_syn_norm()

K = P.K

#xxx

def add_subplot_axes(ax,rect,axisbg='w'):
    fig = plt.gcf()
    box = ax.get_position()
    width = box.width
    height = box.height
    inax_position  = ax.transAxes.transform(rect[0:2])
    transFigure = fig.transFigure.inverted()
    infig_position = transFigure.transform(inax_position)    
    x = infig_position[0]
    y = infig_position[1]
    width *= rect[2]
    height *= rect[3]  # <= Typo was here
    subax = fig.add_axes([x,y,width,height],axisbg=axisbg)
    x_labelsize = subax.get_xticklabels()[0].get_size()
    y_labelsize = subax.get_yticklabels()[0].get_size()
    x_labelsize *= rect[2]**0.5
    y_labelsize *= rect[3]**0.5
    subax.xaxis.set_tick_params(labelsize=x_labelsize)
    subax.yaxis.set_tick_params(labelsize=y_labelsize)

    return subax

maxPCs = 15

f = plt.figure(figsize = (8, 5), dpi = 100)

PCs = np.linspace(1, maxPCs, maxPCs)

eigval_obs = 100. * P.eigVal_obs__k / P.eigVal_obs__k.sum()
eigval_obs_norm = 100. * P.eigVal_obs_norm__k / P.eigVal_obs_norm__k.sum()
eigval_syn_norm = 100. * P.eigVal_syn_norm__k / P.eigVal_syn_norm__k.sum()

plt.plot(PCs, eigval_obs[:maxPCs], 'k+--', label = '$F_{obs}$')
plt.plot(PCs, eigval_obs_norm[:maxPCs], 'k^-', label = '$f_{obs}$')
plt.plot(PCs, eigval_syn_norm[:maxPCs], 'k*-', label = '$f_{syn}$')
plt.legend()
plt.ylim([0, eigval_obs_norm[1] * 1.1])
plt.xlim([1, maxPCs])
plt.xticks(range(1,maxPCs + 1))
plt.title(r'%s - %s' % (K.galaxyName, K.califaID))
plt.xlabel(r'PC')
plt.ylabel(r'Var. [$\%%$]')
plt.grid()
f.savefig('%s/%s-screetest.%s' % (args.outputdir, K.califaID, args.outputimgsuffix))
