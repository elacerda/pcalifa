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


sqN_z = np.sqrt(P.K.N_zone)
a = sqN_z / len(P.l_obs)
ed__l = P.f_obs_norm__zl.std(axis = 0) * sqN_z
ed_syn__l = P.f_syn_norm__zl.std(axis = 0) * sqN_z
edk_mean = np.sqrt(np.abs(P.eigVal_obs_norm__k)).sum() * a
cd__l = np.abs(P.f_obs_norm__zl - P.f_obs_norm__zl.mean(axis = 0)).max(axis = 0)
cd_syn__l = np.abs(P.f_syn_norm__zl - P.f_syn_norm__zl.mean(axis = 0)).max(axis = 0)
k90_obs = ((np.cumsum(P.eigVal_obs_norm__k) / P.eigVal_obs_norm__k.sum()) <= 0.9).sum() + 1
k90_syn = ((np.cumsum(P.eigVal_syn_norm__k) / P.eigVal_syn_norm__k.sum()) <= 0.9).sum() + 1

print '%s & %.5f & %.5f & %.5f & %.5f' % (K.califaID, a, ed__l.mean(), edk_mean, cd__l.mean())
print '%s & %.5f & %.5f & %d & %d & %.5f & %.5f & %.5f & %.5f' % (K.califaID, 
                                                                  (100. * P.eigVal_obs_norm__k[:6] / P.eigVal_obs_norm__k.sum())[:6].sum(), 
                                                                  (100. * P.eigVal_syn_norm__k[:6] / P.eigVal_syn_norm__k.sum())[:6].sum(), 
                                                                  k90_obs, k90_syn, 
                                                                  P.eigVal_obs_norm__k.sum(), P.eigVal_syn_norm__k.sum(),
                                                                  P.eigVal_obs_norm__k.mean(), P.eigVal_syn_norm__k.mean(),
                                                                 )

f, axArr = plt.subplots(3, 1)
f.set_size_inches((7, 9))

ax = axArr[0]
ed_mean = ed__l.mean()
ax.plot(P.l_obs, ed__l, label=r'$f_{obs}$')
ax.plot(P.l_obs, ed_syn__l, label=r'$f_{syn}$')
#x.plot(P.l_obs, cd__l, label=r'$c_\lambda$')
#x.plot(P.l_obs, P.covMat_obs_norm__ll.diagonal(), label=r'$\sigma_\lambda{}^2$')
mask2d = (P.maskQFlag & P.maskEmLines & P.maskLambdaConstrains)[..., np.newaxis] & np.ones_like(K.f_obs, dtype = np.bool)
#ax.fill_between(K.l_obs,
#                np.ma.masked_array(K.f_obs / K.fobs_norm, mask = ~mask2d).max(axis = 1),
#                np.ma.masked_array(K.f_obs / K.fobs_norm, mask = ~mask2d).min(axis = 1),
#                edgecolor='gray', facecolor='lightgray')
ax.set_xlabel(r'$\lambda$')
#ax.set_ylabel(r'f_obs')
ax.set_ylabel(r'$d_\lambda')
ax.grid()
ax.legend()
ax.xaxis.set_minor_locator(mpl.ticker.MaxNLocator(nbins = 35))
ax.text(P.l_obs[1], ed_mean, r'%.2f' % ed_mean,
        fontsize = 12,
        backgroundcolor='w',
        verticalalignment = 'center',
        horizontalalignment = 'center')

for xmin in ax.xaxis.get_minorticklocs():
    ax.axvline(x = xmin, ls = ':', c = 'grey')

ax = axArr[1]
#ax.plot(P.l_obs, ed_syn__l, label=r'$d_\lambda$')
#ax.plot(P.l_obs, cd_syn__l, label=r'$c_\lambda$')
ax.plot(P.l_obs, cd__l, label=r'$f_{obs}$')
ax.plot(P.l_obs, cd_syn__l, label=r'$f_{syn}$')
#mask2d = (P.maskQFlag & P.maskEmLines & P.maskLambdaConstrains)[..., np.newaxis] & np.ones_like(K.f_obs, dtype = np.bool)
#ax.fill_between(K.l_obs,
#                np.ma.masked_array(K.f_syn / K.fobs_norm, mask = ~mask2d).max(axis = 1),
#                np.ma.masked_array(K.f_syn / K.fobs_norm, mask = ~mask2d).min(axis = 1),
#                edgecolor='gray', facecolor='lightgray')
ax.set_xlabel(r'$\lambda$')
ax.set_ylabel(r'$c_\lambda$')
#ax.set_ylabel(r'f_syn')
ax.grid()
ax.legend()
ax.xaxis.set_minor_locator(mpl.ticker.MaxNLocator(nbins = 35))
#ax.text(P.l_obs[1], ed_mean, r'%.2f' % ed_mean,
#        fontsize = 12,
#        backgroundcolor='w',
#        verticalalignment = 'center',
#        horizontalalignment = 'center')

for xmin in ax.xaxis.get_minorticklocs():
    ax.axvline(x = xmin, ls = ':', c = 'grey')

subaxPCs = 15

ax = axArr[2]
subpos = [0.7, 0.5, 0.26, 0.45]
ax.plot(P.eigVal_obs_norm__k)
ax.set_xlim([0, len(P.l_obs)])
ax.set_yscale('log')
ax.set_xlabel(r'autoespectro')
ax.set_ylabel(r'$\log\ \Lambda_k$')
subax = add_subplot_axes(ax, subpos)
subax.plot(P.eigVal_obs_norm__k, 'k+-')
subax.set_ylim([P.eigVal_obs_norm__k[0:subaxPCs].min(), P.eigVal_obs_norm__k[0:subaxPCs].max()])
subax.set_xlim([0, subaxPCs + 1])
subax.set_yscale('log')
f.suptitle('%s' % K.califaID)
f.savefig('%s/%s-dist.%s' % (args.outputdir, K.califaID, args.outputimgsuffix))
