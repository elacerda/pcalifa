#!/usr/bin/python
# -*- coding: utf-8 -*-
from os.path import expanduser
import numpy as np
import matplotlib as mpl
import pcalifa as PCA
from matplotlib import pyplot as plt
from scipy import stats as st
from parser_opt import *

args = parser_args()

print('Output directory: %s' % args.outputDir)

P = PCA.PCAlifa(fitsFile = args.fitsfile, quantilQFlag = 0.95, lc = args.lc)
P.setStarlightMaskFile('%s/CALIFA/Mask.mC' % expanduser('~'))

K = P.K

# xxx

nSpec = 5

mask = P.maskQFlag & P.maskEmLines & P.maskLambdaConstrains
mask2d = (mask[..., np.newaxis] & np.ones_like(K.f_obs, dtype = np.bool))
mf_obs__lz = np.ma.masked_array(K.f_obs, mask = ~mask2d)
mf_obs_sigma__l = mf_obs__lz.std(axis = 1)
mf_obs_norm__lz = np.ma.masked_array(K.f_obs / K.fobs_norm, mask = ~mask2d)
mf_obs_norm_sigma__l = mf_obs_norm__lz.std(axis = 1)
mf_syn_norm__lz = np.ma.masked_array(K.f_syn / K.fobs_norm, mask = ~mask2d)
mf_syn_norm_sigma__l = mf_syn_norm__lz.std(axis = 1)

# xxx

prop = {
    'arr'   : [ K.at_flux__yx, np.log10(K.aZ_flux__yx / 0.019), K.A_V__yx, K.v_0__yx, K.v_d__yx ],
    'label' : [ r'$\langle \log\ t \rangle_L\ [yr]$', r'$\log\ \langle Z \rangle_L\ [Z_\odot]$', r'$A_V\ [mag]$', r'$v_\star\ [km/s]$', r'$\sigma_\star\ [km/s]$' ],
    'name'  : [ 'at_flux', 'aZ_flux', 'AV', 'v0', 'vd' ]
}

f, axArr = plt.subplots(3, 3)
f.set_size_inches((15, 13))

for ax in f.axes:
    ax.set_axis_off()

p5n, p50n, p95n = np.percentile(mf_obs_norm__lz, [5, 50, 95], axis = 1)
p5n[~mask] = np.nan
p50n[~mask] = np.nan
p95n[~mask] = np.nan
mfo_max_norm__l = mf_obs_norm__lz.max(axis = 1)
mfo_min_norm__l = mf_obs_norm__lz.min(axis = 1)

ax = axArr[0, 0]
ax.set_axis_on()
ax.plot(K.l_obs, p95n, 'k-', lw = 0.5)
ax.fill_between(K.l_obs, mfo_max_norm__l, mfo_min_norm__l, edgecolor = 'gray', facecolor = 'lightgray')
ax.set_ylabel(r'$f_{obs}$')
ax.xaxis.set_minor_locator(mpl.ticker.MaxNLocator(nbins = 35))
ax.grid()

p5sn, p50sn, p95sn = np.percentile(mf_syn_norm__lz, [5, 50, 95], axis = 1)
p5sn[~mask] = np.nan
p50sn[~mask] = np.nan
p95sn[~mask] = np.nan
mfs_max_norm__l = mf_syn_norm__lz.max(axis = 1)
mfs_min_norm__l = mf_syn_norm__lz.min(axis = 1)

ax = axArr[0, 1]
ax.set_axis_on()
ax.plot(K.l_obs, p95sn, 'k-', lw = 0.5)
ax.fill_between(K.l_obs, mfs_max_norm__l, mfs_min_norm__l, edgecolor = 'gray', facecolor = 'lightgray')
ax.set_ylabel(r'$f_{syn}\ [10^{-16} erg/s/cm^2/\AA]$')
ax.xaxis.set_minor_locator(mpl.ticker.MaxNLocator(nbins = 35))
ax.grid()

galimg = plt.imread(args.galaxyImgFile)

ax = axArr[0, 2]
ax.set_axis_on()
plt.setp(ax.get_xticklabels(), visible = False)
plt.setp(ax.get_yticklabels(), visible = False)

ax.imshow(galimg)

ax = axArr[1, 0]
ax.set_axis_on()
fobs_norm__yx = K.zoneToYX(K.fobs_norm / (K.zoneArea_pix * 1.e-16), extensive = False)
ax.set_title(r'$\log\ F_{\lambda 5635}\ [10^{-16} erg/s/cm^2/\AA]$')
im = ax.imshow(np.log10(fobs_norm__yx), origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'hot_r')
f.colorbar(ax = ax, mappable = im, use_gridspec = True)

ax = axArr[1, 1]
ax.set_axis_on()
p_i = 0
ax.set_title(prop['label'][p_i])
im = ax.imshow(prop['arr'][p_i], origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'hot_r')
f.colorbar(ax = ax, mappable = im, use_gridspec = True)

ax = axArr[1, 2]
ax.set_axis_on()
p_i = 1
ax.set_title(prop['label'][p_i])
im = ax.imshow(prop['arr'][p_i], origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'hot_r')
f.colorbar(ax = ax, mappable = im, use_gridspec = True)

ax = axArr[2, 0]
ax.set_axis_on()
p_i = 2
ax.set_title(prop['label'][p_i])
im = ax.imshow(prop['arr'][p_i], origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'hot_r')
f.colorbar(ax = ax, mappable = im, use_gridspec = True)

ax = axArr[2, 1]
ax.set_axis_on()
p_i = 3
ax.set_title(prop['label'][p_i])
im = ax.imshow(prop['arr'][p_i], origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'hot_r')
f.colorbar(ax = ax, mappable = im, use_gridspec = True)

ax = axArr[2, 2]
ax.set_axis_on()
p_i = 4
ax.set_title(prop['label'][p_i])
prc = np.percentile(K.v_d, 98.)
print prc
im = ax.imshow(prop['arr'][p_i], origin = 'lower', interpolation = 'nearest', aspect = 'auto', cmap = 'hot_r', vmin = 0, vmax = prc)
f.colorbar(ax = ax, mappable = im, use_gridspec = True)

# plt.suptitle(r'%s - %s' % (K.galaxyName, K.califaID))
f.savefig('%s/%s-apresent.%s' % (args.outputDir, K.califaID, args.outputImgSuffix))
