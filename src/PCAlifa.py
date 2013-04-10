'''
Created on 19/10/2012

@author: lacerda
'''

from pycasso.fitsdatacube import fitsQ3DataCube
import numpy as np
from scipy import linalg
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
import pystarlight.io
import atpy

califaID = 'K0277'
#califaID = 'K0577'
fitsfile = '/home/lacerda/CALIFA/%s/%s_synthesis_eBR_v20_q027.d13c512.ps3b.k1.mC.CCM.Bgsd01.v01.fits' % (califaID, califaID)
maskfile = '/home/lacerda/workspace/PCA/src/Mask.mC'
K = fitsQ3DataCube(fitsfile)
tmax = 20 # numero maximo de eigenvalues
cutZeros__l = True

def PCA(arr, num, axis = -1, arrMean = None):
    if arrMean == None:
        arrMean = arr.mean(axis = axis)

    diff = arr - arrMean
    covMat = np.dot(diff.T, diff) / (num - 1.)
    w, e = linalg.eig(covMat)

    return diff, covMat, arrMean, np.real(w), np.real(e)

def tomogram(I, eigvec, extensive = True):
    t__zk = np.dot(I, eigvec)
    #t__kyx = K.zoneToYX(t__zk.T, extensive = False)
    t__kyx = K.zoneToYX(t__zk.T, extensive = extensive)

    return t__zk, t__kyx

def tomoPlot(tn, l, t, eigvec, eigval, npref):
    fig = plt.figure(figsize = (15, 5))
    gs = gridspec.GridSpec(1, 2, width_ratios = [4, 7])
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])
    ax1.set_title('tomogram %02i' % tn)
    im = ax1.imshow(t[tn, :, :], origin = 'lower', interpolation = 'nearest', aspect = 'auto')
    fig.colorbar(ax = ax1, mappable = im, use_gridspec = True)
    ax2.set_title('eigval %.4e' % eigval[tn])
    ax2.plot(l, eigvec[:, tn])
    plt.tight_layout()
    fig.savefig('%stomo%02i.png' % (npref, tn))
    plt.close()

def buildStarlightMask(maskfile, l_obs):
    t = atpy.Table(maskfile = maskfile, type = 'starlight_mask')

    mask = (l_obs > t[0]['l_up'])

    for i in range(1, len(t)):
        if (t[i]['flag'] == 0.0):
            mask = mask & ((l_obs < t[i]['l_low']) | (l_obs > t[i]['l_up']))

    return mask

def screeTest(eigval, max, npref):
    fig = plt.figure(figsize = (8, 6))
    eigval_norm = eigval / eigval.sum()
    plt.plot(eigval_norm[:max], linestyle = '-', marker = '*')
    plt.ylim([0, eigval_norm[1] * 1.1])
    plt.xticks(range(max))
    plt.title('%s scree test' % npref)
    plt.grid()
    fig.savefig('%sscree.png' % npref)
    plt.close()

#if __name__ == '__main__':
##########################################################################
################################# f_obs ##################################
##########################################################################
#
#    i = 0
#
#    if cutZeros__l is True:
#        while K.f_obs[i, :].max() == 0.:
#            i = i + 1
#
#    l_obs = K.l_obs[i:]
#    f_obs__zl = K.f_obs[i:, :].T
#    I_obs__zl, covMat_obs, ms_obs, eigval_obs__k, eigvec_obs__lk = PCA(f_obs__zl, K.N_zone, 0)
#    t__zk, t__kyx = tomogram(I_obs__zl, eigvec_obs__lk)
#
#    npref = '%s_f_obs_' % K.califaID
#
#    screeTest(eigval_obs__k, tmax, npref)
#
#    for tn in range(tmax):
#        tomoPlot(tn, l_obs, t__kyx, eigvec_obs__lk, eigval_obs__k, npref)
#
##########################################################################
############################# f_obs_norm #################################
##########################################################################
#
#    i = 0
#
#    if cutZeros__l is True:
#        while K.f_obs[i, :].max() == 0.:
#            i = i + 1
#
#    l_obs = K.l_obs[i:]
#    f_obs_norm__zl = (K.f_obs[i:, :] / K.fobs_norm).T
#    I_obs_norm__zl, covMat_obs_norm, ms_obs_norm, eigval_obs_norm__k, eigvec_obs_norm__lk = PCA(f_obs_norm__zl, K.N_zone, 0)
#    t__zk, t__kyx = tomogram(I_obs_norm__zl, eigvec_obs_norm__lk)
#
#    npref = '%s_f_obs_norm_' % K.califaID
#
#    screeTest(eigval_obs_norm__k, tmax, npref)
#
#    for tn in range(tmax):
#        tomoPlot(tn, l_obs, t__kyx, eigvec_obs_norm__lk, eigval_obs_norm__k, npref)
#
##########################################################################
############################### f_syn ####################################
##########################################################################
#
#    i = 0
#
#    if cutZeros__l is True:
#        while K.f_syn[i, :].max() == 0.:
#            i = i + 1
#
#    l_obs = K.l_obs[i:]
#    f_syn__zl = K.f_syn[i:, :].T
#    I_syn__zl, covMat_syn, ms_syn, eigval_syn__k, eigvec_syn__lk = PCA(f_syn__zl, K.N_zone, 0)
#    t__zk, t__kyx = tomogram(I_syn__zl, eigvec_syn__lk)
#
#    npref = '%s_f_syn_' % K.califaID
#
#    screeTest(eigval_syn__k, tmax, npref)
#
#    for tn in range(tmax):
#        tomoPlot(tn, l_obs, t__kyx, eigvec_syn__lk, eigval_syn__k, npref)
#
##########################################################################
############################# f_syn_norm #################################
##########################################################################
#
#    i = 0
#
#    if cutZeros__l is True:
#        while K.f_syn[i, :].max() == 0.:
#            i = i + 1
#
#    l_obs = K.l_obs[i:]
#    f_syn_norm__zl = (K.f_syn[i:, :] / K.fobs_norm).T
#    I_syn_norm__zl, covMat_syn_norm, ms_syn_norm, eigval_syn_norm__k, eigvec_syn_norm__lk = PCA(f_syn_norm__zl, K.N_zone, 0)
#    t__zk, t__kyx = tomogram(I_syn_norm__zl, eigvec_syn_norm__lk)
#
#    npref = '%s_f_syn_norm_' % K.califaID
#
#    screeTest(eigval_syn_norm__k, tmax, npref)
#
#    for tn in range(tmax):
#        tomoPlot(tn, l_obs, t__kyx, eigvec_syn_norm__lk, eigval_syn_norm__k, npref)
#
##########################################################################
############################# f_obs_syn ##################################
##########################################################################
#
#    i = 0
#
#    if cutZeros__l is True:
#        while K.f_obs[i, :].max() == 0.:
#            i = i + 1
#        j = 0
#        while K.f_syn[j, :].max() == 0.:
#            j = j + 1
#        i = max(i, j)
#
#    l_obs = K.l_obs[i:]
#    f_obs_syn__zl = (K.f_obs[i:, :] - K.f_syn[i:, :]).T
#    I_obs_syn__zl, covMat_obs_syn, ms_obs_syn, eigval_obs_syn__k, eigvec_obs_syn__lk = PCA(f_obs_syn__zl, K.N_zone, 0)
#    t__zk, t__kyx = tomogram(I_obs_syn__zl, eigvec_obs_syn__lk)
#
#    npref = '%s_f_obs_syn_' % K.califaID
#
#    screeTest(eigval_obs_syn__k, tmax, npref)
#
#    for tn in range(tmax):
#        tomoPlot(tn, l_obs, t__kyx, eigvec_obs_syn__lk, eigval_obs_syn__k, npref)
#
##########################################################################
########################### f_obs_syn_norm ###############################
##########################################################################
#
#    i = 0
#
#    if cutZeros__l is True:
#        while K.f_obs[i, :].max() == 0.:
#            i = i + 1
#        j = 0
#        while K.f_syn[j, :].max() == 0.:
#            j = j + 1
#        i = max(i, j)
#
#    f_obs_syn_norm__zl = ((K.f_obs[i:, :] - K.f_syn[i:, :]) / K.fobs_norm).T
#    I_obs_syn_norm__zl, covMat_obs_syn_norm, ms_obs_syn_norm, eigval_obs_syn_norm__k, eigvec_obs_syn_norm__lk = PCA(f_obs_syn_norm__zl, K.N_zone, 0)
#    t__zk, t__kyx = tomogram(I_obs_syn_norm__zl, eigvec_obs_syn_norm__lk)
#
#    npref = '%s_f_obs_syn_norm_' % K.califaID
#
#    screeTest(eigval_obs_syn_norm__k, tmax, npref)
#
#    for tn in range(tmax):
#        tomoPlot(tn, l_obs, t__kyx, eigvec_obs_syn_norm__lk, eigval_obs_syn_norm__k, npref)
