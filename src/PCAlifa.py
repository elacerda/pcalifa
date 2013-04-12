'''
Created on 19/10/2012

@author: lacerda
'''

from pycasso.fitsdatacube import fitsQ3DataCube
from matplotlib.ticker import MaxNLocator
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
removeStarlightEmLines = False
flagLinesQuantil = 0.9

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
    ax2.xaxis.set_major_locator(MaxNLocator(20))
    ax2.grid()
    plt.tight_layout()
    fig.savefig('%stomo%02i.png' % (npref, tn))
    plt.close()

def removeStarlightMask(maskfile, l_obs):
    t = atpy.Table(maskfile = maskfile, type = 'starlight_mask')

    mask = (l_obs > t[0]['l_up'])

    for i in range(1, len(t)):
        if (t[i]['weight'] == 0.0):
            mask = mask & ((l_obs < t[i]['l_low']) | (l_obs > t[i]['l_up']))

    return mask

def screeTest(eigval, maxInd, npref):
    fig = plt.figure(figsize = (8, 6))
    eigval_norm = eigval / eigval.sum()
    plt.plot(eigval_norm[:maxInd], linestyle = '-', marker = '*')
    plt.ylim([0, eigval_norm[1] * 1.1])
    plt.xticks(range(maxInd))
    plt.title('%s scree test' % npref)
    plt.grid()
    fig.savefig('%sscree.png' % npref)
    plt.close()

if __name__ == '__main__':
    # removendo linhas flaggeds do califa
    histo = K.f_flag.sum(axis = 1) / K.N_zone
    maskHisto = (histo < flagLinesQuantil)
    mask = maskHisto & (K.l_obs > 3800) & (K.l_obs < 6850)

    if removeStarlightEmLines == True:
        maskEmLines = removeStarlightMask(maskfile, K.l_obs)
        mask = mask & maskEmLines

    l_obs = K.l_obs[mask]
    f_obs__zl = K.f_obs[mask].transpose()
    f_obs_norm__zl = (K.f_obs[mask] / K.fobs_norm).transpose()
    f_syn__zl = K.f_syn[mask].transpose()
    f_syn_norm__zl = (K.f_syn[mask] / K.fobs_norm).transpose()
    f_res__zl = (K.f_obs[mask] - K.f_syn[mask]).transpose()
    f_res_norm__zl = ((K.f_obs[mask] - K.f_syn[mask]) / K.fobs_norm).transpose()

#########################################################################
################################ f_obs ##################################
#########################################################################

    I_obs__zl, covMat_obs, ms_obs, eigval_obs__k, eigvec_obs__lk = PCA(f_obs__zl, K.N_zone, 0)
    t__zk, t__kyx = tomogram(I_obs__zl, eigvec_obs__lk)

    npref = '%s_f_obs_' % K.califaID

    screeTest(eigval_obs__k, tmax, npref)

    for tn in range(tmax):
        tomoPlot(tn, l_obs, t__kyx, eigvec_obs__lk, eigval_obs__k, npref)

########################################################################
########################### f_obs_norm #################################
########################################################################

    I_obs_norm__zl, covMat_obs_norm, ms_obs_norm, eigval_obs_norm__k, eigvec_obs_norm__lk = PCA(f_obs_norm__zl, K.N_zone, 0)
    t__zk, t__kyx = tomogram(I_obs_norm__zl, eigvec_obs_norm__lk)

    npref = '%s_f_obs_norm_' % K.califaID

    screeTest(eigval_obs_norm__k, tmax, npref)

    for tn in range(tmax):
        tomoPlot(tn, l_obs, t__kyx, eigvec_obs_norm__lk, eigval_obs_norm__k, npref)

########################################################################
############################# f_syn ####################################
########################################################################

    I_syn__zl, covMat_syn, ms_syn, eigval_syn__k, eigvec_syn__lk = PCA(f_syn__zl, K.N_zone, 0)
    t__zk, t__kyx = tomogram(I_syn__zl, eigvec_syn__lk)

    npref = '%s_f_syn_' % K.califaID

    screeTest(eigval_syn__k, tmax, npref)

    for tn in range(tmax):
        tomoPlot(tn, l_obs, t__kyx, eigvec_syn__lk, eigval_syn__k, npref)

########################################################################
########################### f_syn_norm #################################
########################################################################

    I_syn_norm__zl, covMat_syn_norm, ms_syn_norm, eigval_syn_norm__k, eigvec_syn_norm__lk = PCA(f_syn_norm__zl, K.N_zone, 0)
    t__zk, t__kyx = tomogram(I_syn_norm__zl, eigvec_syn_norm__lk)

    npref = '%s_f_syn_norm_' % K.califaID

    screeTest(eigval_syn_norm__k, tmax, npref)

    for tn in range(tmax):
        tomoPlot(tn, l_obs, t__kyx, eigvec_syn_norm__lk, eigval_syn_norm__k, npref)

########################################################################
########################### f_res ##################################
########################################################################

    I_res__zl, covMat_res, ms_res, eigval_res__k, eigvec_res__lk = PCA(f_res__zl, K.N_zone, 0)
    t__zk, t__kyx = tomogram(I_res__zl, eigvec_res__lk)

    npref = '%s_f_res_' % K.califaID

    screeTest(eigval_res__k, tmax, npref)

    for tn in range(tmax):
        tomoPlot(tn, l_obs, t__kyx, eigvec_res__lk, eigval_res__k, npref)

########################################################################
######################### f_res_norm ###############################
########################################################################

    I_res_norm__zl, covMat_res_norm, ms_res_norm, eigval_res_norm__k, eigvec_res_norm__lk = PCA(f_res_norm__zl, K.N_zone, 0)
    t__zk, t__kyx = tomogram(I_res_norm__zl, eigvec_res_norm__lk)

    npref = '%s_f_res_norm_' % K.califaID

    screeTest(eigval_res_norm__k, tmax, npref)

    for tn in range(tmax):
        tomoPlot(tn, l_obs, t__kyx, eigvec_res_norm__lk, eigval_res_norm__k, npref)

