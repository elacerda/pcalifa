'''
Created on 19/10/2012

@author: lacerda
'''
import matplotlib
matplotlib.use('agg')

import sys
import numpy as np
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MaxNLocator
from matplotlib import pyplot as plt
import PCAlifa as PCA

fitsDir = sys.argv[1]
califaID = sys.argv[2]
maskfile = '/home/lacerda/workspace/PCA/src/Mask.mC'
flagLinesQuantil = 0.9
remFlaggedLambdas = True
remStarlightEmLines = False
tmax = 20 # numero maximo de eigenvalues

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

def screeTestPlot(eigval, maxInd, npref):
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
    P = PCA.PCAlifa(califaID = califaID,
                    fitsDir = fitsDir,
                    flagLinesQuantil = flagLinesQuantil,
                    remFlaggedLambdas = remFlaggedLambdas,
                    runDefaultPCA = True)

    if (remStarlightEmLines == True):
        P.removeStarlightEmLines(maskfile)

#########################################################################
################################ f_obs ##################################
#########################################################################

    npref = '%s_f_obs_' % P.K.califaID

    screeTestPlot(P.eigval_obs__k, tmax, npref)

    for tn in range(tmax):
        tomoPlot(tn, P.l_obs, P.tomo_obs__kyx, P.eigvec_obs__lk, P.eigval_obs__k, npref)

########################################################################
########################### f_obs_norm #################################
########################################################################

    npref = '%s_f_obs_norm_' % P.K.califaID

    screeTestPlot(P.eigval_obs_norm__k, tmax, npref)

    for tn in range(tmax):
        tomoPlot(tn, P.l_obs, P.tomo_obs_norm__kyx, P.eigvec_obs_norm__lk, P.eigval_obs_norm__k, npref)

########################################################################
############################# f_syn ####################################
########################################################################

    npref = '%s_f_syn_' % P.K.califaID

    screeTestPlot(P.eigval_syn__k, tmax, npref)

    for tn in range(tmax):
        tomoPlot(tn, P.l_obs, P.tomo_syn__kyx, P.eigvec_syn__lk, P.eigval_syn__k, npref)

########################################################################
########################### f_syn_norm #################################
########################################################################

    npref = '%s_f_syn_norm_' % P.K.califaID

    screeTestPlot(P.eigval_syn_norm__k, tmax, npref)

    for tn in range(tmax):
        tomoPlot(tn, P.l_obs, P.tomo_syn_norm__kyx, P.eigvec_syn_norm__lk, P.eigval_syn_norm__k, npref)

########################################################################
############################### f_res ##################################
########################################################################

    npref = '%s_f_res_' % P.K.califaID

    screeTestPlot(P.eigval_res__k, tmax, npref)

    for tn in range(tmax):
        tomoPlot(tn, P.l_obs, P.tomo_res__kyx, P.eigvec_res__lk, P.eigval_res__k, npref)

########################################################################
############################# f_res_norm ###############################
########################################################################

    npref = '%s_f_res_norm_' % P.K.califaID

    screeTestPlot(P.eigval_res_norm__k, tmax, npref)

    for tn in range(tmax):
        tomoPlot(tn, P.l_obs, P.tomo_res_norm__kyx, P.eigvec_res_norm__lk, P.eigval_res_norm__k, npref)
