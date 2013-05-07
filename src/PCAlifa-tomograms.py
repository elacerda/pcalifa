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
import argparse as ap


def parser_args():
    parser = ap.ArgumentParser(description = 'PCAlifa - correlations')
    parser.add_argument('--califaID', '-c',
                        help = 'Califa ID (ex: K0277)',
                        type = str,
                        default = 'K0277')
    parser.add_argument('--fitsDir', '-f',
                        help = 'Califa FITS directory',
                        metavar = 'DIR',
                        type = str,
                        default = '/home/lacerda/CALIFA')
    parser.add_argument('--rSEL', '-S',
                        help = 'Remove Starlight Emission Lines ',
                        metavar = 'MASK FILENAME',
                        type = str,
                        default = False)
    parser.add_argument('--rFL', '-Q',
                        help = 'Remove Flagged Lamdas',
                        metavar = 'QUANTIL',
                        type = float,
                        default = 0.9)
    parser.add_argument('--tmax', '-t',
                        help = 'Number max of eigenvectors to plot',
                        metavar = 'INT',
                        type = int,
                        default = 20)

    return parser.parse_args()


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
    args = parser_args()

    P = PCA.PCAlifa(califaID = args.califaID,
                    fitsDir = args.fitsDir,
                    flagLinesQuantil = args.rFL)

    if args.rSEL:
        P.removeStarlightEmLines(args.rSEL)

    P.runPCA()

#########################################################################
################################ f_obs ##################################
#########################################################################

    npref = '%s_f_obs_' % P.K.califaID

    screeTestPlot(P.eigVal_obs__k, args.tmax, npref)

    for tn in range(args.tmax):
        tomoPlot(tn, P.l_obs, P.tomo_obs__kyx, P.eigVec_obs__lk, P.eigVal_obs__k, npref)

########################################################################
########################### f_obs_norm #################################
########################################################################

    npref = '%s_f_obs_norm_' % P.K.califaID

    screeTestPlot(P.eigVal_obs_norm__k, args.tmax, npref)

    for tn in range(args.tmax):
        tomoPlot(tn, P.l_obs, P.tomo_obs_norm__kyx, P.eigVec_obs_norm__lk, P.eigVal_obs_norm__k, npref)

########################################################################
############################# f_syn ####################################
########################################################################

    npref = '%s_f_syn_' % P.K.califaID

    screeTestPlot(P.eigVal_syn__k, args.tmax, npref)

    for tn in range(args.tmax):
        tomoPlot(tn, P.l_obs, P.tomo_syn__kyx, P.eigVec_syn__lk, P.eigVal_syn__k, npref)

########################################################################
########################### f_syn_norm #################################
########################################################################

    npref = '%s_f_syn_norm_' % P.K.califaID

    screeTestPlot(P.eigVal_syn_norm__k, args.tmax, npref)

    for tn in range(args.tmax):
        tomoPlot(tn, P.l_obs, P.tomo_syn_norm__kyx, P.eigVec_syn_norm__lk, P.eigVal_syn_norm__k, npref)

########################################################################
############################### f_res ##################################
########################################################################

    npref = '%s_f_res_' % P.K.califaID

    screeTestPlot(P.eigVal_res__k, args.tmax, npref)

    for tn in range(args.tmax):
        tomoPlot(tn, P.l_obs, P.tomo_res__kyx, P.eigVec_res__lk, P.eigVal_res__k, npref)

########################################################################
############################# f_res_norm ###############################
########################################################################

    npref = '%s_f_res_norm_' % P.K.califaID

    screeTestPlot(P.eigVal_res_norm__k, args.tmax, npref)

    for tn in range(args.tmax):
        tomoPlot(tn, P.l_obs, P.tomo_res_norm__kyx, P.eigVec_res_norm__lk, P.eigVal_res_norm__k, npref)
