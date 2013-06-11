'''
Created on 10/06/2013

@author: lacerda
'''
import matplotlib
matplotlib.use('agg')

import sys
import numpy as np
import PCAlifa as PCA
import argparse as ap
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MaxNLocator

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
    parser.add_argument('--lc', '-l',
                        help = 'Lambda constrains',
                        metavar = 'LAMBDA',
                        type = int,
                        nargs = 2,
                        default = [3800, 6850])
    parser.add_argument('--tmax', '-t',
                        help = 'Number max of eigenvectors to plot',
                        metavar = 'INT',
                        type = int,
                        default = 20)
    parser.add_argument('--eigv', '-e',
                        help = 'number of eigenvectors used to rebuild the signal',
                        type = int,
                        nargs = '+',
                        default = [1, 2, 3, 4, 5, 6, 10, 20])
    parser.add_argument('--zones', '-z',
                        help = 'The program will rebuild the spectrum of ZONES zones + last zone.',
                        type = int,
                        nargs = '+',
                        default = [0, 10, 20, 100, 200])

    return parser.parse_args()

def parseArrArgs(args, lastZone):
    args.eigvArr = np.array(args.eigv)
    args.zones.append(lastZone)
    args.zonesArr = np.array(args.zones)

def zoneRebuildSpec(iZone, evRebArr, l, O, tomo, eVec, eVal, mean, nPref, resid = False):
    P = PCA.PCAlifa()

    i = 0
    nCols = 2
    nRows = len(evRebArr) / (1. * nCols)

    if nRows > int(nRows):
        nRows = nRows + 1

    nRows = int(nRows)

    fig = plt.figure(figsize = (19.8, 10.8))
    gs = gridspec.GridSpec(nRows, nCols, width_ratios = [1, 1], hspace = 0)

    for i, ne in enumerate(evRebArr):
        eVMask = np.zeros_like(eVal, dtype = np.bool)
        eVMask[:ne] = True

        I_reconstr__il, f_reconstr__il = P.rebuildSpectra(tomo[:, eVMask], eVec[:, eVMask], mean)

        ax = plt.subplot(gs[i])

        P.zoneRebuildSpecAxisPlot(ax, l, O[iZone, :], f_reconstr__il[iZone, :],
                                  eVal, eVec, eVMask, nPref, 7, resid)

        if (i < len(evRebArr) - 2):
            plt.setp(ax.get_xticklabels(), visible = False)
        elif (i == len(evRebArr) - 2):
            i = -1

        i = i + 2

        ax.set_ylabel('using %d eigvec' % ne)

    plt.suptitle('CALIFA ID: %s%sZone %04d' % (nPref[:5], ' ' * 60, iZone))
    plt.tight_layout(pad = 2., w_pad = 0., h_pad = 0.)
    fig.savefig('%szone_%04d.png' % (nPref, iZone))
    plt.close()

if __name__ == '__main__':
    args = parser_args()

    P = PCA.PCAlifa(args.califaID, args.fitsDir, args.rFL, args.lc)

    parseArrArgs(args, P.K.N_zone - 1)
    
    if args.rSEL:
        P.setStarlightMaskFile(args.rSEL)

    P.runPCA()

    npref_f_obs = '%s-f_obs_' % P.K.califaID
    npref_f_obs_norm = '%s-f_obs_norm_' % P.K.califaID
    npref_f_syn = '%s-f_syn_' % P.K.califaID
    npref_f_syn_norm = '%s-f_syn_norm_' % P.K.califaID
    npref_f_res = '%s-f_res_' % P.K.califaID
    npref_r_res_norm = '%s-f_res_norm_' % P.K.califaID
 
#########################################################################
############################## Tomograms ################################

    P.screeTestPlot(P.eigVal_obs__k, args.tmax, npref_f_obs)
    P.screeTestPlot(P.eigVal_obs_norm__k, args.tmax, npref_f_obs_norm)
    P.screeTestPlot(P.eigVal_syn__k, args.tmax, npref_f_syn)
    P.screeTestPlot(P.eigVal_syn_norm__k, args.tmax, npref_f_syn_norm)
    P.screeTestPlot(P.eigVal_res__k, args.tmax, npref_f_res)
    P.screeTestPlot(P.eigVal_res_norm__k, args.tmax, npref_r_res_norm)

    for ti in range(args.tmax):
        P.tomoPlot(P.tomo_obs__kyx, P.l_obs, P.eigVec_obs__lk, P.eigVal_obs__k, ti, npref_f_obs)
        P.tomoPlot(P.tomo_obs_norm__kyx, P.l_obs, P.eigVec_obs_norm__lk, P.eigVal_obs_norm__k, ti, npref_f_obs_norm)
        P.tomoPlot(P.tomo_syn__kyx, P.l_obs, P.eigVec_syn__lk, P.eigVal_syn__k, ti, npref_f_syn)
        P.tomoPlot(P.tomo_syn_norm__kyx, P.l_obs, P.eigVec_syn_norm__lk, P.eigVal_syn_norm__k, ti, npref_f_syn_norm)
        P.tomoPlot(P.tomo_res__kyx, P.l_obs, P.eigVec_res__lk, P.eigVal_res__k, ti, npref_f_res)
        P.tomoPlot(P.tomo_res_norm__kyx, P.l_obs, P.eigVec_res_norm__lk, P.eigVal_res_norm__k, ti, npref_r_res_norm)

#########################################################################
#########################################################################


#########################################################################
########################### Rebuild Spectra #############################

    if (P.K.N_zone > args.zonesArr[-2]):
        for nZone in args.zonesArr:
            zoneRebuildSpec(nZone, args.eigvArr, P.l_obs, P.f_obs__zl, P.tomo_obs__zk, P.eigVec_obs__lk, P.eigVal_obs__k, P.ms_obs__l, npref_f_obs)
            zoneRebuildSpec(nZone, args.eigvArr, P.l_obs, P.f_obs_norm__zl, P.tomo_obs_norm__zk, P.eigVec_obs_norm__lk, P.eigVal_obs_norm__k, P.ms_obs_norm__l, npref_f_obs_norm)
            zoneRebuildSpec(nZone, args.eigvArr, P.l_obs, P.f_syn__zl, P.tomo_syn__zk, P.eigVec_syn__lk, P.eigVal_syn__k, P.ms_syn__l, npref_f_syn)
            zoneRebuildSpec(nZone, args.eigvArr, P.l_obs, P.f_syn_norm__zl, P.tomo_syn_norm__zk, P.eigVec_syn_norm__lk, P.eigVal_syn_norm__k, P.ms_syn_norm__l, npref_f_syn_norm)
            zoneRebuildSpec(nZone, args.eigvArr, P.l_obs, P.f_res__zl, P.tomo_res__zk, P.eigVec_res__lk, P.eigVal_res__k, P.ms_res__l, npref_f_res, resid = True)
            zoneRebuildSpec(nZone, args.eigvArr, P.l_obs, P.f_res_norm__zl, P.tomo_res_norm__zk, P.eigVec_res_norm__lk, P.eigVal_res_norm__k, P.ms_res_norm__l, npref_r_res_norm, resid = True)
    else:
        print "%s N_Zone < %d" % (P.K.califaID, args.zonesArr[-2])
        
#########################################################################
#########################################################################


#########################################################################
############################ Correlations ###############################

    colArr = [
            P.K.at_flux__z,
            P.K.aZ_flux__z / 0.019,
            P.K.A_V,
            P.K.v_0,
            P.K.v_d
    ]

############################### OBS NORM ###############################     

    nRows = 10
    nCols = len(colArr) + 1
    f, axArr = plt.subplots(nRows, nCols)
    f.set_size_inches(19.8, 10.8)

    for i in range(nRows):
        axArr[i, 0].set_ylabel('PC%d' % i)

        for j in range(nCols)[:-1]:
            P.correlationAxisPlot(colArr[j], P.tomo_obs_norm__zk[:, i], axArr[i, j])

        axArr[i, nCols - 1].plot(P.l_obs, P.eigVec_obs_norm__lk[:, i])
        plt.setp(axArr[i, nCols - 1].get_yticklabels(), visible = False)

    f.subplots_adjust(hspace = 0.0)
    f.subplots_adjust(wspace = 0.05)

    plt.setp([a.get_xticklabels() for a in f.axes[:-nCols]], visible = False)
    plt.setp([a.get_yticklabels() for a in f.axes[::nCols]], visible = True)
    axArr[0, 0].set_title('at_flux')
    axArr[0, 1].set_title('aZ_flux')
    axArr[0, 2].set_title('A_V')
    axArr[0, 3].set_title('v_0')
    axArr[0, 4].set_title('v_d')
    axArr[0, 5].set_title('eigenvector')

    plt.suptitle('Correlations PC0 ... PC9 - OBS NORM')
    f.savefig('%s-corre_obs_norm_0-9.png' % P.K.califaID)
    plt.close()

############################### SYN NORM ###############################     
    f, axArr = plt.subplots(nRows, nCols)
    f.set_size_inches(19.8, 10.8)

    for i in range(nRows):
        axArr[i, 0].set_ylabel('PC%d' % i)

        for j in range(nCols)[:-1]:
            P.correlationAxisPlot(colArr[j], P.tomo_syn_norm__zk[:, i], axArr[i, j])

        axArr[i, nCols - 1].plot(P.l_obs, P.eigVec_syn_norm__lk[:, i])
        plt.setp(axArr[i, nCols - 1].get_yticklabels(), visible = False)

    f.subplots_adjust(hspace = 0.0)
    f.subplots_adjust(wspace = 0.05)

    plt.setp([a.get_xticklabels() for a in f.axes[:-nCols]], visible = False)
    plt.setp([a.get_yticklabels() for a in f.axes[::nCols]], visible = True)
    axArr[0, 0].set_title('at_flux')
    axArr[0, 1].set_title('aZ_flux')
    axArr[0, 2].set_title('A_V')
    axArr[0, 3].set_title('v_0')
    axArr[0, 4].set_title('v_d')
    axArr[0, 5].set_title('eigenvector')

    plt.suptitle('Correlations PC0 ... PC9 - SYN NORM')
    f.savefig('%s-corre_syn_norm_0-9.png' % P.K.califaID)
    plt.close()

############################### RES NORM ###############################     
    f, axArr = plt.subplots(nRows, nCols)
    f.set_size_inches(19.8, 10.8)

    for i in range(nRows):
        axArr[i, 0].set_ylabel('PC%d' % i)

        for j in range(nCols)[:-1]:
            P.correlationAxisPlot(colArr[j], P.tomo_res_norm__zk[:, i], axArr[i, j])

        axArr[i, nCols - 1].plot(P.l_obs, P.eigVec_res_norm__lk[:, i])
        plt.setp(axArr[i, nCols - 1].get_yticklabels(), visible = False)

    f.subplots_adjust(hspace = 0.0)
    f.subplots_adjust(wspace = 0.1)

    plt.setp([a.get_xticklabels() for a in f.axes[:-nCols]], visible = False)
    plt.setp([a.get_yticklabels() for a in f.axes[::nCols]], visible = True)
    axArr[0, 0].set_title('at_flux')
    axArr[0, 1].set_title('aZ_flux')
    axArr[0, 2].set_title('A_V')
    axArr[0, 3].set_title('v_0')
    axArr[0, 4].set_title('v_d')
    axArr[0, 5].set_title('eigenvector')

    plt.suptitle('Correlations PC0 ... PC9 - SYN NORM')
    f.savefig('%s-corre_res_norm_0-9.png' % P.K.califaID)
    plt.close()

#########################################################################
#########################################################################


#########################################################################
###################### Population Correlations ##########################

############################### OBS NORM ###############################
    arrInd_1 = range(0, 5)
    arrInd_2 = range(5, 17)
    arrInd_3 = range(17, 27)
    arrInd_4 = range(27, 36)
    arrInd_5 = range(36, len(P.K.ageBase))

    logt = np.log10(P.K.ageBase)
    maskpopx1 = (logt < 7.5)
    maskpopx2 = (logt < 8.5) & (logt >= 7.5)
    maskpopx3 = (logt < 9.5) & (logt >= 8.5)
    maskpopx4 = (logt >= 9.5)

    popxtot = P.K.popx.sum(axis = 1).sum(axis = 0)

    colArr = []

    colArr.append(np.tensordot((P.K.popx.sum(axis = 1))[maskpopx1, :], logt[maskpopx1], (0, 0)) / popxtot)
    colArr.append(np.tensordot((P.K.popx.sum(axis = 1))[maskpopx2, :], logt[maskpopx2], (0, 0)) / popxtot)
    colArr.append(np.tensordot((P.K.popx.sum(axis = 1))[maskpopx3, :], logt[maskpopx3], (0, 0)) / popxtot)
    colArr.append(np.tensordot((P.K.popx.sum(axis = 1))[maskpopx4, :], logt[maskpopx4], (0, 0)) / popxtot)

    nRows = 10
    nCols = len(colArr) + 1
    f, axArr = plt.subplots(nRows, nCols)
    f.set_size_inches(19.8, 10.8)

    for i in range(nRows):
        axArr[i, 0].set_ylabel('PC%d' % i)

        for j in range(nCols)[:-1]:
            P.correlationAxisPlot(colArr[j], P.tomo_obs_norm__zk[:, i], axArr[i, j])

        axArr[i, nCols - 1].plot(P.l_obs, P.eigVec_obs_norm__lk[:, i])
        plt.setp(axArr[i, nCols - 1].get_yticklabels(), visible = False)

    f.subplots_adjust(hspace = 0.0)
    f.subplots_adjust(wspace = 0.05)

    plt.setp([a.get_xticklabels() for a in f.axes[:-nCols]], visible = False)
    plt.setp([a.get_yticklabels() for a in f.axes[::nCols]], visible = True)
    axArr[0, 0].set_title(r'$10^6 to 10^{7.5}yr$')
    axArr[0, 1].set_title(r'$10^{7.5} to 10^{8.5}yr$')
    axArr[0, 2].set_title(r'$10^{8.5} to 10^{9.5}yr$')
    axArr[0, 3].set_title(r'$10^{9.5} to 10^{10.2}yr$')
    axArr[0, 4].set_title('eigenvector')

    plt.suptitle('Correlations PC0 ... PC9 - OBS NORM')
    f.savefig('%s-corre_obs_norm_popx_0-9.png' % P.K.califaID)
    plt.close()

############################### SYN NORM ###############################
    f, axArr = plt.subplots(nRows, nCols)
    f.set_size_inches(19.8, 10.8)

    for i in range(nRows):
        axArr[i, 0].set_ylabel('PC%d' % i)

        for j in range(nCols)[:-1]:
            P.correlationAxisPlot(colArr[j], P.tomo_syn_norm__zk[:, i], axArr[i, j])

        axArr[i, nCols - 1].plot(P.l_obs, P.eigVec_syn_norm__lk[:, i])
        plt.setp(axArr[i, nCols - 1].get_yticklabels(), visible = False)

    f.subplots_adjust(hspace = 0.0)
    f.subplots_adjust(wspace = 0.05)

    plt.setp([a.get_xticklabels() for a in f.axes[:-nCols]], visible = False)
    plt.setp([a.get_yticklabels() for a in f.axes[::nCols]], visible = True)
    axArr[0, 0].set_title(r'$10^6 to 10^{7.5}yr$')
    axArr[0, 1].set_title(r'$10^{7.5} to 10^{8.5}yr$')
    axArr[0, 2].set_title(r'$10^{8.5} to 10^{9.5}yr$')
    axArr[0, 3].set_title(r'$10^{9.5} to 10^{10.2}yr$')
    axArr[0, 4].set_title('eigenvector')

    plt.suptitle('Correlations PC0 ... PC9 - OBS NORM')
    f.savefig('%s-corre_syn_norm_popx_0-9.png' % P.K.califaID)
    plt.close()

############################### RES NORM ###############################
    f, axArr = plt.subplots(nRows, nCols)
    f.set_size_inches(19.8, 10.8)

    for i in range(nRows):
        axArr[i, 0].set_ylabel('PC%d' % i)

        for j in range(nCols)[:-1]:
            P.correlationAxisPlot(colArr[j], P.tomo_res_norm__zk[:, i], axArr[i, j])

        axArr[i, nCols - 1].plot(P.l_obs, P.eigVec_res_norm__lk[:, i])
        plt.setp(axArr[i, nCols - 1].get_yticklabels(), visible = False)

    f.subplots_adjust(hspace = 0.0)
    f.subplots_adjust(wspace = 0.05)

    plt.setp([a.get_xticklabels() for a in f.axes[:-nCols]], visible = False)
    plt.setp([a.get_yticklabels() for a in f.axes[::nCols]], visible = True)
    axArr[0, 0].set_title(r'$10^6 to 10^{7.5}yr$')
    axArr[0, 1].set_title(r'$10^{7.5} to 10^{8.5}yr$')
    axArr[0, 2].set_title(r'$10^{8.5} to 10^{9.5}yr$')
    axArr[0, 3].set_title(r'$10^{9.5} to 10^{10.2}yr$')
    axArr[0, 4].set_title('eigenvector')

    plt.suptitle('Correlations PC0 ... PC9 - OBS NORM')
    f.savefig('%s-corre_res_norm_popx_0-9.png' % P.K.califaID)
    plt.close()

#########################################################################
#########################################################################
