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

plt.rcParams.update({'font.family' : 'Times New Roman',
                     'text.usetex' : False,
                     'backend' : 'ps'
                     })

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
    parser.add_argument('--rebuildSpec', '-R',
                        help = 'Plot the rebuilds --zones spectra.',
                        action = 'store_true',
                        default = False)
    parser.add_argument('--tomograms', '-T',
                        help = 'Plots tomograms.',
                        action = 'store_true',
                        default = False)
    parser.add_argument('--correlat', '-C',
                        help = 'Plots correlations.',
                        action = 'store_true',
                        default = False)
    parser.add_argument('--numcorrepc',
                        help = 'Number of PCs to plot correlations.',
                        type = int,
                        default = 5)
    parser.add_argument('--outputfname',
                        help = 'If provided, this adds OUTPUTFNAME to filenames.',
                        type = str,
                        default = False)


    return parser.parse_args()

def parseArrArgs(args, lastZone):
    args.eigvArr = np.array(args.eigv)
    args.zones.append(lastZone)
    args.zonesArr = np.array(args.zones)

def zoneRebuildSpec(iZone, evRebArr, l, O, tomo, eVec, eVal, mean, nPref, resid):
    P = PCA.PCAlifa()

    i = 0
    nCols = 2
    nRows = len(evRebArr) / (1. * nCols)

    if nRows > int(nRows):
        nRows = nRows + 1

    nRows = int(nRows)

    f, axArr = plt.subplots(nRows, nCols)
    f.set_size_inches(19.2, 10.8)

    k = 0

    for i in range(nRows):
        for j in range(nCols):
            ax = axArr[i, j]
            ne = evRebArr[k]
            ax.set_ylabel('using %d eigvec' % ne)

            eVMask = np.zeros_like(eVal, dtype = np.bool)
            eVMask[:ne] = True

            I_reconstr__il, f_reconstr__il = P.rebuildSpectra(tomo[:, eVMask], eVec[:, eVMask], mean)

            P.zoneRebuildSpecAxisPlot(ax, l, O[iZone, :], f_reconstr__il[iZone, :], eVal, eVec, eVMask, nPref, 10, resid)

            plt.setp(ax.get_yticklabels(), visible = False)
            k = k + 1

    f.subplots_adjust(hspace = 0.0)
    f.subplots_adjust(wspace = 0.15)
    plt.setp([a.get_xticklabels() for a in f.axes[:-nCols]], visible = False)
    plt.setp([a.get_yticklabels() for a in f.axes[::nCols]], visible = True)
    plt.suptitle(r'CALIFA ID: %s%sZone %04d' % (nPref[:5], ' ' * 60, iZone))
    fname = '%szone_%04d.png' % (nPref, iZone)
    f.savefig(fname)
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
    npref_f_res_norm = '%s-f_res_norm_' % P.K.califaID

    if args.outputfname:
        npref_f_obs = '%s%s_' % (npref_f_obs, args.outputfname)
        npref_f_obs_norm = '%s%s_' % (npref_f_obs_norm, args.outputfname)
        npref_f_syn = '%s%s_' % (npref_f_syn, args.outputfname)
        npref_f_syn_norm = '%s%s_' % (npref_f_syn_norm, args.outputfname)
        npref_f_res = '%s%s_' % (npref_f_res, args.outputfname)
        npref_f_res_norm = '%s%s_' % (npref_f_res_norm, args.outputfname)

    #### SANITY CHECK ####
    P.sanityCheck(P.l_obs, P.f_obs__zl[0, :], npref_f_obs)
    P.sanityCheck(P.l_obs, P.f_obs_norm__zl[0, :], npref_f_obs_norm)
    P.sanityCheck(P.l_syn, P.f_syn__zl[0, :], npref_f_syn)
    P.sanityCheck(P.l_syn, P.f_syn_norm__zl[0, :], npref_f_syn_norm)
    P.sanityCheck(P.l_obs, P.f_res__zl[0, :], npref_f_res)
    P.sanityCheck(P.l_obs, P.f_res_norm__zl[0, :], npref_f_res_norm)
    ######################

#########################################################################
############################## Tomograms ################################

    if args.tomograms:
        P.screeTestPlot(P.eigVal_obs__k, args.tmax, npref_f_obs, '%s FOBS' % P.K.califaID)
        P.screeTestPlot(P.eigVal_obs_norm__k, args.tmax, npref_f_obs_norm, '%s FOBS NORM' % P.K.califaID)
        P.screeTestPlot(P.eigVal_syn__k, args.tmax, npref_f_syn, '%s FSYN' % P.K.califaID)
        P.screeTestPlot(P.eigVal_syn_norm__k, args.tmax, npref_f_syn_norm, '%s FSYN NORM' % P.K.califaID)
        P.screeTestPlot(P.eigVal_res__k, args.tmax, npref_f_res, '%s RES' % P.K.califaID)
        P.screeTestPlot(P.eigVal_res_norm__k, args.tmax, npref_f_res_norm, '%s RES NORM' % P.K.califaID)

        for ti in range(args.tmax):
            P.tomoPlot(P.tomo_obs__kyx, P.l_obs, P.eigVec_obs__lk, P.eigVal_obs__k, P.ms_obs__l, ti, npref_f_obs)
            P.tomoPlot(P.tomo_obs_norm__kyx, P.l_obs, P.eigVec_obs_norm__lk, P.eigVal_obs_norm__k, P.ms_obs_norm__l, ti, npref_f_obs_norm)
            P.tomoPlot(P.tomo_syn__kyx, P.l_syn, P.eigVec_syn__lk, P.eigVal_syn__k, P.ms_syn__l, ti, npref_f_syn)
            P.tomoPlot(P.tomo_syn_norm__kyx, P.l_syn, P.eigVec_syn_norm__lk, P.eigVal_syn_norm__k, P.ms_syn_norm__l, ti, npref_f_syn_norm)
            P.tomoPlot(P.tomo_res__kyx, P.l_obs, P.eigVec_res__lk, P.eigVal_res__k, P.ms_res__l, ti, npref_f_res)
            P.tomoPlot(P.tomo_res_norm__kyx, P.l_obs, P.eigVec_res_norm__lk, P.eigVal_res_norm__k, P.ms_res_norm__l, ti, npref_f_res_norm)

#########################################################################
#########################################################################


#########################################################################
########################### Rebuild Spectra #############################
    if args.rebuildSpec:
        if (P.K.N_zone > args.zonesArr[-2]):
            for nZone in args.zonesArr:
                zoneRebuildSpec(nZone, args.eigvArr, P.l_obs,
                                P.f_obs__zl, P.tomo_obs__zk, P.eigVec_obs__lk, P.eigVal_obs__k, P.ms_obs__l, npref_f_obs, False)
                zoneRebuildSpec(nZone, args.eigvArr, P.l_obs,
                                P.f_obs_norm__zl, P.tomo_obs_norm__zk, P.eigVec_obs_norm__lk, P.eigVal_obs_norm__k, P.ms_obs_norm__l, npref_f_obs_norm, False)
                zoneRebuildSpec(nZone, args.eigvArr, P.l_syn,
                                P.f_syn__zl, P.tomo_syn__zk, P.eigVec_syn__lk, P.eigVal_syn__k, P.ms_syn__l, npref_f_syn, False)
                zoneRebuildSpec(nZone, args.eigvArr, P.l_syn,
                                P.f_syn_norm__zl, P.tomo_syn_norm__zk, P.eigVec_syn_norm__lk, P.eigVal_syn_norm__k, P.ms_syn_norm__l, npref_f_syn_norm, False)
                zoneRebuildSpec(nZone, args.eigvArr, P.l_obs,
                                P.f_res__zl, P.tomo_res__zk, P.eigVec_res__lk, P.eigVal_res__k, P.ms_res__l, npref_f_res, True)
                zoneRebuildSpec(nZone, args.eigvArr, P.l_obs,
                                P.f_res_norm__zl, P.tomo_res_norm__zk, P.eigVec_res_norm__lk, P.eigVal_res_norm__k, P.ms_res_norm__l, npref_f_res_norm, True)
        else:
            print "%s N_Zone < %d" % (P.K.califaID, args.zonesArr[-2])
#########################################################################
#########################################################################


#########################################################################
############################ Correlations ###############################
    if args.correlat:
        colArr = [
                P.K.at_flux__z,
                np.log10(P.K.aZ_flux__z / 0.019),
                P.K.A_V,
                P.K.v_0,
                P.K.v_d,
        ]

        colNames = [
            r'$\log\ t\ [yr]$',
            r'$\log\ Z\ [Z_\odot]$',
            r'$A_V\ [mag]$',
            r'$v_\star\ [km/s]$',
            r'$\sigma_\star\ [km/s]$',
            r'eigenvector',
        ]

        nRows = args.numcorrepc
        nCols = len(colArr) + 1

    ############################### OBS NORM ###############################     

        f, axArr = plt.subplots(nRows, nCols)
        f.set_size_inches(19.2, 10.8)

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

        for i in range(nCols):
            axArr[0, i].set_title(colNames[i])

        plt.suptitle(r'Correlations - OBS NORM')
        f.savefig('%s-corre_obs_norm.png' % P.K.califaID)
        plt.close()

    ############################### SYN NORM ###############################     
        f, axArr = plt.subplots(nRows, nCols)
        f.set_size_inches(19.2, 10.8)

        for i in range(nRows):
            axArr[i, 0].set_ylabel('PC%d' % i)

            for j in range(nCols)[:-1]:
                P.correlationAxisPlot(colArr[j], P.tomo_syn_norm__zk[:, i], axArr[i, j])

            axArr[i, nCols - 1].plot(P.l_syn, P.eigVec_syn_norm__lk[:, i])
            plt.setp(axArr[i, nCols - 1].get_yticklabels(), visible = False)

        f.subplots_adjust(hspace = 0.0)
        f.subplots_adjust(wspace = 0.05)

        plt.setp([a.get_xticklabels() for a in f.axes[:-nCols]], visible = False)
        plt.setp([a.get_yticklabels() for a in f.axes[::nCols]], visible = True)

        for i in range(nCols):
            axArr[0, i].set_title(colNames[i])

        plt.suptitle(r'Correlations - SYN NORM')
        f.savefig('%s-corre_syn_norm.png' % P.K.califaID)
        plt.close()

    ############################### RES NORM ###############################     
        f, axArr = plt.subplots(nRows, nCols)
        f.set_size_inches(19.2, 10.8)

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

        for i in range(nCols):
            axArr[0, i].set_title(colNames[i])

        plt.suptitle(r'Correlations - RES NORM')
        f.savefig('%s-corre_res_norm.png' % P.K.califaID)
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

        colNames = [
            r'$10^6\ to\ 10^{7.5}yr$',
            r'$10^{7.5}\ to\ 10^{8.5}yr$',
            r'$10^{8.5}\ to\ 10^{9.5}yr$',
            r'$10^{9.5}\ to\ 10^{10.2}yr$',
            r'eigenvector',
        ]

        nRows = args.numcorrepc
        nCols = len(colArr) + 1
        f, axArr = plt.subplots(nRows, nCols)
        f.set_size_inches(19.2, 10.8)

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

        for i in range(nCols):
            axArr[0, i].set_title(colNames[i])

        plt.suptitle(r'Correlations - OBS NORM')
        f.savefig('%s-corre_obs_norm_popx.png' % P.K.califaID)
        plt.close()

    ############################### SYN NORM ###############################
        f, axArr = plt.subplots(nRows, nCols)
        f.set_size_inches(19.2, 10.8)

        for i in range(nRows):
            axArr[i, 0].set_ylabel('PC%d' % i)

            for j in range(nCols)[:-1]:
                P.correlationAxisPlot(colArr[j], P.tomo_syn_norm__zk[:, i], axArr[i, j])

            axArr[i, nCols - 1].plot(P.l_syn, P.eigVec_syn_norm__lk[:, i])
            plt.setp(axArr[i, nCols - 1].get_yticklabels(), visible = False)

        f.subplots_adjust(hspace = 0.0)
        f.subplots_adjust(wspace = 0.05)

        plt.setp([a.get_xticklabels() for a in f.axes[:-nCols]], visible = False)
        plt.setp([a.get_yticklabels() for a in f.axes[::nCols]], visible = True)

        for i in range(nCols):
            axArr[0, i].set_title(colNames[i])

        plt.suptitle(r'Correlations - SYN NORM')
        f.savefig('%s-corre_syn_norm_popx.png' % P.K.califaID)
        plt.close()

    ############################### RES NORM ###############################
        f, axArr = plt.subplots(nRows, nCols)
        f.set_size_inches(19.2, 10.8)

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

        for i in range(nCols):
            axArr[0, i].set_title(colNames[i])

        plt.suptitle(r'Correlations - RES NORM')
        f.savefig('%s-corre_res_norm_popx.png' % P.K.califaID)
        plt.close()

#########################################################################
#########################################################################


#########################################################################
############################## PCA LOG ##################################
    npref_f_obs = '%s-logf_obs_' % P.K.califaID
    npref_f_obs_norm = '%s-logf_obs_norm_' % P.K.califaID
    npref_f_syn = '%s-logf_syn_' % P.K.califaID
    npref_f_syn_norm = '%s-logf_syn_norm_' % P.K.califaID

    if args.outputfname:
        npref_f_obs = '%s%s_' % (npref_f_obs, args.outputfname)
        npref_f_obs_norm = '%s%s_' % (npref_f_obs_norm, args.outputfname)
        npref_f_syn = '%s%s_' % (npref_f_syn, args.outputfname)
        npref_f_syn_norm = '%s%s_' % (npref_f_syn_norm, args.outputfname)

    I_obs__zl, ms_obs__l, covMat_obs__ll, eigVal_obs__k, eigVec_obs__lk = P.PCA(np.log10(np.abs(P.f_obs__zl)), P.K.N_zone, 0)
    I_obs_norm__zl, ms_obs_norm__l, covMat_obs_norm__ll, eigVal_obs_norm__k, eigVec_obs_norm__lk = P.PCA(np.log10(np.abs(P.f_obs_norm__zl)), P.K.N_zone, 0)
    I_syn__zl, ms_syn__l, covMat_syn__ll, eigVal_syn__k, eigVec_syn__lk = P.PCA(np.log10(np.abs(P.f_syn__zl)), P.K.N_zone, 0)
    I_syn_norm__zl, ms_syn_norm__l, covMat_syn_norm__ll, eigVal_syn_norm__k, eigVec_syn_norm__lk = P.PCA(np.log10(np.abs(P.f_syn_norm__zl)), P.K.N_zone, 0)

    tomo_obs__zk, tomo_obs__kyx = P.tomogram(I_obs__zl, eigVec_obs__lk)
    tomo_obs_norm__zk, tomo_obs_norm__kyx = P.tomogram(I_obs_norm__zl, eigVec_obs_norm__lk)
    tomo_syn__zk, tomo_syn__kyx = P.tomogram(I_syn__zl, eigVec_syn__lk)
    tomo_syn_norm__zk, tomo_syn_norm__kyx = P.tomogram(I_syn_norm__zl, eigVec_syn_norm__lk)

    if args.tomograms:
        P.screeTestPlot(eigVal_obs__k, args.tmax, npref_f_obs, '%s LOG FOBS' % P.K.califaID)
        P.screeTestPlot(eigVal_obs_norm__k, args.tmax, npref_f_obs_norm, '%s LOG FOBS NORM' % P.K.califaID)
        P.screeTestPlot(eigVal_syn__k, args.tmax, npref_f_syn, '%s LOG FSYN' % P.K.califaID)
        P.screeTestPlot(eigVal_syn_norm__k, args.tmax, npref_f_syn_norm, '%s LOG FSYN NORM' % P.K.califaID)

        for ti in range(args.tmax):
            P.tomoPlot(tomo_obs__kyx, P.l_obs, eigVec_obs__lk, eigVal_obs__k, ms_obs__l, ti, npref_f_obs)
            P.tomoPlot(tomo_obs_norm__kyx, P.l_obs, eigVec_obs_norm__lk, eigVal_obs_norm__k, ms_obs_norm__l, ti, npref_f_obs_norm)
            P.tomoPlot(tomo_syn__kyx, P.l_syn, eigVec_syn__lk, eigVal_syn__k, ms_syn__l, ti, npref_f_syn)
            P.tomoPlot(tomo_syn_norm__kyx, P.l_syn, eigVec_syn_norm__lk, eigVal_syn_norm__k, ms_syn_norm__l, ti, npref_f_syn_norm)


#########################################################################
############################ Correlations ###############################

    if args.correlat:
        colArr = [
                P.K.at_flux__z,
                np.log10(P.K.aZ_flux__z / 0.019),
                P.K.A_V,
                P.K.v_0,
                P.K.v_d,
        ]

        colNames = [
            r'$\log\ t[yr]$',
            r'$\log\ Z\ [Z_\odot]$',
            r'$A_V\ [mag]$',
            r'$v_\star\ [km/s]$',
            r'$\sigma_\star\ [km/s]$',
            r'eigenvector',
        ]
    ############################### OBS NORM ###############################     

        nRows = args.numcorrepc
        nCols = len(colArr) + 1
        f, axArr = plt.subplots(nRows, nCols)
        f.set_size_inches(19.2, 10.8)

        for i in range(nRows):
            axArr[i, 0].set_ylabel('PC%d' % i)

            for j in range(nCols)[:-1]:
                P.correlationAxisPlot(colArr[j], tomo_obs_norm__zk[:, i], axArr[i, j])

            axArr[i, nCols - 1].plot(P.l_obs, eigVec_obs_norm__lk[:, i])
            plt.setp(axArr[i, nCols - 1].get_yticklabels(), visible = False)

        f.subplots_adjust(hspace = 0.0)
        f.subplots_adjust(wspace = 0.05)

        plt.setp([a.get_xticklabels() for a in f.axes[:-nCols]], visible = False)
        plt.setp([a.get_yticklabels() for a in f.axes[::nCols]], visible = True)

        for i in range(nCols):
            axArr[0, i].set_title(colNames[i])

        plt.suptitle(r'Correlations - OBS NORM')
        f.savefig('%s-corre_logf_obs_norm.png' % P.K.califaID)
        plt.close()

    ############################### SYN NORM ###############################     
        f, axArr = plt.subplots(nRows, nCols)
        f.set_size_inches(19.2, 10.8)

        for i in range(nRows):
            axArr[i, 0].set_ylabel('PC%d' % i)

            for j in range(nCols)[:-1]:
                P.correlationAxisPlot(colArr[j], tomo_syn_norm__zk[:, i], axArr[i, j])

            axArr[i, nCols - 1].plot(P.l_syn, eigVec_syn_norm__lk[:, i])
            plt.setp(axArr[i, nCols - 1].get_yticklabels(), visible = False)

        f.subplots_adjust(hspace = 0.0)
        f.subplots_adjust(wspace = 0.05)

        plt.setp([a.get_xticklabels() for a in f.axes[:-nCols]], visible = False)
        plt.setp([a.get_yticklabels() for a in f.axes[::nCols]], visible = True)

        for i in range(nCols):
            axArr[0, i].set_title(colNames[i])

        plt.suptitle(r'Correlations - SYN NORM')
        f.savefig('%s-corre_logf_syn_norm.png' % P.K.califaID)
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

        colNames = [
            r'$10^6\ to\ 10^{7.5}yr$',
            r'$10^{7.5}\ to\ 10^{8.5}yr$',
            r'$10^{8.5}\ to\ 10^{9.5}yr$',
            r'$10^{9.5}\ to\ 10^{10.2}yr$',
            r'eigenvector',
        ]

        nRows = args.numcorrepc
        nCols = len(colArr) + 1
        f, axArr = plt.subplots(nRows, nCols)
        f.set_size_inches(19.2, 10.8)

        for i in range(nRows):
            axArr[i, 0].set_ylabel('PC%d' % i)

            for j in range(nCols)[:-1]:
                P.correlationAxisPlot(colArr[j], tomo_obs_norm__zk[:, i], axArr[i, j])

            axArr[i, nCols - 1].plot(P.l_obs, eigVec_obs_norm__lk[:, i])
            plt.setp(axArr[i, nCols - 1].get_yticklabels(), visible = False)

        f.subplots_adjust(hspace = 0.0)
        f.subplots_adjust(wspace = 0.05)

        plt.setp([a.get_xticklabels() for a in f.axes[:-nCols]], visible = False)
        plt.setp([a.get_yticklabels() for a in f.axes[::nCols]], visible = True)

        for i in range(nCols):
            axArr[0, i].set_title(colNames[i])

        plt.suptitle(r'Correlations - OBS NORM')
        f.savefig('%s-corre_logf_obs_norm_popx.png' % P.K.califaID)
        plt.close()

    ############################### SYN NORM ###############################
        f, axArr = plt.subplots(nRows, nCols)
        f.set_size_inches(19.2, 10.8)

        for i in range(nRows):
            axArr[i, 0].set_ylabel('PC%d' % i)

            for j in range(nCols)[:-1]:
                P.correlationAxisPlot(colArr[j], tomo_syn_norm__zk[:, i], axArr[i, j])

            axArr[i, nCols - 1].plot(P.l_syn, eigVec_syn_norm__lk[:, i])
            plt.setp(axArr[i, nCols - 1].get_yticklabels(), visible = False)

        f.subplots_adjust(hspace = 0.0)
        f.subplots_adjust(wspace = 0.05)

        plt.setp([a.get_xticklabels() for a in f.axes[:-nCols]], visible = False)
        plt.setp([a.get_yticklabels() for a in f.axes[::nCols]], visible = True)

        for i in range(nCols):
            axArr[0, i].set_title(colNames[i])

        plt.suptitle(r'Correlations - SYN NORM')
        f.savefig('%s-corre_logf_syn_norm_popx.png' % P.K.califaID)
        plt.close()
