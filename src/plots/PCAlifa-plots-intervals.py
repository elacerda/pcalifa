'''
Created on 17/07/2013

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
from pystarlight import io

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

def pixDevLineFluxRestFrame(l, res__lz, n_zone, model_line_flux_max = 6562.85, search_delta = 30):
    model_line_flux_max = 6562.85
    search_delta = 10.
    llow = model_line_flux_max - search_delta
    lup = model_line_flux_max + search_delta
    lpixlow = np.where(l >= llow)[0][0]
    lpixup = np.where(l <= lup)[0][-1]
    line_max_lpix__z = lpixlow + res__lz[lpixlow:lpixup, :].argmax(axis = 0)
    l_step = l[1] - l[0]
    dim = np.ones((len(l), n_zone)).T
    lmax__z = ((l * dim).T)[line_max_lpix__z][:, 0]
    dl__z = lmax__z - model_line_flux_max
    dpix__z = np.array((dl__z / l_step))

    return dpix__z, dl__z

def setLineFluxRestFrame(l, f_obs__lz, f_syn__lz, n_zone, restFramePixelDeviance__z):
    interp_step = 0.01
    dpix__z = np.array((1. / interp_step * restFramePixelDeviance__z), dtype = np.int)

    l_pix = np.arange(0, len(l), 1)
    l_interp_pix = np.arange(0, len(l), 0.01)

    interp_f_obs__lz = np.zeros((len(l_interp_pix), n_zone))
    interp_f_syn__lz = np.zeros((len(l_interp_pix), n_zone))
    interp_f_obs_notroll__lz = np.zeros((len(l_interp_pix), n_zone))
    interp_f_syn_notroll__lz = np.zeros((len(l_interp_pix), n_zone))

    f_obs_rf__lz = np.zeros((len(l_pix), n_zone))
    f_syn_rf__lz = np.zeros((len(l_pix), n_zone))

    for z in range(n_zone):
        dev = dpix__z[z]
        interp_f_obs_notroll__lz[:, z] = np.interp(l_interp_pix, l_pix, f_obs__lz[:, z])
        interp_f_syn_notroll__lz[:, z] = np.interp(l_interp_pix, l_pix, f_syn__lz[:, z])

        interp_f_obs__lz[:, z] = np.roll(interp_f_obs_notroll__lz[:, z], dev)
        interp_f_syn__lz[:, z] = np.roll(interp_f_syn_notroll__lz[:, z], dev)

        if dev < 0:
            interp_f_obs__lz[dev:, z] = interp_f_obs_notroll__lz[dev:, z].mean()
            interp_f_syn__lz[dev:, z] = interp_f_syn_notroll__lz[dev:, z].mean()

        f_obs_rf__lz[:, z] = np.interp(l_pix, l_interp_pix, interp_f_obs__lz[:, z])
        f_syn_rf__lz[:, z] = np.interp(l_pix, l_interp_pix, interp_f_syn__lz[:, z])

    return f_obs_rf__lz, f_syn_rf__lz

if __name__ == '__main__':
    args = parser_args()

    P = PCA.PCAlifa(args.califaID, args.fitsDir, args.rFL, args.lc)

    parseArrArgs(args, P.K.N_zone - 1)

    ####################################################################################
    ####################################################################################
    ####################################################################################
    ################ PCA just in intervals defined in StarlightMaskFile ################

    P.setStarlightMaskFile(args.rSEL)

    ###### Setting up the mask just for the lambda intervals in StarlighMaskFile
    mask = (P.maskEmLines == False) & (P.maskLambdaConstrains) & (P.maskQFlag)

    l_obs = P.K.l_obs[mask]
    f_obs__zl = P.K.f_obs[mask].transpose()
    f_syn__zl = P.K.f_syn[mask].transpose()
    # to correct some spurious negative fluxes
    f_obs__zl[f_obs__zl < 0] = f_syn__zl[f_obs__zl < 0]

    f_res__zl = f_obs__zl - f_syn__zl
    f_res_norm__zl = ((P.K.f_obs[mask] - P.K.f_syn[mask]) / P.K.fobs_norm).transpose()

    ### PCA RES ###
    I_res__zl, ms_res__l, covMat_res__ll, eigVal_res__k, eigVec_res__lk = P.PCA(f_res__zl, P.K.N_zone, 0)
    ### PCA RES NORM ###
    I_res_norm__zl, ms_res_norm__l, covMat_res_norm__ll, eigVal_res_norm__k, eigVec_res_norm__lk = P.PCA(f_res_norm__zl, P.K.N_zone, 0)
    ### TOMO RES ### 
    tomo_res__zk, tomo_res__kyx = P.tomogram(I_res__zl, eigVec_res__lk)
    ### TOMO RES NORM ### 
    tomo_res_norm__zk, tomo_res_norm__kyx = P.tomogram(I_res_norm__zl, eigVec_res_norm__lk)

    ### PLOTS ###
    npref_f_res = '%s-f_res_intervals_' % P.K.califaID
    npref_f_res_norm = '%s-f_res_norm_intervals_' % P.K.califaID

    if args.outputfname:
        npref_f_res = '%s%s_' % (npref_f_res, args.outputfname)
        npref_f_res_norm = '%s%s_' % (npref_f_res_norm, args.outputfname)

    #### SANITY CHECK ####
    P.sanityCheck(l_obs, f_res__zl[0, :], npref_f_res)
    ######################

    ####################################################################################
    #################################### TOMOGRAMS #####################################  
    ####################################################################################

    if args.tomograms:
        P.screeTestPlot(eigVal_res__k, args.tmax, npref_f_res, '%s RES' % P.K.califaID)
        P.screeTestPlot(eigVal_res_norm__k, args.tmax, npref_f_res_norm, '%s RES NORM' % P.K.califaID)

        for ti in range(args.tmax):
            P.tomoPlot(tomo_res__kyx, l_obs, eigVec_res__lk, eigVal_res__k, ms_res__l, ti, npref_f_res)
            P.tomoPlot(tomo_res_norm__kyx, l_obs, eigVec_res_norm__lk, eigVal_res_norm__k, ms_res_norm__l, ti, npref_f_res_norm)

    ####################################################################################
    ################################ Rebuild Intervals #################################  
    ####################################################################################

    if args.rebuildSpec:
        if (P.K.N_zone > args.zonesArr[-2]):
            for nZone in args.zonesArr:
                zoneRebuildSpec(nZone, args.eigvArr, l_obs,
                                f_res__zl, tomo_res__zk, eigVec_res__lk, eigVal_res__k, ms_res__l, npref_f_res, True)
                zoneRebuildSpec(nZone, args.eigvArr, l_obs,
                                f_res_norm__zl, tomo_res_norm__zk, eigVec_res_norm__lk, eigVal_res_norm__k, ms_res_norm__l, npref_f_res_norm, True)
        else:
            print "%s N_Zone < %d" % (P.K.califaID, args.zonesArr[-2])

    ####################################################################################
    ################################### Correlations ###################################  
    ####################################################################################

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

        ### RES NORM ###
        f, axArr = plt.subplots(nRows, nCols)
        f.set_size_inches(19.2, 10.8)

        for i in range(nRows):
            axArr[i, 0].set_ylabel('PC%d' % i)

            for j in range(nCols)[:-1]:
                P.correlationAxisPlot(colArr[j], tomo_res_norm__zk[:, i], axArr[i, j])

            axArr[i, nCols - 1].plot(l_obs, eigVec_res_norm__lk[:, i])
            plt.setp(axArr[i, nCols - 1].get_yticklabels(), visible = False)

        f.subplots_adjust(hspace = 0.05)
        f.subplots_adjust(wspace = 0.15)

        plt.setp([a.get_xticklabels() for a in f.axes[:-nCols]], visible = False)
        plt.setp([a.get_yticklabels() for a in f.axes[::nCols]], visible = True)

        for i in range(nCols):
            axArr[0, i].set_title(colNames[i])

        plt.suptitle(r'Correlations - intervals - RES NORM')
        f.savefig('%s-corre_res_norm_intervals.png' % P.K.califaID)
        plt.close()
