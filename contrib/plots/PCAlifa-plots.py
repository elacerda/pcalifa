'''
Created on 10/06/2013

@author: lacerda
'''
import matplotlib
matplotlib.use('agg')
import sys
import numpy as np
import pcalifa as PCA
import argparse as ap
from matplotlib import pyplot as plt
from matplotlib import cm
from scipy import stats as st
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MaxNLocator

plt.rcParams.update({'font.family' : 'serif',
                     'text.usetex' : False,
                     })

def parser_args():
    parser = ap.ArgumentParser(description = 'PCAlifa - correlations')
    parser.add_argument('--fitsfile', '-f',
                        help = 'The file must be named KXXXX*.fits',
                        metavar = 'PyCASSO FITS FILE',
                        type = str,
                        default = None)
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
                        help = 'Plot tomograms.',
                        action = 'store_true',
                        default = False)
    parser.add_argument('--correlat', '-C',
                        help = 'Plot correlations.',
                        action = 'store_true',
                        default = False)
    parser.add_argument('--pcalog',
                        help = 'PCA from logfluxes',
                        action = 'store_true',
                        default = False)
    parser.add_argument('--pcainterval',
                        help = 'PCA from different spectral ranges defined by MASKFILE (mandatory use of --rSEL=MASKFILE)',
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
    parser.add_argument('--outputimgsuffix', '-o',
                        help = 'Suffix of image file. Sometimes denote the image type. (Ex.: image.png)',
                        type = str,
                        default = 'png')

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
    fname = '%szone_%04d.%s' % (nPref, iZone, args.outputimgsuffix)
    f.savefig(fname)
    plt.close()


def correlations_PCxPhys_eigv_plot(P, eigv__lk, PCArr__zk, nPC, l_obs, title, fname):
    prop = {
        'arr'   : [ P.K.at_flux__z, np.log10(P.K.aZ_flux__z / 0.019), P.K.A_V, P.K.v_0, P.K.v_d ],
        'label' : [ r'$\log\ t\ [yr]$', r'$\log\ Z\ [Z_\odot]$', r'$A_V\ [mag]$', r'$v_\star\ [km/s]$', r'$\sigma_\star\ [km/s]$', r'eigenvector' ],
    }

    nCols = len(prop['label'])

    f, axArr = plt.subplots(nPC, nCols)
    f.set_size_inches(19.2, 10.8)

    for i in range(nPC):
        axArr[i, 0].set_ylabel('PC%d' % i)
        ymin = PCArr__zk[:, i].min()
        ymax = PCArr__zk[:, i].max()

        for j in range(nCols)[:-1]:
            ax = axArr[i, j]
            ax.set_ylim(ymin, ymax)
            P.correlationAxisPlot(prop['arr'][j], PCArr__zk[:, i], axArr[i, j])

        axArr[i, nCols - 1].plot(l_obs, eigv__lk[:, i])
        plt.setp(axArr[i, nCols - 1].get_yticklabels(), visible = False)

    f.subplots_adjust(hspace = 0.0)
    f.subplots_adjust(wspace = 0.05)

    plt.setp([a.get_xticklabels() for a in f.axes], rotation = 45)
    plt.setp([a.get_xticklabels() for a in f.axes[:-nCols]], visible = False)
    plt.setp([a.get_yticklabels() for a in f.axes[::nCols]], visible = True)

    for i in range(nCols):
        axArr[0, i].set_title(prop['label'][i])

    plt.suptitle(title)

    if (fname):
        f.savefig(fname)
    else:
        plt.show()

def correlations_PCxPC_plot(P, PCArr__zk, nPC, title, fnamepref):
    prop = {
        'arr'   : [ P.K.at_flux__z, np.log10(P.K.aZ_flux__z / 0.019), P.K.A_V, P.K.v_0, P.K.v_d ],
        'label' : [ r'$\log\ t\ [yr]$', r'$\log\ Z\ [Z_\odot]$', r'$A_V\ [mag]$', r'$v_\star\ [km/s]$', r'$\sigma_\star\ [km/s]$' ],
        'fname' : [ 'at_flux', 'aZ_flux', 'AV', 'v0', 'vd' ]
    }

    for p_i, p in enumerate(prop['arr']):
        nRows = nPC - 1
        nCols = nPC - 1
        f, axArr = plt.subplots(nRows, nCols)

        for ax in f.axes:
            ax.set_axis_off()

        fig_width_pt = 1080.
        inches_per_pt = 1.0 / 72.27
        golden_mean = (5 ** 0.5 - 1.0) / 2.0
        fig_width = fig_width_pt * inches_per_pt
        fig_height = fig_width * golden_mean

        f.set_size_inches(fig_width, fig_height)
        f.set_dpi(200)

        for i in range(nRows):
            y_i = i

            ymin = PCArr__zk[:, y_i].min()
            ymax = PCArr__zk[:, y_i].max()

            for j in range(i, nCols):
                ax = axArr[i, j]
                ax.set_axis_on()

                x_i = j + 1

                xmin = PCArr__zk[:, x_i].min()
                xmax = PCArr__zk[:, x_i].max()
                ax.set_xlim(xmin, xmax)
                ax.set_ylim(ymin, ymax)

                x = PCArr__zk[:, x_i]
                y = PCArr__zk[:, y_i]
                z = p

                plt.setp(ax.get_xticklabels(), visible = False)
                plt.setp(ax.get_yticklabels(), visible = False)

                rhoPearson, pvalPearson = st.pearsonr(x, y)
                rhoSpearman, pvalSpearman = st.spearmanr(x, y)
                pTxt = 'p: %.2f' % rhoPearson
                spTxt = 's: %.2f' % rhoSpearman

                if (prop['fname'])[p_i] == 'vd' :
                    prc = np.percentile(p, 98.)
                    im = ax.scatter(x, y, c = z, edgecolor = 'None', s = 3, cmap = cm.jet, vmax = prc)
                else:
                    im = ax.scatter(x, y, c = z, edgecolor = 'None', s = 3, cmap = cm.jet)

                ax.text(0.95, 0.88, pTxt,
                        fontsize = 10, transform = ax.transAxes,
                        horizontalalignment = 'right',
                        verticalalignment = 'center',
                        multialignment = 'right',
                        weight = 'bold')
                ax.text(0.95, 0.72, spTxt,
                        color = 'red',
                        fontsize = 10, transform = ax.transAxes,
                        horizontalalignment = 'right',
                        verticalalignment = 'center',
                        multialignment = 'right',
                        weight = 'bold')
                plt.setp(ax.get_yticklabels(), visible = False)

            axArr[i, i].set_ylabel('PC%d' % i)
            plt.setp(axArr[i, i].get_xticklabels(), visible = True, rotation = 45)
            plt.setp(axArr[i, i].get_yticklabels(), visible = True, rotation = 45)

        f.subplots_adjust(hspace = 0.0)
        f.subplots_adjust(wspace = 0.0)
        f.subplots_adjust(right = 0.8)
        cbar_ax = f.add_axes([0.85, 0.15, 0.05, 0.7])
        cb = f.colorbar(im, cax = cbar_ax)
        cb.set_label(prop['label'][p_i])

        for i in range(nCols):
            axArr[0, i].set_title('PC%d' % (i + 1))

        plt.suptitle(r'%s - %s' % (title, prop['fname'][p_i]))

        if fnamepref:
            f.savefig('%scorre_PCxPC_%s.%s' % (fnamepref, prop['fname'][p_i], args.outputimgsuffix))
        else:
            plt.show()

def correlations_PCxPopulations_eigv(P, eigv__lk, PCArr__zk, nPC, l_obs, title, fname):
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

    popx = {
        'arr' : [],
        'label': [ r'$10^6\ to\ 10^{7.5}yr$', r'$10^{7.5}\ to\ 10^{8.5}yr$', r'$10^{8.5}\ to\ 10^{9.5}yr$', r'$10^{9.5}\ to\ 10^{10.2}yr$', r'eigenvector', ],
    }

    popx['arr'].append(np.tensordot((P.K.popx.sum(axis = 1))[maskpopx1, :], logt[maskpopx1], (0, 0)) / popxtot)
    popx['arr'].append(np.tensordot((P.K.popx.sum(axis = 1))[maskpopx2, :], logt[maskpopx2], (0, 0)) / popxtot)
    popx['arr'].append(np.tensordot((P.K.popx.sum(axis = 1))[maskpopx3, :], logt[maskpopx3], (0, 0)) / popxtot)
    popx['arr'].append(np.tensordot((P.K.popx.sum(axis = 1))[maskpopx4, :], logt[maskpopx4], (0, 0)) / popxtot)

    nCols = len(popx['label'])
    f, axArr = plt.subplots(nPC, nCols)
    f.set_size_inches(19.2, 10.8)

    for i in range(nPC):
        axArr[i, 0].set_ylabel('PC%d' % i)
        ymin = PCArr__zk[:, i].min()
        ymax = PCArr__zk[:, i].max()

        for j in range(nCols)[:-1]:
            axArr[i, j].set_ylim(ymin, ymax)
            P.correlationAxisPlot(popx['arr'][j], PCArr__zk[:, i], axArr[i, j])

        axArr[i, nCols - 1].plot(l_obs, eigv__lk[:, i])
        plt.setp(axArr[i, nCols - 1].get_yticklabels(), visible = False)

    f.subplots_adjust(hspace = 0.0)
    f.subplots_adjust(wspace = 0.05)

    plt.setp([a.get_xticklabels() for a in f.axes], rotation = 45)
    plt.setp([a.get_xticklabels() for a in f.axes[:-nCols]], visible = False)
    plt.setp([a.get_yticklabels() for a in f.axes[::nCols]], visible = True)

    for i in range(nCols):
        axArr[0, i].set_title(popx['label'][i])

    plt.suptitle(title)

    if fname:
        f.savefig(fname)
    else:
        plt.show()

if __name__ == '__main__':
    args = parser_args()

    P = PCA.PCAlifa(args.fitsfile, args.rFL, args.lc)

    parseArrArgs(args, P.K.N_zone - 1)

    if args.rSEL:
        P.setStarlightMaskFile(args.rSEL)

    P.setImgSuffix(args.outputimgsuffix)

    P.runPCA()

    if not args.pcainterval:
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
            P.screeTestPlot(P.eigVal_obs__k, args.tmax, '%s FOBS' % P.K.califaID, npref_f_obs)
            P.screeTestPlot(P.eigVal_obs_norm__k, args.tmax, '%s FOBS NORM' % P.K.califaID, npref_f_obs_norm)
            P.screeTestPlot(P.eigVal_syn__k, args.tmax, '%s FSYN' % P.K.califaID, npref_f_syn)
            P.screeTestPlot(P.eigVal_syn_norm__k, args.tmax, '%s FSYN NORM' % P.K.califaID, npref_f_syn_norm)
            P.screeTestPlot(P.eigVal_res__k, args.tmax, '%s RES' % P.K.califaID, npref_f_res)
            P.screeTestPlot(P.eigVal_res_norm__k, args.tmax, '%s RES NORM' % P.K.califaID, npref_f_res_norm)

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
            nRows = args.numcorrepc

            #########################################################################
            ######################## PC-Phys Correlations ###########################
            ############################### OBS NORM ################################
            correlations_PCxPhys_eigv_plot(P, P.eigVec_obs_norm__lk, P.tomo_obs_norm__zk, nRows, P.l_obs, r'Correlations - OBS NORM', '%scorre.%s' % (npref_f_obs_norm, args.outputimgsuffix))
            ############################### SYN NORM ################################
            correlations_PCxPhys_eigv_plot(P, P.eigVec_syn_norm__lk, P.tomo_syn_norm__zk, nRows, P.l_syn, r'Correlations - SYN NORM', '%scorre.%s' % (npref_f_syn_norm, args.outputimgsuffix))
            ############################### RES NORM ################################
            correlations_PCxPhys_eigv_plot(P, P.eigVec_res_norm__lk, P.tomo_res_norm__zk, nRows, P.l_obs, r'Correlations - RES NORM', '%scorre.%s' % (npref_f_res_norm, args.outputimgsuffix))

            #########################################################################
            ######################### PC-PC Correlations ############################
            ############################### OBS NORM ################################
            correlations_PCxPC_plot(P, P.tomo_obs_norm__zk, nRows, 'Correlations - OBS NORM', npref_f_obs_norm)
            ############################### SYN NORM ################################
            correlations_PCxPC_plot(P, P.tomo_syn_norm__zk, nRows, 'Correlations - SYN NORM', npref_f_syn_norm)
            ############################### RES NORM ################################
            correlations_PCxPC_plot(P, P.tomo_res_norm__zk, nRows, 'Correlations - RES NORM', npref_f_res_norm)

            #########################################################################
            ###################### Population Correlations ##########################
            ############################### OBS NORM ################################
            correlations_PCxPopulations_eigv(P, P.eigVec_obs_norm__lk, P.tomo_obs_norm__zk, nRows, P.l_obs, r'Correlations - OBS NORM', '%scorre_popx.%s' % (npref_f_obs_norm, args.outputimgsuffix))
            ############################### SYN NORM ###############################
            correlations_PCxPopulations_eigv(P, P.eigVec_syn_norm__lk, P.tomo_syn_norm__zk, nRows, P.l_syn, r'Correlations - SYN NORM', '%scorre_popx.%s' % (npref_f_syn_norm, args.outputimgsuffix))
            ############################### RES NORM ###############################
            correlations_PCxPopulations_eigv(P, P.eigVec_res_norm__lk, P.tomo_res_norm__zk, nRows, P.l_obs, r'Correlations - RES NORM', '%scorre_popx.%s' % (npref_f_res_norm, args.outputimgsuffix))

        #########################################################################
        #########################################################################

        #########################################################################
        ############################### PCA LOG #################################
        if args.pcalog:
            npref_f_obs = '%s-logf_obs_' % P.K.califaID
            npref_f_obs_norm = '%s-logf_obs_norm_' % P.K.califaID
            npref_f_syn = '%s-logf_syn_' % P.K.califaID
            npref_f_syn_norm = '%s-logf_syn_norm_' % P.K.califaID

            if args.outputfname:
                npref_f_obs = '%s%s_' % (npref_f_obs, args.outputfname)
                npref_f_obs_norm = '%s%s_' % (npref_f_obs_norm, args.outputfname)
                npref_f_syn = '%s%s_' % (npref_f_syn, args.outputfname)
                npref_f_syn_norm = '%s%s_' % (npref_f_syn_norm, args.outputfname)

            #### SANITY CHECK ####
            P.sanityCheck(P.l_obs, np.log10(np.abs(P.f_obs__zl[0, :])), npref_f_obs)
            P.sanityCheck(P.l_obs, np.log10(np.abs(P.f_obs_norm__zl[0, :])), npref_f_obs_norm)
            P.sanityCheck(P.l_syn, np.log10(np.abs(P.f_syn__zl[0, :])), npref_f_syn)
            P.sanityCheck(P.l_syn, np.log10(np.abs(P.f_syn_norm__zl[0, :])), npref_f_syn_norm)
            ######################

            I_obs__zl, ms_obs__l, covMat_obs__ll, eigVal_obs__k, eigVec_obs__lk = P.PCA(np.log10(np.abs(P.f_obs__zl)), P.K.N_zone, 0)
            I_obs_norm__zl, ms_obs_norm__l, covMat_obs_norm__ll, eigVal_obs_norm__k, eigVec_obs_norm__lk = P.PCA(np.log10(np.abs(P.f_obs_norm__zl)), P.K.N_zone, 0)
            I_syn__zl, ms_syn__l, covMat_syn__ll, eigVal_syn__k, eigVec_syn__lk = P.PCA(np.log10(np.abs(P.f_syn__zl)), P.K.N_zone, 0)
            I_syn_norm__zl, ms_syn_norm__l, covMat_syn_norm__ll, eigVal_syn_norm__k, eigVec_syn_norm__lk = P.PCA(np.log10(np.abs(P.f_syn_norm__zl)), P.K.N_zone, 0)

            tomo_obs__zk, tomo_obs__kyx = P.tomogram(I_obs__zl, eigVec_obs__lk)
            tomo_obs_norm__zk, tomo_obs_norm__kyx = P.tomogram(I_obs_norm__zl, eigVec_obs_norm__lk)
            tomo_syn__zk, tomo_syn__kyx = P.tomogram(I_syn__zl, eigVec_syn__lk)
            tomo_syn_norm__zk, tomo_syn_norm__kyx = P.tomogram(I_syn_norm__zl, eigVec_syn_norm__lk)

            if args.tomograms:
                P.screeTestPlot(eigVal_obs__k, args.tmax, '%s LOG FOBS' % P.K.califaID, npref_f_obs)
                P.screeTestPlot(eigVal_obs_norm__k, args.tmax, '%s LOG FOBS NORM' % P.K.califaID, npref_f_obs_norm)
                P.screeTestPlot(eigVal_syn__k, args.tmax, '%s LOG FSYN' % P.K.califaID, npref_f_syn)
                P.screeTestPlot(eigVal_syn_norm__k, args.tmax, '%s LOG FSYN NORM' % P.K.califaID, npref_f_syn_norm)

                for ti in range(args.tmax):
                    P.tomoPlot(tomo_obs__kyx, P.l_obs, eigVec_obs__lk, eigVal_obs__k, ms_obs__l, ti, npref_f_obs)
                    P.tomoPlot(tomo_obs_norm__kyx, P.l_obs, eigVec_obs_norm__lk, eigVal_obs_norm__k, ms_obs_norm__l, ti, npref_f_obs_norm)
                    P.tomoPlot(tomo_syn__kyx, P.l_syn, eigVec_syn__lk, eigVal_syn__k, ms_syn__l, ti, npref_f_syn)
                    P.tomoPlot(tomo_syn_norm__kyx, P.l_syn, eigVec_syn_norm__lk, eigVal_syn_norm__k, ms_syn_norm__l, ti, npref_f_syn_norm)

            #########################################################################
            ############################# Correlations ##############################
            if args.correlat:
                nRows = args.numcorrepc
                #########################################################################
                ######################## PC-Phys Correlations ###########################
                ############################### OBS NORM ###############################
                correlations_PCxPhys_eigv_plot(P, eigVec_obs_norm__lk, tomo_obs_norm__zk, nRows, P.l_obs, r'Correlations - OBS NORM', '%scorre.%s' % (npref_f_obs_norm, args.outputimgsuffix))
                ############################### SYN NORM ###############################
                correlations_PCxPhys_eigv_plot(P, eigVec_syn_norm__lk, tomo_syn_norm__zk, nRows, P.l_syn, r'Correlations - SYN NORM', '%scorre.%s' % (npref_f_syn_norm, args.outputimgsuffix))

                #########################################################################
                ######################### PC-PC Correlations ############################
                ############################### OBS NORM ################################
                correlations_PCxPC_plot(P, tomo_obs_norm__zk, nRows, 'Correlations - OBS NORM', npref_f_obs_norm)
                ############################### SYN NORM ################################
                correlations_PCxPC_plot(P, tomo_syn_norm__zk, nRows, 'Correlations - SYN NORM', npref_f_syn_norm)

                #########################################################################
                ###################### Population Correlations ##########################
                ############################### OBS NORM ################################
                correlations_PCxPopulations_eigv(P, eigVec_obs_norm__lk, tomo_obs_norm__zk, nRows, P.l_obs, r'Correlations - OBS NORM', '%scorre_popx.%s' % (npref_f_obs_norm, args.outputimgsuffix))
                ############################### SYN NORM ###############################
                correlations_PCxPopulations_eigv(P, eigVec_syn_norm__lk, tomo_syn_norm__zk, nRows, P.l_syn, r'Correlations - SYN NORM', '%scorre_popx.%s' % (npref_f_syn_norm, args.outputimgsuffix))

    #########################################################################
    ######### PCA just in intervals defined in StarlightMaskFile ############
    else:
        ###### Setting up the mask just for the lambda intervals in StarlighMaskFile
        mask = (P.maskEmLines == False) & (P.maskLambdaConstrains) & (P.maskQFlag)

        l_obs = P.K.l_obs[mask]
        f_obs__zl = P.K.f_obs[mask].transpose()
        f_syn__zl = P.K.f_syn[mask].transpose()

        # in order to correct some spurious negative fluxes
        f_obs__zl[f_obs__zl < 0] = f_syn__zl[f_obs__zl < 0]

        f_res__zl = f_obs__zl - f_syn__zl
        f_res_norm__zl = ((P.K.f_obs[mask] - P.K.f_syn[mask]) / P.K.fobs_norm).transpose()

        I_res__zl, ms_res__l, covMat_res__ll, eigVal_res__k, eigVec_res__lk = P.PCA(f_res__zl, P.K.N_zone, 0)
        I_res_norm__zl, ms_res_norm__l, covMat_res_norm__ll, eigVal_res_norm__k, eigVec_res_norm__lk = P.PCA(f_res_norm__zl, P.K.N_zone, 0)
        tomo_res__zk, tomo_res__kyx = P.tomogram(I_res__zl, eigVec_res__lk)
        tomo_res_norm__zk, tomo_res_norm__kyx = P.tomogram(I_res_norm__zl, eigVec_res_norm__lk)

        npref_f_res = '%s-f_res_intervals_' % P.K.califaID
        npref_f_res_norm = '%s-f_res_norm_intervals_' % P.K.califaID

        if args.outputfname:
            npref_f_res = '%s%s_' % (npref_f_res, args.outputfname)
            npref_f_res_norm = '%s%s_' % (npref_f_res_norm, args.outputfname)

        #### SANITY CHECK ####
        P.sanityCheck(l_obs, f_res__zl[0, :], npref_f_res)
        ######################

        #########################################################################
        ############################## Tomograms ################################
        if args.tomograms:
            P.screeTestPlot(eigVal_res__k, args.tmax, '%s RES' % P.K.califaID, npref_f_res)
            P.screeTestPlot(eigVal_res_norm__k, args.tmax, '%s RES NORM' % P.K.califaID, npref_f_res_norm)

            for ti in range(args.tmax):
                P.tomoPlot(tomo_res__kyx, l_obs, eigVec_res__lk, eigVal_res__k, ms_res__l, ti, npref_f_res)
                P.tomoPlot(tomo_res_norm__kyx, l_obs, eigVec_res_norm__lk, eigVal_res_norm__k, ms_res_norm__l, ti, npref_f_res_norm)

        #########################################################################
        ########################## Rebuild Intervals ############################
        if args.rebuildSpec:
            if (P.K.N_zone > args.zonesArr[-2]):
                for nZone in args.zonesArr:
                    zoneRebuildSpec(nZone, args.eigvArr, l_obs,
                                    f_res__zl, tomo_res__zk, eigVec_res__lk, eigVal_res__k, ms_res__l, npref_f_res, True)
                    zoneRebuildSpec(nZone, args.eigvArr, l_obs,
                                    f_res_norm__zl, tomo_res_norm__zk, eigVec_res_norm__lk, eigVal_res_norm__k, ms_res_norm__l, npref_f_res_norm, True)
            else:
                print "%s N_Zone < %d" % (P.K.califaID, args.zonesArr[-2])
        #########################################################################
        ############################ Correlations ###############################
        if args.correlat:
            nRows = args.numcorrepc

            ############################### RES NORM ###############################
            correlations_PCxPhys_eigv_plot(P, eigVec_res_norm__lk, tomo_res_norm__zk, nRows, l_obs, r'Correlations - intervals - RES NORM', '%scorre.%s' % (npref_f_res_norm, args.outputimgsuffix))

            ######################### PC-PC Correlations ############################
            ############################### RES NORM ################################
            correlations_PCxPC_plot(P, tomo_res_norm__zk, nRows, 'Correlations - intervals - RES NORM', npref_f_res_norm)

            ###################### Population Correlations ##########################
            ############################### RES NORM ###############################
            correlations_PCxPopulations_eigv(P, eigVec_res_norm__lk, tomo_res_norm__zk, nRows, l_obs, r'Correlations - intervals- RES NORM', '%scorre_popx.%s' % (npref_f_res_norm, args.outputimgsuffix))
