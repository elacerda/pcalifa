'''
Created on 18/04/2013

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
#fitsDir = '/home/lacerda/CALIFA'
#califaID = 'K0277'

maskfile = '/home/lacerda/workspace/PCA/src/Mask.mC'
flagLinesQuantil = 0.9
remFlaggedLambdas = True
remStarlightEmLines = False
tmax = 20 # numero maximo de eigenvalues

def zoneRebuildPlot(P, nZone, evRebArr, O, tomo, eVec, eVal, mean, nPref, resid = False):
    i = 0
    nCols = 2
    nRows = len(evRebArr) / nCols
    fig = plt.figure(figsize = (19.8, 10.8))
    #plt.rcParams.update({'font.family': 'serif', 'text.usetex': True, 'backend': 'ps'})
    gs = gridspec.GridSpec(nRows, nCols, width_ratios = [1, 1], hspace = 0)

    adev = np.zeros_like(evRebArr)

    for i, ne in enumerate(evRebArr):
        ax = plt.subplot(gs[i])

        I_reb__zl, M = P.rebuildSpectra(tomo, eVec, mean, ne)
        diff = O[zone, :] - M[zone, :]

        sigmaNReb = 0.
        sigmaReb = np.sqrt(eVal[:ne].sum())

        if (ne + 1 != len(eVec)):
            sigmaNReb = np.sqrt(eVal[ne + 1:].sum())

        sigmaRatio = sigmaNReb / sigmaReb

        textStrAdev = ''

        if (resid == True):
            ax.set_ylim([1.1 * O[zone, :].min(), 1.1 * O[zone, :].max()])
        else:
            ax.set_ylim([-0.5 * O[zone, :].mean(), 1.5 * O[zone, :].mean()])
            adev[i] = 100. * (1. / P.K.N_zone) * (np.abs(diff) / O[zone, :]).sum()
            textStrAdev = 'adev =  %.4f %% - ' % adev[i]

        textStr = '%ssigmaReb = %.2e - sigmaNReb = %.2e - ratio = %.4f' % (textStrAdev, sigmaReb, sigmaNReb, sigmaRatio)

        ax.plot(P.l_obs, O[zone, :], label = 'Obs',
                )
#                linewidth = 0.8)
        ax.plot(P.l_obs, M[zone, :], label = 'Mod',
                )
#                linewidth = 0.8)
        ax.plot(P.l_obs, diff * 5., label = 'Res x 5')
        ax.xaxis.set_major_locator(MaxNLocator(20))
        ax.text(0.01, 0.92, textStr,
                fontsize = 10,
                transform = ax.transAxes,
                horizontalalignment = 'left',
                verticalalignment = 'center',
                multialignment = 'left')
        ax.set_ylabel('using %d eigvec' % ne)
        ax.legend(prop = {'size':7})
        ax.grid()

        if (i < len(evRebArr) - 2):
            plt.setp(ax.get_xticklabels(), visible = False)
        elif (i == len(evRebArr) - 2):
            i = -1

        i = i + 2

    plt.suptitle('CALIFA ID: %s%sZone %04d' % (P.K.califaID, ' ' * 60, zone))
    plt.tight_layout(pad = 2., w_pad = 0., h_pad = 0.)
    fig.savefig('%s-zone-%04d.png' % (nPref, zone))
    plt.close()

    if (resid == False):
        fig = plt.figure(figsize = (19.8, 10.8))
        plt.plot(evRebArr, adev, label = 'adev')
        plt.legend()
        plt.xlabel('Number of eigenvectors')
        plt.ylabel('adev')
        fig.savefig('%s-zone-%04d-adev.png' % (nPref, zone))
        plt.close()

if __name__ == '__main__':
    P = PCA.PCAlifa(califaID = califaID,
                    fitsDir = fitsDir,
                    flagLinesQuantil = flagLinesQuantil,
                    remFlaggedLambdas = remFlaggedLambdas,
                    runDefaultPCA = False)

    if (remStarlightEmLines == True):
        P.removeStarlightEmLines(maskfile)

    if (P.K.N_zone > 200):
        zonesRebuild = np.array([0, 10, 20, 100, 200, P.K.N_zone - 1])

        for zone in zonesRebuild:
            eigvecRebuildArr = np.array([1, 2, 3, 4, 5, 6, 10, 20])

            P.PCA_obs()
            P.tomograms_obs()

            nPref = '%s-f_obs' % P.K.califaID
            zoneRebuildPlot(P, zone, eigvecRebuildArr, P.f_obs__zl,
                            P.tomo_obs__zk, P.eigVec_obs__lk,
                            P.eigVal_obs__k, P.ms_obs__l, nPref)

            P.PCA_obs_norm()
            P.tomograms_obs_norm()

            nPref = '%s-f_obs_norm' % P.K.califaID
            zoneRebuildPlot(P, zone, eigvecRebuildArr, P.f_obs_norm__zl,
                            P.tomo_obs_norm__zk, P.eigVec_obs_norm__lk,
                            P.eigVal_obs_norm__k, P.ms_obs_norm__l, nPref)

            P.PCA_syn()
            P.tomograms_syn()

            nPref = '%s-f_syn' % P.K.califaID
            zoneRebuildPlot(P, zone, eigvecRebuildArr, P.f_syn__zl,
                            P.tomo_syn__zk, P.eigVec_syn__lk,
                            P.eigVal_syn__k, P.ms_syn__l, nPref)

            P.PCA_syn_norm()
            P.tomograms_syn_norm()

            nPref = '%s-f_syn_norm' % P.K.califaID
            zoneRebuildPlot(P, zone, eigvecRebuildArr, P.f_syn_norm__zl,
                            P.tomo_syn_norm__zk, P.eigVec_syn_norm__lk,
                            P.eigVal_syn_norm__k, P.ms_syn_norm__l, nPref)

            P.PCA_res()
            P.tomograms_res()

            nPref = '%s-f_res' % P.K.califaID
            zoneRebuildPlot(P, zone, eigvecRebuildArr, P.f_res__zl,
                            P.tomo_res__zk, P.eigVec_res__lk,
                            P.eigVal_res__k, P.ms_res__l, nPref, resid = True)

            P.PCA_res_norm()
            P.tomograms_res_norm()

            nPref = '%s-f_res_norm' % P.K.califaID
            zoneRebuildPlot(P, zone, eigvecRebuildArr, P.f_res_norm__zl,
                            P.tomo_res_norm__zk, P.eigVec_res_norm__lk,
                            P.eigVal_res_norm__k, P.ms_res_norm__l, nPref, resid = True)
    else:
        print "%s N_Zone < 200" % P.K.califaID
