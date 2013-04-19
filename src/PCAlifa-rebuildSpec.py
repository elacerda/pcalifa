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
maskfile = '/home/lacerda/workspace/PCA/src/Mask.mC'
flagLinesQuantil = 0.9
remFlaggedLambdas = True
remStarlightEmLines = False
tmax = 20 # numero maximo de eigenvalues

def zoneRebuildPlot(P, nZone, evRebArr, O, tomo, evec, eval, mean, nPref):
    i = 0
    nCols = 2
    nRows = len(evRebArr) / nCols
    fig = plt.figure(figsize = (19.8, 10.8))
    #plt.rcParams.update({'font.family': 'serif', 'text.usetex': True, 'backend': 'ps'})
    gs = gridspec.GridSpec(nRows, nCols, width_ratios = [1, 1], hspace = 0)

    for ne in evRebArr:
        ax = plt.subplot(gs[i])

        I_reb__zl, M = P.rebuildSpectra(tomo, evec, mean, ne)
        diff = O[zone, :] - M[zone, :]

        adev = 100. * (1. / P.K.N_zone) * (np.abs(diff) / O[zone, :]).sum()
        sigmaNReb = 0.
        sigmaReb = np.sqrt(eval[:ne].sum())

        if (ne + 1 != len(evec)):
            sigmaNReb = np.sqrt(eval[ne + 1:].sum())

        sigmaRatio = sigmaNReb / sigmaReb

        ax.plot(P.l_obs, O[zone, :], label = 'Obs')
        ax.plot(P.l_obs, M[zone, :], label = 'Mod')
        ax.plot(P.l_obs, diff, label = 'Res')
        ax.xaxis.set_major_locator(MaxNLocator(20))
        textStr = 'adev =  %.4f %% - sigmaReb = %.2e - sigmaNReb = %.2e - ratio = %.4f' % (adev, sigmaReb, sigmaNReb, sigmaRatio)
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

if __name__ == '__main__':
    P = PCA.PCAlifa(califaID = califaID,
                    fitsDir = fitsDir,
                    flagLinesQuantil = flagLinesQuantil,
                    remFlaggedLambdas = remFlaggedLambdas,
                    runDefaultPCA = True)

    if (remStarlightEmLines == True):
        P.removeStarlightEmLines(maskfile)

    if (P.K.N_zone > 200):
        zonesRebuild = np.array([0, 10, 20, 100, 200, P.K.N_zone - 1])

        for zone in zonesRebuild:
            eigvecRebuildArr = np.array([1, 2, 3, 4, 5, 6, 10, 20])
            nPref = '%s-f_obs' % P.K.califaID
            zoneRebuildPlot(P, zone, eigvecRebuildArr, P.f_obs__zl,
                            P.tomo_obs__zk, P.eigVec_obs__lk,
                            P.eigVal_obs__k, P.ms_obs, nPref)

            nPref = '%s-f_obs_norm' % P.K.califaID
            zoneRebuildPlot(P, zone, eigvecRebuildArr, P.f_obs_norm__zl,
                            P.tomo_obs_norm__zk, P.eigVec_obs_norm__lk,
                            P.eigVal_obs_norm__k, P.ms_obs_norm, nPref)

            nPref = '%s-f_syn' % P.K.califaID
            zoneRebuildPlot(P, zone, eigvecRebuildArr, P.f_syn__zl,
                            P.tomo_syn__zk, P.eigVec_syn__lk,
                            P.eigVal_syn__k, P.ms_syn, nPref)

            nPref = '%s-f_syn_norm' % P.K.califaID
            zoneRebuildPlot(P, zone, eigvecRebuildArr, P.f_syn_norm__zl,
                            P.tomo_syn_norm__zk, P.eigVec_syn_norm__lk,
                            P.eigVal_syn_norm__k, P.ms_syn_norm, nPref)

            nPref = '%s-f_res' % P.K.califaID
            zoneRebuildPlot(P, zone, eigvecRebuildArr, P.f_res__zl,
                            P.tomo_res__zk, P.eigVec_res__lk,
                            P.eigVal_res__k, P.ms_res, nPref)

            nPref = '%s-f_res' % P.K.califaID
            zoneRebuildPlot(P, zone, eigvecRebuildArr, P.f_res_norm__zl,
                            P.tomo_res_norm__zk, P.eigVec_res_norm__lk,
                            P.eigVal_res_norm__k, P.ms_res_norm, nPref)
