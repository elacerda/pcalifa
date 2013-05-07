'''
Created on 02/05/2013

@author: lacerda
'''
import matplotlib
matplotlib.use('agg')

import sys
import numpy as np
import scipy.stats as st
from matplotlib import pyplot as plt
import PCAlifa as PCA

#fitsDir = '/home/lacerda/CALIFA'
#califaID = 'K0802'
fitsDir = sys.argv[1]
califaID = sys.argv[2]

maskfile = '/home/lacerda/workspace/PCA/src/Mask.mC'
flagLinesQuantil = 0.9
remFlaggedLambdas = True
remStarlightEmLines = True
tmax = 20 # numero maximo de eigenvalues

def corrPlot(x, y, ax):
        rhoPearson, pvalPearson = st.pearsonr(x, y)
        rhoSpearman, pvalSpearman = st.spearmanr(x, y)
        pTxt = 'p: %.2f' % rhoPearson
        spTxt = 's: %.2f' % rhoSpearman
        ax.plot(x, y, '.')
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

if __name__ == '__main__':
    P = PCA.PCAlifa(califaID = califaID,
                    fitsDir = fitsDir,
                    flagLinesQuantil = flagLinesQuantil,
                    remFlaggedLambdas = remFlaggedLambdas)

    if (remStarlightEmLines == True):
        P.removeStarlightEmLines(maskfile)

    P.PCA_obs_norm()
    P.tomograms_obs_norm()
    P.PCA_syn_norm()
    P.tomograms_syn_norm()
    P.PCA_res_norm()
    P.tomograms_res_norm()

    colArr = [
            P.K.at_flux__z,
            P.K.aZ_flux__z / 0.019,
            P.K.A_V,
            P.K.v_0,
            P.K.v_d
    ]

    nRows = 10
    nCols = len(colArr) + 1
    f, axArr = plt.subplots(nRows, nCols)
    f.set_size_inches(19.8, 10.8)

    for i in range(nRows):
        axArr[i, 0].set_ylabel('PC%d' % i)

        for j in range(nCols)[:-1]:
            corrPlot(colArr[j], P.tomo_obs_norm__zk[:, i], axArr[i, j])

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
    f.savefig('%s-corre_obs_norm_0-9.png' % califaID)
    plt.close()

############################### SYN NORM ###############################     
    f, axArr = plt.subplots(nRows, nCols)
    f.set_size_inches(19.8, 10.8)

    for i in range(nRows):
        axArr[i, 0].set_ylabel('PC%d' % i)

        for j in range(nCols)[:-1]:
            corrPlot(colArr[j], P.tomo_syn_norm__zk[:, i], axArr[i, j])

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
    f.savefig('%s-corre_syn_norm_0-9.png' % califaID)
    plt.close()

############################### SYN NORM ###############################     
    f, axArr = plt.subplots(nRows, nCols)
    f.set_size_inches(19.8, 10.8)

    for i in range(nRows):
        axArr[i, 0].set_ylabel('PC%d' % i)

        for j in range(nCols)[:-1]:
            corrPlot(colArr[j], P.tomo_res_norm__zk[:, i], axArr[i, j])

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
    f.savefig('%s-corre_res_norm_0-9.png' % califaID)
    plt.close()

############################### POP CORR ###############################
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
            corrPlot(colArr[j], P.tomo_obs_norm__zk[:, i], axArr[i, j])

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
    f.savefig('%s-corre_obs_norm_popx_0-9.png' % califaID)
    plt.close()

    ############################### SYN NORM ###############################
    f, axArr = plt.subplots(nRows, nCols)
    f.set_size_inches(19.8, 10.8)

    for i in range(nRows):
        axArr[i, 0].set_ylabel('PC%d' % i)

        for j in range(nCols)[:-1]:
            corrPlot(colArr[j], P.tomo_syn_norm__zk[:, i], axArr[i, j])

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
    f.savefig('%s-corre_syn_norm_popx_0-9.png' % califaID)
    plt.close()

    ############################### SYN NORM ###############################
    f, axArr = plt.subplots(nRows, nCols)
    f.set_size_inches(19.8, 10.8)

    for i in range(nRows):
        axArr[i, 0].set_ylabel('PC%d' % i)

        for j in range(nCols)[:-1]:
            corrPlot(colArr[j], P.tomo_res_norm__zk[:, i], axArr[i, j])

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
    f.savefig('%s-corre_res_norm_popx_0-9.png' % califaID)
    plt.close()
