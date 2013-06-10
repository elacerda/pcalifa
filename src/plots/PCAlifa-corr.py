'''
Created on 02/05/2013

@author: lacerda
'''
import matplotlib
matplotlib.use('agg')

import sys
import numpy as np
import PCAlifa as PCA
import argparse as ap
from matplotlib import pyplot as plt

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

    return parser.parse_args()

if __name__ == '__main__':
    args = parser_args()

    P = PCA.PCAlifa(args.califaID, args.fitsDir, args.rFL, args.lc)

    if args.rSEL:
        P.setStarlightMaskFile(args.rSEL)

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

############################### SYN NORM ###############################     
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

    ############################### SYN NORM ###############################
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
