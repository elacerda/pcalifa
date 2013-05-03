'''
Created on 02/05/2013

@author: lacerda
'''
import matplotlib
matplotlib.use('agg')

import sys
import numpy as np
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

def corr(a, b):
    cov = (a - a.mean()) * (b - b.mean())
    adevquad = (a - a.mean()) ** 2.0
    bdevquad = (b - b.mean()) ** 2.0

    c = cov.sum() / np.sqrt(adevquad.sum() * bdevquad.sum())

    return c

def corrPlot(x, y, ax):
        c = corr(x, y)
        label = '%.2f' % c
        ax.plot(x, y, '.', label = label)
        ax.text(0.87, 0.88, label,
                fontsize = 10, transform = ax.transAxes,
                horizontalalignment = 'left',
                verticalalignment = 'center',
                multialignment = 'left',
                weight = 'bold')
        plt.setp(ax.get_yticklabels(), visible = False)

if __name__ == '__main__':
    P = PCA.PCAlifa(califaID = califaID,
                    fitsDir = fitsDir,
                    flagLinesQuantil = flagLinesQuantil,
                    remFlaggedLambdas = remFlaggedLambdas,
                    runDefaultPCA = False)

    if (remStarlightEmLines == True):
        P.removeStarlightEmLines(maskfile)

    P.PCA_obs_norm()
    P.tomograms_obs_norm()

    colArr = [
            P.K.at_flux__z,
            P.K.aZ_flux__z / 0.019,
            P.K.A_V,
            P.K.v_0,
            P.K.v_d
    ]

    nRows = 10
    nCols = len(colArr)
    f, axArr = plt.subplots(nRows, nCols)
    f.set_size_inches(19.8, 10.8)

    for i in range(nRows):
        axArr[i, 0].set_ylabel('PC%d' % i)
        for j in range(nCols):
            corrPlot(colArr[j], P.tomo_obs_norm__zk[:, i], axArr[i, j])

    f.subplots_adjust(hspace = 0.0)
    f.subplots_adjust(wspace = 0.0)

    plt.setp([a.get_xticklabels() for a in f.axes[:-5]], visible = False)
    plt.setp([a.get_yticklabels() for a in f.axes[::5]], visible = True)
    axArr[0, 0].set_title('at_flux')
    axArr[0, 1].set_title('aZ_flux')
    axArr[0, 2].set_title('A_V')
    axArr[0, 3].set_title('v_0')
    axArr[0, 4].set_title('v_d')

    #plt.legend()
    plt.suptitle('Correlations PC0 ... PC9 - OBS NORM')
    #plt.tight_layout(pad = 2., w_pad = 0.5, h_pad = 0.5)
    f.savefig('%s-corre_obs_norm_0-9.png' % califaID)
    plt.close()

    f, axArr = plt.subplots(nRows, nCols)
    f.set_size_inches(19.8, 10.8)

    for i in range(nRows):
        k = i + nRows
        axArr[i, 0].set_ylabel('PC%d' % k)
        for j in range(nCols):
            corrPlot(colArr[j], P.tomo_obs_norm__zk[:, k], axArr[i, j])

    f.subplots_adjust(hspace = 0.0)
    f.subplots_adjust(wspace = 0.0)

    plt.setp([a.get_xticklabels() for a in f.axes[:-5]], visible = False)
    plt.setp([a.get_yticklabels() for a in f.axes[::5]], visible = True)
    axArr[0, 0].set_title('at_flux')
    axArr[0, 1].set_title('aZ_flux')
    axArr[0, 2].set_title('A_V')
    axArr[0, 3].set_title('v_0')
    axArr[0, 4].set_title('v_d')

    #plt.legend()
    plt.suptitle('Correlations PC10 ... PC19 - OBS NORM')
    #plt.tight_layout(pad = 2., w_pad = 0.5, h_pad = 0.5)
    f.savefig('%s-corre_obs_norm_10-19.png' % califaID)
    plt.close()


############################### SYN NORM ###############################     
    P.PCA_syn_norm()
    P.tomograms_syn_norm()

    f, axArr = plt.subplots(nRows, nCols)
    f.set_size_inches(19.8, 10.8)

    for i in range(nRows):
        axArr[i, 0].set_ylabel('PC%d' % i)
        for j in range(nCols):
            corrPlot(colArr[j], P.tomo_syn_norm__zk[:, i], axArr[i, j])

    f.subplots_adjust(hspace = 0.0)
    f.subplots_adjust(wspace = 0.0)

    plt.setp([a.get_xticklabels() for a in f.axes[:-5]], visible = False)
    plt.setp([a.get_yticklabels() for a in f.axes[::5]], visible = True)
    axArr[0, 0].set_title('at_flux')
    axArr[0, 1].set_title('aZ_flux')
    axArr[0, 2].set_title('A_V')
    axArr[0, 3].set_title('v_0')
    axArr[0, 4].set_title('v_d')

    #plt.legend()
    plt.suptitle('Correlations PC0 ... PC9 - SYN NORM')
    #plt.tight_layout(pad = 2., w_pad = 0.5, h_pad = 0.5)
    f.savefig('%s-corre_syn_norm_0-9.png' % califaID)
    plt.close()

    f, axArr = plt.subplots(nRows, nCols)
    f.set_size_inches(19.8, 10.8)

    for i in range(nRows):
        k = i + nRows
        axArr[i, 0].set_ylabel('PC%d' % k)
        for j in range(nCols):
            corrPlot(colArr[j], P.tomo_obs_norm__zk[:, k], axArr[i, j])

    f.subplots_adjust(hspace = 0.0)
    f.subplots_adjust(wspace = 0.0)

    plt.setp([a.get_xticklabels() for a in f.axes[:-5]], visible = False)
    plt.setp([a.get_yticklabels() for a in f.axes[::5]], visible = True)
    axArr[0, 0].set_title('at_flux')
    axArr[0, 1].set_title('aZ_flux')
    axArr[0, 2].set_title('A_V')
    axArr[0, 3].set_title('v_0')
    axArr[0, 4].set_title('v_d')

    #plt.legend()
    plt.suptitle('Correlations PC10 ... PC19 - SYN NORM')
    #plt.tight_layout(pad = 2., w_pad = 0.5, h_pad = 0.5)
    f.savefig('%s-corre_syn_norm_10-19.png' % califaID)
    plt.close()


############################### RES NORM ###############################     
    P.PCA_res_norm()
    P.tomograms_res_norm()

    f, axArr = plt.subplots(nRows, nCols)
    f.set_size_inches(19.8, 10.8)

    for i in range(nRows):
        axArr[i, 0].set_ylabel('PC%d' % i)
        for j in range(nCols):
            corrPlot(colArr[j], P.tomo_res_norm__zk[:, i], axArr[i, j])

    f.subplots_adjust(hspace = 0.0)
    f.subplots_adjust(wspace = 0.0)

    plt.setp([a.get_xticklabels() for a in f.axes[:-5]], visible = False)
    plt.setp([a.get_yticklabels() for a in f.axes[::5]], visible = True)
    axArr[0, 0].set_title('at_flux')
    axArr[0, 1].set_title('aZ_flux')
    axArr[0, 2].set_title('A_V')
    axArr[0, 3].set_title('v_0')
    axArr[0, 4].set_title('v_d')

    #plt.legend()
    plt.suptitle('Correlations PC0 ... PC9 - RES NORM')
    #plt.tight_layout(pad = 2., w_pad = 0.5, h_pad = 0.5)
    f.savefig('%s-corre_res_norm_0-9.png' % califaID)
    plt.close()

    f, axArr = plt.subplots(nRows, nCols)
    f.set_size_inches(19.8, 10.8)

    for i in range(nRows):
        k = i + nRows
        axArr[i, 0].set_ylabel('PC%d' % k)
        for j in range(nCols):
            corrPlot(colArr[j], P.tomo_res_norm__zk[:, k], axArr[i, j])

    f.subplots_adjust(hspace = 0.0)
    f.subplots_adjust(wspace = 0.0)

    plt.setp([a.get_xticklabels() for a in f.axes[:-5]], visible = False)
    plt.setp([a.get_yticklabels() for a in f.axes[::5]], visible = True)
    axArr[0, 0].set_title('at_flux')
    axArr[0, 1].set_title('aZ_flux')
    axArr[0, 2].set_title('A_V')
    axArr[0, 3].set_title('v_0')
    axArr[0, 4].set_title('v_d')

    #plt.legend()
    plt.suptitle('Correlations PC10 ... PC19 - RES NORM')
    #plt.tight_layout(pad = 2., w_pad = 0.5, h_pad = 0.5)
    f.savefig('%s-corre_res_norm_10-19.png' % califaID)
    plt.close()

