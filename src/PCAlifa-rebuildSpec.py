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

#fitsDir = sys.argv[1]
#califaID = sys.argv[2]
fitsDir = '/home/lacerda/CALIFA'
califaID = 'K0277'
maskfile = '/home/lacerda/workspace/PCA/src/Mask.mC'
flagLinesQuantil = 0.9
remFlaggedLambdas = True
remStarlightEmLines = False
tmax = 20 # numero maximo de eigenvalues

if __name__ == '__main__':
    P = PCA.PCAlifa(califaID = califaID,
                    fitsDir = fitsDir,
                    flagLinesQuantil = flagLinesQuantil,
                    remFlaggedLambdas = remFlaggedLambdas,
                    runDefaultPCA = False)

    if (remStarlightEmLines == True):
        P.removeStarlightEmLines(maskfile)

    P.PCA_obs()
    P.tomograms_obs()

    for zone in [0, 10, 20, 100, 200, 500]:
        i = 0

        fig = plt.figure(figsize = (19.8, 10.8))
        gs = gridspec.GridSpec(6, 2, width_ratios = [1, 1])
        #gs.update(hspace = 0.05)

        for n in [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 15, 20]:
            ax = plt.subplot(gs[i])

            I_obs_reb__zl, f_obs_reb__zl = P.rebuildSpectra(P.tomo_obs__zk, P.eigVec_obs__lk, P.ms_obs, n)
            diff = P.f_obs__zl[zone, :] - f_obs_reb__zl[zone, :]
            adev = 100. * (1. / P.K.N_zone) * (np.abs(diff) / P.f_obs__zl[zone, :]).sum()

            ax.plot(P.l_obs, P.f_obs__zl[zone, :], label = 'Obs')
            ax.plot(P.l_obs, f_obs_reb__zl[zone, :], label = 'Mod')
            ax.plot(P.l_obs, diff, label = 'Res')
            ax.xaxis.set_major_locator(MaxNLocator(20))
            ax.set_ylabel('n_eigvec = %d' % n)
            ax.set_title('adev = %.4f' % adev)
            ax.legend(prop = {'size':7})
            #delrms = np.sqrt()
            ax.grid()

            if (i < 10):
                plt.setp(ax.get_xticklabels(), visible = False)

            i = i + 1

        plt.suptitle('CALIFA ID: %s - Zone %d' % (P.K.califaID, zone))
        plt.tight_layout(pad = 2., w_pad = 0.05, h_pad = 0.6)
        fig.savefig('zone-%d.png' % zone)
        plt.close()
