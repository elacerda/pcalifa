'''
Created on 17/04/2013

@author: lacerda
'''
from pycasso.fitsdatacube import fitsQ3DataCube
import pystarlight.io
import numpy as np
import atpy
from scipy import linalg
import scipy.stats as st
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MaxNLocator
from matplotlib import pyplot as plt

fitsDirDefault = '/home/lacerda/CALIFA'

class PCAlifa:
    def __init__(self, califaID = False, fitsDir = fitsDirDefault, flagLinesQuantil = 0.9):
        self.histo = None
        self.tStarlight = None
        self.starlightMaskFile = None
        self.califaID = califaID
        self.initVars()

        if califaID:
            self.flagLinesQuantil = flagLinesQuantil
            self.fitsDir = fitsDir
            self.fitsFile = '%s/%s/%s_synthesis_eBR_v20_q027.d13c512.ps3b.k1.mC.CCM.Bgsd01.v01.fits' % (fitsDir, califaID, califaID)
            self.K = fitsQ3DataCube(self.fitsFile)
            self.mask = np.ones_like(self.K.l_obs, dtype = np.bool)
            self.maskEmLines = self.mask

            self.setVars()

            if flagLinesQuantil:
                self.removeFlaggedLambda(flagLinesQuantil)

    def PCA(self, arr, num, axis = -1, arrMean = None):
        if arrMean == None:
            arrMean = arr.mean(axis = axis)

        diff = arr - arrMean
        covMat = np.dot(diff.T, diff) / (num - 1.)
        w, e = linalg.eigh(covMat)

        S = np.argsort(w)[::-1]
        wS = w[S]
        eS = e[:, S]

        return diff, arrMean, covMat, wS, eS

    def PCA_obs(self):
        self.I_obs__zl, self.ms_obs__l, self.covMat_obs__ll, self.eigVal_obs__k, self.eigVec_obs__lk = self.PCA(self.f_obs__zl, self.K.N_zone, 0)

    def PCA_obs_norm(self):
        self.I_obs_norm__zl, self.ms_obs_norm__l, self.covMat_obs_norm__ll, self.eigVal_obs_norm__k, self.eigVec_obs_norm__lk = self.PCA(self.f_obs_norm__zl, self.K.N_zone, 0)

    def PCA_syn(self):
        self.I_syn__zl, self.ms_syn__l, self.covMat_syn__ll, self.eigVal_syn__k, self.eigVec_syn__lk = self.PCA(self.f_syn__zl, self.K.N_zone, 0)

    def PCA_syn_norm(self):
        self.I_syn_norm__zl, self.ms_syn_norm__l, self.covMat_syn_norm__ll, self.eigVal_syn_norm__k, self.eigVec_syn_norm__lk = self.PCA(self.f_syn_norm__zl, self.K.N_zone, 0)

    def PCA_res(self):
        self.I_res__zl, self.ms_res__l, self.covMat_res__ll, self.eigVal_res__k, self.eigVec_res__lk = self.PCA(self.f_res__zl, self.K.N_zone, 0)

    def PCA_res_norm(self):
        self.I_res_norm__zl, self.ms_res_norm__l, self.covMat_res_norm__ll, self.eigVal_res_norm__k, self.eigVec_res_norm__lk = self.PCA(self.f_res_norm__zl, self.K.N_zone, 0)

    def tomograms_obs(self):
        self.tomo_obs__zk, self.tomo_obs__kyx = self.tomogram(self.I_obs__zl, self.eigVec_obs__lk)

    def tomograms_obs_norm(self):
        self.tomo_obs_norm__zk, self.tomo_obs_norm__kyx = self.tomogram(self.I_obs_norm__zl, self.eigVec_obs_norm__lk)

    def tomograms_syn(self):
        self.tomo_syn__zk, self.tomo_syn__kyx = self.tomogram(self.I_syn__zl, self.eigVec_syn__lk)

    def tomograms_syn_norm(self):
        self.tomo_syn_norm__zk, self.tomo_syn_norm__kyx = self.tomogram(self.I_syn_norm__zl, self.eigVec_syn_norm__lk)

    def tomograms_res(self):
        self.tomo_res__zk, self.tomo_res__kyx = self.tomogram(self.I_res__zl, self.eigVec_res__lk)

    def tomograms_res_norm(self):
        self.tomo_res_norm__zk, self.tomo_res_norm__kyx = self.tomogram(self.I_res_norm__zl, self.eigVec_res_norm__lk)

    def runPCA(self):
        self.PCA_obs()
        self.PCA_obs_norm()
        self.PCA_res()
        self.PCA_res_norm()
        self.PCA_syn()
        self.PCA_syn_norm()

        self.tomograms_obs()
        self.tomograms_obs_norm()
        self.tomograms_res()
        self.tomograms_res_norm()
        self.tomograms_syn()
        self.tomograms_syn_norm()

    def delPCA(self):
        del self.I_obs__zl, self.ms_obs, self.covMat_obs__ll, self.eigVal_obs__k, self.eigVec_obs__lk
        del self.I_obs_norm__zl, self.ms_obs_norm, self.covMat_obs_norm__ll, self.eigVal_obs_norm__k, self.eigVec_obs_norm__lk
        del self.I_syn__zl, self.ms_syn, self.covMat_syn__ll, self.eigVal_syn__k, self.eigVec_syn__lk
        del self.I_syn_norm__zl, self.ms_syn_norm, self.covMat_syn_norm__ll, self.eigVal_syn_norm__k, self.eigVec_syn_norm__lk
        del self.I_res__zl, self.ms_res, self.covMat_res__ll, self.eigVal_res__k, self.eigVec_res__lk
        del self.I_res_norm__zl, self.ms_res_norm, self.covMat_res_norm__ll, self.eigVal_res_norm__k, self.eigVec_res_norm__lk

    def delTomograms(self):
        del self.tomo_obs__zk, self.tomo_obs__kyx
        del self.tomo_obs_norm_zk, self.tomo_obs_norm__kyx
        del self.tomo_syn__zk, self.tomo_syn__kyx
        del self.tomo_syn_norm__zk, self.tomo_syn_norm__kyx
        del self.tomo_res__zk, self.tomo_res__kyx
        del self.tomo_res_norm__zk, self.tomo_res_norm__kyx

    def tomogram(self, I, eigVec, extensive = True):
        t__zk = np.dot(I, eigVec)
        t__kyx = False

        if self.califaID:
            t__kyx = self.K.zoneToYX(t__zk.T, extensive = extensive)

        return t__zk, t__kyx

    def setVars(self):
        self.l_obs = self.K.l_obs[self.mask]
        self.f_obs__zl = self.K.f_obs[self.mask].transpose()
        self.f_obs_norm__zl = (self.K.f_obs[self.mask] / self.K.fobs_norm).transpose()
        self.f_syn__zl = self.K.f_syn[self.mask].transpose()
        self.f_syn_norm__zl = (self.K.f_syn[self.mask] / self.K.fobs_norm).transpose()
        self.f_res__zl = (self.K.f_obs[self.mask] - self.K.f_syn[self.mask]).transpose()
        self.f_res_norm__zl = ((self.K.f_obs[self.mask] - self.K.f_syn[self.mask]) / self.K.fobs_norm).transpose()

    def initVars(self):
        self.l_obs = False
        self.f_obs__zl = False
        self.f_obs_norm__zl = False
        self.f_syn__zl = False
        self.f_syn_norm__zl = False
        self.f_res__zl = False
        self.f_res_norm__zl = False

        self.I_obs__zl = False
        self.ms_obs = False
        self.covMat_obs__ll = False
        self.eigVal_obs__k = False
        self.eigVec_obs__lk = False
        self.tomo_obs__zk = False
        self.tomo_obs__kyx = False

        self.I_obs_norm__zl = False
        self.ms_obs_norm = False
        self.covMat_obs_norm__ll = False
        self.eigVal_obs_norm__k = False
        self.eigVec_obs_norm__lk = False
        self.tomo_obs_norm_zk = False
        self.tomo_obs_norm__kyx = False

        self.I_syn__zl = False
        self.ms_syn = False
        self.covMat_syn__ll = False
        self.eigVal_syn__k = False
        self.eigVec_syn__lk = False
        self.tomo_syn__zk = False
        self.tomo_syn__kyx = False

        self.I_syn_norm__zl = False
        self.ms_syn_norm = False
        self.covMat_syn_norm__ll = False
        self.eigVal_syn_norm__k = False
        self.eigVec_syn_norm__lk = False
        self.tomo_syn_norm__zk = False
        self.tomo_syn_norm__kyx = False

        self.I_res__zl = False
        self.ms_res = False
        self.covMat_res__ll = False
        self.eigVal_res__k = False
        self.eigVec_res__lk = False
        self.tomo_res__zk = False
        self.tomo_res__kyx = False

        self.I_res_norm__zl = False
        self.ms_res_norm = False
        self.covMat_res_norm__ll = False
        self.eigVal_res_norm__k = False
        self.eigVec_res_norm__lk = False
        self.tomo_res_norm__zk = False
        self.tomo_res_norm__kyx = False

    def removeFlaggedLambda(self, quantil):
        self.flagLinesQuantil = quantil
        self.histo = self.K.f_flag.sum(axis = 1) / self.K.N_zone
        self.maskHisto = (self.histo < quantil)
        self.mask = self.maskHisto & (self.K.l_obs > 3800) & (self.K.l_obs < 6850)
        self.setVars()

    def removeStarlightEmLines(self, maskFile):
        self.starlightMaskFile = maskFile
        t = atpy.Table(maskfile = maskFile, type = 'starlight_mask')
        mask = (self.K.l_obs > t[0]['l_up'])

        for i in range(1, len(t)):
            if (t[i]['weight'] == 0.0):
                mask = mask & ((self.K.l_obs < t[i]['l_low']) | (self.K.l_obs > t[i]['l_up']))

        self.tStarlight = t
        self.maskEmLines = mask
        self.mask = self.mask & self.maskEmLines
        self.setVars()

    def rebuildSpectra(self, tomo, eigVec, mean):
        I_rec = np.dot(tomo, eigVec.transpose())
        f_rec = I_rec + mean

        return I_rec, f_rec

    ''' 
    TODO: 
        adicionar P.l_obs na chamada da funcao em PCAlifa-rebuildSpec.py 
    '''
    def zoneRebuildPlot(self, iZone, evRebArr, l, O, tomo, eVec, eVal, mean, nPref, resid = False):
        i = 0
        nCols = 2
        nRows = len(evRebArr) / (1. * nCols)

        if nRows > int(nRows):
            nRows = nRows + 1

        nRows = int(nRows)

        fig = plt.figure(figsize = (19.8, 10.8))
        gs = gridspec.GridSpec(nRows, nCols, width_ratios = [1, 1], hspace = 0)

        adev = np.zeros_like(evRebArr, dtype = np.float64)

        for i, ne in enumerate(evRebArr):
            ax = plt.subplot(gs[i])

            I_reb__zl, M = self.rebuildSpectra(tomo[:, :ne], eVec[:, :ne], mean)
            diff = O[iZone, :] - M[iZone, :]
            sigmaNReb = 0.
            sigmaReb = np.sqrt(eVal[:ne].sum())

            if (ne + 1 != len(eVec)):
                sigmaNReb = np.sqrt(eVal[ne + 1:].sum())

            sigmaRatio = sigmaNReb / sigmaReb
            textStrAdev = ''

            if (resid == True):
                ax.set_ylim([1.1 * O[iZone, :].min(), 1.1 * O[iZone, :].max()])
            else:
                ax.set_ylim([-0.5 * O[iZone, :].mean(), 1.5 * O[iZone, :].mean()])
                adev[i] = 100. * (1. / len(l)) * (np.abs(diff) / O[iZone, :]).sum()
                textStrAdev = 'adev =  %.4f %% - ' % adev[i]

            textStr = '%ssigmaReb = %.2e - sigmaNReb = %.2e - ratio = %.4f' % (textStrAdev, sigmaReb, sigmaNReb, sigmaRatio)

            ax.plot(l, O[iZone, :], label = 'Obs')
            ax.plot(l, M[iZone, :], label = 'Mod')
            ax.plot(l, diff * 5., label = 'Res x 5')
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

        suptitle_txt = ''

        if self.califaID:
            suptitle_txt = 'CALIFA ID: %s%s' % (self.califaID, ' ' * 60)

        plt.suptitle('%sZone %04d' % (suptitle_txt, iZone))
        plt.tight_layout(pad = 2., w_pad = 0., h_pad = 0.)
        fig.savefig('%s-iZone-%04d.png' % (nPref, iZone))
        plt.close()

        if (resid == False):
            fig = plt.figure(figsize = (19.8, 10.8))
            plt.plot(evRebArr, adev, 'o', label = 'adev')
            plt.legend()
            plt.xlabel('Number of eigenvectors')
            plt.ylabel('adev')
            plt.grid()
            fig.savefig('%s-iZone-%04d-adev.png' % (nPref, iZone))
            plt.close()

    def corrPlot(self, x, y, ax):
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

    def tomoPlot(self, tn, l, t, eigvec, eigval, npref):
        fig = plt.figure(figsize = (15, 5))
        gs = gridspec.GridSpec(1, 2, width_ratios = [4, 7])
        ax1 = plt.subplot(gs[0])
        ax2 = plt.subplot(gs[1])
        ax1.set_title('tomogram %02i' % tn)
        im = ax1.imshow(t[tn, :, :], origin = 'lower', interpolation = 'nearest', aspect = 'auto')
        fig.colorbar(ax = ax1, mappable = im, use_gridspec = True)
        ax2.set_title('eigval %.4e' % eigval[tn])
        ax2.plot(l, eigvec[:, tn])
        ax2.xaxis.set_major_locator(MaxNLocator(20))
        ax2.grid()
        plt.tight_layout()
        fig.savefig('%stomo%02i.png' % (npref, tn))
        plt.close()

    def screeTestPlot(self, eigval, maxInd, npref):
        fig = plt.figure(figsize = (8, 6))
        eigval_norm = eigval / eigval.sum()
        plt.plot(eigval_norm[:maxInd], linestyle = '-', marker = '*')
        plt.ylim([0, eigval_norm[1] * 1.1])
        plt.xticks(range(maxInd))
        plt.title('%s scree test' % npref)
        plt.grid()
        fig.savefig('%sscree.png' % npref)
        plt.close()

