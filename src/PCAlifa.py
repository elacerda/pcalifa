'''
Created on 17/04/2013

@author: lacerda
'''
import matplotlib
matplotlib.use('agg')

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
quantilFlagDefault = 0.9
fitsFilenameSuffix = '_synthesis_eBR_v20_q027.d13c512.ps3b.k1.mC.CCM.Bgsd01.v01.fits'

class PCAlifa:
    def __init__(self, califaID = False, fitsDir = fitsDirDefault, quantilQFlag = quantilFlagDefault, lc = []):
        self.califaID = califaID
        self._initVars()

        if califaID:
            self.readCALIFACube(self.califaID, fitsDir, quantilQFlag, lc)

    def readCALIFACube(self, califaID, fitsDir = fitsDirDefault, quantilQFlag = quantilFlagDefault, lc = []):
        self.califaID = califaID

        self._initVars()

        self.quantilQFlag = quantilQFlag
        self.fitsDir = fitsDir
        self.fitsFile = '%s/%s/%s%s' % (fitsDir, califaID, califaID, fitsFilenameSuffix)

        self.K = fitsQ3DataCube(self.fitsFile)

        self.maskEmLines = np.ones_like(self.K.l_obs, dtype = np.bool)
        self.maskQFlag = self.maskEmLines
        self.maskLambdaConstrains = self.maskEmLines

        self._setVars()

        if len(lc):
            self.setLambdaConstrains(lc)

        if quantilQFlag:
            self.setQFlag(quantilQFlag)

    def PCA(self, arr, num, axis = -1, arrMean = False):
        if not arrMean:
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

    def tomogram(self, I, eigVec, extensive = True):
        t__zk = np.dot(I, eigVec)
        t__kyx = False

        if self.califaID:
            t__kyx = self.K.zoneToYX(t__zk.T, extensive = extensive)

        return t__zk, t__kyx

    def _setVars(self):
        mask = self.maskEmLines & self.maskLambdaConstrains & self.maskQFlag
        self.l_obs = self.K.l_obs[mask]
        self.f_obs__zl = self.K.f_obs[mask].transpose()
        self.f_obs_norm__zl = (self.K.f_obs[mask] / self.K.fobs_norm).transpose()
        self.f_syn__zl = self.K.f_syn[mask].transpose()
        self.f_syn_norm__zl = (self.K.f_syn[mask] / self.K.fobs_norm).transpose()
        self.f_res__zl = (self.K.f_obs[mask] - self.K.f_syn[mask]).transpose()
        self.f_res_norm__zl = ((self.K.f_obs[mask] - self.K.f_syn[mask]) / self.K.fobs_norm).transpose()

    def _initVars(self):
        self.maskEmLines = []
        self.maskQFlag = []
        self.maskLambdaConstrains = []

        self.histo = False
        self.tStarlight = False
        self.starlightMaskFile = False

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

    def setLambdaConstrains(self, lc):
        lc = np.array(lc)
        s = np.argsort(lc)
        ldown = lc[s][0]
        lup = lc[s][1]

        self.maskLambdaConstrains = (self.K.l_obs > ldown) & (self.K.l_obs < lup)

        self._setVars()

    def unsetLambdaConstrains(self):
        self.maskLambdaConstrains = np.ones_like(self.K.l_obs, dtype = np.bool)
        self._setVars()

    def setQFlag(self, quantil):
        self.quantilQFlag = quantil
        self.histo = self.K.f_flag.sum(axis = 1) / self.K.N_zone
        self.maskQFlag = (self.histo < quantil)
        self._setVars()

    def unsetQFlag(self):
        self.quantilQFlag = False
        self.histo = False
        self.maskQFlag = np.ones_like(self.K.l_obs, dtype = np.bool)
        self._setVars()

    def setStarlightMaskFile(self, maskFile):
        self.starlightMaskFile = maskFile
        t = atpy.Table(maskfile = maskFile, type = 'starlight_mask')
        mask = (self.K.l_obs > t[0]['l_up'])

        for i in range(1, len(t)):
            if (t[i]['weight'] == 0.0):
                mask = mask & ((self.K.l_obs < t[i]['l_low']) | (self.K.l_obs > t[i]['l_up']))

        self.tStarlight = t
        self.maskEmLines = mask
        self._setVars()

    def unsetStarlightMaskFile(self):
        self.starlightMaskFile = None
        del self.tStarlight
        self.tStarlight = None
        self.maskEmLines = np.ones_like(self.K.l_obs, dtype = np.bool)
        self._setVars()

    def rebuildSpectra(self, tomo, eigVec, mean):
        I_rec = np.dot(tomo, eigVec.transpose())
        f_rec = I_rec + mean

        return I_rec, f_rec

    def zoneRebuildSpecAxisPlot(self, ax, l, O, R, eVal, eVec, eVMask, npref, fontsize = 7, resid = False):
        ''' criando uma string com os eigenvectors usados para reconstruir o cubo'''
        res = O - R
        adev = 100. * (1. / len(l)) * (np.abs(res) / O).sum()

        sigmaNReb = 0.
        sigmaReb = np.sqrt(eVal[eVMask].sum())

        eVMask_not_used = np.asarray([not x for x in eVMask])

        if eVMask_not_used.any():
            sigmaNReb = np.sqrt(eVal[eVMask_not_used].sum())

        sigmaRatio = sigmaNReb / sigmaReb

        textStrAdev = ''

        if resid:
            ax.set_ylim([1.1 * O.min(), 1.1 * O.max()])
        else:
            ax.set_ylim([-0.5 * O.mean(), 1.5 * O.mean()])
            adev = 100. * (1. / len(l)) * (np.abs(res) / O).sum()
            textStrAdev = 'adev =  %.4f %% - ' % adev

        textStr = '%ssigmaReb = %.2e - sigmaNReb = %.2e - ratio = %.4f' % (textStrAdev, sigmaReb, sigmaNReb, sigmaRatio)

        ax.plot(l, O, label = 'Obs')
        ax.plot(l, R, label = 'Mod')
        ax.plot(l, res, label = 'res')
        ax.text(0.01, 0.92, textStr, fontsize = fontsize + 3, transform = ax.transAxes,
                horizontalalignment = 'left', verticalalignment = 'center', multialignment = 'left')
        ax.legend(prop = {'size' : fontsize})
        ax.grid()

    def correlationAxisPlot(self, x, y, ax):
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
        fig.savefig('%stomo_%02i.png' % (npref, tn))
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

