'''
Created on 17/04/2013

@author: lacerda
'''
import matplotlib
matplotlib.use('agg')
import numpy as np
import atpy
from pycasso.fitsdatacube import fitsQ3DataCube
from pystarlight import io
from scipy import linalg
from scipy import stats as st
from matplotlib import gridspec
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator

# plt.rcParams.update({'font.family' : 'Times New Roman',
#                     'text.usetex' : False,
#                     'backend' : 'ps'
#                     })

fitsDirDefault = '/home/lacerda/CALIFA/gal_fits'
quantilFlagDefault = 0.9
fitsFilenameSuffix = '_synthesis_eBR_v20_q027.d13c512.ps3b.k1.mC.CCM.Bgsd01.v01.fits'

class PCAlifa:
    def __init__(self, califaID = False, fitsDir = fitsDirDefault, quantilQFlag = quantilFlagDefault, lc = []):
        self.califaID = califaID
        self._initVars()

        if califaID:
            self.readCALIFACube(self.califaID, fitsDir)

            if len(lc):
                self.setLambdaConstrains(lc)

            self.quantilQFlag = quantilQFlag

            if quantilQFlag:
                self.setQFlag(quantilQFlag)

    def readCALIFACube(self, califaID, fitsDir = fitsDirDefault):
        self.califaID = califaID

        self._initVars()

        self.fitsDir = fitsDir
        self.fitsFile = '%s/%s/%s%s' % (fitsDir, califaID, califaID, fitsFilenameSuffix)

        self.K = fitsQ3DataCube(self.fitsFile)

        self.maskEmLines = np.ones_like(self.K.l_obs, dtype = np.bool)
        self.maskQFlag = self.maskEmLines
        self.maskLambdaConstrains = self.maskEmLines

        self._setVars()

#    def pixDevLineFluxRestFrame(self, modelLineFluxMax = 6562.85, searchRange = 30):
#        self.restFrameLine = modelLineFluxMax
#
#        K = self.K
#
#        modelLineFluxMax = 6562.85
#        searchRange = 10.
#        l_obs = K.l_obs
#        res__lz = (K.f_obs - K.f_syn)
#        llow = modelLineFluxMax - searchRange
#        lup = modelLineFluxMax + searchRange
#        lpixlow = np.where(l_obs >= llow)[0][0]
#        lpixup = np.where(l_obs <= lup)[0][-1]
#        line_max_lpix__z = lpixlow + res__lz[lpixlow:lpixup, :].argmax(axis = 0)
#        l_step = l_obs[1] - l_obs[0]
#        dim = np.ones((len(l_obs), K.N_zone)).T
#        lmax__z = ((l_obs * dim).T)[line_max_lpix__z][:, 0]
#        dpix__z = np.array(((lmax__z - modelLineFluxMax) / l_step))
#        self.restFramePixelDeviance__z = dpix__z
#        self.restFrameLambdaDeviance__z = lmax__z - modelLineFluxMax
#
#    def setLineFluxRestFrame(self):
#        K = self.K
#        interp_step = 0.01
#        dpix__z = np.array((1. / interp_step * self.restFramePixelDeviance__z), dtype = np.int)
#
#        l_obs = K.l_obs
#        l_obs_pix = np.arange(0, len(l_obs), 1)
#        l_interp_pix = np.arange(0, len(l_obs), 0.01)
#
#        interp_f_obs__lz = np.zeros((len(l_interp_pix), K.N_zone))
#        interp_f_syn__lz = np.zeros((len(l_interp_pix), K.N_zone))
#        interp_f_obs_notroll__lz = np.zeros((len(l_interp_pix), K.N_zone))
#        interp_f_syn_notroll__lz = np.zeros((len(l_interp_pix), K.N_zone))
#
#        f_obs_rf__lz = np.zeros((len(l_obs_pix), K.N_zone))
#        f_syn_rf__lz = np.zeros((len(l_obs_pix), K.N_zone))
#
#        for z in range(K.N_zone):
#            dev = dpix__z[z]
#            interp_f_obs_notroll__lz[:, z] = np.interp(l_interp_pix, l_obs_pix, K.f_obs[:, z])
#            interp_f_syn_notroll__lz[:, z] = np.interp(l_interp_pix, l_obs_pix, K.f_syn[:, z])
#
#            interp_f_obs__lz[:, z] = np.roll(interp_f_obs_notroll__lz[:, z], dev)
#            interp_f_syn__lz[:, z] = np.roll(interp_f_syn_notroll__lz[:, z], dev)
#
#            if dev < 0:
#                interp_f_obs__lz[dev:, z] = interp_f_obs_notroll__lz[dev:, z].mean()
#                interp_f_syn__lz[dev:, z] = interp_f_syn_notroll__lz[dev:, z].mean()
#
#            f_obs_rf__lz[:, z] = np.interp(l_obs_pix, l_interp_pix, interp_f_obs__lz[:, z])
#            f_syn_rf__lz[:, z] = np.interp(l_obs_pix, l_interp_pix, interp_f_syn__lz[:, z])
#
#        self.f_obs_rf__zl = f_obs_rf__lz.T
#        self.f_syn_rf__zl = f_syn_rf__lz.T

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

    def tomogram(self, I, eigVec):
        t__zk = np.dot(I, eigVec)
        t__kyx = False

        if self.califaID:
            t__kyx = self.K.zoneToYX(t__zk.T, extensive = False)

        return t__zk, t__kyx

    def _setVars(self):
        mask = self.maskEmLines & self.maskLambdaConstrains & self.maskQFlag
        synMask = self.maskLambdaConstrains & self.maskQFlag

        self.l_obs = self.K.l_obs[mask]
        self.l_syn = self.K.l_obs[synMask]
        self.f_obs__zl = self.K.f_obs[mask].transpose()
        self.f_obs_norm__zl = (self.K.f_obs[mask] / self.K.fobs_norm).transpose()
        self.f_syn__zl = self.K.f_syn[synMask].transpose()
        self.f_syn_norm__zl = (self.K.f_syn[synMask] / self.K.fobs_norm).transpose()
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

        self.maskLambdaConstrains = (self.K.l_obs >= ldown) & (self.K.l_obs <= lup)

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
        mask = (self.K.l_obs >= t[0]['l_up'])

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

    def zoneRebuildSpecAxisPlot(self, ax, l, O, R, eVal, eVec, eVMask, npref, fontsize, resid):
        res = O - R
        sigmaNReb = 0.
        sigmaReb = np.sqrt(eVal[eVMask].sum())
        varReb = eVal[eVMask].sum() / eVal.sum()

        eVMask_not_used = np.asarray([not x for x in eVMask])

        if eVMask_not_used.any():
            sigmaNReb = np.sqrt(eVal[eVMask_not_used].sum())
            varNReb = eVal[eVMask_not_used].sum() / eVal.sum()

        sigmaRatio = sigmaNReb / sigmaReb

        textStrAdev = ''
        res_plot = res
        res_label = 'res'

        if resid:
            ax.set_ylim([1.1 * O.min(), 1.1 * O.max()])
        else:
            res_plot = 5. * res
            res_label = r'$5\ \times$ res'
            ax.set_ylim([-0.5 * O.mean(), 1.5 * O.mean()])
            adev = 100. * (1. / len(l)) * (np.abs(res) / O).sum()
            textStrAdev = 'adev\ =\ %.4f\ \%%' % adev

        textStr = r'$%s \ \ \sigma_{reb}\ =\ %s\ \ \sigma_{Nreb}\ =\ %s\ \ ratio\ =\ %.2f$' % (textStrAdev, self.nToStrSciNot(sigmaReb), self.nToStrSciNot(sigmaNReb), sigmaRatio)

        ax.plot(l, O, label = 'Obs')
        ax.plot(l, R, label = 'Mod')
        ax.plot(l, res_plot, label = res_label)

        ax.text(0.05, 0.92, textStr,
                fontsize = fontsize + 1, transform = ax.transAxes,
                horizontalalignment = 'left',
                verticalalignment = 'center',
                multialignment = 'left')
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

    def tomoPlot(self, t, l, eigvec, eigval, ms, ti, npref):
        tomogram = t[ti, :, :]
        x = l
        y = eigvec[:, ti]
        y2 = ms

        fig = plt.figure(figsize = (15, 5))
        gs = gridspec.GridSpec(1, 2, width_ratios = [4, 7])
        ax1 = plt.subplot(gs[0])
        ax2 = plt.subplot(gs[1])

        ax1.set_title(r'tomogram $%02i$' % ti)
        im = ax1.imshow(tomogram, origin = 'lower', interpolation = 'nearest', aspect = 'auto')
        fig.colorbar(ax = ax1, mappable = im, use_gridspec = True)

        eigval_norm = 100. * eigval[ti] / eigval.sum()
        ax2.set_title(r'Eigenvalue $%s$ ($%.2f\ \%%$)' % (self.nToStrSciNot(eigval[ti]), eigval_norm))
        ax2.plot(x, y)
        ax2.xaxis.set_major_locator(MaxNLocator(20))
        ax2.set_ylabel(r'$PC %02i$' % ti)
        ax2.grid()

        ax3 = ax2.twinx()
        ax3.plot(x, y2, color = '0.55')
        ax3.set_ylabel(r'Mean spetrum')

#        bins = np.arange(0, 3 + 0.1, 0.1)
#        bin_center = (bins[1:] + bins[:-1]) / 2.0
#        pc__HLR = self.K.radialProfile(tomogram, bins)
#        ax3 = plt.axes([0.2, 0.8, .2, .2], axisbg = 'y')
#        ax3.plot(bin_center, pc__HLR)
#        ax3.set_title(r'Radial Profile')
#        ax3.set_xlabel(r'$R_{50}$')

        plt.tight_layout()
        fig.savefig('%stomo_%02i.png' % (npref, ti))
        plt.close()

    def screeTestPlot(self, eigval, maxInd, npref, title):
        fig = plt.figure(figsize = (8, 6))
        eigval_norm = 100. * eigval / eigval.sum()
        plt.plot(eigval_norm[:maxInd], linestyle = '-', marker = '*')
        plt.ylim([0, eigval_norm[1] * 1.1])
        plt.xticks(range(maxInd))
        plt.title(r'%s Scree test' % title)
        plt.xlabel(r'Principal component')
        plt.ylabel(r'Normalized variance [$\%%$]')
        plt.grid()
        fig.savefig('%sscree.png' % npref)
        plt.close()

    def nToStrSciNot(self, n):
        e = np.floor(np.log10(np.abs(n)))
        m = n / 10.**(e)

        if np.int(e) == 0:
            nStr = r'%.2f' % n
        else:
            nStr = r'%.2f \times 10^{%d}' % (m, e)

        return nStr

    def sanityCheck(self, x, y, npref):
        f = plt.figure()
        f.set_size_inches(19.2, 10.8)
        plt.plot(x, y)
        plt.plot(x, y, 'r.')
        plt.ylabel(r'Obs. flux $[erg/s/cm^2/\AA]$')
        plt.xlabel(r'Wavelength $[\AA]$')
        plt.title(r'Mask check (zone 0)')
        plt.grid()
        f.savefig('%ssanity_check.png' % npref)


