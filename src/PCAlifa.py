'''
Created on 17/04/2013

@author: lacerda
'''
from pycasso.fitsdatacube import fitsQ3DataCube
import pystarlight.io
import numpy as np
import atpy
from scipy import linalg

fitsDirDefault = '/home/lacerda/CALIFA'

class PCAlifa:
    def __init__(self, califaID, fitsDir = fitsDirDefault, flagLinesQuantil = 0.9):
        self.histo = None
        self.tStarlight = None
        self.starlightMaskFile = None
        self.califaID = califaID
        self.flagLinesQuantil = flagLinesQuantil
        self.fitsDir = fitsDir
        self.fitsFile = '%s/%s/%s_synthesis_eBR_v20_q027.d13c512.ps3b.k1.mC.CCM.Bgsd01.v01.fits' % (fitsDir, califaID, califaID)
        self.K = fitsQ3DataCube(self.fitsFile)
        self.mask = np.ones_like(self.K.l_obs, dtype = np.bool)
        self.maskEmLines = self.mask

        self.initVars()

        if flagLinesQuantil:
            self.removeFlaggedLambda(flagLinesQuantil)

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

    def PCA(self, arr, num, axis = -1, arrMean = None):
        if arrMean == None:
            arrMean = arr.mean(axis = axis)

        diff = arr - arrMean
        covMat = np.dot(diff.T, diff) / (num - 1.)
        w, e = linalg.eig(covMat)

        S = np.argsort(w)[::-1]
        wS = w[S]
        eS = e[:, S]

        return diff, arrMean, covMat, np.real(wS), np.real(eS)

    def PCAins(self, arr, num, axis = -1, arrMean = None):
        if arrMean == None:
            arrMean = arr.mean(axis = axis)

        diff = arr - arrMean
        covMat = np.dot(diff.T, diff) / (num - 1.)
        w, e = linalg.eig(covMat)

        return diff, arrMean, covMat, w, e

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
        t__kyx = self.K.zoneToYX(t__zk.T, extensive = extensive)

        return t__zk, t__kyx

    def initVars(self):
        self.l_obs = self.K.l_obs[self.mask]
        self.f_obs__zl = self.K.f_obs[self.mask].transpose()
        self.f_obs_norm__zl = (self.K.f_obs[self.mask] / self.K.fobs_norm).transpose()
        self.f_syn__zl = self.K.f_syn[self.mask].transpose()
        self.f_syn_norm__zl = (self.K.f_syn[self.mask] / self.K.fobs_norm).transpose()
        self.f_res__zl = (self.K.f_obs[self.mask] - self.K.f_syn[self.mask]).transpose()
        self.f_res_norm__zl = ((self.K.f_obs[self.mask] - self.K.f_syn[self.mask]) / self.K.fobs_norm).transpose()

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
        self.initVars()

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
        self.initVars()

    def rebuildSpectra(self, tomo, eigVec, mean, ne):
        I_rec = np.dot(tomo[:, :ne], eigVec[:, :ne].transpose())
        f_rec = I_rec + mean

        return I_rec, f_rec
