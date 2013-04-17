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
    def __init__(self, califaID, fitsDir = fitsDirDefault, flagLinesQuantil = 0.9, remFlaggedLambdas = True):
        self.califaID = califaID
        self.flagLinesQuantil = flagLinesQuantil
        self.fitsDir = fitsDir
        self.fitsFile = '%s/%s/%s_synthesis_eBR_v20_q027.d13c512.ps3b.k1.mC.CCM.Bgsd01.v01.fits' % (fitsDir, califaID, califaID)
        self.K = fitsQ3DataCube(self.fitsFile)
        self.histo = None
        self.mask = np.ones_like(self.K.l_obs, dtype = np.bool)
        self.maskEmLines = self.mask
        self.tStarlight = None
        self.remFlaggedLambdas = remFlaggedLambdas

        self.initVars()

        if (self.remFlaggedLambdas == True):
            self.removeFlaggedLambda()

    def PCA(self, arr, num, axis = -1, arrMean = None):
        if arrMean == None:
            arrMean = arr.mean(axis = axis)

        diff = arr - arrMean
        covMat = np.dot(diff.T, diff) / (num - 1.)
        w, e = linalg.eig(covMat)

        S = np.argsort(w)[::-1]
        wS = w[S]
        eS = e[:, S]

        return diff, arrMean, np.real(wS), np.real(eS)

    def tomogram(self, I, eigvec, extensive = True):
        t__zk = np.dot(I, eigvec)
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

    def removeFlaggedLambda(self):
        self.histo = self.K.f_flag.sum(axis = 1) / self.K.N_zone
        self.maskHisto = (self.histo < self.flagLinesQuantil)
        self.mask = self.maskHisto & (self.K.l_obs > 3800) & (self.K.l_obs < 6850)
        self.initVars()

    def removeStarlightEmLines(self, maskfile):
        t = atpy.Table(maskfile = maskfile, type = 'starlight_mask')
        mask = (self.l_obs > t[0]['l_up'])

        for i in range(1, len(t)):
            if ([i]['weight'] == 0.0):
                mask = mask & ((self.l_obs < t[i]['l_low']) | (self.l_obs > t[i]['l_up']))

        self.tStarlight = t
        self.mask = self.mask & mask
        self.initVars(self)
