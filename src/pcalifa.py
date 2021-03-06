'''
Created on 17/04/2013

@author: lacerda
'''
import matplotlib
# matplotlib.use('agg')
import numpy as np
import atpy
import pyfits
from pycasso.fitsdatacube import fitsQ3DataCube
from pystarlight import io
from scipy import linalg
from scipy import stats as st
from matplotlib import cm
from matplotlib import gridspec
from matplotlib import pyplot as plt
from matplotlib.ticker import MaxNLocator

__all__ = [ 'PCAlifa' ]

fitsDirDefault = '/home/lacerda/CALIFA/gal_fits'
quantilFlagDefault = 0.9

plt.rcParams.update({'font.family' : 'serif',
                     'text.usetex' : False,
                     })

def PCA(arr, num, axis = -1, arrMean = False, sort = True):
    if not arrMean or not arrMean.any():
        arrMean = arr.mean(axis = axis)

    diff = arr - arrMean
    covMat = np.dot(diff.T, diff) / (num - 1.)
    w, e = linalg.eigh(covMat)
    wS = w
    eS = e

    if sort:
        S = np.argsort(w)[::-1]
        wS = w[S]
        eS = e[:, S]

    return diff, arrMean, covMat, wS, eS

class PCAlifa:
    # Remove califaID and fitsDir... let only FITS load by fitsFile
    def __init__(self, fitsFile = False, quantilQFlag = 1., lc = []):
        self._initVars()
        self.K = None

        if fitsFile:
            self.fitsFile = fitsFile
            self.readCALIFACube()

            # XXX search for califaID from fitsName
            self.califaID = self.K.califaID

            if len(lc):
                self.setLambdaConstrains(lc)

            if quantilQFlag < 1.:
                self.setQFlag(quantilQFlag)

            self.CALIFACubeInfo()

        self.setImgSuffix('pdf')

    def readCALIFACube(self, fitsFile = False):
        self._initVars()

        if fitsFile:
            self.fitsFile = fitsFile

        self.K = fitsQ3DataCube(self.fitsFile)

        self.maskEmLines = np.ones_like(self.K.l_obs, dtype = np.bool)
        self.maskQFlag = self.maskEmLines
        self.maskLambdaConstrains = self.maskEmLines

        self._setVars()

    def readPCAlifaFits(self, fitsFile):
        hdulist = pyfits.open(fitsFile)
        self.l_obs = hdulist['L_OBS'].data
        self.I__zl = np.ma.masked_array(hdulist['I__ZL'].data, mask = hdulist['I_MASK__ZL'].data)
        self.ms__l = np.ma.masked_array(hdulist['MS__L'].data, mask = hdulist['MS_MASK__L'].data)
        self.covMat__ll = hdulist['COVMAT__LL'].data
        self.eigVal__k = hdulist['EIGVAL__K'].data
        self.eigVec__lk = np.ma.masked_array(hdulist['EIGVEC__LK'].data, mask = hdulist['EIGVEC_MASK__LK'].data)
        self.tomo__zk = np.ma.masked_array(hdulist['TOMO__ZK'].data, mask = hdulist['TOMO_MASK__ZK'].data)
        self.tomo__kyx = np.ma.masked_array(hdulist['TOMO__KYX'].data, mask = hdulist['TOMO_MASK__KYX'].data)

    def savePCAlifaFits(self, fitsFile = False, overwrite = False):
        hdulist = pyfits.HDUList()
        hdulist.append(pyfits.ImageHDU(data = self.l_obs, name = 'l_obs'))
        hdulist.append(pyfits.ImageHDU(data = self.I__zl.data, name = 'I__zl'))
        hdulist.append(pyfits.ImageHDU(data = self.I__zl.mask.astype(int), name = 'I_mask__zl'))
        hdulist.append(pyfits.ImageHDU(data = self.ms__l.data, name = 'ms__l'))
        hdulist.append(pyfits.ImageHDU(data = self.ms__l.mask.astype(int), name = 'ms_mask__l'))
        hdulist.append(pyfits.ImageHDU(data = self.covMat__ll, name = 'covMat__ll'))
        hdulist.append(pyfits.ImageHDU(data = self.eigVal__k, name = 'eigVal__k'))
        hdulist.append(pyfits.ImageHDU(data = self.eigVec__lk.data, name = 'eigVec__lk'))
        hdulist.append(pyfits.ImageHDU(data = self.eigVec__lk.mask.astype(int), name = 'eigVec_mask__lk'))
        hdulist.append(pyfits.ImageHDU(data = self.tomo__zk.data, name = 'tomo__zk'))
        hdulist.append(pyfits.ImageHDU(data = self.tomo__zk.mask.astype(int), name = 'tomo_mask__zk'))
        hdulist.append(pyfits.ImageHDU(data = self.tomo__kyx.data, name = 'tomo__kyx'))
        hdulist.append(pyfits.ImageHDU(data = self.tomo__kyx.mask.astype(int), name = 'tomo_mask__kyx'))

        if not fitsFile:
            fitsFile = self.fitsFile.split('/')[-1].replace('.fits', '_pcalifa.fits')

        print 'Writing to %s...' % fitsFile
        hdulist.writeto(fitsFile, clobber = overwrite)

    def setImgSuffix(self, suffix):
        self.imgSuffix = suffix

    def CALIFACubeInfo(self):
        K = self.K

        print 'Galaxy Name: %s' % K.galaxyName
        print 'CALIFA ID: %s' % K.califaID
        print 'N zones: %d' % K.N_zone
        print 'N X: %d' % K.N_x
        print 'N Y: %d' % K.N_y

        total = K.zoneArea_pix.size
        onepix = np.where(K.zoneArea_pix == 1)[0].size
        morethanten = np.where(K.zoneArea_pix > 10)[0].size
        onepixperc = onepix / np.double(total)
        morethantenperc = morethanten / np.double(total)

        print 'zoneArea pix: %d' % total
        print '1 zones - 1 pix: %d (ratio: %f)' % (onepix, onepixperc)
        print '1 zone > 10 pix: %d (ratio: %f)' % (morethanten, morethantenperc)
        print 'N l_obs: %d, l_ini: %d, l_fin: %d' % (K.l_obs.size, K.l_ini, K.l_fin)
        print 'PCALifa:'
        print 'N l_obs: %d, l_ini: %d, l_fin: %d' % (self.l_obs.size, self.l_obs[0], K.l_obs[-1])
        print 'maskLambdaConstrains'
        print 'N l_obs: %d' % K.l_obs[self.maskLambdaConstrains].size
        print 'maskQFlag'
        print 'N l_obs: %d' % K.l_obs[self.maskQFlag].size
        print 'maskQFlag + lambda constrains'
        print 'N l_obs: %d' % K.l_obs[self.maskQFlag & self.maskLambdaConstrains].size

    def PCA(self, arr, num, axis = -1, arrMean = False, sort = True):
        diff, arrMean, covMat, wS, eS = PCA(arr, num, axis, arrMean, sort)
        arrMeanMasked = np.ma.masked_all((self.K.Nl_obs))
        arrMeanMasked[self.mask] = arrMean

        diffMasked = np.ma.masked_all((self.K.N_zone, self.K.Nl_obs))
        diffMasked[self.mask2d.T] = diff.ravel()

        mask2d = np.repeat(self.mask[..., np.newaxis], len(wS), axis = 1)
        esMasked = np.ma.masked_all((self.K.Nl_obs, len(wS)))
        esMasked[mask2d] = eS.ravel()

        return diffMasked, arrMeanMasked, covMat, wS, esMasked

    def PCA_obs(self):
        self.I__zl, self.ms__l, self.covMat__ll, self.eigVal__k, self.eigVec__lk = self.PCA(self.f_obs__zl, self.K.N_zone, 0)

    def PCA_obs_norm(self):
        self.I__zl, self.ms__l, self.covMat__ll, self.eigVal__k, self.eigVec__lk = self.PCA(self.f_obs_norm__zl, self.K.N_zone, 0)

    def PCA_syn(self):
        self.I__zl, self.ms__l, self.covMat__ll, self.eigVal__k, self.eigVec__lk = self.PCA(self.f_syn__zl, self.K.N_zone, 0)

    def PCA_syn_norm(self):
        self.I__zl, self.ms__l, self.covMat__ll, self.eigVal__k, self.eigVec__lk = self.PCA(self.f_syn_norm__zl, self.K.N_zone, 0)

    def PCA_res(self):
        self.I__zl, self.ms__l, self.covMat__ll, self.eigVal__k, self.eigVec__lk = self.PCA(self.f_res__zl, self.K.N_zone, 0)

    def PCA_res_norm(self):
        self.I__zl, self.ms__l, self.covMat__ll, self.eigVal__k, self.eigVec__lk = self.PCA(self.f_res_norm__zl, self.K.N_zone, 0)

    def tomograms(self):
        self.tomo__zk, self.tomo__kyx = self.tomogram(self.I__zl, self.eigVec__lk)

    def tomogram(self, I, eigVec):
        t__zk = np.ma.dot(I, eigVec)
        t__kyx = False

        if self.K:
            t__kyx = self.K.zoneToYX(t__zk.T, extensive = False)

        return t__zk, t__kyx

    def _setVars(self):
        mask = self.maskEmLines & self.maskLambdaConstrains & self.maskQFlag
        synMask = mask

        self.l_obs = self.K.l_obs
        self.l_syn = self.K.l_obs
        self.f_obs__zl = self.K.f_obs[mask].transpose()
        self.f_obs_norm__zl = (self.K.f_obs[mask] / self.K.fobs_norm).transpose()
        self.f_syn__zl = self.K.f_syn[synMask].transpose()
        self.f_syn_norm__zl = (self.K.f_syn[synMask] / self.K.fobs_norm).transpose()
        self.f_res__zl = (self.K.f_obs[mask] - self.K.f_syn[mask]).transpose()
        self.f_res_norm__zl = ((self.K.f_obs[mask] - self.K.f_syn[mask]) / self.K.fobs_norm).transpose()

        self.mask = mask
        self.mask2d = np.repeat(self.mask[..., np.newaxis], self.K.N_zone, axis = 1)

    def _initVars(self):
        self.maskEmLines = []
        self.maskQFlag = []
        self.maskLambdaConstrains = []

        self.l_obs = False
        self.f_obs__zl = False
        self.f_obs_norm__zl = False
        self.f_syn__zl = False
        self.f_syn_norm__zl = False
        self.f_res__zl = False
        self.f_res_norm__zl = False

        self.I__zl = False
        self.ms__l = False
        self.covMat__ll = False
        self.eigVal__k = False
        self.eigVec__lk = False

    def setLambdaConstrains(self, lc):
        lc = np.array(lc)
        self.lc = lc

        s = np.argsort(lc)
        ldown = lc[s][0]
        lup = lc[s][1]

        self.maskLambdaConstrains = (self.K.l_obs >= ldown) & (self.K.l_obs <= lup)

        self._setVars()

    def unsetLambdaConstrains(self):
        self.lc = False
        self.maskLambdaConstrains = np.ones_like(self.K.l_obs, dtype = np.bool)
        self._setVars()

    def setQFlag(self, quantil):
        K = self.K
        self.quantilQFlag = quantil

        self.histo = (K.f_flag > 0).astype(int).sum(axis = 1) / K.N_zone
        self.maskQFlag = (self.histo < (1.0 - quantil))
        self._setVars()

    def unsetQFlag(self):
        self.quantilQFlag = 1.
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

        print 'maskEmLines'
        print 'N l_obs: %d' % self.K.l_obs[self.maskEmLines].size
        print 'maskEmLines + lambda constrains'
        print 'N l_obs: %d' % self.K.l_obs[self.maskEmLines & self.maskLambdaConstrains].size
        print 'maskEmLines + maskQFlag + lambda constrains'
        print 'N l_obs: %d' % self.K.l_obs[self.maskEmLines & self.maskQFlag & self.maskLambdaConstrains].size

    def unsetStarlightMaskFile(self):
        del self.tStarlight
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
                fxxontsize = fontsize + 1, transform = ax.transAxes,
                horizontalalignment = 'left',
                verticalalignment = 'center',
                multialignment = 'left')
        ax.legend(prop = {'size' : fontsize})
        ax.grid()

    def correlationAxisPlot(self, x, y, ax):
        rhoSpearman, pvalSpearman = st.spearmanr(x, y)
        txt = 's: %.2f' % rhoSpearman

        ax.scatter(x, y, marker = 'o', s = 0.1)

        textbox = dict(boxstyle = 'round', facecolor = 'wheat', alpha = 0.)

        ax.text(0.76, 0.15, txt,
                fontsize = 8,
                transform = ax.transAxes,
                verticalalignment = 'top',
                bbox = textbox)

        plt.setp(ax.get_yticklabels(), visible = False)

    def tomoPlot(self, ti, npref):
        t = self.tomo__kyx
        eigvec = self.eigVec__lk
        eigval = self.eigVal__k
        ms = self.ms__l

        tomogram = t[ti, :, :]
        x = self.l_obs
        y = eigvec[:, ti]
        y2 = ms

        fig = plt.figure(figsize = (15, 5))
        # fig = plt.figure()
        gs = gridspec.GridSpec(1, 2, width_ratios = [4, 7])
        ax1 = plt.subplot(gs[0])
        ax2 = plt.subplot(gs[1])

        ax1.set_title(r'tomogram $%02i$' % ti)
        im = ax1.imshow(tomogram, origin = 'lower', interpolation = 'nearest', aspect = 'auto')
        fig.colorbar(ax = ax1, mappable = im, use_gridspec = True, cmap = cm.jet)

        eigval_norm = 100. * eigval[ti] / eigval.sum()
        ax2.set_title(r'Eigenvalue $%s$ ($%.2f\ \%%$)' % (self.nToStrSciNot(eigval[ti]), eigval_norm))
        ax2.plot(x, y)
        ax2.xaxis.set_major_locator(MaxNLocator(20))
        ax2.set_ylabel(r'$PC %02i$' % ti)
        plt.setp(ax2.get_xticklabels(), rotation = 45)
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

        fig.set_tight_layout(True)

        if npref:
            fig.savefig('%stomo_%02i.%s' % (npref, ti, self.imgSuffix))
        else:
            fig.show()

    def screeTestPlot(self, maxInd, title, npref):
        # set_eps_output_1()
        eigval = self.eigVal__k

        f = plt.figure()
        f.set_size_inches(19.2, 10.8)
        eigval_norm = 100. * eigval / eigval.sum()
        plt.plot(eigval_norm[:maxInd], linestyle = '-', marker = '*')
        plt.ylim([0, eigval_norm[1] * 1.1])
        plt.xticks(range(maxInd))
        plt.title(r'%s Scree test' % title)
        plt.xlabel(r'Principal component')
        plt.ylabel(r'Normalized variance [$\%%$]')
        plt.grid()

        if npref:
            plt.savefig('%sscree.%s' % (npref, self.imgSuffix))
        else:
            plt.show()

    def nToStrSciNot(self, n):
        e = np.floor(np.log10(np.abs(n)))
        m = n / 10.**(e)

        if np.int(e) == 0:
            nStr = r'%.2f' % n
        else:
            nStr = r'%.2f \times 10^{%d}' % (m, e)

        return nStr

    def sanityCheck(self, x, y, npref):
        # set_eps_output_1()
        f = plt.figure()
        f.set_size_inches(19.2, 10.8)
        plt.plot(x, y)
        plt.plot(x, y, 'r.')
        plt.ylabel(r'Obs. flux $[erg/s/cm^2/\AA]$')
        plt.xlabel(r'Wavelength $[\AA]$')
        plt.title(r'Mask check (zone 0)')
        plt.grid()

        if npref:
            f.savefig('%ssanity_check.%s' % (npref, self.imgSuffix))
        else:
            f.show()

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
