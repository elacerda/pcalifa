'''
Created on 17/04/2013

@author: lacerda
'''
from pycasso.fitsdatacube import fitsQ3DataCube
import pystarlight.io
import numpy as np
import atpy
from scipy import linalg

fitsDir = '/home/lacerda/CALIFA'

def loadOneCube(califaID):
    fitsFile = '%s/%s/%s_synthesis_eBR_v20_q027.d13c512.ps3b.k1.mC.CCM.Bgsd01.v01.fits' % (fitsDir, califaID, califaID)
    return fitsQ3DataCube(fitsFile)

def PCA(arr, num, axis = -1, arrMean = None):
    if arrMean == None:
        arrMean = arr.mean(axis = axis)

    diff = arr - arrMean
    covMat = np.dot(diff.T, diff) / (num - 1.)
    w, e = linalg.eig(covMat)

    return diff, covMat, arrMean, np.real(w), np.real(e)

def tomogram(K, I, eigvec, extensive = True):
    t__zk = np.dot(I, eigvec)
    #t__kyx = K.zoneToYX(t__zk.T, extensive = False)
    t__kyx = K.zoneToYX(t__zk.T, extensive = extensive)

    return t__zk, t__kyx

def removeStarlightMask(maskfile, l_obs):
    t = atpy.Table(maskfile = maskfile, type = 'starlight_mask')

    mask = (l_obs > t[0]['l_up'])

    for i in range(1, len(t)):
        if (t[i]['weight'] == 0.0):
            mask = mask & ((l_obs < t[i]['l_low']) | (l_obs > t[i]['l_up']))

    return mask

