'''
Created on 19/10/2012

@author: lacerda
'''

from califa.Q3DataCube import Q3DataCube
import numpy as np
from scipy import linalg

def PCA(arr, axis = 1, arrMean = None):
    if arrMean == None:
        arrMean = arr.mean(axis = axis)

    diff = arr - arrMean.reshape(len(arrMean), 1)
    covMat = 1. / (arr.shape[1] - 1.) * np.dot(diff, diff.T)
    w, e = linalg.eig(covMat)
    return arrMean, w, e

if __name__ == '__main__':
    fitsfile = '/home/lacerda/CALIFA/K0277/K0277_synthesis_eBR_v20_q027.d13c512.ps3b.k1.mC.CCM.Bgsd01.v01.fits'

    c = Q3DataCube(fitsfile)

    ms, eigval, eigvec = PCA(c.f_obs, 1)
