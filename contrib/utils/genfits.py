#!/usr/bin/python
import pcalifa as PCA
from parser_opt import *

args = parser_args()

print 'FITSFILE: %s' % args.fitsfile
# remove flagged lambdas
print 'rFL: %f' % args.rFL
# lambda constrains
print 'lc: %s' % args.lc
# remove STARLIGHT emission lines mask
print 'rSEL: %s' % args.rSEL

P = PCA.PCAlifa(fitsFile = args.fitsfile, quantilQFlag = args.rFL, lc = args.lc)

rSELMask = ''

if args.rSEL:
    rSELMask = '_STMASK'
    P.setStarlightMaskFile(args.rSEL)

P.PCA_obs()
P.tomograms()
outFitsFile = P.fitsFile.split('/')[-1].replace('.fits', '_f_obs_PCA_%.2f%s.fits' % (args.rFL, rSELMask))
P.savePCAlifaFits(outFitsFile)

P.PCA_syn()
P.tomograms()
outFitsFile = P.fitsFile.split('/')[-1].replace('.fits', '_f_syn_PCA_%.2f%s.fits' % (args.rFL, rSELMask))
P.savePCAlifaFits(outFitsFile)

P.PCA_obs_norm()
P.tomograms()
outFitsFile = P.fitsFile.split('/')[-1].replace('.fits', '_f_obs_norm_PCA_%.2f%s.fits' % (args.rFL, rSELMask))
P.savePCAlifaFits(outFitsFile)

P.PCA_syn_norm()
P.tomograms()
outFitsFile = P.fitsFile.split('/')[-1].replace('.fits', '_f_syn_norm_PCA_%.2f%s.fits' % (args.rFL, rSELMask))
P.savePCAlifaFits(outFitsFile)
