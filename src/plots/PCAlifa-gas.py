'''
Created on 01/08/2013

@author: lacerda
'''
import argparse as ap
import numpy as np
import PCAlifa as PCA
from astropy.io import fits

gasFileSuffix = '_synthesis_eBR_v20_q036.d13c512.ps03.k2.mC.CCM.Bgsd61.EML.MC100.fits'

def parser_args():
    parser = ap.ArgumentParser(description = 'PCAlifa - correlations')
    parser.add_argument('--califaID', '-c',
                        help = 'Califa ID (ex: K0277)',
                        type = str,
                        default = 'K0277')
    parser.add_argument('--fitsDir', '-f',
                        help = 'Califa FITS directory',
                        metavar = 'DIR',
                        type = str,
                        default = '/home/lacerda/CALIFA')
    parser.add_argument('--gasFitsDir', '-g',
                        help = 'Califa Gas FITS directory',
                        metavar = 'DIR',
                        type = str,
                        default = '/home/lacerda/CALIFA/rgb-gas')
    parser.add_argument('--tmax', '-t',
                        help = 'Number max of eigenvectors to plot',
                        metavar = 'INT',
                        type = int,
                        default = 20)
    parser.add_argument('--eigv', '-e',
                        help = 'number of eigenvectors used to rebuild the signal',
                        type = int,
                        nargs = '+',
                        default = [1, 2, 3, 4, 5, 6, 10, 20])
    parser.add_argument('--zones', '-z',
                        help = 'The program will rebuild the spectrum of ZONES zones + last zone.',
                        type = int,
                        nargs = '+',
                        default = [0, 10, 20, 100, 200])
    parser.add_argument('--rebuildSpec', '-R',
                        help = 'Plot the rebuilds --zones spectra.',
                        action = 'store_true',
                        default = False)
    parser.add_argument('--tomograms', '-T',
                        help = 'Plots tomograms.',
                        action = 'store_true',
                        default = False)
    parser.add_argument('--correlat', '-C',
                        help = 'Plots correlations.',
                        action = 'store_true',
                        default = False)
    parser.add_argument('--numcorrepc',
                        help = 'Number of PCs to plot correlations.',
                        type = int,
                        default = 5)
    parser.add_argument('--outputfname',
                        help = 'If provided, this adds OUTPUTFNAME to filenames.',
                        type = str,
                        default = False)

    return parser.parse_args()

def parseArrArgs(args, lastZone):
    args.eigvArr = np.array(args.eigv)
    args.zones.append(lastZone)
    args.zonesArr = np.array(args.zones)

if __name__ == '__main__':
    args = parser_args()

    P = PCA.PCAlifa(args.califaID, args.fitsDir, False)

    parseArrArgs(args, P.K.N_zone - 1)
    # This call of PCAlifa is only to use P.K.zoneToYX(), P.K.N_zone and PCA methods;
    hdu = fits.open('%s/%s%s' % (args.gasFitsDir, args.califaID, gasFileSuffix))

    l = []

    n_hdu = len(hdu)
    n_data = n_hdu - 2

    flux__lz = np.zeros((n_data, P.K.N_zone))
    fwhm__lz = np.zeros((n_data, P.K.N_zone))
    ew__lz = np.zeros((n_data, P.K.N_zone))

    # Building the data arrays
    for i in range(n_data):
        j = i + 2
        flux__lz[i] = np.array(object = hdu[j].data['flux'], dtype = hdu[j].data.dtype['flux'])
        fwhm__lz[i] = np.array(object = hdu[j].data['fwhm'], dtype = hdu[j].data.dtype['fwhm'])
        ew__lz[i] = np.array(object = hdu[j].data['EW'], dtype = hdu[j].data.dtype['EW'])

        l.append(hdu[j].header['EXTNAME'])

    l = np.array(l, dtype = np.int)

    I_flux__zl, ms_flux__l, covMat_flux__ll, eigVal_flux__k, eigVec_flux__lk = P.PCA(flux__lz.T, P.K.N_zone, 0)
    I_fwhm__zl, ms_fwhm__l, covMat_fwhm__ll, eigVal_fwhm__k, eigVec_fwhm__lk = P.PCA(fwhm__lz.T, P.K.N_zone, 0)
    I_ew__zl, ms_ew__l, covMat_ew__ll, eigVal_ew__k, eigVec_ew__lk = P.PCA(ew__lz.T, P.K.N_zone, 0)

    tomo_flux__zk, tomo_flux__kyx = P.tomogram(I_flux__zl, eigVec_flux__lk)
    tomo_fwhm__zk, tomo_fwhm__kyx = P.tomogram(I_fwhm__zl, eigVec_fwhm__lk)
    tomo_ew__zk, tomo_ew__kyx = P.tomogram(I_ew__zl, eigVec_ew__lk)

    npref_flux = '%s-flux_' % P.K.califaID
    npref_fwhm = '%s-fwhm_' % P.K.califaID
    npref_ew = '%s-ew_' % P.K.califaID

    if args.outputfname:
        npref_flux = '%s%s_' % (npref_flux, args.outputfname)
        npref_fwhm = '%s%s_' % (npref_fwhm, args.outputfname)
        npref_ew = '%s%s_' % (npref_ew, args.outputfname)

    if args.tomograms:
        P.screeTestPlot(eigVal_flux__k, args.tmax, npref_flux, '%s flux' % P.K.califaID)
        P.screeTestPlot(eigVal_fwhm__k, args.tmax, npref_fwhm, '%s fwhm' % P.K.califaID)
        P.screeTestPlot(eigVal_ew__k, args.tmax, npref_ew, '%s ew' % P.K.califaID)

        for ti in range(args.tmax):
            P.tomoPlot(tomo_flux__kyx, l, eigVec_flux__lk, eigVal_flux__k, ms_flux__l, ti, npref_flux)
            P.tomoPlot(tomo_fwhm__kyx, l, eigVec_fwhm__lk, eigVal_fwhm__k, ms_fwhm__l, ti, npref_fwhm)
            P.tomoPlot(tomo_ew__kyx, l, eigVec_ew__lk, eigVal_ew__k, ms_ew__l, ti, npref_ew)
