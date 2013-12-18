'''
Created on 01/08/2013

@author: lacerda
'''
import argparse as ap
import numpy as np
import PCAlifa as PCA
from astropy.io import fits
from matplotlib import pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MaxNLocator

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
                        default = '/home/lacerda/CALIFA/gal_fits/')
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

def tomoPlot(t, l, eigvec, eigval, ms, ti, npref):
    P = PCA.PCAlifa()
    tomogram = t[ti, :, :]
    x = l
    y = eigvec[:, ti]
    y2 = ms

    N = len(l)
    ind = np.arange(N)
    width = 0.35

    fig = plt.figure(figsize = (15, 5))
    gs = gridspec.GridSpec(1, 2, width_ratios = [4, 7])
    ax1 = plt.subplot(gs[0])
    ax2 = plt.subplot(gs[1])

    ax1.set_title(r'tomogram $%02i$' % ti)
    im = ax1.imshow(tomogram, origin = 'lower', interpolation = 'nearest', aspect = 'auto')
    fig.colorbar(ax = ax1, mappable = im, use_gridspec = True)

    eigval_norm = 100. * eigval[ti] / eigval.sum()
    ax2.set_title(r'Eigenvalue $%s$ ($%.2f\ \%%$)' % (P.nToStrSciNot(eigval[ti]), eigval_norm))
    ax2.bar(ind, y)
    ax2.set_xticks(ind + width)
    ax2.set_xticklabels(l)
    ax2.set_ylabel(r'$PC %02i$' % ti)

#    plt.setp(xticsNames, rotation = 45)
#    ax2.xaxis.set_major_locator(MaxNLocator(len(l)))
#    ax2.grid()

#    ax3 = ax2.twinx()
#    ax3.plot(np.arange(len(y2)), y2, 'o', color = '0.55')
#    ax3.vlines(np.arange(len(y2)), [0], y2)
#    ax3.set_ylabel(r'Mean spetrum')

    plt.setp(ax2, xticklabels = l)
    plt.tight_layout()
    fig.savefig('%stomo_%02i.png' % (npref, ti))
    plt.close()

if __name__ == '__main__':
    args = parser_args()

    # This call of PCAlifa is only to use P.K.zoneToYX(), P.K.N_zone and PCA methods;
    P = PCA.PCAlifa(args.califaID, args.fitsDir, False)

    parseArrArgs(args, P.K.N_zone - 1)

    hdu = fits.open('%s/%s%s' % (args.gasFitsDir, args.califaID, gasFileSuffix))

    l = []

    n_hdu = len(hdu)

    corelines = ['4861', '5007', '6583']
    othlines = ['3727', '6300', '6717', '6731']
    totlines = corelines + othlines

    n_data = len(totlines)

    '''
    Como n_data eh a dimensao de componentes principais e esse numero eh pequeno
    o argumento --tmax eh redefinido para n_data.
    '''
    args.tmax = n_data

    flux__lz = np.zeros((n_data, P.K.N_zone))
    fwhm__lz = np.zeros((n_data, P.K.N_zone))
    ew__lz = np.zeros((n_data, P.K.N_zone))
    eflux__lz = np.zeros((n_data, P.K.N_zone))

    # Building the data arrays
    j = 0
    lines = totlines

    lineindex = {}

    f_Ha = False

    for i in range(n_hdu):
        if not f_Ha and hdu[i].name == '6563':
            f_Ha__z = np.array(object = hdu[i].data['flux'], dtype = hdu[i].data.dtype['flux'])
            ef_Ha__z = np.array(object = hdu[i].data['eflux'], dtype = hdu[i].data.dtype['eflux'])

            f_Ha = True

        for line in lines:
            if hdu[i].name == line:
                flux__lz[j] = np.array(object = hdu[i].data['flux'], dtype = hdu[i].data.dtype['flux'])
                fwhm__lz[j] = np.array(object = hdu[i].data['fwhm'], dtype = hdu[i].data.dtype['fwhm'])
                ew__lz[j] = np.array(object = hdu[i].data['EW'], dtype = hdu[i].data.dtype['EW'])
                eflux__lz[j] = np.array(object = hdu[i].data['eflux'], dtype = hdu[i].data.dtype['eflux'])

                lineindex[line] = j

                l.append(hdu[i].name)
                j = j + 1

    l = np.array(l, dtype = np.int)

    ########### MASKS ###########
    mask = (flux__lz[lineindex['4861']] > 0) & (flux__lz[lineindex['5007']] > 0)
    mask &= (f_Ha__z > 0) & (flux__lz[lineindex['6583']] > 0)
    mask &= (flux__lz[lineindex['3727']] > 0) & (flux__lz[lineindex['6300']] > 0)
    mask &= (flux__lz[lineindex['6717']] > 0) & (flux__lz[lineindex['6731']] > 0)
    mask &= (eflux__lz[lineindex['4861']] > 0)

    maskHa = (f_Ha__z > 0) & (ef_Ha__z > 0)

    SNHa = f_Ha__z / ef_Ha__z
    SNHb = flux__lz[lineindex['4861']] / eflux__lz[lineindex['4861']]
    SNO1 = flux__lz[lineindex['6300']] / eflux__lz[lineindex['6300']]
    SNO3 = flux__lz[lineindex['5007']] / eflux__lz[lineindex['5007']]

    masklines = mask & maskHa

    masklines = mask & maskHa & (SNHb >= 3.0)
    print "Galaxia com %i zonas. %i mascaradas" % (P.K.N_zone, P.K.N_zone - np.where(masklines == True)[0].shape[0])
    masklines__lz = ~(masklines & np.ones(flux__lz.shape, dtype = np.bool))

    ########### MASKED DATA ###########
    flux_m__lz = np.ma.array(flux__lz, mask = masklines__lz, fill_value = np.nan)
    flux_m_normHa__lz = np.ma.array(flux__lz, mask = masklines__lz, fill_value = np.nan)
    fwhm_m__lz = np.ma.array(fwhm__lz, mask = masklines__lz, fill_value = np.nan)
    ew_m__lz = np.ma.array(ew__lz, mask = masklines__lz, fill_value = np.nan)

    flux_m_normHa__lz = flux_m_normHa__lz / f_Ha__z

    ########### PCA ###########
    I_flux__zl, ms_flux__l, covMat_flux__ll, eigVal_flux__k, eigVec_flux__lk = P.PCA(flux_m_normHa__lz.T, P.K.N_zone, 0)
    # I_flux__zl, ms_flux__l, covMat_flux__ll, eigVal_flux__k, eigVec_flux__lk = P.PCA(flux_m__lz.T, P.K.N_zone, 0)
    I_fwhm__zl, ms_fwhm__l, covMat_fwhm__ll, eigVal_fwhm__k, eigVec_fwhm__lk = P.PCA(fwhm_m__lz.T, P.K.N_zone, 0)
    I_ew__zl, ms_ew__l, covMat_ew__ll, eigVal_ew__k, eigVec_ew__lk = P.PCA(ew_m__lz.T, P.K.N_zone, 0)

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
        P.screeTestPlot(eigVal_flux__k, args.tmax, '%s flux' % P.K.califaID, npref_flux)
        P.screeTestPlot(eigVal_fwhm__k, args.tmax, '%s fwhm' % P.K.califaID, npref_fwhm)
        P.screeTestPlot(eigVal_ew__k, args.tmax, '%s ew' % P.K.califaID, npref_ew)

        for ti in range(args.tmax):
            tomoPlot(tomo_flux__kyx, l, eigVec_flux__lk, eigVal_flux__k, ms_flux__l, ti, npref_flux)
            tomoPlot(tomo_fwhm__kyx, l, eigVec_fwhm__lk, eigVal_fwhm__k, ms_fwhm__l, ti, npref_fwhm)
            tomoPlot(tomo_ew__kyx, l, eigVec_ew__lk, eigVal_ew__k, ms_ew__l, ti, npref_ew)

    if args.correlat:
        logO3Hb = np.log10(flux_m__lz[lineindex['5007']] / flux_m__lz[lineindex['4861']])
        logN2Ha = np.log10(flux_m__lz[lineindex['6583']] / f_Ha__z)
        logHaHb = np.log10(f_Ha__z / flux_m__lz[lineindex['4861']])
        logS2S2 = np.log10(flux_m__lz[lineindex['6717']] / flux_m__lz[lineindex['6731']])

        colArr = [
                  logO3Hb,
                  logN2Ha,
                  logO3Hb - logN2Ha,
                  logHaHb,
                  logS2S2,
        ]

        colNames = [
            r'$\log\ F_{[OIII]} / F_{H\beta}$',
            r'$\log\ F_{[NII]} / F_{H\alpha}$',
            r'$\log\ (F_{[OIII]} / F_{H\beta}) / (F_{[NII]} / F_{H\alpha})$',
            r'$\log\ F_{H\alpha} / F_{H\beta}$',
            r'$\log\ F_{6717} / F_{6731}$',
            r'eigenvector',
        ]

        nRows = args.tmax
        nCols = len(colArr) + 1

        ###############################
        ###############################
        ###############################

        f, axArr = plt.subplots(nRows, nCols)
        f.set_size_inches(19.2, 10.8)

        N = len(l)
        ind = np.arange(N)
        width = 0.35

        for i in range(nRows):
            axArr[i, 0].set_ylabel('PC%d' % i)

            for j in range(nCols)[:-1]:
                P.correlationAxisPlot(colArr[j], np.log10(tomo_flux__zk[:, i]), axArr[i, j])

            axArr[i, nCols - 1].bar(ind, eigVec_flux__lk[:, i])
            axArr[i, nCols - 1].set_xticks(ind + width)
            axArr[i, nCols - 1].set_xticklabels(l)
            plt.setp(axArr[i, nCols - 1].get_yticklabels(), visible = False)

        f.subplots_adjust(hspace = 0.0)
        f.subplots_adjust(wspace = 0.05)

        plt.setp([a.get_xticklabels() for a in f.axes], rotation = 45)
        plt.setp([a.get_xticklabels() for a in f.axes[:-nCols]], visible = False)
        plt.setp([a.get_yticklabels() for a in f.axes[::nCols]], visible = True)

        for i in range(nCols):
            axArr[0, i].set_title(colNames[i])

        plt.suptitle(r'Correlations - Flux - %i masked zones')
        f.savefig('%s-corre-%s_flux.png' % (args.outputfname, P.K.califaID))
        plt.close()
