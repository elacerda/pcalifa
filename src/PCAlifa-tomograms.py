'''
Created on 19/10/2012

@author: lacerda
'''
import matplotlib
matplotlib.use('agg')

import sys
import numpy as np
import matplotlib.gridspec as gridspec
from matplotlib.ticker import MaxNLocator
from matplotlib import pyplot as plt
import PCAlifa as PCA
import argparse as ap


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
    parser.add_argument('--rSEL', '-S',
                        help = 'Remove Starlight Emission Lines ',
                        metavar = 'MASK FILENAME',
                        type = str,
                        default = False)
    parser.add_argument('--rFL', '-Q',
                        help = 'Remove Flagged Lamdas',
                        metavar = 'QUANTIL',
                        type = float,
                        default = 0.9)
    parser.add_argument('--tmax', '-t',
                        help = 'Number max of eigenvectors to plot',
                        metavar = 'INT',
                        type = int,
                        default = 20)

    return parser.parse_args()

if __name__ == '__main__':
    args = parser_args()

    P = PCA.PCAlifa(califaID = args.califaID,
                    fitsDir = args.fitsDir,
                    flagLinesQuantil = args.rFL)

    if args.rSEL:
        P.removeStarlightEmLines(args.rSEL)

    P.runPCA()

#########################################################################
################################ f_obs ##################################
#########################################################################

    npref = '%s_f_obs_' % P.K.califaID

    P.screeTestPlot(P.eigVal_obs__k, args.tmax, npref)

    for tn in range(args.tmax):
        P.tomoPlot(tn, P.l_obs, P.tomo_obs__kyx, P.eigVec_obs__lk, P.eigVal_obs__k, npref)

########################################################################
########################### f_obs_norm #################################
########################################################################

    npref = '%s_f_obs_norm_' % P.K.califaID

    P.screeTestPlot(P.eigVal_obs_norm__k, args.tmax, npref)

    for tn in range(args.tmax):
        P.tomoPlot(tn, P.l_obs, P.tomo_obs_norm__kyx, P.eigVec_obs_norm__lk, P.eigVal_obs_norm__k, npref)

########################################################################
############################# f_syn ####################################
########################################################################

    npref = '%s_f_syn_' % P.K.califaID

    P.screeTestPlot(P.eigVal_syn__k, args.tmax, npref)

    for tn in range(args.tmax):
        P.tomoPlot(tn, P.l_obs, P.tomo_syn__kyx, P.eigVec_syn__lk, P.eigVal_syn__k, npref)

########################################################################
########################### f_syn_norm #################################
########################################################################

    npref = '%s_f_syn_norm_' % P.K.califaID

    P.screeTestPlot(P.eigVal_syn_norm__k, args.tmax, npref)

    for tn in range(args.tmax):
        P.tomoPlot(tn, P.l_obs, P.tomo_syn_norm__kyx, P.eigVec_syn_norm__lk, P.eigVal_syn_norm__k, npref)

########################################################################
############################### f_res ##################################
########################################################################

    npref = '%s_f_res_' % P.K.califaID

    P.screeTestPlot(P.eigVal_res__k, args.tmax, npref)

    for tn in range(args.tmax):
        P.tomoPlot(tn, P.l_obs, P.tomo_res__kyx, P.eigVec_res__lk, P.eigVal_res__k, npref)

########################################################################
############################# f_res_norm ###############################
########################################################################

    npref = '%s_f_res_norm_' % P.K.califaID

    P.screeTestPlot(P.eigVal_res_norm__k, args.tmax, npref)

    for tn in range(args.tmax):
        P.tomoPlot(tn, P.l_obs, P.tomo_res_norm__kyx, P.eigVec_res_norm__lk, P.eigVal_res_norm__k, npref)
