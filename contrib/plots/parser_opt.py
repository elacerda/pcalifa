#!/usr/bin/python
# -*- coding: utf-8 -*-
import sys
import argparse as ap
import numpy as np

default = { 
    'rSEL' : '/home/lacerda/CALIFA/Mask.mC',
    'rFL' : 0.95,
    'lc' : [3800, 6850],
    'outputImgSuffix' : 'png',
    'outputDir' : '.',
    'firstTomo' : -1,
    'nTomo' : 5,
    'zone' : 0,
    'eigv' : [1, 2, 3, 4, 5, 6, 10, 20],
}

def parser_args():
    parser = ap.ArgumentParser(description = '%s' % sys.argv[0])
    parser.add_argument('--fitsfile', '-f',
                        help = 'The file must be named KXXXX*.fits',
                        metavar = 'PyCASSO FITS FILE',
                        type = str,
                        default = None)
    parser.add_argument('--lc', '-l',
                        help = 'Lambda constrains',
                        metavar = 'LAMBDA',
                        type = int,
                        nargs = 2,
                        default = default['lc'])
    parser.add_argument('--firstTomo', '-F',
                        help = 'First tomogram to print',
                        metavar = 'NUM',
                        type = int,
                        default = default['firstTomo'])
    parser.add_argument('--nTomo', '-N',
                        help = 'Number of tomograms to print',
                        metavar = 'NUM',
                        type = int,
                        default = default['nTomo'])
    parser.add_argument('--zone', '-z',
                        help = 'Zone number to print or rebuild',
                        metavar = 'NUM',
                        type = int,
                        default = default['zone'])
    parser.add_argument('--galaxyImgFile', '-g',
                        help = 'The image of the galaxy',
                        metavar = 'FILE',
                        type = str,
                        default = None)
    parser.add_argument('--outputImgSuffix', '-o',
                        help = 'Suffix of image file. Sometimes denote the image type. (Ex.: image.png)',
                        type = str,
                        default = default['outputImgSuffix'])
    parser.add_argument('--outputDir', '-d',
                        help = 'Image output directory',
                        metavar = 'DIR',
                        type = str,
                        default = default['outputDir'])
    parser.add_argument('--rSEL', '-S',
                        help = 'Remove Starlight Emission Lines ',
                        metavar = 'MASK FILENAME',
                        type = str,
                        default = default['rSEL'])
    parser.add_argument('--rFL', '-Q',
                        help = 'Remove Flagged Lamdas',
                        metavar = 'QUANTIL',
                        type = float,
                        default = default['rFL'])
    parser.add_argument('--eigv', '-e',
                        help = 'number of eigenvectors used to rebuild the signal',
                        type = int,
                        nargs = '+',
                        default = default['eigv'])

    return parser.parse_args()

def parseArrArgs(args):
    args.eigvArr = np.array(args.eigv)
 
