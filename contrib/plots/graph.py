#!/usr/bin/python
import numpy as np
from matplotlib import pyplot as plt
from matplotlib import cm
from atpy import Table
from parser_opt import *
from pycasso import fitsQ3DataCube

args = parser_args()

# Carregar arquivo FITS com os dados.
K = fitsQ3DataCube(args.fitsfile)

# Converter zonas para imagem.
at_image = K.zoneToYX(K.at_flux__z, extensive = False)

## Desenhar o mapa.
f = plt.figure()
plt.imshow(at_image, origin = 'lower', interpolation = 'nearest')
cb = plt.colorbar()
cb.set_label(r'$\langle \log$ t $\rangle_L$ [anos]')
plt.title(r'%s - %s' % (K.galaxyName, K.califaID))
plt.xlabel(r'Pixels')
f.savefig('%s/%s-at_flux_zone.%s' % (args.outputDir, K.califaID, args.outputImgSuffix))

# Calcular o perfil radial.
bins = np.arange(0, 26, 1)
bin_center = (bins[1:] + bins[:-1]) / 2.0
at_rad = K.radialProfile(at_image, bins, rad_scale = 1.0)
LobsnSD__r = K.radialProfile(K.LobnSD__yx, bin_r=bins, rad_scale=1, mode='sum')
at_flux_LobsnSD__r = K.radialProfile(K.at_flux__yx * K.LobnSD__yx, 
                                     bin_r=bins, rad_scale=1, mode='sum')
at_flux__r = at_flux_LobsnSD__r / LobsnSD__r 

# Desenhar o perfil.
f = plt.figure()
f.set_size_inches(7, 5)
plt.title(r'%s - %s' % (K.galaxyName, K.califaID))
plt.xlabel('radius [arcsec]')
plt.ylabel(r'$\langle \log$ t $\rangle_L$ [anos]')
plt.plot(bin_center, at_rad)
f.savefig('%s/%s-at_flux_radprof.%s' % (args.outputDir, K.califaID, args.outputImgSuffix))
