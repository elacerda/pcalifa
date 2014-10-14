#!/bin/bash
SRCDIR=/Users/lacerda/dev/astro/PCA/src
WORKDIR=/Users/lacerda/CALIFA
STARLIGHTMASKFILE=${WORKDIR}/Mask.mC
#FITSVERS="v20_q036.d13c"
#FITSVERS="v20_q043.d14a"
FITSVERS="px1_q043.d14a"
FITSDIR=${WORKDIR}/gal_fits/${FITSVERS}

GAL_ET="K0864 K0119"
GAL_SP="K0073 K0518 K0008 K0277"
GAL_SP2="K0073 K0518 K0008"
GAL_ME="K0213 K0802"
GAL="$GAL_SP $GAL_ET $GAL_ME"
GAL="K0577"

#FITSUFFIX="_synthesis_eBR_${FITSVERS}512.ps03.k2.mC.CCM.Bgsd61.fits"
FITSUFFIX="_synthesis_eBR_${FITSVERS}512.ps03.k1.mE.CCM.Bgsd6e.fits"

OUTPUTDIR=.
OUTPUTIMGFMT=pdf

for g in $GAL
do
    FITSFILE=${FITSDIR}/${g}${FITSUFFIX}
    GALIMGFILE=${WORKDIR}/images/${g}.jpg

#    time ./tudo.py -f $FITSFILE -d $OUTPUTDIR -o $OUTPUTIMGFMT -g $GALIMGFILE 
    time ./apres_gal.py -f $FITSFILE -d $OUTPUTDIR -o $OUTPUTIMGFMT -g $GALIMGFILE  
    time ./correPCvsPC.py -f $FITSFILE -d $OUTPUTDIR -o $OUTPUTIMGFMT 
    time ./correPCvsPhys.py -f $FITSFILE -d $OUTPUTDIR -o $OUTPUTIMGFMT 
    time ./tomo_obs_norm.py -f $FITSFILE -d $OUTPUTDIR -o $OUTPUTIMGFMT 
    time ./tomo_syn_norm.py -f $FITSFILE -d $OUTPUTDIR -o $OUTPUTIMGFMT  
#    time ./scree-test.py -f $FITSFILE -d $OUTPUTDIR -o $OUTPUTIMGFMT
#
done
