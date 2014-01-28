#!/bin/bash
SRCDIR=/home/lacerda/workspace/PCA/src
WORKDIR=/home/lacerda/CALIFA
STARLIGHTMASKFILE=${WORKDIR}/Mask.mC
FITSDIR=${WORKDIR}/gal_fits

ordem morph: 4 119 8 277 577 213 802
GAL_ET="K0864 K0119"
GAL_SP="K0073 K0518 K0008 K0277"
GAL_ME="K0213 K0802"
GAL="$GAL_SP $GAL_ET $GAL_ME"
#GAL="K0277" 
EIGV="--eigv 1 2 3 4 10 20"
ZONES="--zones 0 10 20 100 200"
TMAX="--tmax 20"
QUBICMASKQUANTIL="--rFL 0.9"
REMSTARLIGHTEMLINES="--rSEL $STARLIGHTMASKFILE"
NUMCORREPCS="--numcorrepc 5"
RUNOPTS="-R -C -T"
RUNOPTS="-C -T"
PCAPLOT=${SRCDIR}/plots/PCAlifa-plots.py
MASKDIR=${SRCDIR}/plots


for g in $GAL
do
     [ ! -d "$g" ] && mkdir -- $g

    runname="complete_spectra"
    dir="${g}/$runname"
    [ ! -d "$dir" ] && mkdir -p -- $dir
    logfile=${dir}/info.log
    echo "python $PCAPLOT --califaID $g -d $FITSDIR -l 3800 6850 $EIGV $ZONES $TMAX $QUBICMASKQUANTIL $NUMCORREPCS $RUNOPTS --pcalog"
    time python $PCAPLOT --califaID $g -d $FITSDIR -l 3800 6850 $EIGV $ZONES $TMAX $QUBICMASKQUANTIL $NUMCORREPCS $RUNOPTS > $logfile
    mv ${g}*.png ${dir}/.

#    runname="noEmLines"
#    dir="${g}/$runname"
#    [ ! -d "$dir" ] && mkdir -p -- $dir
#    logfile=${dir}/info.log
#    echo "python $PCAPLOT --califaID $g -d $FITSDIR -l 3800 6850 $EIGV $ZONES $TMAX $QUBICMASKQUANTIL $REMSTARLIGHTEMLINES $NUMCORREPCS $RUNOPTS --pcalog --outputfname $runname"
#    time python $PCAPLOT --califaID $g -d $FITSDIR -l 3800 6850 $EIGV $ZONES $TMAX $QUBICMASKQUANTIL $REMSTARLIGHTEMLINES $NUMCORREPCS $RUNOPTS --outputfname $runname > $logfile
#    mv ${g}*.png ${dir}/.
#
#    runname="justEmLines"
#    dir="${g}/$runname"
#    [ ! -d "$dir" ] && mkdir -p -- $dir
#    logfile=${dir}/info.log
#    echo "python $PCAPLOT --pcainterval --califaID $g -d $FITSDIR -l 3800 6850 $EIGV $ZONES $TMAX $QUBICMASKQUANTIL $NUMCORREPCS $RUNOPTS --rSEL ${MASKDIR}/mask_lines --outputfname $runname"
#    time python $PCAPLOT --pcainterval --califaID $g -d $FITSDIR -l 3800 6850 $EIGV $ZONES $TMAX $QUBICMASKQUANTIL $NUMCORREPCS $RUNOPTS --rSEL ${MASKDIR}/mask_lines --outputfname $runname > $logfile
#    mv ${g}*.png ${dir}/.
#
#    runname="N2Ha"
#    dir="${g}/$runname"
#    [ ! -d "$dir" ] && mkdir -p -- $dir
#    logfile=${dir}/info.log
#    echo "python $PCAPLOT --pcainterval --califaID $g -d $FITSDIR -l 3800 6850 $EIGV $ZONES $TMAX $QUBICMASKQUANTIL $NUMCORREPCS $RUNOPTS --rSEL ${MASKDIR}/mask_N2Ha --outputfname $runname"
#    time python $PCAPLOT --pcainterval --califaID $g -d $FITSDIR -l 3800 6850 $EIGV $ZONES $TMAX $QUBICMASKQUANTIL $NUMCORREPCS $RUNOPTS --rSEL ${MASKDIR}/mask_N2Ha --outputfname $runname > $logfile
#    mv ${g}*.png ${dir}/.
#
#    runname="O3Hb"
#    dir="${g}/$runname"
#    [ ! -d "$dir" ] && mkdir -p -- $dir
#    logfile=${dir}/info.log
#    echo "python $PCAPLOT --pcainterval --califaID $g -d $FITSDIR -l 3800 6850 $EIGV $ZONES $TMAX $QUBICMASKQUANTIL $NUMCORREPCS $RUNOPTS --rSEL ${MASKDIR}/mask_O3Hb --outputfname $runname"
#    time python $PCAPLOT --pcainterval --califaID $g -d $FITSDIR -l 3800 6850 $EIGV $ZONES $TMAX $QUBICMASKQUANTIL $NUMCORREPCS $RUNOPTS --rSEL ${MASKDIR}/mask_O3Hb --outputfname $runname > $logfile
#    mv ${g}*.png ${dir}/.
#
done
