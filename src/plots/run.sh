#!/bin/bash
SRCDIR=/home/lacerda/workspace/PCA/src
WORKDIR=/home/lacerda/CALIFA
STARLIGHTMASKFILE=${WORKDIR}/Mask.mC
FITSDIR=${WORKDIR}/gal_fits

imgoutput=pdf

opt_runall=false

opt_complete=false
opt_noemlines=true
opt_justemlines=false
opt_n2ha=false
opt_o3hb=false

ll=3850
lu=6800

numcorrepc=5

tmax=20

#ordem morph: 4 119 8 277 577 213 802

GAL_ET="K0864 K0119"
GAL_SP="K0073 K0518 K0008 K0277"
GAL_ME="K0213 K0802"
GAL="$GAL_SP $GAL_ET $GAL_ME"
GAL="K0277"

QUBICMASKQUANTIL="--rFL 0.95"
REMSTARLIGHTEMLINES="--rSEL $STARLIGHTMASKFILE"

LAMBDACONSTRAINS="-l $ll $lu"

CORROPTS="-C --numcorrepc $numcorrepc"
RECOOPTS="-R --zones 0 10 20 100 200 --eigv 1 2 3 4 10 20"
TOMOOPTS="-T --tmax $tmax"

RUNOPTS="$CORROPTS $RECOOPTS $TOMOOPTS"
RUNOPTS="$CORROPTS"

PCAPLOT=${SRCDIR}/plots/PCAlifa-plots.py
MASKDIR=${SRCDIR}/plots

for g in $GAL
do
    [ ! -d "$g" ] && mkdir -- $g
    
    $opt_runall || $opt_complete && {
        runname="complete_spectra"
        dir="${g}/$runname"
        [ ! -d "$dir" ] && mkdir -p -- $dir
        logfile=${dir}/info.log
        echo "python $PCAPLOT --califaID $g -d $FITSDIR $QUBICMASKQUANTIL $LAMBDACONSTRAINS $RUNOPTS --pcalog"
        time python $PCAPLOT --califaID $g -d $FITSDIR $QUBICMASKQUANTIL $LAMBDACONSTRAINS $RUNOPTS --pcalog > $logfile
        mv ${g}*.${imgoutput} ${dir}/.
    }

    $opt_runall || $opt_noemlines && {
        runname="noEmLines"
        dir="${g}/$runname"
        [ ! -d "$dir" ] && mkdir -p -- $dir
        logfile=${dir}/info.log
        echo "python $PCAPLOT --califaID $g -d $FITSDIR $QUBICMASKQUANTIL $LAMBDACONSTRAINS $REMSTARLIGHTEMLINES $RUNOPTS --pcalog --outputfname $runname"
        time python $PCAPLOT --califaID $g -d $FITSDIR $QUBICMASKQUANTIL $LAMBDACONSTRAINS $REMSTARLIGHTEMLINES $RUNOPTS --pcalog --outputfname $runname > $logfile
        mv ${g}*.${imgoutput} ${dir}/.
    }

    $opt_runall || $opt_justemlines && {
        runname="justEmLines"
        dir="${g}/$runname"
        [ ! -d "$dir" ] && mkdir -p -- $dir
        logfile=${dir}/info.log
        echo "python $PCAPLOT --pcainterval --califaID $g -d $FITSDIR $QUBICMASKQUANTIL $LAMBDACONSTRAINS $RUNOPTS --rSEL ${MASKDIR}/mask_lines --outputfname $runname"
        time python $PCAPLOT --pcainterval --califaID $g -d $FITSDIR $QUBICMASKQUANTIL $LAMBDACONSTRAINS $RUNOPTS --rSEL ${MASKDIR}/mask_lines --outputfname $runname > $logfile
        mv ${g}*.${imgoutput} ${dir}/.
    }

    $opt_runall || $opt_n2ha && {
        runname="N2Ha"
        dir="${g}/$runname"
        [ ! -d "$dir" ] && mkdir -p -- $dir
        logfile=${dir}/info.log
        echo "python $PCAPLOT --pcainterval --califaID $g -d $FITSDIR $QUBICMASKQUANTIL $LAMBDACONSTRAINS $RUNOPTS --rSEL ${MASKDIR}/mask_N2Ha --outputfname $runname"
        time python $PCAPLOT --pcainterval --califaID $g -d $FITSDIR $QUBICMASKQUANTIL $LAMBDACONSTRAINS $RUNOPTS --rSEL ${MASKDIR}/mask_N2Ha --outputfname $runname > $logfile
        mv ${g}*.${imgoutput} ${dir}/.
    }

    $opt_runall || $opt_o3hb && {
        runname="O3Hb"
        dir="${g}/$runname"
        [ ! -d "$dir" ] && mkdir -p -- $dir
        logfile=${dir}/info.log
        echo "python $PCAPLOT --pcainterval --califaID $g -d $FITSDIR $QUBICMASKQUANTIL $LAMBDACONSTRAINS $RUNOPTS --rSEL ${MASKDIR}/mask_O3Hb --outputfname $runname"
        time python $PCAPLOT --pcainterval --califaID $g -d $FITSDIR $QUBICMASKQUANTIL $LAMBDACONSTRAINS $RUNOPTS --rSEL ${MASKDIR}/mask_O3Hb --outputfname $runname > $logfile
        mv ${g}*.${imgoutput} ${dir}/.
    }

done
