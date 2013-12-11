#!/bin/bash
SRCDIR=/home/lacerda/workspace/PCA/src
WORKDIR=/home/lacerda/CALIFA
STARLIGHTMASKFILE=${WORKDIR}/Mask.mC
FITSDIR=${WORKDIR}/gal_fits

#ordem morph: 4 119 8 277 577 213 802
#GAL="K0004 K0119 K0008 K0277 K0577 K0213 K0802" 
#GAL="K0277 K0802"
GAL="K0277 K0846"
EIGV="--eigv 1 2 3 4 10 20"
ZONES="--zones 0 10 20 100 200"
TMAX="--tmax 10"
QUBICMASKQUANTIL="--rFL 0.9"
REMSTARLIGHTEMLINES="--rSEL $STARLIGHTMASKFILE"
NUMCORREPCS="--numcorrepc 5"
RUNOPTS="-R -C -T --pcalog"
RUNOPTS="-C"
ts="`date +%s`"

for g in $GAL
do
     [ ! -d "$g" ] && mkdir -- $g

     runname="complete_spectra"
     echo "python PCAlifa-plots.py --califaID $g -f $FITSDIR -l 3800 6850 $EIGV $ZONES $TMAX $QUBICMASKQUANTIL $NUMCORREPCS $RUNOPTS"
     time python PCAlifa-plots.py --califaID $g -f $FITSDIR -l 3800 6850 $EIGV $ZONES $TMAX $QUBICMASKQUANTIL $NUMCORREPCS $RUNOPTS
     dir="${g}/$runname"
     [ ! -d "$dir" ] && mkdir -p -- $dir
     mv ${g}*.png ${dir}/.
     
     runname="noEmLines"
     echo "python PCAlifa-plots.py --califaID $g -f $FITSDIR -l 3800 6850 $EIGV $ZONES $TMAX $QUBICMASKQUANTIL $REMSTARLIGHTEMLINES $NUMCORREPCS $RUNOPTS --outputfname $runname"
     time python PCAlifa-plots.py --califaID $g -f $FITSDIR -l 3800 6850 $EIGV $ZONES $TMAX $QUBICMASKQUANTIL $REMSTARLIGHTEMLINES $NUMCORREPCS $RUNOPTS --outputfname $runname
     dir="${g}/$runname"
     [ ! -d "$dir" ] && mkdir -p -- $dir
     mv ${g}*.png ${dir}/.
     
     runname="justEmLines"
     echo "python PCAlifa-plots-intervals.py --califaID $g -f $FITSDIR -l 3800 6850 $EIGV $ZONES $TMAX $QUBICMASKQUANTIL $NUMCORREPCS $RUNOPTS --rSEL mask_lines --outputfname $runname"
     time python PCAlifa-plots-intervals.py --califaID $g -f $FITSDIR -l 3800 6850 $EIGV $ZONES $TMAX $QUBICMASKQUANTIL $NUMCORREPCS $RUNOPTS --rSEL mask_lines --outputfname $runname
     dir="${g}/$runname"
     [ ! -d "$dir" ] && mkdir -p -- $dir
     mv ${g}*.png ${dir}/.
     
     runname="N2Ha"
     echo "python PCAlifa-plots-intervals.py --califaID $g -f $FITSDIR -l 3800 6850 $EIGV $ZONES $TMAX $QUBICMASKQUANTIL $NUMCORREPCS $RUNOPTS --rSEL mask_N2Ha --outputfname $runname"
     time python PCAlifa-plots-intervals.py --califaID $g -f $FITSDIR -l 3800 6850 $EIGV $ZONES $TMAX $QUBICMASKQUANTIL $NUMCORREPCS $RUNOPTS --rSEL mask_N2Ha --outputfname $runname
     dir="${g}/$runname"
     [ ! -d "$dir" ] && mkdir -p -- $dir
     mv ${g}*.png ${dir}/.
     
     runname="O3Hb"
     echo "python PCAlifa-plots-intervals.py --califaID $g -f $FITSDIR -l 3800 6850 $EIGV $ZONES $TMAX $QUBICMASKQUANTIL $NUMCORREPCS $RUNOPTS --rSEL mask_O3Hb --outputfname $runname"
     time python PCAlifa-plots-intervals.py --califaID $g -f $FITSDIR -l 3800 6850 $EIGV $ZONES $TMAX $QUBICMASKQUANTIL $NUMCORREPCS $RUNOPTS --rSEL mask_O3Hb --outputfname $runname
     dir="${g}/$runname"
     [ ! -d "$dir" ] && mkdir -p -- $dir
     mv ${g}*.png ${dir}/.

     runname="N2Ha-noKin"
     echo "python PCAlifa-plots-noKin.py --califaID $g -f $FITSDIR -l 3800 6850 $EIGV $ZONES $TMAX $QUBICMASKQUANTIL $NUMCORREPCS $RUNOPTS --rSEL mask_N2Ha --outputfname $runname"
     time python PCAlifa-plots-noKin.py --califaID $g -f $FITSDIR -l 3800 6850 $EIGV $ZONES $TMAX $QUBICMASKQUANTIL $NUMCORREPCS $RUNOPTS --rSEL mask_N2Ha --outputfname $runname
     dir="${g}/$runname"
     [ ! -d "$dir" ] && mkdir -p -- $dir
     mv ${g}*.png ${dir}/.

#    runname="gas-$ts"
#    echo "python PCAlifa-gas.py --califaID $g -f $FITSDIR --outputfname $runname"
#    time python PCAlifa-gas.py --califaID $g -f $FITSDIR --outputfname $runname
#    dir="${g}/$runname"
#    [ ! -d "$dir" ] && mkdir -p -- "$dir"
#    mv ${g}*.png ${dir}/.


done
