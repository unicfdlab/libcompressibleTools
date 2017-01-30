#!/bin/bash

source ./libEnv.sh

ls "FoamFourierAnalysis/$FFTW_LIB"
cd FoamFourierAnalysis/$FFTW_LIB
make clean
make distclean
rm -rf fftw.timeStamp
cd ../..
wclean

#
#END-OF-FILE
#

