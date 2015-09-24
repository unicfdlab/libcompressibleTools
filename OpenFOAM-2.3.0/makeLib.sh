#!/bin/bash

source ./libEnv.sh

THIS_DIR=`pwd`

cd $THIS_DIR/FoamFourierAnalysis/$FFTW_LIB

if [ ! -e fftw.timeStamp ]
then

    CFLAGS=-fPIC\\
    CXXFLAGS=-fPIC\\
    ./configure --prefix=$FOAM_USER_LIBBIN/$FFTW_LIB --enable-shared

    make
    make install
    touch fftw.timeStamp

    cd $FOAM_USER_LIBBIN
    libs=`ls $FFTW_LIB/lib/lib*.so*`
    ls $libs

    for lib in $libs
    do
	ln -s $lib
    done
fi

cd $THIS_DIR

wmake libso
#
#END-OF-FILE
#

