#!/bin/sh
# Script/Notes for installing 3rd party software for Obit
# Need to set LD_LIBRARY_PATH PYTHON_PATH
#------------------------------------------------------------------------------
#  Which libraries wanted
doPLPLOT=yes
doCFITSIO=yes
doGLIB=yes
doFFTW=yes
doGSL=yes
doZLIB=yes
doMOTIF=yes
doPYTHON=yes
doWWW=yes
doCURL=yes
doXMLRPC=yes

# Check command line arguments
arg=$1
if test $arg = -without; then
    for x in $@; do
	if test $x = PLPLOT;  then doPLPLOT=no;fi
	if test $x = CFITSIO; then doCFITSIO=no;fi
	if test $x = GLIB;    then doGLIB=no; fi
	if test $x = FFTW;    then doFFTW=no; fi
	if test $x = GSL;     then doGSL=no; fi
	if test $x = ZLIB;    then doZLIB=no; fi
	if test $x = MOTIF;   then doMOTIF=no; fi
	if test $x = PYTHON;  then doPYTHON=no; fi
	if test $x = WWW;     then doWWW=no; fi
	if test $x = CURL;    then doCURL=no; fi
	if test $x = XMLRPC;  then doXMLRPC=no; fi
    done
fi
if test $arg = -help; then
    echo "Build Obit software"
    echo "Third party software may be deselected with -without and any of the following"
    echo "to use a version of the package installed in a standard place"
    echo "PLPLOT - graphics library"
    echo "CFITSIO - FITS library"
    echo "GLIB - GNU extensions to c"
    echo "FFTW - Fourier transform package"
    echo "GSL - GNU Scientific Library"
    echo "ZLIB - Compression library"
    echo "MOTIF - Motif graphics library"
    echo "PYTHON - Python"
    echo "WWW - Internet protocol library"
    echo "CURL - Internet URL library"
    echo "XMLRPC - XMLRPC network protocol library"
fi
echo "doPLPLOT $doPLPLOT"
echo "doCFITSIO $doCFITSIO"
echo "doGLIB $doGLIB"
echo "doFFTW $doFFTW"
echo "doGSL $doGSL"
echo "doZLIB $doZLIB"
echo "doMOTIF $doMOTIF"
echo "doPYTHON $doPYTHON"
echo "doWWW $doWWW"
echo "doCURL $doCURL"
echo "doXMLRPC $doXMLRPC"

# Base address
BASE3=`pwd`; export BASE3

# Set compiler
CC=/usr/bin/gcc; export CC

# Obit base directory
OBIT=$BASE3/../Obit;export OBIT

# Set paths
PYTHONPATH=$OBIT/python;export PYTHONPATH
LD_LIBRARY_PATH=$BASE3/lib; export LD_LIBRARY_PATH
PATH=$PATH:$BASE3/bin;export PATH

# Third party software:
#plplot
if test $doPLPLOT = yes; then
    cd $BASE3
# cleanup
    rm -f -r plplot-5.8.0
    tar xzvf tarballs/plplot-5.8.0.tar.gz
    cd  plplot-5.8.0
    ./configure --prefix=$BASE3/ -enable-java=no --enable-tcl=no \
	--without--python --with-double=no
    make clean all
    make install
    make clean
fi

# cfitsio
if test $doCFITSIO = yes; then
# cleanup
    cd $BASE3
    rm -f -r cfitsio
    tar xzvf tarballs/cfitsio3100.tar.gz
    cd cfitsio
    ./configure --prefix=$BASE3/ 
    make clean all install
    make clean
fi


# glib , needs pkgconfig
if test $doGLIB = yes; then
    cd $BASE3
# pkgconfig
    rm -f -r pkgconfig-0.14.0
    tar xzvf tarballs/pkgconfig-0.14.0.tar.gz
    cd pkgconfig-0.14.0
    ./configure --prefix=$BASE3/
    make clean all install
# cleanup
    cd $BASE3
    rm -f -r glib-2.2.0
    tar xzvf tarballs/glib-2.2.0.tar.gz
    cd glib-2.2.0
    ./configure --prefix=$BASE3/
    make clean all install
# Link headers where they can be found at compile time
#ln -s $BASE3/include/glib-2.0/*.h  $BASE3/include
#ln -s $BASE3/lib/glib-2.0/include/*.h  $BASE3/include
#    make clean
fi

# fftw3
if test $doFFTW = yes; then
    cd $BASE3
# cleanup
    rm -f -r fftw-3.1.2
    tar xzvf tarballs/fftw-3.1.2.tar.gz
    cd fftw-3.1.2
# Use --enable-threads for multithreaded
    ./configure --prefix=$BASE3/ --with-gcc --enable-shared --enable-float
    make clean all install
fi

# gsl
if test $doGSL = yes; then
    cd $BASE3
# cleanup
    rm -f -r gsl-1.6
    tar xzvf tarballs/gsl-1.6.tar.gz
    cd gsl-1.6
    ./configure --prefix=$BASE3 
    make clean all install
    make clean
fi

# zlib
if test $doZLIB = yes; then
    cd $BASE3
# cleanup
    rm -f -r zlib-1.2.3  
    tar xzvf tarballs/zlib-1.2.3.tar.gz
    cd zlib-1.2.3
    ./configure --prefix=$BASE3
    make clean all 
    make install
    make clean
fi

# Open Motif
if test $doMOTIF = yes; then
    cd $BASE3
    # cleanup
    rm -f -r openmotif-2.3.0
    tar xzvf tarballs/openmotif-2.3.0.tar.gz
    cd openmotif-2.3.0
    ./configure --prefix=$BASE3
    make clean all 
    make install
    make clean
fi

# Python
if test $doPYTHON = yes; then
    cd $BASE3
# cleanup
    rm -f -r Python-2.5.1
    tar xzvf tarballs/Python-2.5.1.tgz
    cd Python-2.5.1
    ./configure --prefix=$BASE3 --exec-prefix=$BASE3/../ --enable-shared 
    make clean all 
    make install
    make clean
fi

# WWW (must have zlib installed in a standard location)
if test $doWWW = yes; then
    cd $BASE3
# cleanup
    rm -f -r w3c-libwww-5.4.0
    tar xzvf tarballs/w3c-libwww-5.4.0.tgz
    cd w3c-libwww-5.4.0
    ./configure --prefix=$BASE3 --with-zlib=$BASE3/lib/libz.a
    make clean all 
    make install
    make clean
fi

# curl
if test $doCURL = yes; then
    cd $BASE3
# cleanup
    rm -f -r curl-7.17.0
    tar xzvf tarballs/curl-7.17.0.tar.gz
    cd curl-7.17.0
    ./configure --prefix=$BASE3 --with-zlib=$BASE3/lib
    make clean all install
# Link headers where they can be found at compile time
    ln -s $BASE3/include/curl/*.h  $BASE3/include
    make clean
fi

# xmlrpc 
# Note: this needs libwww installed and libwww-config in the path
# Also needs libcurl
#setenv PATH "$PATH : $BASE3/bin"
if test $doXMLRPC = yes; then
    cd $BASE3
# cleanup
    rm -f -r xmlrpc-c-1.06.18
    tar xzvf tarballs/xmlrpc-c-1.06.18.tgz
    cd xmlrpc-c-1.06.18
# Optionally use ObitSystem www, curl
    curl_opt=" --disable-curl-client "
    if test $doCURL = yes; then
        curl_opt=" CURL_CONFIG=$BASE3/bin/curl-config "
    fi
    libwww_opt=" --disable-libwww-client"
    if test $doWWW = yes; then
        libwww_opt=" LIBWWW_CONFIG=$BASE3/bin/libwww-config "
    fi
    ./configure --prefix=$BASE3 --disable-cplusplus $curl_opt $libwww_opt
    make clean all install
# Patch xmlrpc install bugs
    rm -f $BASE3/include/XmlRpcCpp.h
    install-sh lib/libutil/.libs/libxmlrpc_util.so.3 $BASE3/lib
    install-sh lib/libutil/.libs/libxmlrpc_util.so.3.6.15 $BASE3/lib
    install-sh lib/libutil/.libs/libxmlrpc_util.a $BASE3/lib
    install-sh lib/libutil/.libs/libxmlrpc_util.lai $BASE3/lib
    #install-sh lib/libutil/.libs/libxmlrpc_util.la $BASE3/lib
    make clean
fi




