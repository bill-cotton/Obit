#!/bin/sh
# Script/Notes for building 3rd party software for Obit
# DEBUG version only builds packages if the tarball is not already present.
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
doWWW=no
doCURL=yes
doXMLRPC=yes

# Check command line arguments
arg=$1
if test $arg=-without; then
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
if test $arg=-help; then
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
echo "doCURL $doCURL"
echo "doXMLRPC $doXMLRPC"

# Base address
BASE=`pwd`; export BASE
BASE3=`pwd`/other; export BASE3

# Set compiler
CC=/usr/bin/gcc; export CC

# Obit base directory
OBIT=`pwd`/Obit;export OBIT

# Set paths
PYTHONPATH=$OBIT/python;export PYTHONPATH
LD_LIBRARY_PATH="$BASE3/lib"; export LD_LIBRARY_PATH
echo "LD_LIBRARY_PATH =$LD_LIBRARY_PATH"
PATH=$PATH:$BASE3/bin;export PATH

# Create other directories for third party software
if !(test -d other); then
    mkdir other
fi
if !(test -d other/tarballs); then
    mkdir other/tarballs
fi

# Third party software:
#plplot
cd $BASE
if test $doPLPLOT = yes; then
    plplotdir=plplot-5.8.0
    plplottar=$plplotdir".tar.gz"
    # Copy tarball if necssary
    if !(test -f other/tarballs/$plplottar); then
	wget https://svn.cv.nrao.edu/svn/ObitInstall/other/tarballs/$plplottar
	mv $plplottar other/tarballs/
	cd $BASE3
	rm -f -r $plplotdir
	tar xzvf tarballs/$plplottar
	cd  $plplotdir
	./configure --prefix=$BASE3/ -enable-java=no --enable-tcl=no \
	    --without--python --with-double=no
	make clean all
	make install
	make clean
    fi
fi

# cfitsio
cd $BASE
if test $doCFITSIO = yes; then
# cleanup
    cfitsiodir=cfitsio
    cfitsiotar=$cfitsiodir"3100.tar.gz"
    # Copy tarball if necssary
    if !(test -f other/tarballs/$cfitsiotar); then
	wget https://svn.cv.nrao.edu/svn/ObitInstall/other/tarballs/$cfitsiotar
	mv $cfitsiotar other/tarballs/
	cd $BASE3
	rm -f -r $cfitsiodir
	tar xzvf tarballs/$cfitsiotar
	cd $cfitsiodir
	./configure --prefix=$BASE3/ 
	make clean all install
	make clean
    fi
fi


# glib , needs pkgconfig
cd $BASE
if test $doGLIB = yes; then
    pkgconfigdir=pkgconfig-0.14.0
    pkgconfigtar=$pkgconfigdir".tar.gz"
    # Copy tarball if necssary
    if !(test -f other/tarballs/$pkgconfigtar); then
	wget https://svn.cv.nrao.edu/svn/ObitInstall/other/tarballs/$pkgconfigtar
	mv $pkgconfigtar other/tarballs/
	cd $BASE3
# pkgconfig
	rm -f -r $pkgconfigdir
	tar xzvf tarballs/$pkgconfigtar
	cd $pkgconfigdir
	./configure --prefix=$BASE3/
	make clean all install
    fi
# now glib
    cd $BASE
    glibdir=glib-2.2.0
    glibtar=$glibdir".tar.gz"
    # Copy tarball if necssary
    if !(test -f other/tarballs/$glibtar); then
	wget https://svn.cv.nrao.edu/svn/ObitInstall/other/tarballs/$glibtar
	mv $glibtar other/tarballs/
    cd $BASE3
    rm -f -r $glibdir
    tar xzvf tarballs/$glibtar
    cd $glibdir
    ./configure --prefix=$BASE3/
    make clean all install
    fi
# Link headers where they can be found at compile time
#ln -s $BASE3/include/glib-2.0/*.h  $BASE3/include
#ln -s $BASE3/lib/glib-2.0/include/*.h  $BASE3/include
#    make clean
fi

# fftw3
cd $BASE
if test $doFFTW = yes; then
    fftw3dir=fftw-3.1.2
    fftw3tar=$fftw3dir".tar.gz"
    # Copy tarball if necssary
    if !(test -f other/tarballs/$fftw3tar); then
	wget https://svn.cv.nrao.edu/svn/ObitInstall/other/tarballs/$fftw3tar
	mv $fftw3tar other/tarballs/
	cd $BASE3
# cleanup
	rm -f -r $fftw3dir
	tar xzvf tarballs/$fftw3tar
	cd $fftw3dir
# Use --enable-threads for multithreaded
	./configure --prefix=$BASE3/ --with-gcc --enable-shared --enable-float
	make clean all install
    fi
fi

# gsl
cd $BASE
if test $doGSL = yes; then
    gsldir=gsl-1.6
    gsltar=$gsldir".tar.gz"
    # Copy tarball if necssary
    if !(test -f other/tarballs/$gsltar); then
	wget https://svn.cv.nrao.edu/svn/ObitInstall/other/tarballs/$gsltar
	mv $gsltar other/tarballs/
	cd $BASE3
# cleanup
	rm -f -r $gsldir
	tar xzvf tarballs/$gsltar
	cd $gsldir
	./configure --prefix=$BASE3 
	make clean all install
	make clean
    fi
fi

# zlib
cd $BASE
if test $doZLIB = yes; then
    zlibdir=zlib-1.2.3
    zlibtar=$zlibdir".tar.gz"
    # Copy tarball if necssary
    if !(test -f other/tarballs/$zlibtar); then
	wget https://svn.cv.nrao.edu/svn/ObitInstall/other/tarballs/$zlibtar
	mv $zlibtar other/tarballs/
	cd $BASE3
# cleanup
	rm -f -r $zlibdir
	tar xzvf tarballs/$zlibtar
	cd $zlibdir
	./configure --prefix=$BASE3
	make clean all 
	make install
	make clean
    fi
fi

# Open Motif
cd $BASE
if test $doMOTIF = yes; then
    openmotifdir=openmotif-2.3.0
    openmotiftar=$openmotifdir".tar.gz"
    # Copy tarball if necssary
    if !(test -f other/tarballs/$openmotiftar); then
	wget https://svn.cv.nrao.edu/svn/ObitInstall/other/tarballs/$openmotiftar
	mv $openmotiftar other/tarballs/$openmotiftar
	cd $BASE3
    # cleanup
	rm -f -r $openmotifdir
	tar xzvf tarballs/$openmotiftar
	cd $openmotifdir
	./configure --prefix=$BASE3
	make clean all 
	make install
	make clean
    fi
fi

# Python
cd $BASE
if test $doPYTHON = yes; then
    Pythondir=Python-2.5.1
    Pythontar=$Pythondir".tgz"
    # Copy tarball if necssary
    if !(test -f other/tarballs/$Pythontar); then
	wget https://svn.cv.nrao.edu/svn/ObitInstall/other/tarballs/$Pythontar
	mv $Pythontar other/tarballs/$Pythontar
	cd $BASE3
# cleanup
	rm -f -r $Pythondir
	tar xzvf tarballs/$Pythontar
	cd $Pythondir
	./configure --prefix=$BASE3 --exec-prefix=$BASE3/../ --enable-shared 
	make clean all 
	make install
	make clean
    fi
fi

# WWW (must have zlib installed in a standard location)
# Using WWW is a bad idea 
cd $BASE
if test $doWWW = yes; then
    wwwdir=w3c-libwww-5.4.0
    wwwtar=$wwwdir".tgz"
    # Copy tarball if necssary
    if !(test -f other/tarballs/$wwwtar); then
	wget https://svn.cv.nrao.edu/svn/ObitInstall/other/tarballs/$wwwtar
	mv $wwwtar other/tarballs/$wwwtar
	cd $BASE3
# cleanup
	rm -f -r $wwwdir
	tar xzvf tarballs/$wwwtar
	cd $wwwdir
	./configure --prefix=$BASE3 --with-zlib=$BASE3/lib/libz.a
	make clean all 
	make install
	make clean
    fi
fi

# curl
cd $BASE
if test $doCURL = yes; then
    curldir=curl-7.17.0
    curltar=$curldir".tar.gz"
    # Copy tarball if necssary
    if !(test -f other/tarballs/$curltar); then
	wget https://svn.cv.nrao.edu/svn/ObitInstall/other/tarballs/$curltar
	mv $curltar other/tarballs/$curltar
	cd $BASE3
# cleanup
	rm -f -r $curldir
	tar xzvf tarballs/$curltar
	cd $curldir
	./configure --prefix=$BASE3 --with-zlib=$BASE3/lib
	make clean all install
# Link headers where they can be found at compile time
	ln -s $BASE3/include/curl/*.h  $BASE3/include
	make clean
    fi
fi

# xmlrpc 
# Note: this needs libwww or libcurl installed and libwww-config in the path
# Also needs libcurl
#setenv PATH "$PATH : $BASE3/bin"
cd $BASE
if test $doXMLRPC = yes; then
    xmlrpcdir=xmlrpc-c-1.06.18
    xmlrpctar=$xmlrpcdir".tgz"
    # Copy tarball if necssary
    if !(test -f other/tarballs/$xmlrpctar); then
	wget https://svn.cv.nrao.edu/svn/ObitInstall/other/tarballs/$xmlrpctar
	mv $xmlrpctar other/tarballs/$xmlrpctar
	cd $BASE3
# cleanup
	rm -f -r $xmlrpcdir
	tar xzvf tarballs/$xmlrpctar
	cd $xmlrpcdir
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
fi




