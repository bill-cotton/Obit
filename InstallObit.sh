#!/bin/sh
# Script/Notes for installing Obit
# See InstallObit.sh -help for usage
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
doTHIRD=yes
doObit=yes
doObitView=yes
doObitTalk=yes
doObitSD=yes

# Check command line arguments
arg=$1
if test $arg = "-without"; then
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
	if test $x = THIRD;   then doTHIRD=no; fi
	if test $x = Obit;    then doObit=no; fi
	if test $x = ObitView;  then doObitView=no; fi
	if test $x = ObitTalk;  then doObitTalk=no; fi
	if test $x = ObitSD;  then doObitSD=no; fi
    done
fi

if test $arg = "-help"; then
    echo "Build Obit software"
    echo "Third party software may be deselected with -without and any of the following"
    echo "to use a version of the package installed in a standard place"
    echo "PLPLOT - graphics library"
    echo "CFITSIO - FITS library"
    echo "GLIB - GNU extensions to c"
    echo "FFTW - Fourier transform package(3) "
    echo "GSL - GNU Scientific Library"
    echo "ZLIB - Compression library"
    echo "MOTIF - Motif graphics library "
    echo "PYTHON - Python"
    echo "WWW - Internet protocol library - Best to avoid - DISABLED"
    echo "CURL - Internet URL library"
    echo "XMLRPC - XMLRPC network protocol library"
    echo ""
    echo "Select Obit options"
    echo "THIRD - don't (re) install third party software"
    echo "Obit  - don't (re) install Obit"
    echo "ObitView - don't (re) install ObitView"
    echo "ObitTalk - don't (re) install ObitTalk"
    echo "ObitSD - don't (re) install Obit single dish package"
    exit
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
echo "doTHIRD $doTHIRD"
echo "doObit $doObit"
echo "doObitView $doObitView"
echo "doObitTalk $doObitTalk"
echo "doObitSD $doObitSD"

# Obit system Base address
BASE=`pwd`; export BASE

# Third party software
THIRD=$BASE/other; export THIRD

# Set compiler
CC=/usr/bin/gcc; export CC

# Obit base directory
OBIT=$BASE/ObitSystem/Obit;export OBIT

# Obit single dish base directory
OBITSD=$BASE/ObitSystem/ObitSD;export OBITSD

# Set paths
PYTHONPATH=$OBIT/python;export PYTHONPATH
LD_LIBRARY_PATH=$BASE/other/lib; export LD_LIBRARY_PATH
PATH=$BASE/bin:$PATH;export PATH

# install Third party software
if test $doTHIRD = yes; then
    cd $BASE/other
    ./Install3rdParty.sh $@
    # Add link to python executable if built
    cd $BASE
    if test $doPYTHON  = yes; then 
	rm -f $BASE/bin/python
	ln -s $BASE/other/bin/python $BASE/bin/python
    fi
fi

PLPLOT="--with-plplot=$THIRD"
GSL="--with-gsl=$THIRD"
GLIB="--with-glib-prefix=$THIRD"
FFTW="--with-fftw3=$THIRD"
CFITSIO="--with-cfitsio=$THIRD"
MOTIF="--with-motif=$THIRD"
WWW="--with-www=$THIRD"
CURL="--with-curl=$THIRD"
XMLRPC="--with-xmlrpc=$THIRD"
ZLIB="--with-zlib=$THIRD"
PYTHON="--with-python=$THIRD -with-python-includes=$THIRD/include/python2.7"

# Which ones wanted?
if test $doPLPLOT  = no; then PLPLOT= ; fi
if test $doCFITSIO = no; then CFITSIO= ; fi
if test $doFFTW    = no; then FFTW= ; fi
if test $doGSL     = no; then GSL= ; fi
if test $doGLIB    = no; then GLIB= ; fi
if test $doWWW     = no; then WWW= ; fi
if test $doCURL    = no; then CURL= ; fi
if test $doXMLRPC  = no; then XMLRPC= ; fi
if test $doZLIB    = no; then ZLIB= ; fi
if test $doMOTIF   = no; then MOTIF= ; fi
if test $doPYTHON  = no; then PYTHON= ; fi

# Set LD_LIBRARY_PATH for configure
if test \"x$LD_LIBRARY_PATH\"=\"x\"; then
  LD_LIBRARY_PATH=$LD_LIBRARY_PATH; export LD_LIBRARY_PATH
else
  LD_LIBRARY_PATH=$LD_LIBRARY_PATH:\$LD_LIBRARY_PATH; export LD_LIBRARY_PATH
fi


# Obit
if test $doObit = yes; then
    cd $BASE
# Note: this attempts to use system version of python
    cd ObitSystem/Obit
    echo ./configure --exec_prefix=$BASE --with-obit=$OBIT PATH=$BASE/other/bin:$PATH \
	$PLPLOT $GSL $GLIB $FFTW $CFITSIO $WWW $CURL $XMLRPC $ZLIB $PYTHON \
	OBIT=$OBIT  OBITINSTALL=$BASE LD_LIBRARY_PATH=$LD_LIBRARY_PATH \
	PKG_CONFIG_PATH=$BASE/other/lib/pkgconfig/ LDFLAGS=-Wl,-rpath,$LD_LIBRARY_PATH
    ./configure --exec_prefix=$BASE --with-obit=$OBIT PATH=$BASE/other/bin:$PATH \
	$PLPLOT $GSL $GLIB $FFTW $CFITSIO $WWW $CURL $XMLRPC $ZLIB $PYTHON \
	OBIT=$OBIT  OBITINSTALL=$BASE LD_LIBRARY_PATH=$LD_LIBRARY_PATH \
	PKG_CONFIG_PATH=$BASE/other/lib/pkgconfig/ LDFLAGS=-Wl,-rpath,$LD_LIBRARY_PATH
    make clean all 
#make clean
fi

# ObitView
if test $doObitView = yes; then
    cd $BASE
    cd ObitSystem/ObitView
    echo ./configure --exec_prefix=$BASE --with-obit=$OBIT PATH=$BASE/other/bin:$PATH \
	 $PLPLOT $GSL $GLIB $FFTW $CFITSIO $WWW $CURL $XMLRPC $ZLIB $MOTIF \
	OBIT=$OBIT OBITINSTALL=$BASE LD_LIBRARY_PATH=$LD_LIBRARY_PATH \
	PKG_CONFIG_PATH=$BASE/other/lib/pkgconfig/ LDFLAGS=-L$OBIT/lib
    ./configure --exec_prefix=$BASE --with-obit=$OBIT PATH=$BASE/other/bin:$PATH \
	 $PLPLOT $GSL $GLIB $FFTW $CFITSIO $WWW $CURL $XMLRPC $ZLIB $MOTIF \
	OBIT=$OBIT OBITINSTALL=$BASE LD_LIBRARY_PATH=$LD_LIBRARY_PATH \
	PKG_CONFIG_PATH=$BASE/other/lib/pkgconfig/ LDFLAGS=-L$OBIT/lib
    make clean all 
    make install
#make clean
fi

# ObitTalk
if test $doObitTalk = yes; then
    echo "OBIT $OBIT BASE $BASE THIRD $THIRD"
    cd $BASE
    cd ObitSystem/ObitTalk
    rm -f config.status config.log
    ./configure --prefix=$BASE/opt --exec_prefix=$BASE --with-obit=$OBIT \
	ADDPATH=$THIRD/lib OBIT=$OBIT LD_LIBRARY_PATH=$LD_LIBRARY_PATH \
	PKG_CONFIG_PATH=$BASE/other/lib/pkgconfig/
    make clean all install
    cd python;make all;cd ..
    make install
#make clean
fi

# Obit Single dish
if test $doObitSD = yes; then
    cd $BASE
    cd ObitSystem/ObitSD
    echo ./configure --exec_prefix=$BASE --with-obit=$OBIT PATH=$BASE/other/bin:$PATH \
	$PLPLOT $GSL $GLIB $FFTW $CFITSIO $WWW $CURL $XMLRPC $ZLIB $PYTHON \
	OBIT=$OBIT OBITINSTALL=$BASE LD_LIBRARY_PATH=$LD_LIBRARY_PATH \
	PKG_CONFIG_PATH=$BASE/other/lib/pkgconfig/ LDFLAGS=-L$OBIT/lib
    ./configure --exec_prefix=$BASE --with-obit=$OBIT PATH=$BASE/other/bin:$PATH \
	$PLPLOT $GSL $GLIB $FFTW $CFITSIO $WWW $CURL $XMLRPC $ZLIB $PYTHON \
	OBIT=$OBIT OBITINSTALL=$BASE LD_LIBRARY_PATH=$LD_LIBRARY_PATH \
	PKG_CONFIG_PATH=$BASE/other/lib/pkgconfig/ LDFLAGS=-L$OBIT/lib
    make clean all 
    make install
#make clean
fi

# Write setup scripts
cd $BASE
echo "# c shell version" > setup.csh
echo "# Setup environment to run Obit software" >> setup.csh
echo "setenv OBIT $OBIT" >> setup.csh
echo "setenv OBITSD $OBITSD" >> setup.csh
tstring="\"$LD_LIBRARY_PATH:\$LD_LIBRARY_PATH\""
echo "setenv LD_LIBRARY_PATH $tstring" >> setup.csh
tstring="\"$OBIT/python:$BASE/opt/share/obittalk/python/\""
if test $doObitSD = yes; then
    tstring="\"$OBITSD/python:$OBIT/python:$BASE/opt/share/obittalk/python/\""
fi
echo "setenv PYTHONPATH $tstring" >> setup.csh
tstring="\"$BASE/bin:\$PATH\""
echo "setenv PATH $tstring" >> setup.csh
echo "setenv OBITINSTALL $BASE" >> setup.csh
if test $doPLPLOT  = yes;  then echo "setenv PLPLOT_DRV_DIR $THIRD/lib/plplot5.8.0/drivers" >> setup.csh;fi
chmod +x setup.csh

# Write setup scripts
echo "#!/bin/sh" > setup.sh
echo "# Setup environment to run Obit software" >> setup.sh
echo "OBIT=$OBIT; export OBIT" >> setup.sh
echo "OBITSD=$OBITSD; export OBITSD" >> setup.sh
echo "if test \"x\$LD_LIBRARY_PATH\"=\"x\"; then" >> setup.sh
echo  "  LD_LIBRARY_PATH=$LD_LIBRARY_PATH; export LD_LIBRARY_PATH" >> setup.sh
echo "else" >> setup.sh
tstring="\"$LD_LIBRARY_PATH:\$LD_LIBRARY_PATH\""
echo  "  LD_LIBRARY_PATH=$tstring; export LD_LIBRARY_PATH" >> setup.sh
echo "fi" >> setup.sh
tstring="\"$OBIT/python:$BASE/opt/share/obittalk/python/\""
if test $doObitSD = yes; then
    tstring="\"$OBITSD/python:$OBIT/python:$BASE/opt/share/obittalk/python/\""
fi
echo "PYTHONPATH=$tstring; export PYTHONPATH" >> setup.sh
tstring="\"$BASE/bin:\$PATH\""
echo "PATH=$tstring; export PATH" >> setup.sh
echo "OBITINSTALL=$BASE; export OBITINSTALL" >> setup.sh
if test $doPLPLOT  = yes;  then echo "PLPLOT_DRV_DIR=$THIRD/lib/plplot5.8.0/drivers; export PLPLOT_DRV_DIR" >> setup.sh;fi
chmod +x setup.sh

