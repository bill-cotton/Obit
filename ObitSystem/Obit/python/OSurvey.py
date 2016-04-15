"""
Python utility package for survey catalog browsers
"""

import Obit, InfoList, OErr, OPrinter, Image, UV, Table
import os
# $Id$
#-----------------------------------------------------------------------
#  Copyright (C) 2012-2016
#  Associated Universities, Inc. Washington DC, USA.
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License as
#  published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public
#  License along with this program; if not, write to the Free
#  Software Foundation, Inc., 675 Massachusetts Ave, Cambridge,
#  MA 02139, USA.
#
#  Correspondence about concerning this software should be addressed as follows:
#         Internet email: bcotton@nrao.edu.
#         Postal address: William Cotton
#                         National Radio Astronomy Observatory
#                         520 Edgemont Road
#                         Charlottesville, VA 22903-2475 USA
#-----------------------------------------------------------------------

def PVLPrint (VLTable, image, streamname, err):
    """ 
    Print the contents of a VL (survey catalog) table

    * VLTable    = VL table to print
    * image      = Image to to which VL table is attached
    * streamname = Name of stream to use, "stdout", or "stderr"
    * err        = Python Obit Error/message stack
    * window     = in not None, then edit this OWindow in the server
    """
    ################################################################
    # Checks
    if not Table.PIsA(VLTable):
        raise TypeError,"VLTable MUST be a Python Obit Table"
    Obit.OSurveyVLPrint(VLTable.me, image.me, streamname, err.me)
    # end PVLPrint

def PNVSSPrint (printer, data, err, VLVer=1, first=True, last=True):
    """ 
    Print selected contents of an NVSS (VLA 1.4 GHz D config) catalog

    Returns logical quit to indicate user wants to quit
    * printer    = Printer for output
    * data       = OData to to which VL table is attached
                   will cast from Image or UV types
                   with control parameters on the info:
        Object    string  Name of object
        equinCode long    Epoch code for output, 
                          1=>B1950, 2=>J2000, 3=>Galactic [def 2]
        Fitted    bool    If True give fitted values [def F]
        doraw     bool    If True give raw values from table, else 
                          corrected/deconvolved.[def F]
        RA        double  RA center of search, degrees in equinCode [def 0]
        Dec       double  Dec center of search, degrees in equinCode [def 0]
        Search    double  Search radius in arcsec, <= 0 => all selected. [def 15] 
        Box       double[2] RA and Dec halfwidth of search 
                          rectangle in hr,deg [0,0]
        Silent    double  Half width asec of silent search box. [def 720]
        minFlux   float   Minimum peak flux density. [def 0]
        maxFlux   float   Maximum peak flux density. [def LARGE]
        minPol    float   Minimum percent integrated polarization [def 0]
        minGlat   float   Minimum abs galactic latitude [def any]
        maxGlat   float   Minimum abs galactic latitude [def any]
    * err        = Python Obit Error/message stack
    * VLVer      = VL table version
    * first      = True if this first write to printer
    * last       = True if this is the last write to printer
    """
    ################################################################
    # Checks
    if not OPrinter.PIsA(printer):
        raise TypeError,"printer MUST be a Python Obit printer"
    # cast data if necessary
    if Image.PIsA(data) or UV.PIsA(data):
        ldata = data.cast("ObitData")
    else:
        ldata = data.me
    ret = Obit.OSurveyNVSSPrint(printer.me, ldata, VLVer, first, last, err.me)
    return ret != 0
    # end PNVSSPrint

def PVLSSPrint (printer, data, err, VLVer=1, first=True, last=True):
    """ 
    Print selected contents of an VLSS (VLA 74 MHz B config) catalog

    Returns logical quit to indicate user wants to quit
    * printer    = Printer for output
    * data       = OData to to which VL table is attached
                   will cast from Image or UV types
                   with control parameters on the info:
        Object    string  Name of object
        equinCode long    Epoch code for output, 
                          1=>B1950, 2=>J2000, 3=>Galactic [def 2]
        Fitted    bool    If True give fitted values [def F]
        doraw     bool    If True give raw values from table, else 
                          corrected/deconvolved.[def F]
        RA        double  RA center of search, degrees in equinCode [def 0]
        Dec       double  Dec center of search, degrees in equinCode [def 0]
        Search    double  Search radius in arcsec, <= 0 => all selected. [def 15] 
        Box       double[2] RA and Dec halfwidth of search 
                          rectangle in hr,deg [0,0]
        Silent    double  Half width asec of silent search box. [def 720]
        minFlux   float   Minimum peak flux density. [def 0]
        maxFlux   float   Maximum peak flux density. [def LARGE]
        minGlat   float   Minimum abs galactic latitude [def any]
        maxGlat   float   Minimum abs galactic latitude [def any]
    * err        = Python Obit Error/message stack
    * VLVer      = VL table version
    * first      = True if this first write to printer
    * last       = True if this is the last write to printer
    """
    ################################################################
    # Checks
    if not OPrinter.PIsA(printer):
        raise TypeError,"printer MUST be a Python Obit printer"
    # cast data if necessary
    if Image.PIsA(data) or UV.PIsA(data):
        ldata = data.cast("ObitData")
    else:
        ldata = data.me
    ret = Obit.OSurveyVLSSPrint(printer.me, ldata, VLVer, first, last, err.me)
    return ret != 0
    # end PVLSSPrint

def PGenPrint (printer, data, err, VLVer=1, first=True, last=True):
    """ 
    Print selected contents of a generic VL format catalog

    Returns logical quit to indicate user wants to quit
    * printer    = Printer for output
    * data       = OData to to which VL table is attached
                   will cast from Image or UV types
                   with control parameters on the info:
        Object    string  Name of object
        equinCode long    Epoch code for output, 
                          1=>B1950, 2=>J2000, 3=>Galactic [def 2]
        Fitted    bool    If True give fitted values [def F]
        doraw     bool    If True give raw values from table, else 
                          corrected/deconvolved.[def F]
        RA        double  RA center of search, degrees in equinCode [def 0]
        Dec       double  Dec center of search, degrees in equinCode [def 0]
        Search    double  Search radius in arcsec, <= 0 => all selected. [def 15] 
        Box       double[2] RA and Dec halfwidth of search 
                          rectangle in hr,deg [0,0]
        Silent    double  Half width asec of silent search box. [def 720]
        minFlux   float   Minimum peak flux density. [def 0]
        maxFlux   float   Maximum peak flux density. [def LARGE]
        minGlat   float   Minimum abs galactic latitude [def any]
        maxGlat   float   Minimum abs galactic latitude [def any]
        Calibration Parameters
        fluxScale float Flux density scaling factor, def 1.0
        biasRA    float RA position bias in deg, def 0.0
        biasDec   float Dec position bias in deg, def 0.0
        calRAEr   float Cal component of RA position error (squared), def 0.0
        calDecEr  float Cal component of Dec position error (squared), def 0.0
        ClnBiasAv float Mean CLEAN bias in Jy, def 0.0
        ClnBiasEr float Uncertainty in CLEAN bias in Jy, def 0.0
        calAmpEr  float Cal component of amplitude error as fraction, def 0.03
        calSizeEr float Cal component of amplitude error as fraction, def 0.02
        calPolEr  float Cal component of polarization error, def 0.003
    * err        = Python Obit Error/message stack
    * VLVer      = VL table version
    * first      = True if this first write to printer
    * last       = True if this is the last write to printer
    """
    ################################################################
    # Checks
    if not OPrinter.PIsA(printer):
        raise TypeError,"printer MUST be a Python Obit printer"
    # cast data if necessary
    if Image.PIsA(data) or UV.PIsA(data):
        ldata = data.cast("ObitData")
    else:
        ldata = data.me
    ret = Obit.OSurveyGenPrint(printer.me, ldata, VLVer, first, last, err.me)
    return ret != 0
    # end PGenPrint

def PGetCalParms(beamMaj, beamMin, beamPA, \
                     fluxScale=1.0,biasRA=0.0,biasDec=0.0,calRAEr=0.0,calDecEr=0.0, \
                     ClnBiasAv=0.0,ClnBiasEr=0.0,calAmpEr=0.03,calSizeEr=0.02,calPolEr=0.003):
    """ 
    Collect calibration parameters needed for PGenCorErr

    * beamMaj   CLEAN restoring beam major axis in deg
    * beamMin   CLEAN restoring beam minor axis in deg
    * beamPA    CLEAN restoring beam position angle axis in deg
    * fluxScale Flux density scaling factor
    * biasRA    RA position bias in deg
    * biasDec   Dec position bias in deg
    * calRAEr   Cal component of RA position error (squared)
    * calDecEr  Cal component of Dec position error (squared)
    * ClnBiasAv Mean CLEAN bias in Jy
    * ClnBiasEr Uncertainty in CLEAN bias in Jy
    * calAmpEr  Cal component of amplitude error as fraction
    * calSizeEr Cal component of amplitude error as fraction
    * calPolEr  Cal component of polarization error
    Returns dict with Calibration Parameters:
        fluxScale Flux density scaling factor
        biasRA    RA position bias in deg
        biasDec   Dec position bias in deg
        calRAEr   Cal component of RA position error (squared)
        calDecEr  Cal component of Dec position error (squared)
        ClnBiasAv Mean CLEAN bias in Jy
        ClnBiasEr Uncertainty in CLEAN bias in Jy
        calAmpEr  Cal component of amplitude error as fraction
        calSizeEr Cal component of amplitude error as fraction
        calPolEr  Cal component of polarization error
        beamMaj   CLEAN restoring beam major axis in deg
        beamMin   CLEAN restoring beam minor axis in deg
        beamPA    CLEAN restoring beam position angle axis in deg
     """
    ################################################################
    out = { \
        "fluxScale":fluxScale, \
        "biasRA":biasRA, \
        "biasDec":biasDec, \
        "calRAEr":calRAEr, \
        "calDecEr":calDecEr, \
        "ClnBiasAv":ClnBiasAv, \
        "ClnBiasEr":ClnBiasEr, \
        "calAmpEr":calAmpEr, \
        "calSizeEr":calSizeEr, \
        "calPolEr":calPolEr, \
        "beamMaj":beamMaj, \
        "beamMin":beamMin, \
        "beamPA":beamPA, \
        }
    return out
   # end PGetCalParms

def PGenCorErr(VLrow,calParms):
     """ 
    (re) calibration and error analysis of VL table row

    * VLrow      = VL table row dict from ReadRow
    * calParms   = Calibration parameters dict from PGetCalParms
    Returns dict:
      RA       Right Ascention in Deg
      eRA      Error in Right Ascention in asec
      Dec      Declination in deg
      eDec     Error in Declination in asec
      peak     Peak flux density in Jy
      flux     Integrated flux density in Jy
      eflux    Error in Integrated flux density in Jy
      major    Fitted Gaussian major axis size in deg
      minor    Fitted Gaussian minor axis size in deg
      posang   Fitted Gaussian pposition angle in deg
      pflux    Central polarized flux density in Jy
      epflux   Error in Central polarized flux density in Jy
      pang     Polarization angle deg (fblank if missing)
      epang    Error in polarization angle deg (fblank if missing)
      dmajor   Deconvolved Gaussian major axis in asec
      edmajor  Error in deconvolved Gaussian major axis in asec
      dminor   Deconvolved Gaussian minor axis in asec
      edminor  Error in deconvolved Gaussian minor axis in asec
      dpa      Deconvolved gaussian position angle (deg)
      edpa     Error in deconvolved position angle (deg)
      rflag    Resolution flags[4], 
               0=any resolution, 1=1 axis, 2=major only, 3=minor only
    """
    ################################################################
     out = Obit.OSurveyGenCorErr(VLrow,calParms)
     # convert flags to booleans
     for i in range(0,len(out['rflag'])):
         out['rflag'][i] = out['rflag'][i]!=0
     return out
# end PGenCorErr


import urllib, urllib2
def PWebFetch (url, args, outfile, err):
    """ 
    Generic Fetch a web product

    Sends a http post request to url and sames response to outfile
    Throws exception on error
    * url     = where the request to be sent,
                e.g."http://www.cv.nrao.edu/cgi-bin/postage.pl"
    * args    = dict or arguments, e.g. {"arg1":"arg1"}
    * outfile = Name of the output file, absolute path or relative to CWD
                None => use name from server
    * err     = Python Obit Error/message stack
    """
    ################################################################
    # Package parameters
    encoded_args = urllib.urlencode(args)
    NVSShost    = "http://www.cv.nrao.edu/cgi-bin/postage.pl"
    # fetch
    try:
        request      = urllib2.Request(url)
        response     = urllib2.urlopen(request, encoded_args)
        data         = response.read()
    except Exception, exception:
        print exception
        OErr.PLog(err, OErr.Error, "Request from server failed")
        OErr.printErrMsg(err)
    if outfile == None:     # Name from server?
        outfile = os.path.basename(response.headers['URI'])
    fd = open(outfile,"wb")
    fd.write(data)
    fd.close()
    # Info
    print "Response code =",response.code, response.msg
    print "Response type =",response.headers["Content-Type"]
# end PWebFetch

def PNVSSFetch (RA, Dec, outfile, err, \
                Type = 'image/xfits', Equinox = '2000', \
                ObjName='unnamed', Size = '0.50 0.50', \
                Poltype = 'I', MAPROJ = 'SIN',
                rotate = 0.0, Cells=[15.,15.]):
    """ 
    Postage Stamp server for the NRAO/VLA Sky Survey (1.4 GHz)

    Fetch a postage stamp (AKA cutout) FITS image, jpeg or contour from the NVSS
    If a FITS file, returns an ObitImage object pointing to the result,
    otherwise the output file name
    If something goes wrong with the request or download, the output file
    will likely contain an HTML file explaining the results.
    Throws exception on error
    * RA      = Right Ascension as 'hh mm ss.sss'
    * Dec     = Declination as 'sdd mm ss.ss'
    * outfile = Name of the output file, absolute path or relative to CWD
                None => use name from server
    * err     = Python Obit Error/message stack
    * Type    = Type = "type" of result;
                'image/jpeg' = jpeg image,
                'application/octet-stream' = fits file
                'image/x-fits' = fits image
                (default = /image/jpeg)
    * Equinox = coordinate equinox ('2000' or '1950')
    * ObjName = object name - used for labeling only (optional)
    * Size    = field size (RA,dec) in deg. (default 0.50 0.50)
    * Poltype = total intensity 'I' or 'IQU' (default 'I')
    * MAPROJ  = Map projection type: 'SIN', 'TAN', 'ARC', 'NCP',
                'GLS', 'MER', 'AIT','STG'
    * rotate  = Rotation on sky in deg in sense of from North to East
    * Cells   = pixel spacing in aseconds (default 15. 15.)
    """
    ################################################################
    # Package parameters
    query_args = {'Type'   :Type,
                  'Equinox':Equinox,
                  'ObjName':ObjName,
                  'RA'     :RA,
                  'Dec'    :Dec,
                  'Size'   :Size,
                  'Poltype':Poltype,
                  'MAPROJ' :MAPROJ,
                  'rotate' :str(rotate),
                  'Cells'  :str(Cells[0])+' '+str(Cells[1])}
    NVSSURL    = "http://www.cv.nrao.edu/cgi-bin/postage.pl"
    # fetch
    PWebFetch (NVSSURL, query_args, outfile, err)
    # Get fits image if requested
    if (Type == 'image/x-fits') or (Type == 'application/octet-stream'):
        #print "Get FITS"
        outfits = Image.newPFImage(ObjName, outfile, 0, True, err)
        OErr.printErrMsg(err, "Problem with FITS image, see contents of outfile")
    else:
        return outfile
    return outfits
# end PNVSSFetch

def PVLSSFetch (RA, Dec, outfile, err, \
                Type = 'image/xfits', Equinox = '2000', \
                ObjName='unnamed', Size = '1.0 1.0', \
                MAPROJ = 'SIN', rotate = 0.0, Cells=[20.,20.]):
    """
    Postage Stamp server for the VLA Low-Frequency Sky Survey (74 MHz)
    
    Fetch a postage stamp (AKA cutout) FITS image, jpeg or contour from the VLSS
    If a FITS file, returns an ObitImage object pointing to the result,
    otherwise the output file name
    If something goes wrong with the request or download, the output file
    will likely contain an HTML file explaining the results.
    Throws exception on error
    * RA      = Right Ascension as 'hh mm ss.sss'
    * Dec     = Declination as 'sdd mm ss.ss'
    * outfile = Name of the output file, absolute path or relative to CWD
                None => use name from server
    * err     = Python Obit Error/message stack
    * Type    = Type = "type" of result;
                'image/jpeg' = jpeg image,
                'application/octet-stream' = fits file
                'image/x-fits' = fits image
                (default = /image/jpeg)
    * Equinox = coordinate equinox ('2000' or '1950')
    * ObjName = object name - used for labeling only (optional)
    * Size    = field size (RA,dec) in deg. (default 1.0 1.0)
    * MAPROJ  = Map projection type: 'SIN', 'TAN', 'ARC', 'NCP',
                'GLS', 'MER', 'AIT','STG'
    * rotate  = Rotation on sky in deg in sense of from North to East
    * Cells   = pixel spacing in aseconds (default 20. 20.)
    """
    ################################################################
    # Package parameters
    query_args = {'Type'   :Type,
                  'Equinox':Equinox,
                  'ObjName':ObjName,
                  'RA'     :RA,
                  'Dec'    :Dec,
                  'Size'   :Size,
                  'MAPROJ' :MAPROJ,
                  'rotate' :str(rotate),
                  'Cells'  :str(Cells[0])+' '+str(Cells[1])}
    VLSSURL    = "http://www.cv.nrao.edu/cgi-bin/VLSSpostage.pl"
    # fetch
    PWebFetch (VLSSURL, query_args, outfile, err)
    # Get fits image if requested
    if (Type == 'image/x-fits') or (Type == 'application/octet-stream'):
        outfits = Image.newPFImage(ObjName, outfile, 0, True, err)
        OErr.printErrMsg(err, "Problem with FITS image, see contents of outfile")
    else:
        return outfile
    return outfits
# end PVLSSFetch

