""" Python Obit OTFArrayGeom Obit "On the Fly" array geometry class.
This contains information about the detector array geometry and related
"""
# $Id$
#-----------------------------------------------------------------------
#  Copyright (C) 2009
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
#  Correspondence concerning this software should be addressed as follows:
#         Internet email: bcotton@nrao.edu.
#         Postal address: William Cotton
#                         National Radio Astronomy Observatory
#                         520 Edgemont Road
#                         Charlottesville, VA 22903-2475 USA
#-----------------------------------------------------------------------

# Python shadow class to ObitOTFArrayGeom class
import Obit, string, math

class OTFArrayGeomPtr :
    def __init__(self,this):
        self.this = this
        self.thisown = 0
     #def __del__(self):
     #    if self.thisown == 1 :
     #       # If Obit has been unloaded don't bother
     #       if Obit.__class__ == Obit:
     #           Obit.delete_OTFArrayGeom(self.this)
    def __setattr__(self,name,value):
        if name == "me" :
            Obit.OTFArrayGeom_me_set(self.this,value)
            return
        if name=="Dict":
            return PSetDict(self,value)
        self.__dict__[name] = value
    def __getattr__(self,name):
        if name == "me" : 
            return Obit.OTFArrayGeom_me_get(self.this)
        # Functions to return members
        if name=="List":
            return PGetList(self)
        if name=="Dict":
            return PGetDict(self)
        raise AttributeError,name
    def __repr__(self):
        return "<C OTFArrayGeom instance>"
class OTFArrayGeom(OTFArrayGeomPtr):
    """ Python Obit OTFArrayGeom Obit "On the Fly" array geometry class.
    This contains information about the detector array geometry and related
    """
    def __init__(self,name) :
        self.this = Obit.new_OTFArrayGeom(name)
        self.thisown = 1
    def __del__(self):
        if Obit!=None:
            Obit.delete_OTFArrayGeom(self.this)
            
    def ParAng (self, time, ra, dec):
        """ Get parallactic angle for a given time, direction

        return parallactic angle in degrees 
        self   = Python Obit OTF object
        time   = Time (days)
        ra     = RA (deg)
        dec    = Declination (deg)
        """
        return Obit.OTFArrayGeomParAng (self.me, time, ra, dec)
        # end ParAng 
        
            
    def Elev (self, time, ra, dec):
        """ Get elevation angle for a given time, direction

        return elevation angle in degrees 
        self   = Python Obit OTF object
        time   = Time (days)
        ra     = RA (deg)
        dec    = Declination (deg)
        """
        return Obit.OTFArrayGeomElev (self.me, time, ra, dec)
        # end Elev 
            
    def Coord (self, raPoint, decPoint, rot):
        """ Determine celestial coordinates of each detector in the array

        return array of arrays ([[ra],[dec]] )
        self     = Python Obit OTF object
        raPoint  = RA (deg)
        decPoint = Declination (deg)
        rot      = Array rotation (parallactic angle ) (deg)
        """
        return Obit.OTFArrayGeomCoord(self.me, raPoint, decPoint, rot)
        # end Coord
            
    def Proj (self, raPoint, decPoint, rot, raProj, decProj, Proj=0):
        """ Project the locations of the array detectors onto a specified plane.

        return array of arrays ([[x_offset],[y_offset]] )
        self     = Python Obit OTF object
        raPoint  = Telescope pointing RA (deg)
        decPoint = Telescope pointing Declination (deg)
        rot      = Array rotation (parallactic angle ) (deg)
        raProj   = Central RA of desired projection (deg)
        decProj  = Central Dec of desired projection (deg)
        Proj     = Projection code:
                  0=OBIT_OTF_SIN, 1=OBIT_OTF_ARC, 2=OBIT_OTF_TAN
        """
        return Obit.OTFArrayGeomProj(self.me, raPoint, decPoint, rot, raProj, decProj, Proj)
        # end Proj
    # end class OTFArrayGeom

# Symbolic names for access codes
OBIT_OTF_SIN = 0
OBIT_OTF_ARC = 1
OBIT_OTF_TAN = 2

def PCopy (inAGeom, err):
    """  Copy an OTFArrayGeom

    returns copy of the input OTFArrayGeom
    inAGeom  = Python Obit input OTFArrayGeom
    err     = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inAGeom):
        raise TypeError,"inAGeom MUST be a Python Obit OTFArrayGeom"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be a Python ObitErr"
    #
    out = OTFArrayGeom.OTFArrayGeom("None")
    if err.isErr: # existing error?
        return out
    out.me = Obit.OTFArrayGeomCopy(inAGeom.me, out.me, err.me)
    return out
    # end PCopy

 
def PCopyAGeom (inAGeom, outAGeom, err):
    """  Copies contenst of OTFArrayGeom
    
    inAGeom  = Python Obit input OTFArrayGeom
    outAGeom = Python Obit output OTFArrayGeom
    err     = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inAGeom):
        raise TypeError,"inAGeom MUST be a Python Obit OTFArrayGeom"
    if not PIsA(outAGeom):
        raise TypeError,"outAGeom MUST be a Python Obit OTFArrayGeom"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be a Python ObitErr"
    if err.isErr: # existing error?
        return 
    #
    Obit.OTFArrayGeomCopyAGeom (inAGeom.me, outAGeom.me, err.me)
    # end PCopyAGeom


def PIndex (inAGeom):
    """  Index an OTFArrayGeom

    inAGeom  = Python Obit input OTFArrayGeom
    """
    ################################################################
    # Checks
    if not PIsA(inAGeom):
        raise TypeError,"inAGeom MUST be a Python Obit OTFArrayGeom"
    #
    Obit.OTFArrayGeomIndex (inAGeom.me)
    # end PIndex


def PGetDict(inAGeom):
    """  Convert Obit OTFArrayGeom to a Python Dictionary

    returns contents of an ImageAGeom as a Python Dictionary
    inAGeom  = Python Obit input OTFArrayGeom
    """
    ################################################################
    # Checks
    if not PIsA(inAGeom):
        raise TypeError,"inAGeom MUST be a Python Obit OTFArrayGeom"
    #
    return Obit.OTFArrayGeomGetDict(inAGeom.me)
    # end PGetDict

def PSetDict(inAGeom, inDict):
    """  Convert Python Dictionary to a Obit OTFArrayGeom 

    inAGeom  = Python Obit input OTFArrayGeom to update
    """
    ################################################################
    # Checks
    if not PIsA(inAGeom):
        raise TypeError,"inAGeom MUST be a Python Obit OTFArrayGeom"
    #
    Obit.OTFArrayGeomSetDict(inAGeom.me, inDict)
    # end PSetDict


def PGetList (inAGeom):
    """  Get InfoList from OTFArrayGeom

    returns InfoList
    inAGeom  = Python Obit input OTFArrayGeom
    """
    ################################################################
    # Checks
    if not PIsA(inAGeom):
        raise TypeError,"inAGeom MUST be a Python Obit OTFArrayGeom"
    #
    out    = InfoList.infoList()
    out.me = Obit.InfoListUnref(out.me)
    out.me = Obit.OTFArrayGeomGetList(inAGeom.me)
    return out
    # end PGetList 


def PCorrPoint (azOff, elOff, pa):
    """ Determine the Offset of the pointing position in az and el

    (Not sure this does anything useful)
    return (raPoint,decPoint) as a list
    azOff   =  azimuth ( xcos el) in deg.
    eloff   =  elevation offset in deg
    pa      =  parallactic angle (+any feed rotation)
    """
    ################################################################
    return OTFArrayGeomCorrPoint (azOff, elOff, pa)
    # end PZap

def PIsA (inAG):
    """ Tells if the input really is a Python Obit OTFArrayGeom

    returns true or false (1,0)
    inAG = Python OTFArrayGeom to test
    """
    ################################################################
     # Checks
    if inAG.__class__ != OTFArrayGeom:
        return 0
    #
    return Obit.OTFArrayGeomIsA(inAG.me)
    # end  PIsA

def PHeader (inAG):
    """ Print the contents of a descriptor

    inAG   = Python ImageAGeom to print
    """
    ################################################################
    # Checks
    if not PIsA(inAG):
        raise TypeError,"inAG MUST be a Python Obit ImageAGeom"
    #
    dict = inAG.Dict
    print "Object: %8s" % dict["object"] #"date"
    print "Observed: %8s Telescope:  %8s Created: %8s" % \
          (dict["obsdat"],dict["teles"],dict["date"])
    print "Diameter: %5.1f m  Beamsize %10.5f deg" % \
          (dict["diameter"],dict["beamSize"])
    print " # samples %10d  Sort order = %s" % \
          (dict["nrecord"],dict["isort"])
    print "--------------------------------------------------------------"
    # Columns - not in descriptor!!! :'{
    #print "Column  Format  Units"
    print "Type    Pixels   Coord value     at Pixel     Coord incr   Rotat"
    i = -1
    for ctype in dict["ctype"]:
        i = i+1
        if ctype != "        ":
            # Conversion on some types
            stuff =  PPoslabel (ctype, dict["crval"][i], dict["cdelt"][i])
            print "%8s%6d%16s%11.2f%15s%8.2f" % \
                  (ctype, dict["inaxes"][i], stuff["crval"], dict["crpix"][i], \
                  stuff["cdelt"] , dict["crota"][i])
    print "--------------------------------------------------------------"
    print "Coordinate equinox %6.1f  Coordinate epoch %7.2f" % \
          (dict["equinox"], dict["epoch"])
    print "Observed RA %16s Observed Dec %15s" % \
          (PRA2HMS(dict["obsra"]),  PDec2DMS(dict["obsdec"]))
    #VelDef  = dict["VelDef"]
    #VelDefStr = ["LSR", "Helio", "Observer"]
    #VelType  = dict["VelDef"]
    #VelTypeStr = ["Optical", "radio"]
    #print "Rest freq %12g Vel type: %s,  wrt  %s" % \
    #      (dict["restFreq"], VelDefStr[VelDef-1], VelTypeStr[VelType])
    #print "Alt ref value %12.5g  wrt pixel %8.2f" % \
    #      (dict["altRef"], dict["altCrpix"])
    # end PHeader

def PPoslabel (ctype, crval, cdelt):
    """ Convert a coordinate for display

    returns dict with entries "ctype", "crval", "cdelt"
    giving the relevant strings to display
    ctype   = coordinate type (e.g. "RA---SIN")
    crval   = coordinate value
    cdelt   = coordinate increment
    """
    ################################################################
    out = {"ctype":ctype}
    if ctype[0:2]=="RA":
        out["crval"] = PRA2HMS (crval)
        out["cdelt"] = "%15.6g" % (3600.0*cdelt)
    elif ctype[0:3]=="DEC":
        out["crval"] = PDec2DMS (crval)
        out["cdelt"] = "%15.6g" % (3600.0*cdelt)
    elif ctype[0:6]=="STOKES":
        if crval == 1.0:
            out["crval"] = "      IPol      "
        elif crval == 2.0:
            out["crval"] = "      QPol      "
        elif crval == 3.0:
            out["crval"] = "      UPol      "
        elif crval == 4.0:
            out["crval"] = "      VPol      "
        elif crval == -1.0:
            out["crval"] = "      RPol      "
        elif crval == -2.0:
            out["crval"] = "      LPol      "
        else:
            out["crval"] = "%16.5g" % crval
        out["cdelt"] = "%15.6g" % cdelt
    else:
        out["crval"] = "%16.5g" % crval
        out["cdelt"] = "%15.6g" % cdelt
        
    
    return out
    # end PPoslabel

def PRA2HMS (ra):
    """ Convert a right ascension in degrees to hours, min, seconds

    ra   = Right ascension in deg.
    """
    ################################################################
    p = ra / 15.0
    h = int(p)
    p = (p - h) * 60.0
    m = int(p)
    s = (p - m) * 60.0
    out = "  %2d %2d %8.5f" % (h,m,s)
    return out
    # end PRA2HMS

def PDec2DMS (dec):
    """ Convert a declination in degrees to degrees, min, seconds

    dec  = Declination in deg.
    """
    ################################################################
    p = math.fabs(dec)
    if dec>0.0:
        sgn = " "
    else:
        sgn = "-"
    d = int (p)
    p = (p - d) * 60.0
    m = int(p)
    s = (p - m) * 60.0
    out = "%s%2.2d %2d %7.4f " % (sgn, d,m,s)
    return out
    # end PDec2DMS

def PHMS2RA (rast,  sep=":"):
    """ Convert a right ascension string to degrees

    return RA in degrees
    rast      RA string as "hh:mm:ss.s"
    sep       sympol to use to separate components instead of ":"
    """
    ################################################################
    pp = rast.split(sep)
    if len(pp)>0:
        hour = int(pp[0])
    else:
        hour = 0
    if len(pp)>1:
        min = int(pp[1])
    else:
        min = 0
    if len(pp)>2:
        ssec = float(pp[2])
    else:
        ssec = 0.0
    ra =  hour + min/60.0 + ssec/3600.0
    return ra*15.0
    # end PHMS2RA

def PDMS2Dec (decst, sep=":"):
    """ Convert a declination string to degrees

    Returns dec in deg
    decst     Dec string as "dd:mm:ss.s"
    sep       sympol to use to separate components instead of ":"
    """
    ################################################################
    pp = decst.split(sep)
    if len(pp)>0:
        deg = int(pp[0])
    else:
        deg = 0
    if len(pp)>1:
        min = int(pp[1])
    else:
        min = 0
    if len(pp)>2:
        ssec = float(pp[2])
    else:
        ssec = 0.0
    dec =  abs(deg) + min/60.0 + ssec/3600.0
    if pp[0].find("-") >=0:
        dec = -dec
    return dec
    # end PDec2DMS
