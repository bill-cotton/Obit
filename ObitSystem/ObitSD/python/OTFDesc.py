""" Python Obit OTFDesc Obit "On the Fly" data descriptor class.
This contains information about the observations and the size and 
structure of the data.
"""
# $Id$
#-----------------------------------------------------------------------
#  Copyright (C) 2004-2008
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

# Python shadow class to ObitOTFDesc class
import Obit, string, math

class OTFDescPtr :
    def __init__(self,this):
        self.this = this
        self.thisown = 0
     #def __del__(self):
     #    if self.thisown == 1 :
     #       # If Obit has been unloaded don't bother
     #       if Obit.__class__ == Obit:
     #           Obit.delete_OTFDesc(self.this)
    def __setattr__(self,name,value):
        if name == "me" :
            Obit.OTFDesc_me_set(self.this,value)
            return
        if name=="Dict":
            return PSetDict(self,value)
        self.__dict__[name] = value
    def __getattr__(self,name):
        if name == "me" : 
            return Obit.OTFDesc_me_get(self.this)
        # Functions to return members
        if name=="List":
            return PGetList(self)
        if name=="Dict":
            return PGetDict(self)
        raise AttributeError,name
    def __repr__(self):
        return "<C OTFDesc instance>"
class OTFDesc(OTFDescPtr):
    """ Python Obit OTFDesc Obit "On the Fly" data descriptor class.
    This contains information about the observations and the size and 
    structure of the data.
    """
    def __init__(self,name) :
        self.this = Obit.new_OTFDesc(name)
        self.thisown = 1
    def __del__(self):
        if Obit!=None:
            Obit.delete_OTFDesc(self.this)

def PCopy (inDesc, err):
    """  Copy an OTFDesc

    returns copy of the input OTFDesc
    inDesc  = Python Obit input OTFDesc
    err     = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inDesc):
        raise TypeError,"inDesc MUST be a Python Obit OTFDesc"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be a Python ObitErr"
    #
    out = OTFDesc.OTFDesc("None")
    if err.isErr: # existing error?
        return out
    out.me = Obit.OTFDescCopy(inDesc.me, out.me, err.me)
    return out
    # end PCopy

 
def PCopyDesc (inDesc, outDesc, err):
    """  Copies contenst of OTFDesc
    
    inDesc  = Python Obit input OTFDesc
    outDesc = Python Obit output OTFDesc
    err     = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inDesc):
        raise TypeError,"inDesc MUST be a Python Obit OTFDesc"
    if not PIsA(outDesc):
        raise TypeError,"outDesc MUST be a Python Obit OTFDesc"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be a Python ObitErr"
    if err.isErr: # existing error?
        return 
    #
    Obit.OTFDescCopyDesc (inDesc.me, outDesc.me, err.me)
    # end PCopyDesc


def PIndex (inDesc):
    """  Index an OTFDesc

    inDesc  = Python Obit input OTFDesc
    """
    ################################################################
    # Checks
    if not PIsA(inDesc):
        raise TypeError,"inDesc MUST be a Python Obit OTFDesc"
    #
    Obit.OTFDescIndex (inDesc.me)
    # end PIndex


def PGetDict(inDesc):
    """  Convert Obit OTFDesc to a Python Dictionary

    returns contents of an ImageDesc as a Python Dictionary
    inDesc  = Python Obit input OTFDesc
    """
    ################################################################
    # Checks
    if not PIsA(inDesc):
        raise TypeError,"inDesc MUST be a Python Obit OTFDesc"
    #
    return Obit.OTFDescGetDict(inDesc.me)
    # end PGetDict

def PSetDict(inDesc, inDict):
    """  Convert Python Dictionary to a Obit OTFDesc 

    inDesc  = Python Obit input OTFDesc to update
    """
    ################################################################
    # Checks
    if not PIsA(inDesc):
        raise TypeError,"inDesc MUST be a Python Obit OTFDesc"
    #
    Obit.OTFDescSetDict(inDesc.me, inDict)
    # end PSetDict


def PGetList (inDesc):
    """  Get InfoList from OTFDesc

    returns InfoList
    inDesc  = Python Obit input OTFDesc
    """
    ################################################################
    # Checks
    if not PIsA(inDesc):
        raise TypeError,"inDesc MUST be a Python Obit OTFDesc"
    #
    out    = InfoList.infoList()
    out.me = Obit.InfoListUnref(out.me)
    out.me = Obit.OTFDescGetList(inDesc.me)
    return out
    # end PGetList 


def PIsA (inID):
    """ Tells if the input really is a Python Obit OTFDesc

    returns true or false (1,0)
    inID = Python OTFDesc to test
    """
    ################################################################
     # Checks
    if inID.__class__ != OTFDesc:
        return 0
    #
    return Obit.OTFDescIsA(inID.me)
    # end  PIsA

def PHeader (inID):
    """ Print the contents of a descriptor

    inID   = Python ImageDesc to print
    """
    ################################################################
    # Checks
    if not PIsA(inID):
        raise TypeError,"inID MUST be a Python Obit ImageDesc"
    #
    dict = inID.Dict
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
