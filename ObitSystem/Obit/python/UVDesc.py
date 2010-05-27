""" Python Obit UVDesc class

    This contains information about A UV data set

    Image Members with python interfaces:
    InfoList  - used to pass instructions to processing
                Member List 
    Dict      - (virtual) Python dictionary with contents of descriptor
                Member Dict
"""
# $Id$
#-----------------------------------------------------------------------
#  Copyright (C) 2005,2007,2008
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
 
# Python shadow class to ObitUVDesc class
import Obit, InfoList, OErr, string, math

class UVDescPtr :
    def __init__(self,this):
        self.this = this
    def __setattr__(self,name,value):
        if name == "me" :
            Obit.UVDesc_me_set(self.this,value)
            return
        if name=="Dict":
            return PSetDict(self,value)
        self.__dict__[name] = value
    def __getattr__(self,name):
        if name == "me" : 
            return Obit.UVDesc_me_get(self.this)
        # Functions to return members
        if name=="List":
            return PGetList(self)
        if name=="Dict":
            return PGetDict(self)
        raise AttributeError,str(name)
    def __repr__(self):
        return "<C UVDesc instance>"
class UVDesc(UVDescPtr):
    """ Python Obit Image descriptor class

    """
    def __init__(self, name) :
        self.this = Obit.new_UVDesc(name)
    def __del__(self):
        if Obit!=None:
            Obit.delete_UVDesc(self.this)

def PDefault (name):
    """ Default UVDesc

    returns new UVDesc
    name = optional name for object
    """
    ################################################################
    out = UVDesc("None")
    out.me = Obit.UVDescDefault(name)
    return out
    # end PDefault

def PGetDict (inUD):
    """ Returns the contents of an UVDesc as a Python Dictionary

    returns dictionary
    inUD = Python UVDesc to read
    """
    ################################################################
    # Checks
    if not PIsA(inUD):
        raise TypeError,"inUD MUST be a Python Obit UVDesc"
    #
    return Obit.UVDescGetDict(inUD.me)
    # end PGetDict


def PSetDict (inUD, inDict):
    """ Copies the contents of a Python Dictionaty to an UVDesc

    No type or dimension checking.  Not all values are writeable.
    It's best if this was created by PGetDict.
    inUD   = Python UVDesc to update
    inDict = Python dictionary with values
    """
    ################################################################
    # Checks
    if not PIsA(inUD):
        raise TypeError,"inUD MUST be a Python Obit UVDesc"
    #
    Obit.UVDescSetDict(inUD.me, inDict)
    # end PSetDict

def PGetList (inDesc):
    """  Get InfoList from UVDesc

    returns InfoList
    inDesc  = Python Obit input UVDesc
    """
    ################################################################
    # Checks
    if not PIsA(inDesc):
        raise TypeError,"inDesc MUST be a Python Obit UVDesc"
    #
    out    = InfoList.InfoList()
    out.me = Obit.InfoListUnref(out.me)
    out.me = Obit.UVDescGetList(inDesc.me)
    return out
    # end PGetList 

def PIsA (inUD):
    """ Tells if the input really is a Python Obit UVDesc

    returns true or false (1,0)
    inUD = Python UVDesc to test
    """
    ################################################################
     # Checks
    if inUD.__class__ != UVDesc:
        return 0
    #
    return Obit.UVDescIsA(inUD.me)
    # end  PIsA


def PHeader (inID):
    """ Print the contents of a descriptor

    inID   = Python ImageDesc to print
    """
    ################################################################
    # Checks
    if not PIsA(inID):
        raise TypeError,"inID MUST be a Python Obit UVDesc"
    #
    dict = inID.Dict
    PHeaderDict(dict)
    # end PHeader

def PHeaderDict (dict):
    """ Print the contents of a descriptor as python dict

    dict   = Python ImageDesc to print as python dict
    """
    ################################################################
    print "Object: %8s" % dict["object"] #"date"
    print "Observed: %8s Telescope:  %8s Created: %8s" % \
          (dict["obsdat"],dict["teles"],dict["date"])
    print "Observer: %8s   Instrument: %8s " % \
          (dict["observer"],dict["instrume"])
    print " # visibilities %10d  Sort order = %s" % \
          (dict["nvis"],dict["isort"])
    # Random parameters 
    print "Rand axes: %s %s %s %s %s" % \
          (dict["ptype"][0],dict["ptype"][1],dict["ptype"][2],\
           dict["ptype"][3],dict["ptype"][4])
    if dict["nrparm"] > 5:
        print "           %s %s %s %s %s" % \
              (dict["ptype"][5],dict["ptype"][6],dict["ptype"][7],\
               dict["ptype"][9],dict["ptype"][9])
    if dict["nrparm"] > 10:
        print "           %s %s %s %s " % \
              (dict["ptype"][10],dict["ptype"][11],dict["ptype"][12],\
               dict["ptype"][13])
    print "--------------------------------------------------------------"
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
    if dict["xshift"]!=0.0 or dict["yshift"]!=0.0:
        print "Phase shifted in X %10.3f in Y %10.3f" % \
              (dict["xshift"], dict["yshift"])
    if dict["beamMaj"]>0.0:
        print "Clean Beam %10g x %10g asec, PA %7.1f deg." % \
              (3600.0*dict["beamMaj"], 3600.0*dict["beamMin"], \
               dict["beamPA"])
    VelDef  = dict["VelDef"]
    VelDefStr = ["LSR", "Helio", "Observer"]
    VelType  = dict["VelDef"]
    VelTypeStr = ["Optical", "radio"]
    print "Rest freq %12g Vel type: %s,  wrt  %s" % \
          (dict["restFreq"], VelDefStr[VelDef-1], VelTypeStr[VelType])
    print "Alt ref value %12.5g  wrt pixel %8.2f" % \
          (dict["altRef"], dict["altCrpix"])
    # end PHeaderDict

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

def PHMS2RA (st):
    """ Convert a right ascension in hours, min, seconds to degrees 

    st  = Right ascension as "hh mm ss.ss".
    """
    ################################################################
    p = st.split()
    h = int(p[0])
    m = int(p[1])
    s = float(p[2])
    ra = 15.0 * (float(h) + float(m) / 60.0 + s/3600.0)
    return ra
    # end PHMS2RA

def PDMS2Dec (st):
    """ Convert a declination in degrees, min, seconds to degrees t

    st  = Declination  as "dd mm ss.ss".
    """
    ################################################################
    p = st.split()
    h = abs(int(p[0]))
    m = abs(int(p[1]))
    s = abs(float(p[2]))
    dec = (float(h) + float(m)/60.0 + s/3600.0)
    if st.__contains__("-"):
        dec = -dec
    return dec
    # end PDMS2Dec

def PPos (uv):
    """ Get pointing position from uv data (single source)

     uv   Obit UV data object
     returns [RA,Dec,epoch] in deg
    """
    ################################################################
    # Get descriptor dictionary
    dict = uv.Desc.Dict
    ra  =  dict["crval"][dict["jlocr"]]
    dec = dict["crval"][dict["jlocd"]]
    epoch = dict["equinox"]
    return [ra,dec,epoch]
   # end PPos

