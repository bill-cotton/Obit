# $Id$
#-----------------------------------------------------------------------
#  Copyright (C) 2004,2005
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
 
# Python shadow class to ObitImageDesc class
import Obit, InfoList, OErr, string, math

class ImageDescPtr :
    def __init__(self,this):
        self.this = this
    def __setattr__(self,name,value):
        if name == "me" :
            Obit.ImageDesc_me_set(self.this,value)
            return
        if name=="Dict":
            return PSetDict(self,value)
        self.__dict__[name] = value
    def __getattr__(self,name):
        if name == "me" : 
            return Obit.ImageDesc_me_get(self.this)
        # Functions to return members
        if name=="List":
            return PGetList(self)
        if name=="Dict":
            return PGetDict(self)
        raise AttributeError,str(name)
    def __repr__(self):
        return "<C ImageDesc instance>"
class ImageDesc(ImageDescPtr):
    """ Python Obit Image descriptor class

    This contains information about the observations and the size and 
    coordinates in the image.
    Also included are the current location of the image in an ObitImage
    image buffer and the specified subimaging parameters.

    Image Members with python interfaces:
    InfoList  - used to pass instructions to processing
                Member List (readonly)
    Dict      - (virtual) Python dictionary with contents of descriptor
                Member Dict
    """
    def __init__(self, name) :
        self.this = Obit.new_ImageDesc(name)
    def __del__(self):
        if Obit!=None:
            Obit.delete_ImageDesc(self.this)

def PDefault (name):
    """ Default ImageDesc

    returns new ImageDesc
    name = optional name for object
    """
    ################################################################
    out = ImageDesc("None")
    out.me = Obit.ImageDescDefault(name)
    return out
    # end PDefault

def PCopyDesc (inID, outID, err):
    """ Copy the descriptive information from one descriptor to another

    Structural values not copied.
    inID    = Python Obit ImageDesc for input
    outID   = Python Obit ImageDesc for output
    err     = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inID):
        raise TypeError,"inID MUST be a Python Obit ImageDesc"
    if not PIsA(outID):
        raise TypeError,"outID MUST be a Python Obit ImageDesc"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    Obit.ImageDescCopyDesc (inID.me, outID.me, err.me)
    if err.isErr:
        printErrMsg(err, "Error copying Image descriptor")
    # end PCopyDesc

def POverlap (inID1, inID2, err):
    """ Determine if there is any overlap between images

    Compares ImageDesc objects to see if the associated
    images overlap on the sky.
    Returns True if so sles False
    inID1   = First Python Obit ImageDesc for test
    inID2   = Second Python Obit ImageDesc for test
    err     = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inID1):
        raise TypeError,"inID1 MUST be a Python Obit ImageDesc"
    if not PIsA(inID2):
        raise TypeError,"inID2 MUST be a Python Obit ImageDesc"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    #
    res = Obit.ImageDescOverlap (inID1.me, inID2.me, err.me)
    if err.isErr:
        printErrMsg(err, "Error determining overlap")
    return res!=0
    # end POverlap

def PCvtPixel (inID, inPixel, outID, err):
    """ Return pixel location in outID corresponding to pixel inPixel in inID

    returns array of 2 floats giving pixel position in outID.
    inID    = Python Obit ImageDesc for input
    inPixel = array of floats giving position in image described by inID
              only first 2 used.
    outID   = Python Obit ImageDesc for output
    err     = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inID):
        raise TypeError,"inID MUST be a Python Obit ImageDesc"
    if not PIsA(outID):
        raise TypeError,"outID MUST be a Python Obit ImageDesc"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    if inPixel[0].__class__ != float:
        print "class is" , inPixel[0].__class__
        raise TypeError,"inPixel MUST be float"
    if len(inPixel) < 2:
        raise RuntimeError,"inPixel has fewer then 2 entries"
    #
    outTmp = Obit.ImageDescCvtPixel (inID.me, outID.me, inPixel, err.me)
    if err.isErr:
        printErrMsg(err, "Error converting pixel location")
    out = outTmp[0:2]
    return out
    # end PCvtPixel

def PGetPixel (inID, inPos, err):
    """ Return pixel location in inID corresponding position inPos

    returns array of 2 floats giving pixel position in outID.
    inID    = Python Obit ImageDesc for input
    inPos   = array of floats (RA, dec) in deg
    err     = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inID):
        raise TypeError,"inID MUST be a Python Obit ImageDesc"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    if inPos[0].__class__ != float:
        print "class is" , inPos[0].__class__
        raise TypeError,"inPos MUST be float"
    if len(inPos) < 2:
        raise RuntimeError,"inPos has fewer then 2 entries"
    #
    outTmp = Obit.ImageDescGetPixel (inID.me, inPos, err.me)
    if err.isErr:
        printErrMsg(err, "Error determining pixel")
    out = outTmp[0:2]
    return out
    # end PGetPixel

def PGetPos (inID, inPixel, err):
    """ Return position of pixel inPixel in inID

    returns array of 2 floats giving (RA, Dec) in deg of pixel inPixel
    inID    = Python Obit ImageDesc for input
    inPixel = array of floats giving position in image described by inID
              only first 2 used.
    err     = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inID):
        raise TypeError,"inID MUST be a Python Obit ImageDesc"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    if inPixel[0].__class__ != float:
        print "class is" , inPixel[0].__class__
        raise TypeError,"inPixel MUST be float"
    if len(inPixel) < 2:
        raise RuntimeError,"inPixel has fewer then 2 entries"
    #
    outTmp = Obit.ImageDescGetPos (inID.me, inPixel, err.me)
    if err.isErr:
        printErrMsg(err, "Error converting pixel location to position")
    out = outTmp[0:2]
    return out
    # end PGetPos

def PGetDict (inID):
    """ Returns the contents of an ImageDesc as a Python Dictionary

    returns dictionary
    inID = Python ImageDesc to read
    """
    ################################################################
    # Checks
    if not PIsA(inID):
        raise TypeError,"inID MUST be a Python Obit ImageDesc"
    #
    return Obit.ImageDescGetDict(inID.me)
    # end PGetDict


def PSetDict (inID, inDict):
    """ Copies the contents of a Python Dictionary to an ImageDesc

    No type or dimension checking.  Not all values are writeable.
    It's best if this was created by PGetDict.
    inID   = Python ImageDesc to update
    inDict = Python dictionary with values
    """
    ################################################################
    # Checks
    if not PIsA(inID):
        raise TypeError,"inID MUST be a Python Obit ImageDesc"
    #
    Obit.ImageDescSetDict(inID.me, inDict)
    # end PSetDict

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
    print "Observer: %8s   Instrument: %8s " % \
          (dict["observer"],dict["instrume"])
    print "Minimum = %12.5g  Maximum =  %12.5g %8s" % \
          (dict["minval"],dict["maxval"],dict["bunit"])
    if dict["areBlanks"]:
        print "Image contains magic value blanking"
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
    if dict["niter"]>0:
        print "no. Comp %8d" % dict["niter"]
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
    # end PHeader

def PGetList (inDesc):
    """  Get InfoList from ImageDesc

    returns InfoList
    inDesc  = Python Obit input ImageDesc
    """
    ################################################################
    # Checks
    if not PIsA(inDesc):
        raise TypeError,"inDesc MUST be a Python Obit ImageDesc"
    #
    out    = InfoList.InfoList()
    out.me = Obit.InfoListUnref(out.me)
    out.me = Obit.ImageDescGetList(inDesc.me)
    return out
    # end PGetList 

def PIsA (inID):
    """ Tells if the input really is a Python Obit ImageDesc

    returns true or false (1,0)
    inID = Python ImageDesc to test
    """
    ################################################################
     # Checks
    if inID.__class__ != ImageDesc:
        return 0
    #
    return Obit.ImageDescIsA(inID.me)
    # end  PIsA

def PUnref (inID):
    """ Decrement reference count

    Decrement reference count which will destroy object if it goes to zero
    Python object stays defined.
    inID   = Python ImageDesc object
    """
    ################################################################
     # Checks
    if not PIsA(inID):
        raise TypeError,"inID MUST be a Python Obit ImageDesc"

    inID.me = Obit.ImageDescUnref(inID.me)
    # end PUnref

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
