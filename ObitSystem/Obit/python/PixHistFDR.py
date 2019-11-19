"""
Python Obit PixHistFDR class

This class handles image pixel histograms and calculation of
false detection rates (FDR) of sources.

PixHistFDR Members with python interfaces:
=======  =================================================================
List      used to pass instructions to processing
Sigma     pixel RMS in selected region (from PPixHistFDRHisto)
ImPix     Image pixels FArray
DifHist   Differential histogram FArray
IntHist   Integral histogram FArray
=======  =================================================================
"""
# $Id$
#-----------------------------------------------------------------------
#  Copyright (C) 2011,2019
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

# Obit PixHistFDR
from __future__ import absolute_import
import Obit, _Obit, OErr, Image, FArray, InfoList

# Python shadow class to ObitPixHistFDR class

# class name in C
myClass = "ObitPixHistFDR"
 
#
class PixHistFDR(Obit.PixHistFDR):
    """ Python Obit PixHistFDR class
    
    This class enables fitting models to images.
    Fits models defined in a FitRegion to an image
    
    PixHistFDR Members with python interfaces:
    List      - used to pass instructions to processing 
    """
    def __init__(self, name) :
        super(PixHistFDR, self).__init__()
        Obit.CreatePixHistFDR(self.this, name)
        self.myClass = myClass
    def __del__(self, DeletePixHistFDR=_Obit.DeletePixHistFDR):
        if _Obit!=None:
            DeletePixHistFDR(self.this)
    def __setattr__(self,name,value):
        if name == "me" :
            # Out with the old
            if self.this!=None:
                Obit.PixHistFDRUnref(Obit.PixHistFDR_Get_me(self.this))
            # In with the new
            Obit.PixHistFDR_Set_me(self.this,value)
            return
        if name=="List":
            PSetList(self,value)
            return
        self.__dict__[name] = value
    def __getattr__(self,name):
        if not isinstance(self, PixHistFDR):
            return "Bogus dude",str(self.__class__)
        if name == "me" : 
            return Obit.PixHistFDR_Get_me(self.this)
        # Virtual members
        if name=="List":
            return PGetList(self)
        if name=="Sigma":
            return PSigma(self)
        if name=="ImPix":
            return PImPix(self)
        if name=="DifHist":
            return PDifHist(self)
        if name=="IntHist":
            return PIntHist(self)
        raise AttributeError(str(name))
    def __repr__(self):
        if not isinstance(self, PixHistFDR):
            return "Bogus dude",str(self.__class__)
        return "<C PixHistFDR instance> " + Obit.PixHistFDRGetName(self.me)
    def Histo (self, blc, trc, nbin, range, err):
        """ Computes histograms
        
        self     = object 
        blc      = BLC (1-rel) in image
        trc      = TRC (1-rel) in image
        nbin     = number of bins (should be odd)
        range    = range in units of Image RMS about 0
        err      = Obit error/message stack
        """
        ################################################################
        # Checks
        if not PIsA(self):
            raise TypeError("self MUST be a Python Obit PixHistFDR")
        # make sure blc, trc are long
        lblc=[]; ltrc=[]
        for t in trc:
            ltrc.append(int(t))
        for t in blc:
            lblc.append(int(t))
        Obit.PixHistFDRHisto (self.me, lblc, ltrc, nbin, range, err.me)
        if err.isErr:
            OErr.printErrMsg(err, "Error creating PixHistFDR histogram")
        return
    # end Histo
    
    def Flux (self, maxFDR, err):
        """ Computes flux density at target FDR
        
        returns flux density.
        self     = object 
        maxFDR   = target False detection rate as a fraction
        err      = Obit error/message stack
        """
        ################################################################
        # Checks
        if not PIsA(self):
            raise TypeError("self MUST be a Python Obit PixHistFDR")
        ret =  Obit.PixHistFDRFlux (self.me, maxFDR, err.me)
        if err.isErr:
            OErr.printErrMsg(err, "Error determining PixHistFDR flux")
        return ret
    # end Flux    
    # end class PixHistFDR

def PCreate (name, image, err):
    """ Create the underlying structures of a PixHistFDR

    return object created.
    name      = Name to be given to object
    image     = Obit python Image to use
    err       = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not Image.PIsA(image):
            raise TypeError("image MUST be a Python Obit Image")
    #
    out = PixHistFDR("None");
    out.me = Obit.PixHistFDRCreate(name, image.me, err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error creating PixHistFDR")
    return out;
    # end PCreate

def PCopy (inPixHistFDR, outPixHistFDR, err):
    """ Make a shallow copy of input object.

    Makes structure the same as inPixHistFDR, copies pointers
    inPixHistFDR  = Python PixHistFDR object to copy
    outPixHistFDR = Output Python PixHistFDR object, must be defined
    err         = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inPixHistFDR):
        raise TypeError("inPixHistFDR MUST be a Python Obit PixHistFDR")
    if not PIsA(outPixHistFDR):
        raise TypeError("outPixHistFDR MUST be a Python Obit PixHistFDR")
    if not OErr.OErrIsA(err):
        raise TypeError("err MUST be an OErr")
    #
    Obit.PixHistFDRCopy (inPixHistFDR.me, outPixHistFDR.me, err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error copying PixHistFDR")
    # end PCopy

def PGetList (inPixHistFDR):
    """ Return the member InfoList

    returns InfoList
    inPixHistFDR  = Python PixHistFDR object
    """
    ################################################################
     # Checks
    if not PIsA(inPixHistFDR):
        raise TypeError("inPixHistFDR MUST be a Python Obit PixHistFDR")
    #
    out    = InfoList.InfoList()
    out.me = Obit.PixHistFDRGetList(inPixHistFDR.me)
    return out
    # end PGetList

def PSigma (inPixHistFDR):
    """ Return the pixel RMS

    returns pixel RMS
    inPixHistFDR  = Python PixHistFDR object
    """
    ################################################################
     # Checks
    if not PIsA(inPixHistFDR):
        raise TypeError("inPixHistFDR MUST be a Python Obit PixHistFDR")
    #
    return Obit.PixHistFDRSigma(inPixHistFDR.me)
    # end PSigma

def PImPix (inPixHistFDR):
    """ Return the Image Pixel FArray

    returns Image Pixel FArray
    inPixHistFDR  = Python PixHistFDR object
    """
    ################################################################
     # Checks
    if not PIsA(inPixHistFDR):
        raise TypeError("inPixHistFDR MUST be a Python Obit PixHistFDR")
    #
    out    = FArray.FArray("ImPix")
    out.me = Obit.FArrayUnref(out.me)
    out.me = Obit.PixHistFDRGetImPix(inPixHistFDR.me)
    return out
    # end PImPix 

def PDifHist (inPixHistFDR):
    """ Return the differential pixel histogram

    returns differential pixel histogram
    inPixHistFDR  = Python PixHistFDR object
    """
    ################################################################
     # Checks
    if not PIsA(inPixHistFDR):
        raise TypeError("inPixHistFDR MUST be a Python Obit PixHistFDR")
    #
    out    = FArray.FArray("DifHist")
    out.me = Obit.FArrayUnref(out.me)
    out.me = Obit.PixHistFDRGetHisto(inPixHistFDR.me)
    return out
    # end PifHist

def PIntHist (inPixHistFDR):
    """ Return the integrated pixel histogram

    returns integrated pixel histogram
    inPixHistFDR  = Python PixHistFDR object
    """
    ################################################################
     # Checks
    if not PIsA(inPixHistFDR):
        raise TypeError("inPixHistFDR MUST be a Python Obit PixHistFDR")
    #
    out    = FArray.FArray("IntHist")
    out.me = Obit.FArrayUnref(out.me)
    out.me = Obit.PixHistFDRGetIntHisto(inPixHistFDR.me)
    return out
    # end PIntHist

def PGetName (inPixHistFDR):
    """ Tells Image object name (label)

    returns name as character string
    inPixHistFDR  = Python PixHistFDR object
    """
    ################################################################
     # Checks
    if not PIsA(inPixHistFDR):
        raise TypeError("inPixHistFDR MUST be a Python Obit PixHistFDR")
    #
    return Obit.PixHistFDRGetName(nPixHistFDR.me)
    # end PGetName

def PIsA (inPixHistFDR):
    """ Tells if input really a Python Obit PixHistFDR

    return True, False
    inPixHistFDR   = Python PixHistFDR object
    """
    ################################################################
    # Checks - allow inheritence
    if not str(inPixHistFDR.__class__).startswith("PixHistFDR"):
        return False
    return Obit.PixHistFDRIsA(inPixHistFDR)!=0
    # end PIsA
