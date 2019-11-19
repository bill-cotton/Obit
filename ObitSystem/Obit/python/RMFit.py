""" Python Obit RMFit class

Class for fitting rotation measures to image pixels
This class does least squares fitting of rotation measure (RM in rad m^2) 
and EVPA at the reference lambda^2 in radians.
Either an image cube or a set of single plane images at arbitrary 
frequencies may be fitted.
The result is an image cube with planes, RM, EVPA0 as the planes, 
if doErr requested, the statistical errors of these planes
and the Chi Squared is also returned.
"""
# $Id$
#-----------------------------------------------------------------------
#  Copyright (C) 2013-2019
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

# Obit RMFit
from __future__ import absolute_import
from __future__ import print_function
import Obit, _Obit, OErr, Image, InfoList, BeamShape

# Python shadow class to ObitRMFit class

# class name in C
myClass = "ObitRMFit"
 
#
class RMFit(Obit.RMFit):
    """ Python Obit RMFit class
    
    Class for fitting spectra to image pixels
    
    RMFit Members with python interfaces:
        List      - used to pass instructions to processing
    """
    def __init__(self, name="None", nterm=2) :
        super(RMFit, self).__init__()
        Obit.CreateRMFit (self.this, name, nterm)
        self.myClass = myClass
    def __del__(self, DeleteRMFit=_Obit.DeleteRMFit):
        if _Obit!=None:
            DeleteRMFit(self.this)
    def __setattr__(self,name,value):
        if name == "me" :
            # Out with the old
            if self.this!=None:
                Obit.RMFitUnref(Obit.RMFit_Get_me(self.this))
            # In with the new
            Obit.RMFit_Set_me(self.this,value)
            return
        self.__dict__[name] = value
    def __getattr__(self,name):
        if not isinstance(self, RMFit):
            return  "Bogus dude"+str(self.__class__)
        if name == "me" : 
            return Obit.RMFit_Get_me(self.this)
        if name=="List":
            if not PIsA(self):
                raise TypeError("input MUST be a Python Obit RMFit")
            out    = InfoList.InfoList()
            out.me = Obit.RMFitGetList(self.cast(myClass))
            return out
        raise AttributeError(name)
    def __repr__(self):
        if not isinstance(self, RMFit):
            return  "Bogus dude"+str(self.__class__)
        return "<C RMFit instance> " + Obit.RMFitGetName(self.me)
    
    def Cube (self, inQImage, inUImage, outImage, err):
        """ Fit RMs to an image cube
        
        Fitted rotation measure parameters returned in outImage
        There are two options, if doRMSyn then a straight RM synthesis is done
        resulting in the maximum value.  Better for low SNR data.
        Otherwise a coarse RM synthesis is done followed by a nonlinear least
        squares fit.  Best for high SNR data.
        Can run with multiple threads if enabled:
        OSystem.PAllowThreads(2)  # 2 threads
        self     = RMFit object, parameters on List:
            refLamb2  double scalar Reference lambda^2 for fit [def ref for inQImage]
            minQUSNR  float min. SNR for Q and U pixels [def 3.0]
            minFrac   float min fraction of samples included  [def 0.5]
            doError   boolean scalar If true do error analysis [def False]
            doRMSyn   boolean scalar If true do max RM synthesis [def False]
            maxRMSyn  float scalar max RM to search (rad/m^2) [def ambiguity]
            minRMSyn  float scalar min RM to search (rad/m^2)[def -ambiguity]
                      values outside of [minRMSyn,maxRMSyn] ignored
            delRMSyn  float scalar RM increment for search (rad/m^2)[def 1.0]
            maxChi2   float scalar Max. chi^2 for doRMSyn [def 10]
        inQImage = Q Image cube to be fitted
        inUImage = U Image cube to be fitted, must have same geometry as inQImage
        outImage = Image cube to accept fitted RM, EVPA0.
                   Should be defined but not created.
                   If doRMSyn: plane 1=RM (rad/m^2), 2=EVPA(rad), 3=amp (Jy), 4=Chi^2
                   else
                   Planes 1->nterm are coefficients per pixel (rad m^2, rad)
                   Planes nterm+1->2*nterm are uncertainties in coefficients
                   Plane 2*nterm+1 = Chi squared of fit
        err      = Obit error stack
        """
        ################################################################
        # Checks
        if not PIsA(self):
            raise TypeError("self MUST be a Python Obit RMFit")
        if not inQImage.ImageIsA():
            raise TypeError("inQImage MUST be a Python Obit Image")
        if not inUImage.ImageIsA():
            raise TypeError("inUImage MUST be a Python Obit Image")
        if not outImage.ImageIsA():
            raise TypeError("outImage MUST be a Python Obit Image")
        #
        Obit.RMFitCube(self.me, inQImage.me, inUImage.me,  outImage.me, err.me)
    # end Cube

    def ImArr (self, imQArr, imUArr, outImage, err):
        """ Fit rotation measures to an array of images
        
        Fitted RM parameters returned in outImage
        Can run with multiple threads if enabled:
        OSystem.PAllowThreads(2)  # 2 threads
        self     = RMFit object, parameters on List:
            refLamb2  double scalar Reference lambda^2 for fit [def ref for inQImage]
            minQUSNR  float min. SNR for Q and U pixels [def 3.0]
            minFrac   float min fraction of samples included [def 0.5]
            doError   boolean scalar If true do error analysis [def False]
        imQArr   = Array of 2D images to be fitted
        imUArr   = Array of 2D images to be fitted, same geometries as imQArr
        outImage = Image cube with fitted spectra.
                   Should be defined but not created.
                   Planes 1->nterm are coefficients per pixel
                   Planes nterm+1->2*nterm are uncertainties in coefficients
                   Plane 2*nterm+1 = Chi squared of fit
        err      = Obit error stack
        """
        ################################################################
        # Checks
        if not PIsA(self):
            raise TypeError("self MUST be a Python Obit RMFit")
        if not imArr[0].ImageIsA():
            raise TypeError("imArr[0] MUST be a Python Obit Image")
        if not outImage.ImageIsA():
            raise TypeError("outImage MUST be a Python Obit Image")
        #
        nimage = len(imArr)
        imQArrMe = []
        imUArrMe = []
        for x in imQArr:
            imQArrMe.append(x.me)
        for x in imUArr:
            imUArrMe.append(x.me)

        Obit.RMFitImArr(self.me, nimage, imQArrMe, imUArrMe, outImage.me, err.me)
    # end ImArr
    
def PIsA (inRMFit):
    """ Tells if input really a Python Obit RMFit

    return True, False
    inRMFit   = Python RMFit object
    """
    ################################################################
    if not isinstance(inRMFit, RMFit):
        print("Actually",inRMFit.__class__)
        return False
    # Checks - allow inheritence
    return Obit.RMFitIsA(inRMFit.me)!=0
    # end PIsA

def PCreate (name, nterm=2):
    """ Create the underlying structures of a BeamShape

    return object created.
    name      = Name to be given to object
    nterm     = Number of coefficients, RM, EVPA0 (1 or 2)
    """
    ################################################################
    #
    out = RMFit("None");
    out.me = Obit.RMFitCreate(name, nterm)
    return out;
    # end PCreate

def PSingle (refLamb2, lamb2, qflux, qsigma, uflux, usigma, err, nterm=2):
    """  Fit RM, EVPA0 to Q, U flux measurements
    
    Also does error analysis
    Returns  array of fitter parameters, errors for each and Chi Squares of fit
    refLamb2 = Reference lambda^2 for fit (m^2)
    lamb2    = Array of lambda^2 for fit (m^2)
    qflux    = Array of Q fluxes (Jy) same dim as lamb2
    qsigma   = Array of Q errors (Jy) same dim as lamb2
    uflux    = Array of U fluxes (Jy) same dim as lamb2
    usigma   = Array of U errors (Jy) same dim as lamb2
    err      = Obit error stack
    nterm    = Number of coefficients to fit (1 or 2)
    """
    ################################################################
    #
    nlamb2 = len(lamb2)
    ret = Obit.RMFitSingle(nlamb2, nterm, refLamb2, lamb2, 
                           qflux, qsigma, uflux, usigma, err.me)
    OErr.printErr(err)
    OErr.printErrMsg(err,"Fitting failed")
    return ret
    # end PSingle
