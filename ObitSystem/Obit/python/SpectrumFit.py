""" Python Obit SpectrumFit class

Class for fitting spectra to image pixels
This class does least squares fitting of log(s) as a polynomial in log($\nu$).
Either an image cube or a set of single plane images at arbitrary 
frequencies may be fitted.
Model: S = S_0 exp (alpha*ln(nu/nu_0) + beta*ln(nu/nu_0)**2+...)
The result is an image cube with planes, S_0, alpha, beta,...
as the planes, if doErr requested, the statistical errors of these planes
and the Chi Squared is also returned.
The function ObitSpectrumFitEval will evaluate this fit and return an image
with the flux densities at the desired frequencies.
"""
# $Id$
#-----------------------------------------------------------------------
#  Copyright (C) 2008-2019
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

# Obit SpectrumFit
from __future__ import absolute_import
from __future__ import print_function
import Obit, _Obit, OErr, Image, InfoList

# Python shadow class to ObitSpectrumFit class

# class name in C
myClass = "ObitSpectrumFit"
 

class SpectrumFit(Obit.SpectrumFit):
    """ Python Obit SpectrumFit class
    
    Class for fitting spectra to image pixels
    
    SpectrumFit Members with python interfaces:
        List      - used to pass instructions to processing
    """
    def __init__(self, name="None", nterm=2) :
        super(SpectrumFit, self).__init__()
        Obit.CreateSpectrumFit (self.this, name, nterm)
        self.myClass = myClass
    def __del__(self, DeleteSpectrumFit=_Obit.DeleteSpectrumFit):
        if _Obit!=None:
            DeleteSpectrumFit(self.this)
    def __setattr__(self,name,value):
        if name == "me" :
            # Out with the old
            if self.this!=None:
                Obit.SpectrumFitUnref(Obit.SpectrumFit_Get_me(self.this))
            # In with the new
            Obit.SpectrumFit_Set_me(self.this,value)
            return
        self.__dict__[name] = value
    def __getattr__(self,name):
        if not isinstance(self, SpectrumFit):
            return"Bogus dude"+str(self.__class__)
        if name == "me" : 
            return Obit.SpectrumFit_Get_me(self.this)
        if name=="List":
            if not PIsA(self):
                raise TypeError("input MUST be a Python Obit SpectrumFit")
            out    = InfoList.InfoList()
            out.me = Obit.SpectrumFitGetList(self.me)
            return out
        raise AttributeError(name)
    def __repr__(self):
        if not isinstance(self, SpectrumFit):
            return"Bogus dude"+str(self.__class__)
        return "<C SpectrumFit instance> " + Obit.SpectrumFitGetName(self.me)
    
    def Cube (self, inImage, outImage, err):
        """ Fit spectra to an image cube

        Fitted spectral polynomials returned in outImage
        Can run with multiple threads if enabled:
        OSystem.PAllowThreads(2)  # 2 threads
        Note: for images of the type produced by MFImage (ImageMF class) use
        MFImage.PFitSpec.
        self     = SpectrumFit object, parameters on List:
            refFreq   double scalar Reference frequency for fit [def ref for inImage]
            maxChi2   float scalar Max. Chi Sq for accepting a partial spectrum [def 2.0]
            doError   boolean scalar If true do error analysis [def False]
            doPBCor   boolean scalar If true do primary beam correction. [def False]
            doBrokePow boolean scalar If true do broken power law (3 terms). [def False]
            calFract  float (?,1,1) Calibration error as fraction of flux
                      One per frequency or one for all, def 0.05
            PBmin     float (?,1,1) Minimum beam gain correction
                      One per frequency or one for all, def 0.05,
                      1.0 => no gain corrections
            antSize   float (?,1,1) Antenna diameter (m) for gain corr, 
                      One per frequency or one for all, def 25.0
            corAlpha  float scalar spectral index correction to apply, def = 0.0
        inImage  = Image cube to be fitted
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
            raise TypeError("self MUST be a Python Obit SpectrumFit")
        if not inImage.ImageIsA():
            raise TypeError("inImage MUST be a Python Obit Image")
        if not outImage.ImageIsA():
            raise TypeError("outImage MUST be a Python Obit Image")
        #
        Obit.SpectrumFitCube(self.me, inImage.me, outImage.me, err.me)
    # end Cube

    def ImArr (self, imArr, outImage, err):
        """ Fit spectra to an array of images
        
        Fitted spectral polynomials returned in outImage
        Can run with multiple threads if enabled:
        OSystem.PAllowThreads(2)  # 2 threads
        NB: there is a maximum of 20 images in imArr
        self     = SpectrumFit object, parameters on List:
            refFreq   double scalar Reference frequency for fit [def average of inputs]
            maxChi2   float scalar Max. Chi Sq for accepting a partial spectrum [def 2.0]
            doError   boolean scalar If true do error analysis [def False]
            doPBCor   boolean scalar If true do primary beam correction. [def False]
            doBrokePow boolean scalar If true do broken power law (3 terms). [def False]
            calFract  float (?,1,1) Calibration error as fraction of flux
                      One per frequency or one for all, def 0.05
            PBmin     float (?,1,1) Minimum beam gain correction
                      One per frequency or one for all, def 0.05,
                      1.0 => no gain corrections
            antSize   float (?,1,1) Antenna diameter (m) for gain corr, 
                      One per frequency or one for all, def 25.0
            corAlpha  float scalar spectral index correction to apply, def = 0.0
        imArr    = Array of 2D images to be fitted
                   NB: there is a maximum of 20 images in imArr
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
            raise TypeError("self MUST be a Python Obit SpectrumFit")
        if not imArr[0].ImageIsA():
            raise TypeError("imArr[0] MUST be a Python Obit Image")
        if not outImage.ImageIsA():
            raise TypeError("outImage MUST be a Python Obit Image")
        nimage = len(imArr)
        if nimage>20:
            raise TypeError("Too many input images > 20")
            
        #
        # SWIG SUCKS
        if nimage>0:
            im1 = imArr[0].me
        else:
            im1 = None
        if nimage>1:
            im2 = imArr[1].me
        else:
            im2 = None
        if nimage>2:
            im3 = imArr[2].me
        else:
            im3 = None
        if nimage>3:
            im4 = imArr[3].me
        else:
            im4 = None
        if nimage>4:
            im5 = imArr[4].me
        else:
            im5 = None
        if nimage>5:
            im6 = imArr[5].me
        else:
            im6 = None
        if nimage>6:
            im7 = imArr[6].me
        else:
            im7 = None
        if nimage>7:
            im8 = imArr[7].me
        else:
            im8 = None
        if nimage>8:
            im9 = imArr[8].me
        else:
            im9 = None
        if nimage>9:
            im10 = imArr[9].me
        else:
            im10 = None
        if nimage>10:
            im11 = imArr[10].me
        else:
            im11 = None
        if nimage>11:
            im12 = imArr[11].me
        else:
            im12 = None
        if nimage>12:
            im13 = imArr[12].me
        else:
            im13 = None
        if nimage>13:
            im14 = imArr[13].me
        else:
            im14 = None
        if nimage>14:
            im15 = imArr[14].me
        else:
            im15 = None
        if nimage>15:
            im16 = imArr[15].me
        else:
            im16 = None
        if nimage>16:
            im17 = imArr[16].me
        else:
            im17 = None
        if nimage>17:
            im18 = imArr[17].me
        else:
            im18 = None
        if nimage>18:
            im19 = imArr[18].me
        else:
            im19 = None
        if nimage>19:
            im20 = imArr[19].me
        else:
            im20 = None
        Obit.SpectrumFitImArr(self.me, nimage, 
                              im1, im2, im3, im4, im5, im6, im7, im8,
                              im9, im10, im11, im12, im13, im14, im15,
                              im16, im17, im18, im19, im20,
                              outImage.me, err.me)
    # end ImArr
    
    def Eval (self, inImage, outFreq, outImage, err):
        """ Evaluate spectrum at a given frequency
        
        Returns image at specified frequency
        self     = SpectrumFit object
        inImage  = Spectral coefficient image
                   Planes 1->nterm are coefficients per pixel
                   Planes nterm+1->2*nterm are uncertainties in coefficients
                   Plane 2*nterm+1 = Chi squared of fit
        outFreq  = Output Frequency in Hz
        outImage = Image to write, must be defined but not yet exist.
        err      = Obit error stack
        """
        ################################################################
        # Checks
        if not PIsA(self):
            raise TypeError("self MUST be a Python Obit SpectrumFit")
        if not inImage.ImageIsA():
            raise TypeError("inImage MUST be a Python Obit Image")
        if not outImage.ImageIsA():
            raise TypeError("outImage MUST be a Python Obit Image")
        #
        Obit.SpectrumFitEval(self.me, inImage.me, outFreq, outImage.me, err.me)
    # end Eval
    
    # end class SpectrumFit
    
def PIsA (inSpectrumFit):
    """ Tells if input really a Python Obit SpectrumFit

    return True, False 
    inSpectrumFit   = Python SpectrumFit object
    """
    ################################################################
    if not isinstance(inSpectrumFit, SpectrumFit):
        print("Actually",inSpectrumFit.__class__)
        return False
    # Checks
    return Obit.SpectrumFitIsA(inSpectrumFit.me)!=0
    # end PIsA

def PCreate (name, nterm):
    """ Create the underlying structures of a Spectrum fitter

    return object created.
    name      = Name to be given to object
    nterm     = Number of coefficients of powers of log10(nu)
    """
    ################################################################
    #
    out = SpectrumFit("None");
    out.me = Obit.SpectrumFitCreate(name, nterm)
    return out;
    # end PCreate

def PSingle (nterm, refFreq, freq, flux, sigma, err, doBrokePow=False):
    """  Fit single spectrum to flux measurements
    
    Does error analysis and makes primary beam correction
    Returns  array of fitter parameters, errors for each and Chi Squares of fit
             Initial terms are in Jy, other in log.
    nterm   = Number of coefficients of powers of log(nu) to fit
    refFreq = Reference frequency for fit (Hz)
    freq    = Array of Frequencies (Hz)
    flux    = Array of fluxes (Jy) same dim as freq
    sigma   = Array of errors (Jy) same dim as freq
    err     = Obit error stack
    doBrokePow = If true do broken power law (3 terms)
    """
    ################################################################
    #
    nfreq = len(freq)
    ret = Obit.SpectrumFitSingle(nfreq, nterm, refFreq, freq, flux, sigma, \
                                 doBrokePow, err.me)
    OErr.printErr(err)
    OErr.printErrMsg(err,"Fitting failed")
    return ret
    # end PSingle
