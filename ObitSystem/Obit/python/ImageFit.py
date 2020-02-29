""" Python Obit ImageFit class

This class enables fitting models to images

ImageFit Members with python interfaces:

* List - used to pass instructions to processing 
"""
# $Id$
#-----------------------------------------------------------------------
#  Copyright (C) 2007,2019
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

# Obit ImageFit
from __future__ import absolute_import
from __future__ import print_function
import Obit, _Obit, OErr, Image, ImageDesc, InfoList, FitRegion, FitModel

# Python shadow class to ObitImageFit class

# class name in C
myClass = "ObitImageFit"
 
class ImageFit(Obit.ImageFit):
    """
    Python Obit ImageFit class
    
    This class enables fitting models to images.
    Fits models defined in a FitRegion to an image
    
    ImageFit Members with python interfaces:

    * List - used to pass instructions to processing 
    """
    def __init__(self, name) :
        super(ImageFit, self).__init__()
        Obit.CreateImageFit(self.this, name)
        self.myClass = myClass
    def __del__(self, DeleteImageFit=_Obit.DeleteImageFit):
        if _Obit!=None:
            DeleteImageFit(self.this)
    def __setattr__(self,name,value):
        if name == "me" :
            # Out with the old
            if self.this!=None:
                Obit.ImageFitUnref(Obit.ImageFit_Get_me(self.this))
            # In with the new
            Obit.ImageFit_Set_me(self.this,value)
            return
        if name=="List":
            PSetList(self,value)
            return
        self.__dict__[name] = value
    def __getattr__(self,name):
        if not isinstance(self, ImageFit):
            return "Bogus dude "+str(self.__class__)
        if name == "me" : 
            return Obit.ImageFit_Get_me(self.this)
        # Virtual members
        if name=="List":
            return PGetList(self)
        raise AttributeError(str(name))
    def __repr__(self):
        if not isinstance(self, ImageFit):
            return "Bogus dude "+str(self.__class__)
        return "<C ImageFit instance> " + Obit.ImageFitGetName(self.me)
    def cast(self, toClass):
        """ Casts object pointer to specified class

        This doesn't seem necessary
        * self     = object whose cast pointer is desired
        * toClass  = Class string to cast to
        """
        ################################################################
        # Get pointer with type of this class
        out =  self.me
        out = out.replace(self.myClass, toClass)
        return out
    # end cast
    
    # Input structure for model fitting
    cFitInput={'structure':['Fit',[
        ('fitImage',  'Image to be fitted'),
        ('fitRegion', 'FitRegion to be fitted'),
        ('MaxIter',   'Maximum number of iterations [def. 10 per fitted parameter]'),
        ('prtLv',     'Message level, 0=>none [def 0]'),
        ('PosGuard',  'Distance (cells) from edge to allow center  [def no bound]'),
        ('FluxLow',   'Lower bounds on Flux density [def no bound]'),
        ('FixFlux',   'Fix fluxes to input values [def False]'),
        ('FixPos',    'Fix Positions to input values [def False]'),
        ('GMajUp',    'Major axis upper bound (cells) [def no bound]'),
        ('GMajLow',   'Major axis lower bound (cells) [def no bound]'),
        ('GMinUp',    'Minor axis upper bound (cells) [def no bound]'),
        ('GMinLow',   'Minor axis lower bound (cells) [def no bound]')]],
               # defaults
               'fitImage':None,
               'fitRegion':None,
               'MaxIter':0,
               'prtLv':0,
               'PosGuard':0.0,
               'FluxLow':0.0,
               'FixFlux':False,
               'FixPos':False,
               'GMajUp':1.0e20,
               'GMajLow':0.0,
               'GMinUp':1.0e20,
               'GMinLow':0.0}
    
    def Fit (self, err, input=cFitInput):
        """
        Fit a model to an image
        
        Resultant model left in FitRegion reg

        * inImageFit = Python ImageFit object
        * image      = ObitImage to be fitted
        * reg        = Fit region defining what is to be fitted and initial guess
        * err        = Python Obit Error/message stack
        * input      = input parameter dictionary
        
        Input dictionary entries:

        ========= ================================================================
        fitImage  Image to be fitted
        fitRegion FitRegion to be fitted
        MaxIter   int Maximum number of iterations [def. 10 per fitted parameter]
        prtLv     int Message level, 0=>none [def 0]
        PosGuard  float Distance (cells) from edge to allow center  [def no bound]
        FluxLow   float Lower bounds on Flux density [def no bound]
        FixFlux   bool  Fix flux densities? [def False]
        FixPos    bool  Fix Positions? [def False]
        GMajUp    float Major axis upper bound (cells) [def no bound]
        GMajLow   float Major axis lower bound (cells) [def no bound]
        GMinUp    float Minor axis upper bound (cells) [def no bound]
        GMinLow   float Minor axis lower bound (cells) [def no bound]
        ========= ================================================================
        """

        PFit(self, err, input=input)
        # end Fit
    
    # end class ImageFit

# Allow external access to inputs
FitInput = ImageFit.cFitInput

def input(inputDict):
    """
    Print the contents of an input Dictionary

    * inputDict = Python Dictionary containing the parameters for a routine

    There should be a member of the dictionary ('structure') with a value
    being a list containing:

    1) The name for which the input is intended (string)
    2) a list of tuples consisting of (parameter name, doc string)
       with an entry for each parameter in the dictionary.

    The display of the the inputs dictionary will be in the order of
    the tuples and display the doc string after the value.
    
    An example:
   
    >>> Soln2CalInput={'structure':['Soln2Cal',[('InData','Input OTF'),
                                                ('soln','input soln table version'),
                                                ('oldCal','input cal table version, -1=none'),
                                                ('newCal','output cal table')]],
                       'InData':None, 'soln':0, 'oldCal':-1, 'newCal':0}
    """
    ################################################################
    structure = inputDict['structure']  # Structure information
    print('Inputs for ',structure[0])
    for k,v in structure[1]:
        print('  ',k,' = ',inputDict[k],' : ',v)
        
    # end input

def PCopy (inImageFit, outImageFit, err):
    """
    Make a shallow copy of input object.
    
    Makes structure the same as inImageFit, copies pointers

    * inImageFit  = Python ImageFit object to copy
    * outImageFit = Output Python ImageFit object, must be defined
    * err         = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not PIsA(inImageFit):
        raise TypeError("inImageFit MUST be a Python Obit ImageFit")
    if not PIsA(outImageFit):
        raise TypeError("outImageFit MUST be a Python Obit ImageFit")
    if not OErr.OErrIsA(err):
        raise TypeError("err MUST be an OErr")
    #
    #smi = inImageFit.cast(myClass)   # cast pointer
    #smo = outImageFit.cast(myClass)  # cast pointer
    Obit.ImageFitCopy (inImageFit.me, outImageFit.me, err.me)
    if err.isErr:
        OErr.printErrMsg(err, "Error copying ImageFit")
    # end PCopy

def PGetList (inImageFit):
    """
    Return the member InfoList
    
    returns InfoList

    * inImageFit  = Python ImageFit object
    """
    ################################################################
     # Checks
    if not PIsA(inImageFit):
        raise TypeError("inImageFit MUST be a Python Obit ImageFit")
    #
    #sm = inImageFit.cast(myClass)  # cast pointer
    out    = InfoList.InfoList()
    out.me = Obit.ImageFitGetList(inImageFit.me)
    return out
    # end PGetList

def PCreate (name):
    """
    Create the parameters and underlying structures of a ImageFit.

    * name      = Name to be given to object
                Most control parameters are in InfoList member
    """
    ################################################################
    #
    out = ImageFit(name);
    return out;
    # end PCreate

def PFit (inImageFit, err, input=FitInput):
    """
    Fit a model to an image
    
    Resultant model left in FitRegion reg

    * inImageFit = Python ImageFit object
    * image      = ObitImage to be fitted
    * reg        = Fit region defining what is to be fitted and initial guess
    * err        = Python Obit Error/message stack
    * input      = input parameter dictionary
    
    Input dictionary entries:

    ========= ================================================================
    fitImage  Image to be fitted
    fitRegion FitRegion to be fitted
    MaxIter   int Maximum number of iterations [def. 10 per fitted parameter]
    prtLv     int Message level, 0=>none [def 0]
    PosGuard  float Distance (cells) from edge to allow center  [def no bound]
    FluxLow   float Lower bounds on Flux density [def no bound]
    FixFlux   bool  Fix flux densities? [def False]
    FixPos    bool  Fix Positions? [def False]
    GMajUp    float Major axis upper bound (cells) [def no bound]
    GMajLow   float Major axis lower bound (cells) [def no bound]
    GMinUp    float Minor axis upper bound (cells) [def no bound]
    GMinLow   float Minor axis lower bound (cells) [def no bound]
    ========= ================================================================
    """
    ################################################################
    # Get input parameters
    fitImage   = input["fitImage"]
    fitRegion  = input["fitRegion"]
    # Checks
    if not PIsA(inImageFit):
        raise TypeError("inImageFit MUST be a Python Obit ImageFit")
    if not Image.PIsA(fitImage):
        raise TypeError("fitImage MUST be a Python Obit Image")
    if not FitRegion.PIsA(fitRegion):
        raise TypeError("fitRegion MUST be a Python Obit FitRegion")
    #
    dim = [1,1,1,1,1]
    #
    # Set control values on ImageFit
    inInfo = PGetList(inImageFit)    # 
    InfoList.PAlwaysPutInt    (inInfo, "MaxIter",  dim, [input["MaxIter"]])
    InfoList.PAlwaysPutInt   (inInfo, "prtLv",    dim, [input["prtLv"]])
    InfoList.PAlwaysPutBoolean (inInfo, "FixPos",  dim, [input["FixPos"]])
    InfoList.PAlwaysPutBoolean (inInfo, "FixFlux", dim, [input["FixFlux"]])
    InfoList.PAlwaysPutDouble (inInfo, "PosGuard", dim, [input["PosGuard"]])
    InfoList.PAlwaysPutDouble (inInfo, "FluxLow",  dim, [input["FluxLow"]])
    InfoList.PAlwaysPutDouble (inInfo, "GMajUp",   dim, [input["GMajUp"]])
    InfoList.PAlwaysPutDouble (inInfo, "GMajLow",  dim, [input["GMajLow"]])
    InfoList.PAlwaysPutDouble (inInfo, "GMinUp",   dim, [input["GMinUp"]])
    InfoList.PAlwaysPutDouble (inInfo, "GMinLow",  dim, [input["GMinLow"]])
    #
    # Do operation
    #sm = inImageFit.cast(myClass)  # cast pointer
    ret = Obit.ImageFitFit(inImageFit.me, fitImage.me, fitRegion.me, err.me)
    if ret==0:
        print("Fit converged")
    else:
        print("Fit hit iteration limit")
    if err.isErr:
        OErr.printErrMsg(err, "Error Fitting model")
    # end PFit

def PGetName (inImageFit):
    """
    Tells Image object name (label)
    
    returns name as character string

    * inImageFit  = Python ImageFit object
    """
    ################################################################
     # Checks
    if not PIsA(inImageFit):
        raise TypeError("inImageFit MUST be a Python Obit ImageFit")
    #
    #sm = inImageFit.cast(myClass)  # cast pointer
    return Obit.ImageFitGetName(inImageFit.me)
    # end PGetName

def PIsA (inImageFit):
    """
    Tells if input really a Python Obit ImageFit
    
    return True, False
    * inImageFit   = Python ImageFit object
    """
    ################################################################
    # Checks - allow inheritence
    if not isinstance(inImageFit,ImageFit):
        return False
    return Obit.ImageFitIsA(inImageFit.me)!=0
    # end PIsA
