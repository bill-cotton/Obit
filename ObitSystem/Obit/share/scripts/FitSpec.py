# Fit spectra to a set of images where flux density>minFlux
# Using weighting from AvgWeights.pickle
#exec(open('FitSpec.py').read())
minFlux   = 0.0002 # Minimum flux density for SI
defaultSI = 0.0    # Default spectral index
doPBCor   = False  # Do PB correct, not for Mosaics
antSize   = 24.5   # Antenna size for PB Corr
PBmin     = 0.05   # Min PB gain for correction
doBlank   = False  # Blank total intensity and SI where highest plane blanked?
                   # After any PB Correction
nthreads  = 16     # Number of threads to use
dir       = None   # If given the data directory
                   # Else the value from AvgWeights.pickle
post      = None   # If given, the ending of file names,
                   # Else the value from AvgWeights.pickle
fitsfiles = None   # If given, list of (name root, broadband RMS)
                   # Else the value from AvgWeights.pickle
# Set non defaults here

import Image, OErr, OSystem, FArray
OSystem.PAllowThreads(nthreads)  # with threads
from PipeUtil import SaveObject, FetchObject
import math, pickle

# Get weighting info
stuff = FetchObject("AvgWeights.pickle")
refFreq   = stuff['EffFreq'] # Effective frequency in MHz
wwts      = stuff['weights'] # Subband weights in %
if dir==None:
    dir       = stuff['dir']     # Data directory
if fitsfiles==None:
    fitsfiles = stuff['RMSes']   # List of (file_name root, RMS)
if post==None:
    post = stuff['post']    # Ending of file name

import Obit, ImageMF, Image, FArray, OSystem, OErr
err = OErr.OErr()
OSystem.PAllowThreads(nthreads)  # with threads

for fitsfile in fitsfiles:
    print ("Doing",fitsfile[0])
    inMF  = ImageMF.newPFImageMF('Input',dir+fitsfile[0]+post, 0, True, err)
    ImageMF.PFitSpec(inMF, err, antSize=0., nOrder=1,corAlpha=defaultSI,
                     refFreq=refFreq*1.0e6, Weights=wwts, minFlux=minFlux, 
                     doPBCor=doPBCor, PBmin=PBmin)
    if doBlank:
        print ('Blanking by high freq. bin')
        # Blanking mask
        outIm = Image.newPFImage('Output',dir+fitsfile[0]+post, 0, True, err)
        hiPlane = [outIm.Desc.Dict['inaxes'][2],1,1,1,1]
        if doPBCor:
            # Mask by primary beam correction and high bin
            PBImg = Image.Image('PBImage'); Image.PCloneMem(outIm,PBImg,err)
            # Blank where highest plane blanked
            ImageUtil.PPBImage(outIm,PBImg,err,PBmin,outPlane=hiPlane,antSize=antSize)
            PBImg.GetPlane(None,hiPlane, err); mask1 = PBImg.FArray
            outIm.GetPlane(None,hiPlane, err); FArray.PBlank(mask1, outIm.FArray, mask1)
        else:
            # Mask by highest frequency
            outIm.GetPlane(None,hiPlane, err); mask1 = FArray.PCopy(outIm.FArray, err)
            PBImg = None
        # Plane 1 - total intensity
        outIm.GetPlane(None, [1,1,1,1,1], err)
        FArray.PBlank(outIm.FArray, mask1, outIm.FArray)
        outIm.PutPlane(None, [1,1,1,1,1], err)
        # Plane 2 - Spectral index
        outIm.GetPlane(None, [2,1,1,1,1], err)
        FArray.PBlank(outIm.FArray, mask1, outIm.FArray)
        outIm.PutPlane(None, [2,1,1,1,1], err)
        #DEBUG Image.PFArray2FITS(mask1, "Mask.fits", err, outDisk=0)
        del PBImg, mask1
# End blanking
# end loop
