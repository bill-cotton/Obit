# Fit spectra with errors to a set of images where flux density>minFlux
# Writes new 5 plane cube
# Using weighting from AvgWeights.pickle
#exec(open('FitError.py').read())
minFlux   = 0.0002 # Minimum flux density for SI
defaultSI = 0.0    # Default spectral index
doPBCor   = False  # Do PB correct, not for Mosaics
maxChi2   = 0.01   # fit acceptance level
calFract  = 0.00   # % relative calibration error for error estimates
antSize   = 24.5   # Antenna size for PB Corr
PBmin     = 0.05   # Min PB gain for correction
doBlank   = False  # Blank total intensity and SI where highest plane blanked?
                   # After any PB Correction
nthreads  = 16     # Number of threads to use
post      = None   # If given, the ending of file names,
                   # Else the value from AvgWeights.pickle
fitsfiles = None   # If given, fine name root, RMS
                   # Else the value from AvgWeights.pickle
opost     = "_5Pln.fits" # Post for output
doBlank   = True  # Blank total intensity and SI where highest plane blanked?
import Image, OErr, OSystem, FArray
OSystem.PAllowThreads(nthreads)  # with threads
from PipeUtil import SaveObject, FetchObject
import math, pickle

# Get weighting info
stuff = FetchObject("AvgWeights.pickle")
refFreq   = stuff['EffFreq'] # Effective frequency in MHz
wwts      = stuff['weights'] # Subband weights in %
dir       = stuff['dir']     # Data directory
if fitsfiles==None:
    fitsfiles = stuff['RMSes']   # List of (file_name root, RMS)
if post==None:
    post = stuff['post']    # Ending of spectral cube file name

import Obit, ImageMF, Image, FArray, OSystem, OErr
err = OErr.OErr()
OSystem.PAllowThreads(nthreads)  # with threads

for fitsfile in fitsfiles:
    print ("Doing",fitsfile[0])
    inMF   = ImageMF.newPFImageMF('Input',dir+fitsfile[0]+post, 0, True, err)
    outMF  = ImageMF.newPFImageMF('Output',dir+fitsfile[0]+opost, 0, False, err)
    ImageMF.PFitSpec2(inMF, outMF, err, \
                      corAlpha=defaultSI, doError=True, maxChi2=maxChi2,calFract=calFract, \
                      refFreq=refFreq*1.0e6, Weights=wwts, minFlux=minFlux, doPBCor=doPBCor)
    if doBlank:
        print ('Blanking by high freq. bin')
        # Blanking mask
        outIm = Image.newPFImage('Output',dir+fitsfile[0]+post, 0, True, err)
        hiPlane = [outIm.Desc.Dict['inaxes'][2],1,1,1,1]
        if doPBCor:
            # Mask by primary beam correction and high bin
            PBImg = Image.Image('PBImage'); Image.PCloneMem(outIm,PBImg,err)
            ImageUtil.PPBImage(outIm,PBImg,err,PBmin,outPlane=hiPlane,antSize=antSize)
            PBImg.GetPlane(None,hiPlane, err); mask1 = PBImg.FArray
            outIm.GetPlane(None,hiPlane, err); FArray.PBlank(mask1, outIm.FArray, mask1)
        else:
            # Mask by highest frequency
            outIm.GetPlane(None,hiPlane, err); mask1 = FArray.PCopy(outIm.FArray, err)
            PBImg = None
        # Plane 1 - total intensity
        outIm  = Image.newPFImage('Output',dir+fitsfile[0]+opost, 0, True, err)
        outIm.GetPlane(None, [1,1,1,1,1], err)
        FArray.PBlank(outIm.FArray, mask1, outIm.FArray)
        outIm.PutPlane(None, [1,1,1,1,1], err)
        # Plane 2 - Spectral index
        outIm.GetPlane(None, [2,1,1,1,1], err)
        FArray.PBlank(outIm.FArray, mask1, outIm.FArray)
        outIm.PutPlane(None, [2,1,1,1,1], err)
        Image.PFArray2FITS(mask1, "Mask.fits", err, outDisk=0)
        del PBImg
# end loop
