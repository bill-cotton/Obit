# python/Obit script to calibrate Stokes Vpol SiO maser line cubes
# Vpol is calculated by averaging  the R/L gain ratio determined from
# Stokes Ipol and Vpol images is constant.A

import Image, History, ImageDesc, OSystem, FArray, OErr
from Obit import Bomb

# Init Obit
err=OErr.OErr()
myName = "VPolCal"
ObitSys=OSystem.OSystem (myName, 1, 100, 1, ["../AIPSdata/"], 1, ["../PythonData/"], 1, 0, err)
OErr.printErrMsg(err, "Error with Obit startup")

# Files (FITS)
IPolDisk = 1
IPolFile   = 'U_Ori_E_1.Icube.fits.gz'
VPolDisk = 1
VPolFile   = 'U_Ori_E_1.Vcube.fits'
outDisk = 1
VPolFile   = '00000+00000.PCUBE.gz'
IPolFile   = '00000+00000.PCUBE.gz'
outFile  = '!VPolCal.fits'

# Set clipping factors
RMSFactor = 10.0
MaxFactor = 0.20

# Set data
IImage = Image.newPImage("Input IPol image", IPolFile, IPolDisk, True, err)
VImage = Image.newPImage("Input VPol image", VPolFile, VPolDisk, True, err)
OErr.printErrMsg(err, "Error initializing")

# For debugging
#Bomb()

# Open images
IImage.Open(Image.READONLY, err)
VImage.Open(Image.READONLY, err)
OErr.printErrMsg(err, "Error opening input")

# Get image descriptors
IdescDict = IImage.Desc.Dict
Idim  = IdescDict["inaxes"]            # Dimensionality
VdescDict = VImage.Desc.Dict
Vdim  = VdescDict["inaxes"]            # Dimensionality

# Check that dimensions compatible
if (Idim[0]!=Vdim[0]) | (Idim[1]!=Vdim[1]) | (Idim[2]!=Vdim[2]):
    print "I dim",Idim,"V dim",Vdim
    raise RuntimeError,"I and V images have inconsistent dimensions"

# Get data FArrays
IData = IImage.FArray
VData = VImage.FArray

# Get average ratio - loop over planes
sumX = 0.0; sumWt = 0.0   # Zero accumulators
for plane in range(Idim[2]):
    print "Plane",plane+1
    planeDim = [plane+1,1,1,1,1]  # Only 3D
    IImage.GetPlane(None, planeDim, err)   # Read Ipol plane (IData)
    VImage.GetPlane(None, planeDim, err)   # Read Vpol plane (VData)
    # Get IPol Image plane statistics for clipping
    pos=[0,0]
    Imax = FArray.PMax(IData, pos)
    Irms = IData.RMS
    Imin = max(Imax*MaxFactor,Irms*RMSFactor)
    print "Max",Imax,"RMS",Irms,"clip",Imin
    # Blank noise and below DR limit
    WtData = FArray.PCopy(IData, err)        # Work array copy of IPol
    DivData = FArray.PClone(IData, err)      # Clone of IPol array
    OErr.printErrMsg(err, "Error copying FArray")
    FArray.PClipBlank(WtData, Imin, 1.0e30)
    # Determine mean fractional VPol
    FArray.PDiv(VData, WtData, DivData)
    # Print Mean ratio
    print "Mean ratio",DivData.Mean
    # Multiple fraction times weight
    FArray.PMul(DivData, WtData, DivData)
    # Accumulate
    sumX  = sumX  + DivData.Sum
    sumWt = sumWt + WtData.Sum

del WtData, DivData  # Some cleanup

# Get average instrumental VPol calibration
VCor = sumX/sumWt
print "Average fractional VPol",VCor
# debug VCor = 0.1


# Generate scratch file from VPol
tmpImage  = Image.PScratch(VImage, err)
Image.POpen(tmpImage, Image.WRITEONLY, err)   # Open
OErr.printErrMsg(err, "Error cloning output")
outData = tmpImage.FArray                     # Get data FArray

# Do history
inHistory  = History.History("history", VImage.List, err)
outHistory = History.History("history", tmpImage.List, err)
History.PCopyHeader(inHistory, outHistory, err)
History.POpen(outHistory, History.READWRITE, err)
History.PTimeStamp(outHistory," Start Obit "+ObitSys.pgmName,err)
History.PWriteRec(outHistory,-1,ObitSys.pgmName+" / IPol input = "+IPolFile,err)
History.PWriteRec(outHistory,-1,ObitSys.pgmName+" / VPol input = "+VPolFile,err)
History.PWriteRec(outHistory,-1,ObitSys.pgmName+" / Corrected VPol calibration by "+str(VCor),err)
History.PClose(outHistory, err)
OErr.printErrMsg(err, "Error with history")

VCor = min (0.3, VCor)
# Correct
for plane in range(Idim[2]):
    print "Correcting Plane",plane+1
    planeDim = [plane+1,1,1,1,1]  # Only 3D
    IImage.GetPlane(None, planeDim, err)   # Read Ipol plane (IData)
    VImage.GetPlane(None, planeDim, err)   # Read Vpol plane (VData)
    # Scale IPol by correction and subtract
    FArray.PSMul(IData, VCor)
    FArray.PSub(VData, IData, outData)
    tmpImage.PutPlane(None, planeDim, err) # Write corrected (outData)

# Close images
Image.PClose(IImage, err)
Image.PClose(VImage, err)
Image.PClose(tmpImage, err)
OErr.printErrMsg(err, "Error closing files")

# Copy to quantized integer image with history
outImage = Image.newPImage("Output image", outFile, outDisk, False, err)
Image.PCopyQuantizeFITS (tmpImage, outImage, err, inHistory=outHistory)

# debug
#outImage.Zap(err)  # don't work?
Image.PZap(outImage, err)

# Shutdown Obit
OErr.printErr(err)
del ObitSys

