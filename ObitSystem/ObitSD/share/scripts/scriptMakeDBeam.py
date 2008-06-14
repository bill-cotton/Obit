# Program to generate a dirty beam
import Image, ImageDesc, OSystem, OErr, FArray

# Init Obit
err=OErr.OErr()
ObitSys=OSystem.OSystem ("MakeDBeam", 1, 103, 1, ["None"], 1, ["./FITSdata/"], 1, 0, err)
OErr.printErrMsg(err, "Error with Obit startup")

# Files
disk = 1
outFile  = "!PARDBeam1.5.fits"     # output image

# Create image object
outImage=  Image.newPFImage("Output image", outFile, disk, 0, err)
OErr.printErrMsg(err, "Error initializing")

# define image
RA      = 0.0
Dec     = 0.0
nx      = 51
ny      = 51
xCells  = 2.0 / 3600.0
yCells  = 2.0 / 3600.0
beamsize = 8.0/3600.0  # beamsize in deg.

# Set image parameters in descriptor
imagDesc = outImage.Desc
imagDict = imagDesc.Dict
imagDict["ctype"][0] = "RA---SIN"
imagDict["ctype"][1] = "DEC--SIN"
imagDict["crval"][0] = RA
imagDict["crval"][1] = Dec
imagDict["cdelt"][0] = -abs(xCells) # RA goes backwards
imagDict["cdelt"][1] = yCells
imagDict["naxis"] = 2
imagDict["inaxes"][0] = nx
imagDict["inaxes"][1] = ny
imagDict["crpix"][0] = nx/2+1
imagDict["crpix"][1] = ny/2+1
imagDict["bitpix"] = -32
imagDict["bunit"] = "JY/BEAM"
imagDict["beamMin"] = beamsize
imagDict["beamMaj"] = beamsize
imagDict["beamPA"] = 0.0
# Update descriptor
outImage.Desc.Dict = imagDict

# 
Image.PFullInstantiate (outImage, 2, err)
OErr.printErrMsg(err, "Error creating image")

# Create FArray for image
data = FArray.FArray("Image Data",[nx,ny])

# Add Gaussian
FWHM = beamsize/abs(xCells) # beam size in pixels
FArray.PCGauss2D(data, [nx/2,ny/2], FWHM)

# Write output
Image.PWritePlane (outImage, data, err)

# Say something
print "Wrote beam with size",beamsize*3600.0,"to",outFile

# Shutdown Obit
OErr.printErr(err)
del ObitSys

