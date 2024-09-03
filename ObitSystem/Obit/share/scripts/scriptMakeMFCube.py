# Generate an MFImage-like image cube from FITS images
# Files should be in the current working directory
# Arguments:
# 1) Path of output MFImage style image to create
# 2) Path of a template MFImage style image to clone
# 3) Number of spectral planes for the output image

import sys, OSystem, OErr, InfoList, Image, ImageMF, ImageDesc, History

# Init Obit
err=OErr.OErr()
ObitSys=OSystem.OSystem  ("MakeMFImage", 1, 100, 1, ["None"], 1, ["."], 1, 0, err)
OErr.printErrMsg(err, "Error with Obit startup")

# Get file names, nspec
outFile = sys.argv[1]
inFile = sys.argv[2]
nspec = int(sys.argv[3])
inDisk = 1
outDisk = 1

# Convert files into ImageMFs
inImage = ImageMF.newPFImageMF(inFile, inFile, inDisk, True, err)
OErr.printErrMsg(err, "Error initializing input image")

# output image
outImage = ImageMF.newPFImageMF(outFile, "!"+outFile, outDisk, False, err)
OErr.printErrMsg(err, "Error initializing output image")
if OErr.PIsErr(err):
    OErr.printErr(err); raise RuntimeError("File Error")

# Make sure input if MFImage type
id = inImage.Desc
od = outImage.Desc
if id.Dict['ctype'][2] != "SPECLNMF":
    OErr.PLog(err, OError.Error, "Template image NOT MFImage type")
    OErr.printErr(err); raise RuntimeError("Improper input")

# Modify header adding number of planes
ImageDesc.PCopyDesc(id, od, err)
dict = od.Dict
inaxes = id.Dict["inaxes"]
inaxes[2] = id.List.get("NTERM")[4][0] + nspec
dict["naxis"]  = len(inaxes)
dict["inaxes"] = inaxes
dict["bitpix"] = -32
od.Dict = dict
# Update output Dict.List
od.List.set('NSPEC',nspec)

# Reset Frequency Info in Desc.List
for i in range(1,nspec+1):
    key = 'FREQ%4.4d'%i; od.List.set(key,-1.0)
    key = 'FREH%4.4d'%i; od.List.set(key,-1.0)
    key = 'FREL%4.4d'%i; od.List.set(key,-1.0)

# Create
outImage.Open(ImageMF.WRITEONLY,err)
# Close
outImage.Close(err)
if OErr.PIsErr(err):
    OErr.printErr(err); raise RuntimeError("File Error")

# Add history as Image
outImage = Image.newPFImage(outFile, outFile, outDisk, True, err)
OErr.printErrMsg(err, "Error initializing output history")
outHistory  = History.History("history", outImage.List, err)
z=outHistory.Open(History.WRITEONLY, err)
outHistory.TimeStamp(" Start Obit "+ObitSys.pgmName,err)
z=outHistory.WriteRec(-1,ObitSys.pgmName+" / clone from "+inFile,err)
z=outHistory.WriteRec(-1,ObitSys.pgmName+" / nspec = "+str(nspec),err)
z=outHistory.Close(err)
outImage.UpdateDesc(err)

# Say something
print ("Cloned",inFile,"into",outFile," with ",nspec,"Spectral planes")

# Shutdown Obit
OErr.printErr(err)
OSystem.Shutdown(ObitSys)
