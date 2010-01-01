#testbed for Mosaicing utilities
import Obit, OErr, OSystem, Image, MosaicUtil, FArray

# Init Obit
err=OErr.OErr()
ObitSys=OSystem.OSystem ("Weight", 1, 100, 1, ["../AIPSdata/"], 1, ["../PythonData/"], 1, 0, err)
OErr.printErrMsg(err, "Error with Obit startup")

# For debugging
#print sys.argv
#Obit.Bomb()

# Get file names
tmplFile = "C0000P00.gz"      # Template
tmpFile1 = "testSum1.fits"    # work accumulation
tmpFile2 = "testSum2.fits"    # work accumulation
outFile  = "Comb0000P00.fits" # output
inDisk   = 1
outDisk  = 1

# Convert files into Images
tmplImage  = Image.newPFImage("Template image", tmplFile, inDisk,  True, err)
SumWtImage = Image.newPFImage("Accum image 1", tmpFile1, outDisk, False, err)
SumWt2     = Image.newPFImage("Accum image 2", tmpFile2, outDisk, False, err)
OErr.printErrMsg(err, "Error initializing images")

# Do it
MosaicUtil.PMakeMaster(tmplImage, [1000,1000], SumWtImage, SumWt2, err)
OErr.printErrMsg(err, "Error initializing images")
print "Made accumulation images "

print "Made master images ",tmpFile1,tmpFile2

del tmplImage # cleanup

# Weight them together
inFile=["00000+00000.PCUBE.gz","00000-00260.PCUBE.gz","00000+00260.PCUBE.gz"]
factor = 1.0
for x in inFile:
    print "Accumulate",x
    inImage = Image.newPFImage(x, x, inDisk,  1, err)
    MosaicUtil.PWeightImage(inImage, factor, SumWtImage, SumWt2, err)
    OErr.printErrMsg(err, "Error accumulating image "+x)

# Create output
outImage = Image.newPFImage("Output image", outFile, outDisk, False, err)
Image.PClone(SumWtImage, outImage, err)   # Same structure etc. as SumWtImage
OErr.printErrMsg(err, "Error creating output")

# Normalize
MosaicUtil.PNormalizeImage(SumWtImage, SumWt2, outImage, err, minWt=0.2)
OErr.printErrMsg(err, "Error normalizing output")

# Say something
print "Made accumulation images",tmpFile1, tmpFile2,"from",tmplFile,"to",outFile

# Shutdown Obit
OErr.printErr(err)
del ObitSys
