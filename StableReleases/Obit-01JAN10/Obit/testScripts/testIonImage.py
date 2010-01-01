# test for task IonImage - run from python

import Obit, OErr, OSystem, UV, Image, ObitTask, AIPS, FITS



# On Smeagle
adirs = ["/export/data_1/aips/DATA/SMEAGLE_1",
         "/export/data_1/aips/DATA/SMEAGLE_2", \
         "/export/data_1/aips/DATA/SMEAGLE_3", \
         "/export/data_2/aips/DATA/SMEAGLE_4", \
         "/export/data_2/aips/DATA/SMEAGLE_5", \
         "/export/data_2/aips/DATA/SMEAGLE_6"]
fdirs = ["../testIt","/export/data_1/bcotton/Software.dir/AIPS/FITS"]

# On Gollum
adirs = ["/export/data_1/GOLLUM_1",
         "/export/data_1/GOLLUM_2", \
         "/export/data_1/GOLLUM_3", \
         "/export/data_1/GOLLUM_4", \
         "/export/data_2/GOLLUM_5", \
         "/export/data_2/GOLLUM_6", \
         "/export/data_2/GOLLUM_7", \
         "/export/data_2/GOLLUM_8"]
fdirs = ["../testIt","/export/users/aips/FITS"]

## Init Obit
err=OErr.OErr()
user = 105
ObitSys=OSystem.OSystem ("testIonImage", 1, user, len(adirs), adirs, \
                         len(fdirs), fdirs, True, False, err)
OErr.printErrMsg(err, "Error with Obit startup")

# Setup AIPS, FITS
AIPS.AIPS.userno = user
disk = 0
for ad in adirs:
    disk += 1
    AIPS.AIPS.disks.append(AIPS.AIPSDisk(None, disk, ad))
disk = 0
for fd in fdirs:
    disk += 1
    FITS.FITS.disks.append(FITS.FITSDisk(None, disk, fd))

# Get file names
inDisk        = 1
inFile        = "IonImageTestIn.uvtab"
outDisk       = 1
outFile       = "IonImageTestOut.fits"
out2Disk      = 1
out2File      = "IonImageTestOut.uvtab"
masterDisk    = 1
masterFile    = 'IonImageTestMaster.fits'
masterUVDisk  = 1
masterUVFile  = 'IonImageTestMaster.uvtab'
stokes='IV'
source = '1100+398'

ionimg=ObitTask.ObitTask("IonImage")

# DEBUGGING
#ionimg.debug = True

# Set data files
ionimg.DataType = 'FITS'
ionimg.inFile   = inFile;   ionimg.inDisk   = inDisk
ionimg.outFile  = outFile;  ionimg.outDisk  = outDisk
ionimg.out2File = out2File; ionimg.out2Disk = out2Disk

# Imaging parameters
ionimg.Sources[0] = source
ionimg.FOV        = 5.0
ionimg.Niter      = 5000
ionimg.minFlux    = 0.05
ionimg.autoWindow = True
ionimg.nThreads   = 2
ionimg.xCells     = 20.0;
ionimg.yCells     = 20.0; 
ionimg.WtBox      = 1
ionimg.Beam       = [80.0,80.0,0.0]
ionimg.Robust     = -5.0
ionimg.prtLv      = 2
ionimg.autoCen    = 20.0

# NVSS catalog/calibration
ionimg.OutlierDist = 20.0     # Maximum distance to add outlyers (deg)
ionimg.OutlierFlux = 3.0      # Minimum estimated outlier flux density (Jy)
ionimg.OutlierSI   = -0.7     # Spectral index to estimate flux density
ionimg.OutlierSize = 100      # Size of outlyer field
ionimg.solInt      = 1.0      # Solution interval
ionimg.nZern       = 5        # Number of Zernike terms
ionimg.MinPeak     = 3.0      # Min. acceptable snapshot peak
ionimg.MaxRMS      = 20.0     # Max acceptable seeing
ionimg.dispURL     = "None"   # No image display
ionimg.PeelFlux    = 20.0     # Min. flux for Peeling
ionimg.PeelRefAnt  = 18       # Peeling reference antenna
ionimg.PeelLoop    = 1        # Number of peeling self cal loops
ionimg.PeelSNRMin  = 3.0      # Min peeling self cal SNR
ionimg.PeelSolInt  = 1.0      # Peeling self cal solution interval
ionimg.PeelType    = 'L1'     # Peeling self cal soln type
ionimg.PeelMode    = 'P'      # Peeling self cal soln mode
ionimg.PeelNiter   = 200      # Peeling self cal max CLEAN iterations
ionimg.PeelMinFlux = 0.05     # Peeling self cal mix CLEAN flux

ionimg.logFile="IonImage.log" # Also used for comparison

#ionimg.i
# Run task
ionimg.g

# output image - name source_name+"I"+outFile
outImage  = Image.newPFImage("output image", source+"I"+outFile, outDisk, True, err)

# Compare with master lie [rms diff, max abs diff, max. master]
masterImage  = Image.newPFImage("Master image",   masterFile,  masterDisk,  True, err)
diff = Image.PCompare(outImage, masterImage, err);
line = "Comparison, rel. max. residual "+str(diff[1]/diff[0])+" rel RMS residual "+str(diff[2]/diff[0])
print line
fd = open(ionimg.logFile, "a")
fd.write(line+"\n")
fd.close()

# Compare with master UV data lie
outData     = UV.newPFUV("Output UVData",  source+out2File, out2Disk, True, err)
masterData  = UV.newPFUV("Master UVData", masterUVFile, masterUVDisk, True, err)
diff = UV.PUtilVisCompare(outData, masterData, err);
line = "UV Comparison with master lie, fractional RMS R,I difference "+str(diff);
print line
fd = open(ionimg.logFile, "a")
fd.write(line+"\n")
fd.close()

# Say something
print "IonImage test, input ",inFile,"Clean image",outFile

# Any final word
OErr.printErr(err)

# Shutdown Obit
OErr.printErr(err)
OSystem.Shutdown(ObitSys)
