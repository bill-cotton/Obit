# test for task IonImage, must be run by ObitTalk rather than python
# Tests application of calibration, imaging, self calibration, autoWindow,
# autoCenter, and Beam IonImage correction with systematics limited data
# Needs (smeagle)
#setenv FITS   /export/data_1/users/bcotton/Software.dir/SVN/ObitInstall/ObitSystem/Obit/testIt
#setenv FITS01 /export/data_1/users/bcotton/Software.dir/SVN/ObitInstall/ObitSystem/Obit/testIt
# Needs (gollum)
#setenv FITS   /export/data_1/users/bcotton/Software.dir/SVN/ObitInstall/ObitSystem/Obit/testIt
#setenv FITS01 /export/data_1/users/bcotton/Software.dir/SVN/ObitInstall/ObitSystem/Obit/testIt
# Needs (vino)
#setenv FITS   /export/data_1/users/bcotton/SVN/ObitInstall/ObitSystem/Obit/testIt
#setenv FITS01 /export/data_1/users/bcotton/SVN/ObitInstall/ObitSystem/Obit/testIt

import Obit, OErr, OSystem, UV, Image, ObitTask, FITS

# Init Obit
err=OErr.OErr()
ObitSys=OSystem.OSystem ("testIonImage", 1, 100, 0, ["None"], 1, ["../testIt/"], \
                         True, False, err)
OErr.printErrMsg(err, "Error with Obit startup")

# Get file names
inDisk   = 1
inFile = "IonImageTestIn.uvtab"
outDisk  = 1
outFile = "IonImageTestOutNoPeel.fits"
out2Disk  = 1
out2File = "IonImageTestOut.uvtab"
masterDisk = 1
masterFile  = 'IonImageTestMaster.fits'
masterUVDisk = 1
masterUVFile  = 'IonImageTestMaster.uvtab'
stokes='I'
source = '0900+398'

ii=ObitTask.ObitTask("IonImage")

# DEBUGGING
ii.debug=True

# Set data files
ii.DataType='FITS'
ii.inFile=inFile; ii.inDisk=inDisk
ii.outFile=outFile; ii.outDisk=outDisk
ii.out2File=out2File; ii.out2Disk=out2Disk

# Imaging parameters
ii.doCalib=2; ii.doBand=1;
ii.doCalib=-1; ii.doBand=-1;
ii.Sources[0]=source
ii.FOV = 5.0
ii.Stokes = stokes
ii.Niter=1000
ii.minFlux=0.1
ii.autoWindow=True
ii.autoCen=10.0
ii.PBCor=False; ii.antSize=25.
ii.nThreads = 2
ii.logFile="IonImageNoPeel.log"
ii.doFCal = True
ii.UpdateInt = 0.25
ii.solInt = 1.0
ii.nZern = 5
ii.MinPeak = 4.0
ii.MaxWt = 10.0
ii.doINEdit = True
ii.MaxQual = 1
ii.MaxRMS = 20.0
ii.MinRat = 0.5
ii.FCNiter = 200
ii.FCminFlux = 1.0

# NVSS catalog
ii.OutlierDist = 40.0  # Maximum distance to add outlyers (deg)
ii.OutlierFlux = 3.0   # Minimum estimated outlier flux density (Jy)
ii.OutlierSI   = -0.7  # Spectral index to estimate flux density
ii.OutlierSize = 100    # Size of outlyer field

# Peeling
ii.PeelFlux = 1000.0
ii.PeelLoop = 1
ii.PeelSolInt =  10.0
#ii.dispURL="None"
ii.PeelMinFlux = 0.50
ii.PeelRefAnt     = 21


# Run task
ii.i
l = ii.g

# output image - name has source_name+stokes prepended
outImage  = Image.newPFImage("output image",  source+stokes[0]+outFile, outDisk, True, err)

# Compare with master lie [rms diff, max abs diff, max. master]
masterImage  = Image.newPFImage("Master image",   masterFile,  masterDisk,  True, err)
diff = Image.PCompare(outImage, masterImage, err);
print "Comparison, rel. max. residual",diff[1]/diff[0], " rel RMS residual",diff[2]/diff[0]

# Compare with master UV data lie
outData     = UV.newPFUV("Output UVData",  source+out2File, out2Disk, True, err)
masterData  = UV.newPFUV("Master UVData", masterUVFile, masterUVDisk, True, err)
diff = UV.PUtilVisCompare(outData, masterData, err);
print "UV Comparison with master lie, fractional RMS R,I difference",diff

# Say something
print "IonImage test, input ",inFile,"Clean",outFile

# Any final word
OErr.printErr(err)

# Shutdown Obit
OErr.printErr(err)
OSystem.Shutdown(ObitSys)
