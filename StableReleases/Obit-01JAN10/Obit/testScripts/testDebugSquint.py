# test for task Squint, must be run by ObitTalk rather than python
# Tests application of calibration, imaging, self calibration, autoWindow,
# autoCenter, and Beam Squint correction with systematics limited data
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
ObitSys=OSystem.OSystem ("testSquint", 1, 100, 0, ["None"], 1, ["../testIt/"], \
                         True, False, err)
OErr.printErrMsg(err, "Error with Obit startup")

# Get file names
inDisk   = 1
inFile = "SquintTestIn.uvtab"
outDisk  = 1
outFile = "SquintTestOut.fits"
out2Disk  = 1
out2File = "SquintTestOut.uvtab"
masterDisk = 1
masterFile  = 'SquintTestMaster.fits'
masterUVDisk = 1
masterUVFile  = 'SquintTestMaster.uvtab'
stokes='I'
source = '0319+415H2'

sq=ObitTask.ObitTask("Squint")

# DEBUGGING
sq.debug=True

# Set data files
sq.DataType='FITS'
sq.inFile=inFile; sq.inDisk=inDisk
sq.outFile=outFile; sq.outDisk=outDisk
sq.out2File=out2File; sq.out2Disk=out2Disk

# Imaging parameters
sq.doCalib=2; sq.doBand=1;
sq.Sources[0]=source
sq.FOV = 0.35
sq.Stokes = stokes
sq.Niter=1000
sq.minFlux=0.01
sq.autoWindow=True
sq.autoCen=1.0
sq.PBCor=True; sq.antSize=25.
sq.nThreads = 2
sq.logFile="SquintDebug.log"

# NVSS catalog
sq.OutlierDist = 1.0   # Maximum distance to add outlyers (deg)
sq.OutlierFlux = 0.005 # Minimum estimated outlier flux density (Jy)
sq.OutlierSI   = -0.7  # Spectral index to estimate flux density
sq.OutlierSize = 50    # Size of outlyer field

# Self calibration
sq.maxPSCLoop = 2
sq.minFluxPSC = 0.1
sq.solPInt    = 0.5
sq.solPType   = 'L1'
sq.solPMode   = 'P'
sq.maxASCLoop = 1
sq.minFluxASC = 0.5
sq.solAInt    = 2.0
sq.refAnt     = 21
sq.dispURL    = 'ObitView'

# Run task
l = sq.g

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
print "Squint test, input ",inFile,"Clean",outFile

# Any final word
OErr.printErr(err)

# Shutdown Obit
OErr.printErr(err)
OSystem.Shutdown(ObitSys)
