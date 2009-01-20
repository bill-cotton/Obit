# Program to do basic atmospheric and instrumental calibration
# of GBT Qband test and write calibrated output.
import OTF, OSystem, OErr, InfoList

# Init Obit
err=OErr.OErr()
ObitSys=OSystem.OSystem ("Python", 1, 103, 1, ["None"], 1, ["../PythonData/"], 1, 0, err)
OErr.printErrMsg(err, "Error with Obit startup")

# Files
disk = 1
inFile   = "QbandOTF.fits"         # input OTF data
outFile = "!QbandCalOTF.fits"      # output OTF data

# Set data
inData =OTF.newPOTF("Input data", inFile, disk, 1, err)
OErr.printErrMsg(err, "Error creating input data object")
outData = OTF.newPOTF("Output data", outFile, disk, 0, err)
OTF.PClone(inData, outData, err)     # Same structure etc
OErr.printErrMsg(err, "Error initializing")

# delete any prior calibration tables
OTF.ClearCal(inData,err)

# Select scans
dim = OTF.dim
dim[0] = 2; dim[1] = 1
scans = [108,133]
inInfo = OTF.PGetList(inData)
InfoList.PAlwaysPutInt(inInfo, "SCANS",  dim, scans)
dim[0] = 1
InfoList.PAlwaysPutBoolean (inInfo,"doCalSelect" , dim, [True])
print "inInfo",inInfo.me

# Atmospheric cal
OTF.AtmCalInput["InData"] = inData
OTF.AtmCalInput["solint"] = 120.0           # Solution interval in seconds
OTF.AtmCalInput["tau0"]   = 0.09           # Zenith opacity in nepers - guess
OTF.AtmCalInput["minEl"]  = 5.0            # min elevation deg
OTF.AtmCalInput["aTemp"]  = [0.735,0.722,0.722,0.722,0.722,0.722,0.722,0.722]  # from GC ob
OTF.AtmCalInput["tRx"]    = [18.0,13.8,26.9,21.5, 15.2,13.6,30.5,23.9]  # to get zero about correct
#OTF.AtmCalInput["tRx"]    = [5.4,   3.8,  8.1,6.3, 4.5,4.1, 9.2, 7.2 ]  # to get zero about correct
OTF.AtmCalInput["calJy"]  = [9.6,9.6,9.6,9.6,9.6,9.6,9.6,9.6]  # From 3C286 very crude
OTF.AtmCalInput["calJy"]  = [1.0,1.0,1.0,1.0,1.0,1.0,1.0,1.0]  # From 3C286 very crude
OTF.AtmCal(err, OTF.AtmCalInput)

# Split
OTF.SplitInput["InData"]  = inData
OTF.SplitInput["OutData"] = outData
OTF.SplitInput["average"] = 0   # No averaging in frequency
OTF.SplitInput["gainuse"] = 0   # apply highest gain table
OTF.SplitInput["flagver"] = -1  # no flagging table
OTF.Split(err, OTF.SplitInput)

# Shutdown 
OErr.printErr(err)
print 'Done: wrote calibrated data in',outFile
