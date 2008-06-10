# Program to do basic atmospheric and instrumental calibration
# of GBT Xband daisy test and write calibrated output.
import OTF, OSystem, OErr

# Init Obit
err=OErr.OErr()
ObitSys=OSystem.OSystem ("Python", 1, 103, 1, ["None"], 1, ["../FITSdata/"], 1, 0, err)
OErr.printErrMsg(err, "Error with Obit startup")

# Files
disk = 1
inFile   = "GBTDaisyXOTF.fits"         # input OTF data
outFile = "!GBTDaisyX2OTF.fits"        # output OTF data

# Set data
inData =OTF.newPOTF("Input data", inFile, disk, 1, err)
OErr.printErrMsg(err, "Error creating input data object")
outData = OTF.newPOTF("Output data", outFile, disk, 0, err)
OTF.POTFClone(inData, outData, err)     # Same structure etc
OErr.printErrMsg(err, "Error initializing")

# delete any prior calibration tables
OTF.ClearCal(inData,err)

# Atmospheric cal
OTF.AtmCalInput["InData"] = inData
OTF.AtmCalInput["solint"] = 15.0           # Solution interval in seconds
OTF.AtmCalInput["tau0"]   = 0.009          # Zenith opacity in nepers - guess
OTF.AtmCalInput["minEl"]  = 5.0            # min elevation deg
OTF.AtmCalInput["aTemp"]  = [0.735,0.722]  # from GC ob
OTF.AtmCalInput["tRx"]    = [9.252,9.230]  # to get zero about correct
OTF.AtmCalInput["tRx"]    = [8.900,9.500]  # to get zero about correct
OTF.AtmCalInput["calJY"]  = [1.644,1.676]  # From 3C286 GC obs
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
