# Program to test editing of OTF data from python in Obit
# Test for Richard Prestages Q band data
import OTF, OTFUtil, OErr, OSystem
from Obit import Bomb

# Init Obit
err=OErr.OErr()
ObitSys=OSystem.OSystem  ("Python", 1, 103, 1, ["None"], 1, ["../PythonData/"], 1, 0, err)
OErr.printErrMsg(err, "Error with Obit startup")

# Bomb if needed for debugging
#Bomb()

# Files
disk = 1
inFile   = "QbandOTF.fits"     # input OTF data

# Editing
timerange  = [0.0,1.0e20]   # all times
chans      = [1,0]          # all channels
flagVer    = 1              # Flagging table 1
target     =  "Any"         # all targets
feed       = 2              # Flag microphonic data
stokes     = "100"          # Stokes flag
reason     = "Microphonic"  # reason string

# Set data
inData = OTF.newPOTF("Input data", inFile, disk, 1, err)
OErr.printErrMsg(err, "Error creating input data object")

# Make edit
OTFUtil.PFlag(inData,err,timerange,target,flagVer,chans,stokes,feed,reason)
feed = 4
OTFUtil.PFlag(inData,err,timerange,target,flagVer,chans,stokes,feed,reason)
OErr.printErrMsg(err, "Error writing flag table")

# tell results
print "Flag",inFile,"timerange",timerange,"stokes",stokes,'feed',feed

# Shutdown Obit
OErr.printErr(err)
