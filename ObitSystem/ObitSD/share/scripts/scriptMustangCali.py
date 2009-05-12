# Template script for bulk calibration of Mustang data
# This script applies file global calibrations which do not depend on
# a source model.
# Optionally can update OTF data from the GBT archive
# Calibration tables left on the OFT data file.
import OSystem, OErr, InfoList, History
import OTF, GBTUtil, FITSDir, PARCal, Obit

# Init Obit
err=OErr.OErr()
FITS = ["./FITSdata"]  # Where to put output/scratch files
ObitSys=OSystem.OSystem  ("MustangCali", 1, 1, 1, ["None"], \
                          len(FITS), FITS, True, False, err)
OErr.printErrMsg(err, "Error with Obit startup")

#######################################################################################
# Define parameters

# Root of data directory
DataRoot="/home/gbtdata/AGBT09A_052_01/"
DataRoot = None  # To suppress attempt to update from archive

# Define data
# OTF file
inFile   = "AGBT09A_052_01OTF.fits" # OTF data file
#           ^^^^^^^^^^^^^^^^^^^^^   SET THIS
inDisk   = 0                        # 0 means current working directory
config = "./PAROTF09.cfg"           # Detector configuration file if reading data

# Default Calibration info
flagver    = 1                      # Flag table to apply
CalJy      = [38.5]                 # Cal values in Jy, one for all or per detector
BLInt      = 30.0                   # Baseline filter time in sec
AtmInt     = 20.0                   # Atmospheric filter time in sec

# Table of pointing offsets in time order [time(day) d Xel (asec), d el (asec)]
PointTab=[\
    [0.0, 0.0, 0.0], \
    [1.0, 0.0, 0.0]]
#    *********** SET THIS ********
# Example
#PointTab=[
#    [ 0.35232,    1.30,  -10.34], #   1337-1257, Beam=( 9.25,  8.26,   -8.3), Peak= 8.04 \ 
#    [ 0.43735,    2.09,   -8.41], #   1337-1257, Beam=( 8.95,  8.71,   58.7), Peak= 7.57 \ 
#    [ 0.45643,    2.25,   -8.80]  #   1337-1257, Beam=( 9.01,  8.89,   31.0), Peak= 7.16 \ 
#    ]

# Table of opacities in time order [time(day), zenith opacity(nepers)]
tau0 = [[0.0000,    0.100],       \
        [1000.0000, 0.100]]
#    *********** SET THIS if needed ********

################################# Update data #########################################
# Update OTF with any new scans in archive, 50 msec averaging
inOTF = GBTUtil.UpdateOTF ("PAROTF","Rcvr_PAR",inFile, inDisk, DataRoot, err, \
                           avgTime=0.05, config=config)
OErr.printErrMsg(err, "Error creating/updating input data object")
inInfo = inOTF.List

################################## Initial calibration #############################
PARCal.InitCal(inOTF, ["    "], err,\
               flagver=flagver, CalJy=CalJy, BLInt=BLInt, AtmInt=AtmInt,tau0=tau0, \
               PointTab=PointTab)

##################################### History #########################################
print "Add history"
# Loop over dirty, clean images, output data
outHistory = History.History("out history", inOTF.List, err)

# Add this programs history
outHistory.Open(History.READWRITE, err)
outHistory.TimeStamp(" Start Obit "+ObitSys.pgmName,err)

# Calibration
outHistory.WriteRec(-1,ObitSys.pgmName+" BLInt = "+str(BLInt),err)
outHistory.WriteRec(-1,ObitSys.pgmName+" AtmInt = "+str(AtmInt),err)
outHistory.Close(err)

# Copy history to header
inHistory  = History.History("in history",  inOTF.List, err)
outHistory = History.History("out history", inOTF.List, err)
History.PCopy2Header(inHistory, outHistory, err)
OErr.printErrMsg(err, "Error with history")

# Shutdown Obit
OErr.printErr(err)
del ObitSys
