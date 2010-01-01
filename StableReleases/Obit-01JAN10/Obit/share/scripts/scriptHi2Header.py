# script to copy history from a FITS table to the FITS header
# FITS images only, works in current directory
# Argument:
# 1) Name of input FITS
# example:
# Python scriptHi2Header.py myImage.fits

import sys, Obit, Image, History, OSystem, OErr

# Init Obit
err=OErr.OErr()
ObitSys=OSystem.OSystem ("Hi2Header", 1, 100, 1, ["None"], 1, ["./"], 1, 0, err)
OErr.printErrMsg(err, "Error with Obit startup")

# Files (FITS)
inFile   = sys.argv[1]
inDisk   = 0

# Set data
inImage   = Image.newPFImage("Input image", inFile,  inDisk,   1, err)
OErr.printErrMsg(err, "Error initializing")

# For debugging
#Obit.Bomb()

# Make history
inInfo  = Image.PGetList(inImage)
outInfo = Image.PGetList(inImage)
inHistory  = History.History("history", inInfo, err)
outHistory = History.History("history", outInfo, err)
OErr.printErrMsg(err, "Error initializing history")

History.PCopy2Header(inHistory, outHistory, err)
OErr.printErrMsg(err, "Error copying history to FITS header")

# Say something
print "Copied History table to FITS header for",inFile

# Shutdown Obit
OErr.printErr(err)


