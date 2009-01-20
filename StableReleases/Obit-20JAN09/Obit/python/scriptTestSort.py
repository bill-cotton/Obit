# Test Sorting

import Obit, Image, Table, OSystem, OErr

# Init Obit
err=OErr.OErr()
ObitSys=OSystem.OSystem ("TestSort", 1, 100, 1, ["../AIPSdata/"], 1, ["../PythonData/"], 1, 0, err)
OErr.printErrMsg(err, "Error with Obit startup")

# Files (FITS)
inDisk = 1
inFile   = 'GaussTest.fits'
#outFile  = 'GaussTest.fits'
# AIPS test AIPS image user 100 C346R422    .POLI  .   2,, disk 8/142
#aname = "C346R422"
#aclass = "POLI"
#aseq = 2

# Set data
inImage   = Image.newPFImage("Input image", inFile,  inDisk,   1, err)
#outImage   = Image.newPAImage("Output image", aname, aclass, inDisk, aseq, 1, err)
OErr.printErrMsg(err, "Error initializing")

# For debugging
#Obit.Bomb()

# Make Table
inTable = Image.PImageGetTable (inImage, 3, "AIPS VL", 1, err)
OErr.printErrMsg(err, "Eror making table")

#Sort
Table.PSort (inTable, "DEC(2000)", 0, err)
OErr.printErrMsg(err, "Error sorting")

# Say something
print "Sorted VL table in",inFile

# Shutdown Obit
OErr.printErr(err)
del ObitSys

