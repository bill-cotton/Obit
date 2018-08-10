# python/Obit script to fetch an NVSS postage stamp
# The output image is the input image written as 16 or 32 bit integers
# with a quantization level a given fraction of the minimum RMS in any plane.
# Arguments:
# 1) Name of output"FITS image to be quantized
# 2) RA as string "hh mm ss.ss"
# 3) Dec as string "sdd mm ss.ss"
# All files in ./
# example frorm command line:
# ObitTalk scriptFetchNVSS.py NVSSOut.fits "1 2 3.4" "+45 0 0.0"
#
import sys, Image, OSystem, OSurvey, OErr

# Init Obit
err=OErr.OErr()
ObitSys=OSystem.OSystem ("FetchNVSS", 1, 100, 1, ["./"], 1, ["./"], 1, 0, err)
OErr.printErrMsg(err, "Error with Obit startup")

# Get Arguments
outFile= sys.argv[1]
ra     = sys.argv[2]
dec    = sys.argv[3]
# Parameters
Size   = '0.50 0.50'   # Field of view in deg
Cells  = [15.0, 15.0]  # Cell spacings in asec
OSurvey.PNVSSFetch(ra, dec, outFile, err, Size=Size, Cells=Cells)
# Shutdown Obit
OErr.printErr(err)
del ObitSys

