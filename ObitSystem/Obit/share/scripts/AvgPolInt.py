# python/Obit script to average polarized intensity in a Q/U pair of MFImage FITS files
# Arguments:
# 1) Name of Q MFImage, first plane to be rewritten
# 2) Name of U image, first plane to be rewritten
# 3) nThreads (optional)
# All files are in the current working directory

import sys, Obit, Image, FArray, CArray, FFT, OSystem, OErr, InfoList, History

# Init Obit
err=OErr.OErr()
ObitSys=OSystem.OSystem ("AvgMFPol", 1, 100, 1, ["None"], 1, ["./"], True, False, err)
OErr.printErrMsg(err, "Error with Obit startup")

# Get file names
inQFile = sys.argv[1]
inUFile = sys.argv[2]
if len(sys.argv)>3:
    nThreads = int(sys.argv[3])
    OSystem.PAllowThreads(nThreads)  # Threaded

print (ObitSys.pgmName+":", inQFile,inUFile)
print ("Avg Pol. intensity")

# Convert files into Images
inQ = Image.newPFImage("Q MF Cube", inQFile, 0, True, err)
inU = Image.newPFImage("U MF Cube", inUFile, 0, True, err)
OErr.printErrMsg(err, "Error initializing images")

# Check that MFImage files
if (inQ.Desc.Dict['ctype'][2]!='SPECLNMF') or (inU.Desc.Dict['ctype'][2]!='SPECLNMF'):
    OErr.PLog(err,OErr.Fatal,"Input NOT MFImage files"); OErr.PSet(err)
    exit

# Check the same dimensions
qdim = inQ.Desc.Dict['inaxes'][0:3]
udim = inU.Desc.Dict['inaxes'][0:3]
if (qdim[0]!=udim[0]) or (qdim[1]!=udim[1]) or (qdim[2]!=udim[2]):
    OErr.PLog(err,OErr.Fatal,"Inputs have different dimensions"); OErr.PSet(err)
    exit

# Info
nspec = inQ.Desc.List.Dict['NSPEC'][2][0]
nterm = inQ.Desc.List.Dict['NTERM'][2][0]

# Working arrays
cplane = CArray.CArray('plane', naxis=qdim[0:2])   # Plane complex array
fppol  = FArray.FArray('ppol', naxis=qdim[0:2])    # Plane ppol
faccum = FArray.FArray('accum', naxis=qdim[0:2])   # accumulated weighted ppol
FArray.PFill(faccum, 0.0)                          # Zero accumulator
sumWt = 0.0
# Loop over planes:
for ipln in range (0,nspec):
    pln = [ipln+nterm+1,1,1,1,1]  # Plane in cube 1 rel
    inQ.GetPlane(None, pln, err)
    inU.GetPlane(None, pln, err)
    OErr.printErrMsg(err, "Error reading plane "+str(ipln+1))
    CArray.PComplex(inQ.FArray, inU.FArray, cplane)  # Complex plane
    CArray.PAmp(cplane, fppol)                       # Pol intensity
    rms = (inQ.FArray.RMS+inU.FArray.RMS)*0.5
    if rms>0.0:
        wt = 1.0/rms;                                # Weighting
        sumWt += wt
        FArray.PSMul(fppol, wt)                      # weight
        FArray.PAdd(faccum, fppol, faccum)           # Accumulate
        #print (pln[0],"wt",wt)

# end loop over planes

# Normalize
norm = 1.0/sumWt
FArray.PSMul(faccum, norm)

# Write to first plane
pln = [1,1,1,1,1]
inQ.PutPlane(faccum, pln, err)
inU.PutPlane(faccum, pln, err)
OErr.printErrMsg(err, "Error writing results")

# Do history 
outHistory = History.History("history", inQ.List, err)
outHistory.Open(History.READWRITE, err)
# Only timestamp needed
outHistory.TimeStamp(" Start Obit "+ObitSys.pgmName,err)
outHistory.WriteRec(-1,ObitSys.pgmName+"  / Avg. Pol. intensity",err)
outHistory.WriteRec(-1,ObitSys.pgmName+" inQFile = "+inQFile,err)
outHistory.WriteRec(-1,ObitSys.pgmName+" inUFile = "+inUFile,err)
outHistory.Close(err)
outHistory = History.History("history", inU.List, err)
outHistory.Open(History.READWRITE, err)
outHistory.TimeStamp(" Start Obit "+ObitSys.pgmName,err)
outHistory.WriteRec(-1,ObitSys.pgmName+"  / Avg. Pol. intensity",err)
outHistory.WriteRec(-1,ObitSys.pgmName+" inQFile = "+inQFile,err)
outHistory.WriteRec(-1,ObitSys.pgmName+" inUFile = "+inUFile,err)
outHistory.Close(err)
OErr.printErrMsg(err, "Error with history")

# Shutdown Obit
OErr.printErr(err)
del ObitSys

