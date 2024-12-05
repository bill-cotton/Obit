# Example script to make mosaic of a set of overlapping FITS images
#exec(open("Mosaic.py").read()) # to execute from ObitTalk

# NOTE!!!!!! images MUST have frequency axis in frequency rather than velocity
import os
exec(open("Targ.py").read())          # Target mosaic list
exec(open("PointingList.py").read())  # Pointing list, read to pointings
pointings = []
for p in data:
    pointings.append(p[0])
#################################################################################################
IStream  = 1         # Processing stream
strId    = 'S'+str(IStream)+stktrans   # Stream Id, unique, <= 5 char
restart  = False     # Restarting?
ch0      = 0         # restart first at channel ch0
minWt    = 0.35      # Minimum weight
doGPU    = False     # Use GPU?
minCh    = 1         # Minimum channel number
maxCh    = 16        # Maximum channel number
antSize  = 25.       # Antenna diameter (m)
datadisk  = 0        # disk for data (FITS) 0=CWD
fitsdir   = "./"     # FITS data directory
accumdisk = 1        # AIPS disk for accumulators  1=>RAM disk if setup
accumseq  = 10       # AIPS seq for accumulators
maxd      = 0.6      # How far from pointing center (deg) for overlap

nThreads = 24        # How many threads
prtLv    = 1         # Diagnostic print level, 0, 1, 2

OSystem.PAllowThreads(nThreads) # enable threading
print ("Target=",targets[0][0],"restart=",restart,"ch0=",ch0,"minWt=",minWt)
OErr.PInit(err, taskLog ='Stream'+str(IStream)+'.log')
OErr.PLog(err, OErr.Info, " Start processing "+strId);
exec(open("mosaicBasic.py").read())  # Do it
OErr.PLog(err, OErr.Info, "End processing "+strId);
OErr.printErr(err)
OSystem.Shutdown()

