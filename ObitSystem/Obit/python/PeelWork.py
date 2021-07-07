#----------------------mySource -------------------------------------
source = 'mySource'; gainuse=0;
AUname=source; AUclass='SCal'; AUdisk=2; AUseq=1   # Initial uv data
AIname=source; AIclass='IPScal'; AIdisk=1; AIseq=1 # Self cal image w/ CC table
# First source to peel
ra='01:02:18.71';dec='-75:46:52.6'; seq=1 

# Second source to peel
ra='01:11:33.96';dec='-75:38:10.5'; seq=2; gainuse=-1 
AUclass='UPeel'; AUdisk=2; AUseq=seq

# Third source to peel
ra='01:13:21.34';dec='-75:28:21.2'; seq=3; gainuse=-1; AUseq=seq

# parameters likely to change
oDisk = 1       # output disk
FOV   = 1.2     # Radius of full field of view (deg)
Robust=-1.5     # Briggs robust factor
Niter = 100000  # Maximum number of CLEAN components
minFlux=50.0e-6 # CLEAN depth (Jy)
nthreads=24     # How any threads?
doGPU = False   # Use GPU for model calculation

from PeelScripts import *

# For each pool source do steps 1-4
# May need to clean up uv data files after they are no longer needed

# 1) select other non peel CCs
x=Image.newPAImage('inIm',AIname, AIclass, AIdisk, AIseq, True,err)
rad = ImageDesc.PHMS2RA(ra); decd= ImageDesc.PDMS2Dec(dec)
peelPos = [rad,decd]
# Select within 60" of position.
SelectCC(x,seq,seq+1,60./3600,peelPos,err)

# 2) Subtract others
uvi=UV.newPAUV('inUV',AUname, AUclass, AUdisk, max(1,AUseq-1), True,err)
uvs=UVSub4Peel(uvi,source,x,seq+1,err,nThreads=nthreads,gainUse=gainuse)
uvs.doGPU=True; uvs.outSeq=seq
uvs.g

# 3) Image peel source, up to 4 loops 30 sec SI
us=UV.newPAUV('inUV',AUname, '4Peel', AUdisk, AUseq, True,err)
mm=ImagePeel(uvs,peelPos,err,nThreads=nthreads,maxPSCLoop=4,solPInt=0.5, \
             minFlux=0.0001,seq=seq)
mm.dispURL= "http://localhost:8765/RPC2"
mm.Sources=[source]; mm.Niter=5000; mm.minFluxPSC=0.005; mm.minFluxASC=0.005
mm.Robust=-1.5; mm.autoWin=False; mm.CLEANBox=[-1,5,0,0]
maxAWLoop = 1
addParam(mm,"maxAWLoop", paramVal=maxAWLoop, shortHelp="Max. middle CLEAN loop", \
        longHelp="  maxAWLoop....Max. middle CLEAN Loop count\n"+ \
        "              Override the default behavior for the autoWin middle loop\n"+
        "              if > 0.\n")
mm.doGPU=doGPU
mm.g

# 4) remove peel source
if (seq==1):
    uv = UV.newPAUV('inData', AUname, AUclass, AUdisk, seq, True,err)
else:
    uv = UV.newPAUV('inData', AUname, AUclass, AUdisk, seq-1, True,err)

uvp=UV.newPAUV('inUV', source[0:12], mm.out2Class, mm.out2Disk, mm.out2Seq, True,err)
imp=Image.newPAImage('inIm', source, mm.outClass, mm.out2Disk, mm.outSeq, True,err)
uvpeel=SubPeel(uv,source, imp, uvp, err, doGPU=doGPU, seq=seq, addBack=False,  nThreads=nthreads)

# Image data with no further selfcalibration, may want to fiddle parameters
mf=ObitTask('MFImage')
setname(uvpeel,mf);mf.Sources = [source]; 
mf.doCalib=-1;mf.gainUse=0; mf.flagVer=1; mf.OutlierDist=1.2*FOV; mf.OutlierFlux=0.001
mf.Stokes='I';mf.doBand=-1; mf.BPVer=1; mf.doPol=False; mf.PDVer=2; 
mf.outClass='IPeel';mf.outSeq=1; mf.FOV=FOV; mf.OutlierSize=515; mf.minPatch=500
mf.Niter=Niter; mf.BLFact=1.01; mf.PBCor=False; mf.prtLv=2; mf.autoWindow=True
mf.outDisk=oDisk; mf.out2Disk=oDisk; mf.out2Seq=1; mf.Catalog='AllSkyVZ.FIT'
mf.doGPU=doGPU; mf.nThreads=nthreads; mf.Robust=Robust; mf.maxPixel=1000000
mf.Gain=0.05;mf.autoWindow=True; mf.ccfLim=0.50; mf.minFlux=minFlux;
mf.maxPSCLoop=0; mf.minFluxPSC=0.01; mf.solPMode='P'; mf.solPType='L1';  mf.solPInt=0.50; 
mf.maxASCLoop=0; mf.minFluxASC=1.0; mf.solAMode='A&P'; mf.solAInt=10.0; mf.solAType='L1'
mf.avgPol=True; mf.autoCen=1000.0; mf.noNeg=False;
mf.dispURL= "http://localhost:8765/RPC2" # default ObitView
mf.logFile=source+'.MFImage.log'
maxAWLoop = 1
addParam(mf,"maxAWLoop", paramVal=maxAWLoop, shortHelp="Max. middle CLEAN loop", \
        longHelp="  maxAWLoop....Max. middle CLEAN Loop count\n"+ \
        "              Override the default behavior for the autoWin middle loop\n"+
        "              if > 0.\n")
minFList = [0.00007,0.00005,0.00005, 0.00005] 
addParam(mf,"minFList", paramVal=minFList, \
        shortHelp="minFlux list after SC", \
        longHelp="  minFList....minFluxes to use in IPol after selfcals\n")

mf.g  # run imaging

# Restore peeled sources
xpeel = Image.newPAImage('peel',source[0:12],'IPeel', mf.outDisk, mf.outSeq,True,err)
nseq = seq; CCVer=1
for iSeq in range(1,nseq+1):
    xi = Image.newPAImage('mod',source[0:12], mm.outClass, mm.outDisk, iSeq,True,err)
    RestorePeel(xi, CCVer, xpeel, err)

