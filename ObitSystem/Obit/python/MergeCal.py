# Function to run the equivalent of AIPSish MERGECAL
# Fixes VLBA screwed up calibration tables.

from AIPSTask import AIPSTask

def MergeCal(inUV, err, \
             GCver=1, TYver=1, PCver=1, outVer=0, timeTol=1.0):
    """ Fix up VLBA screwed up calibration

    Merges redundant entries in VLBA GC, TY and PC tables
    Translated from the AIPSish MERGECAL (requires AIPS task TAMRG)
    inUV     = UV data object to fix, MUST be in AIPS format
    err      = Python Obit Error/message stack
    GCver    = Version number of GC table to be fixed.
    TYver    = Version number of TY table to be fixed.
    PCver    = Version number of PC table to be fixed.
    outVer   = Version number of output tables
               (same for all three table types). 
    timeTol  = Tolerance for time comparisons in seconds.
               Records will not be merged if their times differ
               by more than this amount.
    """
    ################################################################
    if inUV.FileType!='AIPS':
        raise RuntimeError,"Can ONLY handle AIPS data"
    # Set up
    tamrg = AIPSTask("tamrg")
    tamrg.inname   = inUV.Aname
    tamrg.inclass  = inUV.Aclass
    tamrg.indisk   = inUV.Disk
    tamrg.inseq    = inUV.Aseq

    # GC table
    tamrg.inext     = "GC"
    tamrg.invers    = GCver
    tamrg.outvers   = outVer
    tamrg.aparm[1:] = [1.,1.,2.,1.,3.,1.,0.,0.,0.,0.]
    tamrg.bparm[1:] = [1.,2.,3.,0.,0.,0.,0.,0.,0.,0.]
    tamrg.cparm[1:] = [0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]
    tamrg.dparm[1:] = [0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]
    tamrg.g
    
    # TY table
    tamrg.inext     = "TY"
    tamrg.invers    = TYver
    tamrg.outvers   = outVer
    tamrg.aparm[1:] = [1.,1.,4.,1.,5.,1.,6.,1.,0.,0.]
    tamrg.bparm[1:] = [1.,3.,4.,5.,6.,0.,0.,0.,0.,0.]
    tamrg.cparm[1:] = [timeTol/(24.*60.*60.),0.,0.,0.,0.,0.,0.,0.,0.,0.]
    tamrg.dparm[1:] = [0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]
    tamrg.g
    
    # PC table
    tamrg.inext     = "PC"
    tamrg.invers    = PCver
    tamrg.outvers   = outVer
    tamrg.aparm[1:] = [1.,1.,4.,1.,5.,1.,6.,1.,0.,0.]
    tamrg.bparm[1:] = [1.,3.,4.,5.,6.,0.,0.,0.,0.,0.]
    tamrg.cparm[1:] = [timeTol/(24.*60.*60.),0.,0.,0.,0.,0.,0.,0.,0.,0.]
    tamrg.dparm[1:] = [0.,0.,0.,0.,0.,0.,0.,0.,0.,0.]
    tamrg.g
    
    # end MergeCal
