""" Functions for calibrating and editing VLA data
"""
import UV, UVDesc, ObitTask, AIPS, OTObit, OErr
def VLAClearCal(uv, err):
    """ Clear previous calibration

    Delete all SN tables, all CL but CL 1
    uv       = UV data object to clear
    err      = Obit error/message stack
    """
    ################################################################
    uv.ZapTable("AIPS SN",-1,err)
    ver = uv.GetHighVer("AIPS CL")
    while (ver>1):
        uv.ZapTable ('AIPS CL', ver, err)
        ver = ver-1
    OErr.printErrMsg(err, "VLAClearCal: Error reseting calibration")
    

    # end VLAClearCal

def VLACal(uv, target, ACal, err, \
           PCal=None, FQid=1, calFlux=None):
    """ Basic Amplitude and phase cal for VLA data

    uv       = UV data object to calibrate
    target   = Target source name
    ACal     = Amp calibrator
    err      = Obit error/message stack
    PCal     = if given, the phase calibrator name
    FQid     = Frequency Id to process
    calFlux  = ACal flux density if given
    """
    ################################################################

    # Run SetJy
    setjy = ObitTask.ObitTask("SetJy")
    OTObit.setname(uv,setjy)
    setjy.Sources=[ACal]
    if FQid:
        setjy.FreqID=FQid
    if calFlux:
        setjy.ZeroFlux=[calFlux]
    else:
        setjy.OPType="CALC"
    setjy.g
    if PCal:
        setjy.ZeroFlux=[0.0,0.0,0.0,0.0]
        setjy.Sources=[PCal]
        setjy.OPType="REJY"
        setjy.g
    # Calib on Amp cal
    calib = ObitTask.ObitTask("Calib")
    OTObit.setname(uv,calib)
    calib.Sources = [ACal]
    calib.flagVer = 1
    calib.ampScalar = True
    calib.doCalib = 2
    calib.solMode="A&P"
    calib.solnVer = 1
    calib.g

    # ClCal CL1 + SN1 => Cl2
    clcal=ObitTask.ObitTask("CLCal")
    OTObit.setname(uv,clcal)
    clcal.solnVer = uv.GetHighVer("AIPS SN")
    clcal.calSour = [ACal]
    clcal.calIn   = uv.GetHighVer("AIPS CL")
    clcal.calOut  = clcal.calIn+1
    clcal.interMode="2PT"
    clcal.FreqID = FQid

    # Calib on phase reference if given
    if PCal:
        calib.Sources = [PCal]
        calib.flagVer = 1
        calib.ampScalar = True
        calib.doCalib = 2
        calib.solMode ="P"
        calib.g
        
        # GetJy to set flux density scale
        getjy = ObitTask.ObitTask("GetJy")
        OTObit.setname(uv,getjy)
        getjy.calSour=[ACal]
        getjy.Sources=[PCal]
        getjy.FreqID = FQid
        getjy.g

        # Set up for CLCal
        clcal.calSour = [PCal]
        
    else:  # Only amp
        
     clcal.g

    # end PCal calibration
    # end VLACal

def VLASplit(uv, target, err, FQid=1, outClass="      "):
    """ Write calibrated data

    uv       = UV data object to clear
    target   = Target source name
    err      = Obit error/message stack
    FQid     = Frequency Id to process
    """
    ################################################################
    split=ObitTask.ObitTask("Split")
    OTObit.setname(uv,split)
    split.Sources=[target]
    split.doCalib = 2
    split.gainUse = 0
    split.flagVer = 1
    split.FreqID = FQid
    split.outClass = outClass
    split.outDisk  = split.inDisk
    split.g
    # end VLAsplit

def VLASetImager (uv, target, outIclass=""):
    """ Setup to run Imager

    return Imager task interface object
    uv       = UV data object to image
    target   = Target source name
    outIclass= output class
    FQid     = Frequency Id to process
    """
    ################################################################
    img = ObitTask("Imager")
    setname(uv,img)
    img.outClass=outIClass
    img.UVTaper=[30.0, 30.0]
    img.UVRange=[0.0,50.0]
    img.FOV=0.05
    img.autoWindow=True
    img.Niter=5000
    img.Gain=0.03
    img.maxPSCloop=3
    img.minFluxPSC=0.5
    img.solPInt=10.0/60.
    img.solPType="L1"
    img.maxASCloop=1
    img.minFluxPSC=1.5
    img.solPInt=10.0
    img.minSNR=3.0
    img.avgPol=True
    img.avgIF=True
    return img
# end VLASetImager
