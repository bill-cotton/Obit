# Clip RM and EVLA planes in RMFit cube by ppol
import Image, FArray, ImageUtil, FArray, OSystem, OErr
import History
# filename = 'J1723+6547_AK_PB.fits'
# inMF = Image.newPFImage('in', filename, 0, True, err)
# OSystem.PAllowThreads(24)
# exec(open('PBCorImageMF.py').read()) 
#RMClip (inMF, err, minPPol=50.0e-6)
RMClip = None; del RMClip
def RMClip (inIm, err, minPPol=0.0):
    """
    Clip RM and EVLA planes in RMFit cube by ppol

    WARNING: This routine modifies the input image;
    ONLY RUN IT ONCE.
    * inIm    = Image to be modified
    * err     = Python Obit Error/message stack
    * minPPol = min PPol (inIm plane 3) in Jy
    """
    # Check
    if not inIm.ImageIsA():
        raise TypeError("input MUST be a Python Obit Image")
    # Image info
    d = inIm.Desc.Dict
    (nx,ny) = d['inaxes'][0:2]
    # Read planes
    RMPix = FArray.FArray('RM', [nx,ny])
    inIm.GetPlane(RMPix, [1,1,1,1,1], err)
    OErr.printErrMsg(err, "Error reading image")
    EVPAPix = FArray.FArray('EVPA', [nx,ny])
    inIm.GetPlane(EVPAPix, [2,1,1,1,1], err)
    PPolPix = FArray.FArray('PPol', [nx,ny])
    inIm.GetPlane(PPolPix, [3,1,1,1,1], err)
    OErr.printErrMsg(err, "Error reading image")
    # Make mask in PPolPix
    FArray.PInClip(PPolPix, -1.0e6, minPPol, FArray.fblank)
    # Blank
    FArray.PBlank(RMPix,   PPolPix, RMPix)
    FArray.PBlank(EVPAPix, PPolPix, EVPAPix)
    # Rewrite
    inIm.PutPlane(RMPix,   [1,1,1,1,1], err)
    inIm.PutPlane(EVPAPix, [2,1,1,1,1], err)
    OErr.printErrMsg(err, "Error writing image")

    # Add history
    outHistory = History.History("history", inIm.List, err)
    outHistory.Open(History.READWRITE, err)
    outHistory.TimeStamp(" Start Obit RMClip",err)
    outHistory.WriteRec(-1,"RMClip"+" minPPol = "+str(minPPol),err)
    outHistory.Close(err)

# End RMClip
