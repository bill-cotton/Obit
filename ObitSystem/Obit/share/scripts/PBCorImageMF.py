# Primary beam correct an ImageMF 
import Image, FArray, ImageMF, ImageUtil, FArray, OSystem, OErr
# filename = 'Abell_194.fits'
# inMF = ImageMF.newPFImageMF('in', filename, 0, True, err)
# OSystem.PAllowThreads(24)
# exec(open('PBCorImageMF.py').read()) 
MFPBCor = None; del MFPBCor
def MFPBCor (inIm, err, antSize=25, minGain=0.05):
    """
    Apply primary beam corrections to an ImageMF

    WARNING: This routine modifies the input image;
    ONLY RUN IT ONCE.
    * inIm    = Image to be modified
    * err     = Python Obit Error/message stack
    * antSize = Antenna diameter in m.
    * minGain = Minimum antenna gain
    """
    # Check
    if not inIm.ImageIsA():
        raise TypeError("input MUST be a Python Obit Image")
    # Image info
    nterm = inIm.Desc.List.Dict['NTERM'][2][0]
    nspec = inIm.Desc.List.Dict['NSPEC'][2][0]
    freqs = []
    for i in range(1,nspec+1):
        key = 'FREQ%4.4d'%i
        freqs.append(inIm.Desc.List.Dict[key][2][0])
    
    # end loop
    # Make scratch image for beam
    beam = Image.Image("PBeam")
    Image.PCloneMem(inIm, beam, err)
    OErr.printErrMsg(err, "Error with scratch beam image")
    # Loop over planes
    for i in range(1,nspec+1):
        # Read plane
        plane=[i+nterm,1,1,1,1]
        Image.PGetPlane (inIm, None, plane, err)
        OErr.printErrMsg(err, "Error reading image")
        # Set frequency for PBCor
        d = inIm.Desc.Dict; d['crval'][2] = freqs[i-1]; inIm.Desc.Dict = d
        # Make PB Image
        ImageUtil.PPBImage(beam, beam, err, minGain=minGain, antSize=antSize, outPlane=plane)
        OErr.printErrMsg(err, "Error making PB image")
        # Divide
        FArray.PDivClip (inIm.FArray, beam.FArray, minGain, inIm.FArray)
        # Rewrite plane
        Image.PPutPlane (inIm, None, plane, err)
        OErr.printErrMsg(err, "Error writing image")
   
    # end loop

# End MFPBCor
