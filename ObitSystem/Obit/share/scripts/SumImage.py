# routine to Average planes in an image
import Image, OErr, FArray, math

def SumImageCube (inCube, outImage, err):
    """
    Average the planes in an image 

    Ignores all blank or zero images
    * inCube   = cube to sum
    * outImage = output average plane, defined but not instantiated
    * err      =  Python Obit Error/message stack
    """
    # Checks
    if not Image.PIsA(inCube):
        raise TypeError,"inCube MUST be a Python Obit Image"
    if not Image.PIsA(outImage):
        raise TypeError,"outImage MUST be a Python Obit Image"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"
    # Clone output image
    inCube.List.set("BLC",[1,1,1])
    inCube.List.set("TRC",[0,0,1])
    inCube.Clone(outImage,err)
    inCube.List.set("TRC",[0,0,0])
    OErr.printErrMsg(err, "Error initializing images")
    # Open files
    inCube.Open(Image.READONLY,err)
    outImage.Open(Image.WRITEONLY,err)
    OErr.printErrMsg(err, "Error opening images")
    nplane = inCube.Desc.Dict["inaxes"][2]
    count = 0
    for i in range (1,nplane+1):
        plane = [i,1,1,1,1]
        inCube.GetPlane(None, plane, err)
        OErr.printErrMsg(err, "Error reading image")
        # Anything here
        rms = inCube.FArray.RawRMS
        if (not math.isnan(rms)) and (rms>0.0):
            count += 1
            FArray.PAdd(inCube.FArray, outImage.FArray, outImage.FArray)
    # end loop
    norm = 1.0 / float(count)
    FArray.PSMul(outImage.FArray, norm)
    plane = [1,1,1,1,1]
    outImage.PutPlane(None, plane, err)
    OErr.printErrMsg(err, "Error writing image")
    inCube.Close(err)
    outImage.Close(err)
    OErr.printErrMsg(err, "Error closing images")
# end SumImageCube

