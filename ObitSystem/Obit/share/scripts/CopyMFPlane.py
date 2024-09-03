# Copy a spectral plane from one MFImage to another

import OErr, InfoList, Image, ImageMF, ImageDesc, History
def CopyMFPlane(inFile, inPlane, outFile, outPlane, err):
    """
    Copy a spectral plane from one FITS MFImage to another

    Files should be in current working directory
    * inFile   = input MFImage file name
    * inPlane  = 1-rel Spectral plane number from input (after NTERM)
    * outFile  = extant output MFImage file name
    * outPlane = 1-rel Spectral plane number for output (after NTERM)
    * err      = Obit error/message object
    """
    # Convert files into Images
    inImage = Image.newPFImage(inFile, inFile, 0, True, err)
    OErr.printErrMsg(err, "Error initializing input image")
    inNterm = inImage.Desc.List.Dict['NTERM'][2][0]
    # output image
    outImage = Image.newPFImage(outFile, outFile, 0, True, err)
    OErr.printErrMsg(err, "Error initializing output image")
    if OErr.PIsErr(err):
        OErr.printErr(err); raise RuntimeError("File Error")
    outNterm = outImage.Desc.List.Dict['NTERM'][2][0]
  
    # Copy plane
    plane = [inPlane+inNterm,1,1,1,1]
    inImage.GetPlane(None, plane, err)
    plane = [outPlane+outNterm,1,1,1,1]
    outImage.PutPlane(inImage.FArray, plane,err)
    if OErr.PIsErr(err):
        OErr.printErr(err); raise RuntimeError("File Error")
    # Open output to update header
    outImage.Open(ImageMF.READWRITE,err)
    # Copy Frequency info
    idl = inImage.Desc.List;
    odl = outImage.Desc.List
    key = 'FREQ%4.4d'%inPlane; val = idl.get(key)[4][0]
    key = 'FREQ%4.4d'%outPlane; odl.set(key,val)
    key = 'FREH%4.4d'%inPlane; val = idl.get(key)[4][0]
    key = 'FREH%4.4d'%outPlane; odl.set(key,val)
    key = 'FREL%4.4d'%inPlane; val = idl.get(key)[4][0]
    key = 'FREL%4.4d'%outPlane; odl.set(key,val)
    # Close
    outImage.Close(err)
    
    # Add history as Image
    outHistory  = History.History("history", outImage.List, err)
    z=outHistory.Open(History.READWRITE, err)
    outHistory.TimeStamp(" Start Obit CopyMFPlane",err)
    z=outHistory.WriteRec(-1,"CopyMFPlane / copy plane "+str(inPlane)+" from "+inFile,err)
    z=outHistory.WriteRec(-1,"CopyMFPlane / to plane "+str(outPlane),err)
    z=outHistory.Close(err)
    outImage.UpdateDesc(err)
    
    # Say something
    print ("Copy plane "+str(inPlane)+" from "+inFile,"to plane "+str(outPlane))
    
# end CopyMFPlane
