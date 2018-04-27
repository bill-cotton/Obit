import History, OErr
#CopyBeamHistory(x,b,err)
def CopyBeamHistory (inIm, outIm, err):
    """
    Copy beam and history from one image to another (beam)

    History may not appear in AIPS header (but is there)
    FITS History is written to "History" table
    * inIm   Input Obit image
    * outIm  Output Obit image (beam)
    * err    Obit Error/message object
    """
    # Copy Beam
    din  = inIm.Desc.Dict
    dout = outIm.Desc.Dict
    dout['beamMaj'] = din['beamMaj']
    dout['beamMin'] = din['beamMin']
    dout['beamPA']  = din['beamPA']
    outIm.Desc.Dict = dout
    outIm.UpdateDesc(err)
    # Copy History
    inHis  = History.History("in",inIm.List,err)
    outHis = History.History("out",outIm.List,err)
    History.PCopy(inHis, outHis,err)
    outIm.UpdateDesc(err)
    OErr.printErr(err)
# end CopyBeamHistory
