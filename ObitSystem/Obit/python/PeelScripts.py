""" 
Routines for peeling a source at a specified position
"""
# $Id$
#-----------------------------------------------------------------------
#  Copyright (C) 2017-2021
#  Associated Universities, Inc. Washington DC, USA.
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License as
#  published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public
#  License along with this program; if not, write to the Free
#  Software Foundation, Inc., 675 Massachusetts Ave, Cambridge,
#  MA 02139, USA.
#
#  Correspondence concerning this software should be addressed as follows:
#         Internet email: bcotton@nrao.edu.
#         Postal address: William Cotton
#                         National Radio Astronomy Observatory
#                         520 Edgemont Road
#                         Charlottesville, VA 22903-2475 USA
#-----------------------------------------------------------------------

from __future__ import absolute_import
import Obit, Image, ImageDesc, SkyGeom, Table, History, OErr
import UV, UVDesc, OSystem, UVSelfCal, FArray
from ObitTask import ObitTask
from OTObit import setname, set2name, setoname
from math import cos, radians
from six.moves import range

def SelectCC(im, inCC, outCC, radius, peelPos, err):
    """
    Select/copy CCs more than radius from  peelPos
    
    This generates a CC table which can be subtracted from the uv data
    and remove all sources but the peel source area.
    * im     = Python Image with CC Tables
    * inCC   = input CC version
    * outCC  = output CC version
    * radius = radius (deg) of zone of exclusion
    * peelPos= [RA, Dec] in deg.
    * err    = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not Image.PIsA(im):
        raise TypeError("im MUST be a Python Obit Image")
    if not OErr.OErrIsA(err):
        raise TypeError("err MUST be an OErr")
    #
    # Geometry
    xref = im.Desc.Dict['crval'][0]
    yref = im.Desc.Dict['crval'][1]
    xrefpix = im.Desc.Dict['crpix'][0]
    yrefpix = im.Desc.Dict['crpix'][1]
    xinc = abs(im.Desc.Dict['cdelt'][0])
    yinc = im.Desc.Dict['cdelt'][1] 
    rot  = im.Desc.Dict['crota'][1] 
    imtype = im.Desc.Dict['ctype'][0][4:]
    # Input CC
    inTab = im.NewTable(Table.READONLY, "AIPS CC", inCC, err)
    if err.isErr:
        OErr.printErrMsg(err, "Error finding input CC Table")
        return
    # Output CC
    nrow = inTab.Desc.Dict['nrow']
    noParms = inTab.Desc.List.Dict['NO_PARMS'][2][0]
    outTab = im.NewTable(Table.WRITEONLY, "AIPS CC", outCC, err, \
                             noParms = noParms)
    if err.isErr:
        OErr.printErrMsg(err, "Error creating output CC Table")
        return
    # Open
    inTab.Open(Table.READONLY, err)
    outTab.Open(Table.WRITEONLY, err)
    if err.isErr:
        OErr.printErrMsg(err, "Error opening CC Tables")
        return
    orow = 1; count = 0; sumf = 0.0
    OErr.PLog(err, OErr.Info, "Excluding:")
    for irow in range(1,nrow+1):
        row = inTab.ReadRow(irow, err)
        # Want this one?
        dx = row['DELTAX'][0]; dy = row['DELTAY'][0]; 
        [ierr,xpos,ypos] = SkyGeom.PWorldPosLM(dx, dy, xref, yref, xinc, yinc, rot, imtype)
        # Small angle approximation
        dra = (xpos-peelPos[0])*cos(radians(xpos))
        delta = ((dra)**2+(ypos-peelPos[1])**2)**0.5
        if delta>radius:
            outTab.WriteRow(orow, row, err); orow += 1
        else:
            #print irow,xpos,ypos
            count += 1; sumf += row['FLUX'][0]
            ras = ImageDesc.PRA2HMS(xpos); decs = ImageDesc.PDec2DMS(ypos)
            OErr.PLog(err, OErr.Info, "%6d %s %s flux= %f"%(irow,ras, decs, row['FLUX'][0]))
    # End loop
    OErr.PLog(err, OErr.Info, "Drop %6d CCs, sum flux= %f"%(count,sumf))
    OErr.printErr(err)
    inTab.Close(err); outTab.Close(err)
    if err.isErr:
        OErr.printErrMsg(err, "Error copying CC Table")
        return
    # end SelectCC

def inSelectCC(im, inCC, outCC, radius, peelPos, err):
    """
    Select/copy CCs within radius from  peelPos
    
    This generates a CC table which can be subtracted from the uv data
    and remove all sources from the peel source area.
    * im     = Python Image with CC Tables
    * inCC   = input CC version
    * outCC  = output CC version
    * radius = radius (deg) of zone of inclusion
    * peelPos= [RA, Dec] in deg.
    * err    = Python Obit Error/message stack
    """
    ################################################################
    # Checks
    if not Image.PIsA(im):
        raise TypeError("im MUST be a Python Obit Image")
    if not OErr.OErrIsA(err):
        raise TypeError("err MUST be an OErr")
    #
    # Geometry
    xref = im.Desc.Dict['crval'][0]
    yref = im.Desc.Dict['crval'][1]
    xrefpix = im.Desc.Dict['crpix'][0]
    yrefpix = im.Desc.Dict['crpix'][1]
    xinc = abs(im.Desc.Dict['cdelt'][0])
    yinc = im.Desc.Dict['cdelt'][1] 
    rot  = im.Desc.Dict['crota'][1] 
    imtype = im.Desc.Dict['ctype'][0][4:]
    # Input CC
    inTab = im.NewTable(Table.READONLY, "AIPS CC", inCC, err)
    if err.isErr:
        OErr.printErrMsg(err, "Error finding input CC Table")
        return
    # Output CC
    nrow = inTab.Desc.Dict['nrow']
    noParms = inTab.Desc.List.Dict['NO_PARMS'][2][0]
    outTab = im.NewTable(Table.WRITEONLY, "AIPS CC", outCC, err, \
                             noParms = noParms)
    if err.isErr:
        OErr.printErrMsg(err, "Error creating output CC Table")
        return
    # Open
    inTab.Open(Table.READONLY, err)
    outTab.Open(Table.WRITEONLY, err)
    if err.isErr:
        OErr.printErrMsg(err, "Error opening CC Tables")
        return
    orow = 1; count = 0; sumf = 0.0
    OErr.PLog(err, OErr.Info, "Excluding:")
    for irow in range(1,nrow+1):
        row = inTab.ReadRow(irow, err)
        # Want this one?
        dx = row['DELTAX'][0]; dy = row['DELTAY'][0]; 
        [ierr,xpos,ypos] = SkyGeom.PWorldPosLM(dx, dy, xref, yref, xinc, yinc, rot, imtype)
        # Small angle approximation
        dra = (xpos-peelPos[0])*cos(radians(xpos))
        delta = ((dra)**2+(ypos-peelPos[1])**2)**0.5
        if delta<=radius:
            outTab.WriteRow(orow, row, err); orow += 1
            count += 1; sumf += row['FLUX'][0]
            ras = ImageDesc.PRA2HMS(xpos); decs = ImageDesc.PDec2DMS(ypos)
            OErr.PLog(err, OErr.Info, "%6d %s %s flux= %f"%(irow,ras, decs, row['FLUX'][0]))
        else:
            pass
            #print irow,xpos,ypos
    # End loop
    OErr.PLog(err, OErr.Info, "Include %6d CCs, sum flux= %f"%(count,sumf))
    OErr.printErr(err)
    inTab.Close(err); outTab.Close(err)
    if err.isErr:
        OErr.printErrMsg(err, "Error copying CC Table")
        return
    # end inSelectCC

def UVSub4Peel(uv, source, im, inCC, err, \
                   nfield=1, doGPU=False, gainUse=1, flagVer=-1, \
                   nThreads=1, noScrat=[0,0,0], taskLog='', debug=False):
    """
    Sets up for subtraction of non peel sources from data
    
    UV data should be have calibration tables from self calibration
    Output data will be on the same disk as the input, seq=1, class='4Peel' and
    with the name of the source (up to 12 char)
    Returns UVSub task object
    * uv        Dataset to be subtracted from
    * source    source name
    * im        Python Image with CC Tables
    * inCC      input CC version, should have had the peel source CCs removed
                using SelectCC
    * err       Python Obit Error/message stack
    * doGPU     Use GPU if available?
    * nfield    Number of facet images
    * gainUse   CL (SN) table to apply, -1=> no cal
    * flagVer   FG table to apply, -1=> no flag
    * noThreads number of threads to use
    * noScrat   AIPS disks not to use for scratch
    * taskLog   Log file
    * debug     If True leave debug Input file in /tmp
    """
    ################################################################
    # Checks
    if not UV.PIsA(uv):
        raise TypeError("uv MUST be a Python Obit UV data")
    if not Image.PIsA(im):
        raise TypeError("im MUST be a Python Obit Image")
    if not OErr.OErrIsA(err):
        raise TypeError("err MUST be an OErr")
    uvsub = ObitTask('UVSub')
    setname(uv,uvsub)
    set2name(im,uvsub)
    uvsub.CCVer = inCC; uvsub.nfield = nfield
    if uv.FileType=='AIPS':
        uvsub.DataType = 'AIPS'
    else:
        uvsub.DataType = 'FITS'
    if im.FileType=='AIPS':
        uvsub.DataType2 = 'AIPS'
    else:
        uvsub.DataType2 = 'FITS'
    uvsub.outDisk = uvsub.inDisk; uvsub.outSeq=1; uvsub.outClass='4Peel'
    uvsub.outName = source[0:12]
    uvsub.nThreads = nThreads; uvsub.noScrat = noScrat
    uvsub.Cmethod = 'DFT'; uvsub.noNeg = False; uvsub.doGPU = doGPU
    uvsub.PBCor = False; uvsub.taskLog = taskLog; uvsub.debug = debug
    uvsub.flagVer = flagVer; uvsub.Sources[0] = source
    if gainUse>=0:
        uvsub.doCalib = 2; uvsub.gainUse = gainUse
    return uvsub
    # end UVSub4Peel

def ImagePeel(uvsub, peelPos, err, \
                  nxy=512, Niter=1000, minFlux=0.001, \
                  maxPSCLoop=2, minFluxPSC=0.01, solPInt=1.0, \
                  minSNR = 3.5, Robust=0.0, doGPU=False, seq=1, \
                  nThreads=1, noScrat=[0,0,0], taskLog='', debug=False):
    """
    Sets up to image subtracted uv data from UVSub4Peel, self calibrate peel source.
    
    Only does A&P self cal
    Returns MFImage task object, output image "IPlMod', uv 'UVPeel', Seq seq
    * uvsub     task object from  UVSub4Peel
    * peelPos   [RA, Dec] in deg of source to peel
    * err       Python Obit Error/message stack
    * nxy       Size in pixels of x,y
    * Niter     max number of iterations
    * minFlux   Min flux density first CLEAN
    * maxPSCLoop max number A&P self cal loops
    * minFluxPSC min peak for self cal
    * solPInt    solution interval for self cal
    * minSNR     min SNR of self cal solutions
    * Robust     Briggs Robust factor
    * seq        Sequence number for output
    * doGPU      Use GPU if available?
    * nThreads  number of threads to use
    * noScrat   AIPS disks not to use for scratch
    * taskLog   Log file
    * debug     If True leave debug Input file in /tmp
    """
    ################################################################
    # Checks
    if not OErr.OErrIsA(err):
        raise TypeError("err MUST be an OErr")
    mfp = ObitTask('MFImage')
    mfp.DataType = 'AIPS'; mfp.inName = uvsub.outName
    mfp.inDisk = uvsub.outDisk;  mfp.inSeq = uvsub.outSeq; 
    mfp.inClass = uvsub.outClass;
    mfp.outDisk = mfp.inDisk; mfp.outSeq=seq;   mfp.outClass='IPlMod'
    mfp.out2Disk = mfp.inDisk; mfp.out2Seq=seq; mfp.out2Class='UVPeel'
    # Shift
    uv = UV.newPAUV('uv',mfp.inName,mfp.inClass,mfp.inDisk,mfp.inSeq, True,  err)
    d = uv.Desc.Dict; ra = d['crval'][4];  dec = d['crval'][5]
    shift = SkyGeom.PShiftXY(ra, dec, d['crota'][1], peelPos[0], peelPos[1])
    mfp.RAShift[0] = 3600*shift[0]; mfp.DecShift[0] = 3600*shift[1]; 
    mfp.nThreads = nThreads; mfp.noScrat = noScrat
    mfp.Cmethod = 'DFT'; mfp.noNeg = False; mfp.doGPU = doGPU
    mfp.PBCor = False; mfp.taskLog = taskLog; mfp.debug = debug
    mfp.prtLv = 2; mfp.PBCor = False; 
    mfp.NField = 1; mfp.nx[0] = nxy; mfp.ny[0] = nxy; 
    mfp.OutlierDist=0.0; mfp.minFlux = minFlux
    mfp.Niter = Niter; mfp.Robust = Robust; mfp.minSNR = minSNR
    mfp.maxPSCLoop = maxPSCLoop; mfp.minFluxPSC = minFluxPSC
    mfp.solPInt = solPInt; mfp.solPType='L1'; mfp.solPMode = 'A&P'
    mfp.maxASCLoop = 0; mfp.minFluxASC = minFluxPSC
    mfp.solAInt = solPInt*2; mfp.solAType='L1'; mfp.solAMode = 'A&P'
    return mfp
    # end ImagePeel

def SubPeel(uv, source, imp, uvp, err, addBack=False, seq=999, \
                flagVer=0, nThreads=1, doGPU=False, noScrat=[0,0,0], taskLog='', debug=False):
    """
    Subtract Peel model w/ solutions, then optionally add back w/o corruptions
    
    UV data should have calibration tables from self calibration
    Output data will be on the same disk as the input, seq=seq, class='PelSub' and
    with name = source (up to 12 char).
    Returns Peel source subtracted/replaced data
    * uv        Dataset with cal tables
                Needs at least the self cal gain table
    * source    source name
    * imp       Peel source model (CC table from ImagePeel)
    * uvp       UV data the result of peel (ImagePeel)
    * err       Python Obit Error/message stack
    * seq       Sequence number for output
    * addBack   Add model back to data w/o corruptions? Not recommended.
    * flagVer   FG table to apply, -1=> no flag
    * nThreads  number of threads to use
    * doGPU     Use GPU if available?
    * noScrat   AIPS disks not to use for scratch
    * taskLog   Log file
    * debug     If True leave debug Input file in /tmp
    """
    ################################################################
    # Checks
    if not UV.PIsA(uv):
        raise TypeError("uv MUST be a Python Obit UV data")
    if not Image.PIsA(imp):
        raise TypeError("imp MUST be a Python Obit Image")
    if not OErr.OErrIsA(err):
        raise TypeError("err MUST be an OErr")
    # Split main data set
    OErr.PLog(err, OErr.Info, "Copy data"); OErr.printErr(err)
    split = ObitTask('Split'); setname(uv,split)
    split.outDisk = split.inDisk; split.outSeq = seq; split.outClass = 'UPeel'
    split.Sources[0] = source; split.flagVer = flagVer
    if uv.GetHighVer('AIPS SN')>0:
        split.doCalib = 2; split.gainUse = 0; split.doCalWt = True
    else:
        split.doCalib = -1
    split.taskLog = taskLog; split.debug = debug
    split.g
    outClass = split.outClass; outDisk = split.outDisk; outSeq = split.outSeq

    # Get data
    OErr.PLog(err, OErr.Info, "Make Peel model with corruptions"); OErr.printErr(err)
    if UV.AExist(source[0:12],outClass,outDisk,outSeq,  err):
        datauv = UV.newPAUV('data',source[0:12],outClass,outDisk,outSeq, True,  err)
    else:
        datauv = UV.newPAUV('data',source[0:8],outClass,outDisk,outSeq, True,  err)
    # Make data set with the model peel source with peel cal applied
    uvsub = ObitTask('UVSub'); setname(uv,uvsub); uvsub.outName = source[0:12]
    uvsub.outDisk = uvsub.inDisk; uvsub.outSeq = 1; uvsub.outClass = 'Model'
    uvsub.Sources[0] = source; uvsub.flagVer = flagVer
    uvsub.doCalib = -1; uvsub.gainUse = 0; 
    set2name(imp,uvsub); uvsub.CCVer=1; uvsub.nfield = 1;
    uvsub.Cmethod = 'DFT'; uvsub.Opcode = 'MODL'; uvsub.PBCor = False;
    uvsub.noScrat = noScrat; uvsub.noNeg = False;
    uvsub.taskLog = taskLog; uvsub.nThreads=nThreads; uvsub.doGPU=doGPU; uvsub.debug = debug
    uvsub.g
 
    # Get model data
    modeluv = UV.newPAUV('model',uvsub.outName,uvsub.outClass,uvsub.outDisk,uvsub.outSeq,True, err)
    # Copy/invert/unblank SN table from peeluv
    hiPeelSN = uvp.GetHighVer('AIPS SN')
    inTab   = max(1, hiPeelSN)
    sntab = uvp.NewTable(Table.READONLY, 'AIPS SN', inTab, err)
    z=UVSelfCal.PInvertSN(sntab, modeluv, 1, True, err)
    # Apply calibration table in subtract
    modeluv.List.set('doCalSelect', True); 
    modeluv.List.set('doCalib', 2); modeluv.List.set('gainUse', 1); 
    modeluv.List.set('passAll', True); 

    # Subtract model from main data
    OErr.PLog(err, OErr.Info, "Subtract Corrupted Peel model from uv data")
    UV.PUtilVisSub(datauv, modeluv, datauv, err)
    OErr.printErr(err)

    # Add model without corrupting calibration
    if addBack:
        OErr.PLog(err, OErr.Info, "Add Peel model without corruptions"); 
        uvsub = ObitTask('UVSub'); setname(datauv,uvsub); setoname(datauv,uvsub)
        uvsub.outSeq = uvsub.inSeq+1; 
        uvsub.Sources[0] = source; uvsub.flagVer = flagVer
        uvsub.doCalib = -1; uvsub.gainUse = 0; 
        set2name(imp,uvsub); uvsub.CCVer=1; uvsub.nfield = 1;
        uvsub.Cmethod = 'DFT'; uvsub.Factor=-1.; uvsub.PBCor = False;
        uvsub.noScrat = noScrat; uvsub.noNeg = False;
        uvsub.taskLog = taskLog; uvsub.nThreads=nThreads; uvsub.doGPU=doGPU; uvsub.debug = debug
        uvsub.g
        outClass = uvsub.outClass; outDisk = uvsub.outDisk; outSeq = uvsub.outSeq
        # end add back
    OErr.printErr(err)
    # Delete model dataset
    if not debug:
        modeluv.Zap(err)

   # final data
    if UV.AExist(source[0:12],outClass,outDisk,outSeq, err):
        datauv2 = UV.newPAUV('data', source[0:12], outClass, outDisk, outSeq, True,  err)
    else:
        datauv2 = UV.newPAUV('data', source[0:8], outClass, outDisk, outSeq, True,  err)
    return datauv2
    # end SubPeel

def RestorePeel(peelMod, CCVer, image, err):
    """
    Restore CCs from one image onto another
    
    If images are ImageMF then multiple planes restored.
    * peelMod   Image with CC table (as Image)
    * CCver     CC version on peelMod to restore
    * image     Output Image to which components to be added
    * err       Python Obit Error/message stack
    * nThreads  number of threads to use
    """
    ################################################################
    import Obit, Image, ImageDesc, SkyGeom, Table, History, OErr
    import UV, UVDesc, OSystem, UVSelfCal, FArray
    # Checks
    if not Image.PIsA(peelMod):
        raise TypeError("uv MUST be a Python Obit Image")
    if not Image.PIsA(image):
        raise TypeError("imp MUST be a Python Obit Image")
    
    ph = peelMod.Desc.Dict; ih = image.Desc.Dict  # Descriptors
    ph_x   = ph['crval'][0]; ph_y  = ph['crval'][1]; ph_rot  = ph['crota'][1];
    ph_dx  = ph['cdelt'][0]; ph_dy = ph['cdelt'][1]; ph_type = ph['ctype'][0][4:];
    ph_ref_x  = ph['crpix'][0]; ph_ref_y = ph['crpix'][1]; 
    ih_x   = ih['crval'][0]; ih_y  = ih['crval'][1]; ih_rot  = ih['crota'][1];
    ih_dx  = ih['cdelt'][0]; ih_dy = ih['cdelt'][1]; ih_type = ih['ctype'][0][4:];
    ih_ref_x  = ih['crpix'][0]; ih_ref_y = ih['crpix'][1]; 
    bmaj   = ih['beamMaj'];  bmin  = ih['beamMin']; bpa = ih['beamPA']; 
    bmaj /= abs(ih_dx); bmin /= abs(ih_dx)   # Beam in pixels (square grid)
    
    #Beam to insert
    beam =  FArray.FArray('beam',naxis=[21,21]); bx = 10.; by=10.
    FArray.PEGauss2D(beam, 1.0, [bx,by], [bmaj, bmin, bpa])
    
    cctab = peelMod.NewTable(Table.READONLY, 'AIPS CC', CCVer, err) # CC Table
    cctab.Open(Table.READONLY, err)
    OErr.printErr(err)
    
    # First plane
    image.GetPlane(None, [1,1,1,1,1],err)
    imArr = image.FArray
    ncc = cctab.Desc.Dict['nrow']  # Number of CCs
    for irow in range(1,ncc+1):
        row  = cctab.ReadRow(irow,err)
        flux = row['FLUX'][0]; 
        dx   = row['DELTAX'][0]; dy = row['DELTAY'][0]; dz = row['DELTAZ'][0];
        [ierr,ra,dec] = SkyGeom.PWorldPosLM(dx,dy,ph_x, ph_y, ph_dx, ph_dy, ph_rot, ph_type)
        # output pixel
        [ierr,xpix,ypix] = SkyGeom.PXYpix(ra, dec, ih_x, ih_y, ih_ref_x, ih_ref_y, ih_dx, ih_dy, ih_rot, ih_type)
        # Need beam if correct location
        tbx = xpix - int(xpix+0.5); tby = ypix - int(ypix+0.5);
        FArray.PFill(beam, 0.0);
        FArray.PEGauss2D(beam, 1.0, [bx+tbx,by+tby], [bmaj, bmin, bpa])
        # Add
        FArray.PShiftAdd(imArr, [int(xpix+0.5)-1,int(ypix+0.5)-1], beam, [int(bx+0.5),int(by+0.5)], flux, imArr)
    # end loop
    image.PutPlane(None, [1,1,1,1,1],err) # rewrite
    
    # MFImage Planes
    if ih['ctype'][2]=='SPECLNMF ':
        nterm = image.Desc.List.Dict['NTERM'][2][0]; nspec = image.Desc.List.Dict['NSPEC'][2][0]
        for ip in range(nterm+1,nterm+nspec+1):
            image.GetPlane(None, [ip,1,1,1,1],err)
            imArr = image.FArray
            for irow in range(1,ncc+1):
                row  = cctab.ReadRow(irow,err)
                flux = row['PARMS'][3+ip-nterm]
                dx   = row['DELTAX'][0]; dy = row['DELTAY'][0]; dz = row['DELTAZ'][0];
                [ierr,ra,dec] = SkyGeom.PWorldPosLM(dx,dy,ph_x, ph_y, ph_dx, ph_dy, ph_rot, ph_type)
                # output pixel
                [ierr,xpix,ypix] = SkyGeom.PXYpix(ra, dec, ih_x, ih_y, ih_ref_x, ih_ref_y, ih_dx, ih_dy, ih_rot, ih_type)
                # Need beam in correct location
                tbx = xpix - int(xpix+0.5); tby = ypix - int(ypix+0.5);
                FArray.PFill(beam, 0.0);
                FArray.PEGauss2D(beam, 1.0, [bx+tbx,by+tby], [bmaj, bmin, bpa])
                # Add
                FArray.PShiftAdd(imArr, [int(xpix+0.5)-1,int(ypix+0.5)-1], beam, [int(bx+0.5),int(by+0.5)], flux, imArr)
            # end loop
            image.PutPlane(None, [ip,1,1,1,1],err) # rewrite
      # end plane loop
    # end MFImage Planes
    # History
    x = image; y=peelMod
    hi = x.History(3,err)
    hi.Open(3,err)
    hi.TimeStamp('Restore Peel', err)
    if x.FileType=='FITS':
        hiCard = y.FileName.strip()+'.'+str(y.Disk)
    else:
        hiCard = y.Aname.strip()+'.'+y.Aclass.strip()+'.'+str(y.Aseq)+'.'+str(y.Disk)
    
    hi.WriteRec(-1,"RestorePeel / file="+hiCard,err)
    hi.Close(err)
# end RestorePeel
