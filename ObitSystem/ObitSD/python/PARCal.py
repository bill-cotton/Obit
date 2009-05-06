""" Mustang calibration Utility package
"""
# $Id$
#-----------------------------------------------------------------------
#  Copyright (C) 2008
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
import Obit, Image, OTF, OTFUtil, GBTUtil, CleanOTF, OTFGetSoln
import ImageFit, FitRegion, FitModel
import math
import OSystem, OErr, ODisplay

def input(inputDict):
    """ Print the contents of an input Dictionary

    inputDict = Python Dictionary containing the parameters for a routine
    There should be a member of the dictionary ('structure') with a value
    being a list containing:
    1) The name for which the input is intended (string)
    2) a list of tuples consisting of (parameter name, doc string)
       with an entry for each parameter in the dictionary.
       The display of the the inputs dictionary will be in the order of
       the tuples and display the doc string after the value.
       An example:
       Soln2CalInput={'structure':['Soln2Cal',[('InData','Input OTF'),
                                               ('soln','input soln table version'),
                                               ('oldCal','input cal table version, -1=none'),
                                               ('newCal','output cal table')]],
                      'InData':None, 'soln':0, 'oldCal':-1, 'newCal':0}
    """
    ################################################################
    structure = inputDict['structure']  # Structure information
    print 'Inputs for ',structure[0]
    for k,v in structure[1]:
        print '  ',k,' = ',inputDict[k],' : ',v
        
    # end input

# Define FitCal input dictionary
FitCalInput={'structure':['FitCal',[('InData','Input OTF'),
                                    ('scrDisk','FITS scrath file disk number'),
                                    ('PSF','Input Beam image for telescope'),
                                    ('scanList','list of lists of target,scans to be processed'),
                                    ('disp','Image display to edit window'),
                                    ('save','If True, keep derived images'),
                                    ('nx','number of pixels in x = RA'),
                                    ('ny','number of pixels in y = de'),
                                    ('xCells','Cell spacing in x (asec)'),
                                    ('yCells','Cell spacing in y (asec)'),
                                    ('minWt','minimum summed weight in gridded image wrt max '),
                                    ('ConvType','Conv. fn type, 0=pillbox,3=Gaussian,4=exp*sinc,5=Sph wave'),
                                    ('ConvParm','Conv. fn parameters'),
                                    ('gainUse','cal. table version, -1=none'),
                                    ('flagVer','flag table version, -1=none'),
                                    ('Niter','Maximum number of CLEAN iterations def[100]'),
                                    ('Patch','Beam patch in pixels [def 100]'),
                                    ('BeamSize','Restoring beam FWHM (deg)'),
                                    ('Gain','CLEAN loop gain def [0.1]'),
                                    ('autoWindow','Automatically set Windows? [False]'),
                                    ('doFilter','Filter out of band noise? [True]'),
                                    ('solType','Soln type, Common, Detector, Both'),
                                    ('solInt','Min Solution interval(sec)'),
                                    ('soln','Array of multiples of solInt for cal')]],
             'InData':None, 'scrDisk':1, 'PSF':None, "scanList":[["All",1]],
             'nx':100, 'ny':100, 'xCells':2.0, 'yCells':2.0, 'minWt':0.01,
             'ConvParm':[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0], 'ConvType':5,
             'gainUse':-1, 'flagVer':-1,'disp':None, 'save':False, 'doFilter':True,
             'Niter':100, 'Patch':100, 'BeamSize':0.0, 'Gain':0.1, 'autoWindow':False,
             'solType':"Both", 'solInt':1.0, 'soln':[3.,1.] }
def FitCal (err, input=FitCalInput):
    """ Image and fit a sequence of calibrator scans

    Does imaging and iterative calibration of a sequence of scans on a calibrator
    (Target).
    Calibration is by imaging  and performing a CLEAN to derive a sky model
    and then subtracting the sky model from the data.  An estimate of the
    residual background is derived from filtering the residuals.
    The number of calibration cycles is determined by the number of entries in
    soln and each cycle may consist of a "Common" calibration, or a "Detector"
    calibration or both.
    Data are imaged as grouped in scanList, scans in each list in scanList
    are imaged together
    Returns an array of dict's, one per list in scanList with
        "Time"   Center time in Days
        "Target" Source Name
        "Peak"   Peak flux density
        "Gauss"  Gaussian parameters (major, minor axis size (asec), PA (deg))
        "elOff"  Offset in elevation in asec
        "azOff"  Offset in azimuth ( actually Xel)
    err     = Python Obit Error/message stack
    input   = input parameter dictionary, for interactive use, the function input
              will display the contents in human readable format.

    Input dictionary entries:
    InData  = Python input OTF to calibrate
    scrDisk = Disk number for scratch files
    PSF     = Image with telescope psf
    scanList= list of lists of target/scans to be imaged together
              list of forn [target, scan_1...]
    disp    = Image display to show final "dirty" maps
    save    = if True, save derived images, else delete
              These will be in scrDisk with names of the form
              target.+scan+.CalImage.fits
    nx      = number of pixels in "x" = RA
    ny      = number of pixels in 'Y' = dec
    xCells  = Cell spacing in x (asec)
    yCells  = Cell spacing in y (asec)
    ConvType= Convolving function Type 0=pillbox,3=Gaussian,4=exp*sinc,5=Sph wave
    ConvParm= Convolving function parameters depends on ConvType
      Type 2 = Sinc, (poor function - don't use)
        Parm[0] = halfwidth in cells,
        Parm[1] = Expansion factor
      Type 3 = Gaussian,
        Parm[0] = halfwidth in cells,[def 3.0]
        Parm[1] = Gaussian with as fraction or raw beam [def 1.0]
      Type 4 = Exp*Sinc
        Parm[0] = halfwidth in cells, [def 2.0]
        Parm[1] = 1/sinc factor (cells) [def 1.55]
        Parm[2] = 1/exp factor (cells) [def 2.52]
        Parm[3] = exp power [def 2.0]
      Type 5 = Spherodial wave 
        Parm[0] = halfwidth in cells [def 3.0]
        Parm[1] = Alpha [def 5.0]
        Parm[2] = Expansion factor [not used]
    gainUse = version number of prior table (Soln or Cal) to apply, -1 is none
    flagVer = version number of flagging table to apply, -1 is none
    Niter       = Maximum number of CLEAN iterations
    Patch       = Beam patch in pixels [def 100]
    BeamSize    = Restoring beam (deg)
    Gain        = CLEAN loop gain
    autoWindow  = True if autoWindow feature desired
    doFilter    = Filter out of band noise?
    solInt      = solution interval (sec)
    solType     = solution type:
       "Common"   solve for common mode.additive effects on timescales longer
               than solInt
       "Detector" Solve for detector additive terms on timescales longer
               than 3 * solInt
       "Both" Both Common and Detector solutions each calibration cycle
    soln        = list of solution intervals, one cycle per interval
                  NO value should be less than solInt
    """
    ################################################################
    # Get input parameters
    dim = [1,1,1,1,1]
    inData  = input["InData"]
    inInfo = inData.List
   # Checks
    if not OTF.PIsA(inData):
        raise TypeError,'Image: Bad input OTF'
    if err.isErr: # existing error?
        return None

    # Set calibration 
    gainUse = input["gainUse"]
    flagVer = input["flagVer"]
    # Default table versions (the Obit routines will do this as well)
    if gainUse == 0:   # Get highest numbered OTFCal table
        gainUse = Obit.OTFGetHighVer(inData.me, "OTFCal")
    if gainUse == 0:   # Doesn't seem to be one, try OTFSoln
        gainUse = Obit.OTFGetHighVer(inData.me, "OTFSoln")
    if gainUse == 0:   # Must not be one
        gainUse = -1
    inInfo.set("doCalSelect", True)       
    inInfo.set("flagVer", flagVer)
    if gainUse>0:
        inInfo.set("doCalib", 1)       
        inInfo.set("gainUse", gainUse)       
        
    # Set imaging/calibration parameters
    scrDisk = input["scrDisk"]
    # Imaging parameters
    OTF.ImageInput["disk"]    = scrDisk 
    OTF.ImageInput["Beam"]    = input["PSF"]
    OTF.ImageInput["xCells"]  = input["xCells"]
    OTF.ImageInput["yCells"]  = input["yCells"]
    OTF.ImageInput["nx"]      = input["nx"]
    OTF.ImageInput["ny"]      = input["ny"]
    OTF.ImageInput["gainUse"] = 0
    OTF.ImageInput["flagVer"] = flagVer
    OTF.ImageInput["minWt"]   = input["minWt"]
    OTF.ImageInput["ConvType"]  = input["ConvType"]
    OTF.ImageInput["ConvParm"]  = input["ConvParm"]
    OTF.ImageInput["doFilter"]  = input["doFilter"]
    
    # Calibration parameters (some reset in loop)
    OTF.ResidCalInput["solType"] = "MultiBeam"
    OTF.ResidCalInput["solInt"]  = input["solInt"]
    OTF.ResidCalInput["gainUse"] = 0
    OTF.ResidCalInput["minEl"]   = -90.0 
    OTF.ResidCalInput["flagVer"] = flagVer
    OTF.ResidCalInput["minFlux"] = 0.0
    
    # CLEAN parameters
    CleanOTF.CleanInput["autoWindow"]= input["autoWindow"]
    CleanOTF.CleanInput["Patch"]    = input["Patch"]
    CleanOTF.CleanInput["Niter"]    = input["Niter"]
    CleanOTF.CleanInput["Gain"]     = input["Gain"]
    CleanOTF.CleanInput["BeamSize"] = input["BeamSize"]
    #DEBUG CleanOTF.CleanInput["disp"]     = input["disp"]
    CleanOTF.CleanInput["minFlux"]  = 0.0
    
    # Reset Soln2Cal parameters for self cal
    OTF.Soln2CalInput["oldCal"]  = 1    
    OTF.Soln2CalInput["newCal"]  = 2

    disp    = input["disp"]  # Image display

    # Initialize output
    results = []
    # Loop over images
    scount = 0;
    for scans in input["scanList"]:
        # Make scratch OTF data
        scrData  = inData.Scratch(err)
        scrInfo  = scrData.List
        
        target = scans[0]
        pos =  GBTUtil.GetTargetPos(inData, target, err)
        ra  = pos[0]                      # ra of center
        dec = pos[1]                      # dec of center
        OTF.ImageInput["ra"]      = ra 
        OTF.ImageInput["dec"]     = dec
        OTF.ImageInput["OutName"] = target+"."+str(scans[1])+".CalImage.fits"
 
        # copy first scan
        inInfo.set("Scans",[scans[1],scans[1]])
        inInfo.set("Targets",target)
        inData.Copy(scrData, err)
        # Concatenate rest
        OTF.ConcatInput["InData"]  = inData
        OTF.ConcatInput["OutData"] = scrData
        for scan in scans[2:]:
            inInfo.set("Targets",target)
            inInfo.set("Scans",[scan,scan])
            OTF.Concat(err, OTF.ConcatInput)
        if err.isErr:
            OErr.printErrMsg(err, "Error creating scratch OTF data")

        # delete any prior calibration tables
        OTF.ClearCal(scrData,err)

       # Dummy cal table
        inter = input["solInt"]/4
        OTFGetSoln.POTFGetDummyCal (scrData, scrData, inter, 1, 1, err)
        if err.isErr:
            OErr.printErrMsg(err, "Error creating dummy cal table")
        scrInfo.set("gainUse",0)
        scrInfo.set("doCalSelect", True)       
        OTF.ImageInput["InData"]     = scrData
        OTF.ResidCalInput["InData"]  = scrData
        OTF.Soln2CalInput["InData"]  = scrData
        OTF.ImageInput["gainUse"]    = 0

        # Create dirty image
        OTF.ImageInput["InData"]  = scrData
        DirtyImg = OTF.makeImage(err, OTF.ImageInput)
        OErr.printErrMsg(err, "Error making initial image")

        # Image for cleaning
        CleanImg = DirtyImg.Scratch(err)

        # Create CleanOTF
        dirtyBeam = input["PSF"]
        CleanObj = CleanOTF.PCreate("Clean", DirtyImg, dirtyBeam, CleanImg, err)
        OErr.printErrMsg(err, "Error creating CLEAN object")
        CleanOTF.CleanInput["CleanOTF"] = CleanObj     # Clean object

        #  CLEAN window
        win = [-1, 20, input["nx"]/2+1, input["ny"]/2+1]
        CleanOTF.PAddWindow(CleanObj, win, err)
        
        # Calibration loop
        count=0
        soln = input["soln"]
        OTF.ImageInput["gainUse"]    = 1     # Which cal table to apply
        OTF.ResidCalInput["gainUse"] = 0     # Prior calibration,
        scrInfo.set("deMode", False)  # Don't remove mode first cycle
        if input["solType"]=="Both":
            for si in soln:
                count = count+1
                # First calibration of pair
                OTF.ResidCalInput["solInt"]  = 3*si
                OTF.ResidCalInput["solType"] = "Offset"
                OTF.Soln2CalInput["oldCal"]  = 1
                OTF.Soln2CalInput["newCal"]  = 2           # (Re)Use 2
                OTF.SelfCal(err, OTF.ImageInput, CleanOTF.CleanInput, OTF.ResidCalInput, OTF.Soln2CalInput)
                OErr.printErrMsg(err, "Error in self cal")
                OTF.ImageInput["gainUse"]    = OTF.Soln2CalInput["newCal"]   # Which cal table to apply
                OTF.ResidCalInput["gainUse"] = OTF.Soln2CalInput["newCal"]   # Prior calibration
                # Second
                OTF.ResidCalInput["solInt"]  = si
                OTF.ResidCalInput["solType"] = "MultiBeam"
                OTF.Soln2CalInput["oldCal"]  = 2
                OTF.Soln2CalInput["newCal"]  = 3           # (Re)Use 3
                OTF.SelfCal(err, OTF.ImageInput, CleanOTF.CleanInput, OTF.ResidCalInput, OTF.Soln2CalInput)
                OTF.ImageInput["gainUse"]    = OTF.Soln2CalInput["newCal"]   # Which cal table to apply
                OTF.ResidCalInput["gainUse"] = 1   # Prior calibration for next cycle
                # Cleanup Soln tables
                scrData.ZapTable("OTFSoln",1,err)
                scrData.ZapTable("OTFSoln",2,err)
                scrInfo.set("deMode", True) # remove mode
        else:   # Only one type
            if input["solType"]=="Common":
                stype = "MultiBeam"
            else:
                stype = "Offset"                
            for si in soln:
                count = count+1
                OTF.ResidCalInput["solInt"]  = si
                OTF.ResidCalInput["solType"] = stype
                OTF.Soln2CalInput["oldCal"]  = 1
                OTF.Soln2CalInput["newCal"]  = 2
                OTF.SelfCal(err, OTF.ImageInput, CleanOTF.CleanInput, OTF.ResidCalInput, OTF.Soln2CalInput)
                OTF.ImageInput["gainUse"]    = OTF.Soln2CalInput["newCal"]   # Which cal table to apply
                OTF.ResidCalInput["gainUse"] = 1   # Prior calibration for next cycle
                # Cleanup Soln tables
                scrData.ZapTable("OTFSoln",1,err)
                scrInfo.set("deMode", True) # remove mode
        # Final image
        DirtyImg = OTF.makeImage(err,  OTF.ImageInput)
        OErr.printErrMsg(err, "Error in final dirty image")

        print "Display final image"
        ODisplay.PImage(disp, DirtyImg, err)
        OErr.printErrMsg(err, "Error displaying final image")
        
        # Do fitting
        # Model fitting
        import ImageFit, FitRegion, FitModel
        
        # Set up model
        fm = FitModel.FitModel(type=FitModel.GaussMod, Peak=1.0, parms=[3.,3.,0.], \
                               DeltaX=8., DeltaY=8.)
        cen=[8,8]
        DirtyImg.Open(Image.READONLY,err)
        DirtyImg.ReadPlane(err)
        fm.Peak = FArray.PMax(DirtyImg.FArray,cen)
        DirtyImg.Close(err)
        corner=[cen[0]-8, cen[1]-8]
        cen=[8,8]
        d = DirtyImg.Desc.Dict
        fr = FitRegion.FitRegion(corner=corner,dim=[17,17],models=[fm])
        
        imf = ImageFit.ImageFit("ImageFit")
        fitinput = ImageFit.FitInput
        fitinput["fitImage"]  = DirtyImg 
        fitinput["fitRegion"] = fr
        fitinput["MaxIter"]   = 100
        fitinput["prtLv"]     = 0
        fitinput["PosGuard"]  = 1.
        print "Initial guess Corner=",corner,"Peak=",fm.Peak
        # Fit
        imf.Fit(err, fitinput)

        # Offsets
        raOff  = (fm.DeltaX+corner[0]-d["crpix"][0])*d["cdelt"][0]
        decOff = (fm.DeltaY+corner[1]-d["crpix"][1])*d["cdelt"][1]

        print "Processing target/scans",scans
        print "Results Peak=",fm.Peak,"dx=",raOff*3600.0,"dy=",decOff*3600.0
        print "Results Center=[",fm.DeltaX+corner[0],",",fm.DeltaY+corner[1],"]"
        print "Results Gauss=",fm.parms

        # Get time, pa of data
        scrInfo.set("doCalSelect", False)
        scrInfo.set("doCalib", -1)
        scrData.Open(OTF.READONLY,err)
        nrow   = scrData.Desc.Dict["nrecord"]
        row    = scrData.ReadRec(err)
        ttime  = row.time
        pa     = row.rot
        scrData.Close(err)
        del row

        print "target",target,"time",ttime*24.0,"par ang",pa
        pa /= 57.296

        # Offsets
        raOff  = (fm.DeltaX+corner[0]-d["crpix"][0])*d["cdelt"][0]
        decOff = (fm.DeltaY+corner[1]-d["crpix"][1])*d["cdelt"][1]

        # Gaussian parameters to asec
        gparm = [fm.parms[0]*abs(d["cdelt"][0])*3600.0, \
                 fm.parms[1]*abs(d["cdelt"][1])*3600.0,
                 fm.parms[2]*57.296]

        # Rotate to az, el
        azOff = raOff*math.cos(pa)  - decOff*math.sin(pa)
        elOff = decOff*math.cos(pa) + raOff*math.sin(pa)
        # to asec
        azOff *= 3600.0
        elOff *= 3600.0

        # Set output
        fitdata = {"Time":ttime, "Peak":fm.Peak, "Gauss":gparm, \
                   "elOff":elOff, "azOff":azOff,"Target":target}
        results.append(fitdata)
    
        # Cleanup
        scrData.Zap(err)
        if not  input["save"]:
            DirtyImg.Zap(err)
        else:
            del DirtyImg
            CleanImg.Zap(err)
        del CleanObj
        scount += 1
        # end loop over images

    return results
    # end FitCal

def InitCal (inData, targets, err, \
             flagver=1, CalJy=[38.5], BLInt=30., AtmInt=20.,tau0=0.1, 
             PointTab=None, prior=None, priorModel = True, PSF=None):
    """ Initial calibration of Mustang (PAR) data

    Any prior calibration tables are removed and then does the following:
      1) Generate an initial Cal table with a time increment of 1/4 of the
         lesser of BLInt and AtmInt
      2) Gain and Weight calibration, conversion to Jy uses CalJy
         CalJy can contain either a single value for all detectors or
         one value per detector.  This is the value of the cal in Jy.
      3) a "Baseline" per detector offsets on timescales longer than BLInt
         are determined and applied
      4) "Atmosphere" calibration determining one offset per detector per scan
         and a common mode offset on time scales longer than AtmInt.
         Opacity corrections are based on tau0 (zenith opacity in nepers)
      5) If PointTab is specified, pointing corrections are applied.
    When the procedure is finished, data cal be calibrated using the highest
    Cal table.
    If prior is given it is a Clean image if the target to be subtracted
    prior to the Baseline and Atmosphere calibration.
    Note: this only makes sense when all targets are covered by prior.
    If this option is used the instrumental PSF must also be provided in PSF.
    
    inData  = OTF data set to be calibrated
    targets = list of target names, empty list = all
    err     = Python Obit Error/message stack
    flagver =
    CalJy   = Array of the cal in Jy.  CalJy can contain either a single value
              for all detectors or one value per detector.  
    BLInt   = Baseline filter shortest timescale in sec
    AtmInt  = Atmospheric filter shortest timescale in sec
    tau0    = zenith opacity in nepers
              this can be either a scalar constant opacity, or a table in the form
              of a list of lists of time (days) and opacity, e.g.:
              # From Ron's weather server from CLEO
              tau0 = [[0.000, 0.102],        \
                      [0.042, 0.100],  \
                      [0.083, 0.100], \
                      [0.250, 0.154]] 
    PointTab= a table of pointing offsets in time order
            [time(day) d Xel (asec), d el (asec)]
            an example:
            PointTab=[[0.06189,-0.38, 0.38],  \
                      [0.09888,-2.77,-1.80],  \
                      [0.13052,-2.10,-0.92],  \
                      [0.16667,-1.86, 1.07],  \
                      [0.19509,-2.32, 1.29]]
           Such a table can be generated from a dataset by CalImage
    prior   = If given, a CLEAN image covering all targets given,
            This model will be subtracted from the data prior to
            "Baseline" and "Atmosphere" calibration
            If this option is used the instrumental PSF must also be provided
            in PSF.
    priorModel  = If prior defined then priorModel is true if the model (CC table)
            is to be used, if False then the image itself is used.
    PSF     = If prior is given, this is the instrumental PSF to use in the
            subtraction.
    """

    # Initial calibration
    # delete any prior calibration tables
    print "Remove previous calibration"
    OTF.ClearCal(inData,err)
    OErr.printErrMsg(err, "Error deleting prior cal tables")

    solInt = min (BLInt, AtmInt)


    # Create an initial dummy table with a interval 1/4 of the shortest
    # Filter type solution interval.
    inter = solInt/4
    print "Create initial calibration table, interval",inter
    OTFGetSoln.POTFGetDummyCal (inData, inData, inter, 1, 1, err)
    OErr.printErrMsg(err, "Error creating initial cal table")

    # Gain/Weight calibration
    print "Gain/Weight  calibration"
    inInfo = inData.List
    inInfo.set("calJy", CalJy)
    inInfo.set("doWate", True) # Do weight calibration
    OTFGetSoln.POTFGetSolnPARGain(inData, inData, err)
    OErr.printErrMsg(err, "Error with Gain/Weight calibration")

    # Update OTF Cal table
    OTF.Soln2CalInput["InData"]  = inData       # Input data object
    OTF.Soln2CalInput["oldCal"]  = 1            # Use initial cal
    OTF.Soln2CalInput["newCal"]  = 2            # New cal table
    OTF.Soln2Cal(err, OTF.Soln2CalInput)        # Apply
    OErr.printErrMsg(err, "Error updating Cal table with Soln")

    ################################## Target specification #################################
    inInfo.set("Targets", targets)       # select only target data
    inInfo.set("Stokes", "    ")         # Set Stokes
    inInfo.set("doCalSelect", True)       
    inInfo.set("flagVer", flagver)
    inInfo.set("gainUse", 0)
    
    ########################### Residual data from prior model #############################
    # Get prior model, if any and compute residual OTF
    if prior!=None:
        print "Using prior model, priorModel=",priorModel
        gainuse = 0
        inInfo.set("gainUse", gainuse)
        inInfo.set("doCalib", 1)
        inInfo.set("flagVer", flagver)
        if priorModel:
            # Use CC table
            resid = OTFUtil.PSubModel(inData, None, prior, PSF, err)
        else:
            # Use Image
            prior.Open(Image.READONLY,err)
            prior.Read(err)
            prior.Close(err)
            resid = OTF.PScratch (inData, err)
            OTFUtil.PSubImage(inData, resid, prior.FArray, prior.Desc, err)
        OErr.printErrMsg(err, "Error with residial data")
    else:
        print "No prior model used"
        resid = inData

    ############################## "Baseline" filter  #############################
    print "Baseline Filter"
    solint = BLInt/86400.0
    inInfo.set("solInt", solint)
    inInfo.set("doCalSelect", True)
    inInfo.set("flagVer", flagver)
    gainuse = 2
    inInfo.set("gainUse", gainuse)
    inInfo.set("doCalib", 1)
    OTFGetSoln.PFilter(resid, inData, err)
    OErr.printErrMsg(err, "Error with Baseline calibration")

    # Soln2Cal parameters for filter cal (most defaulted)
    OTF.Soln2CalInput["InData"]  = inData      # Input data object
    OTF.Soln2CalInput["oldCal"]  = 2           # Use gain cal output
    OTF.Soln2CalInput["newCal"]  = 3           # New cal table
    OTF.Soln2Cal(err, OTF.Soln2CalInput)  # Apply
    OErr.printErrMsg(err, "Error updating Cal table with Soln")
    
    ########################### Residual data from prior model #############################
    # Get prior model, if any and compute residual OTF
    if prior!=None:
        print "Using prior model, priorModel=",priorModel
        gainuse = 0
        inInfo.set("gainUse", gainuse)
        inInfo.set("doCalib", 1)
        inInfo.set("flagVer", flagver)
        if priorModel:
            # Use CC table
            resid = OTFUtil.PSubModel(inData, None, prior, PSF, err)
        else:
            # Use Image
            prior.Open(Image.READONLY,err)
            prior.Read(err)
            prior.Close(err)
            resid = OTF.PScratch (inData, err)
            OTFUtil.PSubImage(inData, resid, prior.FArray, prior.Desc, err)
        OErr.printErrMsg(err, "Error with residial data")
    else:
        print "No prior model used"
        resid = inData

    ############################## Common atmosphere + offset #############################
    print "Common atmosphere removal"
    inInfo.set("Tau0", tau0)                   # Opacity table
    inInfo.set("doCalSelect", True)
    inInfo.set("flagVer", flagver)
    gainuse = 0
    inInfo.set("gainUse", gainuse)
    inInfo.set("doCalib", 1)
    solint = AtmInt/86400.0
    inInfo.set("solInt", solint)
    clipsig = 5.0
    inInfo.set("ClipSig", clipsig)
    plotDet = -10
    inInfo.set("plotDet", plotDet)
    OTFGetSoln.PMBBase(resid, inData, err)
    OErr.printErrMsg(err, "Error with atmosphere calibration")
    
    # Soln2Cal for Atm cal
    OTF.Soln2CalInput["InData"]  = inData      # Input data object
    OTF.Soln2CalInput["oldCal"]  = 3           # Use baseline cal output
    OTF.Soln2CalInput["newCal"]  = 4           # New cal table
    OTF.Soln2Cal(err, OTF.Soln2CalInput)       # Apply
    OErr.printErrMsg(err, "Error updating Cal table with Soln")

    ############################### Pointing correction ##################################
    if PointTab != None:
        print "Apply pointing corrections"
        inInfo.set("POffset", PointTab)
        OTFGetSoln.POTFGetSolnPointTab(inData, inData, err)
        OErr.printErrMsg(err, "Error with pointing corrections")
        
        # Soln2Cal for Point cal
        OTF.Soln2CalInput["InData"]  = inData       # Input data object
        OTF.Soln2CalInput["oldCal"]  = 4           # Use baseline cal output
        OTF.Soln2CalInput["newCal"]  = 5           # New cal table
        OTF.Soln2Cal(err, OTF.Soln2CalInput)       # Apply
        OErr.printErrMsg(err, "Error updating Cal table with Soln")
    

    # End InitCal


import OPlot, OTFRec, FArray
def PlotData (inData, targets, scans, feeds, err, \
              output="None", bgcolor=0, nx=1, ny=1):
    """ Plot selected data

    Plot data in inData selected by targets and scans
    inData  = OTF data set to be plotted, any calibration and editing will be applied
    targets = list of target names, empty list = all
    scans   = Range of scan number, 0's => all
    feeds   = list of feeds to plot
    err     = Python Obit Error/message stack
    output  = name and type of output device:
              "None"  interactive prompt
              "xwin"  X-Window (Xlib)
              "gcw"   Gnome Canvas Widget (interacts with ObitTalk)
              "ps"    PostScript File (monochrome)
              "psc"   PostScript File (color)
              "xfig"  Fig file
              "png"   PNG file
              "jpeg"  JPEG file
              "gif"   GIF file
              "null"  Null device
    bgcolor   = background color index (1-15), symbolic names:
                BLACK, RED(default), YELLOW, GREEN, 
                AQUAMARINE, PINK, WHEAT, GRAY, BROWN,
                BLUE, BLUEVIOLET, CYAN, TURQUOISE,
                MAGENTA, SALMON, WHITE
    nx        = Number of horizontal subpages
    ny        = Number of vertical subpages
    """
    # Set selection
    inInfo = inData.List
    inInfo.set("doCalSelect",True)       # Do selection
    if len(targets)>0:
        inInfo.set("Targets", targets)   # select only target data
    if scans[0]>0 and scans[1]>0:
        lscans = scans
    else:
        lscans = [1,10000000]
    inInfo.set("Scans",lscans)
    if len(feeds)>0 and feeds[0]>0:
        lfeeds = feeds
    else:
        lfeeds = []
        naxis = inData.Desc.Dict["inaxes"]
        nfeed = naxis[1]*max(naxis[2],1)*max(naxis[3],1)*max(naxis[4],1)
        for i in range(0,nfeed):
            lfeeds.append(i+1)
  
    inData.Open(OTF.READONLY,err)
    nrec = inData.Desc.Dict["nrecord"]
    plotdata = [[]]       # time
    datagood = [False]    # data validity
    for j in lfeeds:
        plotdata.append([])
        datagood.append(False)
    # Loop over data
    for i in range(1,nrec+1):
        rec=inData.ReadRec (err)
        # eod of file?
        eof = 'EOF' in rec.__dict__
        if eof:
            break
        OErr.printErrMsg(err, "Error reading OTF data")
        plotdata[0].append(rec.time*24.0)  # Save time
        k = 1
        for j in  lfeeds:
            if rec.data[j-1][1]!=0.0:
                plotdata[k].append(rec.data[j-1][0])
                datagood[k] = rec.data[j-1][0]!=0.0
            else:
                plotdata[k].append(FArray.fblank)
            k += 1
    
    inData.Close(err)

    # Anything selected
    if len(plotdata[0])<=0:
        raise RuntimeError,'No data selected to plot'
    # plot
    plot = OPlot.newOPlot("plot", err, nx=nx, ny=ny, output=output, bgcolor=bgcolor)
    plotInfo = plot.List
    plotInfo.set("XLABEL","Time(hr)")
    plotInfo.set("YLABEL","Flux density(Jy)")
    nplot = len(plotdata)-1
    for i in range(0,nplot):
        plotInfo.set("TITLE","Detector "+str(lfeeds[i]))
        if datagood[i+1] and max(plotdata[i+1])>min(plotdata[i+1]):
            OPlot.PXYPlot(plot, 2, plotdata[0], plotdata[i+1], err)
        # Tolerate error and keep going
        if err.isErr:
            err.Clear()
            print "ERROR occured in plotting detector",i+1
    OPlot.PShow(plot,err)
    OErr.printErrMsg(err, "Error plotting OTF data")
   # end PlotData

import TimeFilter
def PlotPower (inData, targets, scans, feeds, err, \
              output="None", bgcolor=0, nx=1, ny=1):
    """ Plot power spectrum of selected data

    Plot data in inData selected by targets and scans
    inData  = OTF data set to be plotted, any calibration and editing will be applied
    targets = list of target names, empty list = all
    scans   = Range of scan number, 0's => all
    feeds   = list of feeds to plot
    err     = Python Obit Error/message stack
    output  = name and type of output device:
              "None"  interactive prompt
              "xwin"  X-Window (Xlib)
              "gcw"   Gnome Canvas Widget (interacts with ObitTalk)
              "ps"    PostScript File (monochrome)
              "psc"   PostScript File (color)
              "xfig"  Fig file
              "png"   PNG file
              "jpeg"  JPEG file
              "gif"   GIF file
              "null"  Null device
    bgcolor   = background color index (1-15), symbolic names:
                BLACK, RED(default), YELLOW, GREEN, 
                AQUAMARINE, PINK, WHEAT, GRAY, BROWN,
                BLUE, BLUEVIOLET, CYAN, TURQUOISE,
                MAGENTA, SALMON, WHITE
    nx        = Number of horizontal subpages
    ny        = Number of vertical subpages
    """
    # Set selection
    inInfo = inData.List
    inInfo.set("doCalSelect",True)       # Do selection
    if len(targets)>0:
        inInfo.set("Targets", targets)   # select only target data
    if scans[0]>0 and scans[1]>0:
        lscans = scans
    else:
        lscans = [1,10000000]
    inInfo.set("Scans",lscans)
    if len(feeds)>0 and feeds[0]>0:
        lfeeds = feeds
    else:
        lfeeds = []
        naxis = inData.Desc.Dict["inaxes"]
        nfeed = naxis[1]*max(naxis[2],1)*max(naxis[3],1)*max(naxis[4],1)
        for i in range(0,nfeed):
            lfeeds.append(i+1)
  
    inData.Open(OTF.READONLY,err)
    nrec = inData.Desc.Dict["nrecord"]
    plotdata = [[]]       # time
    datagood = [False]    # data validity
    for j in lfeeds:
        plotdata.append([])
        datagood.append(False)
    # Loop over data
    for i in range(1,nrec+1):
        rec=inData.ReadRec (err)
        # eod of file?
        eof = 'EOF' in rec.__dict__
        if eof:
            break
        OErr.printErrMsg(err, "Error reading OTF data")
        plotdata[0].append(rec.time*24.0)  # Save time
        k = 1
        for j in  lfeeds:
            if rec.data[j-1][1]!=0.0:
                plotdata[k].append(rec.data[j-1][0])
                datagood[k] = rec.data[j-1][0]!=0.0
            else:
                plotdata[k].append(FArray.fblank)
            k += 1
    
    inData.Close(err)

    # Create/populate TimeSeries
    ts = TimeFilter.newTimeFilter ("timeSeries", len(plotdata[0]), len(lfeeds))
    k = 1
    nTime = len(plotdata[0])
    dTime = (plotdata[0][nTime-1] - plotdata[0][0]) / nTime
    for j in lfeeds:
        TimeFilter.PGridTime(ts, k-1, dTime, nTime, plotdata[0], plotdata[k])
        k += 1
    # Determine spectra
    TimeFilter.P2Freq (ts)

    # Anything selected
    if len(plotdata[0])<=0:
        raise RuntimeError,'No data selected to plot'
    # plot
    plot = OPlot.newOPlot("plot", err, nx=nx, ny=ny, output=output, bgcolor=bgcolor)
    plotInfo = plot.List
    plotInfo.set("XLABEL","Frequency (Hz)")
    plotInfo.set("YLABEL","Spectral Power")
    plotInfo.set("YOPT", "BCTSL")
    nplot = len(plotdata)-1
    for i in range(0,nplot):
        # Get power spectra
        powerDict = TimeFilter.PGetPower(ts, i)
        # Zero first term
        powerDict["data"][0] = powerDict["data"][1]
        plotInfo.set("TITLE","Detector "+str(lfeeds[i]))
        if datagood[i+1] and max(plotdata[i+1])>min(plotdata[i+1]):
            OPlot.PXYPlot(plot, 2, powerDict["freq"], powerDict["data"], err)
        OErr.printErrMsg(err, "Error plotting OTF data")
    OPlot.PShow(plot,err)
    OErr.printErrMsg(err, "Error plotting OTF data")
   # end PlotPower

import OPlot, OTFRec, FArray, OTFArrayGeom
def PlotElev (inData, targets, scans, feeds, err, \
              output="None", bgcolor=0, nx=1, ny=1):
    """ Plot selected data vs elevation

    Plot data in inData selected by targets and scans
    inData  = OTF data set to be plotted, any calibration and editing will be applied
    targets = list of target names, empty list = all
    scans   = Range of scan number, 0's => all
    feeds   = list of feeds to plot
    err     = Python Obit Error/message stack
    output  = name and type of output device:
              "None"  interactive prompt
              "xwin"  X-Window (Xlib)
              "gcw"   Gnome Canvas Widget (interacts with ObitTalk)
              "ps"    PostScript File (monochrome)
              "psc"   PostScript File (color)
              "xfig"  Fig file
              "png"   PNG file
              "jpeg"  JPEG file
              "gif"   GIF file
              "null"  Null device
    bgcolor   = background color index (1-15), symbolic names:
                BLACK, RED(default), YELLOW, GREEN, 
                AQUAMARINE, PINK, WHEAT, GRAY, BROWN,
                BLUE, BLUEVIOLET, CYAN, TURQUOISE,
                MAGENTA, SALMON, WHITE
    nx        = Number of horizontal subpages
    ny        = Number of vertical subpages
    """
    # Set selection
    inInfo = inData.List
    inInfo.set("doCalSelect",True)       # Do selection
    if len(targets)>0:
        inInfo.set("Targets", targets)   # select only target data
    if scans[0]>0 and scans[1]>0:
        lscans = scans
    else:
        lscans = [1,10000000]
    inInfo.set("Scans",lscans)
    if len(feeds)>0 and feeds[0]>0:
        lfeeds = feeds
    else:
        lfeeds = []
        naxis = inData.Desc.Dict["inaxes"]
        nfeed = naxis[1]*max(naxis[2],1)*max(naxis[3],1)*max(naxis[4],1)
        for i in range(0,nfeed):
            lfeeds.append(i+1)
  
    inData.Open(OTF.READONLY,err)
    ag = inData.ArrayGeom    # Get array geometry
    nrec = inData.Desc.Dict["nrecord"]
    plotdata = [[]]       # time
    datagood = [False]    # data validity
    for j in lfeeds:
        plotdata.append([])
        datagood.append(False)
    # Loop over data
    for i in range(1,nrec+1):
        rec=inData.ReadRec (err)
        # eod of file?
        eof = 'EOF' in rec.__dict__
        if eof:
            break
        OErr.printErrMsg(err, "Error reading OTF data")
        elev = ag.Elev(rec.time, rec.ra, rec.dec)
        plotdata[0].append(elev)  # Save elev
        k = 1
        for j in  lfeeds:
            if rec.data[j-1][1]!=0.0:
                plotdata[k].append(rec.data[j-1][0])
                datagood[k] = rec.data[j-1][0]!=0.0
            else:
                plotdata[k].append(FArray.fblank)
            k += 1
    
    inData.Close(err)

    # Anything selected
    if len(plotdata[0])<=0:
        raise RuntimeError,'No data selected to plot'
    # plot
    plot = OPlot.newOPlot("plot", err, nx=nx, ny=ny, output=output, bgcolor=bgcolor)
    plotInfo = plot.List
    plotInfo.set("XLABEL","Elev(deg)")
    plotInfo.set("YLABEL","Flux density(Jy)")
    nplot = len(plotdata)-1
    for i in range(0,nplot):
        plotInfo.set("TITLE","Detector "+str(lfeeds[i]))
        if datagood[i+1] and max(plotdata[i+1])>min(plotdata[i+1]):
            OPlot.PXYPlot(plot, 2, plotdata[0], plotdata[i+1], err)
        # Tolerate error and keep going
        if err.isErr:
            err.Clear()
            print "ERROR occured in plotting detector",i+1
    OPlot.PShow(plot,err)
    OErr.printErrMsg(err, "Error plotting OTF data")
   # end PlotElev

# Define CleanSkyModel input dictionary
CleanSkyModelInput={'structure':['CleanSkyModel',[('InData','Input OTF'),
                                           ('DirtyName','Dirty image name, None = scratch'),
                                           ('CleanName','Clean image name, None = scratch'),
                                           ('outDisk','FITS file disk number for output'),
                                           ('PSF','Input Beam image for telescope'),
                                           ('scan','scan range to be processed'),
                                           ('target','list of targets to include'),
                                           ('nx','number of pixels in x = RA'),
                                           ('ny','number of pixels in y = de'),
                                           ('xCells','Cell spacing in x (asec)'),
                                           ('yCells','Cell spacing in y (asec)'),
                                           ('minWt','minimum summed weight in gridded image wrt max '),
                                           ('ConvType','Conv. fn type, 0=pillbox,3=Gaussian,4=exp*sinc,5=Sph wave'),
                                           ('ConvParm','Conv. fn parameters'),
                                           ('gainUse','cal. table version, -1=none'),
                                           ('flagVer','flag table version, -1=none'),
                                           ('Niter','Maximum number of CLEAN iterations def[100]'),
                                           ('Patch','Beam patch in pixels [def 100]'),
                                           ('noResid','If True do not include residuals in restored image'),
                                           ('BeamSize','Restoring beam FWHM (deg)'),
                                           ('Window','list of Clean windows, def [[-1,20,50,50]]'),
                                           ('Gain','CLEAN loop gain def [0.1]'),
                                           ('autoWindow','Automatically set Windows? [False]')]],
                    'InData':None, 'outDisk':0, 'PSF':None, "scan":[1,1000000], 'target':["All"],
                    'DirtyName':None, 'CleanName':None, 'Window':[[-1,20,50,50]],
                    'nx':100, 'ny':100, 'xCells':2.0, 'yCells':2.0, 'minWt':0.0001,
                    'ConvParm':[0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0], 'ConvType':5,
                    'gainUse':-1, 'flagVer':-1, 'noResid':False,
                    'Niter':100, 'Patch':100, 'BeamSize':0.0, 'Gain':0.1, 'autoWindow':False}
def CleanSkyModel (err, input=CleanSkyModelInput):
    """ Create a CLEAN image sky model of an OTF data set

    Does imaging and cleaning of specified data.
    Returns an image model of the CLEAN components convolved with PSF
    err     = Python Obit Error/message stack
    input   = input parameter dictionary, for interactive use, the function input
              will display the contents in human readable format.

    Input dictionary entries:
    InData  = Python input OTF to calibrate
    DirtyName = Dirty image name, None = scratch
    CleanName = Clean image name, None = scratch
    outDisk = Disk number for output files
    PSF     = Image with telescope psf
    scans   = scan range to image
    target  = list of targets to include, center will be position of first
    nx      = number of pixels in 'X' = RA
    ny      = number of pixels in 'Y' = dec
    xCells  = Cell spacing in x (asec)
    yCells  = Cell spacing in y (asec)
    ConvType= Convolving function Type 0=pillbox,3=Gaussian,4=exp*sinc,5=Sph wave
    ConvParm= Convolving function parameters depends on ConvType
      Type 2 = Sinc, (poor function - don't use)
        Parm[0] = halfwidth in cells,
        Parm[1] = Expansion factor
      Type 3 = Gaussian,
        Parm[0] = halfwidth in cells,[def 3.0]
        Parm[1] = Gaussian with as fraction or raw beam [def 1.0]
      Type 4 = Exp*Sinc
        Parm[0] = halfwidth in cells, [def 2.0]
        Parm[1] = 1/sinc factor (cells) [def 1.55]
        Parm[2] = 1/exp factor (cells) [def 2.52]
        Parm[3] = exp power [def 2.0]
      Type 5 = Spherodial wave 
        Parm[0] = halfwidth in cells [def 3.0]
        Parm[1] = Alpha [def 5.0]
        Parm[2] = Expansion factor [not used]
    gainUse = version number of prior table (Soln or Cal) to apply, -1 is none
    flagVer = version number of flagging table to apply, -1 is none
    Niter       = Maximum number of CLEAN iterations
    Patch       = Beam patch in pixels [def 100]
    noResid     = If True do not include residuals in restored image
    BeamSize    = Restoring beam (deg)
    Gain        = CLEAN loop gain
    Window      = list of Clean windows
    autoWindow  = True if autoWindow feature wanted.
    """
    ################################################################
    # Get input parameters
    dim = [1,1,1,1,1]
    inData  = input["InData"]
    inInfo = inData.List
    # Checks
    if not OTF.PIsA(inData):
        raise TypeError,'Image: Bad input OTF'
    if err.isErr: # existing error?
        return None

    # Set calibration 
    gainUse = input["gainUse"]
    flagVer = input["flagVer"]
    # Default table versions (the Obit routines will do this as well)
    if gainUse == 0:   # Get highest numbered OTFCal table
        gainUse = Obit.OTFGetHighVer(inData.me, "OTFCal")
    if gainUse == 0:   # Doesn't seem to be one, try OTFSoln
        gainUse = Obit.OTFGetHighVer(inData.me, "OTFSoln")
    if gainUse == 0:   # Must not be one
        gainUse = -1
    inInfo.set("doCalSelect", True)       
    inInfo.set("flagVer", flagVer)
    if gainUse>0:
        inInfo.set("doCalib", 1)       
        inInfo.set("gainUse", gainUse)       
        
    # Get position from OTF
    target = input["target"]
    pos =  GBTUtil.GetTargetPos(inData, target[0], err)
    ra  = pos[0]                      # ra of center
    dec = pos[1]                      # dec of center

    # Set imaging/calibration parameters
    outDisk = input["outDisk"]
    # Imaging parameters
    DirtyName = input["DirtyName"]
    if  DirtyName==None:
        DirtyName = "tmp"+target[0]+"Dirty.fits"
    OTF.ImageInput["InData"]  = inData
    OTF.ImageInput["OutName"] = DirtyName
    OTF.ImageInput["disk"]    = outDisk 
    OTF.ImageInput["Beam"]    = input["PSF"]
    OTF.ImageInput["ra"]      = ra
    OTF.ImageInput["dec"]     = dec
    OTF.ImageInput["xCells"]  = input["xCells"]
    OTF.ImageInput["yCells"]  = input["yCells"]
    OTF.ImageInput["nx"]      = input["nx"]
    OTF.ImageInput["ny"]      = input["ny"]
    OTF.ImageInput["gainUse"] = 0
    OTF.ImageInput["flagVer"] = flagVer
    OTF.ImageInput["minWt"]   = input["minWt"]
    OTF.ImageInput["ConvType"]  = input["ConvType"]
    OTF.ImageInput["ConvParm"]  = input["ConvParm"]
    
    #  Make image
    DirtyImg = OTF.makeImage(err, OTF.ImageInput)
    OErr.printErrMsg(err, "Error making initial image")

    # Clean to get sky model
    cleanFile = input["CleanName"]
    if cleanFile:
        CleanImg = Image.newPFImage("Clean Image", cleanFile, outDisk, False, err)
    else:  # scratch
        CleanImg = DirtyImg.Scratch(err)
    
    CleanObj = CleanOTF.PCreate("Clean", DirtyImg, input["PSF"], CleanImg, err)
    for win in input["Window"]:
        CleanOTF.PAddWindow(CleanObj, win, err)

    # CLEAN parameters
    CleanOTF.CleanInput["CleanOTF"] = CleanObj     # Clean object
    CleanOTF.CleanInput["autoWindow"]= input["autoWindow"]
    CleanOTF.CleanInput["Patch"]    = input["Patch"]
    CleanOTF.CleanInput["Niter"]    = input["Niter"]
    CleanOTF.CleanInput["Gain"]     = input["Gain"]
    CleanOTF.CleanInput["BeamSize"] = input["BeamSize"]
    CleanOTF.CleanInput["noResid"]  = input["noResid"]
    #DEBUG CleanOTF.CleanInput["disp"]     = input["disp"]
    CleanOTF.CleanInput["minFlux"]  = 0.0
    
    resid = CleanObj.Clean                 # Copy image just produced
    Image.PCopy(DirtyImg, resid, err)      # to clean(resid)        
    #input(CleanInp)
    CleanOTF.PClean (err, CleanOTF.CleanInput)  # Do Clean
    OErr.printErrMsg(err, "Error Cleaning")

    # Replace with Clean model
    CCTab    = CleanImg.NewTable(1,"AIPS CC", 0, err)
    OErr.printErrMsg(err, "Error getting CCTable")
    # Use all components
    dim = [1,1,1,1,1]
    CCTab.List.set("BComp", 1)
    CCTab.List.set("EComp", 0)

    # Make model image, CC convolved with beam
    CleanImg.Open(Image.READONLY, err)
    model = CleanImg.ReadPlane(err)
    CleanImg.Close(err)
    model = OTFUtil.PConvBeam (CCTab, CleanObj.Beam, model, err);
    OErr.printErrMsg(err, "Error making model")
    # Write model to clean image
    CleanImg.Open(Image.READWRITE, err)
    CleanImg.WriteFA(model, err)
    CleanImg.Close(err)

    # Delete Dirty Image?
    if input["DirtyName"]==None:
        DirtyImg.Zap(err)

    return CleanImg
    # end CleanSkyModel

