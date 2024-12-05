# Basic mosaicing functions
import Image, ImageDesc, MosaicUtil, OSystem, OErr, History
from ObitTask import ObitTask
from OTObit import imlod, setname
exec(open("CheckPos.py").read())

def StitchF (nameList, stktrans, SumWtImage, SumWt2, err, beam=None, 
             galactic=False, disk=1, restart=0, maxd=0.6):
    """ Accumulate Mosaic looping over planes of a list of cubes 

        nameList  = list of pointing names
        stktrans  = Stokes or transition, e.g. "I"
        SumWtImage= Image*wt accumulator
        SumWt2    = Wt**2 accumulator
        err       = Obit err stack
        beam      = Beam to convolve to if given e.g. [8.0, 8.0, 0.0]
        galactic  = Convert the output to Galactic coordinates?
        disk      = AIPS disk for working files
        restart   = restart channel no. 0-rel
        maxd      = Maximum distance (deg) from pointing center to consider overlap.
    """
    # How many channels? */
    nchan = SumWtImage.Desc.Dict["inaxes"][SumWtImage.Desc.Dict["jlocf"]]
    # Crude position setup
    #   Mosaic center
    mosPos   = SumWtImage.Desc.Dict['crval'][0:2]
    # How close is close enough for proper test?
    mosDelta = abs(SumWtImage.Desc.Dict['cdelt'][0])*SumWtImage.Desc.Dict['inaxes'][0]
    # Set up for beam convolution, input=FITS, output=AIPS
    if beam:
        cvl = ObitTask('Convol'); cvl.DataType='FITS'; cvl.inDisk=0
        cvl.outDType ='AIPS';cvl.outDisk=disk; cvl.outSeq=2; cvl.outClass='Convol'
    # Loop over images
    i = -1;
    acc = None; acct = None
    for name in nameList:
        print (name)
        acc = None; acct = None
        i += 1
        # Check if file exists
        f = fitsdir+name+'_'+stktrans+".fits"
        if not os.path.exists(f):
            # try gzipped
            f = fitsdir+name+'_'+stktrans+".fits.gz"
            if not os.path.exists(f):
                print (f,"does not exist")
                continue
        # Does this overlap the mosaic?
        if CheckPos(name, SumWtImage, cat, galactic, prtLv, err, maxd=maxd):
            acct  = Image.newPFImage("accum", f, 0, True, err)
            try:
                # Get overlap
                if galactic:
                    blc=[1,1,1,1,1,1,1]; trc = acct.Desc.Dict['inaxes']               # Galactic
                else:
                    (blc,trc) = MosaicUtil.PGetOverlap(acct, SumWtImage, err) # Equatorial 
                if (trc[0]>blc[0]) and ((trc[1]>blc[1])):
                    # Images will be rotated due to Equatorial->Galactic conversion,
                    # if there is any overlap, take all of image, workaround for IO problem
                    # Channel selection 
                    blc[0]=1; blc[1]=1;
                    trc[0:2]=acct.Desc.Dict['inaxes'][0:2]
                    blc[2] = minCh; trc[2] = maxCh  # Limit channels
                    # Convolve for common beam?
                    if beam:
                        cvl.inFile = f; cvl.outName = name[0:12]
                        cvl.Beam = beam; cvl.g
                        acc = Image.newPAImage("in",cvl.outName, cvl.outClass, cvl.outDisk, 2, True, err)
                    else:
                        # straight copy FITS to AIPS
                        acc = imlod(f,0,name[0:12],stktrans+strId,disk,1,err)
                    # pointing correction? - Apply offset to reference pixel
                    d = acc.Desc.Dict; 
                    ox = (cat[name][3]/3600.)/d['cdelt'][0]; oy = (cat[name][4]/3600.)/d['cdelt'][1]; 
                    d['crpix'][0] += ox; d['crpix'][1] += oy # "+" since offsetting reference
                    acc.Desc.Dict = d;
                    acc.UpdateDesc(err)
                    # Use first plane to normalize RMS factors
                    Image.PGetPlane (acc, None, [1,1,1,1,1], err)
                    OErr.printErrMsg(err, "Error reading image for "+Image.PGetName(acc))
                    factor = 1.0;
                    RMS = acc.FArray.RMS;    # combined plane RMS
                    maxRMS = 20*RMS          # maximum acceptable plane RMS
                    factor *= cat[name][2]   # any additional factors
                    # DEBUG
                    # DEBUGacc.List.set("BLC",[1,1,1,1,1,1,1]); acc.List.set("BLC",[0,0,0,1,1,1,1]); 
                    #print("StitchF: acc.List before",acc.List.Dict)
                    if prtLv>=1:
                        print ("Accumulate",name,"factor",factor, "maxRMS",maxRMS)
                        print ("overlap", blc, trc)
                    # Accumulate
                    MosaicUtil.PWeightImage(acc, factor, SumWtImage, SumWt2, err, minGain=0.30,
                                            iblc=blc, itrc=trc, restart=restart, planeWt=True,
                                            hwidth=2, doGPU=False, inWtImage=None,
                                            maxRMS=maxRMS, antSize=antSize, prtLv=prtLv)
            except Exception as exception:
                OErr.printErr(err)
                print (exception)
                print ("Failed")
                print ("StitchF: blc",blc,"trc",trc)
                print("acc.List",acc.List.Dict)
                err.Clear()
            if acc:
                acc.Zap(err); del acc; acc = None
            if acct:
                del acct; acct = None
        else:
            print ("Not in mosaic:",name)
            if acct:
                del acct; acct = None
    
    # Cleanup
    if beam:
        del cvl
 # end StitchF

# Loop over targets
for target in targets:
    # Make accumulation images, cells, size[0] x size[1]
    if not restart:
        # make AIPS copy of template file 11 x 11 cells
        f = fitsdir+"tmpl.fits"
        if not os.path.exists(f):
            f = fitsdir+"tmpl.fits.gz"
            if not os.path.exists(f):
                print (f,"does not exist")
                continue
        targ  = Image.newPFImage("target", f, 0, True, err)
        d     = targ.Desc.Dict
        si     = ObitTask("SubImage")
        setname(targ,si)
        # select central 11 x 11 cells
        si.BLC = [int(d['crpix'][0]-5), int(d['crpix'][1]-5), 1]
        si.TRC = [int(d['crpix'][0]+5), int(d['crpix'][1]+5),  d['inaxes'][2]]
        si.TRC[2] = min(si.TRC[2],maxCh)  # Limit channels
        si.BLC[2] = max(si.BLC[2],minCh)
        si.outDisk = accumdisk;  si.outSeq=accumseq; si.outDType="AIPS"
        si.outName =  target[0]; si.outClass = stktrans+"Tp"; 
        si.g
        beam = target[2]  # Beam size (asec, asec, deg)
        tmpl = Image.newPAImage("tmpl", si.outName, si.outClass, si.outDisk, si.outSeq, True, err)
        OErr.printErr(err)     
        # Fiddle header to match mosaic
        d    = tmpl.Desc.Dict
        d["cdelt"][0]=-cells/3600.;   d["cdelt"][1]=cells/3600.      # Cellsize
        d["crval"][0] = target[1][0]; d["crval"][1] = target[1][1];  # Center
        if galactic:
            d["ctype"][0] = 'GLON-SIN';   d["ctype"][1] = 'GLAT-SIN';
        if beam:
            d['beamMaj'] = beam[0]/3600.; d['beamMin'] = beam[1]/3600.;d['beamPA'] = beam[2]
        if stktrans=='Q':
            d["crval"][3] = 2.
        if stktrans=='U':
            d["crval"][3] = 3.
        tmpl.Desc.Dict=d;tmpl.UpdateDesc(err)
        f    = target[0]
        # Accumulators in AIPS
        SumWtImage = Image.newPAImage("WtI", f, stktrans+"WI", accumdisk, accumseq, False, err)
        SumWt2     = Image.newPAImage("Wt2", f, stktrans+"W2", accumdisk, accumseq, False, err)
        print ("Create accumulators")
        MosaicUtil.PMakeMaster(tmpl, size, SumWtImage, SumWt2, err)
        # Put reference pixel at center
        d = SumWtImage.Desc.Dict
        d['crpix'][0] = d['inaxes'][0]/2; d['crpix'][1] = d['inaxes'][1]/2;
        SumWtImage.Desc.Dict = d; SumWtImage.UpdateDesc(err)
        d = SumWt2.Desc.Dict
        d['crpix'][0] = d['inaxes'][0]/2; d['crpix'][1] = d['inaxes'][1]/2;
        SumWt2.Desc.Dict = d; SumWt2.UpdateDesc(err)
    else: # restarting
        f    = target
        SumWtImage = Image.newPAImage("WtI", f, stktrans+"WI", accumdisk, accumseq, True, err)
        SumWt2     = Image.newPAImage("Wt2", f, stktrans+"W2", accumdisk, accumseq, True, err)
        restart = False;ch0 = 0;  # Not again
    # end startup
    
    # accumulate pointings into mosaic accumulators
    OErr.printErr(err)
    OErr.PLog(err,OErr.Info, "Start interpolation"); OErr.printErr(err)
    StitchF (pointings, stktrans, SumWtImage, SumWt2, err, beam=beam, disk=accumdisk,
             galactic=galactic, restart=ch0,maxd=maxd)
    OErr.PLog(err,OErr.Info, "Finish interpolation"); OErr.printErr(err)
    
    # Normalize to FITS
    f = fitsdir+target[0]+stktrans+'_Mosaic.fits'
    print ("Normalize",f)
    mosaic = Image.newPFImage("Mosaic", f, 0, False, err)
    SumWtImage.Clone(mosaic, err)
    MosaicUtil.PNormalizeImage(SumWtImage, SumWt2, mosaic, err, minWt=minWt)
    
    # History
    outHistory  = History.History("history", mosaic.List, err)
    outHistory.Open(History.READWRITE, err)
    outHistory.TimeStamp(" Start Obit Mosaic",err)
    History.PWriteRec(outHistory,-1,"Mosaic   / antSize = "+str(antSize),err)
    History.PWriteRec(outHistory,-1,"Mosaic   / minWt = "+str(minWt),err)
    History.PWriteRec(outHistory,-1,"Mosaic   / minCh = "+str(minCh),err)
    History.PWriteRec(outHistory,-1,"Mosaic   / maxCh = "+str(maxCh),err)
    outHistory.Close(err)

    # Cleanup
    zap(SumWtImage); zap(SumWt2); zap(tmpl)
    del tmpl, SumWtImage, SumWt2, mosaic
# end loop

