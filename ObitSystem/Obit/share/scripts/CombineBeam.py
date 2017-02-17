# Script to combine Holography beams using multiple resolutions and frequency ranges
space = 0.0212 # cell spacing in deg (highest resolution, frequency)
ncell = 101    # number of cells in Az, El (to cover lowest resolution, low frequency)
nFcomb= 2      # Number of frequency bins to combine
stoke = 'LL'   # Stokes
stoke = 'RR'   # Stokes
weight = [0.25, 1.0]  # Relative weights of resolutions
# List of base of input file names, ordered by increasing resolution, freq
dataRoot = [ \
    ['C_Band/HoloCBandLoCoarse', 'C_Band/HoloCBandLoFine'],   \
    ['C_Band/HoloCBandHiCoarse', 'C_Band/HoloCBandHiFine'],   \
]
outRoot = 'C_Band/HoloCBand'  # Output name base

import Image, ImageDesc, Table, OData, FArray, FInterpolate, ImageUtil, OErr, History
fblank = FArray.fblank     # Magic blanking value
# get file names, Image objects
dataFiles = []; dataImag = []; i = 0
for freqRoots in dataRoot:
    dataFiles.append([]); dataImag.append([])
    for root in freqRoots:
        f = root+'.'+stoke+'All.fits'
        dataFiles[i].append(f)
        xf = Image.newPFImage(root,f,0,True,err)
        dataImag[i].append(xf)
        z=Image.PFullInstantiate(xf, Image.READONLY, err)
        OErr.printErrMsg(err, message='Error with input image')
    i += 1  # next freq bin

# end loops

# Clone output from 1st input
f = outRoot+'.'+stoke+'.fits'
outIm = Image.newPFImage(outRoot,f,0,False,err)
dataImag[0][0].Clone(outIm, err)
OErr.printErrMsg(err, message='Error creating output image')
# Redefine
d = outIm.Desc.Dict
d['inaxes'][0] = ncell; d['inaxes'][1] = ncell;
d['inaxes'][3] *= nFcomb; numIF = d['inaxes'][3]; nChan = d['inaxes'][2];
d['cdelt'][0] = space; d['cdelt'][1] = space; 
d['crpix'][0] = (ncell+1)/2.0; d['crpix'][1] = (ncell+1)/2.0; 
outIm.Desc.Dict = d
z=Image.PFullInstantiate(outIm, Image.WRITEONLY, err)
OErr.printErrMsg(err, message='Error creating output image')

# Generate output FQ Table (delete old)
z = outIm.ZapTable('AIPS FQ',1,err)
outFQ = outIm.NewTable(Table.READWRITE,'AIPS FQ',1,err,numIF=numIF)
outFQ.Open(Table.READWRITE, err)
outRefFreq = outIm.Desc.Dict['crval'][2]
# Generate row structure
outRow = {'FRQSEL': [1], 'NumFields': 7, 'Table name': 'AIPS FQ', '_status': [0], \
          'CH WIDTH': numIF*[0.0], 'TOTAL BANDWIDTH':  numIF*[0.0], 'SIDEBAND':numIF*[0], \
          'RXCODE': numIF*['        '], 'BAND': numIF*[0.0], 'IF FREQ': numIF*[0.0], \
          }
# Loop over frequency bins, get IF Info from first resolution
for i in range (0,nFcomb):
    x   = dataImag[i][0]
    d =  x.Desc.Dict
    inRefFreq = d['crval'][2]
    inNIF     = d['inaxes'][3]
    xfq = dataImag[i][0].NewTable(Table.READONLY,"AIPS FQ",1,err)
    xfq.Open(Table.READONLY,err)
    rrr = xfq.ReadRow(1,err)
    for j in range (0,inNIF):
        outRow['IF FREQ'][i*inNIF+j]  = rrr['IF FREQ'][j] + inRefFreq - outRefFreq
        outRow['SIDEBAND'][i*inNIF+j] = rrr['SIDEBAND'][j]
        outRow['CH WIDTH'][i*inNIF+j] = rrr['CH WIDTH'][j]
        outRow['TOTAL BANDWIDTH'][i*inNIF+j] = rrr['TOTAL BANDWIDTH'][j]
    xfq.Close(err); del xfq, rrr

outFQ.WriteRow(1, outRow, err)
outFQ.Close(err)
OErr.printErrMsg(err, message='Error creating output FQ Table')
# Cleanup
del outRow

# Work arrays
work  = FArray.FArray('work', [ncell,ncell])  # work array for interpolated image
work2 = FArray.FArray('work2',[ncell,ncell])  # work array for interpolated image
accWI = FArray.FArray('accWI',[ncell,ncell])  # Sum wt * image
accW  = FArray.FArray('accW',[ncell,ncell])   # Sum wt

# Pixel array images
XPixIm = Image.Image("XPixIm")
YPixIm = Image.Image("YPixIm")
outIm.List.set('BLC',[1,1,1,1,1,1,1])
outIm.List.set('TRC',[0,0,1,1,1,1,1])
outIm.Open(Image.READWRITE,err); outIm.Close(err)
Image.PCloneMem(outIm, XPixIm, err)
Image.PCloneMem(outIm, YPixIm, err)
outIm.List.set('TRC',[0,0,0,0,1,1,1])
outIm.Open(Image.READWRITE,err); outIm.Close(err)

# Interpolator 
fi = FInterpolate.FInterpolate("Interp",dataImag[0][0].FArray, dataImag[0][0].Desc, 2)
half = ncell/2

# Loop over frequency bins
oplane = [1,1,1,1,1]       # Output plane
for i in range (0,nFcomb):
    iplane = [1,1,1,1,1]   # Input plane
    # Loop over IF
    for iIF in range(0,numIF):
        if iIF>=numIF:
            break     #WHAT???  Looping doesn't appear to work correctly
        # Loop over channels
        for iChan in range (0,nChan):
           # loop over resolutions
            ires = 0
            FArray.PFill(accW, 0.0)   # Zero accumulators
            FArray.PFill(accWI, 0.0)  # Zero accumulators
            for x in dataImag[i]:
                numIF = x.Desc.Dict['inaxes'][3];
                # Read input plane
                x.GetPlane(None, iplane, err)
                inArr = x.FArray
                # Input pixels
                ImageUtil.PGetXYPixels(x, outIm, XPixIm, YPixIm, err)
                xpix = XPixIm.FArray; ypix = YPixIm.FArray; 
                # Set interpolator to input image
                FInterpolate.PSetDesc(fi, x.Desc)
                FInterpolate.PReplace(fi, inArr)
                # Interpolate
                for ix in range(0,ncell):
                    for iy in range (0,ncell):
                        pos = [space*(ix-half), space*(iy-half)]
                        pix = ImageDesc.PGetPixel(x.Desc,pos,err)
                        val = FInterpolate.PPixel(fi, pix, err)
                        work.set(val, ix,iy)
                # end interpolation loops
                # Accumulate weight image, zero where image blanked
                FArray.PFill(work2, weight[ires])
                FArray.PBlank(work2, work, work2)
                FArray.PDeblank(work2, 0.0)  # Replace any blanking with 0.0
                FArray.PAdd(accW, work2, accW)
                # Image * wt
                FArray.PDeblank(work, 0.0)  # Replace any blanking with 0.0
                FArray.PSMul(work, weight[ires])
                # Accumulate image*wt
                FArray.PAdd(accWI, work, accWI)
                ires += 1;  # Next resolution
            
            # end resolution loop
            FArray.PDiv(accWI, accW, accWI)  # Normalize to accWI
            # Write output
            outIm.PutPlane(accWI, oplane, err)
            iplane[0] += 1  # Next input plane
            if iplane[0]>nChan:
                iplane[0] = 1; iplane[1] += 1;
            oplane[0] += 1  # Next output plane
            if oplane[0]>nChan:
                oplane[0] = 1; oplane[1] += 1;
        # end channel loop
    # End IF loop
# end frequency bin loop

# copy history from first image
inHistory  = History.History("history", dataImag[0][0].List, err)
outHistory = History.History("history", outIm.List, err)
History.PCopyHeader(inHistory, outHistory, err)
# Add this programs history
outHistory.Open(History.READWRITE, err)
outHistory.TimeStamp("CombineBeam",err)
i = 0
for freqRoots in dataRoot:
    ires = 0
    for root in freqRoots:
        f = root+'.'+stoke+'All.fits'
        outHistory.WriteRec(-1,ObitSys.pgmName+' inFile['+str(i)+'] = '+f,err)
        f = ' weight['+str(i)+'] = '+str(weight[ires])+' \  weight'
        outHistory.WriteRec(-1,ObitSys.pgmName+f,err)
        ires += 1; i += 1
f = outRoot+'.'+stoke+'.fits'
outHistory.WriteRec(-1,ObitSys.pgmName+" outFile = "+f,err)
outHistory.Close(err)
OErr.printErrMsg(err, "Error with history")

# Cleanup
del outIm, work, work2, accWI, accW
