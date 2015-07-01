# Spectrum plotting  utilities
import math, string, Obit, Image, ImageUtil, ImageDesc, FArray, FArrayUtil, InfoList
import Table, OErr
import OPlot

def ParseSpectrum(path, bchan=5, echan=124, lob=5, loe=20, hib=110, hie=124, fact=1.0):
    """ Read Total intensity spectrum from an edited POSSM file

    A linear baseline is fitted and subtracted and the desired channels selected
    First line should contain the values, bchan echan, lob, loe, hib, hie
    where bchan, echan are the desired 1-rel channel numbers
    lob, loe are the channel numbers for the "off" for low end of the spectrum
    hib, hie are the channel numbers for the "off" for high end of the spectrum
    returns a dictionary with entries
       ch   = channel numbers
       vel  = velocity (km/s) coordinate (first col)
       flux = spectral values
       raw  = unclear
    path  = file path
    bchan, echan = the desired 1-rel channel numbers
    lob, loe = the channel numbers for the "off" for low end of the spectrum
    hib, hie = the channel numbers for the "off" for high end of the spectrum
    fact = flux scaling factor
    """
    ################################################################
    input=open(path)
    # Spool through file until line starting with "Channel  IF "
    line = input.readline()
    while (line):
        if line.__contains__("Channel "):
            # Find parts
            parts=line.split()
            ivelo=3
            i= 0
            for p in parts:
                if p=="Velocity":
                    ivelo=i; break
                i = i+1
            iflux=4
            i= 0
            for p in parts:
                if p.startswith("Real"):
                    iflux=i; break
                i = i+1
            break
        line = input.readline()
    
    # parse rest of file
    ich=[]; ivel=[]; ir=[];
    line = input.readline()
    while (line):
        if line.__contains__("FLAGGED"):  # Ignore bad channels
            line = input.readline()
            continue
        parts=line.split()
        ich.append(int(parts[0]))
        ivel.append(float(parts[ivelo]))        # Velocity
        ir.append(float(parts[iflux])*fact)    # Flux
        if int(parts[0])==128:  # Done?
            break
        line = input.readline()
    
    input.close()

    # Fit linear baseline
    sumxh=0.0; sumxl=0.0; sumyh = 0.0; sumyl=0.0; countl=0; counth=0
    for i in range(0,len(ich)):
        if ((ich[i]>=lob) and (ich[i]<=loe)):
            countl = countl+1
            iv = ivel[i]-ivel[0]
            sumxl  = sumxl  + iv
            sumyl  = sumyl  + ir[i]
        if ((ich[i]>=hib) and (ich[i]<=hie)):
            counth = counth+1
            iv = ivel[i]-ivel[0]
            sumxh  = sumxh  + iv
            sumyh  = sumyh  + ir[i]
        
    xlo=sumxl/countl; ylo=sumyl/countl
    xhi=sumxh/counth; yhi=sumyh/counth
    slope = (yhi-ylo)/(xhi-xlo)
    r0 = ylo - (xlo) * slope
    #print "slope",slope,"intercept",r0
    # remove linear baseline
    for i in range(0,len(ich)):
        ir[i] = ir[i] - r0 - slope * (ivel[i]-ivel[0])
    
    # Select
    ch=[]; vel=[];flux=[]
    for i in range(0,len(ich)):
        if ((ich[i]>=bchan) and (ich[i]<=echan)):
            ch.append(ich[i])
            vel.append(ivel[i])
            flux.append(ir[i])
    
    return {'ch':ch,'vel':vel,'flux':flux,'raw':ir}
    # end ParseSpectrum
    
def SumCCSpectrum(path, err, chOff=-9999):
    """ Obtain a spectrum from the CC tables in an FITS image

    Sum all the clean components in each table (one per channel)
    returns a dictionary with entries
       ch   = channel numbers
       vel  = velocity (km/s) coordinate (first col)
       flux = spectral values
    path  = FITS file path, assume FITS disk 0
    chOff = offset to add to channel number in computing Vel
            if not given is determined from image header
    """
    ################################################################
    inImage   = Image.newPFImage("Input image", path, 0, True, err)
    OErr.printErrMsg(err, "Error making image")
    imDict = inImage.Desc.Dict
    ntable = inImage.GetHighVer("AIPS CC")
    inImage.Open(Image.READONLY,err)
    OErr.printErrMsg(err, "Error making image")
    if chOff<-9000:
        d=inImage.Desc.Dict
        chOff = 1-d["crpix"][d["jlocf"]]
    
    # Loop over tables
    ch = []; vel = []; flux=[]
    for i in range(1,ntable+1):
        # Make Table
        inTable = inImage.NewTable (Table.READONLY, "AIPS CC", i, err)
        OErr.printErrMsg(err, "Eror making table")
        inTable.Open(Table.READONLY,err)
        OErr.printErrMsg(err, "Error with table")
        nrow = inTable.Desc.Dict["nrow"]
        sum = 0.0;
        for irow in range (1,nrow+1):
            row=inTable.ReadRow(irow,err)
            OErr.printErrMsg(err, "Error with table")
            sum = sum + row['FLUX'][0]
        inTable.Close(err)
        OErr.printErrMsg(err, "Error with table")
        #print "debug", i, nrow, sum
        # save
        ch.append(i)
        v = 0.001*(imDict['crval'][2]+imDict['cdelt'][2]*(i+chOff-imDict['crpix'][2]))
        vel.append(v)
        flux.append(sum)
    
    inImage.Close(err)
    OErr.printErrMsg(err, "Error closinging image")
    return {'ch':ch,'vel':vel,'flux':flux}
    # end SumCCSpectrum
    

def SumImageSpectrum(path, err, chOff=-9999):
    """ Obtain a spectrum from integral of image planes in an FITS image

    Sum all pixel values and divides by the beam area
    returns a dictionary with entries
       ch   = channel numbers
       vel  = velocity (km/s) coordinate (first col)
       flux = spectral values
    path = FITS file path, assume FITS disk 0
    chOff = offset to add to channel number in computing Vel
            if not given is determined from image header
    """
    ################################################################
    inImage   = Image.newPFImage("Input image", path, 0, True, err)
    OErr.printErrMsg(err, "Error making image")
    imDict = inImage.Desc.Dict
    nplane = imDict["inaxes"][2]
    blc = [1,1,1]; trc=[imDict["inaxes"][0], imDict["inaxes"][1], 1]
    beamMaj = imDict["beamMaj"] / abs(imDict["cdelt"][0])
    beamMin = imDict["beamMin"] / abs(imDict["cdelt"][1])
    barea = 1.1331 * beamMaj * beamMin;
    inImage.Open(1,err)
    OErr.printErrMsg(err, "Error making image")
    if chOff<-9000:
        d=inImage.Desc.Dict
        chOff = 1-d["crpix"][d["jlocf"]]

    # Loop over planes
    ch = []; vel = []; flux=[]
    for i in range(1,nplane+1):
        # Read plane
        blc[2]=i; trc[2]=i
        data = inImage.ReadPlane(err,blc=blc,trc=trc)
        OErr.printErrMsg(err, "Error reading plane")
        s = data.Sum
        
        #print "debug", i, s/barea
        # save
        ch.append(i)
        v = 0.001*(imDict['crval'][2]+imDict['cdelt'][2]*(i+chOff-imDict['crpix'][2]))
        vel.append(v)
        flux.append(s/barea)
    
    inImage.Close(err)
    OErr.printErrMsg(err, "Error closinging image")
    return {'ch':ch,'vel':vel,'flux':flux}
    # end SumImageSpectrum
    

def ImageSpectrum(path, pixel, err, chOff=-9999,bchan=1, echan=0):
    """ Obtain a spectrum from a given pixel in an FITS image

    returns a dictionary with entries
       ch   = channel numbers
       vel  = velocity (km/s) coordinate (first col)
       freq = frequency (Hz)
       flux = spectral values
    path = FITS file path, assume FITS disk 0
    pixel = array of integers giving 1-rel pixel
    chOff = offset to add to channel number in computing Vel
            if not given is determined from image header
    bchan = first 1-rel plane number
    echan = highest plane number 0->all
    """
    ################################################################
    inImage   = Image.newPFImage("Input image", path, 0, True, err)
    OErr.printErrMsg(err, "Error making image")
    imDict = inImage.Desc.Dict
    #print imDict
    nplane = imDict["inaxes"][2]
    blc = [pixel[0],pixel[1],1];
    trc=[pixel[0], imDict["inaxes"][1], 1]
    nchan = echan
    if echan<=0:
        echan = imDict["inaxes"][2]
    inImage.Open(1,err)
    OErr.printErrMsg(err, "Error making image")
    if chOff<-9000:
        d=inImage.Desc.Dict
        chOff = 1-d["crpix"][d["jlocf"]]

    # Loop over planes
    ch = []; vel = []; freq=[]; flux=[]
    for i in range(bchan,echan+1):
        # Read plane
        blc[2]=i; trc[2]=i
        data = inImage.ReadPlane(err,blc=blc,trc=trc)
        OErr.printErrMsg(err, "Error reading plane")
        s = data.get(0,0)
        
        # save
        ch.append(i)
        v = 0.001*(imDict['crval'][2]+imDict['cdelt'][2]*(i+chOff-imDict['crpix'][2]))
        vel.append(v)
        freq.append(v)
        flux.append(s)
        #print "debug", i, s, v
    
    inImage.Close(err)
    OErr.printErrMsg(err, "Error closinging image")
    return {'ch':ch,'vel':vel,'freq':freq,'flux':flux}
    # end ImageSpectrum
    
def OImageSpectrum(inImage, pixel, err, chOff=-9999,bchan=1, echan=0):
    """ Obtain a spectrum from a given pixel in an Obit image

    returns a dictionary with entries
       ch   = channel numbers
       vel  = velocity (km/s) coordinate (first col)
       flux = spectral values
    inImage = Obit image
    pixel   = array of integers giving 1-rel pixel
    chOff   = offset to add to channel number in computing Vel
              if not given is determined from image header
    bchan   = first 1-rel plane number
    echan   = highest plane number 0->all
    """
    ################################################################
    imDict = inImage.Desc.Dict
    # save info
    nx       = imDict["inaxes"][0]
    ny       = imDict["inaxes"][1]
    nplane   = imDict["inaxes"][2]
    altCrpix = imDict["altCrpix"]
    blc = [pixel[0],pixel[1],1];
    trc=[pixel[0], imDict["inaxes"][1], 1]
    trc=[pixel[0], pixel[1], 1]
    if echan<=0:
        echan = imDict["inaxes"][2]
    nchan = echan
    inImage.Open(1,err)
    OErr.printErrMsg(err, "Error making image")
    if chOff<-9000:
        d=inImage.Desc.Dict
        chOff = 1-d["crpix"][d["jlocf"]]

    # Loop over planes
    ch = []; vel = []; flux=[]
    for i in range(bchan,nchan+1):
        # Read plane
        blc[2]=i; trc[2]=i
        data = inImage.ReadPlane(err,blc=blc,trc=trc)
        OErr.printErrMsg(err, "Error reading plane")
        s = data.get(0,0)
        
        #print "debug", i, blc, trc, s
        # save
        ch.append(i)
        v = 0.001*(imDict['crval'][2]+imDict['cdelt'][2]*(i+chOff-imDict['crpix'][2]))
        vel.append(v)
        flux.append(s)
    
    # restore info
    imDict["inaxes"][0] = nx; imDict["inaxes"][1] = ny; imDict["inaxes"][2] = nplane
    imDict["altCrpix"]  = altCrpix
    inImage.Desc.Dict   = imDict
    inImage.Close(err)
    OErr.printErrMsg(err, "Error closinging image")
    return {'ch':ch,'vel':vel,'flux':flux}
    # end OImageSpectrum
    

def OImageMFSpectrum(inImage, pixel, err, bchan=1, echan=0):
    """ Obtain a spectrum from a given pixel in an Obit MF image

    returns a dictionary with entries
       ch      = channel numbers
       freq    = velocity (MHz) coordinate
       flux    = spectral values (mJy)
       refFreq = reference frequency
       fitSpec = fitted spectrum
    inImage = Obit MF image (as from MFImage)
    pixel   = array of integers giving 1-rel pixel
    bchan   = first 1-rel plane number (relative to cube)
    echan   = highest plane number 0->all
    """
    ################################################################
    # Frequency info
    imDict = inImage.Desc.List.Dict
    nterm = imDict["NTERM"][2][0]
    nspec = imDict["NSPEC"][2][0]
    hfreq = []
    for i in range(0,nspec):
        key = "FREQ%4.4d"%(i+1)
        hfreq.append(imDict[key][2][0])
    # Image
    info = inImage.List
    info.set("BLC",[1,1,1,1,1,1,1])
    info.set("TRC",[0,0,0,0,0,0,0])
    inImage.Open(Image.READONLY, err)
    OErr.printErrMsg(err, "Error opening image")
    imDict = inImage.Desc.Dict
    nplane = imDict["inaxes"][2]
    blc = [pixel[0],pixel[1],1];
    trc=[pixel[0], imDict["inaxes"][1], 1]
    trc=[pixel[0], pixel[1], 1]
    nchan = echan
    if echan<=0:
        echan = imDict["inaxes"][2]
        nchan = echan
    refFreq = imDict["crval"][2]

    # Loop over planes reading fittes spectrum
    fitSpec = []
    for i in range(1,nterm+1):
        # Read plane
        blc[2]=i; trc[2]=i
        data = inImage.ReadPlane(err,blc=blc,trc=trc)
        OErr.printErrMsg(err, "Error reading plane")
        fitSpec.append(data.get(0,0))
    
    # Loop over planes reading channel fluxes
    ch = []; freq = []; flux=[]
    for i in range(bchan+nterm,nchan+1):
        # Read plane
        blc[2]=i; trc[2]=i
        data = inImage.ReadPlane(err,blc=blc,trc=trc)
        OErr.printErrMsg(err, "Error reading plane")
        s = data.get(0,0)
        
        #print "debug", i, blc, trc, s
        # save
        ch.append(i)
        freq.append(hfreq[i-nterm-1]*1.0e-9)
        flux.append(s*1.0e3)
    
    inImage.Close(err)
    OErr.printErrMsg(err, "Error closinging image")
    bchan = 1; echan = 0;
    return {'ch':ch,'freq':freq,'flux':flux,'refFreq':refFreq,'fitSpec':fitSpec}
    # end OImageMFSpectrum
    

def PlotInp(spec, spec2, label, file, \
            xlabel="Velocity (km/sec)", ylabel="Flux Density (Jy)",
            xrange=None, yrange=None):
    """ Creates pgPlot input for a spectrum

    accepts spectrum in a dictionary with entries
       vel  = x coordinate
       flux = y coordinate
    spec   = first input spectrum, plotted with solid line and used to scale
    spec2  = second input spectrum, plotted with *, None = ignore
    label  = plot label
    file   = name of pgPlot input file
    xlabel = [optional] x axis label
    ylabel = [optional] y axis label
    xrange = [optional] range of values in x, None = Use range
    yrange = [optional] range of values in y, None = Use range
    """
    ################################################################
    if not xrange:
        xrange=[0.0,0.0]
        mx = max(spec["vel"])
        mn = min(spec["vel"])
        xrange[0] = mn - 0.03*(mx-mn)
        xrange[1] = mx + 0.03*(mx-mn)
        xpos = xrange[0] +0.05*(mx-mn)
    if not yrange:
        yrange=[0.0,0.0]
        mx = max(spec["flux"])
        mn = min(spec["flux"])
        yrange[0] = mn - 0.1*(mx-mn)
        yrange[1] = mx + 0.1*(mx-mn)
        ypos = yrange[1] -0.10*(mx-mn)
    pfile = open(file,"w")
    pfile.write("      \n")
    pfile.write(xlabel+"\n")
    pfile.write(ylabel+"\n")
    pfile.write("    1\n")
    pfile.write("(2F10.5,I5,1X,A)\n")
    ss="%8.3f"%xpos
    xo = string.rjust(ss,10)
    ss="%8.3f"%ypos
    yo = string.rjust(ss,10)
    llabel = label.replace("_"," ")  # no underscore
    pfile.write(xo+yo+"    1 "+llabel+"\n")
    pfile.write("(2F20.5,2I5)\n")
    xmin=string.rjust(str(xrange[0]),20)
    xmax=string.rjust(str(xrange[1]),20)
    ymin=string.rjust(str(yrange[0]),20)
    ymax=string.rjust(str(yrange[1]),20)
    pfile.write(xmin+xmax+ymin+ymax+"\n")
    pfile.write("    2    2    2\n")
    # First spectrum 
    i = 0
    sym = 1
    x = spec["vel"]
    y = spec["flux"]
    for xx in x:
        if y[i]>-99990.0:
            xo = string.rjust(str(xx),20)
            yo = string.rjust(str(y[i]),20)
            lo = string.rjust(str(1),5)
            so = string.rjust(str(sym),5)
            pfile.write(xo+yo+lo+so+"\n")
        i = i+1
    
    # Second spectrum, if given
    if not spec2:
        pfile.close()
        return
    i = 0
    sym = 3
    x = spec2["vel"]
    y = spec2["flux"]
    for xx in x:
        if y[i]>-99990.0:
            xo = string.rjust(str(xx),20)
            yo = string.rjust(str(y[i]),20)
            lo = string.rjust(str(0),5)
            so = string.rjust(str(sym),5)
            pfile.write(xo+yo+lo+so+"\n")
            # error bar
            #yo = string.rjust(str(y[i]-err[i]),20)
            #lo = string.rjust(str(0),5)
            #so = string.rjust(str(1),5)
            #pfile.write(xo+yo+lo+so+"\n")
            #yo = string.rjust(str(y[i]+err[i]),20)
            #lo = string.rjust(str(1),5)
            #so = string.rjust(str(1),5)
            #pfile.write(xo+yo+lo+so+"\n")
        i = i+1
    
    pfile.close()
    # end PlotInp

def doPlot(spec, spec2, label, file="None", \
           xlabel="Velocity (km/s)", ylabel="Flux Density (Jy)", \
           lwidth=1,csize=1,symbol=3):
    """ Plots and displays a spectrum

    accepts spectrum in a dictionary with entries
       vel = vel coordinate
       flux = flux
    spec   = input spectrum to plot with solid line
    spec2  = input spectrum to plot with symbol
    label  = plot label
    file   = name of output postscript file, if "None" ask
             "/xterm" gives xterm display
    xlabel = [optional] x axis label
    ylabel = [optional] y axis label
    lwidth = [optional] line width (integer)
    csize  = [optional] symbol size (integer)
    symbol = [optional] plot symbol code
    """
    ################################################################
    xmin=min(min(spec['vel']),min(spec2['vel']))
    xmax=max(max(spec['vel']),max(spec2['vel']))
    ymin=1.05*min(min(spec['flux']),min(spec2['flux']))
    ymax=1.1*max(max(spec['flux']),max(spec2['flux']))
    # Add "/PS" for postscript files
    if file!="None" and file!="/xterm":
        filename = file+"/ps"
    else:
        filename = file
    # Create plot object
    err=Image.err
    plot=OPlot.newOPlot("plot", err,output=filename)
    # Set labeling on plot InfoList member
    info = plot.List
    dim = [1,1,1,1,1]
    InfoList.PPutFloat   (info, "XMAX", dim, [xmax],err)
    InfoList.PPutFloat   (info, "XMIN", dim, [xmin],err)
    InfoList.PPutFloat   (info, "YMAX", dim, [ymax],err)
    InfoList.PPutFloat   (info, "YMIN", dim, [ymin],err)
    InfoList.PPutInt     (info, "LWIDTH", dim, [lwidth],err)
    InfoList.PPutInt     (info, "CSIZE",  dim, [csize],err)
    dim[0] = max (1,len(label))
    InfoList.PAlwaysPutString  (info, "TITLE", dim, [label])
    dim[0] = max (1,len(xlabel))
    InfoList.PAlwaysPutString  (info, "XLABEL",dim, [xlabel])
    dim[0] =  max (1,len(ylabel))
    InfoList.PAlwaysPutString  (info, "YLABEL",dim, [ylabel])
    x   = spec["vel"]
    y   = spec["flux"]
    OPlot.PXYPlot (plot, -1, x, y, err)
    x   = spec2["vel"]
    y   = spec2["flux"]
    OPlot.PXYOver (plot, symbol, x, y, err)
    OPlot.PShow  (plot, err)
    # end doPlot

def doSimplePlot(spec, label, file="None", \
           xlabel="Velocity (km/s)", ylabel="Flux Density (Jy)", \
           lwidth=1,csize=1,symbol=3):
    """ Plots and displays a simple spectrum

    accepts spectrum in a dictionary with entries
       vel = vel coordinate
       flux = flux
    spec   = input spectrum to plot with symbol
    label  = plot label
    file   = name of output postscript file, if "None" ask
             "/xterm" gives xterm display
    xlabel = [optional] x axis label
    ylabel = [optional] y axis label
    lwidth = [optional] line width (integer)
    csize  = [optional] symbol size (integer)
    symbol = [optional] plot symbol code
    """
    ################################################################
    xmin=0.95*min(spec['vel'])
    xmax=1.05*max(spec['vel'])
    ymin=0.95*min(spec['flux'])
    ymax=1.05*max(spec['flux'])
    # Add "/PS" for postscript files
    if file!="None" and file!="/xterm":
        filename = file+"/ps"
    else:
        filename = file
    # Create plot object
    err=Image.err
    plot=OPlot.newOPlot("plot", err,output=filename)
    # Set labeling on plot InfoList member
    info = plot.List
    dim = [1,1,1,1,1]
    InfoList.PPutFloat   (info, "XMAX", dim, [xmax],err)
    InfoList.PPutFloat   (info, "XMIN", dim, [xmin],err)
    InfoList.PPutFloat   (info, "YMAX", dim, [ymax],err)
    InfoList.PPutFloat   (info, "YMIN", dim, [ymin],err)
    InfoList.PPutInt     (info, "LWIDTH", dim, [lwidth],err)
    InfoList.PPutInt     (info, "CSIZE",  dim, [csize],err)
    dim[0] = max (1,len(label))
    InfoList.PAlwaysPutString  (info, "TITLE", dim, [label])
    dim[0] = max (1,len(xlabel))
    InfoList.PAlwaysPutString  (info, "XLABEL",dim, [xlabel])
    dim[0] =  max (1,len(ylabel))
    InfoList.PAlwaysPutString  (info, "YLABEL",dim, [ylabel])
    x   = spec["vel"]
    y   = spec["flux"]
    OPlot.PXYPlot (plot, symbol, x, y, err)
    OPlot.PShow  (plot, err)
    del plot
    # end doSimplePlot

def doFourPlot(specs, labels, file="None", \
           xlabel="Velocity (km/s)", ylabel="Flux Density (Jy)", \
           lwidth=1,csize=1,symbol=3):
    """ Plots and displays a simple spectrum

    accepts spectrum in a dictionary with entries
       vel = vel coordinate
       flux = flux
    specs  = list of input spectra (up to four) to plot with symbol
    labels = plot label list
    file   = name of output postscript file, if "None" ask
             "/xterm" gives xterm display
    xlabel = [optional] x axis label
    ylabel = [optional] y axis label
    lwidth = [optional] line width (integer)
    csize  = [optional] symbol size (integer)
    symbol = [optional] plot symbol code
    """
    ################################################################
    # Add "/PS" for postscript files
    if file!="None" and file!="/xterm":
        filename = file+"/ps"
    else:
        filename = file
    # Create plot object
    err=Image.err
    plot=OPlot.newOPlot("plot", err,output=filename, nx=2,ny=2)
    # Loop over plots
    for i in range(0,len(specs)):
        spec = specs[i]; label = labels[i]
        xmin=0.95*min(spec['vel'])
        xmax=1.05*max(spec['vel'])
        ymin=0.95*min(spec['flux'])
        ymax=1.05*max(spec['flux'])
    # Set labeling on plot InfoList member
        info = plot.List
        dim = [1,1,1,1,1]
        InfoList.PPutFloat   (info, "XMAX", dim, [xmax],err)
        InfoList.PPutFloat   (info, "XMIN", dim, [xmin],err)
        InfoList.PPutFloat   (info, "YMAX", dim, [ymax],err)
        InfoList.PPutFloat   (info, "YMIN", dim, [ymin],err)
        InfoList.PPutInt     (info, "LWIDTH", dim, [lwidth],err)
        InfoList.PPutInt     (info, "CSIZE",  dim, [csize],err)
        dim[0] = max (1,len(label))
        InfoList.PAlwaysPutString  (info, "TITLE", dim, [label])
        dim[0] = max (1,len(xlabel))
        InfoList.PAlwaysPutString  (info, "XLABEL",dim, [xlabel])
        dim[0] =  max (1,len(ylabel))
        InfoList.PAlwaysPutString  (info, "YLABEL",dim, [ylabel])
        x   = spec["vel"]
        y   = spec["flux"]
        OPlot.PXYPlot (plot, symbol, x, y, err)
        # end loop over plots
    OPlot.PShow  (plot, err)
    del plot
    # end doFourPlot

def doSimpleLogPlot(spec, label, file="None", \
           xlabel="Freq (MHz)", ylabel="Flux Density (mJy)", \
           lwidth=1,csize=1,symbol=3):
    """ Plots and displays a simple log=log spectrum

    accepts spectrum in a dictionary with entries
       freq    = freq coordinate (MHz)
       flux    = flux (mJy)
       refFreq = reference frequency
       fitSpec = fitted spectrum, plotted as line if given
    spec   = input spectrum to plot with symbol
    label  = plot label
    file   = name of output postscript file, if "None" ask
             "/xterm" gives xterm display
    xlabel = [optional] x axis label
    ylabel = [optional] y axis label
    lwidth = [optional] line width (integer)
    csize  = [optional] symbol size (integer)
    symbol = [optional] plot symbol code
    """
    ################################################################
    xmin=0.95*min(spec['freq'])
    xmax=1.05*max(spec['freq'])
    ymin=0.95*min(spec['flux'])
    ymax=1.05*max(spec['flux'])
    # Add "/PS" for postscript files
    if file!="None" and file!="/xterm":
        filename = file+"/ps"
    else:
        filename = file
    # Create plot object
    err=Image.err
    plot=OPlot.newOPlot("plot", err,output=filename)
    # Set labeling on plot InfoList member
    info = plot.List
    dim = [1,1,1,1,1]
    InfoList.PPutFloat   (info, "XMAX", dim, [xmax],err)
    InfoList.PPutFloat   (info, "XMIN", dim, [xmin],err)
    InfoList.PPutFloat   (info, "YMAX", dim, [ymax],err)
    InfoList.PPutFloat   (info, "YMIN", dim, [ymin],err)
    InfoList.PPutInt     (info, "LWIDTH", dim, [lwidth],err)
    InfoList.PPutInt     (info, "CSIZE",  dim, [csize],err)
    dim[0] = max (1,len(label))
    InfoList.PAlwaysPutString  (info, "TITLE", dim, [label])
    dim[0] = max (1,len(xlabel))
    InfoList.PAlwaysPutString  (info, "XLABEL",dim, [xlabel])
    dim[0] =  max (1,len(ylabel))
    InfoList.PAlwaysPutString  (info, "YLABEL",dim, [ylabel])
    xopt = "BCLTS"
    dim[0] =  max (1,len(xopt))
    InfoList.PAlwaysPutString  (info, "XOPT",  dim, [xopt])
    yopt = "BCLTS"
    dim[0] =  max (1,len(yopt))
    InfoList.PAlwaysPutString  (info, "YOPT",  dim, [yopt])
    x   = spec["freq"]
    y   = spec["flux"]
    OPlot.PXYPlot (plot, symbol, x, y, err)
    # Fitted spectrum given?
    if ("fitSpec" in spec) and (spec["fitSpec"]!=None) and (len(spec["fitSpec"])>0):
        refFreq = spec["refFreq"]*1.0e-9
        ss = spec["fitSpec"]
        x = []; y = []
        deltax = xmax - xmin
        for i in range (0,101):
            f = xmin + i*0.01*deltax
            x.append(f)
            ll  = math.log(f/refFreq)
            lll = math.log(f/refFreq)
            arg = 0.0
            for t in ss[1:]:
                arg += t*lll
                lll *= ll
            y.append(ss[0] * 1000.0 * math.exp(arg))
        OPlot.PXYOver (plot, 0, x, y, err)
        if len(ss)>=3:
            speclabel="S#d%5.1f GHz#u=%5.3f mJy, #ga= %5.2f, #gb=%5.3f"%(refFreq,ss[0]*1000.0,ss[1],ss[2])
        elif len(ss)==2:
            speclabel="S#d%5.1f GHz#u=%5.3f mJy, #ga= %5.2f"%(refFreq,ss[0]*1000.0,ss[1])
        OPlot.PRelText(plot, "T", -01.0, 0.3, 0.0, speclabel, err)

        
    labelLogPlot (plot, xmin, xmax, ymin, ymax, err)
    OPlot.PShow  (plot, err)
    del plot
    # end doSimpleLogPlot

def doVPlot(spec, spec2, label, fact2=1.0, file="None", \
           xlabel="Velocity (km/s)", ylabel="Flux Density (Jy)", \
           lwidth=2,csize=1,symbol=3, Stokes="V"):
    """ Plots and displays an I and Stokes Stokes spectrum 

    accepts spectrum in a dictionary with entries
       vel = vel coordinate
       flux = flux
    spec   = IPol spectrum
    spec2  = VPol spectrum to plot with solid line, scaled on return
    fact2  = Factor to multiple times spec2 flux
    label  = plot label
    file   = name of output postscript file, if "None" ask
             "/xterm" gives xterm display
    xlabel = [optional] x axis label
    ylabel = [optional] y axis label
    lwidth = [optional] line width (integer)
    csize  = [optional] symbol size (integer)
    symbol = [optional] plot symbol code
    """
    ################################################################
    # Scale second spectrum
    if fact2!=1.0:
        for i in range(0,len(spec2["flux"])):
            spec2["flux"][i] *= fact2
        
    xmin=min(min(spec['vel']),min(spec2['vel']))
    xmax=max(max(spec['vel']),max(spec2['vel']))
    ymin=1.05*min(min(spec['flux']),min(spec2['flux']))
    ymax=1.1*max(max(spec['flux']),max(spec2['flux']))
    # Add "/PS" for postscript files
    if file!="None" and file!="/xterm":
        filename = file+"/ps"
    else:
        filename = file
    # Create plot object
    err=Image.err
    plot=OPlot.newOPlot("plot", err,output=filename)
    # Set labeling on plot InfoList member
    info = plot.List
    dim = [1,1,1,1,1]
    InfoList.PPutFloat   (info, "XMAX", dim, [xmax],err)
    InfoList.PPutFloat   (info, "XMIN", dim, [xmin],err)
    InfoList.PPutFloat   (info, "YMAX", dim, [ymax],err)
    InfoList.PPutFloat   (info, "YMIN", dim, [ymin],err)
    InfoList.PPutInt     (info, "LWIDTH", dim, [lwidth],err)
    InfoList.PPutInt     (info, "CSIZE",  dim, [csize],err)
    dim[0] = max (1,len(label))
    InfoList.PAlwaysPutString  (info, "TITLE", dim, [label])
    dim[0] = max (1,len(xlabel))
    InfoList.PAlwaysPutString  (info, "XLABEL",dim, [xlabel])
    dim[0] =  max (1,len(ylabel))
    InfoList.PAlwaysPutString  (info, "YLABEL",dim, [ylabel])
    x   = spec["vel"]
    y   = spec["flux"]
    OPlot.PXYPlot (plot, -1, x, y, err)
    x   = spec2["vel"]
    y   = spec2["flux"]
    OPlot.PXYOver (plot, symbol, x, y, err)
    # Label
    tx = xmin + 0.02*(xmax-xmin)
    ty = ymax - 0.05*(ymax-ymin)
    OPlot.PText(plot, tx, ty, 0.0, 0.0, "Line = I", err)
    ty = ymax - 0.09*(ymax-ymin)
    OPlot.PText(plot, tx, ty, 0.0, 0.0, "* = "+Stokes+"*"+str(fact2), err)
    OPlot.PShow  (plot, err)
    # end doVPlot
    
def labelLogPlot(plot, xmin, xmax, ymin, ymax, err, vals=None):
    """ Label log-log plots - plPlot is completely incompetent at this

    plot   = plot to label
    xmin   = xmin of plot, normal units
    xmax   = xmax of plot, normal units
    ymin   = ymin of plot, normal units
    ymax   = ymax of plot, normal units
    vals   = list of strings of values to plot
    """
    ################################################################
    # Default values of labels
    defvals = ["0.002", "0.003", "0.004", "0.005", "0.007", \
               "0.02", "0.03", "0.04", "0.05", "0.07", \
               "0.2", "0.3", "0.4", "0.5", "0.7", \
               "2", "3", "4", "5", "6", "7", "8", "9", \
               "20", "30", "40", "50", "70", \
               "200", "300", "400", "500", "700", \
               "2000", "3000", "4000", "5000", "7000", \
               ]
    if vals==None:
        pltVals = defvals
    else:
        pltVals = vals

    # logs of extrema
    logxmin = math.log10(xmin)
    logxmax = math.log10(xmax)
    logymin = math.log10(ymin)
    logymax = math.log10(ymax)
    
    # X axis
    fjust = 0.5
    for val in pltVals:
        num = float(val)
        if (num>=xmin) and (num<=xmax):
            disp  = len(val)
            coord = (math.log10(num)-logxmin) / (logxmax - logxmin)
            OPlot.PRelText (plot, "B", disp, coord, fjust, val, err)
    # end loop
    # Y axis
    fjust = 0.0
    for val in pltVals:
        num = float(val)
        if (num>=ymin) and (num<=ymax):
            disp  = len(val)
            coord = (math.log10(num)-logymin) / (logymax - logymin)
            OPlot.PRelText (plot, "LV", disp, coord, fjust, val, err)
    # end loop

def FSSpectrum(fstab, entry, err, bchan=1, echan=0):
    """ Obtain a spectrum from a given entry in an FS table

    returns a dictionary with entries
       ch   = channel numbers
       vel  = velocity (km/s) coordinate (first col)
       freq = frequency (Hz)
       flux = spectral values
    fstab = FS table
    entry = FS table entry Id
    bchan = first 1-rel plane number
    echan = highest plane number 0->all
    """
    ################################################################
    fstab.Open(Table.READONLY,err)
    fsrow = fstab.ReadRow(entry, err)
    OErr.printErrMsg(err, "Error reading Table")
    VelRef  = fstab.Desc.List.Dict['VEL_REF'][2][0]
    VelRPix = fstab.Desc.List.Dict['VEL_RPIX'][2][0]
    VelDelt = fstab.Desc.List.Dict['VEL_DEL'][2][0]
    nchan = echan
    if echan<=0:
        echan = fstab.Desc.List.Dict['NO_CH'][2][0]
    bchan = max (0, bchan)
    echan = min (echan, fstab.Desc.List.Dict['NO_CH'][2][0]-1)
    # Loop over channels
    ch = []; vel = []; freq=[]; flux=[]
    for i in range(bchan,echan+1):
        # save
        ch.append(i)
        v = 1.0e-3*(VelRef + (i-VelRPix)*VelDelt)  # Velocity
        s = fsrow['SPECTRUM'][i-1]                 # Spectral value
        vel.append(v)
        freq.append(v)
        flux.append(s)
        #print "debug", i, s, v
    fstab.Close(err)
    OErr.printErrMsg(err, "Error closinging image")
    return {'ch':ch,'vel':vel,'freq':freq,'flux':flux}
    # end FSSpectrum
    
