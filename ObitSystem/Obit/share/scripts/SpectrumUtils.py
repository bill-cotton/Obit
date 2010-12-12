# Mira Sprectrum plotting  utilities
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
    ch = []; vel = []; flux=[]
    for i in range(bchan,nchan+1):
        # Read plane
        blc[2]=i; trc[2]=i
        data = inImage.ReadPlane(err,blc=blc,trc=trc)
        OErr.printErrMsg(err, "Error reading plane")
        s = data.get(0,0)
        
        #print "debug", i, s
        # save
        ch.append(i)
        v = 0.001*(imDict['crval'][2]+imDict['cdelt'][2]*(i+chOff-imDict['crpix'][2]))
        vel.append(v)
        flux.append(s)
    
    inImage.Close(err)
    OErr.printErrMsg(err, "Error closinging image")
    return {'ch':ch,'vel':vel,'flux':flux}
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
    nplane = imDict["inaxes"][2]
    blc = [pixel[0],pixel[1],1];
    trc=[pixel[0], imDict["inaxes"][1], 1]
    trc=[pixel[0], pixel[1], 1]
    nchan = echan
    if echan<=0:
        echan = imDict["inaxes"][2]
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
    
    inImage.Close(err)
    OErr.printErrMsg(err, "Error closinging image")
    return {'ch':ch,'vel':vel,'flux':flux}
    # end OImageSpectrum
    

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
    xmin=min(spec['vel'])
    xmax=max(spec['vel'])
    ymin=1.05*min(spec['flux'])
    ymax=1.1*max(spec['flux'])
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
    # end doSimplePlot

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
