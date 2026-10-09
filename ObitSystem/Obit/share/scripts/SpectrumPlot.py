# Plot a spectrum using matplotlib/numpy
#exec(open('SpectrumPlot.py').read())
boxSpec=None; linRegW=None; SpectrumPlot=None; ArrayMask=None; IsInWindow=None
pointSpec=None
import OSystem, Image, FArray, OWindow

del ArrayMask
def ArrayMask(inIm, winList, err):
    """
    Make masking FArray for inIm from winList

    * inIm     = Obit MFImage-like image
    * winList  = window list of regions in inIm to EXCLUDE
    * err      = Python Obit Error/message stack    
    """
    d=inIm.Desc.Dict
    naxis = d['inaxes'][0:2]
    mask = FArray.FArray('mask', naxis)
    # Fill with 1, not blanked
    FArray.PFill(mask, 1.0)
    blank = FArray.fblank
    # loop over winList
    for w in winList:
        #print (win)
        if w[1]==0: # Rectangle [(2,3),(4,5)]
            swin=[w[2],w[3],w[4],w[5]] # window in mask
            FArray.PRectFill(mask, swin, blank)
        else:       # Circle (rad=2, x=3, y=4]
            swin=[w[2],w[3],w[4]] # window in mask
            FArray.PRoundFill(mask, swin, blank)
    return mask
    # end ArrayMask

del IsInWindow
def IsInWindow(xpix, ypix, winList):
    """
    Determine if a pixel is in the window list
    
    Returns True or False
    * xpix    = x pixel (0-rel)
    * ypix    = y pixel (0-rel)
    * winList = window list, use CleanWindowEdit.py
    """
    # loop over winList
    for w in winList:
        #print (win)
        if w[1]==0: # Rectangle [(2,3),(4,5)]
            if (xpix>=w[2]) and (xpix<=w[4]) and (ypix>=w[3]) and (ypix<=w[5]):
                return True
        else:       # Circle (rad=2, x=3, y=4]
            del2 = (xpix-w[3])**2 + (ypix-w[4])**2 
            if del2<w[2]*w[2]:
                return True
    return False  # No match
# end IsInWindow

del boxSpec
def boxSpec (inIm, err,
             blc=[1,1], trc=[0,0], nThreads=1,  winList=None, debug=False):
    """
    Get box integral spectrum over an image

    Returns {"flux":[flux density], "ferr":[flux density error (Jy)],
             "freqs":[Frequencies (Hz)], "vals":[spectral vals (Jy)],
             "sig":[plane RMS (Jy)]}
    * inIm     = Obit MFImage-like image
    * err      = Python Obit Error/message stack
    * blc      = bottom left corner pixel in each plane (1-rel) 
    * trc      = top right corner pixel (1-rel)
    * nThreads = number of threads allowed
    * winList  = window list of regions in inIm to EXCLUDE
    * debug    = If True, write masked first plane as "debug.fits"
    """
    OSystem.PAllowThreads(nThreads)  # threading

    # Get subband frequencies
    inIm.Open(Image.READONLY,err)
    d=inIm.Desc.List.Dict
    nterm = d['NTERM'][2][0]
    nspec = d['NSPEC'][2][0]
    freqs = []
    for ip in range(nterm+1,nterm+nspec+1):
        key = "FREQ%4.4d"%(ip-nterm)
        freqs.append(d[key][2][0])
        # end get frequencies
    # Beam area
    d=inIm.Desc.Dict
    beamarea = 1.1331*(d["beamMaj"]/abs(d["cdelt"][0])) * \
                   (d["beamMin"]/abs(d["cdelt"][1]))

    # if winList given, convert to a mask
    if winList:
        mask = ArrayMask(inIm, winList, err)
    else:
        mask = None
    # Flux broadband density, error
    beamnorm = 1.0/beamarea
    plane = [1,1,1,1,1]
    inIm.GetPlane(None, plane, err)
    # if mask given, apply
    if mask:
        FArray.PBlank(inIm.FArray, mask,inIm.FArray)
    subarr = FArray.PSubArr(inIm.FArray, blc, trc, err) # Select
    # Debug dump to fits?
    if debug:
        Image.PFArray2FITS(inIm.FArray,"debug.fits",err,outDisk=0)
        print ("Wrote masked plane 1 to debug.fits")
    flux = beamnorm * subarr.Sum
    c    = subarr.Count
    ferr = subarr.RMS * ((c*beamnorm)**0.5)
    
    # Loop over subbands
    vals= []; sig = []
    plane = [1,1,1,1,1]
    for ip in range(nterm+1,nterm+nspec+1):
        plane[0] = ip
        inIm.GetPlane(None, plane, err)
        # if mask given, apply
        if mask:
            FArray.PBlank(inIm.FArray, mask,inIm.FArray)
        subarr = FArray.PSubArr(inIm.FArray, blc, trc, err) # Select
        # if mask given, apply
        v = subarr.Sum
        s = subarr.RMS
        c = subarr.Count
        vals.append(v*beamnorm)
        sig.append(s*((c*beamnorm)**0.5))
    inIm.Close(err)
    OErr.printErrMsg(err,message='Error reading data')
   
    return {"flux":flux, "ferr":ferr, "freqs":freqs, "vals":vals, "sig":sig}
    # end  boxSpec

del pointSpec
def pointSpec (inIm, pos, err, nThreads=1):
    """
    Get the spectrum for a pixel in inIm

    Returns {"flux":[flux density], "ferr":[flux density error (Jy)],
             "freqs":[Frequencies (Hz)], "vals":[spectral vals (Jy)],
             "sig":[plane RMS (Jy)]}
    * inIm     = Obit MFImage-like image
    * pos      = [xpixel, ypixel] 1 rel for spectrum
    * err      = Python Obit Error/message stack
    * nThreads = number of threads allowed
    """
    OSystem.PAllowThreads(nThreads)  # threading

    # Get subband frequencies
    inIm.Open(Image.READONLY,err)
    d=inIm.Desc.List.Dict
    nterm = d['NTERM'][2][0]
    nspec = d['NSPEC'][2][0]
    freqs = []
    for ip in range(nterm+1,nterm+nspec+1):
        key = "FREQ%4.4d"%(ip-nterm)
        freqs.append(d[key][2][0])
        # end get frequencies
    # Flux broadband density, error
    plane = [1,1,1,1,1]
    inIm.GetPlane(None, plane, err)
    flux = inIm.FArray.get(pos[0]-1, pos[1]-1)
    ferr = inIm.FArray.RMS
    
    # Loop over subbands
    vals= []; sig = []
    plane = [1,1,1,1,1]
    for ip in range(nterm+1,nterm+nspec+1):
        plane[0] = ip
        inIm.GetPlane(None, plane, err)
        v = inIm.FArray.get(pos[0]-1, pos[1]-1)
        s = inIm.FArray.RMS
        vals.append(v)
        sig.append(s)
    inIm.Close(err)
    OErr.printErrMsg(err,message='Error reading data')
   
    return {"flux":flux, "ferr":ferr, "freqs":freqs, "vals":vals, "sig":sig}
    # end  pointSpec

del linRegW 
def linRegW (nu, s, e):
    """
    Weighted linear regression for 2 term spectral fitting
    
    returns (flux@freq[0], spectral_index)
    * nu    = list of frequencies 
    * s     = list of flux densities
    * e     = list of flux density errors
    """
    from math import log, exp
    s_w=0; s_x=0.; s_y=0.; s_xx=0.; s_yy=0.; s_xy=0.
    nu_0 = nu[0]; s_0 = s[0]
    lnu = []; ls = []; w = []
    for i in range(0,len(nu)):
        if s[i] and s[i]>0:
            lnu.append(log(nu[i]/nu_0))
            ls.append(log(s[i]))
            w.append(1/e[i])
    n = len (lnu)
    for j in range(0,n):
        s_x += w[j]*lnu[j]; s_y += w[j]*ls[j]; s_xy += w[j]*lnu[j]*ls[j]
        s_xx += w[j]*lnu[j]**2; s_yy += w[j]*ls[j]**2; s_w+=w[j]
    a = (s_w*s_xy-s_x*s_y)/(s_w*s_xx-s_x*s_x)
    b = (s_y-a*s_x)/s_w
    #print (a, b, "S_0=", exp(b), "alpha=", a)
    return (exp(b),a)
# end linRegW

# Plotting stuff
try:
    import matplotlib
    import matplotlib.pyplot as plt
    import matplotlib.ticker as mticker

    class CleanLogFormatter(mticker.LogFormatterMathtext):
        def __call__(self, x, pos=None):
            # If the tick value is exactly 1 (which would be 10^0)
            if x == 1:
                return '1'
            # Otherwise, use Matplotlib's standard math formatting (10^1, 10^2, etc.)
            return super().__call__(x, pos)
        
except Exception as exception:
    print(exception)
    print ("Failed or No matplotlib")
# end CleanLogFormatter
        
# Band by color, range of subbands, color
del SpectrumPlot
def SpectrumPlot (freqs, vals, sig, plotfile, err,
                  title=None, xlabel="log Frequency (MHz)", ylabel="log Flux density (mJy)",
                  doFit=True,refFreq=None,bands=[["Observed",1,0,"red",1.0,0.05]]):
    """
    Make spectrum-like log-log plot
    * freqs   = Frequency array
    * vals    = Array of Spectral values
    * sig     = Array of Gaussian sigmas on vals, or None
    * plotfile= base name of plot file, ".pdf" added
    * err     = Python Obit Error/message stack
    * title   = title of plot, spectral index added if fitted
    * xlabel  = X axis label, should describe freqs
    * ylabel  = Y axis label, should describe vals
    * doFit   = If True, fit spectrum, plot
    * refFreq = If given, the reference frequency for the display of the doFit, fit
                [def freqs[0]]
    * bands   = list of ["name",first_ch, highest-ch, color scale_factor, cal. err]
                for different parts of the spectrum
    """
    import math
    from math import isnan, log
    # Fit Spectrum? Simple linear regression
    # Copy of vals
    vvals = []; 
    for v in vals:
        vvals.append(v)
    # Plot - use matplotlib
    nbands = len(bands)
    nch = len(freqs)
    # Be sure to include all
    if bands and (len(bands)>=1):
        bands[nbands-1][2] = max (bands[nbands-1][2],nch)
    try:
        import matplotlib
        import matplotlib.pyplot as plt
        import matplotlib.ticker as mticker

        fig, ax = plt.subplots()
        ssig = []
        for ib in range(0,nbands):
            #print ("band",bands[ib])
            # Scale
            c0 = max(0,bands[ib][1])-1; c1 = min(nch,bands[ib][2]); #print ("range",c0,c1)
            for ic in range(c0,c1):
                vvals[ic] *= bands[ib][4]
                if sig:  # Add cal error in quadrature
                    ssig.append(((sig[ic]**2)+(bands[ib][5]*vvals[ic]**2))**0.5)
            print ("band",bands[ib],"scaled by", bands[ib][4])
            if sig:
                ax.errorbar(freqs[c0:c1-1], vvals[c0:c1-1], yerr=ssig[c0:c1-1], color=bands[ib][3], fmt='o',label=bands[ib][0])
            else:
                ax.scatter(freqs[c0:c1-1], vvals[c0:c1-1], marker="*", mfc=bands[ib][3],label=bands[ib][0])
        # end nbands loop
        if doFit:
            if refFreq==None:
                refFreq = freqs[0] # Default
            FFitParms = linRegW (freqs, vvals, ssig)  # Ref is freqs[0]
            # reset reference of fit
            FitParms = (FFitParms[0]*(refFreq/freqs[0])**FFitParms[1], FFitParms[1])
            # Plot fitted spectrum
            x0 = min(freqs); y0 = FitParms[0]*(x0/refFreq)**FitParms[1]
            x1 = max(freqs); y1 = FitParms[0]*(x1/refFreq)**FitParms[1]
            ax.plot([x0,x1], [y0,y1],'g',label='Fit')  # Plot fitted spectrum
        # end doFit
        
        if title:
            print ("title",title,"file",plotfile)
            if doFit:
                ax.set_title(r" %s flux=%5.2f,$\alpha$=%5.2f"%(title,FitParms[0],FitParms[1]))
            else:
                ax.set_title(title)
        print ("fit flux=%5.2f,alpha=%5.2f"%(FitParms[0],FitParms[1]))
        ticks=[0.3e3,0.5e3,0.8e3,1.0e3,1.2e3,1.5e3,2.e3,5.0e3,8.0e3,10.0e3,12.0e3,15.0e3,20.0e3,50.0e3]
        ax.set_xlabel(xlabel); ax.set_xscale('log'); ax.set_xticks(ticks,minor=True)
        ax.set_ylabel(ylabel); ax.set_yscale('log')
        # From Dr. Google
        # 1. Force Matplotlib to place a specific number of ticks across your narrow range
        # (Adjust 'numticks' up or down to add or remove labels)
        ax.xaxis.set_major_locator(mticker.LogLocator(subs='all', numticks=6))
        ax.yaxis.set_major_locator(mticker.LogLocator(subs='all', numticks=5))

        # 2. Use standard string formatting ('%g') which dynamically formats 
        # numbers cleanly without scientific notation (e.g. 1, 6, 1000)
        ax.xaxis.set_major_formatter(mticker.FormatStrFormatter('%g'))
        ax.yaxis.set_major_formatter(mticker.FormatStrFormatter('%g'))

        # 3. Turn off the minor formatters completely so they don't leak labels
        ax.xaxis.set_minor_formatter(mticker.NullFormatter())
        ax.yaxis.set_minor_formatter(mticker.NullFormatter())
        # Set plot range
        ax.set_xlim(0.95*min(freqs),1.05*max(freqs))
        ax.legend()
        matplotlib.pyplot.savefig(plotfile+".pdf")
        del ax # cleanup
        # if doFit and sig given, get residual, chi^2
        if doFit and sig:
            sum1=0.0; cnt1=0
            for i in range(0,nch):
                y = FitParms[0]*(freqs[i]/refFreq)**FitParms[1]
                d = vvals[i]-y
                sum1 += (d/ssig[i])**2; cnt1+=1
                # end sums
            chi2 = (sum1/(cnt1-2))**0.5
            print ("reduced chi2 =",chi2)
            return chi2
            # end get RMS residual
        return
    except Exception as exception:
        print(exception)
        print ("Failed or No matplotlib")
# end SpectrumPlot
