# Plot Spectrum at a list of positions in a FITS ImageMF
# On either raw (no PBCor) or PBCorImageMF.py PB corrected

name   = 'Abell3395' # Base of plot name
inFile = 'Abell3395_IPol.fits'; fdisk=0  # Input file on cwd
#   Specify positions (name, ra, dec)
srcpos = [ \
           ('Source_1', '06:26:27.35', '-54:32:51.3'), \
           ('Source_2', '06:26:36.67', '-54:38:49.3'), \
           ('Source_3', '06:26:45.92', '-54:40:25.1'), \
           ('Source_4', '06:26:19.93', '-54:40:18.4'), \
           ('Source_5', '06:26:04.71', '-54:40:14.2'), \
]

# Shouldn't need to change below here
import math
from math import isnan
import FITS, OPlot, Image, ImageDesc, ImageMF, FArray, FInterpolate, ImageUtil, SpectrumFit

# Setup Obit
import OErr, OSystem

adirs = []
fdirs = ['./']

# Init Obit
err=OErr.OErr()
user = 100
ObitSys=OSystem.OSystem ("plotSpectrum", 1, user, len(adirs), adirs, \
                         len(fdirs), fdirs, True, False, err)
OErr.printErrMsg(err, "Error with Obit startup")

fblank = FArray.fblank

# In image
icube = Image.newPFImage('In', inFile, fdisk, True, err)
# Image info
nterm = icube.Desc.List.Dict['NTERM'][2][0]
nspec = icube.Desc.List.Dict['NSPEC'][2][0]
freqs = []
# Interpolator
fi = FInterpolate.FInterpolate('FI', icube.FArray, icube.Desc, 2)
# Get pixel locations
pixels = []; vals = []
for s in srcpos:
    pos = [ImageDesc.PHMS2RA(s[1]), ImageDesc.PDMS2Dec(s[2])]
    pixel = ImageDesc.PGetPixel(icube.Desc, pos, err)
    pixels.append(pixel)
    vals.append([])

# end loop
OErr.printErrMsg(err,message='Finding positions')
# Loop over planes interpolating values
rms = []
for i in range(0,nspec):
    # Read plane
    plane=[i+nterm+1,1,1,1,1]
    Image.PGetPlane (icube, None, plane, err)
    rmss = icube.FArray.RMS
    if not isnan(rmss): # Check validity
        key = 'FREQ%4.4d'%(i+1)
        freqs.append(1.0e-9*icube.Desc.List.Dict[key][2][0]) # GHz
        rms.append(rmss) # Plane RMS
        # Interpolate positions
        for j in range(0,len(pixels)):
            pixel = pixels[j]
            val = FInterpolate.PPixel(fi, pixel, err)
            #print (i,j,val,rms[i],freqs[i])
            if val>0.0:
                vals[j].append(val)
            else:
                vals[j].append(1.0e-6) # fooey

# end loop over planes
OErr.printErrMsg(err,message='Interpolating values')
    
# Plot
pname = name
plot=OPlot.newOPlot("Plot", err,output=pname+".ps/ps")
i = -1
for s in srcpos:
    i += 1
    # Fit Spectrum
    nterm=2; refFreq = 1.0 # 1 GHz
    FitParms = SpectrumFit.PSingle(nterm,refFreq, freqs, vals[i], rms ,err)
    if FitParms[1]==fblank:  # fit failed?
        FitParms[1] = 0.0
    #print ("\nFit",FitParms)
    #print ("vals",vals[i])
    plot.List.set("TITLE"," %s %s %s %s, #ga=%6.3f (%5.3f)"\
                  %(name, s[0], s[1],s[2], FitParms[1], FitParms[3]))
    plot.List.set("YLABEL","Flux density (Jy)")
    plot.List.set("XLABEL","Frequency (GHz)")
    plot.List.set("XOPT","LBCN")
    plot.List.set("YOPT","LBCN")
    ymn = 0.9*10**int(-0.9+math.log10(min(vals[i])))
    plot.List.set("YMIN",ymn)
    plot.List.set("CSIZE",1)
    plot.List.set("SSIZE",3)
    plot.List.set("LWIDTH",3)
    OPlot.PXYErr(plot, 2, freqs, vals[i], rms, err)
    # Plot fitted spectrum
    x0 = min(freqs); y0 = FitParms[0]*(x0/refFreq)**FitParms[1]
    x1 = max(freqs); y1 = FitParms[0]*(x1/refFreq)**FitParms[1]
    OPlot.PDrawLine(plot, x0, y0, x1, y1, err)
    OErr.printErrMsg(err,message='Plotting Spectrum')
# end loop over positions
OPlot.PShow(plot,err)
OErr.printErrMsg(err,message='Plotting Spectrum')

# Shutdown Obit
OErr.printErr(err)
#OSystem.Shutdown(ObitSys)
