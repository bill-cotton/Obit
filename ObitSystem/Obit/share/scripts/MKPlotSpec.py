# Plot Spectrum at a list of positions in a FITS ImageMF
# On either raw (no PBCor) or PBCorImageMF.py PB corrected

name   = 'Abell_194' # Base of plot name
inFile = 'Abell_194.fits'; fdisk=0  # Input file on cwd
#   Specify positions (name, ra, dec)
srcpos = [ \
           ('core',  '01:26:00.597', '-01:20:43.71'), \
           ('src_1', '01:25:38.632', '-01:11:43.34'), \
           ('src_2', '01:25:13.0137','-01:28:02.41'), \
           ('src_5', '01:24:00.795', '-01:05:04.99'), \
           ('src_8', '01:27:02.936' ,'-02:08:19.02'), \
           ('src_9', '01:23:10.24',  '-01:56:24.7'),  \
           ('src_10','01:23:21.18',  '-01:58:36.75'), \
           ('src_11','01:23:27.71',  '-02:04:45.92'), \
]

# Shouldn't need to change below here
import math
from math import isnan
import FITS, OPlot, Image, ImageDesc, ImageMF, FArray, FInterpolate, ImageUtil

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
for i in range(1,nspec+1):
    key = 'FREQ%4.4d'%i
    freqs.append(1.0e-9*icube.Desc.List.Dict[key][2][0]) # GHz
    
# end loop
# Interpolator
fi = FInterpolate.FInterpolate('FI', icube.FArray, icube.Desc, 2)
# Get pixel locations
pixels = []; vals = []
for s in srcpos:
    pos = [ImageDesc.PHMS2RA(s[1]), ImageDesc.PDMS2Dec(s[2])]
    pixel = ImageDesc.PGetPixel(icube.Desc, pos, err)
    pixels.append(pixel)
    vals.append(nspec*[fblank])

# end loop
OErr.printErrMsg(err,message='Finding positions')
# Loop over planes interpolating values
rms = []
for i in range(0,nspec):
    # Read plane
    plane=[i+nterm+1,1,1,1,1]
    Image.PGetPlane (icube, None, plane, err)
    rms.append(icube.FArray.RMS) # Plane RMS
    # Interpolate positions
    for j in range(0,len(pixels)):
        pixel = pixels[j]
        val = FInterpolate.PPixel(fi, pixel, err)
        print (i,j,val,rms[i],freqs[i])
        if val>0.0:
            vals[j][i] = val

# end loop over planes
OErr.printErrMsg(err,message='Interpolating values')
    
# Plot
pname = name
plot=OPlot.newOPlot("Plot", err,output=pname+".ps/ps")
i = -1
for s in srcpos:
    i += 1
    plot.List.set("TITLE"," %s %s %s %s"\
                  %(name, s[0], s[1],s[2]))
    plot.List.set("YLABEL","Flux density (Jy)")
    plot.List.set("XLABEL","Frequency (GHz)")
    plot.List.set("XOPT","LBCN")
    plot.List.set("YOPT","LBCN")
    #plot.List.set("XMAX",1.1*max(freqs))
    #plot.List.set("XMIN",0.9*min(freqs))
    ymn = 0.9*10**int(-0.9+math.log10(min(vals[i])))
    #ymx = ymn * 10.
    print (ymn,min(vals[i]))
    plot.List.set("YMIN",ymn)
    #plot.List.set("YMAX",ymx)
    plot.List.set("CSIZE",1)
    plot.List.set("SSIZE",3)
    plot.List.set("LWIDTH",3)
    OPlot.PXYErr(plot, 2, freqs, vals[i], rms, err)
    OErr.printErrMsg(err,message='Plotting Spectrum')
# end loop over positions
OPlot.PShow(plot,err)
OErr.printErrMsg(err,message='Plotting Spectrum')

# Shutdown Obit
OErr.printErr(err)
#OSystem.Shutdown(ObitSys)
