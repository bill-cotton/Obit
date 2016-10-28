# Determine polarization angle at positions in Q,U,cube
import math
from math import isnan
import OPlot, Image, RMFit

# Setup Obit
import OErr, OSystem

# On Smeagle
adirs = [
    "/export/sdd/bcotton/SDD", \
    "/export/data_2/bcotton/DATA/SMEAGLE_4/" \
]
fdirs = ["/export/sdd/bcotton/FITS"]


# Init Obit
err=OErr.OErr()
user = 100
ObitSys=OSystem.OSystem ("plotRMFit", 1, user, len(adirs), adirs, \
                         len(fdirs), fdirs, True, False, err)
OErr.printErrMsg(err, "Error with Obit startup")

# Setup AIPS, FITS
#AIPS.AIPS.userno = user
#disk = 0
#for ad in adirs:
#    disk += 1
#    AIPS.AIPS.disks.append(AIPS.AIPSDisk(None, disk, ad))
#disk = 0
#for fd in fdirs:
#    disk += 1
#    FITS.FITS.disks.append(FITS.FITSDisk(None, disk, ad))
#


#icube = Image.newPFImage("I","HercA.IAll_L.fits",0,True,err)
#qcube = Image.newPFImage("I","HercA.QAll_L.fits",0,True,err)
#ucube = Image.newPFImage("I","HercA.UAll_L.fits",0,True,err)

icube = Image.newPFImage("I","Results/J2002+006I_Mosaic.fits",0,True,err)
qcube = Image.newPFImage("Q","Results/J2002+006Q_Mosaic.fits",0,True,err)
ucube = Image.newPFImage("U","Results/J2002+006U_Mosaic.fits",0,True,err)

# get lambda sq
d=icube.Desc.List.Dict
clight=2.997924562e8
nterm = d['NTERM'][2][0]
nspec = d['NSPEC'][2][0]
lamb2 = []
for ip in range(nterm+1,nterm+nspec+1):
    key = "FREQ%4.4d"%(ip-nterm)
    frq = d[key][2][0]
    lamb2.append((clight/frq)**2)
    print ip,frq


# Loop over positions
# Initial

pixels=[(3977,3034)] # J2002+006

for pixel in pixels:
    icube.Open(Image.READONLY,err)
    qcube.Open(Image.READONLY,err)
    ucube.Open(Image.READONLY,err)
    iflux=[]; iRMS=[]  # Ipol flux, RMS in plane
    qflux=[]; qRMS=[]  # Qpol flux, RMS in plane
    uflux=[]; uRMS=[]  # Upol flux, RMS in plane
    plane = [1,1,1,1,1]
    for ip in range(nterm+1,nterm+nspec+1):
        plane[0] = ip
        icube.GetPlane(None, plane, err)
        iRMS.append(icube.FArray.RMS)
        iflux.append(icube.FArray.get(pixel[0]-1, pixel[1]-1))
        qcube.GetPlane(None, plane, err)
        qRMS.append(qcube.FArray.RMS)
        qflux.append(qcube.FArray.get(pixel[0]-1, pixel[1]-1))
        ucube.GetPlane(None, plane, err)
        uRMS.append(ucube.FArray.RMS)
        uflux.append(ucube.FArray.get(pixel[0]-1, pixel[1]-1))
    
    icube.Close(err)
    qcube.Close(err)
    ucube.Close(err)

    # Get quanties to plot and errors
    phs=[];  ephs=[]; l2=[]; nsamp = 0
    fpol=[]; efpol=[]
    for ip in range(0,nspec):
        # Data OK?
        if (not isnan(qflux[ip])) and (not isnan(qflux[ip])) and (not isnan(uflux[ip])):
            nsamp += 1
            l2.append(lamb2[ip])
            phs.append(0.5*math.degrees(math.atan2(uflux[ip],qflux[ip])))
            arg    = abs(uflux[ip]/qflux[ip])
            sigarg = arg*((qRMS[ip]**2/qflux[ip]**2)+(uRMS[ip]**2/uflux[ip]**2))**0.5
            ephs.append(abs(math.degrees(sigarg/(1+arg**2))))
            arg = ((uflux[ip]**2+qflux[ip]**2)**0.5)/iflux[ip]
            fpol.append(arg)
            vararg  = (qRMS[ip]**2 +uRMS[ip]**2)
            sigarg  = arg*((vararg/arg**2) + iRMS[ip]**2/iflux[ip]**2)**0.5
            efpol.append(sigarg)
            ip, iflux[ip], (uflux[ip]**2+qflux[ip]**2)**0.5, 0.5*math.degrees(math.atan2(uflux[ip],qflux[ip]))
    
    # Remove phase wraps
    for ip in range(1,nsamp):
        if phs[ip]-phs[ip-1]>130.:
            phs[ip] -=180.
        if phs[ip]-phs[ip-1]<-130.:
            phs[ip] +=180.
    
    # Fit RM
    refLamb2 = lamb2[nspec/2]
    RMParms = RMFit.PSingle(refLamb2, lamb2, qflux, qRMS, uflux, uRMS, err)

    pname = "Poln%4.4d.%4.4d"%(pixel[0],pixel[1])
    print "Plot", pname
    plot=OPlot.newOPlot("Plot", err,output=pname+".ps/ps",ny=2)
    # Plot EVPA
    plot.List.set("TITLE","Pixel %d %d, RM=%6.2f (%5.2f)"\
                  %(pixel[0],pixel[1], RMParms[0], RMParms[2]))
    plot.List.set("YLABEL","EVPA (deg)")
    plot.List.set("CSIZE",1)
    plot.List.set("SSIZE",3)
    plot.List.set("LWIDTH",3)
    OPlot.PXYErr(plot, 2, l2, phs, ephs, err)
    # Plot RM fit - a few times
    noff = 7
    offs = [0.0, 180, -180., 360., -360, 540. -540.]
    for off in offs:
        x1 = lamb2[0] - refLamb2
        y1 = math.degrees(RMParms[1]+RMParms[0]*x1)+off
        x2 = lamb2[nspec-1] - refLamb2
        y2 = math.degrees(RMParms[1]+RMParms[0]*x2)+off
        OPlot.PDrawLine(plot, lamb2[0], y1, lamb2[nspec-1], y2, err)
    # Plot fractional poln
    plot.List.set("TITLE"," ")
    plot.List.set("XLABEL","#gl#u2#d (m#u2#d)")
    plot.List.set("YLABEL","frac. pol.")
    plot.List.set("XMAX",lamb2[0])
    plot.List.set("XMIN",lamb2[nspec-1])
    plot.List.set("YMIN",0.0)
    plot.List.set("YMAX",0.6)
    OPlot.PXYErr(plot, 2, lamb2, fpol, efpol, err)
    PLPLOTOPlot.PShow(plot,err)
# end loop over pixels

# Shutdown Obit
OErr.printErr(err)
OSystem.Shutdown(ObitSys)
