# Determine polarization angle at positions in Q,U,cube, plot
# Vertical dashed line is the reference wavelength^2
import math
from math import isnan
import AIPS, FITS, OPlot, Image, FArray, RMFit

# Setup Obit
import OErr, OSystem

# Cheeta
adirs = [
    "/export/data_1/bcotton/DATA/CHEETA_1/", \
    #"/export/nvme/bcotton/CHEETA_1/", \
    "/export/data_1/bcotton/DATA/CHEETA_2/", \
    "/export/data_1/bcotton/DATA/CHEETA_3/", \
    "/export/data_1/bcotton/DATA/CHEETA_4/", \
]
fdirs = [
    "/export/data_1/bcotton/FITS", \
]
# Smeagle
adirs = [
    "/mnt/ramdisk/RAM1/", \
    "/export/ssd/bcotton/SSD", \
    "/export/raid_1/aips/DATA/SMEAGLE_4", \
    "/export/raid_1/aips/DATA/SMEAGLE_5", \
    "/export/raid_1/aips/DATA/SMEAGLE_6", \
    "/export/raid_2/aips/DATA/SMEAGLE_7", \
    "/export/raid_2/aips/DATA/SMEAGLE_8", \
    "/export/raid_2/aips/DATA/SMEAGLE_9" \
]
fdirs = ["/export/raid_1/bcotton/fits"]


# Init Obit
err=OErr.OErr()
user = 100
ObitSys=OSystem.OSystem ("plotRMFit", 1, user, len(adirs), adirs, \
                         len(fdirs), fdirs, True, False, err)
OErr.printErrMsg(err, "Error with Obit startup")

# Setup AIPS, FITS
AIPS.AIPS.userno = user
disk = 0
for ad in adirs:
    disk += 1
    AIPS.AIPS.disks.append(AIPS.AIPSDisk(None, disk, ad))
disk = 0
for fd in fdirs:
    disk += 1
    FITS.FITS.disks.append(FITS.FITSDisk(None, disk, ad))
OSystem.PAllowThreads(24)  # 24 threads
#
#Perley Butler values, freq (GHz) frac pol, EVPA
PB={"3C48": [(1.050, 0.003, 25.),(1.45, 0.005, 140.), (1.64, 0.007, -5.0)], \
    "3C138":[(1.050, 0.056, -14.),(1.45, 0.075, -11.), (1.64, 0.084, -10.0)], \
    "3C286":[(1.050, 0.086, 33.),(1.45, 0.095, 33.), (1.64, 0.099, 33.)], \
    "J0521+1638_S":[(1.050, 0.056, -14.),(1.45, 0.075, -11.), (1.64, 0.084, -10.0),
                  (1.95,9.0,-10.), (2.45,10.4,-9.), (2.95,10.7,-10.), (3.25,10.0,-10.)], \
    "J1331+3030_S":[(1.050, 0.086, 33.),(1.45, 0.095, 33.), (1.64, 0.099, 33.),
                  (1.95,10.1,33.), (2.45,10.5,33.), (2.95,10.8,-33.), (3.25,10.9,-33.)], \
}
# Do RM plotting
src = [(   '0408-65_S',(459, 459)), \
       ('J0521+1638_S',(456, 457)), \
       ('J0538-4405_S',(465, 465)), \
       ('J1331+3030_S',(464, 464)),]
adisk = 1; aseq=200; iclass='ICals'; qclass='QCals'; uclass='UCals'; post = ''
fblank = FArray.fblank

for s in src:
    field = s[0]
    name = field+post; print ('do',field)
    icube = Image.newPAImage('I', name[0:12],iclass, adisk,aseq,True,err)
    qcube = Image.newPAImage('Q', name[0:12],qclass, adisk,aseq,True,err)
    ucube = Image.newPAImage('U', name[0:12],uclass, adisk,aseq,True,err)
    OErr.printErrMsg(err,message='Finding images')
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
        #print ip,frq
    #
    frq0 = d['FREL0001'][2][0] # ref freq
    lamb20 = (clight/frq0)**2
    print('lambdas',lamb2[0],lamb20)
    pixel = s[1]  # Pixel to plot
    label = field
    icube.Open(Image.READONLY,err)
    qcube.Open(Image.READONLY,err)
    ucube.Open(Image.READONLY,err)
    iflux=[]; iRMS=[]  # Ipol flux, RMS in plane
    qflux=[]; qRMS=[]  # Qpol flux, RMS in plane
    uflux=[]; uRMS=[]  # Upol flux, RMS in plane
    llamb2 = []; nplot = 0
    plane = [1,1,1,1,1]
    for ip in range(nterm+1,nterm+nspec+1):
        plane[0] = ip
        icube.GetPlane(None, plane, err)
        if icube.FArray.get(pixel[0]-1, pixel[1]-1)!=0.0:
            llamb2.append(lamb2[ip-nterm-1])
            nplot += 1  # How many valid data points
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
    OErr.printErrMsg(err,message='Determining Rotation Measure')

    # Get quanties to plot and errors
    phs=[];  ephs=[]; l2=[]; nsamp = 0
    fpol=[]; efpol=[]; qqflux=[]; uuflux=[]; qqRMS=[]; uuRMS=[]
    for ip in range(0,nplot):
        # Data OK?
        #FITS if (not isnan(iflux[ip])) and (not isnan(qflux[ip])) and (not isnan(uflux[ip])):
        if (not (iflux[ip]==fblank) and (not (qflux[ip]==fblank)) and (not (uflux[ip]==fblank)) \
            and (not (iflux[ip]==0.0)) and (not (qflux[ip]==0.0))):
            nsamp += 1
            l2.append(lamb2[ip])
            phs.append(0.5*math.degrees(math.atan2(uflux[ip],qflux[ip])))
            arg    = abs(uflux[ip]/qflux[ip])
            sigarg = arg*((qRMS[ip]**2/qflux[ip]**2)+(uRMS[ip]**2/uflux[ip]**2))**0.5
            ephs.append(abs(math.degrees(sigarg/(1+arg**2))))
            arg = ((uflux[ip]**2+qflux[ip]**2)**0.5)/iflux[ip]
            fpol.append(100.*arg)
            vararg  = (qRMS[ip]**2 +uRMS[ip]**2)
            sigarg  = arg*((vararg/arg**2) + iRMS[ip]**2/iflux[ip]**2)**0.5
            efpol.append(100.*sigarg)
            ip, iflux[ip], (uflux[ip]**2+qflux[ip]**2)**0.5, 0.5*math.degrees(math.atan2(uflux[ip],qflux[ip]))
            qqflux.append(qflux[ip]); qqRMS.append(qRMS[ip])
            uuflux.append(uflux[ip]); uuRMS.append(uRMS[ip])
    
    # Remove phase wraps
    for ip in range(1,nsamp):
        if phs[ip]-phs[ip-1]>130.:
            phs[ip] -=180.
        if phs[ip]-phs[ip-1]<-130.:
            phs[ip] +=180.
        
    # Fit RM
    refLamb2 = lamb2[nspec//2]
    #RMParms = RMFit.PSingle(refLamb2, llamb2, qflux, qRMS, uflux, uRMS, err)
    RMParms = RMFit.PSingle(refLamb2, l2, qqflux, qqRMS, uuflux, uuRMS, err)
    evpa0 =  math.degrees(RMParms[1]+RMParms[0]*(lamb20-refLamb2)) # EVPA @ ref. freq
    eevpa0 = math.degrees(RMParms[3])
    print (field, math.degrees(RMParms[1]), lamb20, evpa0, eevpa0)
    print ('phs',phs)
    # Guess offset (180 ambiguity) from difference in evpa and phs[0]
    pdif =  phs[0]-math.degrees(RMParms[1])
    evpa0 += round(pdif/180)*180
    print (math.degrees(RMParms[1]),phs[0],pdif)

    pname = field+'_'+str(aseq)+'_RMPlot'
    plot=OPlot.newOPlot("Plot", err,output=pname+".ps/ps",ny=2)
    # Plot EVPA
    plot.List.set("TITLE"," %s: RM=%6.2f (%5.2f), EVPA#d_ref#u =%7.2f (%5.2f)"\
                  %(label, RMParms[0], RMParms[2],evpa0,eevpa0))
    plot.List.set("YLABEL","EVPA (deg)")
    plot.List.set("XLABEL","#gl#u2#d (m#u2#d)")
    plot.List.set("XMAX",1.05*lamb20)
    plot.List.set("XMIN",0.95*lamb2[nspec-1])
    avgP = sum(phs)/len(phs)
    #print ("avgP",avgP)
    plot.List.set("YMIN",(min(min(phs)-15,avgP-20.)))
    plot.List.set("YMAX",(max(max(phs)+15,avgP+20.)))
    plot.List.set("CSIZE",1)
    plot.List.set("SSIZE",3)
    plot.List.set("LWIDTH",3)
    OPlot.PXYErr(plot, 2, l2, phs, ephs, err)
    # In Perley&Butler?
    if field in PB:
        for dd in PB[field]:
            x1 = ((clight/(dd[0]*1.0e9))**2)
            y1 = dd[2];
            #print "PB", dd,x1,y1
            OPlot.PDrawSymbol(plot,x1,y1,8,err)
    # Plot RM fit - a few times
    noff = 7
    offs = [0.0, 180, -180., 360., -360, 540. -540.]
    for off in offs:
        x1 = 1.1*lamb2[0] - refLamb2
        y1 = math.degrees(RMParms[1]+RMParms[0]*x1)+off
        x2 = 0.9*lamb2[nspec-1] - refLamb2
        y2 = math.degrees(RMParms[1]+RMParms[0]*x2)+off
        OPlot.PDrawLine(plot, 1.1*lamb2[0], y1, 0.9*lamb2[nspec-1], y2, err)
    OPlot.PSetLineStyle(plot, 3, err)
    OPlot.PDrawLine(plot, lamb20, -2000., lamb20, 2000., err)
    OPlot.PSetLineStyle(plot, 1, err)
    # Plot fractional poln
    plot.List.set("TITLE"," ")
    plot.List.set("XLABEL","#gl#u2#d (m#u2#d)")
    plot.List.set("YLABEL","frac. pol. (%)")
    plot.List.set("YMIN",0.0)
    plot.List.set("YMAX",min(1.20*max(fpol),20.0))
    OPlot.PXYErr(plot, 2, lamb2, fpol, efpol, err)
    # In Perley&Butler?
    if field in PB:
        for dd in PB[field]:
            x1 = ((clight/(dd[0]*1.0e9))**2)
            y1 = 100.*dd[1];
            OPlot.PDrawSymbol(plot,x1,y1,8,err)
    OPlot.PShow(plot,err)
    OErr.printErrMsg(err,message='Plotting Rotation Measure')
# end loop over images

# Shutdown Obit
OErr.printErr(err)
#OSystem.Shutdown(ObitSys)

