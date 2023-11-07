# python/Obit script to plot Faraday rotation measure in a time series of visibilities.
# Fitting uses per IF values of Q and U visibility amplitudes.
# Data should be AIPS FITAB format dataaveraged over baseline (task AvgBL) 
# and over channels in an # IF (task Split, NB doesn't work well on fits UV data)
# Arguments:
# 1) Name of uv data in AIPS/UVTAB format (e.g. (myUVData.uvtab")
# 2) Name of output plot file (e.g. "myRMPlot.ps")


import sys, Obit, UV, Table, RMFit, OErr, OPlot, OSystem
err=OErr.OErr()
ObitSys=OSystem.OSystem ("FitFarRot", 1, 100, 1, ["../AIPSdata/"], 1, ["../PythonData/"], 1, 0, err)
OErr.printErrMsg(err, "Error with Obit startup")

# Get args
UV_File = sys.argv[1]
pfile = sys.argv[2]

uv=UV.newPFUV('uv',UV_File, 0,True,err,nvis=1)
OErr.printErrMsg(err, "Error finding input data")
uv.List.set("doCalSelect",True)
uv.List.set("Stokes","IQUV")
uv.List.set("corrType",0)
z=uv.Open(UV.READCAL, err)
OErr.printErrMsg(err, "Error opening input data")
nvis = uv.Desc.Dict['nvis']
jif   = uv.Desc.Dict["jlocif"]        # Order in data of IF
jf    = uv.Desc.Dict['jlocf']         # Order in data of freq.
nif   = uv.Desc.Dict["inaxes"][jif]   # Number of IFs
js    = uv.Desc.Dict["jlocs"]         # Order in data of Stokes
nstok = uv.Desc.Dict["inaxes"][js]    # Number of Stokes
ifinc = uv.Desc.Dict["incif"]//3   
# IF Frequencies
clight = 2.997924562e8   # Speed of light
fqtab=uv.NewTable(Table.READONLY, "AIPS FQ",1,err)
fqtab.Open(Table.READONLY, err)
fqrow = fqtab.ReadRow(1,err)
OErr.printErrMsg(err, "Error getting IF frequencies data")
freqs = fqrow['IF FREQ']; lamb2 = []; refLamb2=1.0e-6
for i in range(0,len(freqs)):
    freqs[i] += uv.Desc.Dict['crval'][jf]
    lamb2.append(((clight/freqs[i])**2))  # Wavelengths^2

fqtab.Close(err)

x=[]; y=[]; e=[]
for ivis in range(0,nvis):
    try:
        v = uv.ReadVis(err, firstVis=ivis+1)
    except Exception as exception:
        print (exception)
    else:
        pass
    OErr.printErrMsg(err, "Error reading data")
    t = v.time*24 # Hours
    qrms = []; urms = [];  qvis = []; uvis = []; llamb2=[]
    for iif in range(0,nif):
        if v.vis[iif*ifinc+1][1]>0.0:
            qvis.append(v.vis[iif*ifinc+1][0].real)
            qrms.append(0.05); llamb2.append(lamb2[iif])
        if v.vis[iif*ifinc+2][1]>0.0:
            uvis.append(v.vis[iif*ifinc+2][0].real)
            urms.append(0.05)
    # Fit RM
    rmval = RMFit.PSingle(refLamb2, llamb2,qvis, qrms, uvis, urms, err)
    x.append(t); y.append(rmval[0]); e.append(min (0.01,rmval[2]))
    
    # end loop 
    
OErr.printErrMsg(err, "Error reading or fitting data")
plot=OPlot.newOPlot("plot", err, output=pfile+'/ps')
plot.List.set("XLABEL","Time (hrs)")
plot.List.set("YLABEL", "Rotation Measure (rad m#u-2#d)")
plot.List.set("TITLE",UV_File)
#Huh? OPlot.PXYErr(plot, 3, x, y, e, err)
OPlot.PXYPlot(plot, 3, x, y, err)
OPlot.PShow(plot, err)
OErr.printErrMsg(err, "Error plotting data")

# Shutdown Obit
OErr.printErr(err)
del ObitSys
