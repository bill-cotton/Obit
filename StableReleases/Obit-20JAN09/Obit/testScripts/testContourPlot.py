# test Contour plot
import Obit, OTObit, Image, ImageUtil, OSystem, OErr, OPlot, math
# Init Obit
err=OErr.OErr()
ObitSys=OSystem.OSystem ("Plot", 1, 100, 0, ["None"], 1, ["../testIt/"], 1, 0, err)
OErr.printErrMsg(err, "Error with Obit startup")

# Contour plot
file = "AGNVLA.fits"
plotfile = "testCont.ps"
lev=0.1
cntfac = 2.0**0.5
blc =[ 230,119]
trc = [285,178]
img = Image.newPFImage("in image",file, 1, True, err)
info=img.List
info.set("BLC",blc)
info.set("TRC",trc)

cont = OPlot.newOPlot("cont", err, output=plotfile+"/ps")
OPlot.PContour(cont, "Some AGN", img, lev, cntfac, err)
# Text
txt=[(265.0,155.0,"A"), \
     (258.0,155.0,"B"), \
     (251.0,152.0,"C"), \
     (247.0,149.0,"D"), \
     (241.0,146.0,"E"), \
     (244.0,141.0,"F")]
OPlot.PSetCharSize(cont, 1 ,err)
for ss in txt:
    x=-(ss[0]-1.-250.9)*2.0/60.0; y=(ss[1]-1.-150.1)*2.0/60.0;
    OPlot.PText(cont, x, y, 180.0, 0.5, ss[2],err)

OPlot.PShow(cont,err)
print "Coutour plot of",file,"written to",plotfile

# Shutdown Obit
OErr.printErr(err)
OSystem.Shutdown(ObitSys)
