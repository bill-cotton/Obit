#Debug sample script that does a full setup of AIPS crap
# Gollum  user 100 disk 1
from __future__ import absolute_import
from __future__ import print_function
import OErr, OSystem, UV, Obit, Obit, ObitTalkUtil

# On Smeagle
adirs = ["/export/ssd/bcotton/SSD", \
         "/export/raid_1/aips/DATA/SMEAGLE_4"]
fdirs = ["/export/raid_1/bcotton/fits"]

# Init Obit
err=OErr.OErr()
user = 100
ObitSys=OSystem.OSystem ("debug", 1, user, len(adirs), adirs, len(fdirs), fdirs, True, False, err)
OErr.printErrMsg(err, "Error with Obit startup")

print("NewDebug")
import AIPS, FITS, AIPSData, FITSData
from OTObit import newDisplay, tvlod, AMcat
# This shit really bites
AIPS.AIPS.userno = user
disk = 0
for ad in adirs:
    AIPS.AIPS.disks.append(AIPS.AIPSDisk(None, disk, ad))
    disk += 1

disk = 0
for fd in fdirs:
    FITS.FITS.disks.append(FITS.FITSDisk(None, disk, ad))
    disk += 1

# List directories
#ObitTalkUtil.ListAIPSDirs()
#ObitTalkUtil.ListFITSDirs()

print("AIPS.AIPS.disks",AIPS.AIPS.disks)
#DAMN print "AIPSData",help(AIPSData.AIPSCat)
#DAMN AMcat(1)

import Image,FitModel,FitRegion, Obit, ODisplay, OErr
err=OErr.OErr()
from OTObit import newDisplay, tvlod, AMcat, getFITS
AMcat(1)
##disp=ODisplay.ODisplay(8765)
#url = "http://localhost:8765/RPC2" 
#disp = ODisplay.ODisplay("ObitView", url, err) 
#tvlod(x)
x=Image.newPAImage('cc','DemoTest','I',1,1,True,err)
if err.isErr:
    printErr(err)
#x.TVFit(disp,err)
x.Header(err)
#import FitModel
#m=FitModel.FitModel(mtype=1, Peak=0.00054, DeltaX=14, DeltaY=17, parms=[5.0, 5.0, 0.0])
