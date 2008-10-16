Debug sample script that does a full setup of AIPS crap
# Gollum  user 101, disk 3
import OErr, OSystem, UV, AIPS, FITS
import VLACal

# On Gollum
adirs = ["/export/data_1/GOLLUM_1",
         "/export/data_1/GOLLUM_2", \
         "/export/data_1/GOLLUM_3", \
         "/export/data_1/GOLLUM_4", \
         "/export/data_2/GOLLUM_5", \
         "/export/data_2/GOLLUM_6", \
         "/export/data_2/GOLLUM_7", \
         "/export/data_2/GOLLUM_8"]
fdirs = ["/export/users/aips/FITS"]

# Init Obit
err=OErr.OErr()
user = 101
ObitSys=OSystem.OSystem ("debug", 1, user, len(adirs), adirs, len(fdirs), fdirs, True, False, err)
OErr.printErrMsg(err, "Error with Obit startup")

# This shit really bites
AIPS.AIPS.userno = user
disk = 0
for ad in adirs:
    disk += 1
    AIPS.AIPS.disks.append(AIPS.AIPSDisk(None, disk, ad))
disk = 0
for fd in fdirs:
    disk += 1
    FITS.FITS.disks.append(FITS.FITSDisk(None, disk, ad))

# Set uv data 
u3=UV.newPAUV("in","20051116", "K BAND", 3,1,True,err)
u3.Header(err)
VLACal.VLAClearCal(u3,err)
calModel="3C286_K.MODEL"
target="M87"
ACal="1331+305"
PCal="1239+075"

VLACal.VLACal(u3, target, ACal, err, calModel=calModel,calDisk=1)
VLACal.VLASplit(u3, target, err, outClass="IKDarr")


