#Debug sample script that does a full setup of AIPS crap
# Gollum  user 103 disk 1
from __future__ import absolute_import
from __future__ import print_function
import OErr, OSystem, UV, Obit, Obit, ObitTalkUtil

# On Smeagle
adirs = ["/export/ssd/bcotton/SSD", \
         "/export/raid_1/aips/DATA/SMEAGLE_4", \
         "/export/raid_1/aips/DATA/SMEAGLE_5", \
         "/export/raid_1/aips/DATA/SMEAGLE_6"]
fdirs = ["/export/raid_1/bcotton/fits"]

# Init Obit
err=OErr.OErr()
user = 103
ObitSys=OSystem.OSystem ("debug", 1, user, len(adirs), adirs, len(fdirs), fdirs, True, False, err)
OErr.printErrMsg(err, "Error with Obit startup")

print("TableListDebug")
import AIPS, FITS, AIPSData, FITSData, os
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

# setup environment
AIPS_ROOT    = "/home/AIPS/"
AIPS_VERSION = "31DEC19/"
DA00         = "/lustre/cv/projects/bcotton/zuul05/DA00/"
OBIT_EXEC    = "/export/ssd/bcotton/Git/Obit/trunk/ObitSystem/"
ObitTalkUtil.SetEnviron(AIPS_ROOT=AIPS_ROOT, AIPS_VERSION=AIPS_VERSION, \
                        OBIT_EXEC=OBIT_EXEC, DA00=DA00, ARCH="LNX64", \
                        aipsdirs=adirs, fitsdirs=fdirs)
os.environ['NET0'] = '/dev/null'  # Reduce screaming

# Make sure AIPS Tasks enabled
if 'LD_LIBRARY_PATH' in os.environ:
    os.environ['LD_LIBRARY_PATH']+=':'+os.environ['AIPS_ROOT']+os.environ['AIPS_VERSION']+os.environ['ARCH']+'/LIBR/INTELCMP/'
else:
    os.environ['LD_LIBRARY_PATH'] = os.environ['AIPS_ROOT']+os.environ['AIPS_VERSION']+os.environ['ARCH']+'/LIBR/INTELCMP/'


print("AIPS.AIPS.disks",AIPS.AIPS.disks)
#DAMN print "AIPSData",help(AIPSData.AIPSCat)
#DAMN AMcat(1)

import  Obit, UV, TableList, OErr,  EVLACal
err=OErr.OErr()
from OTObit import  AUcat, getname
#AUcat(3)
uv=UV.newPAUV('in','NGC1068A','UVDaKa',3,1,True,err)
if err.isErr:
    OErr.printErr(err)
uv.Header(err)
TableList.PCheck(uv.TableList,err)
if err.isErr:
    OErr.printErr(err)

print("input TL OK")

# spectral plot
scr=EVLACal.EVLASpecPlot(uv,"3C48",[0.,2.0], 1,err)
EVLACal.EVLAWritePlots(scr,1,0,"./testPlot.ps",err,logfile='doofus.log')
scr.Zap(err)

#OSystem.Shutdown()
