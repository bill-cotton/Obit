# Define AIPS  and FITS disks
# On Gollum
adirs = [(None, "/export/data_1/GOLLUM_1"),
         (None, "/export/data_1/GOLLUM_2"),
         (None, "/export/data_1/GOLLUM_3"),
         (None, "/export/data_1/GOLLUM_4"),
         (None, "/export/data_2/GOLLUM_5"),
         (None, "/export/data_2/GOLLUM_6"),
         (None, "/export/data_2/GOLLUM_7"),
         (None, "/export/data_2/GOLLUM_8")]
fdirs = [(None, "/export/users/aips/FITS")]

# On Cheeta
adirs = [ \
          (None, "/export/nvme/bcotton/CHEETA_1/"), \
          (None, "/export/data_1/bcotton/DATA/CHEETA_2/"), \
          (None, "/export/data_1/bcotton/DATA/CHEETA_3/"), \
          (None, "/export/data_1/bcotton/DATA/CHEETA_4/"), \
      ]
fdirs = [(None, "/export/data_1/bcotton/FITS")]
AIPS_ROOT    = "/home/aips/"
AIPS_VERSION = "31DEC19/"
DA00         = "/lustre/cv/projects/bcotton/zuul05/DA00/"
# Define OBIT_EXEC for access to Obit Software 
OBIT_EXEC    = None  # (def /usr/lib/obit/bin)
OBIT_EXEC    = "/export/ssd/bcotton/Git/Obit/trunk/ObitSystem/"


# On Smeagle
adirs = [(None, "/export/raid_1/aips/DATA/SMEAGLE_4"),
         (None, "/export/raid_1/aips/DATA/SMEAGLE_5"), \
         (None, "/export/raid_1/aips/DATA/SMEAGLE_6"), \
         (None, "/export/raid_2/aips/DATA/SMEAGLE_7"), \
         (None, "/export/raid_2/aips/DATA/SMEAGLE_8"), \
         (None, "/export/raid_2/aips/DATA/SMEAGLE_8")]
fdirs = [(None, "/export/raid_1/bcotton/fits")]

############################# Initialize OBIT ##########################################
err     = OErr.OErr()
user    = 103
ffdirs = []
for f in fdirs:
    ffdirs.append(f[1])
aadirs = []
for a in adirs:
    aadirs.append(a[1])
ObitSys = OSystem.OSystem ("Pipeline", 1, user, len(aadirs), aadirs, \
                         len(ffdirs), ffdirs, True, False, err)
OErr.printErrMsg(err, "Error with Obit startup")

# Setup AIPS, FITS
AIPS.userno = user

AIPS_ROOT    = "/home/AIPS/"
AIPS_VERSION = "31DEC19/"
DA00         = "/lustre/cv/projects/bcotton/zuul05/DA00/"
# Define OBIT_EXEC for access to Obit Software 
OBIT_EXEC    = None  # (def /usr/lib/obit/bin)
OBIT_EXEC    = "/export/ssd/bcotton/Git/Obit/trunk/ObitSystem/"

# setup environment
ObitTalkUtil.SetEnviron(AIPS_ROOT=AIPS_ROOT, AIPS_VERSION=AIPS_VERSION, \
                        OBIT_EXEC=OBIT_EXEC, DA00=DA00, ARCH="LNX64", \
                        aipsdirs=adirs, fitsdirs=fdirs)

# Make sure AIPS Tasks enabled
if 'LD_LIBRARY_PATH' in os.environ:
    os.environ['LD_LIBRARY_PATH'] = os.environ['AIPS_ROOT']+os.environ['AIPS_VERSION']+os.environ['ARCH']+'/LIBR/INTELCMP/:'+os.environ['LD_LIBRARY_PATH']
else:
    os.environ['LD_LIBRARY_PATH'] = os.environ['AIPS_ROOT']+os.environ['AIPS_VERSION']+os.environ['ARCH']+'/LIBR/INTELCMP/'
# List directories
#debug ObitTalkUtil.ListAIPSDirs()
#debug ObitTalkUtil.ListFITSDirs()

# Disks to avoid
noScrat     = [0,0,0]      # AIPS disks to avoid 
nThreads    = 16           # Number of threads allowed
disk        =  3           # AIPS disk number
