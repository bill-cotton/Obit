# Define AIPS  and FITS disks
# On Smeagle
adirs = [(None, "/export/data_1/aips/DATA/SMEAGLE_1"),
         (None, "/export/data_1/aips/DATA/SMEAGLE_2"), \
         (None, "/export/data_1/aips/DATA/SMEAGLE_3"), \
         (None, "/export/data_2/aips/DATA/SMEAGLE_4"), \
         (None, "/export/data_2/aips/DATA/SMEAGLE_5"), \
         (None, "/export/data_2/aips/DATA/SMEAGLE_6")]
fdirs = [(None, "/export/data_1/bcotton/Software.dir/AIPS/FITS")]

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

# On Mortibus
adirs = [(None, "/export/data_1/MORTIBUS_1"),
         (None, "/export/data_1/MORTIBUS_2"),
         (None, "/export/data_1/MORTIBUS_3"),
         (None, "/export/data_1/MORTIBUS_4"),
         (None, "/export/data_2/obit/MORTIBUS_5"),
         (None, "/export/data_2/obit/MORTIBUS_6"),
         (None, "/export/data_2/obit/MORTIBUS_7"),
         (None, "/export/data_2/obit/MORTIBUS_8"),
         (None, "/export/data_3/obit/MORTIBUS_9"),
         (None, "/export/data_3/obit/MORTIBUS_10"),
         (None, "/export/data_3/obit/MORTIBUS_11"),
         (None, "/export/data_3/obit/MORTIBUS_12")]
fdirs = [(None, "/export/data_1/FITS"),
         (None, "/export/data_2/FITS"),
         (None, "/export/data_3/FITS")]


# On Panther
adirs = [(None, "/export/data/aips/DATA/PANTHER_1")]
fdirs = [(None, "/export/users/bcotton/Software.dir/AIPS/FITS")]

############################# Initialize OBIT ##########################################
err     = OErr.OErr()
user    = 105
ObitSys = OSystem.OSystem ("Pipeline", 1, user, 0, [" "], \
                         0, [" "], True, False, err)
OErr.printErrMsg(err, "Error with Obit startup")

# Setup AIPS, FITS
AIPS.userno = user

AIPS_ROOT    = "/export/data/aips/"
AIPS_VERSION = "31DEC09/"
DA00         = "/export/data/aips/DA00/PANTHER"
# Define OBIT_EXEC for access to Obit Software 
OBIT_EXEC    = "/export/data_1/obit/ObitInstall/ObitSystem/Obit/"

# setup environment
ObitTalkUtil.SetEnviron(AIPS_ROOT=AIPS_ROOT, AIPS_VERSION=AIPS_VERSION, \
                        OBIT_EXEC=OBIT_EXEC, DA00=DA00, ARCH="LINUX", \
                        aipsdirs=adirs, fitsdirs=fdirs)

# List directories
ObitTalkUtil.ListAIPSDirs()
ObitTalkUtil.ListFITSDirs()

# Disks to avoid
noScrat     = [0]          # AIPS disks to avoid 
nThreads    = 2            # Number of threads allowed
disk        = 1            # AIPS disk number
