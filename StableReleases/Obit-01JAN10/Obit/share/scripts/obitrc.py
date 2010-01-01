# Startup script
print "Executing startup script "
import ObitTalkUtil

###################### Define ###################################
# Define AIPS_ROOT and AIPS_VERSION for access to AIPS Software
AIPS_ROOT    = "/export/data_1/users/aips/"
AIPS_VERSION = "31DEC06/"
DA00         = "/export/data_1/users/aips/DA00/SMEAGLE/"
# Define OBIT_EXEC for access to Obit Software 
OBIT_EXEC    = None  # (def /usr/lib/obit/bin)
OBIT_EXEC    = "/export/data_1/users/bcotton/Software.dir/SVN/ObitInstall/ObitSystem/Obit/"

# Define AIPS directories (URL, disk name)
# URL = None for local disks
aipsdirs = [ \
    (None, "/export/data_1/aips/DATA/SMEAGLE_1"), \
    (None, "/export/data_1/aips/DATA/SMEAGLE_2"), \
    (None, "/export/data_1/aips/DATA/SMEAGLE_3"), \
    (None, "/export/data_2/aips/DATA/SMEAGLE_4"), \
    (None, "/export/data_2/aips/DATA/SMEAGLE_5"), \
    (None, "/export/data_2/aips/DATA/SMEAGLE_6"), \
    (None, "/export/data_2/bcotton/SMEAGLE_7")]

# Define FITS directories (URL, disk name)
# URL = None for local disks
fitsdirs = [ \
    (None, "/export/data_1/users/bcotton/Software.dir/AIPS/FITS")]

# setup environment
ObitTalkUtil.SetEnviron(AIPS_ROOT=AIPS_ROOT, AIPS_VERSION=AIPS_VERSION, \
                        OBIT_EXEC=OBIT_EXEC, ARCH="LINUX", \
                        aipsdirs=aipsdirs, fitsdirs=fitsdirs)

# List directories
ObitTalkUtil.ListAIPSDirs()
ObitTalkUtil.ListFITSDirs()

# Any other customization goes here
