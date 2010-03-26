# Sample Startup script
print "Executing startup script "
import ObitTalkUtil

###################### Define ###################################
# Define AIPS_ROOT and AIPS_VERSION for access to AIPS Software
AIPS_ROOT    = "/export/users/aips"
AIPS_VERSION = "31DEC07"
# Define OBIT_EXEC for access to Obit Software 
OBIT_EXEC    = None  # (def /usr/lib/obit/bin)
# This is the lib/obit directory in the installation
OBIT_EXEC    = "/export/data_1/users/bcotton/ObitBinTest/obit-1.1.166-4/lib/obit"

# Define AIPS directories (URL, disk name)
# URL = None for local disks
aipsdirs = [ \
    (None, "/export/data_1/GOLLUM_1"), \
    (None, "/export/data_1/GOLLUM_2"), \
    (None, "/export/data_1/GOLLUM_3"), \
    (None, "/export/data_1/GOLLUM_4"), \
    (None, "/export/data_2/GOLLUM_5"), \
    (None, "/export/data_2/GOLLUM_6"), \
    (None, "/export/data_2/GOLLUM_7"), \
    (None, "/export/data_2/GOLLUM_8")]

# Define FITS directories (URL, disk name)
# URL = None for local disks
fitsdirs = [ \
    (None, "/export/data_1/users/bcotton/Software.dir/AIPS/FITS")]

# setup environment
ObitTalkUtil.SetEnviron(AIPS_ROOT=AIPS_ROOT, AIPS_VERSION=AIPS_VERSION, \
                        OBIT_EXEC=OBIT_EXEC, \
                        aipsdirs=aipsdirs, fitsdirs=fitsdirs)

# List directories
ObitTalkUtil.ListAIPSDirs()
ObitTalkUtil.ListFITSDirs()

# Any other customization goes here
