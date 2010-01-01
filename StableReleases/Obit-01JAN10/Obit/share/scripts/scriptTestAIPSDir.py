# Test AIPS directory listing

import Obit, OSystem, OErr, AIPSDir, Image, ODisplay

# Init Obit
err=OErr.OErr()
ObitSys=OSystem.OSystem ("Imager", 1, 100, 1, ["/DATA/LUSUS_1/"], 1, ["../PythonData/"], 1, 0, err)
OErr.printErrMsg(err, "Error with Obit startup")

from Obit import Bomb

disk=1
# list directory
AIPSDir.PListDir(disk, err)

a=Image.newPACNO(1,4,True,err)
d=ODisplay.ODisplay("ObitView","ObitView",err)
ODisplay.PImage(d,a,err)



# Shutdown Obit
OErr.printErr(err)
del ObitSys
