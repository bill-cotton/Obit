# The argument, if given, is the data directory, defaults to "../testIt"
import OErr, OSystem, Image, InfoList, FArray, sys
from OErr import Bomb

if len(sys.argv)>=2:
    dataDir = sys.argv[1]
else:
    dataDir = "../testIt/"
err=OErr.OErr()
ObitSys=OSystem.OSystem ("Python", 1, 103, 1, ["None"], 1, [dataDir], 1, 0, err)
OErr.printErrMsg(err, "Error with Obit startup")

# for debugging
#Bomb()

list = InfoList.InfoList();
dim = [1,0,0,0,0]
InfoList.PPutLong (list, "longItem", dim, [321], err)
x=InfoList.PGet (list, "longItem")

x=FArray.FArray("MyFArray",[30,50])
FArray.PSetVal(x,[25,25],3.1415926)
pos=[-1,-1]
result = FArray.PMax(x,pos)
print "FArray max is",result, "should be 3.1415926"
result = FArray.PMin(x,pos)
print "FArray min is",result, "should be 0"

# Try printing Memory allocation
OSystem.PMemPrint()

file = "testPcube.fits"
disk = 1
image = Image.newPFImage("myImage", file, disk, True, err);
OErr.printErr(err)
Image.POpen(image,1,err)
OErr.printErr(err)
Image.PRead(image,err)
OErr.printErr(err)
Image.PClose(image,err)
data = Image.PGetFArray(image);
pos=[-1,-1]
result = FArray.PMax(data,pos)
print "Max in",file,"is",result
result = FArray.PRMS(data)
print "RMS in",file,"is",result

# Shutdown Obit
OErr.printErr(err)
OSystem.Shutdown(ObitSys)
