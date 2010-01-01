# Program to select a section of OTF data
import OTF, OTFUtil, OSystem, OErr, InfoList

# Init Obit
err=OErr.OErr()
ObitSys=OSystem.OSystem ("FitCal", 1, 103, 1, ["None"], 1, ["../PythonData/"], 1, 0, err)
OErr.printErrMsg(err, "Error with Obit startup")

# Files
disk = 1
#inFile   = "GCXbandDay3OTF.fits"         # input OTF data
#inFile   = "GCCbandDay1OTF.fits"         # input OTF data
#inFile   = "GCLband3C286OTF.fits"
inFile = "QbandOTF.fits"

# Set data
inData  = OTF.newPOTF("Input data",  inFile, disk, 1, err)
OErr.printErrMsg(err, "Error initializing")

# Set scans
#GCXbandDay3OTF.fits scan = [3004,3005,3006,3007]
scan = [2007,2008,2009,2010] #GCCbandDay2OTF.fits
scan = [2003,2004,2005,2006] #GCCbandDay2OTF.fits
scan = [1005,1006,1007,1008] #GCCbandDay1OTF.fits
scan = [16,17,18,19]#GCLband3C286OTF.fits"
scan = [80,81,82,83]
tau0 = 0.01 # Xband 
tau0 = 0.008 # Cband
tau0 = 0.006 # Lband
tau0=0.03 #Qband - wild guess
ATemp = [.780, .790]       # Xband day 3
ATemp = [0.67205,0.803175] # C band day 1
ATemp = [0.60122,0.73162] # C band day 1
ATemp = [1.0,1.0] # L band day both- wild guess
OTF.AtmCalInput["aTemp"]  = [0.735,0.722,0.722,0.722,0.722,0.722,0.722,0.722]  # wild guess
dim = OTF.dim
dim[0] = 1; dim[1] = 1
inInfo = OTF.PGetList(inData)
InfoList.PAlwaysPutFloat(inInfo, "TAU0",  dim, [tau0])
dim[0] = 2
InfoList.PAlwaysPutFloat (inInfo,"ATEMP" , dim, ATemp)
dim[0] = 4
InfoList.PAlwaysPutInt(inInfo, "Scan",  dim, scan)

# Do fitting
OTFUtil.PFitCal(inData, -1,  err)
OErr.printErrMsg(err, "Error fitting cal scan")

# Give results
print "Fitting scans",scan,"in",inFile
stuff = InfoList.PGet (inInfo, "TRX")
print "Average TRx = ",stuff[4]
stuff = InfoList.PGet (inInfo, "CALJY")
print "Average CalJy = ",stuff[4]
stuff = InfoList.PGet (inInfo, "Timeoff")
print "Average Timeoff = ",stuff[4]
print "Use opposite sign of Timeoff in DCR2OTR input"

# Shutdown Obit
OErr.printErr(err)

