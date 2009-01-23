""" GBT Obit OTF Utilities

This utility package contains utilities of relevance to the GBT OTF data
Some of these must be run from ObitTalk.
"""
# $Id$
#-----------------------------------------------------------------------
#  Copyright (C) 2007,2008
#  Associated Universities, Inc. Washington DC, USA.
#
#  This program is free software; you can redistribute it and/or
#  modify it under the terms of the GNU General Public License as
#  published by the Free Software Foundation; either version 2 of
#  the License, or (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public
#  License along with this program; if not, write to the Free
#  Software Foundation, Inc., 675 Massachusetts Ave, Cambridge,
#  MA 02139, USA.
#
#  Correspondence concerning this software should be addressed as follows:
#         Internet email: bcotton@nrao.edu.
#         Postal address: William Cotton
#                         National Radio Astronomy Observatory
#                         520 Edgemont Road
#                         Charlottesville, VA 22903-2475 USA
#-----------------------------------------------------------------------
import OTF, Image, Table, ObitTask, FITSDir, OErr

# Routine to update OTF
def UpdateOTF (OTFTask, Rcvr, outFile, outDisk, DataRoot, err, \
               scanList=None, offTime=None, avgTime=None, config=None, \
               scanNo=None, doQuadD=False, doBS=None, dataNorm=None):
    """ Update OTF and return object

    return Python OTF object
    Update an OTF object with scans not yet in the file.
    Read DataRoot/scanLog.fits and determine which scans are available and
    are not already in the output OTF.
    These scans are appended and the output OTF returned.
    Does some data validity checks.

    OTFTask  = Name of Obit task to read GBT archive and write to OTF format
               DCR = "DCROTF"
               PAR = "PAROTF" (Mustang)
               CCB = "CCBOTF"
    Rcvr     = directory name with data for receiver.
               DCR = "DCR"
               PAR = "Rcvr_PAR" (Mustang)
               CCB = "CCB26_40"
    outFile  = Name of output OTF FITS file
    outDisk  = disk number for outFile, 0=> current working directory
    DataRoot = Root of GBT archive for current project
               If None, don't attempt
    err      = Python Obit Error/message stack

    Optional, backend specific values
    scanList = All, list of scan numbers, None=>all
    offTime  = PAR, CCB Offset in sec to be added to time
    avgTime  = PAR Data averaging time in seconds
    config   = PAR  path of configuration file
    scanNo   = PAR, CCB, replace GBT scan number with this value
    doQuadD  = PAR, Use Quadrant Detector offsets
    doBS     = CCB, output beamswitched data
    dataNorm = CCB normalization factors for beamswitched data
    """
    ################################################################
    # Checks
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"

    # Anything to do?
    if not DataRoot:
        outOTF=OTF.newPOTF("output", outFile, outDisk, True, err, nrec=1000)
        return outOTF

    # Scan list 
    if scanList:
        sList = scanList
    else:
        sList = range(1,1000000)
   
    # How much data is already there?
    maxscan=0  # max scan number found so far
    outOTF = None
    if  FITSDir.PExist(outFile, outDisk, err):
        outOTF=OTF.newPOTF("output", outFile, outDisk, True, err)
        indx = outOTF.NewTable(Table.READONLY,"OTFIndex",1,err)
        indx.Open(Table.READONLY,err)
        norow = indx.Desc.Dict["nrow"]
        indxRow=indx.ReadRow(norow,err)
        maxscan = indxRow["SCAN_ID"][0]
        indx.Close(err)
    
    OErr.printErrMsg(err, "Error with Output OTF file")
    print "Highest previous scan is",maxscan

    # Task to read/append
    otf=ObitTask.ObitTask(OTFTask)
    otf.DataRoot = DataRoot
    otf.outOTF   = outFile
    otf.outDisk  = outDisk
    # Optional parameters
    if offTime:
        otf.offTime = offTime
    if avgTime:
        otf.avgTime = avgTime
    if config:
        otf.config = config
    if scanNo:
        otf.scanNo = scanNo
    if doQuadD:
        otf.doQuadD = doQuadD
    if doBS:
        otf.doBS = doBS
    if dataNorm:
        otf.dataNorm = dataNorm

    img=Image.newPFImage("Scan log",DataRoot+"ScanLog.fits",0, True, err)
    scanLog=img.NewTable(Table.READONLY,"ScanLog",1,err)
    scanLog.Open(Table.READONLY,err)
    OErr.printErrMsg(err, "Error with ScanLog")
    norow = scanLog.Desc.Dict["nrow"]
    scnstr=[]    # list of scan numbers
    for irow in range (1,norow+1):
        scan=scanLog.ReadRow(irow,err)
        scanno = scan["SCAN"][0]
        if (scanno>maxscan) and (scanno in sList):
            # Take this one
            maxscan = scanno
            datetime = scan["DATE-OBS"][0]
            # Table ScanLog, multiple entries per scan, need to translate
            # timestring to scan, "-" => "_", "T" => "_"
            scanstr=datetime.replace("-","_").replace("T","_").rstrip()
            # Make sure there is actual data
            test = DataRoot+Rcvr+"/"+scanstr+".fits"
            OK = FITSDir.PExist(test, 0, err) or FITSDir.PExist(test+".gz", 0, err)
            test = DataRoot+"GO/"+scanstr+".fits"
            OK = OK and FITSDir.PExist(test, 0, err) or FITSDir.PExist(test+".gz", 0, err)
            test = DataRoot+"Antenna/"+scanstr+".fits"
            #print "antenna",test
            OK = OK and FITSDir.PExist(test, 0, err) or FITSDir.PExist(test+".gz", 0, err)
            #print "DEBUG scan",scanstr, Rcvr, OK
            if  OK:
                scnstr.append(scanstr)
                # Run PAROTF to append this scan
                otf.Scan=scanstr
                try:
                    #otf.debug=True
                    otf.g
                except:
                    print "failed on scan",scanstr
                else:
                    pass
            
    # Now have list of scans>maxscan in scnstr
    scanLog.Close(err)
    OErr.printErrMsg(err, "Error with ScanLog")

    print "Added Scans",scnstr

    # Output if just created
    if not outOTF:
        outOTF=OTF.newPOTF("output", outFile, outDisk, True, err)
        OErr.printErrMsg(err, "Error with Output OTF file")
    return outOTF
    # end UpdateOTF

# Get target position from OTFTarget table
def GetTargetPos (OTF, Target, err):
    """ Get target position from OTFTarget table

    return [raepo, decepo] in deg
    Loop through target table on OTF looking for Target and return position
    
    OTF      = OTF data file to check
    Target   = Name of target e.g. "MARS"
    """
    ################################################################
    # Checks
    if not OTF.OTFIsA():
        raise TypeError,"OTF MUST be a Python Obit OTF"
    if  type(Target)!=str:
        raise TypeError,"Target MUST be a string"
    if not OErr.OErrIsA(err):
        raise TypeError,"err MUST be an OErr"

    tarTab=OTF.NewTable(Table.READONLY,"OTFTarget",1,err)
    tarTab.Open(Table.READONLY,err)
    OErr.printErrMsg(err, "Error with OTFTarget table")
    norow = tarTab.Desc.Dict["nrow"]
    for irow in range (1,norow+1):
        row=tarTab.ReadRow(irow,err)
        if row["TARGET"][0].strip() == Target.strip():
            out = [row["RAEPO"][0], row["DECEPO"][0]]
            tarTab.Close(err)
            return out
    
    tarTab.Close(err)
    # Didn't find it if it gets here
    OErr.printErrMsg(err, "Error getting Target")
    raise RuntimeError,"Failed to find target "+Target+" in OTFTarget on "+OTF.GetName()
    return [None, None]
    # end GetTargetPos
   

