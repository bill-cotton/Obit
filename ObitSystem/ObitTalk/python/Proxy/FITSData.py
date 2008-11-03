# Copyright (C) 2008 Associated Universities, Inc
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

"""

This module provides the bits and pieces to implement FITSImage and
FITSUVData objects.

"""

# Bits from Obit.
import Obit, OErr, OSystem, ODisplay
import FITSDir
import Image, UV
import TableList
import os

class FITSData:
    def __init__(self, desc=None):
        self.err     = OErr.OErr()
        self.doInit  = False
        self.ObitSys = None
        return

    def exists(self, desc):
        try:
            self._init(desc, verbose=False)
        except OErr.OErr, err:
            OErr.PClear(err)
            return False
        return True

    def _verify(self, desc):
        # Initialize Obit if needed
        if not self.err:
            self.err=OErr.OErr()
        try:
            data = self._init(desc)
        except OErr.OErr, err:
            print err
            #OErr.printErrMsg(err, "FITSData._verify")
            raise RuntimeError, "Cannot open data set %s" % desc['filename']
        #OErr.printErrMsg(self.err, "FITSData._verify")
        return data

    def verify(self, desc):
        data = self._verify(desc)
        return True                # Return something other than None.

    def header(self, desc):
        data = self._verify(desc)
        retval = data.Desc.Dict
        OErr.printErrMsg(self.err, "Error with Descriptor")
        return retval

    def tables(self, desc):
        data = self._verify(desc)
        retval = TableList.PGetList(data.TableList, self.err)
        OErr.printErrMsg(self.err, "Error with Obit Table")
        return retval

    def table_highver(self, desc, type):
        data = self._verify(desc)
        retval = TableList.PGetHigh(data.TableList, type)
        OErr.printErrMsg(self.err, "Error with Obit Table")
        return retval

    def zap(self, desc):
        data = self._verify(desc)
        data.Zap(self.err)
        OErr.printErrMsg(self.err, "Error with Zapping data")
        return True                # Return something other than None.

    def header_table(self, desc, type, version):
        data = self._verify(desc)
        try:
            table = data.NewTable(1, type, version, self.err)
        except OErr.OErr, err:
            OErr.printErrMsg(err, "FITSData.header_table")
            if self.doInit:   # Initialized Obit?
                OSystem.Shutdown(self.ObitSys)
                self.doInit = False
            msg = "Cannot open %s table version %d", (type, version)
            raise RuntimeError, msg
        retval = table.Desc.Dict
        OErr.printErrMsg(self.err, "Error with Obit Table")
        return retval

    def getrow_table(self, desc, type, version, rowno):
        data = self._verify(desc)
        try:
            table = data.NewTable(1, type, version, self.err)
            table.Open(3, self.err)
        except OErr.OErr, err:
            OErr.printErrMsg(err, "FITSData.getrow_table")
            if self.doInit:   # Initialized Obit?
                OSystem.Shutdown(self.ObitSys)
                self.doInit = False
            msg = "Cannot open %s table version %d", (type, version)
            raise RuntimeError, msg
        retval = table.ReadRow(rowno, self.err)
        OErr.printErrMsg(self.err, "Error with Obit Table")
        return retval

    def zap_table(self, desc, type, version):
        data = self._verify(desc)
        try:
            data.ZapTable(type, version, self.err)
            data.UpdateTables(self.err)
        except OErr.OErr, err:
            OErr.printErrMsg(err, "FITSData.zap_table")
            if self.doInit:   # Initialized Obit?
                OSystem.Shutdown(self.ObitSys)
                self.doInit = False
            msg = "Cannot zap %s table version %d", (type, version)
            raise RuntimeError, msg
        OErr.printErrMsg(self.err, "Error with Obit Table")
        return True                # Return something other than None.
    
    pass


class FITSImage(FITSData):
    def _init(self, desc, verbose=True):
        # Open with full path in disk 0
        disk = 0
        path = desc['dirname']+"/"+desc['filename']
        image = Image.newPFImage(desc['filename'], path, disk, 
                                 True, self.err,
                                 verbose = verbose)
        if not image.isOK:  # Exception if something went wrong
            raise OErr.OErr
        OErr.printErrMsg(self.err, "Error with FITSImage")
        return image

    def display(self, desc, url):
        data = self._verify(desc)
        try:
            # Display server object
            disp = ODisplay.ODisplay("ObitView", url, self.err)
            # Send to display
            ODisplay.PImage(disp, data, self.err)
            del disp
        except OErr.OErr, self.err:
            OErr.printErrMsg(err, "FITSImage Display")
            msg = "Cannot display image"
            raise RuntimeError, msg
        OErr.printErrMsg(self.err, "Error with Obit Image display")
        return True           # Return something other than None.

    pass


class FITSUVData(FITSData):
    def _init(self, desc, verbose=True):
        # Open with full path in disk 0
        disk = 0
        path = desc['dirname']+"/"+desc['filename']
        uvdata = UV.newPFUV(desc['filename'], path, disk, 
                            True, self.err,
                                 verbose = verbose)
        if not uvdata.isOK:  # Exception if something went wrong
            raise OErr.OErr
        OErr.printErrMsg(self.err, "Error with FITSUVdata")
        return uvdata

    pass


class FITSCat:
    def __init__(self):
        #self.err = OErr.OErr()
        return

    def cat(self, disk, dir, userno, url, FITSdirs):
        catalog = []
        if dir==None:
            dir = "./"
        if (disk==0):
            listDir = dir
        else:
            listDir = FITSdirs[disk-1]
        flist = os.listdir(listDir)
        slot = 1
        for f in flist:
            catalog.append((slot, f))
            slot += 1
        return catalog

    pass
