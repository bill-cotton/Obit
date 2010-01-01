# Copyright (C) 2005 Joint Institute for VLBI in Europe
# Copyright (C) 2007 Associated Universities, Inc
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

This module provides the bits and pieces to implement AIPSImage and
AIPSUVData objects.

"""

# Bits from Obit.
import Obit, OErr, OSystem, ODisplay
import AIPSDir
import Image, UV
import TableList

class AIPSData:
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
            if self.doInit:   # Initialized Obit?
                OSystem.Shutdown(self.ObitSys)
            return False
        if self.doInit:   # Initialized Obit?
            OSystem.Shutdown(self.ObitSys)
        return True

    def _verify(self, desc):
        # Initialize Obit if needed
        if not OSystem.PIsInit() and desc:
            AIPSdirs = desc["dirs"]
            userno   = desc["userno"]
            popsno   = 1
            if not self.err:
                self.err=OErr.OErr()
            self.ObitSys=OSystem.OSystem ("", popsno, userno, \
                                          len(AIPSdirs), AIPSdirs, \
                                          0, [], True, False, self.err)
            OErr.printErrMsg(self.err, "Error with Obit startup")
            self.doInit = True
        else:
            self.doInit = False
        try:
            data = self._init(desc)
        except OErr.OErr, err:
            print err
            #OErr.printErrMsg(err, "AIPSData._verify")
            if self.doInit:   # Initialized Obit?
                OSystem.Shutdown(self.ObitSys)
            raise RuntimeError, "Cannot open data set %s" % desc['name']
        #OErr.printErrMsg(self.err, "AIPSData._verify")
        return data

    def verify(self, desc):
        data = self._verify(desc)
        if self.doInit:   # Initialized Obit?
            OSystem.Shutdown(self.ObitSys)
            self.doInit = False
        return True                # Return something other than None.

    def header(self, desc):
        data = self._verify(desc)
        retval = data.Desc.Dict
        if self.doInit:   # Initialized Obit?
            OSystem.Shutdown(self.ObitSys)
            self.doInit = False
        OErr.printErrMsg(self.err, "Error with Descriptor")
        return retval

    def clearstat(self, desc, code):
        data = self._verify(desc)
        retval = data.Desc.Dict
        if self.doInit:   # Initialized Obit?
            OSystem.Shutdown(self.ObitSys)
            self.doInit = False
        OErr.printErrMsg(self.err, "Error with Descriptor")
        userno = desc["userno"]
        adisk  = desc['disk']
        cno = data.Acno
        e = Obit.AIPSDirStatus(adisk, userno, cno, code, self.err.me)
        return retval

    def tables(self, desc):
        data = self._verify(desc)
        retval = TableList.PGetList(data.TableList, self.err)
        if self.doInit:   # Initialized Obit?
            OSystem.Shutdown(self.ObitSys)
            self.doInit = False
        OErr.printErrMsg(self.err, "Error with Obit Table")
        return retval

    def table_highver(self, desc, type):
        data = self._verify(desc)
        retval = TableList.PGetHigh(data.TableList, type)
        if self.doInit:   # Initialized Obit?
            OSystem.Shutdown(self.ObitSys)
            self.doInit = False
        OErr.printErrMsg(self.err, "Error with Obit Table")
        return retval

    def zap(self, desc):
        data = self._verify(desc)
        data.Zap(self.err)
        if self.doInit:   # Initialized Obit?
            OSystem.Shutdown(self.ObitSys)
            self.doInit = False
        OErr.printErrMsg(self.err, "Error with Zapping data")
        return True                # Return something other than None.

    def header_table(self, desc, type, version):
        data = self._verify(desc)
        try:
            table = data.NewTable(1, type, version, self.err)
        except OErr.OErr, err:
            OErr.printErrMsg(err, "AIPSData.header_table")
            if self.doInit:   # Initialized Obit?
                OSystem.Shutdown(self.ObitSys)
                self.doInit = False
            msg = "Cannot open %s table version %d", (type, version)
            raise RuntimeError, msg
        retval = table.Desc.Dict
        if self.doInit:   # Initialized Obit?
            OSystem.Shutdown(self.ObitSys)
            self.doInit = False
        OErr.printErrMsg(self.err, "Error with Obit Table")
        return retval

    def getrow_table(self, desc, type, version, rowno):
        data = self._verify(desc)
        try:
            table = data.NewTable(1, type, version, self.err)
            table.Open(3, self.err)
        except OErr.OErr, err:
            OErr.printErrMsg(err, "AIPSData.getrow_table")
            if self.doInit:   # Initialized Obit?
                OSystem.Shutdown(self.ObitSys)
                self.doInit = False
            msg = "Cannot open %s table version %d", (type, version)
            raise RuntimeError, msg
        retval = table.ReadRow(rowno, self.err)
        if self.doInit:   # Initialized Obit?
            OSystem.Shutdown(self.ObitSys)
            self.doInit = False
        OErr.printErrMsg(self.err, "Error with Obit Table")
        return retval

    def zap_table(self, desc, type, version):
        data = self._verify(desc)
        try:
            data.ZapTable(type, version, self.err)
            data.UpdateTables(self.err)
        except OErr.OErr, err:
            OErr.printErrMsg(err, "AIPSData.zap_table")
            if self.doInit:   # Initialized Obit?
                OSystem.Shutdown(self.ObitSys)
                self.doInit = False
            msg = "Cannot zap %s table version %d", (type, version)
            raise RuntimeError, msg
        if self.doInit:   # Initialized Obit?
            OSystem.Shutdown(self.ObitSys)
            self.doInit = False
        OErr.printErrMsg(self.err, "Error with Obit Table")
        return True                # Return something other than None.
    
    pass


class AIPSImage(AIPSData):
    def _init(self, desc, verbose=True):
        userno = OSystem.PGetAIPSuser()
        OSystem.PSetAIPSuser(desc['userno'])
        image = Image.newPAImage(desc['name'], desc['name'], desc['klass'],
                                 desc['disk'], desc['seq'], True, self.err,
                                 verbose = verbose)
        OSystem.PSetAIPSuser(userno)
        if not image.isOK:  # Exception if something went wrong
            raise OErr.OErr
        OErr.printErrMsg(self.err, "Error with AIPSImage")
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
            OErr.printErrMsg(err, "AIPSImage Display")
            if self.doInit:   # Initialized Obit?
                OSystem.Shutdown(self.ObitSys)
            msg = "Cannot display image"
            raise RuntimeError, msg
        if self.doInit:   # Initialized Obit?
            OSystem.Shutdown(self.ObitSys)
            self.doInit = False
        OErr.printErrMsg(self.err, "Error with Obit Image display")
        return True           # Return something other than None.

    pass


class AIPSUVData(AIPSData):
    def _init(self, desc, verbose=True):
        userno = OSystem.PGetAIPSuser()
        OSystem.PSetAIPSuser(desc['userno'])
        uvdata = UV.newPAUV(desc['name'], desc['name'], desc['klass'],
                            desc['disk'], desc['seq'], True, self.err,
                                 verbose = verbose)
        OSystem.PSetAIPSuser(userno)
        if not uvdata.isOK:  # Exception if something went wrong
            raise OErr.OErr
        OErr.printErrMsg(self.err, "Error with AIPSUVdata")
        return uvdata

    pass


class AIPSCat:
    def __init__(self):
        self.err = OErr.OErr()
        #self.sys = OSystem.OSystem("ObitTalk", 1, 1, -1, [], -1, [],
        #                           True, False, self.err)
        return

    def cat(self, disk, userno, url, AIPSdirs):
        # Init Obit if needed
        if not OSystem.PIsInit():
            popsno = 1
            if not self.err:
                self.err=OErr.OErr()
            ObitSys=OSystem.OSystem ("", popsno, userno, \
                                     len(AIPSdirs), AIPSdirs, \
                                     0, [], True, False, self.err)
            OErr.printErrMsg(self.err, "Error with Obit startup")
            doInit = True
        else:
            doInit = False
        _userno = OSystem.PGetAIPSuser()
        OSystem.PSetAIPSuser(userno)

        try:
            num_slots = AIPSDir.PNumber(disk, userno, self.err)
        except OErr.OErr, err:
            OErr.PClear(err)
            if doInit:   # Initialized Obit?
                OSystem.Shutdown(ObitSys)
                self.doInit = False
            return []

        catalog = []
        for slot in xrange(1, num_slots):
            entry = AIPSDir.PInfo(disk, userno, slot, self.err)
            if entry:
                catalog.append((slot, entry))
                pass
            continue
        # Restore Obit to initial state
        OSystem.PSetAIPSuser(_userno)
        if doInit:   # Initialized Obit?
            OSystem.Shutdown(ObitSys)
            self.doInit = False
        OErr.printErrMsg(self.err, "Error with AIPS Catalog")
        return catalog

    pass
