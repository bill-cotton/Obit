# Copyright (C) 2005, 2006 Joint Institute for VLBI in Europe
# Copyright (C) 2007  Associated Universities, Inc
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

This module provides a simple XML-RPC server that can run classic AIPS
and Obit tasks and provide verb-like access to AIPS data on a machine.

"""

from SimpleXMLRPCServer import SimpleXMLRPCServer
from SocketServer import ThreadingMixIn
import sys

# Global AIPS defaults.
import AIPS
import OSystem

# Import AIPS modules.
import Proxy.ObitTask
import Proxy.AIPSTask
import Proxy.AIPSData
import Proxy.ObitScriptP
ObitScript = Proxy.ObitScriptP.ObitScript()

# Use port 8000 or first command line argument
print sys.argv
if len(sys.argv)>1:
    port = sys.argv[1]
else:
    port=8000

# debug diagnostics
#print "Obit Init",OSystem.PIsInit ()
#print "AIPS user no.",OSystem.PGetAIPSuser()


class XMLRPCServer(SimpleXMLRPCServer):
    allow_reuse_address = True
    pass

class ServerFuncs:
    def __init__(self):
        self.ObitTask   = Proxy.ObitTask.ObitTask()
        self.AIPSTask   = Proxy.AIPSTask.AIPSTask()
        self.AIPSUVData = Proxy.AIPSData.AIPSUVData()
        self.AIPSImage  = Proxy.AIPSData.AIPSImage()
        self.AIPSCat    = Proxy.AIPSData.AIPSCat()
        self.ObitScript = Proxy.ObitScriptP.ObitScript()
        return

    def _dispatch(self, name, args):
        # For security reasons, SimpleXMLRPCServer in Python
        # 2.3.5/2.4.1, no longer resolves names with a dot in it.  Se
        # here we explicitly accept names starting with 'AIPS'
        # or 'Obit' and containing a single dot; that should be safe enough. 
        #print "received call ",name,args
        if (name.startswith('AIPS') or name.startswith('Obit')) \
               and name.count('.') == 1:
            name = name.split('.')
            inst = getattr(self, name[0])
            method = getattr(inst, name[1])
            return method(*args)
        msg = "object has no attribute '%s'" % name
        raise AttributeError, msg

    pass                                # class ServerFuncs

server = XMLRPCServer(('', port))
server.register_instance(ServerFuncs())
print "Starting ObitTalkServer on port",port
server.serve_forever()
