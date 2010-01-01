# Copyright (C) 2005 Joint Institute for VLBI in Europe
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

This module provides the AIPSTV class.  This class makes it possible
to control the AIPS TV from Python.

"""

# Generic Python stuff.
import os, socket, struct, time


class AIPSTV(object):
    def __init__(self, host='localhost'):
        self.host = host
        self.port = socket.getservbyname('sssin', 'tcp')
        self._lock_pid = 0
        self._server_pid = 0
        return

    def _send(self, opcode, dat1=0, dat2=0, dat3=0, dat4=0):
        """Send command to the AIPS TV server."""

        s = struct.pack('!6h', opcode, dat1, dat2, dat3, dat4, 0)
        self._socket.send(s)
        if opcode == 12:
            return
        s = self._socket.recv(4)
        (status,) = struct.unpack('!h2x', s)
        if status:
            msg = "AIPS TV returned %d", status
            raise IOError, msg
        return

    def _open(self):
        """Open connection to the AIPS TV server."""

        self._socket = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
        self._socket.connect((self.host, self.port))
        self._send(11)
        return

    def _close(self):
        """Close connection to the AIPS TV server."""

        self._socket.send(struct.pack('!6h', 12, 0, 0, 0, 0, 0))
        self._socket.close()
        return

    def start(self):
        """Start the AIPS TV server."""

        # Check if we already started the AIPS TV.
        if self._lock_pid and self._server_pid:
            raise RuntimeError, "the AIPS TV has already been started"

        # Create an environment for the AIPS TV processes.
        env = os.environ.copy()
        env['TVDEV'] = 'TVDEV01'
        env['TVDEV01'] = 'sssin:localhost'

        # Start the AIPS TV lock daemon.
        file = env['LOAD'] + '/TVSERV.EXE'
        self._lock_pid = os.spawnve(os.P_NOWAIT, file, ['TVSERVER'], env)

        # Start the AIPS TV server.
        file = env['LOAD'] + '/XAS'
        self._server_pid = os.spawnve(os.P_NOWAIT, file, ['XAS'], env)

        # Wait until the TV server has been started.
        time.sleep(2)
        return

    def clear(self):
        """Init the AIPS TV server."""

        self._open()
        self._send(15)
        self._close()
        return

    def kill(self):
        """Close down the AIPS TV server."""

        self._open()
        self._socket.send(struct.pack('!6h', 18, 0, 0, 0, 0, 0))
        self._socket.close()

        # Kill of the zombies.
        waited = False
        if self._lock_pid:
            os.waitpid(self._lock_pid, 0)
            self._lock_pid = 0
            waited = True
            pass
        if self._server_pid:
            os.waitpid(self._server_pid, 0)
            self._server_pid = 0
            waited = True
            pass
        if not waited:
            # Take a nap to avoid confusing users with the output of
            # the dying processes.
            time.sleep(2)

        return

    pass                                # Class AIPSTV
