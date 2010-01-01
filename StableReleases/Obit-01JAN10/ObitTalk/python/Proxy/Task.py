"""

This module provides bits and pieces to implement generic Task proxy
objects.

"""
# Copyright (C) 2005 Joint Institute for VLBI in Europe
# Copyright (C) 2006,2007 Associated Universities, Inc. Washington DC, USA.
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


# Generic Python stuff.
import fcntl, os, pty, select, signal, time

class Task:
    """ Basic server Task class """
    def __init__(self):
        self._pid = {}

    def spawn(self, path, args, env=None):
        """Spawn the task.
         
        There appears to be two calls, the first in this something undocumented
        is done to the terminal and the sceond to actually start the function
        Returns tid
        path   path to executable
        args   arguments to task
        env    environment dictionary
      `  """

        # Fork
        (pid, tid) = pty.fork()
        if pid == 0:   # Actually start function
            try:
                if env:
                    os.execvpe(path, args, env)
                else:
                    os.execv(path, args)
            except Exception, exception:
                print exception
            os._exit(1)
        else:
            # Setup terminal
            fcntl.fcntl(tid, fcntl.F_SETFL, os.O_NONBLOCK)
            self._pid[tid] = pid
            return tid

    def finished(self, tid):
        """Check whether the task has finished.
        
        tid   = Task id in pid table of process
        """
        # Check that tid in bounds
        #if tid>len(self._pid):
        #    return True   # finished one way or another
        time.sleep(0.1)
        return self._pid[tid] == 0

    def messages(self, tid):
        """Return task messages
        
        Returns list of messages
        tid   = Task id in pid table of process
        """

        # Check that tid in bounds
        #if tid>len(self._pid):
        #    return []

        try:
            (iwtd, owtd, ewtd) = select.select([tid], [], [], 0.25)
        except:
            return []
        if tid in iwtd:
            try:
                # hang 10 - don't get ahead of yourself
                #time.sleep(0.010)
                messages = os.read(tid, 1024)
                # Test for task completion here - Macs never fail on read
                (pid, status) = os.waitpid(self._pid[tid], os.WNOHANG)
                if pid:
                    assert(pid == self._pid[tid])
                    if os.WIFEXITED(status) or os.WIFSIGNALED(status):
                        self._pid[tid] = 0
                #return [msg for msg in messages if msg]
                return [messages]
            except:
                # If reading failed, it's (probably) because the child
                # process died.
                (pid, status) = os.waitpid(self._pid[tid], os.WNOHANG)
                if pid:
                    assert(pid == self._pid[tid])
                    if os.WIFEXITED(status) or os.WIFSIGNALED(status):
                        self._pid[tid] = 0
        return []

    def feed(self, tid, banana):
        """Feed the task a  BANANA.
        
        Pass a message to a running task's sdtin
        tid     = Task id in pid table of process
        bananna = text message to pass to task input
        """
        
        os.write(tid, banana)
        pass
    
    def wait(self, tid):
        """Wait for the task to finish.
        
        tid   = Task id in pid table of process
        """

        assert(self.finished(tid))

        del self._pid[tid]
        os.close(tid)
        return

    def abort(self, tid, sig=signal.SIGKILL):
        """Abort the task
        
        Calls abort function for task tid
        None return value
        tid   = Task id in pid table of process to be terminated
        sig   = signal to sent to the task
        """
        #print "DEBUG kill pid, sig",self._pid[tid],sig
        os.kill (self._pid[tid], sig)

        del self._pid[tid]
        os.close(tid)
        return
    # End Class Task


