""" TaskMsgBuffer executes a task asynchronously

"""
# Task message buffer for running tasks asynchronously
# $Id$
#-----------------------------------------------------------------------
#  Copyright (C) 2008
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

# Task message buffer class
import thread, threading, time


class TaskMsgBuffer(threading.Thread):
    """ TaskMsgBuffer executes a task asynchronously
    
    Create a thread and execute a task.
    To use, create object and then call it's start() member
    a call to the wait() function will wait for program to finish
    a call to the Messages() function will return list of messages
    Messages are saved in member messages
    """
    
    def __init__(self, TaskObj):
        """ Create and initialize Task message buffer
        
        TaskObj  = Obit or AIPS or similar task object of the task to execute
        """
        self.myTask        = TaskObj
        self.number        = 1
        self.messages      = []
        self.proxy         = None  
        self.tid           = None
        self.done          = False
        self.Failed        = False
        self.Started       = False
        threading.Thread.__init__(self)
        # end __init__
        
    def run(self):
        """ Execute and manage task
        """
        ################################################################
        TaskBuf           = self
        myTask            = self.myTask
        TaskBuf.done      = False
        TaskBuf.Started   = True
        print "DEBUG TaskMsgBuffer.run"
        (TaskBuf.proxy, TaskBuf.tid) = myTask.spawn()
        # Logging to file?
        if len(myTask.logFile)>0:
            TaskLog = file(myTask.logFile,'a')
        else:
            TaskLog = None
        TaskBuf.Failed = False
        try:
            while not myTask.finished(TaskBuf.proxy, TaskBuf.tid):
                messages = myTask.messages(TaskBuf.proxy, TaskBuf.tid)
                if messages:
                    for message in messages:
                        if type(message)==str:
                            TaskBuf.messages.append(message)
                        else:
                            TaskBuf.messages.append(message[1])
                    if TaskLog:
                        for message in messages:
                            if type(message)==str:
                                x=TaskLog.write('%s\n' % message)
                            else:
                                x=TaskLog.write('%s\n' % message[1])
                        TaskLog.flush()
                continue
        except KeyboardInterrupt, exception:
            myTask.abort(TaskBuf.proxy, TaskBuf.tid)
            raise exception
        except Exception, exception:   # Aborts throw exceptions that get caught here
            TaskBuf.Failed = True
            print exception

            if not TaskBuf.Failed:
                TaskBuf.wait()
            TaskBuf.messages.append("Task Finished")
        else:
            pass
        TaskBuf.done = True
        if TaskLog:
            TaskLog.close()
        # end run
        
    def wait(self):
        """ wait for task to end
        """
        # Wait for it to start first
        while not self.Started:
            time.sleep(0.2)
        time.sleep(0.1)  # Another to be sure
        #print "in wait", self.proxy, self.tid
        if self.done:
            return
        # If the task is aborted, this will throw an exception
        try:
            self.myTask.wait(self.proxy, self.tid)
        except:
            print exception
            pass
        # end wait

    def Messages(self):
        """ Return all messages
        """
        return self.messages
    # end messages

    def abort(self):
        """ Abort task
        """
        # already done?
        if self.done:
            return
        # Abort task
        self.myTask.abort(self.proxy, self.tid)
        self.done = True
        # end abort
    # end TaskMsgBuffer class

