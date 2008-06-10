""" TaskWindow executes a task asynchronously

More thoughts:
- how to extract task object from event
- what happens when task wants to talk to user and ask input?
"""
# Task Window for running tasks asynchronously
# $Id: TaskWindow.py,v 1.4 2006/07/21 14:36:42 bcotton Exp $
#-----------------------------------------------------------------------
#  Copyright (C) 2006
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

# Task window class
import thread, threading, wx, time
import OTWindow


class TaskWindow(threading.Thread):
    """ TaskWindow executes a task asynchronously
    
    Create a thread and execute a task.
    To use, create object and then call it's start() member
    a call to the wait() function will wait for program to finish
    Messages are displayed in a ScrolledText window
    NB: The root Tk window should not be closed until you are through with tasks
    """
    
    def __init__(self, TaskObj):
        """ Create and initialize Task Window
        
        TaskObj  = Obit or AIPS or similar task object of the task to execute
        """
        self.myTask        = TaskObj
        self.number        = 1
        self.MsgWinId      = -1  # wxPython ID of message window
        self.scrollTextId  = -1  # wxPython ID of message text
        self.commandId     = -1  # wxPython ID of command window
        self.statusId      = -1  # wxPython ID of status label
        self.abortButtonId = -1  # wxPython ID of abort/close button
        self.proxy         = None  
        self.tid           = None
        self.done          = False
        self.Failed        = False
        self.Started       = False
        threading.Thread.__init__(self)
        # end __init__
        
    def run(self):
        """ Execute and manage task

        Note: must send info for Message window through events
        """
        ################################################################
        import MsgWin
        TaskWin           = self
        myTask            = self.myTask
        TaskWin.done      = False
        # Hang around until gui is started
        while TaskWin.statusId==-1:
            time.sleep(0.2)
        TaskWin.Started   = True
        OTWindow.CallSetLabel (TaskWin.statusId, "Task Running")
        (TaskWin.proxy, TaskWin.tid) = myTask.spawn()
        # Logging to file?
        if len(myTask.logFile)>0:
            TaskLog = file(myTask.logFile,'a')
        else:
            TaskLog = None
        TaskWin.Failed = False
        try:
            while not myTask.finished(TaskWin.proxy, TaskWin.tid):
                messages = myTask.messages(TaskWin.proxy, TaskWin.tid)
                if messages:
                    OTWindow.CallMessage(self.scrollTextId, messages)
                    if TaskLog:
                        for message in messages:
                            if type(message)==str:
                                x=TaskLog.write('%s\n' % message)
                            else:
                                x=TaskLog.write('%s\n' % message[1])
                        TaskLog.flush()
                continue
        except KeyboardInterrupt, exception:
            myTask.abort(TaskWin.proxy, TaskWin.tid)
            raise exception
        except:   # Aborts throw exceptions that get caught here
            TaskWin.Failed = True
            #print "DEBUG in go - an exception thrown"
        #print "DEBUG in go, finished", myTask.myTask.finished(myTask.proxy, myTask.tid)

        if not TaskWin.Failed:
            TaskWin.wait()
            OTWindow.CallSetLabel(TaskWin.statusId, "Task Finished")
            OTWindow.CallUpdate(TaskWin.MsgWinId);
        else:
            OTWindow.CallSetLabel(TaskWin.statusId, "Task Failed")
            OTWindow.CallUpdate(TaskWin.MsgWinId);
            pass
        TaskWin.done = True
        # Turn Abort button into Close
        OTWindow.CallSetLabel(TaskWin.abortButtonId, "Close")
        OTWindow.CallBind(TaskWin.abortButtonId, MsgWin.delete)
        OTWindow.CallUpdate(TaskWin.MsgWinId);
        if TaskLog:
            TaskLog.close()
        #if AIPS.log:
        #    for message in log:
        #        AIPS.log.write('%s\n' % message)
        #        continue
        #    pass
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
            pass
        # end wait

    def abort(self):
        """ Abort task
        """
        #print "in abort"
        if self.done:
            return
        # Abort task
        self.myTask.abort(self.proxy, self.tid)
        self.done = True
        
        # Update message window
        OTWindow.CallMessage(self.scrollTextId, "Task aborted")
        OTWindow.CallSetLabel(self.statusId, "Task Aborted")
        
        # Turn "Abort" button into "Close"
        OTWindow.CallSetLabel(self.abortButtonId, "Close")
        OTWindow.CallBind(self.abortButtonId, MsgWin.delete)
        # end abort
    # end TaskWindow class

