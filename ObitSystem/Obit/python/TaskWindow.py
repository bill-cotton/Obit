""" TaskWindow executes a task asynchronously, messages in ObitMess Window

ObitMess must be started independently of python
"""
# Task Window for running tasks asynchronously
# $Id$
#-----------------------------------------------------------------------
#  Copyright (C) 2006,2009
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
import thread, threading, time
from xmlrpclib import ServerProxy


class TaskWindow(threading.Thread):
    """ TaskWindow executes a task asynchronously, messages displayed
    
    Create a thread and execute a task.
    To use, create object and then call it's start() member
    a call to the wait() function will wait for program to finish
    Messages are displayed in a ObitMess window
    """
    
    def __init__(self, TaskObj, URL="http://localhost:8777/RPC2"):
        """ Create and initialize Task Window
        
        TaskObj  = Obit or AIPS or similar task object of the task to execute
        """
        self.myTask        = TaskObj
        self.number        = 1
        self.taskID        = -1  # ObitMess task ID
        self.done          = False
        self.proxy         = None  
        self.tid           = None
        self.Failed        = False
        self.Started       = False
        self.url           = URL
        threading.Thread.__init__(self)
        # end __init__
        
    def run(self):
        """ Execute and manage task

        """
        ################################################################
        TaskWin           = self
        myTask            = self.myTask
        TaskWin.done      = False

        # Start ObitMess Window for task myTask._name
        try:
            server = ServerProxy(self.url)
            answer = server.CreateWindow(myTask._name)
            self.taskID = answer["taskID"]
        except Exception, e:
            print "Failed to talk to ObitMess",e
            raise RuntimeError,"Cannot talk to ObitMess - start it "
        
        # Hang around until gui is started
        time.sleep(1.)
        TaskWin.Started   = True
        deadGUI = False  # has GUI died? May want to keep running/logging.
        arg = {"taskID":self.taskID, "status":"Task Running"}
        answer = server.SetStatus(arg)
        deadGUI = deadGUI or (answer['Status']['code']!=0)
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
                    if not deadGUI:
                        msg = ""
                        # Bundle messages
                        for message in messages:
                            if type(message)==str:
                                msg += message+'\n'
                            else:
                                msg += message[1]+'\n'                       
                        arg = {"taskID":self.taskID, "message":msg}
                        answer  = server.DisplayMessage(arg)
                        # loop if busy
                        while answer['Status']['reason']=="Busy":
                            answer  = server.DisplayMessage(arg)
                        deadGUI = deadGUI or (answer['Status']['code']!=0)
                        # Abort request?
                        doAbort = answer["Abort"]
                        # Input (AIPS) request in message?
                        if (msg.__contains__("** press RETURN for more") or
                            msg.__contains__("just hit RETURN to continue ")):
                            answer  = server.UserResponse(self.taskID)
                            # loop if busy
                            while answer['Status']['reason']=="Busy":
                                answer  = server.UserResponse(self.taskID)
                            deadGUI = deadGUI or (answer['Status']['code']!=0)
                            doAbort = doAbort or answer["Abort"]
                            reply   = answer["Result"]
                            if not TaskWin.done:
                                # Feed the task the command
                                myTask.feed(TaskWin.proxy, TaskWin.tid, reply+"\n")
                            else: # Task finished
                                arg = {"taskID":self.taskID, "message":"Task already finished\n"}
                                answer = server.DisplayMessage(arg)
                        # Now abort if requested
                        if doAbort:
                            time.sleep(1.)  # time to shutdown if requested
                            if not myTask.finished(self.proxy, self.tid):
                                # this seems to leave the task in an undead state
                                self.myTask.abort(self.proxy, self.tid)
                            TaskWin.Failed = True
                            TaskWin.done   = True
                        # end if GUI alive
                    if TaskLog:
                        for message in messages:
                            if type(message)==str:
                                x=TaskLog.write('%s\n' % message)
                            else:
                                x=TaskLog.write('%s\n' % message[1])
                        TaskLog.flush()
                continue
        except KeyboardInterrupt, exception:
            print "Something went wrong:",exception
            myTask.abort(TaskWin.proxy, TaskWin.tid)
            raise exception
        except Exception, e:   # Aborts throw exceptions that get caught here
            TaskWin.Failed = True
            TaskWin.done   = True
            #print "An exception was thrown, task aborted:",e

        if not TaskWin.Failed:
            TaskWin.wait()
            arg = {"taskID":self.taskID, "status":"Task Finished"}
            answer = server.SetStatus(arg)
        else:
            arg = {"taskID":self.taskID, "status":"Task Failed"}
            answer = server.SetStatus(arg)
        TaskWin.done = True

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
        server = ServerProxy(self.url)
        arg = {"taskID":self.taskID, "message":"Task Aborted\n"}
        answer = server.DisplayMessage(arg)
        arg = {"taskID":self.taskID, "status":"Task Aborted"}
        answer = server.SetStatus(arg)

        # end abort
    # end TaskWindow class

