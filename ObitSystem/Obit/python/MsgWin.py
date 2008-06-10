""" MsgWin ObitTalk message window class

a MsgWin executes a task for ObitTalk snd displays any messages.
"""
# $Id: MsgWin.py,v 1.1 2006/07/19 14:08:38 bcotton Exp $
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

import wx
import TaskWindow, OTWindow

# event handlers
def save(event):
    """ Save current messages into a file selected by a browser
    """
    tw = event.GetEventObject().TskWin
    fb=wx.FileDialog(None,message="Output Log file")
    x=fb.ShowModal()
    filename=fb.GetFilename()
    if len(filename) > 0:  # did it work?
        file = open(filename, 'w')
        file.write(tw.scrollText.GetValue())
        file.close()
    # end save
    
def delete(event):
    """ delete window
    """
    #print "in delete"
    # Wait for finish
    win = event.GetEventObject()
    tw = win.TskWin
    parent = win.window
    tw.wait()             # wait for task to finish
    x=parent.Destroy()    # kill da wabbit
    #del win
    # end delete
        
def abort(event):
    """ Abort task
    """
    #print "in abort"
    abortButton = event.GetEventObject()
    tw = abortButton.TskWin
    if tw.done:
        tw.scrollText.AppendText(tw.scrollTextId, "Task already finished\n")
        return
    tw.myTask.abort(tw.proxy, tw.tid)
    tw.scrollText.AppendText("Task aborted\n")
    tw.status.SetLabel("Task Aborted")
    tw.done = True
    # Turn Abort button into Close
    tw.abortButton.SetLabel("Close")
    tw.abortButton.Bind(wx.EVT_BUTTON, delete)
    # end abort

def sndcommand(event):
    """ Sent text from the command window to the task
    """
    #print "in delete"
    # Wait for finish
    win = event.GetEventObject()
    cmnd = win.text.GetValue() # read command text
    win.text.SetValue("")      # reset to blank
    #print "DEBUG send task ",cmnd,win.TskWin.done
    if not win.TskWin.done:
        # Feed the task the command
        win.TskWin.myTask.feed(win.TskWin.proxy, win.TskWin.tid, cmnd+"\n")
    else:
        win.TskWin.scrollText.AppendText("Task already finished\n")
   # end sndcommand
        

class MsgWin(wx.Frame):
    """
    Class to manage communications with asynchronous ObitTalk tasks
    """

    def __init__(self, TaskWin):
        """ Create and initialize MsgWin
        
        TaskWin  = Task window  object of the task to execute
        """
        ################################################################
        # Initilize base class
        wx.Frame.__init__(self, OTWindow.gui.app.catcher,
                          title=TaskWin.myTask._name+" Messages",
                          size=(500,400))
        self.myTaskWin        = TaskWin  # Save task window object
        self.number           = 0        # Number of messages
        TaskWin.MsgWinId      = self.GetId();  # ID for communications
        bkg                   = wx.Panel(self)
        scrollText            = wx.TextCtrl(bkg, style=wx.TE_MULTILINE|wx.HSCROLL)
        TaskWin.scrollText   = scrollText
        TaskWin.scrollTextId  = scrollText.GetId();  # ID for communications
        command               = wx.TextCtrl(bkg)
        sendButton            = wx.Button(bkg,label="Send")
        sendButton.text       = command
        sendButton.TskWin     = TaskWin   # allow access when button hit
        sendButton.Bind(wx.EVT_BUTTON, sndcommand)
        TaskWin.commandId     = command.GetId();  # ID for communications
        abortButton           = wx.Button(bkg,label="Abort")
        abortButton.Bind(wx.EVT_BUTTON, abort)
        abortButton.TskWin    = TaskWin   # allow access when button hit
        abortButton.window    = self      # allow access when button hit
        TaskWin.abortButtonId = abortButton.GetId();  # ID for communications
        TaskWin.abortButton   = abortButton
        saveButton            = wx.Button(bkg,label="Save as")
        saveButton.Bind(wx.EVT_BUTTON, save)
        saveButton.TskWin     = TaskWin   # allow access when button hit
        status                = wx.StaticText(bkg, label="status")
        TaskWin.statusId      = status.GetId();  # ID for communications
        TaskWin.status        = status
        # package into window
        hbox                  = wx.BoxSizer()
        x=hbox.Add(status, proportion=1, flag=wx.EXPAND)
        x=hbox.Add(abortButton, proportion=0, flag=wx.LEFT,border=5)
        x=hbox.Add(saveButton, proportion=0, flag=wx.LEFT,border=5)
        hbox2                 = wx.BoxSizer()
        x=hbox2.Add(command, proportion=1, flag=wx.EXPAND)
        x=hbox2.Add(sendButton, proportion=0, flag=wx.LEFT,border=5)
        vbox=wx.BoxSizer(wx.VERTICAL)
        x=vbox.Add(scrollText, proportion=1,
                   flag=wx.EXPAND|wx.LEFT|wx.LEFT|wx.LEFT, border=5)
        x=vbox.Add(hbox2, proportion=0, flag=wx.EXPAND|wx.ALL,border=5)
        x=vbox.Add(hbox, proportion=0, flag=wx.EXPAND|wx.ALL,border=5)
        bkg.SetSizer(vbox)
        #print "DEBUG MsgWinId, scrollTextId, commandId, statusId, abortButtonId",\
        #      TaskWin.MsgWinId, TaskWin.scrollTextId, TaskWin.commandId, \
        #      TaskWin.statusId, TaskWin.abortButtonId
        # end __init__

    def write(self, messages):
        """ Add messaged to message server

        message = message to display, a newline will be added
        """
        ################################################################
        #print "DEBUG write",messages
        if not messages:
            return
        TaskWin   = self.myTaskWin
        for message in messages:
            if type(message)==str:
                self.scrollText.AppendText(message+"\n")
            else:
                self.scrollText.AppendText(message[1]+"\n")
            self.number += 1
        # end write
        
