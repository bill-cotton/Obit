""" OTWindow allows running wxPython widgets

Widgets using wxPython must all be created and run in the same thread.
This class creates a wxPython App in a separate thread and allows starting
new widgets in this same thread.
   New widgets can be created using the functions
   newMsgWin(tw) create message windiow and execute TaskWindow tw
"""
# $Id: OTWindow.py,v 1.2 2007/02/28 19:17:48 bcotton Exp $
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
class OTWindow:
    def start(self):
        """ start the GUI thread
        """
        import  thread
        thread.start_new_thread(self.run, ())

    def run(self):
        """
            Note that OTWindow2 is first imported ***here***.
            This is the second thread.
            OTWindow2  imports wxPython, if we imported it at
            the module level instead of in this function,
            the import would occur in the main thread and
            wxPython would not run correctly in the second thread.
            The wxPython GUI MainLoop is run here (i.e. no return)
        """
        ################################################################
        try:
            import  OTWindow2
            self.app = OTWindow2.OTGUIApp()
            self.started = True
            self.app.MainLoop()
        except TypeError:
            self.app = None
        except Exception, e:
            self.app = None
            #print "DEBUG: oh bugger untrapped exception in OTWindow.run"
            #print e

    def add_MsgWin(self, tw):
        """
        New Task message widget
        
        Send an event to the catcher window in the
        other thread and tell it to create a MsgWin window.
        tw = TaskWindow of task to be run
        """
        ################################################################
        import  OTWindow2, MsgWin

        if self.app:
            evt = OTWindow2.MsgWinEvt()
            evt.evt_type = OTWindow2.EVT_NEW   # Set event type
            evt.evt_tw   = tw                  # add task window
            self.app.catcher.AddPendingEvent(evt);
        else:
            OTWindow2.add_MsgWin(tw)
        # end add_MsgWin

    def SetLabel(self, Id, label):
        """
        Set Widget label
        
        Send an event to the catcher window in the other thread
        and tell it to Setlabel on window Id to label
        Id    = widget Id
        label = New text to label
        """
        ################################################################
        import  OTWindow2

        evt           = OTWindow2.MsgWinEvt()
        evt.evt_type  = OTWindow2.EVT_LABEL   # Set event type
        evt.evt_Id    = Id                    # Set widget Id
        evt.evt_label = label                 # add new label text
        self.app.catcher.AddPendingEvent(evt);
        # end SetLabel

    def Bind(self, Id, handler):
        """
        Set Button event handler 
        
        Send an event to the catcher window in the other thread
        and tell it to rebind the event handler on button Id
        Id      = widget Id
        handler = new event handler
        """
        ################################################################
        import  OTWindow2

        evt             = OTWindow2.MsgWinEvt()
        evt.evt_type    = OTWindow2.EVT_BIND  # Set event type
        evt.evt_Id      = Id                  # Set widget Id
        evt.evt_handler = handler             # add task window
        self.app.catcher.AddPendingEvent(evt);
        # end Bind

    def Update(self, Id):
        """
        Update Widget Id
        
        Send an event to the catcher window in the other thread
        and tell it to refresh the display of widget Id
        Id      = widget Id
        """
        ################################################################
        import  OTWindow2

        evt          = OTWindow2.MsgWinEvt()
        evt.evt_type = OTWindow2.EVT_UPDATE  # Set event type
        evt.evt_Id   = Id                    # Set widget Id
        self.app.catcher.AddPendingEvent(evt);
        # end Update

    def Message(self, Id, message):
        """
        Write messages in TextCtrl
        
        Send an event to the catcher window in the other thread
        and tell it to append message(s) in widget Id
        Id      = widget Id (a TextCtrl)
        message = either a single string or an array of strings
        """
        ################################################################
        import  OTWindow2

        evt          = OTWindow2.MsgWinEvt()
        evt.evt_type = OTWindow2.EVT_MESS # Set event type
        evt.evt_Id   = Id                    # Set widget Id
        evt.evt_mess = message               # add task message(s)
        self.app.catcher.AddPendingEvent(evt);
        # end Message

    # end class OTWindow

# Startup wxPython windowing
gui = OTWindow()
gui.started = False
gui.start()

# Externally callable routine to create a MsgWin (task message window)
def newMsgWin(tw):
    """
    New task message window
    
    Create a new task message window, run the task displaying messages
    and handling communications
    tw = TaskWindow for task to be executed.
    """
    ################################################################
    # Be sure gui thread started
    import  time
    while (not gui.started):
        time.sleep(0.2)
    gui.add_MsgWin(tw)
    # end newMsgWin

def CallSetLabel (Id, label):
    """
    Set label on widget Id

    Id    = widget Id
    label = New text to label
    """
    gui.SetLabel (Id, label)
    # end CallSetLabel

def CallBind (Id, handler):
    """
    Set Button event handler 

    Id      = widget Id
    handler = new event handler
    """
    gui.Bind (Id, handler)
    # end CallBind

def CallUpdate(Id):
    """
    Update Widget Id
    
    Id    = widget Id
    """
    gui.Update (Id)
    # end CallUpdate

def CallMessage (Id, message):
    """
    Set label on widget Id
    
    Id      = widget Id
    message = either a single string or an array of strings
    """
    gui.Message (Id, message)
    # end CallMessage 

   
