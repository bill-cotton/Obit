# OTWindow stuff to run in second (window) thread
import wx, thread, MsgWin

# Using new event binder
wx_EVT_MSGWIN = wx.NewEventType()
EVT_MSGWIN = wx.PyEventBinder(wx_EVT_MSGWIN, 1)

# Event type numbers
EVT_NEW     = 1  # Add MsgWin event
EVT_LABEL   = 2  # SetLabel event
EVT_BIND    = 3  # Bind Handler event
EVT_UPDATE  = 4  # Update event
EVT_MESS    = 5  # Message event

class MsgWinEvt(wx.PyEvent):
    def __init__(self):
        wx.PyEvent.__init__(self)
        self.SetEventType(wx_EVT_MSGWIN)


class HiddenCatcher(wx.Frame):
    """
        The "catcher" frame in the second thread.
        It is invisible.  It's only job is to receive
        Events from the main thread, and create and manipulate
        the appropriate windows.
    """
    def __init__(self):
        wx.Frame.__init__(self, None, -1, '')
        self.Bind(EVT_MSGWIN, self.EvtHandler)

    def EvtHandler(self, evt):
        """
        Event handler for HiddenCatcher class

        The event "type" is given by member evt_type:
        EVT_NEW     = 1  # Add MsgWin event
        EVT_LABEL   = 2  # SetLabel event
        EVT_BIND    = 3  # Bind Handler event
        EVT_UPDATE  = 4  # Update event
        EVT_MESS    = 5  # Message event
        The relevent additional parameters are given by type
        """
        if evt.evt_type==EVT_NEW:
            add_MsgWin(evt.evt_tw)
        elif evt.evt_type==EVT_LABEL:
            widget = self.FindWindowById(evt.evt_Id)
            widget.SetLabel(evt.evt_label)
        elif evt.evt_type==EVT_BIND:   # Better be a button
            widget = self.FindWindowById(evt.evt_Id)
            widget.Bind(wx.EVT_BUTTON, evt.evt_handler)
        elif evt.evt_type==EVT_UPDATE:
            widget = self.FindWindowById(evt.evt_Id)
            widget.Update()
        elif evt.evt_type==EVT_MESS:
            widget = self.FindWindowById(evt.evt_Id)
            for message in evt.evt_mess:
                if type(message)==str:
                    widget.AppendText(message+"\n")
                else:
                    widget.AppendText(message[1]+"\n")
       # end AddMsgWin

#---------------------------------------------------------------------------

class OTGUIApp(wx.App):
    """
    wxApp to run ObitTalk windows

    should run in same thread as GUI windows.    
    """
    def OnInit(self):
        catcher = HiddenCatcher()
        self.catcher = catcher
        return True

#---------------------------------------------------------------------------

def add_MsgWin(tw):
    """
    Make a task message window

    tw = TaskWindow Object
    """
    frame = MsgWin.MsgWin(tw)
    frame.Show(True)
    # end  add_MsgWin


