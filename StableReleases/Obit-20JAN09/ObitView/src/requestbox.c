/* $Id$  */
/* Request (of client) boxes  for ObitView */
/*-----------------------------------------------------------------------
*  Copyright (C) 2005-2008
*  Associated Universities, Inc. Washington DC, USA.
*  This program is free software; you can redistribute it and/or
*  modify it under the terms of the GNU General Public License as
*  published by the Free Software Foundation; either version 2 of
*  the License, or (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*-----------------------------------------------------------------------*/
#include <Xm/Xm.h> 
#include <Xm/DialogS.h> 
#include <Xm/MainW.h> 
#include <Xm/Form.h>
#include <Xm/PushB.h>
#include <Xm/Label.h>
#include <Xm/ToggleB.h>
#include <Xm/RowColumn.h>
#include <Xm/Separator.h>
#include <Xm/MessageB.h>
#include <Xm/TextF.h>
#include <Xm/Text.h>
#include "imagedisp.h"
#include "helpbox.h"
#include "messagebox.h"
#include "requestbox.h"
#include "drawbox.h"
#include "XMLRPCserver.h"
#include "ObitRPC.h"

/**
 *  \file requestbox.c
 * displays Request (of client) dialogs.
 */

/*--------------- file global data ----------------*/
/** is the option box active? */
int RequestBoxActive = 0;

/** Message for edit failed */
gchar *failedMessage = "Edit failed";

/* global structure for things to talk to each other */
typedef struct {
  Widget dialog, fieldlabel, field, OK;
  Widget radioCont, radioAbort, radioQuit, radioNoTv, radioView;
  gshort  xpos, ypos; /* location of the box */
  olong ReqType;  /* request type ObitRPCRequestType */
  olong Field;    /* request field number */
  olong curField; /* currently displayed field */
  olong nfield;   /* total number of fields in mosaic */
  olong cancel;   /* cancel requested if != 0 */
  XtIntervalId timerId;  /* X timer ID for timeout */
  olong setTO;    /* Timeout activated? */
} RequestBoxStuff;
RequestBoxStuff RBdia;

/**
 * Read  field value
 * \param w           widget 
 */
void ReadField (Widget w)
{
  gchar   *value=NULL;
  olong    temp;
  
  /* read value */
  value = XmTextGetString (w);
  if (!value) /* error */
    {MessageShow ("Error reading requested field value");
    return;}
  if (!sscanf (value, "%d", &temp))
    { /* error */
      MessageShow ("Error reading requested field value");
      if (value) XtFree(value); value = NULL;
      return;}
  if (value) XtFree(value);value = NULL;

  /* OK, save */
  temp = MAX (1, MIN (temp, RBdia.nfield));
  RBdia.Field = temp;
  
} /* end ReadField */

/**
 * Callback for Edit Request type
 * Editing requests passed to client
 * \param w      widget activated
 * \param clientData  client data = which  
 *        =0=>Continue, 1=>Abort, 2=>Quit, 3=>NoTV, 4=>View
 * \param callData    call data = state  button state
 */
void EditRequestCB (Widget w,  XtPointer clientData, XtPointer callData)
{
  int which = (int)clientData;
  XmToggleButtonCallbackStruct *state = (XmToggleButtonCallbackStruct*)callData;

  /* Cancel any timeout */
  if (RBdia.setTO) XtRemoveTimeOut(RBdia.timerId);
  RBdia.setTO  = 0;  /* Timeout deactivated */
 
  if (state->set) {
    /*fprintf (stderr,"DEBUG which %d\n",which);*/
    if (which==0) RBdia.ReqType = OBIT_RPC_Request_Continue;
    else if (which==1) RBdia.ReqType = OBIT_RPC_Request_Abort;
    else if (which==2) RBdia.ReqType = OBIT_RPC_Request_Quit;
    else if (which==3) RBdia.ReqType = OBIT_RPC_Request_NoTV;
    else if (which==4) RBdia.ReqType = OBIT_RPC_Request_Field;
  } else return;

} /* end EditRequestCB */

/* button callbacks */
/**
 * Callback for Clear windows button hit 
 * \param w           widget activated
 * \param clientData  client data
 * \param callData    call data
 */
static void ReqClearButCB (Widget w, XtPointer clientData, XtPointer callData)
{
  ClearBox();

  /* Cancel any timeout */
  if (RBdia.setTO) XtRemoveTimeOut(RBdia.timerId);
  RBdia.setTO  = 0;  /* Timeout deactivated */
  
} /* end ReqClearButCB */

/**
 * Callback for Edit Windows button hit 
 * \param w           widget activated
 * \param clientData  client data
 * \param callData    call data
 */
static void ReqEditButCB (Widget w, XtPointer clientData, XtPointer callData)
{
  /* Cancel any timeout */
  if (RBdia.setTO) XtRemoveTimeOut(RBdia.timerId);
  RBdia.setTO  = 0;  /* Timeout deactivated */
  
  /* make box disappear but still exist */
  XtPopdown(RBdia.dialog);
  EditBox ();
  
} /* end ReqEditButCB */

/**
 * Callback for OK button hit 
 * \param w           widget activated
 * \param clientData  client data
 * \param callData    call data
 */
static void ReqOKButCB (Widget w, XtPointer clientData, XtPointer callData)
{
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  
  ReadField (RBdia.field);
  
  /* make it disappear but still exist */
  XtPopdown(RBdia.dialog);

  /* set request on requestList */
  if (requestList==NULL) requestList = newObitInfoList();
  ObitInfoListPut (requestList, "Request", OBIT_long, dim,  
		   (gpointer)&RBdia.ReqType, err);
  if (RBdia.ReqType==OBIT_RPC_Request_Field)
  ObitInfoListPut (requestList, "Field", OBIT_long, dim,  
		   (gpointer)&RBdia.Field, err);
  if (err->error) ObitErrLog(err);  /* any errors */
  
  /* Cancel any timeout */
  if (RBdia.setTO) XtRemoveTimeOut(RBdia.timerId);
  RBdia.setTO  = 0;  /* Timeout deactivated */
  
  /* Release pending XML reply */
  XMLRPCRelease(NULL);
} /* end ReqOKButCB */

/**
 * Callback for Cancel button hit 
 * \param w           widget activated
 * \param clientData  client data
 * \param callData    call data
 */
static void ReqCancelButCB (Widget w, XtPointer clientData, XtPointer callData)
{
  /* make it disappear but still exist */
  XtPopdown(RBdia.dialog);

  /* Want to cancel editing */
  RBdia.cancel = 1;

  /* Cancel any timeout */
  if (RBdia.setTO) XtRemoveTimeOut(RBdia.timerId);
  RBdia.setTO  = 0;  /* Timeout deactivated */
  
  /* Release pending XML reply */
  XMLRPCRelease(failedMessage);
} /* end ReqCancelButCB */

/**
 * Callback for Timeout - activate OK button
 * \param clientData  client data
 * \param timer_ID    X Timer ID
 */
static void TimeOutCB (XtPointer clientData, XtIntervalId* timer_ID)
{
  ReqOKButCB (RBdia.OK, clientData, NULL);
} /* end TimeOutCB */

/**
 * Dialog box for edit Requests
 * When the OK button is hit the request is entered in global ObitInfoList requestList
 * and any pending RPC reply is released.
 * When the Cancel button is hit, any pending RPC reply is released.
 */
void EditRequestBox ()
{
  Widget form=NULL, label1=NULL, sep=NULL;
  Widget radio=NULL;
  Widget OKbutton=NULL, CancelButton=NULL, HelpButton=NULL;
  Widget ClearButton=NULL, EditButton=NULL;
  XmString     label=NULL, Field=NULL, WierdString=NULL;
  XmString     Continue=NULL, Abort=NULL, Quit=NULL, View=NULL, NoTV=NULL;
  gchar        valuestr[61];
  unsigned long interval;
#define REQUESTBOX_WIDTH 130
#define REQUESTBOX_HEIGHT 280
  
  /* field info from global */
  RBdia.curField = image[CurImag].Field;
  RBdia.nfield   = image[CurImag].NField;

  RBdia.cancel = 0;  /* Unset cancel request */
  RBdia.setTO  = 0;  /* No timeout activated */
  RBdia.timerId= 0;  /* TimerID */
    
  /* don't make another one */
  if (RequestBoxActive) {
    if (XtIsRealized (RBdia.dialog))
      XMapRaised (XtDisplay(RBdia.dialog), XtWindow(RBdia.dialog));
    
    /* bring it back where we can see it */
    XtPopup(RBdia.dialog, XtGrabNonexclusive);
    
    /* put it some place reasonable */
    /*  where is parent? */
    XtVaGetValues (Display_shell, XmNx, &RBdia.xpos, XmNy, &RBdia.ypos,  NULL);
    RBdia.ypos += 170;
    if (RBdia.xpos<0) RBdia.xpos += 100;
    XMoveWindow (XtDisplay(Display_shell), XtWindow(RBdia.dialog), 
		 RBdia.xpos, RBdia.ypos);

    /* Set request button to "Continue" */
    if (RBdia.ReqType != OBIT_RPC_Request_Continue) {
      XtVaSetValues (RBdia.radioAbort, XmNset, False, NULL);
      XtVaSetValues (RBdia.radioQuit, XmNset, False, NULL); 
      XtVaSetValues (RBdia.radioNoTv, XmNset, False, NULL); 
      XtVaSetValues (RBdia.radioView, XmNset, False, NULL); 
      XtVaSetValues (RBdia.radioCont, XmNset, True, NULL);
      RBdia.ReqType = OBIT_RPC_Request_Continue;
    }

    /* set value of requested field  */
    g_snprintf (valuestr, 60, "%d", RBdia.curField);
    XtVaSetValues(RBdia.field,
		  XmNvalue,  valuestr,
		  NULL);
  
    /* set field number label field  */
    g_snprintf (valuestr, 60, "Request field of %d", RBdia.nfield);
    WierdString = XmStringCreateSimple (valuestr);
    XtVaSetValues(RBdia.fieldlabel,
		  XmNlabelString,   WierdString,
		  NULL);
    if (WierdString) XmStringFree(WierdString); WierdString = NULL;
  
    /* Timeout enabled? */
    if (ERtime_out>0.0) {
      interval = 1000*MAX(ERtime_out,5.0);
      RBdia.timerId =
	XtAppAddTimeOut (myApp, interval, (XtTimerCallbackProc)TimeOutCB, NULL);
      RBdia.setTO  = 1;  /* Timeout activated */
  }
    return;
  } /* end of update dialog */

  label    = XmStringCreateSimple ("Requests panel");
  Continue = XmStringCreateSimple ("continue");
  Abort    = XmStringCreateSimple ("abort");
  Quit     = XmStringCreateSimple ("quit");
  View     = XmStringCreateSimple ("view field");
  NoTV     = XmStringCreateSimple ("Turn Off TV");

  /* mark as active */
  RequestBoxActive = 1;
  
  RBdia.dialog = XtVaCreatePopupShell ("RequestBox", xmDialogShellWidgetClass, 
				       Display_shell, 
				       XmNautoUnmanage, False,
				       XmNwidth,     REQUESTBOX_WIDTH,
				       XmNheight,    REQUESTBOX_HEIGHT,
				       XmNdeleteResponse, XmDESTROY,
				       NULL);
  
  /* make Form widget to stick things on */
  form = XtVaCreateManagedWidget ("RequestForm", xmFormWidgetClass,
				  RBdia.dialog,
				  XmNautoUnmanage, False,
				  XmNwidth,     REQUESTBOX_WIDTH,
				  XmNheight,    REQUESTBOX_HEIGHT,
				  XmNx,           0,
				  XmNy,           0,
				  NULL);
  
  /* info label widgets */
  label1 = XtVaCreateManagedWidget ("Label1", xmLabelWidgetClass, 
				    form, 
				    XmNwidth,           REQUESTBOX_WIDTH,
				    XmNlabelString,   label,
				    XmNtopAttachment, XmATTACH_FORM,
				    XmNleftAttachment,  XmATTACH_FORM,
				    NULL);
  
  /* Request radio buttons - aint Motif wonderful? */
  radio = 
    XtVaCreateManagedWidget("RadioBox", xmRowColumnWidgetClass,
			    form, 
			    XmNorientation, XmVERTICAL,
			    XmNradioBehavior, True,
			    XmNradioAlwaysOne, True,
			    XmNtopAttachment, XmATTACH_WIDGET,
			    XmNtopWidget,     label1,
			    XmNleftAttachment,  XmATTACH_FORM,
			    NULL);

  
 RBdia.radioCont = 
   XtVaCreateManagedWidget( "Continue",
			    xmToggleButtonWidgetClass, radio,
			    XmNset, True,
			    NULL);
 XtAddCallback (RBdia.radioCont, XmNvalueChangedCallback, EditRequestCB, (XtPointer)0);

 RBdia.radioAbort = 
   XtVaCreateManagedWidget("Abort", xmToggleButtonWidgetClass, radio, NULL);
 XtAddCallback (RBdia.radioAbort, XmNvalueChangedCallback, EditRequestCB, (XtPointer)1);

 RBdia.radioQuit = 
  XtVaCreateManagedWidget("Quit operation", xmToggleButtonWidgetClass, radio, NULL);
 XtAddCallback (RBdia.radioQuit, XmNvalueChangedCallback, EditRequestCB, (XtPointer)2);

 RBdia.radioNoTv = 
  XtVaCreateManagedWidget("Turn off TV", xmToggleButtonWidgetClass, radio, NULL);
 XtAddCallback (RBdia.radioNoTv, XmNvalueChangedCallback, EditRequestCB, (XtPointer)3);

 RBdia.radioView = 
  XtVaCreateManagedWidget("View field:", xmToggleButtonWidgetClass, radio, NULL);
 XtAddCallback (RBdia.radioView, XmNvalueChangedCallback, EditRequestCB, (XtPointer)4);

  /* requested field */
  g_snprintf (valuestr, 60, "Request field of %d", RBdia.nfield);
  Field = XmStringCreateSimple (valuestr);
  RBdia.fieldlabel = XtVaCreateManagedWidget ("FieldRequest", xmLabelWidgetClass,
					    form,
					    XmNwidth,           REQUESTBOX_WIDTH,
					    XmNlabelString,   Field,
					    XmNtopAttachment, XmATTACH_WIDGET,
					    XmNtopWidget,     radio,
					    XmNleftAttachment,  XmATTACH_FORM,
					    NULL);
  
  RBdia.Field = RBdia.curField;
  g_snprintf (valuestr, 60, "%d", RBdia.Field);
  RBdia.field = XtVaCreateManagedWidget ("Field", xmTextFieldWidgetClass,
				       form, 
				       XmNwidth,         REQUESTBOX_WIDTH,
				       XmNvalue,         valuestr,
				       XmNtopAttachment, XmATTACH_WIDGET,
				       XmNtopWidget,     RBdia.fieldlabel,
				       XmNleftAttachment,  XmATTACH_FORM,
				       NULL);
  
  /* separator */
  sep = XtVaCreateManagedWidget ("sep", xmSeparatorWidgetClass,
				 form, 
				 XmNwidth,           REQUESTBOX_WIDTH,
				 XmNtopAttachment,   XmATTACH_WIDGET,
				 XmNtopWidget,        RBdia.field,
				 XmNleftAttachment,  XmATTACH_FORM,
				 NULL);
  
  /* Edit button */
  EditButton = XtVaCreateManagedWidget ("Edit", xmPushButtonWidgetClass, 
					form, 
					XmNtopAttachment, XmATTACH_WIDGET,
					XmNtopWidget,     sep,
					XmNleftAttachment,  XmATTACH_FORM,
					NULL);
  XtAddCallback (EditButton, XmNactivateCallback, ReqEditButCB, NULL);
  
  /* Clear button */
  ClearButton = XtVaCreateManagedWidget ("Clear", xmPushButtonWidgetClass, 
					  form, 
					 XmNtopAttachment, XmATTACH_WIDGET,
					 XmNtopWidget,     RBdia.field,
					 XmNleftAttachment, XmATTACH_WIDGET,
					 XmNleftWidget,    EditButton,
					 NULL);
  XtAddCallback (ClearButton, XmNactivateCallback, ReqClearButCB, NULL);
  /* OK button */
  OKbutton = XtVaCreateManagedWidget (" OK ", xmPushButtonWidgetClass, 
				      form, 
				      XmNbottomAttachment, XmATTACH_FORM,
				      XmNleftAttachment,  XmATTACH_FORM,
				      NULL);
  XtAddCallback (OKbutton, XmNactivateCallback, ReqOKButCB, NULL);
  RBdia.OK = OKbutton;

  /* Cancel button */
  CancelButton = XtVaCreateManagedWidget ("Cancel", xmPushButtonWidgetClass, 
					  form, 
					  XmNbottomAttachment, XmATTACH_FORM,
					  XmNleftAttachment, XmATTACH_WIDGET,
					  XmNleftWidget,     OKbutton,
					  NULL);
  XtAddCallback (CancelButton, XmNactivateCallback, ReqCancelButCB, NULL);
  
  /* Help button */
  HelpButton = XtVaCreateManagedWidget ("Help", xmPushButtonWidgetClass, 
					form, 
					XmNbottomAttachment, XmATTACH_FORM,
					XmNleftAttachment, XmATTACH_WIDGET,
					XmNleftWidget,     CancelButton,
					NULL);
  XtAddCallback (HelpButton, XmNactivateCallback,  HelpBoxTopicCB, 
		 (XtPointer)"Window Editing");
  
  if (label)    XmStringFree(label);    label  = NULL;
  if (Field)    XmStringFree(Field);    Field  = NULL;
  if (Continue) XmStringFree(Continue); Continue = NULL;
  if (Abort   ) XmStringFree(Abort);    Abort= NULL;
  if (Quit)     XmStringFree(Quit);     Quit = NULL;
  if (View)     XmStringFree(View);     View = NULL;
  if (NoTV)     XmStringFree(NoTV);     NoTV = NULL;
  
  /* set it up */
  XtManageChild (RBdia.dialog);
  
  /* put it some place reasonable */
  /*  where is parent? */
  XtVaGetValues (Display_shell,
		 XmNx, &RBdia.xpos,
		 XmNy, &RBdia.ypos,
		 NULL);
  RBdia.ypos += 170;
  RBdia.xpos += 100;
  if (RBdia.xpos<0) RBdia.xpos = 0;
  XMoveWindow (XtDisplay(Display_shell), XtWindow(RBdia.dialog), 
	       RBdia.xpos, RBdia.ypos);

  /* Timeout enabled? */
  if (ERtime_out>0.0) {
    interval = 1000*MAX(ERtime_out,5.0);
    RBdia.timerId =
      XtAppAddTimeOut (myApp, interval, (XtTimerCallbackProc)TimeOutCB, NULL);
      RBdia.setTO  = 1;  /* Timeout activated */
  }
} /* end RequestBox */

