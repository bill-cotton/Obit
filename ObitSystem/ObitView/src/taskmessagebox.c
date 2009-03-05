/* $Id: $  */
/* Task Message box  */
/* adopted from "Power programming Motif" by E. F. Johnson and
   K. Reichard, 1993, MIS Press, New York */
/*-----------------------------------------------------------------------
*  Copyright (C) 2009
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
*  Correspondence concerning ObitView should be addressed as follows:
*         Internet email: bcotton@nrao.edu.
*         Postal address: William Cotton
*                         National Radio Astronomy Observatory
*                         520 Edgemont Road
*                         Charlottesville, VA 22903-2475 USA
*-----------------------------------------------------------------------*/

#include <Xm/Xm.h> 
#include <Xm/Form.h> 
#include <Xm/Label.h> 
#include <Xm/List.h> 
#include <Xm/PushB.h> 
#include <Xm/RowColumn.h> 
#include <Xm/Separator.h> 
#include <Xm/Text.h> 
#include <Xm/TextF.h> 
#include <Xm/DialogS.h> 
#include <Xm/FileSB.h> 
#include <X11/IntrinsicP.h>
#include <glib.h>
#include <stdio.h>
#include <unistd.h>
#include "taskmessagebox.h"
#include "taskmessage.h"

/**
 *  \file taskmessagepbox.c
 * displays task message dialog.
 */

/*--------------- Private data ----------------*/
gchar *SaveAsDir=NULL;

/*---------------Private function prototypes----------------*/
void InitTaskMessageText(void);
void TaskMessageBoxSelectCB(Widget parent, XtPointer clientData, XtPointer callData);
static void unmanage_dialogCB(Widget widget, 
			      XtPointer client_data, 
			      XtPointer call_data);
static void SaveAsCB(Widget widget, 
		     XtPointer client_data, 
		     XtPointer call_data);
static void UserInputCB(Widget widget, 
			XtPointer client_data, 
			XtPointer call_data);
static void AbortCB(Widget widget, 
		    XtPointer client_data, 
		    XtPointer call_data);
void TaskMessageBoxShow(TaskMessageBox* tmb);
/** Delete task message window */
void TaskMessageDelete(glong taskID, glong *status);

/*------------------- Public functions   --------------------*/
/**
 * Create task message box 
 * \param parent           parent widget
 * \param taskname         Task name
 * \param taskID           TaskID
 * \return pointer to newly created dialog widget object structure
 */
TaskMessageBox* TaskMessageBoxMake (Widget parent, gchar *taskname, 
				    glong taskID)
{
  TaskMessageBox *tmb=NULL;
  Arg         args[20];
  Cardinal    n;
  Widget      row, dismiss, saveas, abort;
  XFontStruct *XFont;
  GC          gc;
  int         chei, cwid;
  Dimension   textHei, textWid;
  XmString     label = NULL;

  /* Create Structure */
  tmb = g_malloc0(sizeof(TaskMessageBox));
  tmb->task_dialog = NULL;
  tmb->task_text   = NULL;
  tmb->user_input  = NULL;
  tmb->status      = NULL;
  tmb->taskID      = taskID;
  tmb->doAbort     = FALSE;
  tmb->doUserInput = FALSE;
  
  /* see how big characters are to make boxes to size */
  /* create graphics context for box */
  gc = XCreateGC (XtDisplay (parent), 
		  DefaultRootWindow (XtDisplay(parent)),
		  0, NULL);
  
  /* how big are the characters */
  XFont = XQueryFont(XtDisplay (parent), XGContextFromGC(gc));
  chei = XFont->ascent + XFont->descent + 2;
  cwid = XFont->max_bounds.width;
  
  /* text ~120 char x 30 lines */
  textWid = 120*cwid;
  textHei =  30*chei;
  
  if (gc) XFreeGC(XtDisplay(parent), gc); gc = NULL;
  
  n = 0;
  XtSetArg(args[n], XmNallowResize, True); n++;
  XtSetArg(args[n], XmNtitle, taskname); n++;
  XtSetArg(args[n], XmNautoUnmanage, False); n++;
  tmb->task_dialog  = XmCreateFormDialog(parent, "Messages", args, n);
  
  /* create button area at bottom */
  /* Note: the stuff at the bottom needs to be put in first due to the
     adjustlast policy of the RowColumn widget */
  /* Abort button */
  abort = XtVaCreateManagedWidget("Abort",
				  xmPushButtonWidgetClass, tmb->task_dialog,
				  /*XmNtopAttachment,   XmATTACH_WIDGET,
				    XmNtopWidget,       sep,*/
				  /*XmNleftAttachment,    XmATTACH_FORM,*/
				  XmNrightAttachment,   XmATTACH_FORM,
				  XmNbottomAttachment,  XmATTACH_FORM,
				  NULL);
  /* Add callback */
  XtAddCallback(abort, XmNactivateCallback, 
		(XtCallbackProc)AbortCB, (XtPointer)tmb);

  /* SaveAs button */
  saveas = XtVaCreateManagedWidget("SaveAs",
				   xmPushButtonWidgetClass, tmb->task_dialog,
				   /*XmNtopAttachment,   XmATTACH_WIDGET,
				     XmNtopWidget,       sep,*/
				   /*XmNleftAttachment,    XmATTACH_FORM,*/
				   XmNrightAttachment,   XmATTACH_WIDGET,
				   XmNrightWidget,       abort,
				   XmNbottomAttachment,  XmATTACH_FORM,
				   NULL);
  /* Add callback */
  XtAddCallback(saveas, XmNactivateCallback, 
		(XtCallbackProc)SaveAsCB, (XtPointer)tmb);

  /* Dismiss button */
  dismiss = XtVaCreateManagedWidget("Dismiss",
				    xmPushButtonWidgetClass, tmb->task_dialog,
				    /*XmNtopAttachment,   XmATTACH_WIDGET,
				      XmNtopWidget,       sep
				      XmNleftAttachment,  XmATTACH_FORM,,*/
				    XmNrightAttachment,  XmATTACH_WIDGET,
				    XmNrightWidget,       saveas,
				    XmNbottomAttachment, XmATTACH_FORM,
				    NULL);
  /* Add callback */
  XtAddCallback(dismiss, XmNactivateCallback, 
		(XtCallbackProc)unmanage_dialogCB, (XtPointer)tmb);

  /* Add Status label */
  label  = XmStringCreateSimple ("Task Running");
  tmb->status = XtVaCreateManagedWidget ("statusLabel", xmLabelWidgetClass, 
					 tmb->task_dialog, 
					 XmNwidth,         20,
					 XmNleftAttachment,    XmATTACH_FORM,
					 XmNrightAttachment,   XmATTACH_WIDGET,
					 XmNrightWidget,       dismiss,
					 XmNbottomAttachment,  XmATTACH_FORM,
					 XmNlabelString,   label,
					 NULL);

  /* User input text box */
  tmb->user_input = XtVaCreateManagedWidget("userinput",xmTextFieldWidgetClass,
					    tmb->task_dialog,
					    XmNwidth,           textWid,
					    XmNbottomAttachment,XmATTACH_WIDGET,
					    XmNbottomWidget,    dismiss,
					    XmNleftAttachment,  XmATTACH_FORM,
					    XmNrightAttachment, XmATTACH_FORM,
					    NULL);
  /* Add callback */
  XtAddCallback(tmb->user_input, XmNactivateCallback, 
		(XtCallbackProc)UserInputCB, (XtPointer)tmb);
   
  /* Form for text */
  row = XtVaCreateWidget("row",
			 xmFormWidgetClass, tmb->task_dialog,
			 XmNorientation, XmHORIZONTAL,
			 XmNwidth, textWid+20,
			 XmNresizeHeight, True,
			 XmNresizeWidth,  True,
			 NULL);
  
  /* Create text Widget to display messages */
  n = 0;
  XtSetArg(args[n], XmNeditMode, XmMULTI_LINE_EDIT); n++;
  XtSetArg(args[n], XmNbottomAttachment, XmATTACH_FORM); n++;
  XtSetArg(args[n], XmNtopAttachment,    XmATTACH_FORM); n++;
  XtSetArg(args[n], XmNrightAttachment,  XmATTACH_FORM); n++;
  XtSetArg(args[n], XmNleftAttachment,   XmATTACH_FORM); n++;
  XtSetArg(args[n], XmNheight, textHei); n++; 
  XtSetArg(args[n], XmNwidth,  textWid); n++;
  XtSetArg(args[n], XmNresizeHeight, True); n++;
  XtSetArg(args[n], XmNresizeWidth,  True); n++;
  XtSetArg(args[n], XmNeditable, False); n++;
  XtSetArg(args[n], XmNautoShowCursorPosition, False); n++;
  tmb->task_text = XmCreateScrolledText (row, "task_text", args, n);
  
  /* set up scrolled area attachments */
  XtVaSetValues (row,
		 XmNtopAttachment,    XmATTACH_FORM,
		 XmNleftAttachment,   XmATTACH_FORM,
		 XmNrightAttachment,  XmATTACH_FORM,
		 XmNbottomAttachment, XmATTACH_WIDGET,
		 XmNbottomWidget,     tmb->user_input,
		 NULL);
  
  XtManageChild(tmb->task_text);
  XtManageChild(row);

  /* show the box */
  TaskMessageBoxShow(tmb);

  /* Cleanup */
  if (label)   XmStringFree(label); label   = NULL;

  return tmb;
} /* end  TaskMessageBoxCreate */

/**
 * Display message  
 * Leave display at the bottom if it is currently there
 * \param tmb      TaskMessageBox object
 * \param message  list of messages
 */
void TaskMessageBoxDisplay (TaskMessageBox* tmb, gchar *message)
{
  Arg  args[20];
  Cardinal n;
  int maximum=0, minimum=0, size=0, value=0;
  Widget vScroll=NULL;
  XmTextPosition bottom;
  gboolean showBottom=FALSE;

  if (tmb->task_text == NULL) return; /* sanity check */
  if (message == NULL) return;
  
  /* Where are we now? Look at vertical scrollbar*/
  n = 0;
  XtSetArg (args[n], XmNverticalScrollBar, &vScroll); n++;
  XtGetValues (XtParent(tmb->task_text), args, n);
  if (vScroll) {
    n = 0;
    XtSetArg (args[n], XmNmaximum,    &maximum); n++;
    XtSetArg (args[n], XmNminimum,    &minimum); n++;
    XtSetArg (args[n], XmNvalue,      &value); n++;
    XtSetArg (args[n], XmNsliderSize, &size); n++;
    XtGetValues (vScroll, args, n);
    showBottom  = (value+size)>=maximum;
  } else {
     showBottom = TRUE;
  }

  /* insert text at bottom */
  bottom  = XmTextGetLastPosition (tmb->task_text);
  XmTextInsert(tmb->task_text, bottom, message);

  /* track bottom if there */
  if (showBottom) {
    bottom  = XmTextGetLastPosition (tmb->task_text);
    XmTextShowPosition (tmb->task_text, bottom);
  }

} /* end  TaskMessageBoxDisplay */

/**
 * Setup for  user input - response from UserInputCB
 * \param tmb      TaskMessageBox object
 */
void TaskMessageBoxUserRequest (TaskMessageBox* tmb)
{
  char     *blank="";

  if (tmb->user_input == NULL) return; /* sanity check */

  /* Blank */
  XmTextSetString (tmb->user_input, blank);

  /* Reset status */
  TaskMessageBoxStatus(tmb, "Awaiting user input");

  /* Now awaiting user input */
  tmb->doUserInput = TRUE;
} /* end  TaskMessageBoxUserRequest */

/**
 * Set task status
 * \param tmb      TaskMessageBox object
 * \param status   New status string
 */
void TaskMessageBoxStatus (TaskMessageBox* tmb, gchar *status)
{
  XmString  wierdstring = NULL;
  Arg       wargs[5]; 

  if (tmb->status == NULL) return; /* sanity check */
  if (status == NULL) return;
  
  /* Convert to X string and update */
  wierdstring = XmStringCreateSimple (status);
  XtSetArg (wargs[0], XmNlabelString, wierdstring);
  XtSetValues (tmb->status, wargs, 1);
  if (wierdstring) XmStringFree(wierdstring); wierdstring = NULL;
} /* end  TaskMessageBoxStatus */

/**
 * Delete object and dialog
 * \param tmb      TaskMessageBox object
 */
void TaskMessageBoxDelete (TaskMessageBox* tmb)
{
  if (tmb->task_dialog == NULL) return; /* sanity check */

  /* If there is a pending request for user input, clear it */
  if (tmb->doUserInput) {
    tmb->doUserInput = FALSE;  /* No longer waiting */
    tmb->doAbort     = TRUE;   /* Try to kill it */
    /* Return response */
    TaskMessageUserResponseReply (tmb->taskID, "q");
    /* Give it a chance to happen */
    usleep(500000);  /* 500 msec */
  }

  /* Kill da wabbit */
  if (XtIsManaged(tmb->task_dialog)) XtUnmanageChild(tmb->task_dialog);

  /* free object */
  g_free(tmb);

} /* end  TaskMessageBoxDelete */

/*---------------Private functions -----------------------*/
/**
 * Callback to unmanage dialog
 * \param widget      widget activated
 * \param clientData  client data
 * \param callData    call data, pointer to TaskMessageBox
 */
static void unmanage_dialogCB(Widget widget, 
			      XtPointer client_data, 
			      XtPointer call_data) {
  TaskMessageBox *tmb = (TaskMessageBox*)client_data;
  glong status;

  /* Use TaskMessage for proper cleanup */
  TaskMessageDelete(tmb->taskID, &status);
} /* end unmanage_dialogCB */

/**
 * Save messagebox text to file
 * \param filebox     widget activated
 * \param clientData  client data (TaskMessageBox*)
 * \param callData    call data
 */
static void SaveAsFileOKCB (Widget filebox, 
			    XtPointer client_data, 
			    XtPointer call_data)
{
  TaskMessageBox *tmb=(TaskMessageBox*)client_data;
  gchar *filename, *directory;
  XmFileSelectionBoxCallbackStruct *cbs;
  FILE *file;                /* file pointer for write */
  gchar *text;
  
  cbs = (XmFileSelectionBoxCallbackStruct *) call_data;
  
  /* get file name */
  if (!XmStringGetLtoR (cbs->value, XmSTRING_DEFAULT_CHARSET, &filename))
    return; /* error */

  /* Save contents to filename */
  file = fopen ((const char*)filename,"at");
  if (!file) { /* failed? */
    return; } /* end of open failed */

  /* Write file */
  text = XmTextGetString(tmb->task_text);
  if (text) {
    fprintf (file, "%s\n", text);
    XtFree(text);
  }

  /* Close */
  fclose(file);

  /* save directory name */
  if (!XmStringGetLtoR (cbs->dir, XmSTRING_DEFAULT_CHARSET, &directory))
    return; /* error */
  if (SaveAsDir) g_free(SaveAsDir);
  SaveAsDir = g_strdup(directory);
  if (filename) XtFree(filename); filename = NULL; 
  if (directory) XtFree(directory); directory = NULL;
  
  /* Shazam disappear and die */
#ifndef KEEP_FILE_DIALOG  /* option to let it hang around */
  XtUnmanageChild (filebox);
  XtPopdown (XtParent (filebox));
  XtDestroyWidget(filebox); 
#endif
  
} /* end SaveAsFileOKCB */

/**
 * Cancel save messagebox text to file
 */
static void SaveAsFileCancelCB (Widget filebox, 
				XtPointer client_data, 
				XtPointer call_data)
     /* cancel file selection dialog box */
{
  /* Shazam disappear and die */
  XtUnmanageChild (filebox); 
  XtPopdown (XtParent (filebox)); 
  XtDestroyWidget(filebox); 
} /* end SaveAsFileCancelCB */

/**
 * Callback to save contents in a file
 * \param widget      widget activated
 * \param clientData  client data
 * \param callData    call data, pointer to TaskMessageBox
 */
static void SaveAsCB(Widget widget, 
		     XtPointer client_data, 
		     XtPointer call_data) 
{
  Widget       filebox, help_button;
  XmString     wierdstring = NULL;
  Arg          wargs[5]; 

  /* Create file selection box */
  filebox = (Widget) XmCreateFileSelectionDialog (widget, "file_open", NULL, 0);
  XtAddCallback (filebox, XmNokCallback, SaveAsFileOKCB, client_data);
  XtAddCallback (filebox, XmNcancelCallback, SaveAsFileCancelCB, client_data);
  /* Gray out help button */
  help_button = XmFileSelectionBoxGetChild(filebox, XmDIALOG_HELP_BUTTON);
  XtSetSensitive(help_button, False);
  
  /* set directory if it is defined */
  if (SaveAsDir) {
    wierdstring = XmStringCreateSimple (SaveAsDir);
    XtSetArg (wargs[0], XmNdirectory, wierdstring);
    XtSetValues (filebox, wargs, 1);
    if (wierdstring) XmStringFree(wierdstring); wierdstring = NULL;
  }

  /* Shazam appear: */
  XtManageChild (filebox);
  XtPopup (XtParent (filebox), XtGrabNone);
  /* all the action is in the callback routine */
} /* end SaveAsCB */

/**
 * Callback for abort button
 * \param widget      widget activated
 * \param clientData  client data
 * \param callData    call data, pointer to TaskMessageBox
 */
static void AbortCB(Widget widget, 
		     XtPointer client_data, 
		     XtPointer call_data) 
{
  TaskMessageBox *tmb = (TaskMessageBox*)client_data;
  XmString  wierdstring = NULL;
  Arg       wargs[5]; 

  tmb->doAbort = TRUE;
  
  /* Write in message box */
  TaskMessageBoxDisplay (tmb, "Abort requested\n");

  /* In task status */
  if (tmb->status != NULL) {
    wierdstring = XmStringCreateSimple ("Abort requested");
    XtSetArg (wargs[0], XmNlabelString, wierdstring);
    XtSetValues (tmb->status, wargs, 1);
    if (wierdstring) XmStringFree(wierdstring); wierdstring = NULL;
  }

} /* end AbortCB */

/**
 * Callback for user response ready
 * \param widget      widget activated (enter key)
 * \param clientData  client data (TaskMessageBox*)
 * \param callData    call data
 */
static void UserInputCB(Widget widget, 
			XtPointer client_data, 
			XtPointer call_data) {
  char     *value=NULL, *blank="";
  char      tmp[200];
  TaskMessageBox *tmb = (TaskMessageBox*)client_data;

  /* Read */
  value = XmTextGetString (tmb->user_input);
  
  /* Copy to object */
  strncpy (tmb->userString, value, 199);
  XtFree(value);

  /* Write in message box */
  g_snprintf (tmp, 199, "%s\n", tmb->userString);
  TaskMessageBoxDisplay (tmb, tmp);

  /* Blank */
  XmTextSetString (tmb->user_input, blank);

  /* Reset status */
  TaskMessageBoxStatus(tmb, "Task Running");

  /* Is something waiting for this? */
  if (tmb->doUserInput) {
    tmb->doUserInput = FALSE;  /* No longer waiting */
    /* Return response */
    TaskMessageUserResponseReply (tmb->taskID, tmb->userString);
  }
  
} /* end UserInputCB */

/**
 * Show TaskMessage dialog
 */
void TaskMessageBoxShow(TaskMessageBox* tmb) {
  if (!XtIsManaged(tmb->task_dialog)) {
    XtManageChild(tmb->task_dialog);
  }
} /* end TaskMessageBoxShow */

