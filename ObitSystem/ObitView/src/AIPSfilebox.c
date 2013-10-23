/* $Id$  */
/* AIPS image selection box for ObitView */
/*-----------------------------------------------------------------------
*  Copyright (C) 2013
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
#include <Xm/PanedW.h>
#include <Xm/Separator.h>
#include <Xm/MessageB.h>
#include <Xm/TextF.h>
#include <Xm/Text.h>
#include <Xm/List.h>
#include "imagedisp.h"
#include "helpbox.h"
#include "messagebox.h"
#include "AIPSfilebox.h"
#include "ObitAIPSDir.h"

/**
 *  \file AIPSfilebox.c
 * AIPS image selection box
 */

/*--------------- file global data ----------------*/
/** is the option box active? */
int AIPSFileBoxActive = 0;

/* global structure for things to talk to each other */
AIPSFileBoxStuff AFBdia;

/* button callbacks */
/**
 * Callback for Open button hit 
 * \param w           widget activated
 * \param clientData  client data
 * \param callData    call data
 */
static void AFBOpenButCB (Widget w, XtPointer clientData, XtPointer callData)
{
  gchar   *value=NULL;
  olong   temp;
  
  /* read AIPS directory value */
  value = XmTextFieldGetString (AFBdia.dir);
  if (!value) /* error */
    {MessageShow ("Error reading requested directory");
    return;}
  if (value) XtFree(value);value = NULL;

  /* Save info from dialog */
  FStrngFill (AFBdia.AIPS_dir, value);
  if (value) XtFree(value);value = NULL;
  value =   XmTextFieldGetString (AFBdia.nameT);
  FStrngFill (AFBdia.AName, value);
  if (value) XtFree(value);value = NULL;
  value =   XmTextFieldGetString (AFBdia.classT);
  FStrngFill (AFBdia.AClass, value);
  if (value) XtFree(value);value = NULL;
  value =   XmTextFieldGetString (AFBdia.seqT);
  if (!sscanf (value, "%d", &temp))
    { /* error */
      MessageShow ("Error reading requested sequence");
      if (value) XtFree(value); value = NULL;
      return;}
  AFBdia.ASeq = temp;
  if (value) XtFree(value);value = NULL;
  value =   XmTextFieldGetString (AFBdia.userT);
  if (!sscanf (value, "%d", &temp))
    { /* error */
      MessageShow ("Error reading AIPS user number");
      if (value) XtFree(value); value = NULL;
      return;}
  AFBdia.AUser = temp;
  if (value) XtFree(value);value = NULL;

  /* make it disappear but still exist */
  XtPopdown(AFBdia.dialog);
  AIPSFileBoxActive = 0;   /* No longer active */

} /* end AFBOpenButCB */

/**
 * Callback for Cancel button hit 
 * \param w           widget activated
 * \param clientData  client data
 * \param callData    call data
 */
static void AFBCancelButCB (Widget w, XtPointer clientData, XtPointer callData)
{
  /* make it disappear but still exist */
  XtPopdown(AFBdia.dialog);

  /* Want to cancel editing */
  AFBdia.cancel     = 1;
  AIPSFileBoxActive = 0;   /* No longer active */

} /* end AFBCancelButCB */

/**
 * Callback for change AIPS directory
 * \param text_w,    Text window with directory
 * \param clientData Client data
 * \param callData   AFBdia structure
 */
static void NewADirCB (Widget text_w, XtPointer clientData, XtPointer callData)
{
  gchar   *value=NULL, line[104], time[24];
  olong   disk, ncat, cno, user, temp, position=1;
  ObitAIPSDirCatEntry *entry=NULL;
  XmString item;
 
  /* read value */
  value = XmTextFieldGetString (text_w);
  if (!value) { /* error */
    MessageShow ("Error reading requested directory");
    return;
  }
  /* Just return if blank */
  if ((value[0]==' ') && (value[1]==' ')) 
    {if (value) XtFree(value);return;}

  /* Bail if initial message */
  if (!strncmp (value, "specify AIPS ", 12)) 
    {if (value) XtFree(value); return;}

  /* Better be a full path - feeble test */
  if ((value[0]!='/') && (value[1]!='/')) {
    MessageShow ("You must specify a full path");
    if (value) XtFree(value); return;
  }

  /* OK, save directory name */
  ObitTrimTrail(value);  /* cut trailing blanks */
  /* Add slash to be sure */
  if (value[strlen(value)-1]!='/') sprintf(line,"%s/", value);
  else sprintf(line,"%s", value);
  FStrngFill (AFBdia.AIPS_dir, line);
  if (value) XtFree(value);value = NULL;

  /* Get user number */
  value =   XmTextFieldGetString (AFBdia.userT);
  if (!sscanf (value, "%d", &temp))
    { /* error */
      MessageShow ("Error reading AIPS user number");
      if (value) XtFree(value); value = NULL;
      return;}
  AFBdia.AUser = temp;
  if (value) XtFree(value);value = NULL;

  /* Set directory */
  disk = 1;
  user = AFBdia.AUser;
  ObitAIPSSetDirname(disk, AFBdia.AIPS_dir->sp, err);
  if (err->error) {
      Obit_log_error(err, OBIT_Error, "Error setting AIPS directory"); ObitErrLog(err); return;
  }

  ncat = ObitAIPSDirNumber(disk, user, err);
  if (err->error) {
      Obit_log_error(err, OBIT_Error, "Error finding max. cat entry, possibly bad directory"); 
      ObitErrLog(err); 
      return;
  }
  /* Loop over images */  
  for (cno=1; cno<=ncat; cno++) {
    entry = ObitAIPSDirGetEntry(disk, user, cno, err);
    if (err->error) {
      Obit_log_error(err, OBIT_Error, "Error reading entry"); ObitErrLog(err); return;
    }
    /* exists and an Image */
    if (entry) {
      if (entry->type[0]=='M' && entry->type[1]=='A' && entry->user==user) {
	ObitAIPSDirGetAccess (entry, time); time[20] = 0;
	g_snprintf (line, 100, "%5d %-12.12s.%-6.6s. %4d %s",
		    cno, entry->name, entry->class, entry->seq, time);
	/*fprintf (stderr, "%s\n", line);   DEBUG */
	item = XmStringCreateLocalized (line);
	XmListAddItem (AFBdia.SList, item, position++);
	if (item) XmStringFree (item);item = NULL;
  } /* end if image */
      g_free (entry);
    } /* end if exists */
  } /* end loop reading */

  /* select previous */
  if ((AFBdia.selPos>0) && (AFBdia.selPos<position)) {
    XmListSetPos (AFBdia.SList, AFBdia.selPos);
    XmListSelectPos (AFBdia.SList, AFBdia.selPos, False);
  }
} /* end NewADirCB */

/**
 * Callback for Select AIPS Image
 * \param text_w,    Text window with directory
 * \param clientData Client data
 * \param callData   selection info
 */
static void SelectAImageCB (Widget text_w, XtPointer client_data, XtPointer call_data)
{
  XmListCallbackStruct *cbs = (XmListCallbackStruct *) call_data;
    char *choice, sub[50];

    choice = (char *) XmStringUnparse (cbs->item,
				       XmFONTLIST_DEFAULT_TAG,
				       XmCHARSET_TEXT,
				       XmCHARSET_TEXT,
				       NULL, 0,
				       XmOUTPUT_ALL);
    /* Save which one */
    AFBdia.selPos = cbs->item_position;
    /* Update AIPS file info */
    strncpy(sub, &choice[6], 12); sub[12] = 0;
    XmTextFieldSetString (AFBdia.nameT, sub);
    strncpy(sub, &choice[19], 6); sub[6] = 0;
    XmTextFieldSetString (AFBdia.classT, sub);
    strncpy(sub, &choice[26], 6); sub[6] = 0;
    XmTextFieldSetString (AFBdia.seqT, sub);
    XtFree (choice);
} /* end SelectAImageCB */

/**
 * Dialog box for selecting an AIPS image
 *  w       = widget activated
 *  image   = Image structure to be updates
 */
AIPSFileBoxStuff* AIPSFileBox (Widget w, ImageData *image)
{
  Widget  form=NULL, labelB1=NULL, labelB2=NULL, labelB3=NULL, labelB4=NULL;
  Widget OpenButton=NULL, CancelButton=NULL, HelpButton=NULL, sep=NULL;
  XmString     label=NULL;
  gchar        tstr[64];
#define AIPSFILEBOX_WIDTH 400
#define AIPSFILEBOX_HEIGHT 480
  
  AFBdia.cancel = 0;  /* Unset cancel request */

  /* don't make another one */
  if (AIPSFileBoxActive) {
    if (XtIsRealized (AFBdia.dialog))
      XMapRaised (XtDisplay(AFBdia.dialog), XtWindow(AFBdia.dialog));
    
    /* bring it back where we can see it */
    XtPopup(AFBdia.dialog, XtGrabNonexclusive);
    
    /* put it some place reasonable */
    /*  where is parent? */
    XtVaGetValues (Display_shell, XmNx, &AFBdia.xpos, XmNy, &AFBdia.ypos,  NULL);
    AFBdia.ypos += 170;
    if (AFBdia.xpos<0) AFBdia.xpos += 100;
    XMoveWindow (XtDisplay(Display_shell), XtWindow(AFBdia.dialog), 
		 AFBdia.xpos, AFBdia.ypos);
    
    /* set value of AIPS Directory  */
    XmTextFieldSetString (AFBdia.dir, AFBdia.AIPS_dir->sp);

    /* User number */
    sprintf (tstr, "%d", AFBdia.AUser);
    XmTextFieldSetString (AFBdia.userT, tstr);

    /* Set size */
    XtVaSetValues (AFBdia.dialog,
		   XmNwidth,     AIPSFILEBOX_WIDTH,
		   XmNheight,    AIPSFILEBOX_HEIGHT,
		   NULL);
    /* Fill in images */
    XtVaSetValues (AFBdia.Form,
		   XmNwidth,     AIPSFILEBOX_WIDTH,
		   XmNheight,    AIPSFILEBOX_HEIGHT,
		   NULL);
    XtManageChild (AFBdia.dialog);
    XtManageChild (AFBdia.Form);
    XtManageChild (AFBdia.Open);
    NewADirCB (AFBdia.dir, (XtPointer)&AFBdia, (XtPointer)&AFBdia); 
    return &AFBdia;
  } /* end of update dialog */
  
  /* mark as active */
  AIPSFileBoxActive = 1;
  
  /* Initialize structure */
  if (!AFBdia.AIPS_dir) AFBdia.AIPS_dir = MakeFStrngSize(120);
  if (!AFBdia.AName )   AFBdia.AName    = MakeFStrngSize(12);
  if (!AFBdia.AClass)   AFBdia.AClass   = MakeFStrngSize(6);
  FStrngFill (AFBdia.AIPS_dir, image->AIPSDir->sp);
  FStrngFill (AFBdia.AName,    "            ");
  FStrngFill (AFBdia.AClass,   "      ");
  AFBdia.cno   = 0;
  AFBdia.AUser = image->AIPSuser;
  AFBdia.ASeq  = 0;
  AFBdia.AClass   = MakeFStrngSize(6);
  /* FStrngFill (AFBdia.AIPS_dir, "/export/data_1/GOLLUM_1/"); DEBUG */
  /* AFBdia.AUser = 100;     DEBUG */
  /* Create window */
  AFBdia.dialog = XtVaCreatePopupShell ("AIPS Image Selection", xmDialogShellWidgetClass, 
					Display_shell, 
					XmNautoUnmanage, False,
					XmNwidth,     AIPSFILEBOX_WIDTH,
					XmNheight,    AIPSFILEBOX_HEIGHT,
					XmNdeleteResponse, XmDESTROY,
					NULL);

  /* form to hang things on */
  form =  XmCreateForm (AFBdia.dialog, "AIPSpanel", NULL, 0);
  XtVaSetValues (form,
		 XmNwidth,     AIPSFILEBOX_WIDTH,
		 XmNheight,    AIPSFILEBOX_HEIGHT,
		 XmNx,           0,
		 XmNy,           0,
		 NULL);
  AFBdia.Form = form;

  /* Scrolled list for image selection */
  AFBdia.SList = XmCreateScrolledList (form, "scrolled_list",  NULL, 0);
  XtVaSetValues (AFBdia.SList,             /* List part */
		 XmNwidth,            AIPSFILEBOX_WIDTH,
		 XmNcolumns,          64,
		 XmNvisibleItemCount, 20,
		 XmNselectionPolicy,  XmBROWSE_SELECT,
		 NULL);
  XtVaSetValues (XtParent(AFBdia.SList),    /* Scrolling part */
		 XmNwidth,            AIPSFILEBOX_WIDTH,
		 XmNcolumns,          64,
		 XmNselectionPolicy,  XmBROWSE_SELECT,
		 XmNtopAttachment,    XmATTACH_FORM,
		 XmNrightAttachment,  XmATTACH_FORM,
		 XmNleftAttachment,   XmATTACH_FORM,
		 NULL);
  XtAddCallback (AFBdia.SList, XmNbrowseSelectionCallback, SelectAImageCB,
		 (XtPointer)&AFBdia);

  /* AIPS directory text box */
  AFBdia.dir = XmCreateTextField (form, "AIPS dir", NULL, 0);
  XtVaSetValues (AFBdia.dir,
		 XmNwidth,           AIPSFILEBOX_WIDTH,
		 XmNcolumns,         50,
		 XmNtopAttachment,   XmATTACH_WIDGET,
		 XmNtopWidget,       AFBdia.SList,
		 XmNleftAttachment,  XmATTACH_FORM,
		 NULL);
  XmTextFieldSetString (AFBdia.dir, AFBdia.AIPS_dir->sp);
  XtAddCallback (AFBdia.dir, XmNactivateCallback, NewADirCB,
		 (XtPointer)&AFBdia);

  /* User number */
  labelB1 = XmCreateLabel (form, "AIPS User Number:", NULL, 0);
  XtVaSetValues (labelB1,
		 XmNtopAttachment,  XmATTACH_WIDGET,
		 XmNtopWidget,      AFBdia.dir,
		 XmNleftAttachment, XmATTACH_FORM,
		 XmNheight,         30,
		 NULL);
 AFBdia.userT = XmCreateTextField (form, "user", NULL, 0);
 XtVaSetValues (AFBdia.userT,
		XmNcolumns,         6,
		XmNtopAttachment,   XmATTACH_WIDGET,
		XmNtopWidget,       AFBdia.dir,
		XmNleftAttachment,  XmATTACH_WIDGET,
		XmNleftWidget,      labelB1,
		XmNheight,          30,
		NULL);
  sprintf (tstr, "%d", AFBdia.AUser);
  XmTextFieldSetString (AFBdia.userT, tstr);

  /* Separator */
  sep = XmCreateSeparator (form, "sep", NULL, 0);
  XtVaSetValues (sep,
		 XmNtopAttachment,    XmATTACH_WIDGET,
		 XmNtopWidget,        labelB1,
		 XmNleftAttachment,   XmATTACH_FORM,
		 XmNrightAttachment,   XmATTACH_FORM,
		 NULL);

  /* Image name */
  labelB2 = XmCreateLabel (form, "AIPS Name:", NULL, 0);
  XtVaSetValues (labelB2,
		 XmNtopAttachment,    XmATTACH_WIDGET,
		 XmNtopWidget,        sep,
		 XmNleftAttachment,   XmATTACH_FORM,
		 XmNheight,           30,
		 NULL);
  AFBdia.nameT = XmCreateTextField (form, "class", NULL, 0);
  XtVaSetValues (AFBdia.nameT,
		 XmNcolumns,         12,
		 XmNtopAttachment,   XmATTACH_WIDGET,
		 XmNtopWidget,       sep,
		 XmNleftAttachment,  XmATTACH_WIDGET,
		 XmNleftWidget,      labelB2,
		 XmNheight,          30,
		 NULL);
  XmTextFieldSetString (AFBdia.nameT, AFBdia.AClass->sp);

  /* Image class */
  labelB3 = XmCreateLabel (form, "Class:", NULL, 0);
  XtVaSetValues (labelB3,
		 XmNtopAttachment,  XmATTACH_WIDGET,
		 XmNtopWidget,      sep,
		 XmNleftAttachment,  XmATTACH_WIDGET,
		 XmNleftWidget,      AFBdia.nameT,
		 XmNheight,         30,
		 NULL);
 AFBdia.classT = XmCreateTextField (form, "class", NULL, 0);
 XtVaSetValues (AFBdia.classT,
		XmNcolumns,         6,
		XmNtopAttachment,   XmATTACH_WIDGET,
		XmNtopWidget,       sep,
		XmNleftAttachment,  XmATTACH_WIDGET,
		XmNleftWidget,      labelB3,
		XmNheight,          30,
		NULL);
  XmTextFieldSetString (AFBdia.classT, AFBdia.AClass->sp);

  /* Seq number */
  labelB4 = XmCreateLabel (form, "Seq:", NULL, 0);
  XtVaSetValues (labelB4,
		 XmNtopAttachment,   XmATTACH_WIDGET,
		 XmNtopWidget,       sep,
		 XmNleftAttachment,  XmATTACH_WIDGET,
		 XmNleftWidget,      AFBdia.classT,
		 XmNheight,          30,
		 NULL);
 AFBdia.seqT = XmCreateTextField (form, "seq", NULL, 0);
 XtVaSetValues (AFBdia.seqT,
		XmNcolumns,         6,
		XmNtopAttachment,   XmATTACH_WIDGET,
		XmNtopWidget,       sep,
		XmNleftAttachment,  XmATTACH_WIDGET,
		XmNleftWidget,      labelB4,
		/*XmNrightAttachment, XmATTACH_FORM,*/
		XmNheight,          30,
		NULL);
  sprintf (tstr, "%d", AFBdia.ASeq);
  XmTextFieldSetString (AFBdia.seqT, tstr);

  /* Fill in images */
  NewADirCB (AFBdia.dir, (XtPointer)&AFBdia, (XtPointer)&AFBdia); 

   /* Open button */
  OpenButton = XmCreatePushButton (form, "Open", NULL, 0);
  XtVaSetValues (OpenButton,
		 XmNbottomAttachment, XmATTACH_FORM,
		 XmNleftAttachment,  XmATTACH_FORM,
		 NULL);
  XtAddCallback (OpenButton, XmNactivateCallback, AFBOpenButCB, NULL);
  AFBdia.Open = OpenButton;

  /* Cancel button */
  CancelButton =  XmCreatePushButton (form, "Cancel", NULL, 0);
  XtVaSetValues (CancelButton,
		 XmNbottomAttachment, XmATTACH_FORM,
		 XmNleftAttachment,   XmATTACH_WIDGET,
		 XmNleftWidget,       OpenButton,
		 NULL);
  XtAddCallback (CancelButton, XmNactivateCallback, AFBCancelButCB, NULL);
  AFBdia.Cancel = CancelButton;
 
  /* Help button */
  HelpButton =  XmCreatePushButton (form, "Help", NULL, 0);
  XtVaSetValues (HelpButton,
		 XmNbottomAttachment, XmATTACH_FORM,
		 XmNleftAttachment,   XmATTACH_WIDGET,
		 XmNleftWidget,       CancelButton,
		 NULL);
  /*XtAddCallback (HelpButton, XmNactivateCallback,  HelpBoxTopicCB, 
    (XtPointer)"File/Open AIPS");*/
  AFBdia.Help = HelpButton;

  if (label)    XmStringFree(label);    label  = NULL;
  
  /* set it up */
  XtManageChild (AFBdia.SList);
  XtManageChild (AFBdia.dir);
  XtManageChild (labelB1);       XtManageChild (AFBdia.userT);
  XtManageChild (sep);
  XtManageChild (labelB2);       XtManageChild (AFBdia.nameT);
  XtManageChild (labelB3);       XtManageChild (AFBdia.classT);
  XtManageChild (labelB4);       XtManageChild (AFBdia.seqT);
  XtManageChild (OpenButton);
  XtManageChild (CancelButton);
  XtManageChild (HelpButton);
  XtManageChild (form);
  XtManageChild (AFBdia.dialog);
  
  /* put it some place reasonable */
  /*  where is parent? */
  XtVaGetValues (Display_shell,
		 XmNx, &AFBdia.xpos,
		 XmNy, &AFBdia.ypos,
		 NULL);
  AFBdia.ypos += 170;
  AFBdia.xpos += 100;
  if (AFBdia.xpos<0) AFBdia.xpos = 0;
  XMoveWindow (XtDisplay(Display_shell), XtWindow(AFBdia.dialog), 
	       AFBdia.xpos, AFBdia.ypos);

  return &AFBdia;
} /* end AIPSFileBox */

