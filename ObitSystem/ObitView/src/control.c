/* $Id: control.c,v 1.3 2005/09/01 11:42:39 bcotton Exp $  */
/*    control panel functions for ObitView */
/*-----------------------------------------------------------------------
*  Copyright (C) 1996,2002-2008
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
#include <Xm/DrawingA.h> 
#include <Xm/MainW.h> 
#include <stdlib.h>
#include <stdio.h>
#include <Xm/Scale.h>
#include <Xm/Form.h>
#include <Xm/PushB.h>
#include <Xm/Label.h>
#include <Xm/ToggleB.h>
#include "imagedisp.h"
#include "color.h"

/**
 *  \file control.c
 * image display control functions.
 */

/* callback functions */
/**
 * Callback for contrast slider
 * \param w           widget activated
 * \param clientData  client data
 * \param callData    call data
 */
void ContrastCB (Widget w, XtPointer clientData, XtPointer callData)
{
  ImageDisplay *IDdata = (ImageDisplay *)clientData;
  XmScaleCallbackStruct *call_data = (XmScaleCallbackStruct *) callData;
  
  /* read value of scrollbar */
  IDdata->value[1] = call_data->value;
  
  /* reset color table */
  SetColorTable(IDdata);
} /* end ContrastCB */

/**
 * Callback for brightness slider
 * \param w           widget activated
 * \param clientData  client data
 * \param callData    call data
 */
void BrightnessCB (Widget w, XtPointer clientData, XtPointer callData)
{
  ImageDisplay *IDdata = (ImageDisplay *)clientData;
  XmScaleCallbackStruct *call_data = (XmScaleCallbackStruct *) callData;
  
  /* read value of scrollbar */
  IDdata->value[0] = call_data->value;
  
  /* reset color table */
  SetColorTable(IDdata);
} /* end BrightnessCB */


/**
 * Make brightness/contrast and info panel
 * \param mainWindow  Window on which to place widgets
 * \param data        call data
 */
Widget MakeControl(Widget mainWindow, XtPointer data)
{
  Widget control, BriScroll, ConScroll, toplab, brilab, conlab;
  XmString     topstr = XmStringCreateSimple ("Control panel");
  XmString     bristr = XmStringCreateSimple ("Brightness");
  XmString     constr = XmStringCreateSimple ("Contrast");
  XmString     blankstr = XmStringCreateSimple ("");
  ImageDisplay *IDdata = (ImageDisplay *)data;
  
  /* make Form widget for control/info */
  control = XtVaCreateManagedWidget ("control", xmFormWidgetClass,
				     mainWindow,
				     XmNwidth,           CONTROLWIDTH,
				     XmNheight,          400,
				     XmNtopAttachment,  XmATTACH_FORM,
				     XmNleftAttachment,  XmATTACH_FORM,
				     NULL);
  
  
  /* panel label */
  toplab = XtVaCreateManagedWidget ("CPtoplab", xmLabelWidgetClass, 
				    control, 
				    XmNwidth,           CONTROLWIDTH,
				    XmNlabelString,   topstr,
				    XmNtopAttachment,  XmATTACH_FORM,
				    XmNleftAttachment,  XmATTACH_FORM,
				    NULL);
  
  /* add scroll bar ("scale" in X) controls */
  
  /* Brightness */
  /*  label */
  brilab = XtVaCreateManagedWidget ("CPbrilab", xmLabelWidgetClass, 
				    control, 
				    XmNwidth,           CONTROLWIDTH,
				    XmNlabelString,   bristr,
				    XmNtopAttachment, XmATTACH_WIDGET,
				    XmNtopWidget,     toplab, 
				    XmNleftAttachment,  XmATTACH_FORM,
				    NULL);
  
  BriScroll = XtVaCreateManagedWidget ("briscroll", xmScaleWidgetClass, 
				       control, 
				       XmNwidth,   CONTROLWIDTH,
				       XmNmaximum,          256,
				       XmNminimum,            1,
				       XmNvalue,            128,
				       XmNshowValue,       True,
				       XmNscaleMultiple,      1,
				       XmNorientation, XmHORIZONTAL,
				       XmNprocessingDirection, XmMAX_ON_RIGHT,
				       XmNtopAttachment, XmATTACH_WIDGET,
				       XmNtopWidget,     brilab, 
				       NULL);
  XtAddCallback(BriScroll, XmNvalueChangedCallback, BrightnessCB, data);
  XtAddCallback(BriScroll, XmNdragCallback, BrightnessCB, data);
  IDdata->BriScroll = BriScroll;
  
  /* Contrast */  
  /*  label */
  conlab = XtVaCreateManagedWidget ("CPconlab", xmLabelWidgetClass, 
				    control, 
				    XmNwidth,        CONTROLWIDTH,
				    XmNlabelString,   constr,
				    XmNtopAttachment, XmATTACH_WIDGET,
				    XmNtopWidget,     BriScroll,
				    XmNleftAttachment,  XmATTACH_FORM,
				    NULL);
  
  ConScroll = XtVaCreateManagedWidget ("conscroll", xmScaleWidgetClass, 
				       control,
				       XmNwidth,   CONTROLWIDTH,
				       XmNmaximum,          256,
				       XmNminimum,            1,
				       XmNvalue,            128,
				       XmNscaleMultiple,      1,
				       XmNshowValue,       True,
				       XmNorientation, XmHORIZONTAL,
				       XmNprocessingDirection, XmMAX_ON_RIGHT,
				       XmNtopAttachment, XmATTACH_WIDGET,
				       XmNtopWidget,     conlab,
				       NULL);
  
  XtAddCallback(ConScroll, XmNvalueChangedCallback, ContrastCB, data);
  XtAddCallback(ConScroll, XmNdragCallback, ContrastCB, data);
  IDdata->ConScroll = ConScroll;
  
  /* info label widgets */
  IDdata->Info1 = XtVaCreateManagedWidget ("Info1", xmLabelWidgetClass, 
					   control, 
					   XmNwidth,        CONTROLWIDTH,
					   XmNlabelString,   blankstr,
					   XmNtopAttachment, XmATTACH_WIDGET,
					   XmNtopWidget,     ConScroll,
					   XmNleftAttachment,  XmATTACH_FORM,
					   NULL);
  
  
  IDdata->Info2 = XtVaCreateManagedWidget ("Info2", xmLabelWidgetClass, 
					   control, 
					   XmNwidth,        CONTROLWIDTH,
					   XmNlabelString,   blankstr,
					   XmNtopAttachment, XmATTACH_WIDGET,
					   XmNtopWidget, IDdata->Info1,
					   XmNleftAttachment,  XmATTACH_FORM,
					   NULL);
  
  
  IDdata->Info3 = XtVaCreateManagedWidget ("Info3", xmLabelWidgetClass, 
					   control, 
					   XmNwidth,        CONTROLWIDTH,
					   XmNlabelString,   blankstr,
					   XmNtopAttachment, XmATTACH_WIDGET,
					   XmNtopWidget, IDdata->Info2,
					   XmNleftAttachment,  XmATTACH_FORM,
					   NULL);
  
  
  IDdata->Info4 = XtVaCreateManagedWidget ("Info4", xmLabelWidgetClass, 
					   control, 
					   XmNwidth,         CONTROLWIDTH,
					   XmNlabelString,   blankstr,
					   XmNtopAttachment, XmATTACH_WIDGET,
					   XmNtopWidget, IDdata->Info3,
					   XmNleftAttachment,  XmATTACH_FORM,
					   NULL);
  
  
  IDdata->Info5 = XtVaCreateManagedWidget ("Info5", xmLabelWidgetClass, 
					   control, 
					   XmNwidth,          CONTROLWIDTH,
					   XmNlabelString,   blankstr,
					   XmNtopAttachment, XmATTACH_WIDGET,
					   XmNtopWidget,  IDdata->Info4,
					   XmNleftAttachment,  XmATTACH_FORM,
					   NULL);
  
  
  IDdata->Info6 = XtVaCreateManagedWidget ("Info6", xmLabelWidgetClass, 
					   control, 
					   XmNwidth,         CONTROLWIDTH,
					   XmNlabelString,   blankstr,
					   XmNtopAttachment, XmATTACH_WIDGET,
					   XmNtopWidget,  IDdata->Info5,
					   XmNleftAttachment,  XmATTACH_FORM,
					   NULL);
  
  IDdata->Info7 = XtVaCreateManagedWidget ("Info7", xmLabelWidgetClass, 
					   control, 
					   XmNwidth,         CONTROLWIDTH,
					   XmNlabelString,   blankstr,
					   XmNtopAttachment, XmATTACH_WIDGET,
					   XmNtopWidget,  IDdata->Info6,
					   XmNleftAttachment,  XmATTACH_FORM,
					   NULL);
  
  
  XtManageChild (control);
  
  if (topstr) XmStringFree(topstr); topstr = NULL;
  if (bristr) XmStringFree(bristr); bristr = NULL;
  if (constr) XmStringFree(constr); constr = NULL;
  if (blankstr) XmStringFree(blankstr); blankstr = NULL;
  return control;
} /* end MakeControl */

