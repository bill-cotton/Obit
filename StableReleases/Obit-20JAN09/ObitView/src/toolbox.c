/* $Id$ */
/* Toolbox panel functions for ObitView */
/*-----------------------------------------------------------------------
*  Copyright (C) 1998,2002-2008
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
/*#include <X11/xpm.h> */
#include "imagedisp.h"
#include "menu.h"
#include "infobox.h"
#include "moviebox.h"
#include "blinkbox.h"
#include "optionbox.h"
#include "markpos.h"
#include "lookpos.h"
#include "helpbox.h"
#include "logger.h"

/**
 *  \file tookbox.c
 * displays tookbox (toolbar?).
 */

/*--------------- file global data ----------------*/
/* button pixmaps */
/*#include "pixmap/open.xpm"*/

/* callback functions */
/**
 * Callback for Open
 * \param w           widget activated
 * \param clientData  ImageDisplay
 * \param callData    call data
 */
void TBOpenCB (Widget w, XtPointer clientData, XtPointer callData)
{
  /*    ImageDisplay *IDdata = (ImageDisplay *)clientData;*/
  /*    XmScaleCallbackStruct *call_data = (XmScaleCallbackStruct *) callData;*/
  
  /* read value of scrollbar */
  /*    IDdata->value[1] = call_data->value;*/
  
  /* reset color table */
  /*    SetColorTable(IDdata);*/
} /* end OpenCB */


/**
 * Make toolbox panel
 * \param mainWindow  Main window
 * \param topWidget   Widget onto which to attach
 * \param data        call data
 * \return toolBox widget
 */
Widget MakeToolbox(Widget mainWindow, Widget topWidget, XtPointer data)
{
  Widget       toolbox, tools[5];
  int          iTool = 0;
  XmString     label = XmStringCreateSimple ("Open");
  ImageDisplay *IDdata = (ImageDisplay *)data;
  /*  XpmAttributes *attributes;*/
  /*Drawable     draw = (Drawable)XtWindow(mainWindow);*/
  /*Display      *display = XtDisplay(mainWindow);*/
  
  /* make Form widget for toolbox - same width as the Control */
  toolbox = XtVaCreateManagedWidget ("toolbox", xmFormWidgetClass,
				     mainWindow,
				     XmNwidth,           CONTROLWIDTH,
				     XmNheight,          200,
				     XmNtopAttachment,  XmATTACH_WIDGET,
				     XmNtopWidget,       topWidget,
				     XmNleftAttachment,  XmATTACH_FORM,
				     XmNbottomAttachment,  XmATTACH_FORM,
				     NULL);
  
  /*------------------------- new row ------------------------------*/
  /* Open button */
  /*  XpmCreatePixmapFromData (display, draw, OpenPM, buttonPM, shapemask, attributes);*/
  tools[iTool] = XtVaCreateManagedWidget ("Open", xmPushButtonWidgetClass, 
					  toolbox, 
					  XmNwidth,           50,
					  XmNheight,          50,
					  XmNtopAttachment, XmATTACH_FORM,
					  XmNleftAttachment,  XmATTACH_FORM,
					  /*				    XmNlabelPixmap,  buttonPM,*/
					  NULL);
  XtAddCallback (tools[iTool++], XmNactivateCallback, OpenCB, (XtPointer)IDdata);
  
  /* Preview button */
  /*  XpmCreatePixmapFromData (display, draw, OpenPM, buttonPM, shapemask, attributes);*/
  tools[iTool] = XtVaCreateManagedWidget ("preview", xmPushButtonWidgetClass, 
					  toolbox, 
					  XmNwidth,           50,
					  XmNheight,          50,
					  XmNtopAttachment, XmATTACH_FORM,
					  XmNleftAttachment,  XmATTACH_WIDGET,
					  XmNleftWidget,    tools[iTool-1],
					  /*				    XmNlabelPixmap,  buttonPM,*/
					  NULL);
  XtAddCallback (tools[iTool++], XmNactivateCallback, PreviewCB, (XtPointer)IDdata);
  
  /* Info button */
  /*  XpmCreatePixmapFromData (display, draw, OpenPM, buttonPM, shapemask, attributes);*/
  tools[iTool] = XtVaCreateManagedWidget ("info", xmPushButtonWidgetClass, 
					  toolbox, 
					  XmNwidth,           50,
					  XmNheight,          50,
					  XmNtopAttachment, XmATTACH_FORM,
					  XmNleftAttachment,  XmATTACH_WIDGET,
					  XmNleftWidget,    tools[iTool-1],
					  /*				    XmNlabelPixmap,  buttonPM,*/
					  NULL);
  XtAddCallback (tools[iTool++], XmNactivateCallback, InfoBoxCB, (XtPointer)IDdata);
  
  /*------------------------- new row ------------------------------*/
  /* Movie button */
  /*  XpmCreatePixmapFromData (display, draw, OpenPM, buttonPM, shapemask, attributes);*/
  tools[iTool] = XtVaCreateManagedWidget ("movie", xmPushButtonWidgetClass, 
					  toolbox, 
					  XmNwidth,           50,
					  XmNheight,          50,
					  XmNleftAttachment,  XmATTACH_FORM,
					  XmNtopAttachment, XmATTACH_WIDGET,
					  XmNtopWidget,    tools[iTool-3],
					  /*				    XmNlabelPixmap,  buttonPM,*/
					  NULL);
  XtAddCallback (tools[iTool++], XmNactivateCallback, MovieBoxCB, (XtPointer)IDdata);
  
  /* Blink swap button */
  /*  XpmCreatePixmapFromData (display, draw, OpenPM, buttonPM, shapemask, attributes);*/
  tools[iTool] = XtVaCreateManagedWidget ("swap", xmPushButtonWidgetClass, 
					  toolbox, 
					  XmNwidth,           50,
					  XmNheight,          50,
					  XmNleftAttachment,  XmATTACH_WIDGET,
					  XmNleftWidget,    tools[iTool-1],
					  XmNtopAttachment, XmATTACH_WIDGET,
					  XmNtopWidget,    tools[iTool-1],
					  /*				    XmNlabelPixmap,  buttonPM,*/
					  NULL);
  XtAddCallback (tools[iTool++], XmNactivateCallback, BlinkSwapCB, (XtPointer)IDdata);
  
  /* Blink blink button */
  /*  XpmCreatePixmapFromData (display, draw, OpenPM, buttonPM, shapemask, attributes);*/
  tools[iTool] = XtVaCreateManagedWidget ("blink", xmPushButtonWidgetClass, 
					  toolbox, 
					  XmNwidth,           50,
					  XmNheight,          50,
					  XmNleftAttachment,  XmATTACH_WIDGET,
					  XmNleftWidget,    tools[iTool-1],
					  XmNtopAttachment, XmATTACH_WIDGET,
					  XmNtopWidget,    tools[iTool-3],
					  /*				    XmNlabelPixmap,  buttonPM,*/
					  NULL);
  XtAddCallback (tools[iTool++], XmNactivateCallback, BlinkBlinkCB, (XtPointer)IDdata);
  
  /*------------------------- new row ------------------------------*/
  XtManageChild (toolbox);
  
  /* delete temporary strings */
  if (label) XmStringFree(label); label = NULL;
  
  return toolbox;
} /* end MakeToolbox */

