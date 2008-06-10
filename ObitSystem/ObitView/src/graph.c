/* $Id: graph.c,v 1.5 2005/10/19 19:36:25 bcotton Exp $ */
/*   Graphics overlay plane functions */
/*-----------------------------------------------------------------------
*  Copyright (C) 2005
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
#include <stdio.h> 
#include <string.h> 
#include <math.h>
#include <glib.h>
#include <Xm/DialogS.h> 
#include <Xm/Form.h>
#include <Xm/Label.h>
#include <Xm/PushB.h>
#include "graph.h"
#include "imagedisp.h"
#include "messagebox.h"
#define min(a,b) (((a) < (b)) ? (a) : (b))
#define max(a,b) (((a) > (b)) ? (a) : (b))

/* global structure for things to talk to each other */
/* Graphics overlay color */
typedef struct {
  ImageDisplay *BoxData;
  Widget dialog; /* box */
} GphColBoxStuff;
GphColBoxStuff GphColData;
/* is the GphCol box active? */
int GphColBoxActive = 0;

/**
 *  \file graph.c
 * Manages display overlay graphics
 */
/*---------------Private function prototypes----------------*/
void GphCSetColCB (Widget w, int which, XmToggleButtonCallbackStruct *state);
void GphCSetCancelButCB (Widget w, XtPointer clientData, XtPointer callData);


/*------------------- Public functions -------------------*/
/**
 * Draw line in the graphics overlay of an ImageData
 * Line may be partially off the image.
 * Note: display "y" is reverse of image pixel numbering
 * \param ImgData     ImageData structure pointer
 * \param blc         bottom left corner (0-rel)
 * \param trc         top right corner (0-rel)
 * \param gcolor      color to write, 0 = erase
 */
void GraphDrawLine (ImageData *ImgData, olong blc[2], olong trc[2], olong gcolor)
{
  olong ix, iy, nx, ny, is, ie;
  olong addr, i, n;
  ofloat slope, interc, dx, dy, x, y, d;
  gchar   *gpharray = ImgData->gpharray;
  gchar errMsg[121];

  /* Range check */
  if ((gcolor<0) || (gcolor>MAXGRAPH)) {
    sprintf (errMsg, "GraphDrawLine: color %d out of range [0,%d]", 
	     gcolor, MAXGRAPH);
    MessageShow(errMsg);
    return;
  }
  nx = ImgData->myDesc->inaxes[0];
  ny = ImgData->myDesc->inaxes[1];
  if ((blc[0]<0) || (blc[0]>=nx) || (blc[1]<0) || (blc[1]>=ny) ||
      (trc[0]<0) || (trc[0]>=nx) || (trc[1]<0) || (trc[1]>=ny)) {
    sprintf (errMsg, "GraphDrawLine: line out of range [%d,%d] [%d,%d]", 
	     blc[0], blc[1], trc[0],trc[1]);
    MessageShow(errMsg);
    return;
  }

  dx = (ofloat)(trc[0] - blc[0]);
  dy = (ofloat)(trc[1] - blc[1]);

  /* vertical lines - can go either direction */
  if (fabs(dx) < 0.5) {
    ix = blc[0];
    is = min (blc[1], trc[1]);
    ie = max (blc[1], trc[1]);
    for (iy=is; iy<=ie; iy++) {
      addr = (ny - iy - 1)* nx + ix;
      *(gpharray+addr) = gcolor;
    }
    return;
  }

  /* horizontal lines */
  if (fabs(dy) < 0.5) {
    iy = ny - blc[1] - 1;
    is = min (blc[0], trc[0]);
    ie = max (blc[0], trc[0]);
    for (ix=is; ix<=ie; ix++) {
      addr = iy * nx + ix;
      *(gpharray+addr) = gcolor;
    }
    return;
  }

  /* Line parameters for sloping lines */
  slope = dy / dx;
  interc = ((ofloat)blc[1]) - slope*((ofloat)blc[0]);

  /* determine location every 1/20 pixel in x */
  n = (trc[0] - blc[0]) * 20;
  if (blc[0] < trc[0]) d = 0.05;  /* direction from blc in X */
  else d = -0.05;
  for (i=0; i<n; i++) {
    x =  blc[0] + ((float)i)*d;
    y = interc + x * slope;
    ix = (int)(x + 0.5);
    iy = (int)(y + 0.5);
    /* flip y for display */
    iy = ny - iy - 1;
    /* in image? */
    if ((ix>0) && (ix<nx) && (iy>0) && (iy<ny)) {
      addr = iy * nx + ix;
      *(gpharray+addr) = gcolor;
    }
  } /* end loop over pixels */

} /* end GraphDrawLine */

/**
 * Draw circle in the graphics overlay of an ImageData
 * Line may be partially off the image.
 * Note: display "y" is reverse of image pixel numbering
 * \param ImgData     ImageData structure pointer
 * \param center      center image pixel (0-rel)
 * \param radius      radius in pixels
 * \param gcolor      color to write, 0 = erase
 */
void GraphDrawCircle (ImageData *ImgData, olong center[2], olong radius, olong gcolor)
{
  olong ix, iy, nx, ny;
  olong addr, i, n;
  ofloat x0, y0, x, y, r2, arg;
  gchar   *gpharray = ImgData->gpharray;
  gchar errMsg[121];

  /* Color Range check */
  if ((gcolor<0) || (gcolor>MAXGRAPH)) {
    sprintf (errMsg, "GraphDrawCircle: color %d out of range [0,%d]", 
	     gcolor, MAXGRAPH);
    MessageShow(errMsg);
    return;
  }

  nx = ImgData->myDesc->inaxes[0];
  ny = ImgData->myDesc->inaxes[1];
  if ((center[0]<0) || (center[0]>=nx) || (center[1]<0) || (center[1]>=ny)) {
    sprintf (errMsg, "GraphDrawCircle: center out of range [%d,%d]", 
	     center[0], center[1]);
    MessageShow(errMsg);
    return;
  }

  /* determine location every 1/20 pixel in x */
  r2 = ((ofloat)radius) * ((ofloat)radius);
  x0 = ((ofloat)center[0]);
  y0 = ((ofloat)center[1]);
  n = 2 * radius * 20 + 20;
  for (i=0; i<n; i++) {
    x = ((float)(i-n/2))*0.05;
    ix = (olong)(x0 + x + 0.5);
    arg = r2 - x*x;
    if (arg<0.0) arg = 0.0;
    arg = sqrt (arg);

    /* positive solution */
    y = y0 + arg;
    iy = (olong)(y + 0.5);
    iy = ny - iy - 1;
    /* in image? */
    if ((ix>0) && (ix<nx) && (iy>0) && (iy<ny)) {
      addr = iy * nx + ix;
      *(gpharray+addr) = gcolor;
    }

    /* negative solution */
    y = y0 - arg;
    iy = (olong)(y + 0.5);
    iy = ny - iy - 1;
    /* in image? */
    if ((ix>0) && (ix<nx) && (iy>0) && (iy<ny)) {
      addr = iy * nx + ix;
      *(gpharray+addr) = gcolor;
    }
  } /* end loop over pixels */

}  /* end GraphDrawCircle */

/**
 * Clear graphics overlay
 * \param ImgData     ImageData structure pointer
 */
void GraphClear (ImageData *ImgData)
{
  olong   i, n;
  gchar   *gpharray = ImgData->gpharray;

  if (!gpharray) return;         /* anything to clear? */
  if (!ImgData->myDesc) return;  /* anything to clear? */

  n = (olong)ImgData->myDesc->inaxes[0] * (olong)ImgData->myDesc->inaxes[1];
  for (i=0; i<n; i++) *(gpharray++) = 0;

}  /* end GraphClear */

/**
 * Callback routine to clear graphics overlay
 * \param w           widget activated
 * \param clientData  client data
 * \param callData    call data
 */
void ResetGphCB (Widget w, XtPointer clientData, XtPointer callData)
{
  ImageDisplay *IDdata = (ImageDisplay *)clientData;

  GraphClear(&image[CurImag]);
  /* redraw image */
  PaintImage(IDdata);
} /* end ResetGphCB */

/**
 * Callback for create dialog box for Graphics color
 * \param w           widget activated
 * \param clientData  client data
 * \param callData    call data
 */
void GphCSetCB (Widget parent, XtPointer clientData, XtPointer callData)
{
  int    iGphCol;
  Widget form;
  Widget label1, radio;
  Widget CancelButton;
  XmString  label   = NULL;
  XmString  magenta = NULL;
  XmString  yellow  = NULL;
  XmString  red     = NULL;
  XmString  cyan    = NULL;
  XmString  blue    = NULL;
  XmString  green   = NULL;
  ImageDisplay *IDdata = (ImageDisplay*)clientData;
  
  
  /* register IDdata */
  GphColData.BoxData = IDdata;
  
  /* don't make another one */
  if (GphColBoxActive) {
    if (XtIsRealized (GphColData.dialog))
      XMapRaised (XtDisplay(GphColData.dialog), XtWindow(GphColData.dialog));
    return;
  }
  GphColBoxActive = 1;
  
  label   = XmStringCreateSimple ("Set Color for graphics");
  magenta = XmStringCreateSimple ("Magenta");
  yellow  = XmStringCreateSimple ("Yellow");
  red     = XmStringCreateSimple ("Red");
  cyan    = XmStringCreateSimple ("Cyan");
  blue    = XmStringCreateSimple ("Blue");
  green   = XmStringCreateSimple ("Green");
  
  GphColData.dialog = XtVaCreatePopupShell ("GphColBox", xmDialogShellWidgetClass, 
					    IDdata->shell, 
					    XmNautoUnmanage, False,
					    XmNwidth,     180,
					    XmNheight,   220,
					    XmNdeleteResponse, XmDESTROY,
					    NULL);
  
  /* make Form widget to stick things on */
  form = XtVaCreateManagedWidget ("OptionForm", xmFormWidgetClass,
				  GphColData.dialog,
				  XmNautoUnmanage, False,
				  XmNwidth,     180,
				  XmNheight,    220,
				  XmNx,           0,
				  XmNy,           0,
				  NULL);
  
  /* info label widgets */
  label1 = XtVaCreateManagedWidget ("Label1", xmLabelWidgetClass, 
				    form, 
				    XmNwidth,           180,
				    XmNlabelString,   label,
				    XmNtopAttachment, XmATTACH_FORM,
				    XmNleftAttachment,  XmATTACH_FORM,
				    NULL);
  
  /* GphCol radio buttons */
  iGphCol = IDdata->gcolor-1;
  radio = XmVaCreateSimpleRadioBox(form, "GraphColor", iGphCol, 
				   (XtCallbackProc) GphCSetColCB,
				   XmNwidth,           180,
				   XmNtopAttachment, XmATTACH_WIDGET,
				   XmNtopWidget,     label1,
				   XmNleftAttachment,  XmATTACH_FORM,
				   XmVaRADIOBUTTON, magenta, NULL, NULL, NULL,
				   XmVaRADIOBUTTON, yellow,  NULL, NULL, NULL,
				   XmVaRADIOBUTTON, red,     NULL, NULL, NULL,
				   XmVaRADIOBUTTON, cyan,    NULL, NULL, NULL,
				   XmVaRADIOBUTTON, blue,    NULL, NULL, NULL,
				   XmVaRADIOBUTTON, green,   NULL, NULL, NULL,
				   NULL);
  XtManageChild(radio);
  
  /* Cancel button */
  CancelButton = XtVaCreateManagedWidget ("Cancel", xmPushButtonWidgetClass, 
					  form, 
					  XmNbottomAttachment, XmATTACH_FORM,
					  NULL);
  XtAddCallback (CancelButton, XmNactivateCallback, GphCSetCancelButCB, 
		 (XtPointer)IDdata);
  
  if (label)   XmStringFree(label);   label   = NULL;
  if (magenta) XmStringFree(magenta); magenta = NULL;
  if (yellow)  XmStringFree(yellow);  yellow  = NULL;
  if (red)     XmStringFree(red);     red     = NULL;
  if (cyan)    XmStringFree(cyan);    cyan    = NULL;
  if (blue)    XmStringFree(blue);    blue    = NULL;
  if (green)   XmStringFree(green);   green   = NULL;
  
  /* set it up */
  XtManageChild (GphColData.dialog);
  
} /* end SetEquCB */

/**
 * Callback for Set Graphics Color button 
 * \param w       widget activated
 * \param which   =0=>use header, 1=>J2000 2=> B1950 
 * \param state   Button state
 */
void GphCSetColCB (Widget w, int which, XmToggleButtonCallbackStruct *state)
{
  ImageDisplay *IDdata = GphColData.BoxData;

  if (!state->set) return; /* ignore buttons being cleared */
  
  /* set color */
  IDdata->gcolor = which+1;

  /* done with box */
  XtDestroyWidget (GphColData.dialog);
  GphColBoxActive = 0;
} /* end GphCSetCB */

/* button callbacks */
/**
 * Callback for Cancel button hit
 * \param w           widget activated
 * \param clientData  client data
 * \param callData    call data
 */
void GphCSetCancelButCB (Widget w, XtPointer clientData, XtPointer callData)
{
  if (!GphColBoxActive) return;
  XtDestroyWidget (GphColData.dialog);
  GphColBoxActive = 0;
} /* end GphCSetCancelButCB */

