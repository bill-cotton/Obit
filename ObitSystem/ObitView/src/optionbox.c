/* $Id$  */
/* Option dialog box  for ObitView */
/*-----------------------------------------------------------------------
*  Copyright (C) 1996,2002-2026
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
#include <Xm/Scale.h>
#include <Xm/Form.h>
#include <Xm/PushB.h>
#include <Xm/Label.h>
#include <Xm/ToggleB.h>
#include <Xm/ToggleBG.h>
#include <Xm/RowColumn.h>
#include <Xm/Separator.h>
#include <Xm/MessageB.h>
#include <Xm/TextF.h>
#include <Xm/Text.h>
#include "obitview.h"
#include "imagedisp.h"
#include "messagebox.h"
#include "Image2Pix.h"
#include "XMLRPCserver.h"

/**
 *  \file optionbox.c
 * displays "option" dialog.
 */

/*--------------- file global data ----------------*/
/** is the option box active? */
int OptionBoxActive = 0;

/* global structure for things to talk to each other */
typedef struct {
  ImageDisplay *BoxData;
  Widget dialog, pixran1, pixran2, planelab, data1, data2;/* box, min, max */
  Widget ERtimeout; /* Edit request timeout (sec) */
  Widget plane;            /* Plane Number */
  Widget dim4, dim5, dim6; /* higher dimensions */
  short  xpos, ypos; /* location of the box */
} OptionBoxStuff;
OptionBoxStuff dia;

/**
 * Callback for minimum pixel value
 * \param w           widget activated
 * \param clientData  client data
 * \param callData    call data
 */
void ReadMinCB (Widget w, XtPointer clientData, XtPointer callData)
{
  /*ImageDisplay *IDdata = (ImageDisplay *)clientData;*/
  char     *value=NULL;
  float    temp;
  
  /* read value */
  value = XmTextGetString (w);
  if (!value) /* error */
    {MessageShow ("Error reading minimum pixel value");
    return;}
  if (!sscanf (value, "%e", &temp))
    { /* error */
      MessageShow ("Error reading minimum pixel value");
      if (value) {XtFree(value);} value = NULL;
      return;}
  if (value) {XtFree(value);} value = NULL;
  
  /* OK, save */
  image[CurImag].PixRange[0] = temp;
  
} /* end ReadMinCB */

/**
 * Callback for maximum pixel value
 * \param w           widget activated
 * \param clientData  client data
 * \param callData    call data
 */
void ReadMaxCB (Widget w, XtPointer clientData, XtPointer callData)
{
  /*ImageDisplay *IDdata = (ImageDisplay *)clientData;*/
  char     *value=NULL;
  float    temp;
  
  /* read value */
  value = XmTextGetString (w);
  if (!value) /* error */
    {MessageShow ("Error reading maximum pixel value");
    return;}
  if (!sscanf (value, "%e", &temp))
    { /* error */
      MessageShow ("Error reading maximum pixel value");
      if (value) {XtFree(value);} value = NULL;
      return;}
  if (value) {XtFree(value);} value = NULL;
  
  /* OK, save */
  image[CurImag].PixRange[1] = temp;
  
} /* end ReadMaxCB */

/**
 * Callback for Edit request timeout value
 * \param w           widget activated
 * \param clientData  client data
 * \param callData    call data
 */
void ReadERTimeoutCB (Widget w, XtPointer clientData, XtPointer callData)
{
  /*ImageDisplay *IDdata = (ImageDisplay *)clientData;*/
  char     *value=NULL;
  float    temp;
  
  /* read value */
  value = XmTextGetString (w);
  if (!value) /* error */
    {MessageShow ("Error reading timeout value");
    return;}
  if (!sscanf (value, "%e", &temp))
    { /* error */
      MessageShow ("Error reading timeout value");
      if (value) {XtFree(value);} value = NULL;
      return;}
  if (value) {XtFree(value);} value = NULL;
  
  /* OK, save - lock mutex to avoid collision with inicoming request */
  pthread_mutex_lock(&request_lock); /* lock mutex */
  ERtime_out = temp;
  pthread_mutex_unlock(&request_lock); /* unlock mutex */
  
} /* end ReadERTimeoutCB */

/**
 * Callback for read Plane number
 * \param w           widget activated
 * \param clientData  client data
 * \param callData    call data
 */
void ReadPlaneCB (Widget w, XtPointer clientData, XtPointer callData)
{
  /*ImageDisplay *IDdata = (ImageDisplay *)clientData;*/
  char     *value=NULL;
  int      itemp;
  
  /* read value */
  value = XmTextGetString (w);
  if (!value) /* error */
    {MessageShow ("Error reading plane number");
    return;}
  if (!sscanf (value, "%d", &itemp))
    { /* error */
      MessageShow ("Error reading plane number");
      if (value) {XtFree(value);} value = NULL;
      return;}
  if (value) {XtFree(value);} value = NULL;
  
  /* internally 0 rel; externally 1 rel */
  itemp--;
  
  /* check value */
  if ((image[CurImag].valid) && 
      ((itemp<0) || (itemp>=image[CurImag].myDesc->inaxes[2])))
    { /* error */
      MessageShow ("Error: plane number out of range");
      return;}
  
  /* OK, save value */
  image[CurImag].PlaneNo = itemp;
  
} /* end ReadPlaneCB */

/**
 * Callback for read Dimension 4 number
 * \param w           widget activated
 * \param clientData  client data
 * \param callData    call data
 */
void ReadDim4CB (Widget w, XtPointer clientData, XtPointer callData)
{
  /*ImageDisplay *IDdata = (ImageDisplay *)clientData;*/
  char     *value=NULL;
  int      itemp;
  
  /* read value */
  value = XmTextGetString (w);
  if (!value) /* error */
    {MessageShow ("Error reading dimension 4 number");
    return;}
  if (!sscanf (value, "%d", &itemp))
    { /* error */
      MessageShow ("Error reading dimension 4 number");
      if (value) {XtFree(value);} value = NULL;
      return;}
  if (value) {XtFree(value);} value = NULL;
  
  /* internally 0 rel; externally 1 rel */
  itemp--;
  
  /* check value */
  if ((image[CurImag].valid) && 
      ((itemp<0) || (itemp>=image[CurImag].myDesc->inaxes[3])))
    { /* error */
      MessageShow ("Error: dimension 4 out of range");
      return;}
  
  /* OK, save value */
  image[CurImag].hiDim[0] = itemp;
  
} /* end ReadDim4CB */

/**
 * Callback for read Dimension 5 number
 * \param w           widget activated
 * \param clientData  client data
 * \param callData    call data
 */
void ReadDim5CB (Widget w, XtPointer clientData, XtPointer callData)
{
  /*ImageDisplay *IDdata = (ImageDisplay *)clientData;*/
  char     *value=NULL;
  int      itemp;
  
  /* read value */
  value = XmTextGetString (w);
  if (!value) /* error */
    {MessageShow ("Error reading dimension 5 number");
    return;}
  if (!sscanf (value, "%d", &itemp))
    { /* error */
      MessageShow ("Error reading dimension 5 number");
      if (value) {XtFree(value);} value = NULL;
      return;}
  if (value) {XtFree(value);} value = NULL;
  
  /* internally 0 rel; externally 1 rel */
  itemp--;
  
  /* check value */
  if ((image[CurImag].valid) && 
      ((itemp<0) || (itemp>=image[CurImag].myDesc->inaxes[4])))
    { /* error */
      MessageShow ("Error: dimension 5 out of range");
      return;}
  
  /* OK, save value */
  image[CurImag].hiDim[1] = itemp;
  
} /* end ReadDim5CB */

/**
 * Callback for read Dimension 6 number
 * \param w           widget activated
 * \param clientData  client data
 * \param callData    call data
 */
void ReadDim6CB (Widget w, XtPointer clientData, XtPointer callData)
{
  /*ImageDisplay *IDdata = (ImageDisplay *)clientData;*/
  char     *value=NULL;
  int      itemp;
  
  /* read value */
  value = XmTextGetString (w);
  if (!value) /* error */
    {MessageShow ("Error reading dimension 6 number");
    return;}
  if (!sscanf (value, "%d", &itemp))
    { /* error */
      MessageShow ("Error reading dimension 6 number");
      if (value) {XtFree(value);} value = NULL;
      return;}
  if (value) {XtFree(value);} value = NULL;
  
  /* internally 0 rel; externally 1 rel */
  itemp--;
  
  /* check value */
  if ((image[CurImag].valid) && 
      ((itemp<0) || (itemp>=image[CurImag].myDesc->inaxes[5])))
    { /* error */
      MessageShow ("Error: dimension 6 out of range");
      return;}
  
  /* OK, save value */
  image[CurImag].hiDim[2] = itemp;
  
} /* end ReadDim6CB */

/**
 * Callback for mapping function type
 * \param w      widget activated
 * \param which  =0=>linear, 1=>nonlinear, 2=asinh, 3=>hist. eq.
 * \param state  button state
 */
void MapFnCB (Widget w, int which, XmToggleButtonCallbackStruct *state)
{
  if (state->set) image[CurImag].mapFunc = which;
} /* end MapFnCB */

/* button callbacks */
/**
 * Callback for OK button hit 
 * \param w           widget activated
 * \param clientData  client data
 * \param callData    call data
 */
void OptOKButCB (Widget w, XtPointer clientData, XtPointer callData)
{
  
  ReadMinCB (dia.data1, (XtPointer)dia.BoxData, NULL);
  ReadMaxCB (dia.data2, (XtPointer)dia.BoxData, NULL);
  ReadPlaneCB (dia.plane, (XtPointer)dia.BoxData, NULL);
  ReadDim4CB (dia.dim4, (XtPointer)dia.BoxData, NULL);
  ReadDim5CB (dia.dim5, (XtPointer)dia.BoxData, NULL);
  ReadDim6CB (dia.dim6, (XtPointer)dia.BoxData, NULL);
  ReadERTimeoutCB (dia.ERtimeout, (XtPointer)dia.BoxData, NULL);
  
  /* make it disappear but still exist */
  XtPopdown(dia.dialog);
  
} /* end OptOKButCB */

/**
 * Callback for Cancel button hit 
 * \param w           widget activated
 * \param clientData  client data
 * \param callData    call data
 */
void OptCancelButCB (Widget w, XtPointer clientData, XtPointer callData)
{
  /* make it disappear but still exist */
  XtPopdown(dia.dialog);
  
} /* end OptCancelButCB */

/**
 * Callback for Refresh button hit
 * \param w           widget activated
 * \param clientData  client data
 * \param callData    call data
 */
void OptRefreshButCB (Widget w, XtPointer clientData, XtPointer callData)
{
  int scrollx, scrolly;
  gchar szErrMess[120];

  ReadMinCB (dia.data1, (XtPointer)dia.BoxData, NULL);
  ReadMaxCB (dia.data2, (XtPointer)dia.BoxData, NULL);
  ReadPlaneCB (dia.plane, (XtPointer)dia.BoxData, NULL);
  ReadDim4CB (dia.dim4, (XtPointer)dia.BoxData, NULL);
  ReadDim5CB (dia.dim5, (XtPointer)dia.BoxData, NULL);
  ReadDim6CB (dia.dim6, (XtPointer)dia.BoxData, NULL);
  
  /* reload current file info pixarray */
  if (!image[CurImag].valid) return;  /* tests for valid image */
  if (!image[CurImag].FileName) return;  /* tests for valid name */
  if (!image[CurImag].FileName->sp) return;
  if (!image[CurImag].FileName->length) return;
  
  /* save old scroll */
  scrollx = dia.BoxData->scrollx;
  scrolly = dia.BoxData->scrolly; 

  if (Image2Pix (&image[CurImag], dia.BoxData, 1))
    {/* error */
      g_snprintf (szErrMess, 120, "Error reading FITS file = %s", 
	       image[CurImag].FileName->sp);
      MessageShow(szErrMess);
    }
  
  /* reset scroll */
  dia.BoxData->scrollx = scrollx;
  dia.BoxData->scrolly = scrolly; 
  
  /* reset display */
  PaintImage(dia.BoxData);
  
} /* end OptRefreshButCB */

/**
 * Callback for Reload (from external) button hit
 * \param w           widget activated
 * \param clientData  client data
 * \param callData    call data
 */
void OptReloadButCB (Widget w, XtPointer clientData, XtPointer callData)
{
  int scrollx, scrolly;
  gchar szErrMess[120];

  ReadMinCB (dia.data1, (XtPointer)dia.BoxData, NULL);
  ReadMaxCB (dia.data2, (XtPointer)dia.BoxData, NULL);
  ReadPlaneCB (dia.plane, (XtPointer)dia.BoxData, NULL);
  ReadDim4CB (dia.dim4, (XtPointer)dia.BoxData, NULL);
  ReadDim5CB (dia.dim5, (XtPointer)dia.BoxData, NULL);
  ReadDim6CB (dia.dim6, (XtPointer)dia.BoxData, NULL);
  
  /* reload current file info pixarray */
  if (!image[CurImag].valid) return;  /* tests for valid image */
  if (!image[CurImag].FileName) return;  /* tests for valid name */
  if (!image[CurImag].FileName->sp) return;
  if (!image[CurImag].FileName->length) return;

  /* Force reload */
  image[CurImag].reLoad = TRUE;

  /* save old scroll */
  scrollx = dia.BoxData->scrollx;
  scrolly = dia.BoxData->scrolly; 

  if (Image2Pix (&image[CurImag], dia.BoxData, 1))
    {/* error */
      g_snprintf (szErrMess, 120, "Error reading FITS file = %s", 
	       image[CurImag].FileName->sp);
      MessageShow(szErrMess);
    }
  
  /* reset scroll */
  dia.BoxData->scrollx = scrollx;
  dia.BoxData->scrolly = scrolly; 
  
  /* reset display */
  PaintImage(dia.BoxData);
  
} /* end OptReloadButCB */

/**
 * Callback for create dialog box for Options
 * \param parent      parent widget
 * \param clientData  client data
 * \param callData    call data
 */
void OptionBoxCB (Widget parent, XtPointer clientData, XtPointer callData)
{
  Widget form=NULL, label1=NULL, label1a=NULL, label2=NULL, label3=NULL, label4=NULL;
  Widget radio=NULL, sep=NULL;
  Widget OKbutton=NULL, CancelButton=NULL, RefreshButton=NULL, ReloadButton=NULL;
  XmString     label = NULL, minpix = NULL, maxpix = NULL, pixran = NULL, higher = NULL;
  XmString     labelto = NULL, linear = NULL, nlSqrt = NULL, nlAsinh = NULL, histEq=NULL;
  XmString     plalab = NULL, PixelStr = NULL, wierdstring = NULL;
  char         valuestr[61];
  ImageDisplay *IDdata = (ImageDisplay*)clientData;
  int          start;
  Arg          wargs[5]; 
#define OPTIONBOX_WIDTH 205
#define OPTIONBOX_HEIGHT 350
  
  
  /* register IDdata */
  dia.BoxData = IDdata;
  
  /* don't make another one */
  if (OptionBoxActive) {
    if (XtIsRealized (dia.dialog))
      XMapRaised (XtDisplay(dia.dialog), XtWindow(dia.dialog));
    
    /* bring it back where we can see it */
    XtPopup(dia.dialog, XtGrabNonexclusive);
    
    /* put it some place reasonable */
    /*  where is parent? */
    XtVaGetValues (IDdata->shell, XmNx, &dia.xpos, XmNy, &dia.ypos,  NULL);
    dia.ypos += 160;
    if (dia.xpos<0) dia.xpos = 0;
    XMoveWindow (XtDisplay(IDdata->shell), XtWindow(dia.dialog), 
		 dia.xpos, dia.ypos);
    
    /* set values */
    /* range of values in plane - labels are difficult */
    g_snprintf (valuestr, 60, "  %g %g", image[CurImag].minVal, 
	     image[CurImag].maxVal);
    wierdstring = XmStringCreateSimple (valuestr);
    XtSetArg (wargs[0], XmNlabelString, wierdstring);
    XtSetValues (dia.pixran2, wargs, 1);
    if (wierdstring) {XmStringFree(wierdstring);} wierdstring = NULL;
    
    /* range of planes in image */
    if (image[CurImag].myDesc) 
      g_snprintf (valuestr, 60, "Plane no. (1 -  %d)", 
		  image[CurImag].myDesc->inaxes[2]);
    else g_snprintf (valuestr, 60, "Plane no. (1 - ?)");
    wierdstring = XmStringCreateSimple (valuestr);
    XtSetArg (wargs[0], XmNlabelString, wierdstring);
    XtSetValues (dia.planelab, wargs, 1);
    if (wierdstring) {XmStringFree(wierdstring);} wierdstring = NULL;
    
    /* current min */
    g_snprintf (valuestr, 60, "%f", image[CurImag].PixRange[0]);
    XmTextSetString (dia.data1, valuestr);
    
    /* current max */
    g_snprintf (valuestr, 60, "%f", image[CurImag].PixRange[1]);
    XmTextSetString (dia.data2, valuestr);
    
    /* current plane number */
    g_snprintf (valuestr, 60, "%d", image[CurImag].PlaneNo+1);
    XmTextSetString (dia.plane, valuestr);
    g_snprintf (valuestr, 60, "%d", image[CurImag].hiDim[0]+1);
    XmTextSetString (dia.dim4, valuestr);
    g_snprintf (valuestr, 60, "%d", image[CurImag].hiDim[1]+1);
    XmTextSetString (dia.dim5, valuestr);
    g_snprintf (valuestr, 60, "%d", image[CurImag].hiDim[2]+1);
    XmTextSetString (dia.dim6, valuestr);
    
    /* current ER timeout */
    g_snprintf (valuestr, 60, "%g", ERtime_out);
    XmTextSetString (dia.ERtimeout, valuestr);

    return;
  } /* end of update dialog */
  
  label     = XmStringCreateSimple ("Options panel");
  labelto   = XmStringCreateSimple ("Edit timeout (sec) <0 = forever");
  minpix    = XmStringCreateSimple ("Minimum pixel value");
  maxpix    = XmStringCreateSimple ("Maximum pixel value");
  pixran    = XmStringCreateSimple ("Pixel range in plane");
  linear    = XmStringCreateSimple ("linear");
  nlSqrt    = XmStringCreateSimple ("sqrt");
  nlAsinh   = XmStringCreateSimple ("asinh");
  histEq    = XmStringCreateSimple ("Histogram Equalization");
  higher    = XmStringCreateSimple ("Higher dimensions");
  /* mark as active */
  OptionBoxActive = 1;
  
  dia.dialog = XtVaCreatePopupShell ("OptionBox", xmDialogShellWidgetClass, 
				     IDdata->shell, 
				     XmNautoUnmanage, False,
				     XmNwidth,     (int)(OPTIONBOX_WIDTH*sizeFactor),
				     XmNheight,    (int)(OPTIONBOX_HEIGHT*sizeFactor),
				     XmNdeleteResponse, XmDESTROY,
				     XmNfontList, textFontList, // Set font 
				     NULL);
  
  /* make Form widget to stick things on */
  form = XtVaCreateManagedWidget ("OptionForm", xmFormWidgetClass,
				  dia.dialog,
				  XmNautoUnmanage, False,
				  XmNwidth,     (int)(OPTIONBOX_WIDTH*sizeFactor),
				  XmNheight,    (int)(OPTIONBOX_HEIGHT*sizeFactor),
				  XmNx,           0,
				  XmNy,           0,
				  XmNfontList, textFontList, // Set font 
				  NULL);
  
  /* Box label widgets */
  label1 = XtVaCreateManagedWidget ("Label1", xmLabelWidgetClass, 
				    form, 
				    XmNwidth,            (int)(OPTIONBOX_WIDTH*sizeFactor),
				    XmNlabelString,   label,
				    XmNtopAttachment, XmATTACH_FORM,
				    XmNleftAttachment,  XmATTACH_FORM,
				    XmNfontList, textFontList, // Set font 
				    NULL);
  
  /* Edit request timeout  label widgets */
  label1a = XtVaCreateManagedWidget ("Label1a", xmLabelWidgetClass, 
				     form, 
				     XmNwidth,           (int)(OPTIONBOX_WIDTH*sizeFactor),
				     XmNlabelString,   labelto,
				     XmNtopAttachment, XmATTACH_WIDGET,
				     XmNtopWidget,     label1,
				     XmNleftAttachment,  XmATTACH_FORM,
				     XmNfontList, textFontList, // Set font 
				     NULL);
  
  g_snprintf (valuestr, 60, "%f", ERtime_out);
  dia.ERtimeout = XtVaCreateManagedWidget ("OptionTimeO", xmTextFieldWidgetClass, 
					   form, 
					   XmNwidth,     (int)(OPTIONBOX_WIDTH*sizeFactor),
					   XmNvalue,   valuestr,
					   XmNtopAttachment, XmATTACH_WIDGET,
					   XmNtopWidget,     label1a,
					   XmNleftAttachment,  XmATTACH_FORM,
					   XmNfontList, textFontList, // Set font 
					   NULL);
  /* actual pixel range in plane */
  dia.pixran1 = XtVaCreateManagedWidget ("PixelRange1", xmLabelWidgetClass,
					 form,
					 XmNwidth,       (int)(OPTIONBOX_WIDTH*sizeFactor),
					 XmNlabelString,   pixran,
					 XmNtopAttachment, XmATTACH_WIDGET,
					 XmNtopWidget,     dia.ERtimeout,
					 XmNleftAttachment,  XmATTACH_FORM,
					 XmNfontList, textFontList, // Set font 
					 NULL);
  
  g_snprintf (valuestr, 60, "  %g %g", image[CurImag].minVal, 
	   image[CurImag].maxVal);
  PixelStr = XmStringCreateSimple (valuestr);
  dia.pixran2 = XtVaCreateManagedWidget ("PixelRange2", xmLabelWidgetClass, 
					 form, 
					 XmNwidth,            (int)(OPTIONBOX_WIDTH*sizeFactor),
					 XmNlabelString,   PixelStr,
					 XmNtopAttachment, XmATTACH_WIDGET,
					 XmNtopWidget,     dia.pixran1,
					 XmNleftAttachment,  XmATTACH_FORM,
					 XmNfontList, textFontList, // Set font 
					 NULL);
  /* minimum pixrange */
  label2 = XtVaCreateManagedWidget ("OptionLabel2", xmLabelWidgetClass,
				    form,
				    XmNwidth,       (int)(OPTIONBOX_WIDTH*sizeFactor),
				    XmNlabelString,   minpix,
				    XmNtopAttachment, XmATTACH_WIDGET,
				    XmNtopWidget,     dia.pixran2,
				    XmNleftAttachment,  XmATTACH_FORM,
				    XmNfontList, textFontList, // Set font 
				    NULL);
  
  g_snprintf (valuestr, 60, "%f", image[CurImag].PixRange[0]);
  dia.data1 = XtVaCreateManagedWidget ("OptionData1", xmTextFieldWidgetClass, 
				       form, 
				       XmNwidth,       (int)(OPTIONBOX_WIDTH*sizeFactor),
				       XmNvalue,   valuestr,
				       XmNtopAttachment, XmATTACH_WIDGET,
				       XmNtopWidget,     label2,
				       XmNleftAttachment,  XmATTACH_FORM,
				       XmNfontList, textFontList, // Set font 
				       NULL);
  
  /* maximum pixrange */
  label3 = XtVaCreateManagedWidget ("OptionLabel3", xmLabelWidgetClass,
				    form,
				    XmNwidth,           (int)(OPTIONBOX_WIDTH*sizeFactor),
				    XmNlabelString,   maxpix,
				    XmNtopAttachment, XmATTACH_WIDGET,
				    XmNtopWidget,     dia.data1,
				    XmNleftAttachment,  XmATTACH_FORM,
				    XmNfontList, textFontList, // Set font 
				    NULL);
  
  g_snprintf (valuestr, 60, "%f", image[CurImag].PixRange[1]);
  dia.data2 = XtVaCreateManagedWidget ("OptionData2", xmTextFieldWidgetClass, 
				       form, 
				       XmNwidth,           (int)(OPTIONBOX_WIDTH*sizeFactor),
				       XmNvalue,   valuestr,
				       XmNtopAttachment, XmATTACH_WIDGET,
				       XmNtopWidget,     label3,
				       XmNleftAttachment,  XmATTACH_FORM,
				       XmNfontList, textFontList, // Set font 
				       NULL);
  
  /* Plane number */
  if (image[CurImag].myDesc) 
    g_snprintf (valuestr, 60, "Plane no. (1 -  %d)", 
		image[CurImag].myDesc->inaxes[2]);
  else g_snprintf (valuestr, 60, "Plane no. (1 - ?)");
  plalab = XmStringCreateSimple (valuestr);
  dia.planelab = XtVaCreateManagedWidget ("OptionLabel3", xmLabelWidgetClass,
					  form,
					  XmNwidth,           (int)(OPTIONBOX_WIDTH*sizeFactor),
					  XmNlabelString,   plalab,
					  XmNtopAttachment, XmATTACH_WIDGET,
					  XmNtopWidget,     dia.data2,
					  XmNleftAttachment,  XmATTACH_FORM,
					  XmNfontList, textFontList, // Set font 
					  NULL);
  
  g_snprintf (valuestr, 60, "%d", image[CurImag].PlaneNo+1);
  dia.plane = XtVaCreateManagedWidget ("OptionData2", xmTextFieldWidgetClass, 
				       form, 
				       XmNwidth,           (int)(OPTIONBOX_WIDTH*sizeFactor),
				       XmNvalue,   valuestr,
				       XmNtopAttachment, XmATTACH_WIDGET,
				       XmNtopWidget,     dia.planelab,
				       XmNleftAttachment,  XmATTACH_FORM,
				       XmNfontList, textFontList, // Set font 
				       NULL);
  /* Higher dimensions */
  label4 = XtVaCreateManagedWidget ("OptionLabel4", xmLabelWidgetClass,
				    form,
				    XmNwidth,           OPTIONBOX_WIDTH,
				    XmNlabelString,   higher,
				    XmNtopAttachment, XmATTACH_WIDGET,
				    XmNtopWidget,     dia.plane,
				    XmNleftAttachment,  XmATTACH_FORM,
				    XmNfontList, textFontList, // Set font 
				    NULL);
   g_snprintf (valuestr, 60, "%d", image[CurImag].hiDim[0]+1);
   dia.dim4 = XtVaCreateManagedWidget ("OptionDim4", xmTextFieldWidgetClass, 
				       form, 
				       XmNwidth,           (int)((OPTIONBOX_WIDTH/4)*sizeFactor),
				       XmNvalue,   valuestr,
				       XmNtopAttachment, XmATTACH_WIDGET,
				       XmNtopWidget,     label4,
				       XmNleftAttachment,  XmATTACH_FORM,
				       XmNfontList, textFontList, // Set font 
				       NULL);
   g_snprintf (valuestr, 60, "%d", image[CurImag].hiDim[1]+1);
   dia.dim5 = XtVaCreateManagedWidget ("OptionDim5", xmTextFieldWidgetClass, 
				       form, 
				       XmNwidth,           (int)((OPTIONBOX_WIDTH/4)*sizeFactor),
				       XmNvalue,   valuestr,
				       XmNtopAttachment, XmATTACH_WIDGET,
				       XmNtopWidget,     label4,
				       XmNleftAttachment,  XmATTACH_WIDGET,
				       XmNleftWidget,     dia.dim4,
				       XmNfontList, textFontList, // Set font 
 				       NULL);
  g_snprintf (valuestr, 60, "%d", image[CurImag].hiDim[2]+1);
   dia.dim6 = XtVaCreateManagedWidget ("OptionDim5", xmTextFieldWidgetClass, 
				       form, 
				       XmNwidth,           (int)((OPTIONBOX_WIDTH/3)*sizeFactor),
				       XmNvalue,   valuestr,
				       XmNtopAttachment, XmATTACH_WIDGET,
				       XmNtopWidget,     label4,
				       XmNleftAttachment,  XmATTACH_WIDGET,
				       XmNleftWidget,     dia.dim5,
				       XmNfontList, textFontList, // Set font 
				       NULL);
 /* separator */
  sep = XtVaCreateManagedWidget ("sep", xmSeparatorWidgetClass,
				 form, 
				 XmNwidth,          (int)(OPTIONBOX_WIDTH*sizeFactor),
				 XmNtopAttachment, XmATTACH_WIDGET,
				 XmNtopWidget,     dia.dim5,
				 XmNleftAttachment,  XmATTACH_FORM,
				 NULL);
  /* linear/nonlinear/histogram eq. radio buttons */
  start = image[CurImag].mapFunc;
  radio = XmVaCreateSimpleRadioBox(form, "Mapping_type", start, 
				   (XtCallbackProc)MapFnCB,
				   XmNwidth,           (int)(OPTIONBOX_WIDTH*sizeFactor),
				   XmNtopAttachment, XmATTACH_WIDGET,
				   XmNtopWidget,     sep,
				   XmNleftAttachment,  XmATTACH_FORM,
				   XmVaRADIOBUTTON, linear, NULL, NULL, NULL,
				   XmVaRADIOBUTTON, nlSqrt, NULL,NULL,NULL,
				   XmVaRADIOBUTTON, nlAsinh, NULL,NULL,NULL,
				   XmVaRADIOBUTTON, histEq, NULL,NULL,NULL,
				   XmNfontList, textFontList, // Set font 
				   NULL);
  XtManageChild(radio);
  // set fonts on radio buttons with help from Dr. Google:
  WidgetList children;
  Cardinal numChildren;
  int iii;
  // 2. Fetch the automatically created children array from the container
  XtVaGetValues(radio, 
		XmNchildren, &children, 
		XmNnumChildren, &numChildren, 
		NULL);
  // 3. Loop through and apply your fontList to each child button
  for (iii = 0; iii < numChildren; iii++) {
    // XmVaCreateSimpleRadioBox creates Gadgets by default
    if (XtIsSubclass(children[iii], xmToggleButtonGadgetClass)) {
      XtVaSetValues(children[iii], XmNfontList, textFontList, NULL);
    }
  }
  // 4. Manage the container now that children are altered
  XtManageChild(radio);


  /* OK button */
  OKbutton = XtVaCreateManagedWidget (" OK ", xmPushButtonWidgetClass, 
				      form, 
				      XmNbottomAttachment, XmATTACH_FORM,
				      XmNleftAttachment,  XmATTACH_FORM,
				      XmNfontList, textFontList, // Set font 
				      XmNfontList, textFontList, // Set font 
				      NULL);
  XtAddCallback (OKbutton, XmNactivateCallback, OptOKButCB, (XtPointer)IDdata);
  
  /* Cancel button */
  CancelButton = XtVaCreateManagedWidget ("Cancel", xmPushButtonWidgetClass, 
					  form, 
					  XmNbottomAttachment, XmATTACH_FORM,
					  XmNleftAttachment, XmATTACH_WIDGET,
					  XmNleftWidget,     OKbutton,
					  XmNfontList, textFontList, // Set font 
					  NULL);
  XtAddCallback (CancelButton, XmNactivateCallback, OptCancelButCB, 
		 (XtPointer)IDdata);
  
  /* Refresh button */
  RefreshButton = XtVaCreateManagedWidget ("Refresh", xmPushButtonWidgetClass, 
					   form, 
					   XmNbottomAttachment, XmATTACH_FORM,
					   XmNleftAttachment, XmATTACH_WIDGET,
					   XmNleftWidget,     CancelButton,
					   XmNfontList, textFontList, // Set font 
					   NULL);
  XtAddCallback (RefreshButton, XmNactivateCallback, OptRefreshButCB, 
		 (XtPointer)IDdata);
  
  /* Reload button */
  ReloadButton = XtVaCreateManagedWidget ("Reload", xmPushButtonWidgetClass, 
					  form, 
					  XmNbottomAttachment, XmATTACH_FORM,
					  XmNrightAttachment, XmATTACH_FORM,
					  XmNleftAttachment, XmATTACH_WIDGET,
					  XmNleftWidget,     RefreshButton,
					  XmNfontList, textFontList, // Set font 
					  NULL);
  XtAddCallback (ReloadButton, XmNactivateCallback, OptReloadButCB, 
		 (XtPointer)IDdata);
  
  if (label)     {XmStringFree(label);}     label   = NULL;
  if (labelto)   {XmStringFree(labelto);}   labelto = NULL;
  if (minpix)    {XmStringFree(minpix);}    minpix = NULL;
  if (maxpix)    {XmStringFree(maxpix);}    maxpix = NULL;
  if (pixran)    {XmStringFree(pixran);}    pixran = NULL;
  if (linear)    {XmStringFree(linear);}    linear = NULL;
  if (nlSqrt)    {XmStringFree(nlSqrt);}    nlSqrt = NULL;
  if (nlAsinh)   {XmStringFree(nlAsinh);}   nlAsinh= NULL;
  if (histEq)    {XmStringFree(histEq);}    histEq = NULL;
  if (plalab)    {XmStringFree(plalab);}    plalab = NULL;
  if (PixelStr)  {XmStringFree(PixelStr);}  PixelStr = NULL;
  if (higher)    {XmStringFree(higher);}    higher = NULL;
  
  /* set it up */
  XtManageChild (dia.dialog);
  
  /* put it some place reasonable */
  /*  where is parent? */
  XtVaGetValues (IDdata->shell,
		 XmNx, &dia.xpos,
		 XmNy, &dia.ypos,
		 NULL);
  dia.ypos += 170;
  dia.xpos -= 10;
  if (dia.xpos<0) dia.xpos = 0;
  XMoveWindow (XtDisplay(IDdata->shell), XtWindow(dia.dialog), 
	       dia.xpos, dia.ypos);
} /* end OptionBox */

