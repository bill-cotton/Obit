/* $Id$  */
/* image blink control box  for ObitView */
/*-----------------------------------------------------------------------
*  Copyright (C) 1996,2002-2011
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
#include <Xm/MessageB.h>
#include <Xm/TextF.h>
#include <Xm/Text.h>
#include <time.h>
#include <math.h>
#include "imagedisp.h"
#include "color.h"
#include <unistd.h>

/**
 *  \file blinkbox.c
 * displays "blink" dialog.
 */

/*--------------- file global data ----------------*/
/** is the blink box active? */
int BlinkBoxActive = 0;

/* global structures for things to talk to each other */
typedef struct {
  ImageDisplay *BoxData;   /* image display structure pointer */
  int    oldImag;          /* CurImag at start */
  int    Stop;             /* true if blinking terminated */
  float  Dwell;            /* dwell time on each image */
  Widget dialog;           /* box, start, end planes */
  Widget DwellScroll;      /* Dwell time scroll */
  Widget DwellScrollLabel; /* label for scroll */
  Widget QuitButton;       /* quit button window */
  XtIntervalId timerId;    /* X timer ID for timeout */
} BlinkBoxStuff;

typedef struct { /* define DisplayInfo structure */
  short          value[2];  /* Brightness and contrast scroll bar values */
  unsigned short red[MAXCOLOR], green[MAXCOLOR], blue[MAXCOLOR]; /* colors */
  int            zoom;      /* zoom factor, neg = zoom out */
  int            scrollx;   /* center "x" pixel in display */
  int            scrolly;   /* center "y" pixel in display */
}  BlinkDisplayInfo;

/* file global variables */
BlinkBoxStuff BlinkDia;
BlinkDisplayInfo info[2];
int validInfo[2]={0,0};

/*---------------Private function prototypes----------------*/
void SaveBlinkInfo(ImageDisplay *IDdata);
void RestoreBlinkInfo(ImageDisplay *IDdata);
void AdjustBlinkInfo(void);
int BlinkAreAligned (void);


/**
 * Searches event queue for click of Quit button 
 * \return true if event found
 * THIS IS PROBABLY A BAD THING TO DO - Seems to cause crashes
 */
Boolean BlinkCheckQuit()
{
  XButtonPressedEvent QuitEvent;
  Boolean ret=0;
  
  if (XCheckTypedEvent(XtDisplay(BlinkDia.dialog), ButtonPress, 
		       (XEvent*)&QuitEvent))
    {
      if (QuitEvent.window==XtWindow(BlinkDia.QuitButton)) { /* found it */
      	BlinkDia.Stop = 1; ret=1;
	if (BlinkDia.timerId!=0) XtRemoveTimeOut(BlinkDia.timerId);
	BlinkDia.timerId = 0;
	usleep(250000);  /* Give things a chance to settle - 250 msec */
     }  else {/* it's not yours - put it back */
      	XPutBackEvent (XtDisplay(BlinkDia.dialog), 
		       (XEvent*)&QuitEvent); return 0;
      }
    }
  return ret;
} /* end BlinkCheckQuit */

/**
 * Tell dwell time, writes in box
 */
void BlinkLabelDwell(void)
{
  char      string[30];
  XmString  WierdString = NULL;
  
  /* value label */
  sprintf (string,"    Dwell Time %4.1f sec.", BlinkDia.Dwell);
  WierdString = XmStringCreateSimple (string);
  XtVaSetValues(BlinkDia.DwellScrollLabel, 
		XmNlabelString,   WierdString,
		NULL);
  if (WierdString) XmStringFree(WierdString); WierdString = NULL;
} /* end BlinkLabelDwell */

/**
 * Shutdown blink operation and restore to previous condition 
 */
void BlinkShutDown (void)
{
  /* restore original condition */
  CurImag = BlinkDia.oldImag; /* reset original image */
  RestoreBlinkInfo(BlinkDia.BoxData); /* reset display */
  /* show current image */
  ResetDisplay(BlinkDia.BoxData);
  /* reset color table */
  SetColorTable (BlinkDia.BoxData);
  
  /* I could just die */
  BlinkBoxActive = 0; /* mark as inactive */
  XtDestroyWidget (BlinkDia.dialog);  /* your wish is granted */

} /* end BlinkShutDown */


/**
 * Callback for scroll bar used to select plane
 * \param w           scroll widget activated
 * \param clientData  client data
 * \param callData    call data
 */
void BlinkScrollDwellCB (Widget w, XtPointer clientData, XtPointer callData)
{
  /*ImageDisplay *IDdata = (ImageDisplay *)clientData;*/
  XmScaleCallbackStruct *call_data = (XmScaleCallbackStruct *) callData;
  usleep(250000);  /* Give things a chance to settle - 250 msec */
  /* read value of scrollbar */
  BlinkDia.Dwell = call_data->value * 0.1; /* 10th of sec */
  /* update label */
  BlinkLabelDwell();
} /* end BlinkScrollPlaneCB */

/**
 * Callback for display next image and set next call to itself until quit hit
 * \param w           scroll widget activated
 * \param clientData  client data
 * \param callData    call data
 */
void BlinkPlay (Widget w, XtPointer clientData, XtPointer callData)
{
  unsigned long interval;
  
  /* Done? */
  /*if (BlinkDia.Stop || BlinkCheckQuit()) quit button? */
  if (BlinkDia.Stop) /* quit button? */
    { usleep(250000);  /* Give things a chance to settle - 250 msec */
      /*fprintf (stderr, "Stop blink\n");*/
      BlinkShutDown();
    return;}
  
  /* toggle image */
  if (CurImag) CurImag = 0;
  else CurImag = 1;
  
  /* start timer */
  interval = 1000.0 * BlinkDia.Dwell;
  BlinkDia.timerId = XtAppAddTimeOut(BlinkDia.BoxData->app, interval, 
				     (XtTimerCallbackProc)BlinkPlay, NULL);
  
  BlinkDia.BoxData->showInfo = 0; /* turn off display of pixel information */
  /* show other image */
  ResetDisplay(BlinkDia.BoxData);
  /* reset display info */
  RestoreBlinkInfo(BlinkDia.BoxData); 
  /* sync events - what does this do other than crash */
  XSync (XtDisplay(BlinkDia.dialog), False);
  
  /* reset color table */
  SetColorTable (BlinkDia.BoxData);
} /* end BlinkPlay */

/**
 * Callback for quit button
 * \param w           widget activated
 * \param clientData  client data
 * \param callData    call data
 */
void BlinkQuitButCB (Widget w, XtPointer clientData, XtPointer callData)
{
  /* just mark as stopped and BlinkPlay will shut down */
  BlinkDia.Stop = 1;
  
} /* end BlinkQuitButCB */

/**
 * Callback for swap current and "Blink" store images
 * \param w           widget activated
 * \param clientData  client data
 * \param callData    call data
 */
void BlinkSwapCB (Widget parent, XtPointer clientData, XtPointer callData)
     /* client data is ImageDisplay pointer */
{
  ImageDisplay *IDdata = (ImageDisplay*)clientData;
  
  if (!IDdata) return; /* validity check */
  
  /* register IDdata */
  BlinkDia.BoxData = IDdata;
  
  IDdata->showInfo = 0; /* turn off display of pixel information */
  SaveBlinkInfo(IDdata); /* save display info */
  if (CurImag==0) CurImag = 1; /* swap images */
  else CurImag = 0;
  if (validInfo[CurImag]) RestoreBlinkInfo(IDdata);
  else /* default setup */
    {
      IDdata->zoom = 1;
      IDdata->scrollx = 0;
      IDdata->scrolly = 0;
      IDdata->iXCorn = 0;
      IDdata->iYCorn = 0;
      IDdata->value[0] = 128;
      IDdata->value[1] = 128;
      InitColorTable(IDdata);
    }
  /* show other image */
  ResetDisplay(IDdata);
  /* reset color table */
  SetColorTable (IDdata);
} /* end BlinkSwapCB */

/**
 * Callback for create dialog box for "Blink" display of planes
 *  client data is ImageDisplay pointer
 * \param w           widget activated
 * \param clientData  client data
 * \param callData    call data
 */
void BlinkBlinkCB (Widget parent, XtPointer clientData, XtPointer callData)
{
  Widget       form;
  XmString     WierdString = NULL;
  char         valuestr[30];
  ImageDisplay *IDdata = (ImageDisplay*)clientData;
  short        xpos, ypos;
  
  
  /* register IDdata */
  BlinkDia.BoxData = IDdata;
  
  /* don't make another one */
  if (BlinkBoxActive) {
    if (XtIsRealized (BlinkDia.dialog))
      XMapRaised (XtDisplay(IDdata->shell), XtWindow(BlinkDia.dialog));
    return;
  }
  
  /* mark as active */
  BlinkBoxActive = 1;
  
  /* other initialization */
  BlinkDia.Dwell = 2.0;
  BlinkDia.Stop = 0;
  BlinkDia.timerId = 0;
 
  /* save display states */
  SaveBlinkInfo(IDdata); /* save display info */
  AdjustBlinkInfo(); /* adjust zoom scroll if they are the same on the sky */
  BlinkDia.oldImag = CurImag;
  
  BlinkDia.dialog = XtVaCreatePopupShell ("BlinkBox",
					  xmDialogShellWidgetClass, 
					  IDdata->shell, 
					  XmNautoUnmanage, False,
					  XmNwidth,     200,
					  XmNheight,     80,
					  XmNdeleteResponse, XmDESTROY,
					  NULL);
  
  /* make Form widget to stick things on */
  form = XtVaCreateManagedWidget ("BlinkForm", xmFormWidgetClass,
				  BlinkDia.dialog,
				  XmNautoUnmanage, False,
				  XmNwidth,     200,
				  XmNheight,    80,
				  XmNx,           0,
				  XmNy,           0,
				  NULL);
  
  /* plane scroll */
  BlinkDia.DwellScroll = XtVaCreateManagedWidget ("BlinkDwellScroll", 
						  xmScaleWidgetClass, 
						  form,
						  XmNwidth,            200,
						  XmNheight,           15,
						  XmNmaximum,          100,
						  XmNminimum,           10,
						  XmNvalue,             20,
						  XmNshowValue,       True,
						  XmNorientation, XmHORIZONTAL,
						  XmNprocessingDirection, XmMAX_ON_RIGHT,
						  XmNleftAttachment,  XmATTACH_FORM,
						  XmNtopAttachment, XmATTACH_FORM,
						  NULL);
  
  XtAddCallback(BlinkDia.DwellScroll, XmNvalueChangedCallback, 
                BlinkScrollDwellCB, (XtPointer)IDdata);
  
  /* scroll label */
  sprintf (valuestr, "dummy string");
  WierdString = XmStringCreateSimple (valuestr);
  BlinkDia.DwellScrollLabel = XtVaCreateManagedWidget ("BlinkScrollLabel", 
						       xmLabelWidgetClass, 
						       form, 
						       XmNwidth,           200,
						       XmNtopAttachment, XmATTACH_WIDGET,
						       XmNtopWidget,     BlinkDia.DwellScroll,
						       XmNleftAttachment,  XmATTACH_FORM,
						       NULL);
  if (WierdString) XmStringFree(WierdString); WierdString = NULL;
  BlinkLabelDwell();
  
  /* Quit button */
  BlinkDia.QuitButton = 
    XtVaCreateManagedWidget ("Quit", xmPushButtonWidgetClass, 
			     form, 
			     XmNbottomAttachment, XmATTACH_FORM,
			     XmNleftAttachment,  XmATTACH_FORM,
			     XmNrightAttachment, XmATTACH_FORM,
			     NULL);
  XtAddCallback (BlinkDia.QuitButton, XmNactivateCallback, BlinkQuitButCB, 
		 (XtPointer)IDdata);
  
  /* put it some place reasonable */
  /*  where is parent? */
  XtVaGetValues (IDdata->shell,
		 XmNx, &xpos,
		 XmNy, &ypos,
		 NULL);
  xpos -= 50;
  if (xpos<0) xpos = 0;
  ypos += 160;
  XMoveWindow (XtDisplay(IDdata->shell), XtWindow(BlinkDia.dialog),
	       xpos, ypos);
  
  /* set it up */
  XtManageChild (BlinkDia.dialog);
  
  /* start it going */
  BlinkPlay (BlinkDia.dialog, clientData, callData);
} /* end BlinkBox */

/**
 * Save image display info
 * \param IDdata   Image display
 */
void SaveBlinkInfo(ImageDisplay *IDdata)
{
  int i;
  if (!IDdata) return; /* validity check */
  
  validInfo[CurImag] = 1;
  info[CurImag].value[0] = IDdata->value[0];
  info[CurImag].value[1] = IDdata->value[1];
  if (IDdata->zoom==0) IDdata->zoom = 1;
  info[CurImag].zoom = IDdata->zoom;
  info[CurImag].scrollx = IDdata->scrollx;
  info[CurImag].scrolly = IDdata->scrolly;
  for (i=0;i<MAXCOLOR;i++) {
    info[CurImag].red[i] = IDdata->red[i];
    info[CurImag].blue[i] = IDdata->blue[i];
    info[CurImag].green[i] = IDdata->green[i];}
} /* end SaveBlinkInfo */

/**
 * Restore image display info
 * \param IDdata   Image display
 */
void RestoreBlinkInfo(ImageDisplay *IDdata)
{
  int i;
  if (!IDdata) return; /* validity check */
  
  IDdata->value[0] = info[CurImag].value[0];
  IDdata->value[1] = info[CurImag].value[1];
  if (info[CurImag].zoom==0) info[CurImag].zoom = 1;
  IDdata->zoom = info[CurImag].zoom;
  IDdata->scrollx = info[CurImag].scrollx;
  IDdata->scrolly = info[CurImag].scrolly;
  for (i=0;i<MAXCOLOR;i++) {
    IDdata->red[i] = info[CurImag].red[i];
    IDdata->blue[i] = info[CurImag].blue[i];
    IDdata->green[i] = info[CurImag].green[i];}
} /* end RestoreBlinkInfo */

/**
 * Set zoom and scroll if the two blink images are the same part of the sky
 */
void AdjustBlinkInfo(void)
{ 
  int other;
  
  if (!BlinkAreAligned()) return;
  if (CurImag) other = 0;
  else other = 1;
  info[other].zoom = info[CurImag].zoom;
  info[other].scrollx = info[CurImag].scrollx;
  info[other].scrolly = info[CurImag].scrolly;
}  /* end AdjustBlinkInfo */

/**
 * Routine to determine if the two blink images have coincident pixels
 * \return 1 if so else 0
 */
olong BlinkAreAligned (void)
{  
  olong same;
  odouble arg1, arg2;
  ObitImageDesc *desc1=image[0].myDesc, *desc2=image[1].myDesc;
  
  if (desc1==NULL) return 0;
  if (desc2==NULL) return 0;

  /* this is easy if they are the same file.  */
  same = !strcmp(image[0].FileName->sp, image[1].FileName->sp);
  if (same) return 1;
  if (desc1->inaxes[0]!=desc2->inaxes[0]) return 0;
  if (desc1->inaxes[1]!=desc2->inaxes[1]) return 0;
  arg1 = (desc1->crpix[0]-desc2->crpix[0]);
  arg2 = (desc2->crpix[0]);
  if (fabs(arg1) > 0.1 * fabs(arg2)) return 0;
  arg1 = (desc1->crpix[1]-desc2->crpix[1]);
  arg2 = (desc2->crpix[1]);
  if (fabs(arg1) > 0.1 * fabs(arg2)) return 0;
  if (fabs(desc1->crval[0]-desc2->crval[0]) >
      0.1 * fabs(desc2->crpix[0])) return 0;
  if (fabs(desc1->crval[1]-desc2->crval[1]) >
      0.1 * fabs(desc2->crpix[1])) return 0;
  arg1 = (desc1->crota[1]-desc2->crota[1]);
  if (fabs(arg1) > 1.0) return 0;
  if (strcmp(desc1->ctype[0],desc2->ctype[0])) return 0;
  if (strcmp(desc1->ctype[1],desc2->ctype[1])) return 0;

  return 1;  /* Seems OK  */
} /* end BlinkAreAligned  */

