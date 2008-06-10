/* $Id: moviebox.c,v 1.8 2007/01/11 19:12:34 bcotton Exp $  */
/* Movie control box  for ObitView */
/*-----------------------------------------------------------------------
*  Copyright (C) 1996,1997,1999, 2002-2008-2008
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
#include <Xm/ScrollBar.h>
#include <Xm/Form.h>
#include <Xm/PushB.h>
#include <Xm/Label.h>
#include <Xm/ToggleB.h>
#include <Xm/RowColumn.h>
#include <Xm/Separator.h>
#include <Xm/MessageB.h>
#include <Xm/TextF.h>
#include <Xm/Text.h>
#include <time.h>
#include "imagedisp.h"
#include "Image2Pix.h"
#include "poslabel.h"
#include "messagebox.h"

/**
 *  \file moviebox.c
 * displays "movie" dialog.
 */

/*--------------- file global data ----------------*/
/** is the movie box active? */
int MovieBoxActive = 0;

/* global structure for things to talk to each other */
typedef struct {
  ImageDisplay *BoxData;
  int    StartPlane, EndPlane, CurPlane;
  int    Stop;                /* if true stop button hit */
  float  Dwell;
  Widget dialog, start, end;  /* box, start, end planes */
  Widget dwell;               /* Dwell time */
  Widget PlaneScroll;         /* plane select/display scroll */
  Widget PlaneScrollLabel;    /* label for plane number */
  Widget PlaneLabel;          /* label for plane axis value */
  Widget StopButton;          /* stop button */
  Window window;              /* stop button window */
  short  xpos, ypos;          /* location of the box */
} MovieBoxStuff;
MovieBoxStuff MovieDia;

/**
 * Searches event queue for click of Stop button 
 * \return true if event found
 */
Boolean MovieCheckStop()
{
  XButtonPressedEvent *StopEvent;
  XEvent              evt;
  
  if (XCheckTypedEvent(XtDisplay(MovieDia.dialog), ButtonPress, &evt))
    { StopEvent = (XButtonPressedEvent*)&evt;
    if (StopEvent->window==MovieDia.window) { /* found it */
      MovieDia.Stop = 1; 
      /* clear out queue */
      while (XPending(XtDisplay(MovieDia.dialog)))
	XNextEvent (XtDisplay(MovieDia.dialog), &evt);
      XSync (XtDisplay(MovieDia.dialog), True);
      return 1;}
    else /* it's not yours - put it back */
      {XPutBackEvent (XtDisplay(MovieDia.dialog), &evt); 
      return 0;}
    }
  return 0;
} /* end MovieCheckStop */

/**
 * Callback for get first plane for movie
 * \param w           widget activated
 * \param clientData  client data
 * \param callData    call data
 */
void ReadStartPlaneCB (Widget w, XtPointer clientData, XtPointer callData)
{
  /*ImageDisplay *IDdata = (ImageDisplay *)clientData;*/
  char     *value=NULL;
  int    itemp;
  
  /* read value */
  value = XmTextGetString (w);
  if (!value) /* error */
    {MessageShow ("Error reading start plane number");
    MovieDia.Stop = 1; return;}
  if (!sscanf (value, "%d", &itemp))
    { /* error */
      MessageShow("Error reading start plane number");
      if (value) XtFree(value); value = NULL;
      MovieDia.Stop = 1; return;}
  if (value) XtFree(value); value = NULL;
  
  /* internally 0 rel; externally 1 rel */
  itemp--;
  
  /* check value */
  if ((image[CurImag].valid) && 
      ((itemp<0) || (itemp>=image[CurImag].myDesc->inaxes[2])))
    { /* error */
      MessageShow ("Error: start plane number out of range");
      MovieDia.Stop = 1; return;}
  /* OK, save */
  MovieDia.StartPlane = itemp;
} /* end ReadStartPlaneCB */

/**
 * Callback for get highest plane for movie
 * \param w           widget activated
 * \param clientData  client data
 * \param callData    call data
 */
void ReadEndPlaneCB (Widget w, XtPointer clientData, XtPointer callData)
{
  /*ImageDisplay *IDdata = (ImageDisplay *)clientData;*/
  char     *value=NULL;
  int    itemp;
  
  /* read value */
  value = XmTextGetString (w);
  if (!value) /* error */
    {MessageShow("Error reading end plane number");
    MovieDia.Stop = 1; return;}
  if (!sscanf (value, "%d", &itemp))
    { /* error */
      MessageShow("Error reading end plane number");
      if (value) XtFree(value); value = NULL;
      MovieDia.Stop = 1; return;}
  if (value) XtFree(value); value = NULL;
  
  /* internally 0 rel; externally 1 rel */
  itemp--;
  
  /* check value */
  if ((image[CurImag].valid) && 
      ((itemp<0) || (itemp>=image[CurImag].myDesc->inaxes[2])))
    { /* error */
      MessageShow ("Error: end plane number out of range");
      MovieDia.Stop = 1; return;}
  /* OK, save */
  MovieDia.EndPlane = itemp;
} /* end ReadEndPlaneCB */

/**
 * Callback for get Dwell time
 * \param w           scroll widget activated
 * \param clientData  client data
 * \param callData    call data
 */
void ReadDwellCB (Widget w, XtPointer clientData, XtPointer callData)
{
  /*ImageDisplay *IDdata = (ImageDisplay *)clientData;*/
  char     *value;
  float    temp;
  
  /* read value */
  value = XmTextGetString (w);
  if (!value) /* error */
    {MessageShow ("Error Dwell time");
    MovieDia.Stop = 1; return;}
  if (!sscanf (value, "%e", &temp))
    { /* error */
      MessageShow ("Error Dwell time");
      if (value) XtFree(value); value = NULL;
      MovieDia.Stop = 1; return;}
  if (value) XtFree(value); value = NULL;
  if (MovieDia.Dwell<0.0) MovieDia.Dwell = 0.0;  /* constrain */
  if (MovieDia.Dwell>100.0) MovieDia.Dwell = 100.0;
  /* OK, save */
  MovieDia.Dwell = temp;
} /* end ReadDwellCB */

/**
 * Label Movie current plane
 */
void MovieLabelPlane(void)
{
  gchar axtype[9], label[20], string[30];
  odouble value;
  XmString     WierdString;
  
  /* valid image to load? */
  if (!image[CurImag].valid) return;

  /* value label - the descriptor was subimaged so that the 
   "current" plane is 1 */
  strncpy(axtype, image[CurImag].myDesc->ctype[2], 8); 
  axtype[8] = 0;
  value = image[CurImag].myDesc->crval[2] + 
    (1.0 - image[CurImag].myDesc->crpix[2]) * 
    image[CurImag].myDesc->cdelt[2];  /* assume linear */
  AxisLabel(value, axtype, label);  /* human readable */
  sprintf (string,"    %s",label);
  WierdString = XmStringCreateSimple (string);
  XtVaSetValues(MovieDia.PlaneLabel, 
		XmNlabelString,   WierdString,
		NULL);
  if (WierdString) XmStringFree(WierdString); WierdString = NULL;
} /* end MovieLabelPlane */

/**
 * Load Movie current plane CurPlane 
 */
void MovieLoadPlane (void)
{
  XmString     WierdString;
  gchar         valuestr[30], szErrMess[120];
  olong         xscroll, yscroll;  

  /* valid image to load? */
  if (!image[CurImag].valid) return;

  /* make sure CurPlane in range */
  if (MovieDia.EndPlane<MovieDia.StartPlane) 
    MovieDia.EndPlane = MovieDia.StartPlane;
  if (MovieDia.CurPlane>=image[CurImag].myDesc->inaxes[2]) 
    MovieDia.CurPlane = image[CurImag].myDesc->inaxes[2]-1;
  if (MovieDia.CurPlane<0)  MovieDia.CurPlane = 0;
  
  /* load CurPlane */
  if (!image[CurImag].valid) return;  /* tests for valid image */
  if (!image[CurImag].FileName) return;  /* tests for valid name */
  if (!image[CurImag].FileName->sp) return;
  if (!image[CurImag].FileName->length) return;
  image[CurImag].PlaneNo = MovieDia.CurPlane;   /* Which Plane */
  image[CurImag].reLoad  = TRUE;                /* Force load */
 
  /* save current scroll info which Image2Pix will reset */
  xscroll = MovieDia.BoxData->scrollx;
  yscroll = MovieDia.BoxData->scrolly;
  
  if (Image2Pix (&image[CurImag], MovieDia.BoxData, 0))
    {/* error */
      sprintf (szErrMess, "Error reading FITS file = %s", 
	       image[CurImag].FileName->sp);
      MessageShow (szErrMess); 
      image[CurImag].PlaneNo = 0;  /* Can't read this plane */
      return;
    }
  
  /* restore current scroll info*/
  MovieDia.BoxData->scrollx = xscroll;
  MovieDia.BoxData->scrolly = yscroll;
  
  /*  display */
  PaintImage(MovieDia.BoxData);
  
  /* update labels */
  XtVaSetValues(MovieDia.PlaneScroll, 
		XmNvalue,       MovieDia.CurPlane+1,
		NULL);
  sprintf (valuestr, "    Plane no. %d", MovieDia.CurPlane+1);
  WierdString = XmStringCreateSimple (valuestr);
  XtVaSetValues(MovieDia.PlaneScrollLabel, 
		XmNlabelString,   WierdString,
		NULL);
  if (WierdString) XmStringFree(WierdString); WierdString = NULL;
  MovieLabelPlane();
} /* end MovieLoadPlane */

void MovieScrollPlaneCB (Widget w, XtPointer clientData, XtPointer callData)
{
  /* scroll bar used to select plane */
  /*ImageDisplay *IDdata = (ImageDisplay *)clientData;*/
  XmScrollBarCallbackStruct *call_data = 
    (XmScrollBarCallbackStruct *) callData;
  
  /* read value of scrollbar */
  MovieDia.CurPlane = call_data->value-1; /* 0 rel */
  
  /* Show selected plane */
  MovieLoadPlane();
} /* end MovieScrollPlaneCB */

/**
 *  Display next plane and set next call
 * \param w           widget activated
 * \param clientData  client data
 * \param callData    call data
 */
void MoviePlay (Widget w, XtPointer clientData, XtPointer callData)
{
  unsigned long interval;
  XtIntervalId  timer;
  
  /* valid image to load? */
  if (!image[CurImag].valid) return;

  /* Done? */
  if (MovieDia.Stop) return; /* stop button */
  if (MovieDia.CurPlane>=MovieDia.EndPlane) 
    /*    {MovieDia.Stop = 1; return;}*/ /* done ?*/
    /* just continue looping */
    {MovieDia.CurPlane=MovieDia.StartPlane-1;}
  if (MovieCheckStop()) return;
  MovieDia.CurPlane++; /* next plane */
  
  /* let server catch up */
  /*  XSynchronize(XtDisplay(MovieDia.dialog), True);*/
  
  /* start timer */
  interval = 1000.0 * MovieDia.Dwell;
  timer = XtAppAddTimeOut(MovieDia.BoxData->app, interval, 
			  (XtTimerCallbackProc)MoviePlay, NULL);
  /*  fprintf (stderr,"MoviePlay: events pending = %d\n",
      XPending(XtDisplay(MovieDia.dialog)));*/
  /* show next plane */
  MovieLoadPlane();
  
} /* end MoviePlay */

/* button callbacks */
/**
 * Callback for Play button hit
 * \param w           widget activated
 * \param clientData  client data
 * \param callData    call data
 */
void MoviePlayButCB (Widget w, XtPointer clientData, XtPointer callData)
{
  Display *dpy = XtDisplay(w);
  
  /* ignore if it's already running */
  if (MovieDia.Stop==0) return;
  
  MovieDia.Stop = 0; /* indicate movie in progress */
  
  /* get instructions */
  ReadStartPlaneCB (MovieDia.start, (XtPointer)MovieDia.BoxData, NULL);
  ReadEndPlaneCB (MovieDia.end, (XtPointer)MovieDia.BoxData, NULL);
  ReadDwellCB (MovieDia.dwell, (XtPointer)MovieDia.BoxData, NULL);
  
  /* see if any of the read routines has turned off the movie */
  if (MovieDia.Stop) {
    XBell(dpy, 50); 
    return;}
  
  /* start movie here */
  MovieDia.CurPlane = MovieDia.StartPlane-1;
  
  MoviePlay (w, clientData, callData); /* first plane */
} /* end MoviePlayButCB */

/**
 * Callback for Quit button hit
 * \param w           widget activated
 * \param clientData  client data
 * \param callData    call data
 */
void MovieQuitButCB (Widget w, XtPointer clientData, XtPointer callData)
{
  /* if it's running just stop the movie */
  if (MovieDia.Stop==0) 
    {MovieDia.Stop = 1; return;}
  
  /* make it disappear but still exist
  XtPopdown(MovieDia.dialog); */

  /* Better to delete as tne next image may have a different  no.planes */
  XtDestroyWidget(MovieDia.dialog);   
  MovieDia.dialog = NULL;
  MovieBoxActive = 0;
} /* end MovieQuitButCB */

void MovieStopButCB (Widget w, XtPointer clientData, XtPointer callData)
     /* Stop button hit */
{
  MovieDia.Stop = 1; /* leave marker */
  
} /* end MovieStopButCB */

/**
 * Callback for create dialog box for "Movie" display of planes 
 * \param w           widget activated
 * \param clientData  ImageDisplay pointer
 * \param callData    call data
 */
void MovieBoxCB (Widget parent, XtPointer clientData, XtPointer callData)
{
  Widget form, label1, label2, label3;
  Widget PlayButton, QuitButton;
  XmString     StartLab = NULL;
  XmString     EndLab   = NULL;
  XmString     DwellLab = NULL;
  XmString     WierdString = NULL;
  char         valuestr[30];
  ImageDisplay *IDdata = (ImageDisplay*)clientData;
  int          hiPlane;
  int          value, increment, slider_size, page_increment;
  
  
  /* register IDdata */
  MovieDia.BoxData = IDdata;
  
  /* don't make another one */
  if (MovieBoxActive) {
    if (XtIsRealized (MovieDia.dialog))
      XMapRaised (XtDisplay(IDdata->shell), XtWindow(MovieDia.dialog));
    
    /* bring it back where we can see it */
    XtPopup(MovieDia.dialog, XtGrabNonexclusive);
    
    /* put it some place reasonable */
    /*  where is parent? */
    XtVaGetValues (IDdata->shell, XmNx, &MovieDia.xpos, XmNy, &MovieDia.ypos,
		   NULL);
    MovieDia.ypos += 160;
    if (MovieDia.xpos<0) MovieDia.xpos = 0;
    XMoveWindow (XtDisplay(IDdata->shell), XtWindow(MovieDia.dialog), 
		 MovieDia.xpos, MovieDia.ypos);
    
    /* set values */
    MovieDia.EndPlane = image[CurImag].myDesc->inaxes[2];
    MovieDia.StartPlane = 1;
    MovieDia.CurPlane = 0;
    
    /* start Plane */
    sprintf (valuestr, "%d", MovieDia.StartPlane);
    XmTextSetString (MovieDia.start, valuestr);
    
    /* end Plane */
    sprintf (valuestr, "%d", MovieDia.EndPlane);
    XmTextSetString (MovieDia.end, valuestr);
    
    /* Dwell */
    sprintf (valuestr, "%f", MovieDia.Dwell);
    XmTextSetString (MovieDia.dwell, valuestr);
    
    /* update labels */
    sprintf (valuestr, "    Plane no. %d", MovieDia.CurPlane+1);
    WierdString = XmStringCreateSimple (valuestr);
    XtVaSetValues(MovieDia.PlaneScrollLabel, 
		  XmNlabelString,   WierdString,
		  NULL);
    if (WierdString) XmStringFree(WierdString); WierdString = NULL;
    MovieLabelPlane();
    return;
    
    /* plane scroll */
    hiPlane = MovieDia.EndPlane+1;
    slider_size = hiPlane/10; if (slider_size<1) slider_size = 1;
    value = MovieDia.CurPlane+1;
    increment = 1;
    page_increment = hiPlane/10; if (page_increment<1) page_increment = 1;
    XtVaSetValues(MovieDia.PlaneScroll, 
		  XmNsliderSize,  slider_size,
		  XmNvalue,       MovieDia.CurPlane+1,
		  XmNmaximum,     MovieDia.EndPlane+1+slider_size,
		  NULL);
    XmScrollBarSetValues (MovieDia.PlaneScroll, value, slider_size, 
			  increment, page_increment, True);
    
    XtAddCallback(MovieDia.PlaneScroll, XmNvalueChangedCallback, 
		  MovieScrollPlaneCB, (XtPointer)IDdata);
    
    /* set it up */
    XtManageChild (MovieDia.dialog);
    
  } /* end of update box */
  
  StartLab = XmStringCreateSimple ("Planes");
  EndLab   = XmStringCreateSimple ("to");
  DwellLab = XmStringCreateSimple ("Dwell (sec)");
  /* mark as active */
  MovieBoxActive = 1;
  
  
  /* other initialization */
  MovieDia.CurPlane = image[CurImag].PlaneNo;
  MovieDia.StartPlane = 1;
  MovieDia.EndPlane = image[CurImag].myDesc->inaxes[2];
  if (MovieDia.EndPlane<MovieDia.StartPlane) 
    MovieDia.EndPlane = MovieDia.StartPlane;
  MovieDia.Dwell = 1.0;
  MovieDia.Stop = 1;
  hiPlane = MovieDia.EndPlane+1;
  if (hiPlane == MovieDia.StartPlane) hiPlane++;
  
  MovieDia.dialog = XtVaCreatePopupShell ("MovieBox",
					  xmDialogShellWidgetClass, 
					  IDdata->shell, 
					  XmNautoUnmanage, False,
					  XmNwidth,     150,
					  XmNheight,    230,
					  XmNdeleteResponse, XmDESTROY,
					  NULL);
  
  /* make Form widget to stick things on */
  form = XtVaCreateManagedWidget ("MovieForm", xmFormWidgetClass,
				  MovieDia.dialog,
				  XmNautoUnmanage, False,
				  XmNwidth,     150,
				  XmNheight,    230,
				  XmNx,           0,
				  XmNy,           0,
				  NULL);
  
  /* Play button */
  PlayButton = XtVaCreateManagedWidget (" Play ", xmPushButtonWidgetClass, 
					form, 
					XmNtopAttachment, XmATTACH_FORM,
					XmNleftAttachment,  XmATTACH_FORM,
					NULL);
  XtAddCallback (PlayButton, XmNactivateCallback, MoviePlayButCB, 
		 (XtPointer)IDdata);
  
  /* Stop button */
  MovieDia.StopButton = XtVaCreateManagedWidget (" Stop ", 
						 xmPushButtonWidgetClass, 
						 form, 
						 XmNtopAttachment, XmATTACH_FORM,
						 XmNrightAttachment, XmATTACH_FORM,
						 NULL);
  XtAddCallback (MovieDia.StopButton, XmNactivateCallback, MovieStopButCB, 
		 (XtPointer)IDdata);
  
  /* start Plane */
  sprintf (valuestr, "%d", MovieDia.StartPlane);
  MovieDia.start = XtVaCreateManagedWidget ("MovieStartPlane", 
					    xmTextFieldWidgetClass, 
					    form, 
					    XmNwidth,           75,
					    XmNvalue,   valuestr,
					    XmNtopAttachment, XmATTACH_WIDGET,
					    XmNtopWidget,     MovieDia.StopButton,
					    XmNrightAttachment,  XmATTACH_FORM,
					    NULL);
  label1 = XtVaCreateManagedWidget ("MovieStartLabel", xmLabelWidgetClass, 
				    form, 
				    XmNwidth,           75,
				    XmNheight,          30,
				    XmNlabelString,   StartLab,
				    XmNtopAttachment, XmATTACH_WIDGET,
                                    XmNtopWidget,      PlayButton,
				    XmNleftAttachment,  XmATTACH_FORM,
				    NULL);
  
  /* end Plane */
  sprintf (valuestr, "%d", MovieDia.EndPlane);
  MovieDia.end = XtVaCreateManagedWidget ("MovieEndPlane", 
					  xmTextFieldWidgetClass, 
					  form, 
					  XmNwidth,           75,
					  XmNvalue,   valuestr,
					  XmNtopAttachment, XmATTACH_WIDGET,
					  XmNtopWidget,     MovieDia.start,
					  XmNrightAttachment,  XmATTACH_FORM,
					  NULL);
  label2 = XtVaCreateManagedWidget ("MovieEndLabel", xmLabelWidgetClass, 
				    form, 
				    XmNwidth,           75,
				    XmNheight,          30,
				    XmNlabelString,   EndLab,
				    XmNtopAttachment, XmATTACH_WIDGET,
                                    XmNtopWidget,      MovieDia.start,
				    XmNleftAttachment,  XmATTACH_FORM,
				    NULL);
  
  /* Dwell */
  sprintf (valuestr, "%f", MovieDia.Dwell);
  MovieDia.dwell = XtVaCreateManagedWidget ("MovieDwell", 
					    xmTextFieldWidgetClass, 
					    form, 
					    XmNwidth,           75,
					    XmNvalue,   valuestr,
					    XmNtopAttachment, XmATTACH_WIDGET,
					    XmNtopWidget,     MovieDia.end,
					    XmNrightAttachment,  XmATTACH_FORM,
					    NULL);
  
  label3 = XtVaCreateManagedWidget ("MovieDwellLabel", xmLabelWidgetClass, 
				    form, 
				    XmNwidth,           75,
				    XmNheight,          30,
				    XmNlabelString,   DwellLab,
				    XmNtopAttachment, XmATTACH_WIDGET,
                                    XmNtopWidget,      MovieDia.end,
				    XmNleftAttachment,  XmATTACH_FORM,
				    NULL);
  /* plane scroll */
  slider_size = hiPlane/10; if (slider_size<1) slider_size = 1;
  value = MovieDia.CurPlane+1;
  increment = 1;
  page_increment = hiPlane/10; if (page_increment<1) page_increment = 1;
  MovieDia.PlaneScroll = XtVaCreateManagedWidget ("MoviePlaneScroll", 
						  xmScrollBarWidgetClass, 
						  form,
						  XmNwidth,            150,
						  XmNheight,           20,
						  XmNmaximum,    hiPlane+slider_size,
						  XmNminimum,           1,
						  XmNvalue,       MovieDia.CurPlane+1,
						  XmNshowValue,       True,
						  XmNorientation, XmHORIZONTAL,
						  XmNprocessingDirection, XmMAX_ON_RIGHT,
						  XmNleftAttachment,  XmATTACH_FORM,
						  XmNtopAttachment, XmATTACH_WIDGET,
						  XmNtopWidget,     MovieDia.dwell,
						  XmNtopOffset,            10,
						  NULL);
  XmScrollBarSetValues (MovieDia.PlaneScroll, value, slider_size, 
			increment, page_increment, False);
  
  XtAddCallback(MovieDia.PlaneScroll, XmNvalueChangedCallback, 
                MovieScrollPlaneCB, (XtPointer)IDdata);
  
  /* scroll label */
  sprintf (valuestr, "    Plane no. %d", MovieDia.CurPlane+1);
  WierdString = XmStringCreateSimple (valuestr);
  MovieDia.PlaneScrollLabel = XtVaCreateManagedWidget ("MovieScrollLabel", 
						       xmLabelWidgetClass, 
						       form, 
						       XmNwidth,           150,
						       XmNtopAttachment, XmATTACH_WIDGET,
						       XmNtopWidget,     MovieDia.PlaneScroll,
						       XmNtopOffset,            10,
						       XmNleftAttachment,  XmATTACH_FORM,
						       NULL);
  XtVaSetValues(MovieDia.PlaneScrollLabel, 
		XmNlabelString,   WierdString,
		NULL);
  
  /* plane value label */
  sprintf (valuestr, "Plane value label");
  if (WierdString) XmStringFree(WierdString); WierdString = NULL;
  WierdString = XmStringCreateSimple (valuestr);
  MovieDia.PlaneLabel = XtVaCreateManagedWidget ("MoviePlaneLabel", 
						 xmLabelWidgetClass, 
						 form, 
						 XmNwidth,           150,
						 XmNlabelString,   WierdString,
						 XmNtopAttachment, XmATTACH_WIDGET,
						 XmNtopWidget,  MovieDia.PlaneScrollLabel,
						 XmNleftAttachment,  XmATTACH_FORM,
						 NULL);
  MovieLabelPlane(); /* fill in actual value */
  
  /* Quit button */
  QuitButton = XtVaCreateManagedWidget ("Quit", xmPushButtonWidgetClass, 
					form, 
					XmNbottomAttachment, XmATTACH_FORM,
					XmNrightAttachment, XmATTACH_FORM,
					NULL);
  XtAddCallback (QuitButton, XmNactivateCallback, MovieQuitButCB, 
		 (XtPointer)IDdata);
  
  if (WierdString) XmStringFree(WierdString); WierdString = NULL;
  if (StartLab) XmStringFree(StartLab); StartLab = NULL;
  if (EndLab) XmStringFree(EndLab); EndLab = NULL;
  if (DwellLab) XmStringFree(DwellLab); DwellLab = NULL;
  
  /* set it up */
  XtManageChild (MovieDia.dialog);
  
  /* put it some place reasonable */
  /*  where is parent? */
  XtVaGetValues (IDdata->shell,
		 XmNx, &MovieDia.xpos,
		 XmNy, &MovieDia.ypos,
		 NULL);
  MovieDia.ypos += 160;
  if (MovieDia.xpos<0) MovieDia.xpos = 0;
  XMoveWindow (XtDisplay(IDdata->shell), XtWindow(MovieDia.dialog), 
	       MovieDia.xpos, MovieDia.ypos);
  MovieDia.window = XtWindow(MovieDia.StopButton); /* save Stop Button window */
} /* end MovieBox */

