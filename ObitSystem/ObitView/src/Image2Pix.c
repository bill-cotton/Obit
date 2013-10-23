/* $Id$  */
/* routines to load a FITS file to a ZPixmap for ObitView*/
/*-----------------------------------------------------------------------
*  Copyright (C) 1996-2009
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
#include <unistd.h> 
#include <math.h>
#include <Xm/Xm.h> 
#include <Xm/DrawingA.h> 
#include <Xm/MainW.h> 
#include <Xm/DialogS.h> 
#include <Xm/Form.h>
#include <Xm/PushB.h>
#include <Xm/Label.h>
#include <X11/cursorfont.h> 
#include <pthread.h> 
#include "ObitSystem.h"
#include "ObitImage.h"
#include "ObitAIPS.h"
#include "ObitAIPSDir.h"
#include "ObitMem.h"
#include "imagedisp.h"
#include "messagebox.h"
#include "histo.h"

/**
 *  \file Image2Pix.c
 * utilities to read FITS image convert to ZPixMap.
 */

typedef struct { /* define ReadImage argument */
  ImageData *data;
  ObitErr *err;
}  ReadImageArg;

/*--------------- file global data ----------------*/
Widget   shell;
static Widget   dialog=NULL;
static Widget   line1, CancelButton;
static Boolean  stopped;
gboolean LoadDone;
gboolean ReadFail;
pthread_t readThread;
pthread_attr_t readThread_attr;
ObitImage *workImage;  /* in case load stopps */

/*---------------Private function prototypes----------------*/
void WorkingCursor(gboolean on, gboolean verbose);
Boolean CheckForCancel();
void* ReadImage (void *arg);

/*---------------Public functions ----------------*/

/**
 * Routine to convert an Obitimage to a pixel array with a color table.
 * Obit must have been initialized prior to this call
 * using ObitSystemStartup.  The first AIPS disk should be given
 * a dummy name which will be replaced here using data->AIPSDir.
 * \param image    Image to display
 * \param IDdata   Image Display
 * \param verbose  if true then a progress/cancel box appears
 * \return returns 0 if OK else failed.
 */
olong Image2Pix (ImageData *image, ImageDisplay *IDdata, gboolean verbose)
{
  olong     dims[2];
  unsigned long pasize, yaddr, addr;
  olong      i, j, jj, nx, ny, icol;
  ofloat    val, valmin, valmax;
  ofloat    *valP, blanked, irange, c1, c2, *eq_map=NULL;
  gboolean  newRange, bugOut;
  olong      oldType;
  ReadImageArg RIArg;
  
  /* Magic blanking value */
  blanked = ObitMagicF();
  /* set shell for working box */
  shell = IDdata->shell;

  /* Does the file need to be reloaded?
     if only the pix range or mapping type changed then the 
     image doesn't need to be reloaded */

  if ((image->reLoad) || (!image->valid)) {
    
    image->reLoad = FALSE; /* May not need next time */
    image->valid  = FALSE; /* while reloading */
    /* Start Read in another thread */
    /* Allow immediate cancellation of threads */
    pthread_setcanceltype( PTHREAD_CANCEL_ASYNCHRONOUS, &oldType);
    /* initialize attributes */
    if (pthread_attr_init(&readThread_attr)) {
      MessageShow ("Image2Pix: Failed to initialize attributes for thread");
      /*ObitThreadUnlock (image->thread);*/
      return -1;
    }
    RIArg.data = image;
    RIArg.err  = err;
    LoadDone = FALSE;
    ReadFail = FALSE;
    pthread_create (&readThread, &readThread_attr, ReadImage, 
		    (void*)&RIArg);
    /* restore cancellation type */
    pthread_setcanceltype(oldType, NULL);
    
    /*  start Cancel load dialog box */
    WorkingCursor(TRUE, verbose);
    
    /*  check for disgruntled user - will set LoadDone*/
    if (verbose) {
      /* This will terminate when Load done or user cancels */
      while (1) {
	usleep(250000);  /* sleep 250 msec */
	CheckForCancel();
	if (LoadDone) break;
      }
    }
    
    /* Wait for Load to finish */
    pthread_join (readThread, NULL);
    
    /* end progress message/cancel box if it's still there */
    WorkingCursor(FALSE, FALSE);
  } /* End load image */
  
  /* Did Load work? */
  if ((!image->valid) || err->error || ReadFail) {
    ObitErrLog(err); /* show any error messages on err */
    /* Clean up mess and bail out */
    workImage       = ObitImageUnref (workImage);
    image->myDesc   = ObitImageDescUnref(image->myDesc);
    image->myPixels = ObitImageDescUnref(image->myPixels);
    image->valid    = FALSE;  /* now no valid data */
    /*ObitThreadUnlock (image->thread);*/
    return -1;
  }
  
  /* Lock */
  ObitThreadLock (image->thread);
  
  /* set center scroll */
  IDdata->scrollx  = image->myDesc->inaxes[0] / 2; 
  IDdata->scrolly  = image->myDesc->inaxes[1] / 2; 
  IDdata->showInfo = FALSE; /* turn off display of pixel information */
  
  /*  Setup pixarray/gpharray if needed */
  dims[0] = image->myDesc->inaxes[0];
  dims[1] = image->myDesc->inaxes[1];
  if ((image->pixarray==NULL) || (image->gpharray==NULL) ||
      (image->nxArray != dims[0]) || (image->nyArray != dims[1]))   {
    image->nxArray = dims[0];
    image->nyArray = dims[1];

    /* delete old */
    if (image->pixarray) ObitMemFree ((gchar*) image->pixarray); 
    image->pixarray=NULL;
    if (image->gpharray) ObitMemFree ((gchar*) image->gpharray); 
    image->gpharray=NULL;
    /* create new pix map */
    pasize = dims[0] * dims[1]; /* keep as 8 bit */
    image->pixarray = (gchar*) ObitMemAlloc0 (pasize);
    /* graphics plane */
    image->gpharray = (gchar*) ObitMemAlloc0 (pasize);
  } /* end of (re)build pixarray/gpharray  */
  
  /* specify initial pixel range */
  get_extrema (image->myPixels, &valmax, &valmin);
  image->maxVal = valmax;
  image->minVal = valmin;
  newRange = TRUE;
  if (((image->PixRange[0]!=0.0) || 
       (image->PixRange[1]!=0.0))){ /*  User selected range */
    valmin = image->PixRange[0]; 
    valmax = image->PixRange[1];
    newRange = FALSE;
  }

  /* setup for pixel conversion */
  bugOut = 0;
  if (newRange) { /* find plausible range? */
    bugOut = get_range (image->myPixels, image->mapFunc, &valmax, &valmin);
  }
  if (image->mapFunc==2) { /* histogram equalization */
    bugOut = bugOut || equalize (image->myPixels, &valmax, &valmax, &eq_map);
  }
  if (bugOut) { /* failure in histogram equalization/pixel range */
    MessageShow ( "Error in histogram analysis");
    /* Clean up mess and bail out */
    if (eq_map) g_free(eq_map);  eq_map = NULL;
    image->myDesc   = ObitImageDescUnref(image->myDesc);
    image->myPixels = ObitImageDescUnref(image->myPixels);
    image->valid = FALSE;  /* now no valid data */
    ObitThreadUnlock (image->thread);  /* Unlock */
    return -1;
  }
  
  /* scaling for range of pixel values */
  irange = valmax - valmin;
  if (fabs(irange)<1.0e-25)
    irange = 1.0;
  else
    irange = 1.0 / irange;
  c1 = (MAXCOLOR - 1.0) * irange;
  c2 = valmin * c1 - 0.5;
  
  /* Convert to pixarray */
  valP = image->myPixels->array;  /* pointer in pixel array */
  nx = image->myDesc->inaxes[0]; 
  ny = image->myDesc->inaxes[1]; 
  for (j = 0; j<ny; j++) { /* loop over rows */
    /* pixarray is reversed top-to-bottom wrt image */
    jj = ny - j - 1;
    yaddr = jj * image->myDesc->inaxes[0];
    for (i=0; i<nx; i++) { /* loop over cols */
      if ((i<dims[0]) && (j<dims[1]))
	val = *(valP++);
      else 
	val = blanked;  /* blank fill edge? */
      if (val!=blanked)  { /* Set color index. */
	if (image->mapFunc==1)  { /* square root */
	  if (val<valmin) val=valmin;
	  icol =  (((MAXCOLOR)-1.0) *
		   sqrt(((val-valmin) * irange))+0.5);
	} else if (image->mapFunc==2) { /* histogram equalization */
	  icol = map_pixel (eq_map, val);
	} else  /* Linear */
	  icol = c1 * val - c2;
	if (icol<1) icol = 1;  /* minimum color = 1 */
	if (icol>MAXCOLOR-1) icol = MAXCOLOR-1;
      } else
	icol = 0;  /* blanked pixel */
      /* save pixel value */
      addr = yaddr + i;
      *(image->pixarray+addr) = IDdata->colut[icol];
    } /* end loop over columns */
  } /* end loop over rows */

  /* Unlock */
  ObitThreadUnlock (image->thread);
  if (eq_map) g_free(eq_map);  /* Cleanup */
  return 0;
}  /* end of Image2Pix */

/*---------------Private functions ----------------*/
/**
 * Callback for user hit "cancel" button on load box
 * \param w           scroll widget activated
 * \param clientData  client data
 * \param callData    call data
 */
void StopFITSLoadButCB (Widget dialog, XtPointer clientData, 
			XtPointer callData)
{
  if (!stopped) pthread_cancel(readThread);
  LoadDone = TRUE;
  stopped  = TRUE;
} /* end StopFITSLoad */

/**
 * If verbose create working dialog else change cursor to hourglass 
 * \param on       Want diaglog displayed or not
 * \param verbose  if true then a cancel box appears
 * \param filename Name of fits file being loaded
 */
void WorkingCursor(gboolean on, gboolean verbose)
{
  static olong locked, boxup=FALSE;
  static Cursor cursor=(Cursor)NULL;
  extern Widget shell;
  XSetWindowAttributes attrs;
  Display *dpy = XtDisplay(shell);
  XmString str = NULL;
  olong  newdialog;
  Widget form, line2;
  gchar cstring[121];
  
  /* keep track of how many times this is "locked" */
  on? locked++:locked--;
  if (((locked>1) || (locked ==1)) && on ==0)
    return; /* already locked */
  if (on) stopped = False; /* initialize */
  if (!cursor) /* make sure hourglass cursor initialized */
    cursor = XCreateFontCursor(dpy, XC_watch);
  /* if on and verbose bring up working dialog */
  if (on && verbose)  {
    /* make working dialog box */
    newdialog = !dialog;
    if (newdialog)
      {
	dialog = XtVaCreatePopupShell ("Load Image", 
				       xmDialogShellWidgetClass, 
				       shell, 
				       XmNwidth,     250,
				       XmNheight,    100,
				       XmNdeleteResponse, XmDESTROY,
				       NULL);
	/* make Form widget to stick things on */
	form = XtVaCreateManagedWidget ("WorkingForm", xmFormWidgetClass,
					dialog,
					XmNwidth,     250,
					XmNheight,    100,
					XmNx,           0,
					XmNy,           0,
					NULL);
	
	/* info label widgets */
	g_snprintf (cstring, 120, "Loading Image");
	str = XmStringCreateLocalized (cstring);
	line1 = XtVaCreateManagedWidget ("Line1", xmLabelWidgetClass, 
					 form, 
					 XmNlabelString,     str,
					 XmNtopAttachment,   XmATTACH_FORM,
					 XmNrightAttachment, XmATTACH_FORM,
					 XmNleftAttachment,  XmATTACH_FORM,
					 NULL);
	if (str) XmStringFree(str); str = NULL;
	
	g_snprintf (cstring, 120, "Canceling FITS load behaves poorly");
	str = XmStringCreateLocalized (cstring);
	line2 = XtVaCreateManagedWidget ("Line2", xmLabelWidgetClass, 
					 form, 
					 XmNlabelString,   str,
					 XmNtopAttachment, XmATTACH_WIDGET,
					 XmNtopWidget,     line1,
					 XmNrightAttachment, XmATTACH_FORM,
					 XmNleftAttachment,  XmATTACH_FORM,
					 NULL);
	if (str) XmStringFree(str); str = NULL;
	
	/* Cancel button */
	CancelButton = XtVaCreateManagedWidget ("Cancel", 
						xmPushButtonWidgetClass, 
						form, 
						XmNbottomAttachment,XmATTACH_FORM,
						XmNrightAttachment, XmATTACH_FORM,
						XmNleftAttachment,  XmATTACH_FORM,
						NULL);
	XtAddCallback (CancelButton, XmNactivateCallback, 
		       StopFITSLoadButCB,  NULL);
	/* Shazam appear: */
	XtManageChild (dialog);
	boxup = TRUE;
       /* end create box */
      } else { /* dialog exists - check if done */
	/*  check for disgruntled user - will set LoadDone*/
	CheckForCancel();
	/* check if done */
	if (LoadDone) {  /* Done */
	  /* reset cursor */
	  if (!boxup) {
	    attrs.cursor = on ? cursor : None;
	    XChangeWindowAttributes(dpy, XtWindow(shell), CWCursor, &attrs);
	  }
	  XFlush(dpy);
	  /* kill box if user aborted load - 
	     I haven't a clue as to why this is needed */
	  if (dialog) { 
	    XtDestroyWidget(dialog);
	    dialog = NULL;}
	  
	  boxup = FALSE;
	} /* end of shutdown/reset */
      } /* end check if done */
  } /* end of start up working dialog */

  /* only change cursor */
  else if (on && (!verbose))
    {
      boxup = FALSE;
      attrs.cursor = on ? cursor : None;
      XChangeWindowAttributes(dpy, XtWindow(shell), CWCursor, &attrs);
      XFlush(dpy);
      /* end of change cursor */
      /* done with working dialog */
    } else if (!on) {
      /* reset cursor */
      if (!boxup) {
	attrs.cursor = on ? cursor : None;
	XChangeWindowAttributes(dpy, XtWindow(shell), CWCursor, &attrs);
      }
      XFlush(dpy);
      /* kill box if user aborted load - I haven't a clue as to why this is needed */
      if (dialog) { 
	XtDestroyWidget(dialog);
	dialog = NULL;}
      
      boxup = FALSE;
    } /* end of shutdown/reset */
} /* end WorkingCursor */

/**
 * See if user hit cancel button 
 * \return True if hit
 */
Boolean CheckForCancel()
     /* See if user hit cancel button */
{
  extern Widget shell;
  Display *dpy = XtDisplay(shell);
  Window win = XtWindow(CancelButton);
  XEvent event;
  
  /* Make sure all our requests get to the server */
  XFlush(dpy);
  
  /* take care of pending expose events */
  XmUpdateDisplay(shell);
  
  /* check the event queue for the cancel button */
  while(XCheckMaskEvent(dpy,
			ButtonPressMask | ButtonReleaseMask | 
			ButtonMotionMask | PointerMotionMask |
			KeyPressMask | KeyReleaseMask,
			&event)) { /* caught possibly intresting one */
    if (event.xany.window==win)
      {XtDispatchEvent(&event); /* it's in the cancel button - do it */
      /* Cancel loading thread */
      if (!stopped) {
	pthread_cancel(readThread);
	MessageShow ( "Canceling Image Load");
      }
      LoadDone = TRUE;
      stopped  = TRUE;}
  }
  return stopped;
} /* End CheckForCancel */

/**
 * Routine to read an Obit Image
 * Obit must have been initialized prior to this call
 * using ObitSystemStartup.  The first AIPS disk should be given
 * a dummy name which will be replaced here using data->AIPSDir.
 * \param arg a ReadImageArg pointer with arguments
 * \li data  ImageData structure 
 * \li err   Obit error/message structure
 */
void* ReadImage (void *arg)
{
  ObitImage *image=NULL;
  olong blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  olong disk, cno;
  ReadImageArg *RIArg = (ReadImageArg*)arg;
  ImageData* data = RIArg->data;
  ObitErr *err    = RIArg->err;
  ObitImageDesc *desc=NULL;
  gchar imageName[100];
  gchar AName[13], AClass[7];
  gchar *routine  = "ReadImage";

  /* Remove old structures */
  data->myDesc   = ObitImageDescUnref(data->myDesc);
  data->myPixels = ObitImageDescUnref(data->myPixels);
  data->valid = FALSE;  /* now no valid data */

  if (err->error) {
    LoadDone = TRUE;
    return NULL;  /* previous error? */
  }

  /* Create basic image */
  image = newObitImage("ObitView Image");
  workImage = image;  /* file global */

  /* Set plane */
  blc[2] = trc[2] = data->PlaneNo+1;
  blc[3] = trc[3] = data->hiDim[0]+1;
  blc[4] = trc[4] = data->hiDim[1]+1;
  blc[5] = trc[5] = data->hiDim[2]+1;

  /* Attach to external form by data type */
  if (data->DataType==OBIT_IO_FITS) {   /* FITS */
    /* Assume full path */
    disk = 0;
    ObitImageSetFITS(image,OBIT_IO_byPlane,disk,data->FileName->sp,blc,trc,err);
  } else if (data->DataType==OBIT_IO_AIPS) {  /* AIPS */
    disk = 1;
    ObitSystemSetAIPSuser (data->AIPSuser);
    ObitAIPSSetDirname (disk, data->AIPSDir->sp, err);
    /* Lookup */
    disk = 1;
    strncpy (AName, data->AIPSName->sp, 12);
    AName[12] = 0;
    strncpy (AClass, data->AIPSClass->sp, 6);
    AClass[6] = 0;
    cno = ObitAIPSDirFindCNO(disk, data->AIPSuser, AName, AClass,
			     "MA", data->AIPSseq, err);
   if (cno<=0) {
      LoadDone = TRUE;
      ReadFail = TRUE;
      Obit_log_error(err, OBIT_Error, "Failure looking up input file");
      image = ObitImageUnref(image);  /* Cleanup */
      return NULL;
    }
    
    /* Attach */
    ObitImageSetAIPS(image,OBIT_IO_byPlane,disk,cno,data->AIPSuser,blc,trc,err);

    /* Replace file name */
    g_snprintf (imageName, 100, "%s.%s.%d user %d",  
		AName, AClass, data->AIPSseq, data->AIPSuser);
    FStrngFill (data->FileName, imageName);

  } /* end attaching to external */
  /* Error? */
  if (err->error) {
    LoadDone = TRUE;
    ReadFail = TRUE;
    image = ObitImageUnref(image);  /* Cleanup */
    Obit_traceback_val (err, routine, "LoadImage", NULL);
  }
  
  /* Open */
  ObitImageOpen (image, OBIT_IO_ReadOnly, err);
  /* Error? */
  if (err->error) {
    ReadFail = TRUE;
    LoadDone = TRUE;
    image = ObitImageUnref(image);  /* Cleanup */
    Obit_traceback_val (err, routine, "LoadImage", NULL);
  }

  /* Read */
  ObitImageRead (image, NULL, err);
  /* Error? */
  if (err->error) {
    ReadFail = TRUE;
    LoadDone = TRUE;
    image = ObitImageUnref(image);  /* Cleanup */
    Obit_traceback_val (err, routine, "LoadImage", NULL);
  }

  /* Save descriptor and pixel array */
  data->myDesc   = ObitImageDescRef(image->myDesc);
  data->myPixels = ObitImageDescRef(image->image);
  data->valid = FALSE;   /* Not yet */

  /* Get number of planes from IO descriptor */
  desc = (ObitImageDesc*)image->myIO->myDesc;
  data->myDesc->inaxes[2] = desc->inaxes[2];
  data->myDesc->inaxes[3] = desc->inaxes[3];
  data->myDesc->inaxes[4] = desc->inaxes[4];

  /* Close */
  ObitImageClose (image, err);
  /* Error? */
  if (err->error) {
    ReadFail = TRUE;
    LoadDone = TRUE;
    image = ObitImageUnref(image);  /* Cleanup */
    Obit_traceback_val (err, routine, "LoadImage", NULL);
  }

  /* Cleanup */
  image = ObitImageUnref(image);
  workImage = image;  /* file global */

  data->valid = TRUE;  /* Now have valid data */
  LoadDone = TRUE;
  ReadFail = FALSE;
  return NULL;
} /* end ReadImage */
