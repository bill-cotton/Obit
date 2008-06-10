/* $Id:     */
/*    ObitView: image viewer for Obit */
/* This program requires the Motif library */
/* Cloned from ObitView */
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
*
*  You should have received a copy of the GNU General Public
*  License along with this program; if not, write to the Free
*  Software Foundation, Inc., 675 Massachusetts Ave, Cambridge,
*  MA 02139, USA.
*
*  Correspondence concerning ObitView should be addressed as follows:
*         Internet email: bcotton@nrao.edu.
*         Postal address: William Cotton
*                         National Radio Astronomy Observatory
*                         520 Edgemont Road
*                         Charlottesville, VA 22903-2475 USA
*-----------------------------------------------------------------------*/
#define OBITVIEWMAIN
#include <Xm/Xm.h> 
#include <Xm/DrawingA.h> 
#include <Xm/MainW.h> 
#include <stdlib.h>
#include <stdio.h>
#include <Xm/Separator.h>
#include <Xm/ScrolledW.h>
#include <Xm/Scale.h>
#include <Xm/Form.h>
#include <Xm/PushB.h>
#include <Xm/CascadeB.h>
#include <Xm/RowColumn.h>
#include <Xm/Label.h>
#include <Xm/ToggleB.h>
#include "Obit.h"
#include "obitview.h"
#include "imagedisp.h"
#include "Image2Pix.h"
#include "color.h"
#include "optionbox.h"
#include "moviebox.h"
#include "blinkbox.h"
#include "infobox.h"
#include "control.h"
#include "markpos.h"
#include "cursor.h"
#include "textfile.h"
#include "scrolltext.h"
#include "menu.h"
#include "messagebox.h"
#include "toolbox.h"
#include "drawbox.h"
#include "XMLRPCserver.h"
#include "Obit.h"
#include "ObitSystem.h"

/**
 *  \file ObitView.c
 * ObitView main program
 */

/*---------------Private function prototypes----------------*/
void InitImage (ImageDisplay *IDdata, int narg, char *filename);


/*--------------- file global data ----------------*/
/* stuff for error messages */
olong hwndErr = 1;
gchar szErrMess[120];

/**
 * Main ObitView program
 * Checks command line for option of form "-port xxxx'
 * If found, port xxx is used unless xxx='none'.
 * If xxx='none', no xmlrpc interface is started
 * "-port"  and "xxxx" removed from arguments
 * Default port = 8765
 * Initializes and starts the event loop
 * \param argc  number of command line arguments
 * \param argv  array of command line arguments
 * \return 0 is successful
 */
int main ( int argc, char **argv ) 
{
  Widget       mainWindow, menu, control, form;
  /*Widget       toolbox;*/
  ImageData    image;
  ImageDisplay *IDdata;
  Arg          args[10];
  olong         n, user;
  guint dbug;
  pthread_t    thread[2];
  pthread_attr_t thread_attr[2];
  ObitSystem *mySystem=NULL;
  gchar *AIPSdir[] = {"Dummy"};
  gchar *FITSdir[] = {"."};
  olong           i, j, port;
  gchar          *arg;
  gboolean       noPort=FALSE;

  /* Check if port number given */
  port = -1; 
  for (i=1;i<(argc)-1; i++) {
    arg = argv[i];
    if (strcmp (arg, "-port")==0) {
      /* no port request? */
      if (strcmp (argv[i+1], "none")==0) noPort = TRUE;
      /* port request */
      else port = (olong)strtol(argv[i+1], NULL, 0);
      
      /* remove from arguments */
      for (j=i+2; j<(argc); j++) argv[j-2] = argv[j];
      argc -= 2;
      
    }
    if (port>0) break;  /* found it */
  } /* end loop over arguments */

  /* default port */
  if (port<=0) port = 8765;

  /* initialize attributes */
  for (i=0;i<2;i++) {
    if (pthread_attr_init(&thread_attr[i])) {
      fprintf (stderr, "Failed to initialize attributes for thread %d\n",i);}
  }
  
    /* Start xmlrpc thread if needed */
  if (!noPort) {
    pthread_create (&thread[0], &thread_attr[0], start_myServer, (void*)&port);
  }

  /*  initialize toolkit */
  Display_shell = XtAppInitialize (&myApp, 
				   "ObitView", 
				   NULL, 0, &argc, argv,  NULL, NULL, 0); 
  /* Set initial size */
  n = 0;
  /*  XtSetArg(args[n], XtNwidth, 640);  n++;
      XtSetArg(args[n], XtNheight, 480);  n++; */
  
  /* create main window */
  mainWindow = XtCreateManagedWidget ("mainWindow", xmMainWindowWidgetClass,
				      Display_shell, args, n);
  /* make a form to hang everything on*/
  form = XtVaCreateManagedWidget ("form", xmFormWidgetClass,
				  mainWindow,
				  XmNwidth,     640,
				  XmNheight,    480,
				  XmNx,           0,
				  XmNy,           0,
				  NULL);
  
  
  /* make image display window */
  IDdata = MakeDisplay (form, Display_shell);

  /* Initialize drawing */
  DrawBoxInit (IDdata);

  /* Initialize Obit error/message */
  err = newObitErr();
  /* Set handler so that Obit messages appear in Message Box */
  dbug = g_log_set_handler (NULL,  G_LOG_LEVEL_WARNING | G_LOG_FLAG_FATAL |
			    G_LOG_LEVEL_CRITICAL | G_LOG_LEVEL_MESSAGE |
			    G_LOG_LEVEL_INFO | G_LOG_LEVEL_DEBUG |
			    G_LOG_FLAG_RECURSION, 
			    (GLogFunc)MessageObitHandler, NULL);
   /* Initialize Obit */
  user = 1;
  mySystem = ObitSystemStartup ("", 1, user, 1, AIPSdir, 1, FITSdir, 
				(oint)TRUE, (oint)FALSE, err);
  ObitErrLog(err); /* show any error messages on err */

  /*  make display control panel */
  control = MakeControl (form, (XtPointer)IDdata);
  
  /*  make t panel */
  /* not ready for prime time */
  /*  toolbox = MakeToolbox (form, control, (XtPointer)IDdata); */
  
  /* attach image display to control on form */
  XtVaSetValues(IDdata->display, 
		XmNrightAttachment,  XmATTACH_FORM,
		XmNleftAttachment,  XmATTACH_WIDGET,
		XmNleftWidget,  control,
		XmNtopAttachment,  XmATTACH_FORM,
		XmNbottomAttachment,  XmATTACH_FORM,
		NULL);
  
  /*  make main menu */
  menu = MakeMainMenu (mainWindow, (XtPointer)&image, 
		       (XtPointer)IDdata);
  
  /*  save some widget names */
  XtVaSetValues ( mainWindow, XmNmenuBar, menu, NULL );
  
  /*  Presto - appear on the screen */
  XtRealizeWidget (Display_shell);
  
  /*  create / init color map, cursor */
  SetupColorMap (Display_shell, IDdata);
  
  /* set cursor if possible */
  IDdata->cursor = MakeImageCursor (XtDisplay (IDdata->display), 
				    XtWindow(IDdata->canvas));
  if (IDdata->cursor)
    XDefineCursor (XtDisplay (IDdata->display), XtWindow(IDdata->canvas),
		   IDdata->cursor);

  /* initialize image, display file given as argument */
  InitImage (IDdata, argc, argv[1]);
  
  /*   set display */
  SetDisplay(IDdata);
  
  /* save application context */
  IDdata->app = myApp;
  
  /* DEBUG test Obit error
  Obit_log_error(err, OBIT_InfoErr,"Test Obit Message, dbug =  %d", dbug);
  ObitErrLog(err);
  g_warning("Test g warning");
  g_log (NULL,G_LOG_LEVEL_WARNING,
	 "%s: %s", "SomeLevel", "SomeMessage"); */

  /* redraw message box if it contains any start up messages */
  MessageRefresh();
  
  /* Setup periodic checks for XMLRPC requests */
  XtAppAddWorkProc (myApp, (XtWorkProc)XMLRPCWatcher, IDdata); 

  /* main event loop */
  XtAppMainLoop (myApp);
  
  /* when X-windows terminates quit XMLRPC thread */
  pthread_join(thread[0], NULL);
  pthread_cancel (thread[0]); /* cancel RPC thread */
  pthread_detach(thread[0]);
  /* Shutdown Obit */
  mySystem = ObitSystemShutdown (mySystem);
  
  return 0;
} /* end of main */

/**
 * Initialize image display
 * \param IDdata   Image display to initialize
 * \param narg     number of command line arguments, if >=2 a filename is assumed
 * \param filename Name of FITS file to open.
 */
void InitImage (ImageDisplay *IDdata, olong narg, gchar *filename)
     /* initialize image data structures */
{
  olong i,j,k;
  gchar direct[121];
  
  CurImag     = 0; /* set current image pointer */
  FITS_dir    = NULL;
  mark_dir    = NULL;
  log_dir     = NULL;
  requestList = NULL;
  doLog = 0;   /* position logging initially turned off */
  usr_equinox = -1.0;  /* use equinox of image */
  ERtime_out  = -1.0;  /* Edit request Timeout, init forever */


  /* Initialize structures */
  for (j=0;j<2;j++) {
    image[j].valid     = FALSE;
    image[j].reLoad    = TRUE;
    image[j].pixarray  = NULL;
    image[j].gpharray  = NULL;
    image[j].nxArray   = 0;
    image[j].nyArray   = 0;
    image[j].myDesc    = NULL;
    image[j].myPixels  = NULL;
    image[j].DataType  = OBIT_IO_FITS;
    image[j].FileName  = MakeFStrng("NONE");
    image[j].AIPSName  = MakeFStrng("NONE");
    image[j].AIPSClass = MakeFStrng("NONE");
    image[j].AIPSDir   = MakeFStrng("NONE");
    image[j].AIPSseq   = 0;
    image[j].AIPSuser  = 0;
    image[j].Field     = 0;
    image[j].NField    = 0;
    image[j].maxVal    = 0.0;
    image[j].minVal    = 0.0;
    image[j].iXPixel   = 0;
    image[j].iYPixel   = 0;
    image[j].fXpixel   = 0.0;
    image[j].fYpixel   = 0.0;
    image[j].PixRange[0] = 0.0;
    image[j].PixRange[1] = 0.0;
    image[j].fBpixel     = 0;
    image[j].mapFunc     = 0;
    image[j].PlaneNo     = 0;
  } /* end of loop over ImageData structures */

  /* was a FITS file name passed as an argument? */
  if (narg<2) return;
  
  /* read FITS file to pixmap */
  FStrngFill (image[CurImag].FileName, filename);
  image[CurImag].DataType = OBIT_IO_FITS;    /* FITS */
  image[CurImag].reLoad   = TRUE;            /* Force load */
  /* read FITS file to PixMap */
  if (Image2Pix (&image[CurImag], IDdata, TRUE))
    { /* error */
      g_snprintf (szErrMess, 120, "Error reading FITS file %s", filename);
      MessageShow (szErrMess);
    }
  /* get directory name */
  j = strlen(filename);
  k = 0;
  for (i=0;i<j;i++) if (filename[i]=='/') k = i+1;
  for (i=0;i<k;i++) direct[i]=filename[i]; direct[i]=0;
  FITS_dir = MakeFStrng(direct);
  PaintImage(IDdata); /* draw image */
} /* end of InitImage */

