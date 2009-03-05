/* $Id:  $ */
/*    ObitMess: Task message server for Obit */
/* This program requires the Motif library */
/* Cloned from ObitView */
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
*
*  You should have received a copy of the GNU General Public
*  License along with this program; if not, write to the Free
*  Software Foundation, Inc., 675 Massachusetts Ave, Cambridge,
*  MA 02139, USA.
*
*  Correspondence concerning this software should be addressed as follows:
*         Internet email: bcotton@nrao.edu.
*         Postal address: William Cotton
*                         National Radio Astronomy Observatory
*                         520 Edgemont Road
*                         Charlottesville, VA 22903-2475 USA
*-----------------------------------------------------------------------*/
#define OBITMESSMAIN
#include <Xm/Xm.h> 
#include <Xm/MainW.h> 
#include <stdlib.h>
#include <stdio.h>
#include <Xm/Separator.h>
#include <Xm/ScrolledW.h>
#include <Xm/Label.h>
#include <Xm/Form.h>
#include <Xm/PushB.h>
#include <Xm/PanedW.h>
#include "obitmess.h"
#include "cursor.h"
#include "TMessMenu.h"
#include "taskmessage.h"
#include "messagebox.h"
#include "XMLRPCTaskMessServer.h"
#include "Obit.h"
#include "ObitSystem.h"

/**
 *  \file ObitMess.c
 * ObitMess main program
 */

/*---------------Private function prototypes----------------*/
static void QuitCB(Widget widget, 
		   XtPointer client_data, 
		   XtPointer call_data);

/*--------------- file global data ----------------*/
/* stuff for error messages
olong hwndErr = 1;
gchar szErrMess[120]; */

/**
 * Main ObitMess program
 * Checks command line for option of form "-port xxxx'
 * If found, port xxx is used unless xxx='none'.
 * If xxx='none', no xmlrpc interface is started
 * "-port"  and "xxxx" removed from arguments
 * Default port = 8777
 * Initializes and starts the event loop
 * \param argc  number of command line arguments
 * \param argv  array of command line arguments
 * \return 0 is successful
 */
int main ( int argc, char **argv ) 
{
  Widget       pane, label1, QuitBut;
  XtAppContext app_context;
  Arg          args[10];
  Cardinal     n, xsize, ysize;
  pthread_t    thread[2];
  pthread_attr_t thread_attr[2];
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
  if (port<=0) port = 8777;

  /* initialize attributes */
  for (i=0;i<2;i++) {
    if (pthread_attr_init(&thread_attr[i])) {
      fprintf (stderr, "Failed to initialize attributes for thread %d\n",i);}
  }
  
    /* Start xmlrpc thread if needed */
  if (!noPort) {
    pthread_create (&thread[0], &thread_attr[0], start_myTWServer, (void*)&port);
  }

  /* Define main window */
  xsize = 200;
  ysize = 70;
  n = 0;
  XtSetArg(args[n], XmNwidth,  xsize); n++;
  XtSetArg(args[n], XmNheight, ysize); n++;
  XtSetArg(args[n], XmNx,       50); n++;
  XtSetArg(args[n], XmNy,       50); n++;
  parent = XtAppInitialize(&app_context,
			       "ObitMess",
			       (XrmOptionDescList)NULL,
			       0, &argc, argv,
			       (String*)NULL, args, n);
  /* make Paned-window widget to stick things on */
  n = 0;
  pane = XmCreatePanedWindow (parent,"Pane", args,n);
  
  n = 0;
  label1 = XmCreateLabel(pane, "Obit task message server", args, n);
  /* Manage */
  XtManageChild(label1);

  /* Quit button */
  n = 0;
  QuitBut = XmCreatePushButton (pane, "Quit", args, n);
  /* Add callback */
  XtAddCallback(QuitBut, XmNactivateCallback, 
		(XtCallbackProc)QuitCB, (XtPointer)NULL);
  /* Manage */
  XtManageChild(QuitBut);


  /* Now manage pane */
  XtManageChild(pane);
  /*  Presto - appear on the screen */
  XtRealizeWidget (parent);

  /* Initialize taks windowing */
  TaskMessageInit();
  
  /* Setup periodic checks for XMLRPC requests*/
  XtAppAddWorkProc (app_context, (XtWorkProc)XMLRPCTWWatcher, parent); 
  myApp = app_context;  /* Global version */
  Display_shell = parent;

  /* Initialize Obit error/message */
  err = newObitErr();

  /* Other init */
  xmlStatus    = 0;
  taskID = -1;


  /* main event loop */
  XtAppMainLoop (app_context);
  
  /* when X-windows terminates quit XMLRPC thread */
  pthread_join(thread[0], NULL);
  pthread_cancel (thread[0]); /* cancel RPC thread */
  pthread_detach (thread[0]);
  
  return 0;
} /* end of main */


/**
 * Callback for Quit button - bail out
 * \param w           widget activated
 * \param clientData  client data
 * \param callData    call data
 */
static void QuitCB (Widget w, XtPointer clientData, XtPointer callData)
{
  /* terminate program */
  exit (0);
} /* end QuitCB */


