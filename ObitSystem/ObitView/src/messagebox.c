/* $Id: messagebox.c,v 1.5 2007/01/11 19:11:55 bcotton Exp $ */
/* MessageBox routines for ObitView */
/* uses a ScrollText to display messages from ObitView      */
/*-----------------------------------------------------------------------
*  Copyright (C) 1996-2008
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
#include "scrolltext.h"
#include "messagebox.h"
#include "obitview.h"
#include <stdio.h> 

/**
 *  \file messagebox.c
 * displays "message" dialog.
 */

/*--------------- file global data ----------------*/
/**  MessageBox ScrollText, if non NULL then message box active */
ScrollTextPtr MessScroll = NULL; 

/*---------------Private function prototypes----------------*/
/* shutdown messagebox, arg not used */
void MessageDismiss(XtPointer arg);

/*-----------------public functions--------------*/

/**
 * Display message in message box, creating if necessary
 * \param message  message to display
 */
void MessageShow (char *message)
{
  int next, length, new = 0;
  
  /* new ScrollText? */
  if (!MessScroll)
    {
      /* make ScrollText */
      new = 1;
      MessScroll = ScrollTextMake (Display_shell, "ObitView Messages");
      if (!MessScroll) { /* error, print to stderr */
	fprintf (stderr, message); fprintf (stderr,"\n");
	return;
      } /* end create error */
      /* add dismiss callback */
      MessScroll->DismissProc = (TextFileProc)MessageDismiss;
    } /* end create ScrollText */
  
  /* add text */
  next = MessScroll->num_lines;
  if (next>=MAX_LINE) next = MAX_LINE - 1; /* add to end no matter */
  length = strlen(message);
  MessScroll->lines[next] = (char*)g_malloc(length+1);
  strcpy (MessScroll->lines[next], message);
  next++;
  MessScroll->num_lines = next; /* save number in ScrollText */
  /* make some noise */
  XBell(XtDisplay(MessScroll->Parent), 50); 
  /*  setup*/
  ScrollTextInit (MessScroll); 
  /* Go to bottom */
  ScrollTextBottom (MessScroll); 
  
  if (!new) { /* pop to front */
    if (XtIsRealized (MessScroll->ScrollBox))
      XMapRaised (XtDisplay(MessScroll->Parent), 
		  XtWindow(MessScroll->ScrollBox));
    /* redraw */
    STextExposeCB (MessScroll->ScrollBox, (XtPointer)MessScroll, NULL);
  }
  /* DEBUG
     fprintf(stderr,"MessageShow: %s\n",message); */
  
} /* end MessageShow */

/**
 * g_log handler for Obit error messages
 * \param message  message to display
 */
void MessageObitHandler (const gchar *log_domain,
			 GLogLevelFlags log_level,
			 const gchar *message,
			 gpointer user_data)
{
  MessageShow ((gchar*)message);
} /* end MessageObitHandler */


/**
 * Shutdown messagebox
 *  as this is called when the ScrollText self destructs it does not delete
 *  the ScrollText
 * \param arg   not used
 */
void MessageDismiss(XtPointer arg)
{
  MessScroll = NULL;
} /* end MessageDismiss */

/**
 * redraw message box if it exists
 */
void MessageRefresh (void)
{
  /* does it exist? */
  if (!MessScroll) return;
  
  /*  setup */
  ScrollTextInit (MessScroll);
  if (XtIsRealized (MessScroll->ScrollBox))
    XMapRaised (XtDisplay(MessScroll->Parent), 
		XtWindow(MessScroll->ScrollBox));
  /* redraw */
  STextExposeCB (MessScroll->ScrollBox, (XtPointer)MessScroll, NULL);
} /* end MessageRefresh */
