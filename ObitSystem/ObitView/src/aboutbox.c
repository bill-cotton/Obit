/* $Id$  */
/* about dialog box  for ObitView */
/*-----------------------------------------------------------------------
*  Copyright (C) 1998-2016
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
#include <stdio.h>
#include "obitview.h"
#include "scrolltext.h"
#include <ObitVersion.h>

/**
 *  \file aboutbox.c
 * displays "about" dialog.
 */

/*---------------Private function prototypes----------------*/
void HelpAbout (ScrollTextPtr STextPtr);

/**
 * Callback 
 * \param w           scroll widget activated
 * \param clientData  client data
 * \param callData    call data
 */
void HelpAboutCB (Widget w, XtPointer clientData, XtPointer callData) 
{
  ScrollTextPtr STextPtr;
  /* make ScrollText */
  STextPtr = ScrollTextMake (Display_shell, "About ObitView");
  
  /* copy text */
  if (STextPtr) HelpAbout(STextPtr);
  /* final setup */
  ScrollTextInit (STextPtr);
} /* end HelpAboutCB */

/**
 * Supply About text
 * \param STextPtr Text box to write to
 */
void HelpAbout (ScrollTextPtr STextPtr)  
{
  int loop, next, length;
  char *ver = ObitVersion();
  char label[80];
  g_snprintf (label,79,"ObitView %s Viewer for images in FITS or AIPS format\n", ver);
  g_free(ver);
  char *line[] = {
    "ObitView 1.3 Viewer for images in FITS or AIPS format ",
    "Copyright NRAO/AUI 2005-2016 ",
    " ",
    "   This software is distributed free of charge by NRAO. ",
    "The (USA) National Radio Astronomy Observatory (http://www.nrao.edu/) ",
    "is a facility of the (USA) National Science Foundation operated under ",
    "cooperative agreement by Associated Universities, Inc.  ",
    "(http://www.aui.edu).  ", 
    "Only very limited user support of this software is available.",
    "Suggestions and comments should be sent to Bill Cotton at NRAO ",
    "(bcotton@nrao.edu). NRAO's headquarters address is:  ",
    "National Radio Astronomy Observatory ",
    "520 Edgemont Road, ",
    "Charlottesville, VA 22903, USA",
    " ",
    "This Software is distributed in the hope that it will be useful, but ",
    "WITHOUT ANY WARRANTY; without even the implied warranty of ",
    "MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. ",
    " ",
    "*** FINISHED ***" }; /* end of text */
  line[0] = label;  /* version info */
  loop = 0;
  next = STextPtr->num_lines;
  while (1) { /* loop till done */
    if (!strcmp (line[loop], "*** FINISHED ***")) break; /* Finished */
    if (next>=MAX_LINE) break;
    /* copy */
    length = strlen(line[loop]);
    STextPtr->lines[next] = (char*)g_malloc(length+1);
    strcpy (STextPtr->lines[next], line[loop]);
    loop++; next++;
  } /* end loop loading info */
  STextPtr->num_lines = next; /* save number in ScrollText */
} /* end HelpAbout */


