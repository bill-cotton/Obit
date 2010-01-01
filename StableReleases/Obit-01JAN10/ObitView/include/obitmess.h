/* $Id:  $ */
/* ObitMess header file */
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
*  Correspondence concerning this software should be addressed as follows:
*         Internet email: bcotton@nrao.edu.
*         Postal address: William Cotton
*                         National Radio Astronomy Observatory
*                         520 Edgemont Road
*                         Charlottesville, VA 22903-2475 USA
*-----------------------------------------------------------------------*/
#include <Xm/Xm.h> 
#include <Xm/DrawingA.h> 
#include <Xm/MainW.h> 
#include <stdlib.h>
#include <stdio.h>
#include "glib.h"
#include "FStrng.h"
#include "ObitErr.h"
#include "ObitTypes.h"
#ifndef OBITMESS
#define OBITMESS

/* List of task message windows 
return ID when creating and use it to pick window */

/* global data structures */
#ifdef OBITMESSMAIN
XtAppContext myApp;  /* Program application context */     
Widget parent;       /* highest level widget */
FStrng *LogDir;      /* Task Logfile directory name */
ObitErr *err;        /* Obit Error/message structure */
/* not used but for compatability with ObitView */
FStrng *FITS_dir;    /* FITSfile directory */
FStrng *log_dir;     /* logging directory */
Widget Display_shell;/* highest level widget */
#endif               /* end of declarations for main */

#ifndef OBITMESSMAIN
extern XtAppContext myApp;   /* Program application context */     
extern Widget parent;        /* highest level widget */
extern FStrng *LogDir;       /* Logfile directory name */
extern ObitErr *err;         /* Obit Error/message structure */
#endif /* end of declarations for other files */

#endif /* OBITMESS */ 
