/* $Id$ */
/* ObitView header file */
/*-----------------------------------------------------------------------
*  Copyright (C) 2005-2013
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
#include <Xm/DrawingA.h> 
#include <Xm/MainW.h> 
#include <stdlib.h>
#include <stdio.h>
#include "glib.h"
#include "FStrng.h"
#include "ObitErr.h"
#include "ObitTypes.h"
#include "ObitImageDesc.h"
#include "ObitDConCleanWindow.h"
#include "ObitFArray.h"
#include "ObitInfoList.h"
#ifndef OBITVIEW
#define OBITVIEW
#define MAXCOLOR 128     /* number of colors in display  
			  also defined in histo.h */
#define MAXGRAPH 6       /* Number of graphics colors */


typedef struct {
/* FITS stuff */
  gboolean valid;    /* If true then a valid file is connected */
  gboolean reLoad;   /* If true then a (re) load is needed */
  ObitImageDesc *myDesc;   /* Image descriptor */
  ObitFArray    *myPixels; /* Image pixels */
  ObitDConCleanWindow *edtWindow;  /* Clean window being edited */
  gchar   *pixarray; /* image pixel array, 1 byte per pixel */
                     /* ordered by rows then column */
  gchar   *gpharray; /* image graphs array, 1 byte per pixel */
                     /* ordered by rows then column */
  olong   nxArray, nyArray; /* Dimensions of pixarray and gpharray */
  ObitIOType DataType;/* Type of Image, OBIT_IO_FITS, OBIT_IO_AIPS */
  FStrng *FileName;  /* FITS file name */
  FStrng *AIPSName;  /* AIPS name */
  FStrng *AIPSClass; /* AIPS class */
  FStrng *AIPSDir;   /* AIPS directory name */
  olong AIPSseq;      /* AIPS sequence number */
  olong AIPSuser;     /* AIPS user number */
  olong Field;        /* Field number (1-rel) of mosaic */
  olong NField;       /* Number of fields in mosaic */
  ofloat usr_equinox;/* User specified equinox */
  ofloat maxVal;     /* Maximum pixel value in plane */
  ofloat minVal;     /* Minimum pixel value in plane */
  olong iXPixel;     /* Cursor X pixel (in image) */
  olong iYPixel;     /* Cursor Y pixel (in image) */
  gboolean iFitted;  /* TRUE if position fitted else FALSE */
  ofloat fXpixel;    /* fitted X pixel number 1-rel */
  ofloat fYpixel;    /* fitted Y pixel number 1-rel*/
  ofloat fBpixel;    /* fitted peak brightness */
  ofloat PixRange[2];/* range of pixel values to display */
  olong   mapFunc;   /* mapping: 0=>linear, 1=>sqrt, 2=>histo. Eq. */
  olong PlaneNo;     /* image plane to display 0 rel */
  olong hiDim[3];    /* image dimensions 4-6 to display 0 rel */
} ImageData;

/* global data structures */
#ifdef OBITVIEWMAIN
XtAppContext myApp;  /* Program application context */     
ImageData image[2];  /* images */
gshort  CurImag;     /* which image is current 0 or 1 */
FStrng *FITS_dir;    /* FITSfile directory */
FStrng *mark_dir;    /* Mark position directory */
FStrng *log_dir;     /* logging directory */
gboolean doLog;      /* it true position logging turned on */
ofloat  usr_equinox; /* Equinox desired by user */
Widget Display_shell;/* highest level widget */
ObitInfoList *requestList; /* InfoList for requests to return to client */
ObitErr *err;        /* Obit Error/message structure */
#endif               /* end of declarations for main */

#ifndef OBITVIEWMAIN
extern XtAppContext myApp;  /* Program application context */     
extern ImageData image[2];  /* images */
extern gshort  CurImag;     /* which image is current 0 or 1 */
extern FStrng *FITS_dir;    /* FITSfile directory */
extern FStrng *mark_dir;    /* Mark position directory */
extern FStrng *log_dir;     /* logging directory */
extern gboolean doLog;      /* if true position logging turned on */
extern ofloat  usr_equinox; /* Equinox desired by user */
extern Widget Display_shell;/* highest level widget */
extern ObitInfoList *requestList; /* InfoList for requests to return to client */
extern ObitErr *err;        /* Obit Error/message structure */
#endif /* end of declarations for other files */

#endif /* OBITVIEW */ 
