/* $Id$  */
/* function prototypes for AIPSfilebox.c */
/*-----------------------------------------------------------------------
*  Copyright (C) 2013
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
#ifndef AIPSFILEBOX_H
#define AIPSFILEBOX_H
#include "obitview.h"

/* Communication structure */
typedef struct {
  Widget dialog, Form, Open, Cancel, Help, SList, dir;
  Widget userT, nameT, classT, seqT;
  ImageDisplay *IDdata; /* Image display info */
  olong cno;            /* image catalog number */
  olong AUser;          /* aips user number */
  olong ASeq;           /* image sequence */
  FStrng *AName;        /* Image AIPS name */
  FStrng *AClass;       /* Image AIPS class */
  FStrng *AIPS_dir;     /* Image AIPS Directory name */
  olong cancel;         /* cancel requested if != 0 */
  gshort  xpos, ypos;   /* location of the box */
  int selPos;           /* Element selected */
} AIPSFileBoxStuff;

AIPSFileBoxStuff* AIPSFileBox (Widget w, ImageData *image);
#endif  /* AIPSFILEBOX_H */
