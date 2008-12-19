/* $Id$  */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2004-2008                                          */
/*;  Associated Universities, Inc. Washington DC, USA.                */
/*;  This program is free software; you can redistribute it and/or    */
/*;  modify it under the terms of the GNU General Public License as   */
/*;  published by the Free Software Foundation; either version 2 of   */
/*;  the License, or (at your option) any later version.              */
/*;                                                                   */
/*;  This program is distributed in the hope that it will be useful,  */
/*;  but WITHOUT ANY WARRANTY; without even the implied warranty of   */
/*;  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the    */
/*;  GNU General Public License for more details.                     */
/*;                                                                   */
/*;  You should have received a copy of the GNU General Public        */
/*;  License along with this program; if not, write to the Free       */
/*;  Software Foundation, Inc., 675 Massachusetts Ave, Cambridge,     */
/*;  MA 02139, USA.                                                   */
/*;                                                                   */
/*;  Correspondence concerning Obit should be addressed as follows:   */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
/*  Define the basic components of the ObitImageMosaic structure      */
/*  This is intended to be included in a class structure definition   */
/**
 * \file ObitImageMosaicDef.h
 * ObitImageMosaic structure members for derived classes.
 */
#include "ObitDef.h"  /* Parent class definitions */
/** Threading info member object  */
ObitThread *thread;
/** Linked list of arrays of data.  */
ObitInfoList *info;
/** Number of images in Mosaic */
olong numberImages;
/** Number of images in Fly's eye */
olong nFlyEye;
/** Number of images already initialized */
olong nInit;
/** Image array */
ObitImage **images;
/** Full field image */
ObitImage *FullField;
/** Field of view as radius (deg) */
ofloat FOV;
/** Radius of the usable region of a given tile (cells) */
ofloat Radius;
/** Cell Spacing in X (deg) */
ofloat xCells;
/** Cell Spacing in Y (deg) */
ofloat yCells;
/** Requested size (CLEAN window) of outliers */
olong OutlierSize;
/** Number of pixels in X for each image */
olong *nx;
/** Number of pixels in Y for each image */
olong *ny;
/** Number of planes for each image */
olong *nplane;
/** RA shift (deg) for each image */
ofloat *RAShift;
/** Dec shift (deg) for each image */
ofloat *DecShift;
/** Are image OBIT_IO_AIPS or OBIT_IO_FITS? */
ObitIOType fileType;
/** Name of Mosaic images */
gchar *imName;
/** Class of Mosaic images, 2 char if AIPS */
gchar *imClass;
/** Disk numbers of Mosaic images */
olong *imDisk;
/** Sequence number of Mosaic images */
olong imSeq;
/** Is a full field image desired? */
gboolean doFull;
/** Restoring beam, values in deg*/
ofloat bmaj, bmin, bpa;
