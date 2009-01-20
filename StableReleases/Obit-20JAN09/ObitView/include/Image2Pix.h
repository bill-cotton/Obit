/* $Id$ */
/* Header file for FITS2Pix routine */
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
#ifndef IMAGE2PIX_H 
#define IMAGE2PIX_H 
#include <Xm/Xm.h> 
#include <glib.h> 
#include "obitview.h"

olong Image2Pix (ImageData *image, ImageDisplay *IDdata, gboolean verbose);
#endif /*  IMAGE2PIX_H  */ 
