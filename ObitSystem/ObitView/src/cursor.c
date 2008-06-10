/* $Id$  */
/* creates cursor for ObitView window */
/* tries to use a 16x16 cursor; if this isn't appropriate then the standard */
/* cursor is returned                                                       */
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
#include <X11/X.h> 
#include <Xm/Xm.h> 
#include <Xm/DrawingA.h> 
#include <Xm/MainW.h> 
#include "definecursor.h"

/**
 *  \file cursorc
 * cursor related functions.
 */

/**
 * Make custom cursor for image display window if possible else return NULL
 * \param dosplay     X display
 * \param screen      X screen
 */
Cursor MakeImageCursor (Display *display, Drawable screen)
{
  unsigned int width, height;
  Pixmap csource, cmask;
  int xhot = ImageCursor_x_hot, yhot = ImageCursor_y_hot;
  XColor  cfore, cback;
  Cursor  cursor;
  
  /* see if this display wants a 16x16 cursor, if not return Null */
  XQueryBestCursor (display, screen, ImageCursor_width, ImageCursor_height, 
		    &width, &height);
  if ((width!=ImageCursor_width) || (height!=ImageCursor_height)) return 0;
  
  /* set colors (should do this from resource file)*/
  /* foreground */
  cfore.red=0; cfore.green=65535; cfore.blue = 65535;
  cfore.flags = DoRed | DoGreen | DoBlue;
  /* background */
  cback.red=65535; cback.green=65535; cback.blue = 0;
  cback.flags = DoRed | DoGreen | DoBlue;
  
  /* make pixmaps */
  csource = XCreateBitmapFromData(display, screen, 
				  (const char *)ImageCursor_bits,
				  width, height);
  if (!csource) return 0;
  cmask = XCreateBitmapFromData(display, screen, 
				(const char *)ImageCursorMask_bits,
				width, height);
  if (!cmask) return 0;
  
  /* make cursor */
  cursor = XCreatePixmapCursor (display, csource, cmask, &cfore, &cback, 
				xhot, yhot);
  /* delete pixmaps */
  XFreePixmap(display, csource);
  XFreePixmap(display, cmask);
  return cursor;
} /* end MakeImageCursor */
