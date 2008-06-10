/* $Id$ */
/* imagedisp (image display widget) header file */
/*-----------------------------------------------------------------------
*  Copyright (C) 1996,1999-2008
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
#include "glib.h"
#include "obitview.h"
#ifndef IMAGEDISP_H
#define IMAGEDISP_H
#define CONTROLWIDTH 150 /* width of control panel */
typedef struct {
  /* display information */
  XtAppContext   app;       /* program application context */
  Widget         parent;    /* main program window */
  Widget         shell;     /* highest level form */
  Widget         display;   /* form container widget */
  GC             gc;        /* graphics context */
  Widget         canvas;    /* canvas work widget */
  Widget         goober;    /* widget the size of the full window */
  Cursor         cursor;    /* Cursor for image display, Null = standard */
  XImage         *work;     /* ZPixmap for zooming current screen */
  olong           workNx;    /* X dimension of work ZPixmap  */
  olong           workNy;    /* Y dimension of work ZPixmap  */
  Widget         hscroll;   /* horizonal scroll scale */
  olong           hscr_vis;  /* 1=> hscroll visible, 0=not */
  olong           hscr_max;  /* hscroll maximum value   */
  olong           hscr_min;  /* hscroll minimum value   */
  olong           hscr_half; /* half width of display in x (im. pixels)*/
  Widget         vscroll;   /* vertical scroll scale */
  olong           vscr_vis;  /* 1=> vscroll visible, 0=not */
  olong           vscr_max;  /* hscroll maximum value   */
  olong           vscr_min;  /* hscroll minimum value   */
  olong           vscr_half; /* half width of display in y */
  olong           disp_wid;  /* display width */
  olong           disp_hei;  /* display height */
  
  /* image information */
  olong            zoom;      /* zoom factor, neg = zoom out */
  olong            scrollx;   /* center "x" pixel in display */
  olong            scrolly;   /* center "y" pixel in display */
  olong            iXCorn;    /* x corner pixel in display */
  olong            iYCorn;    /* y corner pixel in display */
  /* color mapping information */
  olong            depth;     /* number of bits in display */
  olong            ncolors;   /* number of colors in color table */
  olong            ngraph;    /* number of graphics overlay colors */
  olong            gcolor;    /* current color for graphics overlay */
  unsigned long  colut[MAXCOLOR+MAXGRAPH]; /* translation table for pixel values 
					      to color index, assumed in order */
  Colormap       cmap;       /* Windows color map */
  unsigned short red[MAXCOLOR+MAXGRAPH], 
    green[MAXCOLOR+MAXGRAPH], blue[MAXCOLOR+MAXGRAPH]; /* colortable */
  gboolean       trueColor;  /*true (1) if using a TrueColor visual */
  unsigned long  coltab[MAXCOLOR+MAXGRAPH]; /* pixel values for TrueColor */
  XColor         colors[MAXCOLOR+MAXGRAPH]; /* Current set of colors used */
  /* Control widgets */
  Widget         BriScroll; /* Brightness control scroll bar */
  Widget         ConScroll; /* Contrast control scroll bar */
  gshort         value[2];  /* Brightness and contrast scroll bar values */
  Widget         Info1;     /* First information line on control panel */
  Widget         Info2;     /* Second information line on control panel */
  Widget         Info3;     /* Third information line on control panel */
  Widget         Info4;     /* Fourth information line on control panel */
  Widget         Info5;     /* Fifth information line on control panel */
  Widget         Info6;     /* sixth information line on control panel */
  Widget         Info7;     /* seventhth information line on control panel */
  gboolean       showInfo;  /* if True show pixel information */
} ImageDisplay;

/* prototypes */
/* create initialize display */
ImageDisplay* MakeDisplay (Widget parent, Widget shell);
/* set display parameters for current image and force expose to redraw */
void ResetDisplay (ImageDisplay *IDdata);
/* set size and location of display after display and */
/* image sizes are known - sets up scroll bars   */
void SetDisplay (ImageDisplay* IDdata);
/* draw image */
void PaintImage (ImageDisplay* IDdata);
/* callbacks */
void FitPosCB (Widget w, XtPointer clientData, XtPointer callData);
void ImageResizeCB (Widget w, XtPointer clientData, XtPointer callData);
void Zoom25CB (Widget w, XtPointer clientData, XtPointer callData);
void Zoom50CB (Widget w, XtPointer clientData, XtPointer callData);
void Zoom100CB (Widget w, XtPointer clientData, XtPointer callData);
void Zoom200CB (Widget w, XtPointer clientData, XtPointer callData);
void Zoom400CB (Widget w, XtPointer clientData, XtPointer callData);
void Zoom800CB (Widget w, XtPointer clientData, XtPointer callData);
void Zoom1600CB (Widget w, XtPointer clientData, XtPointer callData);
void ButtonEH (Widget w, XtPointer data, XEvent *event);

#endif /* IMAGEDISP */ 
