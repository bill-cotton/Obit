/* $Id$  */
/* Display widget for images for ObitView allows zoom and scroll       */
/*-----------------------------------------------------------------------
*  Copyright (C) 1996,1997,1999,2002-2008
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
#include <Xm/ScrollBar.h>
#include <Xm/Form.h>
#include <X11/IntrinsicP.h>
#include <X11/cursorfont.h> 
#include <glib.h> 
#include "imagedisp.h"
#include "logger.h"
#include "menu.h"
#include "messagebox.h"
#include "poslabel.h"
#include "ObitSkyGeom.h"
#define min(a,b) (((a) < (b)) ? (a) : (b))
#define max(a,b) (((a) > (b)) ? (a) : (b))
#define SBAR_WIDTH 16; /* width of scroll bars */

/**
 *  \file imagedisp.c
 * Image Display Definitions.
 */

/*---------------Private function prototypes----------------*/
int  FitPos(ImageDisplay *IDdata);
void UpdateInfo (ImageDisplay *IDdata);
olong pfit (ofloat a[9][9], ofloat *s, ofloat dx[2], ofloat fblank);
int momnt (float ara[9][9], double x, double y, int nx, int ny, 
	   double momar[6],  float fblank);
int momnt2 (float ara[9][9], double *x, double *y, int nx, int ny, 
	   float fblank);
void matvmu (double m[3][3], double vi[3], double vo[3], int n);

/* progress/cancel box in FITS2Pix.c*/
void WorkingCursor(int on, int verbose, char* filename);


/**
 * Display zoomed image
 * \param IDdata   Image display
 * \param iXsize   horizontal (x) size of display area in pixels
 * \param iYsize   vertical (y) size of display area in pixels
 * \param iXpos    x pixel to put in bottom left corner (1-rel)
 * \param iYpos    x pixel to put in bottom left corner
 * \param Width    desired width of the region in the image (pixels)
 * \param Height   desired height of the region in the image (pixels)
 */
void ZoomDisplay (ImageDisplay *IDdata, olong iXsize, olong iYsize, olong iXpos, 
		  olong iYpos, olong Width, olong Height)
{
  olong xinc, yinc, xrep, yrep, ix, iy, iix=0, iiy, i, j;
  olong yaddr, addr=0, yaddr2, addr2, nxin, nxout, ival;
  XImage *work;
  unsigned long cval;
  gchar *pixarray, *gpharray, *work_data;
  
  /*  debug */
  /*
    Dimension cwid,chei;

    XtVaGetValues (IDdata->canvas,
    XmNwidth, &cwid,
    XmNheight, &chei[,
    NULL);
    fprintf (stderr," ZoomDisplay: canvas size: Width %d, Height %d\n",cwid,chei);
    fprintf (stderr," ZoomDisplay: iXpos %d, iYpos %d,  Width %d, Height %d\n",
    iXpos, iYpos, Width, Height);
  */
  
  /* don't bother if there's no image yet */
  if (!IDdata) return;
  
  /* is there a display yet? */
  if (!IDdata->canvas) return;
  if (!XtIsRealized (IDdata->canvas)) return;
  
  /* a few more validity checks */
  if (!image[CurImag].valid) return;
  if (!image[CurImag].pixarray) return;
  
  pixarray = image[CurImag].pixarray;
  gpharray = image[CurImag].gpharray;
  nxin = image[CurImag].myDesc->inaxes[0]; 
  work = IDdata->work;
  nxout = work->bytes_per_line;
  work_data = work->data;
  /* set increment, multiples in pixel array */
  if (IDdata->zoom>=1) /* no zoom or zoom in */
    { xinc = 1;
    xrep = IDdata->zoom;
    yinc = 1;
    yrep = IDdata->zoom;}
  else  /* zoom out */
    { xinc = -IDdata->zoom;
    xrep = 1;
    yinc = -IDdata->zoom;
    yrep = 1;}
  
  
  /* ix, iy are image pixel numbers, iix, iiy are display pixels */
  
  if (IDdata->depth==8) /* 8 bit display goes faster */
    {
      /* copy selected portion of image to display */
      if (IDdata->zoom<=1) /* no zooming or zoom out */
	{iiy = -1;
	for (iy=iYpos; iy<iYpos+Height; iy += yinc)  /* loop over columns */
	  {iiy++; iix = -1;
	  yaddr = iy * nxin; yaddr2 = iiy * nxout;
	  for (ix=iXpos; ix<iXpos+Width; ix += xinc)  /* loop down rows */
	    {iix++;
	    addr = yaddr + ix; addr2 = yaddr2 + iix;
	    /* Graphics or image? */
	    if ((*(gpharray + addr))>0) {
	      ival = (*(gpharray + addr)) + IDdata->ncolors-1;  /* graph */
	      cval = IDdata->coltab[ival];
	    } else cval = *(pixarray + addr);
	    *(work_data + addr2) = cval;
	    }
	  }
	}
      else  /* zooming */
	{iiy = -1;
	for (iy=iYpos; iy<iYpos+Height; iy++)  /* loop over columns */
	  {for (j=0; j<yrep; j++)
	    {iiy++; iix = -1;
	    yaddr = iy * nxin; yaddr2 = iiy * nxout;
	    for (ix=iXpos; ix<iXpos+Width; ix++)  /* loop down rows */
	      {addr = yaddr + ix; 
	      /* Graphics or image? */
	      if ((*(gpharray + addr))>0) {
		ival = (*(gpharray + addr)) + IDdata->ncolors-1;  /* graph */
		cval = IDdata->coltab[ival];
	      } else cval = *(pixarray + addr);
	      for (i=0; i<xrep; i++)
		{iix++;
		addr2 = yaddr2 + iix;
		*(work_data + addr2) = cval;
		}
	      }
	    }
	  }
	}
    }
  else  /* slowly for non 8-bit displays */
    {iiy = -1;
    for (iy=iYpos; iy<iYpos+Height; iy+=yinc)  /* loop over columns */
      {for (j=0; j<yrep; j++)
	{iiy++; iix = -1;
	for (ix=iXpos; ix<iXpos+Width; ix+=xinc)  /* loop down rows */
	  { for (i=0; i<xrep; i++)
	    {iix++;
	    /* Graphics or image? */
	    if ((*(gpharray + addr))>0) 
	      ival = (*(gpharray + addr)) + IDdata->ncolors-1;  /* graph */
	    else ival = *(pixarray + addr);
	    cval = IDdata->coltab[ival];
	    XPutPixel (work, iix, iiy, cval);
	    }
	  }
	}
      }
    }
  
  /* debug */
  /*fprintf (stderr,"ZoomDisplay: max address = %lu\n",maxaddr);*/
  
  /* write it out */
  XPutImage (XtDisplay (IDdata->canvas), XtWindow(IDdata->canvas), 
	     IDdata->gc, IDdata->work, 0, 0, 0, 0,
	     (unsigned int)iix, (unsigned int)iiy);
} /* end ZoomDisplay */

/**
 * Display zoomed 24 bit TrueColor image
 * \param IDdata   Image display
 * \param iXsize   horizontal (x) size of display area in pixels
 * \param iYsize   vertical (y) size of display area in pixels
 * \param iXpos    x pixel to put in bottom left corner (1-rel)
 * \param iYpos    x pixel to put in bottom left corner
 * \param Width    desired width of the region in the image (pixels)
 * \param Height   desired height of the region in the image (pixels)
 */
void ZoomDisplay24 (ImageDisplay *IDdata, int iXsize, int iYsize, int iXpos, 
		    int iYpos, int Width, int Height)
{
  olong xinc, yinc, xrep, yrep, ix, iy, iix=0, iiy, i, j;
  olong yaddr, addr, nxin, nxout, ival;
  unsigned long cval;
  XImage *work;
  gchar *pixarray, *gpharray, *work_data;
  
  /*  debug */
  /*  Dimension cwid,chei;
      XtVaGetValues (IDdata->canvas,
      XmNwidth, &cwid,
      XmNheight, &chei,
      NULL);
      fprintf (stderr," ZoomDisplay24: canvas size: Width %d, Height %d\n",cwid,chei);
  */
  
  /* don't bother if there's no image yet */
  if (!IDdata) return;
  
  /* is there a display yet? */
  if (!IDdata->canvas) return;
  if (!XtIsRealized (IDdata->canvas)) return;
  if (!IDdata->work) return;
  
  /* a few more validity checks */
  if (!image[CurImag].valid) return;
  if (!image[CurImag].pixarray) return;
  
  pixarray = image[CurImag].pixarray;
  gpharray = image[CurImag].gpharray;
  nxin = image[CurImag].myDesc->inaxes[0]; 
  work = IDdata->work;
  nxout = work->bytes_per_line;
  work_data = work->data;
  /* set increment, multiples in pixel array */
  if (IDdata->zoom>=1) /* no zoom or zoom in */
    { xinc = 1;
    xrep = IDdata->zoom;
    yinc = 1;
    yrep = IDdata->zoom;}
  else  /* zoom out */
    { xinc = -IDdata->zoom;
    xrep = 1;
    yinc = -IDdata->zoom;
    yrep = 1;}
  
  
  /* ix, iy are image pixel numbers, iix, iiy are display pixels */
  /* copy selected portion of image to display */
  
  if (IDdata->zoom<=1) /* no zooming or zoom out */
    {iiy = -1;
    for (iy=iYpos; iy<iYpos+Height; iy += yinc)  /* loop over columns */
      {iiy++; iix = -1;
      yaddr = iy * nxin;
      for (ix=iXpos; ix<iXpos+Width; ix += xinc)  /* loop down rows */
	{iix++;
	addr = yaddr + ix;
	/* Graphics or image? */
	if ((*(gpharray + addr))>0) 
	  ival = (*(gpharray + addr)) + IDdata->ncolors-1;  /* graph */
	else ival = (*(pixarray + addr));                   /* Image */
	XPutPixel (work, iix, iiy, IDdata->coltab[ival]);
	}
      }
    }
  else  /* zooming */
    {iiy = -1;
    for (iy=iYpos; iy<iYpos+Height; iy++)  /* loop over columns */
      {for (j=0; j<yrep; j++)
	{iiy++; iix = -1;
	yaddr = iy * nxin;
	for (ix=iXpos; ix<iXpos+Width; ix++)  /* loop down rows */
	  {addr = yaddr + ix; 
	  /* Graphics or image? */
	  if ((*(gpharray + addr))>0) 
	    ival = (*(gpharray + addr)) + IDdata->ncolors-1;  /* graph */
	  else ival = *(pixarray + addr);                     /* Image */
	  cval = IDdata->coltab[ival];
	  for (i=0; i<xrep; i++)
	    {iix++;
	    XPutPixel (work, iix, iiy, cval);
	    }
	  }
	}
      }
    }
  
  /* write it out */
  XPutImage (XtDisplay (IDdata->canvas), XtWindow(IDdata->canvas), 
	     IDdata->gc, IDdata->work, 0, 0, 0, 0,
	     (unsigned int)iix, (unsigned int)iiy);
} /* end ZoomDisplay24 */

/**
 * Set display parameters and ensures a single expose event to repaint the
 * image on the display.
 * \param IDdata   Image display
 */
void ResetDisplay (ImageDisplay *IDdata)
{
  olong old_wid, old_hei;
  
  /* save old display size */
  old_wid = IDdata->disp_wid;
  old_hei = IDdata->disp_hei;
  
  /* reset display parameters*/
  SetDisplay(IDdata);
  
  /* terrible hack - if SetDisplay made the display area larger then this
     will provoke an expose event.  If not then create one here. */
  if ((XtIsRealized (IDdata->canvas)) && (old_wid>=IDdata->disp_wid) && 
      (old_hei>=IDdata->disp_hei)) {
    PaintImage(IDdata); /* redraw */
  }
} /* end ResetDisplay */

/**
 * Sets size and position of the Image display window
 * \param IDdata   Image display
 */
void SetDisplay (ImageDisplay* IDdata) 
{
  olong iXHalf, iYHalf, inXZoom, inYZoom, iScrWid, iScrHei;
  Dimension cwid, chei;
  int it[5];
  olong iZoom;
  gchar *ZP = NULL;
  Display     *dpy = XtDisplay (IDdata->display);
  unsigned int unx, uny;
  unsigned long pasize, ilong;
  olong value, slider_size, increment, page_increment;
  olong sbar_width = SBAR_WIDTH; 
  
  /* update info in control pixel info area*/
  UpdateInfo(IDdata);  
  
  /* don't bother if there's no valid image yet */
  if (!IDdata) return; /* shouldn't happen */
  if (!image[CurImag].valid) { /* set display for no image */
    /* no scroll bars */
    if (XtIsRealized (IDdata->hscroll)) XtUnmapWidget (IDdata->hscroll);
    IDdata->hscr_vis = 0;
    if (XtIsRealized (IDdata->vscroll)) XtUnmapWidget (IDdata->vscroll);
    IDdata->vscr_vis = 0;
    /* no image - give default label */
    XtVaSetValues(IDdata->shell, XmNtitle, "ObitView", NULL);
    /* hide display */
    if (XtIsRealized (IDdata->canvas)) XtUnmapWidget (IDdata->canvas);
    return;}
  
  /*                    maximum display size */
  XtVaGetValues (IDdata->display,
		 XmNwidth,  &it[0],
		 XmNheight, &it[2],
		 NULL);

  cwid = (Dimension)it[0];
  chei = (Dimension)it[2];
  /*     cwid = cwid - CONTROLWIDTH;  leave room for control panel */
  IDdata->disp_wid = cwid;
  IDdata->disp_hei = chei;
  
  iZoom = IDdata->zoom; /* for convinence */
  
  /*         display should be an integral multiple of image pixels */
  if (iZoom>1)
    {IDdata->disp_wid =
       ((IDdata->disp_wid/iZoom)*iZoom);
    IDdata->disp_hei =
      ((IDdata->disp_hei/iZoom)*iZoom);}
  /*                    zoomed image size */
  inXZoom = image[CurImag].myDesc->inaxes[0];
  if (iZoom>1) inXZoom =  inXZoom * iZoom;
  if (iZoom<0) inXZoom = -inXZoom / iZoom;
  inYZoom = image[CurImag].myDesc->inaxes[1];
  if (iZoom>1) inYZoom =  inYZoom * iZoom;
  if (iZoom<0) inYZoom = -inYZoom / iZoom;
  /*                     scroll bar size */
  iScrWid = sbar_width;
  iScrHei = sbar_width;
  /*                      scroll bars needed? (iterate to get it right) */
  if ((IDdata->disp_wid-iScrWid)>=inXZoom) iScrHei = 0;
  if ((IDdata->disp_hei-iScrHei)>=inYZoom) iScrWid = 0;
  if ((IDdata->disp_wid-iScrWid)>=inXZoom) iScrHei = 0;
  if ((IDdata->disp_hei-iScrHei)>=inYZoom) iScrWid = 0;
  if (image[CurImag].valid) /* something in display? */
    /* Display needn't be larger than display+scrollbars */
    /* This sets the image size for no & negative zooms */
    {IDdata->disp_wid=min(IDdata->disp_wid, inXZoom+iScrWid);
    IDdata->disp_hei=min(IDdata->disp_hei,inYZoom+iScrHei);}
  /*                      correct for size of scroll bars */
  IDdata->disp_wid -= iScrWid; 
  if (IDdata->disp_wid<0) IDdata->disp_wid = iScrWid;
  IDdata->disp_hei -= iScrHei; 
  if (IDdata->disp_hei<0) IDdata->disp_hei = iScrHei;
  /*   Display may still be too big */
  if (image[CurImag].valid) /* something in display? */
    {IDdata->disp_wid = min (IDdata->disp_wid, inXZoom);
    IDdata->disp_hei = min (IDdata->disp_hei, inYZoom);}
  else
    {IDdata->disp_wid += iScrWid;  /* no scroll bars if no image */
    IDdata->disp_hei += iScrHei;} 
  /* leave at least the width of the scroll bars (SBAR_WIDTH) around the edge */
  IDdata->disp_wid = min (IDdata->disp_wid, (olong)cwid-sbar_width);
  IDdata->disp_hei = min (IDdata->disp_hei, (olong)chei-sbar_width);
  /*                    display should have an even number  of rows */
  /*     IDdata->disp_hei = 2 * ((IDdata->disp_hei+1)/2); */
  
  /* resize window */
  if (XtIsRealized (IDdata->canvas)) XtMapWidget (IDdata->canvas);
  XtResizeWidget(IDdata->canvas, 
		 (Dimension)IDdata->disp_wid,
		 (Dimension)IDdata->disp_hei, 
		 (Dimension)0);
  
  /*                      Half Size of display in image pixels */
  iXHalf = IDdata->disp_wid/2;
  if (iZoom>1) iXHalf =  iXHalf / iZoom;
  if (iZoom<0) iXHalf = -iXHalf * iZoom;
  iYHalf = IDdata->disp_hei/2;
  if (iZoom>1) iYHalf =  iYHalf / iZoom;
  if (iZoom<0) iYHalf = -iYHalf * iZoom;
  iXHalf = min (iXHalf, image[CurImag].myDesc->inaxes[0]/2);
  iYHalf = min (iYHalf, image[CurImag].myDesc->inaxes[1]/2);
  
  /* setup and center scroll */
  /*     IDdata->scrollx = image[CurImag].myDesc->inaxes[0] / 2; */
  /*     IDdata->scrolly = image[CurImag].myDesc->inaxes[1] / 2; */
  IDdata->hscr_max = image[CurImag].myDesc->inaxes[0]-2*iXHalf;
  IDdata->hscr_min = 1;
  if (IDdata->hscr_max<=IDdata->hscr_min) 
    IDdata->hscr_max = IDdata->hscr_min + 1;
  IDdata->hscr_half = iXHalf;
  value = IDdata->scrollx - iXHalf;
  slider_size = 2 * iXHalf;
  increment = iXHalf / 5; if (increment<1) increment = 1;
  page_increment = 3 * iXHalf / 2; if (page_increment<1) page_increment = 1;
  if (value>IDdata->hscr_max) value = IDdata->hscr_max;
  if (value<IDdata->hscr_min) value = IDdata->hscr_min;
  if (iScrHei)
    {
      /* keep X-Windows from blowing it's tiny mind */
      XmScrollBarSetValues (IDdata->hscroll, 1, 1, 1, 1, False);
      XtVaSetValues(IDdata->hscroll, 
		    XmNheight,    iScrHei,
		    XmNvalue,     value,
		    XmNmaximum,   (Dimension)(IDdata->hscr_max+2*iXHalf),
		    XmNminimum,   (Dimension)(IDdata->hscr_min),
		    NULL);
      XmScrollBarSetValues (IDdata->hscroll, value, slider_size, increment, 
			    page_increment, False);}
  IDdata->vscr_max = image[CurImag].myDesc->inaxes[1]-2*iYHalf;
  IDdata->vscr_min = 1;
  if (IDdata->vscr_max<=IDdata->vscr_min) 
    IDdata->vscr_max = IDdata->vscr_min + 1;
  IDdata->vscr_half = iYHalf;
  value = IDdata->scrolly - iYHalf;
  slider_size = 2 * iYHalf;
  increment = iYHalf / 5; if (increment<1) increment = 1;
  page_increment = 3 * iYHalf / 2; if (page_increment<1) page_increment = 1;
  if (value>IDdata->vscr_max) value = IDdata->vscr_max;
  if (value<IDdata->vscr_min) value = IDdata->vscr_min;
  if (iScrWid)
    {
      /* keep X-Windows from blowing it's tiny mind */
      XmScrollBarSetValues (IDdata->vscroll, 1, 1, 1, 1, False);
      XtVaSetValues(IDdata->vscroll, 
		    XmNwidth,     iScrWid,
		    XmNvalue,     value,
		    XmNmaximum,   (Dimension)(IDdata->vscr_max+2*iYHalf),
		    XmNminimum,   (Dimension)(IDdata->vscr_min),
		    NULL);
      XmScrollBarSetValues (IDdata->vscroll, value, slider_size, increment, 
			    page_increment, False);}
  
  /*     make horizonal scroll bar visible or invisible as necessary */
  /*  iScrHei = 0 => no horizontal scroll */
  if (iScrHei) /* visible */
    {XtMapWidget (IDdata->hscroll);
    IDdata->hscr_vis = 1;}
  else /* invisible */
    {XtUnmapWidget (IDdata->hscroll);
    IDdata->hscr_vis = 0;}
  
  /*     make vertical scroll bar visible or invisible as necessary */
  /*  iScrWid = 0 => no vertical scroll */
  if (iScrWid) /* visible */
    {XtMapWidget (IDdata->vscroll);
    IDdata->vscr_vis = 1;}
  else /* invisible */
    {XtUnmapWidget (IDdata->vscroll);
    IDdata->vscr_vis = 0;}
  
  /* make work ZPixmap for display*/
  /* delete old if necessary*/
  unx = IDdata->disp_wid;
  uny = IDdata->disp_hei;
  if ( !IDdata->work || 
       (unx>IDdata->work->width) || (uny>IDdata->work->height)) 
    {
      if (IDdata->work) 
	{if (IDdata->work->data) g_free(IDdata->work->data);  IDdata->work->data=NULL; 
	IDdata->work->data = NULL;
	XtFree((XtPointer)IDdata->work);
	IDdata->work = NULL;} 
      
      /* create new pix map */
      pasize = (unx+10) * uny * ((IDdata->depth+1)/8);
      /* 24 bit displays are different */
      if (IDdata->depth>12 ) pasize = unx * uny * 4;
      ZP = (char*) g_malloc (pasize);
      IDdata->work = 
	XCreateImage (dpy, 
		      DefaultVisual(dpy, DefaultScreen(dpy)), 
		      IDdata->depth, ZPixmap, 0, 
		      ZP, unx, uny, 32, 0); 
      /* blank it out */
      for (ilong=0; ilong<pasize; ilong++) 
	IDdata->work->data[ilong] = 0;
    } /* end of create new ZPixmap */
} /* end SetDisplay */

/**
 * Redraw the image on the canvas computing zoom and scroll
 * \param IDdata   Image display
 */
void PaintImage (ImageDisplay* IDdata)
{
  olong nXDIB, nYDIB, iZoom, iSrcWid, iSrcHei, iXHalf, iYHalf;
  olong iXSize, iYSize, iXPos, iYPos, iXPage, iYPage;
  Dimension cWidth, cHeight;
  olong i, ch1;
  gchar TitleString[501], *cptr;

 /* blank if there's no image yet */
  if (!IDdata) return; /* this shouldn't happen */
  if (!image[CurImag].valid) {
    if (XtIsRealized (IDdata->canvas)) 
      XClearArea (XtDisplay (IDdata->canvas), XtWindow (IDdata->canvas), 
		  0, 0, 0, 0, TRUE);
    return;
  }
  if (!image[CurImag].pixarray) return;
  
  /* reset window title to file name */
  if (image[CurImag].valid && image[CurImag].FileName) { /* title for image only */
    cptr = image[CurImag].FileName->sp;
    ch1 = 0;
    for (i=0; i<image[CurImag].FileName->length;i++)
      if (cptr[i]=='/') ch1 = i+1;
    if (image[CurImag].myDesc->inaxes[2]>1)
      {g_snprintf (TitleString,500, "%s, plane %d",
		&cptr[ch1],
		image[CurImag].PlaneNo+1);}
    else /* one plane */
      {g_snprintf (TitleString,500,"%s", &cptr[ch1]);}
    XtVaSetValues(IDdata->shell, 
		  XmNtitle,   TitleString,
		  NULL);
  } /* end of label */
  
  cWidth  = IDdata->disp_wid;   /* current display size */
  cHeight = IDdata->disp_hei;
  
  nXDIB = image[CurImag].myDesc->inaxes[0];
  nYDIB = image[CurImag].myDesc->inaxes[1];
  iZoom = IDdata->zoom;
  if (iZoom == 0) iZoom = 1;
  /* size of display in image pixels */
  iXHalf = IDdata->disp_wid;  /* divide by two later */
  if (iZoom>1) iXHalf =  iXHalf / iZoom;
  if (iZoom<0) iXHalf = -iXHalf * iZoom;
  iYHalf = IDdata->disp_hei; /* divide by two later */
  if (iZoom>1) iYHalf =  iYHalf / iZoom;
  if (iZoom<0) iYHalf = -iYHalf * iZoom;
  iXHalf = min (iXHalf, nXDIB);
  iYHalf = min (iYHalf, nYDIB);
  iSrcWid = iXHalf;
  iSrcHei = iYHalf;
  iXHalf = iXHalf / 2;
  iYHalf = iYHalf / 2;
  /* Size of display area */
  iXSize = IDdata->disp_wid;
  iYSize = IDdata->disp_hei;
  iXSize = min (iXSize, nXDIB);
  iYSize = min (iYSize, nYDIB);
  if (iZoom>1) {iYSize=iYSize*iZoom; iXSize=iXSize*iZoom;}

  /*  for negative zooms, iDisp* is set in SetDisplay */
  iXSize = min (iXSize, (olong)cWidth);
  iYSize = min (iYSize, (olong)cHeight);
  iXSize = max (iXSize, 1);
  iYSize = max (iYSize, 1);
  /* "page" size for scrolling */
  iXPage = iSrcWid;
  iYPage = iSrcHei;
  iXPos = IDdata->scrollx - iXHalf;
  iYPos = IDdata->scrolly - iYHalf;
  iXPos = max (iXPos, 0);
  iYPos = max (iYPos, 0);
  iXPos = min (iXPos, image[CurImag].myDesc->inaxes[0]-iXHalf*2);
  iYPos = min (iYPos, image[CurImag].myDesc->inaxes[1]-iYHalf*2);
  IDdata->iXCorn = iXPos;
  IDdata->iYCorn = iYPos;

  /* copy Pixarray to display */
  if (IDdata->trueColor) /* slow, 24 bit color */
    ZoomDisplay24 (IDdata, iXSize, iYSize, iXPos, iYPos, iSrcWid, iSrcHei);
  else /* fast, 8 bit color */
    ZoomDisplay (IDdata, iXSize, iYSize, iXPos, iYPos, iSrcWid, iSrcHei);
  
} /* end PaintImage */

/**
 * Callback for motion in horizonal scroll
 * \param w           widget activated
 * \param clientData  client data
 * \param callData    call data
 */
void hscrollCB (Widget w, XtPointer clientData, XtPointer callData)
{
  ImageDisplay *IDdata = (ImageDisplay *)clientData;
  XmScrollBarCallbackStruct *call_data = 
    (XmScrollBarCallbackStruct *) callData;
  
  /* read value of scrollbar */
  IDdata->scrollx = call_data->value+IDdata->hscr_half;
  PaintImage((ImageDisplay*)clientData); /* redraw */
} /* end hscrollCB */

/**
 * Callback for motion in vertical scroll
 * \param w           scroll widget activated
 * \param clientData  client data
 * \param callData    call data
 */
void vscrollCB (Widget w, XtPointer clientData, XtPointer callData)
{
  ImageDisplay *IDdata = (ImageDisplay *)clientData;
  XmScrollBarCallbackStruct *call_data = 
    (XmScrollBarCallbackStruct *) callData;
  /* read value of scrollbar */
  IDdata->scrolly = call_data->value+IDdata->vscr_half;
  /* reset scroll position and redraw */
  PaintImage((ImageDisplay*)clientData);
} /* end vscrollCB */

/**
 * Callback for redisplay image canvas
 * \param w           canvas widget activated
 * \param clientData  client data
 * \param callData    call data
 */
void   ImageRedispCB (Widget w, XtPointer clientData, XtPointer callData)
{
  ImageDisplay *IDdata = (ImageDisplay *)clientData;
  XmDrawingAreaCallbackStruct *cb = (XmDrawingAreaCallbackStruct *)callData;
  XExposeEvent  *event = (XExposeEvent *) cb->event;
  
  /* NOP unless widget is canvas */
  if (!IDdata) return;
  if (!IDdata->canvas) return;
  if (w!=IDdata->canvas) return;
  
  /*fprintf (stderr,"RedispCB: x %d, y %d, nx %d, ny %d, count %d\n",
    event->x, event->y, event->width, event->height, event->count);*/
  
  /* return if there are more expose events in the queue */
  if (event->count) return;
  
  /* reset display */
  SetDisplay(IDdata);
  /* redraw exposed area from ImageData */
  PaintImage (IDdata);
} /* end of ImageRedispCB */

/**
 * Callback for resize image canvas
 * \param w           canvas widget activated
 * \param clientData  client data
 * \param callData    call data
 */
void   ImageResizeCB (Widget w, XtPointer clientData, XtPointer callData)
{
  ImageDisplay  *IDdata = (ImageDisplay *)clientData;
  
  if (!IDdata) return;
  if (!IDdata->canvas) return;
  
  /*fprintf(stderr,"ImageResizeCB: window resized\n "); debug */
  
  /* reset display */
  ResetDisplay(IDdata);
  /* redraw exposed area from ImageData */
  /*    PaintImage (IDdata); */
} /* end of ImageResizeCB */

/* event handlers */
/**
 * Event handler for button press in image display window
 * \param w      canvas widget activated
 * \param data   event data
 * \param event  event
 */
void ButtonEH (Widget w, XtPointer data, XEvent *event)
{
  olong         ix, iy, isOK=0;
  ofloat       fx, fy, z;
  ImageDisplay *IDdata = (ImageDisplay*)data;
  
  if (!IDdata) return;
  if (!image[CurImag].valid) return;
  
  /* determine image pixel */
  fx = (ofloat)event->xbutton.x;
  fy = (ofloat)event->xbutton.y;
  z = IDdata->zoom;
  if (z>1)  /* zoomed? Take care of sub pixel precision*/
    {fx = (fx - (z*0.5) + 0.5) / z;
    fy = (fy - (z*0.5) + 0.5) / z;}
  if (z<0) 
    {fx = -fx * z; fy = -fy * z; }
  fx = fx + IDdata->iXCorn;
  fy = image[CurImag].myDesc->inaxes[1] - (fy + IDdata->iYCorn) - 1.0;
  ix = (olong)(fx+0.5); iy = (olong)(fy+0.5);
  IDdata->showInfo = 1; /* have something to show */
  if (event->xbutton.button==1) /* left button - display pixel value*/
    { image[CurImag].fXpixel = fx;
    image[CurImag].fYpixel = fy;
    image[CurImag].iXPixel = ix;
    image[CurImag].iYPixel = iy;
    image[CurImag].iFitted = 0; 
    UpdateInfo(IDdata);  /* write info in control info area*/
    isOK = 1;
    } /* end Left button */
  
  if (event->xbutton.button==3) /* right button - fit position */
    { image[CurImag].iXPixel = ix;
    image[CurImag].iYPixel = iy;
    image[CurImag].iFitted = FitPos(IDdata);
    UpdateInfo(IDdata);  /* write info in control info area*/
    isOK = image[CurImag].iFitted==1;
    } /* end Right button */
  if (doLog &&  isOK) DoLogger(); /* log if necessary */
} /* end ButtonEH */

/**
 * Create image display widget
 * \param parent  Parent widget
 * \param shell   Shell widget
 * \return the image display widget
 */
ImageDisplay* MakeDisplay (Widget parent, Widget shell)
{
  ImageDisplay *IDdata;     /* display data structure */
  olong          i;
  Dimension     cwid, chei;
  int           it[5];
  olong          sbar_width = SBAR_WIDTH;
  
  IDdata = (ImageDisplay*) XtMalloc (sizeof (ImageDisplay));
  
  /* set size */
  /*                    maximum display size */
  XtVaGetValues (parent,
		 XmNwidth,  &it[0],
		 XmNheight, &it[2],
		 NULL);
  cwid = (Dimension)it[0];
  chei = (Dimension)it[2];
  cwid = cwid - CONTROLWIDTH; /* leave room for control panel */
  /*  chei = chei - 40; */
  
  IDdata->disp_wid = cwid;
  IDdata->disp_hei = chei;
  
  /* save main window */
  IDdata->parent = parent;
  
  /* set top level form */
  IDdata->shell = shell;
  
  /* make Form widget for display */
  IDdata->display = XtVaCreateManagedWidget ("display", xmFormWidgetClass,
					     parent,
					     /* add extra pixels for scroll bars */
					     XmNwidth,     IDdata->disp_wid+sbar_width,
					     XmNheight,    IDdata->disp_hei+sbar_width,
					     XmNtopAttachment,  XmATTACH_FORM,
					     XmNrightAttachment,  XmATTACH_FORM,
					     NULL);
  
  /* drawing canvas */
  IDdata->canvas = 
    XtVaCreateManagedWidget ("canvas", xmDrawingAreaWidgetClass, 
			     IDdata->display, 
			     /*			     XmNwidth,     IDdata->disp_wid-sbar_width,
						     XmNheight,    IDdata->disp_hei-sbar_width,*/
			     XmNwidth,     1,
			     XmNheight,    1,
			     XmNtopAttachment,  XmATTACH_FORM,
			     XmNleftAttachment,  XmATTACH_FORM,
			     NULL);
  /*   Add callbacks to handle resize and exposures.  */
  XtAddCallback (IDdata->canvas, XmNexposeCallback, ImageRedispCB, 
		 (XtPointer)IDdata);
  
  /* trap button press */
  XtAddEventHandler (IDdata->canvas, ButtonPressMask, FALSE, 
		     (XtEventHandler)ButtonEH, 
		     (XtPointer)IDdata);
  
  /* add horizonal scroll bar ("scale" in X-Windows) */
  
  IDdata->hscroll = 
    XtVaCreateManagedWidget ("hscroll", xmScrollBarWidgetClass, 
			     IDdata->display,
			     XmNwidth,     IDdata->disp_wid,
			     XmNheight,     sbar_width,
			     XmNmaximum,   IDdata->disp_wid,
			     XmNminimum,            1,
			     XmNvalue,              1,
			     XmNshowValue,       False,
			     XmNorientation, XmHORIZONTAL,
			     XmNprocessingDirection, XmMAX_ON_RIGHT,
			     /*				       XmNtopAttachment, XmATTACH_WIDGET,
							       XmNtopWidget,     IDdata->canvas,*/
			     XmNbottomAttachment, XmATTACH_FORM,
			     XmNleftAttachment,  XmATTACH_FORM,
			     XmNrightAttachment,  XmATTACH_FORM,
			     XmNrightOffset,     sbar_width,
			     NULL);
  IDdata->hscr_vis = 1;  /* it's visible */
  /* add call backs */
  XtAddCallback(IDdata->hscroll, XmNvalueChangedCallback, hscrollCB, 
		(XtPointer) IDdata);
  /*  XtAddCallback(IDdata->hscroll, XmNdragCallback, hscrollCB, 
      (XtPointer) IDdata); */
  
  /* add vertical scroll bar */
  
  IDdata->vscroll = 
    XtVaCreateManagedWidget ("vscroll", xmScrollBarWidgetClass, 
			     IDdata->display,
			     XmNheight,     IDdata->disp_hei,
			     XmNwidth,          sbar_width,
			     XmNmaximum,   IDdata->disp_hei,
			     XmNminimum,            1,
			     XmNvalue,              1,
			     XmNshowValue,       False,
			     XmNorientation, XmVERTICAL,
			     XmNprocessingDirection, XmMAX_ON_BOTTOM,
			     XmNtopAttachment,  XmATTACH_FORM,
			     XmNbottomAttachment,  XmATTACH_FORM,
			     XmNbottomOffset,     sbar_width,
			     XmNrightAttachment,  XmATTACH_FORM,
			     /*				       XmNleftAttachment, XmATTACH_WIDGET,
							       XmNleftWidget,   IDdata->canvas,*/
			     NULL);
  IDdata->vscr_vis = 1;  /* it's visible */
  /* add call back */
  XtAddCallback(IDdata->vscroll, XmNvalueChangedCallback, vscrollCB, 
		(XtPointer) IDdata);
  /*  XtAddCallback(IDdata->vscroll, XmNdragCallback, vscrollCB, 
      (XtPointer) */
  
  /* make dummy widget the size of the display which can detect when 
     it's size has changed */
  IDdata->goober = 
    XtVaCreateManagedWidget ("goober", xmDrawingAreaWidgetClass, 
			     IDdata->display, 
			     XmNtopAttachment,  XmATTACH_FORM,
			     XmNleftAttachment,  XmATTACH_FORM,
			     XmNbottomAttachment,  XmATTACH_FORM,
			     XmNrightAttachment,  XmATTACH_FORM,
			     NULL);
  /* go boom  XtUnmapWidget (IDdata->goober);  make invisible */
  XtAddCallback (IDdata->goober, XmNresizeCallback, ImageResizeCB,
		 (XtPointer)IDdata);
  /* number of colors */
  IDdata->ncolors = XDisplayCells (XtDisplay (IDdata->display), 
				   XDefaultScreen (XtDisplay (IDdata->display)));
  if (IDdata->ncolors>MAXCOLOR) IDdata->ncolors = MAXCOLOR;

  /* Number of graphics overlays */
  IDdata->ngraph = MAXGRAPH;
  IDdata->gcolor = 1;  /* initial graphics overlay color */

  /* init color arrays */
  for (i=0; i<MAXCOLOR; i++) 
    {IDdata->colut[i]= i;  /* need ramp for True color */
    IDdata->coltab[i]= i; /* need ramp for True color */
    IDdata->red[i]   = 0;
    IDdata->green[i] = 0;
    IDdata->blue[i]  = 0;}
  
  /* bit planes in display */
  IDdata->depth = XDisplayPlanes (XtDisplay (IDdata->display), 
				  XDefaultScreen (XtDisplay (IDdata->display)));
  /* context */
  IDdata->gc = XCreateGC (XtDisplay (IDdata->display), 
			  DefaultRootWindow (XtDisplay(IDdata->display)),
			  0, NULL );
  /* init image */
  IDdata->work    = NULL;
  IDdata->zoom    = 1;
  IDdata->scrollx = 0;
  IDdata->scrolly = 0;
  IDdata->iXCorn  = 0;
  IDdata->iYCorn  = 0;
  /* no color map yet */
  /*??  IDdata->cmap = 0;*/
  
  /* no slider controls yet */
  IDdata->BriScroll = NULL;
  IDdata->ConScroll = NULL;
  IDdata->value[0]  = 0;
  IDdata->value[1]  = 0;
  IDdata->Info1     = NULL;
  IDdata->Info2     = NULL;
  IDdata->Info3     = NULL;
  IDdata->Info4     = NULL;
  IDdata->Info5     = NULL;
  IDdata->Info6     = NULL;
  IDdata->Info7     = NULL;
  IDdata->showInfo  = 0;
  
  /* return Image display structure */
  return IDdata;
  
} /* end MakeDisplay */

/**
 * Callback routine to fit position 
 * \param w           widget activated
 * \param clientData  client data
 * \param callData    call data
 */
void FitPosCB (Widget w, XtPointer clientData, XtPointer callData)
{
  gboolean     isOK;
  ImageDisplay *IDdata = (ImageDisplay *)clientData;
  
  /* use last position selected */
  image[CurImag].iFitted = FitPos(IDdata);
  UpdateInfo(IDdata);  /* write info in control info area*/
  isOK = image[CurImag].iFitted==1;
  if (doLog &&  isOK) DoLogger(); /* log if necessary */
  FitPos(IDdata);
} /* End FitPosCB */

/**
 * routine to fit a point near the position of the cursor:
 * IDdata->iXPixel, iYPixel
 * Sets Dis[iDs] values fXpixel, fYpixel, fBpixel
 * \param IDdata   Image display
 * \return 1 if worked otherwise 0, -1 => fit failed.
 */
int FitPos(ImageDisplay *IDdata)
{ 
  olong iXcen, iYcen, iX, iY, iXp, iYp;
  ofloat data[9][9], fblank;
  olong pos[7] = {0, 0, 0, 0, 0, 0, 0}, ndim;
  ofloat s, dx[2], pixval=0.0, *pixP;
  Display *dpy = XtDisplay(IDdata->canvas);
 
  /* validity checks */
  if (!image[CurImag].valid) return 0;
  if (!image[CurImag].myDesc) return 0;

  ndim=image[CurImag].myDesc->naxis;
  iXcen = image[CurImag].iXPixel;
  iYcen = image[CurImag].iYPixel;
  fblank = ObitMagicF();
  /* default return values */
  image[CurImag].fBpixel = 0.;
  image[CurImag].fXpixel = 0.;
  image[CurImag].fYpixel = 0.;
  
  if (!image[CurImag].valid) return 0;
  
  /* use 9x9 values around center */
  pos[2] = image[CurImag].PlaneNo;
  for (iY=0; iY<9; iY++) {
    iYp = iYcen + iY - 4;
    pos[1] = iYp;
    for (iX=0; iX<9; iX++) {
      iXp = iXcen + iX - 4;
      pos[0] = iXp;
      /* get value */
      pixP =  ObitFArrayIndex (image[CurImag].myPixels, pos);
      if (pixP) pixval = *pixP;  /* In image? */
      else pixval = fblank;
      if (pixval!=fblank)
	data[iY][iX] = pixval;
      else  /* blank */
	data [iY][iX] = fblank;
    }
  }  /*  end of loop loading image data in array */
  if (pfit (data, &s, dx, fblank))
    {/*fprintf (stderr, "fit failed\n");*/
      XBell(dpy, 50); 
      return -1;}  /* fit failed*/
  image[CurImag].fBpixel = s;
  image[CurImag].fXpixel = (ofloat)iXcen + dx[0];
  image[CurImag].fYpixel = (ofloat)iYcen + dx[1];
  return 1;
} /* end fit pos */

/**
 * Update information about selected pixel in control window
 * \param IDdata   Image display
 */
void UpdateInfo (ImageDisplay *IDdata)
{
  olong          ndim;
  olong          i, posOK;
  XmString      wierdstring = NULL;
  gchar         jerkstring[101];
  Arg           wargs[5]; 
  ofloat        fblank, val, pix[7], *pixP;
  gchar         axtype[3][9], label[3][20];
  odouble       pos[3];
  gboolean      valid;
  olong         ipos[7] = {0, 0, 0, 0, 0, 0, 0};
  Display *dpy = XtDisplay(IDdata->display);
  static Cursor cursor=(Cursor)NULL;
  ObitImageDesc *desc;

  /* check that info defined for current image */
  if (!image[CurImag].valid)  IDdata->showInfo = 0;
  
  /* need to blank display? */
  if (!IDdata->showInfo)
    {
      wierdstring = XmStringCreateSimple ("    ");
      XtSetArg (wargs[0], XmNlabelString, wierdstring);
      XtSetValues (IDdata->Info1, wargs, 1);
      XtSetValues (IDdata->Info2, wargs, 1);
      XtSetValues (IDdata->Info3, wargs, 1);
      XtSetValues (IDdata->Info4, wargs, 1);
      XtSetValues (IDdata->Info5, wargs, 1);
      XtSetValues (IDdata->Info6, wargs, 1);
      XtSetValues (IDdata->Info7, wargs, 1);
      if (wierdstring) XmStringFree(wierdstring); wierdstring = NULL;
      return;
    } /* end blank info area */
  /* return OK if no valid image */
  if (!image[CurImag].valid) return;
  ndim = image[CurImag].myDesc->naxis;  /* Number of dimensions */
  /* make sure hourglass cursor initialized */
  if (!cursor) cursor = XCreateFontCursor(dpy, XC_watch);
  /*  fitted or current */
  if (ndim>2)
    g_snprintf (jerkstring, 100, "(%7.2f,%7.2f, %d)",
	     image[CurImag].fXpixel+1.0, image[CurImag].fYpixel+1.0,
	     image[CurImag].PlaneNo+1);
  else
    g_snprintf (jerkstring, 100, "(%7.2f,%7.2f)",
	     image[CurImag].fXpixel+1.0, image[CurImag].fYpixel+1.0);
  wierdstring = XmStringCreateSimple (jerkstring);
  XtSetArg (wargs[0], XmNlabelString, wierdstring);
  XtSetValues (IDdata->Info1, wargs, 1);
  if (wierdstring) XmStringFree(wierdstring); wierdstring = NULL;
  
  /*  get flux */
  fblank = ObitMagicF();
  /*  fitted or current? */
  if (image[CurImag].iFitted>0) /* fitted values */
    {valid = 1;
    if (image[CurImag].fBpixel==fblank)
      /* blanked */
      g_snprintf (jerkstring, 100, "pixel blanked");
    else
      g_snprintf (jerkstring, 100, "value=%f", image[CurImag].fBpixel);}
  else  /* current position */
    {val = 0;
    for (i=3; i<7; i++) ipos[i] = 0;
    /* Note: Image class pixel numbers are 0 rel; geometry is 1 rel. */
    ipos[0] = image[CurImag].iXPixel; ipos[1] = image[CurImag].iYPixel;
    ipos[2] = image[CurImag].PlaneNo;
    pixP =  ObitFArrayIndex (image[CurImag].myPixels, ipos);
    if (pixP) val = *pixP;
    else val = fblank;
    valid = (ndim>1);
    g_snprintf (jerkstring, 100, "invalid pixel");
    if (valid) {
      if (val==fblank)
	/* blanked */
	g_snprintf (jerkstring, 100, "pixel blanked");
      else
	g_snprintf (jerkstring, 100, "value=%f", val);
    }
    }
  /* write second line (flux density) */
  jerkstring[16] = 0; /* limit size of string */
  wierdstring = XmStringCreateSimple (jerkstring);
  XtSetArg (wargs[0], XmNlabelString, wierdstring);
  XtSetValues (IDdata->Info2, wargs, 1);
  if (wierdstring) XmStringFree(wierdstring); wierdstring = NULL;
  
  /* celestial position */
  /* equinox is line 3 */
  if ((usr_equinox>0.0) && (image[CurImag].myDesc->equinox>0.0))
    g_snprintf(jerkstring,100, "Equinox %6.1f", usr_equinox);
  else if (image[CurImag].myDesc->equinox>0.0)
    g_snprintf(jerkstring,100,"Equinox %6.1f", image[CurImag].myDesc->equinox);
  else
    g_snprintf(jerkstring,100,"Equinox unknown");
  wierdstring = XmStringCreateSimple (jerkstring);
  XtSetArg (wargs[0], XmNlabelString, wierdstring);
  XtSetValues (IDdata->Info3, wargs, 1);
  if (wierdstring) XmStringFree(wierdstring); wierdstring = NULL;
  
  /*  Get position */
  pix[0] = image[CurImag].fXpixel+1.0; 
  pix[1] = image[CurImag].fYpixel+1.0;
  pix[2] = (ofloat)image[CurImag].PlaneNo+1.0;
  if (valid) {
    strncpy(axtype[0], image[CurImag].myDesc->ctype[0], 8);
    strncpy(axtype[1], image[CurImag].myDesc->ctype[1], 8);
    strncpy(axtype[2], image[CurImag].myDesc->ctype[2], 8);
    axtype[0][8] = 0; axtype[1][8] = 0; axtype[2][8] = 0;

    /* Obit position stuff */
    desc = image[CurImag].myDesc;
    posOK = 
      !ObitSkyGeomWorldPos(pix[0], pix[1], 
			   desc->crval[desc->jlocr], desc->crval[desc->jlocd],
			   desc->crpix[desc->jlocr], desc->crpix[desc->jlocd],
			   desc->cdelt[desc->jlocr], desc->cdelt[desc->jlocd],
			   desc->crota[desc->jlocd], &desc->ctype[desc->jlocr][4],
			   &pos[0], &pos[1]);
    /* Need to precess? */
    if ((usr_equinox > 0.0) && (usr_equinox != image[CurImag].myDesc->equinox))
      {if (usr_equinox==1950.0) ObitSkyGeomBtoJ (&pos[0], &pos[1]);
      if (usr_equinox==2000.0) ObitSkyGeomJtoB (&pos[0], &pos[1]);}

    /* Third axis - assume linear - descriptor has been sub imaged -> plane 1*/
    pos[2] =  desc->crval[2] + desc->cdelt[2]*(1.0 - desc->crpix[2]);

    AxisLabel(pos[0], axtype[0], label[0]);  /* human readable */
    if (ndim>=2) AxisLabel(pos[1], axtype[1], label[1]);
    if (ndim>=3) AxisLabel(pos[2], axtype[2], label[2]);
    if (posOK) {  /* valid position */
      /* write fourth line (first axis) */
      wierdstring = 
	XmStringCreateSimple (label[0]);
      XtSetArg (wargs[0], XmNlabelString, wierdstring);
      XtSetValues (IDdata->Info4, wargs, 1);
      if (wierdstring) XmStringFree(wierdstring); wierdstring = NULL;
      
      /* write fifth line (second axis) */
      if (ndim>=2) 
	{wierdstring = 
	   XmStringCreateSimple (label[1]);
	XtSetArg (wargs[0], XmNlabelString, wierdstring);
	XtSetValues (IDdata->Info5, wargs, 1);
	if (wierdstring) XmStringFree(wierdstring); wierdstring = NULL;}
      
      /* write sixth line (third axis) */
      if (ndim>=3) 
	{wierdstring = 
	   XmStringCreateSimple (label[2]);
	XtSetArg (wargs[0], XmNlabelString, wierdstring);
	XtSetValues (IDdata->Info6, wargs, 1);
	if (wierdstring) XmStringFree(wierdstring); wierdstring = NULL;}
    }}
  else {  /* invalid position */
    g_snprintf (jerkstring, 100, "invalid pixel");
    /* write third line (invalid pixel message) */
    wierdstring = 
      XmStringCreateSimple (jerkstring);
    XtSetArg (wargs[0], XmNlabelString, wierdstring);
    XtSetValues (IDdata->Info3, wargs, 1);
    if (wierdstring) XmStringFree(wierdstring); wierdstring = NULL;
  }  /* end of give position on invalid */
  /*  fitted or current? */
  if (image[CurImag].iFitted>0) /* fitted values */
    g_snprintf (jerkstring, 100, "fitted");
  else if (image[CurImag].iFitted<0) /* fit failed */
    g_snprintf (jerkstring, 100, "fit failed!");
  else /* no fit */
    g_snprintf (jerkstring, 100, " ");
  
  /* write seventh line  */
  wierdstring = 
    XmStringCreateSimple (jerkstring);
  XtSetArg (wargs[0], XmNlabelString, wierdstring);
  XtSetValues (IDdata->Info7, wargs, 1);
  if (wierdstring) XmStringFree(wierdstring); wierdstring = NULL;
} /* end UpdateInfo */

/* call backs for zoom */
/**
 * Callback for 25% zoom
 * \param w           scroll widget activated
 * \param clientData  client data
 * \param callData    call data
 */
void Zoom25CB (Widget w, XtPointer clientData, XtPointer callData)
{
  ImageDisplay  *IDdata = (ImageDisplay *)clientData;
  
  /* NOP if same zoom */
  if (IDdata->zoom==-4) return;
  
  IDdata->zoom = -4; /* set new zoom factor */
  
  /* reset display */
  ResetDisplay(IDdata);
  /* mark in menu */
  MenuMarkZoom (0);
} /* end of Zoom25CB */

/**
 * Callback for 50% zoom
 * \param w           scroll widget activated
 * \param clientData  client data
 * \param callData    call data
 */
void Zoom50CB (Widget w, XtPointer clientData, XtPointer callData)
{
  ImageDisplay  *IDdata = (ImageDisplay *)clientData;
  
  /* NOP if same zoom */
  if (IDdata->zoom==-2) return;
  
  IDdata->zoom = -2; /* set new zoom factor */
  
  /* reset display */
  ResetDisplay(IDdata);
  /* mark in menu */
  MenuMarkZoom (1);
} /* end of Zoom50CB */

/**
 * Callback for 100% zoom
 * \param w           scroll widget activated
 * \param clientData  client data
 * \param callData    call data
 */
void Zoom100CB (Widget w, XtPointer clientData, XtPointer callData)
{
  ImageDisplay  *IDdata = (ImageDisplay *)clientData;
  
  /* NOP if same zoom */
  if (IDdata->zoom==1) return;
  
  IDdata->zoom = 1; /* set new zoom factor */
  
  /* reset display */
  ResetDisplay(IDdata);
  /* mark in menu */
  MenuMarkZoom (2);
} /* end of Zoom100CB */

/**
 * Callback for 200% zoom
 * \param w           scroll widget activated
 * \param clientData  client data
 * \param callData    call data
 */
void Zoom200CB (Widget w, XtPointer clientData, XtPointer callData)
{
  ImageDisplay  *IDdata = (ImageDisplay *)clientData;
  
  /* NOP if same zoom */
  if (IDdata->zoom==2) return;
  
  IDdata->zoom = 2; /* set new zoom factor */
  
  /* reset display */
  ResetDisplay(IDdata);
  /* mark in menu */
  MenuMarkZoom (3);
} /* end of Zoom200CB */

/**
 * Callback for 400% zoom
 * \param w           scroll widget activated
 * \param clientData  client data
 * \param callData    call data
 */
void Zoom400CB (Widget w, XtPointer clientData, XtPointer callData)
{
  ImageDisplay  *IDdata = (ImageDisplay *)clientData;
  
  /* NOP if same zoom */
  if (IDdata->zoom==4) return;
  
  IDdata->zoom = 4; /* set new zoom factor */
  
  /* reset display */
  ResetDisplay(IDdata);
  /* mark in menu */
  MenuMarkZoom (4);
} /* end of Zoom400CB */

/**
 * Callback for 800% zoom
 * \param w           scroll widget activated
 * \param clientData  client data
 * \param callData    call data
 */
void Zoom800CB (Widget w, XtPointer clientData, XtPointer callData)
{
  ImageDisplay  *IDdata = (ImageDisplay *)clientData;
  
  /* NOP if same zoom */
  if (IDdata->zoom==8) return;
  
  IDdata->zoom = 8; /* set new zoom factor */
  
  /* reset display */
  ResetDisplay(IDdata);
  /* mark in menu */
  MenuMarkZoom (5);
} /* end of Zoom800CB */

/**
 * Callback for 1600% zoom
 * \param w           scroll widget activated
 * \param clientData  client data
 * \param callData    call data
 */
void Zoom1600CB (Widget w, XtPointer clientData, XtPointer callData)
{
  ImageDisplay  *IDdata = (ImageDisplay *)clientData;
  
  /* NOP if same zoom */
  if (IDdata->zoom==16) return;
  
  IDdata->zoom = 16; /* set new zoom factor */
  
  /* reset display */
  ResetDisplay(IDdata);
  /* mark in menu */
  MenuMarkZoom (6);
} /* end of Zoom1600CB */

int momnt (float ara[9][9], double x, double y, int nx, int ny, 
	   double momar[6], float fblank)
/*-----------------------------------------------------------------------*/
/*  Calculate all 0th, 1st, and 2nd moments of a nx*ny subarray of ara,  */
/*  centered at x,y.  nx and ny should be odd.                           */
/*  Inputs:                                                              */
/*     ara     Input data array                                          */
/*     x       x-center for moment calculation (1-rel)                   */
/*     y       y-center                                                  */
/*     nx      # of points to include in x-direction. nx                 */
/*             should be odd.  The points will be centered               */
/*             about x (rounded)                                         */
/*     ny      # of points in y-direction                                */
/*     fblank  magic blanking value                                      */
/*  Outputs:                                                             */
/*     momar   00,10,20,01,11,02 yx-moments of ara                       */
/*     return  0 => o.k.                                                 */
/*             1 => subarray doesn't fit in main array                   */
/*-----------------------------------------------------------------------*/
{
  long  ind, i, j, nj, indx, indy, iax1, iax2, iay1, iay2;
  float t;
  double s, arg, prod;

  /*      compute loop limits (1-rel)   */
  i = x + 0.5;
  iax1 = i - nx/2;
  iax2 = iax1 + nx - 1;
  i = y + 0.5;
  iay1 = i - ny/2;
  iay2 = iay1+ ny - 1;
  /*      check loop limits             */
  if ((iax1<1) || (iax2>9) || (iay1<1) || (iay2>9)) return 1;
  /*      compute moments               */
  ind = 0;
  for (i = 1; i<=3; i++) {
    nj = 4 - i;
    for (j=1; j<=nj; j++) {
      ind++;
      s = 0.0;
      for (indx=iax1; indx<=iax2; indx++) {
	for (indy=iay1; indy<=iay2; indy++) {
	  t = ara[indy-1][indx-1];
	  if (t!=fblank) {
	    if (i>1) {
	      arg = indx - x;
	      prod = pow (arg,(i-1.0));
	      t *= prod;}     /* ((indx - x)**(i-1)); */
	    if (j>1) {
	      arg = indy - y;
	      prod = pow (arg,(j-1.0));
	      t *= prod;}  /* ((indy - y)**(j-1));*/
	    s += t;}
	}
      }
      momar[ind-1] = s;
    }
  }
    return 0;
  }  /* end of momnt */

int momnt2 (float ara[9][9], double *x, double *y, int nx, int ny, 
	   float fblank)
/*-----------------------------------------------------------------------*/
/*  Compute position from moments                                        */
/*  centered at x,y.  nx and ny should be odd.                           */
/*  Inputs:                                                              */
/*     ara     Input data array                                          */
/*     x       initial x-center for moment calculation (1-rel)           */
/*     y       initial y-center                                          */
/*     nx      # of points to include in x-direction. nx                 */
/*             should be odd.  The points will be centered               */
/*             about x (rounded)                                         */
/*     ny      # of points in y-direction                                */
/*     fblank  magic blanking value                                      */
/*  Outputs:                                                             */
/*     x       fitted x-center for moment calculation (1-rel)            */
/*     y       fitted y-center                                           */
/*     return  0 => o.k.                                                 */
/*             1 => subarray doesn't fit in main array                   */
/*-----------------------------------------------------------------------*/
{
  long  i, indx, indy, iax1, iax2, iay1, iay2;
  float t;
  double sum, sumx, sumy;

  /* compute loop limits (1-rel)   */
  i = *x + 0.5;
  iax1 = i - nx/2;
  iax2 = iax1 + nx - 1;
  i = *y + 0.5;
  iay1 = i - ny/2;
  iay2 = iay1+ ny - 1;

  /*  check loop limits */
  if ((iax1<1) || (iax2>9) || (iay1<1) || (iay2>9)) return 1;

  /* compute moments */
  sum = sumx = sumy = 0.0;
  for (indx=iax1; indx<=iax2; indx++) {
    for (indy=iay1; indy<=iay2; indy++) {
      t = ara[indy-1][indx-1];
      if (t!=fblank) {
	sum  += t;
	sumx += t * (indx - *x);
	sumy += t * (indy - *y);
      }
    } /* end y loop */
  } /* end x loop */

  if (sum<=0.0) return 1;
  /* Correct previous estimation */
  *x += sumx / sum;
  *y += sumy / sum;
  return 0;
}  /* end of momnt2 */

void matvmu (double m[3][3], double vi[3], double vo[3], int n)
/*-----------------------------------------------------------------------*/
/*  Matrix-vector multiplication  vo = vi * m                            */
/*  Inputs:                                                              */
/*     m       Input matrix                                              */
/*     vi      Input vector                                              */
/*     n       Array dimension                                           */
/*  Outputs:                                                             */
/*     vo      Output vector                                             */
/*-----------------------------------------------------------------------*/
{
  int  i, j;
  double s;

  for (i=0; i<n; i++) {
    s = 0.0;
    for (j=0; j<n; j++) s = s + m[j][i] * vi[j];
    vo[i] = s;
  }
}  /* end of matvmu */

int pfit (float a[9][9], float *s, float dx[2], float fblank)
/*--------------------------------------------------------------------*/
/*  Make a parabolic least-squares fit to a 9x9 matrix about the      */
/*  peak value in an array and determine strength and position        */
/*  of maximum.                                                       */
/*   Inputs:                                                          */
/*     a       Data input array[dec][RA]                              */
/*     fblank  Value for blanked pixel                                */
/*  Outputs:                                                          */
/*     dx  Position of max relative to a[5][5]                        */
/*     s   Strength of max                                            */
/*     return  0=OK else failed                                       */
/*--------------------------------------------------------------------*/
{
   double  absmax, minpix, x, y, temp[6], momar[6], d;
   double mat[3][3] = {{0.55555, -0.33333, -0.33333}, {-0.33333, 0.5, 0.0},
		     {-0.33333, 0.0, 0.5}};
   int ix, iy, ixmax, iymax;
   /* default return values */
   dx[0] = 0;
   dx[1] = 0;
   *s = a[4][4];

   /*  find peak in array    */
   absmax = 0.0; ixmax = -1; iymax = -1; minpix=1.0e20;
   for (ix=0; ix<9; ix++) {
     for (iy=0; iy<9; iy++) {
       if ((a[iy][ix]!=fblank) && (fabs(a[iy][ix])>absmax)) 
	 {absmax=fabs(a[iy][ix]); ixmax = ix; iymax = iy;}
       if ((a[iy][ix]!=fblank) && (a[iy][ix]<minpix)) 
	 minpix = a[iy][ix];
     }
   }

   /* Subtract minimum */
   for (ix=0; ix<9; ix++) {
     for (iy=0; iy<9; iy++) {
       if (a[iy][ix]!=fblank) a[iy][ix] -= minpix;
     }
   }
   
   /* check for valid data */
   if ((ixmax<0) || (iymax<0)) return 1;
   /*  00, 01, 02, 10, 11, 20 */
   x = ixmax+1.0;
   y = iymax+1.0;
   /*  default values       */
   dx[0] = x - 5.0;
   dx[1] = y - 5.0;
   *s = a[iymax][ixmax];
   if (momnt (a, x, y, 3, 3, momar, fblank)) return 1;

   /*  multiply matrix * even moms  yields const & quadratic terms */
   temp[0] = momar[0];
   temp[1] = momar[2];
   temp[2] = momar[5];
   matvmu (mat, temp, &temp[3], 3);

   /*  pick up linear & cross term  */
   temp[0] = momar[1] / 6.;
   temp[1] = momar[3] / 6.;
   temp[2] = momar[4] / 4.;
   
   /*  offset of peak */
   d = 4.* temp[4] * temp[5] - (temp[2]*temp[2]);
   if (d==0.0) return 2;
   dx[0] = (temp[2]*temp[0] - 2.*temp[1]*temp[4]) / d;
   dx[1] = (temp[2]*temp[1] - 2.*temp[0]*temp[5]) / d;

   /*  value of peak */
   *s = temp[3] + dx[0]*(temp[1] + dx[0]*temp[5]
			 + dx[1]*temp[2]) + dx[1]*(temp[0]+dx[1]*temp[4]);
   dx[0] += x - 5.0;  /* correct wrt center of input array */
   dx[1] += y - 5.0;

   /* Check that in bounds */
   if ((fabs(dx[0])>5.0) || (fabs(dx[1])>5.0)) {
     /* Use more primitive method */
     if (momnt2(a, &x, &y, 3, 3, fblank)!=0) return 1;
     if (momnt2(a, &x, &y, 3, 3, fblank)!=0) return 1;
     dx[0] = x - 5.0;  /* correct wrt center of input array */
     dx[1] = y - 5.0;
     *s = absmax;
 }
   return 0;  /* OK */
} /* end of pfit */
