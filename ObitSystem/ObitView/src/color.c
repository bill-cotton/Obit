/* $Id$ */
/*    Colormap functions for ObitView */
/*-----------------------------------------------------------------------
*  Copyright (C) 1996,1998,1999, 2002-2026
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
#include <Xm/Scale.h> 
#include <stdlib.h>
#include <stdio.h>
#include "imagedisp.h"
#include "color.h"
#include "messagebox.h"

/**
 *  \file color.c
 * manages display colors
 */
/*---------------Private function prototypes----------------*/
void SetColor1(ImageDisplay* IDdata);
void SetColor2(ImageDisplay *IDdata);
void RevColor(ImageDisplay *IDdata);
void UpdateColor (ImageDisplay *IDdata);
void RGB2Pixel(XColor *colors, ImageDisplay *IDdata);

/*---------------Public functions ----------------*/

/**
 * Callback routine Color contour colors 
 * \param w           widget activated
 * \param clientData  client data
 * \param callData    call data
 */
void ColContourCB (Widget w, XtPointer clientData, XtPointer callData)
{
  ImageDisplay *IDdata = (ImageDisplay *)clientData;
  SetColor1(IDdata);     /* set new colors */
  SetColorTable(IDdata); /* install new colors */
} /* end ColContourCB */

/**
 * Callback routine forPseudo Flame colors
 * \param w           widget activated
 * \param clientData  client data
 * \param callData    call data
 */
void PhlameCB (Widget w, XtPointer clientData, XtPointer callData)
{
  ImageDisplay *IDdata = (ImageDisplay *)clientData;
  SetColor2(IDdata);     /* set new colors */
  SetColorTable(IDdata); /* install new colors */
} /* end PhlameCB */

/**
 * Callback routine for Gray scale
 * \param w           widget activated
 * \param clientData  client data
 * \param callData    call data
 */
void GrayscaleCB (Widget w, XtPointer clientData, XtPointer callData)
{
  ImageDisplay *IDdata = (ImageDisplay *)clientData;
  InitColorTable(IDdata);     /* set new colors */
  SetColorTable(IDdata); /* install new colors */
} /* end GrayscaleCB */

/**
 * Callback routine for Reverse color table c
 * \param w           widget activated
 * \param clientData  client data
 * \param callData    call data
 */
void ReverseColCB (Widget w, XtPointer clientData, XtPointer callData)
{
  ImageDisplay *IDdata = (ImageDisplay *)clientData;
  RevColor(IDdata);     /* set new colors */
  SetColorTable(IDdata); /* install new colors */
} /* end ReverseColCB */

/**
 * Callback routine for reset color table
 * \param w           widget activated
 * \param clientData  client data
 * \param callData    call data
 */
void ResetColCB (Widget w, XtPointer clientData, XtPointer callData)
{
  ImageDisplay *IDdata = (ImageDisplay *)clientData;
  IDdata->value[0] = 128; IDdata->value[1]=128;
  /* reset sliders */
  XmScaleSetValue(IDdata->BriScroll, 128);
  XmScaleSetValue(IDdata->ConScroll, 128);
  InitColorTable(IDdata);     /* set new colors */
  SetColorTable(IDdata);      /* install new colors */
} /* end ResetColCB */

/**
 *  Create and install a color map - fill with grayscale ramp
 *  If an 8 bit PseudoColor visual is available, it is used;
 *  else if a 24 bit TrueColor visual is available it is used
 *  instead of a colormap
 * \param shell    shell widget
 * \param IDdata   Image display
 */
void SetupColorMap (Widget shell, ImageDisplay *IDdata)
{
  olong        i, ncolor;
  Display     *dpy = XtDisplay (shell);
  XColor      *Colors;
  XVisualInfo visual_info;
  
  /* Initialize colors */
  for (i = 0; i < IDdata->ncolors; i++) IDdata->colors[i].pixel = -1;
  
  /* allocate temp color array */
  ncolor = XDisplayCells (dpy, XDefaultScreen (dpy));
  Colors = (XColor*) XtMalloc (sizeof (XColor) * ncolor);
  
  /* debug: */
  /*
    if (!XMatchVisualInfo (dpy, XDefaultScreen (dpy), 24, TrueColor, 
    &visual_info))
    {
    MessageShow ("No TrueColor(24 bit) visual available");
    } else {
    MessageShow ("A TrueColor(24 bit) visual is available, will use it.");
    MessageShow ("Note: an 8 bit PseudoColor visual would run faster."); 
    IDdata->trueColor = 1;
    if (!XAllocColorCells(dpy, DefaultColormap (dpy, DefaultScreen(dpy)), 
    FALSE, NULL, 0, IDdata->colut,  
    IDdata->ncolors)) {
    MessageShow ("Allocating colors failed");
    return;
    }
    for (i = 0; i < IDdata->ncolors; i++) {
    IDdata->coltab[i] = IDdata->colut[i];
    IDdata->colors[i].pixel = IDdata->colut[i];
    }
    if (Colors) {XtFree ((XtPointer) Colors);} Colors = NULL;
    return;
    } 
  */
  
  /* Need true color? Have to repaint the screen.*/
  if (IDdata->trueColor) {
    RGB2Pixel(Colors, IDdata); /* compute pixel values */
    if (Colors) {XtFree ((XtPointer) Colors);} Colors = NULL; /* no longer need this */
    ResetColCB (IDdata->display, (XtPointer)IDdata, NULL);
    return;
  } /* end true color */
  /* end debug*/
  
  
  /* The first time, create a colormap and install it in the
     canvas widget's window. Also set up the window manager
     colormap windows list so both colormaps get installed,
     if the system is capable of handling multiple colormaps. */
  
  
  if (!(IDdata->cmap || IDdata->trueColor)) {
    IDdata->cmap = DefaultColormap (dpy, DefaultScreen(dpy));
    if (!XAllocColorCells(dpy, IDdata->cmap, FALSE, NULL, 0, IDdata->colut, 
			  IDdata->ncolors)) {
      /* could not allocate colors - try a new color table */
      /* first see if there is a pseudocolor visual available */
      
      /* debug */
      /*
	fprintf (stderr,"SetupColorMap: must allocate colormap\n");
      */
      if (!XMatchVisualInfo (dpy, XDefaultScreen (dpy), 8, PseudoColor, 
			     &visual_info)) {
	if (!XMatchVisualInfo (dpy, XDefaultScreen (dpy), 24, TrueColor, 
			       &visual_info))
	  {  /* Oh bother */
	    MessageShow ("No PseudoColor(8 bit) or TrueColor(24 bit) visual available, I cannot cope!");
	    if (Colors) {XtFree ((XtPointer) Colors);} Colors = NULL;
	    return;
	  }
	/* MessageShow ("A TrueColor(24 bit) visual is available, will use it."); */
	/* MessageShow ("Note: an 8 bit PseudoColor visual would run faster."); */
	IDdata->trueColor = 1;  /* remember */
      } /* end 24 bit TrueColor section */
      
      /* create new color map if necesary */
      if (!IDdata->trueColor) {
	MessageShow ("Colormap full, create new one");
	IDdata->cmap = 
	  XCreateColormap (dpy, XtWindow(IDdata->canvas), visual_info.visual, AllocAll);
	if (!IDdata->cmap) {
	  MessageShow ("Could not create colormap");
	  if (Colors) {XtFree ((XtPointer) Colors);} Colors = NULL;
	  return;}
	/* Copy as much of the default color table as possible */
	for (i = 0; i < ncolor; i++) {
	  Colors[i].pixel = i;
	  Colors[i].flags = DoRed | DoGreen | DoBlue;
	}
	XQueryColors (dpy, DefaultColormap(dpy, DefaultScreen(dpy)),
		      Colors, ncolor);
	/* debug */
	/*
	  fprintf (stderr,"SetupColorMap: trueColor= %d cmap=%d \n",
	  IDdata->trueColor, IDdata->cmap);
	*/
	
	XStoreColors (dpy, IDdata->cmap, Colors, ncolor);
	/* install color map */
	XSetWindowColormap (dpy, XtWindow (shell), IDdata->cmap);
	/* let window see changes in colormap */
	/*ugly results      XSelectInput (dpy, XtWindow (shell), ColormapChangeMask);*/
	/* use middle of color map */
	for (i=0; i<IDdata->ncolors;i++) 
	  /* 64 colors       IDdata->colut[i] = i+ncolor/2-IDdata->ncolors;*/
	  /* 128 colors IDdata->colut[i] = i+((3*ncolor)/4)-IDdata->ncolors; */
	  IDdata->colut[i] = i+((5*ncolor)/8)-IDdata->ncolors; /* 256 colors */
      } /* end of setup new virtual colormap */
    } /* end of install virtual colormap */
  } /* end of allocation/installation of colormap */
  
  /* fill with ramp in gray */
  InitColorTable (IDdata);
  
  /* first color always black */
  Colors[0].pixel = IDdata->colut[0];
  Colors[0].flags = DoRed|DoGreen|DoBlue;
  Colors[0].red   = Colors[0].blue = Colors[0].green =  BlackPixel(dpy, 0);
  IDdata->red[0] = Colors[0].red; 
  IDdata->blue[0] = Colors[0].blue; 
  IDdata->green[0] = Colors[0].green;
  
  /* rest of image color table */
  for (i = 0; i < IDdata->ncolors; i++) {
    Colors[i].pixel = IDdata->colut[i];
    Colors[i].flags = DoRed|DoGreen|DoBlue;
    Colors[i].red   = IDdata->red[i];
    Colors[i].blue  = IDdata->blue[i];
    Colors[i].green = IDdata->green[i];
    IDdata->red[i] = Colors[i].red; 
    IDdata->blue[i] = Colors[i].blue; 
    IDdata->green[i] = Colors[i].green;
  } /* end of fill color table loop */
  
  /* allocate colors */
  if (IDdata->trueColor)
    RGB2Pixel (Colors, IDdata);
  else
    XStoreColors (dpy, IDdata->cmap, Colors, IDdata->ncolors);
  
  if (Colors) {XtFree ((XtPointer) Colors);} Colors = NULL;
  
  /* reset sliders */
  IDdata->value[0] = 128; IDdata->value[1]=128;
  if (IDdata->BriScroll) XmScaleSetValue(IDdata->BriScroll, 128);
  if (IDdata->ConScroll) XmScaleSetValue(IDdata->ConScroll, 128);

  /* Init graphics color */
  InitGraphColor(IDdata);
} /* end of SetupColorMap */

/**
 *  Routine to change the color table
 * IDdata->value[0] = Brightness control (0 - 255) 
 * IDdata->value[1] = Contrast control (0- 255) 
 * IDdata->ncolors  = number of image colors (e.g. 128) 
 * IDdata->red      = default red color table (0 - 65535) 
 * IDdata->green    = default green color table (0 - 65535) 
 * IDdata->blue     = default blue color table (0 - 65535) 
 * 
 * The new color table is the entries from the default table with  
 * 2(IDdata->value[0]-128) added to each and  then multiplied by  
 * 1 + 5*((IDdata->value[1]-128)/128); 
 * \param IDdata   Image display
 */
void SetColorTable (ImageDisplay* IDdata)
{
  olong i;
  olong lColor;
  ofloat scale, offset;
  XColor *colors;
  Display     *dpy = XtDisplay (IDdata->canvas);

  /* Effects of brightness and contrast controls */
  scale = 1.0 + 5.0*((IDdata->value[1]-128.0)/128.0);
  offset=2.0*(IDdata->value[0]-128.0); 
  
  /* allocate color array */
  colors = (XColor*) XtMalloc (sizeof (XColor) * (IDdata->ncolors));
  
  /* replace table / update colormap; reserve color 0=(0,0,0) for blanked */
  for (i=0; i<IDdata->ncolors; i++) {
    lColor = scale * ((ofloat)i + offset) + 0.5;
    if (lColor>=IDdata->ncolors) lColor = IDdata->ncolors-1; 
    if (lColor<0) lColor = 0;
    /* black is reserved for blanked */
    if ((i>0) && (lColor==0)) lColor = 1;
    if (i==0) lColor = 0;
    colors[i].pixel = IDdata->colut[i];
    colors[i].flags = DoRed|DoGreen|DoBlue;
    colors[i].red   = IDdata->red[lColor]; 
    colors[i].green = IDdata->green[lColor]; 
    colors[i].blue  = IDdata->blue[lColor]; 
    if (i==0)
      colors[i].red = colors[i].blue = colors[i].green = BlackPixel(dpy, 0);
  }  /* end of loop over color table */
  
  /* Don't bother if no colormap */
  if (!IDdata->cmap) {
    if (colors) {XtFree ((XtPointer) colors);} colors = NULL;
    return;}
  
  /* allocate colors */
  if (IDdata->trueColor) {
    RGB2Pixel (colors, IDdata);
    PaintImage (IDdata);       /* redraw */
  }  else {
    XStoreColors (dpy, IDdata->cmap, colors, IDdata->ncolors);
  }
  
  if (colors) {XtFree ((XtPointer) colors);}  colors = NULL;
  return;
}  /* end of SetColorTable */


/**
 * Set initial values in color table
 * IDdata->ncolors = number of colors (e.g. 128 or 256)
 * IDdata->red = default red color table (0 - 65535) 
 * IDdata->green = default green color table (0 - 65535) 
 * IDdata->blue = default blue color table (0 - 65535) 
 * \param IDdata   Image display
 */
void InitColorTable(ImageDisplay* IDdata)
{
  olong i;
  olong lColor;
  ofloat frac;
  Display     *dpy = XtDisplay (IDdata->canvas);
  
  for (i=0; i<IDdata->ncolors; i++) {
    frac = (ofloat)(i+1) / (ofloat)IDdata->ncolors;
    lColor = frac * 65536.0 + frac * 256.0;
    if (lColor > 65535) lColor = 65535;
    if (i==0) lColor = BlackPixel(dpy, 0);  /* reserve black for blanked */
    IDdata->red[i]   = lColor;
    IDdata->green[i] = lColor;
    IDdata->blue[i]  = lColor;
  }
} /* end InitColorTable */

/**
 * Set color scheme lifted from AIPS
 * IDdata->ncolors = number of colors (e.g. 128)
 * IDdata->red   =  red color table
 * IDdata->green =  green color table
 * IDdata->blue  =  blue color table
 * \param IDdata   Image display
 */
void SetColor1(ImageDisplay* IDdata)
{
  gshort i;
  olong    r, g, b, rt, gt, bt;
#if MAXCOLOR==256   /* 256 color table */
  unsigned char bc_tab[256]={0,0,            /*blue table  */
			     15,15,15,15,15,15,15,15,15,15,15,15,15,15,
			     15,15,15,15,15,15,15,15,15,15,15,15,15,15,   
			     72,72,72,72,72,72,72,72,72,72,72,72,72,72,
			     72,72,72,72,72,72,72,72,72,72,72,72,72,72,
			     127,127,127,127,127,127,127,127,127,127,127,127,127,127,
			     127,127,127,127,127,127,127,127,127,127,127,127,127,127,  
			     203,203,203,203,203,203,203,203,203,203,203,203,203,203,
			     203,203,203,203,203,203,203,203,203,203,203,203,203,203,
			     0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			     0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			     0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0, 
			     0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			     0,0,0,0,0,0,0,0,0,0,0,0,0,0, 0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  unsigned char gc_tab[256]={0,0,         /* green table  */
			     15,15,15,15,15,15,15,15,15,15,15,15,15,15,
			     15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15, 
			     0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			     0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			     0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			     0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			     76,76,76,76,76,76,76,76,76,76,76,76,76,76,
			     76,76,76,76,76,76,76,76,76,76,76,76,76,76,
			     59,59,59,59,59,59,59,59,59,59,59,59,59,59,
			     59,59,59,59,59,59,59,59,59,59,59,59,59,59,  
			     229,229,229,229,229,229,229,229,229,229,229,229,229,229,
			     229,229,229,229,229,229,229,229,229,229,229,229,229,229,
			     255,255,255,255,255,255,255,255,255,255,255,255,255,255,
			     255,255,255,255,255,255,255,255,255,255,255,255,255,255, 
			     89,89,89,89,89,89,89,89,89,89,89,89,89,89,
			     89,89,89,89,89,89,89,89,89,89,89,89,89,89, 
			     0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			     0,0,0,0,0,0,0,0,0,0,0,0,0,0};
  unsigned char rc_tab[256]={0,0,          /*red table  */
                             15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
			     15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,  
			     36,36,36,36,36,36,36,36,36,36,36,36,36,36,
			     36,36,36,36,36,36,36,36,36,36,36,36,36,36,
			     0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			     0,0,0,0,0,0,0,0,0,0,0,0,0,0, 
			     15,15,15,15,15,15,15,15,15,15,15,15,15,15,
			     15,15,15,15,15,15,15,15,15,15,15,15,15,15,
			     0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			     0,0,0,0,0,0,0,0,0,0,0,0,0,0, 
			     0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			     0,0,0,0,0,0,0,0,0,0,0,0,0,0,
			     255,255,255,255,255,255,255,255,255,255,255,255,255,255,
			     255,255,255,255,255,255,255,255,255,255,255,255,255,255,  
			     255,255,255,255,255,255,255,255,255,255,255,255,255,255,
			     255,255,255,255,255,255,255,255,255,255,255,255,255,255,
			     255,255,255,255,255,255,255,255,255,255,255,255,255,255,
			     255,255,255,255,255,255,255,255,255,255,255,255,255,255};
#elif MAXCOLOR==128   /* 128 color table */
  unsigned char bc_tab[128]={0,            /*blue table  */
			     15,15,15,15,15,15,15, 15,15,15,15,15,15,15,15,   
			     72,72,72,72,72,72,72, 72,72,72,72,72,72,72,
			     127,127,127,127,127,127,127, 127,127,127,127,127,127,127,  
			     203,203,203,203,203,203,203, 203,203,203,203,203,203,203,
			     0,0,0,0,0,0,0, 0,0,0,0,0,0,0,
			     0,0,0,0,0,0,0, 0,0,0,0,0,0,0,
			     0,0,0,0,0,0,0, 0,0,0,0,0,0,0, 
			     0,0,0,0,0,0,0, 0,0,0,0,0,0,0,
			     0,0,0,0,0,0,0, 0,0,0,0,0,0,0};
  unsigned char gc_tab[128]={0,         /* green table  */
			     15,15,15,15,15,15,15, 15,15,15,15,15,15,15,15, 
			     0,0,0,0,0,0,0,  0,0,0,0,0,0,0,
			     0,0,0,0,0,0,0,  0,0,0,0,0,0,0,
			     76,76,76,76,76,76,76, 76,76,76,76,76,76,76,
			     59,59,59,59,59,59,59,  59,59,59,59,59,59,59,  
			     229,229,229,229,229,229,229, 229,229,229,229,229,229,229,
			     255,255,255,255,255,255,255, 255,255,255,255,255,255,255, 
			     89,89,89,89,89,89,89, 89,89,89,89,89,89,89, 
			     0,0,0,0,0,0,0,  0,0,0,0,0,0,0};
  unsigned char rc_tab[128]={0,          /*red table  */
			     15,15,15,15,15,15,15,  15, 15,15,15,15,15,15,15,  
			     36,36,36,36,36,36,36, 36,36,36,36,36,36,36,
			     0,0,0,0,0,0,0, 0,0,0,0,0,0,0, 
			     15,15,15,15,15,15,15, 15,15,15,15,15,15,15,
			     0,0,0,0,0,0,0, 0,0,0,0,0,0,0, 
			     0,0,0,0,0,0,0, 0,0,0,0,0,0,0,
			     255,255,255,255,255,255,255, 255,255,255,255,255,255,255,  
			     255,255,255,255,255,255,255, 255,255,255,255,255,255,255,
			     255,255,255,255,255,255,255, 255,255,255,255,255,255,255};
#else /* 64 color table */
  unsigned char bc_tab[64]={0,            /*blue table  */
			    15,15,15,15,15,15,15,   72,72,72,72,72,72,72,
			    127,127,127,127,127,127,127,  203,203,203,203,203,203,203,
			    0,0,0,0,0,0,0,  0,0,0,0,0,0,0,
			    0,0,0,0,0,0,0,  0,0,0,0,0,0,0,
			    0,0,0,0,0,0,0};
  unsigned char gc_tab[64]={0,         /* green table  */
			    15,15,15,15,15,15,15,  0,0,0,0,0,0,0,
			    0,0,0,0,0,0,0,  76,76,76,76,76,76,76,
			    59,59,59,59,59,59,59,  229,229,229,229,229,229,229,
			    255,255,255,255,255,255,255,  89,89,89,89,89,89,89,
			    0,0,0,0,0,0,0};
  unsigned char rc_tab[64]={0,          /*red table  */
			    15,15,15,15,15,15,15,  36,36,36,36,36,36,36,
			    0,0,0,0,0,0,0,  15,15,15,15,15,15,15,
			    0,0,0,0,0,0,0,  0,0,0,0,0,0,0,
			    255,255,255,255,255,255,255,  255,255,255,255,255,255,255,
			    255,255,255,255,255,255,255};  /*255  */
#endif
  for (i=1; i<(short)IDdata->ncolors; i++)
    { rt = rc_tab[i]; gt = gc_tab[i]; bt = bc_tab[i];
    r = rt + rt * 256; g = gt + gt * 256; b = bt + bt * 256;
    if (r>65535) r = 65535;
    if (g>65535) g = 65535;
    if (b>65535) b = 65535;
    IDdata->red[i]   = r;
    IDdata->green[i] = g;
    IDdata->blue[i]  = b;}
}  /* End of SetColor1 */

/**
 * Set color scheme lifted from AIPS (PHLAME) 
 * IDdata->ncolors = number of colors (e.g. 128)
 * IDdata->red   =  red color table
 * IDdata->green =  green color table
 * IDdata->blue  =  blue color table
 * \param IDdata   Image display
 */
void SetColor2(ImageDisplay *IDdata)
{
  gshort i;
  olong    r, g, b, rt, gt, bt;
#if MAXCOLOR==256   /* 256 color table */
  unsigned char bc_tab[]=         /* blue table  */
    {0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,
     0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,
     0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,
     0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,
     0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,
     0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,
     0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,
     0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,
     0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,
     0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,
     0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,
     0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,
     0,0,  26,36,  53,59,  66,72,  79,84,  89,94,  99,103, 108,112,
     117,121, 115,122, 133,137, 140,143, 147,150, 154,157, 161,164, 167,171,
     174,177, 180,183, 186,188, 191,194, 197,199, 202,205, 208,210, 212,215,
     219,221, 224,226, 229,231, 233,235, 238,240, 242,245, 248,250, 252,254};
  unsigned char gc_tab[]=         /* green table  */
    {0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,
     0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,
     0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,
     0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,
     0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,
     0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,
     0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,  31,31,
     0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,  15,15,
     31,34,  38,42,  46,48,  51,54,  57,59,  62,65,  68,70,  72,75,
     77,79,  81,83,  86,87,  89,91,  93,95,  97,99, 101,102, 104,104,
     107,109, 110,112, 114,115, 117,118, 120,121, 123,125, 127,128, 130,130,
     133,134, 135,136, 138,139, 141,142, 144,145, 146,147, 149,150, 151,152,
     154,155, 156,157, 159,160, 162,163, 165,166, 167,169, 170,171, 172,173,
     175,177, 179,180, 182,183, 185,186, 189,190, 192,194, 195,196, 198,200,
     202,203, 205,207, 209,210, 212,214, 215,216, 218,220, 222,223, 225,226,
     228,229, 231,233, 234,235, 237,238, 240,242, 243,244, 246,247, 251,253};
  unsigned char rc_tab[]=          /*red table  */
    {0,0,   0,0,   0,0,   0,0,  0,0,   0,0,    0,0,   0,0,
     0,0,   0,0,   0,0,   0,0,   0,0,   0,0,   0,0,  15,15,
     29,32,  36,39,  43,46,  49,52, 55,57,  59,61,   64,66,  68,70,
     73,75,  77,79,  81,82,  84,86,  88,89,  91,93,  95,96,  98,100, 
     102,103, 105,107, 109,111, 112,113, 114,115, 117,119, 120,121, 123,124,
     126,127, 128,129, 131,132, 134,135, 137,138, 139,141, 142,143, 144,145,
     147,148, 149,150, 152,153, 154,155, 156,157, 158,160, 161,162, 163,164,
     165,166, 167,169, 170,171, 172,173, 174,175, 176,177, 178,179, 180,181,
     183,184, 185,186, 187,188, 189,190, 191,192, 193,194, 195,195, 196,197,
     198,199, 200,201, 202,203, 204,206, 207,208, 209,210, 211,211, 212,213,
     214,215, 216,217, 218,218, 219,220, 221,222, 223,224, 225,225, 226,227,
     228,229, 230,231, 232,232, 233,234, 235,235, 236,237, 238,239, 240,241,
     242,242, 243,244, 245,246, 247,247, 248,248, 249,249, 251,251, 253,253,
     255,255, 255,255, 255,255, 255,255, 255,255, 255,255, 255,255, 255,255,
     255,255, 255,255, 255,255, 255,255, 255,255, 255,255, 255,255, 255,255,
     255,255, 255,255, 255,255, 255,255, 255,255, 255,255, 255,255, 255,255};
#elif MAXCOLOR==128   /* 128 color table */
  unsigned char bc_tab[]=         /* green table  */
    {0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
     0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
     0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
     0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
     0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
     0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
     0,  26,  53,  66,  79,  89,  99, 108, 117, 115, 133, 140, 147, 154, 161, 167,
     174, 180, 186, 191, 197, 202, 208, 212, 219, 224, 229, 233, 238, 242, 248, 252};
  unsigned char gc_tab[]=         /* green table  */
    {0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
     0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
     0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
     0,   0,   0,   0,   0,   0,   0,  31,   0,   0,   0,   0,   0,   0,   0,  15,
     31,  38,  46,  51,  57,  62,  68,  72,  77,  81,  86,  89,  93,  97, 101, 104,
     107, 110, 114, 117, 120, 123, 127, 130, 133, 135, 138, 141, 144, 146, 149, 151,
     154, 156, 159, 162, 165, 167, 170, 172, 175, 179, 182, 185, 189, 192, 195, 198,
     202, 205, 209, 212, 215, 218, 222, 225, 228, 231, 234, 237, 240, 243, 246, 251};
  unsigned char rc_tab[]=          /*red table  */
    {0,   0,   0,   0,  0,   0,    0,   0,   0,   0,   0,   0,   0,   0,   0,  15,
     29,  36,  43,  49, 55,  59,   64,  68,  73,  77,  81,  84,  88,  91,  95,  98, 
     102, 105, 109, 112, 114, 117, 120, 123, 126, 128, 131, 134, 137, 139, 142, 144,
     147, 149, 152, 154, 156, 158, 161, 163, 165, 167, 170, 172, 174, 176, 178, 180,
     183, 185, 187, 189, 191, 193, 195, 196, 198, 200, 202, 204, 207, 209, 211, 212,
     214, 216, 218, 219, 221, 223, 225, 226, 228, 230, 232, 231, 235, 236, 238, 240,
     242, 244, 246, 247, 248, 249, 251, 253, 255, 255, 255, 255, 255, 255, 255, 255,
     255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255};
#else /* 64 color table */
  unsigned char bc_tab[64]=            /*blue table  */
    {0,   0,   0,   0,   0,   0,   0,   0,
     0,   0,   0,   0,   0,   0,   0,   0,
     0,   0,   0,   0,   0,   0,   0,   0,
     0,   0,   0,   0,   0,   0,   0,   0,
     0,   0,   0,   0,   0,   0,   0,   0,
     0,   0,   0,   0,   0,   0,   0,   0,
     0,  53,  79,  99, 117, 133, 147, 161,
     174, 186, 197, 208, 219, 229, 238, 248};
  unsigned char gc_tab[64]=         /* green table  */
    {0,   0,   0,   0,   0,   0,   0,   0,
     0,   0,   0,   0,   0,   0,   0,   0,
     0,   0,   0,   0,   0,   0,   0,   0,
     0,   0,   0,   0,   0,   0,   0,  31,
     46,  57,  68,  77,  86,  93, 101, 107,
     114, 120, 127, 133, 138, 144, 149, 154,
     159, 165, 170, 175, 182, 189, 195, 202,
     209, 215, 222, 228, 234, 240, 246, 251};
  unsigned char rc_tab[64]=          /*red table  */
    {0,   0,   0,   0,   0,   0,   0,  29,
     43,  55,  64,  73,  81,  88,  95, 102,
     109, 114, 120, 126, 131, 137, 142, 147,
     152, 156, 161, 165, 170, 174, 178, 183,
     187, 191, 195, 198, 202, 207, 211, 214,
     218, 221, 225, 228, 232, 235, 238, 243,
     246, 248, 251, 255, 255, 255, 255, 255,
     255, 255, 255, 255, 255, 255, 255, 255};
#endif
  
  for (i=1; i<(short)IDdata->ncolors; i++) 
    { rt = rc_tab[i]; gt = gc_tab[i]; bt = bc_tab[i];
    r = rt + rt * 256; g = gt + gt * 256; b = bt + bt * 256;
    if (r>65535) r = 65535;
    if (g>65535) g = 65535;
    if (b>65535) b = 65535;
    IDdata->red[i] =   r;
    IDdata->green[i] = g;
    IDdata->blue[i] =  b;}
}  /* End of SetColor2 */

/**
 * Reverse color table
 * IDdata->ncolors = number of colors (e.g. 128)
 * IDdata->red =  red color table
 * IDdata->green =  green color table
 * IDdata->blue =  blue color table
 * \param IDdata   Image display
 */
void RevColor(ImageDisplay *IDdata)
{
  gshort i, j;
  unsigned short cTemp;
  for (i=1; i<IDdata->ncolors/2; i++)
    {j = IDdata->ncolors - i;
    cTemp = IDdata->red[i];
    IDdata->red[i] =   IDdata->red[j];
    IDdata->red[j] =   cTemp;
    cTemp = IDdata->green[i];
    IDdata->green[i] =   IDdata->green[j];
    IDdata->green[j] =   cTemp;
    cTemp = IDdata->blue[i];
    IDdata->blue[i] =   IDdata->blue[j];
    IDdata->blue[j] =   cTemp;
    }
}  /* End of RevColor */

/**
 * Update the color map
 * IDdata->value[0] = Brightness control (0 - 255) 
 * IDdata->value[1] = Contrast control (0- 255) 
 * IDdata->ncolors  = number of colors (e.g. 128) 
 * IDdata->red      = default red color table 
 * IDdata->green    = default green color table 
 * IDdata->blue     = default blue color table 
 * The new color table is the entries from the default table with  
 * 2(IDdata->value[0]-128) added to each and  then multiplied by  
 * 1 + 5*((IDdata->value[1]-128)/128); 
 * \param IDdata   Image display
 */
void UpdateColor (ImageDisplay *IDdata)
{
  SetColorTable (IDdata);
}  /* end of UpdateColor */

/**
 * Reallocate colors and get pixel values
 * Convert RGB values to pixel values IDdata->coltab 
 * Used IDdata members:
 * \li colors = previously allocated colors 
 * \li coltab = pixel values for colors 
 *
 * \param colors an array of color values to be used 
 *   in computing IDdata->coltab
 * \param IDdata   Image display
 */
void RGB2Pixel (XColor *colors, ImageDisplay *IDdata)
{
  olong i;
  
  /* Don't bother if no colormap */
  if (!IDdata->cmap) return;
  
  for (i=0; i<IDdata->ncolors; i++) {
    /* deallocate old colors (don't need for static display)*/
    /*    if (IDdata->colors[i].pixel>=0)
	  XFreeColors (XtDisplay (IDdata->canvas), IDdata->cmap, 
	  &IDdata->colors[i].pixel, 
	  1, IDdata->depth);*/
    /* allocate new color */
    XAllocColor (XtDisplay (IDdata->canvas), IDdata->cmap, &colors[i]);
    IDdata->colors[i].pixel = colors[i].pixel;
    IDdata->colors[i].red   = colors[i].red;
    IDdata->colors[i].green = colors[i].green;
    IDdata->colors[i].blue  = colors[i].blue;
    IDdata->colors[i].flags = colors[i].flags;
    IDdata->coltab[i]       = colors[i].pixel;
  }
}  /* end of RGB2Pixel */

void InitGraphColor(ImageDisplay* IDdata)
/* Setup graphics colors */
/* IDdata->graph   = number of graphics colors (e.g. 6) */
/* IDdata->red     = red color table   (0 - 65535) */
/* IDdata->green   = green color table (0 - 65535) */
/* IDdata->blue    = blue color table  (0 - 65535) */
{
   olong i;
   XColor *colors;
   Display     *dpy = XtDisplay (IDdata->canvas);
   /* Graphics colors for 6 colors */
   int rg[] = {1, 1, 1, 0, 0, 0};
   int gg[] = {0, 1, 0, 1, 0, 1};
   int bg[] = {1, 0, 0, 1, 1, 0};

   /* allocate color array */
   colors = (XColor*) XtMalloc (sizeof (XColor) * (IDdata->ngraph));

   /* graphics colors - 6 */
   for (i=0; i<IDdata->ngraph; i++) {
     IDdata->colut[i+IDdata->ncolors]  = IDdata->colut[i+IDdata->ncolors-1]+1;
     IDdata->coltab[i+IDdata->ncolors] = IDdata->colut[i+IDdata->ncolors];
     colors[i].pixel = IDdata->colut[i+IDdata->ncolors];
     colors[i].flags = DoRed|DoGreen|DoBlue;
     colors[i].red   = rg[i]*65535;
     colors[i].green = gg[i]*65535;
     colors[i].blue  = bg[i]*65535;
   }

   /* allocate colors */
   if (IDdata->trueColor) {   /* True color */
     for (i=0; i<IDdata->ngraph; i++) {
       XAllocColor (XtDisplay (IDdata->canvas), IDdata->cmap, &colors[i]);
       IDdata->colors[i+IDdata->ncolors].pixel = colors[i].pixel;
       IDdata->colors[i+IDdata->ncolors].red   = colors[i].red;
       IDdata->colors[i+IDdata->ncolors].green = colors[i].green;
       IDdata->colors[i+IDdata->ncolors].blue  = colors[i].blue;
       IDdata->colors[i+IDdata->ncolors].flags = colors[i].flags;
       IDdata->coltab[i+IDdata->ncolors]       = colors[i].pixel;
     }
   }  else {   /* Pseudo color */
     XStoreColors (dpy, IDdata->cmap, colors, IDdata->ngraph);
   }
   
   /* Cleanup */
   if (colors) {XtFree ((XtPointer) colors);}  colors = NULL;
   
 } /* end InitGraphColor */
