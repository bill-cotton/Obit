/* $Id$ */
/*   Interactively draw boxes */
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
#include <stdio.h> 
#include <string.h> 
#include <math.h>
#include <Xm/Xm.h> 
#include "drawbox.h"
#include "graph.h"
#include "messagebox.h"
#include "requestbox.h"
#define min(a,b) (((a) < (b)) ? (a) : (b))
#define max(a,b) (((a) > (b)) ? (a) : (b))

/* file global structure for things to talk to each other */
typedef struct {
  ImageDisplay *IDdata;          /* Image display */
  ObitDConCleanWindow* myWindow; /* Clean window being edited */
  gboolean modify;               /* True if any window modified */
  gboolean editActive;           /* True currently editing a clean window */
  olong   boxType;                /* Type of box being drawn */
  olong   which;                  /* which component of box to modify */
  olong  box[4];                 /* box specification currently being edited */
  olong  field;                  /* which field in myWindow */
  olong  Id;                     /* which window Id? */
  short  xpos, ypos;             /* location of selected pixel */
} DrawBoxStuff;
DrawBoxStuff DrawData;

/**
 *  \file drawbow.c
 * Interactively sets boxes on an image
 */
/*---------------Private function prototypes----------------*/
static void DrawAll (ImageData *ImgData);
static void DrawAnyBox (ImageData *ImgData, olong *box, olong gcolor);
static void DrawRect (ImageData *ImgData, olong *box, olong gcolor);
static void DrawCirc (ImageData *ImgData, olong *box, olong gcolor);
static void TellWindow (void);
void DrawBoxButEH (Widget w, XtPointer data, XEvent *event);
void DrawBoxKeyEH (Widget w, XtPointer data, XEvent *event);


/*------------------- Public functions -------------------*/
/**
 * Initialize structures for drawing and interactively editing boxes.
 * \param IDdata      ImageDisplay
 */
void DrawBoxInit (ImageDisplay *IDdata)
{

  olong i;
  /* Save Image display */
  DrawData.IDdata = IDdata;

  /* Other initialization */
  DrawData.myWindow   = NULL;
  DrawData.modify     = FALSE;
  DrawData.editActive = FALSE;
  DrawData.boxType    = 0;
  DrawData.field      = 0;
  DrawData.Id         = -1;
  for (i=0; i<4; i++) DrawData.box[i] = 0;

} /* end DrawBoxInit */

/**
 * Draw Clean window on image
 * Box types supported:
 * \li Rectangle, box[0]>0, blc={box[0], box[1]}, trc={box[2], box[3]}
 * \li Circle, box[0]=-10, center={box[2], box[3]}, radius = box[1]
 *
 * \param IDdata      ImageDisplay
 * \param window      Clean window object
 *                    input: initial window
 *                    output: user selected window
 *                    Will be Unrefed when done.
 * \param field       field in window corresponding to image
 * \return returns 0 if OK else failed.
 */
olong DrawBox (ImageDisplay *IDdata, ObitDConCleanWindow* window, 
	      olong field)
{
  olong boxType=0, i;
  char errMsg[121];

  /* Range check */
  if ((IDdata->gcolor<0) || (IDdata->gcolor>MAXGRAPH)) {
    g_snprintf (errMsg, 120, "GraphDrawLine: color %d out of range [0,%d]", 
	     IDdata->gcolor, MAXGRAPH);
    MessageShow(errMsg);
    return -1;
  }

   /* Check image validity */
  if (!image[CurImag].valid) {
    g_snprintf (errMsg, 120, "No image loaded");
    MessageShow(errMsg);
    return -1;
 }
    
 /* Check that window and current image are compatible */
  if ((window->naxis[field-1][0]!=image[CurImag].myDesc->inaxes[0]) || 
      (window->naxis[field-1][1]!=image[CurImag].myDesc->inaxes[1])) {
    g_snprintf (errMsg, 120, "window/image sizes incompatible [ %d, %d] [ %d, %d]", 
	     window->naxis[field-1][0], window->naxis[field-1][1],
	     image[CurImag].myDesc->inaxes[0],
	     image[CurImag].myDesc->inaxes[1]);
    MessageShow(errMsg);
    return -1;
 }
    
  /* Save info where the event handler can get to it */
  for (i=0; i<4; i++) DrawData.box[i] = 0;
  DrawData.boxType  = boxType;
  DrawData.which    = 1;
  DrawData.IDdata   = IDdata;
  DrawData.myWindow = ObitDConCleanWindowUnref(DrawData.myWindow);
  DrawData.myWindow = window;
  DrawData.modify     = FALSE;
  DrawData.editActive = FALSE;
  DrawData.field      = 1;
  DrawData.Id         = -1;

  /* Update image pointer to window */
  image[CurImag].edtWindow = ObitDConCleanWindowUnref(image[CurImag].edtWindow);
  image[CurImag].edtWindow = ObitDConCleanWindowRef(window);

  /* Draw initial */
  DrawAll(&image[CurImag]);

  PaintImage(IDdata); /* redraw image */

  /* Request Box */
  EditRequestBox ();

  return 0;
} /* end DrawBox */

/**
 * Interactively Edit CLEAN boxes
 * Uses parameters saved by DrawBox
 * Box types supported:
 * \li Rectangle, box[0]>0, blc={box[0], box[1]}, trc={box[2], box[3]}
 * \li Circle, box[0]=-10, center={box[2], box[3]}, radius = box[1]
 *
 * Keys for user editing:
 * \li [aA]: switch part of window being edited; blc/trc for rectangle,
 *           center/radius for circle
 * \li [bB]: Delete current box
 * \li [cC]: edit another box (goes into box search mode)
 * \li [dD]: finish editing
 * \li [eE]: new rectangle
 * \li [fF]: new circle
 * \li [fF]: new circle
 * \li [gG]: toggle box/unbox
 *
 * \return returns 0 if OK else failed.
 */
olong EditBox (void)
{
  olong i;
  char errMsg[121];

  /* Check image validity */
  if (!image[CurImag].valid) {
    g_snprintf (errMsg, 120, "No image loaded");
    MessageShow(errMsg);
    return -1;
 }
    
  /* Instructions to user */
  MessageShow(" Begin editing boxes, click mode to specify");
  MessageShow(" [Aa]: Toggle blc/trc center/radius");
  MessageShow(" [Bb]: Delete current box");
  MessageShow(" [Cc]: Edit another box (click to specify)");
  MessageShow(" [Dd]: Exit editing");
  MessageShow(" [Ee]: New rectangle ");
  MessageShow(" [Ff]: New Circle ");
  MessageShow(" [Gg]: Toggle box/unbox");

  /* Save info where the event handler can get to it */
  for (i=0; i<4; i++) DrawData.box[i] = 0;
  DrawData.which      = 1;
  DrawData.modify     = FALSE;
  DrawData.editActive = TRUE;
  DrawData.Id         = -1;

  /* keystrokes in canvas */
  XtInsertEventHandler (DrawData.IDdata->canvas, 
			KeyPressMask, FALSE, 
			(XtEventHandler)DrawBoxKeyEH, 
			(XtPointer)&DrawData, XtListHead);
  /* Trap mouse clicks in canvas */
  XtInsertEventHandler (DrawData.IDdata->canvas, 
			ButtonPressMask|ButtonMotionMask, FALSE, 
			(XtEventHandler)DrawBoxButEH, 
			(XtPointer)&DrawData, XtListHead);
 
  /* Draw initial */
  DrawAll(&image[CurImag]);

  PaintImage(DrawData.IDdata); /* redraw image */

  g_snprintf (errMsg, 120, "editing box  %d", DrawData.Id);
  if (DrawData.Id>0) MessageShow(errMsg);

  return 0;
} /* end EditBox */

/**
 * Delete all windows in current structure
 * Uses parameters saved by DrawBox
 *
 * \return returns 0 if OK else failed.
 */
olong ClearBox (void)
{
  olong i;
  char errMsg[121];

  /* Check image validity */
  if (!image[CurImag].valid) {
    g_snprintf (errMsg, 120, "No image loaded");
    MessageShow(errMsg);
    return -1;
 }
    
  /* Save info where the event handler can get to it */
  for (i=0; i<4; i++) DrawData.box[i] = 0;
  DrawData.which      = 1;
  DrawData.Id         = -1;

  /* clear */
  ObitDConCleanWindowDel (DrawData.myWindow, 1, -1, err);
  if (err->error) ObitErrLog(err);
  DrawData.modify = TRUE;

  /* redraw */
  DrawAll(&image[CurImag]);

  PaintImage(DrawData.IDdata); /* redisplay image */

  return 0;
} /* end ClearBox */

/*------------------- Private functions -------------------*/
/**
 * Draw All boxes and leave the last on in editing mode
 * \param ImgData     ImageData structure pointer
 */
static void DrawAll (ImageData *ImgData)
{
  olong i, iD, nId, color=0;
  olong *win;
  ObitDConCleanWindowType type;

  /* Clear graphics */
  GraphClear (ImgData);

  /* Draw outer box in another color */
  iD=-1;
  if (ObitDConCleanWindowInfo(DrawData.myWindow, DrawData.field, iD, &type, &win, err)) {
    /* Which box type */
    if (type==OBIT_DConCleanWindow_rectangle)  DrawData.boxType = 1; /* rectangle */
    if (type==OBIT_DConCleanWindow_round)      DrawData.boxType = 2; /* circle */
    if (type==OBIT_DConCleanWindow_unrectangle) DrawData.boxType = 3; /* unrectangle */
    if (type==OBIT_DConCleanWindow_unround)     DrawData.boxType = 4; /* uncircle */
    
    /* convert to 0 rel */
    if ((DrawData.boxType==1) || (DrawData.boxType==3)) {
      for (i=0; i<4; i++) DrawData.box[i] = (olong)win[i]-1;
    }
    if ((DrawData.boxType==2) || (DrawData.boxType==4)) {
      DrawData.box[0]=-1; DrawData.box[1]=(olong)win[0]; 
      DrawData.box[2]=(olong)win[1]-1; DrawData.box[3]=(olong)win[2]-1;
    }
    
    /* Draw  */
    color = (DrawData.IDdata->gcolor+2) % 8;
    DrawAnyBox (&image[CurImag], DrawData.box, color);
    
  } /* end if iD present */
  /* Any errors? */
  if (err->error) ObitErrLog(err);

  /* How many potential inner windows? */
  nId = DrawData.myWindow->maxId[DrawData.field-1];
  /* Extract from Obit Object and draw */
  for (iD=1; iD<=nId; iD++) {
    if (ObitDConCleanWindowInfo(DrawData.myWindow, DrawData.field, iD, &type, &win, err)) {
      /* Which box type */
      if (type==OBIT_DConCleanWindow_rectangle)  DrawData.boxType = 1; /* rectangle */
      if (type==OBIT_DConCleanWindow_round)      DrawData.boxType = 2; /* circle */
      if (type==OBIT_DConCleanWindow_unrectangle) DrawData.boxType = 3; /* unrectangle */
      if (type==OBIT_DConCleanWindow_unround)     DrawData.boxType = 4; /* uncircle */
      
      /* convert to 0 rel */
      if (DrawData.boxType==1) {
	color = DrawData.IDdata->gcolor;
	for (i=0; i<4; i++) DrawData.box[i] = (olong)win[i]-1;
      }
      if (DrawData.boxType==2) {
	color = DrawData.IDdata->gcolor;
	DrawData.box[0]=-1; DrawData.box[1]=(olong)win[0]; 
	DrawData.box[2]=(olong)win[1]-1; DrawData.box[3]=(olong)win[2]-1;
      }
      if (DrawData.boxType==3) {
	color = (DrawData.IDdata->gcolor+1) % 8;
	for (i=0; i<4; i++) DrawData.box[i] = (olong)win[i]-1;
      }
      if (DrawData.boxType==4) {
	color = (DrawData.IDdata->gcolor+1) % 8;
	DrawData.box[0]=-1; DrawData.box[1]=(olong)win[0]; 
	DrawData.box[2]=(olong)win[1]-1; DrawData.box[3]=(olong)win[2]-1;
      }
      /* DEBUG
      fprintf (stderr,"Get Box %d %d  %d %d  %d\n",iD,win[0],win[1], win[2],win[3]); */
      
      /* Draw  */
      DrawAnyBox (&image[CurImag], DrawData.box, color);

      /* In case this one is last */
      DrawData.Id = iD;
     
    } /* end if iD present */
  } /* End loop over fields */
  /* Any errors? */
  if (err->error) ObitErrLog(err);

}  /*  end  DrawAll */

/**
 * Draw rectangle in graphics overlay in ImgData
 * \param ImgData     ImageData structure pointer
 * \param box         blc={box[0], box[1]}, trc={box[2], box[3]}
 * \param gcolor      color to write, 0 = erase
 */
static void DrawRect (ImageData *ImgData, olong *box, olong gcolor)
{
  olong blc[2], trc[2];

  /* Draw 4 lines around rectangle */
  blc[0] = box[0]; blc[1] = box[1];
  trc[0] = box[0]; trc[1] = box[3];
  GraphDrawLine (ImgData, blc, trc, gcolor);

  blc[0] = (olong)box[0]; blc[1] = (olong)box[3];
  trc[0] = (olong)box[2]; trc[1] = (olong)box[3];
  GraphDrawLine (ImgData, blc, trc, gcolor);

  blc[0] = (olong)box[2]; blc[1] = (olong)box[3];
  trc[0] = (olong)box[2]; trc[1] = (olong)box[1];
  GraphDrawLine (ImgData, blc, trc, gcolor);

  blc[0] = (olong)box[2]; blc[1] = (olong)box[1];
  trc[0] = (olong)box[0]; trc[1] = (olong)box[1];
  GraphDrawLine (ImgData, blc, trc, gcolor);
}
 /*  end  DrawRect */

/**
 * Draw circle in graphics overlay in ImgData
 * \param ImgData     ImageData structure pointer
 * \param box         center={box[2], box[3]}, radius = box[1]
 * \param gcolor      color to write, 0 = erase
 */
static void DrawCirc (ImageData *ImgData, olong *box, olong gcolor)
{
  olong center[2], radius;

  center[0] = (olong)box[2]; center[1] = (olong)box[3];
  radius =  (olong)box[1]; 
  GraphDrawCircle (ImgData, center, radius, gcolor);

} /*  end DrawCirc */

/**
 * Draw box in graphics overlay in ImgData
 * \param ImgData   ImageData structure pointer
 * \param box       box parameters
 * \param gcolor    color to write, 0 = erase
 */
static void DrawAnyBox (ImageData *ImgData, olong *box, olong gcolor)
{
  if (DrawData.boxType==1) DrawRect (&image[CurImag], DrawData.box, gcolor);
  if (DrawData.boxType==2) DrawCirc (&image[CurImag], DrawData.box, gcolor);
  if (DrawData.boxType==3) DrawRect (&image[CurImag], DrawData.box, gcolor);
  if (DrawData.boxType==4) DrawCirc (&image[CurImag], DrawData.box, gcolor);
}  /*  end  DrawAnyBox */

/**
 * Print current contents of myWindow
 */
static void TellWindow (void)
{
  olong iD, nId;
  olong *win;
  ObitDConCleanWindowType type;
  char errMsg[121];

  /* How many potential windows? */
  nId = DrawData.myWindow->maxId[DrawData.field-1];
  /* Extract from Obit Object and draw */
  for (iD=1; iD<=nId; iD++) {
    if (ObitDConCleanWindowInfo(DrawData.myWindow, DrawData.field, iD, &type, &win, err)) {
      /* Which box type */
      if (type==OBIT_DConCleanWindow_rectangle) 
	g_snprintf (errMsg, 120, "Box %d rectangle, blc = [ %d, %d], trc = [ %d, %d]", 
		 iD, win[0],win[1], win[2], win[3]);
      if (type==OBIT_DConCleanWindow_round)
	g_snprintf (errMsg, 120, "Box %d round, radius =  %d, center = [ %d, %d]", 
		 iD, win[0], win[1], win[2]);
      if (type==OBIT_DConCleanWindow_unrectangle) 
	g_snprintf (errMsg, 120, "Unbox %d rectangle, blc = [ %d, %d], trc = [ %d, %d]", 
		 iD, win[0],win[1], win[2], win[3]);
      if (type==OBIT_DConCleanWindow_unround)
	g_snprintf (errMsg, 120, "Unbox %d round, radius =  %d, center = [ %d, %d]", 
		 iD, win[0], win[1], win[2]);
      MessageShow(errMsg);
   }
  } /* end loop over windows */
  /* Any errors? */
  if (err->error) ObitErrLog(err);
} /* end TellWindow */

/**
 * Event handler for button press in image display window
 * \param w      canvas widget activated
 * \param data   event data
 * \param event  event
 */
void DrawBoxButEH (Widget w, XtPointer data, XEvent *event)
{
  olong         ix, iy, nx, ny, iD, i, which, color;
  olong        *win, pixel[2], toler, twin[4];
  ofloat       fx, fy, z, dist;
  ObitDConCleanWindowType type;
  ImageDisplay *IDdata = DrawData.IDdata;
  char errMsg[121];
  
  if (!IDdata) return;
  if (!image[CurImag].valid) return;
  if (!DrawData.editActive) return;
  
  /* determine image pixel accounting for zoom */
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

  /* Pixel positing (0-rel) in image */
  ix = (olong)(fx+0.5); iy = (olong)(fy+0.5);
  nx = image[CurImag].myDesc->inaxes[0];
  ny = image[CurImag].myDesc->inaxes[1];

  /* Editing or searching? */

  if (DrawData.Id>0) { /* editing */
    /* DEBUG
    fprintf (stderr,"DEBUG type %d box %d  %d %d  %d\n", DrawData.boxType,
	     DrawData.box[0],DrawData.box[1],DrawData.box[2],DrawData.box[3]); */
    /* erase previous */
    if (DrawData.boxType==1)  
      DrawRect (&image[CurImag], DrawData.box, 0);
    else if (DrawData.boxType==2) 
      DrawCirc (&image[CurImag], DrawData.box, 0);
    else if (DrawData.boxType==3)  
      DrawRect (&image[CurImag], DrawData.box, 0);
    else if (DrawData.boxType==4) 
      DrawCirc (&image[CurImag], DrawData.box, 0);
    
    /* Change blc for rectangle or center for circle */
    if ((DrawData.boxType==1) || (DrawData.boxType==3)) {            /* Rectangle */
      if (DrawData.which==1) { 	/* DrawData.which  = 2;*/
	DrawData.box[0] = ix;
	DrawData.box[1] = iy;
	DrawData.modify = TRUE;
	twin[0] = DrawData.box[0]+1; twin[1] = DrawData.box[1]+1;
	twin[2] = DrawData.box[2]+1; twin[3] = DrawData.box[3]+1;
	if (DrawData.boxType==1) 
	  ObitDConCleanWindowUpdate (DrawData.myWindow, DrawData.field,
				     DrawData.Id, OBIT_DConCleanWindow_rectangle,
				     twin, err);
	else
	  ObitDConCleanWindowUpdate (DrawData.myWindow, DrawData.field,
				     DrawData.Id, OBIT_DConCleanWindow_unrectangle,
				     twin, err);
      } else { 	/* DrawData.which  = 1;*/
	DrawData.box[2] = ix;
	DrawData.box[3] = iy;
	DrawData.modify = TRUE;
	twin[0] = DrawData.box[0]+1; twin[1] = DrawData.box[1]+1;
	twin[2] = DrawData.box[2]+1; twin[3] = DrawData.box[3]+1;
	if (DrawData.boxType==1) 
	  ObitDConCleanWindowUpdate (DrawData.myWindow, DrawData.field,
				     DrawData.Id, OBIT_DConCleanWindow_rectangle,
				     twin, err);
	else
	  ObitDConCleanWindowUpdate (DrawData.myWindow, DrawData.field,
				     DrawData.Id, OBIT_DConCleanWindow_unrectangle,
				     twin, err);
      }
    } else if ((DrawData.boxType==2) || (DrawData.boxType==4)) {  /* Circle */
      if (DrawData.which==1) {	/* DrawData.which  = 1;*/
	DrawData.box[2] = ix;
	DrawData.box[3] = iy;
	DrawData.modify = TRUE;
	twin[0] = DrawData.box[1];
	twin[1] = DrawData.box[2]+1; twin[2] = DrawData.box[3]+1;
	if (DrawData.boxType==2)
	  ObitDConCleanWindowUpdate (DrawData.myWindow, DrawData.field,
				     DrawData.Id, OBIT_DConCleanWindow_round,
				     twin, err);
	else
 	  ObitDConCleanWindowUpdate (DrawData.myWindow, DrawData.field,
				     DrawData.Id, OBIT_DConCleanWindow_unround,
				     twin, err);
     /* DEBUG 
      fprintf (stderr,"Set Box %d  %d %d  %d\n",DrawData.Id,twin[0],twin[1], twin[2]);*/

      } else { 	/*DrawData.which  = 2;*/
	dist = (fx -  DrawData.box[2]) * (fx -  DrawData.box[2]) +
	  (fy -  DrawData.box[3]) * (fy -  DrawData.box[3]);
	dist = sqrt (max (0.01, dist));
	DrawData.box[1] = (olong)(dist + 0.5);
	DrawData.modify = TRUE;
	twin[0] = DrawData.box[1];
	twin[1] = DrawData.box[2]+1; twin[2] = DrawData.box[3]+1;
	if (DrawData.boxType==2)
	  ObitDConCleanWindowUpdate (DrawData.myWindow, DrawData.field,
				     DrawData.Id, OBIT_DConCleanWindow_round,
				     twin, err);
	else
	  ObitDConCleanWindowUpdate (DrawData.myWindow, DrawData.field,
				     DrawData.Id, OBIT_DConCleanWindow_unround,
				     twin, err);
      }
    } /* end change box */
    /* Any errors? */
    if (err->error) ObitErrLog(err);
    
    /* Redraw */
    color = DrawData.IDdata->gcolor;
    /* different color for unboxes */
    if ((DrawData.boxType==3) || (DrawData.boxType==4)) 
      color = (DrawData.IDdata->gcolor+1)%8;
    DrawAnyBox (&image[CurImag], DrawData.box, color);
    PaintImage(IDdata); /* redraw image */
    return;
    /* end editing */

    /* Searching for new box to edit? */
  } else if (DrawData.Id<0) {
    pixel[0] = ix; pixel[1] = iy; /* Pixel 1-rel */
    toler = 4; /* how close is needed? */
    iD = ObitDConCleanWindowSearch(DrawData.myWindow, DrawData.field, pixel, 
				   toler, &which, err);
    if (iD>0) { /* find one? */
      if (ObitDConCleanWindowInfo(DrawData.myWindow, DrawData.field, iD, 
				  &type, &win, err)) {
	DrawData.Id    = iD;
	DrawData.which = which;
	/* Which box type */
	if (type==OBIT_DConCleanWindow_rectangle)  DrawData.boxType = 1; /* rectangle */
	if (type==OBIT_DConCleanWindow_round)      DrawData.boxType = 2; /* circle */
	if (type==OBIT_DConCleanWindow_unrectangle)  DrawData.boxType = 3; /* unrectangle */
	if (type==OBIT_DConCleanWindow_unround)      DrawData.boxType = 4; /* uncircle */
	
	/* box 0-rel, win 1-rel */
	if ((DrawData.boxType==1) || (DrawData.boxType==3))
	  for (i=0; i<4; i++) DrawData.box[i] = (olong)(win[i]-1);
	if ((DrawData.boxType==2) || (DrawData.boxType==4)) {
	  DrawData.box[0]=-1; DrawData.box[1]=(olong)win[0]; 
	  DrawData.box[2]=(olong)(win[1]-1); DrawData.box[3]=(olong)(win[2]-1);
	}

	/* tell which */
	g_snprintf (errMsg, 120, "editing box %d", iD);
	MessageShow(errMsg);
 
	/* DEBUG 
	fprintf (stderr,"Get Box %d  %d %d  %d\n",DrawData.Id,win[0],win[1], win[2]);*/

      } /* end if found window */
      /* Any errors? */
      if (err->error) ObitErrLog(err);
      return;
    }
  } /* end of search for new box to edit */
} /* end DrawBoxButEH */

/**
 * Event handler for keystrokes in image display window
 * \li A or a => switch part of box being adjusted
 * \li B or b => Delete current box
 * \li C or c => edit another box (goes into box search mode)
 * \li D or d => Exit box editing mode 
 * \li E or e => New rectangle
 * \li F or f => New circle
 *
 * \param w      canvas widget activated
 * \param data   event data
 * \param event  event
 */
void DrawBoxKeyEH (Widget w, XtPointer data, XEvent *event)
{
  gchar        buffer[11];
  ofloat       fx, fy, z;
  olong        twin[4];
  olong         ix, iy, color;
  char         errMsg[121];
  ImageDisplay *IDdata = DrawData.IDdata;
  
  if (!IDdata) return;
  if (!image[CurImag].valid) return;
  if (!DrawData.editActive) return;


  /* Get which key hit */
  XLookupString ((XKeyEvent*)event, buffer, 10, NULL, NULL);

  /* determine image pixel accounting for zoom */
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

  /* Pixel 0-rel positing in image */
  ix = (olong)(fx+0.5); iy = (olong)(fy+0.5);

  /* If A or a switch which parameter being modified */
  if ((buffer[0]=='a') || (buffer[0]=='A')) {
    DrawData.which = 3 -  DrawData.which;
    return;
  }

  /* if B or b delete box */
  if ((buffer[0]=='b') || (buffer[0]=='B')) {
    /* erase */
    DrawAnyBox (&image[CurImag], DrawData.box, 0);
    PaintImage(IDdata); /* redraw image */
    /* remove from myWindow */
    ObitDConCleanWindowDel (DrawData.myWindow, DrawData.field, DrawData.Id, err);
    DrawData.Id = -1;
    MessageShow(" Box deleted");
    /* Any errors? */
    if (err->error) ObitErrLog(err);
    return;
  }
  
  /* if C or c edit another box */
  if ((buffer[0]=='c') || (buffer[0]=='C')) {
    if (DrawData.Id>0) {
      DrawData.Id = -1;
    }
    MessageShow(" Select another box");
    return;
  }
  
  /* If D or d Exit box editing mode */
  if ((buffer[0]=='d') || (buffer[0]=='D')) {
    /* Handler for  mouse clicks in canvas */
    XtRemoveEventHandler (IDdata->canvas, ButtonPressMask, FALSE, 
			  (XtEventHandler)DrawBoxButEH, 
			  (XtPointer)&DrawData);
    /* keystrokes in canvas */
    XtRemoveEventHandler (IDdata->canvas, KeyPressMask, FALSE, 
			  (XtEventHandler)DrawBoxKeyEH, 
			  (XtPointer)&DrawData);
   
    MessageShow(" End editing boxes");
    DrawData.editActive = FALSE;
    TellWindow();
    /* Check for request */
    EditRequestBox ();
    return;
  }

  /* if E or e new rectangle */
  if ((buffer[0]=='e') || (buffer[0]=='E')) {
    /* Put it where the button clicked */
    DrawData.box[0] = ix-4; DrawData.box[1] = iy-5; 
    DrawData.box[2] = ix+6; DrawData.box[3] = iy+5; 
    twin[0] = DrawData.box[0]+1; twin[1] = DrawData.box[1]+1;
    twin[2] = DrawData.box[2]+1; twin[3] = DrawData.box[3]+1;
    DrawData.boxType = 1;
    DrawData.which   = 1;
    DrawData.Id = ObitDConCleanWindowAdd(DrawData.myWindow, DrawData.field, 
					 OBIT_DConCleanWindow_rectangle,
					 twin, err);
    DrawAnyBox (&image[CurImag], DrawData.box, DrawData.IDdata->gcolor);
    PaintImage(IDdata); /* redraw image */
    /* Any errors? */
    if (err->error) ObitErrLog(err);
    return;
  } /* end Ee */
  
  /* if F or f new circle */
  if ((buffer[0]=='f') || (buffer[0]=='F')) {
    /* Put it where the button clicked */
    DrawData.box[0] = -1;   DrawData.box[1] = 5;
    DrawData.box[2] = ix; DrawData.box[3] = iy; 
    twin[0] = DrawData.box[1];
    twin[1] = DrawData.box[2]+1; twin[2] = DrawData.box[3]+1;
    DrawData.boxType = 2;
    DrawData.which   = 1;
    DrawData.Id = ObitDConCleanWindowAdd(DrawData.myWindow, DrawData.field, 
					 OBIT_DConCleanWindow_round,
					 twin, err);
    DrawAnyBox (&image[CurImag], DrawData.box, DrawData.IDdata->gcolor);
    PaintImage(IDdata); /* redraw image */
    /* Any errors? */
    if (err->error) ObitErrLog(err);
    return;
  } /* end Ff */

  /* if G or g toggle box/unbox */
  if ((buffer[0]=='g') || (buffer[0]=='G')) {
    DrawData.modify = TRUE;
    /* Toggle */
    if (DrawData.boxType==1) {/* rect->unrect */
      DrawData.boxType = 3; 
      twin[0] = DrawData.box[0]+1; twin[1] = DrawData.box[1]+1;
      twin[2] = DrawData.box[2]+1; twin[3] = DrawData.box[3]+1;
      ObitDConCleanWindowUpdate (DrawData.myWindow, DrawData.field,
				 DrawData.Id, OBIT_DConCleanWindow_unrectangle,
				 twin, err);
       g_snprintf (errMsg, 120, "Convert window %d to rect unbox", DrawData.field);
   }
    else if (DrawData.boxType==2) { /* round->unround */
      DrawData.boxType = 4; 
      twin[0] = DrawData.box[1];
      twin[1] = DrawData.box[2]+1; twin[2] = DrawData.box[3]+1;
      ObitDConCleanWindowUpdate (DrawData.myWindow, DrawData.field,
				 DrawData.Id, OBIT_DConCleanWindow_unround,
				 twin, err);
       g_snprintf (errMsg, 120, "Convert window %d to round unbox", DrawData.field);
   }
    else if (DrawData.boxType==3) { /* unrect->rect */
      DrawData.boxType = 1; 
      twin[0] = DrawData.box[0]+1; twin[1] = DrawData.box[1]+1;
      twin[2] = DrawData.box[2]+1; twin[3] = DrawData.box[3]+1;
      ObitDConCleanWindowUpdate (DrawData.myWindow, DrawData.field,
				 DrawData.Id, OBIT_DConCleanWindow_rectangle,
				 twin, err);
      g_snprintf (errMsg, 120, "Convert window %d to rect box", DrawData.field);
    }
    else if (DrawData.boxType==4) {
      DrawData.boxType = 1; /* unround->round */
      twin[0] = DrawData.box[1];
      twin[1] = DrawData.box[2]+1; twin[2] = DrawData.box[3]+1;
      ObitDConCleanWindowUpdate (DrawData.myWindow, DrawData.field,
				 DrawData.Id, OBIT_DConCleanWindow_round,
				 twin, err);
      g_snprintf (errMsg, 120, "Convert window %d to round box", DrawData.field);
    }
    MessageShow(errMsg);
    /* erase */
    DrawAnyBox (&image[CurImag], DrawData.box, 0);
    /* Draw in new color */
    color = DrawData.IDdata->gcolor;
    /* different color for unboxes */
    if ((DrawData.boxType==3) || (DrawData.boxType==4)) 
      color = (DrawData.IDdata->gcolor+1)%8;
    DrawAnyBox (&image[CurImag], DrawData.box, color);
    PaintImage(IDdata); /* redraw image */
    /* Any errors? */
    if (err->error) ObitErrLog(err);
    return;
  } /* end Gg */
  
} /* end DrawBoxKeyEH */

