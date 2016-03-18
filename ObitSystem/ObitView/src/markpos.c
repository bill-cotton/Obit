/* $Id$ */
/* mark position dialog box  for ObitView */
/*-----------------------------------------------------------------------
*  Copyright (C) 1996,1997,1999,2002-2016
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
#include <Xm/Xm.h> 
#include <Xm/DialogS.h> 
#include <Xm/MainW.h> 
#include <Xm/ScrollBar.h>
#include <Xm/Form.h>
#include <Xm/PushB.h>
#include <Xm/Label.h>
#include <Xm/ToggleB.h>
#include <Xm/RowColumn.h>
#include <Xm/Separator.h>
#include <Xm/MessageB.h>
#include <Xm/TextF.h>
#include <Xm/Text.h>
#include <Xm/FileSB.h>
#include "imagedisp.h"
#include "poslabel.h"
#include "markpos.h"
#include "messagebox.h"
#include "graph.h"
#include "ObitSkyGeom.h"
#define min(a,b) (((a) < (b)) ? (a) : (b))
#define max(a,b) (((a) > (b)) ? (a) : (b))

/**
 *  \file positionbox.c
 * displays "position" dialog.
 */

/*--------------- file global data ----------------*/
/** is the position box active? */
olong PositionBoxActive = 0;

/* global structure for things to talk to each other */
typedef struct {
  Widget       dialog, posfilebox, equlab, label1, label2;
  Widget       data1, data2, data3; /* RA, Dec, size text windows */
  ImageDisplay *IDdata;
  olong         markx, marky, iInner, iOuter, bPixel;
  odouble      ra, dec;
} MarkPosStuff;
MarkPosStuff mark;

/**
 * Read RA from text window
 * \param w           text widget
 * \return return 0 if OK, 1 if error
 */
int ReadPosRA (Widget w)
{
  gchar     *value=NULL, ctype[5];
  olong      h, m, toHours, iNumRead, bOK;
  ofloat    s;
  /*Display *dpy = XtDisplay(w);*/
  odouble   ra;
  
  /* read value */
  value = XmTextGetString (w);
  if (!value) /* error */
    {MessageShow ("Error reading RA for Mark Position");
    return 1;}
  /*  first axis in hours or degrees? */
  strncpy (ctype, image[CurImag].myDesc->ctype[0], 4);
  toHours = (!strncmp(ctype, "RA", 2)) || (!strncmp(ctype, "LL", 2));
  iNumRead = sscanf (value, "%d %d %f", &h, &m, &s);
  if (value) XtFree(value); value = NULL;
  bOK = 0;
  mark.bPixel = 0; /* is this a pixel number rather than a real coordinate? */
  /* Deal with just a single value for pixel number */
  if (iNumRead==1) {
    bOK = 1;
    mark.bPixel = 1; /* coordinates are pixel values */
    ra = h;
  }
  if (toHours)
    {if (iNumRead==3) bOK = !hmsra(h, m, s, &ra);}
  else
    if (iNumRead==3) bOK = !dmsdec(h, m, s, &ra);
  if (!bOK)
    { /* error */
      MessageShow ("Error reading RA for Mark Position");
      return 1;}
  
  /* OK, save */
  mark.ra = ra;
  return 0;
} /* end ReadPosRA */

/**
 * Read Dec from text window
 * \param w           text widget
 * \return return 0 if OK, 1 if error
 */
int ReadPosDec (Widget w)
{
  gchar     *value=NULL;
  olong      d, m, iNumRead, bOK, bSouth, i;
  ofloat    s;
  odouble   dec;
  /*Display *dpy = XtDisplay(w);*/
  
  /* read value */
  value = XmTextGetString (w);
  if (!value) /* error */
    {MessageShow("Error reading Dec for Mark Position");
    return 1;}
  bOK = 0;
  iNumRead = sscanf (value, "%d %d %f", &d, &m, &s);
  /* southern declination? look for minus sign*/
  bSouth = 0;
  for (i=0;i<20;i++) 
    {if (value[i]=='-') bSouth=1; 
    if (value[i]==0) break;}
  if (value) XtFree(value); value = NULL;
  /* Deal with just a single value for pixel number */
  if (iNumRead==1) {
    bOK = 1;
    mark.bPixel = 1; /* coordinates are pixel values */
    dec = d;
  }
  if (iNumRead==3) bOK = !dmsdec(d, m, s, &dec);
  if (!bOK)
    { /* error */
      MessageShow("Error reading Dec for Mark Position");
      return 1;}
  
  if (bSouth && (dec>0.0)) dec = -dec; /* make sure declination sign OK */
  
  /* OK, save */
  mark.dec = dec;
  return 0;
} /* end ReadPosDec */

int ReadPosRADec (Widget w)
{
  gchar     *value=NULL, ctype[5];
  olong      h, rm, dd, dm, toHours, iNumRead, bOK, bSouth, i;
  ofloat    rs, ds;
  /*Display *dpy = XtDisplay(w);*/
  odouble   ra, dec;
  
  /* read value */
  value = XmTextGetString (w);
  if (!value) /* error */
    {MessageShow ("Error reading RA for Mark Position");
    return 1;}
  /*  first axis in hours or degrees? */
  strncpy (ctype, image[CurImag].myDesc->ctype[0], 4);
  toHours = (!strncmp(ctype, "RA", 2)) || (!strncmp(ctype, "LL", 2));
  iNumRead = sscanf (value, "%d %d %f %d %d %f", &h, &rm, &rs, &dd, &dm, &ds);
  /* southern declination? look for minus sign*/
  bSouth = 0;
  for (i=0;i<strlen(value);i++) 
    {if (value[i]=='-') bSouth=1; 
    if (value[i]==0) break;}
 if (value) XtFree(value); value = NULL;
  bOK = 0;
  mark.bPixel = 0; /* is this a pixel number rather than a real coordinate? */
  /* Deal with just a single value for pixel number */
  if (iNumRead==2) {
    bOK = 1;
    mark.bPixel = 1; /* coordinates are pixel values */
    ra = h;
  }
  if (toHours)
    {if (iNumRead>=3) bOK = !hmsra(h, rm, rs, &ra);}
  else
    if (iNumRead>=3) bOK = !dmsdec(h, rm, rs, &ra);
  if (!bOK)
    { /* error */
      MessageShow ("Error reading RA for Mark Position");
      return 1;}
  
  /* OK, save */
  mark.ra = ra;

  /* Declination */
  if (iNumRead==2) {
    bOK = 1;
    mark.bPixel = 1; /* coordinates are pixel values */
    dec = rm;
  }
  if (iNumRead>=6) bOK = !dmsdec(dd, dm, ds, &dec);
  if (!bOK)
    { /* error */
      MessageShow("Error reading Dec for Mark Position");
      return 1;}
  
  if (bSouth && (dec>0.0)) dec = -dec; /* make sure declination sign OK */
  
  /* OK, save */
  mark.dec = dec;
  return 0;
} /* end ReadPosRADec */

/**
 *  Read size from text window
 * \param w           text widget
 * \return return 0 if OK, 1 if error
 */
int ReadPosSize (Widget w)
{
  gchar     *value=NULL;
  olong      iInner, iOuter, iNumRead;
  /*Display *dpy = XtDisplay(w);*/
  
  /* read value */
  value = XmTextGetString (w);
  if (!value) /* error */
    {MessageShow("Error reading cross size for Mark Position");
    return 1;}
  iNumRead = sscanf (value, "%d %d", &iInner, &iOuter);
  if (value) XtFree(value); value = NULL;
  if (iNumRead!=2)
    { /* error */
      MessageShow("Error reading cross size for Mark Position");
      return 1;}
  
  /* OK, save */
  mark.iInner = iInner;
  mark.iOuter = iOuter;
  return 0;
} /* end ReadPosSize */


/**
 * Callback for file selection button
 * \param w           widget activated
 * \param clientData  client data
 * \param callData    call data
 */
/* file selection button callbacks */
void PosFileOKCB (Widget filebox, XtPointer clientData, XtPointer callData)
{
  gchar *filename=NULL, *directory=NULL;
  XmFileSelectionBoxCallbackStruct *cbs;
  ImageDisplay  *IDdata;
  /*Display *dpy = XtDisplay(filebox);*/
  
  cbs = (XmFileSelectionBoxCallbackStruct *) callData;
  IDdata = (ImageDisplay *)clientData;
  
  /* get file name */
  if (!XmStringGetLtoR (cbs->value, XmSTRING_DEFAULT_CHARSET, &filename))
    return; /* error */
  
  /* get directory name */
  if (!XmStringGetLtoR (cbs->dir, XmSTRING_DEFAULT_CHARSET, &directory))
    return; /* error */
  if (!mark_dir) mark_dir = MakeFStrng(" ");
  FStrngFill (mark_dir, directory);
  
  /* process file marking positions */
  
  MarkPosFromFile (filename, IDdata);
  
  /*clean up */
  if (filename) XtFree(filename); filename = NULL;
  if (directory) XtFree(directory); directory = NULL;
  
  /* make it disappear but still exist */
  XtPopdown (XtParent (filebox)); 
} /* end FileButCB */

/* button callbacks */
void PosFileCancelCB (Widget filebox, XtPointer clientData, XtPointer callData)
     /* File select cancel button hit- bail out */
{
  /* make it disappear but still exist */
  XtPopdown (XtParent (filebox)); 
} /* end FileButCB */

/* button callbacks */
/**
 * Callback for File button hit - get info from a file 
 * \param w           widget activated
 * \param clientData  client data
 * \param callData    call data
 */
void PosFileButCB (Widget w, XtPointer clientData, XtPointer callData)
{
  ImageDisplay *IDdata;
  XmString     wierdstring;
  Arg          wargs[5]; 
  
  IDdata = (ImageDisplay *)clientData;
  
  mark.posfilebox = (Widget) XmCreateFileSelectionDialog (w, "filebox", 
							  NULL, 0);
  XtAddCallback (mark.posfilebox, XmNokCallback, PosFileOKCB, clientData);
  XtAddCallback (mark.posfilebox, XmNcancelCallback, PosFileCancelCB, 
		 clientData);
  
  /* set directory if it is defined */
  if (mark_dir) {
    wierdstring = XmStringCreateSimple (mark_dir->sp);
    XtSetArg (wargs[0], XmNdirectory, wierdstring);
    XtSetValues (mark.posfilebox, wargs, 1);
    if (wierdstring) XmStringFree(wierdstring); wierdstring = NULL;
  }
  /* Shazam appear: */
  XtManageChild (mark.posfilebox);
  XtPopup (XtParent (mark.posfilebox), XtGrabNone);
  /* all the action is in the callback routine */
  
} /* end PosFileButCB */

/**
 * Callback for Mark button hit - get info and mark
 * \param w           widget activated
 * \param clientData  client data
 * \param callData    call data
 */
void PosOKButCB (Widget w, XtPointer clientData, XtPointer callData)
{
  olong iRet, bOK=0, bDec=0, iX, iY, iScroll;
  ofloat xp, yp;
  gchar ctype[5], szErrMess[120], *value=NULL;
  /*Display *dpy = XtDisplay(w);*/
  ImageDisplay *IDdata = (ImageDisplay*)clientData;
  ObitImageDesc *desc;
  
  /* read selected values  - is Dec given? */
  value = XmTextGetString (mark.data2);
  if (!value) /* error */
    {MessageShow("Error reading Dec for Mark Position");
    return;}
  bDec = strlen(value)>2;
  if (bDec) {  /* Read RA, Dec separately */
  if (ReadPosRA(mark.data1)) return;
  if (ReadPosDec(mark.data2)) return;
  } else { /* Read RA, Dec from RA line */
    if (ReadPosRADec(mark.data1)) return;
  }
  if (ReadPosSize(mark.data3)) return;
  
  /* convert to pixel */
  if (mark.bPixel) { /* position pixel values */
    xp = mark.ra;
    yp = mark.dec;
  } else { /* convert position to pixel */
    strncpy (ctype, &image[CurImag].myDesc->ctype[0][4], 4); 
    ctype[4] = 0;
   desc = image[CurImag].myDesc;
   iRet = ObitSkyGeomXYpix(mark.ra, mark.dec, 
			   desc->crval[desc->jlocr], desc->crval[desc->jlocd],
			   desc->crpix[desc->jlocr], desc->crpix[desc->jlocd],
			   desc->cdelt[desc->jlocr], desc->cdelt[desc->jlocd],
			   desc->crota[desc->jlocd], &desc->ctype[desc->jlocr][4],
			   &xp, &yp);
    bOK = (iRet==0);
  }

  /* "Y" axis upside down */
  yp = (ofloat)image[CurImag].myDesc->inaxes[1] - yp;
  
  /*  check that pixel is in image. */
  bOK = bOK && (xp>0.0) && (yp>0.0) && 
    (xp<=(ofloat)image[CurImag].myDesc->inaxes[0]) &&
    (yp<=(ofloat)image[CurImag].myDesc->inaxes[1]);

  /* everything OK? */
  if (!bOK) 
    {sprintf (szErrMess, "invalid position to mark %9.2f, %9.2f", xp, yp);
    MessageShow (szErrMess);
    return;}
  
  /* mark position in pixel array; convert to 0-rel */
  iX = (olong)(xp - 0.5);
  iY = (olong)(yp - 0.5);
  
  MarkPix (mark.IDdata, iX, iY, mark.iInner, mark.iOuter);
  
  /* for scroll flip y axis back */
  yp = (ofloat)image[CurImag].myDesc->inaxes[1] - yp;
  iY = (olong)(yp - 0.5);
  
  /* set scroll for last position */
  iScroll = image[CurImag].myDesc->inaxes[1] - iY;
  if (IDdata->vscr_vis) 
    {
      IDdata->scrolly = iScroll;
      iScroll -= IDdata->vscr_half;
      iScroll = min (iScroll, IDdata->vscr_max);
      iScroll = max (iScroll, IDdata->vscr_min);
      /* set scrollbar */ 
      XtVaSetValues(IDdata->vscroll,  XmNvalue, iScroll, NULL);
    }
  if (IDdata->hscr_vis) 
    {
      iScroll = iX - IDdata->hscr_half;
      iScroll = min (iScroll, IDdata->hscr_max);
      iScroll = max (iScroll, IDdata->hscr_min);
      IDdata->scrollx = iX;
      /* set scroll bar */
      XtVaSetValues(IDdata->hscroll,  XmNvalue, iScroll, NULL);
    }
  
  /* repaint image */
  PaintImage (mark.IDdata);
  
} /* end PosOKButCB */

/**
 * Callback for Clear button hit - clear position fields
 * \param w           widget activated
 * \param clientData  client data
 * \param callData    call data
 */
void ClearOKButCB (Widget w, XtPointer clientData, XtPointer callData)
{
  gchar valuestr[100];
  /*Display *dpy = XtDisplay(w);*/
  /*ImageDisplay *IDdata = (ImageDisplay*)clientData;*/
  
  sprintf (valuestr, " ");
  /* RA = mark.data1 */
  XmTextSetString (mark.data1, valuestr);
  /* Dec = mark.data2 */
  XmTextSetString (mark.data2, valuestr);
  
} /* end ClearOKButCB */

/**
 * Callback for Cancel button hit
 * \param w           widget activated
 * \param clientData  client data
 * \param callData    call data
 */
void PosCancelButCB (Widget w, XtPointer clientData, XtPointer callData)
{
  /* make it disappear but still exist */
  XtPopdown(mark.dialog);
} /* end PosCancelButCB */

/**
 * Callback for reate dialog box for specifying positions to mark
 * \param w           widget activated
 * \param clientData  client data
 * \param callData    call data
 */
void MarkPosCB (Widget parent, XtPointer clientData, XtPointer callData)
{
  Widget form, toplabel, label3;
  Widget sep1, sep2, sep3;
  Widget FileButton, OKButton, CancelButton, ClearButton;
  XmString     RA=NULL, Dec=NULL,label=NULL, sizestr=NULL, equstr=NULL;
  XmString     wierdstring = NULL;
  Arg          wargs[5]; 
  gchar        valuestr[100], ctype[5];
  olong         h, d, m, toHours;
  gshort       xpos, ypos;
  ofloat       s;
  ImageDisplay *IDdata = (ImageDisplay*)clientData;
  /*Display *dpy = XtDisplay(parent);*/
  
  
  /* validity checks */
  if (!IDdata) return;
  if (!image[CurImag].valid) return;
  
  /* register IDdata */
  mark.IDdata = IDdata;
  
  /* don't make another one */
  if (PositionBoxActive) {
    if (XtIsRealized (mark.dialog))
      XMapRaised (XtDisplay(IDdata->shell), XtWindow(mark.dialog));
    
    /* bring it back where we can see it */
    XtPopup(mark.dialog, XtGrabNonexclusive);
    
    /* put it some place reasonable */
    /*  where is parent? */
    XtVaGetValues (IDdata->shell, XmNx, &xpos, XmNy, &ypos,  NULL);
    ypos += 160;
    if (xpos<0) xpos = 0;
    XMoveWindow (XtDisplay(IDdata->shell), XtWindow(mark.dialog), 
		 xpos, ypos);
    
    /* set values */
    /* equinox - labels are difficult */
    if ((usr_equinox>0.0) && (image[CurImag].myDesc->equinox>0.0))
      sprintf(valuestr,"Equinox %7.1f", usr_equinox);
    else if (image[CurImag].myDesc->equinox>0.0)
      sprintf(valuestr,"Equinox %7.1f", image[CurImag].myDesc->equinox);
    else
      sprintf(valuestr,"Equinox unknown");
    wierdstring = XmStringCreateSimple (valuestr);
    XtSetArg (wargs[0], XmNlabelString, wierdstring);
    XtSetValues (mark.equlab, wargs, 1);
    if (wierdstring) XmStringFree(wierdstring); wierdstring = NULL;
    
    /* set RA, Dec labels from FITS file info */
    strncpy (ctype, image[CurImag].myDesc->ctype[0], 4); ctype[4] = 0;
    if (ctype[2]=='-') ctype[2]=' '; if (ctype[3]=='-') ctype[3]=' ';
    wierdstring = XmStringCreateSimple (ctype);
    XtSetArg (wargs[0], XmNlabelString, wierdstring);
    XtSetValues (mark.label1, wargs, 1);
    if (wierdstring) XmStringFree(wierdstring); wierdstring = NULL;
    strncpy (ctype, image[CurImag].myDesc->ctype[1], 4); ctype[4] = 0;
    if (ctype[2]=='-') ctype[2]=' '; if (ctype[3]=='-') ctype[3]=' ';
    wierdstring = XmStringCreateSimple (ctype);
    XtSetArg (wargs[0], XmNlabelString, wierdstring);
    XtSetValues (mark.label2, wargs, 1);
    if (wierdstring) XmStringFree(wierdstring); wierdstring = NULL;
    
    
    /* Position */
    /*  first axis in hours or degrees? */
    strncpy (ctype, image[CurImag].myDesc->ctype[0], 4);
    toHours = (!strncmp(ctype, "RA", 2)) || (!strncmp(ctype, "LL", 2));
    /* default is reference position */
    mark.ra  = image[CurImag].myDesc->crval[0];
    mark.dec = image[CurImag].myDesc->crval[1];
    
    /* precess if necessary */
    if ((usr_equinox > 0.0) && 
	(image[CurImag].usr_equinox !=
	 image[CurImag].myDesc->equinox))
      {if (image[CurImag].myDesc->equinox==1950.0) 
	ObitSkyGeomBtoJ (&mark.ra, &mark.dec);
      if (image[CurImag].myDesc->equinox==2000.0) 
	ObitSkyGeomJtoB (&mark.ra, &mark.dec);}
    
    /* RA */
    if (toHours)
      {rahms (mark.ra, &h, &m, &s);
      sprintf (valuestr, "%2.2d %2.2d %5f", h, m, s);}
    else
      {decdms (mark.ra, &h, &m, &s);
      sprintf (valuestr, "%3.2d %2.2d %5f", h, m, s);
      }
    XmTextSetString (mark.data1, valuestr);
    
    /* declination */
    decdms (mark.dec, &d, &m, &s);
    if (d<0) d = -d;
    sprintf (valuestr, "%3.2d %2.2d %5f", d, m, s);
    if (mark.dec<0.0) valuestr[0]='-'; /* neg 0 dec. */
    XmTextSetString (mark.data2, valuestr);
    
    return;
  } /* end of update dialog */
  
  label = XmStringCreateSimple ("Position to mark on image");
  sizestr = XmStringCreateSimple ("size");
  /* mark as active */
  PositionBoxActive = 1;
  
  mark.dialog = XtVaCreatePopupShell ("MarkPos", xmDialogShellWidgetClass, 
				      IDdata->shell, 
				      XmNautoUnmanage, False,
				      XmNwidth,     150,
				      XmNheight,    190,
				      XmNdeleteResponse, XmDESTROY,
				      NULL);
  
  /* make Form widget to stick things on */
  form = XtVaCreateManagedWidget ("MarkForm", xmFormWidgetClass,
				  mark.dialog,
				  XmNautoUnmanage, False,
				  XmNwidth,     150,
				  XmNheight,    190,
				  XmNx,           0,
				  XmNy,           0,
				  NULL);
  
  /* info label widgets */
  toplabel = XtVaCreateManagedWidget ("Label", xmLabelWidgetClass, 
				      form, 
				      XmNwidth,           150,
				      XmNlabelString,   label,
				      XmNtopAttachment, XmATTACH_FORM,
				      XmNleftAttachment,  XmATTACH_FORM,
				      NULL);
  
  /* Equinox */
  if ((usr_equinox>0.0) && (image[CurImag].myDesc->equinox>0.0))
    sprintf(valuestr,"Equinox %7.1f", usr_equinox);
  else if (image[CurImag].myDesc->equinox>0.0)
    sprintf(valuestr,"Equinox %7.1f", image[CurImag].myDesc->equinox);
  else
    sprintf(valuestr,"Equinox unknown");
  equstr = XmStringCreateSimple (valuestr);
  mark.equlab = XtVaCreateManagedWidget ("EquLabel", xmLabelWidgetClass, 
					 form, 
					 XmNwidth,           150,
					 XmNlabelString,   equstr,
					 XmNtopAttachment, XmATTACH_WIDGET,
					 XmNtopWidget,     toplabel,
					 XmNleftAttachment,  XmATTACH_FORM,
					 NULL);
  
  /* set labels from FITS file info */
  strncpy (ctype, image[CurImag].myDesc->ctype[0], 4); ctype[4] = 0;
  if (ctype[2]=='-') ctype[2]=' '; if (ctype[3]=='-') ctype[3]=' ';
  RA = XmStringCreateSimple (ctype);
  strncpy (ctype, image[CurImag].myDesc->ctype[1], 4); ctype[4] = 0;
  if (ctype[2]=='-') ctype[2]=' '; if (ctype[3]=='-') ctype[3]=' ';
  Dec = XmStringCreateSimple (ctype);
  
  /* RA */
  mark.label1 = XtVaCreateManagedWidget ("RA", xmLabelWidgetClass,
					 form,
					 XmNheight,    30,
					 XmNlabelString,   RA,
					 XmNtopAttachment, XmATTACH_WIDGET,
					 XmNtopWidget,     mark.equlab,
					 XmNleftAttachment,  XmATTACH_FORM,
					 NULL);
  /*  first axis in hours or degrees? */
  strncpy (ctype, image[CurImag].myDesc->ctype[0], 4);
  toHours = (!strncmp(ctype, "RA", 2)) || (!strncmp(ctype, "LL", 2));
  /* default is reference position */
  mark.ra  = image[CurImag].myDesc->crval[0];
  mark.dec = image[CurImag].myDesc->crval[1];
  
  /* precess if necessary */
  if ((usr_equinox > 0.0) && 
      (image[CurImag].usr_equinox !=
       image[CurImag].myDesc->equinox))
    {if (image[CurImag].myDesc->equinox==1950.0) 
      ObitSkyGeomBtoJ (&mark.ra, &mark.dec);
    if (image[CurImag].myDesc->equinox==2000.0) 
      ObitSkyGeomJtoB (&mark.ra, &mark.dec);}
  
  if (toHours)
    {rahms (mark.ra, &h, &m, &s);
    sprintf (valuestr, "%2.2d %2.2d %5f", h, m, s);}
  else
    {decdms (mark.ra, &h, &m, &s);
    sprintf (valuestr, "%3.2d %2.2d %5f", h, m, s);
    }
  mark.data1 = XtVaCreateManagedWidget ("RA data", xmTextFieldWidgetClass, 
					form, 
					XmNvalue,        valuestr,
					XmNtopAttachment, XmATTACH_WIDGET,
					XmNtopWidget,     mark.equlab,
					XmNleftAttachment, XmATTACH_WIDGET,
					XmNleftWidget,     mark.label1,
					XmNrightAttachment,  XmATTACH_FORM,
					NULL);
  /* separator */
  sep1 = XtVaCreateManagedWidget ("sep1", xmSeparatorWidgetClass,
				  form, 
				  XmNwidth,           150,
				  XmNtopAttachment, XmATTACH_WIDGET,
				  XmNtopWidget,     mark.data1,
				  XmNleftAttachment,  XmATTACH_FORM,
				  NULL);
  /* Dec */
  mark.label2 = XtVaCreateManagedWidget ("Dec", xmLabelWidgetClass,
					 form,
					 XmNheight,    30,
					 XmNlabelString,   Dec,
					 XmNtopAttachment, XmATTACH_WIDGET,
					 XmNtopWidget,     sep1,
					 XmNleftAttachment,  XmATTACH_FORM,
					 NULL);
  
  decdms (mark.dec, &d, &m, &s);
  if (d<0) d = -d;
  sprintf (valuestr, "%3.2d %2.2d %5f", d, m, s);
  if (mark.dec<0.0) valuestr[0]='-'; /* neg 0 dec. */
  mark.data2 = XtVaCreateManagedWidget ("Dec data", xmTextFieldWidgetClass, 
					form, 
					XmNvalue,        valuestr,
					XmNtopAttachment, XmATTACH_WIDGET,
					XmNtopWidget,     sep1,
					XmNleftAttachment, XmATTACH_WIDGET,
					XmNleftWidget,     mark.label2,
					XmNrightAttachment,  XmATTACH_FORM,
					NULL);
  
  /* separator */
  sep2 = XtVaCreateManagedWidget ("sep2", xmSeparatorWidgetClass,
				  form, 
				  XmNwidth,           150,
				  XmNtopAttachment, XmATTACH_WIDGET,
				  XmNtopWidget,     mark.data2,
				  XmNleftAttachment,  XmATTACH_FORM,
				  NULL);
  /* size */
  label3 = XtVaCreateManagedWidget ("pos size", xmLabelWidgetClass,
				    form,
				    XmNheight,    30,
				    XmNlabelString,   sizestr,
				    XmNtopAttachment, XmATTACH_WIDGET,
				    XmNtopWidget,     sep2,
				    XmNleftAttachment,  XmATTACH_FORM,
				    NULL);
  
  mark.iInner = 6;  mark.iOuter = 15;
  sprintf (valuestr, " %d %d", mark.iInner, mark.iOuter);
  mark.data3 = XtVaCreateManagedWidget ("size data", xmTextFieldWidgetClass, 
					form, 
					XmNvalue,        valuestr,
					XmNtopAttachment, XmATTACH_WIDGET,
					XmNtopWidget,     sep2,
					XmNleftAttachment, XmATTACH_WIDGET,
					XmNleftWidget,     label3,
					XmNrightAttachment,  XmATTACH_FORM,
					NULL);
  /* separator */
  sep3 = XtVaCreateManagedWidget ("sep3", xmSeparatorWidgetClass,
				  form, 
				  XmNwidth,           150,
				  XmNtopAttachment, XmATTACH_WIDGET,
				  XmNtopWidget,     mark.data3,
				  XmNleftAttachment,  XmATTACH_FORM,
				  NULL);
  /* File button */
  FileButton = XtVaCreateManagedWidget ("File", xmPushButtonWidgetClass, 
					form, 
					XmNtopAttachment, XmATTACH_WIDGET,
					XmNtopWidget,     sep3,
					XmNleftAttachment,  XmATTACH_FORM,
					NULL);
  XtAddCallback (FileButton, XmNactivateCallback, PosFileButCB, 
		 (XtPointer)IDdata);
  
  /* Cancel button */
  CancelButton = XtVaCreateManagedWidget ("Cancel", xmPushButtonWidgetClass, 
					  form, 
					  XmNtopAttachment, XmATTACH_WIDGET,
					  XmNtopWidget,     sep3,
					  XmNleftAttachment, XmATTACH_WIDGET,
					  XmNleftWidget,     FileButton,
					  NULL);
  XtAddCallback (CancelButton, XmNactivateCallback, PosCancelButCB, 
		 (XtPointer)IDdata);
  
  /* Mark button */
  OKButton = XtVaCreateManagedWidget ("Mark", xmPushButtonWidgetClass, 
				      form, 
				      XmNtopAttachment, XmATTACH_WIDGET,
				      XmNtopWidget,     sep3,
				      XmNleftAttachment, XmATTACH_WIDGET,
				      XmNleftWidget,     CancelButton,
				      XmNrightAttachment, XmATTACH_FORM,
				      NULL);
  XtAddCallback (OKButton, XmNactivateCallback, PosOKButCB, (XtPointer)IDdata);
  
  /* Clear position button */
  ClearButton = XtVaCreateManagedWidget ("Clear", xmPushButtonWidgetClass, 
					 form, 
					 XmNtopAttachment, XmATTACH_WIDGET,
					 XmNtopWidget,     FileButton,
					 XmNbottomAttachment, XmATTACH_FORM,
					 XmNleftAttachment, XmATTACH_FORM,
					 NULL);
  XtAddCallback (ClearButton, XmNactivateCallback, ClearOKButCB, (XtPointer)IDdata);

  if (label) XmStringFree(label); label = NULL;
  if (equstr) XmStringFree(equstr); equstr = NULL;
  if (RA) XmStringFree(RA); RA = NULL;
  if (Dec) XmStringFree(Dec); Dec = NULL;
  if (sizestr) XmStringFree(sizestr); sizestr = NULL;
  
  /* set it up */
  XtManageChild (mark.dialog);
  
  /* put it some place reasonable */
  /*  where is parent? */
  XtVaGetValues (IDdata->shell,
		 XmNx, &xpos,
		 XmNy, &ypos,
		 NULL);
  ypos += 180;
  if (xpos<0) xpos = 0;
  XMoveWindow (XtDisplay(IDdata->shell), XtWindow(mark.dialog), 
	       xpos, ypos);
} /* end MarkPosCB */

/**
 * Routine to read next position from file and return pixel position and 
 * cross size information.
 * \param hFile     file to read, expected open on call and is closed if EOF or error. 
 * \param ImageDesc Image display
 * \param xpix      x pixel of center
 * \param ypix      y pixel of center
 * \param inner     inner size of cross (pixels)
 * \param outer     outer size of cross (pixels)
 * \param iErr      error return; 0=OK else an error occured;
 *        iErr = -1 means the position is out of the image.
 * \return 0 if OK else failed  
 */
int MarkPosRead(FILE *hFile, ImageData *ImageDesc, ofloat *xpix, 
		ofloat *ypix, olong *inner, olong *outer, olong *iErr)
{
  olong h, d, rm, dm, iNumByte, iNumRead, iRet, iInner, iOuter, toHours, i;
  ofloat rs, ds, xp, yp;
  odouble ra, dec;
  gboolean bRAOK, bDecOK, bOK, bSouth;
  gchar szLine[120],szErrMess[120],  ctype[5];
  ObitImageDesc *desc;
 
  if (!fgets (szLine, 120, hFile)) {
    /*   error or EOF */
    if (ferror(hFile)) *iErr = 2;  /* Check if error */
    fclose(hFile);  /* close */
    return 0;}  /*  end of error check */
  iNumByte = strlen(szLine);
  iNumRead = 0; bOK = False; bRAOK = False; bDecOK = False;
  if (iNumByte>10){ /* decode line */
    iNumRead = sscanf (szLine, "%d %d %f %d %d %f %d %d", 
		       &h, &rm, &rs, &d, &dm, &ds, &iInner, &iOuter);
    strncpy (ctype, ImageDesc->myDesc->ctype[0], 4);
    toHours = (!strncmp(ctype, "RA", 2))  || 
      (!strncmp(ctype, "LL", 2));
    if (toHours)
      {if (iNumRead>=3) bRAOK = !hmsra(h, rm, rs, &ra);}
    else
      if (iNumRead>=3) bRAOK = !dmsdec(h, rm, rs, &ra);
    if (iNumRead>=6) bDecOK = !dmsdec(d, dm, ds, &dec);
    
    /* if 2 values then they are pixel values */
    mark.bPixel = 0;
    if (iNumRead==2) {
      mark.bPixel = 1;
      mark.ra = h;
      mark.dec = rm;
    }
    
    /* southern declination? look for minus sign*/
    bSouth = 0;
    for (i=5;i<120;i++) 
      {if (szLine[i]=='-') bSouth=1; 
      if (szLine[i]==0) break; }
    if (bSouth && (dec>0.0)) dec = -dec;
    
    if (bRAOK && bDecOK) {  /* get pixel */
      strncpy (ctype, &ImageDesc->myDesc->ctype[0][4], 4);
      ctype[4] = 0;
      desc = image[CurImag].myDesc;
      iRet = ObitSkyGeomXYpix(ra, dec, 
			      desc->crval[desc->jlocr], desc->crval[desc->jlocd],
			      desc->crpix[desc->jlocr], desc->crpix[desc->jlocd],
			      desc->cdelt[desc->jlocr], desc->cdelt[desc->jlocd],
			      desc->crota[desc->jlocd], &desc->ctype[desc->jlocr][4],
			      &xp, &yp);
      bOK = (iRet==0);
    }  /* end of convert to pixel  */
    
    /* "Y" axis upside down */
    yp = (ofloat)image[CurImag].myDesc->inaxes[1] - yp;
  
   /*  Inner and outer read? */
    bOK = bOK && (iNumRead>=8) && (iInner>=0) && (iInner<100) &&
      (iOuter>0) && (iOuter<100);
  }  /* End of decode valid line */
  
  if (bOK)   /* everything OK? */
    {*xpix = xp; *ypix = yp;
    *inner = iInner; *outer = iOuter;
    /*  check that pixel is in image. */
    if ((xp>0.0) && (yp>0.0) && 
	(xp<(ofloat)(ImageDesc->myDesc->inaxes[0])) &&
	(yp<(ofloat)(ImageDesc->myDesc->inaxes[1]))) *iErr = 0;
    else *iErr = -1;  /* out of image */
    return 1;}
  else  /* bogus dudes */
    {*iErr = 3;
    fclose(hFile);  /* close */
    sprintf (szErrMess, "Error with file entry %s", szLine);
    MessageShow (szErrMess);
    return 0;}
}  /* end MarkPosRead */

/**
 * read file with marking positions and mark image
 * \param filename  Name of file with positions to mark
 * \param IDdata    Image display to update
 */
void MarkPosFromFile (char* filename, ImageDisplay* IDdata)
{
  olong iXcen=0, iYcen=0, iScroll, iMore, iErr, iNumOut=0;
  olong inner_size, outer_size;
  ofloat xpix, ypix;
  gchar szErrMess[120];
  FILE   *hFile;
  
  iErr = 0;
  /* Open file */
  hFile = fopen (filename, "rt");
  if (hFile==NULL) iErr = 1;
  
  /* loop over entries */
  iMore = (olong)(!iErr);
  while (iMore) {
    iMore = MarkPosRead(hFile, &image[CurImag], &xpix, &ypix, &inner_size, 
			&outer_size, &iErr);
    if ((!iMore) || (iErr>0)) break;
    /* determing center pixel */
    if (iErr==0) {  /* Only plot if in image otherwise just count */
      iXcen = (olong)(xpix - 0.5);  /* 0-rel posn */
      iYcen = (olong)(ypix - 0.5);
      MarkPix (IDdata, iXcen, iYcen, inner_size, outer_size);}
    else
      iNumOut++;
  }  /* end of loop over positions */
  /* set scroll for last position */
  if (IDdata->vscr_vis) 
    {
      iScroll = iYcen;
      IDdata->scrolly = iScroll;
      iScroll -= IDdata->vscr_half;
      iScroll = min (iScroll, IDdata->vscr_max);
      iScroll = max (iScroll, IDdata->vscr_min);
      /* set scrollbar */ 
      XtVaSetValues(IDdata->vscroll, XmNvalue, iScroll, NULL);
    }
  
  if (IDdata->hscr_vis) 
    {
      /*                      Take into account the edges of the image. */
      iScroll =iXcen  - IDdata->hscr_half;
      iScroll = min (iScroll, IDdata->hscr_max);
      iScroll = max (iScroll, IDdata->hscr_min);
      IDdata->scrollx = iXcen;
      /* set scroll bar */
      XtVaSetValues(IDdata->hscroll, XmNvalue, iScroll, NULL);
    }
  
  if (iNumOut)  /* Tell if some out of image */
    {sprintf (szErrMess, "%d positions were out of the image", iNumOut);
    MessageShow (szErrMess);}
  if (iErr>0)
    {sprintf (szErrMess, "Error %d with file %s", iErr, filename);
    MessageShow (szErrMess);}
  
  /* redraw image */
  PaintImage(IDdata);
  
} /* end MarkPosFromFile */

#include "drawbox.h"   /* DEBUG */
#include "ObitDConCleanWindow.h"   /* DEBUG */
#include "ObitMem.h"   /* DEBUG */
/**
 * Routine to mark cross on gpharray, iX, iY are 0-rel.
 * \param IDdata   Image display
 * \param iX       x pixel of center 0-rel
 * \param iY       y pixel of center 0-rel
 * \param iInner   number of pixels from center on each arm of cross not to mark
 * \param iOuter   half width of cross in pixels
 */
void MarkPix (ImageDisplay *IDdata, olong iX, olong iY, olong iInner, olong iOuter)
{
  olong iiY, iNx, iNy, i1, i2, i3, i4, j1, j2, j3, j4;
  olong blc[2], trc[2];
  olong box[4], nax[2];
  ObitDConCleanWindow* myWindow;   /* DEBUG */
  
  iNx = image[CurImag].myDesc->inaxes[0];
  iNy = image[CurImag].myDesc->inaxes[1];
  iiY = iNy - iY - 1; /* image "upside down" */

  /* DEBUG - draw circle - test edit box if iInner <= 0 */
  if (iInner <= 0) {
    GraphClear(&image[CurImag]);
    box[0] = iOuter; box[1] = iX; box[2] = iY; box[3] = iOuter;
    nax[0] = iNx; nax[1] = iNy; 
    myWindow = ObitDConCleanWindowCreate1 ("DEBUGCleanBox", nax, err);
    myWindow->nfield = 1;
    myWindow->ndim   = 2;
    ObitDConCleanWindowAdd (myWindow, 1, OBIT_DConCleanWindow_round, box, err);
    DrawBox (IDdata, myWindow, 1);
    return;
  }
  /* done DEBUG */
  
  /*  horizional */
  i1 = iX - iOuter; i1 = max (0, i1);
  i2 = iX - iInner; i2 = max (0, i2);
  i3 = iX + iInner; i3 = min (iNx-1, i3);
  i4 = iX + iOuter; i4 = min (iNx-1, i4);
  blc[0] = i1; blc[1] = iiY;
  trc[0] = i2; trc[1] = iiY;
  GraphDrawLine(&image[CurImag], blc, trc, IDdata->gcolor);
  blc[0] = i3; blc[1] = iiY;
  trc[0] = i4; trc[1] = iiY;
  GraphDrawLine(&image[CurImag], blc, trc, IDdata->gcolor);
  
  /*  vertical */
  j1 = iiY - iOuter; j1 = max (0, j1);
  j2 = iiY - iInner; j2 = max (0, j2);
  j3 = iiY + iInner; j3 = min (iNy-1, j3);
  j4 = iiY + iOuter; j4 = min (iNy-1, j4);
  blc[0] = iX; blc[1] = j1;
  trc[0] = iX; trc[1] = j2;
  GraphDrawLine(&image[CurImag], blc, trc, IDdata->gcolor);
  blc[0] = iX; blc[1] = j3;
  trc[0] = iX; trc[1] = j4;
  GraphDrawLine(&image[CurImag], blc, trc, IDdata->gcolor);
}  /* End of MarkPix */
