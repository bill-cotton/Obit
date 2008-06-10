/* $Id$  */
/* position logging  for ObitView */
/*-----------------------------------------------------------------------
*  Copyright (C) 1996,1997-2008
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
#include "logger.h"
#include "menu.h"
#include "imagedisp.h"
#include "poslabel.h"
#include "messagebox.h"
#include "ObitSkyGeom.h"

/*--------------- file global data ----------------*/
/** Text widget for display */
TextFilePtr LogText = NULL; /* logging TextFile */

/*---------------Private function prototypes----------------*/
/* initialize log file, argument not used */
void LoggerInit (XPointer arg);
/* file selection canceled, arg not used */
void LoggerCancel (XPointer arg);

/*----------------------Public functions---------------------------*/

/**
 * Callback for toggle position logging on/off 
 * \param w           widget activated
 * \param clientData  ImageDisplay
 * \param callData    call data
 */
void LoggerCB (Widget w, XtPointer clientData, XtPointer callData)
{
  ImageDisplay  *IDdata= (ImageDisplay *)clientData;
  /* turn on or off ?*/
  if (doLog)
    { /* turn off */
      if (!LogText) return; /* bail out if no log TextFile */
      doLog = 0;
      MenuMarkLogger (1); /* reset menu item label */
      TextFileKill(LogText); /* delete structures */
      LogText = NULL;
    }
  else
    {/* turn on */
      if (LogText) TextFileKill (LogText); /* shouldn't need this */
      if (!log_dir) log_dir = MakeFStrng("");
      LogText = TextFileMake (IDdata->shell, NULL, log_dir->sp);
      /* ask for file name and initialize file */
      TextFileFind (2, LogText, LoggerInit, LoggerCancel);
      doLog = 1;
      MenuMarkLogger (2); /* reset menu item label */
    }
} /* end LoggerCB */

/**
 * Log position 
 */
void DoLogger(void)
{
  gboolean   valid, posOK;
  olong      ipos[7], ndim, i;
  gchar       flux[20], equistr[7], pixel[23], line[100];
  gchar       axtype[3][9], label[3][30];
  odouble     pos[3];
  ofloat      pix[3], fblank, val, *pixP;
  ImageData  *Image = &image[CurImag];
  ObitImageDesc *desc;
 
  if (!Image) return; /* Validity check */
  if (!Image->valid) return; 
  if (!Image->myDesc) return; /* Validity checks */
  if (!LogText) return; 
  
  /* encode pixel location (convert to 1 relative) */
  sprintf (pixel, "(%7.2f,%7.2f,%4d)", 
	   Image->fXpixel+1.0, Image->fYpixel+1.0, Image->PlaneNo+1);
  
  /*  get flux  */
  ndim = Image->myDesc->naxis;  /* Number of dimensions  */
  fblank =  ObitMagicF();
  /*  fitted or current?  */
  if (Image->iFitted>0) /* fitted values  */
    {valid = 1;
    val = Image->fBpixel;}
  else  /* current position  */
    {val = 0;

    for (i=3; i<7; i++) ipos[i] = 0;
    /* Note: Image class pixel numbers are 0 rel; geometry is 1 rel.  */
    ipos[0] = Image->iXPixel; ipos[1] = Image->iYPixel;
    ipos[2] = Image->PlaneNo;
    /* get value */
    pixP =  ObitFArrayIndex (image[CurImag].myPixels, ipos);
    if (pixP) val  = *pixP; /* In image? */
    else val = fblank;
    valid = (ndim>1) && (val!=fblank);
    } /* end of read current pixel */
  if ((val==fblank) || (!valid)) /* blanked or invalid? - quit*/
    {return;}
  sprintf (flux, "%f",val);
  flux[7] = 0;   /* only 6 char */
  
  /* equinox */
  sprintf (equistr, "  ????");
  if ((usr_equinox>0.0) && (Image->myDesc->equinox>0.0))
    {if (usr_equinox==2000.0) sprintf (equistr, " J2000");
    if(usr_equinox==1950.0) sprintf (equistr, " B1950");}
  if ((usr_equinox<0.0) && (Image->myDesc->equinox>0.0))
    {if (Image->myDesc->equinox==2000.0) sprintf (equistr, " J2000");
    if (Image->myDesc->equinox==1950.0) sprintf (equistr, " B1950");}
  
  /* celestial position  */
  pix[0] = Image->fXpixel+1.0; pix[1] = Image->fYpixel+1.0;
  pix[2] = (float)(Image->PlaneNo)+1.0;
  strncpy(axtype[0], Image->myDesc->ctype[0], 8);
  strncpy(axtype[1], Image->myDesc->ctype[1], 8);
  axtype[0][8] = 0; axtype[1][8] = 0;
  desc = image[CurImag].myDesc;
  posOK = 
    !ObitSkyGeomWorldPos(pix[0], pix[1], 
			 desc->crval[desc->jlocr], desc->crval[desc->jlocd],
			 desc->crpix[desc->jlocr], desc->crpix[desc->jlocd],
			 desc->cdelt[desc->jlocr], desc->cdelt[desc->jlocd],
			 desc->crota[desc->jlocd], &desc->ctype[desc->jlocr][4],
			 &pos[0], &pos[1]);
  AxisLabel(pos[0], axtype[0], label[0]);  /* human readable  */
  label[1][0] = 0;
  if (ndim>=2) AxisLabel(pos[1], axtype[1], label[1]);
  label[0][17]=0; /* only 17 characters */
  label[1][17]=0; /* only 17 characters */
  if (!posOK)   /* valid position? if not quit */
    {return;} /* bail out - logging OK */
  
  /* Open Log File */
  if (TextFileOpen (LogText, 2) != 1) {
    MessageShow ("Error opening logging file");
    TextFileKill(LogText);
    LogText = NULL;
    doLog = 0;
    MenuMarkLogger (1); /* reset menu item label */
    return;
  } /* end of error trap */
  
  /* write output line */
  sprintf (line, "%s %17s,%17s,%6s,%6s,%s", 
	   pixel, label[0], label[1], equistr,flux, Image->FileName->sp);
  if (!TextFileWrite(LogText, line)) {
    MessageShow ("Error writing logging file");
    TextFileClose(LogText);
    TextFileKill(LogText);
    LogText = NULL;
    doLog = 0;
    MenuMarkLogger (1); /* reset menu item label */
    return;
  } /* end of error trap */
  TextFileClose(LogText); /* close file */
  return;
} /* end DoLogger */

/**
 * initialize log file, argument not used, called after file selected
 * \param arg   ignored
 */
void LoggerInit (XPointer arg)
{
  char line[100];
  
  /* Open/create Text File */
  if (TextFileOpen (LogText, 2) != 1) {
    MessageShow ("Error opening logging file");
    TextFileKill(LogText);
    LogText = NULL;
    doLog = 0;
    MenuMarkLogger (1); /* reset menu item label */
    return;
  } /* end of error trap */
  
  /* Add first entry */
  sprintf (line, "ObitView position logging: (pixel) celestial pos., equinox, value, filename");
  if (!TextFileWrite(LogText, line)) {
    MessageShow ("Error writing logging file");
    TextFileClose(LogText);
    TextFileKill(LogText);
    LogText = NULL;
    doLog = 0;
    MenuMarkLogger (1); /* reset menu item label */
    return;
  } /* end of error trap */
  TextFileClose(LogText); /* close file */
  
} /* end LoggerInit */

/**
 * File selection canceled
 * \param arg   ignored
 */
void LoggerCancel (XPointer arg)
{
  /* shut down, bail out */
  /*if (LogText&&LogText->State) TextFileClose(LogText);  close file if open*/
  if (LogText) TextFileKill(LogText); /* delete structure */
  LogText = NULL;
  doLog = 0;
  MenuMarkLogger (1); /* reset menu item label */
} /* end LoggerCancel */
