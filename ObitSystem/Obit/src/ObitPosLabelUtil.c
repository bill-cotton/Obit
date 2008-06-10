/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2006                                               */
/*;  Associated Universities, Inc. Washington DC, USA.                */
/*;                                                                   */
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
/*; Correspondence about this software should be addressed as follows:*/
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/

/* position labeling utilities for ObitView */
#include "ObitPosLabelUtil.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitPosLabelUtil.c
 * ObitPosLabelUtil -  position labeling utilities.
 * Liberally adopted from XFITSview/ObitView
 */

/*----------------------Public functions---------------------------*/

/**
 * Create appropriate character string describing a location
 * \param  pos   location on axis, angles in degrees  
 * \param  axis  label for the axis 
 * \param  label output string, must be large enough for 17 char
 *               contains first 4 characters of axis and the value.
 */
void ObitPosLableUtilAxisLabel(odouble pos, gchar* axis, gchar* label)
{
  odouble temp;
  olong h, m, itemp, i;
  gboolean recog, toHours, isStokes, iStokes;
  ofloat s;
  gchar rectypes[9][5]={"RA  ", "GLON", "ELON", "LL  ",
			"DEC ", "GLAT", "ELAT", "MM  ", "STOK"};
  gchar stoktypes[4][5]={"Ipol", "Qpol", "Upol", "Vpol"};
  gchar minus[2]="-";

  /* convert to hours for "RA" or "LL" only */
  toHours = (!strncmp(axis, rectypes[0], 2))  || 
    (!strncmp(axis, rectypes[3], 2));

  /* is this a Stokes Axis */
  isStokes = !strncmp(axis, rectypes[8], 4);
  
  /* make label (strip minus signs) */
  for (i=0;i<4;i++) label[i] = ' '; label[4] = 0;
  for (i=0;i<4;i++) {if (axis[i]==minus[0]) break; label[i]=axis[i];}
 recog = FALSE;  /* look for recognized position types */
 for (i=0;i<8;i++) recog = recog || (!strcmp (label, rectypes[i]));
 if(recog) { /* conversion for position */
   temp = pos; 
   if (temp>360.0) temp = temp - 360.0;
   if (toHours) 
     {if (temp<0.0) temp = temp + 360.0;
     temp = temp / 15.0; }
   itemp = (olong)temp; h = itemp;
   temp = temp - (odouble)itemp;
   if (temp<0.0) temp = - temp;
   temp = temp*60.0; itemp = (int)temp; m = itemp;
   temp = temp - (odouble)itemp;
   temp = temp*60.0; s = (ofloat)temp;}
 if (recog)
   {if (toHours)  /* display in hours */
     sprintf (&label[4], "%2.2d %2.2d %7.4f", h, m, s);
   else {        /* display in degrees */
     if (h <0) h = -h;
     sprintf (&label[4], "%3.2d %2.2d %6.3f", h, m, s);
     if (pos<0.0) label[4]='-'; /* neg declination */
   }
   }
 else if (isStokes)
   {iStokes = (int)(pos + 0.5) - 1;
   iStokes = MAX (MIN (iStokes, 3), 0);
   sprintf (&label[4], "     %s", stoktypes[iStokes]);}
 else  /* random type */
   sprintf (&label[4], "%13.6g", pos);
} /* End of ObitPosLabelUtilAxisLabel */

/** 
 * Convert RA in degrees to hh mm ss.sss (13 chars)
 * Gives simple display if axislabel not recognized
 * \param ra    RA in degrees 
 * \param rach  label for RA axis  
 * \param rast  [out] Ra string
 */
void ObitPosLabelUtilRA2HMS (odouble ra, gchar* rach, gchar* rast)
{
  odouble temp;
  olong h, m, itemp, i;
  gboolean recog, toHours;
  ofloat s;
  gchar minus[2]="-";
  gchar rectypes[8][5]={"RA  ", "GLON", "ELON", "LL  ",
			"DEC ", "GLAT", "ELAT", "MM  "};
  
  /* convert to hours for "RA" or "LL" only */
  toHours = (!strncmp(rach, rectypes[0], 2))  || 
    (!strncmp(rach, rectypes[3], 2));
  
  temp = ra; 
  if (temp>360.0) temp = temp - 360.0;
  if (toHours) 
    {if (temp<0.0) temp = temp + 360.0;
    temp = temp / 15.0; }
  itemp = (int)temp; h = itemp;
  temp = temp - (odouble)itemp;
  if (temp<0.0) temp = - temp;
  temp = temp*60.0; itemp = (int)temp; m = itemp;
  temp = temp - (odouble)itemp;
  temp = temp*60.0; s = (ofloat)temp;
  
  /* make label - use output as temp  */
 for (i=0;i<4;i++) rast[i] = ' '; rast[4] = 0;
 for (i=0;i<4;i++) {if (rach[i]==minus[0]) break; rast[i]=rach[i];}
  recog = 0;  /* look for recognized types */
  for (i=0;i<8;i++) recog = recog || (!strcmp (rast, rectypes[i]));
  if (recog)
    {if (toHours)  /* display in hours */
      sprintf (rast, "%2.2d %2.2d %7.4f", h, m, s);
    else         /* display in degrees */
      sprintf (rast, "%3.2d %2.2d %6.3f", h, m, s);
    }
  else
    sprintf (rast, "%13.6g", ra);
} /* End of ObitPosLabelUtilRA2HMS */


/**
 * Convert dec in degrees to dd mm ss.sss (14 char)
 * Gives simple display if axislabel not recognized
 * \param dec dec in degrees 
 * \param decch  label for Declination axis 
 * \param decst  [out] dec string
 */
void ObitPosLabelUtilDec2DMS(odouble dec, gchar* decch, gchar* decst)
{
  odouble temp;
  olong d, m, itemp, i;
  gboolean recog;
  ofloat s;
  gchar sign[2]=" ", minus[2]="-";
  gchar rectypes[8][5]={"RA  ", "GLON", "ELON", "LL  ",
			"DEC ", "GLAT", "ELAT", "MM  "};
 
  temp = dec; if (temp<0.0) temp = -temp; itemp = (olong)temp; d = itemp;
  temp = temp - (odouble)itemp;
  temp = temp*60.0; itemp = (olong)temp; m = itemp;
  temp = temp - (odouble)itemp;
  temp = temp*60.0; s = (ofloat)temp;
  if (dec<0.0) sign[0]=minus[0];   /* sign */
  
  /* make label - use output as temp  */
  for (i=0;i<4;i++) decst[i] = ' '; decst[4] = 0;
 for (i=0;i<4;i++) {if (decch[i]==minus[0]) break; decst[i]=decch[i];}
 recog = 0;  /* look for recognized types */
  for (i=0;i<8;i++) recog = recog || (!strcmp (decst, rectypes[i]));
  if (recog)
    sprintf (decst, "%1s%2.2d %2.2d %6.3f", sign, d, m, s);
  else
    sprintf (decst, "%13.6g", dec);
} /* End of ObitPosLabelUtilDec2DMS */

/** 
 * Convert RA in degrees to hours min and seconds
 * \param  ra  RA in degrees 
 * \param  h   [out] hours
 * \param  m   [out] min
 * \param  s   [out] sec 
 */
void ObitPosLabelUtilRAHMS (odouble ra, olong *h, olong *m, ofloat *s)
{
  odouble temp;
  olong itemp;
  
  temp = ra; if (temp<0.0) temp = temp + 360.0;
  if (temp>360.0) temp = temp - 360.0;
  temp =temp / 15.0; itemp = (olong)temp; *h = itemp;
  temp = temp - (odouble)itemp;
  temp = temp*60.0; itemp = (int)temp; *m = itemp;
  temp = temp - (odouble)itemp;
  temp = temp*60.0; *s = (ofloat)temp;
} /* End of ObitPosLabelUtilRAHMS */

/** 
 * Convert dec in degrees to degrees, min, sec
 * \param  dec Declination in degrees 
 * \param  d   [out] degrees
 * \param  m   [out] min
 * \param  s   [out] sec 
 */
void ObitPosLabelUtilDecDMS (odouble dec, olong *d, olong *m, ofloat *s)
{
  odouble temp;
  olong itemp;
  
  temp = dec; itemp = (olong)temp; *d = itemp;
  temp = temp - (odouble)itemp;
  if (temp<0.0) temp = -temp;
  temp = temp*60.0; itemp = (olong)temp; *m = itemp;
  temp = temp - (odouble)itemp;
  temp = temp*60.0; *s = (ofloat)temp;
} /* End of ObitPosLabelUtilDecDMS */

/** 
 * Convert RA in hours min and seconds to degrees
 * \param  h   hours
 * \param  m   min
 * \param  s   sec 
 * \param  ra  [out] RA in degrees 
 * \return 0 if in 0-24 hours else 1
 */
olong ObitPosLabelUtilHMSRA (gint h, olong m, ofloat s, odouble *ra)
{
  *ra = h + m/60.0 + s/3600.0;
  *ra = *ra * 15.0;
  if (*ra<0.0) return 1;
  if (*ra>360.0) return 1;
  return 0;
} /* End of ObitPosLabelUtilHMSRA */

/** 
 * convert dec in degrees, min, sec  to degrees 
 * \param  d   degrees
 * \param  m   min
 * \param  s   sec 
 * \param  ra  [out] RA in degrees 
 * \return 0 if in range +/-360 else 1 
 */
olong ObitPosLabelUtilDMSDec (gint d, olong m, ofloat s, odouble *dec)
{
  olong absdec = d;
  
  if (absdec<0) absdec = -absdec;
  *dec = absdec + m/60.0 + s/3600.0;
  if (d<0) *dec = -(*dec);
  if (*dec<-360.0) return 1;
  if (*dec>360.0) return 1;
  return 0;
} /* End of ObitPosLabelUtilDMSDec */

