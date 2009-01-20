/* $Id$  */
/* position labeling utilities for ObitView */
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
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "poslabel.h"
#define min(a,b) (((a) < (b)) ? (a) : (b))
#define max(a,b) (((a) > (b)) ? (a) : (b))
/**
 * Create appropriate character string describing a location
 * \param  pos   location on axis, angles in degrees  
 * \param  axis  label for the axis 
 * \param  label output string, must be large enough for 17 char
 *               contains first 4 characters of axis and the value.
 */
void AxisLabel(odouble pos, gchar* axis, gchar* label)
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
   iStokes = max (min (iStokes, 3), 0);
   sprintf (&label[4], "     %s", stoktypes[iStokes]);}
 else  /* random type */
   sprintf (&label[4], "%13.6g", pos);
} /* End of AxisLabel */

/**
 * Convert RA in degrees to hh mm ss.sss (17 chars)
 * Gives simple display if axislabel not recognized
 * \param ra    RA in degrees 
 * \param rach  label for RA axis  
 * \param rast  Ra string
 */
void ra2hms(odouble ra, gchar* rach, gchar* rast)
{
  odouble temp;
  olong h, m, itemp, i;
  gboolean recog, toHours;
  ofloat s;
  gchar rectypes[8][5]={"RA  ", "GLON", "ELON", "LL  ",
			"DEC ", "GLAT", "ELAT", "MM  "};
 gchar minus[2]="-";

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

 /* make label (strip minus signs) */
 for (i=0;i<4;i++) rast[i] = ' '; rast[4] = 0;
 for (i=0;i<4;i++) {if (rach[i]==minus[0]) break; rast[i]=rach[i];}
 recog = 0;  /* look for recognized types */
 for (i=0;i<8;i++) recog = recog || (!strcmp (rast, rectypes[i]));
 if (recog)
   {if (toHours)  /* display in hours */
     sprintf (&rast[4], "%2.2d %2.2d %6.4f", h, m, s);
   else         /* display in degrees */
     sprintf (&rast[4], "%3.2d %2.2d %6.3f", h, m, s);
   }
 else
   sprintf (&rast[4], "%13.6g", ra);
} /* End of ra2hms */


/**
 * Convert dec in degrees to Dec dd mm ss.sss (13 char)
 * Gives simple display if axislabel not recognized
 * \param dec dec in degrees 
 * \param decch  label for Declination axis 
 * \param decst  dec string
 */
void dec2dms(odouble dec, gchar* decch, gchar* decst)
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

 /* make label (strip minus signs) */
 for (i=0;i<4;i++) decst[i] = ' '; decst[4] = 0;
 for (i=0;i<4;i++) {if (decch[i]==minus[0]) break; decst[i]=decch[i];}
 recog = 0;  /* look for recognized types */
 for (i=0;i<8;i++) recog = recog || (!strcmp (decst, rectypes[i]));
 if (recog)
   sprintf (&decst[4], "%1s%2.2d %2.2d %6.3f", sign, d, m, s);
 else
   sprintf (&decst[4], "%13.6g", dec);
} /* End of dec2dms */

void rahms(double ra, int *h, int *m, float *s)
/* convert RA in degrees to hours min and seconds */
{double temp;
 int itemp;
 temp = ra; if (temp<0.0) temp = temp + 360.0;
 if (temp>360.0) temp = temp - 360.0;
 temp =temp / 15.0; itemp = (int)temp; *h = itemp;
 temp = temp - (double)itemp;
 temp = temp*60.0; itemp = (int)temp; *m = itemp;
 temp = temp - (double)itemp;
 temp = temp*60.0; *s = (float)temp;
} /* End of rahms */

void decdms(double dec, int *d, int *m, float *s)
/* convert dec in degrees to degrees, min, sec */
{double temp;
 int itemp;
 temp = dec; itemp = (int)temp; *d = itemp;
 temp = temp - (double)itemp;
 if (temp<0.0) temp = -temp;
 temp = temp*60.0; itemp = (int)temp; *m = itemp;
 temp = temp - (double)itemp;
 temp = temp*60.0; *s = (float)temp;
} /* End of decdms */

int hmsra(int h, int m, float s, double *ra)
/* convert RA in hours min and seconds to degrees*/
/* returns 0 if in 0-24 hours else 1 */
{
 *ra = h + m/60.0 + s/3600.0;
 *ra = *ra * 15.0;
 if (*ra<0.0) return 1;
 if (*ra>360.0) return 1;
 return 0;
 } /* End of hmsra */

int dmsdec(int d, int m, float s, double *dec)
/* convert dec in degrees, min, sec  to degrees*/
/* Note: this is also used for longitudes */
/* returns 0 if in range +/-360 else 1 */
{
 int absdec = d;
 if (absdec<0) absdec = -absdec;
 *dec = absdec + m/60.0 + s/3600.0;
 if (d<0) *dec = -(*dec);
 if (*dec<-360.0) return 1;
 if (*dec>360.0) return 1;
 return 0;
} /* End of dmsdec */

