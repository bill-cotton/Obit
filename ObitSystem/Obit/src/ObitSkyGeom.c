/* $Id$      */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2004-2014                                          */
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
/*;Correspondence about this software should be addressed as follows: */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
#include "ObitSkyGeom.h"
#include <math.h>
#include <string.h>
#include <stdio.h>
/*#include "precess.h"*/
/*#include "dsssubs.h"*/
/*#include "zsubs.h"*/


/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitSkyGeom.c
 * ObitSkyGeom function definitions.
 *
 * This file contains celestial coordinate utilities
 */

/*----------------- Macroes ---------------------------*/
#define RAGP  192.25*DG2RAD /* Ra (1950) of north galactic pole(rad) */
#define DECGP  27.40*DG2RAD /* Dec (1950)of north galactic pole (rad)*/
#define LONCP 123.00*DG2RAD /* Longitude of North Celestial pole (rad) */

/*---------------Private function prototypes----------------*/

/*----------------------Public functions---------------------------*/
/**
 * Determine the shift in RA and Dec between two celestial positions.
 * The shift is in (possibly) rotated coordinates.
 * Adopted from the AIPSish SHFTXY.FOR
 * \param ra       Initial Right Ascension in deg.
 * \param dec      Initial declination in deg.
 * \param rotate   Rotation of field, to E from N, deg.
 * \param shiftRA  Shifted Right Ascension in deg.
 * \param shiftDec Shifted declination in deg.
 * \param xShift   (out) Shift from ra to shiftRA in deg.
 * \param yShift   (out) Shift from dec to shiftDec in deg.
 */
void ObitSkyGeomShiftXY (odouble ra, odouble dec, ofloat rotate,
			odouble shiftRA, odouble shiftDec,
			ofloat *xShift, ofloat *yShift)
{
  odouble xxshft, yyshft;
  ofloat maprr;

  maprr = rotate * DG2RAD;

  /* Simple linear shift */
  xxshft = cos (DG2RAD*dec) * (shiftRA - ra);
  yyshft = shiftDec - dec;

  /*  Undo rotation */
  *xShift = (cos (maprr) * xxshft + sin (maprr) * yyshft);
  *yShift = (-sin (maprr) * xxshft + cos (maprr) * yyshft);
} /* end ObitSkyGeomShiftXY */

/**
 * Determine result of a shift applied to a celestial position.
 * The shift is in (possibly) rotated coordinates.
 * Adopted from the AIPSish XYSHFT.FOR
 * \param ra       Initial Right Ascension in deg.
 * \param dec      Initial declination in deg.
 * \param xShift   Shift from ra to shiftRA in deg.
 * \param yShift   Shift from dec to shiftDec in deg.
 * \param rotate   Rotation of field, to E from N, deg.
 * \param shiftRA  (out) Shifted Right Ascension in deg.
 * \param shiftDec (out) Shifted declination in deg.
 */
void ObitSkyGeomXYShift (odouble ra, odouble dec, 
			ofloat xShift, ofloat yShift, ofloat rotate,
			odouble *shiftRA, odouble *shiftDec)
{
  odouble cosDec, xxshft, yyshft;
  ofloat maprr;

  maprr = rotate * DG2RAD;
 
 /*  Undo rotation */
  xxshft = (cos (maprr) * xShift - sin (maprr) * yShift);
  yyshft = (sin (maprr) * xShift + cos (maprr) * yShift);

  /* Simple linear shift */
  *shiftDec = dec + yyshft;
  cosDec = cos (DG2RAD * dec);
  if (fabs(cosDec)<0.001) cosDec = 1.0;  /* trap pole */
  if (cosDec != 0.0) 
    *shiftRA = ra + xxshft / cosDec;
  else
    *shiftRA = ra;

} /* end ObitSkyGeomXYShift */

/**
 * Determine the shift in coordinate reference pixel between 
 * two celestial positions.
 *  The shift is in (possibly) rotated coordinates.
 * Adopted from the AIPSish SHFCRP.FOR
 * \param type     Projection type ("-SIN", "-NCP", otherwise linear, blank = -SIN)
 * \param ra       Initial (reference) Right Ascension in deg.
 * \param dec      Initial (reference) declination in deg.
 * \param rotate   Rotation of field, to E from N, deg.
 * \param xra      Shifted Right Ascension in deg.
 * \param ydec     Shifted declination in deg.
 * \param xshift   (out) Shift from ra to xra in deg.
 * \param yshift   (out) Shift from dec to xDec in deg.
 */
void  
ObitSkyGeomShiftCRP (gchar *type, odouble ra, odouble dec, ofloat rotate,
		     odouble xra, double xdec, 
		     ofloat *xshift, ofloat *yshift)
{
  odouble xxshft, yyshft;
  ofloat maprr;

  maprr = rotate * DG2RAD;
  /* L: SIN, NCP projection*/
  xxshft = cos (DG2RAD*xdec) * sin (DG2RAD*(xra-ra));
  
  /*  M: NCP PROJECTION */
  if (!strncmp(type,"-NCP",4)) {
    yyshft = (cos (DG2RAD*dec) - cos (DG2RAD*xdec) *
	      cos (DG2RAD*(xra-ra))) / sin (DG2RAD*dec);
  /*  M: SIN PROJECTION */
  } else if ((!strncmp(type,"-SIN",4)) ||(!strncmp(type,"    ",4))) {
    yyshft = sin (DG2RAD*xdec) * cos (DG2RAD*dec) - cos
      (DG2RAD*xdec) * sin (DG2RAD*dec) * cos (DG2RAD*(xra-ra));
  } else { /* Something else - do simple "linear" terms */
    xxshft = cos (DG2RAD*dec) * (xra - ra);
    yyshft = xdec - dec;
  }
  /* Undo rotation */
  *xshift = ( cos (maprr) * xxshft + sin (maprr) * yyshft) * RAD2DG;
  *yshift = (-sin (maprr) * xxshft + cos (maprr) * yyshft) * RAD2DG;
} /* end ObitSkyGeomShiftCRP */

/**
 * Finds coordinate shift from RA,DEC to XRA, XDEC and phase terms for
 *  -SIN projection.
 * Adapted from the AIPSish SHISIN.FOR.
 * \param  ra     Initial RA in degrees: reference position
 * \param  dec    Initial Declination in degrees
 * \param  rotate Image rotation in degrees
 * \param  xra    RA of shifted point in degrees
 * \param  xdec   Declination of shifted point in degrees
 * \param  dxyzc  (out) Phase term for Position offsets in x,y,z (2pi turns)
 */
void  
ObitSkyGeomShiftSIN (odouble ra, odouble dec, ofloat rotate,
		    odouble xra, double xdec, ofloat dxyzc[3])
{
  ofloat maprr, xshift, yshift;
  odouble xxshft, yyshft, dzctmp;

  maprr = rotate * DG2RAD;
  
  /*  l, m */
  xxshft = cos (DG2RAD*xdec) * sin (DG2RAD*(xra-ra));
  yyshft = sin (DG2RAD*xdec) * cos (DG2RAD*dec) -
    cos (DG2RAD*xdec) * sin (DG2RAD*dec) * cos (DG2RAD*(xra-ra));
  
  /*  undo rotation */
  xshift = ( cos (maprr) * xxshft + sin (maprr) * yyshft);
  yshift = (-sin (maprr) * xxshft + cos (maprr) * yyshft);
  
  /* Phase shift parameters */
  dxyzc[0] = xshift;
  dxyzc[1] = yshift;
  dzctmp = 1.0 - (dxyzc[0]*dxyzc[0]) - (dxyzc[1]*dxyzc[1]);
  dzctmp = sqrt (dzctmp);
  
  /* Prepare phase calc; mult by 2pi */
  dxyzc[0] = 2.0*G_PI * dxyzc[0];
  dxyzc[1] = 2.0*G_PI * dxyzc[1];
  dxyzc[2] = 2.0*G_PI *(dzctmp - 1.0);
  
} /* end  ObitSkyGeomShiftSIN */

/**
 * Finds coordinate shift from RA,DEC to XRA, XDEC and phase terms for
 *  -NCP projection.
 * Adapted from the AIPSish SHINCP.FOR.
 * \param  ra     Initial RA in degrees: reference position
 * \param  dec    Initial Declination in degrees
 * \param  rotate Image rotation in degrees
 * \param  xra    RA of shifted point in degrees
 * \param  xdec   Declination of shifted point in degrees
 * \param  dxyzc  (out) Phase term for Position offsets in x,y,z
 */
void  
ObitSkyGeomShiftNCP (odouble ra, odouble dec, ofloat rotate,
		    odouble xra, double xdec, ofloat dxyzc[3])
{
  ofloat maprr, xshift, yshift;
  odouble xxshft, yyshft;

  maprr = rotate * DG2RAD;

  /* L, M */
  xxshft = cos (DG2RAD*xdec) * sin (DG2RAD*(xra-ra));
  yyshft = sin (DG2RAD* dec);
  if (yyshft != 0.0) {
    yyshft = (cos (DG2RAD*dec) - cos (DG2RAD*xdec) *
	      cos (DG2RAD*(xra-ra))) / yyshft;
  } else {
    yyshft = xdec - dec;
  }
    
    /*  Undo rotation */
    xshift = ( cos (maprr) * xxshft + sin (maprr) * yyshft);
    yshift = (-sin (maprr) * xxshft + cos (maprr) * yyshft);
    /* Convert to radians */
    dxyzc[0] = xshift;
    dxyzc[1] = yshift;
    
    /* Prepare phase calc; mult by 2pi */
    dxyzc[0] = 2.0*G_PI * dxyzc[0];
    dxyzc[1] = 2.0*G_PI * dxyzc[1];
    dxyzc[2] = 0.0;
} /* end  ObitSkyGeomShiftNCP */

/**
 * Determines the coordinates (raout,decout) corresponding to a 
 * displacement (l,m) given by the direction cosines from coordinates 
 * (ra0,dec0).  the direction cosine l is assumed to be positive to 
 * the east; m is positive to the north.  the routine works for 
 * 4 kinds of projective geometries and for celestial, ecliptic, or 
 * galactic coordinate systems. 
 * this subroutine always uses an accurate computation. 
 * All angles in this subroutine are in radians. 
 * Adapted from the AIPSish NEWPOS.FOR.
 *  \param  Proj   Projection type (Aitoff and  Mercator not supported)
 *  \param  ra0    coordinate reference right ascension (longitude) 
 *  \param  dec0   coordinate reference declination (latitude) 
 *  \param  l      cosine angle of displacement to east 
 *  \param  m      cosine angle of displacement to north 
 *  \param  raout  [out] right ascension or longitude at (l,m) 
 *  \param  decout [out] declination or latitude at (l,m) 
 *  \param  ierr   [out] error condition: 0 = ok, 1 = l,m crazy, 
 *                 2 = bad type,  3 = answer undefined 
 */
void 
ObitSkyGeomNewPos (ObitSkyGeomProj Proj, odouble ra0, odouble dec0, 
		   odouble l, odouble m, odouble *raout, odouble *decout, 
		   olong *ierr)
{
  odouble sins, coss, dect, rat, dt, mg, da, dz, cos0, sin0;
  odouble twopi=6.28318530717959, deps=1.0e-5;
  gboolean OK;
  
  /* default output */
  *decout = 0.00;
  *raout  = 0.00;
  
  /* Check legitimate projection */
  *ierr = 2;
  if ((Proj<OBIT_SkyGeom_SIN) || (Proj>OBIT_SkyGeom_STG)) return;

  /* Check if angles in range */
  *ierr = 1;
  sins = l*l + m*m;
  OK = (sins <= 1.0);
  OK = OK || (Proj==OBIT_SkyGeom_STG);
  OK = OK || (Proj==OBIT_SkyGeom_AIT);
  OK = OK || (Proj==OBIT_SkyGeom_GLS);
  OK = OK || ((Proj==OBIT_SkyGeom_NCP) && (sins < twopi*twopi/4.0));
  if (!OK) return;

  *ierr = 3;
  cos0 = cos(dec0);
  sin0 = sin(dec0);

  /* Branch by type */
  switch (Proj) {
  case OBIT_SkyGeom_SIN:
    /* sin projection */
    coss = sqrt (1.0 - sins);
    dt = sin0 * coss + cos0 * m;
    if (fabs(dt) > 1.0) return;
    dect = asin (dt);
    rat = cos0 * coss - sin0 * m;
    if ((rat == 0.0)  &&  (l == 0.0)) return;
    rat = atan2 (l, rat) + ra0;
    break;
  case OBIT_SkyGeom_TAN:
    /* tan projection */
    dect = cos0 - m * sin0;
    if (dect == 0.0) return;
    rat = ra0 + atan2 (l, dect);
    dect = atan (cos(rat-ra0) * (m * cos0 + sin0) / dect);
    break;
  case OBIT_SkyGeom_ARC:
    /* arc projection */
    sins = sqrt (sins);
    coss = cos (sins);
    if (sins != 0.0) {
      sins = sin (sins) / sins;
    } else {
      sins = 1.0;
    } 
    dt = m * cos0 * sins + sin0 * coss;
    if (fabs(dt) > 1.0) return;
    dect = asin (dt);
    da = coss - dt * sin0;
    dt = l * sins * cos0;
    if ((da == 0.0)  &&  (dt == 0.0)) return;
    rat = ra0 + atan2 (dt, da);
    break;
  case OBIT_SkyGeom_NCP:
    /* WSRT (north pole projection) */
    dect = cos0 - m * sin0;
    if (dect == 0.0) return;
    rat = ra0 + atan2 (l, dect);
    dt = cos (rat-ra0);
    if (dt == 0.0) return;
    dect = dect / dt;
    if (fabs(dect) > 1.0) return;
    dect = acos (MIN(dect, 1.000));
    if (dec0 < 0.0) dect = -dect;
    break;
  case OBIT_SkyGeom_GLS:
    /* global sinusoid */
    dect = dec0 + m;
    if (fabs(dect) > twopi/4.0) return;
    coss = cos (dect);
    if (fabs(l) > twopi*coss/2.0) return;
    rat = ra0;
    if (coss > deps) rat = rat + l / coss;
    break;
  case OBIT_SkyGeom_MER: /* Not supported due to extra parameters */
    /* Mercator */
    *ierr = 2;
    return;
    /*
      rat = l / geomd1(locnum) + ra0;
      if (fabs(rat-ra0) > twopi) return;
      dt = 0.0;
      if (geomd2(locnum) != 0.0) dt = (m + geomd3(locnum)) / geomd2(locnum);
      dt = exp (dt);
      dect = 2.0 * atan (dt) - twopi / 4.0;
    */
    break;
  case OBIT_SkyGeom_AIT:  /* Not supported due to extra parameters */
    /* Aitoff */
    *ierr = 2;
    return;
    /*    rat = ra0;
	  dect = dec0;
	  if ((l != 0.0)  ||  (m != 0.0)) {
	  dz = 4.0 - l*l/(4.0*geomd1(locnum)*geomd1(locnum)) -
	  ((m+geomd3(locnum))/geomd2(locnum))**2;
	  if ((dz > 4.0)  ||  (dz < 2.0)) return;
	  dz = 0.50 * sqrt (dz);
	  dd = (m+geomd3(locnum)) * dz / geomd2(locnum);
	  if (fabs(dd-) > 1.0) return;
	  dd = asin (dd);
	  if (fabs(cos (dd-)) < deps) return;
	  da = l * dz / (2.0 * geomd1(locnum) * cos[dd-1]);
	  if (fabs(da) > 1.0) return;
	  da = asin (da);
	  rat = ra0 + 2.0 * da;
	  dect = dd;
	  if (fabs(dect) > twopi/4.0) return;
	  if (fabs(da) > twopi/4.0) return;
	  }
    */
    break;
  case OBIT_SkyGeom_STG:
    /* stereographic */
    dz = (4.0 - sins) / (4.0 + sins);
    if (fabs(dz) > 1.0) return;
    dect = dz * sin0 + m * cos0 * (1.0+dz) / 2.0;
    if (fabs(dect) > 1.0) return;
    dect = asin (dect);
    rat = cos(dect);
    if (fabs(rat) < deps) return;
    rat = l * (1.0+dz) / (2.0 * rat);
    if (fabs(rat) > 1.0) return;
    rat = asin (rat);
    mg = 1.0 + sin(dect) * sin0 + cos(dect) * cos0 * cos(rat);
    if (fabs(mg) < deps) return;
    mg = 2.0 * (sin(dect) * cos0 - cos(dect) * sin0 * cos(rat)) / mg;
    if (fabs(mg-m) > deps) rat = twopi/2.0 - rat;
    rat = ra0 + rat;
    break;
  default: /* Unknown */
    return;
    break;
  }; /* end switch by projection type */

  /* return: in range ra */
  *raout = rat;
  *decout = dect;
  *ierr = 0;
  if (*raout-ra0 > twopi/2.0)  *raout -= twopi;
  if (*raout-ra0 < -twopi/2.0) *raout += twopi;
} /* end of routine ObitSkyGeomNewPos */ 

/**
 * Determine accurate position for pixel coordinates from IRAF
 * style CD matrix.  
 * Note: xinc, yinc, and rot can be derived from cd1 and cd2 and 
 * should be compatible with them.
 * Taken from FITSview family
 * \param  xpix    x pixel number  (RA or long without rotation)
 * \param  ypix    y pixel number  (dec or lat without rotation)
 * \param  xref    x reference coordinate value (deg)
 * \param  yref    y reference coordinate value (deg)
 * \param  xrefpix x reference pixel
 * \param  yrefpix y reference pixel
 * \param  xinc    x coordinate increment (deg)
 * \param  yinc    y coordinate increment (deg)
 * \param  rot     rotation (deg)  (from N through E)
 * \param  type    projection type code e.g. "-SIN", blank = -SIN
 *                 Does: -SIN, -TAN, -ARC, -NCP, -GLS, -MER, -AIT 
 *                 projections anything else is linear 
 * \param  cd1     first column of CD matrix
 * \param  cd2     second column of CD matrix
 * \param  xpos    [out] x (RA) coordinate (deg)
 * \param  ypos    [out]y (dec) coordinate (deg)
 * \return  0 if successful otherwise: 1 = angle too large for projection;
 */
olong 
ObitSkyGeomCDpos(ofloat xpix, ofloat ypix, odouble xref, odouble yref,
		 ofloat xrefpix, ofloat yrefpix, ofloat xinc, ofloat yinc, 
		 ofloat rot, ofloat cd1[2], ofloat cd2[2], gchar *type, 
		 odouble *xpos, odouble *ypos)
 {
   odouble dx, dy, l, m;

  /*   Offset from ref pixel  */
  dx = (xpix-xrefpix);
  dy = (ypix-yrefpix);

  /* convert to l and m  */
  l = cd1[0]*dx + cd1[1]*dy;
  m = cd2[0]*dx + cd2[1]*dy;

  /* determine position */
  return ObitSkyGeomWorldPosLM (l, m, xref, yref, xinc, yinc, rot, 
				type, xpos, ypos);
}  /* End of ObitSkyGeomCDpos */

/**
 * Determine accurate position for pixel coordinates.
 * Taken from FITSview family
 * \param  xpix    x pixel number  (RA or long without rotation)
 * \param  ypix    y pixel number  (dec or lat without rotation)
 * \param  xref    x reference coordinate value (deg)
 * \param  yref    y reference coordinate value (deg)
 * \param  xrefpix x reference pixel
 * \param  yrefpix y reference pixel
 * \param  xinc    x coordinate increment (deg)
 * \param  yinc    y coordinate increment (deg)
 * \param  rot     rotation (deg)  (from N through E)
 * \param  type    projection type code e.g. "-SIN", blank = -SIN
 *                 Does: -SIN, -TAN, -ARC, -NCP, -GLS, -MER, -AIT 
 *                 projections anything else is linear 
 * \param  xpos    [out] x (RA) coordinate (deg)
 * \param  ypos    [out]y (dec) coordinate (deg)
 * \return  0 if successful otherwise: 1 = angle too large for projection;
 */
olong 
ObitSkyGeomWorldPos(ofloat xpix, ofloat ypix, odouble xref, odouble yref, 
		    ofloat xrefpix, ofloat yrefpix, ofloat xinc, ofloat yinc, 
		    ofloat rot, gchar *type, odouble *xpos, odouble *ypos)
{
  odouble cosr, sinr, dx, dy, temp;
  odouble cond2r=1.745329252e-2;
  
  /*   Offset from ref pixel  */
  dx = (xpix-xrefpix) * xinc;
  dy = (ypix-yrefpix) * yinc;

  /*   Take out rotation  */
  cosr = cos(rot*cond2r);
  sinr = sin(rot*cond2r);
  if (rot!=0.0)
    {temp = dx * cosr - dy * sinr;
    dy = dy * cosr + dx * sinr;
    dx = temp;}

  /* determine position */
  return ObitSkyGeomWorldPosLM (dx, dy, xref, yref, xinc, yinc, rot, type, 
				xpos, ypos);
}  /* End of ObitSkyGeomWorldPos */

/**
 * Determine accurate position for pixel coordinates from offsets 
 * from the reference position.
 * Taken from FITSview family
 * \param  dx      x coordinate offset  (RA or long) 
 * \param  dy      y coordinate offset  (dec or lat)
 * \param  xref    x reference coordinate value (deg)
 * \param  yref    y reference coordinate value (deg
 * \param  xinc    x coordinate increment (deg)
 * \param  yinc    y coordinate increment (deg)
 * \param  rot     rotation (deg)  (from N through E)
 * \param  type    projection type code e.g. "-SIN", blank = -SIN
 *                 Does: -SIN, -TAN, -ARC, -NCP, -GLS, -MER, -AIT 
 *                 projections anything else is linear 
 * \param  xpos    [out] x (RA) coordinate (deg)
 * \param  ypos    [out] y (dec) coordinate (deg)
 * \return  0 if successful otherwise: 1 = angle too large for projection;
 */
olong 
ObitSkyGeomWorldPosLM(odouble dx, odouble dy, odouble xref, odouble yref, 
		      ofloat xinc, ofloat yinc, ofloat rot, gchar *type, 
		      odouble *xpos, odouble *ypos)
{
  odouble cosr, sinr, dz;
  odouble sins, coss, dect=0.0, rat=0.0, dt, l, m, mg, da, dd, cos0, sin0;
  odouble dec0, ra0, decout, raout;
  odouble geo1, geo2, geo3;
  odouble cond2r=1.745329252e-2;
  odouble twopi = 6.28318530717959, deps = 1.0e-5;
  olong   i, itype;
  gchar ctypes[8][5] ={"-SIN","-TAN","-ARC","-NCP", "-GLS", "-MER", "-AIT",
		       "-STG"};

  /*   rotation  */
  cosr = cos(rot*cond2r);
  sinr = sin(rot*cond2r);

  /*  find type  */
  itype = 0;  /* default type is linear */
  for (i=0;i<8;i++) if (!strncmp(type, ctypes[i], 4)) itype = i+1;
  /* Trap blank = -SIN = EVLA/AIPS++ screwup */
  if (!strncmp(type, "    ", 4)) itype = 1;
  
  /* default, linear result for error return  */
  *xpos = xref + dx;
  *ypos = yref + dy;

  /* convert to radians  */
  ra0 = xref * cond2r;
  dec0 = yref * cond2r;
  l = dx * cond2r;
  m = dy * cond2r;
  sins = l*l + m*m;
  decout = 0.0;
  raout = 0.0;
  cos0 = cos(dec0);
  sin0 = sin(dec0);

  /* process by case  */
  switch (itype) {
  case 0:   /* linear */
    rat =  ra0 + l;
    dect = dec0 + m;
    break;
  case 1:   /* -SIN sin*/ 
    if (sins>1.0) return 1;
    coss = sqrt (1.0 - sins);
    dt = sin0 * coss + cos0 * m;
    if ((dt>1.0) || (dt<-1.0)) return 1;
    dect = asin (dt);
    rat = cos0 * coss - sin0 * m;
    if ((rat==0.0) && (l==0.0)) return 1;
    rat = atan2 (l, rat) + ra0;
    break;
  case 2:   /* -TAN tan */
    if (sins>1.0) return 1;
    dect = cos0 - m * sin0;
    if (dect==0.0) return 1;
    rat = ra0 + atan2 (l, dect);
    dect = atan (cos(rat-ra0) * (m * cos0 + sin0) / dect);
    break;
  case 3:   /* -ARC Arc*/
    if (sins>=twopi*twopi/4.0) return 1;
    sins = sqrt(sins);
    coss = cos (sins);
    if (sins!=0.0) sins = sin (sins) / sins;
    else
      sins = 1.0;
    dt = m * cos0 * sins + sin0 * coss;
    if ((dt>1.0) || (dt<-1.0)) return 1;
    dect = asin (dt);
    da = coss - dt * sin0;
    dt = l * sins * cos0;
    if ((da==0.0) && (dt==0.0)) return 1;
    rat = ra0 + atan2 (dt, da);
    break;
  case 4:   /* -NCP North celestial pole*/
    dect = cos0 - m * sin0;
    if (dect==0.0) return 1;
    rat = ra0 + atan2 (l, dect);
    dt = cos (rat-ra0);
    if (dt==0.0) return 1;
    dect = dect / dt;
    if ((dect>1.0) || (dect<-1.0)) return 1;
    dect = acos (MIN(dect, 1.000));
    if (dec0<0.0) dect = -dect;
    break;
  case 5:   /* -GLS global sinusoid */
    dect = dec0 + m;
    if (fabs(dect)>twopi/4.0) return 1;
    coss = cos (dect);
    if (fabs(l)>twopi*coss/2.0) return 1;
    rat = ra0;
    if (coss>deps) rat = rat + l / coss;
    break;
  case 6:   /* -MER mercator*/
    dt = yinc * cosr + xinc * sinr;
    if (dt==0.0) dt = 1.0;
    dy = (yref/2.0 + 45.0) * cond2r;
    dx = dy + dt / 2.0 * cond2r;
    dy = log (tan (dy));
    dx = log (tan (dx));
    geo2 = dt * cond2r / (dx - dy);
    geo3 = geo2 * dy;
    geo1 = cos (yref*cond2r);
    if (geo1<=0.0) geo1 = 1.0;
    rat = l / geo1 + ra0;
    if (fabs(rat - ra0) > twopi) return 1; /* added 10/13/94 DCW/EWG */
    dt = 0.0;
    if (geo2!=0.0) dt = (m + geo3) / geo2;
    dt = exp (dt);
    dect = 2.0 * atan (dt) - twopi / 4.0;
    break;
  case 7:   /* -AIT Aitoff*/
    dt = yinc*cosr + xinc*sinr;
    if (dt==0.0) dt = 1.0;
    dt = dt * cond2r;
    dy = yref * cond2r;
    dx = sin(dy+dt)/sqrt((1.0+cos(dy+dt))/2.0) -
      sin(dy)/sqrt((1.0+cos(dy))/2.0);
    if (dx==0.0) dx = 1.0;
    geo2 = dt / dx;
    dt = xinc*cosr - yinc* sinr;
    if (dt==0.0) dt = 1.0;
    dt = dt * cond2r;
    dx = 2.0 * cos(dy) * sin(dt/2.0);
    if (dx==0.0) dx = 1.0;
    geo1 = dt * sqrt((1.0+cos(dy)*cos(dt/2.0))/2.0) / dx;
    geo3 = geo2 * sin(dy) / sqrt((1.0+cos(dy))/2.0);
    rat = ra0;
    dect = dec0;
    if ((l==0.0) && (m==0.0)) break;
    dz = 4.0 - l*l/(4.0*geo1*geo1) - ((m+geo3)/geo2)*((m+geo3)/geo2) ;
    if ((dz>4.0) || (dz<2.0)) return 1;;
    dz = 0.5 * sqrt (dz);
    dd = (m+geo3) * dz / geo2;
    if (fabs(dd)>1.0) return 1;;
    dd = asin (dd);
    if (fabs(cos(dd))<deps) return 1;;
    da = l * dz / (2.0 * geo1 * cos(dd));
    if (fabs(da)>1.0) return 1;;
    da = asin (da);
    rat = ra0 + 2.0 * da;
    dect = dd;
    break;
  case 8:   /* -STG Sterographic*/
    dz = (4.0 - sins) / (4.0 + sins);
    if (fabs(dz)>1.0) return 1;
    dect = dz * sin0 + m * cos0 * (1.0+dz) / 2.0;
    if (fabs(dect)>1.0) return 1;
    dect = asin (dect);
    rat = cos(dect);
    if (fabs(rat)<deps) return 1;
    rat = l * (1.0+dz) / (2.0 * rat);
    if (fabs(rat)>1.0) return 1;
    rat = asin (rat);
    mg = 1.0 + sin(dect) * sin0 + cos(dect) * cos0 * cos(rat);
    if (fabs(mg)<deps) return 1;
    mg = 2.0 * (sin(dect) * cos0 - cos(dect) * sin0 * cos(rat)) / mg;
    if (fabs(mg-m)>deps) rat = twopi/2.0 - rat;
    rat = ra0 + rat;
    break;
  } /* end switch */
  
  /*  return ra in range  */
  raout = rat;
  decout = dect;
  if (raout-ra0>twopi/2.0) raout = raout - twopi;
  if (raout-ra0<-twopi/2.0) raout = raout + twopi;
  if (raout < 0.0) raout += twopi; /* added by DCW 10/12/94 */
  
  /*  correct units back to degrees  */
  *xpos  = raout  / cond2r;
  *ypos  = decout  / cond2r;
  return 0;
}  /* End of ObitSkyGeomWorldPosLM */

/**
 * Determine pixel for given coordinate
 * Taken from FITSview family
 * \param  xpos    x (RA) coordinate (deg)
 * \param  ypos    y (dec) coordinate (deg)
 * \param  xref    x reference coordinate value (deg)
 * \param  yref    y reference coordinate value (deg)
 * \param  xrefpix x reference pixel
 * \param  yrefpix y reference pixel
 * \param  xinc    x coordinate increment (deg)
 * \param  yinc    y coordinate increment (deg)
 * \param  rot     rotation (deg)  (from N through E)
 * \param  type    projection type code e.g. "-SIN", blank = -SIN
 *                 Does: -SIN, -TAN, -ARC, -NCP, -GLS, -MER, -AIT 
 *                 projections anything else is linear 
 * \param  xpix    [out] x pixel number  (RA or long without rotation)
 * \param  ypix    [out] y pixel number  (dec or lat without rotation)
 * \return  0 if successful otherwise: 1 = angle too large for projection
 *          2 = bad values 
 */
olong 
ObitSkyGeomXYpix(odouble xpos, odouble ypos, odouble xref, odouble yref, 
		 ofloat xrefpix, ofloat yrefpix, ofloat xinc, ofloat yinc, 
		 ofloat rot, gchar *type, ofloat *xpix, ofloat *ypix)
{
  odouble dx, dy, dz, sinr, cosr;
  olong    iret;
  odouble cond2r=1.745329252e-2;
  
  /* get coordinate offset */
  iret = ObitSkyGeomXYPixLM (xpos, ypos, xref, yref, xinc, yinc, rot, type, 
			     &dx, &dy);
  
  /*  Correct for rotation */
  cosr = cos(rot*cond2r);
  sinr = sin(rot*cond2r);
  dz = dx*cosr + dy*sinr;
  dy = dy*cosr - dx*sinr;
  dx = dz;
  
  /*     convert to pixels  */
  *xpix = dx / xinc + xrefpix;
  *ypix = dy / yinc + yrefpix;
  return iret;
}  /* end ObitSkyGeomXYpix */

/**
 * Determine accurate pixel coordinates for an RA and Dec uses IRAF  
 * style CD matrix. 
 * Note: xinc, yinc, and rot can be derived from cd1 and cd2 and 
 * should be compatible with them.
 * Taken from FITSview family
 * \param  xpos    x (RA) coordinate (deg)
 * \param  ypos    y (dec) coordinate (deg)
 * \param  xref    x reference coordinate value (deg)
 * \param  yref    y reference coordinate value (deg)
 * \param  xrefpix x reference pixel
 * \param  yrefpix y reference pixel
 * \param  xinc    x coordinate increment (deg)
 * \param  yinc    y coordinate increment (deg)
 * \param  rot     rotation (deg)  (from N through E)
 * \param  type    projection type code e.g. "-SIN", blank = -SIN
 *                 Does: -SIN, -TAN, -ARC, -NCP, -GLS, -MER, -AIT 
 *                 projections anything else is linear 
 * \param  cd1     first column of CD matrix
 * \param  cd2     second column of CD matrix
 * \param  xpix    [out] x pixel number  (RA or long without rotation)
 * \param  ypix    [out] y pixel number  (dec or lat without rotation)
 * \return  0 if successful otherwise: 1 = angle too large for projection
 *          2 = bad values 
 */
olong 
ObitSkyGeomCDpix(odouble xpos, odouble ypos, odouble xref, odouble yref, 
		 ofloat xrefpix, ofloat yrefpix, ofloat xinc, ofloat yinc, ofloat rot,
		 ofloat icd1[2], ofloat icd2[2], gchar *type, 
		 ofloat *xpix, ofloat *ypix)
{
  odouble l, m, dx, dy;
  olong    iret;
  
  /* get coordinate offset */
  iret = ObitSkyGeomXYPixLM (xpos, ypos, xref, yref, xinc, yinc, rot, type, 
			     &l, &m);
  
  /*  Correct by inverse CD matrix */
  dx = icd1[0]*l + icd1[1]*m;
  dy = icd2[0]*l + icd2[1]*m;
  
  /*     convert to pixels  */
  *xpix = dx + xrefpix;
  *ypix = dy + yrefpix;
  return iret;
}  /* end ObitSkyGeomCDpix */

/**
 * Determine accurate coordinate offsets for an RA and Dec 
 * Taken from FITSview family
 * \param  xpos    x (RA) coordinate (deg)
 * \param  ypos    y (dec) coordinate (deg)
 * \param  xref    x reference coordinate value (deg)
 * \param  yref    y reference coordinate value (deg)
 * \param  xrefpix x reference pixel
 * \param  yrefpix y reference pixel
 * \param  xinc    x coordinate increment (deg)
 * \param  yinc    y coordinate increment (deg)
 * \param  rot     rotation (deg)  (from N through E)
 * \param  type    projection type code e.g. "-SIN", blank = -SIN
 *                 Does: -SIN, -TAN, -ARC, -NCP, -GLS, -MER, -AIT 
 *                 projections anything else is linear 
 * \param  dx      [out] x projected offset  (RA or long without rotation)
 * \param  dy      [out] y projected offset  (dec or lat without rotation)
 * \return  0 if successful otherwise: 1 = angle too large for projection
 *          2 = bad values 
 */
olong 
ObitSkyGeomXYPixLM(odouble xpos, odouble ypos, odouble xref, odouble yref, 
		   ofloat xinc, ofloat yinc, ofloat rot, gchar *type, 
		   odouble *dx, odouble *dy)
{
  odouble dz, r, ra0, dec0, ra, dec, coss, sins, dt, da, dd, sint, ddx, ddy;
  odouble l, m=0.0, geo1, geo2, geo3, sinr, cosr;
  odouble cond2r=1.745329252e-2, deps=1.0e-5, twopi=6.28318530717959;
  olong   i, itype;
  gchar ctypes[8][5] ={"-SIN","-TAN","-ARC","-NCP", "-GLS", "-MER", "-AIT",
		       "-STG"};
  
  /* 0h wrap-around tests added by D.Wells 10/12/94: */
  dt = (xpos - xref);
  if (dt > +180) xpos -= 360;
  if (dt < -180) xpos += 360;
  /* NOTE: changing input argument xpos is OK (call-by-value in C!) */
  
  /* default values - linear */
  *dx = xpos - xref;
  *dy = ypos - yref;
  dz = 0.0;
  /*     check axis increments - bail out if either 0 */
  if ((xinc==0.0) || (yinc==0.0)) {*dx=0.0; *dy=0.0; return 2;}
  
  /*  find type  */
  itype = 0;  /* default type is linear */
  for (i=0;i<8;i++) if (!strncmp(type, ctypes[i], 4)) itype = i+1;
  /* Trap blank = -SIN = EVLA/AIPS++ screwup */
  if (!strncmp(type, "    ", 4)) itype = 1;
  if (itype==0) return 0;  /* done if linear */
  
  /*  rotation */
  r = rot * cond2r;
  cosr = cos (r);
  sinr = sin (r);
  
  /* Non linear position */
  ra0 = xref * cond2r;
  dec0 = yref * cond2r;
  ra = xpos * cond2r;
  dec = ypos * cond2r;
  
  /* compute direction cosine */
  coss = cos (dec);
  sins = sin (dec);
  l = sin(ra-ra0) * coss;
  sint = sins * sin(dec0) + coss * cos(dec0) * cos(ra-ra0);
  /* process by case  */
  switch (itype) {
  case 1:   /* -SIN sin*/ 
    if (sint<0.0) return 1;
    m = sins * cos(dec0) - coss * sin(dec0) * cos(ra-ra0);
    break;
  case 2:   /* -TAN tan */
    if (sint<=0.0) return 1;
    m = sins * sin(dec0) + coss * cos(dec0) * cos(ra-ra0);
    l = l / m;
    m = (sins * cos(dec0) - coss * sin(dec0) * cos(ra-ra0)) / m;
    break;
  case 3:   /* -ARC Arc*/
    m = sins * sin(dec0) + coss * cos(dec0) * cos(ra-ra0);
    if (m<-1.0) m = -1.0;
    if (m>1.0) m = 1.0;
    m = acos (MIN (m,1.000));
    if (m!=0) 
      m = m / sin(m);
    else
      m = 1.0;
    l = l * m;
    m = (sins * cos(dec0) - coss * sin(dec0) * cos(ra-ra0)) * m;
    break;
  case 4:   /* -NCP North celestial pole*/
    if (dec0==0.0) 
      return 1;  /* can't stand the equator */
    else
      m = (cos(dec0) - coss * cos(ra-ra0)) / sin(dec0);
    break;
  case 5:   /* -GLS global sinusoid */
    dt = ra - ra0;
    if (fabs(dec)>twopi/4.0) return 1;
    if (fabs(dec0)>twopi/4.0) return 1;
    m = dec - dec0;
    l = dt * coss;
    break;
  case 6:   /* -MER mercator*/
    dt = yinc * cosr + xinc * sinr;
    if (dt==0.0) dt = 1.0;
    ddy = (yref/2.0 + 45.0) * cond2r;
    ddx = ddy + dt / 2.0 * cond2r;
    ddy = log (tan (ddy));
    ddx = log (tan (ddx));
    geo2 = dt * cond2r / (ddx - ddy);
    geo3 = geo2 * ddy;
    geo1 = cos (yref*cond2r);
    if (geo1<=0.0) geo1 = 1.0;
    dt = ra - ra0;
    l = geo1 * dt;
    dt = dec / 2.0 + twopi / 8.0;
    dt = tan (dt);
    if (dt<deps) return 2;
    m = geo2 * log (dt) - geo3;
    break;
  case 7:   /* -AIT Aitoff*/
    l = 0.0;
    m = 0.0;
    da = (ra - ra0) / 2.0;
    if (fabs(da)>twopi/4.0) return 1;
    dt = yinc*cosr + xinc*sinr;
    if (dt==0.0) dt = 1.0;
    dt = dt * cond2r;
    ddy = yref * cond2r;
    ddx = sin(ddy+dt)/sqrt((1.0+cos(ddy+dt))/2.0) -
      sin(ddy)/sqrt((1.0+cos(ddy))/2.0);
    if (ddx==0.0) ddx = 1.0;
    geo2 = dt / ddx;
    dt = xinc*cosr - yinc* sinr;
    if (dt==0.0) dt = 1.0;
    dt = dt * cond2r;
    ddx = 2.0 * cos(ddy) * sin(dt/2.0);
    if (ddx==0.0) ddx = 1.0;
    geo1 = dt * sqrt((1.0+cos(ddy)*cos(dt/2.0))/2.0) / ddx;
    geo3 = geo2 * sin(ddy) / sqrt((1.0+cos(ddy))/2.0);
    dt = sqrt ((1.0 + cos(dec) * cos(da))/2.0);
    if (fabs(dt)<deps) return 3;
    l = 2.0 * geo1 * cos(dec) * sin(da) / dt;
    m = geo2 * sin(dec) / dt - geo3;
    break;
  case 8:   /* -STG Sterographic*/
    da = ra - ra0;
    if (fabs(dec)>twopi/4.0) return 1;
    dd = 1.0 + sins * sin(dec0) + coss * cos(dec0) * cos(da);
    if (fabs(dd)<deps) return 1;
    dd = 2.0 / dd;
    l = l * dd;
    m = dd * (sins * cos(dec0) - coss * sin(dec0) * cos(da));
    break;
  }  /* end of itype switch */
  
  /*   back to degrees  */
  *dx = l / cond2r;
  *dy = m / cond2r;
  return 0;
}  /* end ObitSkyGeomXYPixLM */

/**
 * Converts B1950 RA, Dec to J2000 
 * Using method on page B42 of The Astronomical Almanac (1990 ed.)
 * Revised 90/05/07 J. J. Condon 
 * Taken from FITSview family
 * \param  ra    in/out Right Ascension in degrees
 * \param  dec   in/out Declination in degrees
 */
void ObitSkyGeomBtoJ (odouble *ra, odouble *dec)
{
  olong    i;
  odouble ra0, dec0, r0ta, rat=0.0, dect;
  odouble sina, cosa, sindd, cosdd, rnorm;
  odouble pi = 3.1415926536;
  odouble a[3] = {-1.62557e-6,-0.31919e-6,-0.13843e-6};
  odouble m1[3] = {0.9999256782, 0.0111820610, 0.0048579479};
  odouble m2[3] = {-0.0111820611, 0.9999374784, -0.0000271474};
  odouble m3[3] = {-0.0048579477, -0.0000271765, +0.9999881997};
  odouble r0[3], r1[3], r[3];
  
  /*    First convert input B1950 Ra, Dec to radians*/
  ra0 = *ra * pi / 180.0;
  dec0 = *dec * pi / 180.0;

  /*    Then convert B1950 RA, Dec to cartesian coordinates*/
  r0[0] = cos (ra0) * cos (dec0);
  r0[1] = sin (ra0) * cos (dec0);
  r0[2] = sin (dec0);

  /*    Remove aberration E-terms      */
  /*    (r0ta = scalar product of r0 and a)*/
  r0ta = r0[0] * a[0] + r0[1] * a[1] + r0[2] * a[2];
  for (i=0; i<3; i++) r1[i] = r0[i] - a[i] + r0ta * r0[i];

  /*    Precess from B1950 to J2000*/
  for (i=0; i<3; i++) r[i] = m1[i] * r1[0] + m2[i] * r1[1] + m3[i] * r1[2];
  
  /*    Convert J2000 Cartesian coordinates to J2000 RA, Dec (radians)*/
  rnorm = sqrt (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
  sindd = r[2] / rnorm;
  dect = asin (sindd);
  cosdd = sqrt (1.0 - sindd * sindd);
  cosa = r[0] / (rnorm * cosdd);
  sina = r[1] / (rnorm * cosdd);
  if (cosa != 0.0) rat = atan (sina / cosa);
  if (cosa == 0.0)
    {if (sina > 0.0) rat = pi / 2.0;
    if (sina < 0.0) rat = 1.50 * pi;}
  /*     Resolve 12h ambiguity of RA*/
  rat = rat * 12.0 / pi;
  if (cosa < 0.0) rat = rat + 12.0;
  if (rat < 0.0) rat = rat + 24.0;
  /*    Convert to degrees */
  *ra = rat * 15.0;
  *dec = dect * 180.0 / pi;
} /* end ObitSkyGeomBtoJ */

/**
 * Converts J2000 RA, Dec to B1950  
 * Using method on page B42 of The Astronomical Almanac (1990 ed.)
 * Revised 90/10/15 J. J. Condon 
 * Taken from FITSview family
 * \param  ra    in/out Right Ascension in degrees
 * \param  dec   in/out Declination in degrees
 */
void ObitSkyGeomJtoB (double *ra, double *dec)
{
  olong    i, iter;
  odouble ra0, dec0, sta, rat=0.0, dect;
  odouble sina, cosa, sindd, cosdd, rnorm, r1norm;
  odouble pi = 3.1415926536;
  odouble a[3] = {-1.62557e-6, -0.31919e-6, -0.13843e-6};
  odouble minv1[3] = {0.9999256795, -0.0111814828, -0.0048590040};
  odouble minv2[3] = {0.0111814828, 0.9999374849, -0.0000271557};
  odouble minv3[3] = {0.0048590039, -0.0000271771, +0.9999881946};
  odouble r0[3], r1[3], r[3], s[3], s1[3];

  /*    First convert input J2000 Ra, Dec to radians*/
  ra0 = *ra * pi / 180.0;
  dec0 = *dec * pi / 180.0;
  /*    Then convert J2000 RA, Dec to cartesian coordinates*/
  r0[0] = cos (ra0) * cos (dec0);
  r0[1] = sin (ra0) * cos (dec0);
  r0[2] = sin (dec0);

  /*    Precess from J2000 to B1950 */
  for (i=0; i<3; i++) r1[i] = minv1[i]*r0[0] + minv2[i]*r0[1] + minv3[i]*r0[2];

  /*    include aberration E-terms      */
  r1norm = sqrt (r1[0] * r1[0] + r1[1] * r1[1] + r1[2] * r1[2]);
  for (i=0; i<3; i++) s1[i] = r1[i] / r1norm;
  for (i=0; i<3; i++) s[i] = s1[i];
  
  /*    Three-Step iteration for r*/
  for (iter=0; iter<3; iter++)
    {
      /*  (sta = scalar product of s and a)*/
      sta = s[0] * a[0] + s[1] * a[1] + s[2] * a[2];
      /*  calculate or recalculate r*/
      for (i=0; i<3; i++) r[i] = s1[i] + a[i]  - sta * s[i];
      rnorm = sqrt (r[0] * r[0] + r[1] * r[1] + r[2] * r[2]);
      /*   calculate or recalculate s*/
      for (i=0; i<3; i++) s[i] = r[i] / rnorm;
    } /* end iteration loop */

  /*       Convert B1950 Cartesian coordinates (r-transpose) */
  /*       to B1950 RA, Dec (radians)*/
  sindd = r[2] / rnorm;
  dect = asin (sindd);
  cosdd = sqrt (1.0 - sindd * sindd);
  cosa = r[0] / (rnorm * cosdd);
  sina = r[1] / (rnorm * cosdd);
  if (cosa != 0.0) rat = atan (sina / cosa);
  if (cosa == 0.0) {
    if (sina > 0.0) rat = pi / 2.0;
    if (sina < 0.0) rat = 1.50 * pi;
  }

  /*      Then convert to deg of dec, hr of ra, */
  /*      resolve 12h ambiguity of RA */
  *dec = dect * 180.0 / pi;
  rat = rat * 12.0 / pi;
  if (cosa < 0.0) rat = rat + 12.0;
  if (rat < 0.0) rat = rat + 24.0;

  /*      Finally convert B1950 RA to degrees*/
  *ra = rat * 15.0;
} /* end ObitSkyGeomJtoB */

/**
 * Converts Convert Equatorial (B1950)to Galactic coordinates
 * Lifted from the AIPSish COORDT
 * \param  RALong    in/out Right Ascension/longitude in degrees
 * \param  DecLat    in/out Declination.latitude in degrees
 */
void ObitSkyGeomEq2Gal (odouble *RALong, odouble *DecLat)
{
  odouble ra, dec, glat, glong;

  /*  Celestial to Galactic */
  dec = *DecLat * DG2RAD;
  ra  = *RALong * DG2RAD;
  glat  = asin (sin(dec)*sin(DECGP) + cos(dec)*cos(DECGP)* cos(ra-RAGP));
  glong = LONCP + atan2 (cos(dec)*sin(RAGP-ra), sin(dec)* cos(DECGP) - 
			cos(dec)*sin(DECGP)*cos(ra-RAGP));
  if (glong >= 2.0*G_PI) glong -= - 2.0*G_PI;
  if (glong < 0.0)       glong += + 2.0*G_PI;
  *RALong = glong * RAD2DG;
  *DecLat = glat * RAD2DG;

} /* end ObitSkyGeomEq2Gal */

/**
 * Converts Convert Galactic to Equatorial (B1950) coordinates
 * Lifted from the AIPSish COORDT
 * \param  RALong    in/out Right Ascension/longitude in degrees
 * \param  DecLat    in/out Declination.latitude in degrees
 */
void ObitSkyGeomGal2Eq (odouble *RALong, odouble *DecLat)
{
  odouble ra, dec, glat, glong;

  /*  Galactic to Celestial */
  glat  = *DecLat * DG2RAD;
  glong = *RALong * DG2RAD;
  dec = asin (sin(glat)*sin(DECGP) + cos(glat)*cos(DECGP)*cos(glong-LONCP));
  ra = RAGP + atan2 (cos(glat)*sin(LONCP-glong), sin(glat)*cos(DECGP) - 
		     cos(glat)*sin(DECGP)*cos(glong-LONCP));
  if (ra >= 2.0*G_PI) ra = ra - 2.0*G_PI;
  if (ra < 0.0)       ra = ra + 2.0*G_PI;

  *RALong = ra * RAD2DG;
  *DecLat = dec * RAD2DG;
} /* end  ObitSkyGeomGal2Eq */

/**
 * Converts Convert Equatorial to Ecliptic coordinates
 * Lifted from the AIPSish COORDT
 * \param  RALong    in/out Right Ascension/longitude in degrees
 * \param  DecLat    in/out Declination.latitude in degrees
 * \param  epoch     Epoch of the coordinates to transform
 */
void ObitSkyGeomEq2Ec (odouble *RALong, odouble *DecLat, ofloat epoch)
{
  odouble ra, dec, elat, elong, dt, eps;
  odouble deps0[] = {23.452294*DG2RAD, -0.0130125*DG2RAD, -1.64e-6*DG2RAD, 5.03e-7*DG2RAD};

  dec = *DecLat * DG2RAD;
  ra  = *RALong * DG2RAD;

  dt = (epoch - 1900.0) / 100.0;
  eps = deps0[0] + dt * (deps0[1] + dt*(deps0[2] + dt*deps0[3]));

  elat  = asin (sin(dec)*cos(eps) - cos(dec)*sin(eps)* sin(ra));
  elong = atan2 (sin(dec)*sin(eps) - cos(dec)*cos(eps)*sin(ra), cos(dec)*cos(ra));
  if (elong >= 2.0*G_PI) elong = elong - 2.0*G_PI;
  if (elong < 0.0)       elong = elong + 2.0*G_PI;
  *RALong = elong * RAD2DG;
  *DecLat = elat * RAD2DG;
} /* end ObitSkyGeomEq2GEc */

/**
 * Converts Convert Ecliptic to Equatorial coordinates
 * Lifted from the AIPSish COORDT
 * \param  RALong    in/out Right Ascension/longitude in degrees
 * \param  DecLat    in/out Declination.latitude in degrees
 * \param  epoch     Epoch of the coordinates to transform
 */
void ObitSkyGeomEc2Eq (odouble *RALong, odouble *DecLat, ofloat epoch)
{
  odouble ra, dec, elat, elong, dt, eps;
  odouble deps0[] = {23.452294*DG2RAD, -0.0130125*DG2RAD, -1.64e-6*DG2RAD, 5.03e-7*DG2RAD};

  elat  = *DecLat * DG2RAD;
  elong = *RALong * DG2RAD;

  dt = (epoch - 1900.0) / 100.0;
  eps = deps0[0] + dt * (deps0[1] + dt*(deps0[2] + dt*deps0[3]));

  dec = asin (sin(elat)*cos(eps) + cos(elat)*sin(eps)*sin(elong));
  ra  = atan2 (cos(elat)*cos(eps)*sin(elong) - sin(elat)*sin(eps), 
	       cos(elat)*cos(elong));

  if (ra >= 2.0*G_PI) ra = ra - 2.0*G_PI;
  if (ra < 0.0)       ra = ra + 2.0*G_PI;

  *RALong = ra * RAD2DG;
  *DecLat = dec * RAD2DG;
} /* end ObitSkyGeomEc2Eq */

/**
 * Converts from celestial coordinates, expressed in terms of  
 * a reference position and a shift from this position to coordinates  
 * in a plane for fitting an Ionospheric phase screen model  
 * consisting of Zernike polynomials.  The output coordinates are  
 * normalized to unity at a 10 deg radius from the reference position.  
 * The coordinates are projected onto a plane tangent to the sky at  
 * the reference position.  
 * Routine translated from the AIPSish ZERGEOM.FOR/RD2ZER 
 * \param ra      Right Ascention of reference position (deg) 
 * \param dec     Declination of reference position (deg) 
 * \param xshift  Shift in X (RA) to desired position (deg) 
 * \param yshift  Shift in Y (Dec) to desired position (deg) 
 * \param xzer    [out] x-coordinate on Zernike plane 
 * \param yzer    [out] y-coordinate on Zernike plane 
 * \param ierr    0 ok, 1 out of range 
 */
void ObitSkyGeomRADec2Zern (odouble ra, odouble dec, 
			    ofloat xshift, ofloat yshift, 
			    ofloat* xzer, ofloat* yzer, olong *ierr) 
{
  odouble dx, dy, a, b, ax, ay, xxinc, coss, sins, sint, dt;
  odouble pi, twopi, zernor;

  /* Constants */
  pi = 3.14159265358979323846e0;
  twopi = 2.0e0*pi;
  /* Zernike normalization to Unit circle (10 deg) */
  zernor = 1.0 / (10.0e0 * DG2RAD);

  /* initial values */
  *xzer = 0.0;
  *yzer = 0.0;
  *ierr = 0;
  
  /* Convert to radians */
  a = DG2RAD * ra;
  b = DG2RAD * dec;
  /* Convert shift to position */
  xxinc = cos (DG2RAD * dec);
  if (xxinc != 0) {
    ax = ra + xshift / xxinc;
  } else {
    ax = ra;
  } 
  ay = dec + yshift;
  ax = ax * DG2RAD;
  ay = ay * DG2RAD;
  
  /* Get projection cosines */
  *ierr = 1;
  if (fabs(ay) > pi/2.0e0) return;
  if (fabs(b) > pi/2.0e0) return;
  *ierr = 0;
  coss = cos (ay);
  sins = sin (ay);
  dt = ax - a;
  if (dt > pi) dt = dt - twopi;
  if (dt < -pi) dt = dt + twopi;
  dx = sin (dt) * coss;
  sint = sins * sin (b) + coss * cos (b) * cos (dt);
  
  /* SIN projection */
  /*      IF (SINT.LT.0.0D0) IERR = 2 */
  if (sint < 0.0e0) {
    *ierr= 1;
    return;
  } 
  dy = sins * cos (b) - coss * sin (b) * cos (dt);

  /* Normalize to Zernike unit sphere */
  *xzer = dx * zernor;
  *yzer = dy * zernor;
  return;
} /* end of routine ObitSkyGeomRADec2Zern */ 

/*----------------------Private functions---------------------------*/
