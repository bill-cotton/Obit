/* $Id$      */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 1996,1997-2008-2008                                */
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

#include "ObitPosition.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitPosition.c
 * ObitPosition module function definitions.
 * Utility functions for coordinates
 * Adopted from XFITSview/FITSview which adapted from AIPS
 */

/*---------------Private function prototypes----------------*/

/*----------------------Public functions---------------------------*/
/**
 * Routine to determine accurate position for pixel coordinates.
 * Does: -SIN, -TAN, -ARC, -NCP, -GLS, -MER, -AIT projections,
 * anything else is linear
 *
 * \param pixel   Pixel number (1-rel) on each axis
 * \param desc    Obit Image Descriptor
 * \param coord   [out] Coordinate value on each axis
 * \return Pointer 0 if OK, else angle too large for projection
 */
olong ObitPositionWorldPos(ofloat pixel[2], ObitImageDesc *desc,
      odouble coord[2])
 {
   odouble cosr, sinr, doff[2], temp;
   odouble cond2r=1.745329252e-2;
   
   /*   Offset from ref pixel  */
   doff[0] = (pixel[0]-desc->crpix[0]) * desc->cdelt[0];
   doff[1] = (pixel[1]-desc->crpix[1]) * desc->cdelt[1];

   /*   Take out rotation  */
   cosr = cos(desc->crota[1]*cond2r);
   sinr = sin(desc->crota[1]*cond2r);
   if (desc->crota[1]!=0.0) {
     temp   = doff[0] * cosr - doff[1] * sinr;
     doff[1] = doff[1] * cosr + doff[0] * sinr;
     doff[0] = temp;
   }
   /* determine position */
   return ObitPositionWorldPosLM (doff, desc, coord);
 }  /* End of ObitPositionWorldPos */

/**
 * Routine to determine accurate position for pixel coordinates from 
 * offsets from the reference position.
 * Does: -SIN, -TAN, -ARC, -NCP, -GLS, -MER, -AIT projections,
 * anything else is linear
 *
 * \param doff    Coordinate offsets on each axis (deg)
 * \param desc    Obit Image Descriptor
 * \param coord   [out] Coordinate on each axis (deg) 
 * \return Pointer 0 if OK, else angle too large for projection
 */
olong ObitPositionWorldPosLM(odouble doff[2], ObitImageDesc *desc,
			    odouble coord[2])
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
   cosr = cos(desc->crota[1]*cond2r);
   sinr = sin(desc->crota[1]*cond2r);
   
   /*  find type  */
   itype = 0;  /* default type is linear */
   for (i=0;i<8;i++) if (!strncmp(&desc->ctype[0][4], ctypes[i], 4)) itype = i+1;

   /* default, linear result for error return  */
   coord[0] = desc->crval[0] + doff[0];
   coord[1] = desc->crval[1] + doff[1];

   /* convert to radians  */
   ra0  = desc->crval[0] * cond2r;
   dec0 = desc->crval[1] * cond2r;
   l = doff[0] * cond2r;
   m = doff[1] * cond2r;
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
     dect = acos (MIN (dect, 1.000));
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
     dt = desc->cdelt[1] * cosr + desc->cdelt[0] * sinr;
     if (dt==0.0) dt = 1.0;
     doff[1] = (desc->crval[1]/2.0 + 45.0) * cond2r;
     doff[0] = doff[1] + dt / 2.0 * cond2r;
     doff[1] = log (tan (doff[1]));
     doff[0] = log (tan (doff[0]));
     geo2 = dt * cond2r / (doff[0] - doff[1]);
     geo3 = geo2 * doff[1];
     geo1 = cos (desc->crval[1]*cond2r);
     if (geo1<=0.0) geo1 = 1.0;
     rat = l / geo1 + ra0;
     if (fabs(rat - ra0) > twopi) return 1; /* added 10/13/94 DCW/EWG */
     dt = 0.0;
     if (geo2!=0.0) dt = (m + geo3) / geo2;
     dt = exp (dt);
     dect = 2.0 * atan (dt) - twopi / 4.0;
     break;

   case 7:   /* -AIT Aitoff*/
     dt = desc->cdelt[1]*cosr + desc->cdelt[0]*sinr;
     if (dt==0.0) dt = 1.0;
     dt = dt * cond2r;
     doff[1] = desc->crval[1] * cond2r;
     doff[0] = sin(doff[1]+dt)/sqrt((1.0+cos(doff[1]+dt))/2.0) -
       sin(doff[1])/sqrt((1.0+cos(doff[1]))/2.0);
     if (doff[0]==0.0) doff[0] = 1.0;
     geo2 = dt / doff[0];
     dt = desc->cdelt[0]*cosr - desc->cdelt[1]* sinr;
     if (dt==0.0) dt = 1.0;
     dt = dt * cond2r;
     doff[0] = 2.0 * cos(doff[1]) * sin(dt/2.0);
     if (doff[0]==0.0) doff[0] = 1.0;
     geo1 = dt * sqrt((1.0+cos(doff[1])*cos(dt/2.0))/2.0) / doff[0];
     geo3 = geo2 * sin(doff[1]) / sqrt((1.0+cos(doff[1]))/2.0);
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

   default:
     fprintf(stderr,"Help\n");
     /* really bad if you get here */
   } /* end of switch by type */

   /*  return ra in range  */
   raout = rat;
   decout = dect;
   if (raout-ra0>twopi/2.0) raout = raout - twopi;
   if (raout-ra0<-twopi/2.0) raout = raout + twopi;
   if (raout < 0.0) raout += twopi; /* added by DCW 10/12/94 */
   
   /*  correct units back to degrees  */
   coord[0]  = raout  / cond2r;
   coord[1]  = decout  / cond2r;
   return 0;
 }  /* End of ObitPositionWorldPosLM */

/**
 * Routine to determine accurate pixel coordinates for an RA and Dec
 * offsets from the reference position.
 * Does: -SIN, -TAN, -ARC, -NCP, -GLS, -MER, -AIT projections,
 * anything else is linear
 *
 * \param coord   Coordinate value on each axis
 * \param desc    Obit Image Descriptor
 * \param pixel   [out] Pixel number (1-rel) on each axis
 * \return Pointer 0 if OK, else 1=angle too large for projection, 
 *   2 = bad values 
 */
olong ObitPositionXYpix(odouble coord[2], ObitImageDesc *desc,
      ofloat pixel[2])
 {
   odouble doff[2], dz, sinr, cosr;
   olong    iret;
   odouble cond2r=1.745329252e-2;
   
   /* get coordinate offset */
   iret = ObitPositionXYpixLM (coord, desc, doff);
   
   /*  Correct for rotation */
   cosr = cos(desc->crota[1]*cond2r);
   sinr = sin(desc->crota[1]*cond2r);
   dz      = doff[0]*cosr + doff[1]*sinr;
   doff[1] = doff[1]*cosr - doff[0]*sinr;
   doff[0] = dz;
   
   /*     convert to pixels  */
   pixel[0] = doff[0] / desc->cdelt[0] + desc->crpix[0];
   pixel[1] = doff[1] / desc->cdelt[1] + desc->crpix[1];
   return iret;
 }  /* end ObitPositionXYpix */

/**
 * Routine to determine accurate coordinate offsets for an RA and Dec
 * offsets from the reference position.
 * Does: -SIN, -TAN, -ARC, -NCP, -GLS, -MER, -AIT projections,
 * anything else is linear
 *
 * \param coord   Coordinate value on each axis
 * \param desc    Obit Image Descriptor
 * \param doff   [out] Pixel number (1-rel) on each axis
 * \return Pointer 0 if OK, else 1=angle too large for projection, 
 *   2 = bad values 
 */
olong ObitPositionXYpixLM(odouble coord[2], ObitImageDesc *desc,
			 odouble doff[2])
{
  odouble dz, r, ra0, dec0, ra, dec, coss, sins, dt, da, dd, sint, ddx, ddy;
  odouble l, m=0.0, geo1, geo2, geo3, sinr, cosr;
  odouble cond2r=1.745329252e-2, deps=1.0e-5, twopi=6.28318530717959;
  odouble xpos, ypos;
  olong   i, itype;
  gchar ctypes[8][5] ={"-SIN","-TAN","-ARC","-NCP", "-GLS", "-MER", "-AIT",
		       "-STG"};
  /* local values */
  xpos = coord[0];
  ypos = coord[1];
  
  /* 0h wrap-around tests added by D.Wells 10/12/94: */
  dt = (xpos - desc->crval[0]);
  if (dt > +180) xpos -= 360;
  if (dt < -180) xpos += 360;
  
  /* default values - linear */
  doff[0] = xpos - desc->crval[0];
  doff[1] = ypos - desc->crval[1];
  dz = 0.0;
  
  /*     check axis increments - bail out if either 0 */
  if ((desc->cdelt[0]==0.0) || (desc->cdelt[1]==0.0)) {doff[0]=0.0; doff[1]=0.0; return 2;}
  
  /*  find type  */
  itype = 0;  /* default type is linear */
  for (i=0;i<8;i++) if (!strncmp(&desc->ctype[0][4], ctypes[i], 4)) itype = i+1;
  if (itype==0) return 0;  /* done if linear */
  
  /*  rotation */
  r = desc->crota[1] * cond2r;
  cosr = cos (r);
  sinr = sin (r);
  
  /* Non linear position */
  ra0  = desc->crval[0] * cond2r;
  dec0 = desc->crval[1] * cond2r;
  ra   = xpos * cond2r;
  dec  = ypos * cond2r;
  
  /* compute direction cosine */
  coss = cos (dec);
  sins = sin (dec);
  l    = sin(ra-ra0) * coss;
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
    m = acos (MIN (m, 1.000));
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
    dt = desc->cdelt[1] * cosr + desc->cdelt[0] * sinr;
    if (dt==0.0) dt = 1.0;
    ddy = (desc->crval[1]/2.0 + 45.0) * cond2r;
    ddx = ddy + dt / 2.0 * cond2r;
    ddy = log (tan (ddy));
    ddx = log (tan (ddx));
    geo2 = dt * cond2r / (ddx - ddy);
    geo3 = geo2 * ddy;
    geo1 = cos (desc->crval[1]*cond2r);
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
    dt = desc->cdelt[1]*cosr + desc->cdelt[0]*sinr;
    if (dt==0.0) dt = 1.0;
    dt = dt * cond2r;
    ddy = desc->crval[1] * cond2r;
    ddx = sin(ddy+dt)/sqrt((1.0+cos(ddy+dt))/2.0) -
      sin(ddy)/sqrt((1.0+cos(ddy))/2.0);
    if (ddx==0.0) ddx = 1.0;
    geo2 = dt / ddx;
    dt = desc->cdelt[0]*cosr - desc->cdelt[1]* sinr;
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
  default:
    /* trouble if you get here */
    fprintf(stderr,"Help\n");
  }  /* end of itype switch */
  
  /*   back to degrees  */
  doff[0] = l / cond2r;
  doff[1] = m / cond2r;
  return 0;
}  /* end ObitPositionXYpixLM */

