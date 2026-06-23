/* $Id: $   */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2026                                               */
/*;  Associated Universities, Inc. Washington DC, USA.                */
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
/*;  Correspondence this software should be addressed as follows:     */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
/*  Define the CUDASkyGeom (GPU version of ObitImageDesc+ObitSkyGeom) */
/* Have to put code in header to get it included with cuda kernal     */
/**
 * \file CUDASkyGeom.h
 * GPU version of ObitImageDesc+ObitSkyGeom utilities
 */
#ifndef CUDASKYGEOM_H 
#define CUDASKYGEOM_H

/*--------------Class definitions-------------------------------------*/
/** CUDASkyGeom Class structures: CUDAImageDesc, CUDAUVDesc */
/* Fooey - need local definition -
   I REALLY DON'T UNDERSTAND why #ifndef CUDAIMAGEDESCDEF_H doesn't work  */
#if HAVE_GPU==1  /* Have a GPU? */
#ifndef CUDAIMAGEDESCDEF_H // Prevent multiple definitions
#define CUDAIMAGEDESCDEF_H
#define IM_MAXDIM 7       /* maximum array dimension */
#define IMLEN_VALUE 41    /* Maximum length of descriptor string value */
#define IMLEN_KEYWORD 21  /* Maximum length of descriptor keyword  */
typedef struct {
#include "CUDAImageDescDef.h"
} CUDAImageDesc;
#else //DEBUG
//CUDAImageDesc *shouldnothappen;
#endif /* CUDAIMAGEDESCDEF_H */
#endif /* HAVE_GPU */

#include "ObitGPUSkyGeomDef.h"  // Ugly

#if IS_CUDA==1  /* CUDA code */
// now in ObitGPUSkyGeomDef.h #include "CUDASkyGeomDef.h"   /* this class definition */

/** Public: Create/initialize CUDAImageDesc structures */
extern "C"
CUDAImageDesc* CUDAImageDescCreate ();

/** Public: ImageDesc Destructor */
extern "C"
void  CUDAImageDescZap (CUDAImageDesc *in);

/** Public: Create/initialize CUDAUVDesc structures */
 extern "C"
CUDAUVDesc*  CUDAUVDescCreate ();

/** Public: UVDesc Destructor */
extern "C"
void CUDAUVDescZap (CUDAUVDesc *in);

/** Public: Coordinate offsets from Reference Position 
    from ObitPosition.c */
/** Public: Coordinate offsets from Reference Position 
    from ObitPosition.c: ObitPositionXYpixLM 
 * Routine to determine accurate coordinate offsets for an RA and Dec
 * offsets from the reference position.
 * Does: -SIN, -TAN, -ARC, -NCP, -GLS, -MER, -AIT projections,
 * anything else is linear
 *
 * \param coord   Coordinate value on each axis [deg]
 * \param desc    CUDA Image Descriptor
 * \param doff   [out] offset [deg] each axis
 * \return Pointer 0 if OK, else 1=angle too large for projection, 
 *   2 = bad values 
 */
extern "C" __device__ 
long CUDAPositionXYpixLM(double coord[2], CUDAImageDesc *desc,
		  	 double doff[2])
{
  double r, ra0, dec0, ra, dec, coss, sins, dt, da, dd, sint, ddx, ddy;
  double l, m=0.0, tmpm, geo1, geo2, geo3, sinr, cosr;
  double cond2r=1.745329252e-2, deps=1.0e-5, twopi=6.28318530717959;
  double xpos, ypos;
  long   i, itype;
  char ctypes[8][5] ={"-SIN","-TAN","-ARC","-NCP", "-GLS", "-MER", "-AIT", "-STG"};
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
  
  /*     check axis increments - bail out if either 0 */
  if ((desc->cdelt[0]==0.0) || (desc->cdelt[1]==0.0)) {doff[0]=0.0; doff[1]=0.0; return 2;}
  
  /*  find type  */
  itype = 0;  /* default type is linear */
  /* Can't use strncmp in device */
  for (i=0;i<8;i++)
    //if (!strncmp(&desc->ctype[0][4], ctypes[i], 4)) itype = i+1;
    if ((desc->ctype[0][4]==ctypes[i][0]) && (desc->ctype[0][5]==ctypes[i][1]) &&
       (desc->ctype[0][6]==ctypes[i][2]) && (desc->ctype[0][7]==ctypes[i][3])) itype = i+1;
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
    tmpm = m;
    if (m>1.0) tmpm = 1.0;
    m = acos (tmpm);
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
    /* trouble if you get here - can't call for help in the device */
    //fprintf(stderr,"Help\n");
    l = 0.0; m = 0.0;
  }  /* end of itype switch */
  
  /*   back to degrees  */
  doff[0] = l / cond2r;
  doff[1] = m / cond2r;
  return 0;
} /* end CUDAPositionXYpixLM */

/** Public: Pixels for  RA, Dec offsets
    from ObitPosition.c:ObitPositionXYpix
 * Routine to determine accurate pixel coordinates for an RA and Dec
 * offsets from the reference position.
 * Does: -SIN, -TAN, -ARC, -NCP, -GLS, -MER, -AIT projections,
 * anything else is linear
 *
 * \param coord   Coordinate value on each axis [deg]
 * \param desc    CUDA Image Descriptor
 * \param rotar   Rotation angle (rad)
 * \param pixel  [out] Pixel number (1-rel) on each axis
 * \return Pointer 0 if OK, else 1=angle too large for projection, 
 *   2 = bad values
 */
extern "C" __device__ 
long CUDAPositionXYpix(double coord[2], CUDAImageDesc *desc, float rotar, 
		       float pixel[2])
{
  double doff[2], dz;
  float sinr, cosr;
  long    iret;
   
   /* get coordinate offset */
   iret = CUDAPositionXYpixLM (coord, desc, doff);
   
   /*  Correct for rotation */
    __sincosf(rotar, &sinr, &cosr);
   dz      = doff[0]*cosr + doff[1]*sinr;
   doff[1] = doff[1]*cosr - doff[0]*sinr;
   doff[0] = dz;
   
   /*     convert to pixels  */
   pixel[0] = doff[0] / desc->cdelt[0] + desc->crpix[0];
   pixel[1] = doff[1] / desc->cdelt[1] + desc->crpix[1];
   return iret;
 } /* end GPUCUDAPositionXYpix */

#else  /* Not CUDA */
/** Public: Create/initialize CUDAImageDesc structures */
CUDAImageDesc* CUDAImageDescCreate ();

/** Public: ImageDesc Destructor */
void  CUDAImageDescZap (CUDAImageDesc *in);

/** Public: Create/initialize CUDAUVDesc structures */
CUDAUVDesc*  CUDAUVDescCreate ();

/** Public: UVDesc Destructor */
void CUDAUVDescZap (CUDAUVDesc *in);

/** Public: Coordinate offsets from Reference Position 
    from ObitPosition.c */
long CUDAPositionXYpixLM(double coord[2], void *desc,
		  	    double doff[2]);

/** Public: Pixels for  RA, Dec offsets
    from ObitPosition.c */
long CUDAPositionXYpix(double coord[2], void *desc,
		          float pixel[2]);
#endif /* IS_CUDA */
#endif  /* CUDASKYGEOM_H */
