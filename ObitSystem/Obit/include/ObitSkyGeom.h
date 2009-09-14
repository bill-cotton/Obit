/* $Id$  */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2004-2009                                          */
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
/*;  Correspondence this software should be addressed as follows:     */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
#ifndef OBITSKYGEOM_H
#define OBITSKYGEOM_H
#include "Obit.h"
#include "ObitErr.h"

/**
 * \file ObitSkyGeom.h
 * Celestial coordinates utility module
 * 
 * This file contains utility functions for celestial coordinates
 */

/*-------------- enumerations -------------------------------------*/
/**
 * \enum obitSkyGeomProj
 * enum for SkyGeom projection types.
 */
enum obitSkyGeomProj {
  /** -SIN  Sin projection */
  OBIT_SkyGeom_SIN = 0, 
  /** -TAN  Tan projection */
  OBIT_SkyGeom_TAN,
  /** -ARC Arc projection */
  OBIT_SkyGeom_ARC, 
  /** -NCP NCP (WSRT) projection */
  OBIT_SkyGeom_NCP, 
  /** -GLS Global sinusoid projection */
  OBIT_SkyGeom_GLS, 
  /** -MER Mercator projection */
  OBIT_SkyGeom_MER, 
  /** -AIT  Aitoff projection */
  OBIT_SkyGeom_AIT, 
  /** -STG Stereographic projection */
  OBIT_SkyGeom_STG 
}; /* end enum obitIOType */

/** typedef for enum for SkyGeom projection types. */
typedef enum obitSkyGeomProj ObitSkyGeomProj;

/*---------------Public functions---------------------------*/
/** Public: Determine shift between two positions */
void ObitSkyGeomShiftXY (odouble ra, odouble dec, ofloat rotate,
			odouble shiftRA, odouble shiftDec,
			ofloat *xShift, ofloat *yShift);

/** Public: Determine result of a shift to a position */
void ObitSkyGeomXYShift (odouble ra, odouble dec, 
			ofloat xShift, ofloat yShift, ofloat rotate,
			odouble *shiftRA, odouble *shiftDec);

/** Public: Get shift of coordinate reference pixel. */
void  
ObitSkyGeomShiftCRP (gchar *type, odouble ra, odouble dec, ofloat rotate,
		     odouble xra, double xdec, 
		     ofloat *xshift, ofloat *yshift);

/** Public: Get shift parameters for -SIN projection. */
void  
ObitSkyGeomShiftSIN (odouble ra, odouble dec, ofloat rotate,
		    odouble xra, double xdec, ofloat dxyzc[3]);

/** Public: Get shift parameters for -NCP projection. */
void  
ObitSkyGeomShiftNCP (odouble ra, odouble dec, ofloat rotate,
		    odouble xra, double xdec, ofloat dxyzc[3]);

/** Public: Returns astronomical coordinates given direction cosines, projection */
void 
ObitSkyGeomNewPos (ObitSkyGeomProj Proj, odouble ra0, odouble dec0, odouble l, odouble m, 
		   odouble *raout, odouble *decout, olong *ierr);

/** Public: accurate position for pixel coordinates */
olong 
ObitSkyGeomWorldPos(ofloat xpix, ofloat ypix, odouble xref, odouble yref, 
		    ofloat xrefpix, ofloat yrefpix, ofloat xinc, ofloat yinc, 
		    ofloat rot, gchar *type, odouble *xpos, odouble *ypos);

/** Public: Position for pixel coordinates from IRAF style CD matrix */
olong 
ObitSkyGeomCDpos(ofloat xpix, ofloat ypix, odouble xref, odouble yref,
		 ofloat xrefpix, ofloat yrefpix, ofloat xinc, ofloat yinc, ofloat rot,
		 ofloat cd1[2], ofloat cd2[2], gchar *type, odouble *xpos, odouble *ypos);

/** Public: Pixel coordinates for an RA and Dec*/
olong 
ObitSkyGeomXYpix(odouble xpos, odouble ypos, odouble xref, odouble yref, 
		 ofloat xrefpix, ofloat yrefpix, ofloat xinc, ofloat yinc, 
		 ofloat rot, gchar *type, ofloat *xpix, ofloat *ypix);

/** Public:pixel coordinates for an RA and Dec from IRAF  style CD matrix. */
olong 
ObitSkyGeomCDpix(odouble xpos, odouble ypos, odouble xref, odouble yref, 
		 ofloat xrefpix, ofloat yrefpix, ofloat xinc, ofloat yinc, ofloat rot,
		 ofloat icd1[2], ofloat icd2[2], gchar *type, 
		 ofloat *xpix, ofloat *ypix);

/** Public: Position for pixel coordinates from  offsets from the reference position.*/
olong 
ObitSkyGeomWorldPosLM(odouble dx, odouble dy, odouble xref, odouble yref, 
		      ofloat xinc, ofloat yinc, ofloat rot, gchar *type, 
		      odouble *xpos, odouble *ypos);

/** Public: Coordinate offsets for an RA and Dec   */
olong 
ObitSkyGeomXYPixLM(odouble xpos, odouble ypos, odouble xref, odouble yref, 
		   ofloat xinc, ofloat yinc, ofloat rot, gchar *type, 
		   odouble *dx, odouble *dy);

/** Public: Precess B1950 to J2000 coordinates  */
void 
ObitSkyGeomBtoJ (odouble *ra, odouble *dec);

/** Public: Precess J2000 to B1950 coordinates */
void 
ObitSkyGeomJtoB (odouble *ra, odouble *dec);

/** Public: Convert Equatorial (B1950) to Galactic coordinates  */
void ObitSkyGeomEq2Gal (odouble *RALong, odouble *DecLat);

/** Public: Convert Galactic to Equatorial (B1950)  */
void ObitSkyGeomGal2Eq (odouble *RALong, odouble *DecLat);

/** Public: Convert Equatorial to Ecliptic coordinates */
void ObitSkyGeomEq2Ec (odouble *RALong, odouble *DecLat, ofloat epoch);

/** Public: Convert Ecliptic to Equatorial */
void ObitSkyGeomEc2Eq (odouble *RALong, odouble *DecLat, ofloat epoch);

/** Public: Projection to Zernike plane */
void ObitSkyGeomRADec2Zern (odouble ra, odouble dec, ofloat xshift, ofloat yshift, 
			    ofloat* xzer, ofloat* yzer, olong *ierr);
#endif  /* OBITSKYGEOM_H */
