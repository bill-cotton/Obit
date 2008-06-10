/* $Id: SkyGeom.inc,v 1.3 2007/10/30 13:04:20 bcotton Exp $  */  
/*--------------------------------------------------------------------*/
/* Swig module description for Image utilities                        */
/*                                                                    */
/*;  Copyright (C) 2007                                               */
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

%{
#include "ObitSkyGeom.h"
#include "ObitPBUtil.h"
%}



%inline %{

/** Public: Determine shift between two positions */
void SkyGeomShiftXY (double ra, double dec, float rotate,
			double shiftRA, double shiftDec,
			float* outFlt1, float* outFlt2) {
  ObitSkyGeomShiftXY ((odouble)ra, (odouble)dec, (ofloat)rotate,
		     (odouble)shiftRA, (odouble)shiftDec,
		      outFlt1, outFlt2);
}

/** Public: Determine result of a shift to a position */
void SkyGeomXYShift (double ra, double dec, 
			float xShift, float yShift, float rotate,
			double* outDbl1, double* outDbl2) {
  ObitSkyGeomXYShift ((odouble)ra, (odouble)dec, 
		     (ofloat)xShift, (ofloat)yShift, (ofloat)rotate,
                      outDbl1, outDbl2);
}

/** Public: Returns astronomical coordinates given direction cosines, projection */
void 
SkyGeomNewPos (char *type, double ra0, double dec0, double l, double m, 
		   double* outDbl1, double* outDbl2, int* outInt1) {
  ObitSkyGeomProj Proj = OBIT_SkyGeom_SIN;
  /* Projection type */
  if (!strncmp (type, "-TAN", 4)) Proj = OBIT_SkyGeom_TAN;
  if (!strncmp (type, "-ARC", 4)) Proj = OBIT_SkyGeom_ARC;
  if (!strncmp (type, "-NCP", 4)) Proj = OBIT_SkyGeom_NCP;
  if (!strncmp (type, "-GLS", 4)) Proj = OBIT_SkyGeom_GLS;
  if (!strncmp (type, "-MER", 4)) Proj = OBIT_SkyGeom_MER;
  if (!strncmp (type, "-AIT", 4)) Proj = OBIT_SkyGeom_AIT;
  if (!strncmp (type, "-STG", 4)) Proj = OBIT_SkyGeom_STG;
  ObitSkyGeomNewPos (Proj, (odouble)ra0, (odouble)dec0, (odouble)l, (odouble)m, 
		     outDbl1, outDbl2, outInt1);
}

/** Public: accurate position for pixel coordinates */
int 
SkyGeomWorldPos(float xpix, float ypix, double xref, double yref, 
		    float xrefpix, float yrefpix, float xinc, float yinc, 
		    float rot, char *type, double *outDbl1, double *outDbl2) {
  return ObitSkyGeomWorldPos((ofloat)xpix, (ofloat)ypix, (odouble)xref, (odouble)yref, 
		    (ofloat)xrefpix, (ofloat)yrefpix, (ofloat)xinc, (ofloat)yinc, 
		    (ofloat)rot, (gchar*)type, outDbl1, outDbl2);
}

/** Public: Position for pixel coordinates from IRAF style CD matrix */
int 
SkyGeomCDpos(float xpix, float ypix, double xref, double yref,
		 float xrefpix, float yrefpix, float xinc, float yinc, float rot,
		 float cd1[2], float cd2[2], char *type, double *outDbl1, double *outDbl2) {
  ofloat lcd1[2], lcd2[2];
  lcd1[0] = cd1[0];
  lcd1[1] = cd1[1];
  lcd2[0] = cd2[0];
  lcd2[1] = cd2[1];
  return ObitSkyGeomCDpos((ofloat)xpix, (ofloat)ypix, (odouble)xref, (odouble)yref,
		 (ofloat)xrefpix, (ofloat)yrefpix, (ofloat)xinc, (ofloat)yinc, (ofloat)rot,
		 lcd1, lcd2, (gchar*)type, outDbl1, outDbl2);
}

/** Public: Pixel coordinates for an RA and Dec*/
int 
SkyGeomXYpix(double xpos, double ypos, double xref, double yref, 
		 float xrefpix, float yrefpix, float xinc, float yinc, 
		 float rot, char *type, float *outFlt1, float *outFlt2){
  return ObitSkyGeomXYpix((odouble)xpos, (odouble)ypos, (odouble)xref, (odouble)yref, 
		 (ofloat)xrefpix, (ofloat)yrefpix, (ofloat)xinc, (ofloat)yinc, 
		 (ofloat)rot, (gchar*)type, outFlt1, outFlt2);
}

/** Public:pixel coordinates for an RA and Dec from IRAF  style CD matrix. */
int 
SkyGeomCDpix(double xpos, double ypos, double xref, double yref, 
		 float xrefpix, float yrefpix, float xinc, float yinc, float rot,
		 float icd1[2], float icd2[2], char *type, 
	 	 float *outFlt1, float *outFlt2) {
  ofloat licd1[2], licd2[2];
  licd1[0] = icd1[0];
  licd1[1] = icd1[1];
  licd2[0] = icd2[0];
  licd2[1] = icd2[1];
  return ObitSkyGeomCDpix((odouble)xpos, (odouble)ypos, (odouble)xref, (odouble)yref, 
		 (ofloat)xrefpix, (ofloat)yrefpix, (ofloat)xinc, (ofloat)yinc, (ofloat)rot,
		 licd1, licd2, (gchar*)type, outFlt1, outFlt2);
}

/** Public: Position for pixel coordinates from  offsets from the reference position.*/
int 
SkyGeomWorldPosLM(double dx, double dy, double xref, double yref, 
		      float xinc, float yinc, float rot, char *type, 
		      double *outDbl1, double *outDbl2) {
  return ObitSkyGeomWorldPosLM((odouble)dx, (odouble)dy, (odouble)xref, (odouble)yref, 
		      (ofloat)xinc, (ofloat)yinc, (ofloat)rot, (gchar*)type, 
		      outDbl1, outDbl2);
}

/** Public: Coordinate offsets for an RA and Dec   */
int 
SkyGeomXYPixLM(double xpos, double ypos, double xref, double yref, 
		   float xinc, float yinc, float rot, char *type, 
		   double *outDbl1, double *outDbl2) {
  return ObitSkyGeomXYPixLM((odouble)xpos, (odouble)ypos, (odouble)xref, (odouble)yref, 
		   (ofloat)xinc, (ofloat)yinc, (ofloat)rot, (gchar*)type, 
		    outDbl1, outDbl2);
}

/** Public: Precess B1950 to J2000 coordinates  */
void 
SkyGeomBtoJ (double *outDbl1, double *outDbl2) {
  ObitSkyGeomBtoJ (outDbl1, outDbl2);
}

/** Public: Precess J2000 to B1950 coordinates */
void 
SkyGeomJtoB (double *outDbl1, double *outDbl2) {
  ObitSkyGeomJtoB (outDbl1, outDbl2);
}

/** Public: Convert Equatorial (B1950) to Galactic coordinates  */
void SkyGeomEq2Gal (double *outDbl1, double *outDbl2) {
  ObitSkyGeomEq2Gal (outDbl1, outDbl2);
}

/** Public: Convert Galactic to Equatorial (B1950)  */
void SkyGeomGal2Eq (double *outDbl1, double *outDbl2) {
  ObitSkyGeomGal2Eq (outDbl1, outDbl2);
}

/** Public: Convert Equatorial to Ecliptic coordinates */
void SkyGeomEq2Ec (double *outDbl1, double *outDbl2, float epoch){
  ObitSkyGeomEq2Ec (outDbl1, outDbl2, (ofloat)epoch);
}

/** Public: Convert Ecliptic to Equatorial */
void SkyGeomEc2Eq (double *outDbl1, double *outDbl2, float epoch) {
  ObitSkyGeomEc2Eq (outDbl1, outDbl2, (ofloat)epoch);
}

/** Public: Projection to Zernike plane */
void SkyGeomRADec2Zern (double ra, double dec, float xshift, float yshift, 
			    float* outFlt1, float* outFlt2, int *outInt1) {
  ObitSkyGeomRADec2Zern ((odouble)ra, (odouble)dec, (ofloat)xshift, (ofloat)yshift, 
			    outFlt1, outFlt2, outInt1);
}

/** Public: Compute VLA beam shape from a fitted polynomial */
float PBUtilPoly (double Angle, double Freq, float pbmin) {
  ofloat out;
  out = ObitPBUtilPoly ((odouble)Angle, (odouble)Freq, (ofloat)pbmin);
  return (float)out;
}

/** Public:  Compute Antenna beam assuming uniform illumination of an antenna */
float PBUtilJinc (double Angle, double Freq, float antSize, float pbmin) {
  ofloat out;
  out = ObitPBUtilJinc ((odouble)Angle, (odouble)Freq, (ofloat)antSize, (ofloat)pbmin);
  return (float)out;
}

/** Public:  Calculates the relative gain at a reference frequency */
float PBUtilRelPB (double Angle, int nfreq, double *Freq, float antSize, 
                   float pbmin, double refFreq) {
  ofloat out;
  out = ObitPBUtilRelPB ((odouble)Angle, (olong)nfreq, (odouble*)Freq, 
        (ofloat)antSize, (ofloat)pbmin, (odouble)refFreq);
  return (float)out;
}

%}
