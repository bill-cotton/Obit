/* $Id: ObitPennArrayUtil.c,v 1.1.1.1 2004/07/19 17:04:45 bcotton Exp $                            */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003                                               */
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

#include <math.h>
#include "ObitPennArrayUtil.h"
#include "ObitOTF.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitPennArrayUtil.c
 * Utility routine definitions.
 */


/*----------------------Public functions---------------------------*/
/**
 * Constructor for default GBT Penn Array ObitOTF.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitOTF* ObitPennArrayCreate (gchar* name)
{
  ObitOTF* out=NULL;
  ObitOTFArrayGeom *geom=NULL;
  olong numberDetect, ncol, nside, i, j, k, ncopy;
  ofloat x, y;
  gchar *Date = "20000919";

  /* Create basic structure */
  out = newObitOTF(name);

   /* create Array geometry with 64 elements */
  numberDetect = 64;
  geom = ObitOTFArrayGeomCreate (numberDetect);
  geom->azOffset = g_realloc(geom->azOffset, numberDetect*sizeof(ofloat));
  geom->elOffset = g_realloc(geom->elOffset, numberDetect*sizeof(ofloat));
 
  /* Fill with square array with 4" separations */
  nside = sqrt(numberDetect);
  k = 0;
  for (j=0; j<nside; j++) {
    y = (j-nside/2) * 4.0 / 3600.0; /* Elevation offset in deg */
    for (i=0; i<nside; i++) {
      x = (i-nside/2) * 4.0 / 3600.0; /* Azimuth offset in deg */
      geom->azOffset[k] = x;
      geom->elOffset[k] = y;
      k++;
    }
  }

  /* Other information - time info for 19 Sept. 2000 (random date) */
  ncopy = strlen (Date);
  for (i=0; i<ncopy; i++) geom->RefDate[i] = Date[i]; geom->RefDate[i]=0;
  geom->TimeSys[0] = 'U'; geom->TimeSys[1] = 'T';geom->TimeSys[2] = 'C';geom->TimeSys[3] = 0;
  geom->TeleX   =  882879.8949; /* 140 ft */
  geom->TeleY   = -4924482.3088;
  geom->TeleZ   =  3944130.6875;
  geom->DegDay  = 3.6098564497330e+02;
  geom->GSTiat0 = 3.5820740467876e+02;
  geom->PolarX  = 1.2269999831915e-02;
  geom->PolarY  = 2.4084000289440e-01;
  geom->ut1Utc  = 1.8335899710655e-01;
  geom->dataUtc = 0.0;
  geom->iatUtc  = 0.0;

  /* Compute some useful terms */
  /* telescope latitude in radians */
  if (fabs(geom->TeleX)<1.0) geom->TeleX = 1.0;
  geom->lat = asin (geom->TeleZ / 
		    sqrt(geom->TeleX*geom->TeleX + 
			 geom->TeleY*geom->TeleY + 
			 geom->TeleZ*geom->TeleZ));
  /* telescope longitude in radians */
  geom->lon = atan2 (geom->TeleY,  geom->TeleX);
  /* LST at iat0 in radians */
  geom->LSTiat0 = geom->GSTiat0*1.74533e-2 + geom->lon;
  /* Earth rotation rate in rad/day */
  geom->RadDay = geom->DegDay*1.74533e-2;
  /* Data - IAT in days */
  geom->dataIat = (geom->dataUtc - geom->iatUtc) / 86400.0;

  /* Attach Array geometry to OTF */
  out->geom = ObitOTFArrayGeomRef(geom);
  geom = ObitOTFArrayGeomUnref(geom);

 /* Initialize default GBT Penn Array Descriptor */
  strncpy (out->myDesc->object, "Simulation", OTFLEN_VALUE);
  strncpy (out->myDesc->teles,  "GBT       ", OTFLEN_VALUE);
  strncpy (out->myDesc->origin, "Obit Simulator", OTFLEN_VALUE);
  out->myDesc->isort[0] = 'T';  /* Time ordered */
  ncol = 0;

  out->myDesc->JDObs = 0.0;
  out->myDesc->epoch = 2000.0;
  out->myDesc->equinox = 2000.0;
  strncpy (out->myDesc->bunit,  "COUNTS  ", OTFLEN_VALUE);
  strncpy (out->myDesc->obsdat, "2000-09-19", OTFLEN_VALUE);
 /* Time */
  strncpy (out->myDesc->colType[ncol], "TIME    ", OTFLEN_KEYWORD);
  strncpy (out->myDesc->colUnit[ncol], "DAYS    ", OTFLEN_VALUE);
  out->myDesc->colRepeat[ncol] = 1;
  ncol++;

  /* Integration time */
  strncpy (out->myDesc->colType[ncol], "TIME_INT", OTFLEN_KEYWORD);
  strncpy (out->myDesc->colUnit[ncol], "DAYS    ", OTFLEN_VALUE);
  out->myDesc->colRepeat[ncol] = 1;
  ncol++;

  /* Target index */
  strncpy (out->myDesc->colType[ncol], "TARGET  ", OTFLEN_KEYWORD);
  strncpy (out->myDesc->colUnit[ncol], "        ", OTFLEN_VALUE);
  out->myDesc->colRepeat[ncol] = 1;
  ncol++;

  /* Scan index */
  strncpy (out->myDesc->colType[ncol], "SCAN    ", OTFLEN_KEYWORD);
  strncpy (out->myDesc->colUnit[ncol], "        ", OTFLEN_VALUE);
  out->myDesc->colRepeat[ncol] = 1;
  ncol++;

  /* Pointing RA */
  strncpy (out->myDesc->colType[ncol], "RA      ", OTFLEN_KEYWORD);
  strncpy (out->myDesc->colUnit[ncol], "DEGREE  ", OTFLEN_VALUE);
  out->myDesc->colRepeat[ncol] = 1;
  ncol++;

  /* Pointing Dec */
  strncpy (out->myDesc->colType[ncol], "DEC     ", OTFLEN_KEYWORD);
  strncpy (out->myDesc->colUnit[ncol], "DEGREE  ", OTFLEN_VALUE);
  out->myDesc->colRepeat[ncol] = 1;
  ncol++;

  /* Rotation of array on sky */
  strncpy (out->myDesc->colType[ncol], "ROTATE  ", OTFLEN_KEYWORD);
  strncpy (out->myDesc->colUnit[ncol], "DEGREE  ", OTFLEN_VALUE);
  out->myDesc->colRepeat[ncol] = 1;
  ncol++;

  /* Cal on? */
  strncpy (out->myDesc->colType[ncol], "CAL     ", OTFLEN_KEYWORD);
  strncpy (out->myDesc->colUnit[ncol], "        ", OTFLEN_VALUE);
  out->myDesc->colRepeat[ncol] = 1;
  ncol++;

  /* Data - MUST be last column */
  out->myDesc->numDesc = ncol;
  strncpy (out->myDesc->colType[ncol], "DATA    ", OTFLEN_KEYWORD);
  strncpy (out->myDesc->colUnit[ncol], "COUNTS  ", OTFLEN_VALUE);
  out->myDesc->colRepeat[ncol] = numberDetect;
  ncol++;

  out->myDesc->ncol = ncol;

  /* Data array descriptors */
  out->myDesc->naxis = 0;
  
  /* Detector axis */
  out->myDesc->inaxes[out->myDesc->naxis] = numberDetect;
  strncpy (out->myDesc->ctype[out->myDesc->naxis], "FEED", OTFLEN_KEYWORD);
  out->myDesc->cdelt[out->myDesc->naxis] = 1.0;
  out->myDesc->crpix[out->myDesc->naxis] = 1.0;
  out->myDesc->crota[out->myDesc->naxis] = 0.0;
  out->myDesc->crval[out->myDesc->naxis] = 1.0;
  out->myDesc->naxis++;
  
  /* Stokes axis */
  out->myDesc->inaxes[out->myDesc->naxis] = 1;
  strncpy (out->myDesc->ctype[out->myDesc->naxis], "STOKES  ", OTFLEN_KEYWORD);
  out->myDesc->cdelt[out->myDesc->naxis] = 1.0;
  out->myDesc->crpix[out->myDesc->naxis] = 1.0;
  out->myDesc->crota[out->myDesc->naxis] = 0.0;
  out->myDesc->crval[out->myDesc->naxis] = 1.0;
  out->myDesc->naxis++;
  
  /* Frequency axis */
  out->myDesc->inaxes[out->myDesc->naxis] = 1;
  strncpy (out->myDesc->ctype[out->myDesc->naxis], "FREQ    ", OTFLEN_KEYWORD);
  out->myDesc->cdelt[out->myDesc->naxis] = 1.0;
  out->myDesc->crpix[out->myDesc->naxis] = 1.0;
  out->myDesc->crota[out->myDesc->naxis] = 0.0;
  out->myDesc->crval[out->myDesc->naxis] = 100.0e9;
  out->myDesc->naxis++;
  
  /* Index the descriptor */
  ObitOTFDescIndex (out->myDesc);

 return out;
} /* end ObitPennArrayCreate */




