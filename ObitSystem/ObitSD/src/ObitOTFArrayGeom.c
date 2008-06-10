/* $Id$*/
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2008                                          */
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

#include "ObitOTFArrayGeom.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitOTFArrayGeom.c
 * GBT/OTF array geometry class function definitions.
 * This class is derived from the Obit base class.
 */

/*--------------- File Global Variables  ----------------*/
/** name of the class defined in this file */
static gchar *myClassName = "ObitOTFArrayGeom";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo structure ObitOTFArrayGeomClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitOTFArrayGeomClassInfo myClassInfo = {FALSE};

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitOTFArrayGeomInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitOTFArrayGeomClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitOTFArrayGeomClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitOTFArrayGeom* newObitOTFArrayGeom (gchar* name)
{
  ObitOTFArrayGeom* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitOTFArrayGeomClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitOTFArrayGeom));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitOTFArrayGeomInit((gpointer)out);

 return out;
} /* end newObitOTFArrayGeom */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitOTFArrayGeomGetClass (void)
{
  return (gconstpointer)&myClassInfo;
} /* end ObitOTFArrayGeomGetClass */

/**
 * Make a deep copy of an ObitOTFArrayGeom.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitOTFArrayGeom* 
ObitOTFArrayGeomCopy  (ObitOTFArrayGeom *in, ObitOTFArrayGeom *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  gboolean oldExist;
  olong i;
  gchar *outName;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return out;
  g_assert (ObitIsA(in, &myClassInfo));
  if (out) g_assert (ObitIsA(out, &myClassInfo));

  /* Create if it doesn't exist */
  oldExist = out!=NULL;
  if (!oldExist) {
    /* derive object name */
    outName = g_strconcat ("Copy: ",in->name,NULL);
    out = newObitOTFArrayGeom(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  out->info = ObitInfoListUnref(out->info);
  out->info = ObitInfoListCopy(in->info);

  /* this class data */
  out->numberDetect = in->numberDetect;
  out->azOffset = g_realloc(out->azOffset, in->numberDetect*sizeof(ofloat));
  out->elOffset = g_realloc(out->elOffset, in->numberDetect*sizeof(ofloat));
  for (i=0; i<in->numberDetect; i++) out->azOffset[i] = in->azOffset[i];
  for (i=0; i<in->numberDetect; i++) out->elOffset[i] = in->elOffset[i];
  for (i=0; i<12; i++) out->RefDate[i] = in->RefDate[i];
  for (i=0; i<4; i++)  out->TimeSys[i] = in->TimeSys[i]; 
  out->TeleX   = in->TeleX;
  out->TeleY   = in->TeleY;
  out->TeleZ   = in->TeleZ;
  out->DegDay  = in->DegDay;
  out->GSTiat0 = in->GSTiat0;
  out->PolarX  = in->PolarX;
  out->PolarY  = in->PolarY;
  out->ut1Utc  = in->ut1Utc;
  out->dataUtc = in->dataUtc;
  out->iatUtc  = in->iatUtc;
 return out;
} /* end ObitOTFArrayGeomCopy */

/**
 * Make a deep copy of an ObitOTFArrayGeom averaging all frequency channels
 * \param in      The object to copy
 * \param inDesc  OTFdescriptor for in
 * \param out     An existing object pointer for output or NULL if none exists.
 * \param outDesc OTFdescriptor for out
 * \param err     Obit error stack object.
 * \return pointer to the new object.
 */
ObitOTFArrayGeom* 
ObitOTFArrayGeomAver  (ObitOTFArrayGeom *in, ObitOTFDesc *inDesc, 
		       ObitOTFArrayGeom *out, ObitOTFDesc *outDesc, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  gboolean oldExist;
  olong i, iRow, oRow;
  olong nstoke, nchan, nfeed;
  olong  istoke, ifeed, ichan, isoff, ifoff, osoff, ofoff;
  gchar *outName;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return out;
  g_assert (ObitOTFArrayGeomIsA(in));
  if (out) g_assert (ObitOTFArrayGeomIsA(out));

  /* Create if it doesn't exist */
  oldExist = out!=NULL;
  if (!oldExist) {
    /* derive object name */
    outName = g_strconcat ("Copy: ",in->name,NULL);
    out = newObitOTFArrayGeom(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = ((ObitClassInfo*)in->ClassInfo)->ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  out->info = ObitInfoListUnref(out->info);
  out->info = ObitInfoListCopy(in->info);

  /* this class data */
  for (i=0; i<12; i++) out->RefDate[i] = in->RefDate[i];
  for (i=0; i<4; i++)  out->TimeSys[i] = in->TimeSys[i]; 
  out->TeleX   = in->TeleX;
  out->TeleY   = in->TeleY;
  out->TeleZ   = in->TeleZ;
  out->DegDay  = in->DegDay;
  out->GSTiat0 = in->GSTiat0;
  out->PolarX  = in->PolarX;
  out->PolarY  = in->PolarY;
  out->ut1Utc  = in->ut1Utc;
  out->dataUtc = in->dataUtc;
  out->iatUtc  = in->iatUtc;
  
  /* Allocate output offset arrays */
  out->numberDetect = 
    outDesc->colRepeat[outDesc->ncol-1] / outDesc->incdatawt; /* number of detectors */
  out->azOffset = g_realloc(out->azOffset, out->numberDetect*sizeof(ofloat));
  out->elOffset = g_realloc(out->elOffset, out->numberDetect*sizeof(ofloat));

  nstoke = inDesc->inaxes[inDesc->jlocs];
  nfeed  = inDesc->inaxes[inDesc->jlocfeed];
  nchan  = outDesc->inaxes[outDesc->jlocf];

  /* Copy selected part of table, entries in same order as data */
  /* loop over data - Stokes outer loop */
  for (istoke=0; istoke<nstoke; istoke++) {
    isoff = istoke * inDesc->incs / inDesc->incdatawt; /* offset in data */
    osoff = istoke * outDesc->incs / outDesc->incdatawt;
    
    /* Loop over feeds */
    for (ifeed=0; ifeed<nfeed; ifeed++) {
      ifoff = isoff + ifeed * inDesc->incfeed / inDesc->incdatawt;
      ofoff = osoff + ifeed * outDesc->incfeed / outDesc->incdatawt;
      
      /* Over output Channel */
      for (ichan=0; ichan<nchan; ichan++) {
	iRow = ifoff + ichan * inDesc->incf / inDesc->incdatawt;  /* Input element number */
	oRow = ofoff + ichan * outDesc->incf / outDesc->incdatawt; /* Output element number */

	/* Copy entry */
	out->azOffset[oRow] = in->azOffset[iRow];
	out->elOffset[oRow] = in->elOffset[iRow];
      }

    } /* end feed loop */
  } /* end Stokes */

  return out;
} /* end ObitOTFArrayGeomAver */

/**
 * Creates an ObitOTFArrayGeom for a specified number of detectors.
 * \param ndetect  Number of detectors in array
 * \return the new object.
 */
ObitOTFArrayGeom* ObitOTFArrayGeomCreate (olong ndetect)
{
  ObitOTFArrayGeom* out;

  /* Create basic structure */
  out = newObitOTFArrayGeom ("Array Geometry");

  /* set members */
  out->numberDetect = ndetect;
  out->azOffset = g_realloc(out->azOffset, out->numberDetect*sizeof(ofloat));
  out->elOffset = g_realloc(out->elOffset, out->numberDetect*sizeof(ofloat));

  return out;
} /* end ObitOTFArrayGeomCreate */

/**
 * Read ObitOTFArrayGeom information from a table
 * \param in    Array geometry to update, if (*in) Null, it is created
 * \param table table to read from
 * \param err   Error stack
 * \return return code OBIT_IO_OK => OK
 */
ObitIOCode ObitOTFArrayGeomRead (ObitOTFArrayGeom **in, ObitTableOTFArrayGeom *table, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitOTFArrayGeom *out=NULL;
  ObitTableOTFArrayGeomRow *row=NULL;
  olong i, numberDetect;
  gchar *routine = "ObitOTFArrayGeomRead";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  if (*in) g_assert (ObitIsA(*in, &myClassInfo));
  g_assert (ObitTableOTFArrayGeomIsA(table));

  /* Open table */
  retCode = ObitTableOTFArrayGeomOpen (table, OBIT_IO_ReadWrite, err);
    if ((retCode != OBIT_IO_OK) || (err->error)) 
      Obit_traceback_val (err, routine, table->name, retCode);

  numberDetect = table->myDesc->nrow;
  /* Does input exist? */
  if (*in!=NULL) {
    out = *in;  /* exists */
  } else {
    /* must create */
    out = ObitOTFArrayGeomCreate (numberDetect);
  }

  /* Allocate arrays */
  out->azOffset = g_realloc(out->azOffset, numberDetect*sizeof(ofloat));
  out->elOffset = g_realloc(out->elOffset, numberDetect*sizeof(ofloat));

  /* Copy header information */
  for (i=0; i<12; i++) out->RefDate[i] = table->RefDate[i]; out->RefDate[i]=0;
  for (i=0; i<4; i++)  out->TimeSys[i] = table->TimeSys[i]; out->TimeSys[i]=0;
  out->TeleX   = table->TeleX;
  out->TeleY   = table->TeleY;
  out->TeleZ   = table->TeleZ;
  out->DegDay  = table->DegDay;
  out->GSTiat0 = table->GSTiat0;
  out->PolarX  = table->PolarX;
  out->PolarY  = table->PolarY;
  out->ut1Utc  = table->ut1Utc;
  out->dataUtc = table->dataUtc;
  out->iatUtc  = table->iatUtc;

  /* Compute some useful terms */
  /* telescope latitude in radians */
  if (fabs(out->TeleX)<1.0) out->TeleX = 1.0;
  out->lat = asin (out->TeleZ / 
		    sqrt(out->TeleX*out->TeleX + 
			 out->TeleY*out->TeleY + 
			 out->TeleZ*out->TeleZ));
  /* telescope longitude in radians */
  out->lon = atan2 (out->TeleY,  out->TeleX);
  /* LST at iat0 in radians */
  out->LSTiat0 = out->GSTiat0*DG2RAD + out->lon;
  /* Earth rotation rate in rad/day */
  out->RadDay = out->DegDay*DG2RAD;
  /* Data - IAT in days */
  out->dataIat = (out->dataUtc - out->iatUtc) / 86400.0;

  /* Create row object */
  row = newObitTableOTFArrayGeomRow (table);

  /* Loop over table */
  for (i=0; i<table->myDesc->nrow; i++) {
    /* Read table row */
    retCode = ObitTableOTFArrayGeomReadRow (table, i+1, row, err);
    if ((retCode != OBIT_IO_OK) || (err->error)) 
      Obit_traceback_val (err, routine, table->name, retCode);

    /* Save data */
    out->azOffset[row->detector-1] = row->azOff;
    out->elOffset[row->detector-1] = row->elOff;
  } /* end loop over rows */

  /* Release row object */
  row = ObitTableOTFArrayGeomRowUnref (row);

  /* Close table */
  retCode = ObitTableOTFArrayGeomClose (table, err);
  if ((retCode != OBIT_IO_OK) || (err->error)) 
    Obit_traceback_val (err, routine, table->name, retCode);
  
/* set pointer to array geometry if created here */
  if (*in==NULL) *in = out;

  return retCode;
} /* end ObitOTFArrayGeomRead */

/**
 * Write ObitOTFArrayGeom information to a table
 * \param in    Array geometry to write
 * \param table table to write to
 * \param err   Error stack
 * \return return code OBIT_IO_OK => OK
 */
ObitIOCode ObitOTFArrayGeomWrite (ObitOTFArrayGeom *in, ObitTableOTFArrayGeom *table, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  olong i, ncopy;
  ObitTableOTFArrayGeomRow *row=NULL;
  gchar *routine = "ObitOTFArrayGeomWrite";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitTableOTFArrayGeomIsA(table));

  /* Open table */
  retCode = ObitTableOTFArrayGeomOpen (table, OBIT_IO_WriteOnly, err);
  if ((retCode != OBIT_IO_OK) || (err->error)) 
    Obit_traceback_val (err, routine, table->name, retCode);

  /* Copy header information */
  ncopy = MIN (12, MAXKEYCHARTABLEOTFArrayGeom);
  for (i=0; i<ncopy; i++) table->RefDate[i] = in->RefDate[i]; table->RefDate[i]=0;
  ncopy = MIN (4, MAXKEYCHARTABLEOTFArrayGeom);
  for (i=0; i<ncopy; i++) table->TimeSys[i] = in->TimeSys[i]; table->TimeSys[i]=0;
  table->TeleX   = in->TeleX;
  table->TeleY   = in->TeleY;
  table->TeleZ   = in->TeleZ;
  table->DegDay  = in->DegDay;
  table->GSTiat0 = in->GSTiat0;
  table->PolarX  = in->PolarX;
  table->PolarY  = in->PolarY;
  table->ut1Utc  = in->ut1Utc;
  table->dataUtc = in->dataUtc;
  table->iatUtc  = in->iatUtc;

  /* Create row object */
  row = newObitTableOTFArrayGeomRow (table);

  /* Loop over table */
  for (i=0; i<in->numberDetect; i++) {
    /* Save data */
    row->detector = i+1;
    row->azOff = in->azOffset[i];
    row->elOff = in->elOffset[i];

    /* Write table row */
    retCode = ObitTableOTFArrayGeomWriteRow (table, i+1, row, err);
    if ((retCode != OBIT_IO_OK) || (err->error)) 
      Obit_traceback_val (err, routine, table->name, retCode);

  } /* end loop over rows */

  /* Release row object */
  row = ObitTableOTFArrayGeomRowUnref (row);

  /* Close table */
  retCode = ObitTableOTFArrayGeomClose (table, err);
  if ((retCode != OBIT_IO_OK) || (err->error)) 
    Obit_traceback_val (err, routine, table->name, retCode);
  
  return retCode;
} /* end ObitOTFArrayGeomWrite */

/**
 * Calculate parallactic angle.
 * \param in    Array geometry to write
 * \param time  Time (days)
 * \param ra    RA (deg)
 * \param dec   Declination (deg)
 * \return parallactic angle in degrees ).
 */
ofloat ObitOTFArrayGeomParAng (ObitOTFArrayGeom *in, ofloat time, ofloat ra, ofloat dec)
{
  ofloat decr, lst, ha, t, PA = 0.0;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* declination in radians */
  decr = dec * DG2RAD;

  /* time in IAT */
  t = time - in->dataIat;

  /* Local siderial time in radians */
  lst = in->LSTiat0  + t*in->RadDay;

  /* Hour angle in radians */
  ha = lst - ra * DG2RAD;

  PA = atan2 (cos(in->lat) * sin(ha), (sin(in->lat)*cos(decr) - cos(in->lat)*sin(decr)*cos(ha)));
  PA *= RAD2DG;
  
  return PA;
} /* end ObitOTFArrayGeomParAng */

/**
 * Calculate Elevation  angle.
 * \param in    Array geometry to write
 * \param time  Time (days)
 * \param ra    RA (deg)
 * \param dec   Declination (deg)
 * \return elevation angle in degrees ).
 */
ofloat ObitOTFArrayGeomElev (ObitOTFArrayGeom *in, ofloat time, ofloat ra, ofloat dec)
{
  ofloat decr, lst, ha, t, El = 0.0;
  odouble darg;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* declination in radians */
  decr = dec * DG2RAD;

  /* time in IAT */
  t = time - in->dataIat;

  /* Local siderial time in radians */
  lst = in->LSTiat0  + t*in->RadDay;

  /* Hour angle in radians */
  ha = lst - ra * DG2RAD;

  darg = sin (in->lat) * sin(decr) + cos (in->lat) * cos(decr) * cos (ha);
  El = (1.570796327 - acos (darg)) *  RAD2DG;

  return El;
} /* end ObitOTFArrayGeomElev */

/**
 * Return the celestial coordinates of each detector in the array.
 * \param in       Array geometry
 * \param raPoint  Telescope pointing RA (deg)
 * \param decPoint Telescope pointing Declination (deg)
 * \param rot      Array rotation (parallactic angle ) (deg)
 * \param x        Array of output RAs (deg).
 *                 allocated in calling routine, must be at least in->numberDetect.
 * \param y        Array of output Decs.
 */
void ObitOTFArrayGeomCoord(ObitOTFArrayGeom *in, ofloat raPoint, ofloat decPoint, ofloat rot,
			   ofloat *x, ofloat *y)
{
  olong i;
  ofloat crot, srot, cdec;

 /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (x!=NULL);
  g_assert (y!=NULL);

  /* First convert to RA, Dec */
  crot = cos(rot*DG2RAD);
  srot = -sin(rot*DG2RAD);

  /* Cosine dec correction */
  cdec =  cos(decPoint*DG2RAD);
  if (cdec<1.0e-5) cdec = 1.0;

  /* Rotate offset into RA, dec and add center*/
  for (i=0; i<in->numberDetect; i++) {
    x[i] = raPoint  + (crot * in->azOffset[i] - srot * in->elOffset[i])/cdec;
    y[i] = decPoint + (crot * in->elOffset[i] + srot * in->azOffset[i]);
  }
} /* end ObitOTFArrayGeomCoord */
 

/**
 * Project the locations of the array detectors onto a specified plane.
 * \param in       Array geometry
 * \param raPoint  Telescope pointing RA (deg)
 * \param decPoint Telescope pointing Declination (deg)
 * \param rot      Array rotation (parallactic angle ) (deg)
 * \param raProj   Central RA of desired projection (deg)
 * \param decProj  Central Dec of desired projection (deg)
 * \param Proj     enum defining which projection
 *                 OBIT_OTF_SIN, OBIT_OTF_ARC, OBIT_OTF_TAN
 * \param x        Array of output "X" offsets (deg) on projected plane
 *                 allocated in calling routine, must be at least in->numberDetect.
 * \param y        Array of output "Y" offsets (deg) on projected plane
 */
void ObitOTFArrayGeomProj(ObitOTFArrayGeom *in, ofloat raPoint, ofloat decPoint, ofloat rot,
			  ofloat raProj, ofloat decProj, ObitOTFProj Proj, ofloat *x, ofloat *y)
{
  olong i;
  ofloat crot, srot;
  odouble xt, yt, xxt, yyt, dt, cdec, sdec, cdec0, sdec0;
  ofloat azOff, elOff;

 /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (x!=NULL);
  g_assert (y!=NULL);

  /* First convert to RA, Dec */
  cdec0 = cos(decProj*DG2RAD);
  sdec0 = sin(decProj*DG2RAD);
  crot = cos(rot*DG2RAD);
  srot = -sin(rot*DG2RAD);

  /* Rotate offset into RA, dec */
  for (i=0; i<in->numberDetect; i++) {
    azOff = in->azOffset[i];
    /* debug - negate az
    azOff = -in->azOffset[i]; */
    elOff = in->elOffset[i];
    x[i] = crot * azOff - srot * elOff;
    y[i] = crot * elOff + srot * azOff;
  }

  /* Project by Projection type */
  switch (Proj) { 
  case OBIT_OTF_SIN:  /* -SIN */
    for (i=0; i<in->numberDetect; i++) {
      /* debug this needed for bug in IDL simulator 
	 xt = (odouble)raPoint  + x[i];
	 yt = ((odouble)decPoint + y[i]) * DG2RAD;
	 cdec = cos(yt);
	 sdec = sin(yt);
	 dt = (xt - raProj) * DG2RAD;
	 x[i] = RAD2DG*(sin(dt) * cdec);
	 y[i] = RAD2DG*(sdec * cdec0 - cdec * sdec0 * cos(dt));*/
      /* debug - this one actually more correct */
      xt = (odouble)raPoint;
      yt = ((odouble)decPoint) * DG2RAD;
      cdec = cos(yt);
      sdec = sin(yt);
      dt = (xt - raProj) * DG2RAD;
      x[i] += RAD2DG*(sin(dt) * cdec);
      y[i] += RAD2DG*(sdec * cdec0 - cdec * sdec0 * cos(dt));
    }
    break;
  case OBIT_OTF_ARC:   /* -ARC */
    for (i=0; i<in->numberDetect; i++) {
      xxt = x[i]; /* save offsets */
      yyt = y[i];
      xt = (odouble)raPoint;
      yt = ((odouble)decPoint) * DG2RAD;
      cdec = cos(yt);
      sdec = sin(yt);
      dt = (xt - raProj) * DG2RAD;
      x[i] = sin(dt) * cdec;
      y[i] = sdec * sdec0 + cdec*cdec0*cos(dt);
      y[i] = MIN (1.0, MAX(-1.0, y[i]));
      y[i] = acos(y[i]);
      if (y[i] != 0.0) 
	y[i] = y[i] / sin(y[i]);
      else
	y[i] = 1.0;
      x[i] = RAD2DG*(x[i] * y[i]);
      y[i] = RAD2DG*((sdec * cdec0 - cdec * sdec0 * cos(dt)) * y[i]);
      x[i] += xxt;  /* add offsets back in */
      y[i] += yyt;
    }
    break;
  case OBIT_OTF_TAN:   /* -TAN */
    for (i=0; i<in->numberDetect; i++) {
      xxt = x[i]; /* save offsets */
      yyt = y[i];
      xt = (odouble)raPoint;
      yt = ((odouble)decPoint) * DG2RAD;
      cdec = cos(yt);
      sdec = sin(yt);
      dt = (xt - raProj) * DG2RAD;
      x[i] = sin(dt) * cdec;
      y[i] = sdec * sdec0 + cdec * cdec0 * cos(dt);
      x[i] = RAD2DG*(x[i] / y[i]);
      y[i] = RAD2DG*((sdec * cdec0 - cdec * sdec0 * cos(dt)) / y[i]);
      x[i] += xxt;  /* add offsets back in */
      y[i] += yyt;
   }
    break;
  default:
    g_assert_not_reached(); /* unknown, barf */
  }; /* end switch to find projection */
  
} /* end ObitOTFArrayGeomProj */

/**
 * Offset the pointing position in az and el
 * \param azOff    azimuth ( xcos el) in deg.
 * \param eloff    elevation offset in deg
 * \param pa       parallactic angle (+any feed rotation)
 * \param raPoint  Telescope pointing RA (deg), returned corrected
 * \param decPoint Telescope pointing Declination (deg) returned corrected
 */
void ObitOTFArrayGeomCorrPoint(ofloat azOff, ofloat elOff, ofloat pa,
			       ofloat *raPoint, ofloat *decPoint)
{
  ofloat ra=0.0, dec=0.0, crot, srot;

  /* anything to do? */
  if ((abs(azOff)<1.0e-20) && (abs(elOff)<1.0e-20)) return;
  crot = cos(pa*DG2RAD);
  srot = sin(pa*DG2RAD);

  /* Correct pointing */
  *raPoint  = ra  + (crot * azOff - srot * elOff);
  *decPoint = dec + (crot * elOff + srot * azOff);
} /* end ObitOTFArrayGeomCorrPoint */
 
/**
 * Initialize global ClassInfo Structure.
 */
void ObitOTFArrayGeomClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitOTFArrayGeomClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitOTFArrayGeomClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitOTFArrayGeomClassInfoDefFn (gpointer inClass)
{
  ObitOTFArrayGeomClassInfo *theClass = (ObitOTFArrayGeomClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitOTFArrayGeomClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitOTFArrayGeomClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitOTFArrayGeomGetClass;
  theClass->ObitClear     = (ObitClearFP)ObitOTFArrayGeomClear;
  theClass->ObitInit      = (ObitInitFP)ObitOTFArrayGeomInit;
  theClass->newObit       = (newObitFP)newObitOTFArrayGeom;
  theClass->ObitCopy      = (ObitCopyFP)ObitOTFArrayGeomCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitOTFArrayGeomCreate = (ObitOTFArrayGeomCreateFP)ObitOTFArrayGeomCreate;

} /* end ObitOTFArrayGeomClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitOTFArrayGeomInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitOTFArrayGeom *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->thread       = newObitThread();
  in->info         = newObitInfoList(); 
  in->azOffset     = NULL;
  in->elOffset     = NULL;
  in->numberDetect = 0;

} /* end ObitOTFArrayGeomInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * For some reason this wasn't build into the GType class.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitOTFArrayGeom* cast to an Obit*.
 */
void ObitOTFArrayGeomClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitOTFArrayGeom *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  in->thread    = ObitThreadUnref(in->thread);
  in->info      = ObitInfoListUnref(in->info);
  if (in->azOffset)  g_free(in->azOffset);  in->azOffset = NULL;
  if (in->elOffset)  g_free(in->elOffset);  in->elOffset = NULL;
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitOTFArrayGeomClear */


