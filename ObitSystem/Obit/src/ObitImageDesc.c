/* $Id$  */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2009                                          */
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
/*;Correspondence about this software should be addressed as follows: */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
#include "Obit.h"
#include "ObitImageDesc.h"
#include "ObitSkyGeom.h"
#include "ObitMem.h"
#include "ObitZernike.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitImageDesc.c
 * ObitImageDesc Obit Image descriptor class definition.
 * This contains information about the observations and the coordinates
 * in the image.
 */

/*--------------- File Global Variables  ----------------*/
/** name of the class defined in this file */
static gchar *myClassName = "ObitImageDesc";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo global structure ObitIOClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitImageDescClassInfo myClassInfo = {FALSE};

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitImageDescInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitImageDescClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitImageDescClassInfoDefFn (gpointer inClass);

/*---------------Public functions---------------------------*/
/**
 * Construct Object.
 * \param name  Optional name for object, NULL = don't use
 * \return pointer to object created.
 */
ObitImageDesc* newObitImageDesc (gchar *name)
{
  ObitImageDesc* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) 
    ObitImageDescClassInit();

  /* allocate structure */
  out = ObitMemAlloc0Name(sizeof(ObitImageDesc), "ImageDesc");

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

 /* set classInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitImageDescInit((gpointer)out);

  return out;
} /* end newObitImageDesc */

/**
 * Returns ClassInfo pointer for the class.
 * Initializes class if needed on first call.
 * \return pointer to the class structure.
 */
gconstpointer ObitImageDescGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) 
    ObitImageDescClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitImageDescGetClass */

/**
 * Copy constructor.
 * The output descriptor will have the size and reference pixel
 * modified to reflect subimaging on the input, i.e. the output
 * descriptor describes the input after being read.
 * \param in Pointer to object to be copied.
 * \param out Pointer to object to be written.  
 *            If NULL then a new structure is created.
 * \param err ObitErr error stack
 * \return Pointer to new object.
 */
ObitImageDesc* ObitImageDescCopy (ObitImageDesc* in, ObitImageDesc* out, 
			  ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  gboolean oldExist;
  gchar *outName;
  olong i,j;
  gchar *routine = "ObitImageDescCopy";

  /* error checks */
  if (err->error) return out;
  Obit_retval_if_fail(ObitIsA(in, &myClassInfo), err, out,
		      "%s: Invalid input Image descriptor for %s", 
		      routine, in->name);
  if (out)  Obit_retval_if_fail(ObitIsA(out, &myClassInfo), err, out,
				"%s: Invalid output Image descriptor for %s", 
				routine, out->name);


  /* Don't bother it they are the same */
  if (in==out) return out;

  /* Create if it doesn't exist */
  oldExist = out!=NULL;
  if (!oldExist) {
    /* derive object name */
    outName = g_strconcat ("Copy: ",in->name,NULL);
    out = newObitImageDesc(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /* initialize/copy */
  out->bitpix  = in->bitpix;
  out->naxis   = in->naxis;
  out->epoch   = in->epoch;
  out->equinox = in->equinox;
  out->maxval  = in->maxval;
  out->minval  = in->minval;
  out->obsra   = in->obsra;
  out->obsdec  = in->obsdec;
  out->plane   = in->plane;
  out->row     = in->row;
  out->IOsize  = in->IOsize;
  out->areBlanks = in->areBlanks;
  out->altCrpix  = in->altCrpix;
  out->altRef    = in->altRef;
  out->restFreq  = in->restFreq;
  out->VelReference = in->VelReference;
  out->VelDef  = in->VelDef;
  out->coordType= in->coordType;
  out->xshift  = in->xshift;
  out->yshift  = in->yshift;
  out->niter   = in->niter;
  out->do3D    = in->do3D;
  out->xPxOff  = in->xPxOff;
  out->yPxOff  = in->yPxOff;
  out->beamMaj = in->beamMaj;
  out->beamMin = in->beamMin;
  out->beamPA  = in->beamPA;
  out->jlocr   = in->jlocr;
  out->jlocd   = in->jlocd;
  out->jlocs   = in->jlocs;
  out->jlocf   = in->jlocf;
  out->jlocif  = in->jlocif;
  for (i=0; i<IMLEN_VALUE; i++) out->object[i] = in->object[i];
  for (i=0; i<IMLEN_VALUE; i++) out->teles[i]  = in->teles[i];
  for (i=0; i<IMLEN_VALUE; i++) out->instrument[i] = in->instrument[i];
  for (i=0; i<IMLEN_VALUE; i++) out->observer[i]   = in->observer[i];
  for (i=0; i<IMLEN_VALUE; i++) out->obsdat[i] = in->obsdat[i];
  for (i=0; i<IMLEN_VALUE; i++) out->origin[i] = in->origin[i];
  for (i=0; i<IMLEN_VALUE; i++) out->date[i]   = in->date[i];
  for (i=0; i<IMLEN_VALUE; i++) out->bunit[i]  = in->bunit[i];

  /* loop over axes */
  for (j=0; j<IM_MAXDIM; j++) {
    out->inaxes[j] = in->inaxes[j];
    out->cdelt[j] = in->cdelt[j];
    out->crota[j] = in->crota[j];
    out->crpix[j] = in->crpix[j];
    out->crval[j] = in->crval[j];
    for (i=0; i<IMLEN_KEYWORD; i++) out->ctype[j][i] = in->ctype[j][i];
  }

  /* Free any existing info members */
  if (in->info!=NULL) {
    if (out->info) out->info = ObitInfoListUnref (out->info); 
    out->info = ObitInfoListCopy (in->info);
  }

  /* Index as well */
  ObitImageDescIndex(out);

  return out;
} /* end ObitImageDescCopy */

/**
 * Copy descriptive material (i.e. things that don't define the structure).
 * \param in  Pointer to object to be copied.
 * \param out Pointer to object to be written.  
 * \param err ObitErr error stack
 */
void ObitImageDescCopyDesc (ObitImageDesc* in, ObitImageDesc* out, 
			    ObitErr *err)
{
  olong i,j;
  gchar *routine = "ObitImageDescCopyDesc";

  /* error checks */
  if (err->error) return;
  Obit_return_if_fail(ObitIsA(in, &myClassInfo), err,
		      "%s: Invalid input Image descriptor for %s", 
		      routine, in->name);
  if (out)  Obit_return_if_fail(ObitIsA(out, &myClassInfo), err,
				"%s: Invalid output Image descriptor for %s", 
				routine, out->name);
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitIsA(out, &myClassInfo));
 
  /* initialize/copy */
  out->epoch   = in->epoch;
  out->equinox = in->equinox;
  out->obsra   = in->obsra;
  out->obsdec  = in->obsdec;
  out->altCrpix = in->altCrpix;
  out->altRef   = in->altRef;
  out->restFreq = in->restFreq;
  out->VelReference = in->VelReference;
  out->VelDef   = in->VelDef;
  out->coordType= in->coordType;
  out->xshift   = in->xshift;
  out->yshift   = in->yshift;
  out->beamMaj  = in->beamMaj;
  out->beamMin  = in->beamMin;
  out->beamPA   = in->beamPA;
  out->niter    = in->niter;
  out->do3D     = in->do3D;
  out->xPxOff   = in->xPxOff;
  out->yPxOff   = in->yPxOff;
  for (i=0; i<IMLEN_VALUE; i++) out->object[i] = in->object[i];
  for (i=0; i<IMLEN_VALUE; i++) out->teles[i]  = in->teles[i];
  for (i=0; i<IMLEN_VALUE; i++) out->instrument[i] = in->instrument[i];
  for (i=0; i<IMLEN_VALUE; i++) out->observer[i]   = in->observer[i];
  for (i=0; i<IMLEN_VALUE; i++) out->obsdat[i] = in->obsdat[i];
  for (i=0; i<IMLEN_VALUE; i++) out->origin[i] = in->origin[i];
  for (i=0; i<IMLEN_VALUE; i++) out->date[i]   = in->date[i];
  for (i=0; i<IMLEN_VALUE; i++) out->bunit[i]  = in->bunit[i];

  /* loop over axes */
  for (j=0; j<IM_MAXDIM; j++) {
    out->cdelt[j]  = in->cdelt[j];
    out->crota[j]  = in->crota[j];
    out->crpix[j]  = in->crpix[j];
    out->crval[j]  = in->crval[j];
    for (i=0; i<IMLEN_KEYWORD; i++) out->ctype[j][i] = in->ctype[j][i];
  }

  /* index output */
  ObitImageDescIndex (out);

  return;
} /* end ObitImageDescCopyDesc */

/**
 * Return default ImageDesc structure.
 * Default is a 1 x 1 floating array
 * \param name  Optional name for object, NULL = don't use
 * \return Pointer to new object.
 */
ObitImageDesc* ObitImageDescDefault (gchar *name)
{
  olong i,j;
  ObitImageDesc* out = NULL;
  gchar blank[IMLEN_VALUE];

  /* Create */
  out = newObitImageDesc(name);

  /* initialize/copy */
  out->bitpix  = -32;
  out->naxis   = 2;
  out->epoch   = 2000.0;
  out->equinox = 2000.0;
  out->maxval  = -1.0e20;
  out->minval  =  1.0e20;
  out->obsra   = 0.0;
  out->obsdec  = 0.0;
  out->plane   = 0;
  out->row     = 0;
  out->IOsize  = OBIT_IO_byPlane;
  out->areBlanks = FALSE;
  out->altCrpix  = 0.0;
  out->altRef    = 0.0;
  out->restFreq  = 0.0;
  out->VelReference = 1;
  out->VelDef  = 0;
  out->VelDef  = OBIT_Equatorial;
  out->xshift  = 0.0;
  out->yshift  = 0.0;
  out->niter   = 0;
  out->do3D    = TRUE;
  out->xPxOff  = 0.0;
  out->yPxOff  = 0.0;
  out->beamMaj = 0.0;
  out->beamMin = 0.0;
  out->beamPA  = 0.0;
  out->jlocr   = 0;
  out->jlocd   = 0;
  out->jlocs   = 0;
  out->jlocf   = 0;
  out->jlocif  = -1;
  for (i=0; i<IMLEN_VALUE; i++) blank[i] = ' ';
  for (i=0; i<IMLEN_VALUE; i++) out->object[i] = blank[i];
  for (i=0; i<IMLEN_VALUE; i++) out->teles[i]  = blank[i];
  for (i=0; i<IMLEN_VALUE; i++) out->instrument[i] = blank[i];
  for (i=0; i<IMLEN_VALUE; i++) out->observer[i]   = blank[i];
  for (i=0; i<IMLEN_VALUE; i++) out->obsdat[i] = blank[i];
  for (i=0; i<IMLEN_VALUE; i++) out->origin[i] = blank[i];
  for (i=0; i<IMLEN_VALUE; i++) out->date[i]   = blank[i];
  for (i=0; i<IMLEN_VALUE; i++) out->bunit[i]  = blank[i];

  /* loop over axes */
  for (j=0; j<IM_MAXDIM; j++) {
    out->inaxes[j] = 1;
    out->cdelt[j] = 0.0;
    out->crota[j] = 0.0;
    out->crpix[j] = 0.0;
    out->crval[j] = 0.0;
    for (i=0; i<IMLEN_KEYWORD; i++) out->ctype[j][i] = blank[i];
  }

  return out;
} /* end ObitImageDescDefault */

/**
 * Set axis order indicators.
 * \param in Pointer to object.
 */
void ObitImageDescIndex (ObitImageDesc* in)
{
  olong i;

  /* error check */
  g_assert (ObitIsA(in, &myClassInfo));

  /* initialize */
  in->jlocr  =  0;
  in->jlocd  =  1;
  in->jlocs  = -1;
  in->jlocf  = -1;
  in->jlocif = -1;

  /* loop over axes looking for labels */
  for (i=0; i<in->naxis; i++) {
    /* Ignore projection codes on celestial positions */
    if (!strncmp (in->ctype[i], "RA--",   4)) {
      in->jlocr = i;
      in->coordType = OBIT_Equatorial;
    }
    if (!strncmp (in->ctype[i], "GLON",   4)) {
      in->jlocr = i;
      in->coordType = OBIT_Galactic;
      in->epoch = in->equinox = 1950.0; /* defined in B1950 */
    }
    if (!strncmp (in->ctype[i], "ELON",   4)) {
      in->jlocr = i;
      in->coordType = OBIT_Ecliptic;
    }
    if (!strncmp (in->ctype[i], "DEC-",   4)) {
      in->jlocd = i;
      in->coordType = OBIT_Equatorial;
    }
    if (!strncmp (in->ctype[i], "GLAT",   4)) {
      in->jlocd = i;
      in->coordType = OBIT_Galactic;
      in->epoch = in->equinox = 1950.0; /* defined in B1950 */
    }
    if (!strncmp (in->ctype[i], "ELAT",   4)) {
      in->jlocd = i;
      in->coordType = OBIT_Ecliptic;
    }
    if (!strncmp (in->ctype[i], "STOKES", 6)) in->jlocs = i;
    if (!strncmp (in->ctype[i], "FREQ",     4)) in->jlocf = i;
    else if (!strncmp (in->ctype[i], "VELO",     4)) in->jlocf = i;
    else if (!strncmp (in->ctype[i], "FELO",     4)) in->jlocf = i;
    else if (!strncmp (in->ctype[i], "FREQ",     4)) in->jlocf = i;
    else if (!strncmp (in->ctype[i], "SPECLOGF", 8)) in->jlocf = i;
    else if (!strncmp (in->ctype[i], "SPECLNMF", 8)) in->jlocf = i;
    if (!strncmp (in->ctype[i], "IF",       2)) in->jlocif = i;
  }
  /* Make sure equinox set if epoch set */
  if ((in->equinox<0.1) && (in->epoch>0.1)) in->equinox = in->epoch;
} /* end ObitImageDescIndex */

/**
 * Determine the pixel location in image described by out corresponding
 * to pixel inPixel in image described by in.
 * \param in       input image descriptor
 * \param out      output image descriptor
 * \param inPixel  Pixel location in input image
 * \param outPixel [out] Pixel location in input image
 * \param err      ObitErr error stack
 * \return TRUE if outPixel is in out, else FALSE
 */
gboolean 
ObitImageDescCvtPixel(ObitImageDesc* in, ObitImageDesc* out, 
		      ofloat *inPixel, ofloat *outPixel, ObitErr *err)
{
  olong bad;
  odouble ra, dec;
  gboolean OK = FALSE;
  gchar *routine = "ObitImageDescCvtPixel";

  /* error checks */
  /*g_assert (ObitErrIsA(err));*/
  if (err->error) return OK;
  /*g_assert (ObitIsA(in, &myClassInfo));*/
  /*g_assert (ObitIsA(out, &myClassInfo));*/
  /*g_assert(inPixel!=NULL);*/
  /*g_assert(outPixel!=NULL);*/

  /* Convert input pixel to position */
  bad = 
    ObitSkyGeomWorldPos(inPixel[0], inPixel[1],
			in->crval[in->jlocr], in->crval[in->jlocd],
			in->crpix[in->jlocr], in->crpix[in->jlocd],
			in->cdelt[in->jlocr], in->cdelt[in->jlocd],
			in->crota[in->jlocd], &in->ctype[in->jlocr][4],
			&ra, &dec);
  if (bad!=0) {
    Obit_log_error(err, OBIT_Error, 
		   "%s: Error %d determining location of pixel in %s", 
		   routine, bad, in->name);
    return OK;
  }

  /* If input not Equatorial - convert */
  if (in->coordType == OBIT_Equatorial) { /* no conversion */
  } else if (in->coordType == OBIT_Galactic) {
    ObitSkyGeomGal2Eq (&ra, &dec);
  } else if (in->coordType == OBIT_Ecliptic) {
    ObitSkyGeomEc2Eq (&ra, &dec, in->epoch);
  }

  /* Does it need to be precessed? */
  if ((in->equinox==2000.0) && (out->equinox==1950.0)) {
    ObitSkyGeomJtoB (&ra, &dec);
  } else if ((in->equinox==1950.0) && (out->equinox==2000.0)) {
    ObitSkyGeomBtoJ (&ra, &dec);
  } else if (in->equinox!=1950.0 && out->equinox!=2000.0) {
    /* Oh bother, I don't know this one */
     Obit_log_error(err, OBIT_Error, 
		   "%s: cannot precess %f to %f", 
		   routine, in->equinox, out->equinox);
    return OK;
 }

  /* If output not Equatorial - convert */
  if (out->coordType == OBIT_Equatorial) { /* no conversion */
  } else if (out->coordType == OBIT_Galactic) {
    ObitSkyGeomEq2Gal (&ra, &dec);
  } else if (out->coordType == OBIT_Ecliptic) {
    ObitSkyGeomEq2Ec (&ra, &dec, out->epoch);
  }

  /* Convert position to output pixel */
  bad = 
    ObitSkyGeomXYpix(ra, dec,
		     out->crval[out->jlocr], out->crval[out->jlocd],
		     out->crpix[out->jlocr], out->crpix[out->jlocd],
		     out->cdelt[out->jlocr], out->cdelt[out->jlocd],
		     out->crota[out->jlocd], &out->ctype[out->jlocr][4],
		     &outPixel[0], &outPixel[1]);
  if (bad!=0) {
    Obit_log_error(err, OBIT_Error, 
		   "%s: Error %d determining location of pixel in %s", 
		   routine, bad, in->name);
    return OK;
  }

  /* Is outPixel in out? or close enough? */
  OK = (outPixel[0]>=0.5) && (outPixel[1]>=0.5) && 
    (outPixel[0]<=(out->inaxes[0]+0.5)) && (outPixel[1]<=(out->inaxes[1]+0.5));
  return OK;
			  
} /* end ObitImageDescCvtPixel */

/**
 * Determine the pixel location in image described by out corresponding
 * to pixel inPixel in image described by in and given a Zernike
 * model of the distortion of in.
 * \param in       input image descriptor
 * \param out      output image descriptor
 * \param nZern    Number of Zernike terms, can handle up to 17
 * \param ZCoef    Array of Zernike coefficients (piston ignored)
 * \param inPixel  Pixel location in input image
 * \param outPixel [out] Pixel location in input image
 * \param err      ObitErr error stack
 * \return TRUE if outPixel is in out, else FALSE
 */
gboolean 
ObitImageDescCvtZern(ObitImageDesc* in, ObitImageDesc* out, 
		     olong nZern, ofloat *ZCoef,
		     ofloat *inPixel, ofloat *outPixel, ObitErr *err)
{
  olong i, bad, ierr;
  odouble ra, dec;
  ofloat x, y, zx, zy, dr, dd;
  gboolean OK = FALSE;
  gchar *routine = "ObitImageDescCvtZern";

  /* error checks */
  if (err->error) return OK;

  /* Convert input pixel to position */
  bad = 
    ObitSkyGeomWorldPos(inPixel[0], inPixel[1],
			in->crval[in->jlocr], in->crval[in->jlocd],
			in->crpix[in->jlocr], in->crpix[in->jlocd],
			in->cdelt[in->jlocr], in->cdelt[in->jlocd],
			in->crota[in->jlocd], &in->ctype[in->jlocr][4],
			&ra, &dec);
  if (bad!=0) {
    Obit_log_error(err, OBIT_Error, 
		   "%s: Error %d determining location of pixel in %s", 
		   routine, bad, in->name);
    return OK;
  }

  /* If input not Equatorial - convert */
  if (in->coordType == OBIT_Equatorial) { /* no conversion */
  } else if (in->coordType == OBIT_Galactic) {
    ObitSkyGeomGal2Eq (&ra, &dec);
  } else if (in->coordType == OBIT_Ecliptic) {
    ObitSkyGeomEc2Eq (&ra, &dec, in->epoch);
  }

  /* Does it need to be precessed? */
  if ((in->equinox==2000.0) && (out->equinox==1950.0)) {
    ObitSkyGeomJtoB (&ra, &dec);
  } else if ((in->equinox==1950.0) && (out->equinox==2000.0)) {
    ObitSkyGeomBtoJ (&ra, &dec);
  } else if (in->equinox!=1950.0 && out->equinox!=2000.0) {
    /* Oh bother, I don't know this one */
     Obit_log_error(err, OBIT_Error, 
		   "%s: cannot precess %f to %f", 
		   routine, in->equinox, out->equinox);
    return OK;
 }

  /* If output not Equatorial - convert */
  if (out->coordType == OBIT_Equatorial) { /* no conversion */
  } else if (out->coordType == OBIT_Galactic) {
    ObitSkyGeomEq2Gal (&ra, &dec);
  } else if (out->coordType == OBIT_Ecliptic) {
    ObitSkyGeomEq2Ec (&ra, &dec, out->epoch);
  }

  /* Get offset from pointing position for Zernike correction */
  ObitSkyGeomShiftXY (in->crval[in->jlocr], in->crval[in->jlocd], 
		      in->crota[in->jlocd], ra, dec, 
		      &x, &y);

  /* Offset on Zernike plane */
  ObitSkyGeomRADec2Zern (in->crval[in->jlocr], in->crval[in->jlocd], 
			 x, y, &zx, &zy, &ierr);
  Obit_retval_if_fail((ierr==0), err, OK,
		      "%s: Error projecting onto Zernike Unit circle", routine);

  /* Evaluate Zernike gradient = pointing offset */
  dr = 0.0;
  dd = 0.0;
  for (i=0; i<nZern; i++) {
    dr += ZCoef[i] * ObitZernikeGradX(i+2, zx, zy);
    dd += ZCoef[i] * ObitZernikeGradY(i+2, zx, zy);
  }
  ra  -= dr;  /* apply correction */
  dec -= dd;


  /* Convert position to output pixel */
  bad = 
    ObitSkyGeomXYpix(ra, dec,
		     out->crval[out->jlocr], out->crval[out->jlocd],
		     out->crpix[out->jlocr], out->crpix[out->jlocd],
		     out->cdelt[out->jlocr], out->cdelt[out->jlocd],
		     out->crota[out->jlocd], &out->ctype[out->jlocr][4],
		     &outPixel[0], &outPixel[1]);
  if (bad!=0) {
    Obit_log_error(err, OBIT_Error, 
		   "%s: Error %d determining location of pixel in %s", 
		   routine, bad, in->name);
    return OK;
  }

  /* Is outPixel in out? or close enough? */
  OK = (outPixel[0]>=0.5) && (outPixel[1]>=0.5) && 
    (outPixel[0]<=(out->inaxes[0]+0.5)) && (outPixel[1]<=(out->inaxes[1]+0.5));
  return OK;
			  
} /* end ObitImageDescCvtZern */

/**
 * Determine the celestial coordinates of a pixel in an image.
 * \param in       input image descriptor
 * \param inPixel  Pixel location in input image
 * \param pos      [out] celestial coordinate  (RA, dec in deg)
 * \param err      ObitErr error stack
 */
void 
ObitImageDescGetPos(ObitImageDesc* in, ofloat *inPixel, 
		    odouble *pos, ObitErr *err)
{
  olong bad; 
  gchar *routine = "ObitImageDescGetPos";

 /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert(inPixel!=NULL);
  g_assert(pos!=NULL);

  /* Convert input pixel to position */
  bad = 
    ObitSkyGeomWorldPos(inPixel[0], inPixel[1],
			in->crval[in->jlocr], in->crval[in->jlocd],
			in->crpix[in->jlocr], in->crpix[in->jlocd],
			in->cdelt[in->jlocr], in->cdelt[in->jlocd],
			in->crota[in->jlocd], &in->ctype[in->jlocr][4],
			&pos[0], &pos[1]);
  if (bad!=0) {
    Obit_log_error(err, OBIT_Error, 
		   "%s: Error %d determining location of pixel in %s", 
		   routine, bad, in->name);
    return;
  }
} /* end ObitImageDescGetPos */

/**
 * Determine the pixel of a celestial coordinate in an image.
 * to pixel inPixel in image described by in.
 * \param in       input image descriptor
 * \param pos      celestial coordinate (RA, dec in deg)
 * \param outPixel [out] Pixel location in input image
 * \param err      ObitErr error stack
 */
void ObitImageDescGetPixel(ObitImageDesc* in, odouble *pos, 
			   ofloat *outPixel, ObitErr *err)

{
  olong bad; 
  gchar *routine = "ObitImageDescGetPixel";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert(pos!=NULL);
  g_assert(outPixel!=NULL);
 
  /* Convert position to output pixel */
  bad = 
    ObitSkyGeomXYpix(pos[0], pos[1],
		     in->crval[in->jlocr], in->crval[in->jlocd],
		     in->crpix[in->jlocr], in->crpix[in->jlocd],
		     in->cdelt[in->jlocr], in->cdelt[in->jlocd],
		     in->crota[in->jlocd], &in->ctype[in->jlocr][4],
		     &outPixel[0], &outPixel[1]);
  if (bad!=0) {
    Obit_log_error(err, OBIT_Error, 
		   "%s: Error %d determining location of pixel in %s", 
		   routine, bad, in->name);
    return;
  }
			  
} /* end ObitImageDescGetPixel */

/**
 * Determine if there is an overlap is the selected regions described by.
 * a pair of image descriptors
 * Test if the difference in ra and dec is less than the sum of 
 * the halfwidths.  Nonlinearities of coordinates ignored but test is 
 * somewhat generous.
 * \param in1      first input image descriptor
 * \param in2      second input image descriptor
 * \param err      ObitErr error stack
 * \return TRUE if there is overlap, else FALSE
 */
gboolean ObitImageDescOverlap(ObitImageDesc *in1, ObitImageDesc *in2, 
			      ObitErr *err)
{
  gboolean out = FALSE;
  odouble ra1, dec1, ra2, dec2, deltaRa, deltaDec;
  ofloat halfx1, halfx2, halfy1, halfy2;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return out;

  /* Positions of centers */
  ra1 = in1->crval[in1->jlocr] + 
    (0.5*in1->inaxes[in1->jlocr] + 1.0 - in1->crpix[in1->jlocr]) * 
    in1->cdelt[in1->jlocr];
  dec1 = in1->crval[in1->jlocd] + 
    (0.5*in1->inaxes[in1->jlocd] + 1.0 - in1->crpix[in1->jlocd]) * 
    in1->cdelt[in1->jlocd];
  ra2 = in2->crval[in2->jlocr] + 
    (0.5*in2->inaxes[in2->jlocr] + 1.0 - in2->crpix[in2->jlocr]) * 
    in2->cdelt[in2->jlocr];
  dec2 = in2->crval[in2->jlocd] + 
    (0.5*in2->inaxes[in2->jlocd] + 1.0 - in2->crpix[in2->jlocd]) * 
    in2->cdelt[in2->jlocd];

  /* Offsets - allow wrap in RA */
  deltaRa  = fabs (ra1-ra2);
  if (deltaRa>180.0) deltaRa  = fabs (deltaRa-360.0);
  deltaRa  = deltaRa * cos(dec1*DG2RAD);
  deltaDec = fabs (dec1-dec2);

  /* Half widths - add a little slop */
  halfx1 = fabs (0.6*in1->inaxes[in1->jlocr]*in1->cdelt[in1->jlocr]);
  halfx2 = fabs (0.6*in2->inaxes[in2->jlocr]*in2->cdelt[in2->jlocr]);
  halfy1 = fabs (0.6*in1->inaxes[in1->jlocd]*in1->cdelt[in1->jlocd]);
  halfy2 = fabs (0.6*in2->inaxes[in2->jlocd]*in2->cdelt[in2->jlocd]);

  /* Are the offsets smaller than the sum of the halfwidths?
     not terrible accurate but it doesn't need to be */
  out = (deltaRa <= (halfx1+halfx2)) && (deltaDec <= (halfy1+halfy2));

  return out;
} /* end ObitImageDescOverlap */

/**
 * Return image rotation angle on sky.
 * \param imDesc Image descriptor
 * \return rotation angle on sky (of u,v,w) in deg.
 */
ofloat ObitImageDescRotate (ObitImageDesc* imDesc)
{
  ofloat rot;

  rot = imDesc->crota[imDesc->jlocd];

  return rot;
} /* end ObitImageDescRotate */

/**
 * Return Observed (pointing) RA and Dec.
 * If observed value given in in (not both 0), these are returned,
 * else the reference position.
 * \param in     Image descriptor
 * \param [out] RAPnt  Observed RA (deg)
 * \param [out] DecPnt Observed Dec (deg)
 */
void ObitImageDescGetPoint(ObitImageDesc *in, odouble *RAPnt, odouble *DecPnt)
{
  /* Use "Observed" position if given */
  *RAPnt   = in->obsra;
  *DecPnt  = in->obsdec;
  if ((fabs(*RAPnt)<1.0e-5) && (fabs(*DecPnt)<1.0e-5)) {
    /* if zeroes - use reference position */
    *RAPnt  = in->crval[in->jlocr];
    *DecPnt = in->crval[in->jlocd];
  }
  /* Make sure RA in range [0,360] */
  if (*RAPnt<0.0)   *RAPnt += 360.0;
  if (*RAPnt>360.0) *RAPnt -= 360.0;
} /* end ObitImageDescGetPoint */

/**
 * Determine angle on sky from a position to the antenna pointing position
 * \param in  Image descriptor
 * \param x   X offset in deg from ref pixel in in
 *            This should be in the same system and equinox as in in.
 * \param y   Y offset in deg from ref pixel in in
 * \param err Error stack, returns if not empty.
 * \return angular distance on the sky (deg)
 */
ofloat ObitImageDescAngle (ObitImageDesc *in, ofloat y, ofloat x)
{
  ofloat dist = 0.0;
  odouble RAPnt, DecPnt, ra, dec, xx, yy, zz;

  /* Get pointing position */
  ObitImageDescGetPoint(in, &RAPnt, &DecPnt);
  RAPnt  *= DG2RAD;
  DecPnt *= DG2RAD;

  /* Convert offset to position */
  ObitSkyGeomXYShift (in->crval[in->jlocr], in->crval[in->jlocd],
		      x, y, ObitImageDescRotate(in), &ra, &dec);

  /* Compute distance */
  xx = DG2RAD * ra;
  yy = DG2RAD * dec;
  zz = sin (yy) * sin (DecPnt) + cos (yy) * cos (DecPnt) * cos (xx-RAPnt);
  zz = MIN (zz, 1.000);
  dist = acos (zz) * RAD2DG;
  return dist;
} /* end ObitImageDescAngle */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitImageDescClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitImageDescClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitImageDescClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitImageDescClassInfoDefFn (gpointer inClass)
{
  ObitImageDescClassInfo *theClass = (ObitImageDescClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitImageDescClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitImageDescClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitImageDescGetClass;
  theClass->newObit       = (newObitFP)newObitImageDesc;
  theClass->ObitCopy      = (ObitCopyFP)ObitImageDescCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitImageDescClear;
  theClass->ObitInit      = (ObitInitFP)ObitImageDescInit;

} /* end ObitImageDescClassDefFn */

/*---------------Private functions--------------------------*/
/**
 * Creates empty member objects, initialize reference count.
 * Does (recursive) initialization of base class members before 
 * this class.
 * \param inn Pointer to the object to initialize.
 */
void ObitImageDescInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitImageDesc *in = inn;
  olong i, j;

  /* error checks */
  g_assert (in != NULL);
  
  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->bitpix = -1;
  in->naxis = -1;
  in->maxval = -1.0e20;
  in->minval =  1.0e20;
  in->plane   = 0;
  in->row     = 0;
  in->areBlanks = FALSE;
  in->do3D      = TRUE;
  in-> xPxOff   = 0.0;
  in-> yPxOff   = 0.0;
  in->info      = newObitInfoList();
  for (i=0; i<IM_MAXDIM; i++) {
    in->inaxes[i] = 0;
    in->cdelt[i]  = 0.0;
    in->crpix[i]  = 0.0;
    in->crota[i]  = 0.0;
     for (j=0; j<IMLEN_KEYWORD; j++) in->ctype[i][j] = ' ';
  }
} /* end ObitImageDescInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 */
void ObitImageDescClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitImageDesc *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* free this class members */
  if (in->info) ObitInfoListUnref (in->info); in->info = NULL;
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);

} /* end ObitImageDescClear */



