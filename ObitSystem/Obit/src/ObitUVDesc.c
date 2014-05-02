/* $Id$      */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2014                                          */
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
#include "ObitUVDesc.h"
#include "ObitTableFQ.h"
#include "ObitTableFQUtil.h"
#include "ObitSkyGeom.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitUVDesc.c
 * ObitUVDesc Obit UV data descriptor class definition.
 * This contains information about the observations and the coordinates
 * in the image.
 */

/*--------------- File Global Variables  ----------------*/
/** name of the class defined in this file */
static gchar *myClassName = "ObitUVDesc";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo global structure ObitIOClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitUVDescClassInfo myClassInfo = {FALSE};

/** truncate double to integer precision */
#ifndef AINT  
#define AINT(x) (odouble)((olong)(x))
#endif

/** Degrees to radians factor */
#ifndef DG2RAD  
#define DG2RAD G_PI / 180.0
#endif

/**  Radians to degrees factor */
#ifndef RAD2DG  
#define RAD2DG 180.0 / G_PI
#endif
/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitUVDescInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitUVDescClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitUVDescClassInfoDefFn (gpointer inClass);

/*---------------Public functions---------------------------*/
/**
 * Construct Object.
 * \return pointer to object created.
 */
ObitUVDesc* newObitUVDesc (gchar *name)
{
  ObitUVDesc *out;

   /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitUVDescClassInit();

  /* allocate structure */
  out = g_malloc0(sizeof(ObitUVDesc));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

 /* set classInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitUVDescInit((gpointer)out);

 return out;
} /* end newObitUVDesc */

/**
 * Returns ClassInfo pointer for the class.
 * Initializes class if needed on first call.
 * \return pointer to the class structure.
 */
gconstpointer ObitUVDescGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitUVDescClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitUVDescGetClass */

/**
 * Copy constructor.
 * The output descriptor will have the size and reference pixel
 * modified to reflect subimaging on the input, i.e. the output
 * descriptor describes the input.
 * Output will have frequency arrays deallocated, if necessary use
 * ObitUVDescGetFreq to rebuild.
 * \param in Pointer to object to be copied.
 * \param out Pointer to object to be written.  
 *            If NULL then a new structure is created.
 * \param err ObitErr error stack
 * \return Pointer to new object.
 */
ObitUVDesc* ObitUVDescCopy (ObitUVDesc* in, ObitUVDesc* out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  gboolean oldExist;
  gchar *outName;
  olong i, j;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return out;
  g_assert (ObitIsA(in, &myClassInfo));
  if (out) g_assert (ObitIsA(out, &myClassInfo));
 
  /* Don't bother it they are the same */
  if (in==out) return out;

  /* Create if it doesn't exist */
  oldExist = out!=NULL;
  if (!oldExist) {
    /* derive object name */
    outName = g_strconcat ("Copy: ",in->name,NULL);
    out = newObitUVDesc(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /* initialize/copy */
  out->access  = in->access;
  out->nvis    = in->nvis;
  out->naxis   = in->naxis;
  out->nrparm  = in->nrparm;
  out->epoch   = in->epoch;
  out->equinox = in->equinox;
  out->maxBL   = in->maxBL;
  out->maxW    = in->maxW;
  out->JDObs   = in->JDObs;
  out->obsra   = in->obsra;
  out->obsdec  = in->obsdec;
  out->firstVis= in->firstVis;
  out->numVisBuff = in->numVisBuff;
  out->altCrpix = in->altCrpix;
  out->altRef   = in->altRef;
  out->restFreq = in->restFreq;
  out->VelReference = in->VelReference;
  out->VelDef   = in->VelDef;
  out->xshift   = in->xshift;
  out->yshift   = in->yshift;
  out->beamMaj  = in->beamMaj;
  out->beamMin  = in->beamMin;
  out->beamPA   = in->beamPA;
  for (i=0; i<UVLEN_VALUE; i++) out->object[i] = in->object[i];
  for (i=0; i<UVLEN_VALUE; i++) out->teles[i]  = in->teles[i];
  for (i=0; i<UVLEN_VALUE; i++) out->instrument[i] = in->instrument[i];
  for (i=0; i<UVLEN_VALUE; i++) out->observer[i]   = in->observer[i];
  for (i=0; i<UVLEN_VALUE; i++) out->obsdat[i] = in->obsdat[i];
  for (i=0; i<UVLEN_VALUE; i++) out->origin[i] = in->origin[i];
  for (i=0; i<UVLEN_VALUE; i++) out->date[i]   = in->date[i];
  for (i=0; i<UVLEN_VALUE; i++) out->bunit[i]  = in->bunit[i];
  for (i=0; i<3; i++)           out->isort[i]  = in->isort[i];

  /* loop over axes */
  for (j=0; j<UV_MAXDIM; j++) {
    out->inaxes[j] = in->inaxes[j];
    out->cdelt[j]  = in->cdelt[j];
    out->crota[j]  = in->crota[j];
    out->crpix[j]  = in->crpix[j];
    out->crval[j]  = in->crval[j];
    for (i=0; i<UVLEN_KEYWORD; i++) out->ctype[j][i] = in->ctype[j][i];
  }

  /* Random parameter labels */
  for (j=0; j<UV_MAX_RANP; j++) {
   for (i=0; i<UVLEN_KEYWORD; i++) out->ptype[j][i] = in->ptype[j][i];
  }

  /* Copy info members */
  out->info = ObitInfoListCopyData (in->info,out->info );
  
  /* Release frequency related arrays - use ObitUVDescGetFreq to rebuild */
  if (out->freqArr)  g_free(out->freqArr);  out->freqArr = NULL;
  if (out->fscale)   g_free(out->fscale);   out->fscale = NULL;
  if (out->freqIF)   g_free(out->freqIF);   out->freqIF = NULL;
  if (out->chIncIF)  g_free(out->chIncIF);  out->chIncIF = NULL;
  if (out->sideband) g_free(out->sideband); out->sideband = NULL;

  /* index output */
  ObitUVDescIndex (out);
  
  return out;
} /* end ObitUVDescCopy */

/**
 * Copy descriptive material (i.e. things that don't define the structure).
 * \param in  Pointer to object to be copied.
 * \param out Pointer to object to be written.  
 * \param err ObitErr error stack
 */
void ObitUVDescCopyDesc (ObitUVDesc* in, ObitUVDesc* out, 
			    ObitErr *err)
{
  olong i,j;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitIsA(out, &myClassInfo));
 
  /* initialize/copy */
  out->epoch   = in->epoch;
  out->equinox = in->equinox;
  out->maxBL   = in->maxBL;
  out->maxW    = in->maxW;
  out->JDObs   = in->JDObs;
  out->obsra   = in->obsra;
  out->obsdec  = in->obsdec;
  out->altCrpix = in->altCrpix;
  out->altRef   = in->altRef;
  out->restFreq = in->restFreq;
  out->VelReference = in->VelReference;
  out->VelDef   = in->VelDef;
  out->xshift   = in->xshift;
  out->yshift   = in->yshift;
  out->beamMaj  = in->beamMaj;
  out->beamMin  = in->beamMin;
  out->beamPA   = in->beamPA;
  for (i=0; i<UVLEN_VALUE; i++) out->object[i] = in->object[i];
  for (i=0; i<UVLEN_VALUE; i++) out->teles[i]  = in->teles[i];
  for (i=0; i<UVLEN_VALUE; i++) out->instrument[i] = in->instrument[i];
  for (i=0; i<UVLEN_VALUE; i++) out->observer[i]   = in->observer[i];
  for (i=0; i<UVLEN_VALUE; i++) out->obsdat[i] = in->obsdat[i];
  for (i=0; i<UVLEN_VALUE; i++) out->origin[i] = in->origin[i];
  for (i=0; i<UVLEN_VALUE; i++) out->date[i]   = in->date[i];
  for (i=0; i<UVLEN_VALUE; i++) out->bunit[i]  = in->bunit[i];
  for (i=0; i<3; i++)           out->isort[i]  = in->isort[i];

  /* loop over axes */
  for (j=0; j<UV_MAXDIM; j++) {
    out->cdelt[j]  = in->cdelt[j];
    out->crota[j]  = in->crota[j];
    out->crpix[j]  = in->crpix[j];
    out->crval[j]  = in->crval[j];
    for (i=0; i<UVLEN_KEYWORD; i++) out->ctype[j][i] = in->ctype[j][i];
  }

  /* Random parameter labels 
  for (j=0; j<UV_MAX_RANP; j++) {
   for (i=0; i<UVLEN_KEYWORD; i++) out->ptype[j][i] = in->ptype[j][i];
  }*/

  if (in->chIncIF)  g_free(in->chIncIF);  in->chIncIF = NULL;
  /* index output */
  ObitUVDescIndex (out);

  /* Copy list */
  ObitInfoListCopyData(in->info, out->info);

  return;
} /* end ObitUVDescCopyDesc */

/**
 * Copy frequency information arrays not copied by Copy
 * \param in  Pointer to object to be copied.
 * \param out Pointer to object to be written.  
 * \param err ObitErr error stack
 */
void ObitUVDescCopyFreq (ObitUVDesc* in, ObitUVDesc* out, 
			 ObitErr *err)
{
  olong nfreq, nif, size, i;

  /* error checks */
  if (err->error) return;

  /* Frequency related arrays - how big? */
  nfreq = 1;
  if (in->jlocf>=0) nfreq = MAX (1, in->inaxes[in->jlocf]);
  nif = 1;
  if (in->jlocif>=0) nif = MAX (1, in->inaxes[in->jlocif]);
  size = nfreq*nif;
  if (size<=0) return;
 
  /* Frequencies */
  if (in->freqArr) {
    if (out->freqArr)
      out->freqArr = g_realloc(out->freqArr, size*sizeof(odouble));
    else
      out->freqArr = g_malloc0(size*sizeof(odouble));
    for (i=0; i<size; i++) out->freqArr[i] = in->freqArr[i];
  }
  
  /* Frequency Scaling */
  if (in->fscale) {
    if (out->fscale)
      out->fscale = g_realloc(out->fscale, size*sizeof(odouble));
    else
      out->fscale = g_malloc0(size*sizeof(odouble));
    for (i=0; i<size; i++) out->fscale[i] = in->fscale[i];
  }
  
  /* IF frequency */
  if (in->freqIF) {
    if (out->freqIF)
      out->freqIF = g_realloc(out->freqIF, size*sizeof(odouble));
    else
      out->freqIF = g_malloc0(size*sizeof(odouble));
    for (i=0; i<nif; i++) out->freqIF[i] = in->freqIF[i];
  }
  
  /* IF channel bandwidths */
  if (in->chIncIF) {
    if (out->chIncIF)
      out->chIncIF = g_realloc(out->chIncIF, size*sizeof(ofloat));
    else
      out->chIncIF = g_malloc0(size*sizeof(ofloat));
    for (i=0; i<nif; i++) out->chIncIF[i] = in->chIncIF[i];
  }
  
  /* IF sidebands */
  if (in->sideband) {
    if (out->sideband)
      out->sideband = g_realloc(out->sideband, size*sizeof(olong));
    else
      out->sideband = g_malloc0(size*sizeof(olong));
    for (i=0; i<nif; i++) out->sideband[i] = in->sideband[i];
  }
  
  return;
} /* end ObitUVDescCopyFreq */

/**
 * Set axis, random parameter order indicators and data increments.
 * \param in Pointer to object.
 */
void ObitUVDescIndex (ObitUVDesc* in)
{
  olong i, size;

  /* error check */
  g_assert (ObitIsA(in, &myClassInfo));

  /* regular axis values */
  size=ObitUVDescRegularIndices(in);

    /* Set increments in visibility array (in units of correlations) */
    in->incs  = 1;
    in->incf  = 1;
    in->incif = 1;
    if (in->jlocs>0) { /* Product of dimensions of previous axes */
      for (i=0; i<in->jlocs; i++) in->incs *= in->inaxes[i];
    }
    if (in->jlocf>0) {
      for (i=0; i<in->jlocf; i++) in->incf *= in->inaxes[i];
    }
    if (in->jlocif>0) {
      for (i=0; i<in->jlocif; i++) in->incif *= in->inaxes[i];
    }

  /* random parameter values */
  /* initialize */
  in->ilocws = -1;
  in->ilocu  = -1;
  in->ilocv  = -1;
  in->ilocw  = -1;
  in->iloct  = -1;
  in->ilocb  = -1;
  in->ilocsu = -1;
  in->ilocfq = -1;
  in->ilocit = -1;
  in->ilocid = -1;

  /* loop over parameters looking for labels */
  for (i=0; i<in->nrparm; i++) {
    /* ignore projection codes */
    if (!strncmp (in->ptype[i], "UU-",      3)) in->ilocu  = i;
    if (!strncmp (in->ptype[i], "VV-",      3)) in->ilocv  = i;
    if (!strncmp (in->ptype[i], "WW-",      3)) in->ilocw  = i;
    if (!strncmp (in->ptype[i], "BASELINE", 8)) in->ilocb  = i;
    if (!strncmp (in->ptype[i], "TIME1",    5)) in->iloct  = i;
    if (!strncmp (in->ptype[i], "TIME",     4)) in->iloct  = i;
    if (!strncmp (in->ptype[i], "DATE",     4)) in->iloct  = i;
    if (!strncmp (in->ptype[i], "SOURCE",   6)) in->ilocsu = i;
    if (!strncmp (in->ptype[i], "FREQSEL",  7)) in->ilocfq = i;
    if (!strncmp (in->ptype[i], "INTTIM",   6)) in->ilocit = i;
    if (!strncmp (in->ptype[i], "CORR-ID",  7)) in->ilocid = i;
    if (!strncmp (in->ptype[i], "WEIGHT",   6)) in->ilocws = i;
  }

  /* number of correlations */
  in->ncorr = size / MAX (1, in->inaxes[in->jlocc]);

  /* trap for FITS compressed uv data for which the descriptors
     are for the data as pairs of shorts */
  if (in->inaxes[in->jlocc]==2) size /= (sizeof(ofloat)/sizeof(gshort));

  /* total size */
  in->lrec = in->nrparm + size;

  /* Set frequency information */
  if (in->jlocf>=0) in->freq = in->crval[in->jlocf];
  else in->freq = 1.0;

  /* Set Julian date if not there */
  if (in->JDObs<1.0) ObitUVDescDate2JD (in->obsdat, &in->JDObs);

  /* Make sure equinox set if epoch set */
  if ((in->equinox<0.1) && (in->epoch>0.1)) in->equinox = in->epoch;
} /* end  ObitUVDescIndex */

/**
 *  Find the indices correspondoning to regular parameters.  
 * \param in Pointer to object.
 * \return the overall number of axes across all regular parameters.
 */
olong ObitUVDescRegularIndices(ObitUVDesc* in)
{
  olong i, size;

  /* initialize */
  in->jlocc  = -1;
  in->jlocs  = -1;
  in->jlocf  = -1;
  in->jlocif = -1;
  in->jlocr  = -1;
  in->jlocd  = -1;

  /* loop over axes looking for labels */
  size = 1;
  for (i=0; i<in->naxis; i++) {
    size *= MAX (1, in->inaxes[i]); /* how big */
    if (!strncmp (in->ctype[i], "COMPLEX", 7)) in->jlocc  = i;
    else if (!strncmp (in->ctype[i], "STOKES",  6)) in->jlocs  = i;
    else if (!strncmp (in->ctype[i], "FREQ",    4)) in->jlocf  = i;
    else if (!strncmp (in->ctype[i], "VELO",    4)) in->jlocf = i;
    else if (!strncmp (in->ctype[i], "FELO",    4)) in->jlocf = i;
    else if (!strncmp (in->ctype[i], "IF",      2)) in->jlocif = i;
    else if (!strncmp (in->ctype[i], "RA",      2)) in->jlocr = i;
    else if (!strncmp (in->ctype[i], "GLON",    4)) in->jlocr = i;
    else if (!strncmp (in->ctype[i], "ELON",    4)) in->jlocr = i;
    else if (!strncmp (in->ctype[i], "DEC",     3)) in->jlocd = i;
    else if (!strncmp (in->ctype[i], "GLAT",    4)) in->jlocd = i;
    else if (!strncmp (in->ctype[i], "ELAT",    4)) in->jlocd = i;
    else
    {
      /* Unknown regular parameter, means trouble. */
      g_warning("Unknown regular parameter %s", in->ctype[i]);
    }
  }  

  return size;
} /* end ObitUVDescRegularIndices */

/**
 * Fills in frequency information in Descriptor from header and FQ table.
 * These are the freqArr, fscale, freqIF sideband, and chIncIF array members.
 * \param in       Descriptor to update.
 * \param fqtab    FQ table with IF frequency information
 *                 Actually an ObitTable but passed as parent class to avoid 
 *                 If NUll no FQ table available - use defaults from header
 *                 circular definitions inf include files.
 * \param SouIFOff if NonNULL, source dependent values to be added to the 
 *                 IF frequencies.  Includes selection by IF.
 * \param err      ObitErr for reporting errors.
 */
void ObitUVDescGetFreq (ObitUVDesc* in, Obit *fqtab, odouble *SouIFOff,
			ObitErr *err)
{
  olong size, i, j, fqid, nfreq;
  oint  *sideBand=NULL, nif, tnif;
  odouble *freqOff=NULL, *freq=NULL;
  ofloat *chBandw=NULL, deltaf=0.0;
  ObitIOCode retCode;
  ObitTableFQ *FQtable = (ObitTableFQ*)fqtab;
  gchar *routine = "ObitUVDescGetFreq";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* release old memory if necessary */
  if (in->freqArr)  g_free(in->freqArr);  in->freqArr = NULL;
  if (in->fscale)   g_free(in->fscale);   in->fscale = NULL;
  if (in->freqIF)   g_free(in->freqIF);   in->freqIF = NULL;
  if (in->chIncIF)  g_free(in->chIncIF);  in->chIncIF = NULL;
  if (in->sideband) g_free(in->sideband); in->sideband = NULL;

  /* how big */
  nfreq = 1;
  if (in->jlocf>=0) nfreq = MAX (1, in->inaxes[in->jlocf]);
  nif = 1;
  if (in->jlocif>=0) nif = MAX (1, in->inaxes[in->jlocif]);
  size = nfreq*nif;
  
  /* allocate */
  in->freqArr = g_malloc0(size*sizeof(odouble));
  in->fscale  = g_malloc0(size*sizeof(ofloat));
  in->freqIF  = g_malloc0(nif*sizeof(odouble));
  in->chIncIF = g_malloc0(nif*sizeof(ofloat));
  in->sideband= g_malloc0(nif*sizeof(olong));

  /* Get information from FQ table if one exists */
  if (FQtable) {
    fqid = 1; /* may need more */
    retCode = ObitTableFQGetInfo (FQtable, fqid, &tnif, &freqOff, &sideBand, 
				  &chBandw, err);
    /* Ignore missing table */
    if (retCode==OBIT_IO_OpenErr) {
      ObitErrClear(err); 
      retCode = OBIT_IO_OK;
      FQtable = NULL;
    }
  }
    if (err->error) Obit_traceback_msg (err, routine, in->name);
  /* Fake FQ info */
   if (FQtable==NULL) {
    freqOff = g_malloc0(nif*sizeof(odouble));
    chBandw = g_malloc0(nif*sizeof(ofloat));
    sideBand = g_malloc0(nif*sizeof(oint));
    for (j=0; j<nif; j++) {
      freqOff[j]  = 0.0;
      chBandw[j]  = in->cdelt[in->jlocf];
      sideBand[j] = 0;
    }
    /* Warn if more than 1 IF and likely to need this information */
    if ((nif>1) && ((in->access==OBIT_IO_ReadCal)||in->access==OBIT_IO_ReadOnly)) {
      Obit_log_error(err, OBIT_InfoWarn, "Missing Frequency (FQ) information");
    }
   }

  /* Make sure reference frequency set */
  in->freq = in->crval[in->jlocf];

  /* Source dependent OFFSets */
  if (SouIFOff!=NULL) {
    for (j=0; j<nif; j++) {
      freqOff[j] += SouIFOff[j];
    }
  }

  /* put frequencies in list in order they occur in data */
  freq = in->freqArr;
  if ((in->jlocf < in->jlocif) || (in->jlocif<=0)) {
    /* Frequency first */
    for (j=0; j<nif; j++) {
      deltaf = chBandw[j];
      if (sideBand[j]<0) deltaf = -fabs(deltaf);  /* Trap lower sideband */
      for (i=0; i<nfreq; i++) {
	*freq = freqOff[j] + in->freq + (i+1.0 - in->crpix[in->jlocf]) * deltaf;
	freq++;
      }
    }
  } else {
    /* IF first */
    for (i=0; i<nfreq; i++) {
      for (j=0; j<nif; j++) {
	if (sideBand[j]<0) deltaf = -fabs(deltaf);  /* Trap lower sideband */
	*freq = freqOff[j] + in->freq + (i+1.0 - in->crpix[in->jlocf]) * deltaf;
	freq++;
      }
    }
  }

  /* Frequency scaling factors to reference frequency */
  for (i=0; i<size; i++) {
    in->fscale[i] = in->freqArr[i] / in->freq;
  }

 /* Save IF frequency info */
  for (j=0; j<nif; j++) {
    in->freqIF[j]  = in->freq + freqOff[j];
    in->chIncIF[j] = chBandw[j];
    in->sideband[j]= sideBand[j];
  }

  /* cleanup */
  if (freqOff)  g_free(freqOff);
  if (sideBand) g_free(sideBand);
  if (chBandw)  g_free(chBandw);

 } /* end ObitUVDescGetFreq */

/**
 * Parse date string as ("yyyy-mm-dd" or "dd/mm/yy" or "yyymmdd")
 * and convert to Julian Date.
 * Algorithm from ACM Algorithm number 199
 * This routine is good from 1 Mar 1900 indefinitely.
 * \param date [in] Date string
 * \param date [out] Julian date.
 */
void ObitUVDescDate2JD (const gchar* date, odouble *JD)
{
  olong ymd[3], n;
  odouble ya, m, d, k, c;
  gchar temp[20];
  
  g_assert (date!=NULL);
  g_assert (JD!=NULL);

  /* Get date */
  if (date[2]=='/') { /* old format 'DD/MM/YY' */
    strncpy (temp, date, 8); temp[8] = 0;
    temp[2] = ' '; temp[5] = ' ';
    n = sscanf (temp, "%2d%3d%3d", &ymd[2], &ymd[1], &ymd[0]);
    /* guess century */
    if (ymd[0]>50) ymd[0]+=1900;
    else ymd[0]+=2000;
  } else if (date[4]=='-') {            /* new format 'YYYY-MM-DD' */
     strncpy (temp, date, 10); temp[10] = 0;
     temp[4] = ' '; temp[7] = ' ';
     n = sscanf (temp, "%4d%3d%3d", &ymd[0], &ymd[1], &ymd[2]);
  } else {  /* assume YYYYMMDD */
    n = sscanf (date, "%4d%2d%2d", &ymd[0], &ymd[1], &ymd[2]);
  }
  /* if something went wrong return bad date */
  if (n!=3) {*JD = -1.0; return;}
  
  /* Convert to Days */
  ya = ymd[0];
  m  = ymd[1];
  d  = ymd[2];
  if (m<=2.0) {
    m  = AINT(m+9.0);
    ya = AINT(ya-1.0);
  } else {
    m  = AINT(m-3.0);
  }
  c = AINT(ya/100.0);
  ya = ya - 100.0 * c;
  k = AINT(146097.0 * c * 0.25) +
    AINT(1461.0 * ya * 0.25) +
    AINT((153.0 * m + 2.0) * 0.2) +
    d;

  /* Following good for > 20th cent. */
  *JD = AINT(k) + 1721118.50;
} /* end ObitUVDescDate2JD */

/**
 * Convert a Julian date to a string in form "yyyy-mm-dd".
 * Apdapted from ASM Algorithm no. 199
 * \param date [in] Julian date.
 * \param date [out] Date string, Must be at least 11 characters.
 */
void ObitUVDescJD2Date (odouble JD, gchar *date)
{
  olong  id, im, iy, ic;
  odouble j, y, d, m;

  g_assert (date!=NULL);

  /* error check */
  if (JD<1.0) {
    g_snprintf (date, 11, "BAD DATE");
    return;
  }
  
  j = AINT (JD + 0.50 - 1721119.0);
  y = AINT ((4.0*j - 1.00) / 146097.0);
  ic = y + 0.00010;
  j = 4.0*j - 1.00 - 146097.0*y;
  d = AINT (j * 0.250);
  j = AINT ((4.00*d + 3.00) / 1461.00);
  d = 4.00*d + 3.00 - 1461.00*j;
  d = AINT ((d+4.00) * 0.250);
  m = AINT ((5.00*d - 3.00) / 153.00);
  id = 5.00*d - 3.00 - 153.00*m;
  id = (id + 5) / 5;
  iy = j + 100*ic;
  im = m;
  if (im < 10) {
    im = im + 3;
  } else {
    im = im - 9;
  iy = iy + 1;
  }

  /* convert to string */
  g_snprintf (date, 11, "%4.4d-%2.2d-%2.2d", iy, im, id);
} /* end ObitUVDescJD2Date */

/**
 * Determine the phase shift parameters for a position shift
 * from uv data to an image.
 * Recognizes "-SIN" and "-NCP" projections.
 * Assumes 3D imaging for -SIN and not for -NCP.
 * Apdapted from AIPS.
 * \param uvDesc UV data descriptor
 * \param imDesc Image descriptor.
 * \param dxyzc  (out) the derived shift parameters.
 * \param err    ObitErr for reporting errors.
 */
void ObitUVDescShiftPhase (ObitUVDesc* uvDesc, 
			   ObitImageDesc* imDesc, 
			   ofloat dxyzc[3], ObitErr *err)
{
  ofloat maprot, inPixel[2];
  odouble  uvra, uvdec, imra, imdec, pos[2];
  gchar *routine = "ObitUVDescShiftPhase";

  /* get positions */
  uvra  = uvDesc->crval[uvDesc->jlocr];
  uvdec = uvDesc->crval[uvDesc->jlocd];
  if (imDesc->do3D) {  /* 3D reference position is at center */
    imra  = imDesc->crval[imDesc->jlocr];
    imdec = imDesc->crval[imDesc->jlocd];
  } else { /* 2D reference is tangent - calculate position of center */
    inPixel[0] = 1.0 + imDesc->inaxes[imDesc->jlocr]*0.5;
    inPixel[1] = 1.0 + imDesc->inaxes[imDesc->jlocd]*0.5;
    ObitImageDescGetPos (imDesc, inPixel, pos, err);
    imra  = pos[0];
    imdec = pos[1];
  }
  maprot = ObitImageDescRotate(imDesc);
  
  /* which projection type? Use projection code on first image axis to decide. */
  if (!strncmp(&imDesc->ctype[0][4], "-NCP", 4)) {
    
    /* correct positions for shift since not 3D imaging ???
    uvra  += uvDesc->xshift;
    uvdec += uvDesc->yshift;
    imra  += imDesc->xshift;
    imdec += imDesc->yshift; */

    /*  -NCP projection */
    ObitSkyGeomShiftNCP (uvra, uvdec, maprot, imra, imdec, dxyzc);
    
  } else if (!strncmp(&imDesc->ctype[0][4], "-SIN", 4)) {
    
    /*  -SIN projection */
    ObitSkyGeomShiftSIN (uvra, uvdec, maprot, imra, imdec, dxyzc);
  } /* end "-SIN" projection */

  else if (!strncmp(&imDesc->ctype[0][4], "    ", 4)) {
    /* Default sine projection */
    ObitSkyGeomShiftSIN (uvra, uvdec, maprot, imra, imdec, dxyzc);

  } else {
      Obit_log_error(err, OBIT_Error, 
		     "%s: Unsupported projection %s",
		     routine,&imDesc->ctype[0][4]);
      return;
  }
  
} /* end ObitUVDescShiftPhase */

/**
 * Determine the phase shift parameters for a position shift
 * from uv data position by a given shift.
 * Recognizes "-SIN" and "-NCP" projections.
 * Assumes 3D imaging for -SIN and not for -NCP.
 * \param uvDesc  UV data descriptor
 * \param xShift  Shift from ra in deg.
 * \param yShift  Shift from dec in deg.
 * \param dxyzc   (out) the derived shift parameters.
 * \param err     ObitErr for reporting errors.
 */
void ObitUVDescShiftPosn (ObitUVDesc* uvDesc, 
			  ofloat xShift, ofloat yShift, 
			  ofloat dxyzc[3], ObitErr *err)
{
  odouble  uvra, uvdec, shra, shdec;
  gchar *routine = "ObitUVDescShiftPosn";
  ofloat uvrot;

  /* get positions - assume 3D imaging */
  uvrot  = ObitUVDescRotate(uvDesc);
  uvra  = uvDesc->crval[uvDesc->jlocr];
  uvdec = uvDesc->crval[uvDesc->jlocd];
  ObitSkyGeomXYShift (uvra, uvdec, xShift, yShift, 0.0, &shra, &shdec);
  
  /* which projection type? Use projection code on first image axis to decide. */
  if (!strncmp(&uvDesc->ptype[uvDesc->ilocu][4], "-NCP", 4)) {
    
    /* correct positions for shift since not 3D imaging */
    uvra  += uvDesc->xshift;
    uvdec += uvDesc->yshift;
    shra  += uvDesc->xshift;
    shdec += uvDesc->yshift;

    /*  -NCP projection */
    ObitSkyGeomShiftNCP (uvra, uvdec, uvrot, shra, shdec, dxyzc);
    
  } else if (!strncmp(&uvDesc->ptype[uvDesc->ilocu][4], "-SIN", 4)) {
    
    /*  -SIN projection */
    ObitSkyGeomShiftSIN (uvra, uvdec, uvrot, shra, shdec, dxyzc);
  } /* end "-SIN" projection */

  else if (!strncmp(&uvDesc->ptype[uvDesc->ilocu][4], "    ", 4)) {
    /* Default sine projection */
    ObitSkyGeomShiftSIN (uvra, uvdec, uvrot, shra, shdec, dxyzc);

  }
  else {
      Obit_log_error(err, OBIT_Error, 
		     "%s: Unsupported projection %s",
		     routine, &uvDesc->ptype[uvDesc->ilocu][4]);
      return;
  }
  
} /* end ObitUVDescShiftPosn */

/**
 * Return uvdata rotation angle on sky.
 * \param uvDesc UV data descriptor
 * \return rotation angle on sky (of u,v,w) in deg.
 */
ofloat ObitUVDescRotate (ObitUVDesc* uvDesc)
{
  ofloat rot;

  rot = uvDesc->crota[uvDesc->jlocd];

  return rot;
} /* end ObitUVDescRotate */

/**
 * Returns re-projection matrices to convert initial ref point (u,v,w)
 * and (X,Y,Z) (phase shift terms) to those at a new tangent point
 * Note that these implement the -SIN projection.  They are not needed
 * for -NCP which should not use 3D imaging.
 * If imDesc->do3D=FALSE the matrices used for 2D imaging are returned:
 * Returns re-projection matrix to convert (u,v,w) to (u',v',w') for
 * shifted gridding in the same geometry as the central facet
 * URot3D(3,1) = -L, URot3D(3,2) = -M, URot3D(3,3) = -1
 * Algorithm adapted  from AIPS PRJMAT.FOR+P2DMAT.FOR (E. W. Greisen author)
 * \param uvDesc  UV data descriptor giving initial position 
 * \param imDesc  Image descriptor giving desired position 
 * \param URot3D  [out] 3D rotation matrix for uv plane
 * \param PRot3D  [out] 3D rotation matrix for image plane
 * \return TRUE if these rotation need to be applied.
 */
gboolean ObitUVDescShift3DMatrix (ObitUVDesc *uvDesc, ObitImageDesc* imDesc,
				  ofloat URot3D[3][3], ofloat PRot3D[3][3])
{
  odouble ra, dec, xra, xdec;
  ofloat urotat, mrotat, pixel[2];
  olong  i, j, k;
  odouble rm[3][3], sa, ca, sd, cd, sd0, cd0, t[3][3], r, x[3][3], ll, mm, dd;
  gboolean isNCP, do3Dmul = FALSE;

  /* error checks */
  g_assert (ObitUVDescIsA(uvDesc));
  g_assert (ObitImageDescIsA(imDesc));
  g_assert (URot3D != NULL);
  g_assert (PRot3D != NULL);

  /* init rotation matrix */
  for (i=0; i<3; i++) 
    for (j=0; j<3; j++) 
      rm[i][j] = 0.0;
  rm[2][2] = 1.0;
  
  /* Get positions, rotations */
  ra   = uvDesc->crval[uvDesc->jlocr];
  dec  = uvDesc->crval[uvDesc->jlocd];
  xra  = imDesc->crval[imDesc->jlocr];
  xdec = imDesc->crval[imDesc->jlocd];
  mrotat = ObitImageDescRotate(imDesc);
  urotat = ObitUVDescRotate(uvDesc);
  
  /* Do these corrections need to be applied? */
  do3Dmul = (mrotat!=0.0) || (urotat!=0.0) ||
    (fabs(xra-ra) > 1.0e-10) || (fabs(xdec-dec) > 1.0e-10) ||
    (fabs(imDesc->xPxOff)>1.0e-20) || (fabs(imDesc->yPxOff)>1.0e-20);
  
  /* Don't want for NCP projection */
  isNCP = !strncmp(&imDesc->ctype[0][4], "-NCP", 4);
  do3Dmul = do3Dmul && (!isNCP);
  
  /* Stuff that depends on imaging type */
  if (imDesc->do3D) {
  
    /*  sin's and cos's */
    sa  = sin (DG2RAD * (xra-ra));
    ca  = cos (DG2RAD * (xra-ra));
    sd  = sin (DG2RAD * xdec);
    cd  = cos (DG2RAD * xdec);
    sd0 = sin (DG2RAD * dec);
    cd0 = cos (DG2RAD * dec);

    /* map + */
    r = DG2RAD * mrotat;
    rm[0][0] = cos (r);
    rm[1][1] = rm[0][0];
    rm[1][0] = -sin (r);
    rm[0][1] = sin (r);
    
    /* forward matrix */
    x[0][0] = ca;
    x[1][0] = -sd * sa;
    x[2][0] = cd * sa;
    x[0][1] = sd0 * sa;
    x[1][1] = cd * cd0 + sd * sd0 * ca;
    x[2][1] = sd * cd0 - cd * sd0 * ca;
    x[0][2] = -cd0 * sa;
    x[1][2] = cd * sd0 - sd * cd0 * ca;
    x[2][2] = sd * sd0 + cd * cd0 * ca;
    
    /* multiply */
    for (i=0; i<3; i++) {
      for (j=0; j<3; j++) {
	t[j][i] = 0.0;
	for (k=0; k<3; k++) {
	  t[j][i] += x[k][i] * rm[j][k];
	}
      }
    }
    
    
    /* uv - */
    r = DG2RAD * urotat;
    rm[0][0] = cos (r);
    rm[1][1] = rm[0][0];
    rm[1][0] = sin (r);
    rm[0][1] = -sin (r);
    
    /* multiply */
    for (i=0; i<3; i++) {
      for (j=0; j<3; j++) {
	r = 0.0;
	for (k=0; k<3; k++) {
	  r += rm[k][i] * t[j][k];
	}
	URot3D[j][i] = r;
      }
    }
    
    /* uv + */
    r = DG2RAD * urotat;
    rm[0][0] = cos (r);
    rm[1][1] = rm[0][0];
    rm[1][0] = -sin (r);
    rm[0][1] = sin (r);
    
    /* backward matrix */
    x[0][0] = ca;
    x[1][0] = sd0 * sa;
    x[2][0] = -cd0 * sa;
    x[0][1] = -sd * sa;
    x[1][1] = cd * cd0 + sd * sd0 * ca;
    x[2][1] = sd0 * cd - cd0 * sd * ca;
    x[0][2] = cd * sa;
    x[1][2] = cd0 * sd - sd0 * cd * ca;
    x[2][2] = sd * sd0 + cd * cd0 * ca;
    
    /* multiply */
    for (i=0; i<3; i++) {
      for (j=0; j<3; j++) {
	t[j][i] = 0.0;
	for (k=0; k<3; k++) {
	  t[j][i] += x[k][i] * rm[j][k];
	}
      }
    }
    
    /* map - */
    r = DG2RAD * mrotat;
    rm[0][0] = cos (r);
    rm[1][1] = rm[0][0];
    rm[1][0] = sin (r);
    rm[0][1] = -sin (r);
    
    /* multiply */
    for (i=0; i<3; i++) {
      for (j=0; j<3; j++) {
	r = 0.0;
	for (k=0; k<3; k++) {
	  r += rm[k][i] * t[j][k];
	}
	PRot3D[j][i] = r;
      }
    }
  } else {  /* 2D */
    /* Get positions */
    pixel[0] = 1.0 + imDesc->inaxes[0]*0.5;
    pixel[1] = 1.0 + imDesc->inaxes[1]*0.5;
    ObitSkyGeomWorldPos(pixel[0], pixel[1],
                        imDesc->crval[imDesc->jlocr], imDesc->crval[imDesc->jlocd],
                        imDesc->crpix[imDesc->jlocr], imDesc->crpix[imDesc->jlocd],
                        imDesc->cdelt[imDesc->jlocr], imDesc->cdelt[imDesc->jlocd],
                        imDesc->crota[imDesc->jlocd], &imDesc->ctype[imDesc->jlocr][4],
                        &xra, &xdec);
  
    /*  sin's and cos's */
    sa  = sin (DG2RAD * (xra-ra));
    ca  = cos (DG2RAD * (xra-ra));
    sd  = sin (DG2RAD * xdec);
    cd  = cos (DG2RAD * xdec);
    sd0 = sin (DG2RAD * dec);
    cd0 = cos (DG2RAD * dec);

    /* Flat coordinates */
    ll = cd * sa;
    mm = sd * cd0 - cd * sd0 * ca;
    dd = sqrt (1.0 - ll*ll - mm*mm);

    /* map + */
    r = DG2RAD * mrotat;
    rm[0][0] = cos (r);
    rm[1][1] = rm[0][0];
    rm[1][0] = -sin (r);
    rm[0][1] = sin (r);

    /* forward matrix - init */
    for (i=0; i<3; i++) 
      for (j=0; j<3; j++) 
	x[i][j] = 0.0;
    x[0][0] =  1.0;
    x[1][1] =  1.0;
    x[2][2] = -1.0;
    x[0][2] = -ll / dd;
    x[1][2] = -mm / dd;

    /* multiply */
    for (i=0; i<3; i++) {
      for (j=0; j<3; j++) {
	t[j][i] = 0.0;
	for (k=0; k<3; k++) {
	  t[j][i] += x[k][i] * rm[j][k];
	}
      }
    }

    /* uv - */
    r = DG2RAD * urotat;
    rm[0][0] = cos (r);
    rm[1][1] = rm[0][0];
    rm[1][0] = sin (r);
    rm[0][1] = -sin (r);
    
    /* multiply */
    for (i=0; i<3; i++) {
      for (j=0; j<3; j++) {
	r = 0.0;
	for (k=0; k<3; k++) {
	  r += rm[k][i] * t[j][k];
	}
	URot3D[j][i] = r;
      }
    }
    
    /* uv + */
    r = DG2RAD * urotat;
    rm[0][0] = cos (r);
    rm[1][1] = rm[0][0];
    rm[1][0] = -sin (r);
    rm[0][1] = sin (r);
    
    /*  x,y,z unchanged by PRot3D */
    for (i=0; i<3; i++) 
      for (j=0; j<3; j++) 
	x[i][j] = 0.0;
    x[0][0] = 1.0;
    x[1][1] = 1.0;
    x[2][2] = 1.0;

    /* multiply */
    for (i=0; i<3; i++) {
      for (j=0; j<3; j++) {
	t[j][i] = 0.0;
	for (k=0; k<3; k++) {
	  t[j][i] += x[k][i] * rm[j][k];
	}
      }
    }

    /* map - */
    r = DG2RAD * mrotat;
    rm[0][0] = cos (r);
    rm[1][1] = rm[0][0];
    rm[1][0] = sin (r);
    rm[0][1] = -sin (r);
    
    /* multiply */
    for (i=0; i<3; i++) {
      for (j=0; j<3; j++) {
	r = 0.0;
	for (k=0; k<3; k++) {
	  r += rm[k][i] * t[j][k];
	}
	PRot3D[j][i] = r;
      }
    }
  } /* end 2D */

  return do3Dmul; /* need to apply? */
} /* end ObitUVDescShift3DMatrix */

/**
 * Returns re-projection matrices to convert initial ref point (u,v,w)
 * and (X,Y,Z) (phase shift terms) to those at a new tangent point
 * Note that these implement the -SIN projection.  They are not needed
 * for -NCP which should not use 3D imaging.
 * Version for a given offset and rotation from the UV data
 * Algorithm adapted  from AIPS. PRJMAT.FOR (E. W. Greisen author)
 * \param uvDesc  UV data descriptor giving initial position 
 * \param shift   RA and Dec offsets(deg) from uv position
 * \param mrotat  rotation in degrees of offset image
 * \param do3D    If true use 3D Imaging, else 2D
 * \param URot3D  [out] 3D rotation matrix for uv plane
 * \param PRot3D  [out] 3D rotation matrix for image plane
 * \return TRUE if these rotation need to be applied.
 */
gboolean ObitUVDescShift3DPos (ObitUVDesc *uvDesc, ofloat shift[2], ofloat mrotat,
			       gboolean do3D, ofloat URot3D[3][3], ofloat PRot3D[3][3])
{
  odouble ra, dec, xra, xdec;
  ofloat urotat;
  olong  i, j, k;
  odouble rm[3][3], sa, ca, sd, cd, sd0, cd0, t[3][3], r, x[3][3], ll, mm, dd;
  gboolean isNCP, do3Dmul = FALSE;

  /* error checks */
  g_assert (ObitUVDescIsA(uvDesc));
  g_assert (URot3D != NULL);
  g_assert (PRot3D != NULL);
  
  /* Get position */
  ra  = uvDesc->crval[uvDesc->jlocr];
  dec = uvDesc->crval[uvDesc->jlocd];
  urotat = ObitUVDescRotate(uvDesc);
  ObitSkyGeomXYShift (ra, dec, shift[0], shift[1], mrotat, &xra, &xdec);
  
  /*  sin's and cos's */
  sa  = sin (DG2RAD * (xra-ra));
  ca  = cos (DG2RAD * (xra-ra));
  sd  = sin (DG2RAD * xdec);
  cd  = cos (DG2RAD * xdec);
  sd0 = sin (DG2RAD * dec);
  cd0 = cos (DG2RAD * dec);

  /* rotation matrix */
  for (i=0; i<3; i++) 
    for (j=0; j<3; j++) 
      rm[i][j] = 0.0;
  rm[2][2] = 1.0;

  /* Do these corrections need to be applied? */
  do3Dmul = (mrotat!=0.0) || (urotat!=0.0) ||
    (fabs(xra-ra) > 1.0e-10) || (fabs(xdec-dec) > 1.0e-10);
  
  /* Don't want for NCP projection */
  isNCP = !strncmp(&uvDesc->ctype[0][4], "-NCP", 4);
  do3Dmul = do3Dmul && (!isNCP);
    
  /* Stuff that depends on imaging type */
  if (do3D) { /* 3D */
    /* map + */
    r = DG2RAD * mrotat;
    rm[0][0] = cos (r);
    rm[1][1] = rm[0][0];
    rm[1][0] = -sin (r);
    rm[0][1] = sin (r);
    
    /* forward matrix */
    x[0][0] = ca;
    x[1][0] = -sd * sa;
    x[2][0] = cd * sa;
    x[0][1] = sd0 * sa;
    x[1][1] = cd * cd0 + sd * sd0 * ca;
    x[2][1] = sd * cd0 - cd * sd0 * ca;
    x[0][2] = -cd0 * sa;
    x[1][2] = cd * sd0 - sd * cd0 * ca;
    x[2][2] = sd * sd0 + cd * cd0 * ca;
    
    /* multiply */
    for (i=0; i<3; i++) {
      for (j=0; j<3; j++) {
	t[j][i] = 0.0;
	for (k=0; k<3; k++) {
	  t[j][i] += x[k][i] * rm[j][k];
	}
      }
    }
    
    /* uv - */
    r = DG2RAD * urotat;
    rm[0][0] = cos (r);
    rm[1][1] = rm[0][0];
    rm[1][0] = sin (r);
    rm[0][1] = -sin (r);
    
    /* multiply */
    for (i=0; i<3; i++) {
      for (j=0; j<3; j++) {
	r = 0.0;
	for (k=0; k<3; k++) {
	  r += rm[k][i] * t[j][k];
	}
	URot3D[j][i] = r;
      }
    }
    
    /* uv + */
    r = DG2RAD * urotat;
    rm[0][0] = cos (r);
    rm[1][1] = rm[0][0];
    rm[1][0] = -sin (r);
    rm[0][1] = sin (r);
    
    /* backward matrix */
    x[0][0] = ca;
    x[1][0] = sd0 * sa;
    x[2][0] = -cd0 * sa;
    x[0][1] = -sd * sa;
    x[1][1] = cd * cd0 + sd * sd0 * ca;
    x[2][1] = sd0 * cd - cd0 * sd * ca;
    x[0][2] = cd * sa;
    x[1][2] = cd0 * sd - sd0 * cd * ca;
    x[2][2] = sd * sd0 + cd * cd0 * ca;
    
    /* multiply */
    for (i=0; i<3; i++) {
      for (j=0; j<3; j++) {
	t[j][i] = 0.0;
	for (k=0; k<3; k++) {
	  t[j][i] += x[k][i] * rm[j][k];
	}
      }
    }
    
    /* map - */
    r = DG2RAD * mrotat;
    rm[0][0] = cos (r);
    rm[1][1] = rm[0][0];
    rm[1][0] = sin (r);
    rm[0][1] = -sin (r);
    
    /* multiply */
    for (i=0; i<3; i++) {
      for (j=0; j<3; j++) {
	r = 0.0;
	for (k=0; k<3; k++) {
	  r += rm[k][i] * t[j][k];
	}
	PRot3D[j][i] = r;
      }
    }
  } else { /* 2D */
    /* Flat coordinates */
    ll = cd * sa;
    mm = sd * cd0 - cd * sd0 * ca;
    dd = sqrt (1.0 - ll*ll - mm*mm);

    /* map + */
    r = DG2RAD * mrotat;
    rm[0][0] = cos (r);
    rm[1][1] = rm[0][0];
    rm[1][0] = -sin (r);
    rm[0][1] = sin (r);

    /* forward matrix - init */
    for (i=0; i<3; i++) 
      for (j=0; j<3; j++) 
	x[i][j] = 0.0;
    x[0][0] =  1.0;
    x[1][1] =  1.0;
    x[2][2] = -1.0;
    x[0][2] = -ll / dd;
    x[1][2] = -mm / dd;

    /* multiply */
    for (i=0; i<3; i++) {
      for (j=0; j<3; j++) {
	t[j][i] = 0.0;
	for (k=0; k<3; k++) {
	  t[j][i] += x[k][i] * rm[j][k];
	}
      }
    }

    /* uv - */
    r = DG2RAD * urotat;
    rm[0][0] = cos (r);
    rm[1][1] = rm[0][0];
    rm[1][0] = sin (r);
    rm[0][1] = -sin (r);
    
    /* multiply */
    for (i=0; i<3; i++) {
      for (j=0; j<3; j++) {
	r = 0.0;
	for (k=0; k<3; k++) {
	  r += rm[k][i] * t[j][k];
	}
	URot3D[j][i] = r;
      }
    }
    
    /* uv + */
    r = DG2RAD * urotat;
    rm[0][0] = cos (r);
    rm[1][1] = rm[0][0];
    rm[1][0] = -sin (r);
    rm[0][1] = sin (r);
    
    /*  x,y,x unchanged by PRot3D */
    for (i=0; i<3; i++) 
      for (j=0; j<3; j++) 
	x[i][j] = 0.0;
    x[0][0] = 1.0;
    x[1][1] = 1.0;
    x[2][2] = 1.0;

    /* multiply */
    for (i=0; i<3; i++) {
      for (j=0; j<3; j++) {
	t[j][i] = 0.0;
	for (k=0; k<3; k++) {
	  t[j][i] += x[k][i] * rm[j][k];
	}
      }
    }

    /* map - */
    r = DG2RAD * mrotat;
    rm[0][0] = cos (r);
    rm[1][1] = rm[0][0];
    rm[1][0] = sin (r);
    rm[0][1] = -sin (r);
    
    /* multiply */
    for (i=0; i<3; i++) {
      for (j=0; j<3; j++) {
	r = 0.0;
	for (k=0; k<3; k++) {
	  r += rm[k][i] * t[j][k];
	}
	PRot3D[j][i] = r;
      }
    }
  } /* end 2D */

  return do3Dmul; /* need to apply? */
} /* end ObitUVDescShift3DMatrix */

/**
 * Set number of vis per IO for a UV data set
 * Limit not to exceed 0.5 GByte
 * \param uvDesc   UV data descriptor , data should have been fully 
 *                 instantiated to get this filled in.
 * \param info     InfoList potentially with nThreads
 * \param nvis     Desired (max) number of visibilities per thread
 * \return number of visibilities per I/O
 */
olong ObitUVDescSetNVis (ObitUVDesc *uvDesc, ObitInfoList* info, olong nvis)
{
  olong nvispio;
  ollong tsize, maxsize;
  olong nThreads, lrec;
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};

  nThreads = 1;
  ObitInfoListGetTest(info, "nThreads", &type, dim, &nThreads);
  lrec = uvDesc->lrec;
  /* Sanity check */
  lrec = MAX (10, lrec);
    
  tsize = nvis * MAX (1, nThreads);                       /* How big is what is desired? */
  maxsize = 500000000 / (lrec*sizeof (ofloat));   /* How big is allowed */
  if (tsize>maxsize) {
    nvispio = (olong)(maxsize);
  } else nvispio = (olong)(tsize);
  
  nvispio = MAX (1, nvispio);   /* At least 1 */

  return nvispio;
} /* end ObitUVDescSetNVis  */


/**
 * Get data projection code
 * \param uvDesc   UV data descriptor
 * \return Projection type, one of
 *  OBIT_SkyGeom_SIN,  OBIT_SkyGeom_TAN, OBIT_SkyGeom_ARC, OBIT_SkyGeom_NCP, 
 *  OBIT_SkyGeom_GLS,  OBIT_SkyGeom_MER, OBIT_SkyGeom_AIT, OBIT_SkyGeom_STG
 */
ObitSkyGeomProj ObitUVDescGetProj (ObitUVDesc *uvDesc)
{
  if (!strncmp(&uvDesc->ptype[uvDesc->ilocu][4], "-SIN", 4)) return OBIT_SkyGeom_SIN;
  if (!strncmp(&uvDesc->ptype[uvDesc->ilocu][4], "-NCP", 4)) return OBIT_SkyGeom_NCP;
  if (!strncmp(&uvDesc->ptype[uvDesc->ilocu][4], "    ", 4)) return OBIT_SkyGeom_SIN;
  if (!strncmp(&uvDesc->ptype[uvDesc->ilocu][4], "-TAN", 4)) return OBIT_SkyGeom_TAN;
  if (!strncmp(&uvDesc->ptype[uvDesc->ilocu][4], "-ARC", 4)) return OBIT_SkyGeom_ARC;
  if (!strncmp(&uvDesc->ptype[uvDesc->ilocu][4], "-GLS", 4)) return OBIT_SkyGeom_GLS;
  if (!strncmp(&uvDesc->ptype[uvDesc->ilocu][4], "-MER", 4)) return OBIT_SkyGeom_MER;
  if (!strncmp(&uvDesc->ptype[uvDesc->ilocu][4], "-AIT", 4)) return OBIT_SkyGeom_AIT;
  if (!strncmp(&uvDesc->ptype[uvDesc->ilocu][4], "-STG", 4)) return OBIT_SkyGeom_STG;
  return OBIT_SkyGeom_SIN;  /* Default */
} /* end ObitUVDescGetProj  */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitUVDescClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitUVDescClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitUVDescClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitUVDescClassInfoDefFn (gpointer inClass)
{
  ObitUVDescClassInfo *theClass = (ObitUVDescClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitUVDescClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitUVDescClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitUVDescGetClass;
  theClass->newObit       = (newObitFP)newObitUVDesc;
  theClass->ObitCopy      = (ObitCopyFP)ObitUVDescCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitUVDescClear;
  theClass->ObitInit      = (ObitInitFP)ObitUVDescInit;

} /* end ObitUVDescClassDefFn */

/*---------------Private functions--------------------------*/
/**
 * Creates empty member objects, initialize reference count.
 * Does (recursive) initialization of base class members before 
 * this class.
 * \param inn Pointer to the object to initialize.
 */
void ObitUVDescInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitUVDesc *in = inn;

  /* error checks */
  g_assert (in != NULL);
  
  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->naxis      = -1;
  in->nrparm     = -1;
  in->firstVis   = 0;
  in->numVisBuff = 0;
  in->maxBL      = -1.0;
  in->maxW       = -1.0;
  in->maxAnt     = 0;
  in->info       = newObitInfoList();
  in->freqArr    = NULL;
  in->fscale     = NULL;
  in->freqIF     = NULL;
  in->chIncIF    = NULL;
  in->sideband   = NULL;
  in->numAnt     = NULL;
} /* end ObitUVDescInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 */
void ObitUVDescClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitUVDesc *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* free this class members */
  if (in->info) ObitInfoListUnref (in->info); in->info = NULL;
  if (in->freqArr)  g_free(in->freqArr);  in->freqArr = NULL;
  if (in->fscale)   g_free(in->fscale);   in->fscale  = NULL;
  if (in->freqIF)   g_free(in->freqIF);   in->freqIF  = NULL;
  if (in->chIncIF)  g_free(in->chIncIF);  in->chIncIF = NULL;
  if (in->sideband) g_free(in->sideband); in->sideband= NULL;
  if (in->numAnt)   g_free(in->numAnt);   in->numAnt  = NULL;
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);

} /* end ObitUVDescClear */

