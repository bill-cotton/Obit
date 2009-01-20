/* $Id$ */
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

#include "ObitPennArrayAtmFit.h"
#include "Obit2DLegendre.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitPennArrayAtmFit.c
 * ObitPennArrayAtmFit class function definitions.
 * Class with atmospheric model fitting routines for the GBT Penn Array
 * This class is derived from the Obit base class.
 */

/*--------------- File Global Variables  ----------------*/
/** name of the class defined in this file */
static gchar *myClassName = "ObitPennArrayAtmFit";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo structure ObitPennArrayAtmFitClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitPennArrayAtmFitClassInfo myClassInfo = {FALSE};

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitPennArrayAtmFitInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitPennArrayAtmFitClear (gpointer in);

/**  Private: Do actual fitting */
static void ObitPennArrayAtmFitdoFit (ObitPennArrayAtmFit *in, const ofloat *data, 
				      olong incs, ofloat *coef);

/**  Private: Filter residuals */
static gboolean 
ObitPennArrayAtmFitFilter(olong numberDetect, ofloat *resid, ofloat clip, 
			  ofloat *data, olong incs);
/** Private: Get average around median value of an array */
static ofloat medianAvg (ofloat *array, olong n);

/** Private: Get average around median value of an array with weights */
static ofloat medianAvgWt (ofloat *array, olong n);

/** Private: Set Class function pointers. */
static void ObitPennArrayAtmFitClassInfoDefFn (gpointer inClass);

/** Private: qsort ofloat comparison */
static int compare_gfloat  (const void* arg1,  const void* arg2);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitPennArrayAtmFit* newObitPennArrayAtmFit (gchar* name)
{
  ObitPennArrayAtmFit* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitPennArrayAtmFitClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitPennArrayAtmFit));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitPennArrayAtmFitInit((gpointer)out);

 return out;
} /* end newObitPennArrayAtmFit */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitPennArrayAtmFitGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitPennArrayAtmFitClassInit();

  return (gconstpointer)&myClassInfo;
} /* end  ObitPennArrayAtmFitGetClass */

/**
 * Construct object from values
 * \param name    Name for object
 * \param geom    Penn Array feed array descriptor
 * \param ncoef   Number of coefficients in polynomial model to fit for.
 * \param err     ObitError stack.
 * \return pointer to the new object.
 */
ObitPennArrayAtmFit* 
ObitPennArrayAtmFitValue (gchar *name, ObitOTFArrayGeom *geom, olong ncoef, ObitErr *err)
{
  olong size;
  olong  iTerm, iDet, iCoef;
  ofloat azMax, elMax, azNorm, elNorm;
  ObitPennArrayAtmFit *out = NULL;

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return NULL;
  g_assert (ObitOTFArrayGeomIsA(geom));

  /* Create basic structure */
  out = newObitPennArrayAtmFit(name);

  /* Attach local information */
  out->ncoef = ncoef;
  out->geom = ObitOTFArrayGeomRef(geom);

  /* Create Legendre polynomial array */
  size = ncoef * geom->numberDetect;
  out->poly = g_malloc0(size*sizeof(ofloat));

  /* Work arrays */
  out->sum2p1 = g_malloc0(ncoef*sizeof(ofloat));
  out->sum2p2 = g_malloc0(ncoef*sizeof(ofloat));
  out->tcoef  = g_malloc0(ncoef*sizeof(ofloat));
  out->resid  = g_malloc0(geom->numberDetect*sizeof(ofloat));

  /* find maximum az, el offset */
  azMax = elMax = -1.0;
  for (iDet=0; iDet<geom->numberDetect; iDet++) {
    azMax = MAX (azMax, fabs(geom->azOffset[iDet]));
    elMax = MAX (elMax, fabs(geom->elOffset[iDet]));
  }

  /* Az, el normalization factors */
  if (azMax>0.0) azNorm = 1.0 / azMax;
  else azNorm = 1.0;
  if (elMax>0.0) elNorm = 1.0 / elMax;
  else elNorm = 1.0;

  /* fill Legendre polynomial terms, Normalize az, el offsets to maximum = 1.0 */
  iTerm = 0;
  for (iDet=0; iDet<geom->numberDetect; iDet++) {
    for (iCoef=0; iCoef<ncoef; iCoef++) {
      out->poly[iTerm++] = Obit2DLegendre (iCoef, azNorm*geom->azOffset[iDet], 
					   elNorm*geom->elOffset[iDet]);
    }
  }

  return out;
} /* end ObitPennArrayAtmFitValue */

/**
 * Fit Legendre polynomial to model data if more than one term fitted.
 * average of 9 values around median value if only one term desired
 * Sky data is passed in data and a model of the atmosphere is returned in coef.
 * High residuals are clipped at the 3 sigma level
 *
 * \param in    The fitting object
 * \param data  Sky brightness values to fit, invalid values are blanked
 * \param incs  increment in data array, >1 => weights
 * \param coef  Fitted Legendre polynomial coefficients.
 *              Array should be at least size of ncoef passed to 
 *              ObitPennArrayAtmFitValue.
 *              Returned as corrections.
 */
void ObitPennArrayAtmFitFit (ObitPennArrayAtmFit *in, ofloat *data, olong incs,
			     ofloat *coef)
{
  gboolean   more = TRUE;
  ofloat clip, value, fblank = ObitMagicF();

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (data != NULL);
  g_assert (coef != NULL);

  /* Single term or polynomial */
  if (in->ncoef==1) { /* median */
    if (incs==1) {
      value = medianAvg(data, in->geom->numberDetect);
    } else { /* with weights */
      value = medianAvgWt(data, in->geom->numberDetect);
    }
    if (value==fblank) coef[0] = fblank;
    else coef[0] = -value;
  } else { /* polynomial */
    /* clip at 3 sigma */
    clip = 3.0;
    
    /* Loop fitting/editing data */
    while (more) {
      
      /* fit with current data */
      ObitPennArrayAtmFitdoFit (in, data, incs, coef);
      
      /* filter out discrepant points */
      more = ObitPennArrayAtmFitFilter(in->geom->numberDetect, in->resid, 
				       clip, data, incs);
    }
  }

} /* end ObitPennArrayAtmFitFit */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitPennArrayAtmFitClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitPennArrayAtmFitClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitPennArrayAtmFitClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitPennArrayAtmFitClassInfoDefFn (gpointer inClass)
{
  ObitPennArrayAtmFitClassInfo *theClass = (ObitPennArrayAtmFitClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitPennArrayAtmFitClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitPennArrayAtmFitClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitPennArrayAtmFitGetClass;
  theClass->ObitClear     = (ObitClearFP)ObitPennArrayAtmFitClear;
  theClass->ObitInit      = (ObitInitFP)ObitPennArrayAtmFitInit;
  theClass->newObit       = (newObitFP)newObitPennArrayAtmFit;
  theClass->ObitCopy      = NULL;
  theClass->ObitClone     = NULL;

} /* end ObitPennArrayAtmFitClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param in Pointer to the object to initialize.
 */
void ObitPennArrayAtmFitInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitPennArrayAtmFit *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->thread         = newObitThread();
  in->info           = newObitInfoList(); 
  in->geom           = NULL;
  in->poly           = NULL;
  in->sum2p1         = NULL;
  in->sum2p2         = NULL;
  in->tcoef          = NULL;
  in->resid          = NULL;
  in->ncoef          = 0;

} /* end ObitPennArrayAtmFitInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * For some reason this wasn't build into the GType class.
 * \param  in Pointer to the object to deallocate.
 *           Actually it should be an ObitPennArrayAtmFit* cast to an Obit*.
 */
void ObitPennArrayAtmFitClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitPennArrayAtmFit *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  in->thread      = ObitThreadUnref(in->thread);
  in->info        = ObitInfoListUnref(in->info);
  in->geom        = ObitOTFArrayGeomUnref(in->geom);
  if (in->poly)   g_free(in->poly); in->poly = NULL;
  if (in->sum2p1) g_free(in->sum2p1); in->sum2p1 = NULL;
  if (in->sum2p2) g_free(in->sum2p2); in->sum2p2 = NULL;
  if (in->tcoef)  g_free(in->tcoef); in->tcoef = NULL;
  if (in->resid)  g_free(in->resid); in->resid = NULL;

  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitPennArrayAtmFitClear */

/**
 * Fit Legendre polynomial to model data.
 * Sky data is passed in data and a model of the atmosphere is returned in coef.
 * Solution uses a relaxation method from Fred Schwab:
 * Pn+1 = Pn + atan2 ((dChi2/dP), (d2Chi2/dP2))
 * for each parameter P where n or n+1 indicates a given iteration,
 * dChi2/dP is the first partial derivative of Chi squared wrt P,
 * d2Chi2/d2P is the second partial derivative of Chi squared wrt P,
 * Chi2 = sum_j ((Dj - sum_i (Pi Lij))**2)
 * Dj  = Brightness of detector j
 * Pi  = Legendre coefficient i [coef(i) below],
 * Lij = Legendre term i for detector j.
 * Residuals from the model are left in in->resid.
 *
 * \param in    The fitting object
 * \param data  Sky brightness values to fit, invalid values are blanked
 * \param incs  increment in data array, >1 => weights
 * \param coef  Fitted Legendre polynomial coefficients.
 *              Array should be at least size of ncoef passed to 
 *              ObitPennArrayAtmFitValue.
 *              Returned as corrections.
 */
static void ObitPennArrayAtmFitdoFit (ObitPennArrayAtmFit *in, const ofloat *data, 
				      olong incs, ofloat *coef)
{
  olong iCoef, iDet, count, nDet;
  ofloat sum, fblank = ObitMagicF();
  olong   ncoef, iter, i, mcoef, rmscnt;
  gboolean   convgd;
  ofloat      delta, tol, test, rms, norm, model, rmsLast, res, pd1, pd2, wx;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (data != NULL);
  g_assert (coef != NULL);

  /* Initialize output */
  ncoef = in->ncoef;
  for (iCoef=0; iCoef<ncoef; iCoef++) coef[iCoef] = 0.0;
  for (iCoef=0; iCoef<ncoef; iCoef++) in->tcoef[iCoef] = 0.0;
  for (i=0; i<in->geom->numberDetect; i++) in->resid[i] = 0.0;

  /* Get zeroth order term = average */
  sum = 0.0;
  count = 0;
  nDet = in->geom->numberDetect;
  for (iDet = 0; iDet<nDet; iDet++) {
    if (data[iDet*incs]!=fblank) {
      sum += data[iDet*incs];
      count ++;
    }
  }

  /* better be some good data */
  if (count>0) {
    coef[0] = sum / count;
  } else { /* bad */
    coef[0] = fblank;
    return;
  }

  
  /* loop over iterations - 
     increase number of coefficients fitted each time up to ncoef */
  rmsLast = 1.0e20;
  for (iter= 1; iter<=100; iter++) { /* loop 600 */
    /* can decide otherwise: */
    convgd = TRUE;

    /* How many coefficients this time? */
    mcoef = MIN (ncoef, iter+1);

    /* zero sums */
    for (iCoef= 0; iCoef<mcoef; iCoef++) {
      in->sum2p1[iCoef] = 0.0;
      in->sum2p2[iCoef] = 1.0e-20;
    }

    /* loop over data doing sums */
    for (iDet = 0; iDet<nDet; iDet++) {
      if (data[iDet*incs]!=fblank) {

	/* current model */
	model = 0.0;
	for (iCoef= 0; iCoef<ncoef; iCoef++) model += coef[iCoef]*in->poly[iDet*ncoef+iCoef];
	
	/* calculate residual */
	res = data[iDet*incs] - model;
	in->resid[iDet] = res;

	/* partial derivatives for detector */
	for (iCoef= 0; iCoef<mcoef; iCoef++) {
	  pd1 = -res * in->poly[iDet*ncoef+iCoef];
	  pd2 = in->poly[iDet*ncoef+iCoef]*in->poly[iDet*ncoef+iCoef];

	  /* sum */
	  in->sum2p1[iCoef] += pd1;
	  in->sum2p2[iCoef] += pd2;
	}
      } 
    }

    /* update solutions - may have to reduce corrections and increase damping 
       until solution improves over previous iteration */
    wx = 2.0; /* damping factor */

    while (wx > 1.0e-10) {
      wx = wx * 0.5;
      
      /* convergence criterion - lower the bar each attempt  */
      tol = 5.0e-6 + iter * 1.0e-5;
      
      norm = 0.0;
      convgd = (mcoef==ncoef);
      for (iCoef= 0; iCoef<mcoef; iCoef++) { /* loop 120 */
	delta = atan2 (in->sum2p1[iCoef], in->sum2p2[iCoef]);
	test = tol * MAX (0.01, abs (coef[iCoef]));
	in->tcoef[iCoef] = coef[iCoef] - wx * delta;
	
	/* convergence test - must be fitting for all desired coefficients */
	convgd = convgd  && (fabs(delta) <= test);
	norm = norm + delta*delta;
      } /* end loop updating test coefficients */
      
      /* Determine rms */
      rms = 0.0;
      rmscnt = 0;
      for (iDet= 0; iDet<nDet; iDet++) {
	if (data[iDet*incs]!=fblank) {
	  /* current model */
	  model = 0.0;
	  for (iCoef= 0; iCoef<ncoef; iCoef++) model += in->tcoef[iCoef]*in->poly[iDet*ncoef+iCoef];
	  
	  /* calculate residual */
	  res = data[iDet*incs] - model;
	  in->resid[iDet] = res;
	  
	  /* residual statistics */
	  rms += res*res;
	  rmscnt++;
	} 
      }
      
      /* Force residuals to decrease */
      if (rmscnt > 0) {
	rms = sqrt (rms/rmscnt);
      } else {
	rms = -1.0;
      } 
      
      /* Are we there yet? */
      if (rms < rmsLast) break;
    } /* end loop decreasing change in solution */

    /* Remember this rms for next iteration */
    rmsLast = rms;

    /* Save values */
    for (iCoef= 0; iCoef<ncoef; iCoef++) coef[iCoef] = in->tcoef[iCoef];
    
    /* converged? */
    if (convgd || (wx < 1.0e-10)) {
      /* Convert to corrections */
     for (iCoef= 0; iCoef<ncoef; iCoef++) coef[iCoef] = -coef[iCoef];
     return;
    }
    
  } /* end of iteration loop */

  /* Convert to corrections */
  for (iCoef= 0; iCoef<ncoef; iCoef++) coef[iCoef] = -coef[iCoef];
  
} /* end ObitPennArrayAtmFitFit */

/**
 * Filter data based on residuals to fit
 * \param numberDetect  Number of elements
 * \param resid         Residuals to model, invalid values are blanked
 * \param clip          Clip level in sigma
 * \param data          Sky brightness values, flagged values blanked
 * \param incs  increment in data array, >1 => weights
 * \return TRUE if some data points removed.
 */
static gboolean 
ObitPennArrayAtmFitFilter(olong numberDetect, ofloat *resid, ofloat clip, 
			  ofloat *data, olong incs)
{
  ofloat sum, rms,fblank = ObitMagicF();
  olong  i, rmscnt;
  gboolean out = FALSE;

  /* get RMS about 0 */
  sum = 0.0; rmscnt = 0;
  for (i=0; i<numberDetect; i++) {
    if (resid[i]!=fblank) {
      sum += (resid[i]) * (resid[i]);
      rmscnt++;
    }
  }
  
  /* need enough data */
  if (rmscnt<12) return out;
  rms = sqrt(sum/rmscnt);
  clip*=rms;

  /* do clipping */
  for (i=0; i<numberDetect; i++) {
    if ((resid[i]!=fblank) && (fabs(resid[i])>clip))  { /* clip? */
      data[i*incs] = fblank;
      if (incs>1) data[i*incs+1] = 0.0;
      out = TRUE;
    }
  }
  

  return out;
} /* end ObitPennArrayAtmFitFilter */

/**
 * Return average of central 50% of values around median
 * Does qsort then returns average
 * \param array array of values, on return will be in ascending order
 * \param n     dimension of array
 * \return      Median/average
 */
static ofloat medianAvg (ofloat *array, olong n)
{
  olong i, cnt, wid, ngood, ind, lo, hi;
  ofloat temp, out, center, fblank = ObitMagicF();

  /* Set blanked data to large value */
  ngood = 0;
  for (i=0; i<n; i++) {
    if (array[i]==fblank) array[i] = 2.0e20;
    else ngood++;
  }
  if (ngood<1) return fblank;  /* Any good data? */

  /* Sort to ascending order */
  qsort ((void*)array, n, sizeof(ofloat), compare_gfloat);

  /* How many good values? */
  cnt = 0;
  for (i=0; i<n; i++) if (array[i*2]<1.0e20) cnt++;
  wid = cnt/4;
  ngood = cnt;
  if (ngood<1) return fblank;

  ind = ngood/2 - 1;
  out = array[ind];

  /* average central 50% values */
  if (ngood>5) {
    temp = 0.0; cnt = 0;
    center = (ngood-1)/2.0;
    lo = (olong)(center - wid + 0.6);
    hi = (olong)(center + wid);
    for (i=lo; i<=hi; i++) {
      if (array[i]<1.0e20) {
	cnt++;
	temp += array[i];
      }
    }
    if (cnt>0) out = temp/cnt;
  }

  /* reblank values */
  for (i=0; i<n; i++) {
    if (array[i]>1.0e20) array[i] = fblank;
  }

  return out;
} /* end medianAvg */

/**
 * Return average of central 50% of  values around median
 * Does qsort then returns average of  Center
 * \param array array of values, on return will be in ascending order
 * \param n     dimension of array/2
 * \return      Median/average
 */
static ofloat medianAvgWt (ofloat *array, olong n)
{
  olong i, wid, cnt, ngood, ind, hi, lo;
  ofloat temp, wt, out, fblank = ObitMagicF();
  ofloat maxWt, clipWt, center;

  /* Clip to top 95% and set to 2.0e20 */
  ngood = 0;
  maxWt = -1.0e20;
  for (i=0; i<n; i++) if (array[i*2]!=fblank) maxWt = MAX(maxWt, array[i*2+1]);
  /* Clip weights and set to large value */
  clipWt = 0.05*maxWt;
  for (i=0; i<n; i++) {
    if ((array[i*2]==fblank) || (array[i*2+1]<clipWt)) {
      array[i*2]   = 2.0e20;
      array[i*2+1] = 0.0;
    } else ngood++;
  }
  if (ngood<1) return fblank;  /* Any good data? */

  /* Sort to ascending order */
  qsort ((void*)array, n, 2*sizeof(ofloat), compare_gfloat);

  /* How many good values? */
  cnt = 0;
  for (i=0; i<n; i++) if (array[i*2]<1.0e20) cnt++;
  wid = cnt/4;
  ngood = cnt;

  ind = (ngood/2)*2;  /* make center pixel even to get value not wt */
  ind -= 2;           /* 0 rel */
  out = array[ind];

  /* average central 50% values */
  if (n>5) {
    temp = 0.0; wt = 0.0;
    center = (ngood-1)/2.0;
    lo = (olong)(center - wid + 0.6);
    hi = (olong)(center + wid);
    ind = ngood/2 - 1;
    for (i=lo; i<=hi; i++) {
      if (array[i*2]<1.0e20) {
	wt += array[i*2+1];
	temp += array[i*2]*array[i*2+1];
      }
    }
    if (wt>0.0) out = temp/wt;
  }

  /* reblank values */
  for (i=0; i<n; i++) {
    if (array[i*2]>1.0e20) array[i*2] = fblank;
  }

  return out;
} /* end medianAvgWt */

/**
 * ofloat comparison of two arguments
 * \param arg1 first value to compare
 * \param arg2 second value to compare
 * \return negative if arg1 is less than arg2, zero if equal
 *  and positive if arg1 is greater than arg2.
 */
static int compare_gfloat  (const void* arg1,  const void* arg2)
{
  int out = 0;
  ofloat larg1, larg2;

  larg1 = *(ofloat*)arg1;
  larg2 = *(ofloat*)arg2;
  if (larg1<larg2) out = -1;
  else if (larg1>larg2) out = 1;
  return out;
} /* end compare_gfloat */
