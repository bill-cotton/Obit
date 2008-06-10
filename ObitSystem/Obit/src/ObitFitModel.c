/* $Id: ObitFitModel.c,v 1.3 2007/10/11 13:35:59 bcotton Exp $        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2006                                               */
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

#include "ObitFitModel.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitFitModel.c
 * ObitFitModel class function definitions.
 * This class is derived from the Obit base class.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitFitModel";

/** Function to obtain parent ClassInfo */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo structure ObitFitModelClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitFitModelClassInfo myClassInfo = {FALSE};

/*--------------- File Global Variables  ----------------*/


/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitFitModelInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitFitModelClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitFitModelClassInfoDefFn (gpointer inClass);

/* Private: Deconvolve Gaussian  */
static olong deconv (ofloat fmaj, ofloat fmin, ofloat fpa, 
		    ofloat cmaj, ofloat cmin, ofloat cpa, 
		    ofloat *rmaj, ofloat *rmin, ofloat *rpa);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitFitModel* newObitFitModel (gchar* name)
{
  ObitFitModel* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitFitModelClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitFitModel));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitFitModelInit((gpointer)out);

 return out;
} /* end newObitFitModel */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitFitModelGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitFitModelClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitFitModelGetClass */

/**
 * Make a deep copy of an ObitFitModel.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitFitModel* ObitFitModelCopy  (ObitFitModel *in, ObitFitModel *out, 
				 ObitErr *err)
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
    out = newObitFitModel(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  if (out->parms)
    out->parms = g_realloc(out->parms, in->nparm*sizeof(ofloat));
  else
    out->parms = g_malloc0(in->nparm*sizeof(ofloat));
  for (i=0; i<in->nparm; i++) out->parms[i] = in->parms[i];
  if (out->eparms)
    out->eparms = g_realloc(out->eparms, in->nparm*sizeof(ofloat));
  else
    out->eparms = g_malloc0(in->nparm*sizeof(ofloat));
  for (i=0; i<in->nparm; i++) out->eparms[i] = in->eparms[i];
  out->nparm  = in->nparm;
  out->type   = in->type;
  out->Peak   = in->Peak;
  out->DeltaX = in->DeltaX;
  out->DeltaY = in->DeltaY;

  return out;
} /* end ObitFitModelCopy */

/**
 * Make a copy of a object but do not copy the actual data
 * This is useful to create an FitModel similar to the input one.
 * \param in  The object to copy
 * \param out An existing object pointer for output, must be defined.
 * \param err Obit error/message stack object.
 */
void ObitFitModelClone  (ObitFitModel *in, ObitFitModel *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitIsA(out, &myClassInfo));

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  out->parms = g_realloc(out->parms, in->nparm*sizeof(ofloat));
  out->nparm = in->nparm;
  out->eparms = g_realloc(out->eparms, in->nparm*sizeof(ofloat));

} /* end ObitFitModelClone */

/**
 * Creates an ObitFitModel 
 * \param name   An optional name for the object.
 * \param type   Model type of the model component
 * \param Peak   Peak density
 * \param DeltaX "X" (RA) offset (deg) of center from reference position
 * \param DeltaY "Y" (Dec) offset (deg) of center from reference position
 * \param nparm  Number of parameters
 * \param parms  Model parameters, type dependent 
 * \return the new object.
 */
ObitFitModel* ObitFitModelCreate (gchar* name, ObitFitModelCompType type, 
				  ofloat Peak, ofloat DeltaX, ofloat DeltaY, 
				  olong nparm, ofloat *parms)
{
  ObitFitModel* out;
  olong i;

  /* Create basic structure */
  out = newObitFitModel (name);

  /* Set values */
  out->parms = g_malloc(nparm*sizeof(ofloat));
  for (i=0; i<nparm; i++) out->parms[i] = parms[i];
  out->eparms = g_malloc(nparm*sizeof(ofloat));
  out->nparm  = nparm;
  out->type   = type;
  out->Peak   = Peak;
  out->DeltaX = DeltaX;
  out->DeltaY = DeltaY;

  /* initialize errors to -1.0 */
  for (i=0; i<nparm; i++) out->eparms[i] =-1.0;
  out->ePeak   = -1.0;
  out->eDeltaX = -1.0;
  out->eDeltaY = -1.0;

  return out;
} /* end ObitFitModelCreate */

/**
 * Subroutine BMVAL deconvoles the fitted beam from the clean beam and  
 * also generates appropriate errors.  
 * Routine translated from the AIPSish VSAD.FOR/BMVAL  
 * \param bmaj    Fitted major axis 
 * \param bmin    Fitted minor axis 
 * \param bpa     Fitted pos. angle (deg) 
 * \param ebmaj   Fitted major axis error 
 * \param ebmin   Fitted minor axis error 
 * \param ebpa    Fitted pos. angle error (deg) 
 * \param cbmaj   Clean beam major axis 
 * \param cbmin   Clean beam minor axis 
 * \param cbpa    Clean beam pos. angle (deg) 
 * \param dgau    Deconvolved (Major, minor, PA) array 
 *                dgau[0][*] = deconvolved
 *                dgau[1][*] = upper bound
 *                dgau[2][*] = lower bound
 * \return   Error return 0-> Can completely deconvolve 
 *           1-> Cannot deconv some limits 
 *           2-> Cannot deconv fitted source 
 */
olong ObitFitModelDeconGau (ofloat bmaj, ofloat bmin, ofloat bpa, 
			   ofloat ebmaj, ofloat ebmin, ofloat ebpa,
			   ofloat cbmaj, ofloat cbmin, ofloat cbpa,
			   ofloat dgau[3][3])
{
  olong   ier, ierr, i, j, k, ic;
  ofloat  b1, b2, b3, delt[3]={-0.7,0.0,0.7}, major, minor, pa, temp;
  
  /* Deconvolve the fit */
  ier = 0;
  dgau[1][0] = 1.0e20;
  dgau[1][1] = 1.0e20;
  dgau[1][2] = 1.0e20;
  dgau[2][0] = -1.0e20;
  dgau[2][1] = -1.0e20;
  dgau[2][2] = -1.0e20;
  ierr = deconv (bmaj, bmin, bpa, cbmaj, cbmin, cbpa, &major, &minor, &pa);

  /* Could not deconvolve */
  if (ierr != 0) {
    ier = 2;
    dgau[0][0] = 0.;
    dgau[1][0] = 0.;
    dgau[2][0] = major;
    dgau[0][1] = 0.;
    dgau[1][1] = 0.;
    dgau[2][1] = major;
    dgau[0][2] = 0.;
    dgau[1][2] = 0.;
    dgau[2][2] = 180.00;

    /* Put in deconvolved size */
  } else {
    dgau[0][0] = major;
    dgau[0][1] = minor;
    dgau[0][2] = fmodf ((pa+720.0), 180.0);
  } 
  
  /* Set up looping */
  ic = 0;
  for (k= 1; k<=3; k++) { /* loop 50 */
    b3 = bpa + delt[k-1] * ebpa;
    for (j= 1; j<=3; j++) { /* loop 50 */
      b2 = bmin + delt[j-1] * ebmin;
      for (i= 1; i<=3; i++) { /* loop 50 */
	b1 = bmaj + delt[i-1] * ebmaj;
	ic = ic + 1;
	ierr = deconv (b1, b2, b3, cbmaj, cbmin, cbpa, &major, &minor, &pa);
	if (ier > 0) goto L40;
	if ((ier == 0)  &&  (ierr == 0)) goto L22;

	/* Hit an impossible deconvolution */
	dgau[1][0] = 0.;
	dgau[1][1] = 0.;
	dgau[1][2] = 0.;
	dgau[2][2] = 180.000;
	ier = 1;
	goto L40;

	/* Get upper and lower bounds. but first look at PA */
      L22:
	pa = fmodf((dgau[0][2] - pa +720.0), 180.0);
	if (pa < 45) goto L24;
	if (pa > 135) goto L26;
	
	/* Switch major, minor axes */
	temp = minor;
	minor = major;
	major = temp;
	pa = pa + dgau[0][2] -90.0;
	goto L30;

      L24:
	pa = pa + dgau[0][2];
	goto L30;

      L26:
	pa = pa + dgau[0][2] - 180.0;
	
	/* Upper and lower bounds */
      L30:            
	dgau[1][0] = MIN (dgau[1][0], major);
	dgau[2][0] = MAX (dgau[2][0], major);
	dgau[1][1] = MIN (dgau[1][1], minor);
	dgau[2][1] = MAX (dgau[2][1], minor);
	dgau[1][2] = MIN (dgau[1][2], pa);
	dgau[2][2] = MAX (dgau[2][2], pa);
	continue;
	/* No conv. look for max */
      L40:            
	dgau[2][0] = MAX (dgau[2][0], major);
	dgau[2][1] = dgau[2][0];
      } /* end loop  L50:  */
    } /* end loop  L50:  */
  } /* end loop  L50:  */
  return ier;
} /* end ObitFitModelDeconGau */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitFitModelClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitFitModelClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitFitModelClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitFitModelClassInfoDefFn (gpointer inClass)
{
  ObitFitModelClassInfo *theClass = (ObitFitModelClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitFitModelClassInit;
  theClass->newObit       = (newObitFP)newObitFitModel;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitFitModelClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitFitModelGetClass;
  theClass->ObitCopy      = (ObitCopyFP)ObitFitModelCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitFitModelClear;
  theClass->ObitInit      = (ObitInitFP)ObitFitModelInit;
  theClass->ObitFitModelCreate = (ObitFitModelCreateFP)ObitFitModelCreate;

} /* end ObitFitModelClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitFitModelInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitFitModel *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->parms  = NULL;
  in->eparms = NULL;
  in->nparm  = 0;
  in->type   = OBIT_FitModel_Unknown;
   } /* end ObitFitModelInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitFitModel* cast to an Obit*.
 */
void ObitFitModelClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitFitModel *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  if (in->parms)  g_free (in->parms);
  if (in->eparms) g_free (in->eparms);
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
   if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
     ParentClass->ObitClear (inn);
  
} /* end ObitFitModelClear */

/**
 * Deconvolves a Gaussian "beam" from a gaussian component.  
 * Routine translated from the AIPSish DECONV.FOR/DECONV  
 * \param fmaj   Fitted major axis 
 * \param fmin   Fitted minor axis 
 * \param fpa    Fitted position angle of major axis 
 * \param cmaj   Point source major axis 
 * \param cmin   Point source minor axis 
 * \param cpa    Point source position angle of major axis 
 * \param rmaj   [out] Real major axis; = 0 => unable to fit 
 * \param rmin   [out] Real minor axis; = 0 => unable to fit 
 * \param rpa    [out] Real position angle of major axis 
 * \return Error return: 0 => ok 
 *         1,2-> # components unable to deconvolve 
 */
static olong deconv (ofloat fmaj, ofloat fmin, ofloat fpa, 
		    ofloat cmaj, ofloat cmin, ofloat cpa, 
		    ofloat *rmaj, ofloat *rmin, ofloat *rpa)
{
  olong   ierr=0;
  ofloat      cmj2, cmn2, fmj2, fmn2, sinc, cosc, rhoc, sigic2, 
    det, rhoa, lfpa, lcpa, konst = 28.647888;

  /* Get useful constants */
  lfpa = fmodf (fpa+900.0, 180.0);
  lcpa = fmodf (cpa+900.0, 180.0);
  cmj2 = cmaj * cmaj;
  cmn2 = cmin * cmin;
  fmj2 = fmaj * fmaj;
  fmn2 = fmin * fmin;
  sinc = (lfpa - lcpa) / konst;
  cosc = cos (sinc);
  sinc = sin (sinc);

  /* Trigonometry now */
  rhoc = (fmj2 - fmn2) * cosc - (cmj2 - cmn2);
  if (rhoc == 0.0) {
    sigic2 = 0.0;
    rhoa = 0.0;
  } else {
    sigic2 = atan((fmj2 - fmn2) * sinc / rhoc);
    rhoa = ((cmj2 - cmn2) - (fmj2 - fmn2) * cosc) / (2.0 * cos (sigic2));
  } 

  (*rpa) = sigic2 * konst + lcpa;
  det = ((fmj2 + fmn2) -(cmj2 + cmn2)) / 2.0;
  (*rmaj) = det - rhoa;
  (*rmin) = det + rhoa;
  ierr = 0;
  if (*rmaj < 0.0) ierr++;
  if (*rmin < 0.0) ierr++;

  /* Swap to get major > minor */
  (*rmaj) = MAX (0.0, *rmaj);
  (*rmin) = MAX (0.0, *rmin);
  (*rmaj) = sqrt (fabs (*rmaj));
  (*rmin) = sqrt (fabs (*rmin));
  if (*rmaj < *rmin) {
    sinc = (*rmaj);
    (*rmaj) = (*rmin);
    (*rmin) = sinc;
    (*rpa) = (*rpa)+90.0;
  } 

  /* Fix up PA */
  (*rpa) = fmodf (*rpa+900.0, 180.0);
  if (*rmaj == 0.0) {
    (*rpa) = 0.0;
  } else if (*rmin == 0.0) {
    if ((fabs(*rpa-lfpa) > 45.0)  &&  (fabs(*rpa-lfpa) < 135.0)) 
      (*rpa) = fmodf (*rpa+450.0, 180.0);
  } 
  return ierr;
} /* end of routine deconv */ 
