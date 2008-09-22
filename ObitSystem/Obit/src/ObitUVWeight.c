/* $Id$    */
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

#include "ObitUVWeight.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitUVWeight.c
 * ObitUVWeight class function definitions.
 * This class is derived from the Obit base class.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitUVWeight";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/** Degrees to radians factor */
#ifndef DG2RAD  
#define DG2RAD G_PI / 180.0
#endif

/**  Radians to degrees factor */
#ifndef RAD2DG  
#define RAD2DG 180.0 / G_PI
#endif

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitUVWeight
 * returns a ObitUVWeight*.
 * in = object to unreference
 */
#define ObitUVWeightUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitUVWeight.
 * returns a ObitUVWeight*.
 * in = object to reference
 */
#define ObitUVWeightRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitUVWeightIsA(in) ObitIsA (in, ObitUVWeightGetClass())

/*--------------- File Global Variables  ----------------*/
/**
 * ClassInfo structure ObitUVWeightClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitUVWeightClassInfo myClassInfo = {FALSE};

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitUVWeightInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitUVWeightClear (gpointer in);

/** Private: Get input parameters from uvdata */
void ObitUVWeightInput (ObitUVWeight *in, ObitUV *uvdata, ObitErr *err);

/** Private: Read uv data accumulating to grid */
void ObitUVWeightReadUV (ObitUVWeight *in, ObitUV *UVin, ObitErr *err);

/** Private: Reweight uv data accumulating to grids */
void ObitUVWeightWtUV (ObitUVWeight *in, ObitUV *UVin, ObitErr *err);

/** Private: Prepare visibility data for gridding weights*/
static void PrepBuffer (ObitUVWeight* in, ObitUV *uvdata);

/** Private: convolve a uv data buffer and sum to grid */
static void GridBuffer (ObitUVWeight* in, ObitUV *uvdata);

/** Private: Process grid */
static void ProcessGrid (ObitUVWeight* in, ObitErr *err);

/** Private: Apply corrections for a buffer of data */
static void WeightBuffer (ObitUVWeight* in, ObitUV *uvdata);

/** Private: Fill convolving function table */
static void ConvFunc (ObitUVWeight* in);

/** Private: Set Class function pointers. */
static void ObitUVWeightClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitUVWeight* newObitUVWeight (gchar* name)
{
  ObitUVWeight* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitUVWeightClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitUVWeight));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitUVWeightInit((gpointer)out);

 return out;
} /* end newObitUVWeight */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitUVWeightGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitUVWeightClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitUVWeightGetClass */

/**
 * Convolves UV weights onto a grid to determine weighting function.
 * Then the data weights are modified by the weighting function.
 * The control parameters are attached to the ObitInfoList member info
 * on uvdata.  See ObitUVWeight class documentation for details
 */
void ObitUVWeightData (ObitUV *uvdata, ObitErr *err)
{
  ObitUVWeight *myWeight = NULL;
  gchar *outName = NULL;
  olong naxis[2];
  gboolean doUnifWt;
  gchar *routine = "ObitUVWeightData";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitUVIsA(uvdata));

  /* create object */
  outName = g_strconcat("UVWeight for: ",uvdata->name, NULL);
  myWeight = newObitUVWeight(outName);
  g_free(outName);

  /* get weighting information */
  ObitUVWeightInput (myWeight, uvdata, err);
  if (err->error) Obit_traceback_msg (err, routine, uvdata->name);

  /* frequency tables if not defined */
  if ((uvdata->myDesc->freqArr==NULL) || (uvdata->myDesc->fscale==NULL)) {
    ObitUVGetFreq (uvdata, err);
    if (err->error) Obit_traceback_msg (err, routine, uvdata->name);
  } /* end setup frequency table */

  /* Are we uniform Weighting? */
  doUnifWt = myWeight->Robust < 7;

  /* Gridding for uniform weighting */
  if (doUnifWt) {
    /* Set convolving function */
    ConvFunc(myWeight);
    
    /* Create weighting grids */
    naxis[0] = 1 + myWeight->nuGrid/2;
    naxis[1] = myWeight->nvGrid;
    myWeight->cntGrid = ObitFArrayCreate ("Count Grid", 2, naxis);
    myWeight->wtGrid  = ObitFArrayCreate ("Weight Grid", 2, naxis);
    
    /* Get/grid weights if uniform weighting */
    ObitUVWeightReadUV (myWeight, uvdata, err);
    if (err->error) Obit_traceback_msg (err, routine, uvdata->name);
    
    /* Process grid */
    ProcessGrid (myWeight, err);
    if (err->error) Obit_traceback_msg (err, routine, uvdata->name);

    /* Informative messages */
     Obit_log_error(err, OBIT_InfoErr, 
		    "Using Robust uniform weighting for %s", uvdata->name);

     /* debug
     fprintf (stderr,"Grid size %d  %d \n",naxis[0],naxis[1]); */
  } else {
    /* Only natural weighting (and taper and power)*/
    myWeight->wtScale = 1.0;
    myWeight->temperance  = 0.0;

    /* Informative messages */
     Obit_log_error(err, OBIT_InfoErr, 
		    "Using natural weighting for %s", uvdata->name);
  } /* End of natural weighting */

  /* Modify Weights */ 
  ObitUVWeightWtUV (myWeight, uvdata, err);
  if (err->error) Obit_traceback_msg (err, routine, uvdata->name);

  /* final diagnostics */
     Obit_log_error(err, OBIT_InfoErr, 
		    "Weighting increased noise by %f for %s", 
		    myWeight->noiseFactor, uvdata->name);

  /* cleanup */
  myWeight = ObitUVWeightUnref(myWeight);
} /* end ObitUVWeightData */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitUVWeightClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitUVWeightClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitUVWeightClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitUVWeightClassInfoDefFn (gpointer inClass)
{
  ObitUVWeightClassInfo *theClass = (ObitUVWeightClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitUVWeightClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitUVWeightClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitUVWeightGetClass;
  theClass->newObit       = (newObitFP)newObitUVWeight;
  theClass->ObitCopy      = NULL;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitUVWeightClear;
  theClass->ObitInit      = (ObitInitFP)ObitUVWeightInit;
  theClass->ObitUVWeightData = (ObitUVWeightDataFP)ObitUVWeightData;

} /* end ObitUVWeightClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitUVWeightInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitUVWeight *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->cntGrid    = NULL;
  in->wtGrid     = NULL;
  in->convfn     = NULL;

} /* end ObitUVWeightInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * For some reason this wasn't build into the GType class.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitUVWeight* cast to an Obit*.
 */
void ObitUVWeightClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitUVWeight *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  in->cntGrid   = ObitFArrayUnref(in->cntGrid);  
  in->wtGrid    = ObitFArrayUnref(in->wtGrid);  
  in->convfn    = ObitFArrayUnref(in->convfn);
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitUVWeightClear */

/**
 * Read input parameters from uvdata.
 * The control parameters are attached to the ObitInfoList member info
 * on uvdata.  See ObitUVWeight Class documentation for details.
 * \param uvdata ObitUV whose data is to be weighted.
 * \param err    ObitErr stack for reporting problems.
 */
void ObitUVWeightInput (ObitUVWeight *in, ObitUV *uvdata, ObitErr *err)
{
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  ofloat temp, xCells, yCells, farr[10], *fptr, sigma2u, sigma2v, cpa, spa;
  olong   itemp, *iptr, nfield;
  gboolean gotIt;
  gchar *routine = "ObitUVWeightInput";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitUVIsA(uvdata));
  g_assert (ObitInfoListIsA(uvdata->info));

  /* Weighting grid size */
  itemp = 0;  
  gotIt = ObitInfoListGetP(uvdata->info, "nuGrid", &type, dim, (gpointer*)&iptr);
  if (gotIt) in->nuGrid = *iptr;
  else { /* better have "nx" */
    if (!ObitInfoListGetP(uvdata->info, "nx", &type, dim, (gpointer*)&iptr)) {
      Obit_log_error(err, OBIT_Error, 
		     "%s: MUST define grid size for %s", routine, uvdata->name);
      return;
    }
    in->nuGrid = *iptr;
  }
  /* Should be odd */
  in->nuGrid = ((in->nuGrid-1)/2) * 2 + 1;

  gotIt = ObitInfoListGetP(uvdata->info, "nvGrid", &type, dim, (gpointer*)&iptr);
  if (gotIt) in->nvGrid = *iptr;
  else { /* better have "ny" */
    if (!ObitInfoListGetP(uvdata->info, "ny", &type, dim, (gpointer*)&iptr)) {
      Obit_log_error(err, OBIT_Error, 
		     "%s: MUST define grid size for %s", routine, uvdata->name);
      return;
    }
    in->nvGrid = *iptr;
  }
  /* Should be odd */
  in->nvGrid = ((in->nvGrid-1)/2) * 2 + 1;
  nfield = dim[0];   /* Number of fields */
  
  /* Image cell size */
  if (!ObitInfoListGetP(uvdata->info, "xCells", &type, dim, (gpointer*)&fptr)) {
    Obit_log_error(err, OBIT_Error, 
		   "%s: MUST define cell spacing for %s", routine, uvdata->name);
    return;
  }
  xCells = (fptr[0])/3600.0;

  if (!ObitInfoListGetP(uvdata->info, "yCells", &type, dim, (gpointer*)&fptr)) {
    Obit_log_error(err, OBIT_Error, 
		   "%s: MUST define cell spacing for %s", routine, uvdata->name);
    return;
  }
  yCells = (fptr[0])/3600.0;

  /* WtBox */
  itemp = 0;
  ObitInfoListGetTest(uvdata->info, "WtBox", &type, dim, &itemp);
  in->WtBox = itemp;
 
  /* WtFunc */
  itemp = 1;
  ObitInfoListGetTest(uvdata->info, "WtFunc", &type, dim, &itemp);
  if (itemp<1) itemp = 1;
  in->WtFunc = itemp;
 
  /* Robust */
  temp = 0.0;
  ObitInfoListGetTest(uvdata->info, "Robust", &type, dim, &temp);
  in->Robust = temp;
 
  /* Weight power */
  temp = 1.0;
  ObitInfoListGetTest(uvdata->info, "WtPower", &type, dim, &temp);
  in->WtPower = fabs (temp);
 
  /* baseline range */
  temp = 1.0e15;
  ObitInfoListGetTest(uvdata->info, "MaxBaseline", &type, dim, &temp);
  in->blmax = temp;
  temp = 0.0;
  ObitInfoListGetTest(uvdata->info, "MinBaseline", &type, dim, &temp);
  in->blmin = temp;

  /* set uv to cells factors */
  in->UScale = in->nuGrid * (DG2RAD * fabs(xCells));
  in->VScale = in->nvGrid * (DG2RAD * fabs(yCells));

  /* Taper [default none] */
  farr[0] = farr[1] = farr[2] = 0.0;
  ObitInfoListGetTest(uvdata->info, "UVTaper", &type, dim, farr);
  if ((farr[0]>0.0) || (farr[1]>0.0)) {
    farr[0] *= 1.0e3;  /* To lambdas */
    farr[1] *= 1.0e3;  /* To lambdas */
    sigma2u = log(0.3)/(farr[0]*farr[0]);
    sigma2v = log(0.3)/(farr[1]*farr[1]);
    cpa = cos(farr[2]*DG2RAD);
    spa = sin(farr[2]*DG2RAD);
    in->sigma1 = (cpa*cpa*sigma2v + spa*spa*sigma2u);
    in->sigma2 = (spa*spa*sigma2v + cpa*cpa*sigma2u);
    in->sigma3 = 2.0*cpa*spa*(sigma2v - sigma2u);
  } else {
    in->sigma1 = 0.0;
    in->sigma2 = 0.0;
    in->sigma3 = 0.0;
  }

} /* end ObitUVWeightInput */

/**
 * Read a UV data object, applying any taper or shift and accumulating to grid.
 * Buffering of data will use the buffers as defined on UVin 
 * ("nVisPIO" in info member).
 * The UVin object will be closed at the termination of this routine.
 * Requires setup by #ObitUVWeightCreate.
 * The gridding information should have been stored in the ObitInfoList on in:
 * \li "UVTaper" OBIT_float [3] = UV taper (maj, min, pa) in kilowavelengths, deg.
 *             Default = no taper.
 * \li "Guardband" OBIT_float scalar = maximum fraction of U or v range allowed in grid.
 *             Default = 0.4.
 * \li "MaxBaseline" OBIT_float scalar = maximum baseline length in wavelengths.
 *             Default = 1.0e15.
 * \li "startChann" OBIT_long scalar = first channel (1-rel) in uv data to grid.
 *             Default = 1.
 * \li "numberChann" OBIT_long scalar = number of channels in uv data to grid.
 *             Default = all.
 * \param in      Object to initialize
 * \param UVin    Uv data object to be gridded.
 *                Should be the same as passed to previous call to 
 *                #ObitUVWeightSetup for input in.
 * \param err     ObitErr stack for reporting problems.
 */
void ObitUVWeightReadUV (ObitUVWeight *in, ObitUV *UVin, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_OK;
  gchar *routine = "ObitUVWeightReadUV";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitUVWeightIsA(in));
  g_assert (ObitUVIsA(UVin));
  g_assert (ObitUVDescIsA(UVin->myDesc));
  g_assert (UVin->myDesc->fscale!=NULL); /* frequency scaling table */

  retCode = ObitUVOpen (UVin, OBIT_IO_ReadOnly, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* loop gridding data */
  while (retCode == OBIT_IO_OK) {

    /* read buffer */
    retCode = ObitUVRead (UVin, NULL, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    
    /* prepare data */
    PrepBuffer (in, UVin);
    
    /* grid */
    GridBuffer (in, UVin);
  } /* end loop reading/gridding data */

  /* Close data */
  retCode = ObitUVClose (UVin, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
} /* end ObitUVWeightReadUV  */

/**
 * Apply any weighting to the uvdata in UVin and rewrite
 * Compute the increase in the noise due to the weighting.
 * \param in      Weighting object.
 * \param UVin    Uv data object to be corrected.
 * \param err     ObitErr stack for reporting problems.
 */
void ObitUVWeightWtUV (ObitUVWeight *in, ObitUV *UVin, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_OK;
  odouble sumInWt, sumOutWt, sumO2IWt, fract;
  olong firstVis;
  gchar *routine = "ObitUVWeightWtUV";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitUVWeightIsA(in));
  g_assert (ObitUVIsA(UVin));
  g_assert (ObitUVDescIsA(UVin->myDesc));
  g_assert (UVin->myDesc->fscale!=NULL); /* frequency scaling table */

  retCode = ObitUVOpen (UVin, OBIT_IO_ReadWrite, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Init noise factor info */
  in->noiseFactor = 1.0;
  in->wtSums[0] = 0.0; /* Sum of input weights      */
  in->wtSums[1] = 0.0; /* Sum of output weights     */
  in->wtSums[2] = 0.0; /* Sum Out_wt*Out_wt / In_wt */ 
  in->numberBad = 0;   /* Number of visibilities outside of the inner 90% */ 

  /* loop correcting data */
  while (retCode == OBIT_IO_OK) {

    /* read buffer */
    retCode = ObitUVRead (UVin, NULL, err);
    if (retCode == OBIT_IO_EOF) break; /* done? */
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    firstVis = UVin->myDesc->firstVis;
    
    /* Apply weighting */
    WeightBuffer (in, UVin);

    /* rewrite buffer */
    retCode = ObitUVWrite (UVin, NULL, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    UVin->myDesc->firstVis = firstVis;  /* reset first vis in buffer */
    ((ObitUVDesc*)UVin->myIO->myDesc)->firstVis = firstVis;
    
  } /* end loop weighting data */

  /* Close data */
  retCode = ObitUVClose (UVin, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Calculate noise increase factor */
  sumInWt  = in->wtSums[0];
  sumOutWt = in->wtSums[1];
  sumO2IWt = in->wtSums[2];

  /* Make sure all data not flagged */
  if (sumOutWt<=0.0) {
    Obit_log_error(err, OBIT_Error, 
		    "ERROR: All data flagged for %s", in->name);
    return;
   }

  if (sumOutWt < 1.0e-9) sumOutWt = 1.0;
  in->noiseFactor = sumO2IWt * sumInWt / (sumOutWt*sumOutWt);
  in->noiseFactor = sqrt (MAX (0.0, in->noiseFactor));

  /* Check for excessive numbers of visibilities flagged for being outside
     the inner 90% of the grid */
  fract = in->numberBad;
  fract = 100.0 * fract / UVin->myDesc->nvis;
  if (fract > 1.0) {
      Obit_log_error(err, OBIT_InfoWarn, 
		    "WARNING: %6.1f percent of data flagged outside of inner 0.9 for %s", 
		     fract, UVin->name);
  }

} /* end ObitUVWeightWtUV  */

/**
 * Prepares a buffer load of visibility data for gridding:
 * \li Convert to cells at the reference frequency.
 * \li enforce guardband - no data near outer edges of grid 
 * \li All data should be converted to the positive V half plane.
 * \param in      Object with grid to accumulate.
 * \param uvdata  Object with uvdata in buffer.
 */
static void PrepBuffer (ObitUVWeight* in, ObitUV *uvdata)
{
  olong ivis, nvis, ifreq, nif, iif, nfreq, loFreq, hiFreq;
  ofloat *u, *v, *w, *vis, *ifvis, *vvis;
  ofloat bl2, blmax2, blmin2, wt, guardu, guardv;
  ObitUVDesc *desc;
  gboolean flip, doFlag, doPower, doOne;

  /* error checks */
  g_assert (ObitUVWeightIsA(in));
  g_assert (ObitUVIsA(uvdata));
  g_assert (uvdata->myDesc != NULL);
  g_assert (uvdata->buffer != NULL);

  /* how much data? */
  desc  = uvdata->myDesc;
  nvis  = desc->numVisBuff;
  if (nvis<=0) return; /* need something */
  nfreq = desc->inaxes[desc->jlocf];
  nif = 1;
  if (desc->jlocif>=0) nif = desc->inaxes[desc->jlocif];
  
  /* range of channels (0-rel) */
  loFreq = 0;
  hiFreq = nfreq-1;

 /* initialize data pointers */
  u   = uvdata->buffer+desc->ilocu;
  v   = uvdata->buffer+desc->ilocv;
  w   = uvdata->buffer+desc->ilocw;
  vis = uvdata->buffer+desc->nrparm;

  /* what needed */
  /* Raising weight to a power? */
  doPower = (fabs (in->WtPower-1.0) > 0.01) && (in->WtPower > 0.01);
  /* Replacing weights with 1.0? */
  doOne = (in->WtPower < 0.01);

 /* Baseline max, min values */
  blmax2 = in->blmax * in->blmax;
  blmin2 = in->blmin * in->blmin;

  /* guardband (90%)in wavelengths */
  guardu = (0.9 * (ofloat)in->wtGrid->naxis[0]) / fabs(in->UScale);
  guardv = (0.9 * ((ofloat)in->wtGrid->naxis[1])/2) / fabs(in->VScale);

  /* Loop over visibilities */
  for (ivis=0; ivis<nvis; ivis++) {

    /* check extrema */
    bl2 = (*u)*(*u) + (*v)*(*v);
    doFlag = ((bl2<blmin2) || (bl2>blmax2));
    /* enforce guardband */
    doFlag = doFlag || ((fabs(*u)>guardu) || (fabs(*v)>guardv));

    /* in the correct half plane? */
    flip = (*u) < 0.0;

    /* loop over IFs */
    ifvis = vis;
    for (iif = 0; iif<nif; iif++) {

      /* loop over frequencies */
      vvis = ifvis;
      for (ifreq = loFreq; ifreq<=hiFreq; ifreq++) {

	/* is this one wanted? */
	if (doFlag)  vvis[2] = 0.0;  /* baseline out of range? */
	
	wt = vvis[2];                /* data weight */
	if (wt <= 0.0) continue;
	
	/* Replacing weights with one? */
	if (doOne) vis[2] = 1.0;

	/* Weights to a power? */
	if (doPower && (vis[2]>0.0)) vis[2] = pow (vis[2], in->WtPower);
	
	vvis += desc->incf; /* visibility pointer */
      } /* end loop over frequencies */
      ifvis += desc->incif; /* visibility pointer */
    } /* Loop over IFs */

    /* Scale u,v to cells at reference frequency */
    if (flip) { /* put in other half plane */
      *u = -((*u) * in->UScale);
      *v = -((*v) * in->VScale);
    } else { /* no flip */
      *u = (*u) * in->UScale;
      *v = (*v) * in->VScale;
    }

    /* update data pointers */
    u += desc->lrec;
    v += desc->lrec;
    w += desc->lrec;
    vis += desc->lrec;
  } /* end loop over visibilities */
} /* end PrepBuffer */

/**
 * Convolves weights in buffer on uvdata onto the grid member of in.
 * Rows in the grid are in U and the data should have all been converted to the 
 * positive U half plane.
 * U, V, should be in cells and data not to be included on the grid should 
 * have zero weight.  Convolution functions must be created.
 * This uses two grids, one for the counts of visibilities in cells and the other
 * for the sum of the weights.  These are needed for Briggs Robust weighting.
 * \param in      Object with grid to accumulate
 * \param uvdata  Object with uv data in buffer, prepared for gridding.
 */
static void GridBuffer (ObitUVWeight* in, ObitUV *uvdata)
{
  olong ivis, nvis, ifreq, nfreq, ncol, iu, iv, icu, icv, lGridRow, lGridCol, itemp;
  olong istok, nstok;
  olong iif, ifq, nif, loFreq, hiFreq, uoff, voff, uuoff=0.0, vvoff, vConvInc, uConvInc;
  ofloat *grid, *ggrid, *cntGrid, *u, *v, *w, *vis, *vvis, *fvis, *ifvis, *wt;
  ofloat *convfnp, weight, rtemp, uf, vf;
  olong fincf, fincif;
  olong pos[] = {0,0,0,0,0};
  ObitUVDesc *desc;

  /* error checks */
  g_assert (ObitUVWeightIsA(in));
  g_assert (ObitUVIsA(uvdata));
  g_assert (in->cntGrid != NULL);
  g_assert (in->wtGrid != NULL);
  g_assert (uvdata->myDesc != NULL);
  g_assert (uvdata->buffer != NULL);

  /* how much data? */
  desc  = uvdata->myDesc;
  nvis  = desc->numVisBuff;
  if (nvis<=0) return; /* need something */
  nfreq = desc->inaxes[desc->jlocf];
  nif = 1;
  if (desc->jlocif>=0) nif = desc->inaxes[desc->jlocif]; 
  nstok = 1;
  if (desc->jlocs>=0) nstok = desc->inaxes[desc->jlocs];
  nstok = MIN (2, nstok);

  /* range of channels (0-rel) */
  loFreq = 0;
  hiFreq = nfreq-1;

  /* Channel and IF increments in frequcncy scaling array */
  fincf  = MAX (1, (desc->incf  / 3) / desc->inaxes[desc->jlocs]);
  fincif = MAX (1, (desc->incif / 3) / desc->inaxes[desc->jlocs]);

 /* initialize data pointers */
  u   = uvdata->buffer+desc->ilocu;
  v   = uvdata->buffer+desc->ilocv;
  w   = uvdata->buffer+desc->ilocw;
  vis = uvdata->buffer+desc->nrparm;

  lGridRow = in->cntGrid->naxis[0]; /* length of row */
  lGridCol = in->cntGrid->naxis[1]; /* length of column */

  /* convolution fn pointer */
  pos[0] = 0;
  convfnp = ObitFArrayIndex (in->convfn, pos);

  /* Loop over visibilities */
  for (ivis=0; ivis<nvis; ivis++) {

    /* loop over IFs */
    ifvis = vis;
    for (iif=0; iif<nif; iif++) {

      /* loop over frequencies */
      fvis = ifvis;
      for (ifreq = loFreq; ifreq<=hiFreq; ifreq++) {
	ifq = iif*fincif + ifreq*fincf;  /* index in IF/freq table */

	/* Scale u,v for frequency (w not used) - already in cells and correct half plane */
	uf = (*u) * desc->fscale[ifq];
	vf = (*v) * desc->fscale[ifq];
	
	/* get center cell */
	if (vf > 0.0) iv = (olong)(vf + 0.5);
	else iv = (olong)(vf - 0.5);
	iu = (olong)(uf + 0.5);

	/* Add this visibility to the count grid */
	pos[0] = iu;
	pos[1] = iv + lGridCol/2;
	cntGrid = ObitFArrayIndex (in->cntGrid, pos); /* pointer in grid */

	/* Check if datum in grid - cntGrid != NULL */
	if  (cntGrid != NULL) {
	  *cntGrid += 1.0;  /* increment count */

	  /* Loop over stokes */
	  vvis = fvis;
	  for (istok=0; istok<nstok; istok++) {
	
	    /* is this one wanted? */
	    wt = vvis + 2; /* data weight */
	    if (*wt <= 0.0) {vvis += desc->incs; continue;}
	    
	    /* weight to grid */
	    weight = (*wt);
	    
	    /* convolve weight onto the weight grid */
	    /* back off half Kernel width */
	    iu -= in->convWidth/2;
	    iv -= in->convWidth/2;
	    
	    /* Starting convolution location, table has in->convNperCell points per cell */
	    /* Determine fraction of the cell to get start location in convolving table. */
	    if (uf > 0.0) itemp = (olong)(uf + 0.5);
	    else itemp = ((olong)(uf - 0.5));
	    rtemp = in->convNperCell*(itemp - (uf) - 0.5);
	    if (rtemp > 0.0) rtemp += 0.5;
	    else rtemp -= 0.5;
	    uoff = (olong)rtemp + in->convNperCell;
	    /* Increment between u convolving */
	    uConvInc = in->convNperCell;
	    
	    /* now v convolving fn */
	    if (vf > 0.0) itemp = (olong)(vf + 0.5);
	    else itemp = ((olong)(vf - 0.5));
	    rtemp = in->convNperCell*(itemp - (vf) - 0.5);
	    if (rtemp > 0.0) rtemp += 0.5;
	    else rtemp -= 0.5;
	    voff = (olong)rtemp + in->convNperCell;
	    /* Increment between v convolving entries */
	    vConvInc = in->convWidth * in->convNperCell + 1;
	    voff = voff * vConvInc;
	    
	    /* if too close to the center, have to break up and do conjugate halves */
	    if (iu > 0) { /* all in same half */
	      ncol = in->convWidth;
	      pos[0] = iu;
	      pos[1] = iv + lGridCol/2;
	      
	    } else { 
	      /* have to split - grid part in conjugate half */
	      pos[0] = -iu;
	      pos[1] = -iv + lGridCol/2;
	      grid = ObitFArrayIndex (in->wtGrid, pos); /* pointer in grid */
	      ncol = -iu;
	      vvoff = voff;
	      for (icv=0; icv<in->convWidth; icv++) {
		uuoff = uoff;
		ggrid  = grid;
		for (icu=0; icu<=ncol; icu++) {
		  *ggrid   += weight * convfnp[uuoff+vvoff];
		  uuoff += uConvInc;  /* Convolution kernel pointer */
		  ggrid -= 1; /* gridding pointer - opposite of normal gridding */
		} /* end inner gridding loop */
		vvoff += vConvInc;  /* Convolution kernel pointer */
		grid -= lGridRow; /* gridding pointer - reverse direction for conjugate */
	      } /* end outer loop */
	      
	      /* set up for rest of grid */
	      ncol = (in->convWidth + iu); /* how many columns left? */
	      iu = 0;      /* by definition  start other half plane at iu=0 */
	      pos[0] = iu; 
	      pos[1] = iv + lGridCol/2;
	      uoff = uuoff - in->convNperCell; /* for other half in u */
	    } /* End of dealing with conjugate portion */
	    
	    /* main loop gridding */
	    grid = ObitFArrayIndex (in->wtGrid, pos); /* pointer in grid */
	    vvoff = voff;
	    for (icv=0; icv<in->convWidth; icv++) {
	      uuoff = uoff;
	      ggrid  = grid;
	      for (icu=0; icu<ncol; icu++) {
		*ggrid   += weight * convfnp[uuoff+vvoff];
		uuoff += uConvInc;  /* Convolution kernel pointer */
		ggrid += 1; /* gridding pointer */
	      } /* end inner gridding loop */
	      vvoff += vConvInc;  /* Convolution kernel pointer */
	      grid += lGridRow; /* gridding pointer */
	    } /* end outer gridding loop */
	    vvis += desc->incs; /* visibility pointer */
	  } /* end of Stokes gridding loop */
	  
	  fvis += desc->incf; /* visibility pointer */
	} /* End of if datum in grid */
      } /* end loop over frequencies */
      ifvis += desc->incif; /* visibility pointer */
    } /* Loop over IFs */

    /* update data pointers */
    u += desc->lrec;
    v += desc->lrec;
    w += desc->lrec;
    vis += desc->lrec;
  } /* end loop over visibilities */
} /* end GridBuffer */

/**
 * Processes grid for Robust weighting.
 * Calculates:
 * \li temperance = factor for Briggs Robust weighting
 * \li wtScale = weighting normalization factor.
 * \param in   Object with grid.
 * \param err  Object for informative messages and errors.
 */
static void ProcessGrid (ObitUVWeight* in, ObitErr *err)
{
  ofloat sumWtCnt, sumCnt;
  /* error checks */
  g_assert (ObitUVWeightIsA(in));
  g_assert (in->cntGrid != NULL);
  g_assert (in->wtGrid != NULL);

  /* Get Sum of Weights*counts */
  sumWtCnt = ObitFArrayDot(in->cntGrid, in->wtGrid);

  /* Get sum of counts */
  sumCnt = ObitFArraySum(in->cntGrid);

  /* Robust temperance value */
  in->wtScale    = sumWtCnt / MAX (1.0, sumCnt);
  in->temperance = in->wtScale * pow (10.0, in->Robust) / 5.0;

  /* If Robust out of range turn it off */
  if ((in->Robust<-7.0) || (in->Robust>7.0)) in->temperance = 0.0;

  /* Weight scaling factor */
  in->wtScale += in->temperance;

} /* end ProcessGrid */

/**
 * Corrects data in buffer using weighting grid.
 * Adds temperance factor in->temperance and multiplies by in->wtScale
 * \param in      Object with grid to accumulate
 * \param uvdata  Object with uv data in buffer, prepared for gridding.
 */
static void WeightBuffer (ObitUVWeight* in, ObitUV *uvdata)
{
  olong ivis, nvis, ifreq, nfreq, iu, iv, lGridRow, lGridCol=0;
  olong istok, nstok;
  olong ifq, iif, nif, loFreq, hiFreq;
  ofloat *grid=NULL, *u, *v, *w, *vis, *vvis, *fvis, *ifvis, *wt;
  ofloat tape, tfact, inWt, outWt, guardu, guardv, uf, vf, minWt;
  ofloat ucell, vcell, uucell, vvcell;
  olong pos[] = {0,0,0,0,0};
  olong fincf, fincif;
  ObitUVDesc *desc;
  gboolean doPower, doOne, doTaper, doUnifWt, doFlag;
  odouble sumInWt, sumOutWt, sumO2IWt,numberBad ;

  /* error checks */
  g_assert (ObitUVWeightIsA(in));
  g_assert (ObitUVIsA(uvdata));
  g_assert (uvdata->myDesc != NULL);
  g_assert (uvdata->buffer != NULL);

  /* initialize weighting sums */
  sumInWt  = in->wtSums[0];
  sumOutWt = in->wtSums[1];
  sumO2IWt = in->wtSums[2];
  numberBad = in->numberBad;

  /* how much data? */
  desc  = uvdata->myDesc;
  nvis  = desc->numVisBuff;
  if (nvis<=0) return; /* need something */
  nfreq = desc->inaxes[desc->jlocf];
  nif = 1;
  if (desc->jlocif>=0) nif = desc->inaxes[desc->jlocif];
  nstok = 1;
  if (desc->jlocs>=0) nstok = desc->inaxes[desc->jlocs];
  nstok = MIN (2, nstok);

  /* Minimum allowed weight */
  minWt = 1.0e-15;
  if (in->temperance>0.0) minWt = 1.0e-7 * in->temperance;
  if (minWt==0.0)  minWt = 1.0e-25;
 
  /* Channel and IF increments in frequency scaling array */
  fincf  = MAX (1, (desc->incf  / 3) / desc->inaxes[desc->jlocs]);
  fincif = MAX (1, (desc->incif / 3) / desc->inaxes[desc->jlocs]);

  /* range of channels (0-rel) */
  loFreq = 0;
  hiFreq = nfreq-1;

 /* initialize data pointers */
  u   = uvdata->buffer+desc->ilocu;
  v   = uvdata->buffer+desc->ilocv;
  w   = uvdata->buffer+desc->ilocw;
  vis = uvdata->buffer+desc->nrparm;

  /* what needed */
  /* Need taper? */
  doTaper = (in->sigma1!=0.0) || (in->sigma2!=0.0);
  /* Raising weight to a power? */
  doPower = (fabs (in->WtPower-1.0) > 0.01) && (in->WtPower > 0.01);
  /* Replacing weights with 1.0? */
  doOne = (in->WtPower < 0.01);
  /* Uniform Weighting */
  doUnifWt = in->Robust < 7;

  if (doUnifWt) {
    g_assert (in->wtGrid != NULL);

    lGridRow = in->wtGrid->naxis[0]; /* length of row */
    lGridCol = in->wtGrid->naxis[1]; /* length of column */
    /* guardband (90%)in wavelengths */
    guardu = (0.9 * (ofloat)in->wtGrid->naxis[0]) / fabs(in->UScale);
    guardv = (0.9 * ((ofloat)in->wtGrid->naxis[1])/2) / fabs(in->VScale);
  } else { /* Not uniform weighting */
    guardu = 1.0e30;
    guardv = 1.0e30;
  }

  /* Loop over visibilities */
  for (ivis=0; ivis<nvis; ivis++) {

    /* enforce guardband */
    doFlag = FALSE;
    if ((fabs(*u)>guardu) || (fabs(*v)>guardv)) {
      doFlag = TRUE;
      numberBad++;
    }

    /* Scale u,v to cells at reference frequency */
    if ((*u) <= 0.0) { /* put in other half plane */
      ucell = -((*u) * in->UScale);
      vcell = -((*v) * in->VScale);
    } else { /* no flip */
      ucell = (*u) * in->UScale;
      vcell = (*v) * in->VScale;
    }
 
    /* loop over IFs */
    ifvis = vis;
    for (iif=0; iif<nif; iif++) {

      /* loop over frequencies */
      fvis = ifvis;
      for (ifreq = loFreq; ifreq<=hiFreq; ifreq++) {
	ifq = iif*fincif + ifreq*fincf;  /* index in IF/freq table */

	  /* Loop over stokes */
	  vvis = fvis;
	  for (istok=0; istok<nstok; istok++) {
	
	    /* is this one wanted? */
	    wt = vvis + 2; /* data weight */
	    if (doFlag) *wt = 0.0;
	    if (*wt <= 0.0) {vvis += desc->incs; continue;}
	    
	    /* Input weight */
	    inWt = *wt;
	    
	    /* Replacing weights with one? */
	    if (doOne) vis[2] = 1.0;
	    
	    /* Weights to a power? */
	    if (doPower && (vis[2]>0.0)) vis[2] = pow (vis[2], in->WtPower);
	    
	    /* Doing uniform weighting? */
	    if (doUnifWt) {
	      /* Scale u,v (cells) for frequency (w not used) */
	      uucell = ucell * desc->fscale[ifq];
	      vvcell = vcell * desc->fscale[ifq];
	      
	      /* get center cell */
	      if (vvcell > 0.0) iv = (olong)(vvcell + 0.5);
	      else iv = (olong)(vvcell - 0.5);
	      iu = (olong)(uucell + 0.5);
	      
	      /* location in grid */
	      pos[0] = iu; 
	      pos[1] = iv + lGridCol/2;
	      grid = ObitFArrayIndex (in->wtGrid, pos); /* pointer in weight grid */
	      
	      /* zero weight if not in grid */
	      if ((grid==NULL) || (fabs(*grid)==0.0)) {
		*wt = 0.0;
	      } else {
		/* Make weighting correction */
		if (*grid >minWt) *wt *= in->wtScale / (*grid + in->temperance);
		/* OLD else    	  *wt *= in->wtScale / (minWt + in->temperance);*/
		else    	  *wt *= in->wtScale;
	      }
	    } /* end of uniform weighting correction */
	    
	    /* apply any taper to the weight. */
	    if (doTaper) {
	      /* Scale u,v (wavelengths) for frequency (w not used) */
	      uf = *u * desc->fscale[ifq];
	      vf = *v * desc->fscale[ifq];
	      
	      tape = ((uf)*(uf)*in->sigma2 + (vf)*(vf)*in->sigma1 + (uf)*(vf)*in->sigma3);
	      if (tape<-14.0) tfact = 0.0; /* underflow */
	      else tfact = exp(tape);
	      *wt *= tfact;
	    }
	    
	    /* Output weight */
	    outWt = *wt;
	    
	    /* Weighting sums for statistics */
	    sumInWt  += inWt;
	    sumOutWt += outWt;
	    sumO2IWt += outWt * outWt / inWt;
	    vvis += desc->incs; /* visibility pointer */
	  } /* end of Stokes gridding loop */
	    
	  fvis += desc->incf; /* visibility pointer */
      } /* end loop over frequencies */
      ifvis += desc->incif; /* visibility pointer */
    } /* Loop over IFs */

    /* update data pointers */
    u += desc->lrec;
    v += desc->lrec;
    w += desc->lrec;
    vis += desc->lrec;
  } /* end loop over visibilities */

  /* save weighting sums */
  in->wtSums[0] = sumInWt;
  in->wtSums[1] = sumOutWt;
  in->wtSums[2] = sumO2IWt;
  in->numberBad = numberBad;

} /* end WeightBuffer */

/**
 * Calculates convolving function and attaches it to in.
 * Makes 2D convolving function tabulated every 1/5 of a cell.
 * \param in      Object with table to init.
 */
static void ConvFunc (ObitUVWeight* in)
{
  ofloat xinc, fwt, fu, fv=0.0, rad, radmax, xi=0.0;
  ofloat u, v, absu, absv=0.0, umax, vmax, *convfnp;
  olong i, j, k, size, lim, bias, naxis[2];
  gboolean round;
  /*gfloat shit[701]; DEBUG */

  /* error checks */
  g_assert (ObitUVWeightIsA(in));

  /* Width of convolving kernel in cells */
  in->convWidth    = 2*in->WtBox+1;  /* number of cells - should be odd */
  in->convNperCell = 5;   /* Number of tabulated points per cell in convfn */
  round = (in->WtFunc > 0) && (in->WtBox>0); /* Round or square function? */

  /* allocate array */
  lim = in->convWidth * in->convNperCell + 1;
  size = lim;
  naxis[0] = size;
  naxis[1] = size;
  in->convfn = ObitFArrayUnref(in->convfn);
  in->convfn = ObitFArrayCreate (in->name, 2L, naxis);

  /* get pointer to memory array */
  naxis[0] = 0; naxis[1] = 0;
  convfnp = ObitFArrayIndex (in->convfn, naxis);

  xinc = 1.0 / (ofloat)in->convNperCell;
  if (in->WtBox==0) {  /* if using a single cell, use everything */
    in->WtFunc = 1;    /* only use pill box for single cell */
    umax = 1.0;
    vmax = 1.0;
  } else { /* multi cell weighting function */
    umax = in->WtBox + 0.5;
    vmax = in->WtBox + 0.5;
  }
  radmax = umax * vmax;
  bias = (in->convNperCell/2) * in->convWidth;

  /*+++++++++++++++++ Pillbox ++++++++++++++++++++++++++++++++++++++++*/
  if (abs(in->WtFunc)==1) {
    /* fill function */
    k = 0;
    for (j=0; j<lim; j++) {
      /* distance in v */
      v = (j-bias) * xinc;
      absv = fabs (v);
      fv =  1.0;
      if (absv == vmax) fv = 0.5;
      else if (absv > vmax) fv = 0.0;
      for (i=0; i<lim; i++) {
	/* distance in u */
	u = (i-bias) * xinc;
	if (round) {
	  /* Use circularly symmetric version */
	  rad = u*u + v*v;
	  fwt = 1.0;
	  if (rad == radmax) fwt = 0.5;
	  else if (rad > radmax)  fwt = 0.0;
	} else {
	  /* use square version */
	  absu = fabs (u);
	  fu = 1.0;
	  if (absu == umax) fu = 0.5;
	  else if (absu > umax)  fu = 0.0;
	  fwt = fu * fv;
	}
	convfnp[k++] = fwt;  /* save it in u major order */
      }
    } 
  } /* end pillbox */
  
  else if (abs(in->WtFunc)==2) {
    /*+++++++++++++++++ Linear ++++++++++++++++++++++++++++++++++++++++*/
    /* fill function */
    xi = 1.0 / ((ofloat)in->WtBox + 1.0);
    k = 0;
    for (j=0; j<lim; j++) {
      /* distance in v */
      v = (j-bias) * xinc;
      absv = fabs (v);
      fv =  1.0 - absv*xi;
      if (absv > vmax) fv = 0.0;
      for (i=0; i<lim; i++) {
	/* distance in u */
	u = (i-bias) * xinc;
	if (round) {
	  /* Use circularly symmetric version */
	  rad = sqrt(u*u + v*v);
	  fwt = 1.0 - rad*xi;
	  if (rad > radmax)  fwt = 0.0;
	} else {
	  /* use square version */
	  absu = fabs (u);
	  fu = 1.0- absu*xi;
	  if (absu > umax)  fu = 0.0;
	  fwt = fu * fv;
	}
	convfnp[k++] = fwt;  /* save it in u major order */
      }
    }  
  } /* end linear */
  
  else if (abs(in->WtFunc)==3) {
    /*+++++++++++++++++ Exponential ++++++++++++++++++++++++++++++++++++++++*/
    /* fill function */
    k = 0;
    for (j=0; j<lim; j++) {
      /* distance in v */
      v = (j-bias) * xinc;
      if (!round) {
	absv = fabs (v);
	fv = exp(-2.0*absv*xi);
	if (absv > vmax) fv = 0.0;
      }
      for (i=0; i<lim; i++) {
	/* distance in u */
	u = (i-bias) * xinc;
	if (round) {
	  /* Use circularly symmetric version */
	  rad = sqrt(u*u + v*v);
	  fwt = exp(-2.0*rad*xi);
	  if (rad > radmax)  fwt = 0.0;
	} else {
	  /* use square version */
	  absu = fabs (u);
	  fu = exp(-2.0*absu*xi);
	  if (absu > umax) fu = 0.0;
	  fwt = fu * fv;
	}
	convfnp[k++] = fwt;  /* save it in u major order */
      }
    }
  } /* end exponential */
  
  else if (abs(in->WtFunc)==4) {
    /*+++++++++++++++++ Gaussian ++++++++++++++++++++++++++++++++++++++++*/
    /* fill function */
    k = 0;
    for (j=0; j<lim; j++) {
      /* distance in v */
      v = (j-bias) * xinc;
      if (!round) {
	absv = fabs (v);
	fv = exp(-4.0*absv*xi*absv*xi);
	if (absv > vmax) fv = 0.0;
      }
      else if (absv >  vmax) fv = 0.0;
      for (i=0; i<lim; i++) {
	/* distance in u */
	u = (i-bias) * xinc;
	if (round) {
	  /* Use circularly symmetric version */
	  rad = sqrt(u*u + v*v);
	  fwt = exp(-4.0*rad*xi*rad*xi);
	  if (rad == radmax) fwt = 0.5;
	  else if (rad > radmax)  fwt = 0.0;
	} else {
	  /* use square version */
	  absu = fabs (u);
	  fu = exp(-4.0*absu*xi*absu*xi);
	  if (absu > umax)  fu = 0.0;
	  fwt = fu * fv;
	}
	convfnp[k++] = fwt;  /* save it in u major order */
      }
    }
  } /* end computing convolving fns */
  else { /* should never get here */
    g_error("Unknown convolving function type %d",in->WtFunc);
  }
} /* end ConvFunc */



