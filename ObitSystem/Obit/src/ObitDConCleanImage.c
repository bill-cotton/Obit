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

#include "ObitDConCleanImage.h"
#include "ObitMem.h"
#include "ObitFFT.h"
#include "ObitTableCCUtil.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitDConCleanImage.c
 * ObitDConCleanImage class function definitions.
 * Image based CLEAN class.
 * This class is derived from the ObitDCon class.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitDConCleanImage";

/** Function to obtain parent ClassInfo - ObitDConClean */
static ObitGetClassFP ObitParentGetClass = ObitDConCleanGetClass;

/**
 * ClassInfo structure ObitDConCleanImageClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitDConCleanImageClassInfo myClassInfo = {FALSE};

/*--------------- File Global Variables  ----------------*/


/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitDConCleanImageInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitDConCleanImageClear (gpointer in);

/** Private: Compute fransfer function. */
static void  GetXfer (ObitDConCleanImage *in, ObitErr *err);

/** Private: Set Class function pointers. */
static void ObitDConCleanImageClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitDConCleanImage* newObitDConCleanImage (gchar* name)
{
  ObitDConCleanImage* out;
  /*gchar *routine = "newObitDConCleanImage";*/

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitDConCleanImageClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitDConCleanImage));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitDConCleanImageInit((gpointer)out);

 return out;
} /* end newObitDConCleanImage */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitDConCleanImageGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitDConCleanImageClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitDConCleanImageGetClass */

/**
 * Make a deep copy of an ObitDConCleanImage.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitDConCleanImage* ObitDConCleanImageCopy  (ObitDConCleanImage *in, 
					     ObitDConCleanImage *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  gboolean oldExist;
  gchar *outName;
  gchar *routine = "ObitDConCleanImageCopy";

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
    out = newObitDConCleanImage(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, out);

  /*  copy this class */
  out->transfer = ObitCArrayUnref(out->transfer);
  out->transfer = ObitCArrayCopy (in->transfer, out->transfer, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, out);
  out->forFFT = ObitFFTUnref(out->forFFT);
  out->forFFT = ObitFFTRef(out->forFFT);
  out->revFFT = ObitFFTUnref(out->revFFT);
  out->revFFT = ObitFFTRef(out->revFFT);
  out->startComp = in->startComp;
  out->endComp   = in->endComp;

  return out;
} /* end ObitDConCleanImageCopy */

/**
 * Make a copy of a object but do not copy the actual data
 * This is useful to create an DConCleanImage similar to the input one.
 * \param in  The object to copy
 * \param out An existing object pointer for output, must be defined.
 * \param err Obit error stack object.
 */
void ObitDConCleanImageClone  (ObitDConCleanImage *in, ObitDConCleanImage *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  gchar *routine = "ObitDConCleanImageClone";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitIsA(out, &myClassInfo));

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /*  copy this class */
  out->transfer = ObitCArrayUnref(out->transfer);
  out->transfer = ObitCArrayRef(in->transfer);
  out->forFFT = ObitFFTUnref(out->forFFT);
  out->forFFT = ObitFFTRef(out->forFFT);
  out->revFFT = ObitFFTUnref(out->revFFT);
  out->revFFT = ObitFFTRef(out->revFFT);
  out->startComp = in->startComp;
  out->endComp   = in->endComp;

} /* end ObitDConCleanImageClone */

/**
 * Creates an ObitDConCleanImage 
 * defined for convenience of derived classes 
 * \param name  An optional name for the object.
 * \param mosaic from which to create object
 * \param err Obit error stack object.
 * \return the new object.
 */
ObitDConCleanImage* ObitDConCleanImageCreate (gchar* name, ObitImageMosaic *mosaic,  
			  ObitErr *err)
{
  olong nfield;
  ObitDConCleanImage* out=NULL;
  gchar *routine = "ObitDConCleanImageCreate";

 /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return out;
  g_assert (ObitImageMosaicIsA(mosaic));

  /* Create basic structure */
  out = newObitDConCleanImage (name);

  /* Save Image Mosaic reference */
  out->mosaic = ObitImageMosaicRef(mosaic);

  /* Window object */
  out->window = ObitDConCleanWindowCreate ("CleanWindow", mosaic, err);
  if (err->error) Obit_traceback_val (err, routine, name, out);

  /* Arrays per field - including those in parent classes */
  nfield =  mosaic->numberImages;
  out->gain    = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean Loop gain");
  out->minFlux = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean minFlux");
  out->factor  = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean factor");
  out->nfield  = nfield;

  return out;
} /* end ObitDConCleanImageCreate */

/**
 * Do deconvolution, uses function on class pointer
 * Does final flatten if FullField member of mosaic member is defined.
 * CLEAN control parameters are in the ObitInfoList member:
 * \li "Niter"   OBIT_long scalar   = Maximum number of CLEAN iterations
 * \li "maxPixel" OBIT_long scalar  = Maximum number of residuals [def 20000]
 * \li "minPatch" OBIT_long scalar  = Minimum beam patch in pixels [def 50]
 * \li "BMAJ"    OBIT_float scalar = Restoring beam major axis (deg)
 * \li "BMIN"    OBIT_float scalar = Restoring beam minor axis (deg)
 * \li "BPA"     OBIT_float scalar = Restoring beam position angle (deg)
 * \li "Beam"    OBIT_float array[3]= (BMAJ, BMIN, BPA) alternate form
 * \li "CCVer"   OBIT_long array  = CLEAN table version per field
 * \li "Gain"    OBIT_float array  = CLEAN loop gain per field
 * \li "minFlux" OBIT_float array  = Minimun flux density (Jy)  per field
 * \li "Factor"  OBIT_float array  = CLEAN depth factor per field
 * \li "Plane"   OBIT_long array    = Plane being processed, 1-rel indices of axes 3-?
 * \param in   The object to deconvolve
 * \param err Obit error stack object.
 */
void ObitDConCleanImageDeconvolve (ObitDCon *inn, ObitErr *err)
{
  ObitDConCleanImage *in;
  gboolean done;
  const ObitDConCleanImageClassInfo *inClass;
  gchar *routine = "ObitDConCleanImageDeconvolve";

  /* Cast input to this type */
  in = (ObitDConCleanImage*)inn;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitDConCleanImageIsA(in));

  inClass = (ObitDConCleanImageClassInfo*)in->ClassInfo; /* class structure */

  /* Get parameters */
  inClass->ObitDConGetParms(inn, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Set transfer function */
  GetXfer (in, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Create Pixel List if needed */
  if (!in->Pixels) {
    in->Pixels = ObitDConCleanPxListCreate("Pixel List", in->mosaic, 
					   in->maxPixel, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
  }

  /* Copy control info to PixelList */
  ObitInfoListCopyData(in->info, in->Pixels->info);

  /* Reset/Init Pixel list*/
  ObitDConCleanPxListReset (in->Pixels, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Loop until Deconvolution done */
  done = FALSE;
  while (!done) {
    in->currentField = 1;  /* Can only do oone */
    /* Get image/beam statistics needed for this cycle */
    inClass->ObitDConCleanPixelStats((ObitDConClean*)in, NULL, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);

    /* Pick components for this major cycle, tells if finished CLEAN */
    done = inClass->ObitDConCleanSelect((ObitDConClean*)in, NULL, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);

    /* Subtract components and make new residual(s) */
    inClass->ObitDConCleanSub((ObitDConClean*)in, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);

  } /* end clean loop */

  /* Restore */
  if (in->doRestore)
    inClass->ObitDConCleanRestore((ObitDConClean*)in, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

} /* end ObitDConCleanImageDeconvolve */

/**
 * Read any base class parameters and then
 * read CLEAN control parameters from the ObitInfoList member:
 * \li "doRestore" OBIT_bool       = Restore image when done? [def TRUE]
 *
 * From Parent classes:
 * \li "Niter"   OBIT_long scalar   = Maximum number of CLEAN iterations
 * \li "maxPixel" OBIT_long scalar  = Maximum number of residuals [def 20000]
 * \li "minPatch" OBIT_long scalar  = Minimum beam patch in pixels [def 100]
 * \li "BMAJ"    OBIT_float scalar = Restoring beam major axis (deg)
 * \li "BMIN"    OBIT_float scalar = Restoring beam minor axis (deg)
 * \li "BPA"     OBIT_float scalar = Restoring beam position angle (deg)
 * \li "Beam"    OBIT_float array[3]= (BMAJ, BMIN, BPA) alternate form
 * \li "CCVer"   OBIT_long array  = CLEAN table version per field
 * \li "Gain"    OBIT_float array  = CLEAN loop gain per field
 * \li "minFlux" OBIT_float array  = Minimun flux density (Jy)  per field
 * \li "Factor"  OBIT_float array  = CLEAN depth factor per field
 * \li "Plane"   OBIT_long array    = Plane being processed, 1-rel indices of axes 3-?
 * \param in  The CLEAN object as base class
 * \param err Obit error stack object.
 */
void  ObitDConCleanImageGetParms (ObitDCon *inn, ObitErr *err)
{
  ObitDConCleanImage *in = (ObitDConCleanImage*)inn;  /* as this class */
  ObitDConClassInfo *ParentClass;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  /*    olong itemp;
	union ObitInfoListEquiv InfoReal; */
  gchar *routine = "ObitDConCleanImageGetParms";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Read any parent class parameters */
  ParentClass = (ObitDConClassInfo*)myClassInfo.ParentClass;
  ParentClass->ObitDConGetParms(inn, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Restore image when done? */
  ObitInfoListGetTest(in->info, "doRestore", &type, dim, &in->doRestore);

} /* end ObitDConCleanImageGetParms */

/**
 * Subtract components and generate new residual image(s).
 * Grids and FFT components, multiplies by transfer function
 * (FFT of beam), reverse FFTs and subtracts from image.
 * \param in   The object to deconvolve
 * \param err Obit error stack object.
 */
void ObitDConCleanImageSub(ObitDConCleanImage *in, ObitErr *err)
{
  ObitIOCode retCode;
  ObitTable *tempTable=NULL;
  ObitTableCC *CCTable = NULL;
  ObitImage *image=NULL;
  ObitImageDesc *imDesc = NULL;
  ObitIOSize IOsize = OBIT_IO_byPlane;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong  blc[IM_MAXDIM], trc[IM_MAXDIM];
  olong i, first, last, ncomp, ver, ndim, naxis[2];
  gchar *tabType = "AIPS CC";
  ofloat gparm[3];
  ObitFArray *grid = NULL;
  ObitCArray *uvGrid = NULL;
  gchar *routine = "ObitDConCleanImageSub";

  /* DEBUG
     ObitFArray *tempFArray = NULL; */
  /* END DEBUG */

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitDConCleanImageIsA(in));

  /* which Image? */
  image = in->mosaic->images[0];
  imDesc = image->myDesc;
  
  /* Full field, correct plane */
  dim[0] = IM_MAXDIM;
  blc[0] = blc[1] = 1;
  for (i=0; i<IM_MAXDIM-2; i++) blc[i+2] = in->plane[i];
  ObitInfoListPut (image->info, "BLC", OBIT_long, dim, blc, err); 
  trc[0] = trc[1] = 0;
  for (i=0; i<IM_MAXDIM-2; i++) trc[i+2] = in->plane[i];
  ObitInfoListPut (image->info, "TRC", OBIT_long, dim, trc, err); 
  dim[0] = 1;
  ObitInfoListPut (image->info, "IOBy", OBIT_long, dim, &IOsize, err);
  
  /* Open Image */
  retCode = ObitImageOpen (image, OBIT_IO_ReadWrite, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_msg (err, routine, in->name);
  
  /* Get CC table */
  ver = in->CCver;
  tempTable = newObitImageTable (image, OBIT_IO_ReadWrite, tabType, &ver, err);
  if ((tempTable==NULL) || (err->error)) Obit_traceback_msg (err, routine, in->name);
  CCTable = ObitTableCCConvert(tempTable);
  tempTable = ObitTableUnref(tempTable);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  in->CCver = ver;  /* save if defaulted (0) */
  
  /* Grid components - do all not already done */
  first = in->startComp;
  last = 0;
  retCode = ObitTableCCUtilGrid (CCTable, 1, &first, &last, FALSE, 1.0, 0.0, 1.0e20,
				 imDesc, &grid, gparm, &ncomp, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_msg (err, routine, in->name);
  
  /* Keep track of what done */
  in->startComp = last + 1;
  in->endComp = last;

  /* Free CC table */
  CCTable = ObitTableCCUnref(CCTable);
  
  /* FFT to image plane */
  ObitFArray2DCenter (grid); /* Swaparoonie to FFT order */
  /* Make Output of FFT if needed */
  ndim = 2;
  naxis[0] = 1+grid->naxis[1]/2; naxis[1] = grid->naxis[0]; 
  if (uvGrid) uvGrid = ObitCArrayRealloc(uvGrid, ndim, naxis);
  else uvGrid = ObitCArrayCreate ("FFT output", ndim, naxis);
  /* FFT */
  ObitFFTR2C (in->forFFT, grid, uvGrid);
  /* Put the center at the center */
  ObitCArray2DCenter (uvGrid);
  
 
  /* DEBUG */
  /* tempFArray = ObitCArrayMakeF(uvGrid);*/  /* Temp FArray */
  /*ObitCArrayImag (uvGrid, tempFArray); */  /* Get imaginary part */
  /*ObitImageUtilArray2Image ("DbuguvGridImag.fits", 1, tempFArray, err);*/
  /* tempFArray = ObitFArrayUnref(tempFArray);*/   /* delete temporary */
  /* END DEBUG */
  
  /* Multiply by Transfer function  */
  ObitCArrayMul (uvGrid, in->transfer, uvGrid);
  
  /* DEBUG */
  /*tempFArray = ObitCArrayMakeF(in->transfer);*/   /* Temp FArray */
  /*ObitCArrayReal (in->transfer, tempFArray);  */   /* Get real part */
  /*ObitImageUtilArray2Image ("DbugTransferReal.fits", 1, tempFArray, err);*/
  /*tempFArray = ObitFArrayUnref(tempFArray);*/   /* delete temporary */
  /* END DEBUG */
  
  /* FFT back to image */
  ObitCArray2DCenter (uvGrid); /* Swaparoonie to FFT order */
  ObitFFTC2R (in->revFFT, uvGrid, grid);
  /* Put the center at the center */
  ObitFArray2DCenter (grid);
  
  /* DEBUG */
  /* ObitImageUtilArray2Image ("DbugConvCC.fits", 1, grid, err);*/  
  /* if (err->error) Obit_traceback_msg (err, routine, in->name);*/  
  /*fprintf(stderr,"After WritingGriddedImage\n"); */
  /* END DEBUG */

 /* read residuals */
  retCode = ObitImageRead (image, image->image->array, err);
  if (err->error) Obit_traceback_msg (err, routine, image->name);
  
  /* Subtract convolved components */
  ObitFArraySub (image->image, grid, grid);
  
  /* Close Image and reopen to get correctly positioned */
  retCode = ObitImageClose (image, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_msg (err, routine, in->name);
  
  /* Open Image */
  retCode = ObitImageOpen (image, OBIT_IO_ReadWrite, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_msg (err, routine, in->name);

  /* rewrite image */
  retCode = ObitImageWrite (image, grid->array, err);
  if (err->error) Obit_traceback_msg (err, routine, image->name);
  
   /* Close Image */
  retCode = ObitImageClose (image, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_msg (err, routine, in->name);
  
 /* Free image memory */
  image->image = ObitFArrayUnref(image->image);
  
  /* Cleanup */
  grid   = ObitFArrayUnref(grid);
  uvGrid = ObitCArrayUnref(uvGrid);

} /* end ObitDConCleanImageSub */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitDConCleanImageClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitDConCleanImageClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitDConCleanImageClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitDConCleanImageClassInfoDefFn (gpointer inClass)
{
  ObitDConCleanImageClassInfo *theClass = (ObitDConCleanImageClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitDConCleanImageClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitDConCleanImageClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitDConCleanImageGetClass;
  theClass->newObit       = (newObitFP)newObitDConCleanImage;
  theClass->ObitCopy      = (ObitCopyFP)ObitDConCleanImageCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitDConCleanImageClear;
  theClass->ObitInit      = (ObitInitFP)ObitDConCleanImageInit;
  theClass->ObitDConGetParms        = (ObitDConGetParmsFP)ObitDConCleanImageGetParms;
  theClass->ObitDConDeconvolve      = (ObitDConDeconvolveFP)ObitDConCleanImageDeconvolve;
  theClass->ObitDConCleanPixelStats = (ObitDConCleanPixelStatsFP)ObitDConCleanPixelStats;
  theClass->ObitDConCleanSelect     = (ObitDConCleanSelectFP)ObitDConCleanSelect;
  theClass->ObitDConCleanSub        = (ObitDConCleanSubFP)ObitDConCleanImageSub;
  theClass->ObitDConCleanRestore    = (ObitDConCleanRestoreFP)ObitDConCleanRestore;
  theClass->ObitDConCleanFlatten    = NULL;  /* Not relevant */

} /* end ObitDConCleanImageClassDefFn */


/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitDConCleanImageInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitDConCleanImage *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->transfer  = NULL;
  in->forFFT    = NULL;
  in->revFFT    = NULL;
  in->doRestore = TRUE;

} /* end ObitDConCleanImageInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitDConCleanImage* cast to an Obit*.
 */
void ObitDConCleanImageClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitDConCleanImage *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  in->transfer = ObitCArrayUnref(in->transfer);
  in->forFFT   = ObitFFTUnref(in->forFFT);
  in->revFFT   = ObitFFTUnref(in->revFFT);
 
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitDConCleanImageClear */

/**
 * Create FFT objects and FFT the Beam to the Transfer function
 * \param in  The object to deconvolve
 * \param err Obit error stack object.
 */
static void  GetXfer (ObitDConCleanImage *in, ObitErr *err)
{
  ObitIOCode retCode;
  ObitImage *Beam = NULL;
  ObitIOSize IOsize = OBIT_IO_byPlane;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong  i, blc[IM_MAXDIM], trc[IM_MAXDIM];
  ofloat norm;
  olong naxis[2], ndim;
  gchar *routine = "GetXfer";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* If transfer function already exists, simply return */
  if (in->transfer) return;

  /* Get Beam */
  Beam = (ObitImage*)(in->mosaic->images[0]->myBeam);

  /* Set output to full image, plane at a time */
  dim[0] = IM_MAXDIM;
  blc[0] = blc[1] = blc[2] = blc[3] = blc[4] = blc[5] = 1;
  trc[0] = trc[1] = trc[2] = trc[3] = trc[4] = trc[5] = 0;
  /* multiplane? */
  if (Beam->myDesc->inaxes[2]>1) {
    for (i=0; i<IM_MAXDIM-2; i++) trc[i+2] = in->plane[i];
    for (i=0; i<IM_MAXDIM-2; i++) blc[i+2] = in->plane[i];
  }
  ObitInfoListPut (Beam->info, "BLC", OBIT_long, dim, blc, err); 
  ObitInfoListPut (Beam->info, "TRC", OBIT_long, dim, trc, err); 
  dim[0] = 1;
  ObitInfoListPut (Beam->info, "IOBy", OBIT_long, dim, &IOsize, err);
 
  retCode = ObitImageOpen (Beam, OBIT_IO_ReadOnly, err);
  if (err->error) Obit_traceback_msg (err, routine, Beam->name);

  retCode = ObitImageRead (Beam, Beam->image->array, err);
  if (err->error) Obit_traceback_msg (err, routine, Beam->name);

  retCode = ObitImageClose (Beam, err);
  if (err->error) Obit_traceback_msg (err, routine, Beam->name);

  /* Make FFTs */
  dim[0] = Beam->image->naxis[0]; dim[1] = Beam->image->naxis[1];
  in->forFFT = newObitFFT("FFT:FTImage", OBIT_FFT_Forward, 
			  OBIT_FFT_HalfComplex, 2, dim);
  in->revFFT = newObitFFT("FFT:FTuv", OBIT_FFT_Reverse, 
			  OBIT_FFT_HalfComplex, 2, dim);

  /* Make transfer */
  naxis[0] =  1+Beam->image->naxis[0]/2; naxis[1] = Beam->image->naxis[1];
  ndim = 2;
  in->transfer = ObitCArrayCreate ("FFT output", ndim, naxis);

  /* FFT */
  ObitFArray2DCenter (Beam->image); /* Swaparoonie to FFT order */
  ObitFFTR2C (in->forFFT, Beam->image, in->transfer);
  /* Put the center at the center */
  ObitCArray2DCenter (in->transfer);

  /* Normalize */
  norm = ((ofloat)Beam->image->naxis[0]) * ((ofloat)Beam->image->naxis[1]);
  norm = 1.0 / norm;
  ObitCArraySMul (in->transfer, norm);

  /* Free Image array? */
  Beam->image = ObitFArrayUnref(Beam->image);

} /* end GetXfer */
