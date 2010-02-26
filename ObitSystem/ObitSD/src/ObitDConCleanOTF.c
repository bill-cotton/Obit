/* $Id$  */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2005-2009                                          */
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

#include "ObitDConCleanOTF.h"
#include "ObitMem.h"
#include "ObitFFT.h"
#include "ObitTableCCUtil.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitDConCleanOTF.c
 * ObitDConCleanOTF class function definitions.
 * Image based CLEAN class.
 * This class is derived from the ObitDCon class.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitDConCleanOTF";

/** Function to obtain parent ClassInfo - ObitDConClean */
static ObitGetClassFP ObitParentGetClass = ObitDConCleanGetClass;

/**
 * ClassInfo structure ObitDConCleanOTFClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitDConCleanOTFClassInfo myClassInfo = {FALSE};

/*--------------- File Global Variables  ----------------*/


/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitDConCleanOTFInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitDConCleanOTFClear (gpointer in);

/** Read Beam into Beam patch */
static void ReadBPOTF (ObitDConCleanOTF* in, ObitErr *err);

/** Private: Set Class function pointers. */
static void ObitDConCleanOTFClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitDConCleanOTF* newObitDConCleanOTF (gchar* name)
{
  ObitDConCleanOTF* out;
  /*gchar *routine = "newObitDConCleanOTF";*/

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitDConCleanOTFClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitDConCleanOTF));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitDConCleanOTFInit((gpointer)out);

 return out;
} /* end newObitDConCleanOTF */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitDConCleanOTFGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitDConCleanOTFClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitDConCleanOTFGetClass */

/**
 * Make a deep copy of an ObitDConCleanOTF.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitDConCleanOTF* ObitDConCleanOTFCopy  (ObitDConCleanOTF *in, 
					     ObitDConCleanOTF *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  gboolean oldExist;
  gchar *outName;
  gchar *routine = "ObitDConCleanOTFCopy";

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
    out = newObitDConCleanOTF(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, out);

  /*  copy this class */
  out->dirty = ObitImageUnref(out->dirty);
  out->dirty = ObitImageRef(in->dirty);
  out->beam  = ObitImageUnref(out->beam);
  out->beam  = ObitImageRef(in->beam);
  out->clean = ObitImageUnref(out->clean);
  out->clean = ObitImageRef(in->clean);

  return out;
} /* end ObitDConCleanOTFCopy */

/**
 * Make a copy of a object but do not copy the actual data
 * This is useful to create an DConCleanOTF similar to the input one.
 * \param in  The object to copy
 * \param out An existing object pointer for output, must be defined.
 * \param err Obit error stack object.
 */
void ObitDConCleanOTFClone  (ObitDConCleanOTF *in, ObitDConCleanOTF *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  gchar *routine = "ObitDConCleanOTFClone";

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
  out->dirty = ObitImageUnref(out->dirty);
  out->dirty = ObitImageRef(in->dirty);
  out->beam  = ObitImageUnref(out->beam);
  out->beam  = ObitImageRef(in->beam);
  out->clean = ObitImageUnref(out->clean);
  out->clean = ObitImageRef(in->clean);

} /* end ObitDConCleanOTFClone */

/**
 * Creates an ObitDConCleanOTF 
 * \param name   An optional name for the object.
 * \param dirty  Dirty image
 * \param beam   Dirty beam, actual raw instrumental response
 * \param clean  Image to receive Clean residual or restored image
 * \param err    Obit error stack object.
 * \return the new object.
 */
ObitDConCleanOTF* 
ObitDConCleanOTFCreate (gchar* name, ObitImage *dirty, ObitImage *beam,  
			ObitImage *clean, ObitErr *err)
{
  olong nfield;
  ObitDConCleanOTF* out=NULL;
  gchar *routine = "ObitDConCleanOTFCreate";

 /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return out;
  g_assert (ObitImageIsA(dirty));
  if (beam) g_assert (ObitImageIsA(beam));
  g_assert (ObitImageIsA(clean));

  /* Create basic structure */
  out = newObitDConCleanOTF (name);

  /* save inputs */
  out->dirty   = ObitImageRef(dirty);
  out->beam    = ObitImageRef(beam);
  out->clean   = ObitImageRef(clean);

  /* Arrays - including those in parent classes */
  nfield =  1;
  out->gain    = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean Loop gain");
  out->minFlux = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean minFlux");
  out->factor  = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean factor");
  out->currentFields = ObitMemAlloc0Name((nfield+1)*sizeof(olong),"Current fields");
  out->BeamPatches   = ObitMemAlloc0Name(nfield*sizeof(ObitFArray*),"Beam patch");
  out->numCurrentField  = 1;
  out->currentFields[0] = 1;
  out->currentFields[1] = 0;

  /* Copy dirty image to clean image and use as residual */
  out->clean = ObitImageCopy (out->dirty, out->clean, err);
  if (err->error) Obit_traceback_val (err, routine, name, out);

  /* Image mosaic for window */
  out->mosaic = newObitImageMosaic (name, 1);
  ObitImageMosaicSetImage (out->mosaic, 0, out->clean, err);
  if (err->error) Obit_traceback_val (err, routine, name, out);

  /* Window object */
  out->window = ObitDConCleanWindowCreate ("CleanWindow", out->mosaic, err);
  if (err->error) Obit_traceback_val (err, routine, name, out);

  return out;
} /* end ObitDConCleanOTFCreate */

/**
 * Do deconvolution, uses function on class pointer
 * Does final flatten if FullField member of mosaic member is defined.
 * CLEAN control parameters are in the ObitInfoList member:
 * \li "Niter"   OBIT_long scalar   = Maximum number of CLEAN iterations
 * \li "maxPixel" OBIT_long scalar  = Maximum number of residuals [def 20000]
 * \li "minPatch" OBIT_long scalar  = Minimum beam patch in pixels [def 50]
 * \li "BMaj"    OBIT_float scalar = Restoring beam major axis (deg)
 * \li "BMin"    OBIT_float scalar = Restoring beam minor axis (deg)
 * \li "BPA"     OBIT_float scalar = Restoring beam position angle (deg)
 * \li "CCVer"   OBIT_long array  = CLEAN table version per field
 * \li "Gain"    OBIT_float array  = CLEAN loop gain per field
 * \li "minFlux" OBIT_float array  = Minimun flux density (Jy)  per field
 * \li "Factor"  OBIT_float array  = CLEAN depth factor per field
 * \li "Plane"   OBIT_long array    = Plane being processed, 1-rel indices of axes 3-?
 * \param in   The object to deconvolve
 * \param err Obit error stack object.
 */
void ObitDConCleanOTFDeconvolve (ObitDCon *inn, ObitErr *err)
{
  ObitDConCleanOTF *in;
  gboolean done;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong numPix;
  olong maxPix;
  const ObitDConCleanOTFClassInfo *inClass;
  gchar *routine = "ObitDConCleanOTFDeconvolve";

  /* Cast input to this type */
  in = (ObitDConCleanOTF*)inn;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitDConCleanOTFIsA(in));

  inClass = (ObitDConCleanOTFClassInfo*)in->ClassInfo; /* class structure */

  /* Get parameters */
  inClass->ObitDConGetParms(inn, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* How many pixels selected? */
  numPix = ObitDConCleanWindowCount (in->window, 1, err);
  numPix += 10; /* for good measure */
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* (Re)Create Pixel List if needed */
  if ((in->Pixels) && (numPix>in->Pixels->maxPixel))  
    in->Pixels = ObitDConCleanPxListUnref(in->Pixels);
  if (!in->Pixels) {
    in->Pixels = ObitDConCleanPxListCreate("Pixel List", in->mosaic, 
					   numPix, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
  }

  /* Copy control info to PixelList */
  ObitInfoListCopyData(in->info, in->Pixels->info);
  dim[0] = dim[1] = 1;
  ObitInfoListAlwaysPut(in->info, "minPatch", OBIT_long, dim, &in->beamPatchSize);
  maxPix = numPix;
  ObitInfoListAlwaysPut(in->info, "maxPixel", OBIT_long, dim, &maxPix);

  /* Reset/Init Pixel list*/
  ObitDConCleanPxListReset (in->Pixels, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Read Beam patch if needed*/
  if (!in->BeamPatches[0]) ReadBPOTF (in, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /*Only one pass deconvolution needed */
  /* Pick components */
  done = inClass->ObitDConCleanSelect((ObitDConClean*)in, NULL, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  
  /* Make residual image */
  inClass->ObitDConCleanSub((ObitDConClean*)in, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Scale CCs? */
  if (in->doScaleCC) 
    inClass->ObitDConCleanOTFScaleCC(in, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  
  /* Restore */
  if (in->doRestore) 
    inClass->ObitDConCleanRestore((ObitDConClean*)in, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

} /* end ObitDConCleanOTFDeconvolve */

/**
 * Read any base class parameters and then
 * read CLEAN control parameters from the ObitInfoList member:
 * From Parent classes:
 * \li "Niter"   OBIT_long scalar   = Maximum number of CLEAN iterations
 * \li "CCVer"   OBIT_long array    = CLEAN table version
 * \li "Gain"    OBIT_float array  = CLEAN loop gain
 * \li "minFlux" OBIT_float array  = Minimum flux density (Jy)
 * \li "Factor"  OBIT_float array  = CLEAN depth factor
 * \li "Plane"   OBIT_long array    = Plane being processed, 1-rel indices of axes 3-?
 * \li "doRestore" OBIT_bool       = Restore image when done? [def TRUE]
 * 
 * This Class:
 * \li "Patch"   OBIT_long scalar   = Beam patch in pixels [def 100]
 * \li "BeamSize"OBIT_float scalar = Restoring beam FWHM (deg)
 * \li "doScaleCC" OBIT_bool scalar = Scale CCs by ratio of dirty to CLEAN beam areas?
 * \param in  The CLEAN object as base class
 * \param err Obit error stack object.
 */
void  ObitDConCleanOTFGetParms (ObitDCon *inn, ObitErr *err)
{
  ObitDConCleanOTF *in = (ObitDConCleanOTF*)inn;  /* as this class */
  ObitDConClassInfo *ParentClass;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  /*  olong itemp;
      union ObitInfoListEquiv InfoReal; */
  gchar *routine = "ObitDConCleanOTFGetParms";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Read any parent class parameters */
  ParentClass = (ObitDConClassInfo*)myClassInfo.ParentClass;
  ParentClass->ObitDConGetParms(inn, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Use Beam size */
  in->cleanSize = 1.0/3600.0;
  ObitInfoListGetTest(in->info, "BeamSize", &type, (gint32*)dim, &in->cleanSize);

  /* beam patch */
  in->beamPatchSize = 100;
  ObitInfoListGetTest(in->info, "Patch", &type, (gint32*)dim, &in->beamPatchSize);

  /* Restore image when done? */
  in->doRestore = TRUE;
  ObitInfoListGetTest(in->info, "doRestore", &type, dim, &in->doRestore);

  /* Scale CCs by ratio of dirty to CLEAN beam areas? */
  in->doScaleCC = TRUE;
  ObitInfoListGetTest(in->info, "doScaleCC", &type, dim, &in->doScaleCC);

} /* end ObitDConCleanOTFGetParms */

/**
 * Make residual image from Cleaned PxList
 * Extract pixel values from the PxList and replace in residual image
 * \param in   The object to deconvolve
 * \param err Obit error stack object.
 */
void ObitDConCleanOTFSub(ObitDConClean *inn, ObitErr *err)
{
  ObitDConCleanOTF *in;
  ObitIOSize IOsize = OBIT_IO_byPlane;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong  i, blc[IM_MAXDIM], trc[IM_MAXDIM];
  olong x, y, nx, pos[2];
  ofloat *pixels=NULL;
  ObitImage *image=NULL;
  ObitFArray *residual=NULL;
  gchar *routine = "ObitDConCleanOTFSub";

  /* Cast input to this type */
  in = (ObitDConCleanOTF*)inn;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitDConCleanOTFIsA(in));

  /* which Image? */
  image = in->clean;
  
  /* Full field, correct plane */
  dim[0] = IM_MAXDIM;
  for (i=0; i<IM_MAXDIM; i++) blc[i] = 1;
  for (i=0; i<IM_MAXDIM; i++) trc[i] = 0;
  for (i=0; i<IM_MAXDIM-2; i++) blc[i+2] = in->plane[i];
  ObitInfoListPut (image->info, "BLC", OBIT_long, dim, blc, err); 
  for (i=0; i<IM_MAXDIM-2; i++) trc[i+2] = in->plane[i];
  ObitInfoListPut (image->info, "TRC", OBIT_long, dim, trc, err); 
  dim[0] = 1;
  ObitInfoListPut (image->info, "IOBy", OBIT_long, dim, &IOsize, err);
  
  /* Read clean ((copy or dirty) image */
  in->clean->extBuffer = FALSE;
  ObitImageOpen (image, OBIT_IO_ReadWrite, err); 
  ObitImageRead (image, NULL, err);
  ObitImageClose (image, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* copy pixels */
  residual = ObitFArrayCopy (in->clean->image, residual, err);

  /* Replace pixels with those from the PxList */
  pos[0] = pos[1] = 0;
  pixels = ObitFArrayIndex(residual, pos);
  nx = in->clean->image->naxis[0];
  for (i=0; i<in->Pixels->nPixel; i++) {
    x = in->Pixels->pixelX[i];
    y = in->Pixels->pixelY[i];
    pixels[x+y*nx] = in->Pixels->pixelFlux[i];
  }

  /* Write residual image (in clean) */
  in->clean->extBuffer = TRUE;
  ObitImageOpen (in->clean, OBIT_IO_ReadWrite, err); 
  ObitImageWrite (in->clean, residual->array, err) ;
  ObitImageClose (in->clean, err);
  in->clean->extBuffer = FALSE;
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Cleanup */
  residual = ObitFArrayUnref(residual);
} /* end ObitDConCleanOTFSub */

/**
 * Restore Clean components to residual
 * Following Parameters used from info member:
 * \li noResid    OBIT_bool scalar If defined and TRUE, set residuals to Zero
 * \li doScale    OBIT_bool scalar If defined and TRUE, Scale residuals by 
 *                ratio of beam areas, def TRUE
 * \li Niter      OBIT_int scalar  Maximum number of iterations [def 200].
 * \li BeamSize   OBIT_float scalar CLEAN restoring beam FWHM in asec [def 3.0].
 *                If < 0 then don't restore.
 *
 * \param in   Object with OTFCLEAN structures.
 * \param err  Error stack
 */
void ObitDConCleanOTFRestore (ObitDConClean *inn, ObitErr *err)
{
  ObitDConCleanOTF *in;
  olong iRow;
  ObitIOSize IOsize = OBIT_IO_byPlane;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong  blc[IM_MAXDIM], trc[IM_MAXDIM];
  olong  i, pos[2], loop, beamCen[2], nrow;
  ofloat DBArea, CBArea, scale, fmax, cleanSize;
  ObitTableCC *CCTable=NULL;
  ObitTableCCRow *row=NULL;
  ObitFArray *residual=NULL, *cleanBeam=NULL;
  gboolean noresid = FALSE, doScale = TRUE;
  gchar *tname;
  gchar *routine = "ObitDConCleanOTFRestore";

  /* Cast input to this type */
  in = (ObitDConCleanOTF*)inn;

 /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Is restoration desired? */
  if (in->cleanSize<0.0) return;

  /* Convert restoring beam size into cells */
  cleanSize = in->cleanSize / (fabs(in->clean->myDesc->cdelt[0]));

  /* Generate clean beam - clone from BeamPatch */
  if (!in->BeamPatches[0]) ReadBPOTF(in, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* center of beam = peak */
  fmax = ObitFArrayMax (in->BeamPatches[0], beamCen);

  /* Copy */
  cleanBeam = ObitFArrayCopy(in->BeamPatches[0], cleanBeam, err);
  if (err->error) Obit_traceback_msg (err, routine, in->BeamPatches[0]->name);

  /* Get dirty beam area */
  DBArea = ObitFArraySum (in->BeamPatches[0]);

  /* Replace with Gaussian */
  ObitFArray2DCGauss (cleanBeam, beamCen, cleanSize);

  /* Get CLEAN beam area */
  CBArea = ObitFArraySum (cleanBeam);

  /* Read residual image */
  /* Full field, correct plane */
  dim[0] = IM_MAXDIM;
  blc[0] = blc[1] = 1;
  for (i=0; i<IM_MAXDIM-2; i++) blc[i+2] = in->plane[i];
  ObitInfoListPut (in->clean->info, "BLC", OBIT_long, dim, blc, err); 
  trc[0] = trc[1] = 0;
  for (i=0; i<IM_MAXDIM-2; i++) trc[i+2] = in->plane[i];
  ObitInfoListPut (in->clean->info, "TRC", OBIT_long, dim, trc, err); 
  dim[0] = 1;
  ObitInfoListPut (in->clean->info, "IOBy", OBIT_long, dim, &IOsize, err);
  in->clean->extBuffer = FALSE;
  ObitImageOpen (in->clean, OBIT_IO_ReadWrite, err); 
  ObitImageRead (in->clean, NULL, err);
  ObitImageClose (in->clean, err); 
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Save clean to residual */
  residual = ObitFArrayCopy (in->clean->image, residual, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  /* Want residuals */
  noresid = FALSE;
  ObitInfoListGetTest(in->info, "noResid", &type, dim, &noresid);
  /* Want to scale residuals? */
  doScale = TRUE;
  ObitInfoListGetTest(in->info, "doScale", &type, dim, &doScale);

  /* Scale residuals by ratio of dirty to clean beam areas to get into
     units of restored components. */
  if (CBArea>0.0) scale = DBArea / CBArea;
  else scale = 1.0;  /* something went wrong */
  if (!doScale) scale = 1.0;  /* Scale */
  if (noresid) scale = 0.0;  /* want residuals? */
  Obit_log_error(err, OBIT_InfoErr, "Scaling residuals by %f", scale);
  ObitFArraySMul (residual, scale);

  /* Initialize CC table on CLEAN image */
  tname = g_strconcat ("Clean table for: ",in->name, NULL);
  CCTable =  newObitTableCCValue (tname, (ObitData*)in->clean, &in->CCver, 
				  OBIT_IO_ReadWrite, 0, err);
  g_free (tname);
  
  ObitTableCCOpen (CCTable, OBIT_IO_ReadWrite, err); 
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  row = newObitTableCCRow(CCTable);

  /* Restoration loop */
  nrow = CCTable->myDesc->nrow;
  Obit_log_error(err, OBIT_InfoErr, "Restoring %d components", nrow); /* debug */
  for (loop=1; loop<=nrow; loop++) {

    /* Read CC table */
    iRow = loop;
    ObitTableCCReadRow (CCTable, iRow, row, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);

    /* Restore to residual */
    pos[0] = -0.5 + (row->DeltaX / in->clean->myDesc->cdelt[0]) + in->clean->myDesc->crpix[0];
    pos[1] = -0.5 + (row->DeltaY / in->clean->myDesc->cdelt[1]) + in->clean->myDesc->crpix[1];
    ObitFArrayShiftAdd (residual, pos, cleanBeam, beamCen, row->Flux, residual);
  } /* End Restoration loop */


  /* Close CC table */
  ObitTableCCClose (CCTable, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Write residual image */
  in->clean->extBuffer = TRUE;
  ObitImageOpen (in->clean, OBIT_IO_ReadWrite, err); 
  /* Update resolution on CLEAN image */
  in->clean->myDesc->beamMaj = in->cleanSize;
  in->clean->myDesc->beamMin = in->cleanSize;
  in->clean->myDesc->beamPA  = 0;
  in->clean->myDesc->niter   = nrow;
  ObitImageWrite (in->clean, residual->array, err) ;
  ObitImageClose (in->clean, err);
  in->clean->extBuffer = FALSE;
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Cleanup */
  CCTable   = ObitTableCCUnref(CCTable);
  row       = ObitTableCCRowUnref(row);
  cleanBeam = ObitFArrayUnref(cleanBeam);
  residual  = ObitFArrayUnref(residual);

} /* end ObitDConCleanOTFRestore */

/**
 * Scale Clean Components by ratio of dirty beam area to clean beam area
 * Following Parameters used from info member:
 *
 * \param in   Object with OTFCLEAN structures.
 * \param err  Error stack
 */
void ObitDConCleanOTFScaleCC (ObitDConCleanOTF *in, ObitErr *err)
{
  ObitFArray *cleanBeam=NULL;
  olong  iRow, loop, beamCen[2], nrow;
  ofloat DBArea, CBArea, scale, fmax, cleanSize;
  ObitTableCC *CCTable=NULL;
  ObitTableCCRow *row=NULL;
  gchar *tname;
  gchar *routine = "ObitDConCleanOTFScaleCC";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Convert restoring beam size into cells */
  cleanSize = in->cleanSize / (fabs(in->clean->myDesc->cdelt[0]));

  /* Generate clean beam - clone from BeamPatch */
  if (!in->BeamPatches[0]) ReadBPOTF(in, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* center of beam = peak */
  fmax = ObitFArrayMax (in->BeamPatches[0], beamCen);

  /* Get dirty beam area */
  DBArea = ObitFArraySum (in->BeamPatches[0]);

  /* Copy Dirty beam */
  cleanBeam = ObitFArrayCopy(in->BeamPatches[0], cleanBeam, err);
  if (err->error) Obit_traceback_msg (err, routine, in->BeamPatches[0]->name);

  /* Replace with CLEAN Gaussian */
  ObitFArray2DCGauss (cleanBeam, beamCen, cleanSize);

  /* Get CLEAN beam area */
  CBArea = ObitFArraySum (cleanBeam);

  /* Scaling factor */
  if (CBArea>0.0) scale = DBArea / CBArea;
  else scale = 1.0;  /* something went wrong */

  /* Initialize CC table on CLEAN image */
  tname = g_strconcat ("Clean table for: ",in->name, NULL);
  CCTable =  newObitTableCCValue (tname, (ObitData*)in->clean, &in->CCver, 
				  OBIT_IO_ReadWrite, 0, err);
  g_free (tname);
  
  ObitTableCCOpen (CCTable, OBIT_IO_ReadWrite, err); 
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  row = newObitTableCCRow(CCTable);

  /* Scaling loop */
  nrow = CCTable->myDesc->nrow;
  for (loop=1; loop<=nrow; loop++) {

    /* Read CC table */
    iRow = loop;
    ObitTableCCReadRow (CCTable, iRow, row, err);
    row->Flux *= scale;
    ObitTableCCWriteRow (CCTable, iRow, row, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
  } /* End Restoration loop */


  /* Close CC table */
  ObitTableCCClose (CCTable, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Cleanup */
  CCTable   = ObitTableCCUnref(CCTable);
  row       = ObitTableCCRowUnref(row);
  cleanBeam = ObitFArrayUnref(cleanBeam);

} /* end ObitDConCleanOTFScaleCC */

/**
 * Select/subtract components from PxList
 * \param in   The object to deconvolve
 * \param pixarray    If NonNULL use instead of the flux densities from the image file.
 *                    Array of ObitFArrays corresponding to fields in in->currentFields 
 * \param err         Obit error stack object.
 * \return TRUE if deconvolution is complete
 */
gboolean ObitDConCleanOTFSelect(ObitDConClean *inn, ObitFArray **pixarray, ObitErr *err)
{
  ObitDConCleanOTF *in;
  gboolean done = FALSE;
  olong numPix;
  olong fields[] = {1,0};
  gchar *routine = "ObitDConCleanOTFSelect";

  /* Cast input to this type */
  in = (ObitDConCleanOTF*)inn;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return done;
  g_assert (ObitDConCleanOTFIsA(in));

   /* How many pixels selected? */
  numPix = ObitDConCleanWindowCount (in->window, 1, err);
  numPix += 10; /* for good measure */
  if (err->error) Obit_traceback_val (err, routine, in->name, done);

  /* (Re)Create Pixel List if needed */
  if ((in->Pixels) && (numPix>in->Pixels->maxPixel))  
    ObitDConCleanPxListResize (in->Pixels, numPix, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, done);
  
 /* Load PxList */
  ObitDConCleanPxListUpdate (in->Pixels, fields, 0, 0.0, in->autoWinFlux,
			     in->window, in->BeamPatches, pixarray, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, done);

  /* Clean */
  ObitDConCleanPxListCLEAN (in->Pixels, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, done);

  return done;
} /* end ObitDConCleanOTFSelect */

/**
 * Get image and beam statistics 
 * \param in   The object to deconvolve
 * \param pixarray    If NonNULL use instead of the flux densities from the image file.
 * \param err Obit error stack object.
 */
void ObitDConCleanOTFPixelStats(ObitDConClean *in, ObitFArray **pixarray, 
				ObitErr *err)
{
  const ObitDConCleanClassInfo *inClass;
  gchar *routine = "ObitDConCleanOTFPixelStats";

 inClass = (ObitDConCleanClassInfo*)in->ClassInfo; /* class structure */

  /* Adjust window if autoWindow */
  if (in->autoWindow) 
    inClass->ObitDConCleanAutoWindow (in, in->currentFields, pixarray, err);
  else 
    in->autoWinFlux = -1.0e20; 
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  return;
} /* end ObitDConCleanPixelStat */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitDConCleanOTFClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitDConCleanOTFClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitDConCleanOTFClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitDConCleanOTFClassInfoDefFn (gpointer inClass)
{
  ObitDConCleanOTFClassInfo *theClass = (ObitDConCleanOTFClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitDConCleanOTFClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitDConCleanOTFClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitDConCleanOTFGetClass;
  theClass->ObitClear     = (ObitClearFP)ObitDConCleanOTFClear;
  theClass->ObitInit      = (ObitInitFP)ObitDConCleanOTFInit;
  theClass->newObit       = (newObitFP)newObitDConCleanOTF;
  theClass->ObitCopy      = (ObitCopyFP)ObitDConCleanOTFCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitDConGetParms        = (ObitDConGetParmsFP)ObitDConCleanOTFGetParms;
  theClass->ObitDConDeconvolve      = (ObitDConDeconvolveFP)ObitDConCleanOTFDeconvolve;
  theClass->ObitDConCleanPixelStats = (ObitDConCleanPixelStatsFP)ObitDConCleanOTFPixelStats;
  theClass->ObitDConCleanSelect     = (ObitDConCleanSelectFP)ObitDConCleanOTFSelect;
  theClass->ObitDConCleanSub        = (ObitDConCleanSubFP)ObitDConCleanOTFSub;
  theClass->ObitDConCleanRestore    = (ObitDConCleanRestoreFP)ObitDConCleanOTFRestore;
  theClass->ObitDConCleanFlatten    = NULL;  /* Not relevant */
  theClass->ObitDConCleanOTFScaleCC = (ObitDConCleanOTFScaleCCFP)ObitDConCleanOTFScaleCC;

} /* end ObitDConCleanOTFClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitDConCleanOTFInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitDConCleanOTF *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->dirty  = NULL;
  in->beam   = NULL;
  in->clean  = NULL;
  in->weight = NULL;
  in->doRestore = TRUE;

} /* end ObitDConCleanOTFInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitDConCleanOTF* cast to an Obit*.
 */
void ObitDConCleanOTFClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitDConCleanOTF *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  in->dirty  = ObitImageUnref(in->dirty);
  in->beam   = ObitImageUnref(in->beam);
  in->clean  = ObitImageUnref(in->clean);
  in->weight = ObitImageUnref(in->weight);
 
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitDConCleanOTFClear */

/**
 * Read Beam patch into BeamPatches[0]
 * The beam patch is symmetric about the center position allowing the 
 * beam itself not to be symmetric.
 * If the beam member is NULL, the BeamPatch is created and filled with the 
 * Gaussian size given in the dirty image header.
 * \param in   The object to deconvolve
 * \param err  Obit error stack object.
 */
static void ReadBPOTF (ObitDConCleanOTF* in, ObitErr *err)
{
  ObitIOCode retCode;
  ObitIOSize IOsize = OBIT_IO_byPlane;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong  blc[IM_MAXDIM], trc[IM_MAXDIM];
  olong  ablc[2], atrc[2], pos[2], naxis[2], cen[2];
  olong icenx, iceny, nx, ny, mxPatch;
  ofloat fmax, beamSize ;
  ObitImage *Beam;
  gchar *routine = "ObitDConCleanOTF:ReadBPOTF";
  
  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitDConCleanIsA(in));

  /* If no Beam defined in in->dirty, create a BeamPatch using resolution in dirty image */
  if (in->dirty->myBeam==NULL) {
    naxis[0] =  naxis[1] = 1 + 2*in->beamPatchSize;
    cen[0] = cen[1] = in->beamPatchSize;
    /* Create */
    in->BeamPatches[0] = ObitFArrayUnref(in->BeamPatches[0]);
    in->BeamPatches[0] = ObitFArrayCreate("BeamPatch", 2, naxis);
    /* Fill with Gaussian with size of that in dirty image */
    beamSize = in->dirty->myDesc->beamMaj / fabs(in->dirty->myDesc->cdelt[0]);
    ObitFArray2DCGauss (in->BeamPatches[0], cen, beamSize);
    return;
 } /* End create Gaussian Beam patch if no image given */
  
  /* Beam image */
  Beam = (ObitImage*)in->dirty->myBeam;
  
  /* Set output to full image, plane at a time */
  dim[0] = IM_MAXDIM;
  blc[0] = blc[1] = blc[2] = blc[3] = blc[4] = blc[5] = 1;
  ObitInfoListPut (Beam->info, "BLC", OBIT_long, dim, blc, err); 
  trc[0] = trc[1] = trc[2] = trc[3] = trc[4] = trc[5] = 0;
  ObitInfoListPut (Beam->info, "TRC", OBIT_long, dim, trc, err); 
  dim[0] = 1;
  ObitInfoListPut (Beam->info, "IOBy", OBIT_long, dim, &IOsize, err);
  
  retCode = ObitImageOpen (Beam, OBIT_IO_ReadOnly, err);
  if (err->error) Obit_traceback_msg (err, routine, Beam->name);

  retCode = ObitImageRead (Beam, Beam->image->array, err);
  if (err->error) Obit_traceback_msg (err, routine, Beam->name);

  /* Replace any blanks with zeroes */
  ObitFArrayDeblank (Beam->image, 0.0);

  /* Compute region of image to take */
  /* Center pixel - make 0-rel */
  nx = Beam->myDesc->inaxes[0];
  ny = Beam->myDesc->inaxes[1];
  /* center = peak */
  fmax = ObitFArrayMax (Beam->image, pos);
  icenx = pos[0];
  iceny = pos[1];

  /* Check that Dirty Map and Beam have the same cell spacing */
  if ((fabs(fabs(in->dirty->myDesc->cdelt[0]) - fabs(Beam->myDesc->cdelt[0])) >
       0.01*fabs(in->dirty->myDesc->cdelt[0])) ||
      (fabs(fabs(in->dirty->myDesc->cdelt[1]) - fabs(Beam->myDesc->cdelt[1])) >
       0.01*fabs(in->dirty->myDesc->cdelt[1]))) {
    Obit_log_error(err, OBIT_Error, "Dirty map and beam have different cell spacings");
    Obit_log_error(err, OBIT_Error, "Map %f %f Beam %f %f", 
		   3600.0*in->dirty->myDesc->cdelt[0], 3600.0*in->dirty->myDesc->cdelt[1],
		   3600.0*Beam->myDesc->cdelt[0],      3600.0*Beam->myDesc->cdelt[1]);
    return;
  }

 /* Check that center pretty close to 1.0 */
  if ((fmax<0.99) || (fmax>1.01)) {
      Obit_log_error(err, OBIT_Error, "%s Beam peak, %f not 1.0", 
		     routine, fmax);
      return;
  }

  /* How big can the patch be? */
  mxPatch = MIN (icenx-1, nx-icenx);
  mxPatch = MIN ( mxPatch, iceny-1);
  mxPatch = MIN ( mxPatch, ny-iceny);
  in->beamPatchSize =  MIN (in->beamPatchSize, mxPatch);

  /* Window as 0-rel */
  ablc[0] = icenx - in->beamPatchSize;
  atrc[0] = icenx + in->beamPatchSize;
  ablc[1] = iceny - in->beamPatchSize;
  atrc[1] = iceny + in->beamPatchSize;

  /* Save Beam patch */
  in->BeamPatches[0] = ObitFArrayUnref(in->BeamPatches[0]);
  in->BeamPatches[0] = ObitFArraySubArr(Beam->image, ablc, atrc, err);
  if (err->error) Obit_traceback_msg (err, routine, Beam->name);

  retCode = ObitImageClose (Beam, err);
  if (err->error) Obit_traceback_msg (err, routine, Beam->name);

  /* Free Image array? */
  Beam->image = ObitFArrayUnref(Beam->image);
} /* end ReadBPOTF */

