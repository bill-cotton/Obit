/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2010-2014                                          */
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
/*;  Correspondence concerning Obit should be addressed as follows:   */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/

#include "ObitDConCleanPxListMF.h"
#include "ObitMem.h"
#include "ObitImageMF.h"
#include "ObitImageMosaicMF.h"
#include "ObitSpectrumFit.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitDConCleanPxListMF.c
 * ObitDConCleanPxListMF class function definitions.
 * This class determines the pixel histogram of an image.
 * This class is derived from the Obit base class.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitDConCleanPxListMF";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitDConCleanPxListGetClass;

/**
 * ClassInfo structure ObitDConCleanPxListMFClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitDConCleanPxListMFClassInfo myClassInfo = {FALSE};

/*--------------- File Global Variables  ----------------*/

/*---------------Private structures----------------*/
/* CLEAN threaded function argument */
typedef struct {
  /* PixelList object */
  ObitDConCleanPxListMF    *PixelList;
  /* First (1-rel) record in otfdata buffer to process this thread */
  olong        first;
  /* Highest (1-rel) record in otfdata buffer to process this thread  */
  olong        last;
  /* thread number , >0 -> no threading  */
  olong        ithread;
  /* Pixel number of peak */
  olong       ipeak;
  /* Combined peak flux to subtract */
  ofloat       combPeak;
  /* Array of nSpec flux densities for component */
  ofloat*       Spec;
  /* Sum per spectral channel for SDI CLEAN */
  ofloat*      sum;
} CLEANFuncArg;

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitDConCleanPxListMFInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitDConCleanPxListMFClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitDConCleanPxListMFClassInfoDefFn (gpointer inClass);

/** Private: Get Spectral beam patches   */
static void GetSpecBeamPatch (ObitDConCleanPxListMF *in, olong ifld, 
			      ObitImageMF* image, ObitErr *err);

/** Private: Threaded CLEAN   */
static gpointer ThreadCLEAN (gpointer arg);

/** Private: Threaded SDI CLEAN */
static gpointer ThreadSDICLEAN (gpointer arg);

/** Private: Make arguments for Threaded CLEAN */
static olong MakeCLEANArgs (ObitDConCleanPxListMF *in, olong maxThread,
			    CLEANFuncArg ***args, ObitErr *err);

/** Private: Delete arguments for Threaded CLEAN */
static void KillCLEANArgs (olong nargs, CLEANFuncArg **args);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitDConCleanPxListMF* newObitDConCleanPxListMF (gchar* name)
{
  ObitDConCleanPxListMF* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitDConCleanPxListMFClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitDConCleanPxListMF));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitDConCleanPxListMFInit((gpointer)out);

 return out;
} /* end newObitDConCleanPxListMF */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitDConCleanPxListMFGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitDConCleanPxListMFClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitDConCleanPxListMFGetClass */

/**
 * Make a deep copy of an ObitDConCleanPxListMF.
 * Note: incomplete copy
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitDConCleanPxListMF* 
ObitDConCleanPxListMFCopy  (ObitDConCleanPxListMF *in, ObitDConCleanPxListMF *out, 
			  ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  gboolean oldExist;
  gchar *outName;
  /*gchar *routine = "ObitDConCleanPxListMFCopy";*/

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
    out = newObitDConCleanPxListMF(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  return out;
} /* end ObitDConCleanPxListMFCopy */

/**
 * Make a copy of a object but do not copy the actual data
 * This is useful to create an DConCleanPxListMF similar to the input one.
 * Note: incomplete copy
 * \param in  The object to copy
 * \param out An existing object pointer for output, must be defined.
 * \param err Obit error stack object.
 */
void ObitDConCleanPxListMFClone  (ObitDConCleanPxListMF *in, ObitDConCleanPxListMF *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  /*gchar *routine = "ObitDConCleanPxListMFClone";*/

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
} /* end ObitDConCleanPxListMFClone */

/**
 * Creates an ObitDConCleanPxListMF
 * Get spectral order from first image on mosaic
 * \param name     An optional name for the object.
 * \param mosaic   Image mosaic to be deconvolved.
 * \param uvdata   Data from which images are being made
 * \param maxPixel Maximum number of pixels allowed (dim. or arrays)
 * \return the new object.
 */
ObitDConCleanPxListMF* 
ObitDConCleanPxListMFCreate (gchar* name, ObitImageMosaic *mosaic,  
			     ObitUV *uvdata, olong maxPixel, ObitErr *err)
{
  ObitDConCleanPxListMF* out=NULL;
  olong i, nfield, maxOrder, nSpec;
  gchar *routine = "ObitDConCleanPxListMFCreate";

  /* error checks */
  if (err->error) return out;
  g_assert (ObitImageMosaicIsA(mosaic));
  /* Check input types */
  Obit_retval_if_fail((ObitImageMosaicMFIsA((ObitImageMosaicMF*)mosaic) && 
		       (ObitImageMFIsA((ObitImageMF*)mosaic->images[0]))), err, 
		      out,
		      "%s: Image mosaic or images not MF", routine);
  
  /* Create basic structure */
  out = newObitDConCleanPxListMF (name);
  
  /* How many spectral orders? - get maximum from mosaic */
  /* There may be different orders in different fields - deal with it */
  maxOrder = 0;
  nSpec    = 0;
  for (i=0; i<mosaic->numberImages; i++) {
    maxOrder = MAX (maxOrder, ((ObitImageMF*)mosaic->images[i])->maxOrder);
    nSpec    = MAX (nSpec, ((ObitImageMF*)mosaic->images[i])->nSpec);
  }
  out->nSpecTerm = maxOrder;
  out->curOrder  = maxOrder;
  out->maxOrder  = maxOrder;
  out->nSpec     = nSpec;

  /* Create arrays including parent class */
  out->maxPixel   = maxPixel;
  out->pixelX     = ObitMemAlloc0Name (maxPixel*sizeof(olong),  "PxListMF X pixel");
  out->pixelY     = ObitMemAlloc0Name (maxPixel*sizeof(olong),  "PxListMF Y pixel");
  out->pixelFld   = ObitMemAlloc0Name (maxPixel*sizeof(gshort), "PxListMF pixel field");
  out->pixelFlux  = ObitMemAlloc0Name (maxPixel*sizeof(ofloat), "PxListMF pixel Flux");
  out->specFreq   = ObitMemAlloc0Name (maxPixel*sizeof(odouble), "PxListMF specFreq");
  out->channelFlux= ObitMemAlloc0Name (maxPixel*sizeof(ofloat*), "PxListMF channelFlux");
  for (i=0; i<nSpec; i++) 
    out->channelFlux[i] = ObitMemAlloc0Name (maxPixel*sizeof(ofloat), "chnFlux");

  /* Save Image Mosaic reference */
  out->mosaic = ObitImageMosaicRef(mosaic);

  /* Per field */
  nfield      = mosaic->numberImages;
  out->nfield = nfield;
  out->iterField  = ObitMemAlloc0Name (nfield*sizeof(olong),  "PxListMF Clean CC count");
  out->CCver      = ObitMemAlloc0Name (nfield*sizeof(olong),  "PxListMF Clean CC version");
  out->fluxField  = ObitMemAlloc0Name (nfield*sizeof(ofloat), "PxListMF Clean Flux");
  out->circGaus   = ObitMemAlloc0Name (nfield*sizeof(ofloat), "PxListMF Gaussian");
  out->gain       = ObitMemAlloc0Name (nfield*sizeof(ofloat), "PxListMF Clean Loop gain");
  out->minFlux    = ObitMemAlloc0Name (nfield*sizeof(ofloat), "PxListMF Clean Mix flux");
  out->factor     = ObitMemAlloc0Name (nfield*sizeof(ofloat), "PxListMF Clean factor");
  out->CCTable    = ObitMemAlloc0Name (nfield*sizeof(ObitTableCC*), "PxListMF CC tables");
  out->CCRow      = ObitMemAlloc0Name (nfield*sizeof(ObitTableCCRow*), "PxListMF CC rows");
  out->BeamPatch   = ObitMemAlloc0Name (nfield*sizeof(ObitFArray*), "PxList BeamPatch");
  for (i=0; i<nfield; i++) {
    out->iterField[i] = 0;
    out->CCver[i]     = 0;
    out->fluxField[i] = 0.0;
    out->circGaus[i]  = mosaic->BeamTaper[i];
    out->gain[i]      = 0.1;
    out->minFlux[i]   = 0.0;
    out->factor[i]    = 0.0;
    out->CCTable[i]   = NULL;
    out->CCRow[i]     = NULL;
  }
  
  /* Frequency info */
  out->refFreq = ((ObitImageMF*)mosaic->images[0])->refFreq; /* Reference from image */

  /* Coarse channel frequencies from first image */
  for (i=0; i<nSpec; i++) 
    out->specFreq[i] = ((ObitImageMF*)mosaic->images[0])->specFreq[i];

  return out;
} /* end ObitDConCleanPxListMFCreate */

/**
 * Read CLEAN control parameters from the ObitInfoList member:
 * \li "Niter"   OBIT_long scalar   = Maximum number of CLEAN iterations
 * \li "Gain"    OBIT_float array  = CLEAN loop gain per field
 *                                   If only one given it is used for all.
 * \li "minFlux" OBIT_float array  = Minimun flux density (Jy)  per field
 *                                   If only one given it is used for all.
 * \li "Factor"  OBIT_float array  = CLEAN depth factor per field
 *                                   If only one given it is used for all.
 * \li "fGauss"  OBIT_float array  = Gaussian size (deg) per field
 * \li "CCVer"   OBIT_long array   = CLEAN table version per field
 *                                   If only one given it is used for all.
 * \li "prtLv"   OBIT_long         = message level  [def 2]
 *               0=none, 1=summary, 2=normal, higher numbers for diagnostics
 * \li "ccfLim"  OBIT_float        = Min. fraction of residual peak to CLEAN to
 *                                   clipped to [0.0,0.9]
 * \param inn The Pixel list CLEAN object
 * \param err Obit error stack object.
 */
void  ObitDConCleanPxListMFGetParms (ObitDConCleanPxList *inn, ObitErr *err)
{
  ObitDConCleanPxListClassInfo *ParentClass;
  ObitDConCleanPxListMF *in = (ObitDConCleanPxListMF*)inn;
  /*union ObitInfoListEquiv InfoReal; */
  gchar *routine = "ObitDConCleanPxListMFGetParms";

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Read any parent class parameters */
  ParentClass = (ObitDConCleanPxListClassInfo*)myClassInfo.ParentClass;
  ParentClass->ObitDConCleanPxListGetParms(inn, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

} /* end ObitDConCleanPxListMFGetParms */

/**
 * Resets CLEAN information
 * Makes sure all potential CC tables are instantiated
 * \param inn         The Pixel List object 
 * \param err Obit error stack object.
 */
void ObitDConCleanPxListMFReset (ObitDConCleanPxList *inn, ObitErr *err)
{
  olong i, ver=0;
  oint noParms;
  ObitDConCleanPxListMFClassInfo *myClass;
  ObitDConCleanPxListMF *in = (ObitDConCleanPxListMF*)inn;
  gchar *routine = " ObitDConCleanPxListMFReset";

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Control info */
  myClass = (ObitDConCleanPxListMFClassInfo*)in->ClassInfo;
  myClass->ObitDConCleanPxListGetParms(inn, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Number of components */
  in->currentIter = 0;
  in->totalFlux   = 0.0;
  for (i=0; i<in->nfield; i++) {
    in->iterField[i] = 0;
    in->fluxField[i] = 0.0;

    /* Get CC table if not defined */
    if (in->CCTable[i]==NULL) {
      /* Are these Gaussians or point? */
      noParms = in->nSpec;     /* Possible spectra */
      /* If adding spectra, also need type parameters even for point */
      if ((in->circGaus[i]>0.0) || (in->nSpec>1)) noParms += 4;
      else noParms = 0;
      ver = MAX (ver, in->CCver[i]);  /* Use last if not defined */
      in->CCTable[i] = ObitTableCCUnref(in->CCTable[i]);  /* Free old */
      in->CCTable[i] = 
	newObitTableCCValue ("Clean Table", (ObitData*)in->mosaic->images[i],
			     &ver, OBIT_IO_WriteOnly, noParms, err);
      if (err->error) Obit_traceback_msg (err, routine, in->name);
      in->CCver[i] = ver;  /* save if defaulted (0) */
      /* Delete old row */
      if (in->CCRow[i]) in->CCRow[i] = ObitTableCCUnref(in->CCRow[i]);
      
    }  /* End create table object */

    /* Reset and instantiate if needed */
    ObitTableClearRows ((ObitTable*)in->CCTable[i], err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
   
    /* Init number of CCs to 1*/
    ObitImageOpen(in->mosaic->images[i], OBIT_IO_ReadWrite, err);
    in->mosaic->images[i]->myDesc->niter = 1;
    in->mosaic->images[i]->myStatus = OBIT_Modified;
    ObitImageClose(in->mosaic->images[i], err);
  } /* end loop over images */

} /* end ObitDConCleanPxListMFReset */

/**
 * Resizes arras as needed
 * \param inn      The Pixel List object 
 * \param maxPixel Maximum number of pixels allowed (dim. or arrays)
 * \param err      Obit error stack object.
 */
void ObitDConCleanPxListMFResize (ObitDConCleanPxList *inn, olong maxPixel, 
				  ObitErr *err)
{
  ObitDConCleanPxListMF *in = (ObitDConCleanPxListMF*)inn;
  olong i;
  /*gchar *routine = " ObitDConCleanPxListMFResize";*/

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  if (maxPixel<in->maxPixel) return;  /* This needed? */

  in->maxPixel   = maxPixel;
  in->pixelX     = ObitMemRealloc (in->pixelX,    maxPixel*sizeof(olong));
  in->pixelY     = ObitMemRealloc (in->pixelY,    maxPixel*sizeof(olong));
  in->pixelFld   = ObitMemRealloc (in->pixelFld,  maxPixel*sizeof(gshort));
  in->pixelFlux  = ObitMemRealloc (in->pixelFlux, maxPixel*sizeof(ofloat));
  for (i=0; i<in->nSpec; i++) 
    in->channelFlux[i]  = ObitMemRealloc (in->channelFlux[i], maxPixel*sizeof(ofloat));

} /* end ObitDConCleanPxListMFResize */

/**
 * Update pixel list and window
 * If the frequency axis has ctype "SPECLOGF" and if "NTERM" exists in the first image 
 * descriptor InfoList and is > 0 and there are at least nterm planes in the image, 
 * then the subsequent planes are used as spectral terms and incorporated into the 
 * PixelList.
 * \param inn         The Pixel List object 
 * \param fields      Which fields? (1-rel) as zero terminated list
 * \param nSkip       Number of residuals to skip between ones accepted
 * \param minFluxLoad Minimum pixel flux density to accept
 * \param autoWinFlux min. residual flux allowed for auto Window
 * \param window      Windows object corresponding to Image Mosaic being CLEANED
 *                    Only pixels inside of the CLEAN window are used.
 * \param BeamPatch   Array of Dirty beam patch to use, for fields in mosaic
 * \param pixarray    If NonNULL use instead of the flux densities from the image file.
 *                    Array of ObitFArrays containing pixels for fields in fields
 *                    Only combined image.
 * \param err         Obit error stack object.
 */
void ObitDConCleanPxListMFUpdate (ObitDConCleanPxList *inn, 
				  olong *fields, olong nSkip,
				  ofloat minFluxLoad,
				  ofloat autoWinFlux,
				  ObitDConCleanWindow *window, 
				  ObitFArray **BeamPatch,
				  ObitFArray **pixarray,
				  ObitErr *err)
{
  ObitImageMF *image=NULL;
  ObitFArray **inFArrays;
  ofloat **sdata=NULL;
  olong i, j, field,  ifld, number, excess, ix, iy, nx, ny, pos[2], skipCnt;
  olong nplanes, iplane, naxis[2], nfield, nCCparms, parmoff;
  olong ip, lox, loy, hix, hiy,xoff, yoff;
  olong plane[5] = {1,1,1,1,1};
  ofloat maxChFlux, *data=NULL, fblank=ObitMagicF();
  gboolean blewIt=FALSE, *mask=NULL, rebuild=FALSE;
  ObitDConCleanPxListMF *in = (ObitDConCleanPxListMF*)inn;
  const ObitDConCleanPxListClassInfo *pxListClass;
  gchar *routine = "ObitDConCleanPxListMFUpdate";

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Save min flux density */
  in->minFluxLoad = minFluxLoad;
  in->autoWinFlux = autoWinFlux;  /* min. allowed for autoWindow */
  
  /* Diagnostics */
  if (in->prtLv>2) {
    Obit_log_error(err, OBIT_InfoErr, "%s field %d min %f auto min %f",
		   routine, fields[0], minFluxLoad, autoWinFlux);
  }

  /* Create inFArrays work array */
  nplanes   = 1 + in->nSpec;
  if (in->nSpec==1) nplanes = 1;  /* normal imaging */
  inFArrays = g_malloc0(nplanes*sizeof(ObitFArray*));
  sdata     = g_malloc0(nplanes*sizeof(ofloat*));

  /* Count number of fields */
  ifld = 0;
  field = fields[ifld];
  nfield = 0;
  while (field>0) {
    nfield++;
    ifld++;
    field = fields[ifld];
  } /* End counting fields */

  /* Create spectral beam array if necessary */
  if (((in->nSpec>1) && (in->nSBeamPatch<(nfield*in->nSpec))) || 
      (in->SBeamPatch==NULL)) {
    /* Get rid of old */
    if (in->SBeamPatch) {
      for (i=0; i<in->nSBeamPatch; i++) {
	in->SBeamPatch[i] = ObitFArrayUnref(in->SBeamPatch[i]);
      }
      in->SBeamPatch = ObitMemFree (in->SBeamPatch);
      in->sigma      = ObitMemFree (in->sigma);
    }
    /* New */
    in->nSBeamPatch = in->mosaic->numberImages * in->nSpec;
    in->SBeamPatch  = ObitMemAlloc0(in->nSBeamPatch*sizeof(ObitFArray*));
    in->sigma       = ObitMemAlloc0(in->nSBeamPatch*sizeof(ofloat));
  } /* end create spectral beam array */

  /* Loop over selected fields */
  ifld = 0;
  field = fields[ifld];
  number = excess = 0;
  skipCnt = 0;
  while (field>0) {
    
    /* Check field */
    if ((field<=0) || (field>in->nfield)) {
      Obit_log_error(err, OBIT_Error,"%s field %d out of range 1- %d in %s",
		     routine, field, in->nfield, in->mosaic->name);
      return;
    }
    
    /* Save Beam patch */
    in->BeamPatch[field-1] = ObitFArrayUnref(in->BeamPatch[field-1]);
    in->BeamPatch[field-1] = ObitFArrayRef(BeamPatch[field-1]);

    /* Which image? */
    image = (ObitImageMF*)in->mosaic->images[field-1];

    /* Get Spectral beam patches */
    if (in->nSpec>1)
      GetSpecBeamPatch (in, field-1, image, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);

    /* May need to rebuild CC Tables if their definitions are wrong 
       only increase record size */
    nfield = in->mosaic->numberImages;
    rebuild = FALSE;
    ip = 0;
    for (i=0; i<nfield; i++) {
      if (in->CCTable[i]->parmsCol>=0)
	nCCparms = in->CCTable[i]->myDesc->dim[in->CCTable[i]->parmsCol][0];
      else nCCparms = 0;
      if (in->circGaus[i]>0.0) parmoff = 4;
      else parmoff = 0;
      if (nCCparms<parmoff+in->nSpec) {
	rebuild = TRUE;
	in->CCTable[i] = ObitTableCCUnref(in->CCTable[i]);
	ObitImageZapTable (in->mosaic->images[i], "AIPS CC", in->CCver[i], err);
	if (err->error) Obit_traceback_msg (err, routine, in->name);
      }
    }
    if (rebuild) {
      Obit_log_error(err, OBIT_InfoWarn,"%s: Rebuilding CC Tables to add spectral terms",
		     routine);
      pxListClass = (ObitDConCleanPxListClassInfo*)in->ClassInfo; 
      pxListClass->ObitDConCleanPxListReset (inn, err);
      if (err->error) Obit_traceback_msg (err, routine, in->name);
    }
    
    /* Store combined flux in inFArrays[0] 
     Use pixarray[ifld] if given, else read */
    inFArrays[0]  = ObitFArrayUnref(inFArrays[0]);
    naxis[0] = image->myDesc->inaxes[0];
    naxis[1] = image->myDesc->inaxes[1];
    if (pixarray && pixarray[ifld]) {  /* Use pixarray */
      inFArrays[0] = ObitFArrayRef(pixarray[ifld]);
    } else {  /* Read from plane 1 */
      inFArrays[0] = ObitFArrayCreate (NULL, 2, naxis);
      plane[0] = 1;
      ObitImageGetPlane ((ObitImage*)image, inFArrays[0]->array, plane, err);
      if (err->error) Obit_traceback_msg (err, routine, image->name);
    }
    /* Read all coarse spectral planes */
    for (iplane=1; iplane<nplanes; iplane++) {
      inFArrays[iplane] = ObitFArrayUnref (inFArrays[iplane]);
      inFArrays[iplane] = ObitFArrayCreate (NULL, 2, naxis);
      /* Coarse spectrum follows spectral planes */
      plane[0] = 1 + in->maxOrder + iplane;
      ObitImageGetPlane ((ObitImage*)image, inFArrays[iplane]->array, plane, err);
      if (err->error) Obit_traceback_msg (err, routine, image->name);
      /* Get RMS */
      in->sigma[(field-1)*in->nSpec+iplane-1] = ObitFArrayRMS(inFArrays[iplane]);
    }
    
    /* pointers to data */
    nx = image->myDesc->inaxes[0];
    ny = image->myDesc->inaxes[1];
    lox = 0; hix = nx;
    loy = 0; hiy = ny;
    pos[0] = 0; pos[1] = loy;
    data = ObitFArrayIndex(inFArrays[0], pos);
    /* subtract closest integer to reference pixel */
    if (image->myDesc->crpix[0]>0.0)  
      xoff = (olong)(image->myDesc->crpix[0]+0.5);
    else xoff = (olong)(image->myDesc->crpix[0]-0.5);
    if (image->myDesc->crpix[1]>0.0)  
      yoff = (olong)(image->myDesc->crpix[1]+0.5);
    else yoff = (olong)(image->myDesc->crpix[1]-0.5);
    xoff--; yoff--; /* To 0 rel */

    /* Spectral planes */
    if (in->nSpec>1) {
      for (j=0; j<in->nSpec; j++) sdata[j] = ObitFArrayIndex(inFArrays[j+1], pos);
    }    

    /* Loop over image saving selected values */
    for (iy=loy; iy<hiy; iy++) {
      /* Get window mask */
      if (ObitDConCleanWindowRow(window, field, iy+1, &mask, err)) {
	for (ix=lox; ix<hix; ix++) {
	  /* Want this one? */
	  if (mask[ix] && (data[ix]!=fblank) && 
	      (fabs(data[ix])>=minFluxLoad)) {
	    if (skipCnt>=nSkip) { /* take this one */
	      skipCnt = 0;
	      if ((number+1)>=in->maxPixel) {
		blewIt = TRUE;   /* Blew arrays */
		excess++;
	      } else { /* OK */
		in->pixelX[number]    = ix - xoff;
		in->pixelY[number]    = iy - yoff;
		in->pixelFld[number]  = field;
		in->pixelFlux[number] = data[ix];
		/* Spectral planes - find max abs */
		maxChFlux = -1.0;
		if (in->nSpec>1) {
		  for (j=0; j<in->nSpec; j++) {
		    if (sdata[j][ix]!=fblank) {
		      in->channelFlux[j][number] = sdata[j][ix];
		      maxChFlux = MAX(fabs(sdata[j][ix]),maxChFlux);
		    } else in->channelFlux[j][number] = 0.0;
		  }
		}
		/* Only accept this one if the combined value is less than the max abs 
		   channel value.  This prevents the convolution from putting components
		   on top of a cell that is exxentially zero in all channels which can cause 
		   CLEAN to get stuck */
		if (fabs(data[ix])<maxChFlux) number++;
	      }
	    }
	  } else { /* skip this one to make them fit */
	    skipCnt++;
	  }
	}
      }
      /* Update pointers to next row */
      data += nx;
      for (j=0; j<in->nSpec; j++) sdata[j] += nx;
    } /* end loop in  y */

    in->nPixel = number;   /* Save actual number */

    /* Free Image data arrays */
    for (iplane=0; iplane<nplanes; iplane++) 
      inFArrays[iplane] = ObitFArrayUnref(inFArrays[iplane]);
    
    /* Cleanup - next field may have different size */
    if ((mask) && (ObitMemValid (mask))) mask = ObitMemFree (mask);

    ifld++;
    field = fields[ifld];
  } /* end loop over fields */

  /* Cleanup */
  if (inFArrays) g_free(inFArrays);
  if (sdata)     g_free(sdata);
  /* Free image buffer */
  image->image = ObitFArrayUnref(((ObitImage*)image)->image);

  /* Give warning if blew arrays */
  if (blewIt) 
    Obit_log_error(err, OBIT_InfoWarn,"%s: Number of pixels exceed arrays by  %d",
		   routine, excess);

} /* end ObitDConCleanPxListMFUpdate */

/**
 * Iteratively perform BGC CLEAN on Pixel lists
 * \param inn   The Pixel list object 
 * \param err   Obit error stack object.
 * \return TRUE if hit limit of niter or min. flux density.
 */
gboolean ObitDConCleanPxListMFCLEAN (ObitDConCleanPxList *inn, ObitErr *err)
{
  gboolean done = FALSE;
  olong iter, field=0,  beamPatch;
  olong i, j, lpatch, irow, xoff, yoff, lastField=-1;
  ofloat minFlux=0.0, lastFlux=0.0, CCmin, atlim, xfac=1.0, xflux;
  ofloat minFluxLoad;
  ofloat subval, peak, ccfLim=0.5;
  odouble totalFlux, *fieldFlux=NULL;
  gchar reason[51];
  ObitTableCCRow *CCRow = NULL;
  ObitImageDesc *desc = NULL;
  ObitIOCode retCode;
  olong ithread, maxThread, nThreads, nCCparms;
  CLEANFuncArg **targs=NULL;
  ObitThreadFunc myThreadFunc;
  olong npix, lopix, hipix, npixPerThread, parmoff;
  gboolean OK = TRUE, *doField=NULL, quit=FALSE;
  ObitDConCleanPxListMF *in = (ObitDConCleanPxListMF*)inn;
  /* Clean loop before writing */
#ifndef MAXBGCLOOP
#define MAXBGCLOOP 20
#endif
  olong ibgc, nbgc, mbgc, iXres, iYres, tipeak[MAXBGCLOOP];
  ofloat *tspec[MAXBGCLOOP], tcombPeak[MAXBGCLOOP];
  gchar *routine = "ObitDConCleanPxListMFCLEAN";
  /* DEBUG */
  olong pos[2]={0,0}, ix, ipeak=0;
  ofloat tmax;

  /* error checks */
  if (err->error) return done;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Check number of residuals - bail if none */
  if (in->nPixel<=0) {
    Obit_log_error(err, OBIT_InfoWarn,"%s NO Residuals to CLEAN in %s",
		   routine, in->name);
    in->complCode = OBIT_CompReasonNoPixel;
    return TRUE;
  }
  minFluxLoad = in->minFluxLoad;

  /* Local arrays */
  doField = g_malloc0(in->nfield*sizeof(gboolean));
  for (i=0; i<in->nfield; i++) doField[i] = FALSE;
  for (ibgc=0; ibgc<MAXBGCLOOP; ibgc++) {
    tipeak[ibgc] = 0;
    tcombPeak[ibgc] = 0.0;
    tspec[ibgc] = g_malloc0(in->nSpec*sizeof(ofloat));
  }

   /* How many components already done? */
  iter = MAX (0, in->currentIter);

  /* Set function */
  myThreadFunc = (ObitThreadFunc)ThreadCLEAN;

  /* Remove any blanks from beam patches */
  lpatch = 0;
  for (i=0; i<in->nfield; i++) {
    if (in->BeamPatch[i]) {
      ObitFArrayDeblank (in->BeamPatch[i], 0.0);
      /* max. beam patch */
      lpatch = MAX (lpatch, in->BeamPatch[i]->naxis[0]);
      /* Loop over spectral channel beams */
      if (in->nSpec>1) {
	for (j=0; j<in->nSpec; j++) 
	  if (in->SBeamPatch[i*in->nSpec+j])  
	      ObitFArrayDeblank (in->SBeamPatch[i*in->nSpec+j], 0.0);
      }
    }
  }

 /* Setup */
  beamPatch     = (lpatch-1)/2;
  CCmin         = 1.0e20;
  atlim         = 0.0;
  in->complCode = OBIT_CompReasonUnknown;  /* No reason for completion yet */

  /* Tell details */
  if (in->prtLv>1) {
    Obit_log_error(err, OBIT_InfoErr,"MF BGC CLEAN: Beam patch = %d cells, min. residual = %g Jy",
		   2*beamPatch, in->minFluxLoad);
    Obit_log_error(err, OBIT_InfoErr," %d residuals loaded ",
		   in->nPixel);
    ObitErrLog(err);  /* Progress Report */
  }

  /* Setup Threading */
  /* Only thread large cases */
  if (in->nPixel>500) maxThread = 1000;
  else maxThread = 1;
  nThreads = MakeCLEANArgs (in, maxThread, &targs, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, TRUE);

  /* Divide up work */
  npix = in->nPixel;
  npixPerThread = npix/nThreads;
  lopix = 1;
  hipix = npixPerThread;
  hipix = MIN (hipix, npix);

  /* Set up thread arguments */
  for (ithread=0; ithread<nThreads; ithread++) {
    if (ithread==(nThreads-1)) hipix = npix;  /* Make sure to do all */
    targs[ithread]->PixelList = in;
    targs[ithread]->first     = lopix;
    targs[ithread]->last      = hipix;
    if (ithread==nThreads) targs[ithread]->last = in->nPixel;
    if (nThreads>1) targs[ithread]->ithread = ithread;
    else targs[ithread]->ithread = -1;
    targs[ithread]->ipeak     = 0;
    targs[ithread]->combPeak  = 0.0;
    /* Update which pix */
    lopix += npixPerThread;
    hipix += npixPerThread;
    hipix = MIN (hipix, npix);
  }

  /* Local accumulators for flux */
  totalFlux = 0.0;
  fieldFlux = g_malloc0(in->nfield*sizeof(odouble));
  for (i=0; i<in->nfield; i++) fieldFlux[i] = 0.0;

  /* CLEAN loop */
  while (!done) {

    /* Set up for inner loop to enhance performance */
    nbgc = MIN (MAXBGCLOOP, in->niter-iter);
    for (ibgc=0; ibgc<nbgc; ibgc++) {   /* loop finding components */
    
      /* Do subtraction/ find next peak  */
      OK = ObitThreadIterator (in->thread, nThreads, 
			       (ObitThreadFunc)myThreadFunc, 
			       (gpointer**)targs);
      
      /* Check for problems */
      if (!OK) Obit_log_error(err, OBIT_Error,"%s: Problem in threading", routine);
      
      /* Get peak info */
      mbgc = 0;
      peak = fabs(targs[0]->combPeak);
      for (ithread=1; ithread<nThreads; ithread++) {
	if (fabs(targs[ithread]->combPeak)>peak) {
	  mbgc = ithread;
	  peak = fabs(targs[ithread]->combPeak);
	}
      }
      tcombPeak[ibgc] = targs[mbgc]->combPeak;
      tipeak[ibgc]    = targs[mbgc]->ipeak;
      for (j=0; j<in->nSpec; j++)
	tspec[ibgc][j] = targs[mbgc]->Spec[j] ;

      /* Set up thread arguments for next subtraction */
      for (ithread=0; ithread<nThreads; ithread++) {
	targs[ithread]->ipeak    = tipeak[ibgc];
	targs[ithread]->combPeak = tcombPeak[ibgc];
	for (j=0; j<in->nSpec; j++)
	  targs[ithread]->Spec[j]= tspec[ibgc][j];
      }
      
     } /* end bgc loop finding components */

    quit = FALSE;
    for (ibgc=0; ibgc<nbgc; ibgc++) {   /* loop saving/testing components */
      /* If first pass set up stopping criteria */
      xflux = tcombPeak[ibgc];
      field = in->pixelFld[tipeak[ibgc]] - 1;
      if (in->resMax < 0.0) {
	in->resMax = MAX (fabs(xflux), 1.0e-10);
	/* lower algorithm limit if fitted values are less that minFluxLoad */
	if (minFluxLoad>in->resMax) minFluxLoad *= 0.90;
      }      
      
      /* Save info */
      xfac     = pow ((minFluxLoad / in->resMax), in->factor[field]);
      ccfLim   = in->resMax*in->ccfLim;  /* Fraction of peak limit */
      lastFlux = in->resMax;
      doField[field] = TRUE;
      minFlux  = in->minFlux[field];
      subval   = tcombPeak[ibgc] * in->gain[field];
      lastFlux = tcombPeak[ibgc];
      iXres    = in->pixelX[tipeak[ibgc]];
      iYres    = in->pixelY[tipeak[ibgc]];
      CCmin    = MIN (CCmin, fabs(tcombPeak[ibgc]));
      
     /* Keep statistics */
      in->iterField[field]++;
      iter++;  /* iteration count */
      fieldFlux[field] += subval;
      totalFlux += subval;
      atlim += xfac / (ofloat)iter;   /* update BGC stopping criterion */
      
      /* Write component to Table */    
      /* Open table if not already open */
      if (in->CCTable[field]->myStatus == OBIT_Inactive) {
	retCode = ObitTableCCOpen (in->CCTable[field], OBIT_IO_ReadWrite, err);
	if ((retCode != OBIT_IO_OK) || (err->error))
	  Obit_traceback_val (err, routine, in->name, done);
      }
      /* Create row if needed */
      if (!in->CCRow[field]) in->CCRow[field] = newObitTableCCRow (in->CCTable[field]);
      CCRow = in->CCRow[field];   /* Get local pointer to  Table Row  */
      lastField = field;
      
      /* Set value */
      desc = in->mosaic->images[field]->myDesc;
      /* correct by offset in ObitDConCleanPxListUpdate */
      if (desc->crpix[0]>0.0)  
	xoff = (olong)(desc->crpix[0]+0.5);
      else xoff = (olong)(desc->crpix[0]-0.5);
      if (desc->crpix[1]>0.0)  
	yoff = (olong)(desc->crpix[1]+0.5);
      else yoff = (olong)(desc->crpix[1]-0.5);
      xoff--; yoff--; /* To 0 rel */
      /* What's in AIPS is a bit more complex and adds field offset from tangent */
      CCRow->DeltaX = (iXres - desc->crpix[0]+1+xoff)*desc->cdelt[0];
      CCRow->DeltaY = (iYres - desc->crpix[1]+1+yoff)*desc->cdelt[1];
      CCRow->Flux   =  subval;
      /* May need Gaussian components */
      parmoff = 0;
      if (in->CCTable[field]->parmsCol>=0)
	nCCparms = in->CCTable[field]->myDesc->dim[in->CCTable[field]->parmsCol][0];
      else nCCparms = 0;
      if (nCCparms>=4) {
	if (in->circGaus[field]>0.0) {
	  CCRow->parms[0] = in->circGaus[field];
	  CCRow->parms[1] = in->circGaus[field];
	  CCRow->parms[2] = 0.0;
	  CCRow->parms[3] = 1;  /* type 1 = Gaussian */
	  parmoff = 4;
	} else if (nCCparms>=(in->nSpecTerm+4)) { /* point */
	  CCRow->parms[0] = 0.0;
	  CCRow->parms[1] = 0.0;
	  CCRow->parms[2] = 0.0;
	  CCRow->parms[3] = 0;  /* type 0 = Point */
	  parmoff = 4;
	}
      } /* end add Gaussian components */
      
      /* May need Spectral components */
      if (nCCparms>=(parmoff+in->nSpec)) {
	
	if (in->nSpec>0) {
	  CCRow->parms[3] += 20.0;  /* mark as also having tabulated spectrum */
	  for (i=0; i<in->nSpec; i++) {
	    CCRow->parms[parmoff+i] = tspec[ibgc][i] * in->gain[field];
	  }
	} /* end add Spectral components */
      } /* end if need spectral component */
      /* DEBUG  */
      if (in->prtLv>5) 
	fprintf (stderr,"Component: fld %2d %5d comb %9.6f flux %9.6f tot %9.6f pos %6d  %6d \n",
		 field+1, tipeak[ibgc], xflux, subval, in->totalFlux+totalFlux, iXres, iYres);
      
      
      irow = in->iterField[field];
      retCode = ObitTableCCWriteRow (in->CCTable[field], irow, CCRow, err);
      if ((retCode != OBIT_IO_OK) || (err->error)) 
	Obit_traceback_val (err, routine, in->name, done);
      
      /* Test various stopping conditions */ 
      /* Are we finished after this? */
      done = done || (iter>=in->niter) || (fabs(xflux)<minFlux);
      if (done) {  /* set completion reason string */
	if (iter>=in->niter)     g_snprintf (reason, 50, "Reached Iter. limit");
	if (iter>=in->niter)     in->complCode = OBIT_CompReasonNiter; 
	if (fabs(xflux)<minFlux) g_snprintf (reason, 50, "Reached min Clean flux");
	if (fabs(xflux)<minFlux) in->complCode = OBIT_CompReasonMinFlux;
	quit = TRUE;
	break;
      } 
      
      /* Diverging? */
      if (fabs (xflux) > 2.0*CCmin) {
	/* Give warning */
	Obit_log_error(err, OBIT_InfoWarn,"%s: Clean has begun to diverge, Stopping",
		       routine);
	g_snprintf (reason, 50, "Solution diverging");
	in->complCode = OBIT_CompReasonDiverge;
	done = TRUE;
	quit = TRUE;
	break;
      }
      
      /* BGC tests to tell if we should quit now */
      if (fabs (xflux) < minFluxLoad * (1.0 + atlim)) {
	g_snprintf (reason, 50, "Reached minimum algorithm flux");
	in->complCode = OBIT_CompReasonBGCLimit;
	/* DEBUG  */
	if (in->prtLv>5)  {
	  tmax = -1.0e20;
	  for (ix=0; ix<in->nPixel; ix++) {
	    if (fabs(in->pixelFlux[ix]) > tmax) {
	      tmax = fabs(in->pixelFlux[ix]);
	      peak = in->pixelFlux[ix];
	      ipeak = ix;
	      pos[0] = in->pixelX[ix];
	      pos[1] = in->pixelY[ix];
	    }
	  }
	  fprintf (stderr, "Quit: %f < %f, minFluxLoad %f real peak %f @ %d %d %d\n",  
		   xflux, minFluxLoad * (1.0 + atlim), minFluxLoad, peak, pos[0], pos[1], ipeak);
	}  /* End DEBUG */
	quit = TRUE;
	break;  /* jump out of CLEAN loop */
      }
      
      /* autoWindow tests to tell if we should quit now */
      if (fabs (xflux) < in->autoWinFlux) {
	g_snprintf (reason, 50, "Reached minimum autoWindow flux");
	in->complCode = OBIT_CompReasonAutoWin;
	/* DEBUG - this to keep SubNewCCs from running until it's fixed 
	   in->complCode = OBIT_CompReasonMinFract;*/
	quit = TRUE;
	break;  /* jump out of CLEAN loop */
      }
      
      /* Deep enough fraction of peak residual */
      if (fabs(xflux)<ccfLim) {
	g_snprintf (reason, 50, "Reached min fract of peak resid");
	in->complCode = OBIT_CompReasonMinFract;
	quit = TRUE;
	break;
      }
    } /* End loop saving, testing components */
    if (done || quit) break;
  } /* end CLEANing loop */
  
  
  /* Keep statistics */
  in->currentIter = iter;
  /* Save accumulators for flux */
  in->totalFlux += totalFlux;
  for (i=0; i<in->nfield; i++) in->fluxField[i] += fieldFlux[i];
  if (fieldFlux) g_free(fieldFlux);

  /* Loop over CC tables closing */
  for (i=0; i<in->nfield; i++) {
    if (in->CCTable[i]->myStatus != OBIT_Inactive) {
      retCode = ObitTableCCClose (in->CCTable[i], err);
      if ((retCode != OBIT_IO_OK) || (err->error))
	Obit_traceback_val (err, routine, in->name, done);
      in->CCRow[i] = ObitTableCCRowUnref(in->CCRow[i]);
    }
  } /* end loop closing tables */
  
    /* Cleanup */
  KillCLEANArgs (nThreads, targs);
  ObitThreadPoolFree (in->thread);

  /* Clear spectral beam patches */
  if (in->SBeamPatch) {
    for (i=0; i<in->nSBeamPatch; i++) {
      in->SBeamPatch[i] = ObitFArrayUnref(in->SBeamPatch[i]);
    }
  }
  
  /* Tell about results */
  if (in->prtLv>1) {
    Obit_log_error(err, OBIT_InfoErr,"Clean stopped because: %s", reason);
    Obit_log_error(err, OBIT_InfoErr,"Min. Flux density %f",
		   xflux);
    if (in->nfield>1) /* Multiple fields? */
      for (i=0; i<in->nfield; i++) {
	if (doField[i]) {
	  Obit_log_error(err, OBIT_InfoErr,"Field %d has %d CCs with %g Jy",
			 i+1, in->iterField[i], in->fluxField[i]);
	}
      }
    
    Obit_log_error(err, OBIT_InfoErr,"Total CLEAN %d CCs with %g Jy",
		   in->currentIter, in->totalFlux);
  }
  /* Keep maximum abs residual */
  in->maxResid = fabs(xflux);

  /* Cleanup */
  if (doField) g_free(doField);
  for (ibgc=0; ibgc<MAXBGCLOOP; ibgc++) {
    if (tspec[ibgc]) g_free(tspec[ibgc]);
  }

  return done;
} /* end ObitDConCleanPxListMFCLEAN */

/**
 * Steer-Dewney-Ito-Greisen CLEAN on Pixel list
 * Lifted from AIPS
 * \param inn   The Pixel list object 
 * \param err   Obit error stack object.
 * \return TRUE if hit limit of niter or min. flux density.
 */
gboolean ObitDConCleanPxListMFSDI (ObitDConCleanPxList *inn, ObitErr *err)
{
  ObitDConCleanPxListMF *in = (ObitDConCleanPxListMF*)inn;
  gboolean drop, done = FALSE;
  olong iter, iresid, field=0, beamPatch, lpatch=0, ipeak, iXres, iYres;
  olong lastField=-1, xoff, yoff, irow, i, j;
  ofloat minFlux=0.0, xflux;
  ofloat minVal=-1.0e20, *sum=NULL, *wt=NULL, *tspec=NULL, mapLim;
  odouble totalFlux, *fieldFlux=NULL;
  gchar reason[51];
  ObitTableCCRow *CCRow = NULL;
  ObitImageDesc *desc = NULL;
  ObitIOCode retCode;
  olong ithread, maxThread, nThreads, nCCparms, ispec, nSpec;
  CLEANFuncArg **targs=NULL;
  olong npix, lopix, hipix, npixPerThread, parmoff;
  gboolean OK = TRUE, *doField=NULL;
  gchar *routine = "ObitDConCleanPxListMFSDI";

  /* error checks */
  if (err->error) return done;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Check number of residuals - bail if none */
  if (in->nPixel<=0) {
    Obit_log_error(err, OBIT_InfoWarn,"%s NO Residuals to CLEAN in %s",
		   routine, in->name);
    in->complCode = OBIT_CompReasonNoPixel;
    return TRUE;
  }

  /* Number of spectral channels */
  nSpec = in->nSpec;
  if ((in->nSpec==1) || (in->curOrder==0)) nSpec = 1;

  /* Local arrays */
  sum     = g_malloc0(nSpec*sizeof(ofloat));
  wt      = g_malloc0(nSpec*sizeof(ofloat));
  tspec   = g_malloc0(nSpec*sizeof(ofloat));
  doField = g_malloc0(in->nfield*sizeof(gboolean));
  for (i=0; i<in->nfield; i++) doField[i] = FALSE;

   /* How many components already done? */
  iter = MAX (0, in->currentIter);

  /* Remove any blanks from beam patches */
  lpatch = 0;
  for (i=0; i<in->nfield; i++) {
    if (in->BeamPatch[i]) {
      ObitFArrayDeblank (in->BeamPatch[i], 0.0);
      /* max. beam patch */
      lpatch = MAX (lpatch, in->BeamPatch[i]->naxis[0]);
      /* Loop over spectral channel beams */
      if (in->nSpec>1) {
	for (j=0; j<in->nSpec; j++) 
	  if (in->SBeamPatch[i*in->nSpec+j])  
	      ObitFArrayDeblank (in->SBeamPatch[i*in->nSpec+j], 0.0);
      }
    }
  }

  /* Setup */
  beamPatch     = (lpatch-1)/2;
  mapLim        = in->minFluxLoad;
  in->complCode = OBIT_CompReasonUnknown;  /* No reason for completion yet */

  /* Adjust data array - amount in excess of mapLim */
  for (iresid=0; iresid<in->nPixel; iresid++) {
    if (in->pixelFlux[iresid]>mapLim) in->pixelFlux[iresid] -= mapLim;
    else if (in->pixelFlux[iresid]<-mapLim) in->pixelFlux[iresid] += mapLim;
    else in->pixelFlux[iresid] = 0.0;
    for (j=0; j<in->nSpec; j++) {
      if (in->channelFlux[j][iresid]>mapLim)  in->channelFlux[j][iresid] -= mapLim;
      else if (in->channelFlux[j][iresid]<-mapLim) in->channelFlux[j][iresid] += mapLim;
      else in->channelFlux[j][iresid] = 0.0;
    }
  } /* end loop over array */

  /* Tell details */
  if (in->prtLv>1) {
    Obit_log_error(err, OBIT_InfoErr,"SDI CLEAN Beam patch = %d cells, CLEAN above = %g Jy",
		   2*beamPatch, in->minFluxLoad);
    Obit_log_error(err, OBIT_InfoErr," %d residuals loaded ",
		   in->nPixel);
    ObitErrLog(err);  /* Progress Report */
  }


  /* Setup Threading */
  /* Only thread large cases */
  if (in->nPixel>500) maxThread = 1000;
  else maxThread = 1;
  nThreads = MakeCLEANArgs (in, maxThread, &targs, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, TRUE);

  /* Divide up work */
  npix = in->nPixel;
  npixPerThread = npix/nThreads;
  lopix = 1;
  hipix = npixPerThread;
  hipix = MIN (hipix, npix);

  /* Set up thread arguments */
  for (ithread=0; ithread<nThreads; ithread++) {
    if (ithread==(nThreads-1)) hipix = npix;  /* Make sure to do all */
    targs[ithread]->PixelList = in;
    targs[ithread]->first     = lopix;
    targs[ithread]->last      = hipix;
    if (nThreads>1) targs[ithread]->ithread = ithread;
    else targs[ithread]->ithread = -1;
    targs[ithread]->ipeak     = 0;
    for (ispec=0; ispec<nSpec; ispec++) targs[ithread]->sum[ispec] = 0.0;
    /* Update which pix */
    lopix += npixPerThread;
    hipix += npixPerThread;
    hipix = MIN (hipix, npix);
  }

  /* Local accumulators for flux */
  totalFlux = 0.0;
  fieldFlux = g_malloc0(in->nfield*sizeof(odouble));
  for (i=0; i<in->nfield; i++) fieldFlux[i] = 0.0;
    
  /* Do outer loop over list */
  for (ipeak=0; ipeak<in->nPixel; ipeak++) {
    xflux   = in->pixelFlux[ipeak];
    for (ispec=0; ispec<nSpec; ispec++) tspec[ispec] = in->channelFlux[ispec][ipeak];
    if (fabs(xflux)<=0.0) continue;
    iXres   = in->pixelX[ipeak];
    iYres   = in->pixelY[ipeak];
    field   = in->pixelFld[ipeak];
    minFlux = in->minFlux[field-1];
    doField[field-1] = TRUE;
    
    /* Determine weight factor for this pixel - 
       dot product of beam and data array */
    /* Set up thread arguments */
    for (ithread=0; ithread<nThreads; ithread++) {
      targs[ithread]->ipeak = ipeak;
      for (ispec=0; ispec<nSpec; ispec++) targs[ithread]->sum[ispec] = 0.0;
    }
    /* operation possibly in threads */
    OK = ObitThreadIterator (in->thread, nThreads, 
			     (ObitThreadFunc)ThreadSDICLEAN, 
			     (gpointer**)targs);

    /* Check for problems */
    if (!OK) Obit_log_error(err, OBIT_Error,"%s: Problem in threading", routine);
    
    /* sum per spectral channel */
    for (ispec=0; ispec<nSpec; ispec++) sum[ispec] = targs[0]->sum[ispec];
    for (ithread=1; ithread<nThreads; ithread++) {
      for (ispec=0; ispec<nSpec; ispec++) sum[ispec] += targs[ithread]->sum[ispec];
    }
    
    /* Weights for this pixel */
    drop = TRUE;
    for (ispec=0; ispec<nSpec; ispec++) {
      if (sum[ispec]!=0.0) wt[ispec] = fabs(tspec[ispec]/sum[ispec]);
      else wt[ispec] = 1.0;
      wt[ispec] = MAX (0.001, MIN (0.5, wt[ispec]));
      tspec[ispec] *= wt[ispec];
      drop = drop && fabs(tspec[ispec])<=0.0;
    }
    minVal = MAX (minVal, fabs(xflux)+mapLim);
    xflux *= wt[0];
    if (drop) continue;
   
    /* Keep statistics */
    in->iterField[field-1]++;
    fieldFlux[field-1] += xflux;
    totalFlux          += xflux;
    
    /* Write component to Table */
    /* Open table if not already open */
    if (in->CCTable[field-1]->myStatus == OBIT_Inactive) {
      retCode = ObitTableCCOpen (in->CCTable[field-1], OBIT_IO_ReadWrite, err);
      if ((retCode != OBIT_IO_OK) || (err->error))
	Obit_traceback_val (err, routine, in->name, done);
    }
    /* Create row if needed */
    if (!in->CCRow[field-1]) in->CCRow[field-1] = newObitTableCCRow (in->CCTable[field-1]);
    CCRow = in->CCRow[field-1];  /* Get local pointer to Table Row  */
    lastField = field;
    
    /* Set value */
    desc = in->mosaic->images[field-1]->myDesc;
    /* correct by offset in ObitDConCleanPxListUpdate */
    if (desc->crpix[0]>0.0)  
      xoff = (olong)(desc->crpix[0]+0.5);
    else xoff = (olong)(desc->crpix[0]-0.5);
    if (desc->crpix[1]>0.0)  
      yoff = (olong)(desc->crpix[1]+0.5);
    else yoff = (olong)(desc->crpix[1]-0.5);
    xoff--; yoff--; /* To 0 rel */
    /* What's in AIPS is a bit more complex and adds field offset from tangent */
    CCRow->DeltaX = (iXres - desc->crpix[0]+1+xoff)*desc->cdelt[0];
    CCRow->DeltaY = (iYres - desc->crpix[1]+1+yoff)*desc->cdelt[1];
    CCRow->Flux   = xflux;
    /* May need Gaussian components */
    parmoff = 0;
    if (in->CCTable[field-1]->parmsCol>=0)
      nCCparms = in->CCTable[field-1]->myDesc->dim[in->CCTable[field-1]->parmsCol][0];
    else nCCparms = 0;
    if (nCCparms>=4) {
      if (in->circGaus[field-1]>0.0) {
	CCRow->parms[0] = in->circGaus[field-1];
	CCRow->parms[1] = in->circGaus[field-1];
	CCRow->parms[2] = 0.0;
	CCRow->parms[3] = 1;  /* type 1 = Gaussian */
	parmoff = 4;
      } else if (nCCparms>=(in->nSpecTerm+4)) { /* point */
	CCRow->parms[0] = 0.0;
	CCRow->parms[1] = 0.0;
	CCRow->parms[2] = 0.0;
	CCRow->parms[3] = 0;  /* type 0 = Point */
	parmoff = 4;
      }
    } /* end add Gaussian components */

    /* May need Spectral components */
    if (nCCparms>=(parmoff+in->nSpec)) {

      if (in->nSpec>0) {
	CCRow->parms[3] += 20.0;  /* mark as also having tabulated spectrum */
	for (ispec=0; ispec<nSpec; ispec++) {
	  CCRow->parms[parmoff+ispec] = tspec[ispec];
	}
      } /* end add Spectral components */
    } /* end if need spectral component */

    /* DEBUG  */
    if (in->prtLv>5) 
      fprintf (stderr,"Component: fld %2d %5d comb %9.6f tot %9.6f pos %6d  %6d S: %f %f %f %f %f %f %f %f %f %f\n",
	       field, ipeak, xflux, in->totalFlux+totalFlux, iXres, iYres,
	       tspec[0], tspec[1], tspec[2], tspec[3], tspec[4],
	       tspec[5], tspec[6], tspec[7], tspec[8], tspec[9]);
    

    irow = in->iterField[field-1];
    retCode = ObitTableCCWriteRow (in->CCTable[field-1], irow, CCRow, err);
    if ((retCode != OBIT_IO_OK) || (err->error)) 
      Obit_traceback_val (err, routine, in->name, done);
    
    /* Test various stopping conditions */ 
    iter++;  /* iteration count */
    if (iter>=in->niter) {
      in->complCode = OBIT_CompReasonNiter; 
      break;
    }
     
  } /* end loop over list */
  
  /* Keep statistics */
  in->currentIter = iter;

  /* Save accumulators for flux */
  in->totalFlux += totalFlux;
  for (i=0; i<in->nfield; i++) in->fluxField[i] += fieldFlux[i];
  if (fieldFlux) g_free(fieldFlux);

  /* Are we finished after this? */
  done = (iter>=in->niter) || (fabs(minVal)<=minFlux);
  g_snprintf (reason, 50, "Completed layer");
  in->complCode = OBIT_CompReasonSDILayer;
  if (iter>=in->niter) {
    g_snprintf (reason, 50, "Reached Iter. limit");
    in->complCode = OBIT_CompReasonNiter; 
  }
  if (fabs(minVal)<=minFlux) {
    g_snprintf (reason, 50, "Reached min. flux density");
    in->complCode = OBIT_CompReasonMinFlux;
  }
  
  /* Loop over CC tables closing */
  for (i=0; i<in->nfield; i++) {
    if (in->CCTable[i]->myStatus != OBIT_Inactive) {
      retCode = ObitTableCCClose (in->CCTable[i], err);
      if ((retCode != OBIT_IO_OK) || (err->error))
	Obit_traceback_val (err, routine, in->name, done);
      in->CCRow[i] = ObitTableCCRowUnref(in->CCRow[i]);
    }
  } /* end loop closing tables */
  
  /* Cleanup */
  KillCLEANArgs (nThreads, targs);
  ObitThreadPoolFree (in->thread);  
 
  /* Tell about results */
  if (in->prtLv>1) {
    Obit_log_error(err, OBIT_InfoErr,"Clean stopped because: %s", reason);
    Obit_log_error(err, OBIT_InfoErr,"%s: Min. Flux density %f",
		   routine, minVal);
    if (in->nfield>1) /* Multiple fields? */
      for (i=0; i<in->nfield; i++) {
	if (doField[i]) {
	  Obit_log_error(err, OBIT_InfoErr,"Field %d has %d CCs with %g Jy",
			 i+1, in->iterField[i], in->fluxField[i]);
	}
      }
    
    Obit_log_error(err, OBIT_InfoErr,"Total CLEAN %d CCs with %g Jy",
		   in->currentIter, in->totalFlux);
  }
  /* Keep maximum abs residual */
  in->maxResid = fabs(in->minFluxLoad);

  /* Cleanup */
  if (doField) g_free(doField);
  if (sum)     g_free(sum);
  if (wt)      g_free(wt);
  if (tspec)   g_free(tspec);

  return done;
} /* end ObitDConCleanPxListMFSDI */

/**
 * Tells results of CLEANing
 * \param inn         The Pixel List object 
 * \param ncomps      Array of total number of components per field, 
 *                    same order as in ObitImageMosaic
 * \param err Obit error stack object.
 */
olong ObitDConCleanPxListMFResult (ObitDConCleanPxList *inn, olong *ncomp,
				   ObitErr *err)
{
  olong out = 0;
  olong i;
  ObitDConCleanPxListMF *in = (ObitDConCleanPxListMF*)inn;

  /* error checks */
  if (err->error) return out;
  g_assert (ObitIsA(in, &myClassInfo));

  out = in->currentIter;
  for (i=0; i<in->nfield; i++) ncomp[i] = in->iterField[i];

  return out;
} /* end ObitDConCleanPxListMFResult */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitDConCleanPxListMFClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitDConCleanPxListMFClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitDConCleanPxListMFClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitDConCleanPxListMFClassInfoDefFn (gpointer inClass)
{
  ObitDConCleanPxListMFClassInfo *theClass = (ObitDConCleanPxListMFClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitDConCleanPxListMFClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitDConCleanPxListMFClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitDConCleanPxListMFGetClass;
  theClass->newObit       = (newObitFP)newObitDConCleanPxListMF;
  theClass->ObitCopy      = (ObitCopyFP)ObitDConCleanPxListMFCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitDConCleanPxListMFClear;
  theClass->ObitInit      = (ObitInitFP)ObitDConCleanPxListMFInit;
  theClass->ObitDConCleanPxListCreate   = 
    (ObitDConCleanPxListCreateFP)ObitDConCleanPxListMFCreate;
  theClass->ObitDConCleanPxListGetParms = 
    (ObitDConCleanPxListGetParmsFP)ObitDConCleanPxListMFGetParms;
  theClass->ObitDConCleanPxListReset    = 
    (ObitDConCleanPxListResetFP)ObitDConCleanPxListMFReset;
  theClass->ObitDConCleanPxListResize   = 
    (ObitDConCleanPxListResizeFP)ObitDConCleanPxListMFResize;
  theClass->ObitDConCleanPxListUpdate   = 
    (ObitDConCleanPxListUpdateFP)ObitDConCleanPxListMFUpdate;
  theClass->ObitDConCleanPxListCLEAN    = 
    (ObitDConCleanPxListCLEANFP)ObitDConCleanPxListMFCLEAN;
  theClass->ObitDConCleanPxListSDI      = 
    (ObitDConCleanPxListSDIFP)ObitDConCleanPxListMFSDI;
  theClass->ObitDConCleanPxListResult   = 
    (ObitDConCleanPxListResultFP)ObitDConCleanPxListMFResult;

} /* end ObitDConCleanPxListMFClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitDConCleanPxListMFInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitDConCleanPxListMF *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->channelFlux = NULL;
  in->specFreq    = NULL;
  in->sigma       = NULL;
  in->SBeamPatch  = NULL;
  in->nSBeamPatch = 0;
  in->maxOrder    = 0;
  in->curOrder    = 0;
  in->nSpec       = 0;
  in->refFreq     = 1.0;

} /* end ObitDConCleanPxListMFInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitDConCleanPxListMF* cast to an Obit*.
 */
void ObitDConCleanPxListMFClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  olong i;
  ObitDConCleanPxListMF *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  if (in->specFreq) in->specFreq = ObitMemFree (in->specFreq);
  if (in->channelFlux) {
    for (i=0; i<in->nSpec; i++) {
      if (in->channelFlux[i]) 
	in->channelFlux[i] = ObitMemFree (in->channelFlux[i]) ; 
    }
    in->channelFlux = ObitMemFree (in->channelFlux);
  }
  
  if (in->SBeamPatch) {
    for (i=0; i<in->nSBeamPatch; i++) {
      in->SBeamPatch[i] = ObitFArrayUnref(in->SBeamPatch[i]);
    }
    in->SBeamPatch = ObitMemFree (in->SBeamPatch);
  }
  in->sigma = ObitMemFree(in->sigma);
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitDConCleanPxListMFClear */

/**
 * Read spectral beam patches for the spectral planes in the beam for image
 * These are places in in->SBeamPatch elements 
 * ifld*in->nSpec to (ifld+1)*in->nSpec-1
 * \param in    The Pixel list object 
 * \param ifld  0-rel field number
 * \param image to extract beam patches
 * \param err   Obit error stack object.
 */
static void GetSpecBeamPatch (ObitDConCleanPxListMF *in, olong ifld, 
			      ObitImageMF* image, ObitErr *err)
{
  ObitImage *theBeam;
  olong ablc[2], atrc[2], pos[2], plane[]={1,1,1,1,1};
  olong i, ip, nx, ny, icenx, iceny, beamPatchSize;
  ofloat fmax, zero=0.0;
  ObitImageClassInfo *imgClass;
  ObitFArray *FAtemp=NULL;
  gboolean bad;
  gchar *routine = "GetSpecBeamPatch";

  /* error checks */
  if (err->error) return;
  Obit_return_if_fail((in->nSBeamPatch>=((ifld+1)*in->nSpec)), err,
 		      "%s: SBeamPatch NOT big enough for %s", routine, in->name);

  imgClass  = (ObitImageClassInfo*)image->ClassInfo;
  theBeam   = imgClass->ObitImageGetBeam((ObitImage*)image, 0, plane, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Get beam patch size same as for combined beam */
  beamPatchSize = in->BeamPatch[ifld]->naxis[0]/2;

  /* Center pixel  */
  nx = theBeam->myDesc->inaxes[0];
  ny = theBeam->myDesc->inaxes[1];
  
  /* Loop over planes */
  ip = ifld*in->nSpec;  /* Initial output element */
  for (i=0; i<in->nSpec; i++) {
    in->SBeamPatch[ip] = ObitFArrayUnref(in->SBeamPatch[ip]);

    /* Read whole plane */
    plane[0] = i+in->maxOrder+2;
    ObitImageGetPlane (theBeam, NULL, plane, err);
    if (err->error) Obit_traceback_msg (err, routine, image->name);
    
    /* Setup to find peak in inner quarter - corners may be high */
    icenx = nx/2;
    iceny = ny/2;
    ablc[0] = icenx - nx/4;
    atrc[0] = icenx + nx/4;
    ablc[1] = iceny - nx/4;
    atrc[1] = iceny + nx/4;

    /* Find center in inner quarter */
    FAtemp = ObitFArraySubArr(theBeam->image, ablc, atrc, err);
    if (err->error) Obit_traceback_msg (err, routine, theBeam->name);
    fmax = ObitFArrayMax (FAtemp, pos);
    FAtemp = ObitFArrayUnref(FAtemp);

    /* Set if Beam OK - peak>0.5 and near center */
    bad = (fmax<0.5) || (fmax>1.01) || 
      (abs(pos[0]+ablc[0]-icenx)>3) || (abs(pos[1]+ablc[1]-iceny)>3);
    if (!bad) {
      icenx = pos[0]+ablc[0];
      iceny = pos[1]+ablc[1];
    }
    
    /* Beam patch window as 0-rel */
    ablc[0] = icenx - beamPatchSize;
    atrc[0] = icenx + beamPatchSize;
    ablc[1] = iceny - beamPatchSize;
    atrc[1] = iceny + beamPatchSize;
    
    /* Save Beam patch */
    in->SBeamPatch[ip] = ObitFArraySubArr(theBeam->image, ablc, atrc, err);
    if (err->error) Obit_traceback_msg (err, routine, theBeam->name);

    /* If beam bad, zero */
    if (bad) ObitFArrayFill (in->SBeamPatch[ip], zero);
    
    ip++;  /* next */
  } /* end loop over planes */
 
  /* Free Image array? */
  theBeam->image = ObitFArrayUnref(theBeam->image);
  
} /* end GetSpecBeamPatch  */

/**
 * Inner loop BGC CLEAN with spectral estimation
 * Subtracts component from residuals using beam patch and
 * returns largest residual in region cleaned
 * Callable as thread
 * Arguments are given in the structure passed as arg
 * \param arg Pointer to CLEANFuncArg argument with elements:
 * \li PixelList PixelList object
 * \li first    First (1-rel) pixel no. to process this thread
 * \li last     Highest (1-rel) pixel no. to process this thread
 * \li ithread  thread number, <0-> no threads
 * \li ipeak    [in/out] Pixel number of peak
 * \li combPeak [in/out] peak flux to subtract 
 * \li Spec     [in/out] Array of nSpec flux densities for component 
 * \return NULL
 */
gpointer ThreadCLEAN (gpointer args)
{
  CLEANFuncArg *largs  = (CLEANFuncArg*)args;
  ObitDConCleanPxListMF *in = largs->PixelList;
  olong  loPix            = largs->first-1;
  olong  hiPix            = largs->last;
  olong  ipeak            = largs->ipeak;
  ofloat combPeak         = largs->combPeak;
  ofloat *Spec            = largs->Spec;

  ofloat xflux, subCval, peak;
  olong j, iresid, nSpec, ispec, iXres, iYres, lpatch, beamPatch, iBeam, field, pos[2];
  ofloat *beam00=NULL, *Sbeam=NULL;

  /* Number of spectral channels */
  nSpec = in->nSpec;
  if ((in->nSpec==1) || (in->curOrder==0)) nSpec = 0;

  /* Do subtraction if peak non zero */
  field  = in->pixelFld[ipeak];
  lpatch = in->BeamPatch[field-1]->naxis[0];
  beamPatch = lpatch/2;
  peak  = -1.0;
  xflux  = 0.0;
  if (fabs(combPeak)>0.0) {
    pos[0]  = pos[1] = 0;
    beam00  = ObitFArrayIndex(in->BeamPatch[field-1], pos); /* Beam patch pointer */
    subCval = combPeak * in->gain[field-1];  /* Combined flux */
    iXres   = in->pixelX[ipeak];
    iYres   = in->pixelY[ipeak];
    for (iresid=loPix; iresid<hiPix; iresid++) {
      /* Is this inside the Beam patch ? */
      if ((abs(in->pixelX[iresid]-iXres) <= beamPatch) && 
	  (abs(in->pixelY[iresid]-iYres) <= beamPatch)) {
	/* Index in beam patch array */
	iBeam = (beamPatch + (in->pixelY[iresid] - iYres)) * lpatch +
	  (beamPatch + (in->pixelX[iresid] - iXres));
	in->pixelFlux[iresid] -= subCval * beam00[iBeam];
	/* Subtract from spectral planes */
	for (ispec=0; ispec<nSpec; ispec++) {
	  Sbeam = in->SBeamPatch[(field-1)*in->nSpec+ispec]->array; /* Beam patch pointer */
	  in->channelFlux[ispec][iresid] -= in->gain[field-1] * Spec[ispec] * Sbeam[iBeam];
	}
      }
      /* Now, Find Peak */
      xflux = fabs(in->pixelFlux[iresid]);
      if (xflux>peak) {
	peak = xflux;
	ipeak = iresid;
      }
    } /* end loop over array */
  } else {  /* Only find peak */
    
    /* Find peak abs residual */
    for (iresid=loPix; iresid<hiPix; iresid++) {
      xflux = fabs(in->pixelFlux[iresid]);
      if (xflux>peak) {
	peak = xflux;
	ipeak = iresid;
      }
    }
  }

  /* Deal with normal imaging */
  if (nSpec<=0) {
    largs->combPeak = in->pixelFlux[ipeak];
    largs->ipeak    = ipeak;
    goto finish;
  }

  /* Return values */
  largs->combPeak = in->pixelFlux[ipeak];
  for (j=0; j<nSpec; j++) largs->Spec[j] = in->channelFlux[j][ipeak];
  largs->ipeak  = ipeak;

  /* Indicate completion */
 finish:
  if (largs->ithread>=0)
    ObitThreadPoolDone (in->thread, (gpointer)&largs->ithread);
  return NULL;
} /* end ThreadCLEAN */

/**
 * Inner loop SDI CLEAN
 * Gets dot product of beam with residuals
 * Callable as thread
 * Arguments are given in the structure passed as arg
 * \param arg Pointer to CLEANFuncArg argument with elements:
 * \li PixelList PixelList object
 * \li first   First (1-rel) pixel no. to process this thread
 * \li last    Highest (1-rel) pixel no. to process this thread
 * \li ithread thread number, <0-> no threads
 * \li sum;    [out] dot product returned
 * \return NULL
 */
gpointer ThreadSDICLEAN (gpointer args)
{
  CLEANFuncArg *largs  = (CLEANFuncArg*)args;
  ObitDConCleanPxListMF *in = largs->PixelList;
  olong  loPix            = largs->first-1;
  olong  hiPix            = largs->last;
  olong  ipeak            = largs->ipeak;
  ofloat *sum             = largs->sum;

  ofloat *beam=NULL;
  olong iresid, nSpec, ispec, iXres, iYres, lpatch, beamPatch, iBeam, field, pos[2];

  /* Number of spectral channels */
  nSpec = in->nSpec;
  if ((in->nSpec==1) || (in->curOrder==0)) nSpec = 0;

  /* Setup */
  for (ispec=0; ispec<nSpec; ispec++) sum[ispec] = 0.0;
  iXres  = in->pixelX[ipeak];
  iYres  = in->pixelY[ipeak];
  field  = in->pixelFld[ipeak];
  lpatch = in->BeamPatch[field-1]->naxis[0];
  beamPatch = (lpatch-1)/2;

  /* Determine weight factor for this pixel - 
     dot product of beam and data array*/
  /* Loop over spectral channels */
  for (ispec = 0; ispec<nSpec; ispec++) {
    pos[0]  = pos[1] = 0;
    beam    = ObitFArrayIndex( in->SBeamPatch[(field-1)*in->nSpec+ispec], pos); 
    for (iresid=loPix; iresid<hiPix; iresid++) {
      /* Is this inside the Beam patch ? */
      if ((abs(in->pixelY[iresid]-iYres) <= beamPatch) && 
	  (abs(in->pixelX[iresid]-iXres) <= beamPatch)) {
	/* Index in beam patch array */
	iBeam = (beamPatch + (in->pixelY[iresid] - iYres)) * lpatch +
	  (beamPatch + (in->pixelX[iresid] - iXres));
	sum[ispec] += in->channelFlux[ispec][iresid] * beam[iBeam];
      }
    } /* end loop over array */
  } /* end loop over channel */

  /* Indicate completion */
  if (largs->ithread>=0)
    ObitThreadPoolDone (in->thread, (gpointer)&largs->ithread);
  return NULL;
} /* end ThreadSDICLEAN */

/**
 * Make arguments for Threaded CLEAN
 * \param in         CLEAN object
 * \param maxThread  Maximum desirable no. threads
 * \param args       Created array of CLEANFuncArg, 
 *                   delete with KillCLEANArgs
 * \return number of elements in args.
 */
static olong MakeCLEANArgs (ObitDConCleanPxListMF *in, olong maxThread,
			    CLEANFuncArg ***args, ObitErr *err)
{
  olong i, nThreads;
  gchar *routine = "MakeCLEANArgs";
  
  /* Setup for threading */
  /* How many threads? */
  nThreads = MAX (1, ObitThreadNumProc(in->thread));
  nThreads = MIN (nThreads, maxThread);

  /* Initialize threadArg array */
  *args = g_malloc0(nThreads*sizeof(CLEANFuncArg*));
  for (i=0; i<nThreads; i++) 
    (*args)[i] = g_malloc0(sizeof(CLEANFuncArg)); 
  
  for (i=0; i<nThreads; i++) {
    (*args)[i]->PixelList = in;
    (*args)[i]->first     = 1;
    (*args)[i]->last      = 0;
    (*args)[i]->ithread   = i;
    (*args)[i]->ipeak     = 0;
    if (in->nSpec>1) (*args)[i]->Spec = g_malloc0(in->nSpec*sizeof(ofloat));
    else (*args)[i]->Spec = NULL;
    if (in->nSpec>1) (*args)[i]->sum = g_malloc0(in->nSpec*sizeof(ofloat));
    else (*args)[i]->sum = g_malloc0(sizeof(ofloat));
    if (err->error) Obit_traceback_val (err, routine, in->name, nThreads);
  }

  return nThreads;
} /*  end MakeCLEANArgs */

/**
 * Delete arguments for Threaded CLEAN
 * \param nargs      number of elements in args.
 * \param args       Array of CLEANFuncArg, type CLEANFuncArg
 */
static void KillCLEANArgs (olong nargs, CLEANFuncArg **args)
{
  olong i;

  if (args==NULL) return;
  for (i=0; i<nargs; i++) {
    if (args[i]) {
      if (args[i]->Spec) g_free(args[i]->Spec);
      if (args[i]->sum)  g_free(args[i]->sum);
      g_free(args[i]);
    }
  }
  g_free(args);
} /*  end KillCLEANArgs */


