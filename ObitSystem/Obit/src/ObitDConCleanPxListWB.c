/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2010                                               */
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

#include "ObitDConCleanPxListWB.h"
#include "ObitMem.h"
#include "ObitImageWB.h"
#include "ObitImageMosaicWB.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitDConCleanPxListWB.c
 * ObitDConCleanPxListWB class function definitions.
 * This class determines the pixel histogram of an image.
 * This class is derived from the Obit base class.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitDConCleanPxListWB";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitDConCleanPxListGetClass;

/**
 * ClassInfo structure ObitDConCleanPxListWBClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitDConCleanPxListWBClassInfo myClassInfo = {FALSE};

/*--------------- File Global Variables  ----------------*/

/*---------------Private structures----------------*/
/* CLEAN threaded function argument */
typedef struct {
  /* PixelList object */
  ObitDConCleanPxListWB    *PixelList;
  /* First (1-rel) record in otfdata buffer to process this thread */
  olong        first;
  /* Highest (1-rel) record in otfdata buffer to process this thread  */
  olong        last;
  /* thread number , >0 -> no threading  */
  olong        ithread;
  /* Pixel number of peak */
  olong       ipeak;
  /* Peak flux to subtract */
  ofloat       peak;
  /* Spectral index of peak */
  ofloat       si;
  /* Spectral curvature of peak */
  ofloat       curve;
  /* 2nd order spectral curvature of peak */
  ofloat       curve2;
} CLEANFuncArg;

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitDConCleanPxListWBInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitDConCleanPxListWBClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitDConCleanPxListWBClassInfoDefFn (gpointer inClass);

/** Private: Threaded CLEAN  for different beam orders */
static gpointer ThreadCLEAN0 (gpointer arg);
static gpointer ThreadCLEAN1 (gpointer arg);
static gpointer ThreadCLEAN2 (gpointer arg);

/** Private: Make arguments for Threaded CLEAN */
static olong MakeCLEANArgs (ObitDConCleanPxListWB *in, olong maxThread,
			    CLEANFuncArg ***args);

/** Private: Delete arguments for Threaded CLEAN */
static void KillCLEANArgs (olong nargs, CLEANFuncArg **args);

/** Private: Estimate exponential alpha from linear spectral index */
static void Linear2Exp1 (ObitDConCleanPxListWB *in, ofloat a, ofloat *alpha);

/** Private: Estimate exponential alpha,beta from linear spectral index, curve */
static void Linear2Exp2 (ObitDConCleanPxListWB *in, ofloat a, ofloat b,
			 ofloat *alpha, ofloat *beta);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitDConCleanPxListWB* newObitDConCleanPxListWB (gchar* name)
{
  ObitDConCleanPxListWB* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitDConCleanPxListWBClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitDConCleanPxListWB));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitDConCleanPxListWBInit((gpointer)out);

 return out;
} /* end newObitDConCleanPxListWB */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitDConCleanPxListWBGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitDConCleanPxListWBClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitDConCleanPxListWBGetClass */

/**
 * Make a deep copy of an ObitDConCleanPxListWB.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitDConCleanPxListWB* 
ObitDConCleanPxListWBCopy  (ObitDConCleanPxListWB *in, ObitDConCleanPxListWB *out, 
			  ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  gboolean oldExist;
  gchar *outName;
  /*gchar *routine = "ObitDConCleanPxListWBCopy";*/

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
    out = newObitDConCleanPxListWB(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  return out;
} /* end ObitDConCleanPxListWBCopy */

/**
 * Make a copy of a object but do not copy the actual data
 * This is useful to create an DConCleanPxListWB similar to the input one.
 * \param in  The object to copy
 * \param out An existing object pointer for output, must be defined.
 * \param err Obit error stack object.
 */
void ObitDConCleanPxListWBClone  (ObitDConCleanPxListWB *in, ObitDConCleanPxListWB *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  /*gchar *routine = "ObitDConCleanPxListWBClone";*/

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
} /* end ObitDConCleanPxListWBClone */

/**
 * Creates an ObitDConCleanPxListWB
 * Get spectral order from first image on mosaic
 * \param name     An optional name for the object.
 * \param mosaic   Image mosaic to be deconvolved.
 * \param uvdata   Data from which images are being made
 * \param maxPixel Maximum number of pixels allowed (dim. or arrays)
 * \return the new object.
 */
ObitDConCleanPxListWB* 
ObitDConCleanPxListWBCreate (gchar* name, ObitImageMosaic *mosaic,  
			     ObitUV *uvdata, olong maxPixel, ObitErr *err)
{
  ObitDConCleanPxListWB* out=NULL;
  olong i, nfield, num, nif, nfreq, maxOrder;
  gchar *routine = "ObitDConCleanPxListWBCreate";

  /* error checks */
  if (err->error) return out;
  g_assert (ObitImageMosaicIsA(mosaic));
  /* Check input types */
  Obit_retval_if_fail((ObitImageMosaicWBIsA((ObitImageMosaicWB*)mosaic) && 
		       (ObitImageWBIsA((ObitImageWB*)mosaic->images[0]))), err, 
		      out,
		      "%s: Image mosaic or images not WB", routine);
  
  /* Create basic structure */
  out = newObitDConCleanPxListWB (name);
  
  /* How many spectral orders? - get maximum from mosaic */
  /* There may be different orders in different fields - deal with it */
  maxOrder = 0;
  for (i=0; i<mosaic->numberImages; i++)
    maxOrder = MAX(maxOrder, ObitImageWBGetBeamOrder(mosaic->images[i]));
  out->nSpecTerm = maxOrder;
  out->curOrder  = maxOrder;
  /* Create arrays including parent class */
  out->maxPixel   = maxPixel;
  out->pixelX     = ObitMemAlloc0Name (maxPixel*sizeof(olong),  "PxListWB X pixel");
  out->pixelY     = ObitMemAlloc0Name (maxPixel*sizeof(olong),  "PxListWB Y pixel");
  out->pixelFld   = ObitMemAlloc0Name (maxPixel*sizeof(gshort), "PxListWB pixel field");
  out->pixelFlux  = ObitMemAlloc0Name (maxPixel*sizeof(ofloat), "PxListWB pixel Flux");
  if (out->nSpecTerm>0) { /* first order */
    out->pixelFlux1 = ObitMemAlloc0Name (maxPixel*sizeof(ofloat), "PxListWB pixel Flux1");
  }
  if (out->nSpecTerm>1) { /* second order */
    out->pixelFlux2 = ObitMemAlloc0Name (maxPixel*sizeof(ofloat), "PxListWB pixel Flux2");
  }
  if (out->nSpecTerm>2) { /* third order */
    out->pixelFlux3 = ObitMemAlloc0Name (maxPixel*sizeof(ofloat), "PxListWB pixel Flux3");
  }

  /* Save Image Mosaic reference */
  out->mosaic = ObitImageMosaicRef(mosaic);

  /* Per field */
  nfield      = mosaic->numberImages;
  out->nfield = nfield;
  out->iterField  = ObitMemAlloc0Name (nfield*sizeof(olong),  "PxListWB Clean CC count");
  out->CCver      = ObitMemAlloc0Name (nfield*sizeof(olong),  "PxListWB Clean CC version");
  out->fluxField  = ObitMemAlloc0Name (nfield*sizeof(ofloat), "PxListWB Clean Flux");
  out->circGaus   = ObitMemAlloc0Name (nfield*sizeof(ofloat), "PxListWB Gaussian");
  out->gain       = ObitMemAlloc0Name (nfield*sizeof(ofloat), "PxListWB Clean Loop gain");
  out->minFlux    = ObitMemAlloc0Name (nfield*sizeof(ofloat), "PxListWB Clean Mix flux");
  out->factor     = ObitMemAlloc0Name (nfield*sizeof(ofloat), "PxListWB Clean factor");
  out->CCTable    = ObitMemAlloc0Name (nfield*sizeof(ObitTableCC*), "PxListWB CC tables");
  /* Beam patches - which depend or spectral order */
  out->BeamPatch00  = ObitMemAlloc0Name (nfield*sizeof(ObitFArray*), "Beam00");
  if (out->nSpecTerm>0) { /* first order */
    out->BeamPatch01  = ObitMemAlloc0Name (nfield*sizeof(ObitFArray*), "Beam01");
    out->BeamPatch10  = ObitMemAlloc0Name (nfield*sizeof(ObitFArray*), "Beam10");
    out->BeamPatch11  = ObitMemAlloc0Name (nfield*sizeof(ObitFArray*), "Beam11");
  }
  if (out->nSpecTerm>1) { /* second order */
    out->BeamPatch02  = ObitMemAlloc0Name (nfield*sizeof(ObitFArray*), "Beam02");
    out->BeamPatch20  = ObitMemAlloc0Name (nfield*sizeof(ObitFArray*), "Beam20");
    out->BeamPatch21  = ObitMemAlloc0Name (nfield*sizeof(ObitFArray*), "Beam21");
    out->BeamPatch12  = ObitMemAlloc0Name (nfield*sizeof(ObitFArray*), "Beam12");
    out->BeamPatch22  = ObitMemAlloc0Name (nfield*sizeof(ObitFArray*), "Beam22");
  }
		 
  for (i=0; i<nfield; i++) {
    out->iterField[i] = 0;
    out->CCver[i]     = 0;
    out->fluxField[i] = 0.0;
    out->circGaus[i]  = 0.0;
    out->gain[i]      = 0.1;
    out->minFlux[i]   = 0.0;
    out->factor[i]    = 0.0;
    out->CCTable[i]   = NULL;
  }
  
  /* Frequency info */
  out->refFreq = ((ObitImageWB*)mosaic->images[0])->refFreq; /* Reference from image */
  /* Find max and min frequencies in uv data */
  out->loFreq  =  1.0e25;
  out->hiFreq  = -1.0e25;
  if (uvdata->myDesc->jlocif>=0) nif = uvdata->myDesc->inaxes[uvdata->myDesc->jlocif];
  else nif = 1;
  nfreq = uvdata->myDesc->inaxes[uvdata->myDesc->jlocf];
  num = nif * nfreq;
  for (i=0; i<num; i++) {
    out->loFreq = MIN (out->loFreq, uvdata->myDesc->freqArr[i]);
    out->hiFreq = MAX (out->hiFreq, uvdata->myDesc->freqArr[i]);
  }

  return out;
} /* end ObitDConCleanPxListWBCreate */

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
void  ObitDConCleanPxListWBGetParms (ObitDConCleanPxList *inn, ObitErr *err)
{
  ObitDConCleanPxListClassInfo *ParentClass;
  ObitDConCleanPxListWB *in = (ObitDConCleanPxListWB*)inn;
  /*union ObitInfoListEquiv InfoReal; */
  gchar *routine = "ObitDConCleanPxListWBGetParms";

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Read any parent class parameters */
  ParentClass = (ObitDConCleanPxListClassInfo*)myClassInfo.ParentClass;
  ParentClass->ObitDConCleanPxListGetParms(inn, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

} /* end ObitDConCleanPxListWBGetParms */

/**
 * Resets CLEAN information
 * Makes sure all potential CC tables are instantiated
 * \param inn         The Pixel List object 
 * \param err Obit error stack object.
 */
void ObitDConCleanPxListWBReset (ObitDConCleanPxList *inn, ObitErr *err)
{
  olong i, ver=0;
  oint noParms;
  ObitDConCleanPxListWBClassInfo *myClass;
  ObitDConCleanPxListWB *in = (ObitDConCleanPxListWB*)inn;
  gchar *routine = " ObitDConCleanPxListWBReset";

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Control info */
  myClass = (ObitDConCleanPxListWBClassInfo*)in->ClassInfo;
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
      noParms = in->nSpecTerm;     /* Possible spectra */
      /* If adding spectra, also need type parameters even for point */
      if ((in->circGaus[i]>0.0) || (in->nSpecTerm>0)) noParms += 4;
      else noParms = 0;
      ver = MAX (ver, in->CCver[i]);  /* Use last if not defined */
      in->CCTable[i] = ObitTableCCUnref(in->CCTable[i]);  /* Free old */
      in->CCTable[i] = 
	newObitTableCCValue ("Clean Table", (ObitData*)in->mosaic->images[i],
			     &ver, OBIT_IO_WriteOnly, noParms, err);
      if (err->error) Obit_traceback_msg (err, routine, in->name);
      in->CCver[i] = ver;  /* save if defaulted (0) */
      
    }  /* End create table object */

    /* Reset and instantiate if needed */
    ObitTableClearRows ((ObitTable*)in->CCTable[i], err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
   
  } /* end loop over images */

} /* end ObitDConCleanPxListWBReset */

/**
 * Resets CLEAN information
 * Makes sure all potential CC tables are instantiated
 * \param inn      The Pixel List object 
 * \param maxPixel Maximum number of pixels allowed (dim. or arrays)
 * \param err      Obit error stack object.
 */
void ObitDConCleanPxListWBResize (ObitDConCleanPxList *inn, olong maxPixel, 
				  ObitErr *err)
{
  ObitDConCleanPxListWB *in = (ObitDConCleanPxListWB*)inn;
  /*gchar *routine = " ObitDConCleanPxListWBResize";*/

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  if (maxPixel<in->maxPixel) return;  /* This needed? */

  in->maxPixel   = maxPixel;
  in->pixelX     = ObitMemRealloc (in->pixelX,    maxPixel*sizeof(olong));
  in->pixelY     = ObitMemRealloc (in->pixelY,    maxPixel*sizeof(olong));
  in->pixelFld   = ObitMemRealloc (in->pixelFld,  maxPixel*sizeof(gshort));
  in->pixelFlux  = ObitMemRealloc (in->pixelFlux, maxPixel*sizeof(ofloat));
  if (in->nSpecTerm>0) { /* first order */
    in->pixelFlux1 = ObitMemRealloc (in->pixelFlux1, maxPixel*sizeof(ofloat));
  }
  if (in->nSpecTerm>1) { /* second order */
    in->pixelFlux2 = ObitMemRealloc (in->pixelFlux2, maxPixel*sizeof(ofloat));
  }
  if (in->nSpecTerm>2) { /* third order */
    in->pixelFlux3 = ObitMemRealloc (in->pixelFlux3, maxPixel*sizeof(ofloat));
  }

} /* end ObitDConCleanPxListWBResize */

/**
 * Update pixel list and window
 * Only consider facet with the highest order number.
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
 *                    Not used for ObitDConCleanPxListWB
 * \param pixarray    If NonNULL use instead of the flux densities from the image file.
 *                    Array of ObitFArrays containing pixels for fields in fields
 *                    1 per spectral order
 * \param err Obit error stack object.
 */
void ObitDConCleanPxListWBUpdate (ObitDConCleanPxList *inn, 
				  olong *fields, olong nSkip,
				  ofloat minFluxLoad,
				  ofloat autoWinFlux,
				  ObitDConCleanWindow *window, 
				  ObitFArray **BeamPatch,
				  ObitFArray **pixarray,
				  ObitErr *err)
{
  ObitImageWB *image=NULL;
  ObitFArray *inFArrays[10]={NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,NULL,};
  olong i, field,  ifld, number, excess, ix, iy, nx, ny, pos[2], skipCnt;
  olong nplanes, iplane, naxis[2], nfield, nCCparms, parmoff;
  olong ip, lox, loy, hix, hiy, order, maxOrder, xoff, yoff;
  ofloat *data=NULL, *data1=NULL, *data2=NULL, *data3=NULL, fblank=ObitMagicF();
  gboolean blewIt=FALSE, *mask=NULL, rebuild=FALSE;
  ObitDConCleanPxListWB *in = (ObitDConCleanPxListWB*)inn;
  const ObitDConCleanPxListClassInfo *pxListClass;
  gchar *routine = "ObitDConCleanPxListWBUpdate";

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
  
  /* There may be different orders in different fields - deal with it */
  maxOrder = 0;
  for (ifld=0; ifld<in->nfield; ifld++) {
    field = fields[ifld];
    if (field <=0) break;
    maxOrder = MAX(maxOrder, ObitImageWBGetBeamOrder(in->mosaic->images[field-1]));
  }
  in->curOrder = maxOrder;  /* current imaging order */

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
    
    /* Which image? */
    image = (ObitImageWB*)in->mosaic->images[field-1];

    /* Order of this image */
    order = ObitImageWBGetBeamOrder(in->mosaic->images[field-1]);
    if (order<maxOrder) continue;  /* bother with this one? */

    /* Need work arrays? */
    if ((image->ResidArr[0]==NULL) || (image->BeamMatx[0][0]==NULL))
      ObitImageWBMakeWork(image, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    
    /* Save Beam patches - 0 order */
    in->BeamPatch00[field-1] = ObitFArrayUnref(in->BeamPatch00[field-1]);
    in->BeamPatch00[field-1] = ObitFArrayRef(image->BeamMatx[0][0]);
    if ((in->nSpecTerm>0) && (order>=1)) { /* first order */
      in->BeamPatch01[field-1] = ObitFArrayUnref(in->BeamPatch01[field-1]);
      in->BeamPatch01[field-1] = ObitFArrayRef(image->BeamMatx[0][1]);
      in->BeamPatch10[field-1] = ObitFArrayUnref(in->BeamPatch10[field-1]);
      in->BeamPatch10[field-1] = ObitFArrayRef(image->BeamMatx[1][0]);
      in->BeamPatch11[field-1] = ObitFArrayUnref(in->BeamPatch11[field-1]);
      in->BeamPatch11[field-1] = ObitFArrayRef(image->BeamMatx[1][1]);
    }
    if ((in->nSpecTerm>1) && (order>=2)) { /* second order */
      in->BeamPatch02[field-1] = ObitFArrayUnref(in->BeamPatch02[field-1]);
      in->BeamPatch02[field-1] = ObitFArrayRef(image->BeamMatx[0][2]);
      in->BeamPatch20[field-1] = ObitFArrayUnref(in->BeamPatch20[field-1]);
      in->BeamPatch20[field-1] = ObitFArrayRef(image->BeamMatx[2][0]);
      in->BeamPatch21[field-1] = ObitFArrayUnref(in->BeamPatch21[field-1]);
      in->BeamPatch21[field-1] = ObitFArrayRef(image->BeamMatx[2][1]);
      in->BeamPatch12[field-1] = ObitFArrayUnref(in->BeamPatch12[field-1]);
      in->BeamPatch12[field-1] = ObitFArrayRef(image->BeamMatx[1][2]);
      in->BeamPatch22[field-1] = ObitFArrayUnref(in->BeamPatch22[field-1]);
      in->BeamPatch22[field-1] = ObitFArrayRef(image->BeamMatx[2][2]);
    }
    
    /* May need to rebuild CC Tables since their definitions are wrong 
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
      if (nCCparms<parmoff+in->nSpecTerm) {
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
    
    /* Loop over planes of flux, spectra (if given) and save in inFArray[*] */
    nplanes = 1 + order;
    naxis[0] = image->myDesc->inaxes[0];
    naxis[1] = image->myDesc->inaxes[1];
    for (iplane=0; iplane<nplanes; iplane++) {
      /* If pixarray given use it instead of the flux plane */
      inFArrays[iplane]  = ObitFArrayUnref(inFArrays[iplane]);
      if (pixarray && pixarray[ip]) {
	inFArrays[iplane]  = ObitFArrayRef(pixarray[ip++]);
      } else {
	inFArrays[iplane]  = ObitFArrayRef(((ObitImageWB*)image)->ResidArr[iplane]);
      }
      
    } /* end loop fetching planes */
    
    /* pointers to data */
    nx = image->myDesc->inaxes[0];
    ny = image->myDesc->inaxes[1];
    lox = nx/4; hix = 3*nx/4;
    loy = ny/4; hiy = 3*ny/4;
    pos[0] = 0; pos[1] = loy;
    data = ObitFArrayIndex(inFArrays[0], pos);
    /* subtract closest interger to reference pixel */
    if (image->myDesc->crpix[0]>0.0)  
      xoff = (olong)(image->myDesc->crpix[0]+0.5);
    else xoff = (olong)(image->myDesc->crpix[0]-0.5);
    if (image->myDesc->crpix[1]>0.0)  
      yoff = (olong)(image->myDesc->crpix[1]+0.5);
    else yoff = (olong)(image->myDesc->crpix[1]-0.5);
    /* Spectral planes */
    if (order>=1)  data1 = ObitFArrayIndex(inFArrays[1], pos);
    if (order>=2)  data2 = ObitFArrayIndex(inFArrays[2], pos);
    if (order>=3)  data3 = ObitFArrayIndex(inFArrays[3], pos);
    
    /* Loop over image saving selected values - only select those in inner quarter.*/
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
		if (order>=1) in->pixelFlux1[number] = data1[ix];
		if (order>=2) in->pixelFlux2[number] = data2[ix];
		if (order>=3) in->pixelFlux3[number] = data3[ix];
		number++;
	      }
	    }
	  } else { /* skip this one to make them fit */
	    skipCnt++;
	  }
	}
      }
      /* Update pointers to next row */
      data += nx;
      if (order>=1) data1 += nx;
      if (order>=2) data2 += nx;
      if (order>=3) data3 += nx;
    } /* end loop in  y */

    in->nPixel = number;   /* Save actual number */

    /* Free Image data arrays */
    for (iplane=0; iplane<nplanes; iplane++) 
      inFArrays[iplane]  = ObitFArrayUnref(inFArrays[iplane]);
    
    /* Cleanup - next field may have different size */
    if ((mask) && (ObitMemValid (mask))) mask = ObitMemFree (mask);

    ifld++;
    field = fields[ifld];
  } /* end loop over fields */

  /* Give warning if blew arrays */
  if (blewIt) 
    Obit_log_error(err, OBIT_InfoWarn,"%s: Number of pixels exceed arrays by  %d",
		   routine, excess);

} /* end ObitDConCleanPxListWBUpdate */

/**
 * Iteratively perform BGC CLEAN al la Sault-Wieringa on Pixel lists
 * \param inn   The Pixel list object 
 * \param err   Obit error stack object.
 * \return TRUE if hit limit of niter or min. flux density.
 */
gboolean ObitDConCleanPxListWBCLEAN (ObitDConCleanPxList *inn, ObitErr *err)
{
  gboolean done = FALSE;
  olong iter, field=0, iXres, iYres, beamPatch, tipeak=0;
  olong i, lpatch, irow, xoff, yoff, lastField=-1;
  ofloat minFlux=0.0, lastFlux=0.0, CCmin, atlim, xfac=1.0, xflux;
  ofloat peak, tpeak, tsi, tcurve, tcurve2;
  ofloat subval, ccfLim=0.5;
  ofloat a, b, alpha=-0.7, beta=0.0;
  odouble totalFlux, *fieldFlux=NULL;
  gchar reason[51];
  ObitTableCCRow *CCRow = NULL;
  ObitImageDesc *desc = NULL;
  ObitIOCode retCode;
  olong ithread, maxThread, nThreads, nCCparms;
  CLEANFuncArg **targs=NULL;
  ObitThreadFunc myThreadFunc;
  olong npix, lopix, hipix, npixPerThread, parmoff;
  gboolean OK = TRUE, *doField=NULL;
  ObitDConCleanPxListWB *in = (ObitDConCleanPxListWB*)inn;
  gchar *routine = "ObitDConCleanPxListWBCLEAN";

  /* error checks */
  if (err->error) return done;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Check number of residuals - bail if none */
  if (in->nPixel<=0) {
    Obit_log_error(err, OBIT_InfoWarn,"%s NO Residuals to CLEAN in %s",
		   routine, in->name);
    return TRUE;
  }

  /* Local arrays */
  doField = g_malloc0(in->nfield*sizeof(gboolean));
  for (i=0; i<in->nfield; i++) doField[i] = FALSE;

   /* How many components already done? */
  iter = MAX (0, in->currentIter);

  /* Set function  by spectral order */
  myThreadFunc = (ObitThreadFunc)ThreadCLEAN0;
  if (in->curOrder==1) myThreadFunc = (ObitThreadFunc)ThreadCLEAN1;
  if (in->curOrder==2) myThreadFunc = (ObitThreadFunc)ThreadCLEAN2;

   /* Remove any blanks from beam patches */
  lpatch = 0;
  for (i=0; i<in->nfield; i++) {
    if (in->BeamPatch00[i]) {
      ObitFArrayDeblank (in->BeamPatch00[i], 0.0);
      /* max. beam patch */
      lpatch = MAX (lpatch, in->BeamPatch00[i]->naxis[0]);
    }
  }

 /* Setup */
  beamPatch     = (lpatch-1)/2;
  CCmin         = 1.0e20;
  atlim         = 0.0;
  in->complCode = OBIT_CompReasonUnknown;  /* No reason for completion yet */

  /* Tell details */
  if (in->prtLv>1) {
    Obit_log_error(err, OBIT_InfoErr,"SW BGC CLEAN: Beam patch = %d cells, min. residual = %g Jy",
		   2*beamPatch, in->minFluxLoad);
    Obit_log_error(err, OBIT_InfoErr," %d residuals loaded ",
		   in->nPixel);
    ObitErrLog(err);  /* Progress Report */
  }

  /* Setup Threading */
  /* Only thread large cases */
  if (in->nPixel>500) maxThread = 1000;
  else maxThread = 1;
  /* Threading doesn't seem to help much
  maxThread = 1; */
  nThreads = MakeCLEANArgs (in, maxThread, &targs);

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
    targs[ithread]->peak      = 0.0;
    targs[ithread]->si        = 0.0;
    targs[ithread]->curve     = 0.0;
    targs[ithread]->curve2    = 0.0;
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
    
    /* Do subtraction/ find next peak  */
    OK = ObitThreadIterator (in->thread, nThreads, 
			     (ObitThreadFunc)myThreadFunc, 
			     (gpointer**)targs);

    /* Check for problems */
    if (!OK) Obit_log_error(err, OBIT_Error,"%s: Problem in threading", routine);
    
    /* Get  next peak info */
    peak  = fabs(targs[0]->peak);
    tpeak  = targs[0]->peak;
    tsi    = targs[0]->si;
    tcurve = targs[0]->curve;
    tcurve2= targs[0]->curve2;
    tipeak = targs[0]->ipeak;
    for (ithread=1; ithread<nThreads; ithread++) {
      if (fabs(targs[ithread]->peak)>peak) {
	peak   = fabs(targs[ithread]->peak);
	tpeak  = targs[ithread]->peak;
	tsi    = targs[ithread]->si;
	tcurve = targs[ithread]->curve;
	tcurve2= targs[ithread]->curve2;
	tipeak = targs[ithread]->ipeak;
      }
    }

    /* If first pass set up stopping criteria */
    xflux = tpeak;
    field   = in->pixelFld[tipeak];
    if (in->resMax < 0.0) {
      in->resMax = MAX (fabs(xflux), 1.0e-10);
    }      
    
    /* Save info */
    xfac = pow ((in->minFluxLoad / in->resMax), in->factor[field-1]);
    ccfLim = in->resMax*in->ccfLim;  /* Fraction of peak limit */
    lastFlux = in->resMax;
    doField[field-1] = TRUE;
    minFlux = in->minFlux[field-1];
    /* Damp any divergence */
    if (tpeak<-fabs(lastFlux)) subval = -fabs(lastFlux) * in->gain[field-1];
    else if (tpeak>fabs(lastFlux)) subval = -fabs(lastFlux) * in->gain[field-1];
    else subval = tpeak * in->gain[field-1];
    subval = tpeak * in->gain[field-1];/* DEBUG - no goes unstable */
    lastFlux = tpeak;
    iXres   = in->pixelX[tipeak];
    iYres   = in->pixelY[tipeak];
    CCmin   = MIN (CCmin, fabs(tpeak));

    /* Set up thread arguments for next subtraction */
    for (ithread=0; ithread<nThreads; ithread++) {
      targs[ithread]->ipeak  = tipeak;
      /*targs[ithread]->peak   = subval/in->gain[field-1]; */
      targs[ithread]->peak   = tpeak;
      targs[ithread]->si     = tsi;
      targs[ithread]->curve  = tcurve;
      targs[ithread]->curve2 = tcurve2;
    }
    
    /* Keep statistics */
    in->iterField[field-1]++;
    iter++;  /* iteration count */
    fieldFlux[field-1] += subval;
    totalFlux += subval;
    atlim += xfac / (ofloat)iter;   /* update BGC stopping criterion */

    /* Write component to Table */    
    /* Open table if not already open */
    if (in->CCTable[field-1]->myStatus == OBIT_Inactive) {
      retCode = ObitTableCCOpen (in->CCTable[field-1], OBIT_IO_ReadWrite, err);
      if ((retCode != OBIT_IO_OK) || (err->error))
	Obit_traceback_val (err, routine, in->name, done);
    }
    
    /* Need Table Row - if different field then it may be different */
    if (field!=lastField) CCRow = ObitTableCCRowUnref(CCRow);
    lastField = field;
    if (!CCRow) CCRow = newObitTableCCRow (in->CCTable[field-1]);
    
    /* Set value */
    desc = in->mosaic->images[field-1]->myDesc;
    /* correct by offset in ObitDConCleanPxListUpdate */
    if (desc->crpix[0]>0.0)  
      xoff = (olong)(desc->crpix[0]+0.5);
    else xoff = (olong)(desc->crpix[0]-0.5);
    if (desc->crpix[1]>0.0)  
      yoff = (olong)(desc->crpix[1]+0.5);
    else yoff = (olong)(desc->crpix[1]-0.5);
    /* What's in AIPS is a bit more complex and adds field offset from tangent */
    CCRow->DeltaX = (iXres - desc->crpix[0]+1+xoff)*desc->cdelt[0];
    CCRow->DeltaY = (iYres - desc->crpix[1]+1+yoff)*desc->cdelt[1];
    CCRow->Flux   =  subval;
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
    if (nCCparms>=(parmoff+in->nSpecTerm)) {

      /* Convert linear to exponsntial coefficients */
      if (in->nSpecTerm==1) {
	beta = 0.0;
	a = tsi/xflux;
	Linear2Exp1 (in, a, &alpha);
      } else if (in->nSpecTerm==2) {
	a = tsi/xflux;
	b = tcurve/xflux;
	Linear2Exp2 (in, a, b, &alpha, &beta);
      }
      if (in->nSpecTerm>0) {
	CCRow->parms[3] += 10.0;  /* mark as also having spectrum */
	if (in->nSpecTerm>0) CCRow->parms[0+parmoff] = alpha;
	if (in->nSpecTerm>1) CCRow->parms[1+parmoff] = beta;
	if (in->nSpecTerm>2) CCRow->parms[2+parmoff] = tcurve2/xflux;
      } /* end add Spectral components */
    } /* end if need spectral component */
    /* DEBUG  */
    fprintf (stderr,"Component: field %d flux %f total %f pos %d  %d si %6.3f curve %6.3f\n",
	     field, subval, in->totalFlux+totalFlux, iXres+1, iYres+1, alpha, beta);
    

    irow = in->iterField[field-1];
    retCode = ObitTableCCWriteRow (in->CCTable[field-1], irow, CCRow, err);
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
      break;
    } 

    /* Diverging? */
    if (fabs (xflux) > 2.0*CCmin) {
      /* Give warning */
      Obit_log_error(err, OBIT_InfoWarn,"%s: Clean has begun to diverge, Stopping",
		     routine);
      g_snprintf (reason, 50, "Solution diverging");
      in->complCode = OBIT_CompReasonDiverge;
      break;
    }
  
    /* BGC tests to tell if we should quit now */
    if (fabs (xflux) < in->minFluxLoad * (1.0 + atlim)) {
      g_snprintf (reason, 50, "Reached minimum algorithm flux");
      in->complCode = OBIT_CompReasonBGCLimit;
      break;  /* jump out of CLEAN loop */
    }
   
    /* autoWindow tests to tell if we should quit now */
    if (fabs (xflux) < in->autoWinFlux) {
      g_snprintf (reason, 50, "Reached minimum autoWindow flux");
      in->complCode = OBIT_CompReasonAutoWin;
      /* DEBUG - this to keep SubNewCCs from running until it's fixed 
      in->complCode = OBIT_CompReasonMinFract;*/
      break;  /* jump out of CLEAN loop */
    }

    /* Deep enough fraction of peak residual */
    if (fabs(xflux)<ccfLim) {
      g_snprintf (reason, 50, "Reached min fract of peak resid");
      in->complCode = OBIT_CompReasonMinFract;
      break;
    }
   
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
    }
  } /* end loop closing tables */
  
    /* Cleanup */
  CCRow = ObitTableCCRowUnref(CCRow);  
  KillCLEANArgs (nThreads, targs);
   
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

  return done;
} /* end ObitDConCleanPxListWBCLEAN */

/**
 * Tells results of CLEANing
 * \param inn         The Pixel List object 
 * \param ncomps      Array of total number of components per field, 
 *                    same order as in ObitImageMosaic
 * \param err Obit error stack object.
 */
olong ObitDConCleanPxListWBResult (ObitDConCleanPxList *inn, olong *ncomp,
				   ObitErr *err)
{
  olong out = 0;
  olong i;
  ObitDConCleanPxListWB *in = (ObitDConCleanPxListWB*)inn;

  /* error checks */
  if (err->error) return out;
  g_assert (ObitIsA(in, &myClassInfo));

  out = in->currentIter;
  for (i=0; i<in->nfield; i++) ncomp[i] = in->iterField[i];

  return out;
} /* end ObitDConCleanPxListWBResult */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitDConCleanPxListWBClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitDConCleanPxListWBClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitDConCleanPxListWBClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitDConCleanPxListWBClassInfoDefFn (gpointer inClass)
{
  ObitDConCleanPxListWBClassInfo *theClass = (ObitDConCleanPxListWBClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitDConCleanPxListWBClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitDConCleanPxListWBClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitDConCleanPxListWBGetClass;
  theClass->newObit       = (newObitFP)newObitDConCleanPxListWB;
  theClass->ObitCopy      = (ObitCopyFP)ObitDConCleanPxListWBCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitDConCleanPxListWBClear;
  theClass->ObitInit      = (ObitInitFP)ObitDConCleanPxListWBInit;
  theClass->ObitDConCleanPxListCreate   = 
    (ObitDConCleanPxListCreateFP)ObitDConCleanPxListWBCreate;
  theClass->ObitDConCleanPxListGetParms = 
    (ObitDConCleanPxListGetParmsFP)ObitDConCleanPxListWBGetParms;
  theClass->ObitDConCleanPxListReset    = 
    (ObitDConCleanPxListResetFP)ObitDConCleanPxListWBReset;
  theClass->ObitDConCleanPxListResize   = 
    (ObitDConCleanPxListResizeFP)ObitDConCleanPxListWBResize;
  theClass->ObitDConCleanPxListUpdate   = 
    (ObitDConCleanPxListUpdateFP)ObitDConCleanPxListWBUpdate;
  theClass->ObitDConCleanPxListCLEAN    = 
    (ObitDConCleanPxListCLEANFP)ObitDConCleanPxListWBCLEAN;
  theClass->ObitDConCleanPxListResult   = 
    (ObitDConCleanPxListResultFP)ObitDConCleanPxListWBResult;

} /* end ObitDConCleanPxListWBClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitDConCleanPxListWBInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitDConCleanPxListWB *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->BeamPatch00 = NULL;
  in->BeamPatch01 = NULL;
  in->BeamPatch02 = NULL;
  in->BeamPatch10 = NULL;
  in->BeamPatch11 = NULL;
  in->BeamPatch12 = NULL;
  in->BeamPatch20 = NULL;
  in->BeamPatch21 = NULL;
  in->BeamPatch22 = NULL;
  in->pixelFlux1  = NULL;
  in->pixelFlux2  = NULL;
  in->pixelFlux3  = NULL;

} /* end ObitDConCleanPxListWBInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitDConCleanPxListWB* cast to an Obit*.
 */
void ObitDConCleanPxListWBClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  olong i;
  ObitDConCleanPxListWB *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  if (in->BeamPatch00) {
    for (i=0; i<in->nfield; i++) {
      in->BeamPatch00[i] = ObitFArrayUnref(in->BeamPatch00[i]);
    }
    g_free(in->BeamPatch00);   in->BeamPatch00 = NULL;
  }
  if (in->BeamPatch01) {
    for (i=0; i<in->nfield; i++) {
      in->BeamPatch01[i] = ObitFArrayUnref(in->BeamPatch01[i]);
    }
    g_free(in->BeamPatch01);   in->BeamPatch01 = NULL;
  }
  if (in->BeamPatch02) {
    for (i=0; i<in->nfield; i++) {
      in->BeamPatch02[i] = ObitFArrayUnref(in->BeamPatch02[i]);
    }
    g_free(in->BeamPatch02);   in->BeamPatch02 = NULL;
  }
  if (in->BeamPatch10) {
    for (i=0; i<in->nfield; i++) {
      in->BeamPatch10[i] = ObitFArrayUnref(in->BeamPatch10[i]);
    }
    g_free(in->BeamPatch10);   in->BeamPatch10 = NULL;
  }
  if (in->BeamPatch11) {
    for (i=0; i<in->nfield; i++) {
      in->BeamPatch11[i] = ObitFArrayUnref(in->BeamPatch11[i]);
    }
    g_free(in->BeamPatch11);   in->BeamPatch11 = NULL;
  }
  if (in->BeamPatch12) {
    for (i=0; i<in->nfield; i++) {
      in->BeamPatch12[i] = ObitFArrayUnref(in->BeamPatch12[i]);
    }
    g_free(in->BeamPatch12);   in->BeamPatch12 = NULL;
  }
  if (in->BeamPatch20) {
    for (i=0; i<in->nfield; i++) {
      in->BeamPatch20[i] = ObitFArrayUnref(in->BeamPatch20[i]);
    }
    g_free(in->BeamPatch20);   in->BeamPatch20 = NULL;
  }
  if (in->BeamPatch21) {
    for (i=0; i<in->nfield; i++) {
      in->BeamPatch21[i] = ObitFArrayUnref(in->BeamPatch21[i]);
    }
    g_free(in->BeamPatch21);   in->BeamPatch21 = NULL;
  }
  if (in->BeamPatch22) {
    for (i=0; i<in->nfield; i++) {
      in->BeamPatch22[i] = ObitFArrayUnref(in->BeamPatch22[i]);
    }
    g_free(in->BeamPatch22);   in->BeamPatch22 = NULL;
  }
  if (in->pixelFlux1) {
    g_free(in->pixelFlux1);   in->pixelFlux1 = NULL;
  }
  if (in->pixelFlux2) {
    g_free(in->pixelFlux2);   in->pixelFlux2 = NULL;
  }
  if (in->pixelFlux3) {
    g_free(in->pixelFlux3);   in->pixelFlux3 = NULL;
  }


  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitDConCleanPxListWBClear */


/**
 * Inner loop BGC CLEAN, 0th order
 * Subtracts component from residuals using beam patch and
 * returns largest residual in region cleaned
 * Callable as thread
 * Arguments are given in the structure passed as arg
 * \param arg Pointer to CLEANFuncArg argument with elements:
 * \li PixelList PixelList object
 * \li first   First (1-rel) pixel no. to process this thread
 * \li last    Highest (1-rel) pixel no. to process this thread
 * \li ithread thread number, <0-> no threads
 * \li ipeak   [in/out] Pixel number of peak
 * \li peak;   [in/out] Peak flux to subtract
 * \return NULL
 */
gpointer ThreadCLEAN0 (gpointer args)
{
  CLEANFuncArg *largs  = (CLEANFuncArg*)args;
  ObitDConCleanPxListWB *in = largs->PixelList;
  olong  loPix            = largs->first-1;
  olong  hiPix            = largs->last;
  olong  ipeak            = largs->ipeak;
  ofloat peak             = largs->peak;

  ofloat xflux, subval;
  olong iresid, iXres, iYres, lpatch, beamPatch, iBeam, field, pos[2];
  ofloat *beam00=NULL;
 
  /* Do subtraction if peak non zero */
  field  = in->pixelFld[ipeak];
  lpatch = in->BeamPatch00[field-1]->naxis[0];
  beamPatch = lpatch/2;
  if (fabs(peak)>0.0) {
    pos[0] = pos[1] = 0;
    beam00 = ObitFArrayIndex(in->BeamPatch00[field-1], pos); /* Beam patch pointer */
    xflux  = in->pixelFlux[ipeak];
    subval = xflux * in->gain[field-1];
    iXres  = in->pixelX[ipeak];
    iYres  = in->pixelY[ipeak];
    for (iresid=loPix; iresid<hiPix; iresid++) {
      /* Is this inside the Beam patch ? */
      if ((abs(in->pixelX[iresid]-iXres) <= beamPatch) && 
	  (abs(in->pixelY[iresid]-iYres) <= beamPatch)) {
	/*&& (in->pixelFld[iresid]==field)) {*/
	/* Index in beam patch array */
	iBeam = (beamPatch + (in->pixelY[iresid] - iYres)) * lpatch +
	  (beamPatch + (in->pixelX[iresid] - iXres));
	in->pixelFlux[iresid] -= subval * beam00[iBeam];
      }
      /* May not be ordered if (in->pixelY[iresid]-iYres > (beamPatch+5)) break; No more in Y? */
    } /* end loop over array */
  }
    
  /* Find peak abs residual */
  peak  = -1.0;
  xflux  = 0.0;
  for (iresid=loPix; iresid<hiPix; iresid++) {
    xflux = fabs(in->pixelFlux[iresid]);
    if (xflux>peak) {
      peak = xflux;
      ipeak = iresid;
    }
  }

  /* Return values */
  largs->peak   = in->pixelFlux[ipeak];
  largs->si     = 0.0;
  largs->curve  = 0.0;
  largs->curve2 = 0.0;
  largs->ipeak  = ipeak;
  
  /* Indicate completion */
  if (largs->ithread>=0)
    ObitThreadPoolDone (in->thread, (gpointer)&largs->ithread);
  return NULL;
} /* end ThreadCLEAN0 */

/**
 * Inner loop BGC CLEAN, 1st order
 * Subtracts component from residuals using beam patch and
 * returns largest residual in region cleaned
 * Callable as thread
 * Arguments are given in the structure passed as arg
 * \param arg Pointer to CLEANFuncArg argument with elements:
 * \li PixelList PixelList object
 * \li first   First (1-rel) pixel no. to process this thread
 * \li last    Highest (1-rel) pixel no. to process this thread
 * \li ithread thread number, <0-> no threads
 * \li ipeak   [in/out] Pixel number of peak
 * \li peak    [in/out] Peak flux to subtract
 * \li si      [in/out] Spectral index (times peak)
 * \return NULL
 */
gpointer ThreadCLEAN1 (gpointer args)
{
  CLEANFuncArg *largs  = (CLEANFuncArg*)args;
  ObitDConCleanPxListWB *in = largs->PixelList;
  olong  loPix            = largs->first-1;
  olong  hiPix            = largs->last;
  olong  ipeak            = largs->ipeak;
  ofloat peak             = largs->peak;
  ofloat si               = largs->si;

  ofloat xflux, subval, subval1;
  olong iresid, iXres, iYres, lpatch, beamPatch, iBeam, field, pos[2];
  ofloat *beam00=NULL, *beam10=NULL, *beam01=NULL, *beam11=NULL;
  odouble iDelta=1.0, A00=0.0, A10=0.0, A01=0.0, A11=0.0;
 
  /* Do subtraction if peak non zero */
  field  = in->pixelFld[ipeak];
  lpatch = in->BeamPatch00[field-1]->naxis[0];
  beamPatch = lpatch/2;
  if (fabs(peak)>0.0) {
    pos[0]  = pos[1] = 0;
    beam00  = ObitFArrayIndex(in->BeamPatch00[field-1], pos); /* Beam patch pointer */
    beam01  = ObitFArrayIndex(in->BeamPatch01[field-1], pos); 
    beam10  = ObitFArrayIndex(in->BeamPatch10[field-1], pos); 
    beam11  = ObitFArrayIndex(in->BeamPatch11[field-1], pos); 
    subval  = peak * in->gain[field-1];
    subval1 = si * in->gain[field-1];
    iXres   = in->pixelX[ipeak];
    iYres   = in->pixelY[ipeak];
    for (iresid=loPix; iresid<hiPix; iresid++) {
      /* Is this inside the Beam patch ? */
      if ((abs(in->pixelX[iresid]-iXres) <= beamPatch) && 
	  (abs(in->pixelY[iresid]-iYres) <= beamPatch)) {
	/*&& (in->pixelFld[iresid]==field)) {*/
	/* Index in beam patch array */
	iBeam = (beamPatch + (in->pixelY[iresid] - iYres)) * lpatch +
	  (beamPatch + (in->pixelX[iresid] - iXres));
	/* Sault-Wieringa values */
	in->pixelFlux[iresid]  -= subval*beam00[iBeam] + subval1*beam10[iBeam];
	in->pixelFlux1[iresid] -= subval*beam01[iBeam] + subval1*beam11[iBeam];
      }
      /* May not be ordered if (in->pixelY[iresid]-iYres > (beamPatch+5)) break; No more in Y? */
    } /* end loop over array */
  }
   
  /* Find peak abs residual */
  peak  = -1.0;
  si     = 0.0;
  xflux  = 0.0;
  field = -1;
  for (iresid=loPix; iresid<hiPix; iresid++) {
    /* Center of new beam? */
    if (field!=in->pixelFld[ipeak]) {
      field  = in->pixelFld[ipeak];
      lpatch = in->BeamPatch00[field-1]->naxis[0];
      beamPatch = lpatch/2;
      pos[0] = lpatch/2; pos[1] = lpatch/2;
      A00 = *ObitFArrayIndex(in->BeamPatch00[field-1], pos);
      A01 = *ObitFArrayIndex(in->BeamPatch01[field-1], pos);
      A10 = *ObitFArrayIndex(in->BeamPatch10[field-1], pos);
      A11 = *ObitFArrayIndex(in->BeamPatch11[field-1], pos);
      iDelta = 1.0 / (A00*A11 - A01*A10);  /* Determinant */
    }
    /* Sault-Wieringa  eq 22 values for max impact */
    xflux = in->pixelFlux[iresid]*in->pixelFlux[iresid]*A11 +
      in->pixelFlux1[iresid]*in->pixelFlux1[iresid]*A00 -
      2.0 * in->pixelFlux[iresid]*in->pixelFlux1[iresid] * A01;
    /*xflux = fabs((A11*in->pixelFlux[iresid]  - A01*in->pixelFlux1[iresid])); max flux */
    /*xflux = in->pixelFlux[ipeak]; DEBUG */
    if (xflux>peak) {
      peak = xflux;
      ipeak = iresid;
    }
  }

  /* Sault-Wieringa values */
  peak = (A11*in->pixelFlux[ipeak]  - A01*in->pixelFlux1[ipeak]) * iDelta;
  si   = (A00*in->pixelFlux1[ipeak] - A10*in->pixelFlux[ipeak])  * iDelta;

  /* Sanity check - resolution will appear to steepen spectrum */
  if ((si/peak)> 2.0)    si =  0.0;
  if ((si/peak)<-10.0)   si =  0.0;
  if ((si/peak)>1.0)     si =  peak;
  if ((si/peak)<-5.)     si = -5.0*peak;
  
  /* Return values */
  largs->peak   = peak;
  largs->si     = si;
  largs->curve  = 0.0;
  largs->curve2 = 0.0;
  largs->ipeak  = ipeak;
  
  /* Indicate completion */
  if (largs->ithread>=0)
    ObitThreadPoolDone (in->thread, (gpointer)&largs->ithread);
  return NULL;
} /* end ThreadCLEAN1 */

/**
 * Inner loop BGC CLEAN - 2nd order, flux, si, curve
 * Subtracts component from residuals using beam patch and
 * returns largest residual in region cleaned
 * Callable as thread
 * Arguments are given in the structure passed as arg
 * \param arg Pointer to CLEANFuncArg argument with elements:
 * \li PixelList PixelList object
 * \li first   First (1-rel) pixel no. to process this thread
 * \li last    Highest (1-rel) pixel no. to process this thread
 * \li ithread thread number, <0-> no threads
 * \li ipeak   [in/out] Pixel number of peak
 * \li peak    [in/out] Peak flux to subtract 
 * \li si      [in/out] Spectral index (times peak)
 * \li curve   [in/out] Spectral curvature (times peak)
 * \return NULL
 */
gpointer ThreadCLEAN2 (gpointer args)
{
  CLEANFuncArg *largs  = (CLEANFuncArg*)args;
  ObitDConCleanPxListWB *in = largs->PixelList;
  olong  loPix            = largs->first-1;
  olong  hiPix            = largs->last;
  olong  ipeak            = largs->ipeak;
  ofloat peak             = largs->peak;
  ofloat si               = largs->si;
  ofloat curve            = largs->curve;

  ofloat xflux, subval, subval1, subval2, iDelta=1.0;
  olong iresid, iXres, iYres, lpatch, beamPatch, iBeam, field, pos[2];
  ofloat *beam00=NULL, *beam10=NULL, *beam01=NULL, *beam11=NULL;
  ofloat *beam02=NULL, *beam20=NULL, *beam21=NULL, *beam12=NULL, *beam22=NULL;
  odouble A00=0.0, A10=0.0, A01=0.0, A11=0.0, A20=0.0, A02=0.0, A21=0.0, A12=0.0, A22=0.0, R0, R1, R2;

  /* Do subtraction if peak non zero */
  field  = in->pixelFld[ipeak];
  lpatch = in->BeamPatch00[field-1]->naxis[0];
  beamPatch = lpatch/2;
  if (fabs(peak)>0.0) {
    pos[0] = pos[1] = 0;
    beam00  = ObitFArrayIndex(in->BeamPatch00[field-1], pos); /* Beam patch pointer */
    beam01  = ObitFArrayIndex(in->BeamPatch01[field-1], pos);
    beam10  = ObitFArrayIndex(in->BeamPatch10[field-1], pos);
    beam11  = ObitFArrayIndex(in->BeamPatch11[field-1], pos);
    beam20  = ObitFArrayIndex(in->BeamPatch20[field-1], pos);
    beam02  = ObitFArrayIndex(in->BeamPatch02[field-1], pos);
    beam12  = ObitFArrayIndex(in->BeamPatch12[field-1], pos);
    beam21  = ObitFArrayIndex(in->BeamPatch21[field-1], pos);
    beam22  = ObitFArrayIndex(in->BeamPatch22[field-1], pos);
    subval  = peak * in->gain[field-1];
    subval1 = si * in->gain[field-1];
    subval2 = curve * in->gain[field-1];
    iXres  = in->pixelX[ipeak];
    iYres  = in->pixelY[ipeak];
    for (iresid=loPix; iresid<hiPix; iresid++) {
      /* Is this inside the Beam patch ? */
      if ((abs(in->pixelX[iresid]-iXres) <= beamPatch) && 
	  (abs(in->pixelY[iresid]-iYres) <= beamPatch)) {
	/* && (in->pixelFld[iresid]==field)) {*/
	/* Index in beam patch array */
	iBeam = (beamPatch + (in->pixelY[iresid] - iYres)) * lpatch +
	  (beamPatch + (in->pixelX[iresid] - iXres));
	in->pixelFlux[iresid]  -= subval*beam00[iBeam] + subval1*beam10[iBeam] + subval2*beam20[iBeam];
	in->pixelFlux1[iresid] -= subval*beam01[iBeam] + subval1*beam11[iBeam] + subval2*beam21[iBeam]; 
	in->pixelFlux2[iresid] -= subval*beam02[iBeam] + subval1*beam12[iBeam] + subval2*beam22[iBeam]; 
      }
      /* May not be ordered if (in->pixelY[iresid]-iYres > (beamPatch+5)) break; No more in Y? */
    } /* end loop over array */
  }
    
  /* Find peak abs residual */
  peak  = -1.0;
  si     = 0.0;
  curve  = 0.0;
  xflux  = 0.0;
  field = -1;
  for (iresid=loPix; iresid<hiPix; iresid++) {
    /* Center of new beam? */
    if (field!=in->pixelFld[ipeak]) {
      field  = in->pixelFld[ipeak];
      lpatch = in->BeamPatch00[field-1]->naxis[0];
      beamPatch = lpatch/2;
      pos[0] = lpatch/2; pos[1] = lpatch/2;
      A00 = *ObitFArrayIndex(in->BeamPatch00[field-1], pos);
      A01 = *ObitFArrayIndex(in->BeamPatch01[field-1], pos);
      A10 = *ObitFArrayIndex(in->BeamPatch10[field-1], pos);
      A11 = *ObitFArrayIndex(in->BeamPatch11[field-1], pos);
      A20 = *ObitFArrayIndex(in->BeamPatch20[field-1], pos);
      A02 = *ObitFArrayIndex(in->BeamPatch02[field-1], pos);
      A21 = *ObitFArrayIndex(in->BeamPatch21[field-1], pos);
      A12 = *ObitFArrayIndex(in->BeamPatch12[field-1], pos);
      A22 = *ObitFArrayIndex(in->BeamPatch22[field-1], pos);
      iDelta = 1.0 / (A02*A02*A11 - 2.0*A01*A02*A12 + A01*A01*A22 + A00*(A12*A12-A11*A22));
    }
    R0 = in->pixelFlux[iresid];
    R1 = in->pixelFlux1[iresid];
    R2 = in->pixelFlux2[iresid];
    /* Max flux 
    xflux = fabs(A12*A12*R0 + A01*A22*R1 - A12*( A02*R1 + A01*R2) + A11*(-A22*R0 + A02*R2));*/
    /* Sault-Wieringa eq 22 values */
    xflux = A12*A12*R0*R0 + (A02*A02-A00*A22)*R1*R1 + A01*A01*R2*R2 +
      2.0*A01*R1*(A22*R0 - A02*R2) - 2.0*A12*(A02*R0*R1 + (A01*R0-A00*R1)*R2) -
      A11*(A22*R0*R0 + R2*(-2.0*A02*R0 + A00*R2)); 
    if (fabs(xflux)>peak) {
      peak = fabs(xflux);
      ipeak = iresid;
    }
  }
  
  /* Sault-Wieringa values */
  R0 = in->pixelFlux[ipeak];
  R1 = in->pixelFlux1[ipeak];
  R2 = in->pixelFlux2[ipeak];
  peak  =  A12*A12*R0 + A01*A22*R1 - A12*( A02*R1 + A01*R2) + A11*(-A22*R0 + A02*R2);
  si    =  A01*A22*R0 + A02*A02*R1 - A02*( A12*R0 + A01*R2) + A00*(-A22*R1 + A12*R2);
  curve = -A01*A12*R0 + A02*(A11*R0 - A01*R1) + A01*A01*R2 + A00*(A12*R1 - A11*R2);

  /* Sanity check resolution will appear to steepen spectrum */
  if ((si/peak)> 2.0)   {si =  0.0; curve=0.0;}
  if ((si/peak)<-10.0)   {si =  0.0; curve=0.0;}
  if ((si/peak)>1.0)     si =  peak;
  if ((si/peak)<-5.)     si = -5.0*peak;

  /*if ((curve/peak)> 0.5) {curve = 0.0; si =  0.0; }
    if ((curve/peak)<-1.0) {curve = 0.0; si =  0.0; } DEBUG*/
  if ((curve/peak)>0.25) curve =  0.25*peak;
  if ((curve/peak)<-0.5) curve = -0.5*peak; 

  /* Return values */
  largs->peak   = peak*iDelta;
  largs->si     = si*iDelta;
  largs->curve  = curve*iDelta;
  largs->curve2 = 0.0;
  largs->ipeak  = ipeak;
  
  /* Indicate completion */
  if (largs->ithread>=0)
    ObitThreadPoolDone (in->thread, (gpointer)&largs->ithread);
  return NULL;
} /* end ThreadCLEAN2 */

/**
 * Make arguments for Threaded CLEAN
 * \param in         OTF with internal buffer to be modified.
 * \param maxThread  Maximum desirable no. threads
 * \param args       Created array of CLEANFuncArg, 
 *                   delete with KillCLEANArgs
 * \return number of elements in args.
 */
static olong MakeCLEANArgs (ObitDConCleanPxListWB *in, olong maxThread,
			    CLEANFuncArg ***args)
{
  olong i, nThreads;

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
    (*args)[i]->peak      = 0.0;
    (*args)[i]->si        = 0.0;
    (*args)[i]->curve     = 0.0;
    (*args)[i]->curve2    = 0.0;
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
    if (args[i]) g_free(args[i]);
  }
  g_free(args);
} /*  end KillCLEANArgs */

/**
 * Calculates equivalent exponential spectral index from linear
 * Averages equivalent at low and high frequencies.
 * \param in      The Pixel List object 
 * \param a       Linear spectral index
 * \param alpha   [out] equivalent exponential term
 */
void Linear2Exp1 (ObitDConCleanPxListWB *in, ofloat a, ofloat *alpha)
{
  ofloat l1, temp1;

  l1 = log (in->loFreq/in->refFreq);

  /* Keep this from blowing up */
  temp1 = MAX (a*l1, -0.8);

  *alpha = log(1.0+temp1)/l1;
  *alpha = a;
  /* DEBUG 
  fprintf (stderr,"Equivalent exponential si %f\n", *alpha);*/
} /* end Linear2Exp1 */

/**
 * Calculates equivalent exponential spectral index and curvature from linear
 * \param in      The Pixel List object 
 * \param a       Linear spectral index
 * \param a       Linear spectral curvature
 * \param alpha   [out] equivalent exponential spectral index term
 * \param alpha   [out] equivalent exponential curvature term
 */
void Linear2Exp2 (ObitDConCleanPxListWB *in, ofloat a, ofloat b, ofloat *alpha, ofloat *beta)
{
  ofloat l1, temp1;

  l1 = log (in->loFreq/in->refFreq);

  /* Keep this from blowing up */
  temp1 = MAX (a*l1, -0.8);

  *alpha = log(1.0+temp1)/l1;
  *alpha = a;
 *beta = b;
  /*  *beta = (log(1.0+temp1+b*l2*l2)*l2 - log(1.0+temp2+b*l2*l2)*l1)/ (l1*l1*l2 - l2*l2*l1);*/
  /* DEBUG
  fprintf (stderr,"Equivalent exponential si %f curve %f \n", *alpha, *beta); */
} /* end Linear2Exp2 */

