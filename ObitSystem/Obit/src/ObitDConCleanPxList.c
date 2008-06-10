/* $Id: ObitDConCleanPxList.c,v 1.21 2008/03/14 12:16:26 bcotton Exp $ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2004-2008                                          */
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

#include "ObitDConCleanPxList.h"
#include "ObitMem.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitDConCleanPxList.c
 * ObitDConCleanPxList class function definitions.
 * This class determines the pixel histogram of an image.
 * This class is derived from the Obit base class.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitDConCleanPxList";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo structure ObitDConCleanPxListClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitDConCleanPxListClassInfo myClassInfo = {FALSE};

/*--------------- File Global Variables  ----------------*/

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitDConCleanPxListInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitDConCleanPxListClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitDConCleanPxListClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitDConCleanPxList* newObitDConCleanPxList (gchar* name)
{
  ObitDConCleanPxList* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitDConCleanPxListClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitDConCleanPxList));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitDConCleanPxListInit((gpointer)out);

 return out;
} /* end newObitDConCleanPxList */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitDConCleanPxListGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitDConCleanPxListClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitDConCleanPxListGetClass */

/**
 * Make a deep copy of an ObitDConCleanPxList.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitDConCleanPxList* 
ObitDConCleanPxListCopy  (ObitDConCleanPxList *in, ObitDConCleanPxList *out, 
			  ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  gboolean oldExist;
  gchar *outName;
  olong i, nfield=0;
  gchar *routine = "ObitDConCleanPxListCopy";

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
    out = newObitDConCleanPxList(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  out->nfield = in->nfield;
  /* Free old */
  out->info = ObitInfoListUnref(out->info);
  out->mosaic = ObitImageMosaicUnref(out->mosaic);
  out->window = ObitDConCleanWindowUnref(out->window);
  out->BeamPatch = ObitFArrayUnref(out->BeamPatch);
  if (out->pixelX)    out->pixelX   =  ObitMemFree (out->pixelX);
  if (out->pixelY)    out->pixelY   =  ObitMemFree (out->pixelY);
  if (out->pixelFld)  out->pixelFld =  ObitMemFree (out->pixelFld);
  if (out->pixelFlux) out->pixelFlux=  ObitMemFree (out->pixelFlux);
  if (out->iterField) out->iterField=  ObitMemFree (out->iterField);
  if (out->CCver)     out->CCver    =  ObitMemFree (out->CCver);
  if (out->fluxField) out->fluxField=  ObitMemFree (out->fluxField);
  if (out->circGaus)  out->circGaus =  ObitMemFree (out->circGaus);
  if (out->gain)      out->gain     =  ObitMemFree (out->gain);
  if (out->minFlux)   out->minFlux  =  ObitMemFree (out->minFlux);
  if (out->factor)    out->factor   =  ObitMemFree (out->factor);
  if (out->CCTable) {
    nfield = in->mosaic->numberImages;
    for (i=0; i<nfield; i++) {
      out->CCTable[i] = ObitTableCCUnref(out->CCTable[i]);
    }
    out->CCTable   =  ObitMemFree (out->CCTable);
  }

  /* Copy new */
  out->info      = ObitInfoListCopy(in->info);
  out->window    = ObitDConCleanWindowCopy(in->window, out->window, err);
  out->mosaic    = ObitImageMosaicCopy(in->mosaic, out->mosaic, err);
  out->BeamPatch = ObitFArrayCopy(in->BeamPatch, out->BeamPatch, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, out);

  /* Create arrays */
  out->maxPixel = in->maxPixel;
  out->pixelX    = ObitMemAlloc0Name (out->maxPixel*sizeof(olong),  "PxList X pixel");
  out->pixelY    = ObitMemAlloc0Name (out->maxPixel*sizeof(olong),  "PxList Y pixel");
  out->pixelFld  = ObitMemAlloc0Name (out->maxPixel*sizeof(gshort), "PxList pixel field");
  out->pixelFlux = ObitMemAlloc0Name (out->maxPixel*sizeof(ofloat), "PxList pixel Flux");
  for (i=0; i<out->maxPixel; i++) {
    out->pixelX[i]    = in->pixelX[i];
    out->pixelY[i]    = in->pixelY[i];
    out->pixelFld[i]  = in->pixelFld[i];
    out->pixelFlux[i] = in->pixelFlux[i];
  }

  /* Per field */
  out->fluxField = ObitMemAlloc0Name (nfield*sizeof(ofloat), "PxList Clean Flux");
  out->circGaus  = ObitMemAlloc0Name (nfield*sizeof(ofloat), "PxList Gaussian");
  out->iterField = ObitMemAlloc0Name (nfield*sizeof(olong),  "PxList Clean CC count");
  out->CCver     = ObitMemAlloc0Name (nfield*sizeof(olong),  "PxList Clean CC version");
  out->gain      = ObitMemAlloc0Name (nfield*sizeof(ofloat), "PxList Clean Loop gain");
  out->minFlux   = ObitMemAlloc0Name (nfield*sizeof(ofloat), "PxList Clean Min flux");
  out->factor    = ObitMemAlloc0Name (nfield*sizeof(ofloat), "PxList Clean factor");
  out->CCTable   = ObitMemAlloc0Name (nfield*sizeof(ObitTableCC*), "PxList CC tables");
  for (i=0; i<nfield; i++) {
    out->fluxField[i] = in->fluxField[i];
    out->circGaus[i]  = in->circGaus[i];
    out->iterField[i] = in->iterField[i];
    out->CCver[i]     = in->CCver[i];
    out->gain[i]      = in->gain[i];
    out->minFlux[i]   = in->minFlux[i];
    out->factor[i]    = in->factor[i];
    if(in->CCTable[i]) out->CCTable[i] = ObitTableCCRef(in->CCTable[i]);
  }

  return out;
} /* end ObitDConCleanPxListCopy */

/**
 * Make a copy of a object but do not copy the actual data
 * This is useful to create an DConCleanPxList similar to the input one.
 * \param in  The object to copy
 * \param out An existing object pointer for output, must be defined.
 * \param err Obit error stack object.
 */
void ObitDConCleanPxListClone  (ObitDConCleanPxList *in, ObitDConCleanPxList *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  olong i, nfield;
  gchar *routine = "ObitDConCleanPxListClone";

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
  out->nfield = in->nfield;
  /* Free old */
  out->info = ObitInfoListUnref(out->info);
  out->mosaic = ObitImageMosaicUnref(out->mosaic);
  out->BeamPatch = ObitFArrayUnref(out->BeamPatch);
  if (out->pixelX)    out->pixelX   =  ObitMemFree (out->pixelX);
  if (out->pixelY)    out->pixelY   =  ObitMemFree (out->pixelY);
  if (out->pixelFld)  out->pixelFld =  ObitMemFree (out->pixelFld);
  if (out->pixelFlux) out->pixelFlux=  ObitMemFree (out->pixelFlux);
  if (out->iterField) out->iterField=  ObitMemFree (out->iterField);
  if (out->CCver)     out->CCver    =  ObitMemFree (out->CCver);
  if (out->fluxField) out->fluxField=  ObitMemFree (out->fluxField);
  if (out->circGaus)  out->circGaus =  ObitMemFree (out->circGaus);
  if (out->gain)      out->gain     =  ObitMemFree (out->gain);
  if (out->minFlux)   out->minFlux  =  ObitMemFree (out->minFlux);
  if (out->factor)    out->factor   =  ObitMemFree (out->factor);
  if (out->CCTable) {
    nfield = in->mosaic->numberImages;
    for (i=0; i<nfield; i++) {
      out->CCTable[i] = ObitTableCCUnref(out->CCTable[i]);
    }
    out->CCTable   =  ObitMemFree (out->CCTable);
  }

  /* Copy new */
  out->info = ObitInfoListCopy(in->info);
  ObitDConCleanWindowClone(in->window, out->window, err);
  out->mosaic = ObitImageMosaicRef(in->mosaic);
  out->BeamPatch = ObitFArrayRef(in->BeamPatch);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Create arrays */
  out->maxPixel = in->maxPixel;
  out->pixelX    = ObitMemAlloc0Name (out->maxPixel*sizeof(olong),  "PxList X pixel");
  out->pixelY    = ObitMemAlloc0Name (out->maxPixel*sizeof(olong),  "PxList Y pixel");
  out->pixelFld  = ObitMemAlloc0Name (out->maxPixel*sizeof(gshort), "PxList pixel field");
  out->pixelFlux = ObitMemAlloc0Name (out->maxPixel*sizeof(ofloat), "PxList pixel Flux");
  for (i=0; i<out->maxPixel; i++) {
    out->pixelX[i]    = in->pixelX[i];
    out->pixelY[i]    = in->pixelY[i];
    out->pixelFld[i]  = in->pixelFld[i];
    out->pixelFlux[i] = in->pixelFlux[i];
  }

  /* Per field */
  nfield = in->mosaic->numberImages;
  out->fluxField = ObitMemAlloc0Name (nfield*sizeof(ofloat), "PxList Clean Flux");
  out->circGaus  = ObitMemAlloc0Name (nfield*sizeof(ofloat), "PxList Gaussian");
  out->iterField = ObitMemAlloc0Name (nfield*sizeof(olong),  "PxList Clean CC count");
  out->CCver     = ObitMemAlloc0Name (nfield*sizeof(olong),  "PxList Clean CC version");
  out->gain      = ObitMemAlloc0Name (nfield*sizeof(ofloat), "PxList Clean Loop gain");
  out->minFlux   = ObitMemAlloc0Name (nfield*sizeof(ofloat), "PxList Clean Min flux");
  out->factor    = ObitMemAlloc0Name (nfield*sizeof(ofloat), "PxList Clean factor");
  out->CCTable   = ObitMemAlloc0Name (nfield*sizeof(ObitTableCC*), "PxList CC tables");
  for (i=0; i<nfield; i++) {
    out->fluxField[i] = in->fluxField[i];
    out->circGaus[i]  = in->circGaus[i];
    out->iterField[i] = in->iterField[i];
    out->CCver[i]     = in->CCver[i];
    out->gain[i]      = in->gain[i];
    out->minFlux[i]   = in->minFlux[i];
    out->factor[i]    = in->factor[i];
    if(in->CCTable[i]) out->CCTable[i] = ObitTableCCRef(in->CCTable[i]);
 }

} /* end ObitDConCleanPxListClone */

/**
 * Creates an ObitDConCleanPxList
 * \param name     An optional name for the object.
 * \param mosaic   Image mosaic to be deconvolved.
 * \param maxPixel Maximum number of pixels allowed (dim. or arrays)
 * \return the new object.
 */
ObitDConCleanPxList* 
ObitDConCleanPxListCreate (gchar* name, ObitImageMosaic *mosaic,  
			   olong maxPixel, ObitErr *err)
{
  ObitDConCleanPxList* out=NULL;
  olong i, nfield;
  /*gchar *routine = "ObitDConCleanPxListCreate";*/

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return out;
  g_assert (ObitImageMosaicIsA(mosaic));

  /* Create basic structure */
  out = newObitDConCleanPxList (name);

  /* Create arrays */
  out->maxPixel = maxPixel;
  out->pixelX    = ObitMemAlloc0Name (maxPixel*sizeof(olong),  "PxList X pixel");
  out->pixelY    = ObitMemAlloc0Name (maxPixel*sizeof(olong),  "PxList Y pixel");
  out->pixelFld  = ObitMemAlloc0Name (maxPixel*sizeof(gshort), "PxList pixel field");
  out->pixelFlux = ObitMemAlloc0Name (maxPixel*sizeof(ofloat), "PxList pixel Flux");

  /* Save Image Mosaic reference */
  out->mosaic = ObitImageMosaicRef(mosaic);

  /* Per field */
  nfield      = mosaic->numberImages;
  out->nfield = nfield;
  out->iterField = ObitMemAlloc0Name (nfield*sizeof(olong),  "PxList Clean CC count");
  out->CCver     = ObitMemAlloc0Name (nfield*sizeof(olong),  "PxList Clean CC version");
  out->fluxField = ObitMemAlloc0Name (nfield*sizeof(ofloat), "PxList Clean Flux");
  out->circGaus  = ObitMemAlloc0Name (nfield*sizeof(ofloat), "PxList Gaussian");
  out->gain      = ObitMemAlloc0Name (nfield*sizeof(ofloat), "PxList Clean Loop gain");
  out->minFlux   = ObitMemAlloc0Name (nfield*sizeof(ofloat), "PxList Clean Mix flux");
  out->factor    = ObitMemAlloc0Name (nfield*sizeof(ofloat), "PxList Clean factor");
  out->CCTable   = ObitMemAlloc0Name (nfield*sizeof(ObitTableCC*), "PxList CC tables");
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
  

  return out;
} /* end ObitDConCleanPxListCreate */

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
 * \li "CCVer"   OBIT_long array    = CLEAN table version per field
 *                                   If only one given it is used for all.
 * \li "prtLv"   OBIT_long          = message level  [def 2]
 *               0=none, 1=summary, 2=normal, higher numbers for diagnostics
 * \param in  The Pixel list CLEAN object
 * \param err Obit error stack object.
 */
void  ObitDConCleanPxListGetParms (ObitDConCleanPxList *in, ObitErr *err)
{
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  olong i, itemp, *ipnt=NULL;
  olong jtemp=0;
  union ObitInfoListEquiv InfoReal; 
  /* gchar *routine = "ObitDConCleanPxListGetParms";*/

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Niter */
  InfoReal.itg = 0;type = OBIT_oint;
  ObitInfoListGetTest(in->info, "Niter", &type, (gint32*)dim, &InfoReal);
  if (type==OBIT_float) itemp = InfoReal.flt + 0.5;
  else itemp = InfoReal.itg;
  in->niter  = itemp;

  /* Loop CC version per field */
  InfoReal.itg = in->CCver[0]; dim[0]=1; type = OBIT_oint;
  if (ObitInfoListGetP(in->info, "CCVer", &type, dim, (gpointer)&ipnt)) {
    /* If only one, use for all */
    if ((in->nfield>1) && (dim[0]==1)) {
      if (type==OBIT_long) jtemp = (olong)ipnt[0];
      else if (type==OBIT_long) jtemp = ((olong*)ipnt)[0];
      else if (type==OBIT_oint) jtemp = ((oint*)ipnt)[0];
      for (i=0; i<in->nfield; i++) in->CCver[i] = jtemp;
    } else { /* read one per field */
       for (i=0; i<dim[0]; i++) {
	 if (type==OBIT_long) jtemp = (olong)ipnt[i];
	 else if (type==OBIT_long) jtemp = ((olong*)ipnt)[i];
	 else if (type==OBIT_oint) jtemp = ((oint*)ipnt)[i];
	 in->CCver[i] = jtemp;
       }
   }
  }

  /* Loop GAIN per field */
  ObitInfoListGetTest(in->info, "Gain", &type, (gint32*)dim, in->gain);
  /* If only one, use for all */
  if ((in->nfield>1) && (dim[0]==1))
    for (i=1; i<in->nfield; i++) in->gain[i] = in->gain[0];

  /* Minimum flux density per field */
  ObitInfoListGetTest(in->info, "minFlux", &type, (gint32*)dim, in->minFlux);
  /* If only one, use for all */
  if ((in->nfield>1) && (dim[0]==1))
    for (i=1; i<in->nfield; i++) in->minFlux[i] = in->minFlux[0];

  /* CLEAN depth factor per field */
  ObitInfoListGetTest(in->info, "Factor", &type, (gint32*)dim, in->factor);
  /* If only one, use for all */
  if ((in->nfield>1) && (dim[0]==1))
    for (i=1; i<in->nfield; i++) in->factor[i] = in->factor[0];

  /* Gaussian convolution per field */
  ObitInfoListGetTest(in->info, "fGauss", &type, (gint32*)dim, in->circGaus);

  /* Print level */
  in->prtLv = 2;  /* default = normal */
  ObitInfoListGetTest(in->info, "prtLv", &type, dim, &in->prtLv);
} /* end ObitDConCleanPxListGetParms */

/**
 * Resets CLEAN information
 * Makes sure all potential CC tables are instantiated
 * \param in          The Pixel List object 
 * \param err Obit error stack object.
 */
void ObitDConCleanPxListReset (ObitDConCleanPxList *in, ObitErr *err)
{
  olong i, ver=0;
  oint noParms;
  gchar *routine = " ObitDConCleanPxListReset";

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Control info */
  ObitDConCleanPxListGetParms(in, err);
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
      if (in->circGaus[i]>0.0) noParms = 4;
      else noParms = 0;
      ver = MAX (ver, in->CCver[i]);  /* Use last if not defined */
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

} /* end ObitDConCleanPxListReset */

/**
 * Resets CLEAN information
 * Makes sure all potential CC tables are instantiated
 * \param in       The Pixel List object 
 * \param maxPixel Maximum number of pixels allowed (dim. or arrays)
 * \param err      Obit error stack object.
 */
void ObitDConCleanPxListResize (ObitDConCleanPxList *in, olong maxPixel, 
				ObitErr *err)
{
  /*gchar *routine = " ObitDConCleanPxListResize";*/

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  if (maxPixel<in->maxPixel) return;  /* This needed? */

  in->maxPixel  = maxPixel;
  in->pixelX    = ObitMemRealloc (in->pixelX,    maxPixel*sizeof(olong));
  in->pixelY    = ObitMemRealloc (in->pixelY,    maxPixel*sizeof(olong));
  in->pixelFld  = ObitMemRealloc (in->pixelFld,  maxPixel*sizeof(gshort));
  in->pixelFlux = ObitMemRealloc (in->pixelFlux, maxPixel*sizeof(ofloat));

} /* end ObitDConCleanPxListResize */

/**
 * Update pixel list and window
 * \param in          The Pixel List object 
 * \param fields      Which fields? (1-rel) as zero terminated list
 * \param nSkip       Number of residuals to skip between ones accepted
 * \param minFluxLoad Minimum pixel flux density to accept
 * \param autoWinFlux min. residual flux allowed for auto Window
 * \param window      Windows object corresponding to Image Mosaic being CLEANED
 *                    Only pixels inside of the CLEAN window are used.
 * \param err Obit error stack object.
 */
void ObitDConCleanPxListUpdate (ObitDConCleanPxList *in, 
				olong *fields, olong nSkip,
				ofloat minFluxLoad,
				ofloat autoWinFlux,
				ObitDConCleanWindow *window, 
				ObitFArray *BeamPatch,
				ObitErr *err)
{
  ObitIOCode retCode;
  ObitImage *image=NULL;
  ObitIOSize IOsize = OBIT_IO_byPlane;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong  blc[IM_MAXDIM], trc[IM_MAXDIM];
  olong i, field,  ifld, number, excess, ix, iy, nx, ny, pos[2], skipCnt;
  ofloat *data, fblank =  ObitMagicF();
  gboolean blewIt=FALSE, *mask=NULL;
  gchar *routine = "ObitDConCleanPxListUpdate";

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
      
  /* Save Beam patch */
  in->BeamPatch = ObitFArrayUnref(in->BeamPatch);
  in->BeamPatch = ObitFArrayRef(BeamPatch);

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
    image = in->mosaic->images[field-1];
    
    /* Set output to full image, plane at a time */
    image->extBuffer = FALSE;
    dim[0] = IM_MAXDIM;
    blc[0] = blc[1] = 1;
    for (i=0; i<IM_MAXDIM-2; i++) blc[i+2] = in->plane[i];
    ObitInfoListPut (image->info, "BLC", OBIT_long, dim, blc, err); 
    trc[0] = trc[1] = 0;
    for (i=0; i<IM_MAXDIM-2; i++) trc[i+2] = in->plane[i];
    ObitInfoListPut (image->info, "TRC", OBIT_long, dim, trc, err); 
    dim[0] = 1;
    ObitInfoListPut (image->info, "IOBy", OBIT_long, dim, &IOsize, err);
    
    retCode = ObitImageOpen (image, OBIT_IO_ReadOnly, err);
    if (err->error) Obit_traceback_msg (err, routine, image->name);
    
    retCode = ObitImageRead (image, image->image->array, err);
    if (err->error) Obit_traceback_msg (err, routine, image->name);
    
    /* pointer to data */
    pos[0] = pos[1] = 0;
    data = ObitFArrayIndex(image->image, pos);
    
    /* Loop over image saving selected values */
    nx = image->myDesc->inaxes[0];
    ny = image->myDesc->inaxes[1];
    for (iy=0; iy<ny; iy++) {
      /* Get window mask */
      if (ObitDConCleanWindowRow(window, field, iy+1, &mask, err)) {
	for (ix=0; ix<nx; ix++) {
	  /* Want this one? */
	  if (mask[ix] && (data[ix]!=fblank) && 
	      (fabs(data[ix])>=minFluxLoad)) {
	    if (skipCnt>=nSkip) { /* take this one */
	      skipCnt = 0;
	      if ((number+1)>=in->maxPixel) {
		blewIt = TRUE;   /* Blew arrays */
		excess++;
	      } else { /* OK */
		in->pixelX[number]    = ix;
		in->pixelY[number]    = iy;
		in->pixelFld[number]  = field;
		in->pixelFlux[number] = data[ix];
		number++;
	      }
	    } else { /* skip this one to make them fit */
	      skipCnt++;
	    }
	  }
        }
      }
      data += nx;
    }
    in->nPixel = number;   /* Save actual number */

    retCode = ObitImageClose (image, err);
    if (err->error) Obit_traceback_msg (err, routine, image->name);
    
    /* Free Image array */
    image->image = ObitFArrayUnref(image->image);
    
    /* Cleanup - next field may have different size */
    if ((mask) && (ObitMemValid (mask))) mask = ObitMemFree (mask);

    ifld++;
    field = fields[ifld];
  } /* end loop over fields */

  /* Give warning if blew arrays */
  if (blewIt) 
    Obit_log_error(err, OBIT_InfoWarn,"%s: Number of pixels exceed arrays by  %d",
		   routine, excess);

} /* end ObitDConCleanPxListUpdate */

/**
 * Iteratively perform BGC CLEAN on Pixel list
 * \param in    The Pixel list object 
 * \param err   Obit error stack object.
 * \return TRUE if hit limit of niter or min. flux density.
 */
gboolean ObitDConCleanPxListCLEAN (ObitDConCleanPxList *in, ObitErr *err)
{
  gboolean done = FALSE;
  olong iter, iresid, ipeak=0, field, iXres, iYres, beamPatch, iBeam, pos[2];
  olong i, lpatch, ix, iy, irow, lastField=-1;
  ofloat peak, minFlux, factor, CCmin, atlim, xfac=1.0, resmax, xflux;
  ofloat subval, *beam=NULL;
  gchar reason[51];
  ObitTableCCRow *CCRow = NULL;
  ObitImageDesc *desc = NULL;
  ObitIOCode retCode;
  gchar *routine = "ObitDConCleanPxListCLEAN";

  /* error checks */
  if (err->error) return done;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Check number of residuals - bail if none */
  if (in->nPixel<=0) {
    Obit_log_error(err, OBIT_InfoWarn,"%s NO Residuals to CLEAN in %s",
		   routine, in->name);
    return TRUE;
  }

   /* How many components already done? */
  iter = MAX (0, in->currentIter);

  /* Setup */
  lpatch = in->BeamPatch->naxis[0];
  beamPatch = (lpatch-1)/2;
  pos[0] = pos[1] = 0;
  beam = ObitFArrayIndex(in->BeamPatch, pos); /* Beam patch pointer */
  CCmin = 1.0e20;
  atlim = 0.0;
  resmax = -1.0e20;

  /* DEBUG
  in->minFluxLoad = 0.01;
  in->niter = 4; */

  /* Tell details */
  if (in->prtLv>1) {
    Obit_log_error(err, OBIT_InfoErr,"Beam patch = %d cells, min. residual = %g Jy",
		   2*beamPatch, in->minFluxLoad);
    Obit_log_error(err, OBIT_InfoErr," %d residuals loaded ",
		   in->nPixel);
  }
  /* CLEAN loop */
  while (!done) {
    
    /* Find peak abs residual */
    peak = -1.0;
    for (iresid=0; iresid<in->nPixel; iresid++) {
      xflux = fabs(in->pixelFlux[iresid]);
      if (xflux>peak) {
	peak = xflux;
	ipeak = iresid;
      }
    }
    
    /* Setup for subtraction */
    field   = in->pixelFld[ipeak];
    minFlux = in->minFlux[field-1];
    factor  = in->factor[field-1];
    xflux = in->pixelFlux[ipeak];
    subval = xflux * in->gain[field-1];
    iXres =  in->pixelX[ipeak];
    iYres =  in->pixelY[ipeak];
    CCmin = MIN (CCmin, peak);
    
    /* If first pass set up stopping criteria */
    if (resmax < 0.0) {
      resmax = MAX (fabs(xflux), 1.0e-10);
      xfac = pow ((in->minFluxLoad / resmax), in->factor[field-1]);
    }
    
    /* Keep statistics */
    in->iterField[field-1]++;
    in->fluxField[field-1] += subval;
    in->totalFlux += subval;
    /* DEBUG 
       fprintf (stderr,"%s field %d flux %f total %f pos %d  %d\n",
       routine,  field, xflux, in->totalFlux, iXres, iYres);*/
    
    /* Do subtraction */
    for (iresid=0; iresid<in->nPixel; iresid++) {
      /* Is this inside the Beam patch ? */
      ix = in->pixelX[iresid];
      iy = in->pixelY[iresid] ;
      if ((in->pixelFld[iresid]==field) &&
	  (abs(ix-iXres) <= beamPatch) && 
	  (abs(iy-iYres) <= beamPatch)) {
	/* Index in beam patch array */
	iBeam = (beamPatch + (iy - iYres)) * lpatch +
	  (beamPatch + (ix - iXres));
	in->pixelFlux[iresid] -= subval * beam[iBeam];
      }
    }

    /* Write component to List */
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
    /* What's in AIPS is a bit more complex and adds field offset from tangent */
    CCRow->DeltaX = (iXres - desc->crpix[0]+1)*desc->cdelt[0];
    CCRow->DeltaY = (iYres - desc->crpix[1]+1)*desc->cdelt[1];
    CCRow->Flux = subval;
    /* May need Gaussian components */
    if ((in->CCTable[field-1])->myDesc->dim[(in->CCTable[field-1])->parmsCol]<=0) {
      if (in->circGaus[field-1]>0.0) {
	CCRow->parms[0] = in->circGaus[field-1];
	CCRow->parms[1] = in->circGaus[field-1];
	CCRow->parms[2] = 0.0;
	CCRow->parms[3] = 1;  /* type 1 = Gaussian */
      } else { /* point */
	CCRow->parms[0] = 0.0;
	CCRow->parms[1] = 0.0;
	CCRow->parms[2] = 0.0;
	CCRow->parms[3] = 0;  /* type 0 = Point */
      }
    }

    irow = in->iterField[field-1];
    retCode = ObitTableCCWriteRow (in->CCTable[field-1], irow, CCRow, err);
    if ((retCode != OBIT_IO_OK) || (err->error)) 
      Obit_traceback_val (err, routine, in->name, done);
    
    /* Test various stopping conditions */ 
    iter++;  /* iteration count */
    /* Are we finished after this? */
    done = done || (iter>=in->niter) || (fabs(xflux)<minFlux);
    if (done) {  /* set completion reason string */
      if (iter>=in->niter)     g_snprintf (reason, 50, "Reached Iter. limit");
      if (fabs(xflux)<minFlux) g_snprintf (reason, 50, "Reached min Clean flux");
      break;
    } 

    /* Diverging? */
    if (fabs (xflux) > 2.0*CCmin) {
      /* Give warning */
      Obit_log_error(err, OBIT_InfoWarn,"%s: Clean has begun to diverge, Stopping",
		     routine);
      g_snprintf (reason, 50, "Solution diverging");
      break;
    }
  
    /* BGC tests to tell if we should quit now */
    if (fabs (xflux) < in->minFluxLoad * (1.0 + atlim)) {
      g_snprintf (reason, 50, "Reached minimum algorithm flux");
      break;  /* jump out of CLEAN loop */
    }
    atlim += xfac / (ofloat)iter;   /* update BGC stopping criterion */
   
    /* autoWindow tests to tell if we should quit now */
    if (fabs (xflux) < in->autoWinFlux) {
      g_snprintf (reason, 50, "Reached minimum autoWindow flux");
      break;  /* jump out of CLEAN loop */
    }
    
  } /* end CLEANing loop */
  

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
  
  /* Keep statistics */
  in->currentIter = iter;
  
  /* Tell about results */
  if (in->prtLv>1) {
    Obit_log_error(err, OBIT_InfoErr,"Clean stopped because: %s", reason);
    Obit_log_error(err, OBIT_InfoErr,"%s: Min. Flux density %f",
		   routine, xflux);
    if (in->nfield>1) /* Multiple fields? */
      Obit_log_error(err, OBIT_InfoErr,"Field %d has %d CCs with %g Jy",
		     field, in->iterField[field-1], 
		     in->fluxField[field-1]);
    
    Obit_log_error(err, OBIT_InfoErr,"Total CLEAN %d CCs with %g Jy",
		   in->currentIter, in->totalFlux);
  }
  /* Keep maximum abs residual */
  in->maxResid = fabs(xflux);

  return done;
} /* end ObitDConCleanPxListCLEAN */

/**
 * Tells results of CLEANing
 * \param in          The Pixel List object 
 * \param ncomps      Array of total number of components per field, 
 *                    same order as in ObitImageMosaic
 * \param err Obit error stack object.
 */
olong ObitDConCleanPxListResult (ObitDConCleanPxList *in, olong *ncomp,
				 ObitErr *err)
{
  olong out = 0;
  olong i;

  /* error checks */
  if (err->error) return out;
  g_assert (ObitIsA(in, &myClassInfo));

  out = in->currentIter;
  for (i=0; i<in->nfield; i++) ncomp[i] = in->iterField[i];

  return out;
} /* end ObitDConCleanPxListResult */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitDConCleanPxListClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitDConCleanPxListClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitDConCleanPxListClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitDConCleanPxListClassInfoDefFn (gpointer inClass)
{
  ObitDConCleanPxListClassInfo *theClass = (ObitDConCleanPxListClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitDConCleanPxListClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitDConCleanPxListClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitDConCleanPxListGetClass;
  theClass->newObit       = (newObitFP)newObitDConCleanPxList;
  theClass->ObitCopy      = (ObitCopyFP)ObitDConCleanPxListCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitDConCleanPxListClear;
  theClass->ObitInit      = (ObitInitFP)ObitDConCleanPxListInit;

} /* end ObitDConCleanPxListClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitDConCleanPxListInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  olong i;
  ObitDConCleanPxList *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->thread    = newObitThread();
  in->info      = newObitInfoList();
  in->mosaic    = NULL;
  in->window    = NULL;
  in->BeamPatch = NULL;
  in->pixelX    = NULL;
  in->pixelY    = NULL;
  in->pixelFld  = NULL;
  in->pixelFlux = NULL;
  in->fluxField = NULL;
  in->fluxField = NULL;
  in->circGaus  = NULL;
  in->CCver     = NULL;
  in->gain      = NULL;
  in->minFlux   = NULL;
  in->factor    = NULL;
  in->CCTable   = NULL;
  for (i=0; i<IM_MAXDIM-2; i++) in->plane[i] = 1;
  in->autoWinFlux  = -1.0e20;


} /* end ObitDConCleanPxListInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitDConCleanPxList* cast to an Obit*.
 */
void ObitDConCleanPxListClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  olong i;
  ObitDConCleanPxList *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  for (i=0; i<in->nfield; i++) 
    in->CCTable[i] = ObitTableCCUnref(in->CCTable[i]);
  in->thread    = ObitThreadUnref(in->thread);
  in->info      = ObitInfoListUnref(in->info);
  in->mosaic    = ObitImageMosaicUnref(in->mosaic);
  in->window    = ObitDConCleanWindowUnref(in->window);
  in->BeamPatch = ObitFArrayUnref(in->BeamPatch);
  if (in->CCTable)   in->CCTable  =  ObitMemFree (in->CCTable);
  if (in->pixelX)    in->pixelX   =  ObitMemFree (in->pixelX);
  if (in->pixelY)    in->pixelY   =  ObitMemFree (in->pixelY);
  if (in->pixelFld)  in->pixelFld =  ObitMemFree (in->pixelFld);
  if (in->pixelFlux) in->pixelFlux=  ObitMemFree (in->pixelFlux);
  if (in->fluxField) in->fluxField=  ObitMemFree (in->fluxField);
  if (in->circGaus)  in->circGaus =  ObitMemFree (in->circGaus);
  if (in->iterField) in->iterField=  ObitMemFree (in->iterField);
  if (in->CCver)     in->CCver    =  ObitMemFree (in->CCver);
  if (in->gain)      in->gain     =  ObitMemFree (in->gain);
  if (in->minFlux)   in->minFlux  =  ObitMemFree (in->minFlux);
  if (in->factor)    in->factor   =  ObitMemFree (in->factor);

  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitDConCleanPxListClear */


