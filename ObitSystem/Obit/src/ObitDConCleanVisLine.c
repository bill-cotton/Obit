/* $Id$  */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2013                                               */
/*;  Associated Universities, Inc. Washington DC, USA.                */
/*;                                                                   */
/*;  This program is free software; you can redistribute it and/or    */
/*;  modify it under the terms of the GNU General Public License as   */
/*;  published by the Free Software Foundation; either version 2 of   */
/*;  the License, or (at your op2tion) any later version.              */
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

#include "ObitImageMosaic.h"
#include "ObitImageUtil.h"
#include "ObitDConClean.h"
#include "ObitDConCleanVis.h"
#include "ObitDConCleanVisLine.h"
#include "ObitMem.h"
#include "ObitFFT.h"
#include "ObitTableUtil.h"
#include "ObitTableCCUtil.h"
#include "ObitSkyGeom.h"
#include "ObitImageWB.h"
#include "ObitSkyModelMF.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitDConCleanVisLine.c
 * ObitDConCleanVisLine class function definitions.
 * Visibility based CLEAN class for spectral line imaging.
 * This class is derived from the ObitDConCleanVis class.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitDConCleanVisLine";

/** Function to obtain parent ClassInfo - ObitDConClean */
static ObitGetClassFP ObitParentGetClass = ObitDConCleanVisGetClass;

/**
 * ClassInfo structure ObitDConCleanVisLineClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitDConCleanVisLineClassInfo myClassInfo = {FALSE};

/** Private: Set Class function pointers. */
static void ObitDConCleanVisLineClassInfoDefFn (gpointer inClass);

/*--------------- File Global Variables  ----------------*/


/*---------------Private structures----------------*/
/* Single channel threaded function argument */
typedef struct {
  /** CLEAN Object */
  ObitDConCleanVisLine *in;
  /** channel number 1-rel */
  olong      chann;
  /** First 1-rel component to subtract in next subtraction 
                 by nfield  */
  olong      *bcomp;
  /** Pixel array for intermediate CLEAN */
  ObitFArray** pixarray;
  /** Was operation finished? */
  gboolean   done;
  /** thread number, <0 -> no threading  */
  olong      ithread;
  /** Obit Thread object */
  ObitThread  *thread;
  /** Obit error stack object */
  ObitErr    *err;
} ChanFuncArg;


/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitDConCleanVisLineInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitDConCleanVisLineClear (gpointer in);

/** Private: (re)make residuals. */
static void  MakeResiduals (ObitDConCleanVis *in, olong *fields, 
			    gboolean doBeam, ObitErr *err);

/** Private: (re)make all residuals. */
static void  MakeAllResiduals (ObitDConCleanVis *in, ObitErr *err);

/** Private: reset sky model. */
static gboolean ResetSkyModel (ObitDConCleanVis *in, ObitErr *err);

/** Private: reset Pixel List. */
static void ResetPixelList (ObitDConCleanVis *in, ObitErr *err);

/* Select components to be subtracted*/
gboolean ObitDConCleanVisLineSelect(ObitDConClean *inn, ObitFArray **pixarray, 
				    ObitErr *err);

/** Private: Low accuracy subtract CLEAN model. */
static void SubNewCCs (ObitDConCleanVis *in, olong *newCC, 
		       ObitFArray **pixarray, ObitErr *err);

/** Private: Make Threaded channel args */
static olong MakeChanFuncArgs (ObitDConCleanVisLine *in, ObitThread *thread,
			       ObitErr *err, ChanFuncArg ***ThreadArgs);

/** Private: Delete Threaded Channel args */
static void KillChanFuncArgs (olong nargs, ChanFuncArg **ThreadArgs);

/** Private: Low accuracy subtract CLEAN model. */
static void SubNewCCs (ObitDConCleanVis *in, olong *newCC, 
		       ObitFArray **pixarray, ObitErr *err);

/** Private: Create/init PxList. */
static void NewPxList (ObitDConCleanVis *in, ObitErr *err);

/** Private: Create/init Pixarray. */
static ObitFArray** NewPxArray (ObitDConCleanVis *in, olong *startCC, ObitErr *err);

/** Private: Delete Pixarray. */
static ObitFArray** KillPxArray (ObitDConCleanVis *in, ObitFArray **pixarray);

/** Private: Delete BeamPatches. */
static void KillBeamPatches (ObitDConCleanVis *in);

/** Private: Find peak brightness. */
static void FindPeak (ObitDConCleanVis *in, ObitErr *err);

/** Private: Pick next field(s) and get Residual image(s) */
static gboolean PickNext2D(ObitDConCleanVis *in, ObitErr *err);
static gboolean PickNext3D(ObitDConCleanVis *in, ObitErr *err);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitDConCleanVisLine* newObitDConCleanVisLine (gchar* name)
{
  ObitDConCleanVisLine* out;
  /*gchar *routine = "newObitDConCleanVisLine";*/

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitDConCleanVisLineClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitDConCleanVisLine));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitDConCleanVisLineInit((gpointer)out);

 return out;
} /* end newObitDConCleanVisLine */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitDConCleanVisLineGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitDConCleanVisLineClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitDConCleanVisLineGetClass */

/**
 * Make a deep copy of an ObitDConCleanVisLine.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitDConCleanVisLine* ObitDConCleanVisLineCopy (ObitDConCleanVisLine *in, 
						ObitDConCleanVisLine *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  gboolean oldExist;
  gchar *outName;
  gchar *routine = "ObitDConCleanVisLineCopy";

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
    out = newObitDConCleanVisLine(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, out);

  /*  copy this class */
  out->nPar = in->nPar;
  return out;
} /* end ObitDConCleanVisLineCopy */

/**
 * Make a copy of a object but do not copy the actual data
 * This is useful to create an DConCleanVisLine similar to the input one.
 * \param in  The object to copy
 * \param out An existing object pointer for output, must be defined.
 * \param err Obit error stack object.
 */
void ObitDConCleanVisLineClone  (ObitDConCleanVisLine *in, ObitDConCleanVisLine *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  gchar *routine = "ObitDConCleanVisLineClone";

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
  out->nPar = in->nPar;
} /* end ObitDConCleanVisLineClone */

/**
 * Creates an ObitDConCleanVisLine for parallel spectral line processing 
 * Initializes threading array (chArgs)
 * \param name   An optional name for the object.
 * \param nPar   Number of parallel channels
 * \param nAvg   Number of channels to be averaged
 * \param uvdata from which to create object, should have all control
                 information defined on info member.
 * \param err    Obit error stack object.
 * \return the new object.
 */
ObitDConCleanVisLine* 
ObitDConCleanVisLineCreate (gchar* name, olong nPar,  olong nAvg, ObitUV *uvdata,  
			    ObitErr *err)
{
  ObitDConCleanVisLine* out=NULL;
  ObitImage *img=NULL;
  olong nfield, i, loChan, hiChan, nchan, nif;
  gboolean freqFirst;
  odouble loFreq, Freq, hiFreq;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  olong bif=1, bchan=1, *BIFSpec=NULL, *EIFSpec=NULL, *BChanSpec=NULL, *EChanSpec=NULL;
  gchar keyword[12];
  gchar *routine = "ObitDConCleanVisLineCreate";

 /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return out;
  g_assert (ObitUVIsA(uvdata));

  /* Create basic structure */
  out = newObitDConCleanVisLine (name);
  out->nPar = nPar;

  /* How many channels to average? */
  dim[0] = dim[1] = dim[2] = dim[3] = dim[4] = 1;
  ObitInfoListAlwaysPut(uvdata->info, "nChAvg", OBIT_long, dim, &nAvg);
 
  /* Create Image Mosaic */
  out->mosaic = ObitImageMosaicCreate(name, uvdata, err);
  /* Define images */
  ObitImageMosaicDefine (out->mosaic, uvdata, TRUE, err);
  if (err->error) Obit_traceback_val (err, routine, name, out);

  /* Set channel selection on images in mosaic */
  BIFSpec   = g_malloc0(nPar*sizeof(olong));
  EIFSpec   = g_malloc0(nPar*sizeof(olong));
  BChanSpec = g_malloc0(nPar*sizeof(olong));
  EChanSpec = g_malloc0(nPar*sizeof(olong));
  ObitInfoListGetTest(uvdata->info, "BIF",   &type, dim, &bif);
  ObitInfoListGetTest(uvdata->info, "BChan", &type, dim, &bchan);
  /* These should be 0-rel */
  bif--; bchan--;
  for (i=0; i<nPar; i++) {
    BIFSpec[i]   = bif;   /* Better be only one */
    EIFSpec[i]   = bif;   /* Better be only one */
    BChanSpec[i] = bchan+i*nAvg;
    EChanSpec[i] = bchan+(i+1)*nAvg - 1;
    EChanSpec[i] = MIN (EChanSpec[i], uvdata->myDesc->inaxes[uvdata->myDesc->jlocf]-1);
  }
  for (i=0; i<out->mosaic->numberImages; i++) {
    img = out->mosaic->images[i];
    dim[0] = dim[1] = dim[2] = dim[3] = dim[4] = 1;
    ObitInfoListAlwaysPut(img->info, "nSpec", OBIT_long, dim, &nPar);
    dim[0] = nPar;
    ObitInfoListAlwaysPut(img->info, "BIFSpec",   OBIT_long, dim, BIFSpec);
    ObitInfoListAlwaysPut(img->info, "EIFSpec",   OBIT_long, dim, EIFSpec);
    ObitInfoListAlwaysPut(img->info, "BChanSpec", OBIT_long, dim, BChanSpec);
    ObitInfoListAlwaysPut(img->info, "EChanSpec", OBIT_long, dim, EChanSpec);
    if (img->myBeam) {   /* Beam if given */
      img = (ObitImage*)img->myBeam;
      dim[0] = dim[1] = dim[2] = dim[3] = dim[4] = 1;
      ObitInfoListAlwaysPut(img->info, "nSpec", OBIT_long, dim, &nPar);
      dim[0] = nPar;
      ObitInfoListAlwaysPut(img->info, "BIFSpec",   OBIT_long, dim, BIFSpec);
      ObitInfoListAlwaysPut(img->info, "EIFSpec",   OBIT_long, dim, EIFSpec);
      ObitInfoListAlwaysPut(img->info, "BChanSpec", OBIT_long, dim, BChanSpec);
      ObitInfoListAlwaysPut(img->info, "EChanSpec", OBIT_long, dim, EChanSpec);
    }
  }
  /* Create UV imager */
  out->imager = (ObitUVImager*)ObitUVImagerCreate2("UVImager", uvdata, out->mosaic, err);
  if (err->error) Obit_traceback_val (err, routine, name, out);

  /* Create SkyModel object */
  if (nPar>1)
    out->skyModel = (ObitSkyModel*)ObitSkyModelMFCreate ("SkyModel", out->mosaic);
  else
    out->skyModel = (ObitSkyModel*)ObitSkyModelCreate ("SkyModel", out->mosaic);

   /* Copy control info to SkyModel */
  ObitInfoListCopyData(uvdata->info, out->skyModel->info);
  
  /* Window object */
  out->window = ObitDConCleanWindowCreate ("CleanWindow", out->mosaic, err);
  if (err->error) Obit_traceback_val (err, routine, name, out);

  /* Arrays per field - including those in parent classes */
  nfield =  out->mosaic->numberImages;
  out->nfield  = nfield;
  out->gain        = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean Loop gain");
  out->minFlux     = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean minFlux");
  out->factor      = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean factor");
  out->quality     = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean quality");
  out->cleanable   = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean cleanable");
  out->fresh       = ObitMemAlloc0Name(nfield*sizeof(gboolean),"Clean fresh");
  out->maxAbsRes   = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean max res");
  out->avgRes      = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean avg res");
  out->imgRMS      = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Image RMS");
  out->imgPeakRMS  = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Image Peak/RMS");
  out->beamPeakRMS = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Beam Peak/RMS");
  out->currentFields = ObitMemAlloc0Name((nfield+3)*sizeof(olong),"Current fields");
  for (i=0; i<nfield; i++) {
    out->maxAbsRes[i]   = -1.0;
    out->avgRes[i]      = -1.0;
    out->quality[i]     = -1.0;
    out->cleanable[i]   = -1.0;
    out->fresh[i]       = FALSE;
    out->imgRMS[i]      = -1.0;
    out->imgPeakRMS[i]  = -1.0;
    out->beamPeakRMS[i] = -1.0;
  }

  /* init array for threading */
  out->nArgs = MakeChanFuncArgs (out, uvdata->thread, err, (ChanFuncArg ***)&out->chArgs);
  
  /* Check that Nargs compatable with nPar */
  Obit_retval_if_fail ((out->nArgs >= out->nPar), err, out, 
		       "%s: Fewer threads arguments %d than parallel channels %d",
		       routine, out->nArgs,  out->nPar);

  /* Set frequency info on first image of mosaic */
  img = out->mosaic->images[0];
  ObitImageOpen(img, OBIT_IO_ReadWrite, err);
  if  (err->error) Obit_traceback_val (err, routine, img->name, out);
  dim[0] = dim[1] = dim[2] = dim[3] = dim[4] = 1;
  ObitInfoListAlwaysPut(img->myDesc->info, "NSPEC", OBIT_long, dim, &nPar);
  nchan = uvdata->myDesc->inaxes[uvdata->myDesc->jlocf];
  if (uvdata->myDesc->jlocif>=0) nchan = uvdata->myDesc->inaxes[uvdata->myDesc->jlocif];
  else nif = 1;

  /* Frequency or IF faster in data */
  freqFirst = (uvdata->myDesc->jlocf<uvdata->myDesc->jlocif) || (uvdata->myDesc->jlocif<0);
  for (i=0; i<nPar; i++) {
    /* Get frequencies */
    if (freqFirst) {
      loChan = BIFSpec[i]*nchan + BChanSpec[i];
      hiChan = EIFSpec[i]*nchan + EChanSpec[i];
    } else {
      loChan = BIFSpec[i] + BChanSpec[i]*nif;
      hiChan = EIFSpec[i] + EChanSpec[i]*nif;
    }
    loFreq = MIN(uvdata->myDesc->freqArr[loChan], uvdata->myDesc->freqArr[hiChan]);
    hiFreq = MAX(uvdata->myDesc->freqArr[loChan], uvdata->myDesc->freqArr[hiChan]);
    Freq = 0.5 * (loFreq + hiFreq);
    sprintf (keyword, "FREQ%4.4d",i+1);
    ObitInfoListAlwaysPut(img->myDesc->info, keyword, OBIT_double, dim, &Freq);
    sprintf (keyword, "FREL%4.4d",i+1);
    ObitInfoListAlwaysPut(img->myDesc->info, keyword, OBIT_double, dim, &loFreq);
    sprintf (keyword, "FREH%4.4d",i+1);
    ObitInfoListAlwaysPut(img->myDesc->info, keyword, OBIT_double, dim, &hiFreq);
  }
  img->myStatus = OBIT_Modified;  /* Force update */
  ObitImageClose(img, err);
  if  (err->error) Obit_traceback_val (err, routine, img->name, out);

  /* Cleanup */
  g_free(BIFSpec);  g_free(EIFSpec); g_free(BChanSpec); g_free(EChanSpec);

  return out;
} /* end ObitDConCleanVisLineCreate */

/**
 * Creates an ObitDConCleanVisLine from optional components 
 * defined for convenience of derived classes 
 * \param name     An optional name for the object.
 * \param uvdata   from which to create object, should have all control
                   information defined on info member.
 * \param imager   Optional ObitUVImager to use, if NULL use default
 *                 Reference "stolen" (i.e. no need to Unref after call)
 * \param skyModel Optional ObitSkyModel to use, if NULL use default
 *                 Reference "stolen" (i.e. no need to Unref after call)
 * \param err      Obit error stack object.
 * \return the new object.
 */
ObitDConCleanVisLine* 
ObitDConCleanVisLineCreate2 (gchar* name, ObitUV *uvdata,  
			 ObitUVImager *imager, ObitSkyModel *skyModel, 
			 ObitErr *err)
{
  olong nfield, i;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitDConCleanVisLine* out=NULL;
  ofloat ftemp;
  gchar *routine = "ObitDConCleanVisLineCreate";

 /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return out;
  g_assert (ObitUVIsA(uvdata));

  /* Create basic structure */
  out = newObitDConCleanVisLine (name);

  /* Use or create UV imager and create its ImageMosaic */
  if (imager==NULL) {
    out->imager = ObitUVImagerCreate("UVImager", uvdata, err);
    if (err->error) Obit_traceback_val (err, routine, name, out);
  } else out->imager = ObitUVImagerRef(imager);

  /* Save uv Mosaic reference */
  out->mosaic = ObitUVImagerGetMosaic(out->imager, err);

  /*  Use or create SkyModel object */
  if (skyModel==NULL) out->skyModel = ObitSkyModelCreate ("SkyModel", out->mosaic);
  else out->skyModel = ObitSkyModelRef(skyModel);

   /* Copy control info to SkyModel */
  ObitInfoListCopyData(uvdata->info, out->skyModel->info);
  /* Disable any value of minFlux */
  ftemp = -1.0e20;
  dim[0] = 1;dim[1] = 1;
  ObitInfoListAlwaysPut (out->skyModel->info, "minFlux", OBIT_float, dim, &ftemp);
  
  /* Window object */
  out->window = ObitDConCleanWindowCreate ("CleanWindow", out->mosaic, err);
  if (err->error) Obit_traceback_val (err, routine, name, out);

  /* Arrays per field - including those in parent classes */
  nfield =  out->mosaic->numberImages;
  out->nfield  = nfield;
  out->gain        = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean Loop gain");
  out->minFlux     = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean minFlux");
  out->factor      = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean factor");
  out->quality     = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean quality");
  out->cleanable   = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean cleanable");
  out->fresh       = ObitMemAlloc0Name(nfield*sizeof(gboolean),"Clean fresh");
  out->maxAbsRes   = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean max res");
  out->avgRes      = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean avg res");
  out->imgRMS      = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Image RMS");
  out->imgPeakRMS  = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Image Peak/RMS");
  out->beamPeakRMS = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Beam Peak/RMS");
  out->currentFields = ObitMemAlloc0Name((nfield+3)*sizeof(olong),"Current fields");
  for (i=0; i<nfield; i++) {
    out->maxAbsRes[i]   = -1.0;
    out->avgRes[i]      = -1.0;
    out->quality[i]     = -1.0;
    out->cleanable[i]   = -1.0;
    out->fresh[i]       = FALSE;
    out->imgRMS[i]      = -1.0;
    out->imgPeakRMS[i]  = -1.0;
    out->beamPeakRMS[i] = -1.0;
 }

  /* init array for threading */
  out->nArgs = MakeChanFuncArgs (out, uvdata->thread, err, (ChanFuncArg ***)&out->chArgs);
  
  /* Check that Nargs compatable with nPar */
  Obit_retval_if_fail ((out->nArgs >= out->nPar), err, out, 
		       "%s: Fewer threads arguments %d than parallel channels %d",
		       routine, out->nArgs,  out->nPar);

  return out;
} /* end ObitDConCleanVisLineCreate2 */

/**
 * Update channel selection 
 * It's unclear how well this works
 * \param in     CLEAN to update
 * \param nPar   Number of parallel channels
 * \param BIF    IF to select
 * \param BChan  First Channel to select
 */
void ObitDConCleanVisLineUpdate (ObitDConCleanVisLine *in, olong nPar, 
				 olong BIF, olong BChan)
{
  ObitImage *img=NULL;
  olong i;
  gint32 dim[MAXINFOELEMDIM];
  olong bif=1, bchan=1, *BIFSpec=NULL, *EIFSpec=NULL, *BChanSpec=NULL, *EChanSpec=NULL;

 /* error checks */
  g_assert (ObitDConCleanVisLineIsA(in));

  /* Set channel selection on images in mosaic */
  BIFSpec   = g_malloc0(nPar*sizeof(olong));
  EIFSpec   = g_malloc0(nPar*sizeof(olong));
  BChanSpec = g_malloc0(nPar*sizeof(olong));
  EChanSpec = g_malloc0(nPar*sizeof(olong));
  bif   = BIF;
  bchan = BChan;
  /* These should be 0-rel */
  bif--; bchan--;
  for (i=0; i<nPar; i++) {
    BIFSpec[i]   = bif;   /* Better be only one */
    EIFSpec[i]   = bif;   /* Better be only one */
    BChanSpec[i] = bchan+i*nPar;
    EChanSpec[i] = bchan+(i+1)*nPar - 1;
  }
  for (i=0; i<in->mosaic->numberImages; i++) {
    img    = in->mosaic->images[i];
    dim[0] = dim[1] = dim[2] = dim[3] = dim[4] = 1;
    ObitInfoListAlwaysPut(img->info, "nSpec", OBIT_long, dim, &nPar);
    dim[0] = nPar;
    ObitInfoListAlwaysPut(img->info, "BIFSpec",   OBIT_long, dim, BIFSpec);
    ObitInfoListAlwaysPut(img->info, "EIFSpec",   OBIT_long, dim, EIFSpec);
    ObitInfoListAlwaysPut(img->info, "BChanSpec", OBIT_long, dim, BChanSpec);
    ObitInfoListAlwaysPut(img->info, "EChanSpec", OBIT_long, dim, EChanSpec);
    if (img->myBeam) {   /* Beam if given */
      img = (ObitImage*)img->myBeam;
      dim[0] = dim[1] = dim[2] = dim[3] = dim[4] = 1;
      ObitInfoListAlwaysPut(img->info, "nSpec", OBIT_long, dim, &nPar);
      dim[0] = nPar;
      ObitInfoListAlwaysPut(img->info, "BIFSpec",   OBIT_long, dim, BIFSpec);
      ObitInfoListAlwaysPut(img->info, "EIFSpec",   OBIT_long, dim, EIFSpec);
      ObitInfoListAlwaysPut(img->info, "BChanSpec", OBIT_long, dim, BChanSpec);
      ObitInfoListAlwaysPut(img->info, "EChanSpec", OBIT_long, dim, EChanSpec);
    }
  }
  /* Cleanup */
  g_free(BIFSpec);  g_free(EIFSpec); g_free(BChanSpec); g_free(EChanSpec);

} /* end ObitDConCleanVisLineUpdate */

/**
 * Set default CLEAN windows in mosaic
 * Loops over channels.
 * If mosaic member  Radius>0 then make round boxes on Fly's eye field
 * with this radius, else use rectangular box including all but outer 5 pixels
 * On outlier fields, use rectangular box of width OutlierSize.
 * If CLEANBox defined in in->info then its contents are used for field 1.
 * Assumes all images in mosaic have descriptors defined.
 * If there are multiple tapers, any window on facet 1 are copied to
 * the corresponding facets at other resolutions.
 * Uses base class function.
 * \param in   The CLEAN object
 * \param err Obit error stack object.
 */
void ObitDConCleanVisLineDefWindow(ObitDConClean *inn, ObitErr *err)
{
  ObitDConCleanVisLine *in;
  ChanFuncArg *inArr=NULL;
  const ObitDConCleanVisClassInfo *ParentClass = myClassInfo.ParentClass;
  olong ichan;
  /*gchar *routine = "ObitDConCleanVisLineDefWindow";*/

  /* Cast input to this type */
  in = (ObitDConCleanVisLine*)inn;

  /* error checks */
  if (err->error) return;
  g_assert (ObitDConCleanVisLineIsA(in));

  /* Set on main */
  ParentClass->ObitDConCleanDefWindow(inn, err);

  /* Loop over channels */
  for (ichan=0; ichan<in->nArgs; ichan++) {
    inArr = ((ChanFuncArg*)in->chArgs[ichan]);

    inArr->done = FALSE;  /* initialize inArr */

    /* Call parent for operation */
    ParentClass->ObitDConCleanDefWindow((ObitDConClean*)inArr->in, err);
  } /* End channel loop */    
} /* end ObitDConCleanVisLineDefWindow */

/**
 * Subtract components from uv data for all channels
 * Generate a fake MFImageTSpec CC table with all the CCs of the planes 
 * concatenated.
 * \param in   The CLEAN object
 * \param err Obit error stack object.
 */
void ObitDConCleanVisLineSub(ObitDConClean *inn, ObitErr *err)
{
  ObitDConCleanVisLine *in;
  const ObitDConCleanVisClassInfo *ParentClass = myClassInfo.ParentClass;
  olong ifld, ichan,  bcomp, ecomp, iCCVer, nCCVer, oCCVer;
  olong *bcopy=NULL, *bc=NULL, *ec=NULL, *cv=NULL;
  gboolean doGaus;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ChanFuncArg *inArr=NULL;
  gchar *routine = "ObitDConCleanVisLineSub";

  /* Cast input to this type */
  in = (ObitDConCleanVisLine*)inn;

  /* error checks */
  if (err->error) return;
  g_assert (ObitDConCleanVisLineIsA(in));

  /* Collect new components into a TSpec CC table */
  bcopy = g_malloc0(in->nArgs*sizeof(olong));
  bc    = g_malloc0(in->nfield*sizeof(olong));
  ec    = g_malloc0(in->nfield*sizeof(olong));
  cv    = g_malloc0(in->nfield*sizeof(olong));
  /* Loop over fields building TSpec CC Tables */
  for (ifld=0; ifld<in->nfield; ifld++) {
    /* Get number of components to use from each */
    for (ichan=0; ichan<in->nPar; ichan++) {
      inArr = ((ChanFuncArg*)in->chArgs[ichan]);
      bcopy[ichan] = inArr->bcomp[ifld];  /* How many already done? */
      inArr->bcomp[ifld] = inArr->in->Pixels->iterField[ifld];  /* for next time */
    }
      
    /* Make combined TSpec CC Table if more than 1 channel*/
    if (in->nPar>1) {
      iCCVer = 1;
      nCCVer = in->nPar;
      oCCVer = iCCVer + in->nArgs;
      doGaus = in->maxBeamTaper>0.0;
      ObitTableCCUtilCombTSpec (in->mosaic->images[ifld], iCCVer, nCCVer, oCCVer,
				bcopy, &bcomp, &ecomp, doGaus, err);
      cv[ifld] = oCCVer;
      bc[ifld] = bcomp;
      ec[ifld] = ecomp;
      in->skyModel->startComp[ifld] = bcomp;
      in->skyModel->endComp[ifld]   = ecomp;
      /* Lower routines too clever */
      in->Pixels->iterField[ifld] = ecomp;
      /* Number of spectral planes */
      dim[0] = dim[1] = dim[2] = dim[3] = 1;
      ObitInfoListAlwaysPut (in->mosaic->images[ifld]->info, "NSPEC", OBIT_long, dim, &in->nPar);
    } else {  /* 1 channel */
      iCCVer = 1;
      in->skyModel->CCver[ifld] = iCCVer;
      oCCVer = iCCVer;
      inArr = ((ChanFuncArg*)in->chArgs[0]);
      cv[ifld] = oCCVer;
      bc[ifld] = inArr->in->skyModel->startComp[ifld];
      ec[ifld] = inArr->in->skyModel->endComp[ifld];
      in->skyModel->startComp[ifld] = inArr->in->skyModel->startComp[ifld];
      in->skyModel->endComp[ifld]   = inArr->in->skyModel->endComp[ifld];
      in->Pixels->iterField[ifld]   = inArr->in->Pixels->iterField[ifld];
   }
  } /* end loop over fields */

  /* Set components on Sky Model */
  dim[0] = in->nfield; dim[1] = dim[2] = dim[3] = 1;
  ObitInfoListAlwaysPut (in->skyModel->info, "CCVer", OBIT_long, dim, cv);
  ObitInfoListAlwaysPut (in->skyModel->info, "BComp", OBIT_long, dim, bc);
  ObitInfoListAlwaysPut (in->skyModel->info, "EComp", OBIT_long, dim, ec);

  in->CCver = oCCVer;

  /* Most work in parent class */
  ParentClass->ObitDConCleanSub(inn, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Mark everything as unfresh */
    for (ichan=0; ichan<in->nPar; ichan++) {
      inArr = ((ChanFuncArg*)in->chArgs[ichan]);
      for (ifld=0; ifld<in->nfield; ifld++)inArr->in->fresh[ifld] = FALSE;
    }


  if (bcopy) g_free(bcopy);
  if (bc)    g_free(bc);
  if (ec)    g_free(ec);

} /* end ObitDConCleanVisLineSub */

/**
 * \param in        The CLEAN object.
 * \param pixarray  If NonNULL use instead of the flux densities from the image file.
 *                  This is an array of ObitFarrays corresponding to fields in
                    in->currentFields
 * \param err       Obit error stack object.
 * \return TRUE if some planes had good data.
 */
gboolean ObitDConCleanVisLineValid(ObitDConCleanVisLine *in)
{
  gboolean allBlank, anyGood;
  olong ipln, npln;
  ObitDConCleanVisLine *tmpClean;

  allBlank = TRUE;
  npln = in->nArgs;
  for (ipln=0; ipln<npln; ipln++) {
    tmpClean = ((ChanFuncArg*)in->chArgs[ipln])->in;
    allBlank = allBlank && (tmpClean->maxAbsRes[0]<=0.0);
  }
  anyGood = !allBlank;
  return anyGood;
} /* end  ObitDConCleanVisLineValid */

/**
 * Get image and beam statistics and prepare to deconvolve
 * If autoWindow option is selected then windows may be added 
 * \param in        The object to deconvolve
 * \param pixarray  If NonNULL use instead of the flux densities from the image file.
 *                  This is an array of ObitFarrays corresponding to fields in
                    in->currentFields
 * \param err       Obit error stack object.
 * \return TRUE if CLEAN should continue
 */
gboolean ObitDConCleanVisLinePixelStats(ObitDConClean *inn, ObitFArray **pixarray, 
					ObitErr *err)
{
  const ObitDConCleanVisClassInfo *ParentClass = myClassInfo.ParentClass;
  ObitDConCleanVisLine *tin=(ObitDConCleanVisLine*)inn, *in;
  gboolean t, newWin=FALSE;
  olong ichan, ifld, best;
  ofloat bestQuality;
  ChanFuncArg *inArr=NULL;
  gchar *routine = "ObitDConCleanVisLinePixelStats";

   /* May need to use upper level Clean for this, use motherShip if inn is a lower level */
  if (tin->chArgs==NULL) in = (ObitDConCleanVisLine*)tin->motherShip;
  else                   in = tin;

  /* error checks */
  if (err->error) return newWin;
  g_assert (ObitDConCleanIsA(in));

  /* Find channel with highest quality field to copy values to main */
  best          = 0;
  bestQuality = 0.0;

  /* Loop over channels */
  for (ichan=0; ichan<in->nPar; ichan++) {
    inArr = ((ChanFuncArg*)in->chArgs[ichan]);
    /* Channel work in parent */
    t = ParentClass->ObitDConCleanPixelStats((ObitDConClean*)inArr->in, 
					     inArr->pixarray, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, newWin);
    newWin = newWin || t;

    /* find best quality */
    for (ifld=0; ifld<in->nfield; ifld++) {
      if (inArr->in->currentFields[ifld] <= 0) break;
      if (inArr->in->quality[inArr->in->currentFields[ifld]-1]>bestQuality) {
	best        = ichan;
	bestQuality = inArr->in->quality[inArr->in->currentFields[ifld]-1];
      }
    }
  } /* end channel loop */
  
  /* Copy info for best channel to main */
  inArr = ((ChanFuncArg*)in->chArgs[best]);
    
  in->peakFlux = inArr->in->peakFlux;
  for (ifld=0; ifld<in->nfield; ifld++) {
    in->currentFields[ifld] = inArr->in->currentFields[ifld];
    if (inArr->in->currentFields[ifld] <= 0) break;
    in->quality[in->currentFields[ifld]-1] = 
      inArr->in->quality[in->currentFields[ifld]-1];
    in->cleanable[in->currentFields[ifld]-1] = 
      inArr->in->cleanable[in->currentFields[ifld]-1];
    in->maxAbsRes[in->currentFields[ifld]-1] = 
      inArr->in->maxAbsRes[in->currentFields[ifld]-1];
    in->minFlux[in->currentFields[ifld]-1] = 
      inArr->in->minFlux[in->currentFields[ifld]-1];
    in->fresh[in->currentFields[ifld]-1] = 
      inArr->in->fresh[in->currentFields[ifld]-1];
    in->beamPeakRMS[in->currentFields[ifld]-1] = 
      inArr->in->beamPeakRMS[in->currentFields[ifld]-1];
    in->imgPeakRMS[in->currentFields[ifld]-1] = 
      inArr->in->imgPeakRMS[in->currentFields[ifld]-1];
    in->avgRes[in->currentFields[ifld]-1] = 
      inArr->in->avgRes[in->currentFields[ifld]-1];
    in->imgRMS[in->currentFields[ifld]-1] = 
      inArr->in->imgRMS[in->currentFields[ifld]-1];
  } /* end loop over fields */
    
  return newWin;
} /* end ObitDConCleanVisLinePixelStats */

/**
 * Get Image statistics to help decide which field is next to process
 * Multiple channel version, some initialization of parameters.
 * The outer window is used to specify valid pixels, 
 * For this version the following are calculated:
 * \li maxAbsRes   Maximum absolute windowed residual value (doBEAM=FALSE)
 * \li avgRes      Average windowed residual value (doBEAM=FALSE)
 * \li imgRMS      Image RMS  (doBeam=FALSE)
 * \li imgPeakRMS  Image Peak/RMS  (doBeam=FALSE)
 * \li beamPeakRMS Beam Peak/RMS  (doBeam = TRUE)
 * \param in    The object to deconvolve
 * \param field Which field? (1-rel) <=0 -> all;
 * \param doBeam If TRUE, do Beam statistics else Image
 * \param err   Obit error stack object.
 */
void ObitDConCleanVisLineImageStats(ObitDConClean *inn, olong field, gboolean doBeam, 
				    ObitErr *err)
{
  const ObitDConCleanVisClassInfo *ParentClass = myClassInfo.ParentClass;
  ObitDConCleanVisLine *in=(ObitDConCleanVisLine*)inn;
  olong ichan, ifld, best, hiFld, loFld, i;
  ofloat bestQuality;
  ChanFuncArg *inArr=NULL;
  gchar *routine = "ObitDConCleanVisLinePixelStats";

  /* error checks */
  if (err->error) return;
  g_assert (ObitDConCleanIsA(in));

  /* Find channel with highest quality field to copy values to main */
  best          = 0;
  bestQuality = 0.0;

  /* Field range */
  if (field<=0) {
    loFld = 0;
    hiFld = in->nfield - 1;
  } else {
    loFld = field - 1;
    hiFld = field - 1;
  }

  /* Loop over channels */
  for (ichan=0; ichan<in->nPar; ichan++) {
    inArr = ((ChanFuncArg*)in->chArgs[ichan]);

    /* Copy control parameters of first field */
    if (field<=1) {
      inArr->in->niter        = in->niter;
      inArr->in->Pixels->niter= in->niter;
      inArr->in->minPatchSize = in->minPatchSize;
      inArr->in->minFluxLoad  = in->minFluxLoad;
      inArr->in->maxPixel     = in->maxPixel;
      inArr->in->bmaj         = in->bmaj;
      inArr->in->bmin         = in->bmin;
      inArr->in->bpa          = in->bpa;
      inArr->in->nfield       = in->nfield ;
      inArr->in->modelMode    = in->modelMode;
      inArr->in->doRestore    = in->doRestore;
      inArr->in->doXRestore   = in->doXRestore;
      inArr->in->doFlatten    = in->doFlatten;
      inArr->in->doWeight     = in->doWeight;
      inArr->in->doRecenter   = in->doRecenter;
      inArr->in->autoWinFlux  = in->autoWinFlux;
      inArr->in->ccfLim       = in->ccfLim;
      inArr->in->SDIGain      = in->SDIGain;
      inArr->in->doSDI        = in->doSDI;
      inArr->in->autoWindow   = in->autoWindow;
      inArr->in->reuseFlux    = in->reuseFlux;
      inArr->in->autoCen      = in->autoCen;
      inArr->in->maxBeamTaper = in->maxBeamTaper;
      inArr->in->minBeamTaper = in->minBeamTaper;
      for (i=0; i<10; i++) inArr->in->MResKnob[i]  = in->MResKnob[i];
    }
    for (ifld=loFld; ifld<hiFld; ifld++) {
      inArr->in->gain[ifld]    = in->gain[ifld];
      inArr->in->minFlux[ifld] = in->minFlux[ifld];
      inArr->in->factor[ifld]  = in->factor[ifld];
    }

    /* Set channel */
    inArr->in->plane[0] = ichan+1;
    inArr->in->nfield   = in->nfield;

    /* Channel work in parent */
    ParentClass->ObitDConCleanImageStats((ObitDConClean*)inArr->in, 
					 field, doBeam, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);

    /* find best quality */
    for (ifld=loFld; ifld<hiFld; ifld++) {
      if (inArr->in->currentFields[ifld] <= 0) break;
      if (inArr->in->quality[ifld]>bestQuality) {
	best        = ichan;
	bestQuality = inArr->in->quality[ifld];
      }
    }
  } /* end channel loop */
  
  /* Copy stats info for best channel to main */
  inArr = ((ChanFuncArg*)in->chArgs[best]);

  for (ifld=loFld; ifld<hiFld; ifld++) {
      if (doBeam) in->beamPeakRMS[ifld] = inArr->in->beamPeakRMS[ifld];
      else {
	in->maxAbsRes[ifld]  = inArr->in->maxAbsRes[ifld];
	in->avgRes[ifld]     = inArr->in->avgRes[ifld];
	in->imgRMS[ifld]     = inArr->in->imgRMS[ifld];
	in->imgPeakRMS[ifld] = inArr->in->imgPeakRMS[ifld];
      }
    } /* end loop over fields */
    
} /* end ObitDConCleanVisLineImageStats */

/**
 * Automatically set a box in window if needed
 * Loops over channels, sets pre CL`uEan no. components
 * autoWindow feature will automatically set CLEAN windows inside 
 * a predefined outer window.  Each cycle the residuals inside the outer 
 * window are searched to the maximum value; if the peak is outside the 
 * inner window and > 5 sigma, a new round box  is added 
 * to the window.  Cleaning in each cycle will stop when the peak residual 
 * drops to the level of the highest value outside the CLEAN window.
 * Facets might have a window added iff 1) it is the facet with the highest 
 * in window pixel, or 2) it has a cleanable, out of window pixel within 30% of the peak.
 * \param in       The object to restore
 * \param fields   Field numbers (1-rel) in ImageMosaic
 *                 zero terminated list, no longer than in->numCurrentField
 * \param pixarray If nonNULL, then use this array of pixels rather than in the ImageMosaic
 *                 Elements are Pixel arrays corresponding to fields in fields [only 1]
 * \param err      Obit error stack object.
 * \return TRUE is a new window added
 */
gboolean ObitDConCleanVisLineAutoWindow(ObitDConClean *inn, olong *fields, ObitFArray **pixarray, 
					ObitErr *err)
{
  gboolean newWin = FALSE, t;
  const ObitDConCleanVisClassInfo *ParentClass = myClassInfo.ParentClass;
  ObitDConCleanVisLine *tin=(ObitDConCleanVisLine*)inn, *in;
  olong ichan, ifld;
  ChanFuncArg *inArr=NULL;
  gchar *routine = "ObitDConCleanVisLinAutoWindowe";
  
  /* May need to use upper level Clean for this, use motherShip if inn is a lower level */
  if (tin->chArgs==NULL) in = (ObitDConCleanVisLine*)tin->motherShip;
  else                   in = tin;

  /* error checks */
  if (err->error) return newWin;
  g_assert (ObitDConCleanIsA(in));

  /* Loop over channels */
  for (ichan=0; ichan<in->nPar; ichan++) {
    inArr = ((ChanFuncArg*)in->chArgs[ichan]);
    inArr->in->window = ObitDConCleanWindowRef(in->window);
    /* Channel work in parent */
    t = ParentClass->ObitDConCleanAutoWindow((ObitDConClean*)inArr->in, 
					     inArr->in->currentFields,
					     inArr->pixarray, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, newWin);
   newWin = newWin || t;
  } /* end channel loop */

  /* save max. number of CLEAN components on main */
  for (ifld=0; ifld<in->nfield; ifld++) in->Pixels->iterField[ifld] = 0;
  for (ifld=0; ifld<in->nfield; ifld++) {
    for (ichan=0; ichan<in->nPar; ichan++) {
     inArr = ((ChanFuncArg*)in->chArgs[ichan]);
     in->Pixels->iterField[ifld] = MAX(in->Pixels->iterField[ifld], 
				       inArr->in->Pixels->iterField[ifld]);
    }
  } /* end field loop */
  
  return newWin;
} /* end ObitDConCleanVisLineAutoWindow */

/**
 * Select components to be subtracted, loops over channels
 * \param in   The object to deconvolve
 * \param pixarray   If NonNULL use instead of the flux densities from the image file.
 *                   Array of ObitFArrays corresponding to fields in in->currentFields 
 * \param err        Obit error stack object.
 * \return TRUE if deconvolution is complete
 */
gboolean ObitDConCleanVisLineSelect(ObitDConClean *inn, ObitFArray **pixarray, 
				    ObitErr *err)
{
  gboolean done=TRUE, t;
  const ObitDConCleanVisClassInfo *ParentClass = myClassInfo.ParentClass;
  ObitDConCleanVisLine *in=(ObitDConCleanVisLine*)inn;
  olong ichan;
  ChanFuncArg *inArr=NULL;
  gchar *routine = "ObitDConCleanSelectLineVis";

  /* error checks */
  if (err->error) return done;
  g_assert (ObitDConCleanIsA(in));

  in->Pixels->complCode = OBIT_CompReasonUnknown;
  in->Pixels->currentIter = 0;
  /* Loop over channels */
  for (ichan=0; ichan<in->nPar; ichan++) {
    inArr = ((ChanFuncArg*)in->chArgs[ichan]);
    ((ObitDConClean*)inArr->in)->CCver = ichan+1;
    /* Channel work in parent */
    if (!inArr->done) {
     if (err->prtLv>=2) Obit_log_error(err, OBIT_InfoErr, 
					"Clean chan %d of %d",
					ichan+1, in->nPar);
     
      t = ParentClass->ObitDConCleanSelect((ObitDConClean*)inArr->in, 
					   inArr->pixarray, err);
    } else t = TRUE;
    done = done && t;
    if (err->error) Obit_traceback_val (err, routine, in->name, done);
    in->Pixels->complCode = MAX (in->Pixels->complCode, inArr->in->Pixels->complCode);
    inArr->done = inArr->done || t;
    /* AutoWin trumps all */
    if ((inArr->in->Pixels->complCode==OBIT_CompReasonAutoWin) || 
	(in->Pixels->complCode==OBIT_CompReasonAutoWin))
      in->Pixels->complCode = OBIT_CompReasonAutoWin;
    /* Maximum number of components */
    in->Pixels->currentIter = MAX (in->Pixels->currentIter, inArr->in->Pixels->currentIter);
  } /* end loop over channels */
 
  return done;
} /* end ObitDConCleanVisLineSelect */

/**
 * Restore components removed from the residual image(s)
 * Loop over channels
 * \param in   The object to restore
 * \param err Obit error stack object.
 */
void ObitDConCleanVisLineRestore(ObitDConClean *inn, ObitErr *err)
{
  const ObitDConCleanVisClassInfo *ParentClass = myClassInfo.ParentClass;
  ObitDConCleanVisLine *in=(ObitDConCleanVisLine*)inn;
  olong ichan;
  ChanFuncArg *inArr=NULL;
  gchar *routine = "ObitDConCleanVisLineRestore";

  /* error checks */
  if (err->error) return;
  g_assert (ObitDConCleanVisLineIsA(in));

  /* Loop over channels */
  for (ichan=0; ichan<in->nPar; ichan++) {
    inArr = ((ChanFuncArg*)in->chArgs[ichan]);

    /* Set channel */
    inArr->in->plane[0] = ichan+1;
    inArr->in->nfield   = in->nfield;
    inArr->in->CCver    = ichan+1;

    /* Channel work in parent */
    ParentClass->ObitDConCleanRestore((ObitDConClean*)inArr->in, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
  } /* end channel loop */

} /* end ObitDConCleanVisLineRestore */

/**
 * Restore components removed from one field but also 
 * appearing in another.  Does brute force convolution.
 * Adopted from the AIPSish QOOP:QCLEAN.FOR(CLOVER)
 * Presumes in->mosaic and image descriptors filled in.
 * \param in   The object to restore
 * \param err Obit error stack object.
 */
void ObitDConCleanVisLineXRestore(ObitDConClean *inn, ObitErr *err)
{
  const ObitDConCleanVisClassInfo *ParentClass = myClassInfo.ParentClass;
  ObitDConCleanVisLine *in=(ObitDConCleanVisLine*)inn;
  ChanFuncArg *inArr=NULL;
  olong ichan;
  gchar *routine = "ObitDConCleanVisLineXRestore";

  /* error checks */
  if (err->error) return;
  g_assert (ObitDConCleanVisLineIsA(in));

  /* Loop over channels */
  for (ichan=0; ichan<in->nPar; ichan++) {
    inArr = ((ChanFuncArg*)in->chArgs[ichan]);

    /* Set channel */
    inArr->in->plane[0] = ichan+1;
    inArr->in->nfield   = in->nfield;
    inArr->in->CCver    = ichan+1;

    /* Channel work in parent */
    ParentClass->ObitDConCleanXRestore((ObitDConClean*)inArr->in, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
  } /* end channel loop */

} /* end ObitDConCleanVisLineXRestore */

/**
 * Flatten multiple facets if needed
 * Does Flatten if FullField member of mosaic member is defined.
 * \param in   The object to deconvolve
 * \param err Obit error stack object.
 */
void ObitDConCleanVisLineFlatten(ObitDConClean *inn, ObitErr *err)
{
  const ObitDConCleanVisClassInfo *ParentClass = myClassInfo.ParentClass;
  ObitDConCleanVisLine *in=(ObitDConCleanVisLine*)inn;
  ChanFuncArg *inArr=NULL;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong ichan;
  gchar *routine = "ObitDConCleanVisLineFlatten";

  /* error checks */
  if (err->error) return;
  g_assert (ObitDConCleanVisLineIsA(in));

  /* Bail if mosaic has no Full Field - nothing to do */
  if (in->mosaic->FullField==NULL) return;

  /* Loop over channels */
  for (ichan=0; ichan<in->nPar; ichan++) {
    inArr = ((ChanFuncArg*)in->chArgs[ichan]);

    /* Set Full field on channel clean */
    if (inArr->in->mosaic->FullField==NULL) 
      inArr->in->mosaic->FullField = ObitImageRef(in->mosaic->FullField);

    /* Set channel */
    inArr->in->plane[0] = ichan+1;
    inArr->in->nfield   = in->nfield;
    ObitInfoListAlwaysPut(inArr->in->mosaic->info, "planeNo", OBIT_long, dim, 
			  &inArr->in->plane[0]);

    /* Channel work in parent */
    ParentClass->ObitDConCleanFlatten((ObitDConClean*)inArr->in, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
  } /* end channel loop */

}/* end ObitDConCleanVisLineFlatten */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitDConCleanVisLineClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitDConCleanVisLineClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitDConCleanVisLineClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitDConCleanVisLineClassInfoDefFn (gpointer inClass)
{
  ObitDConCleanVisLineClassInfo *theClass = (ObitDConCleanVisLineClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitDConCleanVisLineClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitDConCleanVisLineClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitDConCleanVisLineGetClass;
  theClass->newObit       = (newObitFP)newObitDConCleanVisLine;
  theClass->ObitCopy      = (ObitCopyFP)ObitDConCleanVisLineCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitDConCleanVisLineClear;
  theClass->ObitInit      = (ObitInitFP)ObitDConCleanVisLineInit;
  theClass->ObitDConCleanSub = (ObitDConCleanSubFP)ObitDConCleanVisLineSub;
  theClass->ObitDConCleanPixelStats = (ObitDConCleanPixelStatsFP)ObitDConCleanVisLinePixelStats;
  theClass->ObitDConCleanImageStats = (ObitDConCleanImageStatsFP)ObitDConCleanVisLineImageStats;
  theClass->ObitDConCleanSelect  = (ObitDConCleanSelectFP)ObitDConCleanVisLineSelect;
  theClass->ObitDConCleanVisDefWindow = (ObitDConCleanVisDefWindowFP)ObitDConCleanVisLineDefWindow;
  theClass->ObitDConCleanRestore = (ObitDConCleanRestoreFP)ObitDConCleanVisLineRestore;
  theClass->ObitDConCleanFlatten = (ObitDConCleanFlattenFP)ObitDConCleanVisLineFlatten;
  theClass->ObitDConCleanXRestore= (ObitDConCleanXRestoreFP)ObitDConCleanVisLineXRestore;
  theClass->ObitDConCleanAutoWindow = 
    (ObitDConCleanAutoWindowFP)ObitDConCleanVisLineAutoWindow;

  /* Private functions for derived classes */
  theClass->MakeResiduals   = (MakeResidualsFP)MakeResiduals;
  theClass->MakeAllResiduals= (MakeAllResidualsFP)MakeAllResiduals;
  theClass->SubNewCCs       = (SubNewCCsFP)SubNewCCs;
  theClass->NewPxList       = (NewPxListFP)NewPxList;
  theClass->NewPxArray      = (NewPxArrayFP)NewPxArray;
  theClass->KillBeamPatches = (KillBeamPatchesFP)KillBeamPatches;
  theClass->PickNext2D      = (PickNext2DFP)PickNext2D;
  theClass->PickNext3D      = (PickNext3DFP)PickNext3D;
  theClass->KillPxArray     = (KillPxArrayFP)KillPxArray;
  theClass->ResetSkyModel   = (ResetSkyModelFP)ResetSkyModel;
  theClass->ResetPixelList  = (ResetPixelListFP)ResetPixelList;
  theClass->FindPeak        = (FindPeakFP)FindPeak;
} /* end ObitDConCleanVisLineClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitDConCleanVisLineInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitDConCleanVisLine *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->nPar       = 0;
  in->nArgs      = 0;
  in->chArgs     = NULL;
  in->motherShip = NULL;
} /* end ObitDConCleanVisLineInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitDConCleanVisLine* cast to an Obit*.
 */
void ObitDConCleanVisLineClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitDConCleanVisLine *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */

  if (in->nArgs>0)  KillChanFuncArgs (in->nArgs, (ChanFuncArg**)in->chArgs);
  in->chArgs     = NULL;

  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitDConCleanVisLineClear */

/**
 * Make arguments for Threaded Single channel operation
 * Arrays of structures for independent parallel channel cleans.
 * This should only be built on the first (main) CLEAN object.
 * \param in         CLEAN structure, copies made in output
 * \param thread     ObitThread object to be used for function
 * \param err        Obit error stack object.
 * \param ThreadArgs [out] Created array of ChanFuncArg, 
 *                   delete with KillChanFuncArgs
 * \return number of elements in args (number of allowed threads).
 */
static olong MakeChanFuncArgs (ObitDConCleanVisLine *in, ObitThread *thread, 
			       ObitErr *err, ChanFuncArg ***ThreadArgs)
{
  olong i, j, nThreads;

  /* Setup for threading */
  /* How many threads? */
  nThreads = MAX (1, ObitThreadNumProc(thread));

  /* Initialize threadArg array */
  *ThreadArgs = g_malloc0(nThreads*sizeof(ChanFuncArg*));
  for (i=0; i<nThreads; i++) 
    (*ThreadArgs)[i] = g_malloc0(sizeof(ChanFuncArg)); 
  for (i=0; i<nThreads; i++) {
    (*ThreadArgs)[i]->in             = ObitDConCleanVisLineCopy(in, NULL, err);
    (*ThreadArgs)[i]->in->motherShip = (Obit*)in;  /* NEVER Unref this - BAD things will happen */
    (*ThreadArgs)[i]->chann          = i+1;
    (*ThreadArgs)[i]->done           = FALSE;
    (*ThreadArgs)[i]->bcomp          = g_malloc0(in->nfield*sizeof(olong));
    for (j=0; j<in->nfield; j++) (*ThreadArgs)[i]->bcomp[j] = 1;
    (*ThreadArgs)[i]->pixarray       = NULL;
    (*ThreadArgs)[i]->ithread        = i;
    (*ThreadArgs)[i]->thread         = thread;
    (*ThreadArgs)[i]->err            = err;
  }

  return nThreads;
} /*  end MakeChanFuncArgs */

/**
 * Delete arguments for Threaded Single channel operation
 * \param nargs      number of elements in args.
 * \param ThreadArgs Array of ChanFuncArg
 */
static void KillChanFuncArgs (olong nargs, ChanFuncArg **ThreadArgs)
{
  olong i;

  if (ThreadArgs==NULL) return;
  
  ObitThreadPoolFree (ThreadArgs[0]->thread);  /* Free thread pool */
  for (i=0; i<nargs; i++) {
    if (ThreadArgs[i]) {
      g_free((ThreadArgs)[i]->bcomp);
      (ThreadArgs)[i]->pixarray = ObitFArrayUnref((ThreadArgs)[i]->pixarray);
      (ThreadArgs)[i]->in       = ObitDConCleanVisLineUnref((ThreadArgs)[i]->in);
      g_free(ThreadArgs[i]);
    }
  }
  g_free(ThreadArgs);
} /*  end KillChanFuncArgs */

/**
 * Create Pixel list for cleaning, one per parallel argument
 * \param in       The Clean object
 * \param err      Obit error stack object.
 * \return TRUE if attempted, FALSE if cannot do 
 */
static void NewPxList (ObitDConCleanVis *inn, ObitErr *err)
{
  ObitDConCleanVisLine *in;
  ChanFuncArg *inArr=NULL;
  const ObitDConCleanVisClassInfo *ParentClass = myClassInfo.ParentClass;
  olong ichan;
  /*gchar *routine = "ObitDConCleanVis:NewPxList";*/

  /* Cast input to this type */
  in = (ObitDConCleanVisLine*)inn;

   /* error checks */
  if (err->error) return;

  /* Set on main */
  ParentClass->NewPxList(inn, err);
  /* Loop over channels */
  for (ichan=0; ichan<in->nArgs; ichan++) {
    inArr = ((ChanFuncArg*)in->chArgs[ichan]);

    /* Call parent for operation */
    ParentClass->NewPxList((ObitDConCleanVis*)inArr->in, err);
  } /* End channel loop */    

} /* end NewPxList */

/**
 * Create Pixel array for intermediate CLEAN, one per parallel argument
 * \param in       The Clean object
 * \param startCC  [out] Current number of components per field
 * \param err      Obit error stack object.
 * \return array of ObitFArrays for cleaning
 */
static ObitFArray** NewPxArray (ObitDConCleanVis *inn, olong *startCC, 
				ObitErr *err)
{
  ObitDConCleanVisLine *in;
  ChanFuncArg*inArr=NULL;
  const ObitDConCleanVisClassInfo *ParentClass = myClassInfo.ParentClass;
  ObitFArray** pixarray;
  olong ichan;
  /*gchar *routine = "ObitDConCleanVis:NewPxArray";*/

  /* Cast input to this type */
  in = (ObitDConCleanVisLine*)inn;

   /* error checks */
  if (err->error) return NULL;

  /* Main */
  in->plane[0] = 1;
  pixarray = ParentClass->NewPxArray(inn, startCC, err);

  /* Loop over channels */
  for (ichan=0; ichan<in->nPar; ichan++) {
    inArr = ((ChanFuncArg*)in->chArgs[ichan]);

    /* Set channel */
    inArr->in->plane[0] = ichan+1;
    
    /* Call parent for operation */
    inArr->pixarray = ParentClass->NewPxArray((ObitDConCleanVis*)inArr->in, startCC, err);
  } /* End channel loop */    

  return ((ChanFuncArg*)in->chArgs[0])->pixarray;
} /* end NewPxArray */

/**
 * Delete Pixel array for intermediate CLEAN, one per parallel argument
 * \param in       The Clean object
 * \param pixarray Array to delete
 * \return NULL
 */
static ObitFArray** KillPxArray (ObitDConCleanVis *inn,  ObitFArray **pixarray)
{
  ObitDConCleanVisLine *in;
  ChanFuncArg *inArr=NULL;
  const ObitDConCleanVisClassInfo *ParentClass = myClassInfo.ParentClass;
  olong ichan;

  /* Cast input to this type */
  in = (ObitDConCleanVisLine*)inn;

  /* Loop over channels */
  for (ichan=0; ichan<in->nPar; ichan++) {
    inArr = ((ChanFuncArg*)in->chArgs[ichan]);

    /* Call parent for operation */
    inArr->pixarray = ParentClass->KillPxArray((ObitDConCleanVis*)inArr->in, inArr->pixarray);
  } /* End channel loop */    

  return pixarray;
} /* end KillPxArray */

/**
 * Delete Beam Patches on in, one per parallel argument
 * \param in       The Clean object
 */
static void KillBeamPatches (ObitDConCleanVis *inn)
{
  ObitDConCleanVisLine *in = (ObitDConCleanVisLine*)inn;
  ChanFuncArg *inArr=NULL;
  const ObitDConCleanVisClassInfo *ParentClass = myClassInfo.ParentClass;
  olong ichan;

  /* Main */
   ParentClass->KillBeamPatches(inn);

  /* Loop over channels */
  for (ichan=0; ichan<in->nPar; ichan++) {
    inArr = ((ChanFuncArg*)in->chArgs[ichan]);

    /* Call parent for operation */
    ParentClass->KillBeamPatches((ObitDConCleanVis*)inArr->in);
  } /* End channel loop */    
} /* end KillBeamPatches */

/**
 * Low accuracy subtract pixels from image for current CLEAN fields.
 * Loops over channels
 * Loops over fields in in->currentFields
 * Uses BeamPatch to subtract list of components
 * Updates cleanable member on in
 * \param in       The Clean object
 * \param newCC    Start CC -1 per in->mosaic field for subtraction
 * \param pixarray Array of arrays of pixels to subtract, in order of fields in 
 *                 in->currentFields
 * \param err      Obit error stack object.
 * \return TRUE if attempted, FALSE if cannot do 
 */
static void SubNewCCs (ObitDConCleanVis *inn, olong *newCC, ObitFArray **pixarray, 
		       ObitErr *err)
{
  ObitDConCleanVisLine *in = (ObitDConCleanVisLine*)inn;
  ChanFuncArg *inArr=NULL;
  const ObitDConCleanVisClassInfo *ParentClass = myClassInfo.ParentClass;
  olong ichan, *lnewCC=NULL, ifld;
  gchar *routine = "LineSubNewCCs";

  lnewCC   = g_malloc0(in->nfield*sizeof(olong));
 
 /* Loop over channels */
  for (ichan=0; ichan<in->nPar; ichan++) {
    inArr = ((ChanFuncArg*)in->chArgs[ichan]);
    
    /* Prior number of CCs */
    for (ifld=0; ifld<inArr->in->numCurrentField; ifld++) {
      	lnewCC[ifld] = inArr->in->Pixels->iterField[inArr->in->currentFields[ifld]-1];
    }

    /* Call parent for operation */
    ParentClass->SubNewCCs((ObitDConCleanVis*)inArr->in, lnewCC, inArr->pixarray, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
  } /* End channel loop */    
  
  g_free(lnewCC);
} /* end SubNewCCs */

/**
 * 2D version - multiple parallel facets possible, loops over channels
 * Most work done in parent function.
 * Pick next fields to clean and create residual image in mosaic member
 * The selected fields are in->currentFields
 * Looks for "autoCenFlux" value in the first autoWindow image.
 * Adopted from the AIPSish CLBSTF (QOOP:QCLEAN.FOR)
 * \param in   The object to deconvolve
 * \param err Obit error stack object.
 * \return TRUE iff all channels reached minimum flux density or 
 *         max. number  comp. or no fields with SNR>5
 */
static gboolean PickNext2D(ObitDConCleanVis *inn, ObitErr *err)
{
  gboolean done=TRUE, t;
  ObitDConCleanVisLine *in = (ObitDConCleanVisLine*)inn;
  ChanFuncArg *inArr=NULL;
  const ObitDConCleanVisClassInfo *ParentClass = myClassInfo.ParentClass;
  olong ichan, i;
  gchar *routine = "LinePickNext2D";

  /* error checks */
  if (err->error) return done;

  /* Main values reset */
  in->peakFlux = 0.0;
  for (i=0; i<in->nfield; i++) in->minFlux[i] = 0.0;
  
  /* Loop over channels */
  for (ichan=0; ichan<in->nPar; ichan++) {
    inArr = ((ChanFuncArg*)in->chArgs[ichan]);
    /* Call parent for operation */
    t = ParentClass->PickNext2D((ObitDConCleanVis*)inArr->in, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, done);
    done = done && t;

    /* Set some values on main CLEAN */
    in->peakFlux = MAX (in->peakFlux, inArr->in->peakFlux);
    for (i=0; i<in->nfield; i++) in->minFlux[i] = MAX(in->minFlux[i], inArr->in->minFlux[i]);
 
    /* Fields selected */
    if (inArr->in->numCurrentField>in->numCurrentField) {
      in->numCurrentField = inArr->in->numCurrentField;
       for (i=0; i<in->numCurrentField; i++) in->currentFields[i] = inArr->in->currentFields[i];
    }
  } /* End channel loop */    

  return done;
} /* end PickNext2D */

/**
 * 3D version, loops over channels
 * Most work done in parent function.
 * Pick next field to clean.
 * The selected field is in->currentFields[0]
 * \param in   The object to deconvolve
 * \param err Obit error stack object.
 * \return TRUE iff all channels reached minimum flux density or 
 *         max. number  comp. or no fields with SNR>5
 */
static gboolean PickNext3D(ObitDConCleanVis *inn, ObitErr *err)
{
  gboolean done=TRUE, t;
  ObitDConCleanVisLine *in = (ObitDConCleanVisLine*)inn;
  ChanFuncArg *inArr=NULL;
  const ObitDConCleanVisClassInfo *ParentClass = myClassInfo.ParentClass;
  olong ichan, i;
  gchar *routine = "LinePickNext3D";

   /* error checks */
  if (err->error) return done;

  /* Main values reset */
  in->peakFlux = 0.0;
  for (i=0; i<in->nfield; i++) in->minFlux[i] = 0.0;
  
  /* Loop over channels */
  for (ichan=0; ichan<in->nPar; ichan++) {
    inArr = ((ChanFuncArg*)in->chArgs[ichan]);
    /* Call parent for operation */
    t = ParentClass->PickNext3D((ObitDConCleanVis*)inArr->in,err);
    if (err->error) Obit_traceback_val (err, routine, in->name, done);

    done = done && t;
 
    /* Set some values on main CLEAN */
    in->peakFlux = MAX (in->peakFlux, inArr->in->peakFlux);
    for (i=0; i<in->nfield; i++) in->minFlux[i] = MAX(in->minFlux[i], inArr->in->minFlux[i]);
 
    /* Fields selected */
    if (inArr->in->numCurrentField>in->numCurrentField) {
      in->numCurrentField = inArr->in->numCurrentField;
       for (i=0; i<in->numCurrentField; i++) in->currentFields[i] = inArr->in->currentFields[i];
    }
  } /* End channel loop */    

  return done;
} /* end PickNext3D */

/**
 * Reset Sky Model, one per parallel argument
 * If reuseFlux member is > 0 then any components in the sky model above 
 * this level are left and TRUE is returned if any exist.
 * Resets the number of components
 * the reuseFlux option uses the Pixels member and the list of CC Tables
 * which must be fully instantiated.
 * Sets BChan, EChan, BIF, EIF to default
 * \param in   The Clean object
 * \param err Obit error stack object.
 * \return true if components in the sky model need to be subtracted
 */
static gboolean ResetSkyModel(ObitDConCleanVis *inn, ObitErr *err)
{
  gboolean doSub=FALSE, tdoSub;
  ObitDConCleanVisLine *in = (ObitDConCleanVisLine*)inn;
  ChanFuncArg *inArr=NULL;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  const ObitDConCleanVisClassInfo *ParentClass = myClassInfo.ParentClass;
  olong ichan, CCVer;
  gchar *pixelParms[] = { /* parameters on Pixels for CLEAN */
    "Niter", "Gain", "minFlux", "Factor", "fGauss", "ccfLim", "prtLv",
    NULL};
  gchar *routine = "LineResetSkyModel";

  /* error checks */
  if (err->error) return doSub;

  /* Do main */
  doSub = ParentClass->ResetSkyModel(inn, err);
  /* Loop over channels */
  for (ichan=0; ichan<in->nArgs; ichan++) {
    inArr = ((ChanFuncArg*)in->chArgs[ichan]);

    /* Get input CLEAN parameters  */
    ObitInfoListCopyList (in->Pixels->info, inArr->in->Pixels->info, pixelParms);
    if (err->error) Obit_traceback_val (err, routine, in->name, doSub);
    
   /* Set CC Table to use in Pixel Array */
    CCVer = ichan+1;
    ObitInfoListAlwaysPut (inArr->in->Pixels->info, "CCVer", OBIT_long, dim, &CCVer);

    /* Call parent for operation */
   tdoSub =  ParentClass->ResetSkyModel((ObitDConCleanVis*)inArr->in, err);
   if (err->error) Obit_traceback_val (err, routine, in->name, doSub);
   doSub  = doSub || tdoSub;
  } /* End channel loop */    

  return doSub;
} /* end ResetSkyModel */
/**
 * Reset Pixel List for beginning of a CLEAN
 * \param in   The Clean object
 * \param err Obit error stack object.
 */
static void ResetPixelList(ObitDConCleanVis *inn, ObitErr *err)
{
  ObitDConCleanVisLine *in=(ObitDConCleanVisLine*)inn;
  ChanFuncArg *inArr=NULL;
  const ObitDConCleanPxListClassInfo *pxListClass;
  olong ichan, CCVer;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *routine = "LineResetPixelList";

  /* Cast input to this type */
  in = (ObitDConCleanVisLine*)inn;

  /* error checks */
  if (err->error) return;
  
  /* Loop over channels */
  for (ichan=0; ichan<in->nArgs; ichan++) {
    inArr = ((ChanFuncArg*)in->chArgs[ichan]);

    /* Set CC Table to use */
    CCVer = ichan+1;
    ObitInfoListAlwaysPut (inArr->in->Pixels->info, "CCVer", OBIT_long, dim, &CCVer);

    /* PxList class structure */
    pxListClass = 
      (ObitDConCleanPxListClassInfo*)inArr->in->Pixels->ClassInfo; 
    /* Now reset */
    pxListClass->ObitDConCleanPxListReset (inArr->in->Pixels, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
  } /* End channel loop */    

  /* PxList class structure */
  pxListClass = (ObitDConCleanPxListClassInfo*)in->Pixels->ClassInfo; 

  CCVer = 1;   /* DEBUG */
  ObitInfoListAlwaysPut (in->Pixels->info, "CCVer", OBIT_long, dim, &CCVer);
  pxListClass->ObitDConCleanPxListReset (in->Pixels, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
} /* end ResetPixelList */

/**
 * Make selected residual images and get statistics, multiple channels
 * \param in     The Clean object
 * \param fields zero terminated list of field numbers to image
 * \param doBeam If TRUE also make beam
 * \param err    Obit error stack object.
 */
static void  MakeResiduals (ObitDConCleanVis *inn, olong *fields, 
			    gboolean doBeam, ObitErr *err)
{
  ObitDConCleanVisLine *tin=(ObitDConCleanVisLine*)inn, *in;
  gboolean doWeight, doFlatten, found;
  const ObitDConCleanVisLineClassInfo *inClass;
  ObitUVImagerClassInfo *imgClass = (ObitUVImagerClassInfo*)tin->imager->ClassInfo;
  ChanFuncArg *inArr=NULL;
  olong i, ifld, jfld, ichan, field=1, *stale=NULL, best;
  ofloat bestQuality;
  gchar *routine = "LineMakeResidual";

  /* May need to use upper level Clean for this, use motherShip if inn is a lower level */
  if (tin->chArgs==NULL) in = (ObitDConCleanVisLine*)tin->motherShip;
  else                   in = tin;

  inClass = (ObitDConCleanVisLineClassInfo*)in->ClassInfo; /* class structure */

  /* Are residuals fresh? */
  doWeight  = FALSE;
  doFlatten = FALSE;

  /* Make copy of only stale images */
  stale = g_malloc0((in->nfield+1)*sizeof(olong));
  ifld = jfld = 0;
  while(fields[ifld]>0) {
    if (!in->fresh[fields[ifld]-1]) stale[jfld++] = fields[ifld];
    ifld++;
  }
  stale[jfld] = 0;

  /* if none stale then bail out */
  if (stale[0]<=0) goto done;

  /* Make residual images for stale fields */
  imgClass->ObitUVImagerImage (in->imager, stale, doWeight, doBeam, doFlatten, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Find channel with highest quality field to copy values to main */
  best          = 0;
  bestQuality = 0.0;

 /* Statistics per channel */
  for (ichan=0; ichan<in->nPar; ichan++) {
    inArr = ((ChanFuncArg*)in->chArgs[ichan]);
    for (ifld=0; ifld<in->nfield; ifld++) {
      field = stale[ifld];
      if (field<=0) break;
      /* Statistics all done in first call */
      if (ichan==0) {
	/* Get statistics  for image */
	inClass->ObitDConCleanImageStats ((ObitDConClean*)in, field, FALSE, err);
	/* Need Beam statistics? */
	if (doBeam)
	  inClass->ObitDConCleanImageStats ((ObitDConClean*)in, field, TRUE, err);
	if (err->error) Obit_traceback_msg (err, routine, in->name);
      } /* End get image statistics */
      
      /* Quality measure */
      inArr->in->quality[field-1] = 
	ObitDConCleanVisQuality((ObitDConCleanVis*)inArr->in, field, err);
      
      /* Max cleanable flux */
      inArr->in->cleanable[field-1] = 
	inClass->ObitDConCleanVisCleanable((ObitDConCleanVis*)inArr->in, field, NULL, err);
      inArr-> in->fresh[field-1]    = TRUE;
    } /* end statistics loop */
  
    if (err->prtLv>1) ObitErrLog(err);  /* Progress Report */
    else ObitErrClear(err);

    /* For any shifted fields get statistics */
    for (ifld=0; ifld<in->nfield; ifld++) {
      if (in->mosaic->isShift[ifld] > 0) {
	
	/* Has this field reached the min flux? maxAbsRes<minFlux */
	if (inArr->in->maxAbsRes[ifld]<=inArr->in->minFlux[ifld]) continue;
	jfld = in->mosaic->isShift[ifld]-1;
	/* If the corresponding autoCenter field just remade (as well as the shifted
	   version - determine the statistics from the image */
	found = FALSE;
	for (i=0; i<in->nfield; i++) {
	  if (fields[i]==0) break;
	  if ((jfld+1)==stale[i]) {found = TRUE; break;}
	}
	
	if (found) {
	  if (ichan==0) {
	    /* Get statistics  for image */
	    inClass->ObitDConCleanImageStats ((ObitDConClean*)in, ifld+1, FALSE, err);
	  }
	  /* Quality measure */
	  in->quality[ifld] = ObitDConCleanVisQuality((ObitDConCleanVis*)in, ifld+1, err);
	  /* Max cleanable flux */
	  inArr->in->cleanable[ifld]   = 
	    inClass->ObitDConCleanVisCleanable((ObitDConCleanVis*)inArr->in, ifld+1, NULL, err);
	  inArr->in->beamPeakRMS[ifld] = inArr->in->beamPeakRMS[jfld];
	  inArr->in->fresh[ifld]       = TRUE;
	} else { /* Not remade - Use residual */
	  inArr->in->maxAbsRes[ifld]   = inArr->in->maxAbsRes[jfld];
	  inArr->in->quality[jfld]     = inArr->in->maxAbsRes[jfld];  /* Reset Quality to residual */
	  inArr->in->quality[ifld]     = inArr->in->quality[jfld];
	  inArr->in->cleanable[ifld]   = MIN(inArr->in->cleanable[jfld], inArr->in->maxAbsRes[jfld]);
	  inArr->in->beamPeakRMS[ifld] = inArr->in->beamPeakRMS[jfld];
	  inArr->in->fresh[ifld]       = inArr->in->fresh[jfld];
	}
      }
    } /* end loop over fields */
    /* find best quality */
    for (ifld=0; ifld<in->nfield; ifld++) {
      if (inArr->in->currentFields[ifld] <= 0) break;
      if (inArr->in->quality[ifld]>bestQuality) {
	best        = ichan;
	bestQuality = inArr->in->quality[ifld];
      }
    }
  } /* end channel loop */

  /* Copy stats info for best channel to main */
  inArr = ((ChanFuncArg*)in->chArgs[best]);

  for (ifld=0; ifld<in->nfield; ifld++) {
    if (doBeam) in->beamPeakRMS[ifld] = inArr->in->beamPeakRMS[ifld];
    in->maxAbsRes[ifld]  = inArr->in->maxAbsRes[ifld];
    in->avgRes[ifld]     = inArr->in->avgRes[ifld];
    in->imgRMS[ifld]     = inArr->in->imgRMS[ifld];
    in->imgPeakRMS[ifld] = inArr->in->imgPeakRMS[ifld];
    in->beamPeakRMS[ifld]= inArr->in->beamPeakRMS[ifld];
    in->quality[ifld]    = inArr->in->quality[ifld];
    in->cleanable[ifld]  = inArr->in->cleanable[ifld];
    in->fresh[ifld]      = inArr->in->fresh[ifld];
  } /* end loop over fields */
  
 done:   /* finished */
  if(stale) g_free(stale);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
} /* end MakeResiduals */

/**
 * Make all residual images and get statistics
 * Multiple channel version
 * \param in     The Clean object
 * \param err    Obit error stack object.
 */
static void  MakeAllResiduals (ObitDConCleanVis *inn, ObitErr *err)
{
  ObitDConCleanVisLine *tin=(ObitDConCleanVisLine*)inn, *in;
  gboolean doBeam = inn->doBeam;
  gboolean doWeight, doFlatten;
  olong i, ifld, jfld, ichan, best, fields[2]={0,0};
  ChanFuncArg *inArr=NULL;
  ObitUVImagerClassInfo *imgClass = 
    (ObitUVImagerClassInfo*)tin->imager->ClassInfo;
  const ObitDConCleanVisLineClassInfo *inClass;
  ofloat bestQuality;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *routine = "LineMakeAllResiduals";

  if (err->error) return;   /* existing error */
  
  /* Need to use upper level Clean for this, use motherShip if inn is a lower level */
  if (tin->chArgs==NULL) in = (ObitDConCleanVisLine*)tin->motherShip;
  else                   in = tin;

  inClass = (ObitDConCleanVisLineClassInfo*)in->ClassInfo; /* class structure */

  /* Turn off things not needed */
  doWeight  = FALSE;
  doFlatten = FALSE;

  /* Copy prtLv to in->mosaic->info */
  dim[0] = 1;dim[1] = 1;
  ObitInfoListAlwaysPut (in->mosaic->info, "prtLv", OBIT_long, dim, &err->prtLv);

  /* Parallel Image images without needing beam */
  imgClass->ObitUVImagerImage (in->imager, fields,  doWeight, doBeam, doFlatten, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  
   /* Find channel with highest quality field to copy values to main */
  best          = 0;
  bestQuality = 0.0;

  /* Statistics per channel */
  for (ichan=0; ichan<in->nPar; ichan++) {
    inArr = ((ChanFuncArg*)in->chArgs[ichan]);
    if (inArr->in->Pixels->currentIter<0) {
      inArr->done = FALSE;  /* initialize inArr */
    }
    /* Loop over fields getting statistics for Image and Beam */
    for (i=0; i<in->nfield; i++) {
      /* Statistics all done in first call */
      if (ichan==0) {
	inClass->ObitDConCleanImageStats ((ObitDConClean*)in, i+1, FALSE, err);
	if (doBeam)
	  inClass->ObitDConCleanImageStats ((ObitDConClean*)in, i+1, TRUE, err);
      } /* end get statistics */
      if (err->error) Obit_traceback_msg (err, routine, in->name);
      
      /* Quality measure */
      inArr->in->quality[i] = 
	ObitDConCleanVisQuality((ObitDConCleanVis*)inArr->in, i+1, err);
      /* Max cleanable flux */
      inArr->in->cleanable[i] = 
	inClass->ObitDConCleanVisCleanable((ObitDConCleanVis*)inArr->in, i+1, NULL, err);
      inArr->in->fresh[i]     = TRUE;  /* Freshly made image */
      if (err->error) Obit_traceback_msg (err, routine, in->name);
    } /* end loop over field */
    
    /* For any shifted fields get statistics */
    for (ifld=0; ifld<in->mosaic->numberImages; ifld++) {
      if (in->mosaic->isShift[ifld] > 0) {
	jfld = in->mosaic->isShift[ifld]-1;
	/* Get statistics  for image */
	if (ichan==0) 
	  inClass->ObitDConCleanImageStats ((ObitDConClean*)in, ifld+1, FALSE, err);
	/* Quality measure */
	inArr->in->quality[ifld] = 
	  ObitDConCleanVisQuality((ObitDConCleanVis*)inArr->in, ifld+1, err);
	/* Max cleanable flux */
	inArr->in->cleanable[ifld]   = 
	  inClass->ObitDConCleanVisCleanable((ObitDConCleanVis*)inArr->in, ifld+1, NULL, err);
	inArr->in->beamPeakRMS[ifld] = inArr->in->beamPeakRMS[jfld];
	inArr->in->fresh[ifld]       = inArr->in->fresh[jfld];
      }
    } /* end loop over fields */
    /* find best quality */
    for (ifld=0; ifld<in->nfield; ifld++) {
      if (inArr->in->currentFields[ifld] <= 0) break;
      if (inArr->in->quality[ifld]>bestQuality) {
	best        = ichan;
	bestQuality = inArr->in->quality[ifld];
      }
    }
  } /* end channel loop */

  /* Copy stats info for best channel to main */
  inArr = ((ChanFuncArg*)in->chArgs[best]);

  for (ifld=0; ifld<in->nfield; ifld++) {
    if (doBeam) in->beamPeakRMS[ifld] = inArr->in->beamPeakRMS[ifld];
    in->maxAbsRes[ifld]  = inArr->in->maxAbsRes[ifld];
    in->avgRes[ifld]     = inArr->in->avgRes[ifld];
    in->imgRMS[ifld]     = inArr->in->imgRMS[ifld];
    in->imgPeakRMS[ifld] = inArr->in->imgPeakRMS[ifld];
    in->beamPeakRMS[ifld]= inArr->in->beamPeakRMS[ifld];
    in->quality[ifld]    = inArr->in->quality[ifld];
    in->cleanable[ifld]  = inArr->in->cleanable[ifld];
    in->fresh[ifld]      = inArr->in->fresh[ifld];
  } /* end loop over fields */
} /* end MakeAllResiduals */

/**
 * Find maximum brightness using CC tables
 * Set in->peakFlux to the larger of in->peakFlux and the maximum 
 * brightness in the CC tables.
 * \param in      Object of interest.
 * \param err     ObitErr for reporting errors.
 */
static void FindPeak (ObitDConCleanVis *inn, ObitErr *err)
{
  ObitDConCleanVisLine *in = (ObitDConCleanVisLine*)inn;
  const ObitDConCleanVisClassInfo *ParentClass = myClassInfo.ParentClass;
  ChanFuncArg *inArr=NULL;
  olong ichan;
  gchar *routine = "LineFindPeak";

  /* Main */
  in->plane[0] = 1;
  
  /* Loop over channels */
  for (ichan=0; ichan<in->nPar; ichan++) {
    inArr = ((ChanFuncArg*)in->chArgs[ichan]);

    /* Set plane */
    inArr->in->plane[0] = ichan+1;
    inArr->in->nfield   = in->nfield;
    inArr->in->CCver    = ichan+1;

    /* Call parent for operation - ignore if no CLEAN - this uses CCs */
    if (inArr->in->Pixels->currentIter>0) 
      ParentClass->FindPeak((ObitDConCleanVis*)inArr->in, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
  } /* End channel loop */    
} /* end FindPeak */
