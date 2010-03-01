/* $Id: ObitUVImagerMF.c 129 2009-09-27 22:00:59Z bill.cotton $        */
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
/*;Correspondence about this software should be addressed as follows: */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/

#include "ObitUVImagerMF.h"
#include "ObitUVImager.h"
#include "ObitImageMF.h"
#include "ObitUVWeight.h"
#include "ObitImageUtil.h"
#include "ObitImageMosaicMF.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitUVImagerMF.c
 * ObitUVImagerMF class function definitions.
 * This class is derived from the Obit base class.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitUVImagerMF";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitUVImagerGetClass;

/**
 * ClassInfo structure ObitUVImagerMFClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitUVImagerMFClassInfo myClassInfo = {FALSE};

/*--------------- File Global Variables  ----------------*/


/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitUVImagerMFInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitUVImagerMFClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitUVImagerMFClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitUVImagerMF* newObitUVImagerMF (gchar* name)
{
  ObitUVImagerMF* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitUVImagerMFClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitUVImagerMF));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitUVImagerMFInit((gpointer)out);

 return out;
} /* end newObitUVImagerMF */

/**
 * Initializes from ObitInfoList.
 * Initializes class if needed on first call.
 * \param out     the new object.to be initialized
 * \param prefix  If NonNull, string to be added to beginning of inList entry name
 *                "xxx" in the following
 * \param inList  InfoList to extract object information from 
 *      \li "xxxClassType" string UVImager type, "MF" for this class
 *      \li "xxxUVData" prefix for uvdata member, entry with value "None" => doesn't exist
 *      \li "xxxUVWork" prefix for uvwork member, entry with value "None" => doesn't exist
 *      \li "xxxMosaic" prefix for mosaic member, entry with value "None" => doesn't exist
 * \param err     ObitErr for reporting errors.
 */
void ObitUVImagerMFFromInfo (ObitUVImager *out, gchar *prefix, ObitInfoList *inList, 
			     ObitErr *err)
{ 
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *keyword=NULL, *value=NULL;
  gboolean missing;
  gchar *Type = "MF";
  gchar *routine = "ObitUVImagerMFFromInfo";
  
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitUVImagerMFClassInit();

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(out, &myClassInfo));

  /* check class type */
  missing = ObitInfoListGetP(inList, keyword, &type, dim, (gpointer*)&value);
  if ((missing) || (type!=OBIT_string) || (!strncmp(Type,value,dim[0]))) {
    Obit_log_error(err, OBIT_Error,"%s Wrong class type %s!=%s", routine, value, Type);
    return;
  }

} /* end ObitUVImagerMFFromInfo */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitUVImagerMFGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitUVImagerMFClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitUVImagerMFGetClass */

/**
 * Make a deep copy of an ObitUVImagerMF.
 * \param int   The object to copy
 * \param outt  An existing object pointer for output or NULL if none exists.
 * \param err  Obit error stack object.
 * \return pointer to the new object.
 */
ObitUVImager* ObitUVImagerMFCopy  (ObitUVImager *inn, ObitUVImager *outt, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  ObitUVImagerMF *in  = (ObitUVImagerMF*)inn;
  ObitUVImagerMF *out = (ObitUVImagerMF*)outt;
  gboolean oldExist;
  gchar *outName;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return outt;
  g_assert (ObitIsA(in, &myClassInfo));
  if (out) g_assert (ObitIsA(out, &myClassInfo));

  /* Create if it doesn't exist */
  oldExist = out!=NULL;
  if (!oldExist) {
    /* derive object name */
    outName = g_strconcat ("Copy: ",in->name,NULL);
    out = newObitUVImagerMF(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (inn, outt, err);

  /*  copy this class - just pointers */

  return (ObitUVImager*)out;
} /* end ObitUVImagerMFCopy */

/**
 * Make a copy of a object but do not copy the actual data
 * This is useful to create an UVImagerMF similar to the input one.
 * \param inn  The object to copy
 * \param outt  An existing object pointer for output, must be defined.
 * \param err  Obit error stack object.
 */
void ObitUVImagerMFClone  (ObitUVImager *inn, ObitUVImager *outt, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  ObitUVImagerMF *in  = ( ObitUVImagerMF*)inn;
  ObitUVImagerMF *out = ( ObitUVImagerMF*)outt;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitIsA(out, &myClassInfo));

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class - just pointers */

} /* end ObitUVImagerMFClone */

/**
 * Creates an ObitUVImagerMF given an ObitUV with control information.
 * The output ImageMosaic member is created
 * \param name   An optional name for the object.
 * \param order  Spectral imaging order,0=flux,1=si, 2=curve
 * \param maxFBW Maximum IF center fractional bandwidth.
 * \param alpha  Spectral index correction applied to uv data making mosaic
 * \param uvdata ObitUV object with info member containng the output image
 *               specifications and all processing parameters.
 * \li FileType = Underlying file type, OBIT_IO_FITS, OBIT_IO_AIPS
 * \li Name     = Name of image, used as AIPS name or to derive FITS filename
 * \li Class    = Root of class, used as AIPS class or to derive FITS filename
 * \li Seq      = Sequence number
 * \li Disk     = Disk number for underlying files
 * \li FOV      = Field of view (deg) for Mosaic 
 *                If > 0.0 then a mosaic of images will be added to cover this region.
 *                Note: these are in addition to the NField fields added by 
 *                other parameters
 * \li doFull   = if TRUE, create full field image to cover FOV [def. FALSE]
 * \li NField   = Number of fields defined in input,
 *                if unspecified derive from data and FOV
 * \li "xCells" = Cell spacing in X (asec) for all images,
 *                if unspecified derive from data
 * \li "yCells" = Cell spacing in Y (asec) for all images,
 *                if unspecified derive from data
 * \li "BMAJ"   = OBIT_float scalar = Restoring beam major axis (asec)
 *                if = 0 then write fitted value to header
 * \li "BMIN"   = OBIT_float scalar = Restoring beam minor axis (asec)
 * \li "BPA"    = OBIT_float scalar = Restoring beam position angle (deg)
 * \li "Beam"   = OBIT_float [3] = (BMAJ, BMIN, BPA) alternate form
 * \li nx       = Minimum number of cells in X for NField images
 *                if unspecified derive from data
 * \li ny       = Minimum number of cells in Y for NField images
 *                if unspecified derive from data
 * \li RAShift  = Right ascension shift (AIPS convention) for each field
 *                if unspecified derive from FOV and data
 * \li DecShift = Declination for each field
 *                if unspecified derive from FOV and data
 * \li Catalog  =    AIPSVZ format catalog for defining outliers, 
 *                   'None'=don't use [default]
 *                   'Default' = use default catalog.
 *                   Assumed in FITSdata disk 1.
 * \li OutlierDist = Maximum distance (deg) from center to include outlier fields
 *                   from Catalog. [default 1 deg]
 * \li OutlierFlux = Minimum estimated flux density include outlier fields
 *                   from Catalog. [default 0.1 Jy ]
 * \li OutlierSI   = Spectral index to use to convert catalog flux density to observed
 *                   frequency.  [default = -0.75]
 * \li OutlierSize = Width of outlier field in pixels.  [default 50]
 * \param err Obit error stack object.
 * \return the new object.
 */
ObitUVImagerMF* ObitUVImagerMFCreate (gchar* name, olong order,ofloat maxFBW, 
				      ofloat alpha, ObitUV *uvdata,  ObitErr *err)
{
  ObitUVImagerMF* out=NULL;
  gchar *routine = "ObitUVImagerMFCreate";

  /* Error checks */
  if (err->error) return out;
  g_assert(ObitUVIsA(uvdata));

  /* Create basic structure */
  out = newObitUVImagerMF (name);

  /* Save uvdata */
  out->uvdata = ObitUVRef(uvdata);
  out->maxOrder = order;
  out->maxFBW   = maxFBW;

  /* Create output mosaic */
  out->mosaic = 
    (ObitImageMosaic*)ObitImageMosaicMFCreate (name, order, maxFBW, alpha,
					       uvdata, err);
  if (err->error) Obit_traceback_val (err, routine, name, out);

  /* Define images */
  ObitImageMosaicMFDefine (out->mosaic, uvdata, TRUE, err);
  if (err->error) Obit_traceback_val (err, routine, name, out);

  return out;
} /* end ObitUVImagerMFCreate */

/**
 * Creates an ObitUVImagerMF given an ObitUV with control information
 * and a previously existing ImageMosaic
 * The output ImageMosaic member is created
 * \param name   An optional name for the object.
 * \param order  Spectral imaging order,0=flux,1=si, 2=curve
 * \param uvdata ObitUV object with info member containng the output image
 *               specifications and all processing parameters.
 * \param mosaic ImageMosaicMF (as ImageMosaic) to use
 * \param err Obit error stack object.
 * \return the new object.
 */
ObitUVImagerMF* ObitUVImagerMFCreate2 (gchar* name, olong order, ObitUV *uvdata, 
				       ObitImageMosaicMF *mosaic, ObitErr *err)
{
  ObitUVImagerMF* out=NULL;
  gchar *routine = "ObitUVImagerMFCreate2 ";

  /* Error checks */
  if (err->error) return out;
  g_assert(ObitUVIsA(uvdata));

  /* Check input types */
  Obit_retval_if_fail((ObitImageMosaicMFIsA(mosaic) && 
		       (ObitImageMFIsA((ObitImageMF*)mosaic->images[0]))), err, 
		       out,
		       "%s: Image mosaic or images not MF", routine);
		      
  /* Create basic structure */
  out = newObitUVImagerMF (name);

  /* Save uvdata */
  out->uvdata   = ObitUVRef(uvdata);
  out->maxOrder = order;

  /* Save mosaic */
  out->mosaic = ObitImageMosaicRef(mosaic);

  return out;
} /* end ObitUVImagerMFCreate2 */

/**
 * Create shifted 2D imaging facets
 * \param inn    The input UVImager object
 * \param field  Zero terminated list of field numbers to image, 0=> all
 * \param doall  If true then all images were remade
 * \param err    Obit error stack object.
 */
void ObitUVImagerMFShifty (ObitUVImager *inn, olong *field, gboolean doall, 
			   ObitErr *err)
{
  ObitUVImagerMF *in = (ObitUVImagerMF*)inn;
  olong i, j, ofield, ifield;
  ofloat shift[2];
  gboolean found;
  gchar *routine = "ObitUVImagerMFShifty";

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  if (!ObitImageMosaicIsA(in->mosaic)) {
    Obit_log_error(err, OBIT_Error,"%s ImageMosaic not defined in %s", 
		   routine, in->name);
    return;
  }

  /* Loop making any (2D) shifted images */
  for (i=0; i<in->mosaic->numberImages; i++) {
    /* Only needed for 2D */
    if (in->mosaic->images[i]->myDesc->do3D) continue;
    ofield = i;
    ifield = in->mosaic->isShift[i]-1;
    if (ifield >= 0) {
      /* See if input image just made - named in field */
      found = FALSE;
      for (j=0; j<in->mosaic->numberImages; j++) {
	if (field[j]==0) break;
	if (field[j]==(ifield+1)) found = TRUE;
      }
      if (found || doall) {
	shift[0] = in->mosaic->images[ofield]->myDesc->crpix[0] - in->mosaic->images[ifield]->myDesc->crpix[0];
	shift[1] = in->mosaic->images[ofield]->myDesc->crpix[1] - in->mosaic->images[ifield]->myDesc->crpix[1];
	ObitImageUtilMFShift (in->mosaic->images[ifield], in->mosaic->images[ofield], shift, err);
	/* Beam - no shift */
	shift[0] = shift[1] = 0.0;
	ObitImageUtilMFShift ((ObitImage*)in->mosaic->images[ifield]->myBeam, 
			      (ObitImage*)in->mosaic->images[ofield]->myBeam, shift, err);
      }
    } /* end if shifted field */
  } /* End loop making shifted images */
  if (err->error) Obit_traceback_msg (err, routine, in->name);

} /* end ObitUVImageMFShifty */

/**
 * return ImageMosaic member
 * \param inn  The input object
 * \param err  Obit error stack object.
 * \return reference to ImageMosaic.
 */
ObitImageMosaic* ObitUVImagerMFGetMosaic (ObitUVImager *inn, ObitErr *err)
{ 
  ObitUVImagerMF *in  = ( ObitUVImagerMF*)inn;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return NULL;
  g_assert (ObitIsA(in, &myClassInfo));

  return ObitImageMosaicRef(in->mosaic);
} /* end ObitUVImagerMFGetMosaic */

/**
 * Convert structure information to entries in an ObitInfoList
 * \param inn     Object of interest.
 * \param prefix  If NonNull, string to be added to beginning of outList entry name
 *                "xxx" in the following
 * \param outList InfoList to write entries into
 *      \li "xxxClassType" string UVImagerMF type, "Base" for base class
 *      \li "xxxUVData" prefix for uvdata member, entry with value "None" => doesn't exist
 *      \li "xxxUVWork" prefix for uvwork member, entry with value "None" => doesn't exist
 *      \li "xxxMosaic" prefix for mosaic member, entry with value "None" => doesn't exist
 *      \li various weighting parameters
 * \param err     ObitErr for reporting errors.
 */
void ObitUVImagerMFGetInfo (ObitUVImager *inn, gchar *prefix, ObitInfoList *outList, 
			    ObitErr *err)
{ 
  ObitUVImagerMF *in  = ( ObitUVImagerMF*)inn;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *keyword=NULL, *None = "None", *OK="OK", *Type="MF";
  olong i;
  gpointer listPnt;
  gchar *parm[] = 
    {"do3D", "FOV", "doFull", "NField", "xCells", "yCells", "nx", "ny", 
     "RAShift", "DecShift", "Sources", 
     "Catalog",  "OutlierDist", "OutlierFlux", "OutlierSI", "OutlierSize",
     "nuGrid", "nvGrid", "WtBox", "WtFunc", "UVTaper", "Robust", "WtPower",
     NULL};
  gchar *routine = "ObitUVImagerMFGetInfo";

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Class Type */
  if (prefix) keyword = g_strconcat (prefix, "ClassType", NULL);
  else        keyword = g_strdup("ClassType");
  dim[0] = strlen(Type);
  ObitInfoListAlwaysPut(outList, keyword, OBIT_string, dim, Type);
  g_free(keyword);

  /* uv data */
  if (prefix) keyword = g_strconcat (prefix, "UVData", NULL);
  else        keyword = g_strdup("UVData");
  if (in->uvdata) {
    ObitDataGetFileInfo((ObitData*)in->uvdata, keyword, outList, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    dim[0] = strlen(OK);
    ObitInfoListAlwaysPut(outList, keyword, OBIT_string, dim, OK);
  } else {
    dim[0] = strlen(None);
    ObitInfoListAlwaysPut(outList, keyword, OBIT_string, dim, None);
  }
  g_free(keyword);

  /* uv work */
  if (prefix) keyword = g_strconcat (prefix, "UVWork", NULL);
  else        keyword = g_strdup("UVWork");
  if (in->uvdata) {
    ObitDataGetFileInfo((ObitData*)in->uvwork, keyword, outList, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    dim[0] = strlen(OK);
    ObitInfoListAlwaysPut(outList, keyword, OBIT_string, dim, OK);
  } else {
    dim[0] = strlen(None);
    ObitInfoListAlwaysPut(outList, keyword, OBIT_string, dim, None);
  }
  g_free(keyword);

  /* ImageMosaic */
  if (prefix) keyword = g_strconcat (prefix, "Mosaic", NULL);
  else        keyword = g_strdup("Mosaic");
  if (in->mosaic) {
    ObitImageMosaicGetInfo(in->mosaic, keyword, outList, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    dim[0] = strlen(OK);
    ObitInfoListAlwaysPut(outList, keyword, OBIT_string, dim, OK);
  } else {
    dim[0] = strlen(None);
    ObitInfoListAlwaysPut(outList, keyword, OBIT_string, dim, None);
  }
  g_free(keyword);

  /* Weighting/Imaging parameters */
  i = 0;
  while (parm[i]) {
    if (prefix) keyword = g_strconcat (prefix, parm[i], NULL);
    else        keyword = g_strdup(parm[i]);
    if (ObitInfoListGetP(in->uvdata->info, parm[i], &type, dim, (gpointer*)&listPnt)) {
      ObitInfoListAlwaysPut(outList, keyword, type, dim, listPnt);
    }
    i++;
    g_free(keyword);
  }

} /* end ObitUVImagerMFGetInfo */

/**
 * Get number of parallel images
 * Target memory usage is 1 GByte.
 * \param inn     Object of interest.
 * \return the number of parallel images.
 */
olong ObitUVImagerMFGetNumPar (ObitUVImager *inn, ObitErr *err)
{
  ObitUVImagerMF *in  = (ObitUVImagerMF*)inn;
  olong out=8, nSpec;
  odouble lenVis, numVis, imSize, bufSize;

  if (err->error) return out;
  g_assert(ObitUVImagerMFIsA(in));

  /* How big are things? */
  numVis = (odouble)ObitImageUtilBufSize (in->uvdata);  /* Size of buffer */
  lenVis = (odouble)in->uvdata->myDesc->lrec;
  imSize = 2.0 * in->mosaic->images[0]->myDesc->inaxes[0] * 
    in->mosaic->images[0]->myDesc->inaxes[1];  /* Image plane size */

  nSpec = ((ObitImageMF*)in->mosaic->images[0])->nSpec;
  bufSize = numVis*lenVis + imSize*nSpec;  /* Approx memory (words) per parallel image */
  bufSize *= sizeof(ofloat);               /* to bytes */

  out = 1.0e9 / bufSize;  /* How many fit in a gByte? */

  return out;
} /*  end ObitUVImagerMFGetNumPar */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitUVImagerMFClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitUVImagerMFClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitUVImagerMFClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitUVImagerMFClassInfoDefFn (gpointer inClass)
{
  ObitUVImagerMFClassInfo *theClass = (ObitUVImagerMFClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit      = (ObitClassInitFP)ObitUVImagerMFClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitUVImagerMFClassInfoDefFn;
  theClass->ObitGetClass       = (ObitGetClassFP)ObitUVImagerMFGetClass;
  theClass->ObitUVImagerCreate = (ObitUVImagerCreateFP)ObitUVImagerCreate;
  theClass->ObitUVImagerCreate2= (ObitUVImagerCreate2FP)ObitUVImagerMFCreate2;
  theClass->ObitUVImagerGetInfo= (ObitUVImagerGetInfoFP)ObitUVImagerMFGetInfo;
  theClass->ObitUVImagerGetNumPar= (ObitUVImagerGetNumParFP)ObitUVImagerMFGetNumPar;
  theClass->ObitUVImagerShifty = (ObitUVImagerShiftyFP)ObitUVImagerMFShifty;

} /* end ObitUVImagerMFClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitUVImagerMFInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitUVImagerMF *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */

} /* end ObitUVImagerMFInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitUVImagerMF* cast to an Obit*.
 */
void ObitUVImagerMFClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitUVImagerMF *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitUVImagerMFClear */

