/* $Id$        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2010-2021                                          */
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
#include "ObitImageMosaic.h"
#include "ObitImageMosaicMF.h"
#include "ObitUVImagerMFDef.h"
#include "ObitUVImagerMF.h"
#include "ObitUVImager.h"
#include "ObitImageMF.h"
#include "ObitUVWeight.h"
#include "ObitImageUtil.h"
#ifdef __APPLE__
#include <sys/sysctl.h>
#else
#include "sys/sysinfo.h"
#endif
#include "unistd.h"

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
 * \param alphaRefF Reference frequency for alpha
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
 * \li CatDisk  =    FITS disk for Catalog [def 1]
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
				      ofloat alpha, odouble alphaRefF, 
				      ObitUV *uvdata,  ObitErr *err)
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

   /* Open and close uvdata to set descriptor for scratch file */
  ObitUVOpen (out->uvdata, OBIT_IO_ReadCal, err);
  ObitUVClose (out->uvdata, err);
  if (err->error) Obit_traceback_val (err, routine, name, out);

  /* Create scratch uvwork if it doesn't exist */
  if (out->uvwork==NULL) out->uvwork = newObitUVScratch (out->uvdata, err);
  if (err->error) Obit_traceback_val (err, routine, name, out);

  /* Create output mosaic */
  out->mosaic = 
    (ObitImageMosaic*)ObitImageMosaicMFCreate (name, order, maxFBW, alpha,
					       alphaRefF, uvdata, err);
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
 * Attach data for a second polarization
 * \param in  The input object
 * \param uvdata2 ObitUV object with info member containng the output image
 *                specifications and all processing parameters.
 * \param err Obit error stack object.
 */
void ObitUVImagerMFAddPol2 (ObitUVImagerMF *in, ObitUV *uvdata2, ObitErr *err)
{
  /*gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
    ObitInfoType type;*/
  olong i;
  gboolean incompatible; 
  ObitImageMosaicMF *mosaic = (ObitImageMosaicMF*)in->mosaic;
  ObitUVDesc *in1Desc, *in2Desc;
  gchar *routine = "ObitUVImagerAddPol2";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  in->noPolImage = 2;  /* Mark object as having 2 poln */

 /* Save uvdata */
  in->uvdata2 = ObitUVRef(uvdata2);
  
 /* Open and close uvdata2 to set descriptor for scratch file */
  ObitUVOpen (in->uvdata2, OBIT_IO_ReadCal, err);
  ObitUVClose (in->uvdata2, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Check compatability between inUV1, inUV2 */
  in1Desc = in->uvdata->myDesc;
  in2Desc = in->uvdata2->myDesc;
  incompatible = in1Desc->nvis!=in2Desc->nvis;
  incompatible = incompatible || (in1Desc->ncorr!=in2Desc->ncorr);
  incompatible = incompatible || (in1Desc->jlocs!=in2Desc->jlocs);
  incompatible = incompatible || (in1Desc->jlocf!=in2Desc->jlocf);
  incompatible = incompatible || (in1Desc->jlocif!=in2Desc->jlocif);
  incompatible = incompatible || (in1Desc->ilocb!=in2Desc->ilocb);
  if (incompatible) {
     Obit_log_error(err, OBIT_Error,"%s uvdata and uvdata2 have incompatible structures",
		   routine);
      return ;
 }

  /* Create scratch uvwork2 if it doesn't exist */
  if (in->uvwork2==NULL) in->uvwork2 = newObitUVScratch (in->uvdata2, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Create mosaic2 - copy of mosaic */
  in->mosaic2 = 
    (ObitImageMosaic*)ObitImageMosaicMFCreate (in->name, in->maxOrder, in->maxFBW, mosaic->alpha,
					       mosaic->alphaRefF, uvdata2, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  /* Copy any autoCen facets */
  for (i=0; i<in->mosaic->numberImages; i++) {
    if (in->mosaic->isAuto[i]>0)
      ObitImageMosaicMFAddField (in->mosaic2, uvdata2, in->mosaic->nx[i], in->mosaic->ny[i],
				 in->mosaic->nplane[i], in->mosaic->RAShift[i], in->mosaic->DecShift[i], 
				 in->mosaic->isAuto[i], err);
  } /* End add autoCen Facets */

  /* Define images */
  ObitImageMosaicMFDefine (in->mosaic2, uvdata2, TRUE, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

 
} /* end ObitUVImagerMFAddPol2 */

/**
 * Apply weighting to uvdata and write to uvwork member
 * \param in  The input object
 *   The following uvdata info items control behavior:
 *   \li "HalfStoke"   OBIT_boo (1,1,1)   If true, half Stokes are passed in uvwork [def F]
 *                                        else I
 *   \li "FullStoke"   OBIT_boo (1,1,1)   If true, all Stokes are passed in uvwork [def F]
 * \param err Obit error stack object.
 */
void ObitUVImagerMFWeight (ObitUVImager *inn, ObitErr *err)
{
  ObitUVImagerMF *in  = (ObitUVImagerMF*)inn;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  gboolean HalfStoke=FALSE, FullStoke=FALSE, Tr=TRUE;
  gchar *HStokes="HALF", *FStokes="FULL", IStokes[5];
  /* List of control parameters on uvwork */
  gchar *controlList[] = 
    {"FOV", "doFull", "NField", "xCells", "yCells", "nx", "ny", 
     "RAShift", "DecShift", "Sources", "Beam", "targBeam", 
     "Catalog",  "CatDisk", "OutlierDist", "OutlierFlux", "OutlierSI", "OutlierSize",
     "nuGrid", "nvGrid", "WtBox", "WtFunc", "UVTaper", "UVITaper", "Robust", "WtPower",
     "RobustIF", "TaperIF", "MFTaper",
     NULL};
  gchar *routine = "ObitUVImagerWeight";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Get Stokes being imaged */
  strcpy (IStokes, "F   "); 
  ObitInfoListGetTest (in->uvdata->info, "Stokes", &type, dim, IStokes);

  /* Want parallel poln? */
  ObitInfoListGetTest (in->uvdata->info, "HalfStoke", &type, dim, &HalfStoke);
  ObitInfoListGetTest (in->uvdata->info, "FullStoke", &type, dim, &FullStoke);

  /* Want RR, LL (XX,YY) data */
  if (HalfStoke) {
    dim[0] = 4;
    ObitInfoListAlwaysPut (in->uvdata->info, "Stokes", OBIT_string, dim, HStokes);
  }

  /* Want all Stokes correlations? */
  if (FullStoke) {
    dim[0] = 4;
    ObitInfoListAlwaysPut (in->uvdata->info, "Stokes", OBIT_string, dim, FStokes);
  }

  /* Open and close uvdata to set descriptor for scratch file */
  ObitUVOpen (in->uvdata, OBIT_IO_ReadCal, err);
  ObitUVClose (in->uvdata, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Create scratch uvwork if it doesn't exist */
  if (in->uvwork==NULL) in->uvwork = newObitUVScratch (in->uvdata, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
 
  /* Copy/calibrate/select uvdata to uvwork */
  in->uvwork = ObitUVCopy (in->uvdata, in->uvwork, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Copy control info to uvwork */
  ObitInfoListCopyList (in->uvdata->info, in->uvwork->info, controlList);

  /* Weight uvwork */
  ObitUVWeightData (in->uvwork, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Image I */
  dim[0] = 4;
  ObitInfoListAlwaysPut (in->uvwork->info, "Stokes", OBIT_string, dim, IStokes);

  /* Restore original Stokes on uvdata */
  ObitInfoListAlwaysPut (in->uvdata->info, "Stokes", OBIT_string, dim, IStokes);
  dim[0] = 1;
  ObitInfoListAlwaysPut (in->uvwork->info, "doCalSelect", OBIT_bool, dim, &Tr);

  /* Repeat for second polarization (U) if given */
  if (in->noPolImage==2) {
    dim[0] = 4;
    ObitInfoListAlwaysPut (in->uvdata->info, "Stokes", OBIT_string, dim, "U   ");
    
    /* Open and close uvdata2 to set descriptor for scratch file */
    ObitUVOpen (in->uvdata2, OBIT_IO_ReadCal, err);
    ObitUVClose (in->uvdata2, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    
    /* Create scratch uvwork if it doesn't exist */
    if (in->uvwork2==NULL) in->uvwork2 = newObitUVScratch (in->uvdata2, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    
    /* Copy/calibrate/select uvdata to uvwork */
    in->uvwork2 = ObitUVCopy (in->uvdata2, in->uvwork2, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    
    /* Copy control info to uvwork */
    ObitInfoListCopyList (in->uvdata2->info, in->uvwork2->info, controlList);
    
    /* Weight uvwork */
    ObitUVWeightData (in->uvwork2, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);

    /* Image I - just kidding */
    dim[0] = 4;
    ObitInfoListAlwaysPut (in->uvwork2->info, "Stokes", OBIT_string, dim, IStokes);
    
    /* Restore original Stokes on uvdata */
    ObitInfoListAlwaysPut (in->uvdata2->info, "Stokes", OBIT_string, dim, IStokes);
    dim[0] = 1;
    ObitInfoListAlwaysPut (in->uvwork2->info, "doCalSelect", OBIT_bool, dim, &Tr);
  } /* end weight 2nd polarization */

} /* end ObitUVImagerMFWeight */

/**
 * Image data in uvwork2 (secondary poln)  if defined, else uvdata2 writing results in mosaic2.
 * If an autoCenter image is imaged, its shifted version is also shifted.
 * \param in        The input object
 * \param field     zero terminated list of field numbers to image, 0=> all
 * \param doWeight  If TRUE do Weighting ov uv data first
 *                  If TRUE then input data is modified.
 * \param doBeam    If True calculate dirty beams first
 * \param doFlatten If TRUE, flatten images when done
 * \param err       Obit error stack object.
 */
void ObitUVImagerMFImage2 (ObitUVImagerMF *in, olong *field, gboolean doWeight, 
			   gboolean doBeam, gboolean doFlatten, ObitErr *err)
{ 
  ObitUVImager *inn = (ObitUVImager*)in;
  ObitUV *data=NULL;
  ObitImage **imageList=NULL;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitUVDesc *UVDesc;
  ObitImageDesc *imageDesc;
  olong i, j, n, fldCnt, ifield, channel=0, nDo, nLeft, nImage, prtLv, *fldNo=NULL;
  olong NumPar, myAuto;
  ofloat sumwts[2];
  ObitImage *theBeam=NULL;
  gboolean *forceBeam=NULL, needBeam, doall, found;
  ObitUVImagerClassInfo *imgClass = (ObitUVImagerClassInfo*)in->ClassInfo;
  gchar        *dataParms[] = {  /* Imaging info */
    "xShift", "yShift",
    NULL
  };
  gchar *routine = "ObitUVImagerMFImage2";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  data = in->uvwork2;
  if (!data) data = in->uvdata2;
  if (!ObitUVIsA(data)) {
    Obit_log_error(err, OBIT_Error,"%s UV data not defined in %s", routine, data->name);
    return;
  }

  /* Copy imaging info if uvwork2 already defined */
  if (in->uvwork2) 
    ObitInfoListCopyList (in->uvdata2->info, in->uvwork2->info, dataParms);

  /* List of need to force making beam */
  forceBeam = g_malloc0(in->mosaic2->numberImages*sizeof(gboolean));
  fldNo     = g_malloc0(in->mosaic2->numberImages*sizeof(olong));

  /* Count number of fields to image */
  nImage = 0;
  for (i=0; i<in->mosaic2->numberImages; i++) {
    if (field[i]==0) break;
    nImage++;
  }
  /* All wanted? */
  doall =  (nImage <= 0);

  /* get prtLv */
  prtLv = 1;
  if (ObitInfoListGetTest(in->mosaic2->info, "prtLv", &type, dim, &prtLv)) 
    err->prtLv = prtLv;  /* Add to err */

  /* Single or multiple images (including beams) */
  if ((!doBeam) && ((nImage==1) || (in->mosaic2->numberImages==1))) {
    /* single - no beam */
    UVDesc    = data->myDesc;
    ifield = MAX (0, (field[0]-1));
    imageDesc = in->mosaic->images[ifield]->myDesc;
    if (UVDesc->crval[UVDesc->jlocs]>0) 
      imageDesc->crval[imageDesc->jlocs] = UVDesc->crval[UVDesc->jlocs];

    /* reset image max/min */
    imageDesc->maxval    = -1.0e20;
    imageDesc->minval    =  1.0e20;

    /* Is this an anto center image? */    
    myAuto = in->mosaic2->isShift[ifield]-1;  /* Corresponding autoCenter */
    if (myAuto<0) myAuto = ifield; /* No */

    /* Need to force beam? */
    theBeam = (ObitImage*)in->mosaic2->images[myAuto]->myBeam;
    forceBeam[0] = (theBeam==NULL) ||
      !ObitInfoListGetTest(theBeam->info, "SUMWTS", &type, dim, (gpointer)&sumwts);
    needBeam = (doBeam || forceBeam[0]);

    /* Image */
    ObitImageUtilMakeImage (data, in->mosaic2->images[myAuto], channel, 
			    needBeam, doWeight, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
 
    /* Made beam? */
    if (doBeam || forceBeam[0]) {
      /* If no beam size given take this one */
      if (in->mosaic2->bmaj==0.0) {
	in->mosaic2->bmaj = in->mosaic2->images[myAuto]->myDesc->beamMaj;
	in->mosaic2->bmin = in->mosaic2->images[myAuto]->myDesc->beamMin;
	in->mosaic2->bpa  = in->mosaic2->images[myAuto]->myDesc->beamPA;
      } else if (in->mosaic2->bmaj>0.0) { /* beam forced */
	in->mosaic2->images[myAuto]->myDesc->beamMaj = in->mosaic2->bmaj;
	in->mosaic2->images[myAuto]->myDesc->beamMin = in->mosaic2->bmin;
	in->mosaic2->images[myAuto]->myDesc->beamPA  = in->mosaic2->bpa;
	/* Tell if field 1 */
	if (myAuto==0) {
	  Obit_log_error(err, OBIT_InfoErr, 
			 "Using Beamsize %f x %f asec PA=%f",
			 in->mosaic2->bmaj*3600.0, in->mosaic2->bmin*3600.0, 
			 in->mosaic2->bpa);
	}
      }
    }

    goto shifty;
  } /* end single */

  /* Multiple (including beams) - do in parallel */

  NumPar = imgClass->ObitUVImagerGetNumPar(inn, doBeam, err); /* How many to do? */

  /* Get list of images */
  imageList = g_malloc0(in->mosaic2->numberImages*sizeof(ObitImage*));
  if (doall) n = in->mosaic2->numberImages;
  else n = nImage;
  fldCnt = 0;
  for (i=0; i<n; i++) {
    if (doall) ifield = i;
    else ifield = field[i]-1;
    if (in->mosaic2->isShift[ifield]<=0) {  /* Special handling for shifted imaged */
      imageList[fldCnt] = in->mosaic2->images[ifield];
      fldNo[fldCnt++]   = ifield;                 /* Number (0-rel) in mosaic */
    } else if (!doall) { /* Make autoCenter version instead */
      myAuto = in->mosaic2->isShift[ifield]-1;  /* Corresponding autoCenter */
      /* See if already in list, if so ignore */
      found = FALSE;
      for (j=0; j<n; j++) found = found || ((myAuto+1)==field[j]);
      if (!found) {
	imageList[fldCnt] = in->mosaic2->images[myAuto];
	fldNo[fldCnt++]   = myAuto;
	field[i] = myAuto+1;
      }
    }
  }

  /* Loop over fields initializing */
  for (i=0; i<fldCnt; i++) {

    /* Set Stokes */
    UVDesc    = data->myDesc;
    imageDesc = imageList[i]->myDesc;
    if (UVDesc->crval[UVDesc->jlocs]>0) 
      imageDesc->crval[imageDesc->jlocs] = UVDesc->crval[UVDesc->jlocs];

    /* reset image max/min */
    imageDesc->maxval    = -1.0e20;
    imageDesc->minval    =  1.0e20;
 
    /* Is this an anto center image? */    
    myAuto = in->mosaic2->isShift[i]-1;  /* Corresponding autoCenter */
    if (myAuto<0) myAuto = i; /* No */

    /* Need to force beam? */
    theBeam = (ObitImage*)imageList[myAuto]->myBeam;
    forceBeam[i] = (theBeam==NULL) ||
      !ObitInfoListGetTest(theBeam->info, "SUMWTS", &type, dim, (gpointer)&sumwts);

  } /* end loop initializing */

  /* Image - don't do more than NumPar at a time */
  nLeft = fldCnt;
  ifield = 0;
  while (nLeft>0) {
    nDo = MIN (NumPar, nLeft);
    ObitImageUtilMakeImagePar (data, nDo, &imageList[ifield],
			       doBeam, doWeight, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    ifield += nDo;
    nLeft  -= nDo;
    if (prtLv>1) ObitErrLog(err);  /* Progress Report */
    else ObitErrClear(err);
  } /* End loop gridding */

  /* Loop over fields finalizing */
  for (i=0; i<fldCnt; i++) {
    /* If it made a beam check the beam size */
    if (doBeam || forceBeam[i]) {
      /* If no beam size given take this one */
      if (in->mosaic->bmaj==0.0) {
	in->mosaic2->bmaj = imageList[i]->myDesc->beamMaj;
	in->mosaic2->bmin = imageList[i]->myDesc->beamMin;
	in->mosaic2->bpa  = imageList[i]->myDesc->beamPA;
      } else if (in->mosaic2->bmaj>0.0) { /* beam forced */
	imageList[i]->myDesc->beamMaj = in->mosaic2->bmaj;
	imageList[i]->myDesc->beamMin = in->mosaic2->bmin;
	imageList[i]->myDesc->beamPA  = in->mosaic2->bpa;
	/* Tell if field 1 */
	if ((i==0) && (fldNo[0]==1)) {
	  Obit_log_error(err, OBIT_InfoErr, 
			 "Using Beamsize %f x %f asec PA=%f",
			 in->mosaic2->bmaj*3600.0, in->mosaic2->bmin*3600.0, 
			 in->mosaic2->bpa);
	}
      }
    }
  } /* End loop over fields finalizing*/
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Make any 2D shifted images */
 shifty:
  imgClass->ObitUVImagerShifty (inn, field, doall, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Need to flatten? */
  if (doFlatten) ObitUVImagerFlatten (inn, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Cleanup */
  if (forceBeam) g_free(forceBeam);
  if (fldNo) g_free(fldNo);
  if (imageList) g_free(imageList);
} /* end ObitUVImagerMFImage2 */

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
  ObitImageMosaic *mosaic = NULL;
  olong i, j, ofield, ifield;
  ofloat shift[2];
  gboolean found;
  gchar *routine = "ObitUVImagerMFShifty";

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* which mosaic? */
  if (in->whichPol==1) mosaic = in->mosaic;
  else                 mosaic = in->mosaic2;
  if (!ObitImageMosaicIsA(mosaic)) {
    Obit_log_error(err, OBIT_Error,"%s ImageMosaic not defined in %s", 
		   routine, in->name);
    return;
  }

  /* Loop making any (2D) shifted images */
  for (i=0; i<mosaic->numberImages; i++) {
    /* Only needed for 2D */
    if (mosaic->images[i]->myDesc->do3D) continue;
    ofield = i;
    ifield = mosaic->isShift[i]-1;
    if (ifield >= 0) {
      /* See if input image just made - named in field */
      found = FALSE;
      for (j=0; j<mosaic->numberImages; j++) {
	if (field[j]==0) break;
	if (field[j]==(ifield+1)) found = TRUE;
      }
      if (found || doall) {
	shift[0] = mosaic->images[ofield]->myDesc->crpix[0] - mosaic->images[ifield]->myDesc->crpix[0];
	shift[1] = mosaic->images[ofield]->myDesc->crpix[1] - mosaic->images[ifield]->myDesc->crpix[1];
	ObitImageUtilMFShift (mosaic->images[ifield], mosaic->images[ofield], shift, err);
	/* Beam - no shift */
	shift[0] = shift[1] = 0.0;
	ObitImageUtilMFShift ((ObitImage*)mosaic->images[ifield]->myBeam, 
			      (ObitImage*)mosaic->images[ofield]->myBeam, shift, err);
      }
    } /* end if shifted field */
  } /* End loop making shifted images */
  if (err->error) Obit_traceback_msg (err, routine, in->name);

} /* end ObitUVImageMFShifty */

/**
 * return ImageMosaic member using whichPol
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

  if (in->whichPol==1) return ObitImageMosaicRef(in->mosaic);
  else                 return ObitImageMosaicRef(in->mosaic2);
} /* end ObitUVImagerMFGetMosaic */

/**
 * return secondary ImageMosaic member (mosaic2)
 * \param in  The input object
 * \param err Obit error stack object.
 * \return reference to ImageMosaic.
 */
ObitImageMosaic* ObitUVImagerMFGetMosaic2 (ObitUVImagerMF *in, ObitErr *err)
{ 
  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return NULL;
  g_assert (ObitIsA(in, &myClassInfo));

  return ObitImageMosaicMFRef(in->mosaic2);
} /* end ObitUVImagerMFGetMosaic2 */

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
     "RobustIF", "TaperIF", "MFTaper",
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
 * Target memory usage is 0.75 GByte if 32 bit, 3 GByte if 64.
 * \param inn     Object of interest.
 * \return the number of parallel images.
 */
olong ObitUVImagerMFGetNumPar (ObitUVImager *inn, gboolean doBeam, ObitErr *err)
{
  ObitUVImagerMF *in  = (ObitUVImagerMF*)inn;
  olong out=8, nSpec;
  odouble lenVis, numVis, imSize, bufSize, tSize, mSize;
  ObitImage *beam;
  /*long int avphys_pages;*/
  long int phys_pages;
  int pagesize;
  gchar *routine="ObitUVImagerMFGetNumPar";

  if (err->error) return out;
  g_assert(ObitUVImagerMFIsA(in));

#ifdef __APPLE__
  int mib[2];
  int64_t physical_memory;
  size_t length;

  // Get the Physical memory size
  mib[0] = CTL_HW;
  mib[1] = HW_MEMSIZE;
  length = sizeof(physical_memory);
  sysctl(mib, 2, &physical_memory, &length, NULL, 0);  
  mSize = physical_memory;
#else
   /* Inquire as to the amount of available pages of memory */
  /*avphys_pages = get_avphys_pages();*/
  phys_pages   = get_phys_pages();
  pagesize     = getpagesize();
  /* How much to ask for (64-bit) - up to 1/5 total */
  mSize = 0.2*phys_pages*pagesize;
#endif

  /* How big are things? assumes 2nd polarization the same as first */
  numVis = (odouble)ObitImageUtilBufSize (in->uvdata);  /* Size of buffer */
  lenVis = (odouble)in->uvdata->myDesc->lrec;
  imSize = in->mosaic->images[0]->myDesc->inaxes[0] * 
    in->mosaic->images[0]->myDesc->inaxes[1];  /* Image plane size */
  if (doBeam) {
    /* Beam size */
    beam = (ObitImage*)in->mosaic->images[0]->myBeam;
    if (beam!=NULL) imSize += beam->myDesc->inaxes[0] * beam->myDesc->inaxes[1];
    else imSize *= 5;  /* Assume 4X as large (plus map) */
  }

  nSpec = ((ObitImageMF*)in->mosaic->images[0])->nSpec;
  bufSize = numVis*lenVis + imSize*nSpec;  /* Approx memory (words) per parallel image */
  bufSize *= sizeof(ofloat);               /* to bytes */

  if (sizeof(olong*)==4)      tSize = 0.75e9;  /* Use sizeof a pointer type to get 32/64 bit */
  else if (sizeof(olong*)==8) tSize =  mSize;  /*3.0e9; */
  else                        tSize = 1.0e9;  /* Shouldn't happen */
  out = tSize / bufSize;  /* How many fit in a tSize? */

  /* Better be at least 1 */
  Obit_retval_if_fail((out>=1), err, out,
		      "%s: Insufficient memory to make images", routine);
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
  theClass->ObitUVImagerWeight = (ObitUVImagerWeightFP)ObitUVImagerMFWeight;
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
  in->uvdata2 = NULL;
  in->uvwork2 = NULL;
  in->mosaic2 = NULL;
  in->noPolImage = 1;
  in->whichPol   = 1;
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
  in->uvdata2 = ObitUVUnref(in->uvdata2);
  in->uvwork2 = ObitUVUnref(in->uvwork2);
  in->mosaic2 = ObitImageMosaicUnref(in->mosaic2);
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitUVImagerMFClear */

