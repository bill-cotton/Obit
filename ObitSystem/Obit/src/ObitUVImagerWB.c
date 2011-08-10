/* $Id$       */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2010,2011                                          */
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

#include "ObitUVImager.h"
#include "ObitUVImagerWB.h"
#include "ObitImageWB.h"
#include "ObitUVWeight.h"
#include "ObitImageUtil.h"
#include "ObitImageMosaicWB.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitUVImagerWB.c
 * ObitUVImagerWB class function definitions.
 * This class is derived from the Obit base class.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitUVImagerWB";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitUVImagerGetClass;

/**
 * ClassInfo structure ObitUVImagerWBClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitUVImagerWBClassInfo myClassInfo = {FALSE};

/*--------------- File Global Variables  ----------------*/


/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitUVImagerWBInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitUVImagerWBClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitUVImagerWBClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitUVImagerWB* newObitUVImagerWB (gchar* name)
{
  ObitUVImagerWB* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitUVImagerWBClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitUVImagerWB));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitUVImagerWBInit((gpointer)out);

 return out;
} /* end newObitUVImagerWB */

/**
 * Initializes from ObitInfoList.
 * Initializes class if needed on first call.
 * \param out     the new object.to be initialized
 * \param prefix  If NonNull, string to be added to beginning of inList entry name
 *                "xxx" in the following
 * \param inList  InfoList to extract object information from 
 *      \li "xxxClassType" string UVImager type, "WB" for this class
 *      \li "xxxUVData" prefix for uvdata member, entry with value "None" => doesn't exist
 *      \li "xxxUVWork" prefix for uvwork member, entry with value "None" => doesn't exist
 *      \li "xxxMosaic" prefix for mosaic member, entry with value "None" => doesn't exist
 * \param err     ObitErr for reporting errors.
 */
void ObitUVImagerWBFromInfo (ObitUVImager *out, gchar *prefix, ObitInfoList *inList, 
			     ObitErr *err)
{ 
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *keyword=NULL, *value=NULL;
  gboolean missing;
  gchar *Type = "WB";
  gchar *routine = "ObitUVImagerWBFromInfo";
  
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitUVImagerWBClassInit();

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(out, &myClassInfo));

  /* check class type */
  missing = ObitInfoListGetP(inList, keyword, &type, dim, (gpointer*)&value);
  if ((missing) || (type!=OBIT_string) || (!strncmp(Type,value,dim[0]))) {
    Obit_log_error(err, OBIT_Error,"%s Wrong class type %s!=%s", routine, value, Type);
    return;
  }

} /* end ObitUVImagerWBFromInfo */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitUVImagerWBGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitUVImagerWBClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitUVImagerWBGetClass */

/**
 * Make a deep copy of an ObitUVImagerWB.
 * \param int   The object to copy
 * \param outt  An existing object pointer for output or NULL if none exists.
 * \param err  Obit error stack object.
 * \return pointer to the new object.
 */
ObitUVImager* ObitUVImagerWBCopy  (ObitUVImager *inn, ObitUVImager *outt, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  ObitUVImagerWB *in  = (ObitUVImagerWB*)inn;
  ObitUVImagerWB *out = (ObitUVImagerWB*)outt;
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
    out = newObitUVImagerWB(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (inn, outt, err);

  /*  copy this class - just pointers */

  return (ObitUVImager*)out;
} /* end ObitUVImagerWBCopy */

/**
 * Make a copy of a object but do not copy the actual data
 * This is useful to create an UVImagerWB similar to the input one.
 * \param inn  The object to copy
 * \param outt  An existing object pointer for output, must be defined.
 * \param err  Obit error stack object.
 */
void ObitUVImagerWBClone  (ObitUVImager *inn, ObitUVImager *outt, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  ObitUVImagerWB *in  = ( ObitUVImagerWB*)inn;
  ObitUVImagerWB *out = ( ObitUVImagerWB*)outt;

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

} /* end ObitUVImagerWBClone */

/**
 * Creates an ObitUVImagerWB given an ObitUV with control information.
 * The output ImageMosaic member is created
 * \param name   An optional name for the object.
 * \param order  Spectral imaging order,0=flux,1=si, 2=curve
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
ObitUVImagerWB* ObitUVImagerWBCreate (gchar* name, olong order, ObitUV *uvdata, 
				      ObitErr *err)
{
  ObitUVImagerWB* out=NULL;
  gchar *routine = "ObitUVImagerWBCreate";

  /* Error checks */
  if (err->error) return out;
  g_assert(ObitUVIsA(uvdata));

  /* Create basic structure */
  out = newObitUVImagerWB (name);

  /* Save uvdata */
  out->uvdata = ObitUVRef(uvdata);
  out->norder = order;  /* Save order */

  /* Create output mosaic */
  out->mosaic = (ObitImageMosaic*)ObitImageMosaicWBCreate (name, order, uvdata, err);
  if (err->error) Obit_traceback_val (err, routine, name, out);

  /* Define images */
  ObitImageMosaicWBDefine (out->mosaic, uvdata, TRUE, err);
  if (err->error) Obit_traceback_val (err, routine, name, out);

  return out;
} /* end ObitUVImagerWBCreate */

/**
 * Creates an ObitUVImagerWB given an ObitUV with control information
 * and a previously existing ImageMosaic
 * The output ImageMosaic member is created
 * \param name   An optional name for the object.
 * \param order  Spectral imaging order,0=flux,1=si, 2=curve
 * \param uvdata ObitUV object with info member containng the output image
 *               specifications and all processing parameters.
 * \param mosaic ImageMosaicWB (as ImageMosaic) to use
 * \param err Obit error stack object.
 * \return the new object.
 */
ObitUVImagerWB* ObitUVImagerWBCreate2 (gchar* name, olong order, ObitUV *uvdata, 
				       ObitImageMosaic *mosaic, ObitErr *err)
{
  ObitUVImagerWB* out=NULL;
  gchar *routine = "ObitUVImagerWBCreate2 ";

  /* Error checks */
  if (err->error) return out;
  g_assert(ObitUVIsA(uvdata));

  /* Check input types */
  Obit_retval_if_fail((ObitImageMosaicWBIsA((ObitImageMosaicWB*)mosaic) && 
		       (ObitImageWBIsA((ObitImageWB*)mosaic->images[0]))), err, 
		       out,
		       "%s: Image mosaic or images not WB", routine);
		      
  /* Create basic structure */
  out = newObitUVImagerWB (name);

  /* Save uvdata */
  out->uvdata = ObitUVRef(uvdata);
  out->norder = order;  /* Save order */

  /* Save mosaic */
  out->mosaic = ObitImageMosaicRef(mosaic);

  return out;
} /* end ObitUVImagerWBCreate2 */

/**
 * Image data in uvwork if defined, else uvdata writing results in mosaic.
 * If an autoCenter image is imaged, its shifted version is also shifted.
 * \param inn       The input object
 * \param field     zero terminated list of field numbers to image, 0=> all
 * \param doWeight  If TRUE do Weighting ov uv data first
 *                  If TRUE then input data is modified.
 * \param doBeam    If True calculate dirst beams first
 * \param doFlatten If TRUE, flatten images when done
 * \param err       Obit error stack object.
 */
void ObitUVImagerWBImage (ObitUVImager *inn,  olong *field, gboolean doWeight, 
			  gboolean doBeam, gboolean doFlatten, ObitErr *err)
{ 
  ObitUV *data=NULL;
  ObitUVImagerWB *in  = ( ObitUVImagerWB*)inn;
  ObitImage **imageList=NULL;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitUVDesc *UVDesc;
  ObitImageDesc *imageDesc;
  olong i, j, n, fldCnt, ifield, ofield, channel=0, nDo, nLeft, nImage, prtLv, *fldNo=NULL;
  olong NumPar;
  ofloat sumwts, shift[2];
  ObitImage *theBeam=NULL;
  gboolean *forceBeam=NULL, needBeam, doall, found;
  gchar *routine = "ObitUVImagerWBImage";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  data = in->uvwork;
  if (!data) data = in->uvdata;
  if (!ObitUVIsA(data)) {
    Obit_log_error(err, OBIT_Error,"%s UV data not defined in %s", routine, data->name);
    return;
  }

  /* List of need to force making beam */
  forceBeam = g_malloc0(in->mosaic->numberImages*sizeof(gboolean));
  fldNo     = g_malloc0(in->mosaic->numberImages*sizeof(olong));

  /* Count number of fields to image */
  nImage = 0;
  for (i=0; i<in->mosaic->numberImages; i++) {
    if (field[i]==0) break;
    nImage++;
  }
  /* All wanted? */
  doall =  (nImage <= 0);

  /* get prtLv */
  prtLv = 1;
  ObitInfoListGetTest(in->mosaic->info, "prtLv", &type, dim, &prtLv);

  /* Single or multiple images (including beams) */
  if ((!doBeam) && ((nImage==1) || (in->mosaic->numberImages==1))) {
    /* single - no beam */
    UVDesc    = data->myDesc;
    ifield = MAX (0, (field[0]-1));
    imageDesc = in->mosaic->images[ifield]->myDesc;
    if (UVDesc->crval[UVDesc->jlocs]>0) 
      imageDesc->crval[imageDesc->jlocs] = UVDesc->crval[UVDesc->jlocs];

    /* reset image max/min */
    imageDesc->maxval    = -1.0e20;
    imageDesc->minval    =  1.0e20;

    /* Need to force beam? */
    theBeam = (ObitImage*)in->mosaic->images[ifield]->myBeam;
    forceBeam[0] = (theBeam==NULL) ||
      !ObitInfoListGetTest(theBeam->info, "SUMWTS", &type, dim, (gpointer)&sumwts);
    needBeam = (doBeam || forceBeam[0]);

    /* Image */
    ObitImageUtilMakeImage (data, in->mosaic->images[ifield], channel, 
			    needBeam, doWeight, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
 
    /* Made beam? */
    if (doBeam || forceBeam[0]) {
      /* If no beam size given take this one */
      if (in->mosaic->bmaj==0.0) {
	in->mosaic->bmaj = in->mosaic->images[ifield]->myDesc->beamMaj;
	in->mosaic->bmin = in->mosaic->images[ifield]->myDesc->beamMin;
	in->mosaic->bpa  = in->mosaic->images[ifield]->myDesc->beamPA;
      } else if (in->mosaic->bmaj>0.0) { /* beam forced */
	in->mosaic->images[ifield]->myDesc->beamMaj = in->mosaic->bmaj;
	in->mosaic->images[ifield]->myDesc->beamMin = in->mosaic->bmin;
	in->mosaic->images[ifield]->myDesc->beamPA  = in->mosaic->bpa;
	/* Tell if field 1 */
	if (ifield==0) {
	  Obit_log_error(err, OBIT_InfoErr, 
			 "Using Beamsize %f x %f asec PA=%f",
			 in->mosaic->bmaj*3600.0, in->mosaic->bmin*3600.0, 
			 in->mosaic->bpa);
	}
      }
    }

    goto shifty;
  } /* end single */

  /* Multiple (including beams) - do in parallel */

  NumPar = ObitUVImagerWBGetNumPar(inn, err); /* How many to do? */

  /* Get list of images */
  imageList = g_malloc0(in->mosaic->numberImages*sizeof(ObitImage*));
  if (doall) n = in->mosaic->numberImages;
  else n = nImage;
  fldCnt = 0;
  for (i=0; i<n; i++) {
    if (doall) ifield = i;
    else ifield = field[i]-1;
    if (in->mosaic->isShift[ifield]<0) {  /* Special handling for shifted imaged */
      imageList[fldCnt] = in->mosaic->images[ifield];
      fldNo[fldCnt++]   = ifield;                 /* Number (0-rel) in mosaic */
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
 
    /* Need to force beam? */
    theBeam = (ObitImage*)imageList[i]->myBeam;
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
	in->mosaic->bmaj = imageList[i]->myDesc->beamMaj;
	in->mosaic->bmin = imageList[i]->myDesc->beamMin;
	in->mosaic->bpa  = imageList[i]->myDesc->beamPA;
      } else if (in->mosaic->bmaj>0.0) { /* beam forced */
	imageList[i]->myDesc->beamMaj = in->mosaic->bmaj;
	imageList[i]->myDesc->beamMin = in->mosaic->bmin;
	imageList[i]->myDesc->beamPA  = in->mosaic->bpa;
	/* Tell if field 1 */
	if ((i==0) && (field[0]==1)) {
	  Obit_log_error(err, OBIT_InfoErr, 
			 "Using Beamsize %f x %f asec PA=%f",
			 in->mosaic->bmaj*3600.0, in->mosaic->bmin*3600.0, 
			 in->mosaic->bpa);
	}
      }
    }
  } /* End loop over fields finalizing*/
  if (err->error) Obit_traceback_msg (err, routine, in->name);

 shifty:
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
	ObitImageUtilShift (in->mosaic->images[ifield], in->mosaic->images[ofield], shift, err);
	/* Beam - no shift */
	shift[0] = shift[1] = 0.0;
	ObitImageUtilShift ((ObitImage*)in->mosaic->images[ifield]->myBeam, 
			    (ObitImage*)in->mosaic->images[ofield]->myBeam, shift, err);
      }
    } /* end if shifted field */
  } /* End loop making shifted images */
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Need to flatten? */
  if (doFlatten) ObitUVImagerFlatten ((ObitUVImager*)in, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Cleanup */
  if (forceBeam) g_free(forceBeam);
  if (fldNo) g_free(fldNo);
  if (imageList) g_free(imageList);
} /* end ObitUVImagerWBImage */

/**
 * return ImageMosaic member
 * \param inn  The input object
 * \param err  Obit error stack object.
 * \return reference to ImageMosaic.
 */
ObitImageMosaic* ObitUVImagerWBGetMosaic (ObitUVImager *inn, ObitErr *err)
{ 
  ObitUVImagerWB *in  = ( ObitUVImagerWB*)inn;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return NULL;
  g_assert (ObitIsA(in, &myClassInfo));

  return ObitImageMosaicRef(in->mosaic);
} /* end ObitUVImagerWBGetMosaic */

/**
 * Convert structure information to entries in an ObitInfoList
 * \param inn     Object of interest.
 * \param prefix  If NonNull, string to be added to beginning of outList entry name
 *                "xxx" in the following
 * \param outList InfoList to write entries into
 *      \li "xxxClassType" string UVImagerWB type, "Base" for base class
 *      \li "xxxUVData" prefix for uvdata member, entry with value "None" => doesn't exist
 *      \li "xxxUVWork" prefix for uvwork member, entry with value "None" => doesn't exist
 *      \li "xxxMosaic" prefix for mosaic member, entry with value "None" => doesn't exist
 *      \li various weighting parameters
 * \param err     ObitErr for reporting errors.
 */
void ObitUVImagerWBGetInfo (ObitUVImager *inn, gchar *prefix, ObitInfoList *outList, 
			    ObitErr *err)
{ 
  ObitUVImagerWB *in  = ( ObitUVImagerWB*)inn;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *keyword=NULL, *None = "None", *OK="OK", *Type="WB";
  olong i;
  gpointer listPnt;
  gchar *parm[] = 
    {"do3D", "FOV", "doFull", "NField", "xCells", "yCells", "nx", "ny", 
     "RAShift", "DecShift", "Sources", 
     "Catalog", "CatDisk", "OutlierDist", "OutlierFlux", "OutlierSI", "OutlierSize",
     "nuGrid", "nvGrid", "WtBox", "WtFunc", "UVTaper", "Robust", "WtPower",
     NULL};
  gchar *routine = "ObitUVImagerWBGetInfo";

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

} /* end ObitUVImagerWBGetInfo */

/**
 * Get number of parallel images
 * Target memory usage is 0.75 GByte if 32 bit, 3 GByte if 64.
 * \param inn     Object of interest.
 * \return the number of parallel images.
 */
olong ObitUVImagerWBGetNumPar (ObitUVImager *inn, ObitErr *err)
{
  ObitUVImagerWB *in  = (ObitUVImagerWB*)inn;
  olong out=8;
  odouble lenVis, numVis, imSize, bufSize, tSize;
  ObitImage *beam;
  gchar *routine="ObitUVImagerWBGetNumPar";

  if (err->error) return out;
  g_assert(ObitUVImagerWBIsA(in));

  /* How big are things? */
  numVis = (odouble)ObitImageUtilBufSize (in->uvdata);  /* Size of buffer */
  lenVis = (odouble)in->uvdata->myDesc->lrec;
  imSize =  in->mosaic->images[0]->myDesc->inaxes[0] * 
    in->mosaic->images[0]->myDesc->inaxes[1];  /* Image plane size */
  /* Beam size */
  beam = (ObitImage*)in->mosaic->images[0]->myBeam;
  if (beam!=NULL) imSize += beam->myDesc->inaxes[0] * beam->myDesc->inaxes[1];
  else imSize *= 4;  /* Assume 4X as large */

  bufSize = numVis*lenVis + imSize;  /* Approx memory (words) per parallel image */
  bufSize *= sizeof(ofloat);         /* to bytes */

  if (sizeof(olong*)==4)      tSize = 0.75e9;  /* Use sizeof a pointer type to get 32/64 bit */
  else if (sizeof(olong*)==8) tSize = 3.0e9;
  else                        tSize = 1.0e9;  /* Shouldn't happen */
  out = tSize / bufSize;  /* How many fit in a tSize? */

  /* Better be at least 1 */
  Obit_retval_if_fail((out>=1), err, out,
		      "%s: Insufficient memory to make images", routine);

  return out;
} /*  end ObitUVImagerWBGetNumPar */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitUVImagerWBClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitUVImagerWBClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitUVImagerWBClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitUVImagerWBClassInfoDefFn (gpointer inClass)
{
  ObitUVImagerWBClassInfo *theClass = (ObitUVImagerWBClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit      = (ObitClassInitFP)ObitUVImagerWBClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitUVImagerWBClassInfoDefFn;
  theClass->ObitGetClass       = (ObitGetClassFP)ObitUVImagerWBGetClass;
  theClass->ObitUVImagerCreate = (ObitUVImagerCreateFP)ObitUVImagerCreate;
  theClass->ObitUVImagerCreate2= (ObitUVImagerCreate2FP)ObitUVImagerWBCreate2;
  theClass->ObitUVImagerImage  = (ObitUVImagerImageFP)ObitUVImagerWBImage;
  theClass->ObitUVImagerGetInfo= (ObitUVImagerGetInfoFP)ObitUVImagerWBGetInfo;
  theClass->ObitUVImagerGetNumPar= (ObitUVImagerGetNumParFP)ObitUVImagerWBGetNumPar;

} /* end ObitUVImagerWBClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitUVImagerWBInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitUVImagerWB *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */

} /* end ObitUVImagerWBInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitUVImagerWB* cast to an Obit*.
 */
void ObitUVImagerWBClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitUVImagerWB *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitUVImagerWBClear */

