/* $Id$        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2005-2014                                          */
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
#include "ObitUVWeight.h"
#include "ObitImageUtil.h"
#include "ObitUVImagerIon.h"
#include "ObitUVImagerSquint.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitUVImager.c
 * ObitUVImager class function definitions.
 * This class is derived from the Obit base class.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitUVImager";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo structure ObitUVImagerClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitUVImagerClassInfo myClassInfo = {FALSE};

/*--------------- File Global Variables  ----------------*/


/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitUVImagerInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitUVImagerClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitUVImagerClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitUVImager* newObitUVImager (gchar* name)
{
  ObitUVImager* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitUVImagerClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitUVImager));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitUVImagerInit((gpointer)out);

 return out;
} /* end newObitUVImager */

/**
 * Constructor from ObitInfoList.
 * Initializes class if needed on first call.
 * Also works for derived classes.
 * \param prefix  If NonNull, string to be added to beginning of inList entry name
 *                "xxx" in the following
 * \param inList  InfoList to extract object information from 
 *      \li "xxxClassType" string UVImager type, "Base" for base class
 *      \li "xxxUVData" prefix for uvdata member, entry with value "None" => doesn't exist
 *      \li "xxxUVWork" prefix for uvwork member, entry with value "None" => doesn't exist
 *      \li "xxxMosaic" prefix for mosaic member, entry with value "None" => doesn't exist
 *      \li various weighting parameters
 * \param err     ObitErr for reporting errors.
 * \return the new object.
 */
ObitUVImager* ObitUVImagerFromInfo (gchar *prefix, ObitInfoList *inList, 
				    ObitErr *err)
{ 
  ObitUVImager *out = NULL;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *keyword=NULL, *None = "None", *value=NULL, *classType=NULL;
  ObitUV *uvdata=NULL, *uvwork=NULL;
  ObitImageMosaic *mosaic=NULL;
  olong classCnt, i;
  gboolean missing;
  gpointer listPnt;
  gchar *parm[] = 
    {"do3D", "FOV", "doFull", "NField", "xCells", "yCells", "nx", "ny", 
     "RAShift", "DecShift", "Sources", 
     "Catalog", "CatDisk", "OutlierDist", "OutlierFlux", "OutlierSI", "OutlierSize",
     "nuGrid", "nvGrid", "WtBox", "WtFunc", "UVTaper", "Robust", "WtPower",
     "RobustIF", "TaperIF", "MFTaper", "doGPU",
     NULL};
  gchar ctemp[50];
  gchar *routine = "ObitUVImagerFromInfo";

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitUVImagerClassInit();

  /* error checks */
  if (err->error) return out;

  /* check class type */
  if (prefix) keyword = g_strconcat (prefix, "ClassType", NULL);
  else        keyword = g_strdup("ClassType");
  missing = ObitInfoListGetP(inList, keyword, &type, dim, (gpointer*)&classType);
  if ((missing) || (type!=OBIT_string)) {
    Obit_log_error(err, OBIT_Error,"%s No class type", routine);
    return out;
  }
  classCnt = dim[0]; /* How many characters in name */
  g_free(keyword);

  /* uv data */
  if (prefix) keyword = g_strconcat (prefix, "UVData", NULL);
  else        keyword = g_strdup("UVData");
  missing = ObitInfoListGetP(inList, keyword, &type, dim, (gpointer*)&value);
  /* Does it exist? */
  if ((missing) || (type!=OBIT_string) || (!strncmp(None,value,dim[0]))) {
    Obit_log_error(err, OBIT_Error,"%s UV data not defined in %s", routine, keyword);
    return out;
  } else { /* exists*/
    uvdata = (ObitUV*)ObitDataFromFileInfo(keyword, inList, err);
    if (err->error) Obit_traceback_val (err, routine, keyword, out);
  }
  g_free(keyword);

  /* uv data */
  if (prefix) keyword = g_strconcat (prefix, "UVWork", NULL);
  else        keyword = g_strdup("UVWork");
  missing = ObitInfoListGetP(inList, keyword, &type, dim, (gpointer*)&value);
  /* Does it exist? */
  if ((missing) || (type!=OBIT_string) || (!strncmp(None,value,dim[0]))) {
    uvwork = NULL;
  } else { /* exists*/
    uvwork = (ObitUV*)ObitDataFromFileInfo(keyword, inList, err);
    if (err->error) Obit_traceback_val (err, routine, keyword, out);
  }
  g_free(keyword);

  /* ImageMosaic */
  if (prefix) keyword = g_strconcat (prefix, "Mosaic", NULL);
  else        keyword = g_strdup("Mosaic");
  missing = ObitInfoListGetP(inList, keyword, &type, dim, (gpointer*)&value);
  /* Does it exist? */
  if ((missing) || (type!=OBIT_string) || (!strncmp(None,value,dim[0]))) {
    Obit_log_error(err, OBIT_Error,"%s ImageMosaic not defined in %s", 
		   routine, keyword);
    return out;
  } else { /* exists*/
    mosaic = (ObitImageMosaic*)ObitImageMosaicFromInfo(keyword, inList, err);
    if (err->error) Obit_traceback_val (err, routine, keyword, out);
  }
  g_free(keyword);

  /* Create output - by type */
  if (!strncmp("Base", classType, classCnt)) {
    out = ObitUVImagerCreate2(prefix, uvdata, mosaic, err);
  } else if (!strncmp("Squint", classType, classCnt)) {
    out = (ObitUVImager*)ObitUVImagerSquintCreate2(prefix, uvdata, mosaic, err);
    ObitUVImagerSquintFromInfo(out, prefix, inList, err);
  } else if (!strncmp("Ion", classType, classCnt)) {
    out = (ObitUVImager*)ObitUVImagerIonCreate2(prefix, uvdata, mosaic, err);
    ObitUVImagerIonFromInfo(out, prefix, inList, err);
 } else {  /* Assume base and hope for the best */
    out = ObitUVImagerCreate2(prefix, uvdata, mosaic, err);
    /* Note problem in log */
    strncpy (ctemp, classType, MIN (48,classCnt)); ctemp[MIN (49,classCnt+1)] = 0;
    Obit_log_error(err, OBIT_InfoWarn, "%s: Unknown type %s using base class",
		   routine, ctemp);
  }

  /* Weighting/Imaging parameters */
  i = 0;
  while (parm[i]) {
    if (prefix) keyword = g_strconcat (prefix, parm[i], NULL);
    else        keyword = g_strdup(parm[i]);
    if (ObitInfoListGetP(inList, parm[i], &type, dim, (gpointer*)&listPnt)) {
      ObitInfoListAlwaysPut(uvdata->info, keyword, type, dim, listPnt);
    }
    i++;
    g_free(keyword);
  }

  if (err->error) Obit_traceback_val (err, routine, "Output", out);
  out->uvwork = ObitUVRef(uvwork);

  /* cleanup */
  uvdata = ObitUVUnref(uvdata);
  uvwork = ObitUVUnref(uvwork);
  mosaic = ObitImageMosaicUnref(mosaic);

  return out;
} /* end ObitUVImagerFromInfo */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitUVImagerGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitUVImagerClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitUVImagerGetClass */

/**
 * Make a deep copy of an ObitUVImager.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitUVImager* ObitUVImagerCopy  (ObitUVImager *in, ObitUVImager *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  gboolean oldExist;
  gchar *outName;

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
    out = newObitUVImager(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class - just pointers */
  out->uvdata = ObitUVUnref(out->uvdata);
  out->uvwork = ObitUVUnref(out->uvwork);
  out->mosaic = ObitImageMosaicUnref(out->mosaic);
  out->uvdata = ObitUVRef(in->uvdata);
  out->uvwork = ObitUVRef(in->uvwork);
  out->mosaic = ObitImageMosaicRef(in->mosaic);

  return out;
} /* end ObitUVImagerCopy */

/**
 * Make a copy of a object but do not copy the actual data
 * This is useful to create an UVImager similar to the input one.
 * \param in  The object to copy
 * \param out An existing object pointer for output, must be defined.
 * \param err Obit error stack object.
 */
void ObitUVImagerClone  (ObitUVImager *in, ObitUVImager *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;

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
  out->uvdata = ObitUVUnref(out->uvdata);
  out->uvwork = ObitUVUnref(out->uvwork);
  out->mosaic = ObitImageMosaicUnref(out->mosaic);
  out->uvdata = ObitUVRef(in->uvdata);
  out->uvwork = ObitUVRef(in->uvwork);
  out->mosaic = ObitImageMosaicRef(in->mosaic);

} /* end ObitUVImagerClone */

/**
 * Creates an ObitUVImager given an ObitUV with control information.
 * The output ImageMosaic member is created
 * \param name   An optional name for the object.
 * \param uvdata ObitUV object with info member containng the output image
 *               specifications and all processing parameters.
 * \param err Obit error stack object.
 * \return the new object.
 */
ObitUVImager* ObitUVImagerCreate (gchar* name, ObitUV *uvdata, ObitErr *err)
{
  ObitUVImager* out=NULL;
  gchar *routine = "ObitUVImagerCreate";

  /* Error checks */
  if (err->error) return out;
  g_assert(ObitUVIsA(uvdata));

  /* Create basic structure */
  out = newObitUVImager (name);

  /* Save uvdata */
  out->uvdata = ObitUVRef(uvdata);

  /* Create output mosaic */
  out->mosaic = ObitImageMosaicCreate (name, uvdata, err);
  if (err->error) Obit_traceback_val (err, routine, name, out);

  /* Define images */
  ObitImageMosaicDefine (out->mosaic, uvdata, TRUE, err);
  if (err->error) Obit_traceback_val (err, routine, name, out);

  return out;
} /* end ObitUVImagerCreate */

/**
 * Creates an ObitUVImager given an ObitUV with control information
 * and a previously existing ImageMosaic
 * The output ImageMosaic member is created
 * \param name   An optional name for the object.
 * \param uvdata ObitUV object with info member containng the output image
 *               specifications and all processing parameters.
 * \param mosaic ImageMosaic to use
 * \param err Obit error stack object.
 * \return the new object.
 */
ObitUVImager* ObitUVImagerCreate2 (gchar* name, ObitUV *uvdata, 
				   ObitImageMosaic *mosaic, ObitErr *err)
{
  ObitUVImager* out=NULL;

  /* Error checks */
  if (err->error) return out;
  g_assert(ObitUVIsA(uvdata));

  /* Create basic structure */
  out = newObitUVImager (name);

  /* Save uvdata */
  out->uvdata = ObitUVRef(uvdata);

  /* Save mosaic */
  out->mosaic = ObitImageMosaicRef(mosaic);

  return out;
} /* end ObitUVImagerCreate2 */

/**
 * Apply weighting to uvdata and write to uvwork member
 * \param in  The input object
 * \param err Obit error stack object.
 */
void ObitUVImagerWeight (ObitUVImager *in, ObitErr *err)
{
  /* List of control parameters on uvwork */
  gchar *controlList[] = 
    {"do3D", "FOV", "doFull", "NField", "xCells", "yCells", "nx", "ny", 
     "RAShift", "DecShift", "Sources",  "Beam",
     "Catalog", "CatDisk", "OutlierDist", "OutlierFlux", "OutlierSI", "OutlierSize",
     "nuGrid", "nvGrid", "WtBox", "WtFunc", "UVTaper", "Robust", "WtPower",
     "RobustIF", "TaperIF", "MFTaper", "doGPU",
     NULL};
  gchar *routine = "ObitUVImagerWeight";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

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

} /* end ObitUVImagerWeight */

/**
 * Image data in uvwork if defined, else uvdata writing results in mosaic.
 * If an autoCenter image is imaged, its shifted version is also shifted.
 * \param in        The input object
 * \param field     zero terminated list of field numbers to image, 0=> all
 * \param doWeight  If TRUE do Weighting ov uv data first
 *                  If TRUE then input data is modified.
 * \param doBeam    If True calculate dirty beams first
 * \param doFlatten If TRUE, flatten images when done
 * \param err       Obit error stack object.
 */
void ObitUVImagerImage (ObitUVImager *in, olong *field, gboolean doWeight, 
			gboolean doBeam, gboolean doFlatten, ObitErr *err)
{ 
  ObitUV *data=NULL;
  ObitImage **imageList=NULL;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitUVDesc *UVDesc;
  ObitImageDesc *imageDesc;
  olong i, n, fldCnt, ifield, channel=0, nDo, nLeft, nImage, prtLv, *fldNo=NULL;
  olong NumPar;
  ofloat sumwts[2];
  ObitImage *theBeam=NULL;
  gboolean *forceBeam=NULL, needBeam, doall;
  ObitUVImagerClassInfo *imgClass = (ObitUVImagerClassInfo*)in->ClassInfo;
  gchar        *dataParms[] = {  /* Imaging info */
    "xShift", "yShift",
    NULL
  };
  gchar *routine = "ObitUVImagerImage";

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

  /* Copy imaging info if uvwork already defined */
  if (in->uvwork) 
    ObitInfoListCopyList (in->uvdata->info, in->uvwork->info, dataParms);

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
  if (ObitInfoListGetTest(in->mosaic->info, "prtLv", &type, dim, &prtLv)) 
    err->prtLv = prtLv;  /* Add to err */

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

  NumPar = imgClass->ObitUVImagerGetNumPar(in, doBeam, err); /* How many to do? */

  /* Get list of images */
  imageList = g_malloc0(in->mosaic->numberImages*sizeof(ObitImage*));
  if (doall) n = in->mosaic->numberImages;
  else n = nImage;
  fldCnt = 0;
  for (i=0; i<n; i++) {
    if (doall) ifield = i;
    else ifield = field[i]-1;
    if (in->mosaic->isShift[ifield]<=0) {  /* Special handling for shifted imaged */
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

  /* Make any 2D shifted images */
 shifty:
  imgClass->ObitUVImagerShifty (in, field, doall, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Need to flatten? */
  if (doFlatten) ObitUVImagerFlatten (in, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Cleanup */
  if (forceBeam) g_free(forceBeam);
  if (fldNo) g_free(fldNo);
  if (imageList) g_free(imageList);
} /* end ObitUVImagerImage */

/**
 * Create shifted 2D imaging facets
 * \param in     The input UVImager object
 * \param field  Zero terminated list of field numbers to image, 0=> all
 * \param doall  If true then all images were remade
 * \param err    Obit error stack object.
 */
void ObitUVImagerShifty (ObitUVImager *in, olong *field, gboolean doall, 
			 ObitErr *err)
{
  olong i, j, ofield, ifield;
  ofloat shift[2];
  gboolean found;
  gchar *routine = "ObitUVImagerShifty";

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
	ObitImageUtilShift (in->mosaic->images[ifield], in->mosaic->images[ofield], shift, err);
	/* Beam - no shift */
	shift[0] = shift[1] = 0.0;
	ObitImageUtilShift ((ObitImage*)in->mosaic->images[ifield]->myBeam, 
			    (ObitImage*)in->mosaic->images[ofield]->myBeam, shift, err);
      }
    } /* end if shifted field */
  } /* End loop making shifted images */
  if (err->error) Obit_traceback_msg (err, routine, in->name);

} /* end ObitUVImageShifty */

/**
 * Flatten Image Mosaic
 * \param in  The input object
 * \param err Obit error stack object.
 */
void ObitUVImagerFlatten (ObitUVImager *in, ObitErr *err)
{
  gchar *routine = "ObitUVImagerFlatten";
  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  if (!ObitImageMosaicIsA(in->mosaic)) {
    Obit_log_error(err, OBIT_Error,"%s ImageMosaic not defined in %s", 
		   routine, in->name);
    return;
  }


  ObitImageMosaicFlatten (in->mosaic, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
} /* end ObitUVImagerFlatten */

/**
 * return ImageMosaic member
 * \param in  The input object
 * \param err Obit error stack object.
 * \return reference to ImageMosaic.
 */
ObitImageMosaic* ObitUVImagerGetMosaic (ObitUVImager *in, ObitErr *err)
{ 
  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return NULL;
  g_assert (ObitIsA(in, &myClassInfo));

  return ObitImageMosaicRef(in->mosaic);
} /* end ObitUVImagerGetMosaic */

/**
 * Convert structure information to entries in an ObitInfoList
 * \param in      Object of interest.
 * \param prefix  If NonNull, string to be added to beginning of outList entry name
 *                "xxx" in the following
 * \param outList InfoList to write entries into
 *      \li "xxxClassType" string UVImager type, "Base" for base class
 *      \li "xxxUVData" prefix for uvdata member, entry with value "None" => doesn't exist
 *      \li "xxxUVWork" prefix for uvwork member, entry with value "None" => doesn't exist
 *      \li "xxxMosaic" prefix for mosaic member, entry with value "None" => doesn't exist
 *      \li various weighting parameters
 * \param err     ObitErr for reporting errors.
 */
void ObitUVImagerGetInfo (ObitUVImager *in, gchar *prefix, ObitInfoList *outList, 
			  ObitErr *err)
{ 
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *keyword=NULL, *None = "None", *OK="OK", *Type="Base";
  olong i;
  gpointer listPnt;
  gchar *parm[] = 
    {"do3D", "FOV", "doFull", "NField", "xCells", "yCells", "nx", "ny", 
     "RAShift", "DecShift", "Sources", 
     "Catalog", "CatDisk", "OutlierDist", "OutlierFlux", "OutlierSI", "OutlierSize",
     "nuGrid", "nvGrid", "WtBox", "WtFunc", "UVTaper", "Robust", "WtPower",
     "RobustIF", "TaperIF", "MFTaper",
     NULL};
  gchar *routine = "ObitUVImagerGetInfo";

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

} /* end ObitUVImagerGetInfo */

/**
 * Get number of parallel images
 * Target memory usage is 0.75 GByte if 32 bit, 3 GByte if 64.
 * \param in      Object of interest.
 * \param doBeam    If True calculate dirty beams first
 * \param err       Obit error stack object.
 * \return the number of parallel images.
 */
olong ObitUVImagerGetNumPar (ObitUVImager *in, gboolean doBeam, ObitErr *err)
{
  olong out=8;
  odouble lenVis, numVis, imSize, bufSize, tSize;
  ObitImage *beam;
  gchar *routine="ObitUVImagerGetNumPar";

  if (err->error) return out;
  g_assert(ObitUVImagerIsA(in));

  /* How big are things? */
  numVis = (odouble)ObitImageUtilBufSize (in->uvdata);  /* Size of buffer */
  lenVis = (odouble)in->uvdata->myDesc->lrec;
  imSize = in->mosaic->images[0]->myDesc->inaxes[0] * 
    in->mosaic->images[0]->myDesc->inaxes[1];  /* Image plane size */
  /* Beam size */
  if (doBeam) {
    beam = (ObitImage*)in->mosaic->images[0]->myBeam;
    if (beam!=NULL) imSize += beam->myDesc->inaxes[0] * beam->myDesc->inaxes[1];
    else imSize *= 5;  /* Assume 4X as large (plus map) */
  }

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
} /*  end ObitUVImagerGetNumPar */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitUVImagerClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitUVImagerClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitUVImagerClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitUVImagerClassInfoDefFn (gpointer inClass)
{
  ObitUVImagerClassInfo *theClass = (ObitUVImagerClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitUVImagerClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitUVImagerClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitUVImagerGetClass;
  theClass->newObit       = (newObitFP)newObitUVImager;
  theClass->ObitUVImagerFromInfo= (ObitUVImagerFromInfoFP)ObitUVImagerFromInfo;
  theClass->ObitCopy      = (ObitCopyFP)ObitUVImagerCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitUVImagerClear;
  theClass->ObitInit      = (ObitInitFP)ObitUVImagerInit;
  theClass->ObitUVImagerCreate = (ObitUVImagerCreateFP)ObitUVImagerCreate;
  theClass->ObitUVImagerCreate2= (ObitUVImagerCreate2FP)ObitUVImagerCreate2;
  theClass->ObitUVImagerWeight = (ObitUVImagerWeightFP)ObitUVImagerWeight;
  theClass->ObitUVImagerImage  = (ObitUVImagerImageFP)ObitUVImagerImage;
  theClass->ObitUVImagerShifty = (ObitUVImagerShiftyFP)ObitUVImagerShifty;
  theClass->ObitUVImagerFlatten= (ObitUVImagerFlattenFP)ObitUVImagerFlatten;
  theClass->ObitUVImagerGetMosaic = 
    (ObitUVImagerGetMosaicFP)ObitUVImagerGetMosaic;
  theClass->ObitUVImagerGetInfo= (ObitUVImagerGetInfoFP)ObitUVImagerGetInfo;
  theClass->ObitUVImagerGetNumPar= (ObitUVImagerGetNumParFP)ObitUVImagerGetNumPar;

} /* end ObitUVImagerClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitUVImagerInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitUVImager *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->uvdata = NULL;
  in->uvwork = NULL;
  in->mosaic = NULL;

} /* end ObitUVImagerInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitUVImager* cast to an Obit*.
 */
void ObitUVImagerClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitUVImager *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  in->uvdata = ObitUVUnref(in->uvdata);
  in->uvwork = ObitUVUnref(in->uvwork);
  in->mosaic = ObitImageMosaicUnref(in->mosaic);
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitUVImagerClear */

