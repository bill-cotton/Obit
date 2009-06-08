/* $Id$        */
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
    {"FOV", "doFull", "NField", "xCells", "yCells", "nx", "ny", 
     "RAShift", "DecShift", "Sources", 
     "Catalog",  "OutlierDist", "OutlierFlux", "OutlierSI", "OutlierSize",
     "nuGrid", "nvGrid", "WtBox", "WtFunc", "UVTaper", "Robust", "WtPower",
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
    {"FOV", "doFull", "NField", "xCells", "yCells", "nx", "ny", 
     "RAShift", "DecShift", "Sources", 
     "Catalog",  "OutlierDist", "OutlierFlux", "OutlierSI", "OutlierSize",
     "nuGrid", "nvGrid", "WtBox", "WtFunc", "UVTaper", "Robust", "WtPower",
     NULL};
  gchar *routine = "ObitUVImagerWeight";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

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
 * \param in        The input object
 * \param field     Which field (1-rel) to Image, 0=> all
 * \param doWeight  If TRUE do Weighting ov uv data first
 *                  If TRUE then input data is modified.
 * \param doBeam    If True calculate dirst beams first
 * \param doFlatten If TRUE, flatten images when done
 * \param err       Obit error stack object.
 */
void ObitUVImagerImage (ObitUVImager *in,  olong field, gboolean doWeight, 
			gboolean doBeam, gboolean doFlatten, ObitErr *err)
{ 
  ObitUV *data=NULL;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitUVDesc *UVDesc;
  ObitImageDesc *imageDesc;
  olong ifield, hiField, loField, channel=0, nDo, nLeft;
  ofloat sumwts;
  ObitImage *theBeam=NULL;
  gboolean *forceBeam=NULL, needBeam;
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

  /* List of need to force making beam */
  forceBeam = g_malloc0(in->mosaic->numberImages*sizeof(gboolean));

 /* Single or multiple images (including beams) */
  if ((!doBeam) && ((field>0) || (in->mosaic->numberImages==1))) {
    /* single - no beam */
    UVDesc    = data->myDesc;
    ifield = MAX (0, (field-1));
    imageDesc = in->mosaic->images[ifield]->myDesc;
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

    return;
  } /* end single */

  /* Multiple (including beams) - do in parallel */
  loField = 0;
  hiField = in->mosaic->numberImages-1;

  /* Loop over fields initializing */
  for (ifield=loField; ifield<=hiField; ifield++) {
    /* Set Stokes */
    UVDesc    = data->myDesc;
    imageDesc = in->mosaic->images[ifield]->myDesc;
    imageDesc->crval[imageDesc->jlocs] = UVDesc->crval[UVDesc->jlocs];

    /* reset image max/min */
    imageDesc->maxval    = -1.0e20;
    imageDesc->minval    =  1.0e20;
 
    /* Need to force beam? */
    theBeam = (ObitImage*)in->mosaic->images[ifield]->myBeam;
    forceBeam[ifield-loField] = (theBeam==NULL) ||
      !ObitInfoListGetTest(theBeam->info, "SUMWTS", &type, dim, (gpointer)&sumwts);

  } /* end loop initializing */

  /* Image - don't do more than 8 at a time */
  nLeft = in->mosaic->numberImages;
  ifield = 0;
  while (nLeft>0) {
    nDo = MIN (8, nLeft);
    ObitImageUtilMakeImagePar (data, nDo, &in->mosaic->images[ifield],
			       doBeam, doWeight, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    ifield += nDo;
    nLeft  -= nDo;
  } /* End loop gridding */

  /* Loop over fields finalizing */
  for (ifield=loField; ifield<=hiField; ifield++) {
    /* If it made a beam check the beam size */
    if (doBeam || forceBeam[ifield-loField]) {
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
  } /* End loop over fields finalizing*/
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Need to flatten? */
  if (doFlatten) ObitUVImagerFlatten (in, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  /* Cleanup */
  if (forceBeam) g_free(forceBeam);
} /* end ObitUVImagerImage */

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
    {"FOV", "doFull", "NField", "xCells", "yCells", "nx", "ny", 
     "RAShift", "DecShift", "Sources", 
     "Catalog",  "OutlierDist", "OutlierFlux", "OutlierSI", "OutlierSize",
     "nuGrid", "nvGrid", "WtBox", "WtFunc", "UVTaper", "Robust", "WtPower",
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
  theClass->ObitUVImagerFlatten= (ObitUVImagerFlattenFP)ObitUVImagerFlatten;
  theClass->ObitUVImagerGetMosaic = 
    (ObitUVImagerGetMosaicFP)ObitUVImagerGetMosaic;
  theClass->ObitUVImagerGetInfo= (ObitUVImagerGetInfoFP)ObitUVImagerGetInfo;

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

