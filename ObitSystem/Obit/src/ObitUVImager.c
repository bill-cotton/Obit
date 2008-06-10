/* $Id: ObitUVImager.c,v 1.9 2007/05/09 21:18:28 bcotton Exp $        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2005-2008                                          */
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
  ObitUVDesc *UVDesc;
  ObitImageDesc *imageDesc;
  olong ifield, hiField, loField, channel=0;
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

  /* Which field numbers (0-re) */
  if (field>0) {
    loField = field-1;
    hiField = loField;
  } else {  /* All */
    loField = 0;
    hiField = in->mosaic->numberImages-1;
  }

  /* DEBUG 
  fprintf (stderr," %g doBeam %d field %d\n",in->mosaic->bmaj, doBeam, field);*/
  /* Loop over fields Imaging */
  for (ifield=loField; ifield<=hiField; ifield++) {
    /* Set Stokes */
    UVDesc    = data->myDesc;
    imageDesc = in->mosaic->images[ifield]->myDesc;
    imageDesc->crval[imageDesc->jlocs] = UVDesc->crval[UVDesc->jlocs];

    /* reset image max/min */
    imageDesc->maxval    = -1.0e20;
    imageDesc->minval    =  1.0e20;

    /* Image */
    ObitImageUtilMakeImage (data, in->mosaic->images[ifield], channel, 
			    doBeam, doWeight, err);
    /* If it made a beam check the beam size */
    if (doBeam) {
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
  } /* End loop over fields */
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Need to flatten? */
  if (doFlatten) ObitUVImagerFlatten (in, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
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

