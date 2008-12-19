/* $Id$     */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2006-2008                                          */
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

#include "ObitUVImagerSquint.h"
#include "ObitImageUtil.h"
#include "ObitUVWeight.h"
#include "ObitSkyGeom.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitUVImagerSquint.c
 * ObitUVImagerSquint class function definitions.
 * This class is derived from the Obit base class.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitUVImagerSquint";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitUVImagerGetClass;

/**
 * ClassInfo structure ObitUVImagerSquintClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitUVImagerSquintClassInfo myClassInfo = {FALSE};

/*--------------- File Global Variables  ----------------*/


/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitUVImagerSquintInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitUVImagerSquintClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitUVImagerSquintClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitUVImagerSquint* newObitUVImagerSquint (gchar* name)
{
  ObitUVImagerSquint* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitUVImagerSquintClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitUVImagerSquint));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitUVImagerSquintInit((gpointer)out);

 return out;
} /* end newObitUVImagerSquint */

/**
 * Initializes from ObitInfoList.
 * Initializes class if needed on first call.
 * \param out     the new object.to be initialized
 * \param prefix  If NonNull, string to be added to beginning of inList entry name
 *                "xxx" in the following
 * \param inList  InfoList to extract object information from 
 *      \li "xxxClassType" string UVImager type, "Squint" for this class
 *      \li "xxxUVData" prefix for uvdata member, entry with value "None" => doesn't exist
 *      \li "xxxUVWork" prefix for uvwork member, entry with value "None" => doesn't exist
 *      \li "xxxMosaic" prefix for mosaic member, entry with value "None" => doesn't exist
 * \param err     ObitErr for reporting errors.
 */
void ObitUVImagerSquintFromInfo (ObitUVImager *out, gchar *prefix, ObitInfoList *inList, 
				 ObitErr *err)
{ 
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *keyword=NULL, *value=NULL;
  gboolean missing;
  gchar *Type = "Squint";
  gchar *routine = "ObitUVImagerSquintFromInfo";
  
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitUVImagerSquintClassInit();

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(out, &myClassInfo));

  /* check class type */
  missing = ObitInfoListGetP(inList, keyword, &type, dim, (gpointer*)&value);
  if ((missing) || (type!=OBIT_string) || (!strncmp(Type,value,dim[0]))) {
    Obit_log_error(err, OBIT_Error,"%s Wrong class type %s!=%s", routine, value, Type);
    return;
  }

} /* end ObitUVImagerSquintFromInfo */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitUVImagerSquintGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitUVImagerSquintClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitUVImagerSquintGetClass */

/**
 * Make a deep copy of an ObitUVImagerSquint.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitUVImagerSquint* ObitUVImagerSquintCopy  (ObitUVImagerSquint *in, 
					     ObitUVImagerSquint *out, ObitErr *err)
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
    out = newObitUVImagerSquint(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class - just pointers */

  return out;
} /* end ObitUVImagerSquintCopy */

/**
 * Make a copy of a object but do not copy the actual data
 * This is useful to create an UVImagerSquint similar to the input one.
 * \param in  The object to copy
 * \param out An existing object pointer for output, must be defined.
 * \param err Obit error stack object.
 */
void ObitUVImagerSquintClone  (ObitUVImagerSquint *in, ObitUVImagerSquint *out, ObitErr *err)
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

} /* end ObitUVImagerSquintClone */

/**
 * Creates an ObitUVImagerSquint given an ObitUV with control information.
 * The output ImageMosaic member is created
 * \param name   An optional name for the object.
 * \param uvdata ObitUV object with info member containng the output image
 *               specifications and all processing parameters.
 * \param err Obit error stack object.
 * \return the new object.
 */
ObitUVImagerSquint* ObitUVImagerSquintCreate (gchar* name, ObitUV *uvdata, ObitErr *err)
{
  ObitUVImagerSquint* out=NULL;
  gchar *routine = "ObitUVImagerSquintCreate";

  /* Error checks */
  if (err->error) return out;
  g_assert(ObitUVIsA(uvdata));

  /* Create basic structure */
  out = newObitUVImagerSquint (name);

  /* Verify uv data */
  ObitDataFullInstantiate ((ObitData*)uvdata, TRUE, err);
  if (err->error) Obit_traceback_val (err, routine, uvdata->name, out);

  /* Save uvdata */
  out->uvdata = ObitUVRef(uvdata);

  /* Create output mosaic */
  out->mosaic = ObitImageMosaicCreate (name, uvdata, err);
  if (err->error) Obit_traceback_val (err, routine, uvdata->name, out);

  /* Define images */
  ObitImageMosaicDefine (out->mosaic, uvdata, TRUE, err);
  if (err->error) Obit_traceback_val (err, routine, uvdata->name, out);

  return out;
} /* end ObitUVImagerSquintCreate */

/**
 * Creates an ObitUVImagerSquint given an ObitUV with control information
 * and a previously existing ImageMosaic
 * The output ImageMosaic member is created
 * \param name   An optional name for the object.
 * \param uvdata ObitUV object with info member containng the output image
 *               specifications and all processing parameters.
 * \param mosaic ImageMosaic to use
 * \param err Obit error stack object.
 * \return the new object.
 */
ObitUVImagerSquint* ObitUVImagerSquintCreate2 (gchar* name, ObitUV *uvdata, 
					       ObitImageMosaic *mosaic, ObitErr *err)
{
  ObitUVImagerSquint* out=NULL;

  /* Error checks */
  if (err->error) return out;
  g_assert(ObitUVIsA(uvdata));

  /* Create basic structure */
  out = newObitUVImagerSquint (name);

  /* Save uvdata */
  out->uvdata = ObitUVRef(uvdata);

  /* Save mosaic */
  out->mosaic = ObitImageMosaicRef(mosaic);

  return out;
} /* end ObitUVImagerSquintCreate2 */

/**
 * Apply weighting to uvdata and write to uvwork member
 * \param in  The input object
 * \param err Obit error stack object.
 */
void ObitUVImagerSquintWeight (ObitUVImager *in, ObitErr *err)
{
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  gboolean Tr=TRUE;
  gchar *Stokes="HALF", IStokes[5];
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
 
 /* Get Stokes being imaged */
  strncpy (IStokes, "F   ", 4); 
  ObitInfoListGetTest (in->uvdata->info, "Stokes", &type, dim, IStokes);

  /* Want RR, LL data */
  dim[0] = 4;
  ObitInfoListAlwaysPut (in->uvdata->info, "Stokes", OBIT_string, dim, Stokes);

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

} /* end ObitUVImagerSquintWeight */

/**
 * Flatten Image Mosaic
 * \param in  The input object
 * \param err Obit error stack object.
 */
void ObitUVImagerSquintFlatten (ObitUVImager *in, ObitErr *err)
{
  const ObitUVImagerClassInfo *ParentClass;

  /* Call parent class routine */
  ParentClass = ((ObitUVImagerClassInfo*)(in->ClassInfo))->ParentClass;
  ParentClass->ObitUVImagerFlatten(in, err);
} /* end ObitUVImagerSquintFlatten */

/**
 * return ImageMosaic member
 * \param in  The input object
 * \param err Obit error stack object.
 * \return reference to ImageMosaic.
 */
ObitImageMosaic* ObitUVImagerSquintGetMosaic (ObitUVImager *in, ObitErr *err)
{ 
  const ObitUVImagerClassInfo *ParentClass;

  /* Call parent class routine */
  ParentClass = ((ObitUVImagerClassInfo*)(in->ClassInfo))->ParentClass;
  return  ParentClass->ObitUVImagerGetMosaic(in, err);
} /* end ObitUVImagerSquintGetMosaic */

/**
 * Convert structure information to entries in an ObitInfoList
 * \param in      Object of interest.
 * \param prefix  If NonNull, string to be added to beginning of outList entry name
 *                "xxx" in the following
 * \param outList InfoList to write entries into
 *      \li "xxxClassType" string UVImager type, "Squint" for this class
 *      \li "xxxUVData" prefix for uvdata member, entry with value "None" => doesn't exist
 *      \li "xxxUVWork" prefix for uvwork member, entry with value "None" => doesn't exist
 *      \li "xxxMosaic" prefix for mosaic member, entry with value "None" => doesn't exist
 * \param err     ObitErr for reporting errors.
 */
void ObitUVImagerSquintGetInfo (ObitUVImager *inn, gchar *prefix, 
				ObitInfoList *outList, ObitErr *err)
{ 
  ObitUVImagerSquint *in = (ObitUVImagerSquint*)inn;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *keyword=NULL, *Type="Squint";
  gchar *routine = "ObitUVImagerSquintGetInfo";

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Use Base class */
  ObitUVImagerGetInfo(inn, prefix, outList, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* set Class type */
  if (prefix) keyword = g_strconcat (prefix, "ClassType", NULL);
  else        keyword = g_strdup("ClassType");
  dim[0] = strlen(Type);
  ObitInfoListAlwaysPut(outList, keyword, OBIT_string, dim, Type);

} /* end ObitUVImagerSquintGetInfo */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitUVImagerSquintClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitUVImagerSquintClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitUVImagerSquintClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitUVImagerSquintClassInfoDefFn (gpointer inClass)
{
  ObitUVImagerSquintClassInfo *theClass = (ObitUVImagerSquintClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitUVImagerSquintClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitUVImagerSquintClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitUVImagerSquintGetClass;
  theClass->newObit       = (newObitFP)newObitUVImagerSquint;
  theClass->ObitCopy      = (ObitCopyFP)ObitUVImagerSquintCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitUVImagerSquintClear;
  theClass->ObitInit      = (ObitInitFP)ObitUVImagerSquintInit;
  theClass->ObitUVImagerCreate = (ObitUVImagerCreateFP)ObitUVImagerSquintCreate;
  theClass->ObitUVImagerCreate2= (ObitUVImagerCreate2FP)ObitUVImagerSquintCreate2;
  theClass->ObitUVImagerWeight = (ObitUVImagerWeightFP)ObitUVImagerSquintWeight;
  theClass->ObitUVImagerGetInfo= (ObitUVImagerGetInfoFP)ObitUVImagerSquintGetInfo;
  /*theClass->ObitUVImagerImage  = (ObitUVImagerImageFP)ObitUVImagerSquintImage;*/

} /* end ObitUVImagerSquintClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitUVImagerSquintInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitUVImagerSquint *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */

} /* end ObitUVImagerSquintInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitUVImagerSquint* cast to an Obit*.
 */
void ObitUVImagerSquintClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitUVImagerSquint *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitUVImagerSquintClear */

