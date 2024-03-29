/* $Id$   */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2006-2023                                          */
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

#include "ObitUVImagerIon.h"
#include "ObitImageUtil.h"
#include "ObitTableSN.h"
#include "ObitSkyGeom.h"
#include "ObitIoN2SolNTable.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitUVImagerIon.c
 * ObitUVImagerIon class function definitions.
 * This class is derived from the Obit base class.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitUVImagerIon";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitUVImagerGetClass;

/**
 * ClassInfo structure ObitUVImagerIonClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitUVImagerIonClassInfo myClassInfo = {FALSE};

/*--------------- File Global Variables  ----------------*/


/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitUVImagerIonInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitUVImagerIonClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitUVImagerIonClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitUVImagerIon* newObitUVImagerIon (gchar* name)
{
  ObitUVImagerIon* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitUVImagerIonClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitUVImagerIon));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitUVImagerIonInit((gpointer)out);

 return out;
} /* end newObitUVImagerIon */

/**
 * Initializes from ObitInfoList.
 * Initializes class if needed on first call.
 * \param out     the new object.to be initialized
 * \param prefix  If NonNull, string to be added to beginning of inList entry name
 *                "xxx" in the following
 * \param inList  InfoList to extract object information from
 *      \li "xxxClassType" string UVImager type, "Ion" for this class
 *      \li "xxxUVData" prefix for uvdata member, entry with value "None" => doesn't exist
 *      \li "xxxUVWork" prefix for uvwork member, entry with value "None" => doesn't exist
 *      \li "xxxMosaic" prefix for mosaic member, entry with value "None" => doesn't exist
 *      \li "xxxIonTab" prefix for ionTab member, entry with value "None" => doesn't exist
 *      \li "xxxsnVer"  olong SN table version to use, 0> add new
 * \param err     ObitErr for reporting errors.
 * \return the new object.
 */
void ObitUVImagerIonFromInfo (ObitUVImager *out, gchar *prefix, 
			      ObitInfoList *inList, ObitErr *err)
{ 
  ObitUVImagerIon *ionout = NULL;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *keyword=NULL, *None = "None", *value=NULL, *Type = "Ion";;
  gboolean missing;
  gchar *routine = "ObitUVImagerIonFromInfo";

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitUVImagerClassInit();

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(out, &myClassInfo));
 
  /* check class type */
  if (prefix) keyword = g_strconcat (prefix, "ClassType", NULL);
  else        keyword = g_strdup("ClassType");
  missing = ObitInfoListGetP(inList, keyword, &type, dim, (gpointer*)&value);
  if ((missing) || (type!=OBIT_string) || (!strncmp(Type,value,dim[0]))) {
    Obit_log_error(err, OBIT_Error,"%s Wrong class type %s!=%s", routine, value, Type);
    return;
  }
  g_free(keyword);

  /* pointer to this class type */
  ionout = (ObitUVImagerIon*)out;

  /* Add things for this class */
  /* ionTab */
  if (prefix) keyword = g_strconcat (prefix, "IonTab", NULL);
  else        keyword = g_strdup("IonTab");
  missing = ObitInfoListGetP(inList, keyword, &type, dim, (gpointer*)&value);
  /* Does it exist? */
  if ((missing) || (type!=OBIT_string) || (!strncmp(None,value,dim[0]))) {
    Obit_log_error(err, OBIT_Error,"%s Ion model table not defined in %s", routine, keyword);
    return;
  } else { /* exists*/
    ionout->ionTab = (ObitTableNI*)ObitTableFromFileInfo(keyword, inList, err);
    if (err->error) Obit_traceback_msg (err, routine, keyword);
  }
  g_free(keyword);

 /* snVer */
  if (prefix) keyword = g_strconcat (prefix, "snVer", NULL);
  else        keyword = g_strdup("snVer");
  dim[0] = 1;
  ObitInfoListAlwaysPut(inList, keyword, OBIT_long, dim, &ionout->snVer);
  g_free(keyword);

} /* end ObitUVImagerIonFromInfo */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitUVImagerIonGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitUVImagerIonClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitUVImagerIonGetClass */

/**
 * Make a deep copy of an ObitUVImagerIon.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitUVImagerIon* ObitUVImagerIonCopy  (ObitUVImagerIon *in, ObitUVImagerIon *out, ObitErr *err)
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
    out = newObitUVImagerIon(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class - just pointers */
  out->ionTab = ObitTableNIUnref(out->ionTab);
  out->ionTab = ObitTableNIRef(in->ionTab);
  out->snVer  = in->snVer;

  return out;
} /* end ObitUVImagerIonCopy */

/**
 * Make a copy of a object but do not copy the actual data
 * This is useful to create an UVImagerIon similar to the input one.
 * \param in  The object to copy
 * \param out An existing object pointer for output, must be defined.
 * \param err Obit error stack object.
 */
void ObitUVImagerIonClone  (ObitUVImagerIon *in, ObitUVImagerIon *out, ObitErr *err)
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
  out->ionTab = ObitTableNIUnref(out->ionTab);
  out->ionTab = ObitTableNIRef(in->ionTab);
  out->snVer  = in->snVer;

} /* end ObitUVImagerIonClone */

/**
 * Creates an ObitUVImagerIon given an ObitUV with control information.
 * The output ImageMosaic member is created
 * \param name   An optional name for the object.
 * \param uvdata ObitUV object with info member containng the output image
 *               specifications and all processing parameters.
 * \param err Obit error stack object.
 * \return the new object.
 */
ObitUVImagerIon* ObitUVImagerIonCreate (gchar* name, ObitUV *uvdata, ObitErr *err)
{
  ObitUVImagerIon* out=NULL;
  olong NIver, highVer;
  olong itemp;
  oint ncoef=0;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *tabType = "NI", *snType = "SN";
  gchar *routine = "ObitUVImagerIonCreate";

  /* Error checks */
  if (err->error) return out;
  g_assert(ObitUVIsA(uvdata));

  /* Create basic structure */
  out = newObitUVImagerIon (name);

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

  /* Get NI table */
  itemp = 1;
  ObitInfoListGetTest(uvdata->info, "ionVer", &type, dim, &itemp);
  NIver = itemp;
  /* 0=> highest */
  if (NIver<=0) NIver = ObitTableListGetHigh (uvdata->tableList, tabType);
  out->ionTab = newObitTableNIValue ("Ion cal", (ObitData*)uvdata, &NIver, 
				    OBIT_IO_ReadOnly, ncoef, err);
  if (err->error) Obit_traceback_val (err, routine, uvdata->name, out);
  Obit_retval_if_fail ((out->ionTab!=NULL), err, out,
		       "%s: NI table %d not found", routine, NIver);
		       

  /* Current highest SN table number */
  highVer =  ObitTableListGetHigh (uvdata->tableList, snType);
  out->snVer = highVer+1;

  return out;
} /* end ObitUVImagerIonCreate */

/**
 * Creates an ObitUVImagerIon given an ObitUV with control information
 * and a previously existing ImageMosaic
 * The output ImageMosaic member is created
 * \param name   An optional name for the object.
 * \param uvdata ObitUV object with info member containng the output image
 *               specifications and all processing parameters.
 * \param mosaic ImageMosaic to use
 * \param err Obit error stack object.
 * \return the new object.
 */
ObitUVImagerIon* ObitUVImagerIonCreate2 (gchar* name, ObitUV *uvdata, 
					 ObitImageMosaic *mosaic, ObitErr *err)
{
  ObitUVImagerIon* out=NULL;

  /* Error checks */
  if (err->error) return out;
  g_assert(ObitUVIsA(uvdata));

  /* Create basic structure */
  out = newObitUVImagerIon (name);

  /* Save uvdata */
  out->uvdata = ObitUVRef(uvdata);

  /* Save mosaic */
  out->mosaic = ObitImageMosaicRef(mosaic);

  return out;
} /* end ObitUVImagerIonCreate2 */

/**
 * Apply weighting to uvdata and write to uvwork member
 * \param in  The input object
 * \param err Obit error stack object.
 */
void ObitUVImagerIonWeight (ObitUVImager *in, ObitErr *err)
{
  const ObitUVImagerClassInfo *ParentClass;
  gchar *NIlist[] = {"AIPS NI", NULL};
  gchar *routine = "ObitUVImagerIonWeight";

  /* Call parent class routine */
  ParentClass = ((ObitUVImagerClassInfo*)(in->ClassInfo))->ParentClass;
  ParentClass->ObitUVImagerWeight(in, err);

  /* Copy NI Table */
  ObitUVCopyTables (in->uvdata, in->uvwork, NULL, NIlist, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
} /* end ObitUVImagerIonWeight */

/**
 * Image data in uvwork if defined, else uvdata writing results in mosaic.
 * \param in        The input object
 * \param field     zero terminated list of field numbers to image, 0=> all
 * \param doWeight  If TRUE do Weighting ov uv data first
 *                  If TRUE then input data is modified.
 * \param doBeam    If True calculate dirst beams first
 * \param doFlatten If TRUE, flatten images when done
 * \param err       Obit error stack object.
 */
void ObitUVImagerIonImage (ObitUVImager *inn,  olong *field, gboolean doWeight, 
			   gboolean doBeam, gboolean doFlatten, ObitErr *err)
{ 
  ObitUVImagerIon *in=(ObitUVImagerIon*)inn;
  ObitUV *data=NULL;
  ObitImage **imageList=NULL;
  ObitTableSN *snTab=NULL;
  ObitUVDesc *UVDesc;
  ObitImageDesc *imageDesc;
  ObitInfoType type;
  ObitImage *theBeam=NULL;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong itemp, prtLv;
  gboolean Tr=TRUE, Fl=FALSE, needBeam, *forceBeam=NULL, doall, found=FALSE;
  ofloat shift[2], sumwts;
  olong i, j, n, fldCnt, nImage, ifield, ofield=1, channel=0, lver;
  olong *fldNo=NULL, nLeft, nDo, NumPar;
  gchar *tabtype=NULL;
  gchar *routine = "ObitUVImagerIonImage";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Count number of fields to image */
  nImage = 0;
  for (i=0; i<in->mosaic->numberImages; i++) {
    if (field[i]==0) break;
    nImage++;
  }
  /* All wanted? */
  doall =  (nImage <= 0);

  NumPar = ObitUVImagerIonGetNumPar(inn, doBeam, err); /* How many to do? */

  /* List of need to force making beam, keep track of field numbers */
  forceBeam = g_malloc0(in->mosaic->numberImages*sizeof(gboolean));
  fldNo     = g_malloc0(in->mosaic->numberImages*sizeof(olong));
  
  /* Get list of images */
  imageList = g_malloc0(in->mosaic->numberImages*sizeof(ObitImage*));
  if (doall) n = in->mosaic->numberImages;
  else n = nImage;
  fldCnt = 0;
  for (i=0; i<n; i++) {
    if (doall) ifield = i+1;
    else ifield = field[i];
    if (in->mosaic->isShift[ifield-1]) {  /* Special handling for shifted imaged */
      imageList[fldCnt++] = in->mosaic->images[ifield-1];
      fldNo[i] = ifield-1;                 /* Number (0-rel) in mosaic */
      /* Need to force beam? */
      theBeam = (ObitImage*)imageList[i]->myBeam;
      forceBeam[i] = (theBeam==NULL) ||
	!ObitInfoListGetTest(theBeam->info, "SUMWTS", &type, dim, &sumwts);
    }
  }
  
  /* Weighting? */
  if (doWeight) {
    ObitUVImagerIonWeight (inn, err);
    /* Copy NI Table
    ObitUVCopyTables (in->uvdata, in->uvwork, NULL, NIlist, err); */
    if (err->error) Obit_traceback_msg (err, routine, in->name);
  }

  data = in->uvwork;
  if (!data) data = in->uvdata;
  if (!ObitUVIsA(data)) {
    Obit_log_error(err, OBIT_Error,"%s UV data not defined in %s", routine, data->name);
    return;
  }

  /* Get print level */
  prtLv = 0;
  ObitInfoListGetTest(in->uvdata->info, "prtLv", &type, dim, &prtLv);

  /* Image loop */
  nLeft  = fldCnt;
  ifield = 0;
  while (nLeft>0) {
    nDo = MIN (NumPar, nLeft);
    
    /* Set Stokes */
    UVDesc    = data->myDesc;
    imageDesc = imageList[ifield]->myDesc;
    imageDesc->crval[imageDesc->jlocs] = UVDesc->crval[UVDesc->jlocs];

    /* reset image max/min */
    imageDesc->maxval    = -1.0e20;
    imageDesc->minval    =  1.0e20;

    /* Position shift */
    shift[0] = in->mosaic->RAShift[fldNo[ifield]];
    shift[1] = in->mosaic->DecShift[fldNo[ifield]];

    /* Create SN table with values for this field from input NI table */
    snTab = ObitIoN2SolNTableConvert(data, in->ionTab, snTab, shift, err);

    /* Setup to use this SN table */
    dim[0] = 1;
    itemp = snTab->tabVer;
    ObitInfoListAlwaysPut(data->info, "gainUse",     OBIT_long, dim, &itemp);
    ObitInfoListAlwaysPut(data->info, "doCalSelect", OBIT_bool, dim, &Tr);
    itemp = 2;
    ObitInfoListAlwaysPut(data->info, "doCalib",     OBIT_long, dim, &itemp);
 
    /* diagnostics? */
    if (prtLv>=3) {
	Obit_log_error(err, OBIT_InfoErr, 
		       "%s: Field %d Using NI table %d with shift %lf %lf",
		       routine, fldNo[ifield]+1, in->ionTab->tabVer, shift[0], shift[1]);
    }


    /* Need to force beam? */
    needBeam = doBeam || forceBeam[ifield];

    /* Image */
    ObitImageUtilMakeImage (data, imageList[ifield], channel, needBeam, Fl, err);
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
  
  /* Loop making any (2D) shifted images */
  for (i=0; i<in->mosaic->numberImages; i++) {
    /* Only needed for 2D */
    if (in->mosaic->images[i]->myDesc->do3D) continue;
    ofield = i;
    ifield = in->mosaic->isShift[i]-1;
    if (ifield >= 0) {
      /* See if input image just made - named in field */
      found = doall;
      for (j=0; j<in->mosaic->numberImages; j++) {
	if (field[j]==0) break;
	if (field[j]==(ifield+1)) found = TRUE;
      }
      if (found) {
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
  
  /* Turn off calibration */
  dim[0] = 1;
  itemp = 0;
  ObitInfoListAlwaysPut(data->info, "gainUse",     OBIT_long, dim, &itemp);
  ObitInfoListAlwaysPut(data->info, "doCalSelect", OBIT_bool, dim, &Fl);
  itemp = -1;
  ObitInfoListAlwaysPut(data->info, "doCalib",     OBIT_long, dim, &itemp);
  
  if (err->error) goto cleanup;
  
  /* Diagnostics or clear messages */
  if (prtLv>=1)  ObitErrLog(err);
  else ObitErrClear(err);
  
  /* Need to flatten? */
  if (doFlatten) ObitUVImagerFlatten (inn, err);
  if (err->error) goto cleanup;
  
  /* Cleanup*/
 cleanup:
  if (fldNo)     g_free(fldNo);
  if (forceBeam) g_free(forceBeam);
  if (imageList) g_free(imageList);
  if (snTab) {
    tabtype = g_strdup(snTab->tabType); /* Save ephemerial values */
    lver = snTab->tabVer;
    ObitUVZapTable (data, tabtype, lver, err);
    if (tabtype) {g_free (tabtype); tabtype = NULL;}
    ObitUVUpdateTables (data, err);  /* Update disk header */
    snTab = ObitTableSNUnref(snTab);
  }
  if (err->error) Obit_traceback_msg (err, routine, in->name);
} /* end ObitUVImagerIonImage */

/**
 * Flatten Image Mosaic
 * \param in  The input object
 * \param err Obit error stack object.
 */
void ObitUVImagerIonFlatten (ObitUVImager *in, ObitErr *err)
{
  const ObitUVImagerClassInfo *ParentClass;

  /* Call parent class routine */
  ParentClass = ((ObitUVImagerClassInfo*)(in->ClassInfo))->ParentClass;
  ParentClass->ObitUVImagerFlatten(in, err);
} /* end ObitUVImagerIonFlatten */

/**
 * return ImageMosaic member
 * \param in  The input object
 * \param err Obit error stack object.
 * \return reference to ImageMosaic.
 */
ObitImageMosaic* ObitUVImagerIonGetMosaic (ObitUVImager *in, ObitErr *err)
{ 
  const ObitUVImagerClassInfo *ParentClass;

  /* Call parent class routine */
  ParentClass = ((ObitUVImagerClassInfo*)(in->ClassInfo))->ParentClass;
  return  ParentClass->ObitUVImagerGetMosaic(in, err);
} /* end ObitUVImagerIonGetMosaic */

/**
 * Convert structure information to entries in an ObitInfoList
 * \param in      Object of interest.
 * \param prefix  If NonNull, string to be added to beginning of outList entry name
 *                "xxx" in the following
 * \param outList InfoList to write entries into
 *      \li "xxxClassType" string UVImager type, "Ion" for this class
 *      \li "xxxUVData" prefix for uvdata member, entry with value "None" => doesn't exist
 *      \li "xxxUVWork" prefix for uvwork member, entry with value "None" => doesn't exist
 *      \li "xxxMosaic" prefix for mosaic member, entry with value "None" => doesn't exist
 *      \li "xxxIonTab" prefix for ionTab member, entry with value "None" => doesn't exist
 *      \li "xxxsnVer"  olong SN table version to use, 0> add new
 * \param err     ObitErr for reporting errors.
 */
void ObitUVImagerIonGetInfo (ObitUVImager *inn, gchar *prefix, ObitInfoList *outList, 
			     ObitErr *err)
{ 
  ObitUVImagerIon *in = (ObitUVImagerIon*)inn;
   gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *keyword=NULL, *None = "None", *OK="OK", *Type="Ion";
  gchar *routine = "ObitUVImagerIonGetInfo";

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Use Base class */
  ObitUVImagerGetInfo(inn, prefix, outList, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Add things for this class */
  /* ionTab */
  if (prefix) keyword = g_strconcat (prefix, "IonTab", NULL);
  else        keyword = g_strdup("IonTab");
  if (in->ionTab) {
    ObitTableGetFileInfo((ObitTable*)in->ionTab, keyword, outList, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    dim[0] = strlen(OK);
    ObitInfoListAlwaysPut(outList, keyword, OBIT_string, dim, OK);
  } else {
    dim[0] = strlen(None);
    ObitInfoListAlwaysPut(outList, keyword, OBIT_string, dim, None);
  }
  g_free(keyword);

  /* snVer */
  if (prefix) keyword = g_strconcat (prefix, "snVer", NULL);
  else        keyword = g_strdup("snVer");
  dim[0] = 1;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_long, dim, &in->snVer);
  g_free(keyword);

  /* set Class type */
  if (prefix) keyword = g_strconcat (prefix, "ClassType", NULL);
  else        keyword = g_strdup("ClassType");
  dim[0] = strlen(Type);
  ObitInfoListAlwaysPut(outList, keyword, OBIT_string, dim, Type);

} /* end ObitUVImagerIonGetInfo */

/**
 * Get number of parallel images
 * Can only do 1
 * \param in      Object of interest.
 * \param doBeam    If True calculate dirty beams first
 * \param err       Obit error stack object.
 * \return the number of parallel images.
 */
olong ObitUVImagerIonGetNumPar (ObitUVImager *inn, gboolean doBeam, ObitErr *err)
{
  /*ObitUVImagerIon *in = (ObitUVImagerIon*)inn;*/
  olong out=1;

  return out;
} /*  end ObitUVImagerIonGetNumPar */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitUVImagerIonClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitUVImagerIonClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitUVImagerIonClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitUVImagerIonClassInfoDefFn (gpointer inClass)
{
  ObitUVImagerIonClassInfo *theClass = (ObitUVImagerIonClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitUVImagerIonClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitUVImagerIonClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitUVImagerIonGetClass;
  theClass->newObit       = (newObitFP)newObitUVImagerIon;
  theClass->ObitCopy      = (ObitCopyFP)ObitUVImagerIonCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitUVImagerIonClear;
  theClass->ObitInit      = (ObitInitFP)ObitUVImagerIonInit;
  theClass->ObitUVImagerCreate = (ObitUVImagerCreateFP)ObitUVImagerIonCreate;
  theClass->ObitUVImagerCreate2= (ObitUVImagerCreate2FP)ObitUVImagerIonCreate2;
  theClass->ObitUVImagerWeight = (ObitUVImagerWeightFP)ObitUVImagerIonWeight;
  theClass->ObitUVImagerImage  = (ObitUVImagerImageFP)ObitUVImagerIonImage;
  theClass->ObitUVImagerGetInfo= (ObitUVImagerGetInfoFP)ObitUVImagerIonGetInfo;
  theClass->ObitUVImagerGetNumPar= (ObitUVImagerGetNumParFP)ObitUVImagerIonGetNumPar;

} /* end ObitUVImagerIonClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitUVImagerIonInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitUVImagerIon *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->ionTab = NULL;
  in->snVer  = 0;

} /* end ObitUVImagerIonInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitUVImagerIon* cast to an Obit*.
 */
void ObitUVImagerIonClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitUVImagerIon *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  in->ionTab = ObitTableNIUnref(in->ionTab);
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitUVImagerIonClear */

