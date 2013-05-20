/* $Id$  */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2004-2013                                          */
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

#include "ObitImageMosaic.h"
#include "ObitIOImageFITS.h"
#include "ObitIOImageAIPS.h"
#include "ObitSystem.h"
#include "ObitImageUtil.h"
#include "ObitUVUtil.h"
#include "ObitAIPSDir.h"
#include "ObitUV.h"
#include "ObitSkyGeom.h"
#include "ObitTableVZ.h"
#include "ObitTableSUUtil.h"
#include "ObitTableCCUtil.h"
#include "ObitPBUtil.h"
#include "ObitFFT.h"
#include "ObitMem.h"
/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitImageMosaic.c
 * ObitImageMosaic class function definitions.
 * This class is derived from the Obit base class.
 *
 * This class contains an array of associated astronomical images.
 *
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitImageMosaic";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/*--------------- File Global Variables  ----------------*/
/**
 * ClassInfo structure ObitImageMosaicClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitImageMosaicClassInfo myClassInfo = {FALSE};

/*--------------------------- Macroes ---------------------*/
#ifndef MAXFLD         /* Maxum number of fields */
#define MAXFLD 4096
#endif

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitImageMosaicInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitImageMosaicClear (gpointer in);

/** Private: Bubble sort */
 static void 
 Bubble (ofloat *data, olong* indx, olong number, olong direct);

/** Private: Set Class function pointers. */
static void ObitImageMosaicClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name   An optional name for the object.
 * \param number Number of images
 * \return the new object.
 */
ObitImageMosaic* newObitImageMosaic (gchar* name, olong number)
{
  ObitImageMosaic* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitImageMosaicClassInit();

  /* allocate/init structure */
  out = ObitMemAlloc0Name(sizeof(ObitImageMosaic), "ObitImageMosaic");

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  out->numberImages = number;
  ObitImageMosaicInit((gpointer)out);

 return out;
} /* end newObitImageMosaic */

/**
 * Constructor from ObitInfoList.
 * Initializes class if needed on first call.
 * Also works for derived classes.
 * \param prefix  If NonNull, string to be added to beginning of inList entry name
 *                "xxx" in the following
 * \param inList  InfoList to extract object information from 
 *      \li "xxxnumberImages" olong Number of images in Mosaic
 *      \li "xxxnFlyEye"      olong Number of images already initialized
 *      \li "xxxImagennnnn"   string Prefix for each image nnnnn (1-rel) in images
 *      \li "xxxFullField"    string Prefix for Full field image
 *      \li "xxxFOV"          ofloat Field of view as radius (deg)
 *      \li "xxxRadius"       ofloat Radius of the usable region of a given tile (cells)
 *      \li "xxxxCells"       ofloat Cell Spacing in X (deg)
 *      \li "xxxyCells"       ofloat Cell Spacing in Y (deg)
 *      \li "xxxOutlierSize"  olong requested size (CLEAN window) of outliers
 *      \li "xxxnx"           olong* Number of pixels in X for each image
 *      \li "xxxny"           olong* Number of pixels in Y for each image
 *      \li "xxxnplane"       olong* Number of planes for each image
 *      \li "xxxinFlysEye"    gboolean* Is each facet in Fly's Eye? 
 *      \li "xxxFacetNo"      olong* Untapered facet number - associate various taperings
 *      \li "xxxRAShift"      ofloat* RA shift (deg) for each image
 *      \li "xxxDecShift"     ofloat* Dec shift (deg) for each image
 *      \li "xxxfileType"     olong Are image OBIT_IO_AIPS or OBIT_IO_FITS?
 *      \li "xxximName"       string Name of Mosaic images
 *      \li "xxximClass"      string Class of Mosaic images
 *      \li "xxximSeq"        olong Sequence number of Mosaic images
 *      \li "xxximDisk"       olong* Disk numbers of Mosaic images
 *      \li "xxxdoFull"       boolean Is a full field image desired?
 *      \li "xxxbmaj"         ofloat Restoring beam major axis in deg.
 *      \li "xxxbmin"         ofloat Restoring beam minor axis in deg.
 *      \li "xxxbpa"          ofloat Restoring beam PA in deg.
 *      \li "xxxisAuto"       olong Is this an autoCenter field?
 *      \li "xxxisShift"      olong Is this an autoCenter shifted field?   
 *      \li "xxxnumBeamTapes" olong  Number of Beam tapers (elements in BeamTapes) 
 *      \li "xxxBeamTapes"    ofloat*  List of Beam tapers  
 *      \li "xxxBeamTaper"    ofloat*  Beam Taper per image as FWHM in deg      
 * \param err     ObitErr for reporting errors.
 * \return the new object.
 */
ObitImageMosaic* ObitImageMosaicFromInfo (gchar *prefix, ObitInfoList *inList, 
					  ObitErr *err)
{ 
  ObitImageMosaic *out = NULL;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *keyword=NULL, *None = "None", *value=NULL;
  olong i, numberImages, otemp;
  gboolean missing;
  gchar *routine = "ObitImageMosaicFromInfo";

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitImageMosaicClassInit();

  /* error checks */
  if (err->error) return out;

  /* "xxxnumberImages" olong Number of images in Mosaic; */
  if (prefix) keyword = g_strconcat (prefix, "numberImages", NULL);
  else        keyword = g_strdup("numberImages");
  ObitInfoListGetTest(inList, keyword, &type, dim, &numberImages);
  g_free(keyword);

  /* Create output */
  out = newObitImageMosaic (prefix, numberImages);

  /* Copy any InfoList Parameters */
  if (prefix) keyword = g_strconcat (prefix, "Info", NULL);
  else        keyword = g_strdup("Info");
  ObitInfoListCopyWithPrefix (inList, out->info, keyword, TRUE);
  
  /* "xxxnFlyEye"      olong Number of images already initialized */
  if (prefix) keyword = g_strconcat (prefix, "nFlyEye", NULL);
  else        keyword = g_strdup("nFlyEye");
  ObitInfoListGetTest(inList, keyword, &type, dim, &out->nFlyEye);
  g_free(keyword);

  /* "xxxImagennnnn"  string Prefix for each image nnnnn (1-rel) in images */
  for (i=0; i<out->numberImages; i++) {
    /* Prefix for each entry in images */
    if (prefix) sprintf(keyword, "%sImage%5.5d", prefix, i+1);
    else        sprintf(keyword, "Image%5.5d",   i+1);
    missing = ObitInfoListGetP(inList, keyword, &type, dim, (gpointer)&value);
    /* Does it exist? */
    if ((missing) || (type!=OBIT_string) || (!strncmp(None,value,dim[0]))) {
      Obit_log_error(err, OBIT_Error,"%s image %s not defined i", routine, keyword);
      return out;
    } else { /* exists*/
      out->images[i] = (ObitImage*)ObitDataFromFileInfo(keyword, inList, err);
      if (err->error) Obit_traceback_val (err, routine, keyword, out);
    }
    g_free(keyword);
  } /* end loop over images */

  /* "xxxFullField"    string prefix for Full field image */
  if (prefix) keyword = g_strconcat (prefix, "FullField", NULL);
  else        keyword = g_strdup("FullField");
  missing = ObitInfoListGetP(inList, keyword, &type, dim, (gpointer)&value);
  /* Does it exist? */
  if ((missing) || (type!=OBIT_string) || (!strncmp(None,value,dim[0]))) {
    out->FullField = NULL;
  } else { /* exists*/
    out->FullField = (ObitImage*)ObitDataFromFileInfo(keyword, inList, err);
    if (err->error) Obit_traceback_val (err, routine, keyword, out);
  }

  /* "xxxFOV"          ofloat Field of view as radius (deg) */
  if (prefix) keyword = g_strconcat (prefix, "FOV", NULL);
  else        keyword = g_strdup("FOV");
  ObitInfoListGetTest(inList, keyword, &type, dim, &out->FOV);
  g_free(keyword);

  /* "xxxRadius"       ofloat Radius of the usable region of a given tile (cells) */
  if (prefix) keyword = g_strconcat (prefix, "Radius", NULL);
  else        keyword = g_strdup("Radius");
  ObitInfoListGetTest(inList, keyword, &type, dim, &out->Radius);
  g_free(keyword);

  /* "xxxxCells"       ofloat Cell Spacing in X (deg) */
  if (prefix) keyword = g_strconcat (prefix, "xCells", NULL);
  else        keyword = g_strdup("xCells");
  ObitInfoListGetTest(inList, keyword, &type, dim, &out->xCells);
  g_free(keyword);

  /* "xxxyCells"       ofloat Cell Spacing in Y (deg) */
  if (prefix) keyword = g_strconcat (prefix, "yCells", NULL);
  else        keyword = g_strdup("yCells");
  ObitInfoListGetTest(inList, keyword, &type, dim, &out->yCells);
  g_free(keyword);

  /* "xxxOutlierSize"  olong requested size (CLEAN window) of outliers */
  if (prefix) keyword = g_strconcat (prefix, "OutlierSize", NULL);
  else        keyword = g_strdup("OutlierSize");
  ObitInfoListGetTest(inList, keyword, &type, dim, &out->OutlierSize);
  g_free(keyword);

  /* "xxxnx"           olong* Number of pixels in X for each image; */
  if (prefix) keyword = g_strconcat (prefix, "nx", NULL);
  else        keyword = g_strdup("nx");
  ObitInfoListGetTest(inList, keyword, &type, dim, out->nx);
  g_free(keyword);

  /* "xxxny"           olong* Number of pixels in Y for each image */
  if (prefix) keyword = g_strconcat (prefix, "ny", NULL);
  else        keyword = g_strdup("ny");
  ObitInfoListGetTest(inList, keyword, &type, dim, &out->ny);
  g_free(keyword);

  /* "xxxnplane"       olong* Number of planes for each image */
  if (prefix) keyword = g_strconcat (prefix, "nplane", NULL);
  else        keyword = g_strdup("nplane");
  ObitInfoListGetTest(inList, keyword, &type, dim, &out->nplane);
  g_free(keyword);

  /* "xxxinFlysEye"      gboolean* Is each facet in Fly's Eye? */
  if (prefix) keyword = g_strconcat (prefix, "inFlysEye", NULL);
  else        keyword = g_strdup("inFlysEye");
  ObitInfoListGetTest(inList, keyword, &type, dim, &out->inFlysEye);
  g_free(keyword);

  /* "xxxFacetNo"      olong* Untapered facet number - associate various taperings */
  if (prefix) keyword = g_strconcat (prefix, "FacetNo", NULL);
  else        keyword = g_strdup("FacetNo");
  ObitInfoListGetTest(inList, keyword, &type, dim, &out->FacetNo);
  g_free(keyword);

  /* "xxxRAShift"      ofloat* RA shift (deg) for each image */
  if (prefix) keyword = g_strconcat (prefix, "RAShift", NULL);
  else        keyword = g_strdup("RAShift");
  ObitInfoListGetTest(inList, keyword, &type, dim, &out->RAShift);
  g_free(keyword);

  /* "xxxDecShift"     ofloat* Dec shift (deg) for each image */
  if (prefix) keyword = g_strconcat (prefix, "DecShift", NULL);
  else        keyword = g_strdup("DecShift");
  ObitInfoListGetTest(inList, keyword, &type, dim, &out->DecShift);
  g_free(keyword);

  /* "xxxfileType"     olong Are image OBIT_IO_AIPS or OBIT_IO_FITS? */
  if (prefix) keyword = g_strconcat (prefix, "fileType", NULL);
  else        keyword = g_strdup("fileType");
  otemp = 1;
  ObitInfoListGetTest(inList, keyword, &type, dim, &otemp);
  out->fileType = (ObitIOType)otemp;
  g_free(keyword);

  /* "xxximName"       string* Name of Mosaic images */
  if (prefix) keyword = g_strconcat (prefix, "imName", NULL);
  else        keyword = g_strdup("imName");
  ObitInfoListGetTest(inList, keyword, &type, dim, &out->imName);
  g_free(keyword);

  /* "xxximClass"      string* Class of Mosaic images */
  if (prefix) keyword = g_strconcat (prefix, "imClass", NULL);
  else        keyword = g_strdup("imClass");
  ObitInfoListGetTest(inList, keyword, &type, dim, &out->imClass);
  g_free(keyword);

  /* "xxximSeq"        olong Sequence number of Mosaic images */
  if (prefix) keyword = g_strconcat (prefix, "imSeq", NULL);
  else        keyword = g_strdup("imSeq");
  ObitInfoListGetTest(inList, keyword, &type, dim, &out->imSeq);
  g_free(keyword);

  /* "xxximDisk"       olong Disk number of Mosaic images  FIX THIS */
  if (prefix) keyword = g_strconcat (prefix, "imDisk", NULL);
  else        keyword = g_strdup("imDisk");
  ObitInfoListGetTest(inList, keyword, &type, dim, out->imDisk);
  g_free(keyword);

  /* "xxxdoFull"       boolean Is a full field image desired? */
  if (prefix) keyword = g_strconcat (prefix, "doFull", NULL);
  else        keyword = g_strdup("doFull");
  ObitInfoListGetTest(inList, keyword, &type, dim, &out->doFull);
  g_free(keyword);

  /* "xxxbmaj"         ofloat Restoring beam major axis in deg. */
  if (prefix) keyword = g_strconcat (prefix, "bmaj", NULL);
  else        keyword = g_strdup("bmaj");
  ObitInfoListGetTest(inList, keyword, &type, dim, &out->bmaj);
  g_free(keyword);

  /* "xxxbmin"         ofloat Restoring beam minor axis in deg. */
  if (prefix) keyword = g_strconcat (prefix, "bmin", NULL);
  else        keyword = g_strdup("bmin");
  ObitInfoListGetTest(inList, keyword, &type, dim, &out->bmin);
  g_free(keyword);

  /* "xxxbpa"          ofloat Restoring beam PA in deg. */
  if (prefix) keyword = g_strconcat (prefix, "bpa", NULL);
  else        keyword = g_strdup("bpa");
  ObitInfoListGetTest(inList, keyword, &type, dim, &out->bpa);
  g_free(keyword);

  /* "xxxisAuto"       olong Is this an autoCenter field? */
  if (prefix) keyword = g_strconcat (prefix, "isAuto", NULL);
  else        keyword = g_strdup("isAuto");
  ObitInfoListGetTest(inList, keyword, &type, dim, out->isAuto);
  g_free(keyword);

  /* "xxxisShift"       olong Is this an autoCenter shifted field? */
  if (prefix) keyword = g_strconcat (prefix, "isShift", NULL);
  else        keyword = g_strdup("isShift");
  ObitInfoListGetTest(inList, keyword, &type, dim, out->isShift);
  g_free(keyword);

  /* "xxxnumBeamTapes"       olong  Number of Beam tapers (elements in BeamTapes)*/
  if (prefix) keyword = g_strconcat (prefix, "numBeamTapes", NULL);
  else        keyword = g_strdup("numBeamTapes");
  ObitInfoListGetTest(inList, keyword, &type, dim, &out->numBeamTapes);
  g_free(keyword);

  /* "xxxBeamTapes"    ofloat*  List of Beam tapers     */
  /* Add taper array */
  out->BeamTapes = ObitMemAlloc0Name(out->numBeamTapes*sizeof(ofloat),"BeamTapes");
  if (prefix) keyword = g_strconcat (prefix, "BeamTapes", NULL);
  else        keyword = g_strdup("BeamTapes");
  ObitInfoListGetTest(inList, keyword, &type, dim, out->BeamTapes);
  g_free(keyword);

  /* "xxxBeamTaper"      ofloat*  Beam Taper per image as FWHM in deg */
  if (prefix) keyword = g_strconcat (prefix, "BeamTaper", NULL);
  else        keyword = g_strdup("BeamTaper");
  ObitInfoListGetTest(inList, keyword, &type, dim, out->BeamTaper);
  g_free(keyword);

  return out;
} /* end ObitImageMosaicFromInfo */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitImageMosaicGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitImageMosaicClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitGetImageClass */

/**
 * Zap (delete with underlying structures) selected image member(s).
 * Also deletes any associated beam member
 * \param in      The array of images, infoList has following:
 * \li "saveBeam" boolean  If TRUE, save beams, [def FALSE]
 * \param number  The 0-rel image number, -1=> all
 * \param err     Error stack, returns if not empty.
 */
void 
ObitImageMosaicZapImage  (ObitImageMosaic *in, olong number,
			  ObitErr *err)
{
  olong i;
  ObitImage *myBeam=NULL;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gboolean saveBeam = FALSE;
  gchar *routine="ObitImageMosaicZapImage";

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* God save the Beam? */
  ObitInfoListGetTest(in->info, "saveBeam", &type, dim, &saveBeam);

  /* Kill 'em all? */
  if (number==-1) {
    for (i=0; i<in->numberImages; i++) {
      /* First the beam if any */
      if (!saveBeam && in->images[i]->myBeam!=NULL) {
	myBeam = (ObitImage*)in->images[i]->myBeam;
	in->images[i]->myBeam = (Obit*)ObitImageZap(myBeam, err);
      }

      /* then the image */
      in->images[i] = ObitImageZap(in->images[i], err);
      if (err->error) Obit_traceback_msg (err, routine, in->name);
    }
    return;
  }

  /* Check that number in legal range */
  if ((number<0) || (number>in->numberImages)) {
    Obit_log_error(err, OBIT_Error, 
		   "%s: Image number %d out of range [%d %d] in %s", 
		   routine, number, 0, in->numberImages, in->name);
    return;
  }

  /* First the beam if any */
  if (!saveBeam && in->images[number]->myBeam!=NULL) {
    myBeam = (ObitImage*)in->images[number]->myBeam;
    in->images[number]->myBeam = (Obit*)ObitImageZap(myBeam, err);
  }
  /* then the image */
  in->images[number] = ObitImageZap(in->images[number], err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  
  return;
} /* end ObitImageMosaicZapImage */

/**
 * Return pointer to selected image member.
 * Returned reference has been Refed.
 * \param in      The array of images
 * \param number  The 0-rel image number
 * \param err     Error stack, returns if not empty.
 * \return pointer to the selected image. NULL on failure.
 */
ObitImage* 
ObitImageMosaicGetImage  (ObitImageMosaic *in, olong number,
			  ObitErr *err)
{
  ObitImage *out=NULL;
  gchar *routine="ObitImageMosaicGetImage";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return NULL;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Check that number in legal range */
  if ((number<0) || (number>in->numberImages)) {
    Obit_log_error(err, OBIT_Error, 
		   "%s: Image number %d out of range [%d %d] in %s", 
		   routine, number, 0, in->numberImages, in->name);
    return out;
  }

  out = ObitImageRef(in->images[number]);

  return out;
} /* end ObitImageMosaicGetImage */

/**
 * Set a selected image member.
 * \param in      The array of images
 * \param number  The 0-rel image number in array
 * \param image   Image reference to add
 * \param err     Error stack, returns if not empty.
 */
void 
ObitImageMosaicSetImage  (ObitImageMosaic *in, olong number, 
			  ObitImage* image, ObitErr *err)
{
  gchar *routine="ObitImageMosaicSetImage";

  /* error checks */
  g_assert(ObitErrIsA(err));
  g_assert(ObitImageIsA(image));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Check that number in legal range */
  if ((number<0) || (number>in->numberImages)) {
    Obit_log_error(err, OBIT_Error, 
		   "%s: Image number %d out of range [%d %d] in %s", 
		   routine, number, 0, in->numberImages, in->name);
    return;
  }

  /* Only if new and old are different */
  if (image!=in->images[number]) {
    /* Unref any old image*/
    if (in->images[number]) ObitImageUnref(in->images[number]);
    
    in->images[number] = ObitImageRef(image);
  }
} /* end ObitImageMosaicSetImage */


/**
 * Determine the RMS pixel value of image.
 * Ignores outer 5 pixel in image.
 * \param in      The array of images
 * \param number  The 0-rel image number in array
 * \param image   Image reference to add
 * \param err     Error stack, returns if not empty.
 */
ofloat 
ObitImageMosaicGetImageRMS  (ObitImageMosaic *in, olong number,
			     ObitErr *err)
{
  ofloat RMS = -1.0;
  ObitIOSize IOBy;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong blc[IM_MAXDIM], trc[IM_MAXDIM];
  olong i;
  gchar *routine="ObitImageMosaicGetImageRMS";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return RMS;
  g_assert (ObitIsA(in, &myClassInfo));

  for (i=0; i<IM_MAXDIM; i++) blc[i] = 1;
  for (i=0; i<IM_MAXDIM; i++) trc[i] = 0;

 /* Check that number in legal range */
  if ((number<0) || (number>in->numberImages)) {
    Obit_log_error(err, OBIT_Error, 
		   "%s: Image number %d out of range [%d %d] in %s", 
		   routine, number, 0, in->numberImages, in->name);
    return RMS;
  }

  /* Read plane */
  IOBy = OBIT_IO_byPlane;
  dim[0] = 1;
  ObitInfoListPut (in->images[number]->info, "IOBy", OBIT_long, dim, 
		   (gpointer)&IOBy, err);

  /* Set blc, trc to ignore outer 5 pixels */
  blc[0] = 5;
  blc[1] = 5;
  trc[0] = in->images[number]->myDesc->inaxes[0] - 5;
  trc[1] = in->images[number]->myDesc->inaxes[1] - 5;
  dim[0] = 7;
  ObitInfoListPut (in->images[number]->info, "BLC", OBIT_long, dim, blc, err); 
  ObitInfoListPut (in->images[number]->info, "TRC", OBIT_long, dim, trc, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, RMS);
      
 /* Read image */
  ObitImageOpen (in->images[number], OBIT_IO_ReadOnly, err);
  ObitImageRead (in->images[number], NULL, err); /* Read plane */
  ObitImageClose(in->images[number], err);
  if (err->error) Obit_traceback_val (err, routine, in->name, RMS);

  /* Get RMS */
  RMS = ObitFArrayRMS (in->images[number]->image);

  /* Free image buffer */
  in->images[number]->image = ObitFArrayUnref(in->images[number]->image);

  /* Reset blc, trc */
  blc[0] = 1;
  blc[1] = 1;
  trc[0] = 0;
  trc[1] = 0;
  dim[0] = 7;
  ObitInfoListPut (in->images[number]->info, "BLC", OBIT_long, dim, blc, err); 
  ObitInfoListPut (in->images[number]->info, "TRC", OBIT_long, dim, trc, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, RMS);
      
  return RMS;
} /* end ObitImageMosaicGetImageRMS */

/**
 * Return pointer to full field image member.
 * Returned reference has been Refed.
 * \param in      The array of images
 * \param err     Error stack, returns if not empty.
 * \return pointer to the selected image. NULL on failure.
 */
ObitImage* 
ObitImageMosaicGetFullImage  (ObitImageMosaic *in, ObitErr *err)
{
  ObitImage *out=NULL;

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return NULL;
  g_assert (ObitIsA(in, &myClassInfo));

  out = ObitImageRef(in->FullField);

  return out;
} /* end ObitImageMosaicGetFullImage */

/**
 * Set the Full field image member.
 * \param in      The array of images
 * \param image   Image reference to add
 * \param err     Error stack, returns if not empty.
 */
void 
ObitImageMosaicSetFullImage  (ObitImageMosaic *in, 
			      ObitImage* image, ObitErr *err)
{

  /* error checks */
  g_assert(ObitErrIsA(err));
  g_assert(ObitImageIsA(image));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Unref any old image*/
  ObitImageUnref(in->FullField);

  in->FullField = ObitImageRef(image);
} /* end ObitImageMosaicSetFullImage */


/**
 * Make a shallow copy of input object.
 * Parent class members are included but any derived class info is ignored.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Error stack, returns if not empty.
 * \return pointer to the new object.
 */
ObitImageMosaic* ObitImageMosaicCopy (ObitImageMosaic *in, ObitImageMosaic *out, 
				      ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  gboolean oldExist;
  olong i;
  gchar *outName;
  gchar *routine = "ObitImageMosaicCopy";

  /* error checks */
  if (err->error) return NULL;
  g_assert (ObitIsA(in, &myClassInfo));
  if (out) g_assert (ObitIsA(out, &myClassInfo));

  /* Create if it doesn't exist */
  oldExist = out!=NULL;
  if (!oldExist) {
    /* derive object name */
    outName = g_strconcat ("Copy: ",in->name,NULL);
    out = newObitImageMosaic(outName, in->numberImages);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

   /* if Old exists, check that number in legal range */
  if (oldExist) {
    if ((out->numberImages<in->numberImages)) {
      Obit_log_error(err, OBIT_Error, 
		     "%s: Too few images %d in extant output %s", 
		     routine, out->numberImages, out->name);
      return out;
    }

    /* unreference any old members */
    for (i=0; i<out->numberImages; i++)
      out->images[i] = ObitImageUnref(out->images[i]);
  }

  /* Copy other class members */
  out->xCells   = in->xCells;
  out->yCells   = in->yCells;
  out->fileType = in->fileType;
  out->imSeq    = in->imSeq;
  if (out->imName) g_free(out->imName);
  out->imName = g_strdup(in->imName);
  if (out->imClass) g_free(out->imClass);
  out->imClass = g_strdup(in->imClass);
  out->bmaj    = in->bmaj;
  out->bmin    = in->bmin;
  out->bpa     = in->bpa;
  for (i=0; i<in->numberImages; i++) {
    ObitImageMosaicSetImage(out, i, in->images[i], err);
    out->imDisk[i]   = in->imDisk[i];
    out->isAuto[i]   = in->isAuto[i];
    out->isShift[i]  = in->isShift[i];
    out->inFlysEye[i]= in->inFlysEye[i];
    out->FacetNo[i]  = in->FacetNo[i];
    out->nx[i]       = in->nx[i];
    out->ny[i]       = in->ny[i];
    out->nplane[i]   = in->nplane[i];
    out->RAShift[i]  = in->RAShift[i];
    out->DecShift[i] = in->DecShift[i];
    out->BeamTaper[i]= in->BeamTaper[i];
  }
  out->numBeamTapes = in->numBeamTapes;
  if (out->BeamTapes==NULL) 
    out->BeamTapes = ObitMemAllocName(in->numBeamTapes*sizeof(ofloat), "BeamTapes");
  for (i=0; i<out->numBeamTapes; i++) {
    out->BeamTapes[i]= in->BeamTapes[i];
  }
  if (err->error) Obit_traceback_val (err, routine, in->name, out);

  return out;
} /* end ObitImageMosaicCopy */

/**
 * Attach images to their underlying files.
 * For AIPS files:
 * Mosaic images have classes with the first character of the imClass
 * followed by 'M', followed by 4 digits of the field number.  
 * The Beams are the same except that the second character is 'B'.
 * The full field image has class imClass and 'F' as the 6th character.
 *
 * For FITS files:
 * Image classes are imClass+digits of the field
 * Beam classes are the same except the second character is replaced 
 * with a 'B' (unless it already is 'B' in which case 'b' is used).
 * The full field image has class imClass followed by 'Full'
 * 
 * \param in     The object with images,  Details are defined in members:
 * \li numberImage - Number of images in Mosaic
 * \li nInit    - number of images already initialized
 * \li images   - Image array
 * \li fileType - Are image OBIT_IO_AIPS or OBIT_IO_FITS?
 * \li imName   - Name of Mosaic images
 * \li imClass  - imClass
 * \li imSeq    - imSeq
 * \li imDisk   - imDisk one per image
 * \param doBeam  If true, make beam as well.
 * \param err     Error stack, returns if not empty.
 */
void ObitImageMosaicSetFiles  (ObitImageMosaic *in, gboolean doBeam, ObitErr *err) 
{
  olong i, j, user, cno;
  olong blc[IM_MAXDIM] = {1,1,1,1,1};
  olong trc[IM_MAXDIM] = {0,0,0,0,0};
  gboolean exist;
  ObitImage *image=NULL;
  gchar ct, strTemp[48], Aname[13], Aclass[7], Atclass[3], Atype[3] = "MA";
  gchar *routine = "ObitImageMosaicSetFiles";

 /* Create full field image if needed */
  /*???if (in->doFull && (in->FOV>0.0) && (in->numberImages>1) && (in->nInit<=0)) {*/
  if (in->doFull && (in->numberImages>1) && (in->nInit<=0)) {
    image = in->FullField;

    /* AIPS file */
    if (in->fileType==OBIT_IO_AIPS) {
      /* Get AIPS user number */
      user = ObitSystemGetAIPSuser();
      /* set class */
      sprintf(strTemp, "%s", in->imClass);
      strncpy (Aclass, strTemp,    6);  Aclass[6] = 0;
      Aclass[5] = 'F';
      strncpy (Aname,  in->imName, 12); Aname[12] = 0;
      /* allocate */
      cno = ObitAIPSDirAlloc(in->imDisk[0], user, Aname, Aclass, Atype, in->imSeq, &exist, err);
      /* Set info */
      ObitImageSetAIPS(image,OBIT_IO_byPlane,in->imDisk[0],cno,user,blc,trc,err);
      /*DEBUG */
      Obit_log_error(err, OBIT_InfoErr, "Making AIPS image %s %s %d on disk %d cno %d",
		     Aname, Aclass, in->imSeq, in->imDisk[0], cno);
      /* fprintf (stderr,"Making AIPS image %s %s on disk %d cno %d\n",
	 Aname, Aclass, in->imDisk[0], cno);*/

    /* FITS file */
    } else if (in->fileType==OBIT_IO_FITS) {
      /* set filename - derive from field */
      sprintf(strTemp, "MA%s.%sFulldseq%d", in->imName, in->imClass, in->imSeq);
      /* replace blanks with underscores */
      for (j=0; j<strlen(strTemp); j++) if (strTemp[j]==' ') strTemp[j]='_';
      ObitImageSetFITS(image,OBIT_IO_byPlane,in->imDisk[0],strTemp,blc,trc,err);
    }
    if (err->error) Obit_traceback_msg (err, routine, in->name);

    /* fully instantiate */
    ObitImageFullInstantiate (image, FALSE, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);

  } /* end do full field image */

  /* Loop over images */
  for (i=in->nInit; i<in->numberImages; i++) {

    /* Define files */
    image = in->images[i];
    if (in->fileType==OBIT_IO_AIPS) {
      /* Get AIPS user number */
      user = ObitSystemGetAIPSuser();
      /* set class */
      /* Only one character of Class allowed here*/
      Atclass[0] = in->imClass[0]; Atclass[1] = 'M'; Atclass[2] = 0;
      sprintf(strTemp, "%s%4.4d", Atclass, i+1);
      strncpy (Aclass, strTemp,    6);  Aclass[6] = 0;
      strncpy (Aname,  in->imName, 12); Aname[12] = 0;
      /* allocate */
      cno = ObitAIPSDirAlloc(in->imDisk[i], user, Aname, Aclass, Atype, in->imSeq, &exist, err);
      /* Set info */
      ObitImageSetAIPS(image,OBIT_IO_byPlane,in->imDisk[i],cno,user,blc,trc,err);
      /*DEBUG */
      Obit_log_error(err, OBIT_InfoErr, "Making AIPS image %s %s %d on disk %d cno %d",
		     Aname, Aclass, in->imSeq, in->imDisk[i], cno);
      /* fprintf (stderr,"Making AIPS image %s %s on disk %d cno %d\n",
	 Aname, Aclass, in->imDisk[i], cno);*/
    } else if (in->fileType==OBIT_IO_FITS) {
      /* set filename - derive from field */
      Atclass[0] = in->imClass[0]; Atclass[1] =  'M'; Atclass[2] = 0;
      sprintf(strTemp, "MA%s.%s%4.4dseq%d", in->imName, Atclass, i, in->imSeq);
      /* replace blanks with underscores */
      for (j=0; j<strlen(strTemp); j++) if (strTemp[j]==' ') strTemp[j]='_';
      ObitImageSetFITS(image,OBIT_IO_byPlane,in->imDisk[i],strTemp,blc,trc,err);
    }
    /* fully instantiate */
    ObitImageFullInstantiate (image, FALSE, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
 
    /* Doing beams ?*/
    if (doBeam) {
      image = (ObitImage*)in->images[i]->myBeam;
      /* Define files - same except second character of "class" is 'B' */
      if (in->fileType==OBIT_IO_AIPS) {
	/* Get AIPS user number */
	user = ObitSystemGetAIPSuser();
	/* set class */
	Aclass[1] = 'B';
	/* allocate */
	cno = ObitAIPSDirAlloc(in->imDisk[i], user, Aname, Aclass, Atype, in->imSeq, &exist, err);
	/* Set info */
	ObitImageSetAIPS(image,OBIT_IO_byPlane,in->imDisk[i],cno,user,blc,trc,err);
	
      } else if (in->fileType==OBIT_IO_FITS) {
	/* set filename - derive from field - insure uniqueness */
	ct = in->imClass[1];
	if (in->imClass[1] != 'B') in->imClass[1] = 'B';
	else in->imClass[1] = 'b';
	Atclass[0] = in->imClass[0]; Atclass[1] =  in->imClass[1]; Atclass[2] = 0;
	sprintf(strTemp, "MA%s.%s%4.4dseq%d", in->imName,  Atclass, i, in->imSeq);
	/* replace blanks with underscores */
	for (j=0; j<strlen(strTemp); j++) if (strTemp[j]==' ') strTemp[j]='_';
	in->imClass[1] = ct;
	ObitImageSetFITS(image,OBIT_IO_byPlane,in->imDisk[i],strTemp,blc,trc,err);
      }
      /* Open and close to fully instantiate */
      ObitImageOpen(image, OBIT_IO_WriteOnly, err);
      ObitImageClose(image, err);
      if (err->error) Obit_traceback_msg (err, routine, in->name);
    } /* end Beams */
  }    /* end loop over images */

  /* Everything should now be initialized */
  in->nInit = in->numberImages;
 
} /* end  ObitImageMosaicSetFiles */

/**
 * Create an Image Mosaic based on a uv data and parameters attached thereto
 * 
 * \param name    Name to be given to new object
 * \param uvData  The object to create images in,  Details are defined in info members:
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
 * \li numBeamTapes= Number of Beam tapers (elements in BeamTapes) 
 * \li BeamTapes   = List of additional Beam tapers  larger than 0
 * \param err     Error stack, returns if not empty.
 * \return Newly created object.
 */
ObitImageMosaic* ObitImageMosaicCreate (gchar *name, ObitUV *uvData, ObitErr *err)
{
  ObitImageMosaic *out = NULL;
  ObitInfoType type;
  ObitIOType Type;
  ObitIOAccess access;
  gint32 i, nf, nif, nfif, dim[MAXINFOELEMDIM];
  olong Seq, *Disk=NULL, NField, nx[MAXFLD], ny[MAXFLD], catDisk, nDisk;
  olong overlap, imsize, fldsiz[MAXFLD], flqual[MAXFLD];
  ofloat FOV=0.0, xCells, yCells, RAShift[MAXFLD], DecShift[MAXFLD], Tapers[MAXFLD], 
    MaxBL, MaxW, Cells, BeamTapes[100], diam=24.5;
  ofloat *farr, Radius = 0.0, maxScale;
  ofloat shift[2] = {0.0,0.0}, cells[2], bmaj, bmin, bpa, beam[3];
  odouble ra0, dec0, refFreq1, refFreq2;
  gboolean doJ2B, doFull, doCalSelect, inFlysEye[MAXFLD], FacetNo[MAXFLD];
  ofloat equinox, minRad, ratio, OutlierFlux, OutlierDist, OutlierSI;
  olong  itemp, OutlierSize=0, nFlyEye = 0, *iarr=NULL, numBeamTapes=0;
  union ObitInfoListEquiv InfoReal; 
  gchar Catalog[100], Aname[100], Aclass[100];
  gchar *routine = "ObitImageMosaicCreate";

  /* Get inputs with plausible defaults */
  NField = 1;
  ObitInfoListGetTest(uvData->info, "NField",  &type, dim,  &NField);
  ObitInfoListGet(uvData->info, "imFileType",&type, dim,  &Type,     err);
  ObitInfoListGet(uvData->info, "imName",    &type, dim,  Aname,     err);
  ObitInfoListGet(uvData->info, "imClass",   &type, dim,  Aclass,    err);
  ObitInfoListGet(uvData->info, "imSeq",     &type, dim,  &Seq,      err);
  nDisk = MAX(1,NField);
  Disk = g_malloc(nDisk*sizeof(olong));
  ObitInfoListGetP(uvData->info, "imDisk",    &type, dim, (gpointer)&iarr);
  /* If only one disk given - use for all */
  Disk[0] = iarr[0];
  if (dim[0]<=1) {
    for (i=1; i<nDisk; i++) Disk[i] = Disk[0];
  }
  ObitInfoListGet(uvData->info, "FOV",     &type, dim,  &InfoReal, err);
  if (type==OBIT_float) FOV = InfoReal.flt;
  else if (type==OBIT_double)  FOV = InfoReal.dbl;
  xCells = 0.0;
  /*ObitInfoListGetTest(uvData->info, "xCells",  &type, dim,  &xCells);*/
  ObitInfoListGetP(uvData->info, "xCells",  &type, dim,  (gpointer)&farr);
  if (farr!=NULL) xCells = farr[0];
  yCells = 0.0;
  /*ObitInfoListGetTest(uvData->info, "yCells",  &type, dim,  &yCells);*/
  ObitInfoListGetP(uvData->info, "yCells",  &type, dim, (gpointer)&farr);
  if (farr!=NULL) yCells = farr[0];
  doFull = FALSE;
  ObitInfoListGetTest(uvData->info, "doFull", &type, dim,  &doFull);
  sprintf (Catalog, "None");
  ObitInfoListGetTest(uvData->info, "Catalog", &type, dim,  Catalog);
  catDisk = 1;  /* FITS directory for catalog */
  ObitInfoListGetTest(uvData->info, "CatDisk", &type, dim, &catDisk);
  if (err->error) Obit_traceback_val (err, routine, uvData->name, out);

  /* Get array inputs */
  for (i=0; i<MAXFLD; i++) nx[i] = 0;
  ObitInfoListGetTest(uvData->info, "nx",       &type, dim, nx);
  for (i=0; i<MAXFLD; i++) ny[i] = 0;
  ObitInfoListGetTest(uvData->info, "ny",       &type, dim, ny);
  for (i=0; i<MAXFLD; i++) RAShift[i] = 0.0;
  ObitInfoListGetTest(uvData->info, "RAShift",  &type, dim, RAShift);
  for (i=0; i<MAXFLD; i++) DecShift[i] = 0.0;
  ObitInfoListGetTest(uvData->info, "DecShift", &type, dim, DecShift);
  /* Optional Beam */
  bmaj = bmin = bpa = 0.0;
  ObitInfoListGetTest(uvData->info, "BMAJ", &type, dim, &bmaj);
  ObitInfoListGetTest(uvData->info, "BMIN", &type, dim, &bmin);
  ObitInfoListGetTest(uvData->info, "BPA",  &type, dim, &bpa);
  /* Try alternate form - all in beam */
  beam[0] = bmaj; beam[1] = bmin; beam[2] = bpa;
  ObitInfoListGetTest(uvData->info, "Beam",  &type, dim, beam);
  bmaj = beam[0]; bmin = beam[1]; bpa = beam[2];
  /* Beam given in asec - convert to degrees */
  bmaj /= 3600.0; bmin /=3600.0;

  /* Get antenna diameter for non-VLA antennas. */
  if (!strncmp(uvData->myDesc->teles, "KAT-7",5)) diam = 12.0;  /* KAT */

  /* Get extrema - note: this has no selection */
  ObitUVUtilUVWExtrema (uvData, &MaxBL, &MaxW, err);
  if (err->error) Obit_traceback_val (err, routine, uvData->name, out);

  /* Get reference frequency without freq selection */
  refFreq1 = uvData->myDesc->crval[uvData->myDesc->jlocf];

  /* Open with any selection of input */
  doCalSelect = FALSE;
  ObitInfoListGetTest(uvData->info, "doCalSelect", &type, dim, &doCalSelect);
  if (doCalSelect) access = OBIT_IO_ReadCal;
  else access = OBIT_IO_ReadOnly;
  if (doCalSelect) {
    ObitUVOpen (uvData, access, err);
    ObitUVClose (uvData, err);
    if (err->error) Obit_traceback_val (err, routine, uvData->name, out);
  }

  /* Get reference frequency with freq selection */
  refFreq2 = uvData->myDesc->crval[uvData->myDesc->jlocf];

  /* Scale MaxBL, MaxW for selection */
  MaxBL *= refFreq2/refFreq1;
  MaxW  *= refFreq2/refFreq1;

  /* Find maximum uv scaling */
  maxScale = 0.0;
  /* how big is table */
  nf = 1;
  if (uvData->myDesc->jlocf>=0) nf = MAX (1, uvData->myDesc->inaxes[uvData->myDesc->jlocf]);
  nif = 1;
  if (uvData->myDesc->jlocif>=0) nif = MAX (1, uvData->myDesc->inaxes[uvData->myDesc->jlocif]);
  nfif = nf*nif;
  if (uvData->myDesc->fscale) {
    for (i=0; i<nfif; i++) maxScale = MAX (maxScale, uvData->myDesc->fscale[i]);
  } else {
    maxScale = 1.0;
  }
  /* Want MaxBL, MaxW, at highest frequency */
  MaxBL *= maxScale;
  MaxW  *= maxScale;

  /* Suggested cellsize and facet size */
  Cells = 0.0;
  ObitImageUtilImagParm (MaxBL, MaxW, &Cells, &Radius);

 /* Check cell spacing if given */
  if ((fabs(xCells)>0.0) || (fabs(yCells)>0.0)) {
    ratio = fabs(Cells) / fabs(xCells);
    Obit_retval_if_fail(((ratio<10.0) && (ratio>0.1)), err, out,
			"%s: Cellsize seriously bad, suggest %g asec", 
			routine, Cells);
    Radius *= ratio;  /* Correct Radius to actual cell size */
  }

  /* Get cells spacing and maximum undistorted radius from uv data if needed */
  if ((FOV>0.0) || (xCells<=0.0) || (yCells<=0.0)) {
    if (xCells==0.0) xCells = -Cells;
    if (yCells==0.0) yCells =  Cells;
    /* tell about it */
    Obit_log_error(err, OBIT_InfoErr, 
		   "Suggested cell spacing %f, undistorted FOV %f cells for %s",
		   Cells, Radius, uvData->name);
  } else { /* Use given FOV */
    Radius = 0.0;
  }
     
  /* Set fly's eye if needed */
  imsize = (olong)(2.0*Radius + 30.99);
  imsize = MAX (64, imsize);  /* Not too small */
  
  /* Not bigger than FOV */
  imsize = MIN (imsize, ((2.0*3600.0*FOV/fabs(xCells))+10.99));
  overlap = 7;
  cells[0] = xCells; cells[1] = yCells;
  ObitUVGetRADec (uvData, &ra0, &dec0, err);
  if (err->error) Obit_traceback_val (err, routine, uvData->name, out);
  /* Copy nx to fldsiz - Make image sizes FFT friendly */
  for (i=0; i<NField; i++) fldsiz[i] = ObitFFTSuggestSize (MAX(32,nx[i]));
  minRad = 0.0;

  /* Get Fly's eye */
  if (FOV>0.0) minRad = FlyEye(FOV, imsize, cells, overlap, shift, ra0, dec0, 
			       &NField, fldsiz, RAShift, DecShift, flqual, err); 
  if (err->error) Obit_traceback_val (err, routine, uvData->name, out);
  nFlyEye = NField;  /* Number in Fly's eye */
  for (i=0; i<nFlyEye; i++) inFlysEye[i] = TRUE;  /* Mark as in fly's eye */
  for (i=0; i<nFlyEye; i++) FacetNo[i]   = i;  /* Untapered Facet number */
  for (i=nFlyEye; i<MAXFLD; i++) inFlysEye[i] = FALSE; /* Rest not */
  for (i=nFlyEye; i<MAXFLD; i++) FacetNo[i]   = 0;

  /* Add outlyers from catalog if requested */
  /* Blank = default */
  if (!strncmp(Catalog, "    ", 4)) sprintf (Catalog, "Default");
  if (strncmp(Catalog, "None", 4)) {
    /* Set default catalog */
     if (!strncmp(Catalog, "Default", 7)) sprintf (Catalog, "NVSSVZ.FIT");

     /* Get outlier related inputs */
     OutlierDist = 1.0;
     ObitInfoListGetTest(uvData->info, "OutlierDist",  &type, dim,  &OutlierDist);
     OutlierFlux = 0.1;
     ObitInfoListGetTest(uvData->info, "OutlierFlux",  &type, dim,  &OutlierFlux);
     OutlierSI = -0.75;
     ObitInfoListGetTest(uvData->info, "OutlierSI",    &type, dim,  &OutlierSI);
     InfoReal.itg = 50;type = OBIT_float;
     ObitInfoListGetTest(uvData->info, "OutlierSize",  &type, dim,  &InfoReal);
     if (type==OBIT_float) itemp = InfoReal.flt + 0.5;
     else itemp = InfoReal.itg;
     OutlierSize = itemp;

     /* Add to list from catalog */
     equinox = uvData->myDesc->equinox;  /* Clear up confusion in AIPS */
     if (equinox<1.0) equinox = uvData->myDesc->epoch;
     doJ2B = (equinox!=2000.0) ;  /* need to precess? */

     AddOutlier (Catalog, catDisk, minRad, cells, OutlierDist, OutlierFlux, OutlierSI, OutlierSize,
		 ra0, dec0, doJ2B, uvData->myDesc->crval[uvData->myDesc->jlocf], Radius, diam, 
		 &NField, fldsiz, RAShift, DecShift, flqual, err);
     if (err->error) Obit_traceback_val (err, routine, uvData->name, out);
     /* Assign facet numbers */
     for (i=nFlyEye; i<NField; i++) FacetNo[i]   = i;
  } /* end add outliers from catalog */

  /* Make sure some fields defined */
  Obit_retval_if_fail((NField>0), err, out, "%s: NO Fields defined", routine);

  /* Add any additional tapered images */
  AddTapers (uvData, &numBeamTapes, BeamTapes, fabs(cells[0]), Tapers, 
	     &NField, fldsiz, RAShift, DecShift, flqual, inFlysEye, FacetNo, err);
  if (err->error) Obit_traceback_val (err, routine, uvData->name, out);

  /* Copy fldsiz to nx, ny */
  for (i=0; i<NField; i++) nx[i] = ny[i] = fldsiz[i];
  /*for (i=0; i<NField; i++) nx[i] = ny[i] = 1024;  debug */
      
  /* Create output object */
  out = newObitImageMosaic (name, NField);

  /* Add taper array */
  out->BeamTapes = ObitMemAlloc0Name(numBeamTapes*sizeof(ofloat),"BeamTapes");

  /* Set values on out */
  out->fileType= Type;
  Aname[12] = 0; Aclass[6]=0;
  out->imName  = g_strdup(Aname);
  out->imClass = g_strdup(Aclass);
  out->imSeq   = Seq;
  out->xCells  = -fabs(xCells)/3600.0;  /* to Deg*/
  out->yCells  = yCells/3600.0;         /* to deg */
  out->FOV     = FOV;
  out->Radius  = Radius;
  out->doFull  = doFull;
  out->bmaj    = bmaj;
  out->bmin    = bmin;
  out->bpa     = bpa;
  out->nFlyEye = nFlyEye;
  out->OutlierSize  = OutlierSize;
  out->numBeamTapes = numBeamTapes;
  for (i=0; i<numBeamTapes; i++) out->BeamTapes[i] = BeamTapes[i];
  for (i=0; i<NField; i++) {
    if (i<nDisk) out->imDisk[i]   = Disk[i];
    else out->imDisk[i]   = Disk[0];
    out->isAuto[i]   = -1;
    out->isShift[i]  = -1;
    out->nx[i]       = nx[i];
    out->ny[i]       = ny[i];
    out->inFlysEye[i]= inFlysEye[i];
    out->FacetNo[i]  = FacetNo[i];
    out->RAShift[i]  = RAShift[i]/3600.0;   /* to Deg */
    out->DecShift[i] = DecShift[i]/3600.0;  /* to Deg */
    out->BeamTaper[i]= Tapers[i];
  }

  /* Cleanup */
  if (Disk) g_free(Disk);

  return out;
} /* end ObitImageMosaicCreate */

/**
 * Define images in an Image Mosaic from a UV data.
 * Saves BeamTaper (as BeamTapr) in Descriptor infoList
 * \param in     The object to create images in,  Details are defined in members:
 * \li numberImage - Number of images in Mosaic
 * \li nInit    - number of images already initialized
 * \li images   - Image array
 * \li xCells   - Cell Spacing in X (deg)
 * \li yCells   - Cell Spacing in Y (deg)
 * \li nx       - Number of pixels in X for each image
 * \li ny       - Number of pixels in Y for each image
 * \li nplane   - Number of planes for each image
 * \li RAShift  - RA shift (asec) for each image
 * \li DecShift - Dec shift (asec) for each image
 * \li fileType - Are image OBIT_IO_AIPS or OBIT_IO_FITS?
 * \li imName   - Name of Mosaic images
 * \li imClass  - imClass
 * \li imSeq    - imSeq
 * \li imDisk   - imDisk
 * \li isAuto   - if >0 this is an autoCenter field and the value is 
 *                any corresponding 2D shifted (aligned) field, per field
 *                For 3D imaging any positive value indicates an autoCenter field.
 * \li isShift  - if >0 this is an autoCenter shifted field and the value is 
 *                the corresponding 2D autoCenter field, per field
 * \param uvData UV data from which the images are to be made.
 * Imaging information on uvData:
 * \li "nChAvg" OBIT_int (1,1,1) number of channels to average.
 *              This is for spectral line observations and is ignored
 *              if the IF axis on the uv data has more than one IF.
 *              Default is continuum = average all freq/IFs. 0=> all.
 * \li "rotate" OBIT_float (?,1,1) Desired rotation on sky (from N thru E) in deg. [0]
 * \param doBeam  If true, make beam as well.
 * \param err     Error stack, returns if not empty.
 */
void ObitImageMosaicDefine (ObitImageMosaic *in, ObitUV *uvData, gboolean doBeam,
			    ObitErr *err)
{
  olong i, nx, ny;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ofloat *farr=NULL, FOV;
  gboolean doCalSelect, *barr=NULL;
  ObitIOAccess access;
  ObitImage *myBeam;
  ObitImageMosaicClassInfo* mosaicClass;
  gchar *routine = "ObitImageMosaicDefine";

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitUVIsA(uvData));

  /* Set values on uv data */
  Obit_return_if_fail((in->numberImages>=1), err,
		      "%s: No images requested", routine);
  dim[0] = in->numberImages;
  ObitInfoListAlwaysPut (uvData->info, "nx",     OBIT_long,   dim, in->nx);
  ObitInfoListAlwaysPut (uvData->info, "nxBeam", OBIT_long,   dim, in->nx);
  ObitInfoListAlwaysPut (uvData->info, "ny",     OBIT_long,   dim, in->ny);
  ObitInfoListAlwaysPut (uvData->info, "nyBeam", OBIT_long,   dim, in->ny);
  ObitInfoListAlwaysPut (uvData->info, "xShift", OBIT_float,  dim, in->RAShift);
  ObitInfoListAlwaysPut (uvData->info, "yShift", OBIT_float,  dim, in->DecShift);
  farr = ObitMemAllocName(in->numberImages*sizeof(ofloat), routine);
  for (i=0; i<in->numberImages; i++) farr[i] = 3600.0*in->xCells;
  ObitInfoListAlwaysPut (uvData->info, "xCells", OBIT_float, dim, farr);
  for (i=0; i<in->numberImages; i++) farr[i] = 3600.0*in->yCells;
  ObitInfoListAlwaysPut (uvData->info, "yCells", OBIT_float, dim, farr);
  ObitMemFree(farr);
  /* Set doGrid=isAuto<=0 - only meaningful if do3D=FALSE */
  barr = ObitMemAllocName(in->numberImages*sizeof(gboolean), routine);
  for (i=0; i<in->numberImages; i++) barr[i] = in->isAuto[i]<=0;
  ObitInfoListAlwaysPut (uvData->info, "doGrid", OBIT_float, dim, barr);
  ObitMemFree(barr);
 
  /* Make sure UV data descriptor has proper info */
  doCalSelect = FALSE;
  ObitInfoListGetTest(uvData->info, "doCalSelect", &type, (gint32*)dim, &doCalSelect);
  if (doCalSelect) access = OBIT_IO_ReadCal;
  else access = OBIT_IO_ReadOnly;

  /* Open/Close to update descriptor */
  ObitUVOpen (uvData, access, err);
  ObitUVClose (uvData, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Loop over uninitialized images */
  for (i=in->nInit; i<in->numberImages; i++) {
    in->images[i] = ObitImageUtilCreateImage(uvData, i+1, doBeam, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    /* For 2D imaging pick up shift */
    if (!in->images[i]->myDesc->do3D) {
      in->RAShift[i]  = in->images[i]->myDesc->xshift;
      in->DecShift[i] = in->images[i]->myDesc->yshift;
    }
    /* Add BeamTaper to Descriptor InfoList - only 8 char allowed */
    dim[0] = dim[1] = 1;
    ObitInfoListAlwaysPut (in->images[i]->myDesc->info, "BeamTapr", OBIT_float, 
			   dim, &in->BeamTaper[i]);
    myBeam = (ObitImage*)(in->images[i]->myBeam);
    ObitInfoListAlwaysPut (myBeam->myDesc->info, "BeamTapr", OBIT_float, 
			   dim, &in->BeamTaper[i]);
  }    /* end loop over images */

 /* Create full field image if needed */
  if (in->doFull && (in->nInit<=0) && (in->numberImages>1)) { 
    /* Basic structure of field 1 */
    in->FullField = ObitImageUtilCreateImage(uvData, 1, FALSE, err);
    /* Replace name */
    g_free(in->FullField->name);
    in->FullField->name = g_strdup("Flattened Image");
    mosaicClass = (ObitImageMosaicClassInfo*)&myClassInfo;  /* Class pointer */
    FOV = mosaicClass->ObitImageMosaicFOV (in, err);
    /* set size */
    nx = 2 * FOV / fabs (in->xCells);
    ny = 2 * FOV / fabs (in->yCells);
    in->FullField->myDesc->inaxes[0] = nx;
    in->FullField->myDesc->inaxes[1] = ny;
    /* set center pixel must be an integer pixel */
    in->FullField->myDesc->crpix[0] = 1 + nx/2;
    in->FullField->myDesc->crpix[1] = 1 + ny/2;
    /* set shift to zero */
    in->FullField->myDesc->xshift = 0.0;
    in->FullField->myDesc->yshift = 0.0;
  }

  ObitImageMosaicSetFiles (in, doBeam, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
 
} /* end ObitImageMosaicDefine */


/**
 * Add a Field to a mosaic
 * \param in       The object with images to modify
 * \param uvData   UV data from which the images are to be made.
 * \param nx       Number of pixels in X for image, 
 *                 will be converted to next larger good FFT size.
 * \param ny       Number of pixels in Y for image
 * \param nplane   Number of planes for image
 * \param RAShift  RA shift (asec) for image
 * \param DecShift Dec shift (asec) for image
 * \param isAuto   If true, this is an autoCenter image
 * \param err      Error stack, returns if not empty.
 */
void ObitImageMosaicAddField (ObitImageMosaic *in, ObitUV *uvData, 
			      olong nx, olong ny, olong nplane, 
			      ofloat RAShift, ofloat DecShift,
			      gboolean isAuto, ObitErr *err)
{
  ofloat    *ftemp;
  olong     i, *itemp;
  gboolean  *btemp;
  ObitImage **imtemp;
  gchar *routine = "ObitImageMosaicAddField";

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Expand arrays */
  in->nInit = in->numberImages;  /* To ensure new field gets initialized */
  in->numberImages++;

  /* Image array */
  imtemp = ObitMemAlloc0Name(in->numberImages*sizeof(ObitImage*),"ImageMosaic images");
  for (i=0; i<in->nInit; i++) imtemp[i] = in->images[i]; imtemp[i] = NULL;
  in->images = ObitMemFree(in->images);
  in->images = imtemp;

  /* Disk array - add new images on same disk as first */
  itemp  = ObitMemAlloc0Name(in->numberImages*sizeof(olong),"ImageMosaic imDisk");
  for (i=in->nInit; i<in->numberImages; i++) itemp[i] = in->imDisk[0]; 
  for (i=0; i<in->nInit; i++) itemp[i] = in->imDisk[i]; 
  in->imDisk = ObitMemFree(in->imDisk);
  in->imDisk = itemp;
 
  /* isAuto array */
  itemp  = ObitMemAlloc0Name(in->numberImages*sizeof(olong),"ImageMosaic isAuto");
  for (i=0; i<in->nInit; i++) itemp[i] = in->isAuto[i]; 
  for (i=in->nInit; i<in->numberImages; i++) itemp[i] = -1; 
  /* Large dummy value in new element if isAuto */
  if (isAuto) itemp[in->numberImages-1] = in->numberImages+9999999;
  in->isAuto = ObitMemFree(in->isAuto);
  in->isAuto = itemp;
 
  /* isShift array */
  itemp  = ObitMemAlloc0Name(in->numberImages*sizeof(olong),"ImageMosaic isShift");
  for (i=0; i<in->nInit; i++) itemp[i] = in->isShift[i]; 
  for (i=in->nInit; i<in->numberImages; i++) itemp[i] = -1; 
  in->isShift = ObitMemFree(in->isShift);
  in->isShift = itemp;
 
  /* inFlysEye array - any new entries = FALSE */
  btemp  = ObitMemAlloc0Name(in->numberImages*sizeof(gboolean),"ImageMosaic inFlysEye");
  for (i=0; i<in->nInit; i++) btemp[i] = in->inFlysEye[i]; 
  for (i=in->nInit; i<in->numberImages; i++) btemp[i] = FALSE; 
  in->inFlysEye = ObitMemFree(in->inFlysEye);
  in->inFlysEye = btemp;
 
  /* FacetNo array - any new entries = number in mosaic */
  itemp  = ObitMemAlloc0Name(in->numberImages*sizeof(olong),"ImageMosaic FacetNo");
  for (i=0; i<in->nInit; i++) itemp[i] = in->FacetNo[i]; 
  for (i=in->nInit; i<in->numberImages; i++) itemp[i] = i; 
  in->FacetNo = ObitMemFree(in->FacetNo);
  in->FacetNo = itemp;
 
  /* Image size */
  itemp  = ObitMemAlloc0Name(in->numberImages*sizeof(olong),"ImageMosaic nx");
  for (i=0; i<in->nInit; i++) itemp[i] = in->nx[i]; 
  itemp[i] = ObitFFTSuggestSize (nx);
  in->nx = ObitMemFree(in->nx);
  in->nx = itemp;
  itemp  = ObitMemAlloc0Name(in->numberImages*sizeof(olong),"ImageMosaic ny");
  for (i=0; i<in->nInit; i++) itemp[i] = in->ny[i]; 
  itemp[i] = ObitFFTSuggestSize (ny);
  in->ny = ObitMemFree(in->ny);
  in->ny = itemp;
  itemp  = ObitMemAlloc0Name(in->numberImages*sizeof(olong),"ImageMosaic nplane");
  for (i=0; i<in->nInit; i++) itemp[i] = in->nplane[i]; 
  itemp[i] = nplane;
  in->nplane = ObitMemFree(in->nplane);
  in->nplane = itemp;

  /* Shift */
  ftemp  = ObitMemAlloc0Name(in->numberImages*sizeof(ofloat),"ImageMosaic RAShift");
  for (i=0; i<in->nInit; i++) ftemp[i] = in->RAShift[i]; 
  ftemp[i] = RAShift;
  in->RAShift = ObitMemFree(in->RAShift);
  in->RAShift = ftemp;
  ftemp  = ObitMemAlloc0Name(in->numberImages*sizeof(ofloat),"ImageMosaic DecShift");
  for (i=0; i<in->nInit; i++) ftemp[i] = in->DecShift[i]; 
  ftemp[i] = DecShift;
  in->DecShift = ObitMemFree(in->DecShift);
  in->DecShift = ftemp;

  /* BeamTaper - new full resolution */
  ftemp  = ObitMemAlloc0Name(in->numberImages*sizeof(ofloat),"ImageMosaic BeamTaper");
  for (i=0; i<in->nInit; i++) ftemp[i] = in->BeamTaper[i]; 
  ftemp[i] = 0.0;
  in->BeamTaper = ObitMemFree(in->BeamTaper);
  in->BeamTaper = ftemp;

  /* Define image */
  ObitImageMosaicDefine (in, uvData, TRUE, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Create/initialize image */
  ObitImageMosaicSetFiles (in, TRUE, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

} /* end ObitImageMosaicAddField */

/**
 * Project the tiles of a Mosaic to the full field flattened image
 * Ignore facets with BeamTaper>0
 * Routine translated from the AIPSish 4MASS/SUB/FLATEN.FOR/FLATEN  
 * \param in     The object with images
 * Control info on info member:
 * \li "planeNo" OBIT_long (1,1,1) plane (1-rel) in mosiac to flatten [def 1]
 * \param err     Error stack, returns if not empty.
 */
void ObitImageMosaicFlatten (ObitImageMosaic *in, ObitErr *err)
{
  olong   blc[IM_MAXDIM]={1,1,1,1,1}, trc[IM_MAXDIM]={0,0,0,0,0};
  olong i, j, radius, rad, plane[IM_MAXDIM] = {1,1,1,1,1}, hwidth = 2;
  olong *naxis, pos1[IM_MAXDIM], pos2[IM_MAXDIM], planeNo;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitImage *out=NULL, *tout1=NULL, *tout2=NULL;
  ObitFArray *sc2=NULL, *sc1=NULL;
  ObitIOSize IOsize = OBIT_IO_byPlane;
  ofloat xpos1[IM_MAXDIM], xpos2[IM_MAXDIM];
  gboolean overlap;
  gchar *routine = "ObitImageMosaicFlatten";

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Must have output image defined */
  if (in->FullField==NULL) {
    Obit_log_error(err, OBIT_Error, "No flattened image defined");
    return;
  }

  /* plane number */
  planeNo = 1;
  ObitInfoListGetTest(in->info, "planeNo", &type, dim, &planeNo);

  /* Tell user */
  Obit_log_error(err, OBIT_InfoErr, "Flattening flys eye to single image plane %d", planeNo);

  /* Copy beam parameters from first image to flattened */
  ObitImageOpen(in->images[0], OBIT_IO_ReadOnly, err);
  ObitImageOpen(in->FullField, OBIT_IO_ReadWrite, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  in->FullField->myDesc->beamMaj = in->images[0]->myDesc->beamMaj;
  in->FullField->myDesc->beamMin = in->images[0]->myDesc->beamMin;
  in->FullField->myDesc->beamPA  = in->images[0]->myDesc->beamPA;
  in->FullField->myDesc->niter = 1;

  /* Create weight scratch array (temp local pointers for output) */
  sc1 = in->FullField->image;              /* Flattened image FArray */
  out = in->FullField;
  sc2 = newObitFArray("Scratch Weight array");
  ObitFArrayClone (sc1, sc2, err);         /* weight array */
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  ObitImageClose(in->images[0], err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Working, memory only Images */
  tout1 = newObitImage ("Interpolated Image");
  tout2 = newObitImage ("Weight Image");

  /* How big do we want ? */
  /*radius = MAX (in->FOV/(3600.0*in->xCells), in->FOV/(3600.0*in->yCells));*/
  radius = (olong)(0.5 + ObitImageMosaicFOV(in, err)/fabs(in->xCells));
  if (err->error) Obit_traceback_msg (err, routine, in->name);
 
  /* Zero fill accumulations */
  ObitFArrayFill (sc1, 0.0);
  ObitFArrayFill (sc2, 0.0);

  for (i=0; i<IM_MAXDIM; i++) blc[i] = 1;
  for (i=0; i<IM_MAXDIM; i++) trc[i] = 0;
  /* set plane */
  blc[2] = trc[2] = planeNo;

  /* Loop over tiles */
  for (i= 0; i<in->numberImages; i++) { /* loop 500 */

    /* Ignore if tapered */
    if (in->BeamTaper[i]>0.0) continue;

    /* Open full image */
    dim[0] = IM_MAXDIM;
    ObitInfoListPut (in->images[i]->info, "BLC", OBIT_long, dim, blc, err); 
    ObitInfoListPut (in->images[i]->info, "TRC", OBIT_long, dim, trc, err); 
    dim[0] = 1;
    ObitInfoListPut (in->images[i]->info, "IOBy", OBIT_long, dim, &IOsize, err);
    ObitImageOpen(in->images[i], OBIT_IO_ReadOnly, err);
    ObitImageClose(in->images[i], err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);

    /* Input subimage dimension  (trim edges) */
    naxis = in->images[i]->myDesc->inaxes;

    /* Fudge a bit at the edges */
    rad = radius + 3;
    if (rad > ((naxis[0]/2)-2)) rad = (naxis[0]/2) - 3;
    if (rad > ((naxis[1]/2)-2)) rad = (naxis[1]/2) - 3;
    blc[0] = (naxis[0] / 2) - rad;
    blc[0] = MAX (2, blc[0]);
    blc[1] = (naxis[1] / 2) + 1 - rad;
    blc[1] = MAX (2, blc[1]);
    trc[0] = (naxis[0] / 2) + rad;
    trc[0] = MIN (naxis[0]-1, trc[0]);
    trc[1] = (naxis[1] / 2) + 1 + rad;
    trc[1] = MIN (naxis[1]-1, trc[1]);

   /* Open/read sub window of image */
    dim[0] = IM_MAXDIM;
    ObitInfoListPut (in->images[i]->info, "BLC", OBIT_long, dim, blc, err); 
    ObitInfoListPut (in->images[i]->info, "TRC", OBIT_long, dim, trc, err); 
    dim[0] = 1;

    /* Is there some overlap with flattened image? */
    in->images[i]->extBuffer = TRUE;  /* Don't need buffer here */
    ObitImageOpen(in->images[i], OBIT_IO_ReadOnly, err);
    ObitImageClose(in->images[i], err);
    in->images[i]->extBuffer = FALSE;  /* May need buffer later */
    overlap = ObitImageDescOverlap(in->images[i]->myDesc, in->FullField->myDesc, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);

    if (overlap) {    /* Something to do? */
      /* Make or resize Interpolated scratch arrays (memory only images) */
      ObitImageClone2 (in->images[i], out, tout1, err);
      ObitImageClone2 (in->images[i], out, tout2, err);

      /* debug 
      fprintf(stderr,"Interpolate tile %d %d  %d\n",
	      i,tout1->image->naxis[0],tout1->image->naxis[1]);*/
	
      naxis = tout1->myDesc->inaxes; /* How big is output */
  
      /* Interpolate and weight image */
      /* reopen with windowing */
      ObitImageOpen (in->images[i], OBIT_IO_ReadOnly, err);
      ObitImageRead (in->images[i], NULL, err); /* Read plane */
      if (err->error) Obit_traceback_msg (err, routine, in->name);

      ObitImageUtilInterpolateWeight (in->images[i], tout1, tout2, TRUE, rad, 
				      plane, plane, hwidth, err);
      if (err->error) Obit_traceback_msg (err, routine, in->name);

      /* DEBUG 
	 ObitImageUtilArray2Image ("DbugInterp.fits", 1, tout1->image, err);
	 ObitImageUtilArray2Image ("DbugInput.fits", 1, in->images[i]->image, err);*/

      /* Close, deallocate buffer */
      ObitImageClose(in->images[i], err);
      in->images[i]->image = ObitFArrayUnref(in->images[i]->image); /* Free buffer */
      if (err->error) Obit_traceback_msg (err, routine, in->name);

      /* Paste into accumulation images */
      /* Weighted image - first need pixel alignment */
      pos1[0] = tout1->image->naxis[0]/2; /* Use center of tile */
      pos1[1] = tout1->image->naxis[1]/2;
      xpos1[0] = pos1[0];
      xpos1[1] = pos1[1];
      /* Find corresponding pixel in output */
      ObitImageDescCvtPixel (tout1->myDesc, in->FullField->myDesc, xpos1, xpos2, err);
      if (err->error) Obit_traceback_msg (err, routine, in->name);
      pos2[0] = (olong)(xpos2[0]+0.5);
      pos2[1] = (olong)(xpos2[1]+0.5);

      /* accumulate image*weight */
      ObitFArrayShiftAdd (sc1, pos2, tout1->image, pos1, 1.0, sc1);

      /* accumulate weight */
      ObitFArrayShiftAdd (sc2, pos2, tout2->image, pos1, 1.0, sc2);

    } /* end if overlap */
      /* reset window on image */
    dim[0] = IM_MAXDIM;
    for (j=0; j<IM_MAXDIM; j++) {blc[j] = 1; trc[j] = 0;}
    ObitInfoListPut (in->images[i]->info, "BLC", OBIT_long, dim, blc, err); 
    ObitInfoListPut (in->images[i]->info, "TRC", OBIT_long, dim, trc, err); 
    dim[0] = 1;
    /* Open and close to reset descriptor */
    in->images[i]->extBuffer = TRUE;  /* Don't need buffer here */
    ObitImageOpen(in->images[i], OBIT_IO_ReadOnly, err);
    ObitImageClose(in->images[i], err);
    in->images[i]->extBuffer = FALSE;  /* May need buffer later */
 } /* end loop  L500 over input images */

  /* DEBUG
  ObitImageUtilArray2Image ("DbugSumWI.fits", 1, sc1, err);
  ObitImageUtilArray2Image ("DbugSumWW.fits", 1, sc2, err); */

  /* Normalize */
  ObitFArrayDivClip (sc1, sc2, 0.01, out->image);

  /* Write output */
  ObitImageWrite (in->FullField, out->image->array, err); 
  ObitImageClose (in->FullField, err);
  in->FullField->image = ObitFArrayUnref(in->FullField->image); /* Free buffer */

  /* Clear BLC,TRC on images */
  dim[0] = IM_MAXDIM;
  for (i=0; i<IM_MAXDIM; i++) blc[i] = 1;
  for (i=0; i<IM_MAXDIM; i++) trc[i] = 0;
  ObitInfoListPut (in->FullField->info, "BLC", OBIT_long, dim, blc, err); 
  ObitInfoListPut (in->FullField->info, "TRC", OBIT_long, dim, trc, err); 
  for (i=0; i<in->numberImages; i++) {
    ObitInfoListPut (in->images[i]->info, "BLC", OBIT_long, dim, blc, err); 
    ObitInfoListPut (in->images[i]->info, "TRC", OBIT_long, dim, trc, err); 
  }

  /* Cleanup */
  sc2   = ObitFArrayUnref(sc2);
  tout1 = ObitImageUnref(tout1);
  tout2 = ObitImageUnref(tout2);

} /* end ObitImageMosaicFlatten */

/**
 * Return the Field of View of a Mosaic
 * Given value on object unless it is zero and then it works it out from the 
 * array of images and saves it on the object.
 * \param in      The object with images
 * \param err     Error stack, returns if not empty.
 * \return the full width of the field of view in deg.
 */
ofloat ObitImageMosaicFOV (ObitImageMosaic *in, ObitErr *err)
{
  ofloat FOV = 0.0;
  ofloat pixel[2], xshift, yshift, shift, rotate=0.0;
  odouble pos[2], pos0[2];
  olong i, n;
  ObitImageDesc *desc=NULL;
  gchar *routine = "ObitImageMosaicFOV";

  /* error checks */
  if (err->error) return FOV;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Do we have it already? */
  if (in->FOV>0.0) return in->FOV;

  /* No images? We're wasting time */
  if (in->numberImages<=0) return FOV;

  /* If only one, it's the size */
  if (in->numberImages==1) {
    desc = in->images[0]->myDesc;
    in->FOV = MAX ((desc->inaxes[0]*fabs(desc->cdelt[0])), 
		   (desc->inaxes[1]*fabs(desc->cdelt[1])));
    return in->FOV;
  }

  /* Use the center of the first image as the reference position */
  desc = in->images[0]->myDesc;
  FOV = 0.5 * MAX ((desc->inaxes[0]*fabs(desc->cdelt[0])), 
		   (desc->inaxes[1]*fabs(desc->cdelt[1])));
  ObitImageDescGetPos (desc, desc->crpix, pos0, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, FOV);
  
  /* Loop over images, in Fly's eye unless that size is zero, then all */
  n = in->nFlyEye;
  if (n<=0) n = in->numberImages;
  for (i=0; i<n; i++) {
    desc = in->images[i]->myDesc;
    /* All four corners */
    pixel[0] = 1; pixel[1] = 1;
    ObitImageDescGetPos (desc, pixel, pos, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, FOV);
    /* Get separation */
    ObitSkyGeomShiftXY (pos0[0], pos0[1], rotate, pos[0], pos[1], 
			&xshift, &yshift);
    shift = sqrt (xshift*xshift + yshift*yshift);
    FOV = MAX (FOV, shift);

    pixel[0] = desc->inaxes[0]; pixel[1] = desc->inaxes[1];
    ObitImageDescGetPos (desc, pixel, pos, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, FOV);
    /* Get separation */
    ObitSkyGeomShiftXY (pos0[0], pos0[1], rotate, pos[0], pos[1], 
			&xshift, &yshift);
    shift = sqrt (xshift*xshift + yshift*yshift);
    FOV = MAX (FOV, shift);

    pixel[0] = desc->inaxes[0]; pixel[1] = 1;
    ObitImageDescGetPos (desc, pixel, pos, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, FOV);
    /* Get separation */
    ObitSkyGeomShiftXY (pos0[0], pos0[1], rotate, pos[0], pos[1], 
			&xshift, &yshift);
    shift = sqrt (xshift*xshift + yshift*yshift);
    FOV = MAX (FOV, shift);

    pixel[0] = 1; pixel[1] = desc->inaxes[1];
    ObitImageDescGetPos (desc, pixel, pos, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, FOV);
    /* Get separation */
    ObitSkyGeomShiftXY (pos0[0], pos0[1], rotate, pos[0], pos[1], 
			&xshift, &yshift);
    shift = sqrt (xshift*xshift + yshift*yshift);
    FOV = MAX (FOV, shift);
  } /* end loop over images */

  in->FOV = FOV;  /* Save result */

  return in->FOV;
} /* end ObitImageMosaicFOV */

/**
 * Return the maximum Field of View of a Mosaic
 * Searches list of images for which startCC[i]>endCC[i] and returns
 * the maximum distance of any corner from the pointing center.
 * Default is the mosaic FOV (if given)
 * \param in      The object with images
 * \param startCC Array of starting component numbers
 * \param endCC   Array of ending component numbers
 * \param err     Error stack, returns if not empty.
 * \return the full width of the field of view in deg.
 */
ofloat ObitImageMosaicMaxFOV (ObitImageMosaic *in, olong *startCC, olong *endCC, 
			      ObitErr *err)
{
  ofloat FOV = -1.0;
  ofloat pixel[2], xshift, yshift, shift, rotate=0.0;
  odouble pos[2], pos0[2];
  olong i, n;
  ObitImageDesc *desc=NULL;
  gchar *routine = "ObitImageMosaicFOV";

  /* error checks */
  if (err->error) return 0.0;
  g_assert (ObitIsA(in, &myClassInfo));

  /* No images? We're wasting time */
  if (in->numberImages<=0) return 0.0;

  /* Use the pointing of the first image as the reference position */
  desc = in->images[0]->myDesc;
  ObitImageDescGetPoint (desc, &pos0[0], &pos0[1]); 
  
  /* Loop over images */
  n = in->numberImages;
  for (i=0; i<n; i++) {
    /* Pay attention to this one? (components in model?) */
    if (startCC[i]>endCC[i]) continue;

    desc = in->images[i]->myDesc;
    /* All four corners */
    pixel[0] = 1; pixel[1] = 1;
    ObitImageDescGetPos (desc, pixel, pos, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, FOV);
    /* Get separation */
    ObitSkyGeomShiftXY (pos0[0], pos0[1], rotate, pos[0], pos[1], 
			&xshift, &yshift);
    shift = sqrt (xshift*xshift + yshift*yshift);
    FOV = MAX (FOV, shift);

    pixel[0] = desc->inaxes[0]; pixel[1] = desc->inaxes[1];
    ObitImageDescGetPos (desc, pixel, pos, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, FOV);
    /* Get separation */
    ObitSkyGeomShiftXY (pos0[0], pos0[1], rotate, pos[0], pos[1], 
			&xshift, &yshift);
    shift = sqrt (xshift*xshift + yshift*yshift);
    FOV = MAX (FOV, shift);

    pixel[0] = desc->inaxes[0]; pixel[1] = 1;
    ObitImageDescGetPos (desc, pixel, pos, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, FOV);
    /* Get separation */
    ObitSkyGeomShiftXY (pos0[0], pos0[1], rotate, pos[0], pos[1], 
			&xshift, &yshift);
    shift = sqrt (xshift*xshift + yshift*yshift);
    FOV = MAX (FOV, shift);

    pixel[0] = 1; pixel[1] = desc->inaxes[1];
    ObitImageDescGetPos (desc, pixel, pos, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, FOV);
    /* Get separation */
    ObitSkyGeomShiftXY (pos0[0], pos0[1], rotate, pos[0], pos[1], 
			&xshift, &yshift);
    shift = sqrt (xshift*xshift + yshift*yshift);
    FOV = MAX (FOV, shift);
  } /* end loop over images */

  /* Got something */
  if (FOV>=0.0) return FOV;

  /* If previous has failed and FOV given, return it */
  if (in->FOV>0.0) return in->FOV;


  return in->FOV;
} /* end ObitImageMosaicMaxFOV */

/**
 * Finds peak summed component within a given radius.  
 * Also determines offset of peak from the central by a moment analysis.  
 * Routine translated from the AIPSish VLAUTIL.FOR/VLMXCC  
 * \param CCTab     CC table object to filter. 
 * \param nccpos    Number of CC to consider. 
 * \param radius    Radius within which to consider components. (deg) 
 * \param maxcmp    [out] Maximum summed CC flux in search radius. 
 * \param xcenter   [out] delta X (ra) of brightest pixel (deg)
 * \param ycenter   [out] delta Y (dec) of brightest pixel (deg)
 * \param xoff      [out] Offset from central pixel in x (ra) (deg) 
 * \param yoff      [out] Offset from central pixel in y (dec) (deg) 
 * \param err       Error/message stack
 */
void ObitImageMosaicMaxCC (ObitTableCC *CCTab, olong nccpos, ofloat radius, 
			   ofloat* maxcmp, ofloat* xcenter, ofloat* ycenter, 
			   ofloat* xoff, ofloat* yoff, ObitErr* err)  
{
  ObitTableCCRow *CCRow = NULL;
  olong   irow, nrow, lim, ncc, i, j, peakcc;
  ofloat mxdis, sum, xsum, ysum, dis, xcen, ycen;
  ofloat *xpos=NULL, *ypos=NULL, *flux=NULL;
  gchar *routine = "ObitImageMosaicMaxCC";

  /* Error checks */
  if (err->error) return;  /* previous error? */

  /* Default values */
  *maxcmp = 0.0;
  *xoff   = 0.0;
  *yoff   = 0.0;

  /* Anything to do? */
  if (nccpos <= 0) return;

  /* Open table */
  ObitTableCCOpen (CCTab, OBIT_IO_ReadOnly, err);
  if (err->error) goto cleanup;

  /* allocate storage */
  nrow = CCTab->myDesc->nrow;
  if (nrow<1) return;
  xpos = g_malloc0(nrow*sizeof(ofloat));
  ypos = g_malloc0(nrow*sizeof(ofloat));
  flux = g_malloc0(nrow*sizeof(ofloat));

  /* Copy table to arrays */
  CCRow = newObitTableCCRow (CCTab);
  ncc = 0;
  nrow = CCTab->myDesc->nrow;
  for (irow=1; irow<=nrow; irow++) {
    ObitTableCCReadRow (CCTab, irow, CCRow, err);
    if (err->error) goto cleanup;
    if (CCRow->status==-1) continue;  /* Deselected? */
    xpos[ncc]   = CCRow->DeltaX;
    ypos[ncc]   = CCRow->DeltaY;
    flux[ncc++] = CCRow->Flux;
  } /* end loop reading table */

  /* Sum up to nccpos fluxes within radius */
  mxdis = radius * radius;
  *maxcmp = -1.0e20;
  peakcc = 0;
  lim = MIN (MAX (1, nccpos), ncc);
  for (i=0; i<lim; i++) { /* loop 400 */
    sum = 0.0;
    for (j=0; j<lim; j++) { /* loop 300 */
      dis = (xpos[i]-xpos[j])*(xpos[i]-xpos[j]) + (ypos[i]-ypos[j])*(ypos[i]-ypos[j]);
      if (dis <= mxdis) sum += flux[j];
    } /* end loop  L300: */
    if ((sum > *maxcmp) && (flux[i]>=0.2*(*maxcmp))) {
      peakcc = i;
      *maxcmp = sum;
    } 
  } /* end loop  L400: */

  /* Get position offset, use all components. */
  xcen = xpos[peakcc];
  ycen = ypos[peakcc];
  sum  = 0.0;
  xsum = 0.0;
  ysum = 0.0;
  nrow = CCTab->myDesc->nrow;
  
  /* Get first moments about reference pixel (0,0). */
  for (i=0; i<ncc; i++) { /* loop 500 */
    dis = (xpos[i]-xcen)*(xpos[i]-xcen) + (ypos[i]-ycen)*(ypos[i]-ycen);
    if (dis <= mxdis) {
      sum  += flux[i];
      xsum += (xpos[i]-xcen) * flux[i];
      ysum += (ypos[i]-ycen) * flux[i];
    } 
  } /* end loop  L500: */
  if (fabs (sum)  <= 1.0e-8) sum = 1.0;

  /* Get moments */
  *xoff = xsum / sum;
  *yoff = ysum / sum;

  /* Return peak pixel */
  *xcenter = xcen;
  *ycenter = ycen;

  /* Don't accept crazy values */
  if (fabs (*xoff)  >  radius) *xoff = 0.0;
  if (fabs (*yoff)  >  radius) *yoff = 0.0;
  
  /* Close table */
  ObitTableCCClose (CCTab, err);
  if (err->error) goto cleanup;
      
 cleanup:
   /* deallocate storage */
  CCRow = ObitTableCCRowUnref(CCRow);  
  if (xpos) g_free(xpos); 
  if (ypos) g_free(ypos); 
  if (flux) g_free(flux); 
  if (err->error) Obit_traceback_msg (err, routine, CCTab->name);
  
} /* end of routine ObitImageMosaicMaxCC */ 

/**
 * Creates a CC table with the components from this field and 
 * overlapping ones from adjacent fields. 
 * If field is an autoCenter field, only add the corresponding shifted field
 * \param mosaic    Input image mosaic
 * \param field     The field (1-rel) of interest
 * \param CCVer     CC table version of input
 * \param err       Error/message stack
 * \return combined CC table
 */
ObitTableCC* ObitImageMosaicCombineCC (ObitImageMosaic *mosaic, olong field, 
				       olong CCVer, ObitErr* err)  
{
  ObitTableCC *inCC=NULL, *outCC=NULL;
  ObitTableCCRow *CCRow = NULL;
  ObitImageDesc *imDesc1=NULL, *imDesc2=NULL;
  olong   irow, orow, nrow, noParms, ifield, jfield, nx, ny, ver, hiVer;
  ofloat pixel1[2], pixel2[2];
  gboolean overlap;
  gchar *routine = "ObitImageMosaicCombineCC";

  /* Error checks */
  if (err->error) return outCC;  /* previous error? */
  g_assert (ObitIsA(mosaic, &myClassInfo));

  /* Merge CC table on field to new output on same image */
  ver = CCVer;
  noParms = 0;
  ifield = field-1;  /* 0 rel */
  inCC = newObitTableCCValue ("input CC", (ObitData*)mosaic->images[ifield],
				   &ver, OBIT_IO_ReadOnly, noParms, 
				   err);
  /* New output version */
  hiVer = ObitTableListGetHigh (mosaic->images[ifield]->tableList, "AIPS CC");
  ver = hiVer+1;
  outCC = newObitTableCCValue ("Comb CC", (ObitData*)mosaic->images[ifield],
				    &ver, OBIT_IO_ReadWrite, inCC->noParms, 
				    err);
  ObitTableCCUtilMerge (inCC, outCC, err);
  if  (err->error) Obit_traceback_val (err, routine, mosaic->images[ifield]->name, outCC);
  inCC  = ObitTableCCUnref(inCC);

  /* If there is nothing in the table ignore */
  if ((outCC->myDesc->nrow<=0) || (ver==CCVer)) return outCC;

  CCRow = newObitTableCCRow (outCC);

  /* Open output table */
  ObitTableCCOpen (outCC, OBIT_IO_ReadWrite, err);
  if  (err->error) Obit_traceback_val (err, routine, mosaic->images[ifield]->name, outCC);

  /* Loop over other fields */
  imDesc2 = (mosaic->images[ifield])->myDesc;
  nx = imDesc2->inaxes[0];  ny = imDesc2->inaxes[1]; 
  for (jfield = 0; jfield<mosaic->numberImages; jfield++) {
    if (jfield==ifield) continue;
    imDesc1 = (mosaic->images[jfield])->myDesc;
 
    /* If 2D and autoCen ignore all but shifted version */
    if ((!mosaic->images[ifield]->myDesc->do3D) && 
	((mosaic->isAuto[ifield]>0) && (mosaic->isShift[jfield]!=field))) continue;
      
   /* Any overlap? */
    if (ObitImageDescOverlap(imDesc1, imDesc2, err)) {
      CCVer   = 1;
      noParms = 0;
      inCC = newObitTableCCValue ("input CC", (ObitData*)mosaic->images[jfield],
				  &CCVer, OBIT_IO_ReadOnly, noParms, 
				  err);
      /* If no CC table just ignore it */
      if (inCC==NULL)  continue;

      /* Open table */
      ObitTableCCOpen (inCC, OBIT_IO_ReadOnly, err);
      if  (err->error) Obit_traceback_val (err, routine, mosaic->images[ifield]->name, outCC);

      /* Loop over table copying any overlapping CCs */
      nrow = inCC->myDesc->nrow;
      for (irow=1; irow<=nrow; irow++) {
	ObitTableCCReadRow (inCC, irow, CCRow, err);
	if  (err->error) Obit_traceback_val (err, routine, mosaic->images[ifield]->name, outCC);

	/* where in output image? */
	pixel1[0] = imDesc1->crpix[0] + CCRow->DeltaX/imDesc1->cdelt[0];
	pixel1[1] = imDesc1->crpix[1] + CCRow->DeltaY/imDesc1->cdelt[1];
	overlap = ObitImageDescCvtPixel(imDesc1, imDesc2, pixel1, pixel2, err);
	if  (err->error) Obit_traceback_val (err, routine, mosaic->images[jfield]->name, outCC);

	/* Does this one overlap? */
	if (overlap) {

	  /* Update position if 3D (2D is OK) */
	  if (imDesc1->do3D) {
	    CCRow->DeltaX = (pixel2[0]-imDesc1->crpix[0]) * imDesc1->cdelt[0];
	    CCRow->DeltaY = (pixel2[1]-imDesc1->crpix[1]) * imDesc1->cdelt[1];
	  }
	  
	  /* Write row */
	  orow = -1;
	  ObitTableCCWriteRow (outCC, orow, CCRow, err);
	  if  (err->error) Obit_traceback_val (err, routine, mosaic->images[ifield]->name, outCC);
	}
      } /* end loop over table */

      /* Closeup/Cleanup */
      ObitTableCCClose (inCC, err);
      inCC  = ObitTableCCUnref(inCC);
    } /* end if overlap */
  } /* end loop over other fields */

  /* Update number of CCs */
  ObitImageOpen(mosaic->images[ifield], OBIT_IO_ReadWrite, err);
  mosaic->images[ifield]->myDesc->niter = outCC->myDesc->nrow;
  mosaic->images[ifield]->myStatus = OBIT_Modified;
  ObitImageClose(mosaic->images[ifield], err);
  
  /* Close output table */
  ObitTableCCClose (outCC, err);
  CCRow = ObitTableCCRowUnref(CCRow);  

  return outCC;
} /* end of routine ObitImageMosaicCombineCC */ 

/**
 * Zeros all CC fluxes within radius of xcenter, ycenter 
 * \param CCTab     CC table object to filter. 
 * \param nccpos    Number of CC to consider. 
 * \param radius    Radius within which to consider components. (deg) 
 * \param xcenter   delta X (ra) of brightest pixel (deg)
 * \param ycenter   delta Y (dec) of brightest pixel (deg)
 * \param err       Error/message stack
 */
void ObitImageMosaicFlagCC (ObitTableCC *CCTab, olong nccpos, ofloat radius, 
			   ofloat xcenter, ofloat ycenter, ObitErr* err) 
{
  ObitTableCCRow *CCRow = NULL;
  olong   irow, nrow;
  ofloat mxdis, dis;
  gchar *routine = "ObitImageMosaicFlagCC";

  /* Error checks */
  if (err->error) return;  /* previous error? */

  /* Anything to do? */
  if (nccpos <= 0) return;
  mxdis = radius * radius;

  /* Open table */
  ObitTableCCOpen (CCTab, OBIT_IO_ReadWrite, err);
  if (err->error) goto cleanup;

  /* allocate storage */
  nrow = CCTab->myDesc->nrow;
  if (nrow<1) return;

  /* Copy table to arrays */
  CCRow = newObitTableCCRow (CCTab);
  nrow = CCTab->myDesc->nrow;
  for (irow=1; irow<=nrow; irow++) {
    ObitTableCCReadRow (CCTab, irow, CCRow, err);
    if (err->error) goto cleanup;
    if (CCRow->status==-1) continue;  /* Deselected? */
    /* Zero this one? */
    dis = (CCRow->DeltaX-xcenter)*(CCRow->DeltaX-xcenter) + 
      (CCRow->DeltaY-ycenter)*(CCRow->DeltaY-ycenter);
    if (dis <= mxdis) { /* Rewrite? */
      CCRow->Flux = 0.0;
      ObitTableCCWriteRow (CCTab, irow, CCRow, err);
      if (err->error) goto cleanup;
    } 
  } /* end loop checking table */

 cleanup:
   /* deallocate storage */
  CCRow = ObitTableCCRowUnref(CCRow);  
  /* Close table */
  ObitTableCCClose (CCTab, err);
  if (err->error) Obit_traceback_msg (err, routine, CCTab->name);
  
} /* end of routine ObitImageMosaicFlagCC */ 

/**
 * Make a single field ImageMosaic corresponding to the field with
 * the highest summed peak of CCs if any above MinFlux.
 * Includes CCs from overlapping fields both in the test and in  the
 * CC table on the imaga on the output mosaic.
 * Only considers autoCenter Fields, if any.
 * Fields with 1-rel numbers in the zero terminated list ignore are ignored
 * \param mosaic   ImageMosaic to process
 * \param MinFlux  Min. flux density for operation.
 * \param ignore   0 terminated list of 1-rel field numbers to ignore
 * \param field    [out] the 1-rel number of the field copied
 * \param err      Error/message stack
 * return Newly created ImageMosaic or NULL
 */
ObitImageMosaic* ObitImageMosaicMaxField (ObitImageMosaic *mosaic, 
					  ofloat MinFlux, olong *ignore, olong *field,
					  ObitErr* err) 
{
  ObitImageMosaic* out=NULL;
  ObitTableCC *CCTab=NULL, *inCCTable=NULL, *outCCTable=NULL, *saveCCTable=NULL;
  ObitImage   *tmpImage=NULL, *tmpBeam=NULL;
  gint32 dim[MAXINFOELEMDIM];
  ObitInfoType type;
  olong   i, nfield, ifield, itemp, noParms, nccpos;
  olong   maxField, other;
  olong  CCVer;
  ofloat tmax, xcenter, ycenter, xoff, yoff, radius, cells[2], maxCC;
  gboolean forget, doAuto;
  gchar *routine = "ObitImageMosaicMaxField";

  /* Error checks */
  if (err->error) return out;  /* previous error? */
  g_assert(ObitImageMosaicIsA(mosaic));

  /* Number of fields */
  nfield = mosaic->numberImages;

  /* Any autoCenter fields? */
  doAuto = FALSE;
  for (ifield=0; ifield<nfield; ifield++) { 
   if (mosaic->isAuto[ifield]>=0) doAuto = TRUE;
 }

  /* Get cellsize */
  cells[0] =  fabs(mosaic->xCells); cells[1] = fabs(mosaic->yCells);

  /* Consider components within 2.5  cells  */
  radius = 2.5 * cells[0];

  /* CC table(s) */
  itemp = 1;
  ObitInfoListGetTest(mosaic->info, "CCVer", &type, dim, &itemp);
  CCVer = itemp;

  /* Loop over fields */
  maxField = -1;
  maxCC = -1.0e20;
  for (ifield=0; ifield<nfield; ifield++) { /* loop 500 */

    /* Is this field in the ignore list? */
    forget = FALSE;
    i = 0;
    while (ignore[i]>0) {
      forget = forget || (ignore[i++]==(ifield+1));
      if (forget) break;
    }
    if (forget) continue;

    /* Open image in case header needs update */
    ObitImageOpen (mosaic->images[ifield], OBIT_IO_ReadWrite, err);
    if  (err->error) Obit_traceback_val (err, routine, mosaic->images[ifield]->name, out);

    /* Only consider autoCenter Fields? */
    if (doAuto && (mosaic->isAuto[ifield]<0)) goto closem;

    /* Ignore if this is a 2D shifted field  */
    if ((!mosaic->images[ifield]->myDesc->do3D) && (mosaic->isShift[ifield]>=0)) goto closem;

    /* If autoCenter use only the CCs in this image - make copy */
    if (doAuto) {
      CCVer = 1;
      noParms = 0;
      inCCTable = newObitTableCCValue ("Peeled CC", (ObitData*)mosaic->images[ifield],
				       &CCVer, OBIT_IO_ReadOnly, noParms, 
				       err);
      CCVer   = 0;
      CCTab = newObitTableCCValue ("CCTab", (ObitData*)mosaic->images[ifield],
				   &CCVer, OBIT_IO_WriteOnly, noParms, err);
      CCTab     = ObitTableCCCopy (inCCTable, CCTab, err);
      inCCTable = ObitTableCCUnref(inCCTable);

      /* Add any components from the 2D shifted version */
      if ((!mosaic->images[ifield]->myDesc->do3D) && (mosaic->isAuto[ifield]>0)) {
	other = mosaic->isAuto[ifield];
	CCVer = 1;
	noParms = 0;
	inCCTable = newObitTableCCValue ("Shift CC", (ObitData*)mosaic->images[other-1],
					 &CCVer, OBIT_IO_ReadOnly, noParms, 
					 err);
	ObitTableCCUtilAppend (inCCTable, CCTab, 0, 0, err);
	if  (err->error) Obit_traceback_val (err, routine, mosaic->images[other-1]->name, out);
	inCCTable   = ObitTableCCUnref(inCCTable);
     }
    } else {
      /* Not autoCenter - Make CC table object including CCs from any overlapping facets */
      CCTab = ObitImageMosaicCombineCC (mosaic, ifield+1, CCVer, err);
    }
      if  (err->error) Obit_traceback_val (err, routine, mosaic->images[ifield]->name, out);

    /* Determine maximum */
    nccpos = CCTab->myDesc->nrow;
    ObitImageMosaicMaxCC (CCTab, nccpos, radius, &tmax, &xcenter, &ycenter, &xoff, &yoff, err);
    if  (err->error) Obit_traceback_val (err, routine, mosaic->images[ifield]->name, out);

    /* Tell about it */
    if (err->prtLv>=2)
      Obit_log_error(err, OBIT_InfoErr, "Peel Max. field %d is %f", ifield+1,tmax);
    
    /* this one of interest? */
    if ((tmax > MinFlux) && (tmax >= maxCC)) {
      /* Zap a previous, now superceeded file */
      if (saveCCTable && (maxField>=0)) {
	ObitDataZapTable((ObitData*)mosaic->images[maxField], saveCCTable->tabType, 
			 saveCCTable->tabVer, err);
	if  (err->error) Obit_traceback_val (err, routine, mosaic->images[ifield]->name, out);
	saveCCTable = ObitTableCCUnref(saveCCTable);
      }
      maxCC = tmax;
      maxField = ifield;
      saveCCTable = ObitTableCCRef(CCTab);  /* Save table */
    } else {
      
      /* Delete temporary table */
      ObitDataZapTable((ObitData*)mosaic->images[ifield], CCTab->tabType, CCTab->tabVer, err);
      if  (err->error) Obit_traceback_val (err, routine, mosaic->images[ifield]->name, out);
      CCTab = ObitTableCCUnref(CCTab);
    }

    /* Close/update image */
  closem:
    ObitImageClose(mosaic->images[ifield], err);
    if  (err->error) Obit_traceback_val (err, routine, mosaic->images[ifield]->name, out);
  } /* end loop  L500: */

  /* Catch anything? */
  if (maxField<0) return out;

  /* Make output ImageMosaic */
  out = newObitImageMosaic ("Temp mosaic", 1);
  tmpImage = newObitImageScratch (mosaic->images[maxField], err);
  /* Copy image */
  tmpImage = ObitImageCopy(mosaic->images[maxField], tmpImage, err);
  ObitImageMosaicSetImage (out, 0, tmpImage, err);
  if  (err->error) Obit_traceback_val (err, routine, mosaic->images[maxField]->name, out);
  /* Give more sensible name */
  if (tmpImage->name) g_free(tmpImage->name);
  tmpImage->name = g_strdup("Peel Image");

  /* Copy beam */
  tmpBeam = newObitImageScratch ((ObitImage*)mosaic->images[maxField]->myBeam, err);
  tmpBeam = ObitImageCopy((ObitImage*)mosaic->images[maxField]->myBeam, tmpBeam, err);
  ObitImageMosaicSetImage (out, 0, tmpImage, err);
  if  (err->error) Obit_traceback_val (err, routine, mosaic->images[maxField]->name, out);
  tmpImage->myBeam = ObitImageUnref(tmpImage->myBeam);
  tmpImage->myBeam = ObitImageRef(tmpBeam);
  /* Give more sensible name */
  if (tmpBeam->name) g_free(tmpBeam->name);
  tmpBeam->name = g_strdup("Peel Beam");
  tmpImage = ObitImageUnref(tmpImage);
  tmpBeam  = ObitImageUnref(tmpBeam);

  /* Other info */
  out->xCells   = mosaic->xCells;
  out->yCells   = mosaic->yCells;
  out->fileType = mosaic->fileType;
  out->bmaj     = mosaic->bmaj;
  out->bmin     = mosaic->bmin;
  out->bpa      = mosaic->bpa;
  out->nx[0]       = mosaic->nx[maxField];
  out->ny[0]       = mosaic->ny[maxField];
  out->nplane[0]   = mosaic->nplane[maxField];
  out->RAShift[0]  = mosaic->RAShift[maxField];
  out->DecShift[0] = mosaic->DecShift[maxField];

  /* Copy Clean components from saveCCTable */
  CCVer = 1;
  noParms = 0;
  if (saveCCTable) noParms = saveCCTable->noParms;
  outCCTable = newObitTableCCValue ("outCC", (ObitData*)out->images[0],
				    &CCVer, OBIT_IO_ReadWrite, noParms, err);
  outCCTable  = ObitTableCCCopy (saveCCTable, outCCTable, err);
  if  (err->error) Obit_traceback_val (err, routine, mosaic->images[ifield]->name, out);
  outCCTable  = ObitTableCCUnref(outCCTable);
  /* Delete saveCCTable */
  if (saveCCTable && (maxField>=0)) {
    ObitDataZapTable((ObitData*)mosaic->images[maxField], saveCCTable->tabType, 
		     saveCCTable->tabVer, err);
    if  (err->error) Obit_traceback_val (err, routine, mosaic->images[ifield]->name, out);
    saveCCTable = ObitTableCCUnref(saveCCTable);
  }

  *field = maxField+1;  /* Which field is this? */
  return out;
} /* end of routine ObitImageMosaicMaxField */ 

/**
 * Convert structure information to entries in an ObitInfoList
 * \param in      Object of interest.
 * \param prefix  If NonNull, string to be added to beginning of outList entry name
 *                "xxx" in the following
 * \param outList InfoList to write entries into
 *      \li "xxxnumberImages" olong Number of images in Mosaic
 *      \li "xxxnFlyEye"      olong Number of images already initialized
 *      \li "xxxImagennnnn"   string Prefix for each image nnnnn (1-rel) in images
 *      \li "xxxFullField"    string Prefix for Full field image
 *      \li "xxxFOV"          ofloat Field of view as radius (deg)
 *      \li "xxxRadius"       ofloat Radius of the usable region of a given tile (cells)
 *      \li "xxxxCells"       ofloat Cell Spacing in X (deg)
 *      \li "xxxyCells"       ofloat Cell Spacing in Y (deg)
 *      \li "xxxOutlierSize"  olong requested size (CLEAN window) of outliers
 *      \li "xxxnx"           olong* Number of pixels in X for each image
 *      \li "xxxny"           olong* Number of pixels in Y for each image
 *      \li "xxxnplane"       olong* Number of planes for each image
 *      \li "xxxinFlysEye"    gboolean* Is each facet in Fly's Eye? 
 *      \li "xxxFacetNo"      olong* Untapered facet number - associate various taperings
 *      \li "xxxRAShift"      ofloat* RA shift (deg) for each image
 *      \li "xxxDecShift"     ofloat* Dec shift (deg) for each image
 *      \li "xxxfileType"     olong Are image OBIT_IO_AIPS or OBIT_IO_FITS?
 *      \li "xxximName"       string Name of Mosaic images
 *      \li "xxximClass"      string Class of Mosaic images
 *      \li "xxximSeq"        olong Sequence number of Mosaic images
 *      \li "xxximDisk"       olong* Disk number of Mosaic images
 *      \li "xxxdoFull"       boolean Is a full field image desired?
 *      \li "xxxbmaj"         ofloat Restoring beam major axis in deg.
 *      \li "xxxbmin"         ofloat Restoring beam minor axis in deg.
 *      \li "xxxbpa"          ofloat Restoring beam PA in deg.
 *      \li "xxxisAuto"       olong Is this an autoCenter field?
 *      \li "xxxisShift"      olong Is this an autoCenter shifted field?   
 *      \li "xxxnumBeamTapes" olong  Number of Beam tapers (elements in BeamTapes) 
 *      \li "xxxBeamTapes"    ofloat*  List of Beam tapers  
 *      \li "xxxBeamTaper"    ofloat*  Beam Taper per image as FWHM in deg      
 * \param err     ObitErr for reporting errors.
 */
void ObitImageMosaicGetInfo (ObitImageMosaic *in, gchar *prefix, ObitInfoList *outList, 
			     ObitErr *err)
{ 
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *keyword=NULL, *None = "None", *OK="OK";
  olong i, otemp;
  /*gchar *routine = "ObitImageMosaicGetInfo";*/

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

   /* Copy any InfoList keywords */
  if (prefix) keyword = g_strconcat (prefix, "Info", NULL);
  else       keyword = g_strdup("Info");
  ObitInfoListCopyAddPrefix (in->info, outList, keyword);

  /* "xxxnumberImages" olong Number of images in Mosaic; */
  if (prefix) keyword = g_strconcat (prefix, "numberImages", NULL);
  else        keyword = g_strdup("numberImages");
  dim[0] = 1;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_long, dim, &in->numberImages);
  g_free(keyword);

  /* "xxxnFlyEye"      olong Number of images already initialized */
  if (prefix) keyword = g_strconcat (prefix, "nFlyEye", NULL);
  else        keyword = g_strdup("nFlyEye");
  dim[0] = 1;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_long, dim, &in->nFlyEye);
  g_free(keyword);

  /* "xxxImagennnnn"  string Prefix for each image nnnnn (1-rel) in images */
  for (i=0; i<in->numberImages; i++) {
    if (in->images[i]) {
      /* Prefix for each entry in images */
      if (prefix) sprintf(keyword, "%sImage%5.5d", prefix, i+1);
      else        sprintf(keyword, "Image%5.5d", i+1);
      ObitDataGetFileInfo((ObitData*)in->images[i], keyword, outList, err);
      dim[0] = strlen(OK);
      ObitInfoListAlwaysPut(outList, keyword, OBIT_string, dim, OK);
    } else {
      dim[0] = strlen(None);
      ObitInfoListAlwaysPut(outList, keyword, OBIT_string, dim, None);
    }
    g_free(keyword);
  }

  /* "xxxFullField"    string prefix for Full field image */
  if (prefix) keyword = g_strconcat (prefix, "FullField", NULL);
  else        keyword = g_strdup("FullField");
  if (in->FullField) {
    ObitDataGetFileInfo((ObitData*)in->FullField, keyword, outList, err);
    dim[0] = strlen(OK);
    ObitInfoListAlwaysPut(outList, keyword, OBIT_string, dim, OK);
  } else {
    dim[0] = strlen(None);
    ObitInfoListAlwaysPut(outList, keyword, OBIT_string, dim, None);
  }
  g_free(keyword);

  /* "xxxFOV"          ofloat Field of view as radius (deg) */
  if (prefix) keyword = g_strconcat (prefix, "FOV", NULL);
  else        keyword = g_strdup("FOV");
  dim[0] = 1; dim[1] = 1;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_float, dim, &in->FOV);
  g_free(keyword);

  /* "xxxRadius"       ofloat Radius of the usable region of a given tile (cells) */
  if (prefix) keyword = g_strconcat (prefix, "Radius", NULL);
  else        keyword = g_strdup("Radius");
  dim[0] = 1;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_float, dim, &in->Radius);
  g_free(keyword);

  /* "xxxxCells"       ofloat Cell Spacing in X (deg) */
  if (prefix) keyword = g_strconcat (prefix, "xCells", NULL);
  else        keyword = g_strdup("xCells");
  dim[0] = 1;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_float, dim, &in->xCells);
  g_free(keyword);

  /* "xxxyCells"       ofloat Cell Spacing in Y (deg) */
  if (prefix) keyword = g_strconcat (prefix, "yCells", NULL);
  else        keyword = g_strdup("yCells");
  dim[0] = 1;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_float, dim, &in->yCells);
  g_free(keyword);

  /* "xxxOutlierSize"  olong requested size (CLEAN window) of outliers */
  if (prefix) keyword = g_strconcat (prefix, "OutlierSize", NULL);
  else        keyword = g_strdup("OutlierSize");
  dim[0] = 1;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_long, dim, &in->OutlierSize);
  g_free(keyword);

  /* "xxxnx"           olong* Number of pixels in X for each image; */
  if (prefix) keyword = g_strconcat (prefix, "nx", NULL);
  else        keyword = g_strdup("nx");
  dim[0] = in->numberImages; dim[1] = 1;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_long, dim, in->nx);
  g_free(keyword);

  /* "xxxny"           olong* Number of pixels in Y for each image */
  if (prefix) keyword = g_strconcat (prefix, "ny", NULL);
  else        keyword = g_strdup("ny");
  dim[0] = in->numberImages; dim[1] = 1;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_long, dim, &in->ny);
  g_free(keyword);

  /* "xxxnplane"       olong* Number of planes for each image */
  if (prefix) keyword = g_strconcat (prefix, "nplane", NULL);
  else        keyword = g_strdup("nplane");
  dim[0] = in->numberImages; dim[1] = 1;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_long, dim, &in->nplane);
  g_free(keyword);

  /* "xxxinFlysEye"      gboolean* Is each facet in Fly's Eye? */
  if (prefix) keyword = g_strconcat (prefix, "inFlysEye", NULL);
  else        keyword = g_strdup("inFlysEye");
  dim[0] = in->numberImages; dim[1] = 1;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_bool, dim, &in->inFlysEye);
  g_free(keyword);

  /* xxxFacetNo"      olong* Untapered facet number - associate various taperings */
  if (prefix) keyword = g_strconcat (prefix, "FacetNo", NULL);
  else        keyword = g_strdup("FacetNo");
  dim[0] = in->numberImages; dim[1] = 1;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_long, dim, &in->FacetNo);
  g_free(keyword);

  /* "xxxRAShift"      ofloat* RA shift (deg) for each image */
  if (prefix) keyword = g_strconcat (prefix, "RAShift", NULL);
  else        keyword = g_strdup("RAShift");
  dim[0] = in->numberImages; dim[1] = 1;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_float, dim, &in->RAShift);
  g_free(keyword);

  /* "xxxDecShift"     ofloat* Dec shift (deg) for each image */
  if (prefix) keyword = g_strconcat (prefix, "DecShift", NULL);
  else        keyword = g_strdup("DecShift");
  dim[0] = in->numberImages; dim[1] = 1;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_float, dim, &in->DecShift);
  g_free(keyword);

  /* "xxxfileType"     olong Are image OBIT_IO_AIPS or OBIT_IO_FITS? */
  if (prefix) keyword = g_strconcat (prefix, "fileType", NULL);
  else        keyword = g_strdup("fileType");
  dim[0] = 1;
  otemp = (olong)in->fileType;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_long, dim, &otemp);
  g_free(keyword);

  /* "xxximName"       string* Name of Mosaic images */
  if (prefix) keyword = g_strconcat (prefix, "imName", NULL);
  else        keyword = g_strdup("imName");
  dim[0] = strlen(in->imName);
  ObitInfoListAlwaysPut(outList, keyword, OBIT_string, dim, &in->imName);
  g_free(keyword);

  /* "xxximClass"      string* Class of Mosaic images */
  if (prefix) keyword = g_strconcat (prefix, "imClass", NULL);
  else        keyword = g_strdup("imClass");
  dim[0] = strlen(in->imClass);
  ObitInfoListAlwaysPut(outList, keyword, OBIT_string, dim, &in->imClass);
  g_free(keyword);

  /* "xxximSeq"        olong Sequence number of Mosaic images */
  if (prefix) keyword = g_strconcat (prefix, "imSeq", NULL);
  else        keyword = g_strdup("imSeq");
  dim[0] = 1;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_long, dim, &in->imSeq);
  g_free(keyword);

  /* "xxximDisk"       olong Disk number of Mosaic images  */
  if (prefix) keyword = g_strconcat (prefix, "imDisk", NULL);
  else        keyword = g_strdup("imDisk");
  dim[0] = in->numberImages; dim[1] = 1;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_long, dim, in->imDisk);
  g_free(keyword);

  /* "xxxdoFull"       boolean Is a full field image desired? */
  if (prefix) keyword = g_strconcat (prefix, "doFull", NULL);
  else        keyword = g_strdup("doFull");
  dim[0] = 1;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_bool, dim, &in->doFull);
  g_free(keyword);

  /* "xxxbmaj"         ofloat Restoring beam major axis in deg. */
  if (prefix) keyword = g_strconcat (prefix, "bmaj", NULL);
  else        keyword = g_strdup("bmaj");
  dim[0] = 1;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_float, dim, &in->bmaj);
  g_free(keyword);

  /* "xxxbmin"         ofloat Restoring beam minor axis in deg. */
  if (prefix) keyword = g_strconcat (prefix, "bmin", NULL);
  else        keyword = g_strdup("bmin");
  dim[0] = 1;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_float, dim, &in->bmin);
  g_free(keyword);

  /* "xxxbpa"          ofloat Restoring beam PA in deg. */
  if (prefix) keyword = g_strconcat (prefix, "bpa", NULL);
  else        keyword = g_strdup("bpa");
  dim[0] = 1;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_float, dim, &in->bpa);
  g_free(keyword);
 
  /* "xxxisAuto"       olong if >=0 the counterpart of an autoCenter field  */
  if (prefix) keyword = g_strconcat (prefix, "isAuto", NULL);
  else        keyword = g_strdup("isAuto");
  dim[0] = in->numberImages; dim[1] = 1;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_long, dim, in->isAuto);
  g_free(keyword);

  /* "xxxisShift"       olong if >=0 the counterpart of an autoCenter shifted field  */
  if (prefix) keyword = g_strconcat (prefix, "isShift", NULL);
  else        keyword = g_strdup("isShift");
  dim[0] = in->numberImages; dim[1] = 1;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_long, dim, in->isShift);
  g_free(keyword);

  /* "xxxnumBeamTapes"          . */
  if (prefix) keyword = g_strconcat (prefix, "numBeamTapes", NULL);
  else        keyword = g_strdup("numBeamTapes");
  dim[0] = 1;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_long, dim, &in->numBeamTapes);
  g_free(keyword);
 
  /* "xxxBeamTapes"          . */
  if (prefix) keyword = g_strconcat (prefix, "BeamTapes", NULL);
  else        keyword = g_strdup("BeamTapes");
  dim[0] = in->numBeamTapes; dim[1] = 1;
  if (in->BeamTapes)
    ObitInfoListAlwaysPut(outList, keyword, OBIT_float, dim, &in->BeamTapes);
  g_free(keyword);
 
  /* "xxxBeamTaper"          . */
  if (prefix) keyword = g_strconcat (prefix, "BeamTaper", NULL);
  else        keyword = g_strdup("BeamTaper");
  dim[0] = in->numberImages; dim[1] = 1;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_float, dim, &in->BeamTaper);
  g_free(keyword);
 
} /* end ObitImageMosaicGetInfo */

/**
 * Concatenate any CC files from facet images onto the FullField member
 * Uses CC version 1 in all cases.
 * If there are multiple tapered images, the following keywords are added to
 * the output table descriptor list:
 *   \li "NTAPER"   olong Number of tapers
 *   \li "TAPEnnnn" ofloat Taper nnnn in deg.
 * \param mosaic   ImageMosaic to process
 * \param err      Error/message stack
 */
void ObitImageMosaicCopyCC (ObitImageMosaic *in, ObitErr *err)
{
  ObitTableCC *inCCTab=NULL, *outCCTab=NULL;
  olong i, CCVer, noParms, ifield;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar keyword[9];
  gchar *routine = "ObitImageMosaicCopyCC";

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Check if the FullField member exists - give warning and return if not. */
  if ((in->FullField==NULL) && (in->nFlyEye>1)) {
    Obit_log_error(err, OBIT_InfoWarn, 
		   "%s: Mosaic %s has no FullField member", 
		   routine, in->name);
  }

  /* Get first input table to get structure */
  CCVer   = 1;
  noParms = 0;
  inCCTab = newObitTableCCValue ("in CC", (ObitData*)in->images[0],
				  &CCVer, OBIT_IO_ReadOnly, noParms, err);
  if  (err->error) Obit_traceback_msg (err, routine,in->images[0]->name);
  
  /* Create output CC table on in->FullField if exists */
  if (in->FullField) {
    CCVer    = 1;
    if (in->numBeamTapes>1) /* Promote points if tapering */
      noParms  = MAX (4, inCCTab->noParms);
    else
      noParms  = inCCTab->noParms;
    outCCTab = newObitTableCCValue ("out CC", (ObitData*)in->FullField,
				    &CCVer, OBIT_IO_ReadWrite, noParms, err);
    /* Clear any existing rows */
    ObitTableClearRows ((ObitTable*)outCCTab, err);
    if  (err->error) Obit_traceback_msg (err, routine, outCCTab->name);

    /* Copy first */
    ObitTableCCUtilAppend (inCCTab, outCCTab, 0, 0, err);
    if  (err->error) Obit_traceback_msg (err, routine, inCCTab->name);
    inCCTab = ObitTableCCUnref(inCCTab);
  } else { /* Copy others to first image */
    CCVer    = 1;
    noParms  = inCCTab->noParms;;
    outCCTab = newObitTableCCValue ("out CC", (ObitData*)in->images[0],
				    &CCVer, OBIT_IO_ReadWrite, noParms, err);
  }

  /* Loop over other tables */
  for (ifield=1; ifield<in->numberImages; ifield++) {
    CCVer   = 1;
    noParms = 0;
    inCCTab = newObitTableCCValue ("in CC", (ObitData*)in->images[ifield],
				   &CCVer, OBIT_IO_ReadOnly, noParms, err);
    if  (err->error) Obit_traceback_msg (err, routine, in->images[ifield]->name);
    if (inCCTab==NULL) continue;  /* Is there one? */
    ObitTableCCUtilAppend (inCCTab, outCCTab, 0, 0, err);
    if  (err->error) Obit_traceback_msg (err, routine, inCCTab->name);
    inCCTab = ObitTableCCUnref(inCCTab); /* Cleanup */
  } /* end loop copying rest of tables */

  /* Add taper info if needed */
  if (in->numBeamTapes>1) {
    ObitTableOpen ((ObitTable*)outCCTab, OBIT_IO_ReadWrite, err);
    if  (err->error) Obit_traceback_msg (err, routine, outCCTab->name);
    dim[0] = dim[1] = 1;
    /* Add to myIO as well */
    ObitInfoListAlwaysPut (outCCTab->myDesc->info, 
			   "NTAPER", OBIT_long, dim, &in->numBeamTapes);
    ObitInfoListAlwaysPut (((ObitTableDesc*)outCCTab->myIO->myDesc)->info, 
			   "NTAPER", OBIT_long, dim, &in->numBeamTapes);
    for (i=0; i<in->numBeamTapes; i++) {
      sprintf(keyword, "TAPE%4.4d", i+1);
      ObitInfoListAlwaysPut (outCCTab->myDesc->info, keyword, OBIT_float, 
			     dim, &in->BeamTaper[i]);
      ObitInfoListAlwaysPut (((ObitTableDesc*)outCCTab->myIO->myDesc)->info, 
			     keyword, OBIT_float, dim, &in->BeamTaper[i]);
    }
    outCCTab->myStatus = OBIT_Modified;
    ObitTableClose ((ObitTable*)outCCTab, err);
    if  (err->error) Obit_traceback_msg (err, routine, outCCTab->name);
  }

  /* Update number of CCs */
  ObitImageOpen(in->FullField, OBIT_IO_ReadWrite, err);
  in->FullField->myDesc->niter = outCCTab->myDesc->nrow;
  in->FullField->myStatus = OBIT_Modified;
  ObitImageClose(in->FullField, err);
  
  outCCTab = ObitTableCCUnref(outCCTab); /* Cleanup */
} /* end ObitImageMosaicCopyCC */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitImageMosaicClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitImageMosaicClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitImageMosaicClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitImageMosaicClassInfoDefFn (gpointer inClass)
{
  ObitImageMosaicClassInfo *theClass = (ObitImageMosaicClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitImageMosaicClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitImageMosaicClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitImageMosaicGetClass;
  theClass->newObit       = (newObitFP)newObitImageMosaic;
  theClass->ObitCopy      = (ObitCopyFP)ObitImageMosaicCopy;
  theClass->ObitClone     = NULL;  /* Different call */
  theClass->ObitClear     = (ObitClearFP)ObitImageMosaicClear;
  theClass->ObitInit      = (ObitInitFP)ObitImageMosaicInit;
  theClass->ObitImageMosaicFromInfo = (ObitImageMosaicFromInfoFP)ObitImageMosaicFromInfo;
  theClass->ObitImageMosaicGetInfo  = (ObitImageMosaicGetInfoFP)ObitImageMosaicGetInfo;
  theClass->ObitImageMosaicZapImage =
    (ObitImageMosaicZapImageFP)ObitImageMosaicZapImage;
  theClass->ObitImageMosaicGetImage = 
    (ObitImageMosaicGetImageFP)ObitImageMosaicGetImage;
  theClass->ObitImageMosaicSetImage = 
    (ObitImageMosaicSetImageFP)ObitImageMosaicSetImage;
  theClass->ObitImageMosaicGetImageRMS = 
    (ObitImageMosaicGetImageRMSFP)ObitImageMosaicGetImageRMS;
  theClass->ObitImageMosaicGetFullImage = 
    (ObitImageMosaicGetFullImageFP)ObitImageMosaicGetFullImage;
  theClass->ObitImageMosaicSetFullImage = 
    (ObitImageMosaicSetFullImageFP)ObitImageMosaicSetFullImage;
  theClass->ObitImageMosaicSetFiles = 
    (ObitImageMosaicSetFilesFP)ObitImageMosaicSetFiles;
  theClass->ObitImageMosaicCreate = 
    (ObitImageMosaicCreateFP)ObitImageMosaicCreate;
  theClass->ObitImageMosaicDefine = 
    (ObitImageMosaicDefineFP)ObitImageMosaicDefine;
  theClass->ObitImageMosaicFlatten = 
    (ObitImageMosaicFlattenFP)ObitImageMosaicFlatten;
  theClass->ObitImageMosaicFOV = 
    (ObitImageMosaicFOVFP)ObitImageMosaicFOV;
  theClass->ObitImageMosaicMaxFOV = 
    (ObitImageMosaicMaxFOVFP)ObitImageMosaicMaxFOV;
  theClass->ObitImageMosaicMaxCC = 
    (ObitImageMosaicMaxCCFP)ObitImageMosaicMaxCC;
  theClass->ObitImageMosaicCombineCC = 
    (ObitImageMosaicCombineCCFP)ObitImageMosaicCombineCC;
  theClass->ObitImageMosaicFlagCC = 
    (ObitImageMosaicFlagCCFP)ObitImageMosaicFlagCC;
  theClass->ObitImageMosaicAddField = 
    (ObitImageMosaicAddFieldFP)ObitImageMosaicAddField;
  theClass->ObitImageMosaicMaxField = 
    (ObitImageMosaicMaxFieldFP)ObitImageMosaicMaxField;
  theClass->ObitImageMosaicCopyCC = 
    (ObitImageMosaicCopyCCFP)ObitImageMosaicCopyCC;
  /* Private functions for derived classes */
  theClass->FlyEye        = (FlyEyeFP)FlyEye;
  theClass->AddField      = (AddFieldFP)AddField;
  theClass->AddOutlier    = (AddOutlierFP)AddOutlier;
  theClass->AddTapers     = (AddTapersFP)AddTapers;

} /* end ObitImageMosaicClassDefFn */


/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitImageMosaicInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  olong i;
  ObitImageMosaic *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->thread    = newObitThread();
  in->info      = newObitInfoList(); 
  in->FullField = NULL;
  in->doFull    = FALSE;
  in->images    = ObitMemAlloc0Name(in->numberImages*sizeof(ObitImage*),"ImageMosaic images");
  for (i=0; i<in->numberImages; i++) in->images[i] = NULL;
  in->inFlysEye = ObitMemAlloc0Name(in->numberImages*sizeof(gboolean),  "ImageMosaic inFlysEye");
  for (i=0; i<in->numberImages; i++) in->inFlysEye[i] = FALSE;
  in->FacetNo   = ObitMemAlloc0Name(in->numberImages*sizeof(olong),  "ImageMosaic FacetNo");
  in->imDisk    = ObitMemAlloc0Name(in->numberImages*sizeof(olong)," ImageMosaic imDisk");
  in->nx        = ObitMemAlloc0Name(in->numberImages*sizeof(olong), "ImageMosaic nx");
  in->ny        = ObitMemAlloc0Name(in->numberImages*sizeof(olong), "ImageMosaic ny");
  in->nplane    = ObitMemAlloc0Name(in->numberImages*sizeof(olong), "ImageMosaic nplane");
  in->RAShift   = ObitMemAlloc0Name(in->numberImages*sizeof(ofloat),"ImageMosaic RAShift");
  in->DecShift  = ObitMemAlloc0Name(in->numberImages*sizeof(ofloat),"ImageMosaic DecShift");
  in->isAuto    = ObitMemAlloc0Name(in->numberImages*sizeof(olong), "ImageMosaic isAuto");
  in->isShift   = ObitMemAlloc0Name(in->numberImages*sizeof(olong), "ImageMosaic isShift");
  in->BeamTaper = ObitMemAlloc0Name(in->numberImages*sizeof(ofloat),"ImageMosaic BeamTaper");
  for (i=0; i<in->numberImages; i++) in->BeamTaper[i] = 0.0;
  in->imName    = NULL;
  in->imClass   = NULL;
  in->numBeamTapes = 0;
  in->BeamTapes = NULL;
  in->bmaj      = 0.0;
  in->bmin      = 0.0;
  in->bpa       = 0.0;
  in->Radius    = 0.0;
  in->OutlierSize= 0;

} /* end ObitImageMosaicInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitImageMosaic* cast to an Obit*.
 */
void ObitImageMosaicClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitImageMosaic *in = inn;
  olong i;
  
  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  
  /* delete this class members */
  in->thread    = ObitThreadUnref(in->thread);
  in->info      = ObitInfoListUnref(in->info);
  in->FullField = ObitUnref(in->FullField);
  if (in->images) {
    for (i=0; i<in->numberImages; i++)
      in->images[i] = ObitUnref(in->images[i]);
    in->images = ObitMemFree(in->images); 
  }
  if (in->imDisk)    ObitMemFree(in->imDisk); in->imDisk = NULL;
  if (in->nx)        ObitMemFree(in->nx);     in->nx     = NULL;
  if (in->ny)        ObitMemFree(in->ny);     in->ny     = NULL;
  if (in->inFlysEye) ObitMemFree(in->inFlysEye);in->inFlysEye = NULL;
  if (in->FacetNo)   ObitMemFree(in->FacetNo);in->FacetNo = NULL;
  if (in->nplane)    ObitMemFree(in->nplane); in->nplane = NULL;
  if (in->RAShift)   ObitMemFree(in->RAShift);in->RAShift=NULL;
  if (in->DecShift)  ObitMemFree(in->DecShift);in->DecShift=NULL;
  if (in->isAuto)    ObitMemFree(in->isAuto); in->isAuto = NULL;
  if (in->isShift)   ObitMemFree(in->isShift);in->isShift= NULL;
  if (in->BeamTapes) ObitMemFree(in->BeamTapes);in->BeamTapes = NULL;
  if (in->BeamTaper) ObitMemFree(in->BeamTaper);in->BeamTaper = NULL;
  if (in->imName)    g_free(in->imName);      in->imName = NULL;
  if (in->imClass)   g_free(in->imClass);     in->imClass= NULL;
 
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitImageMosaicClear */

/**
 * Set fly's eye fields using a hexagonal pattern of overlapping, circular fields.
 * Translated from the AIPSish FLYEYE in $FOURMASS/SUB/ADDFIELDS.FOR
 * \param  radius  radius of fly's eye in degree 
 * \param  imsize  size of circular regions to be imaged 
 *                 this is the size of the undistorted region. 
 * \param  cells   pixel size (asec) 
 * \param  overlap number pixels of overlap 
 * \param  shift   shift all coordinates 
 * \param  ra0     center ra in degrees 
 * \param  dec0    center dec in degrees 
 * \param  nfield  (output) number fields already added 
 * \param  fldsiz  (output) circular image field sizes (pixels) 
 * \param  rash    (output) ra shifts (asec) 
 * \param  decsh   (output) declination shifts (asec) 
 * \param  flqual  (output) field quality (crowding) code, -1 => ignore 
 * \param  err     Error stack, returns if not empty.
 * \return radius of zone of avoidence for externals (deg)
 */
ofloat
FlyEye (ofloat radius, olong imsize, ofloat cells[2], olong overlap,
	ofloat shift[2], odouble ra0, odouble dec0, 
	olong *nfield, olong *fldsiz, ofloat *rash, ofloat *decsh, olong *flqual,
	ObitErr *err) 
{
  ofloat out = 0.0;
  odouble cosdec, sindec, cosfld, sinfld, ra[MAXFLD], dec[MAXFLD], dx, dy, 
    drad, xra, xdec, xd,  xsh[2], ll, mm, maxrad, ra0r, dec0r, dradd;
  ofloat      dist[MAXFLD], mrotat, rxsh[2];
  gboolean   this1, warn;
  olong   ifield, mfield, ii, jj, indx[MAXFLD], i1, i2, j1, j2, jerr, nfini, mxfld;
  gchar *routine = "FlyEye";

  /* error checks */
  if (err->error) return out;

  for (ii=0; ii<MAXFLD; ii++) dist[ii] = -1.0;
  maxrad = 0.0;
  cosdec = cos (dec0*DG2RAD);
  sindec = sin (dec0*DG2RAD);
  ra0r   = ra0 * DG2RAD;
  dec0r  = dec0 * DG2RAD;
  nfini  = *nfield;
  warn   = FALSE;

  /* shift to degrees */
  xsh[1] = shift[1] / 3.6e3;
  xsh[0] = shift[0] / 3.6e3 / cosdec;

  /* spacing of imaging region in  rad.   */
  drad = 0.5 * AS2RAD * fabs(cells[0]) * (imsize - overlap);
  /* a bit less spacing to allow for sky curvature */
  drad *= 1.0 / (1. + tan(drad));
  dx = 1.5 * drad;
  dy = cos (30.0*DG2RAD) * drad;

  /* field size in deg. */
  dradd = drad * RAD2DG / sqrt (2.0);
  dradd *= 1.1;  /* add a bit to be fill in potential gaps */

  /* collect list of fields to add */
  ifield = *nfield;

  /* maximum number of fields to add, only use inscribed circle in square array. */
  mfield = 2.0 * sqrt (MAXFLD / G_PI);

  /* don't go overboard - only twice  radius */
  mxfld = 0.5 + 2.0 * (radius) / dradd;
  mxfld = MAX (1, mxfld);
  mfield = MIN (mfield, mxfld);
  i2 = mfield;
  i1 = -i2;
  j2 = mfield/2 + 1;
  j1 = -j2;
  this1 = FALSE;

  /* See if one field is enough? */
  if (radius < (0.5*imsize*fabs(cells[0])/3600.0)) {  /* Yes one will do */
    i1 = i2 = 0;
    maxrad = radius;
  } /* end of one field is enough */

 /* Loop over declination */
  for (ii= i1; ii<=i2; ii++) { /* loop 20 */
    mm = ii * dy;

    /* Loop over RA */
    for (jj= i1; jj<=i2; jj++) { /* loop 10 */
      ll = jj * dx;
      /* do every other possibility to get hexagonal pattern */
      this1 = !this1;
      if (this1) {

	/* Convert to RA, Dec (sine projection)  */
	ObitSkyGeomNewPos (OBIT_SkyGeom_SIN, ra0r, dec0r, ll, mm, &xra, &xdec, &jerr);
	if (jerr != 0) {
	  /* angle too large for sine projection. */
	  xra  = ra0 + 2.0 * (radius);
	  xdec = dec0 - 2.0 * (radius);
	  xd   = 1000.0 * (fabs(radius)+dradd);
	} else {
	  /* ok - convert to deg. */
	  xra  = xra*RAD2DG  + xsh[0];
	  xdec = xdec*RAD2DG + xsh[1];
	  /* distance from center */
	  cosfld = cos (xdec*DG2RAD);
	  sinfld = sin (xdec*DG2RAD);
	  xd = sinfld * sindec + cosfld * cosdec * cos (xra*DG2RAD - ra0r);
	  if (xd > 1.0) {
	    xd = 0.0;
	  } else {
	    xd = RAD2DG * acos (MIN (xd, 1.000));
	  } 
	} 

	/* Add field to list if close enough  */
	/* debug 
	fprintf (stderr,"%d %lf %lf %lf\n",ifield,ll,mm,xd);*/
	if (xd < (fabs(radius)+dradd)) {
	  if (ifield < MAXFLD) {
	    ra[ifield]     = xra;
	    dec[ifield]    = xdec;
	    dist[ifield]   = xd;
	    fldsiz[ifield] = imsize;
	    flqual[ifield] = 0;
	    indx[ifield]   = ifield;
	    ifield++;
	    maxrad = MAX (maxrad, xd);
	  } else {
	    warn = TRUE;
	  } 
	} /* end add field that's close enough */
      } /* end pick every other entry to get hexagonal pattern */
    } /* end loop  L10 over RA:  */
  } /* end loop  L20 over declination:  */

  mfield = ifield;

  /* sort on distances */
  Bubble (dist, indx, ifield, 1);

  /* save to output */
  if (radius > 0.0) {
    
    mrotat = 0.0;/* assume no systematic rotation */
    for (ii = 0; ii<mfield; ii++) { /* loop 30 */
      ObitSkyGeomShiftXY (ra0, dec0, mrotat, ra[indx[ii]], dec[indx[ii]], 
			  &rxsh[0], &rxsh[1]);
      AddField (rxsh, dec[indx[ii]], imsize, cells, -1, FALSE, 1.0e20,
		nfield, fldsiz, rash, decsh, flqual, err);
    } /* end loop  L30:  */
  } 

  /* tell how many */
  mfield = *nfield - nfini;
  Obit_log_error(err, OBIT_InfoErr, "%s: added %d flys eye fields", routine, mfield);

  /* If no fields complain and die */
  if ((*nfield)<=0) {
  Obit_log_error(err, OBIT_Error, "%s: NO fields to be imaged ",routine);
  }

  /* too many? */
  if (warn) 
    Obit_log_error(err, OBIT_InfoWarn, 
		   "%s: Exceeded field limit - full field will not be imaged", 
		   routine);

  /* want confusing sources outside  of this */
  out = maxrad;
  return out;
} /* end of routine FlyEye */ 


/**
 * Add a field the the field list if there is room and if it not within 3/4 
 * the radius of another field.
 * Translated from the AIPSish ADFLDX in $FOURMASS/SUB/ADDFIELDS.FOR
 * \param   shift    x,y shift in deg 
 * \param   dec      declination of new field (deg) 
 * \param   imsize   desired size as diameter in pixels of field 
 * \param   cells    cell size in arc sec 
 * \param   qual     input field quality (crowding) code. 
 * \param   check    if true check if this is in another field 
 * \param  minImpact if check then MIN (minImpact, 0.75*fldsize) is the 
 *                   distance in cells to consider a match
 * \param   nfield   (output) number fields already added 
 * \param   fldsiz   (output) circular image field sizes (pixels) 
 * \param   rash     (output) ra shifts (asec) 
 * \param   decsh    (output) declination shifts (asec) 
 * \param   flqual   (output) field quality (crowding) code, -1 => ignore 
 * \param   err      Error stack, returns if not empty.
 * \return 0 if added, 1=matches previous, 2=table full
 */
olong 
AddField (ofloat shift[2], ofloat dec, olong imsize, ofloat cells[2], 
	  olong qual, gboolean check, ofloat minImpact,
	  olong *nfield, olong *fldsiz, ofloat *rash, ofloat *decsh, olong *flqual,
	  ObitErr *err) 
{
  olong   i;
  ofloat      dist, maxdis, cosdec, xsh, ysh, xxsh, yysh, mindist;
  gchar *routine = "AddField";
  
  /* error checks */
  if (err->error) return 1;

  /* is this in a previous field? */
  if (check) {
    cosdec = cos (dec*DG2RAD);
    xsh = shift[0]*3600.0;
    ysh = shift[1]*3600.0;
    for (i = 0; i<(*nfield); i++) { /* loop 100 */
      xxsh = rash[i];
      yysh = decsh[i];
      /* work in pixels**2 */
      dist = (cosdec * (xxsh - xsh) / cells[0]) * 
	(cosdec * (xxsh - xsh) / cells[0]) + 
	((yysh - ysh) / cells[1]) * ((yysh - ysh) / cells[1]);
      /* Minimum distance for a match */
      mindist = MIN (minImpact, 0.5*0.75*fldsiz[i]);
      /* distance in pixel**2 to consider a match */
      maxdis = mindist*mindist;
      if (dist < maxdis) {
	Obit_log_error(err, OBIT_InfoWarn, 
		       "%s: Rejecting additional field at %g %g, in previous field %d", 
		       routine, shift[0], shift[1], i+1);
	return 1;
      } /* end of already in another field */
    } /* end loop  L100: */
  } 

  /* if you got here - it must be ok */
  if ((*nfield) < MAXFLD) {
    (*nfield)++;
    fldsiz[(*nfield)-1] = ObitFFTSuggestSize (imsize); /* Make it FFT friendly */
    rash[(*nfield)-1]   = shift[0]*3600.0;
    decsh[(*nfield)-1]  = shift[1]*3600.0;
    flqual[(*nfield)-1] = qual;
  } else {
    /* can't fit it */
    return 2;
  } 
  return 0;
} /* end of routine  AddField */ 

/**
 * In place bubble sort with tracking of swap indices
 * \param data    data array to sort, returned sorted 
 * \param indx    [out] indices (0-rel) for sorted data.
 * \param number  No. elements in data, indx
 * \param direct  sort ordet, 1=ascending, -1=descending
 */
 static void 
 Bubble (ofloat *data, olong* indx, olong number, olong direct)
{
  olong   i, j, it;
  ofloat      rt;
  gboolean   done;
  
  done = FALSE;
  i = 1;
  
  /* descending order sort */
  if (direct == -1) {
    while ((i < number)  &&  (!done)) {
      done = TRUE;
      for (j= number; j>=i+1; j--) { /* loop 20 */
	if (data[j-1] > data[j-2]) {
	  rt = data[j-1];
	  data[j-1] = data[j-2];
	  data[j-2] = rt;
	  it = indx[j-1];
	  indx[j-1] = indx[j-2];
	  indx[j-2] = it;
	  done = FALSE;
	} 
      } /* end loop  L20:  */;
      i++;
    } /* end while loop */
    
    /* ascending order sort */
  } else {
    while ((i < number)  &&  (!done)) {
      done = TRUE;
      for (j= number; j>=i+1; j--) { /* loop 40 */
	if (data[j-1] < data[j-2]) {
	  rt = data[j-1];
	  data[j-1] = data[j-2];
	  data[j-2] = rt;
	  it = indx[j-1];
	  indx[j-1] = indx[j-2];
	  indx[j-2] = it;
	  done = FALSE;
	} 
      } /* end loop  L40:  */;
      i++;
    } /* end while loop */
  } 

} /* end of routine Bubble */ 

/**
 * Searches Catalog for sources in the desired ra, dec, and flux 
 * range taking into account the estimated single-dish beam 
 * Translated from the AIPSish ADNVSS in $FOURMASS/SUB/ADDFIELDS.FOR
 * \param Catalog      FITS AIPSVZ format catalog file name
 * \param catDisk      FITS sidk number for Catalog
 * \param minRad       Minimum distance from field center to include
 * \param cells        pixel size (asec)
 * \param OutlierDist  Maximum distance (deg) from center - if zero return
 * \param OutlierFlux  Minimum estimated flux density 
 * \param OutlierSI    Spectral index to use to convert catalog flux density
 * \param OutlierSize  Width of outlier field in pixels
 * \param ra0          center ra in degrees 
 * \param dec0         center dec in degrees 
 * \param doJ2B        If True precess catalog from J2000 to B1950
 * \param Freq         Frequency of data in Hz
 * \param minImpact    MIN (minImpact, 0.75*fldsiz) is the distance in cells to 
 *                     consider a match with a previous field
 * \param diam         (input) Antenna diameter (m)
 * \param nfield       (output) number fields already added
 * \param fldsiz       (output) circular image field sizes (pixels) 
 * \param rash         (output) ra shifts (asec)
 * \param decsh        (output) declination shifts (asec) 
 * \param flqual       (output) field quality (crowding) code, -1 => ignore 
 * \param err          Error stack, returns if not empty.
 */
void 
AddOutlier (gchar *Catalog, olong catDisk, ofloat minRad, ofloat cells[2], 
	    ofloat OutlierDist, ofloat OutlierFlux, ofloat OutlierSI, olong OutlierSize,
	    odouble ra0, odouble dec0, gboolean doJ2B, odouble Freq, ofloat minImpact, 
	    ofloat diam,  
	    olong *nfield, olong *fldsiz, ofloat *rash, ofloat *decsh, olong *flqual, 
	    ObitErr *err) 
{
  olong count, imsize, jerr, nfini, qual;
  odouble ra, dec, ra2000, dc2000, dra;
  odouble xx, yy, zz, dist, refreq;
  ofloat radius, minflx, asize, alpha, pbf;
  ofloat flux, scale, xsh[2];
  gboolean wanted, doJinc, warn;
  olong blc[IM_MAXDIM] = {1,1,1,1,1};
  olong trc[IM_MAXDIM] = {0,0,0,0,0};
  olong ver, nrows, irow;
  ObitIOCode retCode;
  ObitImage *VZImage=NULL;
  ObitTableVZ *VZTable=NULL;
  ObitTableVZRow *VZRow=NULL;
  gchar *routine = "AddOutlier";
  
  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;

  /* Really needed? */
  if (OutlierDist<=1.0e-10) return;
  
  /* get control parameters */
  radius = OutlierDist;
  minflx = OutlierFlux;
  alpha  = OutlierSI;
  imsize = OutlierSize;
  if (imsize <= 0) imsize = 50;
  asize = diam;      /* antenna diameter in meters */
  nfini = *nfield;
  warn = TRUE;

  /* get j2000 position to lookup in Catalog in radians */
  ra2000 = ra0;
  dc2000 = dec0;
  if (doJ2B) ObitSkyGeomBtoJ (&ra2000, &dc2000);
  ra2000 *= DG2RAD;
  dc2000 *= DG2RAD;

  /* set defaults. */
  if (radius <= 0.0) radius = 15.0;
  if (asize <= 0.0)  asize  = 25.0;
  if (alpha == 0.0)  alpha  = -0.75;

  /* which beam model to use */
  doJinc = (Freq >= 1.0e9);

  /* Open Catalog (VZ table on an image) */
  VZImage = newObitImage("Catalog image");
  ObitImageSetFITS(VZImage, OBIT_IO_byPlane, catDisk, Catalog, blc, trc, err);

  /* Open to fully instantiate */
  ObitImageOpen(VZImage, OBIT_IO_ReadOnly, err);
  if (err->error) Obit_traceback_msg (err, routine, VZImage->name);

  /* Now get VZ table */
  ver = 1;
  VZTable =  newObitTableVZValue("Catalog table", (ObitData*)VZImage, &ver, 
				 OBIT_IO_ReadOnly, err);
  ObitTableVZOpen(VZTable, OBIT_IO_ReadOnly, err);
  VZRow =  newObitTableVZRow (VZTable);  /* Table row */
  if (err->error) Obit_traceback_msg (err, routine, VZTable->name);

  /* Get table info */
  refreq = VZTable->refFreq;
  nrows  = VZTable->myDesc->nrow;

  /* frequency scaling */
  scale = pow ((Freq / refreq), alpha);

  /* loop through table */
  count = 0;
  warn = FALSE;
  for (irow= 1; irow<=nrows; irow++) { /* loop 500 */
    /* read */
    retCode = ObitTableVZReadRow (VZTable, irow, VZRow, err);
    if (err->error) Obit_traceback_msg (err, routine, VZTable->name);
   
    /* spectral scaling of flux density */
    flux = VZRow->PeakInt * scale;

    /* position, etc */
    ra   = VZRow->Ra2000;
    dec  = VZRow->Dec2000;
    qual = VZRow->Quality;

    /* select (crude) */
    xx = DG2RAD * ra;
    yy = DG2RAD * dec;
    dra = fabs (ra0-ra);
    if (dra>180.0) dra -= 360.0;
    if ((fabs(dec0-dec) <= 1.2*radius) && 
	(fabs(dra)*cos(yy) <= 1.2*radius) &&  
	(flux >= minflx)) {
      /* separation from pointing center */
      zz = sin (yy) * sin (dc2000) + cos (yy) * cos (dc2000) * cos (xx-ra2000);
      zz = MIN (zz, 1.000);
      dist = acos (zz) * RAD2DG;

      /* primary beam correction to flux density */
      if (doJinc) {
	pbf = ObitPBUtilJinc (dist, Freq, asize, 0.05);
      } else {
	pbf = ObitPBUtilPoly (dist, Freq, 0.05);
      } 
      flux *= MAX (0.05, pbf); /* Don't trust below 5% */
      
      /* select (fine) */
      wanted = ((flux >= minflx)  &&  (dist <= radius))  && (dist > minRad);
      if (wanted) {
	if (doJ2B) {  /* precess to 1950 if necessary */
	  ObitSkyGeomJtoB (&ra, &dec);
	} 

	/* get shift needed */
	ObitSkyGeomShiftXY (ra0, dec0, 0.0, ra, dec, &xsh[0], &xsh[1]);

	/* add it if you can */
	jerr = AddField (xsh, dec, imsize, cells, qual, TRUE, minImpact,
			 nfield, fldsiz, rash, decsh, flqual, err);
	warn = warn  ||  (jerr != 0);
       
	/* debug */
	if (jerr==0) Obit_log_error(err, OBIT_InfoErr, "Field %d is %f %f", *nfield, ra, dec);
      } /* end if wanted */ 
    } /* end crude selection */
  } /* end loop over table */

  /* tell how many */
  count = *nfield - nfini;
  Obit_log_error(err, OBIT_InfoErr, "%s: Added %d outlying fields", routine, count);

  /* too many? */
  if ((*nfield)>=MAXFLD) Obit_log_error(err, OBIT_InfoWarn, 
					"%s: Hit field limit - some outliers may be ignored", 
					routine);
  /* Close up */
  ObitImageClose(VZImage, err);
  retCode = ObitTableVZClose(VZTable, err);
  if (err->error) Obit_traceback_msg (err, routine, VZTable->name);
  
  /* clean up */
  VZImage = ObitImageUnref(VZImage);
  VZTable = ObitTableUnref(VZTable);
  VZRow   = ObitTableRowUnref(VZRow);
  
} /* end of routine AddOutlier */ 
/**
 * Checks if additional tapers are requested for images and if so
 * duplicates all previoiusly defined facets
 * \param uvData  The object to create images in,  Details are defined in info members:
 * \li numBeamTapes= Number of Beam tapers (elements in BeamTapes) 
 * \li BeamTapes   = List of Beam tapers  larger than 0
 * \param numBeamTapes (output) Total number of beam tapers including no taper
 * \param BeamTapes    (output) List of additional (nonzero) beam tapers in x pixels
 * \param cells        Image x pixel size in asec.
 * \param Tapers       (output) Beam taper per facet;  input initial facet 
 *                        values are zeroed
 * \param nfield       (in/output) number fields already added
 * \param fldsiz       (in/output) circular image field sizes (pixels) 
 * \param rash         (in/output) ra shifts (asec)
 * \param decsh        (in/output) declination shifts (asec) 
 * \param flqual       (in/output) field quality (crowding) code, -1 => ignore 
 * \param inFlysEye    (in/output) number fields already added
 * \param FacetNo      (in/output) Untapered facet number
 * \param err          Error stack, returns if not empty.
 */
void 
AddTapers (ObitUV *uvData, 
	   olong *numBeamTapes, ofloat *BeamTapes, ofloat cells, ofloat *Tapers, 
	   olong *nfield, olong *fldsiz, ofloat *rash, ofloat *decsh, 
	   olong *flqual, gboolean *inFlysEye, olong *FacetNo, ObitErr *err) 
{
  olong i, j, indx, nf, nt=0;
  gboolean blewIt;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *routine = "AddTapers";

  if (err->error) return;  /* previous error? */

  /* Initial previous facets */
  for (i=0; i<(*nfield); i++) Tapers[i] = 0.0;
  *numBeamTapes = 1;
  BeamTapes[0]  = 0.0;

  /* Get requests for new */
  if (!ObitInfoListGetTest(uvData->info, "numBeamTapes", &type, dim, &nt))
    return;
  *numBeamTapes = nt + 1;
  if (!ObitInfoListGetTest(uvData->info, "BeamTapes", &type, dim, &BeamTapes[1]))
    return;
  
  /* Loop over additional Tapers */
  nf = *nfield;
  blewIt = FALSE;  /* OK */
  for (j=0; j<nt; j++) {
    /* Loop over previously defined images */
    for (i=0; i<nf; i++) {
      indx = *nfield;
      fldsiz[indx]    = fldsiz[i];
      rash[indx]      = rash[i];
      decsh[indx]     = decsh[i];
      flqual[indx]    = flqual[i];
      inFlysEye[indx] = inFlysEye[i];
      FacetNo[indx]   = FacetNo[i];
      Tapers[indx]    = (cells/3600.0)*BeamTapes[j+1];
      if ((*nfield) < MAXFLD) {
	(*nfield)++;  /* Count of total facets */
      } else {  /* Blew array */
	blewIt = TRUE;
      }
    } /* end loop over previously defined image  */
  } /* end loop over additional tapers */

  /* too many? */
  if (blewIt) 
    Obit_log_error(err, OBIT_InfoWarn, 
		   "%s: Exceeded field limit - not all tapers will be included", 
		   routine);

  /* tell how many tapers added */
  Obit_log_error(err, OBIT_InfoErr, "%s: Added %d beam tapers", routine, nt);

} /* end AddTapers */
