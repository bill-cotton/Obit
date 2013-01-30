/* $Id$  */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2010-2013                                          */
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

#include "ObitImageMosaicWB.h"
#include "ObitIOImageFITS.h"
#include "ObitIOImageAIPS.h"
#include "ObitSystem.h"
#include "ObitImageUtil.h"
#include "ObitImageWB.h"
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
 * \file ObitImageMosaicWB.c
 * ObitImageMosaicWB class function definitions.
 * This class is derived from the Obit base class.
 *
 * This class contains an array of associated astronomical images.
 * SW wideband version
 *
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitImageMosaicWB";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitImageMosaicGetClass;

/*--------------- File Global Variables  ----------------*/
/**
 * ClassInfo structure ObitImageMosaicWBClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitImageMosaicWBClassInfo myClassInfo = {FALSE};

/*--------------------------- Macroes ---------------------*/
#ifndef MAXFLD         /* Maxum number of fields */
#define MAXFLD 4096
#endif

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitImageMosaicWBInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitImageMosaicWBClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitImageMosaicWBClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name   An optional name for the object.
 * \param number Number of images
 * \return the new object.
 */
ObitImageMosaicWB* newObitImageMosaicWB (gchar* name, olong number)
{
  ObitImageMosaicWB* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitImageMosaicWBClassInit();

  /* allocate/init structure */
  out = ObitMemAlloc0Name(sizeof(ObitImageMosaicWB), "ObitImageMosaicWB");

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  out->numberImages = number;
  ObitImageMosaicWBInit((gpointer)out);

 return out;
} /* end newObitImageMosaicWB */

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
 *      \li "xxxnorder"       olong maximum beam order [def 1]
 *      \li "xxxnumBeamTapes" olong  Number of Beam tapers (elements in BeamTapes) 
 *      \li "xxxBeamTapes"    ofloat*  List of Beam tapers  
 *      \li "xxxBeamTaper"    ofloat*  Beam Taper per image as FWHM in deg      
 * \param err     ObitErr for reporting errors.
 * \return the new object.
 */
ObitImageMosaicWB* ObitImageMosaicWBFromInfo (gchar *prefix, ObitInfoList *inList, 
					      ObitErr *err)
{ 
  ObitImageMosaicWB *out = NULL;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *keyword=NULL, *None = "None", *value=NULL;
  olong i, numberImages, otemp, norder;
  gboolean missing;
  gchar *routine = "ObitImageMosaicWBFromInfo";

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitImageMosaicWBClassInit();

  /* error checks */
  if (err->error) return out;

  /* "xxxnorder"       olong maximum beam order [def 1] */
  if (prefix) keyword = g_strconcat (prefix, "norder", NULL);
  else        keyword = g_strdup("norder");
  ObitInfoListGetTest(inList, keyword, &type, dim, &norder);
  g_free(keyword);

  /* "xxxnumberImages" olong Number of images in Mosaic; */
  if (prefix) keyword = g_strconcat (prefix, "numberImages", NULL);
  else        keyword = g_strdup("numberImages");
  ObitInfoListGetTest(inList, keyword, &type, dim, &numberImages);
  g_free(keyword);

  /* Create output */
  out = newObitImageMosaicWB (prefix, numberImages);
  out->norder = norder;   /* Max beam order */

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
} /* end ObitImageMosaicWBFromInfo */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitImageMosaicWBGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitImageMosaicWBClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitGetImageMosaicWBClass */

/**
 * Zap (delete with underlying structures) selected image member(s).
 * Also deletes any associated beam member
 * \param inn     The array of images, infoList has following:
 * \li "saveBeam" boolean  If TRUE, save beams, [def FALSE]
 * \param number  The 0-rel image number, -1=> all
 * \param err     Error stack, returns if not empty.
 */
void 
ObitImageMosaicWBZapImage  (ObitImageMosaic *inn, olong number,
			    ObitErr *err)
{
  olong i;
  ObitImageMosaicWB *in = (ObitImageMosaicWB*)inn;
  gchar *routine="ObitImageMosaicWBZapImage";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Kill 'em all? */
  if (number==-1) {
    for (i=0; i<in->numberImages; i++) {
      /* The image and all beams */
      in->images[i] = ObitImageWBZap(in->images[i], err);
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

  /* The image and all beams */
  in->images[number] = ObitImageWBZap(in->images[number], err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  
  return;
} /* end ObitImageMosaicWBZapImage */

/**
 * Make a shallow copy of input object.
 * Parent class members are included but any derived class info is ignored.
 * \param inn  The object to copy
 * \param outt An existing object pointer for output or NULL if none exists.
 * \param err  Error stack, returns if not empty.
 * \return pointer to the new object.
 */
ObitImageMosaic* ObitImageMosaicWBCopy (ObitImageMosaic *inn, ObitImageMosaic *outt, 
					ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  gboolean oldExist;
  olong i;
  gchar *outName;
  ObitImageMosaicWB *in  = (ObitImageMosaicWB*)inn;
  ObitImageMosaicWB *out = (ObitImageMosaicWB*)outt;
  gchar *routine = "ObitImageMosaicWBCopy";

  /* error checks */
  if (err->error) return NULL;
  g_assert (ObitIsA(in, &myClassInfo));
  if (out) g_assert (ObitIsA(out, &myClassInfo));

  /* Create if it doesn't exist */
  oldExist = out!=NULL;
  if (!oldExist) {
    /* derive object name */
    outName = g_strconcat ("Copy: ",in->name,NULL);
    out = newObitImageMosaicWB(outName, in->numberImages);
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
      return outt;
    }

    /* unreference any old members */
    for (i=0; i<out->numberImages; i++)
      out->images[i] = ObitImageUnref(out->images[i]);
  }

  /* Copy other class members */
  out->norder   = in->norder;
  return outt;
} /* end ObitImageMosaicWBCopy */

/**
 * Attach images to their underlying files.
 * For AIPS files:
 * Mosaic images have classes with the first character of the imClass
 * followed by 'M', followed by 4 digits of the field number.  
 * The Beams are the same except that the second character is 'B','C','D'...
 * for the varions beams.
 * The full field image has class imClass and 'F' as the 6th character.
 *
 * For FITS files:
 * Image classes are imClass+digits of the field
 * Beam classes are the same except the second character is replaced 
 * with a 'B' (unless it already is 'B' in which case 'b' is used).
 * Higher order beams use 'C' ('c'), 'D' ('d')...
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
void ObitImageMosaicWBSetFiles  (ObitImageMosaic *inn, gboolean doBeam, ObitErr *err) 
{
  olong i, j, user, cno;
  olong blc[IM_MAXDIM] = {1,1,1,1,1};
  olong trc[IM_MAXDIM] = {0,0,0,0,0};
  gboolean exist;
  ObitImage *image=NULL;
  gchar ct, strTemp[48], Aname[13], Aclass[7], Atclass[3], Atype[3] = "MA";
  ObitImageMosaicWB *in = (ObitImageMosaicWB*)inn;
  gchar beamClass[]="BCDEFGHI", altbeamClass[] = "bcdefghi";
  gchar *routine = "ObitImageMosaicWBSetFiles";

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
 
    /* Doing beams? */
    if (doBeam) {
      image = (ObitImage*)in->images[i]->myBeam;
      /* Define files - same except second character of "class" is 'B'... */
      if (in->fileType==OBIT_IO_AIPS) {
	/* Get AIPS user number */
	user = ObitSystemGetAIPSuser();
	/* set class */
	Aclass[1] = beamClass[0];
	/* allocate */
	cno = ObitAIPSDirAlloc(in->imDisk[i], user, Aname, Aclass, Atype, in->imSeq, &exist, err);
	/* Set info */
	ObitImageSetAIPS(image,OBIT_IO_byPlane,in->imDisk[i],cno,user,blc,trc,err);
	
      } else if (in->fileType==OBIT_IO_FITS) {
	/* set filename - derive from field - insure uniqueness */
	ct = in->imClass[1];
	if (in->imClass[1] != beamClass[0]) in->imClass[1] = beamClass[0];
	else in->imClass[1] = altbeamClass[0];
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
    } /* end Beam */
  }    /* end loop over images */
  
  /* Everything should now be initialized */
  in->nInit = in->numberImages;
 
} /* end  ObitImageMosaicWBSetFiles */

/**
 * Create an Image Mosaic based on a uv data and parameters attached thereto
 * 
 * \param name    Name to be given to new object
 * \param order  Spectral imaging order,0=flux,1=si, 2=curve
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
ObitImageMosaicWB* ObitImageMosaicWBCreate (gchar *name, olong order, ObitUV *uvData, 
					    ObitErr *err)
{
  ObitImageMosaicWB *out = NULL;
  ObitInfoType type;
  ObitIOType Type;
  ObitIOAccess access;
  gint32 i, nf, nif, nfif, dim[MAXINFOELEMDIM];
  olong Seq, *Disk=NULL, NField, nx[MAXFLD], ny[MAXFLD], catDisk, nDisk, numBeamTapes;
  olong overlap, imsize, fldsiz[MAXFLD], flqual[MAXFLD];
  ofloat FOV=0.0, xCells, yCells, RAShift[MAXFLD], DecShift[MAXFLD], Tapers[MAXFLD], 
    MaxBL, MaxW, Cells, BeamTapes[100], diam=24.5;
  ofloat *farr, Radius = 0.0, maxScale;
  ofloat shift[2] = {0.0,0.0}, cells[2], bmaj, bmin, bpa, beam[3];
  odouble ra0, dec0, refFreq1, refFreq2;
  gboolean doJ2B, doFull, doCalSelect, inFlysEye[MAXFLD], FacetNo[MAXFLD];
  ofloat equinox, minRad, ratio, OutlierFlux, OutlierDist, OutlierSI;
  olong  itemp, OutlierSize=0, nFlyEye = 0, *iarr=NULL;
  union ObitInfoListEquiv InfoReal; 
  gchar Catalog[100], Aname[100], Aclass[100];
  ObitImageMosaicClassInfo* mosaicClass;
  gchar *routine = "ObitImageMosaicWBCreate";

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitImageMosaicWBClassInit();

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
  
  /* Not bigger than FOV */
  imsize = MIN (imsize, ((2.0*3600.0*FOV/fabs(xCells))+10.99));
  overlap = 7;
  cells[0] = xCells; cells[1] = yCells;
  ObitUVGetRADec (uvData, &ra0, &dec0, err);
  if (err->error) Obit_traceback_val (err, routine, uvData->name, out);
  /* Copy nx to fldsiz - Make image sizes FFT friendly */
  for (i=0; i<NField; i++) fldsiz[i] = ObitFFTSuggestSize (MAX(32,nx[i]));
  minRad = 0.0;
  mosaicClass = (ObitImageMosaicClassInfo*)&myClassInfo;  /* Class pointer */
  if (FOV>0.0) minRad = mosaicClass->FlyEye(FOV, imsize, cells, overlap, shift, 
					    ra0, dec0, &NField, fldsiz, 
					    RAShift, DecShift, flqual, err); 
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

     mosaicClass->AddOutlier (Catalog, catDisk, minRad, cells, OutlierDist, 
			      OutlierFlux, OutlierSI, OutlierSize,
			      ra0, dec0, doJ2B, 
			      uvData->myDesc->crval[uvData->myDesc->jlocf], 
			      Radius, diam, &NField, fldsiz, RAShift, DecShift, flqual, 
			      err);
     if (err->error) Obit_traceback_val (err, routine, uvData->name, out);
  } /* end add outliers from catalog */

  /* Make sure some fields defined */
  Obit_retval_if_fail((NField>0), err, out, "%s: NO Fields defined", routine);

  /* Add any additional tapered images */
  mosaicClass->AddTapers (uvData, &numBeamTapes, BeamTapes,  fabs(cells[0]), Tapers,
			  &NField, fldsiz, RAShift, DecShift, flqual, inFlysEye, 
			  FacetNo, err);
  if (err->error) Obit_traceback_val (err, routine, uvData->name, out);
  
  /* Create output object */
  out = newObitImageMosaicWB (name, NField);

  /* Add taper array */
  out->BeamTapes = ObitMemAlloc0Name(numBeamTapes*sizeof(ofloat),"BeamTapes");

  /* Set maximum beam order */
  out->norder = order;

  /* Copy fldsiz to nx, ny 
     Double sizes of images and beams - need for image plane convolution */
  for (i=0; i<NField; i++) nx[i] = ny[i] = 2*fldsiz[i];
      
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
} /* end ObitImageMosaicWBCreate */

/**
 * Define images in an Image Mosaic from a UV data.
 * 
 * \param inn     The object to create images in,  Details are defined in members:
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
void ObitImageMosaicWBDefine (ObitImageMosaic *inn, ObitUV *uvData, gboolean doBeam,
			      ObitErr *err)
{
  olong i, nx, ny;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ofloat *farr=NULL, FOV;
  gboolean doCalSelect, *barr=NULL;
  ObitIOAccess access;
  ObitImage *tmpImage=NULL;
  ObitImageMosaicWB *in = (ObitImageMosaicWB*)inn;
  ObitImageMosaicClassInfo* mosaicClass;
  gchar *routine = "ObitImageMosaicWBDefine";

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
  ObitInfoListAlwaysPut (uvData->info, "xShift", OBIT_float, dim, in->RAShift);
  ObitInfoListAlwaysPut (uvData->info, "yShift", OBIT_float, dim, in->DecShift);
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
    /* Create regular ObitImage then convert to ObitImageWB */
    tmpImage =  ObitImageUtilCreateImage(uvData, i+1, doBeam, err);
    in->images[i] = (ObitImage*)ObitImageWBFromImage(tmpImage, in->norder, err);
    tmpImage = ObitImageUnref(tmpImage);
    /* If making a beam convert it as well */
    if (doBeam && in->images[i]->myBeam) {
     tmpImage =  (ObitImage*)ObitImageWBFromImage((ObitImage*)in->images[i]->myBeam, 
						  in->norder, err);
     in->images[i]->myBeam = ObitImageUnref((ObitImage*)in->images[i]->myBeam);
     in->images[i]->myBeam = (Obit*)tmpImage;  /* xfer ownership */
    }
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    /* For 2D imaging pick up shift */
    if (!in->images[i]->myDesc->do3D) {
      in->RAShift[i]  = in->images[i]->myDesc->xshift;
      in->DecShift[i] = in->images[i]->myDesc->yshift;
    }
    /* Add BeamTaper to Descriptor InfoList */
    dim[0] = dim[1] = 1;
    ObitInfoListAlwaysPut (in->images[i]->myDesc->info, "BeamTaper", OBIT_float, dim, 
			   &in->BeamTaper[i]);
  }    /* end loop over images */

 /* Create full field image if needed */
  if (in->doFull && (in->nInit<=0) && (in->numberImages>1)) { 
    /* Basic structure of field 1 */
    /* Create regular ObitImage then convert  to ObitImageWB */
    tmpImage =  ObitImageUtilCreateImage(uvData, 1, FALSE, err);
    in->FullField = (ObitImage*)ObitImageWBFromImage(tmpImage, in->norder, err);
    tmpImage = ObitImageUnref(tmpImage);
    /* Replace name */
    g_free(in->FullField->name);
    in->FullField->name = g_strdup("Flattened Image");
    mosaicClass = (ObitImageMosaicClassInfo*)&myClassInfo;  /* Class pointer */
    FOV = mosaicClass->ObitImageMosaicFOV (inn, err);
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

  ObitImageMosaicWBSetFiles (inn, doBeam, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
 
} /* end ObitImageMosaicWBDefine */


/**
 * Add a Field to a mosaic
 * Wideband version: image doubles in size and converted to ObitImageWB
 * \param inn       The object with images to modify
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
void ObitImageMosaicWBAddField (ObitImageMosaic *inn, ObitUV *uvData, 
				olong nx, olong ny, olong nplane, 
				ofloat RAShift, ofloat DecShift,
				gboolean isAuto, ObitErr *err)
{
  ofloat *ftemp;
  olong     i, *itemp;
  gboolean  *btemp;
  ObitImage **imtemp;
  ObitImageMosaicWB *in = (ObitImageMosaicWB*)inn;
  gchar *routine = "ObitImageMosaicWBAddField";

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Expand arrays */
  in->nInit = in->numberImages;  /* To ensure new field gets initialized */
  in->numberImages++;

  /* Image array */
  imtemp = ObitMemAlloc0Name(in->numberImages*sizeof(ObitImage*),"ImageMosaicWB images");
  for (i=0; i<in->nInit; i++) imtemp[i] = in->images[i]; imtemp[i] = NULL;
  in->images = ObitMemFree(in->images);
  in->images = imtemp;

  /* Disk array - add new images on same disk as first */
  itemp  = ObitMemAlloc0Name(in->numberImages*sizeof(olong),"ImageMosaicWB imDisk");
  for (i=in->nInit; i<in->numberImages; i++) itemp[i] = in->imDisk[0]; 
  for (i=0; i<in->nInit; i++) itemp[i] = in->imDisk[i]; 
  in->imDisk = ObitMemFree(in->imDisk);
  in->imDisk = itemp;
 
  /* isAuto array */
  itemp  = ObitMemAlloc0Name(in->numberImages*sizeof(olong),"ImageMosaicWB isAuto");
  for (i=0; i<in->nInit; i++) itemp[i] = in->isAuto[i]; 
  for (i=in->nInit; i<in->numberImages; i++) itemp[i] = -1; 
  /* Large dummy value in new element if isAuto */
  if (isAuto) itemp[in->numberImages-1] = in->numberImages+9999999;
  in->isAuto = ObitMemFree(in->isAuto);
  in->isAuto = itemp;
 
  /* isShift array */
  itemp  = ObitMemAlloc0Name(in->numberImages*sizeof(olong),"ImageMosaicWB isShift");
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
  itemp  = ObitMemAlloc0Name(in->numberImages*sizeof(olong),"ImageMosaicWB nx");
  for (i=0; i<in->nInit; i++) itemp[i] = in->nx[i]; 
  /* Double but at least 512 */
  nx = MAX (2*nx, 512);
  itemp[i] = ObitFFTSuggestSize (nx);
  in->nx = ObitMemFree(in->nx);
  in->nx = itemp;
  itemp  = ObitMemAlloc0Name(in->numberImages*sizeof(olong),"ImageMosaicWB ny");
  for (i=0; i<in->nInit; i++) itemp[i] = in->ny[i]; 
  /* Double but at least 512 */
  ny = MAX (2*ny, 512);
  itemp[i] = ObitFFTSuggestSize (ny);
  in->ny = ObitMemFree(in->ny);
  in->ny = itemp;
  itemp  = ObitMemAlloc0Name(in->numberImages*sizeof(olong),"ImageMosaicWB nplane");
  for (i=0; i<in->nInit; i++) itemp[i] = in->nplane[i]; 
  itemp[i] = nplane;
  in->nplane = ObitMemFree(in->nplane);
  in->nplane = itemp;

  /* Shift */
  ftemp  = ObitMemAlloc0Name(in->numberImages*sizeof(ofloat),"ImageMosaicWB RAShift");
  for (i=0; i<in->nInit; i++) ftemp[i] = in->RAShift[i]; 
  ftemp[i] = RAShift;
  in->RAShift = ObitMemFree(in->RAShift);
  in->RAShift = ftemp;
  ftemp  = ObitMemAlloc0Name(in->numberImages*sizeof(ofloat),"ImageMosaicWB DecShift");
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
  ObitImageMosaicWBDefine (inn, uvData, TRUE, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Create/initialize image */
  ObitImageMosaicWBSetFiles (inn, TRUE, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

} /* end ObitImageMosaicWBAddField */

/**
 * Project the tiles of a Mosaic to the full field flattened image
 * Ignore facets with BeamTaper>0
 * Routine translated from the AIPSish 4MASS/SUB/FLATEN.FOR/FLATEN  
 * Wideband imaging version, flatten all spectral planes.
 * \param inn     The object with images
 * \param err     Error stack, returns if not empty.
 */
void ObitImageMosaicWBFlatten (ObitImageMosaic *inn, ObitErr *err)
{
  olong   blc[IM_MAXDIM]={1,1,1,1,1}, trc[IM_MAXDIM]={0,0,0,0,0};
  olong nterm, i, j, radius, rad, plane[IM_MAXDIM] = {1,1,1,1,1}, hwidth = 2;
  olong *naxis, pos1[IM_MAXDIM], pos2[IM_MAXDIM];
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1}, iplane;
  ObitImage *out=NULL, *tout1=NULL, *tout2=NULL;
  ObitFArray *sc2=NULL, *sc1=NULL;
  ObitIOSize IOsize = OBIT_IO_byPlane;
  ofloat xpos1[IM_MAXDIM], xpos2[IM_MAXDIM];
  gboolean overlap;
  ObitImageMosaicWB *in = (ObitImageMosaicWB*)inn;
  gchar *routine = "ObitImageMosaicWBFlatten";

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Must have output image defined */
  if (in->FullField==NULL) {
    Obit_log_error(err, OBIT_Error, "No flattened image defined");
    return;
  }

  /* Tell user */
  Obit_log_error(err, OBIT_InfoErr, "Flattening flys eye to single image");

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
  radius = in->Radius;

  /* Loop over spectral planes */
  if (ObitImageWBIsA(in->FullField))
    nterm = MAX(1,(((ObitImageWB*)in->FullField)->order+1));
  else
    nterm = 1;
  for (iplane=0; iplane<nterm; iplane++) {
    
    /* Zero fill accumulations */
    ObitFArrayFill (sc1, 0.0);
    ObitFArrayFill (sc2, 0.0);
    
    /* Set plane */
    for (i=0; i<IM_MAXDIM; i++) blc[i] = 1; blc[2] = iplane+1;
    for (i=0; i<IM_MAXDIM; i++) trc[i] = 0; trc[2] = iplane+1;
    
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
  } /* end loop over planes */

  /* Cleanup */
  ObitImageClose (in->FullField, err);
  in->FullField->image = ObitFArrayUnref(in->FullField->image); /* Free buffer */
  sc2   = ObitFArrayUnref(sc2);
  tout1 = ObitImageUnref(tout1);
  tout2 = ObitImageUnref(tout2);

} /* end ObitImageMosaicWBFlatten */

/**
 * Make a single field ImageMosaicWB corresponding to the field with
 * the highest summed peak of CCs if any above MinFlux.
 * Fields with 1-rel numbers in the zero terminated list ignore are ignored
 * \param inn      ImageMosaicWB to process
 * \param MinFlux  Min. flux density for operation.
 * \param ignore   0 terminated list of 1-rel field numbers to ignore
 * \param field    [out] the 1-rel number of the field copied
 * \param err      Error/message stack
 * return Newly created ImageMosaicWB or NULL
 */
ObitImageMosaicWB* ObitImageMosaicWBMaxField (ObitImageMosaic *inn, 
					      ofloat MinFlux, olong *ignore, olong *field,
					      ObitErr* err) 
{
  ObitImageMosaicWB* out=NULL;
  ObitTableCC *CCTab=NULL, *inCCTable=NULL, *outCCTable=NULL;
  ObitImage   *tmpImage=NULL, *tmpBeam=NULL;
  gint32 dim[MAXINFOELEMDIM];
  ObitInfoType type;
  olong   i, nfield, ifield, itemp, noParms, nccpos;
  olong   maxField;
  olong  CCVer;
  ofloat tmax, xcenter, ycenter, xoff, yoff, radius, cells[2], maxCC;
  gboolean forget;
  ObitImageMosaicWB *mosaic = (ObitImageMosaicWB*)inn;
  gchar *routine = "ObitImageMosaicWBMaxField";

  /* Error checks */
  if (err->error) return out;  /* previous error? */
  g_assert(ObitImageMosaicWBIsA(mosaic));

  /* Number of fields */
  nfield = mosaic->numberImages;

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
    if (forget) break;

    /* Open image in case header needs update */
    ObitImageOpen (mosaic->images[ifield], OBIT_IO_ReadWrite, err);
    if  (err->error) Obit_traceback_val (err, routine, mosaic->images[ifield]->name, out);

    /* Make CC table object */
    noParms = 0;
    CCTab = newObitTableCCValue ("Temp CC", (ObitData*)mosaic->images[ifield],
				 &CCVer, OBIT_IO_ReadWrite, noParms, err);
    if  (err->error) Obit_traceback_val (err, routine, mosaic->images[ifield]->name, out);

    /* Determine maximum */
    nccpos = CCTab->myDesc->nrow;
    ObitImageMosaicMaxCC (CCTab, nccpos, radius, &tmax, &xcenter, &ycenter, &xoff, &yoff, err);
    if  (err->error) Obit_traceback_val (err, routine, mosaic->images[ifield]->name, out);
    
    /* this one of interest? */
    if ((tmax > MinFlux) && (tmax > maxCC)) {
      maxCC = tmax;
      maxField = ifield;
    }

    /* Delete temporary table */
    CCTab = ObitTableCCUnref(CCTab);

    /* Close/update image */
    ObitImageClose(mosaic->images[ifield], err);
    if  (err->error) Obit_traceback_val (err, routine, mosaic->images[ifield]->name, out);
  } /* end loop  L500: */

  /* Catch anything? */
  if (maxField<0) return out;

  /* Make output ImageMosaicWB */
  out = newObitImageMosaicWB ("Temp mosaic", 1);
  tmpImage = newObitImageScratch (mosaic->images[maxField], err);
  /* Copy image */
  tmpImage = (ObitImage*)ObitImageWBCopy((ObitImageWB*)mosaic->images[maxField], (ObitImageWB*)tmpImage, err);
  ObitImageMosaicSetImage ((ObitImageMosaic*)out, 0, tmpImage, err);
  if  (err->error) Obit_traceback_val (err, routine, mosaic->images[maxField]->name, out);
  /* Give more sensible name */
  if (tmpImage->name) g_free(tmpImage->name);
  tmpImage->name = g_strdup("Peel Image");

  /* Copy beam NEED to copy higher order beams  */
  tmpBeam = newObitImageScratch ((ObitImage*)mosaic->images[maxField]->myBeam, err);
  tmpBeam = (ObitImage*)ObitImageWBCopy((ObitImageWB*)mosaic->images[maxField]->myBeam, (ObitImageWB*)tmpBeam, err);
  ObitImageMosaicSetImage ((ObitImageMosaic*)out, 0, tmpImage, err);
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

  /* Copy Clean components */
  CCVer = 1;
  noParms = 0;
  inCCTable = newObitTableCCValue ("Peeled CC", (ObitData*)mosaic->images[maxField],
				   &CCVer, OBIT_IO_ReadOnly, noParms, 
				   err);
  CCVer = 1;
  outCCTable = newObitTableCCValue ("outCC", (ObitData*)out->images[0],
				    &CCVer, OBIT_IO_ReadWrite, inCCTable->noParms, 
				    err);
  outCCTable  = ObitTableCCCopy (inCCTable, outCCTable, err);
  inCCTable   = ObitTableCCUnref(inCCTable);
  outCCTable  = ObitTableCCUnref(outCCTable);

  *field = maxField+1;  /* Which field is this? */
  return out;
} /* end of routine ObitImageMosaicWBMaxField */ 

/**
 * Convert structure information to entries in an ObitInfoList
 * \param inn     Object of interest.
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
 *      \li "xxxnorder"       olong maximum beam order
 * \param err     ObitErr for reporting errors.
 */
void ObitImageMosaicWBGetInfo (ObitImageMosaic *inn, gchar *prefix, ObitInfoList *outList, 
			       ObitErr *err)
{ 
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *keyword=NULL;
  ObitImageMosaicWB *in = (ObitImageMosaicWB*)inn;
  const ObitImageMosaicClassInfo *ParentClass = myClassInfo.ParentClass;
  gchar *routine = "ObitImageMosaicWBGetInfo";

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Get basic info from parent class */
  if (ParentClass->ObitImageMosaicGetInfo!=ObitImageMosaicWBGetInfo)
    ParentClass->ObitImageMosaicGetInfo (inn, prefix, outList, err);
  if  (err->error) Obit_traceback_msg (err, routine, in->name);

  /* "xxxnorder"       olong maximum beam order */
  if (prefix) keyword = g_strconcat (prefix, "norder", NULL);
  else        keyword = g_strdup("norder");
  dim[0] = 1; dim[1] = 1;
  ObitInfoListAlwaysPut(outList, keyword, OBIT_long, dim, &in->norder);
  g_free(keyword);

} /* end ObitImageMosaicWBGetInfo */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitImageMosaicWBClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitImageMosaicWBClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitImageMosaicWBClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitImageMosaicWBClassInfoDefFn (gpointer inClass)
{
  ObitImageMosaicWBClassInfo *theClass = (ObitImageMosaicWBClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitImageMosaicWBClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitImageMosaicWBClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitImageMosaicWBGetClass;
  theClass->newObit       = (newObitFP)newObitImageMosaicWB;
  theClass->ObitCopy      = (ObitCopyFP)ObitImageMosaicWBCopy;
  theClass->ObitClone     = NULL;  /* Different call */
  theClass->ObitClear     = (ObitClearFP)ObitImageMosaicWBClear;
  theClass->ObitInit      = (ObitInitFP)ObitImageMosaicWBInit;
  theClass->ObitImageMosaicFromInfo = (ObitImageMosaicFromInfoFP)ObitImageMosaicWBFromInfo;
  theClass->ObitImageMosaicGetInfo  = (ObitImageMosaicGetInfoFP)ObitImageMosaicWBGetInfo;
  theClass->ObitImageMosaicZapImage =
    (ObitImageMosaicZapImageFP)ObitImageMosaicWBZapImage;
  theClass->ObitImageMosaicSetFiles = 
    (ObitImageMosaicSetFilesFP)ObitImageMosaicWBSetFiles;
  theClass->ObitImageMosaicCreate = 
    (ObitImageMosaicCreateFP)ObitImageMosaicWBCreate;
  theClass->ObitImageMosaicDefine = 
    (ObitImageMosaicDefineFP)ObitImageMosaicWBDefine;
  theClass->ObitImageMosaicFlatten = 
    (ObitImageMosaicFlattenFP)ObitImageMosaicWBFlatten;
  theClass->ObitImageMosaicAddField = 
    (ObitImageMosaicAddFieldFP)ObitImageMosaicWBAddField;
  theClass->ObitImageMosaicMaxField = 
    (ObitImageMosaicMaxFieldFP)ObitImageMosaicWBMaxField;

} /* end ObitImageMosaicWBClassDefFn */


/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitImageMosaicWBInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitImageMosaicWB *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */

} /* end ObitImageMosaicWBInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitImageMosaicWB* cast to an Obit*.
 */
void ObitImageMosaicWBClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitImageMosaicWB *in = inn;
  
  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  
  /* delete this class members */
 
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitImageMosaicWBClear */


