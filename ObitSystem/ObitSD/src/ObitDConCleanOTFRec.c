/* $Id$ */
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

#include "ObitDConCleanOTFRec.h"
#include "ObitMem.h"
#include "ObitTableCCUtil.h"
#include "ObitOTFUtil.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitDConCleanOTFRec.c
 * ObitDConCleanOTFRec class function definitions.
 * Record based OTF CLEAN class.
 * This class is derived from the ObitDCon class.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitDConCleanOTFRec";

/** Function to obtain parent ClassInfo - ObitDConClean */
static ObitGetClassFP ObitParentGetClass = ObitDConCleanGetClass;

/**
 * ClassInfo structure ObitDConCleanOTFRecClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitDConCleanOTFRecClassInfo myClassInfo = {FALSE};

/*--------------- File Global Variables  ----------------*/


/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitDConCleanOTFRecInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitDConCleanOTFRecClear (gpointer in);

/** Read Beam into Beam patch */
static void ReadBP (ObitDConCleanOTFRec* in, ObitErr *err);

/** Private: Set Class function pointers. */
static void ObitDConCleanOTFRecClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitDConCleanOTFRec* newObitDConCleanOTFRec (gchar* name)
{
  ObitDConCleanOTFRec* out;
  /*gchar *routine = "newObitDConCleanOTFRec";*/

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitDConCleanOTFRecClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitDConCleanOTFRec));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitDConCleanOTFRecInit((gpointer)out);

 return out;
} /* end newObitDConCleanOTFRec */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitDConCleanOTFRecGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitDConCleanOTFRecClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitDConCleanOTFRecGetClass */

/**
 * Make a deep copy of an ObitDConCleanOTFRec.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitDConCleanOTFRec* ObitDConCleanOTFRecCopy  (ObitDConCleanOTFRec *in, 
					       ObitDConCleanOTFRec *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  gboolean oldExist;
  gchar *outName;
  gchar *routine = "ObitDConCleanOTFRecCopy";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return out;
  g_assert (ObitIsA(in, &myClassInfo));
  if (out) g_assert (ObitIsA(out, &myClassInfo));

  /* Create if it doesn't exist */
  oldExist = out!=NULL;
  if (!oldExist) {
    /* derive object name */
    outName = g_strconcat ("Copy: ",in->name, NULL);
    out = newObitDConCleanOTFRec(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, out);

  /*  copy this class */
  out->beam  = ObitImageUnref(out->beam);
  if (in->beam) out->beam  = ObitImageRef(in->beam);
  out->clean  = ObitImageUnref(out->clean);
  out->clean  = ObitImageRef(in->clean);
  out->myOTF  = ObitOTFUnref(out->myOTF);
  out->myOTF  = ObitOTFRef(in->myOTF);
  out->scrOTF = ObitOTFUnref(out->scrOTF);
  out->scrOTF = ObitOTFRef(in->scrOTF);
  out->display= ObitDisplayUnref(out->display);
  out->display= ObitDisplayRef(in->display);
  out->doRestore= in->doRestore;
  out->nCCSub   = in->nCCSub;
  out->fracPeak = in->fracPeak;
  out->cleanSize= in->cleanSize;

  return out;
} /* end ObitDConCleanOTFRecCopy */

/**
 * Make a copy of a object but do not copy the actual data
 * This is useful to create an DConCleanOTFRec similar to the input one.
 * \param in  The object to copy
 * \param out An existing object pointer for output, must be defined.
 * \param err Obit error stack object.
 */
void ObitDConCleanOTFRecClone  (ObitDConCleanOTFRec *in, ObitDConCleanOTFRec *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  gchar *routine = "ObitDConCleanOTFRecClone";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitIsA(out, &myClassInfo));

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /*  copy this class */
  out->beam  = ObitImageUnref(out->beam);
  if (in->beam) out->beam  = ObitImageRef(in->beam);
  out->clean = ObitImageUnref(out->clean);
  out->clean = ObitImageRef(in->clean);
  out->myOTF  = ObitOTFUnref(out->myOTF);
  out->myOTF  = ObitOTFRef(in->myOTF);
  out->scrOTF = ObitOTFUnref(out->scrOTF);
  out->scrOTF = ObitOTFRef(in->scrOTF);
  out->display= ObitDisplayUnref(out->display);
  out->display= ObitDisplayRef(in->display);
  out->doRestore= in->doRestore;
  out->nCCSub   = in->nCCSub;
  out->fracPeak = in->fracPeak;
  out->cleanSize= in->cleanSize;

} /* end ObitDConCleanOTFRecClone */

/**
 * Creates an ObitDConCleanOTFRec 
 * \param name   An optional name for the object.
 * \param inOTF  Data to image
 * \li "outName" Obit_string       = Base name of image files, 
 *                                   dirty="Dirty"+outName+".fits"
 *                                   beam = "Beam"+outname+".fits"
 *                                   clean = outname+".fits"
 * \li "outDisk" OBIT_long [1,1,1]  = Disk number for output files. [def 1]
 * \li "Scans"   OBIT_long [2,1,1]  = Range of scan numbers, [def=all]
 * \li "doCalib" OBIT_long (1,1,1)  = >0 -> calibrate,
 * \li "keepCal" OBIT_bool(1,1,1)  = True = keep cal-on data [def TRUE]
 * \li "gainUse" OBIT_long (1,1,1)  = SN/CL table version number, 0-> use highest
 * \li "flagVer" OBIT_long (1,1,1)  = Flag table version, 0-> use highest, <0-> none
 * \li "Stokes"  OBIT_string (4,1,1)= Selected output Stokes parameters:
 *                                   "I", "V", " " -> "I" [def "I"]
 * \li "BChan"   OBIT_long (1,1,1)  = First spectral channel selected. [def all]
 * \li "EChan"   OBIT_long (1,1,1)  = Highest spectral channel selected. [def all]
 * \li "Feeds"   OBIT_long (?,1,1)  = List of selected feed numbers, [def all.]
 * \li "Targets" Obit_string [?,?] = List of target names to include
 * \li "RA"      OBIT_float scalar = Center RA in deg
 * \li "Dec"     OBIT_float scalar = Center Dec in deg
 * \li "nx"      OBIT_long scalar   = Number of cells in RA
 * \li "ny"      OBIT_long scalar   = Number of cells in Declination
 * \li "beamNx"  OBIT_long scalar   = Number of "x" pixels [def 32]
 * \li "beamNy"  OBIT_long scalar   = Number of "y" pixels [def 32]
 * \li "xCells"  OBIT_float scalar = Cell spacing (deg) in RA
 * \li "yCells"  OBIT_float scalar = Cell spacing (deg) in Dec
 * \li "minWt"   OBIT_float scalar = Minimum sum of weights
 *                                   Default = 0.1.
 * \li "Proj"    OBIT_string (4,1,1) Projection string "-SIN", "-ARC", "-TAN"
 *                                   [Default "-SIN"]
 * \li "ConvType"OBIT_long scalar = Convolving function type: [def=3]
 *                  0 = pillbox, 3 = Gaussian, 4 = Exp*Sinc, 5 = Spherodial wave
 * \li "ConvParm"OBIT_float[10] = Convolving function parameters (see below)
 * 
 * Gridding convolution functions:
 * \li 0 = pillbox, 
 * \li 2 = Sinc, 
 *    Parm[0] = halfwidth in cells,
 *    Parm[1] = Expansion factor
 * \li 3 = Gaussian,
 *    Parm[0] = halfwidth in cells,[def 3.0]
 *    Parm[1] = Gaussian with as fraction or raw beam [def 1.0]
 * \li 4 = Exp*Sinc
 *    Parm[0] = halfwidth in cells, [def 2.0]
 *    Parm[1] = 1/sinc factor (cells) [def 1.55]
 *    Parm[2] = 1/exp factor (cells) [def 2.52]
 *    Parm[3] = exp power [def 2.0]
 * \li 5 = Spherodial wave
 *    Parm[0] = halfwidth in cells [def 3.0]
 *    Parm[1] = Alpha [def 5.0]
 *    Parm[2] = Expansion factor [not used]
 * \param err    Obit error stack object.
 * \return the new object.
 */
ObitDConCleanOTFRec* 
ObitDConCleanOTFRecCreate (gchar* name, ObitOTF *inOTF, ObitErr *err)
{
  olong nfield;
  ObitDConCleanOTFRec* out=NULL;
  /*gchar *routine = "ObitDConCleanOTFRecCreate";*/

  /* error checks */
  if (err->error) return out;
  g_assert (ObitOTFIsA(inOTF));

  /* Create basic structure */
  out = newObitDConCleanOTFRec (name);

  /* save inputs */
  out->myOTF   = ObitOTFRef(inOTF);

  /* Arrays - including those in parent classes */
  nfield =  1;
  out->gain    = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean Loop gain");
  out->minFlux = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean minFlux");
  out->factor  = ObitMemAlloc0Name(nfield*sizeof(ofloat),"Clean factor");

  return out;
} /* end ObitDConCleanOTFRecCreate */

/**
 * Do deconvolution, uses function on class pointer
 * Does final flatten if FullField member of mosaic member is defined.
 * CLEAN control parameters are in the ObitInfoList member:
 * \li "fracPeak"OBIT_float scalar = Fraction of residual to CLEAN to 
 *                                   before major cycle  [def 0.75]
 * \li "Niter"   OBIT_long scalar   = Maximum number of CLEAN iterations
 * \li "maxPixel" OBIT_long scalar  = Maximum number of residuals [def 20000]
 * \li "Patch"    OBIT_long scalar  = Minimum beam patch in pixels [def 50]
 * \li "BMaj"    OBIT_float scalar = Restoring beam major axis (deg)
 * \li "BMin"    OBIT_float scalar = Restoring beam minor axis (deg)
 * \li "BPA"     OBIT_float scalar = Restoring beam position angle (deg)
 * \li "CCVer"   OBIT_long array  = CLEAN table version per field
 * \li "Gain"    OBIT_float array  = CLEAN loop gain per field
 * \li "minFlux" OBIT_float array  = Minimun flux density (Jy)  per field
 * \li "Factor"  OBIT_float array  = CLEAN depth factor per field
 * \li "Plane"   OBIT_long array    = Plane being processed, 1-rel indices of axes 3-?
 * \li "dispURL" OBIT_string scalar= URL of display server
 * \li "ConvType" OBIT_long scalar = Convolving function type: [def=3]
 *                0 = pillbox, 3 = Gaussian, 4 = Exp*Sinc, 5 = Spherodial wave
 * \li "ConvParm" OBIT_float[10] = Convolving function parameters
 * \li "deBias"   OBIT_bool scalar = Subtract calibration bias from image? [def False]
 *                Note, this doesn't really work the way you would like
 * \li "deMode"   OBIT_bool scalar = Subtract image mode from image? [def False]
 * \li "minWt"    OBIT_float (1,1,1) Minimum summed gridding convolution weight 
 *                as a fraction of the maximum [def 0.01]
 * \param in   The object to deconvolve
 * \param err Obit error stack object.
 */
void ObitDConCleanOTFRecDeconvolve (ObitDCon *inn, ObitErr *err)
{
  ObitDConCleanOTFRec *in;
  gboolean done;
  olong disk;
  olong jtemp;
  gboolean quit, noClean;
  olong numPix;
  ofloat peakFlux;
    olong maxPix;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong blc[IM_MAXDIM] = {1,1,1,1,1};
  olong trc[IM_MAXDIM] = {0,0,0,0,0};
  gchar *tname, filename[129];
  gchar        *imgParms[] = {  /* Imaging parameters to copy to scratch OTF */
    "RA", "Dec", "nx", "ny" , "beamNx", "beamNy", "xCells" , "yCells", "minWt",
    "Proj", "ConvType", "ConvParm", "deMode", "deBias",
    NULL
  };
  const ObitDConCleanOTFRecClassInfo *inClass;
  gchar *routine = "ObitDConCleanOTFRecDeconvolve";

  /* Cast input to this type */
  in = (ObitDConCleanOTFRec*)inn;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitDConCleanOTFRecIsA(in));

  inClass = (ObitDConCleanOTFRecClassInfo*)in->ClassInfo; /* class structure */

  /* Open inOTF with selection to fully define */
  ObitOTFOpen (in->myOTF, OBIT_IO_ReadCal, err);
  ObitOTFClose (in->myOTF, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Initialize Clean object - Create scratch OTF */
  if (!in->scrOTF) in->scrOTF = newObitOTFScratch (in->myOTF, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Copy OTF to scratch */
  in->scrOTF = ObitOTFCopy (in->myOTF, in->scrOTF, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Copy imaging parameters */
  ObitInfoListCopyList (in->myOTF->info, in->scrOTF ->info, imgParms);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
 
  /* Create Image */
  if (!in->clean) {
    in->clean = ObitOTFUtilCreateImage (in->scrOTF, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
  }
 
  /* define output */
  disk = 1;
  ObitInfoListGetTest(in->myOTF->info, "outDisk", &type, dim, &disk);
  strncpy (filename, "Unnamed",128);
  ObitInfoListGetTest(in->myOTF->info, "outName", &type, dim, filename);
  tname = g_strconcat (filename, ".fits", NULL);
  ObitImageSetFITS(in->clean,OBIT_IO_byPlane, disk, tname, blc, trc, err);
  g_free(tname);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Gridding weight image */
  if (!in->weight) {
    in->weight = ObitOTFUtilCreateImage (in->scrOTF, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    tname = g_strconcat (filename, ".wt", NULL);
    ObitImageSetFITS(in->weight,OBIT_IO_byPlane, disk, tname, blc, trc, err);
    g_free(tname);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
  }
  
  /* Image mosaic for window */
  if (!in->mosaic) {
    in->mosaic = newObitImageMosaic (in->name, 1);
    ObitImageMosaicSetImage (in->mosaic, 0, in->clean, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
  }

  /* Window object */
  if (!in->window) {
    in->window = ObitDConCleanWindowCreate ("CleanWindow", in->mosaic, err);
    ObitDConCleanDefWindow((ObitDConClean*)in, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
  }
  
  /* Form initial image */
  ObitOTFUtilMakeImage (in->scrOTF,  in->clean, TRUE, NULL, in->weight, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  
  /* Get beam */
  if (!in->beam) in->beam = ObitImageRef(in->clean->myBeam);
  
  /* Get parameters */
  inClass->ObitDConGetParms(inn, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* How many pixels selected? */
  numPix = ObitDConCleanWindowCount (in->window, 1, err);
  numPix += 10; /* for good measure */
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* (Re)Create Pixel List if needed */
  if ((in->Pixels) && (numPix>in->Pixels->maxPixel))  
    in->Pixels = ObitDConCleanPxListUnref(in->Pixels);
  if (!in->Pixels) {
    in->Pixels = ObitDConCleanPxListCreate("Pixel List", in->mosaic, 
					   numPix, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
  }

  /* Copy control info to PixelList */
  ObitInfoListCopyData(in->info, in->Pixels->info);
  dim[0] = dim[1] = 1;
  ObitInfoListAlwaysPut(in->info, "Patch", OBIT_long, dim, &in->beamPatchSize);
  maxPix = numPix;
  ObitInfoListAlwaysPut(in->info, "maxPixel", OBIT_long, dim, &maxPix);

  /* Reset/Init Pixel list*/
  ObitDConCleanPxListReset (in->Pixels, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Read Beam patch if needed*/
  if (!in->BeamPatch) ReadBP (in, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Reset Sky Model */
  in->nCCSub = 0;
  
  /* Save actual CC version if not specified */
  if (in->CCver<=0) {
    if (in->Pixels->CCver[0]>0) in->CCver = in->Pixels->CCver[0];
    jtemp = in->CCver;
    ObitInfoListAlwaysPut(in->info, "CCVer", OBIT_long, dim, &jtemp);
  }

  /* Save target minimun flux density as minFlux used to control major cycles */
  in->totalMinFlux = in->minFlux[0];
  peakFlux = fabs(in->clean->myDesc->maxval);  /* Peak in image */

  /* If no CLEANing requested then we're done */
  done = FALSE;
  noClean = in->niter<= 0;
  if (noClean) done = TRUE;

  /* Loop until Deconvolution done */
  /* Display/edit windows if enabled */
  if (in->display) {
    if (noClean)
      quit = ObitDisplayShow (in->display, (Obit*)in->mosaic, NULL, 1, err);
    else
      quit = ObitDisplayShow (in->display, (Obit*)in->mosaic, in->window, 
			      1, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    if (quit) {done=TRUE;}
  }

  while (!done) {

    /* Set minFlux */
    in->minFlux[0] = in->fracPeak * peakFlux;

    /* Get image/beam statistics needed for this cycle */
    inClass->ObitDConCleanPixelStats((ObitDConClean*)in, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);

    /* Pick components */
    done = inClass->ObitDConCleanSelect((ObitDConClean*)in, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
  
    peakFlux = in->Pixels->maxResid;  /* Peak residual in window */

    ObitErrLog(err);  /* Progress Report */

    /* Subtract model, remake residual image */
    inClass->ObitDConCleanSub((ObitDConClean*)in, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    
    /* Display/edit windows if enabled */
    if (in->display) {
      quit = ObitDisplayShow (in->display, (Obit*)in->mosaic, in->window, 
			      1, err);
      if (err->error) Obit_traceback_msg (err, routine, in->name);
      if (quit) {done=TRUE; break;}
    }
  } /* end clean loop */

  /* Restore if needed */
  if (in->doRestore && !noClean) {
    inClass->ObitDConCleanRestore((ObitDConClean*)in, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    
    /* Display restored image if enabled */
    if (in->display) {
      quit = ObitDisplayShow (in->display, (Obit*)in->mosaic, NULL, 
			      1, err);
      if (err->error) Obit_traceback_msg (err, routine, in->name);
    }
  } /* end if restore */
} /* end ObitDConCleanOTFRecDeconvolve */

/**
 * Read any base class parameters and then
 * read CLEAN control parameters from the ObitInfoList member:
 * From Parent classes:
 * \li "Niter"   OBIT_long scalar   = Maximum number of CLEAN iterations
 * \li "CCVer"   OBIT_long array    = CLEAN table version
 * \li "Gain"    OBIT_float array  = CLEAN loop gain
 * \li "minFlux" OBIT_float array  = Minimum flux density (Jy)
 * \li "Factor"  OBIT_float array  = CLEAN depth factor
 * \li "Plane"   OBIT_long array    = Plane being processed, 1-rel indices of axes 3-?
 * \li "doRestore" OBIT_bool       = Restore image when done? [def TRUE]
 * 
 * This Class:
 * \li "Patch"   OBIT_long scalar   = Beam patch in pixels [def 100]
 * \li "BeamSize"OBIT_float scalar = Restoring beam FWHM (deg)
 * \li "fracPeak"OBIT_float scalar = Fraction of residual to CLEAN to before 
 *                                   major cycle [def 0.75]
 *
 * The OTF data file should have the following imaging parameters:
 * \li "outName" Obit_string       = Base name of image files, 
 *                                   dirty="Dirty"+outName+".fits"
 *                                   beam = "Beam"+outname+".fits"
 *                                   clean = outname+".fits"
 * \li "outDisk" OBIT_long [1,1,1]  = Disk number for output files. [def 1]
 * \param in  The CLEAN object as base class
 * \param err Obit error stack object.
 */
void  ObitDConCleanOTFRecGetParms (ObitDCon *inn, ObitErr *err)
{
  ObitDConCleanOTFRec *in = (ObitDConCleanOTFRec*)inn;  /* as this class */
  ObitDConClassInfo *ParentClass;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  gchar *dispURL=NULL, tname[129];
  olong i;
  /*  olong itemp;
      union ObitInfoListEquiv InfoReal; */
  gchar *routine = "ObitDConCleanOTFRecGetParms";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Read any parent class parameters */
  ParentClass = (ObitDConClassInfo*)myClassInfo.ParentClass;
  ParentClass->ObitDConGetParms(inn, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Use Raw Beam size */
  in->cleanSize =in->myOTF->myDesc->beamSize ;
  ObitInfoListGetTest(in->info, "BeamSize", &type, dim, &in->cleanSize);
  if (in->cleanSize<=0.0) in->cleanSize =in->myOTF->myDesc->beamSize;

  /* beam patch */
  in->beamPatchSize = 100;
  ObitInfoListGetTest(in->info, "Patch", &type, dim, &in->beamPatchSize);

  /* Fraction of peak per major cycle */
  in->fracPeak = 0.75;
  ObitInfoListGetTest(in->info, "fracPeak", &type, dim, &in->fracPeak);

  /* Restore image when done? */
  in->doRestore = TRUE;
  ObitInfoListGetTest(in->info, "doRestore", &type, dim, &in->doRestore);

  /* Image display? */
  if (!in->display) {
    ObitInfoListGetP(in->info, "dispURL", &type, dim, (gpointer)&dispURL);
    /* dispURL not NULL terminated */
    if (dispURL) {for (i=0; i<dim[0]; i++) tname[i] = dispURL[i]; tname[i]=0;}
    if (dispURL && (strncmp(tname, "None", 4))) 
      in->display = ObitDisplayCreate("Display", tname, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
  }

} /* end ObitDConCleanOTFRecGetParms */

/**
 * Subtract new components from scratch OTF and remake residual image
 * \param in   The object to deconvolve
 * \param err Obit error stack object.
 */
void ObitDConCleanOTFRecSub(ObitDConClean *inn, ObitErr *err)
{
  ObitDConCleanOTFRec *in;
  ObitTableCC *CCTable=NULL;
  ObitTableCCRow *row=NULL;
  olong pos[2], beamCen[2], ndim, loop, nrow, iRow;
  ofloat FWHM, fmax, saveBeam, bmFact;
  ObitFArray *model=NULL, *RawBeam=NULL;
  gchar *tname;
  gchar *routine = "ObitDConCleanOTFRecSub";

  /* Cast input to this type */
  in = (ObitDConCleanOTFRec*)inn;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitDConCleanOTFRecIsA(in));

  /* Create FArray with components to be subtracted convolved with
     the telescope beam */
  ndim = 2;
  model = ObitFArrayCreate("Model", ndim, in->clean->myDesc->inaxes);

  /* Gaussian the size of raw beam in pixels */
  FWHM = in->myOTF->myDesc->beamSize  / (fabs(in->clean->myDesc->cdelt[0]));
  bmFact = 1.0;

  /* Tweak for CCB = 26.0" at low end of band */
  if (in->myOTF->myDesc->OTFType==OBIT_GBTOTF_CCB) {
    bmFact =  27.0e9 / in->myOTF->myDesc->crval[in->myOTF->myDesc->jlocf];
    FWHM = ((26.0 / 3600.0) / (fabs(in->clean->myDesc->cdelt[0]))) * bmFact;
  }

  /* Clone dirty beam with instrumental resolution 
     center of beam = peak */
  fmax = ObitFArrayMax (in->BeamPatch, beamCen);

  /* Copy */
  RawBeam = ObitFArrayCopy(in->BeamPatch, RawBeam, err);
  if (err->error) Obit_traceback_msg (err, routine, in->beam->name);

   /* Replace with Gaussian */
  ObitFArray2DCGauss (RawBeam, beamCen, FWHM);

  /* Initialize CC table on CLEAN image */
  tname = g_strconcat ("Clean table for: ",in->name, NULL);
  CCTable =  newObitTableCCValue (tname, (ObitData*)in->clean, &in->CCver, 
				  OBIT_IO_ReadWrite, 0, err);
  g_free (tname);
  
  ObitTableCCOpen (CCTable, OBIT_IO_ReadWrite, err); 
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  row = newObitTableCCRow(CCTable);

  /*  loop copy adding components */
  nrow = CCTable->myDesc->nrow;
  for (loop=in->nCCSub+1; loop<=nrow; loop++) {

    /* Read CC table */
    iRow = loop;
    ObitTableCCReadRow (CCTable, iRow, row, err);
    if (err->error) goto cleanup;

    /* Sum to model */
    pos[0] = -0.5 + (row->DeltaX / in->clean->myDesc->cdelt[0]) + in->clean->myDesc->crpix[0];
    pos[1] = -0.5 + (row->DeltaY / in->clean->myDesc->cdelt[1]) + in->clean->myDesc->crpix[1];
    ObitFArrayShiftAdd (model, pos, RawBeam, beamCen, row->Flux, model);
  } /* End component loop */

  /* Close CC table */
  ObitTableCCClose (CCTable, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

 /* Update number of components subtracted */
  in->nCCSub = nrow;

  /* subtract - since the model given has the resolution of the raw beam,
   temporarily change resolution in image header to keep from rescaling model. */
  saveBeam = in->clean->myDesc->beamMaj;
  in->clean->myDesc->beamMaj = in->scrOTF->myDesc->beamSize;
  ObitOTFUtilSubImage(in->scrOTF, in->scrOTF, model, in->clean->myDesc, err);
  in->clean->myDesc->beamMaj = saveBeam;
  if (err->error) goto cleanup;

  /* Form residual image - need to remake beam for proper normalization */
  ObitOTFUtilMakeImage (in->scrOTF, in->clean, TRUE, NULL, in->weight, err);
  if (err->error) goto cleanup;
  
 /* Cleanup */
 cleanup:
  CCTable = ObitTableCCUnref(CCTable);
  row     = ObitTableCCRowUnref(row);
  model   = ObitFArrayUnref(model);
  RawBeam = ObitFArrayUnref(RawBeam);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
}
 /* end ObitDConCleanOTFRecSub */

/**
 * Restore Clean components to residual
 * Following Parameters used from info member:
 * \li noResid    OBIT_bool scalar If defined and TRUE, set residuals to Zero
 * \li Niter      OBIT_int scalar  Maximum number of iterations [def 200].
 * \li BeamSize   OBIT_float scalar CLEAN restoring beam FWHM in asec [def 3.0].
 *                if zero then use any fitted beam size ("fitBeamSize") on the image.
 *                If < 0 then don't restore.
 *
 * \param in   Object with OTFCLEAN structures.
 * \param err  Error stack
 * \return pointer to OTFCLEAN object (ref. count updated).
 */
void ObitDConCleanOTFRecRestore (ObitDConClean *inn, ObitErr *err)
{
  ObitDConCleanOTFRec *in;
  olong iRow;
  ObitIOSize IOsize = OBIT_IO_byPlane;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong  blc[IM_MAXDIM], trc[IM_MAXDIM];
  olong  i, pos[2], loop, beamCen[2], nrow;
  ofloat DBArea, CBArea, scale, fmax, cleanSize, tempSize;
  ObitTableCC *CCTable=NULL;
  ObitTableCCRow *row=NULL;
  ObitFArray *residual=NULL, *cleanBeam=NULL;
  gboolean noresid = FALSE;
  gchar *tname;
  gchar *routine = "ObitDConCleanOTFRecRestore";

  /* Cast input to this type */
  in = (ObitDConCleanOTFRec*)inn;

 /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Is restoration desired? */
  if (in->cleanSize<0.0) return;

  /* restoring size */
  tempSize = in->cleanSize;
  ObitInfoListGetTest(in->info, "BeamSize", &type, dim, &tempSize);
  if (tempSize<=0.0) {
    ObitInfoListGetTest(in->clean->info, "fitBeamSize", &type, dim, &tempSize);
  }
 if (tempSize<=0.0) tempSize =in->myOTF->myDesc->beamSize;

  /* Convert restoring beam size into cells */
  cleanSize = tempSize / (fabs(in->clean->myDesc->cdelt[0]));

 /* Generate clean beam - clone from BeamPatch */
  if (!in->BeamPatch) ReadBP(in, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* center of beam = peak */
  fmax = ObitFArrayMax (in->BeamPatch, beamCen);

  /* Copy */
  cleanBeam = ObitFArrayCopy(in->BeamPatch, cleanBeam, err);
  if (err->error) Obit_traceback_msg (err, routine, in->beam->name);

  /* Get dirty beam area */
  DBArea = ObitFArraySum (in->BeamPatch);

  /* Replace with Gaussian */
  ObitFArray2DCGauss (cleanBeam, beamCen, cleanSize);

  /* Get CLEAN beam area */
  CBArea = ObitFArraySum (cleanBeam);

  /* Read residual image */
  /* Full field, correct plane */
  dim[0] = IM_MAXDIM;
  blc[0] = blc[1] = 1;
  for (i=0; i<IM_MAXDIM-2; i++) blc[i+2] = in->plane[i];
  ObitInfoListPut (in->clean->info, "BLC", OBIT_long, dim, blc, err); 
  trc[0] = trc[1] = 0;
  for (i=0; i<IM_MAXDIM-2; i++) trc[i+2] = in->plane[i];
  ObitInfoListPut (in->clean->info, "TRC", OBIT_long, dim, trc, err); 
  dim[0] = 1;
  ObitInfoListPut (in->clean->info, "IOBy", OBIT_long, dim, &IOsize, err);
  in->clean->extBuffer = FALSE;
  ObitImageOpen (in->clean, OBIT_IO_ReadWrite, err); 
  ObitImageRead (in->clean, NULL, err);
  ObitImageClose (in->clean, err); 
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Save clean to residual */
  residual = ObitFArrayCopy (in->clean->image, residual, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  /* Want residuals */
  noresid = FALSE;
  ObitInfoListGetTest(in->info, "noResid", &type, dim, &noresid);

  /* Scale residuals by ratio of dirty to clean beam areas to get into
     units of restored components. */
  if (CBArea>0.0) scale = DBArea / CBArea;
  else scale = 1.0;  /* something went wrong */
  if (noresid) scale = 0.0;  /* want residuals? */
  Obit_log_error(err, OBIT_InfoErr, "Scaling residuals by %f", scale);
  ObitFArraySMul (residual, scale);

  /* Initialize CC table on CLEAN image */
  tname = g_strconcat ("Clean table for: ",in->name, NULL);
  CCTable =  newObitTableCCValue (tname, (ObitData*)in->clean, &in->CCver, 
				  OBIT_IO_ReadWrite, 0, err);
  g_free (tname);
  
  ObitTableCCOpen (CCTable, OBIT_IO_ReadWrite, err); 
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  row = newObitTableCCRow(CCTable);

  /* Restoration loop */
  nrow = CCTable->myDesc->nrow;
  Obit_log_error(err, OBIT_InfoErr, "Restoring %d components", nrow); /* debug */
  for (loop=1; loop<=nrow; loop++) {

    /* Read CC table */
    iRow = loop;
    ObitTableCCReadRow (CCTable, iRow, row, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);

    /* Restore to residual */
    pos[0] = -0.5 + (row->DeltaX / in->clean->myDesc->cdelt[0]) + in->clean->myDesc->crpix[0];
    pos[1] = -0.5 + (row->DeltaY / in->clean->myDesc->cdelt[1]) + in->clean->myDesc->crpix[1];
    ObitFArrayShiftAdd (residual, pos, cleanBeam, beamCen, row->Flux, residual);
  } /* End Restoration loop */


  /* Close CC table */
  ObitTableCCClose (CCTable, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Write residual image */
  in->clean->extBuffer = TRUE;
  ObitImageOpen (in->clean, OBIT_IO_ReadWrite, err); 
  /* Update resolution on CLEAN image */
  in->clean->myDesc->beamMaj = in->cleanSize;
  in->clean->myDesc->beamMin = in->cleanSize;
  in->clean->myDesc->beamPA  = 0;
  in->clean->myDesc->niter   = nrow;
  ObitImageWrite (in->clean, residual->array, err) ;
  ObitImageClose (in->clean, err);
  in->clean->extBuffer = FALSE;
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Cleanup */
  CCTable   = ObitTableCCUnref(CCTable);
  row       = ObitTableCCRowUnref(row);
  cleanBeam = ObitFArrayUnref(cleanBeam);
  residual  = ObitFArrayUnref(residual);

} /* end ObitDConCleanOTFRecRestore */

/**
 * Select/subtract components from PxList
 * \param in   The object to deconvolve
 * \param err Obit error stack object.
 * \return TRUE if deconvolution is complete
 */
gboolean ObitDConCleanOTFRecSelect(ObitDConClean *inn, ObitErr *err)
{
  ObitDConCleanOTFRec *in;
  gboolean done = FALSE;
  olong numPix;
  olong fields[] = {1,0};
  gchar *routine = "ObitDConCleanOTFRecSelect";

  /* Cast input to this type */
  in = (ObitDConCleanOTFRec*)inn;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return done;
  g_assert (ObitDConCleanOTFRecIsA(in));

  /* How many pixels selected? */
  numPix = ObitDConCleanWindowCount (in->window, 1, err);
  numPix += 10; /* for good measure */
  if (err->error) Obit_traceback_val (err, routine, in->name, done);

  /* (Re)Create Pixel List if needed */
  if ((in->Pixels) && (numPix>in->Pixels->maxPixel))  
    ObitDConCleanPxListResize (in->Pixels, numPix, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, done);
  
  /* Load PxList */
  ObitDConCleanPxListUpdate (in->Pixels, fields, 0, 0.0, in->autoWinFlux,
			     in->window, in->BeamPatch, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, done);

  /* Set min. flux this major cycle */
  in->Pixels->minFlux[0] = in->minFlux[0];

  /* Clean */
  ObitDConCleanPxListCLEAN (in->Pixels, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, done);

  /* Completion test, niter and totalMinFlux */
  done = (in->Pixels->currentIter>=in->niter) || 
    (in->Pixels->maxResid<=in->totalMinFlux);

  return done;
} /* end ObitDConCleanOTFRecSelect */

/**
 * Get image and beam statistics 
 * \param in   The object to deconvolve
 * \param err Obit error stack object.
 */
void ObitDConCleanOTFRecPixelStats(ObitDConClean *in, ObitErr *err)
{
  const ObitDConCleanClassInfo *inClass;
  gchar *routine = "ObitDConCleanOTFRecPixelStats";

 inClass = (ObitDConCleanClassInfo*)in->ClassInfo; /* class structure */

  /* Adjust window if autoWindow */
  if (in->autoWindow) 
    inClass->ObitDConCleanAutoWindow (in, in->currentField, err);
  else 
    in->autoWinFlux = -1.0e20; 
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  return;
} /* end ObitDConCleanPixelStat */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitDConCleanOTFRecClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitDConCleanOTFRecClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitDConCleanOTFRecClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitDConCleanOTFRecClassInfoDefFn (gpointer inClass)
{
  ObitDConCleanOTFRecClassInfo *theClass = (ObitDConCleanOTFRecClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitDConCleanOTFRecClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitDConCleanOTFRecClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitDConCleanOTFRecGetClass;
  theClass->ObitClear     = (ObitClearFP)ObitDConCleanOTFRecClear;
  theClass->ObitInit      = (ObitInitFP)ObitDConCleanOTFRecInit;
  theClass->newObit       = (newObitFP)newObitDConCleanOTFRec;
  theClass->ObitCopy      = (ObitCopyFP)ObitDConCleanOTFRecCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitDConGetParms        = (ObitDConGetParmsFP)ObitDConCleanOTFRecGetParms;
  theClass->ObitDConDeconvolve      = (ObitDConDeconvolveFP)ObitDConCleanOTFRecDeconvolve;
  theClass->ObitDConCleanPixelStats = (ObitDConCleanPixelStatsFP)ObitDConCleanOTFRecPixelStats;
  theClass->ObitDConCleanSelect     = (ObitDConCleanSelectFP)ObitDConCleanOTFRecSelect;
  theClass->ObitDConCleanSub        = (ObitDConCleanSubFP)ObitDConCleanOTFRecSub;
  theClass->ObitDConCleanRestore    = (ObitDConCleanRestoreFP)ObitDConCleanOTFRecRestore;
  theClass->ObitDConCleanFlatten    = NULL;  /* Not relevant */

} /* end ObitDConCleanOTFRecClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitDConCleanOTFRecInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitDConCleanOTFRec *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->myOTF  = NULL;
  in->scrOTF = NULL;
  in->display= NULL;
  in->doRestore = TRUE;
  in->fracPeak  = 0.5;
  in->currentField = 1;;

} /* end ObitDConCleanOTFRecInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitDConCleanOTFRec* cast to an Obit*.
 */
void ObitDConCleanOTFRecClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitDConCleanOTFRec *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  in->myOTF  = ObitOTFUnref(in->myOTF);
  in->scrOTF = ObitOTFUnref(in->scrOTF);
  in->display= ObitDisplayUnref(in->display);

  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitDConCleanOTFRecClear */

/**
 * Read Beam patch into BeamPatch
 * The beam patch is symmetric about the center position allowing the 
 * beam itself not to be symmetric.
 * If the beam member is NULL, the BeamPatch is created and filled with the 
 * Gaussian size given in the clean image header.
 * \param in   The object to deconvolve
 * \param err  Obit error stack object.
 */
static void ReadBP (ObitDConCleanOTFRec* in, ObitErr *err)
{
  ObitIOCode retCode;
  ObitIOSize IOsize = OBIT_IO_byPlane;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong  blc[IM_MAXDIM], trc[IM_MAXDIM];
  olong  ablc[2], atrc[2], pos[2], naxis[2], cen[2];
  olong icenx, iceny, nx, ny, mxPatch;
  ofloat fmax, beamSize ;
  ObitImage *Beam;
  gchar *routine = "ObitDConCleanOTFRec:ReadBP";
  
  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitDConCleanIsA(in));

  /* If no Beam defined create a BeamPatch using resolution in clean image */
  if (in->beam==NULL) {
    naxis[0] =  naxis[1] = 1 + 2*in->beamPatchSize;
    cen[0] = cen[1] = in->beamPatchSize;
    /* Create */
    in->BeamPatch = ObitFArrayUnref(in->BeamPatch);
    in->BeamPatch = ObitFArrayCreate("BeamPatch", 2, naxis);
    /* Fill with Gaussian with size of that in  image */
    beamSize = in->clean->myDesc->beamMaj / fabs(in->clean->myDesc->cdelt[0]);
    ObitFArray2DCGauss (in->BeamPatch, cen, beamSize);
    return;
 } /* End create Gaussian Beam patch if no image given */
  
  /* Beam image */
  Beam = in->beam;
  
  /* Set output to full image, plane at a time */
  dim[0] = IM_MAXDIM;
  blc[0] = blc[1] = blc[2] = blc[3] = blc[4] = blc[5] = 1;
  ObitInfoListPut (Beam->info, "BLC", OBIT_long, dim, blc, err); 
  trc[0] = trc[1] = trc[2] = trc[3] = trc[4] = trc[5] = 0;
  ObitInfoListPut (Beam->info, "TRC", OBIT_long, dim, trc, err); 
  dim[0] = 1;
  ObitInfoListPut (Beam->info, "IOBy", OBIT_long, dim, &IOsize, err);
  
  retCode = ObitImageOpen (Beam, OBIT_IO_ReadOnly, err);
  if (err->error) Obit_traceback_msg (err, routine, Beam->name);

  retCode = ObitImageRead (Beam, Beam->image->array, err);
  if (err->error) Obit_traceback_msg (err, routine, Beam->name);

  /* Replace any blanks with zeroes */
  ObitFArrayDeblank (in->beam->image, 0.0);

  /* Compute region of image to take */
  /* Center pixel - make 0-rel */
  nx = Beam->myDesc->inaxes[0];
  ny = Beam->myDesc->inaxes[1];
  /* center = peak */
  fmax = ObitFArrayMax (Beam->image, pos);
  icenx = pos[0];
  iceny = pos[1];

  /* Check that Clean/residual Map and Beam have the same cell spacing */
  if ((fabs(fabs(in->clean->myDesc->cdelt[0]) - fabs(in->beam->myDesc->cdelt[0])) >
       0.01*fabs(in->clean->myDesc->cdelt[0])) ||
      (fabs(fabs(in->clean->myDesc->cdelt[1]) - fabs(in->beam->myDesc->cdelt[1])) >
       0.01*fabs(in->clean->myDesc->cdelt[1]))) {
    Obit_log_error(err, OBIT_Error, "Clean map and beam have different cell spacings");
    Obit_log_error(err, OBIT_Error, "Map %f %f Beam %f %f", 
		   3600.0*in->clean->myDesc->cdelt[0], 3600.0*in->clean->myDesc->cdelt[1],
		   3600.0*in->beam->myDesc->cdelt[0],  3600.0*in->beam->myDesc->cdelt[1]);
    return;
  }

 /* Check that center pretty close to 1.0 */
  if ((fmax<0.99) || (fmax>1.01)) {
      Obit_log_error(err, OBIT_Error, "%s Beam peak, %f not 1.0", 
		     routine, fmax);
      return;
  }

  /* How big can the patch be? */
  mxPatch = MIN (icenx-1, nx-icenx);
  mxPatch = MIN ( mxPatch, iceny-1);
  mxPatch = MIN ( mxPatch, ny-iceny);
  in->beamPatchSize =  MIN (in->beamPatchSize, mxPatch);

  /* Window as 0-rel */
  ablc[0] = icenx - in->beamPatchSize;
  atrc[0] = icenx + in->beamPatchSize;
  ablc[1] = iceny - in->beamPatchSize;
  atrc[1] = iceny + in->beamPatchSize;

  /* Save Beam patch */
  in->BeamPatch = ObitFArrayUnref(in->BeamPatch);
  in->BeamPatch = ObitFArraySubArr(Beam->image, ablc, atrc, err);
  if (err->error) Obit_traceback_msg (err, routine, Beam->name);

  retCode = ObitImageClose (Beam, err);
  if (err->error) Obit_traceback_msg (err, routine, Beam->name);

  /* Free Image array? */
  Beam->image = ObitFArrayUnref(Beam->image);
} /* end ReadBP */

