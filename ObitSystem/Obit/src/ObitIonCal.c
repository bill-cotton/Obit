/* $Id$        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2006,2008                                          */
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

#include "ObitIonCal.h"
#include "ObitTableVZ.h"
#include "ObitSkyGeom.h"
#include "ObitPBUtil.h"
#include "ObitZernike.h"
#include "ObitUVUtil.h"
#include "ObitDConCleanVis.h"
#include "ObitTableNX.h"
#include "ObitIoN2SolNTable.h"
#include "ObitTableNI.h"
#include "ObitMem.h"
#if HAVE_GSL==1  /* GSL stuff */
#include <gsl/gsl_multifit.h>
#endif /* GSL stuff */

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitIonCal.c
 * ObitIonCal class function definitions.
 * This class is derived from the Obit base class.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitIonCal";

/** Function to obtain parent ClassInfo */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/* Catalog list management */
/** Calibrator list element structure */
typedef struct { 
  /**  catalog celestial position */
  odouble ra, dec;
  /** Offset on Zernike plane (deg) */
  ofloat ZernXY[2];
  /** shift from reference position */
  ofloat shift[2];
  /** expected pixel in reference image */
  ofloat pixel[2];
  /** Estimated catalog flux density */
  ofloat flux;
  /** measured offset from expected position (deg) */
  ofloat offset[2];
  /** measured peak flux density (image units) */
  ofloat peak;
  /** measured integrated flux density (image units)  */
  ofloat fint;
  /** determined weight */
  ofloat wt;
  /** catalog quality code */
  olong qual;
  /** calibrator number (0-rel) */
  olong calNo;
  /** epoch number */
  olong epoch;
}  CalListElem; 
  
/** Private: CalListElem Constructor  */ 
static CalListElem* 
newCalListElem (odouble ra, odouble dec, ofloat ZernXy[2], ofloat shift[2],
		ofloat pixel[2], ofloat flux, ofloat offset[2],
		ofloat peak, ofloat fint,ofloat wt, olong qual,
		gint calNo, olong epoch); 

/**  Private: Update contents of an CalListElem */
static void
CalListElemUpdate (CalListElem *elem, odouble ra, odouble dec, 
		   ofloat ZernXY[2], ofloat shift[2], ofloat pixel[2], 
		   ofloat flux, ofloat offset[2], ofloat peak, ofloat fint, 
		   ofloat wt, olong qual, olong calNo, olong epoch); 

/**  Private: Print contents of an CalListElem */
static void
CalListElemPrint (CalListElem *elem, FILE *file);

/** Private: CalListElem Destructor  */ 
static void freeCalListElem(CalListElem *me); 
  
/** CalList structure. */
typedef struct {
  /** Number of entries */
  olong number;
  /** glib singly linked list */
  GSList* list;
} CalList;

/** Private: CalList structure.constructor. */
static CalList* newCalList (void);
 
/** Private: Add an CalListElem to the CalList  */
static void CalListAdd (CalList *in, CalListElem *elem);

/** Private: Remove a CalListElem from the list. */
static void CalListRemove (CalList *in, CalListElem *elem);

/** Private: Remove all items from the list. */
static void CalListClear (CalList *in);

/** Private: Remove all items from the list. */
static void CalListPrint (CalList *in, FILE *file);

/** Private: destructor. */
static void freeCalList (CalList *in);

/**
 * ClassInfo structure ObitIonCalClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitIonCalClassInfo myClassInfo = {FALSE};

/*--------------- File Global Variables  ----------------*/


/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitIonCalInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitIonCalClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitIonCalClassInfoDefFn (gpointer inClass);

/** Private: Fit positions of calibrator(s) in an image. */
static void CalPosFit (ObitImage *image, ofloat inPixel[2], olong* imsi, 
		       ofloat pixelOff[2], ofloat* souflx, ofloat* souint,
		       ObitErr *err);

/** Private: CLEAN failed - mark observations bad. */
static void CalBadTime (ObitIonCal *in, ObitImageMosaic* mosaic, 
			gint epoch, ObitErr *err);

/** Private: Lookup calibrators in an image from a catalog. */
static void 
FindCalImage (ObitImage *image, gchar *Catalog, olong catDisk, 
	      ofloat OutlierFlux, ofloat OutlierSI, olong qual, ofloat AntDiam,
	      CalList *calList, olong prtLv, ObitErr *err);

/** Private: Fit Zernike model for single time */
static ofloat
IonFit1 (gint nobs, olong* isou, ofloat* x, ofloat* y, 
	 ofloat* dx, ofloat* dy, ofloat* w, 
	 olong* ncoef, ofloat* coef, olong prtLv, ObitErr *err);

/** Private: Edit data for single time */
static gboolean
IonEdit1 (gint nobs, olong* isou, ofloat* x, ofloat* y, 
	  ofloat* dx, ofloat* dy, ofloat* w, 
	  ofloat MaxRMS, olong ncoef, ofloat* coef, olong prtLv, 
	  ObitErr *err);

/** Private: Determine timerange for next segment */
static gboolean NextTime (ObitTableNX *NXtab, ofloat solInt, 
			  ofloat timer[2], olong* suba, ObitErr* err);

/** Private: Get calibrator list from an ImageMosaic */
static void 
FindCalMosaic (ObitIonCal *in, ObitImageMosaic *mosaic, 
	      CalList *calList, olong prtLv, ObitErr *err);

/** Private: Get calibrator info from catalog */
static void 
LookupCals (ObitIonCal *in, ObitImageMosaic *mosaic, CalList *calList, 
	    olong prtLv, ObitErr *err);

/** Private: Filter and fit multiepoch offset measurements */
static ObitTableNI* 
FitAll (ObitIonCal *in, olong ncoef, 
	gint nEpoch, ofloat *timeEpoch, ofloat *timeIEpoch, olong *subaEpoch, 
	odouble refFreq, ofloat *seeing, ObitErr* err); 

/** Private: convert timerange to a printable string. */
static void TR2String (ofloat timer[2], gchar *msgBuf);

/* Private: Fit time sequence of Zernike models */
 static ObitTableNI* 
 doIonFitAll (ObitUV *inUV, ofloat MaxRMS, ofloat MinRat, gboolean doINEdit, 
	      olong ncoef, olong nsou, olong* isou, ofloat* x, ofloat* y, 
	      olong nTime, ofloat* Time, ofloat* TimeI, olong* isuba, 
	      olong n, olong* iTime, ofloat* xoff, ofloat* yoff, 
	      ofloat* flux, ofloat *wt, olong* flqual, ofloat* sint, 
	      olong prtLv, odouble refFreq, ofloat* totRMS, ObitErr* err);

/* Private: Initilize Zernike time series */
static gboolean 
initIonModel (gint nobs, olong nsou, olong nTime, olong maxcoe, olong* isou, 
	      olong* iTime, ofloat*  x, ofloat* y, ofloat* dx, ofloat* dy, 
	      ofloat* w, olong*  flqual, olong* ncoef, gboolean* gotsou, 
	      gboolean* fitsou, ofloat* soffx, ofloat* soffy, ofloat** coef, 
	      olong prtLv, ObitErr* err);

/* Private: Least Squares Fit Zernike time series and source offsets */
static void 
FitIonSeries (gint nobs, olong nsou, olong nTime, olong maxcoe, olong* isou, olong* iTime, 
	      ofloat*  x, ofloat* y, ofloat* dx, ofloat* dy, ofloat* w, olong*  ncoef, 
	      gboolean* gotsou, gboolean* fitsou, ofloat* soffx, ofloat* soffy, 
	      ofloat** coef, ofloat* trms, olong prtLv, ObitErr* err) ;

/* Private: Edit data based on Zernike time series  */
static gboolean 
IonEditSeries (gint nobs, olong nsou, olong nTime, 
	       olong maxcoe, olong* isou, olong* iTime, ofloat*  x, ofloat* y, 
	       ofloat* dx, ofloat* dy, ofloat* w, ofloat MaxRMS, olong* ncoef, 
	       ofloat* soffx, ofloat* soffy, ofloat** coef, olong prtLv, ObitErr* err);

/* Private: Time series editing based on amplitudes */
static gboolean 
IonEditAmp (gint nobs, olong nsou, olong nTime, olong* isou, olong* iTime, 
	    ofloat MinRat, ofloat*  flux, olong* ncoef, ofloat* wt, 
	    olong prtLv, ObitErr* err);

/* Private: Edit data for a given source based on Zernike time series */
static void 
IonEditSource (gint source, olong nobs, olong nTime, olong maxcoe, olong nsou, 
	       olong* isou, olong* iTime, ofloat*  x, ofloat* y, 
	       ofloat* dx, ofloat* dy, ofloat* w, ofloat MaxRMS, olong* ncoef, 
	       ofloat* soffx, ofloat* soffy, ofloat** coef, 
	       ofloat *gdx, ofloat *gdy, olong* stoss, olong prtLv, ObitErr* err);


static void PosImage (ObitImage *image, ofloat pixel[2], ofloat minflux, 
		      ofloat mxcdis, olong FitSize, 
		      ofloat offset[2], ofloat *peak, ofloat *fint, ObitErr *err);
static olong pfit (ofloat a[9][9], ofloat *s, ofloat dx[2], ofloat fblank);
static olong momnt (ofloat ara[9][9], ofloat x, ofloat y, olong nx, olong ny, 
		   ofloat momar[6], ofloat fblank);
static void matvmul (ofloat m[3][3], ofloat vi[3], ofloat vo[3], olong n);
static void rd2zer (odouble ra, odouble dec, ofloat xshift, ofloat yshift, 
		    ofloat* xzer, ofloat* yzer, olong *ierr);
static void zerrot (odouble ra, odouble dec, ofloat xshift, ofloat yshift,
		    ofloat* rotate);
static void  FitZernike (olong nobs, olong* isou, ofloat *x, ofloat *y, 
			 ofloat *dx, ofloat *dy, ofloat *wt, 
			 olong nZern, ofloat *coef);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitIonCal* newObitIonCal (gchar* name)
{
  ObitIonCal* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitIonCalClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitIonCal));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitIonCalInit((gpointer)out);

 return out;
} /* end newObitIonCal */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitIonCalGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitIonCalClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitIonCalGetClass */

/**
 * Make a deep copy of an ObitIonCal.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitIonCal* ObitIonCalCopy  (ObitIonCal *in, ObitIonCal *out, ObitErr *err)
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
    out = newObitIonCal(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  out->myData = ObitUVRef(in->myData);

  return out;
} /* end ObitIonCalCopy */

/**
 * Make a copy of a object but do not copy the actual data
 * This is useful to create an IonCal similar to the input one.
 * \param in  The object to copy
 * \param out An existing object pointer for output, must be defined.
 * \param err Obit error stack object.
 */
void ObitIonCalClone  (ObitIonCal *in, ObitIonCal *out, ObitErr *err)
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

  /*  copy this class */
  out->myData = ObitUVRef(in->myData);

} /* end ObitIonCalClone */

/**
 * Creates an ObitIonCal 
 * \param name  An optional name for the object.
 * \return the new object.
 */
ObitIonCal* ObitIonCalCreate (gchar* name)
{
  ObitIonCal* out;

  /* Create basic structure */
  out = newObitIonCal (name);

  return out;
} /* end ObitIonCalCreate */

/**
 * Attach UV data, unreferences any old data
 * \param in    Object to attach UV data to
 * \param inUV  UV data to attach
 */
void ObitIonCalSetData (ObitIonCal *in, ObitUV *inUV)
{
  in->myData = ObitUVUnref(in->myData);  /* out with the old */
  in->myData = ObitUVRef(inUV);          /* in with the new */
} /* end ObitIonCalSetData */

/**
 * Fills the calList member with a list of calibrator sources 
 * selected from a given catalog whose position are inside the
 * field of view of image.
 * Previous contents of the CalList are cleared
 * \param in   IonCal object 
 *   Control parameters are on the info member.
 * \li Catalog     OBIT_char  (?,1,1)  AIPSVZ format FITS catalog for defining outliers, 
 *                 'Default' or blank = use default catalog.
 * \li catDisk     OBIT_int   (1,1,1) FITS disk for catalog [def 1]
 * \li OutlierFlux OBIT_float (1,1,1)  Minimum estimated flux density include cal. fields
 *                  from Catalog. [default 0.1 Jy ]
 * \li OutlierSI   OBIT_float (1,1,1) Spectral index to use to convert catalog flux 
 *                   density to observed frequency.  [default = -0.75]
 * \li MaxQual     OBIT_int   (1,1,1) Max. cal. quality code [def 1]
 * \li prtLv       OBIT_int   (1,1,1) Print level >=2 => give list [def 0]
 *
 * calList Calibrator list each element of which has:
 * \li ra      position RA (deg) of calibrators
 * \li dec     position Dec (deg) of calibrators
 * \li shift   Offset in field [x,y] on unit Zernike circle of position
 * \li pixel   Expected pixel in reference image
 * \li flux    Estimated catalog flux density 
 * \li offset  Measured offset [x,y] from expected position (deg)
 * \li peak    Measured peak flux density (image units)
 * \li fint    Measured integrated flux density (image units)
 * \li wt      Determined weight
 * \li qual    Catalog quality code
 * \param image   Image object 
 * \param err     Error stack
 */
void ObitIonCalFindImage (ObitIonCal *in, ObitImage* image, ObitErr* err)  
{
  CalList *calList;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  gchar Catalog[257];
  olong catDisk, qual, prtLv;
  ofloat OutlierFlux, OutlierSI, AntDiam;
  gchar *routine = "ObitIonCalFindImage";

  /* Error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return ;  /* previous error? */
  g_assert(ObitIonCalIsA(in));

   /* Clear any entries in calList */
  calList = (CalList*)in->calList;
  CalListClear (calList);

  /* Control information. */
  sprintf (Catalog, "Default");
  ObitInfoListGetTest(in->info, "Catalog", &type, dim,  Catalog);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  if (!strncmp(Catalog, "    ", 4)) sprintf (Catalog, "Default");
  if (!strncmp(Catalog, "Default", 7)) sprintf (Catalog, "NVSSVZ.FIT");
  catDisk = 1;
  ObitInfoListGetTest(in->info, "catDisk",   &type, dim, &catDisk);
  OutlierFlux = 1.0;
  ObitInfoListGetTest(in->info, "OutlierFlux", &type, dim, &OutlierFlux);
  OutlierSI = -0.75;
  ObitInfoListGetTest(in->info, "OutlierSI", &type, dim, &OutlierSI);
  AntDiam = 25.0;
  ObitInfoListGetTest(in->info, "AntDiam", &type, dim, &AntDiam);
  qual = 1;
  ObitInfoListGetTest(in->info, "MaxQual",   &type, dim, &qual);
  prtLv = 0;
  ObitInfoListGetTest(in->info, "prtLv",   &type, dim, &prtLv);

  /* Search Catalog */
  FindCalImage (image, Catalog, catDisk, OutlierFlux, OutlierSI, qual, 
		AntDiam, calList, prtLv, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
} /* end ObitIonCalFindImage */

/**
 * Determines calibrator position offsets from a single image.
 * Given an image containing calibrator sources, return fitted fluxes 
 * and position offsets.  
 * Resultant positions are all referred to a tangent plane at the  
 * pointing center, this is referred to as the Zernike plane as this  
 * plane will be used to fit the phase screen.  The "Zernike Unit  
 * Circle" defines the phase screen.  Source position offsets are in  
 * the X and Y (RA, Dec) as defined by this plane.  
 * Routine translated from the AIPSish CALPOS.FOR/CALPOS  
 * \param in   IonCal object 
 *   Control parameters are on the info member.
 * \li "FitDist"  OBIT_int   (1,1,1) dist, from expected location to search 
 *                                   asec [10 pixels]
 * \li "MinPeak"  OBIT_float (1,1,1) Min. acceptable image peak (Jy) [1.0]
 * \li "MaxDist"  OBIT_float (1,1,1) Max. distance (deg/10) to accept calibrator [1.]
 * \li "MaxWt"    OBIT_float (1,1,1) Max. weight [10.0]
 * \li prtLv      OBIT_int   (1,1,1) Print level >=1 => give fits
 * \param image   Image object 
 * \param calList Calibrator list each element of which has:
 * \li ra      position RA (deg) of calibrators
 * \li dec     position Dec (deg) of calibrators
 * \li shift   Offset in field [x,y] on unit Zernike circle of position
 * \li pixel   Expected pixel in reference image
 * \li flux    Estimated catalog flux density 
 * \li offset  Measured offset [x,y] from expected position (deg)
 * \li peak    Measured peak flux density (image units)
 * \li fint    Measured integrated flux density (image units)
 * \li wt      Determined weight
 * \li qual    Catalog quality code
 * \param err     Error stack
 */
void ObitIonCalPosMul (ObitIonCal *in, ObitImage* image, ObitErr* err)  
{
  olong   i, FitSize, good, imsi[2], epoch=1, prtLv;
  olong blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  ofloat FitDist, flux;
  ofloat mxcdis,maxwt;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  ObitIOSize IOBy;
  gboolean bad;
  ofloat offset[2], peak, fint,wt, dist;
  ofloat area, pixelOff[2];
  CalListElem *elem=NULL, *nelem=NULL;
  GSList  *tmp;
  CalList *calList;
  gchar *routine = "ObitIonCalPosMul";

  /* Error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return ;  /* previous error? */
  g_assert(ObitIonCalIsA(in));

  calList = (CalList*)in->calList;

  /* Control information. */
  maxwt = 10.0;
  ObitInfoListGetTest(in->info, "MaxWt",   &type, dim, &maxwt);
  flux = 1.0;
  ObitInfoListGetTest(in->info, "MinPeak", &type, dim, &flux);
  mxcdis = 1.0;
  ObitInfoListGetTest(in->info, "MaxDist", &type, dim, &mxcdis);
  FitDist = 0.0;
  ObitInfoListGetTest(in->info, "FitDist", &type, dim, &FitDist);
  prtLv = 0;
  ObitInfoListGetTest(in->info, "prtLv",   &type, dim, &prtLv);

  /* Open and read image full plane */
  IOBy = OBIT_IO_byPlane;
  dim[0] = 1;
  ObitInfoListPut (image->info, "IOBy", OBIT_long, dim, &IOBy, err);
  dim[0] = 7;
  for (i=0; i<IM_MAXDIM; i++) {blc[i] = 1; trc[i] = 0;}
  ObitInfoListPut (image->info, "BLC", OBIT_long, dim, blc, err); 
  ObitInfoListPut (image->info, "TRC", OBIT_long, dim, trc, err);
  image->extBuffer = FALSE;  /* Make sure it has buffer */
  if (err->error) Obit_traceback_msg (err, routine, image->name);
  if ((ObitImageOpen (image, OBIT_IO_ReadOnly, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, 
		   "ERROR opening image %s", image->name);
    return;
  }
  ObitImageRead (image, NULL, err);
  if (err->error) Obit_traceback_msg (err, routine, image->name);

  /* Convert FitDist to pixels */
  if (FitDist<=0.0) FitDist = 10.0 * fabs(image->myDesc->cdelt[0]);
  FitSize = (olong)((FitDist / fabs(image->myDesc->cdelt[0]))+0.5);
  FitSize = MAX (10, FitSize);
  imsi[0] = imsi[1] = FitSize;

  /* Loop over sources */
  good = 0;
  tmp = calList->list;
  while (tmp!=NULL) {   /* loop 100 */
    elem = (CalListElem*)tmp->data;
    /* if find epoch > 0 then quit loop */
    if (elem->epoch>0) break;
    bad = FALSE;

    /* fit source */
    CalPosFit (image, elem->pixel, imsi, pixelOff, &peak, &fint, err);
    if (err->error) Obit_traceback_msg (err, routine, image->name);

    /* find anything? */
    if (peak<flux) { /* Nope */
      bad = TRUE;
    } else { /* yes */

      /* Offset */
      offset[0] = pixelOff[0] * image->myDesc->cdelt[0];
      offset[1] = pixelOff[1] * image->myDesc->cdelt[1];

      /* Normalize integrated flux by beam area in pixels */
      area = 1.1331 * image->myDesc->beamMaj * image->myDesc->beamMin /  
	(fabs (image->myDesc->cdelt[0]) * fabs (image->myDesc->cdelt[1])) ;
      if (area  <  0.001) area = 1.0;
      fint /= area;

      /* Within maximum distance? */
      dist = sqrt(elem->ZernXY[0]*elem->ZernXY[0] + elem->ZernXY[1]*elem->ZernXY[1]);
      if (dist <= mxcdis)  wt = MIN (maxwt, peak);
      else wt = 0.0;
      if (bad) wt = 0.0;

      /* Update CalList */
      nelem = newCalListElem (elem->ra, elem->dec, elem->ZernXY, elem->shift, 
			      elem->pixel, elem->flux, offset, peak, fint, wt, 
			      elem->qual, elem->calNo, epoch);
      CalListAdd (calList, nelem);
      good++;  /* This one OK */

      /* diagnostics? */
      if (prtLv>=1) {
	Obit_log_error(err, OBIT_InfoErr, 
		       "Fit E=%5d C=%5d Z= %7.4f %7.4f off=%7.2f %7.2f Peak=%6.2f int=%6.2f qual=%d",
		       epoch, elem->calNo, elem->ZernXY[0], elem->ZernXY[1], 
		       offset[0]*3600.0, offset[1]*3600.0, peak, fint, elem->qual);
      }
    } /* End if detected */

    tmp = g_slist_next(tmp);  /* next item in list */
 } /* end loop  L100: over calibrator*/

  /* Close image */
  ObitImageClose (image, err);
  if (err->error) Obit_traceback_msg (err, routine, image->name);
  image->image = ObitImageUnref(image->image);   /* Free buffer */

  /* tell how many */
  Obit_log_error(err, OBIT_InfoErr, "%s: Fitted %d calibrator offsets", 
		 routine, good);

} /* end of routine ObitIonCalPosMul */ 

/**
 * Fit Zernike model to single epoch image distortion
 * Uses fitted data in calList member.
 * Iteratively edits most discrepant point if needed to get RMS 
 * residual down to MaxRMS.
 * \param in   IonCal object 
 *   Control parameters are on the info member.
 * \li "nZern"  OBIT_int   (1,1,1) Zernike polynomial order requested [def 5]
 * \li "MaxRMS" OBIT_float (1,1,1) Target RMS residual (asec), default 10 asec
 * \li prtLv    OBIT_int   (1,1,1) Print level >=3 => give fitting diagnostics
 *                                 [def. 0]
 * \param epoch   1-rel time index i measurements
 * \param coef    [out] Fitted coefficients
 * \param err     Error stack
 * \return RMS residual in deg, -1 on error
 */
ofloat ObitIonCalFit1 (ObitIonCal *in, olong epoch, ofloat *coef, 
		       ObitErr* err)  
{
  CalListElem *elem=NULL;
  ofloat out = -1.0;
  GSList  *tmp;
  olong nZern, count, badCount, prtLv;
  gboolean doMore;
  ofloat rms, MaxRMS, MaxRMS2;
  ofloat *x=NULL, *y=NULL, *dx=NULL, *dy=NULL, *w=NULL;
  olong *isou=NULL;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  CalList *calList;
  gchar *routine = "ObitIonCalFit1";

  /* Error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return out;  /* previous error? */
  g_assert(ObitIonCalIsA(in));

  calList = (CalList*)in->calList;

  /* Control information. */
  nZern = 5;
  ObitInfoListGetTest(in->info, "nZern",   &type, dim, &nZern);
  /* Coerce to range [2,17] */
  nZern = MAX (2, MIN (17, nZern));
  MaxRMS = 10.0;
  ObitInfoListGetTest(in->info, "MaxRMS",   &type, dim, &MaxRMS);
  MaxRMS /= 3600.0;  /* to deg */
  MaxRMS2 = 2.0 * MaxRMS; /* relax for finding sources */
  prtLv = 0;
  ObitInfoListGetTest(in->info, "prtLv",   &type, dim, &prtLv);

  /* How many measurements at epoch epoch? */
  count = 0;
  /* Loop over calList */
  tmp = calList->list;
  while (tmp!=NULL) {
    elem = (CalListElem*)tmp->data;
    if (elem->epoch==epoch) count++;
    if (elem->epoch>epoch) break; /* done? */
    tmp = g_slist_next(tmp);  /* next item in list */
  } 

  /* allocate memory */
  x  = g_malloc(count*sizeof(ofloat));
  y  = g_malloc(count*sizeof(ofloat));
  dx = g_malloc(count*sizeof(ofloat));
  dy = g_malloc(count*sizeof(ofloat));
  w  = g_malloc(count*sizeof(ofloat));
  isou = g_malloc(count*sizeof(olong));

  /* Copy data to arrays - Loop over calList */
  tmp = calList->list;
  count = 0;
  while (tmp!=NULL) {
    elem = (CalListElem*)tmp->data;
    if (elem->epoch==epoch) {
      /* Unit circle = 10 deg */
      x[count]    = elem->ZernXY[0];
      y[count]    = elem->ZernXY[1];
      dx[count]   = elem->offset[0];
      dy[count]   = elem->offset[1];
      w[count]    = elem->wt;
      isou[count] = count;
      count++;
    }
    if (elem->epoch>epoch) break; /* done? */
    tmp = g_slist_next(tmp);  /* next item in list */
  } 

  /* Fit model */
  doMore = TRUE;
  badCount = 0;
  while (doMore) {
    rms = IonFit1 (count, isou, x, y, dx, dy, w, &nZern, coef, prtLv, err);
    if (rms<MaxRMS2) break;
    doMore = IonEdit1 (count, isou, x, y, dx, dy, w, MaxRMS2, nZern, coef, prtLv, err);
    if (doMore) badCount++;
  }

  /* tell how many flagged */
  Obit_log_error(err, OBIT_InfoErr, "%s: Flagged %d of %d calibrators", 
		 routine, badCount, count);

  Obit_log_error(err, OBIT_InfoErr, "%s: Final RMS %f", 
		 routine, 3600.0*rms);
  /* cleanup */
  if (x)    g_free(x);
  if (y)    g_free(y);
  if (dy)   g_free(dx);
  if (dy)   g_free(dy);
  if (w)    g_free(w);
  if (isou) g_free(isou);

  return rms;
} /* end of routine ObitIonCalFit1 */ 

/**
 * Ionospheric calibration  for a uv data set.
 * Loops over time slices, imaging and deconvolving selected fields.  
 * Then determines position offsets and fits an ionospheric model.  
 * Results are stored in an 'NI' table attached to inUV.
 * Current maximum 1024 epochs.
 * Routine translated from the AIPSish IONCAL.FOR/IONCAL  
 * \param in        IonCal object, must have myData UV data attached.
 *                  If myData is a multisource file, then
 *                  selection (1 source) , calibration and editing controls
 *                  must be set on the info member.
 *   Control parameters on the myData info member.
 * \li Catalog     OBIT_char  (?,1,1)  AIPSVZ format FITS catalog for defining outliers, 
 *                 'Default' or blank = use default catalog.
 * \li catDisk     OBIT_int   (1,1,1)  FITS disk for catalog [def 1]
 * \li OutlierDist OBIT_float (1,1,1)  How far from pointing to add calibrators
 * \li OutlierFlux OBIT_float (1,1,1)  Minimum estimated flux density include outlier fields
 *                                     from Catalog. [default 0.1 Jy ]
 * \li OutlierSI   OBIT_float (1,1,1)  Spectral index to use to convert catalog flux 
 *                                     density to observed frequency.  [default = -0.75]
 * \li Niter       OBIT_int   (1,1,1)  Max. number of components to clean
 * \li minFlux     OBIT_float (1,1,1)  Minimum flux density to CLEAN
 * \li autoWindow  OBIT_boolean (1,1,1)True if autoWindow feature wanted.
 * \li dispURL"    OBIT_string (1,1,1) URL of display server
 *
 *   Control parameters on the in->info member.
 * \li nZern       OBIT_int   (1,1,1)  Zernike polynomial order requested [def 5]
 * \li MaxQual     OBIT_int   (1,1,1)  Max. cal. quality code [def 1]
 * \li prtLv       OBIT_int   (1,1,1)  Print level >=2 => give list [def 0]
 * \li MaxRMS      OBIT_float (1,1,1)  Maximum allowable RMS in arcsec. [def 20]
 * \li MinRat      OBIT_float (1,1,1)  Minimum acceptable ratio to average flux  [def 0.1]
 * \li FitDist     OBIT_int   (1,1,1)  Dist, from expected location to search 
 *                                      asec [10 pixels]
 * \li MinPeak     OBIT_float (1,1,1)  Min. acceptable image peak (Jy) [1.0]
 *                 If not given  OutlierFlux is used
 * \li MaxDist     OBIT_float (1,1,1)  Max. distance (deg/10) to accept calibrator [1.]
 * \li MaxWt       OBIT_float (1,1,1)  Max. weight [10.0]
 * \li doINEdit    OBIT_boolean (1,1,1) If true flag solutions for which the seeing 
 *                                     residual could not be determined or exceeds 
 *                                     MaxRMS [def TRUE] 
 * \li doSN        OBIT_boolean (1,1,1) If true, convert IN table to an SN table 
 *                                     attached to inUV. [def False]
 * \li solInt      OBIT_float  (1,1,1) Solution interval (min). [def 1]
 * \param err    Error code: 0 => ok, -1 => all data flagged 
 */
void ObitIonCaldoCal (ObitIonCal*in, ObitErr* err)
{
  const ObitDConCleanVisClassInfo *ClnClass;
  gboolean doSN;
  ofloat   solInt, maxWt, Flux, mxcDis, OutlierFlux;
  olong     qual, prtLv, cprtLv;
  olong     ncal,  ncoef, suba;
  olong    i, ver;
#define MAXEPOCH 1024   /* Maximum number of epochs */
  olong     nZern, nEpoch=0, subaEpoch[MAXEPOCH];
  ofloat   timeEpoch[MAXEPOCH], timeIEpoch[MAXEPOCH]; 
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  ofloat MaxRMS, off[2], timer[2] = {-1.0e20, 0.0};
  odouble refFreq=0.0;
  ObitUV *inUV = NULL;
  ObitDConCleanVis* myClean=NULL;
  ObitImageMosaic* myMosaic=NULL;
  ObitTableNX *NXtab=NULL;
  ObitTableSN *SNtab=NULL;
  ObitTableNI *NItab=NULL;
  CalList *calList;
  ofloat FOV, oldFOV, Beam[3], seeing, maxGap;
  gboolean badTime, *gotit=NULL, Tr=TRUE, Fl=FALSE;;
  odouble obsra, obsdec;
  gchar msgtxt[101], *Stokes="I   ";
  /* Parameters to copy from inUV to CLEAN object */
  gchar *ParmList[] = {"Niter", "minFlux", "dispURL", "autoWindow", NULL};
  /* Parameters to copy from inUV to IonCal object */
  gchar *IonParmList[] = {"Catalog", "catDisk", "OutlierFlux", "OutlierSI", 
			  NULL};
  gchar *routine = "ObitIonCaldoCal";

  /*----------------------------------------------------------------------- */
  /* Previous error condition? */
  if (err->error) return ;
  /* Error checks */
  inUV = in->myData;   /* data on in */
  g_assert(ObitUVIsA(inUV));
  
  /* Get control parameters */
  qual = 1;
  ObitInfoListGetTest(in->info, "MaxQual", &type, dim, &qual);
  prtLv = 0;
  ObitInfoListGetTest(in->info, "prtLv", &type, dim, &prtLv);
  doSN = FALSE;
  ObitInfoListGetTest(in->info, "doSN", &type, dim, &doSN);
  solInt = 1.0;
  ObitInfoListGetTest(in->info, "solInt", &type, dim, &solInt);
  if (solInt<=0.0) solInt = 1.0;
  maxWt = 10.0;
  ObitInfoListGetTest(in->info, "MaxWt",   &type, dim, &maxWt);
  /* Default MinPeak is &OutlierFlux */
  Flux = 1.0;
  if (!ObitInfoListGetTest(in->info, "MinPeak", &type, dim, &Flux)) {
    OutlierFlux = 1.0;
    ObitInfoListGetTest(in->info, "OutlierFlux", &type, dim, &OutlierFlux);
    Flux = OutlierFlux;
    dim[0] = dim[1] = dim[2] = 1;
    ObitInfoListAlwaysPut(in->info, "MinPeak", OBIT_float, dim, &Flux); 
  }
  mxcDis = 1.0;
  ObitInfoListGetTest(in->info, "MaxDist", &type, dim, &mxcDis);
  MaxRMS = 20.0;
  ObitInfoListGetTest(in->info, "MaxRMS",  &type, dim, &MaxRMS);
  nZern = 5;
  ObitInfoListGetTest(in->info, "nZern",   &type, dim, &nZern);

  /* Make sure inUV indexed */
  maxGap = solInt;  /* Max. time gap in indexing */
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut(inUV->info, "maxGap", OBIT_float, dim, &maxGap);
  ObitUVUtilIndex(inUV, err);
  if (err->error) Obit_traceback_msg (err, routine, inUV->name);

  /* Get NX table */
  ver = 1;
  NXtab = newObitTableNXValue ("IonCal NX table", (ObitData*)inUV, &ver, 
			       OBIT_IO_ReadOnly, err);
  if (err->error) goto cleanup;
      
  /* Set FOV to 0 then restore */
  oldFOV = 0.0;
  ObitInfoListGetTest(inUV->info, "FOV", &type, dim, &oldFOV);
  FOV = 0.0;
  ObitInfoListAlwaysPut (inUV->info, "FOV", type,   dim, &FOV);

  /* Enable data selection by timerange */
  dim[0] = 1;
  ObitInfoListAlwaysPut(inUV->info, "doCalSelect", OBIT_bool, dim, &Tr);
  /* Stokes 'I' */
  dim[0] = 4;
  ObitInfoListAlwaysPut(inUV->info, "Stokes", OBIT_string, dim, Stokes);
  
  /* Create Clean (and calibrator mosaic) */
  myClean = ObitDConCleanVisCreate ("IonCal Clean", inUV, err);
  if (err->error) goto cleanup;
  ClnClass = (ObitDConCleanVisClassInfo*)myClean->ClassInfo; /* class structure */


  /* Reset FOV to previous value */
  ObitInfoListAlwaysPut (inUV->info, "FOV", type,   dim, &oldFOV);

  /* Local pointer to calibrator mosaic */
  myMosaic = ObitUVImagerGetMosaic (myClean->imager, err);
  if (err->error) goto cleanup;

  /* Only summary CLEAN messages */
  cprtLv = 1; dim[0] = dim[1] = 1;
  if (prtLv>=4) cprtLv = prtLv;
  ObitInfoListAlwaysPut (myClean->info, "prtLv", OBIT_long, dim, &cprtLv);

  /* No need to cross restore */
  dim[0] = 1;
  ObitInfoListAlwaysPut(myClean->info, "doXRestore", OBIT_bool, dim, &Fl);

  /* ... or flatten */
  dim[0] = 1;
  ObitInfoListAlwaysPut(myClean->info, "doFlatten", OBIT_bool, dim, &Fl);

  /* Copy CLEAN control parameters */
  ObitInfoListCopyList (inUV->info, myClean->info, ParmList);

  /* Copy IonCal  control parameters */
  ObitInfoListCopyList (inUV->info, in->info, IonParmList);

  /* Default Beam */
  Beam[0] = Beam[1] = Beam[2] = 0.0;
  dim[0] = 3; dim[1] = 1;
  ObitInfoListAlwaysPut (myClean->info, "Beam", OBIT_float, dim, Beam);

  /* Diagnostic info*/
  if (prtLv>0) {
    g_snprintf (msgtxt,100,"calibrating ionosphere, max RMS seeing= %6.2f qual=%5d", 
		MaxRMS, qual);
    Obit_log_error(err, OBIT_InfoWarn, msgtxt);
    
    g_snprintf (msgtxt,100,"observing date = %s", inUV->myDesc->obsdat);
    Obit_log_error(err, OBIT_InfoWarn, msgtxt);
    ObitImageDescGetPoint (myMosaic->images[0]->myDesc, &obsra, &obsdec);
    g_snprintf (msgtxt,100,"pointing ra=%10.6f dec=%10.6f deg", obsra, obsdec);
    Obit_log_error(err, OBIT_InfoWarn, msgtxt);
  }

  ObitErrLog(err); /* show any messages on err */
  
  /* Clear any entries in calList */
  calList = (CalList*)in->calList;
  CalListClear (calList);

  /* Generate calibrator list from mosaic */
  FindCalMosaic (in, myMosaic, calList, prtLv, err);
  if (err->error) goto cleanup;

  /* arrays per calibrator */
  ncal = myMosaic->numberImages;
  gotit = g_malloc(ncal*sizeof(gboolean));  /* calibrator ever found? */
  for (i=0; i<ncal; i++) gotit[i] = FALSE; /* Not found anything yet. */

  /* How many coefficients to solve  for?  */
  ncoef = nZern;

  ObitErrLog(err); /* show any messages on err */
  /* Loop over time slices */
  for (i=1; i<=100000; i++) { /* loop 500 */
    if (!NextTime(NXtab, solInt, timer, &suba, err)) break;
    if (err->error) goto cleanup;

    /* Set timerange, subarray */
    dim[0] = 2; dim[1] = 1;
    ObitInfoListAlwaysPut (inUV->info, "timeRange", OBIT_float, dim, timer);
    dim[0] = 1;
    ObitInfoListAlwaysPut (inUV->info, "Subarray", OBIT_long, dim, &suba);

    /* Use fitted beam */
    dim[0] = 3; dim[1] = 1;
    ObitInfoListAlwaysPut (myClean->info, "Beam", OBIT_float, dim, Beam);

    /* Reset Windows first time */
    if (i==1) ClnClass->ObitDConCleanDefWindow((ObitDConClean*)myClean, err);
    if (err->error) goto cleanup;

    ObitErrLog(err);  /* Make sure message stack shown */

    /* Image/deconvolve this one */
    ClnClass->ObitDConDeconvolve ((ObitDCon*)myClean, err);
    /* Be somewhat tolerant of failures here */
    if (err->error) {
       badTime = TRUE;
       ObitErrClearErr (err); /* Clear error messages and condition */
    } else badTime = FALSE;

    /* Diagnostic info*/
    if (prtLv>0) {
      /* Tell time range */
      TR2String (timer, msgtxt);
      Obit_log_error(err, OBIT_InfoErr, "Timerange %s", msgtxt);
      Obit_log_error(err, OBIT_InfoErr, 
		     "Total CLEAN flux density = %8.1f, resid = %8.1f Jy", 
		     myClean->Pixels->totalFlux, myClean->Pixels->maxResid);
      ObitErrLog(err);
    } 

    /* Save epoch information */
    nEpoch = i;
    if (nEpoch<=MAXEPOCH) {
      subaEpoch[nEpoch-1]  = suba;
      timeEpoch[nEpoch-1]  = (timer[0]+timer[1]) * 0.5;
      timeIEpoch[nEpoch-1] = timer[1] - timer[0];
    } else { /* blew arrays */
      Obit_log_error(err, OBIT_InfoWarn, 
		     "%s: Overflowed internal arrays (%d) stopping IonCal", 
		     routine, MAXEPOCH);
      break;  /* Full, stop calibrating */
    }

    /* save image reference frequency */
    refFreq = myMosaic->images[0]->myDesc->crval[myMosaic->images[0]->myDesc->jlocf];

    /* Do fits to each image in mosaic if CLEAN succeeded */
    if (!badTime) {
      /* Digest CLEAN images */
      ObitIonCalPosMosaic (in, myMosaic, nEpoch, err);
    } else {
      /* Mark entries bad */
      CalBadTime (in, myMosaic, nEpoch, err);
      Obit_log_error(err, OBIT_InfoWarn, "%s: CLEAN Failed this time slice", 
		     routine);
      /* Mark entries bad */
    }
    if (err->error) goto cleanup;

    /* Keep track of epoch times, subarrays */

    /* Diagnostic info*/
    if (prtLv>0) {
      /* Print solutions for this epoch */
      ObitErrLog(err); /* show any messages on err */
    } 
  } /* end loop  L500: over snapshots */

  /* Fit model, write IN table */
  NItab = FitAll (in, ncoef, nEpoch, timeEpoch, timeIEpoch, subaEpoch, 
		  refFreq, &seeing, err);
  
  /* Save seeing */
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut(inUV->info, "seeing", OBIT_float, dim, &seeing);


  /* If requested convert to SN */
  if (doSN) {
    off[0] = off[1] = 0.0;
    SNtab = ObitIoN2SolNTableConvert (inUV, NItab, NULL, off, err);
    if (err->error) goto cleanup;
  } 
 
  /* Cleanup */
 cleanup:  myClean = ObitDConCleanVisUnref(myClean);
  NXtab = ObitTableNXUnref(NXtab);
  SNtab = ObitTableSNUnref(SNtab);
  NItab = ObitTableNIUnref(NItab);
  if (gotit) g_free(gotit);
  ObitImageMosaicZapImage (myMosaic, -1, err);
  myMosaic = ObitImageMosaicUnref (myMosaic);
  if (err->error) Obit_traceback_msg (err, routine, inUV->name);
  
  /* Reset timerange, subarray */
  dim[0] = 2; dim[1] = 1;
  timer[0] = -1.0e20; timer[1] = 1.0e20;
  ObitInfoListAlwaysPut (inUV->info, "timeRange", OBIT_float, dim, timer);
  dim[0] = 1;
  suba = 0;
  ObitInfoListAlwaysPut (inUV->info, "Subarray", OBIT_long, dim, &suba);

} /* end of routine ObitIonCaldoCal */ 

/**
 * Determines calibrator position offsets from images in a mosaic.
 * Asumes given a calList with entries corresponding to entries in an
 * Image mosaic, fit the actual positioins in the mosaic and write new 
 * entries in the CalList with offsets for acceptable fits.
 * Resultant positions are all referred to a tangent plane at the  
 * pointing center, this is referred to as the Zernike plane as this  
 * plane will be used to fit the phase screen.  The "Zernike Unit  
 * Circle" defines the phase screen.  Source position offsets are in  
 * the X and Y (RA, Dec) as defined by this plane.  
 * Routine adapted from the AIPSish CALPOS.FOR/CALPOS  
 * \param in   IonCal object 
 *   Control parameters are on the info member.
 * \li "nZern"  OBIT_int   (1,1,1) Zernike polynomial order requested [def 5]
 * \li FitDist  OBIT_int   (1,1,1) dist, from expected location to search 
 *                                 asec [10 pixels]  
 * \li MinPeak  OBIT_float (1,1,1) Min. acceptable image peak (Jy) [1.0]
 * \li MaxDist  OBIT_float (1,1,1) Max. distance (deg/10) to accept calibrator [1.]
 * \li MaxWt    OBIT_float (1,1,1) Max. weight [10.0]
 * \li MaxQual  OBIT_int   (1,1,1) Max. cal. quality code [def 1]
 * \li prtLv    OBIT_int   (1,1,1) Print level >=1 => give fits
 * \param mosaic   Image mosaic object 
 * \param calList Calibrator list each element of which has:
 * \li ra      position RA (deg) of calibrators
 * \li dec     position Dec (deg) of calibrators
 * \li shift   Offset in field [x,y] on unit Zernike circle of position
 * \li pixel   Expected pixel in reference image
 * \li flux    Estimated catalog flux density 
 * \li offset  Measured offset [x,y] from expected position (deg)
 * \li peak    Measured peak flux density (image units)
 * \li fint    Measured integrated flux density (image units)
 * \li wt      Determined weight
 * \li qual    Catalog quality code
 * \param err     Error stack
 * \param epoch   Epoch number
 */
void ObitIonCalPosMosaic (ObitIonCal *in, ObitImageMosaic* mosaic, 
			  olong epoch, ObitErr* err)  
{
  olong   field, j, good, FitSize=0, prtLv, ncoef, nobs, maxQual, addSize, nZern;
  olong number, total;
  ofloat flux, FitDist, dist;
  ofloat mxcdis, maxwt, MaxRMS, MaxRMS2, rms;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  ofloat offset[2], pixel[2], corr[2], coef[17], peak, fint, dr, dd, wt;
  ObitImage *image=NULL;
  CalListElem *elem=NULL, *nelem=NULL;
  GSList  *tmp;
  CalList *calList;
  gboolean doMore;
  ofloat *x=NULL, *y=NULL, *dx=NULL, *dy=NULL, *w=NULL;
  olong *isou=NULL;
  gchar *routine = "ObitIonCalPosMosaic";

  /* Error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return ;  /* previous error? */
  g_assert(ObitIonCalIsA(in));
  g_assert(ObitImageMosaicIsA(mosaic));

  calList = (CalList*)in->calList;

  /* Control information. */
  maxwt = 10.0;
  ObitInfoListGetTest(in->info, "MaxWt",   &type, dim, &maxwt);
  flux = 1.0;
  ObitInfoListGetTest(in->info, "MinPeak", &type, dim, &flux);
  mxcdis = 1.0;
  ObitInfoListGetTest(in->info, "MaxDist", &type, dim, &mxcdis);
  maxQual = 1;
  ObitInfoListGetTest(in->info, "MaxQual",    &type, dim, &maxQual);
  FitDist = 0.0;
  ObitInfoListGetTest(in->info, "FitDist", &type, dim, &FitDist);
  FitDist /= 3600.0; /* to degrees */
  nZern = 5;
  ObitInfoListGetTest(in->info, "nZern",   &type, dim, &nZern);
  ncoef = nZern;
  prtLv = 0;
  ObitInfoListGetTest(in->info, "prtLv",   &type, dim, &prtLv);
  MaxRMS = 10.0;
  ObitInfoListGetTest(in->info, "MaxRMS",   &type, dim, &MaxRMS);
  MaxRMS /= 3600.0;  /* to deg */
  MaxRMS2 = 2.0 * MaxRMS; /* relax for finding sources */
  /* lenient RMS = 3 pixels
  MaxRMS = 3.0*fabs(mosaic->images[0]->myDesc->cdelt[0]); */

  /* First pass to get basic geometric distortions */
  /* allocate memory */
  nobs = mosaic->numberImages;
  x    = g_malloc0(nobs*sizeof(ofloat));
  y    = g_malloc0(nobs*sizeof(ofloat));
  dx   = g_malloc0(nobs*sizeof(ofloat));
  dy   = g_malloc0(nobs*sizeof(ofloat));
  w    = g_malloc0(nobs*sizeof(ofloat));
  isou = g_malloc0(nobs*sizeof(olong));

  /* Loop getting initial positions */
  nobs = 0;
  for (field=0; field<mosaic->numberImages; field++) {
    image = mosaic->images[field];
 
    /* Convert FitDist to pixels */
    if (FitDist<=0.0) FitDist = 10.0 * fabs(image->myDesc->cdelt[0]) / 3600.0;
    FitSize = (olong)((FitDist / fabs(image->myDesc->cdelt[0]))+0.5);
    FitSize = MAX (10, FitSize);

    /* Find cal info */
    tmp = calList->list;
    while (tmp!=NULL) {   /* loop 100 */
      elem = (CalListElem*)tmp->data;
      if (elem->calNo == (field+1)) break;  /* found it */
      tmp = g_slist_next(tmp);  /* next item in list */
    } /* end loop  L100: over calibrator list */

    PosImage (image, elem->pixel, flux, mxcdis, FitSize, 
	      offset, &peak, &fint, err);
    if (err->error) goto cleanup;
    
    if (peak>flux) {  /* Valid datum? */
      /* Within maximum distance? */
      dist = sqrt(elem->ZernXY[0]*elem->ZernXY[0] + elem->ZernXY[1]*elem->ZernXY[1]);
      if (dist <= mxcdis)  w[nobs] = MIN (maxwt, peak);
      else w[nobs] = 0.0;
      if (elem->qual>maxQual) w[nobs] = 0.0;  /* Only acceptable quality */
      x[nobs]    = elem->ZernXY[0];
      y[nobs]    = elem->ZernXY[1];
      dx[nobs]   = offset[0];
      dy[nobs]   = offset[1];
      isou[nobs] = nobs;
      nobs++;  /* This one OK */
      
      /* diagnostics? */
      if (prtLv>=3) {
	Obit_log_error(err, OBIT_InfoErr, 
		       "Initial E=%5d C=%5d Z= %7.4f %7.4f off=%7.2f %7.2f Peak=%6.2f int=%6.2f qual=%d",
		       epoch, elem->calNo, elem->ZernXY[0], elem->ZernXY[1], 
		       offset[0]*3600.0, offset[1]*3600.0, peak, fint, elem->qual);
      }
    } /* End if detected */
    
  } /* end loop over fields */

  /* Work out Zernike model for field - fit and edit */
  doMore = TRUE;
  while (doMore) {
    rms = IonFit1 (nobs, isou, x, y, dx, dy, w, &ncoef, coef, prtLv, err);
    if (rms<MaxRMS2) break;
    doMore = IonEdit1 (nobs, isou, x, y, dx, dy, w, MaxRMS2, ncoef, coef, prtLv, err);
    if (err->error) goto cleanup;
   }

  /* Was it happy with the results? */
  good = 0;
  for (j=0; j<nobs; j++) if (w[j]>0.0) good++;
  if ((good<=2) || (rms>MaxRMS2)) { /* need at least 3, no? - use no distortion */
    ncoef = 2;
    coef[0] = coef[1] = 0;
    FitSize = MAX (10, FitSize);
  } else { /* OK - set size of fitting region */
    FitSize = 10;
    /* Add an amount propostional to tip/tilt */
    addSize = 0.5 + (0.1 * (sqrt (coef[0]*coef[0]+coef[1]*coef[1]) / 
			    fabs(mosaic->images[0]->myDesc->cdelt[0])) + 
		     /* Add a tern for the rms residual */
		     + rms/fabs(mosaic->images[0]->myDesc->cdelt[0]));
    FitSize += addSize;
  }

  if (prtLv>=2) {
    Obit_log_error(err, OBIT_InfoErr, 
		   "%s: Initial no. cal %d rms %f Search Size %d",
		   routine, good, rms*3600.0, FitSize);
  }

 cleanup:
  if (x)    g_free(x);  x =NULL;
  if (y)    g_free(y);  y =NULL;
  if (dx)   g_free(dx); dx=NULL;
  if (dy)   g_free(dy); dy=NULL;
  if (w)    g_free(w);  w =NULL;
  if (isou) g_free(isou);
  if (err->error) Obit_traceback_msg (err, routine, image->name);

  /* Loop over images getting final positions */
  good = 0;
  for (field=0; field<mosaic->numberImages; field++) {
    image = mosaic->images[field];
 
    /* Find cal info */
    tmp = calList->list;
    while (tmp!=NULL) {   /* loop 100 */
      elem = (CalListElem*)tmp->data;
      if (elem->calNo == (field+1)) break;  /* found it */
      tmp = g_slist_next(tmp);  /* next item in list */
    } /* end loop  L100: over calibrator list */

    /* Work out expected position */
    dr = 0.0;
    dd = 0.0;
    for (j=0; j<ncoef; j++) {
      dr += coef[j]*ObitZernikeGradX(j+2, elem->ZernXY[0], elem->ZernXY[1]);
      dd += coef[j]*ObitZernikeGradY(j+2, elem->ZernXY[0], elem->ZernXY[1]);
    } 
    corr[0] = dr / image->myDesc->cdelt[0];
    corr[1] = dd / image->myDesc->cdelt[1];
    pixel[0] = elem->pixel[0] - corr[0];
    pixel[1] = elem->pixel[1] - corr[1];

    PosImage (image, pixel, flux, mxcdis, FitSize, 
	      offset, &peak, &fint, err);
    if (err->error) Obit_traceback_msg (err, routine, image->name);

    /* Adjust result for quess of position */
    offset[0] += dr;
    offset[1] += dd;
    
    if (peak>flux) {  /* Valid datum? */
      /* Within maximum distance? */
      dist = sqrt(elem->ZernXY[0]*elem->ZernXY[0] + elem->ZernXY[1]*elem->ZernXY[1]);
      if (dist <= mxcdis)  wt = MIN (maxwt, peak);
      else wt = 0.0;
      
      /* Update CalList */
      nelem = newCalListElem (elem->ra, elem->dec, elem->ZernXY, elem->shift, 
			      elem->pixel, elem->flux, offset, peak, fint, wt, 
			      elem->qual, field+1, epoch);
      CalListAdd (calList, nelem);
      good++;  /* This one OK */
      
      /* diagnostics? */
      if (prtLv>=1) {
	Obit_log_error(err, OBIT_InfoErr, 
		       "Fit E=%5d C=%5d Z= %7.4f %7.4f off=%7.2f %7.2f Peak=%6.2f int=%6.2f qual=%d",
		       epoch, elem->calNo, elem->ZernXY[0], elem->ZernXY[1], 
		       offset[0]*3600.0, offset[1]*3600.0, peak, fint, elem->qual);
      }
    } /* End if detected */
    
  } /* end loop over fields */
  
  if (prtLv>=5) {
    /*ObitMemPrint (stdout);DEBUG*/
    ObitMemSummary (&number, &total);
    Obit_log_error(err, OBIT_InfoErr, "%s: Memory use, %d entries %d MByte)", 
		 routine, number, total);
  }
 /* tell how many */
  Obit_log_error(err, OBIT_InfoErr, "%s: Fitted %d calibrator offsets", 
		 routine, good);

  /* If there aren't any - flag time */
  if (good<1) {
    CalBadTime (in, mosaic, epoch, err);
    Obit_log_error(err, OBIT_InfoWarn, "%s: CLEAN Failed this time slice", 
		   routine);
  }
} /* end of routine ObitIonCalPosMosaic */ 

/**
 * Search for a source at the given position return fitted  
 * flux density and offsets from the expected positions.  
 * A source is only accepted if the peak is within   
 * 0.5*min(imsi[0],imsi[1])  pixels of the center.  
 * \param image    Image open and plane read and attached
 * \param pixel    Expected pixel (1-rel) of centroid
 * \param minflux  Min. allowed peak flux density
 * \param mxcdis   Max. distance from center - not used here
 * \param FitSize  Full width in pixels of search window around pixel
 * \param offset   [out] Offset from expected position in cells
 * \param peak     [out] Peak flux density
 * \param fint     [out] integrated flux density
 * \param err      Error stack 
 */
static void PosImage (ObitImage *image, ofloat pixel[2], ofloat minflux, 
		      ofloat mxcdis, olong FitSize, 
		      ofloat offset[2], ofloat *peak, ofloat *fint, ObitErr *err)
{
  olong j, imsi[2];
  olong blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitIOSize IOBy;
  ofloat area, pixelOff[2];
  gchar *routine = "PosImage";

  /* Open and read image full plane */
  IOBy = OBIT_IO_byPlane;
  dim[0] = 1;
  ObitInfoListAlwaysPut (image->info, "IOBy", OBIT_long, dim, &IOBy);
  dim[0] = 7;
  for (j=0; j<IM_MAXDIM; j++) {blc[j] = 1; trc[j] = 0;}
  ObitInfoListAlwaysPut (image->info, "BLC", OBIT_long, dim, blc); 
  ObitInfoListAlwaysPut (image->info, "TRC", OBIT_long, dim, trc);
  image->extBuffer = FALSE;  /* Make sure it has buffer */
  ObitImageOpen (image, OBIT_IO_ReadOnly, err); 
  ObitImageRead (image, NULL, err);
  if (err->error) Obit_traceback_msg (err, routine, image->name);
  
  /* fit source */
  imsi[0] = imsi[1] = FitSize;
  CalPosFit (image, pixel, imsi, pixelOff, peak, fint, err);
  if (err->error) Obit_traceback_msg (err, routine, image->name);
  
  /* Close image */
  ObitImageClose (image, err);
  if (err->error) Obit_traceback_msg (err, routine, image->name);
  image->image = ObitImageUnref(image->image);   /* Free buffer */
  
  /* find anything? */
  if (*peak<minflux) { /* Nope */
    *peak = 0.0;
  } else { /* yes */
    
    /* Offset */
    offset[0] = pixelOff[0] * image->myDesc->cdelt[0];
    offset[1] = pixelOff[1] * image->myDesc->cdelt[1];
    
    /* Normalize integrated flux by beam area in pixels */
    area = 1.1331 * image->myDesc->beamMaj * image->myDesc->beamMin /  
      (fabs (image->myDesc->cdelt[0]) * fabs (image->myDesc->cdelt[1])) ;
    if (area  <  0.001) area = 1.0;
    *fint /= area;
  }
} /* end PosImage */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitIonCalClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitIonCalClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitIonCalClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitIonCalClassInfoDefFn (gpointer inClass)
{
  ObitIonCalClassInfo *theClass = (ObitIonCalClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitIonCalClassInit;
  theClass->newObit       = (newObitFP)newObitIonCal;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitIonCalClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitIonCalGetClass;
  theClass->ObitCopy      = (ObitCopyFP)ObitIonCalCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitIonCalClear;
  theClass->ObitInit      = (ObitInitFP)ObitIonCalInit;
  theClass->ObitIonCalCreate = (ObitIonCalCreateFP)ObitIonCalCreate;

} /* end ObitIonCalClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitIonCalInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitIonCal *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->thread    = newObitThread();
  in->info      = newObitInfoList(); 
  in->calList   = (gpointer)newCalList();
  in->myData    = NULL;

} /* end ObitIonCalInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitIonCal* cast to an Obit*.
 */
void ObitIonCalClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitIonCal *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  in->thread    = ObitThreadUnref(in->thread);
  in->info      = ObitInfoListUnref(in->info);
  freeCalList((CalList*)in->calList); in->calList=NULL;
  in->myData    = ObitUVUnref(in->myData);
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitIonCalClear */

/**
 * Search for a source at the given position return fitted  
 * flux density and offsets from the expected positions.  
 * A source is only accepted if the peak is within   
 * 0.5*min(imsi[0],imsi[1])  pixels of the center.  
 * Routine translated from the AIPSish CALPOS.FOR/FXPFIT  
 * \param image    Image open and plane read and attached
 * \param inPixel  Expected pixel (1-rel) of centroid
 * \param imsi     Full width of the window to search around inPixel
 * \param pixelOff Fitted offset from expected pixel (pixels)
 * \param souflx   Fitted peak flux density 
 *                 Value < 0.0 indicate source not found or 
 *                 otherwise unacceptable. 
 * \param souint   Integrated flux density in 9x9 (not normalized by 
 *                 beam area) 
 * \param err      Error stack 
 */
static void CalPosFit (ObitImage *image, ofloat inPixel[2], olong* imsi, 
		       ofloat pixelOff[2], ofloat* souflx, ofloat* souint,
		       ObitErr *err)
{
  ofloat s, dx[2], dis, maxdis;
  olong   i1, i2, j1, j2, ic, jc;
  olong blc[2], trc[2], pos[7], iX, iY, iXp, iYp, iXcen, iYcen;
  ofloat data[9][9], pixval=0.0, *pixP, sum, fblank= ObitMagicF();
  ObitFArray *fitData;
  gchar *routine = "CalPosFit";

  /* Initial values */
  *souflx = -1.0;
  *souint = -1.0;
  pixelOff[0] = 0.0;
  pixelOff[1] = 0.0;

  /* Find peak, copy data */
  ic = inPixel[0] + 0.5;  /* Closest pixel */
  jc = inPixel[1] + 0.5;
  i1 = ic - imsi[0] / 2 ;
  i2 = ic + (imsi[0] / 2) - 1;
  j1 = jc - imsi[1] / 2 ;
  j2 = jc + (imsi[1] / 2) - 1;
  i1 = MAX (MIN (i1, image->myDesc->inaxes[0]), 1);
  i2 = MIN (MAX (i2, i1), image->myDesc->inaxes[0]);
  j1 = MAX (MIN (j1, image->myDesc->inaxes[1]), 1);
  j2 = MIN (MAX (j2, j1), image->myDesc->inaxes[1]);

  blc[0] = i1-1; blc[1] = j1-1; /* to zero rel */
  trc[0] = i2-1; trc[1] = j2-1;
  fitData = ObitFArraySubArr (image->image, blc, trc, err);
  if (err->error) Obit_traceback_msg (err, routine, image->name);

  /* Peak */
  s = ObitFArrayMax (fitData, pos);
  *souflx = s;

  /* use 9x9 values around center */
  iXcen = pos[0];
  iYcen = pos[1];
  sum = 0.0;
  for (iY=0; iY<9; iY++) {
    iYp = iYcen + iY - 4;
    pos[1] = iYp;
    for (iX=0; iX<9; iX++) {
      iXp = iXcen + iX - 4;
      pos[0] = iXp;
      /* get value */
      pixP =  ObitFArrayIndex (fitData, pos);
      if (pixP) pixval = *pixP;  /* in array? */
      else pixval = fblank;
      if (pixval!=fblank) {
	data[iY][iX] = pixval;
	sum += pixval;
      } else  /* blank */
	data [iY][iX] = fblank;
    }
  }  /*  end of loop loading image data in array */

  /* Cleanup */
  fitData = ObitFArrayUnref(fitData);

  /* Fit peak in data */
  if (pfit (data, &s, dx, fblank)) {
    *souflx = -1.0; /* fit failed */
    *souint = -1.0;
    return;
  } 

  *souflx = s;   /* Fitted peak */
  *souint = sum; /* Integral is sum */

  /* get offset from expected */
  dx[0] += iXcen;  /* offset in fitData */
  dx[1] += iYcen;
  /* revert to 1-rel */
  pixelOff[0] = (-(dx[0] + blc[0] + 1 - inPixel[0]));
  pixelOff[1] = (-(dx[1] + blc[1] + 1 - inPixel[1]));

  /* Is he close enough to the center? */
  maxdis = (0.5 * MIN (imsi[0], imsi[1]))*(0.5 * MIN (imsi[0], imsi[1]));
  dis = pixelOff[0]*pixelOff[0] + pixelOff[1]*pixelOff[1];
  if (dis > maxdis) {  /* Nope */
    *souflx = -1.0;
    *souint = -1.0;
    return;
  } 

} /* end of routine CalPosFit */ 

/**
 * Searches Catalog for calibrators an image and with an estimated flux
 * density in excess of a given limit taking into account the estimated 
 * single-dish beam  
 * Adapted from the AIPSish ADNVSS in $FOURMASS/SUB/ADDFIELDS.FOR
 * \param image        image in question
 * \param Catalog      FITS AIPSVZ format catalog file name
 * \param catDisk      FITS disk number for Catalog
 * \param OutlierFlux  Minimum estimated flux density 
 * \param OutlierSI    Spectral index to use to convert catalog flux density
 * \param qual         Maximum qualifier quality code
 * \param AntDiam      Primary antenna diameter (m) [default 25]
 * \param calList      calibrator list, a linked list of calibrators
 * \param prtLv        Print level >=3 => give list
 * \param err          Error stack, returns if not empty.
 */
static void 
FindCalImage (ObitImage *image, gchar *Catalog, olong catDisk, 
	      ofloat OutlierFlux, ofloat OutlierSI, olong qual, ofloat AntDiam,
	      CalList *calList, olong prtLv, ObitErr *err) 
{
  olong count, cqual;
  odouble Freq, ra, dec, ra2000, dc2000, ra2000d, dc2000d;
  odouble xx, yy, zz, dist, refreq;
  ofloat radius, radRA, minflx, asize, alpha, pbf;
  ofloat flux, scale;
  gboolean doJ2B, doJinc, wanted;
  olong blc[IM_MAXDIM] = {1,1,1,1,1};
  olong trc[IM_MAXDIM] = {0,0,0,0,0};
  olong ver, nrows, irow;
  olong bad, calNo, ierr, epoch=0;
  ofloat shift[2], pixel[2], offset[2]={0.0,0.0}, peak=0.0, fint=0.0, wt=0.0;
  ofloat ZernXY[2];
  ObitIOCode retCode;
  ObitImageDesc *desc=NULL;
  ObitImage *VZImage=NULL;
  ObitTableVZ *VZTable=NULL;
  ObitTableVZRow *VZRow=NULL;
  CalListElem *elem=NULL;
  gchar *routine = "FindCalImage";
  
  /* error checks */
  if (err->error) return;
  g_assert(ObitImageIsA(image));

  /* get control parameters */
  minflx = OutlierFlux;
  alpha  = OutlierSI;
  asize = AntDiam; 

  /* set defaults. */
  if (asize <= 0.0)  asize  = 25.0;
  if (alpha == 0.0)  alpha  = -0.75;


  /* Is image in B1950? */
  desc = image->myDesc;
  doJ2B = (fabs(image->myDesc->equinox-1950.0) < 1.0) ||
    (fabs(desc->epoch-1950.0) < 1.0);

  /* get j2000 position to lookup in Catalog in radians */
  ra2000 = DG2RAD * desc->crval[desc->jlocr];
  dc2000 = DG2RAD * desc->crval[desc->jlocd];;
  if (doJ2B) ObitSkyGeomBtoJ (&ra2000, &dc2000);
  ra2000d = ra2000 * RAD2DG;
  dc2000d = dc2000 * RAD2DG;

  /* set crude search radius = 0.5*sqrt(2)*MAX (x_size, Y_size) */
  radius = 0.5 * sqrt (2.0) * 
    MIN (fabs(desc->cdelt[desc->jlocr]*desc->inaxes[desc->jlocr]),
	 fabs(desc->cdelt[desc->jlocd]*desc->inaxes[desc->jlocd]));
  radRA = radius / cos(dc2000);
  
  /* Observing frequency */
  if (desc->jlocf>=0) Freq = desc->crval[desc->jlocf];
  else Freq = 1.0e9;

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
  for (irow= 1; irow<=nrows; irow++) { /* loop 500 */
    /* read */
    retCode = ObitTableVZReadRow (VZTable, irow, VZRow, err);
    if (err->error) Obit_traceback_msg (err, routine, VZTable->name);
   
    /* spectral scaling of flux density */
    flux = VZRow->PeakInt * scale;

    /* position, etc */
    ra    = VZRow->Ra2000;
    dec   = VZRow->Dec2000;
    cqual = VZRow->Quality;

    /* select (crude) */
    if ((fabs(dc2000d-dec) <= radius)  && (fabs(ra2000d-ra) <= radRA) && 
	(flux >= minflx)) {
      /* separation from pointing center */
      xx = DG2RAD * ra;
      yy = DG2RAD * dec;
      zz = sin (yy) * sin (dc2000) + cos (yy) * cos (dc2000) * cos (xx-ra2000);
      zz = MIN (zz, 1.000);
      dist = acos (zz) * RAD2DG;

      if (dist>radius) continue;

      /* primary beam correction to flux density */
      if (doJinc) {
	pbf = ObitPBUtilJinc (dist, Freq, asize, 0.05);
      } else {
	pbf = ObitPBUtilPoly (dist, Freq, 0.05);
      } 
      flux *=  MAX (0.05, pbf); /* Don't trust below 5% */
      
      /* Convert position to pixel and see if in image*/
      bad = 
	ObitSkyGeomXYpix(ra, dec,
			 desc->crval[desc->jlocr], desc->crval[desc->jlocd],
			 desc->crpix[desc->jlocr], desc->crpix[desc->jlocd],
			 desc->cdelt[desc->jlocr], desc->cdelt[desc->jlocd],
			 desc->crota[desc->jlocd], &desc->ctype[desc->jlocr][4],
			 &pixel[0], &pixel[1]);
      bad = bad || (pixel[0]<1.0) || (pixel[0]>desc->inaxes[0]) || 
	(pixel[1]<1.0) || (pixel[1]>desc->inaxes[1]);
      
      /* select (fine) */
      wanted = ((flux >= minflx)  &&  (dist <= radius))  && (bad==0) && (cqual<=qual);
      if (wanted) {
	if (doJ2B) {  /* precess if necessary */
	  ra  *= DG2RAD;
	  dec *= DG2RAD;
	  ObitSkyGeomBtoJ (&ra, &dec);
	  ra  *= RAD2DG;
	  dec *= RAD2DG;
	} 

	/* get shift needed */
	ObitSkyGeomShiftXY (desc->crval[desc->jlocr], desc->crval[desc->jlocd], 
			    desc->crota[desc->jlocd], ra, dec, 
			    &shift[0], &shift[1]);

	/* Offset on Zernike plane */
	ObitSkyGeomRADec2Zern (desc->crval[desc->jlocr], desc->crval[desc->jlocd], 
			       shift[0], shift[1], &ZernXY[0], &ZernXY[1], &ierr);
	Obit_return_if_fail((ierr==0), err, 
			    "%s: Error projecting onto Zernike Unit circle", routine);

	/* Add to CalList */
	calNo = count;
	elem = newCalListElem (ra, dec, ZernXY, shift, pixel, flux, offset, 
			        peak, fint, wt, cqual, calNo, epoch);
	CalListAdd (calList, elem);
	count++;
      } /* end if wanted */ 
    } /* end crude selection */
  } /* end loop over table */

  /* tell how many */
  Obit_log_error(err, OBIT_InfoErr, "%s: Found %d calibrators", routine, count);

  /* Close up */
  ObitImageClose(VZImage, err);
  retCode = ObitTableVZClose(VZTable, err);
  if (err->error) Obit_traceback_msg (err, routine, VZTable->name);
  VZImage->image = ObitImageUnref(VZImage->image);   /* Free buffer */

  /* Diagnostics? */
  if (prtLv>=3) {
    CalListPrint(calList, stdout); 
  }
  
  /* clean up */
  VZImage = ObitImageUnref(VZImage);
  VZTable = ObitTableUnref(VZTable);
  VZRow   = ObitTableRowUnref(VZRow);
  
} /* end of routine FindCalImage */ 

/**
 * Edit Zernike data.
 * Edit data in attempt to get rms residual under MaxRMS
 * by zeroing weight of most discrepant datum which must have
 * at least 1.5 x the average contribution to the variance.
 * Also checks that there is enough data to overdetermine the model.
 * Routine adapted from the AIPSish IONCAL.FOR/IONEDT 
 * \param nobs Number of observations  
 * \param isou  Source number per obs.  (0-rel)
 * \param x     Offset in field X (RA) on unit Zernike circle  per source  
 * \param y     Offset in field Y (Dec) on unit Zernike circle per source  
 * \param dx    Apparent RA position shifts on the Zernike plane (deg) 
 * \param dy    Apparent dec position shifts on the Zernike plane (deg) 
 * \param w     Weight per source, may be set to zero,
 * \param MaxRMS Target RMS residual (units of  dx, dy)
 * \param ncoef number of coefficients in coef
 * \param coef  Zernike coefficients to use 
 * \param prtLv Print level >=2 tell about editing
 * \param err   Error stack
 * \return True if data was edited 
 */
static gboolean
IonEdit1 (gint nobs, olong* isou, ofloat* x, ofloat* y, ofloat* dx, ofloat* dy, 
	  ofloat* w, ofloat MaxRMS, olong ncoef, ofloat* coef, olong prtLv, 
	  ObitErr* err) 
{
  olong      i, j, is, rmscnt, ibad;
  gboolean out = FALSE;
  ofloat    dr, dd, rx, ry, rms;
  ofloat    *gdx=NULL, *gdy=NULL, var, *ResidV=NULL, baddest;
  gchar *routine = "IonEdit1";

  if (nobs<=0) return FALSE;  /* Something to do? */

  /* Allocate memory */
  gdx    = g_malloc(nobs*ncoef*sizeof(ofloat));
  gdy    = g_malloc(nobs*ncoef*sizeof(ofloat));
  ResidV = g_malloc(nobs*sizeof(ofloat));

  /* Compute Zernike gradients */
  for (i=0; i<nobs; i++) {
    is = isou[i];
    for (j=0; j<ncoef; j++) {
      gdx[j*nobs+i] = ObitZernikeGradX(j+2, x[is], y[is]);
      gdy[j*nobs+i] = ObitZernikeGradY(j+2, x[is], y[is]);
    } /* end loop over coefficients */
  } /* end loop over sources */

  /* get RMS and residuals */
  rms = 0.0;
  rmscnt = 0;
  for (i=0; i<nobs; i++) {
    if (w[i] > 0.0) {
      dr = 0.0;
      dd = 0.0;
      for (j=0; j<ncoef; j++) {
	dr += coef[j]*gdx[j*nobs+i];
	dd += coef[j]*gdy[j*nobs+i];
      } 
      rx = dx[i] - dr;
      ry = dy[i] - dd;
      var = rx*rx + ry*ry;
      ResidV[i] = var;
      rms += var;
      rmscnt++;
    } 
  } 

  /* must have enough data to be over determined */
  if ((2*rmscnt) < (ncoef+1)) goto cleanup;
  var = rms/rmscnt;
  rms = sqrt (var);

  /* Is the current version OK? */
  if (rms <= MaxRMS) goto cleanup;

  /* find the most discrepent */
  ibad = -1;
  baddest = 0.0;
   for (i=0; i<nobs; i++) {
     if ((w[i]>0.0) && (ResidV[i]> baddest)) {
       baddest = ResidV[i];
       ibad = i;
     }
   }

   /* Is the bad one at least 1.5x average variance? */
   if (baddest<(1.5*var)) goto cleanup;

   w[ibad] = 0.0;  /* Flag it */
   out = TRUE;     /* Flagged something */

   /* Diagnostics */
   if (prtLv>=2) {
     Obit_log_error(err, OBIT_InfoErr, 
		    "%s: Flag obs %d variance=%f rms=%f",
		    routine, ibad+1, baddest*3600.0, rms*3600.0);
   }

  /* Deallocate memory */
 cleanup: 
   g_free(ResidV);
   g_free(gdx);
   g_free(gdy);

  return out;
} /* end IonEdit1 */ 

/**
 * Get timerange of next interval  
 *  Routine translated from the AIPSish IONCAL.FOR/ICNXTI  
 * \param NXtab    Index table 
 * \param solInt   Solution interval in min. 
 * \param timer    In: previous timerange (days), -1 => none. 
 *                 Out: next time range 
 * \param suba     [out] Subarray number 
 * \param err    Error code: 0 => ok, -1 = no more. 
 * \return TRUE if there is more data to process
 */
static gboolean NextTime (ObitTableNX *NXtab, ofloat solInt, 
			  ofloat timer[2], olong* suba, ObitErr* err) 
{
  gboolean out = FALSE;
  olong  numNX, npiece;
  olong irow, NXrow=0;
  ofloat si, tbeg, tend, vtend=0, vtbeg=0, ts, te;
  ObitTableNXRow *row;
  gchar *routine = "NextTime";

  /* Error checks */
  if (err->error) return FALSE;  /* previous error? */
  g_assert(ObitTableNXIsA(NXtab));
  if (err->error) return FALSE;  /* previous error? */

  /* Create table row */
  row = newObitTableNXRow (NXtab);

  /* Open table */
  ObitTableNXOpen(NXtab, OBIT_IO_ReadOnly, err);
  if (err->error) goto cleanup;
  numNX = NXtab->myDesc->nrow;

  /* Save initial timerange  */
  tbeg = timer[0];
  tend = timer[1];
  si = solInt / 1440.0;
  /* Find next scan */
  for (irow= 1; irow<= numNX; irow++) { /* loop 200 */
    NXrow = irow;
    ObitTableNXReadRow (NXtab, NXrow, row, err);
    if (err->error) goto cleanup;
    if (row->status<0) continue; /* valid row? */
    
    /* time range of this scan */
    vtbeg = row->Time - 0.5 * row->TimeI;
    vtend = row->Time + 0.5 * row->TimeI;

    /* Ends before previous start? */
    if (vtend < tbeg) continue;

    /* Is previous end in this scan? */
    if ((tend > vtbeg)  &&  (tend < vtend)  && 
	/* Anything left? */
	((tend+0.5*si) < vtend)) break;

    /* Past last scan, take this one */
    if (tend < vtbeg) break;
  } /* end loop  L200: */

  out = TRUE;  /* got here without failing */

  /* Are we done? */
  if ((NXrow == numNX)  &&  ((tend+0.5*si) > vtend)) {
    out = FALSE;
    goto cleanup;
  } 

  /* Divide scan into equal pieces */
  npiece = 0.5 + row->TimeI / si;
  npiece = MAX (1, npiece);
  si = row->TimeI / npiece;

  /* which one is this? */
  ts = MAX (vtbeg, tend);
  te = ts + si;

  /* Set output */
  timer[0] = ts;
  timer[1] = te;
  *suba = MAX (1, row->SubA);

  /* Found it, close table */
 cleanup: ObitTableNXClose(NXtab, err);
  row = ObitTableNXRowUnref (row);
  if (err->error) Obit_traceback_val (err, routine, NXtab->name, out);
  return out;
} /* end of routine NextTime */ 


/**
 * Generate a Calibrator list from an ImageMosaic
 * This assumes that the expected calibrator position is the
 * reference pixel in the image.
 * Calibrators added to CalList as epoch 0;
 * \param in           IonCal object (used for catalog info )
 * \param mosaic       Mosaic to use
 * \param calList      calibrator list, a linked list of calibrators
 * \param prtLv        Print level >=1, give cal info  >=3 => give full list
 * \param err          Error stack, returns if not empty.
 */
static void 
FindCalMosaic (ObitIonCal *in, ObitImageMosaic *mosaic, 
	       CalList *calList, olong prtLv, ObitErr *err)
{
  olong field, calNo, cqual, ierr, epoch=0;
  ofloat  flux, peak=0.0, fint=0.0, wt=0.0;
  ofloat ZernXY[2], shift[2], pixel[2], offset[2]={0.0,0.0};
  odouble ra, dec, raPnt, decPnt;
  ObitImageDesc *desc=NULL;
  CalListElem *elem=NULL;
  gchar *routine = "FindCalMosaic";
  
  /* error checks */
  if (err->error) return;
  g_assert(ObitImageMosaicIsA(mosaic));

  /* Get pointing position for offsets */
  desc = mosaic->images[0]->myDesc;
  ObitImageDescGetPoint (desc, &raPnt, &decPnt);

  /* Loop over images */
  for (field=0; field<mosaic->numberImages; field++) {
    calNo = field+1;

    /* Info from Mosaic */
    desc = mosaic->images[field]->myDesc;
    ra  = desc->crval[desc->jlocr];
    dec = desc->crval[desc->jlocd];
    pixel[0] = desc->crpix[desc->jlocr];
    pixel[1] = desc->crpix[desc->jlocd];

    /*  Some items to be filled in from catalog */
    cqual = -1;
    flux  = -1.0;

    /* get shift needed */
    ObitSkyGeomShiftXY (raPnt, decPnt, desc->crota[desc->jlocd], ra, dec, 
			&shift[0], &shift[1]);
    
    /* Offset on Zernike plane */
    ObitSkyGeomRADec2Zern ( raPnt, decPnt, shift[0], shift[1],
			    &ZernXY[0], &ZernXY[1], &ierr);
    Obit_return_if_fail((ierr==0), err,
			"%s: Error projecting onto Zernike Unit circle", routine);
    
    /* Add to CalList */
    elem = newCalListElem (ra, dec, ZernXY, shift, pixel, flux, offset, 
			   peak, fint, wt, cqual, calNo, epoch);
    CalListAdd (calList, elem);
  } /* end loop over fields */

  /* Get calibrator information from calList */
  LookupCals (in, mosaic, calList, prtLv, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Diagnostics? */
  if (prtLv>=3) {
    CalListPrint(calList, stdout); 
  }
} /* end FindCalMosaic */

/**
 * Looks up calibrators in calList in the catalog specified (or implied)
 * by in->info and fills in catalog information.
 * Adapted from the AIPSish ADNVSS in $FOURMASS/SUB/ADDFIELDS.FOR
 * \param in           IonCal with catalog information
 * \param mosaic       ImageMosaic to describing calList fields
 * \param calList      calibrator list, a linked list of calibrators
 * \param prtLv        Print level >=1 => give list
 * \param err          Error stack, returns if not empty.
 */
static void 
LookupCals (ObitIonCal *in, ObitImageMosaic *mosaic, CalList *calList, 
	    olong prtLv, ObitErr *err)
{
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  gchar Catalog[257];
  olong catDisk;

  olong count, cqual;
  odouble Freq, ra, dec, rar, decr, raPntr, decPntr;
  odouble xx, yy, zz, dist, refreq;
  ofloat radius, radRA, asize, alpha, pbf;
  ofloat equinox, flux, scale;
  gboolean doJ2B, doJinc;
  olong blc[IM_MAXDIM] = {1,1,1,1,1};
  olong trc[IM_MAXDIM] = {0,0,0,0,0};
  olong ver, nrows, irow;
  ObitIOCode retCode;
  ObitImage *VZImage=NULL;
  ObitTableVZ *VZTable=NULL;
  ObitTableVZRow *VZRow=NULL;
  ObitImageDesc *desc=NULL;
  GSList *tmp;
  CalListElem *elem=NULL;
  gchar *routine = "LookupCals";
  
  /* error checks */
  if (err->error) return;
  g_assert(ObitImageMosaicIsA(mosaic));

  /* get control parameters */
  /* Control information. */
  sprintf (Catalog, "Default");
  ObitInfoListGetTest(in->info, "Catalog", &type, dim,  Catalog);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  if (!strncmp(Catalog, "    ", 4)) sprintf (Catalog, "Default");
  if (!strncmp(Catalog, "Default", 7)) sprintf (Catalog, "NVSSVZ.FIT");
  catDisk = 1;
  ObitInfoListGetTest(in->info, "catDisk",   &type, dim, &catDisk);
  alpha = -0.75;
  ObitInfoListGetTest(in->info, "OutlierSI", &type, dim, &alpha);
  asize = 25.0;
  ObitInfoListGetTest(in->info, "AntDiam", &type, dim, &asize);

  /* set defaults. */
  if (asize <= 0.0)  asize  = 25.0;
  if (alpha == 0.0)  alpha  = -0.75;

  /* set crude search radius = 0.5*sqrt(2)*MAX (x_size, Y_size) */
  desc = mosaic->images[0]->myDesc;
  radius = 0.5 * sqrt (2.0) * 
    MIN (fabs(desc->cdelt[desc->jlocr]*desc->inaxes[desc->jlocr]),
	 fabs(desc->cdelt[desc->jlocd]*desc->inaxes[desc->jlocd]));
  radRA = radius / cos(desc->crval[desc->jlocd]*DG2RAD);
  
  /* Get pointing position for offsets */
  ObitImageDescGetPoint (desc, &raPntr, &decPntr);
  raPntr  *= DG2RAD;
  decPntr *= DG2RAD;

  /* Observing frequency */
  if (desc->jlocf>=0) Freq = desc->crval[desc->jlocf];
  else Freq = 1.0e9;

  /* which beam model to use */
  doJinc = (Freq >= 1.0e9);

  /* Need precession? */
  equinox = in->myData->myDesc->equinox;  /* Clear up confusion in AIPS */
  if (equinox<1.0) equinox = in->myData->myDesc->epoch;
  doJ2B = (equinox!=2000.0) ;  /* need to precess? */

  /* Open Catalog (VZ table on an image) */
  VZImage = newObitImage("Catalog image");
  ObitImageSetFITS(VZImage, OBIT_IO_byPlane, catDisk, Catalog, blc, trc, err);
  if (err->error) goto cleanup;

  /* Open to fully instantiate */
  ObitImageOpen(VZImage, OBIT_IO_ReadOnly, err);
  if (err->error) goto cleanup;

  /* Now get VZ table */
  ver = 1;
  VZTable =  newObitTableVZValue("Catalog table", (ObitData*)VZImage, &ver, 
				 OBIT_IO_ReadOnly, err);
  ObitTableVZOpen(VZTable, OBIT_IO_ReadOnly, err);
  VZRow =  newObitTableVZRow (VZTable);  /* Table row */
  if (err->error) goto cleanup;

  /* Get table info */
  refreq = VZTable->refFreq;
  nrows  = VZTable->myDesc->nrow;

  /* frequency scaling */
  scale = pow ((Freq / refreq), alpha);

  /* loop through table */
  count = 0;
  for (irow= 1; irow<=nrows; irow++) { /* loop 500 */
    /* read */
    retCode = ObitTableVZReadRow (VZTable, irow, VZRow, err);
    if (err->error) goto cleanup;
   
    /* spectral scaling of flux density */
    flux = VZRow->PeakInt * scale;

    /* position, etc */
    ra    = VZRow->Ra2000;
    dec   = VZRow->Dec2000;
    /* Need in 1950? */
    if (doJ2B) ObitSkyGeomJtoB (&ra, &dec);
    rar   = ra * DG2RAD;
    decr  = dec * DG2RAD;
    cqual = VZRow->Quality;

    /* Loop over calList */
    tmp = calList->list;
    while (tmp!=NULL) {
      elem = (CalListElem*)tmp->data;
      
      /* select (crude) */
      if ((fabs(dec-elem->dec) <= radius)  && (fabs(ra-elem->ra) <= radRA)) {
	/* separation from pointing center */
	xx = DG2RAD * elem->ra;
	yy = DG2RAD * elem->dec;
	zz = sin (yy) * sin (decr) + cos (yy) * cos (decr) * cos (xx-rar);
	zz = MIN (zz, 1.000);
	dist = acos (zz) * RAD2DG;
	
	/* Must be within a pixel */
	if (dist<fabs(desc->cdelt[desc->jlocr])) {
	  
	  /* Distance from pointing center */
	  zz = sin (yy) * sin (decPntr) + cos (yy) * cos (decPntr) * cos (xx-raPntr);
	  zz = MIN (zz, 1.000);
	  dist = acos (zz) * RAD2DG;

	  /* primary beam correction to flux density */
	  if (doJinc) {
	    pbf = ObitPBUtilJinc (dist, Freq, asize, 0.05);
	  } else {
	    pbf = ObitPBUtilPoly (dist, Freq, 0.05);
	  } 
	  flux *= MAX (0.05, pbf); /* Don't trust below 5% */
	  
	  /* update CalList */
	  CalListElemUpdate (elem, elem->ra, elem->dec, elem->ZernXY, elem->shift, 
			     elem->pixel, flux, elem->offset, 
			     elem->peak, elem->fint, elem->wt, cqual, 
			     elem->calNo, elem->epoch);
	  count++;

	  /* Diagnostics? */
	  if (prtLv>=1) {
	    Obit_log_error(err, OBIT_InfoErr, 
			   "Cal %4d RA=%8.5f Dec=%9.5f Z= %7.4f %7.4f Flux=%6.2f qual=%d",
			   elem->calNo,  elem->ra, elem->dec, elem->ZernXY[0], elem->ZernXY[1], 
			  flux, cqual);
	  }
	} /* End update entry */
      } /* end crude selection */
      tmp = g_slist_next(tmp);  /* next item in list */
    } /* end loop over calList */
  } /* end loop over table */

  /* tell how many */
  Obit_log_error(err, OBIT_InfoErr, "%s: Found %d calibrators", routine, count);

  /* Close up */
  ObitImageClose(VZImage, err);
  retCode = ObitTableVZClose(VZTable, err);
  if (err->error) goto cleanup;
  /* clean up */
 cleanup: VZImage = ObitImageUnref(VZImage);
  VZTable = ObitTableUnref(VZTable);
  VZRow   = ObitTableRowUnref(VZRow);
  if (err->error) Obit_traceback_msg (err, routine, VZTable->name);
  
} /* end of routine LookupCals */ 

/**
 * Filter position offset solutions and fit Zernike models to each epoch.
 * Results written into NI table which is returned.
 * Routine adapted from the AIPSish IONCAL.FOR/IONFIT 
 * \param in   IonCal object 
 *   Control parameters are on the info member.
 * \li "FitDist"  OBIT_int   (1,1,1) dist, from expected location to search 
 *                                   asec [10 pixels]
 * \li "MinPeak"  OBIT_float (1,1,1) Min. acceptable image peak (Jy) [1.0]
 * \li "MaxDist"  OBIT_float (1,1,1) Max. distance (deg/10) to accept calibrator [1.]
 * \li "MaxWt"    OBIT_float (1,1,1) Max. weight [10.0]
 * \li prtLv      OBIT_int   (1,1,1) Print level >=1 => give fits
 * \param calList Calibrator list each element of which has:
 * \li ra      position RA (deg) of calibrators
 * \li dec     position Dec (deg) of calibrators
 * \li shift   Offset in field [x,y] on unit Zernike circle of position
 * \li pixel   Expected pixel in reference image
 * \li flux    Estimated catalog flux density 
 * \li offset  Measured offset [x,y] from expected position (deg)
 * \li peak    Measured peak flux density (image units)
 * \li fint    Measured integrated flux density (image units)
 * \li wt      Determined weight
 * \li qual    Catalog quality code
 * \param ncoef      Maximum number of coefficients to fit
 * \param nEpoch     Number of epochs
 * \param timeEpoch  Center time of each Epoch (day)
 * \param timeIEpoch Time interval of each Epoch (day)
 * \param subaEpoch  Epoch subarray
 * \param refFreq    Reference Frequency for NI table (Hz)
 * \param seeing     [out] RMS residual to final fits (asec).
 * \param err        Error stack
 */
static ObitTableNI* 
FitAll (ObitIonCal *in, olong ncoef, olong nEpoch, 
	ofloat *timeEpoch, ofloat *timeIEpoch, olong *subaEpoch,
	odouble refFreq, ofloat *seeing, ObitErr* err)  
{
  ObitTableNI* out = NULL;
  ObitUV *inUV=NULL;
  CalListElem *elem=NULL;
  GSList  *tmp;
  CalList *calList;
  ofloat maxWt, MaxRMS, MinRat, MinPeak;
  gboolean doINEdit;
  olong prtLv, nobs, nTime, nsou, maxQual, is;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  olong *isou=NULL, *iTime=NULL, *flqual=NULL;
  ofloat *x=NULL, *y=NULL, *xoff=NULL, *yoff=NULL, 
    *flux=NULL, *wt=NULL, *sint=NULL;
  gchar *routine = "FitAll";

  /* Error checks */
  if (err->error) return NULL;  /* previous error? */
  g_assert(ObitIonCalIsA(in));

  *seeing = -1.0;
  calList = (CalList*)in->calList;
  inUV = in->myData;

  /* Control information. */
  maxWt = 10.0;
  ObitInfoListGetTest(in->info, "MaxWt",   &type, dim, &maxWt);
  prtLv = 0;
  ObitInfoListGetTest(in->info, "prtLv",   &type, dim, &prtLv);
  MaxRMS = 10;
  ObitInfoListGetTest(in->info, "MaxRMS", &type, dim, &MaxRMS);
  MaxRMS /= 3600.0;  /* to degrees */
  doINEdit = TRUE;
  ObitInfoListGetTest(in->info, "doINEdit", &type, dim, &doINEdit);
  maxQual = 1;
  ObitInfoListGetTest(in->info, "MaxQual",   &type, dim, &maxQual);
  MinRat = 0.1;
  ObitInfoListGetTest(in->info, "MinRat", &type, dim, &MinRat);
  MinPeak = 1.0;
  ObitInfoListGetTest(in->info, "MinPeak", &type, dim, &MinPeak);

  /* Count number of sources, number of times and number of obs */
  nobs = nTime = nsou = 0;
  tmp = calList->list;
  while (tmp!=NULL) {   /* loop 100 */
    elem = (CalListElem*)tmp->data;
    if (elem->epoch > 0) {
      nobs++;
      nTime = MAX (nTime, elem->epoch);
      nsou  = MAX (nsou,  elem->calNo);
    }
    /* if find epoch > 0 then quit loop */
    tmp = g_slist_next(tmp);  /* next item in list */
 } /* end loop  L100: over calibrator*/

  /* Create arrays */
  isou   = g_malloc0(nobs*sizeof(olong));
  iTime  = g_malloc0(nobs*sizeof(olong));
  flqual = g_malloc0(nobs*sizeof(olong));
  x      = g_malloc0((nsou+1)*sizeof(ofloat));
  y      = g_malloc0((nsou+1)*sizeof(ofloat));
  xoff   = g_malloc0(nobs*sizeof(ofloat));
  yoff   = g_malloc0(nobs*sizeof(ofloat));
  flux   = g_malloc0(nobs*sizeof(ofloat));
  wt     = g_malloc0(nobs*sizeof(ofloat));
  sint   = g_malloc0(nobs*sizeof(ofloat));

  /* Extract data from calList */
  nobs = 0;
  tmp = calList->list;
  while (tmp!=NULL) {   /* loop 100 */
    elem = (CalListElem*)tmp->data;
    if (elem->epoch > 0) {
      isou[nobs]   = elem->calNo-1;
      iTime[nobs]  = elem->epoch-1;
      flqual[nobs] = elem->qual;
      is           = MAX (0, MIN (isou[nobs], nsou));
      x[is]        = elem->ZernXY[0];
      y[is]        = elem->ZernXY[1];
      xoff[nobs]   = elem->offset[0];
      yoff[nobs]   = elem->offset[1];
      flux [nobs]  = elem->peak;
      sint[nobs]   = elem->fint;
      /* Impose maximum Quality */
      if (elem->qual<=maxQual) wt[nobs] = elem->wt;
      else wt[nobs] = 0.0;
      wt[nobs] = MIN (wt[nobs], maxWt);
      /* Impose minimum flux density */
      if (flux[nobs]<MinPeak) wt[nobs] = 0.0;
      nobs++;
    }
    /* if find epoch > 0 then quit loop */
    tmp = g_slist_next(tmp);  /* next item in list */
 } /* end loop  L100: over calibrator*/


  /* do fitting */
  out = doIonFitAll (inUV, MaxRMS, MinRat, doINEdit, 
		     ncoef, nsou, isou, x, y, nEpoch, 
		     timeEpoch, timeIEpoch, subaEpoch, 
		     nobs, iTime, xoff, yoff, flux, wt, flqual, sint, 
		     prtLv, refFreq, seeing, err);

  /* deallocate arrays */
  /*    cleanup:*/
  if (isou) g_free(isou);
  if (iTime) g_free(iTime);
  if (flqual) g_free(flqual);
  if (x) g_free(x);
  if (y) g_free(y);
  if (xoff) g_free(xoff);
  if (yoff) g_free(yoff);
  if (flux) g_free(flux);
  if (wt) g_free(wt);
  if (sint) g_free(sint);
  if (err->error) Obit_traceback_val (err, routine, in->name, out);
  return out;
} /* end of routine FitAll  */ 

/**
 * Weighted fit ionospheric model to data using Zernike polynomials.  
 * Also filters data on basis of RMS residual if RMS exceeds MaxRMS.  
 * Writes values to the output IoN (NI) table.  
 * Ignores sources further than 10 degrees from the pointing  
 * Position are all referred to a tangent plane at the  
 * pointing center, this is referred to as the Zernike plane as this  
 * plane will be used to fit the phase screen.  The "Zernike Unit  
 * Circle" defines the phase screen.  Source position offsets are in  
 * the x and y (RA, Dec) as defined by this plane.  
 * Routine translated from the AIPSish IONCAL.FOR/IONFIT  
 * Input:  
 * \param inUV     UV data to which output NI table is to be attached.
 * \param MaxRMS   Max. allowable RMS before filtering data. (deg) 
 * \param MinRat   Minimum acceptable ratio to average flux 
 * \param doINEdit if true flag solutions for which the seeing residual 
 *                 could not be determined or exceeds MAXRMS 
 * \param ncoef    Number of Zernike coeffients to solve for. 
 * \param nsou     Number of sources 
 * \param isou     Source numbers (0-rel)
 * \param x        X Offset (RA) in field on unit Zernike circle per source  
 * \param y        Y Offset (Dec) in field on unit Zernike circle per source  
 * \param nTime    Number of time intervals 
 * \param Time     Center time of each integration 
 * \param TimeI    Time Interval of each integration 
 * \param isuba    Subarray number per time 
 * \param n        Number of data points (observations)
 * \param iTime    Time interval number per observation.(0-rel)
 * \param xoff     Apparent RA position shifts on the Zernike plane (deg) 
 *                 per observation
 * \param yoff     Apparent dec position shifts on the Zernike plane (deg) 
 *                 per observation
 * \param flux     Array of measured peak flux densities 
 *                 per observation
 * \param wt       Array of weights,  per observation
 * \param flqual   Field quality (crowding) code. Used in determining 
 *                 if the apparent position of a source can be moved. 
 *                 per observation
 * \param sint     Array of measured integrated flux densities 
 *                 per observation
 * \param prtLv    Print level >=1 => give fitting diagnostics
 * \param refFreq    Reference Frequency for NI table (Hz)
 * \param totRMS   [out] Total residual RMS seeing in asec
 * \param err      Error stack
 * \return ObitTableNI into which values were written (ver 1)
 */
 static ObitTableNI* 
 doIonFitAll (ObitUV *inUV, ofloat MaxRMS, ofloat MinRat, gboolean doINEdit, 
	      olong ncoef, olong nsou, olong* isou, ofloat* x, ofloat* y, 
	      olong nTime, ofloat* Time, ofloat* TimeI, olong* isuba, 
	      olong n, olong* iTime, ofloat* xoff, ofloat* yoff, 
	      ofloat* flux, ofloat *wt, olong* flqual, ofloat* sint, 
	      olong prtLv, odouble refFreq, ofloat* totRMS, ObitErr* err) 
{
  ObitTableNI *outTab = NULL;
  ObitTableNIRow *NIrow=NULL;
  olong     i, j, ntoss;
  ofloat   xy, timer[2], fblank = ObitMagicF();
  gboolean redo;
  olong     *mcoef=NULL, ngood, nbad, maxcoef;
  ofloat   **coef=NULL, *soffx=NULL, *soffy=NULL, *rms=NULL;
  gboolean *gotsou=NULL, *fitsou=NULL, allFlag, OK;
  olong ver, rowNo;
  gchar msgtxt[80];
  gchar *routine = "doIonFitAll";

  /* Error checks */
  if (err->error) return outTab;  /* previous error? */
  g_assert(ObitUVIsA(inUV));

  /* Something to do? */
  if ((n<1) || (nsou<1)) return outTab;

  /* allocate arrays */
  mcoef  = g_malloc0(nTime*sizeof(olong));
  rms    = g_malloc0((nTime+1)*sizeof(ofloat));
  soffx  = g_malloc0(nsou*sizeof(ofloat));
  soffy  = g_malloc0(nsou*sizeof(ofloat));
  gotsou = g_malloc0(nsou*sizeof(gboolean));
  fitsou = g_malloc0(nsou*sizeof(gboolean));
  coef   = g_malloc0(nTime*sizeof(ofloat*));
  for (i=0; i<nTime; i++) coef[i] = g_malloc0(ncoef*sizeof(ofloat));

  /* Diagnostics - show input data */
  if (prtLv>=3) {
    /* Source info */ 
       for (i=0; i<nsou; i++) {
       Obit_log_error(err, OBIT_InfoErr, 
       "Source %4d Zernike offset = %10.6f %10.6f", i+1, x[i], y[i]);
       }
    /* Epoch info */
    for (i=0; i<nTime; i++) {
      Obit_log_error(err, OBIT_InfoErr, 
		     "Time %4d Time %10.6f hr dTime%10.6f min Subarray %4d", 
		     i+1, Time[i]*24.0, TimeI[i]*1440.0, isuba[i]);
    }
    /* Measurement data */
    for (i=0; i<n; i++) {
      Obit_log_error(err, OBIT_InfoErr, 
		     "Obs %4d Time %5d Sou %5d Offset %6.2f %6.2f Peak %6.2f Int %6.2f qual %2d", 
		     i+1, iTime[i]+1, isou[i]+1, xoff[i]*3600.0, yoff[i]*3600.0, flux[i], sint[i], flqual[i]);
    }
    ObitErrLog(err);
  } /* end diagnostics */


  allFlag = TRUE;
  for (j=0; j<nsou; j++) { /* loop 50 */
    /* Diagnostics */
    if (prtLv>=3) {
      Obit_log_error(err, OBIT_InfoErr, 
		     "%s: S=%4d x=%10.5f y=%10.5f", routine, j+1, x[j], y[j]);
    }

    /* Flag data for sources further than 10 deg (1 on unit circle) */
    xy = sqrt (x[j]*x[j] + y[j]*y[j]);
    if (xy > 1.0) {
      for (i=0; i<n; i++) { /* loop 40 */
	if (isou[i] == (j+1)) flux[i] = 0.0;
      } /* end loop  L40:  */;
    } 
  } /* end loop  L50:  */

  /* further constraints on weights */
  for (i=0; i<n; i++) { /* loop 60 */
    /* Toss data if integrated value less than a third the peak */
    if (sint[i]  <  (0.35*flux[i])) wt[i] = 0.0;
    /* ... more than thrice the peak */
    if (sint[i]  >  (3.0*flux[i])) wt[i] = 0.0;
  } /* end loop  L60:  */

  /* Init Ionospheric model */
  OK = initIonModel (n, nsou, nTime, ncoef, isou, iTime,  x, y, xoff, yoff, wt, 
		     flqual, mcoef, gotsou, fitsou, soffx, soffy, coef, prtLv, err);
  if (err->error) goto cleanup;
  if (!OK) {allFlag=TRUE; goto AllBad;}  /* Fell flat */

  /* How many sources actually have  data?  */
  ngood = 0;
  for (i=0; i<nsou; i++) { /* loop 70 */
    if (gotsou[i]) ngood++;
  } /* end loop  L70:  */
  
  /* Fit Ionospheric model if enough  data  */
  if (ngood > 3)
    FitIonSeries (n, nsou, nTime, ncoef, isou, iTime,  x, y, xoff, yoff, wt,  
		  mcoef, gotsou, fitsou, soffx, soffy, coef, rms, prtLv, err);
  if (err->error) goto cleanup;

  /* If only one source dummy rms at 1 asec */
  if (nsou==1) {
    for (i=0; i<nTime; i++) rms[i] = 1.0/ 3600.0;
  }

  /* Edit Ionospheric model */
  if (doINEdit  &&  (ngood > 3)) {
        redo = IonEditSeries (n, nsou, nTime, ncoef, isou, iTime,  x, y, 
			      xoff, yoff, wt, MaxRMS, mcoef, 
			      soffx, soffy, coef, prtLv, err);
	if (err->error) goto cleanup; 

	/* Edit By amplitude */
        redo = redo || IonEditAmp (n, nsou, nTime, isou, iTime, MinRat,  flux, mcoef, wt, 
				   prtLv, err);
	if (err->error) goto cleanup;
	
	/* Redo Ionospheric model if data edited */
	if (redo) {
	  OK = initIonModel (n, nsou, nTime, ncoef, isou, iTime,  x, y, xoff, yoff, 
			     wt, flqual, mcoef, gotsou, fitsou, soffx, soffy, coef, 
			     prtLv, err);
	  if (err->error) goto cleanup;
	  if (!OK) {allFlag = TRUE; goto AllBad;}
	} 
	
	/* Refit Ionospheric model */
	FitIonSeries (n, nsou, nTime, ncoef, isou, iTime,  x, y, xoff, yoff, 
		      wt,  mcoef, gotsou, fitsou, soffx, soffy, coef, rms, prtLv, err);
	if (err->error) goto cleanup;
  } /* end editing */

  /* Only take solutions with the maximum number of terms */
  maxcoef = -1;
  ntoss = nbad = 0;
  for (j=0; j<nTime; j++) maxcoef = MAX (maxcoef, mcoef[j]);
 
  /* Create outputNI table */
  ver = 1;
  outTab = newObitTableNIValue ("IonCal NI table", (ObitData*)inUV, &ver, 
			       OBIT_IO_WriteOnly, ncoef, err);
  if (err->error) goto cleanup;
  NIrow = newObitTableNIRow (outTab);       /* Row structure */
  NIrow->antNo  = 0;  /* No antenna specific */
  NIrow->SourId = 0;  /* no source ID */
  
  /* Clear any existing rows */
  ObitTableClearRows ((ObitTable*)outTab, err);
  if (err->error) goto cleanup;

  /* Open output NI table */
  ObitTableNIOpen (outTab, OBIT_IO_WriteOnly, err);
  ObitTableNISetRow (outTab, NIrow, err);   /* Attach row for writing */
  if (err->error) goto cleanup;

  /* Header keywords */
  outTab->heightIon = 1.0e7;  /* Height of ionosphere = large */ 
  outTab->refFreq = refFreq;  /* image reference freq */

  /* Loop over times writing results */
  for (i=0; i<nTime; i++) { /* loop 200 */
    /* Set IN table entry bad if mcoef[i] = 0 */
    if ((mcoef[i]<maxcoef) && (mcoef[i]>0)) ntoss++;
    if ((mcoef[i] >= maxcoef) && (rms[i] > 0.0)) {
      NIrow->weight = 1.0 / MAX (0.001, rms[i]);
    } else {
      nbad++;  /* Count the bad */
      NIrow->weight = 0.0;
      rms[i] = -1.0;
      for (j=0; j<ncoef; j++) coef[i][j] = fblank;
    } 

    /* Fill row structure */
    NIrow->Time  = Time[i];
    NIrow->TimeI = TimeI[i];
    NIrow->SubA  = isuba[i];
    for (j=0; j<ncoef; j++) NIrow->coef[j] = coef[i][j];

    /* Write */
    rowNo = -1;
    ObitTableNIWriteRow (outTab, rowNo, NIrow, err);
    if (err->error) goto cleanup;


    /* Tell about it */
    if (prtLv>=1) {
      timer[0] = Time[i]-0.5*TimeI[i]; timer[1] = Time[i]+0.5*TimeI[i]; 
      TR2String (timer, msgtxt);
      Obit_log_error(err, OBIT_InfoErr, 
		     "%5d Time %s  RMS bias corr resid = %7.1f asec",i+1, msgtxt,rms[i]*3600.0);
      /* Give coefficients if fitted */
      if (NIrow->weight<=0.0) {
	Obit_log_error(err, OBIT_InfoWarn, "This interval to be ignored");
      } else {  /* OK */
	if (maxcoef==2)
	  Obit_log_error(err, OBIT_InfoErr, " %10.2f %10.2f ",
			 3600.0*coef[i][0], 3600.0*coef[i][1]);
	else
	  Obit_log_error(err, OBIT_InfoErr, " %10.2f %10.2f %10.2f %10.2f %10.2f",
			 3600.0*coef[i][0], 3600.0*coef[i][1], 3600.0*coef[i][2],
			 3600.0*coef[i][3], 3600.0*coef[i][4]);
	if (maxcoef>5) {
	  Obit_log_error(err, OBIT_InfoErr, " %10.2f %10.2f %10.2f %10.2f %10.2f",
			 3600.0*coef[i][5], 3600.0*coef[i][6], 3600.0*coef[i][7],
			 3600.0*coef[i][8], 3600.0*coef[i][9]);
	}
	if (maxcoef>10) {
	  Obit_log_error(err, OBIT_InfoErr, " %10.2f %10.2f %10.2f %10.2f %10.2f",
			 3600.0*coef[i][10], 3600.0*coef[i][11], 3600.0*coef[i][12],
			 3600.0*coef[i][13], 3600.0*coef[i][14]);
	}
	if (maxcoef>15) {
	  Obit_log_error(err, OBIT_InfoErr, " %10.2f %10.2f %10.2f %10.2f %10.2f",
			 3600.0*coef[i][15], 3600.0*coef[i][16], 3600.0*coef[i][17],
			 3600.0*coef[i][18], 3600.0*coef[i][19]);
	}
      }
    } /* end diagnostics */
  } /* end loop  L200: */

  /* Close NI table */
  ObitTableNIClose (outTab, err);
  if (err->error) goto cleanup;

  /* Total RMS */
  (*totRMS) = 3600.0 * rms[nTime];

  /* Tell about it */
  if (prtLv>=1) {
    if (ntoss>0)
      Obit_log_error(err, OBIT_InfoErr, "%s: Flagged %d of %d intervals for too few Zernike coef.",
		     routine, ntoss, nTime);
    if (nbad>0)
      Obit_log_error(err, OBIT_InfoErr, "%s: Flagged total %d of %d intervals",
		     routine, nbad, nTime);
    Obit_log_error(err, OBIT_InfoErr, "Total corrected rms residual = %7.1f asec",
		   *totRMS);
  
    /* Source offsets */
    for (i=0; i<nsou; i++) { /* loop 240 */
      if (gotsou[i]  &&  fitsou[i]) {
	Obit_log_error(err, OBIT_InfoErr, "Source %5d offset = %7.1f%7.1f ",
		       i+1, 3600.0*soffx[i], 3600.0*soffy[i]);
      } 
    } /* end loop  L240: */
  } /* end diagnostics */

 goto cleanup;

    /* Solution failed */
 AllBad:
  if (allFlag) Obit_log_error(err, OBIT_Error, 
			      "%s: All ionospheric solutions bad", routine);

  /* cleanup - deallocate arrays */
 cleanup:
  NIrow = ObitTableNIRowUnref (NIrow);
  /* deallocate arrays */
  if (mcoef)  g_free(mcoef);
  if (rms)    g_free(rms);
  if (soffx)  g_free(soffx);
  if (soffy)  g_free(soffy);
  if (gotsou) g_free(gotsou);
  if (fitsou) g_free(fitsou);
  if (coef)   {
    for (i=0; i<nTime; i++) g_free(coef[i]);
    g_free(coef);
  }
  if (err->error) Obit_traceback_val (err, routine, inUV->name, outTab);

  return outTab;
} /* end of routine doIonFitAll */ 

/**
 * Initialize Zernike model to apparent position offsets for a time  
 * series allowing for a systematic offset on the apparent source  
 * positions.   
 * Each time is analysed with fixed source positions and the soffx and  
 * soffy are set based on the average residual.  A comparison of the  
 * average residual with the total rms residual is used to decide if an   
 * offset is to be fitted to each source (fitsou). 
 * Offsets are only allowed for flqual values > 0 
 * Routine translated from the AIPSish IONCAL.FOR/INCINI  
 * \param nobs     Number of observations 
 * \param nsou     Number of sources 
 * \param nTime    Number of times   
 * \param maxcoe   Leading dimension of COEF 
 * \param isou     Source number per obs.  (0-rel)
 * \param iTime    Time interval number per observation.
 * \param x        Offset in field X (RA) on unit Zernike circle /source
 * \param y        Offset in field Y (Dec) on unit Zernike circle /source
 * \param dx       Apparent RA position shifts on the Zernike plane (deg) /obs
 * \param dy       Apparent dec position shifts on the Zernike plane (deg) /obs
 * \param w        Weight of observations  /obs
 * \param flqual   Field quality (crowding) code. Used in determining 
 *                 if the apparent position of a source can be moved. 
 * \param ncoef    Number of coefficents fitted per time 
 * \param gotsou   If true, there is data for this souce 
 * \param fitsou   If true fit offset correction for each source 
 * \param soffx    Initial guess for source offset in x=ra (deg) 
 * \param soffy    Initial guess for source offset in y=dec (deg) 
 * \param coef     Initial guess for Zernike coefficients (term,time) 
 * \param prtLv    Print level >=1 => give fitting diagnostics
 * \param err      Error/message stack 
 * \return TRUE if data OK, else too little data
 */
static gboolean 
initIonModel (gint nobs, olong nsou, olong nTime, olong maxcoe, olong* isou, 
	      olong* iTime, ofloat*  x, ofloat* y, ofloat* dx, ofloat* dy, 
	      ofloat* w, olong*  flqual, olong* ncoef, gboolean* gotsou, 
	      gboolean* fitsou, ofloat* soffx, ofloat* soffy, ofloat** coef, 
	      olong prtLv, ObitErr* err) 
{
  gboolean out=FALSE;
  olong   i, j, k, it, is, itim, rmscnt, numobs, numprm, itlast, itb, ite, js, ntgood;
  olong   qual;
  ofloat rms, dr, dd, rx, ry, rr;
  ofloat *sum2p1=NULL, *sum2p2=NULL, *sum3p1=NULL, *sum3p2=NULL;
  ofloat *tcoef=NULL, *gdx=NULL, *gdy=NULL;
  gchar *routine = "initIonModel";

    /* Error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return out;  /* previous error? */

  /* Initial guess (0's) */
  for (itim=0; itim<nTime; itim++) { /* loop 30 */
    ncoef[itim] = 0;
    for (i=0; i<maxcoe; i++) { /* loop 20 */
      coef[itim][i] = 0.0;
    } /* end loop  L20:  */
  } /* end loop  L30:  */

  /* if only one source just use tip and tilt */
  if (nsou==1) {
    for (itim=0; itim<nTime; itim++) { /* loop 30 */
      if (w[itim]>0.0) {
	gotsou[0]     = TRUE;
	ncoef[itim]   = 2;
	coef[itim][0] = dx[itim];
	coef[itim][1] = dy[itim];
	out = TRUE;
      }
    }
    return out;
  } /* end 1 source */

  /* Allocate work arrays */
  tcoef  = g_malloc(maxcoe*sizeof(ofloat));
  gdx    = g_malloc(nsou*maxcoe*sizeof(ofloat));
  gdy    = g_malloc(nsou*maxcoe*sizeof(ofloat));
  sum2p1 = g_malloc(nsou*sizeof(ofloat));
  sum2p2 = g_malloc(nsou*sizeof(ofloat));
  sum3p1 = g_malloc(nsou*sizeof(ofloat));
  sum3p2 = g_malloc(nsou*sizeof(ofloat));

  for (is=0; is<nsou; is++) { /* loop 40 */
    soffx[is] = 0.0;
    soffy[is] = 0.0;
    gotsou[is] = FALSE ;
    fitsou[is] = FALSE ;
    sum2p1[is] = 0.0;
    sum2p2[is] = 1.0e-20;
    sum3p1[is] = 0.0;
    sum3p2[is] = 1.0e-20;
  } /* end loop  L40:  */

  /* How many coeffients can actually  be fitted? */
  for (i=0; i<nobs; i++) { /* loop 60 */
    it = iTime[i];
    is = isou[i];
    if (w[i] > 0.0) {
      ncoef[it]++;
      gotsou[is] = TRUE;
    } 
  } /* end loop  L60:  */

  for (itim=0; itim<nTime; itim++) { /* loop 70 */
    ncoef[itim] = MIN (maxcoe, ncoef[itim]*2);
    if (ncoef[itim] < 5) ncoef[itim] =  MIN (2, ncoef[itim]);
  } /* end loop  L70:  */;
  
  /* Compute Zernike gradients */
  for (i=0; i<nsou; i++) {
    for (j=0; j<maxcoe; j++) {
      gdx[j*nsou+i] = ObitZernikeGradX(j+2, x[i], y[i]);
      gdy[j*nsou+i] = ObitZernikeGradY(j+2, x[i], y[i]);
    } /* end loop over coefficients */
  } /* end loop over sources */

  rms = 0.0;
  rmscnt = 0;
  numobs = 0;
  numprm = 0;
  itlast = 0;
  itb = 0;
  ntgood = 0;
  /* Loop over data fitting each time */
  for (i=0; i<nobs; i++) { /* loop 200 */
    it = iTime[i];
    if ((it > 0)  &&  (ncoef[it] > 0)) {
      is = isou[i];
      ite = i-1;
      /* New time? */
      if (it != itlast) {
	IonFit1 (i-itb, &isou[itb], x, y, &dx[itb], &dy[itb],& w[itb], 
		 &ncoef[itlast], &coef[itlast][0], prtLv, err);
	numprm = numprm + ncoef[itlast];
	if (ncoef[itlast] > 0) ntgood++;
	/* If fewer than 4 sources don't bother */
	if (nsou >= 4) {
	  
	  /* Debug Diagnostics */
	  if ((prtLv>=3) && (ncoef[itlast]<5)) {
	    Obit_log_error(err, OBIT_InfoErr, 
			   "%s: Bad time=%4d ncoef=%4d numobs=%5d numprm=%5d W= %10.3f %10.3f %10.3f %10.3f %10.3f ",
			   routine, itlast, ncoef[itlast], numobs, numprm, 
			   w[itb], w[itb+1], w[itb+2], w[itb+3], w[itb+4]);
	  }
	  
	  /* Get residuals */
	  for (j= itb; j<=ite; j++) { /* loop 150 */
	    if (w[j] > 0.0) {
	      js = isou[j];
	      /* current model */
	      dr = 0.0;
	      dd = 0.0;
	      for (k=0; k<ncoef[itlast]; k++) { /* loop 140 */
		/* model offset */
		dr += coef[itlast][k]*gdx[k*nsou+js];
		dd += coef[itlast][k]*gdy[k*nsou+js];
	      } /* end loop  L140: */;
	      /* Calculate residuals */
	      rx = dx[j] - dr;
	      ry = dy[j] - dd;
	      /* Residual statistics */

	      /* Debug Diagnostics */
	      if (prtLv>=3) {
		Obit_log_error(err, OBIT_InfoErr, 
			       "%s: S=%4d E= %5d resid %10.3f %10.3f ncoef=%4d",
			       routine, js+1, itlast+1, rx*3600.0, ry*3600.0, ncoef[itlast]);
		if (fabs(rx)>300.0) 
		  Obit_log_error(err, OBIT_InfoErr, 
				 "  coef: %10.3f %10.3f %10.3f %10.3f %10.3f ",
				 coef[itlast][0], coef[itlast][1], coef[itlast][2], 
				 coef[itlast][3], coef[itlast][4]);
	      }
	      rms += rx*rx + ry*ry;
	      rmscnt++;
	      /* Accululate sums */
	      sum2p1[js] += rx * w[j];
	      sum2p2[js] += w[j];
	      sum3p1[js] += ry * w[j];
	      sum3p2[js] += w[j];
	      numobs += 2;
	    } /* end if positive weight */ 
	  } /* end loop  L150: */
	} /* end if too few source */
	
	/* For next time interval */
	itb = i;
	itlast = it;
      } 
    }
  } /* end loop  L200: */
  

  /* Last time */
  itlast = nTime-1;
  ite = nobs-1;
  IonFit1 (i-itb,  &isou[itb], x, y, &dx[itb], &dy[itb], &w[itb], 
	   &ncoef[itlast], &coef[itlast][0], prtLv, err);
  numprm += ncoef[itlast];
  if (ncoef[itlast] > 0) ntgood++;
  
  /* If fewer than 4 sources don't  bother */
  if (nsou >= 4) {
    /* Get residuals */
    for (j= itb; j<=ite; j++) { /* loop 250 */
      if (w[j] > 0.0) {
	js = isou[j];
	/* current model */
	dr = 0.0;
	dd = 0.0;
	for (k=0; k<ncoef[itlast]; k++) { /* loop 240 */
	  /* model offset */
	  dr = dr + coef[itlast][k]*gdx[k*nsou+js];
	  dd = dd + coef[itlast][k]*gdy[k*nsou+js];
	} /* end loop  L240: */;
	/* Calculate residuals */
	rx = dx[j] - dr;
	ry = dy[j] - dd;
	rms += rx*rx + ry*ry;
	rmscnt++;
	/* Accululate sums */
	sum2p1[js] += rx * w[j];
	sum2p2[js] += w[j];
	sum3p1[js] += ry * w[j];
	sum3p2[js] += w[j];
	numobs += 2;
      } 
    } /* end loop  L250: */;
  } /* end if enough data */
  
  /* If too few sources skip */
  if (nsou < 4) goto cleanup;
  /* Bail out if no data */
  if (ntgood < 2) 
    Obit_log_error(err, OBIT_Error, 
		   "%s: Insufficient data for ionospheric model", routine);
  
  /* RMS */
  if (numobs > numprm)  rms = sqrt ((rms/rmscnt) * (numobs / (numobs - numprm)));
  else rms = -1.0;
 
  /* Tell about it */
  if (prtLv>=1) {
    Obit_log_error(err, OBIT_InfoErr, 
	        "%s: initial rms residual is %17.5f asec", routine, rms*3600.0);
  }

  /* Decide which sources to fit and  initial values */
  for (is=0; is<nsou; is++) { /* loop 300 */
    if (gotsou[is]) {
      rx = sum2p1[is] / sum2p2[is];
      ry = sum3p1[is] / sum3p2[is];
      rr = sqrt (rx*rx + ry*ry);

      /* Diagnostics */
      if (prtLv>=3) {
	Obit_log_error(err, OBIT_InfoErr, 
		       "%s: source %4d x_off %10.4f y_off %10.4f bad %d", 
		       routine, is+1, rx*3600.0, ry*3600.0, (rr>(0.5*rms)));
      }

      /* Find qual */
      qual = -1;
      for (k=0; k<nobs; k++) {
	if (is == isou[k]) {
	  qual = flqual[k];
	  break;
	}
      }

      /* Fit if residual > 0.5*RMS  and quality worse than 0 */
      if ((rr > (0.5*rms))  &&  (qual > 0)) {
	/*debug            IF ((RR.GT.(0.5*RMS))) THEN */
	soffx[is] = rx;
	soffy[is] = ry;
	fitsou[is] = TRUE;
      } 
    } 
  } /* end loop  L300: */
  out = TRUE; /* OK */
  
  /* Deallocate work arrays */
 cleanup:
  if (tcoef)  g_free(tcoef);
  if (gdx)    g_free(gdx);
  if (gdy)    g_free(gdy);
  if (sum2p1) g_free(sum2p1);
  if (sum2p2) g_free(sum2p2);
  if (sum3p1) g_free(sum3p1);
  if (sum3p2) g_free(sum3p2);
  return out;
} /* end of routine initIonModel */ 

/**
 * Final Ionospheric model for each time, uses fitted source offsets
 * Routine translated from the AIPSish IONCAL.FOR/INCINI  
 * \param nobs     Number of observations 
 * \param nsou     Number of sources 
 * \param nTime    Number of times   
 * \param maxcoe   Leading dimension of COEF 
 * \param isou     Source number per obs.  (0-rel)
 * \param iTime    Time interval number per observation.
 * \param x        Offset in field X (RA) on unit Zernike circle /source
 * \param y        Offset in field Y (Dec) on unit Zernike circle /source
 * \param dx       Apparent RA position shifts on the Zernike plane (deg) /obs
 * \param dy       Apparent dec position shifts on the Zernike plane (deg) /obs
 * \param w        Weight of observations  /obs
 * \param ncoef    Number of coefficents fitted per time 
 * \param fitsou   If true fit offset correction for each source 
 * \param soffx    Source offset in x=ra (deg) 
 * \param soffy    Source offset in y=dec (deg) 
 * \param coef     Initial guess for Zernike coefficients (term,time) , 
 *                 updated on return
 * \param prtLv    Print level >=1 => give fitting diagnostics
 * \param err      Error/message stack 
 * \return TRUE if data OK, else too little data
 */
static void
finalIonModel (gint nobs, olong nsou, olong nTime, olong maxcoe, olong* isou, 
	       olong* iTime, ofloat*  x, ofloat* y, ofloat* dx, ofloat* dy, 
	       ofloat* w, olong* ncoef, gboolean* fitsou, 
	       ofloat* soffx, ofloat* soffy, ofloat** coef, 
	       olong prtLv, ObitErr* err) 
{
  olong   i, j, k, it, is, itim, rmscnt, numobs, numprm, itlast, itb, ite, js, ntgood;
  ofloat rms, dr, dd, rx, ry;
  ofloat *tcoef=NULL, *gdx=NULL, *gdy=NULL;
  gchar *routine = "finalIonModel";

    /* Error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;  /* previous error? */

  /* if only one source just use tip and tilt */
  if (nsou==1) {
    for (itim=0; itim<nTime; itim++) { /* loop 30 */
      if (w[itim]>0.0) {
	ncoef[itim]   = 2;
	coef[itim][0] = dx[itim];
	coef[itim][1] = dy[itim];
      }
    }
  } /* end 1 source */

  /* Allocate work arrays */
  tcoef  = g_malloc(maxcoe*sizeof(ofloat));
  gdx    = g_malloc(nsou*maxcoe*sizeof(ofloat));
  gdy    = g_malloc(nsou*maxcoe*sizeof(ofloat));

  /* Compute Zernike gradients */
  for (i=0; i<nsou; i++) {
    for (j=0; j<maxcoe; j++) {
      gdx[j*nsou+i] = ObitZernikeGradX(j+2, x[i], y[i]);
      gdy[j*nsou+i] = ObitZernikeGradY(j+2, x[i], y[i]);
    } /* end loop over coefficients */
  } /* end loop over sources */

  rms = 0.0;
  rmscnt = 0;
  numobs = 0;
  numprm = 0;
  itlast = 0;
  itb = 0;
  ntgood = 0;
  /* Loop over data fitting each time */
  for (i=0; i<nobs; i++) { /* loop 200 */
    it = iTime[i];
    if ((it > 0)  &&  (ncoef[it] > 0)) {
      is = isou[i];
      ite = i-1;
      /* New time? */
      if (it != itlast) {
	IonFit1 (i-itb, &isou[itb], x, y, &dx[itb], &dy[itb],& w[itb], 
		 &ncoef[itlast], &coef[itlast][0], prtLv, err);
	numprm += ncoef[itlast];
	if (ncoef[itlast] > 0) ntgood++;
	/* If fewer than 4 sources don't bother */
	if (nsou >= 4) {
	  
	  /* Debug Diagnostics */
	  if ((prtLv>=3) && (ncoef[itlast]<5)) {
	    Obit_log_error(err, OBIT_InfoErr, 
			   "%s: Bad time=%4d ncoef=%4d numobs=%5d numprm=%5d W= %10.3f %10.3f %10.3f %10.3f %10.3f ",
			   routine, itlast, ncoef[itlast], numobs, numprm, 
			   w[itb], w[itb+1], w[itb+2], w[itb+3], w[itb+4]);
	  }
	  
	  /* Get residuals */
	  for (j= itb; j<=ite; j++) { /* loop 150 */
	    if (w[j] > 0.0) {
	      js = isou[j];
	      /* current model */
	      dr = 0.0;
	      dd = 0.0;
	      for (k=0; k<ncoef[itlast]; k++) { /* loop 140 */
		/* model offset */
		dr += coef[itlast][k]*gdx[k*nsou+js];
		dd += coef[itlast][k]*gdy[k*nsou+js];
	      } /* end loop  L140: */;
	      /* Calculate residuals */
	      rx = dx[j] - dr;
	      ry = dy[j] - dd;
	      /* Residual statistics */

	      /* Debug Diagnostics */
	      if (prtLv>=3) {
		Obit_log_error(err, OBIT_InfoErr, 
			       "%s: S=%4d E= %5d resid %10.3f %10.3f ncoef=%4d",
			       routine, js+1, itlast+1, rx*3600.0, ry*3600.0, ncoef[itlast]);
		if (fabs(rx)>300.0) 
		  Obit_log_error(err, OBIT_InfoErr, 
				 "  coef: %10.3f %10.3f %10.3f %10.3f %10.3f ",
				 coef[itlast][0], coef[itlast][1], coef[itlast][2], 
				 coef[itlast][3], coef[itlast][4]);
	      }
	      rms += rx*rx + ry*ry;
	      rmscnt++;
	      numobs += 2;
	    } /* end if positive weight */ 
	  } /* end loop  L150: */
	} /* end if too few source */
	
	/* For next time interval */
	itb = i;
	itlast = it;
      } 
    }
  } /* end loop  L200: */
  

  /* Last time */
  itlast = nTime-1;
  ite = nobs-1;
  IonFit1 (i-itb,  &isou[itb], x, y, &dx[itb], &dy[itb], &w[itb], 
	   &ncoef[itlast], &coef[itlast][0], prtLv, err);
  numprm += ncoef[itlast];
  if (ncoef[itlast] > 0) ntgood++;
  
  /* If fewer than 4 sources don't  bother */
  if (nsou >= 4) {
    /* Get residuals */
    for (j= itb; j<=ite; j++) { /* loop 250 */
      if (w[j] > 0.0) {
	js = isou[j];
	/* current model */
	dr = 0.0;
	dd = 0.0;
	for (k=0; k<ncoef[itlast]; k++) { /* loop 240 */
	  /* model offset */
	  dr = dr + coef[itlast][k]*gdx[k*nsou+js];
	  dd = dd + coef[itlast][k]*gdy[k*nsou+js];
	} /* end loop  L240: */;
	/* Calculate residuals */
	rx = dx[j] - dr;
	ry = dy[j] - dd;
	rms += rx*rx + ry*ry;
	rmscnt++;
	numobs += 2;
      } 
    } /* end loop  L250: */;
  } /* end if enough data */
  
  /* If too few sources skip */
  if (nsou < 4) goto cleanup;
  /* Bail out if no data */
  if (ntgood < 2) 
    Obit_log_error(err, OBIT_Error, 
		   "%s: Insufficient data for ionospheric model", routine);
  
  /* RMS */
  if (numobs > numprm)  rms = sqrt ((rms/rmscnt) * (numobs / (numobs - numprm)));
  else rms = -1.0;
 
  /* Tell about it */
  if (prtLv>=1) {
    Obit_log_error(err, OBIT_InfoErr, 
	        "%s: initial rms residual is %17.5f asec", routine, rms*3600.0);
  }

  /* Deallocate work arrays */
 cleanup:
  if (tcoef)  g_free(tcoef);
  if (gdx)    g_free(gdx);
  if (gdy)    g_free(gdy);
  return;
} /* end of routine finalIonModel */ 



/**
 * Fit Zernike model to apparent position offsets for a time series  
 * allowing for a systematic offset on the apparent source positions.  
 * Solution uses a relaxation method from Fred Schwab:  
 * Pn+1 = Pn + atan2 (dChi2/dP), (d2Chi2/dP2))  
 * for each parameter P where n or n+1 indicates a given iteration,   
 * dChi2/dP is the first partial derivative of Chi squared wrt P,  
 * d2Chi2/d2P is the second partial derivative of Chi squared wrt P,  
 * Chi2 = Sum (w Abs(dxj - oxk - Sum (Pi Gxik))**2) +   
 * Sum (w Abs(dyj - oyk - Sum (Pi Gyik))**2)  
 * Must be initialized by a call to initIonModel.
 * Routine translated from the AIPSish IONCAL.FOR/IOALFT  
 *  
 * dxj = apparent offset in x of obs j  
 * dyj = apparent offset in y of obs j  
 * oxk = correction to source k offset in x, k = f(j)  
 * oyk = correction to source k offset in y, k = f(j)  
 * Pi  = Zernike coefficient i [COEF(I) below],  
 * Gxik = Zernike x gradient term i for source k        
 * Gyik = Zernike y gradient term i for source k     
 * \param nobs      Number of observations 
 * \param nsou      Number of sources  
 * \param nTime     Number of times    
 * \param maxcoe    Max. number of coefficients
 * \param isou      Source number per obs.  (0-rel)
 * \param iTime     Time interval number per observation
 * \param x         Offset in field X (RA) on unit Zernike circle /source
 * \param y         Offset in field Y (Dec) on unit Zernike circle /source
 * \param dx        Apparent RA position shifts on the Zernike plane (deg) 
 * \param dy        Apparent dec position shifts on the Zernike plane 
 * \param w         Weights of observations 
 * \param ncoef     Number of coefficents fitted per time 
 * \param gotsou    If true, there is data for this souce 
 * \param fitsou    If true fit offset correctsion for each source 
 * \param soffx     Correction to source offset in x=ra (deg) 
 *                  input: initial guess, output: fitted values 
 * \param soffy     Correction to source offset in y=dec (deg) 
 *                  input: initial guess, output: fitted values 
 * \param coef      Zernike coefficients fitted (term,time) 
 *                  input: initial guess, output: fitted values 
 * \param trms      RMS residual for each time interval (deg) 
 *                  value [nTime+] is the total corrected RMS 
 *                  Values corrected for the number of parameters.
 *                  -1 => not determined. 
 * \param prtLv     Print level >=1 => give fitting diagnostics
 * \param err       Error stack
 */
static void 
FitIonSeries (gint nobs, olong nsou, olong nTime, olong maxcoe, olong* isou, olong* iTime, 
	      ofloat*  x, ofloat* y, ofloat* dx, ofloat* dy, ofloat* w, olong*  ncoef, 
	      gboolean* gotsou, gboolean* fitsou, ofloat* soffx, ofloat* soffy, 
	      ofloat** coef, ofloat* trms, olong prtLv, ObitErr* err) 
{
  olong   i, j, it, is, itim, rmscnt, numobs, numprm, iter;
  gboolean   convgd, OK;
  ofloat      rms=0.0, dr, dd, delta, tol, norm=0.0, test, rx, ry, pd1, pd2, wx, rmslst ;
  ofloat **sum1p1=NULL, **sum1p2=NULL, *sum2p1=NULL, *sum2p2=NULL, *sum3p1=NULL, *sum3p2=NULL;
  ofloat **tcoef=NULL, *tsoffx=NULL, *tsoffy=NULL, *gdx=NULL, *gdy=NULL, *sum1=NULL, *sum2=NULL;
  /* Penalty terms by order of Zernike 1-4 */
  ofloat pen[]={0.001, 0.001,
		0.01, 0.01, 0.01,
		0.03, 0.03, 0.03, 0.03, 0.03, 
		0.05, 0.05, 0.05,0.05,  0.05, 0.05, 0.05, 
		0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08};
  ofloat *sumWt=NULL;
  gchar *routine = "FitIonSeries";

  /* Error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return ;  /* previous error? */

  /* Allocate work arrays */
  tsoffx = g_malloc(nsou*sizeof(ofloat));
  tsoffy = g_malloc(nsou*sizeof(ofloat));
  gdx    = g_malloc(nsou*maxcoe*sizeof(ofloat));
  gdy    = g_malloc(nsou*maxcoe*sizeof(ofloat));
  sum2p1 = g_malloc(nsou*sizeof(ofloat));
  sum2p2 = g_malloc(nsou*sizeof(ofloat));
  sum3p1 = g_malloc(nsou*sizeof(ofloat));
  sum3p2 = g_malloc(nsou*sizeof(ofloat));
  sum1   = g_malloc(nTime*sizeof(ofloat));
  sum2   = g_malloc(nTime*sizeof(ofloat));
  sum1p1 = g_malloc(nTime*sizeof(ofloat*));
  sum1p2 = g_malloc(nTime*sizeof(ofloat*));
  tcoef  = g_malloc(nTime*sizeof(ofloat*));
  sumWt  = g_malloc(nTime*sizeof(ofloat));
  for (i=0; i<nTime; i++) {
    sum1p1[i] = g_malloc(maxcoe*sizeof(ofloat));
    sum1p2[i] = g_malloc(maxcoe*sizeof(ofloat));
    tcoef[i]  = g_malloc(maxcoe*sizeof(ofloat));
  }

  /* Compute Zernike gradients */
  for (i=0; i<nsou; i++) {
    for (j=0; j<maxcoe; j++) {
      gdx[j*nsou+i] = ObitZernikeGradX(j+2, x[i], y[i]);
      gdy[j*nsou+i] = ObitZernikeGradY(j+2, x[i], y[i]);
    } /* end loop over coefficients */
  } /* end loop over sources */

  rmslst = 1.0e20;
  for (i=0; i<nTime; i++) trms[i] = 1.0;

  /* Loop over iterations */
  for (iter=1; iter<=200; iter++) { /* loop 600 */

    convgd = TRUE;  /* Can decide otherwise: */

    /* Zero sums */
    for (j=0; j<maxcoe; j++) { /* loop 130 */
      for (itim=0; itim<nTime; itim++) { /* loop 120 */
	sum1p1[itim][j] = 0.0;
	sum1p2[itim][j] = 1.0e-20;
      } /* end loop  L120: */;
    } /* end loop  L130: */
    for (j=0; j<nsou; j++) {
      sum2p1[j] = 0.0;
      sum2p2[j] = 1.0e-20;
      sum3p1[j] = 0.0;
      sum3p2[j] = 1.0e-20;
    } 
    
    /* Loop over data doing sums */
    for (it=0; it<nTime; it++) sumWt[it] = 0.0;
    for (i=0; i<nobs; i++) { /* loop 200 */
      if (w[i] > 0.0) {
	it = iTime[i];
	is = isou[i];
	sumWt[it] += w[i];  /* Sum weights per time*/
	/* current model */
	dr = 0.0;
	dd = 0.0;
	for (j=0; j<ncoef[it]; j++) { /* loop 140 */
	  /* model offset */
	  dr += coef[it][j]*gdx[j*nsou+is];
	  dd += coef[it][j]*gdy[j*nsou+is];
	} /* end loop  L140: */;
	
	/* Calculate residuals */
	if (fitsou[is]) {
	  rx = dx[i] - soffx[is] - dr;
	  ry = dy[i] - soffy[is] - dd;
	} else {
	  rx = dx[i] - dr;
	  ry = dy[i] - dd;
	}

	/* Partial derivatives for coef */
	for (j=0; j<ncoef[it]; j++) { /* loop 150 */
	  pd1 = -rx * gdx[j*nsou+is] - ry * gdy[j*nsou+is];
	  pd2 = gdx[j*nsou+is]*gdx[j*nsou+is] + gdy[j*nsou+is]*gdy[j*nsou+is];

	  /* Sum */
	  sum1p1[it][j] += w[i] * pd1;
	  sum1p2[it][j] += w[i] * pd2;
	} /* end loop  L150: */;
	
	/* Fitting this source? */
	if (fitsou[is]) {
	  /* Partial derivatives for soffx */
	  pd1 = -rx;
	  pd2 = 1.0;
	  /* Sum */
	  sum2p1[is] += w[i] * pd1;
	  sum2p2[is] += w[i] * pd2;
	  /* Partial derivatives for soffy */
	  pd1 = -ry;
	  pd2 = 1.0;
	  /* Sum */
	  sum3p1[is] += w[i] * pd1;
	  sum3p2[is] += w[i] * pd2;
	} 
      } 
    } /* end loop  L200: */
    
    /* Penalty functions for non zero higher order Zernike coefs */
    for (it=0; it<nTime; it++) { /* loop 230 */
      for (j=0; j<ncoef[it]; j++) { 
	sum1p1[it][j] += coef[it][j] * sumWt[it]*pen[j];
	sum1p2[it][j] += sumWt[it]*pen[j];
      }
    }

    /* Update solutions */
    wx = 1.6;
    OK = FALSE;
    while (!OK) {
      wx = wx * 0.5;
      /* don't loop forever */
      if (wx < 1.0e-10) goto converged;

      /* Convergence criterion - lower the bar  */
      tol = 5.0e-6 + iter * 1.0e-5;
      norm = 0.0;
      numobs = 0;
      numprm = 0;

      /* Coefficients */
      for (itim=0; itim<nTime; itim++) { /* loop 230 */
	if (sumWt[itim] > 0.0) {
	  for (j=0; j<ncoef[itim]; j++) { /* loop 220 */
	    numprm++;
	    delta = atan2 (sum1p1[itim][j], sum1p2[itim][j]);
	    test = tol;
	    tcoef[itim][j] = coef[itim][j] - wx * delta;

	    /* DEBUG
	       if (fabs(tcoef[itim][j])>1.0) {
	       fprintf (stdout, "Lookout time %d parm %d coef %f part %f %f\n",
	       itim, j,tcoef[itim][j], sum1p1[itim][j], sum1p2[itim][j]);
	       }  End DEBUG */

	    /* Convergence test */
	    convgd = convgd && (fabs(delta) <= test);
	    norm += delta*delta;
	  } /* end loop  L220: */;
	} /* end valid data for interval */
	else {  /* don't change coefs */
	  for (j=0; j<ncoef[itim]; j++) tcoef[itim][j] = coef[itim][j];
	}
      } /* end loop  L230: */;

      /* Position offsets */
      for (is=0; is<nsou; is++) { /* loop 250 */
	if (gotsou[is]  &&  fitsou[is]) {

	  /* X offset */
	  delta = atan2 (sum2p1[is], sum2p2[is]);
	  test = tol ;
	  tsoffx[is] = soffx[is] - wx * delta;

	  /* Convergence test */
	  convgd = convgd && (fabs(delta) <= test);
	  norm +=delta*delta;

	  /* Y offset */
	  delta = atan2 (sum3p1[is], sum3p2[is]);
	  test = tol ;
	  tsoffy[is] = soffy[is] - wx * delta;

	  /* Convergence test */
	  convgd = convgd  &&  (fabs(delta) <= test);
	  numprm += 2;
	  norm += delta*delta;

	  /* Debug Diagnostics */
	  if (prtLv>=4) {
	    Obit_log_error(err, OBIT_InfoErr, 
			   "%s: iter %d source %4d off_x=%8.2f off_y=%8.2f delta_y=%8.2f",
			   routine, iter+1, is+1, soffx[is]*3600.0, soffy[is]*3600.0, delta*3600.0);
	  }
	  
	} else {
	  tsoffy[is] = 0.0;
	  tsoffx[is] = 0.0;
	} 
      } /* end loop  L250: */

      /* Determine RMS */
      rms    = 0.0;
      rmscnt = 0;

      /* Time RMS sums */
      for (j=0; j<nTime; j++) { /* loop 310 */
	sum1[j] = 0.0;
	sum2[j] = 1.0e-20;
      } /* end loop  L310: */;

      for (i=0; i<nobs; i++) { /* loop 340 */
	if (w[i] > 0.0) {
	  it = iTime[i];
	  is = isou[i];
	  numobs += 2;
	  /* current model */
	  dr = 0.0;
	  dd = 0.0;
	  for (j=0; j<ncoef[it]; j++) { /* loop 330 */
	    /* model offset */
	    dr += tcoef[it][j]*gdx[j*nsou+is];
	    dd += tcoef[it][j]*gdy[j*nsou+is];
	  } /* end loop  L330: */;

	  /* Calculate residuals */
	  if (fitsou[is]) {
	    rx = dx[i] - tsoffx[is] - dr;
	    ry = dy[i] - tsoffy[is] - dd;
	  } else {
	    rx = dx[i] - dr;
	    ry = dy[i] - dd;
	  }

	  /* Residual statistics */
	  rms += rx*rx + ry*ry;       /* Total */
	  rmscnt++;
	  sum1[it] += rx*rx + ry*ry;  /* Per time */
	  sum2[it] += 1.0;         /* NOTE: this was wrong in AIPS */
	  /* DEBUG
	  if (it==4)
	    fprintf (stdout,"Src %5d residX %10.5f residY %10.5f model %10.5f%10.5f \n",
		     is+1, rx*3600.0, ry*3600.0, dr*3600.0, dd*3600.0); */
	} 
      } /* end loop  L340: */

      /* Force residuals to decrease */
      if (numobs > numprm) rms = sqrt ((rms/rmscnt) * (numobs / (numobs - numprm)));
      else rms = -1.0;
      OK =  (rms <= rmslst);
      rmslst = rms;
    } /* end of loop seeking improvement of fit */

    /* Save values */
    trms[nTime] = rms;  /* total RMS */

    for (j=0; j<nTime; j++) { /* loop 360 */
      /* Any data and enought deg. of freedom to fit? */
      if ((sumWt[j]>0.0) && ((2.0*sum2[j]) > (ncoef[j]+1))) {
	trms[j] = sqrt ((sum1[j] / sum2[j]) * 
			(2.0*sum2[j] / (2.0*sum2[j] - ncoef[j]))) ;
	  /* DEBUG 
	  if (j==4)
	    fprintf (stdout, "RMS %10.5f\n",trms[j]*3600.0);*/
      } else {
	trms[j] = -1.0;
      }

      /* DEBUG
	 if ((fabs(tcoef[j][0])>1.0) || (fabs(tcoef[j][1])>1.0)) {
	 fprintf (stdout, "Bother iter %d time %d coef %f %f\n",
	 it, j+1, tcoef[j][0], tcoef[j][1]);
	 }  End DEBUG */
      
      for (i=0; i<ncoef[j]; i++) { /* loop 350 */
	coef[j][i] = tcoef[j][i];
      } /* end loop  L350: */;
    } /* end loop  L360: */;
    
    for (i=0; i<nsou; i++) { /* loop 370 */
      soffx[i] = tsoffx[i];
      soffy[i] = tsoffy[i];
    } /* end loop  L370: */;

    /* Debug Diagnostics */
    if ((prtLv>=3) && (iter<=2)) {
      /*j = nTime-1;
	Obit_log_error(err, OBIT_InfoErr, 
	"iter %4d coef %10.5f %10.5f %10.5f %10.5f %10.5f ",
	iter, 3600.0*coef[j][0], 3600.0*coef[j][1], 3600.0*coef[j][2], 
	3600.0*coef[j][3], 3600.0*coef[j][4]);*/
      Obit_log_error(err, OBIT_InfoErr, 
		     "iter %4d rms=%10.5f norm=%12.5g", iter, rms*3600.0, norm);
    }
   
    /* Converged? */
    if (convgd) break;
  } /* end loop iteration L600: */;

    /* Converged */
  converged:

  /* Tell about it */
  if (prtLv>=2) {
    Obit_log_error(err, OBIT_InfoErr, 
		   "%s: No. iteration %5d RMS resid=%17.5f norm=%12.5g", 
		   routine, iter+1, rms*3600.0, norm);
  }

  /* Deallocate work arrays */
  if (tsoffx) g_free(tsoffx);
  if (tsoffy) g_free(tsoffy);
  if (sumWt)  g_free(sumWt);
  if (gdx)    g_free(gdx);
  if (gdy)    g_free(gdy);
  if (sum2p1) g_free(sum2p1);
  if (sum2p2) g_free(sum2p2);
  if (sum3p1) g_free(sum3p1);
  if (sum3p2) g_free(sum3p2);
  if (sum1)   g_free(sum1);
  if (sum2)   g_free(sum2);
  if (tcoef) {
    for (i=0; i<nTime; i++) {
      if (tcoef[i]) g_free(tcoef[i]);
    } 
    g_free(tcoef);
  }
  if (sum1p1){
    for (i=0; i<nTime; i++) {
      if (sum1p1[i]) g_free(sum1p1[i]);
    } 
    g_free(sum1p1);
  }
  if (sum1p2){
    for (i=0; i<nTime; i++) {
      if (sum1p2[i]) g_free(sum1p2[i]);
    } 
    g_free(sum1p2);
  }
} /* end of routine FitIonSeries */ 

/**
 * Edit ionospheric data based on model fit.  
 * Times with either insufficient data to determine a decent solution  
 * (4 sources) or excessive RMS residuals will be flagged by setting  
 * the corresponding ncoef to zero.  
 * Ant intervals with fewer coefficients in the model than the maximum
 * will be flagged.
 * Routine translated from the AIPSish IONCAL.FOR/IONEDT  
 * \param nobs      Number of observations 
 * \param nsou      Number of sources  
 * \param nTime     Number of times  
 * \param maxcoe    Maximum number of coefficients
 * \param isou      Source number per obs.  (0-rel)
 * \param iTime     Time interval number per observation.
 * \param x         Offset in field X (RA) on unit Zernike circle /source
 * \param y         Offset in field Y (Dec) on unit Zernike circle  /source
 * \param dx        Apparent RA position shifts on the Zernike plane (deg)
 * \param dy        Apparent dec position shifts on the Zernike plane  (deg)
 * \param w         Weight of observations 
 *                  On output, set to zero to remove datum 
 * \param MaxRMS    Maximum acceptable RMS residual (deg) 
 * \param ncoef     Number of coefficents fitted per time 
 *                  On output, set to -1 to remove time. 
 * \param soffx     Source offset in x=ra (deg) 
 * \param soffy     Source offset in y=dec (deg) 
 * \param coef     Zernike coefficients (term,time) 
 * \param prtLv    Print level >=1 => give fitting diagnostics
 * \param err      Error/message stack 
 * \return TRUE if data edited.
 */
static gboolean 
IonEditSeries (gint nobs, olong nsou, olong nTime, 
	gint maxcoe, olong* isou, olong* iTime, ofloat*  x, ofloat* y, 
	ofloat* dx, ofloat* dy, ofloat* w, ofloat MaxRMS, olong* ncoef, 
	ofloat* soffx, ofloat* soffy, ofloat** coef, olong prtLv, ObitErr* err) 
{
  gboolean out=FALSE;
  olong   i, j, it, itlast, itb, ite, count, ntoss, stoss, numt;
  gboolean  doMore;
  ofloat    rms;
  ofloat *lx=NULL, *ly=NULL, *gdx=NULL, *gdy=NULL;
  gchar *routine = "IonEditSeries";

  /* Error checks */
  if (err->error) return out;  /* previous error? */

  /* Allocate work arrays */
  gdx    = g_malloc(nsou*maxcoe*sizeof(ofloat));
  gdy    = g_malloc(nsou*maxcoe*sizeof(ofloat));
  lx     = g_malloc(nsou*sizeof(ofloat));
  ly     = g_malloc(nsou*sizeof(ofloat));

  /* Compute Zernike gradients */
  for (i=0; i<nsou; i++) {
    for (j=0; j<maxcoe; j++) {
      gdx[j*nsou+i] = ObitZernikeGradX(j+2, x[i], y[i]);
      gdy[j*nsou+i] = ObitZernikeGradY(j+2, x[i], y[i]);
    } /* end loop over coefficients */
  } /* end loop over sources */

  numt  = 0;
  ntoss = 0;
  stoss = 0;
  itlast = iTime[0];
  itb = 0;
  /* Loop over data examining each time */
  for (i=0; i<nobs; i++) { /* loop 200 */
    it = iTime[i];
    ite = i-1;
    
    /* New time? Last? */
    if ((it != itlast)  ||  (i == (nobs-1))) {
      if (i == (nobs-1)) ite = i;  /* Last? */
      /* Debug Diagnostics */
      if (prtLv>=3) {
	Obit_log_error(err, OBIT_InfoErr, 
		       "%s: time %4d coef %10.5f %10.5f %10.5f %10.5f %10.5f ",
		       routine,itlast+1, 3600.0*coef[itlast][0], 3600.0*coef[itlast][1], 
		       3600.0*coef[itlast][2], 3600.0*coef[itlast][3], 
		       3600.0*coef[itlast][4]);
      }
      
      /* Begin editing loop */
      /* Copy actual source Zernike offsets to lx, ly */
      count = 0;
      for (i=itb; i<=ite; i++) {
	lx[count] = x[isou[i]];
	ly[count] = y[isou[i]];
	count++;
      }
      
      /* Fit model/edit loop */
      doMore = TRUE;
      while (doMore) {
	rms = IonFit1 (count, &isou[itb], x, y, &dx[itb], &dy[itb], &w[itb], 
		       &ncoef[itlast], &coef[itlast][0], prtLv, err);
	if (rms<MaxRMS) break;  /* Reached goal? */
	doMore = IonEdit1 (count, &isou[itb], x, y, &dx[itb], &dy[itb], &w[itb],  
			   MaxRMS, ncoef[itlast], &coef[itlast][0], prtLv, err);
	if (doMore) stoss++;
      }

      
      /* If RMS still excessive toss  interval */
      if (rms > MaxRMS) {
	if (prtLv>=2) {
	Obit_log_error(err, OBIT_InfoErr, 
		       "%s: Flagged time %d, excessive rms %f > %f",
		       routine, itlast+1, rms*3600.0, MaxRMS*3600.0);
	}
	ncoef[itlast] = -1;
	ntoss++;
	/* Reject data */
	for (j= itb; j<=ite; j++) { /* loop 180 */
	  w[j] = 0.0;
	} /* end loop  L180: */
      } 
      
      /* For next time interval */
      itb = i;
      itlast = it;
      numt = MAX (numt, itlast);
    } /* end if new time */ 
  } /* end loop  L200: */

  /* Filter each source */
  for (i=0; i<nsou; i++) { /* loop 10 */
    /* Have adjusted the positions with no source offsets */
    soffx[i] = 0.0;
    soffy[i] = 0.0;
    IonEditSource (i, nobs, nTime, maxcoe, nsou, isou, iTime,  x, y, dx, dy, w, 
		   MaxRMS, ncoef, soffx, soffy, coef, gdx, gdy, &stoss, 
		   prtLv, err);
  } /* end loop  L10:  */

  /* Tell about it */
  if (prtLv>=1) {
    Obit_log_error(err, OBIT_InfoErr, 
		   "%s: rejected %5d of %5d intervals, %6d single obs.", 
		   routine, ntoss, numt+1, stoss);
  }
  
  /* Do anything? */
  out = ((ntoss > 0)  ||  (stoss > 0));

  /* Deallocate work arrays */
  if (lx)    g_free(lx);
  if (ly)    g_free(ly);
  if (gdx)   g_free(gdx);
  if (gdy)   g_free(gdy);
  return out;
} /* end of routine IonEditSeries */ 


/**
 * Edit times if average measured flux is less that MinRat wrt average  
 * Times will be flagged by setting the corresponding ncoef to zero.  
 * Routine translated from the AIPSish IONCAL.FOR/IONAED  
 * \param nobs      Number of observations 
 * \param nsou      Number of sources  
 * \param nTime     Number of times  
 * \param isou      Source number per obs.  (0-rel)
 * \param iTime     Time number per obs. 
 * \param MinRat    Minimum acceptable ratio to average flux 
 * \param  flux     Peak flux density of observations 
 * \param ncoef     Number of coefficents fitted per time 
 *                  On output, set to -1 to remove time. 
 * \param wt        Weights of observations, set to 0 if flagged 
 * \param prtLv     Print level >=1 => give fitting diagnostics
 * \param err       Error/message stack 
 * \return TRUE if data edited.
 */
static gboolean 
IonEditAmp (gint nobs, olong nsou, olong nTime, olong* isou, olong* iTime, 
	    ofloat MinRat, ofloat*  flux, olong* ncoef, ofloat* wt, 
	    olong prtLv, ObitErr* err) 
{
  gboolean out = FALSE;
  olong   *cntflx, i, j, isss, ittt, ilast, drop, total, ib, ie ;
  gboolean bad;
  ofloat   *sumflx, sumt1, sumt2;
  gchar flg[11];
  gchar *routine = "IonEditAmp";

  /* Error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return out;  /* previous error? */

  /* allocate arrays */
  cntflx  = g_malloc0(nsou*sizeof(olong));
  sumflx  = g_malloc0(nsou*sizeof(ofloat));

  if (prtLv>=1) {
    Obit_log_error(err, OBIT_InfoErr, 
		   "%s: filtering solutions with min. flux ratio %10.5f", 
		   routine,MinRat );
  }

  drop = 0;
  total = 0;
  /* Average source flux densities */
  for (i=0; i<nsou; i++) cntflx[i] = 0;
  for (i=0; i<nsou; i++) sumflx[i] = 0.0;
  for (i=0; i<nobs; i++) { /* loop 20 */
    if (flux[i] > 0) {
      isss = isou[i];
      cntflx[isss]++;
      sumflx[isss] += flux[i];
    } 
  } /* end loop  L20:  */

  /* Normalize */
  for (i=0; i<nsou; i++) { /* loop 40 */
    if (cntflx[i] > 0) sumflx[i] /= cntflx[i];
  } /* end loop  L40:  */;

  /* Loop over times summing measured and averaged fluxes */
  ilast = 0;
  ib = 0;
  ie = 0;
  sumt1 = 0.0;
  sumt2 = 0.0;
  for (i=0; i<nobs; i++) { /* loop 100 */
    if (iTime[i] > ilast) {
      total++;
      /* New time - check results */
      if (sumt2 > 0.0) sumt1 /= sumt2;
      bad = sumt1  <  MinRat;
      ittt = iTime[i];
      if (bad) {
	/* Flag time/data */
	ncoef[ittt] = -1;
	for (j= ib; j<= ie; j++) { /* loop 80 */
	  flux[j] = 0.0;
	  wt[j]   = 0.0;
	} /* end loop  L80:  */;
	drop++;
      } 
      /*   Debug  Diagnostics*/
      strcpy (flg, " ");
      if (bad) strcpy (flg, "flagged");
      if (prtLv>=3) {
	Obit_log_error(err, OBIT_InfoErr, 
		       "%s: Time %4d average flux= %10.2f %s", 
		       routine,ilast, sumt1, flg);
      }
      sumt1 = 0.0;
      sumt2 = 0.0;
      ilast = ittt;
      ib = i;
    } 
    ie = i;

    /* Sum Fluxes and average values for  source.  */
    if (flux[i] > 0.0) {
      isss = isou[i];
      sumt1 += flux[i];
      sumt2 += sumflx[isss];
    } 
  } /* end loop  L100: */

  /* Deal with last average */
  if (sumt2 > 0.0) sumt1 /= sumt2;
  bad = sumt1  <  MinRat;
  if (bad) {
    /* Flag time/data */
    ncoef[ilast] = -1;
    for (j= ib; j<= ie; j++) { /* loop 120 */
      wt[j] = 0.0;
      flux[j] = 0.0;
    } /* end loop  L120: */;
    drop++;
  }

  /*   Debug Diagnostics */
  strcpy (flg, " ");
  if (bad) strcpy (flg, "flagged");
  if (prtLv>=3) {
    Obit_log_error(err, OBIT_InfoErr, 
		   "%s: Time %4d average flux= %10.2f %s", 
		   routine, ilast, sumt1, flg);
  }

  /* Tell results */
  if (prtLv>=1) {
    Obit_log_error(err, OBIT_InfoErr, 
		   "%s: reject %4d of %4d times due to low peak", 
		   routine, drop, total);
  }

    if (cntflx) g_free(cntflx);
    if (sumflx) g_free(sumflx);
    out = (drop >0);  /* Do anything? */
    return out;
} /* end of routine IonEditAmp */ 

/**
 * Edit ionospheric data based on model fit for a given source.  
 * Routine translated from the AIPSish IONCAL.FOR/INESOU  
 * \param source    Source number to edit 
 * \param nobs      Number of observations 
 * \param nTime     Number of times 
 * \param maxcoe    Maximum number of coefficients
 * \param nsou      Number of sources
 * \param isou      Source number per obs.  (0-rel)
 * \param iTime     Time interval number per observation.
 * \param x         Offset in field X (RA) on unit Zernike circle /source
 * \param y         Offset in field Y (Dec) on unit Zernike circle  /source
 * \param dx        Apparent RA position shifts on the Zernike plane (deg) 
 * \param dy        Apparent dec position shifts on the Zernike plane (deg) 
 * \param w         Weight of observations 
 *                  On output, set to zero to remove datum 
 * \param MaxRMS    Maximum acceptable RMS residual 
 * \param ncoef     Number of coefficents fitted per time 
 *                  On output, set to -1 to remove time. 
 * \param soffx     source offset in x=ra (deg) 
 * \param soffy     source offset in y=dec (deg) 
 * \param coef      Zernike coefficients (term,time) 
 * \param gdx       Zernike gradients in X per source, coef
 * \param gdy       Zernike gradients in Y per source, coef
 * \param stoss     [in/out] Number of observations rejected 
 * \param prtLv     Print level >=1 => give fitting diagnostics
 * \param err       Error stack
 */
static void IonEditSource (gint source, olong nobs, olong nTime, olong maxcoe, olong nsou, 
			   olong* isou, olong* iTime, ofloat*  x, ofloat* y, 
			   ofloat* dx, ofloat* dy, ofloat* w, ofloat MaxRMS, olong* ncoef, 
			   ofloat* soffx, ofloat* soffy, ofloat** coef, 
			   ofloat *gdx, ofloat *gdy, olong* stoss, olong prtLv, ObitErr* err) 
{
  olong   i, j, k, it, rmscnt, js, numt, *ktime, count, i1, i2, wid, idt, nt, *iobs;
  ofloat  rms, dr, dd, sum;
  ofloat  *resid, rx, ry, test, rmst;
  gchar *routine = "IonEditSource";  

  /* allocate arrays */
  resid  = g_malloc0(nTime*sizeof(ofloat));
  ktime  = g_malloc0(nTime*sizeof(olong));
  iobs   = g_malloc0(nTime*sizeof(olong));

  numt   = 0;
  rms    = 0.0;
  rmscnt = 0;
  /* Loop over data examining each time */
  for (i=0; i<nobs; i++) { /* loop 200 */
    if ((isou[i] == source)  &&  (w[i] > 0.0)) {
      it = iTime[i];
      ktime[numt] = iTime[i];
      iobs[numt]  = i;
      js          = isou[i];
      /* current model */
      dr = 0.0;
      dd = 0.0;
      for (k=0; k<ncoef[it]; k++) { /* loop 140 */
	/* model offset */
	dr += coef[it][k]*gdx[k*nsou+js];
	dd += coef[it][k]*gdy[k*nsou+js];
      } /* end loop  L140: */

      /* Calculate residuals */
      rx = dx[i] - soffx[js] - dr;
      ry = dy[i] - soffy[js] - dd;
      resid[numt] = sqrt (rx*rx + ry*ry);

      /* Residual statistics */
      rms += rx*rx + ry*ry;
      rmscnt++;
      numt++;  /* How many found? */
    } 
  } /* end loop  L200: */;

  /* Enough on this source to bother? */
  if (numt <= 5) goto cleanup;

  /* Source RMS */
  rms = sqrt (rms/rmscnt);

  /* Can't be more than target */
  rmst = MIN (rms, MaxRMS);

  /* Loop comparing with 7 point running mean. */
  nt = 0;
  wid = MIN (3, numt/2);
  /* (Best leave this loop 1-rel indexing) */
  for (i= 1; i<=numt; i++) { /* loop 300 */
    sum = 0.0;
    count = 0;
    i1 = i - wid;
    i2 = i + wid;
    if (i1 < 1) {
      i2 = MIN (numt, 2*wid+1);
      i1 = 1;
    } 
    if (i2 > numt) {
      i2 = MIN (numt, i2);
      i1 = MAX (1, i2-2*wid);
    }

    /* Sum points within WID of current  integration but excluding it. */
    for (j= i1; j<=i2; j++) { /* loop 210 */
      idt = abs (ktime[i-1]-ktime[j-1]);
      if ((idt > 0)  &&  (idt <= wid)) {
	sum += resid[j-1];
	count++;
      } 
    } /* end loop  L210: */

    /* Need at least 2 to compare with */
    if (count > 1) {
      test = resid[i-1] - (sum / count);
      /* Toss it if difference greater than 2 sigma. */
      if (test > 2.0*rmst) {
	(*stoss)++;
	nt++;
	w[iobs[i-1]] = 0.0;
      } 
    } else {
      /* Too little to compare - toss */
      (*stoss)++;
      nt++;
      w[iobs[i-1]] = 0.0;
    } 
  } /* end loop  L300: */

  /* Tell results */
  if (prtLv>=1) {
    Obit_log_error(err, OBIT_InfoErr, 
		   "%s: Source %4d rejected %4d of %4d RMS = %10.3f", 
		   routine, source+1, nt, numt, 3600.0*rms);
  }

  /* deallocate arrays */
 cleanup:
  if (resid) g_free(resid);
  if (ktime) g_free(ktime);
  if (iobs)  g_free(iobs);

} /* end of routine IonEditSource */ 

/**
 * Fit Zernike model to apparent position offsets.  
 * Solution uses gsl fitting
 * dxj = apparent offset in x of obs j  
 * dyj = apparent offset in y of obs j  
 * Pi  = Zernike coefficient i [COEF(I) below],  
 * Gxij = Zernike x gradient term i for obs j        
 * Gyij = Zernike y gradient term i for obs j        
 * \param nobs Number of observations  
 * \param isou  Source number per obs.  (0-rel)
 * \param x     Offset in field X (RA) on unit Zernike circle  per source  
 * \param y     Offset in field Y (Dec) on unit Zernike circle per source  
 * \param dx    Apparent RA position shifts on the Zernike plane (deg) 
 * \param dy    Apparent dec position shifts on the Zernike plane (deg) 
 * \param w     Weight per source.
 * \param ncoef [in] Maximum number of coefficients
 *              [out] actual number
 * \param coef  [out] Zernike coefficients fitted 
 * \param prtLv Print level >=3 give diagnostics
 * \param err   Error stack
 * \return RMS residual (deg) -1 on failure
 */
static ofloat
IonFit1 (gint nobs, olong* isou, ofloat* x, ofloat* y, 
	 ofloat* dx, ofloat* dy, ofloat* w, 
	 olong* ncoef, ofloat* coef, olong prtLv, ObitErr* err)
{
  olong      i, j, is, rmscnt, mcoef;
  ofloat out = -1.0;
  ofloat    rms=0.0, dr, dd, rx, ry;
  gchar *routine = "IonFit1";
  
  /* initial values 0 */
  for (i=0; i<*ncoef; i++) coef[i] = 0.0;

  /* How many coeffients can actually  be fitted? */
  rmscnt = 0;
  for (i=0; i<nobs; i++) { /* loop 10 */
    if (w[i] > 0.0) rmscnt++;
  } /* end loop  L10:  */;
  *ncoef = MIN (*ncoef, rmscnt*2);
  if (*ncoef < 17) (*ncoef) = MIN (10, *ncoef);
  if (*ncoef < 10) (*ncoef) = MIN (5, *ncoef);
  if (*ncoef < 5)  (*ncoef) = MIN (2, *ncoef);
  /* Better have something to work  with */
  if (*ncoef  <=  0) return out;
  mcoef = *ncoef;

  /* Use gsl fit */
  FitZernike ((olong)nobs, isou, x, y, dx, dy, w, (olong)mcoef, coef);

  /* Get RMS */
  rms = 0.0;
  rmscnt = 0;
  for (i=0; i<nobs; i++) {
    if (w[i] > 0.0) {
      is = isou[i];
      dr = 0.0;
      dd = 0.0;
      for (j=0; j<mcoef; j++) {
	dr += coef[j]*ObitZernikeGradX(j+2, x[is], y[is]);
	dd += coef[j]*ObitZernikeGradY(j+2, x[is], y[is]);
      } 
      rx = dx[i] - dr;
      ry = dy[i] - dd;
      rms += rx*rx + ry*ry;
      rmscnt++;
      /* Diagnostics */
      if (prtLv>=3) {
	Obit_log_error(err, OBIT_InfoErr, 
		       "%s: obs %5d o=%7.2f %7.2f m=%7.2f %7.2f r=%7.2f %7.2f w=%6.2f",
		       routine, i+1, dx[i]*3600.0, dy[i]*3600.0, dr*3600.0, dd*3600.0, 
		       rx*3600.0, ry*3600.0, w[i]);
      }   /* End Diagnostics */
    } 
  } /* end print loop */ 
  rms = sqrt (rms/rmscnt);
  Obit_log_error(err, OBIT_InfoErr, 
		   "%s: RMS residual=%9.2f asec",routine, 3600.0*rms);
  ObitErrLog(err); /* show any messages on err */

  out = rms;
  return out;
} /* end IonFit1 */ 

/**
 * Make a parabolic least-squares fit to a 9x9 matrix about the
 * peak value in an array and determine strength and position
 * of maximum.
 * 0.5*min(imsi[0],imsi[1])  pixels of the center.  
 * Routine translated from the AIPSish CALPOS.FOR/FXPFIT  
 * \param a       Data input array[dec][RA] 
 * \param dx      Position of max relative to a[5][5]
 * \param s       Strength of max 
 * \param fblank  Value for blanked pixel  
 * \return 0=OK else failed  
 */
static olong pfit (ofloat a[9][9], ofloat *s, ofloat dx[2], ofloat fblank)
{
  float  absmax, x, y, temp[6], momar[6], d;
  float mat[3][3] = {{0.55555, -0.33333, -0.33333}, {-0.33333, 0.5, 0.0},
		     {-0.33333, 0.0, 0.5}};
  int ix, iy, ixmax, iymax;
  /* default return values */
  dx[0] = 0;
  dx[1] = 0;
  *s = a[3][3];
  /*  find peak in array    */
  absmax = 0.0; ixmax = -1; iymax = -1;
  for (ix=0; ix<9; ix++) {
    for (iy=0; iy<9; iy++) {
      if ((a[iy][ix]!=fblank) && (fabs(a[iy][ix])>absmax)) 
	{absmax=fabs(a[iy][ix]); ixmax = ix; iymax = iy;}
    }
  }
  /* check for valid data */
  if ((ixmax<0) || (iymax<0)) return 1;
  /*  00, 01, 02, 10, 11, 20 */
  x = ixmax+1;
  y = iymax+1;
  /*  default values       */
  dx[0] = x - 5.0;
  dx[1] = y - 5.0;
  *s = a[iymax][ixmax];
  if (momnt (a, x, y, 3, 3, momar, fblank)) return 1;

  /*  multiply matrix * even moms  yields const & quadratic terms */
  temp[0] = momar[0];
  temp[1] = momar[2];
  temp[2] = momar[5];
  matvmul (mat, temp, &temp[3], 3);

  /*  pick up linear & cross term  */
  temp[0] = momar[1] / 6.;
  temp[1] = momar[3] / 6.;
  temp[2] = momar[4] / 4.;

  /*  offset of peak */
  d = 4.* temp[4] * temp[5] - (temp[2]*temp[2]);
  if (d==0.0) return 2;
  dx[0] = (temp[2]*temp[0] - 2.*temp[1]*temp[4]) / d;
  dx[1] = (temp[2]*temp[1] - 2.*temp[0]*temp[5]) / d;
  /*  value of peak */
  *s = temp[3] + dx[0]*(temp[1] + dx[0]*temp[5]
			+ dx[1]*temp[2]) + dx[1]*(temp[0]+dx[1]*temp[4]);
  dx[0] = dx[0] + x - 5.0;  /* correct wrt center of input array */
  dx[1] = dx[1] + y - 5.0;
  return 0;
} /* end of pfit */

/**
 * Calculate all 0th, 1st, and 2nd moments of a nx*ny subarray of ara
 * centered at x,y.  nx and ny should be odd.
 * \param ara     Input data array 
 * \param x       x-center for moment calculation (1-rel)
 * \param nx      # of points to include in x-direction. nx  
 *                should be odd.  
 *                The points will be centered about x (rounded)
 * \param ny      # of points in y-direction
 * \param momar   00,10,20,01,11,02 yx-moments of ara  
 * \param fblank  Value for blanked pixel  
 * \return 0=OK else failed  1 => subarray doesn't fit in main array 
 */
static olong 
momnt (ofloat ara[9][9], ofloat x, ofloat y, olong nx, olong ny, 
       ofloat momar[6], ofloat fblank)
  {
    olong  ind, i, j, k, nj, indx, indy, iax1, iax2, iay1, iay2;
    ofloat s, t, arg, prod;
    
    /*      compute loop limits (1-rel)   */
    i = x + 0.5;
    iax1 = i - nx/2;
    iax2 = iax1 + nx - 1;
    i = y + 0.5;
    iay1 = i - ny/2;
    iay2 = iay1+ ny - 1;

    /*      check loop limits             */
    if ((iax1<1) || (iax2>9) || (iay1<1) || (iay2>9)) return 1;
    /*      compute moments               */
    ind = 0;
    for (i = 1; i<=3; i++) {
      nj = 4 - i;
      for (j=1; j<=nj; j++) {
	ind = ind + 1;
	s = 0.0;
	for (indx=iax1; indx<=iax2; indx++) {
	  for (indy=iay1; indy<=iay2; indy++) {
	    t = ara[indy-1][indx-1];
	    if (t!=fblank) {
	      if (i>1) {
		prod = 1.0; 
		arg = indx - x;
		for (k=1; k<i; k++) prod *= arg;
		t = t * prod;}     /* ((indx - x)**(i-1)); */
	      if (j>1) {
		prod = 1.0; 
		arg = indy - y;
		for (k=1; k<j; k++) prod *= arg;
		t = t * prod;}  /* ((indy - y)**(j-1));*/
	      s = s + t;}
	  }
	}
	momar[ind-1] = s;
      }
    }
    return 0;
  }  /* end of momnt */

/**
 *  Matrix-vector multiplication  vo = vi * m  
 * \param  m      Input matrix      
 * \param vi      Input vector 
 * \param n       Array dimension  
 * \param vo      [out] Output vector
 */
static void matvmul (ofloat m[3][3], ofloat vi[3], ofloat vo[3], olong n)
{
  int  i, j;
  float s;
  
  for (i=0; i<n; i++) {
    s = 0.0;
    for (j=0; j<n; j++) s = s + m[j][i] * vi[j];
    vo[i] = s;
  }
}  /* end of matvmul */


/**
 * Converts from celestial coordinates, expressed in terms of  
 * a reference position and a shift from this position to coordinates  
 * in a plane for fitting an Ionospheric phase screen model  
 * consisting of Zernike polynomials.  The output coordinates are  
 * normalized to unity at a 10 deg radius from the reference position.  
 * The coordinates are projected onto a plane tangent to the sky at  
 * the reference position.  
 * Routine translated from the AIPSish ZERGEOM.FOR/RD2ZER 
 * \param ra      Right Ascention of reference position (deg) 
 * \param dec     Declination of reference position (deg) 
 * \param xshift  Shift in X (RA) to desired position (deg) 
 * \param yshift  Shift in Y (Dec) to desired position (deg) 
 * \param xzer    [out] x-coordinate on Zernike plane 
 * \param yzer    [out]  y-coordinate on Zernike plane 
 * \param ierr    0 ok, 1 out of range 
 */
static void rd2zer (odouble ra, odouble dec, ofloat xshift, ofloat yshift, 
		    ofloat* xzer, ofloat* yzer, olong *ierr) 
{
  odouble dx, dy, a, b, ax, ay, xxinc, coss, sins, sint, dt;
  odouble pi, twopi, zernor;

  /* Constants */
  pi = 3.14159265358979323846e0;
  twopi = 2.0e0*pi;
  /* Zernike normalization to Unit circle (10 deg) */
  zernor = 1.0 / (10.0e0 * DG2RAD);

  /* initial values */
  *xzer = 0.0;
  *yzer = 0.0;
  *ierr = 0;
  
  /* Convert to radians */
  a = DG2RAD * ra;
  b = DG2RAD * dec;
  /* Convert shift to position */
  xxinc = cos (DG2RAD * dec);
  if (xxinc != 0) {
    ax = ra + xshift / xxinc;
  } else {
    ax = ra;
  } 
  ay = dec + yshift;
  ax = ax * DG2RAD;
  ay = ay * DG2RAD;
  
  /* Get projection cosines */
  *ierr = 1;
  if (fabs(ay) > pi/2.0e0) return;
  if (fabs(b) > pi/2.0e0) return;
  *ierr = 0;
  coss = cos (ay);
  sins = sin (ay);
  dt = ax - a;
  if (dt > pi) dt = dt - twopi;
  if (dt < -pi) dt = dt + twopi;
  dx = sin (dt) * coss;
  sint = sins * sin (b) + coss * cos (b) * cos (dt);
  
  /* SIN projection */
  /*      IF (SINT.LT.0.0D0) IERR = 2 */
  if (sint < 0.0e0) {
    *ierr= 1;
    return;
  } 
  dy = sins * cos (b) - coss * sin (b) * cos (dt);

  /* Normalize to Zernike unit sphere */
  *xzer = dx * zernor;
  *yzer = dy * zernor;
  return;
} /* end of routine rd2zer */ 

/**
 * Given an offset from a given field center, return the rotation of  
 * the coordinates of the offset wrt the field center.  
 * This is just the difference in RA.  
 * Routine translated from the AIPSish ZERGEOM.FOR/ZERROT  
 * \param ra   Initial RA in degrees: ref point 
 * \param dec   Initial Declination in degrees: ref point 
 * \param xshift   Shift in "X" in degrees - simple shifts not -SIN 
 * \param yshift   Shift in "Y" in degrees 
 * \param rotate   [out] Image rotation in degrees 
 */
static void zerrot (odouble ra, odouble dec, ofloat xshift, ofloat yshift, 
		    ofloat* rotate) 
{
  odouble xxinc, xra, xdec;
  
  /* simple shifts */
  xdec = dec + yshift;
  xxinc = cos (DG2RAD * dec);
  if (xxinc != 0) {
    xra = ra + xshift / xxinc;
  } else {
    xra = ra;
       } 
  /* Difference in RA */
  (*rotate) = (ra - xra);
  /* Fold to range (-180,+180) */
  if (*rotate > +180.0e0) (*rotate) -= 360.0e0;
  if (*rotate < -180.0e0) (*rotate) += 360.0e0;
} /* end of routine zerrot */ 

/**
 * Convert a timerange as time in days to a printable string
 * \param    timer  Beginning time, end time in days
 * \msgBuff  Human readable string as "dd/hh:mm:ss.s-dd/hh:mm:ss.s"
 *           must be allocated at least 30 characters
 */
static void TR2String (ofloat timer[2], gchar *msgBuf)
{
  ofloat rtemp, rt1, rt2;
  olong   id1, id2, it1, it2, it3, it4;

  id1 = timer[0];
  rtemp = 24.0 * (timer[0] - id1);
  it1 = rtemp;
  it1 = MIN (23, it1);
  rtemp = (rtemp - it1)*60.0;
  it2 = rtemp;
  it2 = MIN (59, it2);
  rt1 = (rtemp - it2)*60.0;
  id2 = timer[1];
  rtemp = 24.0 * (timer[1] - id2);
  it3 = rtemp;
  it3 = MIN (23, it3);
  rtemp = (rtemp - it3)*60.0;
  it4 = rtemp;
  it4 = MIN (59, it4);
  rt2 = (rtemp - it4)*60.0;
  g_snprintf (msgBuf, 30, "%2.2d/%2.2d:%2.2d:%5.2f-%2.2d/%2.2d:%2.2d:%5.2f",
	      id1, it1, it2, rt1, id2, it3, it4, rt2);
} /* end of routine TR2String   */ 

/*                 CalListElem functions             */

/**
 * CalListElem Constructor
 * \param ra      Catalog celestial position RA (deg)
 * \param dec     Catalog celestial position RA (deg)
 * \param ZernXY  Offset on Zernike plane (deg)
 * \param shift   Shift from reference position
 * \param pixel   Expected pixel in reference image
 * \param flux    Estimated catalog flux density
 * \param offset  Measured offset from expected position (pixels)
 * \param peak    Measured peak flux density (image units)
 * \param fint    Measured integrated flux density (image units)
 * \param wt      Determined weight
 * \param qual    Catalog quality code 
 * \param calNo   Calibrator number
 * \param epoch   Epoch number
 * \return the new  object.
 */
static CalListElem* 
newCalListElem (odouble ra, odouble dec, ofloat ZernXY[2], ofloat shift[2],
		ofloat pixel[2], ofloat flux, ofloat offset[2],
		ofloat peak, ofloat fint, ofloat wt,  olong qual,
		gint calNo, olong epoch)
{
  CalListElem *out=NULL;

  out = g_malloc0(sizeof(CalListElem));
  out->ra        = ra;
  out->dec       = dec;
  out->ZernXY[0] = ZernXY[0];
  out->ZernXY[1] = ZernXY[1];
  out->shift[0]  = shift[0];
  out->shift[1]  = shift[1];
  out->pixel[0]  = pixel[0];
  out->pixel[1]  = pixel[1];
  out->flux      = flux;
  out->offset[0] = offset[0];
  out->offset[1] = offset[1];
  out->peak      = peak;
  out->fint      = fint;
  out->wt        = wt;
  out->qual      = qual;
  out->calNo     = calNo;
  out->epoch     = epoch;
  return out;
} /* end newCalListElem */

/**
 * Update CalListElem
 * \param elem    CalListElem to update
 * \param ra      Catalog celestial position RA (deg)
 * \param dec     Catalog celestial position Dec (deg)
 * \param ZernXY  Offset on Zernike plane (deg)
 * \param shift   Shift from reference position
 * \param pixel   Expected pixel in reference image
 * \param flux    Estimated catalog flux density
 * \param offset  Measured offset from expected position (pixels)
 * \param peak    Measured peak flux density (image units)
 * \param fint    Measured integrated flux density (image units)
 * \param wt      Determined weight
 * \param qual    Catalog quality code 
 * \param calNo   Calibrator number
 * \param epoch   Epoch number
 * \return the new  object.
 */
static void
CalListElemUpdate (CalListElem *elem, odouble ra, odouble dec, 
		   ofloat ZernXY[2], ofloat shift[2], ofloat pixel[2], 
		   ofloat flux, ofloat offset[2], ofloat peak, ofloat fint, 
		   ofloat wt, olong qual, olong calNo, olong epoch)
{
  elem->ra        = ra;
  elem->dec       = dec;
  elem->ZernXY[0] = ZernXY[0];
  elem->ZernXY[1] = ZernXY[1];
  elem->shift[0]  = shift[0];
  elem->shift[1]  = shift[1];
  elem->pixel[0]  = pixel[0];
  elem->pixel[1]  = pixel[1];
  elem->flux      = flux;
  elem->offset[0] = offset[0];
  elem->offset[1] = offset[1];
  elem->peak      = peak;
  elem->fint      = fint;
  elem->wt        = wt;
  elem->qual      = qual;
  elem->calNo     = calNo;
  elem->epoch     = epoch;
} /* end CalListElemUpdate */

/**
 * Print contents 
 * \param in Object to print
 * \param file  FILE* to write to
 */
static void CalListElemPrint (CalListElem *in, FILE *file)
{
  if (!in) return;

  fprintf (file, "Cal. no.=%d, epoch=%d\n",in->calNo, in->epoch);
  fprintf (file, "RA=%lf, Dec=%lf\n",in->ra, in->dec);
  fprintf (file, "Offset (deg) on Zernike plane [%f,%f]\n",
	   in->ZernXY[0],in->ZernXY[1]);
  fprintf (file, "Position shift (deg) [%f,%f]\n",
	   in->shift[0],in->shift[1]);
  fprintf (file, "Expected pixel [%f,%f]\n",
	   in->pixel[0],in->pixel[1]);
  fprintf (file, "Expected flux density %f, quality %d\n",
	   in->flux,in->qual);
  fprintf (file, "Measured position offset (asec) [%f,%f]\n",
	   in->offset[0]*3600.0,in->offset[1]*3600.0);
  fprintf (file, "Measured peak %f, integ. %f, wt %f\n",
	   in->peak,in->fint, in->wt);

} /* end CalListElemPrint */

/**
 * Destructor 
 * \param in Object to delete
 */
static void freeCalListElem (CalListElem *in)
{
  if (in) g_free(in);
} /* end freeCalListElem */



/*  CalList functions */
/**
 * CalList Constructor
 * \return the new  CalList structure.
 */
static CalList* newCalList (void)
{
  CalList *out=NULL;

  out = g_malloc0(sizeof(CalList));
  out->number = 0;
  out->list   = NULL;
  return out;
} /* end newCalListElem */

/**
 * Attach elem to list in
 * \param in   list to add elem to
 * \param elem the element to add.
 */
static void CalListAdd (CalList *in, CalListElem *elem)
{
  /* link to list */
  in->list = g_slist_append (in->list, elem);
  in->number++;
} /* end CalListAdd */

/**
 * Remove elem from list in
 * \param in   list to remove elem from
 * \param elem the element to remove.
 */
static void CalListRemove (CalList *in, CalListElem *elem)
{
  /* remove from list */
  in->list = g_slist_remove(in->list, elem);
  in->number--; /* keep count */  
} /* end CalListRemove  */

/**
 * Remove all elements from list in
 * \param in   list to remove elem from
 * \param elem the element to remove.
 */
static void CalListClear (CalList *in)
{
  GSList *tmp;

  if (in==NULL) return;  /* Does it exist? */
  if (in->list==NULL) return;  /* Anything in it? */

  /* loop through list deleting elements */
  tmp = in->list;
  while (tmp!=NULL) {
    if (tmp->data) freeCalListElem(tmp->data);
    tmp = g_slist_next(tmp);
  }

  /* delete members  */
  g_slist_free(in->list);
  in->list = NULL;
  in->number = 0;

} /* end CalListClear  */

/**
 * Print all elements in list in to file
 * \param in   list to remove elem from
 * \param elem the element to remove.
 */
static void CalListPrint (CalList *in, FILE *file)
{
  GSList *tmp;

  if (in==NULL) return;  /* Does it exist? */
  if (in->list==NULL) return;  /* Anything in it? */

  fprintf (file, "Listing Of CalList\n");

  /* loop through list deleting elements */
  tmp = in->list;
  while (tmp!=NULL) {
    if (tmp->data) CalListElemPrint(tmp->data, file);
    tmp = g_slist_next(tmp);
  }

} /* end CalListPrint  */

/**
 * Destructor 
 * \param in Object to delete
 */
static void freeCalList (CalList *in)
{
  /* Clear List */
  CalListClear(in);

  /* delete object */
  g_free(in);
} /* end freeCalListElem */

/**
 * Fit Zernike polynomial to position offset pairs
 * Use gsl package.
 * \param nobs    Number of data measurements
 * \param isou    Source number per obs.  (0-rel)
 * \param x       Zernike X on unit circle per source
 * \param y       Zernike Y on unit circle per source
 * \param dx      X offset (deg)
 * \param dy      X offset (deg)
 * \param wt      Data weights
 * \param nZern   Number of Zernike coefficients to fit
 * \param coef    [out] Zernike coefficients
 */
void  FitZernike (olong nobs, olong* isou, ofloat *x, ofloat *y, 
		  ofloat *dx, ofloat *dy, ofloat *wt, 
		  olong nZern, ofloat *coef)
{
#if HAVE_GSL==1  /* GSL stuff */
  olong i, j, k, is, good, p=nZern, npen;
  double xi, chisq, sumwt, penWt;
  gsl_matrix *X, *cov;
  gsl_vector *yy, *w, *c;
  gsl_multifit_linear_workspace *work;
  /* Penalty terms by order of Zernike 1-4 */
  ofloat pen[]={0.001, 0.001,
		0.01, 0.01, 0.01,
		0.03, 0.03, 0.03, 0.03, 0.03, 
		0.05, 0.05, 0.05,0.05,  0.05, 0.05, 0.05, 
		0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08, 0.08};

  /* Only use good data (2 per measurement) */
  good = 0;sumwt = 0.0;
  for (i=0; i<nobs; i++) if (wt[i]>0.0) 
    {good += 2; sumwt += wt[i];}

  /* Number of penalty terms */
  npen = nZern;

  /* allocate arrays */
  X    = gsl_matrix_alloc(good+npen, p);
  yy   = gsl_vector_alloc(good+npen);
  w    = gsl_vector_alloc(good+npen);
  c    = gsl_vector_alloc(p);
  cov  = gsl_matrix_alloc(p, p);
  work = gsl_multifit_linear_alloc (good+npen, p);

  /* set data */
  k = 0;
  for (i=0; i<nobs; i++) {
    if (wt[i]>0.0) {
      is = isou[i];
      /* X offset */
      gsl_vector_set(yy, k, dx[i]);
      gsl_vector_set(w, k, wt[i]);
      for (j=0; j<p; j++) {
	xi = ObitZernikeGradX(j+2, x[is], y[is]);
	gsl_matrix_set(X, k, j, xi);
      }
      k++; 
      /* Y offset */
      gsl_vector_set(yy, k, dy[i]);
      gsl_vector_set(w, k, wt[i]);
      for (j=0; j<p; j++) {
	xi = ObitZernikeGradY(j+2, x[is], y[is]);
	gsl_matrix_set(X, k, j, xi);
      }
      k++; 
    }
  }

  /* Penalty terms - try to get coefficients as small as possible */
  penWt = sumwt / good;  /* Penalty weight data dependency */
  for (i=0; i<npen; i++) {
    gsl_vector_set(yy, k, 0.0);
    gsl_vector_set(w, k, penWt*pen[i]);
    for (j=0; j<p; j++) {
      if (j==i) xi = 1.0;
      else xi = 0.0;
      gsl_matrix_set(X, k, j, xi);
    }
    k++; 
  } /* end penalty loop */

  /* Fit */
  gsl_multifit_wlinear (X, w, yy, c, cov, &chisq, work);

  /* get results */
  for (j=0; j<nZern; j++) coef[j] = gsl_vector_get(c, j);

  /* Deallocate arrays */
  gsl_matrix_free(X);
  gsl_vector_free(yy);
  gsl_vector_free(w);
  gsl_vector_free(c);
  gsl_matrix_free(cov);
  gsl_multifit_linear_free (work);
#else  /* No GSL - stubb */
  g_error ("FitZernike: GSL not available - cannot do fit");
#endif /* GSL stuff */
} /* end FitZernike */

/**
 * Called when CLEAN failed to mark first 5 entries as failed.
 * \param in   IonCal object 
 * \param mosaic   Image mosaic object 
 * \param calList Calibrator list each element of which has:
 * \li ra      position RA (deg) of calibrators
 * \li dec     position Dec (deg) of calibrators
 * \li shift   Offset in field [x,y] on unit Zernike circle of position
 * \li pixel   Expected pixel in reference image
 * \li flux    Estimated catalog flux density 
 * \li offset  Measured offset [x,y] from expected position (deg)
 * \li peak    set to 0.0
 * \li fint    set to 0.0
 * \li wt      set to 0.0
 * \li qual    Catalog quality code
 * \param err     Error stack
 * \param epoch   Epoch number
 */
static void CalBadTime (ObitIonCal *in, ObitImageMosaic* mosaic, 
			gint epoch, ObitErr* err)  
{
  olong   field;
  ofloat offset[2], peak, fint, wt;
  ObitImage *image=NULL;
  CalListElem *elem=NULL, *nelem=NULL;
  GSList  *tmp;
  CalList *calList;
  /*gchar *routine = "CalBadTime";*/

  /* Error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return ;  /* previous error? */
  g_assert(ObitIonCalIsA(in));
  g_assert(ObitImageMosaicIsA(mosaic));

  calList = (CalList*)in->calList;

  /* do up to 5 fields Loop over images  */
  for (field=0; field<MIN (5, mosaic->numberImages); field++) {
    image = mosaic->images[field];
 
    /* Find cal info */
    tmp = calList->list;
    while (tmp!=NULL) {   /* loop 100 */
      elem = (CalListElem*)tmp->data;
      if (elem->calNo == (field+1)) break;  /* found it */
      tmp = g_slist_next(tmp);  /* next item in list */
    } /* end loop  L100: over calibrator list */

    /* zero it out */
    offset[0] = offset[1] = 0.0;
    peak = 0.0;
    fint = 0.0;
    wt   = 0.0;

    /* Update CalList */
    nelem = newCalListElem (elem->ra, elem->dec, elem->ZernXY, elem->shift, 
			    elem->pixel, elem->flux, offset, peak, fint, wt, 
			    elem->qual, field+1, epoch);
    CalListAdd (calList, nelem);
  } /* end loop over fields */
  
} /* end of routine CalBadTime */ 

