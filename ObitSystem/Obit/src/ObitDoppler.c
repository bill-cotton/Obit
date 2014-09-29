/* $Id$          */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2012-2014                                          */
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
/*; Correspondence about this software should be addressed as follows:*/
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
#include "ObitDoppler.h"
#include "ObitTableANUtil.h"
#include "ObitTableSUUtil.h"
#include "ObitPrecess.h"
#ifndef VELIGHT
#define VELIGHT 2.997924562e8
#endif

/** Name of the class defined in this file */
static gchar *myClassName = "ObitDoppler";

/** Function to obtain parent ClassInfo */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo structure ObitDopplerClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitDopplerClassInfo myClassInfo = {FALSE};

/**
 * \file ObitDoppler.c
 * ObitDoppler  Source Doppler shifting
 */
/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
static void  ObitDopplerInit  (gpointer in);

/** Private: Deallocate members. */
static void  ObitDopplerClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitDopplerClassInfoDefFn (gpointer inClass);

/** Private: Calculates the velocity component of the observer . */
odouble VelLSR (odouble ra, odouble dec, olong year, olong doy, 
		odouble ut, odouble x, odouble y, odouble z, odouble* vsun);  
/** Private: coordinate transformation */
static void coordTran (odouble ao, odouble bo, odouble ap, odouble bp, odouble
		a1, odouble b1, odouble* a2, odouble* b2); 
/** Private: Calculates the geodetic latitude and logitude*/
static void GeoDll (odouble x, odouble y, odouble z, 
	     odouble *lat, odouble  *lng);
/** Private: */
static void rdmove (olong  nyri, olong nyrf, olong mo, olong nda, 
	     odouble ra, odouble d, 
	     odouble* delr, odouble* deld, odouble* dc ); 
/** Private: Compute JD for beginning of a given year */
static olong julDay (olong nyr);
/** Frequency shift single autocorrelation function */
static void ACShift (ObitCArray  *Spectrum,  ObitCArray  *Work,  
		     ObitFFT *FFTFor, ObitFFT *FFTRev, 
		     olong sideband, olong nchan, ofloat shift, olong doSmo);
/** Frequency shift single crosscorrelation function */
static void CCShift (ObitCArray *Spectrum,  ObitCArray *Work,
		     ObitFFT *FFTFor, ObitFFT *FFTRev, 
		     olong sideband, olong nchan, ofloat shift);
/** How many channels to shift? */
static ofloat CalcShift (ObitDoppler *doppler, olong ant1, olong ant2, 
			 olong souId, ofloat time, ObitErr *err);
/** Update SU table for multi source */
static void UpdateSUTable(ObitDoppler *doppler, ObitUV *inData, ObitErr *err);
/** Get default reference channel */
static ofloat defaultRefChan(ObitDoppler *doppler, ObitErr *err);
/** Days to human string */
static void day2dhms(ofloat time, gchar *timeString);
/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitDoppler* newObitDoppler (gchar* name)
{
  ObitDoppler* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitDopplerClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitDoppler));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitDopplerInit((gpointer)out);

 return out;
} /* end newObitDoppler */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitDopplerGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitDopplerClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitDopplerGetClass */

/**
 * Make a deep copy of an ObitDoppler.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitDoppler* ObitDopplerCopy  (ObitDoppler *in, ObitDoppler *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  gboolean oldExist;
  gchar *outName;

  /* error checks */
  if (err->error) return out;
  g_assert (ObitIsA(in, &myClassInfo));
  if (out) g_assert (ObitIsA(out, &myClassInfo));

  /* Create if it doesn't exist */
  oldExist = out!=NULL;
  if (!oldExist) {
    /* derive object name */
    outName = g_strconcat ("Copy: ",in->name,NULL);
    out = newObitDoppler(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  out->uvdata   = ObitUVUnref(out->uvdata);
  out->uvdata   = ObitUVRef(in->uvdata);
  out->Spectrum = ObitCArrayCopy(in->Spectrum, out->Spectrum, err);
  out->Work     = ObitCArrayCopy(in->Work, out->Work, err);

  return out;
} /* end ObitDopplerCopy */

/**
 * Make a copy of a object but do not copy the actual data
 * This is useful to create an Doppler similar to the input one.
 * \param in  The object to copy
 * \param out An existing object pointer for output, must be defined.
 * \param err Obit error stack object.
 */
void ObitDopplerClone  (ObitDoppler *in, ObitDoppler *out, ObitErr *err)
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
  out->uvdata   = ObitUVUnref(out->uvdata);
  out->uvdata   = ObitUVRef(in->uvdata);
  out->Spectrum = ObitCArrayCopy(in->Spectrum, out->Spectrum, err);
  out->Work     = ObitCArrayCopy(in->Work, out->Work, err);

} /* end ObitDopplerClone */

/**
 * Creates an ObitDoppler 
 * \param name  An optional name for the object.
 * \return the new object.
 */
ObitDoppler* ObitDopplerCreate (gchar* name, ObitUV *uvdata, ObitErr *err)
{
  ObitDoppler* out=NULL;
  ObitTableAN *ANTab=NULL;
  ObitTableSU *SUTab=NULL;
  gboolean doCalSelect;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  ObitIOAccess access;
  ObitIOCode iretCode;
  ObitUVDesc *inDesc;
  olong ndim=1, naxis[1], iANver, iSUver;
  gchar *routine = "ObitDopplerCreate";

  if (err->error) return out;  /* Existing error? */

  /* Create basic structure */
  out = newObitDoppler (name);

  /* Get data */
  out->uvdata = ObitUVRef(uvdata);

  /* Open/close to see how big things are */
    /* Selection of input? */
  doCalSelect = FALSE;
  ObitInfoListGetTest(uvdata->info, "doCalSelect", &type, dim, &doCalSelect);
  if (doCalSelect) access = OBIT_IO_ReadCal;
  else access = OBIT_IO_ReadOnly;

  /* test open to fully instantiate input and see if it's OK */
  iretCode = ObitUVOpen (uvdata, access, err);
  iretCode = ObitUVClose (uvdata, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) 
    Obit_traceback_val (err, routine, uvdata->name, out);
  inDesc = uvdata->myDesc;

  /* Sizes of things */
  out->nchan = inDesc->inaxes[inDesc->jlocf];
  if (inDesc->jlocif>=0) out->nif = inDesc->inaxes[inDesc->jlocif];
  else                   out->nif = 1;
  if (inDesc->jlocs>=0) out->npoln = inDesc->inaxes[inDesc->jlocs];
  else                  out->npoln= 1;

  /* Create work arrays - need double sideband */
  naxis[0] = 2*ObitFFTSuggestSize(out->nchan);
  out->Spectrum = ObitCArrayCreate("Spectrum", ndim, naxis);
  out->Work     = ObitCArrayCreate("Work", ndim, naxis);

  /* Create FFTs - pad as necessary */
  naxis[0] = 2*ObitFFTSuggestSize (out->nchan);
  out->FFTFor = newObitFFT ("For", OBIT_FFT_Forward, OBIT_FFT_FullComplex,
			     ndim, naxis);
  out->FFTRev = newObitFFT ("Rev", OBIT_FFT_Reverse, OBIT_FFT_FullComplex,
			     ndim, naxis);

  /* Antenna list */
  iANver = 1;
  ANTab = 
    newObitTableANValue (out->name, (ObitData*)uvdata, &iANver, OBIT_IO_ReadOnly, 0, 0, 0, err);
  out->antList = ObitTableANGetList(ANTab, err);
  if (err->error) Obit_traceback_val (err, routine, out->name, out);
  ANTab = ObitTableANUnref(ANTab);

  /* Source */
  out->source = newObitSource("src");
  out->source->SourID = -999;

  /* Source list */
  iSUver = 1;
  SUTab = 
    newObitTableSUValue (out->name, (ObitData*)uvdata, &iSUver, OBIT_IO_ReadOnly, 0, err);
  out->sourceList = ObitTableSUGetList(SUTab, err);
  if (err->error) Obit_traceback_val (err, routine, out->name, out);
  SUTab = ObitTableSUUnref(SUTab);

  return out;
} /* end ObitDopplerCreate */

/**
 * Calculates frequency corresponding to a given LSR velocity
 * for a given line, direction, observer location and date.
 *  
 * This version takes into account components of the observer's  
 * motion due to the rotation of the earth, the revolution of the  
 * earth-moon barycenter about the sun, and the motion of the earth's  
 * centre about the earth-moon barycenter.  The perturbations of the  
 * earth's orbit due to planets are neglected.  The absolute precision  
 * of this version of dop is about 0.004 km/sec, but since the dominant  
 * error term is slowly varying the relative error will be considerably  
 * less for times up to a week or so.  
 *  
 * The routine omits the planetary perturbations on the earth's orbit.  
 * They amount to about 0.003 km/sec and are thought to be the largest  
 * contribution to the error in the velocity.  
 * \param rest  Rest frequency of line
 * \param vlsr  Desired LSR velocity (km/s)
 * \param ra    Right Ascension of date (deg) 
 * \param dec   Declination of date (deg) 
 * \param year  Year of observation 
 * \param doy   Day number of observation 
 * \param ut    UT of datum (days) 
 * \param x     X Geocentric station coord (m)
 * \param y     Y Geocentric station coord (m)  left hand?
 * \param z     Z Geocentric station coord (m)
 * \return      Frequency corresponding the vlsr (Hz)
 */
odouble ObitDopplerFreqLSR (odouble rest, ofloat vlsr, 
			    odouble ra, odouble dec, 
			    olong year, olong doy, odouble ut, 
			    odouble x, odouble y, odouble z)
{
  odouble flsr, vsun, vel;

  /* Velocity of observer in direction of source */
  vel = VelLSR (ra*DG2RAD, dec*DG2RAD, year, doy, ut, x, y, z, &vsun);

  /* Cooresponding line frequency */
  flsr = rest * (VELIGHT/(VELIGHT+vel+(vlsr*1.0e3)));
  return flsr;

} /* end  ObitDopplerFreqLSR */

/**
 * Correct a UV data set for the earth's motion.
 * This version takes into account components of the observer's  
 * motion due to the rotation of the earth, the revolution of the  
 * earth-moon barycenter about the sun, and the motion of the earth's  
 * centre about the earth-moon barycenter.  The perturbations of the  
 * earth's orbit due to planets are neglected.  The absolute precision  
 * of this version is about 0.004 km/sec, but since the dominant  
 * error term is slowly varying the relative error will be considerably  
 * less for times up to a week or so.  
 * The routine omits the planetary perturbations on the earth's orbit.  
 * They amount to about 0.003 km/sec and are thought to be the largest  
 * contribution to the error in the velocity.  
 * Velocity information is updated in output header and SU table if multisource).
 * \param inUV     Input uv data to correct, 
 *                 Any request for calibration, editing and selection honored
 *  Control parameter on info
 * \li JDref    OBIT_double (1) JD of reference velocity.
 * \li RestFreq OBIT_double (1) Rest frequency (Hz) of line
 * \li VelLSR   OBIT_double (1) Desired LSR velocity (m/s) [def 0]
 * \li refChan  Obit_float  (1) Desired channel for VelLSR [def nx/2+1]
 * \param scratch  True if scratch file desired, will be same type as inUV.
 * \param outUV    If not scratch, then the previously defined output file
 *                 May be NULL for scratch only
 *                 If it exists and scratch, it will be Unrefed
 * \param err      Error stack, returns if not empty.
 * \return the frequency corrected ObitUV.
 */
ObitUV* ObitDopplerCVel (ObitUV *inUV, gboolean scratch, ObitUV *outUV, 
			 ObitErr *err)
{
  ObitIOCode iretCode, oretCode;
  gboolean doCalSelect;
  gchar *exclude[]={"AIPS CL", "AIPS SN", "AIPS FG", "AIPS CQ", "AIPS WX",
		    "AIPS AT", "AIPS CT", "AIPS OB", "AIPS IM", "AIPS MC",
		    "AIPS PC", "AIPS NX", "AIPS TY", "AIPS GC", "AIPS HI",
		    "AIPS PL", "AIPS NI", "AIPS BP", "AIPS OF", "AIPS PS",
		    "AIPS FQ", "AIPS SU", "AIPS AN", "AIPS PD",
		    NULL};
  ObitDoppler *doppler=NULL;
  gchar *sourceInclude[] = {"AIPS SU", NULL};
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  union ObitInfoListEquiv InfoReal; 
  olong numAnt, numBL;
  ObitIOAccess access;
  ObitUVDesc *inDesc, *outDesc;
  olong suba, lastSourceID, curSourceID, lastSubA;
  gchar *today=NULL;
  ofloat cbase,  *tmpVis=NULL, uvwScale=1.0;
  ofloat *inBuffer, *outBuffer, *Spectrum, shift;
  olong j, ant1, ant2,ivis=0, iindx=0, oindx=0, NPIO, oldNPIO;
  olong iif, nif, istok, nstok, ichan, nchan, tndx;
  olong incif, incf, incs, oincif, oincf, oincs, off, toff;
  olong ifoff, soff, foff, oifoff, osoff, ofoff, nchanp, sideband;
  gboolean done, gotOne;
  gchar *routine = "ObitDopplerCVel";
 
  /* error checks */
  if (err->error) return outUV;
  g_assert (ObitUVIsA(inUV));
  if (!scratch && (outUV==NULL)) {
    Obit_log_error(err, OBIT_Error,"%s Output MUST be defined for non scratch files",
		   routine);
      return outUV;
  }

  /* Selection/calibration/editing of input? */
  doCalSelect = FALSE;
  ObitInfoListGetTest(inUV->info, "doCalSelect", &type, dim, &doCalSelect);
  if (doCalSelect) access = OBIT_IO_ReadCal;
  else access = OBIT_IO_ReadOnly;

   /* Copy number of records per IO to output */
  ObitInfoListGet (outUV->info, "nVisPIO", &type, dim,  &oldNPIO, err);
  ObitInfoListGet (inUV->info, "nVisPIO", &type, dim,  &NPIO, err);
  ObitInfoListAlwaysPut (outUV->info, "nVisPIO",  type, dim,  &NPIO);

  /* test open to fully instantiate input and see if it's OK */
  iretCode = ObitUVOpen (inUV, access, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine, inUV->name, outUV);

  /* Create scratch? */
  if (scratch) {
    if (outUV) outUV = ObitUVUnref(outUV);
    outUV = newObitUVScratch (inUV, err);
  } else { /* non scratch output must exist - clone from inUV */
    ObitUVClone (inUV, outUV, err);
  }
  if (err->error) Obit_traceback_val (err, routine, inUV->name, inUV);

  /* copy Descriptor */
  outUV->myDesc = ObitUVDescCopy(inUV->myDesc, outUV->myDesc, err);
  inBuffer = inUV->buffer;  /* Local copy of buffer pointer */
 
  /* Output creation date today */
  today = ObitToday();
  strncpy (outUV->myDesc->date, today, UVLEN_VALUE-1);
  if (today) g_free(today);
  
  suba    = 1;
  numAnt  = inUV->myDesc->numAnt[suba-1];/* actually highest antenna number */
  numBL   = (numAnt*(numAnt+1))/2;  /* Include auto correlations */

  /* test open output */
  oretCode = ObitUVOpen (outUV, OBIT_IO_WriteOnly, err);
  /* If this didn't work try OBIT_IO_ReadWrite */
  if ((oretCode!=OBIT_IO_OK) || (err->error)) {
    ObitErrClear(err);
    oretCode = ObitUVOpen (outUV, OBIT_IO_ReadWrite, err);
  }
  /* if it didn't work bail out */
  if ((oretCode!=OBIT_IO_OK) || (err->error)) goto cleanup;

  /* Get descriptors */
  inDesc  = inUV->myDesc;
  outDesc = outUV->myDesc;

  /* sizes of things */
  nchan = inDesc->inaxes[inDesc->jlocf];
  if (inDesc->jlocif>=0) nif = inDesc->inaxes[inDesc->jlocif];
  else                   nif = 1;
  if (inDesc->jlocs>=0) nstok = inDesc->inaxes[inDesc->jlocs];
  else                  nstok= 1;

   /* get increments  */
  incs   = inDesc->incs;
  incf   = inDesc->incf;
  incif  = inDesc->incif;
  oincs  = outDesc->incs;
  oincf  = outDesc->incf;
  oincif = outDesc->incif;

  /* Copy tables before data */
  iretCode = ObitUVCopyTables (inUV, outUV, exclude, NULL, err);
  /* If multisource out then copy SU table, multiple sources selected or
   sources deselected suggest MS out */
  if ((inUV->mySel->numberSourcesList>1) || (!inUV->mySel->selectSources))
  iretCode = ObitUVCopyTables (inUV, outUV, NULL, sourceInclude, err);
  if (err->error) goto cleanup;

  /* reset to beginning of uv data */
  iretCode = ObitIOSet (inUV->myIO,  inUV->info, err);
  oretCode = ObitIOSet (outUV->myIO, outUV->info, err);
  if (err->error) goto cleanup;

  /* Close and reopen input to init calibration which will have been disturbed 
     by the table copy */
  iretCode = ObitUVClose (inUV, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) goto cleanup;

  /* create Doppler object */
  doppler = ObitDopplerCreate ("Doppler", inUV, err);
  if (err->error) goto cleanup;
  /* Work array pointer */
  Spectrum = doppler->Spectrum->array;
  nchanp   = doppler->Spectrum->naxis[0]/2;  /* How big FFT? double sideband array */

  iretCode = ObitUVOpen (inUV, access, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) goto cleanup;

  /* Initialize things */
  lastSourceID = -1;
  curSourceID  = 0;
  outDesc->numVisBuff = 0;
  inBuffer  = inUV->buffer;   /* Local copy of buffer pointer */
  outBuffer = outUV->buffer;  /* Local copy of buffer pointer */

  /* Get control from inData - accept float or double */
  doppler->JDref = -1.0; InfoReal.dbl = doppler->JDref; type = OBIT_double;
  ObitInfoListGetTest(inUV->info, "JDref",    &type, dim, &InfoReal);
  if (type==OBIT_double)     doppler->JDref = InfoReal.dbl;
  else if (type==OBIT_float) doppler->JDref = (odouble)InfoReal.flt;

  doppler->RestFreq = -999.9; InfoReal.dbl = doppler->RestFreq; type = OBIT_double;
  ObitInfoListGetTest(inUV->info, "RestFreq", &type, dim, &InfoReal);
  if (type==OBIT_double)     doppler->RestFreq = InfoReal.dbl;
  else if (type==OBIT_float) doppler->RestFreq = (odouble)InfoReal.flt;

  doppler->VelLSR = -9999.0; InfoReal.flt = doppler->VelLSR; type = OBIT_float;
  ObitInfoListGetTest(inUV->info, "VelLSR",   &type, dim, &InfoReal);
  if (type==OBIT_float)       doppler->VelLSR = InfoReal.flt;
  else if (type==OBIT_double) doppler->VelLSR = (ofloat)InfoReal.dbl;
  doppler->VelLSR *= 1.0e-3;   /* to km/sec */

  doppler->refChan  = -9999.0; InfoReal.flt = doppler->refChan; type = OBIT_float;
  ObitInfoListGetTest(inUV->info, "refChan",  &type, dim, &InfoReal);
  if (type==OBIT_float)       doppler->refChan = InfoReal.flt;
  else if (type==OBIT_double) doppler->refChan = (ofloat)InfoReal.dbl;

  /* Defaults */
  if (doppler->JDref<0.0)       doppler->JDref = inDesc->JDObs;  /* Observation date */
  if (doppler->RestFreq<0.0)    doppler->RestFreq = inDesc->restFreq;
  if (doppler->RestFreq<0.0)    doppler->RestFreq = inDesc->freq;
  if (doppler->VelLSR<-9998.0)  doppler->VelLSR = 0.0;
  if (doppler->refChan<-9998.0) doppler->refChan = defaultRefChan(doppler, err);
  if (err->error) goto cleanup;
  if (doppler->refChan<-9998.0) doppler->refChan = inDesc->altCrpix;
  if (doppler->refChan==0.0)    doppler->refChan = 1+inDesc->inaxes[inDesc->jlocf]/2;
  /* Save */
  dim[0] = dim[1] = dim[2] = dim[3] = dim[4] = 1;
  ObitInfoListAlwaysPut(outUV->info, "refChan", OBIT_float, dim, &doppler->refChan);

  /* Tell reference channel */
  Obit_log_error(err, OBIT_InfoErr, 
		 "Using reference channel %f", doppler->refChan);
  ObitErrLog(err); /* Show messages */


  /* Update output header with velocity info */
  inDesc->restFreq  = doppler->RestFreq;
  inDesc->altRef    = doppler->VelLSR*1000.0;
  inDesc->altCrpix  = doppler->refChan;
  outDesc->restFreq = doppler->RestFreq;
  outDesc->altRef   = doppler->VelLSR*1000.0;
  outDesc->altCrpix = doppler->refChan;

  /* Create work vis vuffer */
  tmpVis = g_malloc0(inDesc->lrec*sizeof(ofloat));

  /* Loop over visibilities */
  done   = FALSE;
  gotOne = FALSE;

  /* Next buffer load */
  while ((iretCode==OBIT_IO_OK) && (oretCode==OBIT_IO_OK)) {
    if ((!gotOne) || (inUV->myDesc->numVisBuff<=0)) { /* need to read new record? */
      if (doCalSelect) iretCode = ObitUVReadSelect (inUV, inUV->buffer, err);
      else iretCode = ObitUVRead (inUV, inUV->buffer, err);
    }
    if (err->error) goto cleanup;

    /* Are we there yet??? */
    done = (inDesc->firstVis >= inDesc->nvis) || (iretCode==OBIT_IO_EOF);
    if (done) break; 

    /* Make sure valid data found */
    if (inUV->myDesc->numVisBuff<=0) continue;
    
    /* How many */
    outDesc->numVisBuff = inDesc->numVisBuff;
    
    /* loop over visibilities in buffer */
    for (ivis=0; ivis<inDesc->numVisBuff; ivis++) { 
      /* Copy random parameters */
      iindx = ivis*inDesc->lrec;
      oindx = ivis*outDesc->lrec;

	/* Are we there yet??? */
	done = (inDesc->firstVis >= inDesc->nvis) || 
	  (iretCode==OBIT_IO_EOF);
	if (done) goto done;

	/* Copy random parameters */
	for (j=0; j<inDesc->nrparm; j++) outBuffer[oindx+j] = inBuffer[iindx+j];

	/* Now process this visibility */
	cbase = inBuffer[iindx+inDesc->ilocb]; /* Baseline */
	ant1 = (cbase / 256.0) + 0.001;
	ant2 = (cbase - ant1 * 256) + 0.001;
	lastSubA = (olong)(100.0 * (cbase -  ant1 * 256 - ant2) + 0.5);
	/* source ID */
	if (inDesc->ilocsu>=0) curSourceID = inBuffer[iindx+inDesc->ilocsu];
	/* Loop over data in visibility, copy to Spectrum to shift 
	 then loop over components copy to outbuffer */
	/* loop over IF */
	for (iif=0; iif<nif; iif++) {
	  ifoff  = iif * incif;
	  oifoff = iif * oincif;
	  /* Loop over polarization */
	  for (istok=0; istok<nstok; istok++) {
	    soff  = ifoff +  istok * incs;
	    osoff = oifoff + istok * oincs;
	    /* Zero fill */
	    for (ichan=0; ichan<nchanp; ichan++) {
	      Spectrum[ichan*2] = 0.0; Spectrum[ichan*2+1] = 0.0;
	    }
	    /* Loop over frequency channel */
	    tndx = 0;
	    for (ichan=0; ichan<nchan; ichan++) {
	      foff = iindx + inUV->myDesc->nrparm + soff + ichan * incf;
	      /* to Spectrum */
	      if (inBuffer[foff+2]>0.0) {
		Spectrum[tndx++] = inBuffer[foff+0];
		Spectrum[tndx++] = inBuffer[foff+1];
	      } else tndx += 2;
	      foff += 3;
	    } /* end frequency loop */
	    /* Shift -  How many channels to shift */
	    shift    = CalcShift (doppler, ant1, ant2, curSourceID,
				  inBuffer[iindx+inDesc->iloct], err);
	    /* Approximate correction in channel number */
	    if (shift>0) off = (olong)(shift+0.5);
	    else         off = (olong)(shift-0.5);
	    /* shift = 1.0;  DEBUG */
	    if (err->error) goto cleanup;

	    /* uvw scaling to first IF, poln */
	    if ((iif==0) && (istok==0)) {
	      uvwScale = (outDesc->freq - outDesc->chIncIF[0]*shift) / 
		outDesc->freq;
	    }

	    sideband = 1;  /* Sideband - always effectively upper */
	    if (ant1==ant2) {
	      /* autocorrelation */
	      ACShift (doppler->Spectrum, doppler->Work,
		       doppler->FFTFor, doppler->FFTRev, 
		       sideband, nchanp, shift, FALSE); 
	    } else {
	      /* crosscorrelation */
	      CCShift (doppler->Spectrum, doppler->Work,
		       doppler->FFTFor, doppler->FFTRev, 
		       sideband, nchanp, shift); 
	    }
	    /* Copy to output buffer */
	    /* Loop over frequency channel */
	    tndx = 0;
	    for (ichan=0; ichan<nchan; ichan++) {
	      foff  = iindx + inUV->myDesc->nrparm  + osoff + ichan * incf;
	      ofoff = oindx + outUV->myDesc->nrparm + osoff + ichan * oincf;
	      toff = foff + (off*incf) + 2;  /* roughly align new and old */
	      if ((inBuffer[toff]>0.0) && (ichan>off) && ((ichan+off)<nchan)) {
		outBuffer[ofoff+0] = Spectrum[tndx++];
		outBuffer[ofoff+1] = Spectrum[tndx++];
		outBuffer[ofoff+2] = inBuffer[toff];  /* more or less */
		/* end valid data */
	      } else {
		tndx += 2;
		outBuffer[ofoff+0] = 0.0;
		outBuffer[ofoff+1] = 0.0;
		outBuffer[ofoff+2] = 0.0;
	      }
	      ofoff += 3;
	      foff  += 3;
	    }/* end frequency loop */
	  } /* end poln loop */
	} /* end IF loop */
    } /* end loop processing buffer of input data */

    /* scale uvw */
    outBuffer[oindx+inDesc->ilocu] *= uvwScale;
    outBuffer[oindx+inDesc->ilocv] *= uvwScale;
    outBuffer[oindx+inDesc->ilocw] *= uvwScale;
       
    /* Write it out */
    oretCode = ObitUVWrite (outUV, outUV->buffer, err);
    if (err->error) goto cleanup;
  } /* End loop over input file */
  
  /* End of processing */
 done:
  /* check for errors */
  if ((iretCode > OBIT_IO_EOF) || (oretCode > OBIT_IO_EOF) ||
      (err->error)) goto cleanup;
  
  /* Restore no vis per read in output */
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut (outUV->info, "nVisPIO", OBIT_long, dim, &oldNPIO);

  /* Update SU table for multisource */
  UpdateSUTable (doppler, outUV, err);
  if (err->error) goto cleanup;

  /* Cleanup */
 cleanup:
  if (tmpVis) g_free(tmpVis); tmpVis = NULL;
  /* close files */
  iretCode = ObitUVClose (inUV, err);
  oretCode = ObitUVClose (outUV, err);
  if ((oretCode!=OBIT_IO_OK) || (iretCode!=OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, outUV->name, outUV);
  
  return outUV;
} /* end ObitDopplerCVel */

/**
 * Convert a Julian date to year, doy and time
 * Apdapted from ASM Algorithm no. 199
 * \param JD Julian date.
 * \param year  year
 * \param doy   Day of year
 * \param ut    Time in days since doy on year of JD
 */
void ObitDopplerJD2Date (odouble JD, olong *year, olong *doy, ofloat *ut)
{
  olong  id, im, iy, ic, it;
  ollong tt;
  odouble j, y, d, m;
  /*              1  2  3  4  5  6  7  8  9 10 11 12 */
  /*olong *months = [31,28,31,30,31,30,31,31,30,31,30,31];*/
  olong months[] = {0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334};

  /* error check */
  if (JD<1.0) {
    *year = -1;
    *doy  = -1;
    *ut   = -1.0;
    return;
  }
  
  j = (olong) (JD + 0.50 - 1721119.0);
  y = (olong) ((4.0*j - 1.00) / 146097.0);
  ic = y + 0.00010;
  j = 4.0*j - 1.00 - 146097.0*y;
  d = (olong) (j * 0.250);
  j = (olong) ((4.00*d + 3.00) / 1461.00);
  d = 4.00*d + 3.00 - 1461.00*j;
  d = (olong) ((d+4.00) * 0.250);
  m = (olong) ((5.00*d - 3.00) / 153.00);
  id = 5.00*d - 3.00 - 153.00*m;
  id = (id + 5) / 5;
  iy = j + 100*ic;
  im = m;
  if (im < 10) {
    im = im + 3;
  } else {
    im = im - 9;
  iy = iy + 1;
  }

  /* get day of year */
  *year = iy;
  it    = months[im-1] + id;
  /* Yeap year? */
  if ((((*year)%4)==0) && (im>2)) it++;
  *doy = it;             /* Day of year */
  tt = (ollong)(JD-0.5);
  *ut  = JD - tt;        /* time */

} /* end ObitDopplerJD2Date */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitDopplerClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitDopplerClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitDopplerClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitDopplerClassInfoDefFn (gpointer inClass)
{
  ObitDopplerClassInfo *theClass = (ObitDopplerClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitDopplerClassInit;
  theClass->newObit       = (newObitFP)newObitDoppler;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitDopplerClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitDopplerGetClass;
  theClass->ObitCopy      = (ObitCopyFP)ObitDopplerCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitDopplerClear;
  theClass->ObitInit      = (ObitInitFP)ObitDopplerInit;
  theClass->ObitDopplerCreate = (ObitDopplerCreateFP)ObitDopplerCreate;

} /* end ObitDopplerClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitDopplerInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitDoppler *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->thread   = newObitThread();
  in->info     = newObitInfoList(); 
  in->uvdata   = NULL;
  in->Spectrum = NULL;
  in->Work     = NULL;
  in->FFTFor   = NULL;
  in->FFTRev   = NULL;
  in->antList  = NULL;
  in->sourceList  = NULL;
  in->source   = NULL;
  in->year     = -999;
  in->doy      = -999;
  in->JDref    = -9999.0;
  in->RestFreq = -9999.0;
  in->VelLSR   = -9999.0;
  in->refChan  = -9999.0;
} /* end ObitDopplerInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitDoppler* cast to an Obit*.
 */
void ObitDopplerClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitDoppler *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  in->thread   = ObitThreadUnref(in->thread);
  in->info     = ObitInfoListUnref(in->info); 
  in->uvdata   = ObitUVUnref(in->uvdata);
  in->Spectrum = ObitCArrayUnref(in->Spectrum);
  in->Work     = ObitCArrayUnref(in->Work);
  in->FFTFor   = ObitFFTUnref(in->FFTFor);
  in->FFTRev   = ObitFFTUnref(in->FFTRev);
  in->antList  = ObitAntennaListUnref(in->antList);
  in->sourceList  = ObitSourceListUnref(in->sourceList);
  in->source   = ObitSourceUnref(in->source);
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitDopplerClear */

/**
 * Calculates the velocity component of the observer  
 * with respect to the local standard of rest as projected onto  
 * a line specified by the right ascension and declination (epoch=  
 * date) for a specified time.  
 * The location of the observer is specified by the geocentric  
 * coordinates X, Y and Z (in metres)  
 *  
 * This version takes into account components of the observer's  
 * motion due to the rotation of the earth, the revolution of the  
 * earth-moon barycenter about the sun, and the motion of the earth's  
 * centre about the earth-moon barycenter.  The perturbations of the  
 * earth's orbit due to planets are neglected.  The absolute precision  
 * of this version of dop is about 0.004 km/sec, but since the dominant  
 * error term is slowly varying the relative error will be considerably  
 * less for times up to a week or so.  
 *  
 * The routine omits the planetary perturbations on the earth's orbit.  
 * They amount to about 0.003 km/sec and are thought to be the largest  
 * contribution to the error in the velocity.  
 * Routine translated from the AIPSish DOPLR.FOR/DOPLR  
 * \param ra    Right Ascension of date (radians) 
 * \param dec   Declination of date (radians) 
 * \param year  Year of observation 
 * \param doy   Day number of observation 
 * \param ut    UT of datum (days) 
 * \param x     X Geocentric station coord (m)
 * \param y     Y Geocentric station coord (m)  left hand?
 * \param z     Z Geocentric station coord (m)
 * \param vsun  [out] The component of the sun's motion 
 *              with respect to the LSR as projected onto the line of
 *              sight to the source. (m/s) 
 * \return      The total velocity component of the observer's
 *              motion  with respect to the LSR as projected 
 *              onto the line of sight to the source. (m/s) 
 */
odouble VelLSR (odouble ra, odouble dec, olong year, olong doy,
		odouble ut, odouble x, odouble y, odouble z,
		odouble* vsun) 
{
  odouble   zro, aaa, dd, dela, deldd, xo, yo, zo, lng,
    vt, ra1, pie, cc, cs, s, cat, wlong, du, tu, utda, smd, t,
    start, c1, gst, xlst, rho, vrho, dlat, vobs, am, e, ai, vs, xlam,
    alam, aa, an, hop, v, omga, omgar, amon, gamp, pim, em, olamm,
    aim, amm, vsm, alamm, anm, aam, hopm, vmon, beta, dc, along, algm,
    betam;

  /* Vel of Sun wrt LSR is 20 km/s towards RA 18 hrs Dec +30 deg */
  zro = 0.0;
  aaa = 18.0 * G_PI / 12.0;
  dd  = 30.0 * G_PI / 180.0;

  /* precesses this direction  to date */
  rdmove( 1900, year, 1, doy, aaa, dd, &dela, &deldd, &dc);
  aaa += dela;
  dd  += deldd;

  /* this velocity is converted to cartesian components */
  xo = 20.0 * cos (aaa) * cos (dd);
  yo = 20.0 * sin (aaa) * cos (dd);
  zo = 20.0 * sin (dd);

  /* ra1 = RA in days */
  ra1 = (ra * (12.0 / G_PI) ) / 24.0;

  /* CC, CS, and S are the direction cosines corresponding to ra and dec */
  cc = cos (dec) * cos (ra);
  cs = cos (dec) * sin (ra);
  s  = sin (dec);

  /* vsun is the projection onto the line of sight to the  star of the sun's motion wrt 
     to the LSR (km/sec) */
  (*vsun) = -xo * cc - yo * cs - zo * s;

  /* determine geodetic lat & long. */
  if ((x == 0.0)  &&  (y == 0.0)  &&  (z == 0.0)) {
    cat = 0.0;
    lng = 0.0;
  } else {
    GeoDll (x, y, z, &cat, &lng);
  } 
  wlong = lng / (2.0 * G_PI);

  /* Time calculations
     The epoch is 1900 jan 0.5 UT  = julian day 2415020.0
     DU is the time from epoch to jan 0.0 of the current year (days) */
  du = (julDay(year) - 2415020) - 0.5;

  /* tu is du converted to julian  centuries */
  tu = du / 36525.0;

  /* utda is the gmt from  jan 0.0 to the present (days) */
  utda = doy + ut;

  /* smd is the time from the epoch to the present (days) */
  smd = du + utda;

  /* t is smd converted to julian centuries */
  t = smd / 36525.0;

  /* start is the greenwich mean sidereal time on jan 0.0 (days) */
  start = (6.0 + 38.0 / 60.0 + (45.836 + 129.1794 +
				8640184.542 * (tu-0.7) +
				0.0929 * tu*tu) / 3600.0) / 24.0; 
  /* c1 is the conversion factor from  solar time to sidereal time */
  c1 = 0.997269566414;

  /* gst is the greenwich mean sidereal time (days) */
  gst = start + utda / c1;

  /* xlst is the local mean sidereal  time (from jan 0) (days) */
  xlst = gst - wlong;
  xlst -= (olong)xlst;

  /* Observer's motion wrt earth's centre 
     rho is the radius vector from the earth's center to the observer */
  rho = sqrt (x*x + y*y + z*z);

  /* vrho is corresponding circular velocity (meters/sidereal day) */
  vrho = 2.0e0 * G_PI * rho;

  /* converted to kilometers/sec */
  vrho =  vrho / 24.0e3 / 3600.0 / c1;

  /* reduction of geodetic latitude  to geocentric latitude (arcsecs) */
  dlat = -(11.0 * 60.0 + 32.7430) * sin (2.0 * cat) + 1.1633
    * sin (4.0 * cat) -0.0026e0 * sin (6.0 * cat);

  /* convert cat to geocentric lat  (radians) */
  cat += dlat * G_PI / 3600.0 / 180.0e0;

  /* vobs is the projection onto the line of sight to the star of the 
     velocity of the observer wrt the  earth's center (km/sec) */
  vobs = vrho * cos (cat) * cos (dec) * sin (2.0 * G_PI * (xlst - ra1));
  
  /* Earth's orbit about sun  am is the mean anomaly 
     (of the earth's orbit) (radians) */
  am = (358.47583 + 0.9856002670 * smd-0.000150 * t*t -
	0.000003 * t*t*t) * G_PI / 180.0; 
  
  /* pie is the mean longitude of  perihelion (radians) */
  pie = (101.22083 + 0.0000470684 * smd + 0.000453 * t*t +
	 0.000003 * t*t*t) * G_PI / 180.0;

  /* e is the eccentricity of the orbit (dimensionless) */
  e = 0.01675104 - 0.00004180 * t - 0.000000126 * t*t;
  
  /* ai is the mean obliquity of the ecliptic (radians) */
  ai = (23.452294 - 0.0130125 * t - 0.00000164 * t*t +
	0.000000503 * t*t*t) * G_PI / 180.0; 

  /* vs is the true anomaly (approximate formula) (radians) */
  vs = am + (2.0 * e - 0.25 * e*e*e) * sin (am) + 1.25 * e*e *
    sin (2.0 * am) + 13.0 / 12.0 * e*e*e * sin (3.0 * am);

  /* xlam is the true longitude of the earth as seen from the sun (radians) */
  xlam = pie + vs;

  /* alam is the true longitude of the sun as seen from the earth (radians) */
  alam = xlam + G_PI;

  /* beta is the latitude of the star (radians). 
     along is the longitude of the star (radians) */
  coordTran (zro, zro, (-G_PI/2.0), (G_PI/2.0 - ai), ra, dec,
	     &along, &beta);

  /* aa is the semi-major axis of  the earth's orbit (km) */
  aa = 149598500.0;

  /* an is the mean angular rate of the earth about the sun (radians/day) */
  an = 2.0 * G_PI / 365.2564;

  /* hop is h/p from smart = the component of the earth's velocity 
     perpendicular to the radius vector (km/day) */
  hop = an * aa / sqrt (1.0 - e*e);

  /* converted to km/sec */
  hop /= 86400.0;

  /* v is the projection onto the line of sight to the star of  the
     velocity of the earth-moon barycenter with respect to the sun (km/sec) */
  v = -hop * cos (beta) * (sin(alam - along) - e * sin (pie - along));

  /* Calculate moon's orbit around the earth-moon barycentre,. 
     omga is the longitude of the mean ascending node of the lunar pole (degrees) */
  omga = 259.183275 - 0.0529539222 * smd + 0.002078 * t*t + 0.000002 * t*t*t;

  /* omgar is omga in radians */
  omgar = omga * G_PI / 180.0;

  /* amon is omga plus the mean lunar longitude of the moon 
     (degrees - should be 13.1763965268) */
  amon = 270.434164 + 13.176396527 * smd - 0.001133 * t*t +
    0.0000019 * t*t*t;

  /* gamp (gamma-prime) is omga plus the lunar longitude of lunar perigee */
  gamp = 334.329556 + 0.1114040803 * smd - 0.010325 * t*t -
    0.000012 * t*t*t;

  /* pim is the mean lunar longitude of lunar perigee (in radians) */
  pim = (gamp - omga) * G_PI / 180.0;

  /* em is the eccentricity of the  lunar orbit */
  em = 0.054900489;

  /* olamm is the mean lunar long. of the moon (in radians) */
  olamm = (amon - omga) * G_PI / 180.0;

  /* aim is the inclination of the lunar orbit to the ecliptic (radians) */
  aim = 5.1453964 * G_PI / 180.0;

  /* amm is the approximate mean anomaly (radians) */
  amm = olamm - pim;

/* vsm is the true anomaly (approx formula) (radians) */
  vsm = amm + (2.0 * em - 0.25 * em*em*em) * sin (amm) + 1.25 *
    em*em * sin (2.0 *amm) + 13.0 / 12.0 * em*em*em * sin (3.0 * amm);  

  /* alamm is the true lunar longitude of the moon (radians) */
  alamm = pim + vsm;

  /* anm is the mean angular rate of the lunar rotation (radians/day) */
  anm = 2.0 * G_PI / 27.321661;

  /* aam is the semi-major axis of the lunar obrit (kilometers) */
  aam = 60.2665 * 6378.388;

  /* betam is the lunar latitude of the star (rads). 
     algm is the  lunar longitude of the star. */
  coordTran (omgar, zro, (omgar - G_PI / 2.0), (G_PI / 2.0 - aim),
	     along, beta, &algm, &betam);
  
  /* hopm is h/p from smart = the component of the lunar velocity
     perpendicular to the radius vector (km/day) */
  hopm = anm * aam / sqrt (1.0 - em*em);

  /* converted to km/sec */
  hopm /=  86400.0;

/* vmon is the projection onto the line of sight to the star of the
   velocity of the earth's center with respect to the earth-moon 
   barycenter (km/s) 
   (the 81.30 is the ratio of the earth's mass to the moon's mass) */
  vmon = -hopm / 81.30 * cos (betam) * (sin(alamm - algm) - em * sin (pim - algm));

  /* vt is observer's velocity in direction of ra, dec */
  vt = v + (*vsun) + vmon + vobs;
  vt *= 1000.0;      /* to m/s */
  (*vsun) *= 1000.0; /* to m/s */
  return vt;
} /* end of routine VelLSR */ 

/**
 * Routine to convert the longitude-like (A1) and latitude-like (B1)  
 * coordinates of a point on a sphere into the corresponding  
 * coordinates (A2,B2) in a different coordinate system that is  
 * specified by the coordinates of its origin (AO, BO) and its north  
 * pole (AP, BP) in the original coordinate system.  
 * The range of A2 will be from -pi to pi.  
 * Routine translated from the AIPSish DOPLR.FOR/CDTRN  
 * \param ao   longitude of origin (rad)
 * \param bo   latitude of origin (rad)
 * \param ap   longitude of pole (rad)
 * \param bp   latitude of pole (rad)
 * \param a1   input longitude (rad)
 * \param b1   input latitude (rad)
 * \param a2   [out] output longitude (rad)
 * \param b2   [out] output latitude (rad)
 */
static void coordTran (odouble ao, odouble bo, odouble ap, odouble bp,
		odouble a1, odouble b1, odouble* a2, odouble* b2)  
{
  odouble sbo, cbo, sbp, cbp, sb1, cb1, sb2, cb2, saa, caa, sbb, cbb,
    sa2, ca2, ta2o2;  
  
  sbo = sin (bo);
  cbo = cos (bo);
  sbp = sin (bp);
  cbp = cos (bp);
  sb1 = sin (b1);
  cb1 = cos (b1);
  
  sb2 = sbp * sb1 + cbp * cb1 * cos (ap - a1);
  (*b2) = asin(sb2);
  cb2 = cos ((*b2));
  
  saa = sin (ap - a1) * cb1 / cb2;
  caa = (sb1 - sb2 * sbp) / (cb2 * cbp);

  cbb = sbo / cbp;
  sbb = sin (ap - ao) * cbo;

  sa2 = saa * cbb - caa * sbb;
  ca2 = caa * cbb + saa * sbb;
  if (ca2<=0.0) ta2o2 = (1.0 - ca2) / sa2;
  else          ta2o2 = sa2 / (1.0 + ca2);
  (*a2) = 2.0 * atan(ta2o2);

return;
} /* end of routine coordTran */ 

/**
 * Calculates the geodetic latitude and logitude of a
 * point on the earth's surface from its geocentric coordinates.
 * Routine translated from the AIPSish GEODLL 
 * \param x     X (m)
 * \param y     Y (m)
 * \param z     Z (m)
 * \param lat   [out] Geodetic latitude (rads)
 * \param lng   [out] Geodetic longitude (rads)
 */
static void GeoDll (odouble x, odouble y, odouble z, 
	     odouble *lat, odouble  *lng)
{
  odouble  gclat, dlat, c1=-3.358513e-3, c2=5.6398e-6, c3=-0.01261e-6;

  *lng = atan2 (y, x);
  gclat = asin ( z / sqrt (x*x + y*y + z*z) );
  dlat =  c1 * sin ( 2.0 * gclat )
        + c2 * sin ( 4.0 * gclat )
        + c3 * sin ( 6.0 * gclat );
  *lat = gclat - dlat;
} /* end GeoDll */

/**
 * Calculates the corrections in RA and Dec  
 * to be added to mean coordinates for epoch yeari.  
 * If the day-number is known, use it for nda and set mo=1  
 * Move also calculates the equation of the equinoxes (dc)  
 * which may be added to the mean sidereal time to give the  
 * apparent sidereal time.  
 * delr and deld contain corrections for precession,  
 * aberration, and some terms of nutation.  
 * Routine translated from the AIPSish DOPLR.FOR/RDMOVE  
 * \param yeari    Epoch of coords 
 * \param yearf    Year of observation 
 * \param mo       Month of observation 
 * \param nda      Day of observation 
 * \param ra       RA of epoch (rads) 
 * \param d        Dec of epoch (rads) 
 * \param delr     [out] RA correction (rads) 
 * \param deld     [out] Dec correction (rads) 
 * \param dc       [out] Equation of the equinoxes (mins of time) 
 */
static void rdmove (olong  yeari, olong yearf, olong mo, olong nda, 
	     odouble ra, odouble d, 
	     odouble* delr, odouble* deld, odouble* dc ) 
{
  odouble    snd, csd, tnd,csr, snr, al, to, t, zetao, z, theta, am,
    an, alam, snl, csl, omega, arg, dlong, doblq;  
  
  snd = sin (d);
  csd = cos (d);
  tnd = snd / csd;
  
  csr = cos (ra);
  snr = sin (ra);

  /* al is an approximate day number  (i.e. the number of days since 
     January of the year yearf). */
  al = 30 * (mo - 1) + nda;

  /* to is the time from 1900 to ayri (centuries), t is the time  
     from ayri to date (yearf,mo,nda) (centuries) 
     (365.2421988 is the number of ephemeris days in a tropical year) */
  to = (yeari - 1900) / 100.0;
  t = ((yearf - yeari) + al / 365.2421988) / 100.0;

  /* zetao is a precessional angle */
  zetao = (2304.250 + 1.396 * to) * t + 0.302 * t*t + 0.018 * t*t*t;

  /* ditto for z */
  z = zetao + 0.791 * t*t;

  /* and theta */
  theta = (2004.682 - 0.853 * to) * t - 0.426 * t*t - 0.042 * t*t*t;

  /* am and an are the m and n  precessional numbers */
  am = (zetao + z) * 4.848136811e-6;
  an = theta * 4.848136811e-6;

  /* alam is an approximate mean longitude for sun (radians) */
  alam = (0.985647 * al + 278.5) * 0.0174532925;
  snl = sin (alam);
  csl = cos (alam);

  /* delr is the annual aberration  term in RA (radians) */
  (*delr) = -9.92413605e-5 * (snl * snr + 0.91745051 * csl * csr) / csd;

  /* plus precession terms + am + an * snr * tnd;  
     deld is ditto above in dec */
  (*deld) = -9.92413605e-5 * (snl * csr * snd - 0.91745051 * csl
			      * snr * snd + 0.39784993 * csl * csd) +
    
    an * csr;  
  
  /* calculate the nutation (approx) omega is the angle of the 
     first term of nutation */
  omega = 259.183275 - 1934.142 * ( to + t );

  /* arg is omega converted to rads */
  arg = omega * 0.0174532925;

  /* dlong is the nutation in longitude (delta-psi) (radians) */
  dlong = -8.3597e-5 * sin (arg);

  /* doblq is the nutation in obliquity (delta-epsilon) (rads) */
  doblq = 4.4678e-5 * cos (arg);

  /* add nutation in ra into delr */
  (*delr) += dlong * ( 0.91745051 + 0.39784993 * snr * tnd ) - 
    csr * tnd * doblq; 

  /* and dec. */
  (*deld) += 0.39784993 * csr * dlong + snr * doblq;
  
  /* dc is the equation of the  equinoxes (minutes of time) */
  (*dc) = dlong * 210.264169;
  return;
} /* end of routine rdmove */ 

/**
 *   This function computes the Julian day number at 12 hrs ut on 
 *   January 0 of the year year (Gregorian calendar). Julda is an 
 *   integer because of this definition.  For example, 
 *   julda = 2439856 for year = 1968. 
 * Translated from the AIPSIsh JULDA 
 * \param  year  Year
 * \return Julian date of begining of year
 */
static olong julDay (year) {
  olong    julda, yearm1, ic;
  
  yearm1 = year-1;
  ic = yearm1/100;
  julda = 1721425 + 365 * yearm1 + yearm1 / 4 - ic + ic / 4;
  return julda;
} /* end of routine julDay */ 

/**
 * Frequency shift the polarization spectra for an antenna/IF
 * Use Fourier shift theorem
 * Adopted from AIPS TPSHFT.FOR
 * \param Spectrum Complex spectrum
 * \param Work     Work CArray the size of Spectrum
 * \param FFTFor   FFT object for forward transform
 * \param FFTRev   FFT object for reverse transform
 * \param sideband which sideband (1=upper, -1 = lower)
 * \param nchan    number of channels in Spectrum
 * \param shift    Number of channels to shift
 * \param doSmo    0=> no smoothing, 1=>Hanning.
 */
static void ACShift (ObitCArray  *Spectrum,  ObitCArray  *Work,  
		     ObitFFT *FFTFor, ObitFFT *FFTRev, 
		     olong sideband, olong nchan, ofloat shift, olong doSmo)
{
  ofloat dela, rfact, arg, xre, xim, norm, *spectrum, *work;
  olong   nfrq, nfrq2, i, n2, jf, jbin, ntrans, fftdir;
  olong pos[2] = {0, 0};
  
  nfrq = nchan;
  rfact = 1.0;
  spectrum = ObitCArrayIndex (Spectrum, pos);
  work     = ObitCArrayIndex (Work, pos);

  /* reflect spectrum */
  nfrq2  = nfrq * 2;
  ntrans = nfrq;
  fftdir = -1;
  ObitFFTC2C (FFTFor, Spectrum, Work);

  /* determine shift parms */
  dela  = -2.0*G_PI * shift / nfrq;

  /* shift AC spectrum. */
  n2 = nfrq / 2;
  for (i= 0; i< nfrq; i++) { /* loop 200 */
    if (i+1 <= n2) {
      jf = i;
    } else {
      jf = -nfrq + i + 1;
    } 
    jbin = jf + n2;

    /* Hanning? */
    if (doSmo == 1) rfact = 0.5*(1.0-cos(2.0*G_PI*jbin/(nfrq-1)));

    /* Multiple by phase ramp to shift */
    arg = dela * jf;
    xre = work[2+i];
    xim = work[2*i+1];
    work[2*i]   = rfact * (cos (arg) * xre - sin (arg) * xim);
    work[2*i+1] = rfact * (sin (arg) * xre + cos (arg) * xim);
  } /* end loop   L200 */;

  /* transform back to spectrum */
  fftdir = -fftdir;
  ObitFFTC2C (FFTRev, Work, Spectrum);

  /* normalize */
  norm = 1.0 / (ofloat)nfrq;
  for (i= 0; i< nfrq2; i++) { /* loop 300 */
    spectrum[i] *= norm;
  } /* end loop   L300 */;

  /* form real only */
  for (i= 0; i< nfrq; i++) { /* loop 400 */
    spectrum[2*i]   = sqrt(spectrum[2*i]*spectrum[2*i] + spectrum[2*i+1]*spectrum[2*i+1]);
    spectrum[2*i+1] = 0.0;
  } /* end loop   L400 */;

} /* end ACShift */

/**
 * Frequency shift the spectra for a baseline/IF
 * Adopted from AIPS XCSHNQ.FOR, NQTOCF.FOR, NQTOCS.FOR
 * Use Fourier shift theorem
 * \param Spectrum Complex spectrum
 * \param Work     Work CArray the size of Spectrum
 * \param FFTFor   FFT object for forward transform
 * \param FFTRev   FFT object for reverse transform
 * \param sideband which sideband (1=upper, -1 = lower)
 * \param nchan    number of channels in Spectrum
 * \param shift    Number of channels to shift
 */
static void CCShift (ObitCArray *Spectrum,  ObitCArray *Work,
		     ObitFFT *FFTFor, ObitFFT *FFTRev, 
		     olong sideband, olong nchan, ofloat shift)
{
  ofloat  del, del1, cd, sd, c, s, store, norm, temp1, temp2, *spectrum, *work;
  olong   nfrq, nxcf, kstart, kstop, k, kk, ll, fftdir, i;
  olong pos[2] = {0, 0};

  nfrq = nchan;
  
  /* fft to double sideband xc-function */
  nxcf = nfrq * 2;
  spectrum = ObitCArrayIndex (Spectrum, pos);
  work     = ObitCArrayIndex (Work, pos);
  
  /* fill lower sideband array slots with zeroes */
  kstart = nfrq + 1;
  kstop  = nxcf;
  for (k=kstart-1; k<kstop; k++) { /* loop 10 */
    spectrum[2*k]   = 0.0;
    spectrum[2*k+1] = 0.0;
  }

  /* transform to xcf */
  fftdir = -sideband;
  /* call fourg (Spectrum, nxcf, fftdir, work); */
  /* Direction depends on sideband */
  if (sideband>0) 
    ObitFFTC2C (FFTFor, Spectrum, Work);
  else
    ObitFFTC2C (FFTRev, Spectrum, Work);
 
  /* flip Spectrum around to center correlation function  in first half of array */
  kstop = nfrq;
  norm = 1.0 / (ofloat)nxcf;
  for (k=0; k< kstop; k++) { /* loop 20 */
    kk = nxcf - kstop + k;
    ll = nxcf - kstop + k;
    temp1 = work[2*k];
    temp2 = work[2*k+1];
    work[2*k]    = work[2*kk] * norm;
    work[2*k+1]  = work[2*kk+1] * norm;
    work[2*ll]   = temp1 * norm;
    work[2*ll+1] = temp2 * norm;
  }
  
  /* minus sign in following equation for "del" is because increasing channel # 
     corresponds to decreasing lag  values. */
  del  = -1.0 * sideband * 2.0 * G_PI * shift / nxcf;
  del1 = -0.5 *  nxcf  * del;
  cd   = cos( del );
  sd   = sin( del );
  c    = cos( del1 );
  s    = sin( del1 );
  
  /* shift xc spectrum.  at this point the array data  should contain a 
     double sideband (2 x nchan complex) correlation function with the zero delay  
     in the nfrq+1'th channel. Multiply by phase ramp. */
  for (i= 0; i< nxcf; i++) { /* loop 50 */
    store       = work[i*2];
    work[i*2]   = work[i*2]*c - work[i*2+1]*s;
    work[i*2+1] =     store*s + work[i*2+1]*c;
    store = c;
    c     = c*cd - s*sd;
    s     = store*sd + s*cd;
  } /* end loop    L50 */;
  
  /* Rearrange the correlation function so that the center channel (nxcf/2+1) winds up 
     in the first element of the array to be  transformed.  
     The complex correlation function to be transformed should be (nfrq*2) points long. */
  for (i= 0; i<nfrq; i++) { 
    k = i + nfrq;
    temp1       = work[2*i];
    temp2       = work[2*i+1];
    work[2*i]   = work[2*k];
    work[2*i+1] = work[2*k+1];
    work[2*k]   = temp1;
    work[2*k+1] = temp2;
  }

  /* fft back to spectrum */
  /* fftdir determines which sideband will end  up in first half of Spectrum */
  fftdir = sideband;
  /* Direction depends on sideband */
  if (sideband>0) 
    ObitFFTC2C (FFTRev, Work, Spectrum);
  else
    ObitFFTC2C (FFTFor, Work, Spectrum);
} /* end CCShift */

/**
 * Calculate number of channels to shift
 * \param doppler  ObitDoppler object
 * \param ant1     1st antenna number
 * \param ant2     2nd  antenna number
 * \param souId    Source ID, 0=single source
 * \param time     Time in days since 0h on reference
 * \param err      Error stack, returns if not empty.
 * \return  Number of channels to shift
 */
static ofloat CalcShift (ObitDoppler *doppler, olong ant1, olong ant2, 
			 olong souId, ofloat time, ObitErr *err)
{
  ofloat shift=0.0;
  ofloat relTime;
  odouble appFreq1, appFreq2, appFreq;
  ObitUVDesc *uvDesc = doppler->uvdata->myDesc;
  gboolean newSource=FALSE;
  olong isou, i;
  gchar *source, timeStr[20];
  gchar *routine = "CalcShift";

  /* Need to calculate year and doy of reference? */
  if (doppler->year<0) {
    ObitDopplerJD2Date (uvDesc->JDObs, &doppler->year, &doppler->doy, &relTime);
  }

  /* New Source? */
  if (souId!=doppler->source->SourID) {
    doppler->source->SourID = souId;
    newSource = TRUE;
    /* Single source? get from header */
    if (souId==0) {
      source = uvDesc->object;
      ObitUVGetRADec (doppler->uvdata, 
		      &doppler->source->RAMean, &doppler->source->DecMean, 
		      err);
      if (err->error) Obit_traceback_val (err, routine, doppler->name, shift);
      /* Compute apparent position */
      ObitPrecessUVJPrecessApp (uvDesc, doppler->source);
   } else { /* Multi source */
     /* Lookup in source list */
     isou = -1;
     for (i=0; i<doppler->sourceList->number; i++) {
       if (doppler->sourceList->SUlist[i]->SourID==souId) 
	 {isou = i; break;}
     }
     /* Found it? */
     Obit_retval_if_fail((isou>=0), err, shift,
			 "%s: Could not find source id %d", 
			 routine, souId);
     /* Get positions  */
     source = doppler->sourceList->SUlist[isou]->SourceName;
     doppler->source->RAMean  = doppler->sourceList->SUlist[isou]->RAMean;
     doppler->source->DecMean = doppler->sourceList->SUlist[isou]->DecMean;
     doppler->source->RAApp   = doppler->sourceList->SUlist[isou]->RAApp;
     doppler->source->DecApp  = doppler->sourceList->SUlist[isou]->DecApp;
     /* Update source list with RestFreq, LSRVel */
     for (i=0; i<doppler->sourceList->SUlist[isou]->numIF; i++) {
       doppler->sourceList->SUlist[isou]->LSRVel[i]   = doppler->VelLSR*1000.0;
       doppler->sourceList->SUlist[isou]->RestFreq[i] = doppler->RestFreq;
     }
    }
  } /* end new Source */

  /* Get time relative to reference Date */
  relTime = time;

  /* DEBUG
  appFreq1 = VelLSR (doppler->source->RAApp*DG2RAD, 
		     doppler->source->DecApp*DG2RAD,
		     doppler->year, doppler->doy, (odouble)relTime,	     
		     doppler->antList->ANlist[ant1-1]->AntXYZ[0],
		     doppler->antList->ANlist[ant1-1]->AntXYZ[1],
		     doppler->antList->ANlist[ant1-1]->AntXYZ[2], &appFreq2); */

  /* Get apparent frequency of line at given velocity at this time */
  /* Antenna 1 */
  appFreq1 = ObitDopplerFreqLSR (doppler->RestFreq, doppler->VelLSR, 
				 doppler->source->RAApp, 
				 doppler->source->DecApp,  
				 doppler->year, doppler->doy, (odouble)relTime, 
				 doppler->antList->ANlist[ant1-1]->AntXYZ[0],
				 doppler->antList->ANlist[ant1-1]->AntXYZ[1],
				 doppler->antList->ANlist[ant1-1]->AntXYZ[2]);
				 /* Antenna 2 */
  appFreq2 = ObitDopplerFreqLSR (doppler->RestFreq, doppler->VelLSR, 
				 doppler->source->RAApp, 
				 doppler->source->DecApp,  
				 doppler->year, doppler->doy, (odouble)relTime, 
				 doppler->antList->ANlist[ant2-1]->AntXYZ[0],
				 doppler->antList->ANlist[ant2-1]->AntXYZ[1],
				 doppler->antList->ANlist[ant2-1]->AntXYZ[2]);

  /* Average */
  appFreq = 0.5 * appFreq1 + 0.5 * appFreq2;

  /* How many channels? CHECK - not sure about this */
  shift = -(appFreq - (uvDesc->crval[uvDesc->jlocf] + 
		       (uvDesc->altCrpix-uvDesc->crpix[uvDesc->jlocf]) * 
		       uvDesc->cdelt[uvDesc->jlocf])) /
    uvDesc->cdelt[uvDesc->jlocf];

  /* Tell shift */
  if (newSource) {
    day2dhms (time, timeStr);
    Obit_log_error(err, OBIT_InfoErr, "%s %s, shift %f ch", 
		   source, timeStr, shift);
    ObitErrLog(err); /* Show messages */
  }

  /* Sanity check */
  Obit_retval_if_fail((shift<uvDesc->inaxes[uvDesc->jlocf]), err, shift,
		      "%s: Excessive shift %f ch, check RestFreq", 
		      routine, shift);
  return shift;
} /* end  CalcShift */


/**
 * For multisource inData update velocity information
 * \param in      The object with SourveList to use
 * \param inData  UV data with SU table to update
 * \param err     Obit error stack object.
 */
static void UpdateSUTable  (ObitDoppler *in, ObitUV *inData, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitTableSU *SUTab=NULL;
  ObitTableSURow *row=NULL;
  olong i, isou, numIF, iRow, iSUver;
  gchar *routine ="UpdateSUTable";

  /* error checks */
  if (err->error) return;

  /* Check if multisource */
  if (inData->myDesc->ilocsu<0) return;
  if ((in->sourceList==NULL) || (in->sourceList->SUlist[0]==NULL)) return;
  if (in->sourceList->number<=0) return;

  /* Get Table */
  iSUver = 1;
  numIF  = in->sourceList->SUlist[0]->numIF;
  SUTab  = 
    newObitTableSUValue (in->name, (ObitData*)inData, &iSUver, OBIT_IO_ReadWrite, numIF, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  if (SUTab==NULL) return;

  /* Open input table */
  retCode = ObitTableSUOpen (SUTab, OBIT_IO_ReadWrite, err);
  if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;

  /* Set row */
  row  = newObitTableSURow (SUTab);
  ObitTableSUSetRow (SUTab, row, err);
  if (err->error) goto cleanup;
 
 /* Loop over table copying selected data */
  for (iRow=1; iRow<=SUTab->myDesc->nrow; iRow++) {
    /* Read */
    retCode = ObitTableSUReadRow (SUTab, iRow, row, err);
    if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;

    /* Lookup in source list */
     isou = -1;
     for (i=0; i<in->sourceList->number; i++) {
       if (in->sourceList->SUlist[i]->SourID==row->SourID) 
	 {isou = i; break;}
     }
     /* if not found - skip */
     if (isou<0) continue;
     /* if no velocity info found - skip */
     if (in->sourceList->SUlist[isou]->RestFreq[0]<=0.0) continue;

     /* Update */
     for (i=0; i<numIF; i++) {
       row->LSRVel[i]   = in->sourceList->SUlist[isou]->LSRVel[i];
       row->RestFreq[i] = in->sourceList->SUlist[isou]->RestFreq[i];
     }

    /* Rewrite */
    retCode = ObitTableSUWriteRow (SUTab, iRow, row, err);
    if ((retCode != OBIT_IO_OK) || (err->error))goto cleanup;
  } /* end loop over table */

  /* Close up/cleanup */
 cleanup:
  retCode = ObitTableSUClose (SUTab, err);
  SUTab   = ObitTableSUUnref(SUTab);
  row     = ObitTableSURowUnref(row);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_msg (err, routine, in->name);

} /* end UpdateSUTable */

/**
 * Calculate default reference channel 
 * Pick channel for which the desired LSR velocity is the closest at 0h
 * on the reference day.
 * \param doppler  ObitDoppler object
 * \param err      Error stack, returns if not empty.
 * \return  default reference channel
 */
static ofloat defaultRefChan(ObitDoppler *doppler, ObitErr *err)
{
  ofloat refChan = 0.0;
  ofloat ut;
  odouble appFreq, shift, cenFreq;
  olong year, doy, ishift, i, iant;
  ObitUVDesc *uvDesc = doppler->uvdata->myDesc;
  gchar *routine = "defaultRefChan";

  /* Need position - Single source? get from header */
  if (uvDesc->ilocsu<0) {
    ObitUVGetRADec (doppler->uvdata, 
		      &doppler->source->RAMean, &doppler->source->DecMean, 
		      err);
      if (err->error) Obit_traceback_val (err, routine, doppler->name, refChan);
      /* Compute apparent position */
      ObitPrecessUVJPrecessApp (uvDesc, doppler->source);
   } else { /* Multi source */
     /* Use first in source list */
     doppler->source->RAMean  = doppler->sourceList->SUlist[0]->RAMean;
     doppler->source->DecMean = doppler->sourceList->SUlist[0]->DecMean;
     doppler->source->RAApp   = doppler->sourceList->SUlist[0]->RAApp;
     doppler->source->DecApp  = doppler->sourceList->SUlist[0]->DecApp;
    }

  /* Day */
  ObitDopplerJD2Date (doppler->JDref, &year, &doy, &ut);

  /* Find antenna with position */
  iant = 0;
  for (i=0; i<doppler->antList->number; i++) {
    if (fabs(doppler->antList->ANlist[i]->AntXYZ[0])>100.0) 
      {iant=i; break;}
  }

  /* Get apparent frequency of line at given velocity at this time, first antenna */
  appFreq = ObitDopplerFreqLSR (doppler->RestFreq, doppler->VelLSR, 
				doppler->source->RAApp, doppler->source->DecApp,  
				year, doy, ut, 
				doppler->antList->ANlist[iant]->AntXYZ[0],
				doppler->antList->ANlist[iant]->AntXYZ[1],
				doppler->antList->ANlist[iant]->AntXYZ[2]);
  /* How many channel shifted from center frequency? */
  cenFreq = uvDesc->crval[uvDesc->jlocf] + 
    ((uvDesc->inaxes[uvDesc->jlocf]/2) - uvDesc->crpix[uvDesc->jlocf]) *
    uvDesc->cdelt[uvDesc->jlocf];
  shift = (appFreq-cenFreq) / uvDesc->cdelt[uvDesc->jlocf];

  /* Minimal shift from center */
  if (shift>0.0) ishift = (olong)(shift+0.5);
  else           ishift = (olong)(shift-0.5);
  refChan = uvDesc->crpix[uvDesc->jlocf] + (uvDesc->inaxes[uvDesc->jlocf]/2) - 1 + ishift;

  return refChan;
} /* end defaultRefChan */

/** 
 * Convert Time in days to a human readable form "dd/hh:mm:ss.s"
 * \param time  Time in days
 * \param timeString [out] time as string, should be >16 char
 */
void day2dhms(ofloat time, gchar *timeString)
{
  olong day, thour, tmin;
  ofloat ttim, ssec;

  /* Trap bad times */
  if ((time<-100.0) || (time>1000.0)) {
    sprintf (timeString, "Bad time");
    return;
  }

  day   = (olong)(time);
  ttim  = 24.0*(time - day);
  thour = MIN ((olong)(ttim), 23);
  ttim  = 60.0*(ttim - thour);
  tmin  = MIN ((olong)(ttim), 59);
  ssec  = 60.0*(ttim - tmin);
  /* avoid silliness */
  if (ssec>59.951) {
    tmin++;
    ssec = 0.0;
  }
  if (tmin>=60) {
    thour++;
    tmin -= 60;
  }
  if (thour>23) {
    day++;
    thour -= 24;
  }
  sprintf (timeString, "%2.2d/%2.2d:%2.2d:%4.1f", day, thour, tmin, ssec);
  /* Zero fill seconds */
  if (timeString[9]==' ') timeString[9] = '0';
} /* end day2dhms */

