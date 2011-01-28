/* $Id$        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2010                                               */
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

#include "ObitOpacity.h"
#include "ObitPrecess.h"
#include "ObitTableANUtil.h"
#include "ObitTableSUUtil.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitOpacity.c
 * ObitOpacity class function definitions.
 * This class is derived from the Obit base class.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitOpacity";

/** Function to obtain parent ClassInfo */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo structure ObitOpacityClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitOpacityClassInfo myClassInfo = {FALSE};

/*--------------- File Global Variables  ----------------*/


/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitOpacityInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitOpacityClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitOpacityClassInfoDefFn (gpointer inClass);

/** Private: Startup. */
static void StartOpacityCalc (ObitOpacity *in, ObitErr *err);

/** Private: Calculate saturation water pressure. */
static ofloat Saturate (ofloat Temp, ofloat Press);

/** Private: Calculate zenith opacity from a combination of weather 
    and seasonal data. */
static ofloat ZenOpacity (odouble Freq, ofloat temp, ofloat DP, ofloat press, 
			  ofloat WXWeight, odouble mjd, ofloat time, 
			  ofloat *nH2O);

/** Private: Scale K band opacity to another frequency. */
static ofloat ScaleKOpac (odouble Freq, ofloat TauK, ofloat Elev);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitOpacity* newObitOpacity (gchar* name)
{
  ObitOpacity* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitOpacityClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitOpacity));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitOpacityInit((gpointer)out);

 return out;
} /* end newObitOpacity */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitOpacityGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitOpacityClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitOpacityGetClass */

/**
 * Make a deep copy of an ObitOpacity.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitOpacity* ObitOpacityCopy  (ObitOpacity *in, ObitOpacity *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  gboolean oldExist;
  gchar *outName;
  olong i;
  gchar *routine = "ObitOpacityCopy";

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
    out = newObitOpacity(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  /* Clear old */
  out->myData  = ObitUVUnref(out->myData);
  out->weather = ObitWeatherUnref(out->weather);
  if (out->antList) {
    for (i=0; i<out->numSubA; i++) 
      out->antList[i] = ObitAntennaListUnref(out->antList[i]);
    g_free(out->antList);
  }
  out->souList   = ObitSourceListUnref(out->souList);
  out->oneSource = ObitSourceUnref(out->oneSource);

  /*  clone this class */
  out->myData = ObitUVRef(in->myData);
  out->weather = ObitWeatherCopy(in->weather, out->weather, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, out);
  out->WeaWt   = in->WeaWt;
  out->numSubA = in->numSubA;
  if (out->numSubA>0) {
    out->antList = g_malloc0(out->numSubA*sizeof(ObitAntennaList*));
    for (i=0; i<out->numSubA; i++) 
      out->antList[i] = ObitAntennaListRef(in->antList[i]);
  }
  if (in->souList)   out->souList   = ObitSourceListRef(in->souList);
  if (in->oneSource) out->oneSource = ObitSourceRef(in->oneSource);

  return out;
} /* end ObitOpacityCopy */

/**
 * Make a copy of a object but do not copy the actual data
 * This is useful to create an Opacity similar to the input one.
 * \param in  The object to copy
 * \param out An existing object pointer for output, must be defined.
 * \param err Obit error stack object.
 */
void ObitOpacityClone  (ObitOpacity *in, ObitOpacity *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  olong i;
  gchar *routine = "ObitOpacityClone";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitIsA(out, &myClassInfo));

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /* Clear old */
  out->myData  = ObitUVUnref(out->myData);
  out->weather = ObitWeatherUnref(out->weather);
  if (out->antList) {
    for (i=0; i<out->numSubA; i++) 
      out->antList[i] = ObitAntennaListUnref(out->antList[i]);
    g_free(out->antList);
  }
  out->souList   = ObitSourceListUnref(out->souList);
  out->oneSource = ObitSourceUnref(out->oneSource);

  /*  clone this class */
  out->myData = ObitUVRef(in->myData);
  ObitWeatherClone(in->weather, out->weather, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  out->WeaWt   = in->WeaWt;
  out->numSubA = in->numSubA;
  if (out->numSubA>0) {
    out->antList = g_malloc0(out->numSubA*sizeof(ObitAntennaList*));
    for (i=0; i<out->numSubA; i++) 
      out->antList[i] = ObitAntennaListRef(in->antList[i]);
  }
  if (in->souList)   out->souList   = ObitSourceListRef(in->souList);
  if (in->oneSource) out->oneSource = ObitSourceRef(in->oneSource);

} /* end ObitOpacityClone */

/**
 * Creates an ObitOpacity 
 * \param name  An optional name for the object.
 * \param weather  optional weather interpolator
 * \param WeaWt    Fractional weight for weather info vs seasonal model
 *                 set to zero if weather==NULL 
 * \return the new object.
 */
ObitOpacity* ObitOpacityCreate (gchar* name, ObitUV *inData)
{
  ObitOpacity* out;

  /* Create basic structure */
  out = newObitOpacity (name);

  out->myData = ObitUVRef(inData);

  return out;
} /* end ObitOpacityCreate */

/**
 *  Calculate opacities at a set of frequencies for a source at a given
 *  time viewed by a given antenna.
 * \param in       Opacity object
 *       Control parameters:
 *   \li WXWeight OBIT_float scalar fractional weight given measured weather
 *                data vs a seasonal model.  0=>ignore weather, 1.0=>only WX. 
                  default = 1.0 (ignored if no WX table in uv data).
 * \param time     Time(days) wrt AntList reference
 * \param nfreq    Number of frequencies
 * \param freqs    List of frequencies for opacities (Hz)
 * \param Ant      Which antenna (1-rel)
 * \param SubA     Which subarray (1-rel)
 * \param Sou      Source Id
 * \param opac     [out] Opacity gain factor at freqs at observed elevation
 * \param err      Obit Error/message stack.
 */
void ObitOpacityCalc (ObitOpacity *in, ofloat time, olong nfreq, odouble *freqs,
		      olong Ant, olong SubA, olong Sou, ofloat *opac, ObitErr *err)
{
  ofloat elev, temp=0.0, press=0.0, DP=0.0, nH2O, WXWeight, stdOpac;
  ofloat windDir, windVel, wvrH2O, ions;
  odouble KBand = 23.0e9;
  ObitSource *source;
  olong i;
  ofloat fblank = ObitMagicF();
  gchar *routine = "ObitOpacityCalc";
 
 /* error checks */
  if (err->error) return;
  g_assert (ObitOpacityIsA(in));
  g_assert (freqs!=NULL);
  g_assert (opac!=NULL);

  /* Have I been properly initialized? */
  if (!in->amInit) StartOpacityCalc (in, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Multi or single source? */
  if (in->souList) source = in->souList->SUlist[Sou-1];
  else source = in->oneSource;

  /* Get elevation */
  elev = ObitAntennaListElev(in->antList[SubA-1], Ant, time, source);

  /* Get weather info */
  if ((in->WeaWt>0.0) && (in->weather))
    ObitWeatherReport(in->weather, time, Ant, SubA,  &temp, &DP, &press, 
		      &windDir, &windVel, &wvrH2O, &ions, err);
  else DP = temp = 0.0;
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* check for blanked */
  WXWeight = in->WeaWt;
  if ((temp==fblank) || (DP==fblank) || (press==fblank)) WXWeight = 0.0;

  /* Loop over frequencies */
  for (i=0; i<nfreq; i++) {
    /* K band opacity */
    stdOpac = ZenOpacity(KBand, temp, DP, press, WXWeight, 
			 in->myData->myDesc->JDObs, time,
			 &nH2O);
    /* Scale to observed frequency */
    opac[i] = ScaleKOpac (freqs[i], stdOpac, elev);
  }

} /* end ObitOpacityCalc */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitOpacityClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitOpacityClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitOpacityClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitOpacityClassInfoDefFn (gpointer inClass)
{
  ObitOpacityClassInfo *theClass = (ObitOpacityClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitOpacityClassInit;
  theClass->newObit       = (newObitFP)newObitOpacity;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitOpacityClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitOpacityGetClass;
  theClass->ObitCopy      = (ObitCopyFP)ObitOpacityCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitOpacityClear;
  theClass->ObitInit      = (ObitInitFP)ObitOpacityInit;
  theClass->ObitOpacityCreate = (ObitOpacityCreateFP)ObitOpacityCreate;
  theClass->ObitOpacityCalc   = (ObitOpacityCalcFP)ObitOpacityCalc;

} /* end ObitOpacityClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitOpacityInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitOpacity *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->amInit    = FALSE;
  in->myData    = NULL;
  in->weather   = NULL;
  in->WeaWt     = 0.0;
  in->antList   = NULL;
  in->numSubA   = 0;
  in->souList   = NULL;
  in->oneSource = NULL;

} /* end ObitOpacityInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitOpacity* cast to an Obit*.
 */
void ObitOpacityClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitOpacity *in = inn;
  olong i;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  in->weather   = ObitWeatherUnref(in->weather);
  in->myData    = ObitUVUnref(in->myData);
  in->souList   = ObitSourceListUnref(in->souList);
  in->oneSource = ObitSourceUnref(in->oneSource);
  if (in->antList) {
    for (i=0; i<in->numSubA; i++)
      in->antList[i]   = ObitAntennaListUnref(in->antList[i]);
    g_free(in->antList);
  }
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitOpacityClear */

/**
 *  Startup
 *  \li Determine if WX table use desired
 *  \li Initialize weather interpolator if needed
 *  \li Get Antenna list(s)
 *  \li Get Source list if multi source or,
 *  \li Set oneSource for a single source file
 * \param in       Opacity object
 * \param err      Obit Error/message stack.
 */
static void StartOpacityCalc (ObitOpacity *in, ObitErr *err)
{
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong        ver, iver, nrow, highVer;
  olong        numIF=0, numOrb=0, numPCal=0;
  ObitTableAN  *ANTable=NULL;
  ObitTableWX  *WXTable=NULL;
  ObitTableSU  *SUTable=NULL;
  gchar *routine = "ObitOpacity:StartOpacityCalc ";
 
  /* error checks */
  if (err->error) return;

  in->amInit = TRUE;  /* This call initializes */

  /* Is use of the WX table requested? */
  in->WeaWt = 1.0;
  ObitInfoListGetTest(in->myData->info, "WXWeight", &type, dim, &in->WeaWt);

  /* Initialize weather interpolator if requested */
  if (in->WeaWt>0.0) {

    /* Have WX table? */
    highVer = ObitTableListGetHigh (in->myData->tableList, "AIPS WX");
    if (highVer<1) in->WeaWt = 0.0;   /* No option here */
    
    /* Is anything in the table? */
    if (in->WeaWt>0.0) {
      ver = 1;
      WXTable = 
	newObitTableWXValue (in->name, (ObitData*)in->myData, &ver, OBIT_IO_ReadOnly, err);
      ObitTableWXOpen (WXTable, OBIT_IO_ReadOnly, err);
      nrow = WXTable->myDesc->nrow;
      ObitTableWXClose (WXTable, err);
      if (err->error) Obit_traceback_msg (err, routine, in->name);
      if (nrow<=0) in->WeaWt = 0.0;  /* Anything there? */

      if (in->WeaWt>0.0) {  /* Keep trying? */
	in->weather = ObitWeatherCreate("Weather", WXTable, in->myData, err);
	if (err->error) Obit_traceback_msg (err, routine, in->name);
     }
    }
    WXTable = ObitUnref(WXTable);   /* Cleanup */
  } /* end initialize weather interpolator */

  /* Antenna lists - how name subarrays? */
  highVer = ObitTableListGetHigh (in->myData->tableList, "AIPS AN");
  in->numSubA = highVer;
  in->antList = g_malloc0(in->numSubA*sizeof(ObitAntennaList*));

  /* Loop over AN tables (subarrays) forming antenna lists*/
  for (iver=1; iver<=in->numSubA; iver++) {
    ver = iver;
    ANTable = newObitTableANValue ("AN table", (ObitData*)in->myData, 
				   &ver, OBIT_IO_ReadOnly, numIF, numOrb, numPCal, err);
    if (ANTable==NULL) Obit_log_error(err, OBIT_Error, "ERROR with AN table");
    in->antList[iver-1] = ObitTableANGetList (ANTable, err);
    if (err->error) Obit_traceback_msg (err, routine, in->myData->name);
    
    ANTable = ObitTableANUnref(ANTable);   /* Cleanup */
  } /* End loop over subarrays */

  /* Is in->myData a single or multi source file? */
  highVer = ObitTableListGetHigh (in->myData->tableList, "AIPS SU");

  if (highVer>0) {
    /* Multisource, get Source List */
  /* Convert SU table into Source List */
    ver = 1;
    if (in->myData->myDesc->jlocif>=0) 
      numIF = in->myData->myDesc->inaxes[in->myData->myDesc->jlocif];
    else numIF = 1;
    SUTable = newObitTableSUValue (in->myData->name, (ObitData*)in->myData, 
				   &ver, numIF, OBIT_IO_ReadOnly, err);
    if (SUTable) in->souList = ObitTableSUGetList (SUTable, err);
    if (err->error) Obit_traceback_msg (err, routine, in->myData->name);
    SUTable = ObitTableSUUnref(SUTable);
  } else {
    /* single source - create oneSource */
    in->oneSource = newObitSource("Single");
    strncpy (in->oneSource->SourceName, in->myData->myDesc->object, 
	     MIN(20,UVLEN_VALUE));
    in->oneSource->equinox = in->myData->myDesc->equinox;
    in->oneSource->RAMean  = in->myData->myDesc->crval[in->myData->myDesc->jlocr];
    in->oneSource->DecMean = in->myData->myDesc->crval[in->myData->myDesc->jlocd];
    /* Compute apparent position */
    ObitPrecessUVJPrecessApp (in->myData->myDesc, in->oneSource);
  }

} /* end StartOpacityCalc */

/**
 * Calculate zenith opacity from a combination of weather 
 * and seasonal data.
 * From the AIPSIsh OPACTY.FOR
 *   Function to calculate day of year from month, day, year from Meeus'
 *   book on Astronomical Algorithms.
 * \param  Freqs     Frequency (Hz)
 * \param  temp      Surface ambient temp (C)
 * \param  DP        Surface dewpoint temp (C)
 * \param  press     Surface pressure in millibars
 * \param  WXWeight  Weight to give weather (0 -> 1), 1-WT given seasonal
 * \param  JD        Current Julian day
 * \param  time      Average IAT after JD (day)
 * \param  nH2O     [out] H2O column density  (requires WT > 0)
 * \return  zenith opacity
 */
static ofloat ZenOpacity (odouble Freq, ofloat temp, ofloat DP, ofloat press, 
			  ofloat WXWeight, odouble JD, ofloat time, 
			  ofloat *nH2O)
{
  ofloat Opac    = 0.0;
  ofloat kk     = 1.380658e-23;
  ofloat mvap   = 18 * 1.6749286e-27;
  ofloat watden = 1000.0;
  ofloat vaph   = 1500.0;    /* vaph is the scale height (m) of H2O vapor */
  ofloat konst, pvap, mdoy, wopac, sopac;
  olong doy, month, day, year, k2, iband, ic;
  odouble j, y, d, m;
  gboolean isLeap;
  ofloat C[24]= { 
    0.0,  0.0,  0.0,           /* 4 band  */
    0.0,  0.0,  0.0,           /* P-band  */
    1.0,  0.0,  0.0,           /* L-band  */
    1.0,  0.0,  0.0,           /* C-band  */
    1.0,  0.0,  0.0,           /* X-band  */
    1.19, 0.0578, -2.486E-4,   /* Ku-band */
    1.77, 0.9060, -1.385E-4,   /* K-band  */
    5.86, 0.1712,  9.201E-4    /* Q-band  */ 
  };
  ofloat D[24]= { 
    0.0,  0.0,  0.0,           /* 4 band  */
    0.0,  0.0,  0.0,           /* P-band  */
    1.0,  0.0,  0.0,           /* L-band  */
    1.0,  0.0,  0.0,           /* C-band  */
    1.0,  0.0,  0.0,           /* X-band  */
    2.19, -0.0085, 2.010E-5,   /* Ku-band */
    22.12, -0.1783, 4.402E-4,  /* K-band  */
    9.70, -0.0357, 9.402E-5    /* Q-band  */ 
  };

  /* Work out band */
  iband = 0;                     /* 4 band  */
  if (Freq>=200.0e6)  iband=1;   /* P-band  */
  if (Freq>= 1.0e9)   iband=2;   /* L-band  */
  if (Freq>= 3.0e9)   iband=3;   /* C-band  */
  if (Freq>= 7.0e9)   iband=4;   /* X-band  */
  if (Freq>=10.0e9)   iband=5;   /* Ku-band */
  if (Freq>=20.0e9)   iband=6;   /* K-band  */
  if (Freq>=36.0e9)   iband=7;   /* Q-band  */
 
 konst  = 1.0e5*mvap*vaph/(watden*kk);

  /*  weather opacity */
   if (WXWeight>0.0) {
    /* No, I know this looks wrong but Brian Butler claims it's right */
    pvap = Saturate (DP, press);
    (*nH2O) = konst * pvap / (temp + 273.15);
    wopac = (C[iband*3] + C[iband*3+1]*(*nH2O) + C[iband*3+2]*(*nH2O)*(*nH2O)) / 100.0;
  }  else {
    wopac = 0.0;
  }
  
  /*  seasonal opacity  */
  if (WXWeight<1.0) {
    (*nH2O) = 0.0;   /* No can do */
    /* Convert JD to date */
    j = (olong) (JD + 0.50 - 1721119.0);
    y = (olong) ((4.0*j - 1.00) / 146097.0);
    ic = y + 0.00010;
    j = 4.0*j - 1.00 - 146097.0*y;
    d = (olong) (j * 0.250);
    j = (olong) ((4.00*d + 3.00) / 1461.00);
    d = 4.00*d + 3.00 - 1461.00*j;
    d = (olong) ((d+4.00) * 0.250);
    m = (olong) ((5.00*d - 3.00) / 153.00);
    day = 5.00*d - 3.00 - 153.00*m;
    day = (day + 5) / 5;
    year = j + 100*ic;
    month = m;
    if (month < 10) {
      month += 3;
    } else {
      month -= 9;
      year += 1;
    }
    
    isLeap = ((year%4)==0) && (((year%100)!=0) || ((year%400)==0));
    if (isLeap)  k2 = 1;
    else         k2 = 2;
    
    doy = ((275 * month) / 9) - k2 * ((month + 9) / 12) + day - 30;
    mdoy = doy;
    if (mdoy>=200.0) mdoy = mdoy - 365.0;
    mdoy += 165.0;
    mdoy += time/24.0;
    sopac = (D[iband*3] + D[iband*3+1]*mdoy + D[iband*3+2]*mdoy*mdoy) / 100.0;
  } else {
    sopac = 0.0;
  }
  
  /* Put it all together */
  Opac = MAX (0.0, MIN (1.0, WXWeight)) * wopac + MAX (0.0, MIN (1.0, 1.0-WXWeight)) * sopac;
  
  return Opac;
} /* end ZenOpacity */

			    /**
			     * Function to calculate saturation vapor pressure over water or ice,
 *   in millibar, given temperature (Temp) in degrees C, and pressure
 *  (Press) in millibar.  
 * From the AIPSIsh SATPRS.FOR
 *   Taken from:
 *     Buck, A.L., New Equations for Computing Vapor Pressure and
 *       Enhancement Factor, J. Appl. Met., v.20, pp. 1527-1532, 1981
 *  who references:
 *     Wexler, A., Vapor Pressure Formulation for Water in the Range
 *        0C to 100C - A Revision, J. Res. Natl. Bur. Stand., v.80A,
 *        pp. 775 ff, 1976
 *  and
 *     Wexler, A., Vapor Pressure Formulation for Ice, J. Res. Natl.
 *       Bur. Stand., v.81A, pp. 5-20, 1977
 *  for the "exact" formulations - which are reworks of the Goff-Gratch
 *  formulation:
 *     Goff, J.A., and S. Gratch, Low-pressure Properties of Water from
 *       -160F to 212F. Trans. Am. Soc. Heat. Vent. Eng., v. 52, 95-121,
 *       1946
 *  note that i use the fw5 and fi5 coefficients for the "enhancement
 *  factor" from Table 3 of Buck.
 *  By Bryan Butler, NRAO
 * \param  Temp  Temperature (C)
 * \param  Press Atmospheric pressure (millibar)
 * \return vapor pressure of saturated water 
 */
static ofloat Saturate (ofloat Temp, ofloat Press)
{
  odouble satp = 0.0;
  odouble theta, theta2, theta3, theta4, ew, ei, fw, fi;
  odouble aow=4.1e-4, aoi=4.8e-4, bow=3.48e-6, boi=3.47e-6, 
    cow=7.4e-10, coi=5.9e-10, dow=30.6, doi=23.8, eow=-3.8e-2, eoi= -3.1e-2;

  theta  = Temp + 273.15;
  theta2 = theta*theta;
  theta3 = theta*theta2;
  theta4 = theta*theta3;

  /*   !!! water !!!  */
  if (Temp>0.01) {
    ew = (-2991.2729 / theta2) + (-6017.0128 / theta) +
      (18.87643854) + (-0.028354721 * theta) +(0.17838301e-4 * theta2) +
      (-0.84150417e-9 * theta3) + (0.44412543e-12 * theta4) +
      (2.858487 * log(theta));
    ew = 0.01 * exp (ew);
    fw = 1.0 + aow + Press * 
      (bow + cow * (Temp + dow + eow * Press)*(Temp + dow + eow * Press));
    satp = ew * fw;
  } else {
    /*     !!! ice !!!  */
    ei = (-5865.3696 / theta) + (22.241033) + (0.013749042 * theta) +
      (-0.34031775e-4 * theta2) + (0.26967687e-7 * theta3) + (0.6918651 * log(theta));
    ei = 0.01 * exp (ei);
    fi = 1.0 + aoi + Press * 
      (boi + coi * (Temp + doi + eoi * Press)*(Temp + doi + eoi * Press));
    satp = ei * fi;
  }
  return satp;
} /* end Saturate */

/**
 * Scale K band opacity to other frequencies at a given elevation
 * Only valid to 1-50 GHz
 * From the AIPSIsh KBOPAC.FOR
 *  Function derived by Josh Marvil and Frazer Owen (2010)
 * \param  nFreq   Number of frequencies in Freqs  
 * \param  Freq    Frequency (Hz)
 * \param  TauK    (radio) K band zenith opacity (nepers)
 * \param  Elev    Elevation (rad)
 * \return  Opac Scaled opacity at Freqs, at observed elevation
 */
static ofloat ScaleKOpac (odouble Freq , ofloat TauK, ofloat Elev)
{
  ofloat Opac=0.0;
  olong j;
  ofloat z1, z2, zopac, pwv, d;
  /* Values for 1.0 - 50.0 GHz interval of 0.25 GHz */
  ofloat TauFun[2*197] = {
    5.272, 0.000,   5.464, 0.000,   5.579, 0.001,
    5.654, 0.002,   5.707, 0.003,   5.748, 0.005,
    5.780, 0.006,   5.808, 0.008,   5.832, 0.009,
    5.855, 0.011,   5.876, 0.013,   5.897, 0.016,
    5.918, 0.018,   5.939, 0.021,   5.960, 0.023,
    5.981, 0.026,   6.003, 0.030,   6.026, 0.033,
    6.049, 0.037,   6.073, 0.040,   6.098, 0.044,
    6.123, 0.048,   6.150, 0.052,   6.178, 0.057,
    6.206, 0.061,   6.235, 0.066,   6.266, 0.071,
    6.297, 0.077,   6.330, 0.082,   6.364, 0.088,
    6.398, 0.094,   6.435, 0.100,   6.473, 0.107,
    6.522, 0.114,   6.550, 0.122,   6.590, 0.129,
    6.635, 0.137,   6.726, 0.146,   6.722, 0.155,
    6.766, 0.164,   6.835, 0.174,   6.864, 0.185,
    6.905, 0.196,   6.954, 0.207,   7.003, 0.220,
    7.055, 0.233,   7.109, 0.247,   7.164, 0.262,
    7.220, 0.277,   7.277, 0.294,   7.336, 0.312,
    7.396, 0.331,   7.458, 0.352,   7.521, 0.374,
    7.586, 0.398,   7.661, 0.423,   7.726, 0.451,
    7.785, 0.482,   7.854, 0.516,   7.924, 0.553,
    8.000, 0.593,   8.081, 0.638,   8.144, 0.688,
    8.219, 0.744,   8.296, 0.806,   8.375, 0.876,
    8.456, 0.955,   8.540, 1.045,   8.634, 1.148,
    8.724, 1.265,   8.819, 1.400,   8.920, 1.557,
    9.028, 1.738,   9.147, 1.949,   9.281, 2.195,
    9.435, 2.483,   9.615, 2.820,   9.830, 3.214,
    10.090, 3.670,  10.404, 4.194,  10.776, 4.786,
    11.204, 5.436,  11.669, 6.117,  12.129, 6.778,
    12.523, 7.328,  12.782, 7.638,  12.872, 7.625,
    12.815, 7.351,  12.667, 6.924,  12.491, 6.426,
    12.335, 5.914,  12.256, 5.420,  12.146, 4.961,
    12.134, 4.546,  12.171, 4.176,  12.249, 3.849,
    12.362, 3.562,  12.503, 3.312,  12.689, 3.093,
    12.859, 2.902,  13.010, 2.735,  13.210, 2.589,
    13.419, 2.461,  13.638, 2.349,  13.864, 2.250,
    14.097, 2.163,  14.342, 2.087,  14.586, 2.019,
    14.838, 1.960,  15.098, 1.907,  15.367, 1.861,
    15.653, 1.820,  16.030, 1.784,  16.215, 1.752,
    16.506, 1.724,  16.820, 1.700,  17.345, 1.679,
    17.567, 1.660,  17.879, 1.644,  18.096, 1.630,
    18.423, 1.619,  18.769, 1.609,  19.127, 1.601,
    19.495, 1.595,  19.871, 1.590,  20.262, 1.587,
    20.663, 1.585,  21.078, 1.584,  21.501, 1.584,
    21.939, 1.585,  22.390, 1.587,  22.853, 1.590,
    23.256, 1.594,  23.747, 1.598,  24.254, 1.604,
    24.777, 1.610,  25.317, 1.616,  25.875, 1.624,
    26.451, 1.631,  27.058, 1.640,  26.713, 1.649,
    27.065, 1.658,  27.674, 1.668,  28.319, 1.678,
    28.991, 1.689,  29.692, 1.700,  30.431, 1.712,
    31.355, 1.724,  31.152, 1.736,  31.939, 1.749,
    32.756, 1.762,  33.606, 1.776,  34.491, 1.789,
    35.414, 1.803,  36.384, 1.818,  37.363, 1.832,
    38.973, 1.847,  40.074, 1.862,  41.222, 1.878,
    42.420, 1.894,  43.671, 1.910,  44.978, 1.926,
    46.345, 1.942,  47.775, 1.959,  52.349, 1.976,
    53.994, 1.993,  55.726, 2.011,  57.631, 2.028,
    59.460, 2.046,  61.421, 2.064,  63.577, 2.083,
    65.854, 2.101,  68.610, 2.120,  71.062, 2.139,
    73.653, 2.158,  76.416, 2.178,  79.303, 2.197,
    82.334, 2.217,  85.573, 2.237,  89.014, 2.257,
    92.743, 2.277,  96.640, 2.298, 100.796, 2.319,
    105.234, 2.340, 109.982, 2.361, 115.072, 2.382,
    120.538, 2.404, 126.423, 2.426, 132.773, 2.448,
    139.643, 2.470, 147.095, 2.492, 155.207, 2.515,
    164.067, 2.538, 173.788, 2.561, 184.509, 2.584,
    196.415, 2.608, 209.874, 2.632
  };

  pwv = -1.710 + 136.47 * TauK;
  j = Freq /0.25e9 + 0.01;
  j = j - 3;
  j = MAX (0, MIN (195, j));  /* Clip to range */
  d = Freq/0.25e9 - (j+3);
  z1 = TauFun[j*2] + TauFun[j*2+1] * pwv;
  z2 = TauFun[(j+1)*2] + TauFun[(j+1)*2+1] * pwv;
  zopac = (z1 + d * (z2 - z1)) / 1000.0;
  Opac = zopac / sin (Elev);
  Opac = sqrt (exp (Opac));

  return Opac;
} /* end ScaleKOpac */
