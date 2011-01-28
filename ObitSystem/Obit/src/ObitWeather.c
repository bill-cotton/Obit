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

#include "ObitWeather.h"
#include "ObitTableANUtil.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitWeather.c
 * ObitWeather class function definitions.
 * This class is derived from the Obit base class.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitWeather";

/** Function to obtain parent ClassInfo */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo structure ObitWeatherClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitWeatherClassInfo myClassInfo = {FALSE};

/*--------------- File Global Variables  ----------------*/


/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitWeatherInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitWeatherClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitWeatherClassInfoDefFn (gpointer inClass);

/** Private: Initialize for a subarray. */
static void InitSubarray (ObitWeather *in, olong subarray, ObitErr *err);

/** Private: Update time of data. */
static void UpdateTime (ObitWeather *in, ofloat time, olong Ant, 
			ObitErr *err);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitWeather* newObitWeather (gchar* name)
{
  ObitWeather* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitWeatherClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitWeather));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitWeatherInit((gpointer)out);

 return out;
} /* end newObitWeather */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitWeatherGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitWeatherClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitWeatherGetClass */

/**
 * Make a deep copy of an ObitWeather.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitWeather* ObitWeatherCopy  (ObitWeather *in, ObitWeather *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  gboolean oldExist;
  olong i;
  gchar *outName;
  gchar *routine = "ObitWeatherCopy";

  /* error checks */
  if (err->error) return out;
  g_assert (ObitIsA(in, &myClassInfo));
  if (out) g_assert (ObitIsA(out, &myClassInfo));

  /* Create if it doesn't exist */
  oldExist = out!=NULL;
  if (!oldExist) {
    /* derive object name */
    outName = g_strconcat ("Copy: ",in->name,NULL);
    out = newObitWeather(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  /* Clear old structures */
  out->WXTable = ObitTableWXUnref(out->WXTable);
  if (out->antList) {
    for (i=0; i<out->numSubA; i++)
      out->antList[i]   = ObitAntennaListUnref(out->antList[i]);
    g_free(out->antList);
  }
  if (out->priorRowNo)  g_free(out->priorRowNo);
  if (out->followRowNo) g_free(out->followRowNo);
  if (out->priorRow){
    for (i=0; i<out->numAnt; i++) {
      ObitTableWXRowUnref(out->priorRow[i]);
    }
    g_free(out->priorRow);
  }
  if (out->followRow){
    for (i=0; i<out->numAnt; i++) {
      ObitTableWXRowUnref(out->followRow[i]);
    }
    g_free(out->followRow);
  }

  /* Create new  */
  out->WXTable  = ObitTableWXRef(in->WXTable);
  out->numSubA  = in->numSubA;
  out->antList  = g_malloc0(out->numSubA*sizeof(ObitAntennaList*));
  for (i=0; i<out->numSubA ; i++) {
    out->antList[i] = ObitAntennaListRef(in->antList[i]);
  }
  out->SubA        = in->SubA;
  out->numAnt      = out->antList[out->SubA-1]->number;
  out->priorRowNo  = g_malloc0(out->numAnt*sizeof(olong));
  out->followRowNo = g_malloc0(out->numAnt*sizeof(olong));
  out->priorRow    = g_malloc0(out->numAnt*sizeof(ObitTableWXRow*));
  out->followRow   = g_malloc0(out->numAnt*sizeof(ObitTableWXRow*));
  for (i=0; i<out->numAnt; i++) {
    out->priorRow[i]  = newObitTableWXRow(out->WXTable);
    out->followRow[i] = newObitTableWXRow(out->WXTable);
    /* Set time to large negative value to indicate uninitialized */
    out->priorRow[i]->Time    = -1.0e20;
    out->followRow[i]->Time   = -1.0e20;
    out->priorRowNo[i]        = -1;
    out->followRowNo[i]       = -1;
  }

  /* Init for Subarray out->SubA */
  InitSubarray(out, out->SubA, err);
  if (err->error) Obit_traceback_val (err, routine, out->name, out);

  return out;
} /* end ObitWeatherCopy */

/**
 * Make a copy of a object but do not copy the actual data
 * This is useful to create an Weather similar to the input one.
 * \param in  The object to copy
 * \param out An existing object pointer for output, must be defined.
 * \param err Obit error stack object.
 */
void ObitWeatherClone  (ObitWeather *in, ObitWeather *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  olong i;
  gchar *routine = "ObitWeatherClone";

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitIsA(out, &myClassInfo));

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  /* Clear old structures */
  out->WXTable = ObitTableWXUnref(out->WXTable);
  if (out->antList) {
    for (i=0; i<out->numSubA; i++)
      out->antList[i]   = ObitAntennaListUnref(out->antList[i]);
    g_free(out->antList);
  }
  if (out->priorRowNo)  g_free(out->priorRowNo);
  if (out->followRowNo) g_free(out->followRowNo);
  if (out->priorRow){
    for (i=0; i<out->numAnt; i++) {
      ObitTableWXRowUnref(out->priorRow[i]);
    }
    g_free(out->priorRow);
  }
  if (out->followRow){
    for (i=0; i<out->numAnt; i++) {
      ObitTableWXRowUnref(out->followRow[i]);
    }
    g_free(out->followRow);
  }

  /* Create new  */
  out->WXTable  = ObitTableWXRef(in->WXTable);
  out->numSubA  = in->numSubA;
  out->antList  = g_malloc0(out->numSubA*sizeof(ObitAntennaList*));
  for (i=0; i<out->numSubA ; i++) {
    out->antList[i] = ObitAntennaListRef(in->antList[i]);
  }
  out->SubA        = in->SubA;
  out->numAnt      = out->antList[out->SubA-1]->number;
  out->priorRowNo  = g_malloc0(out->numAnt*sizeof(olong));
  out->followRowNo = g_malloc0(out->numAnt*sizeof(olong));
  out->priorRow    = g_malloc0(out->numAnt*sizeof(ObitTableWXRow*));
  out->followRow   = g_malloc0(out->numAnt*sizeof(ObitTableWXRow*));
  for (i=0; i<out->numAnt; i++) {
    out->priorRow[i]  = newObitTableWXRow(out->WXTable);
    out->followRow[i] = newObitTableWXRow(out->WXTable);
    /* Set time to large negative value to indicate uninitialized */
    out->priorRow[i]->Time    = -1.0e20;
    out->followRow[i]->Time   = -1.0e20;
    out->priorRowNo[i]        = -1;
    out->followRowNo[i]       = -1;
  }

  /* Init for Subarray out->SubA */
  InitSubarray(out, out->SubA, err);
  if (err->error) Obit_traceback_msg (err, routine, out->name);

} /* end ObitWeatherClone */

/**
 * Creates an ObitWeather 
 * \param name    An optional name for the object.
 * \param WXTable Weather table to interpolate
 * \param UVData  UV data used to get antenna/array information
 * \param err     Obit error stack object.
 * \return the new object.
 */
ObitWeather* ObitWeatherCreate (gchar* name, ObitTableWX *WXTable, 
				ObitUV *UVData, ObitErr *err)
{
  ObitWeather* out;
  ObitTableAN  *ANTable=NULL;
  olong highVer, ver, iver;
  olong numIF=0, numOrb=0, numPCal=0;
  gchar *routine = "ObitWeatherCreate";

  /* Create basic structure */
  out = newObitWeather (name);

  /* Save table */
  out->WXTable = ObitTableWXRef(WXTable);

  /* Get antenna/subarray information */
  /* Antenna lists - how name subarrays? */
  highVer = ObitTableListGetHigh (UVData->tableList, "AIPS AN");
  out->numSubA = highVer;
  out->antList = g_malloc0(out->numSubA*sizeof(ObitAntennaList*));

  /* Loop over AN tables (subarrays) forming antenna lists*/
  for (iver=1; iver<=out->numSubA; iver++) {
    ver = iver;
    ANTable = newObitTableANValue ("AN table", (ObitData*)UVData, 
				   &ver, OBIT_IO_ReadOnly, numIF, numOrb, numPCal, err);
    if (ANTable==NULL) Obit_log_error(err, OBIT_Error, "ERROR with AN table");
    out->antList[iver-1] = ObitTableANGetList (ANTable, err);
    if (err->error) Obit_traceback_val (err, routine, UVData->name, out);
    
    ANTable = ObitTableANUnref(ANTable);   /* Cleanup */
  } /* End loop over subarrays */

  /* Initialize for subarray 1 */
  InitSubarray(out, 1, err);
  if (err->error) Obit_traceback_val (err, routine, UVData->name, out);

  return out;
} /* end ObitWeatherCreate */

/**
 *  Interpolate Weather values at a given time
 *  Calls are expected in time order by antenna.
 * \param in      Weather interpolator
 * \param time    Desired time (day)
 * \param ant     Desired antenna number
 * \param suba    Desired subarray
 * \param temp    [out] Temperature (C), may be fblanked
 * \param DP      [out] Dewpoint (C), may be fblanked
 * \param press   [out] Pressure (millibar), may be fblanked
 * \param windDir [out] Wind direction azimuth deg, may be fblanked
 * \param windVel [out] Wind velocity m/s, may be fblanked
 * \param wvrH2O  [out] Water vapor radiometer (??), may be fblanked
 * \param ions    [out] Ionospheric electron something (??), may be fblanked
 * \param err     Obit error stack object.
 */
void ObitWeatherReport (ObitWeather *in, ofloat time, olong ant, olong suba,
			ofloat *temp, ofloat *DP, ofloat *press, 
			ofloat *windDir, ofloat *windVel, 
			ofloat *wvrH2O, ofloat *ions, 
			ObitErr *err)
{
  ofloat w1, w2, ww1, ww2, diff;
  ofloat fblank = ObitMagicF();
  olong iant = ant-1;
  gchar *routine = "ObitWeatherReport";
  
  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  
  /* New subarray? */
  if (suba!=in->SubA) InitSubarray (in, suba, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Check/update times */
  UpdateTime (in, time, ant, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Set weights */
  /* Prior only? */
  if ((in->priorRow[iant]->Time>-1.0e3) && (in->followRow[iant]->Time<-1.0e3)) {
    w1 = 1.0;
    w2 = 0.0;
    /* If both before use follow */
  } else if ((in->priorRow[iant]->Time<time) && (in->followRow[iant]->Time<time)) { 
    w1 = 0.0;
    w2 = 1.0;
    /* If both after use prior */
  } else if ((in->priorRow[iant]->Time>time) && (in->followRow[iant]->Time>time)) { 
    w1 = 1.0;
    w2 = 0.0;
    /* Follow only - shouldn't happed but just in case */
  } else if ((in->priorRow[iant]->Time<-1.0e3) && (in->followRow[iant]->Time>-1.0e3)) {
    w1 = 0.0;
    w2 = 1.0;
  } else {  /* Both - interpolate */
    diff = in->followRow[iant]->Time-in->priorRow[iant]->Time;
    if (diff==0.0) diff = 1.0;
    w1 = (in->followRow[iant]->Time-time) / diff;
    w2 = 1.0 - w1;
  }

  /* Output data may be blanked - modify weights */
  /* temp    Temperature (C) = temperature*/
  ww1 = w1; ww2 = w2;
  if (in->priorRow[iant]->temperature==fblank)  ww1 = 0.0;
  if (in->followRow[iant]->temperature==fblank) ww2 = 0.0;
  if ((ww1+ww2)>0.0) *temp = ww1*in->priorRow[iant]->temperature + 
		       ww2*in->followRow[iant]->temperature;
  else *temp = fblank;

  /* DP      Dewpoint (C) = dewpoint */
  ww1 = w1; ww2 = w2;
  if (in->priorRow[iant]->dewpoint==fblank)  ww1 = 0.0;
  if (in->followRow[iant]->dewpoint==fblank) ww2 = 0.0;
  if ((ww1+ww2)>0.0) *DP = ww1*in->priorRow[iant]->dewpoint + 
		       ww2*in->followRow[iant]->dewpoint;
  else *DP = fblank;

  /* press   Pressure (millibar) = pressure */
  ww1 = w1; ww2 = w2;
  if (in->priorRow[iant]->pressure==fblank)  ww1 = 0.0;
  if (in->followRow[iant]->pressure==fblank) ww2 = 0.0;
  if ((ww1+ww2)>0.0) *press = ww1*in->priorRow[iant]->pressure + 
		       ww2*in->followRow[iant]->pressure;
  else *press = fblank;

  /* windDir Wind direction azimuth deg = windDirection */
  ww1 = w1; ww2 = w2;
  if (in->priorRow[iant]->windDirection==fblank)  ww1 = 0.0;
  if (in->followRow[iant]->windDirection==fblank) ww2 = 0.0;
  if ((ww1+ww2)>0.0) *windDir = ww1*in->priorRow[iant]->windDirection + 
		       ww2*in->followRow[iant]->windDirection;
  else *windDir = fblank;

  /* windVel Wind velocity m/s = windVelocity */
  ww1 = w1; ww2 = w2;
  if (in->priorRow[iant]->windVelocity==fblank)  ww1 = 0.0;
  if (in->followRow[iant]->windVelocity==fblank) ww2 = 0.0;
  if ((ww1+ww2)>0.0) *windVel = ww1*in->priorRow[iant]->windVelocity + 
		       ww2*in->followRow[iant]->windVelocity;
  else *windVel = fblank;

  /* wvrH2O  Water vapor radiometer (??) = wvrH2O */
  ww1 = w1; ww2 = w2;
  if (in->priorRow[iant]->wvrH2O==fblank)  ww1 = 0.0;
  if (in->followRow[iant]->wvrH2O==fblank) ww2 = 0.0;
  if ((ww1+ww2)>0.0) *wvrH2O = ww1*in->priorRow[iant]->wvrH2O + 
		       ww2*in->followRow[iant]->wvrH2O;
  else *wvrH2O = fblank;

  /* ions    Ionospheric electron something (??) = onosElectron */
  ww1 = w1; ww2 = w2;
  if (in->priorRow[iant]->onosElectron==fblank)  ww1 = 0.0;
  if (in->followRow[iant]->onosElectron==fblank) ww2 = 0.0;
  if ((ww1+ww2)>0.0) (*ions) = ww1*in->priorRow[iant]->onosElectron + 
		       ww2*in->followRow[iant]->onosElectron;
  else (*ions) = fblank;


} /* end ObitWeatherReport */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitWeatherClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitWeatherClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitWeatherClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitWeatherClassInfoDefFn (gpointer inClass)
{
  ObitWeatherClassInfo *theClass = (ObitWeatherClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitWeatherClassInit;
  theClass->newObit       = (newObitFP)newObitWeather;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitWeatherClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitWeatherGetClass;
  theClass->ObitCopy      = (ObitCopyFP)ObitWeatherCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitWeatherClear;
  theClass->ObitInit      = (ObitInitFP)ObitWeatherInit;
  theClass->ObitWeatherCreate = (ObitWeatherCreateFP)ObitWeatherCreate;

} /* end ObitWeatherClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitWeatherInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitWeather *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->WXTable     = NULL;
  in->numSubA     =  0;
  in->antList     = NULL;
  in->SubA        = 0;
  in->numAnt      = 0;
  in->priorRowNo  = NULL;
  in->followRowNo = NULL;
  in->priorRow    = NULL;
  in->followRow   = NULL;

} /* end ObitWeatherInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *             Actually it should be an ObitWeather* cast to an Obit*.
 */
void ObitWeatherClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitWeather *in = inn;
  olong i;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  in->WXTable = ObitTableWXUnref(in->WXTable);
  if (in->antList) {
    for (i=0; i<in->numSubA; i++)
      in->antList[i]   = ObitAntennaListUnref(in->antList[i]);
    g_free(in->antList);
  }
  if (in->priorRowNo)  g_free(in->priorRowNo);
  if (in->followRowNo) g_free(in->followRowNo);
  if (in->priorRow){
    for (i=0; i<in->numAnt; i++) {
      ObitTableWXRowUnref(in->priorRow[i]);
    }
    g_free(in->priorRow);
  }
  if (in->followRow){
    for (i=0; i<in->numAnt; i++) {
      ObitTableWXRowUnref(in->followRow[i]);
    }
    g_free(in->followRow);
  }

  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitWeatherClear */

/**
 * Initialize object for a given subarray
 * Created row structures and reads prior and follow for each antenna.
 * \param in   The object to init
 * \param SubA Which (1-rel) subarray?
 * \param err  Obit error stack object.
 */
static void InitSubarray (ObitWeather *in, olong SubA, ObitErr *err)
{
  olong i, iRow, ant, start;
  gboolean found;
  gchar *routine = "ObitWeather:InitSubarray";

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Clear old structures */
  if (in->priorRowNo)  g_free(in->priorRowNo);
  if (in->followRowNo) g_free(in->followRowNo);
  if (in->priorRow){
    for (i=0; i<in->numAnt; i++) {
      ObitTableWXRowUnref(in->priorRow[i]);
    }
    g_free(in->priorRow);
  }
  if (in->followRow){
    for (i=0; i<in->numAnt; i++) {
      ObitTableWXRowUnref(in->followRow[i]);
    }
    g_free(in->followRow);
  }

  /* Create new for this subarray */
  in->SubA        = SubA;
  in->numAnt      = in->antList[SubA-1]->number;
  in->priorRowNo  = g_malloc0(in->numAnt*sizeof(olong));
  in->followRowNo = g_malloc0(in->numAnt*sizeof(olong));
  in->priorRow    = g_malloc0(in->numAnt*sizeof(ObitTableWXRow*));
  in->followRow   = g_malloc0(in->numAnt*sizeof(ObitTableWXRow*));
  for (i=0; i<in->numAnt; i++) {
    in->priorRow[i]  = newObitTableWXRow(in->WXTable);
    in->followRow[i] = newObitTableWXRow(in->WXTable);
    /* Set time to large negative value to indicate uninitialized */
    in->priorRow[i]->Time    = -1.0e20;
    in->followRow[i]->Time   = -1.0e20;
    in->priorRowNo[i]        = -1;
    in->followRowNo[i]       = -1;
  }
  /* Open table */
  ObitTableWXOpen (in->WXTable, OBIT_IO_ReadOnly, err);
  if (err->error) Obit_traceback_msg (err, routine, in->WXTable->name);
  
  /* Loop through antennas finding first occurance */
  for (i=0; i<in->numAnt; i++) {
    /* Loop over table */
    ant = i + 1;
    found = FALSE;
    for (iRow=1; iRow<=in->WXTable->myDesc->nrow; iRow++) {
      ObitTableWXReadRow (in->WXTable, iRow, in->priorRow[i], err);
      if (err->error) Obit_traceback_msg (err, routine, in->WXTable->name);
      in->priorRowNo[i] = iRow;   /* Save row number */
      /* This the one? antenna 0 matches any */
      if (((in->priorRow[i]->antNo==ant) || (in->priorRow[i]->antNo==0)) && 
	  (in->priorRow[i]->SubA==in->SubA))
	{found=TRUE;break;}
    } /* end loop over table */
    /* Was it found? If not mark */
    if (!found) {
      in->priorRow[i]->Time    = -1.0e20;
      in->priorRowNo[i]        = -1;
    }
  } /* End loop finding first occurance */

  /* Loop through antennas finding second occurance */
  for (i=0; i<in->numAnt; i++) {
    /* Was first occurance found? */
    if (in->priorRowNo[i]>0) { 
      /* Loop over table */
      ant = i + 1;
      found = FALSE;
      start = MIN (in->priorRowNo[i]+1, in->WXTable->myDesc->nrow);
      for (iRow=start; iRow<=in->WXTable->myDesc->nrow; iRow++) {
	ObitTableWXReadRow (in->WXTable, iRow, in->followRow[i], err);
	if (err->error) Obit_traceback_msg (err, routine, in->WXTable->name);
	in->followRowNo[i] = iRow;   /* Save row number */
	/* This the one? antenna 0 matches any */
	if (((in->followRow[i]->antNo==ant) || (in->followRow[i]->antNo==0)) && 
	    (in->followRow[i]->SubA==in->SubA))
	  {found=TRUE;break;}
      } /* end loop over table */
    } /* End if first found */
    /* Was it found? If not mark */
    if (!found) {
      in->followRow[i]->Time    = -1.0e20;
      in->followRowNo[i]        = -1;
    }
  } /* End loop finding second occurance */
  
  /* Close table */
  ObitTableWXClose (in->WXTable, err);
  if (err->error) Obit_traceback_msg (err, routine, in->WXTable->name);
	
} /*  end InitSubarray */

/**
 * Update times of row data if time beyond follow.
 * For times after the last occurance of a given antenna/subarray
 * the followRow->Time=-1.0e-20
 * If no entries were found at all for an antenna, priorRow->Time=-1.0e-20
 * \param in   The object to update
 * \param time Desired time (day)
 * \param Ant  Desired (1-rel) Antenna number
 * \param err  Obit error stack object.
 */
static void UpdateTime (ObitWeather *in, ofloat time, olong Ant, 
			ObitErr *err)
{
  olong iRow, ant, start;
  gboolean found, OK, opened=FALSE;
  gchar *routine = "ObitWeather:UpdateTime";

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Update antenna if necessary */
  ant = Ant-1;
  if ((in->followRowNo[ant]>0) &&                         /* Follow set */
      (in->followRowNo[ant]<in->WXTable->myDesc->nrow) && /* before end */
      (time>in->followRow[ant]->Time)) {
    /* Open table */
    if (!opened)
      ObitTableWXOpen (in->WXTable, OBIT_IO_ReadOnly, err);
    if (err->error) Obit_traceback_msg (err, routine, in->WXTable->name);
    opened = TRUE;
    OK     = FALSE;
    /* loop until follow is past time or the end of the table */
    while(!OK) {
      /* Update - shuffle pointers */
      in->priorRow[ant]   = ObitTableWXRowUnref(in->priorRow[ant]);
      in->priorRow[ant]   = ObitTableWXRowRef(in->followRow[ant]);
      in->priorRowNo[ant] = in->followRowNo[ant];
      
      /* Loop over table looking for next occurance of ant+1 */
      start = MIN (in->priorRowNo[ant]+1, in->WXTable->myDesc->nrow);
      for (iRow=start; iRow<=in->WXTable->myDesc->nrow; iRow++) {
	ObitTableWXReadRow (in->WXTable, iRow, in->followRow[ant], err);
	if (err->error) Obit_traceback_msg (err, routine, in->WXTable->name);
	in->followRowNo[ant] = iRow;   /* Save row number */
	/* This the one? match antenna (0==any) and subarray */
	if (((in->followRow[ant]->antNo==ant) || (in->followRow[ant]->antNo==0)) && 
	    (in->followRow[ant]->SubA==in->SubA))
	  {found=TRUE;break;}
      } /* end loop over table */
	/* Find next? */
      if (!found) {
	in->followRow[ant]->Time = -1.0e20;
      }
      /* This one OK or past end */
      OK = ((time<=in->followRow[ant]->Time) || (iRow>=in->WXTable->myDesc->nrow));
    } /* end !OK loop */
  } /* end needs update */
  
  /* Close table */
  if (opened)
    ObitTableWXClose (in->WXTable, err);
  if (err->error) Obit_traceback_msg (err, routine, in->WXTable->name);
} /*  end UpdateTime */

