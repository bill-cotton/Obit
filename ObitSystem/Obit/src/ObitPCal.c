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

#include "ObitPCal.h"
#include "ObitTableANUtil.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitPCal.c
 * ObitPCal class function definitions.
 * This class is derived from the Obit base class.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitPCal";

/** Function to obtain parent ClassInfo */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo structure ObitPCalClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitPCalClassInfo myClassInfo = {FALSE};

/*--------------- File Global Variables  ----------------*/


/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitPCalInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitPCalClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitPCalClassInfoDefFn (gpointer inClass);

/** Private: Initialize for a subarray. */
static void InitSubarray (ObitPCal *in, olong subarray, ObitErr *err);

/** Private: Update time of data. */
static void UpdateTime (ObitPCal *in, ofloat time, olong Ant, 
			ObitErr *err);

/** Private: Create independent PC Table row */
static ObitTablePCRow* MakePCRow (ObitTablePC *PCTab);

/** Private: Delete independent PC Table row */
static ObitTablePCRow* KillPCRow (ObitTablePCRow *PCRow);

/** Private: Copy PC Table row */
static void CopyPCRow (ObitTablePCRow* PCTabIn, ObitTablePCRow* PCCTabOut, 
		       olong numIF, olong numBand);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitPCal* newObitPCal (gchar* name)
{
  ObitPCal* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitPCalClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitPCal));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitPCalInit((gpointer)out);

 return out;
} /* end newObitPCal */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitPCalGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitPCalClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitPCalGetClass */

/**
 * Make a deep copy of an ObitPCal.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitPCal* ObitPCalCopy  (ObitPCal *in, ObitPCal *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  gboolean oldExist;
  olong i;
  gchar *outName;
  gchar *routine = "ObitPCalCopy";

  /* error checks */
  if (err->error) return out;
  g_assert (ObitIsA(in, &myClassInfo));
  if (out) g_assert (ObitIsA(out, &myClassInfo));

  /* Create if it doesn't exist */
  oldExist = out!=NULL;
  if (!oldExist) {
    /* derive object name */
    outName = g_strconcat ("Copy: ",in->name,NULL);
    out = newObitPCal(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  /* Clear old structures */
  out->PCTable = ObitTablePCUnref(out->PCTable);
  if (out->antList) {
    for (i=0; i<out->numSubA; i++)
      out->antList[i]   = ObitAntennaListUnref(out->antList[i]);
    g_free(out->antList);
  }
  if (out->priorRowNo)  g_free(out->priorRowNo);
  if (out->followRowNo) g_free(out->followRowNo);
  if (out->priorRow){
    for (i=0; i<out->numAnt; i++) {
      KillPCRow(out->priorRow[i]);
    }
    g_free(out->priorRow);
  }
  if (out->followRow){
    for (i=0; i<out->numAnt; i++) {
      KillPCRow(out->followRow[i]);
    }
    g_free(out->followRow);
  }

  /* Create new  */
  out->PCTable  = ObitTablePCRef(in->PCTable);
  out->numSubA  = in->numSubA;
  out->antList  = g_malloc0(out->numSubA*sizeof(ObitAntennaList*));
  for (i=0; i<out->numSubA ; i++) {
    out->antList[i] = ObitAntennaListRef(in->antList[i]);
  }
  out->SubA        = in->SubA;
  out->numIF       = in->numIF;
  out->numTone     = in->numTone;
  out->numPoln     = in->numPoln;
  out->numAnt      = out->antList[out->SubA-1]->number;
  out->priorRowNo  = g_malloc0(out->numAnt*sizeof(olong));
  out->followRowNo = g_malloc0(out->numAnt*sizeof(olong));
  out->priorRow    = g_malloc0(out->numAnt*sizeof(ObitTablePCRow*));
  out->followRow   = g_malloc0(out->numAnt*sizeof(ObitTablePCRow*));
  for (i=0; i<out->numAnt; i++) {
    out->priorRow[i]  = MakePCRow(out->PCTable);
    out->followRow[i] = MakePCRow(out->PCTable);
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
} /* end ObitPCalCopy */

/**
 * Make a copy of a object but do not copy the actual data
 * This is useful to create an PCal similar to the input one.
 * \param in  The object to copy
 * \param out An existing object pointer for output, must be defined.
 * \param err Obit error stack object.
 */
void ObitPCalClone  (ObitPCal *in, ObitPCal *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  olong i;
  gchar *routine = "ObitPCalClone";

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
  out->PCTable = ObitTablePCUnref(out->PCTable);
  if (out->antList) {
    for (i=0; i<out->numSubA; i++)
      out->antList[i]   = ObitAntennaListUnref(out->antList[i]);
    g_free(out->antList);
  }
  if (out->priorRowNo)  g_free(out->priorRowNo);
  if (out->followRowNo) g_free(out->followRowNo);
  if (out->priorRow){
    for (i=0; i<out->numAnt; i++) {
      KillPCRow(out->priorRow[i]);
    }
    g_free(out->priorRow);
  }
  if (out->followRow){
    for (i=0; i<out->numAnt; i++) {
      KillPCRow(out->followRow[i]);
    }
    g_free(out->followRow);
  }

  /* Create new  */
  out->PCTable  = ObitTablePCRef(in->PCTable);
  out->numSubA  = in->numSubA;
  out->numIF    = in->numIF;
  out->numTone  = in->numTone;
  out->numPoln  = in->numPoln;
  out->antList  = g_malloc0(out->numSubA*sizeof(ObitAntennaList*));
  for (i=0; i<out->numSubA ; i++) {
    out->antList[i] = ObitAntennaListRef(in->antList[i]);
  }
  out->SubA        = in->SubA;
  out->numAnt      = out->antList[out->SubA-1]->number;
  out->priorRowNo  = g_malloc0(out->numAnt*sizeof(olong));
  out->followRowNo = g_malloc0(out->numAnt*sizeof(olong));
  out->priorRow    = g_malloc0(out->numAnt*sizeof(ObitTablePCRow*));
  out->followRow   = g_malloc0(out->numAnt*sizeof(ObitTablePCRow*));
  for (i=0; i<out->numAnt; i++) {
    out->priorRow[i]  = MakePCRow(out->PCTable);
    out->followRow[i] = MakePCRow(out->PCTable);
    /* Set time to large negative value to indicate uninitialized */
    out->priorRow[i]->Time    = -1.0e20;
    out->followRow[i]->Time   = -1.0e20;
    out->priorRowNo[i]        = -1;
    out->followRowNo[i]       = -1;
  }

  /* Init for Subarray out->SubA */
  InitSubarray(out, out->SubA, err);
  if (err->error) Obit_traceback_msg (err, routine, out->name);

} /* end ObitPCalClone */

/**
 * Creates an ObitPCal 
 * \param name    An optional name for the object.
 * \param PCTable PCal table to interpolate
 * \param UVData  UV data used to get antenna/array information
 * \param err     Obit error stack object.
 * \return the new object.
 */
ObitPCal* ObitPCalCreate (gchar* name, ObitTablePC *PCTable, 
			  ObitUV *UVData, ObitErr *err)
{
  ObitPCal* out;
  ObitTableAN  *ANTable=NULL;
  olong highVer, ver, iver;
  olong numOrb=0, numPCal=0;
  gchar *routine = "ObitPCalCreate";

  /* Create basic structure */
  out = newObitPCal (name);

  /* Save table */
  out->PCTable = ObitTablePCRef(PCTable);

  /* Get antenna/subarray information */
  /* Antenna lists - how name subarrays? */
  highVer = ObitTableListGetHigh (UVData->tableList, "AIPS AN");
  out->numSubA = highVer;
  out->antList = g_malloc0(out->numSubA*sizeof(ObitAntennaList*));

  /* How many IFs */
  if (UVData->myDesc->jlocif>=0)
    out->numIF = UVData->myDesc->inaxes[UVData->myDesc->jlocif];
  else
    out->numIF = 1;

  /* How many Polns */
  out->numPoln = UVData->myDesc->inaxes[UVData->myDesc->jlocs];
  
  /* How many Tones */
  out->numTone = PCTable->numTones;

  /* Loop over AN tables (subarrays) forming antenna lists*/
  for (iver=1; iver<=out->numSubA; iver++) {
    ver = iver;
    ANTable = newObitTableANValue ("AN table", (ObitData*)UVData, 
				   &ver, OBIT_IO_ReadOnly, out->numIF, numOrb, numPCal, err);
    if (ANTable==NULL) Obit_log_error(err, OBIT_Error, "ERROR with AN table");
    out->antList[iver-1] = ObitTableANGetList (ANTable, err);
    if (err->error) Obit_traceback_val (err, routine, UVData->name, out);
    
    ANTable = ObitTableANUnref(ANTable);   /* Cleanup */
  } /* End loop over subarrays */

  /* Initialize for subarray 1 */
  InitSubarray(out, 1, err);
  if (err->error) Obit_traceback_val (err, routine, UVData->name, out);

  return out;
} /* end ObitPCalCreate */

/**
 *  Interpolate PCal values at a given time
 *  Calls are expected in time order by antenna.
 * \param in      PCal interpolator
 * \param time    Desired time (day)
 * \param ant     Desired antenna number
 * \param suba    Desired subarray
 * \param CableCal [out] CABLE_CAL round() trip, seconds
 * \param Freq1    [out] Array of tone frequencies (Hz), (numTones, numBands) poln 1
 * \param PCal1    [out] Array of tone phases, (numTones, numBands) poln 1, 
 *                       may be fblanked, must be dimensioned at least numIF*numTone
 * \param Freq2    [out] Array of tone frequencies, poln 2
 * \param PCal2    [out] Array of tone phases, poln 2, may be fblanked, 
 * \param err      Obit error stack object.
 */
void ObitPCalReport (ObitPCal *in, ofloat time, olong ant, olong suba,
		     odouble *CableCal,
		     odouble *Freq1, ofloat *PCal1, 
		     odouble *Freq2, ofloat *PCal2, 
		     ObitErr *err)
{
  ofloat w1, w2, ww1, ww2, diff, dt1, dt2, ph1, ph2;
  ofloat fblank = ObitMagicF();
  olong i, iant = ant-1;
  gchar *routine = "ObitPCalReport";
  
  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  
  /* New subarray? */
  if (suba!=in->SubA) InitSubarray (in, suba, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Check/update times */
  UpdateTime (in, time, ant, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  dt1 = (in->priorRow[iant]->Time - time) * 86400.0;  /* In seconds */
  dt2 = (time - in->followRow[iant]->Time) * 86400.0;
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

  /* Cable cal */
  *CableCal = w1*in->priorRow[iant]->CableCal + w2*in->followRow[iant]->CableCal;

  /* Output data may be blanked - modify weights */
  /* First poln */
  for (i=0; i<in->numIF*in->numTone; i++) {
    /* Assume frequencies the same - otherwise this make no sense */
    Freq1[i] = in->priorRow[iant]->PCFreq1[i];

    ww1 = w1; ww2 = w2;
    if (in->priorRow[iant]->PCReal1[i]==fblank)  ww1 = 0.0;
    if (in->followRow[iant]->PCReal1[i]==fblank) ww2 = 0.0;
    if ((ww1+ww2)>0.0) {
      ph1 = atan2(in->priorRow[iant]->PCImag1[i], in->priorRow[iant]->PCReal1[i]);
      ph2 = atan2(in->followRow[iant]->PCImag1[i], in->followRow[iant]->PCReal1[i]);
      /* Add rate terms */
      ph1 += (dt1*in->priorRow[iant]->PCRate1[i])  * 2*G_PI * Freq1[i];
      ph2 -= (dt2*in->followRow[iant]->PCRate1[i]) * 2*G_PI * Freq1[i];
      PCal1[i] = ww1*ph1 + ww2*ph2;
    }
    else PCal1[i] = fblank;
  } /* end IF loop poln 1 */

  /* Second  poln if present */
  if (in->numPoln>1) {
    for (i=0; i<in->numIF*in->numTone; i++) {
     Freq2[i] = in->priorRow[iant]->PCFreq2[i];
     ww1 = w1; ww2 = w2;
      if (in->priorRow[iant]->PCReal2[i]==fblank)  ww1 = 0.0;
      if (in->followRow[iant]->PCReal2[i]==fblank) ww2 = 0.0;
      if ((ww1+ww2)>0.0) {
	ph1 = atan2(in->priorRow[iant]->PCImag2[i], in->priorRow[iant]->PCReal2[i]);
	ph2 = atan2(in->followRow[iant]->PCImag2[i], in->followRow[iant]->PCReal2[i]);
	/* Add rate terms */
	ph1 += (dt1*in->priorRow[iant]->PCRate2[i])  * 2*G_PI * Freq1[i];
	ph2 -= (dt2*in->followRow[iant]->PCRate2[i]) * 2*G_PI * Freq1[i];
	PCal2[i] = ww1*ph1 + ww2*ph2;
    }
      else PCal2[i] = fblank;
    } /* end IF loop poln 1 */
  }  /* end seconf poln */

} /* end ObitPCalReport */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitPCalClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitPCalClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitPCalClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitPCalClassInfoDefFn (gpointer inClass)
{
  ObitPCalClassInfo *theClass = (ObitPCalClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitPCalClassInit;
  theClass->newObit       = (newObitFP)newObitPCal;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitPCalClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitPCalGetClass;
  theClass->ObitCopy      = (ObitCopyFP)ObitPCalCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitPCalClear;
  theClass->ObitInit      = (ObitInitFP)ObitPCalInit;
  theClass->ObitPCalCreate = (ObitPCalCreateFP)ObitPCalCreate;

} /* end ObitPCalClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitPCalInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitPCal *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->PCTable     = NULL;
  in->PCRow       = NULL;
  in->numSubA     =  0;
  in->antList     = NULL;
  in->SubA        = 0;
  in->numAnt      = 0;
  in->priorRowNo  = NULL;
  in->followRowNo = NULL;
  in->priorRow    = NULL;
  in->followRow   = NULL;

} /* end ObitPCalInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *             Actually it should be an ObitPCal* cast to an Obit*.
 */
void ObitPCalClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitPCal *in = inn;
  olong i;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  in->PCTable = ObitTablePCUnref(in->PCTable);
  if (in->antList) {
    for (i=0; i<in->numSubA; i++)
      in->antList[i]   = ObitAntennaListUnref(in->antList[i]);
    g_free(in->antList);
  }
  if (in->priorRowNo)  g_free(in->priorRowNo);
  if (in->followRowNo) g_free(in->followRowNo);
  if (in->priorRow){
    for (i=0; i<in->numAnt; i++) {
      KillPCRow(in->priorRow[i]);
    }
    g_free(in->priorRow);
  }
  if (in->followRow){
    for (i=0; i<in->numAnt; i++) {
      KillPCRow(in->followRow[i]);
    }
    g_free(in->followRow);
  }

  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitPCalClear */

/**
 * Initialize object for a given subarray
 * Created row structures and reads prior and follow for each antenna.
 * \param in   The object to init
 * \param SubA Which (1-rel) subarray?
 * \param err  Obit error stack object.
 */
static void InitSubarray (ObitPCal *in, olong SubA, ObitErr *err)
{
  olong i, iRow, ant, start;
  gboolean found;
  gchar *routine = "ObitPCal:InitSubarray";

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Clear old structures */
  if (in->priorRowNo)  g_free(in->priorRowNo);
  if (in->followRowNo) g_free(in->followRowNo);
  if (in->priorRow){
    for (i=0; i<in->numAnt; i++) {
      KillPCRow(in->priorRow[i]);
    }
    g_free(in->priorRow);
  }
  if (in->followRow){
    for (i=0; i<in->numAnt; i++) {
      KillPCRow(in->followRow[i]);
    }
    g_free(in->followRow);
  }

  /* Create new for this subarray */
  in->SubA        = SubA;
  in->numAnt      = in->antList[SubA-1]->number;
  in->priorRowNo  = g_malloc0(in->numAnt*sizeof(olong));
  in->followRowNo = g_malloc0(in->numAnt*sizeof(olong));
  in->priorRow    = g_malloc0(in->numAnt*sizeof(ObitTablePCRow*));
  in->followRow   = g_malloc0(in->numAnt*sizeof(ObitTablePCRow*));
  for (i=0; i<in->numAnt; i++) {
    in->priorRow[i]  = MakePCRow(in->PCTable);
    in->followRow[i] = MakePCRow(in->PCTable);
    /* Set time to large negative value to indicate uninitialized */
    in->priorRow[i]->Time    = -1.0e20;
    in->followRow[i]->Time   = -1.0e20;
    in->priorRowNo[i]        = -1;
    in->followRowNo[i]       = -1;
  }
  if (in->PCRow==NULL) in->PCRow = MakePCRow(in->PCTable);

  /* Open table */
  ObitTablePCOpen (in->PCTable, OBIT_IO_ReadOnly, err);
  if (err->error) Obit_traceback_msg (err, routine, in->PCTable->name);
  
  /* Loop through antennas finding first occurance */
  for (i=0; i<in->numAnt; i++) {
    /* Loop over table */
    ant = i + 1;
    found = FALSE;
    for (iRow=1; iRow<=in->PCTable->myDesc->nrow; iRow++) {
      ObitTablePCReadRow (in->PCTable, iRow, in->PCRow, err);
      CopyPCRow (in->PCRow, in->priorRow[i], in->numIF, in->numTone);  /* Save values */
      if (err->error) Obit_traceback_msg (err, routine, in->PCTable->name);
      in->priorRowNo[i] = iRow;   /* Save row number */
      /* This the one? */
      if ((in->priorRow[i]->antennaNo==ant) && (in->priorRow[i]->Array==in->SubA)) 
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
      start = MIN (in->priorRowNo[i]+1, in->PCTable->myDesc->nrow);
      for (iRow=start; iRow<=in->PCTable->myDesc->nrow; iRow++) {
	ObitTablePCReadRow (in->PCTable, iRow, in->PCRow, err);
	CopyPCRow (in->PCRow, in->followRow[i], in->numIF, in->numTone);  /* Save values */
	if (err->error) Obit_traceback_msg (err, routine, in->PCTable->name);
	in->followRowNo[i] = iRow;   /* Save row number */
	/* This the one? */
	if ((in->followRow[i]->antennaNo==ant) && (in->followRow[i]->Array==in->SubA)) 
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
  ObitTablePCClose (in->PCTable, err);
  if (err->error) Obit_traceback_msg (err, routine, in->PCTable->name);
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
static void UpdateTime (ObitPCal *in, ofloat time, olong Ant, 
			ObitErr *err)
{
  olong i, iRow, ant, start;
  gboolean found, OK, opened=FALSE;
  gchar *routine = "ObitPCal:UpdateTime";

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Update antenna if necessary */
  ant = Ant-1;
  if ((in->followRowNo[ant]>0) &&                         /* Follow set */
      (in->followRowNo[ant]<in->PCTable->myDesc->nrow) && /* before end */
      (time>in->followRow[ant]->Time)) {
    /* Open table */
    if (!opened)
      ObitTablePCOpen (in->PCTable, OBIT_IO_ReadOnly, err);
    if (err->error) Obit_traceback_msg (err, routine, in->PCTable->name);
    opened = TRUE;
    OK     = FALSE;
    /* loop until follow is past time or the end of the table */
    while(!OK) {
      /* Update - shuffle data */
      CopyPCRow(in->followRow[ant], in->priorRow[ant], in->numIF, in->numTone);
      in->priorRowNo[ant] = in->followRowNo[ant];
      
      /* Loop over table looking for next occurance of ant+1 */
      start = MIN (in->priorRowNo[ant]+1, in->PCTable->myDesc->nrow);
      for (iRow=start; iRow<=in->PCTable->myDesc->nrow; iRow++) {
	ObitTablePCReadRow (in->PCTable, iRow, in->PCRow, err);
	CopyPCRow (in->PCRow, in->followRow[ant], in->numIF, in->numTone);  /* Save values */
	if (err->error) Obit_traceback_msg (err, routine, in->PCTable->name);
	in->followRowNo[ant] = iRow;   /* Save row number */
	/* This the one? match antenna and subarray */
	if ((in->followRow[ant]->antennaNo==(ant+1)) && (in->followRow[ant]->Array==in->SubA)) 
	  {found=TRUE;break;}
      } /* end loop over table */
	/* Find next? */
      if (!found) {
	in->followRow[ant]->Time = -1.0e20;
      }
      /* This one OK or past end */
      OK = ((time<=in->followRow[ant]->Time) || (iRow>=in->PCTable->myDesc->nrow));
    } /* end !OK loop */
  } /* end needs update */
  
  /* Close table */
  if (opened)
    ObitTablePCClose (in->PCTable, err);
  if (err->error) Obit_traceback_msg (err, routine, in->PCTable->name);
} /*  end UpdateTime */

/**
 * Independent TablePC Row Constructor.
 * CANNOT be used for actual IO, only for storage
 * \param PCTab  Table for new structure
 * \return the new Table PC Row object.
 */
static ObitTablePCRow* MakePCRow (ObitTablePC *PCTab) 
{
  ObitTablePCRow *out=NULL;
  olong numIF, size;

  numIF = PCTab->numBand;
  size  = numIF*PCTab->numTones;
  out   = newObitTablePCRow(PCTab);

  /* Replace pointers to buffer with real arrays */
  out->State1  = g_malloc0(numIF*4*sizeof(ofloat));
  out->PCFreq1 = g_malloc0(size*sizeof(odouble));
  out->PCReal1 = g_malloc0(size*sizeof(ofloat));
  out->PCImag1 = g_malloc0(size*sizeof(ofloat));
  out->PCRate1 = g_malloc0(size*sizeof(ofloat));
  if (PCTab->numPol>1) {
    out->State2  = g_malloc0(numIF*4*sizeof(ofloat));
    out->PCFreq2 = g_malloc0(size*sizeof(odouble));
    out->PCReal2 = g_malloc0(size*sizeof(ofloat));
    out->PCImag2 = g_malloc0(size*sizeof(ofloat));
    out->PCRate2 = g_malloc0(size*sizeof(ofloat));
  }
  return out;
} /* end MakePCRow  */

/**
 * Independent TablePC Row Destructor
 * \param PCRow  Row structure to delete;
 * \return NULL
 */
static ObitTablePCRow* KillPCRow (ObitTablePCRow *PCRow) 
{
  if (PCRow==NULL) return NULL;
  if (PCRow->State1)  g_free(PCRow->State1);
  if (PCRow->PCFreq1) g_free(PCRow->PCFreq1);
  if (PCRow->PCReal1) g_free(PCRow->PCReal1);
  if (PCRow->PCImag1) g_free(PCRow->PCImag1);
  if (PCRow->PCRate1) g_free(PCRow->PCRate1);
  if (PCRow->State2)  g_free(PCRow->State2);
  if (PCRow->PCFreq2) g_free(PCRow->PCFreq2);
  if (PCRow->PCReal2) g_free(PCRow->PCReal2);
  if (PCRow->PCImag2) g_free(PCRow->PCImag2);
  if (PCRow->PCRate2) g_free(PCRow->PCRate2);
  return ObitTablePCRowUnref(PCRow);
} /* end  KillPCRow */

/**
 * Copy TablePC Row 
 * \param PCRowIn  Input Row structure
 * \param PCRowOut Output Row structure
 * \param numIF    Number of bands
 * \peram numTone  Number of tones per band
 */
static void CopyPCRow (ObitTablePCRow* PCRowIn, ObitTablePCRow* PCRowOut, 
		       olong numIF, olong numTone)
{
  olong i, numCopy;

  PCRowOut->Time      = PCRowIn->Time;
  PCRowOut->TimeI     = PCRowIn->TimeI;
  PCRowOut->SourID    = PCRowIn->SourID;
  PCRowOut->antennaNo = PCRowIn->antennaNo;
  PCRowOut->Array     = PCRowIn->Array;
  PCRowOut->FreqID    = PCRowIn->FreqID;
  PCRowOut->CableCal  = PCRowIn->CableCal;
  for (i=0; i<numIF*4; i++) PCRowOut->State1[i] = PCRowIn->State1[i];
  numCopy = numIF * numTone;
  for (i=0; i<numCopy; i++) PCRowOut->PCFreq1[i] = PCRowIn->PCFreq1[i];
  for (i=0; i<numCopy; i++) PCRowOut->PCReal1[i] = PCRowIn->PCReal1[i];
  for (i=0; i<numCopy; i++) PCRowOut->PCImag1[i] = PCRowIn->PCImag1[i];
  for (i=0; i<numCopy; i++) PCRowOut->PCRate1[i] = PCRowIn->PCRate1[i];
  if (PCRowIn->State2) for (i=0; i<numIF*4; i++)  PCRowOut->State2[i]  = PCRowIn->State2[i];
  if (PCRowIn->PCFreq2) for (i=0; i<numCopy; i++) PCRowOut->PCFreq2[i] = PCRowIn->PCFreq2[i];
  if (PCRowIn->PCReal2) for (i=0; i<numCopy; i++) PCRowOut->PCReal2[i] = PCRowIn->PCReal2[i];
  if (PCRowIn->PCImag2) for (i=0; i<numCopy; i++) PCRowOut->PCImag2[i] = PCRowIn->PCImag2[i];
  if (PCRowIn->PCRate2) for (i=0; i<numCopy; i++) PCRowOut->PCRate2[i] = PCRowIn->PCRate2[i];
} /* end  CopyPCRow */

