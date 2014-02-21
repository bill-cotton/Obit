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

#include "ObitSwPower.h"
#include "ObitTableANUtil.h"
#include "ObitTableUtil.h"
#include "ObitUVSoln.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitSwPower.c
 * ObitSwPower class function definitions.
 * This class is derived from the Obit base class.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitSwPower";

/** Function to obtain parent ClassInfo */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo structure ObitSwPowerClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitSwPowerClassInfo myClassInfo = {FALSE};

/*--------------- File Global Variables  ----------------*/


/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitSwPowerInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitSwPowerClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitSwPowerClassInfoDefFn (gpointer inClass);

/** Private: Initialize for a subarray. */
static void InitSubarray (ObitSwPower *in, olong subarray, ObitErr *err);

/** Private: Update time of data. */
static void UpdateTime (ObitSwPower *in, ofloat time, olong Ant, 
			ObitErr *err);

/** Private: Update time of data. */
static void UpdateTime (ObitSwPower *in, ofloat time, olong Ant, 
			ObitErr *err);

/** Private: Create independent SY Table row */
static ObitTableSYRow* MakeSYRow (ObitTableSY *SYTab);

/** Private: Delete independent SY Table row */
static ObitTableSYRow* KillSYRow (ObitTableSYRow *SYRow);

/** Private: Copy SY Table row */
static void CopySYRow (ObitTableSYRow* SYTabIn, ObitTableSYRow* SyTabOut, 
		       olong numIF);

/** Private:  Determine number of times */
static olong SYCountTime (ObitTableSY *SYTab, olong isub, ObitErr* err);

/** Private:  Smooth SY table */
void 
SYSmooth (ObitTableSY *SYTab, gchar* smoFunc, gchar* smoType, ofloat alpha, 
	  ofloat *smoParm, olong iif, olong sub, 
	  olong nxt, ofloat* work1, ofloat* work2, gboolean doBlank, ObitErr* err);

/** Private: Generic smoothing */
void SYsmoIt (gchar* smmeth, ofloat width, ofloat alpha, 
	      ofloat* x, ofloat *t, ofloat *w, olong n, 
	      ofloat* xs, ofloat* ws, ofloat *wrk1, ofloat *wrk2, 
	      gboolean doBlank);
/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitSwPower* newObitSwPower (gchar* name)
{
  ObitSwPower* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitSwPowerClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitSwPower));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitSwPowerInit((gpointer)out);

 return out;
} /* end newObitSwPower */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitSwPowerGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitSwPowerClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitSwPowerGetClass */

/**
 * Make a deep copy of an ObitSwPower.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitSwPower* ObitSwPowerCopy  (ObitSwPower *in, ObitSwPower *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  gboolean oldExist;
  olong i;
  gchar *outName;
  gchar *routine = "ObitSwPowerCopy";

  /* error checks */
  if (err->error) return out;
  g_assert (ObitIsA(in, &myClassInfo));
  if (out) g_assert (ObitIsA(out, &myClassInfo));

  /* Create if it doesn't exist */
  oldExist = out!=NULL;
  if (!oldExist) {
    /* derive object name */
    outName = g_strconcat ("Copy: ",in->name,NULL);
    out = newObitSwPower(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  /* Clear old structures */
  out->SYTable = ObitTableSYUnref(out->SYTable);
  if (out->antList) {
    for (i=0; i<out->numSubA; i++)
      out->antList[i]   = ObitAntennaListUnref(out->antList[i]);
    g_free(out->antList);
  }
  if (out->priorRowNo)  g_free(out->priorRowNo);
  if (out->followRowNo) g_free(out->followRowNo);
  if (out->priorRow){
    for (i=0; i<out->numAnt; i++) {
      KillSYRow(out->priorRow[i]);
    }
    g_free(out->priorRow);
  }
  if (out->followRow){
    for (i=0; i<out->numAnt; i++) {
      KillSYRow(out->followRow[i]);
    }
    g_free(out->followRow);
  }

  /* Create new  */
  out->SYTable  = ObitTableSYRef(in->SYTable);
  out->numSubA  = in->numSubA;
  out->antList  = g_malloc0(out->numSubA*sizeof(ObitAntennaList*));
  for (i=0; i<out->numSubA ; i++) {
    out->antList[i] = ObitAntennaListRef(in->antList[i]);
  }
  out->SubA        = in->SubA;
  out->numIF       = in->numIF;
  out->numPoln     = in->numPoln;
  out->numAnt      = out->antList[out->SubA-1]->number;
  out->priorRowNo  = g_malloc0(out->numAnt*sizeof(olong));
  out->followRowNo = g_malloc0(out->numAnt*sizeof(olong));
  out->priorRow    = g_malloc0(out->numAnt*sizeof(ObitTableSYRow*));
  out->followRow   = g_malloc0(out->numAnt*sizeof(ObitTableSYRow*));
  for (i=0; i<out->numAnt; i++) {
    out->priorRow[i]  = newObitTableSYRow(out->SYTable);
    out->followRow[i] = newObitTableSYRow(out->SYTable);
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
} /* end ObitSwPowerCopy */

/**
 * Make a copy of a object but do not copy the actual data
 * This is useful to create an SwPower similar to the input one.
 * \param in  The object to copy
 * \param out An existing object pointer for output, must be defined.
 * \param err Obit error stack object.
 */
void ObitSwPowerClone  (ObitSwPower *in, ObitSwPower *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  olong i;
  gchar *routine = "ObitSwPowerClone";

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
  out->SYTable = ObitTableSYUnref(out->SYTable);
  if (out->antList) {
    for (i=0; i<out->numSubA; i++)
      out->antList[i]   = ObitAntennaListUnref(out->antList[i]);
    g_free(out->antList);
  }
  if (out->priorRowNo)  g_free(out->priorRowNo);
  if (out->followRowNo) g_free(out->followRowNo);
  if (out->priorRow){
    for (i=0; i<out->numAnt; i++) {
      KillSYRow(out->priorRow[i]);
    }
    g_free(out->priorRow);
  }
  if (out->followRow){
    for (i=0; i<out->numAnt; i++) {
      KillSYRow(out->followRow[i]);
    }
    g_free(out->followRow);
  }

  /* Create new  */
  out->SYTable  = ObitTableSYRef(in->SYTable);
  out->numSubA  = in->numSubA;
  out->numIF    = in->numIF;
  out->numPoln  = in->numPoln;
  out->antList  = g_malloc0(out->numSubA*sizeof(ObitAntennaList*));
  for (i=0; i<out->numSubA ; i++) {
    out->antList[i] = ObitAntennaListRef(in->antList[i]);
  }
  out->SubA        = in->SubA;
  out->numAnt      = out->antList[out->SubA-1]->number;
  out->priorRowNo  = g_malloc0(out->numAnt*sizeof(olong));
  out->followRowNo = g_malloc0(out->numAnt*sizeof(olong));
  out->priorRow    = g_malloc0(out->numAnt*sizeof(ObitTableSYRow*));
  out->followRow   = g_malloc0(out->numAnt*sizeof(ObitTableSYRow*));
  for (i=0; i<out->numAnt; i++) {
    out->priorRow[i]  = newObitTableSYRow(out->SYTable);
    out->followRow[i] = newObitTableSYRow(out->SYTable);
    /* Set time to large negative value to indicate uninitialized */
    out->priorRow[i]->Time    = -1.0e20;
    out->followRow[i]->Time   = -1.0e20;
    out->priorRowNo[i]        = -1;
    out->followRowNo[i]       = -1;
  }

  /* Init for Subarray out->SubA */
  InitSubarray(out, out->SubA, err);
  if (err->error) Obit_traceback_msg (err, routine, out->name);

} /* end ObitSwPowerClone */

/**
 * Creates an ObitSwPower 
 * If info item "calInt" on UVData exists and is > 0, then SYTable
 * is copied to a new table and smoothed to calInt (MWF, alpha=0.5)
 * \param name    An optional name for the object.
 * \param SYTable SwPower table to interpolate
 * \param UVData  UV data used to get antenna/array information, info items:
 * \li calInt   OBIT_float (1) Calibration interval in sec [def 30 sec]
 * \li doSmoo   OBIT_boo   (1) Smooth SY table to calInt?  [def TRUE]
 * \param err     Obit error stack object.
 * \return the new object.
 */
ObitSwPower* ObitSwPowerCreate (gchar* name, ObitTableSY *SYTable, 
				ObitUV *UVData, ObitErr *err)
{
  ObitSwPower* out;
  ObitTableAN  *ANTable=NULL;
  gint32 dim[MAXINFOELEMDIM];
  ObitInfoType type;
  olong highVer, ver, iver, isuba;
  ofloat calInt=30.0, smoParm[3];
  olong numOrb=0, numPCal=0;
  gboolean doSmoo=TRUE;
  gchar *routine = "ObitSwPowerCreate";

  /* Create basic structure */
  out = newObitSwPower (name);

  /* Copy/smooth input? */
  ObitInfoListGetTest(UVData->info, "calInt",  &type, dim, &calInt);
  ObitInfoListGetTest(UVData->info, "doSmoo",  &type, dim, &doSmoo);
  if ((calInt>0.0) && doSmoo) {
    /* Copy/smooth */
    highVer = ObitTableListGetHigh (UVData->tableList, "AIPS SY");
    ver     = highVer+1;
    out->SYTable = newObitTableSYValue("SYTable", (ObitData*)UVData, &ver, OBIT_IO_ReadWrite,  
				       SYTable->nIF, SYTable->nPol, err);
    out->SYTable = ObitTableSYCopy (SYTable, out->SYTable, err);
    if (err->error) Obit_traceback_val (err, routine, UVData->name, out);
 
    /* Tell */
    Obit_log_error(err, OBIT_InfoErr, "Smooth AIPS SY %d to %d, smooth time=%f", 
		   SYTable->tabVer, out->SYTable->tabVer, calInt);
    
    
    dim[0] = 3; dim[1] = dim[2] = 1;
    smoParm[0] = calInt/3600.0; smoParm[1] = smoParm[2] = 0.0;
    ObitInfoListAlwaysPut(out->SYTable->info, "smoParm", OBIT_float, dim, &smoParm);
    isuba = 1;
    ObitSwPowerSYSmo (out->SYTable, isuba, err);
    if (err->error) Obit_traceback_val (err, routine, UVData->name, out);
  } else {
    /* Save table */
    out->SYTable = ObitTableSYRef(SYTable);
  }

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
} /* end ObitSwPowerCreate */

/**
 *  Interpolate SwPower values at a given time
 *  Calls are expected in time order by antenna.
 * \param in      SwPower interpolator
 * \param time    Desired time (day)
 * \param ant     Desired antenna number
 * \param suba    Desired subarray
 * \param PwrDif1 [out] (P_on-P_off)*G  Poln # 1, 1 per IF
 *                      must be dimensioned at least numIF
 * \param PwrSum1 [out] (P_on+P_off)*G  Poln # 1, 1 per IF
 * \param Gain1   [out] Post switched power gain Poln # 1, 1 per IF
 * \param PwrDif2 [out] (P_on-P_off)*G  Poln # 1, 1 per IF
 * \param PwrSum2 [out] (P_on+P_off)*G  Poln # 1, 1 per IF
 * \param Gain2   [out] Post switched power gain Poln # 1, 1 per IF
 * \param err     Obit error stack object.
 */
void ObitSwPowerReport (ObitSwPower *in, ofloat time, olong ant, olong suba,
			ofloat *PwrDif1, ofloat *PwrSum1, ofloat *Gain1,
			ofloat *PwrDif2, ofloat *PwrSum2, ofloat *Gain2,
			ObitErr *err)
{
  ofloat w1, w2, ww1, ww2, diff;
  ofloat fblank = ObitMagicF();
  olong i, iant = ant-1;
  gchar *routine = "ObitSwPowerReport";
  
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
  /* Difference */
  /* First poln */
  for (i=0; i<in->numIF; i++) {
    ww1 = w1; ww2 = w2;
    if (in->priorRow[iant]->PwrDif1[i]==fblank)  ww1 = 0.0;
    if (in->followRow[iant]->PwrDif1[i]==fblank) ww2 = 0.0;
    if ((ww1+ww2)>0.0) PwrDif1[i] = ww1*in->priorRow[iant]->PwrDif1[i] + 
			 ww2*in->followRow[iant]->PwrDif1[i];
    else PwrDif1[i] = fblank;
  } /* end IF loop poln 1 */

  /* Second  poln if present */
  if (in->numPoln>1) {
    for (i=0; i<in->numIF; i++) {
      ww1 = w1; ww2 = w2;
      if (in->priorRow[iant]->PwrDif2[i]==fblank)  ww1 = 0.0;
      if (in->followRow[iant]->PwrDif2[i]==fblank) ww2 = 0.0;
      if ((ww1+ww2)>0.0) PwrDif2[i] = ww1*in->priorRow[iant]->PwrDif2[i] + 
			   ww2*in->followRow[iant]->PwrDif2[i];
      else PwrDif2[i] = fblank;
    } /* end IF loop poln 1 */
  }  /* end second poln */

  /* Sum */
  /* First poln */
  for (i=0; i<in->numIF; i++) {
    ww1 = w1; ww2 = w2;
    if (in->priorRow[iant]->PwrSum1[i]==fblank)  ww1 = 0.0;
    if (in->followRow[iant]->PwrSum1[i]==fblank) ww2 = 0.0;
    if ((ww1+ww2)>0.0) PwrSum1[i] = ww1*in->priorRow[iant]->PwrSum1[i] + 
			 ww2*in->followRow[iant]->PwrSum1[i];
    else PwrSum1[i] = fblank;
  } /* end IF loop poln 1 */

  /* Second  poln if present */
  if (in->numPoln>1) {
    for (i=0; i<in->numIF; i++) {
      ww1 = w1; ww2 = w2;
      if (in->priorRow[iant]->PwrSum2[i]==fblank)  ww1 = 0.0;
      if (in->followRow[iant]->PwrSum2[i]==fblank) ww2 = 0.0;
      if ((ww1+ww2)>0.0) PwrSum2[i] = ww1*in->priorRow[iant]->PwrSum2[i] + 
			   ww2*in->followRow[iant]->PwrSum2[i];
      else PwrSum2[i] = fblank;
    } /* end IF loop poln 1 */
  }  /* end second poln */

  /* Gain */
  /* First poln */
  for (i=0; i<in->numIF; i++) {
    ww1 = w1; ww2 = w2;
    if (in->priorRow[iant]->Gain1[i]==fblank)  ww1 = 0.0;
    if (in->followRow[iant]->Gain1[i]==fblank) ww2 = 0.0;
    if ((ww1+ww2)>0.0) Gain1[i] = ww1*in->priorRow[iant]->Gain1[i] + 
			 ww2*in->followRow[iant]->Gain1[i];
    else Gain1[i] = fblank;
  } /* end IF loop poln 1 */

  /* Second  poln if present */
  if (in->numPoln>1) {
    for (i=0; i<in->numIF; i++) {
      ww1 = w1; ww2 = w2;
      if (in->priorRow[iant]->Gain2[i]==fblank)  ww1 = 0.0;
      if (in->followRow[iant]->Gain2[i]==fblank) ww2 = 0.0;
      if ((ww1+ww2)>0.0) Gain2[i] = ww1*in->priorRow[iant]->Gain2[i] + 
			   ww2*in->followRow[iant]->Gain2[i];
      else Gain2[i] = fblank;
    } /* end IF loop poln 1 */
  }  /* end second poln */

} /* end ObitSwPowerReport */

/**
 * Smooths an SY Table in time.
 * Flagged values are optionally interpolated.
 * Controls on SYTab:
 * \li smoFunc   OBIT_string (4,1,1) smoothing function 
 *               'MWF' (median window filter), 
 *               "GAUS' (Gaussian) else "BOX", [def "MWF"]
 * \li alpha     OBIT_float (1,1,1) Alpha factor for MWF (0->box, 1->pure MWF) 
 *               [def 0.5]
 * \li smoParm   OBIT_float (3,1,1) Amplitude smoothing time in hr. [def 0.0]
 *               PwrDif, PwrSum, Gain
 * \li smoType   OBIT_string (4,1,1) Data to be smoothed
 *               "DIF ", "ALL ", "    " = "ALL " [def "DIF "]
 * \li doBlank   OBIT_bool (1,1,1) Replace blanked values with interpolated? [def true]
 * \param SYTab  SY table object 
 * \param isuba  Desired subarray, 0=> 1 
 * \param err    Error/message stack, returns if error.
 */
void ObitSwPowerSYSmo (ObitTableSY *SYTab, olong isuba, ObitErr* err) 
{
  ObitIOCode retCode;
  gint32 dim[MAXINFOELEMDIM];
  ObitInfoType type;
  gchar  smtype[5], smfunc[5];
  olong   i, iif, isub, numant, numpol, numif;
  gboolean doBlank;
  ofloat alpha, smparm[5];
  olong   mxtime, ntmp;
  ofloat *work1=NULL, *work2=NULL;
  gchar *routine = "ObitSwPowerSYSmo";

  /* Error checks */
  if (err->error) return;  /* previous error? */
  g_assert(ObitTableSYIsA(SYTab));

  /* Control Info */
  strcpy (smfunc, "MWF ");dim[0] = strlen(smfunc);
  ObitInfoListGetTest(SYTab->info, "smoFunc", &type, dim, smfunc);
  smfunc[dim[0]] = 0;
  strcpy (smtype, "DIF ");dim[0] = strlen(smtype);
  ObitInfoListGetTest(SYTab->info, "smoType", &type, dim, smtype);
  smtype[dim[0]] = 0;
  if (!strncmp (smtype, "    ",4)) strcpy (smtype, "ALL ");
  alpha = 0.5;
  ObitInfoListGetTest(SYTab->info, "alpha", &type, dim, &alpha);

  /* Smoothing times */
  for (i=0; i<5; i++) smparm[i] = 0.0;
  ObitInfoListGetTest(SYTab->info, "smoParm", &type, dim, smparm);
  for (i=0; i<5; i++) smparm[i] /= 24.0;  /* To days */
  doBlank = FALSE;
  ObitInfoListGetTest(SYTab->info, "doBlank", &type, dim, &doBlank);

  /* Subarray */
  isub = MAX (1, isuba);
  
  /* Sort to time-antenna order - just to be sure */
  ObitTableUtilSort2f ((ObitTable*)SYTab, "TIME    ", 1, FALSE, "ANTENNA NO.", 
		       1, FALSE, err);
  if (err->error) Obit_traceback_msg (err, routine, SYTab->name);

  /* Count number of times appearing in table */
  mxtime = SYCountTime(SYTab, isuba, err);
  if (err->error) Obit_traceback_msg (err, routine, SYTab->name);
  mxtime += 100;  /* Fudge a bit on the number of times */

  /* Sort to antenna-time order */
  ObitTableUtilSort2f ((ObitTable*)SYTab, "ANTENNA NO.", 1, FALSE, "TIME  ", 
		       1, FALSE, err);
  if (err->error) Obit_traceback_msg (err, routine, SYTab->name);

  /* Open table */
  retCode = ObitTableSYOpen (SYTab, OBIT_IO_ReadWrite, err);
  if (err->error) Obit_traceback_msg (err, routine, SYTab->name);
 
  /* Get descriptive info */
  numif  = SYTab->nIF;
  numpol = SYTab->nPol;
  numant = SYTab->nAnt;

  /* Primary work array */
  work1 = g_malloc(10*mxtime*sizeof(ofloat));
  
  /* How many secondary work arrays */
  if (!strncmp(smfunc, "GAUS", 4)) ntmp = 3;
  else if (!strncmp(smfunc, "MWF", 3)) ntmp = 4;
  else ntmp = 2;
  work2 = g_malloc(ntmp*mxtime*sizeof(ofloat));
  
  for (iif=0; iif<numif; iif++) { /* loop over IFs 300 */
    /* Smooth this IF */
    SYSmooth (SYTab, smfunc, smtype, alpha, smparm, iif, isub, 
	      mxtime, work1, work2, doBlank, err);
    if (err->error) goto cleanup;
  } /* end loop  L300: */
  if (work1) g_free(work1); work1 = NULL;
  if (work2) g_free(work2); work2 = NULL;

  /* Return to time-antenna order */
  ObitTableUtilSort2f ((ObitTable*)SYTab, "TIME  ", 1, FALSE, "ANTENNA NO.", 
		       1, FALSE, err);
  if (err->error) Obit_traceback_msg (err, routine, SYTab->name);

  /* Close output table */
cleanup: 
  retCode = ObitTableSYClose (SYTab, err);
  
  /* Cleanup */
  if (work1) g_free(work1);
  if (work2) g_free(work2);
  if (err->error) Obit_traceback_msg (err, routine, SYTab->name);
} /* end of routine ObitSwPowerSYSmo */ 

/**
 * Initialize global ClassInfo Structure.
 */
void ObitSwPowerClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitSwPowerClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitSwPowerClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitSwPowerClassInfoDefFn (gpointer inClass)
{
  ObitSwPowerClassInfo *theClass = (ObitSwPowerClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitSwPowerClassInit;
  theClass->newObit       = (newObitFP)newObitSwPower;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitSwPowerClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitSwPowerGetClass;
  theClass->ObitCopy      = (ObitCopyFP)ObitSwPowerCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitSwPowerClear;
  theClass->ObitInit      = (ObitInitFP)ObitSwPowerInit;
  theClass->ObitSwPowerCreate = (ObitSwPowerCreateFP)ObitSwPowerCreate;

} /* end ObitSwPowerClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitSwPowerInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitSwPower *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->SYTable     = NULL;
  in->SYRow       = NULL;
  in->numSubA     =  0;
  in->antList     = NULL;
  in->SubA        = 0;
  in->numAnt      = 0;
  in->priorRowNo  = NULL;
  in->followRowNo = NULL;
  in->priorRow    = NULL;
  in->followRow   = NULL;

} /* end ObitSwPowerInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *             Actually it should be an ObitSwPower* cast to an Obit*.
 */
void ObitSwPowerClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitSwPower *in = inn;
  olong i;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  in->SYTable = ObitTableSYUnref(in->SYTable);
  in->SYRow   = ObitTableSYRowUnref(in->SYRow);
  if (in->antList) {
    for (i=0; i<in->numSubA; i++)
      in->antList[i]   = ObitAntennaListUnref(in->antList[i]);
    g_free(in->antList);
  }
  if (in->priorRowNo)  g_free(in->priorRowNo);
  if (in->followRowNo) g_free(in->followRowNo);
  if (in->priorRow){
    for (i=0; i<in->numAnt; i++) {
      KillSYRow (in->priorRow[i]);
    }
    g_free(in->priorRow);
  }
  if (in->followRow){
    for (i=0; i<in->numAnt; i++) {
      KillSYRow (in->followRow[i]);
    }
    g_free(in->followRow);
  }

  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitSwPowerClear */

/**
 * Initialize object for a given subarray
 * Created row structures and reads prior and follow for each antenna.
 * \param in   The object to init
 * \param SubA Which (1-rel) subarray?
 * \param err  Obit error stack object.
 */
static void InitSubarray (ObitSwPower *in, olong SubA, ObitErr *err)
{
  olong i, iRow, ant, start;
  gboolean found=FALSE;
  gchar *routine = "ObitSwPower:InitSubarray";

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Clear old structures */
  if (in->priorRowNo)  g_free(in->priorRowNo);
  if (in->followRowNo) g_free(in->followRowNo);
  if (in->priorRow){
    for (i=0; i<in->numAnt; i++) {
      KillSYRow(in->priorRow[i]);
    }
    g_free(in->priorRow);
  }
  if (in->followRow){
    for (i=0; i<in->numAnt; i++) {
      KillSYRow(in->followRow[i]);
    }
    g_free(in->followRow);
  }

  /* Create new for this subarray */
  in->SubA        = SubA;
  in->numAnt      = in->antList[SubA-1]->number;
  in->priorRowNo  = g_malloc0(in->numAnt*sizeof(olong));
  in->followRowNo = g_malloc0(in->numAnt*sizeof(olong));
  in->priorRow    = g_malloc0(in->numAnt*sizeof(ObitTableSYRow*));
  in->followRow   = g_malloc0(in->numAnt*sizeof(ObitTableSYRow*));
  for (i=0; i<in->numAnt; i++) {
    in->priorRow[i]  = MakeSYRow(in->SYTable);
    in->followRow[i] = MakeSYRow(in->SYTable);
    /* Set time to large negative value to indicate uninitialized */
    in->priorRow[i]->Time    = -1.0e20;
    in->followRow[i]->Time   = -1.0e20;
    in->priorRowNo[i]        = -1;
    in->followRowNo[i]       = -1;
  }
  /* Open table */
  ObitTableSYOpen (in->SYTable, OBIT_IO_ReadOnly, err);
  if (err->error) Obit_traceback_msg (err, routine, in->SYTable->name);

  /* Create work row if needed */
  if (in->SYRow==NULL) in->SYRow = newObitTableSYRow(in->SYTable);
  
  /* Loop through antennas finding first occurance */
  for (i=0; i<in->numAnt; i++) {
    /* Loop over table */
    ant = i + 1;
    found = FALSE;
    for (iRow=1; iRow<=in->SYTable->myDesc->nrow; iRow++) {
      ObitTableSYReadRow (in->SYTable, iRow, in->SYRow, err);
      CopySYRow (in->SYRow, in->priorRow[i], in->numIF);  /* Save values */
      if (err->error) Obit_traceback_msg (err, routine, in->SYTable->name);
      in->priorRowNo[i] = iRow;   /* Save row number */
      /* This the one? */
      if ((in->priorRow[i]->antennaNo==ant) && (in->priorRow[i]->SubA==in->SubA)) 
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
      start = MIN (in->priorRowNo[i]+1, in->SYTable->myDesc->nrow);
      for (iRow=start; iRow<=in->SYTable->myDesc->nrow; iRow++) {
	ObitTableSYReadRow (in->SYTable, iRow, in->SYRow, err);
	CopySYRow (in->SYRow, in->followRow[i], in->numIF);  /* Save values */
	if (err->error) Obit_traceback_msg (err, routine, in->SYTable->name);
	in->followRowNo[i] = iRow;   /* Save row number */
	/* This the one? */
	if ((in->followRow[i]->antennaNo==ant) && (in->followRow[i]->SubA==in->SubA)) 
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
  ObitTableSYClose (in->SYTable, err);
  if (err->error) Obit_traceback_msg (err, routine, in->SYTable->name);
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
static void UpdateTime (ObitSwPower *in, ofloat time, olong Ant, 
			ObitErr *err)
{
  olong iRow, ant, start;
  gboolean found, OK, opened=FALSE;
  gchar *routine = "ObitSwPower:UpdateTime";

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Update antenna if necessary */
  ant = Ant-1;
  if ((in->followRowNo[ant]>0) &&                         /* Follow set */
      (in->followRowNo[ant]<in->SYTable->myDesc->nrow) && /* before end */
      (time>in->followRow[ant]->Time)) {
    /* Open table */
    if (!opened)
      ObitTableSYOpen (in->SYTable, OBIT_IO_ReadOnly, err);
    if (err->error) Obit_traceback_msg (err, routine, in->SYTable->name);
    opened = TRUE;
    OK     = FALSE;
    /* loop until follow is past time or the end of the table */
    while(!OK) {
      /* Update - shuffle data */
      CopySYRow(in->followRow[ant], in->priorRow[ant], in->numIF);
      in->priorRowNo[ant] = in->followRowNo[ant];
      
      /* Loop over table looking for next occurance of ant+1 */
      start = MIN (in->priorRowNo[ant]+1, in->SYTable->myDesc->nrow);
      for (iRow=start; iRow<=in->SYTable->myDesc->nrow; iRow++) {
	ObitTableSYReadRow (in->SYTable, iRow, in->SYRow, err);
	CopySYRow (in->SYRow, in->followRow[ant], in->numIF);  /* Save values */
	if (err->error) Obit_traceback_msg (err, routine, in->SYTable->name);
	in->followRowNo[ant] = iRow;   /* Save row number */
	/* This the one? match antenna and subarray */
	if ((in->followRow[ant]->antennaNo==(ant+1)) && (in->followRow[ant]->SubA==in->SubA)) 
	  {found=TRUE; break;}
      } /* end loop over table */
	/* Find next? */
      if (!found) {
	in->followRow[ant]->Time = -1.0e20;
      }
      /* This one OK or past end */
      OK = ((time<=in->followRow[ant]->Time) || (iRow>=in->SYTable->myDesc->nrow));
    } /* end !OK loop */
  } /* end needs update */
  
  /* Close table */
  if (opened)
    ObitTableSYClose (in->SYTable, err);
  if (err->error) Obit_traceback_msg (err, routine, in->SYTable->name);
} /*  end UpdateTime */

/**
 * Independent TableSY Row Constructor.
 * CANNOT be used for actual IO, only for storage
 * \param SYTab  Table for new structure
 * \return the new Table SY Row object.
 */
static ObitTableSYRow* MakeSYRow (ObitTableSY *SYTab) 
{
  ObitTableSYRow *out=NULL;
  olong numIF;

  numIF = SYTab->nIF;
  out = newObitTableSYRow(SYTab);

  /* Replace pointers to buffer with real arrays */
  out->PwrDif1 = g_malloc0(numIF*sizeof(ofloat));
  out->PwrSum1 = g_malloc0(numIF*sizeof(ofloat));
  out->Gain1   = g_malloc0(numIF*sizeof(ofloat));
  out->PwrDif2 = g_malloc0(numIF*sizeof(ofloat));
  out->PwrSum2 = g_malloc0(numIF*sizeof(ofloat));
  out->Gain2   = g_malloc0(numIF*sizeof(ofloat));
  return out;
} /* end MakeSYRow  */

/**
 * Independent TableSY Row Destructor
 * \param SYRow  Row structure to delete;
 * \return NULL
 */
static ObitTableSYRow* KillSYRow (ObitTableSYRow *SYRow) 
{
  if (SYRow==NULL) return NULL;
  if (SYRow->PwrDif1)  g_free(SYRow->PwrDif1);
  if (SYRow->PwrSum1)  g_free(SYRow->PwrSum1);
  if (SYRow->Gain1)    g_free(SYRow->Gain1);
  if (SYRow->PwrDif2)  g_free(SYRow->PwrDif2);
  if (SYRow->PwrSum2)  g_free(SYRow->PwrSum2);
  if (SYRow->Gain2)    g_free(SYRow->Gain2);
  return ObitTableSYRowUnref(SYRow);
} /* end  KillSYRow */

/**
 * Copy TableSY Row 
 * \param SYRowIn  Input Row structure
 * \param SYRowOut Output Row structure
 * \param numIF    Dimension of arrays
 */
static void CopySYRow (ObitTableSYRow* SYRowIn, ObitTableSYRow* SYRowOut, 
		       olong numIF)
{
  olong i;

  SYRowOut->Time      = SYRowIn->Time;
  SYRowOut->TimeI     = SYRowIn->TimeI;
  SYRowOut->SourID    = SYRowIn->SourID;
  SYRowOut->antennaNo = SYRowIn->antennaNo;
  SYRowOut->SubA      = SYRowIn->SubA;
  SYRowOut->FreqID    = SYRowIn->FreqID;
  for (i=0; i<numIF; i++) SYRowOut->PwrDif1[i] = SYRowIn->PwrDif1[i];
  for (i=0; i<numIF; i++) SYRowOut->PwrSum1[i] = SYRowIn->PwrSum1[i];
  for (i=0; i<numIF; i++) SYRowOut->Gain1[i]   = SYRowIn->Gain1[i];
  if (SYRowIn->PwrDif2) for (i=0; i<numIF; i++) SYRowOut->PwrDif2[i] = SYRowIn->PwrDif2[i];
  if (SYRowIn->PwrSum2) for (i=0; i<numIF; i++) SYRowOut->PwrSum2[i] = SYRowIn->PwrSum2[i];
  if (SYRowIn->Gain2)   for (i=0; i<numIF; i++) SYRowOut->Gain2[i]   = SYRowIn->Gain2[i];
} /* end  CopySYRow */

/**
 * Routine to smooth values in an open SY table.  
 * All poln present are smoothed but only one IF.
 * \param SYTab  SY table object; must be opened/closed externally
 * \param smoFunc  Smoothing function: 'MWF', 'GAUS', else BOX 
 * \param smoType  Type of data to smooth
 *                "DIF ", "ALL ", "    " = "ALL " [def "DIF "]
 * \param alpha    Alpha clip for MWF (0 -> box, 1 -> pure MWF) 
 * \param smoParm  Smoothing time in days for:
 *                 PwrDif, PwrSum, Gain
 *                 0=>fill in for blanked only. 
 * \param iif      Desired IF (0-rel)
 * \param sub      Desired subarray (1-rel)
 * \param nxt      Number of times allowed in wrk 
 * \param work1    Work buffer (nxt*16) 
 * \param work2    Work area >= (nxt*m)  (m=2 BOX, 3 GAUS, 4 MWF) 
 * \param doBlank  replace blanked values with interpolated values?
 * \param err      Error/message stack, returns if error.
 */
void 
SYSmooth (ObitTableSY *SYTab, gchar* smoFunc, gchar* smoType, ofloat alpha, 
	  ofloat *smoParm, olong iif, olong sub, 
	  olong nxt, ofloat* work1, ofloat* work2, gboolean doBlank, ObitErr* err) 
{
  olong   loopa, numtim, ant, numrec, nleft, isnrno=0, itime, 
    n1good, n2good, i, numif, numpol, numant;
  olong loopr, fstrec, save;
  ofloat    stdif=0.0, stsum=0.0, stgain, weight, fblank =  ObitMagicF();
  gboolean  need2, dodif, dosum=FALSE, dogain=FALSE;
  odouble   timoff=0.0;
  ObitIOCode retCode;
  ObitTableSYRow *row=NULL;
  gchar *routine = "ObitSwPower:SYSmooth";
  
  /* Error checks */
  if (err->error) return ;  /* previous error? */
  g_assert(ObitTableSYIsA(SYTab));
  
  /* Get number of records in table */
  numrec = SYTab->myDesc->nrow;
  if (numrec <= 0) return;   /* bail if nothing */
  
  /* Get descriptive info */
  numif  = SYTab->nIF;
  numpol = SYTab->nPol;
  numant = SYTab->nAnt;
  
 
  /* Are there 2 polarizations? */
  need2 = numpol>1;

  /* Only dif ? */
  if (!strncmp(smoType, "DIF ",4)) {
    dodif   = TRUE;
    dosum   = FALSE;
    dogain  = FALSE;
    stdif   = smoParm[0];
    stsum   = 0.0;
    stgain  = 0.0;

    /* All? */
  } else if (!strncmp(smoType, "ALL ",4)) {
    dodif   = TRUE;
    dosum   = TRUE;
    dogain  = TRUE;
    stdif   = smoParm[0];
    stsum   = smoParm[1];
    stgain  = smoParm[2];
  }

  /* Create Row */
  row = newObitTableSYRow (SYTab);
  /* Attach row to output buffer */
  ObitTableSYSetRow (SYTab, row, err);
  fstrec = 0;  /* Record number read in table */
  
  /* Loop over antenna */
  for (loopa= 1; loopa<=numant; loopa++) { /* loop 600 */
    ant = loopa;
    /* Set pointers, counters */
    numtim = 0;
    nleft = SYTab->myDesc->nrow - fstrec;  /* How many rows? */
    n1good = 0;
    n2good = 0;
    /* Loop in time, reading */
    for (loopr=1; loopr<=nleft; loopr++) { /* loop 100 */
      isnrno = fstrec + loopr;
      retCode = ObitTableSYReadRow (SYTab, isnrno, row, err);
      if (err->error) Obit_traceback_msg (err, routine, SYTab->name);
      if (row->status<0) continue;  /* Skip deselected record */
      
      /* Finished antenna? */
      if (row->antennaNo < ant) continue; /* Shouldn't happen */
      if (row->antennaNo > ant) break;

      /* Want this record */
      if ((row->SubA == sub)  &&  (row->antennaNo == ant)) {

	/* Put in buffer */
	if (numtim >= nxt) {
	  Obit_log_error(err, OBIT_Error, 
			 "%s: Exceed time limit of %d for %s", 
			 routine, nxt, SYTab->name);
	  row = ObitTableSYRowUnref(row); /* delete row object */
	  return;
	} 
	/* Work1 usage :
	   0 = dif pol 1
	   1 = sum pol 1
	   2 = gain pol 1
	   3 = Weight pol 1
	   4 = dif pol 2
	   5 = sum pol 2
	   6 = gain pol 2
	   7 = Weight pol 2
	   8 = Time(day) relative to first
	   9 = row number
	*/
	
	if (numtim == 0) timoff = row->Time;  /* First time */
	work1[8*nxt+numtim] = row->Time - timoff;
	work1[9*nxt+numtim] = (ofloat)isnrno;
	/*USE WEIGHT TO BLANK CRAZIES;*/
	weight = 1.0;
	/* First polarization */
	if (row->PwrDif1[iif]!=fblank) {
	  work1[0*nxt+numtim]  = row->PwrDif1[iif];
	  work1[1*nxt+numtim]  = row->PwrSum1[iif];
	  work1[2*nxt+numtim]  = row->Gain1[iif];
	  work1[3*nxt+numtim]  = weight;
	  n1good = n1good + 1;
	} else {
	  work1[0*nxt+numtim]  = fblank;
	  work1[1*nxt+numtim]  = fblank;
	  work1[2*nxt+numtim]  = fblank;
	  work1[3*nxt+numtim]  = fblank;
	}
	
	if (need2) {  /* Second polarization */
	  /*USE WEIGHT TO BLANK CRAZIES;*/
	  weight = 1.0;
	if (row->PwrDif2[iif]!=fblank) {
	    work1[4*nxt+numtim]  = row->PwrDif2[iif];
	    work1[5*nxt+numtim]  = row->PwrSum2[iif];
	    work1[6*nxt+numtim]  = row->Gain2[iif];
	    work1[7*nxt+numtim]  = weight;
	    n2good = n2good + 1;
	  } else {
	    work1[4*nxt+numtim]  = fblank;
	    work1[5*nxt+numtim]  = fblank;
	    work1[6*nxt+numtim]  = fblank;
	    work1[7*nxt+numtim]  = fblank;
	  } 
	} /* end second polarization */
      } /* end if want record */ 
      numtim++;   /* count times */
    } /* end loop  L100: */

    save = isnrno - 1; /* How far did we get? */
    if (numtim <= 0) goto endAnt;  /* Catch anything? */
    
    /* Smooth as requested */
    if (n1good > 0) { /* First polarization */
      if (dodif) {  /* Diff */
	SYsmoIt (smoFunc, stdif, alpha, &work1[8*nxt], &work1[0*nxt], &work1[3*nxt], numtim, 
		 &work2[0*nxt], &work2[1*nxt], &work2[2*nxt], &work2[3*nxt], doBlank);
	/* Copy back */
	for (i=0; i<numtim; i++) work1[0*nxt+i] = work2[i];
	/* Save deblanked weights if no other smoothing */
	if (doBlank && (!(dosum||dogain))) for (i=0; i<numtim; i++) work1[3*nxt+i] = work2[1*nxt+i];
      }
      if (dosum) {  /* Sum */
	SYsmoIt (smoFunc, stsum, alpha, &work1[8*nxt], &work1[1*nxt], &work1[3*nxt], numtim, 
		 &work2[0*nxt], &work2[1*nxt], &work2[2*nxt], &work2[3*nxt], doBlank);
	/* Copy back */
	for (i=0; i<numtim; i++) work1[1*nxt+i] = work2[i];
	/* Save deblanked weights if no other smoothing */
	if (doBlank && (!dogain)) for (i=0; i<numtim; i++) work1[3*nxt+i] = work2[1*nxt+i];
      }
      if (dogain) {  /* Gain */
	SYsmoIt (smoFunc, stsum, alpha, &work1[8*nxt], &work1[2*nxt], &work1[3*nxt], numtim, 
		 &work2[0*nxt], &work2[1*nxt], &work2[2*nxt], &work2[3*nxt], doBlank);
	/* Copy back */
	for (i=0; i<numtim; i++) work1[2*nxt+i] = work2[i];
	/* Save deblanked weights  */
	if (doBlank) for (i=0; i<numtim; i++) work1[3*nxt+i] = work2[1*nxt+i];
      }
    } /* end first polarization */
    
    if (n2good > 0) {  /* Second polarization */
      if (dodif) {  /* Diff */
	SYsmoIt (smoFunc, stdif, alpha, &work1[8*nxt], &work1[4*nxt], &work1[7*nxt], numtim, 
		 &work2[0*nxt], &work2[1*nxt], &work2[2*nxt], &work2[3*nxt], doBlank);
	/* Copy back */
	for (i=0; i<numtim; i++) work1[4*nxt+i] = work2[i];
	/* Save deblanked weights if no other smoothing */
	if (doBlank && (!(dosum||dogain))) for (i=0; i<numtim; i++) work1[7*nxt+i] = work2[1*nxt+i];
      }
      if (dosum) {  /* Sum */
	SYsmoIt (smoFunc, stsum, alpha, &work1[8*nxt], &work1[5*nxt], &work1[7*nxt], numtim, 
		 &work2[0*nxt], &work2[1*nxt], &work2[2*nxt], &work2[3*nxt], doBlank);
	/* Copy back */
	for (i=0; i<numtim; i++) work1[5*nxt+i] = work2[i];
	/* Save deblanked weights if no other smoothing */
	if (doBlank && (!dogain)) for (i=0; i<numtim; i++) work1[7*nxt+i] = work2[1*nxt+i];
      }
      if (dogain) {  /* Gain */
	SYsmoIt (smoFunc, stsum, alpha, &work1[8*nxt], &work1[6*nxt], &work1[7*nxt], numtim, 
		 &work2[0*nxt], &work2[1*nxt], &work2[2*nxt], &work2[3*nxt], doBlank);
	/* Copy back */
	for (i=0; i<numtim; i++) work1[6*nxt+i] = work2[i];
	/* Save deblanked weights  */
	if (doBlank) for (i=0; i<numtim; i++) work1[7*nxt+i] = work2[1*nxt+i];
      }
    } /* end second polarization */
    
    /* Replace with smoothed values */
    for (itime=0; itime<numtim; itime++) { /* loop 200 */
      isnrno = (olong)(work1[9*nxt+itime]+0.5);
      retCode = ObitTableSYReadRow (SYTab, isnrno, row, err);
      if (err->error) Obit_traceback_msg (err, routine, SYTab->name);
      if (row->status<0) continue;  /* Skip deselected record */
      
      /* Update */
      /* weights zero rather than fblank */
      if (work1[3*nxt+itime]==fblank) work1[3*nxt+itime] = 0.0;
      if (work1[3*nxt+itime]>0.0) {
	row->PwrDif1[iif]   = work1[0*nxt+itime];
	row->PwrSum1[iif]   = work1[1*nxt+itime];
	row->Gain1[iif]     = work1[2*nxt+itime];
      } else {  /* Datum bad */
	row->PwrDif1[iif]   = fblank;
	row->PwrSum1[iif]   = fblank;
	row->Gain1[iif]     = fblank;
      }
      if (need2) {
	/* weights zero rather than fblank */
	if (work1[7*nxt+itime]==fblank) work1[7*nxt+itime] = 0.0;
	if (work1[7*nxt+itime] > 0.0) {
	  row->PwrDif2[iif]   = work1[4*nxt+itime];
	  row->PwrSum2[iif]   = work1[5*nxt+itime];
	  row->Gain2[iif]     = work1[6*nxt+itime];
	} else {  /* Datum bad */
	  row->PwrDif2[iif]   = fblank;
	  row->PwrSum2[iif]   = fblank;
	  row->Gain2[iif]     = fblank;
	}
      }
      
      /* Rewrite record */
      retCode = ObitTableSYWriteRow (SYTab, isnrno, row, err);
      if (err->error) Obit_traceback_msg (err, routine, SYTab->name);
    } /* end loop rewriting smoothed solutions L200: */;
    /* First SY number of next antenna */
    
    /* End of antenna loop */
  endAnt: fstrec = save;
  } /* end loop over antennas  L600: */;

  row = ObitTableSYRowUnref(row); /* delete row object */
} /* end of routine SYSmooth */ 

/**
 * Routine to determine number of times in an open SY table.  
 * \param SYTab    SY table object 
 * \param isub     Subarray number, 0=>1 
 * \param err      Error/message stack, returns if error.
 * \return number of times
 */
static olong SYCountTime (ObitTableSY *SYTab, olong isub, ObitErr* err) 
{
  olong  loop, sub, count=0;
  ofloat lastTime;
  ObitTableSYRow *row=NULL;
  gchar *routine = "ObitSwPower:SYCountTime";

  /* Error checks */
  if (err->error) return count;  /* previous error? */
  g_assert(ObitTableSYIsA(SYTab));

  /* Subarray */
  sub = MAX (1, isub);
  lastTime = -1.0e20;
  count = 0;
  
  /* Open table */
  ObitTableSYOpen (SYTab, OBIT_IO_ReadOnly, err);
  if (err->error) Obit_traceback_val (err, routine, SYTab->name, count);

  /* Create Row */
  row = newObitTableSYRow (SYTab);
  /* Loop through table */
  for (loop=1; loop<=SYTab->myDesc->nrow; loop++) { /* loop 20 */

    ObitTableSYReadRow (SYTab, loop, row, err);
    if (err->error) break;
    if (row->status<0) continue;  /* Skip deselected record */

    /* Right subarray? */
    if ((row->SubA!=sub) && (row->SubA>0)) continue;

    /* Count times - only allow epsilon time difference */
    if (row->Time>(lastTime+0.0005*row->TimeI)) {
      lastTime = row->Time;
      count++;
    }
  } /* end loop  L20:  */

  ObitTableSYClose (SYTab, err);
  row = ObitTableSYRowUnref(row); /* delete row object */
  if (err->error) Obit_traceback_val (err, routine, SYTab->name, count);

  return count;
} /* end of routine timeCount */ 

/**
  * Routine to call appropriate smoothing routine.  Magic value blanking  
 * is supported.  
 * Routine adopted from the AIPSish 
 * 31DEC02/APL/PGM/NOTST/SNSMO.FOR/SNSMSM  
 * \param smmeth  Method 'BOX','MWF', 'GAUS', unknown = 'BOX' 
 * \param width   Smoothing time (days) 
 * \param alpha   0 -> 1 = pure boxcar -> pure MWF (ALPHA of the 
 *                data samples are discarded and the rest averaged). 
 * \param x       Abscissas of points to be smoothed in increasing 
 *                order 
 * \param y       Values to be smoothed. 
 * \param w       Weights of data. 
 * \param n       Number of points to smooth. 
 * \param ys      [out] Smoothed values. 
 * \param ws      [out] Smoothed weights 
 * \param yor     Scratch 
 * \param wor     Scratch 
 * \param doBlank replace blanked values with interpolated values.
 */
void 
SYsmoIt (gchar* smmeth, ofloat width, ofloat alpha, 
	 ofloat* x, ofloat *y, ofloat *w, olong n, 
	 ofloat* ys, ofloat* ws, ofloat *wrk1, ofloat *wrk2, gboolean doBlank) 
{
  /* Any work to do? */
  if (n <= 0) return;

  /* Smooth */
  if (!strncmp (smmeth, "BOX",3)) {
    ObitUVSolnSmooBox (width, x, y, w, n, ys, ws, doBlank);
  } else if (!strncmp (smmeth, "MWF",3)) {
    ObitUVSolnSmooMWF (width, alpha, x, y, w, n, ys, ws, wrk1, wrk2, doBlank);
  } else if (!strncmp (smmeth, "Gaus",4)) {
    ObitUVSolnSmooGauss (width, x, y, w, n, ys, ws, wrk1, doBlank);
  } else { /* Default "BOX" */
    ObitUVSolnSmooBox (width, x, y, w, n, ys, ws, doBlank);
  }
} /* end of routine SYsmoIt */ 
