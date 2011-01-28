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

#include "ObitTsys.h"
#include "ObitTableANUtil.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitTsys.c
 * ObitTsys class function definitions.
 * This class is derived from the Obit base class.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitTsys";

/** Function to obtain parent ClassInfo */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo structure ObitTsysClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitTsysClassInfo myClassInfo = {FALSE};

/*--------------- File Global Variables  ----------------*/


/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitTsysInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitTsysClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitTsysClassInfoDefFn (gpointer inClass);

/** Private: Initialize for a subarray. */
static void InitSubarray (ObitTsys *in, olong subarray, ObitErr *err);

/** Private: Update time of data. */
static void UpdateTime (ObitTsys *in, ofloat time, olong Ant, 
			ObitErr *err);

/** Private: Create independent TY Table row */
static ObitTableTYRow* MakeTYRow (ObitTableTY *TYTab);

/** Private: Delete independent TY Table row */
static ObitTableTYRow* KillTYRow (ObitTableTYRow *TYRow);

/** Private: Copy TY Table row */
static void CopyTYRow (ObitTableTYRow* TYTabIn, ObitTableTYRow* SyTabOut, 
		       olong numIF);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitTsys* newObitTsys (gchar* name)
{
  ObitTsys* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitTsysClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitTsys));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitTsysInit((gpointer)out);

 return out;
} /* end newObitTsys */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitTsysGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitTsysClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitTsysGetClass */

/**
 * Make a deep copy of an ObitTsys.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitTsys* ObitTsysCopy  (ObitTsys *in, ObitTsys *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  gboolean oldExist;
  olong i;
  gchar *outName;
  gchar *routine = "ObitTsysCopy";

  /* error checks */
  if (err->error) return out;
  g_assert (ObitIsA(in, &myClassInfo));
  if (out) g_assert (ObitIsA(out, &myClassInfo));

  /* Create if it doesn't exist */
  oldExist = out!=NULL;
  if (!oldExist) {
    /* derive object name */
    outName = g_strconcat ("Copy: ",in->name,NULL);
    out = newObitTsys(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  /* Clear old structures */
  out->TYTable = ObitTableTYUnref(out->TYTable);
  if (out->antList) {
    for (i=0; i<out->numSubA; i++)
      out->antList[i]   = ObitAntennaListUnref(out->antList[i]);
    g_free(out->antList);
  }
  if (out->priorRowNo)  g_free(out->priorRowNo);
  if (out->followRowNo) g_free(out->followRowNo);
  if (out->priorRow){
    for (i=0; i<out->numAnt; i++) {
      KillTYRow(out->priorRow[i]);
    }
    g_free(out->priorRow);
  }
  if (out->followRow){
    for (i=0; i<out->numAnt; i++) {
      KillTYRow(out->followRow[i]);
    }
    g_free(out->followRow);
  }

  /* Create new  */
  out->TYTable  = ObitTableTYRef(in->TYTable);
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
  out->priorRow    = g_malloc0(out->numAnt*sizeof(ObitTableTYRow*));
  out->followRow   = g_malloc0(out->numAnt*sizeof(ObitTableTYRow*));
  for (i=0; i<out->numAnt; i++) {
    out->priorRow[i]  = MakeTYRow(out->TYTable);
    out->followRow[i] = MakeTYRow(out->TYTable);
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
} /* end ObitTsysCopy */

/**
 * Make a copy of a object but do not copy the actual data
 * This is useful to create an Tsys similar to the input one.
 * \param in  The object to copy
 * \param out An existing object pointer for output, must be defined.
 * \param err Obit error stack object.
 */
void ObitTsysClone  (ObitTsys *in, ObitTsys *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  olong i;
  gchar *routine = "ObitTsysClone";

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
  out->TYTable = ObitTableTYUnref(out->TYTable);
  if (out->antList) {
    for (i=0; i<out->numSubA; i++)
      out->antList[i]   = ObitAntennaListUnref(out->antList[i]);
    g_free(out->antList);
  }
  if (out->priorRowNo)  g_free(out->priorRowNo);
  if (out->followRowNo) g_free(out->followRowNo);
  if (out->priorRow){
    for (i=0; i<out->numAnt; i++) {
      KillTYRow(out->priorRow[i]);
    }
    g_free(out->priorRow);
  }
  if (out->followRow){
    for (i=0; i<out->numAnt; i++) {
      KillTYRow(out->followRow[i]);
    }
    g_free(out->followRow);
  }

  /* Create new  */
  out->TYTable  = ObitTableTYRef(in->TYTable);
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
  out->priorRow    = g_malloc0(out->numAnt*sizeof(ObitTableTYRow*));
  out->followRow   = g_malloc0(out->numAnt*sizeof(ObitTableTYRow*));
  for (i=0; i<out->numAnt; i++) {
    out->priorRow[i]  = MakeTYRow(out->TYTable);
    out->followRow[i] = MakeTYRow(out->TYTable);
    /* Set time to large negative value to indicate uninitialized */
    out->priorRow[i]->Time    = -1.0e20;
    out->followRow[i]->Time   = -1.0e20;
    out->priorRowNo[i]        = -1;
    out->followRowNo[i]       = -1;
  }

  /* Init for Subarray out->SubA */
  InitSubarray(out, out->SubA, err);
  if (err->error) Obit_traceback_msg (err, routine, out->name);

} /* end ObitTsysClone */

/**
 * Creates an ObitTsys 
 * \param name    An optional name for the object.
 * \param TYTable Tsys table to interpolate
 * \param UVData  UV data used to get antenna/array information
 * \param err     Obit error stack object.
 * \return the new object.
 */
ObitTsys* ObitTsysCreate (gchar* name, ObitTableTY *TYTable, 
				ObitUV *UVData, ObitErr *err)
{
  ObitTsys* out;
  ObitTableAN  *ANTable=NULL;
  olong highVer, ver, iver;
  olong numOrb=0, numPCal=0;
  gchar *routine = "ObitTsysCreate";

  /* Create basic structure */
  out = newObitTsys (name);

  /* Save table */
  out->TYTable = ObitTableTYRef(TYTable);

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
} /* end ObitTsysCreate */

/**
 *  Interpolate Tsys values at a given time
 *  Calls are expected in time order by antenna.
 * \param in      Tsys interpolator
 * \param time    Desired time (day)
 * \param ant     Desired antenna number
 * \param suba    Desired subarray
 * \param TSys1   [out] Tsys pol 1 (K), may be fblanked, 
 *                      must be dimensioned at least numIF
 * \param TSys2   [out] Tsys pol 2 (K), may be fblanked
 * \param err     Obit error stack object.
 */
void ObitTsysReport (ObitTsys *in, ofloat time, olong ant, olong suba,
		     ofloat *Tsys1, ofloat *Tsys2,
		     ObitErr *err)
{
  ofloat w1, w2, ww1, ww2, diff;
  ofloat fblank = ObitMagicF();
  olong i, iant = ant-1;
  gchar *routine = "ObitTsysReport";
  
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
  /* First poln */
  for (i=0; i<in->numIF; i++) {
    ww1 = w1; ww2 = w2;
    if (in->priorRow[iant]->Tsys1[i]==fblank)  ww1 = 0.0;
    if (in->followRow[iant]->Tsys1[i]==fblank) ww2 = 0.0;
    if ((ww1+ww2)>0.0) Tsys1[i] = ww1*in->priorRow[iant]->Tsys1[i] + 
			 ww2*in->followRow[iant]->Tsys1[i];
    else Tsys1[i] = fblank;
  } /* end IF loop poln 1 */

  /* Second  poln if present */
  if (in->numPoln>1) {
    for (i=0; i<in->numIF; i++) {
      ww1 = w1; ww2 = w2;
      if (in->priorRow[iant]->Tsys2[i]==fblank)  ww1 = 0.0;
      if (in->followRow[iant]->Tsys2[i]==fblank) ww2 = 0.0;
      if ((ww1+ww2)>0.0) Tsys2[i] = ww1*in->priorRow[iant]->Tsys2[i] + 
			   ww2*in->followRow[iant]->Tsys2[i];
      else Tsys2[i] = fblank;
    } /* end IF loop poln 1 */
  }  /* end seconf poln */

} /* end ObitTsysReport */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitTsysClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitTsysClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitTsysClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitTsysClassInfoDefFn (gpointer inClass)
{
  ObitTsysClassInfo *theClass = (ObitTsysClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitTsysClassInit;
  theClass->newObit       = (newObitFP)newObitTsys;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitTsysClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitTsysGetClass;
  theClass->ObitCopy      = (ObitCopyFP)ObitTsysCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitTsysClear;
  theClass->ObitInit      = (ObitInitFP)ObitTsysInit;
  theClass->ObitTsysCreate = (ObitTsysCreateFP)ObitTsysCreate;

} /* end ObitTsysClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitTsysInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitTsys *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->TYTable     = NULL;
  in->TYRow       = NULL;
  in->numSubA     =  0;
  in->antList     = NULL;
  in->SubA        = 0;
  in->numAnt      = 0;
  in->priorRowNo  = NULL;
  in->followRowNo = NULL;
  in->priorRow    = NULL;
  in->followRow   = NULL;

} /* end ObitTsysInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *             Actually it should be an ObitTsys* cast to an Obit*.
 */
void ObitTsysClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitTsys *in = inn;
  olong i;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  in->TYTable = ObitTableTYUnref(in->TYTable);
  if (in->antList) {
    for (i=0; i<in->numSubA; i++)
      in->antList[i]   = ObitAntennaListUnref(in->antList[i]);
    g_free(in->antList);
  }
  if (in->priorRowNo)  g_free(in->priorRowNo);
  if (in->followRowNo) g_free(in->followRowNo);
  if (in->priorRow){
    for (i=0; i<in->numAnt; i++) {
      KillTYRow(in->priorRow[i]);
    }
    g_free(in->priorRow);
  }
  if (in->followRow){
    for (i=0; i<in->numAnt; i++) {
      KillTYRow(in->followRow[i]);
    }
    g_free(in->followRow);
  }

  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitTsysClear */

/**
 * Initialize object for a given subarray
 * Created row structures and reads prior and follow for each antenna.
 * \param in   The object to init
 * \param SubA Which (1-rel) subarray?
 * \param err  Obit error stack object.
 */
static void InitSubarray (ObitTsys *in, olong SubA, ObitErr *err)
{
  olong i, iRow, ant, start;
  gboolean found;
  gchar *routine = "ObitTsys:InitSubarray";

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Clear old structures */
  if (in->priorRowNo)  g_free(in->priorRowNo);
  if (in->followRowNo) g_free(in->followRowNo);
  if (in->priorRow){
    for (i=0; i<in->numAnt; i++) {
      KillTYRow(in->priorRow[i]);
    }
    g_free(in->priorRow);
  }
  if (in->followRow){
    for (i=0; i<in->numAnt; i++) {
      KillTYRow(in->followRow[i]);
    }
    g_free(in->followRow);
  }

  /* Create new for this subarray */
  in->SubA        = SubA;
  in->numAnt      = in->antList[SubA-1]->number;
  in->priorRowNo  = g_malloc0(in->numAnt*sizeof(olong));
  in->followRowNo = g_malloc0(in->numAnt*sizeof(olong));
  in->priorRow    = g_malloc0(in->numAnt*sizeof(ObitTableTYRow*));
  in->followRow   = g_malloc0(in->numAnt*sizeof(ObitTableTYRow*));
  for (i=0; i<in->numAnt; i++) {
    in->priorRow[i]  = MakeTYRow(in->TYTable);
    in->followRow[i] = MakeTYRow(in->TYTable);
    /* Set time to large negative value to indicate uninitialized */
    in->priorRow[i]->Time    = -1.0e20;
    in->followRow[i]->Time   = -1.0e20;
    in->priorRowNo[i]        = -1;
    in->followRowNo[i]       = -1;
  }
  /* Open table */
  ObitTableTYOpen (in->TYTable, OBIT_IO_ReadOnly, err);
  if (err->error) Obit_traceback_msg (err, routine, in->TYTable->name);
  
  /* Loop through antennas finding first occurance */
  for (i=0; i<in->numAnt; i++) {
    /* Loop over table */
    ant = i + 1;
    found = FALSE;
    for (iRow=1; iRow<=in->TYTable->myDesc->nrow; iRow++) {
      ObitTableTYReadRow (in->TYTable, iRow, in->TYRow, err);
      CopyTYRow (in->TYRow, in->priorRow[i], in->numIF);  /* Save values */
      if (err->error) Obit_traceback_msg (err, routine, in->TYTable->name);
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
      start = MIN (in->priorRowNo[i]+1, in->TYTable->myDesc->nrow);
      for (iRow=start; iRow<=in->TYTable->myDesc->nrow; iRow++) {
	ObitTableTYReadRow (in->TYTable, iRow, in->TYRow, err);
	CopyTYRow (in->TYRow, in->followRow[i], in->numIF);  /* Save values */
	if (err->error) Obit_traceback_msg (err, routine, in->TYTable->name);
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
  ObitTableTYClose (in->TYTable, err);
  if (err->error) Obit_traceback_msg (err, routine, in->TYTable->name);
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
static void UpdateTime (ObitTsys *in, ofloat time, olong Ant, 
			ObitErr *err)
{
  olong i, iRow, ant, start;
  gboolean found, OK, opened=FALSE;
  gchar *routine = "ObitTsys:UpdateTime";

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Update antenna if necessary */
  ant = Ant-1;
  if ((in->followRowNo[ant]>0) &&                         /* Follow set */
      (in->followRowNo[ant]<in->TYTable->myDesc->nrow) && /* before end */
      (time>in->followRow[ant]->Time)) {
    /* Open table */
    if (!opened)
      ObitTableTYOpen (in->TYTable, OBIT_IO_ReadOnly, err);
    if (err->error) Obit_traceback_msg (err, routine, in->TYTable->name);
    opened = TRUE;
    OK     = FALSE;
    /* loop until follow is past time or the end of the table */
    while(!OK) {
      /* Update - shuffle data */
      CopyTYRow(in->followRow[ant], in->priorRow[ant], in->numIF);
      in->priorRowNo[ant] = in->followRowNo[ant];
      
      /* Loop over table looking for next occurance of ant+1 */
      start = MIN (in->priorRowNo[i]+1, in->TYTable->myDesc->nrow);
      for (iRow=start; iRow<=in->TYTable->myDesc->nrow; iRow++) {
	ObitTableTYReadRow (in->TYTable, iRow, in->TYRow, err);
	CopyTYRow (in->TYRow, in->followRow[i], in->numIF);  /* Save values */
	if (err->error) Obit_traceback_msg (err, routine, in->TYTable->name);
	in->followRowNo[i] = iRow;   /* Save row number */
	/* This the one? match antenna and subarray */
	if ((in->followRow[i]->antennaNo==(ant+1)) && (in->followRow[i]->SubA==in->SubA)) 
	  {found=TRUE;break;}
      } /* end loop over table */
	/* Find next? */
      if (!found) {
	in->followRow[ant]->Time = -1.0e20;
      }
      /* This one OK or past end */
      OK = ((time<=in->followRow[ant]->Time) || (iRow>=in->TYTable->myDesc->nrow));
    } /* end !OK loop */
  } /* end needs update */
  
  /* Close table */
  if (opened)
    ObitTableTYClose (in->TYTable, err);
  if (err->error) Obit_traceback_msg (err, routine, in->TYTable->name);
} /*  end UpdateTime */

/**
 * Independent TableTY Row Constructor.
 * CANNOT be used for actual IO, only for storage
 * \param TYTab  Table for new structure
 * \return the new Table TY Row object.
 */
static ObitTableTYRow* MakeTYRow (ObitTableTY *TYTab) 
{
  ObitTableTYRow *out=NULL;
  olong numIF;

  numIF = TYTab->numIF;
  out = newObitTableTYRow(TYTab);

  /* Replace pointers to buffer with real arrays */
  out->Tsys1 = g_malloc0(numIF*sizeof(ofloat));
  out->Tant1 = g_malloc0(numIF*sizeof(ofloat));
  out->Tsys2 = g_malloc0(numIF*sizeof(ofloat));
  out->Tant2 = g_malloc0(numIF*sizeof(ofloat));
  return out;
} /* end MakeTYRow  */

/**
 * Independent TableTY Row Destructor
 * \param TYRow  Row structure to delete;
 * \return NULL
 */
static ObitTableTYRow* KillTYRow (ObitTableTYRow *TYRow) 
{
  if (TYRow==NULL) return NULL;
  if (TYRow->Tsys1)  g_free(TYRow->Tsys1);
  if (TYRow->Tant1)  g_free(TYRow->Tant1);
  if (TYRow->Tsys2)  g_free(TYRow->Tsys2);
  if (TYRow->Tant2)  g_free(TYRow->Tant2);
  return ObitTableTYRowUnref(TYRow);
} /* end  KillTYRow */

/**
 * Copy TableTY Row 
 * \param TYRowIn  Input Row structure
 * \param TYRowOut Output Row structure
 * \param numIF    Dimension of arrays
 */
static void CopyTYRow (ObitTableTYRow* TYRowIn, ObitTableTYRow* TYRowOut, 
		       olong numIF)
{
  olong i;

  TYRowOut->Time      = TYRowIn->Time;
  TYRowOut->TimeI     = TYRowIn->TimeI;
  TYRowOut->SourID    = TYRowIn->SourID;
  TYRowOut->antennaNo = TYRowIn->antennaNo;
  TYRowOut->SubA      = TYRowIn->SubA;
  TYRowOut->FreqID    = TYRowIn->FreqID;
  for (i=0; i<numIF; i++) TYRowOut->Tsys1[i] = TYRowIn->Tsys1[i];
  for (i=0; i<numIF; i++) TYRowOut->Tant1[i] = TYRowIn->Tant1[i];
  if (TYRowIn->Tsys2) for (i=0; i<numIF; i++) TYRowOut->Tsys2[i] = TYRowIn->Tsys2[i];
  if (TYRowIn->Tant2) for (i=0; i<numIF; i++) TYRowOut->Tant2[i] = TYRowIn->Tant2[i];
} /* end  CopyTYRow */

