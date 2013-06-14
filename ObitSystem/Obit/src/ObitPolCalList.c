/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2012-2013                                          */
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

#include "ObitPolCalList.h"
#include "ObitTablePD.h"

/*-------------- Obit: Merx mollis mortibus nuper ------------*/
/**
 * \file ObitPolCalList.c
 * ObitPolCalList class function definitions.
 *
 * This is a list of polarization calibration parameters derived 
 * from an AIPS PD table.
 */

/*------------------- File Global Variables - ----------------*/
/** name of the class defined in this file */
static gchar *myClassName = "ObitPolCalList";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo global structure ObitPolCalListClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitPolCalListClassInfo myClassInfo = {FALSE};

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitPolCalListInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitPolCalListClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitPolCalListClassInfoDefFn (gpointer inClass);


/*----------------------Public functions---------------------------*/
/**
 * Basic Constructor.
 * Initializes class if needed on first call.
 * \return the new object.
 */
ObitPolCalList* newObitPolCalList (gchar* name)
{
  ObitPolCalList* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitPolCalListClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitPolCalList));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

 /* set classInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitPolCalListInit((gpointer)out);

  return out;
} /* end newObitPolCalList */

/**
 * Returns ClassInfo pointer for the class.
 * Initializes class if needed on first call.
 * \return pointer to the class structure.
 */
gconstpointer ObitPolCalListGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitPolCalListClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitGetIOClass */

/**
 * Creates an ObitPolCalList of a given size.
 * \param name       An optional name for the object.
 * \param PDTabl     PD Table as basic Obit
 * \param err        ObitError stack.
 * \return the new object.
 */
ObitPolCalList* ObitPolCalListCreate (gchar* name, Obit *PDTabl, 
				      ObitErr *err)
{
  ObitPolCalList* out=NULL;
  ObitIOCode retCode;
  ObitTablePD *PDTab = (ObitTablePD*)PDTabl;
  ObitTablePDRow *PDRow=NULL;
  olong i, iRow, ia, iif, ich, indx, jndx, kndx;
  gboolean haveRLPhase;
  gchar *routine = "ObitPolCalListCreate";

  if (err->error) return out;
  /* Create basic structure */
  out = newObitPolCalList (name);

  /* Does PD table exist? */
  Obit_retval_if_fail((PDTab!=NULL), err, out,
		      "%s: Could not find PD table for %s", 
		      routine, PDTab->name);
  
  /* Open input table */
  retCode = ObitTablePDOpen (PDTab, OBIT_IO_ReadOnly, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, PDTab->name, out);
  
  /* Set row */
  PDRow  = newObitTablePDRow (PDTab);

  /* Sizes of things */
  out->numAnt    = PDTab->numAnt;
  out->numPoln   = PDTab->numPol;
  out->numIF     = PDTab->numIF;
  out->numChan   = PDTab->numChan;
  out->numPCal   = PDTab->numPol*2;
  out->polRefAnt = -1;

  /* create arrays */
  out->ANlist      = g_malloc0 (PDTab->numAnt*sizeof(ofloat*));
  out->RLPhaseDiff = g_malloc0 (PDTab->numIF*PDTab->numChan*sizeof(ofloat));

  /* Create elements of ObitPolCal list */
  for (i=0; i<PDTab->numAnt; i++) {
    out->ANlist[i] = g_malloc0(out->numPCal*out->numIF*out->numChan*sizeof(ofloat));
  }

  /* Polarization parameterization type  */
  out->polType = ObitPolCalListGetPolType (PDTab->polType);

  haveRLPhase = PDTab->myDesc->repeat[PDTab->RLPhaseCol]>0;
  /* Loop over table */
  for (iRow=1; iRow<=PDTab->myDesc->nrow; iRow++) {
    retCode = ObitTablePDReadRow (PDTab, iRow, PDRow, err);
    if ((retCode != OBIT_IO_OK) || (err->error))
      Obit_traceback_val (err, routine, PDTab->name, out);
    if (PDRow->status==-1) continue;

    if (out->polRefAnt<0) out->polRefAnt = PDRow->RefAnt;  /* Get first refAnt */

    ia = PDRow->antNo - 1;
    kndx = jndx = indx = 0;
    for (iif=0; iif<out->numIF; iif++) {
      for (ich=0; ich<out->numChan; ich++) {
	/* See if RL (XY) Phase difference given */
	if (haveRLPhase) 
	  out->RLPhaseDiff[jndx++] = PDRow->RLPhase[kndx];
	else 
	  out->RLPhaseDiff[jndx++] = 0.0;
	out->ANlist[ia][indx++] = PDRow->Real1[kndx];
	out->ANlist[ia][indx++] = PDRow->Imag1[kndx];
	out->ANlist[ia][indx++] = PDRow->Real2[kndx];
	out->ANlist[ia][indx++] = PDRow->Imag2[kndx++];
      } /* end Channel loop */
    } /* end IF loop */
  } /* End loop over table */

  /* Close table */
  retCode = ObitTablePDClose (PDTab, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, PDTab->name, out);
  
  /* release row object */
  PDRow  = ObitTablePDRowUnref(PDRow);

  if (out->polRefAnt<0) out->polRefAnt = 1;  /* Default reference antenna */
  return out;
} /* end ObitPolCalListCreate */

/**
 * Return corresponding type code, recognizes:
 * \li 'ORI-ELP ' Orientation-ellipticity (OBIT_UVPoln_ELORI)
 * \li 'APPROX  ' R/L Linear D term approximation (OBIT_UVPoln_Approx)
 * \li 'VLBI    ' R/L Linear D term approximation for resolved sources (OBIT_UVPoln_VLBI)
 * \li 'X-Y LIN ' X/Y Linear D term approximation (OBIT_UVPoln_XYLin)
 * \param type Polarization calibration type name
 * \return code, OBIT_UVPoln_Unknown -> not recognized, 
 *               OBIT_UVPoln_NoCal   -> no code (no poln cal)
 */
ObitUVPolCalType ObitPolCalListGetPolType (gchar* type)
{
  ObitUVPolCalType out = OBIT_UVPoln_Unknown;

  if (!strncmp (type, "    ",    4)) out = OBIT_UVPoln_NoCal;
  if (!strncmp (type, "ORI-ELP", 7)) out = OBIT_UVPoln_ELORI;
  if (!strncmp (type, "APPROX",  6)) out = OBIT_UVPoln_Approx;
  if (!strncmp (type, "VLBI",    4)) out = OBIT_UVPoln_VLBI;
  if (!strncmp (type, "X-Y LIN", 7)) out = OBIT_UVPoln_XYLin;
  return out;
} /* end ObitPolCalListGetPolType */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitPolCalListClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitPolCalListClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitPolCalListClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitPolCalListClassInfoDefFn (gpointer inClass)
{
  ObitPolCalListClassInfo *theClass = (ObitPolCalListClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitPolCalListClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitPolCalListClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitPolCalListGetClass;
  theClass->newObit       = (newObitFP)newObitPolCalList;
  theClass->ObitCopy      = NULL;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitPolCalListClear;
  theClass->ObitInit      = (ObitInitFP)ObitPolCalListInit;

} /* end ObitPolCalListClassDefFn */


/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Does (recursive) initialization of base class members before 
 * this class.
 * \param inn Pointer to the object to initialize.
 */
void ObitPolCalListInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitPolCalList *in = inn;

  /* error checks */
  g_assert (in != NULL);
  
  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->ANlist      = NULL;
  in->RLPhaseDiff = NULL;
  in->numAnt      = 0;
} /* end ObitPolCalListInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 */
void ObitPolCalListClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitPolCalList *in = inn;
  olong i;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* free this class members */
  /* loop through list copying elements */
  for (i=0; i<in->numAnt; i++) 
    if (in->ANlist[i]) g_free(in->ANlist[i]);

  /* delete members */
  if (in->ANlist) g_free(in->ANlist); in->ANlist = NULL;
  if (in->RLPhaseDiff) g_free(in->RLPhaseDiff); in->RLPhaseDiff = NULL;
 
 /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);

} /* end ObitPolCalListClear */


  
