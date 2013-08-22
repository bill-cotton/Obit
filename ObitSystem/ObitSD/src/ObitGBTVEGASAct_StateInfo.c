/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2013                                               */
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
/*;  Correspondence this software should be addressed as follows:     */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/

#include "ObitGBTVEGASAct_StateInfo.h"
#include "ObitTableGBTVEGASACT_STATE.h"

/*-------------- Obit: Merx mollis mortibus nuper ------------*/
/**
 * \file ObitGBTVEGASAct_StateInfo.c
 * ObitGBTVEGASAct_StateInfo class function definitions.
 *
 * This is a list of associated tables.
 */

/*------------------- File Global Variables - ----------------*/
/** name of the class defined in this file */
static gchar *myClassName = "ObitGBTVEGASAct_StateInfo";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo global structure ObitGBTVEGASAct_StateInfoClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitGBTVEGASAct_StateInfoClassInfo myClassInfo = {FALSE};

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitGBTVEGASAct_StateInfoInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitGBTVEGASAct_StateInfoClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitGBTVEGASAct_StateInfoClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Basic Constructor.
 * Initializes class if needed on first call.
 * \param name  A name for the object
 * \return the new object.
 */
ObitGBTVEGASAct_StateInfo* newObitGBTVEGASAct_StateInfo (gchar* name)
{
  ObitGBTVEGASAct_StateInfo* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitGBTVEGASAct_StateInfoClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitGBTVEGASAct_StateInfo));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

 /* set classInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitGBTVEGASAct_StateInfoInit((gpointer)out);

  return out;
} /* end newObitGBTVEGASAct_StateInfo */

/**
 * Constructor from values.
 * \param name  A name for the object
 * \param DataRoot Root of data directory.
 * \param scan     Date/time scan name (e.g. "2003_05_05_05:32:56")
 * \param err      Obit error stack object.
 * \return the new object.
 */
ObitGBTVEGASAct_StateInfo* 
newObitGBTVEGASAct_StateInfoValue (gchar *name, gchar *DataRoot, gchar *scan, ObitErr *err)
{
  ObitGBTVEGASAct_StateInfo* out=NULL;
  ObitTableGBTVEGASACT_STATE     *VEGASAct_Statetable=NULL;
  ObitTableGBTVEGASACT_STATERow  *VEGASAct_Staterow=NULL;
  ObitIOCode retCode;
  gchar *tab, FullFile[128];
  olong irow, disk=0;
  olong ver, nrow, iout;
  gchar *routine = "newObitGBTVEGASAct_StateInfoValue";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return out;
  g_assert (scan!=NULL);

  /* Create basic object */
  out = newObitGBTVEGASAct_StateInfo(name);

  /* Fill in values */
  /* Get VEGAS Act_State setup information from the VEGASAct_State table */
  /* get full file name */
  sprintf (FullFile,"%s/VEGAS/%sA.fits", DataRoot, scan);

  /* Create table structure */
  VEGASAct_Statetable = newObitTableGBTVEGASACT_STATE("VEGASAct_State");

  /* Setup */
  tab = "ACT_STATE";
  ver = 1; 
  nrow = 1;
  ObitTableSetFITS(VEGASAct_Statetable,disk,FullFile,tab,ver,nrow,err);
  
  /* Open */
  retCode = ObitTableGBTVEGASACT_STATEOpen (VEGASAct_Statetable, OBIT_IO_ReadOnly, err);
  if (err->error) return out;
  
  /* Create Row structure */
  VEGASAct_Staterow = newObitTableGBTVEGASACT_STATERow (VEGASAct_Statetable);

  /* Loop over table */
  iout = 0;
  for (irow = 1; irow<=VEGASAct_Statetable->myDesc->nrow; irow++) {
    retCode = ObitTableGBTVEGASACT_STATEReadRow (VEGASAct_Statetable, irow, VEGASAct_Staterow, err);
    if (err->error) return out;

    out->isigref1[iout] =  VEGASAct_Staterow->isigref1 != 0;
    out->isigref2[iout] =  VEGASAct_Staterow->isigref2 != 0;
    out->ical[iout]     =  VEGASAct_Staterow->ical != 0;
    out->esigref1[iout] =  VEGASAct_Staterow->esigref1 != 0;
    out->esigref2[iout] =  VEGASAct_Staterow->esigref2 != 0;
    out->ecal[iout]     =  VEGASAct_Staterow->ecal != 0;
    iout++;
    /* Done?  */
    if (iout>=MAXNUMVEGASACT_STATE) {
      Obit_log_error(err, OBIT_Error, "%s ERROR exceed limit %d VEGASAct_States", 
		     routine, MAXNUMVEGASACT_STATE);
      return out;
    }
    
  } /* end loop over table */
  out->nVEGASAct_State = iout;

  /* Close */
  retCode = ObitTableGBTVEGASACT_STATEClose (VEGASAct_Statetable, err);
  if (err->error) return out;
  
  /* Cleanup */
  VEGASAct_Statetable = ObitTableGBTVEGASACT_STATEUnref(VEGASAct_Statetable);
  VEGASAct_Staterow   = ObitTableGBTVEGASACT_STATEUnref(VEGASAct_Staterow);

  return out;
} /* end newObitGBTVEGASAct_StateInfoValue */

/**
 * Returns ClassInfo pointer for the class.
 * Initializes class if needed on first call.
 * \return pointer to the class structure.
 */
gconstpointer ObitGBTVEGASAct_StateInfoGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitGBTVEGASAct_StateInfoClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitGetIOClass */

/**
 * Make a copy of a object.
 * Parent class members are included but any derived class info is ignored.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitGBTVEGASAct_StateInfo* ObitGBTVEGASAct_StateInfoCopy  (ObitGBTVEGASAct_StateInfo *in, ObitGBTVEGASAct_StateInfo *out, 
				   ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  gboolean oldExist;
  olong i;
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
    out = newObitGBTVEGASAct_StateInfo(outName);
    if (outName) g_free(outName); outName = NULL;
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  out->nVEGASAct_State = in->nVEGASAct_State;
  for (i=0; i<in->nVEGASAct_State; i++) {
    out->isigref1[i] =  in->isigref1[i];
    out->isigref2[i] =  in->isigref2[i];
    out->ical[i]     =  in->ical[i];
    out->esigref1[i] =  in->esigref1[i];
    out->esigref2[i] =  in->esigref2[i];
    out->ecal[i]     =  in->ecal[i];
  }

  return out;
} /* end ObitGBTVEGASAct_StateInfoCopy */


/**
 * Initialize global ClassInfo Structure.
 */
void ObitGBTVEGASAct_StateInfoClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitGBTVEGASAct_StateInfoClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitGBTVEGASAct_StateInfoClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitGBTVEGASAct_StateInfoClassInfoDefFn (gpointer inClass)
{
  ObitGBTVEGASAct_StateInfoClassInfo *theClass = (ObitGBTVEGASAct_StateInfoClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitGBTVEGASAct_StateInfoClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitGBTVEGASAct_StateInfoClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitGBTVEGASAct_StateInfoGetClass;
  theClass->ObitClear     = (ObitClearFP)ObitGBTVEGASAct_StateInfoClear;
  theClass->ObitInit      = (ObitInitFP)ObitGBTVEGASAct_StateInfoInit;
  theClass->newObit       = (newObitFP)newObitGBTVEGASAct_StateInfo;
  theClass->ObitCopy      = (ObitCopyFP)ObitGBTVEGASAct_StateInfoCopy;
  theClass->ObitClone     = NULL;
 
} /* end ObitGBTVEGASAct_StateInfoClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Does (recursive) initialization of base class members before 
 * this class.
 * \param inn Pointer to the object to initialize.
 */
void ObitGBTVEGASAct_StateInfoInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  olong i;
  ObitGBTVEGASAct_StateInfo *in = inn;

  /* error checks */
  g_assert (in != NULL);
  
  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->nVEGASAct_State = 0;
  for (i=0; i<MAXNUMVEGASACT_STATE; i++) {
    in->isigref1[i] = FALSE;
    in->isigref2[i] = FALSE;
    in->ical[i]     = FALSE;
    in->esigref1[i] = FALSE;
    in->esigref2[i] = FALSE;
    in->ecal[i]     = FALSE;
  }

} /* end ObitGBTVEGASAct_StateInfoInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 */
void ObitGBTVEGASAct_StateInfoClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitGBTVEGASAct_StateInfo *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* free this class members */
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitGBTVEGASAct_StateInfoClear */

