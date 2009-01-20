/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2004-2008                                          */
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

#include "ObitGBTDCRStateInfo.h"
#include "ObitTableGBTDCRSTATE.h"
#include "ObitTableGBTSPDATA.h"

/*-------------- Obit: Merx mollis mortibus nuper ------------*/
/**
 * \file ObitGBTDCRStateInfo.c
 * ObitGBTDCRStateInfo class function definitions.
 *
 * This is a list of associated tables.
 */

/*------------------- File Global Variables - ----------------*/
/** name of the class defined in this file */
static gchar *myClassName = "ObitGBTDCRStateInfo";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo global structure ObitGBTDCRStateInfoClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitGBTDCRStateInfoClassInfo myClassInfo = {FALSE};

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitGBTDCRStateInfoInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitGBTDCRStateInfoClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitGBTDCRStateInfoClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Basic Constructor.
 * Initializes class if needed on first call.
 * \param name  A name for the object
 * \return the new object.
 */
ObitGBTDCRStateInfo* newObitGBTDCRStateInfo (gchar* name)
{
  ObitGBTDCRStateInfo* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitGBTDCRStateInfoClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitGBTDCRStateInfo));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

 /* set classInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitGBTDCRStateInfoInit((gpointer)out);

  return out;
} /* end newObitGBTDCRStateInfo */

/**
 * Constructor from values.
 * \param name  A name for the object
 * \param disk     Obit FITS disk number
 * \param scan     Date/time scan name (e.g. "2003_05_05_05:32:56")
 * \param err      Obit error stack object.
 * \return the new object.
 */
ObitGBTDCRStateInfo* 
newObitGBTDCRStateInfoValue (gchar *name, olong disk, gchar *scan, ObitErr *err)
{
  ObitGBTDCRStateInfo* out=NULL;
  ObitTableGBTDCRSTATE     *DCRStatetable=NULL;
  ObitTableGBTDCRSTATERow  *DCRStaterow=NULL;
  ObitIOCode retCode;
  gchar *tab, FullFile[128];
  olong irow;
  olong ver, nrow, iout;
  gchar *routine = "newObitGBTDCRStateInfoValue";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return out;
  g_assert (scan!=NULL);

  /* Create basic object */
  out = newObitGBTDCRStateInfo(name);

  /* Fill in values */
  /* Get DCR State setup information from the DCRState table */
  /* get full file name */
  sprintf (FullFile,"DCR/%s.fits", scan);

  /* Create table structure */
  DCRStatetable = newObitTableGBTDCRSTATE("DCRState");

  /* Setup */
  tab = "STATE";
  ver = 1; 
  nrow = 1;
  ObitTableSetFITS(DCRStatetable,disk,FullFile,tab,ver,nrow,err);
  
  /* Open */
  retCode = ObitTableGBTDCRSTATEOpen (DCRStatetable, OBIT_IO_ReadOnly, err);
  if (err->error) return out;
  
  /* Create Row structure */
  DCRStaterow = newObitTableGBTDCRSTATERow (DCRStatetable);

  /* Loop over table */
  iout = 0;
  for (irow = 1; irow<=DCRStatetable->myDesc->nrow; irow++) {
    retCode = ObitTableGBTDCRSTATEReadRow (DCRStatetable, irow, DCRStaterow, err);
    if (err->error) return out;

    out->blanktim[iout] =  DCRStaterow->blanktim;
    out->phasetim[iout] =  DCRStaterow->phasetim;
    out->sigref[iout]   =  DCRStaterow->sigref;
    out->cal[iout]      =  DCRStaterow->cal;
    iout++;
    /* Done?  */
    if (iout>=MAXNUMDCRSTATE) {
      Obit_log_error(err, OBIT_Error, "%s ERROR exceed limit %d DCRStates", 
		     routine, MAXNUMDCRSTATE);
      return out;
    }
    
  } /* end loop over table */
  out->nDCRState = iout;

  /* Close */
  retCode = ObitTableGBTDCRSTATEClose (DCRStatetable, err);
  if (err->error) return out;
  
  /* Cleanup */
  DCRStatetable = ObitTableGBTDCRSTATEUnref(DCRStatetable);
  DCRStaterow   = ObitTableGBTDCRSTATEUnref(DCRStaterow);

  return out;
} /* end newObitGBTDCRStateInfoValue */

/**
 * Returns ClassInfo pointer for the class.
 * Initializes class if needed on first call.
 * \return pointer to the class structure.
 */
gconstpointer ObitGBTDCRStateInfoGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitGBTDCRStateInfoClassInit();

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
ObitGBTDCRStateInfo* ObitGBTDCRStateInfoCopy  (ObitGBTDCRStateInfo *in, ObitGBTDCRStateInfo *out, 
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
    out = newObitGBTDCRStateInfo(outName);
    if (outName) g_free(outName); outName = NULL;
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  out->nDCRState = in->nDCRState;
  for (i=0; i<in->nDCRState; i++) {
    out->blanktim[i]  = in->blanktim[i];
    out->phasetim[i]  = in->phasetim[i];
    out->sigref[i]    = in->sigref[i];
    out->cal[i]       = in->cal[i];
  }

  return out;
} /* end ObitGBTDCRStateInfoCopy */


/**
 * Initialize global ClassInfo Structure.
 */
void ObitGBTDCRStateInfoClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitGBTDCRStateInfoClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitGBTDCRStateInfoClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitGBTDCRStateInfoClassInfoDefFn (gpointer inClass)
{
  ObitGBTDCRStateInfoClassInfo *theClass = (ObitGBTDCRStateInfoClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitGBTDCRStateInfoClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitGBTDCRStateInfoClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitGBTDCRStateInfoGetClass;
  theClass->ObitClear     = (ObitClearFP)ObitGBTDCRStateInfoClear;
  theClass->ObitInit      = (ObitInitFP)ObitGBTDCRStateInfoInit;
  theClass->newObit       = (newObitFP)newObitGBTDCRStateInfo;
  theClass->ObitCopy      = (ObitCopyFP)ObitGBTDCRStateInfoCopy;
  theClass->ObitClone     = NULL;
 
} /* end ObitGBTDCRStateInfoClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Does (recursive) initialization of base class members before 
 * this class.
 * \param inn Pointer to the object to initialize.
 */
void ObitGBTDCRStateInfoInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  olong i;
  ObitGBTDCRStateInfo *in = inn;

  /* error checks */
  g_assert (in != NULL);
  
  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->nDCRState = 0;
  for (i=0; i<MAXNUMDCRSTATE; i++) {
    in->blanktim[i]  = 0.0;
    in->phasetim[i]  = 0.0;
    in->sigref[i]    = 0;
    in->cal[i]       = 0;
  }

} /* end ObitGBTDCRStateInfoInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 */
void ObitGBTDCRStateInfoClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitGBTDCRStateInfo *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* free this class members */
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitGBTDCRStateInfoClear */

