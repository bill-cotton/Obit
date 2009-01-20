/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2006                                               */
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

#include "ObitGBTCCBStateInfo.h"
#include "ObitTableGBTCCBSTATE.h"

/*-------------- Obit: Merx mollis mortibus nuper ------------*/
/**
 * \file ObitGBTCCBStateInfo.c
 * ObitGBTCCBStateInfo class function definitions.
 *
 * This is a list of associated tables.
 */

/*------------------- File Global Variables - ----------------*/
/** name of the class defined in this file */
static gchar *myClassName = "ObitGBTCCBStateInfo";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo global structure ObitGBTCCBStateInfoClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitGBTCCBStateInfoClassInfo myClassInfo = {FALSE};

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitGBTCCBStateInfoInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitGBTCCBStateInfoClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitGBTCCBStateInfoClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Basic Constructor.
 * Initializes class if needed on first call.
 * \param name  A name for the object
 * \return the new object.
 */
ObitGBTCCBStateInfo* newObitGBTCCBStateInfo (gchar* name)
{
  ObitGBTCCBStateInfo* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitGBTCCBStateInfoClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitGBTCCBStateInfo));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

 /* set classInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitGBTCCBStateInfoInit((gpointer)out);

  return out;
} /* end newObitGBTCCBStateInfo */

/**
 * Constructor from values.
 * \param name  A name for the object
 * \param disk     Obit FITS disk number
 * \param scan     Date/time scan name (e.g. "2003_05_05_05:32:56")
 * \param err      Obit error stack object.
 * \return the new object.
 */
ObitGBTCCBStateInfo* 
newObitGBTCCBStateInfoValue (gchar *name, olong disk, gchar *scan, ObitErr *err)
{
  ObitGBTCCBStateInfo* out=NULL;
  ObitTableGBTCCBSTATE     *CCBStatetable=NULL;
  ObitTableGBTCCBSTATERow  *CCBStaterow=NULL;
  ObitIOCode retCode;
  gchar *tab, FullFile[128];
  olong irow;
  olong ver, nrow, iout;
  gchar *routine = "newObitGBTCCBStateInfoValue";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return out;
  g_assert (scan!=NULL);

  /* Create basic object */
  out = newObitGBTCCBStateInfo(name);

  /* Fill in values */
  /* Get CCB State setup information from the CCBState table */
  /* get full file name */
  sprintf (FullFile,"CCB26_40/%s.fits", scan);

  /* Create table structure */
  CCBStatetable = newObitTableGBTCCBSTATE("CCBState");

  /* Setup */
  tab = "CCBSTATE";
  ver = 1; 
  nrow = 1;
  ObitTableSetFITS(CCBStatetable,disk,FullFile,tab,ver,nrow,err);
  
  /* Open */
  retCode = ObitTableGBTCCBSTATEOpen (CCBStatetable, OBIT_IO_ReadOnly, err);
  if (err->error) return out;
  
  /* Create Row structure */
  CCBStaterow = newObitTableGBTCCBSTATERow (CCBStatetable);

  /* Loop over table */
  iout = 0;
  for (irow = 1; irow<=CCBStatetable->myDesc->nrow; irow++) {
    retCode = ObitTableGBTCCBSTATEReadRow (CCBStatetable, irow, CCBStaterow, err);
    if (err->error) return out;

    out->phia[iout] =  CCBStaterow->phia;
    out->phib[iout] =  CCBStaterow->phib;
    iout++;
    /* Done?  */
    if (iout>=MAXNUMCCBSTATE) {
      Obit_log_error(err, OBIT_Error, "%s ERROR exceed limit %d CCBStates", 
		     routine, MAXNUMCCBSTATE);
      return out;
    }
    
  } /* end loop over table */
  out->nCCBState = iout;

  /* Close */
  retCode = ObitTableGBTCCBSTATEClose (CCBStatetable, err);
  if (err->error) return out;
  
  /* Cleanup */
  CCBStatetable = ObitTableGBTCCBSTATEUnref(CCBStatetable);
  CCBStaterow   = ObitTableGBTCCBSTATEUnref(CCBStaterow);

  return out;
} /* end newObitGBTCCBStateInfoValue */

/**
 * Returns ClassInfo pointer for the class.
 * Initializes class if needed on first call.
 * \return pointer to the class structure.
 */
gconstpointer ObitGBTCCBStateInfoGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitGBTCCBStateInfoClassInit();

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
ObitGBTCCBStateInfo* ObitGBTCCBStateInfoCopy  (ObitGBTCCBStateInfo *in, ObitGBTCCBStateInfo *out, 
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
    out = newObitGBTCCBStateInfo(outName);
    if (outName) g_free(outName); outName = NULL;
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  out->nCCBState = in->nCCBState;
  for (i=0; i<in->nCCBState; i++) {
    out->phia[i]  = in->phia[i];
    out->phib[i]  = in->phib[i];
  }

  return out;
} /* end ObitGBTCCBStateInfoCopy */


/**
 * Initialize global ClassInfo Structure.
 */
void ObitGBTCCBStateInfoClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitGBTCCBStateInfoClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitGBTCCBStateInfoClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitGBTCCBStateInfoClassInfoDefFn (gpointer inClass)
{
  ObitGBTCCBStateInfoClassInfo *theClass = (ObitGBTCCBStateInfoClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitGBTCCBStateInfoClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitGBTCCBStateInfoClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitGBTCCBStateInfoGetClass;
  theClass->ObitClear     = (ObitClearFP)ObitGBTCCBStateInfoClear;
  theClass->ObitInit      = (ObitInitFP)ObitGBTCCBStateInfoInit;
  theClass->newObit       = (newObitFP)newObitGBTCCBStateInfo;
  theClass->ObitCopy      = (ObitCopyFP)ObitGBTCCBStateInfoCopy;
  theClass->ObitClone     = NULL;
 
} /* end ObitGBTCCBStateInfoClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Does (recursive) initialization of base class members before 
 * this class.
 * \param inn Pointer to the object to initialize.
 */
void ObitGBTCCBStateInfoInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  olong i;
  ObitGBTCCBStateInfo *in = inn;

  /* error checks */
  g_assert (in != NULL);
  
  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->nCCBState = 0;
  for (i=0; i<MAXNUMCCBSTATE; i++) {
    in->phia[i]  = 0;
    in->phib[i]  = 0;
  }

} /* end ObitGBTCCBStateInfoInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 */
void ObitGBTCCBStateInfoClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitGBTCCBStateInfo *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* free this class members */
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitGBTCCBStateInfoClear */

