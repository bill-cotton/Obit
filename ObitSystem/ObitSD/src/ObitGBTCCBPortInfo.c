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

#include "ObitGBTCCBPortInfo.h"
#include "ObitTableGBTCCBPORT.h"

/*-------------- Obit: Merx mollis mortibus nuper ------------*/
/**
 * \file ObitGBTCCBPortInfo.c
 * ObitGBTCCBPortInfo class function definitions.
 *
 * This is a list of associated tables.
 */

/*------------------- File Global Variables - ----------------*/
/** name of the class defined in this file */
static gchar *myClassName = "ObitGBTCCBPortInfo";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo global structure ObitGBTCCBPortInfoClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitGBTCCBPortInfoClassInfo myClassInfo = {FALSE};

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitGBTCCBPortInfoInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitGBTCCBPortInfoClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitGBTCCBPortInfoClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Basic Constructor.
 * Initializes class if needed on first call.
 * \param name  A name for the object
 * \return the new object.
 */
ObitGBTCCBPortInfo* newObitGBTCCBPortInfo (gchar* name)
{
  ObitGBTCCBPortInfo* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitGBTCCBPortInfoClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitGBTCCBPortInfo));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

 /* set classInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitGBTCCBPortInfoInit((gpointer)out);

  return out;
} /* end newObitGBTCCBPortInfo */

/**
 * Constructor from values.
 * \param name  A name for the object
 * \param disk     Obit FITS disk number
 * \param scan     Date/time scan name (e.g. "2003_05_05_05:32:56")
 * \param err      Obit error stack object.
 * \return the new object.
 */
ObitGBTCCBPortInfo* 
newObitGBTCCBPortInfoValue (gchar *name, olong disk, gchar *scan, ObitErr *err)
{
  ObitGBTCCBPortInfo* out=NULL;
  ObitTableGBTCCBPORT     *CCBPorttable=NULL;
  ObitTableGBTCCBPORTRow  *CCBPortrow=NULL;
  ObitIOCode retCode;
  gchar *tab, FullFile[128];
  olong irow;
  olong ver, nrow, iout;
  gchar *routine = "newObitGBTCCBPortInfoValue";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return out;
  g_assert (scan!=NULL);

  /* Create basic object */
  out = newObitGBTCCBPortInfo(name);

  /* Fill in values */
  /* Get CCB Port setup information from the CCBPort table */
  /* get full file name */
  sprintf (FullFile,"CCB26_40/%s.fits", scan);

  /* Create table structure */
  CCBPorttable = newObitTableGBTCCBPORT("CCBPort");

  /* Setup */
  tab = "PORT";
  ver = 1; 
  nrow = 1;
  ObitTableSetFITS(CCBPorttable,disk,FullFile,tab,ver,nrow,err);
  
  /* Open */
  retCode = ObitTableGBTCCBPORTOpen (CCBPorttable, OBIT_IO_ReadOnly, err);
  if (err->error) return out;
  
  /* Create Row structure */
  CCBPortrow = newObitTableGBTCCBPORTRow (CCBPorttable);

  /* Loop over table */
  iout = 0;
  for (irow = 1; irow<=CCBPorttable->myDesc->nrow; irow++) {
    retCode = ObitTableGBTCCBPORTReadRow (CCBPorttable, irow, CCBPortrow, err);
    if (err->error) return out;

    out->bank[iout]  =  CCBPortrow->bank;
    out->port[iout]  =  CCBPortrow->port;
    out->slave[iout] =  CCBPortrow->slave;
    iout++;
    /* Done?  */
    if (iout>=MAXNUMCCBPORT) {
      Obit_log_error(err, OBIT_Error, "%s ERROR exceed limit %d CCBPorts", 
		     routine, MAXNUMCCBPORT);
      return out;
    }
    
  } /* end loop over table */
  out->nCCBPort = iout;

  /* Close */
  retCode = ObitTableGBTCCBPORTClose (CCBPorttable, err);
  if (err->error) return out;
  
  /* Cleanup */
  CCBPorttable = ObitTableGBTCCBPORTUnref(CCBPorttable);
  CCBPortrow   = ObitTableGBTCCBPORTUnref(CCBPortrow);

  return out;
} /* end newObitGBTCCBPortInfoValue */

/**
 * Returns ClassInfo pointer for the class.
 * Initializes class if needed on first call.
 * \return pointer to the class structure.
 */
gconstpointer ObitGBTCCBPortInfoGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitGBTCCBPortInfoClassInit();

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
ObitGBTCCBPortInfo* ObitGBTCCBPortInfoCopy  (ObitGBTCCBPortInfo *in, ObitGBTCCBPortInfo *out, 
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
    out = newObitGBTCCBPortInfo(outName);
    if (outName) g_free(outName); outName = NULL;
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  out->nCCBPort = in->nCCBPort;
  for (i=0; i<in->nCCBPort; i++) {
    out->bank[i]  = in->bank[i];
    out->port[i]  = in->port[i];
    out->slave[i] = in->slave[i];
  }

  return out;
} /* end ObitGBTCCBPortInfoCopy */


/**
 * Initialize global ClassInfo Structure.
 */
void ObitGBTCCBPortInfoClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitGBTCCBPortInfoClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitGBTCCBPortInfoClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitGBTCCBPortInfoClassInfoDefFn (gpointer inClass)
{
  ObitGBTCCBPortInfoClassInfo *theClass = (ObitGBTCCBPortInfoClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitGBTCCBPortInfoClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitGBTCCBPortInfoClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitGBTCCBPortInfoGetClass;
  theClass->ObitClear     = (ObitClearFP)ObitGBTCCBPortInfoClear;
  theClass->ObitInit      = (ObitInitFP)ObitGBTCCBPortInfoInit;
  theClass->newObit       = (newObitFP)newObitGBTCCBPortInfo;
  theClass->ObitCopy      = (ObitCopyFP)ObitGBTCCBPortInfoCopy;
  theClass->ObitClone     = NULL;
 
} /* end ObitGBTCCBPortInfoClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Does (recursive) initialization of base class members before 
 * this class.
 * \param inn Pointer to the object to initialize.
 */
void ObitGBTCCBPortInfoInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  olong i;
  ObitGBTCCBPortInfo *in = inn;

  /* error checks */
  g_assert (in != NULL);
  
  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->nCCBPort = 0;
  for (i=0; i<MAXNUMCCBPORT; i++) {
    in->bank[i]  = 0;
    in->port[i]  = 0;
    in->slave[i] = 0;
  }

} /* end ObitGBTCCBPortInfoInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 */
void ObitGBTCCBPortInfoClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitGBTCCBPortInfo *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* free this class members */
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitGBTCCBPortInfoClear */

