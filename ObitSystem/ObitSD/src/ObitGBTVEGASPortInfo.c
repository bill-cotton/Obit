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

#include "ObitGBTVEGASPortInfo.h"
#include "ObitTableGBTVEGASPORT.h"

/*-------------- Obit: Merx mollis mortibus nuper ------------*/
/**
 * \file ObitGBTVEGASPortInfo.c
 * ObitGBTVEGASPortInfo class function definitions.
 *
 * This is a list of associated tables.
 */

/*------------------- File Global Variables - ----------------*/
/** name of the class defined in this file */
static gchar *myClassName = "ObitGBTVEGASPortInfo";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo global structure ObitGBTVEGASPortInfoClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitGBTVEGASPortInfoClassInfo myClassInfo = {FALSE};

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitGBTVEGASPortInfoInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitGBTVEGASPortInfoClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitGBTVEGASPortInfoClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Basic Constructor.
 * Initializes class if needed on first call.
 * \param name  A name for the object
 * \return the new object.
 */
ObitGBTVEGASPortInfo* newObitGBTVEGASPortInfo (gchar* name)
{
  ObitGBTVEGASPortInfo* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitGBTVEGASPortInfoClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitGBTVEGASPortInfo));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

 /* set classInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitGBTVEGASPortInfoInit((gpointer)out);

  return out;
} /* end newObitGBTVEGASPortInfo */

/**
 * Constructor from values.
 * \param name  A name for the object
 * \param DataRoot Root of data directory.
 * \param scan     Date/time scan name (e.g. "2003_05_05_05:32:56")
 * \param err      Obit error stack object.
 * \return the new object.
 */
ObitGBTVEGASPortInfo* 
newObitGBTVEGASPortInfoValue (gchar *name, gchar *DataRoot, gchar *scan, ObitErr *err)
{
  ObitGBTVEGASPortInfo* out=NULL;
  ObitTableGBTVEGASPORT     *VEGASPorttable=NULL;
  ObitTableGBTVEGASPORTRow  *VEGASPortrow=NULL;
  ObitIOCode retCode;
  gchar *tab, FullFile[128];
  olong irow, disk=0;
  olong ver, nrow, iout;
  gchar *routine = "newObitGBTVEGASPortInfoValue";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return out;
  g_assert (scan!=NULL);

  /* Create basic object */
  out = newObitGBTVEGASPortInfo(name);

  /* Fill in values */
  /* Get VEGAS Port setup information from the VEGASPort table */
  /* get full file name */
  sprintf (FullFile,"%s/VEGAS/%sA.fits", DataRoot, scan);

  /* Create table structure */
  VEGASPorttable = newObitTableGBTVEGASPORT("VEGASPort");

  /* Setup */
  tab = "PORT";
  ver = 1; 
  nrow = 1;
  ObitTableSetFITS(VEGASPorttable,disk,FullFile,tab,ver,nrow,err);
  
  /* Open */
  retCode = ObitTableGBTVEGASPORTOpen (VEGASPorttable, OBIT_IO_ReadOnly, err);
  if (err->error) return out;
  
  /* Create Row structure */
  VEGASPortrow = newObitTableGBTVEGASPORTRow (VEGASPorttable);

  /* Loop over table */
  iout = 0;
  for (irow = 1; irow<=VEGASPorttable->myDesc->nrow; irow++) {
    retCode = ObitTableGBTVEGASPORTReadRow (VEGASPorttable, irow, VEGASPortrow, err);
    if (err->error) return out;

    out->bank[iout]    =  VEGASPortrow->bank;
    out->port[iout]    =  VEGASPortrow->port;
    out->measpwr[iout] =  VEGASPortrow->measpwr;
    out->tone[iout] = g_strndup (VEGASPortrow->tone, 5);
    iout++;
    /* Done?  */
    if (iout>=MAXNUMVEGASPORT) {
      Obit_log_error(err, OBIT_Error, "%s ERROR exceed limit %d VEGASPorts", 
		     routine, MAXNUMVEGASPORT);
      return out;
    }
    
  } /* end loop over table */
  out->nVEGASPort = iout;

  /* Close */
  retCode = ObitTableGBTVEGASPORTClose (VEGASPorttable, err);
  if (err->error) return out;
  
  /* Cleanup */
  VEGASPorttable = ObitTableGBTVEGASPORTUnref(VEGASPorttable);
  VEGASPortrow   = ObitTableGBTVEGASPORTUnref(VEGASPortrow);

  return out;
} /* end newObitGBTVEGASPortInfoValue */

/**
 * Returns ClassInfo pointer for the class.
 * Initializes class if needed on first call.
 * \return pointer to the class structure.
 */
gconstpointer ObitGBTVEGASPortInfoGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitGBTVEGASPortInfoClassInit();

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
ObitGBTVEGASPortInfo* ObitGBTVEGASPortInfoCopy  (ObitGBTVEGASPortInfo *in, ObitGBTVEGASPortInfo *out, 
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
    out = newObitGBTVEGASPortInfo(outName);
    if (outName) g_free(outName); outName = NULL;
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  out->nVEGASPort = in->nVEGASPort;
  for (i=0; i<in->nVEGASPort; i++) {
    out->bank[i]  = in->bank[i];
    out->port[i]  = in->port[i];
    out->measpwr[i]  = in->measpwr[i];
    out->tone[i] = g_strndup (in->tone[i], 5);
  }

  return out;
} /* end ObitGBTVEGASPortInfoCopy */


/**
 * Initialize global ClassInfo Structure.
 */
void ObitGBTVEGASPortInfoClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitGBTVEGASPortInfoClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitGBTVEGASPortInfoClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitGBTVEGASPortInfoClassInfoDefFn (gpointer inClass)
{
  ObitGBTVEGASPortInfoClassInfo *theClass = (ObitGBTVEGASPortInfoClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitGBTVEGASPortInfoClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitGBTVEGASPortInfoClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitGBTVEGASPortInfoGetClass;
  theClass->ObitClear     = (ObitClearFP)ObitGBTVEGASPortInfoClear;
  theClass->ObitInit      = (ObitInitFP)ObitGBTVEGASPortInfoInit;
  theClass->newObit       = (newObitFP)newObitGBTVEGASPortInfo;
  theClass->ObitCopy      = (ObitCopyFP)ObitGBTVEGASPortInfoCopy;
  theClass->ObitClone     = NULL;
 
} /* end ObitGBTVEGASPortInfoClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Does (recursive) initialization of base class members before 
 * this class.
 * \param inn Pointer to the object to initialize.
 */
void ObitGBTVEGASPortInfoInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  olong i;
  ObitGBTVEGASPortInfo *in = inn;

  /* error checks */
  g_assert (in != NULL);
  
  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->nVEGASPort = 0;
  for (i=0; i<MAXNUMVEGASPORT; i++) {
    in->bank[i]    = 0;
    in->port[i]    = 0;
    in->measpwr[i] = 0.0;
    in->tone[i]    = NULL;
  }

} /* end ObitGBTVEGASPortInfoInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 */
void ObitGBTVEGASPortInfoClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitGBTVEGASPortInfo *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* free this class members */
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitGBTVEGASPortInfoClear */

