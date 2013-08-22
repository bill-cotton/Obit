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

#include "ObitGBTVEGASSamplerInfo.h"
#include "ObitTableGBTVEGASSAMPLER.h"

/*-------------- Obit: Merx mollis mortibus nuper ------------*/
/**
 * \file ObitGBTVEGASSamplerInfo.c
 * ObitGBTVEGASSamplerInfo class function definitions.
 *
 * This is a list of associated tables.
 */

/*------------------- File Global Variables - ----------------*/
/** name of the class defined in this file */
static gchar *myClassName = "ObitGBTVEGASSamplerInfo";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo global structure ObitGBTVEGASSamplerInfoClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitGBTVEGASSamplerInfoClassInfo myClassInfo = {FALSE};

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitGBTVEGASSamplerInfoInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitGBTVEGASSamplerInfoClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitGBTVEGASSamplerInfoClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Basic Constructor.
 * Initializes class if needed on first call.
 * \param name  A name for the object
 * \return the new object.
 */
ObitGBTVEGASSamplerInfo* newObitGBTVEGASSamplerInfo (gchar* name)
{
  ObitGBTVEGASSamplerInfo* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitGBTVEGASSamplerInfoClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitGBTVEGASSamplerInfo));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

 /* set classInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitGBTVEGASSamplerInfoInit((gpointer)out);

  return out;
} /* end newObitGBTVEGASSamplerInfo */

/**
 * Constructor from values.
 * \param name  A name for the object
 * \param DataRoot Root of data directory.
 * \param scan     Date/time scan name (e.g. "2003_05_05_05:32:56")
 * \param err      Obit error stack object.
 * \return the new object.
 */
ObitGBTVEGASSamplerInfo* 
newObitGBTVEGASSamplerInfoValue (gchar *name, gchar *DataRoot, gchar *scan, ObitErr *err)
{
  ObitGBTVEGASSamplerInfo* out=NULL;
  ObitTableGBTVEGASSAMPLER     *VEGASSamplertable=NULL;
  ObitTableGBTVEGASSAMPLERRow  *VEGASSamplerrow=NULL;
  ObitIOCode retCode;
  gchar *tab, FullFile[128];
  olong i, irow, disk=0;
  olong ver, nrow, iout;
  gchar *routine = "newObitGBTVEGASSamplerInfoValue";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return out;
  g_assert (scan!=NULL);

  /* Create basic object */
  out = newObitGBTVEGASSamplerInfo(name);

  /* Fill in values */
  /* Get VEGAS Sampler setup information from the VEGASSampler table */
  /* get full file name */
  sprintf (FullFile,"%s/VEGAS/%sA.fits", DataRoot, scan);

  /* Create table structure */
  VEGASSamplertable = newObitTableGBTVEGASSAMPLER("VEGASSampler");

  /* Setup */
  tab = "SAMPLER";
  ver = 1; 
  nrow = 1;
  ObitTableSetFITS(VEGASSamplertable,disk,FullFile,tab,ver,nrow,err);
  
  /* Open */
  retCode = ObitTableGBTVEGASSAMPLEROpen (VEGASSamplertable, OBIT_IO_ReadOnly, err);
  if (err->error) return out;
  
  /* Create Row structure */
  VEGASSamplerrow = newObitTableGBTVEGASSAMPLERRow (VEGASSamplertable);

  /* Loop over table */
  iout = 0;
  for (irow = 1; irow<=VEGASSamplertable->myDesc->nrow; irow++) {
    retCode = ObitTableGBTVEGASSAMPLERReadRow (VEGASSamplertable, irow, VEGASSamplerrow, err);
    if (err->error) return out;

    out->bank_a[iout]   = VEGASSamplerrow->bank_a;
    out->port_a[iout]   = VEGASSamplerrow->port_a;
    out->bank_b[iout]   = VEGASSamplerrow->bank_b;
    out->port_b[iout]   = VEGASSamplerrow->port_b;
    out->subband[iout]  = VEGASSamplerrow->subband;
    out->crval1[iout]   = VEGASSamplerrow->crval1;
    out->cdelt1[iout]   = VEGASSamplerrow->cdelt1;
    out->freqres[iout]  = VEGASSamplerrow->freqres;
    for (i=0; i<4; i++) out->datatype[iout][i] = VEGASSamplerrow->datatype[i];
    iout++;
    /* Done?  */
    if (iout>=MAXNUMVEGASSAMPLER) {
      Obit_log_error(err, OBIT_Error, "%s ERROR exceed limit %d VEGASSamplers", 
		     routine, MAXNUMVEGASSAMPLER);
      return out;
    }
    
  } /* end loop over table */
  out->nVEGASSampler = iout;

  /* Close */
  retCode = ObitTableGBTVEGASSAMPLERClose (VEGASSamplertable, err);
  if (err->error) return out;
  
  /* Cleanup */
  VEGASSamplertable = ObitTableGBTVEGASSAMPLERUnref(VEGASSamplertable);
  VEGASSamplerrow   = ObitTableGBTVEGASSAMPLERUnref(VEGASSamplerrow);

  return out;
} /* end newObitGBTVEGASSamplerInfoValue */

/**
 * Returns ClassInfo pointer for the class.
 * Initializes class if needed on first call.
 * \return pointer to the class structure.
 */
gconstpointer ObitGBTVEGASSamplerInfoGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitGBTVEGASSamplerInfoClassInit();

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
ObitGBTVEGASSamplerInfo* ObitGBTVEGASSamplerInfoCopy  (ObitGBTVEGASSamplerInfo *in, ObitGBTVEGASSamplerInfo *out, 
				   ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  gboolean oldExist;
  olong i, j;
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
    out = newObitGBTVEGASSamplerInfo(outName);
    if (outName) g_free(outName); outName = NULL;
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  out->nVEGASSampler = in->nVEGASSampler;
  for (i=0; i<in->nVEGASSampler; i++) {
    out->bank_a[i]   = in->bank_a[i];
    out->port_a[i]   = in->port_a[i];
    out->bank_b[i]   = in->bank_b[i];
    out->port_b[i]   = in->port_b[i];
    out->subband[i]  = in->subband[i];
    out->crval1[i]   = in->crval1[i];
    out->cdelt1[i]   = in->cdelt1[i];
    out->freqres[i]  = in->freqres[i];
    for (j=0; j<4; j++) out->datatype[i][j] = in->datatype[i][j];
  }

  return out;
} /* end ObitGBTVEGASSamplerInfoCopy */


/**
 * Initialize global ClassInfo Structure.
 */
void ObitGBTVEGASSamplerInfoClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitGBTVEGASSamplerInfoClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitGBTVEGASSamplerInfoClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitGBTVEGASSamplerInfoClassInfoDefFn (gpointer inClass)
{
  ObitGBTVEGASSamplerInfoClassInfo *theClass = (ObitGBTVEGASSamplerInfoClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitGBTVEGASSamplerInfoClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitGBTVEGASSamplerInfoClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitGBTVEGASSamplerInfoGetClass;
  theClass->ObitClear     = (ObitClearFP)ObitGBTVEGASSamplerInfoClear;
  theClass->ObitInit      = (ObitInitFP)ObitGBTVEGASSamplerInfoInit;
  theClass->newObit       = (newObitFP)newObitGBTVEGASSamplerInfo;
  theClass->ObitCopy      = (ObitCopyFP)ObitGBTVEGASSamplerInfoCopy;
  theClass->ObitClone     = NULL;
 
} /* end ObitGBTVEGASSamplerInfoClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Does (recursive) initialization of base class members before 
 * this class.
 * \param inn Pointer to the object to initialize.
 */
void ObitGBTVEGASSamplerInfoInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  olong i, j;
  ObitGBTVEGASSamplerInfo *in = inn;

  /* error checks */
  g_assert (in != NULL);
  
  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->nVEGASSampler = 0;
  for (i=0; i<MAXNUMVEGASSAMPLER; i++) {
    in->bank_a[i]   = ' ';
    in->port_a[i]   = ' ';
    in->bank_b[i]   = ' ';
    in->port_b[i]   = ' ';
    in->subband[i]  = 0;
    in->crval1[i]   = 0.0;
    in->cdelt1[i]   = 0.0;
    in->freqres[i]  = 0.0;
    for (j=0; j<4; j++) in->datatype[i][j] = ' ';
  }

} /* end ObitGBTVEGASSamplerInfoInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 */
void ObitGBTVEGASSamplerInfoClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitGBTVEGASSamplerInfo *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* free this class members */
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitGBTVEGASSamplerInfoClear */

