/* $Id$   */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2008                                          */
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

#include "ObitGBTIFInfo.h"
#include "ObitTableGBTIF.h"
#include "ObitTableGBTSPDATA.h"

/*-------------- Obit: Merx mollis mortibus nuper ------------*/
/**
 * \file ObitGBTIFInfo.c
 * ObitGBTIFInfo class function definitions.
 *
 * This is a list of associated tables.
 */

/*------------------- File Global Variables - ----------------*/
/** name of the class defined in this file */
static gchar *myClassName = "ObitGBTIFInfo";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo global structure ObitGBTIFInfoClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitGBTIFInfoClassInfo myClassInfo = {FALSE};

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitGBTIFInfoInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitGBTIFInfoClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitGBTIFInfoClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Basic Constructor.
 * Initializes class if needed on first call.
 * \param name  A name for the object
 * \return the new object.
 */
ObitGBTIFInfo* newObitGBTIFInfo (gchar* name)
{
  ObitGBTIFInfo* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitGBTIFInfoClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitGBTIFInfo));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

 /* set classInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitGBTIFInfoInit((gpointer)out);

  return out;
} /* end newObitGBTIFInfo */

/**
 * Constructor from values.
 * \param name  A name for the object
 * \param backend  Name of the backend for which the information is desired
 *                 (e.g. "DCR", "SpectralProcessor")
 * \param disk     Obit FITS disk number
 * \param scan     Date/time scan name (e.g. "2003_05_05_05:32:56")
 * \param err      Obit error stack object.
 * \return the new object.
 */
ObitGBTIFInfo* 
newObitGBTIFInfoValue (gchar *name, gchar *backend, olong disk, gchar *scan, ObitErr *err)
{
  ObitGBTIFInfo* out=NULL;
  ObitTableGBTSPDATA *DATAtable=NULL;
  ObitTableGBTIF     *IFtable=NULL;
  ObitTableGBTIFRow  *IFrow=NULL;
  ObitIOCode retCode;
  gchar *tab, FullFile[128];
  gboolean isCCB=FALSE;
  olong irow;
  olong ver, nrow, iout;
  gchar *routine = "newObitGBTIFInfoValue";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return out;
  g_assert (backend!=NULL);
  g_assert (scan!=NULL);

  /* Create basic object */
  out = newObitGBTIFInfo(name);

  /* Fill in values */
  out->Backend = g_strdup(backend);

  /* if the Spectral Processor, get number of channels from the 
     SPDATA first dimension */
  if (!strncmp (backend, "SpectralProcessor", 17)) {

    /* get full file name */
    sprintf (FullFile,"SpectralProcessor/%s.fits", scan);

    /* Create table structure */
    DATAtable = newObitTableGBTSPDATA("Data");

    /* Setup */
    tab = "DATA";
    ver = 1; 
    nrow = 1;
    ObitTableSetFITS(DATAtable,disk,FullFile,tab,ver,nrow,err);

    /* Open */
    retCode = ObitTableGBTSPDATAOpen (DATAtable, OBIT_IO_ReadOnly, err);
    if (err->error) return out;

    /* Get number of channels = first dimension of Data column */
    out->nchan = DATAtable->myDesc->dim[DATAtable->dataCol][0];
    
    /* Close */
    retCode = ObitTableGBTSPDATAClose (DATAtable, err);
    if (err->error) return out;
    
    /* Cleanup */
    DATAtable = ObitTableGBTSPDATAUnref(DATAtable);

    /* If CCB nchan = 4 */
  } else if (!strncmp (backend, "CCB26_40", 8)) { 
    isCCB = TRUE;
    out->nchan = 4;
  } else { /* Something else, only one channel */
    out->nchan = 1;
  } /* End of getting number of channels */

  /* Get IF setup information from the IF table */
  /* get full file name */
  sprintf (FullFile,"IF/%s.fits", scan);

  /* Create table structure */
  IFtable = newObitTableGBTIF("IF");

  /* Setup */
  tab = "IF";
  ver = 1; 
  nrow = 1;
  ObitTableSetFITS(IFtable,disk,FullFile,tab,ver,nrow,err);
  
  /* Open */
  retCode = ObitTableGBTIFOpen (IFtable, OBIT_IO_ReadOnly, err);
  if (err->error) return out;
  
  /* Create Row structure */
  IFrow = newObitTableGBTIFRow (IFtable);

  /* Loop over table */
  iout = 0;
  for (irow = 1; irow<=IFtable->myDesc->nrow; irow++) {
    retCode = ObitTableGBTIFReadRow (IFtable, irow, IFrow, err);
    if (err->error) return out;

    /* Want this one? */
    if (!strncmp (backend, IFrow->backend, strlen(backend))) {
      out->poln[iout]  = IFrow->polarize[0];
      if (isCCB) out->delta[iout] = IFrow->bandwdth;
      else out->delta[iout] = IFrow->bandwdth/out->nchan;
      /* Lower sideband? */
      if (IFrow->sideband[0]=='L') out->delta[iout] = -out->delta[iout];
      out->refPixel[iout]     = 1;
      out->refFrequency[iout] = IFrow->CenterSky;
      out->bank[iout][0] = IFrow->bank[0];
      out->bank[iout][1] = IFrow->bank[1];
      out->port[iout]    = IFrow->port;
      out->feed[iout]    = IFrow->feed;
      out->srfeed1[iout] = IFrow->srfeed1;
      out->srfeed2[iout] = IFrow->srfeed2;
      iout++;
    }
    /* Done? */
    if (iout>=MAXNUMIF) {
      Obit_log_error(err, OBIT_Error, "%s ERROR exceep limit %d IFs", routine, MAXNUMIF);
      return out;
    }
    
  } /* end loop over table */
  out->nIF = iout;

  /* Close */
  retCode = ObitTableGBTIFClose (IFtable, err);
  if (err->error) return out;
  
  /* Cleanup */
  IFtable = ObitTableGBTIFUnref(IFtable);
  IFrow   = ObitTableGBTIFUnref(IFrow);

  return out;
} /* end newObitGBTIFInfoValue */

/**
 * Returns ClassInfo pointer for the class.
 * Initializes class if needed on first call.
 * \return pointer to the class structure.
 */
gconstpointer ObitGBTIFInfoGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitGBTIFInfoClassInit();

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
ObitGBTIFInfo* ObitGBTIFInfoCopy  (ObitGBTIFInfo *in, ObitGBTIFInfo *out, 
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
    out = newObitGBTIFInfo(outName);
    if (outName) g_free(outName); outName = NULL;
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  out->nchan = in->nchan;
  out->nIF   = in->nIF;
  for (i=0; i<in->nIF; i++) {
    out->poln[i]         = in->poln[i];
    out->delta[i]        = in->delta[i];
    out->refPixel[i]     = in->refPixel[i];
    out->refFrequency[i] = in->refFrequency[i];
    out->feed[i]         = in->feed[i];
    out->srfeed1[i]      = in->srfeed1[i];
    out->srfeed2[i]      = in->srfeed2[i];
  }

  if (out->Backend) g_free(out->Backend);
  out->Backend = g_strdup(in->Backend);

  return out;
} /* end ObitGBTIFInfoCopy */


/**
 * Initialize global ClassInfo Structure.
 */
void ObitGBTIFInfoClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitGBTIFInfoClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitGBTIFInfoClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitGBTIFInfoClassInfoDefFn (gpointer inClass)
{
  ObitGBTIFInfoClassInfo *theClass = (ObitGBTIFInfoClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitGBTIFInfoClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitGBTIFInfoClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitGBTIFInfoGetClass;
  theClass->ObitClear     = (ObitClearFP)ObitGBTIFInfoClear;
  theClass->ObitInit      = (ObitInitFP)ObitGBTIFInfoInit;
  theClass->newObit       = (newObitFP)newObitGBTIFInfo;
  theClass->ObitCopy      = (ObitCopyFP)ObitGBTIFInfoCopy;
  theClass->ObitClone     = NULL;

} /* end ObitGBTIFInfoClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Does (recursive) initialization of base class members before 
 * this class.
 * \param inn Pointer to the object to initialize.
 */
void ObitGBTIFInfoInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  olong i;
  ObitGBTIFInfo *in = inn;

  /* error checks */
  g_assert (in != NULL);
  
  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->Backend = NULL;
  in->nIF = 0;
  for (i=0; i<MAXNUMIF; i++) {
    in->poln[i]         = ' ';
    in->delta[i]        = 0.0;
    in->refPixel[i]     = 0.0;
    in->refFrequency[i] = 0.0;
    in->feed[i]         = 0;
    in->srfeed1[i]      = 0;
    in->srfeed2[i]      = 0;
  }

} /* end ObitGBTIFInfoInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 */
void ObitGBTIFInfoClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitGBTIFInfo *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* free this class members */
  if (in->Backend) g_free(in->Backend);
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitGBTIFInfoClear */

