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

#include "ObitGBTBeamOffInfo.h"
#include "ObitTableGBTBEAM_OFFSETS.h"
#include "ObitTableGBTSPDATA.h"

/*-------------- Obit: Merx mollis mortibus nuper ------------*/
/**
 * \file ObitGBTBeamOffInfo.c
 * ObitGBTBeamOffInfo class function definitions.
 *
 * This is a list of associated tables.
 */

/*------------------- File Global Variables - ----------------*/
/** name of the class defined in this file */
static gchar *myClassName = "ObitGBTBeamOffInfo";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo global structure ObitGBTBeamOffInfoClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitGBTBeamOffInfoClassInfo myClassInfo = {FALSE};

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitGBTBeamOffInfoInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitGBTBeamOffInfoClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitGBTBeamOffInfoClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Basic Constructor.
 * Initializes class if needed on first call.
 * \param name  A name for the object
 * \return the new object.
 */
ObitGBTBeamOffInfo* newObitGBTBeamOffInfo (gchar* name)
{
  ObitGBTBeamOffInfo* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitGBTBeamOffInfoClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitGBTBeamOffInfo));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

 /* set classInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitGBTBeamOffInfoInit((gpointer)out);

  return out;
} /* end newObitGBTBeamOffInfo */

/**
 * Constructor from values.
 * \param name  A name for the object
 * \param disk     Obit FITS disk number
 * \param scan     Date/time scan name (e.g. "2003_05_05_05:32:56")
 * \param err      Obit error stack object.
 * \return the new object.
 */
ObitGBTBeamOffInfo* 
newObitGBTBeamOffInfoValue (gchar *name, olong disk, gchar *scan, ObitErr *err)
{
  ObitGBTBeamOffInfo* out=NULL;
  ObitTableGBTBEAM_OFFSETS     *BeamOfftable=NULL;
  ObitTableGBTBEAM_OFFSETSRow  *BeamOffrow=NULL;
  ObitIOCode retCode;
  gchar *tab, FullFile[128];
  olong j, irow;
  olong ver, nrow, iout;
  gchar *routine = "newObitGBTBeamOffInfoValue";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return out;
  g_assert (scan!=NULL);

  /* Create basic object */
  out = newObitGBTBeamOffInfo(name);

  /* Get BeamOff setup information from the BeamOff table */
  /* get full file name */
  sprintf (FullFile,"Antenna/%s.fits", scan);

  /* Create table structure */
  BeamOfftable = newObitTableGBTBEAM_OFFSETS("BEAM_OFFSETS");

  /* Setup */
  tab = "BEAM_OFFSETS";
  ver = 1; 
  nrow = 1;
  ObitTableSetFITS(BeamOfftable,disk,FullFile,tab,ver,nrow,err);
  
  /* Open */
  retCode = ObitTableGBTBEAM_OFFSETSOpen (BeamOfftable, OBIT_IO_ReadOnly, err);
  if (err->error) return out;
  
  /* Create Row structure */
  BeamOffrow = newObitTableGBTBEAM_OFFSETSRow (BeamOfftable);

  /* Loop over table */
  iout = 0;
  for (irow = 1; irow<=BeamOfftable->myDesc->nrow; irow++) {
    retCode = ObitTableGBTBEAM_OFFSETSReadRow (BeamOfftable, irow, BeamOffrow, err);
    if (err->error) return out;

    for (j=0; j<8; j++) out->BeamName[iout][j] = BeamOffrow->Name[j];
    out->xeloff[iout]  = BeamOffrow->xeloff;
    out->eloff[iout]   = BeamOffrow->eloff;
    out->srfeed1[iout] = BeamOffrow->srfeed1;
    out->srfeed2[iout] = BeamOffrow->srfeed2;
    iout++;

    /* Done? */
    if (iout>=MAXNUMBEAMOFF) {
      Obit_log_error(err, OBIT_Error, "%s ERROR exceed limit %d BeamOffs", 
		     routine, MAXNUMBEAMOFF);
      return out;
    }
    
  } /* end loop over table */
  out->nBeamOff = iout;

  /* Close */
  retCode = ObitTableGBTBEAM_OFFSETSClose (BeamOfftable, err);
  if (err->error) return out;
  
  /* Cleanup */
  BeamOfftable = ObitTableGBTBEAM_OFFSETSUnref(BeamOfftable);
  BeamOffrow   = ObitTableGBTBEAM_OFFSETSUnref(BeamOffrow);

  return out;
} /* end newObitGBTBeamOffInfoValue */

/**
 * Returns ClassInfo pointer for the class.
 * Initializes class if needed on first call.
 * \return pointer to the class structure.
 */
gconstpointer ObitGBTBeamOffInfoGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitGBTBeamOffInfoClassInit();

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
ObitGBTBeamOffInfo* ObitGBTBeamOffInfoCopy  (ObitGBTBeamOffInfo *in, ObitGBTBeamOffInfo *out, 
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
    out = newObitGBTBeamOffInfo(outName);
    if (outName) g_free(outName); outName = NULL;
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  out->nBeamOff = in->nBeamOff;
  for (i=0; i<in->nBeamOff; i++) {
    for (j=0; j<8;j++) out->BeamName[i][j] = in->BeamName[i][j];
    out->xeloff[i]  = in->xeloff[i];
    out->eloff[i]   = in->eloff[i];
    out->srfeed1[i] = in->srfeed1[i];
    out->srfeed2[i] = in->srfeed2[i];
  }

  return out;
} /* end ObitGBTBeamOffInfoCopy */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitGBTBeamOffInfoClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitGBTBeamOffInfoClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitGBTBeamOffInfoClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitGBTBeamOffInfoClassInfoDefFn (gpointer inClass)
{
  ObitGBTBeamOffInfoClassInfo *theClass = (ObitGBTBeamOffInfoClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitGBTBeamOffInfoClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitGBTBeamOffInfoClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitGBTBeamOffInfoGetClass;
  theClass->ObitClear     = (ObitClearFP)ObitGBTBeamOffInfoClear;
  theClass->ObitInit      = (ObitInitFP)ObitGBTBeamOffInfoInit;
  theClass->newObit       = (newObitFP)newObitGBTBeamOffInfo;
  theClass->ObitCopy      = (ObitCopyFP)ObitGBTBeamOffInfoCopy;
  theClass->ObitClone     = NULL;
 
} /* end ObitGBTBeamOffInfoClassDefFn */


/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Does (recursive) initialization of base class members before 
 * this class.
 * \param inn Pointer to the object to initialize.
 */
void ObitGBTBeamOffInfoInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  olong i, j;
  ObitGBTBeamOffInfo *in = inn;

  /* error checks */
  g_assert (in != NULL);
  
  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->nBeamOff = 0;
  for (i=0; i<MAXNUMBEAMOFF; i++) {
    for (j=0; j<8;j++) in->BeamName[i][j] = ' ';
    in->xeloff[i]  = 0.0;
    in->eloff[i]   = 0.0;
    in->srfeed1[i] = 0;
    in->srfeed2[i] = 0;
  }

} /* end ObitGBTBeamOffInfoInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 */
void ObitGBTBeamOffInfoClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitGBTBeamOffInfo *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitGBTBeamOffInfoClear */

