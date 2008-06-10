/* $Id: ObitOTFSkyModel.c,v 1.5 2005/12/19 00:37:10 bcotton Exp $ */
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

#include "ObitOTFSkyModel.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitOTFSkyModel.c
 * GBT/OTF Sky model class function definitions.
 * This class is derived from the Obit base class.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitOTFSkyModel";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/*--------------- File Global Variables  ----------------*/
/**
 * ClassInfo structure ObitOTFSkyModelClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitOTFSkyModelClassInfo myClassInfo = {FALSE};

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitOTFSkyModelInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitOTFSkyModelClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitOTFSkyModelClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitOTFSkyModel* newObitOTFSkyModel (gchar* name)
{
  ObitOTFSkyModel* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitOTFSkyModelClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitOTFSkyModel));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitOTFSkyModelInit((gpointer)out);

 return out;
} /* end newObitOTFSkyModel */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitOTFSkyModelGetClass (void)
{
  return (gconstpointer)&myClassInfo;
} /* end ObitOTFSkyModelGetClass */

/**
 * Make a deep copy of an ObitOTFSkyModel.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitOTFSkyModel* 
ObitOTFSkyModelCopy  (ObitOTFSkyModel *in, ObitOTFSkyModel *out, ObitErr *err)
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
    out = newObitOTFSkyModel(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  out->info = ObitInfoListUnref(out->info);
  out->info = ObitInfoListCopy(in->info);

  /* this class data */
  out->numberComp = in->numberComp;
  out->RAOffset   = g_realloc(out->RAOffset, in->numberComp*sizeof(ofloat));
  out->DecOffset  = g_realloc(out->DecOffset, in->numberComp*sizeof(ofloat));
  out->flux       = g_realloc(out->flux, in->numberComp*sizeof(ofloat));
  for (i=0; i<in->numberComp; i++) out->RAOffset[i] = in->RAOffset[i];
  for (i=0; i<in->numberComp; i++) out->DecOffset[i] = in->DecOffset[i];
  for (i=0; i<in->numberComp; i++) out->flux[i] = in->flux[i];
  out->RACenter = in->RACenter;
  out->DecCenter = in->DecCenter;
  out->proj      = in->proj;
  return out;
} /* end ObitOTFSkyModelCopy */

/**
 * Creates an ObitOTFSkyModel for a specified number of components.
 * \param ncomp  Number of components
 * \return the new object.
 */
ObitOTFSkyModel* ObitOTFSkyModelCreate (olong ncomp)
{
  ObitOTFSkyModel* out;

  /* Create basic structure */
  out = newObitOTFSkyModel ("Sky Model");

  /* set members */
  out->numberComp = ncomp;
  out->RAOffset  = g_realloc(out->RAOffset,  out->numberComp*sizeof(ofloat));
  out->DecOffset = g_realloc(out->DecOffset, out->numberComp*sizeof(ofloat));
  out->flux      = g_realloc(out->flux, out->numberComp*sizeof(ofloat));

  return out;
} /* end ObitOTFSkyModelCreate */

/**
 * Read ObitOTFSkyModel information from a table
 * \param in    Sky model to update, if Null, it is created
 * \param table table to read from
 * \param err   Error stack
 * \return return code OBIT_IO_OK => OK
 */
ObitIOCode ObitOTFSkyModelRead (ObitOTFSkyModel **in, ObitTableSkyModel *table, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitOTFSkyModel *out=NULL;
  ObitTableSkyModelRow *row=NULL;
  olong i, numberComp=0;
  gchar *routine = "ObitOTFSkyModelRead";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitTableSkyModelIsA(table));

  /* Open table */
  retCode = ObitTableSkyModelOpen (table, OBIT_IO_ReadWrite, err);
    if ((retCode != OBIT_IO_OK) || (err->error)) 
      Obit_traceback_val (err, routine, table->name, retCode);

  /* Does input exist? */
  if (in!=NULL) {
    out = *in;  /* exists */
    numberComp = (*in)->numberComp;
  } else {
    /* must create */
    numberComp = table->myDesc->nrow;
    out = ObitOTFSkyModelCreate (numberComp);
  }

  /* Allocate arrays */
  out->RAOffset  = g_realloc(out->RAOffset,  numberComp*sizeof(ofloat));
  out->DecOffset = g_realloc(out->DecOffset, numberComp*sizeof(ofloat));
  out->flux      = g_realloc(out->flux, numberComp*sizeof(ofloat));

  /* Copy header information */
  out->RACenter  = table->RA;
  out->DecCenter = table->Dec;
  out->proj      = ObitOTFSkyModelProj (table->Proj);

  /* Create row object */
  row = newObitTableSkyModelRow (table);

  /* Loop over table */
  for (i=0; i<table->myDesc->nrow; i++) {
    /* Read table row */
    retCode = ObitTableSkyModelReadRow (table, i+1, row, err);
    if ((retCode != OBIT_IO_OK) || (err->error)) 
      Obit_traceback_val (err, routine, table->name, retCode);

    /* Save data */
    out->RAOffset[i]  = row->RAOff;
    out->DecOffset[i] = row->DecOff;
    out->flux[i]      = row->Flux;
 } /* end loop over rows */

  /* Release row object */
  row = ObitTableSkyModelRowUnref (row);

  /* Close table */
  retCode = ObitTableSkyModelClose (table, err);
  if ((retCode != OBIT_IO_OK) || (err->error)) 
    Obit_traceback_val (err, routine, table->name, retCode);
  
/* set pointer to array geometry if created here */
  if (in==NULL) in = &out;

  return retCode;
} /* end ObitOTFSkyModelRead */

/**
 * Write ObitOTFSkyModel information to a table
 * \param in    Array geometry to write
 * \param table table to write to
 * \param err   Error stack
 * \return return code OBIT_IO_OK => OK
 */
ObitIOCode ObitOTFSkyModelWrite (ObitOTFSkyModel *in, ObitTableSkyModel *table, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  olong i;
  ObitTableSkyModelRow *row=NULL;
  gchar *routine = "ObitOTFSkyModelWrite";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitTableSkyModelIsA(table));

  /* Open table */
  retCode = ObitTableSkyModelOpen (table, OBIT_IO_WriteOnly, err);
  if ((retCode != OBIT_IO_OK) || (err->error)) 
    Obit_traceback_val (err, routine, table->name, retCode);

  /* Copy header information */
  table->RA  = in->RACenter;
  table->Dec = in->DecCenter;
  if (in->proj == OBIT_OTF_SIN) strncpy (table->Proj, "-SIN", 4);
  if (in->proj == OBIT_OTF_ARC) strncpy (table->Proj, "-ARC", 4);
  if (in->proj == OBIT_OTF_TAN) strncpy (table->Proj, "-TAN", 4);

  /* Create row object */
  row = newObitTableSkyModelRow (table);

  /* Loop over table */
  for (i=0; i<in->numberComp; i++) {
    /* Save data */
    row->RAOff  = in->RAOffset[i];
    row->DecOff = in->DecOffset[i];
    row->Flux   = in->flux[i];

    /* Write table row */
    retCode = ObitTableSkyModelWriteRow (table, i+1, row, err);
    if ((retCode != OBIT_IO_OK) || (err->error)) 
      Obit_traceback_val (err, routine, table->name, retCode);

  } /* end loop over rows */

  /* Release row object */
  row = ObitTableSkyModelRowUnref (row);

  /* Close table */
  retCode = ObitTableSkyModelClose (table, err);
  if ((retCode != OBIT_IO_OK) || (err->error)) 
    Obit_traceback_val (err, routine, table->name, retCode);
  
  return retCode;
} /* end ObitOTFSkyModelWrite */

/**
 * Determine Projection type code
 * Recognizes '-SIN', '-ARC', '-TAN'
 * default is -SIN.
 * \param string string code to test.
 * \param table table to write to
 * \param err   Error stack
 * \return return code OBIT_IO_OK => OK
 */
ObitOTFProj ObitOTFSkyModelProj (gchar *string)
{
  ObitOTFProj out;
 
  /* error checks */
  g_assert (string!=NULL);

  out  = OBIT_OTF_SIN; /* Default projection */
  if (!strncmp (string, "-ARC", 4)) out = OBIT_OTF_ARC;
  if (!strncmp (string, "-TAN", 4)) out = OBIT_OTF_TAN;

  return out;
}

/**
 * Initialize global ClassInfo Structure.
 */
void ObitOTFSkyModelClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitOTFSkyModelClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitOTFSkyModelClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitOTFSkyModelClassInfoDefFn (gpointer inClass)
{
  ObitOTFSkyModelClassInfo *theClass = (ObitOTFSkyModelClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitOTFSkyModelClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitOTFSkyModelClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitOTFSkyModelGetClass;
  theClass->ObitClear     = (ObitClearFP)ObitOTFSkyModelClear;
  theClass->ObitInit      = (ObitInitFP)ObitOTFSkyModelInit;
  theClass->newObit       = (newObitFP)newObitOTFSkyModel;
  theClass->ObitCopy      = (ObitCopyFP)ObitOTFSkyModelCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitOTFSkyModelCreate = (ObitOTFSkyModelCreateFP)ObitOTFSkyModelCreate;

} /* end ObitOTFSkyModelClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitOTFSkyModelInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitOTFSkyModel *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->thread       = newObitThread();
  in->info         = newObitInfoList(); 
  in->RAOffset     = NULL;
  in->DecOffset    = NULL;
  in->flux         = NULL;
  in->numberComp = 0;

} /* end ObitOTFSkyModelInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * For some reason this wasn't build into the GType class.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitOTFSkyModel* cast to an Obit*.
 */
void ObitOTFSkyModelClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitOTFSkyModel *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  in->thread    = ObitThreadUnref(in->thread);
  in->info      = ObitInfoListUnref(in->info);
  if (in->RAOffset)  g_free(in->RAOffset);  in->RAOffset = NULL;
  if (in->DecOffset) g_free(in->DecOffset);  in->DecOffset = NULL;
  if (in->flux)      g_free(in->flux);  in->flux = NULL;
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitOTFSkyModelClear */


