/* $Id$    */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2009                                          */
/*;  Associated Universities, Inc. Washington DC, USA.                */
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
/*; Correspondence about this software should be addressed as follows:*/
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
#include "Obit.h"
#include "ObitImageSel.h"
#include "ObitMem.h"
#include "ObitTableFQUtil.h"

/*-------------- Obit: Merx mollis mortibus nuper ------------*/
/**
 * \file ObitImageSel.c
 * ObitImageSel ObitImage selector class definition.
 *
 * This contains information about data selection.
 */

/*--------------- File Global Variables  ----------------*/
/** name of the class defined in this file */
static gchar *myClassName = "ObitImageSel";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo global structure ObitIOClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitImageSelClassInfo myClassInfo = {FALSE};

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitImageSelInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitImageSelClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitImageSelClassInfoDefFn (gpointer inClass);

/*---------------Public functions---------------------------*/
/**
 * Construct Object.
 * \return pointer to object created.
 */
ObitImageSel* newObitImageSel (gchar *name)
{
  ObitImageSel* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) 
    ObitImageSelClassInit();

  /* allocate structure */
  out = ObitMemAlloc0Name(sizeof(ObitImageSel), "ObitImageSel");

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

 /* set classInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitImageSelInit((gpointer)out);

  return out;
} /* end newObitImageSel */

/**
 * Returns ClassInfo pointer for the class.
 * Initializes class if needed on first call.
 * \return pointer to the class structure.
 */
gconstpointer ObitImageSelGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) 
    ObitImageSelClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitImageSelGetClass */

/**
 * Copy constructor.
 * \param in Pointer to object to be copied.
 * \param out Pointer to object to be written.  
 *            If NULL then a new structure is created.
 * \param err ObitErr error stack
 * \return Pointer to new object.
 */
ObitImageSel* 
ObitImageSelCopy (ObitImageSel* in, ObitImageSel* out, 
		  ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  gboolean oldExist;
  gchar *outName;
  olong i;

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
    out = newObitImageSel(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /* copy */
  for (i=0; i<IM_MAXDIM; i++) {
    out->blc[i] = in->blc[i];
    out->trc[i] = in->trc[i];
  }

  return out;
} /* end ObitImageSelCopy */

/**
 * Creates/resizes image I/O buffer
 * for data transfers as described by data members.
 * \param buffer Preexisting buffer or NULL if none
 * \param desc   Pointer input descriptor.
 * \param sel    Image selector.
 * \return size in floats needed for I/O.
 */
ObitFArray* 
ObitImageSelBuffer (ObitFArray *buffer, ObitImageDesc* desc, 
		    ObitImageSel* sel)
{
  ObitFArray *out = buffer;
  olong ndim=0, naxis[2];

  /* error checks */
  if (desc==NULL) return out; 
  g_assert (ObitIsA(desc, ObitImageDescGetClass()));
  g_assert (ObitIsA(sel, &myClassInfo));

  /* make sure defaults filled in */
  ObitImageSelDefault (desc, sel);

  /* determine to size needed */  
  /* big enough for a subimaged row or plane */
  if (desc->IOsize==OBIT_IO_byRow) {
    ndim = 1;
    naxis[0] = sel->trc[0] - sel->blc[0] + 1;
  } else if (desc->IOsize==OBIT_IO_byPlane) {
    ndim = 2;
    naxis[0] = sel->trc[0] - sel->blc[0] + 1;
    naxis[1] = sel->trc[1] - sel->blc[1] + 1;
  }

  /* Create out is none exists */
  if (out==NULL) {
     out = ObitFArrayCreate ("Image Buffer", ndim, naxis);
  } else { /* resize */
    out = ObitFArrayRealloc (out, ndim, naxis);
  }

  return out;
} /* end ObitImageSelBufferSize */

/**
 * minimum inaxes is 1, 
 * blc default 1 but must be between 1 and max(trc,inaxes)
 * trc defaults to inaxes but must be between blc and inaxes.
 * Also indexes structure.
 * \param in Pointer to descriptor.
 * \param sel Image selector, blc, trc members changed if needed.
 */
void ObitImageSelDefault (ObitImageDesc* in, ObitImageSel* sel)
{
  olong i;

  /* error checks */
  g_assert (ObitIsA(in, ObitImageDescGetClass()));
  g_assert (ObitIsA(sel, &myClassInfo));

  /* minimum axis dimension = 1 */
  for (i=0; i<IM_MAXDIM; i++) in->inaxes[i] = MAX (1,in->inaxes[i]);

  /* default blc = 1 */
  for (i=0; i<IM_MAXDIM; i++) 
    if (sel->blc[i] <= 0) sel->blc[i] = 1;
  for (i=0; i<IM_MAXDIM; i++) 
    if (sel->blc[i] > in->inaxes[i]) sel->blc[i] = in->inaxes[i];

  /* default trc = inaxes */
  for (i=0; i<IM_MAXDIM; i++) 
    if (sel->trc[i] <=0) sel->trc[i] = in->inaxes[i];
  for (i=0; i<IM_MAXDIM; i++) 
    if (sel->trc[i] > in->inaxes[i]) sel->trc[i] = in->inaxes[i];

  /* Index as well */
  ObitImageDescIndex(in);

  /* Patch AIPS++ Bug */
  if (in->jlocr>=0) {
    if (in->crval[in->jlocr]<0.0) in->crval[in->jlocr] += 360.0;
  }
  if (in->obsra<0.0) in->obsra += 360.0;
} /* end ObitImageSelDefault */

/**
 * Apply selection criteria to input descriptor to derive output.
 * \param in Pointer to input descriptor.
 * \param sel Image selector, blc, trc members changed if needed.
 * \param out Pointer to output descriptor.
 * \param err Obit error stack
 */
void ObitImageSelSetDesc (ObitImageDesc* in, ObitImageSel* sel,
			  ObitImageDesc* out, ObitErr *err)
{
  olong i;
  gchar *routine = "ObitImageSelSetDesc";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, ObitImageDescGetClass()));
  g_assert (ObitIsA(sel, &myClassInfo));
  g_assert (ObitIsA(out, ObitImageDescGetClass()));

  /* make sure defaults filled in */
  ObitImageSelDefault (in, sel);

  /* copy most values */
  ObitImageDescCopy (in, out, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  out->naxis = in->naxis;
  /* apply sel->blc, sel->trc */
  for (i=0; i<in->naxis; i++) {
    /* number of axes */
    out->inaxes[i] = sel->trc[i] - sel->blc[i] + 1;
    /* reference pixel */
    out->crpix[i] = in->crpix[i] - sel->blc[i] + 1;
  }

  /* Decrement alternate reference pixel for any selection in frequency */
  if (out->jlocf>=0) {
    out->altCrpix = out->altCrpix - sel->blc[out->jlocf] + 1;
  }

} /* end ObitImageSelSetDesc */

/**
 * Correct frequency for any selection by IF
 * \param in      Pointer to input descriptor.
 * \param sel     Image selector
 * \param inImage Pointer to image (as ObitData*)
 * \param err     Obit error stack
 */
void ObitImageSelSetIF (ObitImageDesc* in, ObitImageSel* sel,
			ObitData* inImage, ObitErr *err)
{
  olong highVer, fqid=1, nif=1;
  oint  *sideBand;
  ofloat *chBandW;
  odouble *freqOff;
  ObitTableFQ *FQTab=NULL;
  ObitErr *lerr=NULL;
  gchar *routine = "ObitImageSelSetIF";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitImageDescIsA(in));
  g_assert (ObitIsA(sel, &myClassInfo));
  g_assert (ObitDataIsA(inImage));

  /* If selection by IF correct frequency */
  if ((in->jlocif>=0) && (sel->blc[in->jlocif]>1)) {
    /* See if FQ table given */
    highVer = ObitTableListGetHigh (inImage->tableList, "AIPS FQ");
    if (highVer>=1) {
      lerr = newObitErr();
      highVer = 1;
      FQTab = newObitTableFQValue (inImage->name, inImage, &highVer, 
				   OBIT_IO_ReadOnly,  nif, lerr);
      ObitTableFQGetInfo (FQTab, fqid, &nif, &freqOff, &sideBand, &chBandW, lerr);
      /* Find it? */
      if (lerr->error) 
	Obit_log_error(err, OBIT_InfoWarn, "%s: Could not find/read FQ table", routine);

      lerr = ObitErrUnref(lerr);        /* Cleanup */
      FQTab = ObitTableFQUnref(FQTab);

      /* Update frequency */
      in->crval[in->jlocf] += freqOff[sel->blc[in->jlocif]-1];
	
   } /* end FQ table exists */
  }

} /* end ObitImageSelSetIF */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitImageSelClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitImageSelClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitImageSelClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitImageSelClassInfoDefFn (gpointer inClass)
{
  ObitImageSelClassInfo *theClass = (ObitImageSelClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitImageSelClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitImageSelClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitImageSelGetClass;
  theClass->newObit       = (newObitFP)newObitImageSel;
  theClass->ObitCopy      = (ObitCopyFP)ObitImageSelCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitImageSelClear;
  theClass->ObitInit      = (ObitInitFP)ObitImageSelInit;

} /* end ObitImageSelClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Does (recursive) initialization of base class members before 
 * this class.
 * \param inn Pointer to the object to initialize.
 */
void ObitImageSelInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitImageSel *in = inn;
  olong i;

  /* error checks */
  g_assert (in != NULL);
  
  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  for (i=0; i<IM_MAXDIM; i++) {
    in->blc[i] = 1;
    in->trc[i] = 0;
  }

} /* end ObitImageSelInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 */
void ObitImageSelClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitImageSel *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* free this class members */
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);

} /* end ObitImageSelClear */



