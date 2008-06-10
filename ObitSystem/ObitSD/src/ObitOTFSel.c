/* $Id: ObitOTFSel.c,v 1.11 2008/02/28 15:22:01 bcotton Exp $      */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2008                                          */
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
/*;Correspondence about this software should be addressed as follows: */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
#include "Obit.h"
#include "ObitOTFSel.h"
#include "ObitTableOTFIndex.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitOTFSel.c
 * ObitOTFSel Obit GBT/OTF data selector class definition.
 * This contains information about data selection and calibration.
 */

/*--------------- File Global Variables  ----------------*/
/** name of the class defined in this file */
static gchar *myClassName = "ObitOTFSel";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo global structure ObitIOClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitOTFSelClassInfo myClassInfo = {FALSE};

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitOTFSelInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitOTFSelClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitOTFSelClassInfoDefFn (gpointer inClass);

/*---------------Public functions---------------------------*/
/**
 * Construct Object.
 * \return pointer to object created.
 */
ObitOTFSel* newObitOTFSel (gchar *name)
{
  ObitOTFSel* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitOTFSelClassInit();

  /* allocate structure */
  out = g_malloc0(sizeof(ObitOTFSel));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

 /* set classInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitOTFSelInit((gpointer)out);

  return out;
} /* end newObitOTFSel */

/**
 * Returns ClassInfo pointer for the class.
 * Initializes class if needed on first call.
 * \return pointer to the class structure.
 */
gconstpointer ObitOTFSelGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) 
    ObitOTFSelClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitOTFSelGetClass */

/**
 * Copy constructor.
 * \param in Pointer to object to be copied.
 * \param out Pointer to object to be written.  
 *            If NULL then a new structure is created.
 * \param err ObitErr error stack
 * \return Pointer to new object.
 */
ObitOTFSel* ObitOTFSelCopy (ObitOTFSel* in, ObitOTFSel* out, 
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
    out = newObitOTFSel(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /* This class members */
  out->FileType     = in->FileType;
  out->nRecPIO      = in->nRecPIO;
  out->doCalSelect  = in->doCalSelect;
  out->calVersion   = in->calVersion;
  out->doCal        = in->doCal;
  out->doIndex      = in->doIndex;
  out->transPol     = in->transPol;
  out->doFlag       = in->doFlag;
  out->FGversion    = in->FGversion;
  out->numRow       = in->numRow;
  out->LastRowRead  = in->LastRowRead;
  out->scanFirstRec = in->scanFirstRec;
  out->scanLastRec  = in->scanLastRec;
  out->numberPoln   = in->numberPoln;
  out->jincs        = in->jincs;
  out->startChann   = in->startChann;
  out->numberChann  = in->numberChann;
  out->jincf        = in->jincf;
  out->numberFeed   = in->numberFeed;
  out->jincfeed     = in->jincfeed;
  out->timeRange[0] = in->timeRange[0];
  out->timeRange[1] = in->timeRange[1];
  out->scans[0]     = in->scans[0];
  out->scans[1]     = in->scans[1];
  for (i=0; i<5; i++) out->Stokes[i] = in->Stokes[i];
  out->keepCal      = in->keepCal;
  out->replCal      = in->replCal;

  /* Selected feed list */
  if ((in->feeds!=NULL) && (in->numberFeedList>0)) {
    if (out->feeds) g_free(out->feeds);
    out->feeds = g_malloc (in->numberFeedList*sizeof(olong));
    out->numberFeedList = in->numberFeedList;
    for (i=0; i<in->numberFeedList; i++) out->feeds[i] = in->feeds[i];
  }

  /* Selected target list */
  if ((in->targets!=NULL) && (in->numberTargetList>0)) {
    if (out->targets) g_free(out->targets);
    out->targets = g_malloc (in->numberTargetList*sizeof(olong));
    out->numberTargetList = in->numberTargetList;
    for (i=0; i<in->numberTargetList; i++) out->targets[i] = in->targets[i];
  }

  return out;
} /* end ObitOTFSelCopy */

/**
 * Determines how large a buffer (in floats) is needed
 * for data transfers as described by data members.
 * \param desc Pointer input descriptor.
 * \param sel OTF selector.
 * \return size in floats needed for I/O.
 */
olong ObitOTFSelBufferSize (ObitOTFDesc* desc, ObitOTFSel* sel)
{
  olong size = 0;

  /* error checks */
  if (desc==NULL) return size; 
  g_assert (ObitIsA(desc, ObitOTFDescGetClass()));
  g_assert (ObitIsA(sel, &myClassInfo));

  /* make sure defaults filled in */
  ObitOTFSelDefault (desc, sel);

  /* size of uncompressed vis * number of vis */
  size = desc->lrec * MAX (20, sel->nRecPIO);

  return size;
} /* end ObitOTFSelBufferSize */

/**
 * Not much to do
 * \param in Pointer to descriptor.
 * \param sel OTF selector, output vis descriptor changed if needed.
 */
void ObitOTFSelDefault (ObitOTFDesc* in, ObitOTFSel* sel)
{

  /* error checks */
  g_assert (ObitIsA(in, ObitOTFDescGetClass()));
  g_assert (ObitIsA(sel, &myClassInfo));

  /* Index descriptor as well */
  ObitOTFDescIndex(in);
} /* end ObitOTFSelDefault */

/**
 * Derive the descriptor for data being written; 
 * Doesn't really change.
 * \param in Pointer to input descriptor, this describes the data
 *           as they appear in memory.
 * \param sel OTF selector
 * \param out Pointer to output descriptor, describing form on disk.
 * \param err Obit error stack
 */
void ObitOTFSelGetDesc (ObitOTFDesc* in, ObitOTFSel* sel,
			  ObitOTFDesc* out, ObitErr *err)
{
  gchar *routine = "ObitOTFSelGetDesc";
  
  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, ObitOTFDescGetClass()));
  g_assert (ObitIsA(sel, &myClassInfo));
  g_assert (ObitIsA(out, ObitOTFDescGetClass()));

  /* make sure defaults filled in */
  ObitOTFSelDefault (in, sel);

  /* copy most values */
  ObitOTFDescCopy (in, out, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* make sure defaults, indices filled in */
  ObitOTFSelDefault (in, sel);
  ObitOTFSelDefault (out, sel);

} /* end ObitOTFSelGetDesc */

/**
 * Apply selection criteria to input descriptor to derive output.
 * Note: descriptor modification due to data selection is mostly done in
 * ObitOTFCalApply.
 * Also sets any previously undefined values on sel.
 * \param in  Pointer to input descriptor, this describes the data
 *            as they appear on disk (possibly compressed).
 * \param sel OTF selector, blc, trc members changed if needed.
 * \param out Pointer to output descriptor, this describes the data 
 *            after any processing when read.
 * \param err Obit error stack
 */
void ObitOTFSelSetDesc (ObitOTFDesc* in, ObitOTFSel* sel,
			  ObitOTFDesc* out, ObitErr *err)
{
  gchar *routine = "ObitOTFSelSetDesc";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, ObitOTFDescGetClass()));
  g_assert (ObitIsA(sel, &myClassInfo));
  g_assert (ObitIsA(out, ObitOTFDescGetClass()));

  /* make sure defaults filled in */
  ObitOTFSelDefault (in, sel);

  /* copy most values */
  ObitOTFDescCopy (in, out, err);
  if (err->error) /* add traceback, return on error */
      Obit_traceback_msg (err, routine, in->name);

  /* make sure defaults, indices filled in */
  ObitOTFSelDefault (in, sel);
  ObitOTFSelDefault (out, sel);

  /* set data increments as float. */
  sel->jincs    = in->incs;
  sel->jincf    = in->incf;
  sel->jincfeed = in->incfeed;

  /* Frequency Selection */
  if (sel->numberChann<=0) sel->numberChann = in->inaxes[in->jlocf];
  if (sel->startChann<=0) sel->startChann = 1;
  sel->numberChann = MIN (sel->numberChann, in->inaxes[in->jlocf]);

  /* Number of Feeds */
  sel->numberFeed = in->inaxes[in->jlocfeed];

} /* end ObitOTFSelSetDesc */

/**
 * See if an OTFIndex table exists and if so initialize it to use in deciding
 * which visibilities to read.
 * \param  in      Pointer to the object.
 * \param  desc    OTF descriptor from IO where the next visibility to
 *                 read and the number will be stored.
 * \param  err     Error stack
 * \return TRUE is finished, else FALSE
 */
void ObitOTFSelNextInit (ObitOTFSel *in, ObitOTFDesc *desc, ObitErr *err)
{
  ObitIOCode retCode;
  gchar *routine="ObitOTFSelNextInit";
 
  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitOTFSelIsA(in));
  g_assert (ObitOTFDescIsA(desc));

  /* Open Index table  */
  retCode = 
    ObitTableOTFIndexOpen ((ObitTableOTFIndex*)(in->IndexTable), OBIT_IO_ReadWrite,
			   err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_msg (err, routine, in->name);
  in->numRow = ((ObitTableOTFIndex*)in->IndexTable)->myDesc->nrow;
  in->LastRowRead = 0;

  /* Initialize */
  in->scanFirstRec = -1;
  in->scanLastRec  = -1;

  /* Don't bother if table empty */
  if (in->numRow<=0) {
    in->doIndex = FALSE;
    /* Close table  */
    retCode = ObitTableOTFIndexClose ((ObitTableOTFIndex*)in->IndexTable, err);
    if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
      Obit_traceback_msg (err, routine, in->name);
    return;
}
  
  /* Create row structure */
  in->IndexTableRow = 
    (Obit*)newObitTableOTFIndexRow((ObitTableOTFIndex*)(in->IndexTable));

  /* We want indexing */
  in->doIndex = TRUE; 
  return;
} /* end ObitOTFSelNextInit */

/**
 * Uses selector member to decide which records to read next.
 * If doIndex is TRUE, then use index table.
 * \param  in      Pointer to the object.
 * \param  desc    OTF descriptor from IO where the next visibility to
 *                 read and the number will be stored.
 * \param  err     Error stack
 * \return TRUE if finished, else FALSE
 */
gboolean ObitOTFSelNext (ObitOTFSel *in, ObitOTFDesc *desc, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitTableOTFIndexRow *row;
  olong nleft;
  gboolean gotIt, done = FALSE;
  gchar *routine = "ObitOTFSelNext";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return done;
  g_assert (ObitOTFSelIsA(in));
  g_assert (ObitOTFDescIsA(desc));
  
  /* Is this controlled by an index? */
  if (in->doIndex) {
    /* Which records are wanted? */
    if (desc->firstRec < 1) { /* first read? */
      desc->firstRec = 1;
    } else { /* subsequent reads */
      desc->firstRec += in->numRecRead;
    }

     /* Need a new scan? */
    if (desc->firstRec>in->scanLastRec) {
      
      /* Read index file until a selected scan found */
      row = (ObitTableOTFIndexRow*)in->IndexTableRow;
      gotIt = FALSE;
      while ((!gotIt) && (in->LastRowRead<in->numRow)) {
	in->LastRowRead++;
	retCode = ObitTableOTFIndexReadRow ((ObitTableOTFIndex*)in->IndexTable, 
				      in->LastRowRead, row, err);
	if (err->error) Obit_traceback_val (err, routine, in->name, retCode);
	
	/* Is this one wanted? */
	/* Any overlap with time range? */
	if (row->Time+row->TimeI < in->timeRange[0]) continue;
	if (row->Time > in->timeRange[1]) continue;

	/* A selected scan? */
	if ((row->ScanID<in->scans[0]) || (row->ScanID>in->scans[1])) continue;
	
	/* Requested Target? */
	gotIt = ObitOTFSelWantTarget(in, row->TargetID);
      } /* end loop over file */
      
	/* Save info if scan found */
      if (gotIt) {
	in->scanFirstRec = row->StartRec;
	in->scanLastRec  = row->EndRec;
	desc->firstRec   = in->scanFirstRec; /* beginning of scan */
      } else {
	/* No more scans - Must be finished */
	done = TRUE;
	return done;
      }
    } /* end get new scan */
    
    /* how many to attempt to read? */
    in->numRecRead = in->nRecPIO;
    
    /* but not more than all of scan */
    nleft = in->scanLastRec - desc->firstRec + 1;
    in->numRecRead  = MIN (nleft, in->numRecRead);
    in->numRecRead  = MAX (0, in->numRecRead);
    done = (nleft<=0);

    /* And no more than in table -in case Index bad */
    if (desc->firstRec+in->numRecRead > desc->nrecord) {
      nleft = desc->nrecord - desc->firstRec + 1;
      in->numRecRead = MIN (nleft, in->numRecRead);
      in->numRecRead = MAX (0, in->numRecRead);
      done = (nleft<=0);
    }
    
  } else {
    /* not indexed */
    /* Which records are wanted? */
    if (desc->firstRec < 1) { /* first read? */
      desc->firstRec = 1;
    } else { /* subsequent reads */
      desc->firstRec += in->nRecPIO;
    }
    
    /* how many? */
    in->numRecRead = in->nRecPIO;
    
    /* but not more than all */
    nleft = desc->nrecord - desc->firstRec + 1;
    in->numRecRead = MIN (nleft, in->numRecRead);
    in->numRecRead = MAX (0, in->numRecRead);
    done = (nleft<=0);
  }
  
  /* Don't try to read past EOF */
  done = done || (desc->firstRec>desc->nrecord);

  return done;
} /* end ObitOTFSelNext */

/**
 * Close Index table if open
 * \param  in   Pointer to the Selector.
 * \param  err  Error stack
 * \return TRUE is finished, else FALSE
 */
void ObitOTFSelShutdown (ObitOTFSel *in, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  gchar *routine = "ObitOTFSelShutdown";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitOTFSelIsA(in));

  /* Anything to do? */
  if (!in->doIndex) return;

 /* Close table  */
  retCode = ObitTableOTFIndexClose ((ObitTableOTFIndex*)in->IndexTable, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Release row structure  */
  in->IndexTableRow = ObitTableOTFIndexRowUnref(in->IndexTableRow);
} /* end ObitOTFSelShutdown */


/**
 * Determine if a given target is selected.
 * \param sel    OTF selector.
 * \param TargetID Target ID to be tested, <=0 => any
 * \return TRUE if target selected.
 */
gboolean ObitOTFSelWantTarget (ObitOTFSel* sel, olong TargetID)
{
  olong i;

  /* error checks */
  g_assert (ObitIsA(sel, &myClassInfo));

  /* TargetID<=0 matches any */
  if (TargetID<=0) return TRUE;

  /* If array is null - everything selected */
  if (sel->targets == NULL) return TRUE;

  /* ditto no entries */
  if (sel->numberTargetList == 0) return TRUE;

  /* ditto first entry <= 0 */
  if (sel->targets[0]<= 0) return TRUE;

  for (i=0; i<sel->numberTargetList; i++) {
    if (TargetID==sel->targets[i]) return TRUE;
  }
  /* didn't make the cut */
  return FALSE;

} /* end ObitOTFSelWantTarget */

/**
 * Determine if a given Feed is selected.
 * \param sel    OTF selector.
 * \param feed   feed id to test (0-rel)
 * \return TRUE if feed selected.
 */
gboolean ObitOTFSelWantFeed (ObitOTFSel* sel, olong feed)
{
  olong i, lfeed=feed+1;

  /* error checks */
  g_assert (ObitIsA(sel, &myClassInfo));

  /* If array is null - everything selected */
  if (sel->feeds == NULL) return TRUE;

  /* ditto no entries */
  if (sel->numberFeedList == 0) return TRUE;

  /* ditto first entry zero */
  if (sel->feeds[0] == 0) return TRUE;

  for (i=0; i<sel->numberFeedList; i++) {
    if (lfeed==sel->feeds[i]) return TRUE;
  }
  /* didn't make the cut */
  return FALSE;

} /* end ObitOTFSelWantFeed */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitOTFSelClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitOTFSelClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitOTFSelClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitOTFSelClassInfoDefFn (gpointer inClass)
{
  ObitOTFSelClassInfo *theClass = (ObitOTFSelClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitOTFSelClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitOTFSelClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitOTFSelGetClass;
  theClass->ObitClear     = (ObitClearFP)ObitOTFSelClear;
  theClass->ObitInit      = (ObitInitFP)ObitOTFSelInit;
  theClass->newObit       = (newObitFP)newObitOTFSel;
  theClass->ObitCopy      = (ObitCopyFP)ObitOTFSelCopy;
  theClass->ObitClone     = NULL;

} /* end ObitOTFSelClassDefFn */

/*---------------Private functions--------------------------*/
/**
 * Creates empty member objects, initialize reference count.
 * Does (recursive) initialization of base class members before 
 * this class.
 * \param inn Pointer to the object to initialize.
 */
void ObitOTFSelInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitOTFSel *in = inn;

  /* error checks */
  g_assert (in != NULL);
  
  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->nRecPIO       = 1;
  in->doCalSelect   = FALSE;
  in->feeds   = NULL;
  in->numberFeedList = 0;
  in->targets = NULL;
  in->numberTargetList = 0;
  in->IndexTable    = NULL;
  in->IndexTableRow = NULL;
  in->keepCal       = TRUE;
  in->replCal       = FALSE;
  
} /* end ObitOTFSelInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 */
void ObitOTFSelClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitOTFSel *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* free this class members */
  if (in->feeds)   g_free(in->feeds);   in->feeds = NULL;
  if (in->targets) g_free(in->targets); in->targets = NULL;
  in->IndexTable    = ObitTableOTFIndexUnref((ObitTableOTFIndex*)in->IndexTable);
  in->IndexTableRow = ObitTableOTFIndexRowUnref((ObitTableOTFIndexRow*)in->IndexTableRow);
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);

} /* end ObitOTFSelClear */


