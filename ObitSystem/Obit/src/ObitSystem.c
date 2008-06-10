/* $Id$      */
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
#include <errno.h>
#include "Obit.h"
#include "ObitSystem.h"
#include "ObitAIPS.h"
#include "ObitFITS.h"
#include "ObitImage.h"
#include "ObitTable.h"
#include "ObitUV.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitSystem.c
 * ObitSystem Class definition file.
 */

/*--------------- File Global Variables  ----------------*/
/** name of the class defined in this file */
static gchar *myClassName = "ObitSystem";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo global structure ObitIOClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitSystemClassInfo myClassInfo = {FALSE};

/** Pointer to Class information structure */
static ObitSystem *mySystemInfo = NULL;

/*--------------- File Structure definitions  ----------------*/
/** An element of a scratchList */
typedef struct { 
  /** Obit pointer */
  Obit *item;
}  scratchListElem; 

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitSystemInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitSystemClear (gpointer in);

/** Private: Create a scratchListElem. */
static scratchListElem*
newscratchListElem (Obit *item);

/** Private: Delete a scratchListElem. */
static void freescratchListElem (scratchListElem *in);

/** Private: Add an scratchListElem to the list. */
static void scratchListAdd (ObitSystem *in, scratchListElem *elem);

/** Private: Remove an scratchListElem from the list. */
static void scratchListRemove (ObitSystem *in, scratchListElem *elem);

/** Private: Find item in a list */
static scratchListElem*  scratchListFind(ObitSystem *in, Obit *item);

/** Private: Set Class function pointers. */
static void ObitSystemClassInfoDefFn (gpointer inClass);

/*---------------Public functions---------------------------*/
/**
 * Initialize ObitSystem information
 * \li Initialize AIPS disk information if any
 * \li Initialize FITS disk information if any
 * \li Initialize Scratch file list.
 * \param pgmName        Name of program (max 5 char if AIPS)
 *                       If NULL or empty (zero length), no start/shutdown 
 *                       messages
 * \param pgmNumber      Version number of program (e.g. POPS number).
 * \param AIPSuser       AIPS user number if using AIPS files
 * \param numberAIPSdisk Number of AIPS disks
 *                       If 0, no AIPS files.
 * \param AIPSdir        List of AIPS directory names
 *                       If NULL none.
 * \param numberFITSdisk Number of FITS disks
 *                       If 0, no FITS files.
 * \param FITSdir        List of FITS directory names
 *                       If NULL none.
 * \param F_TRUE         Value of Fortran TRUE (used in Fortran interface)
 * \param F_FALSE        Value of Fortran FALSE
 * \param err            Obit error stack for any error messages.
 * \return pointer to object created.
 */
ObitSystem* 
ObitSystemStartup (gchar *pgmName, olong pgmNumber, 
		   olong AIPSuser,
		   olong numberAIPSdisk, gchar* AIPSdir[], 
		   olong numberFITSdisk, gchar* FITSdir[], 
		   oint F_TRUE, oint F_FALSE, ObitErr *err)
{
  ObitSystem* out;

  /* Init system error flag */
  errno = 0;

  /* initialize AIPS stuff */
  ObitAIPSClassInit (numberAIPSdisk, AIPSdir, F_TRUE, F_FALSE); 

  /* initialize FITS disk stuff */
  ObitFITSClassInit (numberFITSdisk, FITSdir); 

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitSystemClassInit();

  /* allocate structure */
  out = g_malloc0(sizeof(ObitSystem));

 /* set classInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitSystemInit((gpointer)out);

  /* initialize values */
  out->name      = g_strdup("ObitSystem");
  if (pgmName) out->pgmName = g_strdup(pgmName);
  else out->pgmName = g_strdup("Nameless");
  out->pgmNumber = pgmNumber;
  out->AIPSuser  = AIPSuser;
  out->numberAIPSdisk = numberAIPSdisk;
  out->numberFITSdisk = numberFITSdisk;

  /* make global pointer */
  mySystemInfo = (ObitSystem*)ObitRef(out);

  /* Save error/message stack object */
  out->err = ObitErrRef (err);

  /* Startup message if program name given */
  if ((strlen(out->pgmName)>0) && strncmp (out->pgmName, "NameLess", 8))
    Obit_log_error(out->err, OBIT_InfoErr, "%s Begins", out->pgmName);
  ObitErrTimeStamp(out->err);  /* Add Timestamp */
  ObitErrLog(out->err);

  return out;
} /*  end ObitSystemStartup */

/**
 * Returns ClassInfo pointer for the class.
 * Initializes class if needed on first call.
 * \return pointer to the class structure.
 */
gconstpointer ObitSystemGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) 
    ObitSystemClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitSystemGetClass */

/**
 * Shutdown Obit system
 * \li delete any remaining scratch objects
 * \param in   Object to shutdown
 * \return NULL pointer for object destroyed..
 */
ObitSystem* ObitSystemShutdown (ObitSystem* in)
{
   scratchListElem *elem=NULL;
   ObitClassInfo *myClass;
   Obit *tst;
   GSList *tmp;

   /* ignore if I haven't been started */
   if (!mySystemInfo) return NULL;

  /* Lock object aginst other threads */
  ObitThreadLock(mySystemInfo->thread);

  /* loop through scratchList  Unrefing elements */
  tmp = in->scratchList;
  while (tmp!=NULL) {
    elem = (scratchListElem*)tmp->data;
    /* Make sure still valid */
    if (elem) {
      tst = (Obit*)elem->item;
      myClass = (ObitClassInfo*)tst->ClassInfo;
      if (myClass  && myClass->hasScratch) { /* Scratch forms allowed? */
	while (tst && (tst->ReferenceCount>0)){ 
	  /* Zap it - this may modify in->scratchList */
	  tst = ObitUnref(tst);
	}
      }
    }

    /* List not changed, remove this entry anyway */
    if (tmp==in->scratchList)  scratchListRemove (in, elem);

    tmp = in->scratchList; /* go to the new(?) head of the list */
  } /* end loop over scratch List */

  /* Shutdown AIPS */
  ObitAIPSShutdown();

  /* Shutdown FITS */
  ObitFITSShutdown();

  ObitThreadUnlock(mySystemInfo->thread);

  /* Shutdown message if program name given */
  if ((strlen(in->pgmName)>0) && strncmp (in->pgmName, "NameLess", 8))
    Obit_log_error(in->err, OBIT_InfoErr, "%s Ends", in->pgmName);
  ObitErrTimeStamp(in->err);  /* Add Timestamp */
  ObitErrLog(in->err);

  /* delete object */
  mySystemInfo = ObitUnref(mySystemInfo );
  return ObitUnref(in);
} /*  end ObitSystemShutdown */

/**
 * Get disk name, etc information for a new scratch file.
 * \param FileType File system type (OBIT_IO_FITS, OBIT_IO_AIPS);
 * \param type     File type ("MA", "UV" for AIPS, not used for FITS)
 * \param info     ObitInfoList to write assignments to
 * \param err      Obit error stack for any error messages.
 */
void ObitSystemGetScratch (ObitIOType FileType,gchar *type,
			   ObitInfoList *info, ObitErr *err)
{
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong i;
  gchar *routine = "ObitSystemGetScratch";

 /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return ;  /* previous error? */
  g_assert (ObitInfoListIsA(info));
  Obit_return_if_fail((mySystemInfo!=NULL), err,
		      "%s: Obit not initialized",routine);

  /* Lock object aginst other threads */
  ObitThreadLock(mySystemInfo->thread);

  /* update scratch file and last disk numbers */
  mySystemInfo->numberScratch++;
  mySystemInfo->lastDisk++; /* spread scratch files out */

  /* save file type on info */
  ObitInfoListPut (info, "FileType", OBIT_long, dim, (gpointer)&FileType,   err);
  if (err->error) Obit_traceback_msg (err, routine, "System");

  if (FileType == OBIT_IO_FITS) {          /* FITS file */
    /* Keep disk number in bounds */
    if (mySystemInfo->lastDisk>mySystemInfo->numberFITSdisk)
      mySystemInfo->lastDisk = 1;

    /* Assign */
    ObitFITSAssign (mySystemInfo->pgmName, mySystemInfo->pgmNumber,
		    mySystemInfo->lastDisk, mySystemInfo->numberScratch, 
		    info, err);
    if (err->error) Obit_traceback_msg (err, routine, "System");

  } else if (FileType == OBIT_IO_AIPS) {  /* AIPS file */
    /* Check noScrat status for disk */
    for (i=0; i<mySystemInfo->numberAIPSdisk; i++) {
      /* Keep disk number in bounds */
      if (mySystemInfo->lastDisk>mySystemInfo->numberAIPSdisk)
	mySystemInfo->lastDisk = 1;
      if (!ObitAIPSisNoScrat(mySystemInfo->lastDisk)) break;
      /* Keep going */
      mySystemInfo->lastDisk++; /* spread scratch files out */
    } /* End of loop checking */

    /* Make sure some allowed */
    Obit_return_if_fail ((!ObitAIPSisNoScrat(mySystemInfo->lastDisk)), err,
			 "%s: No scratch files allowed on any AIPS disk",
			 routine);

    /* Assign */
    ObitAIPSAssign (mySystemInfo->pgmName, mySystemInfo->pgmNumber,
		    type, mySystemInfo->AIPSuser,
		    mySystemInfo->lastDisk, mySystemInfo->numberScratch, 
 		    info, err);
    if (err->error) Obit_traceback_msg (err, routine, "System");

  } else { /* should never get here */
    g_assert_not_reached(); /* barf and die */
  }

  /* Unlock object */
  ObitThreadUnlock(mySystemInfo->thread);
} /*  end ObitSystemGetScratch */

/**
 * Add a scratch object to the list.
 * \param in   Object (ObitUV or ObitImage) to add to list
 * \param err  Obit error stack for any error messages.
 */
void ObitSystemAddScratch (Obit *in, ObitErr *err)
{
   scratchListElem *elem=NULL;
   ObitClassInfo *myClass;
   gchar *routine="ObitSystemAddScratch";

 /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return ;  /* previous error? */
  Obit_return_if_fail((mySystemInfo!=NULL), err,
		      "%s: Obit not initialized",routine);

  /* Are scratch version allowed */
  myClass = (ObitClassInfo*)in->ClassInfo;
  g_assert (myClass  && myClass->hasScratch);

  /* Lock object aginst other threads */
  ObitThreadLock(mySystemInfo->thread);

  /* Add to list */
  elem = newscratchListElem (in);
  scratchListAdd (mySystemInfo, elem);

  /* Unlock object */
  ObitThreadUnlock(mySystemInfo->thread);
} /*  end ObitSystemAddScratch */

/**
 * Free a scratch object from the list.
 * \param in   Object (ObitUV or ObitImage) to remove from list
 * \param err  Obit error stack for any error messages.
 */
void ObitSystemFreeScratch (Obit *in, ObitErr *err)
{
   scratchListElem *elem=NULL;
   ObitClassInfo *myClass;
   gchar *routine="ObitSystemFreeScratch";

 /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return ;  /* previous error? */
  Obit_return_if_fail((mySystemInfo!=NULL), err,
		      "%s: Obit not initialized",routine);

  /* Are scratch version allowed */
  myClass = (ObitClassInfo*)in->ClassInfo;
  g_assert (myClass  && myClass->hasScratch);

  /* Lock object aginst other threads */
  ObitThreadLock(mySystemInfo->thread);

  /* remove from list */
  elem = scratchListFind (mySystemInfo, in);
  if (elem) scratchListRemove (mySystemInfo, elem);
  
  /* Unlock object */
  ObitThreadUnlock(mySystemInfo->thread);
} /*  end ObitSystemFreeScratch */

/**
 * Tell if Obit is initialized
 * \return TRUE if Obit initialized, else FALSE
 */
gboolean ObitSystemIsInit (void)
{
  gboolean isInit = mySystemInfo!=NULL;
  return isInit;
} /* end ObitSystemIsInit */

/**
 * Tell program name
 * \return Program name
 */
gchar* ObitSystemGetPgmName (void)
{
  if (mySystemInfo==NULL) {
    g_warning ("Obit not initialized");
    return NULL;
  }
  return mySystemInfo->pgmName;
} /* end ObitSystemGetPgmName */

/**
 * Set program name
 * \param pgmName Program name
 */
void ObitSystemSetPgmName (gchar *pgmName)
{
  if (mySystemInfo==NULL) {
    g_warning ("Obit not initialized");
    return;
  }

  if (mySystemInfo->pgmName) g_free(mySystemInfo->pgmName);
  mySystemInfo->pgmName = g_strdup(pgmName);
} /* end ObitSystemSetPgmName */

/**
 * Tell program number
 * \return Program number
 */
olong ObitSystemGetPgmNumber (void)
{
  if (mySystemInfo==NULL) {
    g_warning ("Obit not initialized");
    return 0;
  }
  return mySystemInfo->pgmNumber;
} /* end ObitSystemGetPgmNumber */

/**
 * Set program number
 * \param pgmNumber Program number
 */
void ObitSystemSetPgmNumber (gint pgmNumber)
{
  if (mySystemInfo==NULL) {
    g_warning ("Obit not initialized");
    return;
  }

  mySystemInfo->pgmNumber = pgmNumber;
} /* end ObitSystemSetPgmName */

/**
 * Tell AIPS user number
 * \return AIPS user number
 */
olong ObitSystemGetAIPSuser (void)
{
  if (mySystemInfo==NULL) {
    g_warning ("Obit not initialized");
    return 0;
  }
  return mySystemInfo->AIPSuser;
} /* end ObitSystemGetAIPSuser */

/**
 * Set AIPS User ID member
 * \param AIPSuser       AIPS user number
 */
void ObitSystemSetAIPSuser (gint AIPSuser)
{
  if (mySystemInfo==NULL) {
    g_warning ("Obit not initialized");
    return;
  }
   mySystemInfo->AIPSuser  = AIPSuser;
} /* end ObitSystemSetAIPSuser */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitSystemClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitSystemClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitSystemClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitSystemClassInfoDefFn (gpointer inClass)
{
  ObitSystemClassInfo *theClass = (ObitSystemClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitSystemClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitSystemClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitSystemGetClass;
  theClass->newObit       = NULL;
  theClass->ObitCopy      = NULL;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitSystemClear;
  theClass->ObitInit      = (ObitInitFP)ObitSystemInit;
  theClass->ObitSystemStartup      = 
    (ObitSystemStartupFP)ObitSystemStartup;
  theClass->ObitSystemShutdown     = 
    (ObitSystemShutdownFP)ObitSystemShutdown;
  theClass->ObitSystemGetScratch   = 
    (ObitSystemGetScratchFP)ObitSystemGetScratch;
  theClass->ObitSystemAddScratch   = 
    (ObitSystemAddScratchFP)ObitSystemAddScratch;
  theClass->ObitSystemFreeScratch  = 
    (ObitSystemFreeScratchFP)ObitSystemFreeScratch;

} /* end ObitSystemClassDefFn */

/*---------------Private functions--------------------------*/
/**
 * Creates empty member objects, initialize reference count.
 * Does (recursive) initialization of base class members before 
 * this class.
 * \param inn Pointer to the object to initialize.
 */
void ObitSystemInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitSystem *in = inn;

  /* error checks */
  g_assert (in != NULL);
  
  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->thread        = newObitThread();
  in->pgmName       = NULL;
  in->pgmNumber     = -1;
  in->numberScratch = 0;
  in->lastDisk      = 0;
  in->number        = 0;
  in->scratchList   = NULL;
  in->err           = NULL;

} /* end ObitSystemInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 */
void ObitSystemClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitSystem *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* free this class members */
   in->thread    = ObitThreadUnref(in->thread);
   if (in->pgmName) g_free(in->pgmName);
   g_slist_free(in->scratchList);
   ObitErrUnref(in->err);

  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);

} /* end ObitSystemClear */


/**
 * scratchListElem Constructor
 * \param item  Data
 * \return the new  object.
 */
static scratchListElem* newscratchListElem (Obit *item)
{
  scratchListElem *out=NULL;

  out = g_malloc0(sizeof(scratchListElem));
  out->item = item;
  return out;
} /* end newscratchListElem */

/**
 * Destructor 
 * \param in Object to delete
 */
static void freescratchListElem (scratchListElem *in)
{
  if (in) g_free(in);
} /* end freescratchListElem */

/**
 * Attach elem to list in
 * \param in   list to add elem to
 * \param elem the element to add.
 */
static void scratchListAdd (ObitSystem *in, scratchListElem *elem)
{
  scratchListElem *tmp;

  /* is it already there? */
  tmp = scratchListFind(in, elem->item);
  if (tmp!=NULL) {/* Yes - don't add again */
    freescratchListElem(elem); /* delete object */
    return; 
  }

  /* link to list */
  in->scratchList = g_slist_append (in->scratchList, elem);
  in->number++;

} /* end scratchListAdd */

/**
 * Remove elem from list in
 * \param in   list to remove elem from
 * \param elem the element to remove.
 */
static void scratchListRemove (ObitSystem *in, scratchListElem *elem)
{
 scratchListElem *tmp;

 if (elem==NULL) return;  /* anything to do? */
  /* is it there? */
  tmp = scratchListFind(in, elem->item);
  if (tmp==NULL) {/* No - return */
    return; 
  }
 
 /* remove from list */
  in->scratchList = g_slist_remove(in->scratchList, elem);
  in->number--; /* keep count */

} /* end scratchListRemove  */

/**
 * Find an item in list in
 * \param in   list to search
 * \param item to search for
 * \return pointer to element containing item, NULL if not found.
 */
static scratchListElem* scratchListFind (ObitSystem *in, Obit *item)
{
  GSList *tmp;
  scratchListElem *out = NULL;

  /* loop through list testing elements */
  tmp = in->scratchList;
  while (tmp!=NULL) {
    out = (scratchListElem*)tmp->data;
    /* check if this is a match, to the pointers point to the same thing? */
    if (item==out->item) return out;
    tmp = g_slist_next(tmp);
  }

  return NULL; /* didn't find */
} /* end scratchListFind */

