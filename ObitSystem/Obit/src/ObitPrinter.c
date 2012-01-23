/* $Id$        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2009,2012                                          */
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
/*;Correspondence about this software should be addressed as follows: */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/

#include "ObitPrinter.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitPrinter.c
 * ObitPrinter class function definitions.
 * This class is derived from the Obit base class.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitPrinter";

/** Function to obtain parent ClassInfo */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo structure ObitPrinterClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitPrinterClassInfo myClassInfo = {FALSE};

/*--------------- File Global Variables  ----------------*/


/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitPrinterInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitPrinterClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitPrinterClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitPrinter* newObitPrinter (gchar* name)
{
  ObitPrinter* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitPrinterClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitPrinter));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitPrinterInit((gpointer)out);

 return out;
} /* end newObitPrinter */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitPrinterGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitPrinterClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitPrinterGetClass */

/**
 * Make a deep copy of an ObitPrinter.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitPrinter* ObitPrinterCopy  (ObitPrinter *in, ObitPrinter *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  gboolean oldExist;
  gchar *outName;

  /* error checks */
  if (err->error) return out;
  g_assert (ObitIsA(in, &myClassInfo));
  if (out) g_assert (ObitIsA(out, &myClassInfo));

  /* Create if it doesn't exist */
  oldExist = out!=NULL;
  if (!oldExist) {
    /* derive object name */
    outName = g_strconcat ("Copy: ",in->name,NULL);
    out = newObitPrinter(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  out->isInteractive = in->isInteractive;
  out->outFile       = in->outFile;
  out->LinesPerPage  = in->LinesPerPage;
  if (out->outFileName) g_free(out->outFileName); out->outFileName=NULL;
  if (in->outFileName) out->outFileName = g_strdup(in->outFileName);
  if (out->Title1) g_free(out->Title1); out->Title1 = NULL;
  if (in->Title1) out->Title1 = g_strdup(in->Title1);
  if (out->Title2) g_free(out->Title2); out->Title2 = NULL;
  if (in->Title2) out->Title2 = g_strdup(in->Title2);

  return out;
} /* end ObitPrinterCopy */

/**
 * Make a copy of a object but do not copy the actual data
 * This is useful to create an Printer similar to the input one.
 * \param in  The object to copy
 * \param out An existing object pointer for output, must be defined.
 * \param err Obit error stack object.
 */
void ObitPrinterClone  (ObitPrinter *in, ObitPrinter *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitIsA(out, &myClassInfo));

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  out->isInteractive = in->isInteractive;
  out->outFile       = in->outFile;
  out->LinesPerPage  = in->LinesPerPage;
  if (out->outFileName) g_free(out->outFileName); out->outFileName=NULL;
  if (in->outFileName) out->outFileName = g_strdup(in->outFileName);
  if (out->Title1) g_free(out->Title1); out->Title1 = NULL;
  if (in->Title1) out->Title1 = g_strdup(in->Title1);
  if (out->Title2) g_free(out->Title2); out->Title2 = NULL;
  if (in->Title2) out->Title2 = g_strdup(in->Title2);

} /* end ObitPrinterClone */

/**
 * Creates an ObitPrinter 
 * \param name  An optional name for the object.
 * \param isInteractive If true, output is interactive and written
 *                      to outStream, answers to questions read from
                        stdin.
 * \param outStream     Output non file device, e.g. stdout  
 * \param fileName      Name of output ascii file 
 *                      ignored for interactive
 *                      If NULL, all output will be to outStream
 * \return the new object.
 */
ObitPrinter* ObitPrinterCreate (gchar* name, gboolean isInteractive,
				FILE *outStream, gchar *fileName)
{
  ObitPrinter* out;

  /* Create basic structure */
  out = newObitPrinter (name);

  /* Save entries */
  out->isInteractive = isInteractive;
  out->outStream     = outStream;
  if (fileName) out->outFileName = g_strdup(fileName);
  out->outFile   = newObitFile(name);
  out->myStatus  = OBIT_Inactive;

  return out;
} /* end ObitPrinterCreate */

/**
 * Open an ObitPrinter 
 * \param in            The ObitPrinter object to Open
 *                      info list contains:
 *   "pageLimit"  OBIT_long  Maximum number of pages of output def [50]
 * \param LinesPerPage  Number of lines per page, if <=0:
 *                      24 if interactive, 80 otherwise
 * \param Title1        First title string for each page
 * \param Title2        Second title string for each page
 * \param err           Obit error stack object.
 */
void ObitPrinterOpen (ObitPrinter *in, olong LinesPerPage, gchar *Title1, 
		      gchar *Title2, ObitErr *err)
{
  ObitInfoType type;
  union ObitInfoListEquiv InfoReal; 
  gint32       dim[MAXINFOELEMDIM];

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Page Limit from info */
  if (ObitInfoListGetTest(in->info, "pageLimit", &type, dim, &InfoReal)) {
    if (type==OBIT_oint) in->pageLimit = InfoReal.itg;
    else                 in->pageLimit = InfoReal.otg;
    in->pageLimit = MAX (1, in->pageLimit);
  }

  /* Save values - lines per page */
  if (LinesPerPage<=0) {
    if (in->isInteractive) in->LinesPerPage = 24;
    else                   in->LinesPerPage = 80;
  } else in->LinesPerPage = LinesPerPage;
  in->Title1 = g_strdup(Title1);
  ObitTrimTrail(in->Title1);
  in->Title2 = g_strdup(Title2);
  ObitTrimTrail(in->Title2);

  /* Open file if noninteractive */
  if (!in->isInteractive && in->outFileName) {
    ObitFileOpen(in->outFile, in->outFileName, OBIT_IO_ReadWrite, 
		 OBIT_IO_Text, 0, err);
    if (err->error) 
      Obit_traceback_msg (err, "Error in ObitFileOpen", in->name);
 }

  in->myStatus = OBIT_Active;  /* seems to be OK */
  in->nPages = 0;
  in->nLines = 0;
} /* end ObitPrinterOpen */

/**
 * Add a line to a print output
 * Add page headers if needed and asks user of interactive input
 * if he/she wishes to quit.
 * \param in            The object to write
 * \param line          Line of text to add
 * \param quit          [out] If TRUE, user has had enough
 * \param err            Obit error stack object.
 */
void ObitPrinterWrite (ObitPrinter *in, gchar *line, gboolean *quit, 
		       ObitErr *err)
{
  gchar lline[2000];
  gchar *routine = "ObitPrinterWrite";

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (line != NULL);
  /* A previous error condition? */
  if (in->myStatus==OBIT_ErrorExist) return;
  Obit_return_if_fail((in->myStatus==OBIT_Active), err,
 		      "%s: Printer object %s not opened", routine, in->name);

  *quit = FALSE;  /* Until determined otherwise */

  ObitTrimTrail (line);  /* no trailing blanks */
  /* Write file if noninteractive */
  if (!in->isInteractive && in->outFileName) {
    sprintf (lline, "%s\n", line);
    ObitFileWriteLine(in->outFile, lline, err);
    if (err->error) {
      in->myStatus = OBIT_ErrorExist;
      Obit_traceback_msg (err, "Error in ObitFileWriteLine", in->name);
    }
    /* Line statistics */
    in->nLines++;
    /* Time for a new page? */
    if (in->nLines>=in->LinesPerPage) ObitPrinterNewPage (in, quit, err);
  } else { /* interactive */
    fprintf (in->outStream, "%s\n", line);
    /* Line statistics */
    in->nLines++;
    /* Time for a new page? */
    if (in->nLines>=in->LinesPerPage) ObitPrinterNewPage (in, quit, err);
  } /* end interactive */
  /* Error? */
  if (err->error) {
    in->myStatus = OBIT_ErrorExist;
    Obit_traceback_msg (err, "Error in ObitFileWriteLine", in->name);
  }
} /* end ObitPrinterWrite */

/**
 * Close a printer
 * \param in            The object to write
 * \param err            Obit error stack object.
 */
void ObitPrinterClose (ObitPrinter *in, ObitErr *err)
{
  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  if (in->myStatus==OBIT_ErrorExist) return;
  if (in->myStatus!=OBIT_Active) return;

  /* Close file if noninteractive */
  if (!in->isInteractive && in->outFile) {
    ObitFileFlush (in->outFile, err);
    ObitFileClose (in->outFile, err);
    if (err->error) {
      in->myStatus = OBIT_ErrorExist;
      Obit_traceback_msg (err, "Error in ObitFileClose", in->name);
    }
  } else { /* Flush interactive */
    fflush (in->outStream);
    fclose (in->outStream);
  }

  in->myStatus  = OBIT_Inactive;
} /* end ObitPrinterClose */

/**
 * Update titles
 * \param in       The object to update
 * \param Title1   First title string for each page
 * \param Title2   Second title string for each page
 * \param err      Obit error stack object.
 */
void ObitPrinterSetTitle (ObitPrinter *in, gchar *Title1, 
			  gchar *Title2, ObitErr *err)
{
  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  if (in->Title1) g_free(in->Title1); in->Title1 = NULL;
  if (Title1) in->Title1 = g_strdup(Title1);
  if (in->Title2) g_free(in->Title2);in->Title2 = NULL;
  if (Title2) in->Title2  = g_strdup(Title2);
} /* end ObitPrinterSetTitle */

/**
 * Start a new page on the output
 * Titles will be added if defined (non-NULL)
 * For interactive use, the used will be asked whether to continue or quit.
 * \param in       The object to update
 * \param quit     [out] If TRUE, user has had enough
 * \param err      Obit error stack object.
 */
void ObitPrinterNewPage (ObitPrinter *in, gboolean *quit, ObitErr *err)
{
  gchar lline[138];
  gchar *routine = "ObitPrinterNewPage";

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  /* A previous error condition? */
  if (in->myStatus==OBIT_ErrorExist) return;
  Obit_return_if_fail((in->myStatus==OBIT_Active), err,
 		      "%s: Printer object %s not open", routine, in->name);

  *quit = FALSE;  /* Until determined otherwise */

  /* New page */
  in->nLines = 0;
  in->nPages++;
  if (!in->isInteractive && in->outFileName) {  /* Non interactive */
    sprintf (lline,"                                Page %d\n",in->nPages);
    in->nLines++;
    ObitFileWriteLine(in->outFile, lline, err);
    if (in->Title1) {
      sprintf (lline, "%s\n", in->Title1);
      ObitFileWriteLine(in->outFile, lline, err);
      in->nLines++;
    }
    if (in->Title2) {
      sprintf (lline, "%s\n", in->Title2);
      ObitFileWriteLine(in->outFile, lline, err);
      in->nLines++;
    }
    if (err->error) {
      in->myStatus = OBIT_ErrorExist;
      Obit_traceback_msg (err, "Error in ObitFileWriteLine", in->name);
    }
  } else if (in->isInteractive) { /* interactive */
    /* Does the user want to quit */
    fprintf (in->outStream, "Type Q to stop, just hit RETURN to continue\n");
    /* Read response from stdin */
    fgets (lline, 100, stdin);
    if ((lline[0]=='q') || (lline[0]=='Q')) {
      /* Bail */
      *quit = TRUE;
      in->myStatus  = OBIT_Inactive;
      return;
    }
    if (in->Title1) {
      fprintf (in->outStream, "%s\n", in->Title1);
      in->nLines++;
    }
    if (in->Title2) {
      fprintf (in->outStream, "%s\n", in->Title2);
      in->nLines++;
    }
    /* end of interactive */
  } else {  /* non interactive stream */
    fprintf (in->outStream,"                                                Page %d\n",in->nPages);
    if (in->Title1) {
      fprintf (in->outStream, "%s\n", in->Title1);
      in->nLines++;
    }
    if (in->Title2) {
      fprintf (in->outStream, "%s\n", in->Title2);
      in->nLines++;
    }
  }

  /* Check Page Limit */
  if (in->nPages>=in->pageLimit) {
    if (!in->isInteractive && in->outFileName) {  /* Non interactive */
      sprintf (lline,"Hit page limit %d\n",in->pageLimit);
      in->nLines++;
      ObitFileWriteLine(in->outFile, lline, err);
    } else {  /* interactive or non interactive stream */
      fprintf (in->outStream,"Hit page limit %d\n",in->pageLimit);
    }
    /* Force Quit */
    *quit = TRUE;
    in->myStatus  = OBIT_Inactive;
  }
  
} /* end ObitPrinterNewPage */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitPrinterClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitPrinterClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitPrinterClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitPrinterClassInfoDefFn (gpointer inClass)
{
  ObitPrinterClassInfo *theClass = (ObitPrinterClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitPrinterClassInit;
  theClass->newObit       = (newObitFP)newObitPrinter;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitPrinterClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitPrinterGetClass;
  theClass->ObitCopy      = (ObitCopyFP)ObitPrinterCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitPrinterClear;
  theClass->ObitInit      = (ObitInitFP)ObitPrinterInit;
  theClass->ObitPrinterCreate = (ObitPrinterCreateFP)ObitPrinterCreate;
  theClass->ObitPrinterOpen   = (ObitPrinterOpenFP)ObitPrinterOpen;
  theClass->ObitPrinterWrite  = (ObitPrinterWriteFP)ObitPrinterWrite;
  theClass->ObitPrinterClose  = (ObitPrinterCloseFP)ObitPrinterClose;
  theClass->ObitPrinterSetTitle = 
    (ObitPrinterSetTitleFP)ObitPrinterSetTitle;
  theClass->ObitPrinterNewPage = 
    (ObitPrinterNewPageFP)ObitPrinterNewPage;

} /* end ObitPrinterClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitPrinterInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitPrinter *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->thread      = newObitThread();
  in->info        = newObitInfoList(); 
  in->outFile     = NULL;
  in->outStream   = NULL;
  in->outFileName = NULL;
  in->Title1      = NULL;
  in->Title2      = NULL;
  in->pageLimit   = 50;

} /* end ObitPrinterInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitPrinter* cast to an Obit*.
 */
void ObitPrinterClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitPrinter *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  in->thread    = ObitThreadUnref(in->thread);
  in->info      = ObitInfoListUnref(in->info);
  in->outFile   = ObitFileUnref(in->outFile);
  if (in->outFileName) g_free(in->outFileName);
  if (in->Title1)      g_free(in->Title1);
  if (in->Title2)      g_free(in->Title2);
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitPrinterClear */

