/* $Id$            */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2002-2008                                          */
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
/*; Correspondence about this software should be addressed as follows:*/
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/

#include <time.h>
#include "Obit.h"
#include "ObitMem.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 *  \file Obit.c
 * Obit class function definitions.
 */

/* Documentation for doxygen main page */
/**
 *  \mainpage Obit Classes
 * Obit uses a class derivation scheme that doxygen does not understand
 * so some care is needed in interpreting this documentation.
 * Class hierarchies are generally noted in the names of modules, i.e.
 * Obit is the base class from which (almost) all others are derived.
 * Obit class derivation is by means of nested include files; each class has an
 * include file for the data members and for the class function pointers.
 * These include files include the corresponding includes of their parent class.
 * 
 * The most effective use of this documentation page is to use the "File
 * List" function and select the header file (.h file) for the given
 * class, this gives the functional interface to the class.
 * Class data members can be viewed using the links to the class
 * structure on this page or use the "Class List" function on the main
 * page and select the desired class.
 * 
 * Each class has a "ClassInfo" structure with class specific
 * information, mostly, function pointers.
 * Each instance of a class has a number of data members, including a
 * pointer to the ClassInfo structure of its class.
 */

/*--------------- file global data ----------------*/
/** name of the class defined in this file */
static gchar *myClassName = "Obit";

/**
 * ClassInfo structure ObitClassInfo.
 * This structure is used by class objects to access class functions.
 */
ObitClassInfo myClassInfo = {FALSE};

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
static void  ObitInit  (gpointer in);

/** Private: Deallocate members. */
static void  ObitClear (gpointer in);

/*----------------------Public functions---------------------------*/
/**
 *  Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
Obit* newObit (gchar* name)
{
  Obit* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(Obit));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitInit((gpointer)out);

  return out;
} /* end newObit */

/**
 * Returns ClassInfo pointer for the class.
 * This method MUST be included in each derived class to ensure
 * proper linking and class initialization.
 * Initializes class if needed on first call.
 * \return pointer to the class structure.
 */
gconstpointer ObitGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitGetClass */

/**
 * Make a deep copy of input object.
 * Copies are made of complex members such as files; these will be copied
 * applying whatever selection is associated with the input.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Error stack, returns if not empty.
 * \return pointer to the new (existing) object.
 */
Obit* ObitCopy (Obit *in, Obit *out, ObitErr *err)
{
  gboolean oldExist;

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return NULL;
  g_assert (ObitIsA(in, &myClassInfo));
  if (out) g_assert (ObitIsA(out, &myClassInfo));

  /* Create if it doesn't exist */
   oldExist = out!=NULL;
   if (!out) out = newObit(NULL);

  /* deep copy this class members if newly created */
   if (!oldExist) {
   }

  return out;
} /* end ObitCopy */

/**
 * Make a shallow copy of an object.
 * The result will have pointers to the more complex members.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \return pointer to the new object.
 */
Obit* ObitClone  (Obit *in, Obit *out)
{

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  if (out) g_assert (ObitIsA(out, &myClassInfo));

  /* Create if it doesn't exist */
  if (!out) out = newObit(NULL);

  return out;
} /* end ObitClone */

/**
 * To reference a pointer, incrementing reference count and returning 
 * the pointer.
 * This function should always be used to copy pointers as this 
 * will ensure a proper reference count.
 * Should also work for derived classes
 * \param in Pointer to object to link, if Null, just return.
 * \return the pointer to in.
 */
gpointer ObitRef (gpointer in)
{
  gpointer out = in;
  if (in==NULL) return in;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* increment reference count */
  ((Obit*)in)->ReferenceCount++;

  return out;
} /* end ObitRef */

/**
 * Always use this function to dismiss an object as it will
 * ensure that the object is only deleted when there are no more 
 * pointers to it.
 * if the input pointer is NULL, the reference count is already <=0 
 * or the object is not a valid Obit Object, it simply returns.
 * \param  in Pointer to object to unreference.
 * \return NULL pointer.
 */
gpointer ObitUnref (gpointer inn)
{
  ObitClassInfo *MyClass;
  Obit *in = (Obit*)inn;

  /* Ignore NULL pointers */
  if (in==NULL) return NULL;

  /* Do nothing if not valid Obit Object */
  if (in->ObitId != OBIT_ID) return NULL;

  /* Do nothing if reference count already <=0 */
  if (in->ReferenceCount<=0)  return NULL;

  /* error check */
  g_assert (ObitIsA(in, &myClassInfo));

  /* decrement ref count */
  in->ReferenceCount--;

  /* Deallocate if it goes negative */
  if (in->ReferenceCount<1) {

    /* Free members */
    /* Who is the guy? */
    MyClass = in->ClassInfo;
    if ((MyClass!=NULL) && (MyClass->ObitClear!=NULL)) 
      MyClass->ObitClear (in);

    /* free structure - may be ObitMem allocation */
    if (ObitMemValid (in)) ObitMemFree (in);
    else g_free (in);
  } /* end deallocate */

  return NULL; /* new value for pointer */
} /* end ObitUnref */

/**
 * Determine if the input object is a member of the class whose
 * ClassInfo is pointed to by class,  or of a derived class.
 * Should also work for derived classes.
 * \param in Pointer to object to test.
 * \param class Pointer to ClassInfo structure of the class to be
 *              tested.
 * \return TRUE if member of class or a derived class, else FALSE.
 */
gboolean ObitIsA (gpointer in, gconstpointer class)
{
  const ObitClassInfo *inClass;

  /* if either input is null then it's not a match */
  if (in==NULL) return FALSE;
  if (class==NULL) return FALSE;

  /* Check ObitString */
  if (((Obit*)in)->ObitId != OBIT_ID) return FALSE;

  /* Loop back through the inheritance comparing ClassInfo 
     pointers */
  inClass = ((Obit*)in)->ClassInfo;
  while (inClass!=NULL) {
    if (inClass==class) return TRUE;  /* Found it */
    inClass = inClass->ParentClass;
  }
  return FALSE; /* if it got here it must not match */
} /* end ObitIsA  */

/**
 * Returns Obit magic blanking float value
 * This is adopted from AIPS and correcponds to the string 'INDE'
 * \return float magic value
 */
ofloat ObitMagicF (void)
{
  static union FBLANKequiv {
    gchar string[4];
    ofloat fblank;
  } FBLANK;
  FBLANK.string[0] = 'I'; 
  FBLANK.string[1] = 'N'; 
  FBLANK.string[2] = 'D'; 
  FBLANK.string[3] = 'E'; 
  
  return FBLANK.fblank;
} /* end ObitMagicF */

/**
 * Trims trailing blanks from string
 * Inserts Null after last non blank, at least one character always left.
 * \param str String to trim
 */
void ObitTrimTrail (gchar *str)
{
  olong i;
  for (i=strlen(str); i>0; i--) {
    if (str[i]==' ') str[i] = 0;
    if (str[i]!=0) break;
  }
} /* end ObitTrimTrail */

/**
 * Compare two strings allowing for various combinations of
 * blank filling and NULL termination.  
 * Blanks past the last non blank, non-NULL are considered insignificant
 * \param str1   First string to compare
 * \param str2   Second string to compare
 * \param maxlen Maximum number of characters to compare
 * \return True if all significant characters match, else False
 */
gboolean ObitStrCmp (gchar *str1, gchar *str2, olong maxlen)
{
  gboolean out=TRUE;
  olong i, last1, last2;

  /* Find last significant character */
  last1 = maxlen-1;
  for (i=maxlen-1; i>=0; i--) {
    if ((str1[i]==0) || (str1[i]==' ')) last1 = i-1;
  }
  last2 = maxlen-1;
  for (i=maxlen-1; i>=0; i--) {
    if ((str2[i]==0) || (str2[i]==' ')) last2 = i-1;
  }
  
  /* Compare */
  if (last1!=last2) return FALSE;
  for (i=0; i<last1; i++)
    if (str1[i]!=str2[i]) { out=FALSE; break;}
  return out;
} /* end ObitStrCmp */

/**
 * Returns a string with Today's date as yyyy-mm-dd
 * \return data string, should be g_freeed when done.
 */
gchar* ObitToday (void)
{
  gchar *out=NULL;
  struct tm *lp;
  time_t clock;
  olong datea[3];
 
  /* Get time since 00:00:00 GMT, Jan. 1, 1970 in seconds. */
  time (&clock);

  /* Convert to  broken-down time. */
  lp = localtime (&clock);

  /* to local arrays */
  datea[0] = lp->tm_year;
  if (datea[0]<1000) datea[0] += 1900; /* full year */
  datea[1] = lp->tm_mon+1; /* For some bizzare reason, month is 0-rel */
  datea[2] = lp->tm_mday;

  /* Create output */
  out = g_malloc(14);
  g_snprintf (out, 12, "%4.4d-%2.2d-%2.2d", datea[0], datea[1], datea[2]);

  return out;
} /* end ObitToday */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitClassInit (void)
{
  ObitClassInfo *ParentClass;

  if (myClassInfo.initialized) return;  /* only once */

  /* Initialize (recursively) parent class  */
  ParentClass = NULL; /* Base class */

  /* Class info. */
  myClassInfo.hasScratch    = FALSE; /* No scratch files */
  myClassInfo.ClassName     = g_strdup(myClassName);
  myClassInfo.ParentClass   = ParentClass;
  
  /* Set function pointers */
  ObitClassInfoDefFn ((gpointer)&myClassInfo);

  myClassInfo.initialized = TRUE; /* Now initialized */
} /* end ObitClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 * \param inClass   Pointer to ClassInfo structure of the class to be
 *                  filled.
 * \param callClass Pointer to ClassInfo of calling class
 */
void ObitClassInfoDefFn (gpointer inClass)
{
  ObitClassInfo *theClass = (ObitClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of in */
  g_assert (ObitInfoIsA(theClass, &myClassInfo));

  /* Initialize (recursively) parent class first */
  if (ParentClass!=NULL) 
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitGetClass;
  theClass->newObit       = (newObitFP)newObit;
  theClass->ObitCopy      = (ObitCopyFP)ObitCopy;
  theClass->ObitClone     = (ObitCloneFP)ObitClone;
  theClass->ObitRef       = (ObitRefFP)ObitRef;
  theClass->ObitUnref     = (ObitUnrefFP)ObitUnref;
  theClass->ObitIsA       = (ObitIsAFP)ObitIsA;
  theClass->ObitClear     = (ObitClearFP)ObitClear;
  theClass->ObitInit      = (ObitInitFP)ObitInit;
} /* end ObitClassDefFn */

/**
 * Determine if the class object is this or a derived class.
 * \param in Pointer to object to test.
 * \param class Pointer to ClassInfo structure of the class to be
 *              tested.
 * \return TRUE if test or a derived class, else FALSE.
 */
gboolean ObitInfoIsA (ObitClassInfo* class, ObitClassInfo* type)
{
  const ObitClassInfo *parent;

  /* if either input is null then it's not a match */
  if (class==NULL) return FALSE;
  if (type==NULL)  return FALSE;

  if (class==type) return TRUE;  /* Same */

  /* Loop back through the inheritance comparing ClassInfo 
     pointers */
  parent = class->ParentClass;
  while (parent!=NULL) {
    if (parent==type) return TRUE;  /* Found it */
    parent = parent->ParentClass;
  }
  return FALSE; /* if it got here it must not match */
} /* end ObitInfoIsA  */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Recursively initialized base classes before this class.
 * \param in Pointer to the object to initialize.
 */
static void ObitInit  (gpointer inn)
{
  Obit *in = inn;

  /* error checks */
  g_assert (in != NULL);
  
  /* set members in this class */
  in->ObitId = OBIT_ID;
  in->ReferenceCount = 1;

} /* end ObitInit */

/**
 * Deallocates member objects.
 * \param  in Pointer to the object to deallocate.
 */
static void ObitClear (gpointer inn)
{
  Obit *in = inn;

  /* error check */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  if (in->name){ g_free (in->name); in->name = NULL;}

} /* end ObitClear */


