/* $Id: ObitAIPSObject.c,v 1.3 2007/08/31 17:24:03 bcotton Exp $ */
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
/*; Correspondence about this software should be addressed as follows: */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
#include <sys/types.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "ObitAIPS.h"
#include "ObitAIPSObject.h"

/*-------- ObitIO: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitAIPSObject.c
 * ObitAIPSObject module function definitions.
 */

/*------------------ Structures -----------------------------*/
  /** Maximum number of AIPS Classes */
#ifndef MAXCLS
#define MAXCLS 20
#endif

  /** Maximum number of virtual keywords per class */
#ifndef MAXVKW
#define MAXVKW 50
#endif

  /** Structure to contain virtual keywords */
  typedef struct {
    /** Have I been initialized? */
    gboolean initialized;
    /** Number of defined virtual keywords. */
    oint   VKWNum[MAXCLS];
    /** Name of virtual keyword  */
    gchar              VKWTableKey[MAXCLS][MAXVKW][8];
    /**   category:  1 = in fixed portion of catalog header, pointer is
      *               pointer into type dependent array.  D values must be
      *               copied from R array
      *               2 = in keyword/value portion of catalog header, some
      *               restrictions apply (not more than 2 words of data).
      *               3 = Special derived keywords read access only.  Pointer
      *               specifies a class specific function.
      */
    oint               VKWTableCat[MAXCLS][MAXVKW];
    /** Offset to catalog header entry or function, -1=? not found  */
    oint               VKWTableOff[MAXCLS][MAXVKW];
    /**  data type: 1,2,3,4,5 for D, R, C, I, L data types of associated  data. */
    ObitAIPSObjectType VKWTableType[MAXCLS][MAXVKW];
    /** Size of first dimension (length of string)  */
    oint               VKWTableDm1[MAXCLS][MAXVKW];
    /** Size of second dimension  */
    oint               VKWTableDm2[MAXCLS][MAXVKW];
    /** Class name for virtual keyword. */
    gchar  VKWCls[MAXCLS][8];
  } ObitAIPSObjData;
/*-----------------File Globals ---------------------------*/

/**
 * Virtual keyword data 
 */
static ObitAIPSObjData myVKData = {FALSE};;


/*---------------Private function prototypes----------------*/
/** Private: Initialize virtual keywords for inputs class. */
static void INvini (oint *ierr);

/** Private: Initialize virtual keywords for image class */
static void IMvini (oint *ierr);

/** Private: Initialize virtual keywords for uv data class */
static void UVvini (oint *ierr);

/** Private: Initialize virtual keywords for table class */
static void Tabvini (oint *ierr);

/*---------------AIPS function prototypes----------------*/
void zointd_ ();
void zocrob_ (gchar *name, gchar *class, oint *objnum, oint *ierr);
void zodeob_ (oint *objnum);
void zofnob_ (gchar *name, oint *objnum);
void zocpob_ (gchar *namein, gchar *nameout, oint *ierr);
void zofnle_ (oint *objnum, gchar *keyword, oint *ierr);
void zoinle_ (oint *objnum, gchar *keyword, oint *type, oint *ndim,
              oint *dim, oint *ierr);
void zostdt_ (oint *objnum, gchar *keyword, oint *type, oint *ndim,
              oint *dim, gchar *data, oint *ierr);
void zofedt_ (oint *objnum, gchar *keyword, oint *type, oint *ndim,
              oint *dim, gchar *data, oint *ierr);
void zostct_ (oint *objnum, oint *catblk, oint *ierr);
void zofect_ (oint *objnum, oint *catblk, oint *ierr);
void zofenm_ (oint *objnum, gchar *name, gchar *class, oint *ierr);
/*---------------Public functions---------------------------*/

/**
 * Initialize AIPS Object manager
 * Assumes initialization in AIPS of AIPS structures 
 * \param ierr return code, 0=>OK
 */
void ObitAIPSObjectOBinit (oint *ierr)
{
  oint   loop, i;
  
  *ierr = 0;

  /* init class information. */
  myVKData.initialized = TRUE;

  /* virtual keywords. */
  for (loop= 0; loop<MAXCLS; loop++) {
    /* number of keywords */
    myVKData.VKWNum[loop] = 0;

    /* class names */
    for (i=0; i<8; i++) myVKData.VKWCls[loop][i] = ' ';
  } 

  /* input class (none) */
  INvini (ierr);
  if (*ierr != 0) return;

  /* image class */
  IMvini (ierr);
  if (*ierr != 0) return;

  /* uv data class */
  UVvini (ierr);
  if (*ierr != 0) return;

  /* Table class */
  Tabvini (ierr);
  if (*ierr != 0) return;

} /* end ObitAIPSObjectOBinit */

/**
 * Add a Catalog header keyword to the virtual keyword list
 * \param class   object class
 * \param keyword keyword to add
 * \param type    data type OBIT_AIPSObjectInt, OBIT_AIPSObjectRe, OBIT_AIPSObjectDP, OBIT_AIPSObjectCar
 * \param ierr return code, 0=>OK
 */
void ObitAIPSObjectOBvhkw (AIPSObjClass class, AIPSKey keyword, ObitAIPSObjectType type, oint *ierr)
{
  oint   clasno, i, j, loop;
  gchar Cname[9], Keyword[9];
  gboolean OK;
  
  *ierr = 2;
  
   /* Internally keywords are 8, blank filled characters - convert input */
  for (loop=0; loop<8; loop++) Keyword[loop] = ' '; Keyword[loop] = 0;
  for (loop=0; loop<strlen(keyword); loop++) Keyword[loop] = keyword[loop];
  
   /* Internally class names are 8, blank filled characters - convert input */
  for (loop=0; loop<8; loop++) Cname[loop] = ' '; Cname[loop] = 0;
  for (loop=0; loop<strlen(class); loop++) Cname[loop] = class[loop];
  
  /* find class number */
  clasno = -1;
  for (i= 0; i<MAXCLS; i++) {
    if (!strncmp(myVKData.VKWCls[i], Cname, 8)) {
      clasno = i;
      *ierr = 0;
      break;
    }
  }
  
  /* didn't find? */
  if (clasno<0) g_error ("Unknown Class %s", class);
  
  /* is it already defined? */
  j = myVKData.VKWNum[clasno];
  OK = FALSE;
  for (i= 0; i<j; i++) { 
    if (!strncmp(myVKData.VKWTableKey[clasno][i], Keyword, 8)) {
      OK = TRUE;
    }
  } 
  
  /* is there room in table? */
  if (myVKData.VKWNum[clasno] >= MAXVKW) {
    g_error ("Virtual keyword table full for class %s", class);
    
    /* add entry to table */
  } else {
    myVKData.VKWNum[clasno] = myVKData.VKWNum[clasno] + 1;
    i = myVKData.VKWNum[clasno] - 1;
    strncpy (myVKData.VKWTableKey[clasno][i], Keyword, 8);
  } 
  
  /* insert new parameters */
  myVKData.VKWTableCat[clasno][i] = 2;
  myVKData.VKWTableOff[clasno][i] = 1;
  myVKData.VKWTableType[clasno][i] = type;
  myVKData.VKWTableDm1[clasno][i] = 1;
  /* Character strings always 8 char */
  if (type == OBIT_AIPSObjectCar) myVKData.VKWTableDm1[clasno][i] = 8;
  myVKData.VKWTableDm2[clasno][i] = 1;
  
  *ierr = 0;
} /* end ObitAIPSObjectOBvhkw */

/**
 * See if keyword is an object dependent, virtual keyword.
 * \param objnum object slot number
 * \param keywrd keyword
 * \param keypnt index in tables  , set if keywrd found
 * \param ierr return code, 0=>found, -1=not found, >0 => error
 * \param err Error stack
 */
void ObitAIPSObjectOBkeyv (oint objnum, AIPSKey keywrd, oint *keypnt, olong *ierr, ObitErr *err)
{
  oint   loop, clasno;
  gchar cname[9], Keyword[9];
  gchar *routine = "ObitAIPSObjectOBkeyv";

  /* Internally keywords are 8, blank filled characters - convert input */
  for (loop=0; loop<8; loop++) Keyword[loop] = ' '; Keyword[loop] = 0;
  for (loop=0; loop<strlen(keywrd); loop++) Keyword[loop] = keywrd[loop];
  
  /* lookup class number (clasno) */
  ObitAIPSObjectOBclass (objnum, &clasno, cname, err);
  if (err->error) {
    *ierr = 1; /* Error */
    Obit_traceback_msg (err, routine, "AIPS");
  }
  
  /* look for keyword in table */
  *ierr = -1; /* not yet found */
  for (loop= 1; loop<=myVKData.VKWNum[clasno]; loop++) {
    *keypnt = loop;
    if (!strncmp (myVKData.VKWTableKey[clasno][loop], Keyword, 8)) {
      *ierr = 0;  /* found */
      return; /* found */
    }
  } 
  
  /* not found */
} /* end ObitAIPSObjectOBkeyv */

/**
 * Return description of a specified real (non-virtual) keyword
 * \param objnum object slot number
 * \param keywrd keyword
 * \param type   data type OBIT_AIPSObjectInt, OBIT_AIPSObjectRe, 
 *               OBIT_AIPSObjectDP, OBIT_AIPSObjectCar
 * \param dim    dimensionality of value (valuec)
 * \param err ObitErr for reporting errors.
 * \return TRUE if found, else FALSE
 */
gboolean ObitAIPSObjectOBinfo (oint objnum, AIPSKey keywrd, ObitAIPSObjectType *type, 
			   AIPSKeyDim dim, ObitErr *err)
{
  gboolean Found = FALSE;
  oint   ondim, loop, keypnt, clasno, ierr=0;
  gchar Keyword[9], tname[33], cname[9], tclass[9];
  gchar *routine = "ObitAIPSObjectOBinfo";
  
  /* error tests */
  g_assert(ObitErrIsA(err)); 
  if (err->error) return Found;  /* existing error condition */
  
  /* Internally keywords are 8, blank filled characters - convert input */
  for (loop=0; loop<8; loop++) Keyword[loop] = ' '; Keyword[loop] = 0;
  for (loop=0; loop<strlen(keywrd); loop++) Keyword[loop] = keywrd[loop];

 /* virtual keyword? */
  ObitAIPSObjectOBkeyv (objnum, Keyword, &keypnt, &ierr, err);
  if (ierr==0) { /* found virtual keyword  - set value */
    /* find class number */
    ObitAIPSObjectOBclass (objnum, &clasno, cname, err);
    if (err->error) Obit_traceback_val (err, routine, "AIPS", Found);
    
    /* get keyword type, dimension */
    *type  = myVKData.VKWTableType[clasno][keypnt];
    dim[0] = myVKData.VKWTableDm1[clasno][keypnt];
    dim[1] = myVKData.VKWTableDm2[clasno][keypnt];
    Found = TRUE;
    return Found;
  } /* end virtual keyword */
  
  /* Try real keyword */
  zoinle_ (&objnum, Keyword, (oint*)type, &ondim, dim, &ierr);
  Found = ierr==0;
  if (ierr == 4) {
    Obit_log_error(err, OBIT_Error, "Bad object number" );
  } else if (ierr == 5) {
    Obit_log_error(err, OBIT_Error, "Object does not exist" );
  } 

  /* error message */
  if (err->error) { 
    /* Look up object name */
    strcpy (tname, "unknown object");
    zofenm_ (&objnum, tname, tclass, &ierr);
    /* trailing NULLs */
    tname[32] = 0; tclass[8] = 0;
    Obit_log_error(err, OBIT_Error, "Problem with object %s keyword %s", tname, keywrd);
  }
  
  return Found;
} /* end ObitAIPSObjectOBinfo */

/**
 * Fetch the value (array) for a specified real (non-virtual) keyword
 * \param objnum object slot number
 * \param keywrd keyword
 * \param type   data type OBIT_AIPSObjectInt, OBIT_AIPSObjectRe, 
 *               OBIT_AIPSObjectDP, OBIT_AIPSObjectCar
 * \param dim    dimensionality of value (valuec)
 * \param value  numeric data value (array)
 * \param valuec character data array
 * \param err ObitErr for reporting errors.
 */
void ObitAIPSObjectOBrget (oint objnum, AIPSKey keywrd, ObitAIPSObjectType *type, 
			   AIPSKeyDim dim, gpointer value, gchar *valuec, 
			   ObitErr *err)
{
  oint   otype, ondim, loop, ierr=0;
  AIPSKeyDim odim;
  gchar Keyword[9], tname[33], tclass[9];
  
  /* error tests */
  g_assert(ObitErrIsA(err)); 
  if (err->error) return;  /* existing error condition */
  
  /* Internally keywords are 8, blank filled characters - convert input */
  for (loop=0; loop<8; loop++) Keyword[loop] = ' '; Keyword[loop] = 0;
  for (loop=0; loop<strlen(keywrd); loop++) Keyword[loop] = keywrd[loop];

  /* find type, check existance */
  zoinle_ (&objnum, Keyword, &otype, &ondim, odim, &ierr);
  if (ierr == 1) {
    Obit_log_error(err, OBIT_Error, " Keyword %s not found", Keyword);
  } else if (ierr == 4) {
    Obit_log_error(err, OBIT_Error, "Bad object number" );
  } else if (ierr == 5) {
    Obit_log_error(err, OBIT_Error, "Object does not exist" );
  } 
  
  /* fetch value if found */
  if (!err->error) {
    /* numeric or character different */
    if (otype == OBIT_AIPSObjectCar) { /* character */
      zofedt_ (&objnum, Keyword, (oint*)type, &ondim, dim, valuec, &ierr);
    } else { /* numeric */
      zofedt_  (&objnum, Keyword, (oint*)type, &ondim, dim, (gchar*)value, &ierr);
    } 
  }

  /* error test */
  if (ierr!=0) Obit_log_error(err, OBIT_Error, "Error fetching value" );
  
  /* error message */
 if (err->error) { 
   /* Look up object name */
   strcpy (tname, "unknown object");
   zofenm_ (&objnum, tname, tclass, &ierr);
   /* trailing NULLs */
   tname[32] = 0; tclass[8] = 0;
   Obit_log_error(err, OBIT_Error, "Problem with object %s keyword %s", tname, keywrd);
 }
} /* end ObitAIPSObjectOBrget */

/**
 * Associate an object slot with an object name
 * \param name  Name of object to create
 * \param class Class name of new object
 * \param err ObitErr for reporting errors.
 */
void ObitAIPSObjectOBcrea (AIPSObj name, AIPSObjClass class, ObitErr *err)
{
  oint   loop, objnum, clsnum, ierr=0;
  gchar Class[9];
  gboolean OK=FALSE;

  /* error tests */
  g_assert(ObitErrIsA(err)); 
  if (err->error) return;  /* existing error condition */
  
   /* Internally class names are 8, blank filled characters - convert input */
  for (loop=0; loop<8; loop++) Class[loop] = ' '; Class[loop] = 0;
  for (loop=0; loop<strlen(class); loop++) Class[loop] = class[loop];
  
  /* create object storage */
  zocrob_ (name, Class, &objnum, &ierr);
  if (ierr == 1) {
    Obit_log_error(err, OBIT_Error, "AIPS Object table full");
  } else if (ierr == 2) {
    Obit_log_error(err, OBIT_Error, "AIPS Memory allocation failed");
  } else  if (ierr != 0) {
    Obit_log_error(err, OBIT_Error, "AIPS object creation failed");
  } 

  if (!err->error) {
    /* is class already defined? */
    OK = FALSE;
    clsnum = -1;
    for (loop= 0; loop<MAXCLS; loop++) {
     if (!strncmp(myVKData.VKWCls[loop], Class, 8)) {
       clsnum = loop;
       OK = TRUE;
       break;
     }
    } 
  }
  if (!OK) {
    /* no, find one to assign. */
    clsnum = -1;
    for (loop= 0; loop<MAXCLS; loop++) {
      clsnum = loop;
      if (!strncmp(myVKData.VKWCls[loop], "        ", 8)) {
	/* set class name */
	strncpy (myVKData.VKWCls[clsnum], Class, 8);
	break;
      } 
    }
    
    if (clsnum<0) {
      /* no more classes available */
      Obit_log_error(err, OBIT_Error, "No more AIPS classes can be defined");
    }
  } /* end new class */

  /* error */
  if (err->error) {
    Obit_log_error(err, OBIT_Error, "Problem creating %s", name);
  }
} /* end ObitAIPSObjectOBcrea */

/**
 * Free the object slot associated with an object.
 * \param name Name of object to free
 * \param err ObitErr for reporting errors.
 */
void ObitAIPSObjectOBfree (AIPSObj name, ObitErr *err)
{
  oint   objnum;
  
  /* error tests */
  g_assert(ObitErrIsA(err)); 
  if (err->error) return;  /* existing error condition */
  
  /* find slot number */
  ObitAIPSObjectOBname (name, &objnum, err);
  if (err->error) return;

  /* destroy */
  zodeob_ (&objnum);
  
} /* end ObitAIPSObjectOBfree */

/**
 *  Look up the object slot number of object with name "name"
 * \param name   Name of object to lookup
 * \param objnum object slot number
 * \param err ObitErr for reporting errors.
 */
void  ObitAIPSObjectOBname (AIPSObj name, oint *objnum, ObitErr *err)
{
  oint len, loop;
  gchar Name[33];

  /* error tests */
  g_assert(ObitErrIsA(err)); 
  if (err->error) return;  /* existing error condition */
  
   /* Internally names are 32 char, blank filled characters - convert input */
  for (loop=0; loop<32; loop++) Name[loop] = ' '; Name[loop] = 0;
  len = MIN (32, strlen(name));
  for (loop=0; loop<len; loop++) Name[loop] = name[loop];
  
  /* look up */
  zofnob_ (Name, objnum);
  if (*objnum > 0) return;
  
  /* not found */
  Obit_log_error(err, OBIT_Error, "AIPS object %s not found", name);
  *objnum = -1;
} /* end  ObitAIPSObjectOBname */

/**
 * Look up the class number and name of object number objnum. 
 * \param objnum  object slot number
 * \param classno Class number
 * \param name    name of class of object
 * \param err ObitErr for reporting errors.
 */
void  ObitAIPSObjectOBclass (oint objnum, oint *classno, AIPSObjClass name, ObitErr *err)
{
  olong   loop, clasno=0, ierr=0;
  gchar tname[33], tclass[9];
  gboolean OK;

  /* error tests */
  g_assert(ObitErrIsA(err)); 
  if (err->error) return;  /* existing error condition */
  
  /* get object class name */
  strcpy (tname, "unknown object");
  zofenm_ (&objnum, tname, tclass, &ierr);
  /* trailing NULLs */
  tname[32] = 0; tclass[8] = 0;
  if (ierr == 4) {
    Obit_log_error(err, OBIT_Error, "Bad AIPS object number");
  } else if (ierr == 5) {
    Obit_log_error(err, OBIT_Error, "AIPS object does not exist");
  } else if (ierr != 0) {
    Obit_log_error(err, OBIT_Error, "Error finding AIPS class");
  } 

   if (!err->error) {
     /* find class in Tables */
     OK = FALSE;
     clasno = -1;
     for (loop= 0; loop<MAXCLS; loop++) {
      if (!strncmp(myVKData.VKWCls[loop], tclass, 8)) {
	clasno = loop;
	OK = TRUE;
	break;
      }
     }

     if (clasno<0) {
       /* not found */
       ierr = 1;
       clasno = -1;
       Obit_log_error(err, OBIT_Error, "AIPS Class %s not defined", tclass);
     }
   }

   /* Set return class name, number */
   strncpy (name, tclass, 8);
   *classno = clasno;

   /* error */
   if (err->error) {
     Obit_log_error(err, OBIT_Error, "Problem with %s ", tname);
   }
} /* end  ObitAIPSObjectOBclass */

/**
 * Save an entry in an object creating it if necessary.
 * \param objnum object slot number
 * \param keywrd keyword
 * \param type   data type OBIT_AIPSObjectInt, OBIT_AIPSObjectRe, OBIT_AIPSObjectDP, OBIT_AIPSObjectCar
 * \param dim    dimensionality of value (valuec)
 * \param value  numeric data value (array)
 * \param valuec character data array
 * \param err ObitErr for reporting errors.
 */
void ObitAIPSObjectOBput (oint objnum, AIPSKey keywrd, ObitAIPSObjectType type, 
			  AIPSKeyDim dim, gpointer value, gchar *valuec, 
			  ObitErr *err)
{
  oint   nval, keypnt, clasno, catpnt, dim1, dim2, ndim, nwdpdp, loop, ierr = 0;
  union ObitAIPSCatEquiv catblk;
  gchar cname[9], tname[33], tclass[9], Keyword[9];
  gchar *routine = "ObitAIPSObjectOBput";

 /* error tests */
  g_assert(ObitErrIsA(err)); 
  if (err->error) return;  /* existing error condition */

   /* Internally keywords are 8, blank filled characters - convert input */
  for (loop=0; loop<8; loop++) Keyword[loop] = ' '; Keyword[loop] = 0;
  for (loop=0; loop<strlen(keywrd); loop++) Keyword[loop] = keywrd[loop];
  
 /* size of value */
  dim1 = dim[0];
  dim2 = dim[1];
  nwdpdp = sizeof(odouble) / sizeof(ofloat); /* how big is a double? */
  nval = dim1 * dim2;
  if (type == OBIT_AIPSObjectDP) nval = nval * nwdpdp;

  /* virtual keyword? */
  ObitAIPSObjectOBkeyv (objnum, Keyword, &keypnt, &ierr, err);

  if (ierr==0) { /* found virtual keyword - get value */

    /* get catblk */
    zofect_ (&objnum, (oint*)&catblk, &ierr);
    if (ierr != 0) {
      Obit_log_error(err, OBIT_Error, "Problem %d reading AIPS catalog header", ierr);
    } 
   
    /* find class number */
    ObitAIPSObjectOBclass (objnum, &clasno, cname, err);
    if (err->error) Obit_traceback_msg (err, routine, "AIPS");

    /* catalog header? */
    if (myVKData.VKWTableCat[clasno][keypnt] == 1) {

      /* check type */
      if (type != myVKData.VKWTableType[clasno][keypnt]) {
	ierr = 2;
	Obit_log_error(err, OBIT_Error, "Wrong AIPS keyword data type");
      } 
      /* check dimensions */
      if ((dim1 != myVKData.VKWTableDm1[clasno][keypnt])  || 
	  (dim2 != myVKData.VKWTableDm2[clasno][keypnt])) {
	ierr = 2;
	Obit_log_error(err, OBIT_Error, "Wrong AIPS keyword value dimensionality");
      } 
    
      /* catalog header pointer */
      catpnt = myVKData.VKWTableOff[clasno][keypnt];

      /* Check offset, -1 => not found */
      if (catpnt<0) {
	ierr = 2;
	Obit_log_error(err, OBIT_Error, "Keyword %s NOT in header", Keyword);
      }
	
      /* adjust for double */
      if (type == OBIT_AIPSObjectDP) catpnt = (catpnt-1) * nwdpdp + 1;
	

      /* Save if no problem */
      if (!err->error) {
	/* get data: character */
	if (type == OBIT_AIPSObjectCar) {
	  g_memmove ((gchar*)catblk.itg[catpnt], valuec, nval);
	  
	  /* other - numeric */
	} else { /* other */
	  g_memmove ((gchar*)catblk.itg[catpnt], (gchar*)value, nval*sizeof(ofloat));
	} 
	
	/* store catblk */
	zostct_ (&objnum, (oint*)&catblk, &ierr);
	if (ierr != 0) {
	  Obit_log_error(err, OBIT_Error, "Problem %d writing AIPS catalog header", ierr);
	} 
      }
      
      /* catalog keyword? */
    } else if (myVKData.VKWTableCat[clasno][keypnt] == 2) {
      
      g_error ("Write keyword/value in header +++++ NYI+++++");
      
      /* derived keyword? */
    } else if (myVKData.VKWTableCat[clasno][keypnt] == 3) {
      
      g_error ("Cannot write AIPS derived keyword");
      
      /* unknown category */
    } else {
      g_error ("Unknown AIPS virtual keyword category");
    } 
    return;

  } else { /* no virtual keywords - try real keyword */
    ndim = 2;
    if (type == OBIT_AIPSObjectCar) { /* Character */
      zostdt_ (&objnum, Keyword, (oint*)&type, &ndim, dim, valuec, &ierr);
    } else { /* numeric */
      zostdt_ (&objnum, Keyword, (oint*)&type, &ndim, dim, (gchar*)value, &ierr);
    } 
    if (ierr == 1) {
      Obit_log_error(err, OBIT_Error, "%s: Bad AIPS type", routine);
    } else if (ierr == 2) {
      Obit_log_error(err, OBIT_Error, "%s: Bad number of dimension ", routine);
    } else if (ierr == 3) {
      Obit_log_error(err, OBIT_Error, "%s: Bad dimension array", routine);
    } else if (ierr == 4) {
      Obit_log_error(err, OBIT_Error, "%s: Bad AIPS object number ", routine);
    } else if (ierr == 5) {
      Obit_log_error(err, OBIT_Error, "%s: Object does not exist ", routine);
    } else if (ierr == 6) {
      Obit_log_error(err, OBIT_Error, "%s: Error allocating keyword storage", routine);
    } 
  } 
  if (!err->error) return;

  /* error message */
  if (err->error) { 
    /* Look up object name */
    strcpy (tname, "unknown object");
    zofenm_ (&objnum, tname, tclass, &ierr);
    /* trailing NULLs */
    tname[32] = 0; tclass[8] = 0;
    Obit_log_error(err, OBIT_Error, "%s: Problem with object %s keyword %s", routine,tname, Keyword);
  }
} /* end ObitAIPSObjectOBput */

/**
 *  Fetch the value (array) for a specified keyword
 * \param objnum  object slot number
 * \param keywrd  keyword
 * \param type    data type OBIT_AIPSObjectInt, OBIT_AIPSObjectRe, OBIT_AIPSObjectDP, OBIT_AIPSObjectCar
 * \param dim     dimensionality of value (valuec)
 * \param value   numeric data value (array)
 * \param valuec  character data array
 * \param err ObitErr for reporting errors.
 */
void ObitAIPSObjectOBget (oint objnum, AIPSKey keywrd, ObitAIPSObjectType *type, 
			  AIPSKeyDim dim,  gpointer value, gchar *valuec, 
			  ObitErr *err)
{
  oint   nval, keypnt, clasno, disk, cno, catpnt, dim1, dim2, ierr;
  oint loop, i, nwdpdp;
  gboolean  all;
  gchar cname[9], tname[33], tclass[9], Keyword[9];
  union ObitAIPSCatEquiv catblk;
  gchar *routine = "ObitAIPSObjectOBget";

 /* error tests */
  g_assert(ObitErrIsA(err)); 
  if (err->error) return;  /* existing error condition */
  
   /* Internally keywords are 8, blank filled characters - convert input */
  for (loop=0; loop<8; loop++) Keyword[loop] = ' '; Keyword[loop] = 0;
  for (loop=0; loop<strlen(keywrd); loop++) Keyword[loop] = keywrd[loop];
  
 nwdpdp = sizeof(odouble) / sizeof(ofloat); /* how big is a double? */
  ierr = 0;
 
  /* virtual keyword? */
  ObitAIPSObjectOBkeyv (objnum, Keyword, &keypnt, &ierr, err);

  if (ierr==0) { /* found virtual keyword  - set value */
    /* find class number */
    ObitAIPSObjectOBclass (objnum, &clasno, cname, err);
    if (err->error) Obit_traceback_msg (err, routine, "AIPS");
    
    /* get keyword type, dimension */
    *type = myVKData.VKWTableType[clasno][keypnt];
    dim1 = myVKData.VKWTableDm1[clasno][keypnt];
    dim2 = myVKData.VKWTableDm2[clasno][keypnt];
    
    /* size of value */
    nval = dim1 * dim2;
    if (*type == OBIT_AIPSObjectDP) nval = nval * nwdpdp;
    
    /* get catblk */
    zofect_ (&objnum, (oint*)catblk.itg, &ierr);
    if (ierr != 0) {
      Obit_log_error(err, OBIT_Error, "Problem %d reading AIPS catalog header", ierr);
    } 
  
    /* Fetch value if no error */
    if (!err->error)  { /* Read CATBLK OK */
      
      /* catalog header? */
      if (myVKData.VKWTableCat[clasno][keypnt] == 1) {
	
	/* catalog header pointer */
	catpnt = myVKData.VKWTableOff[clasno][keypnt];
	
	/* Check offset, -1 => not found */
	if (catpnt<0) {
	  ierr = 2;
	  Obit_log_error(err, OBIT_Error, "Keyword %s NOT in header", Keyword);
	}
	
	/* adjust for double */
	if (*type == OBIT_AIPSObjectDP) catpnt = (catpnt-1) * nwdpdp + 1;
	
      /* Get if no problem */
	if (!err->error) {
	  /* get data: character */
	  if (*type == OBIT_AIPSObjectCar) {
	    g_memmove (valuec, (gchar*)&catblk.itg[catpnt], nval);
	    
	    /* other - numeric */
	  } else { /* other */
	    g_memmove ((gchar*)value, (gchar*)&catblk.itg[catpnt], nval*sizeof(ofloat));
	  } 
	}
      
	/* catalog keyword? */
      } else if (myVKData.VKWTableCat[clasno][keypnt] == 2) {
	
	/* get class. */
	strcpy (tname, "unknown object");
	zofenm_ (&objnum, tname, tclass, &ierr);
	/* trailing NULLs */
	tname[32] = 0; tclass[8] = 0;
	
	/* get disk and cno number. */
	ObitAIPSObjectOBdskc (tname, &disk, &cno, err);
	if (err->error) Obit_traceback_msg (err, routine, "AIPS");
	
	/* read keyword/value in header +++++ NYI+++++ */
	g_error ("read keyword/value in header +++++ NYI+++++");
	
	/* ????????????????????????????????????????????????????????
	   to do:
	   
	   1) decide if uvdata or image from class
	   2) get descriptor from catblk
	   3) get value from infoList
	   4) copy to value or valuec
	*/
	
	/* derived keyword? */
      } else if (myVKData.VKWTableCat[clasno][keypnt] == 3) {
	dim[0] = 1;
	dim[1] = 1;
	dim[2] = 1;
	dim[3] = 1;
	dim[4] = 1;
	/* read derived keywod +++++ NYI+++++ */
	g_error ("read derived keyword +++++ NYI+++++");
	
	/* ????????????????????????????????????????????????????????
	   to do:
	   
	   1) decide if uvdata or image from class
	   2) get descriptor from catblk
	   3) get value from descriptor
	   4) copy to value or valuec
	*/
	
      } else {/* unknown category */
	g_error ("unknown virtual keyword category");
      } 
      
      /* save dimensionality */
      dim[0] = dim1;
      dim[1] = dim2;
      dim[2] = 1;
      dim[3] = 1;
      dim[4] = 1;
      return;
    } /* End of if no error reading CATBLK */
      
  } else { /* not virtual keyword, try real keyword. */
    ObitAIPSObjectOBrget (objnum, Keyword, type, dim, value, valuec, err);
    if (err->error) {
      Obit_log_error(err, OBIT_Error, "Problem accessing keyword %s", Keyword);
    } 
  }

  /* Patch AIPS bug past the highest dimension is a 0 dimension followed by
     random? values - replace with ones */
  all = FALSE;
  for (i=1; i<MAXAIPSOBJDIM; i++) {
    all = all || (dim[i]==0);
    if (all) dim[i] = 1;
  }
  
  if (!err->error) return;
  
  /* error message */
  if (err->error) { 
    /* Look up object name */
    strcpy (tname, "unknown object");
    zofenm_ (&objnum, tname, tclass, &ierr);
    /* trailing NULLs */
    tname[32] = 0; tclass[8] = 0;
    Obit_log_error(err, OBIT_Error, "%s: Problem with object %s keyword %s", routine, tname, Keyword);
  }
} /* end ObitAIPSObjectOBget */

/**
 * Return Disk and slot information for object.
 * \param name  Name of object to lookup
 * \param disk  AIPS disk number
 * \param cno   AIPS catalog slot number
 * \param err ObitErr for reporting errors.
 */
void ObitAIPSObjectOBdskc (AIPSObj name, oint *disk, oint *cno, ObitErr *err)
{
  oint   objnum;
  ObitAIPSObjectType type;
  gchar cdummy[1];
  AIPSKeyDim dim;
  gchar *routine = "ObitAIPSObjectOBdskc";

  /* find slot number */
  ObitAIPSObjectOBname (name, &objnum, err);
  if (err->error) Obit_traceback_msg (err, routine, "AIPS");

  /* find disk  and slot */
  ObitAIPSObjectOBrget (objnum, "DISK", &type, dim, (gpointer)disk, cdummy, err);
  if (err->error) Obit_traceback_msg (err, routine, "AIPS");
  
  /* find disk  and slot */
  ObitAIPSObjectOBrget (objnum, "CNO", &type, dim, (gpointer)cno, cdummy, err);
  if (err->error) Obit_traceback_msg (err, routine, "AIPS");
  
} /* end ObitAIPSObjectOBdskc */

/**
 * Return catalog header record for an object.
 * \param name Name of object
 * \param cat  AIPS catalog header block
 * \param err ObitErr for reporting errors.
 */
void ObitAIPSObjectOBhget (AIPSObj name, union ObitAIPSCatEquiv cat, ObitErr *err)
{
  olong  objnum, ierr=0;
  gchar *routine = "ObitAIPSObjectOBhget";

  /* error tests */
  g_assert(ObitErrIsA(err)); 
  if (err->error) return;  /* existing error condition */
  
  /* find slot number */
  ObitAIPSObjectOBname (name, &objnum, err);
  if (err->error) Obit_traceback_msg (err, routine, "AIPS");

  /* copy header */
  zofect_ (&objnum, (oint*)&cat.itg, &ierr);
  if (ierr!=0) {
    Obit_log_error(err, OBIT_Error, "%s: Problem %d saving AIPS CATBLK", routine, ierr);
  }
} /* end ObitAIPSObjectOBhget */

/**
 * Store catalog header record for an object
 * \param name Name of object
 * \param cat  AIPS catalog header block
 * \param err ObitErr for reporting errors.
 */
void ObitAIPSObjectOBhput (AIPSObj name, union ObitAIPSCatEquiv cat, ObitErr *err)
{
  olong   objnum, ierr=0;
  gchar *routine = "ObitAIPSObjectOBhget";

  /* error tests */
  g_assert(ObitErrIsA(err)); 
  if (err->error) return;  /* existing error condition */
  
  /* find slot number */
  ObitAIPSObjectOBname (name, &objnum, err);
  if (err->error) Obit_traceback_msg (err, routine, "AIPS");

  /* copy header */
  zostct_ (&objnum, (oint*)&cat.itg, &ierr);
  if (ierr!=0) {
    Obit_log_error(err, OBIT_Error, "%s: Problem %d saving AIPS CATBLK", routine, ierr);
  }
} /* end ObitAIPSObjectOBhput */

/**
 * Copies one object to another
 * \param namein Input AIPS object name
 * \param namout Output AIPS object name
 * \param err ObitErr for reporting errors.
 */
void ObitAIPSObjectOBcopy (AIPSObj namein, AIPSObj namout, ObitErr *err)
{
  oint ierr=0;
  gchar *routine = "ObitAIPSObjectOBcopy";

  /* error tests */
  g_assert(ObitErrIsA(err)); 
  if (err->error) return;  /* existing error condition */
  
  /* Copy AIPS object */
  zocpob_ (namein, namout, &ierr);
  if (ierr == 1) {
    Obit_log_error(err, OBIT_Error, "%s: Problem with input AIPS object", routine);
  } else  if (ierr == 2) {
    Obit_log_error(err, OBIT_Error, "%s: Problem creating new AIPS object", routine);
  } else if (ierr == 3) {
    Obit_log_error(err, OBIT_Error, "%s: Problem  with copying linked list", routine);
  } 
  if (!err->error) return;

  /* error */
  Obit_log_error(err, OBIT_Error, "%s: Error copying %s to %s", routine, namein, namout);
} /* end ObitAIPSObjectOBcopy */

/*---------------Private functions---------------------------*/


/**
 * Initialize virtual keywords for INPUTS class
 * \param ierr return code, 0=>OK
 */
static void INvini (oint *ierr)
{
  oint   loop, clsnum;
 
  *ierr = 0;

  /* find next class number */
  for (loop= 1; loop<=MAXCLS; loop++) {
    clsnum = loop;
    if (!strncmp(myVKData.VKWCls[loop-1], "        ", 8)) break;
  } 
  
  /* initialize to no virtual  keywords. */
  myVKData.VKWNum[clsnum-1] = 0;
  
  /* set class name */
  strncpy (myVKData.VKWCls[clsnum-1], "INPUTS  ", 8);
} /* end INvini */

/**
 * Initialize virtual keywords for IMAGE class 
 * \param ierr return code, 0=>OK
 */
static void IMvini (oint *ierr)
{
#define NOIMKW 36
  oint   loop, clsnum, numvkw, len, i;
 
  /* Keyword names */
  gchar *DataKey[NOIMKW] = {
    "OBJECT", "TELESCOP", "INSTRUME", "OBSERVER",
    "DATE-OBS", "DATE-MAP", "BUNIT", "NDIM", "NAXIS", "CTYPE",
    "CRVAL", "CDELT", "CRPIX", "CROTA", "EPOCH", "DATAMAX",
    "DATAMIN", "PRODUCT", "NITER", "BMAJ", "BMIN", "BPA", "VELREF",
    "ALTRVAL", "ALTRPIX", "OBSRA", "OBSDEC", "RESTFREQ", "XSHIFT",
    "YSHIFT", "NAMCLSTY", "IMSEQ", "USERNO", "EXTYPE", "EXTVER", "BLANK"};

  /* category (1 for catalog stuff */
  oint DataCat[NOIMKW] = {
    1, 1, 1, 1, 
    1, 1, 1, 1, 1, 1, 
    1, 1, 1, 1, 1, 1, 
    1, 1, 1, 1, 1, 1, 1, 
    1, 1, 1, 1, 1, 1, 
    1, 1, 1, 1, 1, 1, 1};

  /* data type 1=D, 2=R, 3=C, 4=I  */
  ObitAIPSObjectType DataType[NOIMKW] = {
    OBIT_AIPSObjectCar, OBIT_AIPSObjectCar, OBIT_AIPSObjectCar, OBIT_AIPSObjectCar, 
    OBIT_AIPSObjectCar, OBIT_AIPSObjectCar, OBIT_AIPSObjectCar, OBIT_AIPSObjectInt, 
    OBIT_AIPSObjectInt, OBIT_AIPSObjectCar, OBIT_AIPSObjectDP,  OBIT_AIPSObjectRe, 
    OBIT_AIPSObjectRe,  OBIT_AIPSObjectRe,  OBIT_AIPSObjectRe,  OBIT_AIPSObjectRe, 
    OBIT_AIPSObjectRe,  OBIT_AIPSObjectInt, OBIT_AIPSObjectInt, OBIT_AIPSObjectRe, 
    OBIT_AIPSObjectRe,  OBIT_AIPSObjectRe,  OBIT_AIPSObjectInt, OBIT_AIPSObjectDP, 
    OBIT_AIPSObjectRe,  OBIT_AIPSObjectDP,  OBIT_AIPSObjectDP,  OBIT_AIPSObjectDP, 
    OBIT_AIPSObjectRe,  OBIT_AIPSObjectRe,  OBIT_AIPSObjectCar, OBIT_AIPSObjectInt, 
    OBIT_AIPSObjectInt, OBIT_AIPSObjectCar, OBIT_AIPSObjectInt, OBIT_AIPSObjectRe };

  /* 1st dimension (number of dim, extensions hard coded) */
  oint  DataDm1[NOIMKW] = {
    8, 8, 8, 8, 8, 8, 8, 1, 7, 7, 7, 7, 7, 7, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 20, 1, 1, 2, 50, 1};
  /* 2nd dimension */
  oint  DataDm2[NOIMKW] = {
    1, 1, 1, 1, 1, 1, 1, 1, 1, 7, 1, 1, 1, 1, 1, 1, 1, 1,
    1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 50, 1, 1};

  *ierr = 0;

  /* find next class number */
  for (loop= 0; loop<MAXCLS; loop++) {
    clsnum = loop;
    if (!strncmp(myVKData.VKWCls[loop], "        ", 8)) break;
  } 

  /* set class name */
  strncpy (myVKData.VKWCls[clsnum], "IMAGE   ", 8);
  
  /* initialize virtual keywords. */
  numvkw = 0;
  for (loop= 0; loop<NOIMKW; loop++) {
    /* Blank fill keywords to 8 characters */
    for (i=0; i<8; i++) myVKData.VKWTableKey[clsnum][numvkw][i] = ' ';
    len = MIN (8, strlen(DataKey[loop]));
    for (i=0; i<len; i++) myVKData.VKWTableKey[clsnum][numvkw][i] = DataKey[loop][i];
    myVKData.VKWTableCat[clsnum][numvkw]  = DataCat[loop];
    myVKData.VKWTableType[clsnum][numvkw] = DataType[loop];
    myVKData.VKWTableDm1[clsnum][numvkw]  = DataDm1[loop];
    myVKData.VKWTableDm2[clsnum][numvkw]  = DataDm2[loop];
    myVKData.VKWTableOff[clsnum][numvkw]  = ObitAIPSCatOffset(DataKey[loop]); /* Offset in catalog header */
    
    numvkw = numvkw + 1;
  }

  /* how many virtual keywords */
  myVKData.VKWNum[clsnum] = numvkw;

} /* end IMvini */

/**
 * Initialize virtual keywords for Table  class
 * Setup same at UV data as this is more general than images and
 * most tables are attached to UV data.
 * \param ierr return code, 0=>OK
 */
static void Tabvini (oint *ierr)
{
#define NOTABKW 38
  oint   loop, clsnum, numvkw, len, i;
 
  /* Keyword names */
  gchar *DataKey[NOTABKW] = {
    "OBJECT", "TELESCOP", "INSTRUME", "OBSERVER",
    "DATE-OBS", "DATE-MAP", "BUNIT", "NDIM", "NAXIS", "CTYPE",
    "CRVAL", "CDELT", "CRPIX", "CROTA", "EPOCH", "GCOUNT",
    "NRPARM", "PTYPE", "SORTORD", "VELREF", "ALTRVAL", "ALTRPIX",
    "OBSRA", "OBSDEC", "RESTFREQ", "XSHIFT", "YSHIFT", "NAMCLSTY",
    "IMSEQ", "USERNO", "EXTYPE", "EXTVER", "LREC", "NCORR",
    "TYPEUVD", "BMAJ", "BMIN", "BPA"};

  /* category (1 for catalog stuff */
  oint DataCat[NOTABKW] = {
    1, 1, 1, 1, 
    1, 1, 1, 1, 1, 1, 
    1, 1, 1, 1, 1, 1, 
    1, 1, 1, 1, 1, 1, 1, 
    1, 1, 1, 1, 1, 1, 
    1, 1, 1, 1, 1, 1, 1, 1, 1};

  /* data type 1=D, 2=R, 3=C, 4=I  */
  ObitAIPSObjectType DataType[NOTABKW] = {
    OBIT_AIPSObjectCar, OBIT_AIPSObjectCar, OBIT_AIPSObjectCar, OBIT_AIPSObjectCar, 
    OBIT_AIPSObjectCar, OBIT_AIPSObjectCar, OBIT_AIPSObjectCar, OBIT_AIPSObjectInt, 
    OBIT_AIPSObjectInt, OBIT_AIPSObjectCar, OBIT_AIPSObjectDP,  OBIT_AIPSObjectRe, 
    OBIT_AIPSObjectRe,  OBIT_AIPSObjectRe,  OBIT_AIPSObjectRe,  OBIT_AIPSObjectInt, 
    OBIT_AIPSObjectInt, OBIT_AIPSObjectCar, OBIT_AIPSObjectCar, OBIT_AIPSObjectInt, 
    OBIT_AIPSObjectDP,  OBIT_AIPSObjectRe,  OBIT_AIPSObjectDP,  OBIT_AIPSObjectDP, 
    OBIT_AIPSObjectDP,  OBIT_AIPSObjectRe,  OBIT_AIPSObjectRe,  OBIT_AIPSObjectCar, 
    OBIT_AIPSObjectInt, OBIT_AIPSObjectInt, OBIT_AIPSObjectCar, OBIT_AIPSObjectInt, 
    OBIT_AIPSObjectInt, OBIT_AIPSObjectInt, OBIT_AIPSObjectCar, OBIT_AIPSObjectRe, 
    OBIT_AIPSObjectRe,  OBIT_AIPSObjectRe};

  /* 1st dimension (number of dim, extensions hard coded) */
  oint  DataDm1[NOTABKW] = {
      8, 8, 8, 8, 8, 8, 8, 1, 7, 8, 7, 7, 7, 7, 1, 1, 1, 8,
      2, 1, 1, 1, 1, 1, 1, 1, 1, 20, 1, 1, 2, 50, 1, 1, 2, 1, 1, 1};
  /* 2nd dimension */
  oint  DataDm2[NOTABKW] = {
  1, 1, 1, 1, 1, 1, 1, 1, 1, 7, 1, 1, 1, 1, 1, 1, 1,
  14, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 50, 1, 1, 1, 1, 1, 1, 1};

  *ierr = 0;

  /* find next class number */
  for (loop= 0; loop<MAXCLS; loop++) {
    clsnum = loop;
    if (!strncmp(myVKData.VKWCls[loop], "        ", 8)) break;
  } 

  /* set class name */
  strncpy (myVKData.VKWCls[clsnum], "TABLE   ", 8);
  
  /* initialize virtual keywords. */
  numvkw = 0;
  for (loop= 0; loop<NOTABKW; loop++) {
    /* Blank fill keywords to 8 characters */
    for (i=0; i<8; i++) myVKData.VKWTableKey[clsnum][numvkw][i] = ' ';
    len = MIN (8, strlen(DataKey[loop]));
    for (i=0; i<len; i++) myVKData.VKWTableKey[clsnum][numvkw][i] = DataKey[loop][i];
    myVKData.VKWTableCat[clsnum][numvkw]  = DataCat[loop];
    myVKData.VKWTableType[clsnum][numvkw] = DataType[loop];
    myVKData.VKWTableDm1[clsnum][numvkw]  = DataDm1[loop];
    myVKData.VKWTableDm2[clsnum][numvkw]  = DataDm2[loop];
    myVKData.VKWTableOff[clsnum][numvkw]  = ObitAIPSCatOffset(DataKey[loop]); /* Offset in catalog header */
    
    numvkw = numvkw + 1;
  }

  /* how many virtual keywords */
  myVKData.VKWNum[clsnum] = numvkw;
} /* end Tabvini */

/**
 * Initialize virtual keywords for UV data  class
 * \param ierr return code, 0=>OK
 */
static void UVvini (oint *ierr)
{
#define NOUVKW 38
  oint   loop, clsnum, numvkw, len, i;
 
  /* Keyword names */
  gchar *DataKey[NOUVKW] = {
    "OBJECT", "TELESCOP", "INSTRUME", "OBSERVER",
    "DATE-OBS", "DATE-MAP", "BUNIT", "NDIM", "NAXIS", "CTYPE",
    "CRVAL", "CDELT", "CRPIX", "CROTA", "EPOCH", "GCOUNT",
    "NRPARM", "PTYPE", "SORTORD", "VELREF", "ALTRVAL", "ALTRPIX",
    "OBSRA", "OBSDEC", "RESTFREQ", "XSHIFT", "YSHIFT", "NAMCLSTY",
    "IMSEQ", "USERNO", "EXTYPE", "EXTVER", "LREC", "NCORR",
    "TYPEUVD", "BMAJ", "BMIN", "BPA"};

  /* category (1 for catalog stuff */
  oint DataCat[NOUVKW] = {
    1, 1, 1, 1, 
    1, 1, 1, 1, 1, 1, 
    1, 1, 1, 1, 1, 1, 
    1, 1, 1, 1, 1, 1, 1, 
    1, 1, 1, 1, 1, 1, 
    1, 1, 1, 1, 1, 1, 1, 1, 1};

  /* data type 1=D, 2=R, 3=C, 4=I  */
  ObitAIPSObjectType DataType[NOUVKW] = {
    OBIT_AIPSObjectCar, OBIT_AIPSObjectCar, OBIT_AIPSObjectCar, OBIT_AIPSObjectCar, 
    OBIT_AIPSObjectCar, OBIT_AIPSObjectCar, OBIT_AIPSObjectCar, OBIT_AIPSObjectInt, 
    OBIT_AIPSObjectInt, OBIT_AIPSObjectCar, OBIT_AIPSObjectDP,  OBIT_AIPSObjectRe, 
    OBIT_AIPSObjectRe,  OBIT_AIPSObjectRe,  OBIT_AIPSObjectRe,  OBIT_AIPSObjectInt, 
    OBIT_AIPSObjectInt, OBIT_AIPSObjectCar, OBIT_AIPSObjectCar, OBIT_AIPSObjectInt, 
    OBIT_AIPSObjectDP,  OBIT_AIPSObjectRe,  OBIT_AIPSObjectDP,  OBIT_AIPSObjectDP, 
    OBIT_AIPSObjectDP,  OBIT_AIPSObjectRe,  OBIT_AIPSObjectRe,  OBIT_AIPSObjectCar, 
    OBIT_AIPSObjectInt, OBIT_AIPSObjectInt, OBIT_AIPSObjectCar, OBIT_AIPSObjectInt, 
    OBIT_AIPSObjectInt, OBIT_AIPSObjectInt, OBIT_AIPSObjectCar, OBIT_AIPSObjectRe, 
    OBIT_AIPSObjectRe,  OBIT_AIPSObjectRe};

  /* 1st dimension (number of dim, extensions hard coded) */
  oint  DataDm1[NOUVKW] = {
      8, 8, 8, 8, 8, 8, 8, 1, 7, 8, 7, 7, 7, 7, 1, 1, 1, 8,
      2, 1, 1, 1, 1, 1, 1, 1, 1, 20, 1, 1, 2, 50, 1, 1, 2, 1, 1, 1};
  /* 2nd dimension */
  oint  DataDm2[NOUVKW] = {
  1, 1, 1, 1, 1, 1, 1, 1, 1, 7, 1, 1, 1, 1, 1, 1, 1,
  14, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 50, 1, 1, 1, 1, 1, 1, 1};

  *ierr = 0;

  /* find next class number */
  for (loop= 0; loop<MAXCLS; loop++) {
    clsnum = loop;
    if (!strncmp(myVKData.VKWCls[loop], "        ", 8)) break;
  } 

  /* set class name */
  strncpy (myVKData.VKWCls[clsnum], "UVDATA  ", 8);
  
  /* initialize virtual keywords. */
  numvkw = 0;
  for (loop= 0; loop<NOUVKW; loop++) {
    /* Blank fill keywords to 8 characters */
    for (i=0; i<8; i++) myVKData.VKWTableKey[clsnum][numvkw][i] = ' ';
    len = MIN (8, strlen(DataKey[loop]));
    for (i=0; i<len; i++) myVKData.VKWTableKey[clsnum][numvkw][i] = DataKey[loop][i];
    myVKData.VKWTableCat[clsnum][numvkw]  = DataCat[loop];
    myVKData.VKWTableType[clsnum][numvkw] = DataType[loop];
    myVKData.VKWTableDm1[clsnum][numvkw]  = DataDm1[loop];
    myVKData.VKWTableDm2[clsnum][numvkw]  = DataDm2[loop];
    myVKData.VKWTableOff[clsnum][numvkw]  = ObitAIPSCatOffset(DataKey[loop]); /* Offset in catalog header */
    
    numvkw = numvkw + 1;
  }

  /* how many virtual keywords */
  myVKData.VKWNum[clsnum] = numvkw;
} /* end UVvini */

