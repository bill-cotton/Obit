/* $Id$     */
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
/*; Correspondence about this software should be addressed as follows:*/
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
#include "Obit.h"
#include "ObitOTFDesc.h"
#include "ObitTableFQ.h"
#include "ObitTableFQUtil.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitOTFDesc.c
 * ObitOTFDesc Obit GBT/OTF data descriptor class definition.
 * This contains information about the observations and the coordinates
 * in the image.
 */

/*--------------- File Global Variables  ----------------*/
/** name of the class defined in this file */
static gchar *myClassName = "ObitOTFDesc";

/**
 * ClassInfo global structure ObitIOClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitOTFDescClassInfo myClassInfo = {FALSE};

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/** truncate double to integer precision */
#ifndef AINT  
#define AINT(x) (odouble)((olong)(x))
#endif

/** Degrees to radians factor */
#ifndef DG2RAD  
#define DG2RAD G_PI / 180.0
#endif

/**  Radians to degrees factor */
#ifndef RAD2DG  
#define RAD2DG 180.0 / G_PI
#endif
/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitOTFDescInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitOTFDescClear (gpointer in);

/* Convert date to JD */
void ObitOTFDescDate2JD (const gchar* date, odouble *JD);

/** Private: Set Class function pointers. */
static void ObitOTFDescClassInfoDefFn (gpointer inClass);

/*---------------Public functions---------------------------*/
/**
 * Construct Object.
 * \return pointer to object created.
 */
ObitOTFDesc* newObitOTFDesc (gchar *name)
{
  ObitOTFDesc *out;

   /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitOTFDescClassInit();

  /* allocate structure */
  out = g_malloc0(sizeof(ObitOTFDesc));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

 /* set classInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitOTFDescInit((gpointer)out);

 return out;
} /* end newObitOTFDesc */

/**
 * Returns ClassInfo pointer for the class.
 * Initializes class if needed on first call.
 * \return pointer to the class structure.
 */
gconstpointer ObitOTFDescGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitOTFDescClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitOTFDescGetClass */

/**
 * Copy constructor.
 * The output descriptor will have the size and reference pixel
 * modified to reflect selection on the input, i.e. the output
 * descriptor describes the input.
 * \param in Pointer to object to be copied.
 * \param out Pointer to object to be written.  
 *            If NULL then a new structure is created.
 * \param err ObitErr error stack
 * \return Pointer to new object.
 */
ObitOTFDesc* ObitOTFDescCopy (ObitOTFDesc* in, ObitOTFDesc* out, 
			    ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  gboolean oldExist;
  gchar *outName;
  olong i, j;

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
    out = newObitOTFDesc(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /* initialize/copy */
  out->access  = in->access;
  out->nrecord = in->nrecord;
  out->naxis   = in->naxis;
  out->numDesc = in->numDesc;
  out->epoch   = in->epoch;
  out->equinox = in->equinox;
  out->JDObs   = in->JDObs;
  out->obsra   = in->obsra;
  out->obsdec  = in->obsdec;
  out->firstRec= in->firstRec;
  out->numRecBuff = in->numRecBuff;
  out->ncol    = in->ncol;
  out->numDesc = in->numDesc;
  out->beamSize= in->beamSize;
  out->diameter= in->diameter;
  out->OTFType = in->OTFType;
  for (i=0; i<OTFLEN_VALUE; i++) out->object[i] = in->object[i];
  for (i=0; i<OTFLEN_VALUE; i++) out->teles[i]  = in->teles[i];
  for (i=0; i<OTFLEN_VALUE; i++) out->obsdat[i] = in->obsdat[i];
  for (i=0; i<OTFLEN_VALUE; i++) out->origin[i] = in->origin[i];
  for (i=0; i<OTFLEN_VALUE; i++) out->date[i]   = in->date[i];
  for (i=0; i<OTFLEN_VALUE; i++) out->bunit[i]  = in->bunit[i];
  for (i=0; i<3; i++)            out->isort[i]  = in->isort[i];
  for (i=0; i<OTF_MAX_COL; i++) {
    out->colRepeat[i]  = in->colRepeat[i];
    for (j=0; j<OTFLEN_KEYWORD; j++) out->colType[i][j] = in->colType[i][j];
    for (j=0; j<OTFLEN_VALUE; j++)   out->colUnit[i][j] = in->colUnit[i][j];
  }

 /* loop over data axes */
  for (j=0; j<OTF_MAXDIM; j++) {
    out->inaxes[j] = in->inaxes[j];
    out->cdelt[j]  = in->cdelt[j];
    out->crota[j]  = in->crota[j];
    out->crpix[j]  = in->crpix[j];
    out->crval[j]  = in->crval[j];
    for (i=0; i<OTFLEN_KEYWORD; i++) out->ctype[j][i] = in->ctype[j][i];
  }

   /* Free any existing info members */
  if (in->info!=NULL) {
    if (out->info) out->info = ObitInfoListUnref (out->info); 
    out->info = ObitInfoListCopy (in->info);
  }

  /* index output */
  ObitOTFDescIndex (out);

  return out;
} /* end ObitOTFDescCopy */

/**
 * Copy descriptive material (i.e. things that don't define the structure).
 * \param in  Pointer to object to be copied.
 * \param out Pointer to object to be written.  
 * \param err ObitErr error stack
 */
void ObitOTFDescCopyDesc (ObitOTFDesc* in, ObitOTFDesc* out, 
			    ObitErr *err)
{
  olong i, j;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitIsA(out, &myClassInfo));
 
  /* initialize/copy */
  out->epoch   = in->epoch;
  out->equinox = in->equinox;
  out->JDObs   = in->JDObs;
  out->obsra   = in->obsra;
  out->obsdec  = in->obsdec;
  out->OTFType = in->OTFType;
  for (i=0; i<OTFLEN_VALUE; i++) out->object[i] = in->object[i];
  for (i=0; i<OTFLEN_VALUE; i++) out->teles[i]  = in->teles[i];
  for (i=0; i<OTFLEN_VALUE; i++) out->obsdat[i] = in->obsdat[i];
  for (i=0; i<OTFLEN_VALUE; i++) out->origin[i] = in->origin[i];
  for (i=0; i<OTFLEN_VALUE; i++) out->date[i]   = in->date[i];
  for (i=0; i<OTFLEN_VALUE; i++) out->bunit[i]  = in->bunit[i];
  for (i=0; i<3; i++)           out->isort[i]  = in->isort[i];
  for (i=0; i<OTF_MAX_COL; i++) {
    for (j=0; j<OTFLEN_KEYWORD; j++) out->colType[i][j] = in->colType[i][j];
    for (j=0; j<OTFLEN_VALUE; j++)   out->colUnit[i][j] = in->colUnit[i][j];
  }


  /* index output */
  ObitOTFDescIndex (out);

  return;
} /* end ObitOTFDescCopyDesc */

/**
 * Define indices for data
 * \param in Pointer to object.
 */
void ObitOTFDescIndex (ObitOTFDesc* in)
{
  olong numDesc, size, i;

  /* error check */
  g_assert (ObitIsA(in, &myClassInfo));

  /* data axis values */
  /* initialize */
  in->jlocdatawt = -1;
  in->jlocfeed   = -1;
  in->jlocs      = -1;
  in->jlocf      = -1;
  in->jlocstate  = -1;

  /* loop over axes looking for labels */
  size = 1;
  for (i=0; i<in->naxis; i++) {
    size *= MAX (1, in->inaxes[i]); /* how big */
    if (!strncmp (in->ctype[i], "DATAWT",  6)) in->jlocdatawt = i;
    if (!strncmp (in->ctype[i], "FEED",    4)) in->jlocfeed   = i;
    if (!strncmp (in->ctype[i], "STOKES",  6)) in->jlocs      = i;
    if (!strncmp (in->ctype[i], "FREQ",    4)) in->jlocf      = i;
    if (!strncmp (in->ctype[i], "STATE",   5)) in->jlocstate  = i;
  }
    /* Set increments in data array */
    in->incdatawt= 1;
    if (in->jlocdatawt==0) in->incdatawt= 2;  /* Data-Wt length 2 if present */
    in->incfeed  = 1;
    in->incs     = 1;
    in->incf     = 1;
    in->incstate = 1;
    if (in->jlocdatawt>0) { /* Oh Shit */
      g_error ("DATAWT axis if present MUST be first");
    }
    if (in->jlocfeed>0) { /* Product of dimensions of previous axes */
      for (i=0; i<in->jlocfeed; i++) in->incfeed *= in->inaxes[i];
    }
    if (in->jlocs>0) {
      for (i=0; i<in->jlocs; i++) in->incs *= in->inaxes[i];
    }
    if (in->jlocf>0) {
      for (i=0; i<in->jlocf; i++) in->incf *= in->inaxes[i];
    }
    if (in->jlocstate>0) {
      for (i=0; i<in->jlocstate; i++) in->incstate *= in->inaxes[i];
    }
  /* initialize data indices */
  in->iloct    = -1;
  in->ilocti   = -1;
  in->iloctar  = -1;
  in->ilocscan = -1;
  in->ilocra   = -1;
  in->ilocdec  = -1;
  in->ilocrot  = -1;
  in->iloccal  = -1;
  in->ilocdata = -1;

  /* loop over columns looking for labels */
  size = 0;
  numDesc = 0;
  for (i=0; i<in->ncol; i++) {
    size += MAX (1, in->colRepeat[i]); /* how big */
    if (!strncmp (in->colType[i], "TIME    ", 8)) {numDesc++; in->iloct   = i;}
    if (!strncmp (in->colType[i], "TIME_INT", 8)) {numDesc++; in->ilocti  = i;}
    if (!strncmp (in->colType[i], "TARGET  ", 8)) {numDesc++; in->iloctar = i;}
    if (!strncmp (in->colType[i], "SCAN    ", 8)) {numDesc++; in->ilocscan= i;}
    if (!strncmp (in->colType[i], "RA      ", 8)) {numDesc++; in->ilocra  = i;}
    if (!strncmp (in->colType[i], "DEC     ", 8)) {numDesc++; in->ilocdec = i;}
    if (!strncmp (in->colType[i], "ROTATE  ", 8)) {numDesc++; in->ilocrot = i;}
    if (!strncmp (in->colType[i], "CAL     ", 8)) {numDesc++; in->iloccal = i;}
    if (!strncmp (in->colType[i], "DATA    ", 8)) in->ilocdata = i;
  }

  /* total size in floats */
  in->lrec = size;
  
  /* Number of Descriptive parameters */
  in->numDesc = numDesc;

  /* Set Julian date if not there */
  if (in->JDObs<1.0) ObitOTFDescDate2JD (in->obsdat, &in->JDObs);
} /* end  ObitOTFDescIndex */

/**
 * Parse date string as ("yyyy-mm-dd" or "dd/mm/yy")
 * and convert to Julian Date.
 * Algorithm from ACM Algorithm number 199
 * This routine is good from 1 Mar 1900 indefinitely.
 * \param date [in] Date string
 * \param date [out] Julian date.
 */
void ObitOTFDescDate2JD (const gchar* date, odouble *JD)
{
  olong ymd[3], n;
  odouble ya, m, d, k, c;
  gchar temp[20];
  
  g_assert (date!=NULL);
  g_assert (JD!=NULL);

  /* Get date */
  if (date[2]=='/') { /* old format 'DD/MM/YY' */
    strncpy (temp, date, 8); temp[8] = 0;
    temp[2] = ' '; temp[5] = ' ';
    n = sscanf (temp, "%2d%3d%3d", &ymd[2], &ymd[1], &ymd[0]);
    /* guess century */
    if (ymd[0]>50) ymd[0]+=1900;
    else ymd[0]+=2000;
  } else {            /* new format 'YYYY-MM-DD' */
     strncpy (temp, date, 10); temp[10] = 0;
     temp[4] = ' '; temp[7] = ' ';
   n = sscanf (temp, "%4d%3d%3d", &ymd[0], &ymd[1], &ymd[2]);
  }
  /* if something went wrong return bad date */
  if (n!=3) {*JD = -1.0; return;}
  
  /* Convert to Days */
  ya = ymd[0];
  m  = ymd[1];
  d  = ymd[2];
  if (m<=2.0) {
    m  = AINT(m+9.0);
    ya = AINT(ya-1.0);
  } else {
    m  = AINT(m-3.0);
  }
  c = AINT(ya/100.0);
  ya = ya - 100.0 * c;
  k = AINT(146097.0 * c * 0.25) +
    AINT(1461.0 * ya * 0.25) +
    AINT((153.0 * m + 2.0) * 0.2) +
    d;

  /* Following good for > 20th cent. */
  *JD = AINT(k) + 1721118.50;
} /* end ObitOTFDescDate2JD */

/**
 * Convert a Julian date to a string in form "yyyy-mm-dd".
 * Apdapted from ASM Algorithm no. 199
 * \param date [in] Julian date.
 * \param date [out] Date string, Must be at least 11 characters.
 */
void ObitOTFDescJD2Date (odouble JD, gchar *date)
{
  olong  id, im, iy, ic;
  odouble j, y, d, m;

  g_assert (date!=NULL);

  /* error check */
  if (JD<1.0) {
    g_snprintf (date, 11, "BAD DATE");
    return;
  }
  
  j = AINT (JD + 0.50 - 1721119.0);
  y = AINT ((4.0*j - 1.00) / 146097.0);
  ic = y + 0.00010;
  j = 4.0*j - 1.00 - 146097.0*y;
  d = AINT (j * 0.250);
  j = AINT ((4.00*d + 3.00) / 1461.00);
  d = 4.00*d + 3.00 - 1461.00*j;
  d = AINT ((d+4.00) * 0.250);
  m = AINT ((5.00*d - 3.00) / 153.00);
  id = 5.00*d - 3.00 - 153.00*m;
  id = (id + 5) / 5;
  iy = j + 100*ic;
  im = m;
  if (im < 10) {
    im = im + 3;
  } else {
    im = im - 9;
  iy = iy + 1;
  }

  /* convert to string */
  g_snprintf (date, 11, "%4.4d-%2.2d-%2.2d", iy, im, id);
} /* end ObitOTFDescJD2Date */

/**
 * Convert OTF type to a descriptive string
 * Possibilities:\\
 * ``Unknown'',  ``DCR'': GBT DCR, ``SP'': GBT Spectral processor,
 * ``CCB'':CalTech Continuum Backend, ``PAR'':Penn Array Receiver
 * \param OTFType [in] enum
 * \param TString [out] readable representation (8 char+NULL)
 */
void ObitOTFDescType2String (ObitGBTOTFType OTFType, gchar *TString)
{
  if (OTFType==OBIT_GBTOTF_DCR)      strncpy (TString, "DCR     ",9);
  else if (OTFType==OBIT_GBTOTF_SP)  strncpy (TString, "SP      ",9);
  else if (OTFType==OBIT_GBTOTF_CCB) strncpy (TString, "CCB     ",9);
  else if (OTFType==OBIT_GBTOTF_PAR) strncpy (TString, "PAR     ",9);
  else strncpy (TString, "Unknown ",9);
} /* end  ObitOTFDescType2String */

/**
 * Convert a descriptive string to an OTF type
 * Possibilities:\\
 * ``Unknown'',  ``DCR'': GBT DCR, ``SP'': GBT Spectral processor,
 * ``CCB'':CalTech Continuum Backend, ``PAR'':Penn Array Receiver
 * \param OTFType [in] enum
 * \param TString [out] readable representation (8 char)
 */
ObitGBTOTFType ObitOTFDescString2Type (gchar *TString) {
  if (!strncmp(TString, "DCR",3)) return OBIT_GBTOTF_DCR;
  else if (!strncmp(TString, "SP",2)) return OBIT_GBTOTF_SP;
  else if (!strncmp(TString, "CCB",3)) return OBIT_GBTOTF_CCB;
  else if (!strncmp(TString, "PAR",3)) return OBIT_GBTOTF_PAR;
  return OBIT_GBTOTF_Unknown;
} /* end ObitOTFDescString2Type */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitOTFDescClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitOTFDescClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitOTFDescClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitOTFDescClassInfoDefFn (gpointer inClass)
{
  ObitOTFDescClassInfo *theClass = (ObitOTFDescClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitOTFDescClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitOTFDescClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitOTFDescGetClass;
  theClass->ObitClear     = (ObitClearFP)ObitOTFDescClear;
  theClass->ObitInit      = (ObitInitFP)ObitOTFDescInit;
  theClass->newObit       = (newObitFP)newObitOTFDesc;
  theClass->ObitCopy      = (ObitCopyFP)ObitOTFDescCopy;
  theClass->ObitClone     = NULL;

} /* end ObitOTFDescClassDefFn */

/*---------------Private functions--------------------------*/
/**
 * Creates empty member objects, initialize reference count.
 * Does (recursive) initialization of base class members before 
 * this class.
 * \param inn Pointer to the object to initialize.
 */
void ObitOTFDescInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitOTFDesc *in = inn;

  /* error checks */
  g_assert (in != NULL);
  
  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->ncol       = -1;
  in->nrecord    = 0;
  in->firstRec   = 0;
  in->numRecBuff = 0;
  in->info       = newObitInfoList();
  in->ncol       = 0;
  in->numDesc    = 0;
  in->OTFType    = OBIT_GBTOTF_Unknown;
} /* end ObitOTFDescInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 */
void ObitOTFDescClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitOTFDesc *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* free this class members */
  if (in->info) ObitInfoListUnref (in->info); in->info = NULL;
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);

} /* end ObitOTFDescClear */

