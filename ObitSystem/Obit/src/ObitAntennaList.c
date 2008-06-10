/* $Id$ */
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
/*;Correspondence about this software should be addressed as follows: */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/

#include "ObitAntennaList.h"

/*-------------- Obit: Merx mollis mortibus nuper ------------*/
/**
 * \file ObitAntennaList.c
 * ObitAntennaList class function definitions.
 *
 * This is a list of associated tables.
 */

/*------------------- File Global Variables - ----------------*/
/** name of the class defined in this file */
static gchar *myClassName = "ObitAntennaList";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo global structure ObitAntennaListClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitAntennaListClassInfo myClassInfo = {FALSE};

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitAntennaListInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitAntennaListClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitAntennaListClassInfoDefFn (gpointer inClass);


/*----------------------Public functions---------------------------*/
/**
 * Basic Constructor.
 * Initializes class if needed on first call.
 * \return the new object.
 */
ObitAntennaList* newObitAntennaList (gchar* name)
{
  ObitAntennaList* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitAntennaListClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitAntennaList));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

 /* set classInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitAntennaListInit((gpointer)out);

  return out;
} /* end newObitAntennaList */

/**
 * Returns ClassInfo pointer for the class.
 * Initializes class if needed on first call.
 * \return pointer to the class structure.
 */
gconstpointer ObitAntennaListGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitAntennaListClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitGetIOClass */

/**
 * Creates an ObitAntennaList of a given size.
 * \param name       An optional name for the object.
 * \param nant       Number of antennas (actually, the highest Antenna ID).
 * \param numpolcal  Number of Polarization calibration parameters
 * \return the new object.
 */
ObitAntennaList* ObitAntennaListCreate (gchar* name, olong nant, olong numpolcal)
{
  ObitAntennaList* out;
  olong i;
  gchar sname[31];

  /* Create basic structure */
  out = newObitAntennaList (name);

  /* create data array if wanted */
  if (nant<0) return out;

  /* Save information */
  out->number = nant;

  /* create array */
  out->ANlist = g_malloc0 (nant*sizeof(ObitAntenna*));

  /* Create elements of ObitAntenna list */
  for (i=0; i<out->number; i++) {
    g_snprintf (sname, 30, "Antenna  %d", i+1);
    out->ANlist[i] = newObitAntenna (sname);
    /* Polarization calibration arrays */
    out->ANlist[i]->FeedAPCal = g_malloc0(numpolcal*sizeof(ofloat));
    out->ANlist[i]->FeedBPCal = g_malloc0(numpolcal*sizeof(ofloat));
  }

  return out;
} /* end ObitAntennaListCreate */

/**
 * Make a copy of a object.
 * Parent class members are included but any derived class info is ignored.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitAntennaList* ObitAntennaListCopy  (ObitAntennaList *in, ObitAntennaList *out, 
				       ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  gboolean oldExist;
  gchar *outName;
  olong i;
  gchar *routine = "ObitAntennaListCopy";
  
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
    out = newObitAntennaList(outName);
    if (outName) g_free(outName); outName = NULL;
  }
  
  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);
  
  /*  copy this class */
  /* number of elements */
  out->number = in->number;
  
  /* Copy general information */
  out->polType   = in->polType;
  out->GSTIAT0   = in->GSTIAT0;
  out->RotRate   = in->RotRate;
  out->JD        = in->JD;
  out->ut1Utc    = in->ut1Utc;
  out->dataUtc   = in->dataUtc;
  out->numPoln   = in->numPoln;
  out->numPCal   = in->numPCal;
  out->polRefAnt = in->polRefAnt;
  out->numIF     = in->numIF;
  out->FreqID    = in->FreqID;
  for (i=0; i<3; i++) out->ArrayXYZ[i] = in->ArrayXYZ[i];
  for (i=0; i<2; i++) out->PolarXY[i]  = in->PolarXY[i];
  for (i=0; i<12; i++) out->TimSys[i]  = in->TimSys[i];
  for (i=0; i<12; i++) out->ArrName[i] = in->ArrName[i];
  /* Polarization R-L Phases */
  out->RLPhaseDiff = g_realloc (out->RLPhaseDiff, out->numIF*sizeof(ofloat));
  for (i=0; i<out->numIF; i++) out->RLPhaseDiff[i] = in->RLPhaseDiff[i];
  
  /* Reallocate list if needed */
  out->ANlist = g_realloc (out->ANlist, out->number*sizeof(ObitAntenna*));
  
  /* loop through list copying elements */
  for (i=0; i<out->number; i++) 
    out->ANlist[i] = ObitAntennaCopy (in->ANlist[i], out->ANlist[i], err);
  
  if (err->error) Obit_traceback_val (err, routine, in->name, out);
  
  return out;
} /* end ObitAntennaListCopy */

/**
 * Return corresponding type code, recognizes:
 * \li 'ORI-ELP ' Orientation-ellipticity (OBIT_UVPoln_ELORI)
 * \li 'APPROX  ' R/L Linear D term approximation (OBIT_UVPoln_Approx)
 * \li 'VLBI    ' R/L Linear D term approximation for resolved sources (OBIT_UVPoln_VLBI)
 * \li 'X-Y LIN ' X/Y Linear D term approximation (OBIT_UVPoln_XYLin)
 * \param type Polarization calibration type name
 * \return code, OBIT_UVPoln_Unknown -> not recognized, 
 *               OBIT_UVPoln_NoCal   -> no code (no poln cal)
 */
ObitUVPolCalType ObitAntennaListGetPolType (gchar* type)
{
  ObitUVPolCalType out = OBIT_UVPoln_Unknown;

  if (!strncmp (type, "    ",    4)) out = OBIT_UVPoln_NoCal;
  if (!strncmp (type, "ORI-ELP", 7)) out = OBIT_UVPoln_ELORI;
  if (!strncmp (type, "APPROX",  6)) out = OBIT_UVPoln_Approx;
  if (!strncmp (type, "VLBI",    4)) out = OBIT_UVPoln_VLBI;
  if (!strncmp (type, "X-Y LIN", 7)) out = OBIT_UVPoln_XYLin;
  return out;
} /* end ObitAntennaListGetPolType */

/**
 * Return source elevation in radians for a given antenna, time, and source.
 * From the AIPSish SOUELV.FOR
 * \param inAList  Antenna List
 * \param ant      Antenna number (1-rel)
 * \param time     Time in Days
 * \param Source   Source structure
 * \return elevation in radians
 */
ofloat ObitAntennaListElev (ObitAntennaList *inAList, olong ant, 
			    ofloat time, ObitSource *Source)
{
  ofloat elev = 0.0;
  ofloat decr, gst, lst, ha, along, alat;
  olong i, iant;
  odouble darg;

  /* Find antenna in list */
  iant = 0;
  for (i=0; i<inAList->number; i++) {
    if (inAList->ANlist[i]->AntID==ant) {iant = i; break;}
  }

  /* antenna location */
  alat  = inAList->ANlist[iant]->AntLat;
  along = inAList->ANlist[iant]->AntLong;

  /* declination in radians */
  decr = Source->DecApp * DG2RAD;

  /* Greenwich siderial time in radians */
  gst = inAList->GSTIAT0  + time*inAList->RotRate - inAList->dataIat;

  /* Local  siderial time */
  lst = gst + along;

  /* Hour angle in radians */
  ha = lst - Source->RAApp * DG2RAD;

  darg = sin (alat) * sin(decr) + cos (alat) * cos(decr) * cos (ha);
  elev = (1.570796327 - acos (MIN (darg, 1.000)));

  return elev;
} /* end ObitAntennaListElev */

/**
 * Return source azimuth in radians for a given antenna, time, and source.
 * From the AIPSish SOUELV.FOR
 * \param inAList  Antenna List
 * \param ant      Antenna number (1-rel)
 * \param time     Time in Days
 * \param Source   Source structure
 * \return azimuth in radians
 */
ofloat ObitAntennaListAz (ObitAntennaList *inAList, olong ant, 
			  ofloat time, ObitSource *Source)
{
  ofloat az = 0.0;
  ofloat decr, gst, lst, ha, along, alat;
  olong i, iant;
  odouble darg, darg2, daz;

  /* Find antenna in list */
  iant = 0;
  for (i=0; i<inAList->number; i++) {
    if (inAList->ANlist[i]->AntID==ant) {iant = i; break;}
  }

  /* antenna location */
  alat  = inAList->ANlist[iant]->AntLat;
  along = inAList->ANlist[iant]->AntLong;

  /* declination in radians */
  decr = Source->DecApp * DG2RAD;

  /* Greenwich siderial time in radians */
  gst = inAList->GSTIAT0  + time*inAList->RotRate - inAList->dataIat;

  /* Local  siderial time */
  lst = gst + along;

  /* Hour angle in radians */
  ha = lst - Source->RAApp * DG2RAD;

  /* Compute azimuth */
  darg  = sin(decr)*cos(alat) - cos(decr)*sin(alat)*cos(ha);
  darg2 = cos(decr) * sin(ha);
  daz = atan2 (darg, darg2);
  daz = fmod (daz, (2.0*G_PI));
  if (daz<0.0) daz += 2.0*G_PI;
  az = (ofloat)daz;

  return az;
} /* end ObitAntennaListAz */

/**
 * Return source Parallactic angle in radians for a given antenna, time, and source.
 * \param inAList  Antenna List
 * \param ant      Antenna number (1-rel)
 * \param time     Time in Days
 * \param Source   Source structure
 * \return paralactic angle in radians
 */
ofloat ObitAntennaListParAng (ObitAntennaList *inAList, olong ant, 
			      ofloat time, ObitSource *Source)
{
  ofloat parAng = 0.0;
  ofloat decr, gst, lst, ha, along, alat;
  olong i, iant;

  /* Find antenna in list */
  iant = 0;
  for (i=0; i<inAList->number; i++) {
    if (inAList->ANlist[i]->AntID==ant) {iant = i; break;}
  }

  /* antenna location */
  alat  = inAList->ANlist[iant]->AntLat;
  along = inAList->ANlist[iant]->AntLong;

  /* declination in radians */
  decr = Source->DecApp * DG2RAD;

  /* Greenwich siderial time in radians */
  gst = inAList->GSTIAT0  + time*inAList->RotRate - inAList->dataIat;

  /* Local  siderial time */
  lst = gst + along;

  /* Hour angle in radians */
  ha = lst - Source->RAApp * DG2RAD;

  parAng = atan2 (cos(alat) * sin(ha), 
		  (sin(alat)*cos(decr) - cos(alat)*sin(decr)*cos(ha)));
  return parAng;
} /* end ObitAntennaListParAng */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitAntennaListClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitAntennaListClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitAntennaListClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitAntennaListClassInfoDefFn (gpointer inClass)
{
  ObitAntennaListClassInfo *theClass = (ObitAntennaListClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitAntennaListClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitAntennaListClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitAntennaListGetClass;
  theClass->newObit       = (newObitFP)newObitAntennaList;
  theClass->ObitCopy      = (ObitCopyFP)ObitAntennaListCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitAntennaListClear;
  theClass->ObitInit      = (ObitInitFP)ObitAntennaListInit;

} /* end ObitAntennaListClassDefFn */


/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Does (recursive) initialization of base class members before 
 * this class.
 * \param inn Pointer to the object to initialize.
 */
void ObitAntennaListInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitAntennaList *in = inn;

  /* error checks */
  g_assert (in != NULL);
  
  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->ANlist      = NULL;
  in->RLPhaseDiff = NULL;
  in->number = 0;
} /* end ObitAntennaListInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 */
void ObitAntennaListClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitAntennaList *in = inn;
  olong i;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* free this class members */
  /* loop through list copying elements */
  for (i=0; i<in->number; i++) 
    in->ANlist[i] = ObitAntennaUnref (in->ANlist[i]);

  /* delete members */
  if (in->ANlist) g_free(in->ANlist); in->ANlist = NULL;
  if (in->RLPhaseDiff) g_free(in->RLPhaseDiff); in->RLPhaseDiff = NULL;
 
 /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);

} /* end ObitAntennaListClear */


  
