/* $Id$        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2013                                               */
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

#include "ObitUVWCalc.h"
#include "ObitTableSUUtil.h"
#include "ObitTableANUtil.h"
#include "ObitPrecess.h"
#ifndef VELIGHT
#define VELIGHT 2.997924562e8
#endif

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitUVWCalc.c
 * ObitUVWCalc class function definitions.
 * This class is derived from the Obit base class.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitUVWCalc";

/** Function to obtain parent ClassInfo */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo structure ObitUVWCalcClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitUVWCalcClassInfo myClassInfo = {FALSE};

/*--------------- File Global Variables  ----------------*/


/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitUVWCalcInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitUVWCalcClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitUVWCalcClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitUVWCalc* newObitUVWCalc (gchar* name)
{
  ObitUVWCalc* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitUVWCalcClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitUVWCalc));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitUVWCalcInit((gpointer)out);

 return out;
} /* end newObitUVWCalc */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitUVWCalcGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitUVWCalcClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitUVWCalcGetClass */

/**
 * Make a deep copy of an ObitUVWCalc.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitUVWCalc* ObitUVWCalcCopy  (ObitUVWCalc *in, ObitUVWCalc *out, ObitErr *err)
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
    out = newObitUVWCalc(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  out->myData   = ObitUVCopy (in->myData, out->myData, err);
  out->SouList  = ObitSourceListCopy (in->SouList, out->SouList, err);
  out->mySource = ObitSourceCopy (in->mySource, out->mySource, err);
  out->AntList  = ObitAntennaListCopy (in->AntList, out->AntList, err);
  out->maxAnt   = in->maxAnt;
  if (out->antIndex) g_free(out->antIndex);
  out->antIndex = g_malloc0((out->maxAnt+2)*sizeof(olong));
  for (i=0; i<out->maxAnt; i++) out->antIndex[i] = in->antIndex[i];

  return out;
} /* end ObitUVWCalcCopy */

/**
 * Make a copy of a object but do not copy the actual data
 * This is useful to create an UVWCalc similar to the input one.
 * \param in  The object to copy
 * \param out An existing object pointer for output, must be defined.
 * \param err Obit error stack object.
 */
void ObitUVWCalcClone  (ObitUVWCalc *in, ObitUVWCalc *out, ObitErr *err)
{
  olong i;
  const ObitClassInfo *ParentClass;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitIsA(out, &myClassInfo));

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  Clone this class 
      Delete any old */
  out->myData   = ObitUVUnref (out->myData);
  out->SouList  = ObitSourceListUnref (out->SouList);
  out->mySource = ObitSourceUnref (out->mySource);
  out->AntList  = ObitAntennaListUnref (out->AntList);
  if (out->antIndex) g_free(out->antIndex);
  /* reference new */
  out->myData   = ObitUVRef (in->myData);
  out->SouList  = ObitSourceListRef (in->SouList);
  out->mySource = ObitSourceRef (in->mySource);
  out->AntList  = ObitAntennaListRef (in->AntList);
  out->maxAnt   = in->maxAnt;
  out->antIndex = g_malloc0((out->maxAnt+2)*sizeof(olong));
  for (i=0; i<out->maxAnt; i++) out->antIndex[i] = in->antIndex[i];

} /* end ObitUVWCalcClone */

/**
 * Creates an ObitUVWCalc 
 * Initializes various structures and modifies the AntennaList
 * \param name  An optional name for the object.
 * \param inUV  ObitUV for which u,v,ws are desired
 * \param err   Obit error stack object.
 * \return the new object.
 */
ObitUVWCalc* ObitUVWCalcCreate (gchar* name, ObitUV *inUV, ObitErr *err)
{
  ObitUVWCalc  *out=NULL;
  ObitTableSU  *SUTable=NULL;
  ObitTableAN  *ANTable=NULL;
  olong        i, ver, numPCal, numOrb, numIF, cnt;
  odouble      lambda, xx, yy, zz, Long;
  gchar *routine = "ObitUVWCalcCreate";
  
  /* error checks */
  if (err->error) return out;
  
  /* Create basic structure */
  out = newObitUVWCalc (name);

  /* Fully instantiate UV data if needed */
  ObitUVFullInstantiate (inUV, TRUE, err);
  if (err->error) Obit_traceback_val (err, routine, inUV->name, out);

  /* Save uv data */
  out->myData = ObitUVRef(inUV);

  /* Convert SU table into Source List if there is an SourceID parameter */
  if (inUV->myDesc->ilocsu>=0) {
    ver = 1;
    numIF = 0;
    SUTable = newObitTableSUValue (inUV->name, (ObitData*)inUV, &ver, numIF, 
				   OBIT_IO_ReadOnly, err);
    if (SUTable) out->SouList = ObitTableSUGetList (SUTable, err);
    if (err->error) Obit_traceback_val (err, routine, inUV->name, out);
    SUTable = ObitTableSUUnref(SUTable);
  } else {  /* Create a source object */
    out->mySource = newObitSource("Single");
    strncpy (out->mySource->SourceName, inUV->myDesc->object, MIN(20,UVLEN_VALUE));
    out->mySource->equinox = inUV->myDesc->equinox;
    out->mySource->RAMean  = inUV->myDesc->crval[inUV->myDesc->jlocr];
    out->mySource->DecMean = inUV->myDesc->crval[inUV->myDesc->jlocd];
    /* Compute apparent position */
    ObitPrecessUVJPrecessApp (inUV->myDesc, out->mySource);
  }

  /* Convert AN table into AntennaList */
  ver = 1;
  numPCal  = 0;
  numOrb   = 0;
  ANTable = newObitTableANValue (inUV->name, (ObitData*)inUV, &ver, 
				 numIF, numOrb, numPCal, OBIT_IO_ReadOnly, err);
  if (ANTable) out->AntList = ObitTableANGetList (ANTable, err);
  if (err->error) Obit_traceback_val (err, routine, inUV->name, out);
  ANTable = ObitTableANUnref(ANTable);

  /* maximum antenna number */
  out->maxAnt = 0;
  for (i=0; i<out->AntList->number; i++) out->maxAnt = MAX (out->maxAnt, out->AntList->ANlist[i]->AntID);

  /* (1-rel) Index into list */
  out->antIndex = g_malloc0((out->maxAnt+2)*sizeof(olong));
  for (i=0; i<out->maxAnt; i++) out->antIndex[i] = -1;
  for (i=0; i<out->AntList->number; i++) out->antIndex[out->AntList->ANlist[i]->AntID] = i;


  /* Average Longitude */
  cnt  = 0;
  Long = 0.0;
  for (i=0; i<out->maxAnt; i++) {
    out->AntList->ANlist[i]->AntLong = atan2 (out->AntList->ANlist[i]->AntXYZ[1], 
					      out->AntList->ANlist[i]->AntXYZ[0]);
    if (fabs(out->AntList->ANlist[i]->AntXYZ[1])>1.0) {
      cnt++;
      Long += out->AntList->ANlist[i]->AntLong;
    }
  }
  Long /= cnt;
  Long = -Long; /* Other way for rotation */

  /* Modify AntennaList for UVW calculations 
     Convert positions to lambda - rotate to frame of the array */
  lambda = VELIGHT/inUV->myDesc->freq;  /* Reference wavelength */
  for (i=0; i<out->AntList->number; i++) {
    xx = out->AntList->ANlist[i]->AntXYZ[0] / lambda;
    yy = out->AntList->ANlist[i]->AntXYZ[1] / lambda;
    zz = out->AntList->ANlist[i]->AntXYZ[2] / lambda;
    out->AntList->ANlist[i]->AntXYZ[0] = xx*cos(Long) - yy*sin(Long);
    out->AntList->ANlist[i]->AntXYZ[1] = xx*sin(Long) + yy*cos(Long);
    out->AntList->ANlist[i]->AntXYZ[2] = zz;
  }

  return out;
} /* end ObitUVWCalcCreate */

/**
 * Calculate the u,v,w for a given source, time, baseline
 * This is useful to create an UVWCalc similar to the input one.
 * \param in  The object with antenna/source information
 * \param time Time (days) wrt reference day for subA
 * \param SId  Source identifier
 * \param subA Subarray number, only one supported for now
 * \param ant1 First antenna number of baseline
 * \param ant2 Second antenna number of baseline
 * \param uvw  [out] [u,v,w] array in wavelengths at reference frequency
 * \param err Obit error stack object.
 */
void ObitUVWCalcUVW (ObitUVWCalc *in, ofloat time, olong SId,
		     olong subA, olong ant1, olong ant2, ofloat *uvw, 
		     ObitErr *err)
{
  ofloat bl[3], u, v, sinha, cosha, vw;
  olong i, ia1, ia2;
  ObitSource *source=NULL;
  odouble RAR, AntLst, HrAng=0.0, BLLong, dRa, RAOff, DecOff;
  gchar *routine = "ObitUVWCalcUVW";

  /* error checks */
  if (err->error) return;
  Obit_return_if_fail((subA<=1), err,
		      "%s: Only one subarray supported", routine);
  Obit_return_if_fail(((ant1>=0) && (ant1<=in->maxAnt)), err,
		      "%s: Ant1 %d out of bounds", routine, ant1);
  Obit_return_if_fail(((ant2>=0) && (ant2<=in->maxAnt)), err,
		      "%s: Ant2 %d out of bounds", routine, ant2);

  /* Initialize output */
  uvw[0] = uvw[1] = uvw[2] = 0.0;

  /* Done for auto correlations */
  if (ant1==ant2) return;

  /* Antenna indices */
  ia1 = in->antIndex[ant1];
  ia2 = in->antIndex[ant2];

  /* New source?  */
  if (in->curSID!=SId) {
    in->curSID = SId;
    /* Find in Source List */
    if (in->SouList) {
      in->curSource = NULL;
      for (i=0; i<in->SouList->number; i++) {
	if (in->SouList->SUlist[i]->SourID==SId) {
	  in->curSource = in->SouList->SUlist[i];
	  break;}
      }
    } else { /* Single source */
      in->curSource = in->mySource;
    }
    
    /* Check */
    Obit_return_if_fail((in->curSource!=NULL), err,
			"%s: Source ID %d not found", routine, SId);
    
    /* Get rotation to get v north at standard epoch - 
       precess posn. 10" north */
    source = newObitSource("Temp");
    source->RAMean  = in->curSource->RAMean;
    source->DecMean = in->curSource->DecMean + 10.0/3600.0;
    /* Compute apparent position */
    ObitPrecessUVJPrecessApp (in->myData->myDesc, source);
    RAOff  = source->RAApp;
    DecOff = source->DecApp;
    /* uvrot = rotation to north */
    dRa = (RAOff-in->curSource->RAApp) * cos(DG2RAD*source->DecApp);
    source = ObitSourceUnref(source);
    in->uvrot = -(ofloat)atan2(dRa, DecOff-in->curSource->DecApp);
    in->cuvrot = cos(in->uvrot);
    in->suvrot = sin(in->uvrot);
    in->cDec   = cos(in->curSource->DecApp*DG2RAD);
    in->sDec   = sin(in->curSource->DecApp*DG2RAD);
  } /* end new source */
    
  /* Average baseline longitude  */
  BLLong = 0.5*(in->AntList->ANlist[ia1]->AntLong + in->AntList->ANlist[ia2]->AntLong);
  
  bl[0] = in->AntList->ANlist[ia1]->AntXYZ[0] - in->AntList->ANlist[ia2]->AntXYZ[0];
  bl[1] = in->AntList->ANlist[ia1]->AntXYZ[1] - in->AntList->ANlist[ia2]->AntXYZ[1];
  bl[2] = in->AntList->ANlist[ia1]->AntXYZ[2] - in->AntList->ANlist[ia2]->AntXYZ[2];
  
  /* LST and hour angle (radians) */
  AntLst = in->AntList->GSTIAT0 + BLLong + time*in->AntList->RotRate;
  RAR    = in->curSource->RAApp*DG2RAD; /* RA in radians */
  HrAng  = AntLst - RAR;                /* Hour angle */
  
  /* Compute uvw - short baseline approximation */
  cosha  = cos (HrAng);
  sinha  = sin (HrAng);
  vw     =  bl[0]*cosha - bl[1]*sinha;
  uvw[0] =  bl[0]*sinha + bl[1]*cosha;
  uvw[1] = -vw*in->sDec + bl[2]*in->cDec;
  uvw[2] =  vw*in->cDec + bl[2]*in->sDec;
  
  /* Rotate in u-v plane to north of standard epoch */
  u = uvw[0];
  v = uvw[1];
  uvw[0] = u*in->cuvrot - v*in->suvrot;
  uvw[1] = v*in->cuvrot + u*in->suvrot;
  
} /* end ObitUVWCalcUVW */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitUVWCalcClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitUVWCalcClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitUVWCalcClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitUVWCalcClassInfoDefFn (gpointer inClass)
{
  ObitUVWCalcClassInfo *theClass = (ObitUVWCalcClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitUVWCalcClassInit;
  theClass->newObit       = (newObitFP)newObitUVWCalc;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitUVWCalcClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitUVWCalcGetClass;
  theClass->ObitCopy      = (ObitCopyFP)ObitUVWCalcCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitUVWCalcClear;
  theClass->ObitInit      = (ObitInitFP)ObitUVWCalcInit;
  theClass->ObitUVWCalcCreate = (ObitUVWCalcCreateFP)ObitUVWCalcCreate;
  theClass->ObitUVWCalcUVW    = (ObitUVWCalcUVWFP)ObitUVWCalcUVW;
} /* end ObitUVWCalcClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitUVWCalcInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitUVWCalc *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->myData   = NULL;
  in->AntList  = NULL;
  in->SouList  = NULL;
  in->mySource = NULL;
  in->curSource= NULL;
  in->antIndex = NULL;
  in->maxAnt   = -1;
  in->curSID   = -1;
  in->uvrot    = 0.0;

} /* end ObitUVWCalcInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitUVWCalc* cast to an Obit*.
 */
void ObitUVWCalcClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitUVWCalc *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  in->myData   = ObitUVUnref (in->myData);
  in->AntList  = ObitAntennaListUnref (in->AntList);
  in->SouList  = ObitSourceListUnref (in->SouList);
  in->mySource = ObitSourceUnref (in->mySource);
  if (in->antIndex) g_free(in->antIndex);
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitUVWCalcClear */

