/* $Id$        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2012                                               */
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

#include "ObitSourceEphemerus.h"
#include "ObitPrecess.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitSourceEphemerus.c
 * ObitSourceEphemerus class function definitions.
 * This class is derived from the Obit base class.
 * The ObitSourceEphemerus has position information of moving sources.
*/

/** name of the class defined in this file */
static gchar *myClassName = "ObitSourceEphemerus";

/** Function to obtain parent ClassInfo */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo structure ObitSourceEphemerusClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitSourceEphemerusClassInfo myClassInfo = {FALSE};

/*--------------- File Global Variables  ----------------*/


/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitSourceEphemerusInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitSourceEphemerusClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitSourceEphemerusClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitSourceEphemerus* newObitSourceEphemerus (gchar* name)
{
  ObitSourceEphemerus* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitSourceEphemerusClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitSourceEphemerus));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitSourceEphemerusInit((gpointer)out);

 return out;
} /* end newObitSourceEphemerus */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitSourceEphemerusGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitSourceEphemerusClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitSourceEphemerusGetClass */

/**
 * Make a deep copy of an ObitSourceEphemerus.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitSourceEphemerus* ObitSourceEphemerusCopy  (ObitSourceEphemerus *in, 
					       ObitSourceEphemerus *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  gboolean oldExist;
  gchar *outName;
  olong i, j, n;

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
    out = newObitSourceEphemerus(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /* Cleanup any old */
  if (out->nsrc>0) {
    for (i=0; i<out->nsrc; i++) {
      if (out->RADeriv)
	if (out->RADeriv[i])   {g_free(out->RADeriv[i]);   out->RADeriv[i]=NULL;}
      if (out->DecDeriv)
	if (out->DecDeriv[i])  {g_free(out->DecDeriv[i]);  out->DecDeriv[i]=NULL;}
      if (out->DistDeriv)
	if (out->DistDeriv[i]) {g_free(out->DistDeriv[i]); out->DistDeriv[i]=NULL;}
    }
    if (out->SID)          {g_free(out->SID);          out->SID=NULL;}
    if (out->refTime)      {g_free(out->refTime);      out->refTime=NULL;}
    if (out->RARef)        {g_free(out->RARef);        out->RARef=NULL;}
    if (out->numRADeriv)   {g_free(out->numRADeriv);   out->numRADeriv=NULL;}
    if (out->RADeriv)      {g_free(out->RADeriv);      out->RADeriv=NULL;}
    if (out->DecRef)       {g_free(out->DecRef);       out->DecRef=NULL;}
    if (out->numDecDeriv)  {g_free(out->numDecDeriv);  out->numDecDeriv=NULL;}
    if (out->DecDeriv)     {g_free(out->DecDeriv);     out->DecDeriv=NULL;}
    if (out->distRef)      {g_free(out->distRef);      out->distRef=NULL;}
    if (out->numDistDeriv) {g_free(out->numDistDeriv); out->numDistDeriv=NULL;}
    if (out->DistDeriv)    {g_free(out->DistDeriv);    out->DistDeriv=NULL;}
    out->nsrc = 0;
  } /* end cleanup old */

  /*  copy this class -  Allocate arrays */
  out->SID          = g_malloc0(in->nsrc*sizeof(olong));
  out->refTime      = g_malloc0(in->nsrc*sizeof(odouble));
  out->RARef        = g_malloc0(in->nsrc*sizeof(odouble));
  out->numRADeriv   = g_malloc0(in->nsrc*sizeof(olong));
  out->RADeriv      = g_malloc0(in->nsrc*sizeof(odouble*));
  out->DecRef       = g_malloc0(in->nsrc*sizeof(odouble));
  out->numDecDeriv  = g_malloc0(in->nsrc*sizeof(olong));
  out->DecDeriv     = g_malloc0(in->nsrc*sizeof(odouble*));
  out->distRef      = g_malloc0(in->nsrc*sizeof(odouble));
  out->numDistDeriv = g_malloc0(in->nsrc*sizeof(olong));
  out->DistDeriv    = g_malloc0(in->nsrc*sizeof(odouble*));
  out->nsrc         = 0;

  /* Loop over array */
  for (i=0; i<in->nsrc; i++) {
      n = in->numRADeriv[in->nsrc];
      out->SID[in->nsrc] = in->SID[in->nsrc];
      out->refTime[in->nsrc]      = 
      out->RARef[in->nsrc]        = 
      out->numRADeriv[in->nsrc]   = n;
      out->RADeriv[in->nsrc]      = g_malloc0(n*sizeof(odouble));
      for (j=1; j<n; j++) out->RADeriv[in->nsrc][j] = in->RADeriv[in->nsrc][j];
      out->DecRef[in->nsrc]       = in->DecRef[in->nsrc];
      out->numDecDeriv[in->nsrc]  = n;
      out->DecDeriv[in->nsrc]     = g_malloc0(n*sizeof(odouble));
      for (j=1; j<n; j++) out->DecDeriv[in->nsrc][j] = in->DecDeriv[in->nsrc][j];
      out->distRef[in->nsrc]      = 0.0;
      out->numDistDeriv[in->nsrc] = 0;
      out->DistDeriv[in->nsrc]    = NULL;
      in->nsrc++;
    }
  
  return out;
} /* end ObitSourceEphemerusCopy */

/**
 * Make a copy of a object but do not copy the actual data
 * This is useful to create an SourceEphemerus similar to the input one.
 * \param in  The object to copy
 * \param out An existing object pointer for output, must be defined.
 * \param err Obit error stack object.
 */
void ObitSourceEphemerusClone  (ObitSourceEphemerus *in, 
				ObitSourceEphemerus *out, ObitErr *err)
{
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

  /*  copy this class */
  ObitSourceEphemerusCopy (in, out, err);

} /* end ObitSourceEphemerusClone */

/**
 * Creates an ObitSourceEphemerus 
 * \param name  An optional name for the object.
 * \return the new object.
 */
ObitSourceEphemerus* ObitSourceEphemerusCreate (gchar* name)
{
  ObitSourceEphemerus* out;

  /* Create basic structure */
  out = newObitSourceEphemerus (name);

  /* Source structure for precession */
  out->source = newObitSource("EPhem");

  return out;
} /* end ObitSourceEphemerusCreate */

/**
 * Initialize an ObitSourceEphemerusSetup given an ObitASDM
 * \param in      The object to setup
 * \param SDM     SDM with info
 * \param updtime How often to update position (day)
 * \param uvDesc  UV data descriptor
 * \param err     Obit error stack object.
 */
/** Public: Setup given an ObitASDM. */
void ObitSourceEphemerusSetup (ObitSourceEphemerus *in, ObitSDMData *SDM,
			       odouble updtime, ObitUVDesc *uvDesc,
			       ObitErr *err)
{
  olong i, j, n, count;
  ASDMSourceArray*   SourceArray=NULL;
  ASDMFieldTable     *Tab;
  gchar *routine = "ObitSourceEphemerusSetup";

  /* Cleanup any old */
  if (in->nsrc>0) {
    for (i=0; i< in->nsrc; i++) {
      if (in->RADeriv)
	if (in->RADeriv[i])   {g_free(in->RADeriv[i]);   in->RADeriv[i]=NULL;}
      if (in->DecDeriv)
	if (in->DecDeriv[i])  {g_free(in->DecDeriv[i]);  in->DecDeriv[i]=NULL;}
      if (in->DistDeriv)
	if (in->DistDeriv[i]) {g_free(in->DistDeriv[i]); in->DistDeriv[i]=NULL;}
    }
    if (in->SID)          {g_free(in->SID);          in->SID=NULL;}
    if (in->refTime)      {g_free(in->refTime);      in->refTime=NULL;}
    if (in->RARef)        {g_free(in->RARef);        in->RARef=NULL;}
    if (in->numRADeriv)   {g_free(in->numRADeriv);   in->numRADeriv=NULL;}
    if (in->RADeriv)      {g_free(in->RADeriv);      in->RADeriv=NULL;}
    if (in->DecRef)       {g_free(in->DecRef);       in->DecRef=NULL;}
    if (in->numDecDeriv)  {g_free(in->numDecDeriv);  in->numDecDeriv=NULL;}
    if (in->DecDeriv)     {g_free(in->DecDeriv);     in->DecDeriv=NULL;}
    if (in->distRef)      {g_free(in->distRef);      in->distRef=NULL;}
    if (in->numDistDeriv) {g_free(in->numDistDeriv); in->numDistDeriv=NULL;}
    if (in->DistDeriv)    {g_free(in->DistDeriv);    in->DistDeriv=NULL;}
    in->nsrc = 0;
  } /* end cleanup old */

  /* Check if any sources have derivatives on position */
  Tab   = SDM->FieldTab;
  count = 0;
  for (i=0; i<Tab->nrows; i++) {
    if (Tab->rows[i]->numPoly>1) count++;
  }

  /* If none then it's easy */
  if (count<=0) return;

  /* Extract info from ASDM for source number */
  SourceArray = ObitSDMDataGetSourceArray(SDM);
  Obit_return_if_fail((SourceArray), err,
		      "%s: Could not extract Source info from ASDM", 
		      routine);

  in->uvDesc    = ObitUVDescRef(uvDesc);    /* Save data descriptor */
  in->updtime   = updtime;   /* Save update interval */
  in->lastTime  = -1.0e-10;  /* Force update */
  in->validTime = -1.0e-10;

  /* Allocate arrays */
  in->SID          = g_malloc0(count*sizeof(olong));
  in->refTime      = g_malloc0(count*sizeof(odouble));
  in->RARef        = g_malloc0(count*sizeof(odouble));
  in->numRADeriv   = g_malloc0(count*sizeof(olong));
  in->RADeriv      = g_malloc0(count*sizeof(odouble*));
  in->DecRef       = g_malloc0(count*sizeof(odouble));
  in->numDecDeriv  = g_malloc0(count*sizeof(olong));
  in->DecDeriv     = g_malloc0(count*sizeof(odouble*));
  in->distRef      = g_malloc0(count*sizeof(odouble));
  in->numDistDeriv = g_malloc0(count*sizeof(olong));
  in->DistDeriv    = g_malloc0(count*sizeof(odouble*));
  in->nsrc         = 0;

  /* Loop over array */
  count = 0;
  for (i=0; i<Tab->nrows; i++) {
    if (Tab->rows[i]->numPoly>1) {
      n = Tab->rows[i]->numPoly;
      /* Setup/fill */
      for (j=0; j<SourceArray->nsou; j++) {
	if (SourceArray->sou[j]->sourceId==Tab->rows[i]->sourceId)
	  in->SID[count] = SourceArray->sou[j]->sourceNo;
      }
      in->refTime[count]      = Tab->rows[i]->time - SDM->refJD;
      in->RARef[count]        = Tab->rows[i]->referenceDir[0];
      in->numRADeriv[count]   = n;
      in->RADeriv[count]      = g_malloc0(n*sizeof(odouble));
      for (j=1; j<n; j++) in->RADeriv[count][j] = Tab->rows[i]->referenceDir[j*2];
      in->DecRef[count]       = Tab->rows[i]->referenceDir[1];
      in->numDecDeriv[count]  = n;
      in->DecDeriv[count]     = g_malloc0(n*sizeof(odouble));
      for (j=1; j<n; j++) in->DecDeriv[count][j] = Tab->rows[i]->referenceDir[1+j*2];
      in->distRef[count]      = 0.0;
      in->numDistDeriv[count] = 0;
      in->DistDeriv[count]    = NULL;
      count++;
    }
  }
  in->nsrc = count;

  /* Cleanup */
  SourceArray = ObitSDMDataKillSourceArray(SourceArray);
} /* end ObitSourceEphemerusSetup  */

/**
 * Check if a source ID is present and if so return position.
 * Positions in ASDM are mean epoch J2000 and are precessed to 
 * apparent of date.
 * \param in     The object to check
 * \param srcID  Desired source ID
 * \param time   Time (days)
 * \param RA     [out] Apparent RA of date (rad)
 * \param Dec    [out] Apparent Dec of date (rad)
 * \param dist   [out] Distance (m), 0 if not known
 * \param uvroe  [out] rotation angle (rad) for uv to align v with J2000 north.
 * \return TRUE if source found, else FALSE
 */
gboolean 
ObitSourceEphemerusCheckSource(ObitSourceEphemerus *in, olong srcID,
			       ofloat time, odouble *RA, odouble *Dec, 
			       odouble *dist, ofloat *uvrot)
{
  gboolean found=FALSE;
  olong isrc, i;
  odouble dtime, sum, fact, ra, dec, RAOff, DecOff, dRa, RAApp, DecApp;

  /* Initialize output */
  *RA    = 0.0;
  *Dec   = 0.0;
  *dist  = 0.0;
  *uvrot = 0.0;  

  /* Anything? */
  if (in->nsrc<=0) return found;

  /* Loop over list */
  for (isrc=0; isrc<in->nsrc; isrc++) {
    if (srcID==in->SID[isrc]) {found=TRUE; break;}
  }
  if (!found) return found;

  /* Must have it - is it a new source? */
  if (srcID!=in->lastSrc) {
    in->validTime = -1.0e10;  /* force update */
    in->lastSrc   = srcID;
  }

  /* Current position still valid? */
  if (time<in->validTime) {
    *RA    = in->lastRA;
    *Dec   = in->lastDec;
    *dist  = in->lastDist;
    *uvrot = in->lastUVrot;
    return found;
 } /* end still valid */

  /* Need to update */
  in->lastTime = time;
  dtime = in->lastTime - in->refTime[isrc];
  /* RA */
  fact = dtime;
  sum  = in->RARef[isrc];
  for (i=0; i<in->numRADeriv[isrc]; i++) {
    sum += fact*in->RADeriv[isrc][i];
    fact *= dtime;
  }
  ra = sum*RAD2DG;  /* Convert RA to deg */
 
  /* Dec */
  fact = dtime;
  sum  = in->DecRef[isrc];
  for (i=0; i<in->numDecDeriv[isrc]; i++) {
    sum += fact*in->DecDeriv[isrc][i];
    fact *= dtime;
  }
  dec = sum*RAD2DG; /* Convert Dec to deg */

  /* Distance */
  fact = dtime;
  sum  = in->distRef[isrc];
  for (i=0; i<in->numDistDeriv[isrc]; i++) {
    sum += fact*in->DistDeriv[isrc][i];
    fact *= dtime;
  }
  in->lastDist = sum;

  /* Precess to apparent */
  in->source->RAMean  = ra;
  in->source->DecMean = dec;
  /* Compute apparent position */
  ObitPrecessUVJPrecessApp (in->uvDesc, in->source);
  RAApp  = in->source->RAApp;
  DecApp = in->source->DecApp;
  in->lastRA  = DG2RAD * RAApp;
  in->lastDec = DG2RAD * DecApp;

  /* Get rotation angle to align v with north at standard epoch - 
     precess posn. 10" north */
  in->source->RAMean  = ra;
  in->source->DecMean = dec + 1.0/3600.0;
  /* Compute apparent position */
  ObitPrecessUVJPrecessApp (in->uvDesc, in->source);
  RAOff  = in->source->RAApp;
  DecOff = in->source->DecApp;
  /* uvrot global = rotation to north */
  dRa = (RAOff-RAApp) * cos (DG2RAD*DecApp);
  in->lastUVrot = -(ofloat)atan2(dRa, DecOff-DecApp);

  /* Save output */
  *RA    = in->lastRA;
  *Dec   = in->lastDec;
  *uvrot = in->lastUVrot;
  *dist  = in->lastDist;

  /* Update validity time  */
  in->validTime = in->lastTime + in->updtime;

  return found;
} /* end ObitSourceEphemerusCheckSource  */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitSourceEphemerusClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitSourceEphemerusClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitSourceEphemerusClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitSourceEphemerusClassInfoDefFn (gpointer inClass)
{
  ObitSourceEphemerusClassInfo *theClass = (ObitSourceEphemerusClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitSourceEphemerusClassInit;
  theClass->newObit       = (newObitFP)newObitSourceEphemerus;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitSourceEphemerusClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitSourceEphemerusGetClass;
  theClass->ObitCopy      = (ObitCopyFP)ObitSourceEphemerusCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitSourceEphemerusClear;
  theClass->ObitInit      = (ObitInitFP)ObitSourceEphemerusInit;
  theClass->ObitSourceEphemerusCreate = (ObitSourceEphemerusCreateFP)ObitSourceEphemerusCreate;

} /* end ObitSourceEphemerusClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitSourceEphemerusInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitSourceEphemerus *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->nsrc         = 0;
  in->SID          = NULL;
  in->refTime      = NULL;
  in->RARef        = NULL;
  in->numRADeriv   = NULL;
  in->RADeriv      = NULL;
  in->DecRef       = NULL;
  in->numDecDeriv  = NULL;
  in->DecDeriv     = NULL;
  in->distRef      = NULL;
  in->numDistDeriv = NULL;
  in->DistDeriv    = NULL;
  in->source       = NULL;
  in->lastSrc      = -999;
  in->lastTime     = -1.0e-10;  /* Force update */
  in->validTime    = -1.0e-10;

} /* end ObitSourceEphemerusInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitSourceEphemerus* cast to an Obit*.
 */
void ObitSourceEphemerusClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  olong i;
  ObitSourceEphemerus *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
   if (in->nsrc>0) {
    for (i=0; i< in->nsrc; i++) {
      if (in->RADeriv)
	if (in->RADeriv[i])   {g_free(in->RADeriv[i]);   in->RADeriv[i]=NULL;}
      if (in->DecDeriv)
	if (in->DecDeriv[i])  {g_free(in->DecDeriv[i]);  in->DecDeriv[i]=NULL;}
      if (in->DistDeriv)
	if (in->DistDeriv[i]) {g_free(in->DistDeriv[i]); in->DistDeriv[i]=NULL;}
    }
    if (in->SID)          {g_free(in->SID);          in->SID=NULL;}
    if (in->refTime)      {g_free(in->refTime);      in->refTime=NULL;}
    if (in->RARef)        {g_free(in->RARef);        in->RARef=NULL;}
    if (in->numRADeriv)   {g_free(in->numRADeriv);   in->numRADeriv=NULL;}
    if (in->RADeriv)      {g_free(in->RADeriv);      in->RADeriv=NULL;}
    if (in->DecRef)       {g_free(in->DecRef);       in->DecRef=NULL;}
    if (in->numDecDeriv)  {g_free(in->numDecDeriv);  in->numDecDeriv=NULL;}
    if (in->DecDeriv)     {g_free(in->DecDeriv);     in->DecDeriv=NULL;}
    if (in->distRef)      {g_free(in->distRef);      in->distRef=NULL;}
    if (in->numDistDeriv) {g_free(in->numDistDeriv); in->numDistDeriv=NULL;}
    if (in->DistDeriv)    {g_free(in->DistDeriv);    in->DistDeriv=NULL;}
    in->nsrc = 0;
  } /* end cleanup old */

   in->source = ObitSourceUnref(in->source);
   in->uvDesc = ObitUVDescUnref(in->uvDesc);

  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitSourceEphemerusClear */

