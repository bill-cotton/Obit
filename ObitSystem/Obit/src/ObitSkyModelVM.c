/* $Id: ObitSkyModelVM.c,v 1.9 2007/09/07 12:29:52 bcotton Exp $    */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2006-2008                                          */
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

#include "ObitTableCCUtil.h"
#include "ObitSkyModelVM.h"
#include "ObitFFT.h"
#include "ObitUVUtil.h"
#include "ObitImageUtil.h"
#include "ObitPBUtil.h"
#include "ObitMem.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitSkyModelVM.c
 * ObitSkyModelVM class function definitions.
 * This is a virtual ObitSkyModel class for implementing incorporation of 
 * temporarily and/or spatially variable effects.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitSkyModelVM";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitSkyModelGetClass;

/**
 * ClassInfo structure ObitSkyModelVMClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitSkyModelVMClassInfo myClassInfo = {FALSE};

/*---------------Private function prototypes----------------*/
/** Private: FT by DFT, may be overridden in derived class */
void ObitSkyModelVMFTDFT (ObitSkyModelVM *in, olong field, ObitUV *uvdata, ObitErr *err);

/** Private: Initialize newly instantiated object. */
static void  ObitSkyModelVMInit  (gpointer in);

/** Private: Deallocate members. */
static void  ObitSkyModelVMClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitSkyModelVMClassInfoDefFn (gpointer inClass);

/** Private: Chose model type */
void  ObitSkyModelVMChose (ObitSkyModel* in, ObitUV* uvdata);

/** Private: Load point model */
gboolean ObitSkyModelVMLoadPoint (ObitSkyModel *in, ObitUV *uvdata, ObitErr *err);

/** Private: Load Components< */
gboolean ObitSkyModelVMLoadComps (ObitSkyModel *in, olong n, ObitUV *uvdata, 
				  ObitErr *err);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitSkyModelVM* newObitSkyModelVM (gchar* name)
{
  ObitSkyModelVM* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitSkyModelVMClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitSkyModelVM));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitSkyModelVMInit((gpointer)out);

 return out;
} /* end newObitSkyModelVM */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitSkyModelVMGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitSkyModelVMClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitSkyModelVMGetClass */

/**
 * Make a deep copy of an ObitSkyModelVM.
 * Since this is a virtual class it only calls parent class.
 * \param inn  The object to copy
 * \param outt An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitSkyModel* ObitSkyModelVMCopy (ObitSkyModel *inn, ObitSkyModel *outt, 
				  ObitErr *err)
{
  const ObitSkyModelClassInfo *ParentClass;
  ObitSkyModelVM *in  = (ObitSkyModelVM*)inn;
  ObitSkyModelVM *out = (ObitSkyModelVM*)outt;
  /*gchar *routine = "ObitSkyModelVMCopy";*/

  /* Copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL) && 
	    /* Don't call yourself */
	    (ParentClass!=(const ObitSkyModelClassInfo*)&myClassInfo));
  out = (ObitSkyModelVM*)ParentClass->ObitCopy (in, out, err);

   /* Copy times */
  out->endVMModelTime = in->endVMModelTime;
  out->curVMModelTime = in->curVMModelTime;

  /* Component array */
  if (in->VMComps) out->VMComps = ObitFArrayCopy(in->VMComps, out->VMComps, err);

  return (ObitSkyModel*)out;
} /* end ObitSkyModelVMCopy */

/**
 * Creates an ObitSkyModelVM 
 * This is a virtual class, calling this routine is an error
 * \param name  An optional name for the object.
 * \param mosaic ObitImageMosaic giving one or more images/CC tables
 * \return the new object.
 */
ObitSkyModelVM* ObitSkyModelVMCreate (gchar* name, ObitImageMosaic* mosaic)
{
  ObitSkyModelVM* out = NULL;
  gchar *routine = "ObitSkyModelVMCreate";

  /* should never be called - complain */ 
  g_error("%s: Virtual routine stubbed", routine);
  return out;
} /* end ObitSkyModelVMCreate */

/**
 * Initializes an ObitSkyModelVM 
 * Checks range of requested components.
 * Recurses through inheritence heirarchy
 * \param inn  SkyModel to initialize
 * \param uvdata  uv data being modeled.
 * \param err Obit error stack object.
 * \return the new object.
 */
void ObitSkyModelVMInitMod (ObitSkyModel* inn, ObitUV *uvdata, ObitErr *err)
{
  ObitSkyModelVM *in = (ObitSkyModelVM*)inn;
  ObitTable *tempTable=NULL;
  ObitTableCC *CCTable = NULL;
  ObitIOCode retCode;
  olong ver, i;
  gchar *tabType = "AIPS CC";
  gchar *routine = "ObitSkyModelVMInitMod";

  /* initialize Base class */
  ObitSkyModelInitMod(inn, uvdata, err);

  /* Any initialization this class */
  /* Reset time of current model */
  in->endVMModelTime = in->curVMModelTime = -1.0e20;

  /* Is a point model being used? */
  if ((in->pointFlux!=0.0) && (in->mosaic->numberImages <= 1)) return;

  /* Check start and end component numbers */
  for (i=0; i<in->mosaic->numberImages; i++) {
    /* Get CC table */
    ver = in->CCver[i];
    tempTable = newObitImageTable (in->mosaic->images[i], OBIT_IO_ReadOnly, 
				   tabType, &ver, err);
    if (tempTable==NULL) {
      Obit_log_error(err, OBIT_Error,"%s: Error finding %s Table %d for %s", 
		     routine, tabType, ver,  in->mosaic->images[i]->name);
      return;
    }
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    CCTable = ObitTableCCConvert(tempTable);
    tempTable = ObitTableUnref(tempTable);
    if (err->error) Obit_traceback_msg (err, routine, in->name);

    /* Open */
    retCode = ObitTableCCOpen (CCTable, OBIT_IO_ReadOnly, err);
    if ((retCode != OBIT_IO_OK) || (err->error))
      Obit_traceback_msg (err, routine, in->name);

    /* How many? */
    in->endComp[i] = MIN (in->endComp[i], CCTable->myDesc->nrow);
    if (in->endComp[i]<=0) in->endComp[i] = CCTable->myDesc->nrow;
    in->startComp[i] = MIN (in->startComp[i], CCTable->myDesc->nrow+1);
    in->startComp[i] = MAX (in->startComp[i], 1);

    /* Close */
    retCode = ObitTableCCClose (CCTable, err);
    if ((retCode != OBIT_IO_OK) || (err->error))
      Obit_traceback_msg (err, routine, in->name);
    CCTable = ObitTableUnref(CCTable);
  } /* End loop over images */

} /* end ObitSkyModelVMInitMod */

/**
 * Any shutdown operations needed for a model
 * \param inn  SkyModel to initialize
 * \param uvdata  uv data being modeled.
 * \param err Obit error stack object.
 */
void ObitSkyModelVMShutDownMod (ObitSkyModel* inn, ObitUV *uvdata, ObitErr *err)
{

  /* Call parent shutdown */
  ObitSkyModelShutDownMod(inn, uvdata, err);

} /* end ObitSkyModelVMShutDownMod */

/**
 * Initializes an ObitSkyModel for a pass through data in time order
 * Reset times for curent model.
 * \param inn  SkyModel to initialize
 * \param err Obit error stack object.
 */
void ObitSkyModelVMInitModel (ObitSkyModel* inn, ObitErr *err)
{
  ObitSkyModelVM *in = (ObitSkyModelVM*)inn;

  /* Reset time of current model */
  in->endVMModelTime = -1.0e20;
  in->curVMModelTime = -1.0e20;
} /* end ObitSkyModelVMInitModel */

/**
 * Update VM model with time or spatial modifications to model
 * This is a virtual routine that needs to be over ridden in derived classes
 * \param in      SkyModelVM 
 * \param time    current time (d)
 * \param suba    0-rel subarray number
 * \param uvdata  uv data being modeled.
 * \param err Obit error stack object.
 */
void ObitSkyModelVMUpdateModel (ObitSkyModelVM *in, ofloat time, olong suba,
				ObitUV *uvdata, ObitErr *err)
{
  gchar *routine = "ObitSkyModelVMUpdateModel";

  /* should never be called - complain */ 
  Obit_log_error(err, OBIT_Error,"%s: Virtual routine stubbed", routine);
} /* end ObitSkyModelVMUpdateModel */


/**
 * Initialize global ClassInfo Structure.
 */
void ObitSkyModelVMClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitSkyModelVMClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitSkyModelVMClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitSkyModelVMClassInfoDefFn (gpointer inClass)
{
  ObitSkyModelVMClassInfo *theClass = (ObitSkyModelVMClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitSkyModelVMClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitSkyModelVMClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitSkyModelVMGetClass;
  theClass->newObit       = (newObitFP)newObitSkyModelVM;
  theClass->ObitCopy      = NULL;  /* Virtual */
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitSkyModelVMClear;
  theClass->ObitInit      = (ObitInitFP)ObitSkyModelVMInit;
  theClass->ObitSkyModelCreate  = (ObitSkyModelCreateFP)ObitSkyModelVMCreate;
  theClass->ObitSkyModelInitMod = (ObitSkyModelInitModFP)ObitSkyModelVMInitMod;
  theClass->ObitSkyModelShutDownMod= (ObitSkyModelShutDownModFP)ObitSkyModelVMShutDownMod;
  theClass->ObitSkyModelInitModel= (ObitSkyModelInitModelFP)ObitSkyModelVMInitModel;
  theClass->ObitSkyModelLoadPoint = (ObitSkyModelLoadPointFP)ObitSkyModelVMLoadPoint;
  theClass->ObitSkyModelLoadComps = (ObitSkyModelLoadCompsFP)ObitSkyModelVMLoadComps;
  theClass->ObitSkyModelLoadImage = NULL; /* unsupported */
  theClass->ObitSkyModelFTDFT     = (ObitSkyModelFTDFTFP)ObitSkyModelVMFTDFT;
  theClass->ObitSkyModelGetInput  = (ObitSkyModelGetInputFP)ObitSkyModelVMGetInput;
  theClass->ObitSkyModelChose     = (ObitSkyModelChoseFP)ObitSkyModelVMChose;
  theClass->ObitSkyModelVMUpdateModel=
    (ObitSkyModelVMUpdateModelFP)ObitSkyModelVMUpdateModel;

} /* end ObitSkyModelVMClassDefFn */


/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitSkyModelVMInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitSkyModelVM *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->VMComps        = NULL;
  in->endVMModelTime = -1.0e20;
  in->curVMModelTime = -1.0e20;
} /* end ObitSkyModelVMInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * For some reason this wasn't build into the GType class.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitSkyModelVM* cast to an Obit*.
 */
void ObitSkyModelVMClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitSkyModelVM *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  in->VMComps = ObitFArrayUnref(in->VMComps);
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitSkyModelVMClear */

/**
 * Get input parameters from info member
 * \param inn Pointer to the ObitSkyModelVM .
 * \param err Obit error stack object.
 */
void  ObitSkyModelVMGetInput (ObitSkyModel* inn, ObitErr *err)
{
  ObitSkyModelVM *in = (ObitSkyModelVM*)inn;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  union ObitInfoListEquiv InfoReal; 
  gchar *routine = "ObitSkyModelVMGetInput";

  /* error checks */
  if (err->error) return;

  /* Call base class version */
  ObitSkyModelGetInput (inn, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Model update interval - default 1 min */
  if (in->updateInterval<=0.0) in->updateInterval = 1.0;
  InfoReal.flt = in->updateInterval; type = OBIT_float;
  ObitInfoListGetTest(in->info, "UpdateInt", &type, dim, &InfoReal);
  in->updateInterval = InfoReal.flt;
  /* Convert to days */
  in->updateInterval /=  1400.0;

} /* end ObitSkyModelVMGetInput */

/**
 * Decide which method is the most appropriate to calculate the FT of a model
 * Sets currentMode member function
 * Only DFT supported
 * \param in     Pointer to theObitSkyModel .
 * \param uvdata UV data set
 */
void  ObitSkyModelVMChose (ObitSkyModel* in, ObitUV* uvdata) 
{
  in->currentMode = OBIT_SkyModel_DFT;
  return;
} /* end ObitSkyModelVMChose */


/**
 * Load point model into in comps member.
 * Multiplies by factor member.
 * This function may be overridden in a derived class and 
 * should always be called by its function pointer.
 * Adapted from the AIPSish QNOT:VISDFT
 * Output is in member comps with a single row, the entries are
 * \li Field (0)
 * \li CC DeltaX
 * \li CC DeltaY
 * \li Amplitude (Jy)
 * \li -2*pi*x (radians)
 * \li -2*pi*y (radians)
 * \li -2*pi*z (radians)
 * \param inn  SkyModel 
 * \param uvdata UV data set to model
 * \param err Obit error stack object.
 * \return TRUE iff this image produced a valid model (i.e. had some CCs).
 */
gboolean ObitSkyModelVMLoadPoint (ObitSkyModel *inn, ObitUV *uvdata, ObitErr *err)
{
  ObitSkyModelVM *in = (ObitSkyModelVM*)inn;
  gboolean gotSome = FALSE;
  ObitIOCode retCode = OBIT_IO_SpecErr;
  olong ndim, naxis[2];
  ofloat *table, const2, ccrot, ssrot, cpa, spa, uvrot, xmaj, xmin;
  ofloat dxyzc[3], xxoff, yyoff, zzoff;
  gchar *routine = "ObitSkyModelVMLoadPoint";
  
  /* error checks */
  if (err->error) return gotSome;
  retCode = OBIT_IO_OK;

  gotSome = (in->pointFlux * in->factor!=0.0);  /* Non zero model? */
  if (!gotSome) return gotSome;

  /* Get position phase shift parameters */
  ObitUVDescShiftPosn(uvdata->myDesc, -in->pointXOff, -in->pointYOff, 
		       dxyzc, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, gotSome);
 
  /* Set field center offsets */
  uvrot  = ObitUVDescRotate(uvdata->myDesc);
  ssrot = sin (DG2RAD * uvrot);
  ccrot = cos (DG2RAD * uvrot);
  xxoff = dxyzc[0] * ccrot + dxyzc[1] * ssrot;
  yyoff = dxyzc[1] * ccrot - dxyzc[0] * ssrot;
  zzoff = dxyzc[2];
   
  in->modType = in->pointParms[3] + 0.5;  /* Model component type */

  /* (re)allocate structure */
  ndim = 2;
  naxis[0] = 7; naxis[1] = 1;
  if (in->pointParms[3]==1) naxis[0] += 3; /* Gaussian */
  if (in->pointParms[3]==2) naxis[0] += 2; /* Uniform sphere */
  if (in->comps!=NULL) in->comps = ObitFArrayRealloc(in->comps, ndim, naxis);
  else in->comps = ObitFArrayCreate("Components", ndim, naxis);

  /* Fill values */
  naxis[0] = 0; naxis[1]=0;
  table = ObitFArrayIndex(in->comps, naxis);
  table[0] = 0;
  table[1] = in->pointXOff;
  table[2] = in->pointYOff;
  table[3] = in->pointFlux * in->factor;
  table[4] = xxoff;
  table[5] = yyoff;
  table[6] = zzoff;

  /* Gaussian */
  if (in->modType==OBIT_SkyModel_GaussMod) {
    /* const2 converts FWHM(deg) to coefficients for u*u, v*v, u*v */
    const2 = DG2RAD * (G_PI / 1.17741022) * sqrt (0.5) * 2.77777778e-4;
    cpa = cos (DG2RAD * in->pointParms[2]);
    spa = sin (DG2RAD * in->pointParms[2]);
    xmaj = in->pointParms[0] * const2;
    xmin = in->pointParms[1] * const2;
    table[7] = -(((cpa * xmaj)*(cpa * xmaj)) + (spa * xmin)*(spa * xmin));
    table[8] = -(((spa * xmaj)*(spa * xmaj)) + (cpa * xmin)*(cpa * xmin));
    table[9] = -2.0 *  cpa * spa * (xmaj*xmaj - xmin*xmin);
  }

  /* Uniform sphere */
  if (in->modType==OBIT_SkyModel_USphereMod) {
    table[0] = 3.0 * in->pointFlux * in->factor;
    table[7] = in->pointParms[1]  * 0.109662271 * 2.7777778e-4;
    table[8] = 0.1;
 }
  /* One component in model */
  in->numComp = 1;

  return gotSome;
} /* end ObitSkyModelVMLoadPoint */

/**
 * Load components model into in comps member.
 * Multiplies by factor member.
 * This function may be overridden in a derived class and 
 * should always be called by its function pointer.
 * Adapted from the AIPSish QNOT:IONDFT
 * \param inn  SkyModel 
 * \param n   Image number on mosaic, if -1 load all images
 * \param uvdata UV data set to model
 * \param err Obit error stack object.
 * Output is in member comps, the entries are
 * \li Field (0-rel)
 * \li CC DeltaX
 * \li CC DeltaY
 * \li Amplitude (Jy)
 * \li -2*pi*x (radians)
 * \li -2*pi*y (radians)
 * \li -2*pi*z (radians)
 * \li Other model parameters depending on model type
 * \return TRUE iff this image produced a valid model (i.e. had some CCs).
 */
gboolean ObitSkyModelVMLoadComps (ObitSkyModel *inn, olong n, ObitUV *uvdata, 
				  ObitErr *err)
{
  ObitSkyModelVM *in = (ObitSkyModelVM*)inn;
  gboolean gotSome = FALSE;
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitTable *tempTable=NULL;
  ObitTableCC *CCTable = NULL;
  ObitTableCCRow *CCRow = NULL;
  ObitImageDesc *imDesc=NULL;
  ObitUVDesc *uvDesc=NULL;
  ObitFArray *CompArr=NULL;
  olong warray, larray;
  ofloat *array, parms[20], range[2];
  olong ver, i, j, hi, lo, count, ncomp, startComp, endComp, irow, lrec;
  olong outCCVer, ndim, naxis[2];
  ofloat *table, xxoff, yyoff, zzoff;
  ofloat konst, konst2, xyz[3], xp[3], umat[3][3], pmat[3][3];
  ofloat ccrot, ssrot, xpoff, ypoff, maprot, uvrot;
  ofloat dxyzc[3], cpa, spa, xmaj, xmin;
  gboolean doCheck=FALSE, want, do3Dmul;
  gchar *tabType = "AIPS CC";
  gchar *routine = "ObitSkyModelVMLoadComps";
  
  /* error checks */
  if (err->error) return gotSome;

  /* Don't bother if no components requested */
  if ((n>=0) && (in->startComp[n]>in->endComp[n])) return gotSome;

  /* Uv descriptor */
  uvDesc = uvdata->myDesc;

  konst = DG2RAD * 2.0 * G_PI;
  /* konst2 converts FWHM(deg) to coefficients for u*u, v*v, u*v */
  konst2 = DG2RAD * (G_PI / 1.17741022) * sqrt (0.5) * 2.77777778e-4;

  /* Loop over images counting CCs */
  count = 0;
  in->modType = OBIT_SkyModel_Unknown; /* Model type not known */
  if (in->mosaic) {lo = 0; hi = in->mosaic->numberImages-1;}
  else {lo = 0; hi = 0;}
  if (n>=0) {lo = n; hi = n;}
  for (i=lo; i<=hi; i++) {

    /* Expect anything in this table? */
    if ((in->startComp[i]>in->endComp[i]) || (in->endComp[i]<=0)) continue;

    /* Get CC table */
    ver = in->CCver[i];
    tempTable = newObitImageTable (in->mosaic->images[i],OBIT_IO_ReadOnly, 
				   tabType, &ver, err);
    if ((tempTable==NULL) || (err->error)) 
      Obit_traceback_val (err, routine, in->name, retCode);
    CCTable = ObitTableCCConvert(tempTable);
    tempTable = ObitTableUnref(tempTable);
    if (err->error) Obit_traceback_val (err, routine, in->name, retCode);

    /* Open */
    retCode = ObitTableCCOpen (CCTable, OBIT_IO_ReadOnly, err);
    if ((retCode != OBIT_IO_OK) || (err->error))
      Obit_traceback_val (err, routine, in->name, retCode);

    /* How many? */
    endComp = in->endComp[i];
    if (endComp<=0) endComp = CCTable->myDesc->nrow;
    count += MIN(CCTable->myDesc->nrow, endComp) - MAX(1, in->startComp[i]) + 1;

    /* Get model type in first with components */
    /* If only 3 col, or parmsCol 0 size then this is a point model */
    if ((CCTable->myDesc->nfield==3) || 
	(CCTable->myDesc->dim[CCTable->parmsCol]<=0))
      in->modType = OBIT_SkyModel_PointMod;
    if ((in->modType == OBIT_SkyModel_Unknown) && (in->startComp[i]<=endComp)) {
      /* Create table row */
      CCRow = newObitTableCCRow (CCTable);
      /* Read first */
      irow = in->startComp[i];
      retCode = ObitTableCCReadRow (CCTable, irow, CCRow, err);
      if ((retCode != OBIT_IO_OK) || (err->error)) 
	Obit_traceback_val (err, routine, in->name, retCode);

      /* Get model type */
      in->modType = CCRow->parms[3] + 0.5;
      /* Release table row */
      CCRow = ObitTableCCRowUnref (CCRow);
    }

    /* Do we need to check model type */
    doCheck = (CCTable->myDesc->nfield>4) && (CCTable->parmsCol>=0) && 
      (CCTable->myDesc->dim[CCTable->parmsCol][0]>=3);
    
    /* Close */
    retCode = ObitTableCCClose (CCTable, err);
    if ((retCode != OBIT_IO_OK) || (err->error))
      Obit_traceback_val (err, routine, in->name, retCode);

    /* release table  */
    CCTable = ObitTableCCUnref (CCTable);
  } /* end loop counting CCs */

  /* (re)allocate structure */
  ndim = 2;
  naxis[0] = 7; naxis[1] = count;
  if (in->modType==OBIT_SkyModel_GaussMod) naxis[0] += 3; /* Gaussian */
  if (in->modType==OBIT_SkyModel_USphereMod) naxis[0] += 2; /* Uniform sphere */
  if (in->comps!=NULL) in->comps = ObitFArrayRealloc(in->comps, ndim, naxis);
  else in->comps = ObitFArrayCreate("Components", ndim, naxis);
  lrec = naxis[0]; /* Save size of entry */
  /* Get pointer */
  naxis[0] = 0; naxis[1]=0; 
  table = ObitFArrayIndex(in->comps, naxis);

  /* Loop over images loading CCs */
  ncomp = 0;
  for (i=lo; i<=hi; i++) {

    /* Anything to do? */
    if (in->endComp[i]<in->startComp[i]) continue;

    /* Get CC table */
    outCCVer = 0;
    ver = in->CCver[i];
    startComp = in->startComp[i];
    endComp = in->endComp[i];
    range[0] = in->minDFT;  /* Range of merged fluxes for DFT */
    range[1] = 1.0e20;
    CCTable = ObitSkyModelgetPBCCTab (inn, uvdata, (olong)i, &ver, &outCCVer, 
				      &startComp, &endComp, range, err); 
    if (err->error) Obit_traceback_val (err, routine, in->name, retCode);
    
    /* Save values of highest comp - probably bad */
    if (outCCVer==0) {
      /* no translation of table */
      /*??? in->endComp[i] = endComp; */
    } else {
      /* Translated table with only selected values */
      /*??? in->endComp[i] = in->startComp[i] + endComp-startComp; */
    }
    
    /* Anything to do? */
    if (endComp<startComp) {
      /* No - free up tables */
      /* if outCCver>0 then the CCtable is temporary - Zap */
      if (outCCVer>0) {
	CCTable = ObitTableCCUnref (CCTable);
	ObitImageZapTable(in->mosaic->images[i], tabType, outCCVer, err);
      /* else simply release table  */
      } else CCTable = ObitTableCCUnref (CCTable);
      if (err->error) Obit_traceback_val (err, routine, in->name, retCode);
      continue;
    }

    /* Field specific stuff */
    imDesc = in->mosaic->images[i]->myDesc; /* Image descriptor */
    /*  Set field center offsets. */
    maprot = ObitImageDescRotate(imDesc);
    uvrot  = ObitUVDescRotate(uvDesc);
    ssrot = sin (DG2RAD * (uvrot - maprot));
    ccrot = cos (DG2RAD * (uvrot - maprot));

    /* Get position phase shift parameters */
    ObitUVDescShiftPhase(uvDesc, imDesc, dxyzc, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, retCode);
    
    /*    Get reference pixel offsets from tangent point */
    /* xpoff = (imDesc->crpix[imDesc->jlocr] - 
       (imDesc->inaxes[imDesc->jlocr]/2)) * 
       imDesc->cdelt[imDesc->jlocr];
       ypoff = (imDesc->crpix[imDesc->jlocd] - 
       (imDesc->inaxes[imDesc->jlocd]/2) - 1) * 
       imDesc->cdelt[imDesc->jlocd];*/
    /* These should always be zero for 3D imaging? */
    xpoff = 0.0;
    ypoff = 0.0;
    
    /* Set field center offsets */
    xxoff = dxyzc[0] * ccrot + dxyzc[1] * ssrot;
    yyoff = dxyzc[1] * ccrot - dxyzc[0] * ssrot;
    zzoff = dxyzc[2];

    /* 3D rotation matrix if needed */
    if (in->do3D) {
      do3Dmul = ObitUVDescShift3DMatrix (uvDesc, imDesc, umat, pmat);
    } else {do3Dmul = FALSE;}
    
    /* Convert table to merged array */
    CompArr = ObitTableCCUtilMergeSel (CCTable, startComp, endComp, parms, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, retCode);
    /* entries 0=flux, 1= deltaX 2=deltaY per merged CC, other parameters in parms */
    naxis[0] = 0; naxis[1]=0; 
    array = ObitFArrayIndex(CompArr, naxis);
    warray = CompArr->naxis[0];
    larray = CompArr->naxis[1];
   
    /* loop over CCs */
    for (j=0; j<larray; j++) {
 
     /* Only down to first negative? */
      if (in->noNeg && (array[0]<=0.0)) break;

      /* Do we want this one? Only accept components of in->modType */
      want = !doCheck;
      want = want && (fabs(array[0])>0.0);
      want = want && (array[0]>in->minFlux);
      want = want && (ncomp<count);  /* don't overflow */
      if (want) {

	/* Field number */
	table[0] = i;
	table[1] = array[1];
	table[2] = array[2];

	/* Point */
	table[3] = array[0] * in->factor;
	xp[0] = (array[1] + xpoff) * konst;
	xp[1] = (array[2] + ypoff) * konst;
	xp[2] = 0.0;
	if (do3Dmul) {
	  xyz[0] = xp[0]*umat[0][0] + xp[1]*umat[1][0];
	  xyz[1] = xp[0]*umat[0][1] + xp[1]*umat[1][1];
	  xyz[2] = xp[0]*umat[0][2] + xp[1]*umat[1][2];
	  /* PRJMUL (2, XP, UMAT, XYZ); */
	} else {  /* no 3D */
	  xyz[0] = ccrot * xp[0] + ssrot * xp[1];
	  xyz[1] = ccrot * xp[1] - ssrot * xp[0];
	  xyz[2] = 0.0;
 	}
	table[4] = xyz[0] + xxoff;
	table[5] = xyz[1] + yyoff;
	table[6] = xyz[2] + zzoff;

	/* Gaussian */
	if (in->modType==OBIT_SkyModel_GaussMod) {
	  cpa = cos (DG2RAD * parms[2]);
	  spa = sin (DG2RAD * parms[2]);
	  xmaj = parms[0] * konst2;
	  xmin = parms[1] * konst2;
	  table[7] = -(((cpa * xmaj)*(cpa * xmaj)) + (spa * xmin)*(spa * xmin));
	  table[8] = -(((spa * xmaj)*(spa * xmaj)) + (cpa * xmin)*(cpa * xmin));
	  table[9] = -2.0 *  cpa * spa * (xmaj*xmaj - xmin*xmin);
	}

	/* Uniform sphere */
	if (in->modType==OBIT_SkyModel_USphereMod) {
	  table[3] = 3.0 * array[0] * in->factor;
	  table[7] = CCRow->parms[1]  * 0.109662271 * 2.7777778e-4;
	  table[8] = 0.1;
	}
	    
	/* Update */
	table += lrec;
	array += warray;
	ncomp++;
      } /* End only desired */
    } /* end loop over components */

    /* Delete merged CC array */
    CompArr = ObitFArrayUnref(CompArr);

    /* if outCCver>0 then the CCtable is temporary - Zap */
    if (outCCVer>0) {
      CCTable = ObitTableCCUnref (CCTable);
      ObitImageZapTable(in->mosaic->images[i], tabType, outCCVer, err);
    /* else simply release table  */
    } else CCTable = ObitTableCCUnref (CCTable);
    if (err->error) Obit_traceback_val (err, routine, in->name, retCode);

    /* Release table row */
    CCRow = ObitTableCCRowUnref (CCRow);
    
  } /* end loop loading CCs */

  /* Zero any extra entries in table. */
  for (i=ncomp; i<count; i++) {
    /* Zero entry */
    table[0] = -10.0;
    table[1] = 0.0;
    table[2] = 0.0;
    table[3] = 0.0;
    table[4] = 0.0;
    table[5] = 0.0;
    table[6] = 0.0;
    table += lrec;  /* Update pointer */
  } /* end loop zeroing extra components */

  in->numComp = ncomp;

  /* Find anything */
  gotSome = ncomp>0;

  return gotSome;
} /* end ObitSkyModelVMLoadComps */

/**
 * Do Fourier transform using a DFT for a buffer of data.
 * Version for time/spatial dependent effects.
 * If doDivide member is true then FT of model is divided into the data,
 * If doReplace member is true then FT of model replaces the data,
 * else, it is subtracted.
 * If doFlip member is true the Fourier transform is multiplied by sqrt(-1)
 * (for Stokes RL and LR)
 * After the AIPSish QXXPTS, QPTDIV and friends
 * This function may be overridden in a derived class and 
 * should always be called by its function pointer.
 * \param in     SkyModelVM with model components loaded (ObitSkyModelLoad)
 * \param field  Field number being processed (-1 => all)
 * \param uvdata UV data set to model and subtract from current buffer
 * \param err Obit error stack object.
 * \return return code, OBIT_IO_OK=> OK
 */
void ObitSkyModelVMFTDFT (ObitSkyModelVM *in, olong field, ObitUV *uvdata, ObitErr *err)
{
  olong iVis, iIF, iChannel, iStoke, iComp, lcomp, ncomp;
  olong lrec, nrparm, naxis[2];
  olong startPoln, numberPoln, jincs, startChannel, numberChannel;
  olong jincf, startIF, numberIF, jincif, kincf, kincif;
  olong offset, offsetChannel, offsetIF;
  olong ilocu, ilocv, ilocw, iloct, suba, itemp;
  ofloat *visData, *ccData, *data, *fscale;
  ofloat modReal, modImag;
  ofloat amp, arg, freq2, freqFact, wt=0.0, temp;
  odouble phase, tx, ty, tz, sumReal, sumImag, *freqArr;
  const ObitSkyModelVMClassInfo 
    *myClass=(const ObitSkyModelVMClassInfo*)in->ClassInfo;
  gchar *routine = "ObitSkyModelVMFTDFT";

  /* error checks - assume most done at higher level */
  if (err->error) return;

  /* Create buffer for VM model components if needbe */
  if (in->VMComps==NULL) 
    in->VMComps = ObitFArrayCopy (in->comps, in->VMComps, err);
  /* Resize if necessary */
  if (!ObitFArrayIsCompatable(in->comps, in->VMComps))
    ObitFArrayClone (in->comps, in->VMComps, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  
  /* Visibility pointers */
  ilocu =  uvdata->myDesc->ilocu;
  ilocv =  uvdata->myDesc->ilocv;
  ilocw =  uvdata->myDesc->ilocw;
  iloct =  uvdata->myDesc->iloct;

  /* Set channel, IF and Stokes ranges (to 0-rel)*/
  startIF  = in->startIFPB-1;
  numberIF = MAX (1, in->numberIFPB);
  jincif   = uvdata->myDesc->incif;
  startChannel  = in->startChannelPB-1;
  numberChannel = MAX (1, in->numberChannelPB);
  jincf         = uvdata->myDesc->incf;
  startPoln  = in->startPoln-1;
  numberPoln = in->numberPoln;
  jincs      = uvdata->myDesc->incs;  /* increment in real array */
  /* Increments in frequency tables */
  if (uvdata->myDesc->jlocf<uvdata->myDesc->jlocif) { /* freq before IF */
    kincf = 1;
    kincif = uvdata->myDesc->inaxes[uvdata->myDesc->jlocf];
  } else { /* IF beforefreq  */
    kincif = 1;
    kincf = uvdata->myDesc->inaxes[uvdata->myDesc->jlocif];
  }

  /* Get pointer for components */
  naxis[0] = 0; naxis[1] = 0; 
  data = ObitFArrayIndex(in->VMComps, naxis);
  lcomp = in->VMComps->naxis[0];  /* Length of row in comp table */
  ncomp = in->numComp;            /* Number of components */

  /* Get pointer for frequency correction tables */
  fscale = uvdata->myDesc->fscale;
  freqArr = uvdata->myDesc->freqArr;

  /* Loop over vis in buffer */
  visData = uvdata->buffer;         /* Buffer pointer */
  lrec = uvdata->myDesc->lrec;      /* Length of record */
  nrparm = uvdata->myDesc->nrparm;  /* Words of "random parameters" */
  for (iVis=0; iVis<uvdata->myDesc->numVisBuff; iVis++) {

    /* Is current model still valid? */
    if (visData[iloct] > in->endVMModelTime) {
      /* Subarray 0-rel */
      itemp = (olong)visData[uvdata->myDesc->ilocb];
      suba = 100.0 * (visData[uvdata->myDesc->ilocb]-itemp) + 0.5; 
      /* Update */
      myClass->ObitSkyModelVMUpdateModel (in, visData[iloct], suba, uvdata, err);
    }
    if (err->error) Obit_traceback_msg (err, routine, in->name);

    /* Loop over IFs */
    for (iIF=startIF; iIF<startIF+numberIF; iIF++) {
      offsetIF = nrparm + iIF*jincif; 
      for (iChannel=startChannel; iChannel<startChannel+numberChannel; iChannel++) {
	offsetChannel = offsetIF + iChannel*jincf; 
	freqFact = fscale[iIF*kincif + iChannel*kincf];  /* Frequency scaling factor */
	freq2 = freqArr[iIF*kincif + iChannel*kincf];    /* Frequency squared */
	freq2 *= freq2;

	/* Sum over components */
	/* Table values 0=Amp, 1=-2*pi*x, 2=-2*pi*y, 3=-2*pi*z */
	sumReal = sumImag = 0.0;
	ccData = data;
	
	/* Sum by model type */
	switch (in->modType) {
	case OBIT_SkyModel_PointMod:     /* Point */
	  /* From the AIPSish QXXPTS.FOR  */
	  for (iComp=0; iComp<ncomp; iComp++) {
	    tx = ccData[1]*(odouble)visData[ilocu];
	    ty = ccData[2]*(odouble)visData[ilocv];
	    tz = ccData[3]*(odouble)visData[ilocw];
	    phase = freqFact * (tx + ty + tz);
	    sumReal += ccData[0]*cos(phase);
	    sumImag += ccData[0]*sin(phase);
	    /* DEBUG
	       if (iVis==0) {
	       fprintf (stderr,"%s: comp %d real %f imag %f phase %f sum %f %f \n",
	       routine, iComp, ccData[0]*cos(phase), ccData[0]*sin(phase),
	       57.296*phase, sumReal, sumImag); 
	       } */
	    ccData += lcomp;  /* update pointer */
	  } /* end loop over components */
	  break;
	case OBIT_SkyModel_GaussMod:     /* Gaussian on sky */
	  /* From the AIPSish QGASUB.FOR  */
	  for (iComp=0; iComp<ncomp; iComp++) {
	    arg = freq2 * (ccData[4]*visData[ilocu]*visData[ilocu] +
			   ccData[5]*visData[ilocv]*visData[ilocv] +
			   ccData[6]*visData[ilocu]*visData[ilocv]);
	    if (arg>1.0e-5) amp = ccData[0] * exp (arg);
	    else amp = ccData[0];
	    tx = ccData[1]*(odouble)visData[ilocu];
	    ty = ccData[2]*(odouble)visData[ilocv];
	    tz = ccData[3]*(odouble)visData[ilocw];
	    phase = freqFact * (tx + ty + tz);
	    sumReal += amp*cos(phase);
	    sumImag += amp*sin(phase);
	    ccData += lcomp;  /* update pointer */
	  }  /* end loop over components */
	  break;
	case OBIT_SkyModel_USphereMod:    /* Uniform sphere */
	  /* From the AIPSish QSPSUB.FOR  */
	  for (iComp=0; iComp<ncomp; iComp++) {
	    arg = freqFact * sqrt(visData[ilocu]*visData[ilocu] +
				  visData[ilocv]*visData[ilocv]);
	    arg = MAX (arg, 0.1);
	    amp = ccData[0] * ((sin(arg)/(arg*arg*arg)) - cos(arg)/(arg*arg));
	    tx = ccData[1]*(odouble)visData[ilocu];
	    ty = ccData[2]*(odouble)visData[ilocv];
	    tz = ccData[3]*(odouble)visData[ilocw];
	    phase = freqFact * (tx + ty + tz);
	    sumReal += amp*cos(phase);
	    sumImag += amp*sin(phase);
	    ccData += lcomp;  /* update pointer */
	  }  /* end loop over components */
	  break;
	default:
	  Obit_log_error(err, OBIT_Error,"%s Unknown Comp model type %d in %s",
			 routine, in->modType, in->name);
	  return;
	}; /* end switch by model type */

	/* Need to multiply model by sqrt(-1)? */
	if (in->doFlip) {
	  modReal = -(ofloat)sumImag;
	  modImag =  (ofloat)sumReal;
	} else {
	  modReal =  (ofloat)sumReal;
	  modImag =  (ofloat)sumImag;
	}
	
	/* Dividing? */
	if (in->doDivide) {
	  /* Divide model - also correct weight */
	  wt = modReal * modReal + modImag * modImag;
	  modReal /= wt;
	  modImag /= wt;
	  wt = sqrt (wt);
	}

	/* Stokes Loop */
	for (iStoke=startPoln; iStoke<startPoln+numberPoln; iStoke++) {
	  offset = offsetChannel + iStoke*jincs; /* Visibility offset */

	  /* Ignore blanked data */
	  if ((visData[offset+2]<=0.0) && !in->doReplace) continue;
	  
 	  /* Apply model to data */
	  if (in->doDivide) {
	    temp = modReal * visData[offset] + modImag * visData[offset+1];
	    visData[offset+1] = modReal * visData[offset+1] - modImag * visData[offset];
	    visData[offset]   = temp;
	    visData[offset+2] *= wt;  /* correct weight */
	  } else if (in->doReplace) {  /* replace data with model */
	    visData[offset]   = modReal;
	    visData[offset+1] = modImag;
	  } else {
	    /* Subtract model */
	    visData[offset]   -= modReal;
	    visData[offset+1] -= modImag;
	  }

	  /* Factor for next Stokes */
	  modReal *= in->stokFactor;
	  modImag *= in->stokFactor;
	  
	  offset += jincs;
	} /* end loop over Stokes */
	  offsetChannel += jincf;
      } /* end loop over Channel */
 	  offsetIF += jincif;
   } /* end loop over IF */

    visData += lrec; /* Update vis pointer */
  } /* end loop over visibilities */

} /* ObitSkyModelVMFTDFT */

