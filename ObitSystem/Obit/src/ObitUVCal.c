/* $Id$       */
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
/*;  Correspondence this software should be addressed as follows:     */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/

#include "ObitUVCal.h"
#include "ObitUVDesc.h"
#include "ObitUVSel.h"
#include "ObitUVCalCalibrate.h"
#include "ObitUVCalSelect.h"
#include "ObitUVCalFlag.h"
#include "ObitUVCalBandpass.h"
#include "ObitUVCalBaseline.h"
#include "ObitUVCalPolarization.h"
#include "ObitTableANUtil.h"
#include "ObitTableSUUtil.h"
#include "ObitPrecess.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitUVCal.c
 * ObitUVCal class function definitions.
 * This class is derived from the Obit base class.
 */

/*--------------- File Global Variables  ----------------*/
/** name of the class defined in this file */
static gchar *myClassName = "ObitUVCal";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo structure ObitUVCalClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitUVCalClassInfo myClassInfo = {FALSE};

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitUVCalInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitUVCalClear (gpointer in);

/** Public: Init Smoothing */
static void 
ObitUVCalSmoothInit (ObitUVCal *in, ObitUVSel *sel, ObitUVDesc *desc, 
		     ObitErr *err);

/** Private: Set Class function pointers. */
static void ObitUVCalClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitUVCal* newObitUVCal (gchar* name)
{
  ObitUVCal* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitUVCalClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitUVCal));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitUVCalInit((gpointer)out);

 return out;
} /* end newObitUVCal */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitUVCalGetClass (void)
{
  return (gconstpointer)&myClassInfo;
} /* end  ObitUVCalGetClass */

/**
 * Make a deep copy of input object.
 * Copies are made of complex members including disk files; these 
 * will be copied applying whatever selection is associated with the input.
 * Parent class members are included but any derived class info is ignored.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Error stack, returns if not empty.
 * \return pointer to the new object.
 */
ObitUVCal* ObitUVCalCopy (ObitUVCal *in, ObitUVCal *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  gboolean oldExist;
  gchar *outName;

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return NULL;
  g_assert (ObitIsA(in, &myClassInfo));
  if (out) g_assert (ObitIsA(out, &myClassInfo));

  /* Create if it doesn't exist */
  oldExist = out!=NULL;
  if (!oldExist) {
    /* derive object name */
    outName = g_strconcat ("Copy: ",in->name,NULL);
    out = newObitUVCal(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /* copy/set this classes additions only if out newly created */
  if (!oldExist) {
  }
  
  return out;
} /* end ObitUVCalCopy */

/**
 * Make a shallow copy of a object.
 * The result will have pointers to the more complex members.
 * Parent class members are included but any derived class info is ignored.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \return pointer to the new object.
 */
ObitUVCal* ObitUVCalClone  (ObitUVCal *in, ObitUVCal *out)
{
  const ObitClassInfo *myClass, *ParentClass;
  gboolean oldExist;
  gchar *outName;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  if (out) g_assert (ObitIsA(out, &myClassInfo));

  /* Create if it doesn't exist */
  oldExist = out!=NULL;
  if (!oldExist) {
    /* derive object name */
    outName = g_strconcat ("Clone: ",in->name,NULL);
    out = newObitUVCal(outName);
    g_free(outName);
  }

  /* shallow copy any parent class members */
   myClass     = in->ClassInfo;
   ParentClass = myClass->ParentClass;
   g_assert ((ParentClass!=NULL) && (ParentClass->ObitClone!=NULL));
   ParentClass->ObitClone ((Obit*)in, (Obit*)out);

   if (!oldExist) { /* only copy ObitInfoList if just created */
     out->info = ObitInfoListUnref(out->info);
     out->info = ObitInfoListRef(in->info);
   }
     
   /* copy/set this classes additions */

  return out;
} /* end ObitUVCalClone */

/**
 * Creates necessary structures reading what ever calibration
 * files are needed.
 * Output descriptor modified to reflect data selection.
 * \param in      Object CFto initialize.
 * \param sel     Data selector.
 * \param inDesc  Input  data descriptor.
 * \param outDesc Output data descriptor (after transformations/selection).
 * \param err     ObitError stack.
 */
void ObitUVCalStart (ObitUVCal *in, ObitUVSel *sel, ObitUVDesc *inDesc, 
		     ObitUVDesc *outDesc, ObitErr *err)
{
  olong i, sid;
  gchar *routine = "ObitUVCalStart";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitIsA(inDesc, ObitUVDescGetClass()));
  g_assert (ObitIsA(outDesc, ObitUVDescGetClass()));
  g_assert (ObitIsA(sel, ObitUVSelGetClass()));

  /* Data MUST be in Time order */
  if (inDesc->isort[0] != 'T') {
    Obit_log_error(err, OBIT_Error, 
		   "ERROR: Data MUST be time ordered to calibrate/edit %s", 
		   in->name);
    return;
 }

  /* Reference descriptor */
  in->myDesc = ObitUVDescCopy(inDesc,in->myDesc, err);
  /* Needs update if data compressed */
  if (in->myDesc->inaxes[0]==1) {
    in->myDesc->inaxes[0] *= 3;
    in->myDesc->incs  *= 3;
    in->myDesc->incf  *= 3;
    in->myDesc->incif *= 3;
    in->myDesc->lrec   = sel->lrecUC;
    in->myDesc->nrparm = sel->nrparmUC;
  }
  /* Could be 2 words per vis */
  if (in->myDesc->inaxes[0]==2) {
    in->myDesc->inaxes[0] = 3*(in->myDesc->inaxes[0]/2);
    in->myDesc->incs  = 3*(in->myDesc->incs/2);
    in->myDesc->incf  = 3*(in->myDesc->incf/2);
    in->myDesc->incif = 3*(in->myDesc->incif/2);
    in->myDesc->lrec   = sel->lrecUC;
    in->myDesc->nrparm = sel->nrparmUC;
  }

  /* Copy Selector */
  in->mySel = ObitUVSelCopy(sel, in->mySel, err);

  /* initialize in */
  in->numSubA   = inDesc->numSubA;
  in->numStok   = inDesc->inaxes[inDesc->jlocs];
  in->numChan   = inDesc->inaxes[inDesc->jlocf];
  in->numChan   = sel->numberChann;
  in->bChan     = sel->startChann;
  in->eChan     = sel->startChann + sel->numberChann - 1;
  in->numIF     = inDesc->inaxes[inDesc->jlocif];
  in->numIF     = sel->numberIF;
  in->bIF       = sel->startIF;
  in->eIF       = sel->startIF + sel->numberIF - 1;
  in->dropSubA  = sel->dropSubA;

  /* Get source information - is there a source table, or get from header? */
  if (in->SUTable) { /* Read table */
    in->sourceList = ObitTableSUGetList ((ObitTableSU*)in->SUTable, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    /* If only one source, copy name to output descriptor */
    if ((sel->numberSourcesList==1) && (sel->sources!=NULL)) {
      sid = sel->sources[0];
      if ((sid>0) && (sid<=in->sourceList->number))
	strncpy (outDesc->object, in->sourceList->SUlist[sid-1]->SourceName, 17);
    }
  } else { /* Get needed info from UVdescriptor */
    in->sourceList = ObitSourceListCreate ("Single Source", 1);
    in->sourceList->SUlist[0]->SourID = -1;
    in->sourceList->SUlist[0]->RAMean  = inDesc->crval[inDesc->jlocr];
    in->sourceList->SUlist[0]->DecMean = inDesc->crval[inDesc->jlocd];
    /* Precess position */
    ObitPrecessUVJPrecessApp (inDesc, in->sourceList->SUlist[0]);
  }

  /* If only one source selected make sure no "SOURCE" 
     random parameter is written */
  if (((sel->numberSourcesList==1) || (in->SUTable==NULL)) && (outDesc->ilocsu>=0))
    strncpy (outDesc->ptype[outDesc->ilocsu], "REMOVED ", UVLEN_KEYWORD); 

  /* Create AntennaLists */
  in->antennaLists = g_malloc0(in->numSubA*sizeof(ObitAntennaList*));

  /* Read Polarization calibration from AN tables  */
  if ((sel->doPolCal)  || (sel->doBPCal)) {
    for (i=0; i<in->numSubA; i++) {
      in->antennaLists[i] = ObitTableANGetList ((ObitTableAN*)in->ANTables[i], err);
      if (err->error) Obit_traceback_msg (err, routine, in->name);
      
      /* Make sure polarization cal present */
      if (sel->doPolCal) {
	if (in->antennaLists[i]->polType == OBIT_UVPoln_NoCal)
	  Obit_log_error(err, OBIT_Error, "No polarization Cal info for %s", in->name);
	if (in->antennaLists[i]->polType == OBIT_UVPoln_Unknown)
	  Obit_log_error(err, OBIT_Error, "Unknown polarization Cal info for %s", in->name);
      }
    } /* End loop over subarrays */
  }
  if (err->error) return; /* Bail out if trouble */

  /* Initialize Flagging if needed */
  ObitUVCalFlagInit (in, sel, inDesc, err);
  
  /* Initialize Spectral smoothing if needed */
  ObitUVCalSmoothInit(in, sel, inDesc, err);
  
  /* Initialize baseline calibration if needed */
  ObitUVCalBaselineInit(in, sel, inDesc, err);
  
  /* Initialize  amp/phase/ delay/ rate calibration if needed */
  ObitUVCalCalibrateInit(in, sel, inDesc, err);
  
  /* Initialize bandpass calibration */
  ObitUVCalBandpassInit(in, sel, inDesc, err);
  
  /* Initialize polarization calibration */
  ObitUVCalPolarizationInit(in, sel, inDesc, err);
  
  /* Initialize Stokes translation/data selection */
  ObitUVCalSelectInit(in, sel, inDesc, outDesc, err);

  /* If calibrating then the units are Jy (AIPS needs 8 char) */
  if (in->doCal) strncpy (outDesc->bunit, "JY      ", UVLEN_VALUE); 

} /* end UVCalStart */

/**
 * Apply calibration, editing and selection.
 * Data should be uncompressed before calling.
 * \param in   Calibration Object.
 * \param recIn Input raw data array for a single visibility, 
 *        random parameters then visibility array. 
 * \param recOut Output calibrated, edited,transformed data.
 * \param err  ObitError stack.
 * \return TRUE if some of the data is valid, FALSE if none.
 */
gboolean  ObitUVCalApply (ObitUVCal *in, ofloat *recIn, 
			  ofloat *recOut, ObitErr *err)
{
  gboolean OK;
  ofloat *visIn, *visOut, time, scl;
  olong  i,nrparm, ant1, ant2;
  ObitUVDesc *desc;
  ObitUVSel *sel;

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return FALSE;
  g_assert (recIn!=NULL);
  g_assert (recOut!=NULL);

  /* Copy random parameters */
  desc = in->myDesc;
  nrparm = in->mySel->nrparmUC;
  for (i=0; i<nrparm; i++)  recOut[i] = recIn[i];

  /* Set visibility pointers */
  sel = in->mySel;
  visIn  = &recIn[nrparm];
  visOut = &recOut[nrparm];

  /* Get time and baseline */
  time = recIn[desc->iloct];
  ant1 = (olong)(recIn[desc->ilocb]/256);
  ant2 = (olong)(recIn[desc->ilocb] - ant1*256);

  /* Remove subarray info? */
  if (in->dropSubA) recIn[desc->ilocb] = ant2 + ant1*256;

  /* Is this visibility wanted? */
  OK = ObitUVCalWant (in, time, ant1, ant2, recIn, visIn, err);
  if (!OK) return OK;

  /* Apply calibration */
  if (in->doFlag) ObitUVCalFlag (in, time, ant1, ant2, recIn, visIn, err);
  if (in->doSmo)  ObitUVCalSmooth (in, time, ant1, ant2, recIn, visIn, err);
  if (in->doBL)   ObitUVCalBaseline(in, time, ant1, ant2, recIn, visIn, err);
  if (in->doCal)  ObitUVCalCalibrate(in, time, ant1, ant2, recIn, visIn, err);
  if (in->doBP)   ObitUVCalBandpass(in, time, ant1, ant2, recIn, visIn, err);
  if (in->doPol)  ObitUVCalPolarization(in, time, ant1, ant2, recIn, visIn, err);
  OK = ObitUVCalSelect(in, recIn, visIn, visOut,err);

  /* If selecting in IF, scale U,V.W */
  if (sel->startIF>1) {
    scl = desc->freqIF[sel->startIF-1] / desc->freq;
    recOut[desc->ilocu] *= scl;
    recOut[desc->ilocv] *= scl;
    recOut[desc->ilocw] *= scl;
  }

  /* Pass everything? */
  OK =  OK || sel->passAll;
  return OK;
} /* end ObitUVCalApply */

/**
 * Deletes structures and shuts down any open I/O
 * \param in   Calibration Object.
 * \param err  ObitError stack.
 */
void ObitUVCalShutdown (ObitUVCal *in, ObitErr *err)
{
  olong i;

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* return if inactive */

  /* Shutdown as needed */
  if (in->doFlag) ObitUVCalFlagShutdown (in, err);
  if (in->doBL)   ObitUVCalBaselineShutdown(in, err);
  if (in->doCal)  ObitUVCalCalibrateShutdown(in,  err);
  if (in->doBP)   ObitUVCalBandpassShutdown(in, err);
  if (in->doPol)  ObitUVCalPolarizationShutdown(in, err);

  /* Unref tables */
  in->SUTable = ObitTableUnref((ObitTable*)in->SUTable);
  in->BLTable = ObitTableUnref((ObitTable*)in->BLTable);
  in->BPTable = ObitTableUnref((ObitTable*)in->BPTable);
  in->FGTable = ObitTableUnref((ObitTable*)in->FGTable);
  in->CLTable = ObitTableUnref((ObitTable*)in->CLTable);
  in->SNTable = ObitTableUnref((ObitTable*)in->SNTable);
  in->CQTable = ObitTableUnref((ObitTable*)in->CQTable);
  if (in->ANTables) {
    for (i=0; i<in->numANTable; i++) {
      in->ANTables[i] = ObitTableUnref((ObitTable*)in->ANTables[i]);
    }
    g_free(in->ANTables); in->ANTables = NULL;
  }
} /* end ObitUVCalShutdown */

/**
 * Determine if data meets selection criteria
 * Largely adopted from AIPS DATGET.FOR
 * \param in    Calibration Object.
 * \param time  Time of datum
 * \param ant1  first antenna number of baseline
 * \param ant2  second antanna of baseline.
 * \param RP    Random parameter array
 * \param visIn 1 visibility as an array of floats
 * \param err   ObitError stack.
 */
gboolean ObitUVCalWant (ObitUVCal *in, ofloat time, olong ant1, olong ant2, 
			ofloat *RP, ofloat *visIn, ObitErr *err)
{
  gboolean OK = TRUE;
  ObitUVDesc *desc;
  ObitUVSel *sel;
  ofloat uvmin2, uvmax2, uvdis2;
  olong   kbase, FQID, SourID, iSubA;
   

   /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return FALSE;
  g_assert (ObitUVCalIsA(in));

  /* local pointers for structures */
  desc = in->myDesc;
  sel  = in->mySel;

  /* Is time in desired range? */
  if ((time < sel->timeRange[0]) || (time > sel->timeRange[1])) return FALSE;

  /* are these antennas desired? */
  if (!ObitUVSelWantAnt(sel, ant1)) return FALSE;
  if (!ObitUVSelWantAnt(sel, ant2)) return FALSE;

  /* Selection by auto/cross correlation */
  if (sel->corrType != 1) {
    /* Cross only */
    if ((sel->corrType==0) && (ant1==ant2)) return FALSE;
    /* Auto only */
    if ((sel->corrType==2) && (ant1!=ant2)) return FALSE;
  }

  /* Baseline and subarray number in data */
  kbase = (olong)RP[desc->ilocb];
  iSubA = 1 + (olong)(100.0*(RP[desc->ilocb] -(ofloat)kbase) + 0.1);
  /* Wanted? */
  if ((sel->SubA > 0) && (sel->SubA != iSubA)) return FALSE;

  /* Data FQ id */
  if (desc->ilocfq >= 0) FQID = RP[desc->ilocfq] + 0.1;
  else  FQID = 0;
  /* Wanted? */
  if ((sel->FreqID > 0) && (sel->FreqID != FQID)) return FALSE;

  /* Source ID */
  if (desc->ilocsu >= 0) SourID = RP[desc->ilocsu] + 0.1;
  else SourID = 0;
  /* Is this source wanted */
  if (!ObitUVSelWantSour(sel, SourID)) return FALSE;

  /* Check UV range */
  uvmin2 = sel->UVRange[0] * sel->UVRange[0];
  uvmax2 = sel->UVRange[1] * sel->UVRange[1];
  uvdis2 = RP[desc->ilocu]*RP[desc->ilocu] + RP[desc->ilocv]*RP[desc->ilocv];
  if ((uvdis2 < uvmin2) || (uvdis2 > uvmax2)) return FALSE;

  return OK;
} /* end ObitUVCalWant */

/**
 * Smooth data in frequency 
 * Adapted from AIPS SMOSP.FOR
 * \param in    Calibration Object.
 * \param time  Time of datum
 * \param ant1  first antenna number of baseline
 * \param ant2  second antanna of baseline.
 * \param RP    Random parameters array.
 * \param visIn 1 visibility as an array of floats
 * \param err   ObitError stack.
 */
void ObitUVCalSmooth (ObitUVCal *in, float time, olong ant1, olong ant2, 
		      ofloat *RP, ofloat *visIn, ObitErr *err)
{
  ofloat *vis;
  ObitUVDesc *desc;
  ObitUVSel  *sel;
  olong   i, j, j1, j2, l, ioff, ipol, iif, ifrq, kpol, indx, suprad, inxinc;
  ofloat  s, w, fblank = ObitMagicF();
  
  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitUVCalIsA(in));
  
  /* Pointer to visibility data portion of record */
  desc = in->myDesc;
  sel  = in->mySel;
  vis = &visIn[desc->nrparm];

  /* half width of convolution kernal */
  suprad = in->SmoothWidth;

  /* increment for ??? */
  inxinc = desc->incf;
  /* loop over IFs */
  for (iif= in->bIF; iif<=in->eIF; iif++) { /* loop 100 */
    ioff = (iif-1) * desc->incif;
    /* loop over polzns */
    for (ipol= 1; ipol<=in->numStok; ipol++) { /* loop 90 */
      kpol = (ipol-1) * desc->incs;

      /* loop over real/imaginary */
      for (i= 1; i<=2; i++) { /* loop 80 */

	/* copy data to temp array */
	indx = (ioff + kpol) + (in->bChanSmo-1)*desc->incf;

	for (ifrq= in->bChanSmo; ifrq<=in->eChanSmo; ifrq++) { /* loop 10 */
	  if (vis[indx+3-i] <= 0.0) { /* Flagged? Note: there is an error here in AIPS */
	    in->SmoothWork[ifrq-1] = fblank;
	  } else {
	    in->SmoothWork[ifrq-1] = vis[indx+i-1];
	  } 
	  indx = indx + inxinc;
	} /* end loop  L10:  */;

	/* convolve the data */
	indx = (ioff + kpol) + (in->bChan-1)*desc->incf;
	for (ifrq= in->bChan; ifrq<=in->eChan; ifrq++) { /* loop 30 */
	  j1 = MAX (ifrq - suprad, in->bChanSmo);
	  j2 = MIN (ifrq + suprad, in->eChanSmo);
	  s = 0.0;
	  w = 0.0;
	  for (j= j1; j<=j2; j++) { /* loop 20 */
	    if (in->SmoothWork[j-1] != fblank) {
	      l = abs(ifrq-j);
	      s = in->SmoothWork[j-1] * in->SmoothConvFn[l] + s;
	      w = in->SmoothConvFn[l] + w;
	    } 
	  } /* end loop  L20:  */;

	  /* result of smoothing */
	  if (w > 0.0) { /* good */
	    vis[indx+i-1] = s / w;
	  } else { /* bad */
	    vis[indx+i-1] = 0.0;
	    vis[indx+3-i] = 0.0; /* Note: there was also a bug in AIPS here */
	  } 
	  indx = indx + inxinc;
	} /* end loop  L30:  */;
      } /* end real/imag loop  L80:  */;
    } /* end poln loop  L90:  */;
  } /* end IF loop  L100: */;
} /* end ObitUVCalSmooth */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitUVCalClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitUVCalClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitUVCalClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitUVCalClassInfoDefFn (gpointer inClass)
{
  ObitUVCalClassInfo *theClass = (ObitUVCalClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitUVCalClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitUVCalClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitUVCalGetClass;
  theClass->newObit       = (newObitFP)newObitUVCal;
  theClass->ObitCopy      = (ObitCopyFP)ObitUVCalCopy;
  theClass->ObitClone     = (ObitCloneFP)ObitUVCalClone;
  theClass->ObitClear     = (ObitClearFP)ObitUVCalClear;
  theClass->ObitInit      = (ObitInitFP)ObitUVCalInit;
  theClass->ObitUVCalStart= (ObitUVCalStartFP)ObitUVCalStart;
  theClass->ObitUVCalApply= (ObitUVCalApplyFP)ObitUVCalApply;
  theClass->ObitUVCalShutdown = (ObitUVCalShutdownFP)ObitUVCalShutdown;

} /* end ObitUVCalClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param in Pointer to the object to initialize.
 */
void ObitUVCalInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitUVCal *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->thread      = newObitThread();
  in->info        = newObitInfoList(); 
  in->myStatus    = OBIT_Inactive;
  in->myDesc      = NULL;
  in->mySel       = newObitUVSel(in->name);
  in->doFlag      = FALSE;
  in->doSmo       = FALSE;
  in->doCal       = FALSE;
  in->doBP        = FALSE;
  in->doPol       = FALSE;
  in->flag        = NULL;
  in->baselineCal = NULL;
  in->ampPhaseCal = NULL;
  in->bandpassCal = NULL;
  in->polnCal     = NULL;
  in->ANTables    = NULL;
  in->numANTable  = 0;
  in->BLTable     = NULL;
  in->BPTable     = NULL;
  in->CQTable     = NULL;
  in->CLTable     = NULL;
  in->FGTable     = NULL;
  in->SNTable     = NULL;
  in->SmoothConvFn= NULL;
  in->SmoothWork  = NULL;
  in->sourceList  = NULL;
  in->antennaLists= NULL;

} /* end ObitUVCalInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * For some reason this wasn't build into the GType class.
 * \param  in Pointer to the object to deallocate.
 *           Actually it should be an ObitUVCal* cast to an Obit*.
 */
void ObitUVCalClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitUVCal *in = inn;
  olong i;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  in->thread      = ObitThreadUnref(in->thread);
  in->info        = ObitInfoListUnref(in->info);
  in->myDesc      = ObitUVDescUnref(in->myDesc);
  in->mySel       = ObitUVSelUnref(in->mySel);
  in->flag        = ObitUVCalFlagSUnref(in->flag);
  in->baselineCal = ObitUVCalBaselineSUnref(in->baselineCal);
  in->ampPhaseCal = ObitUVCalCalibrateSUnref(in->ampPhaseCal);
  in->bandpassCal = ObitUVCalBandpassSUnref(in->bandpassCal);
  in->polnCal     = ObitUVCalPolarizationSUnref(in->polnCal);
  in->BLTable     = ObitUnref(in->BLTable);
  in->BPTable     = ObitUnref(in->BPTable);
  in->CLTable     = ObitUnref(in->CLTable);
  in->FGTable     = ObitUnref(in->FGTable);
  in->SNTable     = ObitUnref(in->SNTable);
  in->CQTable     = ObitUnref(in->CQTable);
  if (in->ANTables) {
    for (i=0; i<in->numANTable; i++) in->ANTables[i] = ObitUnref(in->ANTables[i]);
    if (in->ANTables) g_free(in->ANTables);
  }
  if (in->SmoothConvFn) g_free(in->SmoothConvFn);
  if (in->SmoothWork) g_free(in->SmoothWork);
  in->sourceList = ObitSourceListUnref(in->sourceList);
  if (in->antennaLists) {
    for (i=0; i<in->numSubA; i++) 
      in->antennaLists[i] = ObitAntennaListUnref(in->antennaLists[i]);
    g_free(in->antennaLists);
  }

  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitUVCalClear */

/**
 * Initialize structures for Spectral smoothing.
 * Convolving function parameters 
 * \li (1)  = type of function, 0 => no smoothing, 1 => Hanning, 2 => Gaussian
 *            3 => boxcar, 4 => sin(x)/x.
 * \li (2) = width of function in channels
 * \li (3) = support of function in channels 
 * Adapted from the AIPS SETSM.FOR
 * \param in   Flag Object.
 * \param sel  Data selector.
 * \param desc Data descriptor.
 * \param err  ObitError stack.
 */
void ObitUVCalSmoothInit (ObitUVCal *in, ObitUVSel *sel, ObitUVDesc *desc, 
		    ObitErr *err)
{
  olong   i, n, lspect, iType, suprad;
  ofloat  fx, x, w;
  ofloat  widths[4] = {4.0, 2.0, 2.0, 3.0};
  ofloat  sups[4] = {1.0, 3.0, 1.0, 4.0};
  
  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitUVCalIsA(in));
  
  /* Are we doing it? -  to be sure check parameters */
  in->doSmo = FALSE;
  iType = sel->smooth[0] + 0.5;
  if (iType <= 0) return;
  in->doSmo = TRUE;

  /* Copy Selector information */
  in->bChan       = sel->startChann;
  in->eChan       = sel->startChann + sel->numberChann - 1;
  in->bIF         = sel->startIF;
  in->eIF         = sel->startIF + sel->numberIF - 1;
  in->smooth[0]   = sel->smooth[0];
  in->smooth[1]   = sel->smooth[1];
  in->smooth[2]   = sel->smooth[2];

  /* Copy descriptor information */
  in->numIF       = desc->inaxes[desc->jlocif];
  in->numChan     = desc->inaxes[desc->jlocf];
  in->numStok     = desc->inaxes[desc->jlocs];
  in->numSubA     = desc->numSubA;

  /* Size of convolution function */
  suprad = 0.1 + in->smooth[2] / 2.0;
  in->SmoothWidth = suprad;

  /* Allocate smoothing arrays */
  if (in->SmoothConvFn) g_free(in->SmoothConvFn);
  if (in->SmoothWork)   g_free(in->SmoothWork);
  in->SmoothConvFn = g_malloc0(in->SmoothWidth*sizeof(ofloat));
  in->SmoothWork   = g_malloc0(in->numChan*sizeof(ofloat));
  
  /* Check ranges */
  if (iType > 4) iType = 1;
  in->smooth[0] = iType;

  lspect = MAX (12, in->numChan);
  if ((in->smooth[1] < 0.5)  ||  (in->smooth[1] > lspect/3)) in->smooth[1] = widths[iType-1];
  if ((in->smooth[2] > 4.*sups[iType-1]*in->smooth[1])  || 
      (in->smooth[2] < in->smooth[1]))in->smooth[2] = sups[iType-1] * in->smooth[1];
  in->smooth[2] = 2.0 * suprad + 1.0;

  /* channel ranges for in->smoothing */
  in->bChanSmo = MAX (1, in->bChan - suprad);
  in->eChanSmo = MIN (in->numChan, in->eChan + suprad);

  /* initialize Smoothing colvolution function */
  for (i=0; i<in->SmoothWidth; i++) in->SmoothConvFn[i] = 0.0;

  /* Set up for computing function */
  n = 1 + suprad;
  fx = 2.0 / in->smooth[1];
  in->SmoothConvFn[0] = 1.0;

  /* compute look-up tables */
  w = in->SmoothConvFn[0];
  switch (iType) {
    case 1:  /* Hanning smooth */
    for (i= 2; i<=n; i++) { /* loop 20 */
      x = i - 1.0;
      in->SmoothConvFn[i-1] = MAX (0.0, 1.0-fx*x);
      w = w + 2 * in->SmoothConvFn[i-1];
    } /* end loop  L20:  */;
    break;
    
  case 2: /* Gaussian smooth */
    fx = -log(2.0) * fx * fx;
    for (i= 2; i<=n; i++) { /* loop 30 */
      x = i - 1.0;
      in->SmoothConvFn[i-1] = exp (fx * x * x);
      w = w + 2 * in->SmoothConvFn[i-1];
    } /* end loop  L30:  */;
    break;
    
  case 3:   /* boxcar smooth */
    fx = 1.0 / fx;
    for (i= 2; i<=n; i++) { /* loop 40 */
      x = i - 1.0;
      if (x < fx) {
	in->SmoothConvFn[i-1] = 1.0;
      } else if (x == fx) {
	in->SmoothConvFn[i-1] = 0.5;
      } 
      w = w + 2 * in->SmoothConvFn[i-1];
      } /* end loop  L40:  */;
      break;
      
  case 4:   /* sinc smooth */
    fx = 3.14159 * fx;
    for (i= 2; i<=n; i++) { /* loop 50 */
      x = (i - 1.0) * fx;
      in->SmoothConvFn[i-1] = sin(x) / x;
      w = w + 2 * in->SmoothConvFn[i-1];
    } /* end loop  L50:  */;
    break;
    
  default:
    g_assert_not_reached(); /* unknown, barf */
  }; /* end switch */
  
  /* normalize integral */
  if (w <= 0.0) w = 1.0;
  for (i= 0; i<n; i++)  in->SmoothConvFn[i] /=  w;
  
} /*  end ObitUVCalFlagInit */




