/* $Id$      */
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

#include "ObitOTFCal.h"
#include "ObitOTFDesc.h"
#include "ObitOTFSel.h"
#include "Obit2DLegendre.h"
#include "ObitTableOTFCal.h"
#include "ObitTableOTFSoln.h"
#include "ObitTableOTFBP.h"
#include "ObitOTFCalFlag.h"
#include "ObitOTFCalBandpass.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitOTFCal.c
 * ObitOTFCal class function definitions.
 * Calibration class for GBT/OTF data
 * This class is derived from the Obit base class.
 */

/*--------------- File Global Variables  ----------------*/
/** name of the class defined in this file */
static gchar *myClassName = "ObitOTFCal";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo structure ObitOTFCalClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitOTFCalClassInfo myClassInfo = {FALSE};

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitOTFCalInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitOTFCalClear (gpointer in);

/** Private: Update calibration arrays. */
static void ObitOTFCalUpdate (ObitOTFCal *in, ofloat time, ObitErr *err);

/** Private:  Read calibration for a new time into the internal arrays. */
static void ObitOTFCalNewTime (ObitOTFCal *in, ofloat time, ObitErr *err);

/** Private: Modify output descriptor for selection and setup */
static void ObitOTFCalSelectInit (ObitOTFCal *in, ObitOTFSel *sel, 
				  ObitOTFDesc *inDesc, 
				  ObitOTFDesc *outDesc, ObitErr *err);

/** Private: Select/translate data */
static gboolean ObitOTFCalSelect (ObitOTFCal *in, ofloat *recIn, ofloat *recOut, 
				  ObitErr *err);
/** Private: Set Class function pointers. */
static void ObitOTFCalClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitOTFCal* newObitOTFCal (gchar* name)
{
  ObitOTFCal* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitOTFCalClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitOTFCal));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitOTFCalInit((gpointer)out);

 return out;
} /* end newObitOTFCal */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitOTFCalGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitOTFCalClassInit();

  return (gconstpointer)&myClassInfo;
} /* end  ObitOTFCalGetClass */

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
ObitOTFCal* ObitOTFCalCopy (ObitOTFCal *in, ObitOTFCal *out, ObitErr *err)
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
    out = newObitOTFCal(outName);
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
} /* end ObitOTFCalCopy */

/**
 * Make a shallow copy of a object.
 * The result will have pointers to the more complex members.
 * Parent class members are included but any derived class info is ignored.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \return pointer to the new object.
 */
ObitOTFCal* ObitOTFCalClone  (ObitOTFCal *in, ObitOTFCal *out)
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
    out = newObitOTFCal(outName);
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
} /* end ObitOTFCalClone */

/**
 * Initialize Calibrator
 * Output descriptor modified to reflect data selection.
 * \param in      Object to initialize.
 * \param sel     Data selector.
 * \param inDesc  Input  data descriptor.
 * \param outDesc Output data descriptor (after transformations/selection).
 * \param err     ObitError stack.
 */
void ObitOTFCalStart (ObitOTFCal *in, ObitOTFSel *sel, ObitOTFDesc *inDesc, 
		      ObitOTFArrayGeom *geom, ObitOTFDesc *outDesc, ObitErr *err)
{
  ObitIOCode retCode;
  ofloat azMax, elMax, azNorm, elNorm;
  olong iDet, iCoef, iTerm;
  gchar *routine="ObitOTFCalStart";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitIsA(inDesc, ObitOTFDescGetClass()));
  g_assert (ObitIsA(outDesc, ObitOTFDescGetClass()));
  g_assert (ObitIsA(sel, ObitOTFSelGetClass()));

  /* Data MUST be in Time order */
  if (inDesc->isort[0] != 'T') {
    Obit_log_error(err, OBIT_Error, 
		   "ERROR: Data MUST be time ordered to calibrate/edit %s", 
		   in->name);
    return;
 }

  /* Reference descriptor */
  in->myDesc = ObitOTFDescUnref(in->myDesc);
  in->myDesc = ObitOTFDescRef(inDesc);

  /* Copy Selector */
  in->mySel = ObitOTFSelCopy(sel, in->mySel, err);

  /* initialize in */
  in->numStok   = inDesc->inaxes[inDesc->jlocs];
  in->numChan   = inDesc->inaxes[inDesc->jlocf];
  in->numFeed   = inDesc->inaxes[inDesc->jlocfeed];
  in->bChan     = sel->startChann;
  in->eChan     = sel->startChann + sel->numberChann - 1;
  in->doFlag    = sel->doFlag;
 

  /* Initialize calibration */
  if (sel->doCal) {
    /* Open calibration table, create row structure, get ndetect, npoly  */
    if (in->doSolnTable) { /* Soln */
      retCode = 
	ObitTableOTFSolnOpen ((ObitTableOTFSoln*)(in->SolnTable), OBIT_IO_ReadWrite, err);
      if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
	Obit_traceback_msg (err, routine, in->name);
      in->SolnTableRow = (Obit*)newObitTableOTFSolnRow((ObitTableOTFSoln*)(in->SolnTable));
      in->numDet =  ((ObitTableOTFSoln*)in->SolnTable)->numDet;
      in->numPoly = ((ObitTableOTFSoln*)in->SolnTable)->numPoly;
      in->numRow =  ((ObitTableOTFSoln*)in->SolnTable)->myDesc->nrow;
    } else {  /* Cal */
      retCode = 
	ObitTableOTFCalOpen ((ObitTableOTFCal*)(in->CalTable), OBIT_IO_ReadWrite, err);
      if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
	Obit_traceback_msg (err, routine, in->name);
      in->CalTableRow = (Obit*)newObitTableOTFCalRow((ObitTableOTFCal*)(in->CalTable));
      in->numDet  = ((ObitTableOTFCal*)in->CalTable)->numDet;
      in->numPoly = ((ObitTableOTFCal*)in->CalTable)->numPoly;
      in->numRow  = ((ObitTableOTFCal*)in->CalTable)->myDesc->nrow;
    }
    
    /* Allocate calibration arrays */
    in->CalApplyCal    = g_malloc0(in->numDet*sizeof(ofloat));
    in->CalApplyAdd    = g_malloc0(in->numDet*sizeof(ofloat));
    in->CalApplyMult   = g_malloc0(in->numDet*sizeof(ofloat));
    in->CalApplyWt     = g_malloc0(in->numDet*sizeof(ofloat));
    in->CalApplyPoly   = g_malloc0(in->numPoly*sizeof(ofloat));
    in->CalPriorCal    = g_malloc0(in->numDet*sizeof(ofloat));
    in->CalPriorAdd    = g_malloc0(in->numDet*sizeof(ofloat));
    in->CalPriorMult   = g_malloc0(in->numDet*sizeof(ofloat));
    in->CalPriorWt     = g_malloc0(in->numDet*sizeof(ofloat));
    in->CalPriorPoly   = g_malloc0(in->numPoly*sizeof(ofloat));
    in->CalFollowCal   = g_malloc0(in->numDet*sizeof(ofloat));
    in->CalFollowAdd   = g_malloc0(in->numDet*sizeof(ofloat));
    in->CalFollowMult  = g_malloc0(in->numDet*sizeof(ofloat));
    in->CalFollowWt    = g_malloc0(in->numDet*sizeof(ofloat));
    in->CalFollowPoly  = g_malloc0(in->numPoly*sizeof(ofloat));
    in->poly           = g_malloc0(in->numPoly*in->numDet*sizeof(ofloat));
    
    /* Initial times to trigger update of calibration arrays */
    in->CalTime       = -1.0e20;
    in->PriorCalTime  = -1.0e20;
    in->FollowCalTime = -1.0e20;
    
    /* find maximum az, el offset */
    azMax = elMax = -1.0;
    for (iDet=0; iDet<in->numDet; iDet++) {
      azMax = MAX (azMax, fabs(geom->azOffset[iDet]));
      elMax = MAX (elMax, fabs(geom->elOffset[iDet]));
    }
    
    /* Setup for Legendre polynomial atmospheric model */
    /* Az, el normalization factors */
    if (azMax>0.0) azNorm = 1.0 / azMax;
    else azNorm = 1.0;
    if (elMax>0.0) elNorm = 1.0 / elMax;
    else elNorm = 1.0;
    
    /* fill Legendre polynomial terms, Normalize az, el offsets to maximum = 1.0 */
    iTerm = 0;
    for (iDet=0; iDet<in->numDet; iDet++) {
      for (iCoef=0; iCoef<in->numPoly; iCoef++) {
	in->poly[iTerm++] = Obit2DLegendre (iCoef, azNorm*geom->azOffset[iDet], 
					    elNorm*geom->elOffset[iDet]);
      }
    }
  } else {
    /* Not calibrating - set some values */
    in->numDet = in->numChan * in->numStok * in->numFeed;
  } /* end initialize calibration */

  /* Initialize bandpass calibration */
  if (sel->doBP) {
    ObitOTFCalBandpassInit (in, sel, inDesc, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
  } /* end BP cal setup */

  /* Must have some detectors */
  if (in->numDet<=0) {
    Obit_log_error(err, OBIT_Error, 
		   "ERROR: NO detectors in %s", 
		   in->name);
    return;
  }

  /* more arrays */
  in->WantDetect     = g_malloc0(in->numDet*sizeof(gboolean));

  /* Initialize Stokes translation/data selection */
  ObitOTFCalSelectInit(in, sel, inDesc, outDesc, err);

  /* initialize flagging if requested */
  if (sel->doFlag) ObitOTFCalFlagInit (in, sel, inDesc, err);

  /* now active */
  in->myStatus = OBIT_Active;

} /* end ObitOTFCalStart */

/**
 * Apply calibration and selection
 * \param in      Calibration Object.
 * \param recIn   input data array for a single record
 * \param recOut  output data array for a single record, calibrated selected and translated.
 * \param err  ObitError stack.
 * \return TRUE if some of the data is valid, FALSE if none.
 */
gboolean  ObitOTFCalApply (ObitOTFCal *in, ofloat *recIn, ofloat *recOut, ObitErr *err)
{
  gboolean OK;
  ofloat *data, time, corr, dRa, dDec, fblank = ObitMagicF();
  olong  i, j, scan, incdatawt;
  olong nfeed, nstok, nchan, nscal, ifeed, istok, ichan;
  ofloat  add, mult, cal;
  gboolean doDataWt;
  ObitOTFDesc *desc;
  ObitOTFSel *sel;
  gchar *routine = "ObitOTFCalApply";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return FALSE;
  g_assert (recIn!=NULL);
  g_assert (recOut!=NULL);

  /* Is the data selected? If not return FALSE. */

  desc = in->myDesc;
  incdatawt = desc->incdatawt; /* increment in data-wt axis */
  doDataWt = incdatawt>1;      /* Have Data-Wt axis? */

  /* Set data pointers */
  sel = in->mySel;
  data  = &recIn[desc->ilocdata];

  /* see if new time - update cal. */
  time = recIn[desc->iloct];
  if (sel->doCal && (time > in->CalTime)) {
    ObitOTFCalUpdate(in, time, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, FALSE);
  }

  /* Is this in the desired time range? */
  if ((time<sel->timeRange[0]) || (time>sel->timeRange[1])) return FALSE;

  /* Is this in a selected scan? */
  scan = (olong)(recIn[desc->ilocscan]+0.5);
  if ((scan<sel->scans[0]) || (scan>sel->scans[1])) return FALSE;

  /* Apply calibration if needed */
  if (sel->doCal) {
    /* position offset - convert az, el offset to RA, dec 
     azOff is actually in Xel */
    dRa  = in->CalApplyAzoff*cos (DG2RAD*recIn[desc->ilocrot]) + 
      in->CalApplyEloff*sin (DG2RAD*recIn[desc->ilocrot]);
    dDec = in->CalApplyEloff*cos (DG2RAD*recIn[desc->ilocrot]) - 
      in->CalApplyAzoff*sin (DG2RAD*recIn[desc->ilocrot]);
    recIn[desc->ilocra]  += dRa;
    recIn[desc->ilocdec] += dDec;

    /* Gain and offset */
    OK = FALSE;
    /* VEGAS data different 1 or 2 polarizations, 1 or more channels, 1 or more feeds, possibly with XPol
       cal data presumed in order, stokes0, stokes2 per feed */
    if (desc->OTFType==OBIT_GBTOTF_VEGAS) {
      /* How many detectors? */
      nfeed = desc->inaxes[desc->jlocfeed];
      nstok = desc->inaxes[desc->jlocs];
      nscal = MIN (2, nstok);              /* Number of stokes in cal table */
      nchan = desc->inaxes[desc->jlocf];
      /* Loop over feed */
      for (ifeed=0; ifeed<nfeed; ifeed++) {
	/* Loop over Stokes */
	for (istok=0; istok<nstok; istok++) {
	  /* Loop over Channel - apply same correction to all channels */
	  /* Depends on istok */
	  if (istok==0) {  /* XX */
	    cal = in->CalApplyCal[ifeed*nscal];
	    add = in->CalApplyAdd[ifeed*nscal];
	    mult = in->CalApplyMult[ifeed*nscal];
	  } else if (istok==1) { /* YY */
	    cal = in->CalApplyCal[ifeed*nscal+1];
	    add = in->CalApplyAdd[ifeed*nscal+1];
	    mult = in->CalApplyMult[ifeed*nscal+1];
	  } else if (istok>=2) { /* Real(XY) or Imag(XY) */
	    cal = 0.0;
	    add = 0.0;
	    mult = sqrt(fabs(in->CalApplyMult[ifeed*nscal]*in->CalApplyMult[ifeed*nscal+1]));
	  }
	  for (ichan=0; ichan<nchan; ichan++) {
	    if ((*data!=fblank) && (add!=fblank) && (mult!=fblank)) {
	      /* Subtract cal value if on */
	      if (fabs(recIn[desc->iloccal])>0.0) *data -= cal ;
	      
	      /* Replace data with cal value? */
	      if (sel->replCal) *data = 0.0;
	      
	      /* Additive and multiplicative terms */
	      *data = mult * (*data + add);

	      /* polynomial atmosphere term */
	      corr = 0.0;
	      for (j=0; j<in->numPoly; j++) corr += in->poly[i*in->numPoly+j]*in->CalApplyPoly[j];
	      *data += corr;
	      
	      /* Weight */
	      if (doDataWt) *(data+1) *= in->CalApplyWt[i];
	      
	      OK = TRUE;
	    } else *data = fblank; /* End if valid */
	    data += incdatawt; /* next channel */
	  } /* end channel loop */
	} /* end Stokes loop */
      } /* end feed loop */
    } else {  /* Non VEGAS */
      for (i=0; i<in->numDet; i++) {
	if ((*data!=fblank) && (in->CalApplyAdd[i] != fblank) && (in->CalApplyCal[i] != fblank)) {
	  /* Subtract cal value if on */
	  if (fabs(recIn[desc->iloccal])>0.0) *data -= in->CalApplyCal[i];
	  
	  /* Replace data with cal value? */
	  if (sel->replCal) *data = 0.0;
	  
	  /* Additive and multiplicative terms */
	  *data = in->CalApplyMult[i] * (*data + in->CalApplyAdd[i]);
	  
	  /* polynomial atmosphere term */
	  corr = 0.0;
	  for (j=0; j<in->numPoly; j++) corr += in->poly[i*in->numPoly+j]*in->CalApplyPoly[j];
	  *data += corr;
	  
	  /* Weight */
	  if (doDataWt) *(data+1) *= in->CalApplyWt[i];
	  
	  OK = TRUE;
	} else {
	  *data = fblank;
	  if (doDataWt) *(data+1) = 0.0;
	}
	
	data += incdatawt; /* next detector */
      }
    } /* end of non VEGAS calibration */
  } /* end apply calibration */

  /* Bandpass calibration if requested */
  if (sel->doBP) ObitOTFCalBandpass (in, time, recIn, data, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, FALSE);

  /* Flagging if requested */
  if (sel->doFlag) ObitOTFCalFlag (in, time, recIn, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, FALSE);

  /* Copy descriptive data to output */
  for (i=0; i<desc->numDesc; i++) recOut[i] = recIn[i];

  /* Check cal-on */
  OK = (sel->keepCal || (recIn[desc->iloccal]==0.0));

  /* Select data to output */
  OK = OK && ObitOTFCalSelect (in, &recIn[desc->ilocdata], &recOut[desc->ilocdata], err);
  
  return OK;
  } /* end ObitOTFCalApply */

/**
 * Deletes structures and shuts down any open I/O
 * \param in   Calibration Object.
 * \param err  ObitError stack.
 * \return NULL pointer for deallocated ObitOTFCal
 */
ObitOTFCal* ObitOTFCalShutdown (ObitOTFCal *in, ObitErr *err)
{
  ObitIOCode retCode;
  gchar *routine="ObitOTFCalShutdown";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return in;
  g_assert (ObitIsA(in, &myClassInfo));

  /* return if inactive */
  if (in->myStatus==OBIT_Inactive) return ObitOTFCalUnref(in);

  /* Shutdown Calibration as needed */
  if (in->mySel->doCal) {
    /* Close calibration table */
    if (in->doSolnTable && (((ObitTableOTFSoln*)(in->SolnTable))->myStatus==OBIT_Active)) { /* Soln */
      retCode = 
	ObitTableOTFSolnClose ((ObitTableOTFSoln*)(in->SolnTable), err);
      if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
	Obit_traceback_val (err, routine, in->name, in);
      in->SolnTableRow = (Obit*)ObitTableOTFSolnRowUnref((ObitTableOTFSoln*)(in->SolnTable));
      
    } else if (((ObitTableOTFCal*)(in->CalTable))->myStatus==OBIT_Active) {  /* Cal */
      retCode = 
	ObitTableOTFCalClose ((ObitTableOTFCal*)(in->CalTable), err);
      if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
	Obit_traceback_val (err, routine, in->name, in);
      in->CalTableRow = (Obit*)ObitTableOTFCalRowUnref((ObitTableOTFCal*)(in->CalTable));
    }
  }

  /* Shutdown flagging if needed */
  if (in->mySel->doFlag) ObitOTFCalFlagShutdown (in, err);

  return ObitOTFCalUnref(in);
} /* end ObitOTFCalShutdown */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitOTFCalClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitOTFCalClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitOTFCalClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitOTFCalClassInfoDefFn (gpointer inClass)
{
  ObitOTFCalClassInfo *theClass = (ObitOTFCalClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitOTFCalClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitOTFCalClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitOTFCalGetClass;
  theClass->ObitClear     = (ObitClearFP)ObitOTFCalClear;
  theClass->ObitInit      = (ObitInitFP)ObitOTFCalInit;
  theClass->newObit       = (newObitFP)newObitOTFCal;
  theClass->ObitCopy      = (ObitCopyFP)ObitOTFCalCopy;
  theClass->ObitClone     = (ObitCloneFP)ObitOTFCalClone;
  theClass->ObitOTFCalStart= (ObitOTFCalStartFP)ObitOTFCalStart;
  theClass->ObitOTFCalApply= (ObitOTFCalApplyFP)ObitOTFCalApply;
  theClass->ObitOTFCalShutdown = (ObitOTFCalShutdownFP)ObitOTFCalShutdown;

} /* end ObitOTFCalClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param in Pointer to the object to initialize.
 */
void ObitOTFCalInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitOTFCal *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->thread         = newObitThread();
  in->info           = newObitInfoList(); 
  in->myStatus       = OBIT_Inactive;
  in->myDesc         = NULL;
  in->mySel          = newObitOTFSel(in->name);
  in->CalTable       = NULL;
  in->CalTableRow    = NULL;
  in->BPCalTable     = NULL;
  in->BPCalTableRow  = NULL;
  in->SolnTable      = NULL;
  in->SolnTableRow   = NULL;
  in->bandpassCal    = NULL;
  in->CalApplyCal    = NULL;
  in->CalApplyAdd    = NULL;
  in->CalApplyMult   = NULL;
  in->CalApplyWt     = NULL;
  in->CalApplyPoly   = NULL;
  in->CalPriorCal    = NULL;
  in->CalPriorAdd    = NULL;
  in->CalPriorMult   = NULL;
  in->CalPriorWt     = NULL;
  in->CalPriorPoly   = NULL;
  in->CalFollowCal   = NULL;
  in->CalFollowAdd   = NULL;
  in->CalFollowMult  = NULL;
  in->CalFollowWt    = NULL;
  in->CalFollowPoly  = NULL;
  in->poly           = NULL;
  in->WantDetect     = NULL;
} /* end ObitOTFCalInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * For some reason this wasn't build into the GType class.
 * \param  in Pointer to the object to deallocate.
 *           Actually it should be an ObitOTFCal* cast to an Obit*.
 */
void ObitOTFCalClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitOTFCal *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  in->thread      = ObitThreadUnref(in->thread);
  in->info        = ObitInfoListUnref(in->info);
  in->myDesc      = ObitOTFDescUnref(in->myDesc);
  in->mySel       = ObitOTFSelUnref(in->mySel);
  in->CalTable    = ObitTableOTFCalUnref((ObitTableOTFCal*)in->CalTable);
  in->CalTableRow = ObitTableOTFCalRowUnref((ObitTableOTFCalRow*)in->CalTableRow);
  in->BPCalTable    = ObitTableOTFBPUnref((ObitTableOTFBP*)in->BPCalTable);
  in->BPCalTableRow = ObitTableOTFBPRowUnref((ObitTableOTFBPRow*)in->BPCalTableRow);
  in->SolnTable   = ObitTableOTFSolnUnref((ObitTableOTFSoln*)in->SolnTable);
  in->SolnTableRow= ObitTableOTFSolnRowUnref((ObitTableOTFSolnRow*)in->SolnTableRow);
  in->bandpassCal = ObitOTFCalBandpassSUnref(in->bandpassCal);
  if (in->CalApplyCal)    g_free(in->CalApplyCal);   in->CalApplyCal    = NULL;
  if (in->CalApplyAdd)    g_free(in->CalApplyAdd);   in->CalApplyAdd    = NULL;
  if (in->CalApplyMult)   g_free(in->CalApplyMult);  in->CalApplyMult   = NULL;
  if (in->CalApplyWt)     g_free(in->CalApplyWt);    in->CalApplyWt     = NULL;
  if (in->CalApplyPoly)   g_free(in->CalApplyPoly);  in->CalApplyPoly   = NULL;
  if (in->CalPriorCal)    g_free(in->CalPriorCal);   in->CalPriorCal    = NULL;
  if (in->CalPriorAdd)    g_free(in->CalPriorAdd);   in->CalPriorAdd    = NULL;
  if (in->CalPriorMult)   g_free(in->CalPriorMult);  in->CalPriorMult   = NULL;
  if (in->CalPriorWt)     g_free(in->CalPriorWt);    in->CalPriorWt     = NULL;
  if (in->CalPriorPoly)   g_free(in->CalPriorPoly);  in->CalPriorPoly   = NULL;
  if (in->CalFollowAdd)   g_free(in->CalFollowAdd);  in->CalFollowAdd   = NULL;
  if (in->CalFollowCal)   g_free(in->CalFollowCal);  in->CalFollowCal   = NULL;
  if (in->CalFollowMult)  g_free(in->CalFollowMult); in->CalFollowMult  = NULL;
  if (in->CalFollowWt)    g_free(in->CalFollowWt);   in->CalFollowWt    = NULL;
  if (in->CalFollowPoly)  g_free(in->CalFollowPoly); in->CalFollowPoly  = NULL;
  if (in->poly)           g_free(in->poly);          in->poly           = NULL;
  if (in->WantDetect)     g_free(in->WantDetect);    in->WantDetect     = NULL;

  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitOTFCalClear */

/**
 * Update ObitOTFCal calibration tables for time time.
 * The current table is interpolated between the previous and following
 * sets of solutions.
 * If a new set of entries is needed from the Soln/Cal table they are read.
 * \param in   Calibrate Object.
 * \param time desired time in days
 * \param err  Error stack for messages and errors.
 */
static void ObitOTFCalUpdate (ObitOTFCal *in, ofloat time, ObitErr *err)
{
  ofloat  delta, wtt1, wtt2, wt1, wt2, fblank = ObitMagicF();
  gboolean bad, newcal;
  olong i;
  gchar *routine="ObitOTFCalUpdate";
 
 
  /* see if time for new table entry */
  if ((in->LastRowRead < in->numRow)  &&  (time > in->FollowCalTime)) {
    ObitOTFCalNewTime (in, time, err);
    if (err->error) Obit_traceback_msg (err, routine, "unspecified");
    newcal = TRUE;
  } else {
    newcal = FALSE;
  }
  /* debug */
  if (in->LastRowRead >= in->numRow) {
    newcal = FALSE;
  }


  /* see if calibration needs update; every 0.1 of solution interval. */
  delta = (time - in->CalTime);  
  if ((!newcal) &&  (delta <= 0.1*(in->FollowCalTime-in->PriorCalTime))) return;

  /* interpolate current calibration to time */
  in->CalTime = time;

  /* set interpolation weights proportional to time difference. */
  wtt1 = 0.0;
  if (time <= in->FollowCalTime) {
    if (in->FollowCalTime > in->PriorCalTime) 
      wtt1 = (in->FollowCalTime - time) / (in->FollowCalTime - in->PriorCalTime);
  } 
  wtt2 = 1.0 - wtt1;

  /* position calibration */
  in->CalApplyAzoff  = wtt1 * in->CalPriorAzoff + wtt2 * in->CalFollowAzoff;
  in->CalApplyEloff = wtt1 * in->CalPriorEloff+ wtt2 * in->CalFollowEloff;

  /* Interpolate additive and multiplicative terms */
  for (i=0; i<in->numDet;  i++) {
    wt1 = wtt1;
    wt2 = wtt2;
    bad = FALSE;
  
    /* Check for flagging */
    if (in->CalPriorAdd[i] == fblank) {
      wt1 = 0.0;
      wt2 = 1.0;
	}
    /* Check for flagging */
    if (in->CalFollowAdd[i] == fblank) {
      bad =  (wt1 <= 0.0);
      wt1 = 1.0;
      wt2 = 0.0;
    }

    /* make sure solutions valid */
    if (in->CalPriorCal[i]==fblank) {
      wt1 = 0.0;
      bad = bad || (wt2 <= 0.0);
    }
    if (in->CalFollowCal[i]==fblank) {
      wt2 = 0.0;
      bad = bad || (wt1 <= 0.0);
    }
    if (in->CalPriorAdd[i]==fblank) {
      wt1 = 0.0;
      bad = bad || (wt2 <= 0.0);
    }
    if (in->CalFollowAdd[i]==fblank) {
      wt2 = 0.0;
      bad = bad || (wt1 <= 0.0);
    }

    if (!bad) {
      in->CalApplyCal[i]  = wt1 *in->CalPriorCal[i]  + wt2 * in->CalFollowCal[i];
      in->CalApplyAdd[i]  = wt1 *in->CalPriorAdd[i]  + wt2 * in->CalFollowAdd[i];
      in->CalApplyMult[i] = wt1 *in->CalPriorMult[i] + wt2 * in->CalFollowMult[i];
      in->CalApplyWt[i]   = wt1 *in->CalPriorWt[i]   + wt2 * in->CalFollowWt[i];
    } else { /* bad - flag */
      in->CalApplyAdd[i]  = fblank;
      in->CalApplyMult[i] = fblank;
      in->CalApplyWt[i]   = 0.0;
    }
  } /* end additive and multiplicative loop */
  
  
  /* Interpolate Polynomial terms */
  for (i=0; i<in->numPoly;  i++) {
    wt1 = wtt1;
    wt2 = wtt2;
    bad = FALSE;
    
    /* Check for flagging */
    if (in->CalPriorPoly[i] == fblank) {
      wt1 = 0.0;
      wt2 = 1.0;
    }
    /* Check for flagging */
    if (in->CalFollowPoly[i] == fblank) {
      bad =  (wt1 <= 0.0);
      wt1 = 1.0;
      wt2 = 0.0;
    }
    if (!bad) {
      in->CalApplyPoly[i] = wt1 *in->CalPriorPoly[i]  + wt2 * in->CalFollowPoly[i];
    } else {
      in->CalApplyPoly[i] = 0.0;
    }
  } /* end polynomial loop */
  /* debug 
     if ((time>2.026e-4) && (time<2.029e-4)) {
     fprintf (stderr, "debug time %e prior %e follow %e wt %f %f\n",
     time, in->PriorCalTime, in->FollowCalTime, wtt1, wtt2);
     fprintf (stderr, "poly %e\n",in->CalApplyPoly[0]);
     }*/

} /* end ObitOTFCalUpdate */
    
/**
 * Read calibration for next time from Cal or Soln table.
 * \param in   Calibrate Object.
 * \param time desired time in days
 * \param in   Error stack for messages and errors.
 */
static void ObitOTFCalNewTime (ObitOTFCal *in, ofloat time, ObitErr *err)
{
  ObitIOCode retCode;
  ofloat fblank = ObitMagicF();
  olong i, j;
  olong  irow, limit;
  gboolean want, first;
  ObitTableOTFSoln *OTFSolnTable = NULL;
  ObitTableOTFSolnRow *OTFSolnTableRow = NULL;
  ObitTableOTFCal *OTFCalTable = NULL;
  ObitTableOTFCalRow *OTFCalTableRow = NULL;
  gchar *routine="ObitOTFCalNewTime";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;

  
  /* initialize Prior and Following arrays if first call */
  if (in->LastRowRead <= 0) {
    first = TRUE;
    for (i=0; i<in->numDet;  i++) in->CalPriorCal[i]  = fblank;
    for (i=0; i<in->numDet;  i++) in->CalPriorAdd[i]  = fblank;
    for (i=0; i<in->numDet;  i++) in->CalPriorMult[i] = fblank;
    for (i=0; i<in->numDet;  i++) in->CalPriorWt[i]   = fblank;
    for (i=0; i<in->numPoly; i++) in->CalPriorPoly[i] = fblank;
    in->CalPriorAzoff  = 0.0;
    in->CalPriorEloff = 0.0;
    in->PriorCalTime = -1.0e20;
    for (i=0; i<in->numDet;  i++) in->CalFollowCal[i]  = fblank;
    for (i=0; i<in->numDet;  i++) in->CalFollowAdd[i]  = fblank;
    for (i=0; i<in->numDet;  i++) in->CalFollowMult[i] = fblank;
    for (i=0; i<in->numDet;  i++) in->CalFollowWt[i]   = fblank;
    for (i=0; i<in->numPoly; i++) in->CalFollowPoly[i] = fblank;
    in->CalFollowAzoff  = 0.0;
    in->CalFollowEloff = 0.0;
    in->FollowCalTime = -1.0e20;
  } /* end of initialize on first call */
  else {
    
    /* Shuffle data from Following to Prior  */
    first = FALSE;
    for (i=0; i<in->numDet;  i++) in->CalPriorCal[i]  = in->CalFollowCal[i];
    for (i=0; i<in->numDet;  i++) in->CalPriorAdd[i]  = in->CalFollowAdd[i];
    for (i=0; i<in->numDet;  i++) in->CalPriorMult[i] = in->CalFollowMult[i];
    for (i=0; i<in->numDet;  i++) in->CalPriorWt[i]   = in->CalFollowWt[i];
    for (i=0; i<in->numPoly; i++) in->CalPriorPoly[i] = in->CalFollowPoly[i];
    in->CalPriorAzoff  = in->CalFollowAzoff;
    in->CalPriorEloff = in->CalFollowEloff;
    in->PriorCalTime   = in->FollowCalTime;
  }

  /* Handle OTFSoln/OTFCal table separately but do the same things for each */
  if (in->doSolnTable) {
    /* OTFSoln table  - set local pointers */
    OTFSolnTable = (ObitTableOTFSoln*)in->SolnTable;
    OTFSolnTableRow = (ObitTableOTFSolnRow*)in->SolnTableRow;
    
    /* Read until after selected time. */
    limit = MAX (1, in->LastRowRead+1);
    for (j= limit; j<=in->numRow; j++) {
      irow = j;
      retCode = ObitTableOTFSolnReadRow (OTFSolnTable, irow, OTFSolnTableRow, err);
      if (err->error) Obit_traceback_msg (err, routine, "Cal(OTFSoln) table");
      in->LastRowRead = irow;
      if (OTFSolnTableRow->status < 0) continue; /* entry flagged? */
      
      want = TRUE;
      /* skip if not wanted */
      if (!want) continue;
      
      /* save */
      in->FollowCalTime = OTFSolnTableRow->Time;
      for (i=0; i<in->numDet;  i++) in->CalFollowCal[i]  = OTFSolnTableRow->cal[i];
      for (i=0; i<in->numDet;  i++) in->CalFollowAdd[i]  = OTFSolnTableRow->add[i];
      for (i=0; i<in->numDet;  i++) in->CalFollowMult[i] = OTFSolnTableRow->mult[i];
      for (i=0; i<in->numDet;  i++) in->CalFollowWt[i]   = OTFSolnTableRow->wt[i];
      for (i=0; i<in->numPoly; i++) in->CalFollowPoly[i] = OTFSolnTableRow->poly[i];
      in->CalFollowAzoff  = OTFSolnTableRow->dAz;
      in->CalFollowEloff = OTFSolnTableRow->dEl;

      /* Is this after the desired time? */
      if (OTFSolnTableRow->Time > time) break;

      /* Shuffle to Prior */
      first = FALSE;
      for (i=0; i<in->numDet;  i++) in->CalPriorCal[i]  = in->CalFollowCal[i];
      for (i=0; i<in->numDet;  i++) in->CalPriorAdd[i]  = in->CalFollowAdd[i];
      for (i=0; i<in->numDet;  i++) in->CalPriorMult[i] = in->CalFollowMult[i];
      for (i=0; i<in->numDet;  i++) in->CalPriorWt[i]   = in->CalFollowWt[i];
      for (i=0; i<in->numPoly; i++) in->CalPriorPoly[i] = in->CalFollowPoly[i];
      in->CalPriorAzoff  = in->CalFollowAzoff;
      in->CalPriorEloff = in->CalFollowEloff;
      in->PriorCalTime   = in->FollowCalTime;
      
    } /* End loop over table */
	
  } else {
    /* OTFCal table - set local pointers */
    OTFCalTable = (ObitTableOTFCal*)in->CalTable;
    OTFCalTableRow = (ObitTableOTFCalRow*)in->CalTableRow;
    
    /* Read until after selected time. */
    limit = MAX (1, in->LastRowRead+1);
    for (j= limit; j<=in->numRow; j++) {
      irow = j;
      retCode = ObitTableOTFCalReadRow (OTFCalTable, irow, OTFCalTableRow, err);
      if (err->error) Obit_traceback_msg (err, routine, "Cal(OTFCal) Table");
      in->LastRowRead = irow;
      if (OTFCalTableRow->status < 0) continue; /* entry flagged? */
      
      want = TRUE;
      /* skip if not wanted */
      if (!want) continue;
      
      /* save */
      in->FollowCalTime = OTFCalTableRow->Time;
      for (i=0; i<in->numDet;  i++) in->CalFollowCal[i]  = OTFCalTableRow->cal[i];
      for (i=0; i<in->numDet;  i++) in->CalFollowAdd[i]  = OTFCalTableRow->add[i];
      for (i=0; i<in->numDet;  i++) in->CalFollowMult[i] = OTFCalTableRow->mult[i];
      for (i=0; i<in->numDet;  i++) in->CalFollowWt[i]   = OTFCalTableRow->wt[i];
      for (i=0; i<in->numPoly; i++) in->CalFollowPoly[i] = OTFCalTableRow->poly[i];
      in->CalFollowAzoff  = OTFCalTableRow->dAz;
      in->CalFollowEloff = OTFCalTableRow->dEl;

      /* Is this after the desired time? */
      if (OTFCalTableRow->Time > time) break;

      /* Shuffle to Prior */
      first = FALSE;
      for (i=0; i<in->numDet;  i++) in->CalPriorCal[i]  = in->CalFollowCal[i];
      for (i=0; i<in->numDet;  i++) in->CalPriorAdd[i]  = in->CalFollowAdd[i];
      for (i=0; i<in->numDet;  i++) in->CalPriorMult[i] = in->CalFollowMult[i];
      for (i=0; i<in->numDet;  i++) in->CalPriorWt[i]   = in->CalFollowWt[i];
      for (i=0; i<in->numPoly; i++) in->CalPriorPoly[i] = in->CalFollowPoly[i];
      in->CalPriorAzoff  = in->CalFollowAzoff;
      in->CalPriorEloff = in->CalFollowEloff;
      in->PriorCalTime   = in->FollowCalTime;
      
    } /* end loop over table entries  */;

  } /* end choice between Soln and Cal */
  
      /* if first call and Prior not filled, Shuffle to Prior */
  if (first) {
    for (i=0; i<in->numDet;  i++) in->CalPriorCal[i]  = in->CalFollowCal[i];
    for (i=0; i<in->numDet;  i++) in->CalPriorAdd[i]  = in->CalFollowAdd[i];
    for (i=0; i<in->numDet;  i++) in->CalPriorMult[i] = in->CalFollowMult[i];
    for (i=0; i<in->numDet;  i++) in->CalPriorWt[i]   = in->CalFollowWt[i];
    for (i=0; i<in->numPoly; i++) in->CalPriorPoly[i] = in->CalFollowPoly[i];
    in->CalPriorAzoff  = in->CalFollowAzoff;
    in->CalPriorEloff = in->CalFollowEloff;
    in->PriorCalTime  = in->FollowCalTime - 2.0/86400.0;
  }

  /* just to be sure something rational in times */
  if (in->PriorCalTime < -1000.0)  in->PriorCalTime  = time - 2.0/86400.0;
  if (in->FollowCalTime > 10000.0) in->FollowCalTime = time + 2.0/86400.0;
  
} /* end ObitOTFCalNewTime */

/**
 * Initialize structures for data selection.
 * Output descriptor modified to reflect data selection.
 * Selection by feed not yet implemented.
 * \param in   Object to initialize.
 * \li in->PolMode Stokes conversion type, 0=>none, 1=I, 2=V 
 * \param sel     Data selector.
 * \param inDesc  Input  data descriptor.
 * \param outDesc Output data descriptor (after transformations/selection).
 * \param err     ObitError stack.
 */
static void ObitOTFCalSelectInit (ObitOTFCal *in, ObitOTFSel *sel, 
				  ObitOTFDesc *inDesc, 
				  ObitOTFDesc *outDesc, ObitErr *err)
{
  olong   i, nmode=3;
  olong   pmode, nstok=1;
  gchar  chmode[3][5] = {"    ", "I   ", "V   "}; /* Only derived types are I, V */
  odouble crval=0.0;
  ofloat  cdelt=1.0;

  /* Copy Selector */
  in->mySel = ObitOTFSelCopy(sel, in->mySel, err);
  
  /* Compute translation parameters, find poln. mode */
  pmode = -1;
  for (i= 0; i<nmode; i++) 
    if (!strncmp(sel->Stokes, chmode[i], 4)) pmode = i;
	
  /* Unrecognized Stokes' */
  if (pmode < 0) {
    Obit_log_error(err, OBIT_Error, 
		   "Unknown Stokes request %s for %s", sel->Stokes, in->name);
    return;
  } 

  /* Save mode */
  in->PolMode = pmode;
  
  /* linear polarized data (x-y) can only be converted to I */
  if ((inDesc->crval[inDesc->jlocs] <= -4.5) && (pmode==2)){
    Obit_log_error(err, OBIT_Error, 
		   "CANNOT convert linear to circular poln for %s", in->name);
    return;
  }

  /* branch by pmode */
  switch (pmode) {
  case 0: /* no conversion    */
    nstok = outDesc->inaxes[outDesc->jlocs];
    crval = outDesc->crval[outDesc->jlocs]; /* Stokes index */
    cdelt = outDesc->cdelt[outDesc->jlocs]; /* stokes increment */
    break;

  case 1: /* I    */
    nstok = 1;
    crval = 1.0; /* Stokes index */
    cdelt = 1.0; /* stokes increment */
    break;
    
  case 2: /* V    */
    nstok = 1;
    crval = 4.0; /* Stokes index */
    cdelt = 1.0; /* stokes increment */
    break;
    

  default: /* should never get here */
    g_assert_not_reached();
    break;
  }; /* end switch on poln mode */

  /* How many polarizations out */
  in->mySel->numberPoln  = nstok;

  /* Feed (Detector) selection */
  for (i=0; i<in->numDet; i++) 
    in->WantDetect[i] = ObitOTFSelWantFeed(in->mySel, i);

  /* Update output descriptor */
  /* Stokes selection */
  outDesc->inaxes[outDesc->jlocs] = nstok;
  outDesc->crpix[outDesc->jlocs] = 1.0;
  outDesc->crval[outDesc->jlocs] = crval;
  outDesc->cdelt[outDesc->jlocs] = cdelt;

  /* Channel selection */
  outDesc->inaxes[outDesc->jlocf] = sel->numberChann;
  outDesc->crpix[outDesc->jlocf]  = 1.0;
  outDesc->crval[outDesc->jlocf]  = 
    inDesc->crval[inDesc->jlocf] + 
    (sel->startChann-inDesc->crpix[inDesc->jlocf]) * 
    inDesc->cdelt[inDesc->jlocf];

  /* Reindex output descriptor */
  ObitOTFDescIndex (outDesc);

} /* end ObitOTFCalSelectInit */


/**
 * Select data and translate to desired Stokes parameter
 * \param in     Calibration Object.
 * \li in->PolMode Stokes conversion type, 0=>none, 1=I, 2=V 
 * \param recIn  input data as an array of floats
 * \param recOut output data as an array of floats
 * \param err    ObitError stack.
 * \returns TRUE if at least some data is valid
 */
static gboolean ObitOTFCalSelect (ObitOTFCal *in, ofloat *recIn, ofloat *recOut, 
				  ObitErr *err)
{
  olong   i, lc, lf, ioff;
  gboolean   good, wantFeed;
  olong incdatawt;
  ObitOTFSel *sel = in->mySel;
  ObitOTFDesc *desc = in->myDesc;
  ofloat fblank = ObitMagicF();

  /* existing error? */
  if (err->error) return FALSE;
  
  good = FALSE;
  incdatawt = desc->incdatawt; /* increment in data-wt axis */

  /* if in->PolMode=0 and no channel selection, 
     no translation is to be done, just copy */
  if ((in->PolMode==0) && (in->numChan==sel->numberChann)) {
    /* Loop over data */
    for (lc=0; lc<in->numDet; lc++) { 
      if (!in->WantDetect[lc]) {
	recIn[lc*incdatawt]   = fblank;  /* Unwanted Feed? */
	recIn[lc*incdatawt+1] = 0.;  /* Unwanted Feed? */
      }
      recOut[lc*incdatawt]   = recIn[lc*incdatawt];
      recOut[lc*incdatawt+1] = recIn[lc*incdatawt+1];
      good = good || (recOut[lc]!=fblank); /* Check if valid */
    }
    return good;
  }

  /* loop checking and copying records. */
  /* Loop over Feed */
  for (lf=0; lf<sel->numberFeed; lf++) { 

    /* This feed wanted? */
    wantFeed  = in->WantDetect[lf];

    /* record offset in input */
    ioff = (sel->startChann - 2) * desc->incf + lf*desc->incfeed;

    /* Loop over channel */
    for (lc=0; lc<sel->numberChann; lc++) { 
      ioff += desc->incf;
      i = lf*desc->incfeed + lc*desc->incf; /* output record number */
     
      if (!wantFeed) {
	recIn[ioff]   = fblank;  /* Unwanted Feed? */
	recIn[ioff+1] = 0.0;  
      }

      /* branch by in->PolMode */
      switch (in->PolMode) {
      case 0: /* no conversion  - just copy  */
	recOut[i]   = recIn[ioff];
	recOut[i+1] = recIn[ioff+1];
	good = good || (recOut[i]!=fblank); /* Check if valid */
	i += desc->incs;
	if (in->numStok>1) { /* Two polarizations? */
	  recOut[i]   = recIn[ioff+desc->incs];
	  recOut[i+1] = recIn[ioff+desc->incs+1];
	  good = good || (recOut[i]!=fblank); /* Check if valid */
	  i += desc->incs;
	}
	break;
	
      case 1: /* I    =(R+L)/2 */
	if ((recIn[ioff]!=fblank) && (recIn[ioff+desc->incs]!=fblank)) {
	  recOut[i]   = 0.5 * (recIn[ioff] + recIn[ioff+desc->incs]);
	  recOut[i+1] = 0.5 * (recIn[ioff+1] + recIn[ioff+desc->incs+1]);
	}
	good = good || (recOut[i]!=fblank); /* Check if valid */
	i += desc->incs;
	break;
	
      case 2: /* V    = (R-L)/2 */
	if ((recIn[ioff]!=fblank) && (recIn[ioff+desc->incs]!=fblank))
	  recOut[i]   = 0.5 * (recIn[ioff] - recIn[ioff+desc->incs]);
	  recOut[i+1] = 0.5 * (recIn[ioff+1] - recIn[ioff+desc->incs+1]);
	good = good || (recOut[i]!=fblank); /* Check if valid */
	i += desc->incs;
	break;
	
      default: /* should never get here */
	g_assert_not_reached();
	break;
      }; /* end switch on poln mode */
      
    } /* end loop over channel  */;
  } /* end loop over Feed */
  return good;
} /* end ObitOTFCalSelect */

