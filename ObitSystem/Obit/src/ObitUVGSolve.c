/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2006-2010                                          */
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

#include "ObitUVSel.h"
#include "ObitUVGSolve.h"
#include "ObitMem.h"
#include "ObitTableUtil.h"
#include "ObitUVUtil.h"
#include "ObitUVSoln.h"
#include "ObitTableANUtil.h"
#include "ObitTableSUUtil.h"
#include "ObitPrecess.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitUVGSolve.c
 * ObitUVGSolve class function definitions.
 * This class enables self calibration of ObitUV data sets
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitUVGSolve";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo structure ObitUVGSolveClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitUVGSolveClassInfo myClassInfo = {FALSE};

/*--------------- File Global Variables  ----------------*/

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitUVGSolveInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitUVGSolveClear (gpointer in);

/** Private:  Read and average next solution interval */
static gboolean 
NextAvg (ObitUV *inUV, ofloat interv, 
	 gboolean avgif, gboolean avgpol, gboolean ampScalar, ofloat* antwt, 
	 ofloat uvrang[2], ofloat wtuv, olong numAnt, olong numFreq, 
         olong numIF, olong numPol, odouble* timec, ofloat* timei, olong* sid, 
         olong* fqid, ofloat* vis, olong *ant1, olong *ant2, olong *nextVisBuf, 
	 ObitErr* err);

/** Private: Solve for Gains for a solution interval */
static void 
doSolve (ofloat* vobs, olong *ant1, olong *ant2, olong numAnt, olong numIF, 
	 olong numPol, olong refant, gboolean avgif, gboolean avgpol, gboolean dol1, 
	 olong mode, olong minno, ofloat snrmin, olong prtlv, 
	 ofloat* creal, ofloat* cimag, ofloat* cwt, olong* refan, 
	 gboolean* gotant, ofloat *gain, ofloat *snr, olong *count, ObitErr *err);

/** Private: Determine SNR of solution */
static void   
calcSNR (ofloat* vobs, olong *ant1, olong *ant2, olong numBL, olong numAnt, 
	 ofloat* gain, ofloat* snr, ofloat closer[2][2], ofloat snrmin, odouble time, 
	 olong iif, olong ist, olong* count, olong prtlv, gchar* prtsou, ObitErr *err);  

/** Private: Gain Soln: Compute least squares gains */
static void 
gainCalc (ofloat* vobs, olong *ant1, olong *ant2, olong numBL, olong numAnt, olong refant, 
	  olong mode, olong minno, ofloat* g, olong* nref, olong prtlv, 
	  olong* ierr, ObitErr* err);

/** Private: Gain Soln: Does L1 solution for gains  */
static void 
gainCalcL1 (ofloat* vobs, olong *ant1, olong *ant2, olong numBL, olong numAnt, olong refant, 
	    olong mode, olong minno, ofloat* g, olong* nref, olong prtlv, 
	    olong* ierr, ObitErr* err);

/** Private: Set Class function pointers. */
static void ObitUVGSolveClassInfoDefFn (gpointer inClass);

/* Set Antenna and source lists */
static void SetLists (ObitUVGSolve *in, ObitUV *inUV, olong suba, 
		      ObitErr* err);
/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitUVGSolve* newObitUVGSolve (gchar* name)
{
  ObitUVGSolve* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitUVGSolveClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitUVGSolve));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitUVGSolveInit((gpointer)out);

 return out;
} /* end newObitUVGSolve */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitUVGSolveGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitUVGSolveClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitUVGSolveGetClass */

/**
 * Make a deep copy of an ObitUVGSolve.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitUVGSolve* ObitUVGSolveCopy  (ObitUVGSolve *in, ObitUVGSolve *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  gboolean oldExist;
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
    out = newObitUVGSolve(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  out->AList = ObitAntennaListCopy(in->AList, out->AList, err);
  out->SList = ObitSourceListCopy(in->SList, out->SList, err);
  out->UVFullRange[0] = in->UVFullRange[0];
  out->UVFullRange[1] = in->UVFullRange[1];

  return out;
} /* end ObitUVGSolveCopy */

/**
 * Make a copy of a object but do not copy the actual data
 * This is useful to create an UVGSolve similar to the input one.
 * \param in  The object to copy
 * \param out An existing object pointer for output, must be defined.
 * \param err Obit error stack object.
 */
void ObitUVGSolveClone  (ObitUVGSolve *in, ObitUVGSolve *out, ObitErr *err)
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
  out->AList = ObitAntennaListCopy(in->AList, out->AList, err);
  out->SList = ObitSourceListCopy(in->SList,  out->SList, err);
  out->UVFullRange[0] = in->UVFullRange[0];
  out->UVFullRange[1] = in->UVFullRange[1];

} /* end ObitUVGSolveClone */

/**
 * Creates an ObitUVGSolve 
 * \param name      An optional name for the object.
 * \return the new object.
 */
ObitUVGSolve* ObitUVGSolveCreate (gchar* name)
{
  ObitUVGSolve* out;

  /* Create basic structure */
  out = newObitUVGSolve (name);

  return out;
} /* end ObitUVGSolveCreate */

/**
 * Determine phase or amp & phase calibration for an UV dataset divided by
 * a source model.
 * If the output table previously exists, deselect any entries corresponding to 
 * data selected on inUV (to be fitted).
 * Routine translated from the AIPSish UVUTIL.FOR/SLFCAL 
 * \param in      Input self cal object. 
 * If only one source fitted in->curSource, has it's ID
 * Control parameters are on the info member.
 * \li "solnVer" OBIT_int   (1,1,1) Solution (SN) table to write; 0=> create new.
 * \li "subA"    OBIT_int   (1,1,1) Selected subarray (default 1)
 * \li "solInt"  OBIT_float (1,1,1) Solution interval (min). (default scan)
 * \li "refAnt"  OBIT_int   (1,1,1) Ref ant to use. (default 1)
 * \li "avgPol"  OBIT_bool  (1,1,1) True if RR and LL to be averaged (false)
 * \li "avgIF"   OBIT_bool  (1,1,1) True if all IFs to be averaged (false)
 * \li "ampScalar" OBIT_bool  (1,1,1) True ampscalar averaging of data wanted
 * \li "minSNR"  OBIT_float (1,1,1) Minimum acceptable SNR (5)
 * \li "doMGM"   OBIT_bool  (1,1,1) True then find the mean gain modulus (false)
 * \li "elevMGM" OBIT_float (1,1,1) Min. elevation to include in mean gain modulus
 * \li "solType" OBIT_string (4,1,1 Solution type '  ', 'L1',  (' ')
 * \li "solMode" OBIT_string (4,1,1 Solution mode: 'A&P', 'P', 'P!A', 'GCON' ('P')
 * \li "minNo"   OBIT_int   (1,1,1) Min. no. antennas. (default 4)
 * \li "antWt"   OBIT_float (*,1,1) Antenna weights. (default 1.0)
 * \li "UVR_Full"OBIT_float (2,1,1) Range of baseline lengths with full weight
 *                                  (lamda). If none is given then 
 *                                  derive one if possible.
 * \li "WtUV"    OBIT_float (1,1,1) Weight outside of UVRANG. (default 1.0)
 * \li "minOK"   OBIT_float (1,1,1) Minimum fraction of valid solutions (def 0.1)
 * \li "prtLv"   OBIT_int   (1,1,1) Print level (default no print)
 *
 * On output the following are set
 * \li "FractOK"    OBIT_float (1,1,1) The fraction of solutions which are OK
 * \param inUV   Input UV data. 
 * \param outUV  UV with which the output  SN is to be associated
 * \param sel    UV selector describing data fitted, used to deselect any 
 *               extant entries in output SN table.
 * \param err    Error/message stack, returns if error.
 * \return Pointer to the newly created SN object which is associated with outUV.
 */
ObitTableSN* ObitUVGSolveCal (ObitUVGSolve *in, ObitUV *inUV, ObitUV *outUV, 
			      ObitUVSel *sel, ObitErr *err)
{
  ObitTableSN *outSoln=NULL;
  ObitTableSNRow *row=NULL;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  olong i, k, nif, npoln, iAnt, numBL, SNver;
  olong nextVisBuf, cntmgm, iSNRow, numFreq, cntGood=0, cntPoss=0, cntBad=0;
  oint numPol, numIF, numAnt, suba, refant, minno, prtlv, mode;
  ofloat solInt, snrmin, uvrang[2], wtuv, summgm, FractOK, minOK=0.1;
  olong itemp, kday, khr, kmn, ksec, *refAntUse=NULL;
  olong *ant1=NULL, *ant2=NULL, *count=NULL;
  ofloat *antwt=NULL, *creal=NULL, *cimag=NULL, *cwt=NULL;
  ofloat *avgVis=NULL, *gain=NULL, *snr=NULL;
  gboolean avgpol, avgif, domgm, dol1, ampScalar, empty;
  gboolean done, good, oldSN, *gotAnt;
  gchar soltyp[5], solmod[5];
  odouble timec=0.0, timex;
  olong sid, fqid=0, sourid=0, lenEntry;
  ofloat timei=0.0, elevmgm, elev=90.0/57.296;
  ObitIOCode retCode;
  gchar *tname, *ModeStr[] = {"A&P", "P", "P!A"};
  gchar *routine = "ObitUVGSolveCal";
  
  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return outSoln;
  g_assert (ObitUVGSolveIsA(in));
  g_assert (ObitUVIsA(inUV));
  g_assert (ObitUVIsA(outUV));
  
  /* Get Solution interval from inUV and copy to SubScanTime*/
  solInt = 0.0;
  ObitInfoListGetTest(in->info, "solInt", &type, dim, &solInt);
  solInt /= 1440.0;  /* Convert to days */
  ObitInfoListAlwaysPut(inUV->info, "SubScanTime", OBIT_float, dim, &solInt);

  /* Min allowable OK fraction */
  minOK = 0.1;
  ObitInfoListGetTest(in->info, "minOK", &type, dim, &minOK);

  /* open UV data  */
  retCode = ObitUVOpen (inUV, OBIT_IO_ReadCal, err);
  if (err->error) Obit_traceback_val (err, routine, inUV->name, outSoln);
  
  /* Update frequency tables on inUV */
  if (!inUV->myDesc->freqArr) ObitUVGetFreq (inUV, err);
  if (err->error) Obit_traceback_val (err, routine, inUV->name, outSoln);
  /* Need array information */
  if (!inUV->myDesc->numAnt)   ObitUVGetSubA (inUV, err);
  if (err->error) Obit_traceback_val (err, routine, inUV->name, outSoln);
  
  /* Create output - version requested? */
  itemp = 0;
  ObitInfoListGetTest(in->info, "solnVer", &type, dim, &itemp);
  SNver = itemp;
  oldSN = SNver > 0;  /* Does SN table already exist? */
  /* 0=> make new */
  itemp = ObitTableListGetHigh (outUV->tableList, "AIPS SN") + 1;
  if (SNver<=0) SNver = itemp;
  tname = g_strconcat ("SN Calibration for: ", outUV->name, NULL);
  if (inUV->myDesc->jlocs>=0)
    numPol = MIN (2, inUV->myDesc->inaxes[inUV->myDesc->jlocs]);
  else numPol = 1;
  if (inUV->myDesc->jlocif>=0)
    numIF  = inUV->myDesc->inaxes[inUV->myDesc->jlocif];
  else numIF  = 1;
  outSoln = newObitTableSNValue(tname, (ObitData*)outUV, &SNver, OBIT_IO_ReadWrite, 
				numPol, numIF, err);
  g_free (tname);
  if (err->error) Obit_traceback_val (err, routine, inUV->name, outSoln);
  if (inUV->myDesc->jlocf>=0)
    numFreq  = inUV->myDesc->inaxes[inUV->myDesc->jlocf];
  else numFreq  = 1;

  /* If SN table previously existed, deselect values about to be redetermined. 
     get information from selector on inUV */
  if (oldSN) ObitUVSolnDeselSN (outSoln, sel->SubA, sel->FreqID, 
				sel->numberAntList, sel->ants, 
				sel->numberSourcesList, sel->sources, 
				sel->timeRange, err);
  if (err->error) Obit_traceback_val (err, routine, inUV->name, outSoln);
 
  /* Init mean gain modulus statistics */
  summgm = 0.0;
  cntmgm = 0;
  
  /* Which subarray? */
  suba = 1;
  ObitInfoListGetTest(in->info, "subA", &type, dim, &suba);
  /* Can only do one */
  Obit_retval_if_fail((suba>0 && suba<=inUV->myDesc->numSubA), err, outSoln,
		      "%s: MUST specify a single subarray for %s", 
		      routine, inUV->name);
  
  /* Create arrays */
  numAnt = inUV->myDesc->numAnt[suba-1];
  gain   = g_malloc0(2*numAnt*sizeof(ofloat));
  snr    = g_malloc0(numAnt*sizeof(ofloat));
  count  = g_malloc0(numAnt*sizeof(olong));
  antwt  = g_malloc0(numAnt*sizeof(ofloat));
  creal  = g_malloc0(numAnt*numIF*numPol*sizeof(ofloat));
  cimag  = g_malloc0(numAnt*numIF*numPol*sizeof(ofloat));
  cwt    = g_malloc0(numAnt*numIF*numPol*sizeof(ofloat));
  gotAnt = g_malloc0(numAnt*sizeof(gboolean));
  refAntUse = g_malloc0(numIF*numPol*sizeof(olong));
  numBL  = (numAnt * (numAnt-1)) / 2;
  ant1   = g_malloc0(numBL*sizeof(olong)+5);
  ant2   = g_malloc0(numBL*sizeof(olong)+5);
  
  /* Get parameters from inUV */
  refant = 1;
  ObitInfoListGetTest(in->info, "refAnt", &type, dim, &refant);
  avgpol = FALSE;
  ObitInfoListGetTest(in->info, "avgPol", &type, dim, &avgpol);
  if (numPol<=1) avgpol = FALSE;  /* Don't average if only one */
  avgif = FALSE;
  ObitInfoListGetTest(in->info, "avgIF",  &type, dim, &avgif);
  if (numIF<=1) avgif = FALSE;  /* Don't average if only one */
  ampScalar = FALSE;
  ObitInfoListGetTest(in->info, "ampScalar",  &type, dim, &ampScalar);
  snrmin = 5.0;
  ObitInfoListGetTest(in->info, "minSNR", &type, dim, &snrmin);
  domgm = FALSE;
  ObitInfoListGetTest(in->info, "doMGM", &type, dim, &domgm);
  elevmgm = 0.0;
  ObitInfoListGetTest(in->info, "elevMGM", &type, dim, &elevmgm);
  elevmgm /= 57.296;  /* to radians */
  soltyp[0] = soltyp[1] = soltyp[2] = soltyp[3] = ' '; soltyp[4] = 0;
  ObitInfoListGetTest(in->info, "solType", &type, dim, soltyp);
  solmod[0] = solmod[1] = solmod[2] = solmod[3] = ' '; solmod[4] = 0;
  ObitInfoListGetTest(in->info, "solMode", &type, dim, solmod);
  minno = 4;
  ObitInfoListGetTest(in->info, "minNo",  &type, dim, &minno);
  for (i=0; i<numAnt; i++) antwt[i] = 1.0;
  ObitInfoListGetTest(in->info, "antWt",  &type, dim, &antwt[0]);
  uvrang[0] = 0.0; uvrang[1] = 1.0e15;
  if (!ObitInfoListGetTest(in->info, "UVR_Full", &type, dim, &uvrang[0])
      || (uvrang[1]<=uvrang[0])) {
    /* If no explicit uv range given, use everything */
    uvrang[0] = 0.0;
    uvrang[1] = 1.0e15;
  }  /* end derive uv range */
  wtuv = 1.0;
  ObitInfoListGetTest(in->info, "WtUV", &type, dim, &wtuv);
  prtlv = 0;
  ObitInfoListGetTest(in->info, "prtLv", &type, dim, &prtlv);
  
  /* Digest SOLMODE and SOLTYPE */
  dol1 = !strncmp(soltyp, "L1",2);
  mode = 1;
  if (!strncmp(solmod, "A&P",3))  mode = 0;
  if (solmod[0]=='P')             mode = 1;
  if (!strncmp(solmod, "P!A", 3)) mode = 2;
  if (prtlv>=1) {
    Obit_log_error(err, OBIT_InfoErr, "Self calibrate in %s mode", ModeStr[mode]);
    ObitErrLog(err);
  }
  
  sid = in->curSource;  /* In case only one */
  /* If averaging gain modulus get antenna and source lists */
  if (domgm && (elevmgm>0.0)) SetLists (in, inUV, sel->SubA, err);
 
  /* Averaging of data? */
  nif = numIF;
  if (avgif) nif = 1;
  npoln = numPol;
  if (avgpol) npoln = 1;
  
  /* Allocate average visibility array */
  lenEntry = 4;
  avgVis  = g_malloc0(lenEntry*nif*npoln*numBL*sizeof(ofloat));
  
  /* Open output table */
  retCode = ObitTableSNOpen (outSoln, OBIT_IO_ReadWrite, err);
  if (err->error) goto cleanup;
  /* Anything already there? */
  empty = outSoln->myDesc->nrow==0;
  if (empty) {  /* Init if empty */
    outSoln->numAnt = numAnt;  /* Number of antennas */
    outSoln->mGMod  = 1.0;     /* initial mean gain modulus */
  }
  
  /* Create Row */
  row = newObitTableSNRow (outSoln);
  
  /* Attach row to output buffer */
  ObitTableSNSetRow (outSoln, row, err);
  if (err->error) goto cleanup;
  
  /* Initialize solution row */
  row->Time   = 0.0; 
  row->TimeI  = 0.0; 
  row->SourID = sid; 
  row->antNo  = 0; 
  row->SubA   = 0; 
  row->FreqID = 0; 
  row->IFR    = 0.0; 
  row->NodeNo = 0; 
  row->MBDelay1 = 0.0; 
  for (i=0; i<numIF; i++) {
    row->Real1[i]   = 0.0; 
    row->Imag1[i]   = 0.0; 
    row->Delay1[i]  = 0.0; 
    row->Rate1[i]   = 0.0; 
    row->Weight1[i] = 0.0; 
    row->RefAnt1[i] = 0; 
  }
  if (numPol>1) {
    row->MBDelay2 = 0.0; 
    for (i=0; i<numIF; i++) {
      row->Real2[i]   = 0.0; 
      row->Imag2[i]   = 0.0; 
      row->Delay2[i]  = 0.0; 
      row->Rate2[i]   = 0.0; 
      row->Weight2[i] = 0.0; 
      row->RefAnt2[i] = 0; 
    }
  }
  
  /* Loop until done */
  done = FALSE;
  nextVisBuf = -1;
  while (!done) {
    /* Read and average next solution interval */
    done =  NextAvg (inUV, solInt, avgif, avgpol, ampScalar, antwt, uvrang, wtuv, 
		     (olong)numAnt, numFreq, (olong)numIF, (olong)numPol, 
		     &timec, &timei, &sid, &fqid, avgVis, ant1, ant2, &nextVisBuf, 
		     err);
    if (err->error) goto cleanup;

    /* Done? 
    if (done) break; still have OK data */
    
    /* Write time if requested */
    if (prtlv >= 4) {
      kday = timec;
      timex = (timec - kday) * 24.;
      khr = timex;
      timex = (timex - khr) * 60.;
      kmn = timex;
      timex = (timex - kmn) * 60.;
      ksec = timex + 0.5;
      Obit_log_error(err, OBIT_InfoErr, " time=%d/ %3d %3d %3d,", 
		     kday, khr, kmn, ksec);
    } 
    
    /* Do solutions */
    doSolve (avgVis, ant1, ant2, numAnt, numIF, numPol, refant, avgif, avgpol, 
	     dol1, mode, minno, snrmin, prtlv, creal, cimag, cwt, refAntUse, 
	     gotAnt, gain, snr, count, err);

    /* Messages */
    if (prtlv>1) ObitErrLog(err);
    
    /* How many possible solutions */
    for (iAnt= 0; iAnt<numAnt; iAnt++) if (gotAnt[iAnt]) cntPoss += numIF*numPol;
    
    /* MGM with elevation limit? look up source index */
    if (domgm && (elevmgm>0.0)) {
      sourid = 0;
      for (k=0; k<in->SList->number; k++) {
	if (in->SList->SUlist[k]->SourID==sid) {sourid=k; break;}
      }
    }

    /* Write solutions to SN table */
    /* Common values */
    row->Time   = timec; 
    row->TimeI  = timei; 
    if (sid>0) row->SourID = sid; 
    row->SubA   = suba; 
    row->FreqID = fqid; 
    
    /* Loop over antennas */
    iSNRow = -1;
    for (iAnt= 0; iAnt<numAnt; iAnt++) {
      /* Need elevation? */
      if (domgm && (elevmgm>0.0)) 
	elev = ObitAntennaListElev (in->AList, iAnt+1, timec, 
				    in->SList->SUlist[sourid]);
      if (gotAnt[iAnt]) {
	good = FALSE; /* Until proven */
	row->antNo  = iAnt+1; 
	for (i=0; i<numIF; i++) {
	  row->Real1[i]   = creal[iAnt+i*numAnt]; 
	  row->Imag1[i]   = cimag[iAnt+i*numAnt]; 
	  row->Weight1[i] = cwt[iAnt+i*numAnt]; 
	  row->RefAnt1[i] = refAntUse[i]; 
	  if (cwt[iAnt+i*numAnt]>0.0) {good = TRUE; cntGood++;}
	  if (cwt[iAnt+i*numAnt]<=0.0) cntBad++;    /* DEBUG */
	  /* If good sum mean gain modulus */
	  if (domgm && (cwt[iAnt+i*numAnt]>0.0) && (elev>elevmgm)) {
	    cntmgm++;
	    summgm += sqrt(row->Real1[i]*row->Real1[i]+
			   row->Imag1[i]*row->Imag1[i]);
	  }
	}
	if (numPol>1) {
	  row->MBDelay2 = 0.0; 
	  for (i=0; i<numIF; i++) {
	    row->Real2[i]   = creal[iAnt+(i+numIF)*numAnt]; 
	    row->Imag2[i]   = cimag[iAnt+(i+numIF)*numAnt]; 
	    row->Weight2[i] = cwt[iAnt+(i+numIF)*numAnt]; 
	    row->RefAnt2[i] = refAntUse[numIF+i];
	    if (cwt[iAnt+(i+numIF)*numAnt]>0.0) {good = TRUE; cntGood++;}
	    if (cwt[iAnt+(i+numIF)*numAnt]<=0.0) cntBad++;    /* DEBUG */
	    /* If good sum mean gain modulus */
	    if (domgm && (cwt[iAnt+(i+numIF)*numAnt]>0.0) && (elev>elevmgm)) {
	      cntmgm++;
	      summgm += sqrt(row->Real2[i]*row->Real2[i]+
			     row->Imag2[i]*row->Imag2[i]);
	    }
	  }
	}
	/* DEBUG
	if ((row->Weight1[0]>1000.0) || (row->Weight1[1]>1000.0) || 
	    (row->Weight2[0]>1000.0) || (row->Weight2[1]>1000.0)) {
	  Obit_log_error(err, OBIT_InfoErr, "BINGO time %f and %d wt %f %f %f %f", 
			 row->Time, row->antNo, row->Weight1[0], row->Weight1[1], 
			 row->Weight2[0], row->Weight2[1]);
 	} */
	retCode = ObitTableSNWriteRow (outSoln, iSNRow, row, err);
	if (err->error) goto cleanup;
      } /* end if gotant */
    }
  } /* end loop processing data */
  
  /* Close input uv data */
  ObitUVClose (inUV, err);
  if (err->error) goto cleanup;
  
  /* Save mean gain modulus if requested */
  if (domgm  &&  (cntmgm > 0)  &&  ((mode == 0) || (mode == 3))) {
    /* Update */
    outSoln->mGMod  = summgm / cntmgm;
  }
  
  /* If suba>1 or not empty at start mark as unsorted */
  if ((suba == 1) && empty) {  
    /* IF subarray 1, new table the sort is time, antenna */
    outSoln->myDesc->sort[0] = outSoln->TimeCol+1;
    outSoln->myDesc->sort[1] = outSoln->antNoCol+1; 
  } else { /* otherwise unsorted */
    outSoln->myDesc->sort[0] = 0;
    outSoln->myDesc->sort[1] = 0;
  }

  /* Close output table */
  retCode = ObitTableSNClose (outSoln, err);
  if (err->error) goto cleanup;
  
  /* Give success rate */
  if (prtlv>=3) {
    Obit_log_error(err, OBIT_InfoErr, " %d of %d possible solutions found",
		   cntGood, cntPoss);
  }
  FractOK = (ofloat)cntGood / (ofloat)cntPoss;
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut(in->info, "FractOK", OBIT_float, dim, &FractOK);
  /* DEBUG Obit_log_error(err, OBIT_InfoErr, " %d of %d possible solutions bad",
     cntBad, cntPoss); */
  
  /* Require at least minOK */
  if (!(cntGood >= minOK*cntPoss)) {
    Obit_log_error(err, OBIT_Error,
		   "%s: TOO FEW Successful selfcal solutions for  %s", 
		   routine, inUV->name);
    goto cleanup;
  }
  
 goto cleanup; /* Cleanup */
  cleanup: row = ObitTableSNUnref(row);
  if (gain)   g_free(gain);
  if (snr)    g_free(snr);
  if (count)  g_free(count);
  if (antwt)  g_free(antwt);
  if (creal)  g_free(creal);
  if (cimag)  g_free(cimag);
  if (cwt)    g_free(cwt);
  if (gotAnt) g_free(gotAnt);
  if (avgVis) g_free(avgVis);
  if (ant1)   g_free(ant1);
  if (ant2)   g_free(ant2);
  if (refAntUse) g_free(refAntUse);
  if (err->error) Obit_traceback_val (err, routine, inUV->name, outSoln);
  
  return outSoln;
} /* end ObitUVGSolveCal */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitUVGSolveClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitUVGSolveClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitUVGSolveClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitUVGSolveClassInfoDefFn (gpointer inClass)
{
  ObitUVGSolveClassInfo *theClass = (ObitUVGSolveClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitUVGSolveClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitUVGSolveClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitUVGSolveGetClass;
  theClass->newObit       = (newObitFP)newObitUVGSolve;
  theClass->ObitCopy      = (ObitCopyFP)ObitUVGSolveCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitUVGSolveClear;
  theClass->ObitInit      = (ObitInitFP)ObitUVGSolveInit;
  theClass->ObitUVGSolveCreate = (ObitUVGSolveCreateFP)ObitUVGSolveCreate;
  theClass->ObitUVGSolveCal    = (ObitUVGSolveCalFP)ObitUVGSolveCal;

} /* end ObitUVGSolveClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitUVGSolveInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitUVGSolve *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->thread    = newObitThread();
  in->info      = newObitInfoList(); 
  in->curSource = 0;
  in->AList     = NULL;
  in->SList     = NULL;
  in->UVFullRange[0] = 0.0;
  in->UVFullRange[1] = 0.0;

} /* end ObitUVGSolveInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitUVGSolve* cast to an Obit*.
 */
void ObitUVGSolveClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitUVGSolve *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  in->thread   = ObitThreadUnref(in->thread);
  in->info     = ObitInfoListUnref(in->info);
  in->AList    = ObitAntennaListUnref(in->AList);
  in->SList    = ObitSourceListUnref(in->SList);
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitUVGSolveClear */

/**
 * Average next solution interval
 *  Data is averaged until one of several conditions is met:
 *  \li the time exceeds the initial time plus the specified interval.
 *  \li the source id (if present) changes
 *  \li the FQ id (if present) changes. 
 * Routine translated from the AIPSish UVUTIL.FOR/NXTAVG
 * \param inUV    Input UV data. 
 * \param interv  Desired time interval of average in days. 
 *                Will use value from inUV->mySel.ObitUVSelSubScan
                  if available, 0=> scan average
 * \param avgif   If true average in IF 
 * \param avgpol  If true average in polarization (Stokes 'I')
 * \param antwt   Extra weights to antennas (>=0 => 1.0)
 * \param uvrang  Range of baseline lengths with full weight (lamda). 
 *                0s => all baselines 
 * \param wtuv    Weight outside of UVRANG. (No default)
 * \param numAnt  Highest antenna number (NOT number of antennas)
 * \param numIF   Maximum number of Frequencies 
 * \param numIF   Maximum number of IFs 
 * \param numPol  Maximum number of polarizations
 * \param timec   [out] Center time of observations (days)
 * \param timei   [out] Actual time interval (days)
 * \param sid     [out] Source Id if present else -1. 
 * \param fqid    [out] FQ id if present else -1.
 * \param avgVis  [out] [4,baseline,IF,pol]
 * \param ant1    [out] Array of first antenna numbers (1-rel) on a baseline
 * \param ant2    [out] Array of second antenna numbers (1-rel) on a baseline
 * \param nextVisBuf [in/out] next vis (0-rel) to read in buffer, 
 *                   -1=> read, -999 done
 * \param err    Error/message stack, returns if error.
 * \return TRUE if all data read, else FALSE
 */
static gboolean 
NextAvg (ObitUV* inUV, ofloat interv, 
	 gboolean avgif, gboolean avgpol, gboolean ampScalar, ofloat* antwt, 
	 ofloat uvrang[2], ofloat wtuv, olong numAnt, olong numFreq, olong numIF, 
	 olong numPol, odouble* timec, ofloat* timei, olong* sid, olong* fqid, 
	 ofloat* avgVis, olong *ant1, olong *ant2, olong *nextVisBuf, ObitErr* err)
{
  gboolean done=FALSE;
  ObitIOCode retCode= OBIT_IO_OK;
  ofloat linterv, ctime, cbase, stime, ltime=0, weight, wt, bl, *visPnt, temp;
  ofloat tmp1, tmp2;
  olong i, j, csid, cfqid, a1, a2, *blLookup=NULL, blIndex, visIndex, lenEntry;
  olong iFreq, iIF, iStok, jBL, jIF, jStok, jncs, jncif, mPol, mIF, offset, numBL;
  odouble timeSum;
  olong timeCount, accumIndex;
  olong maxAllow; /* DEBUG */
  /* olong maxShit=0; DEBUG */
  gboolean gotData;
  gchar *routine = "ObitUVGSolve:NextAvg";
  
  /* Error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return done;  /* previous error? */
  g_assert(ObitUVIsA(inUV));
  
  /* Did we hit the end last time? */
  if (*nextVisBuf==-999) return TRUE;

  /* Get actual average time from Selector if available */
  temp = ObitUVSelSubScan (inUV->mySel);
  if ((temp>0.0)  && (temp<1.0)) linterv = 1.01*temp;
  else linterv = 0.99*interv;
  /* If still zero use one second */
  if (linterv <=0.0) linterv = 1.0 / 86400.0;

  /* Numbers of poln and IFs in accumulator */
  if (avgpol) mPol = 1;
  else mPol = numPol;
  if (avgif) mIF = 1;
  else mIF = numIF;

  /* Increments in accumulator */  
  lenEntry = 4;  /* Length of accumulator */
  jncif = lenEntry*(numAnt*(numAnt-1))/2;
  if (avgpol) jncs  = jncif;
  else jncs  = jncif * mIF;
  
  maxAllow = lenEntry*((numAnt*(numAnt-1))/2)*mPol*mIF; /* DEBUG */
  
  /* Baseline lookup table */
  blLookup = g_malloc0(numAnt*sizeof(olong));
  blLookup[0] = 0;
  for (i=1; i<numAnt; i++) blLookup[i] = blLookup[i-1] + numAnt-i;
  
  /* Fill in antenna numbers in ant* arrays - assume a1<a2 */
  jBL = 0;
  for (i=0; i<numAnt; i++) {
    for (j=i+1; j<numAnt; j++) {
      ant1[jBL] = i+1;
      ant2[jBL] = j+1;
      jBL++;
    }
  }
  
  /* Zero accumulations */
  numBL = (numAnt*(numAnt-1)) / 2;
  for (jBL=0; jBL<numBL; jBL++) {         /* Loop over baseline */
    for (jStok=0; jStok<mPol; jStok++) {  /* Loop over Pol */
      for (jIF=0; jIF<mIF; jIF++) {       /* Loop over IF */
	/* Accumulator index */
	accumIndex = jBL*lenEntry + jStok*jncs + jIF*jncif;
	for (i=0; i<lenEntry; i++) avgVis[accumIndex+i] = 0.0;
      } /* end loop over IF */
    } /* end loop over Pol */
  } /* end loop over baseline */
  *fqid = -1;
  stime = -1.0e20;
  timeSum = 0.0;
  timeCount = 0;
  
  /* Loop reading and averaging data */
  while (!done) {
    
    /* Need to read? */
    if ((*nextVisBuf<0) || (*nextVisBuf>=inUV->myDesc->numVisBuff)) {
      gotData = FALSE;
      while (!gotData) {
	retCode = ObitUVReadSelect (inUV, NULL, err);
	if (err->error) goto cleanup;
	/* Finished? */
	if (retCode==OBIT_IO_EOF) {
	  done = TRUE;
	  *nextVisBuf = -999;
	  break;
	}
	/* Make sure there was actual data */
	gotData = inUV->myDesc->numVisBuff>0;
      }
      *nextVisBuf = 0;
    } 

    /* Get actual average time from Selector if available */
    temp = ObitUVSelSubScan (inUV->mySel);
    if ((temp>0.0)  && (temp<1.0)) linterv = 1.01*temp;
    else linterv = 0.99*interv;
    /* If still zero use one second */
    if (linterv <=0.0) linterv = 1.0 / 86400.0; 

    /* Visibility pointer */
    visPnt = inUV->buffer + (*nextVisBuf) * inUV->myDesc->lrec;
    /* Vis data */
    ctime = visPnt[inUV->myDesc->iloct]; /* Time */
    if (stime<-1000.0) stime = ctime;  /* Set start time from first vis */
    if (inUV->myDesc->ilocsu>=0) csid  = (olong)visPnt[inUV->myDesc->ilocsu]; /* Source */
    else csid  = 0;
    if (timeCount==0) *sid = csid;  /* Set output source id first vis */
    if (inUV->myDesc->ilocfq>=0) cfqid = (olong)visPnt[inUV->myDesc->ilocfq]; /* FQid */
    else cfqid  = 0;
    if (*fqid<0) *fqid = cfqid;  /* Set output fq id first vis */
    cbase = visPnt[inUV->myDesc->ilocb]; /* Baseline */
    
    /* Is this integration done? */
    if ((*sid!=csid) || (*fqid!=cfqid) || (ctime-stime>=linterv)) break;
    
    /* Sum time */
    timeSum += ctime;
    timeCount++;
    ltime = ctime;   /* Last time in accumulation */
    
    /* crack Baseline */
    a1 = (cbase / 256.0) + 0.001;
    a2 = (cbase - a1 * 256) + 0.001;

    /* Check that antenna numbers in range */
    if (!((a1>0) && (a2>0) && (a1<=numAnt)&& (a2<=numAnt))) {
     Obit_log_error(err, OBIT_Error,
		    "%s: Bad antenna number, %d or %d not in [1, %d] in %s", 
		    routine, a1, a2, numAnt, inUV->name);
     goto cleanup;
    }

    /* Ignore autocorrelations */
    if (a1==a2) {(*nextVisBuf)++; continue;}
   
    /* Set extra weighting factors */
    weight = antwt[a1-1]*antwt[a2-1];
    
    /* Check baseline length - posibly modify weight */
    bl = sqrt (visPnt[inUV->myDesc->ilocu]*visPnt[inUV->myDesc->ilocu] + 
	       visPnt[inUV->myDesc->ilocv]*visPnt[inUV->myDesc->ilocv]);
    if ((bl<uvrang[0]) || (bl>uvrang[1])) weight *= wtuv;
    
    /* Baseline index this assumes a1<a2 always */
    blIndex =  blLookup[a1-1] + a2-a1-1;
    
    /* Accumulate */
    offset = lenEntry*blIndex; /* Offset in accumulator */
    /* Loop over polarization */
    for (iStok = 0; iStok<numPol; iStok++) {
      jStok = iStok;      /* Stokes in accumulation */
      if (avgpol) jStok = 0;
      /* Loop over IF */
      for (iIF = 0; iIF<numIF; iIF++) {
	jIF = iIF;        /* IF in accumulation */
	if (avgif) jIF = 0;
	
	/* Loop over frequency */
	for (iFreq = 0; iFreq<numFreq; iFreq++) {

	  /* Visibity index */
	  visIndex = inUV->myDesc->nrparm + iStok*inUV->myDesc->incs + 
	    iIF*inUV->myDesc->incif + iFreq*inUV->myDesc->incf;
	  
	  /* Accumulator index */
	  accumIndex = offset + jStok*jncs + jIF*jncif;
	  
	  /* Accumulate */
	  if (visPnt[visIndex+2]>0.0) {
	    wt = weight * visPnt[visIndex+2];
	    /* maxShit = MAX (maxShit, accumIndex);  DEBUG */
	    /* DEBUGif (accumIndex>maxAllow) { 
	       fprintf (stderr,"bad accum %d a1 %d a2 %d jStok %d jIF %d\n", 
	       accumIndex, a1, a2, jStok, jIF);
	       }  end DEBUG */
	    avgVis[accumIndex]   += wt *  visPnt[visIndex];
	    avgVis[accumIndex+1] += wt * visPnt[visIndex+1];
	    avgVis[accumIndex+2] += wt;
	    if (ampScalar) { /* Ampscalar averaging? Accumulate amplitude */
	      avgVis[accumIndex+3] += wt * sqrt (visPnt[visIndex]*visPnt[visIndex] +
						 visPnt[visIndex+1]*visPnt[visIndex+1]);
	    }
	  } /* End vis valid */
	} /* Loop over Freq */
      } /* Loop over IF */
    } /* Loop over Stokes */
    
    (*nextVisBuf)++;  /* Increment vis being processed */
  } /* End loop summing */
  
  /* Normalize by sum of weights */
  for (jBL=0; jBL<numBL; jBL++) {         /* Loop over baseline */
    for (jStok=0; jStok<mPol; jStok++) {  /* Loop over Pol */
      for (jIF=0; jIF<mIF; jIF++) {       /* Loop over IF */
	/* Accumulator index */
	accumIndex = jBL*lenEntry + jStok*jncs + jIF*jncif;
	if (avgVis[accumIndex+2]>0) {
	  avgVis[accumIndex]  /= avgVis[accumIndex+2];
	  avgVis[accumIndex+1]/= avgVis[accumIndex+2];
	  /* if ampScalar replace vector amp with scalar */
	  if (ampScalar) { 
	    tmp1 = sqrt (avgVis[accumIndex]  *avgVis[accumIndex] + 
			 avgVis[accumIndex+1]*avgVis[accumIndex+1]) *
	      avgVis[accumIndex+2];
	    if (tmp1>1.0e-20) tmp2 = avgVis[accumIndex+3] / tmp1;
	    else tmp2 = 1.0;
	    avgVis[accumIndex]   *= tmp2;
	    avgVis[accumIndex+1] *= tmp2;
	  }
	}
      } /* end loop over IF */
    } /* end loop over Pol */
  } /* end loop over baseline */
  
  /* Average time and interval for output */
  if (timeCount>0) *timec = timeSum/timeCount;
  else *timec = 0.0;
  *timei = (ltime - stime);
  
  /* DEBUG
  fprintf (stderr,"max accumIndex  %d\n", maxShit); */

  /* Cleanup */
  cleanup: if (blLookup) g_free(blLookup);
  if (err->error) Obit_traceback_val (err, routine, inUV->name, done);
  return done;
} /* end NextAvg */

/**
 *   Determine antenna gains.
 *   Does least squares solutions for phase and optionally 
 *   amplitude.  Three methods are available, "normal", "L1" and 
 *   "amplitude constrained". 
 *      All frequencies in each IF and at all times are assumed to have 
 *   been averaged.  If avgif is true then the data is assumed to have 
 *   been averaged in frequency and the solution found for the first IF 
 *   is copied to all IFs. 
 *   If avgpol  is true then the data is assumed to have been averaged in 
 *   polarization  and the data is copied to the second 
 * Routine translated from the AIPSish UVUTIL.FOR/SLFPA
 * \param vobs    Array of observed data [4,baseline,IF,pol]
 * \param ant1    Array of first antenna numbers (1-rel) on a baseline
 * \param ant2    Array of second antenna numbers (1-rel) on a baseline
 * \param numAnt  Maximum antenna number (NOT number of antennas)
 * \param numIF   Maximum number of IFs 
 * \param numPol  Maximum number of polarizations
 * \param numif   Number of IFs 
 * \param numpol  Number of Poln. in initial data (in SN)
 * \param refant  Reference antenna to use.
 * \param avgif   If true average in IF 
 * \param avgpol  If true average in polarization (Stokes 'I')
 * \param dol1    If true, use L1 solution
 * \param mode    Solution mode; 0= full gain, 1=phase  2=phase(ignore amp), 
 *                3=full, constrain amplitude.
 * \param minno   Minimum number of antannas allowed 
 * \param snrmin  Minimum SNR allowed.
 * \param prtlv   Print level, .ge. 4 gives some print.
 * \param creal   [out] (ant,if,pol) Real part of solution
 * \param cimag   [out] (ant,if,pol) Imag part of solution 
 * \param cwt     [out] (ant,if,pol) Weights = SNR
 * \param refan   [out] (if,pol,)     Reference antennas used
 * \param gotant  [out] (ant) If true corresponding antenna has data.
 * \param gain    Work array (2,ant)
 * \param snr     Work array (ant)
 * \param count   Work array, count of vis per antenna
 * \param err    Error/message stack, returns if error.
 */
static void 
doSolve (ofloat* vobs, olong *ant1, olong *ant2, olong numAnt, olong numIF, 
	 olong numPol, olong refant, gboolean avgif, 
	 gboolean avgpol, gboolean dol1, olong mode, olong minno, ofloat snrmin, 
	 olong prtlv, float* creal, ofloat* cimag, ofloat* cwt, olong* refan, 
	 gboolean* gotant, ofloat *gain, ofloat *snr, olong *count, ObitErr *err)
{
  olong iif, ist, iant, iBL, numBL, BLIndex, mPol, mIF;
  olong lprtlv, iref, ierr, lenEntry=4;
  ofloat amp, time=0.0, fblank =  ObitMagicF();
  /* No printout from CLBSNR */
  ofloat closer[2][2] = {{1.0e10,1.0e10},{1.0e10,1.0e10}};
  gchar prtsou[17] = "                ";
  
  /* Error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;  /* previous error? */
  
  /* Initialize output */
  /* (Loop) over Stokes' type */
  for (ist=0; ist<numPol; ist++) {
    /* Loop over IF */
    for (iif=0; iif<numIF; iif++) {
      refan[iif+ist*numIF] = 0;
      for (iant=0; iant<numAnt; iant++) {
	creal[iant+(iif+ist*numIF)*numAnt] = fblank;
	cimag[iant+(iif+ist*numIF)*numAnt] = fblank;
	cwt[iant+(iif+ist*numIF)*numAnt]   = 0.0;
      } /* end antenna loop */
    } /* end IF loop */
  } /* End stokes loop */
  /* No antennas found yet */
  for (iant=0; iant<numAnt; iant++) gotant[iant] = FALSE;
  
  numBL = (numAnt*(numAnt-1)) / 2;  /* Number of baselines */
  
  /* (Loop) over Stokes' type */
  BLIndex = 0;
  mPol = numPol;
  if (avgpol) mPol = 1;
  mIF  = numIF;
  if (avgif)  mIF  = 1;
  for (ist=0; ist<mPol; ist++) { /* loop 600 */
    /* Loop over IF */
    for (iif=0; iif<mIF; iif++) { /* loop 500 */
      /* Check if antennas have any data */
      for (iBL=0; iBL<numBL; iBL++) {
	if (vobs[BLIndex+iBL*lenEntry+2]>0.0) {
	  gotant[ant1[iBL]-1] = TRUE;
	  gotant[ant2[iBL]-1] = TRUE;	  
	}
      }
      iref = refant;
      
      /* Do solution */
      
      if (dol1) {/* L1 solution */
	gainCalcL1 (&vobs[BLIndex], ant1, ant2, numBL, numAnt,
		    refant, mode, minno, gain, &iref, prtlv, &ierr, err);
	
      } else { /* Normal */
	gainCalc (&vobs[BLIndex], ant1, ant2, numBL, numAnt,
		  refant, mode, minno, gain, &iref, prtlv, &ierr, err);
      } 
      
      if (ierr != 0) { /* Solution failed */
	for (iant=0; iant<numAnt; iant++) cwt[iant+(iif+ist*numIF)*numAnt] = -1.0;
	BLIndex += (+numBL)*lenEntry; continue;
      } 
      
      /* Convert amplitude to 1/amp * to correct data. */
      for (iant=0; iant<numAnt; iant++) { /* loop 110 */
	amp = gain[iant*2]*gain[iant*2] + gain[iant*2+1]*gain[iant*2+1];
	if (amp < 1.0e-20) amp = 1.0;
	gain[iant*2]   = gain[iant*2]   / amp;
	gain[iant*2+1] = gain[iant*2+1] / amp;
      } /* end loop  L110: */;
      
      /* Compute SNRs */
      lprtlv = prtlv;
      if (lprtlv > 3)  lprtlv += 3;
      calcSNR (&vobs[BLIndex], ant1, ant2, numBL, 
	       numAnt, gain, snr, closer, snrmin, time, iif, ist, count, 
	       lprtlv, prtsou, err);
      
      /* Save results */
      for (iant=0; iant<numAnt; iant++) { /* loop 150 */
	if (snr[iant] > snrmin) {
	  creal[iant+(iif+ist*numIF)*numAnt] = gain[iant*2];
	  cimag[iant+(iif+ist*numIF)*numAnt] = gain[iant*2+1];
	  cwt[iant+(iif+ist*numIF)*numAnt]   = snr[iant];
	} 
      } /* end loop  L150: */;
      refan[iif+ist*numIF] = iref;  /* Actual reference antenna */
      
      /* Averaging in poln? */
      if (avgpol) {
	for (iant=0; iant<numAnt; iant++) { /* loop 200 */
	  creal[iant+(iif+(ist+1)*numIF)*numAnt] = creal[iant+(iif+ist*numIF)*numAnt];
	  cimag[iant+(iif+(ist+1)*numIF)*numAnt] = cimag[iant+(iif+ist*numIF)*numAnt];
	  cwt[iant+(iif+(ist+1)*numIF)*numAnt]   = cwt[iant+(iif+ist*numIF)*numAnt];
	} /* end loop  L200: */;
	refan[iif+numIF] = refan[iif];
      } 
      BLIndex += (+numBL)*lenEntry;  /* Index into data array for this IF/Poln */
    } /* end loop over IF  L500: */;
  } /* end loop over Stokes' type L600: */;

  /* If averaging in IF copy soln. */
  if (avgif) {
    ist = 0;
    for (iif= 1; iif<numIF; iif++) { /* loop 620 */
      for (iant=0; iant<numAnt; iant++) { /* loop 610 */
	creal[iant+iif*numAnt] = creal[iant];
	cimag[iant+iif*numAnt] = cimag[iant];
	cwt[iant+iif*numAnt]   = cwt[iant];
	if (numPol>1) { /* Other poln */
	  ist = 1;
	  creal[iant+(iif+ist*numIF)*numAnt] = creal[iant+numIF*numAnt];
	  cimag[iant+(iif+ist*numIF)*numAnt] = cimag[iant+numIF*numAnt];
	  cwt[iant+(iif+ist*numIF)*numAnt]   = cwt[iant+numIF*numAnt];
	}
      } /* end loop  L610: */;
      refan[iif] = refan[0];
      if (numPol>1) refan[iif+numIF] = refan[numIF];
    } /* end loop  L620: */;
  } 
} /* end doSolve */

/**
 * calcSNR computes antenna based signal-to-noise ratioes (SNR) based  
 * on phase residuals from a gain model.  The approximation used is  
 * that SNR = 1.0 / RMS phase residual.  
 * If the values in closer are greater than 0 then any values exceeding 
 * these limits will be printed under control of prtlv  
 * Does a weighted solution for the SNR.  If insufficient data is  
 * present to compute an RMS but at least one observation exists, the  
 * SNR is set to 6.0.  
 * Routine translated from the AIPSish CLBSNR.FOR/CLBSNR
 * \param vobs    Normalized visibility [4,baseline,IF,pol]
 *                Zero value assumed invalid on input. 
 *                On return real part is replaced with the phase difference. 
 * \param ant1    Array of first antenna numbers. 
 * \param ant2    Array of second antenna numbers. 
 * \param numBL   Number of observations (baselines) 
 * \param numAnt  Number of antennas. 
 * \param gain    Antenna gains to be applied (real, imaginary) 
 * \param snr     [out] Signal to noise ratio for each antenna. 
 *                If an antenna is used in fewer than 3 baselines 
 *                the returned value is SNRMIN + 1.0 
 *                clipped at 1000.0
 * \param closer  (i,j) If (j=1/2) amplitude/phase closure errors 
 *                exceed these limits on average (i=1) or individually (i=2) 
 *                they will be printed 
 * \param snrmin  Minimum SNR allowed. 
 * \param time    Time in days of data, used for labeling. 
 * \param iif     IF number used for labeling only.
 * \param ist     Stokes parameter of soln. 1-5 => R,L,R,L,I. 
 * \param count   A work array used for the counts for each 
 *                antenna, must be at least MAXANT in size. 
 * \param prtlv   Print level: 0 none, 2 statistics of failures, 
 *                4 individual failures, 5 the antenna SNRs 
 * \param prtsou  Current source name. 
 * \param err     Error/message stack, returns if error.
 */
static void 
calcSNR (ofloat* vobs, olong *ant1, olong *ant2, olong numBL, 
	 olong numAnt, ofloat* gain, ofloat* snr, ofloat closer[2][2], 
	 ofloat snrmin, odouble time, olong iif, olong ist, olong* count, 
	 olong prtlv, gchar* prtsou, ObitErr *err) 
{
  gchar pol[5][5] = {"Rpol","Lpol","Rpol","Lpol","Ipol"};
  olong   loop, ii, jj, nprt, id, ih, im, ls, blprt[3][3], ne,lenEntry=4;
  gboolean   doclos, msgdun, docls1, docls2;
  ofloat      zr, zi, zzr, zzi, prtsnr, tmtemp, phsq, ae, pe, blrprt[3];
  ofloat *sumwt=NULL;
  olong   *error=NULL;
  gchar msgtxt[81];
  
  /* Error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;  /* previous error? */
  
  /* What messages are wanted? */
  docls1 = (closer[0][0] * closer[1][0] > 1.0e-20)  &&  
    (closer[0][0] * closer[1][0] < 1.0e20)  &&  (prtlv > 3);
  docls2 = ((closer[0][1] * closer[1][1] > 1.0e-20)  &&  
	    (closer[0][1] * closer[1][1] < 1.0e20))  &&  (prtlv > 3);
  doclos = (docls1)  ||  (docls2);
  msgdun = FALSE;
  
  /* Label for any closure errors: */
  if (doclos) {
    id = time;
    tmtemp = (time - id) * 24.0;
    ih = tmtemp;
    tmtemp = (tmtemp - ih) * 60.0;
    im = tmtemp;
    tmtemp = (tmtemp - im) * 60.0;
    ls = tmtemp;
    if (ls == 60) {
      ls = 0;
      im = im + 1;
      if (im == 60)  {
	im = 0;
	ih = ih + 1;
	if (ih == 24) {
	  ih = 0;
	  id = id + 1;
	} 
      } 
    } 
    g_snprintf (msgtxt,80,"closure errors at %3d/%2.2d:%2.2d:%2.2d %s IF %3d %s", 
		id, ih, im, ls, prtsou, iif, pol[ist-1]);
  } 

  /* Create work arrays */
  sumwt = g_malloc0(numAnt*sizeof(ofloat));
  error = g_malloc0(numAnt*sizeof(olong));
 
  /* Zero sums, counts etc. */
  nprt = 0;
  for (loop=0; loop<numAnt; loop++) { /* loop 10 */
    /* ?? snr[loop]   = 1.0e-6; */
    snr[loop]   = 0.0;
    sumwt[loop] = 0.0;
    count[loop] = 0;
    error[loop] = 0;
  } /* end loop  L10:  */
  
  /* Determine phase residuals. */
  ne = 0;
  ae = 0.0;
  pe = 0.0;
  for (loop=0; loop<numBL; loop++) { /* loop 30 */
    if (vobs[loop*lenEntry+2]>1.0e-20) {
      ii = ant1[loop]-1;
      jj = ant2[loop]-1;
      zr = gain[ii*2] * gain[jj*2] + gain[ii*2+1] * gain[jj*2+1];
      zi = gain[ii*2] * gain[jj*2+1] - gain[ii*2+1] * gain[jj*2];
      zzr = vobs[loop*lenEntry] * zr - vobs[loop*lenEntry+1] * zi;
      zzi = vobs[loop*lenEntry] * zi + vobs[loop*lenEntry+1] * zr;
      vobs[loop*lenEntry+1] = sqrt (zzr*zzr + zzi*zzi);
      if (vobs[loop*lenEntry+1]*vobs[loop*lenEntry+2] > 1.0e-20) {
	vobs[loop*lenEntry] = atan2 (zzi, zzr);
	if (docls1) {
	  ne = ne + 1;
	  pe = pe + fabs (vobs[loop*lenEntry]);
	  ae = ae + fabs (log10 (fabs (vobs[loop*lenEntry+1]+1.e-20)));
	} 
      } 
    }
  } /* end loop  L30:  */
  
  /* Statistical failure */
  if ((docls1)  &&  (ne > 0)) {
    pe = pe / ne;
    ae = ae / ne;
    ae = pow (10.0, ae) - 1.0;
    /* print header, message */
    if ((ae > closer[0][0])  ||  (pe > closer[1][0])) {
      Obit_log_error(err, OBIT_InfoErr, "%s", msgtxt);
      msgdun = TRUE;
      ae = ae * 100.0;
      pe = pe * 57.296;
      Obit_log_error(err, OBIT_InfoErr, "Average closure error %10.3f %8.2f d", 
		     ae, pe);
    } 
  } 
  
  /* Sum square residuals */
  ne = 0;
  for (loop=0; loop<numBL; loop++) { /* loop 50 */
    if (vobs[loop*lenEntry+1]*vobs[loop*lenEntry+2] >= 1.0e-20) {
      phsq = vobs[loop*lenEntry] * vobs[loop*lenEntry] * vobs[loop*lenEntry+2];
      ii = ant1[loop]-1;
      count[ii]++;
      snr[ii] += phsq;
      sumwt[ii] += vobs[loop*lenEntry+2];
      jj = ant2[loop]-1;
      count[jj]++;
      snr[jj] += phsq;
      sumwt[jj] += vobs[loop*lenEntry+2];
      
      /* check closure error */
      if (docls2) {
	pe = fabs (vobs[loop*lenEntry]);
	ae = vobs[loop*lenEntry+1];
	if ((ae < 1.0)  &&  (ae > 0.0)) ae = 1.0 / ae;
	ae = ae - 1.0;
	if ((ae > closer[0][1])  ||  (pe > closer[1][1])) {
	  error[ii]++;
	  error[jj]++;
	  ne = ne + 1;
	  /* individual messages */
	  if (prtlv >= 5) {
	    /* Print header message */
	    if (!msgdun) {
              Obit_log_error(err, OBIT_InfoErr, "%s", msgtxt);
	      msgdun = TRUE;
	    } 
	    /* Flush buffer if full */
	    if (nprt >= 5) {
              Obit_log_error(err, OBIT_InfoErr, 
		   "%4.2d-%2.2d %7.1f %5d %4.2d-%2.2d %7.1f %5d %4.2d-%2.2d %7.1f %5d ", 
			     blprt[0][0], blprt[0][1],  blrprt[0], blprt[0][2],
			     blprt[1][0], blprt[1][1],  blrprt[1], blprt[1][2],
			     blprt[2][0], blprt[2][1],  blrprt[1], blprt[3][2]);
	      nprt = 0;
	    } 
	    /* New entry */
	    blprt[nprt][0] = ii;
	    blprt[nprt][1] = jj;
	    blprt[nprt][2] = pe * 57.296 + 0.5;
	    blrprt[nprt] = MAX (-9999., MIN (9999., ae*100.));
	    nprt = nprt + 1;
	  } 
	} 
      } 
    } 
  } /* end loop  L50:  */;
  
  /* Flush closure message buffer */
  if (nprt > 0) {
    Obit_log_error(err, OBIT_InfoErr, 
		   "%4.2d-%2.2d %7.1f %5d %4.2d-%2.2d %7.1f %5d %4.2d-%2.2d %7.1f %5d ", 
		   blprt[0][0], blprt[0][1],  blrprt[0], blprt[0][2],
		   blprt[1][0], blprt[1][1],  blrprt[1], blprt[1][2],
		   blprt[2][0], blprt[2][1],  blrprt[1], blprt[3][2]);
    nprt = 0;
  } 
  
  /* summary by antenna */
  if (ne > 0) {
    if (!msgdun) {
      Obit_log_error(err, OBIT_InfoErr, "%s", msgtxt);
      msgdun = TRUE;
    } 
    for (loop=0; loop<numAnt; loop++) { /* loop 60 */
      if (error[loop] > 0) {
	Obit_log_error(err, OBIT_InfoErr, "Antenna %2d had %6d excess closure errors", 
		       loop+1, error[loop]);
      } 
    } /* end loop  L60:  */;
  } 

  /* Convert to SNRs */
  for (loop=0; loop<numAnt; loop++) { /* loop 70 */
    if (count[loop] <= 0) snr[loop] = 0.0;
    /* For insufficient data (but  some) use SNRMIN + 1.0 */
    if ((count[loop] >= 1)  &&  (count[loop] <= 2)) snr[loop] = snrmin + 1.0;
    if ((count[loop] >= 3)  &&  (snr[loop] > 1.0e-20)) {
      snr[loop] = 1.0 / sqrt (snr[loop] / ((count[loop]-1) * sumwt[loop]));
      snr[loop] = MIN (1000.0, snr[loop]);
      /* DEBUG Solution failed
	 if (snr[loop]< snrmin) {
	 Obit_log_error(err, OBIT_InfoErr, "Antenna %2d had %6ld counts %f snr", 
	 loop+1, count[loop], snr[loop]);
	 } */
    }
    
    /* Print result if desired. */
    if ((count[loop] >= 1)  &&  (prtlv >= 5)) {
      prtsnr = MIN (9999.999, snr[loop]);
      Obit_log_error(err, OBIT_InfoErr, "antenna(%2d)  %3d obs, snr = %10.3f", 
		     loop+1, count[loop], prtsnr);
    }
  } /* end loop  L70:  */;

  /* Cleanup */
  if (sumwt) g_free(sumwt);
  if (error) g_free(error);
} /* end of routine calcSNR */ 

/**
 * Computes antenna gains given visibilities divided by the  
 * model values.  
 * Routine translated from the AIPSish GCALC.FOR/GCALC  
 * \param vobs    Complex normalized visibility. 
 * \param is      Array of first antenna numbers. 
 * \param js      Array of second antenna numbers. 
 * \param wt      Array of visibility weights. 
 * \param numBL   Number of observations. 
 * \param numAnt  Maximum antenna number (NOT number of antennas)
 * \param refant  Desired reference antenna. 
 * \param mode    Solution mode: 0 = full gain soln. 
 *                1 = phase only keep ampl. info. 
 *                2 = phase only discard ampl. info. 
 * \param minno   Min. number of antennas acceptable. 
 * \param g       Complex antenna gains to be applied. 
 * \param nref    Reference antenna used. 
 * \param prtlv   Print flag,    0=none, 4=soln, 5=data plus soln
 * \param ierr    Return error code 0 = OK, 1 = no valid data, 2=didn't converge, 
 *                3 = too few antennas
 * \param err     Error/message stack, returns if error.
 */
static void 
gainCalc (ofloat* vobs, olong *ant1, olong *ant2, olong numBL, olong numAnt,
	  olong refant, olong mode, olong minno, ofloat *g, olong* nref,
	  olong prtlv, olong *ierr, ObitErr* err) 
{
  ofloat gng[2], z[2], zr, zi, amp, gd, ph, qq, rms, s, sumwt, tol, w=0.0, x, xx, xxx, 
    yy, wx, arg, xfail, xfail1, xfail2;
  olong  i, ii, it, j, jj, k, nt, ntd, ntm1, itmax, nfail, iref, lenEntry=4;
  gboolean   convgd;
  ofloat *gn=NULL, *glast=NULL, *swt=NULL, tempR, tempI;
  gchar  msgtxt[81];

  /* Error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return ;  /* previous error? */

  itmax = 60;
  tol = 5.0e-6;
  nt = 0;
  ntd = 0;
  sumwt = 0.0;
  *ierr = 0;

  /* Any data? */
  if (numBL <= 0) {*ierr = 1; return; } 

  /* Create work arrays */
  gn    = g_malloc0(2*numAnt*sizeof(ofloat));
  glast = g_malloc0(2*numAnt*sizeof(ofloat));
  swt   = g_malloc0(numAnt*sizeof(ofloat));
  for (i=0; i<numAnt; i++) swt[i] = 0.0;

  /* Print data if nec. */
  if (prtlv >= 5) {
    for (k=0; k<numBL; k++) { /* loop 20 */
      if (vobs[k*lenEntry+2] > 1.0e-20) {
	amp = sqrt (vobs[k*lenEntry]*vobs[k*lenEntry] + 
		    vobs[k*lenEntry+1]*vobs[k*lenEntry+1]);
	ph = 57.296 * atan2 (vobs[k*lenEntry+1], vobs[k*lenEntry]+1.0e-20);
	g_snprintf (msgtxt,80," %4d amp =%12.3e phase =%9.2f  %3d  %3d %12.3e", 
		    k+1, amp, ph, ant1[k], ant2[k], vobs[k*lenEntry+2]);
	Obit_log_error(err, OBIT_InfoErr, msgtxt);
      } 
    } /* end loop  L20:  */;
  } /* end if print */
 
  /* Find which antennas have data  and normalize to unit amplitude if requested. */
  for (k=0; k<numBL; k++) { /* loop 40 */
    i = ant1[k]-1;
    j = ant2[k]-1;
    nt = MAX (nt, ant2[k]);
    xxx = sqrt ((vobs[k*lenEntry]*vobs[k*lenEntry]) + 
		(vobs[k*lenEntry+1]*vobs[k*lenEntry+1]));
    if (xxx <= 1.0e-20) vobs[k*lenEntry+2] = 0.0;
    if (vobs[k*lenEntry+2] > 1.0e-20) {
      swt[i] += vobs[k*lenEntry+2];
      swt[j] += vobs[k*lenEntry+2];
      sumwt  += vobs[k*lenEntry+2];
      if (mode == 2) {  /* normalize visibility */
	vobs[k*lenEntry]   /=  xxx;
	vobs[k*lenEntry+1] /=  xxx;
      } 
    } 
  } /* end loop  L40:  */;

  /* Initialize solutions */
  for (i=0; i<numAnt; i++) { /* loop 50 */
    g[i*2]       = 1.0;
    g[i*2+1]     = 0.0;
    glast[i*2]   = 1.0;
    glast[i*2+1] = 0.0;
    if (swt[i] > 1.0e-20) ntd++;
  } /* end loop  L50:  */;

  /* Any data? */
  if (sumwt <= 1.0e-20) {
    *ierr = 1; 
    goto cleanup;
  } 

  /* Pick reference ant. if the input one not usable */
  (*nref) = refant;
  if (((*nref < 1)  ||  (*nref > numAnt)  ||  (swt[(*nref)-1] == 0.0))) { /* goto L100;*/
    /* Find antenna with highest summed weight */
    xx = 0.0;
    for (i=0; i<nt; i++) { /* loop 90 */
      if (swt[i] > xx) { /*goto L90;*/
	xx = swt[i];
	*nref = i+1;
      }
    } /* end loop  L90:  */;
  } /* L100: */

  ntm1 = nt - 1;
  /* Too few antennas? */
  if (ntd < minno) {
    *ierr = 3; 
    goto cleanup;
  } 

  /* Print statistics */
  if (prtlv >= 5) {
    /* Sum chi squares */
    s = 0.0;
    for (k=0; k<numBL; k++) { /* loop 160 */
      if (vobs[k*lenEntry+2]>0.0) {
	ii = ant1[k];
	jj = ant2[k];
	z[0] = g[jj*2];
	z[1] = - g[jj*2+1];
	zr = g[ii*2] * z[0] - g[ii*2+1] * z[1];
	zi = g[ii*2] * z[1] + g[ii*2+1] * z[0];
	z[0] = vobs[k*lenEntry]   - zr;
	z[1] = vobs[k*lenEntry+1] - zi;
	qq = z[0] * z[0] + z[1] * z[1];
	s = s + vobs[k*lenEntry+2] * qq;
      }
    } /* end loop  L160: */;
    rms = sqrt (s/sumwt);
    it = 0;
    g_snprintf (msgtxt, 80, "iter= %5d s=%15.5e rms=%15.5e", it, s, rms);
    Obit_log_error(err, OBIT_InfoErr, msgtxt);
  } /* end print */
 
  /* Begin solution iteration */
  for (it=1; it<=itmax; it++) { /* loop 310 */
    if (it > 15) tol = (it - 10) * 1.0e-6;

    /* Full amp and phase solution */
    if (mode <= 0) { /* goto L230;*/
      /* Following for amplitude and phase solution. */
      if (it == 1) w = 0.5;
      if (it == 2) w = 0.75;
      if (it > 2)  w = 0.9;
      if (it > 10) w = 0.5;
      if (ntd <= 6) w = 0.25;
      for (i=0; i<nt; i++) { /* loop over antennas 220 */
	if (swt[i] == 0.0) continue; /* no valid data?  goto L220;*/

	gng[0] = 0.0;  /* zero accumulators */
	gng[1] = 0.0;
	gd = 0.0;

	/* Loop over data looking for antenna i */
	for (k=0; k<numBL; k++) {
	  if (vobs[k*lenEntry+2]<0.0) continue;
	  if ((i+1) == ant1[k]) {
	    /* i as first antenna */
	    jj = ant2[k]-1;
	    z[0] = g[jj*2];
	    z[1] = g[jj*2+1];
	    qq = z[0] * z[0] + z[1] * z[1];
	    gd     += vobs[k*lenEntry+2] * qq;
	    gng[0] += vobs[k*lenEntry+2] * (g[jj*2] * vobs[k*lenEntry]   - g[jj*2+1] * vobs[k*lenEntry+1]);
	    gng[1] += vobs[k*lenEntry+2] * (g[jj*2] * vobs[k*lenEntry+1] + g[jj*2+1] * vobs[k*lenEntry]);

	    /* i as second antenna: */
	  } else if ((i+1) == ant2[k]) {
	    ii = ant1[k]-1;
	    z[0] = g[ii*2];
	    z[1] = g[ii*2+1];
	    qq = z[0] * z[0] + z[1] * z[1];
	    gd     += vobs[k*lenEntry+2] * qq;
	    gng[0] += vobs[k*lenEntry+2] * (g[ii*2]   * vobs[k*lenEntry] + g[ii*2+1] * vobs[k*lenEntry+1]);
	    gng[1] += vobs[k*lenEntry+2] * (g[ii*2+1] * vobs[k*lenEntry] - g[ii*2]   * vobs[k*lenEntry+1]);
	  }
	} /* end loop of data */

	/* corrections to gain for antenna i */
	g[i*2]   += w * (gng[0] / gd - g[i*2]);
	g[i*2+1] += w * (gng[1] / gd - g[i*2+1]);
      } /* end loop over data  L220: */

      /* End amplitude and phase */

    } else { /* Phase only L230 */
      /* Following for phase only. */
      for (i=0; i<nt; i++) { /* loop 240 */
	gn[i*2]   = 0.0;
	gn[i*2+1] = 0.0;
      } /* end loop  L240: */;

      /* Damping factor */
      w = .8;
      if (ntd <= 6) w = .25;
      for (k=0; k<numBL; k++) { /* loop over data 250 */
	if (vobs[k*lenEntry+2] > 0.0) {
	  ii = ant1[k]-1;
	  jj = ant2[k]-1;
	  /*               Z = G(II) * CONJG (G(JJ)) * CONJG (WORK(K)) */
	  zr = (g[ii*2]   * g[jj*2] + g[ii*2+1] * g[jj*2+1]);
	  zi = (g[ii*2+1] * g[jj*2] - g[ii*2]   * g[jj*2+1]);
	  /*               GN(II) = GN(II) + Z */
	  tempR = vobs[k*lenEntry+2] * (zr * vobs[k*lenEntry] + zi * vobs[k*lenEntry+1]);
	  tempI = vobs[k*lenEntry+2] * (zi * vobs[k*lenEntry] - zr * vobs[k*lenEntry+1]);
	  gn[ii*2]   += tempR;
	  gn[ii*2+1] += tempI;
	  /*               GN(JJ) = GN(JJ) + CONJG (Z) */
	  gn[jj*2]   += tempR;
	  gn[jj*2+1] -= tempI;
	} /* end of if valid */
      } /* end loop  L250: */

      /* Update solutions */
      for (i=0; i<nt; i++) { /* loop over antennas 260 */
	wx = w;
	if (swt[i] <= 0.0) wx = 0.0;
	/*               XX = ATAN2 (AIMAG (G(I)), REAL (G(I))) */
	xx = atan2 (g[i*2+1], g[i*2]);
	/*               YY = ATAN2 (AIMAG (GN(I)), REAL (GN(I))) */
	yy = atan2 (gn[i*2+1], gn[i*2]+1.0e-30);
	arg = xx - wx * yy;
	g[i*2]   = cos (arg);
	g[i*2+1] = sin (arg);
      } /* end loop  L260: */;
    } /* end phase only */

    /* Convergence test */
    convgd = TRUE; /* L265: */
    nfail = 0;
    xfail = 0.0;
    for (i=0; i<nt; i++) { /* loop 280 */
      if ((swt[i] > 0.0) && convgd) { /*goto L275;*/
      /* Use different criteria before and after  30 iterations */
	if (it <= 30) { /* before 30 goto L270;*/
	  x = tol * (5.0e-7 + sqrt (g[i*2]*g[i*2] + g[i*2+1]*g[i*2+1]));
	  if (fabs (g[i*2]  -glast[i*2])   > x) convgd = FALSE;
	  if (fabs (g[i*2+1]-glast[i*2+1]) > x) convgd = FALSE;

	} else { /* after 30  L270: */
	  xfail1 = fabs (g[i*2]   - glast[i*2]);
	  xfail2 = fabs (g[i*2+1] - glast[i*2+1]);
	  xfail1 = MAX (xfail1, xfail2);
	  x = tol * (5.0e-7 + sqrt (g[i*2]*g[i*2]+g[i*2+1]*g[i*2+1]));
	  if (xfail1 > x) { /*goto L275;*/
	    convgd = FALSE;
	    nfail = nfail + 1;
	    if ((xfail1/x) >= xfail) { /* goto L275; */
	      xfail = xfail1 / x;
	    }
	  }
	}
      } /* end if got data and still "converged" L275: */

      /* Save as last solution */
      glast[i*2]   = g[i*2];
      glast[i*2+1] = g[i*2+1];
    } /* end loop  L280: */;

    if (prtlv >= 5) {
      /* Print statistics */
      s = 0.0;
      for (k=0; k<numBL; k++) { /* loop 290 */
	if (vobs[k*lenEntry+2]>0) {
	  ii = ant1[k]-1;
	  jj = ant2[k]-1;
	  z[0] = g[jj*2];
	  z[1] = - g[jj*2+1];
	  zr = g[ii*2] * z[0] - g[ii*2+1] * z[1];
	  zi = g[ii*2] * z[1] + g[ii*2+1] * z[0];
	  z[0] = vobs[k*lenEntry] - zr;
	  z[1] = vobs[k*lenEntry+1] - zi;
	  qq = z[0] * z[0] + z[1] * z[1];
	  s += vobs[k*lenEntry+2] * qq;
	}
      } /* end loop  L290: */;
      rms = sqrt (s/sumwt);
      g_snprintf (msgtxt,80,"iter= %5d s=%15.5e rms=%15.5e", it, s, rms);
      Obit_log_error(err, OBIT_InfoErr, msgtxt);
    } /* end print statistics */ 

    if (convgd) break;   /* Converged?  goto L400;*/
    if ((it > 30) && (xfail <= 10.0)) break;   /* Close enough goto L400; */
    if (it >= itmax) { /* goto L310; */
      /* Didn't converge */
      *ierr = 2;
      goto cleanup;
    }
  }    /* End of iteration loop */

  /* Refer solutions to reference antenna. */
  iref = *nref-1;
  xxx = 1.0 / sqrt (g[iref*2]*g[iref*2] + g[iref*2+1]*g[iref*2+1]);
  for (i=0; i<nt; i++) { /* loop 600 */
    if ((i != iref)  &&  (swt[i] > 1.0e-20)) {
      /*            G(I) = G(I) / G(NREF) * CABS (G(NREF)) */
      zr = g[i*2]   * g[iref*2] + g[i*2+1] * g[iref*2+1];
      zi = g[i*2+1] * g[iref*2] - g[i*2]   * g[iref*2+1];
      g[i*2]   = xxx * zr;
      g[i*2+1] = xxx * zi;
    } 
  } /* end loop  L600: */;

  /* adjust reference antenna gain */
  g[iref*2] = sqrt (g[iref*2]*g[iref*2] + g[iref*2+1]*g[iref*2+1]);
  g[iref*2+1] = 0.0;

  if (prtlv >= 4) {
    for (i=0; i<nt; i++) { /* loop 610 */
      if (swt[i] > 1.0e-20) {
	/* Print results. */
	amp = sqrt (g[i*2]*g[i*2] + g[i*2+1]*g[i*2+1]);
	ph = 57.2958 * atan2 (g[i*2+1], g[i*2]);
	g_snprintf (msgtxt,80,"ant=  %5d amp=%12.5f phase=%12.2f", 
		    i+1, amp, ph);
	Obit_log_error(err, OBIT_InfoErr, msgtxt);
      } 
    } /* end loop  L610: */;
  } /* end of print */ 

  *ierr = 0; /* ok */
  /* Cleanup */
  cleanup: if (gn)    g_free(gn); 
  if (glast) g_free(glast);
  if (swt)   g_free(swt);
} /* end of routine gainCalc */ 


/**
 * Does an L1 type gain solution similar to gainCalc.  
 * Routine translated from the AIPSish GCALC1.FOR/GCALC1  
 * \param vobs    Complex normalized visibility. 
 * \param is      Array of first antenna numbers. 
 * \param js      Array of second antenna numbers. 
 * \param wt      Array of visibility weights. 
 * \param numBL   Number of observations. 
 * \param numAnt  Maximum antenna number (NOT number of antennas)
 * \param refant  Desired reference antenna. 
 * \param mode    Solution mode: 0 = full gain soln. 
 *                1 = phase only keep ampl. info. 
 *                2 = phase only discard ampl. info. 
 * \param minno   Min. number of antennas acceptable. 
 * \param g       Complex antenna gains to be applied. 
 * \param nref    Reference antenna used. 
 * \param prtlv   Print flag,    0=none, 4=soln, 5=data plus soln
 * \param ierr    Return error code 0 = OK, 1 = no valid data, 3 = too few antennas
 * \param err     Error/message stack, returns if error.
 */
static void 
gainCalcL1 (ofloat* vobs, olong *ant1, olong *ant2, olong numBL, 
	    olong numAnt, olong refant, olong mode, olong minno, ofloat *g, olong* nref, 
	    olong prtlv, olong *ierr, ObitErr* err) 
{
  olong   k, nt, ntd, ie, ne, it, itmax, i, ii, j, jj, iref, lenEntry=4;
  ofloat  amp, ph, xlamb, pterm;
  odouble gd, t, eps, f, tol, s, qq, w=0.0, xx, yy, x, xrr, xrd, xrz, xrtemp, xrgng, 
    xir, xid, xiz, xitemp, xigng, d1, d2, f1, f2, rmean, t1, t2, sumwt, xxx;
  gboolean   convgd;
  ofloat  epss[3] = {5.0e-3, 5.0e-4, 5.0e-5};
  ofloat  *swt=NULL, *p1=NULL, *p2=NULL;
  odouble *xrdg=NULL, *xidg=NULL, *xrglas=NULL, *xiglas=NULL;
  gchar msgtxt[81];

  /* Error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;  /* previous error? */

  /* Anything to do? */
  if (numBL <= 0) {*ierr=1; return;}

  /* Create work arrays */
  swt    = g_malloc0(numAnt*sizeof(ofloat));
  p1     = g_malloc0(numAnt*sizeof(ofloat));
  p2     = g_malloc0(numAnt*sizeof(ofloat));
  xrdg   = g_malloc0(numAnt*sizeof(odouble));
  xidg   = g_malloc0(numAnt*sizeof(odouble));
  xrglas = g_malloc0(numAnt*sizeof(odouble));
  xiglas = g_malloc0(numAnt*sizeof(odouble));

  /* Init values */
  xlamb = 0.0;
  pterm = 0.0;
  ne = 3;
  itmax = 60;
  tol = 5.0e-5;
  nt = 0;
  ntd = 0;
  sumwt = 0.0;
  for (i= 1; i<=numAnt; i++) { /* loop 10 */
    swt[i-1] = 0.;
  } /* end loop  L10:  */;

  /* Find which antennas have data and normalize to unit amplitude if requested. */
  for (k=0; k<numBL; k++) { /* loop 40 */
    /* Data dump if requested */
    if ((prtlv >= 5)  &&  (vobs[k*lenEntry+2] > 1.0e-20)) {
      amp = sqrt (vobs[k*lenEntry]*vobs[k*lenEntry] + vobs[k*lenEntry+1]*vobs[k*lenEntry+1]);
      ph = 57.296 * atan2 (vobs[k*lenEntry+1], vobs[k*lenEntry]+1.0e-20);
      /* DEBUG g_snprintf (msgtxt,80," %4d amp =%12.3e phase =%9.3f  %3d  %3d %12.3e", 
		  k+1, amp, ph, ant1[k], ant2[k], vobs[k*lenEntry+2]); */
      
      g_snprintf (msgtxt,80," %4d real =%12.5e imag =%12.5e  %3d  %3d %12.5e", 
		  k+1, vobs[k*lenEntry], vobs[k*lenEntry+1], ant1[k], ant2[k], vobs[k*lenEntry+2]);
      Obit_log_error(err, OBIT_InfoErr, msgtxt);
    } 
    if (vobs[k*lenEntry+2] > 0.0) {
      i = ant1[k]-1;
      j = ant2[k]-1;
      nt = MAX (nt, ant2[k]);
      xxx = sqrt ((vobs[k*lenEntry]*vobs[k*lenEntry]) + (vobs[k*lenEntry+1]*vobs[k*lenEntry+1]));
      if (xxx > 1.0e-20) {
	swt[i] += vobs[k*lenEntry+2];
	swt[j] += vobs[k*lenEntry+2];
	sumwt  += vobs[k*lenEntry+2];
	/* Normalize by amplitude if only solving for phase without ampl. */
	if (mode == 2) {
	  vobs[k*lenEntry]   /= xxx;
	  vobs[k*lenEntry+1] /= xxx;
	} 
      }
    } 
  } /* end loop  L40:  */

  /* Initialize solutions */
  for (i=0; i<numAnt; i++) { /* loop 50 */
    g[i*2]    = 1.0;
    g[i*2+1]  = 0.0;
    xrdg[i]   = 1.0;
    xidg[i]   = 0.0;
    xrglas[i] = 1.0;
    xiglas[i] = 0.0;
    if (swt[i] > 0.0) ntd++;
  } /* end loop  L50:  */;

  /* Any data found? */
  if (sumwt == 0.0) {
    /* No data */
    *ierr = 1;
    goto cleanup;
  }

  /* Pick reference ant. if the input one not usable */
  (*nref) = refant;
  if (((*nref < 1)  ||  (*nref > numAnt)  ||  (swt[(*nref)-1] <= 0.0))) { /* goto L100;*/
    /* no good, find antenna with highest summed weight */
    xx = 0.0;
    for (i=0; i<nt; i++) { /* loop 90 */
      if (swt[i] > xx) { /*goto L90;*/
	xx = swt[i];
	*nref = i+1;
      }
    } /* end loop  L90:  */;
  } /* L100: */

  if (ntd < minno) { /* goto L930;*/
    /* Too few antennas */
    *ierr = 3;
    goto cleanup;
  }

  /* L1 outer loop */
  for (ie= 1; ie<=ne; ie++) { /* loop 320 */
    eps = epss[ie-1];
    /* Print rms residual if requested */
    if (prtlv >= 5) { /*goto L170;*/
      s = 0.0;
      for (k=0; k<numBL; k++) { /* loop 160 */
	if (vobs[k*lenEntry+2] > 0.0) { /*goto L160;*/
	  i = ant1[k]-1;
	  j = ant2[k]-1;
	  xrz = vobs[k*lenEntry]   - (xrdg[i] * xrdg[j] + xidg[i] * xidg[j]);
	  xiz = vobs[k*lenEntry+1] - (xrdg[j] * xidg[i] - xrdg[i] * xidg[j]);
	  qq = xrz * xrz + xiz * xiz;
	  s += vobs[k*lenEntry+2] * sqrt (qq + eps);
	  /* DEBUG */
	  Obit_log_error(err, OBIT_InfoErr, " %5d  %5d  %5d %15.5g %15.5g %15.5g %15.5g",
			 k+1,i+1,j+1,xrz,xiz,qq,s);
	} /* end of if data valid */
      } /* end loop  L160: */;
      rmean = s / sumwt;
      it = 0;
      g_snprintf (msgtxt,80,"iter= %5d s=%15.5e rmean=%15.5e", it, s, rmean);
      Obit_log_error(err, OBIT_InfoErr, msgtxt);
    } /* end of if print */

    /* Inner Solution loop */
    for (it= 1; it<=itmax; it++) { /* loop  */
      if (it > 15) tol = (it-10) * 1.0e-5;
      if (mode<=0) {
	/* Section to solve for full  complex gains: */
	if (it == 1) w = 0.5;
	if (it == 2) w = 1.5;
	if (it == 3) w = 1.5;
	if (ntd <= 6) w = .25;
	for (i=0; i<nt; i++) { /* loop 220 */
	  if (swt[i] > 0.0) { /*goto L220;*/
	    xrgng = 0.0;
	    xigng = 0.0;
	    gd = 0.0;
	    for (k=0; k<numBL; k++) { /* loop 210 */
	      if (vobs[k*lenEntry+2] > 0.0) { /*goto L210; */
		ii = ant1[k]-1;
		jj = ant2[k]-1;
		if ((swt[ii] == 0.) || (swt[jj] == 0.)) continue; /* goto L210;*/
		if (i == ii) {
		  /* i as first antenna */
		  xrr = vobs[k*lenEntry]   - (xrdg[ii] * xrdg[jj] + xidg[ii] * xidg[jj]);
		  xir = vobs[k*lenEntry+1] - (xrdg[jj] * xidg[ii] - xrdg[ii] * xidg[jj]);
		  t = vobs[k*lenEntry+2] / sqrt (xrr * xrr + xir * xir + eps);
		  xrz = xrdg[jj];
		  xiz = xidg[jj];
		  gd    += t * (xrz * xrz + xiz * xiz);
		  xrgng += t * (xrdg[jj] * vobs[k*lenEntry]   - xidg[jj] * vobs[k*lenEntry+1]);
		  xigng += t * (xrdg[jj] * vobs[k*lenEntry+1] + vobs[k*lenEntry] * xidg[jj]);

		} else if (i == jj) {
		  /* i as second antenna: */
		  xrr = vobs[k*lenEntry]   - (xrdg[ii] * xrdg[jj] + xidg[ii] * xidg[jj]);
		  xir = vobs[k*lenEntry+1] - (xrdg[jj] * xidg[ii] - xrdg[ii] * xidg[jj]);
		  t = vobs[k*lenEntry+2] / sqrt (xrr * xrr + xir * xir + eps);
		  xrz = xrdg[ii];
		  xiz = xidg[ii];
		  gd    += t * (xrz * xrz + xiz * xiz);
		  xrgng += t * (xrdg[ii]  * vobs[k*lenEntry] + xidg[ii] * vobs[k*lenEntry+1]);
		  xigng += t * (vobs[k*lenEntry] * xidg[ii]  - xrdg[ii] * vobs[k*lenEntry+1]);
		} /* end i as second antenna */
	      } /* end if baseline valid data */
	    } /* end loop  L210: */;

	    xrdg[i] += w * (xrgng / gd - xrdg[i]);
	    xidg[i] += w * (xigng / gd - xidg[i]);
	  } /* end if antenna has valid data*/
	} /* end loop  L220: */;
	
      } else {
	/* Phase only solutions */
	for (i=0; i<nt; i++) { /* loop  240*/
	  p1[i] = 0.;
	  p2[i] = 0.;
	} /* end loop  L240: */;
	
	w = 0.5;        /* damping */
	if (it <= 3)  w = 0.25;
	if (ntd <= 6) w = 0.25;
	for (k=0; k<numBL; k++) { /* loop 250 */
	  if (vobs[k*lenEntry+2]>0.0) {
	    i = ant1[k]-1;
	    j = ant2[k]-1;
	    xrd = xrdg[i] * xrdg[j] + xidg[i] * xidg[j];
	    xid = xrdg[j] * xidg[i] - xrdg[i] * xidg[j];
	    xrr = vobs[k*lenEntry]   - xrd;
	    xir = vobs[k*lenEntry+1] - xid;
	    xrtemp = xrd;
	    xitemp = xid;
	    xrd = xrtemp    * vobs[k*lenEntry] + xitemp * vobs[k*lenEntry+1];
	    xid = vobs[k*lenEntry] * xitemp    - xrtemp * vobs[k*lenEntry+1];
	    qq = xrr * xrr + xir * xir;
	    f = qq + eps;
	    f1 = sqrt (f);
	    f2 = f * f1;
	    d1 = 2.0 * xid;
	    d2 = 2.0 * xrd;
	    t1 = 0.5 * vobs[k*lenEntry+2] / f1 * d1;
	    t2 = 0.5 * vobs[k*lenEntry+2] * ((d2 / f1) -0.5 * (d1 * d1) / f2);
	    p1[i] += t1;
	    p1[j] -= t1;
	    p2[i] += t2;
	    p2[j] += t2;
	  } /* end if valid data */
	} /* end loop  L250: */

	/* Update solutions */
	for (i=0; i<nt; i++) { /* loop 260 */
	  if (swt[i] > 0.0) { /*goto L260;*/
	    xx = atan2 (xidg[i], xrdg[i]);
	    yy = atan2 (p1[i], p2[i]);
	    xrdg[i] = cos (xx - w * yy);
	    xidg[i] = sin (xx - w * yy);
	  } /* end if valid antenna data */
	} /* end loop  L260: */
      } /* end Phase only section */
	

      /* Convergence test */
      convgd = TRUE; /* L270: */
      for (i=0; i<nt; i++) { /* loop 280 */
	if ((swt[i] > 0.0)  &&  convgd) { /*goto L275;*/
	  x = tol * (5.0e-7 + sqrt (xrdg[i] * xrdg[i] + xidg[i] * xidg[i]));
	  if (fabs (xrdg[i]-xrglas[i]) > x) convgd = FALSE;
	  if (fabs (xidg[i]-xiglas[i]) > x) convgd = FALSE;
	} /* end if valid and still "converged" L275: */
	xrglas[i] = xrdg[i];
	xiglas[i] = xidg[i];
      } /* end loop  L280: */
      
      /* Print iteration value is requested */
      if (prtlv >= 5) { /*goto L300; */
	s = 0.0;
	for (k=0; k<numBL; k++) { /* loop 290 */
	  if (vobs[k*lenEntry+2] > 0.0) { /*goto L290;*/
	    i = ant1[k]-1;
	    j = ant2[k]-1;
	    xrz = vobs[k*lenEntry]   - (xrdg[i] * xrdg[j] + xidg[i] * xidg[j]);
	    xiz = vobs[k*lenEntry+1] - (xrdg[j] * xidg[i] - xrdg[i] * xidg[j]);
	    qq = xrz * xrz + xiz * xiz;
	    s += vobs[k*lenEntry+2] * sqrt (qq + eps);
	  } /* end if valid */
	} /* end loop  L290: */;
	rmean = s / sumwt;
	g_snprintf (msgtxt,80,"iter= %5d s=%15.5e rmean=%15.5e", it, s, rmean);
	Obit_log_error(err, OBIT_InfoErr, msgtxt);
      } /* end print */

      /* Inner loop converged? */
      if (convgd) break;  /* goto L320; L300: */
      
    } /* end inner solution loop L310: */;
  } /* end outer solution loop  L320: */;

  /* Refer solutions to reference antenna. */
  iref = (*nref)-1;
  g[iref*2]   = xrdg[iref];
  g[iref*2+1] = xidg[iref];
  xx = 1.0 / sqrt (g[iref*2]*g[iref*2] + g[iref*2+1]*g[iref*2+1]);
  for (i=0; i<nt; i++) { /* loop 600 */
    g[i*2]  = xrdg[i];
    g[i*2+1] = xidg[i];
    if ((i != iref)  &&  (swt[i] > 0.0)) { /*goto L600;*/
      xrr = g[i*2]   * g[iref*2] + g[i*2+1] * g[iref*2+1];
      xir = g[i*2+1] * g[iref*2] - g[i*2]   * g[iref*2+1];
      g[i*2]   = xx * xrr;
      g[i*2+1] = xx * xir;
    } /* end if valid */
  } /* end loop  L600: */;

  /* Adjust reference solution */
  g[iref*2]   = sqrt (g[iref*2]*g[iref*2] + g[iref*2+1]*g[iref*2+1]);
  g[iref*2+1] = 0.0;
 
  /* Print final results if requested */
  if (prtlv >= 5) {
    for (i=0; i<nt; i++) { /* loop 610 */
      if (swt[i] > 0.0) { /*goto L610;*/
	amp = sqrt (g[i*2]*g[i*2] + g[i*2+1]*g[i*2+1]);
	ph = 57.2958 * atan2 (g[i*2+1], g[i*2]);
	g_snprintf (msgtxt,80,"ant=  %5d amp=%12.5f phase=%12.2f", 
		    i+1, amp, ph);
	Obit_log_error(err, OBIT_InfoErr, msgtxt);
      } /* end if valid */
    } /* end loop  L610: */;
  } /* end print */

  *ierr = 0;  /* Seems to have worked */

  /* Cleanup */
  cleanup: if (swt) g_free(swt);
  if (p1)     g_free(p1);
  if (p2)     g_free(p2);
  if (xrdg)   g_free(xrdg);
  if (xidg)   g_free(xidg);
  if (xrglas) g_free(xrglas);
  if (xiglas) g_free(xiglas);
} /* end of routine gainCalcL1 */ 

/* Set Antenna and source lists */
/**
 * Make a deep copy of an ObitUVGSolve.
 * If only one source listed set in->curSource
 * \param in   The object to update lists
 * \param inUV Data object
 * \param suba Subarray
 * \param err  Obit error stack object.
 * \return pointer to the new object.
 */
static void SetLists (ObitUVGSolve *in, ObitUV *inUV, olong suba, ObitErr* err)
{
  olong iver;
  ObitTableSU *SUTable=NULL;
  ObitTableAN *ANTable=NULL;
  iver = 1;
  gchar *routine = "SetLists";

  /* Source list */
  in->curSource = 0;
  SUTable = newObitTableSUValue (in->name, (ObitData*)inUV, &iver, 
				 OBIT_IO_ReadOnly, 0, err);
  if (SUTable && (inUV->mySel->selectSources) && inUV->mySel->sources) {
    in->SList = ObitTableSUGetList (SUTable, err);
    if (err->error) Obit_traceback_msg (err, routine, SUTable->name);
    if (inUV->mySel->sources) in->curSource = inUV->mySel->sources[0]; 
  } else {  /* Use position from header */
    in->SList = ObitSourceListCreate ("SList", 1);
    in->SList->SUlist[0]->equinox = inUV->myDesc->equinox;
    in->SList->SUlist[0]->RAMean  = inUV->myDesc->crval[inUV->myDesc->jlocr];
    in->SList->SUlist[0]->DecMean = inUV->myDesc->crval[inUV->myDesc->jlocd];
    /* Compute apparent position */
    ObitPrecessUVJPrecessApp (inUV->myDesc, in->SList->SUlist[0]);
    in->curSource = 0;  /* This the only one and no source id in data */
  }

  SUTable = ObitTableSUUnref(SUTable);   /* Done with table */

  /* Antenna list */
  iver = MAX (1, suba);
  ANTable = newObitTableANValue (in->name, (ObitData*)inUV, &iver, 
				 OBIT_IO_ReadOnly, 0, 0, err);
  in->AList = ObitTableANGetList (ANTable, err);
  ANTable = ObitTableANUnref(ANTable);   /* Done with table */
  if (err->error) Obit_traceback_msg (err, routine, ANTable->name);
} /* end SetLists */
