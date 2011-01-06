/* $Id: ObitUVGSolveWB.c 163 2010-03-01 15:05:52Z bill.cotton $ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2010                                               */
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
#include "ObitUVGSolveWB.h"
#include "ObitMem.h"
#include "ObitTableUtil.h"
#include "ObitUVUtil.h"
#include "ObitUVSoln.h"
#include "ObitTableANUtil.h"
#include "ObitTableSUUtil.h"
#include "ObitPrecess.h"
#include "ObitSinCos.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitUVGSolveWB.c
 * ObitUVGSolveWB class function definitions.
 * This class enables wideband self calibration of ObitUV data sets
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitUVGSolveWB";

/** Function to obtain parent ClassInfo - ObitUVGSolve */
static ObitGetClassFP ObitParentGetClass = ObitUVGSolveGetClass;

/**
 * ClassInfo structure ObitUVGSolveWBClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitUVGSolveWBClassInfo myClassInfo = {FALSE};

/*--------------- File Global Variables  ----------------*/

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitUVGSolveWBInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitUVGSolveWBClear (gpointer in);

/** Private:  Initialize fitter */
static void 
SetupFitter (ObitUVGSolveWB *in, ObitUV *inUV, ObitErr* err);

/** Private:  Read and average next solution interval */
static gboolean 
NextAvgWB (ObitUV *inUV, ofloat interv, ofloat* antwt, 
	   ofloat uvrang[2], ofloat wtuv, ScanData *scanData, 
	   olong *nextVisBuf, ObitErr* err);

/** Private: initial (FFT) solutions for a solution interval */
static void 
initSolve (ObitUVGSolveWB *in, ObitErr *err);

/** Private: final (LS) solutions for a solution interval */
static void 
finalSolve (ObitUVGSolveWB *in, ObitErr *err);

/** Private: find best reference antenna */
static olong 
GetRefAnt (ScanData *scanData, ObitErr *err);

/** Private: initial (FFT) solution for one antenna  */
static gboolean 
initAntSolve (ObitUVGSolveWB *in, olong iAnt, olong refAnt, ObitErr *err);

/** Private: Stack baselines for one antenna  */
static void 
stackAnt (ObitUVGSolveWB *in, olong iAnt, ObitErr *err);

/** Private: Determine SNR of solution */
static void   
calcSNR (ObitUVGSolveWB *in, ofloat snrmin, gboolean *gotAnt, olong prtlv, 
	 ObitErr *err);  

/** Private: Initial guess for big least squares */
static gboolean
fitIFData (ObitUVGSolveWB *in, olong iAnt, olong iPoln, 
	   ofloat *phase, ofloat *delay, ofloat *disp, ObitErr* err);

/** Private: Set Class function pointers. */
static void ObitUVGSolveWBClassInfoDefFn (gpointer inClass);

/* Set Antenna and source lists */
static void SetLists (ObitUVGSolveWB *in, ObitUV *inUV, olong suba, 
		      ObitErr* err);

/** Private: Make Baseline data structure */
static BaselineData* MakeBLData (olong ant1, olong ant2, 
				 olong numIF, olong numPoln,
				 gboolean avgPoln, olong numFreq);

/** Private: Clear Baseline data structure */
static void ClearBLData (BaselineData *in, gboolean avgPoln);

/** Private: Delete Baseline data structure */
static BaselineData* KillBLData (BaselineData *in, gboolean avgPoln);

/** Private: Make Scan data structure */
static ScanData* MakeScanData (ObitUV *inUV, gboolean avgPoln);

/** Private: Clear Scan data structure */
static void ClearScanData (ScanData *in);

/** Private: Delete Scan data structure */
static ScanData* KillScanData (ScanData *in);

/** Private: Find peak in FFT delay function */
static void  FindFFTPeak (ObitUVGSolveWB *in, ofloat *ppos, ofloat *pval);

/** Private: Zero elements of array in1 where array in2 is blanked */
void ObitFArrayZeroBlank (ObitFArray* in1, ObitFArray* in2, ObitFArray* out);

/** Private: Fill last half of CArray with cos/sin of phases in in */
void ObitCArrayFSinCos (ObitFArray* in, ObitCArray* out);

/** Private: Multiply a CArray by an FArray,  second adjusted to end of first*/
void ObitCArrayFMulEnd (ObitCArray* Cin, ObitFArray* Fin, ObitCArray* out);

/** Private: Add elements of two FArrays, second adjusted to end of first  */
void ObitFArrayAddEnd (ObitFArray* in1, ObitFArray* in2, ObitFArray* out);

/** Private: Harmonic sum of elements of two FArrays, second adjusted to end of first  */
void ObitFArrayHarmAddEnd (ObitFArray* in1, ObitFArray* in2, ObitFArray* out);

/** Private: Add square of sum of elements of two FArrays, 
    second adjusted to end of first  */
void ObitFArrayAddEnd2 (ObitFArray* in1, ObitFArray* in2, ObitFArray* out, olong *count);

/** Private: Accumulate square of harmonic sum of elements of two FArrays, 
    second adjusted to end of first  */
void ObitFArrayHarmAccEnd2 (ObitFArray* in1, ObitFArray* in2, ObitFArray* out,
			    olong *count);

/** Privaste: Calculate SNR from decorrelation */
static ofloat 
decorSNR (ofloat amp, ofloat sumw, ofloat sumww, ofloat tfact, olong count);

#if HAVE_GSL==1  /* GSL stuff */
/** Private Fringe fitting function calculating model  */
static int fringeFitFunc (const gsl_vector *coef, void *params, gsl_vector *f);
/** Private Fringe fitting function calculating Jacobean  */
static int fringeFitJacob (const gsl_vector *coef, void *params, gsl_matrix *J);
/** Private Fringe fitting function calculating model and Jacobean  */
static int fringeFitFuncJacob (const gsl_vector *coef, void *params, gsl_vector *f, 
			       gsl_matrix *J);
/** Private Coarse Fringe fitting function calculating model  */
static int coarseFitFunc (const gsl_vector *coef, void *params, gsl_vector *f);
/** Private CoarseFringe fitting function calculating Jacobean  */
static int coarseFitJacob (const gsl_vector *coef, void *params, gsl_matrix *J);
/** Private CoarseFringe fitting function calculating model and Jacobean  */
static int coarseFitFuncJacob (const gsl_vector *coef, void *params, gsl_vector *f, 
			       gsl_matrix *J);
#endif /* GSL stuff */


/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitUVGSolveWB* newObitUVGSolveWB (gchar* name)
{
  ObitUVGSolveWB* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitUVGSolveWBClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitUVGSolveWB));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitUVGSolveWBInit((gpointer)out);

 return out;
} /* end newObitUVGSolveWB */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitUVGSolveWBGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitUVGSolveWBClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitUVGSolveWBGetClass */

/**
 * Make a deep copy of an ObitUVGSolveWB.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitUVGSolveWB* ObitUVGSolveWBCopy  (ObitUVGSolveWB *in, ObitUVGSolveWB *out, ObitErr *err)
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
    out = newObitUVGSolveWB(outName);
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
} /* end ObitUVGSolveWBCopy */

/**
 * Make a copy of a object but do not copy the actual data
 * This is useful to create an UVGSolveWB similar to the input one.
 * \param in  The object to copy
 * \param out An existing object pointer for output, must be defined.
 * \param err Obit error stack object.
 */
void ObitUVGSolveWBClone  (ObitUVGSolveWB *in, ObitUVGSolveWB *out, ObitErr *err)
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

} /* end ObitUVGSolveWBClone */

/**
 * Creates an ObitUVGSolveWB 
 * \param name      An optional name for the object.
 * \return the new object.
 */
ObitUVGSolveWB* ObitUVGSolveWBCreate (gchar* name)
{
  ObitUVGSolveWB* out;

  /* Create basic structure */
  out = newObitUVGSolveWB (name);

  return out;
} /* end ObitUVGSolveWBCreate */

/**
 * Determine phase or amp & phase calibration for an UV dataset divided by
 * a source model.
 * If the output table previously exists, deselect any entries corresponding to 
 * data selected on inUV (to be fitted).
 * Routine translated from the AIPSish UVUTIL.FOR/SLFCAL 
 * \param in      Input gain fitter object. 
 * If only one source fitted in->curSource, has it's ID
 * Control parameters are on the info member.
 * \li "solnVer" OBIT_int   (1,1,1) Solution (SN) table to write; 0=> create new.
 * \li "subA"    OBIT_int   (1,1,1) Selected subarray (default 1)
 * \li "solInt"  OBIT_float (1,1,1) Solution interval (min). (default scan)
 * \li "refAnt"  OBIT_int   (1,1,1) Single Ref ant to use. (default most common)
 * \li "refAnts" OBIT_int   (?,1,1) list of Ref ants to use. (default (most common))
 * \li "avgPol"  OBIT_bool  (1,1,1) True if RR and LL to be averaged (false)
 * \li "avgIF"   OBIT_bool  (1,1,1) True if all IFs to be averaged (false)
 *                                  otherwise individual IF fits.
 * \li "doTwo"   OBIT_bool  (1,1,1) Use 2 BL combinations as well as 1 BL. (true)
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
ObitTableSN* ObitUVGSolveWBCal (ObitUVGSolveWB *in, ObitUV *inUV, ObitUV *outUV, 
				ObitUVSel *sel, ObitErr *err)
{
  ObitTableSN *outSoln=NULL;
  ObitTableSNRow *row=NULL;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  olong i, nif, npoln, iAnt, numBL, SNver;
  olong nextVisBuf, iSNRow, numFreq, cntGood=0, cntPoss=0, cntBad=0;
  oint numPol, numIF, numAnt, suba, minno, prtlv, mode;
  ofloat solInt, snrmin, uvrang[2], wtuv, FractOK, minOK=0.1;
  olong itemp, kday, khr, kmn, ksec;
  ofloat *antwt=NULL, fblank = ObitMagicF();
  gboolean avgpol, dol1, empty;
  gboolean done, good, oldSN, *gotAnt;
  gchar soltyp[5], solmod[5];
  odouble timex;
  olong sid;
  ObitIOCode retCode;
  gchar *tname, *ModeStr[] = {"A&P", "P", "P!A"};
  gchar *routine = "ObitUVGSolveWBCal";
  
  /* error checks */
  if (err->error) return outSoln;
  g_assert (ObitUVGSolveWBIsA(in));
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
  
  /* Create output SN table - version requested? */
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
 
  /* Which subarray? */
  suba = 1;
  ObitInfoListGetTest(in->info, "subA", &type, dim, &suba);
  /* Can only do one */
  Obit_retval_if_fail((suba>0 && suba<=inUV->myDesc->numSubA), err, outSoln,
		      "%s: MUST specify a single subarray for %s", 
		      routine, inUV->name);
  
  /* Create arrays */
  numAnt = inUV->myDesc->numAnt[suba-1];
  antwt  = g_malloc0(numAnt*sizeof(ofloat));
  gotAnt = g_malloc0(numAnt*sizeof(gboolean));
  numBL  = (numAnt * (numAnt-1)) / 2;
  
  /* Get parameters from inUV */
  avgpol = FALSE;
  ObitInfoListGetTest(in->info, "avgPol", &type, dim, &avgpol);
  if (numPol<=1) avgpol = FALSE;  /* Don't average if only one */
  snrmin = 5.0;
  ObitInfoListGetTest(in->info, "minSNR", &type, dim, &snrmin);
  in->minSNR = snrmin;
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
  mode = 2;  /* Delay in P!A mode */
  if (prtlv>=1) {
    Obit_log_error(err, OBIT_InfoErr, "Calibrate delay in %s mode", ModeStr[mode]);
    ObitErrLog(err);
  }

  /* Initialize object */
  SetupFitter (in, inUV, err);
  if (err->error) goto cleanup;
  
  sid = in->curSource;  /* In case only one */
  /* Get antenna and source lists */
  SetLists (in, inUV, sel->SubA, err);
 
  /* Averaging of data? */
  nif = numIF;
  npoln = numPol;
  if (avgpol) npoln = 1;
  
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
  row->SubA   = 1; 
  row->FreqID = 1; 
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
    done =  NextAvgWB (inUV, solInt, antwt, uvrang, wtuv, 
		       in->scanData, &nextVisBuf, err);
    if (err->error) goto cleanup;

    /* Done? 
    if (done) break; still have OK data */
    
    /* Write time if requested */
    if (prtlv >= 3) {
      kday = in->scanData->timec;
      timex = (in->scanData->timec - kday) * 24.;
      khr = timex;
      timex = (timex - khr) * 60.;
      kmn = timex;
      timex = (timex - kmn) * 60.;
      ksec = timex + 0.5;
      Obit_log_error(err, OBIT_InfoErr, " time=%d/ %3d %3d %3d ", 
		     kday, khr, kmn, ksec);
    } 
    
    /* Do solutions */
    initSolve (in, err);
    finalSolve (in, err);
    /* Get SNR of solutions */
    calcSNR (in, snrmin, gotAnt, prtlv, err);
    /* Messages */
    if (prtlv>1) ObitErrLog(err);
    
    /* How many possible solutions */
    for (iAnt= 0; iAnt<numAnt; iAnt++) if (gotAnt[iAnt]) cntPoss += numIF*numPol;
    
    /* Write solutions to SN table */
    /* Common values */
    row->Time   = in->scanData->timec; 
    row->TimeI  = in->scanData->timei; 
    if (in->scanData->sid>0) row->SourID = in->scanData->sid; 
    row->SubA   = MAX (1, in->scanData->suba); 
    row->FreqID = in->scanData->fqid; 
    
    /* Loop over antennas - write as corrections rather than solutions */
    iSNRow = -1;
    for (iAnt= 0; iAnt<numAnt; iAnt++) {
      if (gotAnt[iAnt]) {
	good = FALSE; /* Until proven */
	row->antNo  = iAnt+1; 
	if (in->scanData->avgIF)
	    row->MBDelay1 = -in->antDelay[iAnt*numIF*numPol]; 
	for (i=0; i<numIF; i++) {
	  row->Real1[i]   =  in->antGain[(iAnt*numIF*numPol+i)*2];
	  row->Imag1[i]   =  in->antGain[(iAnt*numIF*numPol+i)*2+1];
	  row->Delay1[i]  =  in->antDelay[iAnt*numIF*numPol+i];
	  row->Weight1[i] =  in->antWeight[iAnt*numIF*numPol+i];
	  if (row->Weight1[i]<=0.0) {
	    row->Real1[i]   = fblank;
	    row->Imag1[i]   = fblank;
	    row->Delay1[i]  = fblank;
	  }
	  row->RefAnt1[i] = in->refAnt;
	  if (row->Weight1[i]>0.0) {good = TRUE; cntGood++;}
	  if (row->Weight1[i]<=0.0) cntBad++;    /* DEBUG */
	}
	if (numPol>1) {
	  if (in->scanData->avgIF)
	    row->MBDelay2 = -in->antDelay[iAnt*numIF*numPol+numIF]; 
	  for (i=0; i<numIF; i++) {
	    row->Real2[i]   =  in->antGain[(iAnt*numIF*numPol+i+numIF)*2];
	    row->Imag2[i]   =  in->antGain[(iAnt*numIF*numPol+i+numIF)*2+1];
	    row->Delay2[i]  =  in->antDelay[iAnt*numIF*numPol+i+numIF];
	    row->Weight2[i] =  in->antWeight[iAnt*numIF*numPol+i+numIF];
	    if (row->Weight2[i]<=0.0) {
	      row->Real2[i]   = fblank;
	      row->Imag2[i]   = fblank;
	      row->Delay2[i]  = fblank;
	    }
	    row->RefAnt2[i] = in->refAnt;
	    if (row->Weight2[i]>0.0) {good = TRUE; cntGood++;}
	    if (row->Weight2[i]<=0.0) cntBad++; 
	  }
	}
	retCode = ObitTableSNWriteRow (outSoln, iSNRow, row, err);
	if (err->error) goto cleanup;
      } /* end if gotant */
    }
  } /* end loop processing data */
  
  /* Close input uv data */
  ObitUVClose (inUV, err);
  if (err->error) goto cleanup;
  
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
  if (prtlv>=2) {
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
 cleanup: 
  row = ObitTableSNUnref(row);
  if (antwt)  g_free(antwt);
  if (gotAnt) g_free(gotAnt);
 
  /* Shutdown  */
  in->scanData = KillScanData (in->scanData);
   
 if (err->error) Obit_traceback_val (err, routine, inUV->name, outSoln);
  
  return outSoln;
} /* end ObitUVGSolveWBCal */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitUVGSolveWBClassInit (void)
{
  ofloat phase=0.5, cp, sp;
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitUVGSolveWBClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */

  /* Init Sine/Cosine calculator - just to be sure about threading */
  ObitSinCosCalc(phase, &sp, &cp);
 
} /* end ObitUVGSolveWBClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitUVGSolveWBClassInfoDefFn (gpointer inClass)
{
  ObitUVGSolveWBClassInfo *theClass = (ObitUVGSolveWBClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitUVGSolveWBClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitUVGSolveWBClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitUVGSolveWBGetClass;
  theClass->newObit       = (newObitFP)newObitUVGSolveWB;
  theClass->ObitCopy      = (ObitCopyFP)ObitUVGSolveWBCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitUVGSolveWBClear;
  theClass->ObitInit      = (ObitInitFP)ObitUVGSolveWBInit;
  theClass->ObitUVGSolveCreate = (ObitUVGSolveCreateFP)ObitUVGSolveWBCreate;
  theClass->ObitUVGSolveCal    = (ObitUVGSolveCalFP)ObitUVGSolveWBCal;

} /* end ObitUVGSolveWBClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitUVGSolveWBInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitUVGSolveWB *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->scanData    = NULL;
  in->antGain     = NULL;
  in->antDelay    = NULL;
  in->antDisp     = NULL;
  in->antWeight   = NULL;
  in->refAnts     = NULL;
  in->cWork1      = NULL;
  in->cWork2      = NULL;
  in->cWork3      = NULL;
  in->FFTFitArray = NULL;
  in->fWork1      = NULL;
  in->fWork2      = NULL;
  in->fWork3      = NULL;
  in->fWorkWt2    = NULL;
  in->myFFT       = NULL;
  in->myCInterp   = NULL;
  in->FFTOverSamp = 4;
#if HAVE_GSL==1  /* GSL stuff */
  in->myCoarseSolver   = NULL;
  in->coarse_coef_init = NULL;
  in->coarseData       = NULL;
  in->coarseWt         = NULL;
  in->coarseFreq       = NULL;
#endif /* GSL stuff */


} /* end ObitUVGSolveWBInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *         Actually it should be an ObitUVGSolveWB* cast to an Obit*.
 */
void ObitUVGSolveWBClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitUVGSolveWB *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  in->scanData     = KillScanData (in->scanData);
  in->cWork1       = ObitCArrayUnref(in->cWork1);
  in->cWork2       = ObitCArrayUnref(in->cWork2);
  in->cWork3       = ObitCArrayUnref(in->cWork3);
  in->FFTFitArray  = ObitCArrayUnref(in->FFTFitArray);
  in->fWork1       = ObitFArrayUnref(in->fWork1);
  in->fWork2       = ObitFArrayUnref(in->fWork2);
  in->fWork3       = ObitFArrayUnref(in->fWork3);
  in->fWorkWt2     = ObitFArrayUnref(in->fWorkWt2);
  in->myFFT        = ObitFFTUnref(in->myFFT);
  in->myCInterp    = ObitCInterpolateUnref(in->myCInterp);
  if (in->antGain)   g_free(in->antGain);
  if (in->antDelay)  g_free(in->antDelay);
  if (in->antDisp)   g_free(in->antDisp);
  if (in->antWeight) g_free(in->antWeight);
  if (in->refAnts)   g_free(in->refAnts);
#if HAVE_GSL==1  /* GSL stuff */
  if (in->myCoarseSolver)   gsl_multifit_fdfsolver_free(in->myCoarseSolver);
  if (in->myCoarseFunc)     g_free(in->myCoarseFunc);
#endif /* GSL stuff */
  if (in->coarse_coef_init) g_free(in->coarse_coef_init);
  if (in->coarseData)       g_free(in->coarseData);
  if (in->coarseWt)         g_free(in->coarseWt);
  if (in->coarseFreq)       g_free(in->coarseFreq);

  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitUVGSolveWBClear */

/** Initialize fitter 
 * \param in      Input gain fitter object. 
 * \param inUV    Input UV data. 
 * \param err    Error/message stack, returns if error.
 */
static void 
SetupFitter (ObitUVGSolveWB *in, ObitUV *inUV, ObitErr* err)
{
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  olong suba, naxis[1], i, *refAnts;
  gboolean avgpol, avgIF;
#if HAVE_GSL==1  /* GSL stuff */
  const gsl_multifit_fdfsolver_type *T;
  olong ndata, ncoef, refChan;
#endif /* GSL stuff */

  /* error checks */
  if (err->error) return;

 /* Get parameters from inUV */
  if (inUV->myDesc->jlocs>=0)
    in->numPoln    = MIN (2, inUV->myDesc->inaxes[inUV->myDesc->jlocs]);
  else in->numPoln = 1;
  if (inUV->myDesc->jlocif>=0)
    in->numIF  = inUV->myDesc->inaxes[inUV->myDesc->jlocif];
  else in->numIF  = 1;
  in->prtLv = 0;
  ObitInfoListGetTest(in->info, "prtLv", &type, dim, &in->prtLv);
  suba = 1;
  ObitInfoListGetTest(in->info, "subA", &type, dim, &suba);
  in->maxAnt  = inUV->myDesc->numAnt[suba-1];
  avgpol = FALSE;
  ObitInfoListGetTest(in->info, "avgPol", &type, dim, &avgpol);
  if (in->numPoln<=1) avgpol = FALSE;  /* Don't average if only one */
  avgIF = FALSE;
  ObitInfoListGetTest(in->info, "avgIF", &type, dim, &avgIF);
  if (in->numIF<=1) avgIF = FALSE;  /* Don't average if only one */

  /* Generate scanData structure */
  in->scanData = MakeScanData (inUV, avgpol);

  in->scanData->avgPoln = avgpol;
  in->scanData->avgIF   = avgIF;
  in->scanData->doTwo   = TRUE;
  ObitInfoListGetTest(in->info, "doTwo", &type, dim, &in->scanData->doTwo);

  /* Setup for interpolation of FFT result */
  in->FFTOverSamp = 4;   /* Over sampling factor for FFT search */
  naxis[0] = in->FFTOverSamp*ObitFFTSuggestSize(in->scanData->BLData[0]->numFreq);
  in->FFTFitArray = ObitCArrayCreate("FFTFit work", 1, naxis);
  in->myCInterp   = newObitCInterpolateCreate ("FFT Search", in->FFTFitArray, NULL, 
					       1.0, 1.0, 2, 2, err);

  /* Reference antenna(s) */
  in->refAnt  = -1;
  ObitInfoListGetTest(in->info, "refAnt", &type, dim, &in->refAnt);
  /* Determine reference antenna is none given */
  if (in->refAnt<=0) 
    in->refAnt = GetRefAnt(in->scanData, err);
  /* Check for list - zero terminated list */
  if (ObitInfoListGetP(in->info, "refAnts",  &type, dim, (gpointer)&refAnts)) {
    in->refAnts = g_malloc0((dim[0]+1)*sizeof(olong));
    for (i=0; i<dim[0]; i++) in->refAnts[i] = refAnts[i]; in->refAnts[i] = 0;
    /* make sure at least one non zero */
    if (in->refAnts[0]==0) in->refAnts[0] = in->refAnt;
  } else { /* Use single antenna list */
    in->refAnts = g_malloc0(2*sizeof(olong));
    in->refAnts[0] = in->refAnt; in->refAnts[1] = 0; 
  }

  /* Setup for coarse IF result fitting */
#if HAVE_GSL==1  /* GSL stuff */
  ndata = in->numIF;  /* Number of data points */
  ncoef = 2;          /* Fitting phase and delay */
  T = gsl_multifit_fdfsolver_lmsder;
  /* T = gsl_multifit_fdfsolver_lmder; Faster */
  in->myCoarseSolver   = gsl_multifit_fdfsolver_alloc(T, ndata, ncoef);
  in->coarse_coef_init =  g_malloc(ncoef*sizeof(double));
  in->ncoefCoarse      = ncoef;
  in->ndataCoarse      = ndata;
  /* Create /fill function structure */
  in->myCoarseFunc      = g_malloc0(sizeof(gsl_multifit_function_fdf));
  in->myCoarseFunc->f   = &coarseFitFunc;      /* Compute function */
  in->myCoarseFunc->df  = &coarseFitJacob;     /* Compute Jacobian (derivative matrix) */
  in->myCoarseFunc->fdf = &coarseFitFuncJacob; /* Compute both function and derivatives */
  in->myCoarseFunc->n   = ndata;               /* Number of data points */
  in->myCoarseFunc->p   = ncoef;               /* number of parameters */
  in->myCoarseFunc->params = in;               /* Object */
  in->coarseData       = g_malloc0(ndata*sizeof(ofloat));  /* Data array */
  in->coarseWt         = g_malloc0(ndata*sizeof(ofloat));  /* Weight array */
  in->coarseFreq       = g_malloc0(ndata*sizeof(ofloat));  /* Frequency offset array */
  refChan = (olong)(in->scanData->refChan+0.5) - 1; /* Ref. channel */
  refChan = MAX (0, refChan);
  refChan = MIN (refChan, in->scanData->BLData[0]->numFreq-1);
  for (i=0; i<ndata; i++) {
    in->coarseData[i] = 0.0;
    in->coarseWt[i]   = 1.0;
    in->coarseFreq[i] = in->scanData->BLData[0]->dFreq[i][refChan]; 
  }
#endif /* GSL stuff */

} /*  End SetupFitter */
/**
 * Average next solution interval
 *  Data is averaged until one of several conditions is met:
 *  \li the time exceeds the initial time plus the specified interval.
 *  \li the source id (if present) changes
 *  \li the FQ id (if present) changes. 
 * \param inUV    Input UV data. 
 * \param interv  Desired time interval of average in days. 
 *                Will use value from inUV->mySel.ObitUVSelSubScan
                  if available, 0=> scan average
 * \param antwt   Extra weights to antennas (>=0 => 1.0)
 * \param uvrang  Range of baseline lengths with full weight (lamda). 
 *                0s => all baselines 
 * \param wtuv    Weight outside of UVRANG. (No default)
 * \param scanData Data structure for averaged data
 * \param nextVisBuf [in/out] next vis (0-rel) to read in buffer, 
 *                   -1=> read, -999 done
 * \param err    Error/message stack, returns if error.
 * \return TRUE if all data read, else FALSE
 */
static gboolean 
NextAvgWB (ObitUV* inUV, ofloat interv, ofloat* antwt, 
	   ofloat uvrang[2], ofloat wtuv,  ScanData *scanData, olong *nextVisBuf, 
	   ObitErr* err)
{
  gboolean done=FALSE;
  BaselineData *BLData;
  ObitIOCode retCode= OBIT_IO_OK;
  ofloat linterv, ctime, cbase, stime, ltime=0, weight, wt, bl, *visPnt, temp;
  ofloat *wtArray, *visArray, *phArray, fblank = ObitMagicF();
  olong i, ip, csid, cfqid, a1, a2, *blLookup=NULL, blIndex, visIndex;
  olong iFreq, iIF, iStok, jBL, jIF, jStok, mPol, mIF, numBL, suba;
  olong numAnt, numFreq, numIF, numPol;
  odouble timeSum;
  olong timeCount;
  gboolean avgpol, gotData;
  gchar *routine = "ObitUVGSolveWB:NextAvg";
  
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

  /* Numbers of things */
  numAnt  = scanData->maxAnt;
  numBL   = scanData->numBase;
  numFreq = scanData->BLData[0]->numFreq;
  numIF   = scanData->BLData[0]->numIF;
  numPol  = scanData->BLData[0]->numPoln;
  avgpol = scanData->avgPoln;

  /* Numbers of poln and IFs in accumulator */
  if (avgpol) mPol = 1;
  else mPol = numPol;
  mIF = numIF;  /* Never average IF */

  /* Baseline lookup table - assumes  a1 < a2 */
  blLookup    = g_malloc0(numAnt*sizeof(olong));
  blLookup[0] = 0;
  for (i=1; i<numAnt; i++) blLookup[i] = blLookup[i-1] + numAnt-i;
  
  /* Zero accumulations */
  ClearScanData (scanData);

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
    
    /* Visibility pointer */
    visPnt = inUV->buffer + (*nextVisBuf) * inUV->myDesc->lrec;
    /* Vis data */
    ctime = visPnt[inUV->myDesc->iloct]; /* Time */
    if (stime<-1000.0) stime = ctime;    /* Set start time from first vis */
    if (inUV->myDesc->ilocsu>=0) csid  = (olong)visPnt[inUV->myDesc->ilocsu]; /* Source */
    else csid  = 0;
    if (timeCount==0) scanData->sid = csid;  /* Set output source id first vis */
    if (inUV->myDesc->ilocfq>=0) cfqid = (olong)visPnt[inUV->myDesc->ilocfq]; /* FQid */
    else cfqid  = 0;
    if (scanData->fqid<0) scanData->fqid = cfqid;  /* Set output fq id first vis */
    
    /* Is this integration done? */
    if ((scanData->sid!=csid) || (scanData->fqid!=cfqid) || (ctime-stime>=linterv)) break;
    
    /* Sum time */
    timeSum += ctime;
    timeCount++;
    ltime = ctime;   /* Last time in accumulation */
    
    /* crack Baseline */
    cbase = visPnt[inUV->myDesc->ilocb]; /* Baseline */
    suba  = 1 + 0.01*(visPnt[inUV->myDesc->ilocb] - (olong)cbase);
    a1 = (cbase / 256.0) + 0.001;
    a2 = (cbase - a1 * 256) + 0.001;
    if (timeCount==0) scanData->suba = suba;  /* Set subarray first vis */

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
    BLData = scanData->BLData[blIndex];  /* Baseline data structure */

    /* Check consistency with BaselineData */
    if (!((a1==BLData->ant1) && (a2==BLData->ant2))) {
     Obit_log_error(err, OBIT_Error,
		    "%s: Antenna numbers inconsistent, %d != %d or %d != %d", 
		    routine, a1, a2, BLData->ant1, BLData->ant2);
     goto cleanup;
    }

    /* Accumulate data */
    /* Loop over polarization */
    ip = 0;
    for (iStok = 0; iStok<numPol; iStok++) {
      if (avgpol) ip = 0;  /* use first poln */
      /* Loop over IF */
      for (iIF = 0; iIF<numIF; iIF++) {
	
	/* Get accumulators */
	visArray = BLData->visArray[ip]->array;
	wtArray  = BLData->wtArray[ip]->array;

	/* Loop over frequency */
	for (iFreq = 0; iFreq<numFreq; iFreq++) {

	  /* Visibity index */
	  visIndex = inUV->myDesc->nrparm + iStok*inUV->myDesc->incs + 
	    iIF*inUV->myDesc->incif + iFreq*inUV->myDesc->incf;
	  
	  /* Accumulate phases and weights */
	  if (visPnt[visIndex+2]>0.0) {
	    wt = MAX (0.0, weight * visPnt[visIndex+2]);
	    visArray[2*iFreq]   += wt * visPnt[visIndex];
	    visArray[1+2*iFreq] += wt * visPnt[visIndex+1];
	    wtArray[iFreq]      += wt;
	    BLData->WtPolnIF[ip]+= wt;
	  } /* End vis valid */
	} /* Loop over Freq */
	ip++;
     } /* Loop over IF */
    } /* Loop over Stokes */
    
    (*nextVisBuf)++;  /* Increment vis being processed */
  } /* End loop summing */
  
  /* Normalize by sum of weights */
  for (jBL=0; jBL<numBL; jBL++) {         /* Loop over baseline */
    BLData = scanData->BLData[jBL];       /* Baseline data structure */
    ip = 0;
    for (jStok=0; jStok<mPol; jStok++) {  /* Loop over Pol */
      for (jIF=0; jIF<mIF; jIF++) {       /* Loop over IF */
	/* Accumulators */
	visArray = BLData->visArray[ip]->array;
	phArray  = BLData->phArray[ip]->array;
	wtArray  = BLData->wtArray[ip]->array;
	/* Loop over frequency */
	for (iFreq = 0; iFreq<numFreq; iFreq++) {
	  if (wtArray[iFreq]>0) {
	    visArray[2*iFreq]   /= wtArray[iFreq];
	    visArray[1+2*iFreq] /= wtArray[iFreq];
	    phArray[iFreq] = atan2(visArray[1+2*iFreq], visArray[2*iFreq]);
	  } else {
	    visArray[2*iFreq]   = fblank;
	    visArray[1+2*iFreq] = fblank;
	    phArray[iFreq]      = fblank;
	  }
 	} /* Loop over Freq */
	ip++;
      } /* end loop over IF */
    } /* end loop over Pol */
  } /* end loop over baseline */
  
  /* Average time and interval for output */
  if (timeCount>0) scanData->timec = timeSum/timeCount;
  else scanData->timec = 0.0;
  scanData->timei = (ltime - stime);
  
  /* Cleanup */
  cleanup: if (blLookup) g_free(blLookup);
  if (err->error) Obit_traceback_val (err, routine, inUV->name, done);
  return done;
} /* end NextAvgWB */

/**
 * calcSNR computes antenna based signal-to-noise ratioes (SNR) based  
 * on phase residuals from a gain model.  The approximation used is  
 * that SNR = 1.0 / RMS phase residual.  
 * Does a weighted solution for the SNR.  If insufficient data is  
 * present to compute an RMS but at least one observation exists, the  
 * SNR is set to 6.0.  
 * If the previous weight is larger that the calculated one, 
 * the previous is used.
 * \param in      Solver data structure
 * \param snrmin  Minimum SNR allowed. 
 * \param count   A work array used for the counts for each 
 *                antenna, must be at least MAXANT in size. 
 * \param gotAnt  Array of flags if antenna found
 * \param prtlv   Print level: 0 none, >5 the antenna SNRs 
 * \param err     Error/message stack, returns if error.
 */
static void 
calcSNR (ObitUVGSolveWB *in, ofloat snrmin, gboolean *gotAnt, olong prtlv, 
	 ObitErr *err) 
{
  BaselineData *BLData;
  olong   iAnt1, iAnt2, loop, nprt, *count=NULL, nval, offset, ip, jp;
  ofloat  zr, zi, zzr, zzi, prtsnr;
  ofloat  *wtArray, *phArray, phase1, phase2, phase, amp2, ph1, ph2;
  ofloat *sumwt=NULL, *snr=NULL;
  olong  npoln, ipoln, iIF, iFreq, iAnt, indx, *error=NULL;

  /* Error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;  /* previous error? */

    /* Create work arrays */
  sumwt = g_malloc0(in->scanData->maxAnt*sizeof(ofloat));
  snr   = g_malloc0(in->scanData->maxAnt*sizeof(ofloat));
  error = g_malloc0(in->scanData->maxAnt*sizeof(olong));
  count = g_malloc(in->scanData->maxAnt*sizeof(olong));
 
  /* How many poln? */
  if (in->scanData->avgPoln) npoln = 1;
  else  npoln = in->numPoln;
  nval = in->numIF * in->numPoln;  /* solution entries per antenna in output */

  for (loop=0; loop<in->scanData->maxAnt; loop++) gotAnt[loop]= FALSE;

  /* Loop over polarization */
  for (ipoln=0; ipoln<npoln; ipoln++) {
    
    /* Zero sums, counts etc. */
    nprt = 0;
    for (loop=0; loop<in->scanData->maxAnt; loop++) { /* loop 10 */
      snr[loop]   = 0.0;
      sumwt[loop] = 0.0;
      count[loop] = 0;
      error[loop] = 0;
    } /* end loop  L10:  */
    
    /* Determine phase residuals. */
    for (loop=0; loop<in->scanData->numBase; loop++) { /* loop 30 */
      BLData = in->scanData->BLData[loop];       /* Baseline data structure */

      /* IF 0 phases - higher IFs have delay added */
      iAnt1 = BLData->ant1-1;
      iAnt2 = BLData->ant2-1;
      indx = (iAnt1*nval + ipoln)*2;
      phase1 = atan2(in->antGain[indx+1], in->antGain[indx]);
      indx = (iAnt2*nval + ipoln)*2;
      phase2 = atan2(in->antGain[indx+1], in->antGain[indx]);

     /* Loop over IF */
      ip = ipoln*in->numIF;  /* data product (IF/poln) indicator */
      for (iIF=0; iIF<in->numIF; iIF++) {
	/* Get accumulators */
	phArray  = BLData->phArray[ip]->array;
	wtArray  = BLData->wtArray[ip]->array;
	/* Loop over frequency */
	for (iFreq=0; iFreq<BLData->numFreq; iFreq++) {
	  if (wtArray[iFreq]>1.0e-20) {
	    gotAnt[iAnt1] = TRUE;  /* Some data for antenna */
	    gotAnt[iAnt2] = TRUE;
	    indx = (iAnt1*nval + iIF*npoln + ipoln);
	    ph1  = phase1 + in->antDelay[indx] * in->scanData->BLData[loop]->dFreq[iIF][iFreq];
	    indx = (iAnt2*nval + iIF*npoln + ipoln);
	    ph2  = phase1 + in->antDelay[indx] * in->scanData->BLData[loop]->dFreq[iIF][iFreq];
	    zr   = cos(ph1)*cos(ph2) - sin(ph1)*sin(ph2);
	    zi   = cos(ph1)*sin(ph2) + sin(ph1)*cos(ph2);
	    zzr  = cos(phArray[iFreq])*zr - sin(phArray[iFreq])*zi;
	    zzi  = cos(phArray[iFreq])*zi + sin(phArray[iFreq])*zr;
	    amp2 = zzr*zzr + zzi*zzi;
	    if (fabs(amp2)>0.0) {
	      phase = atan2 (zzi, zzr);
	      /* counts/sums of residuals squared */
	      count[iAnt1]++;
	      count[iAnt2]++;
	      snr[iAnt1]   += phase*phase*wtArray[iFreq];
	      snr[iAnt2]   += phase*phase*wtArray[iFreq];
	      sumwt[iAnt1] += wtArray[iFreq];
	      sumwt[iAnt2] += wtArray[iFreq];
	    }
	  } 
	} /* End loop over frequency */
	ip++;  /* data product (IF/poln) indicator */
      } /* End loop over IF */
    } /* end loop  L30:  */
    
    /* Convert to SNRs */
    for (loop=0; loop<in->maxAnt; loop++) { /* loop 70 */
      if (count[loop] <= 0) snr[loop] = 0.0;
      /* For insufficient data (but  some) use snrmin + 1.0 */
      if ((count[loop] >= 1)  &&  (count[loop] <= 2)) snr[loop] = snrmin + 1.0;
      if ((count[loop] >= 3)  &&  (snr[loop] > 1.0e-20)) {
	snr[loop] = 1.0 / sqrt (snr[loop] / ((count[loop]-1) * sumwt[loop]));
	snr[loop] = MIN (1000.0, snr[loop]);
	/* Make sure at least snrmin */
	if (snr[loop]<snrmin) snr[loop] = 0.0;
      }
      
      /* Save to antWeight if larger - set to zero if previous <snrmin */
      /* Loop over IF */
      offset = loop*nval;
      ip = 0 + ipoln*in->numIF;
      for (iIF=0; iIF<in->numIF; iIF++) {
	if (in->antWeight[offset+ip]>snrmin) {
	  /*in->antWeight[offset+ip] = MAX (in->antWeight[offset+ip], snr[loop]); 
	    else
	    in->antWeight[offset+ip] = 0.0;DEBUG*/
	}
	if (in->antWeight[offset+ip]<snrmin) in->antWeight[offset+ip] = 0.0;
	ip++;
      }
      
      /* Print result if desired. */
      if ((count[loop] >= 1)  &&  (prtlv >= 5)) {
	prtsnr = MIN (9999.999, snr[loop]);
	Obit_log_error(err, OBIT_InfoErr, "antenna(%2d)  %3d obs, snr = %10.3f", 
		       loop+1, count[loop], prtsnr);
      }
    } /* end loop  L70:  */;
    
  } /* End loop over poln */

  /* If averaging poln solutions, copy results */
  if ((in->scanData->avgPoln) && (in->scanData->BLData[0]->numPoln>1)) {
    nval = in->numIF * in->numPoln;  /* entries per antenna */
    
    for (iAnt=1; iAnt<=in->maxAnt; iAnt++) { /* Loop over antennas */
      offset = (iAnt-1)*nval;          /* This antennas offset on solutions */
      ip = in->scanData->BLData[0]->numIF;
      jp = 0;
      /* Loop over IF */
      for (iIF=0; iIF<in->scanData->BLData[0]->numIF; iIF++) {
	in->antWeight[offset+ip] = in->antWeight[offset+jp];
	jp++;  ip++;
      } /* End IF Loop */
    } /* end antenna loop */
  } /* End copy solutions to other polarization if averaging */

  /* Cleanup */
  if (sumwt) g_free(sumwt);
  if (snr)   g_free(snr);
  if (error) g_free(error);
  if (count) g_free(count);
} /* end of routine calcSNR */ 

/**
 * Set Antenna and source lists
 * If only one source listed set in->curSource
 * \param in   The object to update lists
 * \param inUV Data object
 * \param suba Subarray
 * \param err  Obit error stack object.
 * \return pointer to the new object.
 */
static void SetLists (ObitUVGSolveWB *in, ObitUV *inUV, olong suba, ObitErr* err)
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
				 OBIT_IO_ReadOnly, 0, 0, 0, err);
  in->AList = ObitTableANGetList (ANTable, err);
  ANTable = ObitTableANUnref(ANTable);   /* Done with table */
  if (err->error) Obit_traceback_msg (err, routine, ANTable->name);
} /* end SetLists */

/**
 * Initial guess at solutions based on FFT on baselines to reference antenna
 * Gets stacked antenna phases weithted by SNR of fit.
 * \param in   The structure with data and solutions.
 * \param err  Error/message stack, returns if error.
 */
static void 
initSolve (ObitUVGSolveWB *in, ObitErr *err)
{
  olong iAnt, i, nval, offset, iref, refAnt;
  gboolean good, test;

  if (err->error) return;  /* Prior error? */

  nval = in->numIF * in->numPoln;  /* entries per antenna */

  /* Create solution structures if needed */
  if (in->antGain==NULL)   in->antGain   = g_malloc(in->maxAnt*nval*2*sizeof(ofloat));
  if (in->antDelay==NULL)  in->antDelay  = g_malloc(in->maxAnt*nval*sizeof(ofloat));
  if (in->antDisp==NULL)   in->antDisp   = g_malloc(in->maxAnt*nval*sizeof(ofloat));
  if (in->antWeight==NULL) in->antWeight = g_malloc(in->maxAnt*nval*sizeof(ofloat));

  /* Loop over list of possible reference antennas */
  iref = 0;
  good = FALSE;
  while (!good && (in->refAnts[iref]>0)) {
    refAnt = in->refAnts[iref];

    /* Zero test reference antenna values */
    /* Use 0 phase, derivatives for reference antenna */
    offset = (refAnt-1) * nval;
    for (i=0; i<nval; i++) {
      in->antGain[2*(offset+i)]   = 1.0;
      in->antGain[2*(offset+i)+1] = 0.0;
      in->antDelay[offset+i ]     = 0.0;
      in->antDisp[offset+i]       = 0.0;
      in->antWeight[offset+i ]    = in->minSNR+1;
    }
    /* Loop over antennas making one and two baseline combinations
       between each antenna and the reference.  Then FFT to get
       initial gain, delay and rate from the amp peak of the FFT. */
    for (iAnt=1; iAnt<=in->maxAnt; iAnt++) {
      if (iAnt!=refAnt) {
	/* not reference ant - solve */
	test = initAntSolve (in, iAnt, refAnt, err);
	good = good || test;
      } /* end solve */
    } /* end loop over antennas */
    in->refAnt = in->refAnts[iref];  /* In case this is it */
    if (good) break;
    iref++;    /* Try next */
  } /* end loop over reference antennas */

  /* Now stack data using weights from fitting */
  for (iAnt=1; iAnt<=in->maxAnt; iAnt++) {
    if (iAnt==in->refAnt) continue;
      stackAnt (in, iAnt, err);
  }

} /* end initSolve */

/**
 * Final least squares solutions.
 * Fits single phase, delay, and rate per poln/IF
 * Start off with solutions from initSolve
 * Requires GSL for least squares solution.
 * Model calsulationd potentially multi threaded
 * \param in   The structure with data and solutions.
 * \param err  Error/message stack, returns if error.
 */
static void 
finalSolve (ObitUVGSolveWB *in, ObitErr *err)
{
#if HAVE_GSL==1  /* GSL stuff */
  const gsl_multifit_fdfsolver_type *T;
  gsl_multifit_fdfsolver* solver;
  gsl_multifit_function_fdf func;
  gsl_vector_view coef;
  int status, iter, ncoef, ndata, iAnt, iIF, jIF, iPoln, ip, jp, kp;
  olong loIF, hiIF;
  olong nval, offset, np, iFreq;
  double *coef_init=NULL;
  ofloat phase, delay, disp, IFphase, gnorm, twopi=2.0*G_PI;
  ofloat fblank = ObitMagicF();
  gboolean noData; 
  olong nPoln;
  /*gchar *routine = "finalSolve";*/

  if (err->error) return;  /* Prior error? */
  /* Number of coefficients to solve for 
     phase, delay, for n-1 antennas */
  np = 2;    /* number of parameters per antenna */
  ncoef = np; 
  nval  = in->numIF * in->numPoln;  /* solution entries per antenna */
  if (in->scanData->avgIF) {
    ndata = in->scanData->BLData[0]->numIF * in->scanData->BLData[0]->numFreq;  
    loIF = 0;
    hiIF = 1;
  } else { /* IF at a time */
    ndata = in->scanData->BLData[0]->numFreq;  
    loIF = 0;
    hiIF = in->numIF;
  }
 
  T = gsl_multifit_fdfsolver_lmsder;
  /* T = gsl_multifit_fdfsolver_lmder; Faster */
  solver = gsl_multifit_fdfsolver_alloc(T, ndata, ncoef);

  /* Create /fill function structure */
  func.f   = &fringeFitFunc;      /* Compute function */
  func.df  = &fringeFitJacob;     /* Compute Jacobian (derivative matrix) */
  func.fdf = &fringeFitFuncJacob; /* Compute both function and derivatives */
  func.n = ndata;                 /* Number of data points */
  func.p = ncoef;                 /* number of parameters */
  func.params = in;               /* Object */

  /* Zero initial reference antenna weights */	
  nPoln = in->scanData->BLData[0]->numPoln;
  for (iPoln=0; iPoln<nPoln; iPoln++) {
      for (jIF=loIF; jIF<hiIF; jIF++) {
	ip     = iPoln*in->numIF + jIF;
	in->antWeight[(in->refAnt-1)*nval+ip] = 0.0;
      }
  }

  /* Outer loop over polarization */
  if (in->scanData->avgPoln) nPoln = 1;
  for (iPoln=0; iPoln<nPoln; iPoln++) {
    
    /* Loop over non reference antennas averaging IF values */
    for (iAnt=1; iAnt<=in->maxAnt; iAnt++) {
      if (iAnt==in->refAnt) continue;  /* Not reference antenna */

      /* Loop over IF */
      for (jIF=loIF; jIF<hiIF; jIF++) {
      
	/* Does it have fringes? */
	offset = (iAnt-1)*nval;
	ip     = iPoln*in->numIF;
	if (!in->scanData->avgIF) offset += jIF;
	if (in->antWeight[offset+ip]<in->minSNR) continue;

	/* initial guesses */
	noData = FALSE;
	jp = 0;
	if (coef_init==NULL) coef_init = g_malloc(ncoef*sizeof(double));
	if (in->scanData->avgIF) {
	  /* Get initial value from fitting IF solutions - skip if no data */
	  if (!fitIFData (in, iAnt, iPoln, &phase, &delay, &disp, err)) {
	    noData = TRUE;
	    continue;
	  } /* end no data */
	  /* Set initial values */
	  coef_init[jp++] = (double)phase;
	  coef_init[jp++] = (double)delay;
	} else { /* single IF */
	  ip = iPoln*in->numIF;
	  offset = (iAnt-1)*nval + jIF;
	  coef_init[jp++] = (double)atan2(in->antGain[2*(offset+ip)+1], in->antGain[2*(offset+ip)]);
	  coef_init[jp++] = (double)in->antDelay[offset+ip]*1.0e9;
	}

	/* Create coeffient vector */
	coef = gsl_vector_view_array (coef_init, ncoef);
      
	/* tell what's being solved for */
	in->scanData->curAnt  = iAnt;
	in->scanData->curPoln = iPoln;
	in->scanData->curIF   = jIF;
	in->scanData->avgIF   = in->scanData->avgIF;
	in->scanData->numChan = in->scanData->BLData[0]->numFreq;
	
	/* Set up solver */
	gsl_multifit_fdfsolver_set(solver, &func, &coef.vector);
	
	/* tell initial values */
	if (in->prtLv>=4) {
	  Obit_log_error(err, OBIT_InfoErr, 
			 "ant %d poln %d IF %d initial guess phase %lf delay %lf  RMS %f", 
			 iAnt, iPoln, jIF, coef_init[0]*57.296, coef_init[1], in->scanData->RMSRes*57.296);
	  ObitErrLog(err);
	}

	/* ready to rumble */
	iter = 0;
	do {
	  iter++;
	  status = gsl_multifit_fdfsolver_iterate(solver);
	  if ((status!=GSL_CONTINUE) && (status!=GSL_SUCCESS)) {/* problem? */
	    /*Obit_log_error(err, OBIT_InfoWarn, "%s: Solver status %d %s", 
	      routine,status,gsl_strerror(status));*/
	    /*break;*/
	  }
	  /* convergence test */
	  status = gsl_multifit_test_delta(solver->dx, solver->x, 1.0e-5, 1.0e-5);
	  
	  /* Diagnostic messages */
	  if (in->prtLv>=4) {
	    gnorm = (ofloat)gsl_blas_dnrm2(solver->f);
	    phase  = (ofloat)gsl_vector_get(solver->x, 0);
	    delay  = (ofloat)gsl_vector_get(solver->x, 1);
	    Obit_log_error(err, OBIT_InfoErr, 
			   "ant %d poln %d IF %d Iteration %d phase %f delay %f gradient norm %f RMS %f", 
			   iAnt, iPoln, jIF, iter, phase*57.296, delay, gnorm, in->scanData->RMSRes*57.296);
	    ObitErrLog(err);
	  }
	}
	while ((status == GSL_CONTINUE) && (iter< 100));
	
	/* Get fitted values */
	jp = 0;
	iFreq = in->scanData->BLData[0]->numFreq/2;  /* Center or edge??? */
	iFreq = 0;
	
	phase  = (ofloat)gsl_vector_get(solver->x, jp++);
	phase  = fmodf(phase, twopi);
	delay  = (ofloat)gsl_vector_get(solver->x, jp++)*1.0e-9;
	offset = (iAnt-1)*nval;
	ip = kp = iPoln*in->numIF;
	/* Loop over if averaging IF */
	if (in->scanData->avgIF) {
	  for (iIF=0; iIF<in->numIF; iIF++) {
	    /* Phase of this IF */
	    IFphase = phase + delay * twopi * in->scanData->BLData[0]->dFreq[iIF][iFreq];
	    /* Have data? */
	    if (noData) {
	      in->antGain[2*(offset+ip)]   =  fblank;
	      in->antGain[2*(offset+ip)+1] =  fblank;
	      in->antDelay[offset+kp]      =  fblank;
	    } else { /* OK */
	      in->antGain[2*(offset+ip)]   =  cos(IFphase);
	      in->antGain[2*(offset+ip)+1] =  sin(IFphase);
	      in->antDelay[offset+kp]      =  delay;
	      in->antWeight[offset+kp]     =  1.0/in->scanData->RMSRes;
	      /* Reference antenna SNR the highest of any */
	      in->antWeight[(in->refAnt-1)*nval+kp] = 
		MAX (in->antWeight[(in->refAnt-1)*nval+kp],in->antWeight[offset+kp] );
	    }
	    ip++;
	    kp++;
	  } /* end inner IF loop */
	} else { 
	    offset += jIF;
	    /* Have data? */
	    if (noData) {
	      in->antGain[2*(offset+ip)]   =  fblank;
	      in->antGain[2*(offset+ip)+1] =  fblank;
	      in->antDelay[offset+kp]      =  fblank;
	    } else { /* OK */
	      in->antGain[2*(offset+ip)]   =  cos(phase);
	      in->antGain[2*(offset+ip)+1] =  sin(phase);
	      in->antDelay[offset+kp]      =  delay;
	      in->antWeight[offset+kp]     =  1.0/in->scanData->RMSRes;
	      /* Reference antenna SNR the highest of any */
	      in->antWeight[(in->refAnt-1)*nval+kp+jIF] = 
		MAX (in->antWeight[(in->refAnt-1)*nval+kp+jIF],in->antWeight[offset+kp] );
	    }	  
	}  /* end single IF */
      } /* end outer IF loop */
    } /* end iAnt loop */
  } /* end poln loop */
  
  /* If averaging poln solutions, copy results */
  if ((in->scanData->avgPoln) && (in->scanData->BLData[0]->numPoln>1)) {
    nval = in->numIF * in->numPoln;  /* entries per antenna */
    
    for (iAnt=1; iAnt<=in->maxAnt; iAnt++) { /* Loop over antennas */
      offset = (iAnt-1)*nval;          /* This antennas offset on solutions */
      ip = in->scanData->BLData[0]->numIF;
      jp = 0;
      /* Loop over IF */
      for (iIF=0; iIF<in->scanData->BLData[0]->numIF; iIF++) {
	in->antGain[2*(offset+ip)]   = in->antGain[2*(offset+jp)];
	in->antGain[2*(offset+ip)+1] = in->antGain[2*(offset+jp)+1];
	in->antDelay[offset+ip]      = in->antDelay[offset+jp];
	in->antDisp[offset+ip]       = in->antDisp[offset+jp];
	in->antWeight[offset+ip]     = in->antWeight[offset+jp];
	jp++;  ip++;
      } /* End IF Loop */
    } /* end antenna loop */
  } /* End copy solutions to other polarization if averaging */

  /* cleanup */
  gsl_multifit_fdfsolver_free(solver);
  if (coef_init) g_free(coef_init);
  return;

#else  /* No GSL - stubb */
  gchar *routine = "finalSolve";
  Obit_log_error(err, OBIT_Error, 
		 "%s: GSL not available - cannot do fit", 
		     routine);
  return;
#endif /* GSL stuff */

} /* end finalSolve */

/**
 * Make initial guess for big least squares from IF solutions for one antenna
 * \param in      Solver data structure
 * \param iAnt    Antenna number
 * \param iPoln   Polarization
 * \param phase   [out] fitted phase at reference freq
 * \param delay   [out] fitted group delay in nsec
 * \param disp    [out] fitted dispersion
 * \param err    Error/message stack, returns if error.
 * \return True if some data
 */
static gboolean
fitIFData (ObitUVGSolveWB *in, olong iAnt, olong iPoln, 
	   ofloat *phase, ofloat *delay, ofloat *disp, ObitErr* err)
{
  gboolean out = FALSE;
  olong ip, jp, i, kndx, offset, iIF, nval;
  ofloat sumPhs, sumDelay, modDly, modPhs, sumD1, dph, dly, twopi=2.0*G_PI;
  olong count, count1, iturn;
#if HAVE_GSL==1  /* GSL stuff */
  gsl_vector_view coef;
  int status, iter;
#endif /* GSL stuff */
  /*gchar *routine = "fitIFData";*/

  /* Error check */
  if (err->error) return out;

  /* Check if any data */
  kndx = (iAnt-1) + iPoln*in->scanData->maxAnt;
  sumD1 = ObitFArraySum (in->scanData->antStackWt[kndx]);
  out = sumD1 > 0.0;
  if (!out) return out;

  ip = iPoln*in->numIF;
  sumDelay = 0.0; count = 0;
  nval = in->numIF * in->numPoln;  /* entries per antenna */
  offset = (iAnt-1)*nval;
  /* Loop over IF */
  for (iIF=0; iIF<in->numIF; iIF++) {
    sumDelay  += in->antDelay[offset+ip++];
    count++;
  } /* end IF loop */
  
  /* Initialize output */
  ip = iPoln*in->numIF;
  *phase = atan2(in->antGain[2*(offset+ip)+1], in->antGain[2*(offset+ip)]);
  *delay = 1.0e9*sumDelay/count;  /* In nsec */
  *disp  = 0.0;
  
  /* With only one frequency, this is the best you can do */
  if (in->ndataCoarse<=1) return out;

  /* Set data */
  for (i=0; i<in->ndataCoarse; i++) {
    in->coarseData[i] = atan2(in->antGain[2*(offset+ip)+1], in->antGain[2*(offset+ip)]);
    ip++;
  }
  
  /* Refine delay using First differences */
  modDly = *delay;
  count1 = 0; sumD1 = 0.0;
  for (i=1; i<in->ndataCoarse; i++) {
    count1++;
    /* Measured phase difference */
    dph = (in->coarseData[i] - in->coarseData[i-1]);
    /* Model phase difference */
    modPhs = twopi*modDly*(in->coarseFreq[i] - in->coarseFreq[i-1]);
    /* Add turns of phase to best agreement with prior model */
    iturn = (olong)(0.5 + ((modPhs - dph) / twopi));
    dph += iturn*twopi;
    /* Nearest half turn */
    if (dph>G_PI) dph -= twopi;
    else if (dph<-G_PI) dph += twopi;
    dly = dph / (twopi*(in->coarseFreq[i] - in->coarseFreq[i-1]));
    sumD1 += dly;
  }
  *delay = sumD1/count1;

  /* Rerefine delay using differences by two */
  if (in->ndataCoarse>2) {
    modDly = *delay;
    count1 = 0; sumD1 = 0.0;
    for (i=2; i<in->ndataCoarse; i++) {
      count1++;
      /* Measured phase difference */
      dph = (in->coarseData[i] - in->coarseData[i-2]);
      /* Model phase difference */
      modPhs = twopi*modDly*(in->coarseFreq[i] - in->coarseFreq[i-2]);
      /* Add turns of phase to best agreement with prior model */
      iturn = (olong)(0.5 + ((modPhs - dph) / twopi));
      dph += iturn*twopi;
      /* Nearest half turn */
      if (dph>G_PI) dph -= twopi;
      else if (dph<-G_PI) dph += twopi;
      dly = dph / (twopi*(in->coarseFreq[i] - in->coarseFreq[i-2]));
      sumD1 += dly;
    }
    *delay = sumD1/count1;
  }

  /* Refine phase using delay estimate */
  modDly = *delay;
  sumPhs = 0.0; count = 0;
  for (i=0; i<in->ndataCoarse; i++) {
    modPhs = twopi*modDly*in->coarseFreq[i];
    dph = in->coarseData[i] - modPhs;
    /* Nearest half turn */
    if (dph>G_PI) dph -= twopi;
    else if (dph<-G_PI) dph += twopi;
    sumPhs += dph;
    count++;
  }
  *phase = sumPhs / count;

#if HAVE_GSL==1  /* GSL stuff */

  jp = 0;
  in->coarse_coef_init[jp++] = (double)*phase;
  in->coarse_coef_init[jp++] = (double)*delay;
  coef = gsl_vector_view_array (in->coarse_coef_init, in->ncoefCoarse);
  
  /* Set up solver */
  gsl_multifit_fdfsolver_set(in->myCoarseSolver, in->myCoarseFunc, &coef.vector);
    
    /* ready to rumble */
    iter = 0;
    do {
      iter++;
      status = gsl_multifit_fdfsolver_iterate(in->myCoarseSolver);
      if ((status!=GSL_CONTINUE) && (status!=GSL_SUCCESS)) {/* problem? */
	/*Obit_log_error(err, OBIT_Error, "%s: Solver status %d %s", 
	  routine,status,gsl_strerror(status));
	  break; */
      }
      /* convergence test */
      status = gsl_multifit_test_delta(in->myCoarseSolver->dx, in->myCoarseSolver->x, 1.0e-4, 1.0e-4);
    }
    while ((status == GSL_CONTINUE) && (iter< 100));
    
    /* Get fitted values */
    jp = 0;
    *phase  = (ofloat)gsl_vector_get(in->myCoarseSolver->x, jp++);
    *delay  = (ofloat)gsl_vector_get(in->myCoarseSolver->x, jp++);
#endif /* GSL stuff */
 
    /* Flip sign */
    *phase = -(*phase);

    return out;
}  /* End fitIFData */

/**
 * Determine the antenna with the highest weight in the current data
 * \param in   The structure with scan data, returns -1 if not fully defined.
 * \param err  Error/message stack, returns if error.
 * \return (1-rel) antenna number with highest summed weight
 */
static olong 
GetRefAnt (ScanData *scanData, ObitErr *err) {
  olong refAnt=-1;
  ofloat *sum=NULL, maxSum;
  olong iBase, iIF, iPoln, i, ant1, ant2, nIF, nPoln, ip, best;

  if (err->error) return refAnt;
  if (scanData==NULL) return refAnt;
  if (scanData->BLData==NULL) return refAnt;

  /* Create sum array */
  sum = g_malloc0((scanData->maxAnt+1)*sizeof(ofloat));

  /* Sum weights for antennas (use as 1-rel to avoid confusion) */
  for (iBase=0; iBase<scanData->numBase; iBase++) {  /* Baseline loop */
    ant1  = scanData->BLData[iBase]->ant1;
    ant2  = scanData->BLData[iBase]->ant2;
    nIF   = scanData->BLData[iBase]->numIF;
    nPoln = scanData->BLData[iBase]->numPoln;
    ip    = 0;
    /* Loop over Poln */
    for (iPoln=0; iPoln<nPoln; iPoln++) {
    /* Loop over IF */
      for (iIF=0; iIF<nIF; iIF++) {
	sum[ant1] += MAX (0.0, scanData->BLData[iBase]->WtPolnIF[ip]);
	sum[ant2] += MAX (0.0, scanData->BLData[iBase]->WtPolnIF[ip]);
	ip++;
      } /* end IF loop */
    } /* end poln loop */
  } /* end loop over baselines accumulating */

  /* Find best */
  best = -1; maxSum = -1.0e20;
  for (i=1; i<=scanData->maxAnt; i++) {
    if (sum[i]>maxSum) {
      maxSum = sum[i];
      best   = i;
    }
  }
  refAnt = best;

  /* Cleanup */
  if (sum) g_free(sum);
  return refAnt;
} /* end GetRefAnt */

/**
 * Initial guess at solutions based on FFT on baselines to reference antenna
 * Accumulates one and two baseline connections between iAnt and the 
 * reference antenna and FFTs to determine the value at the peak.
 * This is used to determine gain, delay and rate per IF/Poln
 * \param in     The structure with data and solutions.
 * \param iAnt   which Antenna (1-rel)
 * \param refAnt reference Antenna (1-rel)
 * \param err    Error/message stack, returns if error.
 * \return TRUE if fringes found on at least some IF/Poln
 */
static gboolean 
initAntSolve (ObitUVGSolveWB *in, olong iAnt, olong refAnt, ObitErr *err)
{
  olong naxis[1], maxis[1], iBase, jBase, iPoln, jAnt, iIF, ip, jp, nPoln;
  olong iBase1, iBase2, offset, nval, icend, refChan, nFreq2, k, kndx, kp;
  olong count;
  ofloat cmplx[2] = {0.0,0.0}, delTau, iNorm, pval[2], ppos, dph, dre, dim, tre, tim;
  ofloat peak, fblank = ObitMagicF();
  ofloat amp, sumw, sumww, ph;
  gboolean good=FALSE, OK=FALSE;
  if (err->error) return OK;  /* Prior error? */

  /* Need to create work arrays? */
  naxis[0] = in->scanData->BLData[0]->numFreq;
  iNorm = 1.0 / naxis[0];   /* Normalization factor */
  /* Make FFT friendly, then quadruple to zero pad */
  naxis[0] = in->FFTOverSamp*ObitFFTSuggestSize(naxis[0]);

  nFreq2 =  in->scanData->BLData[0]->numFreq/2;     /* Center frequency channel */
  refChan = (olong)(in->scanData->refChan+0.5) - 1; /* Ref. channel */

  if (in->cWork1==NULL) in->cWork1 = ObitCArrayCreate ("CWork1", 1, naxis);
  else in->cWork1 = ObitCArrayRealloc (in->cWork1, 1, naxis);
  if (in->cWork2==NULL) in->cWork2 = ObitCArrayCreate ("CWork2", 1, naxis);
  else in->cWork2 = ObitCArrayRealloc (in->cWork2, 1, naxis);
  if (in->cWork3==NULL) in->cWork3 = ObitCArrayCreate ("CWork2", 1, naxis);
  else in->cWork3 = ObitCArrayRealloc (in->cWork3, 1, naxis);
  if (in->fWork1==NULL) in->fWork1 = ObitFArrayCreate ("FWork1", 1, naxis);
  else in->fWork1 = ObitFArrayRealloc (in->fWork1, 1, naxis);
  if (in->fWork2==NULL) in->fWork2 = ObitFArrayCreate ("FWork2", 1, naxis);
  else in->fWorkWt2 = ObitFArrayRealloc (in->fWorkWt2, 1, naxis);
  if (in->fWorkWt2==NULL) in->fWorkWt2 = ObitFArrayCreate ("FWorkWt2", 1, naxis);
  else in->fWork2 = ObitFArrayRealloc (in->fWork2, 1, naxis);
  maxis[0] = in->scanData->BLData[0]->numFreq;
  if (in->fWork3==NULL) in->fWork3 = ObitFArrayCreate ("FWork3", 1, maxis);
  else in->fWork3 = ObitFArrayRealloc (in->fWork3, 1, maxis);

  /* Stuff for finding solution */
  icend = naxis[0]/2;  /* Center delay pixel */
  delTau  = 1.0/(naxis[0]*in->scanData->freqAvg);   /* Delay increment */

  /* Need to create FFT? */
  if (in->myFFT==NULL) 
    in->myFFT = 
      newObitFFT ("Fringe", OBIT_FFT_Forward, OBIT_FFT_FullComplex, 1, naxis);
  
  nval = in->numIF * in->numPoln;  /* entries per antenna */
  offset = (iAnt-1)*nval;          /* This antennas offset on solutions */

  /* How many polarizations to solve for? */
  if (in->scanData->avgPoln)  nPoln = 1;
  else nPoln = in->scanData->BLData[0]->numPoln;

  /* Loop over poln */
  ip = 0;
  for (iPoln=0; iPoln<nPoln; iPoln++) {
    /* Loop over IF */
    for (iIF=0; iIF<in->scanData->BLData[0]->numIF; iIF++) {

      /* Zero accumulators (cWorkn, fWorkn) */
      ObitCArrayFill(in->cWork1,   cmplx);
      ObitFArrayFill(in->fWork1,   0.0);
      ObitCArrayFill(in->cWork2,   cmplx);
      ObitFArrayFill(in->fWork2,   0.0);
      ObitCArrayFill(in->cWork3,   cmplx);
      ObitFArrayFill(in->fWork3,   0.0);
      ObitFArrayFill(in->fWorkWt2, 0.0);
      count = 0;

      /* One baseline combination */
      iBase = -1;
      for (jBase=0; jBase<in->scanData->numBase; jBase++) {
	if (((in->scanData->BLData[jBase]->ant1==iAnt) && 
	     (in->scanData->BLData[jBase]->ant2==refAnt)) ||
	    ((in->scanData->BLData[jBase]->ant2==iAnt) && 
	     (in->scanData->BLData[jBase]->ant1==refAnt))) {
	  /* found it  */
	  iBase = jBase;
	  break;
	}
      }

      /* find it? Any data? */
      if ((iBase>=0) && (in->scanData->BLData[iBase]->WtPolnIF[ip]>0.0)) { 
	/* Promote phases to complex, zero pad */
	ObitCArrayFSinCos (in->scanData->BLData[iBase]->phArray[ip], in->cWork2);
	/* Conjugate if reference antenna has lower number */
	if (refAnt<iAnt) ObitCArrayConjg(in->cWork2);
	/* Add to end of accumulator */
	ObitCArrayAdd (in->cWork1,  in->cWork2, in->cWork1);
	/* Multiply by Weights */
	ObitCArrayFMulEnd (in->cWork1, in->scanData->BLData[iBase]->wtArray[ip],  in->cWork1);
	/* Accumulate weights */
	ObitFArrayAddEnd (in->fWork1,  in->scanData->BLData[iBase]->wtArray[ip],  in->fWork1);
	/* Accumulate weights squares, count */
	ObitFArrayAddEnd2 (in->fWorkWt2,  in->scanData->BLData[iBase]->wtArray[ip],  
			   in->fWorkWt2,  &count);
	/* Double the weights */
        ObitCArraySMul (in->cWork1, 2.0);
        ObitFArraySMul (in->fWork1, 2.0);
        ObitFArraySMul (in->fWorkWt2, 4.0);

	/* Zero weights where data blanked */
	ObitFArrayZeroBlank (in->fWork1, in->scanData->BLData[iBase]->phArray[ip], in->fWork1);
	ObitFArrayZeroBlank (in->fWorkWt2, in->scanData->BLData[iBase]->phArray[ip], in->fWorkWt2);
      }

      /* Two baseline combinations - search for antenna with data on both baselines 
	 to iAnt and refAnt */
      if (in->scanData->doTwo) {
	iBase1 = -1; iBase2 = -1;
	for (jAnt = 1; jAnt<=in->scanData->maxAnt; jAnt++) {
	  /* Ignore iAnt and refAnt */
	  if ((jAnt==iAnt) || (jAnt==refAnt)) continue;
	  
	  /* First find baseline 1 between iAnt and jAnt */
	  for (jBase=0; jBase<in->scanData->numBase; jBase++) {
	    if (((in->scanData->BLData[jBase]->ant1==iAnt) && 
		 (in->scanData->BLData[jBase]->ant2==jAnt)) ||
		((in->scanData->BLData[jBase]->ant2==iAnt) && 
		 (in->scanData->BLData[jBase]->ant1==jAnt))) {
	      /* found it  */
	      iBase1 = jBase;
	      break;
	    }
	  }
	  
	  /* Next find baseline 2 between jAnt and refAnt */
	  for (jBase=0; jBase<in->scanData->numBase; jBase++) {
	    if (((in->scanData->BLData[jBase]->ant1==refAnt) && 
		 (in->scanData->BLData[jBase]->ant2==jAnt)) ||
		((in->scanData->BLData[jBase]->ant2==refAnt) && 
		 (in->scanData->BLData[jBase]->ant1==jAnt))) {
	      /* found it  */
	      iBase2 = jBase;
	      break;
	    }
	  }
	  /* find them? Any data on both? */
	  if (((iBase1>=0) && (in->scanData->BLData[iBase1]->WtPolnIF[ip]>0.0)) &&
	      ((iBase2>=0) && (in->scanData->BLData[iBase2]->WtPolnIF[ip]>0.0))) { 
	    /* Get difference in fWork3 -> cWork2 */
	    /* want baseline 1 - baseline 2 ((ant-other)-(ref-other)) 
	     baseline data will have ant2>ant1 - may need to flip sign */
	    if (refAnt>jAnt) {  /* Negate baseline 2 */
	      if (jAnt<iAnt) {      /* Also negate baseline 1 */
		ObitFArraySub (in->scanData->BLData[iBase2]->phArray[ip], 
			       in->scanData->BLData[iBase1]->phArray[ip], in->fWork3);
	      } else                /* baseline 1 OK */
		ObitFArrayAdd (in->scanData->BLData[iBase1]->phArray[ip], 
			       in->scanData->BLData[iBase2]->phArray[ip], in->fWork3);
	    } else {   /* Baseline 2 OK */
	      if (jAnt<iAnt) {  /* negate baseline 1 */
		ObitFArrayAdd (in->scanData->BLData[iBase2]->phArray[ip], 
			       in->scanData->BLData[iBase1]->phArray[ip], in->fWork3);
		ObitFArrayNeg (in->fWork3);
	      } else {          /* Also baseline 1 OK */
		ObitFArraySub (in->scanData->BLData[iBase1]->phArray[ip], 
			       in->scanData->BLData[iBase2]->phArray[ip], in->fWork3);
	      }
	    } /* end ref vs other */
	    
	    /* Promote phases to complex, zero pad */
	    ObitCArrayFSinCos (in->fWork3, in->cWork2);
	    
	    /* Get harmonic sum Weight array for this combination fWork2 */
	    ObitFArrayHarmAddEnd (in->scanData->BLData[iBase1]->wtArray[ip], 
				  in->scanData->BLData[iBase2]->wtArray[ip], 
				  in->fWork2);
	    
	    /* Accumulate weights squares, count */
	    ObitFArrayHarmAccEnd2 (in->scanData->BLData[iBase1]->wtArray[ip], 
				   in->scanData->BLData[iBase2]->wtArray[ip],
				   in->fWorkWt2, &count);
	    
	    /* Zero weights where data blanked */
	    ObitFArrayZeroBlank (in->fWork2, in->fWork3, in->fWork2);

	    /* Multiply by Weights */
	    ObitCArrayFMul (in->cWork2, in->fWork2 , in->cWork2);
	    
	    /* Accumulate to cWork1, fWork1 */
	    ObitCArrayAdd (in->cWork1, in->cWork2, in->cWork1);
	    ObitFArrayAdd (in->fWork1, in->fWork2, in->fWork1);
	    
	  } /* end add two baseline combination */
	} /* end loop over secondary antennas */
      } /* end if doTwo */

      /* Check for good data */
      good = FALSE;
      kndx = (iAnt-1) + iPoln*in->scanData->maxAnt;
      kp = iIF*in->scanData->BLData[0]->numFreq;
      for (k=0; k<in->scanData->antStackWt[kndx]->naxis[0]; k++) {
	good = good || (in->fWork1->array[k]>0.0);  /* Any good data? */
	kp++;
      }

      /* Zero solutions */
      in->antGain[2*(offset+ip)]   = 0.0;
      in->antGain[2*(offset+ip)+1] = 0.0;
      in->antDelay[offset+ip]      = 0.0;
      in->antDisp[offset+ip]       = 0.0; 
      in->antWeight[offset+ip]     = 0.0;

      /* Any valid data? */
      if (good) {
	/* get weight sums */
	sumw  = ObitFArraySum(in->fWork1);
	sumww = ObitFArraySum(in->fWorkWt2);

	/* Replace zero weights by 1.0 */
	ObitFArrayInClip (in->fWork1, -1.0e-6, 1.0e-6, 1.0);
	
	/* Normalize by sum of weights */
	ObitCArrayFDiv (in->cWork1, in->fWork1, in->cWork1);
	
	/* FFT cWork1 to cWork2 */
	ObitFFTC2C (in->myFFT, in->cWork1, in->cWork2);
	
	/* Get location and value at peak, uses cWork2, fWork2 */
	FindFFTPeak (in, &ppos, pval);
	if ((pval[0]!=fblank) && (pval[1]))
	  peak = sqrt (pval[0]*iNorm*pval[0]*iNorm + pval[1]*iNorm*pval[1]*iNorm);
	else 
	  peak = 0.0;
        /* This is REALLY needed to circumvent a gcc bug */
        if (err->error) fprintf (stderr,"GCC Bug workaround\n");

	/* Get SNR */
	amp   = peak;
	in->antWeight[offset+ip] = decorSNR (amp, sumw, sumww, 1.0, count);
	/* Use the best found as the reference antenna value */
	in->antWeight[(refAnt-1)*nval+ip] = 
	  MAX (in->antWeight[(refAnt-1)*nval+ip], in->antWeight[offset+ip]);
	
	/* Check for OK solution */
	if (in->antWeight[offset+ip]>in->minSNR) {
	  
	  /* Gain is peak value - measured at center of IF */
	  ph = atan2(pval[1], pval[0]);
	  in->antGain[2*(offset+ip)]   = cos(ph);
	  in->antGain[2*(offset+ip)+1] = sin(ph);
	  OK = TRUE;  /* At least something worked */
	  
	  /* Delay and rate from position of peak - 
	     if less than 1/2 of a cell from zero - call it zero 
	     WHY ?? - disable with minus */
	  if (fabs(ppos-icend)<-0.125*in->FFTOverSamp) in->antDelay[offset+ip] = 0.0;
	  else in->antDelay[offset+ip]      = (ppos-icend)*delTau;
	  
	  in->antDisp[offset+ip]       = 0.0;  /* No dispersion here */

	  /* Correct phase to reference channel */
	  dph = -2.0*G_PI*in->antDelay[offset+ip] * 
	    (in->scanData->BLData[0]->dFreq[iIF][nFreq2] - 
	     in->scanData->BLData[0]->dFreq[iIF][refChan]);
	  ObitSinCosCalc (dph, &dim, &dre); /* Sine/cosine */
	  tre = in->antGain[2*(offset+ip)];
	  tim = in->antGain[2*(offset+ip)+1];
	  in->antGain[2*(offset+ip)]   = tre*dre - tim*dim;
	  in->antGain[2*(offset+ip)+1] = tre*dim + tim*dre;
	  
	  /* Adjust for baseline direction - no already done  */
	  if (refAnt>iAnt) {
	    /* Flip phases
	    in->antGain[2*(offset+ip)+1] = -in->antGain[2*(offset+ip)+1];
	    in->antDelay[offset+ip]      = -in->antDelay[offset+ip]; */
	  }

	} /* end sensible fit */
	/* Tell result if requested */
	if (in->prtLv>=3) {
	  Obit_log_error(err, OBIT_InfoErr, 
			 "ant %2d poln %d IF %2d  delay %9.3f nsec SNR %8.1f", 
			 iAnt, iPoln, iIF, in->antDelay[offset+ip]*1.0e9, 
			 in->antWeight[offset+ip]);
	  ObitErrLog(err);
	}
	
      } /* end if good data */
      ip++;
    } /* end loop over IF */
  } /* end loop over Poln */

  /* If averaging poln solutions, copy results */
  if ((in->scanData->avgPoln) && (in->scanData->BLData[0]->numPoln>1)) {
    ip = in->scanData->BLData[0]->numIF;
    jp = 0;
    /* Loop over IF */
    for (iIF=0; iIF<in->scanData->BLData[0]->numIF; iIF++) {
      in->antGain[2*(offset+ip)]   = in->antGain[2*(offset+jp)];
      in->antGain[2*(offset+ip)+1] = in->antGain[2*(offset+jp)+1];
      in->antDelay[offset+ip]      = in->antDelay[offset+jp];
      in->antDisp[offset+ip]       = in->antDisp[offset+jp];
      in->antWeight[offset+ip]     = in->antWeight[offset+jp];
      /* Use the best found as the reference antenna value */
      in->antWeight[(refAnt-1)*nval+ip] = 
	MAX (in->antWeight[(refAnt-1)*nval+ip], in->antWeight[offset+jp]);
     jp++;  ip++;
    }
  }
  return OK;
} /* end  initAntSolve */

/**
 * Stack one and two baseline combination phases between 
 * an antenna and the reference.
 * \param in   The structure with data and solutions.
 * \param iAnt which Antenna (1-rel)
 * \param err  Error/message stack, returns if error.
 */
static void 
stackAnt (ObitUVGSolveWB *in, olong iAnt, ObitErr *err)
{
  olong naxis[1], maxis[1], iBase, jBase, iPoln, jAnt, iIF, ip, nPoln;
  olong iBase1, iBase2, offset, nval, k, kndx, kp;
  ofloat cmplx[2] = {0.0,0.0};
  ofloat antWt;
  gboolean good;
  if (err->error) return;  /* Prior error? */

  /* Need to create work arrays? */
  naxis[0] = in->scanData->BLData[0]->numFreq;
  /* Make FFT friendly, then quadruple to zero pad */
  naxis[0] = in->FFTOverSamp*ObitFFTSuggestSize(naxis[0]);

  if (in->cWork1==NULL) in->cWork1 = ObitCArrayCreate ("CWork1", 1, naxis);
  else in->cWork1 = ObitCArrayRealloc (in->cWork1, 1, naxis);
  if (in->cWork2==NULL) in->cWork2 = ObitCArrayCreate ("CWork2", 1, naxis);
  else in->cWork2 = ObitCArrayRealloc (in->cWork2, 1, naxis);
  if (in->cWork3==NULL) in->cWork3 = ObitCArrayCreate ("CWork2", 1, naxis);
  else in->cWork3 = ObitCArrayRealloc (in->cWork3, 1, naxis);
  if (in->fWork1==NULL) in->fWork1 = ObitFArrayCreate ("FWork1", 1, naxis);
  else in->fWork1 = ObitFArrayRealloc (in->fWork1, 1, naxis);
  if (in->fWork2==NULL) in->fWork2 = ObitFArrayCreate ("FWork2", 1, naxis);
  else in->fWorkWt2 = ObitFArrayRealloc (in->fWorkWt2, 1, naxis);
  if (in->fWorkWt2==NULL) in->fWorkWt2 = ObitFArrayCreate ("FWorkWt2", 1, naxis);
  else in->fWork2 = ObitFArrayRealloc (in->fWork2, 1, naxis);
  maxis[0] = in->scanData->BLData[0]->numFreq;
  if (in->fWork3==NULL) in->fWork3 = ObitFArrayCreate ("FWork3", 1, maxis);
  else in->fWork3 = ObitFArrayRealloc (in->fWork3, 1, maxis);

  nval = in->numIF * in->numPoln;  /* entries per antenna */
  offset = (iAnt-1)*nval;          /* This antennas offset on solutions */

  /* How many polarizations? */
  if (in->scanData->avgPoln)  nPoln = 1;
  else nPoln = in->scanData->BLData[0]->numPoln;

  /* Loop over poln */
  ip = 0;
  for (iPoln=0; iPoln<nPoln; iPoln++) {
    /* Loop over IF */
    for (iIF=0; iIF<in->scanData->BLData[0]->numIF; iIF++) {

      /* Zero accumulators (cWorkn, fWorkn) */
      ObitCArrayFill(in->cWork1,   cmplx);
      ObitFArrayFill(in->fWork1,   0.0);
      ObitCArrayFill(in->cWork2,   cmplx);
      ObitFArrayFill(in->fWork2,   0.0);
      ObitCArrayFill(in->cWork3,   cmplx);
      ObitFArrayFill(in->fWork3,   0.0);
      ObitFArrayFill(in->fWorkWt2, 0.0);

      /* One baseline combination */
      iBase = -1;
      for (jBase=0; jBase<in->scanData->numBase; jBase++) {
	if (((in->scanData->BLData[jBase]->ant1==iAnt) && 
	     (in->scanData->BLData[jBase]->ant2==in->refAnt)) ||
	    ((in->scanData->BLData[jBase]->ant2==iAnt) && 
	     (in->scanData->BLData[jBase]->ant1==in->refAnt))) {
	  /* found it 
	     Where there Fringes on this baseline? */
	  if ((in->antWeight[(in->refAnt-1)*nval+ip]>in->minSNR) && 
	      (in->antWeight[(iAnt-1)*nval+ip]>in->minSNR))
	    iBase = jBase;
	  break;
	}
      }

      /* find it? Any data? */
      if ((iBase>=0) && (in->scanData->BLData[iBase]->WtPolnIF[ip]>0.0)) { 
	/* Promote phases to complex, zero pad */
	ObitCArrayFSinCos (in->scanData->BLData[iBase]->phArray[ip], in->cWork2);
	/* Conjugate if reference antenna has lower number */
	if (in->refAnt<iAnt) ObitCArrayConjg(in->cWork2);
	/* Add to end of accumulator */
	ObitCArrayAdd (in->cWork1,  in->cWork2, in->cWork1);
	/* Get data weights */
	ObitFArrayAddEnd (in->fWork1,  in->scanData->BLData[iBase]->wtArray[ip],  in->fWork1);
	/* Multiply data weight by product of antenna weights * 2 */
	/*antWt = 2.0 * sqrt(in->antWeight[(in->refAnt-1)*nval+ip] * 
	  in->antWeight[(iAnt-1)*nval+ip]);
	  Works better without weighting */
	antWt = 2.0;
	ObitFArraySMul (in->fWork1, antWt);
	
	/* Zero weights where data blanked */
	ObitFArrayZeroBlank (in->fWork1, in->scanData->BLData[iBase]->phArray[ip], in->fWork1);

	/* Multiply by  weights */
	ObitCArrayFMulEnd (in->cWork1, in->fWork1,  in->cWork1);
      } /* end of any data on baseline */

      /* Two baseline combinations - search for antenna with data on both baselines 
	 to iAnt and in->refAnt */
      if (in->scanData->doTwo) {
	iBase1 = -1; iBase2 = -1;
	for (jAnt = 1; jAnt<=in->scanData->maxAnt; jAnt++) {
	  /* Ignore iAnt and refAnt */
	  if ((jAnt==iAnt) || (jAnt==in->refAnt)) continue;
	  
	  /* First find baseline 1 between iAnt and jAnt */
	  for (jBase=0; jBase<in->scanData->numBase; jBase++) {
	    if (((in->scanData->BLData[jBase]->ant1==iAnt) && 
		 (in->scanData->BLData[jBase]->ant2==jAnt)) ||
		((in->scanData->BLData[jBase]->ant2==iAnt) && 
		 (in->scanData->BLData[jBase]->ant1==jAnt))) {
	      /* found it  
		 Where there Fringes on this baseline? */
	      if ((in->antWeight[(iAnt-1)*nval+ip]>in->minSNR) && 
		  (in->antWeight[(jAnt-1)*nval+ip]>in->minSNR))
		iBase1 = jBase;
	      break;
	    }
	  }
	  
	  /* Next find baseline 2 between jAnt and in->refAnt */
	  for (jBase=0; jBase<in->scanData->numBase; jBase++) {
	    if (((in->scanData->BLData[jBase]->ant1==in->refAnt) && 
		 (in->scanData->BLData[jBase]->ant2==jAnt)) ||
		((in->scanData->BLData[jBase]->ant2==in->refAnt) && 
		 (in->scanData->BLData[jBase]->ant1==jAnt))) {
	      /* found it  
		 Where there Fringes on this baseline? */
	      if ((in->antWeight[(in->refAnt-1)*nval+ip]>in->minSNR) && 
		  (in->antWeight[(jAnt-1)*nval+ip]>in->minSNR))
		iBase2 = jBase;
	      break;
	    }
	  }
	  /* find them? Any data on both? */
	  if (((iBase1>=0) && (in->scanData->BLData[iBase1]->WtPolnIF[ip]>0.0)) &&
	      ((iBase2>=0) && (in->scanData->BLData[iBase2]->WtPolnIF[ip]>0.0))) { 
	    /* Get difference in fWork3 -> cWork2 */
	    /* want baseline 1 - baseline 2 ((ant-other)-(ref-other)) 
	     baseline data will have ant2>ant1 - may need to flip sign */
	    if (in->refAnt>jAnt) {  /* Negate baseline 2 */
	      if (jAnt<iAnt) {      /* Also negate baseline 1 */
		ObitFArraySub (in->scanData->BLData[iBase2]->phArray[ip], 
			       in->scanData->BLData[iBase1]->phArray[ip], in->fWork3);
	      } else                /* baseline 1 OK */
		ObitFArrayAdd (in->scanData->BLData[iBase1]->phArray[ip], 
			       in->scanData->BLData[iBase2]->phArray[ip], in->fWork3);
	    } else {   /* Baseline 2 OK */
	      if (jAnt<iAnt) {  /* negate baseline 1 */
		ObitFArrayAdd (in->scanData->BLData[iBase2]->phArray[ip], 
			       in->scanData->BLData[iBase1]->phArray[ip], in->fWork3);
		ObitFArrayNeg (in->fWork3);
	      } else {          /* Also baseline 1 OK */
		ObitFArraySub (in->scanData->BLData[iBase1]->phArray[ip], 
			       in->scanData->BLData[iBase2]->phArray[ip], in->fWork3);
	      }
	    } /* end ref vs other */
	    
	    /* Promote phases to complex, zero pad */
	    ObitCArrayFSinCos (in->fWork3, in->cWork2);
	    
	    /* Baseline 1 weights - multiply by product of antenna weights 
	       antWt = sqrt(in->antWeight[(jAnt-1)*nval+ip] * in->antWeight[(iAnt-1)*nval+ip]);*/
	    antWt = 1.0;
	    ObitFArrayFill(in->fWork2,   0.0);
	    ObitFArrayAddEnd (in->fWork2,in->scanData->BLData[iBase1]->wtArray[ip], in->fWork2);
	    /* Multiply data weight by product of antenna weights * 2 */
	    ObitFArraySMul (in->fWork2, antWt);

	    /* Baseline 2 weights - multiply by product of antenna weights 
	       antWt = sqrt(in->antWeight[(in->refAnt-1)*nval+ip] * in->antWeight[(jAnt-1)*nval+ip]);*/
	    antWt = 1.0;
	    ObitFArrayFill(in->fWorkWt2, 0.0);
	    ObitFArrayAddEnd (in->fWorkWt2,in->scanData->BLData[iBase1]->wtArray[ip], in->fWorkWt2);
	    /* Multiply data weight by product of antenna weights * 2 */
	    ObitFArraySMul (in->fWork2, antWt);

	    /* Get harmonic sum Weight array for this combination fWork2 */
	    ObitFArrayHarmAddEnd (in->fWork2, in->fWorkWt2, in->fWork2);

	    /* Zero weights where data blanked */
	    ObitFArrayZeroBlank (in->fWork2, in->fWork3, in->fWork2);

	    /* Multiply phasors by Weights */
	    ObitCArrayFMul (in->cWork2, in->fWork2, in->cWork2);
	    
	    /* Accumulate to cWork1, fWork1 */
	    ObitCArrayAdd (in->cWork1, in->cWork2, in->cWork1);
	    ObitFArrayAdd (in->fWork1, in->fWork2, in->fWork1);
	    
	  } /* end add two baseline combination */
	} /* end loop over secondary antennas */
      } /* end if doTwo */

      /* Copy stacked phases/weights to in->antStackPh */
      good = FALSE;
      kndx = (iAnt-1) + iPoln*in->scanData->maxAnt;
      kp = iIF*in->scanData->BLData[0]->numFreq;
      for (k=0; k<in->scanData->antStackWt[kndx]->naxis[0]; k++) {
	in->scanData->antStackPh[kndx]->array[kp] = 
	  atan2(in->cWork1->array[1+2*k], in->cWork1->array[2*k]);
	in->scanData->antStackWt[kndx]->array[kp] = in->fWork1->array[k];
	good = good || (in->fWork1->array[k]>0.0);  /* Any good data? */
	kp++;
      }
      ip++; /* next Poln/IF */
    } /* End IF Loop */
  } /* end Poln loop */
} /* end  stackAnt */

/**
 * Create Baseline structure
 * \param ant1    First antenna number (1-rel)
 * \param ant2    Second antenna number (1-rel)
 * \param numIF   Number of IFs
 * \param numPoln Number of poln 
 * \param avgPoln If TRUE, agerage polarisations (one poln data array)
 * \param numFreq Number of frequency averages per IF
 * \return pointer to the new structure.
 */
static BaselineData* MakeBLData (olong ant1, olong ant2, 
				 olong numIF, olong numPoln,
				 gboolean avgPoln, olong numFreq)
{
  BaselineData *out=NULL;
  olong iIF, iPoln, iFreq, ip, naxis[1], npoln;

  out          = g_malloc0(sizeof(BaselineData));
  out->ant1    = ant1;
  out->ant2    = ant2;
  out->numIF   = numIF;
  out->numPoln = numPoln;
  out->numFreq = numFreq;
  if (avgPoln) npoln = 1;
  else  npoln = numPoln;
  out->IFindex    = g_malloc0(numIF*npoln*sizeof(olong));
  out->Stokeindex = g_malloc0(numIF*npoln*sizeof(olong));
  out->WtPolnIF   = g_malloc0(numIF*npoln*sizeof(ofloat));
  out->wtArray    = g_malloc0(numIF*npoln*sizeof(ObitFArray*));
  out->phArray    = g_malloc0(numIF*npoln*sizeof(ObitFArray*));
  out->visArray   = g_malloc0(numIF*npoln*sizeof(ObitCArray*));
  out->dFreq      = g_malloc0(numIF*sizeof(ofloat*));

  /* Populate arrays */
  ip = 0;
  naxis[0] = numFreq;
  for (iPoln=0; iPoln<npoln; iPoln++) {
    for (iIF=0; iIF<numIF; iIF++) {
      out->IFindex[ip]    = iIF;
      out->Stokeindex[ip] = iPoln;
      out->wtArray[ip]    = ObitFArrayCreate("wtData",  1, naxis);
      out->phArray[ip]    = ObitFArrayCreate("phData",  1, naxis);
      out->visArray[ip]   = ObitCArrayCreate("visData",  1, naxis);
      ip++;
    } /* end IF loop */
  } /* end Poln loop */


  for (iIF=0; iIF<numIF; iIF++) {
    out->dFreq[iIF] = g_malloc0(numFreq*sizeof(ofloat));
    for (iFreq=0; iFreq<numFreq; iFreq++) out->dFreq[iIF][iFreq] = 0.0;
  }

  return out;
} /* end MakeBLData */

/**
 * Clear (zero) Baseline data structure
 * \param in   The structure to modify
 * \param avgPoln If TRUE, average polarisations (one poln data array)
 * \return NULL pointer
 */
static void ClearBLData (BaselineData *in, olong avgPoln)
{
  olong iIF, iPoln,ip, npoln;
  ofloat cmplx0[2]={0.0,0.0};


  if (in==NULL) return; /* anything there? */
  
  /* Zero arrays */
  if (avgPoln) npoln = 1;
  else  npoln = in->numPoln;
  ip = 0;
  for (iPoln=0; iPoln<npoln; iPoln++) {
    for (iIF=0; iIF<in->numIF; iIF++) {
      ObitFArrayFill(in->wtArray[ip],  0.0);
      ObitFArrayFill(in->phArray[ip],  0.0);
      ObitCArrayFill(in->visArray[ip], cmplx0);
      in->WtPolnIF[ip] = 0.0;
      ip++;
    } /* end IF loop */
  } /* end Poln loop */

} /* end ClearBLData */

/**
 * Delete Baseline data structure
 * \param in   The structure to delete
 * \param avgPoln If TRUE, average polarisations (one poln data array)
 * \return NULL pointer
 */
static BaselineData* KillBLData (BaselineData *in, gboolean avgPoln)
{
  olong iIF, iPoln, npoln, ip;

  if (in==NULL) return NULL; /* anything there? */
  if (in->dFreq) {
    for (iIF=0; iIF<in->numIF; iIF++) {
      if (in->dFreq[iIF]) g_free(in->dFreq[iIF]);
    }
    g_free(in->dFreq); in->dFreq = NULL;
  }
  
  /* Delete wt arrays */
  ip = 0;
  if (avgPoln) npoln = 1;
  else  npoln = in->numPoln;
  if (in->wtArray) {
    for (iPoln=0; iPoln<npoln; iPoln++) {
      for (iIF=0; iIF<in->numIF; iIF++) {
	in->wtArray[ip]   = ObitFArrayUnref(in->wtArray[ip]);
	ip++;
      } /* end IF loop */
    } /* end Poln loop */
    g_free(in->wtArray);
  } /* end if wtArray */

  /* Delete phase arrays */
  ip = 0;
  if (in->phArray) {
    for (iPoln=0; iPoln<npoln; iPoln++) {
      for (iIF=0; iIF<in->numIF; iIF++) {
	in->phArray[ip]   = ObitFArrayUnref(in->phArray[ip]);
	ip++;
      } /* end IF loop */
    } /* end Poln loop */
    g_free(in->phArray);
  } /* end if phArray */

  /* Delete vis arrays */
  ip = 0;
  if (in->visArray) {
    for (iPoln=0; iPoln<npoln; iPoln++) {
      for (iIF=0; iIF<in->numIF; iIF++) {
	in->visArray[ip]   = ObitCArrayUnref(in->visArray[ip]);
	ip++;
      } /* end IF loop */
    } /* end Poln loop */
    g_free(in->visArray);
  } /* end if visArray */

  if (in->IFindex)    g_free(in->IFindex);
  if (in->Stokeindex) g_free(in->Stokeindex);
  if (in->WtPolnIF)   g_free(in->WtPolnIF);
  g_free(in);
  return NULL;
} /* end KillBLData */

/**
 * Make Scan data structure and all baseline data structures
 * \param inUV        Data being processed.
 * \param FFTOverSamp FFT oversampling ratio
 * \return pointer to the new structure.
 */
static ScanData* MakeScanData (ObitUV *inUV, gboolean avgPoln)
{
  ScanData *out=NULL;
  olong ant1, ant2, ip, iIF, iFreq, ifq, numIF, numFreq, numPoln, maxAnt, i, n;
  olong fincf, fincif, naxis[1];
  ofloat refChan;
  odouble refFreq;

  /* Get info from uv data */
  if (inUV->myDesc->jlocs>=0)
    numPoln = MIN (2, inUV->myDesc->inaxes[inUV->myDesc->jlocs]);
  else numPoln = 1;
  if (inUV->myDesc->jlocif>=0)
    numIF  = inUV->myDesc->inaxes[inUV->myDesc->jlocif];
  else numIF  = 1;
  if (inUV->myDesc->jlocf>=0)
    numFreq  = inUV->myDesc->inaxes[inUV->myDesc->jlocf];
  else numFreq  = 1;
  maxAnt = inUV->myDesc->numAnt[0];
  refFreq = inUV->myDesc->crval[inUV->myDesc->jlocf];  /* Reference fequency */
  refChan = inUV->myDesc->crpix[inUV->myDesc->jlocf];  /* Reference channel */
  /* Increments in the frequency table */
  fincf  = MAX (1, (inUV->myDesc->incf  / 3) / inUV->myDesc->inaxes[inUV->myDesc->jlocs]);
  fincif = MAX (1, (inUV->myDesc->incif / 3) / inUV->myDesc->inaxes[inUV->myDesc->jlocs]);

  /* Create output */
  out           = g_malloc0(sizeof(ScanData));
  out->maxAnt   = maxAnt;
  out->numBase  = maxAnt*(maxAnt-1)/2;
  out->suba     = -1;
  out->sid      = -1;
  out->fqid     = -1;
  out->timeAvg  = 0.0;
  out->freqAvg  = fabs(inUV->myDesc->cdelt[inUV->myDesc->jlocf]);
  out->refChan  = refChan;
  out->timec    = 0.0;
  out->timei    = 0.0;

  out->BLData   = g_malloc0(out->numBase*sizeof(BaselineData*));

  ip = 0;
  for (ant1=1; ant1<maxAnt; ant1++) {
    for (ant2=ant1+1; ant2<=maxAnt; ant2++) {
      out->BLData[ip] = MakeBLData (ant1, ant2, numIF, numPoln, avgPoln, numFreq);

      /* Set frequency offset information */
      for (iIF=0; iIF<numIF; iIF++) {
	/* Frequencies */
	for (iFreq=0; iFreq<numFreq; iFreq++) {
	  ifq = iIF*fincif + iFreq*fincf;  /* index in IF/freq table */
	  (out->BLData[ip])->dFreq[iIF][iFreq] = 1.0e-9*(inUV->myDesc->freqArr[ifq] - refFreq);
	}
      }
      ip++;
    } /* end loop over ant2 */
  } /* end loop over ant1 */

  /* Stacked antenna data - first frequencies */
  naxis[0] = out->BLData[0]->numFreq*out->BLData[0]->numIF;
  out->antFreqOff = ObitFArrayCreate("Stacked FreqO", 1, naxis);
  ip = 0;
  for (iIF=0; iIF<out->BLData[0]->numIF; iIF++) {
    for (iFreq=0; iFreq<out->BLData[0]->numFreq; iFreq++) {
      out->antFreqOff->array[ip++] = out->BLData[0]->dFreq[iIF][iFreq];
    }
  }

  /* Then Stacked phase arrays  */
  n = out->maxAnt*out->BLData[0]->numPoln;
  out->antStackPh = g_malloc0(n*sizeof(ObitFArray*));
  out->antStackWt = g_malloc0(n*sizeof(ObitFArray*));
  for (i=0; i<n; i++) {
    out->antStackPh[i] = ObitFArrayCreate("Stacked Phase",  1, naxis);
    out->antStackWt[i] = ObitFArrayCreate("Stacked Weight", 1, naxis);
  }
 
  return out;
} /* end MakeScanData */

/**
 * Clear (zero) Scan data structure and all baseline structs
 * \param in   The structure to delete
 * \return NULL pointer
 */
static void ClearScanData (ScanData *in)
{
  olong i, n;
  if (in==NULL) return; /* anything there? */

  in->suba     = -1;
  in->sid      = -1;
  in->fqid     = -1;

  if (in->BLData) {
    for (i=0; i<in->numBase; i++) {
      ClearBLData (in->BLData[i], in->avgPoln);
    }
  }

  /* Stacked phases */
  n = in->maxAnt*in->BLData[0]->numPoln;
  for (i=0; i<n; i++) {
    ObitFArrayFill(in->antStackPh[i], 0.0);
    ObitFArrayFill(in->antStackWt[i], 0.0);
  }
 
} /* end ClearScanData */

/**
 * Delete Scan data structure and all baseline structs
 * \param in   The structure to delete
 * \return NULL pointer
 */
static ScanData* KillScanData (ScanData *in)
{
  olong i, n;
  if (in==NULL) return NULL; /* anything there? */

  /* Stacked phases */
  in->antFreqOff = ObitFArrayUnref(in->antFreqOff);
  n = in->maxAnt*in->BLData[0]->numPoln;
  for (i=0; i<n; i++) {
    in->antStackPh[i] = ObitFArrayUnref(in->antStackPh[i]);
    in->antStackWt[i] = ObitFArrayUnref(in->antStackWt[i]);
  }
  g_free(in->antStackPh);
  g_free(in->antStackWt);

  if (in->BLData) {
    for (i=0; i<in->numBase; i++) {
      in->BLData[i] = KillBLData (in->BLData[i], in->avgPoln);
    }
    g_free(in->BLData);
  }
  g_free(in);
  return NULL;
} /* end KillScanData */


#if HAVE_GSL==1  /* GSL stuff */
/**
 * Fringe fitting model calculating routine for gsl least squares fitter
 * \param coef   Coefficient array per antenna (phase, delay, dispersion)
 * \param params Data structure (ObitUVGSolveWB)
 * \param f      [out] function residuals
 * \returns GSL completion code
 */
static int fringeFitFunc (const gsl_vector *coef, void *params, gsl_vector *f)
{
  ObitUVGSolveWB *in  = (ObitUVGSolveWB*)params;
  ScanData     *data  = in->scanData; 
  olong i, ndata, kndx, lo, hi, cnt;
  odouble sum, freq0=0.0;
  ofloat lcoef[10], *vdata, *vwt, *vfreq, resid, phase, twopi=2.0*G_PI;

  /* fill local array of coeffieients */
  for (i=0; i<in->ncoefCoarse; i++) 
    lcoef[i] = (ofloat)gsl_vector_get(coef, i);

  /* DEBUG
  if(data->curAnt==22) {
    fprintf (stdout, "Model IF %d %f %f \n", data->curIF, lcoef[0], lcoef[1]);
  } */

  /* Loop over data */
  kndx = (data->curAnt-1) + data->curPoln*data->maxAnt;
  vdata = data->antStackPh[kndx]->array;
  vwt   = data->antStackWt[kndx]->array;
  vfreq = data->antFreqOff->array;
  ndata = data->antStackPh[kndx]->naxis[0];
  sum   = 0.0;   /* RMS residual accumulator */
  cnt    = 0;
  /* Set range of data values */
  if (data->avgIF) {  /* all at once */
    lo = 0;
    hi = ndata;
  } else {   /* Only one IF */
    lo = data->curIF*data->numChan;
    hi = lo + data->numChan;
  }
  if (!data->avgIF) freq0 = vfreq[lo];   /* First frequency of IF */
  for (i=lo; i<hi; i++) {
    phase = lcoef[0] + twopi*lcoef[1]*(vfreq[i]-freq0);
    resid = phase- vdata[i];
    /* Take out any whole turns */
    resid = fmodf(resid, twopi);
    /* Nearest half turn */
    if (resid>G_PI) resid -= twopi;
    else if (resid<-G_PI) resid += twopi;

    /* DEBUG 
    if(data->curAnt==22) {
      fprintf (stdout, " %d, %f %f %f  %f %f \n", 
	       i, phase*57.296,  vdata[i]*57.296, resid*57.296, vwt[i], vfreq[i]-freq0);
    }*/

    /* RMS resid */
    sum += resid*resid;

    /* Weight */
    resid *= vwt[i];
    if (vwt[i]>0.0) cnt++;
    
    /* Save residual */
    gsl_vector_set(f, i-lo, resid);	
  } /* end loop over data */

  /* Save RMS residual */
  data->RMSRes = (ofloat)sqrt(sum)/(hi-lo);
  /* DEBUG
  if(data->curAnt==22) {
    fprintf (stdout, " RMS resid  %f \n", data->RMSRes*57.296);
  } */

  return GSL_SUCCESS;
} /* end  fringeFitFunc */

/**
 * Fringe fitting function calculating Jacobean for gsl least squares fitter
 * This is the partial derivative of the residuals matrix
 * \param coef   Coefficient array (phase, delay, per non reference antenna)
 * \param params Data structure  (ObitUVGSolveWB)
 * \param J      [out] Jacobean values
 * \returns GSL completion code
 */
static int fringeFitJacob (const gsl_vector *coef, void *params, gsl_matrix *J)
{
  ObitUVGSolveWB *in = (ObitUVGSolveWB*)params;
  ScanData     *data  = in->scanData; 
  olong i, ndata, kndx, lo, hi;
  odouble freq0=0.0;
  ofloat lcoef[10], *vdata, *vwt, *vfreq, resid, phase, part1, part2, twopi=2.0*G_PI;
 
  /* Fill local array of coeffieients */
  for (i=0; i<in->ncoefCoarse; i++) 
    lcoef[i] = (ofloat)gsl_vector_get(coef, i);

  /* Loop over data */
  kndx = (data->curAnt-1) + data->curPoln*data->maxAnt;
  vdata = data->antStackPh[kndx]->array;
  vwt   = data->antStackWt[kndx]->array;
  vfreq = data->antFreqOff->array;
  ndata = data->antStackPh[kndx]->naxis[0];
  /* Set range of data values */
  if (data->avgIF) {  /* all at once */
    lo = 0;
    hi = ndata;
  } else {   /* Only one IF */
    lo = data->curIF*data->numChan;
    hi = lo + data->numChan;
  }
  if (!data->avgIF) freq0 = vfreq[lo];   /* First frequency of IF */
  for (i=lo; i<hi; i++) {
    phase = lcoef[0] + twopi*lcoef[1]*(vfreq[i]-freq0);
    resid = phase - vdata[i];
    /* Take out any whole turns */
    resid = fmodf(resid, twopi);
    /* Nearest half turn */
    if (resid>G_PI) resid -= twopi;
    else if (resid<-G_PI) resid += twopi;
    
    /* Weight */
    resid *= vwt[i];
    
    /* compute partials */
    part1 = vwt[i];
    part2 = vwt[i]*twopi*(vfreq[i]-freq0);
    /* DEBUG 
    part1 = resid;
    part2 = resid*twopi*(vfreq[i]-freq0);*/
    
    /* Update Jacobean */
    gsl_matrix_set (J, i-lo, 0,  part1);
    gsl_matrix_set (J, i-lo, 1,  part2);
  } /* end loop over data */

  return GSL_SUCCESS;

} /* end  fringeFitJacob */


/**
 * Compute both function and derivatives for  gsl least squares fitter
 * \param coef   Coefficient array per antenna (phase, delay, dispersion)
 * \param params Data structure
 * \param f      [out] function residuals
 * \param J      [out] Jacobean
 * \returns GSL completion code
 */
static int fringeFitFuncJacob (const gsl_vector *coef, void *params, gsl_vector *f, 
			       gsl_matrix *J)
{
  fringeFitFunc (coef, params, f);
  fringeFitJacob(coef, params, J);
  return GSL_SUCCESS;
} /* end  fringeFitFuncJacob */

/**
 * Fringe fitting model calculating routine for coarse IF least squares fitter
 * \param coef   Coefficient array per antenna (phase, delay, dispersion)
 * \param params Data structure (ObitUVGSolveWB)
 * \param f      [out] function residuals
 * \returns GSL completion code
 */
static int coarseFitFunc (const gsl_vector *coef, void *params, gsl_vector *f)
{
  ObitUVGSolveWB *in = (ObitUVGSolveWB*)params;
  olong i, ndata;
  ofloat lcoef[10], *vdata, *vwt, *vfreq, resid, phase, twopi=2.0*G_PI;

  /* Create/fill local array of coeffieients */
  for (i=0; i<in->ncoefCoarse; i++) 
    lcoef[i] = (ofloat)gsl_vector_get(coef, i);

  /* Loop over data */
  vdata = in->coarseData;
  vwt   = in->coarseWt;
  vfreq = in->coarseFreq;
  ndata = in->ndataCoarse;
  for (i=0; i<ndata; i++) {
    phase = lcoef[0] + twopi*lcoef[1]*vfreq[i];
    resid = phase- vdata[i];
    /* Take out any whole turns */
    resid = fmodf(resid, twopi);
    /* Nearest half turn */
    if (resid>G_PI) resid -= twopi;
    else if (resid<-G_PI) resid += twopi;
    
    /* Weight */
    resid *= vwt[i];
    
    /* Save residual */
    gsl_vector_set(f, i, resid);	
  } /* end loop over data */

  return GSL_SUCCESS;
} /* end  coarseFitFunc */

/**
 * Fringe fitting function calculating Jacobean for coarse IF least squares fitter
 * This is the partial derivative of the residuals matrix
 * \param coef   Coefficient array (phase, delay)
 * \param params Data structure  (ObitUVGSolveWB)
 * \param J      [out] Jacobean values
 * \returns GSL completion code
 */
static int coarseFitJacob (const gsl_vector *coef, void *params, gsl_matrix *J)
{
  ObitUVGSolveWB *in = (ObitUVGSolveWB*)params;
  olong i, ndata;
  ofloat lcoef[10], *vdata, *vwt, *vfreq, resid, phase, part1, part2, twopi=2.0*G_PI;
 
  /* Fill local array of coeffieients */
  for (i=0; i<in->ncoefCoarse; i++) 
    lcoef[i] = (ofloat)gsl_vector_get(coef, i);

  /* Loop over data */
  vdata = in->coarseData;
  vwt   = in->coarseWt;
  vfreq = in->coarseFreq;
  ndata = in->ndataCoarse;
  for (i=0; i<ndata; i++) {
    phase = lcoef[0] + twopi*lcoef[1]*vfreq[i];
    resid = phase - vdata[i];
    /* Take out any whole turns */
    resid = fmodf(resid, twopi);
    /* Nearest half turn */
    if (resid>G_PI) resid -= twopi;
    else if (resid<-G_PI) resid += twopi;
    
    /* Weight */
    resid *= vwt[i];
    
    /* compute partials*/
    part1 = vwt[i];
    part2 = vwt[i]*twopi*vfreq[i]; 
    /* DEBUG
    part1 = resid;
    part2 = resid*twopi*(vfreq[i]); */
    
    /* Update Jacobean */
    gsl_matrix_set (J, i, 0, part1);
    gsl_matrix_set (J, i, 1, part2);
  } /* end loop over data */

  
  return GSL_SUCCESS;
} /* end  coarseFitJacob */


/**
 * Compute both function and derivatives for coarse IF least squares fitter
 * \param coef   Coefficient array per antenna (phase, delay, dispersion)
 * \param params Data structure
 * \param f      [out] function residuals
 * \param J      [out] Jacobean
 * \returns GSL completion code
 */
static int coarseFitFuncJacob (const gsl_vector *coef, void *params, gsl_vector *f, 
			       gsl_matrix *J)
{
  coarseFitFunc (coef, params, f);
  coarseFitJacob(coef, params, J);
  return GSL_SUCCESS;
} /* end  coarseFitFuncJacob */
#endif /* GSL stuff */

/**
 *  Add corresponding elements of two arrays, second adjusted to beginning 
 *  of first
 *  out = in1 + in2,  if either is blanked the result is blanked
 * Only works for 1D
 * \param in1  Input object with data
 * \param in2  Input object with data
 * \param out  Output array (may be an input array).
 */
void ObitFArrayAddEnd (ObitFArray* in1, ObitFArray* in2, ObitFArray* out)
{
  olong i, j;
  ofloat fblank = ObitMagicF();

  if (in1->arraySize==in2->arraySize) {   /* Inputs same size */
    j = 0;
    
    for (i=0; i<in2->arraySize; i++) {
      if ((in1->array[j]!=fblank) && (in2->array[i]!=fblank)) 
	out->array[j] = in1->array[i] + in2->array[j];
      else out->array[j] = fblank;
      j++;
    }
  } else {  /* Inputs different size */
    j = 0; 
    
    for (i=0; i<in2->arraySize; i++) {
      if ((in1->array[j]!=fblank) && (in2->array[i]!=fblank)) 
	out->array[j] = in1->array[j] + in2->array[i];
      else out->array[j] = fblank;
      j++;
    }
  }
} /* end ObitFArrayAddEnd */

/**
 *  Add harmonic sum of corresponding elements of two arrays, 
 * second adjusted to beginning of first
 *  out = 1/(1/in1 + 1/in2),  if either is blanked the result is blanked
 * Only works for 1D
 * \param in1  Input object with data
 * \param in2  Input object with data
 * \param out  Output array (may be an input array).
 */
void ObitFArrayHarmAddEnd (ObitFArray* in1, ObitFArray* in2, ObitFArray* out)
{
  olong i, j;
  ofloat harm;
  ofloat fblank = ObitMagicF();
  
  if (in1->arraySize==in2->arraySize) {   /* Inputs same size */
    j = 0;
    
    for (i=0; i<in2->arraySize; i++) {
      if ((in1->array[j]!=fblank) && (in2->array[i]!=fblank) && 
	  (in1->array[j]!=0.0) && (in2->array[i]!=0.0)) {
	harm = 1.0 / ((1.0/in1->array[j]) + (1.0/in2->array[i]));
	out->array[j] = harm;
      }
      else out->array[j] = fblank;
      j++;
    }
  } else {  /* Inputs different size */
    j = 0; 
    
    for (i=0; i<in2->arraySize; i++) {
      if ((in1->array[j]!=fblank) && (in2->array[i]!=fblank) && 
	  (in1->array[j]!=0.0) && (in2->array[i]!=0.0)) {
	harm = 1.0 / ((1.0/in1->array[j]) + (1.0/in2->array[i]));
	out->array[j] = harm;
      }
      else out->array[j] = fblank;
      j++;
    }
  }
} /* end ObitFArrayHarmAddEnd */

/**
 *  Add square of sum of corresponding elements of two arrays, 
 *  second adjusted to beginning of first
 *  out = in1 + in2,  if either is blanked the result is blanked
 * Only works for 1D
 * \param in1  Input object with data
 * \param in2  Input object with data
 * \param out  Output array (may be an input array).
 * \param count incremented by number of valid data.
 */
void ObitFArrayAddEnd2 (ObitFArray* in1, ObitFArray* in2, ObitFArray* out, olong *count)
{
  olong i, j;
  ofloat fblank = ObitMagicF();
  
  if (in1->arraySize==in2->arraySize) {   /* Inputs same size */
    j = 0;
    
    for (i=0; i<in2->arraySize; i++) {
      if ((in1->array[j]!=fblank) && (in2->array[i]!=fblank)) {
	out->array[j] = (in1->array[i] + in2->array[j])*(in1->array[i] + in2->array[j]);
	(*count)++;
      } else out->array[j] = fblank;    
      j++;
    }
  } else {  /* Inputs different size */
    j = 0; 
    
    for (i=0; i<in2->arraySize; i++) {
      if ((in1->array[j]!=fblank) && (in2->array[i]!=fblank)) {
	out->array[j] = (in1->array[j] + in2->array[i]) * (in1->array[j] + in2->array[i]);
	(*count) ++;
      }  else out->array[j] = fblank;
      j++;
    }
  }
} /* end ObitFArrayAddEnd2 */

/**
 *  Aaccumulate square of harmonic sum of corresponding elements of two arrays, 
 *  second adjusted to beginning of first
 *  out = 1/(1/in1 + 1/in2),  if either is blanked the result is blanked
 * Only works for 1D
 * \param in1  Input object with data
 * \param in2  Input object with data
 * \param out  Output array (may be an input array).
 * \param count incremented by number of valid data.
 */
void ObitFArrayHarmAccEnd2 (ObitFArray* in1, ObitFArray* in2, ObitFArray* out,
			   olong *count)
{
  olong i, j;
  ofloat harm;
  ofloat fblank = ObitMagicF();

  if (in1->arraySize==in2->arraySize) {   /* Inputs same size */
    j = 0;
   
    for (i=0; i<in2->arraySize; i++) {
      if ((in1->array[j]!=fblank) && (in2->array[i]!=fblank) && 
	  (in1->array[j]!=0.0) && (in2->array[i]!=0.0) &&
	  (out->array[j]!=fblank)) {
	harm = 1.0 / ((1.0/in1->array[j]) + (1.0/in2->array[i]));
	out->array[j] += harm*harm;
	(*count)++;
      }
      j++;
    }
  } else {  /* Inputs different size */
    j = 0; 
    
    for (i=0; i<in2->arraySize; i++) {
      if ((in1->array[j]!=fblank) && (in2->array[i]!=fblank) && 
	  (in1->array[j]!=0.0) && (in2->array[i]!=0.0) &&
	  (out->array[j]!=fblank)) {
	harm = 1.0 / ((1.0/in1->array[j]) + (1.0/in2->array[i]));
	out->array[j] *= harm*harm;
	(*count)++;
      }
      j++;
    }
  }
} /* end ObitFArrayHarmAddEnd2 */

/**
 *  Multiply the elements of a CArray by the elements of an FArray, 
 *  Fin adjusted to the beginning of Cin
 *    out = Cin * Fin
 * Only works for 1D
 * \param Cin  Input CArray
 * \param Fin  Input FArray
 * \param out  Output CArray
 */
void ObitCArrayFMulEnd (ObitCArray* Cin, ObitFArray* Fin, ObitCArray* out)
{
  olong i, j;
  ofloat trc, tic, tr;

  /* error checks */
  g_assert (ObitCArrayIsA(Cin));
  g_assert (ObitFArrayIsA(Fin));
  g_assert (ObitCArrayIsCompatable(Cin, out));

  /* Multiply */
  j = 0; 
  for (i=0; i<Fin->arraySize; i++) {
    trc = Cin->array[j];
    tic = Cin->array[j+1];
    tr  = Fin->array[i];
    out->array[j]   = trc * tr;
    out->array[j+1] = tic * tr;
    j += 2;
  }
}  /* end ObitCArrayFMulEnd */


/**
 * Determines location and value of the complex function in in->cWork2
 * Shuffles the center to the center of in->cWork2, 
 * determines the location of the peak and interpolates the value.
 * Uses in->fWork2 as a scratch array.
 * \param in      Input gain fitter object. 
 * \param ppos    [out] Position of peak (0-rel cells)
 * \param pval    [out] Complex value of delay function at peak
 */
static void  
FindFFTPeak (ObitUVGSolveWB *in, ofloat *ppos, ofloat pval[2])
{
  olong i, n, pos[0], i1, i2;
  ofloat *inData, *outData, pmax, sum, wt, pixel;

  /* Shuffle data to interpolator array */
  n = in->cWork2->arraySize/2;
  inData  = in->cWork2->array;
  outData = in->FFTFitArray->array;
  for (i=0; i<2*n; i++) outData[2*n+i] = inData[i];
  for (i=0; i<2*n; i++) outData[i]     = inData[2*n+i];

  /* Get amplitudes in in->fwork2 */
  ObitCArrayAmp (in->FFTFitArray, in->fWork2);

  /* Find peak */
  pmax = ObitFArrayMax (in->fWork2, pos);

  /* Centroid (up to) closest 5 points */
  i1 = MAX (0, pos[0]-2);
  i2 = MIN (in->fWork2->arraySize-1,  pos[0]+2);
  sum = wt = 0.0; 
  for (i=i1; i<=i2; i++) {
    wt  += in->fWork2->array[i];
    sum += i*in->fWork2->array[i];
  }
  if (fabs(wt)>0.0) {
    pixel = sum / wt;
    *ppos = pixel;
  } else { /* No data */
    *ppos = 0.0;
  }

  /* Interpolate */
  ObitCInterpolate1D (in->myCInterp, pixel, pval);
} /* end FindFFTPeak  */

/**
 *  Zero elements of array in1 where array in2 is blanked
 *  out = in1 or zero where in2 is blank
 * Only operates for the number of elements in in1 (WATCH out - special case)
 * \param in1  Input object with data
 * \param in2  Input object with blanking
 * \param out  Output array (may be an input array).
 */
void ObitFArrayZeroBlank (ObitFArray* in1, ObitFArray* in2, ObitFArray* out)
{
  olong i;
  ofloat fblank = ObitMagicF();

   /* error checks */

  for (i=0; i<in1->arraySize; i++) {
    if (in2->array[i]!=fblank)
      out->array[i] = in1->array[i];
    else out->array[i] = 0.0;
  }
} /* end ObitFArrayZeroBlank */

/**
 *  Fill the beginning of a CArray with (cos,sin) of phases in in.
 *  out = complex(cos(in), sin(in))
 * Allows magic value blanking
 * Only works for 1D
 * \param in   Input FArray
 * \param out  Output CArray
 */
void ObitCArrayFSinCos (ObitFArray* in, ObitCArray* out)
{
  olong i, j, k, n;
  ofloat stemp[16], ctemp[16], fblank = ObitMagicF();

  /* j = 2*(out->arraySize - in->arraySize);  Where does in start in out? */
  j = 0; 
  for (i=0; i<in->arraySize; i+=16) {
    /* use fast sin/cos in blocks of 16 */
    n = MIN(16, in->arraySize-i);
    ObitSinCosVec (n, &in->array[i], stemp, ctemp);
    for (k=0; k<n; k++){
      if (in->array[i+k]!=fblank) {
	out->array[j]   = ctemp[k];
	out->array[j+1] = stemp[k];
      } else {
	out->array[j]   = 0.0;
	out->array[j+1] = 0.0;
      }
      j += 2;
    }
  }
} /* end ObitCArrayFSinCos */

/**
 *  Approximate SNR from the peak amplitude of the FFT 
 *  of a set of phasors. The phasors are intended to be 
 *  stacked baselines with combinations arrprximating that 
 *  of a given antenna to a given reference antenna.
 *  This approximation breaks down for large values of SNR (>50).
 *  From the AIPSish FRING.FOR/FRNSR2
 * \param  amp   Amplitude of coherent sum
 * \param  amp   Normalized amplitude of coherent sum
 * \param  sumw  Sum of weights used in stacking
 * \param  sumww Sum of weights squared
 * \param  tfact Weight modification array; used if unequal 
 *               integration times in the data.  Use 1.0
 * \param  count Count of data points added to array FFTed.
 * \return Approximation to SNR
 */
static ofloat 
decorSNR (ofloat amp, ofloat sumw, ofloat sumww, ofloat tfact, olong count)
{
  ofloat snr = 0.0;

  amp = MIN (amp, 0.999);   /* Disallow correlations > 1.0 */
  if ((amp> 0.0) && (sumw>0.0)) {
    snr = tan(1.570796*amp);
    snr = pow(snr,1.163);
    snr *= sqrt(sumw/sqrt(sumww/count));
  }
  /*  Adjust for unequal integration times in the data 
  tfact = MAX(1.0, tfact);
  snr /= sqrt (tfact);*/
  return snr;
} /*  end decorSNR */
