/* $Id:  $  */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2009                                               */
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

#include "ObitUVRFIXize.h"
#include "ObitTableUtil.h"
#include "ObitTableSNUtil.h"
#include "ObitUVUtil.h"
#include "ObitUVSoln.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitUVRFIXize.c
 * ObitUVRFIXize class function definitions.
 * This clas supports RFI estimation and removel from UV data
 * This class is derived from the Obit base class.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitUVRFIXize";

/** Function to obtain parent ClassInfo */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo structure ObitUVRFIXizeClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitUVRFIXizeClassInfo myClassInfo = {FALSE};

/*--------------- File Global Variables  ----------------*/


/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitUVRFIXizeInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitUVRFIXizeClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitUVRFIXizeClassInfoDefFn (gpointer inClass);

/** Private: Get fringe rates from UVSoln object. */
static void getAntFR (ObitUVRFIXize *in, ObitUVSoln *SNSoln, ofloat time, 
		      ofloat *frRate, ObitErr* err);

/** Private: Get delays and fringe rates from UVSoln object. */
static void getAntDelFR (ObitUVRFIXize *in, ObitUVSoln *SNSoln, ofloat time, 
			 ofloat *Real, ofloat *Imag, ofloat *delay, ofloat *frRate, 
			 ObitErr* err);

/** Private:  Read calibration for a new time into the internal arrays. */
static void ObitUVRFIXizeNewTime (ObitUVRFIXize *in, ofloat time,
				  ObitErr *err);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitUVRFIXize* newObitUVRFIXize (gchar* name)
{
  ObitUVRFIXize* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitUVRFIXizeClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitUVRFIXize));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitUVRFIXizeInit((gpointer)out);

 return out;
} /* end newObitUVRFIXize */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitUVRFIXizeGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitUVRFIXizeClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitUVRFIXizeGetClass */

/**
 * Creates an ObitUVRFIXize 
 * \param name  An optional name for the object.
 * \param inUV  UV data with input SN table
 * \return the new object.
 */
ObitUVRFIXize* ObitUVRFIXizeCreate (gchar* name, ObitUV *inUV, 
				    ObitUV *residUV, ObitUV *outUV)
{
  ObitUVRFIXize* out;

  /* Create basic structure */
  out = newObitUVRFIXize (name);

  /* save input uv data pointers */
  out->myUV    = ObitUVRef(inUV);
  out->residUV = ObitUVRef(residUV);
  out->outUV   = ObitUVRef(outUV);

  /* Info */
  out->numAnt = inUV->myDesc->numAnt[0];
  out->numBL  = out->numAnt*out->numAnt/2;
  if (inUV->myDesc->jlocif>=0)
    out->numIF = inUV->myDesc->inaxes[inUV->myDesc->jlocif];
  else
    out->numIF = 1;

  return out;
} /* end ObitUVRFIXizeCreate */

/**
 * Counterrotate and average residual data 
 * \param in     UVRFIXize object
 * Control parameter on info element:
 * \li "solInt"    OBIT_float (1,1,1) Counter rotated SN table interval [def 1 min]
 * \li "timeInt"   OBIT_float (1,1,1) Data integration time in sec [def 10 sec].
 * \li "timeAvg"   OBIT_float  (1,1,1) Time interval over which to average 
 *                 (min) [def = 1 min.]
 *                 NB: this should be at least 2 integrations.
 * \li "doInvert"  OBIT_bool (1,1,1) If TRUE invert solution [def FALSE];
 * \param err    Error/message stack, returns if error.
 */
void ObitUVRFIXizeCounterRot (ObitUVRFIXize *in, ObitErr* err) 
{
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  olong ver, itemp;
  ofloat solInt=60, timeInt=10.0, timeAvg=0.95;
  gboolean Tr=TRUE, invert=FALSE;
  gchar *routine = "ObitUVRFIXizeCounterRot";
  
  /* Error checks */
  if (err->error) return ;  /* previous error? */
  g_assert (ObitUVRFIXizeIsA(in));

  /* Check compatability - allow no output */
  Obit_return_if_fail (((in->myUV->myDesc->lrec==in->outUV->myDesc->lrec) || (in->outUV==NULL)), err, 
		       "%s Input and output data records incompatable", routine);  

  /* Get/set solution table control */
  ObitInfoListGetTest(in->info, "solInt",  &type, dim, &solInt);
  ObitInfoListGetTest(in->info, "timeInt", &type, dim, &timeInt);
  ObitInfoListGetTest(in->info, "doInvert", &type, dim, &invert);
  dim[0] = dim[1] = dim[2] = dim[3] = dim[4] = 1;
  ObitInfoListAlwaysPut(in->residUV->info, "solInt",  OBIT_float, dim, &solInt);
  ObitInfoListAlwaysPut(in->residUV->info, "timeInt", OBIT_float, dim, &timeInt);
  ObitInfoListAlwaysPut(in->residUV->info, "doInvert", OBIT_bool, dim, &invert);

  /* Create zero fringe rate solution table */
  ver    = 0;
  in->SNTable = ObitTableSNGetZeroFR (in->residUV, in->residUV, ver, err);
  if (err->error) Obit_traceback_msg (err, routine, in->residUV->name);

  /* Get averaging time */
  ObitInfoListGetTest(in->info, "timeAvg", &type, dim, &timeAvg);

  /* Want calibration, averaging */
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut(in->residUV->info, "doCalSelect", OBIT_bool, dim, &Tr);
  itemp = 2;
  ObitInfoListAlwaysPut(in->residUV->info, "doCalib", OBIT_long, dim, &itemp);
  itemp = in->SNTable->tabVer;
  ObitInfoListAlwaysPut(in->residUV->info, "gainUse", OBIT_long, dim, &itemp);
  ObitInfoListAlwaysPut(in->residUV->info, "timeAvg", OBIT_float, dim, &timeAvg);

  /* Apply, writing averaged counter rotated data */
  in->RFIUV = ObitUVUtilAvgT (in->residUV, TRUE, in->residUV, err);
  if (err->error) Obit_traceback_msg (err, routine, in->residUV->name);

  /* Save averaging time in days */
  in->AvgTime = timeAvg/1440.0;  

} /* end of routine ObitUVRFIXizeCounterRot */ 

/**
 * Filter Counterrotated/averaged residual data
 * Reads in->filterUV and in->SNTable, decided which data hase RFI and rewrites
 * RFIUV either zeroing, flagging or passing the data depending on whether 
 * \li Zero: The data appears not to contain RFI
 * \li Flag: It is not possibbe to separate sky signal and RFI
 * \li Pass: the data contains RFI.
 * \param in     UVRFIXize object
 * Control parameter on info element:
 * \li "minRFI"   OBIT_float (1,1,1) Minimum RFI amplitude (Jy) to remove [def 50]
 * \li "timeAvg"  OBIT_float  (1,1,1) Time interval over whic RFIUV averaged 
 *                (min) [def = 1 min.]
 * \li "timeInt"  OBIT_float (1,1,1) Data integration time in sec [def 10 sec].
 * \li "maxRot"   OBIT_float (1,1,1) Max. fringe rotation (turns) in an integration 
 *                  to estimate RFI[def 2.0].
 * \li "minRot"   OBIT_float (1,1,1) Min. fringe rotation (turns) in an integration,
 *                  data with smaller values flagged if > minRFI [def 2.0].
 * \param err    Error/message stack, returns if error.
 */
void ObitUVRFIXizeFilter (ObitUVRFIXize *in, ObitErr* err) 
{
  ObitIOCode retCode;
  gboolean done, zero, blank;
  ofloat *frRate=NULL, fr, fblank = ObitMagicF();
  ofloat timeInt=10.0, timeAvg=0.95, minamp2=50.0, amp2;
  ofloat maxRot=2.0, minRot = 0.25;
  odouble freq;
  olong i, jndx, indx, ant1, ant2;
  olong incs, incif, incf,  nif, nfreq, nstok, iif, ichan, istok;
  ObitUVDesc *inDesc;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM]={1,1,1,1,1};
  gchar *interMode="AMBG";
  gchar *routine = "ObitUVRFIXizeFilter";
  
  /* Error checks */
  if (err->error) return ;  /* previous error? */
  g_assert (ObitUVRFIXizeIsA(in));

  /* Make sure in->RFIUV available */
  Obit_return_if_fail ((in->RFIUV), err, 
		       "%s Averaged counterrotated data not available (run Counterrot)", 
		       routine);  

  /* Get averaging time */
  ObitInfoListGetTest(in->info, "timeAvg", &type, dim, &timeAvg);
  timeAvg *= 60.0;  /* to seconds */
  /* Data integration time */
  ObitInfoListGetTest(in->info, "timeInt", &type, dim, &timeInt);
  /* Minimum amplitude to subtract */
  ObitInfoListGetTest(in->info, "minRFI",  &type, dim, &minamp2);
  /* Want square for test */
  minamp2 = minamp2*minamp2;
  /* Minimum rotation */
  ObitInfoListGetTest(in->info, "minRot",  &type, dim, &minRot);
  /* Maximum rotation */
  ObitInfoListGetTest(in->info, "maxRot",  &type, dim, &maxRot);

  /* Local pointers */
  inDesc  = in->RFIUV->myDesc;

  /* How many IFs, freq, poln */
  nfreq = inDesc->inaxes[inDesc->jlocf];
  if (inDesc->jlocf>=0) nif = inDesc->inaxes[inDesc->jlocif];
  else nif = 1;
  if (inDesc->jlocs>=0) nstok = inDesc->inaxes[inDesc->jlocs];
  else nstok = 1;

  /* Data always expanded to 3 words per vis */
  incs  = 3 * inDesc->incs  / inDesc->inaxes[0];
  incf  = 3 * inDesc->incf  / inDesc->inaxes[0];
  incif = 3 * inDesc->incif / inDesc->inaxes[0];

  /* Create uvSoln object */
  in->SNSoln = ObitUVSolnCreate(in->name, in->residUV);
  dim[0] = strlen(interMode);
  ObitInfoListAlwaysPut(in->SNSoln->info,"interMode", OBIT_string, dim, interMode);
  ObitUVSolnStartUp (in->SNSoln, err);
  if (err->error) Obit_traceback_msg (err, routine, in->RFIUV->name);

  /* Fringe rate array */
  frRate = g_malloc0(in->numAnt*sizeof(float));

 /* Open Input Data */
  retCode = ObitUVOpen (in->RFIUV, OBIT_IO_ReadWrite, err);
  if ((retCode != OBIT_IO_OK) || (err->error>0)) 
    Obit_traceback_msg (err, routine, in->RFIUV->name);

  /* Loop over data */
  done = (retCode != OBIT_IO_OK);
  while (!done) {
    
    /* read buffer */
    retCode = ObitUVRead (in->RFIUV, NULL, err);
    if (err->error) Obit_traceback_msg (err, routine, in->RFIUV->name);

    done = (retCode == OBIT_IO_EOF); /* done? */
    if (done) break;

    /* Modify data */
    for (i=0; i<inDesc->numVisBuff; i++) { /* loop over visibilities */
      jndx = i*inDesc->lrec;
      /* Update zero fringe rate table if necessary */
      if (in->RFIUV->buffer[jndx+inDesc->iloct]> in->SNSoln->CalTime){
	getAntFR (in, in->SNSoln, in->RFIUV->buffer[jndx+inDesc->iloct], frRate, err);
	if (err->error) Obit_traceback_msg (err, routine, in->SNSoln->name);
      }
      /* Which baseline? */
      ant1 = (in->RFIUV->buffer[jndx+inDesc->ilocb] / 256.0) + 0.001;
      ant2 = (in->RFIUV->buffer[jndx+inDesc->ilocb] - ant1 * 256) + 0.001;

      /* loop over IF */
      for (iif=0; iif<nif; iif++) {
	/* Loop over frequency channel */
	for (ichan=0; ichan<nfreq; ichan++) { 
	  /* Baseline fringe rate in Hz for this channel */
	  if ((frRate[ant1-1]!=fblank) && (frRate[ant2-1]!=fblank)) {
	    freq = inDesc->freqIF[iif] + 
	      (ichan+1.0 - inDesc->crpix[inDesc->jlocf]) * inDesc->cdelt[inDesc->jlocf];
	    fr   = (frRate[ant1-1] - frRate[ant2-1]) * freq;
	  } else fr = fblank;
	  /* Loop over Stokes */
	  for (istok=0; istok<nstok; istok++) {
	    indx = jndx + inDesc->nrparm + incs*istok + incf*ichan + incif*iif;
	    blank = FALSE;
	    zero  = FALSE;
	    if ((in->RFIUV->buffer[indx+2]>0.0) && (fr!=fblank)) {
	      /* What to do? */
	      amp2 = in->RFIUV->buffer[indx+0]*in->RFIUV->buffer[indx+0] +
		in->RFIUV->buffer[indx+1]*in->RFIUV->buffer[indx+1];
	      zero = amp2 < minamp2;  /* Amp below RFI threshold */
	      /* If more than maxRot turns in an integration ignore */
	      zero = zero || (fabs(timeInt*fr)>maxRot);
	      /* Need at least minRot turn of phase in data averaging */
	      if (fabs(timeAvg*fr)<minRot) {  
		blank = TRUE;
		/*zero = FALSE;*/
	      }
	      if (zero) {
		in->RFIUV->buffer[indx+0] = 0.0;
		in->RFIUV->buffer[indx+1] = 0.0;
	      } else if (blank) {
		in->RFIUV->buffer[indx+2] = 0.0;
	      }
	    } /* end if correlation valid */
	  } /* end loop over polarization */
	} /* end loop over frequency */
      } /* end loop over IF */
    } /* End loop over buffer */
    
    /* Rewrite buffer */
    retCode = ObitUVWrite (in->RFIUV, NULL, err);
    if (err->error) Obit_traceback_msg (err, routine, in->RFIUV->name);
    /* Reset file offset to keep from causing problem in next red */
    in->RFIUV->myDesc->firstVis -= in->RFIUV->myDesc->numVisBuff;
    ((ObitUVDesc*)in->RFIUV->myIO->myDesc)->firstVis = /* oh bother */
      in->RFIUV->myDesc->firstVis;
  } /* End loop over data */
  
  /* Close up */
  retCode = ObitUVClose(in->RFIUV, err);
  
  /* Cleanup */
  ObitUVSolnShutDown (in->SNSoln, err);
  if (err->error) Obit_traceback_msg (err, routine, in->RFIUV->name);
  in->SNSoln = ObitUVSolnUnref(in->SNSoln);
  in->SNTableRow = ObitTableSNRowUnref(in->SNTableRow);
  if (frRate) g_free(frRate);
  
} /* end of routine ObitUVRFIXizeFilter */ 

/**
 * Remove estimated RFI from data 
 * Read RFIUV data and apply inverse of zero fringe rate table phases
 * \param in     UVRFIXize object
 * Output parameters on info element:
 * \li "fractFlag"  OBIT_float (1,1,1) Fraction of data flagged
 * \li "fractMod"   OBIT_float (1,1,1) Fraction of data modified
 * \param err    Error/message stack, returns if error.
 */
void ObitUVRFIXizeCorrect (ObitUVRFIXize *in, ObitErr* err) 
{
  ObitIOCode retCode;
  gboolean doCalSelect=FALSE, done, blank, want=FALSE, someGood;
  ObitInfoType type;
  ObitIOAccess access;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong i, j, jndx, iindx, oindx, ant1, ant2;
  olong incs, incif, incf,  nif, nchan, nstok, iif, ichan, istok;
  olong BLindx, BLoff;
  odouble totalVis=1.0e-10; /* Total number of output visibilities sampled */
  odouble numFlag=0;        /* Number of visibilities flagged */
  odouble numMod=0;         /* Number of visibilities modified */
  ofloat fractFlag, fractMod;
  ObitUVDesc *inDesc, *outDesc;
  gchar *routine = "ObitUVRFIXizeCorrect";
  
  /* Error checks */
  if (err->error) return ;  /* previous error? */
  g_assert (ObitUVRFIXizeIsA(in));
  
  /* Make sure in->RFIUV available */
  Obit_return_if_fail ((in->RFIUV), err, 
		       "%s Averaged counterrotated data not available (run Counterrot)", 
		       routine);  

  /* Local pointers */
  inDesc  = in->myUV->myDesc;
  outDesc = in->outUV->myDesc;
  
  /* Check compatability */
  Obit_return_if_fail ((inDesc->lrec==outDesc->lrec), err, 
		       "%s Input and output data records incompatable", routine);  
  Obit_return_if_fail ((inDesc->lrec==in->RFIUV->myDesc->lrec), err, 
		       "%s Input and RFI model data records incompatable", routine);  

  /* How many IFs, freq, poln */
  nchan = inDesc->inaxes[inDesc->jlocf];
  if (inDesc->jlocf>=0) nif = inDesc->inaxes[inDesc->jlocif];
  else nif = 1;
  if (inDesc->jlocs>=0) nstok = inDesc->inaxes[inDesc->jlocs];
  else nstok = 1;
  
  /* Data always expanded to 3 words per vis */
  incs  = 3 * inDesc->incs  / inDesc->inaxes[0];
  incf  = 3 * inDesc->incf  / inDesc->inaxes[0];
  incif = 3 * inDesc->incif / inDesc->inaxes[0];
  
  /* Selection of input? */
  doCalSelect = TRUE;
  ObitInfoListGetTest(in->myUV->info, "doCalSelect", &type, dim, &doCalSelect);
  if (doCalSelect) access = OBIT_IO_ReadCal;
  else access = OBIT_IO_ReadOnly;

  /* Start up interpolator in estimate of RFI */
  ObitUVRFIXizeFetchStartUp (in, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  
  /* Open Data */
  retCode = ObitUVOpen (in->myUV, access, err);
  if ((retCode != OBIT_IO_OK) || (err->error>0)) 
    Obit_traceback_msg (err, routine, in->myUV->name);
  retCode = ObitUVOpen (in->outUV, OBIT_IO_ReadWrite, err);
  if ((retCode != OBIT_IO_OK) || (err->error>0)) 
    Obit_traceback_msg (err, routine, in->outUV->name);

  /* Loop over data */
  done = (retCode != OBIT_IO_OK);
  while (!done) {
    
    /* read buffer */
    retCode = ObitUVRead (in->myUV, NULL, err);
    if (err->error) Obit_traceback_msg (err, routine, in->myUV->name);
    
    done = (retCode == OBIT_IO_EOF); /* done? */
    if (done) break;
    
    /* Modify data */
    outDesc->numVisBuff = 0;  /* No data yet */
    for (i=0; i<inDesc->numVisBuff; i++) { /* loop over visibilities */
      jndx = i*inDesc->lrec;
      /* Update RFI model if necessary */
      if (in->myUV->buffer[jndx+inDesc->iloct] > in->VisTime){
	want = ObitUVRFIXizeFetch (in, in->myUV->buffer[jndx+inDesc->iloct], err);
	if (err->error) Obit_traceback_msg (err, routine, in->myUV->name);
      }

      /* This one selected? */
      if (!want) continue;

      /* Which baseline? */
      ant1 = (in->myUV->buffer[jndx+inDesc->ilocb] / 256.0) + 0.001;
      ant2 = (in->myUV->buffer[jndx+inDesc->ilocb] - ant1 * 256) + 0.001;
      /* Baseline index this assumes a1<a2 always - ignore auto correlations */
      if (ant1!=ant2) {
	BLoff = (in->blLookup[ant1-1] + ant2-ant1-1) * in->myUV->myDesc->lrec
	  + in->myUV->myDesc->nrparm;
      } else {
	continue;
      }

      /* Copy random parameters */
      someGood = FALSE;
      iindx = jndx;
      oindx = (outDesc->numVisBuff)*outDesc->lrec;
      for (j=0; j<outDesc->nrparm; j++) {
	in->outUV->buffer[oindx+j] = in->myUV->buffer[iindx+j];
      }
      /* Zero output */
      for (j=outDesc->nrparm; j<outDesc->lrec; j++) in->outUV->buffer[oindx+j] = 0.0;
      
      
      /* loop over IF */
      for (iif=0; iif<nif; iif++) {
	/* Loop over frequency channel */
	for (ichan=0; ichan<nchan; ichan++) { 
	  /* Loop over Stokes */
	  for (istok=0; istok<nstok; istok++) {
	    iindx = jndx + inDesc->nrparm + incs*istok + incf*ichan + incif*iif;
	    oindx = (outDesc->numVisBuff)*outDesc->lrec + outDesc->nrparm +
	      incs*istok + incf*ichan + incif*iif;
	    BLindx = BLoff + incs*istok + incf*ichan + incif*iif;
	    if (in->myUV->buffer[iindx+2]>0.0) {
	      totalVis++;   /* Count samples */
	      /* What to do?  is RFI model blanked? */
	      blank = in->VisApply[BLindx+2]<=0.0; 
	      if (blank) { /* Can't make correction - flag */
		numFlag++;   /* Count no. flagged */
		in->outUV->buffer[oindx+0] = 0.0;
		in->outUV->buffer[oindx+1] = 0.0;
		in->outUV->buffer[oindx+2] = 0.0;
	      } else { /* subtract correction */
		if ((in->VisApply[BLindx+0]!=0) && 
		    (in->VisApply[BLindx+1]!=0.0)) numMod++;   /* Count no. modified */
		someGood = TRUE;
		in->outUV->buffer[oindx+0] = 
		  in->myUV->buffer[iindx+0] - in->VisApply[BLindx+0];
		in->outUV->buffer[oindx+1] = 
		  in->myUV->buffer[iindx+1] - in->VisApply[BLindx+1];
		in->outUV->buffer[oindx+2] = in->myUV->buffer[iindx+2];
		/* DEBUG - just replace  
		in->outUV->buffer[oindx+0] = in->VisApply[BLindx+0];
		in->outUV->buffer[oindx+1] = in->VisApply[BLindx+1];
		in->outUV->buffer[oindx+2] = in->VisApply[BLindx+2];*/
		/* End DEBUG */
	    }
	    } /* end if correlation valid */
	  } /* end loop over polarization */
	} /* end loop over frequency */
      } /* end loop over IF */

      if (someGood) outDesc->numVisBuff++;  /* Count output */
    } /* End loop over buffer */
    
    /* Rewrite buffer */
    retCode = ObitUVWrite (in->outUV, NULL, err);
    if (err->error) Obit_traceback_msg (err, routine, in->outUV->name);
    
  } /* End loop over data */

  /* Close up */
  retCode = ObitUVClose(in->myUV, err);
  retCode = ObitUVClose(in->outUV, err);
  if (err->error) Obit_traceback_msg (err, routine, in->outUV->name);

  /* Give report */
  fractFlag = numFlag / totalVis;
  fractMod  = numMod  / totalVis;
  Obit_log_error(err, OBIT_InfoErr, 
		 "Flagged %5.1f percent, modified %5.1f percent",
		 100.0*fractFlag, 100.0*fractMod);
  dim[0] = dim[1] = dim[2] = dim[3] = dim[4] = 1;
  ObitInfoListAlwaysPut(in->info, "fractFlag", OBIT_float, dim, &fractFlag);
  ObitInfoListAlwaysPut(in->info, "fractMod",  OBIT_float, dim, &fractMod);

  /* Cleanup */
  ObitUVRFIXizeFetchShutDown (in, err);
  if (err->error) Obit_traceback_msg (err, routine, in->outUV->name);
} /* end of routine ObitUVRFIXizeCorrect */ 

/**
 * Initialize structures RFI estimation interpolation .
 * \param in   Solution Object.
 * \param err  ObitError stack.
 */
void ObitUVRFIXizeFetchStartUp (ObitUVRFIXize *in, ObitErr *err)
{
  ObitUVDesc *desc = in->RFIUV->myDesc;
  ofloat fblank = ObitMagicF();
  olong i, j, jndx, ant1, ant2;
  gint32 dim[MAXINFOELEMDIM]={1,1,1,1,1};
  gchar *interMode="AMBG";
  gchar *routine="ObitUVRFIXizeFetchStartUp";

  /* error checks */
  if (err->error) return;
  g_assert (ObitUVRFIXizeIsA(in));

  /* Open RFI Estimation UV data */
  ObitUVOpen(in->RFIUV, OBIT_IO_ReadOnly, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Copy descriptor information */
  in->numAnt    = desc->maxAnt;
  in->numBL     = in->numAnt*in->numAnt/2;
  in->numVis    = desc->nvis;

  /* Create baseline index lookup */
  in->blLookup = g_malloc0(in->numAnt*sizeof(olong));
  in->blLookup[0] = 0;
  for (i=1; i<in->numAnt; i++) {
    in->blLookup[i] = in->blLookup[i-1] + in->numAnt-i;
  }

  /* Nothing read yet */
  in->LastVisRead = 200000000;  /* Force read */
  in->NextVisRead = 1;          /* Start at beginning */

  /* In case restarting */
  if (in->PriorVisTime)  g_free (in->PriorVisTime);  in->PriorVisTime  = NULL;
  if (in->FollowVisTime) g_free (in->FollowVisTime); in->FollowVisTime = NULL;
  if (in->ApplyVisTime)  g_free (in->ApplyVisTime);  in->ApplyVisTime  = NULL;
  if (in->PriorVisNum)   g_free (in->PriorVisNum);   in->PriorVisNum   = NULL;
  if (in->FollowVisNum)  g_free (in->FollowVisNum);  in->FollowVisNum  = NULL;
  if (in->ApplyVisNum)   g_free (in->ApplyVisNum);   in->ApplyVisNum   = NULL;
  if (in->VisApply)      g_free (in->VisApply);      in->VisApply      = NULL;
  if (in->VisPrior)      g_free (in->VisPrior);      in->VisPrior      = NULL;
  if (in->VisFollow)     g_free (in->VisFollow);     in->VisFollow     = NULL;
  if (in->AntReal)       g_free (in->AntReal);       in->AntReal       = NULL;
  if (in->AntImag)       g_free (in->AntImag);       in->AntImag       = NULL;
  if (in->AntDelay)      g_free (in->AntDelay);      in->AntDelay      = NULL;
  if (in->AntRate)       g_free (in->AntRate);       in->AntRate       = NULL;

  /* Allocate arrays */
  in->lenVisArrayEntry = desc->lrec; /* length of vis array entry */
  in->VisApply     = g_malloc(desc->lrec*in->numBL*sizeof(ofloat));
  in->VisPrior     = g_malloc(desc->lrec*in->numBL*sizeof(ofloat));
  in->VisFollow    = g_malloc(desc->lrec*in->numBL*sizeof(ofloat));
  in->PriorVisTime = g_malloc(in->numBL*sizeof(ofloat));
  in->FollowVisTime= g_malloc(in->numBL*sizeof(ofloat));
  in->ApplyVisTime = g_malloc(in->numBL*sizeof(ofloat));
  in->PriorVisNum  = g_malloc(in->numBL*sizeof(olong));
  in->FollowVisNum = g_malloc(in->numBL*sizeof(olong));
  in->ApplyVisNum  = g_malloc(in->numBL*sizeof(olong));
  in->AntReal      = g_malloc(in->numAnt*in->numIF*sizeof(ofloat));
  in->AntImag      = g_malloc(in->numAnt*in->numIF*sizeof(ofloat));
  in->AntDelay     = g_malloc(in->numAnt*sizeof(ofloat));
  in->AntRate      = g_malloc(in->numAnt*sizeof(ofloat));

  /* Initial times to trigger update */
  in->VisTime    = -1.0e20;
  in->ReadTime   = -1.0e20;
  in->PriorTime  = -1.0e20;
  in->FollowTime = -1.0e20;
  
  /* initialize Prior and Following arrays */
  /* Loop over antenna 1 */
  for (ant1=0; ant1<in->numAnt-1; ant1++) {
    in->AntRate[ant1] = fblank;
    in->AntDelay[ant1]= fblank;
    /* Loop over antenna 2 */
    for (ant2=ant1+1; ant2<in->numAnt; ant2++) {
      jndx = in->blLookup[ant1] + ant2-ant1-1;
      in->PriorVisTime[jndx]  = -10000.0;
      in->FollowVisTime[jndx] = -10000.0;
      in->ApplyVisTime[jndx]  = -10000.0;
      in->PriorVisNum[jndx]   = -1;
      in->FollowVisNum[jndx]  = -1;
      in->ApplyVisNum[jndx]   = -1;
      for (j=0; j<in->lenVisArrayEntry; j++) in->VisPrior[jndx+j]  = 0.0;
      for (j=0; j<in->lenVisArrayEntry; j++) in->VisFollow[jndx+j] = 0.0;
    } /* end loop over ant2 */
  } /* end loop over ant1 */
  in->AntRate[ant1] = fblank;
  in->AntDelay[ant1]= fblank;
  for (i=0; i<in->numAnt*in->numIF; i++) {
    in->AntReal[i] = fblank;
    in->AntImag[i] = fblank;
  }
  
  /* Create/start uvSoln object */
  in->SNSoln = ObitUVSolnCreate(in->name, in->residUV);
  dim[0] = strlen(interMode);
  ObitInfoListAlwaysPut(in->SNSoln->info,"interMode", OBIT_string, dim, interMode);
  ObitUVSolnStartUp (in->SNSoln, err);
  if (err->error) Obit_traceback_msg (err, routine, in->RFIUV->name);

  /* Setup SN table rows if needed */
  if (in->SNTableRow==NULL) in->SNTableRow = newObitTableSNRow(in->SNSoln->SNTable);
} /*  end ObitUVRFIXizeFetchStartUp */

/**
 * Get RFI estimation for a given time;
 * Interpolates entries closest in time
 * Values for time are in array in->VisApply which is an array
 * of full visibility records, one per baseline.
 * The first baseline for antenna n (baseline n - n+1) starts in
 * slot number in->blLookup[n-1].
 * Apply inverse of zero fringe rate calibration phases.
 * \param in    UVRFIXize Object.
 * \param time  Desired time (days)
 * \param err   ObitError stack.
 * \return TRUE if source, time etc. selected, else FALSE.
 */
gboolean ObitUVRFIXizeFetch (ObitUVRFIXize *in, ofloat time, ObitErr *err)
{
  gboolean out = FALSE;
  gboolean wanted, doPrior, good1, good2;
  olong  k, j, jndx, indx, kndx, BLindx, ant1, ant2;
  olong incs, incif, incf,  nif, nchan, nstok, iif, ichan, istok;
  ofloat gr1, gi1, dly1, rat1, gr2, gi2, dly2, rat2, tr, ti, gr, gi;
  ofloat dt1, dt2, dt, wt1, wt2, amp1, amp2, phase1, phase2, amp, phase=0.0;
  ofloat freqFact, dc, dcr, dci, ggr, ggi;
  ofloat fblank = ObitMagicF();
  gdouble twopi = 2.0* G_PI;
  ObitUVDesc *desc;
  ObitUVSel *sel;
  gchar *routine="ObitUVRFIXizeFetch";

  /* error checks */
  if (err->error) return out;

   /* Save time */
  in->VisTime = time;

 /* local pointers for structures */
  desc = in->myUV->myDesc;
  sel  = in->myUV->mySel;

  /* How many IFs, freq, poln */
  nchan = desc->inaxes[desc->jlocf];
  if (desc->jlocf>=0) nif = desc->inaxes[desc->jlocif];
  else nif = 1;
  if (desc->jlocs>=0) nstok = desc->inaxes[desc->jlocs];
  else nstok = 1;

  /* Data always expanded to 3 words per vis */
  incs  = 3 * desc->incs  / desc->inaxes[0];
  incf  = 3 * desc->incf  / desc->inaxes[0];
  incif = 3 * desc->incif / desc->inaxes[0];

  /* Check if this data wanted */
  wanted = (time>=sel->timeRange[0]) && (time<=sel->timeRange[1]);
  if (!wanted) return wanted;
  out = wanted;

  /* see if new time - update arrays. */
  if (time > in->ReadTime) {
    ObitUVRFIXizeNewTime (in, time, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, out);
  }

  /* Rate frequency factor */
  freqFact = twopi * 86400.0;

  /* Update zero delay/fringe rate table if necessary 
   leaves results in in->AntReal, in->AntImag, in->AntRate and in->AntDelay */
  if (time > in->SNSoln->CalTime){
    getAntDelFR (in, in->SNSoln, time, 
		 in->AntReal, in->AntImag, in->AntDelay, in->AntRate, err);
    if (err->error) Obit_traceback_val (err, routine, in->SNSoln->name, out);
  }

  /* Loop over antenna 1 */
  for (ant1=0; ant1<in->numAnt-1; ant1++) {
    /* Loop over antenna 2 */
    for (ant2=ant1+1; ant2<in->numAnt; ant2++) {
      jndx = in->blLookup[ant1] + ant2-ant1-1;
      BLindx = jndx*desc->lrec;

      /* set interpolation weights proportional to time difference. */
      dt  = in->FollowVisTime[jndx] - in->PriorVisTime[jndx] + 1.0e-20;
      dt1 = time - in->PriorVisTime[jndx];
      dt2 = time - in->FollowVisTime[jndx];
      wt1 = 0.0;
      if ((time < in->FollowVisTime[jndx]) && (dt>1.0e-19)) wt1 = -dt2 / dt;
      /* Before Prior time? */
      if ((time < in->PriorVisTime[jndx]) || (in->FollowVisTime[jndx]<-1000.0)) wt1 = 1.0;
      /* After Follow time? */
      if ((time > in->FollowVisTime[jndx]) && (in->FollowVisTime[jndx]<-1000.0)) wt1 = 0.0;
      wt2 = 1.0 - wt1;

      /* Is prior or follow closer in time */
      doPrior = fabs(time-in->PriorVisTime[jndx]) <= fabs(time-in->FollowVisTime[jndx]);

      /* Copy closest visibility random parameters */
      if (doPrior) {
	in->ApplyVisTime[jndx] = in->PriorVisTime[jndx];
	in->ApplyVisNum[jndx]  = in->PriorVisNum[jndx];
	for (j=0; j<desc->nrparm; j++) 
	  in->VisApply[BLindx+j] = in->VisPrior[BLindx+j];
      } else { /* Follow closer */
	in->ApplyVisTime[jndx] = in->FollowVisTime[jndx];
	in->ApplyVisNum[jndx]  = in->FollowVisNum[jndx];
	for (j=0; j<desc->nrparm; j++) 
	  in->VisApply[BLindx+j] = in->VisFollow[BLindx+j];
      }

      /* Baseline phase calibration info */
      dly1 = in->AntDelay[ant1];
      dly2 = in->AntDelay[ant2];
      rat1 = in->AntRate[ant1];
      rat2 = in->AntRate[ant2];
      
      /* Interpolate visibilities */
      /* Loop over IF */
      for (iif=0; iif<in->numIF; iif++) {
	/* Loop over poln */
	for (istok=0; istok<nstok; istok++) {
	  /* Loop over channel */
	  for (ichan=0; ichan<nchan; ichan++) {
	    /* Index for this datum */
	    kndx = BLindx + desc->nrparm + incs*istok + incf*ichan + incif*iif;
	    
	    /*for (j=desc->nrparm; j<desc->lrec; j+=3) {*/
	    good1 = (in->PriorVisTime[jndx]>-100.0)  && (in->VisPrior[kndx+2]>0.0)  
	      && (wt1>0.0) && (rat1!=fblank);
	    good2 = (in->FollowVisTime[jndx]>-100.0) && (in->VisFollow[kndx+2]>0.0) 
	      && (wt2>0.0) && (rat2!=fblank);
	    
	    /* Prior good, follow bad */
	    if ((good1) && (!good2)) {
	      in->VisApply[kndx]   = in->VisPrior[kndx];
	      in->VisApply[kndx+1] = in->VisPrior[kndx+1];
	      in->VisApply[kndx+2] = in->VisPrior[kndx+2];
	      
	      /* Prior bad, follow good */
	    } else if ((!good1) && (good2)) {
	      in->VisApply[kndx]   = in->VisFollow[kndx];
	      in->VisApply[kndx+1] = in->VisFollow[kndx+1];
	      in->VisApply[kndx+2] = in->VisFollow[kndx+2];
	      
	      /* Both bad - flag */
	    } else if ((!good1) && (!good2)) {
	      in->VisApply[kndx]   = 0.0;
	      in->VisApply[kndx+1] = 0.0;
	      in->VisApply[kndx+2] = 0.0;
	      
	      /* Both good, interpolate in amp phase */
	    } else {
	      amp1 = sqrt (in->VisPrior[kndx]*in->VisPrior[kndx] + 
			   in->VisPrior[kndx+1]*in->VisPrior[kndx+1]);
	      amp2 = sqrt (in->VisFollow[kndx]*in->VisFollow[kndx] + 
			   in->VisFollow[kndx+1]*in->VisFollow[kndx+1]);
	      phase1 = atan2 (in->VisPrior[kndx+1],  in->VisPrior[kndx]+1.0e-20);
	      phase2 = atan2 (in->VisFollow[kndx+1], in->VisFollow[kndx]+1.0e-20);
	      /* Keep phase constant when interpolating with one zero. */
	      if ((amp1!=0.) && (amp2==0.)) {
		phase = phase1;
	      } else if ((amp1==0.) && (amp2!=0.)) {
		phase = phase2;
	      } else {
		/* Force within half turn */
		if ((phase2-phase1)>=G_PI) phase2 -= twopi;
		if ((phase2-phase1)<-G_PI) phase2 += twopi;
		phase = wt1*phase1 + wt2*phase2;
	      }
	      amp   = wt1*amp1 + wt2*amp2;
	      in->VisApply[kndx]   = amp * cos(phase);
	      in->VisApply[kndx+1] = amp * sin(phase);
	      in->VisApply[kndx+2] = in->VisPrior[kndx+2]*wt1 + in->VisFollow[kndx+2]*wt2;
	    }
	  } /* end loop over channel */
	} /* end poln loop */
      } /* end loop over IF */

      /* Calibrate visibility in VisApply[jndx] */
      /* Valid cal? */
      if ((dly1!=fblank) && (dly2!=fblank)) {
	/* Loop over IF */
	for (iif=0; iif<in->numIF; iif++) {
	  indx = ant1*in->numIF+iif;
	  gr1 = in->AntReal[indx];
	  gi1 = in->AntImag[indx];
	  indx = ant2*in->numIF+iif;
	  gr2 = in->AntReal[indx];
	  gi2 = in->AntImag[indx];
	  gr = gr1*gr2 + gi1*gi2;
	  gi = gr2*gi1 - gr1*gi2;

	  /* Loop over poln */
	  for (istok=0; istok<nstok; istok++) {
	    /* Loop over channel */
	    for (ichan=0; ichan<nchan; ichan++) {
	      kndx = BLindx + desc->nrparm + incs*istok + incf*ichan + incif*iif;
	      /* Delay correction - num. turns = freq diff from ref channel * delay */
	      dc = -twopi * (dly1-dly2) * desc->chIncIF[iif] * 
		(ichan - desc->crpix[desc->jlocf] + 1.0);
	      dcr = cos(dc);
	      dci = sin(dc);
	      ggr = gr*dcr + gi*dci;
	      ggi = gr*dci - gi*dcr;
	      tr = in->VisApply[kndx+0];
	      ti = in->VisApply[kndx+1];
	      in->VisApply[kndx+0] = tr*ggr - ti*ggi;
	      in->VisApply[kndx+1] = ggr*ti + tr*ggi;
	    } /* end loop over channel */
	  } /* end poln loop */
	} /* end IF loop */
      } else { /* bad cal - flag */
	for (k=0; k<desc->ncorr; k++) {
	  kndx = BLindx + desc->nrparm + k*3;
	  in->VisApply[kndx+2] = 0.0;
	}
      } /* end bad cal */
    } /* end loop over antenna 2 */
  } /* end loop over antenna 1 */

  return out;
} /* end ObitUVRFIXizeFetch */


/**
 * Shutdown Solution interpolation.
 * Close any open file and destroy structures.
 * \param in   Calibrate Object.
 * \param err  ObitError stack.
 */
void ObitUVRFIXizeFetchShutDown(ObitUVRFIXize *in, ObitErr *err)
{
  ObitIOCode retCode;
  gchar *routine="ObitUVRFIXizeFetchShutDown";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitUVRFIXizeIsA(in));

  /* Close UV data  */
  retCode = ObitUVClose (in->RFIUV, err);
  if ((retCode!=OBIT_IO_OK) || (err->error))
    Obit_traceback_msg (err, routine, in->name);

  ObitUVSolnShutDown (in->SNSoln, err);
  if (err->error) Obit_traceback_msg (err, routine, in->RFIUV->name);
  in->SNSoln     = ObitUVSolnUnref(in->SNSoln);
  in->SNTableRow = ObitTableSNRowUnref(in->SNTableRow);
  if (in->PriorVisTime)  g_free (in->PriorVisTime);  in->PriorVisTime  = NULL;
  if (in->FollowVisTime) g_free (in->FollowVisTime); in->FollowVisTime = NULL;
  if (in->ApplyVisTime)  g_free (in->ApplyVisTime);  in->ApplyVisTime  = NULL;
  if (in->PriorVisNum)   g_free (in->PriorVisNum);   in->PriorVisNum   = NULL;
  if (in->FollowVisNum)  g_free (in->FollowVisNum);  in->FollowVisNum  = NULL;
  if (in->ApplyVisNum)   g_free (in->ApplyVisNum);   in->ApplyVisNum   = NULL;
  if (in->VisApply)      g_free (in->VisApply);      in->VisApply      = NULL;
  if (in->VisPrior)      g_free (in->VisPrior);      in->VisPrior      = NULL;
  if (in->VisFollow)     g_free (in->VisFollow);     in->VisFollow     = NULL;
  if (in->AntReal)       g_free (in->AntReal);       in->AntReal       = NULL;
  if (in->AntImag)       g_free (in->AntImag);       in->AntImag       = NULL;
  if (in->AntDelay)      g_free (in->AntDelay);      in->AntDelay      = NULL;
  if (in->AntRate)       g_free (in->AntRate);       in->AntRate       = NULL;

} /*  end ObitUVRFIXizeFetchShutDown */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitUVRFIXizeClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitUVRFIXizeClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitUVRFIXizeClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitUVRFIXizeClassInfoDefFn (gpointer inClass)
{
  ObitUVRFIXizeClassInfo *theClass = (ObitUVRFIXizeClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitUVRFIXizeClassInit;
  theClass->newObit       = (newObitFP)newObitUVRFIXize;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitUVRFIXizeClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitUVRFIXizeGetClass;
  theClass->ObitCopy      = NULL;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitUVRFIXizeClear;
  theClass->ObitInit      = (ObitInitFP)ObitUVRFIXizeInit;
  theClass->ObitUVRFIXizeCounterRot   = (ObitUVRFIXizeCounterRotFP)ObitUVRFIXizeCounterRot;
  theClass->ObitUVRFIXizeFilter       = (ObitUVRFIXizeFilterFP)ObitUVRFIXizeFilter;
  theClass->ObitUVRFIXizeCorrect      = (ObitUVRFIXizeCorrectFP)ObitUVRFIXizeCorrect;
  theClass->ObitUVRFIXizeCreate       = (ObitUVRFIXizeCreateFP)ObitUVRFIXizeCreate;
  theClass->ObitUVRFIXizeFetchStartUp = 
    (ObitUVRFIXizeFetchStartUpFP)ObitUVRFIXizeFetchStartUp;
  theClass->ObitUVRFIXizeFetch        = (ObitUVRFIXizeFetchFP)ObitUVRFIXizeFetch;
  theClass->ObitUVRFIXizeFetchShutDown = 
    (ObitUVRFIXizeFetchShutDownFP)ObitUVRFIXizeFetchShutDown;

} /* end ObitUVRFIXizeClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitUVRFIXizeInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitUVRFIXize *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->thread        = newObitThread();
  in->info          = newObitInfoList(); 
  in->myUV          = NULL;
  in->SNSoln        = NULL;
  in->SNTable       = NULL;
  in->SNTableRow    = NULL;
  in->PriorVisTime  = NULL;
  in->FollowVisTime = NULL;
  in->ApplyVisTime  = NULL;
  in->PriorVisNum   = NULL;
  in->FollowVisNum  = NULL;
  in->ApplyVisNum   = NULL;
  in->VisApply      = NULL;
  in->VisPrior      = NULL;
  in->VisFollow     = NULL;
  in->AntReal       = NULL;
  in->AntImag       = NULL;
  in->AntDelay      = NULL;
  in->AntRate       = NULL;
  in->blLookup      = NULL;
  in->LastVisRead   = -1;
  in->NextVisRead   = -1;
  in->numAnt        = 0;
  in->CurSourID     = -1;
  in->PriorSourID   = -1;
  in->FollowSourID  = -1;
  in->PriorTime     = -1.0;
  in->FollowTime    = -1.0;
  in->VisTime       = -1.0;
  in->ReadTime      = -1.0;

} /* end ObitUVRFIXizeInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitUVRFIXize* cast to an Obit*.
 */
void ObitUVRFIXizeClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitUVRFIXize *in  = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  in->thread    = ObitThreadUnref(in->thread);
  in->info      = ObitInfoListUnref(in->info);

  in->SNSoln        = ObitUVSolnUnref(in->SNSoln);
  in->SNTable       = ObitTableSNUnref(in->SNTable);
  in->SNTableRow    = ObitTableSNRowUnref(in->SNTableRow);
  in->myUV          = ObitUVUnref(in->myUV);
  in->RFIUV         = ObitUVUnref(in->RFIUV);
  in->outUV         = ObitUVUnref(in->outUV);
  if (in->PriorVisTime)  g_free (in->PriorVisTime);  in->PriorVisTime  = NULL;
  if (in->FollowVisTime) g_free (in->FollowVisTime); in->FollowVisTime = NULL;
  if (in->ApplyVisTime)  g_free (in->ApplyVisTime);  in->ApplyVisTime  = NULL;
  if (in->PriorVisNum)   g_free (in->PriorVisNum);   in->PriorVisNum   = NULL;
  if (in->FollowVisNum)  g_free (in->FollowVisNum);  in->FollowVisNum  = NULL;
  if (in->ApplyVisNum)   g_free (in->ApplyVisNum);   in->ApplyVisNum   = NULL;
  if (in->VisApply)      g_free (in->VisApply);      in->VisApply      = NULL;
  if (in->VisPrior)      g_free (in->VisPrior);      in->VisPrior      = NULL;
  if (in->VisFollow)     g_free (in->VisFollow);     in->VisFollow     = NULL;
  if (in->AntReal)       g_free (in->AntReal);       in->AntReal       = NULL;
  if (in->AntImag)       g_free (in->AntImag);       in->AntImag       = NULL;
  if (in->AntDelay)      g_free (in->AntDelay);      in->AntDelay      = NULL;
  if (in->AntRate)       g_free (in->AntRate);       in->AntRate       = NULL;
  if (in->blLookup)      g_free (in->blLookup);      in->blLookup      = NULL;
 
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitUVRFIXizeClear */

/**
 * Fill array of antenna fringe rates from a ObitUVSoln at a given time
 * \param in      UVRFIXize object
 * \param SNSoln  ObitUVSoln from which to extract fringe rates
 * \param time    desired time in days
 * \param frRate  [out] Fringe rates per antenna in sec/sec
 *                blanked if not given.
 * \param err     Error stack for messages and errors.
 */
static void getAntFR (ObitUVRFIXize *in, ObitUVSoln *SNSoln, ofloat time, 
		      ofloat *frRate, ObitErr* err)
{
  ofloat fblank = ObitMagicF();
  olong iant;
  gchar *routine = "ObitUVRFIXize:getAntFR";

  /* Setup SN table row if needed */
  if (in->SNTableRow==NULL) in->SNTableRow = newObitTableSNRow(SNSoln->SNTable);

  /* Set time */
  in->SNTableRow->Time = time;

  /* see if new time - update tables. */
  if (in->SNTableRow->Time > SNSoln->CalTime) {
    ObitUVSolnUpdate (SNSoln, in->SNTableRow->Time, in->SNTableRow->SourID, err);
    if (err->error) Obit_traceback_msg (err, routine, SNSoln->name);
  }

  /* Loop over antennas */
  for (iant=0; iant<SNSoln->numAnt; iant++) {
    in->SNTableRow->antNo = iant+1;
    if (ObitUVSolnGetSN (SNSoln, in->SNTableRow, err)) {
      /* Should all be the same */
      frRate[iant] = in->SNTableRow->Rate1[0];
    } else { /* not there */
      frRate[iant] = fblank;
    }
    if (err->error) Obit_traceback_msg (err, routine, SNSoln->name);
  } /* end loop over antennas */
} /* end getAntFR */

/**
 * Fill array of antenna fringe rates from a ObitUVSoln at a given time
 * Returns inverse of calibration in table (amplitudes=1.0)
 * \param in      UVRFIXize object
 * \param SNSoln  ObitUVSoln from which to extract fringe rates
 * \param time    desired time in days
 * \param Real    [out] real part of gain per IF and antenna
 * \param Imag    [out] imaginary part of gain per IF and antenna
 * \param delay   [out] Group delay per antenna in sec
 *                blanked if not given.
 * \param frRate  [out] Fringe rates per antenna in sec/sec
 *                blanked if not given.
 * \param err     Error stack for messages and errors.
 */
static void getAntDelFR (ObitUVRFIXize *in, ObitUVSoln *SNSoln, ofloat time, 
			 ofloat *Real, ofloat *Imag, ofloat *delay, ofloat *frRate,  
			 ObitErr* err)
{
  ofloat fblank = ObitMagicF();
  olong iant, iif, jndx;
  gchar *routine = "ObitUVRFIXize:getAntDelFR";

  /* Setup SN table row if needed */
  if (in->SNTableRow==NULL) in->SNTableRow = newObitTableSNRow(SNSoln->SNTable);

  /* Set time */
  in->SNTableRow->Time = time;

  /* see if new time - update tables. */
  if (in->SNTableRow->Time > SNSoln->CalTime) {
    ObitUVSolnUpdate (SNSoln, in->SNTableRow->Time, in->SNTableRow->SourID, err);
    if (err->error) Obit_traceback_msg (err, routine, SNSoln->name);
  }

  /* Loop over antennas */
  for (iant=0; iant<SNSoln->numAnt; iant++) {
    in->SNTableRow->antNo = iant+1;
    if (ObitUVSolnGetSN (SNSoln, in->SNTableRow, err)) {
      /* Should all be the same */
      delay[iant]  = -in->SNTableRow->Delay1[0];
      frRate[iant] = -in->SNTableRow->Rate1[0];
      for (iif=0; iif<in->numIF; iif++) {
	jndx = iant*in->numIF + iif;
	Real[jndx] =  in->SNTableRow->Real1[iif];
	Imag[jndx] = -in->SNTableRow->Imag1[iif];
      }
    } else { /* not there */
      delay[iant]  = fblank;
      frRate[iant] = fblank;
      for (iif=0; iif<in->numIF; iif++) {
	jndx = iant*in->numIF + iif;
	Real[jndx] = fblank;
	Imag[jndx] = fblank;
      }
    }
    if (err->error) Obit_traceback_msg (err, routine, SNSoln->name);
  } /* end loop over antennas */
} /* end getAntDelFR */

/**
 * Read visibilities for next time from in->RFIUV.
 * Updates following time and shuffles to prior.
 * \param in   UVRFIXize object.
 * \param time Desired time in days
 * \param err  Error stack for messages and errors.
 */
static void ObitUVRFIXizeNewTime (ObitUVRFIXize *in, ofloat time,
				  ObitErr *err)
{
  ObitIOCode retCode;
  ofloat sum, endTime;
  olong  j, ant1, ant2, count, indx, jndx, BLoff;
  gboolean done, readAll=FALSE;
  ObitUVDesc *inDesc = in->RFIUV->myDesc;
  gchar *routine="ObitUVRFIXizeNewTime";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;

  /* Save time */
  in->ReadTime = time;

  /* Shuffle data from Following to Prior if time exceeded */
  /* Loop over antenna 1 */
  for (ant1=0; ant1<in->numAnt-1; ant1++) {
    /* Loop over antenna 2 */
    for (ant2=ant1+1; ant2<in->numAnt; ant2++) {
      jndx = in->blLookup[ant1] + ant2-ant1-1;
      if ((time > in->FollowVisTime[jndx])  &&  (in->FollowVisTime[jndx] > -100.)) {
	in->PriorVisTime[jndx] = in->FollowVisTime[jndx];
	in->PriorVisNum[jndx]  = in->FollowVisNum[jndx];
	for (j=0; j<in->lenVisArrayEntry; j++) in->VisPrior[jndx+j] = in->VisFollow[jndx+j];
	/*in->FollowVisTime[jndx] = -10000.0;
	  in->FollowVisNum[jndx ] = -1;*/
      } /* end if shuffle */
    } /* end ant 2 loop */
  } /* end ant 1 loop */

  /* Is the desired next datum (in->NextVisRead in the current buffer or later? */
  in->LastVisRead = in->NextVisRead;
  if (in->NextVisRead<inDesc->firstVis) {
    /* No - Restart IO */
    ObitUVIOReset (in->RFIUV, in->NextVisRead, err);
    if (err->error) Obit_traceback_msg (err, routine, in->RFIUV->name);
    in->LastVisRead = 200000000;  /* Force read */
  }

  endTime = time + 1.2*in->AvgTime; /* Don't search further than this time */

  /* Read  filling in data  until after endTime. */
  readAll = in->NextVisRead>=inDesc->nvis;  /* At end? */
  done = readAll;  /* If all read, leave as is */
  while (!done) {

    /* Need another buffer load? */
    if (in->LastVisRead>=(inDesc->firstVis+inDesc->numVisBuff-1)) {
      retCode = ObitUVRead(in->RFIUV, NULL, err);
      if (err->error) Obit_traceback_msg (err, routine, in->RFIUV->name);
      in->LastVisRead = inDesc->firstVis;
      readAll = retCode==OBIT_IO_EOF; /* Read whole file? */
      done = readAll;
    }
    /* Is time past end of search window? */
    indx  = (in->LastVisRead-inDesc->firstVis) * inDesc->lrec;
    if (in->RFIUV->buffer[indx+inDesc->iloct] >= endTime) break;

    /* Which baseline? */
    ant1  = (in->RFIUV->buffer[indx+inDesc->ilocb] / 256.0) + 0.001;
    ant2  = (in->RFIUV->buffer[indx+inDesc->ilocb] - ant1 * 256) + 0.001;

    /* Trap bad data */
    if ((ant1<1) || (ant2<1) || (ant1>in->numAnt) || (ant2>in->numAnt)) goto skip;

    jndx  = in->blLookup[ant1-1] + ant2-ant1-1;
    BLoff = jndx * inDesc->lrec;

    /* If follow time already passed time ignore */
    if ((in->RFIUV->buffer[indx+inDesc->iloct]>in->FollowVisTime[jndx]) &&
	(in->FollowVisTime[jndx]>time)) goto skip;

    /* data time after target? */
    if (in->RFIUV->buffer[indx+inDesc->iloct] >= time) { 
      /* After target time update Follow */
      /* new following entry - copy to prior */
      if (in->FollowVisTime[jndx]<=time) {
	in->PriorVisTime[jndx] = in->FollowVisTime[jndx];
	in->PriorVisNum[jndx]  = in->FollowVisNum[jndx];
	for (j=0; j<in->lenVisArrayEntry; j++) in->VisPrior[BLoff+j] = in->VisFollow[BLoff+j];
      }

      /* fill in new following values */
      in->FollowVisTime[jndx] = in->RFIUV->buffer[indx+inDesc->iloct];
      in->FollowVisNum[jndx ] = in->LastVisRead;
      for (j=0; j<in->lenVisArrayEntry; j++) in->VisFollow[BLoff+j] = in->RFIUV->buffer[indx+j];

      /* if Prior entry not valid copy following */
      if (in->PriorVisTime[jndx] <= -100.) {
	in->PriorVisTime[jndx] = in->FollowVisTime[jndx];
	in->PriorVisNum[jndx]  = in->FollowVisNum[jndx];
	for (j=0; j<in->lenVisArrayEntry; j++) in->VisPrior[BLoff+j] = in->VisFollow[BLoff+j];
      }
    } else {  /* Prior to target time, update Prior */
      
      in->PriorVisTime[jndx] = in->RFIUV->buffer[indx+inDesc->iloct];
      in->PriorVisNum[jndx]  = in->LastVisRead;
      for (j=0; j<in->lenVisArrayEntry; j++) in->VisPrior[BLoff+j] = in->RFIUV->buffer[indx+j];
      /* No valid Follow time - use prior */
      in->FollowVisTime[jndx] = in->PriorVisTime[jndx];
      in->FollowVisNum[jndx ] = in->PriorVisNum[jndx];
      for (j=0; j<in->lenVisArrayEntry; j++) in->VisFollow[BLoff+j] = in->VisPrior[BLoff+j];
   } /* end update follow or prior */
  skip:
    in->LastVisRead++;  /* Number in buffer processed */
  } /* end loop over file  */
  
  /* Average prior and follow times */
  sum = 0.0; count = 0;
  in->NextVisRead = 1;  /* Restart looking at last prior */
  /* Loop over antenna 1 */
  for (ant1=0; ant1<in->numAnt-1; ant1++) {
    /* Loop over antenna 2 */
    for (ant2=ant1+1; ant2<in->numAnt; ant2++) {
      jndx = in->blLookup[ant1] + ant2-ant1-1;
      /* Find last valid prior */
      in->NextVisRead = MAX (in->NextVisRead, in->PriorVisNum[jndx]);
      if (in->PriorVisTime[jndx] > 0.0) {
	sum += in->PriorVisTime[jndx];
	count++;
      }
    } /* end ant 2 loop */
  } /* end ant 1 loop */
  if (count>0) in->PriorTime = sum/count;
  
  sum = 0.0; count = 0;
  /* Loop over antenna 1 */
  for (ant1=0; ant1<in->numAnt-1; ant1++) {
    /* Loop over antenna 2 */
    for (ant2=ant1+1; ant2<in->numAnt; ant2++) {
      jndx = in->blLookup[ant1] + ant2-ant1-1;
      if (in->FollowVisTime[jndx] > 0.0) {
	sum += in->FollowVisTime[jndx];
	count++;
      }
    } /* end ant 2 loop */
  } /* end ant 1 loop */
  if (count>0) in->FollowTime = sum/count;
 
  /* just to be sure something rational in times */
  if (in->PriorTime < -1000.0)  in->PriorTime  = time - 2.0/86400.0;
  if (in->FollowTime > 10000.0) in->FollowTime = time + 2.0/86400.0;
  
} /* end ObitUVRFIXizeNewTime */

