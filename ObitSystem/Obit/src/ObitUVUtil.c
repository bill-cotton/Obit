/* $Id$   */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2004-2022                                          */
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

#include "ObitUVDesc.h"
#include "ObitUVUtil.h"
#include "ObitTableSUUtil.h"
#include "ObitTableNXUtil.h"
#include "ObitTablePSUtil.h"
#include "ObitTableFQUtil.h"
#include "ObitTableANUtil.h"
#include "ObitTableFG.h"
#include "ObitPrecess.h"
#include "ObitUVWCalc.h"
#include "ObitUVSortBuffer.h"
#if HAVE_GSL==1  /* GSL stuff */
#include <gsl/gsl_randist.h>
#endif /* HAVE_GSL */
/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitUVUtil.c
 * ObitUVUtil module function definitions.
 */

/*---------------Private function prototypes----------------*/
/** Modify output descriptor from effects of frequency averaging */
static ofloat AvgFSetDesc (ObitUVDesc *inDesc, ObitUVDesc *outDesc, 
			   olong NumChAvg, olong *ChanSel, gboolean doAvgAll, 
			   olong *corChan, olong *corIF, olong *corStok, gboolean *corMask,
			   ObitErr *err);

/** Modify output descriptor from effects of frequency bloating */
static ofloat BloatFSetDesc (ObitUVDesc *inDesc, ObitUVDesc *outDesc, 
			     olong nBloat, ObitErr *err);

/** Average visibility in frequency */
static void AvgFAver (ObitUVDesc *inDesc, ObitUVDesc *outDesc, 
		      olong NumChAvg, olong *ChanSel, gboolean doAvgAll, 
		      olong *corChan, olong *corIF, olong *corStok, gboolean *corMask,
		      ofloat *inBuffer, ofloat *outBuffer, ofloat *work, 
		      ObitErr *err);

/** Smooth visibility in frequency */
static void SmooF (ObitUVDesc *inDesc, ObitUVDesc *outDesc, olong NumChSmo, 
		   olong *corChan, olong *corIF, olong *corStok, gboolean *corMask,
		   ofloat *inBuffer, ofloat *outBuffer, ofloat *work, 
		   ObitErr *err);

/** Hann visibility in frequency */
static void Hann (ObitUVDesc *inDesc, ObitUVDesc *outDesc, gboolean doDescm,
		  ofloat *inBuffer, ofloat *outBuffer, ofloat *work, 
		  ObitErr *err);

/** Duplicate channels */
static void Bloat (ObitUVDesc *inDesc, ObitUVDesc *outDesc, olong nBloat,
		   gboolean unHann, ofloat *inBuffer, ofloat *outBuffer, 
		   ObitErr *err);

/** Copy selected channels */
static void FreqSel (ObitUVDesc *inDesc, ObitUVDesc *outDesc, 
		     olong BChan, olong EChan, olong chinc,
		     olong BIF, olong EIF,
		     ofloat *inBuffer, ofloat *outBuffer);

/** Update FQ table for averaging */
static void FQSel (ObitUV *inUV, olong chAvg, olong fqid, ObitErr *err);

/** Low accuracy inverse Sinc function */
static ofloat InvSinc(ofloat arg);
/*----------------------Public functions---------------------------*/

/**
 * Find maximum baseline length and W in a data set.
 * Imaging parameters are on the inUV info member as arrays for a number 
 * of fields.
 * \param inUV     Input uv data. 
 * \param MaxBL    Output maximum baseline length (sqrt(u*u+v*v))
 * \param MaxW     Output Max abs(w) in data.
 * \param err      Error stack, returns if not empty.
 */
void ObitUVUtilUVWExtrema (ObitUV *inUV, ofloat *MaxBL, ofloat *MaxW, ObitErr *err)
{
  ObitIOCode retCode=OBIT_IO_SpecErr;
  olong i;
  ofloat bl, maxbl, maxw, *u, *v, *w;
  gchar *routine = "ObitUVUtilUVWExtrema";
 
   /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitUVIsA(inUV));

  /* See if values are already in the descriptor */
  if ((inUV->myDesc->maxBL>0.0) && (inUV->myDesc->maxW>0.0)) {
    *MaxBL = inUV->myDesc->maxBL;
    *MaxW  = inUV->myDesc->maxW;
    return;
  }

  /* Open uv data if not already open */
  if (inUV->myStatus==OBIT_Inactive) {
    retCode = ObitUVOpen (inUV, OBIT_IO_ReadOnly, err);
    if (err->error) Obit_traceback_msg (err, routine, inUV->name);
  }

  /* Loop through data */
  maxbl = maxw = 0.0;
  while (retCode==OBIT_IO_OK) {
    /* read buffer full */
    retCode = ObitUVRead (inUV, NULL, err);
    if (err->error) Obit_traceback_msg (err, routine, inUV->name);
    
    /* initialize data pointers */
    u   = inUV->buffer+inUV->myDesc->ilocu;
    v   = inUV->buffer+inUV->myDesc->ilocv;
    w   = inUV->buffer+inUV->myDesc->ilocw;
    for (i=0; i<inUV->myDesc->numVisBuff; i++) { /* loop over buffer */
      
      /* Get statistics */
      bl = sqrt ((*u)*(*u) + (*v)*(*v));
      maxbl = MAX (maxbl, bl);
      maxw = MAX (fabs(*w), maxw);
      
      /* update data pointers */
      u += inUV->myDesc->lrec;
      v += inUV->myDesc->lrec;
      w += inUV->myDesc->lrec;
    } /* end loop over buffer */
  } /* end loop over file */
  
    /* Close */
  retCode = ObitUVClose (inUV, err);
  if (err->error) Obit_traceback_msg (err, routine, inUV->name);

  /* Save values in the descriptor */
  inUV->myDesc->maxBL = maxbl;
  inUV->myDesc->maxW  = maxw;

  *MaxBL = maxbl;
  *MaxW  = maxw;
  
} /* end ObitUVUtilUVWExtrema */

/**
 * Make copy of ObitUV with the visibilities zeroed and weights 1.
 * \param inUV     Input uv data to copy. 
 * \param scratch  True if scratch file desired, will be same type as inUV.
 * \param outUV    If not scratch, then the previously defined output file
 *                 May be NULL for scratch only
 *                 If it exists and scratch, it will be Unrefed
 * \param err      Error stack, returns if not empty.
 * \return the zeroed ObitUV.
 */
ObitUV* ObitUVUtilCopyZero (ObitUV *inUV, gboolean scratch, ObitUV *outUV, 
			    ObitErr *err)
{
  ObitIOCode iretCode, oretCode;
  gboolean doCalSelect;
  gchar *exclude[]={"AIPS CL","AIPS SN","AIPS FG","AIPS CQ","AIPS WX",
		    "AIPS AT","AIPS CT","AIPS OB","AIPS IM","AIPS MC",
		    "AIPS PC","AIPS NX","AIPS TY","AIPS GC","AIPS HI",
		    "AIPS PL","AIPS NI","AIPS SY","AIPS PT","AIPS OT",
		    NULL};
  gchar *sourceInclude[] = {"AIPS SU", NULL};
  olong i, j, indx;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  ObitIOAccess access;
  ObitUVDesc *inDesc, *outDesc;
  gchar *today=NULL;
  gchar *routine = "ObitUVUtilCopyZero";
 
  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return outUV;
  g_assert (ObitUVIsA(inUV));
  if (!scratch && (outUV==NULL)) {
    Obit_log_error(err, OBIT_Error,"%s Output MUST be defined for non scratch files",
		   routine);
      return outUV;
  }

  /* Create scratch? */
  if (scratch) {
    if (outUV) outUV = ObitUVUnref(outUV);
    outUV = newObitUVScratch (inUV, err);
  } else { /* non scratch output must exist - clone from inUV */
   outUV->myDesc = ObitUVDescCopy (inUV->myDesc, outUV->myDesc, err);
    ObitUVClone (inUV, outUV, err);
  }
  if (err->error) Obit_traceback_val (err, routine, inUV->name, inUV);

  /* Selection of input? */
  doCalSelect = FALSE;
  ObitInfoListGetTest(inUV->info, "doCalSelect", &type, (gint32*)dim, &doCalSelect);
  if (doCalSelect) access = OBIT_IO_ReadCal;
  else access = OBIT_IO_ReadOnly;

  /* test open to fully instantiate input and see if it's OK */
  iretCode = ObitUVOpen (inUV, access, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine, inUV->name, outUV);

  /* copy Descriptor */
  outUV->myDesc = ObitUVDescCopy(inUV->myDesc, outUV->myDesc, err);

  /* Creation date today */
  today = ObitToday();
  strncpy (outUV->myDesc->date, today, UVLEN_VALUE-1);
  if (today) g_free(today);
  
  /* use same data buffer on input and output 
     so don't assign buffer for output */
  if (outUV->buffer) ObitIOFreeBuffer(outUV->buffer); /* free existing */
  outUV->buffer = NULL;
  outUV->bufferSize = -1;

  /* test open output */
  oretCode = ObitUVOpen (outUV, OBIT_IO_WriteOnly, err);
  /* If this didn't work try OBIT_IO_ReadWrite */
  if ((oretCode!=OBIT_IO_OK) || (err->error)) {
    ObitErrClear(err);
    oretCode = ObitUVOpen (outUV, OBIT_IO_ReadWrite, err);
  }
  /* if it didn't work bail out */
  if ((oretCode!=OBIT_IO_OK) || (err->error)) {
    /* unset output buffer (may be multiply deallocated) */
    outUV->buffer = NULL;
    outUV->bufferSize = 0;
    Obit_traceback_val (err, routine, outUV->name, outUV);
  }

  /* iretCode = ObitUVClose (inUV, err); DEBUG */
  /* Copy tables before data */
  iretCode = ObitUVCopyTables (inUV, outUV, exclude, NULL, err);
  /* If multisource out then copy SU table, multiple sources selected or
   sources deselected suggest MS out */
  if ((inUV->mySel->numberSourcesList>1) || (!inUV->mySel->selectSources))
  iretCode = ObitUVCopyTables (inUV, outUV, NULL, sourceInclude, err);
  if (err->error) {
    outUV->buffer = NULL;
    outUV->bufferSize = 0;
    Obit_traceback_val (err, routine, inUV->name, outUV);
  }

  /* reset to beginning of uv data */
  iretCode = ObitIOSet (inUV->myIO,  inUV->info, err);
  oretCode = ObitIOSet (outUV->myIO, outUV->info, err);
  if (err->error) Obit_traceback_val (err, routine,inUV->name, outUV);

  /* Close and reopen input to init calibration which will have been disturbed 
     by the table copy */
  iretCode = ObitUVClose (inUV, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine,inUV->name, outUV);

  iretCode = ObitUVOpen (inUV, access, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine,inUV->name, outUV);
  outUV->buffer = inUV->buffer;

  /* Get descriptors */
  inDesc  = inUV->myDesc;
  outDesc = outUV->myDesc;

  /* we're in business, copy, zero data, set weight to 1 */
  while ((iretCode==OBIT_IO_OK) && (oretCode==OBIT_IO_OK)) {
    if (doCalSelect) iretCode = ObitUVReadSelect (inUV, inUV->buffer, err);
    else iretCode = ObitUVRead (inUV, inUV->buffer, err);
    if (iretCode!=OBIT_IO_OK) break;
   /* How many */
    outDesc->numVisBuff = inDesc->numVisBuff;

    /* Modify data */
    for (i=0; i<inDesc->numVisBuff; i++) { /* loop over visibilities */
      indx = i*inDesc->lrec + inDesc->nrparm;
      for (j=0; j<inDesc->ncorr; j++) { /* loop over correlations */
	inUV->buffer[indx]   = 0.0;
	inUV->buffer[indx+1] = 0.0;
	inUV->buffer[indx+2] = 1.0;
	indx += inDesc->inaxes[0];
      } /* end loop over correlations */
    } /* end loop over visibilities */

    /* Write */
    oretCode = ObitUVWrite (outUV, inUV->buffer, err);
  } /* end loop processing data */
  
  /* check for errors */
  if ((iretCode > OBIT_IO_EOF) || (oretCode > OBIT_IO_EOF) ||
      (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine,inUV->name, outUV);
  
  /* unset input buffer (may be multiply deallocated ;'{ ) */
  outUV->buffer = NULL;
  outUV->bufferSize = 0;
  
  /* close files */
  iretCode = ObitUVClose (inUV, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) 
    Obit_traceback_val (err, routine, inUV->name, outUV);
  
  oretCode = ObitUVClose (outUV, err);
  if ((oretCode!=OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, outUV->name, outUV);
  
  return outUV;
} /* end ObitUVUtilCopyZero */

#if HAVE_AVX512==1
#include <immintrin.h>
# define ALIGN32_BEG
# define ALIGN32_END __attribute__((packed,aligned(32)))
/* Union allowing c interface */
typedef __m512  V16SF; // vector of 16 float (avx)
typedef __m512i V16SI; // vector of 16 int   (avx)
typedef __mmask16 MASK16; // vector of 16 mask  (avx)
typedef ALIGN32_BEG union {
  float f[16];
  int   i[16];
  V16SF   v;
} ALIGN32_END CV16SF;
typedef ALIGN32_BEG union {
  float f[16];
  int   i[16];
  V16SI   v;
} ALIGN32_END IV16SF;
#endif
/**
 * Divide the visibilities in one ObitUV by those in another
 * outUV = inUV1 / inUV2
 * \param inUV1    Input uv data numerator, no calibration/selection
 * \param inUV2    Input uv data denominator, no calibration/selection
 *                 inUV2 should have the same structure, no. vis etc
 *                 as inUV1.
 * \param outUV    Previously defined output, may be the same as inUV1
 * \param err      Error stack, returns if not empty.
 */
void ObitUVUtilVisDivide (ObitUV *inUV1, ObitUV *inUV2, ObitUV *outUV, 
			  ObitErr *err)
{
  ObitIOCode iretCode, oretCode;
  gchar *exclude[]={"AIPS CL","AIPS SN","AIPS FG","AIPS CQ","AIPS WX",
		    "AIPS AT","AIPS CT","AIPS OB","AIPS IM","AIPS MC",
		    "AIPS PC","AIPS NX","AIPS TY","AIPS GC","AIPS HI",
		    "AIPS PL","AIPS NI","AIPS SY","AIPS PT","AIPS OT",
		    NULL};
  gchar *sourceInclude[] = {"AIPS SU", NULL};
  olong i, j, indx, firstVis;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong NPIO, jlast=0;
  ofloat work[3];
  gboolean incompatible, same, btemp;
  ObitUVDesc *in1Desc, *in2Desc, *outDesc;
  gchar *today=NULL;
  gchar *routine = "ObitUVUtilVisDivide";

#if HAVE_AVX512==1
  CV16SF vr1, vi1, vw1, vr2, vi2, vw2, vz, v1, vwk1, vwk2, vwk3, vt;
  IV16SF vindxr, vindxi, vindxw;
  MASK16 msk, msk1, msk2;
  olong  jndx;
  vz.v     = _mm512_set1_ps(0.0);  /* vector of zeroes */
  v1.v     = _mm512_set1_ps(1.0);  /* vector of ones */
  vindxr.v = _mm512_set_epi32 (45,42,39,36,33,30,27,24,21,18,15,12, 9,6,3,0); /* Reals*/
  vindxi.v = _mm512_set_epi32 (46,43,40,37,34,32,28,25,22,19,16,13,10,7,4,1); /* Imags */
  vindxw.v = _mm512_set_epi32 (47,44,41,38,35,32,29,26,23,20,17,14,11,8,5,2); /* Weight */
#endif

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitUVIsA(inUV1));
  g_assert (ObitUVIsA(inUV2));
  g_assert (ObitUVIsA(outUV));

  /* Are input1 and output the same file? */
  same = ObitUVSame(inUV1, outUV, err);
  if (err->error) Obit_traceback_msg (err, routine, inUV1->name);

  /* test open to fully instantiate input and see if it's OK */
  iretCode = ObitUVOpen (inUV1, OBIT_IO_ReadWrite, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_msg (err, routine,inUV1->name);

  /* copy Descriptor */
  outUV->myDesc = ObitUVDescCopy(inUV1->myDesc, outUV->myDesc, err);

  /* Creation date today */
  today = ObitToday();
  strncpy (outUV->myDesc->date, today, UVLEN_VALUE-1);
  if (today) g_free(today);
 
  /* use same data buffer on input and output.
     If multiple passes are made the input files will be closed
     which deallocates the buffer, use output buffer.
     so free input buffer */
  if (!same) {
    /* use same data buffer on input 1 and output 
       so don't assign buffer for output */
    if (outUV->buffer) ObitIOFreeBuffer(outUV->buffer); /* free existing */
    outUV->buffer = NULL;
    outUV->bufferSize = -1;
  }

   /* Copy number of records per IO to second input */
  ObitInfoListGet (inUV1->info, "nVisPIO", &type, dim,  (gpointer)&NPIO, err);
  ObitInfoListPut (inUV2->info, "nVisPIO",  type, dim,  (gpointer)&NPIO, err);
  if (err->error) Obit_traceback_msg (err, routine, inUV1->name);

  /* Open second input */
  iretCode = ObitUVOpen (inUV2, OBIT_IO_ReadWrite, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) 
    Obit_traceback_msg (err, routine,inUV2->name);

  /* Get input descriptors */
  in1Desc = inUV1->myDesc;
  in2Desc = inUV2->myDesc;

  /* Check compatability between inUV1, inUV2 */
  incompatible = in1Desc->nvis!=in2Desc->nvis;
  incompatible = incompatible || (in1Desc->ncorr!=in2Desc->ncorr);
  incompatible = incompatible || (in1Desc->jlocs!=in2Desc->jlocs);
  incompatible = incompatible || (in1Desc->jlocf!=in2Desc->jlocf);
  incompatible = incompatible || (in1Desc->jlocif!=in2Desc->jlocif);
  incompatible = incompatible || (in1Desc->ilocb!=in2Desc->ilocb);
  if (incompatible) {
     Obit_log_error(err, OBIT_Error,"%s inUV1 and inUV2 have incompatible structures",
		   routine);
      return ;
 }

  /* Look at all data */
  btemp = TRUE;
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut(inUV1->info, "passAll", OBIT_bool, dim, &btemp);
  ObitInfoListAlwaysPut(inUV2->info, "passAll", OBIT_bool, dim, &btemp);

 /* test open output */
  oretCode = ObitUVOpen (outUV, OBIT_IO_WriteOnly, err);
  /* If this didn't work try OBIT_IO_ReadWrite */
  if ((oretCode!=OBIT_IO_OK) || (err->error)) {
    ObitErrClear(err);
    oretCode = ObitUVOpen (outUV, OBIT_IO_ReadWrite, err);
  }
  /* if it didn't work bail out */
  if ((oretCode!=OBIT_IO_OK) || (err->error)) {
    /* unset output buffer (may be multiply deallocated) */
    outUV->buffer = NULL;
    outUV->bufferSize = 0;
    Obit_traceback_msg (err, routine, outUV->name);
  }

  /* Copy tables before data if in1 and out are not the same */
  if (!ObitUVSame (inUV1, outUV, err)) {
    iretCode = ObitUVCopyTables (inUV1, outUV, exclude, NULL, err);
    /* If multisource out then copy SU table, multiple sources selected or
       sources deselected suggest MS out */
    if ((inUV1->mySel->numberSourcesList>1) || (!inUV1->mySel->selectSources))
      iretCode = ObitUVCopyTables (inUV1, outUV, NULL, sourceInclude, err);
    if (err->error) {
      outUV->buffer = NULL;
      outUV->bufferSize = 0;
      Obit_traceback_msg (err, routine, inUV1->name);
    }
  }
  if (err->error) Obit_traceback_msg (err, routine, inUV1->name);

  /* reset to beginning of uv data */
  iretCode = ObitUVIOSet (inUV1, err);
  oretCode = ObitUVIOSet (outUV, err);
  if (err->error) Obit_traceback_msg (err, routine,inUV1->name);

  /* Close and reopen input to init calibration which will have been disturbed 
     by the table copy */
  iretCode = ObitUVClose (inUV1, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error))
    Obit_traceback_msg (err, routine,inUV1->name);
  iretCode = ObitUVOpen (inUV1, OBIT_IO_ReadWrite, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error))
    Obit_traceback_msg (err, routine,inUV1->name);
  outUV->buffer = inUV1->buffer;

  outDesc = outUV->myDesc;   /* Get output descriptor */

  /* we're in business, divide */
  while ((iretCode==OBIT_IO_OK) && (oretCode==OBIT_IO_OK)) {
    /* Read first input */
    iretCode = ObitUVRead (inUV1, inUV1->buffer, err);
    if (iretCode!=OBIT_IO_OK) break;
    /* Read second input */
    iretCode = ObitUVRead (inUV2, inUV2->buffer, err);
    if (iretCode!=OBIT_IO_OK) break;
    /* How many */
    outDesc->numVisBuff = in1Desc->numVisBuff;
    firstVis = in1Desc->firstVis;

    /* compatability check */
    incompatible = in1Desc->numVisBuff!=in2Desc->numVisBuff;
    if (incompatible) break;

    /* Modify data */
    for (i=0; i<in1Desc->numVisBuff; i++) { /* loop over visibilities */
      /* compatability check - check time and baseline or antenna code */
      indx = i*in1Desc->lrec;
      incompatible = 
	inUV1->buffer[indx+in1Desc->iloct] !=inUV2->buffer[indx+in2Desc->iloct];
      if (in1Desc->ilocb>=0) {
	incompatible = incompatible ||
	  inUV1->buffer[indx+in1Desc->ilocb] !=inUV2->buffer[indx+in2Desc->ilocb];
      }
      if (in1Desc->iloca1>=0) {
	incompatible = incompatible ||
	  inUV1->buffer[indx+in1Desc->iloca1]!=inUV2->buffer[indx+in2Desc->iloca1] ||
	  inUV1->buffer[indx+in1Desc->iloca2]!=inUV2->buffer[indx+in2Desc->iloca2];
      }
      if (incompatible) break;

#if HAVE_AVX512==1  /* Vector */
      for (j=0; j<in1Desc->ncorr-16; j+=16) { /* loop over correlations */
        jndx   = i*in1Desc->lrec + in1Desc->nrparm + in1Desc->inaxes[0]*j;
	vr1.v  = _mm512_i32gather_ps (vindxr.v, &inUV1->buffer[jndx], 4);/* Load Reals 1 */
	vi1.v  = _mm512_i32gather_ps (vindxi.v, &inUV1->buffer[jndx], 4);/* Load Imaginaries 1 */
	vw1.v  = _mm512_i32gather_ps (vindxw.v, &inUV1->buffer[jndx], 4);/* Load Weights 1 */
	vr2.v  = _mm512_i32gather_ps (vindxr.v, &inUV2->buffer[jndx], 4);/* Load Reals 2 */
	vi2.v  = _mm512_i32gather_ps (vindxi.v, &inUV2->buffer[jndx], 4);/* Load Imaginaries 2 */
	vw2.v  = _mm512_i32gather_ps (vindxw.v, &inUV2->buffer[jndx], 4);/* Load Weights 2 */
	msk1   = _mm512_cmp_ps_mask(vw1.v, vz.v, _CMP_LE_OQ);            /* find invalid data 1 */
	msk2   = _mm512_cmp_ps_mask(vw2.v, vz.v, _CMP_LE_OQ);            /* find invalid data 2 */
	msk    = _mm512_kxnor (msk1, msk2);                              /* Combined valid */
	/* Complex Divide */
	vwk2.v = _mm512_mul_ps(vr2.v, vr2.v);
	vwk3.v = _mm512_mul_ps(vi2.v, vi2.v);
	vwk1.v = _mm512_add_ps(vwk2.v, vwk3.v);               /* r2^2 + i2^2 */
	vt.v   = _mm512_sqrt_ps(vwk1.v);
	vw1.v  = _mm512_mul_ps(vw1.v, vt.v);                  /* Weight */
	msk1   = _mm512_cmp_ps_mask(vwk1.v, vz.v, _CMP_NEQ_UQ);/* trap zero */
	vwk1.v = _mm512_mask_rcp14_ps(v1.v, msk1, vwk1.v);    /* approximate reciprocal */ 
	vwk2.v = _mm512_mul_ps(vr1.v, vwk1.v);
	vwk3.v = _mm512_mul_ps(vi1.v, vwk1.v);
	vt.v   = _mm512_mul_ps(vwk2.v, vr2.v);
	vr1.v  = _mm512_mul_ps(vwk3.v, vi2.v);
	vr1.v  = _mm512_add_ps(vt.v, vr1.v);
	vt.v   = _mm512_mul_ps(vwk3.v, vr2.v);
	vi1.v  = _mm512_mul_ps(vwk2.v, vi2.v);
	vi1.v  = _mm512_sub_ps(vt.v, vi1.v);
	/* Weight */
	/* Save output */
	_mm512_mask_i32scatter_ps (&inUV1->buffer[jndx], msk, vindxr.v, vr1.v, 4); /* save valid real */
	_mm512_mask_i32scatter_ps (&inUV1->buffer[jndx], msk, vindxi.v, vi1.v, 4); /* save valid imaginary */
	_mm512_mask_i32scatter_ps (&inUV1->buffer[jndx], msk, vindxw.v, vw1.v, 4); /* save valid weight */
      } /* end loop over correlations */
      jlast = j;  /* How far did it get */
#else /* Scalar */
      jlast = 0;  /* Do all */
#endif
      indx = i*in1Desc->lrec + in1Desc->nrparm + jlast*in1Desc->inaxes[0];
      for (j=jlast; j<in1Desc->ncorr; j++) { /* loop over correlations */
	/* Divide */
	ObitUVWtCpxDivide ((&inUV1->buffer[indx]), 
			   (&inUV2->buffer[indx]), 
			   (&inUV1->buffer[indx]), work);
	indx += in1Desc->inaxes[0];
      } /* end loop over correlations */
      if (incompatible) break;
    } /* end loop over visibilities */
    
    /* Write */
    oretCode = ObitUVWrite (outUV, inUV1->buffer, err);
    if (same) {
      outUV->myDesc->firstVis = firstVis;
      ((ObitUVDesc*)(outUV->myIO->myDesc))->firstVis = firstVis;
    }
  } /* end loop processing data */
  
  /* Check for incompatibility */
  if (incompatible) {
    Obit_log_error(err, OBIT_Error,"%s inUV1 and inUV2 have incompatible contents",
		   routine);
    return;
  }
    
  /* check for errors */
  if ((iretCode > OBIT_IO_EOF) || (oretCode > OBIT_IO_EOF) ||
      (err->error)) /* add traceback,return */
    Obit_traceback_msg (err, routine,inUV1->name);
    
  /* unset output buffer (may be multiply deallocated ;'{ ) */
  outUV->buffer = NULL;
  outUV->bufferSize = 0;
  
  /* close files */
  iretCode = ObitUVClose (inUV1, err);
  if (err->error) Obit_traceback_msg (err, routine, inUV1->name);
  
  iretCode = ObitUVClose (inUV2, err);
  if (err->error) Obit_traceback_msg (err, routine, inUV2->name);
  
  oretCode = ObitUVClose (outUV, err);
  if (err->error) Obit_traceback_msg (err, routine, outUV->name);
  /* In case */
  outUV->buffer     = NULL;
  outUV->bufferSize = 0;
 
  /* Reset passAll */
  btemp = FALSE;
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut(inUV1->info, "passAll", OBIT_bool, dim, &btemp);
  ObitInfoListAlwaysPut(inUV2->info, "passAll", OBIT_bool, dim, &btemp);

} /* end ObitUVUtilVisDivide */

/**
 * Divide the cross pol visibilities in one ObitUV by I pol (avg parallel pol)
 * \param inUV     Input uv data, if info member KeepSou TRUE, copy SU table.
 * \param outUV    Previously defined output, may be the same as inUV
 * \param err      Error stack, returns if not empty.
 */
void ObitUVUtilXPolDivide (ObitUV *inUV, ObitUV *outUV, ObitErr *err)
{
  ObitIOCode iretCode, oretCode;
  gchar *exclude[]={"AIPS CL","AIPS SN","AIPS FG","AIPS CQ","AIPS WX",
		    "AIPS AT","AIPS CT","AIPS OB","AIPS IM","AIPS MC",
		    "AIPS PC","AIPS NX","AIPS TY","AIPS GC","AIPS HI",
		    "AIPS PL","AIPS NI","AIPS SY","AIPS PT","AIPS OT",
		    NULL};
  gchar *sourceInclude[] = {"AIPS SU", NULL};
  olong i, indx, jndx, firstVis;
  olong iif, nif, ichan, nchan, istok, nstok, incs, incf, incif;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ofloat work[3], Ivis[3], wt, Ireal, Iimag, Iwt;
  gboolean KeepSou, same;
  ObitUVDesc *inDesc, *outDesc;
  gchar *today=NULL;
  gchar *routine = "ObitUVUtilXPolDivide";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitUVIsA(inUV));
  g_assert (ObitUVIsA(outUV));

  /* Are input1 and output the same file? */
  same = ObitUVSame(inUV, outUV, err);
  if (err->error) Obit_traceback_msg (err, routine, inUV->name);

  /* Keep source table? */
  KeepSou = FALSE;
  ObitInfoListGetTest(inUV->info, "KeepSou", &type, (gint32*)dim, &KeepSou);

  /* test open to fully instantiate input and see if it's OK */
  iretCode = ObitUVOpen (inUV, OBIT_IO_ReadWrite, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_msg (err, routine, inUV->name);

  /* Get input descriptor */
  inDesc = inUV->myDesc;

   /* Set up for parsing data */
  nchan = inDesc->inaxes[inDesc->jlocf];
  if (inDesc->jlocif>=0) nif = inDesc->inaxes[inDesc->jlocif];
  else nif = 1;
  if (inDesc->jlocs>=0) nstok = inDesc->inaxes[inDesc->jlocs];
  else nstok = 1;
  /* get increments */
  incs  = inDesc->incs;
  incf  = inDesc->incf;
  incif = inDesc->incif;

  /* Make sure at least 4 Stokes correlations */
  Obit_return_if_fail ((nstok>=4), err,
		       "%s: MUST have at least 4 Stokes, have  %d",  
		       routine, nstok);  

  /* copy Descriptor */
  outUV->myDesc = ObitUVDescCopy(inUV->myDesc, outUV->myDesc, err);

  /* Creation date today */
  today = ObitToday();
  strncpy (outUV->myDesc->date, today, UVLEN_VALUE-1);
  if (today) g_free(today);
 
  /* use same data buffer on input and output.
     If multiple passes are made the input files will be closed
     which deallocates the buffer, use output buffer.
     so free input buffer */
  if (!same) {
    /* use same data buffer on input 1 and output 
       so don't assign buffer for output */
    if (outUV->buffer) ObitIOFreeBuffer(outUV->buffer); /* free existing */
    outUV->buffer = NULL;
    outUV->bufferSize = -1;
  }

  /* test open output */
  oretCode = ObitUVOpen (outUV, OBIT_IO_WriteOnly, err);
  /* If this didn't work try OBIT_IO_ReadWrite */
  if ((oretCode!=OBIT_IO_OK) || (err->error)) {
    ObitErrClear(err);
    oretCode = ObitUVOpen (outUV, OBIT_IO_ReadWrite, err);
  }
  /* if it didn't work bail out */
  if ((oretCode!=OBIT_IO_OK) || (err->error)) {
    /* unset output buffer (may be multiply deallocated) */
    outUV->buffer = NULL;
    outUV->bufferSize = 0;
    Obit_traceback_msg (err, routine, outUV->name);
  }

  /* Copy tables before data if in1 and out are not the same */
  if (!ObitUVSame (inUV, outUV, err)) {
    iretCode = ObitUVCopyTables (inUV, outUV, exclude, NULL, err);
    /* If multisource out then copy SU table, multiple sources selected or
       sources deselected suggest MS out or KeepSou == TRUE */
    if ((inUV->mySel->numberSourcesList>1) || (!inUV->mySel->selectSources) ||
	KeepSou)
      iretCode = ObitUVCopyTables (inUV, outUV, NULL, sourceInclude, err);
    if (err->error) {
      outUV->buffer = NULL;
      outUV->bufferSize = 0;
      Obit_traceback_msg (err, routine, inUV->name);
    }
  }
  if (err->error) Obit_traceback_msg (err, routine, inUV->name);

  /* reset to beginning of uv data */
  iretCode = ObitUVIOSet (inUV, err);
  oretCode = ObitUVIOSet (outUV, err);
  if (err->error) Obit_traceback_msg (err, routine,inUV->name);

  /* Close and reopen input to init calibration which will have been disturbed 
     by the table copy */
  iretCode = ObitUVClose (inUV, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error))
    Obit_traceback_msg (err, routine,inUV->name);
  iretCode = ObitUVOpen (inUV, OBIT_IO_ReadWrite, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error))
    Obit_traceback_msg (err, routine,inUV->name);
  outUV->buffer = inUV->buffer;

  outDesc = outUV->myDesc;   /* Get output descriptor */

  /* we're in business, divide */
  while ((iretCode==OBIT_IO_OK) && (oretCode==OBIT_IO_OK)) {
    /* Read input */
    iretCode = ObitUVRead (inUV, inUV->buffer, err);
    if (iretCode!=OBIT_IO_OK) break;
    /* How many */
    outDesc->numVisBuff = inDesc->numVisBuff;
    firstVis = inDesc->firstVis;

    /* Modify data */
    for (i=0; i<inDesc->numVisBuff; i++) { /* loop over visibilities */
      indx = inDesc->nrparm + i*inDesc->lrec;
      /* loop over IF */
      for (iif=0; iif<nif; iif++) {
	/* Loop over frequency channel */
	for (ichan=0; ichan<nchan; ichan++) { /* loop 60 */
	  jndx = indx + iif*incif + ichan*incf;
	  /* Get Stokes I - average parallel hands */
	  wt = inUV->buffer[jndx+2];
	  if (wt>0.0) {
	    Ireal = wt*inUV->buffer[jndx];
	    Iimag = wt*inUV->buffer[jndx+1];
	    Iwt   = wt;
	  } else {
	    Ireal = Iimag = Iwt = 0.0;
	  }
	  wt = inUV->buffer[jndx+5];
	  if (wt>0.0) {
	    Ireal += wt*inUV->buffer[jndx+3];
	    Iimag += wt*inUV->buffer[jndx+4];
	    Iwt   += wt;
	  }
	  if (Iwt > 0.0) {
	    Ivis[0] = Ireal/Iwt;
	    Ivis[1] = Iimag/Iwt;
	    Ivis[2] = Iwt;
	  } else {
	    Ivis[0] = Ivis[1] = Ivis[2] = 0.0;
	  }
	  /* Loop over cross polarization */
	  for (istok=2; istok<nstok; istok++) {
	    jndx = indx + iif*incif + ichan*incf + istok*incs;
	    /* Divide */
	    ObitUVWtCpxDivide ((&inUV->buffer[jndx]), Ivis,
			       (&inUV->buffer[jndx]), work);
	  } /* end loop over Stokes */
	} /* end loop over channels */
      } /* end loop over IFs */
    } /* end loop over visibilities */
    
    /* Write */
    oretCode = ObitUVWrite (outUV, inUV->buffer, err);
    if (same) {
      outUV->myDesc->firstVis = firstVis;
      ((ObitUVDesc*)(outUV->myIO->myDesc))->firstVis = firstVis;
    }
  } /* end loop processing data */
  
  /* check for errors */
  if ((iretCode > OBIT_IO_EOF) || (oretCode > OBIT_IO_EOF) ||
      (err->error)) /* add traceback,return */
    Obit_traceback_msg (err, routine,inUV->name);
    
  /* unset output buffer (may be multiply deallocated ;'{ ) */
  outUV->buffer = NULL;
  outUV->bufferSize = 0;
  
  /* close files */
  iretCode = ObitUVClose (inUV, err);
  if (err->error) Obit_traceback_msg (err, routine, inUV->name);
  
  oretCode = ObitUVClose (outUV, err);
  if (err->error) Obit_traceback_msg (err, routine, outUV->name);
  /* In case */
  outUV->buffer     = NULL;
  outUV->bufferSize = 0;
 
} /* end ObitUVUtilXPolDivide */

/**
 * Subtract the visibilities in one ObitUV from those in another
 * outUV = inUV1 - inUV2
 * \param inUV1    First input uv data, no calibration/selection
 * \param inUV2    Second input uv data, calibration/selection allowed
 *                 inUV2 should have the same structure, no. vis etc
 *                 as inUV1.
 *                 On info member (in addition to calibration info):
 * \li "Factor"    OBIT_float (1,1,1) Vis scaling factor (def 1.0)
 * \param outUV    Previously defined output, may be the same as inUV1
 * \param err      Error stack, returns on error
 */
void ObitUVUtilVisSub (ObitUV *inUV1, ObitUV *inUV2, ObitUV *outUV, 
		       ObitErr *err)
{
  ObitIOCode iretCode, oretCode;
  gchar *exclude[]={"AIPS CL","AIPS SN","AIPS FG","AIPS CQ","AIPS WX",
		    "AIPS AT","AIPS CT","AIPS OB","AIPS IM","AIPS MC",
		    "AIPS PC","AIPS NX","AIPS TY","AIPS GC","AIPS HI",
		    "AIPS PL","AIPS NI","AIPS SY","AIPS PT","AIPS OT",
		    NULL};
  gchar *sourceInclude[] = {"AIPS SU", NULL};
  olong i, j, indx, firstVis;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong NPIO;
  ofloat factor=1.0;
  gboolean incompatible, same, doCalSelect, btemp;
  ObitUVDesc *in1Desc, *in2Desc, *outDesc;
  ObitIOAccess access;
  gchar *today=NULL;
  gchar *routine = "ObitUVUtilVisSub";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitUVIsA(inUV1));
  g_assert (ObitUVIsA(inUV2));
  g_assert (ObitUVIsA(outUV));

   ObitInfoListGetTest (inUV2->info, "Factor",  &type, dim, &factor);
 /* Are input1 and output the same file? */
  same = ObitUVSame(inUV1, outUV, err);
  if (err->error) Obit_traceback_msg (err, routine, inUV1->name);

  /* test open to fully instantiate input and see if it's OK */
  iretCode = ObitUVOpen (inUV1, OBIT_IO_ReadWrite, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_msg (err, routine,inUV1->name);

  /* copy Descriptor */
  outUV->myDesc = ObitUVDescCopy(inUV1->myDesc, outUV->myDesc, err);

  /* Creation date today */
  today = ObitToday();
  strncpy (outUV->myDesc->date, today, UVLEN_VALUE-1);
  if (today) g_free(today);
 
  /* use same data buffer on input and output.
     If multiple passes are made the input files will be closed
     which deallocates the buffer, use output buffer.
     so free input buffer */
  if (!same) {
    /* use same data buffer on input 1 and output 
       so don't assign buffer for output */
    if (outUV->buffer) ObitIOFreeBuffer(outUV->buffer); /* free existing */
    outUV->buffer = NULL;
    outUV->bufferSize = -1;
  }

   /* Copy number of records per IO to second input */
  ObitInfoListGet (inUV1->info, "nVisPIO", &type, dim,  (gpointer)&NPIO, err);
  ObitInfoListPut (inUV2->info, "nVisPIO",  type, dim,  (gpointer)&NPIO, err);
  if (err->error) Obit_traceback_msg (err, routine, inUV1->name);

  /* Selection of second input? */
  doCalSelect = FALSE;
  ObitInfoListGetTest(inUV2->info, "doCalSelect", &type, dim, &doCalSelect);
  if (doCalSelect) access = OBIT_IO_ReadCal;
  else access = OBIT_IO_ReadOnly;

  /* Open second input */
  iretCode = ObitUVOpen (inUV2, access, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) 
    Obit_traceback_msg (err, routine,inUV2->name);

  /* Get input descriptors */
  in1Desc = inUV1->myDesc;
  in2Desc = inUV2->myDesc;

  /* Check compatability between inUV1, inUV2 */
  incompatible = in1Desc->nvis!=in2Desc->nvis;
  incompatible = incompatible || (in1Desc->ncorr!=in2Desc->ncorr);
  incompatible = incompatible || (in1Desc->jlocs!=in2Desc->jlocs);
  incompatible = incompatible || (in1Desc->jlocf!=in2Desc->jlocf);
  incompatible = incompatible || (in1Desc->jlocif!=in2Desc->jlocif);
  incompatible = incompatible || (in1Desc->ilocb!=in2Desc->ilocb);
  if (incompatible) {
    if (in1Desc->nvis!=in2Desc->nvis)
      	    Obit_log_error(err, OBIT_InfoWarn, 
			   "nvis %d != %d", in1Desc->nvis, in2Desc->nvis);
    if (in1Desc->ncorr!=in2Desc->ncorr)
      	    Obit_log_error(err, OBIT_InfoWarn, 
			   "ncorr %d != %d", in1Desc->ncorr, in2Desc->ncorr);
    if (in1Desc->jlocs!=in2Desc->jlocs) Obit_log_error(err, OBIT_InfoWarn, "Stokes axes differ");
    if (in1Desc->jlocf!=in2Desc->jlocf) Obit_log_error(err, OBIT_InfoWarn, "Frequency axes differ");
    if (in1Desc->jlocif!=in2Desc->jlocif) Obit_log_error(err, OBIT_InfoWarn, "IF axes differ");
     Obit_log_error(err, OBIT_Error,"%s: %s and %s have incompatible structures",
		    routine, inUV1->name, inUV2->name);
      return ;
  }

  /* Look at all data */
  btemp = TRUE;
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut(inUV1->info, "passAll", OBIT_bool, dim, &btemp);
  ObitInfoListAlwaysPut(inUV2->info, "passAll", OBIT_bool, dim, &btemp);
  
  /* test open output */
  oretCode = ObitUVOpen (outUV, OBIT_IO_WriteOnly, err);
  /* If this didn't work try OBIT_IO_ReadWrite */
  if ((oretCode!=OBIT_IO_OK) || (err->error)) {
    ObitErrClear(err);
    oretCode = ObitUVOpen (outUV, OBIT_IO_ReadWrite, err);
  }
  /* if it didn't work bail out */
  if ((oretCode!=OBIT_IO_OK) || (err->error)) {
    /* unset output buffer (may be multiply deallocated) */
    outUV->buffer = NULL;
    outUV->bufferSize = 0;
    Obit_traceback_msg (err, routine, outUV->name);
  }

  /* Copy tables before data if in1 and out are not the same */
  if (!ObitUVSame (inUV1, outUV, err)) {
    iretCode = ObitUVCopyTables (inUV1, outUV, exclude, NULL, err);
    /* If multisource out then copy SU table, multiple sources selected or
       sources deselected suggest MS out */
    if ((inUV1->mySel->numberSourcesList>1) || (!inUV1->mySel->selectSources))
      iretCode = ObitUVCopyTables (inUV1, outUV, NULL, sourceInclude, err);
    if (err->error) {
      outUV->buffer = NULL;
      outUV->bufferSize = 0;
      Obit_traceback_msg (err, routine, inUV1->name);
    }
  }
  if (err->error) Obit_traceback_msg (err, routine, inUV1->name);

  /* reset to beginning of uv data */
  iretCode = ObitUVIOSet (inUV1, err);
  oretCode = ObitUVIOSet (outUV, err);
  if (err->error) Obit_traceback_msg (err, routine,inUV1->name);

  /* Close and reopen input to init calibration which will have been disturbed 
     by the table copy */
  iretCode = ObitUVClose (inUV1, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error))
    Obit_traceback_msg (err, routine,inUV1->name);
  iretCode = ObitUVOpen (inUV1, OBIT_IO_ReadWrite, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error))
    Obit_traceback_msg (err, routine,inUV1->name);
  outUV->buffer = inUV1->buffer;

  outDesc = outUV->myDesc;   /* Get output descriptor */

  /* we're in business, subtract */
  while ((iretCode==OBIT_IO_OK) && (oretCode==OBIT_IO_OK)) {
    /* Read first input */
    iretCode = ObitUVRead (inUV1, inUV1->buffer, err);
    if (iretCode!=OBIT_IO_OK) break;
    /* Read second input */
    if (doCalSelect) iretCode = ObitUVReadSelect (inUV2, inUV2->buffer, err);
    else iretCode = ObitUVRead (inUV2, inUV2->buffer, err);
    if (iretCode!=OBIT_IO_OK) break;
    /* How many */
    outDesc->numVisBuff = in1Desc->numVisBuff;
    firstVis = in1Desc->firstVis;

    /* compatability check */
    incompatible = in1Desc->numVisBuff!=in2Desc->numVisBuff;
    if (incompatible) break;

    /* Modify data */
    for (i=0; i<in1Desc->numVisBuff; i++) { /* loop over visibilities */
      /* compatability check - check time and baseline or antenna code */
      indx = i*in1Desc->lrec ;
      incompatible = 
	inUV1->buffer[indx+in1Desc->iloct] !=inUV2->buffer[indx+in2Desc->iloct];
      if (in1Desc->ilocb>=0) {
	incompatible = incompatible ||
	  inUV1->buffer[indx+in1Desc->ilocb] !=inUV2->buffer[indx+in2Desc->ilocb];
      }
      if (in1Desc->iloca1>=0) {
	incompatible = incompatible ||
	  inUV1->buffer[indx+in1Desc->iloca1]!=inUV2->buffer[indx+in2Desc->iloca1] ||
	  inUV1->buffer[indx+in1Desc->iloca2]!=inUV2->buffer[indx+in2Desc->iloca2];
      }
      if (incompatible) break;

      indx += in1Desc->nrparm;
      for (j=0; j<in1Desc->ncorr; j++) { /* loop over correlations */
	/* subtract */
	if ((inUV1->buffer[indx+2]>0.0) && (inUV2->buffer[indx+2]>0.0)) {
	  inUV1->buffer[indx]   -= factor * inUV2->buffer[indx];
	  inUV1->buffer[indx+1] -= factor * inUV2->buffer[indx+1];
	  /* this blanks poln data } else {
	    inUV1->buffer[indx]   = 0.0;
	    inUV1->buffer[indx+1] = 0.0;
	    inUV1->buffer[indx+2] = 0.0;*/
	}
	indx += in1Desc->inaxes[0];
      } /* end loop over correlations */
      if (incompatible) break;
    } /* end loop over visibilities */
    
    /* Write */
    oretCode = ObitUVWrite (outUV, inUV1->buffer, err);
    /* suppress vis number update if rewriting the same file */
    if (same) {
      outUV->myDesc->firstVis = firstVis;
      ((ObitUVDesc*)(outUV->myIO->myDesc))->firstVis = firstVis;
    }
 } /* end loop processing data */
  
  /* Check for incompatibility */
  if (incompatible) {
    Obit_log_error(err, OBIT_Error,"%s inUV1 and inUV2 have incompatible contents",
		   routine);
    return;
  }
    
  /* check for errors */
  if ((iretCode > OBIT_IO_EOF) || (oretCode > OBIT_IO_EOF) ||
      (err->error)) /* add traceback,return */
    Obit_traceback_msg (err, routine,inUV1->name);
    
  /* unset output buffer (may be multiply deallocated ;'{ ) */
  outUV->buffer = NULL;
  outUV->bufferSize = 0;
  
  /* Reset passAll */
  btemp = FALSE;
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut(inUV1->info, "passAll", OBIT_bool, dim, &btemp);
  ObitInfoListAlwaysPut(inUV2->info, "passAll", OBIT_bool, dim, &btemp);
  
  /* close files */
  iretCode = ObitUVClose (inUV1, err);
  if (err->error) Obit_traceback_msg (err, routine, inUV1->name);
  
  iretCode = ObitUVClose (inUV2, err);
  if (err->error) Obit_traceback_msg (err, routine, inUV2->name);
  
  oretCode = ObitUVClose (outUV, err);
  if (err->error) Obit_traceback_msg (err, routine, outUV->name);
  
} /* end ObitUVUtilVisSub */

/**
 * Subtract the 1st visibility in one ObitUV from those in another
 * outUV = inUV1[*] - inUV2[0]
 * \param inUV1    First input uv data, no calibration/selection
 * \param inUV2    Second input uv data, calibration/selection allowed
 *                 inUV2 should have the same structure, only the 1st
 *                 visibility used and subtracted from all in inUV1
 * \param outUV    Previously defined output, may be the same as inUV1
 * \param err      Error stack, returns on error
 */
void ObitUVUtilVisSub1 (ObitUV *inUV1, ObitUV *inUV2, ObitUV *outUV, 
		        ObitErr *err)
{
  ObitIOCode iretCode, oretCode;
  gchar *exclude[]={"AIPS CL","AIPS SN","AIPS FG","AIPS CQ","AIPS WX",
		    "AIPS AT","AIPS CT","AIPS OB","AIPS IM","AIPS MC",
		    "AIPS PC","AIPS NX","AIPS TY","AIPS GC","AIPS HI",
		    "AIPS PL","AIPS NI","AIPS SY","AIPS PT","AIPS OT",
		    NULL};
  gchar *sourceInclude[] = {"AIPS SU", NULL};
  olong i, j, indx, kndx, firstVis;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong NPIO;
  gboolean incompatible, same, doCalSelect, btemp;
  ObitUVDesc *in1Desc, *in2Desc, *outDesc;
  ObitIOAccess access;
  gchar *today=NULL;
  gchar *routine = "ObitUVUtilVisSub1";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitUVIsA(inUV1));
  g_assert (ObitUVIsA(inUV2));
  g_assert (ObitUVIsA(outUV));

  /* Are input1 and output the same file? */
  same = ObitUVSame(inUV1, outUV, err);
  if (err->error) Obit_traceback_msg (err, routine, inUV1->name);

  /* test open to fully instantiate input and see if it's OK */
  iretCode = ObitUVOpen (inUV1, OBIT_IO_ReadWrite, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_msg (err, routine,inUV1->name);

  /* copy Descriptor */
  outUV->myDesc = ObitUVDescCopy(inUV1->myDesc, outUV->myDesc, err);

  /* Creation date today */
  today = ObitToday();
  strncpy (outUV->myDesc->date, today, UVLEN_VALUE-1);
  if (today) g_free(today);
 
  /* use same data buffer on input and output.
     If multiple passes are made the input files will be closed
     which deallocates the buffer, use output buffer.
     so free input buffer */
  if (!same) {
    /* use same data buffer on input 1 and output 
       so don't assign buffer for output */
    if (outUV->buffer) ObitIOFreeBuffer(outUV->buffer); /* free existing */
    outUV->buffer = NULL;
    outUV->bufferSize = -1;
  }

   /* Copy number of records per IO to second input */
  ObitInfoListGet (inUV1->info, "nVisPIO", &type, dim,  (gpointer)&NPIO, err);
  ObitInfoListPut (inUV2->info, "nVisPIO",  type, dim,  (gpointer)&NPIO, err);
  if (err->error) Obit_traceback_msg (err, routine, inUV1->name);

  /* Selection of second input? */
  doCalSelect = FALSE;
  ObitInfoListGetTest(inUV2->info, "doCalSelect", &type, dim, &doCalSelect);
  if (doCalSelect) access = OBIT_IO_ReadCal;
  else access = OBIT_IO_ReadOnly;

  /* Open second input */
  iretCode = ObitUVOpen (inUV2, access, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) 
    Obit_traceback_msg (err, routine,inUV2->name);

  /* Get input descriptors */
  in1Desc = inUV1->myDesc;
  in2Desc = inUV2->myDesc;

  /* Check compatability between inUV1, inUV2 */
  incompatible = (in1Desc->ncorr!=in2Desc->ncorr);
  incompatible = incompatible || (in1Desc->jlocs!=in2Desc->jlocs);
  incompatible = incompatible || (in1Desc->jlocf!=in2Desc->jlocf);
  incompatible = incompatible || (in1Desc->jlocif!=in2Desc->jlocif);
  incompatible = incompatible || (in1Desc->ilocb!=in2Desc->ilocb);
  if (incompatible) {
     Obit_log_error(err, OBIT_Error,"%s inUV1 and inUV2 have incompatible structures",
		   routine);
      return ;
  }

  /* Look at all data */
  btemp = TRUE;
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut(inUV1->info, "passAll", OBIT_bool, dim, &btemp);
  ObitInfoListAlwaysPut(inUV2->info, "passAll", OBIT_bool, dim, &btemp);
  
  /* test open output */
  oretCode = ObitUVOpen (outUV, OBIT_IO_WriteOnly, err);
  /* If this didn't work try OBIT_IO_ReadWrite */
  if ((oretCode!=OBIT_IO_OK) || (err->error)) {
    ObitErrClear(err);
    oretCode = ObitUVOpen (outUV, OBIT_IO_ReadWrite, err);
  }
  /* if it didn't work bail out */
  if ((oretCode!=OBIT_IO_OK) || (err->error)) {
    /* unset output buffer (may be multiply deallocated) */
    outUV->buffer = NULL;
    outUV->bufferSize = 0;
    Obit_traceback_msg (err, routine, outUV->name);
  }

  /* Copy tables before data if in1 and out are not the same */
  if (!ObitUVSame (inUV1, outUV, err)) {
    iretCode = ObitUVCopyTables (inUV1, outUV, exclude, NULL, err);
    /* If multisource out then copy SU table, multiple sources selected or
       sources deselected suggest MS out */
    if ((inUV1->mySel->numberSourcesList>1) || (!inUV1->mySel->selectSources))
      iretCode = ObitUVCopyTables (inUV1, outUV, NULL, sourceInclude, err);
    if (err->error) {
      outUV->buffer = NULL;
      outUV->bufferSize = 0;
      Obit_traceback_msg (err, routine, inUV1->name);
    }
  }
  if (err->error) Obit_traceback_msg (err, routine, inUV1->name);

  /* reset to beginning of uv data */
  iretCode = ObitUVIOSet (inUV1, err);
  oretCode = ObitUVIOSet (outUV, err);
  if (err->error) Obit_traceback_msg (err, routine,inUV1->name);

  /* Close and reopen input to init calibration which will have been disturbed 
     by the table copy */
  iretCode = ObitUVClose (inUV1, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error))
    Obit_traceback_msg (err, routine,inUV1->name);
  iretCode = ObitUVOpen (inUV1, OBIT_IO_ReadWrite, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error))
    Obit_traceback_msg (err, routine,inUV1->name);
  outUV->buffer = inUV1->buffer;

  outDesc = outUV->myDesc;   /* Get output descriptor */

  /* Read vis from second input */
  if (doCalSelect) iretCode = ObitUVReadSelect (inUV2, inUV2->buffer, err);
  else iretCode = ObitUVRead (inUV2, inUV2->buffer, err);
  if (err->error) Obit_traceback_msg (err, routine,inUV2->name);

  /* we're in business, subtract */
  while ((iretCode==OBIT_IO_OK) && (oretCode==OBIT_IO_OK)) {
    /* Read first input */
    iretCode = ObitUVRead (inUV1, inUV1->buffer, err);
    if (iretCode!=OBIT_IO_OK) break;
    /* How many */
    outDesc->numVisBuff = in1Desc->numVisBuff;
    firstVis = in1Desc->firstVis;

    /* Modify data */
    for (i=0; i<in1Desc->numVisBuff; i++) { /* loop over visibilities */
      indx = i*in1Desc->lrec+ in1Desc->nrparm;
      kndx = in2Desc->nrparm;
      for (j=0; j<in1Desc->ncorr; j++) { /* loop over correlations */
	/* subtract */
	if ((inUV1->buffer[indx+2]>0.0) && (inUV2->buffer[kndx+2]>0.0)) {
	  inUV1->buffer[indx]   -= inUV2->buffer[kndx];
	  inUV1->buffer[indx+1] -= inUV2->buffer[kndx+1];
	  /* this blanks poln data } else {
	    inUV1->buffer[indx]   = 0.0;
	    inUV1->buffer[indx+1] = 0.0;
	    inUV1->buffer[indx+2] = 0.0;*/
	}
	indx += in1Desc->inaxes[0];
	kndx += in2Desc->inaxes[0];
      } /* end loop over correlations */
      if (incompatible) break;
    } /* end loop over visibilities */
    
    /* Write */
    oretCode = ObitUVWrite (outUV, inUV1->buffer, err);
    /* suppress vis number update if rewriting the same file */
    if (same) {
      outUV->myDesc->firstVis = firstVis;
      ((ObitUVDesc*)(outUV->myIO->myDesc))->firstVis = firstVis;
    }
 } /* end loop processing data */
  
  /* Check for incompatibility */
  if (incompatible) {
    Obit_log_error(err, OBIT_Error,"%s inUV1 and inUV2 have incompatible contents",
		   routine);
    return;
  }
    
  /* check for errors */
  if ((iretCode > OBIT_IO_EOF) || (oretCode > OBIT_IO_EOF) ||
      (err->error)) /* add traceback,return */
    Obit_traceback_msg (err, routine,inUV1->name);
    
  /* unset output buffer (may be multiply deallocated ;'{ ) */
  outUV->buffer = NULL;
  outUV->bufferSize = 0;
  
  /* Reset passAll */
  btemp = FALSE;
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut(inUV1->info, "passAll", OBIT_bool, dim, &btemp);
  ObitInfoListAlwaysPut(inUV2->info, "passAll", OBIT_bool, dim, &btemp);
  
  /* close files */
  iretCode = ObitUVClose (inUV1, err);
  if (err->error) Obit_traceback_msg (err, routine, inUV1->name);
  
  iretCode = ObitUVClose (inUV2, err);
  if (err->error) Obit_traceback_msg (err, routine, inUV2->name);
  
  oretCode = ObitUVClose (outUV, err);
  if (err->error) Obit_traceback_msg (err, routine, outUV->name);
  
} /* end ObitUVUtilVisSub1 */

/**
 * Compare the visibilities in one ObitUV with those in another.
 * Return the RMS of the real and imaginary differences 
 * divided by the amplitude of inUV2.
 * Only valid visibilities compared, zero amplitudes ignored.
 * \param inUV1    Input uv data numerator, no calibration/selection
 *  Control parameter on info
 * \li printRat OBIT_float (1) If given and >0.0 then tell about entries
 *              with a real or imaginary ratio > printRat
 * \param inUV2    Input uv data denominator, no calibration/selection
 *                 inUV2 should have the same structure, no. vis etc
 *                 as inUV1.
 * \param err      Error stack, returns if not empty.
 * \return RMS of the real and imaginary differences /amplitude, 
 *        -1 => no data compared or other error.
 */
ofloat ObitUVUtilVisCompare (ObitUV *inUV1, ObitUV *inUV2, ObitErr *err)
{
  ObitIOCode iretCode;
  olong i, j, indx, jndx, vscnt;
  ollong count, vNo;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong NPIO;
  ofloat amp, rms = -1.0, rrat, irat, printRat=-1.0;
  odouble sum;
  gboolean incompatible;
  ObitUVDesc *in1Desc, *in2Desc;
  gchar *routine = "ObitUVUtilVisCompare";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return rms;
  g_assert (ObitUVIsA(inUV1));
  g_assert (ObitUVIsA(inUV2));

  /* Diagnostics? */
  ObitInfoListGetTest(inUV1->info, "printRat", &type, dim, &printRat);

  /* test open to fully instantiate input and see if it's OK */
  iretCode = ObitUVOpen (inUV1, OBIT_IO_ReadWrite, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine,inUV1->name, rms);

  /* Copy number of records per IO to second input */
  ObitInfoListGet (inUV1->info, "nVisPIO", &type, dim,  (gpointer)&NPIO, err);
  ObitInfoListPut (inUV2->info, "nVisPIO",  type, dim,  (gpointer)&NPIO, err);
  if (err->error) Obit_traceback_val (err, routine, inUV1->name, rms);

  /* Open second input */
  iretCode = ObitUVOpen (inUV2, OBIT_IO_ReadWrite, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) 
    Obit_traceback_val (err, routine,inUV2->name, rms);

  /* Get input descriptors */
  in1Desc = inUV1->myDesc;
  in2Desc = inUV2->myDesc;

  /* Check compatability between inUV1, inUV2 */
  incompatible = in1Desc->nvis!=in2Desc->nvis;
  incompatible = incompatible || (in1Desc->ncorr!=in2Desc->ncorr);
  incompatible = incompatible || (in1Desc->jlocs!=in2Desc->jlocs);
  incompatible = incompatible || (in1Desc->jlocf!=in2Desc->jlocf);
  incompatible = incompatible || (in1Desc->jlocif!=in2Desc->jlocif);
  incompatible = incompatible || (in1Desc->ilocb!=in2Desc->ilocb);
  if (incompatible) {
     Obit_log_error(err, OBIT_Error,"%s inUV1 and inUV2 have incompatible structures",
		   routine);
      return rms;
 }

  /* we're in business, loop comparing */
  count = 0;
  vscnt = 0;
  sum = 0.0;
  while (iretCode==OBIT_IO_OK) {
    /* Read first input */
    iretCode = ObitUVRead (inUV1, inUV1->buffer, err);
    if (iretCode!=OBIT_IO_OK) break;
    /* Read second input */
    iretCode = ObitUVRead (inUV2, inUV2->buffer, err);
    if (iretCode!=OBIT_IO_OK) break;

    /* compatability check */
    incompatible = in1Desc->numVisBuff!=in2Desc->numVisBuff;
    if (incompatible) break;

    /* Compare data */
    for (i=0; i<in1Desc->numVisBuff; i++) { /* loop over visibilities */
      vscnt++;
      /* compatability check - check time and baseline or antenna code */
      indx = i*in1Desc->lrec ;
      jndx = i*in2Desc->lrec ;
      incompatible = 
	inUV1->buffer[indx+in1Desc->iloct] !=inUV2->buffer[jndx+in2Desc->iloct];
      if (in1Desc->ilocb>=0) {
	incompatible = incompatible ||
	  inUV1->buffer[indx+in1Desc->ilocb] !=inUV2->buffer[jndx+in2Desc->ilocb];
      }
      if (in1Desc->iloca1>=0) {
	incompatible = incompatible ||
	  inUV1->buffer[indx+in1Desc->iloca1]!=inUV2->buffer[jndx+in2Desc->iloca1] ||
	  inUV1->buffer[indx+in1Desc->iloca2]!=inUV2->buffer[jndx+in2Desc->iloca2];
      }
      if (incompatible) {
	vNo = indx+in1Desc->firstVis + i;  /* Which visibility */
	if (inUV1->buffer[indx+in1Desc->iloct]!=inUV2->buffer[jndx+in2Desc->iloct])
	  Obit_log_error(err, OBIT_Error, "Incompatible Times %f != %f @ vis %ld", 
			 inUV1->buffer[indx+in1Desc->iloct], inUV2->buffer[jndx+in2Desc->iloct], vNo);
	if ((in1Desc->ilocb>=0) && (inUV1->buffer[indx+in1Desc->ilocb]!=inUV2->buffer[jndx+in2Desc->ilocb]))
	  Obit_log_error(err, OBIT_Error, "Incompatible Baselines %f != %f @ vis %ld", 
			 inUV1->buffer[indx+in1Desc->ilocb], inUV2->buffer[jndx+in2Desc->ilocb], vNo);
	if ((in1Desc->iloca2>=0) && (inUV1->buffer[indx+in1Desc->iloca1]!=inUV2->buffer[jndx+in2Desc->iloca1]))
	  Obit_log_error(err, OBIT_Error, "Incompatible Antenna 1 %f != %f @ vis %ld", 
			 inUV1->buffer[indx+in1Desc->ilocb], inUV2->buffer[jndx+in2Desc->ilocb], vNo);
	if ((in1Desc->iloca2>=0) && (inUV1->buffer[indx+in1Desc->iloca2]!=inUV2->buffer[jndx+in2Desc->iloca2]))
	  Obit_log_error(err, OBIT_Error, "Incompatible Antenna 2 %f != %f @ vis %ld", 
			 inUV1->buffer[indx+in1Desc->ilocb], inUV2->buffer[jndx+in2Desc->ilocb], vNo);
	break;
      }

      indx += in1Desc->nrparm;
      jndx += in2Desc->nrparm;
      for (j=0; j<in1Desc->ncorr; j++) { /* loop over correlations */
	/* Statistics  */
	amp = inUV2->buffer[jndx]*inUV2->buffer[jndx] + 
	  inUV2->buffer[jndx+1]*inUV2->buffer[jndx+1];
	if ((inUV1->buffer[indx+2]>0.0) && (inUV2->buffer[indx+2]>0.0) && (amp>0.0)) {
	  amp = sqrt(amp);
	  rrat = ((inUV1->buffer[indx] - inUV2->buffer[jndx]) *
		  (inUV1->buffer[indx] - inUV2->buffer[jndx])) / amp;
	  irat = ((inUV1->buffer[indx+1] - inUV2->buffer[jndx+1]) *
		  (inUV1->buffer[indx+1] - inUV2->buffer[jndx+1])) / amp;
	  sum += rrat + irat;
	  count += 2;
	  /* Diagnostics */
	  if ((printRat>0.0) && ((rrat>printRat) ||  (irat>printRat))) {
	    Obit_log_error(err, OBIT_InfoWarn, 
			   "High ratio vis %d corr %d ratio %f %f amp %f vis %f %f - %f %f",
			   vscnt, j+1, rrat, irat, amp, inUV1->buffer[indx], inUV1->buffer[indx+1],
			   inUV2->buffer[jndx], inUV2->buffer[jndx+1]);
	  }
	}
	indx += in1Desc->inaxes[0];
	jndx += in2Desc->inaxes[0];
      } /* end loop over correlations */
      if (incompatible) break;
    } /* end loop over visibilities */
  } /* end loop processing data */
  
  /* Check for incompatibility */
  if (incompatible) {
    Obit_log_error(err, OBIT_Error,"%s inUV1 and inUV2 have incompatible contents",
		   routine);
    return rms;
  }
    
  /* check for errors */
  if ((iretCode > OBIT_IO_EOF) || (err->error))
    Obit_traceback_val (err, routine,inUV1->name, rms);
    
  /* close files */
  iretCode = ObitUVClose (inUV1, err);
  if (err->error) Obit_traceback_val (err, routine, inUV1->name, rms);
  
  iretCode = ObitUVClose (inUV2, err);
  if (err->error) Obit_traceback_val (err, routine, inUV2->name, rms);

  /* Get RMS to return */
  if (count>0) rms = sqrt(sum / count);
  else rms = -1.0;

  return rms;
  
} /* end ObitUVUtilVisCompare */

/**
 * Reads the UV and rewrites its Index (AIPS NX) table
 * \param inUV    Input UV data. 
 * Control parameters are on the info member.
 * \li "maxScan"  OBIT_float (1,1,1) max. scan time, min. [def LARGE]
 * \li "maxGap"   OBIT_float (1,1,1) max. time gap in scan, min. [def LARGE]
 * \param err      Error stack, returns if not empty.
 */
void ObitUVUtilIndex (ObitUV *inUV, ObitErr *err)
{ 
  ObitIOCode retCode;
  ObitTableNX* table;
  ObitTableNXRow* row=NULL;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  olong num, i, iRow, ver, startVis=1, curVis, itemp, jtemp;
  olong suba, lastSubA=0, source, lastSource=0, fqid, lastFQID=0;
  ofloat maxScan, maxGap, *vis;
  odouble startTime=0.0, endTime=0.0, lastTime = -1.0e20; 
  gchar *routine = "ObitUVUtilIndex";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitUVIsA(inUV));

  /* Get maximum scan length */
  maxScan = 1.0e20;
  ObitInfoListGetTest(inUV->info, "maxScan", &type, dim, &maxScan);
  maxScan /= 1440.0; /* to days */

  /* Get maximum scan gap */
  maxGap = 1.0e20;
  ObitInfoListGetTest(inUV->info, "maxGap", &type, dim, &maxGap);
  maxGap /= 1440.0; /* to days */

  /* Open UV */
  inUV->bufferSize = MAX (0, inUV->bufferSize); /* Need buffer */
  retCode = ObitUVOpen (inUV, OBIT_IO_ReadWrite, err);
  if ((retCode != OBIT_IO_OK) || (err->error>0)) 
    Obit_traceback_msg (err, routine, inUV->name);
  lastSubA   = -1000; /* initialize subarray number */
  lastFQID   = -1000; /* initialize FQ Id */
  lastSource = -1000; /* initialize source number */
  curVis     = 0;     /* visibility counter */

  /* create Index table object */
  ver = 1;
  table = newObitTableNXValue ("Index table", (ObitData*)inUV, &ver, 
			       OBIT_IO_ReadWrite, err);
  if (err->error) Obit_traceback_msg (err, routine, inUV->name);
  
  /* Clear existing rows */
  ObitTableClearRows ((ObitTable*)table, err);
  if (err->error) goto cleanup;

  /* Open Index table */
  if ((ObitTableNXOpen (table, OBIT_IO_WriteOnly, err)
       != OBIT_IO_OK) || (err->error)) goto cleanup;

  /* Create Index Row */
  row = newObitTableNXRow (table);

  /* initialize row */
  row->SourID   = 0;
  row->SubA     = 0;
  row->Time     = 0.0;
  row->TimeI    = 0.0;
  row->StartVis = -1;
  row->EndVis   = -1;
  row->FreqID   = 0;
  fqid          = 1;
  source        = 1;

  /* attach to table buffer */
  ObitTableNXSetRow (table, row, err);
  if (err->error)  goto cleanup;

  /* Write at beginning of UVIndex Table */
  iRow = 0;
  table->myDesc->nrow = 0; /* ignore any previous entries */

  /* Loop over UV */
  while (retCode==OBIT_IO_OK) {
    retCode = ObitUVRead (inUV, inUV->buffer, err);
    /*if (retCode!=OBIT_IO_OK) fprintf(stderr, "retcode=%d\n", retCode);*/
    /* EOF is OK */
    if (retCode==OBIT_IO_EOF) ObitErrClear(err);
    if (retCode!=OBIT_IO_OK) break;
    if (err->error) goto cleanup;

    /* How many */
    num = inUV->myDesc->numVisBuff;
    
    /* Visibility pointer */
    vis = inUV->buffer;

    /* initialize on first visibility */
    if (curVis<=0) {
      startVis   = 1;
      startTime  = vis[inUV->myDesc->iloct];
      lastTime   = vis[inUV->myDesc->iloct];
    }
    
    /* Loop over buffer */
    for (i=0; i<num; i++) {
      curVis++; /* Current UV visibility number */

      /* require some time change for further checks */
      if (vis[inUV->myDesc->iloct]==lastTime) goto endloop;
      /* Subarray number */
      ObitUVDescGetAnts(inUV->myDesc, vis, &itemp, &jtemp, &suba);
      if (inUV->myDesc->ilocfq>=0) 
	fqid   = vis[inUV->myDesc->ilocfq] + 0.5;  /* Which FQ Id */
      if (inUV->myDesc->ilocsu>=0) 
	source = vis[inUV->myDesc->ilocsu]  + 0.5; /* Source number */
       /* Initialize? */
      if (lastFQID<=0)   lastFQID   = fqid;
      if (lastSubA<=0)   lastSubA   = suba;
      if (lastSource<=0) lastSource = source;
      
      
      /* New scan - change of source, fqid, or reach time limit */
      if (((suba!=lastSubA)     && (lastSubA>0)) || 
	  ((source!=lastSource) && (lastSource>0)) ||
	  ((fqid!=lastFQID)     && (lastFQID>0)) ||
	  ((vis[inUV->myDesc->iloct]-lastTime) > maxGap) ||
	  ((vis[inUV->myDesc->iloct]-startTime)>maxScan))
	{  /* Write index */

	/* Fill index record */
	if (lastSource>0) row->SourID   = lastSource;
	if (lastFQID>0)   row->FreqID   = lastFQID;
	if (lastSubA>0)   row->SubA     = lastSubA;
	row->Time     = 0.5 * (startTime + endTime);
	row->TimeI    = (endTime - startTime);
	row->StartVis = startVis;
	row->EndVis   = curVis-1;
	
	/* Write Index table */
	iRow++;
	if ((ObitTableNXWriteRow (table, iRow, row, err)
	     != OBIT_IO_OK) || (err->error>0)) goto cleanup;

	/* Initialize next suba */
	lastSubA   = suba;
	lastFQID   = fqid;
	lastSource = source;
	startVis   = curVis;
	startTime  = vis[inUV->myDesc->iloct];
	endTime    = vis[inUV->myDesc->iloct];
	
      } /* end of if new scan */

      endloop: endTime    = vis[inUV->myDesc->iloct];   /* potential end time */
      lastTime   = vis[inUV->myDesc->iloct];   /* time of last vis */
      vis += inUV->myDesc->lrec;     /* Update data visibility pointer */
    } /* end loop over buffer */

  } /* End loop over UV */

  /* Last Scan */
  if (lastSource>0) row->SourID   = lastSource;
  if (lastFQID>0)   row->FreqID   = lastFQID;
  if (lastSubA>0)   row->SubA     = lastSubA;
  row->Time     = 0.5 * (startTime + endTime);
  row->TimeI    = (endTime - startTime);
  row->StartVis = startVis;
  row->EndVis   = curVis;
  
  /* Write Index table */
  iRow++;
  if ((ObitTableNXWriteRow (table, iRow, row, err)
       != OBIT_IO_OK) || (err->error>0)) goto cleanup; 
  
  /* Close Index table */
  cleanup: if ((ObitTableNXClose (table, err) 
       != OBIT_IO_OK) || (err->error>0)) 
    Obit_traceback_msg (err, routine, table->name);

  /* Cleanup */
  row = ObitTableNXRowUnref(row);
  table = ObitTableNXUnref(table);
  if (err->error) Obit_traceback_msg (err, routine, inUV->name);

  /* Close UV */
  retCode = ObitUVClose (inUV, err);
  if ((retCode != OBIT_IO_OK) || (err->error>0)) 
    Obit_traceback_msg (err, routine, inUV->name);
} /* end ObitUVUtilIndex */


/**
 * Return a SourceList containing selected sources
 * Uses following criteria:
 * \li initially select all source in SU table
 *     If no SU table and info entry "Sources" is given, and is 
 *        the same to 8 characters to header "OBJECT", it is used.
 *     Else  header "OBJECT"
 * \li Any explicit selection in Sources parameter in UV data selector
 * \li selection by CalCode
 * \li selection to exclude any entries marked done in the PS table
 *     if Infolist item 'doPS' is True. 
 * \li if time range given and an NX table exists, only select 
 *     sources with data  in this timerange.
 * \param inUV   UV data with all selection criteria, other controls:
 * \li "doPS" OBIT_bool (1,1,1) If true and PS 1 exists, exclude [def False]
 * \li "souCode" OBIT_string (4,1,1) Source Cal code desired, '    ' => any code selected
 *                                   '*   ' => any non blank code (calibrators only)
 *                                   '-CAL' => blank codes only (no calibrators)
 *                                   def = '    '
 *   any sources marked done in PS table 1.
 * \param *err   ObitErr error stack.
 * \return requested ObitSourceList
 */
ObitSourceList* ObitUVUtilWhichSources (ObitUV *inUV, ObitErr *err) 
{
  olong iver, i, j, count;
  ObitTableSU *SUTable=NULL;
  ObitTableNX *NXTable=NULL;
  ObitTablePS *PSTable=NULL;
  ObitSourceList *out=NULL, *tsList=NULL;
  ObitInfoType type;
  olong theRow;
  gint32 dim[MAXINFOELEMDIM];
  gboolean want, allCal, allNCal, doTime, doPS, *good=NULL, OK;
  gchar souCode[5], *strTemp;
  gchar *routine = "ObitUVUtilWhichSources";

  /* Full Source list */
  iver = 1;
  SUTable = newObitTableSUValue (inUV->name, (ObitData*)inUV, &iver, 
				 OBIT_IO_ReadOnly, 0, err);
  if (err->error) Obit_traceback_val (err, routine, inUV->name, out);
  if (SUTable) {
    tsList = ObitTableSUGetList (SUTable, err);
    if (err->error) Obit_traceback_val (err, routine, inUV->name, out);
  } else {  /* Use position /name from header */
    out = ObitSourceListCreate ("SList", 1);
    /* Use first name in "Sources" if given, otherwise "object" from Header */
    OK = ObitInfoListGetP (inUV->info, "Sources", &type, dim, (gpointer)&strTemp);
    OK = OK && (strTemp!=NULL) && strncmp (strTemp, "    ", 4);
    OK = OK && !strncmp (strTemp, inUV->myDesc->object, 8);  /* And same as Object */
    if (OK) strncpy (out->SUlist[0]->SourceName, strTemp, MIN(16,UVLEN_VALUE));
    else    strncpy (out->SUlist[0]->SourceName, inUV->myDesc->object, MIN(20,UVLEN_VALUE));
    out->SUlist[0]->equinox = inUV->myDesc->equinox;
    out->SUlist[0]->RAMean  = inUV->myDesc->crval[inUV->myDesc->jlocr];
    out->SUlist[0]->DecMean = inUV->myDesc->crval[inUV->myDesc->jlocd];
    /* Compute apparent position */
    ObitPrecessUVJPrecessApp (inUV->myDesc, out->SUlist[0]);
    return out;
  }

  SUTable = ObitTableSUUnref(SUTable);   /* Done with table */

  /* check selection and calcode */
  souCode[0] = souCode[1] = souCode[2] = souCode[3] = ' ';souCode[4] = 0;
  ObitInfoListGetTest(inUV->info, "souCode", &type, dim, souCode);
  doPS = FALSE;
  ObitInfoListGetTest(inUV->info, "doPS",    &type, dim, &doPS);
  allCal  = !strncmp (souCode, "-CAL", 4); /* all calibrators */
  allNCal = !strncmp (souCode, "*   ", 4); /* All non calibrators */

  /* Need other tables? */
  if (doPS) {
    iver = 1;
    PSTable = newObitTablePSValue (inUV->name, (ObitData*)inUV, &iver, 
				   OBIT_IO_ReadOnly, err);
    if (!PSTable) doPS = FALSE;
    if (doPS) ObitTablePSOpen (PSTable, OBIT_IO_ReadOnly, err);
  }

  doTime = TRUE;
  if (doTime) {
    iver = 1;
    NXTable = newObitTableNXValue (inUV->name, (ObitData*)inUV, &iver, 
				   OBIT_IO_ReadOnly, err);
    if (!NXTable) doTime = FALSE;
    if (doTime) ObitTableNXOpen (NXTable, OBIT_IO_ReadOnly, err);
  }

  /* Keep track of the good ones */
  good = g_malloc0(tsList->number*sizeof(gboolean));
  for (i=0; i<tsList->number; i++) good[i] = TRUE;
  
  /* Loop over list */
  for (i=0; i<tsList->number; i++) {
    /* Explicitly stated */
    want = ObitUVSelWantSour(inUV->mySel, tsList->SUlist[i]->SourID);
    
    if (allCal || allNCal) {
      /* Calibrator/non calibrator */
      want = want  &&
	((allCal  && !strncmp(tsList->SUlist[i]->CalCode, "    ", 4)) ||
	 (allNCal &&  strncmp(tsList->SUlist[i]->CalCode, "    ", 4)));
    }
    
    /* Timerange */
    if (doTime) {
      want = want && 
	ObitTableNXWantSour (NXTable, tsList->SUlist[i]->SourID,
			     inUV->mySel->timeRange, err);
    }
    if (err->error) Obit_traceback_val (err, routine, inUV->name, out);
    
    /* Already done in PS table? */
    if (doPS) {
      want = want && 
	ObitTablePSWantSour (PSTable, tsList->SUlist[i]->SourceName, 
			     &theRow, err);
    }
    if (err->error) Obit_traceback_val (err, routine, inUV->name, out);
    
    /* Want this one? */
    if (!want) good[i] = FALSE;
  } /* End loop over list */

  /* Close tables */
  if (doPS) ObitTablePSClose (PSTable,  err);
  if (doTime) ObitTableNXClose (NXTable, err);
  if (NXTable) NXTable = ObitTableNXUnref(NXTable);   /* Done with table */
   if (PSTable) PSTable = ObitTablePSUnref(PSTable);  /* Done with table */
  
  
  /* create output with only valid entries 
     Count valid */
  count = 0;
  for (i=0; i<tsList->number; i++) if (good[i]) count++;
  out = ObitSourceListCreate ("Desired source list", count);
  
  /* Copy */
  j = 0;
  for (i=0; i<tsList->number; i++) {
    if (good[i]) {
      out->SUlist[j] = ObitSourceCopy(tsList->SUlist[i], out->SUlist[j], err);
      j++;
    }
  }
  
   tsList =  ObitSourceListUnref (tsList);  /* free old */
  if (good) g_free(good);
  return out;
} /* end ObitUVUtilWhichSources */

/**
 * Hanning smooth the data inObitUV.
 * \param inUV     Input uv data to average, 
 *                 Any request for calibration, editing and selection honored
 * Control parameters are on the info member.
 * \li "doDescm"  OBIT_bool (1,1,1) Descimate data after Hanning [def TRUE]
 *              
 * \param scratch  True if scratch file desired, will be same type as inUV.
 * \param outUV    If not scratch, then the previously defined output file
 *                 May be NULL for scratch only
 *                 If it exists and scratch, it will be Unrefed
 * \param err      Error stack, returns if not empty.
 * \return the frequency averaged ObitUV.
 */
ObitUV* ObitUVUtilHann (ObitUV *inUV, gboolean scratch, ObitUV *outUV, 
			ObitErr *err)
{
  ObitIOCode iretCode, oretCode;
  gboolean doCalSelect, doDescm;
  gchar *exclude[]={"AIPS CL", "AIPS SN", "AIPS FG", "AIPS CQ", "AIPS WX",
		    "AIPS AT", "AIPS CT", "AIPS OB", "AIPS IM", "AIPS MC",
		    "AIPS PC", "AIPS NX", "AIPS TY", "AIPS GC", "AIPS HI",
		    "AIPS PL", "AIPS NI", "AIPS BP", "AIPS OF", "AIPS PS",
		    "AIPS FQ", "AIPS SU", "AIPS AN", "AIPS PD", "AIPS OT",
		    NULL};
  gchar *sourceInclude[] = {"AIPS SU", NULL};
  olong NumChAvg, i, j, indx, jndx;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  ObitIOAccess access;
  ObitUVDesc *inDesc, *outDesc;
  gchar *today=NULL;
  ofloat *work=NULL, scale;
  gchar *routine = "ObitUVUtilHann";
 
  /* error checks */
  if (err->error) return outUV;
  g_assert (ObitUVIsA(inUV));
  if (!scratch && (outUV==NULL)) {
    Obit_log_error(err, OBIT_Error,"%s Output MUST be defined for non scratch files",
		   routine);
      return outUV;
  }

  /* Selection/calibration/editing of input? */
  doCalSelect = FALSE;
  ObitInfoListGetTest(inUV->info, "doCalSelect", &type, (gint32*)dim, &doCalSelect);
  if (doCalSelect) access = OBIT_IO_ReadCal;
  else access = OBIT_IO_ReadOnly;

  /* test open to fully instantiate input and see if it's OK */
  iretCode = ObitUVOpen (inUV, access, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine, inUV->name, outUV);

  /* Create scratch? */
  if (scratch) {
    if (outUV) outUV = ObitUVUnref(outUV);
    outUV = newObitUVScratch (inUV, err);
  } else { /* non scratch output must exist - clone from inUV */
    outUV->myDesc = ObitUVDescCopy (inUV->myDesc, outUV->myDesc, err);
    /*ObitUVClone (inUV, outUV, err);*/
  }
  if (err->error) Obit_traceback_val (err, routine, inUV->name, inUV);

  /* copy Descriptor */
  outUV->myDesc = ObitUVDescCopy(inUV->myDesc, outUV->myDesc, err);
 
  /* Descimate output? */
  doDescm = TRUE;
  ObitInfoListGetTest(inUV->info, "doDescm", &type, dim, &doDescm);

  /* Effectively averaging 2 channels if doDescm */
  if (doDescm)  NumChAvg = 2;
  else          NumChAvg = 1;
  dim[0] = dim[1] = dim[2] = dim[3] = dim[4] = 1;
  ObitInfoListAlwaysPut(inUV->info, "NumChAvg", OBIT_long, dim, &NumChAvg);
 
  /* Creation date today */
  today = ObitToday();
  strncpy (outUV->myDesc->date, today, UVLEN_VALUE-1);
  if (today) g_free(today);
  
  /* Get descriptors */
  inDesc  = inUV->myDesc;
  outDesc = outUV->myDesc;

  /* Create work array for averaging */
  work = g_malloc(2*inDesc->lrec*sizeof(ofloat));

  /* Modify descriptor for affects of averaging, get u,v,w scaling */
  ObitUVGetFreq (inUV, err);   /* Make sure frequencies updated */
  if (err->error) goto cleanup;
  /* Update output Descriptor */
  if (doDescm) {/* No descimate - basically moving up one (input) channel */
    outDesc->crpix[outDesc->jlocf]  = 0.5+inDesc->crpix[inDesc->jlocf]/2.0;
    outDesc->cdelt[outDesc->jlocf]  = inDesc->cdelt[inDesc->jlocf]*2;
    outDesc->inaxes[outDesc->jlocf] = inDesc->inaxes[inDesc->jlocf]/2;
  } else {  /* No descimate */
    outDesc->crpix[outDesc->jlocf]  = inDesc->crpix[inDesc->jlocf];
  }
  scale = 1.0;  /* Haven't really changed the frequency */
 
  /* If descimating, last channel incomplete - drop */
  if (doDescm) outDesc->inaxes[outDesc->jlocf]--;

  /* test open output */
  oretCode = ObitUVOpen (outUV, OBIT_IO_WriteOnly, err);
  /* If this didn't work try OBIT_IO_ReadWrite */
  if ((oretCode!=OBIT_IO_OK) || (err->error)) {
    ObitErrClear(err);
    oretCode = ObitUVOpen (outUV, OBIT_IO_ReadWrite, err);
  }
  /* if it didn't work bail out */
  if ((oretCode!=OBIT_IO_OK) || (err->error)) goto cleanup;

  /* Copy tables before data */
  iretCode = ObitUVCopyTables (inUV, outUV, exclude, NULL, err);
  /* If multisource out then copy SU table, multiple sources selected or
   sources deselected suggest MS out */
  if ((inUV->mySel->numberSourcesList>1) || (!inUV->mySel->selectSources))
  iretCode = ObitUVCopyTables (inUV, outUV, NULL, sourceInclude, err);
  /* FQ table selection */
  iretCode = ObitTableFQSelect (inUV, outUV, NULL, 0.0, err);
  /* Correct FQ table for averaging 
     FQSel (outUV, NumChAvg, 1, err); done in FQSelect(?) */
  if (err->error) goto cleanup;

  /* reset to beginning of uv data */
  iretCode = ObitIOSet (inUV->myIO,  inUV->info, err);
  oretCode = ObitIOSet (outUV->myIO, outUV->info, err);
  if (err->error) goto cleanup;

  /* Close and reopen input to init calibration which will have been disturbed 
     by the table copy */
  iretCode = ObitUVClose (inUV, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) goto cleanup;

  iretCode = ObitUVOpen (inUV, access, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) goto cleanup;

  /* we're in business, average data */
  while ((iretCode==OBIT_IO_OK) && (oretCode==OBIT_IO_OK)) {
    if (doCalSelect) iretCode = ObitUVReadSelect (inUV, inUV->buffer, err);
    else iretCode = ObitUVRead (inUV, inUV->buffer, err);
    if (iretCode!=OBIT_IO_OK) break;
    /* How many */
    outDesc->numVisBuff = inDesc->numVisBuff;

    /* Modify data */
    for (i=0; i<inDesc->numVisBuff; i++) { /* loop over visibilities */
      /* Copy random parameters */
      indx = i*inDesc->lrec;
      jndx = i*outDesc->lrec;
      for (j=0; j<inDesc->nrparm; j++) 
	outUV->buffer[jndx+j] =  inUV->buffer[indx+j];

      /* Scale u,v,w for new reference frequency */
      outUV->buffer[jndx+outDesc->ilocu] *= scale;
      outUV->buffer[jndx+outDesc->ilocv] *= scale;
      outUV->buffer[jndx+outDesc->ilocw] *= scale;

      /* Smooth data */
      indx += inDesc->nrparm;
      jndx += outDesc->nrparm;
      /* Average data */
      Hann (inUV->myDesc, outUV->myDesc, doDescm, 
	    &inUV->buffer[indx], &outUV->buffer[jndx], work, err);
      if (err->error) goto cleanup;
    } /* end loop over visibilities */

    /* Write */
    oretCode = ObitUVWrite (outUV, outUV->buffer, err);
    if (err->error) goto cleanup;
  } /* end loop processing data */
  
  /* check for errors */
  if ((iretCode > OBIT_IO_EOF) || (oretCode > OBIT_IO_EOF) ||
      (err->error)) goto cleanup;

  /* Cleanup */
 cleanup:
  if (work)  {g_free(work);}    work    = NULL;
  
  /* close files */
  iretCode = ObitUVClose (inUV, err);
  oretCode = ObitUVClose (outUV, err);
  if ((iretCode!=OBIT_IO_OK) || (oretCode!=OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, outUV->name, outUV);
  
  return outUV;
} /* end ObitUVUtilHann */

/**
 * Duplicate channels in input to output.
 * \param inUV     Input uv data to average, 
 *                 Any request for calibration, editing and selection honored
 * Control parameters are on the info member.
 * \li "unHann"  OBIT_bool (1,1,1) Undo previous Hanning? [def FALSE]
 * \li "nBloat"  OBIT_int  (1,1,1) Number of duplicates of each channel [def 2]
 *              
 * \param scratch  True if scratch file desired, will be same type as inUV.
 * \param outUV    If not scratch, then the previously defined output file
 *                 May be NULL for scratch only
 *                 If it exists and scratch, it will be Unrefed
 * \param err      Error stack, returns if not empty.
 * \return the frequency averaged ObitUV.
 */
ObitUV* ObitUVUtilBloat (ObitUV *inUV, gboolean scratch, ObitUV *outUV, 
			ObitErr *err)
{
  ObitIOCode iretCode, oretCode;
  gboolean doCalSelect, unHann;
  gchar *exclude[]={"AIPS CL", "AIPS SN", "AIPS FG", "AIPS CQ", "AIPS WX",
		    "AIPS AT", "AIPS CT", "AIPS OB", "AIPS IM", "AIPS MC",
		    "AIPS PC", "AIPS NX", "AIPS TY", "AIPS GC", "AIPS HI",
		    "AIPS PL", "AIPS NI", "AIPS BP", "AIPS OF", "AIPS PS",
		    "AIPS FQ", "AIPS SU", "AIPS AN", "AIPS PD", "AIPS SY",
		    "AIPS PT", "AIPS OT",
		    NULL};
  gchar *sourceInclude[] = {"AIPS SU", NULL};
  olong i, j, indx, jndx;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  ObitIOAccess access;
  ObitUVDesc *inDesc, *outDesc;
  gchar *today=NULL;
  olong nBloat=2;
  ofloat scale;
  gchar *routine = "ObitUVUtilBloat";
 
  /* error checks */
  if (err->error) return outUV;
  g_assert (ObitUVIsA(inUV));
  if (!scratch && (outUV==NULL)) {
    Obit_log_error(err, OBIT_Error,"%s Output MUST be defined for non scratch files",
		   routine);
      return outUV;
  }

  /* Selection/calibration/editing of input? */
  doCalSelect = FALSE;
  ObitInfoListGetTest(inUV->info, "doCalSelect", &type, (gint32*)dim, &doCalSelect);
  if (doCalSelect) access = OBIT_IO_ReadCal;
  else access = OBIT_IO_ReadOnly;

  /* test open to fully instantiate input and see if it's OK */
  iretCode = ObitUVOpen (inUV, access, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine, inUV->name, outUV);

  /* Create scratch? */
  if (scratch) {
    if (outUV) outUV = ObitUVUnref(outUV);
    outUV = newObitUVScratch (inUV, err);
  } else { /* non scratch output must exist - clone from inUV */
    outUV->myDesc = ObitUVDescCopy (inUV->myDesc, outUV->myDesc, err);
    /*ObitUVClone (inUV, outUV, err);*/
  }
  if (err->error) Obit_traceback_val (err, routine, inUV->name, inUV);

  /* copy Descriptor */
  outUV->myDesc = ObitUVDescCopy(inUV->myDesc, outUV->myDesc, err);
 
  /* How much bloat output? */
  nBloat = 2;
  ObitInfoListGetTest(inUV->info, "nBloat", &type, dim, &nBloat);
  unHann = FALSE;
  ObitInfoListGetTest(inUV->info, "unHann", &type, dim, &unHann);
  if (unHann) nBloat = 2;

  /* Creation date today */
  today = ObitToday();
  strncpy (outUV->myDesc->date, today, UVLEN_VALUE-1);
  if (today) g_free(today);
  
  /* Get descriptors */
  inDesc  = inUV->myDesc;
  outDesc = outUV->myDesc;

  /* Modify descriptor for affects of averaging, get u,v,w scaling */
  ObitUVGetFreq (inUV, err);   /* Make sure frequencies updated */
  scale = BloatFSetDesc (inDesc, outDesc, nBloat, err);
  if (err->error) goto cleanup;

  /* test open output */
  oretCode = ObitUVOpen (outUV, OBIT_IO_WriteOnly, err);
  /* If this didn't work try OBIT_IO_ReadWrite */
  if ((oretCode!=OBIT_IO_OK) || (err->error)) {
    ObitErrClear(err);
    oretCode = ObitUVOpen (outUV, OBIT_IO_ReadWrite, err);
  }
  /* if it didn't work bail out */
  if ((oretCode!=OBIT_IO_OK) || (err->error)) goto cleanup;

  /* Copy tables before data */
  iretCode = ObitUVCopyTables (inUV, outUV, exclude, NULL, err);
  /* If multisource out then copy SU table, multiple sources selected or
   sources deselected suggest MS out */
  if ((inUV->mySel->numberSourcesList>1) || (!inUV->mySel->selectSources))
  iretCode = ObitUVCopyTables (inUV, outUV, NULL, sourceInclude, err);
  /* FQ table selection */
  iretCode = ObitTableFQSelect (inUV, outUV, NULL, 0.0, err);
  /* Correct FQ table for averaging 
     FQSel (outUV, NumChAvg, 1, err); done in FQSelect(?) */
  if (err->error) goto cleanup;

  /* reset to beginning of uv data */
  iretCode = ObitIOSet (inUV->myIO,  inUV->info, err);
  oretCode = ObitIOSet (outUV->myIO, outUV->info, err);
  if (err->error) goto cleanup;

  /* Close and reopen input to init calibration which will have been disturbed 
     by the table copy */
  iretCode = ObitUVClose (inUV, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) goto cleanup;

  iretCode = ObitUVOpen (inUV, access, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) goto cleanup;

  /* we're in business, average data */
  while ((iretCode==OBIT_IO_OK) && (oretCode==OBIT_IO_OK)) {
    if (doCalSelect) iretCode = ObitUVReadSelect (inUV, inUV->buffer, err);
    else iretCode = ObitUVRead (inUV, inUV->buffer, err);
    if (iretCode!=OBIT_IO_OK) break;
    /* How many */
    outDesc->numVisBuff = inDesc->numVisBuff;

    /* Modify data */
    for (i=0; i<inDesc->numVisBuff; i++) { /* loop over visibilities */
      /* Copy random parameters */
      indx = i*inDesc->lrec;
      jndx = i*outDesc->lrec;
      for (j=0; j<inDesc->nrparm; j++) 
	outUV->buffer[jndx+j] =  inUV->buffer[indx+j];

      /* Scale u,v,w for new reference frequency */
      outUV->buffer[jndx+outDesc->ilocu] *= scale;
      outUV->buffer[jndx+outDesc->ilocv] *= scale;
      outUV->buffer[jndx+outDesc->ilocw] *= scale;

      /* Smooth data */
      indx += inDesc->nrparm;
      jndx += outDesc->nrparm;
      /* duplicate channels */
      Bloat (inUV->myDesc, outUV->myDesc, nBloat, unHann, 
	    &inUV->buffer[indx], &outUV->buffer[jndx], err);
      if (err->error) goto cleanup;
    } /* end loop over visibilities */

    /* Write */
    oretCode = ObitUVWrite (outUV, outUV->buffer, err);
    if (err->error) goto cleanup;
  } /* end loop processing data */
  
  /* check for errors */
  if ((iretCode > OBIT_IO_EOF) || (oretCode > OBIT_IO_EOF) ||
      (err->error)) goto cleanup;

  /* Cleanup */
 cleanup:
  /* close files */
  iretCode = ObitUVClose (inUV, err);
  oretCode = ObitUVClose (outUV, err);
  if ((iretCode!=OBIT_IO_OK) || (oretCode!=OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, outUV->name, outUV);
  
  return outUV;
} /* end ObitUVUtilBloat */

/**
 * Spectrally average the data inObitUV.
 * \param inUV     Input uv data to average, 
 *                 Any request for calibration, editing and selection honored
 * Control parameters on info element of inUV:
 * \li "NumChAvg" OBIT_long scalar Number of channels to average, [def. = all]
 * \li "doAvgAll" OBIT_bool Scalar, if TRUE then average all channels and 
 *                IF. default = FALSE
 * \li "ChanSel"  OBIT_int (4,*) Groups of channels to consider (relative to
 *                channels & IFs selected by BChan, EChan, BIF, EIF)
 *                (start, end, increment, IF) where start and end at the 
 *                beginning and ending channel numbers (1-rel) of the group
 *                to be included, increment is the increment between
 *                selected channels and IF is the IF number (1-rel)
 *                default increment is 1, IF=0 means all IF.
 *                The list of groups is terminated by a start <=0
 *                Default is all channels in each IF.
 * \li "noScale" OBIT_bool Scalar, if TRUE then do NOT scale u,v,w for new 
 *               frequency. Used for holography data.  def FALSE
 *              
 * \param scratch  True if scratch file desired, will be same type as inUV.
 * \param outUV    If not scratch, then the previously defined output file
 *                 May be NULL for scratch only
 *                 If it exists and scratch, it will be Unrefed
 * \param err      Error stack, returns if not empty.
 * \return the frequency averaged ObitUV.
 */
ObitUV* ObitUVUtilAvgF (ObitUV *inUV, gboolean scratch, ObitUV *outUV, 
			ObitErr *err)
{
  ObitIOCode iretCode, oretCode;
  gboolean doCalSelect;
  gchar *exclude[]={"AIPS CL", "AIPS SN", "AIPS FG", "AIPS CQ", "AIPS WX",
		    "AIPS AT", "AIPS CT", "AIPS OB", "AIPS IM", "AIPS MC",
		    "AIPS PC", "AIPS NX", "AIPS TY", "AIPS GC", "AIPS HI",
		    "AIPS PL", "AIPS NI", "AIPS BP", "AIPS OF", "AIPS PS",
		    "AIPS FQ", "AIPS SU", "AIPS AN", "AIPS PD", "AIPS SY",
		    "AIPS PT", "AIPS OT",
		    NULL};
  gchar *sourceInclude[] = {"AIPS SU", NULL};
  olong i, j, indx, jndx;
  olong *corChan=NULL, *corIF=NULL, *corStok=NULL;
  gboolean *corMask=NULL;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  ObitIOAccess access;
  ObitUVDesc *inDesc, *outDesc;
  gchar *today=NULL;
  ofloat *work=NULL, scale;
  olong NumChAvg, *ChanSel=NULL;
  gboolean doAvgAll, noScale;
  olong defSel[] = {1,-10,1,0, 0,0,0,0};
  gchar *routine = "ObitUVUtilAvgF";
 
  /* error checks */
  if (err->error) return outUV;
  g_assert (ObitUVIsA(inUV));
  if (!scratch && (outUV==NULL)) {
    Obit_log_error(err, OBIT_Error,"%s Output MUST be defined for non scratch files",
		   routine);
      return outUV;
  }

  /* Get Parameters */
  NumChAvg = -1;
  ObitInfoListGetTest(inUV->info, "NumChAvg", &type, dim, &NumChAvg);
  doAvgAll = FALSE;
  ObitInfoListGetTest(inUV->info, "doAvgAll", &type, dim, &doAvgAll);
  noScale = FALSE;
  ObitInfoListGetTest(inUV->info, "noScale", &type, dim, &noScale);
  ChanSel = NULL;
  if (!ObitInfoListGetP(inUV->info, "ChanSel", &type, dim, (gpointer)&ChanSel)) {
    ChanSel = defSel;  /* Use default = channels 1 => n */
  }
  /* ChanSel all zero? => default */
  if ((ChanSel[0]<=0) && (ChanSel[1]<=0) && (ChanSel[2]<=0) && (ChanSel[3]<=0)) {
    ChanSel = defSel;  /* Use default = channels 1 => n */
  }

  /* Selection/calibration/editing of input? */
  doCalSelect = FALSE;
  ObitInfoListGetTest(inUV->info, "doCalSelect", &type, (gint32*)dim, &doCalSelect);
  if (doCalSelect) access = OBIT_IO_ReadCal;
  else access = OBIT_IO_ReadOnly;

  /* test open to fully instantiate input and see if it's OK */
  iretCode = ObitUVOpen (inUV, access, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine, inUV->name, outUV);

  /* Create scratch? */
  if (scratch) {
    if (outUV) outUV = ObitUVUnref(outUV);
    outUV = newObitUVScratch (inUV, err);
  } else { /* non scratch output must exist - clone from inUV */
    outUV->myDesc = ObitUVDescCopy (inUV->myDesc, outUV->myDesc, err);
    /*ObitUVClone (inUV, outUV, err);*/
  }
  if (err->error) Obit_traceback_val (err, routine, inUV->name, inUV);

  /* copy Descriptor */
  outUV->myDesc = ObitUVDescCopy(inUV->myDesc, outUV->myDesc, err);
 
  /* Default frequency group selection?  Average all */
  if (ChanSel == defSel) defSel[1] = inUV->myDesc->inaxes[inUV->myDesc->jlocf];

  /* Default number of channels to average = all */
  if ((NumChAvg<=0) || doAvgAll) NumChAvg = inUV->myDesc->inaxes[inUV->myDesc->jlocf];
  NumChAvg = MIN (NumChAvg, inUV->myDesc->inaxes[inUV->myDesc->jlocf]);
 
  /* Creation date today */
  today = ObitToday();
  strncpy (outUV->myDesc->date, today, UVLEN_VALUE-1);
  if (today) g_free(today);
  
  /* Get descriptors */
  inDesc  = inUV->myDesc;
  outDesc = outUV->myDesc;

  /* Create work array for averaging */
  work = g_malloc(2*inDesc->lrec*sizeof(ofloat));
  /* Work arrays defining data */
  corChan = g_malloc(inDesc->ncorr*sizeof(olong));
  corIF   = g_malloc(inDesc->ncorr*sizeof(olong));
  corStok = g_malloc(inDesc->ncorr*sizeof(olong));
  corMask = g_malloc(inDesc->ncorr*sizeof(gboolean));

  /* Modify descriptor for affects of averaging, get u,v,w scaling */
  ObitUVGetFreq (inUV, err);   /* Make sure frequencies updated */
  scale = AvgFSetDesc (inDesc, outDesc, NumChAvg, ChanSel, doAvgAll, 
		       corChan, corIF, corStok, corMask, err);
  if (err->error) goto cleanup;
  /* Scale u,v,w? */
  if (noScale) scale = 1.0;

  /* test open output */
  oretCode = ObitUVOpen (outUV, OBIT_IO_WriteOnly, err);
  /* If this didn't work try OBIT_IO_ReadWrite */
  if ((oretCode!=OBIT_IO_OK) || (err->error)) {
    ObitErrClear(err);
    oretCode = ObitUVOpen (outUV, OBIT_IO_ReadWrite, err);
  }
  /* if it didn't work bail out */
  if ((oretCode!=OBIT_IO_OK) || (err->error)) goto cleanup;

  /* Copy tables before data */
  iretCode = ObitUVCopyTables (inUV, outUV, exclude, NULL, err);
  /* If multisource out then copy SU table, multiple sources selected or
   sources deselected suggest MS out */
  if ((inUV->mySel->numberSourcesList>1) || (!inUV->mySel->selectSources))
  iretCode = ObitUVCopyTables (inUV, outUV, NULL, sourceInclude, err);
  /* FQ table selection */
  iretCode = ObitTableFQSelect (inUV, outUV, NULL, 0.0, err);
  /* Correct FQ table for averaging 
     FQSel (outUV, NumChAvg, 1, err); done in FQSelect(?) */
  if (err->error) goto cleanup;

  /* reset to beginning of uv data */
  iretCode = ObitIOSet (inUV->myIO,  inUV->info, err);
  oretCode = ObitIOSet (outUV->myIO, outUV->info, err);
  if (err->error) goto cleanup;

  /* Close and reopen input to init calibration which will have been disturbed 
     by the table copy */
  iretCode = ObitUVClose (inUV, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) goto cleanup;

  iretCode = ObitUVOpen (inUV, access, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) goto cleanup;

  /* we're in business, average data */
  while ((iretCode==OBIT_IO_OK) && (oretCode==OBIT_IO_OK)) {
    if (doCalSelect) iretCode = ObitUVReadSelect (inUV, inUV->buffer, err);
    else iretCode = ObitUVRead (inUV, inUV->buffer, err);
    if (iretCode!=OBIT_IO_OK) break;
    /* How many */
    outDesc->numVisBuff = inDesc->numVisBuff;

    /* Modify data */
    for (i=0; i<inDesc->numVisBuff; i++) { /* loop over visibilities */
      /* Copy random parameters */
      indx = i*inDesc->lrec;
      jndx = i*outDesc->lrec;
      for (j=0; j<inDesc->nrparm; j++) 
	outUV->buffer[jndx+j] =  inUV->buffer[indx+j];

      /* Scale u,v,w for new reference frequency */
      outUV->buffer[jndx+outDesc->ilocu] *= scale;
      outUV->buffer[jndx+outDesc->ilocv] *= scale;
      outUV->buffer[jndx+outDesc->ilocw] *= scale;

      /* Average data */
      indx += inDesc->nrparm;
      jndx += outDesc->nrparm;
      /* Average data */
      AvgFAver (inUV->myDesc, outUV->myDesc, NumChAvg, ChanSel, doAvgAll, 
		corChan, corIF, corStok, corMask,
		&inUV->buffer[indx], &outUV->buffer[jndx], work, err);
      if (err->error) goto cleanup;
    } /* end loop over visibilities */

    /* Write */
    oretCode = ObitUVWrite (outUV, outUV->buffer, err);
    if (err->error) goto cleanup;
  } /* end loop processing data */
  
  /* check for errors */
  if ((iretCode > OBIT_IO_EOF) || (oretCode > OBIT_IO_EOF) ||
      (err->error)) goto cleanup;

  /* Cleanup */
 cleanup:
  if (work) {g_free(work);}       work    = NULL;
  if (corChan) {g_free(corChan);} corChan = NULL;
  if (corIF) {g_free(corIF);}     corIF   = NULL;
  if (corStok) {g_free(corStok);} corStok = NULL;
  if (corMask) {g_free(corMask);} corMask = NULL;
  
  /* close files */
  iretCode = ObitUVClose (inUV, err);
  oretCode = ObitUVClose (outUV, err);
  if ((iretCode!=OBIT_IO_OK) || (oretCode!=OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, outUV->name, outUV);
  
  return outUV;
} /* end ObitUVUtilAvgF */

/**
 * Temporally average the data inObitUV.
 * \param inUV     Input uv data to average, 
 *                 Any request for calibration, editing and selection honored
 * Control parameter on info element of inUV:
 * \li "timeAvg"   OBIT_float  (1,1,1) Time interval over which to average 
 *                 (min) [def = 1 min.]
 *                 NB: this should be at least 2 integrations.
 *              
 * \param scratch  True if scratch file desired, will be same type as inUV.
 * \param outUV    If not scratch, then the previously defined output file
 *                 May be NULL for scratch only
 *                 If it exists and scratch, it will be Unrefed
 * \param err      Error stack, returns if not empty.
 * \return the frequency averaged ObitUV.
 */
ObitUV* ObitUVUtilAvgT (ObitUV *inUV, gboolean scratch, ObitUV *outUV, 
			ObitErr *err)
{
  ObitIOCode iretCode, oretCode;
  gboolean doCalSelect;
  gchar *exclude[]={"AIPS CL", "AIPS SN", "AIPS FG", "AIPS CQ", "AIPS WX",
		    "AIPS AT", "AIPS CT", "AIPS OB", "AIPS IM", "AIPS MC",
		    "AIPS PC", "AIPS NX", "AIPS TY", "AIPS GC", "AIPS HI",
		    "AIPS PL", "AIPS NI", "AIPS BP", "AIPS OF", "AIPS PS",
		    "AIPS FQ", "AIPS SU", "AIPS AN", "AIPS PD", "AIPS SY",
		    "AIPS PT", "AIPS OT",
		    NULL};
  gchar *sourceInclude[] = {"AIPS SU", NULL};
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  olong ncorr, nrparm, numAnt, jtemp;
  ollong lltmp, i, j, numBL, jndx, indx, blindx, iindx=0, nvis;
  ollong *blLookup=NULL;
  ObitIOAccess access;
  ObitUVDesc *inDesc, *outDesc;
  ObitUVSortBuffer *outBuffer=NULL;
  olong suba, lastSourceID, curSourceID, lastSubA;
  gchar *today=NULL;
  ofloat timeAvg, curTime, startTime, endTime;
  ofloat *accVis=NULL, *accRP=NULL, *ttVis=NULL;
  ofloat *inBuffer;
  olong ant1, ant2;
  olong ivis=0, NPIO;
  gboolean done, gotOne;
  gchar *routine = "ObitUVUtilAvgT";
 
  /* error checks */
  if (err->error) return outUV;
  g_assert (ObitUVIsA(inUV));
  if (!scratch && (outUV==NULL)) {
    Obit_log_error(err, OBIT_Error,"%s Output MUST be defined for non scratch files",
		   routine);
      return outUV;
  }

  /* Get Parameter - Time interval */
  timeAvg = 1.0;  /* default 1 min */
  ObitInfoListGetTest(inUV->info, "timeAvg", &type, dim, &timeAvg);
  if (timeAvg<=(0.01/60.0)) timeAvg = 1.0;
  timeAvg /= 1440.0;  /* convert to days */

  /* Selection/calibration/editing of input? */
  doCalSelect = FALSE;
  ObitInfoListGetTest(inUV->info, "doCalSelect", &type, dim, &doCalSelect);
  if (doCalSelect) access = OBIT_IO_ReadCal;
  else access = OBIT_IO_ReadOnly;

  /* test open to fully instantiate input and see if it's OK */
  iretCode = ObitUVOpen (inUV, access, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine, inUV->name, outUV);

  /* Create scratch? */
  if (scratch) {
    if (outUV) outUV = ObitUVUnref(outUV);
    outUV = newObitUVScratch (inUV, err);
  } else { /* non scratch output must exist - clone from inUV */
    ObitUVClone (inUV, outUV, err);
  }
  if (err->error) Obit_traceback_val (err, routine, inUV->name, inUV);

  /* copy Descriptor */
  outUV->myDesc = ObitUVDescCopy(inUV->myDesc, outUV->myDesc, err);
  inBuffer = inUV->buffer;  /* Local copy of buffer pointer */
 
  /* Output creation date today */
  today = ObitToday();
  strncpy (outUV->myDesc->date, today, UVLEN_VALUE-1);
  if (today) g_free(today);
  
  /* Set number of output vis per read to twice number of baselines */
  suba    = 1;
  numAnt  = inUV->myDesc->numAnt[suba-1];/* actually highest antenna number */
  /* Better be some */
  Obit_retval_if_fail ((numAnt>1), err, outUV,
		       "%s Number of antennas NOT in descriptor",  
		       routine);  
  numBL   = (((ollong)numAnt)*(numAnt+1))/2;  /* Include auto correlations */
  NPIO = 1;
  ObitInfoListGetTest(inUV->info, "nVisPIO", &type, dim, &NPIO);
  jtemp = (olong)(2 * numBL);  /* Might cause trouble if numBL VERY large */
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut(outUV->info, "nVisPIO", OBIT_long, dim, &jtemp);

  /* test open output */
  oretCode = ObitUVOpen (outUV, OBIT_IO_WriteOnly, err);
  /* If this didn't work try OBIT_IO_ReadWrite */
  if ((oretCode!=OBIT_IO_OK) || (err->error)) {
    ObitErrClear(err);
    oretCode = ObitUVOpen (outUV, OBIT_IO_ReadWrite, err);
  }
  /* if it didn't work bail out */
  if ((oretCode!=OBIT_IO_OK) || (err->error)) goto cleanup;

  /* Get descriptors */
  inDesc  = inUV->myDesc;
  outDesc = outUV->myDesc;

  /* Create work arrays for averaging */
  ncorr   = inUV->myDesc->ncorr;
  nrparm  = inUV->myDesc->nrparm;
  lltmp   = 4*numBL*ncorr*sizeof(ofloat);
  accVis  = g_malloc0(lltmp);     /* Vis */
  lltmp   = numBL*(nrparm+1)*sizeof(ofloat);
  accRP   = g_malloc0(lltmp);  /* Rand. parm */
  ttVis   = g_malloc0(( inUV->myDesc->lrec+5)*sizeof(ofloat)); /* Temp Vis */

  /* Baseline lookup table */
  blLookup = g_malloc0 (numAnt* sizeof(ollong));
  blLookup[0] = 0;
  /* Include autocorr */
  for (i=1; i<numAnt; i++) blLookup[i] = blLookup[i-1] + numAnt-i+1; 

  /* Create sort buffer */
  /* Make sort buffer big enough for four copies of each baseline */
  nvis = 4 * numBL;  
  nvis = MIN (nvis, inUV->myDesc->nvis);
  outBuffer = ObitUVSortBufferCreate ("Buffer", outUV, nvis, err);
  if (err->error) goto cleanup;

  /* Copy tables before data */
  iretCode = ObitUVCopyTables (inUV, outUV, exclude, NULL, err);
  /* If multisource out then copy SU table, multiple sources selected or
   sources deselected suggest MS out */
  if ((inUV->mySel->numberSourcesList>1) || (!inUV->mySel->selectSources))
  iretCode = ObitUVCopyTables (inUV, outUV, NULL, sourceInclude, err);
  if (err->error) goto cleanup;

  /* reset to beginning of uv data */
  iretCode = ObitIOSet (inUV->myIO,  inUV->info, err);
  oretCode = ObitIOSet (outUV->myIO, outUV->info, err);
  if (err->error) goto cleanup;

  /* Close and reopen input to init calibration which will have been disturbed 
     by the table copy */
  iretCode = ObitUVClose (inUV, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) goto cleanup;

  iretCode = ObitUVOpen (inUV, access, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) goto cleanup;

  /* Initialize things */
  startTime = -1.0e20;
  endTime   =  1.0e20;
  lastSourceID = -1;
  curSourceID  = 0;
  outDesc->numVisBuff = 0;
  inBuffer  = inUV->buffer;   /* Local copy of buffer pointer */

  /* Loop over intervals */
  done   = FALSE;
  gotOne = FALSE;

  /* we're in business, average data */
  while ((iretCode==OBIT_IO_OK) && (oretCode==OBIT_IO_OK)) {
    if ((!gotOne) || (inUV->myDesc->numVisBuff<=0)) { /* need to read new record? */
      if (doCalSelect) iretCode = ObitUVReadSelect (inUV, inUV->buffer, err);
      else iretCode = ObitUVRead (inUV, inUV->buffer, err);
      if (err->error) goto cleanup;
    }

    /* Are we there yet??? */
    done = (inDesc->firstVis >= inDesc->nvis) || (iretCode==OBIT_IO_EOF);
    if (done && (startTime>0.0)) goto process; /* Final? */

    /* Make sure valid data found */
    if (inUV->myDesc->numVisBuff<=0) continue;

    /* loop over visibilities */
    for (ivis=0; ivis<inDesc->numVisBuff; ivis++) { 
      /* Copy random parameters */
      iindx = ivis*inDesc->lrec;

      gotOne = FALSE;
      
      curTime = inBuffer[iindx+inDesc->iloct]; /* Time */
      if (inDesc->ilocsu>=0) curSourceID = inBuffer[iindx+inDesc->ilocsu];
      if (startTime < -1000.0) {  /* Set time window etc. if needed */
	startTime = curTime;
	endTime   = startTime + timeAvg;
	lastSourceID = curSourceID;
      }

      /* Still in current interval/source? */
      if ((curTime<endTime) && (curSourceID == lastSourceID) && 
	  (inDesc->firstVis<=inDesc->nvis) && (iretCode==OBIT_IO_OK)) {
	/* accumulate */
	ObitUVDescGetAnts(inUV->myDesc, &inBuffer[iindx], &ant1, &ant2, &lastSubA);
	/* Check antenna number */
	Obit_retval_if_fail ((ant2<=numAnt), err, outUV, 
			     "%s Antenna 2=%d > max %d", routine, ant2, numAnt);  
	/* Baseline index this assumes a1<=a2 always */
	blindx =  blLookup[ant1-1] + ant2-ant1;
	
	/* Accumulate RP
	   (1,*) =  count 
	   (2...,*) =  Random parameters, sum u, v, w, time, int. */
	jndx = blindx*(1+nrparm);
	accRP[jndx]++;
	for (i=0; i<nrparm; i++) { 
	  /* Sum known parameters to average */
	  if ((i==inDesc->ilocu) || (i==inDesc->ilocv) || (i==inDesc->ilocw) ||
	      (i==inDesc->iloct) || (i==inDesc->ilocit)) {
	    accRP[jndx+i+1] += inBuffer[iindx+i];
	  } else { /* merely keep the rest */
	    accRP[jndx+i+1]  = inBuffer[iindx+i];
	  }
	} /* end loop over parameters */
	/* Accumulate Vis
	   (1,*) =  count 
	   (2,*) =  sum Real
	   (3,*) =  sum Imag
	   (4,*) =  Sum Wt     */
	indx = iindx+inDesc->nrparm; /* offset of start of vis data */
	for (i=0; i<ncorr; i++) {
	  if (inBuffer[indx+2] > 0.0) {
	    jndx = i*4 + blindx*4*ncorr;
	    accVis[jndx]   += 1.0;
	    accVis[jndx+1] += inBuffer[indx];
	    accVis[jndx+2] += inBuffer[indx+1];
	    accVis[jndx+3] += inBuffer[indx+2];
	  } 
	  indx += 3;
	} /* end loop over correlations */;
      } else {  /* process interval */
	
      process:
	/* Now may have the next record in the IO Buffer */
	if ((iretCode==OBIT_IO_OK) && (ivis<(inDesc->numVisBuff-1))) gotOne = TRUE;
	    
	/* Loop over baselines writing average */
	for (blindx=0; blindx<numBL; blindx++) {

	  /* Anything this baseline? */
	  jndx = blindx*(1+nrparm);
	  if (accRP[jndx]>0.0) {
	    /* Average u, v, w, time random parameters */
	    indx = 0;
	    for (i=0; i<nrparm; i++) { 
	      /* Average known parameters */
	      if ((i==inDesc->ilocu) || (i==inDesc->ilocv) || (i==inDesc->ilocw) ||
		  (i==inDesc->iloct)) 
		accRP[jndx+i+1] /= accRP[jndx];
	      /* Copy to output buffer */
	      ttVis[indx++] = accRP[jndx+i+1];
	    } /* End random parameter loop */

	    /* Average vis data */
	    for (j=0; j<ncorr; j++) {
	      jndx = j*4 + blindx*4*ncorr;
	      if (accVis[jndx]>0.0) {
		accVis[jndx+1] /= accVis[jndx];
		accVis[jndx+2] /= accVis[jndx];
	      }
	      /* Copy to output buffer */
	      ttVis[indx++] = accVis[jndx+1];
	      ttVis[indx++] = accVis[jndx+2];
	      ttVis[indx++] = accVis[jndx+3];
	    } /* end loop over correlators */

	    /* Copy to Sort Buffer (Sorts and writes when full) */
	    ObitUVSortBufferAddVis(outBuffer, ttVis, endTime, err);
	    if (err->error) goto cleanup;
	    /* Write output one vis at a time
	    outDesc->numVisBuff = 1;
	    oretCode = ObitUVWrite (outUV, outBuffer, err);
	    if (err->error) goto cleanup; */
	  } /* End any data this baseline */
	} /* end loop over baselines */
	
	/* Flush Sort Buffer */
	ObitUVSortBufferFlush (outBuffer, err);
	if (err->error) Obit_traceback_val (err, routine, outUV->name, outUV);

	/* Are we there yet??? */
	done = (inDesc->firstVis >= inDesc->nvis) || 
	  (iretCode==OBIT_IO_EOF);
	if (done) goto done;

	/* Reinitialize things */
	startTime = -1.0e20;
	endTime   =  1.0e20;
	for (i=0; i<4*ncorr*numBL; i++)    accVis[i] = 0.0;
	for (i=0; i<(nrparm+1)*numBL; i++) accRP[i] = 0.0;

	/* Now accumulate this visibility */
	ObitUVDescGetAnts(inUV->myDesc, &inBuffer[iindx], &ant1, &ant2, &lastSubA);
	/* Baseline index this assumes a1<=a2 always */
	/* Check antenna number */
	Obit_retval_if_fail ((ant2<=numAnt), err, outUV, 
			     "%s Antenna 2=%d > max %d", routine, ant2, numAnt);  
	blindx =  blLookup[ant1-1] + ant2-ant1;
	
	/* Accumulate RP
	   (1,*) =  count 
	   (2...,*) =  Random parameters, sum u, v, w, time, int. */
	jndx = blindx*(1+nrparm);
	accRP[jndx]++;
	for (i=0; i<nrparm; i++) { 
	  /* Sum known parameters to average */
	  if ((i==inDesc->ilocu) || (i==inDesc->ilocv) || (i==inDesc->ilocw) ||
	      (i==inDesc->iloct) || (i==inDesc->ilocit)) {
	    accRP[jndx+i+1] += inBuffer[iindx+i];
	  } else { /* merely keep the rest */
	    accRP[jndx+i+1]  = inBuffer[iindx+i];
	  }
	} /* end loop over parameters */
	/* Accumulate Vis
	   (1,*) =  count 
	   (2,*) =  sum Real
	   (3,*) =  sum Imag
	   (4,*) =  Sum Wt     */
	indx = iindx+inDesc->nrparm; /* offset of start of vis data */
	for (i=0; i<ncorr; i++) {
	  if (inBuffer[indx+2] > 0.0) {
	    jndx = i*4 + blindx*4*ncorr;
	    accVis[jndx]   += 1.0;
	    accVis[jndx+1] += inBuffer[indx];
	    accVis[jndx+2] += inBuffer[indx+1];
	    accVis[jndx+3] += inBuffer[indx+2];
	  } 
	  indx += 3;
	} /* end loop over correlations */;
      } /* end process interval */
      
    } /* end loop processing buffer of input data */
  } /* End loop over input file */
  
  /* End of processing */
 done:
  /* check for errors */
  if ((iretCode > OBIT_IO_EOF) || (oretCode > OBIT_IO_EOF) ||
      (err->error)) goto cleanup;
  
  /* Restore no vis per read in output */
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut (outUV->info, "nVisPIO", OBIT_long, dim, &NPIO);

  /* Cleanup */
 cleanup:
  if (accVis)   {g_free(accVis);}   accVis   = NULL;
  if (ttVis)    {g_free(ttVis);}    ttVis    = NULL;
  if (accRP)    {g_free(accRP);}    accRP    = NULL;
  if (blLookup) {g_free(blLookup);} blLookup = NULL;
  outBuffer = ObitUVSortBufferUnref(outBuffer);

  /* close files */
  iretCode = ObitUVClose (inUV, err);
  oretCode = ObitUVClose (outUV, err);
  if ((oretCode!=OBIT_IO_OK) || (iretCode!=OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, outUV->name, outUV);
  
  return outUV;
} /* end ObitUVUtilAvgT */

/**
 * Average all visibilities in inUV, write single vis to outUV.
 * For single source data single vis, for multisource data, one vis per scan.
 * Data labeled as baseline 1-2
 * \param inUV     Input uv data to average, 
 *                 Any request for calibration, editing and selection honored
 * \param scratch  True if scratch file desired, will be same type as inUV.
 * \param outUV    If not scratch, then the previously defined output file
 *                 May be NULL for scratch only
 *                 If it exists and scratch, it will be Unrefed
 * \param err      Error stack, returns if not empty.
 * \return the frequency averaged ObitUV.
 */
ObitUV* ObitUVUtilAvg2One (ObitUV *inUV, gboolean scratch, ObitUV *outUV, 
		  	   ObitErr *err)
{
  ObitIOCode iretCode, oretCode;
  gboolean doCalSelect;
  gchar *exclude[]={"AIPS CL", "AIPS SN", "AIPS FG", "AIPS CQ", "AIPS WX",
		    "AIPS AT", "AIPS CT", "AIPS OB", "AIPS IM", "AIPS MC",
		    "AIPS PC", "AIPS NX", "AIPS TY", "AIPS GC", "AIPS HI",
		    "AIPS PL", "AIPS NI", "AIPS BP", "AIPS OF", "AIPS PS",
		    "AIPS FQ", "AIPS SU", "AIPS AN", "AIPS PD", "AIPS SY",
		    "AIPS PT", "AIPS OT",
		    NULL};
  gchar *sourceInclude[] = {"AIPS SU", NULL};
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  olong ncorr, nrparm, numAnt, jtemp;
  ollong lltmp, i, j, numBL, jndx, indx, blindx, iindx=0;
  ObitIOAccess access;
  ObitUVDesc *inDesc, *outDesc;
  olong suba, lastSourceID, curSourceID;
  gchar *today=NULL;
  ofloat curTime, startTime, endTime;
  ofloat *accVis=NULL, *accRP=NULL, *ttVis=NULL;
  ofloat *inBuffer;
  olong ivis=0, NPIO;
  gboolean done, gotOne;
  gchar *routine = "ObitUVUtilAvgT";
 
  /* error checks */
  if (err->error) return outUV;
  g_assert (ObitUVIsA(inUV));
  if (!scratch && (outUV==NULL)) {
    Obit_log_error(err, OBIT_Error,"%s Output MUST be defined for non scratch files",
		   routine);
      return outUV;
  }

  /* Selection/calibration/editing of input? */
  doCalSelect = FALSE;
  ObitInfoListGetTest(inUV->info, "doCalSelect", &type, dim, &doCalSelect);
  if (doCalSelect) access = OBIT_IO_ReadCal;
  else access = OBIT_IO_ReadOnly;

  /* test open to fully instantiate input and see if it's OK */
  iretCode = ObitUVOpen (inUV, access, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine, inUV->name, outUV);

  /* Create scratch? */
  if (scratch) {
    if (outUV) outUV = ObitUVUnref(outUV);
    outUV = newObitUVScratch (inUV, err);
  } else { /* non scratch output must exist - clone from inUV */
    ObitUVClone (inUV, outUV, err);
  }
  if (err->error) Obit_traceback_val (err, routine, inUV->name, inUV);

  /* copy Descriptor */
  outUV->myDesc = ObitUVDescCopy(inUV->myDesc, outUV->myDesc, err);
  inBuffer = inUV->buffer;  /* Local copy of buffer pointer */
 
  /* Output creation date today */
  today = ObitToday();
  strncpy (outUV->myDesc->date, today, UVLEN_VALUE-1);
  if (today) g_free(today);
  
  /* Set number of output vis per read to twice number of baselines */
  suba    = 1;
  numAnt  = inUV->myDesc->numAnt[suba-1];/* actually highest antenna number */
  /* Better be some */
  Obit_retval_if_fail ((numAnt>1), err, outUV,
		       "%s Number of antennas NOT in descriptor",  
		       routine);  
  numBL   = 1;  /* Single output visibility */
  NPIO = 1;
  ObitInfoListGetTest(inUV->info, "nVisPIO", &type, dim, &NPIO);
  jtemp = (olong)(numBL); 
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut(outUV->info, "nVisPIO", OBIT_long, dim, &jtemp);

  /* test open output */
  oretCode = ObitUVOpen (outUV, OBIT_IO_WriteOnly, err);
  /* If this didn't work try OBIT_IO_ReadWrite */
  if ((oretCode!=OBIT_IO_OK) || (err->error)) {
    ObitErrClear(err);
    oretCode = ObitUVOpen (outUV, OBIT_IO_ReadWrite, err);
  }
  /* if it didn't work bail out */
  if ((oretCode!=OBIT_IO_OK) || (err->error)) goto cleanup;

  /* Get descriptors */
  inDesc  = inUV->myDesc;
  outDesc = outUV->myDesc;

  /* Create work arrays for averaging */
  ncorr   = inUV->myDesc->ncorr;
  nrparm  = inUV->myDesc->nrparm;
  lltmp   = 4*numBL*ncorr*sizeof(ofloat);
  accVis  = g_malloc0(lltmp);     /* Vis */
  lltmp   = numBL*(nrparm+1)*sizeof(ofloat);
  accRP   = g_malloc0(lltmp);  /* Rand. parm */
  ttVis   = g_malloc0(( inUV->myDesc->lrec+5)*sizeof(ofloat)); /* Temp Vis */

  /* Copy tables before data */
  iretCode = ObitUVCopyTables (inUV, outUV, exclude, NULL, err);
  /* If multisource out then copy SU table, multiple sources selected or
   sources deselected suggest MS out */
  if ((inUV->mySel->numberSourcesList>1) || (!inUV->mySel->selectSources))
  iretCode = ObitUVCopyTables (inUV, outUV, NULL, sourceInclude, err);
  if (err->error) goto cleanup;

  /* reset to beginning of uv data */
  iretCode = ObitIOSet (inUV->myIO,  inUV->info, err);
  oretCode = ObitIOSet (outUV->myIO, outUV->info, err);
  if (err->error) goto cleanup;

  /* Close and reopen input to init calibration which will have been disturbed 
     by the table copy */
  iretCode = ObitUVClose (inUV, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) goto cleanup;

  iretCode = ObitUVOpen (inUV, access, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) goto cleanup;

  /* Initialize things */
  startTime = -1.0e20;
  endTime   =  1.0e20;
  lastSourceID = -1;
  curSourceID  = 0;
  outDesc->numVisBuff = 0;
  inBuffer  = inUV->buffer;   /* Local copy of buffer pointer */

  /* Loop over intervals */
  done   = FALSE;
  gotOne = FALSE;

  /* we're in business, average data */
  while ((iretCode==OBIT_IO_OK) && (oretCode==OBIT_IO_OK)) {
    if ((!gotOne) || (inUV->myDesc->numVisBuff<=0)) { /* need to read new record? */
      if (doCalSelect) iretCode = ObitUVReadSelect (inUV, inUV->buffer, err);
      else iretCode = ObitUVRead (inUV, inUV->buffer, err);
    }

    /* Are we there yet??? */
    done = (inDesc->firstVis >= inDesc->nvis) || (iretCode==OBIT_IO_EOF);
    if (done && (startTime>0.0)) goto process; /* Final? */

    /* Make sure valid data found */
    if (inUV->myDesc->numVisBuff<=0) continue;

    /* loop over visibilities */
    for (ivis=0; ivis<inDesc->numVisBuff; ivis++) { 
      /* Copy random parameters */
      iindx = ivis*inDesc->lrec;

      gotOne = FALSE;
      
      curTime = inBuffer[iindx+inDesc->iloct]; /* Time */
      if (inDesc->ilocsu>=0) curSourceID = inBuffer[iindx+inDesc->ilocsu];
      if (startTime < -1000.0) {  /* Set time window etc. if needed */
	startTime = curTime;
	endTime   = startTime + 1000.0;
	lastSourceID = curSourceID;
      }

      /* Still in current interval/source? */
      if ((curTime<endTime) && (curSourceID == lastSourceID) && 
	  (inDesc->firstVis<=inDesc->nvis) && (iretCode==OBIT_IO_OK)) {
	/* accumulate all to one vis */
	blindx =  0;
	
	/* Accumulate RP
	   (1,*) =  count 
	   (2...,*) =  Random parameters, sum u, v, w, time, int. */
	jndx = blindx*(1+nrparm);
	accRP[jndx]++;
	for (i=0; i<nrparm; i++) { 
	  /* zero u,v,w */
	  if ((i==inDesc->ilocu) || (i==inDesc->ilocv) || (i==inDesc->ilocw))
	      inBuffer[iindx+i] = 0.0;
	  /* Sum known parameters to average */
	  if ((i==inDesc->iloct) || (i==inDesc->ilocit)) {
	    accRP[jndx+i+1] += inBuffer[iindx+i];
	  } else { /* merely keep the rest */
	    accRP[jndx+i+1]  = inBuffer[iindx+i];
	  }
	} /* end loop over parameters */
	/* Accumulate Vis
	   (1,*) =  count 
	   (2,*) =  sum Real*wt
	   (3,*) =  sum Imag*wt
	   (4,*) =  Sum Wt     */
	indx = iindx+inDesc->nrparm; /* offset of start of vis data */
	for (i=0; i<ncorr; i++) {
	  if (inBuffer[indx+2] > 0.0) {
	    jndx = i*4 + blindx*4*ncorr;
	    accVis[jndx]   += 1.0;
	    accVis[jndx+1] += inBuffer[indx]*inBuffer[indx+2];
	    accVis[jndx+2] += inBuffer[indx+1]*inBuffer[indx+2];
	    accVis[jndx+3] += inBuffer[indx+2];
	  } 
	  indx += 3;
	} /* end loop over correlations */;
      } else {  /* process interval */
	
      process:
	/* Now may have the next record in the IO Buffer */
	if ((iretCode==OBIT_IO_OK) && (ivis<(inDesc->numVisBuff-1))) gotOne = TRUE;
	    
	/* Loop over baselines writing average */
	for (blindx=0; blindx<numBL; blindx++) {

	  /* Anything this baseline? */
	  jndx = blindx*(1+nrparm);
	  if (accRP[jndx]>0.0) {
	    /* Average u, v, w, time random parameters */
	    indx = 0;
	    for (i=0; i<nrparm; i++) { 
	      /* Average known parameters */
	      if ((i==inDesc->ilocu) || (i==inDesc->ilocv) || (i==inDesc->ilocw) ||
		  (i==inDesc->iloct)) 
		accRP[jndx+i+1] /= accRP[jndx];
	      /* Copy to output buffer */
	      ttVis[indx++] = accRP[jndx+i+1];
	    } /* End random parameter loop */
	    /* Set baseline to 1-2 */
	    ObitUVDescSetAnts(outDesc, ttVis, 1, 2, 1);

	    /* Average vis data */
	    for (j=0; j<ncorr; j++) {
	      jndx = j*4 + blindx*4*ncorr;
	      if (accVis[jndx]>0.0) {
		accVis[jndx+1] /= accVis[jndx+3];
		accVis[jndx+2] /= accVis[jndx+3];
	      }
	      /* Copy to output buffer */
	      ttVis[indx++] = accVis[jndx+1];
	      ttVis[indx++] = accVis[jndx+2];
	      ttVis[indx++] = accVis[jndx+3];
	    } /* end loop over correlators */

	    /* Write output one vis at a time */
	    outDesc->numVisBuff = 1;
	    oretCode = ObitUVWrite (outUV, ttVis, err);
	    if (err->error) goto cleanup;
	  } /* End any data this baseline */
	} /* end loop over baselines */
	
	/* Are we there yet??? */
	done = (inDesc->firstVis >= inDesc->nvis) || 
	  (iretCode==OBIT_IO_EOF);
	if (done) goto done;

	/* Reinitialize things */
	startTime = -1.0e20;
	endTime   =  1.0e20;
	for (i=0; i<4*ncorr*numBL; i++)    accVis[i] = 0.0;
	for (i=0; i<(nrparm+1)*numBL; i++) accRP[i] = 0.0;

	/* Now accumulate this visibility */
	blindx = 0;
	
	/* Accumulate RP
	   (1,*) =  count 
	   (2...,*) =  Random parameters, sum u, v, w, time, int. */
	jndx = blindx*(1+nrparm);
	accRP[jndx]++;
	for (i=0; i<nrparm; i++) { 
	  /* Sum known parameters to average */
	  if ((i==inDesc->ilocu) || (i==inDesc->ilocv) || (i==inDesc->ilocw) ||
	      (i==inDesc->iloct) || (i==inDesc->ilocit)) {
	    accRP[jndx+i+1] += inBuffer[iindx+i];
	  } else { /* merely keep the rest */
	    accRP[jndx+i+1]  = inBuffer[iindx+i];
	  }
	} /* end loop over parameters */
	/* Accumulate Vis
	   (1,*) =  count 
	   (2,*) =  sum Real
	   (3,*) =  sum Imag
	   (4,*) =  Sum Wt     */
	indx = iindx+inDesc->nrparm; /* offset of start of vis data */
	for (i=0; i<ncorr; i++) {
	  if (inBuffer[indx+2] > 0.0) {
	    jndx = i*4 + blindx*4*ncorr;
	    accVis[jndx]   += 1.0;
	    accVis[jndx+1] += inBuffer[indx]*inBuffer[indx+2];
	    accVis[jndx+2] += inBuffer[indx+1]*inBuffer[indx+2];
	    accVis[jndx+3] += inBuffer[indx+2];
	  } 
	  indx += 3;
	} /* end loop over correlations */;
      } /* end process interval */
      
    } /* end loop processing buffer of input data */
  } /* End loop over input file */
  
  /* End of processing */
 done:
  /* check for errors */
  if ((iretCode > OBIT_IO_EOF) || (oretCode > OBIT_IO_EOF) ||
      (err->error)) goto cleanup;
  
  /* Restore no vis per read in output */
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut (outUV->info, "nVisPIO", OBIT_long, dim, &NPIO);

  /* Cleanup */
 cleanup:
  if (accVis)   {g_free(accVis);}   accVis   = NULL;
  if (ttVis)    {g_free(ttVis);}    ttVis    = NULL;
  if (accRP)    {g_free(accRP);}    accRP    = NULL;

  /* close files */
  iretCode = ObitUVClose (inUV, err);
  oretCode = ObitUVClose (outUV, err);
  if ((oretCode!=OBIT_IO_OK) || (iretCode!=OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, outUV->name, outUV);
  
  return outUV;
} /* end ObitUVUtilAvg2One */

/**
 * Spectrally smooth the data in ObitUV.
 * \param inUV     Input uv data to average, 
 *                 Any request for calibration, editing and selection honored
 * Control parameters on info element of inUV:
 * \li "NumChSmo" OBIT_long scalar Number of channels to average, [def. = 3]
 * \param scratch  True if scratch file desired, will be same type as inUV.
 * \param outUV    If not scratch, then the previously defined output file
 *                 May be NULL for scratch only
 *                 If it exists and scratch, it will be Unrefed
 * \param err      Error stack, returns if not empty.
 * \return the frequency averaged ObitUV.
 */
ObitUV* ObitUVUtilSmoF (ObitUV *inUV, gboolean scratch, ObitUV *outUV, 
			ObitErr *err)
{
  ObitIOCode iretCode, oretCode;
  gboolean doCalSelect;
  gchar *exclude[]={"AIPS CL", "AIPS SN", "AIPS FG", "AIPS CQ", "AIPS WX",
		    "AIPS AT", "AIPS CT", "AIPS OB", "AIPS IM", "AIPS MC",
		    "AIPS PC", "AIPS NX", "AIPS TY", "AIPS GC", "AIPS HI",
		    "AIPS PL", "AIPS NI", "AIPS BP", "AIPS OF", "AIPS PS",
		    "AIPS FQ", "AIPS SU", "AIPS AN", "AIPS PD", "AIPS SY",
		    "AIPS PT", "AIPS OT",
		    NULL};
  gchar *sourceInclude[] = {"AIPS SU", NULL};
  olong i, j, indx, jndx;
  olong *corChan=NULL, *corIF=NULL, *corStok=NULL;
  gboolean *corMask=NULL;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  ObitIOAccess access;
  ObitUVDesc *inDesc, *outDesc;
  gchar *today=NULL;
  ofloat *work=NULL;
  olong NumChSmo;
  gchar *routine = "ObitUVUtilSmoF";
 
  /* error checks */
  if (err->error) return outUV;
  g_assert (ObitUVIsA(inUV));
  if (!scratch && (outUV==NULL)) {
    Obit_log_error(err, OBIT_Error,"%s Output MUST be defined for non scratch files",
		   routine);
      return outUV;
  }

  /* Get Parameters */
  NumChSmo = 3;
  ObitInfoListGetTest(inUV->info, "NumChSmo", &type, dim, &NumChSmo);
  /* Make sure odd */
  Obit_retval_if_fail (((1+2*(NumChSmo/2))==NumChSmo), err, outUV,
		       "%s NumChSmo MUST be odd not %d",  
		       routine,NumChSmo);  

  /* Selection/calibration/editing of input? */
  doCalSelect = FALSE;
  ObitInfoListGetTest(inUV->info, "doCalSelect", &type, (gint32*)dim, &doCalSelect);
  if (doCalSelect) access = OBIT_IO_ReadCal;
  else access = OBIT_IO_ReadOnly;

  /* test open to fully instantiate input and see if it's OK */
  iretCode = ObitUVOpen (inUV, access, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine, inUV->name, outUV);

  /* Create scratch? */
  if (scratch) {
    if (outUV) outUV = ObitUVUnref(outUV);
    outUV = newObitUVScratch (inUV, err);
  } else { /* non scratch output must exist - clone from inUV */
    outUV->myDesc = ObitUVDescCopy (inUV->myDesc, outUV->myDesc, err);
    /*ObitUVClone (inUV, outUV, err);*/
  }
  if (err->error) Obit_traceback_val (err, routine, inUV->name, inUV);

  /* copy Descriptor */
  outUV->myDesc = ObitUVDescCopy(inUV->myDesc, outUV->myDesc, err);
 
  /* Don't smooth more channels than exist */
  NumChSmo = MIN (NumChSmo, inUV->myDesc->inaxes[inUV->myDesc->jlocf]);
 
  /* Creation date today */
  today = ObitToday();
  strncpy (outUV->myDesc->date, today, UVLEN_VALUE-1);
  if (today) g_free(today);
  
  /* Get descriptors */
  inDesc  = inUV->myDesc;
  outDesc = outUV->myDesc;

  /* Create work array for smoothing */
  work = g_malloc(2*inDesc->lrec*sizeof(ofloat));
  /* Work arrays defining data */
  corChan = g_malloc(inDesc->ncorr*sizeof(olong));
  corIF   = g_malloc(inDesc->ncorr*sizeof(olong));
  corStok = g_malloc(inDesc->ncorr*sizeof(olong));
  corMask = g_malloc(inDesc->ncorr*sizeof(gboolean));

  /* test open output */
  oretCode = ObitUVOpen (outUV, OBIT_IO_WriteOnly, err);
  /* If this didn't work try OBIT_IO_ReadWrite */
  if ((oretCode!=OBIT_IO_OK) || (err->error)) {
    ObitErrClear(err);
    oretCode = ObitUVOpen (outUV, OBIT_IO_ReadWrite, err);
  }
  /* if it didn't work bail out */
  if ((oretCode!=OBIT_IO_OK) || (err->error)) goto cleanup;

  /* Copy tables before data */
  iretCode = ObitUVCopyTables (inUV, outUV, exclude, NULL, err);
  /* If multisource out then copy SU table, multiple sources selected or
   sources deselected suggest MS out */
  if ((inUV->mySel->numberSourcesList>1) || (!inUV->mySel->selectSources))
  iretCode = ObitUVCopyTables (inUV, outUV, NULL, sourceInclude, err);
  /* FQ table selection */
  iretCode = ObitTableFQSelect (inUV, outUV, NULL, 0.0, err);
  if (err->error) goto cleanup;

  /* reset to beginning of uv data */
  iretCode = ObitIOSet (inUV->myIO,  inUV->info, err);
  oretCode = ObitIOSet (outUV->myIO, outUV->info, err);
  if (err->error) goto cleanup;

  /* Close and reopen input to init calibration which will have been disturbed 
     by the table copy */
  iretCode = ObitUVClose (inUV, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) goto cleanup;

  iretCode = ObitUVOpen (inUV, access, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) goto cleanup;

  /* we're in business, average data */
  while ((iretCode==OBIT_IO_OK) && (oretCode==OBIT_IO_OK)) {
    if (doCalSelect) iretCode = ObitUVReadSelect (inUV, inUV->buffer, err);
    else iretCode = ObitUVRead (inUV, inUV->buffer, err);
    if (iretCode!=OBIT_IO_OK) break;
    /* How many */
    outDesc->numVisBuff = inDesc->numVisBuff;

    /* Modify data */
    for (i=0; i<inDesc->numVisBuff; i++) { /* loop over visibilities */
      /* Copy random parameters */
      indx = i*inDesc->lrec;
      jndx = i*outDesc->lrec;
      for (j=0; j<inDesc->nrparm; j++) 
	outUV->buffer[jndx+j] =  inUV->buffer[indx+j];

      /* Average data */
      indx += inDesc->nrparm;
      jndx += outDesc->nrparm;
      /* Average data */
      SmooF (inDesc, outDesc, NumChSmo, corChan, corIF, corStok, corMask,
	     &inUV->buffer[indx], &outUV->buffer[jndx], work, err);
      if (err->error) goto cleanup;
    } /* end loop over visibilities */

    /* Write */
    oretCode = ObitUVWrite (outUV, outUV->buffer, err);
    if (err->error) goto cleanup;
  } /* end loop processing data */
  
  /* check for errors */
  if ((iretCode > OBIT_IO_EOF) || (oretCode > OBIT_IO_EOF) ||
      (err->error)) goto cleanup;

  /* Cleanup */
 cleanup:
  if (work) {g_free(work);}       work    = NULL;
  if (corChan) {g_free(corChan);} corChan = NULL;
  if (corIF) {g_free(corIF);}     corIF   = NULL;
  if (corStok) {g_free(corStok);} corStok = NULL;
  if (corMask) {g_free(corMask);} corMask = NULL;
  
  /* close files */
  iretCode = ObitUVClose (inUV, err);
  oretCode = ObitUVClose (outUV, err);
  if ((iretCode!=OBIT_IO_OK) || (oretCode!=OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, outUV->name, outUV);
  
  return outUV;
} /* end ObitUVUtilSmoF */

/**
 * Temporally average in a baseline dependent fashion the data in a ObitUV.
 * Also optionally average in frequency.
 * Time average UV data with averaging times depending
 * on time and baseline.  The averaging time is the greater of
 * maxInt and the time it takes for time smearing to reduce the 
 * visibility amplitude by maxFact.
 * \param inUV     Input uv data to average, 
 *                 Any request for calibration, editing and selection honored
 * Control parameters on info element of inUV:
 * \li "FOV"      OBIT_float  (1,1,1) Field of view (radius, deg)
 * \li "maxInt"   OBIT_float  (1,1,1) Maximum integration (min)
 * \li "maxFact"  OBIT_float  (1,1,1) Maximum time smearing factor
 * \li "NumChAvg" OBIT_long scalar Number of channels to average, [def. = all]
 * \li "doAvgAll" OBIT_bool Scalar, if TRUE then average all channels and 
 *                IF. default = FALSE
 * \li "ChanSel"  OBIT_int (4,*) Groups of channels to consider (relative to
 *                channels & IFs selected by BChan, EChan, BIF, EIF)
 *                (start, end, increment, IF) where start and end at the 
 *                beginning and ending channel numbers (1-rel) of the group
 *                to be included, increment is the increment between
 *                selected channels and IF is the IF number (1-rel)
 *                default increment is 1, IF=0 means all IF.
 *                The list of groups is terminated by a start <=0
 *                Default is all channels in each IF.
 *              
 * \param scratch  True if scratch file desired, will be same type as inUV.
 * \param outUV    If not scratch, then the previously defined output file
 *                 May be NULL for scratch only
 *                 If it exists and scratch, it will be Unrefed
 * \param err      Error stack, returns if not empty.
 * \return the frequency averaged ObitUV.
 */
ObitUV* ObitUVUtilBlAvgTF (ObitUV *inUV, gboolean scratch, ObitUV *outUV, 
			  ObitErr *err)
{
  ObitIOCode iretCode, oretCode;
  gboolean doCalSelect;
  gchar *exclude[]={"AIPS CL", "AIPS SN", "AIPS FG", "AIPS CQ", "AIPS WX",
		    "AIPS AT", "AIPS CT", "AIPS OB", "AIPS IM", "AIPS MC",
		    "AIPS PC", "AIPS NX", "AIPS TY", "AIPS GC", "AIPS HI",
		    "AIPS PL", "AIPS NI", "AIPS BP", "AIPS OF", "AIPS PS",
		    "AIPS FQ", "AIPS SU", "AIPS AN", "AIPS PD", "AIPS SY",
		    "AIPS PT", "AIPS OT", 
		    NULL};
  gchar *sourceInclude[] = {"AIPS SU", NULL};
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  olong ncorr, nrparm, numAnt;
  ollong lltmp, numBL, i, j, jndx, indx, blLo, blHi, nvis, ivis=0, iindx=0;
  ollong blindx=0, *blLookup=NULL;
  ObitIOAccess access;
  ObitUVDesc *inDesc, *outDesc;
  olong suba, lastSourceID, curSourceID, lastSubA;
  gchar *today=NULL;
  ofloat curTime=-1.0e20, startTime;
  ofloat *accVis=NULL, *accRP=NULL, *lsBlTime=NULL, *stBlTime=NULL, *stBlU=NULL, *stBlV=NULL;
  ofloat *tVis=NULL, *ttVis=NULL;
  ofloat *inBuffer;
  ObitUVSortBuffer *outBuffer=NULL;
  ofloat FOV, maxTime, maxInt, maxFact, UVDist2, maxUVDist2;
  olong ant1=1, ant2=2;
  olong NPIO, itemp, count=0;
  gboolean done, gotOne, doAllBl, sameInteg;
  ofloat *work=NULL, scale=1.0;
  olong NumChAvg, *ChanSel=NULL;
  gboolean doAvgAll, doAvgFreq;
  olong defSel[] = {1,1000000000,1,0, 0,0,0,0};
  olong *corChan=NULL, *corIF=NULL, *corStok=NULL;
  gboolean *corMask=NULL;
  gchar *routine = "ObitUVUtilAvgTF";
 
  /* error checks */
  if (err->error) return outUV;
  g_assert (ObitUVIsA(inUV));
  if (!scratch && (outUV==NULL)) {
    Obit_log_error(err, OBIT_Error,"%s Output MUST be defined for non scratch files",
		   routine);
      return outUV;
  }

  /* Give report */  
  Obit_log_error(err, OBIT_InfoErr, 
		 "Doing baseline dependent time averaging");
  ObitErrLog(err); 

  /* Get Parameters - radius of field of view */
  FOV = 20.0/60.0;  /* default 20 amin */
  ObitInfoListGetTest(inUV->info, "FOV", &type, dim, &FOV);
  /* to radians */
  FOV = MAX(FOV,1.0e-4) * DG2RAD;

  /* max. integration default 1 min */
  maxInt = 1.0;  
  ObitInfoListGetTest(inUV->info, "maxInt", &type, dim, &maxInt);
  if (maxInt<=(1.0e-2/60.0)) maxInt = 1.0;
  maxInt /= 1440.0;  /* convert to days */

  /* max amplitude loss default 1.01 */
  maxFact = 1.01;
  ObitInfoListGetTest(inUV->info, "maxFact", &type, dim, &maxFact);
  if (maxFact<0.99)  maxFact = 1.01;
  maxFact = MIN(MAX(maxFact,1.0), 10.0);

  /* Maximum UV distance squared to allow */
  maxUVDist2 = (InvSinc(1.0/maxFact) / FOV);
  maxUVDist2 = maxUVDist2*maxUVDist2; /* Square */

  /* Get Frequency Parameters */
  NumChAvg = 0;
  ObitInfoListGetTest(inUV->info, "NumChAvg", &type, dim, &NumChAvg);
  NumChAvg = MAX(1, NumChAvg);
  doAvgAll = FALSE;
  ObitInfoListGetTest(inUV->info, "doAvgAll", &type, dim, &doAvgAll);
  ChanSel = NULL;
  if (!ObitInfoListGetP(inUV->info, "ChanSel", &type, dim, (gpointer)&ChanSel)) {
    ChanSel = defSel;  /* Use default = channels 1 => n */
  }
  /* ChanSel all zero? => default */
  if ((ChanSel[0]<=0) && (ChanSel[1]<=0) && (ChanSel[2]<=0) && (ChanSel[3]<=0)) {
    ChanSel = defSel;  /* Use default = channels 1 => n */
  }

  /* Averaging in frequency? */
  doAvgFreq = (NumChAvg>1) || doAvgAll;

  /* Selection/calibration/editing of input? */
  doCalSelect = FALSE;
  ObitInfoListGetTest(inUV->info, "doCalSelect", &type, dim, &doCalSelect);
  if (doCalSelect) access = OBIT_IO_ReadCal;
  else             access = OBIT_IO_ReadOnly;

  /* test open to fully instantiate input and see if it's OK */
  iretCode = ObitUVOpen (inUV, access, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine, inUV->name, outUV);

  /* Is scratch? */
  if (scratch) {
    if (outUV) outUV = ObitUVUnref(outUV);
    outUV = newObitUVScratch (inUV, err);
  } else { /* non scratch output must exist - clone from inUV */
    outUV->myDesc = ObitUVDescCopy (inUV->myDesc, outUV->myDesc, err);
    ObitUVClone (inUV, outUV, err);
  }
  if (err->error) Obit_traceback_val (err, routine, inUV->name, outUV);

  inDesc  = inUV->myDesc;
  /* Create work array for frequency averaging */
  if (doAvgFreq) {
    work = g_malloc(2*inDesc->lrec*sizeof(ofloat));
    /* Work arrays defining data */
    corChan = g_malloc(inDesc->ncorr*sizeof(olong));
    corIF   = g_malloc(inDesc->ncorr*sizeof(olong));
    corStok = g_malloc(inDesc->ncorr*sizeof(olong));
    corMask = g_malloc(inDesc->ncorr*sizeof(gboolean));

    /* Modify descriptor for affects of frequency averaging, get u,v,w scaling */
    ObitUVGetFreq (inUV, err);   /* Make sure frequencies updated */
    scale = AvgFSetDesc (inUV->myDesc, outUV->myDesc, NumChAvg, ChanSel, doAvgAll, 
			 corChan, corIF, corStok, corMask, err);
    if (err->error) goto cleanup;
    
  } else { /* Only time averaging */   
    /* copy Descriptor */
    outUV->myDesc = ObitUVDescCopy(inUV->myDesc, outUV->myDesc, err);
  }

  /* Get Parameters - radius of field of view */
  FOV = 20.0/60.0;  /* default 20 amin */
  ObitInfoListGetTest(inUV->info, "FOV", &type, dim, &FOV);
  /* to radians */
  FOV = MAX(FOV,1.0e-4) * DG2RAD;

  /* max. integration default 1 min */
  maxInt = 1.0;  
  ObitInfoListGetTest(inUV->info, "maxInt", &type, dim, &maxInt);
  if (maxInt<=(1.0e-2/60.0)) maxInt = 1.0;
  maxInt /= 1440.0;  /* convert to days */

  /* max amplitude loss default 1.01 */
  maxFact = 1.01;
  ObitInfoListGetTest(inUV->info, "maxFact", &type, dim, &maxFact);
  if (maxFact<0.99)  maxFact = 1.01;
  maxFact = MIN(MAX(maxFact,1.0), 10.0);

  /* Maximum UV distance squared to allow */
  maxUVDist2 = (InvSinc(1.0/maxFact) / FOV);
  maxUVDist2 = maxUVDist2*maxUVDist2; /* Square */

  /* Get Frequency Parameters */
  NumChAvg = 0;
  ObitInfoListGetTest(inUV->info, "NumChAvg", &type, dim, &NumChAvg);
  NumChAvg = MAX(1, NumChAvg);
  doAvgAll = FALSE;
  ObitInfoListGetTest(inUV->info, "doAvgAll", &type, dim, &doAvgAll);
  ChanSel = NULL;
  if (!ObitInfoListGetP(inUV->info, "ChanSel", &type, dim, (gpointer)&ChanSel)) {
    ChanSel = defSel;  /* Use default = channels 1 => n */
  }
  /* ChanSel all zero? => default */
  if ((ChanSel[0]<=0) && (ChanSel[1]<=0) && (ChanSel[2]<=0) && (ChanSel[3]<=0)) {
    ChanSel = defSel;  /* Use default = channels 1 => n */
  }

  /* Averaging in frequency? */
  doAvgFreq = (NumChAvg>1) || doAvgAll;

  /* Selection/calibration/editing of input? */
  doCalSelect = FALSE;
  ObitInfoListGetTest(inUV->info, "doCalSelect", &type, dim, &doCalSelect);
  if (doCalSelect) access = OBIT_IO_ReadCal;
  else             access = OBIT_IO_ReadOnly;

  /* test open to fully instantiate input and see if it's OK */
  iretCode = ObitUVOpen (inUV, access, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine, inUV->name, outUV);

  /* Is scratch? */
  if (scratch) {
    if (outUV) outUV = ObitUVUnref(outUV);
    outUV = newObitUVScratch (inUV, err);
  } else { /* non scratch output must exist - clone from inUV */
    outUV->myDesc = ObitUVDescCopy (inUV->myDesc, outUV->myDesc, err);
    ObitUVClone (inUV, outUV, err);
  }
  if (err->error) Obit_traceback_val (err, routine, inUV->name, outUV);

  inDesc  = inUV->myDesc;
  /* Create work array for frequency averaging */
  if (doAvgFreq) {
    work = g_malloc(2*inDesc->lrec*sizeof(ofloat));
    /* Work arrays defining data */
    corChan = g_malloc(inDesc->ncorr*sizeof(olong));
    corIF   = g_malloc(inDesc->ncorr*sizeof(olong));
    corStok = g_malloc(inDesc->ncorr*sizeof(olong));
    corMask = g_malloc(inDesc->ncorr*sizeof(gboolean));

    /* Modify descriptor for affects of frequency averaging, get u,v,w scaling */
    ObitUVGetFreq (inUV, err);   /* Make sure frequencies updated */
    scale = AvgFSetDesc (inUV->myDesc, outUV->myDesc, NumChAvg, ChanSel, doAvgAll, 
			 corChan, corIF, corStok, corMask, err);
    if (err->error) goto cleanup;
    
  } else { /* Only time averaging */   
    /* copy Descriptor */
    outUV->myDesc = ObitUVDescCopy(inUV->myDesc, outUV->myDesc, err);
  }

  /* Add integration time if not present */
  if (outUV->myDesc->ilocit<0) {
    strncpy (outUV->myDesc->ptype[outUV->myDesc->nrparm], 
	     "INTTIM", UVLEN_KEYWORD-1);
    outUV->myDesc->ilocit = outUV->myDesc->nrparm;
    outUV->myDesc->nrparm++;
  }

  /* Output creation date today */
  today = ObitToday();
  strncpy (outUV->myDesc->date, today, UVLEN_VALUE-1);
  if (today) g_free(today);
  
  /* Set number of output vis per write  */
  NPIO = 1;
  ObitInfoListGetTest(inUV->info, "nVisPIO", &type, dim, &NPIO);
  itemp = 1000;  /* Internal IO buffer */
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut(outUV->info, "nVisPIO", OBIT_long, dim, &itemp);

  /* test open output */
  oretCode = ObitUVOpen (outUV, OBIT_IO_WriteOnly, err);
  /* If this didn't work try OBIT_IO_ReadWrite */
  if ((oretCode!=OBIT_IO_OK) || (err->error)) {
    ObitErrClear(err);
    oretCode = ObitUVOpen (outUV, OBIT_IO_ReadWrite, err);
  }
  /* if it didn't work bail out */
  if ((oretCode!=OBIT_IO_OK) || (err->error)) goto cleanup;

  /* Create sort buffer - Size depends on OS */
  if (sizeof(olong*)==4) {  /* 32 bit OS */
    /* Make sort buffer big  ~ 0.5 Gbyte */
    nvis = 500000000 / (outUV->myDesc->lrec*sizeof(ofloat));  
  } else if (sizeof(olong*)==8) {  /* 64 bit OS */
    /* Make sort buffer big  ~ 4 Gbyte */
    nvis = 4000000000 / (outUV->myDesc->lrec*sizeof(ofloat));  
  } else nvis = 2000000000 / (outUV->myDesc->lrec*sizeof(ofloat)); 
  nvis = MIN (nvis, inUV->myDesc->nvis);
  outBuffer = ObitUVSortBufferCreate ("Buffer", outUV, nvis, err);
  if (err->error) goto cleanup;

  /* Get descriptors */
  inDesc  = inUV->myDesc;
  outDesc = outUV->myDesc;

  /* Create work arrays for time averaging */
  suba    = 1;
  numAnt  = inDesc->numAnt[suba-1]; /* actually highest antenna number */
  /* Better be some */
  Obit_retval_if_fail ((numAnt>1), err, outUV,
		       "%s Number of antennas NOT in descriptor",  
		       routine);  
  numBL   = (((ollong)numAnt)*(numAnt+1))/2;  /* Include auto correlations */
  ncorr   = inDesc->ncorr;
  nrparm  = inDesc->nrparm;
  lltmp = 4*numBL*ncorr*sizeof(ofloat);
  accVis  = g_malloc0(lltmp);                           /* Vis accumulator */
  tVis    = g_malloc0((inDesc->lrec+5)*sizeof(ofloat)); /* Temp Vis */
  ttVis   = g_malloc0((inDesc->lrec+5)*sizeof(ofloat)); /* Temp Vis */
  lltmp   = numBL*(nrparm+1)*sizeof(ofloat);
  accRP   = g_malloc0(lltmp);   /* Rand. parm */
  lltmp   = numBL*sizeof(ofloat);
  stBlTime= g_malloc0(lltmp);   /* Baseline start time */
  lsBlTime= g_malloc0(lltmp);   /* Baseline last time */
  stBlU   = g_malloc0(lltmp);   /* Baseline start U */
  stBlV   = g_malloc0(lltmp);   /* Baseline start V */

  /* Baseline lookup table */
  blLookup = g_malloc0 (numAnt* sizeof(ollong));
  blLookup[0] = 0;
  /* Include autocorr */
  for (i=1; i<numAnt; i++) blLookup[i] = blLookup[i-1] + numAnt-i+1; 

 /* Copy tables before data */
  iretCode = ObitUVCopyTables (inUV, outUV, exclude, NULL, err);
  /* If multisource out then copy SU table, multiple sources selected or
   sources deselected suggest MS out */
  if ((inUV->mySel->numberSourcesList>1) || (!inUV->mySel->selectSources))
    iretCode = ObitUVCopyTables (inUV, outUV, NULL, sourceInclude, err);
  if (err->error) goto cleanup;

  /* reset to beginning of uv data */
  iretCode = ObitIOSet (inUV->myIO,  inUV->info, err);
  oretCode = ObitIOSet (outUV->myIO, outUV->info, err);
  if (err->error) goto cleanup;

  /* Close and reopen input to init calibration which will have been disturbed 
     by the table copy */
  iretCode = ObitUVClose (inUV, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) goto cleanup;

  iretCode = ObitUVOpen (inUV, access, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) goto cleanup;

  /* Initialize things */
  startTime    = -1.0e20;
  lastSourceID = -1;
  curSourceID  = 0;
  outDesc->numVisBuff = 0;
  inBuffer     = inUV->buffer;   /* Local copy of buffer pointer */

  /* Loop over intervals */
  done   = FALSE;
  gotOne = FALSE;

  /* we're in business, average data */
  while ((iretCode==OBIT_IO_OK) && (oretCode==OBIT_IO_OK)) {
    if ((!gotOne) || (inUV->myDesc->numVisBuff<=0)) { /* need to read new record? */
      if (doCalSelect) iretCode = ObitUVReadSelect (inUV, inUV->buffer, err);
      else iretCode = ObitUVRead (inUV, inUV->buffer, err);
    }

    /* Are we there yet??? */
    done = (inDesc->firstVis >= inDesc->nvis) || (iretCode==OBIT_IO_EOF);
    if (done && (startTime>0.0)) {doAllBl=TRUE; goto process;} /* Final? */

    /* Make sure valid data found */
    if (inUV->myDesc->numVisBuff<=0) continue;

    /* loop over visibilities */
    for (ivis=0; ivis<inDesc->numVisBuff; ivis++) { 

      gotOne = FALSE;
      
      /* Which data is this? */
      iindx = ivis*inDesc->lrec;
      ObitUVDescGetAnts(inUV->myDesc, &inBuffer[iindx], &ant1, &ant2, &lastSubA);
      /* Check antenna number */
      Obit_retval_if_fail ((ant2<=numAnt), err, outUV,
			   "%s Antenna 2=%d > max %d", routine, ant2, numAnt);  
      /* Baseline index this assumes a1<=a2 always */
      blindx =  blLookup[ant1-1] + ant2-ant1;
      blindx = MAX (0, MIN (blindx, numBL-1));
      curTime = inBuffer[iindx+inDesc->iloct]; /* Time */
      if (inDesc->ilocsu>=0) curSourceID = inBuffer[iindx+inDesc->ilocsu];

      /* Set time window etc. if needed */
      if (startTime < -1000.0) {  
	startTime    = curTime;
	lastSourceID = curSourceID;
      }

      /* If end of data, new scan, source, etc., finish all accumulations */
      doAllBl = 
	(curSourceID != lastSourceID) ||        /* Same source */
	(inDesc->firstVis>inDesc->nvis) ||      /* Not end of data */
	(iretCode!=OBIT_IO_OK);                 /* Not end of data */
      lastSourceID = curSourceID;
      
      /* Reset baseline start on first accumulation */
      jndx = blindx*(1+nrparm);
      if (accRP[jndx]<1.0) { 
	stBlTime[blindx] = inBuffer[iindx+inDesc->iloct];
	stBlU[blindx]    = inBuffer[iindx+inDesc->ilocu];
	stBlV[blindx]    = inBuffer[iindx+inDesc->ilocv];
      }
	
      /* Compute square of UV distance since start of integration */
      UVDist2 = 
	(inBuffer[iindx+inDesc->ilocu]-stBlU[blindx])*(inBuffer[iindx+inDesc->ilocu]-stBlU[blindx]) +
	(inBuffer[iindx+inDesc->ilocv]-stBlV[blindx])*(inBuffer[iindx+inDesc->ilocv]-stBlV[blindx]);

      /* Still in current baseline integration? */
      sameInteg = ((!doAllBl) &&                           /* Not end of scan or data */
		   ((curTime-stBlTime[blindx])<maxInt) &&  /* Max. integration */
		   (UVDist2<maxUVDist2));                  /* Max. smearing */
      if (!sameInteg) {  /* Write */
	
      process:
	/* Now may have the next record in the IO Buffer */
	if ((iretCode==OBIT_IO_OK) && (ivis<(inDesc->numVisBuff-1))) gotOne = TRUE;
	/* Finish this or all baselines? */
	if (doAllBl) {
	  blLo = 0;
	  blHi = numBL-1;
	} else { /* only this one */
	  blLo = blindx;
	  blHi = blindx;
	}
	/* Loop over this or all baselines */
	for (blindx=blLo; blindx<=blHi; blindx++) { 
	  /* Anything this baseline? */
	  jndx = blindx*(1+nrparm);
	  if (accRP[jndx]>0.0) {
	    /* Average u, v, w, time random parameters */
	    indx = 0;
	    for (i=0; i<nrparm; i++) { 
	      /* Average known parameters */
	      if ((i==inDesc->ilocu) || (i==inDesc->ilocv) || (i==inDesc->ilocw) ||
		  (i==inDesc->iloct)) 
		accRP[jndx+i+1] /= accRP[jndx];
	      /* Copy to output buffer */
	      ttVis[indx++] = accRP[jndx+i+1];
	    } /* End random parameter loop */
	    
	    /* Average vis data in time */
	    indx = outDesc->nrparm;
	    for (j=0; j<ncorr; j++) {
	      jndx = j*4 + blindx*4*ncorr;
	      if (accVis[jndx]>0.0) {
		accVis[jndx+1] /= accVis[jndx];
		accVis[jndx+2] /= accVis[jndx];
	      }
	      /* Copy to tempory vis */
	      tVis[indx++] = accVis[jndx+1];
	      tVis[indx++] = accVis[jndx+2];
	      tVis[indx++] = accVis[jndx+3];
	      
	    } /* end loop over correlators */
	    
	    if (doAvgFreq) {
	      /* Average data in frequency to output buffer */
	      AvgFAver (inDesc, outDesc, NumChAvg, ChanSel, doAvgAll, 
			corChan, corIF, corStok, corMask,
			&tVis[outDesc->nrparm], &ttVis[outDesc->nrparm], work, err);
	      if (err->error) goto cleanup;

	      /* Scale u,v,w for new reference frequency */
	      ttVis[outDesc->ilocu] *= scale;
	      ttVis[outDesc->ilocv] *= scale;
	      ttVis[outDesc->ilocw] *= scale;
	    } else { /* only time averaging */
	      /* Copy to output buffer */
	      indx = outDesc->nrparm;
	      jndx = outDesc->nrparm;
	      for (j=0; j<ncorr; j++) {
		ttVis[indx++] = tVis[jndx++];
		ttVis[indx++] = tVis[jndx++];
		ttVis[indx++] = tVis[jndx++];
	      }
	    }
	    
	    /* Set integration time */
	    if (inDesc->ilocit>=0)
	      ttVis[outDesc->ilocit] = 
		MAX (lsBlTime[blindx]-stBlTime[blindx],ttVis[outDesc->ilocit]);
	    else
	      ttVis[outDesc->ilocit] = lsBlTime[blindx]-stBlTime[blindx];
	    
	    /* Copy to Sort Buffer (Sorts and writes when full) */
	    count++;
	    maxTime = curTime - 0.6*maxInt;
	    ObitUVSortBufferAddVis(outBuffer, ttVis, maxTime, err);
	    if (err->error) goto cleanup;
	    
	    /* Reinitialize baseline */
	    jndx = blindx*(1+nrparm);
	    for (i=0; i<=nrparm; i++) accRP[jndx+i]  = 0.0;
	    jndx = blindx*4*ncorr;
	    for (i=0; i<4*ncorr; i++) accVis[jndx+i] = 0.0;
	    lsBlTime[blindx] = 0.0;
	    stBlTime[blindx] = 0.0;
	  } /* End any data this baseline */
	  
	} /* end loop over baseline */
      }  /* end process interval */
      
      /* accumulate */
      blindx =  blLookup[ant1-1] + ant2-ant1;
      jndx = blindx*(1+nrparm);
      if (accRP[jndx]<1.0) { /* starting conditions */
	stBlTime[blindx] = inBuffer[iindx+inDesc->iloct];
	stBlU[blindx]    = inBuffer[iindx+inDesc->ilocu];
	stBlV[blindx]    = inBuffer[iindx+inDesc->ilocv];
      }
      lsBlTime[blindx] = inBuffer[iindx+inDesc->iloct]; /* Highest time */
      /* Accumulate RP
	 (1,*)    =  count 
	 (2...,*) =  Random parameters, sum u, v, w, time, int. */
      accRP[jndx]++;
      for (i=0; i<nrparm; i++) { 
	/* Sum known parameters to average */
	if ((i==inDesc->ilocu) || (i==inDesc->ilocv) || (i==inDesc->ilocw) ||
	    (i==inDesc->iloct) || (i==inDesc->ilocit)) {
	  accRP[jndx+i+1] += inBuffer[iindx+i];
	} else { /* merely keep the rest */
	  accRP[jndx+i+1]  = inBuffer[iindx+i];
	}
      } /* end loop over parameters */
	/* Accumulate Vis
	   (1,*) =  count 
	   (2,*) =  sum Real
	   (3,*) =  sum Imag
	   (4,*) =  Sum Wt     */
      indx = iindx+inDesc->nrparm; /* offset of start of vis data */
      for (i=0; i<ncorr; i++) {
	if (inBuffer[indx+2] > 0.0) {
	  jndx = i*4 + blindx*4*ncorr;
	  accVis[jndx]   += 1.0;
	  accVis[jndx+1] += inBuffer[indx];
	  accVis[jndx+2] += inBuffer[indx+1];
	  accVis[jndx+3] += inBuffer[indx+2];
	} 
	indx += 3;
      } /* end loop over correlations */;
      
    } /* end loop processing buffer of input data */
  } /* End loop over input file */
  
  /* End of processing */

  /* check for errors */
  if ((iretCode > OBIT_IO_EOF) || (oretCode > OBIT_IO_EOF) ||
      (err->error)) goto cleanup;
  
  /* Cleanup */
 cleanup:
  if (accVis)   {g_free(accVis);}   accVis   = NULL;
  if (accRP)    {g_free(accRP);}    accRP    = NULL;
  if (tVis)     {g_free(tVis);}     tVis     = NULL;
  if (ttVis)    {g_free(ttVis);}    ttVis    = NULL;
  if (blLookup) {g_free(blLookup);} blLookup = NULL;
  if (lsBlTime) {g_free(lsBlTime);} lsBlTime = NULL;
  if (stBlTime) {g_free(stBlTime);} stBlTime = NULL;
  if (stBlU)    {g_free(stBlU);}    stBlU    = NULL;
  if (stBlV)    {g_free(stBlV);}    stBlV    = NULL;
  if (work)     {g_free(work);}     work     = NULL;
  if (corChan)  {g_free(corChan);}  corChan  = NULL;
  if (corIF)    {g_free(corIF);}    corIF    = NULL;
  if (corStok)  {g_free(corStok);}  corStok  = NULL;
  if (corMask)  {g_free(corMask);}  corMask  = NULL;
  
  /* Flush Sort Buffer */
  ObitUVSortBufferFlush (outBuffer, err);
  if (err->error) Obit_traceback_val (err, routine, outUV->name, outUV);
  outBuffer = ObitUVSortBufferUnref(outBuffer);
 
  /* close files */
  iretCode = ObitUVClose (inUV, err);
  oretCode = ObitUVClose (outUV, err);
  if ((oretCode!=OBIT_IO_OK) || (iretCode!=OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, outUV->name, outUV);

  /* Restore no vis per read in output */
  dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut (outUV->info, "nVisPIO", OBIT_long, dim, &NPIO);

  /* Give report */  
  Obit_log_error(err, OBIT_InfoErr, 
		 "Wrote %d averaged visibilities",count);
  ObitErrLog(err); 

  return outUV;
} /* end ObitUVUtilBlAvgTF */

/**
 * Count number of good correlations per time interval
 * \param inData    Input UV data, data selections, if any, applied
 * \param timeInt   Size of time interval in days, max. 500 intervals
 *                  If data from a new source is found a new interval 
 *                  is started.
 * \param err       Error stack, returns if not empty.
 * \return ObitInfoList with entries:
 * \li "numTime" OBIT_int [1] Number of time intervals
 * \li "numCorr" OBIT_int [1] Number of Correlations per vis
 * \li "Count"   OBIT_int [?] Count of good correlations per interval
 * \li "Bad"     OBIT_int [?] Count of bad correlations per interval
 * \li "Source"  OBIT_int [?] Source ID per interval
 * \li "LST"     OBIT_float [?] Average LST (days) per interval
 *               -1000.0 => no data.
 */
ObitInfoList* ObitUVUtilCount (ObitUV *inUV, ofloat timeInt, ObitErr *err)
{
  ObitInfoList *outList = NULL;
  ObitTableAN *ANTable=NULL;
  ObitIOCode iretCode;
  gboolean doCalSelect;
  olong i, ver, ivis;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitIOAccess access;
  ObitUVDesc *inDesc;
  ofloat *inBuffer;
  olong numTime, ncorr, indx, iindx, lastSourceID, curSourceID;
  gboolean gotOne, done, isVLA;
  odouble GSTiat0, DegDay, ArrayX, ArrayY;
  ofloat dataIat, ArrLong;
  ofloat startTime, endTime, curTime;
  ollong visCnt[500], goodCnt[500], badCnt[500], timeCnt[500], timeSou[500];
  ofloat timeSum[500];
  odouble dtemp[500];
  gchar *routine = "ObitUVUtilCount";

  /* error checks */
  if (err->error) return outList;
  g_assert (ObitUVIsA(inUV));

  /* Initialize sums */
  numTime = 0;
  for (i=0; i<500; i++) {
    visCnt[i] = goodCnt[i] = badCnt[i] = timeCnt[i] = timeSou[i] = 0;
    timeSum[i] = 0.0;
  }

  timeInt /= 1440.0;  /* timeInt to days */

  /* Selection of input? */
  doCalSelect = FALSE;
  ObitInfoListGetTest(inUV->info, "doCalSelect", &type, dim, &doCalSelect);
  if (doCalSelect) access = OBIT_IO_ReadCal;
  else access = OBIT_IO_ReadOnly;

  /* Open input */
  iretCode = ObitUVOpen (inUV, access, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine, inUV->name, outList);
  inDesc  = inUV->myDesc;  /* Get descriptor */

  /* Initialize things */
  startTime = -1.0e20;
  endTime   =  1.0e20;
  lastSourceID = -1;
  curSourceID  = 0;
  inBuffer = inUV->buffer;
  ncorr = inDesc->ncorr;

  /* Loop over intervals */
  done   = FALSE;
  gotOne = FALSE;

  /* we're in business, average data */
  while (iretCode==OBIT_IO_OK)  {
    if ((!gotOne) || (inDesc->numVisBuff<=0)) { /* need to read new record? */
      if (doCalSelect) iretCode = ObitUVReadSelect (inUV, inUV->buffer, err);
      else iretCode = ObitUVRead (inUV, inUV->buffer, err);
      if (iretCode > OBIT_IO_EOF) goto cleanup;
    }

    /* Are we there yet??? */
    done = (inDesc->firstVis >= inDesc->nvis) || (iretCode==OBIT_IO_EOF);
    if (done && (startTime>0.0)) goto process; /* Final? */

    /* Make sure valid data found */
    if (inUV->myDesc->numVisBuff<=0) continue;
    iindx = 0;

    /* loop over visibilities in buffer */
    for (ivis=0; ivis<inDesc->numVisBuff; ivis++) { 
      gotOne = FALSE;

      visCnt[numTime]++;  /* Count vis */
      curTime = inBuffer[iindx+inDesc->iloct]; /* Time */
      if (inDesc->ilocsu>=0) curSourceID = inBuffer[iindx+inDesc->ilocsu];
      if (startTime < -1000.0) {  /* Set time window etc. if needed */
	startTime = curTime;
	endTime   = startTime + timeInt;
	lastSourceID = curSourceID;
      }

      /* Still in current interval */
      if ((curTime<endTime) && (curSourceID == lastSourceID) && 
	  (inDesc->firstVis<=inDesc->nvis) && (iretCode==OBIT_IO_OK)) {

	/* sums */
	timeSou[numTime] = curSourceID;
	timeCnt[numTime]++;
	timeSum[numTime] += curTime;
	indx = iindx+inDesc->nrparm; /* offset of start of vis data */
	for (i=0; i<ncorr; i++) {
	  if (inBuffer[indx+2] > 0.0) goodCnt[numTime]++;
	  else badCnt[numTime]++;
	  indx += 3;
	} /* end loop over correlations */;
      } else {
	numTime++;  /* new interval */
	startTime = -1.0e20;
	endTime   =  1.0e20;
	/* Check over run */
	Obit_retval_if_fail ((numTime<100), err, outList,
		       "%s Too many time intervals %d",  
			     routine, numTime);  
      }
      iindx += inDesc->lrec;
    } /* end loop processing buffer of input data */
  } /* End loop over input file */
  
 process:
  /* Create output */
  outList  = newObitInfoList();

  /* Get time information */
  ver = 1;
  ANTable = newObitTableANValue (inUV->name, (ObitData*)inUV, &ver, 
				 OBIT_IO_ReadOnly, 0, 0, 0, err);
  GSTiat0 = ANTable->GSTiat0;
  DegDay  = ANTable->DegDay;
  ArrayX  = ANTable->ArrayX;
  ArrayY  = ANTable->ArrayY;
  if (!strncmp (ANTable->TimeSys, "IAT", 3)) {
    dataIat = 0.0;  /* in IAT */
  } else {  /* Assume UTC */
    dataIat = ANTable->dataUtc/86400.0;  /* in days */
  }
  isVLA  = !strncmp(ANTable->ArrName, "VLA     ", 8);
  ANTable = ObitTableANUnref(ANTable);   /* Done with table */
  if (err->error) Obit_traceback_val (err, routine, ANTable->name, outList);
  /* Need longitude */
  if (isVLA) ArrLong = 1.878283678;
  else ArrLong = atan2(ArrayY, ArrayX);
  ArrLong *= RAD2DG;

  /* Average times and convert to LST */
  numTime++;
  for (i=0; i<numTime; i++) {
    if (timeCnt[i]>0) {
      timeSum[i] /= timeCnt[i];
      /* To LST in deg */
      timeSum[i] = ((timeSum[i]-dataIat)*DegDay) + GSTiat0 + ArrLong;
      timeSum[i] /= 360.0;  /* Back to days */
    } else timeSum[i] = -1000.0;
  } /* end loop over time */

  /* save values */
  ObitInfoListAlwaysPut (outList, "numTime", OBIT_long, dim, &numTime);
  ObitInfoListAlwaysPut (outList, "numCorr", OBIT_long, dim, &ncorr);
  dim[0] = numTime;
  ObitInfoListAlwaysPut (outList, "Source", OBIT_long,     dim, timeSou);
  ObitInfoListAlwaysPut (outList, "LST",    OBIT_float,    dim, timeSum);
  for (i=0; i<numTime; i++) dtemp[i] = (odouble)goodCnt[i];
  ObitInfoListAlwaysPut (outList, "Count",  OBIT_double,   dim, dtemp);
  for (i=0; i<numTime; i++) dtemp[i] = (odouble)badCnt[i];
  ObitInfoListAlwaysPut (outList, "Bad",    OBIT_double,   dim, dtemp);
  for (i=0; i<numTime; i++) dtemp[i] = (odouble)visCnt[i];
  ObitInfoListAlwaysPut (outList, "Vis",    OBIT_double,   dim, dtemp);
  /* End of processing */

  /* Cleanup */
 cleanup:
  /* close file */
  iretCode = ObitUVClose (inUV, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, inUV->name, outList);
  
  return outList;
} /* end ObitUVUtilCount  */

/**
 * Copy blocks of channels from one UV to a set of output UVs..
 * All selected IFs and polarizations are also copied, there must be
 * an integran number of IFs per output.
 * \param inUV     Input uv data to average, 
 *                 Any request for calibration, editing and selection honored
 * \param nOut     Number of outout images
 * \param outUV    Array of previously defined but not yet instantiated 
 *                 (never opened) UV data objects to receive 1/nOut 
 *                 of the data channels in inUV.
 * \param err      Error stack, returns if not empty.
 */
void ObitUVUtilSplitCh (ObitUV *inUV, olong nOut, ObitUV **outUV, 
			ObitErr *err)
{
  ObitIOCode iretCode=OBIT_IO_SpecErr, oretCode=OBIT_IO_SpecErr;
  gboolean doCalSelect;
  gchar *exclude[]={"AIPS CL", "AIPS SN", "AIPS FG", "AIPS CQ", "AIPS WX",
		    "AIPS AT", "AIPS CT", "AIPS OB", "AIPS IM", "AIPS MC",
		    "AIPS PC", "AIPS NX", "AIPS TY", "AIPS GC", "AIPS HI",
		    "AIPS PL", "AIPS NI", "AIPS BP", "AIPS OF", "AIPS PS",
		    "AIPS FQ", "AIPS SU", "AIPS AN", "AIPS PD", "AIPS SY",
		    "AIPS PT", "AIPS OT",
		    NULL};
  gchar *sourceInclude[] = {"AIPS SU", NULL};
  olong *BChan=NULL, *numChan=NULL, *BIF=NULL, *numIF=NULL;
  olong chinc=1, nchan, nif, nchOut, NPIO;
  olong i, j, indx, jndx, ivis, nIFperOut, oldNumberIF, oldStartIF;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM]={1,1,1,1,1};
  ObitIOAccess access;
  ObitUVDesc *inDesc, *outDesc;
  gchar *today=NULL;
  ofloat *scale = NULL;
  gchar *routine = "ObitUVUtilSplitCh";
 
  /* error checks */
  if (err->error) return;
  g_assert (ObitUVIsA(inUV));
  for (i=0; i<nOut; i++) {
    if (!ObitUVIsA(outUV[i])) {
      Obit_log_error(err, OBIT_Error,"%s Output %d MUST be defined for non scratch files",
		     routine, i);
      return;
    }
  }

  /* Selection/calibration/editing of input? */
  doCalSelect = FALSE;
  ObitInfoListGetTest(inUV->info, "doCalSelect", &type, dim, &doCalSelect);
  if (doCalSelect) access = OBIT_IO_ReadCal;
  else access = OBIT_IO_ReadOnly;

  /* test open to fully instantiate input and see if it's OK */
  ObitUVFullInstantiate (inUV, TRUE, err);
  iretCode = ObitUVOpen (inUV, access, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error))
    Obit_traceback_msg (err, routine, inUV->name);

  /* Get descriptor */
  inDesc  = inUV->myDesc;

  /* Allocate work arrays */
  scale   = g_malloc(nOut*sizeof(ofloat));
  BChan   = g_malloc(nOut*sizeof(olong));
  numChan = g_malloc(nOut*sizeof(olong));
  BIF     = g_malloc(nOut*sizeof(olong));
  numIF   = g_malloc(nOut*sizeof(olong));

  /* Divvy up channels - must be an integral number of IFs per output -
     or only 1 */
  nchan = inDesc->inaxes[inDesc->jlocf];
  if (inDesc->jlocif>=0) nif = inDesc->inaxes[inDesc->jlocif];
  else                   nif = 1;
  nIFperOut = (glong) (0.9999 + (nif / (ofloat)nOut));
  if ((fabs(((ofloat)nIFperOut)-(nif / (ofloat)nOut))>0.001) && (nif>1)) {
    Obit_log_error(err, OBIT_Error,"%s Not an equal number of IFs per output",
		   routine);
    return;
  }
  if (nOut>(nchan* nIFperOut)) {
    Obit_log_error(err, OBIT_Error,"%s Fewer channels, %d than output files %d",
		   routine, nchan, nOut);
    return;
  }
  nchOut = (glong) (0.999 + ((nchan*nif) / (ofloat)nOut)); 
  nchOut = MAX (1, nchOut);
  nchOut = MIN (nchOut, nchan);
  for (i=0; i<nOut; i++) {
    BIF[i]     = 1 + i*nIFperOut;
    numIF[i]   = nIFperOut;
    BChan[i]   = 1 + i*nchOut - (BIF[i]-1)*nchan;
    numChan[i] = MIN (nchOut, nchan-BChan[i]+1);
  }

  /* Set up output UV data */
  for (i=0; i<nOut; i++) {
 
    /* copy Descriptor */
    outUV[i]->myDesc = ObitUVDescCopy(inUV->myDesc, outUV[i]->myDesc, err);
    
    /* Creation date today */
    today = ObitToday();
    strncpy (outUV[i]->myDesc->date, today, UVLEN_VALUE-1);
    if (today) g_free(today);
    
    /* Get descriptor */
    outDesc = outUV[i]->myDesc;

    /* Set output frequency info */
    outDesc->crval[outDesc->jlocf]  = inDesc->crval[outDesc->jlocf] + 
      (BChan[i]-inDesc->crpix[inDesc->jlocf]) * inDesc->cdelt[inDesc->jlocf] +
      (inDesc->freqIF[BIF[i]-1] - inDesc->freqIF[0]);
    outDesc->inaxes[outDesc->jlocf] = numChan[i];
    /*outDesc->crpix[outDesc->jlocf]  = 1.0;*/
    outDesc->cdelt[outDesc->jlocf]  = inDesc->cdelt[inDesc->jlocf] * chinc;
    /* UVW scaling parameter */
    scale[i] = outDesc->crval[outDesc->jlocf]/inDesc->freq;
    /* IFs */
    if (outDesc->jlocif>=0) {
      outDesc->crval[outDesc->jlocif]  = inDesc->crval[outDesc->jlocif] + 
	(BIF[i]-inDesc->crpix[inDesc->jlocif]) * inDesc->cdelt[inDesc->jlocif];
      outDesc->inaxes[outDesc->jlocif] = numIF[i];
      outDesc->crpix[outDesc->jlocif]  = 1.0;
      outDesc->cdelt[outDesc->jlocif]  = inDesc->cdelt[inDesc->jlocif] * chinc;
    }
    /* Alternate frequency/vel */
    outDesc->altCrpix = inDesc->altCrpix - (BChan[i] + 1.0)/chinc;
    outDesc->altRef   = inDesc->altRef;

    /* Copy number of records per IO to output */
    ObitInfoListGet (inUV->info, "nVisPIO", &type, dim, (gpointer)&NPIO, err);
    ObitInfoListAlwaysPut (outUV[i]->info, "nVisPIO", type, dim, (gpointer)&NPIO);
    if (err->error) goto cleanup;

    /* test open output */
    oretCode = ObitUVOpen (outUV[i], OBIT_IO_WriteOnly, err);
    /* If this didn't work try OBIT_IO_ReadWrite */
    if ((oretCode!=OBIT_IO_OK) || (err->error)) {
      ObitErrClear(err);
      oretCode = ObitUVOpen (outUV[i], OBIT_IO_ReadWrite, err);
    }
    /* if it didn't work bail out */
    if ((oretCode!=OBIT_IO_OK) || (err->error)) goto cleanup;

    /* Copy tables before data */
    iretCode = ObitUVCopyTables (inUV, outUV[i], exclude, NULL, err);
    /* If multisource out then copy SU table, multiple sources selected or
       sources deselected suggest MS out */
    if ((inUV->mySel->numberSourcesList>1) || (!inUV->mySel->selectSources))
      iretCode = ObitUVCopyTables (inUV, outUV[i], NULL, sourceInclude, err);
    /* Fiddle IF selection */
    oldNumberIF = inUV->mySel->numberIF;
    oldStartIF  = inUV->mySel->startIF;
    inUV->mySel->numberIF = nIFperOut;
    inUV->mySel->startIF  = BIF[i];
    ObitTableFQSelect (inUV, outUV[i], NULL, 0.0, err);
    /* reset IF selection */
    inUV->mySel->numberIF = oldNumberIF;
    inUV->mySel->startIF  = oldStartIF;
    if (err->error) goto cleanup;
   
    /* reset to beginning of uv data */
    iretCode = ObitIOSet (inUV->myIO,  inUV->info, err);
    oretCode = ObitIOSet (outUV[i]->myIO, outUV[i]->info, err);
    if (err->error) goto cleanup;

    /* Close and reopen input to init calibration which will have been disturbed 
       by the table copy */
    iretCode = ObitUVClose (inUV, err);
    if ((iretCode!=OBIT_IO_OK) || (err->error)) goto cleanup;
    
    iretCode = ObitUVOpen (inUV, access, err);
    if ((iretCode!=OBIT_IO_OK) || (err->error)) goto cleanup;
    
  } /* end loop over output files */

  /* we're in business, copy data */
  while ((iretCode==OBIT_IO_OK) && (oretCode==OBIT_IO_OK)) {
    if (doCalSelect) iretCode = ObitUVReadSelect (inUV, inUV->buffer, err);
    else iretCode = ObitUVRead (inUV, inUV->buffer, err);
    if (iretCode!=OBIT_IO_OK) break;

    /* How many */
    for (i=0; i<nOut; i++) {
      outUV[i]->myDesc->numVisBuff = inDesc->numVisBuff;
    }

    /* Copy data */
    for (ivis=0; ivis<inDesc->numVisBuff; ivis++) { /* loop over visibilities */
      /* Copy random parameters */
      indx = ivis*inDesc->lrec;
      for (i=0; i<nOut; i++) {
	jndx = ivis*outUV[i]->myDesc->lrec;
	for (j=0; j<inDesc->nrparm; j++) 
	  outUV[i]->buffer[jndx+j] =  inUV->buffer[indx+j];
 
	/* Scale u,v,w for new reference frequency */
	outUV[i]->buffer[jndx+inDesc->ilocu] *= scale[i];
	outUV[i]->buffer[jndx+inDesc->ilocv] *= scale[i];
	outUV[i]->buffer[jndx+inDesc->ilocw] *= scale[i];
     }

      /* Copy visibility */
      indx += inDesc->nrparm;
      for (i=0; i<nOut; i++) {
	jndx = ivis*outUV[i]->myDesc->lrec + outUV[i]->myDesc->nrparm;
	FreqSel (inUV->myDesc, outUV[i]->myDesc, 
		 BChan[i], BChan[i]+numChan[i]-1, chinc, BIF[i], BIF[i]+nIFperOut-1,
		 &inUV->buffer[indx], &outUV[i]->buffer[jndx]);
      }
    } /* end loop over visibilities */

    /* Write outputs */
   for (i=0; i<nOut; i++) {
     oretCode = ObitUVWrite (outUV[i], outUV[i]->buffer, err);
     if (err->error) goto cleanup;
   }
  } /* end loop processing data */
  
  /* check for errors */
  if ((iretCode > OBIT_IO_EOF) || (oretCode > OBIT_IO_EOF) ||
      (err->error)) goto cleanup;

  /* Cleanup */
 cleanup:
  if (scale)   g_free(scale);
  if (BChan)   g_free(BChan);
  if (numChan) g_free(numChan);
  if (BIF)     g_free(BIF);
  if (numIF)   g_free(numIF);
 
  /* close files */
  iretCode = ObitUVClose (inUV, err);
  for (i=0; i<nOut; i++) {
    oretCode = ObitUVClose (outUV[i], err);
    if ((iretCode!=OBIT_IO_OK) || (oretCode!=OBIT_IO_OK) || (err->error))
      Obit_traceback_msg (err, routine, inUV->name);
  }
  if (err->error) Obit_traceback_msg (err, routine, inUV->name);
  
  return;
} /* end ObitUVUtilSplitCh */

/**
 * Add Gaussian noise to an UV
 * out = in*scale +  noise (sigma), real, imag, each vis
 * Note: This uses the GSL random number generator, if this is not available
 * then only the scaling is done.
 * \param inUV  Input UV 
 * \param outUV Output UV, must already be defined but may be inUV
 * \param scale  scaling factor for data
 * \param sigma  Standard deviation of Gaussian noise .
 * \param err    Error stack
 */
void ObitUVUtilNoise(ObitUV *inUV, ObitUV *outUV, ofloat scale, ofloat sigma, 
		     ObitErr *err)
{
  ObitIOCode retCode;
  gboolean doCalSelect, done, same;
  ObitInfoType type;
  ObitIOAccess access, oaccess;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ofloat val, fblank = ObitMagicF();
  odouble dsigma = sigma;
  olong i, j, indx, NPIO, firstVis;
  ObitUVDesc *inDesc, *outDesc;
  /* Don't copy Cal and Soln or data or flag tables */
  gchar *exclude[]={"AIPS CL","AIPS SN","AIPS FG","AIPS CQ","AIPS WX",
		    "AIPS AT","AIPS CT","AIPS OB","AIPS IM","AIPS MC",
		    "AIPS PC","AIPS NX","AIPS TY","AIPS GC","AIPS HI",
		    "AIPS PL","AIPS NI","AIPS SY","AIPS PT","AIPS OT",
		    NULL};
#if HAVE_GSL==1  /* GSL stuff */
  gsl_rng *ran=NULL;
#endif /* HAVE_GSL */
  gchar *routine = "ObitUVUtilNoise";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitUVIsA(inUV));
  g_assert (ObitUVIsA(outUV));

  /* Are input and output the same file? */
  same = ObitUVSame(inUV, outUV, err);
  if (err->error) Obit_traceback_msg (err, routine, inUV->name);

  /* Local pointers */
  inDesc  = inUV->myDesc;
  outDesc = outUV->myDesc;

  /* Calibration wanted? */ 
  doCalSelect = FALSE;
  ObitInfoListGetTest(inUV->info, "doCalSelect", &type, (gint32*)dim, &doCalSelect);
  if (doCalSelect) access = OBIT_IO_ReadCal;
  else access = OBIT_IO_ReadWrite;

 /* Open Input Data */
  retCode = ObitUVOpen (inUV, access, err);
  if ((retCode != OBIT_IO_OK) || (err->error>0)) 
    Obit_traceback_msg (err, routine, inUV->name);

  /* use same data buffer on input and output.
     If multiple passes are made the input files will be closed
     which deallocates the buffer, use output buffer.
     so free input buffer */
  if (!same) {
    /* use same data buffer on input 1 and output 
       so don't assign buffer for output */
    if (outUV->buffer) ObitIOFreeBuffer(outUV->buffer); /* free existing */
    outUV->buffer = NULL;
    outUV->bufferSize = -1;
  }

   /* Copy number of records per IO to output */
  ObitInfoListGet (inUV->info, "nVisPIO", &type, dim,   (gpointer)&NPIO, err);
  ObitInfoListPut (outUV->info, "nVisPIO",  type, dim,  (gpointer)&NPIO, err);
  if (err->error) Obit_traceback_msg (err, routine, inUV->name);

  /* Open Output Data */
  if (same) oaccess = OBIT_IO_ReadWrite;
  else      oaccess = OBIT_IO_WriteOnly;
  retCode = ObitUVOpen (outUV, oaccess, err) ;
  if ((retCode != OBIT_IO_OK) || (err->error>0)) {
    outUV->buffer = NULL; /* remove pointer to inUV buffer */
    outUV->bufferSize = 0;
    Obit_traceback_msg (err, routine, outUV->name);
  }

  /* Copy tables before data */
  if (!same) {
    retCode = ObitUVCopyTables (inUV, outUV, exclude, NULL, err);
    if (err->error) {/* add traceback,return */
      outUV->buffer = NULL;
      outUV->bufferSize = 0;
      Obit_traceback_msg (err, routine, inUV->name);
    }
    
    /* Close and reopen input to init calibration which will have been disturbed 
       by the table copy */
    retCode = ObitUVClose (inUV, err);
    if (err->error) {
      outUV->buffer = NULL; outUV->bufferSize = 0;
      Obit_traceback_msg (err, routine, inUV->name);
    }
    
    retCode = ObitUVOpen (inUV, access, err);
    if ((retCode != OBIT_IO_OK) || (err->error>0)) {
      outUV->buffer = NULL; outUV->bufferSize = 0;
      Obit_traceback_msg (err, routine, inUV->name);
    }
  } /* end if not same */
  outUV->buffer = inUV->buffer;

  /* Init random number generator */
#if HAVE_GSL==1  /* GSL stuff */
  ran = gsl_rng_alloc(gsl_rng_default);
#endif /* HAVE_GSL */
  
  /* Loop over data */
  done = (retCode != OBIT_IO_OK);
  while (!done) {
    
    /* read buffer */
    retCode = ObitUVRead (inUV, NULL, err);
    if (err->error) {
      outUV->buffer = NULL; outUV->bufferSize = 0;
      Obit_traceback_msg (err, routine, inUV->name);
    }
    done = (retCode == OBIT_IO_EOF); /* done? */
    if (done) break;

    /* How many? */
    outDesc->numVisBuff = inDesc->numVisBuff;

    /* Modify data */
    for (i=0; i<inDesc->numVisBuff; i++) { /* loop over visibilities */
      indx = i*inDesc->lrec + inDesc->nrparm;
      for (j=0; j<inDesc->ncorr; j++) { /* loop over correlations */
	if (inUV->buffer[indx]!=fblank) {
	  val = inUV->buffer[indx]*scale;
#if HAVE_GSL==1  /* GSL stuff */
	  val += (ofloat)gsl_ran_gaussian (ran, dsigma);
#endif /* HAVE_GSL */
	  inUV->buffer[indx]  = val;
	}
	if (inUV->buffer[indx+1]!=fblank) {
	  val = inUV->buffer[indx+1]*scale;
#if HAVE_GSL==1  /* GSL stuff */
	  val += (ofloat)gsl_ran_gaussian (ran, dsigma);
#endif /* HAVE_GSL */
	  inUV->buffer[indx+1]  = val;
	}
	indx += inDesc->inaxes[0];
      } /* end loop over correlations */
    } /* end loop over visibilities */

    
    /* Write buffer - if same as input, fiddle first vis value */
    firstVis = outDesc->firstVis;
    retCode = ObitUVWrite (outUV, NULL, err);
    if (same) {
      outDesc->firstVis = firstVis;
      ((ObitUVDesc*)(outUV->myIO->myDesc))->firstVis = firstVis;
    }
    if (err->error) {
      outUV->buffer = NULL; outUV->bufferSize = 0;
      Obit_traceback_msg (err, routine, outUV->name);
    }
  } /* end loop over data */
  
  /* unset output buffer (may be multiply deallocated ;'{ ) */
  outUV->buffer = NULL;
  outUV->bufferSize = 0;
  
  /* Free random number generator */
#if HAVE_GSL==1  /* GSL stuff */
  gsl_rng_free(ran);
#endif /* HAVE_GSL */
  
  /* Close input */
  retCode = ObitUVClose (inUV, err);
  if (err->error) Obit_traceback_msg (err, routine, inUV->name);
  
  /* Close output */
  retCode = ObitUVClose (outUV, err);
  if (err->error) Obit_traceback_msg (err, routine, outUV->name);

} /* end ObitUVUtilNoise */

/**
 * Adds flagging entry to associated flag table
 * Input values on inUV
 * \li "flagVer"   OBIT_Int (1,1,1) Flagging table version, default = 1
 * \li "subA"      OBIT_Int (1,1,1) Subarray, default = 0
 * \li "freqID"    OBIT_Int (1,1,1) Frequency ID, default = 0
 * \li "timeRange" OBIT_float (2,1,1) Start and stop times to flag (days) def, 0s=all
 * \li "Chans"     OBIT_Int (2,1,1) First and highest channels to flag (1-rel), def, 0=>all
 * \li "IFs"       OBIT_Int (2,1,1) First and highest IF to flag (1-rel), def, 0=>all
 * \li "Ants"      OBIT_Int (2,1,1) first and second antenna  numbers for a baseline, 0=$>$all
 * \li "Source"    OBIT_string (?,1,1) Name of source, def, "Any" => all
 * \li "Stokes"    OBIT_string (?,1,1) Stokes to flag, def " " = flag all
 *                 "FFFF"  where F is '1' to flag corresponding Stokes, '0' not.
 *                 Stokes order 'R', 'L', 'RL' 'LR' or 'X', 'Y', 'XY', 'YX'
 * \li "Reason"    OBIT_string (?,1,1) reason string for flagging (max. 24 char).
 * \param inUV   Input UV data
 * \param err     Error stack, returns if not empty.
 * \return IO return code, OBIT_IO_OK = OK
 */
ObitIOCode ObitUVUtilFlag (ObitUV *inUV, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  oint iarr[2];
  olong flagVer, subA, freqID, chans[2], ifs[2], ants[2], SouID;
  ofloat timerange[2];
  gchar source[49], stokes[10], reason[49];
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  ObitTableFG    *FlagTable=NULL;
  ObitTableFGRow *FlagRow=NULL;
  ObitTableSU    *SourceTable=NULL;
  gchar *tname, souCode[5];
  olong ver, iRow, Qual, Number;
  gboolean xselect;
  oint numIF;
  gchar *routine = "ObitUVUtilFlag";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitUVIsA(inUV));

  /* Get parameters */
  flagVer = 1;
  ObitInfoListGetTest(inUV->info, "flagVer", &type, dim, &flagVer);

  subA = 0;
  ObitInfoListGetTest(inUV->info, "subA", &type, dim, &subA);

  freqID = 0;
  ObitInfoListGetTest(inUV->info, "freqID", &type, dim, &freqID);

  chans[0] = 1; chans[0] = 0;
  ObitInfoListGetTest(inUV->info, "Chans", &type, dim, chans);

  ifs[0] = 1; ifs[0] = 0;
  ObitInfoListGetTest(inUV->info, "IFs", &type, dim, ifs);

  ants[0] = 1; ants[0] = 0;
  ObitInfoListGetTest(inUV->info, "Ants", &type, dim, ants);

  timerange[0] = -1.0e20;  timerange[1] = 1.0e20;
  ObitInfoListGetTest(inUV->info, "timeRange", &type, dim, timerange);
  /* default */
  if ((timerange[0]==0.0) && (timerange[1]==0.0)) {
    timerange[0] = -1.0e20;
    timerange[1] =  1.0e20;
  }

  g_snprintf (source, 48, "Any");
  ObitInfoListGetTest(inUV->info, "Source", &type, dim, source);
  source[dim[0]] = 0;   /* terminate */

  g_snprintf (stokes, 9, " ");
  ObitInfoListGetTest(inUV->info, "Stokes", &type, dim, stokes);
  stokes[dim[0]] = 0;   /* terminate */

  g_snprintf (reason, 48, " ");
  ObitInfoListGetTest(inUV->info, "Reason", &type, dim, reason);
  reason[dim[0]] = 0;   /* terminate */

  /* Open/close input UV to fully instantiate */
  retCode = ObitUVOpen (inUV, OBIT_IO_ReadWrite, err);
  if (err->error) Obit_traceback_val (err, routine, inUV->name, retCode);
  
  /* Close */
  retCode = ObitUVClose (inUV, err);
  if (err->error) Obit_traceback_val (err, routine, inUV->name, retCode);

  /* Look up Source number if needed */
  if (strncmp ("Any", source, 3)) {
    /* Instantiate/Create Source Table */
    retCode = OBIT_IO_ReadErr;
    tname = g_strconcat ("SU table for: ",inUV->name, NULL);
    ver = 1;
    if (inUV->myDesc->jlocif>=0) numIF = inUV->myDesc->inaxes[inUV->myDesc->jlocif];
    else numIF = 1;
    SourceTable = newObitTableSUValue(tname, (ObitData*)inUV, &ver, 
				      OBIT_IO_ReadWrite, numIF, err);
    dim[0] = strlen(source);
    dim[1] = 1;
    Number = 1;
    Qual   = -1;
    sprintf (souCode, "    ");
    ObitTableSULookup (SourceTable, dim, source, Qual, souCode, iarr, 
		       &xselect, &Number, err);
    if (err->error) Obit_traceback_val (err, routine, inUV->name, retCode);
    SouID = iarr[0];  /* Source ID */
    SourceTable = ObitTableSUUnref(SourceTable);
    g_free (tname);
  } else { /* flag all sources */
    SouID = 0;
  }

  /* Instantiate/Create output Flag Table */
  tname = g_strconcat ("FG table for: ", inUV->name, NULL);
  ver = flagVer;
  FlagTable = newObitTableFGValue(tname, (ObitData*)inUV, &ver, 
				       OBIT_IO_ReadWrite, err);
  if (err->error) Obit_traceback_val (err, routine, inUV->name, retCode);
  g_free (tname);

  /* Open table */
  retCode = ObitTableFGOpen (FlagTable, OBIT_IO_ReadWrite, err);
  if (err->error) Obit_traceback_val (err, routine, inUV->name, retCode);

  /* Create Table Row */
  FlagRow = newObitTableFGRow (FlagTable);
  
  /* Attach  row to output buffer */
  ObitTableFGSetRow (FlagTable, FlagRow, err);
  if (err->error) Obit_traceback_val (err, routine, inUV->name, retCode);

  /* If there are entries in the table, mark it unsorted */
  if (FlagTable->myDesc->nrow>0) 
    {FlagTable->myDesc->sort[0]=0; FlagTable->myDesc->sort[1]=0;}
  
  /* Fill in Flag row */
  FlagRow->SourID = SouID;
  FlagRow->SubA   = subA;
  FlagRow->freqID = freqID;
  FlagRow->TimeRange[0] = timerange[0];
  FlagRow->TimeRange[1] = timerange[1];
  FlagRow->ants[0]  = ants[0];
  FlagRow->ants[1]  = ants[1];
  FlagRow->chans[0] = chans[0];
  FlagRow->chans[1] = chans[1];
  FlagRow->ifs[0]   = ifs[0];
  FlagRow->ifs[1]   = ifs[1];
  FlagRow->pFlags[0] = 0;
  strncpy (FlagRow->reason, reason, 24);
  if (stokes[0]==' ') {
    FlagRow->pFlags[0] = 15;
  } else {
    if (stokes[0]!='0') FlagRow->pFlags[0] += 1;
    if (stokes[1]!='0') FlagRow->pFlags[0] += 2;
    if (stokes[2]!='0') FlagRow->pFlags[0] += 4;
    if (stokes[3]!='0') FlagRow->pFlags[0] += 8;
  }
  
  /* write row */
  iRow = FlagTable->myDesc->nrow+1;
  retCode =  ObitTableFGWriteRow (FlagTable, iRow, FlagRow, err);
  if (err->error) Obit_traceback_val (err, routine, inUV->name, retCode);
   
  /* Close Flag table */
  retCode =  ObitTableFGClose (FlagTable, err);
  if (err->error) Obit_traceback_val (err, routine, inUV->name, retCode);

  /* Cleanup */
  FlagTable = ObitTableFGUnref(FlagTable);
  FlagRow   = ObitTableFGRowUnref(FlagRow);

  return retCode;
} /* end ObitUVUtilFlag */

/**
 * Make copy of ObitUV with the u,v,w terms calculated
 * \param inUV     Input uv data to copy. 
 * \param outUV    If not scratch, then the previously defined output file
 *                 May be NULL for scratch only
 *                 If it exists and scratch, it will be Unrefed
 * \param err      Error stack, returns if not empty.
 */
void ObitUVUtilCalcUVW (ObitUV *inUV, ObitUV *outUV,  ObitErr *err)
{
  ObitIOCode iretCode, oretCode;
  gboolean doCalSelect;
  gchar *exclude[]={"AIPS CL","AIPS SN","AIPS FG","AIPS CQ","AIPS WX",
		    "AIPS AT","AIPS CT","AIPS OB","AIPS IM","AIPS MC",
		    "AIPS PC","AIPS NX","AIPS TY","AIPS GC","AIPS HI",
		    "AIPS PL","AIPS NI","AIPS OT",
		    NULL};
  gchar *sourceInclude[] = {"AIPS SU", NULL};
  ObitUVWCalc *uvwCalc=NULL;
  olong i, indx, SId, subA, ant1, ant2;
  ofloat uvw[3];
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  ObitIOAccess access;
  ObitUVDesc *inDesc, *outDesc;
  gchar *today=NULL;
  gchar *routine = "ObitUVUtilCopyZero";
 
  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitUVIsA(inUV));
  if (outUV==NULL) {
    Obit_log_error(err, OBIT_Error,"%s Output MUST be defined",
		   routine);
      return;
  }

  /* Selection of input? */
  doCalSelect = FALSE;
  ObitInfoListGetTest(inUV->info, "doCalSelect", &type, (gint32*)dim, &doCalSelect);
  if (doCalSelect) access = OBIT_IO_ReadCal;
  else access = OBIT_IO_ReadOnly;

  /* test open to fully instantiate input and see if it's OK */
  iretCode = ObitUVOpen (inUV, access, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_msg (err, routine, inUV->name);

  /* copy Descriptor */
  outUV->myDesc = ObitUVDescCopy(inUV->myDesc, outUV->myDesc, err);

  /* Creation date today */
  today = ObitToday();
  strncpy (outUV->myDesc->date, today, UVLEN_VALUE-1);
  if (today) g_free(today);
  
  /* use same data buffer on input and output 
     so don't assign buffer for output */
  if (outUV->buffer) ObitIOFreeBuffer(outUV->buffer); /* free existing */
  outUV->buffer = NULL;
  outUV->bufferSize = -1;

  /* test open output */
  oretCode = ObitUVOpen (outUV, OBIT_IO_WriteOnly, err);
  /* If this didn't work try OBIT_IO_ReadWrite */
  if ((oretCode!=OBIT_IO_OK) || (err->error)) {
    ObitErrClear(err);
    oretCode = ObitUVOpen (outUV, OBIT_IO_ReadWrite, err);
  }
  /* if it didn't work bail out */
  if ((oretCode!=OBIT_IO_OK) || (err->error)) {
    /* unset output buffer (may be multiply deallocated) */
    outUV->buffer = NULL;
    outUV->bufferSize = 0;
    Obit_traceback_msg (err, routine, outUV->name);
  }

  /* iretCode = ObitUVClose (inUV, err); DEBUG */
  /* Copy tables before data */
  iretCode = ObitUVCopyTables (inUV, outUV, exclude, NULL, err);
  /* If multisource out then copy SU table, multiple sources selected or
   sources deselected suggest MS out */
  if ((inUV->mySel->numberSourcesList>1) || (!inUV->mySel->selectSources))
  iretCode = ObitUVCopyTables (inUV, outUV, NULL, sourceInclude, err);
  if (err->error) {
    outUV->buffer = NULL;
    outUV->bufferSize = 0;
    Obit_traceback_msg (err, routine, inUV->name);
  }

  /* reset to beginning of uv data */
  iretCode = ObitIOSet (inUV->myIO,  inUV->info, err);
  oretCode = ObitIOSet (outUV->myIO, outUV->info, err);
  if (err->error) Obit_traceback_msg (err, routine,inUV->name);

  /* Close and reopen input to init calibration which will have been disturbed 
     by the table copy */
  iretCode = ObitUVClose (inUV, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error))
    Obit_traceback_msg (err, routine,inUV->name);

  iretCode = ObitUVOpen (inUV, access, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) 
    Obit_traceback_msg (err, routine,inUV->name);
  outUV->buffer = inUV->buffer;

  /* Get descriptors */
  inDesc  = inUV->myDesc;
  outDesc = outUV->myDesc;
  SId = 0;   /* In case single source */

  uvwCalc = ObitUVWCalcCreate("UVWCalc", outUV, err);
  if (err->error) Obit_traceback_msg (err, routine, outUV->name);

  /* we're in business, copy, recompute u,v,w */
  while ((iretCode==OBIT_IO_OK) && (oretCode==OBIT_IO_OK)) {
    if (doCalSelect) iretCode = ObitUVReadSelect (inUV, inUV->buffer, err);
    else iretCode = ObitUVRead (inUV, inUV->buffer, err);
    if (iretCode!=OBIT_IO_OK) break;
   /* How many */
    outDesc->numVisBuff = inDesc->numVisBuff;

    /* Modify data */
    for (i=0; i<inDesc->numVisBuff; i++) { /* loop over visibilities */
      indx = i*inDesc->lrec;
      if (inUV->myDesc->ilocsu>=0) 
	SId = (olong)(inUV->buffer[indx+inUV->myDesc->ilocsu] + 0.5);
      ObitUVDescGetAnts(inUV->myDesc, &inUV->buffer[indx], &ant1, &ant2, &subA);
      ObitUVWCalcUVW(uvwCalc, inUV->buffer[indx+inDesc->iloct], SId, 
		     subA, ant1, ant2, uvw, err);
      inUV->buffer[indx+inDesc->ilocu] = uvw[0];
      inUV->buffer[indx+inDesc->ilocv] = uvw[1];
      inUV->buffer[indx+inDesc->ilocw] = uvw[2];
    } /* end loop over visibilities */

    /* Write */
    oretCode = ObitUVWrite (outUV, inUV->buffer, err);
    if (err->error) {
      uvwCalc = ObitUVWCalcUnref(uvwCalc);
      Obit_traceback_msg (err, routine,inUV->name);
    }
  } /* end loop processing data */
  
  /* check for errors */
  if ((iretCode > OBIT_IO_EOF) || (oretCode > OBIT_IO_EOF) ||
      (err->error)) /* add traceback,return */
    Obit_traceback_msg (err, routine,inUV->name);
  
  /* unset input buffer (may be multiply deallocated ;'{ ) */
  outUV->buffer = NULL;
  outUV->bufferSize = 0;
  
  uvwCalc = ObitUVWCalcUnref(uvwCalc);  /* Cleanup */

  /* close files */
  iretCode = ObitUVClose (inUV, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) 
    Obit_traceback_msg (err, routine, inUV->name);
  
  oretCode = ObitUVClose (outUV, err);
  if ((oretCode!=OBIT_IO_OK) || (err->error))
    Obit_traceback_msg (err, routine, outUV->name);
  
  return;
} /* end ObitUVUtilCalcUVW */

/**
 * Compute low precision  visibility uvw
 * \param b    Baseline vector
 * \param dec  Pointing declination (rad)
 * \param ha   Pointing hour angle (rad)
 * \param uvw  [out] baseline u,v,w in units of b
 */
void ObitUVUtilUVW(const ofloat b[3], odouble dec, ofloat ha, ofloat uvw[3])
{
  odouble cosdec, sindec;
  ofloat     sinha, cosha, vw;

  cosdec = cos (dec);
  sindec = sin (dec);
  cosha = cos (ha);
  sinha = sin (ha);
  vw = b[0]*cosha - b[1]*sinha;
  uvw[0] =  b[0]*sinha + b[1]*cosha;
  uvw[1] = -vw*sindec + b[2]*cosdec;
  uvw[2] =  vw*cosdec + b[2]*sindec;
} /* end  ObitUVUtilUVW */

/**
 * Append the contents of one UV onto the end of another
 * \param inUV  Input UV 
 * \param outUV Output UV, must already be defined
 * \param err    Error stack
 */
void ObitUVUtilAppend(ObitUV *inUV, ObitUV *outUV, ObitErr *err)
{
  ObitIOCode retCode;
  gboolean doCalSelect, done;
  ObitInfoType type;
  ObitIOAccess access;
  gboolean incompatible;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitUVDesc *inDesc, *outDesc;
  olong inNPIO, outNPIO, NPIO;
  gchar *routine = "ObitUVUtilAppend";

  /* error checks */
  if (err->error) return;
  g_assert (ObitUVIsA(inUV));
  g_assert (ObitUVIsA(outUV));

  /* Get input descriptors */
  inDesc  = inUV->myDesc;
  outDesc = outUV->myDesc;

  /* Check compatability between inUV, outUV */
  incompatible = (inDesc->ncorr!=outDesc->ncorr);
  incompatible = incompatible || (inDesc->jlocs!=outDesc->jlocs);
  incompatible = incompatible || (inDesc->jlocf!=outDesc->jlocf);
  incompatible = incompatible || (inDesc->jlocif!=outDesc->jlocif);
  incompatible = incompatible || (inDesc->ilocb!=outDesc->ilocb);
  if (incompatible) {
     Obit_log_error(err, OBIT_Error,"%s inUV and outUV have incompatible structures",
		   routine);
      return ;
 }
  /* Calibration wanted? */ 
  doCalSelect = FALSE;
  ObitInfoListGetTest(inUV->info, "doCalSelect", &type, dim, &doCalSelect);
  if (doCalSelect) access = OBIT_IO_ReadCal;
  else access = OBIT_IO_ReadWrite;

  /* Set number of vis per I/O */
  inNPIO = 1000;
  ObitInfoListGetTest (inUV->info, "nVisPIO", &type, dim, &inNPIO);
  outNPIO = 1000;
  ObitInfoListGetTest (outUV->info, "nVisPIO", &type, dim, &outNPIO);
  NPIO = 1000; dim[0] = dim[1] = dim[2] = 1;
  ObitInfoListAlwaysPut (inUV->info,  "nVisPIO", OBIT_long, dim,  &NPIO);
  ObitInfoListAlwaysPut (outUV->info, "nVisPIO", OBIT_long, dim,  &NPIO);

  /* Open Input Data */
  retCode = ObitUVOpen (inUV, access, err);
  if ((retCode != OBIT_IO_OK) || (err->error>0)) 
    Obit_traceback_msg (err, routine, inUV->name);
  
  /* use same data buffer on input and output 
     so don't assign buffer for output */
  if (outUV->buffer) ObitIOFreeBuffer(outUV->buffer); /* free existing */
  outUV->buffer     = inUV->buffer;
  outUV->bufferSize = inUV->bufferSize;

  /* Open Output Data */
  retCode = ObitUVOpen (outUV, OBIT_IO_ReadWrite, err) ;
  if ((retCode != OBIT_IO_OK) || (err->error>0)) {
    outUV->buffer = NULL; /* remove pointer to inUV buffer */
    outUV->bufferSize = 0;
    Obit_traceback_msg (err, routine, outUV->name);
  }
  outDesc->firstVis = outDesc->nvis+1; /* Write to end */

  /* Loop over data */
  done = (retCode != OBIT_IO_OK);
  while (!done) {
    
    /* read buffer */
    retCode = ObitUVRead (inUV, NULL, err);
    if (err->error) {
      outUV->buffer = NULL; outUV->bufferSize = 0;
      Obit_traceback_msg (err, routine, inUV->name);
    }
    done = (retCode == OBIT_IO_EOF); /* done? */
    if (done) break;

    /* How many? */
    outDesc->numVisBuff = inDesc->numVisBuff;
   
    /* Write buffer */
    retCode = ObitUVWrite (outUV, NULL, err);
    if (err->error) {
      outUV->buffer = NULL; outUV->bufferSize = 0;
      Obit_traceback_msg (err, routine, outUV->name);
    }
  } /* end loop over data */
  
  /* unset output buffer (may be multiply deallocated ;'{ ) */
  outUV->buffer = NULL;
  outUV->bufferSize = 0;
  
  /* Close input */
  retCode = ObitUVClose (inUV, err);
  if (err->error) Obit_traceback_msg (err, routine, inUV->name);
  
  /* Close output */
  retCode = ObitUVClose (outUV, err);
  if (err->error) Obit_traceback_msg (err, routine, outUV->name);

  /* Reset number of vis per I/O */
  ObitInfoListAlwaysPut (inUV->info,  "nVisPIO", OBIT_long, dim,  &inNPIO);
  ObitInfoListAlwaysPut (outUV->info, "nVisPIO", OBIT_long, dim,  &outNPIO);

} /* end ObitUVUtilAppend */

#ifndef VELIGHT
#define VELIGHT 2.997924562e8
#endif /* VELIGHT */
/**
 *  How many channels can I average
 * \param inUV    Input UV 
 * \param maxFact Maximum allowed bandwith smearing amplitude loss
 * \param FOV     Radius of desired FOV (deg)
 * \param err     Error stack
 */
olong ObitUVUtilNchAvg(ObitUV *inUV, ofloat maxFact, ofloat FOV, ObitErr *err)
{
  olong out=1;
  ObitIOCode iretCode;
  ObitUVDesc *inDesc;
  ObitTableAN *ANTable=NULL;
  ObitAntennaList **AntList=NULL;
  olong i, j, numSubA, iANver, numIF, numOrb, numPCal;
  ofloat chBW, maxBL, BL, fact, beta, tau;
  gchar *routine = "ObitUVUtilNchAv";

  /* error checks */
  if (err->error) return out;
  g_assert (ObitUVIsA(inUV));

  /* test open to fully instantiate input and see if it's OK */
  ObitUVFullInstantiate (inUV, TRUE, err);
  iretCode = ObitUVOpen (inUV, OBIT_IO_ReadCal, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, inUV->name, out);

  /* Get descriptor */
  inDesc = inUV->myDesc;

  /* Channel bandwidth */
  chBW = inDesc->cdelt[inDesc->jlocf];

  /* Antenna List 
     How many AN tables (no. subarrays)?  */
  numSubA = ObitTableListGetHigh (inUV->tableList, "AIPS AN");
  AntList = g_malloc0(numSubA*sizeof(ObitAntennaList*));
  maxBL = 0.0;

  /* Loop over AN tables (subarrays) */
  for (iANver=1; iANver<=numSubA; iANver++) {
    numOrb   = 0;
    numPCal  = 0;
    numIF    = 0;

    ANTable = newObitTableANValue ("AN table", (ObitData*)inUV, 
				   &iANver, OBIT_IO_ReadOnly, numIF, numOrb, numPCal, err);
    if (ANTable==NULL) Obit_log_error(err, OBIT_Error, "ERROR with AN table");
    AntList[iANver-1] = ObitTableANGetList (ANTable, err);
    if (err->error) Obit_traceback_val (err, routine, inUV->name, out);
    
    /* Cleanup */
    ANTable = ObitTableANUnref(ANTable);

    /* Find maximum baseline */
    for (i=0; i<AntList[iANver-1]->number-1; i++) {
      /* Antenna in array? */
      if ((fabs(AntList[iANver-1]->ANlist[i]->AntXYZ[0]<1.0)) &&
	   (fabs(AntList[iANver-1]->ANlist[i]->AntXYZ[1]<1.0)) &&
	  (fabs(AntList[iANver-1]->ANlist[i]->AntXYZ[2]<1.0))) continue;
      for (j=i+1; j<AntList[iANver-1]->number; j++) {
	/* Antenna in array? */
	if ((fabs(AntList[iANver-1]->ANlist[j]->AntXYZ[0]<1.0)) &&
	    (fabs(AntList[iANver-1]->ANlist[j]->AntXYZ[1]<1.0)) &&
	    (fabs(AntList[iANver-1]->ANlist[j]->AntXYZ[2]<1.0))) continue;
	BL = 
	  (AntList[iANver-1]->ANlist[i]->AntXYZ[0]-AntList[iANver-1]->ANlist[j]->AntXYZ[0]) *
	  (AntList[iANver-1]->ANlist[i]->AntXYZ[0]-AntList[iANver-1]->ANlist[j]->AntXYZ[0]) + 
	  (AntList[iANver-1]->ANlist[i]->AntXYZ[1]-AntList[iANver-1]->ANlist[j]->AntXYZ[1]) *
	  (AntList[iANver-1]->ANlist[i]->AntXYZ[1]-AntList[iANver-1]->ANlist[j]->AntXYZ[1]) + 
	  (AntList[iANver-1]->ANlist[i]->AntXYZ[2]-AntList[iANver-1]->ANlist[j]->AntXYZ[2]) *
	  (AntList[iANver-1]->ANlist[i]->AntXYZ[2]-AntList[iANver-1]->ANlist[j]->AntXYZ[2]);
	maxBL = MAX (maxBL, BL);
      }
    }
  } /* End loop over subarrays */
  maxBL = sqrt(maxBL);
  /* Cleanup */
  if (AntList) {
    for (i=0; i<numSubA; i++) {
      AntList[i] = ObitAntennaListUnref(AntList[i]);
    }
    g_free(AntList);
  }

  /* Close up */
  iretCode = ObitUVClose (inUV, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) 
    Obit_traceback_val (err, routine, inUV->name, out);

  /* Calculate number of channels to average */
  tau = sin(FOV*DG2RAD) * maxBL / VELIGHT;  /* Maximum delay across FOV */
  for (i=2; i<=inDesc->inaxes[inDesc->jlocf]; i++) {
    beta = i * chBW;
    fact = (G_PI*beta*tau) / fabs(sin(G_PI*beta*tau));
    if (fact>maxFact) break;
    out = i;
  }
  
  return out;
} /* end ObitUVUtilNchAvg */

/*----------------------Private functions---------------------------*/
/**
 * Count channels after selection and averaging and modify descriptor.
 * Also determines arrays giving the channel, IF and stokes of each correlator
 * \param inDesc   Input UV descriptor
 * \param outDesc  Output UV descriptor to be modified
 * \param NumChAvg Number of channels to average
 * \param ChanSel  Groups of channels/IF in input to be averaged together
 *                 (start, end, increment, IF) where start and end at the 
 *                 beginning and ending channel numbers (1-rel) of the group
 *                 to be averaged together, increment is the increment between
 *                 selected channels and IF is the IF number (1-rel)
 *                 ChanSel is ignored if NumChAvg is given and > 0.
 *                 default increment is 1, IF=0 means all IF.
 *                 The list of groups is terminated by a start <=0
 * \param doAvgAll If true all channels and IFs to be averaged together
 * \param corChan  [out] 0-rel output channel numbers
 *                 Should be externally allocated to at least the number of correlators
 * \param corIF    [out] 0-rel output IF numbers
 *                 Should be externally allocated to at least the number of correlators
 * \param corStok  [out] 0-rel output Stokes parameter code.
 *                 Should be externally allocated to at least the number of correlators
 * \param corMask  [out] Array of masks per correlator for selected channels/IFs
 *                 Should be externally allocated to at least the number of correlators
 * \param err      Error stack, returns if not empty.
 * \return scaling factor for U,V,W to new reference frequency
 */
static ofloat AvgFSetDesc (ObitUVDesc *inDesc, ObitUVDesc *outDesc, 
			   olong NumChAvg, olong *ChanSel, gboolean doAvgAll, 
			   olong *corChan, olong *corIF, olong *corStok, gboolean *corMask,
			   ObitErr *err)
{
  ofloat scale = 1.0;
  olong i, ii, count, count2, NumIFAvg, numChan = 1, numIF = 1;
  olong ichan, iif, istok, nchan, nif, nstok;
  olong incs, incf, incif, ioff, lfoff, soff;
  olong *selTemp;
  olong *tcorChan=NULL, *tcorIF=NULL;
  gboolean more, match;
  odouble oldFreq, sum, sum2;

  /* error checks */
  if (err->error) return scale;
  g_assert(ObitUVDescIsA(inDesc));
  g_assert(ObitUVDescIsA(outDesc));

  /* Set up for parsing data */
  nchan = inDesc->inaxes[inDesc->jlocf];
  if (inDesc->jlocif>=0) nif = inDesc->inaxes[inDesc->jlocif];
  else nif = 1;
  if (inDesc->jlocs>=0) nstok = inDesc->inaxes[inDesc->jlocs];
  else nstok = 1;

  numChan = 1 + (nchan-1) / MAX (1, NumChAvg);  /* Number of output channels */
  /* IF averaging */
  NumIFAvg   = 1;
  if (doAvgAll) NumIFAvg   = nif;

  /* Enforce ChanSel Limits */
  selTemp = ChanSel;
  more = selTemp[1]>0;
  while (more) {
    if (selTemp[0]<1)       selTemp[0] = 1;
    if (selTemp[1]>nchan)   selTemp[1] = nchan;
    if (selTemp[2]<1)       selTemp[2] = 1;
    if (selTemp[3]<0)       selTemp[3] = 0;
    if (selTemp[3]>nif)     selTemp[3] = nif;
    selTemp += 4;
    more = selTemp[1]>0;
  } /* end loop enforcing limits */

  /* get increments (one word per correlator) */
  incs  = inDesc->incs  / inDesc->inaxes[0];
  incf  = inDesc->incf  / inDesc->inaxes[0];
  incif = inDesc->incif / inDesc->inaxes[0];

  /* temp arrays */
  tcorChan = g_malloc0(inDesc->ncorr*sizeof(olong));
  tcorIF   = g_malloc0(inDesc->ncorr*sizeof(olong));

  /* loop over IF */
  for (iif=0; iif<nif; iif++) {
    lfoff = iif * incif;
    ioff  = lfoff;
    
    /* Loop over frequency channel */
    for (ichan=0; ichan<nchan; ichan++) { /* loop 60 */
      soff = ioff;

      /* Loop over polarization */
      for (istok=0; istok<nstok; istok++) {
	corChan[soff]  = (ichan/NumChAvg); 
	corIF[soff]    = (iif/NumIFAvg); 
	corStok[soff]  = istok;
	corMask[soff]  = FALSE;
	tcorChan[soff] = (ichan); 
	tcorIF[soff]   = (iif); 

	/* Is this correlator/IF in ChanSel? */
	selTemp = ChanSel;
	more = selTemp[1]>0;
	while (more) {
	  /* In IF range? */
	  if ((selTemp[3]<=0) || (selTemp[3]==(iif+1))) {
	    if ((selTemp[0]<=(ichan+1)) && (selTemp[1]>=(ichan+1))) {
	      /* Desired channel? */
	      match = (selTemp[2]<=1) || (((ichan-selTemp[0])%(selTemp[2]))==0);
	      corMask[soff] = match;
	    } /* end channel range */
	  } /* end IF range */
	  selTemp += 4;
	  more = selTemp[1]>0;
	}
	soff += incs;
      } /* end loop over stokes */
      ioff  += incf;     
    } /* end loop over channel */
  } /* end loop over IF */

  /* Averaging all channels/IF? */
  if (doAvgAll) {
    outDesc->inaxes[outDesc->jlocf] = nchan;
    /* Average all frequencies */
    sum  = 0.0; count  = 0;
    sum2 = 0.0; count2 = 0;
    for (i=0; i<inDesc->ncorr; i++) {
      ii = tcorChan[i] + tcorIF[i]*nchan;
      sum2 += inDesc->freqArr[ii];
      count2++;
      if (corMask[i]) {
	sum += inDesc->freqArr[ii];
	count++;
      }
    } /* end loop over correlators */

    if (count>0)
      outDesc->crval[outDesc->jlocf] = sum / (ofloat)count;
    else if (count2>0)
      outDesc->crval[outDesc->jlocf] = sum2 / (ofloat)count2;
    outDesc->inaxes[outDesc->jlocf] = numChan;
    outDesc->crpix[outDesc->jlocf] = 1.0;
    outDesc->cdelt[outDesc->jlocf] = 
      inDesc->cdelt[inDesc->jlocf] * inDesc->inaxes[inDesc->jlocf];
    if (outDesc->jlocif>=0) {
      outDesc->inaxes[outDesc->jlocif] = numIF;
      outDesc->crval[outDesc->jlocif] = 1.0;
      outDesc->crpix[outDesc->jlocif] = 1.0;
      outDesc->cdelt[outDesc->jlocif] = 1.0;
    }
    /* Alternate frequency/vel */
    outDesc->altCrpix = 1.0;
    outDesc->altRef   = inDesc->altRef;

    return scale;
  } /* End doAvgAll */

  /* Averaging channels  - IFs unaffected */
  numChan = 1 + ((inDesc->inaxes[inDesc->jlocf]-1) / NumChAvg);
  /* Get frequency of first average */
  sum  = 0.0; count  = 0;
  sum2 = 0.0; count2 = 0;
  ii = -1;
  for (i=0; i<inDesc->ncorr; i++) {
    if ((corChan[i]>0) || (corStok[i]>0) || (corIF[i]>0)) continue;  /* want? */
    ii++;
    sum2 += inDesc->freqArr[ii];
    count2++;
    if (corMask[i]) {
      sum += inDesc->freqArr[ii];      count++;
    }
  } /* end loop over correlators */

  oldFreq = outDesc->crval[outDesc->jlocf]; /* Save old reference freq */
  outDesc->inaxes[outDesc->jlocf] = numChan;
  if (count>0)
    outDesc->crval[outDesc->jlocf]  = sum / (ofloat)count;
  else if (count2>0)
    outDesc->crval[outDesc->jlocf]  = sum2 / (ofloat)count2;
  outDesc->crpix[outDesc->jlocf]  = 1.0;
  outDesc->cdelt[outDesc->jlocf]  = inDesc->cdelt[inDesc->jlocf] * NumChAvg;

  /* Frequency scaling */
  scale = outDesc->crval[outDesc->jlocf]/oldFreq;

  /* Alternate frequency/vel */
  outDesc->altCrpix = 1.0 + (outDesc->altCrpix-1.0)/NumChAvg;
  outDesc->altRef   = inDesc->altRef;

  /* Cleanup */
  if (tcorChan) g_free(tcorChan);
  if (tcorIF) g_free(tcorIF);

  return scale;
  
} /* end AvgFSetDesc */

/**
 * Modify descriptor for duplicatiing channels
 * \param inDesc   Input UV descriptor
 * \param outDesc  Output UV descriptor to be modified, assumed mostly copied.
 * \param nBloat   Number of time to duplicate channels 
 * \param err      Error stack, returns if not empty.
 * \return scaling factor for U,V,W to new reference frequency
 */
static ofloat BloatFSetDesc (ObitUVDesc *inDesc, ObitUVDesc *outDesc, 
			     olong nBloat, ObitErr *err)
{
  olong nchan;
  ofloat scale = 1.0;

  /* error checks */
  if (err->error) return scale;
  g_assert(ObitUVDescIsA(inDesc));
  g_assert(ObitUVDescIsA(outDesc));

  /* Modify */
  nchan = (inDesc->inaxes[inDesc->jlocf]+nBloat-1) * nBloat;
  outDesc->inaxes[outDesc->jlocf] = nchan;
  outDesc->cdelt[outDesc->jlocf]  = inDesc->cdelt[inDesc->jlocf] / nBloat;
  outDesc->crval[outDesc->jlocf]  = inDesc->crval[inDesc->jlocf] - 
    0.5 * outDesc->cdelt[outDesc->jlocf];
  scale = outDesc->crval[outDesc->jlocf] / inDesc->crval[inDesc->jlocf];

  /* Alternate frequency/vel */
  outDesc->altCrpix = 1.0 + (outDesc->altCrpix-1.0)*nBloat;
  outDesc->altRef   = inDesc->altRef;

  return scale;
  
} /* end BloatFSetDesc */

/**
 * Average visibility in frequency
 * \param inDesc    Input UV descriptor
 * \param outDesc   Output UV descriptor to be modified
 * \param NumChAvg  Number of channels to average
 * \param ChanSel   Groups of channels/IF in input to be averaged together
 *                  (start, end, increment, IF) where start and end at the 
 *                  beginning and ending channel numbers (1-rel) of the group
 *                  to be averaged together, increment is the increment between
 *                  selected channels and IF is the IF number (1-rel)
 *                  ChanSel is ignored if NumChAvg is given and > 0.
 *                  default increment is 1, IF=0 means all IF.
 *                  The list of groups is terminated by a start <=0
 *                  corMask actually used to select.
 * \param doAvgAll  If true all channels and IFs to be averaged together
 * \param corChan  0-rel output channel numbers
 * \param corIF    0-rel output IF numbers
 * \param corStok  0-rel output Stokes parameter code.
 * \param corMask  Array of masks per correlator for selected channels/IFs
 *                 TRUE => select
 * \param inBuffer  Input buffer (data matrix)
 * \param outBuffer Output buffer (data matrix)
 * \param work      Work array twice the size of the output visibility
 * \param err       Error stack, returns if not empty.
 */
static void AvgFAver (ObitUVDesc *inDesc, ObitUVDesc *outDesc, 
		      olong NumChAvg, olong *ChanSel, gboolean doAvgAll, 
		      olong *corChan, olong *corIF, olong *corStok, gboolean *corMask,
		      ofloat *inBuffer, ofloat *outBuffer, ofloat *work, 
		      ObitErr *err)
{
  olong i, n, indx, jndx, jf, jif, js, nochan, noif;
  olong nchan, nif, nstok, incs, incf, incif;

  /* error checks */
  if (err->error) return;

  nchan = inDesc->inaxes[inDesc->jlocf];
  if (inDesc->jlocif>=0) nif = inDesc->inaxes[inDesc->jlocif];
  else nif = 1;
  if (inDesc->jlocs>=0) nstok = inDesc->inaxes[inDesc->jlocs];
  else nstok = 1;
  incs  = outDesc->incs;
  incf  = outDesc->incf;
  incif = outDesc->incif;
  nochan = outDesc->inaxes[outDesc->jlocf];
  if (outDesc->jlocif>=0) noif = outDesc->inaxes[outDesc->jlocif];
  else noif = 1;

  /* Zero work accumulator */
  n = 4 * inDesc->ncorr;
  for (i=0; i<n; i++) work[i] = 0.0;

  /* Accumulate order, channel, IF, poln */
  indx = 0;
  for (i=0; i<inDesc->ncorr; i++) {
    jndx = 4*(corChan[i] + corIF[i]*nchan + corStok[i]*nchan*nif);
    if ((inBuffer[indx+2]>0.0) && corMask[i]) {
      work[jndx]   += inBuffer[indx];
      work[jndx+1] += inBuffer[indx+1];
      work[jndx+2] += inBuffer[indx+2];
      work[jndx+3] += 1.0;
    }
    indx += inDesc->inaxes[0];
  } /* end accumulation loop */
  
  /* Normalize to output */
  /* Loop over Stokes */
  for (js=0; js<nstok; js++) {  /* Stokes loop */
    for (jif=0; jif<noif; jif++) {  /* IF loop */
      for (jf=0; jf<nochan; jf++) {  /* Frequency loop */
	jndx = 4*(jf + jif*nchan + js*nchan*nif);
	indx = js*incs + jif*incif + jf*incf;
	if (work[jndx+3]>0.0) {
	  outBuffer[indx]   = work[jndx]   / work[jndx+3];
	  outBuffer[indx+1] = work[jndx+1] / work[jndx+3];
	  outBuffer[indx+2] = work[jndx+2];
	} else {
	  outBuffer[indx]   = 0.0;
	  outBuffer[indx+1] = 0.0;
	  outBuffer[indx+2] = 0.0;
	}
      } /* end Frequency loop */
    } /* end IF loop */
  } /* end Stokes loop */

} /* end AvgFAver */

/**
 * Frequency boxcar smooth a visibility
 * \param inDesc   Input UV descriptor
 * \param outDesc  Output UV descriptor to be modified
 * \param NumChSMO Number of channels to smooth, should be odd
 * \param corChan  0-rel output channel numbers
 * \param corIF    0-rel output IF numbers
 * \param corStok  0-rel output Stokes parameter code.
 * \param corMask  Array of masks per correlator for selected channels/IFs
 *                 TRUE => select
 * \param inBuffer  Input buffer (data matrix)
 * \param outBuffer Output buffer (data matrix)
 * \param work      Work array twice the size of the output visibility
 * \param err       Error stack, returns if not empty.
 */
static void SmooF (ObitUVDesc *inDesc, ObitUVDesc *outDesc, olong NumChSmo,
		   olong *corChan, olong *corIF, olong *corStok, gboolean *corMask,
		   ofloat *inBuffer, ofloat *outBuffer, ofloat *work, 
		   ObitErr *err)
{
  olong i, j, half, n, indx, jndx, kndx, jf, jif, js, nochan, noif, off;
  olong nchan, nif, nstok, iincs, iincf, iincif, incs, incf, incif;

  /* error checks */
  if (err->error) return;

  half   = NumChSmo/2;
  iincs  = inDesc->incs;
  iincf  = inDesc->incf;
  iincif = inDesc->incif;
  nchan = inDesc->inaxes[inDesc->jlocf];
  if (inDesc->jlocif>=0) nif = inDesc->inaxes[inDesc->jlocif];
  else nif = 1;
  if (inDesc->jlocs>=0) nstok = inDesc->inaxes[inDesc->jlocs];
  else nstok = 1;
  incs  = outDesc->incs;
  incf  = outDesc->incf;
  incif = outDesc->incif;
  nochan = outDesc->inaxes[outDesc->jlocf];
  if (outDesc->jlocif>=0) noif = outDesc->inaxes[outDesc->jlocif];
  else noif = 1;

  /* Zero work accumulator */
  n = 4 * inDesc->ncorr;
  for (i=0; i<n; i++) work[i] = 0.0;

  /* Accumulate order, channel, IF, poln */
  for (js=0; js<nstok; js++) {  /* Stokes loop */
    for (jif=0; jif<nif; jif++) {  /* IF loop */
      for (jf=half; jf<nchan-half+1; jf++) {  /* Frequency loop */
	indx = js*iincs + jif*iincif + jf*iincf;
	jndx = 4*(jf + jif*nchan + js*nchan*nif);
	off = -half*iincf;
	/* Sums for smoothing */
	for (j=0; j<NumChSmo; j++) {
	  if (inBuffer[indx+off+2]>0.0) {  /* Valid? */
	    work[jndx]   += inBuffer[indx+off];
	    work[jndx+1] += inBuffer[indx+off+1];
	    work[jndx+2] += inBuffer[indx+off+2];
	    work[jndx+3] += 1.0;
	    off += iincf;
	  } /* end data valid */
	} /* end smoothing loop */
      } /* end freq loop */
      /* Copy end channels */
      jndx = 4*(half + jif*nchan + js*nchan*nif);
      for (jf=0; jf<half; jf++) {
	kndx = 4*(jf + jif*nchan + js*nchan*nif);
	work[kndx]   = work[jndx];
	work[kndx+1] = work[jndx+1];
	work[kndx+2] = work[jndx+2];
	work[kndx+3] = work[jndx+3];
      }
      jndx = 4*(nchan-half-1 + jif*nchan + js*nchan*nif);
      for (jf=half; jf>0; jf--) {
	kndx = 4*(nchan-jf + jif*nchan + js*nchan*nif);
	work[kndx]   = work[jndx];
	work[kndx+1] = work[jndx+1];
	work[kndx+2] = work[jndx+2];
	work[kndx+3] = work[jndx+3];
      }
    } /* end IF loop */
  } /* end stokes loop */

  /* Normalize to output */
  /* Loop over Stokes */
  for (js=0; js<nstok; js++) {  /* Stokes loop */
    for (jif=0; jif<noif; jif++) {  /* IF loop */
      for (jf=0; jf<nochan; jf++) {  /* Frequency loop */
	jndx = 4*(jf + jif*nchan + js*nchan*nif);
	indx = js*incs + jif*incif + jf*incf;
	/* Check for zero data */
	if ((work[jndx+3]>0.0) && !(( work[jndx]==0.0) && (work[jndx+1]==0.0))) {
	  outBuffer[indx]   = work[jndx]   / work[jndx+3];
	  outBuffer[indx+1] = work[jndx+1] / work[jndx+3];
	  outBuffer[indx+2] = work[jndx+2];
	} else {
	  outBuffer[indx]   = 0.0;
	  outBuffer[indx+1] = 0.0;
	  outBuffer[indx+2] = 0.0;
	}
      } /* end Frequency loop */
    } /* end IF loop */
  } /* end Stokes loop */

} /* end SmooF */

/**
 * Hanning smooth a visibility
 * \param inDesc   Input UV descriptor
 * \param outDesc  Output UV descriptor to be modified
 * \param doDescm  If TRUE drop every other channel
 * \param inBuffer  Input buffer (data matrix)
 * \param outBuffer Output buffer (data matrix)
 * \param work      Work array twice the size of the output visibility
 * \param err       Error stack, returns if not empty.
 */
static void Hann (ObitUVDesc *inDesc, ObitUVDesc *outDesc, gboolean doDescm,
		  ofloat *inBuffer, ofloat *outBuffer, ofloat *work, 
		  ObitErr *err)
{
  olong i, n, indx, jndx, jf, jif, js, nochan, noif, off, wincf;
  olong nchan, nif, nstok, iincs, iincf, iincif, incs, incf, incif;

  /* error checks */
  if (err->error) return;

  iincs  = inDesc->incs;
  iincf  = inDesc->incf;
  iincif = inDesc->incif;
  nchan = inDesc->inaxes[inDesc->jlocf];
  if (inDesc->jlocif>=0) nif = inDesc->inaxes[inDesc->jlocif];
  else nif = 1;
  if (inDesc->jlocs>=0) nstok = inDesc->inaxes[inDesc->jlocs];
  else nstok = 1;
  incs  = outDesc->incs;
  incf  = outDesc->incf;
  incif = outDesc->incif;
  nochan = outDesc->inaxes[outDesc->jlocf];
  if (outDesc->jlocif>=0) noif = outDesc->inaxes[outDesc->jlocif];
  else noif = 1;

  /* Zero work accumulator */
  n = 4 * inDesc->ncorr;
  for (i=0; i<n; i++) work[i] = 0.0;

  /* Accumulate order, channel, IF, poln */
  for (js=0; js<nstok; js++) {  /* Stokes loop */
    for (jif=0; jif<nif; jif++) {  /* IF loop */
      for (jf=0; jf<nchan; jf++) {  /* Frequency loop */
	indx = js*iincs + jif*iincif + jf*iincf;
	jndx = 4*(jf + jif*nchan + js*nchan*nif);
	off = -iincf;
	if ((jf>0) && (inBuffer[indx+off+2]>0.0)) {
	  /* Prior channel 0.25 wt */
	  work[jndx]   += 0.25*inBuffer[indx+off];
	  work[jndx+1] += 0.25*inBuffer[indx+off+1];
	  work[jndx+2] += 0.25*inBuffer[indx+off+2];
	  work[jndx+3] += 0.25;
	}
	/* Center channel 0.5 wt */
	if (inBuffer[indx+2]>0.0) {
	  work[jndx]   += 0.5*inBuffer[indx];
	  work[jndx+1] += 0.5*inBuffer[indx+1];
	  work[jndx+2] += 0.5*inBuffer[indx+2];
	  work[jndx+3] += 0.5;
	}
	/* Follow channel 0.25 wt */
	off = iincf;
	if ((jf<(nchan-1)) && (inBuffer[indx+off+2]>0.0)) {
	  work[jndx]   += 0.25*inBuffer[indx+off];
	  work[jndx+1] += 0.25*inBuffer[indx+off+1];
	  work[jndx+2] += 0.25*inBuffer[indx+off+2];
	  work[jndx+3] += 0.25;
	}
      } /* end freq loop */
    } /* end IF loop */
  } /* end stokes loop */

  /* Descimating output? */
  if (doDescm) wincf = 2;
  else         wincf = 1;
  
  /* Normalize to output */
  /* Loop over Stokes */
  for (js=0; js<nstok; js++) {  /* Stokes loop */
    for (jif=0; jif<noif; jif++) {  /* IF loop */
      for (jf=0; jf<nochan; jf++) {  /* Frequency loop */
	jndx = 4*(jf*wincf + jif*nchan + js*nchan*nif);
	indx = js*incs + jif*incif + jf*incf;
	/* Check for zero data */
	if ((work[jndx+3]>0.0) && !(( work[jndx]==0.0) && (work[jndx+1]==0.0))) {
	  outBuffer[indx]   = work[jndx]   / work[jndx+3];
	  outBuffer[indx+1] = work[jndx+1] / work[jndx+3];
	  outBuffer[indx+2] = work[jndx+2];
	} else {
	  outBuffer[indx]   = 0.0;
	  outBuffer[indx+1] = 0.0;
	  outBuffer[indx+2] = 0.0;
	}
      } /* end Frequency loop */
    } /* end IF loop */
  } /* end Stokes loop */

} /* end Hann */

/**
 * Duplicate channels 
 * \param inDesc   Input UV descriptor
 * \param outDesc  Output UV descriptor to be modified
 * \param nBloat   Number of duplicates of each input channel
 * \param unHann   If true, undo prior Hanning
 * \param inBuffer  Input buffer (data matrix)
 * \param outBuffer Output buffer (data matrix)
 * \param err       Error stack, returns if not empty.
 */
static void Bloat (ObitUVDesc *inDesc, ObitUVDesc *outDesc, olong nBloat,
		   gboolean unHann, ofloat *inBuffer, ofloat *outBuffer, 
		   ObitErr *err)
{
  olong j, indx, jndx, jf, jif, js, nochan, noif;
  olong nchan, nstok, iincs, iincf, iincif, incs, incf, incif;
  ofloat WtFact;

  /* error checks */
  if (err->error) return;

  WtFact = 1.0 / nBloat;     /* Weight factor */
  if (unHann) WtFact = 1.0;  /* Undoing Hanning? */
  iincs  = inDesc->incs;
  iincf  = inDesc->incf;
  iincif = inDesc->incif;
  nchan = inDesc->inaxes[inDesc->jlocf];
  if (inDesc->jlocs>=0) nstok = inDesc->inaxes[inDesc->jlocs];
  else nstok = 1;
  incs  = outDesc->incs;
  incf  = outDesc->incf;
  incif = outDesc->incif;
  nochan = outDesc->inaxes[outDesc->jlocf];
  if (outDesc->jlocif>=0) noif = outDesc->inaxes[outDesc->jlocif];
  else noif = 1;

  /* Duplicate to output */
  for (js=0; js<nstok; js++) {  /* Stokes loop */
    for (jif=0; jif<noif; jif++) {  /* IF loop */
      for (jf=0; jf<nchan; jf++) {  /* Frequency loop */
	indx = js*iincs + jif*iincif + jf*iincf;
	for (j=0; j<nBloat; j++) {
	  jndx = (jf*nBloat+j)*incf + jif*incif + js*incs;
	  /* Check for zero data, final channels */
	  if ((inBuffer[indx+2]>0.0) && ((jf*nBloat+j)<nochan)) {
	    outBuffer[jndx]   = inBuffer[indx];
	    outBuffer[jndx+1] = inBuffer[indx+1];
	    outBuffer[jndx+2] = inBuffer[indx+2] * WtFact;
	  } else {
	    outBuffer[jndx]   = 0.0;
	    outBuffer[jndx+1] = 0.0;
	    outBuffer[jndx+2] = 0.0;
	  }
	} /* end duplication */
      } /* end Frequency loop */
    } /* end IF loop */
  } /* end Stokes loop */
} /* end Bloat */

/**
 * Select visibility in frequency
 * \param inDesc    Input UV descriptor
 * \param outDesc   Output UV descriptor
 * \param BChan     First (1-rel) input channel selected
 * \param EChan     Highest (1-rel) channel selected
 * \param chinc     Increment of channels selected
 * \param BIF       First IF (1-rel) selected
 * \param EIF       Highest IF selected
 * \param inBuffer  Input buffer (data matrix)
 * \param outBuffer Output buffer (data matrix)
 */
static void FreqSel (ObitUVDesc *inDesc, ObitUVDesc *outDesc, 
		     olong BChan, olong EChan, olong chinc,
		     olong BIF, olong EIF,
		     ofloat *inBuffer, ofloat *outBuffer)
{
  olong indx, jndx, jf, jif, js, ojf, ojif;
  olong nstok, iincs, iincf, iincif, oincs, oincf, oincif;
  olong bchan=BChan-1, echan=EChan-1, bif=BIF-1, eif=EIF-1;

  /* Data increments */
  nstok  = outDesc->inaxes[outDesc->jlocs];
  iincs  = outDesc->incs;
  iincf  = outDesc->incf;
  iincif = outDesc->incif;
  oincs  = outDesc->incs;
  oincf  = outDesc->incf;
  oincif = outDesc->incif;

  /* Copy to output */
  /* Loop over Stokes */
  for (js=0; js<nstok; js++) {  /* Stokes loop */
    for (jif=bif; jif<=eif; jif++) {  /* IF loop */
      for (jf=bchan; jf<=echan; jf+=chinc) {  /* Frequency loop */
	jndx = js*iincs + jif*iincif + jf*iincf;
	ojf  = jf  - bchan;
	ojif = jif - bif;
	indx = js*oincs + ojif*oincif + ojf*oincf;
	outBuffer[indx]   = inBuffer[jndx];
	outBuffer[indx+1] = inBuffer[jndx+1];
	outBuffer[indx+2] = inBuffer[jndx+2];
      } /* end Frequency loop */
    } /* end IF loop */
  } /* end Stokes loop */

} /* end FreqSel */

/** 
 * Update FQ table for averaging 
 * \param inUV   Input UV data
 * \param chAvg  Number of channels averaged
 * \param fqid   Desired FQ ID to update
 * \param err    Error stack, returns if not empty.
 */
static void FQSel (ObitUV *inUV, olong chAvg, olong fqid, ObitErr *err)
{
  ObitTableFQ    *inTab=NULL;
  olong iFQver, highFQver;
  oint numIF;
  olong i, nif;
  odouble *freqOff=NULL;
  ofloat *chBandw=NULL;
  oint *sideBand=NULL;
  gchar *FQType = "AIPS FQ";
  gchar *routine = "ObitUVUtil:FQSel";

  /* error checks */
  if (err->error) return;
  g_assert (ObitUVIsA(inUV));

  /* How many FQ tables  */
  highFQver = ObitTableListGetHigh (inUV->tableList, FQType);

  /* Are there any? */
  if (highFQver <= 0) return;

  /* Should only be one FQ table */
  iFQver = 1;
  if (inUV->myDesc->jlocif>=0) 
    nif = inUV->myDesc->inaxes[inUV->myDesc->jlocif];
  else
    nif = 1;

  /* Get input table */
  numIF = 0;
  inTab = 
    newObitTableFQValue (inUV->name, (ObitData*)inUV, &iFQver, OBIT_IO_ReadOnly, 
			 numIF, err);
  if (err->error) Obit_traceback_msg (err, routine, inUV->name);
  /* Find it? */
   Obit_return_if_fail(((inTab!=NULL) || (nif<=1)), err,
		      "%s: Could not find FQ table for %s %d IFs", 
		      routine, inUV->name, nif);

   /* Fetch values */
   ObitTableFQGetInfo (inTab, fqid, &nif, &freqOff, &sideBand, &chBandw, err);
   if (err->error) Obit_traceback_msg (err, routine, inTab->name);

   /* Update channel widths */
   for (i=0; i<nif; i++) chBandw[i] *= (ofloat)chAvg;

   /* Save values */
   ObitTableFQPutInfo (inTab, fqid, nif, freqOff, sideBand, chBandw, err);
   if (err->error) Obit_traceback_msg (err, routine, inTab->name);

   /* Cleanup */
   if (freqOff)  g_free(freqOff);
   if (sideBand) g_free(sideBand);
   if (chBandw)  g_free(chBandw);
   inTab = ObitTableFQUnref(inTab);
  
 } /* end FQSel */

/**
 * Low accuracy inverse Sinc function
 * \param arg Argument
 * \return angle in radians
 */
static ofloat InvSinc(ofloat arg)
{
  odouble x0, x1, a;
  olong i, n=1000;

  /* Some iterations of Newton-Raphson starting with the value near 1.0 */
  x1 = 0.001;
  for (i=0; i<n; i++) {
    x0 = x1;
    a = x0 * G_PI;
    x1 = x0 - ((sin(a)/a) - arg) / ((a*cos(a) - G_PI*sin(a))/(a*a));
    if (fabs(x1-x0)<1.0e-6) break;  /* Convergence test */
  }

  return (ofloat)x1;
} /* end InvSinc */
