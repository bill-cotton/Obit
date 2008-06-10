/* $Id: ObitUVUtil.c,v 1.28 2008/03/05 14:21:35 bcotton Exp $   */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2004-2008                                          */
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

#include "ObitUVUtil.h"
#include "ObitTableSUUtil.h"
#include "ObitTableNXUtil.h"
#include "ObitTablePSUtil.h"
#include "ObitTableAN.h"
#include "ObitPrecess.h"
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

/** Average visibility in frequency */
static void AvgFAver (ObitUVDesc *inDesc, ObitUVDesc *outDesc, 
		      olong NumChAvg, olong *ChanSel, gboolean doAvgAll, 
		      olong *corChan, olong *corIF, olong *corStok, gboolean *corMask,
		      ofloat *inBuffer, ofloat *outBuffer, ofloat *work, 
		      ObitErr *err);
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
		    "AIPS PL", "AIPS NI",
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
  outUV->buffer = inUV->buffer;
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
		    "AIPS PL", "AIPS NI",
		    NULL};
  gchar *sourceInclude[] = {"AIPS SU", NULL};
  olong i, j, indx;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong NPIO;
  ofloat work[3];
  gboolean incompatible, same;
  ObitUVDesc *in1Desc, *in2Desc, *outDesc;
  gchar *today=NULL;
  gchar *routine = "ObitUVUtilVisDivide";

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
    outUV->buffer = inUV1->buffer;
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
  if (incompatible) {
     Obit_log_error(err, OBIT_Error,"%s inUV1 and inUV2 have incompatible structures",
		   routine);
      return ;
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

    /* compatability check */
    incompatible = in1Desc->numVisBuff!=in2Desc->numVisBuff;
    if (incompatible) break;

    /* Modify data */
    for (i=0; i<in1Desc->numVisBuff; i++) { /* loop over visibilities */
      /* compatability check - check time and baseline code */
      indx = i*in1Desc->lrec ;
      incompatible = 
	inUV1->buffer[indx+in1Desc->iloct]!=inUV1->buffer[indx+in2Desc->iloct] ||
	inUV1->buffer[indx+in1Desc->ilocb]!=inUV1->buffer[indx+in2Desc->ilocb];
      if (incompatible) break;

      indx += in1Desc->nrparm;
      for (j=0; j<in1Desc->ncorr; j++) { /* loop over correlations */
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
  
} /* end ObitUVUtilVisDivide */

/**
 * Subtract the visibilities in one ObitUV from those in another
 * outUV = inUV1 - inUV2
 * \param inUV1    First input uv data, no calibration/selection
 * \param inUV2    Second input uv data, calibration/selection allowed
 *                 inUV2 should have the same structure, no. vis etc
 *                 as inUV1.
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
		    "AIPS PL", "AIPS NI",
		    NULL};
  gchar *sourceInclude[] = {"AIPS SU", NULL};
  olong i, j, indx, firstVis;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong NPIO;
  gboolean incompatible, same, doCalSelect;
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
    outUV->buffer = inUV1->buffer;
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
  if (incompatible) {
     Obit_log_error(err, OBIT_Error,"%s inUV1 and inUV2 have incompatible structures",
		   routine);
      return ;
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
      /* compatability check - check time and baseline code */
      indx = i*in1Desc->lrec ;
      incompatible = 
	inUV1->buffer[indx+in1Desc->iloct]!=inUV2->buffer[indx+in2Desc->iloct] ||
	inUV1->buffer[indx+in1Desc->ilocb]!=inUV2->buffer[indx+in2Desc->ilocb];
      if (incompatible) break;

      indx += in1Desc->nrparm;
      for (j=0; j<in1Desc->ncorr; j++) { /* loop over correlations */
	/* subtract */
	if ((inUV1->buffer[indx+2]>0.0) && (inUV2->buffer[indx+2]>0.0)) {
	  inUV1->buffer[indx]   -= inUV2->buffer[indx];
	  inUV1->buffer[indx+1] -= inUV2->buffer[indx+1];
	} else {
	  inUV1->buffer[indx]   = 0.0;
	  inUV1->buffer[indx+1] = 0.0;
	  inUV1->buffer[indx+2] = 0.0;
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
  
  /* close files */
  iretCode = ObitUVClose (inUV1, err);
  if (err->error) Obit_traceback_msg (err, routine, inUV1->name);
  
  iretCode = ObitUVClose (inUV2, err);
  if (err->error) Obit_traceback_msg (err, routine, inUV2->name);
  
  oretCode = ObitUVClose (outUV, err);
  if (err->error) Obit_traceback_msg (err, routine, outUV->name);
  
} /* end ObitUVUtilVisSub */

/**
 * Compare the visibilities in one ObitUV with those in another.
 * Return the RMS of the real and imaginary differences 
 * divided by the amplitude of inUV2.
 * Only valid visibilities compared, zero amplitudes ignored.
 * \param inUV1    Input uv data numerator, no calibration/selection
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
  olong i, j, indx, count;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong NPIO;
  ofloat amp, rms = -1.0;
  odouble sum;
  gboolean incompatible;
  ObitUVDesc *in1Desc, *in2Desc;
  gchar *routine = "ObitUVUtilVisCompare";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return rms;
  g_assert (ObitUVIsA(inUV1));
  g_assert (ObitUVIsA(inUV2));

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
  if (incompatible) {
     Obit_log_error(err, OBIT_Error,"%s inUV1 and inUV2 have incompatible structures",
		   routine);
      return rms;
 }

  /* we're in business, loop comparing */
  count = 0;
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
      /* compatability check - check time and baseline code */
      indx = i*in1Desc->lrec ;
      incompatible = 
	inUV1->buffer[indx+in1Desc->iloct]!=inUV1->buffer[indx+in2Desc->iloct] ||
	inUV1->buffer[indx+in1Desc->ilocb]!=inUV1->buffer[indx+in2Desc->ilocb];
      if (incompatible) break;

      indx += in1Desc->nrparm;
      for (j=0; j<in1Desc->ncorr; j++) { /* loop over correlations */
	/* Statistics  */
	amp = inUV2->buffer[indx]*inUV2->buffer[indx] + 
	  inUV2->buffer[indx+1]*inUV2->buffer[indx+1];
	if ((inUV1->buffer[indx+2]>0.0) && (inUV2->buffer[indx+2]>0.0) && (amp>0.0)) {
	  amp = sqrt(amp);
	  sum += ((inUV1->buffer[indx] - inUV2->buffer[indx]) + 
	    (inUV1->buffer[indx] - inUV2->buffer[indx])) / amp;
	  sum += ((inUV1->buffer[indx+1] - inUV2->buffer[indx+1]) + 
	    (inUV1->buffer[indx+1] - inUV2->buffer[indx+1])) / amp;
	  count += 2;
	}
	indx += in1Desc->inaxes[0];
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
  if (count>0) rms = sum / count;
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
  olong num, i, lrec, iRow, ver, startVis=1, curVis, itemp;
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
  lrec = inUV->myDesc->lrec;  /* Size of visibility */
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

      itemp  = vis[inUV->myDesc->ilocb];           /* Subarray number */
      suba = (100.0*(vis[inUV->myDesc->ilocb]-itemp)) + 1.5;
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
 *     If no SU table, a single entry is returned.
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
  gboolean want, allCal, allNCal, doTime, doPS, *good=NULL;
  gchar souCode[5];
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
    strncpy (out->SUlist[0]->SourceName, inUV->myDesc->object, MIN(20,UVLEN_VALUE));
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
  gchar *exclude[]={"AIPS CL","AIPS SN","AIPS FG","AIPS CQ","AIPS WX",
		    "AIPS AT","AIPS CT","AIPS OB","AIPS IM","AIPS MC",
		    "AIPS PC","AIPS NX","AIPS TY","AIPS GC","AIPS HI",
		    "AIPS PL","AIPS NI","AIPS BP","AIPS OF","AIPS PS",
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
  gboolean doAvgAll;
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
  ChanSel = NULL;
  if (!ObitInfoListGetP(inUV->info, "ChanSel", &type, dim, (gpointer)&ChanSel)) {
    ChanSel = defSel;  /* Use default = channels 1 => n */
  }
  /* ChanSel all zero? => default */
  if ((ChanSel[0]<=0) && (ChanSel[1]<=0) && (ChanSel[2]<=0) && (ChanSel[3]<=0)) {
    ChanSel = defSel;  /* Use default = channels 1 => n */
  }


  /* Create scratch? */
  if (scratch) {
    if (outUV) outUV = ObitUVUnref(outUV);
    outUV = newObitUVScratch (inUV, err);
  } else { /* non scratch output must exist - clone from inUV */
    ObitUVClone (inUV, outUV, err);
  }
  if (err->error) Obit_traceback_val (err, routine, inUV->name, inUV);

  /* Selection/calibration/editing of input? */
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
  scale = AvgFSetDesc (inDesc, outDesc, NumChAvg, ChanSel, doAvgAll, 
		       corChan, corIF, corStok, corMask, err);
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
  if (work) g_free(work);       work    = NULL;
  if (corChan) g_free(corChan); corChan = NULL;
  if (corIF) g_free(corIF);     corIF   = NULL;
  if (corStok) g_free(corStok); corStok = NULL;
  if (corMask) g_free(corMask); corMask = NULL;
  
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
  gchar *exclude[]={"AIPS CL","AIPS SN","AIPS FG","AIPS CQ","AIPS WX",
		    "AIPS AT","AIPS CT","AIPS OB","AIPS IM","AIPS MC",
		    "AIPS PC","AIPS NX","AIPS TY","AIPS GC","AIPS HI",
		    "AIPS PL", "AIPS NI","AIPS BP","AIPS OF","AIPS PS",
		    NULL};
  gchar *sourceInclude[] = {"AIPS SU", NULL};
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  olong ncorr, nrparm, numAnt, numBL;
  ObitIOAccess access;
  ObitUVDesc *inDesc, *outDesc;
  olong suba, lastSourceID, curSourceID, lastSubA, lastFQID=-1;
  gchar *today=NULL;
  ofloat cbase, timeAvg, curTime, startTime, endTime, lastTime=-1.0;
  ofloat *accVis=NULL, *accRP=NULL;
  ofloat *inBuffer, *outBuffer;
  olong ant1, ant2, blindx, *blLookup=NULL;
  olong i, j, ivis=0, iindx, jndx, indx, NPIO, itemp;
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
  if (timeAvg<=(1.0/60.0)) timeAvg = 1.0;
  timeAvg /= 1440.0;  /* convert to days */

  /* Create scratch? */
  if (scratch) {
    if (outUV) outUV = ObitUVUnref(outUV);
    outUV = newObitUVScratch (inUV, err);
  } else { /* non scratch output must exist - clone from inUV */
    ObitUVClone (inUV, outUV, err);
  }
  if (err->error) Obit_traceback_val (err, routine, inUV->name, inUV);

  /* Selection/calibration/editing of input? */
  doCalSelect = FALSE;
  ObitInfoListGetTest(inUV->info, "doCalSelect", &type, dim, &doCalSelect);
  if (doCalSelect) access = OBIT_IO_ReadCal;
  else access = OBIT_IO_ReadOnly;

  /* test open to fully instantiate input and see if it's OK */
  iretCode = ObitUVOpen (inUV, access, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_val (err, routine, inUV->name, outUV);

  /* copy Descriptor */
  outUV->myDesc = ObitUVDescCopy(inUV->myDesc, outUV->myDesc, err);
  inBuffer = inUV->buffer;  /* Local copy of buffer pointer */
 
  /* Output creation date today */
  today = ObitToday();
  strncpy (outUV->myDesc->date, today, UVLEN_VALUE-1);
  if (today) g_free(today);
  
 /* Set number of output vis per read to 1 */
  NPIO = 1;
  ObitInfoListGetTest(inUV->info, "nVisPIO", &type, dim, &NPIO);
  itemp = 1;
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

 /* Get descriptors */
  inDesc  = inUV->myDesc;
  outDesc = outUV->myDesc;

  /* Create work arrays for averaging */
  suba    = 1;
  numAnt  = inUV->myDesc->numAnt[suba-1];/* actually highest antenna number */
  numBL   = (numAnt*(numAnt+1))/2;  /* Include auto correlations */
  ncorr   = inUV->myDesc->ncorr;
  nrparm  = inUV->myDesc->nrparm;
  accVis  = g_malloc0(4*numBL*ncorr*sizeof(ofloat));     /* Vis */
  accRP   = g_malloc0(numBL*(nrparm+1)*sizeof(ofloat));  /* Rand. parm */

  /* Baseline lookup table */
  blLookup = g_malloc0 (numAnt* sizeof(olong));
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
  startTime = -1.0e20;
  endTime   =  1.0e20;
  lastSourceID = -1;
  curSourceID  = 0;
  outDesc->numVisBuff = 0;
  outBuffer = outUV->buffer;  /* Local copy of buffer pointer */
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
	lastTime  = curTime;
	endTime   = startTime + timeAvg;
	lastSourceID = curSourceID;
      }

      /* Still in current interval/source? */
      if ((curTime<endTime) && (curSourceID == lastSourceID) && 
	  (inDesc->firstVis<=inDesc->nvis) && (iretCode==OBIT_IO_OK)) {
	/* accumulate */
	cbase = inBuffer[iindx+inUV->myDesc->ilocb]; /* Baseline */
	ant1 = (cbase / 256.0) + 0.001;
	ant2 = (cbase - ant1 * 256) + 0.001;
	lastSubA = (olong)(100.0 * (cbase -  ant1 * 256 - ant2) + 0.5);
	/* Baseline index this assumes a1<=a2 always */
	blindx =  blLookup[ant1-1] + ant2-ant1;
	if (inDesc->ilocfq>=0) lastFQID = (olong)(inBuffer[iindx+inDesc->ilocfq]+0.5);
	else lastFQID = 0;
	lastTime = curTime;
	
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
	      outBuffer[indx++] = accRP[jndx+i+1];
	    } /* End random parameter loop */

	    /* Average vis data */
	    for (j=0; j<ncorr; j++) {
	      jndx = j*4 + blindx*4*ncorr;
	      if (accVis[jndx]>0.0) {
		accVis[jndx+1] /= accVis[jndx];
		accVis[jndx+2] /= accVis[jndx];
	      }
	      /* Copy to output buffer */
	      outBuffer[indx++] = accVis[jndx+1];
	      outBuffer[indx++] = accVis[jndx+2];
	      outBuffer[indx++] = accVis[jndx+3];
	    } /* end loop over correlators */

	    /* Write output one vis at a time */
	    outDesc->numVisBuff = 1;
	    oretCode = ObitUVWrite (outUV, outBuffer, err);
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
  if (accVis)   g_free(accVis);   accVis   = NULL;
  if (accRP)    g_free(accRP);    accRP    = NULL;
  if (blLookup) g_free(blLookup); blLookup = NULL;
  
  /* close files */
  iretCode = ObitUVClose (inUV, err);
  oretCode = ObitUVClose (outUV, err);
  if ((oretCode!=OBIT_IO_OK) || (iretCode!=OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, outUV->name, outUV);
  
  return outUV;
} /* end ObitUVUtilAvgT */

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
  ofloat lastTime=-1.0, *inBuffer;
  olong ncorr, indx, iindx, lastSourceID, curSourceID;
  gboolean gotOne, done, isVLA;
  odouble GSTiat0, DegDay, ArrayX, ArrayY, ArrayZ;
  ofloat dataIat, ArrLong;
  ofloat startTime, endTime, curTime;
  olong numTime, goodCnt[500], badCnt[500], timeCnt[500], timeSou[500];
  ofloat timeSum[500];
  gchar *routine = "ObitUVUtilCount";

  /* error checks */
  if (err->error) return outList;
  g_assert (ObitUVIsA(inUV));

  /* Initialize sums */
  numTime = 0;
  for (i=0; i<100; i++) {
    goodCnt[i] = badCnt[i] = timeCnt[i] = timeSou[i] = 0;
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

      curTime = inBuffer[iindx+inDesc->iloct]; /* Time */
      if (inDesc->ilocsu>=0) curSourceID = inBuffer[iindx+inDesc->ilocsu];
      if (startTime < -1000.0) {  /* Set time window etc. if needed */
	startTime = curTime;
	lastTime  = curTime;
	endTime   = startTime + timeInt;
	lastSourceID = curSourceID;
      }

      /* Still in current interval */
      if ((curTime<endTime) && (curSourceID == lastSourceID) && 
	  (inDesc->firstVis<=inDesc->nvis) && (iretCode==OBIT_IO_OK)) {
	lastTime = curTime;

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
				 OBIT_IO_ReadOnly, 0, 0, err);
  GSTiat0 = ANTable->GSTiat0;
  DegDay  = ANTable->DegDay;
  ArrayX  = ANTable->ArrayX;
  ArrayY  = ANTable->ArrayY;
  ArrayZ  = ANTable->ArrayZ;
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
  ObitInfoListAlwaysPut (outList, "Source", OBIT_long,   dim, timeSou);
  ObitInfoListAlwaysPut (outList, "Count",  OBIT_long,   dim, goodCnt);
  ObitInfoListAlwaysPut (outList, "Bad",    OBIT_long,   dim, badCnt);
  ObitInfoListAlwaysPut (outList, "LST",    OBIT_float, dim, timeSum);
  /* End of processing */

  /* Cleanup */
 cleanup:
  /* close file */
  iretCode = ObitUVClose (inUV, err);
  if ((iretCode!=OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, inUV->name, outList);
  
  return outList;
} /* end ObitUVUtilAvgT  */

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
  if (inDesc->jlocf>=0) nif = inDesc->inaxes[inDesc->jlocif];
  else nif = 1;
  if (inDesc->jlocs>=0) nstok = inDesc->inaxes[inDesc->jlocs];
  else nstok = 1;

  numChan = nchan / NumChAvg;  /* Number of output channels */
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
  numChan = inDesc->inaxes[inDesc->jlocf] / NumChAvg;
  /* Get frequency of first average */
  sum = 0.0; count = 0;
    sum2 = 0.0; count2 = 0;
  for (i=0; i<inDesc->ncorr; i++) {
    if (corChan[i]>=NumChAvg) break;  /* finished? */
    ii = corChan[i] + corIF[i]*nchan;
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
  if (inDesc->jlocf>=0) nif = inDesc->inaxes[inDesc->jlocif];
  else nif = 1;
  if (inDesc->jlocs>=0) nstok = inDesc->inaxes[inDesc->jlocs];
  else nstok = 1;
  incs  = outDesc->incs;
  incf  = outDesc->incf;
  incif = outDesc->incif;
  nochan = outDesc->inaxes[outDesc->jlocf];
  if (outDesc->jlocf>=0) noif = outDesc->inaxes[outDesc->jlocif];
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
