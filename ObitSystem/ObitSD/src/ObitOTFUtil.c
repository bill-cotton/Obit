/* $Id$  */
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

#include <time.h>
#include <gsl/gsl_randist.h>
#include "ObitThread.h"
#include "ObitOTFUtil.h"
#include "ObitOTFGetSoln.h"
#include "ObitTableOTFIndex.h"
#include "ObitTableOTFTargetUtil.h"
#include "ObitImageUtil.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitOTFUtil.c
 * ObitOTF class utility function definitions.
 */

/*---------------Private structures----------------*/
/* SubImage threaded function argument */
typedef struct {
  /* OTF data set to model and subtract from current buffer */
  ObitOTF       *otfdata;
  /* First (1-rel) record in otfdata buffer to process this thread */
  olong        first;
  /* Highest (1-rel) record in otfdata buffer to process this thread  */
  olong        last;
  /* thread number  */
  olong        ithread;
  /* Obit error stack object */
  ObitErr      *err;
  /* Image Interpolator */
  ObitFInterpolate *Interp;
  /* Scaling factor for model */
  ofloat factor;
  /* Work arrays the size of ndetect */
  ofloat *xpos, *ypos;
} SubImageFuncArg;

/*---------------Private function prototypes----------------*/
/** Private: Subtract an image interpolator from a buffer of data. */
void ObitOTFUtilSubImageBuff (ObitOTF *in, ObitFInterpolate *image, ofloat factor, 
			      olong nThread, SubImageFuncArg **args, ObitErr *err);

/** Private: Convert an ObitOTFDesc to an ObitImageDesc */
static void 
ObitOTFUtilOTF2ImageDesc(ObitOTFDesc *OTFDesc, ObitImageDesc *imageDesc, 
			 gchar *Proj);

/** Private: Get Date string for current date */
static void ObitOTFUtilCurDate (gchar *date, olong len);

/** Private: Threaded OTFUtilSubImageBuff */
static gpointer ThreadOTFUtilSubImageBuff (gpointer arg);

/** Private: Make arguments for Threaded OTFUtilSubImageBuff */
static glong MakeOTFUtilSubImageArgs (ObitOTF *in, ObitErr *err, 
				      SubImageFuncArg ***args);

/** Private: Delete arguments for Threaded OTFUtilSubImageBuff */
static void KillOTFUtilSubImageArgs (olong nargs, SubImageFuncArg **args);

#define MAXSAMPLE 10000   /* Maximum number of samples in a scan */

/*----------------------Public functions---------------------------*/
/**
 * Subtract a 2D ObitFArray from an OTF
 * \param inOTF  Input OTF 
 * \param outOTF Output OTF, must already be defined
 * \param image  Image plane to subtract
 * \param desc   Image descriptor for image
 * \param err    Error stack
 */
void ObitOTFUtilSubImage(ObitOTF *inOTF, ObitOTF *outOTF, ObitFArray *image, 
			 ObitImageDesc *desc, ObitErr *err) 
{
  ObitFInterpolate *imageInt=NULL;
  ObitIOCode retCode;
  gboolean doCalSelect, done, same;
  olong firstRec;
  ObitInfoType type;
  ObitIOAccess access;
  gint32 dim[MAXINFOELEMDIM];
  ofloat scale = 1.0;
  olong NPIO, nThreads;
  SubImageFuncArg **targs=NULL;
  /* Don't copy Cal and Soln or data or flag tables */
  gchar *exclude[]={"OTFSoln", "OTFCal", "OTFScanData", "OTFFlag", NULL};
  gchar *routine = "ObitOTFUtilSubImage";

    /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitFArrayIsA(image));
  g_assert (ObitImageDescIsA(desc));
  g_assert (ObitOTFIsA(inOTF));
  g_assert (ObitOTFIsA(outOTF));

  /* An interpolator for the image */
  imageInt = newObitFInterpolateCreate (image->name, image, desc, 2);

  /* are input and utput the same? */
  same = ObitOTFSame (inOTF, outOTF, err);
  
  /* Calibration wanted? */ 
  doCalSelect = FALSE;
  ObitInfoListGetTest(inOTF->info, "doCalSelect", &type, (gint32*)dim, 
		      &doCalSelect);
  if (doCalSelect) access = OBIT_IO_ReadCal;
  else access = OBIT_IO_ReadWrite;

  /* Make sure NPIO the same */
  NPIO = 1000;
  dim[0] = 1;
  ObitInfoListGetTest (inOTF->info, "nRecPIO", &type, dim,  &NPIO);
  ObitInfoListAlwaysPut (inOTF->info, "nRecPIO", OBIT_long, dim,  &NPIO);
  ObitInfoListAlwaysPut (outOTF->info, "nRecPIO", OBIT_long, dim,  &NPIO);

 /* Open Input Data */
  retCode = ObitOTFOpen (inOTF, access, err);
  if ((retCode != OBIT_IO_OK) || (err->error>0)) goto cleanup;

  if (!same) {
  /* use same data buffer on input and output 
     so don't assign buffer for output */
    if (outOTF->buffer) ObitIOFreeBuffer(outOTF->buffer); /* free existing */
    outOTF->buffer     = inOTF->buffer;
    outOTF->bufferSize = inOTF->bufferSize;

   /* Open Output Data */
    retCode = ObitOTFOpen (outOTF, OBIT_IO_WriteOnly, err) ;
    if ((retCode != OBIT_IO_OK) || (err->error>0)) goto cleanup;
  }

  /* Copy tables before data if they are different */
  if (!same) {
    retCode = ObitOTFCopyTables (inOTF, outOTF, exclude, NULL, err);
    if (err->error) goto cleanup;
  }

  /* Close and reopen input to init calibration which will have 
     been disturbed by the table copy */
  retCode = ObitOTFClose (inOTF, err);
  if (err->error) goto cleanup;

  retCode = ObitOTFOpen (inOTF, access, err);
  if ((retCode != OBIT_IO_OK) || (err->error>0)) goto cleanup;

  /* Scaling if needed to convert from image brightness to data 
     square of ratio of beam areas */
  scale = inOTF->myDesc->beamSize / desc->beamMaj;
  scale = scale * scale;
  /*fprintf (stderr,"Scaling image to data by %f\n",scale); *//*debug */

  /* Setup Threading */
  nThreads = MakeOTFUtilSubImageArgs (inOTF, err, &targs);

  /* Loop over data */
  done = (retCode != OBIT_IO_OK);
  while (!done) {

    /* read buffer */#
    retCode = ObitOTFRead (inOTF, NULL, err);
    if (err->error) goto cleanup;
    done = (retCode == OBIT_IO_EOF); /* done? */
    if (done) break;
    firstRec = inOTF->myDesc->firstRec;

     /* How many? */
     outOTF->myDesc->numRecBuff = inOTF->myDesc->numRecBuff;

     if (outOTF->myDesc->numRecBuff>0) {
       /* Subtract image from this buffer */
       ObitOTFUtilSubImageBuff (outOTF, imageInt, scale, nThreads, targs, err);
       
       /* Write buffer */
       retCode = ObitOTFWrite (outOTF, NULL, err);
     }
     if (err->error) goto cleanup;

     /* suppress vis number update if rewriting the same file */
     if (same) {
       outOTF->myDesc->firstRec = firstRec;
       ((ObitOTFDesc*)(outOTF->myIO->myDesc))->firstRec = firstRec;
     }

  } /* end scan loop */
  
 cleanup:
  /* unset output buffer (may be multiply deallocated ;'{ ) */
  if (!same) {
    outOTF->buffer = NULL;
    outOTF->bufferSize = 0;
    retCode = ObitOTFClose (outOTF, err); /* Close output */
  } 

  /* Close input */
  retCode = ObitOTFClose (inOTF, err);

  /* cleanup */
  imageInt = ObitFInterpolateUnref(imageInt);
  if (err->error) Obit_traceback_msg (err, routine, inOTF->name);
  /* Shutdown Threading */
  KillOTFUtilSubImageArgs (nThreads, targs);

} /* end ObitOTFUtilSubImage */

/**
 * Replace the data in an OTF with the model values from FArray image
 * \param inOTF  Input OTF 
 * \param outOTF Output OTF, must already be defined
 * \param image  Image plane to subtract
 * \param desc   Image descriptor for image
 * \param err    Error stack
 */
void ObitOTFUtilModelImage(ObitOTF *inOTF, ObitOTF *outOTF, ObitFArray *image, 
			   ObitImageDesc *desc, ObitErr *err) 
{
  ObitFInterpolate *imageInt=NULL;
  ObitIOCode retCode;
  gboolean doCalSelect, done, same;
  olong firstRec;
  ObitInfoType type;
  ObitIOAccess access;
  gint32 dim[MAXINFOELEMDIM];
  ofloat scale = 1.0;
  olong NPIO;
  /* Don't copy Cal and Soln or data or flag tables */
  gchar *exclude[]={"OTFSoln", "OTFCal", "OTFScanData", "OTFFlag", NULL};
  gchar *routine = "ObitOTFUtilModelImage";

    /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitFArrayIsA(image));
  g_assert (ObitImageDescIsA(desc));
  g_assert (ObitOTFIsA(inOTF));
  g_assert (ObitOTFIsA(outOTF));

  /* An interpolator for the image */
  imageInt = newObitFInterpolateCreate (image->name, image, desc, 2);

  /* are input and utput the same? */
  same = ObitOTFSame (inOTF, outOTF, err);
  
  /* Calibration wanted? */ 
  doCalSelect = FALSE;
  ObitInfoListGetTest(inOTF->info, "doCalSelect", &type, dim, &doCalSelect);
  if (doCalSelect) access = OBIT_IO_ReadCal;
  else access = OBIT_IO_ReadWrite;

  /* Make sure NPIO the same */
  NPIO = 1000;
  dim[0] = 1;
  ObitInfoListGetTest (inOTF->info, "nRecPIO", &type, dim,  &NPIO);
  ObitInfoListAlwaysPut (inOTF->info, "nRecPIO", OBIT_long, dim,  &NPIO);
  ObitInfoListAlwaysPut (outOTF->info, "nRecPIO", OBIT_long, dim,  &NPIO);

 /* Open Input Data */
  retCode = ObitOTFOpen (inOTF, access, err);
  if ((retCode != OBIT_IO_OK) || (err->error>0)) goto cleanup;

  if (!same) {
  /* use same data buffer on input and output 
     so don't assign buffer for output */
    if (outOTF->buffer) ObitIOFreeBuffer(outOTF->buffer); /* free existing */
    outOTF->buffer     = inOTF->buffer;
    outOTF->bufferSize = inOTF->bufferSize;

   /* Open Output Data */
    retCode = ObitOTFOpen (outOTF, OBIT_IO_WriteOnly, err) ;
    if ((retCode != OBIT_IO_OK) || (err->error>0)) goto cleanup;
  }

  /* Copy tables before data if they are different */
  if (!same) {
    retCode = ObitOTFCopyTables (inOTF, outOTF, exclude, NULL, err);
    if (err->error) goto cleanup;
  }

  /* Close and reopen input to init calibration which will have 
     been disturbed by the table copy */
  retCode = ObitOTFClose (inOTF, err);
  if (err->error) goto cleanup;

  retCode = ObitOTFOpen (inOTF, access, err);
  if ((retCode != OBIT_IO_OK) || (err->error>0)) goto cleanup;

  /* Scaling if needed to convert from image brightness to data 
     square of ratio of beam areas */
  scale = inOTF->myDesc->beamSize / desc->beamMaj;
  scale = scale * scale;
  /*fprintf (stderr,"Scaling image to data by %f\n",scale); *//*debug */

  /* Loop over data */
  done = (retCode != OBIT_IO_OK);
  while (!done) {

    /* read buffer */#
    retCode = ObitOTFRead (inOTF, NULL, err);
    if (err->error) goto cleanup;
    done = (retCode == OBIT_IO_EOF); /* done? */
    if (done) break;
    firstRec = inOTF->myDesc->firstRec;

     /* How many? */
     outOTF->myDesc->numRecBuff = inOTF->myDesc->numRecBuff;

     if (outOTF->myDesc->numRecBuff>0) {
       /* Replace data with model in this buffer */
       ObitOTFUtilModelImageBuff (outOTF, imageInt, scale, err);
       
       /* Write buffer */
       retCode = ObitOTFWrite (outOTF, NULL, err);
     }
     if (err->error) goto cleanup;

     /* suppress vis number update if rewriting the same file */
     if (same) {
       outOTF->myDesc->firstRec = firstRec;
       ((ObitOTFDesc*)(outOTF->myIO->myDesc))->firstRec = firstRec;
     }

  } /* end scan loop */
  
 cleanup:
  /* unset output buffer (may be multiply deallocated ;'{ ) */
  if (!same) {
    outOTF->buffer = NULL;
    outOTF->bufferSize = 0;
    retCode = ObitOTFClose (outOTF, err); /* Close output */
  } 

  /* Close input */
  retCode = ObitOTFClose (inOTF, err);

  /* cleanup */
  imageInt = ObitFInterpolateUnref(imageInt);
  if (err->error) Obit_traceback_msg (err, routine, inOTF->name);

} /* end ObitOTFUtilModelImage */

/**
 * Multiply all data in an OTF by a scale factor and add an offset
 * out = in*scale + offset
 * \param inOTF  Input OTF 
 * \param outOTF Output OTF, must already be defined
 * \param scale  scaling factor for data
 * \param offset additive term to be applied to all data.
 * \param err    Error stack
 */
void ObitOTFUtilScale(ObitOTF *inOTF, ObitOTF *outOTF, ofloat scale, ofloat offset,
		      ObitErr *err)
{
  ObitIOCode retCode;
  gboolean doCalSelect, done;
  ObitInfoType type;
  ObitIOAccess access;
  gint32 dim[MAXINFOELEMDIM];
  ofloat *data, fblank = ObitMagicF();
  olong i, j, ndetect, ndata;
  olong incdatawt;
  ObitOTFDesc* desc;
 /* Don't copy Cal and Soln or data or flag tables */
  gchar *exclude[]={"OTFSoln", "OTFCal", "OTFScanData", "OTFFlag", NULL};
  gchar *routine = "ObitOTFUtilScale";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitOTFIsA(inOTF));
  g_assert (ObitOTFIsA(outOTF));

  /* Local pointers */
  desc = inOTF->myDesc;
  incdatawt = desc->incdatawt; /* increment in data-wt axis */

  /* Calibration wanted? */ 
  doCalSelect = FALSE;
  ObitInfoListGetTest(inOTF->info, "doCalSelect", &type, (gint32*)dim, &doCalSelect);
  if (doCalSelect) access = OBIT_IO_ReadCal;
  else access = OBIT_IO_ReadWrite;

 /* Open Input Data */
  retCode = ObitOTFOpen (inOTF, access, err);
  if ((retCode != OBIT_IO_OK) || (err->error>0)) 
    Obit_traceback_msg (err, routine, inOTF->name);

  /* use same data buffer on input and output 
     so don't assign buffer for output */
  if (outOTF->buffer) ObitIOFreeBuffer(outOTF->buffer); /* free existing */
  outOTF->buffer     = inOTF->buffer;
  outOTF->bufferSize = inOTF->bufferSize;

 /* Open Output Data */
  retCode = ObitOTFOpen (outOTF, OBIT_IO_WriteOnly, err) ;
  if ((retCode != OBIT_IO_OK) || (err->error>0)) {
    outOTF->buffer = NULL; /* remove pointer to inOTF buffer */
    outOTF->bufferSize = 0;
    Obit_traceback_msg (err, routine, outOTF->name);
  }

 /* Copy tables before data */
  retCode = ObitOTFCopyTables (inOTF, outOTF, exclude, NULL, err);
  if (err->error) {/* add traceback,return */
    outOTF->buffer = NULL;
    outOTF->bufferSize = 0;
    Obit_traceback_msg (err, routine, inOTF->name);
  }

  /* Close and reopen input to init calibration which will have been disturbed 
     by the table copy */
  retCode = ObitOTFClose (inOTF, err);
  if (err->error) {
    outOTF->buffer = NULL; outOTF->bufferSize = 0;
    Obit_traceback_msg (err, routine, inOTF->name);
  }

  retCode = ObitOTFOpen (inOTF, access, err);
  if ((retCode != OBIT_IO_OK) || (err->error>0)) {
    outOTF->buffer = NULL; outOTF->bufferSize = 0;
    Obit_traceback_msg (err, routine, inOTF->name);
  }

  ndetect = inOTF->geom->numberDetect;    /* How many detectors */

  /* Loop over data */
  done = (retCode != OBIT_IO_OK);
  while (!done) {

    /* read buffer */
    retCode = ObitOTFRead (inOTF, NULL, err);
    if (err->error) {
      outOTF->buffer = NULL; outOTF->bufferSize = 0;
      Obit_traceback_msg (err, routine, inOTF->name);
    }
    done = (retCode == OBIT_IO_EOF); /* done? */
    if (done) break;

     /* How many? */
     ndata   = inOTF->myDesc->numRecBuff;
     outOTF->myDesc->numRecBuff = ndata;

     if (ndata>0) {
       /* Modify data from this buffer */
       
       data = inOTF->buffer;  /* data pointer */

       /* Loop over buffer */
       for (i=0; i<ndata; i++) {
	 /* Loop over detectors */
	 for (j=0; j<ndetect; j++) {
	   /* Modify */
	   if (data[desc->ilocdata+j*incdatawt]!=fblank) {
	     data[desc->ilocdata+j*incdatawt] *= scale;
	     data[desc->ilocdata+j*incdatawt] += offset;
	   }	   
	 } /* end loop over detectors */
	 data += desc->lrec; /* update buffer pointer */
       } /* end  Loop over buffer */  
 
       /* Write buffer */
       retCode = ObitOTFWrite (outOTF, NULL, err);
     }
     if (err->error) {
       outOTF->buffer = NULL; outOTF->bufferSize = 0;
       Obit_traceback_msg (err, routine, outOTF->name);
     }

  } /* end data loop */
  
  /* unset output buffer (may be multiply deallocated ;'{ ) */
  outOTF->buffer = NULL;
  outOTF->bufferSize = 0;

  /* Close input */
  retCode = ObitOTFClose (inOTF, err);
  if (err->error) Obit_traceback_msg (err, routine, inOTF->name);

  /* Close output */
  retCode = ObitOTFClose (outOTF, err);
  if (err->error) Obit_traceback_msg (err, routine, outOTF->name);

} /* end ObitOTFUtilScale */

/**
 * Add Gaussian noise to an OTF
 * out = in*scale +offset +  noise (sigma)
 * \param inOTF  Input OTF 
 * \param outOTF Output OTF, must already be defined
 * \param scale  scaling factor for data
 * \param offset additive term to be applied to all data.
 * \param sigma  Standard deviation of Gaussian noise .
 * \param err    Error stack
 */
void ObitOTFUtilNoise(ObitOTF *inOTF, ObitOTF *outOTF, ofloat scale, ofloat offset,
		      ofloat sigma, ObitErr *err)
{
  ObitIOCode retCode;
  gboolean doCalSelect, done;
  ObitInfoType type;
  ObitIOAccess access;
  gint32 dim[MAXINFOELEMDIM];
  ofloat *data, fblank = ObitMagicF();
  odouble dsigma = sigma;
  olong i, j, ndetect, ndata;
  olong incdatawt;
  ObitOTFDesc* desc;
  /* Don't copy Cal and Soln or data or flag tables */
  gchar *exclude[]={"OTFSoln", "OTFCal", "OTFScanData", "OTFFlag", NULL};
  gsl_rng *ran=NULL;
  gchar *routine = "ObitOTFUtilNoise";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitOTFIsA(inOTF));
  g_assert (ObitOTFIsA(outOTF));

  /* Local pointers */
  desc = inOTF->myDesc;
  incdatawt = desc->incdatawt; /* increment in data-wt axis */

  /* Calibration wanted? */ 
  doCalSelect = FALSE;
  ObitInfoListGetTest(inOTF->info, "doCalSelect", &type, (gint32*)dim, &doCalSelect);
  if (doCalSelect) access = OBIT_IO_ReadCal;
  else access = OBIT_IO_ReadWrite;

 /* Open Input Data */
  retCode = ObitOTFOpen (inOTF, access, err);
  if ((retCode != OBIT_IO_OK) || (err->error>0)) 
    Obit_traceback_msg (err, routine, inOTF->name);

  /* use same data buffer on input and output 
     so don't assign buffer for output */
  if (outOTF->buffer) ObitIOFreeBuffer(outOTF->buffer); /* free existing */
  outOTF->buffer     = inOTF->buffer;
  outOTF->bufferSize = inOTF->bufferSize;

  /* Open Output Data */
  retCode = ObitOTFOpen (outOTF, OBIT_IO_WriteOnly, err) ;
  if ((retCode != OBIT_IO_OK) || (err->error>0)) {
    outOTF->buffer = NULL; /* remove pointer to inOTF buffer */
    outOTF->bufferSize = 0;
    Obit_traceback_msg (err, routine, outOTF->name);
  }

  /* Copy tables before data */
  retCode = ObitOTFCopyTables (inOTF, outOTF, exclude, NULL, err);
  if (err->error) {/* add traceback,return */
    outOTF->buffer = NULL;
    outOTF->bufferSize = 0;
    Obit_traceback_msg (err, routine, inOTF->name);
  }

  /* Close and reopen input to init calibration which will have been disturbed 
     by the table copy */
  retCode = ObitOTFClose (inOTF, err);
  if (err->error) {
    outOTF->buffer = NULL; outOTF->bufferSize = 0;
    Obit_traceback_msg (err, routine, inOTF->name);
  }
  
  retCode = ObitOTFOpen (inOTF, access, err);
  if ((retCode != OBIT_IO_OK) || (err->error>0)) {
    outOTF->buffer = NULL; outOTF->bufferSize = 0;
    Obit_traceback_msg (err, routine, inOTF->name);
  }

  ndetect = inOTF->geom->numberDetect;    /* How many detectors */

  /* Init random number generator */
  ran = gsl_rng_alloc(gsl_rng_taus);
  
  /* Loop over data */
  done = (retCode != OBIT_IO_OK);
  while (!done) {
    
    /* read buffer */
    retCode = ObitOTFRead (inOTF, NULL, err);
    if (err->error) {
      outOTF->buffer = NULL; outOTF->bufferSize = 0;
      Obit_traceback_msg (err, routine, inOTF->name);
    }
    done = (retCode == OBIT_IO_EOF); /* done? */
    if (done) break;

     /* How many? */
     ndata   = inOTF->myDesc->numRecBuff;
     outOTF->myDesc->numRecBuff = ndata;

     if (ndata>0) {
       /* Modify data from this buffer */
       
       data = inOTF->buffer;  /* data pointer */

       /* Loop over buffer */
       for (i=0; i<ndata; i++) {
	 /* Loop over detectors */
	 for (j=0; j<ndetect; j++) {
	   /* Modify */
	   if (data[desc->ilocdata+j*incdatawt]!=fblank) {
	     data[desc->ilocdata+j*incdatawt] *= scale;
	     data[desc->ilocdata+j*incdatawt] += offset + 
	       (ofloat)gsl_ran_gaussian (ran, dsigma);
	   }	   
	 } /* end loop over detectors */
	 data += desc->lrec; /* update buffer pointer */
       } /* end  Loop over buffer */  
 
       /* Write buffer */
       retCode = ObitOTFWrite (outOTF, NULL, err);
     }
     if (err->error) {
       outOTF->buffer = NULL; outOTF->bufferSize = 0;
       Obit_traceback_msg (err, routine, outOTF->name);
     }

  } /* end data loop */
  
  /* unset output buffer (may be multiply deallocated ;'{ ) */
  outOTF->buffer = NULL;
  outOTF->bufferSize = 0;

  /* Free random number generator */
  gsl_rng_free(ran);
  
  /* Close input */
  retCode = ObitOTFClose (inOTF, err);
  if (err->error) Obit_traceback_msg (err, routine, inOTF->name);

  /* Close output */
  retCode = ObitOTFClose (outOTF, err);
  if (err->error) Obit_traceback_msg (err, routine, outOTF->name);

} /* end ObitOTFUtilNoise */

/**
 * Subtract a sky model from a buffer full of OTF data.
 * \param in         OTF with internal buffer to be modified.
 * \param sky        OTFSkyModel to subtract.
 * \param factor     Scaling factor for sky model.
 */
void ObitOTFUtilSubSkyModelBuff (ObitOTF *in, ObitOTFSkyModel *sky, ofloat factor)
{
  ofloat *xpos, *ypos, *data, dist, maxdist, arg, gf1, gf2, sigma;
  olong incdatawt;
  olong  ndetect, ndata, ncomp, i, j, k;
  ObitOTFDesc* desc;
  ObitOTFArrayGeom* geom;

  /* Error checks */
  g_assert (ObitOTFIsA(in));
  g_assert (ObitOTFSkyModelIsA(sky));

  /* Local pointers */
  desc = in->myDesc;
  geom = in->geom;
  data = in->buffer;
  ndetect = geom->numberDetect;    /* How many detectors */
  ndata   = desc->numRecBuff;      /* How many data records */
  ncomp   = sky->numberComp;        /* How many components */
  incdatawt = desc->incdatawt; /* increment in data-wt axis */

  /* Allocate temporary arrays */
  /* Projected data locations */
  xpos = g_malloc0(ndetect*sizeof(float));
  ypos = g_malloc0(ndetect*sizeof(float));

  /* How close (deg^2) to consider 3 * beam size */
  maxdist = (desc->beamSize * 3.0) * (desc->beamSize * 3.0);

  /* Gaussian factors */
  sigma = desc->beamSize/2.35;
  gf1 = 1.0 / (2.0 * sigma * sigma);
  /* Normalize the flux integral rather than the peak */
  gf2 = factor / sqrt(2.0 * G_PI);
  /* debug - normalize peak */
  gf2 = -1.0;

  /* Loop over data records */
  for (i=0; i<ndata; i++) {

    /* Get Sky locations of the data */
    ObitOTFArrayGeomProj (geom, data[desc->ilocra], data[desc->ilocdec], 
			  data[desc->ilocrot], sky->RACenter, sky->DecCenter, 
			  sky->proj, xpos, ypos);
    
    /* Loop over the array */
    for (j=0; j<ndetect; j++) {

      /* loop over components in Sky model */
      for (k=0; k<ncomp; k++) {
	/* how far is the measurement from the component (distance^2) */
	dist = (xpos[j]-sky->RAOffset[k]) * (xpos[j]-sky->RAOffset[k]) + 
	  (ypos[j]-sky->DecOffset[k]) * (ypos[j]-sky->DecOffset[k]);

	/* Is this one close enough? */
	if (dist<maxdist) {
	  arg = -dist * gf1;
	  data[desc->ilocdata+j*incdatawt] -= gf2 * exp (arg) * sky->flux[k];  /* subtract */
	}
      } /* end loop over components */
    } /* end loop over array */
    data += desc->lrec; /* update buffer pointer */
  } /* end loop over buffer */

  /* cleanup */
  if (xpos) g_free(xpos);
  if (ypos) g_free(ypos);
} /* end ObitOTFUtilSubSkyModelBuff */

/**
 * Subtract the values in an image from a buffer full of OTF data.
 * For CCB beamswitched data (OTFType=OBIT_GBTOTF_CCB, no. States=1)
 * If threading has been enabled by a call to ObitThreadAllowThreads 
 * this routine will divide the buffer up amount the number of processors
 * returned by ObitThreadNumProc.
 * \param in         OTF with internal buffer to be modified.
 * \param image      Image interpolator
 * \param factor     Scaling factor for model.
 * \param nThreads   Number of elements in args
 * \param args       Threaded function argument structs
 * \param err        Error stack
 */
void ObitOTFUtilSubImageBuff (ObitOTF *in, ObitFInterpolate *image, 
			      ofloat factor, olong nThreads, SubImageFuncArg **args, 
			      ObitErr *err)
{
  olong i, nrec, lorec, hirec, nTh, nrecPerThread;
  gboolean OK = TRUE;
  gchar *routine = "ObitOTFUtilSubImageBuff";

  /* error checks - assume most done at higher level */
  if (err->error) return;

  /* Divide up work */
  nrec = in->myDesc->numRecBuff;
  nrecPerThread = nrec/nThreads;
  nTh = nThreads;
  if (nrec<100) {nrecPerThread = nrec; nTh = 1;}
  lorec = 1;
  hirec = nrecPerThread;
  hirec = MIN (hirec, nrec);

  /* Set up thread arguments */
  for (i=0; i<nTh; i++) {
    if (i==(nTh-1)) hirec = nrec;  /* Make sure do all */
    args[i]->otfdata = in;
    args[i]->first  = lorec;
    args[i]->last   = hirec;
    if (i==0) args[i]->Interp = ObitFInterpolateRef(image);
    else args[i]->Interp = ObitFInterpolateClone(image, NULL);
    args[i]->factor = factor;
    /* Update which rec */
    lorec += nrecPerThread;
    hirec += nrecPerThread;
    hirec = MIN (hirec, nrec);
  }

  /* Do operation */
  OK = ObitThreadIterator (in->thread, nTh, 
			   (ObitThreadFunc)ThreadOTFUtilSubImageBuff, 
			   (gpointer**)args);

  /* Check for problems */
  if (!OK) Obit_log_error(err, OBIT_Error,"%s: Problem in threading", routine);
} /* end ObitOTFUtilSubImageBuff */

/**
 * Replace data with the values in an image in a buffer full of OTF data.
 * For CCB beamswitched data (OTFType=OBIT_GBTOTF_CCB, no. States=1)
 * \param in         OTF with internal buffer to be modified.
 * \param image      Image interpolator
 * \param factor     Scaling factor for model.
 * \param err        Error stack
 */
void ObitOTFUtilModelImageBuff (ObitOTF *in, ObitFInterpolate *image, 
				ofloat factor, ObitErr *err)
{
  ofloat *xpos, *ypos, *data, ffact, value, RACenter, DecCenter;
  ObitOTFProj proj;
  odouble coord[IM_MAXDIM];
  olong  ndetect, ndata, i, j;
  olong incfeed, itemp,incdatawt ;
  gboolean CCBBS, isRef;
  ObitOTFDesc* desc;
  ObitOTFArrayGeom* geom;
  ofloat fblank = ObitMagicF();
  gchar Proj[5];

  /* Error checks */
  g_assert (ObitOTFIsA(in));
  g_assert (ObitFInterpolateIsA(image));

  /* Local pointers */
  desc = in->myDesc;
  geom = in->geom;
  data = in->buffer;
  ndetect = geom->numberDetect;    /* How many detectors */
  ndata   = desc->numRecBuff;      /* How many data records */
  if (in->myDesc->jlocfeed>=0) 
    incfeed = in->myDesc->incfeed / in->myDesc->incdatawt;
  else incfeed = 1;  /* This is probably a bad sign */

  /* Get Model center */
  RACenter  = image->myDesc->crval[0];
  DecCenter = image->myDesc->crval[1];
  strncpy (Proj, &image->myDesc->ctype[0][4], 5);
  proj      = ObitOTFSkyModelProj (Proj);
  incdatawt = desc->incdatawt; /* increment in data-wt axis */

  /* Is this CCB beamswitched data? */
  CCBBS = (in->myDesc->OTFType==OBIT_GBTOTF_CCB) &&
    (in->myDesc->jlocstate>=0) &&
    (in->myDesc->inaxes[in->myDesc->jlocstate]==1);

  /* Allocate temporary arrays */
  /* Projected data locations */
  xpos = g_malloc0(ndetect*sizeof(float));
  ypos = g_malloc0(ndetect*sizeof(float));
  ffact = factor;

  /* Loop over data records */
  for (i=0; i<ndata; i++) {

    /* Get Sky locations of the data projected onto the image */
    ObitOTFArrayGeomProj(geom, data[desc->ilocra], data[desc->ilocdec], 
			 data[desc->ilocrot], RACenter, DecCenter, proj, 
			 xpos, ypos);

    /* Loop over the array */
    for (j=0; j<ndetect; j++) {
      
     /* Interpolate - use coordinates on a flat plane at the tangent point of the image 
	 xpos and ypos are offsets from the image center projected onto the plane 
	 of the image, ObitFInterpolateOffset does a linear approximation. */
       coord[0] = xpos[j];  /* + RACenter; */
       coord[1] = ypos[j];  /* + DecCenter; */
       value = ObitFInterpolateOffset (image, coord, err); 

      /* Replace */
      if ((value!=fblank) && (data[desc->ilocdata+j*incdatawt]!=fblank))
	data[desc->ilocdata+j*incdatawt] = value * ffact;
      else {
	data[desc->ilocdata+j*incdatawt] = fblank;
      }

      /* For CCB beamswitched data add value corresponding to this feeds
       reference beam */
      if (CCBBS) {
	/* is this the reference or signal beam? */
	itemp = j / incfeed;
	/* itemp odd means reference beam */
	isRef = itemp != 2*(itemp/2);
	/* Use feed offset for feed 1 if this is reference, else feed 2 */
	if (isRef) itemp = 0;
	else itemp = incfeed;
	/* interpolate */
	coord[0] = xpos[itemp]; /* + RACenter; */
	coord[1] = ypos[itemp]; /* + DecCenter; */
	value = ObitFInterpolateOffset (image, coord, err);
       /* Replace */
       if ((value!=fblank) && (data[desc->ilocdata+j*incdatawt]!=fblank))
	 data[desc->ilocdata+j*incdatawt] = value * ffact;
       else {
	 data[desc->ilocdata+j*incdatawt] = fblank;
       }
      }  /* end adding to reference beam position */
      
    } /* end loop over array */
    data += desc->lrec; /* update buffer pointer */
  } /* end loop over buffer */
  
  /* cleanup */
  if (xpos) g_free(xpos);
  if (ypos) g_free(ypos);
} /* end ObitOTFUtilModelImageBuff */

/**
 * Create basic ObitImage structure and fill out descriptor.
 * Imaging parameters are on the inOTF info member.
 * \li "nx"     OBIT_int (1,1,1) Dimension of image in RA [no default].
 * \li "ny"     OBIT_int (1,1,1) Dimension of image in declination[no default]
 * \li "RA"     OBIT_float (1,1,1) Right Ascension of center of image
 *                                 Default is observed position center in inOTF 
 * \li "Dec"    OBIT_float (1,1,1) Declination of center of image
 *                                 Default is observed position center in inOTF 
 * \li "xCells" OBIT_float (1,1,1) X (=RA) cell spacing in asec [no default]
 * \li "yCells" OBIT_float (1,1,1) Y (=dec) cell spacing in asec [no default]
 * \li "Proj"   OBIT_string (4,1,1) Projection string "-SIN", "-ARC", "-TAN"
 *                         [Default "-SIN"]
 * \param inOTF     Input OTF data. 
 * \param err      Error stack, returns if not empty.
 * \return Pointer to the newly created ObitImage.
 */
ObitImage* ObitOTFUtilCreateImage (ObitOTF *inOTF, ObitErr *err)
{
  ObitImage *outImage=NULL;
  gchar outName[121], Proj[8];
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  olong itemp[10], nx, ny;
  ofloat ftemp[10], xCells, yCells, RA, Dec;
  gchar *routine = "ObitOTFUtilCreateImage";
 
   /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return outImage;
  g_assert (ObitOTFIsA(inOTF));

  /* open/close OTF data to fully instantiate if not already open */
  if (inOTF->myStatus==OBIT_Inactive) {
    ObitOTFFullInstantiate (inOTF, TRUE, err);
    if (err->error) Obit_traceback_val (err, routine, inOTF->name, outImage);
  }

  /* Create output */
  g_snprintf (outName, 120, "%s",inOTF->name);
  outImage = newObitImage(outName);

  /* Get parameters for image */
  /* Image size */
  itemp[0] = -1;
  if (!ObitInfoListGetTest(inOTF->info, "nx", &type, dim, itemp))
    Obit_log_error(err, OBIT_Error, "%s: %s MUST define nx", routine, inOTF->name);
  nx = itemp[0];
 
  if (!ObitInfoListGetTest(inOTF->info, "ny", &type, dim, itemp))
    Obit_log_error(err, OBIT_Error, "%s: %s MUST define ny", routine, inOTF->name);
  ny = itemp[0];

  /* Cell Spacing */
  ftemp[0] = 0.0;
  if (!ObitInfoListGetTest(inOTF->info, "xCells", &type, dim, ftemp))
    Obit_log_error(err, OBIT_Error, "%s: %s MUST define xCells", routine, inOTF->name);
  xCells = ftemp[0];

  if (!ObitInfoListGetTest(inOTF->info, "yCells", &type, dim, ftemp)) 
    Obit_log_error(err, OBIT_Error, "%s: %s MUST define yCells", routine, inOTF->name);
  yCells = ftemp[0];

  /* Center Default is observed position center in inOTF */
  ftemp[0] = inOTF->myDesc->obsra;
  ObitInfoListGetTest(inOTF->info, "RA", &type, dim, ftemp);
  RA = ftemp[0];

  ftemp[0] = inOTF->myDesc->obsdec;
  ObitInfoListGetTest(inOTF->info, "Dec", &type, dim, ftemp);
  Dec = ftemp[0];
 
  /* Default projection is -SIN */
  Proj[0] = '-'; Proj[1] = 'S'; Proj[2] = 'I'; Proj[3] = 'N'; Proj[4] = 0; 
  ObitInfoListGetTest(inOTF->info, "Proj", &type, dim, (gpointer*)&Proj);
  Proj[4] = 0;
  
  /* bail out if an error so far */
  if (err->error) return outImage;

  /* Set values on descriptor(s) */
  outImage->myDesc->xshift = 0.0;
  outImage->myDesc->yshift = 0.0;
  outImage->myDesc->crota[0] = 0.0;
  outImage->myDesc->crota[1] = 0.0;
  outImage->myDesc->cdelt[0] = xCells / 3600.0;
  outImage->myDesc->cdelt[1] = yCells / 3600.0;
  outImage->myDesc->inaxes[0] = nx;
  outImage->myDesc->inaxes[1] = ny;
  outImage->myDesc->crval[0] = RA;
  outImage->myDesc->crval[1] = Dec;

  /* Fill in descriptor */
  ObitOTFUtilOTF2ImageDesc (inOTF->myDesc, outImage->myDesc, Proj);

  return outImage;
} /* end ObitOTFUtilCreateImage */

/**
 * Convolves data onto a grid, accumulated and normalizes.
 * Imaging parameters are on the inOTF info member.
 * \li "ConvType"   OBIT_long scalar = Convolving function type: [def=3]
 *                  0 = pillbox, 3 = Gaussian, 4 = Exp*Sinc, 5 = Spherodial wave
 * \li "ConvParm"  OBIT_float[10] = Convolving function parameters
 * \li "deBias"    OBIT_bool scalar = Subtract calibration bias from image? [def False]
 *                 Note, this doesn't really work the way you would like
 * \li "deMode"    OBIT_bool scalar = Subtract image mode from image? [def False]
 * \li "minWt"     OBIT_float (1,1,1) Minimum summed gridding convolution weight 
 *                 as a fraction of the maximum [def 0.01]
 * \li "doScale"   OBIT_bool scalar If true, convolve/scale beam [def TRUE]
 * \li "doFilter"  OBIT_bool scalar If true, filter out of band noise[def TRUE]
 * \param inOTF    Input OTF data. 
 * \param outImage Image to be written.  Must be previously instantiated.
 * \param doBeam   If TRUE also make convolved beam.  
 *                 Will make the myBeam member of outImage.
 * \param Beam     If non NULL use as instrumental response beam 
 * \param Wt       If non NULL write weight array to.
 * \param err      Error stack, returns if not empty.
 */
void ObitOTFUtilMakeImage (ObitOTF *inOTF, ObitImage *outImage, gboolean doBeam, 
			   ObitImage *Beam, ObitImage *Wt, ObitErr *err)
{
  ObitIOSize IOBy;
  gchar *outName=NULL;
  ObitOTFGrid *myGrid=NULL;
  ObitFArray *biasArray=NULL;  
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong blc[IM_MAXDIM] = {1,1,1,1,1};
  olong trc[IM_MAXDIM] = {0,0,0,0,0};
  ofloat parms[10], xparms[10], mode, radius;
  olong i, nparms, cType, tcType;
  gboolean deBias, replCal, tbool, deMode, doFilter;
  gchar *parmList[] = {"minWt","Clip","beamNx","beamNy","doScale",NULL};
  gchar *routine = "ObitOTFUtilMakeImage";

   /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitOTFIsA(inOTF));
  g_assert (ObitImageIsA(outImage));

 /* Make Gridding object */
  outName = g_strconcat ("OTFGrid for: ", inOTF->name, NULL);
  myGrid = newObitOTFGrid(outName);
  g_free(outName);

  /* Copy control parameters to grid object */
  ObitInfoListCopyList (inOTF->info, myGrid->info, parmList);

  /* Now make image */
  /* Set blc, trc */
  blc[0] = 1;
  blc[1] = 1;
  blc[2] = 1;
  trc[0] = outImage->myDesc->inaxes[0];
  trc[1] = outImage->myDesc->inaxes[0];
  trc[2] = 1;
      
  IOBy = OBIT_IO_byPlane;
  dim[0] = 1;
  ObitInfoListPut (outImage->info, "IOBy", OBIT_long, dim, (gpointer)&IOBy, err);
  dim[0] = 7;
  ObitInfoListPut (outImage->info, "BLC", OBIT_long, dim, (gpointer)blc, err); 
  ObitInfoListPut (outImage->info, "TRC", OBIT_long, dim, (gpointer)trc, err);
  
  /*  Open image to get descriptor */
  if ((ObitImageOpen (outImage, OBIT_IO_WriteOnly, err) 
       != OBIT_IO_OK) || (err->error>0))
    Obit_log_error(err, OBIT_Error, "ERROR opening image %s", outImage->name);
  if (err->error) return;
  
  /* reset max, min */
  outImage->myDesc->minval =  1.0e20;
  outImage->myDesc->maxval = -1.0e20;
  ((ObitImageDesc*)outImage->myIO->myDesc)->minval =  1.0e20;
  ((ObitImageDesc*)outImage->myIO->myDesc)->maxval = -1.0e20;
  
  /* Gaussian beam */
  outImage->myDesc->beamMaj = inOTF->myDesc->beamSize;
  outImage->myDesc->beamMin = inOTF->myDesc->beamSize;
  outImage->myDesc->beamPA  = 0.0;
  
  /* Gridding setup */
  ObitOTFGridSetup (myGrid, inOTF, outImage->myDesc, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);

  /* Beam if requested */
  if (doBeam) ObitOTFGridMakeBeam (myGrid, outImage, Beam, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);
  
  /* Grid */
  ObitOTFGridReadOTF (myGrid, inOTF, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);
  
  /* Normalize */
  ObitOTFGridNorm(myGrid, outImage->image, outImage->myDesc, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);

   /* reset beam size on IO descriptor for the effects of convolution */
  ((ObitImageDesc*)outImage->myIO->myDesc)->beamMaj = outImage->myDesc->beamMaj;
  ((ObitImageDesc*)outImage->myIO->myDesc)->beamMin = outImage->myDesc->beamMin;
  
  /* Remove image Mode */
  deMode = FALSE;
  ObitInfoListGetTest(inOTF->info, "deMode", &type, dim, &deMode);
  if (deMode) {
     /* remove mode of image */
    mode = ObitFArrayMode(outImage->image);
    ObitFArraySAdd(outImage->image, -mode);
    /* Message */
    Obit_log_error(err, OBIT_InfoErr, "Subtract image Mode %g",mode);
  }

  /* Debias? */
  deBias = FALSE;
  ObitInfoListGetTest(inOTF->info, "deBias", &type, dim, &deBias);
  if (deBias) {
    /* Save image array */
    biasArray = ObitFArrayCopy (outImage->image, biasArray, err);
    if (err->error) Obit_traceback_msg (err, routine, outImage->name);

    /* This time replace data with cal value */
    replCal = FALSE;
    ObitInfoListGetTest(inOTF->info, "replCal", &type, dim, &replCal);
    tbool = TRUE;
    ObitInfoListAlwaysPut(inOTF->info, "replCal", OBIT_bool, dim, &tbool);

    /* Use large Gaussian convolving fn twice beam size*/
    cType = 3;
    ObitInfoListGetTest(inOTF->info, "ConvType", &type, dim, &cType);
    tcType = 3;
    ObitInfoListAlwaysPut(inOTF->info, "ConvType", OBIT_long, dim, &tcType);
    for (i=0; i<10; i++) parms[i] = 0.0;  /* Default 0.0 */
    dim[0] = 1;
    ObitInfoListGetTest(inOTF->info, "ConvParm", &type, dim, parms);
    nparms = dim[0];
    for (i=0; i<10; i++) xparms[i] = parms[i]; 
    xparms[0] = 7.0; xparms[1] = 2.0;
    ObitInfoListAlwaysPut(inOTF->info, "ConvParm", OBIT_float, dim, xparms);
  
    /* Gridding setup */
    ObitOTFGridSetup (myGrid, inOTF, outImage->myDesc, err);
    if (err->error) Obit_traceback_msg (err, routine, outImage->name);

    /* reGrid */
    ObitOTFGridReadOTF (myGrid, inOTF, err);
    if (err->error) Obit_traceback_msg (err, routine, outImage->name);
    
    /* Normalize  bias image */
    ObitOTFGridNorm(myGrid, outImage->image, outImage->myDesc, err);
    if (err->error) Obit_traceback_msg (err, routine, outImage->name);

    /* Subtract bias image array */
    ObitFArraySub(biasArray, outImage->image, outImage->image);

    /* remove mode of image */
    mode = ObitFArrayMode(outImage->image);
    ObitFArraySAdd(outImage->image, -mode);

    /* Restore the state of things */
    dim[0] = 1;
    ObitInfoListAlwaysPut(inOTF->info, "replCal", OBIT_bool, dim, &replCal);
    ObitInfoListAlwaysPut(inOTF->info, "ConvType", OBIT_long, dim, &cType);
    dim[0] = nparms;
    ObitInfoListAlwaysPut(inOTF->info, "ConvParm", OBIT_float, dim, parms);
    if (biasArray) ObitFArrayUnref(biasArray);
  } /* end debias */

 /* Write image */
  ObitImageWrite (outImage, NULL, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);
  
  /* tell Max/Min */
  Obit_log_error(err, OBIT_InfoErr, 
		 "Image max %g, min %g for %s", 
		 outImage->myDesc->maxval, outImage->myDesc->minval, 
		 outImage->name);
  
  /* Close Image */
  ObitImageClose (outImage, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);

  /* Filter if requested */
  doFilter = TRUE;
  ObitInfoListGetTest(inOTF->info, "doFilter", &type, dim, &doFilter);
  radius = 0.5* inOTF->myDesc->diameter;
  if (radius<=0.0) radius = 50.0;  /* Default = GBT */
  if (doFilter) 
    ObitImageUtilUVFilter(outImage, outImage, radius, err);
  if (err->error) Obit_traceback_msg (err, routine, outImage->name);

  /* Save Weight image? */
  if (Wt!=NULL) {
    ObitImageClone (outImage, Wt, err);   /* Looks like outImage */
    ObitImageOpen (Wt, OBIT_IO_WriteOnly, err);
    /* reset max, min */
    Wt->myDesc->minval =  1.0e20;
    Wt->myDesc->maxval = -1.0e20;
    ((ObitImageDesc*)Wt->myIO->myDesc)->minval =  1.0e20;
    ((ObitImageDesc*)Wt->myIO->myDesc)->maxval = -1.0e20;
    ObitImageWrite (Wt, myGrid->gridWt->array, err);
    ObitImageClose (Wt, err);
    if (err->error) Obit_traceback_msg (err, routine, Wt->name);
  }

  /* Free myGrid */
  myGrid = ObitOTFGridUnref(myGrid);
  
}  /* end ObitOTFUtilMakeImage */

/**
 * Reads the OTF and rewrites its OTFIndex table
 * \param inOTF    Input OTF data. 
 * \param err      Error stack, returns if not empty.
 */
void ObitOTFUtilIndex (ObitOTF *inOTF, ObitErr *err)
{ 
  ObitIOCode retCode;
  ObitTableOTFIndex* table;
  ObitTableOTFIndexRow* row;
  olong num, i, lrec, iRow, ver, lastscan, iscan, target, lastTarget=0, startRec=1, curRec;
  ofloat *rec;
  odouble startTime=0.0, endTime=0.0; 
  gchar *routine = "ObitOTFUtilIndex";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitOTFIsA(inOTF));

  /* Open OTF */
  retCode = ObitOTFOpen (inOTF, OBIT_IO_ReadWrite, err);
  if ((retCode != OBIT_IO_OK) || (err->error>0)) 
    Obit_traceback_msg (err, routine, inOTF->name);
  lrec = inOTF->myDesc->lrec;  /* Size of record */
  lastscan = -1000; /* initialize scan number */
  curRec = 0;       /* record counter */

  /* create Index table object */
  ver = 1;
  table = newObitTableOTFIndexValue ("Index table", (ObitData*)inOTF, &ver, 
				     OBIT_IO_ReadWrite, err);
  if (err->error) Obit_traceback_msg (err, routine, inOTF->name);

  /* Open Index table */
  if ((ObitTableOTFIndexOpen (table, OBIT_IO_ReadWrite, err) 
       != OBIT_IO_OK) || (err->error))  { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR opening output OTFIndex table");
    return;
  }

  /* Create Index Row */
  row = newObitTableOTFIndexRow (table);

  /* initialize row */
  row->ScanID = 0;
  row->TargetID = 0;
  row->Time = 0.0;
  row->TimeI = 0.0;
  row->StartRec = -1;
  row->EndRec = -1;

  /* attach to table buffer */
  ObitTableOTFIndexSetRow (table, row, err);
  if (err->error) Obit_traceback_msg (err, routine, inOTF->name);

  /* Write at beginning of OTFIndex Table */
  iRow = 0;
  table->myDesc->nrow = 0; /* ignore any previous entries */

  /* Loop over OTF */
  while (retCode==OBIT_IO_OK) {
    retCode = ObitOTFRead (inOTF, inOTF->buffer, err);
    if (retCode!=OBIT_IO_OK) break;

    /* How many */
    num = inOTF->myDesc->numRecBuff;
    
    /* Record pointer */
    rec = inOTF->buffer;

    /* initialize on first record */
    if (curRec<=0) {
      startRec   = 1;
      startTime  = rec[inOTF->myDesc->iloct];
    }
    
    /* Loop over buffer */
    for (i=0; i<num; i++) {
      
      iscan  = rec[inOTF->myDesc->ilocscan] + 0.5; /* Which scan number */
      target = rec[inOTF->myDesc->iloctar] + 0.5;  /* Target number */
      curRec++; /* Current OTF record number */
      /* Initialize? */
      if (lastscan<=0)   lastscan = iscan;
      if (lastTarget<=0) lastTarget = target;
      
      
      /* New scan? */
      if ((iscan!=lastscan) && (lastscan>0)) {  /* Write index */

	/* Record values */
	row->ScanID   = lastscan;
	row->TargetID = lastTarget;
	row->Time     = 0.5 * (startTime + endTime);
	row->TimeI    = (endTime - startTime);
	row->StartRec = startRec;
	row->EndRec   = curRec-1;
	
	/* Write OTFIndex table */
	iRow++;
	if ((ObitTableOTFIndexWriteRow (table, iRow, row, err)
	     != OBIT_IO_OK) || (err->error>0)) { 
	  Obit_log_error(err, OBIT_Error, "ERROR writing OTFIndex Table file");
	  return;
	}

	/* Initialize next scan */
	lastscan   = iscan;
	lastTarget = target;
	startRec   = curRec;
	startTime  = rec[inOTF->myDesc->iloct];
	endTime    = rec[inOTF->myDesc->iloct];
	
      } /* end of if new scan */

      endTime    = rec[inOTF->myDesc->iloct]; /* potential end time */
      rec += inOTF->myDesc->lrec;     /* Update data record pointer */
    } /* end loop over buffer */

  } /* End loop over OTF */

  /* Last Scan */
  /* Record values */
  row->ScanID   = lastscan;
  row->TargetID = lastTarget;
  row->Time     = 0.5 * (startTime + endTime);
  row->TimeI    = (endTime - startTime);
  row->StartRec = startRec;
  row->EndRec   = curRec-1;
  
  /* Write OTFIndex table */
  iRow++;
  if ((ObitTableOTFIndexWriteRow (table, iRow, row, err)
       != OBIT_IO_OK) || (err->error>0)) { 
    Obit_log_error(err, OBIT_Error, "ERROR writing OTFIndex Table file");
    return;
  }

 /* Close OTFIndex table */
  if ((ObitTableOTFIndexClose (table, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR closing output OTFIndex Table file");
    return;
  }

  /* Cleanup */
  row = ObitTableOTFIndexRowUnref(row);
  table = ObitTableOTFIndexUnref(table);

  /* Close OTF */
  retCode = ObitOTFClose (inOTF, err);
  if ((retCode != OBIT_IO_OK) || (err->error>0)) 
    Obit_traceback_msg (err, routine, inOTF->name);
} /* end ObitOTFUtilIndex */

/**
 * Differences the Ons and Offs in a beamswitched nodding scan
 * Output values on inOTF
 * \li "OnOff"    OBIT_float (*,1,1) Differences of On-Off for each detector
 *                in the same order as defined in the data.  For beamswitched
 *                data this will be twice the source strength.
 * \param inOTF    Input OTF data. Applies any calibration specified
 *                 Target position must be in OTFTarget table.
 * \param scan     Scan number
 * \param err      Error stack, returns if not empty.
 */
void ObitOTFUtilDiffNod (ObitOTF *inOTF, olong scan, ObitErr *err)
{
  ObitIOCode retCode;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitTableOTFTarget* targetTable=NULL;
  ofloat *avgOff=NULL, *avgOn=NULL, *feedRA=NULL, *feedDec=NULL;
  olong   *cntOn=NULL, *cntOff=NULL;
  olong i, state, ndetect, incfeed, itemp, iDet;
  olong scans[2], doCal, incdatawt, targID=0;
  olong ver;
  gboolean doCalSelect, isCal, isRef, gotTarinfo=FALSE;
  odouble RACal, DecCal, dra, ddec, val;
  ofloat *rec, FluxCal, fblank = ObitMagicF();
  gchar *routine = "ObitOTFUtilDiffNod";

   /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitOTFIsA(inOTF));

  /* How many detectors? */
  ndetect = inOTF->geom->numberDetect;

  /* Create arrays */
  avgOff  = g_malloc0(ndetect*sizeof(ofloat));
  avgOn   = g_malloc0(ndetect*sizeof(ofloat));
  cntOff  = g_malloc0(ndetect*sizeof(olong));
  cntOn   = g_malloc0(ndetect*sizeof(olong));
  feedRA  = g_malloc(ndetect*sizeof(ofloat));
  feedDec = g_malloc(ndetect*sizeof(ofloat));
  for (i=0; i<ndetect; i++) avgOff[i] = avgOn[i] = 0.0;
  for (i=0; i<ndetect; i++) cntOff[i] = cntOn[i] = 0;

  /* Select scan on input */
  doCalSelect = TRUE;
  dim[0] = 1;
  ObitInfoListAlwaysPut(inOTF->info, "doCalSelect", OBIT_bool, dim, &doCalSelect);
  doCal = 1;
  ObitInfoListAlwaysPut(inOTF->info, "doCalib", OBIT_bool, dim, &doCal);
  scans[0] = scan; scans[1] = scan;
  dim[0] = 2;
  ObitInfoListAlwaysPut(inOTF->info, "Scans", OBIT_long, dim, scans);
  incdatawt = inOTF->myDesc->incdatawt; /* increment in data-wt axis */

   /* open OTF data to fully instantiate  */
  retCode = ObitOTFOpen (inOTF, OBIT_IO_ReadCal, err);
  if (err->error) goto cleanup;

  /* loop reading data */
  retCode = OBIT_IO_OK;
  while (retCode == OBIT_IO_OK) {

    /* read buffer */
    retCode = ObitOTFReadSelect (inOTF, NULL, err);
    if (err->error) goto cleanup;
    if (retCode==OBIT_IO_EOF) break; /* done? */

    /* Record pointer */
    rec = inOTF->buffer;
  
    /* Feed increment in data */
    if (inOTF->myDesc->jlocfeed>=0) 
      incfeed = inOTF->myDesc->incfeed / inOTF->myDesc->incdatawt;
    else incfeed = 1;  /* This is probably a bad sign */

    /* Loop over buffer */
       for (i=0; i<inOTF->myDesc->numRecBuff; i++) {

	 /* Get source info on first record */
	 if (!gotTarinfo) {
	   gotTarinfo = TRUE;
	   /* Get position from table */
	   ver = 1;
	   targID     = (olong)rec[inOTF->myDesc->iloctar];
	   targetTable = 
	     newObitTableOTFTargetValue ("TargetTable", (ObitData*)inOTF, &ver, OBIT_IO_ReadWrite, 
					 err);
	   ObitTableOTFTargetGetSource (targetTable, targID, &RACal, &DecCal, &FluxCal, err);
	   targetTable = ObitTableOTFTargetUnref(targetTable);
	   if (err->error) goto cleanup;
	   
	   /* Make sure there is something */
	   if ((RACal==0.0) || (DecCal==0.0)) {
	     Obit_log_error(err, OBIT_Error, "%s: MISSING Calibrator info for %d %lf %lf in %s", 
			    routine, targID, RACal, DecCal, inOTF->name);
	     goto cleanup;
	   }

	 } /* End of get target info & create arrays */

	 /* Get feed positions */
	 ObitOTFArrayGeomCoord (inOTF->geom,  rec[inOTF->myDesc->ilocra], 
				rec[inOTF->myDesc->ilocdec], rec[inOTF->myDesc->ilocrot], 
				feedRA, feedDec);
	 /* Three states here, sig beam on source (state=1), ref beam on source (state=-1) 
	    or neither (state=0) */
	 state = 0;
	 /* Close to sig beam?  */
	 itemp = 0;
	 dra  = feedRA[itemp]  - RACal;
	 ddec = feedDec[itemp] - DecCal;
	 if (sqrt(dra*dra+ddec*ddec) < 0.3*inOTF->myDesc->beamSize) state = 1;

	 if (!state) {
	   /* Close to reference beam? */
	   itemp = inOTF->myDesc->incfeed;
	   dra  = feedRA[itemp]  - RACal;
	   ddec = feedDec[itemp] - DecCal;
	   if (sqrt(dra*dra+ddec*ddec) < 0.3*inOTF->myDesc->beamSize) state = -1;
	 }

	 /* sum values if on sig or ref position and cal off */
	 isCal = rec[inOTF->myDesc->iloccal]!=0.0; /* Cal on? */
	 if ((!isCal) && (state!=0)) {
	   for (iDet=0; iDet<=ndetect; iDet++) {
	     val = rec[inOTF->myDesc->ilocdata+iDet*incdatawt];
	     if (val==fblank) continue;

	     /* Is this a sig or ref beam */
	     itemp = iDet / incfeed;
	     /* itemp odd is reference beam */
	     isRef = itemp != 2*(itemp/2);

	     if (isRef) { /* reference beam */
	       if (state<0)      {avgOn[iDet]  += val; cntOn[iDet]++;}
	       else if (state>0) {avgOff[iDet] += val; cntOff[iDet]++;}
	     } else { /* signal beam */
	       if (state>0)      {avgOn[iDet]  += val; cntOn[iDet]++;}
	       else if (state<0) {avgOff[iDet] += val; cntOff[iDet]++;}
	     }
	   }
	 }

	 rec += inOTF->myDesc->lrec; /* Data record pointer */	 
       } /* end loop over buffer load */
  } /* end loop reading data */ 
  
  /* Close data */
  retCode = ObitOTFClose (inOTF, err);
  if (err->error) goto cleanup;

  /* Get On-Off */
  for (i=0; i<ndetect; i++) {
    /* Average */
    if (cntOn[i]>0) avgOn[i] /= cntOn[i];
    else avgOn[i] = fblank;
    if (cntOff[i]>0) avgOff[i] /= cntOff[i];
    else avgOff[i] = fblank;

    if ((avgOff[i]!=fblank) && (avgOn[i]!=fblank)) {
       avgOn[i] = (avgOn[i] - avgOff[i]);
    } else {
       avgOn[i] = fblank;
    }
  } /* end loop over detectors */

  /* Save values on inOTF */
  dim[0] = ndetect; dim[1] = 1;
  ObitInfoListPut (inOTF->info, "OnOff",  OBIT_float, dim,  avgOn, err);

  /* Cleanup */
 cleanup:
  if (avgOn)   g_free(avgOn);
  if (avgOff)  g_free(avgOff);
  if (cntOn)   g_free(cntOn);
  if (cntOff)  g_free(cntOff);
  if (feedRA)  g_free(feedRA);
  if (feedDec) g_free(feedDec);
 
}  /* end ObitOTFUtilDiffNod */

/**
 * Create an image and fill the descriptor values for an image cube
 * based on  the descriptor for a single plane and for the uv data 
 * creating the image.
 * This should be called before the image is Opened or instantiated.
 * \param inDesc    Input Image Descriptor.
 * \param UVDesc    Input UV Descriptor.
 * \param outDesc   Output Image Descriptor 
 * \param Stokes    Stokes parameter of image ' '=>'I', (I, Q, U, V, R, L)
 * \param bchan     first (1-rel) channel in UVDesc
 * \param echan     highest (1-rel) channel in UVDesc
 * \param incr      channel increment in input
 * \param nchavg    How many uv channels to average per image channel.
 *                  Ignored if uv data has multiple IFs.
 */
void 
ObitOTFUtilMakeCube (ObitImageDesc *inDesc, ObitOTFDesc *OTFDesc, 
		       ObitImageDesc *outDesc, 
		       gchar *Stokes, olong bchan, olong echan, olong incr, ObitErr *err)
{
  olong numberChann;
  gchar *name;
  gchar *routine = "ObitOTFUtilMakeCube";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitImageDescIsA(inDesc));
  g_assert (ObitOTFDescIsA(OTFDesc));

  /* Save output name */
  if (outDesc->name) name = g_strdup (outDesc->name);
  else  name = g_strdup ("Descriptor");

  /* Most info from inDesc */
  outDesc = ObitImageDescCopy (inDesc, outDesc, err);
  if (err->error) Obit_traceback_msg (err, routine, inDesc->name);

  /* restore name */
  if (outDesc->name) g_free(outDesc->name);
  outDesc->name = name;

  /* Set number of channels */
  numberChann = MIN (echan, OTFDesc->inaxes[OTFDesc->jlocf]) - MAX (1, bchan) + 1;
  outDesc->inaxes[outDesc->jlocf] = MAX (1, numberChann / MAX (1, incr));

  /* Stokes parameter */
  if ((Stokes[0]=='I') || (Stokes[0]==' ')) outDesc->crval[outDesc->jlocs] = 1.0;
  else if (Stokes[0]=='Q') outDesc->crval[outDesc->jlocs] =  2.0;
  else if (Stokes[0]=='U') outDesc->crval[outDesc->jlocs] =  3.0;
  else if (Stokes[0]=='V') outDesc->crval[outDesc->jlocs] =  4.0;
  else if (Stokes[0]=='R') outDesc->crval[outDesc->jlocs] = -1.0;
  else if (Stokes[0]=='L') outDesc->crval[outDesc->jlocs] = -1.0;

  /* reset image max/min */
  outDesc->maxval    = -1.0e20;
  outDesc->minval    =  1.0e20;

  return;
} /* end ObitOTFUtilMakeCube */

/**
 * Convolve a set of Clean components with a beam image
 * \param CCTab      CC Table, following parameters on infoList
 * \li  "BComp" OBIT_int (1,1,1) Start CC to use, 1-rel [def 1 ]
 * \li  "EComp" OBIT_int (1,1,1) Highest CC to use, 1-rel [def to end ]
 * \param Beam       Beam image to convolve with CCs
 * \param Template   Template for output array
 * \param err        Obit Error stack
 * \return An ObitFArray whose size is that of Template, spacing is that of Beam
 *         with the CCs in CCTab convolved with Beam and the (0,0) position is
 *         (nx/2,ny/2) (0-rel)
 */
ObitFArray* ObitOTFUtilConvBeam (ObitTableCC *CCTab, ObitImage *Beam, 
				 ObitFArray *Template, ObitErr *err)
{
  ObitFArray *out = NULL;
  ObitFArray *beamArray = NULL;
  olong iRow;
  ObitIOSize IOsize = OBIT_IO_byPlane;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong  *iptr, bcomp, ecomp;
  olong blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  olong  pos[2], loop, beamCen[2], nrow, nx, ny, xcen, ycen;
  ofloat xdelt, ydelt;
  ObitTableCCRow *row=NULL;
  gchar *routine = "ObitOTFUtilConvBeam";

 /* error checks */
  if (err->error) return out;
  g_assert (ObitTableCCIsA(CCTab));
  g_assert (ObitImageIsA(Beam));
  g_assert (ObitFArrayIsA(Template));

  /* Make output array */
  out = ObitFArrayCreate(routine, Template->ndim, Template->naxis);
  /* Size and center */
  nx = (olong)Template->naxis[0];
  xcen = nx/2;
  ny = (olong)Template->naxis[1];
  ycen = ny/2;

  /* Range of CCs - start CC number */
  if (ObitInfoListGetP(CCTab->info, "BComp",  &type, dim, (gpointer)&iptr)) {
    bcomp = iptr[0];
  } else bcomp = 1;

  /* End CC number */
  if (ObitInfoListGetP(CCTab->info, "EComp",  &type, dim, (gpointer)&iptr)) {
    ecomp = iptr[0];
  } else ecomp = 0;

  /* Read beam image - Full field */
  dim[0] = IM_MAXDIM;
  ObitInfoListPut (Beam->info, "BLC", OBIT_long, dim, blc, err); 
  ObitInfoListPut (Beam->info, "TRC", OBIT_long, dim, trc, err); 
  dim[0] = 1;
  ObitInfoListPut (Beam->info, "IOBy", OBIT_long, dim, &IOsize, err);
  Beam->extBuffer = FALSE;
  ObitImageOpen  (Beam, OBIT_IO_ReadOnly, err); 
  ObitImageRead  (Beam, NULL, err);
  ObitImageClose (Beam, err); 
  if (err->error) Obit_traceback_val (err, routine, Beam->name, out);
  beamArray = Beam->image;  /* FArray with beam */
  beamCen[0] = (olong)(Beam->myDesc->crpix[0] + 0.5);
  beamCen[1] = (olong)(Beam->myDesc->crpix[1] + 0.5);
  xdelt = 1.0 /  Beam->myDesc->cdelt[0];
  ydelt = 1.0 /  Beam->myDesc->cdelt[1];

  /* Open CC table */
  ObitTableCCOpen (CCTab, OBIT_IO_ReadOnly, err); 
  if (err->error) Obit_traceback_val (err, routine, Beam->name, out);
  row = newObitTableCCRow(CCTab);

  /* Restoration loop */
  nrow = CCTab->myDesc->nrow;
  if (ecomp<1) ecomp = nrow;
  Obit_log_error(err, OBIT_InfoErr, "%s: Using %d components", 
		 routine, ecomp-bcomp+1); 
  for (loop=bcomp; loop<=ecomp; loop++) {

    /* Read CC table */
    iRow = loop;
    ObitTableCCReadRow (CCTab, iRow, row, err);
    if (err->error) Obit_traceback_val (err, routine, CCTab->name, out);

    /* Restore to residual */
    pos[0] = (olong)(xcen + (row->DeltaX * xdelt) + 1.5);
    pos[1] = (olong)(ycen + (row->DeltaY * ydelt) + 1.5);
    ObitFArrayShiftAdd (out, pos, beamArray, beamCen, row->Flux, out);
  } /* End Restoration loop */


  /* Close CC table */
  ObitTableCCClose (CCTab, err);
  if (err->error) Obit_traceback_val (err, routine, CCTab->name, out);
  Beam->image = ObitFArrayUnref(Beam->image);

  /* Cleanup */
  row        = ObitTableCCRowUnref(row);
  return out;
} /* end ObitOTFUtilConvBeam */

/*----------------------Private functions---------------------------*/

/**
 * Fill in an image Descriptor from a OTF Descriptor.
 * Needs size, cell spacing and center filled in
 * Information about the first two axes other than the type an 
 * coordinate value need to be set separately.
 * to get the final position correct.
 * \param OTFDesc    Input OTF Descriptor.
 * \param imageDesc  Output image Descriptor
 * \param Proj       Projection code.
 */
static void 
ObitOTFUtilOTF2ImageDesc(ObitOTFDesc *OTFDesc, ObitImageDesc *imageDesc,
			 gchar *Proj)
{
  olong i, iaxis;
  gchar *st1;

  /* error checks */
  g_assert (ObitOTFDescIsA(OTFDesc));
  g_assert (ObitImageDescIsA(imageDesc));
  
  /* Be sure OTF descriptor is indexed */
  ObitOTFDescIndex(OTFDesc);

  /* loop over axes */

  iaxis = 0;
  /* RA axis, pos, inaxes, cdelt, crota, xshift set else where */
  /* Form label string */
  st1 = imageDesc->ctype[iaxis];
  st1[0] = 'R'; st1[1] = 'A'; st1[2] = '-'; st1[3] = '-'; 
  for (i=0; i<4; i++)  st1[i+4]=Proj[i];
  st1[9] = 0;

  /* Reference pixel */
  imageDesc->crpix[iaxis] = 1.0 + imageDesc->inaxes[iaxis] / 2.0;

  /* Dec axis, pos, inaxes, cdelt, crota, xshift set else where */
  iaxis++;
  /* Form label string */
  st1 = imageDesc->ctype[iaxis];
  st1[0] = 'D'; st1[1] = 'E'; st1[2] = 'C'; st1[3] = '-'; 
  for (i=0; i<4; i++)  st1[i+4]=Proj[i];
  st1[9] = 0;

  /* Reference pixel */
  imageDesc->crpix[iaxis] = 1.0 + imageDesc->inaxes[iaxis] / 2.0;

  /* Frequency Axis */
  iaxis++;
  /* Initially set for continuum */
  strncpy (imageDesc->ctype[iaxis], "FREQ    ", IMLEN_KEYWORD-1);
  imageDesc->inaxes[iaxis] = 1;  /* Only one for continuum */
  imageDesc->crpix[iaxis] = 1.0; /* reference pixel */
  imageDesc->crota[iaxis] = 0.0; /* no possible meaning */
  imageDesc->cdelt[iaxis] = OTFDesc->cdelt[OTFDesc->jlocf];
  if (OTFDesc->jlocf>=0)
    imageDesc->crval[iaxis] = OTFDesc->crval[OTFDesc->jlocf];
  else
    imageDesc->crval[iaxis] = 1.0e9; /* unknown frequency */

  /* Stokes Axis */
  iaxis++;
  imageDesc->inaxes[iaxis] = 1;  /* Only one */
  strncpy (imageDesc->ctype[iaxis], "STOKES  ", IMLEN_KEYWORD-1);
  imageDesc->crval[iaxis] = 1.0;
  imageDesc->crpix[iaxis] = 1.0; /* reference pixel */
  imageDesc->cdelt[iaxis] = 1.0; /* coordinate increment */
  imageDesc->crota[iaxis] = 0.0; /* no possible meaning */

  /* Total number of axes */
  imageDesc->naxis = iaxis+1;

  /* Copy information not directly related to an axis */
  /* Strings */
  strncpy (imageDesc->object, OTFDesc->object, IMLEN_VALUE-1);
  strncpy (imageDesc->teles,  OTFDesc->teles,  IMLEN_VALUE-1);
  strncpy (imageDesc->origin, OTFDesc->origin, IMLEN_VALUE-1);
  strncpy (imageDesc->bunit,  "JY/BEAM ",     IMLEN_VALUE-1);
  /* Set current date */
  ObitOTFUtilCurDate (imageDesc->date, IMLEN_VALUE-1);

  /* Observing date */
  if (OTFDesc->JDObs>1.0) ObitOTFDescJD2Date (OTFDesc->JDObs, imageDesc->obsdat);

  imageDesc->epoch        = OTFDesc->epoch;
  imageDesc->obsra        = OTFDesc->obsra;
  imageDesc->obsdec       = OTFDesc->obsdec;

  /* initialize some values */
  imageDesc->areBlanks = FALSE;
  imageDesc->niter     = 0;
  imageDesc->maxval    = -1.0e20;
  imageDesc->minval    =  1.0e20;
  imageDesc->bitpix    = -32;
  imageDesc->beamMaj   = OTFDesc->beamSize;
  imageDesc->beamMin   = OTFDesc->beamSize;
  imageDesc->beamPA    = 0.0;

  /* Index Image descriptor */
  ObitImageDescIndex(imageDesc);

} /* end ObitOTFUtilOTF2ImageDesc */

/**
 * Fills an existing character array with the string for the current date.
 * \param date Character string to accept the string (10 char+null)
 * \param len  Actual length of date (should be at least 11)
 */
static void ObitOTFUtilCurDate (gchar *date, olong len)
{
  struct tm *lp;
  time_t clock;
  
  /* Get time since 00:00:00 GMT, Jan. 1, 1970 in seconds. */
  time (&clock);
  
  /* Convert to  broken-down time. */
  lp = localtime (&clock);
  
  /* Full year */
  if (lp->tm_year<1000)  lp->tm_year += 1900; 
  lp->tm_mon++; /* Month 0-rel rest 1-rel */

  /* to output */
  g_snprintf (date, len, "%4.4d-%2.2d-%2.2d",
	      lp->tm_year, lp->tm_mon, lp->tm_mday);
} /* end ObitOTFUtilCurDate */

/**
 * Subtract the values in an image from a portion of a buffer of OTF data.
 * For CCB beamswitched data (OTFType=OBIT_GBTOTF_CCB, no. States=1)
 * Callable as thread
 * Arguments are given in the structure passed as arg
 * \param arg Pointer to SubImageFuncArg argument with elements:
 * \li sky     OTFSkyModel 
 * \li otfdata OTF data set to model and subtract from current buffer
 * \li first  First (1-rel) rec in otfdata buffer to process this thread
 * \li last   Highest (1-rel) rec inotfdata buffer to process this thread
 * \li ithread thread number
 * \li err Obit error stack object.
 * \li factor  Scaling factor for sky model
 * \li Interp  Image Interpolator
 * \li xpos, ypos, float arrays the size of ndetect
 * \return NULL
 */
gpointer ThreadOTFUtilSubImageBuff (gpointer args)
{
  SubImageFuncArg *largs  = (SubImageFuncArg*)args;
  ObitOTF *in             = largs->otfdata;
  olong loRec             = largs->first-1;
  olong hiRec             = largs->last;
  ofloat factor           = largs->factor;
  ObitFInterpolate *image = largs->Interp;
  ofloat *xpos            = largs->xpos;
  ofloat *ypos            = largs->ypos;
  ObitErr *err            = largs->err;

  ofloat *data, ffact, value, RACenter, DecCenter;
  ObitOTFProj proj;
  odouble coord[IM_MAXDIM];
  olong  ndetect, ndata, i, j;
  olong incfeed, itemp,incdatawt ;
  gboolean CCBBS, isRef;
  ObitOTFDesc* desc;
  ObitOTFArrayGeom* geom;
  ofloat fblank = ObitMagicF();
  gchar Proj[5];

  /* Error checks */
  if (err->error) goto finish;

  /* Local pointers */
  desc = in->myDesc;
  geom = in->geom;
  data = in->buffer+loRec*desc->lrec;  /* Appropriate offset in buffer */
  ndetect = geom->numberDetect;        /* How many detectors */
  ndata   = desc->numRecBuff;          /* How many data records */
  if (in->myDesc->jlocfeed>=0) 
    incfeed = in->myDesc->incfeed / in->myDesc->incdatawt;
  else incfeed = 1;  /* This is probably a bad sign */

  /* Get Model center */
  RACenter  = image->myDesc->crval[0];
  DecCenter = image->myDesc->crval[1];
  strncpy (Proj, &image->myDesc->ctype[0][4], 5);
  proj      = ObitOTFSkyModelProj (Proj);
  incdatawt = desc->incdatawt; /* increment in data-wt axis */

  /* Is this CCB beamswitched data? */
  CCBBS = (in->myDesc->OTFType==OBIT_GBTOTF_CCB) &&
    (in->myDesc->jlocstate>=0) &&
    (in->myDesc->inaxes[in->myDesc->jlocstate]==1);

  ffact = factor;

  /* Loop over data records */
  for (i=loRec; i<hiRec; i++) {

    /* Get Sky locations of the data projected onto the image */
    ObitOTFArrayGeomProj(geom, data[desc->ilocra], data[desc->ilocdec], 
			 data[desc->ilocrot], RACenter, DecCenter, proj, 
			 xpos, ypos);

    /* Loop over the array */
    for (j=0; j<ndetect; j++) {
      
     /* Interpolate - use coordinates on a flat plane at the tangent point of the image 
	 xpos and ypos are offsets from the image center projected onto the plane 
	 of the image, ObitFInterpolateOffset does a linear approximation. */
       coord[0] = xpos[j];  /* + RACenter; */
       coord[1] = ypos[j];  /* + DecCenter; */
       value = ObitFInterpolateOffset (image, coord, err); 

       /* debug 
       if ((fabs(data[desc->ilocdata+j]-value*ffact)>1.0)  && (j<=3)) {
	 fprintf (stderr,"debugSig: val %g pos %lf %lf in data %f %f pos %f %f\n", 
		 value, coord[0]*3600.0, coord[1]*3600.0, 
		 data[desc->ilocdata+j*incdatawt], (data[desc->ilocdata+j*incdatawt]-value*ffact),
		 data[desc->ilocra], data[desc->ilocdec]);
       } */
       /* debug  
       if ((fabs(coord[0])<0.0014) && (fabs(coord[1])<0.0014) && (j<=3)) {
	 fprintf (stderr,"debugSig: val %g pos %lf %lf in data %f %f pos %f %f\n", 
		 value, coord[0]*3600.0, coord[1]*3600.0, 
		 data[desc->ilocdata+j*incdatawt], (data[desc->ilocdata+j*incdatawt]-value*ffact),
		 data[desc->ilocra], data[desc->ilocdec]);
      }*/
      /* debug  
      if (value>0.2) {
	fprintf (stderr,"debug: val %g pos %lf %lf in data %f %f pos %f %f\n", 
		 value, coord[0]*3600.0, coord[1]*3600.0, 
		 data[desc->ilocdata], factor*(data[desc->ilocdata]-value),
		 data[desc->ilocra], data[desc->ilocdec]);
      }*/

      /* Subtract */
      if ((value!=fblank) && (data[desc->ilocdata+j*incdatawt]!=fblank))
	data[desc->ilocdata+j*incdatawt] -= value * ffact;
      else {
	data[desc->ilocdata+j*incdatawt] = fblank;
      }

      /* For CCB beamswitched data add value corresponding to this feeds
       reference beam */
      if (CCBBS) {
	/* is this the reference or signal beam? */
	itemp = j / incfeed;
	/* itemp odd means reference beam */
	isRef = itemp != 2*(itemp/2);
	/* Use feed offset for feed 1 if this is reference, else feed 2 */
	if (isRef) itemp = 0;
	else itemp = incfeed;
	/* interpolate */
	coord[0] = xpos[itemp]; /* + RACenter; */
	coord[1] = ypos[itemp]; /* + DecCenter; */
	value = ObitFInterpolateOffset (image, coord, err);
       /* debug 
       if ((fabs(data[desc->ilocdata+j*incdatawt]+value*ffact)>1.0)  && (j<=3)) {
	 fprintf (stderr,"debugRef: val %g pos %lf %lf in data %f %f pos %f %f\n", 
		 value, coord[0]*3600.0, coord[1]*3600.0, 
		 data[desc->ilocdata+j*incdatawt], (data[desc->ilocdata+j*incdatawt]+value*ffact),
		 data[desc->ilocra], data[desc->ilocdec]);
       } */
       /* debug 
       if ((fabs(coord[0])<0.0014) && (fabs(coord[1])<0.0014) && (j<=3)) {
	 fprintf (stderr,"debugRef: val %g pos %lf %lf in data %f %f pos %f %f\n", 
		 value, coord[0]*3600.0, coord[1]*3600.0, 
		 data[desc->ilocdata+j*incdatawt], (data[desc->ilocdata+j*incdatawt]+value*ffact),
		 data[desc->ilocra], data[desc->ilocdec]);
       } */
       /* Add */
       if ((value!=fblank) && (data[desc->ilocdata+j*incdatawt]!=fblank))
	 data[desc->ilocdata+j*incdatawt] += value * ffact;
       else {
	 data[desc->ilocdata+j*incdatawt] = fblank;
       }
      }  /* end adding to reference beam position */
      
    } /* end loop over array */
    data += desc->lrec; /* update buffer pointer */
  } /* end loop over buffer */
  
  /* Indicate completion */
 finish: ObitThreadPoolDone (in->thread, (gpointer)&largs->ithread);
  return NULL;
} /* end ThreadOTFUtilSubImageBuff */

/**
 * Make arguments for Threaded OTFUtilSubImageBuff
 * \param in         OTF with internal buffer to be modified.
 * \param sky        OTFSkyModel to subtract.
 * \param factor     Scaling factor for sky model.
 * \param err        Obit error stack object.
 * \param args       Created array of SubImageFuncArg, 
 *                   delete with KillOTFUtilSubImageArgs
 * \return number of elements in args.
 */
static glong MakeOTFUtilSubImageArgs (ObitOTF *in, ObitErr *err, 
				      SubImageFuncArg ***args)
{
  olong i, nThreads, ndetect;

  /* Setup for threading */
  /* How many threads? */
  nThreads = MAX (1, ObitThreadNumProc(in->thread));

  /* Initialize threadArg array */
  *args = g_malloc0(nThreads*sizeof(SubImageFuncArg*));
  for (i=0; i<nThreads; i++) 
    (*args)[i] = g_malloc0(sizeof(SubImageFuncArg)); 
  
  ndetect = in->geom->numberDetect;    /* How many detectors */

  for (i=0; i<nThreads; i++) {
    (*args)[i]->otfdata = in;
    (*args)[i]->ithread = i;
    (*args)[i]->err     = err;
    (*args)[i]->Interp  = NULL;
    (*args)[i]->factor  = 1.0;
    (*args)[i]->xpos    = g_malloc0(ndetect*sizeof(ofloat));
    (*args)[i]->ypos    = g_malloc0(ndetect*sizeof(ofloat));
  }

  return nThreads;
} /*  end MakeOTFUtilSubImageArgs */

/**
 * Delete arguments for Threaded OTFUtilSubImageBuff
 * \param nargs      number of elements in args.
 * \param args       Array of SubImageFuncArg, type SubImageFuncArg
 */
static void KillOTFUtilSubImageArgs (olong nargs, SubImageFuncArg **args)
{
  olong i;

  if (args==NULL) return;
  for (i=0; i<nargs; i++) {
    if (args[i]) {
      if (args[i]->Interp) ObitFInterpolateUnref(args[i]->Interp);
      if (args[i]->xpos) g_free (args[i]->xpos);
      if (args[i]->ypos) g_free (args[i]->ypos);
      g_free(args[i]);
    }
  }
  g_free(args);
} /*  end KillOTFUtilSubImageArgs */
