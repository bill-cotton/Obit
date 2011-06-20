/* $Id$        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2011                                               */
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

#include "ObitPixHisto.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitPixHisto.c
 * ObitPixHisto class function definitions.
 * This class is derived from the Obit base class.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitPixHisto";

/** Function to obtain parent ClassInfo */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo structure ObitPixHistoClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitPixHistoClassInfo myClassInfo = {FALSE};

/*--------------- File Global Variables  ----------------*/


/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitPixHistoInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitPixHistoClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitPixHistoClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitPixHisto* newObitPixHisto (gchar* name)
{
  ObitPixHisto* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitPixHistoClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitPixHisto));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitPixHistoInit((gpointer)out);

 return out;
} /* end newObitPixHisto */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitPixHistoGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitPixHistoClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitPixHistoGetClass */

/**
 * Make a deep copy of an ObitPixHisto.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitPixHisto* ObitPixHistoCopy  (ObitPixHisto *in, ObitPixHisto *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  gboolean oldExist;
  olong i;
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
    out = newObitPixHisto(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  /* Delete any old */
  out->image    = ObitImageUnref(out->image);
  out->imagePix = ObitFArrayUnref(out->imagePix);
  out->histo    = ObitFArrayUnref(out->histo);
  out->intHisto = ObitFArrayUnref(out->intHisto);

  /*  copy this class */
  out->image    = ObitImageCopy (in->image,    out->image,    err);
  out->imagePix = ObitFArrayCopy(in->imagePix, out->imagePix, err);
  out->histo    = ObitFArrayCopy(in->histo,    out->histo,    err);
  out->intHisto = ObitFArrayCopy(in->intHisto, out->intHisto, err);
  for (i=0; i<MAXFARRAYDIM; i++) {
    out->blc[i] = in->blc[i];
    out->trc[i] = in->trc[i];
  }

  return out;
} /* end ObitPixHistoCopy */

/**
 * Make a copy of a object but do not copy the actual data
 * This is useful to create an PixHisto similar to the input one.
 * \param in  The object to copy
 * \param out An existing object pointer for output, must be defined.
 * \param err Obit error stack object.
 */
void ObitPixHistoClone  (ObitPixHisto *in, ObitPixHisto *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  olong i;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitIsA(out, &myClassInfo));

  /* shallow copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitClone ((Obit*)in, (Obit*)out);

  /* Delete any old */
  out->image    = ObitImageUnref (out->image);
  out->imagePix = ObitFArrayUnref(out->imagePix);
  out->histo    = ObitFArrayUnref(out->histo);
  out->intHisto = ObitFArrayUnref(out->intHisto);

  /*  copy this class */
  out->image    = ObitImageRef (in->image);
  out->imagePix = ObitFArrayRef(in->imagePix);
  out->histo    = ObitFArrayRef(in->histo);
  out->intHisto = ObitFArrayRef(in->intHisto);
  for (i=0; i<MAXFARRAYDIM; i++) {
    out->blc[i] = in->blc[i];
    out->trc[i] = in->trc[i];
  }

} /* end ObitPixHistoClone */

/**
 * Creates an ObitPixHisto 
 * \param name  An optional name for the object.
 * \param image Image to use, image plane selectioin is by any values 
 *              of BLC and TRC on info
 * \param err   Obit error/message stack.
 * \return the new object.
 */
ObitPixHisto* ObitPixHistoCreate (gchar* name, ObitImage *image, 
				  ObitErr *err)
{
  ObitPixHisto* out;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong tblc[IM_MAXDIM], blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong ttrc[IM_MAXDIM], trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  olong i;
  gchar *routine=" ObitPixHistoCreate";

  /* error checks */
  if (err->error) return NULL;
  g_assert (ObitImageIsA(image));

  /* Create basic structure */
  out = newObitPixHisto (name);

  /* Save basic image */
  out->image = ObitImageRef(image);

  /* Get window - use for plane */
  ObitInfoListGetTest(image->info, "BLC", &type, dim, blc);
  ObitInfoListGetTest(image->info, "TRC", &type, dim, trc);

  /* Set selection */
  tblc[0] = 1; tblc[1] = 1; for (i=2; i<IM_MAXDIM; i++) tblc[i] = blc[i];
  ttrc[0] = 0; ttrc[1] = 0; for (i=2; i<IM_MAXDIM; i++) ttrc[i] = trc[i];
  ObitInfoListAlwaysPut(image->info, "BLC", OBIT_long, dim, tblc);
  ObitInfoListAlwaysPut(image->info, "TRC", OBIT_long, dim, ttrc);

  /* Open image */
  if ((ObitImageOpen (image, OBIT_IO_ReadOnly, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s: ERROR opening image %s", 
		   routine, image->name);
    return out;
  }

  /* Read input plane */
  if ((ObitImageRead (image, NULL , err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "%s: ERROR reading image %s", 
		   routine, image->name);
    return out;
  }

  /* Copy Pixel array */
  out->imagePix =  ObitFArrayCopy(image->image, out->imagePix, err);
  if (err->error) Obit_traceback_val (err, routine, image->name, out);
  
  /* close input */  
  ObitImageClose (image, err) ;
  if (err->error) Obit_traceback_val (err, routine, image->name, out);

  /* Reset BLC, TRC */
  ObitInfoListAlwaysPut(image->info, "BLC", OBIT_long, dim, blc);
  ObitInfoListAlwaysPut(image->info, "TRC", OBIT_long, dim, trc);

  /* Free image buffer if not memory resident */
  if (image->mySel->FileType!=OBIT_IO_MEM) 
    image->image = ObitFArrayUnref(image->image);

  return out;
} /* end ObitPixHistoCreate */

/**
 * Calculates a histogram for a region in the PixHisto image pixel array
 * Does both differential and integral histograms.
 * \param in    PixHisto object
 * \param blc   Bottom left corner (1-rel) in imagePix
 *              only 2 elements used.
 * \param trc   Top right corner (1-rel) in imagePix
 * \param nbin  Number of bins in the histogram
 * \param range Range of values around zero in units of the RMS
 * \param err   Obit error/message stack.
 */
void ObitPixHistoHisto ( ObitPixHisto *in, 
			 olong *blc, olong *trc,
			 olong nbin, ofloat range,
			 ObitErr *err)
{
  ofloat *hi, *hd, under=0.0, over=0.0, total=0.0;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  olong i, j, lblc[2], ltrc[2];
  odouble sum;
  ObitFArray *tempFA=NULL;
  gchar *routine=" ObitPixHistoHisto";

  /* error checks */
  if (err->error) return;
  Obit_return_if_fail(ObitFArrayIsA(in->imagePix), err,
		      "%s:No image pixel array for %s ", 
		      routine, in->name);


  /* Copy pixels to temporary array */
  lblc[0] = blc[0]-1; lblc[1] = blc[1]-1;
  ltrc[0] = trc[0]-1; ltrc[1] = trc[1]-1;
  tempFA = ObitFArraySubArr (in->imagePix, lblc, ltrc, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
  in->blc[0] = blc[0]; in->blc[1] = blc[1];
  in->trc[0] = trc[0]; in->trc[1] = trc[1]; 

  /* Get RMS */
  in->sigma = ObitFArrayRMS(tempFA);

  /* Center, width info */
  in->cenx = (blc[0]+trc[0])/2;
  in->ceny = (blc[1]+trc[1])/2;
  in->FDRsize = trc[0] - in->cenx;  /* some assumptions here */

  /* histogram */
  in->histo = ObitFArrayUnref(in->histo);
  in->histo = ObitFArrayHisto(tempFA, nbin, 
			      -range*in->sigma, range*in->sigma);
  tempFA = ObitFArrayUnref(tempFA);

  /* Integrated histogram */
  in->intHisto = ObitFArrayUnref(in->intHisto);
  in->intHisto =  ObitFArrayCopy(in->histo, in->intHisto, err); /* Copy */
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Underflows, overflows, total */
  ObitInfoListGetTest(in->histo->info, "Under",  &type, dim, &under);
  ObitInfoListGetTest(in->histo->info, "Over",   &type, dim, &over);
  ObitInfoListGetTest(in->histo->info, "Total",  &type, dim, &total);

  /* Integrate - positive */
  hi = in->intHisto->array;
  hd = in->histo->array;
  for (j=nbin/2+1; j<nbin; j++) {
    sum = (odouble)over;
    for (i=j; i<nbin; i++) sum += hd[i];
    hi[j] = (ofloat)sum;
  }

  /* negative */
  for (j=nbin/2-1; j>=0; j--) {
    sum = (odouble)under;
    for (i=0; i<=j; i++) sum += hd[i];
    hi[j] = (ofloat)sum;
  }
  /* Center = total */
  hi[nbin/2] = total;

} /* end ObitPixHistoHisto */

/**
 *  Calculate minimum flux density for a given FDR
 * If all else fails, the returned value is 5 * sigma.
 * \param in      PixHisto object
 * \param maxFDR  Fractional acceptable False Detection Rate, 0=>0.5
 * \param err     Obit error/message stack.
 * \return Flux density level with expected false detection rate of maxFDR
 */
ofloat ObitPixHistoFDRFlux (ObitPixHisto *in, ofloat maxFDR,
			    ObitErr *err)
{
  ofloat out=in->sigma*5;
  olong i, j, nbin=0, cen;
  ofloat amin=0.0, amax=0.0, cell, fdr, flux, *hi, test=0.0;
  ofloat x1, y1, x2, y2;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *routine=" ObitPixHistoFDRFlux";

  /* error checks */
  if (err->error) return out;
  Obit_retval_if_fail(ObitFArrayIsA(in->histo), err, out,
		      "%s:No pixel histogram for %s ", 
		      routine, in->name);

  /* Get histogram info */
  ObitInfoListGetTest(in->histo->info, "nHisto",  &type, dim, &nbin);
  ObitInfoListGetTest(in->histo->info, "Min",     &type, dim, &amin);
  ObitInfoListGetTest(in->histo->info, "Max",     &type, dim, &amax);
  cell = (amax - amin)/(nbin-1.);
  cen = nbin/2+1;

  if (err->prtLv>=2) 
    Obit_log_error(err, OBIT_InfoErr, 
		   "Histogram analysis: Cen: %d %d, hw: %d, cell: %f",
		   in->cenx, in->ceny, in->FDRsize, cell);
  hi = in->intHisto->array;
  j = nbin/2;
  x1 = y1 = 0.0;
  for (i=nbin/2+1; i<nbin; i++) {
    flux = (i-cen+0.5) * cell;
    j--;
    fdr = 1.0 - (hi[i]-hi[j])/hi[i];
    x2 = flux;
    y2 = fdr;
    /* Diagnostics? */
    if (err->prtLv>=2) {
      Obit_log_error(err, 
		     OBIT_InfoErr, "%2.2d %2.2d hist: %8.0f flux: %6.4f FDR: %6.4f", 
		     i, j, hi[i], flux, fdr);
    }
    /* Cross the desired value? */
    if (fdr<maxFDR) {
      test = x1 + cell * ((y1-maxFDR) / MAX(0.0001,(y1-y2)));
      /* sanity check */
      if ((test>in->sigma) && (test<in->sigma*10)) out  = test; 
      break;
    }
    x1 = x2; y1 = y2;
  }

  /* diagnostics? */
  if (err->prtLv>=1) {
   Obit_log_error(err, 
		  OBIT_InfoErr, "FDR lim: %8.5f, range: %8.5f, sigma: %8.5f", 
		  out, amax, in->sigma);
  }

  return out;
} /* end ObitPixHistoFDRFlux */


/**
 * Initialize global ClassInfo Structure.
 */
void ObitPixHistoClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitPixHistoClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitPixHistoClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitPixHistoClassInfoDefFn (gpointer inClass)
{
  ObitPixHistoClassInfo *theClass = (ObitPixHistoClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitPixHistoClassInit;
  theClass->newObit       = (newObitFP)newObitPixHisto;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitPixHistoClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitPixHistoGetClass;
  theClass->ObitCopy      = (ObitCopyFP)ObitPixHistoCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitPixHistoClear;
  theClass->ObitInit      = (ObitInitFP)ObitPixHistoInit;
  theClass->ObitPixHistoCreate  = (ObitPixHistoCreateFP)ObitPixHistoCreate;
  theClass->ObitPixHistoHisto   = (ObitPixHistoHistoFP)ObitPixHistoHisto;
  theClass->ObitPixHistoFDRFlux = (ObitPixHistoFDRFluxFP)ObitPixHistoFDRFlux;
} /* end ObitPixHistoClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitPixHistoInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitPixHisto *in = inn;
  olong i;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->info      = newObitInfoList(); 
  in->image     = NULL;
  in->imagePix  = NULL;
  in->histo     = NULL;
  in->intHisto  = NULL;
  for (i=0; i<MAXFARRAYDIM; i++) {
    in->blc[i] = 1;
    in->trc[i] = 0;
  }
} /* end ObitPixHistoInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitPixHisto* cast to an Obit*.
 */
void ObitPixHistoClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitPixHisto *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  in->info      = ObitInfoListUnref(in->info);
  in->image     = ObitImageUnref(in->image);
  in->imagePix  = ObitFArrayUnref(in->imagePix);
  in->histo     = ObitFArrayUnref(in->histo);
  in->intHisto  = ObitFArrayUnref(in->intHisto);
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitPixHistoClear */

