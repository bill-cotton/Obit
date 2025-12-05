/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2025                                               */
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

#include "ObitSpectrumInterp.h"
#include "ObitPosition.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitSpectrumInterp.c
 * ObitSpectrumInterp class function definitions.
 * This class is derived from the Obit base class.
 * This class supports interpolation in a spectrum using 
 * Lagrange interpolation.
 */

/*--------------- File Global Variables  ----------------*/
/** name of the class defined in this file */
static gchar *myClassName = "ObitSpectrumInterp";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo structure ObitSpectrumInterpClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitSpectrumInterpClassInfo myClassInfo = {FALSE};


/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitSpectrumInterpInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitSpectrumInterpClear (gpointer in);

/** Private: Set convolution kernal */
static void SetConvKernal (odouble Freq, olong nFreq, odouble *Freqs, ofloat *chWt,
			   olong hwidth, ofloat *denom, olong *Start,
			   odouble *cenFreq, ofloat *Kernal);

/** Private: Set denominator */
static void SetDenom (olong nFreq, olong hwidth, olong Start,
		      odouble *Freqs, ofloat *denom);

/** Private: Set Class function pointers. */
static void ObitSpectrumInterpClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitSpectrumInterp* newObitSpectrumInterp (gchar* name)
{
  ObitSpectrumInterp* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitSpectrumInterpClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitSpectrumInterp));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitSpectrumInterpInit((gpointer)out);

 return out;
} /* end newObitSpectrumInterp */

/**
 * Constructor.
 * NOTE: Use GHz as frequencies, else blow floating point range
 * Initializes class if needed on first call.
 * \param name     An optional name for the object.
 * \param nFreq    Number of frequencies in spectrum
 * \param Freqs    Array of frequencies in spectrum
 * \param chWt     Channel weights, NULL => all 1.0
 * \param Spectrum Pointer to Spectrum 
 * \param hwidth   Half width of interpolation kernal (range [1,3] allowed).
                   out of range =>3
 * \return the new object.
 */
ObitSpectrumInterp* newObitSpectrumInterpCreate (gchar* name, olong nFreq,
						 odouble *Freqs, ofloat *Spectrum,
						 ofloat *chWt, olong hwidth)
{
  ObitSpectrumInterp* out;
  olong i;

  /* Create/init output structure */
  out = newObitSpectrumInterp (name);

  out->nFreq = nFreq;  /* How many channels */
  if ((hwidth>0)&&(hwidth<4)) out->hwidth = hwidth;
  else                        out->hwidth = 3;
  /* Attach pointer */
  out->Spectrum = Spectrum;

  /* Create work Arrays - bigger than necessary */
  out->Freqs  = g_malloc0(nFreq*sizeof(odouble));
  out->chWt   = g_malloc0(nFreq*sizeof(ofloat));
  out->Kernal = g_malloc0(nFreq*sizeof(ofloat));
  out->denom  = g_malloc0(nFreq*sizeof(ofloat));

  /* Copy Data */
  for (i=0; i<nFreq; i++) {
    out->Freqs[i] = Freqs[i];
    if (chWt) out->chWt[i] = chWt[i];
    else      out->chWt[i] = 1.0;
  }

  /* Init Lagrangian denominators for hwidth, Start */
  out->Start = 0;
  SetDenom (nFreq, hwidth, out->Start, Freqs, out->denom);
  return out; 
} /* end newObitSpectrumInterpCreate */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitSpectrumInterpGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitSpectrumInterpClassInit();

  return (gconstpointer)&myClassInfo;
} /* end  ObitSpectrumInterpGetClass */

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
ObitSpectrumInterp* ObitSpectrumInterpCopy (ObitSpectrumInterp *in, ObitSpectrumInterp *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  gboolean oldExist;
  olong i;
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
    out = newObitSpectrumInterp(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /* copy/set this classes additions */
  out->nFreq      = in->nFreq;
  out->Start      = in->Start;
  out->Freq       = in->Freq;
  out->cenFreq    = in->cenFreq;
  out->value      = in->value;
  out->hwidth     = in->hwidth;
  /* copy pointer */
  out->Spectrum = in->Spectrum;
  /* Copy arrays */
  if (out->Kernal) g_free(out->Kernal);
  out->Kernal  = g_malloc0(in->nFreq*sizeof(ofloat));
  for (i=0; i<in->nFreq; i++) out->Kernal[i] = in->Kernal[i];
  if (out->denom) g_free(out->denom);
  out->Freqs  = g_malloc0(in->nFreq*sizeof(ofloat));
  for (i=0; i<in->nFreq; i++) out->Freqs[i] = in->denom[i];
  if (out->Freqs) g_free(out->Freqs);
  out->Freqs  = g_malloc0(in->nFreq*sizeof(odouble));
  for (i=0; i<in->nFreq; i++) out->Freqs[i] = in->Freqs[i];
  if (out->chWt) g_free(out->chWt);
  out->chWt  = g_malloc0(in->nFreq*sizeof(ofloat));
  for (i=0; i<in->nFreq; i++) out->chWt[i] = in->chWt[i];
  
  return out;
} /* end ObitSpectrumInterpCopy */

/**
 * Make a shallow copy of a object.
 * The result will have pointers to the more complex members.
 * Parent class members are included but any derived class info is ignored.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \return pointer to the new object.
 */
ObitSpectrumInterp* ObitSpectrumInterpClone  (ObitSpectrumInterp *in, ObitSpectrumInterp *out)
{
  const ObitClassInfo *myClass, *ParentClass;
  gboolean oldExist;
  olong i;
  gchar *outName;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  if (out) g_assert (ObitIsA(out, &myClassInfo));

  /* Create if it doesn't exist */
  oldExist = out!=NULL;
  if (!oldExist) {
    /* derive object name */
    outName = g_strconcat ("Clone: ",in->name,NULL);
    out = newObitSpectrumInterp(outName);
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
  out->nFreq      = in->nFreq;
  out->Start      = in->Start;
  out->Freq       = in->Freq;
  out->cenFreq    = in->cenFreq;
  out->value      = in->value;
  out->hwidth   = in->hwidth;
  /* copy pointers*/
  out->Spectrum = in->Spectrum;
  /* Copy arrays */
  if (out->Kernal) g_free(out->Kernal);
  out->Kernal  = g_malloc0(in->nFreq*sizeof(ofloat));
  for (i=0; i<in->nFreq; i++) out->Kernal[i] = in->Kernal[i];
  if (out->denom) g_free(out->denom);
  out->denom  = g_malloc0(in->nFreq*sizeof(ofloat));
  for (i=0; i<in->nFreq; i++) out->denom[i] = in->denom[i];
   if (out->Freqs) g_free(out->Freqs);
  out->Freqs  = g_malloc0(in->nFreq*sizeof(odouble));
  for (i=0; i<in->nFreq; i++) out->Freqs[i] = in->Freqs[i];
  if (out->chWt) g_free(out->chWt);
  out->chWt  = g_malloc0(in->nFreq*sizeof(ofloat));
  for (i=0; i<in->nFreq; i++) out->chWt[i] = in->chWt[i];
  
  return out;
} /* end ObitSpectrumInterpClone */

/**
 * Replace the Spectrum member to be interpolated.
 * \param in        The object to update
 * \param Spectrum  Pointer to spectrum
 */
void ObitSpectrumInterpReplace  (ObitSpectrumInterp *in, ofloat *Spectrum)
{
  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  in->Spectrum = Spectrum;
} /* end ObitSpectrumInterpReplace */

/**
 * Set channel weights based on a Spectrum
 * Weight set to 1 whre spectrum not 0 or fblank
 * \param in        The object to update
 * \param Spectrum  Pointer to spectrum
 */
void ObitSpectrumInterpSetWt  (ObitSpectrumInterp *in, ofloat *Spectrum)
{
  ofloat fblank =  ObitMagicF();
  olong i;
  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  for (i=0; i<in->nFreq; i++) {
    if ((Spectrum[i]!=0.0) && (Spectrum[i]!=fblank)) in->chWt[i] = 1.0;
    else                                             in->chWt[i] = 0.0;
  }
  in->Spectrum = Spectrum;
} /* end ObitSpectrumInterpReplace */

/**
 * Interpolate value at requested cell in a spectrum
 * Values out of the frequency range are extrapolations of the first/last two
 * \param in    The object to interpolate
 * \param Freq  The frequency at which the spectral value is desired.
 * \return value, magic blanked if invalid
 */
ofloat ObitSpectrumInterpFreq (ObitSpectrumInterp *in, odouble Freq)
{
  ofloat fblank =  ObitMagicF();
  ofloat value, sum, sumwt, sum2, wt;
  ofloat *Kernal, *data, *chWt;
  olong i, iwid;
  olong indx, Start;

  /* Special handling outside of range of Freqs */
  if (Freq<in->Freqs[0]) {                   /* Before frequency range */
    if  ((in->chWt[0]==0.0) || (in->Spectrum[1]==0.0))  /* No initial spectral value */
      {in->value = fblank; return in->value;} 
    value = in->Spectrum[0];
    if ((in->chWt[1]>0.0) && (in->Spectrum[1]!=0.0))
      value -= (in->Spectrum[0]-in->Spectrum[1]) *
	((Freq-in->Freqs[0])/(in->Freqs[1]-in->Freqs[0]));
    in->value = value; return in->value;
  } else if (Freq>in->Freqs[in->nFreq-1]) {  /* After frequency range */
    if  ((in->chWt[in->nFreq-1]==0.0) || (in->Spectrum[in->nFreq-1]==0.0))  /* No final spectral value */
      {in->value = fblank; return in->value;} 
    value = in->Spectrum[in->nFreq-1];
    if ((in->chWt[in->nFreq-1]>0.0) && (in->Spectrum[in->nFreq-1]!=0.0))
      value -= (in->Spectrum[in->nFreq-1]-in->Spectrum[in->nFreq-2]) *
	((Freq-in->Freqs[in->nFreq-1])/(in->Freqs[in->nFreq-2]-in->Freqs[in->nFreq-1]));
   in->value = value; return in->value;
  } else {                                   /* In frequency range */

    /* Update convolving kernal as needed */
    if (Freq != in->Freq) {
      in->Freq = Freq;
      SetConvKernal (Freq, in->nFreq, in->Freqs, in->chWt, in->hwidth,
		     in->denom, &in->Start, &in->cenFreq, in->Kernal);
    }
    
    /* Local versions of things */
    Kernal  = in->Kernal;
    chWt    = in->chWt;
    data    = in->Spectrum;
    Start   = in->Start;
    
    /* Zero sums */
    sum   = 0.0;
    sum2  = 0.0;
    sumwt = 0.0;
    
    /* Loop over data summing values times convolving weights */
    iwid = 2*in->hwidth+1;
    for (i=0; i<iwid; i++) {
      indx = Start + i;
      if (data[indx] != fblank) {
	wt = Kernal[i]*chWt[indx];
	sumwt += wt;
	sum   += data[indx] * wt;
	sum2  += Kernal[i];
      } 
    }
    
    /* normalize sum if not excessive blanking */
    if (sum2 > 0.20) {
      in->value = sum / sumwt;
      return in->value;
    } 
    
    /* No luck - return fblank */
    return fblank;
  } /* end in freq range */
}
/* end  ObitSpectrumInterpFreq */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitSpectrumInterpClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitSpectrumInterpClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitSpectrumInterpClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitSpectrumInterpClassInfoDefFn (gpointer inClass)
{
  ObitSpectrumInterpClassInfo *theClass = (ObitSpectrumInterpClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitSpectrumInterpClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitSpectrumInterpClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitSpectrumInterpGetClass;
  theClass->newObit       = (newObitFP)newObitSpectrumInterp;
  theClass->ObitCopy      = (ObitCopyFP)ObitSpectrumInterpCopy;
  theClass->ObitClone     = (ObitCloneFP)ObitSpectrumInterpClone;
  theClass->ObitClear     = (ObitClearFP)ObitSpectrumInterpClear;
  theClass->ObitInit      = (ObitInitFP)ObitSpectrumInterpInit;

} /* end ObitSpectrumInterpClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param in Pointer to the object to initialize.
 */
void ObitSpectrumInterpInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitSpectrumInterp *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->info        = newObitInfoList();
  in->nFreq       = -1;
  in->hwidth      = 0;
  in->Start       = -100;
  in->Freq        = 0.0;
  in->value       = 0.0;
  in->Freqs       = NULL;
  in->chWt        = NULL;
  in->Spectrum    = NULL;
  in->Kernal      = NULL;
  in->denom       = NULL;

} /* end ObitSpectrumInterpInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * For some reason this wasn't build into the GType class.
 * \param  in Pointer to the object to deallocate.
 *           Actually it should be an ObitSpectrumInterp* cast to an Obit*.
 */
void ObitSpectrumInterpClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitSpectrumInterp *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  in->info        = ObitInfoListUnref(in->info);
  in->Spectrum    = NULL;
  if (in->Freqs)  {g_free (in->Freqs);  in->Freqs=NULL;}
  if (in->chWt)   {g_free (in->chWt);   in->chWt=NULL;}
  if (in->Kernal) {g_free (in->Kernal); in->Kernal=NULL;}
  if (in->denom)  {g_free (in->denom);  in->denom=NULL;}

  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitSpectrumInterpClear */

/**
 * Set Lagrangian denominator.
 * \param  nFreq  Number of frequencies 
 * \param  hwidth Half width of convolution kernal
 * \param  Start first (0-rel) element in arrays
 * \param  Freqs Array of frequencies
 * \param  denom  [out] Reciprocals of Lagrangian denominators
 */
static void SetDenom (olong nFreq, olong hwidth, olong Start,
		      odouble *Freqs, ofloat *denom)
{
  ofloat prod;
  olong i, j, iwid;

   /* Init Lagrangian denominators for hwidth, start
   here assumed = nFreq */
  iwid = 1 + (2*hwidth);
  for (j=0; j<iwid; j++) {
    prod = 1.0;
    for (i=0; i<iwid; i++) {
      if (i != j) prod *= (ofloat)(Freqs[j+Start] - Freqs[i+Start]);
    } 
    denom[j] = 1.0 / prod;
  } 
  return;
}  /* end SetDenom */

/**
 * Set Lagrangian interpolation kernal taking into account ends of the grid.
 * \param  Freq   Desired Frequency of spectral point
 * \param  nFreq  Number of pixels on axis 
 * \param  Freqs  Frequency Array
 * \param  chWt   Channel weights
 * \param  hwidth Half width of convolution kernal
 * \param  denom  [out] Reciprocals of Lagrangian denominators
 * \param  Start  [out] first pixel in array for convolution
 * \param  cenFreq[out] Central (closest) Frequency
 * \param  Kernal [out] convolving kernal.
 */
static void SetConvKernal (double Freq, olong nFreq, odouble *Freqs, ofloat *chWt,
			   olong hwidth, ofloat *denom, olong *Start,
			   odouble *cenFreq, ofloat *Kernal)
{
  ofloat prod;
  olong i, j, iwid, newStart, best, lstart;
  odouble delta, F1, F2, F3, F4,F5, F6, F7, xx;

  /* Find closest frequency / start channel */
  iwid = 2*hwidth + 1;
  delta = fabs(Freq-Freqs[0]); best = 0;
  for (i=1; i<nFreq; i++) {
    if (fabs(Freq-Freqs[i])<delta) {
      delta = fabs(Freq-Freqs[i]); best = i;
    }
  }
  *cenFreq = Freqs[best];
  newStart = best - hwidth;
  /* Deal with edge cases */
  if (newStart<0) newStart = 0;
  if (newStart+iwid>nFreq) newStart = nFreq - iwid;
    
  if (newStart!=*Start) {
    *Start = newStart;
    /* Redo denom as needed */
    SetDenom (nFreq, hwidth, *Start, Freqs, denom);
  }

  xx = Freq-*cenFreq;  /* Relative frequency to interpolate */
  lstart = *Start;    /* Start frequency */
  F1 = xx-(Freqs[lstart+0]-*cenFreq);
  F2 = xx-(Freqs[lstart+1]-*cenFreq);
  F3 = xx-(Freqs[lstart+2]-*cenFreq);
  /* compute interpolating kernal by hwidth */
  switch (iwid) {
  case 3:
    Kernal[0] = denom[0]*(F2)*(F3);
    Kernal[1] = denom[1]*(F1)*(F3);
    Kernal[2] = denom[2]*(F1)*(F2);
    break;
  case 5:
    F4 = xx-(Freqs[lstart+3]-*cenFreq);
    F5 = xx-(Freqs[lstart+4]-*cenFreq);
    Kernal[0] = denom[0]*(F2)*(F3)*(F4)*(F5);
    Kernal[1] = denom[1]*(F1)*(F3)*(F4)*(F5);
    Kernal[2] = denom[2]*(F1)*(F2)*(F4)*(F5);
    Kernal[3] = denom[3]*(F1)*(F2)*(F3)*(F5);
    Kernal[4] = denom[4]*(F1)*(F2)*(F3)*(F4);
    break;
  case 7:
    F4 = xx-(Freqs[lstart+3]-*cenFreq);
    F5 = xx-(Freqs[lstart+4]-*cenFreq);
    F6 = xx-(Freqs[lstart+5]-*cenFreq);
    F7 = xx-(Freqs[lstart+6]-*cenFreq);
    Kernal[0] = denom[0]*(F2)*(F3)*(F4)*(F5)*(F6)*(F7);
    Kernal[1] = denom[1]*(F1)*(F3)*(F4)*(F5)*(F6)*(F7);
    Kernal[2] = denom[2]*(F1)*(F2)*(F4)*(F5)*(F6)*(F7);
    Kernal[3] = denom[3]*(F1)*(F2)*(F3)*(F5)*(F6)*(F7);
    Kernal[4] = denom[4]*(F1)*(F2)*(F3)*(F4)*(F6)*(F7);
    Kernal[5] = denom[5]*(F1)*(F2)*(F3)*(F4)*(F5)*(F7);
    Kernal[6] = denom[6]*(F1)*(F2)*(F3)*(F4)*(F5)*(F6);
    break;
  default:
    /* Anything else here */
    for (j= 0; j<iwid; j++) {
      prod = denom[j];
      for (i= 0; i<iwid; i++) {
	if (i != j) prod *= (xx - Freqs[lstart+i]);
      };
      Kernal[j] = prod;
    }
  }; /* end switch */
  return;
}  /* end SetConvKernal */

/**
 * Set Lagrangian interpolation kernal taking into account ends of the grid.
 * \param  Freq   Desired Frequency of spectral point
 * \param  nFreq  Number of pixels on axis 
 * \param  Freqs  Frequency Array
 * \param  chWt   Channel weights
 * \param  hwidth Half width of convolution kernal
 * \param  denom  [out] Reciprocals of Lagrangian denominators
 * \param  Start  [out] first pixel in array for convolution
 * \param  cenFreq[out] Central (closest) Frequency
 * \param  Kernal [out] convolving kernal.
 */
static void SetConvKernalY (double Freq, olong nFreq, odouble *Freqs, ofloat *chWt,
			   olong hwidth, ofloat *denom, olong *Start,
			   odouble *cenFreq, ofloat *Kernal)
{
  ofloat prod;
  olong i, j, iwid, newStart, best, lstart;
  odouble delta, F1, F2, F3, F4,F5, F6, F7, xx;

  /* Find closest frequency / start channel */
  iwid = 2*hwidth + 1;
  delta = fabs(Freq-Freqs[0]); best = 0;
  for (i=1; i<nFreq; i++) {
    if (fabs(Freq-Freqs[i])<delta) {
      delta = fabs(Freq-Freqs[i]); best = i;
    }
  }
  *cenFreq = Freqs[best];
  newStart = best - hwidth;
  /* Deal with edge cases */
  if (newStart<0) newStart = 0;
  if (newStart+iwid>nFreq) newStart = nFreq - iwid;
    
  if (newStart!=*Start) {
    *Start = newStart;
    /* Redo denom as needed */
    SetDenom (nFreq, hwidth, *Start, Freqs, denom);
  }

  xx = Freq-*cenFreq;  /* Relative frequency to interpolate */
  lstart = *Start;    /* Start frequency */
  F1 = Freqs[lstart+0]-*cenFreq;
  F2 = Freqs[lstart+1]-*cenFreq;
  F3 = Freqs[lstart+2]-*cenFreq;
  /* compute interpolating kernal by hwidth */
  switch (iwid) {
  case 3:
    Kernal[0] = denom[0]*(xx-F2)*(xx-F3);
    Kernal[1] = denom[1]*(xx-F1)*(xx-F3);
    Kernal[2] = denom[2]*(xx-F1)*(xx-F2);
    break;
  case 5:
    F4 = Freqs[lstart+3]-*cenFreq;
    F5 = Freqs[lstart+4]-*cenFreq;
    Kernal[0] = denom[0]*(xx-F2)*(xx-F3)*(xx-F4)*(xx-F5);
    Kernal[1] = denom[1]*(xx-F1)*(xx-F3)*(xx-F4)*(xx-F5);
    Kernal[2] = denom[2]*(xx-F1)*(xx-F2)*(xx-F4)*(xx-F5);
    Kernal[3] = denom[3]*(xx-F1)*(xx-F2)*(xx-F3)*(xx-F5);
    Kernal[4] = denom[4]*(xx-F1)*(xx-F2)*(xx-F3)*(xx-F4);
    break;
  case 7:
    F4 = Freqs[lstart+3]-*cenFreq;
    F5 = Freqs[lstart+4]-*cenFreq;
    F6 = Freqs[lstart+5]-*cenFreq;
    F7 = Freqs[lstart+6]-*cenFreq;
    Kernal[0] = denom[0]*(xx-F2)*(xx-F3)*(xx-F4)*(xx-F5)*(xx-F6)*(xx-F7);
    Kernal[1] = denom[1]*(xx-F1)*(xx-F3)*(xx-F4)*(xx-F5)*(xx-F6)*(xx-F7);
    Kernal[2] = denom[2]*(xx-F1)*(xx-F2)*(xx-F4)*(xx-F5)*(xx-F6)*(xx-F7);
    Kernal[3] = denom[3]*(xx-F1)*(xx-F2)*(xx-F3)*(xx-F5)*(xx-F6)*(xx-F7);
    Kernal[4] = denom[4]*(xx-F1)*(xx-F2)*(xx-F3)*(xx-F4)*(xx-F6)*(xx-F7);
    Kernal[5] = denom[5]*(xx-F1)*(xx-F2)*(xx-F3)*(xx-F4)*(xx-F5)*(xx-F7);
    Kernal[6] = denom[6]*(xx-F1)*(xx-F2)*(xx-F3)*(xx-F4)*(xx-F5)*(xx-F6);
    break;
  default:
    /* Anything else here */
    for (j= 0; j<iwid; j++) {
      prod = denom[j];
      for (i= 0; i<iwid; i++) {
	if (i != j) prod *= (xx - Freqs[lstart+i]);
      };
      Kernal[j] = prod;
    }
  }; /* end switch */
  return;
}  /* end SetConvKernalY */

