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

#include "ObitSpectrumMF.h"
#include "ObitSpectrumFit.h"
#include "ObitPosition.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitSpectrumMF.c
 * ObitSpectrumMF class function definitions.
 * This class is derived from the Obit base class.
 * This class  evaluates an MFImage-like tabulated spectrum.
 * Lagrange interpolation.
 */

/*--------------- File Global Variables  ----------------*/
/** name of the class defined in this file */
static gchar *myClassName = "ObitSpectrumMF";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo structure ObitSpectrumMFClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitSpectrumMFClassInfo myClassInfo = {FALSE};


/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitSpectrumMFInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitSpectrumMFClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitSpectrumMFClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitSpectrumMF* newObitSpectrumMF (gchar* name)
{
  ObitSpectrumMF* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitSpectrumMFClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitSpectrumMF));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitSpectrumMFInit((gpointer)out);

 return out;
} /* end newObitSpectrumMF */

/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name     An optional name for the object.
 * \param nFreq    Number of frequencies (subbands) in spectrum
 * \param Freqs    Array of frequencies in spectrum 
 * \param nTerm    Number of terms in spectral fit
 * \param refFreq  Reference frequency used in Spectral fitting
 * \return the new object.
 */
ObitSpectrumMF*
newObitSpectrumMFCreate (gchar* name, olong nFreq, odouble *Freqs,
			 odouble refFreq, olong nTerm)
{
  ObitSpectrumMF* out;
  olong i;

  /* Create/init output structure */
  out = newObitSpectrumMF (name);

  out->nFreq   = nFreq;   /* How many subbands */
  out->nTerm   = nTerm;   /* How many terms in fit */
  out->refFreq = refFreq; /* Reference frequency*/

  /* CreatewFrequency Array */
  out->Freqs  = g_malloc0(nFreq*sizeof(odouble));

  /* Copy Data */
  for (i=0; i<nFreq; i++)  out->Freqs[i] = Freqs[i];
 
  return out; 
} /* end newObitSpectrumMFCreate */

/**
 * Update info 
 * \param in       The object to update
 * \param nFreq    Number of frequencies (subbands) in spectrum
 * \param Freqs    Array of frequencies in spectrum 
 *                 updated if nonNULL
 * \param nTerm    Number of terms in spectral fit
 * \param refFreq  Reference frequency used in Spectral fitting
 */
void ObitSpectrumMFUpdate (ObitSpectrumMF *in, olong nFreq, odouble *Freqs,
			   odouble refFreq, olong nTerm)
{
  olong i;

  /*  Always update refFreq, nTerm */
  in->refFreq = refFreq; /* Reference frequency*/
  in->nTerm   = nTerm;   /* How many terms in fit */

  /* Freqs given? */
  if (Freqs!=NULL) {
    if ((in->Freqs) && (in->nFreq!=nFreq))  /* Delete old if wrong size */
      {g_free(in->Freqs); in->Freqs = NULL;}
    if (in->Freqs==NULL) in->Freqs = g_malloc0(nFreq*sizeof(odouble));
    in->nFreq   = nFreq;   /* How many subbands */
    /* Copy Data */
    for (i=0; i<nFreq; i++) in->Freqs[i] = Freqs[i];
  } /* end update Freqs */
} /* end ObitSpectrumMFUpDate */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitSpectrumMFGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitSpectrumMFClassInit();

  return (gconstpointer)&myClassInfo;
} /* end  ObitSpectrumMFGetClass */

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
ObitSpectrumMF*
ObitSpectrumMFCopy (ObitSpectrumMF *in, ObitSpectrumMF *out, ObitErr *err)
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
    out = newObitSpectrumMF(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /* copy/set this classes additions */
  out->nFreq      = in->nFreq;
  out->nTerm      = in->nTerm;
  out->refFreq    = in->refFreq;
  /* Copy arrays */
  if (out->Freqs) g_free(out->Freqs);
  out->Freqs  = g_malloc0(in->nFreq*sizeof(odouble));
  for (i=0; i<in->nFreq; i++) out->Freqs[i] = in->Freqs[i];
  
  return out;
} /* end ObitSpectrumMFCopy */

/**
 * Make a shallow copy of a object.
 * The result will have pointers to the more complex members.
 * Parent class members are included but any derived class info is ignored.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \return pointer to the new object.
 */
ObitSpectrumMF* ObitSpectrumMFClone  (ObitSpectrumMF *in, ObitSpectrumMF *out)
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
    out = newObitSpectrumMF(outName);
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
  out->nTerm      = in->nTerm;
  out->refFreq    = in->refFreq;

  /* Copy arrays */
  if (out->Freqs) g_free(out->Freqs);
  out->Freqs  = g_malloc0(in->nFreq*sizeof(odouble));
  for (i=0; i<in->nFreq; i++) out->Freqs[i] = in->Freqs[i];
  
  return out;
} /* end ObitSpectrumMFClone */

/**
 * Interpolate value at requested cell in a spectrum
 * Values out of the frequency range are extrapolations of the first/last two
 * \param in       The object to interpolate
 * \param chFreq   The frequency at which the spectral value is desired.
 * \param sbno     Closest subband number (0-rel)
 * \param Spectrum Subband spectrum; 
 * \param spFit    Spectral fit parameters from ObitSpectrumFitSingle(Arg)
 *                 to subband spectrum
 * \return spectral value, if the  nearest subband is 0, the spectrum evaluates at the 
 *         channel frequency is returned.
 */
ofloat ObitSpectrumMFEval (ObitSpectrumMF *in, odouble chFreq,
			   olong sbno, ofloat *Spectrum, ofloat* spFit)
{
  ofloat fblank =  ObitMagicF();
  ofloat flux = fblank;
  ofloat flux1, flux2, ratio;

  flux1 = ObitSpectrumEval (in->nTerm, in->refFreq, spFit,  in->Freqs[sbno]);  /* Subband Freq */
  flux2 = ObitSpectrumEval (in->nTerm, in->refFreq, spFit,  chFreq);           /* channel freq */
  if (flux1>0.0) ratio = flux2/flux1;
  else           ratio = 1.0;
  if (Spectrum[sbno]!=0.0) flux = ratio*Spectrum[sbno];
  else                     flux = flux2;
  return flux;
} /* end  ObitSpectrumMFEval */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitSpectrumMFClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitSpectrumMFClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitSpectrumMFClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitSpectrumMFClassInfoDefFn (gpointer inClass)
{
  ObitSpectrumMFClassInfo *theClass = (ObitSpectrumMFClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitSpectrumMFClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitSpectrumMFClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitSpectrumMFGetClass;
  theClass->newObit       = (newObitFP)newObitSpectrumMF;
  theClass->ObitCopy      = (ObitCopyFP)ObitSpectrumMFCopy;
  theClass->ObitClone     = (ObitCloneFP)ObitSpectrumMFClone;
  theClass->ObitClear     = (ObitClearFP)ObitSpectrumMFClear;
  theClass->ObitInit      = (ObitInitFP)ObitSpectrumMFInit;

} /* end ObitSpectrumMFClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param in Pointer to the object to initialize.
 */
void ObitSpectrumMFInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitSpectrumMF *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->info        = newObitInfoList();
  in->nFreq       = -1;
  in->nTerm       = -1;
  in->Freqs       = NULL;
  in->refFreq     = -1.0;
} /* end ObitSpectrumMFInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * For some reason this wasn't build into the GType class.
 * \param  in Pointer to the object to deallocate.
 *           Actually it should be an ObitSpectrumMF* cast to an Obit*.
 */
void ObitSpectrumMFClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitSpectrumMF *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  in->info        = ObitInfoListUnref(in->info);
  if (in->Freqs)  {g_free (in->Freqs);  in->Freqs=NULL;}

  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitSpectrumMFClear */
