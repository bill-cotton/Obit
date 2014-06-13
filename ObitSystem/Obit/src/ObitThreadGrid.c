/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2014                                               */
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

#include "ObitThread.h"
#include "ObitThreadGrid.h"
#include "ObitSinCos.h"
#include "ObitExp.h"
#include "ObitUVGridMF.h"
#include "ObitUVGridWB.h"
#include <stdio.h>
#include <math.h>

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitThreadGrid.c
 * ObitThreadGrid class function definitions.
 * This class is derived from the Obit base class.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitThreadGrid";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/*--------------- File Global Variables  ----------------*/
/**
 * ClassInfo structure ObitThreadGridClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitThreadGridClassInfo myClassInfo = {FALSE};

/*---------------Private structures----------------*/
/* threaded function argument */
typedef struct {
  /* facet to grid */
  olong        facet;
  /* vrow to flip conjugate sides */
  olong        vrow;
  /* highest vis (0-rel) */
  olong        hivis;
  /* lowest vis (0-rel) */
  olong        lovis;
  /* thread number, >0 -> no threading  */
  olong        ithread;
  /* grid i (0-rel) or n for image/beam  */
  olong        iGpI, nGpI;
  /* (SW) Beam order 0=normal */
  olong        beamOrd;
  /* Channel selection  */
  olong bChan, eChan;
  /* uv data buffer */
  ofloat       *data;
  /* uv grid */
  ofloat       *grid;
  /* gridding info */
  ObitThreadGridInfo *gridInfo;
  /* output merged/swapped grid */
  ObitCArray *outGrid;
  /* thread */
  ObitThread *thread;
  /* Arrays gridding tapers per uv data channel for forcing the beams */
  ofloat *sigma1, *sigma2, *sigma3;
  /* work arrays */
  ofloat  *fwork1, *fwork2, *fwork3, **cnvfnu, **cnvfnv;
  olong   *iuarr, *ivarr;
  /* debugging array
  ofloat *debug; */
} GridFuncArg;

/*----------------------Private functions---------------------------*/
/** Private: Initialize newly instantiated object. */
void  ObitThreadGridInit  (gpointer in);
/** Private: Deallocate members. */
void  ObitThreadGridClear (gpointer in);
/** Private: Set Class function pointers. */
static void ObitThreadGridClassInfoDefFn (gpointer inClass);
/** Base uv grid setup */
void ObitThreadGridSetupBase (ObitThreadGrid *in, ObitUV *UVin,
			      olong nPar, ObitUVGrid **UVGrids, 
			      olong nThreads, ObitErr *err);
/** Multifrequency uv grid setup */
void ObitThreadGridSetupMF (ObitThreadGrid *in, ObitUV *UVin,
			    olong nPar, ObitUVGrid **UVGrids, 
			    olong nThreads, ObitErr *err);
/** Wideband (Sault-Weirenga) */
void ObitThreadGridSetupWB (ObitThreadGrid *in, ObitUV *UVin,
			    olong nPar, ObitUVGrid **UVGrids, 
			    olong nThreads, ObitErr *err);
/** Private: Grid a data buffer. */
static gpointer ThreadGrid (gpointer args);
/** Private: Add neg u to conjugate cells */
static gpointer ThreadFlip (gpointer args);
/** Private: swap/merge grids */
static gpointer ThreadMerge (gpointer args);
//** Private: pre data for gridding routine */
void fast_prep_grid(olong ivis, GridFuncArg *args);
/** Private: inner gridding routine */
void fast_grid(ofloat *grid, ofloat vis[2], olong iu, olong iv, olong lrow, 
	       olong nconv, ofloat *cu, ofloat *cv);
/** Private: inner 7x7 (AVX) gridding routine */
void fast_grid7(ofloat *grid, ofloat vis[2], olong iu, olong iv, olong lrow, 
		olong nconv, ofloat *cu, ofloat *cv);
/* Round */
static olong _lroundf(ofloat in)
{
  olong out;
  if (in>=0.0) out = (olong)(in+0.5);
  else         out = (olong)(in-0.5);
  return out;
}
/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitThreadGrid* newObitThreadGrid (gchar* name)
{
  ObitThreadGrid* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitThreadGridClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitThreadGrid));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitThreadGridInit((gpointer)out);

 return out;
} /* end newObitThreadGrid */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitThreadGridGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitThreadGridClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitThreadGridGetClass */

/**
 * Prepares for gridding uv data of the type described by UVGrids
 * \param in       Object to initialize
 * \param UVin     Data to be gridded
 * \param nPar     Number of parallel grids 
 * \param UVGrids  Array of ObitUVGrid object to be gridded.
 * \param nThreads Number of threads to use
 * \param err      ObitErr stack for reporting problems.
 */
void ObitThreadGridSetup (ObitThreadGrid *in, ObitUV *UVin,
			  olong nPar, ObitUVGrid **UVGrids, 
			  olong nThreads, ObitErr *err)
{
 gchar *routine="ObitThreadGridSetup";

 /* Branch by gridder type */
 if (ObitUVGridMFIsA(UVGrids[0])) {
   /* Multi Frequency version */
   ObitThreadGridSetupMF (in, UVin, nPar, UVGrids, nThreads, err);
 } else if (ObitUVGridWBIsA(UVGrids[0])) {
   /* Wideband (SW) version */
   ObitThreadGridSetupWB (in, UVin, nPar, UVGrids, nThreads, err);
 } else {
   /* Base version */
   ObitThreadGridSetupBase (in, UVin, nPar, UVGrids, nThreads, err);
 }
 if (err->error) Obit_traceback_msg (err, routine, in->name);
}  /* end ObitThreadGridSetup */

/**
 * Prepares for gridding uv data of the type described by UVGrids
 * Base gridder version
 * \param in       Object to initialize
 * \param UVin     Data to be gridded
 * \param nPar     Number of parallel grids 
 * \param UVGrids  Array of ObitUVGrid object to be gridded.
 * \param nThreads Number of threads to use
 * \param err      ObitErr stack for reporting problems.
 */
void ObitThreadGridSetupBase (ObitThreadGrid *in, ObitUV *UVin,
			      olong nPar, ObitUVGrid **UVGrids, 
			      olong nThreads, ObitErr *err)
{
  ObitUVDesc *uvDesc;
  GridFuncArg **funcarg;
  olong i, j, k, iGrid, iTh, wrksize, nGrid, nGpI, size;
  olong nvis, nvisPth, nChan;
  ofloat noRot[] = {1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0};
  /*gchar *routine="ObitThreadGridSetup";*/

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitThreadGridIsA(in));
  g_assert (ObitUVGridIsA(UVGrids[0]));

   /* Init sincos, exp */
  ObitSinCosInit ();
  ObitExpInit();

  uvDesc   = UVin->myDesc;     /* UV descriptor */
  in->UVin = ObitUVRef(UVin);  /* Save UV data */

  /* How many grids per image - make best use of threads */
  if (nPar>nThreads) nGpI = 1;
  else               nGpI = (olong)(0.999+((ofloat)nThreads)/nPar);
  /* Try to split grid over at least 4 threads */
  if (nThreads>=4) nGpI = MAX (4, nGpI);
  /* How many grids?  */
  nGrid     = nGpI*nPar;
  in->nGrid = nGrid;   /* no. facets gridded = ngrid */
  in->nGpI  = nGpI;

  /* How many vis, per thread? */
  nvis    = uvDesc->numVisBuff;
  nvisPth = nvis/nGpI;
  /* How many channel/IFs? */
  nChan      = uvDesc->inaxes[uvDesc->jlocf] * uvDesc->inaxes[uvDesc->jlocif];

  /* structure arrays */
  in->GridInfo              = g_malloc0(sizeof(ObitThreadGridInfo));
  in->GridInfo->guardu      = g_malloc0(nGrid*sizeof(ofloat));
  in->GridInfo->guardv      = g_malloc0(nGrid*sizeof(ofloat));
  in->GridInfo->maxBL       = g_malloc0(nGrid*sizeof(ofloat));
  in->GridInfo->minBL       = g_malloc0(nGrid*sizeof(ofloat));
  in->GridInfo->nx          = g_malloc0(nGrid*sizeof(olong));
  in->GridInfo->ny          = g_malloc0(nGrid*sizeof(olong));
  in->GridInfo->uscale      = g_malloc0(nGrid*sizeof(ofloat));
  in->GridInfo->vscale      = g_malloc0(nGrid*sizeof(ofloat));
  in->GridInfo->BeamTaperUV = g_malloc0(nGrid*sizeof(ofloat));
  in->GridInfo->shift       = g_malloc0(nGrid*3*sizeof(ofloat));
  in->GridInfo->rotUV       = g_malloc0(nGrid*9*sizeof(ofloat));
  in->GridInfo->isBeam      = g_malloc0(nGrid*sizeof(gboolean));
  in->GridInfo->thArgs      = g_malloc0(nGrid*sizeof(gpointer));
  /* Fill in global values */
  in->nGrid                  = nGrid;
  in->nGpI                   = nGpI;
  in->GridInfo->nfacet       = nGrid;
  in->GridInfo->nThreads     = nThreads;
  in->GridInfo->nchan        = nChan;
  in->GridInfo->nrparm       = uvDesc->nrparm;
  in->GridInfo->nvis         = nvis;
  in->GridInfo->lenvis       = uvDesc->lrec;
  in->GridInfo->freqArr      = uvDesc->fscale;
  in->GridInfo->convWidth    = UVGrids[0]->convWidth;
  in->GridInfo->convNperCell = UVGrids[0]->convNperCell;
  in->GridInfo->convfn       = UVGrids[0]->convgfn->array;
  in->GridInfo->thread       = ObitThreadRef(UVGrids[0]->thread);

  /* Initialize Function args */
  funcarg    = (GridFuncArg**)in->GridInfo->thArgs;
  wrksize    = nChan; /* Number of channels */
  for (iTh=0; iTh<nGrid; iTh++) {
    funcarg[iTh] = g_malloc(sizeof(GridFuncArg));
    funcarg[iTh]->ithread  = iTh;
    funcarg[iTh]->thread   = ObitThreadRef(UVGrids[0]->thread);
    funcarg[iTh]->vrow     = 0;
    funcarg[iTh]->facet    = iTh;
    funcarg[iTh]->data     = UVin->buffer;
    funcarg[iTh]->gridInfo = in->GridInfo;
    funcarg[iTh]->sigma1   = NULL;
    funcarg[iTh]->sigma2   = NULL;
    funcarg[iTh]->sigma3   = NULL;
    funcarg[iTh]->fwork1   = g_malloc0(2*wrksize*sizeof(ofloat));
    funcarg[iTh]->fwork2   = g_malloc0(2*wrksize*sizeof(ofloat));
    funcarg[iTh]->fwork3   = g_malloc0(wrksize*sizeof(ofloat));
    funcarg[iTh]->cnvfnu   = g_malloc0(wrksize*sizeof(ofloat*));
    funcarg[iTh]->cnvfnv   = g_malloc0(wrksize*sizeof(ofloat*));
    funcarg[iTh]->iuarr    = g_malloc0(wrksize*sizeof(olong));
    funcarg[iTh]->ivarr    = g_malloc0(wrksize*sizeof(olong));
    funcarg[iTh]->beamOrd  = 0;
   } /* End loop creating Thread args */

  /* order grids first by beam/non beam, then in order of UVGrids */
  iGrid = 0;
  for (i=0; i<nPar; i++) {
    if (UVGrids[i]->doBeam) {
      /* add gridding facets */
      for (j=0; j<nGpI; j++) {
	in->GridInfo->isBeam[iGrid]       = TRUE;         /* Is a beam */
	in->GridInfo->BeamTaperUV[iGrid]  = UVGrids[i]->BeamTaperUV;
	in->GridInfo->uscale[iGrid]       = UVGrids[i]->UScale;
	in->GridInfo->vscale[iGrid]       = UVGrids[i]->VScale;
	in->GridInfo->shift[iGrid*3+0]    = UVGrids[i]->dxc;
	in->GridInfo->shift[iGrid*3+1]    = UVGrids[i]->dyc;
	in->GridInfo->shift[iGrid*3+2]    = UVGrids[i]->dzc;
	in->GridInfo->nx[iGrid]           = UVGrids[i]->nxBeam;
	in->GridInfo->ny[iGrid]           = UVGrids[i]->nyBeam;
	in->GridInfo->guardu[iGrid]       = UVGrids[i]->guardband*UVGrids[i]->nxBeam;
	in->GridInfo->guardv[iGrid]       = UVGrids[i]->guardband*UVGrids[i]->nyBeam;
	in->GridInfo->maxBL[iGrid]        = UVGrids[i]->blmax;
	in->GridInfo->minBL[iGrid]        = UVGrids[i]->blmin;
	/* u,v,w rotation if do3Dmul */
	if (UVGrids[i]->do3Dmul) {
	  for (k=0; k<3; k++) in->GridInfo->rotUV[iGrid*9+k  ] = UVGrids[i]->URot3D[0][k];
	  for (k=0; k<3; k++) in->GridInfo->rotUV[iGrid*9+3+k] = UVGrids[i]->URot3D[1][k];
	  for (k=0; k<3; k++) in->GridInfo->rotUV[iGrid*9+6+k] = UVGrids[i]->URot3D[2][k];
	} else {  /* identity matrix */
	  for (k=0; k<9; k++) in->GridInfo->rotUV[iGrid*9+k] = noRot[k];
	}
	/* Thread arg stuff */
	funcarg[iGrid]->iGpI   = j;       /* grid number for this beam */
	funcarg[iGrid]->nGpI   = nGpI;    /* number of grids for this beam */
	funcarg[iGrid]->facet  = iGrid;	  /* facet (beam) number */
	funcarg[iGrid]->lovis  = j*nvisPth;	
	funcarg[iGrid]->hivis  = (j+1)*nvisPth;	
	size = 2 * (1 + UVGrids[i]->convWidth/2 + UVGrids[i]->nxBeam/2) * 
	  UVGrids[i]->nyBeam;
	funcarg[iGrid]->grid = g_malloc0(size*sizeof(ofloat));
	funcarg[iGrid]->outGrid = ObitCArrayRef(UVGrids[i]->grid);
	funcarg[iGrid]->beamOrd = 0;  /* (SW) beam order */
	/* Channel selection */
	funcarg[iGrid]->bChan = 0;
	funcarg[iGrid]->eChan = nChan;
	iGrid++;
      } /* end loop over grids per image */
      funcarg[iGrid-1]->hivis = MAX(funcarg[iGrid-1]->hivis, nvis);
    } /* end if beam */
  } /* End loop adding beams */
 
 /* Loop over images */
  for (i=0; i<nPar; i++) {
    if (!UVGrids[i]->doBeam) {
      /* add gridding facets */
      for (j=0; j<nGpI; j++) {
	in->GridInfo->isBeam[iGrid]       = FALSE;         /* Not a beam */
	in->GridInfo->BeamTaperUV[iGrid]  = UVGrids[i]->BeamTaperUV;
	in->GridInfo->uscale[iGrid]       = UVGrids[i]->UScale;
	in->GridInfo->vscale[iGrid]       = UVGrids[i]->VScale;
	in->GridInfo->shift[iGrid*3+0]    = UVGrids[i]->dxc;
	in->GridInfo->shift[iGrid*3+1]    = UVGrids[i]->dyc;
	in->GridInfo->shift[iGrid*3+2]    = UVGrids[i]->dzc;
	in->GridInfo->nx[iGrid]           = UVGrids[i]->nxImage;
	in->GridInfo->ny[iGrid]           = UVGrids[i]->nyImage;
	in->GridInfo->guardu[iGrid]       = UVGrids[i]->guardband*UVGrids[i]->nxImage;
	in->GridInfo->guardv[iGrid]       = UVGrids[i]->guardband*UVGrids[i]->nyImage;
	in->GridInfo->maxBL[iGrid]        = UVGrids[i]->blmax;
	in->GridInfo->minBL[iGrid]        = UVGrids[i]->blmin;
	/* u,v,w rotation if do3Dmul */
	if (UVGrids[i]->do3Dmul) {
	  for (k=0; k<3; k++) in->GridInfo->rotUV[iGrid*9+k  ] = UVGrids[i]->URot3D[0][k];
	  for (k=0; k<3; k++) in->GridInfo->rotUV[iGrid*9+3+k] = UVGrids[i]->URot3D[1][k];
	  for (k=0; k<3; k++) in->GridInfo->rotUV[iGrid*9+6+k] = UVGrids[i]->URot3D[2][k];
	} else {  /* identity matrix */
	  for (k=0; k<9; k++) in->GridInfo->rotUV[iGrid*9+k] = noRot[k];
	}
	/* Thread arg stuff */
	funcarg[iGrid]->iGpI   = j;       /* grid number for this image */
	funcarg[iGrid]->nGpI   = nGpI;    /* number of grids for this image */
	funcarg[iGrid]->facet  = iGrid;	  /* facet (image) number */
	funcarg[iGrid]->lovis  = j*nvisPth;	
	funcarg[iGrid]->hivis  = (j+1)*nvisPth;
	size = 2 * (1 + UVGrids[i]->convWidth/2 + UVGrids[i]->nxImage/2) * 
	  UVGrids[i]->nyImage;
	funcarg[iGrid]->grid = g_malloc0(size*sizeof(ofloat));
	funcarg[iGrid]->outGrid = ObitCArrayRef(UVGrids[i]->grid);
	/* Channel selection */
	funcarg[iGrid]->bChan = 0;
	funcarg[iGrid]->eChan = nChan;
	iGrid++;
      } /* end loop over grids per image */
      funcarg[iGrid-1]->hivis = MAX(funcarg[iGrid-1]->hivis, nvis);
    } /* end if image */
  } /* End loop adding image */
 
}  /* end ObitThreadGridSetupBase */

/**
 * Prepares for gridding uv data of the type described by UVGrids
 * Multifrequency gridder version
 * \param in       Object to initialize
 * \param UVin     Data to be gridded
 * \param nPar     Number of parallel grids 
 * \param UVGrids  Array of ObitUVGrid object to be gridded.
 * \param nThreads Number of threads to use
 * \param err      ObitErr stack for reporting problems.
 */
void ObitThreadGridSetupMF (ObitThreadGrid *in, ObitUV *UVin,
			    olong nPar, ObitUVGrid **UVGridss, 
			    olong nThreads, ObitErr *err)
{
  ObitUVDesc *uvDesc;
  ObitUVGridMF **UVGrids = (ObitUVGridMF**)UVGridss;
  GridFuncArg **funcarg;
  olong i, j, k, iGrid, iTh, wrksize, nGrid, nGpI, size;
  olong nvis, nvisPth, nChan, bCh, eCh, iSpec;
  ofloat noRot[] = {1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0};
  /*gchar *routine="ObitThreadGridSetup";*/

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitThreadGridIsA(in));
  g_assert (ObitUVGridIsA(UVGrids[0]));

   /* Init sincos, exp */
  ObitSinCosInit ();
  ObitExpInit();

  uvDesc   = UVin->myDesc;     /* UV descriptor */
  in->UVin = ObitUVRef(UVin);  /* Save UV data */

  /* How many grids per image - make best use of threads */
  if ((nPar*UVGrids[0]->nSpec)>nThreads) nGpI = 1;
  else               nGpI = (olong)(0.999+((ofloat)nThreads)/
				    (nPar*UVGrids[0]->nSpec));
  /* Try to split grid over at least 4 threads */
  if (nThreads>=4) nGpI = MAX (4/UVGrids[0]->nSpec, nGpI);
  /* How many grids?  one per coarse frequency planes */
  nGrid     = nGpI*nPar*UVGrids[0]->nSpec;
  in->nGrid = nGrid;   /* no. facets gridded = ngrid */
  in->nGpI  = nGpI;

  /* How many vis, per thread? */
  nvis    = uvDesc->numVisBuff;
  nvisPth = nvis/nGpI;
  /* How many channels */
  nChan      = uvDesc->inaxes[uvDesc->jlocf];

  /* structure arrays */
  in->GridInfo              = g_malloc0(sizeof(ObitThreadGridInfo));
  in->GridInfo->guardu      = g_malloc0(nGrid*sizeof(ofloat));
  in->GridInfo->guardv      = g_malloc0(nGrid*sizeof(ofloat));
  in->GridInfo->maxBL       = g_malloc0(nGrid*sizeof(ofloat));
  in->GridInfo->minBL       = g_malloc0(nGrid*sizeof(ofloat));
  in->GridInfo->nx          = g_malloc0(nGrid*sizeof(olong));
  in->GridInfo->ny          = g_malloc0(nGrid*sizeof(olong));
  in->GridInfo->uscale      = g_malloc0(nGrid*sizeof(ofloat));
  in->GridInfo->vscale      = g_malloc0(nGrid*sizeof(ofloat));
  in->GridInfo->BeamTaperUV = g_malloc0(nGrid*sizeof(ofloat));
  in->GridInfo->shift       = g_malloc0(nGrid*3*sizeof(ofloat));
  in->GridInfo->rotUV       = g_malloc0(nGrid*9*sizeof(ofloat));
  in->GridInfo->isBeam      = g_malloc0(nGrid*sizeof(gboolean));
  in->GridInfo->thArgs      = g_malloc0(nGrid*sizeof(gpointer));
  /* Fill in global values */
  in->nGrid                  = nGrid;
  in->nGpI                   = nGpI;
  in->GridInfo->nfacet       = nGrid;
  in->GridInfo->nThreads     = nThreads;
  in->GridInfo->nchan        = nChan * uvDesc->inaxes[uvDesc->jlocif];
  in->GridInfo->nrparm       = uvDesc->nrparm;
  in->GridInfo->nvis         = nvis;
  in->GridInfo->lenvis       = uvDesc->lrec;
  in->GridInfo->freqArr      = uvDesc->fscale;
  in->GridInfo->convWidth    = UVGrids[0]->convWidth;
  in->GridInfo->convNperCell = UVGrids[0]->convNperCell;
  in->GridInfo->convfn       = UVGrids[0]->convgfn->array;
  in->GridInfo->thread       = ObitThreadRef(UVGrids[0]->thread);

  /* Initialize Function args */
  funcarg    = (GridFuncArg**)in->GridInfo->thArgs;
  wrksize    = nChan* uvDesc->inaxes[uvDesc->jlocif]; /* Number of channels/IF  */
  for (iTh=0; iTh<nGrid; iTh++) {
    funcarg[iTh] = g_malloc(sizeof(GridFuncArg));
    funcarg[iTh]->ithread  = iTh;
    funcarg[iTh]->thread   = ObitThreadRef(UVGrids[0]->thread);
    funcarg[iTh]->vrow     = 0;
    funcarg[iTh]->facet    = iTh;
    funcarg[iTh]->data     = UVin->buffer;
    funcarg[iTh]->gridInfo = in->GridInfo;
    funcarg[iTh]->fwork1   = g_malloc0(2*wrksize*sizeof(ofloat));
    funcarg[iTh]->fwork2   = g_malloc0(2*wrksize*sizeof(ofloat));
    funcarg[iTh]->fwork3   = g_malloc0(wrksize*sizeof(ofloat));
    funcarg[iTh]->cnvfnu   = g_malloc0(wrksize*sizeof(ofloat*));
    funcarg[iTh]->cnvfnv   = g_malloc0(wrksize*sizeof(ofloat*));
    funcarg[iTh]->iuarr    = g_malloc0(wrksize*sizeof(olong));
    funcarg[iTh]->ivarr    = g_malloc0(wrksize*sizeof(olong));
    funcarg[iTh]->beamOrd  = 0;
   } /* End loop creating Thread args */

  /* order grids first by beam/non beam, then in order of UVGrids */
  iGrid = 0;
  for (i=0; i<nPar; i++) {
    if (UVGrids[i]->doBeam) {
      /* add gridding facets */
      for (j=0; j<nGpI; j++) {
	/* Loop over coarse frequency bin grids  */
	for (iSpec=0; iSpec<UVGrids[i]->nSpec; iSpec++) {
	  in->GridInfo->isBeam[iGrid]       = TRUE;         /* Is a beam */
	  in->GridInfo->BeamTaperUV[iGrid]  = UVGrids[i]->BeamTaperUV;
	  in->GridInfo->uscale[iGrid]       = UVGrids[i]->UScale;
	  in->GridInfo->vscale[iGrid]       = UVGrids[i]->VScale;
	  in->GridInfo->shift[iGrid*3+0]    = UVGrids[i]->dxc;
	  in->GridInfo->shift[iGrid*3+1]    = UVGrids[i]->dyc;
	  in->GridInfo->shift[iGrid*3+2]    = UVGrids[i]->dzc;
	  in->GridInfo->nx[iGrid]           = UVGrids[i]->nxBeam;
	  in->GridInfo->ny[iGrid]           = UVGrids[i]->nyBeam;
	  in->GridInfo->guardu[iGrid]       = UVGrids[i]->guardband*UVGrids[i]->nxBeam;
	  in->GridInfo->guardv[iGrid]       = UVGrids[i]->guardband*UVGrids[i]->nyBeam;
	  in->GridInfo->maxBL[iGrid]        = UVGrids[i]->blmax;
	  in->GridInfo->minBL[iGrid]        = UVGrids[i]->blmin;
	  /* u,v,w rotation if do3Dmul */
	  if (UVGrids[i]->do3Dmul) {
	    for (k=0; k<3; k++) in->GridInfo->rotUV[iGrid*9+k  ] = UVGrids[i]->URot3D[0][k];
	    for (k=0; k<3; k++) in->GridInfo->rotUV[iGrid*9+3+k] = UVGrids[i]->URot3D[1][k];
	    for (k=0; k<3; k++) in->GridInfo->rotUV[iGrid*9+6+k] = UVGrids[i]->URot3D[2][k];
	  } else {  /* identity matrix */
	    for (k=0; k<9; k++) in->GridInfo->rotUV[iGrid*9+k] = noRot[k];
	  }
	  /* Thread arg stuff */
	  funcarg[iGrid]->iGpI   = j;       /* grid number for this beam */
	  funcarg[iGrid]->nGpI   = nGpI;    /* number of grids for this beam */
	  funcarg[iGrid]->facet  = iGrid;	  /* facet (beam) number */
	  funcarg[iGrid]->lovis  = j*nvisPth;	
	  funcarg[iGrid]->hivis  = (j+1)*nvisPth;	
	  size = 2 * (1 + UVGrids[i]->convWidth/2 + UVGrids[i]->nxBeam/2) * 
	    UVGrids[i]->nyBeam;
	  funcarg[iGrid]->grid    = g_malloc0(size*sizeof(ofloat));
	  funcarg[iGrid]->outGrid = ObitCArrayRef(UVGrids[i]->grids[iSpec]);
	  funcarg[iGrid]->beamOrd = 0;  /* (SW) beam order */
	  /* Channel selection */
	  bCh = UVGrids[i]->BIFSpec[iSpec]*nChan + UVGrids[i]->BChanSpec[iSpec];
	  eCh = UVGrids[i]->EIFSpec[iSpec]*nChan + UVGrids[i]->EChanSpec[iSpec];
	  funcarg[iGrid]->bChan = bCh;
	  funcarg[iGrid]->eChan = eCh;
	  /* beam Tapering parameters */
	  funcarg[iGrid]->sigma1   = UVGrids[i]->sigma1;
	  funcarg[iGrid]->sigma2   = UVGrids[i]->sigma2;
	  funcarg[iGrid]->sigma3   = UVGrids[i]->sigma3;
	  iGrid++;
	} /* end loop over coarse frequency bin grids */
      } /* end loop over grids per image */
      funcarg[iGrid-1]->hivis = MAX(funcarg[iGrid-1]->hivis, nvis);
    } /* end if beam */
  } /* End loop adding beams */
 
  /* Loop over images */
  for (i=0; i<nPar; i++) {
    if (!UVGrids[i]->doBeam) {
      /* add gridding facets */
      for (j=0; j<nGpI; j++) {
	/* Loop over coarse frequency bin grids  */
	for (iSpec=0; iSpec<UVGrids[i]->nSpec; iSpec++) {
	  in->GridInfo->isBeam[iGrid]       = FALSE;         /* Not a beam */
	  in->GridInfo->BeamTaperUV[iGrid]  = UVGrids[i]->BeamTaperUV;
	  in->GridInfo->uscale[iGrid]       = UVGrids[i]->UScale;
	  in->GridInfo->vscale[iGrid]       = UVGrids[i]->VScale;
	  in->GridInfo->shift[iGrid*3+0]    = UVGrids[i]->dxc;
	  in->GridInfo->shift[iGrid*3+1]    = UVGrids[i]->dyc;
	  in->GridInfo->shift[iGrid*3+2]    = UVGrids[i]->dzc;
	  in->GridInfo->nx[iGrid]           = UVGrids[i]->nxImage;
	  in->GridInfo->ny[iGrid]           = UVGrids[i]->nyImage;
	  in->GridInfo->guardu[iGrid]       = UVGrids[i]->guardband*UVGrids[i]->nxImage;
	  in->GridInfo->guardv[iGrid]       = UVGrids[i]->guardband*UVGrids[i]->nyImage;
	  in->GridInfo->maxBL[iGrid]        = UVGrids[i]->blmax;
	  in->GridInfo->minBL[iGrid]        = UVGrids[i]->blmin;
	  /* u,v,w rotation if do3Dmul */
	  if (UVGrids[i]->do3Dmul) {
	    for (k=0; k<3; k++) in->GridInfo->rotUV[iGrid*9+k  ] = UVGrids[i]->URot3D[0][k];
	    for (k=0; k<3; k++) in->GridInfo->rotUV[iGrid*9+3+k] = UVGrids[i]->URot3D[1][k];
	    for (k=0; k<3; k++) in->GridInfo->rotUV[iGrid*9+6+k] = UVGrids[i]->URot3D[2][k];
	  } else {  /* identity matrix */
	    for (k=0; k<9; k++) in->GridInfo->rotUV[iGrid*9+k] = noRot[k];
	  }
	  /* Thread arg stuff */
	  funcarg[iGrid]->iGpI   = j;       /* grid number for this image */
	  funcarg[iGrid]->nGpI   = nGpI;    /* number of grids for this image */
	  funcarg[iGrid]->facet  = iGrid;	  /* facet (image) number */
	  funcarg[iGrid]->lovis  = j*nvisPth;	
	  funcarg[iGrid]->hivis  = (j+1)*nvisPth;
	  size = 2 * (1 + UVGrids[i]->convWidth/2 + UVGrids[i]->nxImage/2) * 
	    UVGrids[i]->nyImage;
	  funcarg[iGrid]->grid = g_malloc0(size*sizeof(ofloat));
	  funcarg[iGrid]->outGrid = ObitCArrayRef(UVGrids[i]->grids[iSpec]);
	  /* Channel selection */
	  bCh = UVGrids[i]->BIFSpec[iSpec]*nChan + UVGrids[i]->BChanSpec[iSpec];
	  eCh = UVGrids[i]->EIFSpec[iSpec]*nChan + UVGrids[i]->EChanSpec[iSpec];
	  funcarg[iGrid]->bChan = bCh;
	  funcarg[iGrid]->eChan = eCh;
	  /* beam Tapering parameters */
	  funcarg[iGrid]->sigma1   = UVGrids[i]->sigma1;
	  funcarg[iGrid]->sigma2   = UVGrids[i]->sigma2;
	  funcarg[iGrid]->sigma3   = UVGrids[i]->sigma3;
	  iGrid++;
	} /* end loop over coarse frequency bin grids */
      } /* end loop over grids per image */
      funcarg[iGrid-1]->hivis = MAX(funcarg[iGrid-1]->hivis, nvis);
    } /* end if image */
  } /* End loop adding image */
}  /* end ObitThreadGridSetupMF */

/**
 * Prepares for gridding uv data of the type described by UVGrids
 * WB (SW) gridder version
 * \param in       Object to initialize
 * \param UVin     Data to be gridded
 * \param nPar     Number of parallel grids 
 * \param UVGrids  Array of ObitUVGrid object to be gridded.
 * \param nThreads Number of threads to use
 * \param err      ObitErr stack for reporting problems.
 */
void ObitThreadGridSetupWB (ObitThreadGrid *in, ObitUV *UVin,
			    olong nPar, ObitUVGrid **UVGridss, 
			    olong nThreads, ObitErr *err)
{
  ObitUVGridWB **UVGrids = (ObitUVGridWB**)UVGridss;
  ObitUVDesc *uvDesc;
  GridFuncArg **funcarg;
  olong i, j, k, iGrid, iTh, wrksize, nGrid, nGpI, size;
  olong nvis, nvisPth, nChan;
  ofloat noRot[] = {1.0,0.0,0.0, 0.0,1.0,0.0, 0.0,0.0,1.0};
  /*gchar *routine="ObitThreadGridSetup";*/

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitThreadGridIsA(in));
  g_assert (ObitUVGridIsA(UVGrids[0]));

   /* Init sincos, exp */
  ObitSinCosInit ();
  ObitExpInit();

  uvDesc   = UVin->myDesc;     /* UV descriptor */
  in->UVin = ObitUVRef(UVin);  /* Save UV data */

  /* How many grids per image - make best use of threads */
  if (nPar>nThreads) nGpI = 1;
  else               nGpI = (olong)(0.999+((ofloat)nThreads)/nPar);
  /* Try to split grid over at least 4 threads */
  if (nThreads>=4) nGpI = MAX (4, nGpI);
  /* How many grids?  */
  nGrid     = nGpI*nPar;
  in->nGrid = nGrid;   /* no. facets gridded = ngrid */
  in->nGpI  = nGpI;

  /* How many vis, per thread? */
  nvis    = uvDesc->numVisBuff;
  nvisPth = nvis/nGpI;
  /* How many channel/IFs? */
  nChan      = uvDesc->inaxes[uvDesc->jlocf] * uvDesc->inaxes[uvDesc->jlocif];

  /* structure arrays */
  in->GridInfo              = g_malloc0(sizeof(ObitThreadGridInfo));
  in->GridInfo->guardu      = g_malloc0(nGrid*sizeof(ofloat));
  in->GridInfo->guardv      = g_malloc0(nGrid*sizeof(ofloat));
  in->GridInfo->maxBL       = g_malloc0(nGrid*sizeof(ofloat));
  in->GridInfo->minBL       = g_malloc0(nGrid*sizeof(ofloat));
  in->GridInfo->nx          = g_malloc0(nGrid*sizeof(olong));
  in->GridInfo->ny          = g_malloc0(nGrid*sizeof(olong));
  in->GridInfo->uscale      = g_malloc0(nGrid*sizeof(ofloat));
  in->GridInfo->vscale      = g_malloc0(nGrid*sizeof(ofloat));
  in->GridInfo->BeamTaperUV = g_malloc0(nGrid*sizeof(ofloat));
  in->GridInfo->shift       = g_malloc0(nGrid*3*sizeof(ofloat));
  in->GridInfo->rotUV       = g_malloc0(nGrid*9*sizeof(ofloat));
  in->GridInfo->isBeam      = g_malloc0(nGrid*sizeof(gboolean));
  in->GridInfo->thArgs      = g_malloc0(nGrid*sizeof(gpointer));
  /* Fill in global values */
  in->nGrid                  = nGrid;
  in->nGpI                   = nGpI;
  in->GridInfo->nfacet       = nGrid;
  in->GridInfo->nThreads     = nThreads;
  in->GridInfo->nchan        = nChan;
  in->GridInfo->nrparm       = uvDesc->nrparm;
  in->GridInfo->nvis         = nvis;
  in->GridInfo->lenvis       = uvDesc->lrec;
  in->GridInfo->freqArr      = uvDesc->fscale;
  in->GridInfo->convWidth    = UVGrids[0]->convWidth;
  in->GridInfo->convNperCell = UVGrids[0]->convNperCell;
  in->GridInfo->convfn       = UVGrids[0]->convgfn->array;
  in->GridInfo->thread       = ObitThreadRef(UVGrids[0]->thread);

  /* Initialize Function args */
  funcarg    = (GridFuncArg**)in->GridInfo->thArgs;
  wrksize    = nChan; /* Number of channels */
  for (iTh=0; iTh<nGrid; iTh++) {
    funcarg[iTh] = g_malloc(sizeof(GridFuncArg));
    funcarg[iTh]->ithread  = iTh;
    funcarg[iTh]->thread   = ObitThreadRef(UVGrids[0]->thread);
    funcarg[iTh]->vrow     = 0;
    funcarg[iTh]->facet    = iTh;
    funcarg[iTh]->data     = UVin->buffer;
    funcarg[iTh]->gridInfo = in->GridInfo;
    funcarg[iTh]->sigma1   = NULL;
    funcarg[iTh]->sigma2   = NULL;
    funcarg[iTh]->sigma3   = NULL;
    funcarg[iTh]->fwork1   = g_malloc0(2*wrksize*sizeof(ofloat));
    funcarg[iTh]->fwork2   = g_malloc0(2*wrksize*sizeof(ofloat));
    funcarg[iTh]->fwork3   = g_malloc0(wrksize*sizeof(ofloat));
    funcarg[iTh]->cnvfnu   = g_malloc0(wrksize*sizeof(ofloat*));
    funcarg[iTh]->cnvfnv   = g_malloc0(wrksize*sizeof(ofloat*));
    funcarg[iTh]->iuarr    = g_malloc0(wrksize*sizeof(olong));
    funcarg[iTh]->ivarr    = g_malloc0(wrksize*sizeof(olong));
    funcarg[iTh]->beamOrd  = 0;
   } /* End loop creating Thread args */

  /* order grids first by beam/non beam, then in order of UVGrids */
  iGrid = 0;
  for (i=0; i<nPar; i++) {
    if (UVGrids[i]->doBeam) {
      /* add gridding facets */
      for (j=0; j<nGpI; j++) {
	in->GridInfo->isBeam[iGrid]       = TRUE;         /* Is a beam */
	in->GridInfo->BeamTaperUV[iGrid]  = UVGrids[i]->BeamTaperUV;
	in->GridInfo->uscale[iGrid]       = UVGrids[i]->UScale;
	in->GridInfo->vscale[iGrid]       = UVGrids[i]->VScale;
	in->GridInfo->shift[iGrid*3+0]    = UVGrids[i]->dxc;
	in->GridInfo->shift[iGrid*3+1]    = UVGrids[i]->dyc;
	in->GridInfo->shift[iGrid*3+2]    = UVGrids[i]->dzc;
	in->GridInfo->nx[iGrid]           = UVGrids[i]->nxBeam;
	in->GridInfo->ny[iGrid]           = UVGrids[i]->nyBeam;
	in->GridInfo->guardu[iGrid]       = UVGrids[i]->guardband*UVGrids[i]->nxBeam;
	in->GridInfo->guardv[iGrid]       = UVGrids[i]->guardband*UVGrids[i]->nyBeam;
	in->GridInfo->maxBL[iGrid]        = UVGrids[i]->blmax;
	in->GridInfo->minBL[iGrid]        = UVGrids[i]->blmin;
	/* u,v,w rotation if do3Dmul */
	if (UVGrids[i]->do3Dmul) {
	  for (k=0; k<3; k++) in->GridInfo->rotUV[iGrid*9+k  ] = UVGrids[i]->URot3D[0][k];
	  for (k=0; k<3; k++) in->GridInfo->rotUV[iGrid*9+3+k] = UVGrids[i]->URot3D[1][k];
	  for (k=0; k<3; k++) in->GridInfo->rotUV[iGrid*9+6+k] = UVGrids[i]->URot3D[2][k];
	} else {  /* identity matrix */
	  for (k=0; k<9; k++) in->GridInfo->rotUV[iGrid*9+k] = noRot[k];
	}
	/* Thread arg stuff */
	funcarg[iGrid]->iGpI   = j;       /* grid number for this beam */
	funcarg[iGrid]->nGpI   = nGpI;    /* number of grids for this beam */
	funcarg[iGrid]->facet  = iGrid;	  /* facet (beam) number */
	funcarg[iGrid]->lovis  = j*nvisPth;	
	funcarg[iGrid]->hivis  = (j+1)*nvisPth;	
	size = 2 * (1 + UVGrids[i]->convWidth/2 + UVGrids[i]->nxBeam/2) * 
	  UVGrids[i]->nyBeam;
	funcarg[iGrid]->grid = g_malloc0(size*sizeof(ofloat));
	funcarg[iGrid]->outGrid = ObitCArrayRef(UVGrids[i]->grid);
	funcarg[iGrid]->beamOrd = UVGrids[i]->order;  /* (SW) beam order */
	/* Channel selection - all */
	funcarg[iGrid]->bChan = 0;
	funcarg[iGrid]->eChan = nChan;
	iGrid++;
      } /* end loop over grids per image */
      funcarg[iGrid-1]->hivis = MAX(funcarg[iGrid-1]->hivis, nvis);
    } /* end if beam */
  } /* End loop adding beams */
 
 /* Loop over images */
  for (i=0; i<nPar; i++) {
    if (!UVGrids[i]->doBeam) {
      /* add gridding facets */
      for (j=0; j<nGpI; j++) {
	in->GridInfo->isBeam[iGrid]       = FALSE;         /* Not a beam */
	in->GridInfo->BeamTaperUV[iGrid]  = UVGrids[i]->BeamTaperUV;
	in->GridInfo->uscale[iGrid]       = UVGrids[i]->UScale;
	in->GridInfo->vscale[iGrid]       = UVGrids[i]->VScale;
	in->GridInfo->shift[iGrid*3+0]    = UVGrids[i]->dxc;
	in->GridInfo->shift[iGrid*3+1]    = UVGrids[i]->dyc;
	in->GridInfo->shift[iGrid*3+2]    = UVGrids[i]->dzc;
	in->GridInfo->nx[iGrid]           = UVGrids[i]->nxImage;
	in->GridInfo->ny[iGrid]           = UVGrids[i]->nyImage;
	in->GridInfo->guardu[iGrid]       = UVGrids[i]->guardband*UVGrids[i]->nxImage;
	in->GridInfo->guardv[iGrid]       = UVGrids[i]->guardband*UVGrids[i]->nyImage;
	in->GridInfo->maxBL[iGrid]        = UVGrids[i]->blmax;
	in->GridInfo->minBL[iGrid]        = UVGrids[i]->blmin;
	/* u,v,w rotation if do3Dmul */
	if (UVGrids[i]->do3Dmul) {
	  for (k=0; k<3; k++) in->GridInfo->rotUV[iGrid*9+k  ] = UVGrids[i]->URot3D[0][k];
	  for (k=0; k<3; k++) in->GridInfo->rotUV[iGrid*9+3+k] = UVGrids[i]->URot3D[1][k];
	  for (k=0; k<3; k++) in->GridInfo->rotUV[iGrid*9+6+k] = UVGrids[i]->URot3D[2][k];
	} else {  /* identity matrix */
	  for (k=0; k<9; k++) in->GridInfo->rotUV[iGrid*9+k] = noRot[k];
	}
	/* Thread arg stuff */
	funcarg[iGrid]->iGpI   = j;       /* grid number for this image */
	funcarg[iGrid]->nGpI   = nGpI;    /* number of grids for this image */
	funcarg[iGrid]->facet  = iGrid;	  /* facet (image) number */
	funcarg[iGrid]->lovis  = j*nvisPth;	
	funcarg[iGrid]->hivis  = (j+1)*nvisPth;
	size = 2 * (1 + UVGrids[i]->convWidth/2 + UVGrids[i]->nxImage/2) * 
	  UVGrids[i]->nyImage;
	funcarg[iGrid]->grid = g_malloc0(size*sizeof(ofloat));
	funcarg[iGrid]->outGrid = ObitCArrayRef(UVGrids[i]->grid);
	/* Channel selection */
	funcarg[iGrid]->bChan = 0;
	funcarg[iGrid]->eChan = nChan;
	iGrid++;
      } /* end loop over grids per image */
      funcarg[iGrid-1]->hivis = MAX(funcarg[iGrid-1]->hivis, nvis);
    } /* end if image */
  } /* End loop adding image */
 
}  /* end ObitThreadGridSetupWB */

/**
 * Grid UV data with threads per facet
 * call ObitThreadPoolFree (thread) after last call.
 * \param grids      Thread grid object
 */
void ObitThreadGridGrid (ObitThreadGrid *grids)
{
  olong iLoop, nLeft, nvis, nvisPth, nDo, j;
  olong nThreads        = grids->GridInfo->nThreads;
  olong nGrid           = grids->nGrid;
  ObitThread *thread    = grids->GridInfo->thread;
  GridFuncArg **funcarg = (GridFuncArg**)grids->GridInfo->thArgs;
  gboolean OK;

  /* Anything to do? */
  if (grids->UVin->myDesc->numVisBuff<=0) return;

  /* Did the number of visibilities change? */
  if (grids->GridInfo->nvis!=grids->UVin->myDesc->numVisBuff) {
    /* Reset on funcarg */
    nvis                  = grids->UVin->myDesc->numVisBuff;
    nvisPth               = nvis/grids->nGpI;   /* Vis per thread */
    grids->GridInfo->nvis = nvis;
    for (iLoop=0; iLoop<grids->nGrid; iLoop++) {
      j = funcarg[iLoop]->iGpI;
      funcarg[iLoop]->lovis  = j*nvisPth;	
      funcarg[iLoop]->hivis  = (j+1)*nvisPth;
      if (j==(funcarg[iLoop]->nGpI-1)) funcarg[iLoop]->hivis = nvis;
    }
  }

  nLeft = nGrid;
  for (iLoop=0; iLoop<nGrid; iLoop+=nThreads) {
    nDo = MIN(nThreads, nLeft);
    if (nDo==1) funcarg[iLoop]->ithread = -1;  /* Only one? */
    OK = ObitThreadIterator (thread, nDo, ThreadGrid, 
			     (gpointer **)&funcarg[iLoop]);
    if ((nDo==1) && (nThreads>1)) funcarg[iLoop]->ithread = iLoop;  /* reset if needed */
    nLeft -= nThreads;
  }

  /* Check for problems
     if (!OK) Obit_log_error(err, OBIT_Error,"%s: Problem in threading", routine);*/
} /* end ObitThreadGridGrid */

/**
 * Flip UV data with threads per facet
 * Add cells in negative u half plane to conjugate cell
 * \param grids      Thread grid object
 */
void  ObitThreadGridFlip (ObitThreadGrid *grids)
{
  olong iLoop, nLeft, nDo;
  olong nThreads        = grids->GridInfo->nThreads;
  olong nGrid           = grids->nGrid;
  ObitThread *thread    = grids->GridInfo->thread;
  GridFuncArg **funcarg = (GridFuncArg**)grids->GridInfo->thArgs;
  gboolean OK;

  nLeft = nGrid;
  for (iLoop=0; iLoop<nGrid; iLoop+=nThreads) {
    nDo = MIN(nThreads, nLeft);
    if (nDo==1) funcarg[iLoop]->ithread = -1;  /* Only one? */
    OK = ObitThreadIterator (thread, nDo, ThreadFlip, 
			     (gpointer **)&funcarg[iLoop]);
    if ((nDo==1) && (nThreads>1)) funcarg[iLoop]->ithread = iLoop;  /* reset if needed */
    nLeft -= nThreads;
  }

  ObitThreadPoolFree (thread);  /* Free thread pool */

  /* Check for problems
     if (!OK) Obit_log_error(err, OBIT_Error,"%s: Problem in threading", routine);*/
  /* Really should clean up */
} /* end ObitThreadGridFlip */


/**
 * Merge/transpose one or more grids for a given image/beam
 * leaving in center-at-edges for FFT
 * The order of the threading objects should be grids for a given
 * image/beam are contigious and all have grids->nGpI.
 * \param grids      Thread grid object
*/
void ObitThreadGridMerge (ObitThreadGrid *grids)
{
  olong iLoop, jLoop, ia, nLeft, nDo, nIt, nTh;
  olong nThreads        = grids->GridInfo->nThreads;
  olong nGrid           = grids->nGrid;
  ObitThread *thread    = grids->GridInfo->thread;
  GridFuncArg **funcarg = (GridFuncArg**)grids->GridInfo->thArgs;
  GridFuncArg **tmparg  = NULL;
  gboolean OK;

  /* Temporary array for arguments */
  nDo = grids->nGrid/grids->nGpI;  /* how many actual images/beams */
  tmparg = g_malloc0(nDo*sizeof(GridFuncArg*));

  /* Loop over grids per image, accumulating one at a time */
  for (jLoop=0; jLoop<grids->nGpI; jLoop++) {
    /* Get array of arguments */
    ia = 0; nDo = 0;
    for (iLoop=0; iLoop<nGrid; iLoop++) {
      /* Include if iGpI = jLoop */
      if (funcarg[iLoop]->iGpI==jLoop) {nDo++; tmparg[ia++] = funcarg[iLoop];}
    }
    /* No more than actual number of image/beam at a time */
    nTh = MIN (nThreads, grids->nGrid/grids->nGpI);
    /* Loop over this set */
    nLeft = nDo;
    for (iLoop=0; iLoop<nDo; iLoop+=nTh) {
     nIt = MIN(nTh, nLeft);
     if (nIt==1) tmparg[iLoop]->ithread = -1;  /* Only one? */
     OK = ObitThreadIterator (thread, nIt, ThreadMerge, 
			       (gpointer **)&tmparg[iLoop]);
     if ((nDo==1) && (nThreads>1)) tmparg[iLoop]->ithread = iLoop;  /* reset if needed */
      nLeft -= nIt;
    }
  } /* end jLoop */

  ObitThreadPoolFree (thread);  /* Free thread pool */

  /* Cleanup */
  if (tmparg) g_free(tmparg);

  /* Check for problems
     if (!OK) Obit_log_error(err, OBIT_Error,"%s: Problem in threading", routine);*/
} /* end ObitThreadGridMerge */


/**
 * Grid UV data with threads per facet
 * Add cells in negative u half plane to conjugate cell
 * \param args  Threading argument
 */
static gpointer ThreadGrid (gpointer args)
{
  GridFuncArg *largs    = (GridFuncArg*)args;
  ObitThreadGridInfo *gridInfo = largs->gridInfo;
  olong ifacet          = largs->facet;
  ofloat *grid          = largs->grid;
  olong  bChan          = largs->bChan;
  olong  eChan          = largs->eChan;

  olong halfWidth    = gridInfo->convWidth/2;         // half width of convolution kernal
  olong fullWidth    = gridInfo->convWidth;           // full width of convolution kernal
  olong kvis, ichan, lrow;

  lrow  = 2*(1 + gridInfo->nx[ifacet]/2 + halfWidth);  // length of grid row in floats
  eChan = MAX (eChan, bChan+1);  /* At least 1 channel */
  /* Loop over vis */
  for (kvis=largs->lovis; kvis<largs->hivis; kvis++) {
    /* Prep grid */
    fast_prep_grid (kvis, largs);
    /* Grid channels */
    for (ichan=bChan; ichan<eChan; ichan++) {
      fast_grid(grid, &largs->fwork1[ichan*2], largs->iuarr[ichan], largs->ivarr[ichan], 
		lrow, fullWidth, largs->cnvfnu[ichan], largs->cnvfnv[ichan]);
    }  /* end channel loop */
  } /* end vis loop */
  /* Done
 done: */
  if (largs->ithread>=0)
    ObitThreadPoolDone (largs->thread, (gpointer)&largs->ithread);
  
  return NULL;
} /* end ThreadGrid */

/**
 * Flip/add conjugate rows in grid in one thread
 * Add cells in negative u half plane to conjugate cell
 * \param args  Threading argument
 */
static gpointer ThreadFlip (gpointer args)
{
  GridFuncArg *largs = (GridFuncArg*)args;
  ObitThreadGridInfo *gridInfo = largs->gridInfo;
  olong ifacet          = largs->facet;
  ofloat *grid          = largs->grid;
  olong halfWidth       = gridInfo->convWidth/2;    /* half width of convolution kernal */
  ofloat *gxi, *gxo, *gci, *gco, xxo[2], cjo[2], xxi[2], cji[2];
  olong ny, ny2, vrow, iu, vc, lrow;

  lrow = 2*(1 + gridInfo->nx[ifacet]/2 + halfWidth);  /* length of grid row */
  ny   = gridInfo->ny[ifacet];
  ny2  = ny/2;
  if (ny<0) goto done;
  /* loop over rows, center v = ny/2 (0 rel) */
  for (vrow=1; vrow<=ny2; vrow++) {
    vc = ny - vrow;               /* conjugate row number */
    /* loop over u columns */
    gci = grid + vc*lrow   + (halfWidth)*2;
    gxo = grid + vrow*lrow + (halfWidth)*2;
    gxi = gxo; gco = gci;
    for (iu=0; iu<=halfWidth; iu++) {
      /* Read initial values from both sides of both rows */
      cji[0] = gci[0]; cji[1] = gci[1];
      cjo[0] = gco[0]; cjo[1] = gco[1];
      xxi[0] = gxi[0]; xxi[1] = gxi[1];
      xxo[0] = gxo[0]; xxo[1] = gxo[1];
      /* update both row and conjugate row */
      gxo[0] = xxo[0] + cji[0];
      gxo[1] = xxo[1] - cji[1];
      gco[0] = cjo[0] + xxi[0];
      gco[1] = cjo[1] - xxi[1];
	/* hard core debug 
      if ((fabs(gxo[1]>0.00001) || (fabs(gxo[1]>0.00001)) {
	  fprintf (stdout,"dbg %d %15.7f %15.7f  \n",
		   iu, gxo[1], gco[1]);
	}*/
      gxo += 2; gco += 2; gxi -= 2; gci -= 2;
    } /* end u loop */
  } /* end v loop */
 
 /* Done */
 done:
  if (largs->ithread>=0)
    ObitThreadPoolDone (largs->thread, (gpointer)&largs->ithread);
  
  return NULL;
} /* end ThreadFlip */

/**
 * shuffle/merge grids for an image/beam in one thread
 *  output zeroed on first
 * \param args  Threading argument
 */
static gpointer ThreadMerge (gpointer args)
{
  GridFuncArg *largs = (GridFuncArg*)args;
  ObitThreadGridInfo *gridInfo = largs->gridInfo;
  olong ifacet        = largs->facet;
  olong halfWidth     = gridInfo->convWidth/2;    /* half width of convolution kernal */
  ofloat *grid        = largs->grid;
  ofloat *outGrid     = largs->outGrid->array;
  ofloat *gi, *go, czero[] = {0.0,0.0};
  olong ny, vrow, iu, vs, ilrow, olrow, ncopy;

  /* Zero output on first */
  if (largs->iGpI==0) ObitCArrayFill(largs->outGrid, czero);

  ilrow = 2*(1 + gridInfo->nx[ifacet]/2 + halfWidth);  /* length of input grid row */
  olrow = 2*(1 + gridInfo->nx[ifacet]/2);              /* length of output grid row */
  ncopy = olrow*sizeof(ofloat);
  ny   = gridInfo->ny[ifacet];
  if (ny<0) goto done;

  /* loop over half v rows */
  vs = (ny/2);         /* swap row number */
  for (vrow=0; vrow<ny/2; vrow++) {
    /* loop over u columns */
    gi = grid + vrow*ilrow + halfWidth*2;
    go = outGrid + vs*olrow;
    for (iu=0; iu<olrow; iu+=2) {
      go[0] += gi[0];
      go[1] += gi[1];
      go += 2; gi += 2;
    }  /*end u loop */
    vs++;
  } /* end v loop */

  /* Loop over other half */
  vs = 0;
  for (vrow=ny/2; vrow<ny; vrow++) {
    /* loop over u columns */
    gi = grid    + vrow*ilrow + halfWidth*2;
    go = outGrid + vs*olrow;
    for (iu=0; iu<olrow; iu+=2) {
      go[0] += gi[0];
      go[1] += gi[1];
      go += 2; gi += 2;
    }  /* end u loop */
    vs++;
  } /* end v loop */

 /* Done */
 done:
  if (largs->ithread>=0)
    ObitThreadPoolDone (largs->thread, (gpointer)&largs->ithread);
  
  return NULL;
} /* end ThreadMerge */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitThreadGridClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitThreadGridClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitThreadGridClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitThreadGridClassInfoDefFn (gpointer inClass)
{
  ObitThreadGridClassInfo *theClass = (ObitThreadGridClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitThreadGridClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitThreadGridClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitThreadGridGetClass;
  theClass->newObit       = (newObitFP)newObitThreadGrid;
  theClass->ObitCopy      = NULL;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitThreadGridClear;
  theClass->ObitInit      = (ObitInitFP)ObitThreadGridInit;
  theClass->ObitThreadGridSetup = (ObitThreadGridSetupFP)ObitThreadGridSetup;
  theClass->ObitThreadGridGrid  = (ObitThreadGridGridFP)ObitThreadGridGrid;
  theClass->ObitThreadGridFlip  = (ObitThreadGridFlipFP)ObitThreadGridFlip;
  theClass->ObitThreadGridMerge = (ObitThreadGridMergeFP)ObitThreadGridMerge;

} /* end ObitThreadGridClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitThreadGridInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitThreadGrid *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->thread       = NULL;
  in->info         = newObitInfoList(); 
  in->myStatus     = OBIT_Inactive;
  in->GridInfo     = NULL;
  in->UVin         = NULL;
} /* end ObitThreadGridInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * For some reason this wasn't build into the GType class.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitThreadGrid* cast to an Obit*.
 */
void ObitThreadGridClear (gpointer inn)
{
  olong i;
  ObitClassInfo *ParentClass;
  GridFuncArg **funcarg;
  ObitThreadGrid *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  in->thread    = ObitThreadUnref(in->thread);
  in->info      = ObitInfoListUnref(in->info);
  in->UVin      = ObitUVUnref(in->UVin);
  if (in->GridInfo) {
    if (in->GridInfo->guardu) g_free(in->GridInfo->guardu);
    if (in->GridInfo->guardv) g_free(in->GridInfo->guardv);
    if (in->GridInfo->maxBL)  g_free(in->GridInfo->maxBL);
    if (in->GridInfo->minBL)  g_free(in->GridInfo->minBL);
    if (in->GridInfo->nx)     g_free(in->GridInfo->nx);
    if (in->GridInfo->ny)     g_free(in->GridInfo->ny);
    if (in->GridInfo->uscale) g_free(in->GridInfo->uscale);
    if (in->GridInfo->vscale) g_free(in->GridInfo->vscale);
    if (in->GridInfo->BeamTaperUV) g_free(in->GridInfo->BeamTaperUV);
    if (in->GridInfo->shift)  g_free(in->GridInfo->shift);
    if (in->GridInfo->rotUV)  g_free(in->GridInfo->rotUV);
    if (in->GridInfo->isBeam) g_free(in->GridInfo->isBeam);
    funcarg = (GridFuncArg**)in->GridInfo->thArgs;
    if (funcarg) {
      for (i=0; i<in->nGrid; i++) {
	if (funcarg[i]) {
	  funcarg[i]->gridInfo = ObitThreadGridUnref(funcarg[i]->gridInfo);
	  funcarg[i]->thread   = ObitThreadUnref(funcarg[i]->thread);
	  if (funcarg[i]->fwork1) g_free(funcarg[i]->fwork1);
	  if (funcarg[i]->fwork2) g_free(funcarg[i]->fwork2);
	  if (funcarg[i]->fwork3) g_free(funcarg[i]->fwork3);
	  if (funcarg[i]->cnvfnu) g_free(funcarg[i]->cnvfnu);
	  if (funcarg[i]->cnvfnv) g_free(funcarg[i]->cnvfnv);
	  if (funcarg[i]->iuarr)  g_free(funcarg[i]->iuarr);
	  if (funcarg[i]->ivarr)  g_free(funcarg[i]->ivarr);
	  if (funcarg[i]->grid)   g_free(funcarg[i]->grid);
	  ObitCArrayUnref(funcarg[i]->outGrid);
	}
	/*g_free(funcarg);??*/
      }
    }
    g_free(in->GridInfo);
  } /* end if in->GridInfo */
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitThreadGridClear */

/** 
 * c version first pass at preparing to grid one visibility/facet
 * \param kvis   visibiliy number (0-rel)
 * \param args   Threaded gridding function argument
 * Saves values:
 * \li   fwork1    Visibility array as (r,i)
 * \li   cnvfnu    Address of u convolving vector
 * \li   cnvfnv    Address of v convolving vector
 * \li   iuarr     first u cell (0-rel)
 * \li   ivarr     first v cell (0-rel)
*/
void fast_prep_grid(olong kvis, GridFuncArg *args)  
{
  ObitThreadGridInfo *gridInfo = args->gridInfo;
  olong ifacet          = args->facet;
  ofloat *vis_in        = args->data;
  olong  bChan          = args->bChan;
  olong  eChan          = args->eChan;
  olong  beamOrd        = args->beamOrd;

  olong halfWidth    = gridInfo->convWidth/2;         // half width of convolution kernal
  olong fullWidth    = gridInfo->convWidth;           // full width of convolution kernal
  olong convNperCell = gridInfo->convNperCell;        // resolution of kernal
  olong ichan, ivis, jvis, lrow, halfv, it, iphase, iu, iv;
  ofloat *rot        = &gridInfo->rotUV[ifacet*9];
  ofloat *shift      = &gridInfo->shift[ifacet*3];
  ofloat *convfn     = gridInfo->convfn;
  ofloat u,v,w, uu, vv, ww, maxBL2, minBL2, bmTaper, BL2, fact;
  ofloat vr, vi, vw, vvr, vvi, phase, guardu, guardv, ftemp;
  ofloat c, s, phaseSign, freqFact;
  gdouble dshift[3], dphase, doTape;
  gboolean want;
  static const gdouble twopi = 2*G_PI;
  static const gdouble itwopi = 1.0/(2*G_PI);

  lrow  = 2*(1 + gridInfo->nx[ifacet]/2 + halfWidth);  // length of grid row in floats
  eChan = MAX (eChan, bChan+1);  /* At least 1 channel */
  halfv = gridInfo->ny[ifacet]/2;
  maxBL2  = gridInfo->maxBL[ifacet]*gridInfo->maxBL[ifacet];
  minBL2  = gridInfo->minBL[ifacet]*gridInfo->minBL[ifacet];
  bmTaper = gridInfo->BeamTaperUV[ifacet];
  doTape  = (fabs(bmTaper)>0.0) || ((args->sigma1!=NULL) && (args->sigma1[0]!=0.0));

  /* guardband in wavelengths */
  guardu = fabs(gridInfo->guardu[ifacet] / gridInfo->uscale[ifacet]);
  guardv = fabs(gridInfo->guardv[ifacet] / gridInfo->vscale[ifacet]);

  dshift[0] = (gdouble)shift[0];  dshift[1] = (gdouble)shift[1];  dshift[2] = (gdouble)shift[2];
  ivis = kvis * gridInfo->lenvis; /* beginning of visibility */
  /*  Assume random parameters start with u,v,w */
  u = vis_in[ivis];
  v = vis_in[ivis+1];
  w = vis_in[ivis+2];
  /* rotate u,v,w for facet */
  uu = u*rot[0] + v*rot[1] + w*rot[2];
  vv = u*rot[3] + v*rot[4] + w*rot[5];
  ww = u*rot[6] + v*rot[7] + w*rot[8];
  /* Only gridding half plane, need to flip to other side? */
  if (uu<=0.0) {
    phaseSign = -1.0;
  } else { /* no flip */
    phaseSign = 1.0;
  }
  /* loop over channels for position shift */
  for (ichan=bChan; ichan<eChan; ichan++) {
    freqFact = phaseSign * gridInfo->freqArr[ichan];
    u = uu * freqFact;  // Scale u,v,w to channel
    v = vv * freqFact;
    w = ww * freqFact;
    jvis = ivis + gridInfo->nrparm + ichan*3;
    vvr = vis_in[jvis];
    vvi = phaseSign*vis_in[jvis+1];   /* Conjugate if neg u */
    vw  = vis_in[jvis+2];
    /* Data valid? positive weight and within guardband */
    want = ((vw>0.) && (u<guardu) && (fabs(v)<guardv));
    /* Baseline limits */
    BL2 = u*u + v*v;
    if ((maxBL2>0) && want) want = want && maxBL2>BL2;
    if ((minBL2>0) && want) want = want && minBL2<BL2;
    if (want) {
      /* If this a beam - don't bother shifting - replace data with (wt*1,0) */
      if (!gridInfo->isBeam[ifacet]) {
	/* real part of vis */
	/* position shift in double, reduce range */
	dphase =  (u*dshift[0] + v*dshift[1] + w*dshift[2]);
	iphase = (olong)(dphase*itwopi);
	phase  = (ofloat)(dphase - iphase * twopi);
	ObitSinCosCalc (phase, &s, &c);
	vr = c*vvr - s*vvi;
	vi = s*vvr + c*vvi; 
      } else {
	vr = 1.0; vi = 0.0;
      }
      /* Tapering? */
      if (doTape) {
	/* Beam taper? */
	/* Other (MF) tapering */
	if ((args->sigma1!=NULL) && (args->sigma1[ichan]!=0.0)) {
	  ftemp = u*u*(args->sigma2[ichan]+bmTaper) + 
	          v*v*(args->sigma1[ichan]+bmTaper) + 
	          u*v*(args->sigma3[ichan]);
	  fact = ObitExpCalc (ftemp);
	  vw  *= fact;
	}  /* end Other (MF) tapering */
	else if (fabs(bmTaper)>0.0) {
	  fact = ObitExpCalc(bmTaper*BL2);
	  vw  *= fact;
	}  /* end beam Taper */
      } /* end any taper */
      /* SW Beam order */
      switch (beamOrd) {
      case 0:   /* Dirty beam */
	break;
      case 1:   /* first order ln(nu/nu0) */
	ftemp = log(gridInfo->freqArr[ichan]);
	vw *= ftemp;
	break;
      case 2:   /* second order ln(nu/nu0)**2 */
	ftemp = log(gridInfo->freqArr[ichan]);
	vw *= ftemp*ftemp;
	break;
      case 3:   /* third order ln(nu/nu0)**3 */
	ftemp = log(gridInfo->freqArr[ichan]);
	vw *= ftemp*ftemp*ftemp;
	break;
      default:
	break;
      }; /* end beamOrd switch */
      /* weighted data - save in fwork1 */
      args->fwork1[2*ichan]   = vr * vw;
      args->fwork1[2*ichan+1] = vi * vw;
      /* convert to u,v cells */
      u *= gridInfo->uscale[ifacet];
      v *= gridInfo->vscale[ifacet];
      iu = (olong)(u+0.5);   /* to grid cell number */
      iv = _lroundf(v);
      /* start in convolution function in blocks of fullWidth */
      it = convNperCell + _lroundf(convNperCell * (iu - u - 0.5));
      args->cnvfnu[ichan] = convfn + fullWidth * it;
      it = convNperCell + _lroundf(convNperCell * (iv - v - 0.5));
      args->cnvfnv[ichan] = convfn + fullWidth * it;
      args->iuarr[ichan]  = iu;
      args->ivarr[ichan]  = iv - halfWidth + halfv;
    } /* end if valid */
    else {  /* Invalid data */
      args->fwork1[2*ichan]   = 0.0;
      args->fwork1[2*ichan+1] = 0.0;
      args->iuarr[ichan]      = halfWidth; 
      args->ivarr[ichan]      = halfv;
      args->cnvfnu[ichan]     = convfn;
      args->cnvfnv[ichan]     = convfn ;
    } /* end invalid data */
  } /* end channel loop */
} /* end fast_prep_grid */

/** AVX implementation 8 floats in parallel */
#if HAVE_AVX==1
#include <immintrin.h>
/* gcc or icc */
# define ALIGN32_BEG
# define ALIGN32_END __attribute__((aligned(32)))

typedef __m256  v8sf;
typedef __m256d v4df;
typedef __m128  v4sf;
typedef __m256i v8si;
typedef __m128i v4si;
/* Unions allowing c interface */
typedef ALIGN32_BEG union {
  float f[8];
  int   i[8];
  v8sf  v;
} ALIGN32_END V8SF;
typedef ALIGN32_BEG union {
  float f[4];
  int   i[4];
  v4sf  v;
} ALIGN32_END V4SF;
typedef ALIGN32_BEG union {
  double    f[4];
  long long i[4];
  v4df      v;
} ALIGN32_END V4DF;

static const v4si _MASK3 = {0xffffffffffffffff, 0x00000000ffffffff};
static const v4si _MASK4 = {0xffffffffffffffff, 0xffffffffffffffff};
static const v8si _mask3 = {0xffffffffffffffff, 0x00000000ffffffff, 0x0000000000000000, 0x0000000000000000};
static const v8si _mask6 = {0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0x0000000000000000};
static const v8si _mask7 = {0xffffffffffffffff, 0xffffffffffffffff, 0xffffffffffffffff, 0x00000000ffffffff};
static const v8sf _half  = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5}; /* 0.5 vector */
static const v8sf _mhalf = {-0.5, -0.5,- 0.5, -0.5, -0.5, -0.5, -0.5, -0.5}; /* -0.5 vector */
static const v8sf _one   =  {1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0, 1.0}; /* 1 vector */

/** 
 * SOMETHING IS WRONG IN THIS VERSION - using c version
 * Get aliased version of strong source at origin
 * AVX version first pass at preparing to grid one visibility/facet
 * \param kvis   visibiliy number (0-rel)
 * \param args   Threaded gridding function argument
 * Saves values:
 * \li   fwork1    Visibility array as (r,i)
 * \li   cnvfnu    Address of u convolving vector
 * \li   cnvfnv    Address of v convolving vector
 * \li   iuarr     first u cell (0-rel)
 * \li   ivarr     first v cell (0-rel)
*/
void fast_prep_gridAVX(olong kvis, GridFuncArg *args)  
{
  ObitThreadGridInfo *gridInfo = args->gridInfo;
  olong ifacet          = args->facet;
  ofloat *vis_in        = args->data;
  olong  bChan          = args->bChan;
  olong  eChan          = args->eChan;
  olong  beamOrd        = args->beamOrd;

  olong halfWidth    = gridInfo->convWidth/2;         // half width of convolution kernal
  olong fullWidth    = gridInfo->convWidth;           // full width of convolution kernal
  olong convNperCell = gridInfo->convNperCell;        // resolution of kernal
  ofloat *rot        = &gridInfo->rotUV[ifacet*9];
  ofloat *shift      = &gridInfo->shift[ifacet*3];
  ofloat *convfn     = gridInfo->convfn;
  ofloat *saveu      = args->fwork1;   /* aliases for temp storage */
  ofloat *savev      = args->fwork1+gridInfo->nchan;
  ofloat *savew      = args->fwork2;
  ofloat *saveph     = args->fwork1;   /* aliases for temp storage */
  ofloat *saves      = args->fwork2;
  ofloat *savec      = args->fwork3;
  ofloat *valid      = args->fwork2+gridInfo->nchan;
  olong ichan, jchan, kchan, ivis, jvis, lrow, halfv, it, iphase;
  olong ndone, iu, iv;
  ofloat u, v, w, uu, vv, ww, iscaleu, iscalev, ftemp, farr[8];
  ofloat vr, vi, tr, ti, tw, maxBL2, minBL2, bmTaper, BL2, fact;
  ofloat c, s, phaseSign, freqFact;
  gboolean doTape;
  gdouble dshift[3], dphase;
  static const ofloat  twopi = 2*G_PI;
  static const gdouble itwopi = 1.0/(2*G_PI);

  v8sf vuvw, vrot1, vrot2, vt, vt1, vt2, vt3, vt4, vru, vrv, vrw, vf, vvr, vvi, vw, vpr, vpi;
  v8sf vtr, vti, vvalid, vguardu, vguardv, vmguardv, vBL2;
  v4df du, dv, dw, ds1, ds2, ds3, dt1, dt2, dt3;
  V8SF vout, vin, chu, chv, chw;
  v4sf v4t, v4out ;

  lrow  = 2*(1 + gridInfo->nx[ifacet]/2 + halfWidth);  // length of grid row in floats
  eChan = MAX (eChan, bChan+1);  /* At least 1 channel */
  halfv = gridInfo->ny[ifacet]/2;
  iscaleu = 1.0 / gridInfo->uscale[ifacet];
  iscalev = 1.0 / gridInfo->vscale[ifacet];
  maxBL2  = gridInfo->maxBL[ifacet]*gridInfo->maxBL[ifacet];
  minBL2  = gridInfo->minBL[ifacet]*gridInfo->minBL[ifacet];
  bmTaper = gridInfo->BeamTaperUV[ifacet];
  doTape  = (fabs(bmTaper)>0.0) || ((args->sigma1!=NULL) && (args->sigma1[0]!=0.0));

  ivis = kvis * gridInfo->lenvis; /* beginning of visibility */
  /*  Assume random parameters start with u,v,w */
  vuvw  =  _mm256_maskload_ps(&vis_in[ivis], _mask3); /* Only 3 */
  v4t   =  _mm_maskload_ps (&vis_in[ivis], _MASK3);
  vuvw  =  _mm256_insertf128_ps (vuvw, v4t, 1);   /* Copy in second 4 */
  vrot1 =  _mm256_maskload_ps(&rot[0], _mask3);
  v4t   =  _mm_maskload_ps (&rot[3],  _MASK3);
  vrot1 =  _mm256_insertf128_ps (vrot1, v4t, 1);   /* 2nd row in second 4 */
  vrot2 =  _mm256_maskload_ps(&rot[6], _mask3);
  vout.v=  _mm256_dp_ps(vuvw, vrot1, 0xff);
  uu    = vout.f[0];
  vv    = vout.f[4];
  vout.v=  _mm256_dp_ps(vuvw, vrot2, 0xff);
  ww    = vout.f[0];
  /* Only gridding half plane, need to flip to other side? */
  if (uu<=0.0) {
    phaseSign = -1.0;
  } else { /* no flip */
    phaseSign = 1.0;
  }
  vru = _mm256_set1_ps (uu*gridInfo->uscale[ifacet]);  /* Promote ref u,v,w to vector in cells */
  vrv = _mm256_set1_ps (vv*gridInfo->vscale[ifacet]);
  vrw = _mm256_set1_ps (ww*gridInfo->uscale[ifacet]);  /* use u scale for w */

  /* loop over channels, scale by freq, save, u,v,w get center, convfn*/
  vguardu  = _mm256_set1_ps ( gridInfo->guardu[ifacet]);
  vguardv  = _mm256_set1_ps ( gridInfo->guardv[ifacet]);
  vmguardv = _mm256_set1_ps (-gridInfo->guardv[ifacet]);
  ndone = bChan;
  for (ichan=bChan; ichan<eChan; ichan+=8) {
    if ((ichan+8)>=eChan) break;  /* must have all 8 */
    ndone += 8;  /* how many channels done */
    vt = _mm256_set1_ps (phaseSign);
    vf = _mm256_loadu_ps(&gridInfo->freqArr[ichan]);
    vt = _mm256_mul_ps(vt, vf);
    chu.v = _mm256_mul_ps(vru, vt);  /* Scale u,v,w to cells at channel */
    _mm256_storeu_ps(&saveu[ichan], chu.v); /* save u in cells at channel */
    chv.v = _mm256_mul_ps(vrv, vt);
    _mm256_storeu_ps(&savev[ichan], chv.v); 
    chw.v = _mm256_mul_ps(vrw, vt);
    _mm256_storeu_ps(&savew[ichan], chw.v); 
    /* Check guardband */
    vtr = _mm256_cmp_ps (chu.v, vguardu, _CMP_LT_OQ);
    vti = _mm256_cmp_ps (chv.v, vguardv, _CMP_LT_OQ);
    vvalid = _mm256_and_ps (vtr, vti);
    vti = _mm256_cmp_ps (chv.v, vmguardv, _CMP_GT_OQ);
    vvalid = _mm256_and_ps (vvalid, vti);
    vtr = _mm256_setzero_ps ();
    /* Check baseline length - include in vvalid */
    vt   = _mm256_set1_ps (iscaleu);
    vt1  = _mm256_mul_ps(chu.v, vt);   /* u in wavelengths */
    vt   = _mm256_set1_ps (iscalev);
    vt2  = _mm256_mul_ps(chv.v, vt);   /* v in to wavelengths */
    vt3  = _mm256_mul_ps(vt1, vt2);    /* u*v  */
    vt1  = _mm256_mul_ps(vt1, vt1);    /* u^2 */
    vt2  = _mm256_mul_ps(vt2, vt2);    /* v^2 */
    vBL2 = _mm256_add_ps(vt1, vt2);    /* Baseline**2 */
    if (maxBL2>0.0) {   /* Check max */
      vt  = _mm256_set1_ps (maxBL2);
      vt4 = _mm256_cmp_ps (vBL2, vt, _CMP_LT_OQ);
      vvalid = _mm256_and_ps (vvalid, vt4);
    }
    if (minBL2>0.0) {
      vt  = _mm256_set1_ps (minBL2);
      vt4 = _mm256_cmp_ps (vBL2, vt, _CMP_GT_OQ);
      vvalid = _mm256_and_ps (vvalid, vt4);
    }/* end max/min BL */
    vvalid = _mm256_blendv_ps (vtr, _one, vvalid);  /* 1/0 mask for valid */
    /* Tapering? */
    if (doTape) {
      /* Other (MF) tapering + bmTaper*/
      if ((args->sigma1!=NULL) && (args->sigma1[ichan]!=0.0)) {
	vf    = _mm256_set1_ps (bmTaper);
	vt    = _mm256_loadu_ps(&args->sigma2[ichan]);
	vt    = _mm256_add_ps (vt, vf);   /* add bmTaper */
	vt1   = _mm256_mul_ps (vt1, vt);  /* sigma2 * u * u */
	vt    = _mm256_loadu_ps(&args->sigma1[ichan]);
	vt    = _mm256_add_ps (vt, vf);   /* add bmTaper */
	vt2   = _mm256_mul_ps (vt2, vt);  /* sigma1 * v * v */
	vt    = _mm256_loadu_ps(&args->sigma3[ichan]);
	vt3   = _mm256_mul_ps (vt3, vt);  /* sigma3 * u * v */
	vin.v = _mm256_add_ps (vt1, vt2);
	vin.v = _mm256_add_ps (vin.v, vt3);
	ObitExpVec(8, vin.f, vout.f); 
	vvalid = _mm256_mul_ps(vout.v, vvalid);
      }  /* end Other (MF) tapering */
      /* Beam taper only if given multiply by vvalid */
      else if (fabs(bmTaper)>0.0) {
	vt    = _mm256_set1_ps (bmTaper);
	vin.v = _mm256_mul_ps (vBL2, vt);
	ObitExpVec(8, vin.f, vout.f); 
	vvalid = _mm256_mul_ps(vout.v, vvalid);
      } /* end beam Taper */
    } /* end any tapering */
    /* Save valid*taper */
    _mm256_storeu_ps(&valid[ichan], vvalid); 
    /* Scale u, v to cells */
    /*chu.v = _mm256_mul_ps(chu.v, vuscale);*/
    /*chv.v = _mm256_mul_ps(chv.v, vvscale);*/
    /* Round u (always non negative) */
    /*NO chu.v = _mm256_add_ps(chu.v, _half);
         chu.v = _mm256_floor_ps(chu.v);*/
    /* AVX doesn't do integers */
    kchan = 0;
    for (jchan=ichan; jchan<ichan+8; jchan++) {
      u = chu.f[kchan]; v = chv.f[kchan]; kchan++;
      /* Flagged data? */
      if (valid[jchan]==0.0) {
	args->iuarr[jchan]  = 0;
	args->ivarr[jchan]  = halfv;
	args->cnvfnv[jchan] = convfn;
	args->cnvfnv[jchan] = convfn;
      } else { /* OK */
	iu = (olong)(u+0.5);   /* to grid cell number */
	iv = _lroundf(v);
	args->iuarr[jchan] = iu; /* u cell including extra halfWidth cells */
	args->ivarr[jchan] = iv + halfv - halfWidth;
	it = convNperCell + _lroundf(convNperCell * (iu - u - 0.5));
	args->cnvfnu[jchan] = convfn + fullWidth * it;
	it = convNperCell + _lroundf(convNperCell * (iv - v - 0.5));
	args->cnvfnv[jchan] = convfn + fullWidth * it;
      }
    } /* end integer channel loop */
  } /* end loop over blocks of 8 channels */
  /* do whatever is left scalar */
  for (ichan=ndone; ichan<eChan; ichan++) {
    freqFact = phaseSign * gridInfo->freqArr[ichan];
    u = uu * freqFact * gridInfo->uscale[ifacet];
    v = vv * freqFact * gridInfo->vscale[ifacet];
    w = ww * freqFact * gridInfo->uscale[ifacet];
   /* Guardband */
    if ((u<gridInfo->guardu[ifacet]) && (fabs(v)<gridInfo->guardv[ifacet])) 
      valid[ichan] = 1.0;
    else valid[ichan] = 0.0;
    /* MAX/MIN BL */
    BL2 = u*u + v*v;
    if ((maxBL2>0) && (BL2>maxBL2))  valid[ichan] = 0.0;
    if ((minBL2>0) && (BL2<minBL2))  valid[ichan] = 0.0;
    /* mask u,v to get out of bound data */
    if (valid[ichan]==0.0) {u = v = halfWidth;}
    /* Save u,v,w */
    saveu[ichan] = u;
    savev[ichan] = v;
    savew[ichan] = w;
    iu = (olong)(u+0.5);   /* to grid cell number */
    iv = _lroundf(v);
    args->iuarr[ichan] = iu; /* u cell including extra halfWidth cells */
    args->ivarr[ichan] = iv + halfv - halfWidth;
    it = convNperCell + _lroundf(convNperCell * (iu - u - 0.5));
    args->cnvfnu[ichan] = convfn + fullWidth * it;
    it = convNperCell + _lroundf(convNperCell * (iv - v - 0.5));
    args->cnvfnv[ichan] = convfn + fullWidth * it;
    /* Tapering? */
    if (doTape) {
      freqFact = phaseSign * gridInfo->freqArr[ichan];
      u = uu * freqFact;   /* u, v in wavelengths at channel */
      v = vv * freqFact;
     /* Other (MF) tapering */
      if ((args->sigma1!=NULL) && (args->sigma1[ichan]!=0.0)) {
	ftemp = u*u*(args->sigma2[ichan]+bmTaper) + 
	        v*v*(args->sigma1[ichan]+bmTaper) + 
	        u*v*args->sigma3[ichan];
	fact = ObitExpCalc (ftemp);
	valid[ichan] *= fact;
      }  /* end Other (MF) tapering */
      /* Beam taper if given multiply by valid */
      else if (fabs(bmTaper)>0.0) {
	fact = ObitExpCalc (bmTaper*BL2);
	valid[ichan] *= fact;
      }  /* end beam Taper */
    } /* end any tapering */
  } /* End remainder loop */

  /* If this a beam - don't bother shifting - replace data with (wt*1,0) */
  if (gridInfo->isBeam[ifacet]) {
    ndone = bChan;
    for (ichan=bChan; ichan<eChan; ichan+=8) {
      if ((ichan+8)>=eChan) break;  /* must have all 8 */
      ndone += 8;  /* how many channels done */
      /* data weight */
      jvis = ivis + gridInfo->nrparm + ichan*3;
      vw  = _mm256_set_ps(vis_in[jvis+23], vis_in[jvis+20], vis_in[jvis+17], vis_in[jvis+14],
			  vis_in[jvis+11], vis_in[jvis+8],  vis_in[jvis+5],  vis_in[jvis+2]);
      vt2 = _mm256_setzero_ps();
      /* Validity*taper mask */
      vvalid = _mm256_loadu_ps(&valid[ichan]);
      /* Modify weight for validity*taper */
      vw  = _mm256_mul_ps(vw, vvalid);
      /* Make vw nonnegative */
      vout.v  = _mm256_max_ps (vw, vt2);
      /* SW Beam order */
      switch (beamOrd) {
      case 0:   /* Dirty beam */
	break;
      case 1:   /* first order ln(nu/nu0) */
	kchan = 0;
	  for (jchan=ichan; jchan<ichan+8; jchan++) {
	    ftemp = log(gridInfo->freqArr[jchan]);
	    farr[kchan++] = ftemp;
	  }
	vt     = _mm256_loadu_ps(farr);
	vout.v = _mm256_mul_ps(vt, vout.v);
	break;
      case 2:   /* second order ln(nu/nu0)**2 */
	kchan = 0;
	  for (jchan=ichan; jchan<ichan+8; jchan++) {
	    ftemp = log(gridInfo->freqArr[jchan]);
	    farr[kchan++] = ftemp*ftemp;
	  }
	vt     = _mm256_loadu_ps(farr);
	vout.v = _mm256_mul_ps(vt, vout.v);
	break;
      case 3:   /* third order ln(nu/nu0)**3 */
	kchan = 0;
	  for (jchan=ichan; jchan<ichan+8; jchan++) {
	    ftemp = log(gridInfo->freqArr[jchan]);
	    farr[kchan++] = ftemp*ftemp*ftemp;
	  }
	vt     = _mm256_loadu_ps(farr);
	vout.v = _mm256_mul_ps(vt, vout.v);
	break;
      default:
	break;
      }; /* end beamOrd switch */
      /* Save weighted data in fwork1 */
      /* Damn - do it the hard way */
      /*vt  = _mm256_unpacklo_ps(vw, vt2);*/ /* Interleave first 4 r/i pairs */
      /*_mm256_storeu_ps(&args->fwork1[2*ichan], vt); */
      /*vt  = _mm256_unpackhi_ps(vw, vt2);*/ /* Interleave second 4 r/i pairs */
      /*_mm256_storeu_ps(&args->fwork1[2*ichan+8], vt); */
      for (jchan=0; jchan<8; jchan++) {
	args->fwork1[2*(ichan+jchan)]   = vout.f[jchan];
	args->fwork1[2*(ichan+jchan)+1] = 0.0;
      }
      
    } /* end loop over blocks of 8 channels */
    /* do whatever is left scalar */
    for (ichan=ndone; ichan<eChan; ichan++) {
      /* Weight, enforce guardband */
      jvis = ivis + gridInfo->nrparm + ichan*3;
      tw = MAX (0.0, vis_in[jvis+2]*valid[ichan]);
      /* SW Beam order */
      switch (beamOrd) {
      case 0:   /* Dirty beam */
	break;
      case 1:   /* first order ln(nu/nu0) */
	ftemp = log(gridInfo->freqArr[ichan]);
	tw *= ftemp;
	break;
      case 2:   /* second order ln(nu/nu0)**2 */
	ftemp = log(gridInfo->freqArr[ichan]);
	tw *= ftemp*ftemp;
	break;
      case 3:   /* third order ln(nu/nu0)**3 */
	ftemp = log(gridInfo->freqArr[ichan]);
	tw *= ftemp*ftemp*ftemp;
	break;
      default:
	break;
      }; /* end beamOrd switch */
      /* weighted data - save in fwork1 */
      args->fwork1[2*ichan]   = tw;
      args->fwork1[2*ichan+1] = 0.0;
    } /* end final channels */

  } else { /* not a beam */
    /* shift to double vector - scale to turns,cells */
    dshift[0] = (gdouble)shift[0] * (itwopi * iscaleu); 
    dshift[1] = (gdouble)shift[1] * (itwopi * iscalev); 
    dshift[2] = (gdouble)shift[2] * (itwopi * iscaleu); 
    ds1 = _mm256_set1_pd(dshift[0]);
    ds2 = _mm256_set1_pd(dshift[1]);
    ds3 = _mm256_set1_pd(dshift[2]);
    
    /* Channel loop shifting vis, first get phase shift, calc in double */
    /* No double AVX dot product */
    ndone = bChan;
    for (ichan=bChan; ichan<eChan; ichan+=4) {
      if ((ichan+4)>=eChan) break;  /* must have all 4 */
      ndone += 4;
      v4t = _mm_maskload_ps(&saveu[ichan], _MASK4);
      du  = _mm256_cvtps_pd(v4t);
      v4t = _mm_maskload_ps(&savev[ichan], _MASK4);
      dv  = _mm256_cvtps_pd(v4t);
      v4t = _mm_maskload_ps(&savew[ichan], _MASK4);
      dw  = _mm256_cvtps_pd(v4t);
      dt1     = _mm256_mul_pd(du, ds1);
      dt2     = _mm256_mul_pd(dv, ds2);
      dt3     = _mm256_mul_pd(dw, ds3);
      dt1     = _mm256_add_pd(dt1, dt2);
      dt1     = _mm256_add_pd(dt1, dt3);       /* u*s0+v*s1+w*s2 in turns */
      dt2     = _mm256_floor_pd(dt1);          /* subtract 2 n pi */
      dt1     = _mm256_sub_pd(dt1, dt2);
      dt2     = _mm256_set1_pd (twopi);
      dt1     = _mm256_mul_pd(dt1, dt2);       /* to radians */
      v4out   = _mm256_cvtpd_ps(dt1);          /* to float */
      _mm_maskstore_ps(&saveph[ichan], _MASK4, v4out); /* save */
    } /* end loop over blocks of 4 channels */
    /* do whatever is left scalar */
    for (ichan=ndone; ichan<eChan; ichan++) {
      u = saveu[ichan];  v = savev[ichan];  w = savev[ichan]; 
      dphase =  (u*dshift[0] + v*dshift[1] + w*dshift[2]);
      iphase = (olong)(dphase);
      saveph[ichan] = (ofloat)(dphase - iphase * twopi);
    } /* end final channels */
    
    /* Get sine/cosine of phase */
    ObitSinCosVec(gridInfo->nchan, saveph, saves, savec);
  
    /* Channel loop shifting vis, first get phase shift, calc in double */
    ndone = bChan;
    for (ichan=bChan; ichan<eChan; ichan+=8) {
      if ((ichan+8)>=eChan) break;  /* must have all 8 */
      ndone += 8;  /* how many channels done */
      
      /* Load next 8 vis to vectors */
      jvis = ivis + gridInfo->nrparm + ichan*3;
      vvr = _mm256_set_ps(vis_in[jvis+21], vis_in[jvis+18], vis_in[jvis+15], vis_in[jvis+12],
			  vis_in[jvis+9],  vis_in[jvis+6],  vis_in[jvis+3],  vis_in[jvis+0]);
      vvi = _mm256_set_ps(vis_in[jvis+22], vis_in[jvis+19], vis_in[jvis+16], vis_in[jvis+13],
			  vis_in[jvis+10], vis_in[jvis+7],  vis_in[jvis+4],  vis_in[jvis+1]);
      vw  = _mm256_set_ps(vis_in[jvis+23], vis_in[jvis+20], vis_in[jvis+17], vis_in[jvis+14],
			  vis_in[jvis+11], vis_in[jvis+8],  vis_in[jvis+5],  vis_in[jvis+2]);
      /* conjugate for negative u */
      vt  = _mm256_set1_ps(phaseSign);
      vvi = _mm256_mul_ps(vt, vvi);
      
      /* Corresponding cosine, sine of phase */
      vpr = _mm256_loadu_ps(&savec[ichan]);
      vpi = _mm256_loadu_ps(&saves[ichan]);
      /* Validity*taper mask */
      vvalid = _mm256_loadu_ps(&valid[ichan]);
      /* Modify weight for validity*taper */
      vw  = _mm256_mul_ps(vw, vvalid);
      /* Make nonnegative */
      vt  = _mm256_setzero_ps();
      vw  = _mm256_max_ps (vw, vt);
      /* Rotate vis */
      vt  = _mm256_mul_ps(vvr, vpr);
      vt1 = _mm256_mul_ps(vvi, vpi);
      vt1 = _mm256_sub_ps(vt, vt1);   /* Real part */
      vt  = _mm256_mul_ps(vvr, vpi);
      vt2 = _mm256_mul_ps(vvi, vpr);
      vt2 = _mm256_add_ps(vt, vt2);   /* Imaginary part */
      /* weight data */
      vin.v  = _mm256_mul_ps(vt1, vw);
      vout.v = _mm256_mul_ps(vt2, vw);
      /* Save weighted data in fwork1 */
       /* Damn - do it the hard way */
      /*vt  = _mm256_unpacklo_ps(vt1, vt2);*/ /* Interleave first 4 r/i pairs */
      /*_mm256_storeu_ps(&args->fwork1[2*ichan], vt); */
      /*vt  = _mm256_unpackhi_ps(vt1, vt2);*/ /* Interleave second 4 r/i pairs */
      /*_mm256_storeu_ps(&args->fwork1[2*ichan+8], vt); */
      for (jchan=0; jchan<8; jchan++) {
	args->fwork1[2*(ichan+jchan)]   = vin.f[jchan];
	args->fwork1[2*(ichan+jchan)+1] = vout.f[jchan];
      }
    } /* end loop over blocks of 8 channels */
    /* do whatever is left scalar */
    for (ichan=ndone; ichan<eChan; ichan++) {
      c = savec[ichan]; s = saves[ichan];
      jvis = ivis + gridInfo->nrparm + ichan*3;
      tr = vis_in[jvis]; ti = phaseSign*vis_in[jvis+1]; 
      /* Weight, enforce guardband */
      tw = MAX (0.0, vis_in[jvis+2]*valid[ichan]);
      vr = c*tr - s*ti;
      vi = s*tr + c*ti;
      /* weighted data - save in fwork1 */
      args->fwork1[2*ichan]   = vr * tw;
      args->fwork1[2*ichan+1] = vi * tw;
    } /* end final channels */
  } /* end not a beam */
} /* end fast_prep_grid */

/** 
 * Fast AVX  gridding  only 7x7 implemented
 * \param grid  base of visibility grid
 * \param vis   weighted visibility (r,i)
 * \param iu    u col. number of start of convolution kernal (as ofloat)
 * \param iv    v row number of start of convolution kernal
 * \param lrow  length of row in ofloats
 * /param nconv dimension in u,v, of convolution kernal (better be 7)
 * \param cu    start of u convolution
 * \param cv    start of v convolution
*/
void fast_grid(ofloat *grid, ofloat vis[2], olong iu, olong iv, olong lrow, 
		olong nconv,  ofloat *cu, ofloat *cv)  
{
  fast_grid7 (grid, vis, iu, iv, lrow, nconv, cu, cv);
} /* end fast_grid */
/** 
 * Fast AVX (8) 7x7 gridding
 * \param grid  base of visibility grid
 * \param vis   weighted visibility (r,i)
 * \param iu    u col. number of start of convolution kernal (as ofloat)
 * \param iv    v row number of start of convolution kernal
 * \param lrow  length of row in ofloats
 * /param nconv dimension in u,v, of convolution kernal (better be 7)
 * \param cu    start of u convolution
 * \param cv    start of v convolution
*/
void fast_grid7(ofloat *grid, ofloat vis[2], olong iu, olong iv, olong lrow, 
		olong nconv,  ofloat *cu, ofloat *cv)  
{
  v8sf vt1, vt2, vcu1, vcu2, vcv, vs, vt, vgrid, vconv;
  V8SF vtmp; 
  int jv, addr, ia;

  if (nconv<=0) return;
  if ((vis[0]==0.0) && (vis[1]==0.0)) return;
  /* Load convolution kernals, copy pairs of values to first and second halves */
  /* Have to do this the hard way */
  vtmp.v = _mm256_loadu_ps(cu);
  vcu1   = _mm256_set_ps(vtmp.f[3], vtmp.f[3], vtmp.f[2], vtmp.f[2], vtmp.f[1], vtmp.f[1], vtmp.f[0], vtmp.f[0]);
  vcu2   = _mm256_set_ps( 0.0,       0.0,      vtmp.f[6], vtmp.f[6], vtmp.f[5], vtmp.f[5], vtmp.f[4], vtmp.f[4]);

  /* Create first and second parts of vis vector x conv u */
  vt1 = _mm256_set1_ps(vis[0]);   /* Real */
  vt2 = _mm256_set1_ps(vis[1]);   /* Imaginary */
  vs  = _mm256_blend_ps(vt1, vt2, 0xaa);   /* Vis, pairs of r,i */
  ia = iv*lrow + iu*2;
  /* Loop in v */
  for (jv=0; jv<nconv; jv++) {
    addr = ia;  /* Address of grid */
    /* Load first half of grid */
    vgrid = _mm256_loadu_ps(&grid[addr]);
    vcv   = _mm256_set1_ps(cv[jv]);
    vconv = _mm256_mul_ps(vcu1, vcv);  /* Convolution function */
    /* Multiply vis by convolution fn */
    vt    = _mm256_mul_ps(vs, vconv);
    /* Update first half of grid */
    vgrid = _mm256_add_ps(vt, vgrid);
    _mm256_storeu_ps(&grid[addr], vgrid);
    /* Second part of grid */
    addr += 8;
    /*vgrid = _mm256_maskload_ps(&grid[addr], _mask6);*/ /*only 6 */
    vgrid = _mm256_loadu_ps(&grid[addr]);
    /* Multiply vis by convolution fn */
    vconv = _mm256_mul_ps(vcu2, vcv);
    vt    = _mm256_mul_ps(vs, vconv);
    /* Update second half of grid */
    vgrid = _mm256_add_ps(vt, vgrid);
   /* _mm256_maskstore_ps(&grid[addr], _mask6, vgrid); */ /*  Only 6 */
    _mm256_storeu_ps(&grid[addr], vgrid);
   ia += lrow;
  } /* end loop in v */

  /* _mm_empty();  wait for operations to finish */
} /* end fast_grid7 */
#else  /* C only versions */
/** 
 * c version grid for any nconv
 * \param grid  base of visibility grid
 * \param vis   visibility (r,i)
 * \param iu    u col. (0-rel) number of start of convolution kernal (as ofloat)
 * \param iv    v row (0-rel) number of start of convolution kernal
 * \param lrow  length of row in ofloats
 * /param nconv dimension in u,v, of convolution kernal
 * \param cu    start of u convolution
 * \param cv    start of v convolution
*/
void fast_grid(ofloat *grid, ofloat vis[2], olong iu, olong iv, olong lrow, 
	       olong nconv, ofloat *cu, ofloat *cv)  
{
  olong ju, jv, addr;
  ofloat cvv, cfn;
  if (nconv<=0) return;
  /* hard core debug 
  if ((ivis==36)&&(ichan==1)&&(facet==0)) {
    fprintf (stdout,"dbg %d %d %d %d %d  %15.5f %15.5f %15.5f   \n", 
	     iu, iv, facet, ivis, ichan, vis[0], cu[0], cv[0]);
    fprintf (stdout,"    %15.5f %15.5f %15.5f %15.5f %15.5f %15.5f %15.5f %15.5f   \n", 
	     cv[0], cv[1], cv[2], cv[3], cv[4], cv[5], cv[6], cv[7]);
  }*/
  if ((vis[0]==0.0) && (vis[1]==0.0)) return;
  for (jv=0; jv<nconv; jv++) {
    cvv = cv[jv];
    if (cvv!=0.0) {
      addr = (iv+jv)*lrow + iu*2;
      for (ju=0; ju<nconv; ju++) {
	cfn = cu[ju]*cvv;
	grid[addr++] += vis[0] * cfn;
	grid[addr++] += vis[1] * cfn;
     }  /* end u loop */
    }
  }  /* end v loop */
} /* end fast_grid */
#endif /* end c version */
