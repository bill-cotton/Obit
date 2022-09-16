/* $Id$        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2014-2022                                          */
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

#include "ObitGPUGrid.h"
#include "ObitUVGrid.h"
#include "ObitUVGridMF.h"
#include "ObitImageMF.h"
#include "ObitImageWB.h"
#if HAVE_GPU==1  /* Compiled with GPU?*/
#include "ObitCUDAUtil.h"
size_t ObitCUDAUtilMemory(int cuda_device); /* Huh???*/
void ObitCUDAUtilHostAny2GPU(void *GPU, void *host, int memsize, int* stream);
#include "ObitCUDAGridInfoDef.h"
#endif  /* Compiled with  GPU?*/

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitGPUGrid.c
 * ObitGPUGrid class function definitions.
 * This class is derived from the Obit base class.
 * The ObitGPUGrid class uses GPUs to calculate sky models.
 * CUDA is used for the actual GPU routines but CUDA is REALLY primitive
 * so very basic implementation is in ObitCUDAGrid.cu and friends.
 * Portions of the class are in CUDA and are only implemented if the
 * compiler option -DHAVE_GPU=1 is used.  Some portions also need to 
 * variable IS_CUDA=1 to be set in the calling routines.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitGPUGrid";

/** Function to obtain parent ClassInfo */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo structure ObitGPUGridClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitGPUGridClassInfo myClassInfo = {FALSE};

/*--------------- File Global Variables  ----------------*/


/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitGPUGridInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitGPUGridClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitGPUGridClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitGPUGrid* newObitGPUGrid (gchar* name)
{
  ObitGPUGrid* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitGPUGridClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitGPUGrid));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitGPUGridInit((gpointer)out);

 return out;
} /* end newObitGPUGrid */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitGPUGridGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitGPUGridClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitGPUGridGetClass */

/**
 * Make a deep copy of an ObitGPUGrid.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitGPUGrid* ObitGPUGridCopy  (ObitGPUGrid *in, ObitGPUGrid *out, ObitErr *err)
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
    out = newObitGPUGrid(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */

  return out;
} /* end ObitGPUGridCopy */

/**
 * Make a copy of a object but do not copy the actual data
 * This is useful to create an GPUGrid similar to the input one.
 * \param in  The object to copy
 * \param out An existing object pointer for output, must be defined.
 * \param err Obit error stack object.
 */
void ObitGPUGridClone  (ObitGPUGrid *in, ObitGPUGrid *out, ObitErr *err)
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

} /* end ObitGPUGridClone */

/**
 * Creates an ObitGPUGrid 
 * \param name    An optional name for the object.
 * \param nfacet  Total number of facets
 * \param iimage  output Image as Obit
 * \param UVin    Input data
 * \return the new object.
 */
ObitGPUGrid* ObitGPUGridCreate (gchar* name, olong nfacet, Obit *iimage, ObitUV *UVin)
{
  ObitGPUGrid* out;
  ObitImage *image = (ObitImage*)iimage;
  ObitUVDesc *uvDesc = UVin->myDesc;
  ObitImageDesc *imDesc = image->myDesc;
  olong nchan, nif, nVisPIO, lenvis, nplane;
  /* ObitInfoType type;
     gint32 dim[MAXINFOELEMDIM];*/
 
  /* Create basic structure */
  out = newObitGPUGrid (name);
  out->cudaInfo = g_malloc0(sizeof(CUDAGridInfo));

  /* Number of vis per IO */
  nVisPIO = MIN(GPU_NVISPIO, uvDesc->nvis);
  /*ObitInfoListGetTest (UVin->info, "nVisPIO", &type, dim, &nVisPIO);*/
  /*  ObitInfoListAlwaysPut (UVin->info, "nVisPIO", type, dim, &nVisPIO);*/
  nchan = uvDesc->inaxes[uvDesc->jlocf];  /* number of channels */
  lenvis = uvDesc->lrec;
  if (uvDesc->jlocif>0) nif = uvDesc->inaxes[uvDesc->jlocif];  /* number of IFs */
  else                  nif = 1;

  /* How many planes? */
  if (ObitImageWBIsA(image)) {
       nplane = 1;  /* Note GPU gridding doesn't work for this */
   } else if (ObitImageMFIsA(image)) {
      nplane = ((ObitImageMF*)image)->nSpec;
    } else if (imDesc->inaxes[imDesc->jlocf]>1) {
      nplane = 1;  /* Note GPU gridding may not work for this */
    } else {  /* Basic image */
      nplane = 1;
    }
 
  /* Init - allocate structures */
#if HAVE_GPU==1  /* Compiled with GPU?*/
  ObitCUDAGridAlloc(out->cudaInfo, (long)nfacet, (long)nchan, (long)nif, 
		    (long)nVisPIO, (long)lenvis, (long)nplane);
#endif  /* Compiled with  GPU?*/

  return out;
} /* end ObitGPUGridCreate */

/**
 * Set GPU to use, return GPU memory

 * \param cuda_device Which GPU?
 * \return  GPU memory in bytes
 */
ollong ObitGPUGridSetGPU (int cuda_device)
{
  ollong out=0;
#if HAVE_GPU==1  /* compiled with GPU?*/
  ObitCUDASetGPU (cuda_device);
  out = (ollong)ObitCUDAUtilMemory(cuda_device);
#endif /* compiled with GPU?*/
  return out;
} /* end ObitGPUGridSetGPU*/

/**
 * Copies gridding info to GPU (CUDA) structure, sets up structures.
 * Copy to GPU
 * Must be called per facet.
 * The object UVin will be opened during this call if it is not already open.
 * The output image should describe the center, size and grid spacing of the desired
 * image.
 * \param in       Object to initialize
 * \param uvgrid   ObitUVGrid with parameters as (Obit*)
 * \param chDone   Array of nplane flags indicating which channels are done.
 *                 If NULL then no channels are "done".
 * \param ifacet   0-rel Facet number, initialize on 0, only even if doBeam
 * \param nfacet   Number of facets
 * \param nplane   Number of planes in grid
 * \param UVin     Uv data object to be gridded.
 * \param imagee   Image to be gridded (as Obit*)
 *                 Descriptor infoList entry "BeamTapr" gives any additional
 *                 tapering in degrees.
 * \param doBeam   TRUE if also make Beams,image in alternating entries
 * \param err      ObitErr stack for reporting problems.
 */
void ObitGPUGridSetGPUStruct (ObitGPUGrid *in, Obit *uvgrid, gboolean *chDone,
			      olong ifacet, olong nfacet, olong nplane, 
			      ObitUV *UVin, Obit *imagee, gboolean doBeam, ObitErr *err)
{
  olong nVisPIO, oldnVisPIO;
  ObitUVDesc *uvDesc = UVin->myDesc;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
#if HAVE_GPU==1  /* compiled with GPU?*/

  ObitImageDesc *imDesc = ((ObitImage*)imagee)->myDesc;
  ObitIOAccess access;
  gboolean doCalSelect;
  olong i, j, n, nchan, nx, ny, nchanIF, nSpec, last;
  odouble *specFreq=NULL, *specFreqLo=NULL, *specFreqHi=NULL;
  gchar keyword[24];
  ObitUVGrid *imGrid=(ObitUVGrid*)((ObitImage*)imagee)->myGrid, *uvGrid=(ObitUVGrid*)uvgrid;
  ObitUVGridMF *uvGridMF=(ObitUVGridMF*)uvgrid;
  CUDAGridInfo* gpuInfo;
  FacetInfo** facetInfo;
  GridInfo* gridInfo;
#endif /* compiled with GPU?*/
  gchar *routine="ObitGPUGridSetGPUStruct";

  /* Number of vis per IO */
  nVisPIO = MIN(GPU_NVISPIO, uvDesc->nvis);
  ObitInfoListGetTest (UVin->info, "nVisPIO", &type, dim, &nVisPIO);
  oldnVisPIO = nVisPIO;

  /* to be sure */
#if HAVE_GPU==1  /* compiled with GPU?*/
  /* Init */
  nchan = uvDesc->inaxes[uvDesc->jlocf];
  nchanIF = nchan * uvDesc->inaxes[uvDesc->jlocif];
  if (ifacet==0) {  /* create/allocate */
    gpuInfo = in->cudaInfo;  /* pointer to host structure */
    gpuInfo->gpu_memory = ObitCUDAUtilMemory((int)in->cuda_device);
    gridInfo = gpuInfo->h_gridInfo;   /* Host Info common to all grids */
    gridInfo->nfacet     = nfacet;
    gpuInfo->d_base = (void*)ObitCUDAUtilAllocGPU(sizeof (CUDAGridInfo));

    gridInfo->nplane     = nplane;
    gridInfo->nchan      = nchanIF;  /* Lower levels don't care about IF */
    gridInfo->nif        = 1;
    gridInfo->nstok      = uvDesc->inaxes[uvDesc->jlocs];
    gridInfo->nrparm     = uvDesc->nrparm;
    gridInfo->lenvis     = uvDesc->lrec;
    /* From ObitUVGrid */
    gridInfo->convWidth  = uvGrid->convWidth;
    gridInfo->convNperCell = uvGrid->convNperCell;
    gridInfo->guardu     = 0.4; /* Fraction of grid */
    gridInfo->guardv     = 0.4;
    gridInfo->rotate     = uvGrid->rotate;
    for (i=0; i<nchanIF; i++) gridInfo->h_freqArr[i] = (float)uvDesc->fscale[i];
    n = uvGrid->convWidth * uvGrid->convNperCell;
    for (i=0; i<n; i++) gridInfo->h_convfn[i] = uvGrid->convfn->array[i];

    facetInfo = gpuInfo->h_facetInfo;   /* Array of per facet structs */;
    gridInfo->prtLv = err->prtLv; /* Print level */

    /* For ImageMF set freqPlane */
    if (ObitImageMFIsA((ObitImage*)imagee)) {
      /* get binning from image header */
      nSpec = -1;
      ObitInfoListGetTest(imDesc->info, "NSPEC", &type, dim, &nSpec);
      specFreq   = (odouble*)g_malloc0(nSpec*sizeof(odouble));
      specFreqLo = (odouble*)g_malloc0(nSpec*sizeof(odouble));
      specFreqHi = (odouble*)g_malloc0(nSpec*sizeof(odouble));
      for (i=0; i<nSpec; i++) {
	specFreq[i] = 1.0;
	sprintf (keyword, "FREQ%4.4d",i+1);
	ObitInfoListGetTest(imDesc->info, keyword, &type, dim, &specFreq[i]);
	sprintf (keyword, "FREL%4.4d",i+1);
	ObitInfoListGetTest(imDesc->info, keyword, &type, dim, &specFreqLo[i]);
	sprintf (keyword, "FREH%4.4d",i+1);
	ObitInfoListGetTest(imDesc->info, keyword, &type, dim, &specFreqHi[i]);
      }
      last = 0;
      for (i=0; i<nchanIF; i++) {
	gridInfo->h_freqPlane[i] = last;  /* For first */
	for (j=0; j<nSpec; j++) {
	  if (uvDesc->freqArr[i] >= specFreqLo[j]) gridInfo->h_freqPlane[i] = j;
	}
	last = gridInfo->h_freqPlane[i];
	gridInfo->h_freqPlane[nchanIF-1] = nSpec-1; /* Highest */
      }
      /* Cleanup */
      if (specFreq)   g_free(specFreq);
      if (specFreqLo) g_free(specFreqLo);
      if (specFreqHi) g_free(specFreqHi);
      /* Taper sigmas */
      for (i=0; i<nchanIF; i++) {
	gridInfo->h_sigma1[i] = uvGridMF->sigma1[i];
	gridInfo->h_sigma2[i] = uvGridMF->sigma2[i];
	gridInfo->h_sigma3[i] = uvGridMF->sigma3[i];
      }
    } else if (nplane<=1) { /* single continuum plane */
      for (i=0; i<nchanIF; i++) gridInfo->h_freqPlane[i] = 0;
      for (i=0; i<nchanIF; i++) 
	gridInfo->h_sigma1[i] = gridInfo->h_sigma2[i] = gridInfo->h_sigma3[i] = 0;
   } else {  /* assume line image */
      for (i=0; i<nchanIF; i++) gridInfo->h_freqPlane[i] = i;
      for (i=0; i<nchanIF; i++) 
	gridInfo->h_sigma1[i] = gridInfo->h_sigma2[i] = gridInfo->h_sigma3[i] = 0;
    }
      /* Are any of the grids already done? Give invalid plane number.  CHECK*/
    if (chDone) {
      for (i=0; i<nchanIF; i++) {
	if (chDone[i]) gridInfo->h_freqPlane[i] = nplane+100;
      }
    } /* end channel done flag */
  
  } /* end create/allocate */

  gpuInfo   = in->cudaInfo;           /* pointer to host structure */
  facetInfo = gpuInfo->h_facetInfo;   /* Array of per facet structs */;
  gridInfo  = gpuInfo->h_gridInfo;    /* Host Info common to all grids */
    
  /* Grid size 
     Grids have half convolving fn width extra x cells */
  nx = imDesc->inaxes[0]; ny = imDesc->inaxes[1]; 
  facetInfo[ifacet]->sizeGrid = 2*(nx/2+1+gridInfo->convWidth/2)*ny;
  /* Plus some slop */
  facetInfo[ifacet]->sizeGrid += nx;

  facetInfo[ifacet]->uscale     = uvGrid->UScale; /* facet (nx,ny) specific */
  facetInfo[ifacet]->vscale     = uvGrid->VScale;
  facetInfo[ifacet]->wscale     = uvGrid->WScale;
  /* If doBeam alternating beam/image pairs
     per facet stuff */
  facetInfo[ifacet]->nx = nx;
  facetInfo[ifacet]->ny = ny;
  if (!doBeam) { /* don't shift beam */
    facetInfo[ifacet]->shift[0] = imGrid->dxc;
    facetInfo[ifacet]->shift[1] = imGrid->dyc;
    facetInfo[ifacet]->shift[2] = imGrid->dzc;
  }
  facetInfo[ifacet]->rotUV[0] = uvGrid->URot3D[0][0];
  facetInfo[ifacet]->rotUV[1] = uvGrid->URot3D[0][1];
  facetInfo[ifacet]->rotUV[2] = uvGrid->URot3D[0][2];
  facetInfo[ifacet]->rotUV[3] = uvGrid->URot3D[1][0];
  facetInfo[ifacet]->rotUV[4] = uvGrid->URot3D[1][1];
  facetInfo[ifacet]->rotUV[5] = uvGrid->URot3D[1][2];
  facetInfo[ifacet]->rotUV[6] = uvGrid->URot3D[2][0];
  facetInfo[ifacet]->rotUV[7] = uvGrid->URot3D[2][1];
  facetInfo[ifacet]->rotUV[8] = uvGrid->URot3D[2][2];
  facetInfo[ifacet]->doBeam = doBeam;
  facetInfo[ifacet]->maxBL   = uvGrid->blmax;
  if (uvGrid->blmax<=0.0) facetInfo[ifacet]->maxBL = 1.0e12;
  facetInfo[ifacet]->minBL   = uvGrid->blmin;
  if (ObitImageMFIsA((ObitImage*)imagee)) 
    facetInfo[ifacet]->bmTaper = uvGrid->BeamTaperUV;
  /* reallocate data buffer on UVin to locked memory
    first reset number of vis per IO */
  if (ifacet==0) {
    ObitUVClose (UVin, err);
    gridInfo->oldnVisPIO = oldnVisPIO;  /* Save */
    /*nVisPIO = MIN(MIN(10000,oldnVisPIO),  uvDesc->nvis);
      dim[0] = dim[1] = dim[2] = dim[3] = 1;
      ObitInfoListAlwaysPut (UVin->info, "nVisPIO", OBIT_long, dim,  &nVisPIO);*/
    gridInfo->nvis = nVisPIO;  /* How many per read? */
    doCalSelect = FALSE;
    ObitInfoListGetTest(UVin->info, "doCalSelect", &type, dim, &doCalSelect);
    if (doCalSelect) access = OBIT_IO_ReadCal;
    else access = OBIT_IO_ReadOnly;   
    ObitUVOpen (UVin, access, err);
    if (err->error) Obit_traceback_msg (err, routine, UVin->name);
  }
  /* Copy to GPU ofter last facet */
if (ifacet>=(nfacet-1)) ObitCUDAGridInit(gpuInfo);
#endif /* HAVE_GPU */

} /* end  ObitGPUGridSetGPUStruct */


/**
 * Copy a bufferload of visibilities to the locked memory
 * \param in      Input Grid object
 * \param uvdata  UV data being operated on
 * \param err     Obit error stack object.
 */
void ObitGPUGrid2GPU (ObitGPUGrid *in, ObitUV *uvdata, ObitErr *err)
{
#if HAVE_GPU==1  /* Compiled with GPU?*/
  size_t memsize = uvdata->myDesc->numVisBuff * uvdata->myDesc->lrec * sizeof(float);  // size of vis buffer
  memcpy(in->cudaInfo->h_gridInfo->h_data_in, uvdata->buffer, memsize);
  /* How much actually read? */
  in->cudaInfo->h_gridInfo->nvis = uvdata->myDesc->numVisBuff;
#endif  /* Compiled with  GPU?*/
} /* end ObitGPUGrid2GPU */

/**
 * Grid a bufferload of visibilities
 * \param in      Input object
 * \param err     Obit error stack object.
 */
void ObitGPUGridGrid (ObitGPUGrid *in, ObitErr *err)
{
#if HAVE_GPU==1  /* Compiled with GPU?*/
  /* CALL CUDA */
  ObitCUDAGridGrid (in->cudaInfo);
#endif  /* Compiled with  GPU?*/
} /* end ObitGPUGridGrid */

/**
 * Fold negative u to positive
 * \param in      Input object
 * \param err     Obit error stack object.
 */
void ObitGPUGridFlip (ObitGPUGrid *in, ObitErr *err)
{
#if HAVE_GPU==1  /* Compiled with GPU?*/
  /* CALL CUDA */
  ObitCUDAGridFlip (in->cudaInfo);
#endif  /* Compiled with  GPU?*/
} /* end ObitGPUGridFlip */

/**
 * Copy grids to FFT buffers via locked memory
 * \param in      Input Grid object
 * \param ifacet  0-rel facet number
 * \param grid    array of output Grid/FFT array
 * \param uvdata  UV data being operated on
 * \param err     Obit error stack object.
 */
void ObitGPUGrid2Host (ObitGPUGrid *in, olong ifacet, ObitCArray **grid, ObitErr *err)
{
#if HAVE_GPU==1  /* Compiled with GPU?*/
  size_t memsize = grid[0]->naxis[0] * 2 * sizeof(float);  // size of Grid row in bytes
  olong iplane, iy, lirow, lorow, halfWidth, indx;
  /* Copy to host memory - loop over planes extracting to grid*/
  ObitCUDAGrid2CPU (in->cudaInfo, (long)ifacet);
  halfWidth  = in->cudaInfo->h_gridInfo->convWidth/2;  /* convolving fn halfwidth (width odd) */
  lorow = grid[0]->naxis[0];
  lirow = lorow + halfWidth;
  for (iplane=0; iplane<in->cudaInfo->h_gridInfo->nplane; iplane++) {
    /* copy to FFT CArray buffers */
    /* Trim rows by half width*2 floats */
    indx = iplane * in->cudaInfo->h_facetInfo[ifacet]->sizeGrid;
    for (iy=0; iy<grid[0]->naxis[1]; iy++) {
      memcpy(&grid[iplane]->array[2*iy*lorow], 
	     &in->cudaInfo->h_facetInfo[ifacet]->h_grid[indx+2*(iy*lirow+halfWidth)], memsize);
    } /* end row loop */
    /* Swaparonie */
    ObitCArray2DCenter(grid[iplane]);
  } /* end plane loop */
#endif  /* Compiled with  GPU?*/
} /* end ObitGPUGrid2Host */

/**
 * Shutdown  ObitGPUGrid 
 * \param in      Iniput object
 * \param uvdata  UV data being operated on
 * \param err     Obit error stack object.
 */
void ObitGPUGridShutdown (ObitGPUGrid *in, ObitUV *uvdata, ObitErr *err)
{
  if (err->error) return;
  if (in->cudaInfo==NULL) return;
#if HAVE_GPU==1  /* Compiled with GPU?*/
  /*int oldnVisPIO;
    gint32 dim[MAXINFOELEMDIM];
    oldnVisPIO = in->cudaInfo->h_gridInfo->oldnVisPIO;*/
  // Free resources GPU and locked memory
  ObitCUDAGridShutdown(in->cudaInfo);
  g_free(in->cudaInfo); in->cudaInfo = NULL;

  /* Replace old nVisPIO
  dim[0] = dim[1] = dim[2] = dim[3] = 1;
  ObitInfoListAlwaysPut (uvdata->info, "nVisPIO", OBIT_long, dim,  
			 &oldnVisPIO); */
#endif  /* Compiled with  GPU?*/
} /*end ObitGPUGridShutdown */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitGPUGridClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitGPUGridClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitGPUGridClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitGPUGridClassInfoDefFn (gpointer inClass)
{
  ObitGPUGridClassInfo *theClass = (ObitGPUGridClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitGPUGridClassInit;
  theClass->newObit       = (newObitFP)newObitGPUGrid;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitGPUGridClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitGPUGridGetClass;
  theClass->ObitCopy      = (ObitCopyFP)ObitGPUGridCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitGPUGridClear;
  theClass->ObitInit      = (ObitInitFP)ObitGPUGridInit;
  theClass->ObitGPUGridCreate = (ObitGPUGridCreateFP)ObitGPUGridCreate;
  theClass->ObitGPUGridSetGPU = (ObitGPUGridSetGPUFP)ObitGPUGridSetGPU;
  theClass->ObitGPUGridSetGPUStruct = (ObitGPUGridSetGPUStructFP)ObitGPUGridSetGPUStruct;
  /*theClass->ObitGPUGridInitGPU = (ObitGPUGridInitGPUFP)ObitGPUGridInitGPU;*/
  theClass->ObitGPUGrid2GPU = (ObitGPUGrid2GPUFP)ObitGPUGrid2GPU;
  theClass->ObitGPUGridGrid = (ObitGPUGridGridFP)ObitGPUGridGrid;
  theClass->ObitGPUGridFlip = (ObitGPUGridFlipFP)ObitGPUGridFlip;
  theClass->ObitGPUGrid2Host = (ObitGPUGrid2HostFP)ObitGPUGrid2Host;
  theClass->ObitGPUGridShutdown = (ObitGPUGridShutdownFP)ObitGPUGridShutdown;
} /* end ObitGPUGridClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitGPUGridInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitGPUGrid *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->cudaInfo = NULL;
  /* Not here
  in->cudaInfo = g_malloc0(sizeof(CUDAGridInfo));
  in->cudaInfo->h_gridInfo = NULL;
  in->cudaInfo->d_gridInfo = NULL;
  in->cudaInfo->h_facetInfo = NULL;
  in->cudaInfo->d_facetInfo = NULL;*/
} /* end ObitGPUGridInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitGPUGrid* cast to an Obit*.
 */
void ObitGPUGridClear (gpointer inn)
{
  ObitClassInfo *ParentClass;

#if HAVE_GPU==1  /* compiled with  GPU?*/
  ObitGPUGrid *in = inn;
  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  /* CUDA related things deleted in ObitGPUGridShutdown */
   
  if (in->cudaInfo) {g_free(in->cudaInfo);} in->cudaInfo = NULL;
#endif  /* Compiled with  GPU?*/
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitGPUGridClear */


