/* $Id: $   */
/* HIDE c (esp.glib) structures from cuda */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2026                                               */
/*;  Associated Universities, Inc. Washington DC, USA.                */
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
/*  Define the CUDA SkyGeom utilities                                 */
/**
 * \file ObitGPUSkyGeom.c
 * implementation of limited GPU version of ObitSkyGeom class
 */
#include "ObitImageDesc.h"
#include "ObitUVDesc.h"
#include "ObitGPUSkyGeom.h"
#include "ObitCUDAUtil.h"

#if HAVE_GPU==1  /* GPU? */
/**
 * Creates a CUDAImageDesc
 * \return the new object in locked host memory
 */
CUDAImageDesc* ObitGPUImageDescCreate ()
{
  CUDAImageDesc* out=NULL;
  int memsize;

  memsize = sizeof(CUDAImageDesc);
  out = (CUDAImageDesc*)ObitCUDAUtilAllocHost(memsize);
  //out = (CUDAImageDesc*)malloc(sizeof(CUDAImageDesc));
  return out;
} /* end ObitGPUImageDescCreate */

/**
 * Destroys a CUDAImageDesc
 * Called from c
 * \param in   Object to delete
 */
void ObitGPUImageDescZap (CUDAImageDesc *in)
{
  if (in) ObitCUDAUtilFreeHost((float*)in);	
} /* end ObitGPUImageDescZap */

/**
 * Creates a CUDAUVDesc
 * \return the new object in locked host memory
 */
CUDAUVDesc* ObitGPUUVDescCreate ()
{
  CUDAUVDesc* out=NULL;
  int memsize;

  memsize = sizeof(CUDAImageDesc);
  out = (CUDAUVDesc*)ObitCUDAUtilAllocHost(memsize);
 // out = (CUDAUVDesc*)malloc(sizeof(CUDAUVDesc));
  return out;
} /* end ObitGPUUVDescCreate */

/**
 * Destroys a CUDAUVDesc
 * Called from c
 * \param in   Object to delete
 */
void ObitGPUUVDescZap (CUDAUVDesc *in)
{
  if (in) ObitCUDAUtilFreeHost((float*)in);	
} /* end ObitGPUUVDescZap */

/** Public: Host to Device Image Descriptor conversion
    Create host resident version of a device Image Descriptor
  * \param  in
*/
CUDAImageDesc* ObitGPUSkyGeomImageH2D (ObitImageDesc *in)
{
  CUDAImageDesc* out=NULL;
  int i,j;
  out = ObitGPUImageDescCreate ();
  out->bitpix   = (long)in->bitpix;
  out->naxis   = (long)in->naxis;
  for (i=0; i<IM_MAXDIM; i++) out->inaxes[i] = (long)in->inaxes[i];
  for (i=0; i<IM_MAXDIM; i++) {
    for (j=0; j<IMLEN_KEYWORD; j++) out->ctype[i][j] = in->ctype[i][j];
  }
  for (i=0; i<IM_MAXDIM; i++) out->cdelt[i] = (float)in->cdelt[i];
  for (i=0; i<IM_MAXDIM; i++) out->crpix[i] = (float)in->crpix[i];
  for (i=0; i<IM_MAXDIM; i++) out->crota[i] = (float)in->crota[i];
  for (i=0; i<IM_MAXDIM; i++) out->crval[i] = (double)in->crval[i];
  out->altCrpix = (float)in->altCrpix;
  out->xshift   = (float)in->xshift;
  out->yshift   = (float)in->yshift;
  out->beamMaj  = (float)in->beamMaj;
  out->beamMin  = (float)in->beamMin;
  out->beamPA   = (float)in->beamPA;
  for (i=0; i<IMLEN_VALUE; i++) out->object[i]     = in->object[i];
  for (i=0; i<IMLEN_VALUE; i++) out->teles[i]      = in->teles[i];
  for (i=0; i<IMLEN_VALUE; i++) out->instrument[i] = in->instrument[i];
  for (i=0; i<IMLEN_VALUE; i++) out->observer[i]   = in->observer[i];
  for (i=0; i<IMLEN_VALUE; i++) out->obsdat[i]     = in->obsdat[i];
  for (i=0; i<IMLEN_VALUE; i++) out->date[i]       = in->date[i];
  for (i=0; i<IMLEN_VALUE; i++) out->origin[i]     = in->origin[i];
  for (i=0; i<IMLEN_VALUE; i++) out->bunit[i]      = in->bunit[i];
  out->epoch   = (float)in->epoch;
  out->equinox = (float)in->equinox;
  out->obsra   = (double)in->obsra;
  out->obsdec  = (double)in->obsdec;
  out->altRef   = (double)in->altRef;
  out->restFreq = (double)in->restFreq;
  out->minval   = (float)in->minval;
  out->maxval   = (float)in->maxval;
  if (in->areBlanks) out->areBlanks = 1;
  else               out->areBlanks = 0;
  out->niter   = (long)in->niter;
  out->VelReference   = (long)in->VelReference;
  out->VelDef   = (long)in->VelDef;
  /* fooey  out->coordType= in->coordType;*/
  out->row      = (long)in->row;
  out->plane    = (long)in->plane;
  out->plane4   = (long)in->plane4;
  out->plane5   = (long)in->plane5;
  out->plane6   = (long)in->plane6;
  out->plane7   = (long)in->plane7;
  out->jlocr    = (long)in->jlocr;
  out->jlocd    = (long)in->jlocd;
  out->jlocs    = (long)in->jlocs;
  out->jlocf    = (long)in->jlocf;
  out->jlocif   = (long)in->jlocif;
  if (in->do3D) out->do3D = 1;
  else          out->do3D = 0;
  out->xPxOff   = (float)in->xPxOff;
  out->yPxOff   = (float)in->yPxOff;
  return out;
} /* end ObitGPUSkyGeomImageH2D */

/** Public: Host to Device UV Descriptor conversion
    Create host resident version of a device UV Descriptor
*/
CUDAUVDesc* ObitGPUSkyGeomUVH2D (ObitUVDesc *in)
{
  CUDAUVDesc* out=NULL;
  int i, j;
  out = ObitGPUUVDescCreate ();
  out->nvis      = (long)in->nvis;
  out->lrec      = (long)in->lrec;
  out->nrparm    = (long)in->nrparm;
  out->firstVis  = (long)in->firstVis;
  out->numVisBuff= (long)in->numVisBuff;
  out->naxis     = (long)in->naxis;
  out->numSubA   = (long)in->numSubA;
  out->maxAnt    = (long)in->maxAnt;
  for (i=0; i<UV_MAXDIM; i++) out->inaxes[i] = (long)in->inaxes[i];
  for (i=0; i<UV_MAXDIM; i++) {
    for (j=0; j<UVLEN_KEYWORD; j++) out->ctype[i][j] = in->ctype[i][j];
  }
  for (i=0; i<UV_MAX_RANP; i++) {
    for (j=0; j<UVLEN_KEYWORD; j++) out->ptype[i][j] = in->ptype[i][j];
  }
  for (i=0; i<UV_MAXDIM; i++) out->cdelt[i] = (float)in->cdelt[i];
  for (i=0; i<UV_MAXDIM; i++) out->crpix[i] = (float)in->crpix[i];
  for (i=0; i<UV_MAXDIM; i++) out->crota[i] = (float)in->crota[i];
  for (i=0; i<UV_MAXDIM; i++) out->crval[i] = (double)in->crval[i];
  out->altCrpix = (float)in->altCrpix;
  out->xshift   = (float)in->xshift;
  out->yshift   = (float)in->yshift;
  out->beamMaj  = (float)in->beamMaj;
  out->beamMin  = (float)in->beamMin;
  out->beamPA   = (float)in->beamPA;
  out->DeltaTime= (float)in->DeltaTime;
  for (i=0; i<UVLEN_VALUE; i++) out->object[i]     = in->object[i];
  for (i=0; i<UVLEN_VALUE; i++) out->teles[i]      = in->teles[i];
  for (i=0; i<UVLEN_VALUE; i++) out->instrument[i] = in->instrument[i];
  for (i=0; i<UVLEN_VALUE; i++) out->observer[i]   = in->observer[i];
  for (i=0; i<UVLEN_VALUE; i++) out->obsdat[i]     = in->obsdat[i];
  for (i=0; i<UVLEN_VALUE; i++) out->date[i]       = in->date[i];
  for (i=0; i<UVLEN_VALUE; i++) out->origin[i]     = in->origin[i];
  for (i=0; i<UVLEN_VALUE; i++) out->bunit[i]      = in->bunit[i];
  out->epoch   = (float)in->epoch;
  out->equinox = (float)in->equinox;
  out->maxBL   = (float)in->maxBL;
  out->maxW    = (float)in->maxW;
  out->JDObs   = (double)in->JDObs;
  out->obsra   = (double)in->obsra;
  out->obsdec  = (double)in->obsdec;
  out->altRef  = (double)in->altRef;
  out->restFreq= (double)in->restFreq;
  out->freq    = (double)in->freq;
  out->VelReference   = (long)in->VelReference;
  out->VelDef  = (long)in->VelDef;
  out->ilocu   = (long)in->ilocu;
  out->ilocv   = (long)in->ilocv;
  out->ilocw   = (long)in->ilocw;
  out->iloct   = (long)in->iloct;
  out->ilocb   = (long)in->ilocb;
  out->iloca1  = (long)in->iloca1;
  out->iloca2  = (long)in->iloca2;
  out->ilocsa  = (long)in->ilocsa;
  out->ilocsu  = (long)in->ilocsu;
  out->ilocfq  = (long)in->ilocfq;
  out->ilocit  = (long)in->ilocit;
  out->ilocid  = (long)in->ilocid;
  out->ilocws  = (long)in->ilocws;
  out->jlocc   = (long)in->jlocc;
  out->jlocs   = (long)in->jlocs;
  out->jlocf   = (long)in->jlocf;
  out->jlocr   = (long)in->jlocr;
  out->jlocd   = (long)in->jlocd;
  out->jlocif  = (long)in->jlocif;
  out->incs    = (long)in->incs;
  out->incf    = (long)in->incf;
  out->incif   = (long)in->incif;
  out->ncorr   = (long)in->ncorr;
  out->isort[0] = in->isort[0]; out->isort[1] = in->isort[1]; out->isort[2] = in->isort[2]; 
  return out;
} /* end ObitGPUSkyGeomUVH2D */

#else  /* Not GPU */

/**
 * Creates a CUDAImageDesc
 * \return the new object in locked host memory
 */
CUDAImageDesc* ObitGPUImageDescCreate ()
{
  return NULL;
} /* end ObitGPUImageDescCreate */

/**
 * Destroys a CUDAImageDesc
 * Called from c
 * \param in   Object to delete
 */
void ObitGPUImageDescZap (CUDAImageDesc *in)
{
  return;
} /* end ObitGPUImageDescZap */

/**
 * Creates a CUDAUVDesc
 * \return the new object in locked host memory
 */
CUDAUVDesc* ObitGPUUVDescCreate ()
{
  return NULL;
} /* end ObitGPUUVDescCreate */

/**
 * Destroys a CUDAUVDesc
 * Called from c
 * \param in   Object to delete
 */
void ObitGPUUVDescZap (CUDAUVDesc *in)
{
  return;
} /* end ObitGPUUVDescZap */

/** Public: Host to Device Image Descriptor conversion
    Create host resident version of a device Image Descriptor
  * \param  in
*/
CUDAImageDesc* ObitGPUSkyGeomImageH2D (ObitImageDesc *in)
{
  return NULL;
} /* end ObitGPUSkyGeomImageH2D*/

/** Public: Host to Device UV Descriptor conversion
    Create host resident version of a device UV Descriptor
*/
CUDAUVDesc* ObitGPUSkyGeomUVH2D (ObitUVDesc *in)
{
  return NULL;
} /* end ObitGPUSkyGeomUVH2D */

#endif /* HAVE_GPU */
