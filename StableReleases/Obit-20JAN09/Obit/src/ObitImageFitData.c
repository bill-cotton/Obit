/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2006-2008                                          */
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

#include "ObitImageFitData.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitImageFitData.c
 * ObitImageFitData class function definitions.
 * This class is derived from the Obit base class.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitImageFitData";

/** Function to obtain parent ClassInfo */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo structure ObitImageFitDataClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitImageFitDataClassInfo myClassInfo = {FALSE};

/*--------------- File Global Variables  ----------------*/


/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitImageFitDataInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitImageFitDataClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitImageFitDataClassInfoDefFn (gpointer inClass);

/** Private: Evaluation function for solver. */
static void fxdvd (gpointer ddata, odouble *p, odouble *f, odouble *grad, olong iflag);
/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitImageFitData* newObitImageFitData (gchar* name)
{
  ObitImageFitData* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitImageFitDataClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitImageFitData));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitImageFitDataInit((gpointer)out);

 return out;
} /* end newObitImageFitData */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitImageFitDataGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitImageFitDataClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitImageFitDataGetClass */

/**
 * Creates an ObitImageFitData from an ObitFitRegion
 * \param name   An optional name for the object.
 * \param reg    ObitFitRegion to copy
 * \param bounds InfoList giving bounds
 * \li "PosGuard" OBIT_double (1,1,1) Distance from edge to allow center
 *                [def no bound]
 * \li "FluxLow"  OBIT_double (1,1,1) Lower bounds on Flux density [no bound]
 * \li "GMajUp"   OBIT_double (1,1,1) Major axis upper bound [no bound]
 * \li "GMajLow"  OBIT_double (1,1,1) Major axis lower bound [no bound]
 * \li "GMinUp"   OBIT_double (1,1,1) Minor axis upper bound [no bound]
 * \li "GMinLow"  OBIT_double (1,1,1) Minor axis lower bound [no bound]
 * \param image  Image with pixel data, (with selection used for reg)
 * \param err    Obit Error/mesage stack
 * \return the new object.
 */
ObitImageFitData* 
ObitImageFitDataCreate (gchar* name, ObitFitRegion *reg, 
			ObitInfoList *bounds, ObitImage *image, ObitErr *err)
{
  ObitImageFitData* out=NULL;
  gint32       dim[MAXINFOELEMDIM];
  ObitInfoType type;
  olong blc[2], trc[2],pos[2] ;
  ofloat cells, rmax, rmin, fblank=ObitMagicF();
  odouble dblank;
  gboolean doPosGuard, doFluxLow, doGMajUp, doGMajLow, doGMinUp, doGMinLow;
  odouble PosGuard, FluxLow, GMajUp, GMajLow, GMinUp, GMinLow;
  olong i, j, k, nparm, nvar;
  gchar *routine = "ObitImageFitDataCreate";

  /* error checks */
  if (err->error) return out;
  g_assert (ObitFitRegionIsA(reg));
  g_assert (ObitImageIsA(image));

  /* Bounds info */
  doPosGuard = ObitInfoListGetTest(bounds, "PosGuard", &type, dim, &PosGuard);
  doFluxLow  = ObitInfoListGetTest(bounds, "FluxLow",  &type, dim, &FluxLow);
  doGMajUp   = ObitInfoListGetTest(bounds, "GMajUp",   &type, dim, &GMajUp); 
  doGMajLow  = ObitInfoListGetTest(bounds, "GMajLow",  &type, dim, &GMajLow);
  doGMinUp   = ObitInfoListGetTest(bounds, "GMinUp",   &type, dim, &GMinUp);
  doGMinLow  = ObitInfoListGetTest(bounds, "GMinLow",  &type, dim, &GMinLow);

  /* Avoid crazy values */
  GMajLow = MAX(1.0e-20, GMajLow); doGMajLow = TRUE;
  GMinLow = MAX(1.0e-20, GMinLow); doGMinLow = TRUE;

  /* Create basic structure */
  out = newObitImageFitData (name);

  /* Read image if not already done */
  if (!image->image || (ObitFArraySum(image->image)==0.0)) {
    ObitImageOpen(image, OBIT_IO_ReadOnly,err);
    ObitImageRead (image, NULL, err);
    ObitImageClose (image, err);
    if (err->error) Obit_traceback_val (err, routine, image->name, out);
  } /* end read image */

  /* Pixel data */
  blc[0] = reg->corner[0]; 
  blc[1] = reg->corner[1];
  trc[0] = reg->corner[0] + reg->dim[0] - 1; 
  trc[1] = reg->corner[1] + reg->dim[1] - 1;
  out->pixels = ObitFArraySubArr (image->image, blc, trc, err);
  out->resids = ObitFArraySubArr (image->image, blc, trc, err);
  out->irms   = ObitFArrayRMS(image->image);

  /* Beam info */
  if (image->myDesc->beamMaj>0.0) {
    out->beam[0] = image->myDesc->beamMaj/fabs(image->myDesc->cdelt[1]);
    out->beam[1] = image->myDesc->beamMin/fabs(image->myDesc->cdelt[1]);
    out->beam[2] = image->myDesc->beamPA*DG2RAD;
  } else {
    out->beam[0] = 4.0;
    out->beam[1] = 4.0;
    out->beam[2] = 0.0;
  }


  /* Check that something to fit */
  Obit_retval_if_fail((ObitFArrayRMS(out->pixels)>0.0), err, out,
		       "%s: Nothing to fit, RMS=0",routine);

  /* Restoring beam area in pixels */
  cells = sqrt (fabs(image->myDesc->cdelt[0]*image->myDesc->cdelt[1]));
  if (image->myDesc->beamMaj>0.0) {
    out-> beamarea = 1.1331 * (image->myDesc->beamMaj/cells) * 
      (image->myDesc->beamMin/cells);
  } else out->beamarea = 1.0;

  /* Fiddling to get fitting to work better
     first scale such that peak value is 5.0 */
  rmax = ObitFArrayMax (out->pixels, pos);
  rmin = ObitFArrayMin (out->pixels, pos);
  out->rscale = 5.0 / rmax;
  if (out->rscale==0.0) out->rscale = 1.0;
  ObitFArraySMul(out->pixels, out->rscale);

  /* count things in reg - for now fit all - only do Gaussians */
  nparm = 0;
  nvar = 0;
  for (i=0; i<reg->nmodel; i++) {
    if (reg->models[i]->type == OBIT_FitModel_GaussMod) nparm += 6;
    if (reg->models[i]->type == OBIT_FitModel_GaussMod) nvar += 6;
  }

  /* Check that something to fit */
  Obit_retval_if_fail((nparm>0), err, out, "%s: Nothing to fit, nparm=0",routine);

  /* Set values */
  dblank  = (odouble)fblank;
  out->fx = (ObitImageFitDataFuncFP) fxdvd;  /* evaluation function */
  out->ncomp  = reg->nmodel;
  out->nparm  = nparm;
  out->nvar   = nvar;

  /* Allocate arrays */
  out->np    = g_malloc0(out->ncomp*sizeof(olong));
  out->type  = g_malloc0(out->ncomp*sizeof(olong));
  out->ptemp = g_malloc0(out->nparm*sizeof(odouble));
  out->ivar  = g_malloc0(out->nparm*sizeof(olong));
  out->jvar  = g_malloc0(out->nparm*sizeof(olong));
  out->p     = g_malloc0(out->ncomp*sizeof(odouble*));
  out->dp    = g_malloc0(out->ncomp*sizeof(odouble*));
  out->pf    = g_malloc0(out->ncomp*sizeof(odouble*));
  out->pu    = g_malloc0(out->ncomp*sizeof(odouble*));
  out->pl    = g_malloc0(out->ncomp*sizeof(odouble*));
  out->e     = g_malloc0(out->ncomp*sizeof(odouble*));
  for (i=0; i<out->ncomp; i++) {
    if (reg->models[i]->type == OBIT_FitModel_GaussMod) out->np[i] = 6;
    else out->np[i] = 1;
    out->p[i]  = g_malloc0(out->np[i]*sizeof(odouble));
    out->dp[i] = g_malloc0(out->np[i]*sizeof(odouble));
    out->pf[i] = g_malloc0(out->np[i]*sizeof(odouble));
    out->pu[i] = g_malloc0(out->np[i]*sizeof(odouble));
    out->pl[i] = g_malloc0(out->np[i]*sizeof(odouble));
    out->e[i]  = g_malloc0(out->np[i]*sizeof(odouble));
  }

  /* Initialize */
  k = 0;
  for (i=0; i<out->ncomp; i++) {
    out->type[i] = reg->models[i]->type;
    for (j=0; j<out->np[i]; j++) {
      out->ivar[k] = i;  /* keeps track of fitted variables */
      out->jvar[k] = j;
      out->p[i][j] = 0.0;
      out->dp[i][j] = 0.001;
      out->pf[i][j] = TRUE;
      out->pu[i][j] = dblank;
      out->pl[i][j] = dblank;
      out->e[i][j]  = 1.0;
      k++;
    }
  } /* end loop over components */

  /* Fill Values - keep away from lower bound */
  for (i=0; i<out->ncomp; i++) {
    if (out->type[i] == OBIT_FitModel_GaussMod) {
      out->p[i][0]  = out->rscale * reg->models[i]->Peak;
      out->p[i][1]  = reg->models[i]->DeltaX;
      out->p[i][2]  = reg->models[i]->DeltaY;
      /* internal and external different maj, min */
      out->p[i][4]  = reg->models[i]->parms[0];
      out->p[i][3]  = reg->models[i]->parms[1]; /* break degeneracy */
      out->p[i][5]  = reg->models[i]->parms[2];
      out->e[i][0]  = reg->models[i]->ePeak;
      out->e[i][1]  = reg->models[i]->eDeltaX;
      out->e[i][2]  = reg->models[i]->eDeltaY;
      out->e[i][4]  = reg->models[i]->eparms[0];
      out->e[i][3]  = reg->models[i]->eparms[1];
      out->e[i][5]  = reg->models[i]->eparms[2];
      /* solution bounds */
      if (doFluxLow) {
	out->pl[i][0] = out->rscale * FluxLow;
	out->p[i][0]  = MAX (out->p[i][0], out->pl[i][0]+out->dp[i][0]);
      }
      if (doPosGuard) {
	out->pl[i][1]  = PosGuard;
	out->pl[i][2]  = PosGuard;
	out->pu[i][1]  = reg->dim[0] - 1 - PosGuard;
	out->pu[i][2]  = reg->dim[1] - 1 - PosGuard;
      }
      if (doGMajUp)  out->pu[i][3] = GMajUp;
      if (doGMajLow) {
	out->pl[i][3] = GMajLow;
	out->p[i][3]  = MAX (out->p[i][3], out->pl[i][3]+out->dp[i][3]);
      }
      if (doGMinUp)  out->pu[i][4] = GMinUp;
      if (doGMinLow) {
	out->pl[i][4] = GMinLow;
 	out->p[i][4]  = MAX (out->p[i][4], out->pl[i][4]+out->dp[i][4]);
     }
      /* Break possible axis degeneracy */
      if (fabs(out->p[i][3]-out->p[i][4])<0.1*out->p[i][3]) {
	out->p[i][3]  *= 1.1;
      }
    } /* end if Gaussian */
  } /* end loop over components */

  return out;
} /* end ObitImageFitDataCreate */
/**
 * copy values from an ObitImageFitData to an ObitFitRegion
 * \param data   Fit data to copy
 * \param reg    ObitFitRegion to update
 */
void ObitImageFitData2Reg (ObitImageFitData* data, ObitFitRegion *reg)
{
  olong i, pos[2];
  ofloat rcorr, corr, temp;

  /* Correction to errors for number of data points */
  rcorr = 1.0 / data->rscale;
  corr  = 1.0 / ((ofloat)data->pixels->naxis[0] * (ofloat)data->pixels->naxis[1]);
  corr  = rcorr*sqrt(corr);

  /* Fill Values */
  for (i=0; i<data->ncomp; i++) {
    if (data->type[i] == OBIT_FitModel_GaussMod) {
      reg->models[i]->Peak      = data->p[i][0] * rcorr;
      reg->models[i]->DeltaX    = data->p[i][1];
      reg->models[i]->DeltaY    = data->p[i][2];
      /* internal and external different maj, min */
      reg->models[i]->parms[0]  = data->p[i][4];
      reg->models[i]->parms[1]  = data->p[i][3];
      reg->models[i]->parms[2]  = data->p[i][5];
      /* Correct major and minor axes */
      if (reg->models[i]->parms[0]<reg->models[i]->parms[1]) {
	temp = reg->models[i]->parms[1];
	reg->models[i]->parms[1] = reg->models[i]->parms[0];
	reg->models[i]->parms[0] = temp;
	reg->models[i]->parms[2] -= 90.0*DG2RAD;
      }
      /* Errors */
      ObitImageFitDataGaussErr (reg->models[i]->Peak, reg->models[i]->parms[0], reg->models[i]->parms[1],
				reg->models[i]->parms[2], data->irms, data->beam,
				&reg->models[i]->ePeak, &reg->models[i]->eDeltaX, &reg->models[i]->eDeltaY,
				&reg->models[i]->eparms[0], &reg->models[i]->eparms[1],
				&reg->models[i]->eparms[2]);
    }
  } /* end loop over components */

  /* stats */
  reg->peakResid = rcorr*ObitFArrayMaxAbs (data->resids, pos);
  reg->RMSResid  = rcorr*ObitFArrayRMS0 (data->resids);
  reg->fluxResid = rcorr*ObitFArraySum (data->resids)/data->beamarea;
} /* end ObitImageFitData2Reg */


/**
 * Routine to determine errors in Gaussian parameters
 * Routine addapted from the AIPSish CORERR.FOR
 * \param peak      Peak Ipol (Jy/beam) 
 * \param major     Fitted major axis size (pixel)
 * \param minor     Fitted minor axis size  (pixel)
 * \param posang    Fitted PA (rad)
 * \param irms      RMS (sigma) in Ipol. 
 * \param beam      Restoring beam major, minor axes and position angle 
 *                  (pixel, pixel, rad) 
 * \param epeak     [out] Error in peak (Jy/beam)
 * \param errra     [out] Error of Right Ascension (pixel)
 * \param errdec    [out] Error of Declination (pixel)
 * \param errmaj    [out] Error of major axis size (pixel)
 * \param errmin    [out] Error of Minor axis size (pixel)
 * \param errpa     [out] Error of position angle (rad)
 */
void ObitImageFitDataGaussErr (ofloat peak, ofloat major, ofloat minor, ofloat posang, 
			       ofloat irms, ofloat* beam, 
			       ofloat* epeak, ofloat* errra, ofloat* errdec, 
			       ofloat* errmaj, ofloat* errmin, ofloat* errpa)
{
  ofloat bemrat, snr, snramp, snrmaj, snrmin, sinc, cosc,
    tmaj, tmin, errx2, erry2, psf[3];

  psf[0] = beam[0];
  psf[1] = beam[1];
  psf[2] = beam[2];

  /* Beam ratio */
  bemrat = (major/psf[0]) * (minor/psf[1]);

  /* Trap under size beams */
  bemrat = MAX (bemrat, 1.0);

  /* Effective SNRs^2 to account for correlated noise. */
  snr = peak / irms;

  /* SNR**2 for amplitude errors */
  snramp = ((major*minor)/(4.0*psf[0]*psf[0])) * 
    (pow(1.0+pow(psf[0]/major,2.0),1.5)) * 
    (pow(1.0+pow(psf[1]/minor,2.0),1.5)) * snr * snr;

  /* SNR**2 for major axis error */
  snrmaj = ((major*minor)/(4.0*psf[0]*psf[0])) * 
    (pow(1.0+pow(psf[0]/major,2.0),2.5)) * 
    (pow(1.0+pow(psf[1]/minor,2.0),0.5)) * snr * snr;

  /* SNR**2 for minor axis/PA errors */
  snrmin = ((major*minor)/(4.0*psf[0]*psf[0])) * 
    (pow(1.0+pow(psf[0]/major,2.0),0.5)) * 
    (pow(1.0+pow(psf[1]/minor,2.0),2.5)) * snr * snr;

  /* Flux errors */
  (*epeak) = sqrt ((2.0 * peak*peak / snramp));

  /* Errors */
  sinc = sin (posang);
  cosc = cos (posang); 

  /* Trap under size beams */
  tmaj = MAX (major, psf[1]);
  tmin = MAX (minor, psf[1]);

  /* Axis sizes include */
  (*errmaj) = sqrt (((2.0 * tmaj*tmaj) / snrmaj));
  (*errmin) = sqrt (((2.0 * tmin*tmin) / snrmin));
  (*errpa)  = sqrt (pow(tmaj*tmin/(tmaj*tmaj-tmin*tmin+1.0e-20),2.0) 
		    * 4.0 / snrmin);
  (*errpa) = MIN ((*errpa), G_PI*0.5);

  /* Position errors */
  errx2 = 2.0 * (2.0 * tmaj*tmaj) / (8.0 * log (2.0) * snrmaj);
  erry2 = 2.0 * (2.0 * tmin*tmin) / (8.0 * log (2.0) * snrmin);

  (*errra)  = sqrt (errx2*sinc*sinc + erry2*cosc*cosc);
  (*errdec) = sqrt (erry2*sinc*sinc + errx2*cosc*cosc);
} /* end of routine ObitImageFitDataGaussErr */ 

/**
 * Initialize global ClassInfo Structure.
 */
void ObitImageFitDataClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitImageFitDataClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitImageFitDataClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitImageFitDataClassInfoDefFn (gpointer inClass)
{
  ObitImageFitDataClassInfo *theClass = (ObitImageFitDataClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitImageFitDataClassInit;
  theClass->newObit       = (newObitFP)newObitImageFitData;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitImageFitDataClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitImageFitDataGetClass;
  theClass->ObitCopy      = NULL;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitImageFitDataClear;
  theClass->ObitInit      = (ObitInitFP)ObitImageFitDataInit;
  theClass->ObitImageFitDataCreate = (ObitImageFitDataCreateFP)ObitImageFitDataCreate;

} /* end ObitImageFitDataClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitImageFitDataInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitImageFitData *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->pixels = NULL;
  in->resids = NULL;
  in->ivar   = NULL;
  in->jvar   = NULL;
  in->ptemp  = NULL;
  in->np     = NULL;
  in->p      = NULL;
  in->dp     = NULL;
  in->pf     = NULL;
  in->pu     = NULL;
  in->pl     = NULL;
  in->e      = NULL;

} /* end ObitImageFitDataInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitImageFitData* cast to an Obit*.
 */
void ObitImageFitDataClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  olong i;
  ObitImageFitData *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  in->pixels = ObitFArrayUnref(in->pixels);
  in->resids = ObitFArrayUnref(in->resids);
  if (in->ivar) g_free(in->ivar);
  if (in->jvar) g_free(in->jvar);
  if (in->ptemp) g_free(in->ptemp);
  if (in->np) g_free(in->np);
  if (in->p) {
    for (i=0; i<in->ncomp; i++) if (in->p[i]) g_free(in->p[i]);
    if (in->p) g_free(in->p);
  }
  if (in->dp) {
    for (i=0; i<in->ncomp; i++) if (in->dp[i]) g_free(in->dp[i]);
    if (in->dp) g_free(in->dp);
  }
  if (in->pf) {
    for (i=0; i<in->ncomp; i++) if (in->pf[i]) g_free(in->pf[i]);
    if (in->pf) g_free(in->pf);
  }
  if (in->pu) {
    for (i=0; i<in->ncomp; i++) if (in->pu[i]) g_free(in->pu[i]);
    if (in->pu) g_free(in->pu);
  }
  if (in->pl) {
    for (i=0; i<in->ncomp; i++) if (in->pl[i]) g_free(in->pl[i]);
    if (in->pl) g_free(in->pl);
  }
  if (in->e) {
    for (i=0; i<in->ncomp; i++) if (in->e[i]) g_free(in->e[i]);
    if (in->e) g_free(in->e);
  }
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitImageFitDataClear */

/**
 * Given the vector P of solution parameters, this routine computes  
 * the value of the chi-squared function F (a sum of squared residuals),  
 * and, optionally, the gradient, GRAD, of F w.r.t. P.  When IFLAG=1,  
 * only F is computed.  Otherwise F and GRAD are both computed.  Note  
 * that P is to contain only the parameters which are being solved for  
 * --- not the parameters that are to be held fixed.  This subroutine is  
 * called by the minimization routine DVDMIN.  
 * Only OBIT_FitModel_GaussMod models currently handled.
 *  
 * Additionally, the residuals (model minus data) are stored in the  
 * data object for use outside the minimization  routine proper.  
 * (The minimization routine DVDMIN doesn't need to know the residuals, 
 * it only needs F and GRAD).  
 * Routine adopted from the AIPSish VSAD.FOR/FXDVD  
 *  
 * \param data    Data object
 * \param p       Vector of least-squares solution pararameters. 
 * \param f       [out] The value of the average chi-squared function corresponding 
 *                to the given P. 
 * \param grad    [out] The gradient of the chi-squared function.  I.e., 
 *                grad[i] = derivative of f w.r.t. p[i]. 
 * \param iflag   iflag=1 ==> compute just F, 
 *                iflag !=1 ==> compute both F and GRAD. 
 */
static void fxdvd (gpointer ddata, odouble *p, odouble *f, odouble* grad, olong iflag) 
{
  olong   nk=0, npts, i, j, k, ix1, iy1, ntot, l;
  ofloat sth2, cth2, s2th, c2th, mj, mn, va, vb, vc, vd, x, y, 
    x2, y2, xy, tworfv, twocon, fv, g4c, g5c, tmp, ex;
  olong pos[2];
  olong nvar, ncomp, nparm, nx, ny;
  ofloat fblank=ObitMagicF();
  ofloat *pixels=NULL, *resid=NULL;
  ofloat con = 2.772589;
  odouble *ptemp, **lp, chisq, dblank, arg;
  /* odouble dbug[6]; DEBUG */
  ObitImageFitData *data =  (ObitImageFitData*)ddata;

  /* Test data for type */
  g_assert (ObitImageFitDataIsA(data));

  /* pointers to pixel data, residuals */
  pos[0] = pos[1] = 0;
  pixels = ObitFArrayIndex (data->pixels, pos);
  resid  = ObitFArrayIndex (data->resids, pos);
  dblank = (odouble)fblank;

  /* Local values */
  nvar  = data->nvar;
  ncomp = data->ncomp;
  nparm = data->nparm;
  nx = data->pixels->naxis[0];
  ny = data->pixels->naxis[1];
  ptemp = data->ptemp;
  lp    = data->p;  /* local pointer to data version of parameters */

  twocon = 2.0*con;
  npts = nx*ny;
  for (k=0; k<npts; k++) { /* loop 10 */
    if (pixels[k] != fblank) resid[k] = -pixels[k];
    else pixels[k] = fblank;
  } /* end loop  L10:  */

  ntot = nk;
  for (i=0; i<nvar; i++) { /* loop 20 */
    /* Enforce hard limits */
    if (data->pl[data->ivar[i]][data->jvar[i]]!=dblank)
      p[i] = MAX (p[i], data->pl[data->ivar[i]][data->jvar[i]]);
    if (data->pu[data->ivar[i]][data->jvar[i]]!=dblank)
      p[i] = MIN (p[i], data->pu[data->ivar[i]][data->jvar[i]]);
    lp[data->ivar[i]][data->jvar[i]] = p[i];
  } /* end loop  L20:  */

  /*  For the Ith Gaussian component,
      p[i][0] = the peak amplitude of the component, 
      p[i][1] = x-position (pixels), 
      p[i][2] = y-position (pixels), 
      p[i][3] = major axis FWHM, 
      p[i][4] = minor axis FWHM, 
      and, p[i][5] = position angle of the major axis, normally measured from 
                     North through East. (rad) */

  ix1 = 1;
  iy1 = 1;
  k = 0;
  chisq = 0.0;
  for (i=0; i<ncomp; i++) { /* loop 80 */
    if (data->type[i] == OBIT_FitModel_GaussMod) {
      tmp  = sin (lp[i][5]);
      sth2 = tmp*tmp;
      tmp  = cos (lp[i][5]);
      cth2 = tmp*tmp;
      s2th = -sin (2.0*lp[i][5]);
      c2th =  cos (2.0*lp[i][5]);
      mj = lp[i][3]*lp[i][3]/con;
      mn = lp[i][4]*lp[i][4]/con;
      va = cth2/mj+sth2/mn;
      vb = sth2/mj+cth2/mn;
      vc = s2th*((1.0/mn)-(1.0/mj));
      for (l=0; l<npts; l++) { /* loop 70 */
	if (pixels[l] != fblank) {
	  x = ix1 + (l % nx) - lp[i][1];
	  y = iy1 + (olong)(l/nx) - lp[i][2];
	  x = ((va*x + vc*y)*x + vb*y*y);
	  /* limit accuracy to 10**-4 to save time */
	  if (x < 9.2) {
	    fv = lp[i][0]*exp(-x);
	  } else {
	    fv = 0.0;
	  } 
	  resid[l] += fv;
	}
      } /* end loop  L70:  */
      k += data->np[i];
    } /* end if Gaussian*/
  } /* end loop  L80:  */

  /* Get Chi squared */
  chisq = 0.0;
  for (k=0; k<npts; k++) { /* loop 10 */
    if (resid[k] != fblank) chisq += resid[k]*resid[k];
  } /* end loop  L10:  */

  /* Returned objective f is Chi squared*/
  *f = chisq;

  /* get gradient if requested */
  if (iflag != 1) {
    for (i=0; i<nparm; i++) { /* loop 90 */
      ptemp[i] = 0.0;
    } /* end loop  L90:  */
    
    k = 0;
    for (i=0; i<ncomp; i++) { /* loop 140 */
      if (data->type[i] == OBIT_FitModel_GaussMod) {
	tmp  = sin (lp[i][5]);
	sth2 = tmp*tmp;
	tmp  = cos (lp[i][5]);
	cth2 = tmp*tmp;
	s2th = -sin (2.0*lp[i][5]);
	c2th =  cos (2.0*lp[i][5]);
	mj = lp[i][3]*lp[i][3]/con;
	mn = lp[i][4]*lp[i][4]/con;
	g4c = twocon/(lp[i][3]*lp[i][3]*lp[i][3]);
	g5c = twocon/(lp[i][4]*lp[i][4]*lp[i][4]);
	va = cth2/mj + sth2/mn;
	vb = sth2/mj + cth2/mn;
	vc = s2th*(1.0/mn - 1.0/mj);
	vd = c2th*(1.0/mn - 1.0/mj);
	for (l=0; l<npts; l++) { /* loop 130 */
	  if (pixels[l] != fblank) {
	    x = ix1 + (l % nx) - lp[i][1];
	    y = iy1 + (olong)(l/nx) - lp[i][2];
	    x2 = x*x;
	    y2 = y*y;
	    xy = x*y;
	    arg = va*x2 + vb*y2 + vc*xy;
	    if (arg < 9.2) ex = exp(-(arg));
	    else ex = 0.0;
	    fv = lp[i][0] * ex;
	    tworfv = 2.0 * resid[l] * fv;
	    if (data->pf[i][0]) ptemp[k]   += 2.0 * resid[l] * ex;
	    if (data->pf[i][1]) ptemp[k+1] += tworfv * (2.0*x*va + y*vc);
	    if (data->pf[i][2]) ptemp[k+2] += tworfv * (2.0*y*vb + x*vc);
	    if (data->pf[i][3]) ptemp[k+3] += tworfv * g4c*(x2*cth2 + y2*sth2 - xy*s2th);
	    if (data->pf[i][4]) ptemp[k+4] += tworfv * g5c*(x2*sth2 + y2*cth2 + xy*s2th);
	    if (data->pf[i][5]) ptemp[k+5] += tworfv * (vc*(x2-y2) + 2.0*vd*xy);
	  } 
	} /* end loop  L130: */
	k += data->np[i];
      } /* end if Gaussian */
    } /* end loop  L140: */
    
    /* copy gradient terms for fitted parameters */
    k = 0;
    l = 0;
    for (i=0; i<ncomp; i++) { /* loop 160 */
      if (data->type[i] == OBIT_FitModel_GaussMod) {
	for (j=0; j<data->np[i]; j++) { /* loop 150 */
	  if (data->pf[i][j]) {
	    grad[k] = ptemp[l];
	    k++;
	  } 
	  l++;
	} /* end loop  L150: */
      } /* end if Gaussian */
    } /* end loop  L160: */
  } /* end of if need gradients */
  
  /* DEBUG
  fprintf (stdout, " ChiSq %lg penalty  %lg\n", chisq, penalty);
  fprintf (stdout, " p= %lf %lf %lf %lf %lf %lf \n", 
	   lp[0][0], lp[0][1], lp[0][2], lp[0][3], lp[0][4], lp[0][5]);
  if (iflag != 1) {
    fprintf (stdout, " grad= %lf %lf %lf %lf %lf %lf\n", 
	     grad[0], grad[1], grad[2], grad[3], grad[4], grad[5]);
    fprintf (stdout, " gpen= %lf %lf %lf %lf %lf %lf\n", 
	     dbug[0], dbug[1], dbug[2], dbug[3], dbug[4], dbug[5]);
  } */
  return;
} /* end of routine fxdvd */ 
