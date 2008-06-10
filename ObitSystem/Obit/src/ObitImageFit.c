/* $Id$        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2006                                               */
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
#include "gsl/gsl_cblas.h"
#include "ObitImageFit.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitImageFit.c
 * ObitImageFit class function definitions.
 * This class is derived from the Obit base class.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitImageFit";

/** Function to obtain parent ClassInfo */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo structure ObitImageFitClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitImageFitClassInfo myClassInfo = {FALSE};

/*--------------- File Global Variables  ----------------*/


/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitImageFitInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitImageFitClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitImageFitClassInfoDefFn (gpointer inClass);

/** Private: Do fitting. */
static void 
dvdmin (ObitImageFitDataFuncFP fx, gpointer data,  odouble* xi, odouble *xerr, 
	gint n, odouble eps, olong itmax, odouble *fopt, odouble *gnopt, 
	gint *ier, olong npr, ObitErr *err);
static odouble dmachx(gint job);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitImageFit* newObitImageFit (gchar* name)
{
  ObitImageFit* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitImageFitClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitImageFit));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitImageFitInit((gpointer)out);

 return out;
} /* end newObitImageFit */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitImageFitGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitImageFitClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitImageFitGetClass */

/**
 * Make a deep copy of an ObitImageFit.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitImageFit* ObitImageFitCopy  (ObitImageFit *in, ObitImageFit *out, 
				 ObitErr *err)
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
    out = newObitImageFit(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  out->image = ObitImageUnref(out->image);
  out->image = ObitImageRef(in->image);
  out->reg   = ObitFitRegionUnref(out->reg);
  out->reg   = ObitFitRegionRef(in->reg);
  out->data  = ObitImageFitDataUnref(out->data);
  out->data  = ObitImageFitDataRef(in->data);

  return out;
} /* end ObitImageFitCopy */

/**
 * Make a copy of a object but do not copy the actual data
 * This is useful to create an ImageFit similar to the input one.
 * \param in  The object to copy
 * \param out An existing object pointer for output, must be defined.
 * \param err Obit error stack object.
 */
void ObitImageFitClone  (ObitImageFit *in, ObitImageFit *out, ObitErr *err)
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
  out->image = ObitImageUnref(out->image);
  out->image = ObitImageRef(in->image);
  out->reg   = ObitFitRegionUnref(out->reg);
  out->reg   = ObitFitRegionRef(in->reg);
  out->data  = ObitImageFitDataUnref(out->data);
  out->data  = ObitImageFitDataRef(in->data);

} /* end ObitImageFitClone */

/**
 * Creates an ObitImageFit 
 * \param name   An optional name for the object.
 * \return the new object.
 */
ObitImageFit* ObitImageFitCreate (gchar* name)
{
  ObitImageFit* out;

  /* Create basic structure */
  out = newObitImageFit (name);

  return out;
} /* end ObitImageFitCreate */

/**
 * Does fit
 * \param in    Fitter to use, control info in info:
 * \li "MaxIter"  OBIT_int (1,1,1) Maximum number of iterations
 *                defaults to 10 per fitted parameter
 * \li "prtLv"    OBIT_int (1,1,1) Message level, 0=>none [def 0]
 * \li "PosGuard" OBIT_double (1,1,1) Distance from edge to allow center
 *                [def no bound]
 * \li "FluxLow"  OBIT_double (1,1,1) Lower bounds on Flux density [no bound]
 * \li "GMajUp"   OBIT_double (1,1,1) Major axis upper bound [no bound]
 * \li "GMajLow"  OBIT_double (1,1,1) Major axis lower bound [no bound]
 * \li "GMinUp"   OBIT_double (1,1,1) Minor axis upper bound [no bound]
 * \li "GMinLow"  OBIT_double (1,1,1) Minor axis lower bound [no bound]
 * \param image Image to fit, should have pixel array attached and
 *              opened as described in reg (BLC, TRC).
 * \param reg   Region in image to fit, on input the initial guess
 *              on output, the fitted result
 * \param err   Obit error/message stack.
 * \return a completion code: 0=> converged, 1=> max. iterations.
 */
olong ObitImageFitFit (ObitImageFit* in,  ObitImage *image, 
		      ObitFitRegion* reg, ObitErr *err)
{
  olong i, j, itmax, npr, nvar, fst, ier = -1;
  odouble eps, fopt, gnopt, *xi=NULL, *xerr=NULL;
  odouble dblank;
  ofloat fblank=ObitMagicF();
  gint32       dim[MAXINFOELEMDIM];
  ObitInfoType type;
  gchar *routine = "ObitImageFitFit";

  /* error checks */
  if (err->error) return ier;
  g_assert (ObitImageFitIsA(in));
  g_assert (ObitImageIsA(image));
  g_assert (ObitFitRegionIsA(reg));

  /* Save input */
  in->image = ObitImageUnref(in->image);
  in->image = ObitImageRef(image);
  in->reg   = ObitFitRegionUnref(in->reg);
  in->reg   = ObitFitRegionRef(reg);

  /* Create ObitImageFitData (delete old) 
   this copies info */
  in->data = ObitImageFitDataUnref(in->data);
  in->data = ObitImageFitDataCreate ("data", reg, in->info, image, err);
  if (err->error) goto cleanup;

  /* Setup for fitting */
  dblank = (odouble)fblank;
  eps  = 1.0e-8;
  xi   = g_malloc0 (in->data->nparm*sizeof(odouble));
  xerr = g_malloc0 (in->data->nparm*sizeof(odouble));
  nvar = 0;
  for (i=0; i<in->data->ncomp; i++) {
    if (in->data->type[i] == OBIT_FitModel_GaussMod) {
      fst = nvar;
      for (j=0; j<in->data->np[i]; j++) {
	if (in->data->pf[i][j]) {  /* Fitting this one? */
	  xi[nvar]   = in->data->p[i][j];
	  /* Make sure in range */
	  if (in->data->pl[i][j]!=dblank) 
	    xi[nvar] = MAX (xi[nvar], in->data->pl[i][j]);
	  if (in->data->pu[i][j]!=dblank) 
	    xi[nvar] = MIN (xi[nvar], in->data->pu[i][j]);
	  xerr[nvar] = fabs(in->data->e[i][j]);
	  nvar++;
	} /* end if fitting */
      } /* end loop over parameters */
   } /* end if Gaussian */
  } /* end loop over components */

  /* Control info */
  npr = 0;
  ObitInfoListGetTest(in->info, "prtLv", &type, dim, &npr);
  itmax = 10 * nvar;
  ObitInfoListGetTest(in->info, "MaxIter", &type, dim, &itmax);

  /* Do fit */
  dvdmin (in->data->fx, in->data, xi, xerr, nvar, 
	  eps, itmax, &fopt, &gnopt, &ier, npr, err);
  if (err->error) goto cleanup;


  /* Save errors */
  nvar = 0;
  for (i=0; i<in->data->ncomp; i++) {
    if (in->data->type[i] == OBIT_FitModel_GaussMod) {
      for (j=0; j<in->data->np[i]; j++) {
	if (in->data->pf[i][j]) {  /* Fitting this one? */
	  in->data->e[i][j] = xerr[nvar];
	  nvar++;
	} /* end if fitting */
      } /* end loop over parameters */
    } /* end if Gaussian */
  }

  /* Copy results back to region */
  ObitImageFitData2Reg (in->data, reg);

  /* RMS residual if requested */
  if (npr>0) {
    fopt = ObitFArrayRMS0(in->data->resids);
    fopt /= in->data->rscale;
    Obit_log_error(err, OBIT_InfoErr, "RMS residual = %g", fopt);
    ObitErrLog(err); 
 }

  /* Cleanup */
 cleanup:
  if (xi)   g_free(xi);
  if (xerr) g_free(xerr);
  if (err->error) Obit_traceback_val (err, routine, in->name, ier);
  return ier;
} /* end ObitImageFitFit  */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitImageFitClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitImageFitClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitImageFitClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitImageFitClassInfoDefFn (gpointer inClass)
{
  ObitImageFitClassInfo *theClass = (ObitImageFitClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitImageFitClassInit;
  theClass->newObit       = (newObitFP)newObitImageFit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitImageFitClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitImageFitGetClass;
  theClass->ObitCopy      = (ObitCopyFP)ObitImageFitCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitImageFitClear;
  theClass->ObitInit      = (ObitInitFP)ObitImageFitInit;
  theClass->ObitImageFitCreate = (ObitImageFitCreateFP)ObitImageFitCreate;
  theClass->ObitImageFitFit    = (ObitImageFitFitFP)ObitImageFitFit;

} /* end ObitImageFitClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitImageFitInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitImageFit *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->image = NULL;
  in->reg   = NULL;
  in->data  = NULL;
  in->info  = newObitInfoList();

} /* end ObitImageFitInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitImageFit* cast to an Obit*.
 */
void ObitImageFitClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitImageFit *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  if (in->info) ObitInfoListUnref (in->info); in->info = NULL;
  in->image = ObitUnref(in->image);
  in->reg   = ObitUnref(in->reg);
  in->data  = ObitUnref(in->data);
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitImageFitClear */

/**
 * This is a translated from a Fortran implementation of Davidon's 
 * optimally conditioned  variable metric (quasi-Newton) method for 
 * function minimization.  It  is based on the algorithm given in 
 * W. C. Davidon:  Optimally conditioned optimization algorithms 
 * without line searches,   Mathematical Programming, vol. 9 (1975) 
 * pp. 1-30.  One should refer   to that reference for the algorithmic 
 * details.  
 *   Here, the steps of   the algorithm which are delineated by COMMENT 
 * lines correspond to  
 * the numbered steps in Davidon's paper.  The user must supply a  
 * subroutine fx to calculate the objective function and its gradient  
 * at a given point.  The objective function F is assumed to be a  
 * real-valued function of N real variables.  Here, 0 is assumed to be  
 * a lower bound for F.  If F can assume negative values, Step 2 must  
 * be modified in one of two different ways, depending on whether a  
 * lower bound is known (see Davidon for details).  
 * Routine translated from the AIPSish VSAD.FOR/DVDMIN  
 *  Remarks: 
 *  0) This works best when all the parameters are of order 1 
 *     (or at least the same order)
 *  1) This algorithm can be used for under-determined problems. 
 *  2) It maintains an approximation, in factored form J*transpose(J), 
 *     to the inverse Hessian of F.  At each iteration, a rank two update 
 *     is added to this approximation.  This approximation remains posi- 
 *     tive definite throughout the iteration.  In cases where an un- 
 *     known, say the Ith unknown, is ill-determined, ERR(I) will be 
 *     finite on exit from this routine. So, in least-squares applica- 
 *     tions, the error estimates for ill-determined parameters are like- 
 *     ly to be too small. 
 *  2.5) In the case of an under-determined problem (i.e., when the 
 *     Hessian matrix is singular) J*transpose(J) is a non-singular
 *     matrix whose inverse is close to the Hessian matrix.
 *  3) Furthermore, in cases where an excellent initial guess is supplied
 *     by the user, DVDMIN is likely to converge before it has iterated
 *     long enough to get a good approximation to the inverse Hessian.
 *     (Understand that it is trying to estimate this second-order in-
 *     formation only from the first-order information that is supplied
 *     by FX.)  So, in least-squares applications, when convergence oc-
 *     curs in just a couple of iterations, the derived error estimates
 *     may be inaccurate.
 *  4) Another Fortran implementation is given in the technical report
 *     by W. C. Davidon and L. Nazareth:  DRVOCR - A Fortran implementa-
 *     tion of Davidon's optimally conditioned method, Argonne National
 *     Lab., Applied Math. Div. Technical Memo. No. 306, August 1977.
 *  5) Comparisons of Davidon's algorithm with other quasi-Newton mini-
 *     mization algorithms are given in  J. N. Lyness:  A bench mark
 *     experiment for minimization algorithms, Math. of Computation,
 *     vol. 33 (1979) pp. 249-264.  This algorithm compares quite favor-
 *     ably with others, including the routine QNMDER of Gill et al.,
 *     and the Harwell Library routine VA13AD.
 *  6) Argonne Lab.'s MINPACK routines (non-proprietary) or NAG Library
 *     routines (proprietary) could be used in place of DVDMIN.  They
 *     would provide somewhat more flexibility.  They're a bit more con-
 *     servative (and therefore more robust, but perhaps less efficient).
 * \param fx         A user-supplied function of the type DVDFuncFP:
 *                   fx (data, x, f, g, k) which is used to calculate the 
 *                   value of the objective function f at x and, optionally, 
 *                   the gradient g of f at x.  When K=1, fx 
 *                   need only compute f.  When k=2, both f and g are 
 *                   required. 
 * \param data       Structure pointer passed to fx
 * \param xi         The user-supplied initial guess is replaced by 
 *                   the location of the best minimum found by the al- 
 *                   gorithm. 
 * \param xerr       The initial estimate supplied by the user is re- 
 *                   placed by an estimate of the square roots of the 
 *                   diagonal elements of the Hessian matrix evaluated 
 *                   at the best minimum found.  In least-squares ap- 
 *                   plications, assuming that F is the sum of squared 
 *                   residuals, estimates of the standard errors of 
 *                   the unknowns can be obtained by multiplying ERR 
 *                   by the r.m.s. residual. 
 * \param n          The number of unknowns. 
 * \param eps        A small positive number used in tests to set a 
 *                   lower bound on the squared Euclidean norm of 
 *                   vectors considered significantly different from 
 *                   0.  EPS is used in the convergence test.  Usually 
 *                   setting EPS in the range 10**(-12) to 10**(-8) is 
 *                   reasonable.  Very close to a minimum, the algo- 
 *                   rithm generally exhibits a quadratic rate of con- 
 *                   vergence, so setting EPS a few orders of magni- 
 *                   tude too small usually is not too costly. 
 * \param itmax      The maximum number of iterations.  On average, a 
 *                   few evaluations of F and slightly more than one 
 *                   evaluation of G are required at each iteration. 
 * \param fopt       [out] The value of F evaluated at the location of the 
 *                   best minimum that was found. 
 * \param gnopt      [out]  The Euclidean norm of the gradient of the objec- 
 *                   tive function, evaluated at the location of the 
 *                   best minimum that was found. 
 * \param ier        [out] An error flag.  When IER=0, convergence was 
 *                   achieved in ITMAX or fewer iterations; otherwise not. 
 * \param npr        A print flag.  When NPR=0, there is no printout; 
 *                   for NPR=1, the value of F and the Euclidean norm 
 *                   of G, both evaluated at the location of the best 
 *                   minimum found so far, are printed at each itera- 
 *                   tion; for NPR=2, the latter information, together 
 *                   with the location of the best minimum, is printed 
 *                   at each iteration. 
 * \param err        Obit error/message stack
 */
static void
dvdmin (ObitImageFitDataFuncFP fx, gpointer data,  odouble *xi, odouble *xerr, 
	gint n, odouble eps, olong itmax, odouble *fopt, odouble *gnopt, 
	gint *ier, olong npr, ObitErr *err) 
{
  olong   nf, ng, it, i, j, l, iflag, lcount;
  odouble lambda, msq, mu, nsq, nu, **xj=NULL, *x0=NULL, *x=NULL, *k0=NULL, *k=NULL, 
    *s=NULL, *gg=NULL, *m=NULL, *p=NULL, *q=NULL, *wun=NULL, *ax=NULL, tinyc, f, gn, 
    f0p, xx,fp, b0, utu, uts, b, gamma, f0, delta, a, c, alf, t1, t2, t3, t4, t5, t6,
    qtk0, tmp;

  /* existing error */
  if (err->error) return;
  *ier = 0;

  /* create arrays */
  xj = g_malloc0(n*sizeof(odouble*));
  for (i=0; i<n; i++)  xj[i] = g_malloc0(n*sizeof(odouble));
  x0  = g_malloc0(n*sizeof(odouble));
  x   = g_malloc0(n*sizeof(odouble));
  k0  = g_malloc0(n*sizeof(odouble));
  k   = g_malloc0(n*sizeof(odouble));
  s   = g_malloc0(n*sizeof(odouble));
  gg  = g_malloc0(n*sizeof(odouble));
  m   = g_malloc0(n*sizeof(odouble));
  p   = g_malloc0(n*sizeof(odouble));
  q   = g_malloc0(n*sizeof(odouble));
  wun = g_malloc0(n*sizeof(odouble));
  ax  = g_malloc0(n*sizeof(odouble));

  /* If diagnostics requested,give header */
  if (npr>0) 
    Obit_log_error(err, OBIT_InfoErr, 
		   "   dvdmin     debugging information    dvdmin");
  
  /* Initialization: */
  tinyc = 1.0e-3*sqrt (dmachx(1));
  nf = 1;
  ng = 1;
  gn = 0.0;
  it = -1;
  for (i=0; i<n; i++) { /* loop 20 */
    x[i]  = xi[i];
    x0[i] = xi[i];
    for (j=0; j<n; j++) { /* loop 10 */
      xj[j][i] = 0.0;
    } /* end loop  L10:  */
    xj[i][i] = xerr[i];
  } /* end loop  L20:  */

  /* evaluate */
  iflag = 2;
  fx (data, x, &f, gg, iflag);
  f0 = f;
  for (i=0; i<n; i++) { /* loop 40 */
    for (j=0; j<n; j++) { /* loop 30 */
      ax[j] = xj[i][j];
    } /* end loop  L30:  */
    wun[i] = cblas_ddot (n, ax, 1, gg, 1);
    k0[i]  = wun[i];
  } /* end loop  L40:  */

  /* Step 1: Begining of iteration */
 L100:  
  it++; /* iteration count */
  gn = cblas_dnrm2 (n, gg, 1);

  /* print iteration values */
  if (npr >= 1) {
    Obit_log_error(err, OBIT_InfoErr, 
		   "iteration # %4d   f=%16.8g, gradient=%16.8g",
		   it, f0, gn);
    Obit_log_error(err, OBIT_InfoErr, "           parameters,     gradients:");
    for (i=0; i<n; i++) 
      Obit_log_error(err, OBIT_InfoErr, "%4d  %16.8g  %16.8g", 
		     i, x0[i], gg[i]);
    ObitErrLog(err); 
  } /* end if print */

  /* Hit iteration limit? */
  if (it > itmax) {*ier = 1; goto L900;}

  for (i=0; i<n; i++) { /* loop 120 */
    s[i] = -k0[i];
  } /* end loop  L120: */

  f0p = cblas_ddot (n, k0, 1, s, 1);
  lambda = 2.0e0;
  if (4.0e0*f0 < -f0p) {
    xx = -4.0e0*f0 / f0p;
    for (i=0; i<n; i++) { /* loop 130 */
      s[i] *= xx;
    } /* end loop  L130: */
    f0p = -4.0 * f0;
  }

  /* Step 2: */
  lcount = 0;  /* Counter for 200 loop */
 L200:
  if (lcount > itmax) {*ier = 1; goto L900;}
  lcount++;
  /* Update solution */  
   /* fprintf (stdout, " delta ");DEBUG */
  for (i=0; i<n; i++) { /* loop 220 */
    for (j=0; j<n; j++) { /* loop 210 */
      ax[j] = xj[j][i];
    } /* end loop  L210: */
    x[i] = x0[i] + cblas_ddot (n, ax, 1, s, 1);
  /* fprintf (stdout, " %g ",x[i]-x0[i]); DEBUG */
  } /* end loop  L220: */
  /* fprintf (stdout, " \n"); DEBUG */

  /* Converged? */
  if (-f0p < eps) {*ier = 0;  goto L900;}
      
  /* evaluate */
  iflag = 1;
  fx (data, x, &f, gg, iflag);
  nf++;
  if (f >= f0) { 
    /* didn't improve - try smaller increment */
    for (i=0; i<n; i++) { /* loop 240 */
      s[i] *= 0.5;
    } /* end loop  L240: */
    f0p *= 0.5;
    lambda = 0.5;
    goto L200;
  }
  
  /* Step 3: */
  /*L300: */
  /* evaluate */
  iflag = 2;
  fx (data, x, &f, gg, iflag);
  nf++;
  ng++;
  for (i=0; i<n; i++) { /* loop 320 */
    for (j=0; j<n; j++) { /* loop 310 */
      ax[j] = xj[i][j];
    } /* end loop  L310: */
    k[i]  = cblas_ddot (n, ax, 1, gg, 1);
    m[i]  = s[i] + k0[i] - k[i];
    k0[i] = k[i];
    x0[i] = x[i];
  } /* end loop  L320: */
      
  fp = cblas_ddot (n, k, 1, s, 1);
  b0 = fp - f0p;
  f0 = f;
  f0p = fp;
  if (b0 < eps) {
    for (i=0; i<n; i++) { /* loop 330 */
      s[i] *= lambda;
    } /* end loop  L330: */;
    f0p *= lambda;
    goto L200;
  }
      
      /* Step 4: */
  /* L400:  */
  tmp = cblas_dnrm2 (n, m, 1);
  msq = tmp*tmp;
  if (msq < eps) goto L100;
  nu = cblas_ddot (n, m, 1, s, 1);
  mu = nu - msq;
  xx = cblas_ddot (n, m, 1, wun, 1) / msq;
  for (i=0; i<n; i++) { /* loop 410 */
    wun[i] += -xx * m[i];
  } /* end loop  L410: */

  tmp = cblas_dnrm2 (n, wun, 1);
  utu = tmp*tmp;
  xx = cblas_ddot (n, m, 1, wun, 1);
  if ((xx < tinyc)  ||  ((1.0e3*xx)*(1.0e3*xx) < msq*utu)) goto L450;
  for (i=0; i<n; i++) { /* loop 420 */
    wun[i] = 0.0;
  } /* end loop  L420: */
  nsq = 0.0;
  goto L500;
      
  /* Step 4A: */
 L450:     
  uts = cblas_ddot (n, wun, 1, s, 1);
  xx = uts / utu;
  for (i=0; i<n; i++) { /* loop 460 */
    wun[i] *= xx;
  } /* end loop  L460: */
  nsq = uts * xx;
  
  /* Step 5: */
 L500:
  xx = nu / msq;
  b = nsq + mu * xx;
  if (b >= eps) goto L600;
  for (i=0; i<n; i++) { /* loop 510 */
    wun[i] = s[i] - xx * m[i];
  } /* end loop  L510: */
  nsq = b0 - mu * xx;
  b = b0;
      
  /* Step 6: */
 L600:
  if (mu*nu >= msq*nsq) {
    gamma = 0.0;
    delta = sqrt (nu/mu);
    goto L700;
  }
  
      /* Step 6A: */
  /* L650:  */
  a = b - mu;
  c = b + nu;
  gamma = sqrt ((1.0-mu*nu/(msq*nsq))/(a*b));
  delta = sqrt (c/a);
  if (c < a) gamma = -gamma;
  
  /* Step 7: */
 L700:  
  xx = nsq * gamma;
  alf = nu + mu * delta + msq * xx;
  t1 = delta - xx;
  t2 = gamma * nu;
  t3 = (1.0+xx) / alf;
  t4 = -gamma * mu / alf;
  xx = mu*nu/alf;
  t5 = nsq * (1.0e0 + gamma*xx);
  t6 = -(1.0e0+delta) * xx;
  for (i=0; i<n; i++) { /* loop 710 */
    p[i]   = t1*m[i] + t2*wun[i];
    q[i]   = t3*m[i] + t4*wun[i];
    wun[i] = t5*m[i] + t6*wun[i];
  } /* end loop  L710: */

  qtk0 = cblas_ddot (n, q, 1, k0, 1);
  for (i=0; i<n; i++) { /* loop 730 */
    k0[i] += qtk0*p[i];
    for (l=0; l<n; l++) { /* loop 720 */
      ax[l] = xj[l][i];
    } /* end loop  L720: */
    
    xx = cblas_ddot (n, ax, 1, q, 1);
    for (j=0; j<n; j++) { /* loop 730A */
      xj[j][i] +=  xx * p[j];
    } /* end loop  L730A: */
  } /* end loop  L730: */
  
  if (nsq > 0.0) goto L100;
  for (i=0; i<n; i++) { /* loop 740 */
    wun[i] = k0[i];
  } /* end loop  L740: */
  goto L100;

  /* Done Exit: */
 L900:
  for (i=0; i<n; i++) { /* loop  */
	xi[i] = x0[i];
	for (j=0; j<n; j++) { /* loop 910 */
	  ax[j] = xj[j][i];
	} /* end loop  L910: */
	xerr[i] = cblas_dnrm2(n, ax, 1);
  } /* end loop  L920: */
  
  /* return values */
  *fopt  = f0;
  *gnopt = gn;

  /* Display of results wanted? */
  if (npr > 0) {
    if (*ier == 0) Obit_log_error(err, OBIT_InfoErr, "***  convergence achieved.");
    if (*ier == 1) Obit_log_error(err, OBIT_InfoErr, "***  maximum number of iterations reached.");
    Obit_log_error(err, OBIT_InfoErr, "%d function evaluations and %d gradient evaluations.", nf, ng);
    Obit_log_error(err, OBIT_InfoErr, "Solution  objective function %lg, grad. norm %lg:",f0, gn);
    Obit_log_error(err, OBIT_InfoErr, "Solution    parameters,      errors:");
    for (i=0; i<n; i++) 
      Obit_log_error(err, OBIT_InfoErr, "%4d  %16.8g %16.8g", i, xi[i], sqrt(xerr[i]*f0));
    ObitErrLog(err); 
  }

  /* cleanup:*/
  /* delete  arrays */
  if (xj) {
    for (i=0; i<n; i++)  if (xj[i]) g_free(xj[i]);
    g_free(xj);
  }
  if (x0)  g_free(x0);
  if (x)   g_free(x);
  if (k0)  g_free(k0);
  if (k)   g_free(k);
  if (s)   g_free(s);
  if (gg)  g_free(gg);
  if (m)   g_free(m);
  if (p)   g_free(p);
  if (q)   g_free(q);
  if (wun) g_free(wun);
  if (ax)  g_free(ax);
  } /* end of routine dvdmin */ 

/**
 * Machine dependent constants 
 * Adopted from Linpack DMACH
 * \param job:
 *   1 => return a very small number epsilon, just greater than zero
 *   2 => return a very small number, somewhat larger than epsilon
 *   3 => return a very large number
 * \return the requested value
 */ 
static odouble dmachx (gint job)
{
  odouble eps, teps, tiny, huge, s;
 
  /* Epsilon */
  teps = 1.0e0; eps = 0.0;
  s = 1.0 + teps;
  while (s>1.0) {
    eps = teps;
    teps *= 0.5;
    s = 1.0 + teps;
  }

  /* Tiny  */
  s = 1.0;
  while (s!=0.0) {
    tiny = s;;
    s /= 16.0;
  }
  tiny = 100.0*(tiny / eps);

  /* Huge */
  huge = 1.0 / tiny;
  
  if (job == 1) return eps;
  if (job == 2) return tiny;
  return huge;
} /* end of routine dmachx */ 
