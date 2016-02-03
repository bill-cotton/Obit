/* $Id$        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2008,2016                                          */
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

#include "ObitBeamShape.h"
#include "ObitPBUtil.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitBeamShape.c
 * ObitBeamShape class function definitions.
 * This class is derived from the Obit base class.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitBeamShape";

/** Function to obtain parent ClassInfo */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo structure ObitBeamShapeClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitBeamShapeClassInfo myClassInfo = {FALSE};

/*--------------- File Global Variables  ----------------*/


/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitBeamShapeInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitBeamShapeClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitBeamShapeClassInfoDefFn (gpointer inClass);

/** Private: Check/set tabulated beam. */
static void FindTabBeam (ObitBeamShape *in);

/** Private: Interpolate tabulated beam. */
static ofloat GetTabBeam (ObitBeamShape *in, odouble Angle);
/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitBeamShape* newObitBeamShape (gchar* name)
{
  ObitBeamShape* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitBeamShapeClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitBeamShape));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitBeamShapeInit((gpointer)out);

 return out;
} /* end newObitBeamShape */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitBeamShapeGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitBeamShapeClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitBeamShapeGetClass */

/**
 * Make a deep copy of an ObitBeamShape.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitBeamShape* ObitBeamShapeCopy  (ObitBeamShape *in, ObitBeamShape *out, ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  gboolean oldExist;
  gchar *outName;

  /* error checks */
  if (err->error) return out;
  g_assert (ObitIsA(in, &myClassInfo));
  if (out) g_assert (ObitIsA(out, &myClassInfo));

  /* Create if it doesn't exist */
  oldExist = out!=NULL;
  if (!oldExist) {
    /* derive object name */
    outName = g_strconcat ("Copy: ",in->name,NULL);
    out = newObitBeamShape(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /*  copy this class */
  out->myDesc  = ObitImageRef(in->myDesc);
  out->pbmin   = in->pbmin;
  out->antSize = in->antSize;
  out->doGain  = in->doGain;
  out->refFreq = in->refFreq;
  out->doJinc  = in->doJinc;
  ObitImageDescGetPoint(out->myDesc, &out->raPnt, &out->decPnt) ;
  return out;
} /* end ObitBeamShapeCopy */

/**
 * Make a copy of a object but do not copy the actual data
 * This is useful to create an BeamShape similar to the input one.
 * \param in  The object to copy
 * \param out An existing object pointer for output, must be defined.
 * \param err Obit error stack object.
 */
void ObitBeamShapeClone  (ObitBeamShape *in, ObitBeamShape *out, ObitErr *err)
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
  out->myDesc  = ObitImageRef(in->myDesc);
  out->pbmin   = in->pbmin;
  out->antSize = in->antSize;
  out->doGain  = in->doGain;
  out->refFreq = in->refFreq;
  out->doJinc  = in->doJinc;
  ObitImageDescGetPoint(out->myDesc, &out->raPnt, &out->decPnt) ;
} /* end ObitBeamShapeClone */

/**
 * Creates an ObitBeamShape 
 * \param name    An optional name for the object.
 * \param image   Image for which beam shape is desired
 *                Control on info member:
 * \li doTab      If TRUE use tabulated beam if available [def FALSE]
 * \param pbmin   Minimum gain, lower values will be clipped at this value
 * \param antSize Size of Antenna in (m)
 * \param doGain  If true gain wanted, else gain set to 1.0
 * \return the new object.
 */
ObitBeamShape* ObitBeamShapeCreate (gchar* name, ObitImage *image, 
				    ofloat pbmin, ofloat antSize, 
				    gboolean doGain)
{
  ObitBeamShape* out;
  gboolean doTab=FALSE;
  gint32       dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;

  /* Create basic structure */
  out = newObitBeamShape (name);
  out->myDesc  = ObitImageRef(image->myDesc);
  out->pbmin   = pbmin;
  out->antSize = antSize;
  out->doGain  = doGain;
  out->refFreq = out->myDesc->crval[out->myDesc->jlocf] + 
    out->myDesc->cdelt[out->myDesc->jlocf] * 
    (out->myDesc->crpix[out->myDesc->jlocf] - out->myDesc->plane);
  out->doJinc  = out->refFreq >= 1.0e9;
  ObitImageDescGetPoint(out->myDesc, &out->raPnt, &out->decPnt) ;
  out->raPnt  *= DG2RAD;  /* to radians */
  out->decPnt *= DG2RAD;  /* to radians */
  /* tabulated beam? */
  ObitInfoListGetTest(image->info, "doTab", &type, dim, &doTab);
  if (doTab) FindTabBeam(out);  /* Use if available */
  return out;
} /* end ObitBeamShapeCreate */

/**
 * Calculate gain in a given direction.
 * Simple function of distance from pointing center.
 * \param in     the BeamShape object
 * \param ra     RA (deg) of direction for gain
 * \param dec    Declination (deg) of direction for gain
 * \param parAng Parallactic angle (rad) NYI
 * \return Gain
 */
ofloat ObitBeamShapeGain (ObitBeamShape *in, odouble ra, odouble dec, 
			  ofloat parAng)
{
  ofloat gain=1.0;
  odouble xx, yy, zz, Angle;

  /* Compute gain? */
  if (!in->doGain) return gain;

  /* get offset */
  xx = DG2RAD * (ra);
  yy = DG2RAD * (dec);
  zz = sin (yy) * sin (in->decPnt) + cos (yy) * cos (in->decPnt) * cos (xx-in->raPnt);
  zz = MIN (zz, 1.000);
  Angle = acos (zz) * RAD2DG;

  /* Compute */
  if (in->doTab)  gain = GetTabBeam (in, Angle);
  else if (in->doJinc) gain = ObitPBUtilJinc(Angle, in->refFreq, in->antSize, in->pbmin);
  else                 gain = ObitPBUtilPoly(Angle, in->refFreq, in->pbmin);
  return gain;
} /* end ObitBeamShapeGain */

/**
 * Calculate gain in a given offset from a symmetric beam shape.
 * Simple function of distance from pointing center.
 * \param in     the BeamShape object
 * \param Angle  Angular distance (deg) from pointing center
 * \return power gain
 */
ofloat ObitBeamShapeGainSym (ObitBeamShape *in, odouble Angle)
{
  ofloat gain=1.0;

  /* Compute gain? */
  if (!in->doGain) return gain;

  /* Compute */
  if (in->doTab)       gain = GetTabBeam (in, Angle);
  else if (in->doJinc) gain = ObitPBUtilJinc(Angle, in->refFreq, in->antSize, in->pbmin);
  else                 gain = ObitPBUtilPoly(Angle, in->refFreq, in->pbmin);
  return gain;
} /* end ObitBeamShapeGainSym */

/**
 * Calculate angular distance from pointing center
 * \param in     the BeamShape object
 * \param ra     RA (deg) of direction for angle
 * \param dec    Declination (deg) of direction for angle
 * \param parAng Parallactic angle (rad) NYI
 * \return Angle from pointing center (deg)
 */
odouble ObitBeamShapeAngle (ObitBeamShape *in, odouble ra, odouble dec, 
			    ofloat parAng)
{
  odouble xx, yy, zz, Angle=0.0;

  /* get offset */
  xx = DG2RAD * (ra);
  yy = DG2RAD * (dec);
  zz = sin (yy) * sin (in->decPnt) + cos (yy) * cos (in->decPnt) * cos (xx-in->raPnt);
  zz = MIN (zz, 1.000);
  Angle = acos (zz) * RAD2DG;
  return Angle;
} /* end ObitBeamShapeAngle */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitBeamShapeClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitBeamShapeClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitBeamShapeClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitBeamShapeClassInfoDefFn (gpointer inClass)
{
  ObitBeamShapeClassInfo *theClass = (ObitBeamShapeClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitBeamShapeClassInit;
  theClass->newObit       = (newObitFP)newObitBeamShape;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitBeamShapeClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitBeamShapeGetClass;
  theClass->ObitCopy      = (ObitCopyFP)ObitBeamShapeCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitBeamShapeClear;
  theClass->ObitInit      = (ObitInitFP)ObitBeamShapeInit;
  theClass->ObitBeamShapeCreate = (ObitBeamShapeCreateFP)ObitBeamShapeCreate;
  theClass->ObitBeamShapeGain   = (ObitBeamShapeGainFP)ObitBeamShapeGain;
  theClass->ObitBeamShapeGainSym= (ObitBeamShapeGainSymFP)ObitBeamShapeGainSym;
  theClass->ObitBeamShapeAngle  = (ObitBeamShapeAngleFP)ObitBeamShapeAngle;

} /* end ObitBeamShapeClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitBeamShapeInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitBeamShape *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->myDesc     = NULL;
  in->myFI       = NULL;
  in->pbmin      = 0.0;
  in->antSize    = 0.0;
  in->icellSize  = 0.0;
  in->doGain     = FALSE;
  in->raPnt      = 0.0;
  in->decPnt     = 0.0;
  in->refFreq    = 0.0;
  in->itabRefFreq= 0.0;
  in->doJinc     = FALSE;
  in->doTab      = FALSE;
} /* end ObitBeamShapeInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitBeamShape* cast to an Obit*.
 */
void ObitBeamShapeClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitBeamShape *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  in->myDesc  = ObitImageDescUnref(in->myDesc);
  in->myFI    = ObitFInterpolateUnref(in->myFI);

  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitBeamShapeClear */

/**
 * Check if a tabulated beam is available and if so enable it
 * Currently implemented: VLA P, S Band
 * \param in     the BeamShape object
 */
static void FindTabBeam (ObitBeamShape *in)
{
  olong i, j, itab, ndim=1, naxis[] = {1};
  ObitFArray *tFA = NULL;
  olong   ntab   = 2;            /* How many? */
  olong   hwidth = 2;            /* Interpolation half width */
  /*                    P      S */
  odouble minFreq[] = {200.0e6, 1.8e9, };  /* Lower bound of bands */
  odouble maxFreq[] = {500.0e6, 4.2e9, };  /* Upper  bound of bands */
  odouble refFreq[] = {340.0e6, 3.0e9, };  /* Tabulated reference freq */

  /* Have a tabulated beam? */
  itab = -1;  /* Which beam */
  for (i=0; i<ntab; i++) {
    if ((in->refFreq>=minFreq[i]) && (in->refFreq<=maxFreq[i])) {
      itab = i;
      break;
    }
  }
  /* Find one? */
  if (itab<0) return; /* Nope */
  /* Fill in details */
  switch (itab) {
  case 0:     /* P band */
    in->itabRefFreq = 1.0/refFreq[itab]; /* 1/tabulated ref freq */
    in->icellSize   = 3600.0/250.0;     /* 1/Tabulated cell spacing */
    olong  Pncell   = 300;
    /* 10 Nov2015, P Band 300 entries  cell 250.0/3600, refFreq  = 340.0e6 */
    ofloat Pbeam[]  = {    /* Fitted to beam */
      1.000000, 0.998248, 0.993048, 0.984058, 0.970796, 0.955100, 0.939748, 0.912643, 0.886224,
      0.867500, 0.830706, 0.793851, 0.771507, 0.747092, 0.722266, 0.695144, 0.666975, 0.638199,
      0.608038, 0.577563, 0.546267, 0.514786, 0.483709, 0.451885, 0.421084, 0.391206, 0.361235,
      0.332097, 0.303311, 0.275288, 0.248056, 0.220403, 0.193054, 0.166338, 0.138971, 0.112269,
      0.086202, 0.061056, 0.035924, 0.012253,-0.011110,-0.033139,-0.056696,-0.074705,-0.091836,
     -0.108188,-0.119736,-0.130505,-0.138099,-0.144492,-0.148338,-0.152569,-0.154990,-0.157243,
     -0.158610,-0.159265,-0.159365,-0.158356,-0.156740,-0.154839,-0.153190,-0.151059,-0.148089,
     -0.144850,-0.141322,-0.137659,-0.133187,-0.128910,-0.124340,-0.119538,-0.114051,-0.108139,
     -0.103371,-0.098704,-0.093910,-0.089438,-0.084935,-0.080243,-0.075727,-0.071550,-0.067402,
     -0.062761,-0.058984,-0.055482,-0.052036,-0.048765,-0.045855,-0.043183,-0.041011,-0.039075,
     -0.037319,-0.035561,-0.033798,-0.032735,-0.031881,-0.031131,-0.030240,-0.029644,-0.029222,
     -0.028921,-0.028558,-0.028655,-0.028869,-0.029189,-0.029516,-0.030054,-0.030669,-0.031355,
     -0.032283,-0.033223,-0.034139,-0.034954,-0.035681,-0.036417,-0.037149,-0.037820,-0.038435,
     -0.038967,-0.039466,-0.039883,-0.040226,-0.040512,-0.040651,-0.040659,-0.040359,-0.040056,
     -0.039284,-0.038972,-0.038596,-0.038105,-0.037470,-0.036376,-0.035456,-0.033998,-0.032957,
     -0.031710,-0.030480,-0.029148,-0.027525,-0.025789,-0.024095,-0.022282,-0.020765,-0.019174,
     -0.017573,-0.015583,-0.013725,-0.012100,-0.010378,-0.008902,-0.007481,-0.006150,-0.004723,
     -0.003377,-0.002155,-0.001147,-0.000536, 0.000203, 0.000865, 0.001340, 0.001923, 0.002231,
      0.002408, 0.002378, 0.002177, 0.001899, 0.001426, 0.001066, 0.000533, 0.000073,-0.000631,
     -0.001448,-0.002298,-0.003168,-0.003954,-0.004920,-0.005864,-0.006699,-0.007743,-0.008659,
     -0.009452,-0.010132,-0.010690,-0.011086,-0.011309,-0.011851,-0.012474,-0.012585,-0.012583,
     -0.012258,-0.011698,-0.010891,-0.010301,-0.009587,-0.008786,-0.007685,-0.006226,-0.004800,
     -0.003321,-0.002036,-0.000653, 0.000758, 0.002152, 0.003543, 0.004861, 0.006130, 0.007524,
      0.008955, 0.010038, 0.011105, 0.013919, 0.014813, 0.015765, 0.016394, 0.016900, 0.017360,
      0.017642, 0.018007, 0.018134, 0.018125, 0.018098, 0.017966, 0.017623, 0.017367, 0.017176,
      0.016799, 0.016424, 0.015927, 0.015253, 0.014508, 0.013456, 0.012689, 0.011830, 0.011033,
      0.010307, 0.009657, 0.009023, 0.008565, 0.007511, 0.006803, 0.006060, 0.005279, 0.004655,
      0.004088, 0.003540, 0.003182, 0.003064, 0.002707, 0.002566, 0.002285, 0.002189, 0.002264,
      0.002423, 0.002526, 0.002785, 0.003154, 0.003435, 0.003552, 0.004070, 0.005006, 0.005364,
      0.005837, 0.006182, 0.006740, 0.007079, 0.007322, 0.007627, 0.007587, 0.007766, 0.007941,
      0.008522, 0.008887, 0.008494, 0.008300, 0.008574, 0.008664, 0.008907, 0.009141, 0.009420,
      0.009731, 0.009579, 0.009681, 0.009798, 0.008858, 0.007708, 0.007665, 0.007506, 0.007470,
      0.007620, 0.007608, 0.007569, 0.006601, 0.006494, 0.006344, 0.005570, 0.006879, 0.006602,
      0.007726, 0.007875};
    naxis[0] = Pncell;   /* Create/fill FArray with beam shape for interpolator */
    tFA = ObitFArrayCreate("TempFA", ndim, naxis);
    for (j=0; j<Pncell; j++) tFA->array[j] = Pbeam[j]*Pbeam[j];  /* Voltage to power */
    in->myFI = newObitFInterpolateCreate ("Interp", tFA, in->myDesc, hwidth);
    tFA = ObitFArrayUnref(tFA);    /* Unreference temp FA */
    break;
  case 1:     /* S band */
    in->itabRefFreq = 1.0/refFreq[itab]; /* 1/tabulated ref freq */
    in->icellSize   = 3600.0/30.0;       /* 1/Tabulated cell spacing */
    olong  Sncell  = 250;
    /* 15 Dec 15 S Band 250 entries, cell 30.0/3600, refFreq  = 3.0e9 */
    ofloat Sbeam[] = {    /* Fitted to voltage beam */
       1.000000, 0.997801, 0.991191, 0.979534, 0.963917, 0.945477, 0.925469, 0.902924,
       0.877989, 0.850618, 0.820960, 0.788785, 0.754457, 0.718072, 0.679851, 0.639915,
       0.598547, 0.556108, 0.512954, 0.469143, 0.425074, 0.381071, 0.337314, 0.293631,
       0.251205, 0.210151, 0.170330, 0.131907, 0.095075, 0.060002, 0.026808,-0.004369,
      -0.033414,-0.060111,-0.084383,-0.106189,-0.125516,-0.142382,-0.156810,-0.168836,
      -0.178508,-0.185931,-0.191124,-0.194186,-0.195177,-0.194241,-0.191455,-0.186927,
      -0.180851,-0.173404,-0.164744,-0.155074,-0.144559,-0.133349,-0.121636,-0.109558,
      -0.097257,-0.084910,-0.072642,-0.060605,-0.048917,-0.037715,-0.027097,-0.017189,
      -0.008081, 0.000135, 0.007436, 0.013753, 0.019071, 0.023365, 0.026618, 0.028834,
       0.030028, 0.030255, 0.029542, 0.027926, 0.025465, 0.022204, 0.018204, 0.013505,
       0.008183, 0.002350,-0.003892,-0.010405,-0.017047,-0.023758,-0.030407,-0.036969,
      -0.043382,-0.049580,-0.055471,-0.060975,-0.066024,-0.070577,-0.074561,-0.078006,
      -0.080809,-0.082940,-0.084391,-0.085216,-0.085304,-0.084722,-0.083528,-0.081607,
      -0.079001,-0.075808,-0.072112,-0.067903,-0.063302,-0.058277,-0.052895,-0.047243,
      -0.041358,-0.035287,-0.029134,-0.022834,-0.016636,-0.010580,-0.004612, 0.001007,
       0.006283, 0.011364, 0.015945, 0.020120, 0.024015, 0.027299, 0.030071, 0.032507,
       0.034261, 0.035475, 0.036158, 0.036473, 0.036126, 0.035281, 0.034121, 0.032342,
       0.030147, 0.027706, 0.024814, 0.021625, 0.018355, 0.014773, 0.011045, 0.007359,
       0.003429,-0.000425,-0.004357,-0.008201,-0.011823,-0.015388,-0.018852,-0.021985,
      -0.024917,-0.027582,-0.029917,-0.031863,-0.033594,-0.034882,-0.035838,-0.036543,
      -0.036923,-0.036918,-0.036655,-0.036032,-0.035119,-0.034013,-0.032658,-0.031008,
      -0.029239,-0.027303,-0.025200,-0.022939,-0.020680,-0.018365,-0.015953,-0.013662,
      -0.011393,-0.009033,-0.006836,-0.004746,-0.002665,-0.000855, 0.000758, 0.002372,
       0.003605, 0.004623, 0.005435, 0.006260, 0.006732, 0.007053, 0.007470, 0.007507,
       0.007187, 0.006721, 0.006158, 0.005485, 0.004737, 0.003936, 0.003081, 0.002209,
       0.001331, 0.000459,-0.000365,-0.001158,-0.001882,-0.002531,-0.003091,-0.003569,
      -0.003950,-0.004220,-0.004365,-0.004399,-0.004327,-0.004131,-0.003800,-0.003388,
      -0.002845,-0.002230,-0.001516,-0.000734, 0.000109, 0.001006, 0.001937, 0.002879,
       0.003840, 0.004793, 0.005728, 0.006657, 0.007614, 0.008493, 0.009317, 0.010081,
       0.010751, 0.011335, 0.011832, 0.012213, 0.012516, 0.012726, 0.012827, 0.012838,
       0.012723, 0.012484, 0.012118, 0.011647, 0.011070, 0.010442, 0.009695, 0.008852,
       0.007935, 0.006959};
    naxis[0] = Sncell;   /* Create/fill FArray with beam shape for interpolator */
    tFA = ObitFArrayCreate("TempFA", ndim, naxis);
    for (j=0; j<Sncell; j++) tFA->array[j] = Sbeam[j]*Sbeam[j];  /*  Voltage to power */
    in->myFI = newObitFInterpolateCreate ("Interp", tFA, in->myDesc, hwidth);
    tFA = ObitFArrayUnref(tFA);    /* Unreference temp FA */
   break;
  default:    /* Doh */
    return;
  }; /* end switch */
  in->doTab = TRUE;    /* Must be OK if it gets here */
} /* end FindTabBeam */

/**
 * Check if a tabulated beam is available and if so enable it
 * Currently implemented: VLA S Band
 * Returns pbmin if out of tabulated beam shape
 * \param in  the BeamShape object
 * \param in  Angle from pointing in degrees
 * \return  beam power gain
 */
static ofloat GetTabBeam (ObitBeamShape *in, odouble Angle)
{
  ofloat pixel;
  pixel = 1.0+Angle*in->icellSize*in->itabRefFreq*in->refFreq;
  if (pixel<=1.1) return 1.0;                /* Trap center */
  if (pixel>in->myFI->nx) return in->pbmin;  /* Beyond tabulation */
  return ObitFInterpolate1D(in->myFI, pixel);
} /* end GetTabBeam */
