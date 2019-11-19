/* $Id$        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2008-2019                                          */
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
#include "ObitImageMF.h"

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

/** Private: Check/set MeerKat tabulated beam. */
static void MeerKATTabBeam (ObitBeamShape *in);

/** Private: Check/set VLITE tabulated beam. */
static void FindVLITEBeam (ObitBeamShape *in);
 
/** Private: MeerKAT beam. */
static ofloat GetMKBeam (ObitBeamShape *in, odouble Angle);

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
  out->doJinc  = in->doJinc;
  out->doTab   = in->doTab;
  out->doVLITE = in->doVLITE;
  out->doMeerKAT = in->doMeerKAT;
  out->beamAng = in->beamAng;
  out->refFreq = in->refFreq;
  out->itabRefFreq = in->itabRefFreq;
  out->icellSize   = in->icellSize;
  if (in->myFI) out->myFI = ObitFInterpolateCopy(in->myFI,  out->myFI, err);
  ObitImageDescGetPoint(out->myDesc, &out->raPnt, &out->decPnt) ;
  out->raPnt      *= DG2RAD;
  out->decPnt     *= DG2RAD;
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
  out->doTab   = in->doTab;
  out->doVLITE = in->doVLITE;
  out->doMeerKAT = in->doMeerKAT;
  out->beamAng = in->beamAng;
  out->refFreq = in->refFreq;
  out->itabRefFreq = in->itabRefFreq;
  out->icellSize   = in->icellSize;
  out->myFI        = ObitFInterpolateCopy(in->myFI,  out->myFI, err);
  ObitImageDescGetPoint(out->myDesc, &out->raPnt, &out->decPnt) ;
  out->raPnt      *= DG2RAD;
  out->decPnt     *= DG2RAD;
} /* end ObitBeamShapeClone */

/**
 * Creates an ObitBeamShape 
 * \param name    An optional name for the object.
 * \param image   Image for which beam shape is desired
 *                Control on info member:
 * \li doTab      If TRUE use tabulated beam if available [def FALSE]
 *                Traps MeerKAT case, uses Tabulated beam
 *                Traps VLITE case, uses Tabulated beam
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
  gboolean doTab=FALSE, isMeerKAT=FALSE, doVLITE=FALSE;
  gint32   dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;

  /* Create basic structure */
  out = newObitBeamShape (name);
  out->myDesc  = ObitImageRef(image->myDesc);
  out->pbmin   = pbmin;
  out->antSize = antSize;
  out->beamAng = -1.0;
  out->doGain  = doGain;
  out->refFreq = ObitImageMFGetPlaneFreq(image);
  out->doJinc  = out->refFreq >= 1.0e9;
  ObitImageDescGetPoint(out->myDesc, &out->raPnt, &out->decPnt) ;
  out->raPnt  *= DG2RAD;  /* to radians */
  out->decPnt *= DG2RAD;  /* to radians */
  /* tabulated beam? */
  ObitInfoListGetTest(image->info, "doTab", &type, dim, &doTab);
  isMeerKAT = !strncmp(image->myDesc->teles, "MeerKAT",7); /* MeerKAT */
  out->doMeerKAT = isMeerKAT;
  ObitInfoListGetTest(image->info, "doVLITE", &type, dim, &doVLITE); /* VLITE */
  doVLITE = doVLITE || !strncmp(image->myDesc->instrument, "VLITE",5);
   /* if (isMeerKAT)   MeerKATTabBeam(out);  Always use for MeerKAT */
  if (doVLITE) FindVLITEBeam(out);   /* Use VLITE beam  */
  else if (doTab)   FindTabBeam(out);     /* Use standard if available */

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
  gboolean doKAT=FALSE;

  /* Compute gain? */
  if (!in->doGain) return gain;
  doKAT = !strncmp(in->myDesc->teles, "KAT-7",5); /* Kat-7 */

  /* Compute */
  if (in->doTab)          gain = GetTabBeam (in, Angle);
  else if (in->doVLITE)   gain = GetTabBeam (in, Angle);
  else if (in->doMeerKAT) gain = GetMKBeam (in, Angle);
  else if (doKAT)         gain = ObitPBUtilKAT7 (Angle, in->refFreq, 0.0);
  else if (in->doJinc)    gain = ObitPBUtilJinc(Angle, in->refFreq, in->antSize, in->pbmin);
  else                    gain = ObitPBUtilPoly(Angle, in->refFreq, in->pbmin);
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
 * Currently implemented: VLA P, S Band, C Band
 * \param in     the BeamShape object
 */
static void FindTabBeam (ObitBeamShape *in)
{
  olong i, j, itab, ndim=1, naxis[] = {1};
  ObitFArray *tFA = NULL;
  olong   ntab   = 4;            /* How many tables? */
  olong   hwidth = 2;            /* Interpolation half width */
  /*                      P      S       C     X*/
  odouble minFreq[] = {200.0e6, 1.8e9, 2.9e9, 7.90e9, };  /* Lower bound of bands */
  odouble maxFreq[] = {500.0e6, 3.8e9, 7.8e9, 11.8e9, };  /* Upper  bound of bands */
  odouble refFreq[] = {340.0e6, 3.0e9, 6.0e9, 10.0e9, };  /* Tabulated reference freq */

  /* Have a tabulated beam? */
  in->doTab = FALSE;
  itab = -1;  /* Which beam - pick hightest that matches */
  for (i=0; i<ntab; i++) {
    if ((in->refFreq>=minFreq[i]) && (in->refFreq<=maxFreq[i])) {
      itab = i;
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
    in->icellSize   = 3600.0/30.000000;
    olong  Sncell  = 271;
    /* 26 Mar 17 S Band 271 entries, cell 30.0/3600, refFreq  = 3.0e9 */
    ofloat Sbeam[] = {    /* Fitted to voltage beam */
     0.999916,   0.998041,   0.992403,   0.983417,   0.970690,   0.954605,   0.934347,   0.911324,
      0.884635,   0.856180,   0.825812,   0.793804,   0.760268,   0.725417,   0.689476,   0.652862,
      0.615359,   0.576879,   0.537562,   0.497731,   0.457156,   0.416678,   0.376476,   0.336644,
      0.297254,   0.257873,   0.220638,   0.184037,   0.148640,   0.114083,   0.080538,   0.047727,
      0.015873,  -0.014357,  -0.042714,  -0.068528,  -0.091635,  -0.111966,  -0.129744,  -0.145239,
     -0.158594,  -0.169904,  -0.179086,  -0.186182,  -0.191374,  -0.194697,  -0.196225,  -0.196027,
     -0.194237,  -0.190898,  -0.186032,  -0.179821,  -0.172408,  -0.163978,  -0.154683,  -0.144632,
     -0.134019,  -0.122921,  -0.111480,  -0.099815,  -0.088034,  -0.076280,  -0.064696,  -0.053410,
     -0.042514,  -0.032083,  -0.022240,  -0.013096,  -0.004734,   0.002770,   0.009357,   0.015028,
      0.019776,   0.023592,   0.026536,   0.028587,   0.029759,   0.030112,   0.029639,   0.028348,
      0.026270,   0.023397,   0.019699,   0.015186,   0.009870,   0.003825,  -0.002823,  -0.009938,
     -0.017318,  -0.024711,  -0.031917,  -0.038710,  -0.045045,  -0.050931,  -0.056429,  -0.061579,
     -0.066376,  -0.070815,  -0.074835,  -0.078350,  -0.081389,  -0.083807,  -0.085635,  -0.086936,
     -0.087541,  -0.087519,  -0.086935,  -0.085700,  -0.083826,  -0.081410,  -0.078385,  -0.074852,
     -0.070854,  -0.066469,  -0.061740,  -0.056791,  -0.051530,  -0.046004,  -0.040204,  -0.034114,
     -0.027919,  -0.021592,  -0.015040,  -0.008719,  -0.002578,   0.003244,   0.008796,   0.013744,
      0.018268,   0.022574,   0.026286,   0.029612,   0.032737,   0.035167,   0.037087,   0.038687,
      0.039559,   0.039934,   0.039820,   0.039373,   0.038304,   0.036761,   0.034935,   0.032771,
      0.030191,   0.027111,   0.023793,   0.020436,   0.016751,   0.012917,   0.009107,   0.004987,
      0.000755,  -0.003402,  -0.007737,  -0.012027,  -0.015965,  -0.019601,  -0.023190,  -0.026405,
     -0.029020,  -0.031540,  -0.033613,  -0.035416,  -0.037009,  -0.038289,  -0.039113,  -0.039709,
     -0.039968,  -0.039834,  -0.039469,  -0.038713,  -0.037749,  -0.036427,  -0.034956,  -0.033265,
     -0.031260,  -0.029195,  -0.027000,  -0.024571,  -0.022189,  -0.019780,  -0.017147,  -0.014622,
     -0.012086,  -0.009569,  -0.006783,  -0.004316,  -0.001946,   0.000559,   0.002634,   0.004491,
      0.006123,   0.007209,   0.008059,   0.008661,   0.009044,   0.009228,   0.009216,   0.009071,
      0.008723,   0.008216,   0.007570,   0.006781,   0.005950,   0.005064,   0.004132,   0.003156,
      0.002169,   0.001178,   0.000197,  -0.000751,  -0.001665,  -0.002525,  -0.003298,  -0.003972,
     -0.004541,  -0.004987,  -0.005301,  -0.005507,  -0.005512,  -0.005331,  -0.005064,  -0.004702,
     -0.004187,  -0.003618,  -0.002920,  -0.002114,  -0.001251,  -0.000301,   0.000719,   0.001817,
      0.002941,   0.004111,   0.005328,   0.006492,   0.007619,   0.008720,   0.009696,   0.010666,
      0.011536,   0.012327,   0.012986,   0.013538,   0.013977,   0.014347,   0.014608,   0.014759,
      0.014807,   0.014793,   0.014674,   0.014425,   0.014052,   0.013554,   0.012892,   0.012219,
      0.011340,   0.010320,   0.009199,   0.007984,   0.006676,   0.005444,   0.003929,   0.002249,
      0.001276,  -0.000303,  -0.001747,  -0.003157,  -0.004322,  -0.005987,   0.000000};
    naxis[0] = Sncell;   /* Create/fill FArray with beam shape for interpolator */
    tFA = ObitFArrayCreate("TempFA", ndim, naxis);
    for (j=0; j<Sncell; j++) tFA->array[j] = Sbeam[j]*Sbeam[j];  /*  Voltage to power */
    in->myFI = newObitFInterpolateCreate ("Interp", tFA, in->myDesc, hwidth);
    tFA = ObitFArrayUnref(tFA);    /* Unreference temp FA */
   break;
  case 2:     /* C band */
    in->itabRefFreq = 1.0/refFreq[itab]; /* 1/tabulated ref freq */
    /* 28 Mar 17 C Band 296 entries, cell 15.0/3600, refFreq  = 6.0e9 */
    in->icellSize   = 3600.0/15.000000;
    olong  Cncell  = 296;
    ofloat Cbeam[] = {    /* Fitted to voltage beam */
      1.000000,   0.998789,   0.994080,   0.986285,   0.975384,   0.961374,   0.944553,   0.924910,
      0.902562,   0.877663,   0.850386,   0.820879,   0.789437,   0.756136,   0.721146,   0.684763,
      0.647141,   0.608446,   0.568897,   0.528680,   0.488029,   0.447100,   0.406164,   0.365446,
      0.325158,   0.285469,   0.246629,   0.208829,   0.172168,   0.136825,   0.102849,   0.070317,
      0.039407,   0.010228,  -0.017136,  -0.042427,  -0.065456,  -0.086199,  -0.104568,  -0.120588,
     -0.134339,  -0.145913,  -0.155426,  -0.163000,  -0.168658,  -0.172445,  -0.174470,  -0.174832,
     -0.173636,  -0.170935,  -0.166869,  -0.161607,  -0.155257,  -0.148005,  -0.139902,  -0.131095,
     -0.121637,  -0.111607,  -0.101219,  -0.090550,  -0.079782,  -0.069028,  -0.058428,  -0.048158,
     -0.038355,  -0.029210,  -0.020800,  -0.013157,  -0.006300,  -0.000255,   0.004986,   0.009438,
      0.013087,   0.015924,   0.017935,   0.019112,   0.019470,   0.019012,   0.017770,   0.015805,
      0.013110,   0.009778,   0.005789,   0.001262,  -0.003725,  -0.009085,  -0.014755,  -0.020604,
     -0.026551,  -0.032540,  -0.038530,  -0.044399,  -0.050128,  -0.055642,  -0.060897,  -0.065805,
     -0.070355,  -0.074453,  -0.078075,  -0.081191,  -0.083775,  -0.085815,  -0.087299,  -0.088217,
     -0.088577,  -0.088379,  -0.087639,  -0.086368,  -0.084577,  -0.082278,  -0.079505,  -0.076293,
     -0.072643,  -0.068613,  -0.064241,  -0.059551,  -0.054626,  -0.049495,  -0.044240,  -0.038911,
     -0.033604,  -0.028336,  -0.023192,  -0.018217,  -0.013461,  -0.008930,  -0.004665,  -0.000688,
      0.002973,   0.006294,   0.009277,   0.011864,   0.014055,   0.015842,   0.017217,   0.018164,
      0.018701,   0.019037,   0.018972,   0.018683,   0.017900,   0.016961,   0.015543,   0.013832,
      0.012005,   0.009745,   0.007464,   0.004897,   0.002173,  -0.000486,  -0.003350,  -0.006102,
     -0.009004,  -0.011703,  -0.014448,  -0.017078,  -0.019304,  -0.021631,  -0.023808,  -0.025588,
     -0.027300,  -0.028618,  -0.029814,  -0.030739,  -0.031292,  -0.031673,  -0.031669,  -0.031487,
     -0.031044,  -0.030205,  -0.029338,  -0.027926,  -0.026508,  -0.024915,  -0.023063,  -0.021161,
     -0.019012,  -0.016874,  -0.014643,  -0.012640,  -0.010466,  -0.008312,  -0.006204,  -0.004161,
     -0.002204,  -0.000364,   0.001382,   0.002985,   0.004264,   0.005533,   0.006616,   0.007360,
      0.008067,   0.008618,   0.008977,   0.009145,   0.009142,   0.008944,   0.008635,   0.008167,
      0.007543,   0.006811,   0.005868,   0.004774,   0.003739,   0.002787,   0.001715,   0.000592,
     -0.000547,  -0.001565,  -0.002836,  -0.004057,  -0.005121,  -0.006138,  -0.007065,  -0.007837,
     -0.008524,  -0.009049,  -0.009474,  -0.009756,  -0.009883,  -0.009882,  -0.009763,  -0.009492,
     -0.009052,  -0.008515,  -0.007905,  -0.007152,  -0.006318,  -0.005464,  -0.004476,  -0.003449,
     -0.002450,  -0.001346,  -0.000213,   0.000831,   0.001945,   0.002987,   0.004041,   0.005036,
      0.005910,   0.006812,   0.007536,   0.008197,   0.008741,   0.009200,   0.009532,   0.009745,
      0.009920,   0.009936,   0.009823,   0.009626,   0.009296,   0.008842,   0.008358,   0.007524,
      0.006567,   0.005737,   0.004643,   0.003468,   0.002529,   0.001292,   0.000357,  -0.000909,
     -0.002171,  -0.003078,  -0.003966,  -0.004995,  -0.005696,  -0.006633,  -0.007247,  -0.007678,
     -0.008085,  -0.008197,  -0.008363,  -0.008267,  -0.008043,  -0.008026,  -0.007602,  -0.007138,
     -0.006989,  -0.006362,  -0.006213,  -0.005359,  -0.004543,  -0.004222,  -0.003133,  -0.002014,
     -0.001619,   0.000667,   0.001032,   0.002247,   0.003623,   0.003929,   0.004971,   0.000000};
    naxis[0] = Cncell;   /* Create/fill FArray with beam shape for interpolator */
    tFA = ObitFArrayCreate("TempFA", ndim, naxis);
    for (j=0; j<Cncell; j++) tFA->array[j] = Cbeam[j]*Cbeam[j];  /*  Voltage to power */
    in->myFI = newObitFInterpolateCreate ("Interp", tFA, in->myDesc, hwidth);
    tFA = ObitFArrayUnref(tFA);    /* Unreference temp FA */
   break;
  case 3:     /* X band */
    in->itabRefFreq = 1.0/refFreq[itab]; /* 1/tabulated ref freq */
    /* 5 April 2017 fittes to beam cubes */
    in->icellSize   = 3600.0/6.000000;
    olong  Xncell  = 396;
    ofloat Xbeam[] = {    /* Fitted to voltage beam */
      1.001476,   1.000759,   0.998609,   0.995033,   0.989962,   0.983529,   0.975584,   0.966563,
      0.956258,   0.944680,   0.931881,   0.917895,   0.902752,   0.886478,   0.869062,   0.850571,
      0.831144,   0.810802,   0.789539,   0.767361,   0.744272,   0.720472,   0.695998,   0.671065,
      0.645587,   0.619629,   0.593261,   0.566494,   0.539517,   0.512347,   0.485009,   0.457614,
      0.430212,   0.402867,   0.375610,   0.348578,   0.321758,   0.295294,   0.269221,   0.243516,
      0.218281,   0.193576,   0.169430,   0.145923,   0.122978,   0.100713,   0.079126,   0.058239,
      0.038053,   0.018605,  -0.000079,  -0.017964,  -0.034976,  -0.051145,  -0.066409,  -0.080668,
     -0.093925,  -0.106139,  -0.117305,  -0.127369,  -0.136390,  -0.144379,  -0.151377,  -0.157434,
     -0.162577,  -0.166840,  -0.170234,  -0.172813,  -0.174620,  -0.175688,  -0.176023,  -0.175680,
     -0.174706,  -0.173106,  -0.170894,  -0.168202,  -0.165011,  -0.161307,  -0.157176,  -0.152607,
     -0.147557,  -0.142081,  -0.136239,  -0.129999,  -0.123437,  -0.116617,  -0.109574,  -0.102373,
     -0.095007,  -0.087537,  -0.080043,  -0.072616,  -0.065248,  -0.058088,  -0.051082,  -0.044310,
     -0.037795,  -0.031549,  -0.025613,  -0.020024,  -0.014806,  -0.009935,  -0.005428,  -0.001278,
      0.002482,   0.005907,   0.008959,   0.011615,   0.013870,   0.015787,   0.017339,   0.018517,
      0.019332,   0.019770,   0.019874,   0.019621,   0.019040,   0.018121,   0.016863,   0.015321,
      0.013509,   0.011401,   0.009026,   0.006425,   0.003626,   0.000544,  -0.002632,  -0.005983,
     -0.009427,  -0.013017,  -0.016689,  -0.020431,  -0.024234,  -0.028042,  -0.031836,  -0.035658,
     -0.039423,  -0.043136,  -0.046810,  -0.050338,  -0.053782,  -0.057091,  -0.060254,  -0.063274,
     -0.066109,  -0.068757,  -0.071207,  -0.073451,  -0.075482,  -0.077280,  -0.078852,  -0.080188,
     -0.081284,  -0.082127,  -0.082721,  -0.083097,  -0.083219,  -0.083109,  -0.082755,  -0.082168,
     -0.081339,  -0.080304,  -0.079066,  -0.077597,  -0.075941,  -0.074069,  -0.072021,  -0.069779,
     -0.067396,  -0.064853,  -0.062152,  -0.059315,  -0.056350,  -0.053291,  -0.050147,  -0.046936,
     -0.043704,  -0.040395,  -0.037138,  -0.033885,  -0.030660,  -0.027473,  -0.024368,  -0.021300,
     -0.018288,  -0.015407,  -0.012582,  -0.009858,  -0.007227,  -0.004706,  -0.002324,  -0.000031,
      0.002130,   0.004195,   0.006064,   0.007698,   0.009200,   0.010474,   0.011585,   0.012511,
      0.013271,   0.013848,   0.014257,   0.014497,   0.014572,   0.014496,   0.014254,   0.013872,
      0.013339,   0.012666,   0.011850,   0.011138,   0.010188,   0.009110,   0.008158,   0.006793,
      0.005538,   0.004206,   0.002612,   0.001258,  -0.000344,  -0.001992,  -0.003338,  -0.004921,
     -0.006434,  -0.007965,  -0.009713,  -0.011143,  -0.012771,  -0.014138,  -0.015626,  -0.017154,
     -0.018403,  -0.019971,  -0.021214,  -0.022396,  -0.023723,  -0.024601,  -0.025517,  -0.026438,
     -0.026992,  -0.027827,  -0.028283,  -0.028620,  -0.029106,  -0.028892,  -0.028978,  -0.028939,
     -0.028456,  -0.028155,  -0.027889,  -0.027506,  -0.026996,  -0.026378,  -0.025654,  -0.024651,
     -0.023704,  -0.022665,  -0.021541,  -0.020337,  -0.019068,  -0.017735,  -0.016342,  -0.014711,
     -0.013237,  -0.011738,  -0.010232,  -0.008693,  -0.007146,  -0.005621,  -0.003841,  -0.002284,
     -0.000776,   0.000724,   0.002214,   0.003644,   0.005050,   0.006943,   0.008500,   0.009509,
      0.010734,   0.011582,   0.012358,   0.013070,   0.013660,   0.014161,   0.014593,   0.014920,
      0.015144,   0.015294,   0.015353,   0.015321,   0.015224,   0.015049,   0.014793,   0.014488,
      0.014118,   0.013675,   0.013162,   0.012633,   0.012059,   0.011410,   0.010753,   0.010072,
      0.009353,   0.008613,   0.007871,   0.007135,   0.006358,   0.005583,   0.004854,   0.004142,
      0.003400,   0.002720,   0.002090,   0.001477,   0.000833,   0.000251,  -0.000275,  -0.000793,
     -0.001264,  -0.001688,  -0.002069,  -0.002376,  -0.002637,  -0.002845,  -0.003002,  -0.003141,
     -0.003198,  -0.003199,  -0.003160,  -0.003053,  -0.002947,  -0.002763,  -0.002559,  -0.002280,
     -0.001979,  -0.001663,  -0.001295,  -0.000888,  -0.000464,   0.000002,   0.000461,   0.000901,
      0.001400,   0.001907,   0.002405,   0.002921,   0.003388,   0.003905,   0.004401,   0.004919,
      0.005420,   0.005887,   0.006330,   0.006758,   0.007116,   0.007456,   0.007774,   0.008073,
      0.008310,   0.008506,   0.008645,   0.008722,   0.008751,   0.008732,   0.008692,   0.008552,
      0.008372,   0.008116,   0.007792,   0.007536,   0.007127,   0.006715,   0.006261,   0.005781,
      0.005441,   0.004954,   0.004442,   0.003903,   0.003344,   0.002783,   0.002483,   0.001928,
      0.001376,   0.000854,   0.000265,   0.000035,  -0.000537,  -0.001076,  -0.001992,  -0.002808,
     -0.002904,  -0.003105,  -0.003284,   0.000000};
    naxis[0] = Xncell;   /* Create/fill FArray with beam shape for interpolator */
    tFA = ObitFArrayCreate("TempFA", ndim, naxis);
    for (j=0; j<Xncell; j++) tFA->array[j] = Xbeam[j]*Xbeam[j];  /*  Voltage to power */
    in->myFI = newObitFInterpolateCreate ("Interp", tFA, in->myDesc, hwidth);
    tFA = ObitFArrayUnref(tFA);    /* Unreference temp FA */
  default:    /* Doh */
    return;
  }; /* end switch */
  in->doTab = TRUE;    /* Must be OK if it gets here */
} /* end FindTabBeam */

/**
 * Set MeerKAT tabulated beam
 * \param in     the BeamShape object
 */
static void MeerKATTabBeam (ObitBeamShape *in)
{
  olong i, j, itab, ndim=1, naxis[] = {1};
  ObitFArray *tFA = NULL;
  olong   ntab   = 1;            /* How many? */
  olong   hwidth = 2;            /* Interpolation half width */
  /*                    */
  odouble minFreq[] = { 800.0e6};  /* Lower bound of bands */
  odouble maxFreq[] = {2000.0e6};  /* Upper  bound of bands */
  odouble refFreq[] = {1350.0e6};  /* Tabulated reference freq */

  /* Have a tabulated beam? */
  in->doTab = FALSE;
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
  case 0:     /* 1.35 GHz beam */
    in->itabRefFreq = 1.0/refFreq[itab]; /* 1/tabulated ref freq */
    in->icellSize   = 3600.0/30.0;       /* 1/Tabulated cell spacing */
    /* 26 Jul 16 S MeerKAT 1.3-1.4 GHz, cell 30.0/3600, refFreq  = 1.35e9 */
    olong  MKncell  = 175;
    ofloat MKbeam[] = {    /* Fitted to voltage beam */
      1.000000,0.999920,0.999679,0.999278,0.998717,0.997996,0.997116,0.996075,0.994876,0.993518,
      0.992002,0.990328,0.988497,0.986509,0.984365,0.982066,0.979612,0.977006,0.974246,0.971334,
      0.968272,0.965059,0.961697,0.958187,0.954532,0.950731,0.946786,0.942698,0.938469,0.934100,
      0.929591,0.924946,0.920166,0.915252,0.910205,0.905028,0.899722,0.894288,.888729, 0.883045,
      0.877241,0.871318,0.865278,0.859118,0.852851,0.846466,0.839977,0.833374,0.826673,0.819866,
      0.812954,0.805952,0.798846,0.791648,0.784367,0.776988,0.769525,0.761978,0.754352,0.746645,
      0.738862,0.731004,0.723075,0.715081,0.707017,0.698888,0.690705,0.682462,0.674161,0.665812,
      0.657410,0.648964,0.640473,0.631937,0.623368,0.614762,0.606119,0.597448,0.588752,0.580035, 
      0.571291,0.562529,0.553751,0.544967,0.536166,0.527358,0.518550,0.509741,0.500926,0.492125,
      0.483325,0.474535,0.465762,0.457003,0.448264,0.439542,0.430849,0.422179,0.413538,0.404930,
      0.396351,0.387820,0.379324,0.370871,0.362465,0.354103,0.345796,0.337535,0.329336,0.321190,
      0.313106,0.305086,0.297133,0.289241,0.281423,0.273674,0.265988,0.258370,0.250832,0.243372,
      0.235992,0.228700,0.221495,0.214373,0.207335,0.200393,0.193529,0.186762,0.180058,0.173477,
      0.166996,0.160628,0.154340,0.148136,0.142154,0.136344,0.130664,0.125113,0.119700,0.114493,
      0.109589,0.105007,0.100781,0.096646,0.092559,0.088636,0.085108,0.081382,0.077876,0.073577,
      0.069953,0.066347,0.062494,0.058545,0.055008,0.051969,0.047982,0.044959,0.042908,0.041673,
      0.039321,0.038724,0.038647,0.039247,0.040620,0.042692,0.045443,0.048629,0.049098,0.052326,
      0.056978,0.052184,0.054663,0.056593,0.059868};
    naxis[0] = MKncell;   /* Create/fill FArray with beam shape for interpolator */
    tFA = ObitFArrayCreate("TempFA", ndim, naxis);
    for (j=0; j<MKncell; j++) tFA->array[j] = MKbeam[j]*MKbeam[j];  /*  Voltage to power */
    in->myFI = newObitFInterpolateCreate ("Interp", tFA, in->myDesc, hwidth);
    tFA = ObitFArrayUnref(tFA);    /* Unreference temp FA */
    break;
  default:    /* Doh */
    return;
  }; /* end switch */
  in->doTab = TRUE;    /* Must be OK if it gets here */
} /* end MeerKATabBeam */

/**
 * Set MeerKAT tabulated beam
 * Add Tabulated Power Beam for VLITE use on PBCorr only
 *  Might be valid for eVLA Pband PBCorr as well
 *  Should not be used for normalizing asymmetric beam correctiosn
 *  \param in   the BeamShape object
 */
static void FindVLITEBeam (ObitBeamShape *in)
{
   olong j, ndim=1, naxis[] = {1};
   ObitFArray *tFA = NULL;
   /* olong   ntab   = 1;            How many? */
   olong   hwidth = 2;            /* Interpolation half width */
  /*                    P       */
  /*odouble minFreq[] = {320.0e6};   Lower bound of bands */
  /*odouble maxFreq[] = {364.0e6};   Upper  bound of bands */
  odouble refFreq[] = {340.85e6};  /* Tabulated reference freq */
  in->itabRefFreq = 1.0/refFreq[0]; /* 1/tabulated ref freq */
  in->icellSize   = 3600.0/60.0;     /* 1/Tabulated cell spacing */
  olong  Pncell   = 305;
  /* 6 June 2016  tabulated entries at 60 arcsec from fits to VLITE Data  */
  ofloat Pbeam[]  = {     /* extrapolated beyond 4 degrees */
    1.000000, 0.999878, 0.999510, 0.998898, 0.998042, 0.996943, 0.995600, 0.994017, 
    0.992192, 0.990128, 0.987827, 0.985289, 0.982518, 0.979514, 0.976280, 0.972819, 
    0.969132, 0.965223, 0.961094, 0.956748, 0.952189, 0.947420, 0.942443, 0.937263, 
    0.931883, 0.926307, 0.920539, 0.914583, 0.908443, 0.902123, 0.895627, 0.888961,
    0.882128, 0.875133, 0.867981, 0.860677, 0.853225, 0.845630, 0.837898, 0.830033, 
    0.822040, 0.813925, 0.805693, 0.797349, 0.788898, 0.780345, 0.771696, 0.762956, 
    0.754130, 0.745223, 0.736242, 0.727190, 0.718074, 0.708899, 0.699669, 0.690391, 
    0.681068, 0.671707, 0.662312, 0.652889, 0.643442, 0.633977, 0.624498, 0.615009, 
    0.605517, 0.596025, 0.586539, 0.577061, 0.567598, 0.558154, 0.548732, 0.539337, 
    0.529973, 0.520216, 0.511762, 0.503513, 0.495461, 0.487594, 0.479860, 0.472211, 
    0.464602, 0.457008, 0.449412, 0.441815, 0.434249, 0.426702, 0.419148, 0.411588, 
    0.404046, 0.396562, 0.389173, 0.381904, 0.374788, 0.367862, 0.361166, 0.354740, 
    0.348596, 0.342691, 0.336999, 0.331508, 0.326200, 0.321041, 0.315993, 0.310997, 
    0.306006, 0.301010, 0.296015, 0.291027, 0.286066, 0.281138, 0.276259, 0.271462, 
    0.266785, 0.262279, 0.257999, 0.253991, 0.250247, 0.246743, 0.243471, 0.240415, 
    0.237547, 0.234840, 0.232247, 0.229724, 0.227250, 0.224829, 0.222460, 0.220124, 
    0.217801, 0.215473, 0.213142, 0.210840, 0.208615, 0.206543, 0.204695, 0.203071, 
    0.201650, 0.200431, 0.199443, 0.198684, 0.198089, 0.197593, 0.197155, 0.196772, 
    0.196469, 0.196265, 0.196165, 0.196158, 0.196211, 0.196283, 0.196370, 0.196456, 
    0.196508, 0.196493, 0.196381, 0.196145, 0.195762, 0.195209, 0.194474, 0.193558, 
    0.192453, 0.191135, 0.189614, 0.187899, 0.186036, 0.184078, 0.182076, 0.180053, 
    0.178012, 0.175960, 0.173921, 0.171933, 0.170017, 0.168171, 0.166373, 0.164585, 
    0.162761, 0.160912, 0.159070, 0.157270, 0.155546, 0.153908, 0.152345, 0.150844, 
    0.149398, 0.147996, 0.146621, 0.145267, 0.143926, 0.142593, 0.141266, 0.139946, 
    0.138640, 0.137347, 0.136058, 0.134761, 0.133446, 0.132103, 0.130722, 0.129305, 
    0.127858, 0.126383, 0.124879, 0.123340, 0.121757, 0.120120, 0.118419, 0.116645, 
    0.114804, 0.112903, 0.110950, 0.108962, 0.106957, 0.104948, 0.102953, 0.100984, 
    0.099052, 0.097169, 0.095345, 0.093592, 0.091916, 0.090322, 0.088818, 0.087409, 
    0.086097, 0.084879, 0.083751, 0.082710, 0.081753, 0.080875, 0.080070, 0.079334, 
    0.078660, 0.078042, 0.077474, 0.076950, 0.076461, 0.075996, 0.075544, 0.075095, 
    0.074642, 0.074184, 0.073720, 0.073249, 0.072772, 0.072287, 0.071795, 0.071295, 
    0.070788, 0.070276, 0.069758, 0.069235, 0.068709, 0.068179, 0.067645, 0.067108, 
    0.066569, 0.066028, 0.065484, 0.064940, 0.064394, 0.063847, 0.063301, 0.062754, 
    0.062207, 0.061660, 0.061113, 0.060566, 0.060020, 0.059473, 0.058926, 0.057779,
    0.057232, 0.056685, 0.056138, 0.055591, 0.055044, 0.054497, 0.053950, 0.053403, 
    0.052856, 0.052309, 0.051762, 0.051215, 0.050668, 0.050121, 0.049574, 0.049027,
    0.048480, 0.047933, 0.047386, 0.046839, 0.046292, 0.045745, 0.045198, 0.044651,
    0.044104, 0.043557, 0.043010, 0.042463, 0.041916, 0.041369, 0.040822, 0.040275,
    0.039728};
  naxis[0] = Pncell;   /* Create/fill FArray with beam shape for interpolator */
  tFA = ObitFArrayCreate("TempFA", ndim, naxis);
  for (j=0; j<Pncell; j++) tFA->array[j] = Pbeam[j];  
  in->myFI = newObitFInterpolateCreate ("Interp", tFA, in->myDesc, hwidth);
  tFA = ObitFArrayUnref(tFA);    /* Unreference temp FA */
    
  in->doVLITE = TRUE;    /* Must be OK if it gets here */
} /* end FindVLITEBeam */

/**
 * Check if a tabulated beam is available and if so enable it
 * Currently implemented: VLA S Band, P Band
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

/**
/ * MeerKAT beam
 *  Calculate cosine beam shape (Condon & Ransom, Essential Radio Astronomy eq 3.95)
 * Compute beamAng if <0
 * \param in  the BeamShape object
 * \param in  Angle from pointing in degrees
 * \return  beam power gain
 */
static ofloat GetMKBeam (ObitBeamShape *in, odouble Angle)
{
  ofloat gain=1.0, rhor, div;
  if (Angle<=0.0) return gain;
  if (in->beamAng<0.0) in->beamAng  = (57.5/60.0) * (1.5e9/in->refFreq);
  /*in->beamAng  = (57.5/60.0) * (1.5e9/in->refFreq);*/
  rhor = 1.18896*Angle/in->beamAng;
  div = (1.-4.*(rhor*rhor));
  if (fabs(div)<1.0e-5) div = 1.0e-5; /* Stop zero divides */
  gain = (cos(G_PI*rhor)/div);

  return gain*gain;
} /* end GetTabBeam */
