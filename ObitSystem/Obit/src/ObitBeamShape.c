/* $Id$        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2008,2015                                          */
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
 * Currently implemented: VLA S Band
 * \param in     the BeamShape object
 */
static void FindTabBeam (ObitBeamShape *in)
{
  olong i, j, itab, ndim=1, naxis[] = {1};
  ObitFArray *tFA = NULL;
  olong   ntab   = 2;            /* How many? */
  olong   hwidth = 2;            /* Interpolation half width */
  /*                    P      S */
  odouble minFreq[] = {200.0e6,   1.8e9, };  /* Lower bound of bands */
  odouble maxFreq[] = {500.0e6,   4.2e9, };  /* Upper  bound of bands */
  odouble refFreq[] = {225.875e6, 3.0e9, };  /* Tabulated reference freq */

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
    in->icellSize   = 3600.0/1000.0;     /* 1/Tabulated cell spacing */
    olong  Pncell  = 100;
    ofloat Pbeam[] = {    /* Fitted to beam */
         1.000000, 0.995992,  0.989139,  0.919887,  0.872861,  0.811718,  0.740684,
         0.666480, 0.585123,  0.498273,  0.414439,  0.334170,  0.257083,  0.183989,
	 0.110730, 0.037010, -0.031470, -0.092954, -0.135057, -0.165120, -0.179161,
	-0.188929,-0.191828, -0.190135, -0.184020, -0.174683, -0.163303, -0.149757,
	-0.136103,-0.122890, -0.110372, -0.098280, -0.087212, -0.077239, -0.069038,
	-0.063323,-0.059285, -0.056697, -0.055298, -0.055924, -0.056855, -0.058437,
	-0.060505,-0.062230, -0.063892, -0.065237, -0.065988, -0.065608, -0.063919,
	-0.061482,-0.058363, -0.055099, -0.050779, -0.046468, -0.041654, -0.036583,
	-0.031563,-0.026960, -0.022766, -0.019660, -0.017064, -0.015663, -0.014495,
	-0.014512,-0.015131, -0.016242, -0.017579, -0.019229, -0.020716, -0.022451,
	-0.023891,-0.024455, -0.024304, -0.022449, -0.020602, -0.016906, -0.013246,
	-0.009531,-0.005468, -0.001181,  0.002665,  0.006218,  0.009146,  0.011199,
	 0.012479, 0.013290,  0.013251,  0.012276,  0.010742,  0.009091,  0.007457,
	 0.005781, 0.004403,  0.003223,  0.002550,  0.002138,  0.002078,  0.002208,
	 0.002334, 0.002794};
    naxis[0] = Pncell;   /* Create/fill FArray with beam shape for interpolator */
    tFA = ObitFArrayCreate("TempFA", ndim, naxis);
    for (j=0; j<Pncell; j++) tFA->array[j] = Pbeam[j]*Pbeam[j];  /* Voltage to power */
    in->myFI = newObitFInterpolateCreate ("Interp", tFA, in->myDesc, hwidth);
    tFA = ObitFArrayUnref(tFA);    /* Unreference temp FA */
    break;
  case 1:     /* S band */
    in->itabRefFreq = 1.0/refFreq[itab]; /* 1/tabulated ref freq */
    in->icellSize   = 3600.0/30.0;       /* 1/Tabulated cell spacing */
    olong  Sncell  = 260;
    /* 30 Oct 2015, center S Band 260 entries  cell 30.0/3600, refFreq  = 3.0e9 */
    ofloat Sbeam[] = {    /* Fitted to voltage beam */
     1.000000,  0.998166,  0.992630,  0.983293,  0.970421,  0.954320,  0.935404,  0.913524,
     0.888822,  0.861374,  0.831416,  0.799036,  0.764535,  0.728045,  0.689913,  0.650275,
     0.609362,  0.567518,  0.524953,  0.481886,  0.438518,  0.395194,  0.352099,  0.308980,
     0.266747,  0.225901,  0.186161,  0.147768,  0.110880,  0.075563,  0.041967,  0.010247,
    -0.019428, -0.046835, -0.071868, -0.094454, -0.114670, -0.132424, -0.147865, -0.160953,
    -0.171729, -0.180223, -0.186506, -0.190641, -0.192712, -0.192838, -0.191127, -0.187693,
    -0.182687, -0.176257, -0.168558, -0.159757, -0.150029, -0.139540, -0.128448, -0.116952,
    -0.105134, -0.093179, -0.081205, -0.069350, -0.057769, -0.046559, -0.035830, -0.025724,
    -0.016359, -0.007807, -0.000134,  0.006612,  0.012402,  0.017222,  0.021069,  0.023921,
     0.025796,  0.026710,  0.026704,  0.025830,  0.024132,  0.021637,  0.018316,  0.014265,
     0.009485,  0.004062, -0.001861, -0.008130, -0.014574, -0.021038, -0.027397, -0.033664,
    -0.039820, -0.045824, -0.051647, -0.057171, -0.062322, -0.067032, -0.071234, -0.074890,
    -0.077957, -0.080561, -0.082522, -0.083701, -0.084332, -0.084191, -0.083518, -0.082184,
    -0.080141, -0.077564, -0.074378, -0.070700, -0.066562, -0.062019, -0.057122, -0.051870,
    -0.046327, -0.040691, -0.034717, -0.028573, -0.022528, -0.016300, -0.010206, -0.004643,
     0.000901,  0.005777,  0.010663,  0.015230,  0.019094,  0.022909,  0.026280,  0.028825,
     0.031233,  0.033134,  0.034181,  0.035043,  0.035392,  0.034940,  0.034293,  0.033187,
     0.031302,  0.029391,  0.027123,  0.024157,  0.021317,  0.017869,  0.014654,  0.011084,
     0.007436,  0.003847,  0.000021, -0.003935, -0.007603, -0.011327, -0.014734, -0.017899,
    -0.021275, -0.023962, -0.026598, -0.028873, -0.030939, -0.032493, -0.034100, -0.035061,
    -0.035822, -0.036300, -0.036455, -0.036332, -0.035764, -0.035022, -0.034104, -0.032824,
    -0.031269, -0.029762, -0.027853, -0.025813, -0.023907, -0.021643, -0.019305, -0.017294,
    -0.014936, -0.012591, -0.010618, -0.008327, -0.006108, -0.004353, -0.002350, -0.000847,
     0.000853,  0.002446,  0.003319,  0.004526,  0.005449,  0.005822,  0.006501,  0.006966,
     0.006886,  0.007098,  0.006748,  0.006236,  0.005676,  0.004979,  0.004272,  0.003450,
     0.002656,  0.001794,  0.000992,  0.000221, -0.000653, -0.001358, -0.002045, -0.002618,
    -0.003103, -0.003549, -0.003854, -0.004102, -0.004201, -0.004221, -0.004101, -0.003907,
    -0.003580, -0.003182, -0.002662, -0.002068, -0.001380, -0.000637,  0.000166,  0.001023,
     0.001935,  0.002858,  0.003828,  0.004772,  0.005751,  0.006665,  0.007557,  0.008394,
     0.009230,  0.009963,  0.010625,  0.011226,  0.011731,  0.012186,  0.012507,  0.012755,
     0.012876,  0.012943,  0.012859,  0.012656,  0.012367,  0.011903,  0.011415,  0.010884,
     0.010101,  0.009424,  0.008422,  0.007649,  0.006481,  0.005636,  0.004270,  0.003321,
     0.002067,  0.001079, -0.000353, -0.001331};
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
