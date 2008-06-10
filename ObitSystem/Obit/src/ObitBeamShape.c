/* $Id: ObitBeamShape.c,v 1.1 2008/05/06 13:20:14 bcotton Exp $        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2008                                               */
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
  if (in->doJinc) gain = ObitPBUtilJinc(Angle, in->refFreq, in->antSize, in->pbmin);
  else            gain = ObitPBUtilPoly(Angle, in->refFreq, in->pbmin);
  return gain;
} /* end ObitBeamShapeGain */

/**
 * Calculate gain in a given offset from a symmetric beam shape.
 * Simple function of distance from pointing center.
 * \param in     the BeamShape object
 * \param Angle  Angular distance (deg) from pointing center
 * \return gain
 */
ofloat ObitBeamShapeGainSym (ObitBeamShape *in, odouble Angle)
{
  ofloat gain=1.0;

  /* Compute gain? */
  if (!in->doGain) return gain;

  /* Compute */
  if (in->doJinc) gain = ObitPBUtilJinc(Angle, in->refFreq, in->antSize, in->pbmin);
  else            gain = ObitPBUtilPoly(Angle, in->refFreq, in->pbmin);
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
  in->myDesc  = NULL;
  in->pbmin   = 0.0;
  in->antSize = 0.0;
  in->doGain  = FALSE;
  in->raPnt   = 0.0;
  in->decPnt  = 0.0;
  in->refFreq = 0.0;
  in->doJinc  = FALSE;
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

  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitBeamShapeClear */

