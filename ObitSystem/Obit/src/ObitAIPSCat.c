/* $Id$  */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2012                                          */
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
#include <sys/types.h>
#include <time.h>
#include <math.h>
#include <string.h>
#include "ObitImageDesc.h"
#include "ObitAIPS.h"
#include "ObitAIPSCat.h"
#include "ObitTableDesc.h"
#include "ObitFile.h"
#include "ObitSystem.h"

/*-------- ObitIO: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitAIPSCat.c
 * ObitAIPSCat class function definitions.
 */

/*-----------------File Globals ---------------------------*/
/** number of bytes per "sector" in ancient aipsish */
static olong AIPS_NBPS = 256*sizeof(AIPSint); 

/*------------------ Structures -----------------------------*/
/** Structure to contain information about the structure of
 *  the AIPS header file.
 * Named beginning with a K are offsets with the
 * second character giving the data type:
 * \li I = integer, offset in type AIPSint,
 * \li H = Hollerith (character - not null terminated), 
 *         offset in type gfloat.
 * \li R = float, offset in type gfloat.
 * \li D = double, offset in type gdouble.
 */
typedef struct {
    /** H(2) Source name (8 char) */
    olong  KHOBJ;
    /** H(2)  Telescope, i.e., 'VLA' (8 char) */
    olong  KHTEL;
    /** H(2)  Instrument  e.g., receiver or correlator (8 char) */
    olong  KHINS;
    /** H(2)  Observer name *(8 char) */
    olong  KHOBS;
    /**  H(2)  Observation date in format 'YYYYMMDD'(8 char) */
    olong  KHDOB;
    /**  H(2)  Date map created in format 'YYYYMMDD'(8 char) */
    olong  KHDMP;
    /** H(2)   Map units, i.e., 'JY/BEAM ' (8 char) */
    olong  KHBUN;
    /** H(14)  Random Parameter types (8 char each) */
    olong  KHPTP;
    /** H(7)    Coordinate type, i.e., 'RA---SIN' */
    olong  KHCTP;
    /** R(7)    Coordinate value increment along axis */
    olong  KRCIC;
    /** R(7)    Coordinate Reference Pixel */
    olong  KRCRP;
    /** R(7)    Coordinate Rotation Angles */
    olong  KRCRT;
    /** R       Epoch of coordinates (years) (actually equinox) */
    olong  KREPO;
    /**  R      Float value of data maximum*/
    olong  KRDMX;
    /** R       Float value of data minimum */
    olong  KRDMN;
    /** R       Value of indeterminate pixel */
    olong  KRBLK;
    /** H       Image name (12 characters) */
    olong  KHIMN;
    /** H       Image class (6 characters) */
    olong  KHIMC;
    /** H       Map physical type (i.e., 'MA','UV') (2 char) */
    olong  KHPTY;
    /** R       Beam major axis in degrees */
    olong  KRBMJ;
    /** R       Beam minor axis in degrees */
    olong  KRBMN;
    /** R       Beam position angle in degrees */
    olong  KRBPA;
    /** R       Alternate ref pixel location (frequency or velocity) */
    olong  KRARP;
    /** R       Offset in X (rotated RA) of phase center */
    olong  KRXSH;
    /** R       Offset in Y (rotated Dec) from tangent pt. */
    olong  KRYSH;
    /**         Name Character offset in HOLLERITH string=1 */
    olong  KHIMNO;
    /**         Class Character offset in HOLLERITH string=13 */
    olong  KHIMCO;
    /**         Type Character offset in HOLLERITH=19 */
    olong  KHPTYO;
    /** D(7)    Coordinate value at reference pixel */
    olong  KDCRV;
    /** D       Antenna pointing Right Ascension */
    olong  KDORA;
    /** D       Antenna pointing Declination */
    olong  KDODE;
    /** D       Rest frequency of line (Hz) */
    olong  KDRST;
    /** D       Alternate ref pixel value (frequency or velocity)*/
    olong  KDARV;
    /** I       Number of random par. groups. (# visibilities) */
    olong  KIGCN;
    /** I       # clean iterations */
    olong  KINIT;
    /**         Max. number of labeled random parameters */
    olong  KIPTPN;
    /**         Max. number of axes */
    olong  KICTPN;
    /**         Max number of extension files */
    olong  KIEXTN;
    /** I       Number of random parameters */
    olong  KIPCN;
    /** I       Number of coordinate axes */
    olong  KIDIM;
    /** I(7)    Number of pixels on each axis */
    olong  KINAX;
    /** I       Image sequence no. */
    olong  KIIMS;
    /** I       Image user ID number */
    olong  KIIMU;
    /** I       Clean map type: 1-4 => normal,
     *  components, residual, points.  
     *  For uv data this word contains a two character sort order code. 
     */
    olong  KITYP;
    /** I       Velocity reference frame: 1-3
     *  => LSR, Helio, Observer +  256 if radio definition. 
     */
    olong  KIALT;
    /** H(?)   Names of extension file types (2 char) 1 per float(?) */
    olong  KHEXT;
    /** I(?)   Number of versions of corresponding extension file.*/
    olong  KIVER;
    /** UV weight normalization */
    olong  KRWTN;
    /** I 3D type, 1=>DO3D=FALSE; 2=>DO3D=TRUE */
    olong  KIITY;
    /**  R "X" pixel offset */
    olong  KRXPO;
    /**  R "Y" pixel offset */
    olong  KRYPO;
    /** Reserved  */
    olong  KIRES;
    /** Reserved */
  olong  KIRESN;
} ObitAIPSCatDHDR;

/* AIPS magic blanking value */
static union FBLANKequiv {
  gchar string[4];
  float fblank;
} FBLANK;

/*-----------------File Globals ---------------------------*/
/**
 * ObitAIPSCatDHDR global structure ObitClassInfo.
 * This structure is used to parse/write AIPS header files.
 */
static ObitAIPSCatDHDR myDHDR;

/*---------------Private function prototypes----------------*/
/** Private: Update access time  */
static void
ObitAIPSCatUpdateAccess(AIPSint *dateTime);

/** Private: Copy string replacing non printing chars with blanks  */
static void CopyDeeNull (gchar *out, gchar *in, olong nchar);

/*---------------Public functions---------------------------*/
/**
 * Converts a buffer filled with the contents of an AIPS header file
 * into an ObitImageDesc.
 * \param desc The Image descriptor to be updated.
 * \param buffer A buffer of AIPS_NBPS char with the contents of the AIPS
 *               header file.
 * \param err    ObitErr structure to report errors.
 */
void ObitAIPSCatImageGetDesc (ObitImageDesc *desc, gchar *buffer, 
			      ObitErr *err)
{
  olong   i, j, ndim;
  AIPSint *header  = (AIPSint*)buffer;
  ofloat  *fheader = (ofloat*)buffer;
  odouble *dheader = (odouble*)buffer;
  gchar   temp[9];

  /* error tests */
  g_assert (ObitIsA(desc, ObitImageDescGetClass()));
  g_assert(ObitErrIsA(err)); 
  g_assert(buffer!=NULL); 

  /* maximum number of dimensions to copy */
  ndim = IM_MAXDIM;

  /* bitpix */
  desc->bitpix = -32; /* AIPS only does floats */

  /* Number of axes. */
  desc->naxis = header[myDHDR.KIDIM];

  /* Dimensions of axes. */
  for (i=0; i<ndim; i++) desc->inaxes[i] = header[myDHDR.KINAX+i];

  /* WCS labels for each dimension of array. */
  for (i=0; i<ndim; i++) {
    g_memmove (desc->ctype[i],  &header[myDHDR.KHCTP+i*2], 8);
    desc->ctype[i][9] = 0; /* null terminate */
  }

  /* Axis coordinate increments. */
  for (i=0; i<ndim; i++) desc->cdelt[i] = fheader[myDHDR.KRCIC+i]; 

  /* Axis reference pixels (1-rel) */
  for (i=0; i<ndim; i++) desc->crpix[i] = fheader[myDHDR.KRCRP+i]; 

  /* Axis rotation angles (deg) */
  for (i=0; i<ndim; i++) desc->crota[i] = fheader[myDHDR.KRCRT+i]; 

  /* Axis coordinate values at reference pixel. */
  for (i=0; i<ndim; i++) desc->crval[i] = dheader[myDHDR.KDCRV+i]; 

  /* Name of object. */
  g_memmove (desc->object,  &header[myDHDR.KHOBJ], 8);
  desc->object[9] = 0; /* null terminate */

  /* Name of telescope making observation. */
  g_memmove (desc->teles,  &header[myDHDR.KHTEL], 8);
  desc->teles[9] = 0; /* null terminate */

  /* Name of instrument making observation. */
  g_memmove (desc->instrument,  &header[myDHDR.KHINS], 8);
  desc->instrument[9] = 0; /* null terminate */

  /* Name of observer. */
  g_memmove (desc->observer,  &header[myDHDR.KHOBS], 8);
  desc->observer[9] = 0; /* null terminate */

  /* Observing date as yyyy-mm-dd ,  AIPS is YYYYMMDD */
  g_memmove (temp,  &header[myDHDR.KHDOB], 8);
  j = 0;
  for (i=0; i<4; i++) desc->obsdat[i] = temp[j++]; 
  desc->obsdat[4] = '-';
  for (i=5; i<7; i++) desc->obsdat[i] = temp[j++];
  desc->obsdat[7] = '-';
  for (i=8; i<10; i++) desc->obsdat[i] = temp[j++];
  desc->obsdat[10] = 0; /* null terminate */

  /* Image creation date as yyyy-mm-dd,  AIPS is YYYYMMDD */
  g_memmove (temp,  &header[myDHDR.KHDMP], 8);
  j = 0;
  for (i=0; i<4; i++) desc->date[i] = temp[j++]; 
  desc->date[4] = '-';
  for (i=5; i<7; i++) desc->date[i] = temp[j++];
  desc->date[7] = '-';
  for (i=8; i<10; i++) desc->date[i] = temp[j++];
  desc->date[10] = 0; /* null terminate */

  /* Origin (software) of image. */
  g_memmove (desc->origin,  "Obit    ", 8);
  desc->origin[9] = 0; /* null terminate */

  /* Units of image */
  g_memmove (desc->bunit,  &header[myDHDR.KHBUN], 8);
  desc->bunit[9] = 0; /* null terminate */

  /* Epoch (years) of celestial coordinates,
   *  This is sometimes confused with equinox.
   */
  desc->epoch = fheader[myDHDR.KREPO];

  /* Mean Equinox of celestial coordinates (e.g. 2000) */
  desc->equinox = fheader[myDHDR.KREPO];

  /* Observed RA (deg) */
  desc->obsra = dheader[myDHDR.KDORA];

  /* Observed Dec (deg) */
  desc->obsdec = dheader[myDHDR.KDODE];

  /* maximum value in image */
  desc->maxval = fheader[myDHDR.KRDMX];

  /* minimum value in image */
  desc->minval = fheader[myDHDR.KRDMN];

  /* Are there blanked pixels in the image? */
  desc->areBlanks = (fheader[myDHDR.KRBLK] == FBLANK.fblank);

  /* x shift */
  desc->xshift = fheader[myDHDR.KRXSH];

  /* y shift */
  desc->yshift = fheader[myDHDR.KRYSH];

  /* Rest frequency of line (Hz) */
  desc->restFreq = dheader[myDHDR.KDRST];

  /* Alternate ref pixel location (frequency or velocity) */
  desc->altCrpix = fheader[myDHDR.KRARP];

  /* Alternate ref pixel value (frequency or velocity) */
  desc->altRef = dheader[myDHDR.KDARV];

  /* Beam major axis in degrees */
  desc->beamMaj = fheader[myDHDR.KRBMJ];
    
  /* Beam minor axis in degrees */
  desc->beamMin = fheader[myDHDR.KRBMN];

  /* 3D stuff - anything but explicit 2D is 3D */
  desc->do3D   = header[myDHDR.KIITY] != 1;
  desc->xPxOff = fheader[myDHDR.KRXPO];
  desc->yPxOff = fheader[myDHDR.KRYPO];
    
  /* Beam position angle in degrees */
  desc->beamPA = fheader[myDHDR.KRBPA];
    
  /* # clean iterations */
  desc->niter = header[myDHDR.KINIT];

  /* Velocity reference frame  */
  desc->VelDef = header[myDHDR.KIALT] / 256;
  desc->VelReference = header[myDHDR.KIALT] - 256*desc->VelDef;

} /* end ObitAIPSImageGetDesc */

/**
 * Converts an ObitImageDesc to entries in a buffer filled with the 
 * contents of an AIPS header file.
 * \param desc The Image descriptor to be updated.
 * \param buffer A buffer of AIPS_NBPS char with the contents of the AIPS
 *               header file.
 * \param init   If true, initialize the contents of buffer before
 *               converting.
 * \param dirEntry Structure with catalog directory with file name information.
 * \param err    ObitErr structure to report errors.
 */
void ObitAIPSCatImageSetDesc (ObitImageDesc *desc, gchar *buffer, 
			      gboolean init, ObitAIPSDirCatEntry* dirEntry,
			      ObitErr *err)
{
  olong   i, j, ndim;
  AIPSint *header  = (AIPSint*)buffer;
  ofloat  *fheader = (ofloat*)buffer;
  odouble *dheader = (odouble*)buffer;
  gchar   *cp, temp[13];
  gboolean doInit;
  gchar   *blank = "                    ";

  /* error tests */
  g_assert (ObitIsA(desc, ObitImageDescGetClass()));
  g_assert(ObitErrIsA(err)); 
  g_assert(buffer!=NULL); 

  /* Need to init header? init or old header has no axes defined */
  doInit = init || (header[myDHDR.KIDIM]<=0);

  /* maximum number of dimensions to copy */
  if (doInit)
    ndim = IM_MAXDIM;
  else
    ndim = MAX(IM_MAXDIM, header[myDHDR.KICTPN]);
  ndim = MIN (7, ndim);

  /* initialize first? */
  if (doInit) { /* yup */
    /* zero fill */
    for (i=0; i<AIPS_NBPS; i++) buffer[i] = 0;

    /* blank axis type values */
    for (i=0; i<header[myDHDR.KICTPN]; i++) 
      g_memmove (&header[myDHDR.KHCTP+i*2], blank,  8);

    /* blank random parameter values */
    for (i=0; i<header[myDHDR.KIPTPN]; i++) 
      g_memmove (&header[myDHDR.KHPTP+i*2], blank,  8);

    /* blank name class... */
    g_memmove (&header[myDHDR.KHIMN], blank, 20);

    /* blank observers name */
    g_memmove (&header[myDHDR.KHOBS], blank, 8);

    /* blank instrument name */
    g_memmove (&header[myDHDR.KHINS], blank, 8);

    /* Image sequence no. */
    header[myDHDR.KIIMS] = 0;
    
    /* Image user ID number */
    header[myDHDR.KIIMU] = 0;
    
    /* x shift */
    fheader[myDHDR.KRXSH] = 0.0;

    /* y shift */
    fheader[myDHDR.KRYSH] = 0.0;

    /* Rest frequency of line (Hz) */
    dheader[myDHDR.KDRST] = 0.0;

    /* Alternate ref pixel value (frequency or velocity) */
    dheader[myDHDR.KDARV] = 0.0;

    /* Alternate ref pixel location (frequency or velocity) */
    fheader[myDHDR.KRARP] = 0.0;

    /* Number of random parameters (# visibilities) */
    header[myDHDR.KIPCN] = 0;
    
    /* Beam major axis in degrees */
    fheader[myDHDR.KRBMJ] = 0.0;
    
    /* Beam minor axis in degrees */
    fheader[myDHDR.KRBMN] = 0.0;
    
    /* Beam position angle in degrees */
    header[myDHDR.KRBPA] = 0.0;
    
    /* # clean iterations */
    header[myDHDR.KINIT] = 0;
    
    /*  Clean map type */
    header[myDHDR.KITYP] = 1;
    
    /* Velocity reference frame  */
    header[myDHDR.KIALT] = 1;

    /* UV weight normalization */
    fheader[myDHDR.KRWTN] = 0.0;

    /* set extension files */
    for (i=0; i<myDHDR.KIEXTN; i++) {
      /* two character code of extension file. */
      g_memmove (&header[myDHDR.KHEXT+i], blank, 2);
      /* Number of versions of corresponding extension file. */
      header[myDHDR.KIVER+i] = 0;
    }
  } /* end of init */

  /* copy to header */

  /* 3D type   */
  if (desc->do3D) header[myDHDR.KIITY] = 2;
  else header[myDHDR.KIITY] = 1;
  
  /* Pixel offsets */
  fheader[myDHDR.KRXPO] = desc->xPxOff;
  fheader[myDHDR.KRYPO] = desc->yPxOff;

  /* AIPS naming info from catalog directory entry */
  cp = (gchar*)&header[myDHDR.KHIMN];
  CopyDeeNull (cp, dirEntry->name, 12);
  CopyDeeNull (cp+12, dirEntry->class, 6);
  CopyDeeNull (cp+18, dirEntry->type, 2);

  /* Image sequence no. */
  header[myDHDR.KIIMS] = dirEntry->seq;
  
  /* Image user ID number */
  header[myDHDR.KIIMU] = dirEntry->user; 
    
  /* Number of axes. */
  header[myDHDR.KIDIM] = desc->naxis;
    
  /* Dimensions of axes. */
  for (i=0; i<ndim; i++) header[myDHDR.KINAX+i] = desc->inaxes[i];
    
  /* WCS labels for each dimension of array. */
  for (i=0; i<ndim; i++) {
    CopyDeeNull ((gchar*)&header[myDHDR.KHCTP+i*2], desc->ctype[i],  8);
  }

  /* Axis coordinate increments. */
  for (i=0; i<ndim; i++) fheader[myDHDR.KRCIC+i] = desc->cdelt[i]; 

  /* Axis reference pixels (1-rel) */
  for (i=0; i<ndim; i++) fheader[myDHDR.KRCRP+i] = desc->crpix[i]; 

  /* Axis rotation angles (deg) */
  for (i=0; i<ndim; i++) fheader[myDHDR.KRCRT+i] = desc->crota[i]; 

  /* Axis coordinate values at reference pixel. */
  for (i=0; i<ndim; i++) dheader[myDHDR.KDCRV+i] = desc->crval[i]; 

  /* Name of object. */
  CopyDeeNull ((gchar*)&header[myDHDR.KHOBJ], desc->object,  8);

  /* Name of telescope making observation. */
  CopyDeeNull ((gchar*)&header[myDHDR.KHTEL], desc->teles,  8);

  /* Name of instrument making observation. */
  CopyDeeNull ((gchar*)&header[myDHDR.KHINS], desc->instrument,  8);

  /* Name of observer. */
  CopyDeeNull ((gchar*)&header[myDHDR.KHOBS], desc->observer,  8);

  /* Observing date as yyyy-mm-dd,  AIPS is YYYMMDD */
  j = 0;
  for (i=0; i<4; i++) temp[j++]  = desc->obsdat[i]; 
  for (i=5; i<7; i++) temp[j++]  = desc->obsdat[i];
  for (i=8; i<10; i++) temp[j++] = desc->obsdat[i];
  CopyDeeNull ((gchar*)&header[myDHDR.KHDOB], temp,  8);

  /* Image creation date as yyyy-mm-dd,  AIPS is YYYMMDD */
  j = 0;
  for (i=0; i<4; i++) temp[j++]  = desc->date[i]; 
  for (i=5; i<7; i++) temp[j++]  = desc->date[i];
  for (i=8; i<10; i++) temp[j++] = desc->date[i];
  CopyDeeNull ((gchar*)&header[myDHDR.KHDMP], temp,  8);

  /* Units of image */
  CopyDeeNull ((gchar*)&header[myDHDR.KHBUN], desc->bunit,  8);

  /* Epoch (years) of celestial coordinates,
   *  This is sometimes confused with equinox.
   */
   fheader[myDHDR.KREPO]= desc->epoch;

  /* Mean Equinox of celestial coordinates (e.g. 2000) */
   fheader[myDHDR.KREPO] = desc->equinox;

  /* Observed RA (deg) */
   dheader[myDHDR.KDORA] = desc->obsra;

  /* Observed Dec (deg) */
   dheader[myDHDR.KDODE] = desc->obsdec;

  /* maximum value in image */
   fheader[myDHDR.KRDMX] = desc->maxval;

  /* minimum value in image */
   fheader[myDHDR.KRDMN] = desc->minval;

   /* Are there blanked pixels in the image? */
   if (desc->areBlanks) {
     fheader[myDHDR.KRBLK] = ObitMagicF();
  } else
     fheader[myDHDR.KRBLK] = 0.0;
   
   /* x shift */
   fheader[myDHDR.KRXSH] = desc->xshift;

   /* y shift */
   fheader[myDHDR.KRYSH] = desc->yshift;

   /* Rest frequency of line (Hz) */
   dheader[myDHDR.KDRST] = desc->restFreq;

   /* Alternate ref pixel location (frequency or velocity) */
   fheader[myDHDR.KRARP] = desc->altCrpix;

   /* Alternate ref pixel value (frequency or velocity) */
   dheader[myDHDR.KDARV] = desc->altRef;

   /* Beam major axis in degrees */
   fheader[myDHDR.KRBMJ] = desc->beamMaj;
    
   /* Beam minor axis in degrees */
   fheader[myDHDR.KRBMN] = desc->beamMin;
    
   /* Beam position angle in degrees */
   fheader[myDHDR.KRBPA] = desc->beamPA;
    
   /* # clean iterations */
   header[myDHDR.KINIT] = desc->niter;

   /* Velocity reference frame  */
   header[myDHDR.KIALT] = desc->VelReference + 256*desc->VelDef;
} /* end ObitAIPSCatImageSetDesc */

/**
 * Converts a buffer filled with the contents of an AIPS header file
 * into an ObitUVDesc.
 * \param desc The UV descriptor to be updated.
 * \param buffer A buffer of AIPS_NBPS char with teh contents of the AIPS
 *               header file.
 * \param err    ObitErr structure to report errors.
 */
void ObitAIPSCatUVGetDesc (ObitUVDesc *desc, gchar *buffer, 
			      ObitErr *err)
{
  olong   i, j, ndim;
  AIPSint *header  = (AIPSint*)buffer;
  ofloat  *fheader = (ofloat*)buffer;
  odouble *dheader = (odouble*)buffer;
  gchar   temp[20];

  /* error tests */
  g_assert (ObitIsA(desc, ObitUVDescGetClass()));
  g_assert(ObitErrIsA(err)); 
  g_assert(buffer!=NULL); 

  /* maximum number of dimensions to copy */
  ndim = IM_MAXDIM;

  /* Number of axes. */
  desc->naxis = header[myDHDR.KIDIM];

  /* Dimensions of axes. */
  for (i=0; i<ndim; i++) desc->inaxes[i] = header[myDHDR.KINAX+i];

  /* WCS labels for each dimension of array. */
  for (i=0; i<ndim; i++) {
    g_memmove (desc->ctype[i],  &header[myDHDR.KHCTP+i*2], 8);
    desc->ctype[i][9] = 0; /* null terminate */
  }

  /* number of random parameters */
  desc->nrparm = header[myDHDR.KIPCN];
  desc->nrparm = MIN (14, desc->nrparm);

  /* number of visibilities */
  desc->nvis = header[myDHDR.KIGCN];

  /* WCS labels for random parameters. */
  for (i=0; i<desc->nrparm; i++) {
    g_memmove (desc->ptype[i],  &header[myDHDR.KHPTP+i*2], 8);
    desc->ctype[i][9] = 0; /* null terminate */
  }

  /*  Sort order */
  g_memmove (desc->isort, &header[myDHDR.KITYP], 2);
  desc->isort[2] = 0;
    
  /* Axis coordinate increments. */
  for (i=0; i<ndim; i++) desc->cdelt[i] = fheader[myDHDR.KRCIC+i]; 

  /* Axis reference pixels (1-rel) */
  for (i=0; i<ndim; i++) desc->crpix[i] = fheader[myDHDR.KRCRP+i]; 

  /* Axis rotation angles (deg) */
  for (i=0; i<ndim; i++) desc->crota[i] = fheader[myDHDR.KRCRT+i]; 

  /* Axis coordinate values at reference pixel. */
  for (i=0; i<ndim; i++) desc->crval[i] = dheader[myDHDR.KDCRV+i]; 

  /* Name of object. */
  g_memmove (desc->object,  &header[myDHDR.KHOBJ], 8);
  desc->object[9] = 0; /* null terminate */

  /* Name of telescope making observation. */
  g_memmove (desc->teles,  &header[myDHDR.KHTEL], 8);
  desc->teles[9] = 0; /* null terminate */

  /* Name of instrument making observation. */
  g_memmove (desc->instrument,  &header[myDHDR.KHINS], 8);
  desc->instrument[9] = 0; /* null terminate */

  /* Name of observer. */
  g_memmove (desc->observer,  &header[myDHDR.KHOBS], 8);
  desc->observer[9] = 0; /* null terminate */

  /* Observing date as yyyy-mm-dd ,  AIPS is YYYMMDD */
  g_memmove (temp,  &header[myDHDR.KHDOB], 8);
  j = 0;
  for (i=0; i<4; i++) desc->obsdat[i] = temp[j++]; 
  desc->obsdat[4] = '-';
  for (i=5; i<7; i++) desc->obsdat[i] = temp[j++];
  desc->obsdat[7] = '-';
  for (i=8; i<10; i++) desc->obsdat[i] = temp[j++];
  desc->obsdat[10] = 0; /* null terminate */

  /* Image creation date as yyyy-mm-dd,  AIPS is YYYMMDD */
  g_memmove (temp,  &header[myDHDR.KHDMP], 8);
  j = 0;
  for (i=0; i<4; i++) desc->date[i] = temp[j++]; 
  desc->date[4] = '-';
  for (i=5; i<7; i++) desc->date[i] = temp[j++];
  desc->date[7] = '-';
  for (i=8; i<10; i++) desc->date[i] = temp[j++];
  desc->date[10] = 0; /* null terminate */

  /* Origin (software) of image. */
  g_memmove (desc->origin,  "Obit    ", 8);
  desc->origin[9] = 0; /* null terminate */

  /* Units of image */
  g_memmove (desc->bunit,  &header[myDHDR.KHBUN], 8);
  desc->bunit[9] = 0; /* null terminate */

  /* Epoch (years) of celestial coordinates,
   *  This is sometimes confused with equinox.
   */
  desc->epoch = fheader[myDHDR.KREPO];

  /* Mean Equinox of celestial coordinates (e.g. 2000) */
  desc->equinox = fheader[myDHDR.KREPO];

  /* Observed RA (deg) */
  desc->obsra = dheader[myDHDR.KDORA];

  /* Observed Dec (deg) */
  desc->obsdec = dheader[myDHDR.KDODE];

  /* x shift */
  desc->xshift = fheader[myDHDR.KRXSH];

  /* y shift */
  desc->yshift = fheader[myDHDR.KRYSH];

  /* Rest frequency of line (Hz) */
  desc->restFreq = dheader[myDHDR.KDRST];

  /* Alternate ref pixel location (frequency or velocity) */
  desc->altCrpix = fheader[myDHDR.KRARP];

  /* Alternate ref pixel value (frequency or velocity) */
  desc->altRef = dheader[myDHDR.KDARV];

  /* Beam major axis in degrees */
  desc->beamMaj = fheader[myDHDR.KRBMJ];
    
  /* Beam minor axis in degrees */
  desc->beamMin = fheader[myDHDR.KRBMN];

  /* Beam position angle in degrees */
  desc->beamPA = fheader[myDHDR.KRBPA];
    
  /* Velocity reference frame  */
  desc->VelDef = header[myDHDR.KIALT] / 256;
  desc->VelReference = header[myDHDR.KIALT] - 256*desc->VelDef;
} /* end ObitAIPSUVGetDesc */

/**
 * Converts an ObitUVDesc to entries in a buffer filled with the 
 * contents of an AIPS header file.
 * \param desc The UV descriptor to be updated.
 * \param buffer A buffer of AIPS_NBPS char with the contents of the AIPS
 *               header file.
 * \param init   If true, initialize the contents of buffer before
 *               converting.
 * \param dirEntry Structure with catalog directory wile file name information.
 * \param err    ObitErr structure to report errors.
 */
void ObitAIPSCatUVSetDesc (ObitUVDesc *desc, gchar *buffer, 
			      gboolean init, ObitAIPSDirCatEntry* dirEntry,
			      ObitErr *err)
{
  olong   i, j, ndim;
  AIPSint *header  = (AIPSint*)buffer;
  ofloat  *fheader = (ofloat*)buffer;
  odouble *dheader = (odouble*)buffer;
  gchar   *cp, temp[13];
  gchar   *blank = "                    ";

  /* error tests */
  g_assert (ObitIsA(desc, ObitUVDescGetClass()));
  g_assert(ObitErrIsA(err)); 
  g_assert(buffer!=NULL); 

  /* maximum number of dimensions to copy */
  if (init)
    ndim = IM_MAXDIM;
  else
    ndim = MAX(IM_MAXDIM, header[myDHDR.KICTPN]);
  ndim = MIN (7, ndim);

  /* initialize first? */
  if (init) { /* yup */
    /* zero fill */
    for (i=0; i<AIPS_NBPS; i++) buffer[i] = 0;

    /* blank axis type values */
    for (i=0; i<header[myDHDR.KICTPN]; i++) 
      g_memmove (&header[myDHDR.KHCTP+i*2], blank,  8);

    /* blank random parameter values */
    for (i=0; i<header[myDHDR.KIPTPN]; i++) 
      g_memmove (&header[myDHDR.KHPTP+i*2], blank,  8);

    /* blank name class... */
    g_memmove (&header[myDHDR.KHIMN], blank, 20);

    /* blank observers name */
    g_memmove (&header[myDHDR.KHOBS], blank, 8);

    /* blank instrument name */
    g_memmove (&header[myDHDR.KHINS], blank, 8);

    /* Image sequence no. */
    header[myDHDR.KIIMS] = 0;
    
    /* Image user ID number */
    header[myDHDR.KIIMU] = 0;
    
    /* x shift */
    fheader[myDHDR.KRXSH] = 0.0;

    /* y shift */
    fheader[myDHDR.KRYSH] = 0.0;

    /* Rest frequency of line (Hz) */
    dheader[myDHDR.KDRST] = 0.0;

    /* Alternate ref pixel value (frequency or velocity) */
    dheader[myDHDR.KDARV] = 0.0;

    /* Alternate ref pixel location (frequency or velocity) */
    fheader[myDHDR.KRARP] = 0.0;

    /* Number of random parameters (# visibilities) */
    header[myDHDR.KIPCN] = 0;
    
    /* Beam major axis in degrees */
    fheader[myDHDR.KRBMJ] = 0.0;
    
    /* Beam minor axis in degrees */
    fheader[myDHDR.KRBMN] = 0.0;
    
    /* Beam position angle in degrees */
    header[myDHDR.KRBPA] = 0.0;
    
    /* # clean iterations */
    header[myDHDR.KINIT] = 0;
    
    /*  Sort order */
    g_memmove (&header[myDHDR.KITYP], blank, 2);
    
    /* Velocity reference frame  */
    header[myDHDR.KIALT] = 1;

    /* UV weight normalization */
    fheader[myDHDR.KRWTN] = 0.0;

    /* set extension files */
    for (i=0; i<myDHDR.KIEXTN; i++) {
      /* two character code of extension file. */
      g_memmove (&header[myDHDR.KHEXT+i], blank, 2);
      /* Number of versions of corresponding extension file. */
      header[myDHDR.KIVER+i] = 0;
    }
  } /* end of init */

  /* copy to header */

  /* AIPS naming info from catalog directory entry */
  cp = (gchar*)&header[myDHDR.KHIMN];
  CopyDeeNull (cp, dirEntry->name, 12);
  CopyDeeNull (cp+12, dirEntry->class, 6);
  CopyDeeNull (cp+18, dirEntry->type, 2);

  /* Image sequence no. */
  header[myDHDR.KIIMS] = dirEntry->seq;
  
  /* Image user ID number */
  header[myDHDR.KIIMU] = dirEntry->user; 
    
  /*  Sort order */
  g_memmove (&header[myDHDR.KITYP], desc->isort, 2);
    
  /* Number of axes. */
  header[myDHDR.KIDIM] = desc->naxis;

  /* Dimensions of axes. */
  for (i=0; i<ndim; i++) header[myDHDR.KINAX+i] = desc->inaxes[i];

  /* WCS labels for each dimension of array. */
  for (i=0; i<ndim; i++) {
    CopyDeeNull ((gchar*)&header[myDHDR.KHCTP+i*2], desc->ctype[i],  8);
  }

  /* number of random parameters */
  header[myDHDR.KIPCN] = desc->nrparm;

  /* number of visibilities */
  header[myDHDR.KIGCN] = desc->nvis;

  /* WCS labels for random partameters. */
  for (i=0; i<desc->nrparm; i++) {
    CopyDeeNull ((gchar*)&header[myDHDR.KHPTP+i*2], desc->ptype[i],  8);
  }

  /* Axis coordinate increments. */
  for (i=0; i<ndim; i++) fheader[myDHDR.KRCIC+i] = desc->cdelt[i]; 

  /* Axis reference pixels (1-rel) */
  for (i=0; i<ndim; i++) fheader[myDHDR.KRCRP+i] = desc->crpix[i]; 

  /* Axis rotation angles (deg) */
  for (i=0; i<ndim; i++) fheader[myDHDR.KRCRT+i] = desc->crota[i]; 

  /* Axis coordinate values at reference pixel. */
  for (i=0; i<ndim; i++) dheader[myDHDR.KDCRV+i] = desc->crval[i]; 

  /* Name of object. */
  CopyDeeNull ((gchar*)&header[myDHDR.KHOBJ], desc->object,  8);

  /* Name of telescope making observation. */
  CopyDeeNull ((gchar*)&header[myDHDR.KHTEL], desc->teles,  8);

  /* Name of instrument making observation. */
  CopyDeeNull ((gchar*)&header[myDHDR.KHINS], desc->instrument,  8);

  /* Name of observer. */
  CopyDeeNull ((gchar*)&header[myDHDR.KHOBS], desc->observer,  8);

  /* Observing date as yyyy-mm-dd,  AIPS is YYYMMDD */
  j = 0;
  for (i=0; i<4; i++) temp[j++]  = desc->obsdat[i]; 
  for (i=5; i<7; i++) temp[j++]  = desc->obsdat[i];
  for (i=8; i<10; i++) temp[j++] = desc->obsdat[i];
  CopyDeeNull ((gchar*)&header[myDHDR.KHDOB], temp,  8);

  /* Image creation date as yyyy-mm-dd,  AIPS is YYYMMDD */
  j = 0;
  for (i=0; i<4; i++) temp[j++]  = desc->date[i]; 
  for (i=5; i<7; i++) temp[j++]  = desc->date[i];
  for (i=8; i<10; i++) temp[j++] = desc->date[i];
  CopyDeeNull ((gchar*)&header[myDHDR.KHDMP], temp,  8);

  /* Units of image */
  CopyDeeNull ((gchar*)&header[myDHDR.KHBUN], desc->bunit,  8);

  /* Epoch (years) of celestial coordinates,
   *  This is sometimes confused with equinox.
   */
   fheader[myDHDR.KREPO]= desc->epoch;

  /* Mean Equinox of celestial coordinates (e.g. 2000) */
   fheader[myDHDR.KREPO] = desc->equinox;

  /* Observed RA (deg) */
   dheader[myDHDR.KDORA] = desc->obsra;

  /* Observed Dec (deg) */
   dheader[myDHDR.KDODE] = desc->obsdec;

   /* x shift */
   fheader[myDHDR.KRXSH] = desc->xshift;

   /* y shift */
   fheader[myDHDR.KRYSH] = desc->yshift;

   /* Rest frequency of line (Hz) */
   dheader[myDHDR.KDRST] = desc->restFreq;

   /* Alternate ref pixel location (frequency or velocity) */
   fheader[myDHDR.KRARP] = desc->altCrpix;

   /* Alternate ref pixel value (frequency or velocity) */
   dheader[myDHDR.KDARV] = desc->altRef;

   /* Beam major axis in degrees */
   fheader[myDHDR.KRBMJ] = desc->beamMaj;
    
   /* Beam minor axis in degrees */
   fheader[myDHDR.KRBMN] = desc->beamMin;
    
   /* Beam position angle in degrees */
   fheader[myDHDR.KRBPA] = desc->beamPA;
    
   /* Velocity reference frame  */
   header[myDHDR.KIALT] = desc->VelReference + 256*desc->VelDef;
} /* end ObitAIPSCatUVSetDesc */

/**
 * Copies table information from a buffer filled with the contents of an 
 * AIPS header file into an ObitTableList.
 * \param tableList The ObitTableList to be updated.
 *                  Entries not already present in tableList will be entered
 *                  with a NULL pointer for the Table.
 * \param buffer    A buffer of AIPS_NBPS char with the contents of the AIPS
 *                  header file.
 * \param user      AIPS user number.
 * \param disk      AIPS disk number.
 * \param cno       AIPS catalog slot number.
 * \param err       ObitErr structure to report errors.
 */
void ObitAIPSCatGetTable (ObitTableList *tableList, gchar *buffer, 
			  olong user, olong disk, olong cno, ObitErr *err)
{
  olong   j, k, ver;
  AIPSint *header  = (AIPSint*)buffer;
  gboolean exist;
  gchar name[8] = {"AIPS XX"}, *tableFileName = NULL;
  gchar *routine = "ObitAIPSCatGetTable";

  /* error tests */
  g_assert (ObitTableListIsA(tableList));
  g_assert(ObitErrIsA(err)); 
  if (err->error) return;  /* existing error condition */
  g_assert(buffer!=NULL); 

  /* loop over AIPS header */
  for (j = 0; j<myDHDR.KIEXTN; j++) {
    /* table name */
    g_memmove (&name[5], (gchar*)&header[myDHDR.KHEXT+j], 2);

    if (header[myDHDR.KIVER+j] > 0) { /* occupied? */
      /* Crazy value? */
      if (header[myDHDR.KIVER+j]>4096) header[myDHDR.KIVER+j] = 1;
      /* Loop over possibilities */
      for (k = 1; k<=header[myDHDR.KIVER+j]; k++) {
	/* not all things in AIPS header are tables 
	   - but pretend they are */
	/* Does it actually exist? */
	tableFileName =  
	  ObitAIPSFilename (OBIT_AIPS_Table, disk, cno, user, &name[5], k, err);
	exist = ObitFileExist(tableFileName, err);
	if (tableFileName) g_free (tableFileName);
	if (err->error) Obit_traceback_msg (err, routine, tableList->name);
	
	/* update TableList */
	if (exist ) {
	  ver = k;
	  ObitTableListPut (tableList, name, &ver, NULL, err);
	  if (err->error) Obit_traceback_msg (err, routine, tableList->name);
	}
      } /* end loop over possibilities */
    }
  } /* end loop over AIPS entries */
} /* end ObitAIPSCatGetTable */

/**
 * Copies an ObitTableList to entries in a buffer filled with the 
 * contents of an AIPS header file.
 * For AIPS tables, names must be of the form "AIPS XX" where XX is the
 * table type.
 * \param tableList The ObitTableList to be updated.
 * \param buffer    A buffer of AIPS_NBPS char with the contents of the AIPS
 *                  header file.
 * \param err       ObitErr structure to report errors.
 */
void ObitAIPSCatSetTable (ObitTableList *tableList, gchar *buffer, 
			  ObitErr *err)
{
  olong   i, j;
  AIPSint *header  = (AIPSint*)buffer;
  gchar   *name, *stemp;
  olong version, freeOne=-1;
  ObitTable *table;
  gboolean found;
  gchar *routine = "ObitAIPSCatSetTable";

  /* error tests */
  g_assert (ObitTableListIsA(tableList));
  g_assert(ObitErrIsA(err)); 
  if (err->error) return;  /* existing error condition */
  g_assert(buffer!=NULL); 

  /* Clear all previous table entries from Catalog header */
  for (j = 0; j<myDHDR.KIEXTN; j++) {
    stemp = (gchar*)&header[myDHDR.KHEXT+j];
    stemp[0] = ' '; stemp[1] = ' '; 
    header[myDHDR.KIVER+j] = 0;
  }

  /* Loop over entries in ObitTableList */
  for (i=1; i<=tableList->number; i++) {
    found = FALSE;
    freeOne=-1;
    ObitTableListGetNumber(tableList, i, &name, &version, &table, err);
    if (err->error) Obit_traceback_msg (err, routine, tableList->name);
    /* This a valid table? */
    if ((ObitTableIsA(table) && (table->ReferenceCount<1)) || 
	(version<=-9)) continue;  
    table = ObitUnref(table); /* Don't need */

    /* Ignore "History" tables */
    if (!strncmp (name, "History", 7)) {
      ObitTableListRemove (tableList, name, version); /* Remove */
      continue;
    }

    /* name better be in AIPSish */
    if (strncmp (name, "AIPS ", 5)) {
	Obit_log_error(err, OBIT_Error, 
		       "Table name %s not in AIPSish for %s", 
		       name, tableList->name);
	if (name) g_free(name); name = NULL;   /* clean up */
	return;
    }
    
    /* See if it already exists in header */
    for (j = 0; j<myDHDR.KIEXTN; j++) {
      /* Check for the first free header */
      if ((header[myDHDR.KIVER+j]<=0) &&(freeOne<0)) freeOne = j;

      /* does this match? */
      if (!strncmp (&name[5], (gchar*)&header[myDHDR.KHEXT+j], 2)) {
	header[myDHDR.KIVER+j] = MAX (version, header[myDHDR.KIVER+j]);
	found = TRUE;
	break;
      }
    }
    if (!found) { /* not found, add new one */
      /* make sure there is room */
      if (freeOne < 0) {
	Obit_log_error(err, OBIT_Error, 
		       "No space for more Tables in AIPS header for %s", 
		       tableList->name);
	if (name) g_free(name); name = NULL;  /* clean up */
	return;
      }

      /* make entry in header */
      header[myDHDR.KIVER+freeOne] = version;
      CopyDeeNull ((gchar*)&header[myDHDR.KHEXT+freeOne], &name[5], 2);
    }

    if (name) g_free(name); name = NULL;   /* clean up */
  } /* end loop over ObitTableList */

} /* end ObitAIPSCatSetTable */

/**
 * Converts a AIPS table header records into an ObitTableDesc.
 * Column labels, units and any keyword value pairs are not dealt with;
 * this routine needs to be coordinated with 
 * #ObitIOTableAIPS::ObitIOTableAIPSReadDescriptor.
 * \param desc         The Table descriptor to be updated.
 * \param tabType      AIPS two character type code
 * \param tabVer       AIPS Table version.
 * \param controlBlock AIPS table control block.
 * \param record       Second "record" with size and type info.
 * \param err          ObitErr structure to report errors.
 */
void 
ObitAIPSCatTableGetDesc (ObitTableDesc *desc, 
			 gchar tabType[3], olong tabVer,
			 AIPSint controlBlock[256],
			 AIPSint record[256], ObitErr *err)
{
  olong i, j, itemp, atype, alen, physical, off;
  ObitInfoType otype;

  /* error tests */
  g_assert (ObitIsA(desc, ObitTableDescGetClass()));
  g_assert(ObitErrIsA(err)); 
  g_assert(controlBlock!=NULL); 
  g_assert(record!=NULL); 

  if (desc->dim)
    for (i=0; i<desc->nfield; i++) {
      g_free(desc->dim[i]); desc->dim[i] = NULL;
    }

  /* Number of rows */
  desc->nrow = controlBlock[4];

  /* Number of columns */
  desc->nfield = controlBlock[9] + 1;

  /* beginning of row data in bytes */
  desc->startData = AIPS_NBPS * (controlBlock[49]-1);

  /* Sort order  */
  /* Going AIPS documentation appears to be in error */
  itemp = abs (controlBlock[42]) % 256;
  if (abs (controlBlock[42])>256) itemp +=256; /* Abs value? */
  desc->sort[0] = itemp;
  if (controlBlock[42] < 0) desc->sort[0] = -desc->sort[0];
  itemp = abs (controlBlock[43]) % 256;
  if (abs (controlBlock[43])>256) itemp +=256; /* Abs value? */
  desc->sort[1] = itemp;
  if (controlBlock[43] < 0) desc->sort[1] = -desc->sort[1];

  /* Table name - this was totally botched in AIPS - must use "AIPS XX" */
  if (desc->TableName) g_free(desc->TableName);
  desc->TableName = g_strndup ("AIPS XX", 7);
  desc->TableName[5] = tabType[0];
  desc->TableName[6] = tabType[1];

  /* table version number */
  desc->version = tabVer;

  /* Initialize structure arrays */
  desc->offset = g_realloc(desc->offset, desc->nfield*sizeof(olong));
  desc->order  = g_realloc(desc->order, desc->nfield*sizeof(olong));
  desc->type   = g_realloc(desc->type, desc->nfield*sizeof(ObitInfoType));
  desc->dim    = g_realloc(desc->dim, desc->nfield*sizeof(gint32*));
  for (i=0; i<desc->nfield; i++) {
    desc->dim[i] = g_malloc0(MAXINFOELEMDIM*sizeof(gint32));
  }

  /* parse types and sizes array to get structure of table */
  for (i=0; i<desc->nfield-1; i++) {

    /* Parse type and length */
    atype = record[128+i];
    alen  = atype / 10;
    atype = atype - 10 * alen;

    /* Convert to Obit type */
    otype = OBIT_oint;
    if (atype==1)      otype = OBIT_double;
    else if (atype==2) otype = OBIT_float;
    else if (atype==3) otype = OBIT_string;
    else if (atype==4) otype = OBIT_oint;
    else if (atype==5) otype = OBIT_bool;
    else if (atype==7) otype = OBIT_bits;
    else if (atype==9) otype = OBIT_oint;

    /* AIPS has real problems with character strings, 
       always grouped in multiples of 4 */
    if (atype==3) alen = 4 + 4*((alen-1) / 4);

    desc->type[i]   = otype;
    desc->dim[i][0] = alen;
    /* AIPS can only do 1-D */
    for (j=1; j<MAXINFOELEMDIM; j++) desc->dim[i][j] = 1;
   }

  /* Add "_status" column - this corresponds to the AIPS select column
     but with slightly different meanings. */
  desc->type[desc->nfield-1] = OBIT_int;
  for (j=0; j<MAXINFOELEMDIM; j++) desc->dim[desc->nfield-1][j] = 1;

  /* index */
  ObitTableDescIndex (desc);

  /* sanity check to be sure AIPS notion of structure is the same as Obit's */
  /* Check correspondance between logical and physical */
  for (i=0; i<desc->nfield; i++) {
    physical = controlBlock[128+i];
    /* test Obit Offset against AIPS */
    if (desc->order[i] != controlBlock[128+i]) {
      Obit_log_error(err, OBIT_Error, 
          "Obit and AIPS table order inconsistency for column %d: %d %d", 
	   i+1, desc->order[i], controlBlock[128+i]);
    }
  } /* End of consistency check */

  /* Check offsets in data records */
  for (i=0; i<desc->nfield; i++) {
    /* String are special - stored in floats in AIPS */
    if (desc->type[i]==OBIT_string) off = 1+desc->offset[i]/4;
    else off = 1+desc->offset[i];
    /* test Obit Offset against AIPS */
    if (off != record[i]) {
      Obit_log_error(err, OBIT_Error, 
          "Obit and AIPS table offset inconsistency for column %d: %d %d", 
	   i+1, 1+off, record[i]);
    }
  } /* End of consistency check */

  /* If Errors tell about table */
  if (err->error)
    Obit_log_error(err, OBIT_Error, 
		   "Errors in AIPS table %s version %d for %s", 
		   tabType, tabVer, desc->name);
  

} /* end ObitAIPSCatTableGetDesc */

/**
 * Converts an ObitTableDesc to entries in buffers filled with the 
 * contents of an AIPS table header.
 * \param desc The Table descriptor to be updated.
 * \param init   If true, initialize the contents of buffers before
 *               converting.
 * \param tabType      AIPS two character type code
 * \param tabVer       AIPS Table version.
 * \param controlBlock AIPS table control block.
 * \param record       Second "record" with size and type info.
 * \param err    ObitErr structure to report errors.
 */
void ObitAIPSCatTableSetDesc (ObitTableDesc *desc, gboolean init, 
			      gchar tabType[3], olong tabVer,
			      AIPSint controlBlock[256],
			      AIPSint record[256], ObitErr *err)
{
  olong i, len, more, itemp, ll;
  olong nrps, nspr=0, ns, nr, off;
  olong AIPStype;
  gchar  tt[4], *pgmName;
  gchar   *blank = "                                                  ";
  gchar *routine = " ObitAIPSCatTableSetDesc";

  /* error tests */
  g_assert (ObitIsA(desc, ObitTableDescGetClass()));
  g_assert(ObitErrIsA(err)); 
  g_assert(controlBlock!=NULL); 
  g_assert(record!=NULL); 
  /* Has descriptor been completed? */
  Obit_return_if_fail((desc->lrowIO > 0), err, 
		      "%s: Table Descriptor not fully defined", routine);

  /* initialize first? */
  if (init) { /* yup */
    /* zero fill blocks */
    for (i=0; i<256; i++) controlBlock[i] = 0;
    for (i=0; i<256; i++) record[i] = 0;

    /* Strings in controlBlock */
    g_memmove ((gchar*)&controlBlock[16], blank, 48);
    g_memmove ((gchar*)&controlBlock[28], blank,  8);
    g_memmove ((gchar*)&controlBlock[38], blank, 48);
    g_memmove ((gchar*)&controlBlock[53], "*AIPS TABLE*", 12);
    g_memmove ((gchar*)&controlBlock[100], blank, 48);
    g_memmove ((gchar*)&controlBlock[112], blank, 48);
    g_memmove ((gchar*)&controlBlock[124], blank, 16);

    /* Fill in title = AIPS XX */
    g_memmove ((gchar*)&controlBlock[100], "AIPS", 4);
    tt[0] = ' '; tt[1] = ' '; tt[2] = tabType[0]; tt[3] = tabType[1];
    g_memmove ((gchar*)&controlBlock[101], tt, 4);

    /* Creator info */
    ObitAIPSCatUpdateAccess(&controlBlock[10]); /* date/time */
    pgmName = ObitSystemGetPgmName ();
    ll = MIN (12, strlen(pgmName));
    g_memmove ((gchar*)&controlBlock[28], pgmName, ll); /* Creator */

    /* Number of logical records for expansion */
    controlBlock[41] = 100;

    /* "Record" for column types and data pointers */
    controlBlock[44] = 2;

    /* "Record" for row selection strings  */
    controlBlock[45] = 3;

    /* "Record" for beginning of titles  */
    controlBlock[46] = 5;

    /* "Record" for beginning of units  */
    controlBlock[47] = controlBlock[46] + 1 + ((desc->nfield-2)/(256/6));

    /* "Record" for beginning of keywords  */
    controlBlock[48] = controlBlock[47] + 1 + ((desc->nfield-2)/(256/2));

    /* "Record" for beginning of data  */
    controlBlock[49] = controlBlock[48] + 1 + ((desc->info->number-1)/(256/5));

    /* What is this */
    controlBlock[50] = 0;

    /* Maximum number of keywords allowed */
    /* how many more will fit? */
    len = 1 + ((desc->info->number-1)/(256/5)); /* how many blocks allocated? */
    len = len * (256/5); /* convert to number of keywords */
    controlBlock[51] = len;

    /* FITS ASCII? */
    controlBlock[59] = 1;

    /* Selection strings? */
    controlBlock[60] = 0;

    /* Next availableSelection strings */
    controlBlock[61] = 1;

    /* Lookup table */
    for (i=0; i<desc->nfield; i++) {
      controlBlock[128+i] = desc->order[i];
    }

    /* Data offsets */
    for (i=0; i<desc->nfield; i++) {
      off = 1+desc->offset[i];
      /* Strings are always different */
      if (desc->type[i]==OBIT_string) off = 1 + (desc->offset[i]/4);
      record[i] = off;
    }

    /* Data types and sizes */
    for (i=0; i<desc->nfield; i++) {
      AIPStype = 4;
      if (desc->type[i]==OBIT_double) AIPStype = 1;
      if (desc->type[i]==OBIT_float)  AIPStype = 2;
      if (desc->type[i]==OBIT_string) AIPStype = 3;
      if (desc->type[i]==OBIT_oint)   AIPStype = 4;
      if (desc->type[i]==OBIT_long)   AIPStype = 4;
      if (desc->type[i]==OBIT_long)    AIPStype = 4;
      if (desc->type[i]==OBIT_bool)   AIPStype = 5;
      if (desc->type[i]==OBIT_bits)   AIPStype = 7;
      if (i==desc->nfield-1) AIPStype = 9;  /* selection */
      record[128+i] = AIPStype + 10 * desc->repeat[i];
    } /* end loop over fields */
    
    /* beginning of row data in bytes */
    desc->startData = AIPS_NBPS * (controlBlock[49]-1);
  } /* end of init section */

  /* Update things that might have changed */

  /* how many sectors per row? */
  nrps = AIPS_NBPS / desc->lrowIO; /* # rows per sector */
  if (nrps > 0) { /* multiple rows per sector */
    /* How many whole sectors */
    ns = desc->nrow / nrps;
    /* How many rows in this sector? */
    nr = desc->nrow % nrps;
  } else { /* multiple sectors per row */
    /* Number of sectors per row */
    nspr = 1 + (desc->lrowIO-1) / AIPS_NBPS;
    /* How many whole sectors */
    ns = desc->nrow * nspr;
    /* How many rows in this sector? */
    nr = 0;
  }

  /* Number of 256 AIPSint blocks in file */
  len = (desc->startData / (AIPS_NBPS)) + ns;
  if (nr>0) len++;
  controlBlock[0] = len;

  /* Maximum number of rows */
  /* how many more will fit? */
  more = 0;
  if (nrps > 0) more = nrps - nr;
  controlBlock[2] = more + desc->nrow;

  /* Current number */
  controlBlock[4] = desc->nrow;  
  /* Number of AIPSints per row including selection */
  controlBlock[7] = desc->lrowIO / sizeof(AIPSint);

  /* Rows per sector or sectors per row */
  if (desc->lrow<AIPS_NBPS) 
    controlBlock[8] = nrps;
  else /* sectors per row */
    controlBlock[8] = -nspr;

  /* number of columns per row not including selection */
  controlBlock[9] = desc->nfield - 1;

  /* Last access */
  ObitAIPSCatUpdateAccess(&controlBlock[32]);
  pgmName = ObitSystemGetPgmName ();
  ll = MIN (12, strlen(pgmName));
  g_memmove ((gchar*)&controlBlock[38], pgmName, ll); /* Creator */

  /* sort order (logical column numbers */
  itemp = abs (desc->sort[0]) % 256;
  if (abs (desc->sort[0])>256) itemp += 256;/* Abs value? */
  if (desc->sort[0]>0) controlBlock[42] = itemp;
  else controlBlock[42] = -itemp;
  /* second key */
  itemp = abs (desc->sort[1]) % 256;
  if (abs (desc->sort[1])>256) itemp += 256;/* Abs value? */
  if (desc->sort[1]>0) controlBlock[43] = itemp;
  else controlBlock[43] = -itemp;

  /* how many keywords now exist? */
  controlBlock[52] = desc->info->number;
   
} /* end ObitAIPSCatTableSetDesc */

/**
 * Sets values about sizes and locations of pieces in an AIPS
 * header file.
 * Uses 0 relative values rather than AIPS 1-rel.
 * This is more-or-less AIPS routine $APLSUB/VHDRIN.FOR 
 */
void ObitAIPSCatInitDHDR(void) {
  olong PI, PH, PR, PD, NWDPDP, INC4C8=2;
  ofloat X, EPS=0.01;

  /* Make sure AIPSint and ofloat the same size */
  g_assert(sizeof(AIPSint)==sizeof(ofloat));

  /* initialize AIPS magic blanking value float equiv of 'INDE' */
  FBLANK.fblank = ObitMagicF();
  
  /* number of words (float) per double */
  NWDPDP = sizeof(double) / sizeof(float);

  /* set address increments */
  PH = 1;

  /* array sizes */
  myDHDR.KIPTPN = 14;
  myDHDR.KICTPN = 7;
  myDHDR.KIEXTN = 50;/* from PHDR.INC */

  /* character of 8 */
  myDHDR.KHOBJ = PH - 1;
  myDHDR.KHTEL = PH +   INC4C8 - 1;
  myDHDR.KHINS = PH + 2*INC4C8 - 1;
  myDHDR.KHOBS = PH + 3*INC4C8 - 1;
  myDHDR.KHDOB = PH + 4*INC4C8 - 1;
  myDHDR.KHDMP = PH + 5*INC4C8 - 1;
  myDHDR.KHBUN = PH + 6*INC4C8 - 1;
  myDHDR.KHPTP = PH + 7*INC4C8 - 1;
  myDHDR.KHCTP = PH + (7+myDHDR.KIPTPN)*INC4C8 - 1;
  PH = PH + INC4C8 * (7 + myDHDR.KIPTPN + myDHDR.KICTPN);
  PI = PH;
  
  /* reals: double precision */
  X = (PI - 1.0) / NWDPDP + 1.0;
  PD = X + EPS;
  if (ABS(PD-X)> EPS) {
    PD = X + 1.0 - EPS;
    PI = (PD-1)*NWDPDP + 1;
  }
  myDHDR.KDCRV = PD - 1;
  PI = PI + NWDPDP * myDHDR.KICTPN;

  /* reals: single precision */
  PR = PI;
  myDHDR.KRCIC = PR - 1;
  myDHDR.KRCRP = PR +   myDHDR.KICTPN - 1;
  myDHDR.KRCRT = PR + 2*myDHDR.KICTPN - 1;
  PR = PR + 3 * myDHDR.KICTPN;
  myDHDR.KREPO = PR - 1;
  myDHDR.KRDMX = PR + 1.0 - 1;
  myDHDR.KRDMN = PR + 2.0 - 1;
  myDHDR.KRBLK = PR + 3.0 - 1;
  PR = PR + 4;
  PI = PR;

  /* GCN now I, other integers */
  myDHDR.KIGCN = PI - 1;
  myDHDR.KIPCN = PI + 1 - 1;
  myDHDR.KIDIM = PI + 2 - 1;
  myDHDR.KINAX = PI + 3 - 1;
  myDHDR.KIIMS = PI + (3+myDHDR.KICTPN) - 1;
  PI = PI + (4+myDHDR.KICTPN);

  /* name \\ class \\ phystype */
  /* character(20) */
  PH = PI;
  myDHDR.KHIMN = PH - 1;
  myDHDR.KHIMC = myDHDR.KHIMN - 1;
  myDHDR.KHPTY = myDHDR.KHIMN - 1;
  myDHDR.KHIMNO = 1 - 1;
  myDHDR.KHIMCO = 13 - 1;
  myDHDR.KHPTYO = 19 - 1;
  PI = PH + 5;

  /* user #, clean parms */
  myDHDR.KIIMU = PI - 1;
  myDHDR.KIALT = PI + 1 - 1;
  PI = PI + 2;
  PR = PI;
  myDHDR.KRBMJ = PR - 1;
  myDHDR.KRBMN = PR + 1 - 1;
  myDHDR.KRBPA = PR + 2 - 1;
  PI = PI + 3;

  /* iterations */
  myDHDR.KINIT = PI - 1;
  myDHDR.KITYP = PI + 1 - 1;
  PI = PI + 2;

  /* Pointing position, alt freq */
  X = (PI - 1.0) / NWDPDP + 1.0;
  PD = X + EPS;
  if (fabs(X-PD) > EPS) {
    PD = X + 1.0 - EPS;
    PI = (PD-1) * NWDPDP + 1;
  }
  myDHDR.KDORA = PD - 1;
  myDHDR.KDODE = PD + 1 - 1;
  myDHDR.KDRST = PD + 2 - 1;
  myDHDR.KDARV = PD + 3 - 1;
  PD = PD + 4;
  PI = (PD - 1) * NWDPDP + 1;
  PR = PI;
  myDHDR.KRARP = PR - 1;
  myDHDR.KRXSH = PR + 1 - 1;
  myDHDR.KRYSH = PR + 2 - 1;
  PI = PI + 3;

  /* extension files */
  myDHDR.KHEXT = PI - 1;
  myDHDR.KIVER = PI + myDHDR.KIEXTN - 1;
  PI = PI + 2*myDHDR.KIEXTN;
  /* G&C wcs removed from main */
  /* header */
  /* myDHDR.KRCOK = PI */
  /* myDHDR.KDLON = PI / NWDPDP + 1 */
  /* myDHDR.KDPRJ = myDHDR.KDLON + 1 */
  /* PD = myDHDR.KDPRJ + 9 */
  /* PI = (PD-1) * NWDPDP + 1 */
  /* myDHDR.KRPCM = PI */
  /* UV weight normalization */
  myDHDR.KRWTN = PI;
  /* Imaging type, xpoff, ypoff */
  myDHDR.KIITY = PI + 1;
  myDHDR.KRXPO = PI + 2;
  myDHDR.KRYPO = PI + 3;

  /* residual space */
  myDHDR.KIRES = PI + 4;
  myDHDR.KIRESN = 257 - myDHDR.KIRES;

} /* end  ObitAIPSCatInitDHDR */

/**
 * Returns offset in AIPS catalog header for given keywords 
 * \param keyword name of keyword.
 * \return 0-rel offset in AIPS catalog header, -1 => not found
 */
olong ObitAIPSCatOffset (gchar *keyword)
{
  olong ncomp, off = 0;

  ncomp = MIN (8, strlen(keyword));

  if (!strncmp(keyword, "OBJECT  ", ncomp))      off = myDHDR.KHOBJ;
  else if (!strncmp(keyword, "TELESCOP", ncomp)) off = myDHDR.KHTEL;
  else if (!strncmp(keyword, "INSTRUME", ncomp)) off = myDHDR.KHINS;
  else if (!strncmp(keyword, "OBSERVER", ncomp)) off = myDHDR.KHOBS;
  else if (!strncmp(keyword, "DATE-OBS", ncomp)) off = myDHDR.KHDOB;
  else if (!strncmp(keyword, "DATE-MAP", ncomp)) off = myDHDR.KHDMP;
  else if (!strncmp(keyword, "BUNIT   ", ncomp)) off = myDHDR.KHBUN;
  else if (!strncmp(keyword, "SORTORD ", ncomp)) off = myDHDR.KITYP;
  else if (!strncmp(keyword, "NDIM    ", ncomp)) off = myDHDR.KIDIM;
  else if (!strncmp(keyword, "NAXIS   ", ncomp)) off = myDHDR.KINAX;
  else if (!strncmp(keyword, "GCOUNT  ", ncomp)) off = myDHDR.KIGCN;
  else if (!strncmp(keyword, "NRPARM  ", ncomp)) off = myDHDR.KIPCN;
  else if (!strncmp(keyword, "TYPEUVD ", ncomp)) off = 3; /* UV data type */
  else if (!strncmp(keyword, "CTYPE   ", ncomp)) off = myDHDR.KHCTP;
  else if (!strncmp(keyword, "PTYPE   ", ncomp)) off = myDHDR.KHPTP;
  else if (!strncmp(keyword, "CRVAL   ", ncomp)) off = myDHDR.KDCRV;
  else if (!strncmp(keyword, "CDELT   ", ncomp)) off = myDHDR.KRCIC;
  else if (!strncmp(keyword, "CRPIX   ", ncomp)) off = myDHDR.KRCRP;
  else if (!strncmp(keyword, "CROTA   ", ncomp)) off = myDHDR.KRCRT;
  else if (!strncmp(keyword, "EPOCH   ", ncomp)) off = myDHDR.KREPO;
  else if (!strncmp(keyword, "DATAMAX ", ncomp)) off = myDHDR.KRDMX;
  else if (!strncmp(keyword, "DATAMIN ", ncomp)) off = myDHDR.KRDMN;
  else if (!strncmp(keyword, "PRODUCT ", ncomp)) off = myDHDR.KITYP;
  else if (!strncmp(keyword, "NITER   ", ncomp)) off = myDHDR.KINIT;
  else if (!strncmp(keyword, "BMAJ    ", ncomp)) off = myDHDR.KRBMJ;
  else if (!strncmp(keyword, "BMIN    ", ncomp)) off = myDHDR.KRBMN;
  else if (!strncmp(keyword, "BPA     ", ncomp)) off = myDHDR.KRBPA;
  else if (!strncmp(keyword, "VELREF  ", ncomp)) off = myDHDR.KIALT;
  else if (!strncmp(keyword, "ALTRVAL ", ncomp)) off = myDHDR.KDARV;
  else if (!strncmp(keyword, "ALTRPIX ", ncomp)) off = myDHDR.KRARP;
  else if (!strncmp(keyword, "OBSRA   ", ncomp)) off = myDHDR.KDORA;
  else if (!strncmp(keyword, "OBSDEC  ", ncomp)) off = myDHDR.KDODE;
  else if (!strncmp(keyword, "RESTFREQ", ncomp)) off = myDHDR.KDRST;
  else if (!strncmp(keyword, "XSHIFT  ", ncomp)) off = myDHDR.KRXSH;
  else if (!strncmp(keyword, "YSHIFT  ", ncomp)) off = myDHDR.KRYSH;
  else if (!strncmp(keyword, "NAMCLSTY", ncomp)) off = myDHDR.KHIMN;
  else if (!strncmp(keyword, "IMSEQ   ", ncomp)) off = myDHDR.KIIMS;
  else if (!strncmp(keyword, "USERNO  ", ncomp)) off = myDHDR.KIIMU;
  else if (!strncmp(keyword, "EXTYPE  ", ncomp)) off = myDHDR.KHEXT;
  else if (!strncmp(keyword, "EXTVER  ", ncomp)) off = myDHDR.KIVER;
  else if (!strncmp(keyword, "BLANK   ", ncomp)) off = myDHDR.KRBLK;
  else off = -1; /* unknown */

  return off; 
} /*  end ObitAIPSCatOffset */

/**
 * Create and write dummy AIPS header
 * \param disk   disk number.
 * \param user   AIPS user number.
 * \param Aname  AIPS name.
 * \param Aclass AIPS class.
 * \param Atype  AIPS file type (MA, UV, SC).
 * \param seq    AIPS sequence number.
 * \param cno    AIPS catalog slot number
 * \param err    Obit error stack.
 */
void ObitAIPSCatDummy (olong disk, olong user, 
		       gchar Aname[13], gchar Aclass[7], gchar Atype[3], 
		       olong seq, olong cno, ObitErr *err)
{
  ObitIOCode retCode;
  ObitImageDesc* imDesc  = NULL;
  ObitUVDesc*    uvDesc  = NULL;
  gchar *HeaderFile;
  gsize size;
  olong i;
  ObitFilePos wantPos;
  AIPSint buffer[256];
  ObitFile *myFile=NULL;
  ObitAIPSDirCatEntry *dirEntry = NULL;
  gchar *routine = "ObitAIPSCatDummy";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  if (cno<=0)   return;
  if (disk<=0) return;

  for (i=0; i<256; i++) buffer[i] = 0; /* clear buffer */

  /* Set file name */
  HeaderFile = ObitAIPSFilename (OBIT_AIPS_Header, disk, cno, 
				 user, NULL, 0, err);
  if (err->error) Obit_traceback_msg (err, routine, Aname);

  /* open Header file */
  myFile = newObitFile(Aname);
  size = 256 * sizeof(AIPSint);
  if (ObitFileOpen (myFile, HeaderFile, OBIT_IO_ReadWrite, 
		     OBIT_IO_Binary, size, err) ||
      (err->error)) /* add traceback on error */
    Obit_traceback_msg (err, routine, Aname);
  g_free(HeaderFile); HeaderFile = NULL;  /* cleanup */

  /* Get catalog entry */
  retCode = OBIT_IO_ReadErr;
  /* Get catalog descriptor */
  dirEntry = ObitAIPSDirGetEntry(disk, user, cno, err);
  if (err->error) Obit_traceback_msg (err, routine, Aname);

  /* Create dummy header by type */
  if ((Atype[0]=='U') && (Atype[1]=='V')) { /* UV data */
    uvDesc = newObitUVDesc(NULL); /* Create dummy header */
    uvDesc->lrec   = 8;
    uvDesc->nrparm = 5;
    uvDesc->naxis  = 1;
    uvDesc->inaxes[0] = 3;
    /* do conversion */
    ObitAIPSCatUVSetDesc (uvDesc, (gchar*)buffer, !myFile->exist, 
			  dirEntry, err);
    if (err->error) Obit_traceback_msg (err, routine, Aname);
    uvDesc = ObitUVDescUnref(uvDesc);
  } else { /* other than UV data */
    imDesc = newObitImageDesc(NULL); /* Create dummy header */
    imDesc->naxis     = 0;
    imDesc->inaxes[0] = 2;
    imDesc->inaxes[1] = 2;
   /* do conversion */
    ObitAIPSCatImageSetDesc (imDesc, (gchar*)buffer, !myFile->exist, 
			     dirEntry, err);
    if (err->error) Obit_traceback_msg (err, routine, Aname);
    imDesc = ObitImageDescUnref(imDesc);
  } /* End write header into buffer */

  /* Now write it */
  size = 256 * sizeof(AIPSint);           /* transfer size in bytes */
  wantPos = 0; /* File location */
  retCode = ObitFileWrite (myFile, wantPos, size, (gchar*)buffer, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback on error */
    Obit_traceback_msg (err, routine, Aname);

  /* Write first dummy keyword/value record */
  for (i=0; i<256; i++) buffer[i] = 0; /* clear buffer */
  wantPos = 256 * sizeof(AIPSint); /* File location */
  /* write block */
  retCode = ObitFileWrite (myFile, wantPos, size, (gchar*)buffer, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback on error */
    Obit_traceback_msg (err, routine, Aname);
  
  /* flush/close file */
  retCode = ObitFileClose (myFile, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback on error */
    Obit_traceback_msg (err, routine, Aname);
  
  /* cleanup */
  myFile = ObitFileUnref(myFile);
  if (dirEntry) g_free(dirEntry); /* free catalog directory entry */

} /* end ObitAIPSCatDummy */

/**
 * Chane AIPS name in header
 * \param disk     Disk number
 * \param user     AIPS user number.
 * \param cno      Slot number to be renamed
 * \param newName  New AIPS name (12 characters)
 * \param newClass New AIPS Class (6 characters)
 * \param newSeq   New AIPS sequence
 * \param err      Obit error stack.
 */
void ObitAIPSCatRename(olong disk, olong user,  olong cno, gchar *newName, 
		      gchar *newClass, olong newSeq, ObitErr *err)
{
  ObitIOCode retCode;
  gchar *HeaderFile;
  gsize size;
  ObitFilePos wantPos;
  AIPSint buffer[256];
  ObitFile *myFile=NULL;
  gchar *cp;
  gchar *routine = "ObitAIPSCatRename";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  if (cno<=0)   return;
  if (disk<=0) return;

  /* Set file name */
  HeaderFile = ObitAIPSFilename (OBIT_AIPS_Header, disk, cno, 
				 user, NULL, 0, err);
  if (err->error) Obit_traceback_msg (err, routine, newName);

  /* open Header file */
  myFile = newObitFile(newName);
  size = 256 * sizeof(AIPSint);
  if (ObitFileOpen (myFile, HeaderFile, OBIT_IO_ReadWrite, 
		     OBIT_IO_Binary, size, err) ||
      (err->error)) /* add traceback on error */
    Obit_traceback_msg (err, routine, newName);
  g_free(HeaderFile); HeaderFile = NULL;  /* cleanup */

  /* Read header record */
  size = 256 * sizeof(AIPSint);           /* transfer size in bytes */
  wantPos = 0; /* File location */
  retCode = ObitFileRead (myFile, wantPos, size, (gchar*)buffer, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback on error */
    Obit_traceback_msg (err, routine, newName);

  /* reset AIPS naming info  */
  cp = (gchar*)&buffer[myDHDR.KHIMN];
  CopyDeeNull (cp, newName, 12);
  CopyDeeNull (cp+12, newClass, 6);
  /* Image sequence no. */
  buffer[myDHDR.KIIMS] = newSeq;

  /* Now write it */
  size = 256 * sizeof(AIPSint);           /* transfer size in bytes */
  wantPos = 0; /* File location */
  retCode = ObitFileWrite (myFile, wantPos, size, (gchar*)buffer, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback on error */
    Obit_traceback_msg (err, routine, newName);

  /* flush/close file */
  retCode = ObitFileClose (myFile, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback on error */
    Obit_traceback_msg (err, routine, newName);
  
  /* cleanup */
  myFile = ObitFileUnref(myFile);

} /* end ObitAIPSCatRename */

/*---------------Private functions---------------------------*/

/**
 * Returns date and time in an AIPSint array.
 * \param dateTime date and time as yyyy, mm, dd, hh, mm, ss.
 */
static void
ObitAIPSCatUpdateAccess(AIPSint *dateTime)
{
  struct tm *lp;
  time_t clock;

  /* Get time since 00:00:00 GMT, Jan. 1, 1970 in seconds. */
  time (&clock);

  /* Convert to  broken-down time. */
  lp = localtime (&clock);
  lp->tm_mon++; /* For some bizzare reason, month is 0-rel */

  /* to output */
  dateTime[0] = lp->tm_year;
  if (dateTime[0]<1000)  dateTime[0] += 1900; /* full year */
  dateTime[1] = MAX (1, lp->tm_mon);
  dateTime[2] = lp->tm_mday;
  dateTime[3] = lp->tm_hour;
  dateTime[4] = lp->tm_min;
  dateTime[5] = lp->tm_sec;
} /* end ObitAIPSCatUpdateAccess */

/**
 * Copies string replacing any non printing characters with blank.
 * \param out    destination Character array 
 * \param in     source character array
 * \param nchar  Number of characters to test
 */
static void CopyDeeNull (gchar *out, gchar *in, olong nchar)
{
  olong i;

  for (i=0; i<nchar; i++) {
    if (!g_ascii_isprint(in[i])) out[i] = ' ';
    else out[i] = in[i];
  }
} /* end  DeeNull */

