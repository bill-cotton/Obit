/* $Id$     */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2002-2011                                          */
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
/*;  Correspondence this software should be addressed as follows:     */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
#ifndef OBITTYPES_H 
#define OBITTYPES_H 
#include <glib.h>

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitTypes.h
 * Obit type definition file.
 *
 * This defines various enumerations and types to be used in Obit.
 */

/*-------------------------- Macroes  -----------------------------------*/
/** Degrees to radians factor */
#ifndef DG2RAD  
#define DG2RAD G_PI / 180.0
#endif

/**  Radians to degrees factor */
#ifndef RAD2DG  
#define RAD2DG 180.0 / G_PI
#endif

/**  Radians to arcsecond factor */
#ifndef RAD2AS  
#define RAD2AS 3600.0 * RAD2DG
#endif

/**  Arcsecond  to Radians factor */
#ifndef AS2RAD  
#define AS2RAD DG2RAD / 3600.0
#endif
/*-------------- type definitions----------------------------------*/
/** Typedef for equivalent to FORTRAN INTEGER - should be size of float */
/* if sizeof (gint32) == sizeof (gfloat) */
typedef gint32 oint;

/* if sizeof(glong) == sizeof(gfloat)
   typedef glong oint;*/

/** Typedef for Obit 32 integers */
typedef gint32 olong;

/** Typedef for Obit 64 bit integers */
typedef gint64 ollong;

/** Typedef for Obit floats */
typedef gfloat ofloat;


/** Typedef for Obit doubles */
typedef gdouble odouble;

/*-------------------- unions -------------------------------------*/
/** Equivalence for obtaining a single element of an uncertain type
   from InfoList */
  union ObitInfoListEquiv { 
    olong   otg;
    oint    itg;
    ofloat  flt;
    odouble dbl;
  };

/*-------------- enumerations -------------------------------------*/
/**
 * \enum obitIOStatus
 * enum for object status.
 * This specifies the status of the connection a disk resident data.
 */
enum obitIOStatus {
  /** Unconnected to disk resident data, no information defined */
  OBIT_Inactive, 
  /** Descriptor created but not connected to a disk file(s) */
  OBIT_Defined,  
  /** Fully defined and connected to disk resident structures */
  OBIT_Active,
  /** Fully defined and I/O active, buffers modified. */
  OBIT_Modified,
  /** An error has been encountered. */
  OBIT_ErrorExist
}; /* end enum obitIOStatus */
/** typedef for enum for ObitIO object status. */
typedef enum obitIOStatus ObitIOStatus;

/**
 * \enum obitIOType
 * enum for I/O file type.
 * This specifies the type of underlying data file.
 */
enum obitIOType {
  /** FITS */
  OBIT_IO_FITS=0, 
  /** AIPS */
  OBIT_IO_AIPS,
  /** Memory only */
  OBIT_IO_MEM,
  /** Binary data */
  OBIT_IO_Binary,
  /** Text file */
  OBIT_IO_Text
}; /* end enum obitIOType */

/** typedef for enum for ObitIO object status. */
typedef enum obitIOType ObitIOType;

/**
 * \enum obitIOAccess
 * enum for I/O access type.
 * This specifies Read or read/write access
 */
enum obitIOAccess {
  /** None (closed) */
  OBIT_IO_None, 
  /** Read Only */
  OBIT_IO_ReadOnly, 
  /** Read with selection, calibration, editing (UV data only) */
  OBIT_IO_ReadCal, 
  /** Write only  */
  OBIT_IO_WriteOnly,
  /** Read and/or write */
  OBIT_IO_ReadWrite
}; /* end enum obitIOAccess */

/** typedef for enum for ObitIO object status. */
typedef enum obitIOAccess ObitIOAccess;

/**
 * \enum obitIOSize
 * enum for I/O data transfer size
 * This specifies Read or read/write access
 */
enum obitIOSize {
  /** Transfer by image row */
  OBIT_IO_byRow,
  /** Transfer by image plane */
  OBIT_IO_byPlane
}; /* end enum obitIOSize */

/** typedef for enum for ObitIO Transfer block size. */
typedef enum obitIOSize ObitIOSize;

/**
 * \enum obitIOCode
 * enum for I/O return codes.
 */
enum obitIOCode {
  /** OK */
  OBIT_IO_OK=0, 
  /** Hit end of file - all higher values represent real errors */
  OBIT_IO_EOF, 
  /** Specification error */
  OBIT_IO_SpecErr,
  /** Open error  */
  OBIT_IO_OpenErr,
  /** Close error  */
  OBIT_IO_CloseErr,
  /** Read Error */
  OBIT_IO_ReadErr,
  /** Write Error */
  OBIT_IO_WriteErr,
  /** Delete error */
  OBIT_IO_DeleteErr,
  /** Item not found error */
  OBIT_IO_NotFoundErr
}; /* end enum obitIOCode */

/** typedef for enum for ObitIO object status. */
typedef enum obitIOCode ObitIOCode;

/**
 * \enum obitInfoType
 * enum for ObitInfoElem element data types.
 * For integers use OBIT_long
 */
enum obitInfoType {
  /** 8-bit signed byte */
  OBIT_byte = 0,  
  /** 16 bit signed integer */
  OBIT_short, 
  /** signed integer */
  OBIT_int,   
  /** Equivalent of FORTRAN INTEGER */
  OBIT_oint,   
  /** 32 bit signed integer */
  OBIT_long,  
  /** 8-bit unsigned byte */
  OBIT_ubyte, 
  /** 16 bit unsigned integer  */
  OBIT_ushort,
  /**  unsigned integer */
  OBIT_uint,  
  /** 32 bit unsigned signed integer */
  OBIT_ulong, 
  /** single precision float */
  OBIT_float, 
  /** double precision float */
  OBIT_double,
  /** single precision complex */
  OBIT_complex,
  /** double precision complex */
  OBIT_dcomplex,
  /** character */
  OBIT_string,
  /** boolean */
  OBIT_bool,
  /** Bit array packed into gints */
  OBIT_bits
}; /* end enum obitInfoType */

/** typedef for enum for ObitInfoElem element data types. */
typedef enum obitInfoType ObitInfoType;


#endif /* OBITTYPES_H */ 
