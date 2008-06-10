/* $Id: ObitTableDef.h,v 1.4 2006/01/24 13:25:09 bcotton Exp $  */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2008                                          */
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
/*  Define the basic components of the ObitTable structure            */
/*  This is intended to be included in a class structure definition   */
/**
 * \file ObitTableDef.h
 * ObitTable structure members for derived classes.
 */
#include "ObitDef.h"  /* Parent class definitions */
/** Threading info member object  */
ObitThread *thread;
/** Linked list of arrays of data.  */
ObitInfoList *info;
/** I/O status */
ObitIOStatus myStatus;
/** I/O interface to access disk resident data */
ObitIO* myIO;
/** Obit Table data Descriptor */
ObitTableDesc* myDesc;
/** Obit Table data Selector */
ObitTableSel* mySel;
/** Table buffer */
ofloat *buffer;
/** Table buffer size in floats */
olong bufferSize;
/** Table type (name of type, e.g. "AIPS AN") */
gchar *tabType;
/** Table version number(1-rel) */
olong tabVer;
/** This table attached to host (actually) ObitData 
this is a secret reference - don't use Ref or Unref functions */
Obit *myHost;
