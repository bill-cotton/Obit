/* $Id: ObitTableListClassDef.h,v 1.1.1.1 2004/07/19 16:42:41 bcotton Exp $                            */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2002-2003                                          */
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
/* Define the basic components of the ObitTableList InfoClass structure*/
/* This is intended to be included in a classInfo structure definition*/
/**
 * \file ObitTableListClassDef.h
 * ObitTableListClassDef ClassInfo structure members for derived classes.
 */
#include "ObitClassDef.h"  /* Parent class ClassInfo definition file */
/** Function pointer to Add item to list */
ObitTableListPutFP ObitTableListPut;
/** Function pointer to obtain item by name */
ObitTableListGetFP ObitTableListGet;
/** Function pointer to obtain item by number */
ObitTableListGetNumberFP ObitTableListGetNumber;
/** Function pointer to remove item from list */
ObitTableListRemoveFP ObitTableListRemove;
