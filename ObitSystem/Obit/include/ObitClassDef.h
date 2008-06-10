/* $Id$  */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2002-2008                                          */
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
/*  Define the basic components of the Obit InfoClass structure       */
/* This is intended to be included in a classInfo structure definition*/
/** Have I been initialized? */
gboolean initialized;
/** Are disk resident "scratch" objects of this class possible? */
gboolean hasScratch;
/** Name of class ("Obit") */
gchar* ClassName;
/** Pointer to parent class ClassInfo, Null if none. */
gconstpointer ParentClass;
/** Function pointer to Class initializer */
ObitClassInitFP ObitClassInit;
/** Function pointer to newObit. */
newObitFP newObit;
/** Function pointer to GetClass. */
ObitGetClassFP ObitGetClass;
/** Function pointer to ClassInfoDefFn. */
ObitClassInfoDefFnFP ObitClassInfoDefFn;
/** Function pointer to shallow copy constructor. */
ObitCopyFP ObitCopy;
/** Function pointer to deep copy constructor. */
ObitCloneFP ObitClone;
/** Function pointer to Object Ref. */
ObitRefFP ObitRef;
/** Function pointer to Object Unref. */
ObitUnrefFP ObitUnref;
/** Function pointer to test if a class member. */
ObitIsAFP ObitIsA;
/* private functions - added in first call to newObit */
/** Function pointer to deallocation function. */
ObitClearFP ObitClear;
/** Function pointer to object initializer. */
ObitInitFP ObitInit;
