/* $Id$ */
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
#ifndef OBITAIPSCAT_H 
#define OBITAIPSCAT_H 
#include <glib.h>
#include "ObitErr.h"
#include "ObitImageDesc.h"
#include "ObitUVDesc.h"
#include "ObitTableDesc.h"
#include "ObitAIPS.h"
#include "ObitAIPSDir.h"
#include "ObitTableList.h"

/*-------- ObitIO: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitAIPSCat.h
 * ObitAIPSCat  class definition.
 *
 * This is a Utility class to handle the interface with the AIPS catalog
 * system.
 * Also handles some of the details of AIPS table headers.
 * This class is non-derivable and has no public instances.
 */
/*---------------Public functions---------------------------*/
/** Public: Convert AIPS header file data to an ObitImageDesc. */
void ObitAIPSCatImageGetDesc (ObitImageDesc *desc, gchar *buffer, 
			      ObitErr *err);

/** Public: Convert an ObitImageDesc to AIPS header file data. */
void ObitAIPSCatImageSetDesc (ObitImageDesc *desc, gchar *buffer, 
			      gboolean init, ObitAIPSDirCatEntry* dirEntry,
			      ObitErr *err);

/** Public: Convert AIPS header file data to an ObitUVDesc. */
void ObitAIPSCatUVGetDesc (ObitUVDesc *desc, gchar *buffer, 
			      ObitErr *err);

/** Public: Convert an ObitUVDesc to AIPS header file data. */
void ObitAIPSCatUVSetDesc (ObitUVDesc *desc, gchar *buffer, 
			      gboolean init, ObitAIPSDirCatEntry* dirEntry,
			      ObitErr *err);

/** Public: Copy AIPS  header to an ObitTableList. */
void ObitAIPSCatGetTable (ObitTableList *tableList, gchar *buffer, 
			  olong user, olong disk, olong cno, ObitErr *err);

/** Public: Copy an ObitTableList to AIPS table header. */
void ObitAIPSCatSetTable (ObitTableList *tableList, gchar *buffer, 
			  ObitErr *err);

/** Public: Convert AIPS Table header to an ObitTableDesc. */
void ObitAIPSCatTableGetDesc (ObitTableDesc *desc, 
			      gchar tabType[3], olong tabVer,
			      AIPSint controlBlock[256],
			      AIPSint record[256], ObitErr *err);

/** Public: Convert an ObitTableDesc to AIPS table header. */
void ObitAIPSCatTableSetDesc (ObitTableDesc *desc, gboolean init, 
			      gchar tabType[3], olong tabVer,
			      AIPSint controlBlock[256],
			      AIPSint record[256], ObitErr *err);

/** Public: Initialize AIPS Header info object */
void ObitAIPSCatInitDHDR(void);

/** Public: Return catalog header offset for by keyword */
olong ObitAIPSCatOffset (gchar *keyword);

/** Public: Write dummy AIPS header */
void ObitAIPSCatDummy (gint disk, olong user, 
		       gchar Aname[13], gchar Aclass[7], gchar Atype[3], 
		       olong seq, olong cno, ObitErr *err);

/** Public: Rename cataloged name, class, seq */
void ObitAIPSCatRename(gint disk, olong user,  olong cno, gchar *newName, 
		      gchar *newClass, olong newSeq, ObitErr *err);

#endif /* OBITAIPSCAT_H */ 

