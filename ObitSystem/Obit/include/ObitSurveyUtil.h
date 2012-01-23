/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2012                                               */
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
/*; Correspondence about this software should be addressed as follows:*/
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
#ifndef OBITSURVEYUTIL_H 
#define OBITSURVEYUTIL_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitTableVL.h"
#include "ObitTableVZ.h"
#include "ObitFitRegionList.h"
#include "ObitPrinter.h"
/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitSurveyUtil.h
 * ObitSurvey class utility routine definition.
 */

/*---------------Public functions---------------------------*/
/** Public: Print all contents of a VL table */
void ObitSurveyUtilVLPrint (ObitTableVL *in, ObitImage *image, FILE  *prtFile, 
			    ObitErr *err);
/** Public: Printed selected NVSS survey catalog entries */
gboolean ObitSurveyNVSSPrint (ObitPrinter *printer, ObitData *data, olong VLVer, 
			      gboolean first, gboolean last, ObitErr* err);

/** Public: Printed selected VLSS survey catalog entries */
gboolean ObitSurveyVLSSPrint (ObitPrinter *printer, ObitData *data, olong VLVer, 
			      gboolean first, gboolean last, ObitErr* err);
#endif /* OBITSURVEYUTIL_H */ 
