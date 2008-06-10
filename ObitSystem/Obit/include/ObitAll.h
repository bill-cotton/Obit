/* $Id: ObitAll.h,v 1.2 2005/10/06 20:22:54 bcotton Exp $         */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003                                               */
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
#ifndef OBITALL_H 
#define OBITALL_H 

/**
 * \file ObitAll.h
 * Obit class definition.
 *
 * This file should be used in the main program to define the
 * Obit class structures.
 * Also, ObitInitAll should be called to initialize
 * in order to properly define the global class structures.
 */
/*-------- Obit: Merx mollis mortibus nuper ------------------*/

#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <glib.h>
#include "ObitSystem.h"
#include "ObitErr.h"
#include "ObitThread.h"
#include "ObitInfoList.h"
#include "ObitImage.h"
#include "ObitIO.h"
#include "ObitIOImageAIPS.h"
#include "ObitIOImageFITS.h"
#include "ObitAIPS.h"
#include "ObitAIPSDir.h"
#include "ObitUV.h"
#include "ObitUVWeight.h"
#include "ObitTable.h"
#include "ObitTableFQ.h"
#include "ObitTableFQUtil.h"
#include "ObitTableAN.h"
#include "ObitTableANUtil.h"
#include "ObitTableSU.h"
#include "ObitTableSUUtil.h"
#include "ObitTableSN.h"
#include "ObitTableBL.h"
#include "ObitTableBP.h"
#include "ObitTableNX.h"
#include "ObitTableNI.h"
#include "ObitImageUtil.h"
#include "ObitIoN2SolNTable.h"
#endif /* OBITALL_H */ 
