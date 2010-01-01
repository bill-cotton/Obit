/* $Id:      */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2005-2008                                          */
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
#ifndef OBITRPCUTIL_H 
#define OBITRPCUTIL_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitXML.h"
#include "ObitRPC.h"
#include "ObitFile.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitRPCUtil.h
 * ObitRPC Remote Procedure Call Utility package for ObitRPC class
 *
 */

/*---------------Public functions---------------------------*/
/** Public: Copy a file over an XMLRPC connection */
void ObitRPCUtilFileSend (ObitFile *file, ObitRPC *client, 
			  gchar *serverURL, ObitErr *err);

#endif /* OBITFRPCUTIL_H */ 
