/* $Id: xml.h,v 1.3 2005/11/29 18:02:32 bcotton Exp $  */
/* xml routines for ObitView */
/*-----------------------------------------------------------------------
*  Copyright (C) 2005
*  Associated Universities, Inc. Washington DC, USA.
*  This program is free software; you can redistribute it and/or
*  modify it under the terms of the GNU General Public License as
*  published by the Free Software Foundation; either version 2 of
*  the License, or (at your option) any later version.
*
*  This program is distributed in the hope that it will be useful,
*  but WITHOUT ANY WARRANTY; without even the implied warranty of
*  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*  GNU General Public License for more details.
*  Correspondence concerning ObitView should be addressed as follows:
*         Internet email: bcotton@nrao.edu.
*         Postal address: William Cotton
*                         National Radio Astronomy Observatory
*                         520 Edgemont Road
*                         Charlottesville, VA 22903-2475 USA
*-----------------------------------------------------------------------*/
#include <pthread.h>
/*#include <xmlrpc-c/base.h>*/
#include <xmlrpc.h>
#include <glib.h>
#include "ObitDConCleanWindow.h"

#ifndef XML_H
#define XML_H
/*------------------ Structures -----------------------------*/


/*----------------- Globals ---------------------------*/

/*---------------Public functions---------------------------*/
xmlrpc_value* Window2XML (xmlrpc_env *env, ObitDConCleanWindow *window, 
			  olong field, ObitErr *err);
ObitDConCleanWindow* XML2Window (xmlrpc_env *env, xmlrpc_value* xmlwindow,
				 ObitErr *err);
xmlrpc_value* 
FileInfo2XML (xmlrpc_env *env, 
	      ObitIOType Type, gchar *Name,
	      gchar *AClass, gchar *ADir, olong ASeq, olong AUser,
	      ObitErr *err);

gint
XML2FileInfo (xmlrpc_env *env, xmlrpc_value* xmlinfo,
	      ObitIOType *Type, gchar **Name,
	      gchar **AClass, gchar **ADir, olong *ASeq, olong *AUser,
	      ObitErr *err);
#endif /* end XML_H */


