/* $Id$        */
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
/*;Correspondence about this software should be addressed as follows: */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/

#include <unistd.h> 
#include <stdio.h>    /* i/o include */
#include <stdlib.h>   /* i/o include */
#include <signal.h>
#include "ObitXML.h"
#include "ObitMem.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitXML.c
 * ObitXML class function definitions.
 * This class is derived from the Obit base class.
 * The XML implementation is based on xmlrpc
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitXML";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo structure ObitXMLClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitXMLClassInfo myClassInfo = {FALSE};

/*--------------- File Global Variables  ----------------*/


/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitXMLInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitXMLClear (gpointer in);

/** Convert InfoList to XML */
static void encodeInfoList (ObitInfoList *desc, 
			    ObitXMLEnv *envP, ObitXMLValue *parmP, 
			    ObitErr *err);

/** Convert ObitXMLValue to InfoList */
static void decodeInfoList (ObitXMLEnv *envP, ObitXMLValue *parmP, 
			    ObitInfoList **out,  ObitErr *err);

/** Private: Set Class function pointers. */
static void ObitXMLClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitXML* newObitXML (gchar* name)
{
  ObitXML* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitXMLClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitXML));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitXMLInit((gpointer)out);

  return out;
} /* end newObitXML */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitXMLGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitXMLClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitXMLGetClass */

/**
 * Create an XML from PRC call arguments received to be sent from client end
 * type = OBIT_XML_RPCCallArg
 * \param func         Name of function to be called
 * \param argList      InfoList with call arguments
 * \param err     Obit Error message
 * \return new ObitXML object suitable to be passed to ObitRPCCall
 */
ObitXML* 
ObitXMLSetCallArg (gchar* func, ObitInfoList *argList, ObitErr *err)
{
  ObitXML  *out=NULL;
  xmlrpc_value *v;
  gchar *routine = "ObitXMLSetCallArg";

  /* error checks */
  if (err->error) return out;

  /* initialize structure */
  out = newObitXML(routine);
  out->type = OBIT_XML_RPCCallArg;
  out->name = g_strdup(routine);
  out->func = g_strdup(func);
  out->parmP = xmlrpc_struct_new(&out->envP);

   /* ObitXML type */
  v = xmlrpc_build_value(&out->envP, "i", (xmlrpc_int32)OBIT_XML_InfoList);
  xmlrpc_struct_set_value(&out->envP, out->parmP, "XMLType", v);
  xmlrpc_DECREF(v);
 
  /* number of entries */
  v = xmlrpc_build_value(&out->envP, "i", (xmlrpc_int32)argList->number);
  xmlrpc_struct_set_value(&out->envP, out->parmP, "number", v);
  xmlrpc_DECREF(v);
  XMLRPC_FAIL_IF_FAULT(&out->envP);

  /* Extract from argList and encode */
  encodeInfoList (argList, (ObitXMLEnv*)&out->envP, (ObitXMLValue*)out->parmP, err);
  if (err->error) Obit_traceback_val (err, routine, out->name, out);  

  /* Make sure everything is cool */
 cleanup:
  if (out->envP.fault_occurred) {
    Obit_log_error(err, OBIT_Error, "%s:XML-RPC Fault: %s (%d)",
		   routine, out->envP.fault_string, out->envP.fault_code);
  }
  
  return out;
} /* end ObitXMLSetCallArg */

/**
 * Convert call arguments received on server end to an InfoList
 * type = OBIT_XML_RPCCallArg
 * \param envP         Environment, structure copied
 * \param paramArrayP  Parameter array, structure NOT copied
 * \param err     Obit Error message
 * \return new ObitInfoList object
 */
ObitInfoList* 
ObitXMLGetCallArg (ObitXMLEnv * const envP, ObitXMLValue * const paramArrayP,
		   ObitErr *err)
{
  ObitInfoList  *out=newObitInfoList();
  xmlrpc_value *v;
  gchar *routine = "ObitXMLGetCallArg";

  /* error checks */
  if (err->error) return out;

  /* Data comes packed into an array */
  xmlrpc_decompose_value(envP, paramArrayP, "(S)", &v);

  /* Convert v to the InfoList */
  decodeInfoList (envP, (ObitXMLValue*)v, &out, err);
  if (err->error) Obit_traceback_val (err, routine, routine, out);  
  xmlrpc_DECREF(v);
  XMLRPC_FAIL_IF_FAULT(envP);
 
  /* Make sure everything is cool */
 cleanup:
  if (envP->fault_occurred) {
    Obit_log_error(err, OBIT_Error, "%s:XML-RPC Fault: %s (%d)",
		   routine, envP->fault_string, envP->fault_code);
  }
  
  return out;
} /* end ObitXMLGetCallArg */

/**
 * Give the XML Environment
 * \param in    XML object
 * \return pointer to XML Environment
 */
ObitXMLEnv ObitXMLGetEnv (ObitXML* in)
{
  return (ObitXMLEnv)in->envP;
} /* end ObitXMLGetEnv */

/**
 * Give the XML Value
 * \param in    XML object
 * \return pointer to XML Value
 */
ObitXMLValue* ObitXMLGetValue (ObitXML* in)
{
  xmlrpc_INCREF(in->parmP);
  return (ObitXMLValue*)in->parmP;
} /* end ObitXMLGetValue */


/**
 * Convert ping data to XML
 * Intended for RPC function "ping"
 * A ping call is just passed an arbitrary gint, use 42
 * type = OBIT_XML_Ping
 * \param err     Obit Error message
 * \return new ObitXML object
 */
ObitXML* 
ObitXMLPing2XML (ObitErr *err)
{
  ObitXML  *out=NULL;
  gchar *routine = "ObitXMLPing2XML";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return out;

  /* initial structure */
  out = newObitXML(routine);
  out->type = OBIT_XML_Ping;
  out->func = g_strdup("ping");
  
  out->parmP = xmlrpc_build_value(&out->envP, "i", (xmlrpc_int32)42);
  
  /* Make sure everything is cool */
  if (out->envP.fault_occurred) {
    Obit_log_error(err, OBIT_Error, "%s:XML-RPC Fault: %s (%d)",
		   routine, out->envP.fault_string, out->envP.fault_code);
  }
  
  return out;
} /* end ObitXMLPing2XML */

/**
 * Convert an ObitXML to ping
 * \param xml     ObitXML Object from which to extract information
 * \param err     Obit Error message
 * \return the random gint
 */
olong
ObitXMLXML2Ping (ObitXML* xml, ObitErr *err)
{
  xmlrpc_int32 out = -1;
  gchar *routine = "ObitXMLXML2Ping";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return (olong)out;
  Obit_retval_if_fail((xml->type == OBIT_XML_Ping), err, out,
		      "%s: xml wrong type %d != %d", 
		      routine, xml->type, OBIT_XML_Ping);

  xmlrpc_decompose_value(&xml->envP, xml->parmP, "i", &out);

  /* Make sure everything is cool */
  if (xml->envP.fault_occurred) {
    Obit_log_error(err, OBIT_Error, "XML-RPC Fault: %s (%d)",
		   xml->envP.fault_string, xml->envP.fault_code);
    return (olong)out;
  }

  return (olong)out;
} /* end ObitXMLXML2Ping */

/**
 * Convert ObitInfoList to XML
 * Note: Can only translate structures with only data types
 * directly translatable into XML.
 * type = OBIT_XML_InfoList
 * \param list    List to convert
 * \param err     Obit Error message
 * \return  new ObitXML object (Unref when done)
 */
ObitXML* 
ObitXMLInfoList2XML (ObitInfoList *list, ObitErr *err)
{
  ObitXML  *out=NULL;
  xmlrpc_value *v;
  gchar *routine = "ObitXMLInfoList2XML";
  
  /* error checks */
  if (err->error) return out;
  g_assert (ObitInfoListIsA(list));

  /* initial structure */
  out = newObitXML(routine);
  out->parmP = xmlrpc_struct_new(&out->envP);
  out->type = OBIT_XML_InfoList;
  out->func = g_strdup("undefined");

  /* ObitXML type */
  v = xmlrpc_build_value(&out->envP, "i", (xmlrpc_int32)OBIT_XML_InfoList);
  xmlrpc_struct_set_value(&out->envP, out->parmP, "XMLType", v);
  xmlrpc_DECREF(v);
 
  /* number of entries */
  v = xmlrpc_build_value(&out->envP, "i", (xmlrpc_int32)list->number);
  xmlrpc_struct_set_value(&out->envP, out->parmP, "number", v);
  xmlrpc_DECREF(v);
  XMLRPC_FAIL_IF_FAULT(&out->envP);

  /* Extract from Obit Object and encode */
  encodeInfoList (list, (ObitXMLEnv*)&out->envP, (ObitXMLValue*)out->parmP, err);
  if (err->error) Obit_traceback_val (err, routine, out->name, out);  

    /* Make sure everything is cool */
 cleanup:
  if (out->envP.fault_occurred) {
    Obit_log_error(err, OBIT_Error, "%s:XML-RPC Fault: %s (%d)",
		   routine, out->envP.fault_string, out->envP.fault_code);
  }
  
  return out;
} /* end ObitXMLInfoList2XML */

/**
 * Convert XML struct to ObitInfoList
 * Works on relatively arbitrary XML but requires 
 * all data in a given entry to be the same type.
 * \param xml     Object to convert
 * \param err     Obit Error message
 * \return  new ObitInfoList (Unref when done)
 */
ObitInfoList*
ObitXMLXML2InfoList (ObitXML *xml, ObitErr *err)
{
  ObitInfoList  *out=NULL;
  ObitInfoType infoType=0;
  xmlrpc_int32 xmlInt;
  xmlrpc_double xmlDouble;
  xmlrpc_bool xmlBool;
  xmlrpc_type xmlType, xmlTypeA=0;
  const char* xmlChar, **xmlCArray=NULL;
  olong size, maxchar = 0;
  gint32 dim[MAXINFOELEMDIM];
  gpointer data=NULL;
  xmlrpc_value *v=NULL, *k=NULL, *e=NULL;
  olong i, j, jj, num, number;
  gchar *key;
  gchar *routine = "ObitXMLXML2InfoList";
  
  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return out;
  g_assert (ObitXMLIsA(xml));

  /* Be sure xml->parmP a struct */
  xmlType = xmlrpc_value_type (xml->parmP);
  Obit_retval_if_fail(((xml->type == OBIT_XML_InfoList) || 
		       (xml->type == OBIT_XML_RPCCallArg)), err, out,
		      "%s: xml wrong type %d != %d", 
		      routine, xml->type, OBIT_XML_InfoList);
  
  /* initial structure */
  out = newObitInfoList();

  /* How many entries? */
  number = xmlrpc_struct_size (&xml->envP, xml->parmP);
  XMLRPC_FAIL_IF_FAULT(&xml->envP);

  /* Extract from XML Object and copy to InfoList */
  for (i=0; i<number; i++) {
    xmlrpc_struct_read_member (&xml->envP, xml->parmP, (unsigned int)i,
			       &k, &v);
      XMLRPC_FAIL_IF_FAULT(&xml->envP);

    /* Get key (- name) */
    xmlrpc_read_string (&xml->envP, k, &xmlChar);
    key = (gchar*)xmlChar;

    /* Handle by type */
    dim[0] = dim[1] = dim[2] = dim[3] = dim[4] = 1;
    xmlType = xmlrpc_value_type (v);
    switch (xmlType) { 
    case XMLRPC_TYPE_INT :
      /* scalar long */
      infoType = OBIT_long;
      data = g_malloc0(sizeof(olong));
      xmlrpc_read_int (&xml->envP, v, &xmlInt);
      XMLRPC_FAIL_IF_FAULT(&xml->envP);
      ((olong*)data)[0] = (olong)xmlInt;
      break;

    case XMLRPC_TYPE_BOOL:
      /* scalar boolean */
      infoType = OBIT_bool;
      data = g_malloc0(sizeof(gboolean));
      xmlrpc_read_bool (&xml->envP, v, &xmlBool);
      XMLRPC_FAIL_IF_FAULT(&xml->envP);
      ((gboolean*)data)[0] = (gboolean)xmlBool;
      break;

    case XMLRPC_TYPE_DOUBLE:
      /* scalar double */
      infoType = OBIT_double;
      data = g_malloc0(sizeof(odouble));
      xmlrpc_read_double (&xml->envP, v, &xmlDouble);
      XMLRPC_FAIL_IF_FAULT(&xml->envP);
      ((odouble*)data)[0] = (odouble)xmlDouble;
      break;

    case XMLRPC_TYPE_STRING:
      /* single string */
      infoType = OBIT_string;
      xmlrpc_read_string (&xml->envP, v, &xmlChar);
      XMLRPC_FAIL_IF_FAULT(&xml->envP);
      data = g_strdup (xmlChar);
      dim[0] = strlen (xmlChar);
      break;

   case XMLRPC_TYPE_ARRAY:
     /* Loop through array - all MUST be the same data type */
     num = xmlrpc_struct_size (&xml->envP, v);
     dim[0] = num;

     for (j=0; j<num; j++) {
       e = xmlrpc_array_get_item (&xml->envP, v, (int)j);

       /* Check data type */
       if (j==0) xmlTypeA = xmlrpc_value_type (e);
       else  Obit_retval_if_fail((xmlTypeA==xmlrpc_value_type (e)), err, out,
		      "%s: All data in array NOT the same datatype", routine);

       switch (xmlTypeA) {
       case XMLRPC_TYPE_INT:
	 infoType = OBIT_long;
	 if (j==0) data = g_malloc0(num*sizeof(olong));
	 xmlrpc_read_int (&xml->envP, e, &xmlInt);
	 XMLRPC_FAIL_IF_FAULT(&xml->envP);
	 ((olong*)data)[j] = (olong)xmlInt;
	 break;
       case XMLRPC_TYPE_BOOL:
	 infoType = OBIT_bool;
	 if (j==0) data = g_malloc0(num*sizeof(gboolean));
	 xmlrpc_read_bool (&xml->envP, e, &xmlBool);
	 XMLRPC_FAIL_IF_FAULT(&xml->envP);
	 ((gboolean*)data)[j] = (gboolean)xmlBool;
	 break;
       case XMLRPC_TYPE_DOUBLE:
	 infoType = OBIT_double;
	 if (j==0) data = g_malloc0(num*sizeof(odouble));
	 xmlrpc_read_double (&xml->envP, e, &xmlDouble);
	 XMLRPC_FAIL_IF_FAULT(&xml->envP);
	 ((odouble*)data)[j] = (gboolean)xmlDouble;
	 break;
       case XMLRPC_TYPE_STRING:
	 infoType = OBIT_string;
	 /* Strings are a bitch */
	 if (j==0) xmlCArray = g_malloc0(num*sizeof(char*));
	 xmlrpc_read_string (&xml->envP, v, &xmlChar);
	 XMLRPC_FAIL_IF_FAULT(&xml->envP);
	 maxchar = MAX (maxchar, strlen (xmlChar));
	 xmlCArray[j] = g_strdup(xmlChar);
	 break;
       default:
	 g_assert_not_reached(); /* unknown, barf */
       }; /* end switch element type */
     } /* end loop over array */

     if (infoType == OBIT_string) {
       /* put all the strings together without NULLs */
       dim[0] = maxchar; dim[1] = num;
       size = dim[0] * dim[1];
       data = g_malloc0(size*sizeof(char)+5);
       for (j=0; j<size; j++) ((gchar*)data)[j] = ' ';
       for (j=0; j<num; j++) {
	 for (jj=0; jj<maxchar; jj++) {
	   ((gchar*)data)[j*maxchar+jj] = xmlCArray[j][jj];
	   g_free((gchar*)xmlCArray[j]);
	 }
       }
     } /* end of glue strings together */

    default:
      g_assert_not_reached(); /* unknown, barf */
    }; /* end switch basic type */
    XMLRPC_FAIL_IF_FAULT(&xml->envP);

    /* Add to ObitInfoList */
    xmlrpc_DECREF(v);
    if (infoType!=OBIT_string)
      ObitInfoListPut(out, key, infoType, dim, data, err);
    else /* Strings - blah */
       ObitInfoListPut(out, key, infoType, dim, data, err);
     
    g_free(data);
    if (err->error)Obit_traceback_val (err, routine, "writing ObitInfoList", out);
  } /* end loop over elements */

    /* Make sure everything is cool */
 cleanup:
  if (xml->envP.fault_occurred) {
    Obit_log_error(err, OBIT_Error, "%s:XML-RPC Fault: %s (%d)",
		   routine, xml->envP.fault_string, xml->envP.fault_code);
  }
  
  return out;
} /* end ObitXMLXML2InfoList */

/**
 * Convert XML created from an ObitInfoList to ObitInfoList
 * \param xml     Object to convert
 * \param err     Obit Error message
 * \return  new ObitInfoList (Unref when done)
 */
ObitInfoList*
ObitXMLXMLInfo2List (ObitXML *xml, ObitErr *err)
{
  ObitInfoList  *out=NULL;
  xmlrpc_int32 XMLType;
  xmlrpc_type xmlType;
  gchar *routine = "ObitXMLXMLInfo2List";
  
  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return out;
  g_assert (ObitXMLIsA(xml));
  Obit_retval_if_fail(((xml->type == OBIT_XML_InfoList) || 
		       (xml->type == OBIT_XML_RPCCallArg)), err, out,
		      "%s: xml wrong type %d != %d", 
		      routine, xml->type, OBIT_XML_InfoList);

  /* Be sure xml->parmP a struct of correct type */
  xmlType = xmlrpc_value_type (xml->parmP);
  if (xmlType==XMLRPC_TYPE_STRUCT) {
    if (xmlrpc_struct_has_key(&xml->envP,xml->parmP,"XMLType")) {
      xmlrpc_decompose_value(&xml->envP, xml->parmP, "{s:i,*}",
			     "XMLType", &XMLType);
	XMLRPC_FAIL_IF_FAULT(&xml->envP);
	Obit_retval_if_fail((XMLType == OBIT_XML_InfoList), err, out,
			    "%s: xml wrong type %d != %d", 
			    routine, XMLType, OBIT_XML_InfoList);
   }
  } else { /* Not a struct */
    Obit_log_error(err, OBIT_Error, "%s: input XML NOT a STRUCT", routine);
    return out;
 }
  
  /* Convert */
  decodeInfoList ((ObitXMLEnv*)&xml->envP, (ObitXMLValue*)xml->parmP, &out, err);
  if (err->error) Obit_traceback_val (err, routine, xml->name, out);  

    /* Make sure everything is cool */
 cleanup:
  if (xml->envP.fault_occurred) {
    Obit_log_error(err, OBIT_Error, "%s:XML-RPC Fault: %s (%d)",
		   routine, xml->envP.fault_string, xml->envP.fault_code);
  }
  return out;
} /* end ObitXMLXMLInfo2List */

/**
 * Convert Image definitions information to XML
 * Intended for RPC function "loadImage"
 * type = OBIT_XML_LoadImage
 * \param Type    Image type, OBIT_IO_FITS or OBIT_IO_AIPS
 * \param Name    FITS file path or AIPS Name (12 char)
 * \param AClass  AIPS class (6 char), may be NULL for FITS
 * \param ADir    Path to AIPS directory, may be NULL for FITS
 * \param ASeq    AIPS image sequence number
 * \param AUser   AIPS User number
 * \param Field   Field number (1-rel)
 * \param NField  Total number of fields in mosaic
 * \param err     Obit Error message
 * \return  new ObitXML object
 */
ObitXML* 
ObitXMLFileInfo2XML (ObitIOType Type, gchar *Name,
		     gchar *AClass, gchar *ADir, olong ASeq, olong AUser,
		     olong Field, olong NField, 
		     ObitErr *err)
{
  ObitXML *out=NULL;
  gchar *nix="None";
  gchar *name, *aclass, *adir;
  gchar *routine = "ObitXMLFileInfo2XML";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return out;
  g_assert(Name!=NULL);

   /* initial structure */
  out = newObitXML(routine);
  out->type = OBIT_XML_LoadImage;
  out->func = g_strdup("loadImage");
  
 /* Deal with values not given */
  name = Name;
  if (AClass) aclass = AClass;
  else aclass = nix;
  if (ADir) adir = ADir;
  else adir = nix;
  
  out->parmP = xmlrpc_build_value(&out->envP, "{s:i,s:i,s:i,s:i,s:i,s:i,s:s,s:s,s:s}",
				  "XMLType",  (xmlrpc_int32)OBIT_XML_LoadImage,
				  "DataType", (xmlrpc_int32)Type,
				  "Field",    (xmlrpc_int32)Field,
				  "NField",   (xmlrpc_int32)NField,
				  "ASeq",     (xmlrpc_int32)ASeq,
				  "AUser",    (xmlrpc_int32)AUser,
				  "Name",     name,
				  "AClass",   aclass,
				  "ADir",     adir);
  
  /* Make sure everything is cool */
  if (out->envP.fault_occurred) {
    Obit_log_error(err, OBIT_Error, "%s: XML-RPC Fault: %s (%d)",
		  routine,  out->envP.fault_string, out->envP.fault_code);
  }

  return out;
} /* end ObitXMLFileInfo2XML */

/**
 * Convert an ObitXML to an Image definition
 * \param xml     ObitXML Object from which to extract information
 * \param Type    [out]Image type, OBIT_IO_FITS or OBIT_IO_AIPS
 * \param Name    [out]FITS file path or AIPS Name (12 char)
 *                will be allocated, g_free when done
 * \param AClass  [out]AIPS class (6 char), 
 *                if non NULL will be allocated, g_free when done
 * \param ADir    [out]Path to AIPS directory
 *                if non NULL will be allocated, g_free when done
 * \param Aseq    [out]AIPS image sequence number
 * \param Auser   [out]AIPS User number
 * \param err     Obit Error message
 * \param Field   [out]Field number (1-rel)
 * \param NField  [out]Total number of fields in mosaic
*/
void
ObitXMLXML2FileInfo (ObitXML* xml, ObitIOType *Type, gchar **Name,
		     gchar **AClass, gchar **ADir, olong *ASeq, olong *AUser,
		     olong *Field, olong *NField, 
		     ObitErr *err)
{
  xmlrpc_int32 type, aseq, auser, field, nfield;
  gchar *name, *aclass, *adir;
  gchar *routine = "ObitXMLXML2FileInfo";
  
  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  Obit_return_if_fail((xml->type == OBIT_XML_LoadImage), err, 
		      "%s: xml wrong type %d != %d", 
		      routine, xml->type, OBIT_XML_LoadImage);

  /* parse information */
  xmlrpc_decompose_value(&xml->envP, xml->parmP, 
			 "{s:i,s:i,s:i,s:i,s:i,s:s,s:s,s:s,*}",
			 "DataType", &type,
			 "Field",    &field,
			 "NField",   &nfield,
			 "ASeq",     &aseq,
			 "AUser",    &auser,
			 "Name",     &name,
			 "AClass",   &aclass,
			 "ADir",     &adir);
  
  /* Make sure everything is cool */
  if (xml->envP.fault_occurred) {
    Obit_log_error(err, OBIT_Error, "%s:XML-RPC Fault: %s (%d)",
		   routine, xml->envP.fault_string, xml->envP.fault_code);
    return;
  }

  /* Copy to output */
  *Type   = (ObitIOType)type;
  *ASeq   = (olong)aseq;
  *AUser  = (olong)auser;
  *Field  = (olong)field;
  *NField = (olong)nfield;
  *Name   = g_strdup(name);
  if (AClass) *AClass = g_strdup(aclass);
  if (ADir) *ADir     = g_strdup(adir);
  return;
} /* end ObitXMLXML2FileInfo */

/**
 * Convert one field of an ObitDConCleanWindow to XML
 * Intended for RPC function "editWindow"
 * type = OBIT_XML_EditWindow
 * \param window  Object from which to extract information
 * \param field   Which field? (1-rel)
 * \param err     Obit Error message
 * \return ObitXML object
 */
ObitXML* ObitXMLWindow2XML (ObitDConCleanWindow *window, olong field, 
			    ObitErr *err)
{ 
  ObitXML  *out=NULL;
  xmlrpc_value *v;
  olong *win, iD, nId;
  gchar bname[21];
  ObitDConCleanWindowType type;
  gchar *routine = "ObitXMLWindow2XML";
  
  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return out;
  g_assert (ObitDConCleanWindowIsA(window));

  /* initial structure */
  out = newObitXML(routine);
  out->parmP = xmlrpc_struct_new(&out->envP);
  out->type = OBIT_XML_EditWindow;
  out->func = g_strdup("editWindow");

  /* ObitXML type */
  v = xmlrpc_build_value(&out->envP, "i", (xmlrpc_int32)OBIT_XML_EditWindow);
  xmlrpc_struct_set_value(&out->envP, out->parmP, "XMLType", v);
  xmlrpc_DECREF(v);
 
  /* nfield = only 1 allowed */
  v = xmlrpc_build_value(&out->envP, "i", (xmlrpc_int32)1);
  xmlrpc_struct_set_value(&out->envP, out->parmP, "nfield", v);
  xmlrpc_DECREF(v);
 
  /* ndim */
  v = xmlrpc_build_value(&out->envP, "i", (xmlrpc_int32)window->ndim);
  xmlrpc_struct_set_value(&out->envP, out->parmP, "ndim", v);
  xmlrpc_DECREF(v);
  XMLRPC_FAIL_IF_FAULT(&out->envP);

  /* naxis */
  v = xmlrpc_build_value(&out->envP, "(ii)", 
			 (xmlrpc_int32)window->naxis[field-1][0],
			 (xmlrpc_int32)window->naxis[field-1][1]);
  xmlrpc_struct_set_value(&out->envP, out->parmP, "naxis", v);
  xmlrpc_DECREF(v);
  XMLRPC_FAIL_IF_FAULT(&out->envP);

  /* maxId */
  v = xmlrpc_build_value(&out->envP, "i", (xmlrpc_int32)window->maxId[field-1]);
  xmlrpc_struct_set_value(&out->envP, out->parmP, "maxId", v);
  xmlrpc_DECREF(v);
  XMLRPC_FAIL_IF_FAULT(&out->envP);

  /* Outer Window */
  iD=-1;
  if (ObitDConCleanWindowInfo(window, field, iD, &type, &win, err)) {
    /* Which box type */
    if ((type==OBIT_DConCleanWindow_rectangle) || 
	(type==OBIT_DConCleanWindow_unrectangle))
      v = xmlrpc_build_value(&out->envP, "(iiiii)", 
			     (xmlrpc_int32)type,
			     (xmlrpc_int32)win[0],(xmlrpc_int32)win[1],
			     (xmlrpc_int32)win[2],(xmlrpc_int32)win[3]);
    if ((type==OBIT_DConCleanWindow_round) ||
	(type==OBIT_DConCleanWindow_unround))
      v = xmlrpc_build_value(&out->envP, "(iiii)", 
			     (xmlrpc_int32)type,
			     (xmlrpc_int32)win[0],(xmlrpc_int32)win[1],
			     (xmlrpc_int32)win[2]);
    XMLRPC_FAIL_IF_FAULT(&out->envP);
    sprintf (bname, "outerWin");
    xmlrpc_struct_set_value(&out->envP, out->parmP, bname, v);
    xmlrpc_DECREF(v);
    XMLRPC_FAIL_IF_FAULT(&out->envP);
  }  /* end ifouter window */

  /* How many potential inner windows? */
  nId = window->maxId[field-1];
  /* Extract from Obit Object and encode */
  for (iD=1; iD<=nId; iD++) {
    if (ObitDConCleanWindowInfo(window, field, iD, &type, &win, err)) {
      /* Which box type */
      if ((type==OBIT_DConCleanWindow_rectangle) || 
	  (type==OBIT_DConCleanWindow_unrectangle))
	v = xmlrpc_build_value(&out->envP, "(iiiii)", 
			       (xmlrpc_int32)type,
			       (xmlrpc_int32)win[0],(xmlrpc_int32)win[1],
			       (xmlrpc_int32)win[2],(xmlrpc_int32)win[3]);
      if ((type==OBIT_DConCleanWindow_round) ||
	  (type==OBIT_DConCleanWindow_unround))
	v = xmlrpc_build_value(&out->envP, "(iiii)", 
			       (xmlrpc_int32)type,
			       (xmlrpc_int32)win[0],(xmlrpc_int32)win[1],
			       (xmlrpc_int32)win[2]);
      XMLRPC_FAIL_IF_FAULT(&out->envP);
      sprintf (bname, "win%4.4d", iD);
      xmlrpc_struct_set_value(&out->envP, out->parmP, bname, v);
      xmlrpc_DECREF(v);
      XMLRPC_FAIL_IF_FAULT(&out->envP);

    }  /* end if valid */
  } /* end Loop over windows */
  
  /* Make sure everything is cool */
  cleanup:
  if (out->envP.fault_occurred) {
    Obit_log_error(err, OBIT_Error, "%s:XML-RPC Fault: %s (%d)",
		   routine, out->envP.fault_string, out->envP.fault_code);
  }
  
  return out;
} /* end ObitXMLWindow2XML */

/**
 * Convert an ObitXML to an ObitDConCleanWindow
 * \param env        xmlrpc environment
 * \param xmlwindow  XML Object from which to extract information
 * \param err        Obit Error message
 * \return Obit Window structure
 */
ObitDConCleanWindow* ObitXMLXML2Window (ObitXML* xml, ObitErr *err)
{ 
  ObitDConCleanWindow *out=NULL;
  olong iD, nId, nax[2], box[4], field=1;
  gchar bname[21];
  olong indx, i;
  xmlrpc_int32 nfield, ndim, naxis[2], maxId, bType;
  xmlrpc_value *v, *vv;
  xmlrpc_type xmlType;
  ObitDConCleanWindowType type=0;
  gchar *routine = "ObitXMLXML2Window";
 
  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return out;
  /* Could be parsing the reply from an RPC Call */
  Obit_retval_if_fail(((xml->type == OBIT_XML_EditWindow) || 
		       (xml->type == OBIT_XML_Reply)), 
		      err, out, "%s: xml wrong type %d != %d", 
		      routine, xml->type, OBIT_XML_EditWindow);

  /* Be sure xml->parmP a struct of correct type */
  xmlType = xmlrpc_value_type (xml->parmP);
  Obit_retval_if_fail((xmlType==XMLRPC_TYPE_STRUCT),
		      err, out, "%s: return NOT a struct", routine);
  Obit_retval_if_fail((xmlrpc_struct_has_key(&xml->envP,xml->parmP,"nfield")),
		      err, out, "%s: xml wrong type", routine);
  XMLRPC_FAIL_IF_FAULT(&xml->envP);
 
  /* parse header information */
  xmlrpc_decompose_value(&xml->envP, xml->parmP, "{s:i,s:i,s:(ii),s:i,*}",
			 "nfield", &nfield,
			 "ndim",   &ndim,
			 "naxis",  &naxis[0], &naxis[1],
			 "maxId",  &maxId);
  XMLRPC_FAIL_IF_FAULT(&xml->envP);
  nax[0] = (olong)naxis[0];
  nax[1] = (olong)naxis[1];

  /* initial structure */
  out = ObitDConCleanWindowCreate1("CleanWindow", nax, err);

  out->nfield      = (olong)nfield;
  out->ndim        = (olong)ndim;
  out->maxId[0]    = 0;

  /* Outer window */
  sprintf (bname, "outerWin");
  indx = xmlrpc_struct_has_key (&xml->envP, xml->parmP, bname);
  XMLRPC_FAIL_IF_FAULT(&xml->envP);
  if (indx!=-1) { /* it exists */
    xmlrpc_struct_find_value (&xml->envP, xml->parmP, bname, &v);
    XMLRPC_FAIL_IF_FAULT(&xml->envP);
    /* xmlrpc_struct_has_key may lie */
    if (v) {
      /* Size depends on type = first word */
      xmlrpc_array_read_item (&xml->envP, v, 0, &vv);
      XMLRPC_FAIL_IF_FAULT(&xml->envP);
      xmlrpc_decompose_value(&xml->envP, vv, "i", &bType);
      xmlrpc_DECREF(vv);
      XMLRPC_FAIL_IF_FAULT(&xml->envP);
      /* by type */
      if ((bType==OBIT_DConCleanWindow_rectangle) || 
	  (bType==OBIT_DConCleanWindow_unrectangle)) { /* rectangle - 4 integers */
	type = bType;
	for (i=1; i<=4; i++) {
	  xmlrpc_array_read_item (&xml->envP, v, i, &vv);
	  XMLRPC_FAIL_IF_FAULT(&xml->envP);
	  xmlrpc_decompose_value(&xml->envP, vv, "i", &bType);
	  box[i-1] = (olong)bType;
	  xmlrpc_DECREF(vv);
	  XMLRPC_FAIL_IF_FAULT(&xml->envP);
	}
      } else if ((bType==OBIT_DConCleanWindow_round) || 
		 (bType==OBIT_DConCleanWindow_unround)) { /* round - 3 integers */
	type = bType;
	for (i=1; i<=3; i++) {
	  xmlrpc_array_read_item (&xml->envP, v, i, &vv);
	  XMLRPC_FAIL_IF_FAULT(&xml->envP);
	  xmlrpc_decompose_value(&xml->envP, vv, "i", &bType);
	  box[i-1] = (olong)bType;
	  xmlrpc_DECREF(vv);
	  XMLRPC_FAIL_IF_FAULT(&xml->envP);
	}
      }
      /* Add it */
      ObitDConCleanWindowOuter (out, field, type, box, err);
    } /* end if it really exists */
  } /* end if exists */
  if (xml->envP.fault_occurred) goto cleanup;  /* something go wrong? */

  /* How many potential inner windows? */
  nId = maxId;
  type = OBIT_DConCleanWindow_rectangle;
  /* Loop over potential entries */
  for (iD=1; iD<=nId; iD++) {
    sprintf (bname, "win%4.4d", iD);
    indx = xmlrpc_struct_has_key (&xml->envP, xml->parmP, bname);
    XMLRPC_FAIL_IF_FAULT(&xml->envP);
    if (indx!=-1) { /* it exists */
      xmlrpc_struct_find_value (&xml->envP, xml->parmP, bname, &v);
      XMLRPC_FAIL_IF_FAULT(&xml->envP);
      /* xmlrpc_struct_has_key may lie */
      if (v) {
	/* Size depends on type = first word */
	xmlrpc_array_read_item (&xml->envP, v, 0, &vv);
	XMLRPC_FAIL_IF_FAULT(&xml->envP);
	xmlrpc_decompose_value(&xml->envP, vv, "i", &bType);
	xmlrpc_DECREF(vv);
	XMLRPC_FAIL_IF_FAULT(&xml->envP);
	/* by type */
	if ((bType==OBIT_DConCleanWindow_rectangle) || 
	    (bType==OBIT_DConCleanWindow_unrectangle)) { /* rectangle - 4 integers */
	  type = bType;
	  for (i=1; i<=4; i++) {
	    xmlrpc_array_read_item (&xml->envP, v, i, &vv);
	    XMLRPC_FAIL_IF_FAULT(&xml->envP);
	    xmlrpc_decompose_value(&xml->envP, vv, "i", &bType);
	    box[i-1] = (olong)bType;
	    xmlrpc_DECREF(vv);
	    XMLRPC_FAIL_IF_FAULT(&xml->envP);
	  }
	} else if ((bType==OBIT_DConCleanWindow_round) || 
		   (bType==OBIT_DConCleanWindow_unround)) { /* round - 3 integers */
	  type = bType;
	  for (i=1; i<=3; i++) {
	    xmlrpc_array_read_item (&xml->envP, v, i, &vv);
	    XMLRPC_FAIL_IF_FAULT(&xml->envP);
	    xmlrpc_decompose_value(&xml->envP, vv, "i", &bType);
	    box[i-1] = (olong)bType;
	    xmlrpc_DECREF(vv);
	    XMLRPC_FAIL_IF_FAULT(&xml->envP);
	  }
	}
	/* Add it */
	ObitDConCleanWindowAdd (out, field, type, box, err);
      } /* end if it really exists */
    } /* end if exists */
    if (xml->envP.fault_occurred) break;  /* something go wrong? */
  } /* end loop */
      

  /* Make sure everything is cool */
  cleanup:
  if (xml->envP.fault_occurred) {
    Obit_log_error(err, OBIT_Error, "XML-RPC Fault: %s (%d)",
		   xml->envP.fault_string, xml->envP.fault_code);
  }

  return out;
} /* end ObitXMLXML2Window */

/**
 * Convert Binary blob with description in InfoList to XML
 * NB: no function name is set on output.
 * type = OBIT_XML_BinBlob
 * \param blob    Binary blob to include
 * \param desc    Description of blob including:
 * \li size  long scalar    Size of blob in bytes
 * \param err     Obit Error message
 * \return  new ObitXML, blob in member "BlobData"
 */
ObitXML* 
ObitXMLBlob2XML (gpointer blob, ObitInfoList *desc, ObitErr *err)
{
  ObitXML  *out=NULL;
  xmlrpc_value *v;
  ObitInfoType infoType;
  gint32 dim[MAXINFOELEMDIM];
  olong size;
  gchar *routine = "ObitXMLBlob2XML";
  
  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return out;

  /* initial structure */
  out = newObitXML(routine);
  out->parmP = xmlrpc_struct_new(&out->envP);
  out->type = OBIT_XML_BinBlob;
  out->func = g_strdup("undefined");

  /* ObitXML type */
  v = xmlrpc_build_value(&out->envP, "i", (xmlrpc_int32)OBIT_XML_BinBlob);
  xmlrpc_struct_set_value(&out->envP, out->parmP, "XMLType", v);
  xmlrpc_DECREF(v);
 
  /* Copy description */
  encodeInfoList (desc, (ObitXMLEnv*)&out->envP, (ObitXMLValue*)out->parmP, err);
  if (err->error) Obit_traceback_val (err, routine, out->name, out);  

  /* get size of blob */
  ObitInfoListGet(desc, "size", &infoType, dim, &size, err);
  if (err->error) Obit_traceback_val (err, routine, out->name, out);  

  /* encode blob */
  v = xmlrpc_base64_new(&out->envP, size, blob);
  xmlrpc_struct_set_value(&out->envP, out->parmP, "BlobData", v);
  xmlrpc_DECREF(v);

  /* Make sure everything is cool */
  if (out->envP.fault_occurred) {
    Obit_log_error(err, OBIT_Error, "%s: XML-RPC Fault: %s (%d)",
		  routine,  out->envP.fault_string, out->envP.fault_code);
  }

  return out;
} /* end  ObitXMLBlob2XML */

/**
 * Convert XML to Binary blob with description in InfoList
 * \param xml     Object to convert
 * \param desc    [out] Description of output blob including:
 * \li size  olong scalar  Size of blob in bytes
 * \param err     Obit Error message
 * \return  new output binary blob (ObitMemFree when done)
 */
gpointer 
ObitXMLXML2Blob (ObitXML *xml, ObitInfoList **desc, ObitErr *err)
{
  gpointer out=NULL;
  const unsigned char* xmlrpc_block=NULL;
  olong size;
  ObitInfoType infoType;
  gint32 dim[MAXINFOELEMDIM];
  xmlrpc_value *v;
  xmlrpc_type xmlType;
  size_t usize;
  gchar *routine = "ObitXMLXML2Blob";
 
  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return out;
  /* Could be parsing the reply from an RPC Call */
  Obit_retval_if_fail(((xml->type == OBIT_XML_BinBlob) || 
		       (xml->type == OBIT_XML_Reply)),
		      err, out, "%s: xml wrong type %d != %d", 
		      routine, xml->type, OBIT_XML_BinBlob);

  /* Be sure xml->parmP a struct of correct type */
  xmlType = xmlrpc_value_type (xml->parmP);
  Obit_retval_if_fail((xmlType==XMLRPC_TYPE_STRUCT),
		      err, out, "%s: return NOT a struct", routine);
  XMLRPC_FAIL_IF_FAULT(&xml->envP);
 
  /* Parse descriptive info */
  decodeInfoList ((ObitXMLEnv*)&xml->envP, (ObitXMLValue*)xml->parmP, desc, err);
  if (err->error) Obit_traceback_val (err, routine, xml->name, out); 

  /* How big is the blob */
  ObitInfoListGet(*desc, "size", &infoType, dim, &size, err);
  if (err->error) Obit_traceback_val (err, routine, xml->name, out);  

  /* allocate blob */
  out = ObitMemAllocName ((gulong)size, "XML:BinBlob");
  Obit_retval_if_fail((out!=NULL),
		      err, out, "%s: memory allocation for %d failed", routine, size);

  /* Parse blob */
  xmlrpc_struct_find_value (&xml->envP, xml->parmP,  "BlobData", &v);
  if (v) xmlrpc_read_base64 (&xml->envP, v, &usize, &xmlrpc_block);
  /* Make sure everything is cool */
  if (xml->envP.fault_occurred || (xmlrpc_block==NULL)) {
    Obit_log_error(err, OBIT_Error, "%s: XML-RPC Fault: %s (%d)",
		  routine,  xml->envP.fault_string, xml->envP.fault_code);
  }
  
  /* Copy to output */
  memcpy(out, xmlrpc_block, usize);
  free ((gchar*)xmlrpc_block);

    /* Make sure everything is cool */
 cleanup:
  if (xml->envP.fault_occurred || (xmlrpc_block==NULL)) {
    Obit_log_error(err, OBIT_Error, "%s: XML-RPC Fault: %s (%d)",
		  routine,  xml->envP.fault_string, xml->envP.fault_code);
  }
  
  return out;

} /* end ObitXMLXML2Blob */

/**
 * Create an ObitXML for the return value from an RPC call
 * \param name  Name string for new object
 * \param parmP xmlrpc_value (as gpointer) as xml data for new object
 *              This passes control of parmP to the new object
 * \param err     Obit Error message
 * \return new ObitXML object
 */
ObitXML* 
ObitXMLReturn (gchar *name, gpointer parmP, ObitErr *err)
{
  ObitXML  *out=NULL;
  xmlrpc_int32 XMLType;
  xmlrpc_type xmlType;
  gchar *routine = "ObitXMLReturn";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return out;
  g_assert (parmP!=NULL);

  /* initial structure */
  out = newObitXML(name);
  out->type  = OBIT_XML_Reply;
  out->func  = g_strdup(routine);
  out->parmP = (xmlrpc_value*)parmP;
  xmlrpc_INCREF(out->parmP);

  /* See if xml type in xml */
  /* Be sure out->parmP a struct of correct type */
  xmlType = xmlrpc_value_type (out->parmP);
  if (xmlType==XMLRPC_TYPE_STRUCT) {
    if (xmlrpc_struct_has_key(&out->envP,out->parmP,"Result")) {
      /* DEBUG */
      out->type  = OBIT_XML_Reply;
    }
    if (xmlrpc_struct_has_key(&out->envP,out->parmP,"XMLType")) {
      xmlrpc_decompose_value(&out->envP, out->parmP, "{s:i,*}",
			     "XMLType", &XMLType);
	/* save type on object */
	out->type  = (ObitXMLType)XMLType;
    }
  }
  
  
  /* Make sure everything is cool */
  if (out->envP.fault_occurred) {
    Obit_log_error(err, OBIT_Error, "%s:XML-RPC Fault: %s (%d)",
		   routine, out->envP.fault_string, out->envP.fault_code);
  }
  
  return out;
} /* end ObitXMLReturn */

/**
 * Create an ObitXML for the return value from an server RPC call
 * \param list    Info list to be written as xml entry "Result"
 * \param err     Obit Error message
 * \return new ObitXML object
 */
ObitXML* 
ObitXMLServerReturn (ObitInfoList *list, ObitErr *err)
{
  ObitXML  *out=NULL;
  xmlrpc_value *v;
  gchar *routine = "ObitXMLServerReturn";

  /* error checks */
  if (err->error) return out;
  g_assert (ObitInfoListIsA(list));

  /* initial structure */
  out = newObitXML(routine);
  out->type  = OBIT_XML_InfoList;
  out->func  = g_strdup(routine);
  out->parmP = xmlrpc_struct_new(&out->envP);

  /* ObitXML type */
  v = xmlrpc_build_value(&out->envP, "i", (xmlrpc_int32)OBIT_XML_InfoList);
  xmlrpc_struct_set_value(&out->envP, out->parmP, "XMLType", v);
  xmlrpc_DECREF(v);

  /* Extract from Obit Object and encode */
  v = xmlrpc_struct_new(&out->envP);
  encodeInfoList (list, (ObitXMLEnv*)&out->envP, (ObitXMLValue*)v, err);
  if (err->error) Obit_traceback_val (err, routine, out->name, out);  
  
  /* Add to out */
  xmlrpc_struct_set_value(&out->envP, out->parmP, "Result", v);
  xmlrpc_DECREF(v);

  /* Make sure everything is cool */
  if (out->envP.fault_occurred) {
    Obit_log_error(err, OBIT_Error, "%s:XML-RPC Fault: %s (%d)",
		   routine, out->envP.fault_string, out->envP.fault_code);
  }
  
  return out;
} /* end ObitXMLServerReturn */

/**
 * Extract InfoList from "Result" in XML object from ObitRPCCall
 * This to be used on the client side to extract the returned data.
 * \param xml   Object from which result to be extracted
 * \param err   Obit Error message
 * \return new  Info 
 */
ObitInfoList* 
ObitXMLGetServerResult (ObitXML *xml,  ObitErr *err)
{
  ObitInfoList  *out=NULL;
  xmlrpc_type xmlType;
  gchar *routine = "ObitXMLGetServerResult";

  /* error checks */
  if (err->error) return out;
  g_assert (ObitXMLIsA(xml));
  
    /* Be sure xml->parmP a struct */
  xmlType = xmlrpc_value_type (xml->parmP);
  Obit_retval_if_fail(((xml->type == OBIT_XML_InfoList) || 
		       (xml->type == OBIT_XML_Reply)), err, out,
		      "%s: xml wrong type %d != %d", 
		      routine, xml->type, OBIT_XML_InfoList);
  
  /* initial structure */
  out = newObitInfoList();

  decodeInfoList (&xml->envP, xml->parmP, &out, err);
  XMLRPC_FAIL_IF_FAULT(&xml->envP);
  if (err->error) Obit_traceback_val (err, routine, xml->name, out); 

  /* Make sure everything is cool */
 cleanup:
  if (xml->envP.fault_occurred) {
    Obit_log_error(err, OBIT_Error, "%s:XML-RPC Fault: %s (%d)",
		   routine, xml->envP.fault_string, xml->envP.fault_code);
  }
  
  return out;
} /* end ObitXMLGetServerResult */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitXMLClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitXMLClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitXMLClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitXMLClassInfoDefFn (gpointer inClass)
{
  ObitXMLClassInfo *theClass = (ObitXMLClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitXMLClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitXMLClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitXMLGetClass;
  theClass->newObit       = (newObitFP)newObitXML;
  theClass->ObitCopy      = NULL;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitXMLClear;
  theClass->ObitInit      = (ObitInitFP)ObitXMLInit;
  
} /* end ObitXMLClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitXMLInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitXML *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);
  
  /* set members in this class */
  in->thread = newObitThread();
  /* Initialize our error-handling environment. */
  xmlrpc_env_init(&in->envP);
  in->parmP = NULL;
  in->func  = NULL;
  
} /* end ObitXMLInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitXML* cast to an Obit*.
 */
void ObitXMLClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitXML *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  in->thread  = ObitThreadUnref(in->thread);
  xmlrpc_env_clean(&in->envP);
  xmlrpc_DECREF(in->parmP);
  if (in->func) g_free(in->func); in->func = NULL;
 
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitXMLClear */

/**
 * Low level convert ObitInfoList to XML
 * Note: Can only translate structures with only data types
 * directly translatable into XML.
 * \param list    List to convert
 * \param envP    XML environment
 * \param parmP   XML value
 * \param err     Obit Error message
 */
static void encodeInfoList (ObitInfoList *list, 
			    ObitXMLEnv *envP, ObitXMLValue *parmP,
			    ObitErr *err)
{
  gchar *nameP;
  ObitInfoType infoType;
  gint32 dim[MAXINFOELEMDIM];
  olong i, j, num;
  gpointer data;
  xmlrpc_value *v, *d, *e;
  gchar bname[21];
  gboolean badType=FALSE;
  gshort *shortP;
  gboolean *booleanP;
  olong   *intP;
  olong  *longP;
  oint   *ointP;
  ofloat *floatP;
  odouble *doubleP;
  gchar *routine = "ObitXML:encodeInfoList";

  /* error checks */
  if (err->error) return;

  /* number of entries */
  v = xmlrpc_build_value(envP, "i", (xmlrpc_int32)list->number);
  xmlrpc_struct_set_value(envP, parmP, "number", v);
  xmlrpc_DECREF(v);
  XMLRPC_FAIL_IF_FAULT(envP);

  /* Extract from Obit Object and encode */
  for (i=1; i<=list->number; i++) {
    if (ObitInfoListGetNumberP(list, (olong)i, &nameP, &infoType, dim, &data)) {
      /* Header for entry */
      v = xmlrpc_build_value(envP, "{s:s,s:i,s:(iiiii)}", 
			     "name",nameP,
			     "type",(xmlrpc_int32)infoType,
			     "dim",(xmlrpc_int32)dim[0],(xmlrpc_int32)dim[1],
			     (xmlrpc_int32)dim[2],(xmlrpc_int32)dim[3], (xmlrpc_int32)dim[4]);
      XMLRPC_FAIL_IF_FAULT(envP);

      /* How much data? */
      num = MAX (1, dim[0]);
      for (j=1; j<MAXINFOELEMDIM; j++) num *= MAX (1,dim[j]);

      /* create array to add */
      d = xmlrpc_array_new(envP);
      XMLRPC_FAIL_IF_FAULT(envP);

      /* Add data by type */
      switch (infoType) { 
      case OBIT_byte:
	badType = FALSE;
	break;
      case OBIT_short:
	shortP = (gshort*)data;
	for (j=0; j<num; j++) {
	  e = xmlrpc_build_value(envP, "i", (xmlrpc_int32)(*shortP++));
	  xmlrpc_array_append_item (envP, d, e);
	  xmlrpc_DECREF(e);
	}
	break;
      case OBIT_int:
	intP = (olong*)data;
	for (j=0; j<num; j++) {
	  e = xmlrpc_build_value(envP, "i", (xmlrpc_int32)(*intP++));
	  xmlrpc_array_append_item (envP, d, e);
	  xmlrpc_DECREF(e);
	}
	break;
      case OBIT_oint:
	ointP = (oint*)data;
	for (j=0; j<num; j++) {
	  e = xmlrpc_build_value(envP, "i", (xmlrpc_int32)(*ointP++));
	  xmlrpc_array_append_item (envP, d, e);
	  xmlrpc_DECREF(e);
	}
	break;
      case OBIT_long:
	longP = (olong*)data;
	for (j=0; j<num; j++) {
	  e = xmlrpc_build_value(envP, "i", (xmlrpc_int32)(*longP++));
	  xmlrpc_array_append_item (envP, d, e);
	  xmlrpc_DECREF(e);
	}
	break;
      case OBIT_ubyte:
	badType = FALSE;
	break;
      case OBIT_ushort:
	badType = FALSE;
	break;
      case OBIT_uint:
	badType = FALSE;
	break;
      case OBIT_ulong:
	badType = FALSE;
	break;
      case OBIT_float:
	floatP = (ofloat*)data;
	for (j=0; j<num; j++) {
	  e = xmlrpc_build_value(envP, "d", (xmlrpc_double)(*floatP++));
	  xmlrpc_array_append_item (envP, d, e);
	  xmlrpc_DECREF(e);
	}
	break;
      case OBIT_double:
	doubleP = (odouble*)data;
	for (j=0; j<num; j++) {
	  e = xmlrpc_build_value(envP, "d", (xmlrpc_double)(*doubleP++));
	  xmlrpc_array_append_item (envP, d, e);
	  xmlrpc_DECREF(e);
	}
	break;
      case OBIT_complex:
	floatP = (ofloat*)data;
	for (j=0; j<num*2-1; j++) {
	  e = xmlrpc_build_value(envP, "d", (xmlrpc_double)(*floatP++));
	  xmlrpc_array_append_item (envP, d, e);
	  xmlrpc_DECREF(e);
	}
	break;
      case OBIT_dcomplex:
	doubleP = (odouble*)data;
	for (j=0; j<num*2-1; j++) {
	  e = xmlrpc_build_value(envP, "d", (xmlrpc_double)(*doubleP++));
	  xmlrpc_array_append_item (envP, d, e);
	  xmlrpc_DECREF(e);
	}
	break;
      case OBIT_string:
	e = xmlrpc_build_value(envP, "s#", (char*)data, num);
	xmlrpc_array_append_item (envP, d, e);
	xmlrpc_DECREF(e);
	break;
      case OBIT_bool:
	booleanP = (gboolean*)data;
	for (j=0; j<num; j++) {
	  e = xmlrpc_build_value(envP, "b", (xmlrpc_bool)(*booleanP++));
	  xmlrpc_array_append_item (envP, d, e);
	  xmlrpc_DECREF(e);
	}
	break;
      case OBIT_bits:
	badType = FALSE;
      default:
	g_assert_not_reached(); /* unknown, barf */
      }; /* end switch  */
      XMLRPC_FAIL_IF_FAULT(envP);

      /* Add data */
      xmlrpc_struct_set_value(envP, v, "data", d);
      xmlrpc_DECREF(d);
      XMLRPC_FAIL_IF_FAULT(envP);

      /* Add entry to output */
      sprintf (bname, "entry%4.4d", i);
      xmlrpc_struct_set_value(envP, parmP, bname, v);
      xmlrpc_DECREF(v);
      XMLRPC_FAIL_IF_FAULT(envP);
    } /* end if found */

    /* Unsupported data type? */
    Obit_return_if_fail((!badType), err,
		      "%s: Unsupported InfoList type in xml  %d", 
		      routine, infoType);

  } /* end loop over list */
  
  /* Make sure everything is cool */
  cleanup:
  if (envP->fault_occurred) {
    Obit_log_error(err, OBIT_Error, "%s:XML-RPC Fault: %s (%d)",
		   routine, envP->fault_string, envP->fault_code);
  }
  
  
 } /*  end encodeInfoList */

/**
 * Low level convert XML to ObitInfoList
 * \param envP    XML environment
 * \param parmP   XML value
 * \param out     [out]List to accept values (created here)
 * \param err     Obit Error message
 */
static void decodeInfoList (ObitXMLEnv *envP, ObitXMLValue *parmP, 
			    ObitInfoList **out, ObitErr *err)
{
  gchar *name;
  ObitInfoType infoType;
  xmlrpc_int32 i32temp, i32infoType=0, i32dim[MAXINFOELEMDIM];
  xmlrpc_double xmlDouble;
  xmlrpc_bool xmlBool;
  const char *xmlChar;
  gboolean badType=FALSE;
  gshort *shortP;
  gboolean *booleanP;
  olong   *intP;
  olong  *longP;
  oint   *ointP;
  ofloat *floatP;
  odouble *doubleP;
  gint32 dim[MAXINFOELEMDIM];
  gpointer data=NULL;
  xmlrpc_value *v=NULL, *d=NULL, *e=NULL;
  olong i, j, num, indx, number;
  gchar bname[21];
  gchar *routine = "ObitXML:decodeInfoList";
  
  /* error checks */
  if (err->error) return;

  /* How many entries? All may not be from an InfoList */
  number = 0;
  if (xmlrpc_struct_has_key(envP,parmP,"number")) {
    xmlrpc_decompose_value(envP, parmP, "{s:i,*}",
			   "number", &i32temp);
    XMLRPC_FAIL_IF_FAULT(envP);
    number = i32temp;
   }

  /* initial structure */
  *out = newObitInfoList();

  /* Extract from XML Object and copy to InfoList */
  for (i=1; i<=number; i++) {
    sprintf (bname, "entry%4.4d", i);
    indx = xmlrpc_struct_has_key(envP,parmP,bname);
    if (indx!=-1) { /* it exists */
      xmlrpc_struct_find_value (envP, parmP, bname, &v);
      XMLRPC_FAIL_IF_FAULT(envP);
      /* xmlrpc_struct_has_key may lie */
      if (v) {
	/* Header */
	xmlrpc_decompose_value(envP, v, "{s:s,s:i,s:(iiiii),*}", 
			       "name", &xmlChar,
			       "type", &i32infoType,
			       "dim", &i32dim[0],&i32dim[1],&i32dim[2],&i32dim[3],&i32dim[4]);
	XMLRPC_FAIL_IF_FAULT(envP);
	name = (gchar*)xmlChar;
	infoType = (ObitInfoType)i32infoType;
	dim[0] = (olong)i32dim[0]; dim[1] = (olong)i32dim[1]; dim[2] = (olong)i32dim[2]; 
	dim[3] = (olong)i32dim[3]; dim[4] = (olong)i32dim[4];
	
	/* How much data? */
	num = MAX (1, dim[0]);
	for (j=1; j<MAXINFOELEMDIM; j++) num *= MAX (1,dim[j]);

	/* Get data array */
	xmlrpc_struct_find_value (envP, v, "data", &d);
	
	/* Get data by type */
	switch (infoType) { 
	case OBIT_byte:
	  badType = FALSE;
	  break;
	case OBIT_short:
	  data = g_malloc (num*sizeof(gshort));
	  shortP = (gshort*)data;
	  for (j=0; j<num; j++) {
	    e = xmlrpc_array_get_item(envP, d, (int)j);
	    XMLRPC_FAIL_IF_FAULT(envP);
	    xmlrpc_read_int (envP, e, &i32temp);
	    XMLRPC_FAIL_IF_FAULT(envP);
	    (*shortP++) = (gshort)i32temp;
	    xmlrpc_DECREF(e);
	  }
	  break;
	case OBIT_int:
	  data = g_malloc (num*sizeof(olong));
	  intP = (olong*)data;
	  for (j=0; j<num; j++) {
	    e = xmlrpc_array_get_item(envP, d, (int)j);
	    XMLRPC_FAIL_IF_FAULT(envP);
	    xmlrpc_read_int (envP, e, &i32temp);
	    XMLRPC_FAIL_IF_FAULT(envP);
	    (*intP++) = (olong)i32temp;
	    xmlrpc_DECREF(e);
	  }
	  break;
	case OBIT_oint:
	  data = g_malloc (num*sizeof(oint));
	  ointP = (oint*)data;
	  for (j=0; j<num; j++) {
	    e = xmlrpc_array_get_item(envP, d, (int)j);
	    XMLRPC_FAIL_IF_FAULT(envP);
	    xmlrpc_read_int (envP, e, &i32temp);
	    XMLRPC_FAIL_IF_FAULT(envP);
	    (*ointP++) = (olong)i32temp;
	    xmlrpc_DECREF(e);
	  }
	  break;
	case OBIT_long:
	  data = g_malloc (num*sizeof(olong));
	  longP = (olong*)data;
	  for (j=0; j<num; j++) {
	    e = xmlrpc_array_get_item(envP, d, (int)j);
	    XMLRPC_FAIL_IF_FAULT(envP);
	    xmlrpc_read_int (envP, e, &i32temp);
	    XMLRPC_FAIL_IF_FAULT(envP);
	    (*longP++) = (olong)i32temp;
	    xmlrpc_DECREF(e);
	  }
	  break;
	case OBIT_ubyte:
	  badType = FALSE;
	  break;
	case OBIT_ushort:
	  badType = FALSE;
	  break;
	case OBIT_uint:
	  badType = FALSE;
	  break;
	case OBIT_ulong:
	  badType = FALSE;
	  break;
	case OBIT_float:
	  data = g_malloc (num*sizeof(ofloat));
	  floatP = (ofloat*)data;
	  for (j=0; j<num; j++) {
	    e = xmlrpc_array_get_item(envP, d, (int)j);
	    XMLRPC_FAIL_IF_FAULT(envP);
	    xmlrpc_read_double (envP, e, &xmlDouble);
	    XMLRPC_FAIL_IF_FAULT(envP);
	    (*floatP++) = (ofloat)xmlDouble;
	    xmlrpc_DECREF(e);
	  }
	  break;
	case OBIT_double:
	  data = g_malloc (num*sizeof(odouble));
	  doubleP = (odouble*)data;
	  for (j=0; j<num; j++) {
	    e = xmlrpc_array_get_item(envP, d, (int)j);
	    XMLRPC_FAIL_IF_FAULT(envP);
	    xmlrpc_read_double (envP, e, &xmlDouble);
	    XMLRPC_FAIL_IF_FAULT(envP);
	    (*doubleP++) = (odouble)xmlDouble;
	    xmlrpc_DECREF(e);
	  }
	  break;
	case OBIT_complex:
	  data = g_malloc (2*num*sizeof(ofloat));
	  floatP = (ofloat*)data;
	  for (j=0; j<num*2-1; j++) {
	    e = xmlrpc_array_get_item(envP, d, (int)j);
	    XMLRPC_FAIL_IF_FAULT(envP);
	    xmlrpc_read_double (envP, e, &xmlDouble);
	    XMLRPC_FAIL_IF_FAULT(envP);
	    (*floatP++) = (ofloat)xmlDouble;
	    xmlrpc_DECREF(e);
	  }
	  break;
	case OBIT_dcomplex:
	  data = g_malloc (2*num*sizeof(ofloat));
	  doubleP = (odouble*)data;
	  for (j=0; j<num*2-1; j++) {
	    e = xmlrpc_array_get_item(envP, d, (int)j);
	    XMLRPC_FAIL_IF_FAULT(envP);
	    xmlrpc_read_double (envP, e, &xmlDouble);
	    XMLRPC_FAIL_IF_FAULT(envP);
	    (*doubleP++) = (odouble)xmlDouble;
	    xmlrpc_DECREF(e);
	  }
	  break;
	case OBIT_string:
	  data = g_malloc (num*sizeof(gchar)+3);
	  /* all one big happy string */
	  e = xmlrpc_array_get_item(envP, d, (int)0);
	  xmlrpc_read_string (envP, e, &xmlChar);
	  strncpy ((char*)data, xmlChar, num);
	  XMLRPC_FAIL_IF_FAULT(envP);
	  xmlrpc_DECREF(e);
	  break;
	case OBIT_bool:
	  data = g_malloc (num*sizeof(gboolean));
	  booleanP = (gboolean*)data;
	  for (j=0; j<num; j++) {
	    e = xmlrpc_array_get_item(envP, d, (int)j);
	    XMLRPC_FAIL_IF_FAULT(envP);
	    xmlrpc_read_bool (envP, e, &xmlBool);
	    XMLRPC_FAIL_IF_FAULT(envP);
	    (*booleanP++) = (olong)xmlBool;
	    xmlrpc_DECREF(e);
	  }
	  break;
	case OBIT_bits:
	  badType = FALSE;
	  break;
	default:
	  g_assert_not_reached(); /* unknown, barf */
	}; /* end switch  */
	XMLRPC_FAIL_IF_FAULT(envP);
	
	/* Unsupported data type? */
	Obit_return_if_fail((!badType), err,
			    "%s: Unsupported InfoList type in xml  %d", 
			    routine, infoType);
	
	/* Add to ObitInfoList */
	xmlrpc_DECREF(v);
	ObitInfoListPut(*out, name, infoType, dim, data, err);
	g_free(data);
	if (err->error)Obit_traceback_msg (err, routine, "writing ObitInfoList");
	
      } /* end if really found */
    } /* end if found */
  } /* end loop over list */
    
  /* Make sure everything is cool */
 cleanup:
  if (envP->fault_occurred) {
    Obit_log_error(err, OBIT_Error, "%s:XML-RPC Fault: %s (%d)",
		   routine, envP->fault_string, envP->fault_code);
  }
  
} /*  end decodeInfoList */


