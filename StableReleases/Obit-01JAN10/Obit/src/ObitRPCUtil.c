/* $Id:    */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2005                                               */
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

#include "ObitRPCUtil.h"
#include "ObitMem.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitRPCUtil.c
 * ObitRPC class utility function definitions.
 */

/*---------------Private function prototypes----------------*/
/** Private: .
static void  ObitRPCUtilNada  (void); */


/*----------------------Public functions---------------------------*/
/**
 * Send a File across an established XMLRPC connection
 * \param file       File to be sent
 * \param client     RPC client to use
 * \param serverURL  URL of service, e.g. "http://localhost:8765/RPC2"
 * \param err        Obit Error message
 */
void ObitRPCUtilFileSend (ObitFile *file, ObitRPC *client, 
			  gchar *serverURL, ObitErr *err)
{
  olong chunkSize=64*1024;  /* Size of pieces to send - there is a limit */
  olong numChunk, Chunk;
  ObitFilePos fileSize, filePos;
  ObitIOCode iretCode = OBIT_IO_SpecErr;
  ObitXML *xml, *reply=NULL;
  ObitInfoList *desc=NULL, *status=NULL;
  ObitInfoType infoType;
  gint32 dim[MAXINFOELEMDIM];
  gpointer buffer=NULL;
  olong retCode;
  gchar *reason=NULL, *fileName=NULL;
  gchar *routine = "ObitRPCCUtilFileSend";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitFileIsA(file));
  g_assert (ObitRPCIsA(client));

  /* Check that file exists */
  Obit_return_if_fail(((file->fileName!=NULL) && (ObitFileExist(file->fileName, err))), 
		      err, "%s: file not defined or doesn't exist: %s", 
		      routine, file->fileName );
  

  /* Open and Get file info */ 
  fileName = g_strdup(file->fileName);
  fileSize = ObitFileSize(fileName, err);
  iretCode = ObitFileOpen (file, fileName,OBIT_IO_ReadOnly, OBIT_IO_Binary, 
			  chunkSize, err);
  if (err->error) Obit_traceback_msg (err, routine, file->name);
  g_free(fileName);

  /* Get just basic file name without path */
  fileName = ObitFileName(file->fileName);

  /* Other info */
  chunkSize = file->blockSize;
  numChunk = (olong)(((ofloat)fileSize / (ofloat)chunkSize) + 0.99999);
  numChunk = MAX (1, numChunk);

  /* position at beginning */
  filePos = 0;
  Chunk   = 1;

  /* Descriptive info */
  desc = newObitInfoList();
  dim[0] = dim[1] = dim[2] = dim[3] = dim[4] = 1;
  /* File name */
  dim[0] = strlen (fileName);
  ObitInfoListPut(desc, "FileName", OBIT_string, dim, fileName, err);
  /* Chunk size */
  dim[0] = 1;
  ObitInfoListPut(desc, "size", OBIT_long, dim, &chunkSize, err);
  /* Number of chunks */
  dim[0] = 1;
  ObitInfoListPut(desc, "numChunk", OBIT_long, dim, &numChunk, err);
  if (err->error) Obit_traceback_msg (err, routine, file->name);

  /* Allocate buffer */
  buffer = ObitMemAllocName (chunkSize, "ObitRPCUtilIOBuffer");

  /* Loop reading */
  while (iretCode==OBIT_IO_OK) {

    /* read next chunk,  No more than rest of file */
    chunkSize = MIN (chunkSize, (fileSize-filePos));
    iretCode = ObitFileRead(file, filePos, chunkSize, buffer, err);
    if (err->error) {
      ObitMemFree(buffer);
      Obit_traceback_msg (err, routine, file->name);
    }
    filePos += chunkSize;

    /* Labeling for this chunk */
    dim[0] = 1;
    ObitInfoListPut(desc, "Chunk", OBIT_long, dim, &Chunk, err);
    ObitInfoListPut(desc, "size", OBIT_long, dim, &chunkSize, err);
    if (err->error) {
      ObitMemFree(buffer);
      Obit_traceback_msg (err, routine, file->name);
    }

    /* Package it into an XML parsel */
    xml = ObitXMLBlob2XML (buffer, desc, err);
    if (err->error) {
      ObitMemFree(buffer);
      Obit_traceback_msg (err, routine, file->name);
    }

    /* Give xml a method name */
    if (xml->func) g_free( xml->func);
    xml->func = g_strdup("copyFile");

    /* Shove it across the connection */
    reply = ObitRPCCall (client, serverURL, xml, &status, NULL, err);
    if (err->error) {
      ObitMemFree(buffer);
      xml = ObitXMLUnref(xml);
      Obit_traceback_msg (err, routine, file->name);
    }

    /* Check Status */
    retCode = -1;
    if (status) {
      ObitInfoListGet (status, "code", &infoType, dim, (gpointer)&retCode, err);
      ObitInfoListGetP (status, "reason", &infoType, dim, (gpointer*)&reason);
    }
    if (err->error) {
      ObitMemFree(buffer);
      xml = ObitXMLUnref(xml);
      Obit_traceback_msg (err, routine, file->name);
    }
    
    /* Did it work? */
    if (retCode!=0) {
      Obit_log_error(err, OBIT_InfoWarn, "%s: Could not talk to Server, code %d",
		     routine, retCode);
      Obit_log_error(err, OBIT_InfoWarn, "   because: %s",reason);
      ObitMemFree(buffer);
      xml = ObitXMLUnref(xml);
     return;
    }

    status  = ObitInfoListUnref(status);
    reply   = ObitXMLUnref(reply);
    xml     = ObitXMLUnref(xml);
    Chunk++;

    /* Are we there yet? */
    if (Chunk>numChunk) break;
  } /* end loop over file */

  g_free(fileName);
  ObitMemFree(buffer);      /* free buffer */
  ObitInfoListUnref(desc);   /* Free Info list */

  iretCode = ObitFileClose (file, err);
  if (err->error) Obit_traceback_msg (err, routine, file->name);
  
} /* end  ObitRPCUtilFileSend */

