/* $Id: XMLRPCserver.c,v 1.15 2006/07/02 22:48:13 bcotton Exp $ */
/* XMLRPC server for ObitView */
/* Much of this material directly adapted from xmlrpc-c-1.2/examples */
/*-----------------------------------------------------------------------
*  Copyright (C) 2005-2008
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
#include <unistd.h> 
#include <errno.h>    /* global C error include */
#include <stdio.h>    /* i/o include */
#include <stdlib.h>   /* i/o include */
#include <signal.h>
#include <xmlrpc_client.h>
#include <xmlrpc_server.h>
#include <xmlrpc_server_abyss.h>
/*#include <xmlrpc-c/base.h>
  #include <xmlrpc-c/server.h>
  #include <xmlrpc-c/server_abyss.h>*/
#include "XMLRPCserver.h"
#include "ObitRPC.h"
#include "obitview.h"
#include "imagedisp.h"
#include "drawbox.h"
#include "Image2Pix.h"
#include "messagebox.h"
#include "xml.h"
#include "ObitDConCleanWindow.h"
#include "ObitMem.h"
#include "ObitFile.h"
#include "ObitInfoList.h"
#include <Xm/Xm.h> 

/*---------------Private function prototypes----------------*/
static void server_loop(gint port);
static xmlrpc_value*
ping(xmlrpc_env *  const envP, 
     xmlrpc_value * const paramArrayP,
     void *const userData);
static xmlrpc_value *
loadFITS(xmlrpc_env *   const envP, 
	 xmlrpc_value * const paramArrayP,
	 void *         const userData);

static xmlrpc_value *
loadImage(xmlrpc_env *   const envP, 
	   xmlrpc_value * const paramArrayP,
	   void *         const userData);

static xmlrpc_value *
editWindow(xmlrpc_env *   const envP, 
	   xmlrpc_value * const paramArrayP,
	   void *         const userData);

static xmlrpc_value *
copyFile(xmlrpc_env *   const envP, 
	 xmlrpc_value * const paramArrayP,
	 void *         const userData);


/*----------------- Globals ---------------------------*/
/** Structure for File information */
typedef struct {
  /** Image type, OBIT_IO_FITS or OBIT_IO_AIPS */
  ObitIOType Type;
  /** AIPS image sequence number  */
  olong ASeq;
  /** AIPS User number */
  olong AUser;
  /** FITS file path or AIPS Name (12 char) */
  gchar *Name;
  /** AIPS class (6 char), may be "None" for FITS */
  gchar *AClass;
  /** Path to AIPS directory, may be "None" for FITS */
  gchar *ADir;
  /** Field number */
  olong Field;
  /** Total Number of fields */
  olong NField;
} FileInfo;
FileInfo myFileInfo;  /* Static structure */

/* Are we busy doing something? */
gboolean WeAreBusy;

/*---------------Public functions ----------------*/
/**
 * start server thread
 * \param port  TCP port to watch
 */
void* start_myServer (void* port)
{

  /* initialization */
  WeAreBusy= FALSE;
  receive_flag = 0;
  return_flag  = 0;
  xmlrpcFunc   = XMLRPC_Inactive;
  if (pthread_mutex_init (&request_lock, NULL)) {
    fprintf (stderr, "Failed to initialize Request lock \n");
    return NULL;
  }

  /* Init LoadImage interface */
  myFileInfo.Name   = NULL;
  myFileInfo.AClass = NULL;
  myFileInfo.ADir   = NULL;
  myFileInfo.Type   = -1;
  myFileInfo.ASeq   = 0;
  myFileInfo.AUser  = 0;
  myFileInfo.Field  = 0;
  myFileInfo.NField = 0;

  /* Start looping */
  server_loop (*(olong*)port);

  return NULL;
} /* end start_server */

/**
 * start server and begin looping
 * \param port  TCP port to watch
 */
static void server_loop(gint port) 
{
  ObitRPC *server=NULL;
  ObitErr *err=NULL;

  /* Create Server */
  err = newObitErr();
  server = ObitRPCCreateServer("ObitView", err);

  ObitRPCAddMethod (server, "ping",       &ping,       NULL, err);
  ObitRPCAddMethod (server, "loadFITS",   &loadFITS,   NULL, err);
  ObitRPCAddMethod (server, "loadImage",  &loadImage,  NULL, err);
  ObitRPCAddMethod (server, "editWindow", &editWindow, NULL, err);
  ObitRPCAddMethod (server, "copyFile",   &copyFile,   NULL, err);

  /* Loop forever 'neath the streets of Boston */
  ObitRPCServerLoop(server, port, "/tmp/xmlrpc_log");

} /* end server_loop */

/**
 * Function to run in X thread to act on requests.
 * Communication with the XMLRPC thread is through globals.
 * (take out loop/timer in ObitView )
 * \li request_lock - mutex to protect global values
 * \li xmlrpcFunc   - code of function called
 * \li xmlrpcData   - data pointer
 * \li receive_flag - boolean flag indicating a request
 * \li return_flag  - boolean flag indicating a return value is ready
 */
void* XMLRPCWatcher (XtPointer clientData) 
{
  ImageDisplay  *IDdata;
  static gchar *IAm = "ObitView";
  gchar szErrMess[120];

  /*while(1) {*/

  IDdata = (ImageDisplay *)clientData;

  if (receive_flag) {
    pthread_mutex_lock(&request_lock); /* lock mutex */
    receive_flag=0;  /* Clear receive flag */
    requestList = ObitInfoListUnref(requestList);  /* make sure any old requests cleared */
    
    /* Branch by function */
    switch (xmlrpcFunc) {
    case XMLRPC_Inactive:   /* Called by mistake */
      break;
    case XMLRPC_Ping:   /* Ping no argument */
      xmlrpcData = IAm;
      return_flag=1;   /* Set return flag */
      break;
    case XMLRPC_LoadFits: /* Load FITS file, name passed */
      /*fprintf (stderr, "DEBUG Load FITS file %s\n", (char*)xmlrpcData);*/
      /* read FITS file to pixmap */
      FStrngFill (image[CurImag].FileName, (char*)xmlrpcData);
      image[CurImag].DataType = OBIT_IO_FITS;    /* FITS */
      image[CurImag].reLoad   = TRUE;            /* Force load */
      if (Image2Pix (&image[CurImag], IDdata, TRUE)) {
	  /* error */
	  sprintf (szErrMess, "Error reading FITS file = %s", (char*)xmlrpcData);
	  MessageShow (szErrMess);
	  xmlrpcData = "Load failed";  /* return message */
	} else 
	  xmlrpcData = "Loaded File";  /* return message */

      /* reset display */
      ResetDisplay(IDdata);
      return_flag=1;   /* Set return flag */
  
      break;
     case XMLRPC_LoadImage: /* Load Image */
       /*fprintf (stderr, "DEBUG LoadImage T %d S %d U %d name %s\n",*/
       /* myFileInfo.Type, myFileInfo.ASeq, */
       /* myFileInfo.AUser, myFileInfo.Name);*/
       /* if (myFileInfo.Type==OBIT_IO_AIPS) {*/
       /*  fprintf (stderr, "Class %s Dir %s\n",*/
       /*  myFileInfo.AClass, myFileInfo.ADir);*/
       /*}*/
       /* read file to pixmap */
       /* Get instructions */
       image[CurImag].DataType = myFileInfo.Type;
       image[CurImag].AIPSseq  = myFileInfo.ASeq;
       image[CurImag].AIPSuser = myFileInfo.AUser;
       image[CurImag].Field    = myFileInfo.Field;
       image[CurImag].NField   = myFileInfo.NField;
       FStrngFill (image[CurImag].FileName,  myFileInfo.Name);
       FStrngFill (image[CurImag].AIPSName,  myFileInfo.Name);
       FStrngFill (image[CurImag].AIPSClass, myFileInfo.AClass);
       FStrngFill (image[CurImag].AIPSDir,   myFileInfo.ADir);
       image[CurImag].reLoad   = TRUE;            /* Force load */
       if (Image2Pix (&image[CurImag], IDdata, TRUE)) {
	  /* error */
	  sprintf (szErrMess, "Error reading Image = %s", myFileInfo.Name);
	  MessageShow (szErrMess);
	  xmlrpcData = "Load failed";  /* return message */
	} else 
	  xmlrpcData = "Loaded File";  /* return message */

      /* reset display */
      ResetDisplay(IDdata);
      return_flag=1;   /* Set return flag */
  
      break;
   case XMLRPC_EditWindow:   /* Edit Window - returns the window */
      if (DrawBox (IDdata, (ObitDConCleanWindow*)xmlrpcData, 1)) {
	  /* error */
	  sprintf (szErrMess, "Error editing Window");
	  MessageShow (szErrMess);
	  xmlrpcData = "Edit failed";  /* return message */
	  return_flag=1;   /* Set return flag */
	} else 
	  return_flag=0;   /* clear return flag */
      break;
   case XMLRPC_FileCopy:   /* File Copy - nothing really to return */
     /* The X thread doesn't really have anything to do here */
     xmlrpcData = IAm;  /* In case anybody cares */
     return_flag=0;     /* clear return flag */
     break;
    default:
      fprintf (stderr, "Unknown function %d\n", xmlrpcFunc);
    }
    pthread_mutex_unlock(&request_lock); /* unlock mutex */
  }
  
  /* Let's not burn too many cycles in this event loop */
  usleep(250000);  /* 250 msec */
  /* } end of loop forever */
  return NULL;
} /* end XMLRPCWatcher */


/**
 * Set return ready flag
 * \param data  if nonNULL use as xmlrpcData
 */
void XMLRPCRelease (void *data)
{  
  /* Something going on? */
  if (xmlrpcFunc == XMLRPC_Inactive) return;
  pthread_mutex_lock(&request_lock);   /* lock mutex */
  return_flag = 1;  /* Send result if requested via XMLRPC */
  if (data) xmlrpcData = data;
  pthread_mutex_unlock(&request_lock); /* unlock mutex */
  xmlrpcFunc = XMLRPC_Inactive;

} /* end XMLRPCRelease */

/**
 * Set return ready flag
 * \param valueP  pointer to the item to be returned
 */
void XMLRPCSetReturn (void* valueP)
{  
  /* Something going on? */
  if (xmlrpcFunc == XMLRPC_Inactive) return;
  pthread_mutex_lock(&request_lock);   /* lock mutex */
  xmlrpcData = valueP;                 /* Set return value */
  pthread_mutex_unlock(&request_lock); /* unlock mutex */
  xmlrpcFunc = XMLRPC_Inactive;

} /* end XMLRPCSetReturn */

/*------------------- Private functions -------------------*/

/**
 * Handle ping function call 
 * \param envP          xmlrpc environment
 * \param paramArrayP   call argument as xml
 * \param userData      ignored
 * \return xml return value for function
 */
static xmlrpc_value *
ping(xmlrpc_env *   const envP, 
     xmlrpc_value * const paramArrayP,
     void *         const userData) 
{  
  xmlrpc_int x;

  /* Are we waiting on something? */
  if (WeAreBusy)
      return xmlrpc_build_value(envP, "{s:{s:i,s:s},s:s}", 
				"Status","code", (xmlrpc_int32)1,"reason","Busy",
				"Result", "Busy");
  WeAreBusy = TRUE; /* busy now */
  
  /* Parse our argument array. */
  xmlrpc_decompose_value(envP, paramArrayP, "(i)", &x);
  if (envP->fault_occurred)
    return NULL;

  /* Pass request to X thread */
  pthread_mutex_lock(&request_lock); /* lock mutex */
  receive_flag = 1;          /* now have a request */
  xmlrpcFunc = XMLRPC_Ping;  /* Function called */
  pthread_mutex_unlock(&request_lock); /* unlock mutex */

  /* Wait for results */
  while(1) {
    if (return_flag) {
      pthread_mutex_lock(&request_lock); /* lock mutex */
      return_flag = 0; /* clear flag */
      pthread_mutex_unlock(&request_lock); /* unlock mutex */
     
      /* Return our result. */
      WeAreBusy = FALSE; /* no longer busy */
      xmlrpcFunc = XMLRPC_Inactive;
      return xmlrpc_build_value(envP, "{s:{s:i,s:s},s:s}", 
				"Status","code", (xmlrpc_int32)0,"reason","OK",
				"Result", (char*)xmlrpcData);
    }
    /* Let's not burn too many cycles in this event loop */
    usleep(250000); /* 250 msec */
  }

} /* end ping */

/**
 * Handle loadFITS function call 
 * \param envP          xmlrpc environment
 * \param paramArrayP   call argument as xml
 * \param userData      ignored
 * \return xml return value for function
 */
static xmlrpc_value *
loadFITS(xmlrpc_env *   const envP, 
	 xmlrpc_value * const paramArrayP,
	 void *         const userData)
{  
  gchar *filename;
  
  /* Are we waiting on something? */
  if (WeAreBusy)
      return xmlrpc_build_value(envP, "{s:{s:i,s:s},s:s}", 
				"Status","code", (xmlrpc_int32)1,"reason","Busy",
 				"Result", "Busy");
   WeAreBusy = TRUE; /* busy now */
 
  /*fprintf (stderr, "DEBUG in loadFITS\n");*/

  /* Parse our argument array. */
  xmlrpc_decompose_value(envP, paramArrayP, "(s)", &filename);
  if (envP->fault_occurred)
    return NULL;

  /* Pass request to X thread */
  pthread_mutex_lock(&request_lock); /* lock mutex */
  receive_flag = 1;              /* now have a request */
  xmlrpcFunc = XMLRPC_LoadFits;  /* Function called */
  xmlrpcData = filename;         /* pass file name */
  pthread_mutex_unlock(&request_lock); /* unlock mutex */

  /* Wait for results */
  while(1) {
    if (return_flag) {
      pthread_mutex_lock(&request_lock); /* lock mutex */
      return_flag = 0; /* clear flag */
      pthread_mutex_unlock(&request_lock); /* unlock mutex */
      g_free (filename);   /* Cleanup */
     
      /* Return our result. */
      WeAreBusy = FALSE; /* no longer busy */
      xmlrpcFunc = XMLRPC_Inactive;
      /* It work? */
      if (!image[CurImag].valid)
	return xmlrpc_build_value(envP, "{s:{s:i,s:s},s:s}", 
				  "Status","code", (xmlrpc_int32)1,"reason","Failed",
				  "Result", (char*)xmlrpcData);
      else /*return xmlrpc_build_value(envP, "s", "Boo");*/
	  return xmlrpc_build_value(envP, "{s:{s:i,s:s},s:s}", 
				  "Status","code", (xmlrpc_int32)0,"reason","OK",
				  "Result", (char*)xmlrpcData);
    }
    /* Let's not burn too many cycles in this event loop */
    usleep(250000); /* 250 msec */
  }

} /* end loadFITS */

/**
 * Handle loadImage function call 
 * \param envP          xmlrpc environment
 * \param paramArrayP   call argument as xml
 * \param userData      ignored
 * \return xml return value for function
 */
static xmlrpc_value *
loadImage(xmlrpc_env *   const envP, 
	 xmlrpc_value * const paramArrayP,
	 void *         const userData)
{  
  xmlrpc_value *strt, *req;
  ObitXML *arg;
  
  /* Are we waiting on something? */
  if (WeAreBusy)
      return xmlrpc_build_value(envP, "{s:{s:i,s:s},s:s}", 
				"Status","code", (xmlrpc_int32)1,"reason","Busy",
				"Result", "Busy");
  WeAreBusy = TRUE; /* busy now */
 
  /*fprintf (stderr, "DEBUG in loadImage\n");*/

  /* Parse our argument structure. */
  xmlrpc_decompose_value(envP, paramArrayP, "(S)", &strt);
  arg =  ObitXMLReturn ("loadImage", strt, err);  /* Convert to ObitXML */
  /* Cleanup old versions */
  if (myFileInfo.Name)   g_free(myFileInfo.Name);   myFileInfo.Name=NULL;
  if (myFileInfo.AClass) g_free(myFileInfo.AClass); myFileInfo.AClass=NULL;
  if (myFileInfo.ADir)   g_free(myFileInfo.ADir);   myFileInfo.ADir=NULL;
  ObitXMLXML2FileInfo(arg, &myFileInfo.Type, &myFileInfo.Name, 
		      &myFileInfo.AClass, &myFileInfo.ADir, &myFileInfo.ASeq, 
		      &myFileInfo.AUser,&myFileInfo.Field,&myFileInfo.NField,err);
  ObitXMLUnref(arg);
  /* DEBUG 
     fprintf (stderr,"Name %s Class %s\n", myFileInfo.Name, myFileInfo.AClass);*/
  if (envP->fault_occurred) {
    ObitErrLog (err);
    WeAreBusy = FALSE; /* no longer busy */
    xmlrpcFunc = XMLRPC_Inactive;
    return xmlrpc_build_value(envP, "{s:{s:i,s:s},s:s}", 
			      "Status","code", 
			      (xmlrpc_int32)500,"reason","Invalid argument",
			      "Result", (char*)xmlrpcData);

  }

  /* Pass request to X thread */
  pthread_mutex_lock(&request_lock);  /* lock mutex */
  receive_flag = 1;                   /* now have a request */
  xmlrpcFunc = XMLRPC_LoadImage;      /* Function called */
  xmlrpcData = NULL;                  /* arguments passed through global myFileInfo */
  pthread_mutex_unlock(&request_lock);/* unlock mutex */

  /* Wait for results */
  while(1) {
    if (return_flag) {
      pthread_mutex_lock(&request_lock); /* lock mutex */
      return_flag = 0; /* clear flag */
      pthread_mutex_unlock(&request_lock); /* unlock mutex */
     
      /* Return our result.- string set in XMLRPCWatcher */
      WeAreBusy = FALSE; /* no longer busy */
      /* It work? */
      if (image[CurImag].valid) { /* OK */
	/* Make and return any request object */
	if (requestList) {
	  /* Make ObitXML but delete all but xmlrpc_value part */
	  arg = ObitXMLInfoList2XML (requestList, err);
	  req = arg->parmP;
	  xmlrpc_INCREF(req);
	  ObitXMLUnref(arg);
	  requestList = ObitInfoListUnref(requestList);
	  return xmlrpc_build_value(envP, "{s:{s:i,s:s},s:V,s:V}", 
				    "Status","code", (xmlrpc_int32)0,"reason","OK",
				    "Result", (char*)xmlrpcData,
				    "Request", req);
	} else  /* OK but no request */
	  return xmlrpc_build_value(envP, "{s:{s:i,s:s},s:s}", 
				    "Status","code", (xmlrpc_int32)0,"reason","OK",
				    "Result", (char*)xmlrpcData);
      }
      else /* FAILED */  /* return xmlrpc_build_value(envP, "s", "Boo");*/
	return xmlrpc_build_value(envP, "{s:{s:i,s:s},s:s}", 
				  "Status","code", (xmlrpc_int32)1,"reason","Failed",
				  "Result", (char*)xmlrpcData);
    }
    /* Let's not burn too many cycles in this event loop */
    usleep(250000); /* 250 msec */
  }

} /* end loadImage */

/**
 * Handle editWindow function call 
 * \param envP          xmlrpc environment
 * \param paramArrayP   call argument as xml
 * \param userData      ignored
 * \return xml return value for function
 */
static xmlrpc_value *
editWindow(xmlrpc_env *   const envP, 
	   xmlrpc_value * const paramArrayP,
	   void *         const userData) 
{  
  ObitDConCleanWindow *window;
  xmlrpc_value *strt, *req;
  ObitXML *arg;

   /* Are we waiting on something? */
  if (WeAreBusy)
    return xmlrpc_build_value(envP, "{s:{s:i,s:s},s:s}", 
			      "Status","code", (xmlrpc_int32)1,"reason","Busy",
			      "Result", "Busy");
   WeAreBusy = TRUE; /* busy now */
 
 /*fprintf (stderr, "DEBUG in editWindow\n");*/

  /* Data comes packed into an array */
  xmlrpc_decompose_value(envP, paramArrayP, "(S)", &strt);
  arg =  ObitXMLReturn ("editWindow", strt, err);  /* Convert to ObitXML */
  window = ObitXMLXML2Window (arg, err);
  ObitXMLUnref(arg);
  if (envP->fault_occurred) {
    ObitErrLog (err);
    WeAreBusy = FALSE; /* no longer busy */
    xmlrpcFunc = XMLRPC_Inactive;
    return xmlrpc_build_value(envP, "{s:{s:i,s:s},s:s}", 
			      "Status","code", 
			      (xmlrpc_int32)500,"reason","Invalid argument",
			      "Result", (char*)xmlrpcData);
  }
 
  /* Pass request to X thread */
  pthread_mutex_lock(&request_lock);   /* lock mutex */
  receive_flag = 1;                    /* now have a request */
  xmlrpcFunc = XMLRPC_EditWindow;      /* Function called */
  xmlrpcData = window;                 /* pass window */
  pthread_mutex_unlock(&request_lock); /* unlock mutex */

  /* Wait for results */
  while(1) {
    if (return_flag) {
      pthread_mutex_lock(&request_lock);   /* lock mutex */
      return_flag = 0; /* clear flag */
      pthread_mutex_unlock(&request_lock); /* unlock mutex */

      /* Return our result. */
      /* Did it work? */
      if (strncmp (xmlrpcData, "Edit failed", 11)) {
	/* Make ObitXML but delete all but xmlrpc_value part */
	arg = ObitXMLWindow2XML (window, 1, err);
	strt = arg->parmP;
	xmlrpc_INCREF(strt);
	ObitXMLUnref(arg);

	/* Make and return any request object */
	if (requestList) {
	  /* Make ObitXML but delete all but xmlrpc_value part */
	  arg = ObitXMLInfoList2XML (requestList, err);
	  req = arg->parmP;
	  xmlrpc_INCREF(req);
	  ObitXMLUnref(arg);
	  requestList = ObitInfoListUnref(requestList);
	  xmlrpcFunc = XMLRPC_Inactive;
	  WeAreBusy = FALSE; /* no longer busy */
	  return xmlrpc_build_value(envP, "{s:{s:i,s:s},s:V,s:V}", 
				    "Status","code", (xmlrpc_int32)0,"reason","OK",
				    "Result", strt,
				    "Request", req);
	} else {
	  xmlrpcFunc = XMLRPC_Inactive;
	  WeAreBusy = FALSE; /* no longer busy */
	  return xmlrpc_build_value(envP, "{s:{s:i,s:s},s:V}", 
				    "Status","code", (xmlrpc_int32)0,"reason","OK",
				    "Result", strt);
	}
      } else { /* Failed */
	xmlrpcFunc = XMLRPC_Inactive;
	WeAreBusy = FALSE; /* no longer busy */
	return xmlrpc_build_value(envP, "{s:{s:i,s:s},s:s}", 
				"Status","code", (xmlrpc_int32)600,"reason","Failed",
				"Result", (char*)xmlrpcData);
      }
    }
    /* Let's not burn too many cycles in this event loop */
    usleep(250000); /* 250 msec */
  }

} /* end editWindow */

/**
 * Handle copyFile function call 
 * This routine does all of the work and sends a return directly.
 * \param envP          xmlrpc environment
 * \param paramArrayP   call argument as xml
 * \param userData      ignored
 * \return xml return value for function
 */
static xmlrpc_value *
copyFile(xmlrpc_env *   const envP, 
	   xmlrpc_value * const paramArrayP,
	   void *         const userData) 
{  
  ObitFile *file=NULL;
  olong numChunk, Chunk, chunkSize;
  ObitIOCode retCode;
  xmlrpc_value *strt;
  gpointer fileData;
  gchar *fileName=NULL;
  ObitInfoList *desc;
  ObitInfoType infoType;
  gint32 dim[MAXINFOELEMDIM];
  ObitXML *arg;

   /* Are we waiting on something? Need to see if this is a continuation of last ?*/
  if (WeAreBusy)
    return xmlrpc_build_value(envP, "{s:{s:i,s:s},s:s}", 
			      "Status","code", (xmlrpc_int32)1,"reason","Busy",
			      "Result", "Busy");
   WeAreBusy = TRUE; /* busy now */
 
   /*fprintf (stderr, "DEBUG in copyFile\n");*/

  /* Data comes packed into an array */
  xmlrpc_decompose_value(envP, paramArrayP, "(S)", &strt);
  arg =  ObitXMLReturn ("copyFile", strt, err); /* Convert to ObitXML */
  fileData = ObitXMLXML2Blob (arg, &desc, err);
  ObitXMLUnref(arg);
  if (envP->fault_occurred || err->error) {
    ObitErrLog (err);
    xmlrpcFunc = XMLRPC_Inactive;
    WeAreBusy = FALSE; /* no longer busy */
    return xmlrpc_build_value(envP, "{s:{s:i,s:s},s:s}", 
			      "Status","code", 
			      (xmlrpc_int32)500,"reason","Invalid argument",
			      "Result", "Failed");
  }

  /* Get descriptive data */
  ObitInfoListGetP(desc, "FileName", &infoType, dim, (gpointer*)&fileName);
  ObitInfoListGet(desc, "size",      &infoType, dim, &chunkSize, err);
  ObitInfoListGet(desc, "numChunk",  &infoType, dim, &numChunk,  err);
  ObitInfoListGet(desc, "Chunk",     &infoType, dim, &Chunk,     err);

  /* Create File Object */
  file = newObitFile("FileCopy");

  /* If this is the first chunk and it exists, delete it first */
  if (Chunk==1) {
    file->fileName = g_strdup(fileName);
    file = ObitFileZap (file, err);
    file = newObitFile("FileCopy");
  }

  /* Open it */
  retCode = ObitFileOpen (file, fileName, OBIT_IO_ReadWrite, OBIT_IO_Binary, 
			  chunkSize, err);
  /* Position at end if not first */
  if (Chunk>1) retCode = ObitFileEnd (file, err);
  /* Write */
  retCode = ObitFileWrite (file, -1, chunkSize, fileData,  err);
  retCode = ObitFileClose (file, err);

  /* Cleanup */
  ObitMemFree(fileData);
  ObitFileUnref(file);
  ObitInfoListUnref(desc);

  /*  Did an error occur? */
  if (err->error) {
    ObitErrLog (err);
    xmlrpcFunc = XMLRPC_Inactive;
    WeAreBusy = FALSE; /* no longer busy */
    return xmlrpc_build_value(envP, "{s:{s:i,s:s},s:s}", 
			      "Status","code", 
			      (xmlrpc_int32)500,"reason","Error writing file",
			      "Result", "Write Failed");
  }
  
  
  /* It all seems to have worked */
  xmlrpcFunc = XMLRPC_Inactive;
  WeAreBusy = FALSE; /* no longer busy */
  return xmlrpc_build_value(envP, "{s:{s:i,s:s},s:s}", 
			    "Status","code", (xmlrpc_int32)0,"reason","OK",
			    "Result", "OK");

} /* end copyFile */

