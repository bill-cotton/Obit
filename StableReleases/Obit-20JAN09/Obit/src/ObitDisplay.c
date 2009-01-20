/* $Id$    */
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

#include "ObitAIPSDir.h"
#include "ObitFITS.h"
#include "ObitDisplay.h"
#include "ObitImageMosaic.h"
#include "ObitImageUtil.h"
#include "ObitRPCUtil.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitDisplay.c
 * ObitDisplay class function definitions.
 * This class is derived from the Obit base class.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitDisplay";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo structure ObitDisplayClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitDisplayClassInfo myClassInfo = {FALSE};

/*--------------- File Global Variables  ----------------*/


/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitDisplayInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitDisplayClear (gpointer in);

/** Private: ping server. */
static gboolean pingServer (ObitDisplay* display, ObitErr *err);

/** Private: Load Image. */
static gboolean LoadImage (ObitDisplay* display, ObitImage *image, 
			   gchar *fitsFile, 
			   olong field, olong nfield, ObitInfoList **request, 
			   ObitErr *err);

/** Private: Copy Image to FITS file on remote display. */
static ObitImage* CopyImage (ObitDisplay* display, ObitImage *image, 
			     gchar *filename, ObitErr *err);

/** Private: Delete temporary file */
static void DelTemp (char *filename, ObitErr *err);

/** Private: Edit Window. */
static gboolean EditWindow (ObitDisplay* display, ObitDConCleanWindow *window, 
			    olong field, ObitInfoList **request, ObitErr *err);

/**  Private: Turn Image info into an ObitXML */
static ObitXML* getImageInfo (ObitImage *image, gchar *fitsFile, 
			      olong field, olong nfield,
			      ObitErr *err);

/** Private: Set Class function pointers. */
static void ObitDisplayClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitDisplay* newObitDisplay (gchar* name)
{
  ObitDisplay* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitDisplayClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitDisplay));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitDisplayInit((gpointer)out);

 return out;
} /* end newObitDisplay */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitDisplayGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitDisplayClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitDisplayGetClass */

/**
 * Creates an ObitDisplay
 * \param name      An optional name for the object.
 * \param ServerURL URL of display server, 
 *        NULL defaults to "http://localhost:8765/RPC2"
 * \param err       Obit Error message
 * \return the new object.
 */
ObitDisplay* ObitDisplayCreate (gchar* name, gchar *ServerURL, ObitErr *err)
{
  ObitDisplay* out;
  gchar *routine = "ObitDisplayCreate";

  /* Only allow one */
  Obit_retval_if_fail((myClassInfo.numberDisplay<1), err, NULL,
		      "%s: ONLY One XMLDisplay  allowed", routine);
  myClassInfo.numberDisplay++;

  /* Create basic structure */
  out = newObitDisplay (name);
  out->client = ObitRPCCreateClient (name, err);
  if ((ServerURL!=NULL) && (strncmp(ServerURL,"      ",6)) && 
      (strncmp(ServerURL,"ObitView",8)))
    out->ServerURL = g_strdup(ServerURL);
  else if ((ServerURL!=NULL) && (!strncmp(ServerURL,"ObitView",8)))
    out->ServerURL = g_strdup("http://localhost:8765/RPC2");
  else out->ServerURL = g_strdup("http://localhost:8765/RPC2");
  
  return out;
} /* end ObitDisplayCreate */

/**
 * Display an image with optional CLEAN window editing
 * Either a single image or an image mosaic can be passed.
 * For a mosaic, the user can request other images from the mosaic.
 * If the display is remote, the image is copied as a gzipped FITS file.
 * \param display    ObitDisplay object
 * \param image      ObitImage or Image Mosaic
 * \param window     if nonNULL window corresponding to image
 *                   possibly edited by user. 
 *                   This MUST correspond to image.
 * \param field      If image= an ImageMosaic then this is the 
 *                   1-rel field number
 * \param err        Obit Error message
 * \return TRUE if user wants to quit
 */
gboolean ObitDisplayShow (ObitDisplay* display, Obit *image, 
			  ObitDConCleanWindow *window, 
			  olong field, ObitErr *err)
{
  gboolean out = FALSE;
  ObitImage *curImage=NULL, *tmpImage=NULL;
  ObitImageDesc *desc=NULL;
  ObitImageMosaic *mosaic=NULL;
  gboolean isMosaic, doEdit=FALSE;
  gint32 dim[MAXINFOELEMDIM];
  ObitInfoType infoType;
  ObitInfoList *request=NULL;
  olong req, ifield, jfield, oldError, nfield=1;
  gchar *tmpfilename = "ObitDisplayFITS.fits.gz";
  gchar *tmpFullpath=NULL, *tmpdir = "/tmp/";
  gchar *routine = "ObitDisplayShow";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return out;
  g_assert (ObitDisplayIsA(display));
  if (display->turnedOff) return out;  /* display turned off? */
  if (window) g_assert (ObitDConCleanWindowIsA(window));
  g_assert (ObitImageIsA(image) || ObitImageMosaicIsA(image));

  isMosaic = ObitImageMosaicIsA(image);
  if (isMosaic) {
    mosaic = (ObitImageMosaic*)image;
    nfield = mosaic->numberImages;

    /* If window given check for compatability */
    if (window) {
      doEdit = TRUE;
      Obit_retval_if_fail((mosaic->numberImages==window->nfield), err, out,
			  "%s: Mosaic and windows incompatible, No. field %d !=  %d",
			  routine, mosaic->numberImages, window->nfield);
    }
  } else {  /* Image */
    curImage = (ObitImage*)image;
    doEdit = window!=NULL;
  }
  
  /* Progress Report to show and clear pending messages */
  ObitErrLog(err);

 /* Check server availability (ping) */
  if (!pingServer(display, err)) {
    /* works but server temporarily unavailable */
    Obit_log_error(err, OBIT_InfoWarn, "%s: Display unavailable",
		   routine);
    return out;
  }
  
  /* Loop until no more requests */
  ifield = field;
  while (1) {
    
    /* check field number */
    Obit_retval_if_fail(((ifield>0) && (ifield<=nfield)), err, out,
			"%s: field %d out of range [1,%d]",
			routine, ifield, nfield);
    if (mosaic) curImage = mosaic->images[ifield-1];

    /* Load image - if display is not local make FITS image and pass */
    if (!strncmp (display->ServerURL, "http://localhost",  16)) {
      if (!LoadImage (display, curImage, NULL, ifield, nfield, &request, err)) {
	/* Failed */
	Obit_log_error(err, OBIT_InfoWarn, "%s: Image load failed",
		       routine);
	goto cleanup;
      }
    } else {
      /* Copy to remote */
      tmpImage = ObitImageUnref(tmpImage);
      tmpFullpath =  g_strconcat (tmpdir, tmpfilename, NULL);
      tmpImage = CopyImage (display, curImage, tmpFullpath, err);
      if (err->error) Obit_traceback_val (err, routine, display->name, out);

      /* Load image, it will have the same name on the remote end */
      if (!LoadImage (display, tmpImage, tmpfilename, ifield, nfield, 
		      &request, err)) {
	/* Failed */
	Obit_log_error(err, OBIT_InfoWarn, "%s: Image load failed",
		       routine);
	goto cleanup;
      }
      curImage = tmpImage;
    } /* end remote copy */

    
    /* Request? */
    req = OBIT_Request_Continue;
    if (request) {
      if (ObitInfoListGetTest (request, "Request", &infoType, dim, &req)) {
	/*fprintf (stdout, "Request from edit %d\n",req);*/
	if (ObitInfoListGetTest (request, "Field", &infoType, dim, &jfield)) {
	  /*fprintf (stdout, "Request field %d\n",field);*/
	}
      }
      request = ObitInfoListUnref(request);
    } /* End request handling */

    /* Edit windows if given */
    if (window && doEdit) {
      
      /* Check compatibility  
	 Descriptor on IO more reliable*/
      desc = (ObitImageDesc*)curImage->myIO->myDesc;
      Obit_retval_if_fail(((desc->inaxes[0]==window->naxis[ifield-1][0]) && 
			   (desc->inaxes[1]==window->naxis[ifield-1][1])), err, out,
			  "%s: Image and window incompatible [ %d, %d] [ %d, %d]",
			  routine, desc->inaxes[0], desc->inaxes[1], 
			  window->naxis[ifield-1][0], window->naxis[ifield-1][1]);
      
      /* User message */
      Obit_log_error(err, OBIT_InfoErr, "Edit boxes, continue, quit, abort?");
      ObitErrLog(err); /* show any  messages */

      /* Do edit */
      if (!EditWindow (display, window, ifield, &request, err)) {
	/* Failed */
	Obit_log_error(err, OBIT_InfoWarn, "%s: Image load failed",
		       routine);
	goto cleanup;
      }
      
      /* Request? */
      if (request) {
	if (ObitInfoListGetTest (request, "Request", &infoType, dim, &req)) {
	  /*fprintf (stdout, "Request from edit %d\n",req);*/
	  if (ObitInfoListGetTest (request, "Field", &infoType, dim, &jfield)) {
	    /*fprintf (stdout, "Request field %d\n",jfield);*/
	  }
	}
	request = ObitInfoListUnref(request);
      }
    } /* end editing windows */
    
    /* Handle any request */
    if (req==OBIT_Request_Abort) {
      Obit_log_error(err, OBIT_Error, "User requests Abort");
      goto cleanup;
    } else if (req==OBIT_Request_Quit) {
      Obit_log_error(err, OBIT_InfoWarn, "User requests Quit");
      out = TRUE;
      goto cleanup;
    } else if (req==OBIT_Request_Edit) {
      ifield = jfield;
      doEdit = window!=NULL;
    } else if (req==OBIT_Request_View) {
      ifield = jfield;
      doEdit = TRUE;
    } else if (req==OBIT_Request_NoTV) {
      ObitDisplayTurnOff (display);
      Obit_log_error(err, OBIT_InfoWarn, "User requests no further displays");
      break;
    } else if (req==OBIT_Request_Continue) {
      break;
    } else break;/* No request - exit loop */
    
  } /* end loop over requests */

  /* cleanup temporary files */
 cleanup:
  tmpImage = ObitImageUnref(tmpImage);
  oldError = err->error;
  err->error = 0;   /* trickery to get DelTemp to work on error */
  DelTemp(tmpFullpath, err);  /* Zap file - ImageZap can't */
  err->error = oldError;
  return out;
} /* end ObitDisplayShow */

/**
 * Enable display (default initial condition)
 * \param display    ObitDisplay object
 */
void ObitDisplayTurnOn (ObitDisplay* display)
{
  /* error check */
  g_assert (ObitDisplayIsA(display));

  display->turnedOff = FALSE;
} /* end ObitDisplayTurnOn */

/**
 * Disable display 
 * \param display    ObitDisplay object
 */
void ObitDisplayTurnOff (ObitDisplay* display)
{
  /* error check */
  g_assert (ObitDisplayIsA(display));

  display->turnedOff = TRUE;
} /* end ObitDisplayTurnOff */


/**
 * Initialize global ClassInfo Structure.
 */
void ObitDisplayClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitDisplayClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitDisplayClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitDisplayClassInfoDefFn (gpointer inClass)
{
  ObitDisplayClassInfo *theClass = (ObitDisplayClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitDisplayClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitDisplayClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitDisplayGetClass;
  theClass->newObit       = (newObitFP)newObitDisplay;
  theClass->ObitCopy      = NULL;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitDisplayClear;
  theClass->ObitInit      = (ObitInitFP)ObitDisplayInit;
  theClass->ObitDisplayCreate = 
    (ObitDisplayCreateFP)ObitDisplayCreate;
  theClass->ObitDisplayShow = 
    (ObitDisplayShowFP)ObitDisplayShow;
  theClass->ObitDisplayTurnOn = 
    (ObitDisplayTurnOnFP)ObitDisplayTurnOn;
  theClass->ObitDisplayTurnOff = 
    (ObitDisplayTurnOffFP)ObitDisplayTurnOff;
  theClass->numberDisplay  = 0;

} /* end ObitDisplayClassDefFn */


/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitDisplayInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitDisplay *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->thread    = newObitThread();
  in->turnedOff = FALSE;
  in->ServerURL = NULL;
} /* end ObitDisplayInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitDisplay* cast to an Obit*.
 */
void ObitDisplayClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitDisplay *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  in->thread  = ObitThreadUnref(in->thread);
  in->client  = ObitRPCUnref(in->client);
  if (in->ServerURL) g_free(in->ServerURL); in->ServerURL = NULL;
  myClassInfo.numberDisplay--;

  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitDisplayClear */

/**
 * Ping server
 * Returns TRUE is the ping was successful and the server is available.
 * A communitations failure will result in the display being "turned Off"
 * and err set
 * \param display   ObitDisplay object
 * \param err       Obit Error message
 * \return TRUE if available
 */
static gboolean pingServer (ObitDisplay* display, ObitErr *err)
{
  gboolean out = FALSE;
  ObitXML *xml=NULL, *result=NULL;
  ObitInfoList *status=NULL, *request=NULL;
  gint32 dim[MAXINFOELEMDIM];
  ObitInfoType infoType;
  olong retCode;
  gchar *reason=NULL;
  gchar *routine = "ObitDisplay:pingServer";

  /* existing error? */
  if (err->error) return out;

  xml = ObitXMLPing2XML(err);
  result = ObitRPCCall (display->client, display->ServerURL, xml, &status, &request, err);
  /* If something goes wrong with communications, turn off */
  if (err->error) {
    ObitDisplayTurnOff (display);
    ObitErrClearErr (err);  /* ignore error */
    goto cleanup;
  }

  /* Check Status */
  retCode = -1;
  if (status) {
    ObitInfoListGet (status, "code", &infoType, dim, (gpointer)&retCode, err);
    ObitInfoListGetP (status, "reason", &infoType, dim, (gpointer*)&reason);
  }
  if (err->error) {
    ObitDisplayTurnOff (display);
    goto cleanup;
  }

  /* Did it work? */
  if (retCode!=0) {
    Obit_log_error(err, OBIT_InfoWarn, "%s: Could not talk to Display, code %d",
		   routine, retCode);
    Obit_log_error(err, OBIT_InfoWarn, "   because: %s",reason);
    goto cleanup;
  }
  /* Must be OK */
  out = TRUE;

  /* Cleanup from load Image */
 cleanup:
  status  = ObitInfoListUnref(status);
  request = ObitInfoListUnref(request);
  xml     = ObitXMLUnref(xml);
  result  = ObitXMLUnref(result);
  if (err->error) Obit_traceback_val (err, routine, display->name, out);

  return out;
} /* end pingServer */

/**
 * Send request to display server to load an image
 * The image must be visible as described to the server.
 * May return a request for further client action. 
 * Returns TRUE is the ping was successful and the server is available.
 * A communitations failure will result in the display being "turned Off"
 * and err set
 * \param display   ObitDisplay object
 * \param image     Image to display (1st plane)
 * \param fitsFile  Name of FITS image to load as seen by Display server
 *                  ignored if NULL
 * \param field     Field number
 * \param nfield    Total number of images in mosaic (1 if single)
 * \param request   [out] request details if any
 * \param err       Obit Error message
 * \return TRUE if successful
 */
static gboolean LoadImage (ObitDisplay* display, ObitImage *image,  
			   gchar *fitsFile, 
			   olong field, olong nfield, ObitInfoList **request, 
			   ObitErr *err)
{
  gboolean out = FALSE;
  ObitXML *xml=NULL, *result=NULL;
  ObitInfoList *status=NULL;
  gint32 dim[MAXINFOELEMDIM];
  ObitInfoType infoType;
  olong retCode;
  gchar *reason=NULL;
  gchar *routine = "ObitDisplay:LoadImage";

  /* existing error? */
  if (err->error) return out;

  /* Get image info  */
  xml = getImageInfo (image, fitsFile, field, nfield, err);

  /* Load Image */
  result = ObitRPCCall (display->client, display->ServerURL, xml, &status, request, err);
  /* If something goes wrong with communications, turn off */
  if (err->error) {
    ObitDisplayTurnOff (display);
    ObitErrClearErr (err);  /* ignore error */
    goto cleanup;
  }

  /* Check Status */
  retCode = -1;
  if (status) {
    ObitInfoListGet (status, "code", &infoType, dim, &retCode, err);
    ObitInfoListGetP (status, "reason", &infoType, dim, (gpointer*)&reason);
  }
  if (err->error) {
    ObitDisplayTurnOff (display);
    goto cleanup;
  }

  /* Did it work? */
  if (retCode!=0) {
    Obit_log_error(err, OBIT_InfoWarn, "%s: Could not talk to Display, code %d",
		   routine, retCode);
    Obit_log_error(err, OBIT_InfoWarn, "   because: %s",reason);
    goto cleanup;
  }

  /* Cleanup from load Image */
 cleanup:
  status  = ObitInfoListUnref(status);
  xml     = ObitXMLUnref(xml);
  result  = ObitXMLUnref(result);
  if (err->error) Obit_traceback_val (err, routine, display->name, out);

  /* Must be OK */
  out = TRUE;

  return out;
} /* end LoadImage */

/**
 * Copy image to a gzipped FITS image on the remote display machine
 * Create local gzipped file name filename  quantized at 0.3 times the RMS.  
 * Copy to /tmp on remote host.
 * \param display   ObitDisplay object
 * \param image     Image to display (1st plane)
 * \param filename  Name of file to be copied
 * \param err       Obit Error message
 * \return True if successful
 */
static ObitImage* CopyImage (ObitDisplay* display, ObitImage *image,
			     gchar *filename, ObitErr *err)
{
  ObitImage *outFITS=NULL;
  ObitIOCode retCode;
  ObitFile *file=NULL;
  gint32 dim[MAXINFOELEMDIM] = {1, 1, 1, 1, 1};
  olong disk=0;
  olong chunk;
  ofloat factor;
  gchar *routine = "ObitDisplay:CopyImage";

  /* existing error? */
  if (err->error) return outFITS;

  /* Get rid of any old version */
  DelTemp (filename, err);
  if (err->error) Obit_traceback_val (err, routine, filename, outFITS);

  /* Quantize image to temporary file */
  dim[0] = dim[1] = dim[2] = 1;
  factor = 0.3;
  ObitInfoListAlwaysPut(image->info, "factor", OBIT_float, dim, &factor);
  outFITS = ObitImageUtilQuanFITS (image, filename, disk, err);
  if (err->error) Obit_traceback_val (err, routine, display->name, outFITS);

  /* Create file object */
  file =  newObitFile(filename);
  chunk = 4196;  /* size of chunk to copy */
  retCode = ObitFileOpen (file, filename, OBIT_IO_ReadOnly, OBIT_IO_Binary, 
			  chunk, err);
  if (err->error) Obit_traceback_val (err, routine, display->name, outFITS);
  retCode = ObitFileClose (file, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) 
    Obit_traceback_val (err, routine, filename, outFITS);

 /* Copy to remote */
  ObitRPCUtilFileSend (file, display->client, display->ServerURL, err);
  if (err->error) Obit_traceback_val (err, routine, display->name, outFITS);

  ObitFileUnref (file);  /* cleanup */

  return outFITS;
} /* end CopyImage  */

/** * Delete file named filename
 * \param filename  Name of file to be deleted 
 * \param err       Obit Error message
 */
static void DelTemp (char *filename, ObitErr *err)
{
  ObitFile *file=NULL;
  olong chunk;
  olong retCode;
  gchar *tfilename;
  gchar *routine = "ObitDisplay:DelTemp";

  /* existing error? */
  if (err->error) return;

  /* Anything to do? */
  if (filename==NULL) return;

  /* Drop any initial '!' */
  if (filename[0]=='!') tfilename = &filename[1];
  else tfilename = filename;

  /* Does it exist? */
  if (!ObitFileExist (filename, err)) return;

  file =  newObitFile(filename);
  chunk = 4196;  /* size of chunk */
  retCode = ObitFileOpen (file, tfilename, OBIT_IO_ReadOnly, OBIT_IO_Binary, 
			  chunk, err);
  if (err->error) Obit_traceback_msg (err, routine, filename);
  retCode = ObitFileClose (file, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) 
    Obit_traceback_msg (err, routine, filename);

  ObitFileZap (file, err);  /* Delete file */
  if (err->error)  Obit_traceback_msg (err, routine, filename);
} /* end DelTemp  */

/**
 * Send selected window to display server and allow user editing
 * Field window information is replaces in window with edited version.
 * May return a request for further client action. 
 * Returns TRUE is the edit was successful and the server is available.
 * A communitations failure will result in the display being "turned Off"
 * and err set
 * \param display   ObitDisplay object
 * \param image     Image to display (1st plane)
 * \param field     Field number (1-rel) to edit
 * \param request   [out] request details if any
 * \param err       Obit Error message
 * \return TRUE if successful
 */
static gboolean EditWindow (ObitDisplay* display, ObitDConCleanWindow *window, 
			   olong field, ObitInfoList **request, ObitErr *err)
{
  gboolean out = FALSE;
  ObitXML *xml=NULL, *result=NULL;
  ObitDConCleanWindow *newWindow=NULL;
  ObitInfoList *status=NULL;
  gint32 dim[MAXINFOELEMDIM] = {1, 1, 1, 1, 1};
  ObitInfoType infoType;
  olong retCode;
  gchar *reason=NULL;
  gchar *routine = "ObitDisplay:EditWindow";

  /* existing error? */
  if (err->error) return out;

  /* Convert window to ObitXML */
  xml = ObitXMLWindow2XML (window, (olong)field, err);

  /* Load Image */
  result = ObitRPCCall (display->client, display->ServerURL, xml, &status, request, err);
  /* If something goes wrong with communications, turn off */
  if (err->error) {
    ObitDisplayTurnOff (display);
    ObitErrClearErr (err);  /* ignore error */
    goto cleanup;
  }

  /* Check Status */
  retCode = -1;
  if (status) {
    ObitInfoListGet (status, "code", &infoType, dim, &retCode, err);
    ObitInfoListGetP (status, "reason", &infoType, dim, (gpointer*)&reason);
  }
  if (err->error) {
    ObitDisplayTurnOff (display);
    goto cleanup;
  }

  /* Did it work? */
  if (retCode!=0) {
    Obit_log_error(err, OBIT_InfoWarn, "%s: Could not talk to Display, code %d",
		   routine, retCode);
    Obit_log_error(err, OBIT_InfoWarn, "   because: %s",reason);
    goto cleanup;
  }

  /* Replace edited version in original window */
  newWindow = ObitXMLXML2Window (result, err);
  ObitDConCleanWindowReplaceField (newWindow, 1, window, field, err);
  if (err->error) goto cleanup;

  /* Cleanup from editing */
 cleanup:
  newWindow = ObitDConCleanWindowUnref(newWindow);
  status    = ObitInfoListUnref(status);
  xml       = ObitXMLUnref(xml);
  result    = ObitXMLUnref(result);
  if (err->error) Obit_traceback_val (err, routine, display->name, out);

  /* Must be OK */
  out = TRUE;

  return out;
} /* end EditWindow  */

/**
 * Send selected window to display server and allow user editing
 * \param image     Image to be described, ignored if fitsFile given
 * \param fitsFile  Name of FITS image to load as seen by Display server
 *                  ignored if NULL
 * \param field     Field number
 * \param nfield    Total number of images in mosaic (1 if single)
 * \param err       Obit Error message
 * \return ObitXML description of the file to send display server 
 */
static ObitXML* getImageInfo (ObitImage *image, gchar *fitsFile, 
			      olong field, olong nfield,
			      ObitErr *err)
{
  ObitXML *out = NULL;
  ObitIOType FileType;
  ObitInfoType type;
  gint32 dim[MAXINFOELEMDIM];
  olong i, disk, UserId=0, CNO, ASeq=0;
  gchar *name=NULL, *filenameP, filename[201], *AClass=NULL, *ADir=NULL;
  ObitAIPSDirCatEntry *catEntry=NULL;
  gchar *routine = "ObitDisplay:getImageInfo";

  /* existing error? */
  if (err->error) return out;


  /* Name as seen by Server given? */
  if (fitsFile) {
    disk = 0;
    name = g_strdup(fitsFile);
    FileType = OBIT_IO_FITS;
    /* Otherwise get  info by type from image */
  } else {
    
    /* Get FileType */
    if (!ObitInfoListGet(image->info, "FileType", &type, dim, 
			 (gpointer)&FileType, err)) {
      Obit_log_error(err, OBIT_Error, 
		     "%s: entry FileType not in InfoList Object %s",
		     routine, image->name);
      return out;
    }
    
    /* FITS */
    if (FileType==OBIT_IO_FITS) { 
      if(!ObitInfoListGet(image->info, "Disk", &type, dim, &disk, err))
	Obit_traceback_val (err, routine, image->name, out);
      
      if(!ObitInfoListGetP(image->info, "FileName", &type, dim, 
			   (gpointer*)&filenameP)) 
	Obit_traceback_val (err, routine, image->name, out);
      for (i=0; i<MIN (201, dim[0]); i++) filename[i] = filenameP[i];
      filename[i] = 0;
      
      /* Full name */
      name = ObitFITSFilename (disk, filename, err);
      if (err->error) goto cleanup;
      
      /* AIPS image */
    } else if (FileType==OBIT_IO_AIPS) { /* AIPS file */
      if(!ObitInfoListGet(image->info, "Disk", &type, dim, &disk, err)) 
	Obit_traceback_val (err, routine, image->name, out);
      
      if(!ObitInfoListGet(image->info, "User", &type, dim, &UserId, err)) 
	Obit_traceback_val (err, routine, image->name, out);
      
      if(!ObitInfoListGet(image->info, "CNO", &type, dim, &CNO, err))
	Obit_traceback_val (err, routine, image->name, out);
      
      ADir = ObitAIPSDirname(disk, err);
      if (err->error) Obit_traceback_val (err, routine, image->name, out);
      
      catEntry = ObitAIPSDirGetEntry(disk, UserId, CNO, err);
      if (err->error) Obit_traceback_val (err, routine, image->name, out);
      
      name = g_malloc(13*sizeof(gchar));
      for (i=0; i<12; i++) name[i] = catEntry->name[i]; name[i] = 0;
      AClass = g_malloc(7*sizeof(gchar));
      for (i=0; i<6; i++) AClass[i] = catEntry->class[i]; AClass[i] = 0;
      ASeq = catEntry->seq;
      
    } else if (FileType==OBIT_IO_MEM) {  /* Memory resident only */
      /* Can't cope with this one */
      Obit_log_error(err, OBIT_InfoErr, 
		     "%s: Can't send memory resident images to display",
		     routine);
      return out;
    }
  } /* end of get info from image */

  /* Construct XML */
  out = ObitXMLFileInfo2XML (FileType, name, AClass, ADir, ASeq, UserId, 
			     field, nfield, err);

  /* Cleanup */
 cleanup:
  if (name)     g_free(name);
  if (AClass)   g_free(AClass);
  if (catEntry) g_free(catEntry);
  if (err->error) Obit_traceback_val (err, routine, routine, out);

  return out;
} /* end getImageInfo */
