/* $Id$        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2005-2016                                          */
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
#ifndef OBITDISPLAY_H 
#define OBITDISPLAY_H 

#include "Obit.h"
#include "ObitErr.h"
#include "ObitRPC.h"
#include "ObitDConCleanWindow.h"

/*-------- Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitDisplay.h
 * ObitDisplay Image Display class
 *
 * This class is derived from the #Obit class.
 *
 * This class communicates with the Image Display server
 * The implementation is based on ObitView, only one allowed
 * and return
 *
 * \section ObitDisplay Display Interface
 * The ObitDisplay class is the client interface between applications 
 * software and the image display server which runs asynchronously, 
 * perhaps on another computer.  The server is specified by a URL in
 * the #ObitDisplayCreate call and an example of a local server running on 
 * port 8765 is "http://localhost:8765/RPC2"
 * The actual communications uses the  #ObitRPC class which communicates 
 * using internet (http) protocols and xml to package the information 
 * communicated.  
 * The actual communication is through #ObitRPCCall.
 * In Obit, xml is encapsulated in the #ObitXML class.  
 * In ObitView the server is managed in XMLRPCserver.
 *
 *    The communication model is the client/server model where the client makes 
 * a call with a single (xml) argument and the server performs its service and 
 * returns a single (xml) reply. (Since xml is very flexible the requirement of 
 * a single argument and reply is not a limitation).  
 * In this model, the interaction between client and server is stateless, i.e.,
 * makes no assumption about any previous or future interactions, although both
 * the client and do have state.  
 * In this implementation, the response from the server may include a request for 
 * a further action by the client.  Examples of this are to display and/or edit 
 * the CLEAN window for another field of a mosaic or to abort the program.
 *
 * The returned xml from the server contains up to three components:
 * \li Result
 * This will be the #ObitXML returned by the #ObitRPCCall and depends on the
 * Call type.  There is an #ObitXML function to translate this to the 
 * call-specific useful form.
 * \li Status
 * This is returned as an #ObitInfoList by #ObitRPCCall and gives the status 
 * (success or failure) of the call.  See below for details.
 * \li Request [optional]
 * The request section, if present, is returned as an #ObitInfoList by #ObitRPCCall .
 * This describes a request by the server (generally user input) for an action by
 * the client and is described further below.
 *
 * \subsection Status Returned
 * The status argument returned by #ObitRPCCall is an #ObitInfoList translation of 
 * Status portion of the server response.
 * There will be two entries:
 * \li code (olong) which is a return code, 0=OK everything else indicates an error
 * \li reason (gchar*) which is a reason for the code.  
 * "OK" is given for success, "Busy" if the server is still processing another 
 * request (code=1).  Other values may also appear.
 *
 * \subsection Request Returned
 * The request argument returned by #ObitRPCCall is an #ObitInfoList translation of 
 * Request portion of the server response.  
 * This is a copy of an ObitInfoList generated by the server and the details depend 
 * on the exact request.  The Request code is the "Request" entry (olong) in the 
 * ObitInfoList(#ObitDisplayRequest enum defined for convienence) and the  
 * following values are defined:
 * \li 0 (OBIT_Request_Continue) Continue program (no request for action) , 
 * \li 1 (OBIT_Request_Abort)    Abort program, current results are assumed of no value
 * \li 2 (OBIT_Request_Quit)     Quit, graceful shutdown with current results saved
 * \li 3 (OBIT_Request_Edit)     Edit, send window for editing, must match current image
 * \li 4 (OBIT_Request_View)     View, send another field in the same ObitImageMosaic 
 *     and possibly corresponding window to edit.  
 *     This request will also include member "Field" (olong) which is the 1-rel 
 *     field number requested.
 *
 * \section Server Functions
 *    The display server currently supported is ObitView which has the following
 * callable functions:
 * \li ping
 * Simple function to test communications.
 * \li loadFITS (not used here)
 * Sends the name of a FITS image to display.
 * Name should be the full path as seen from the server and may be an internet URL.
 * \li loadImage
 * Load any image type accessable by Obit in the server
 * \li editWindow
 * Passed an ObitDConCleanWindow for a single image, the server overlays it on
 * the image and allows interactive editing, returning the edited window structure.
 *
 * \subsection ping call arguments and return
 * The argument passed to the ping call is generated by #ObitXMLPing2XML 
 * although no actual information is passed to the server (a random integer is used).
 * The return value is a string containing the name of the server (e.g. "ObitView").
 * The Status return value from #ObitRPCCall gives the availability of the server 
 * (may be "busy"). If there was a communications failure (e.g. no server or network
 * connection or bad xml) then an error will be entered on err.
 *
 * \subsection loadFITS call arguments and return
 * This call is not used by Obit.  
 * The argument is a single string giving the full path or URL of a FITS image to load.
 *
 * \subsection loadImage call arguments and return
 * This call is used to display an image accessable by Obit running in the server.
 * The argument is generated by #ObitXMLFileInfo2XML which takes a description of the
 * file as seen by the server. 
 * Note: the FITS path may be an internet URL but only local AIPS images are
 * accessable.  The return value is a string of either "Loaded File" or "Load failed"
 * but the succes of the function is better determined from the Status return.
 *
 * \subsection editWindow call arguments and return
 * This call sends the #ObitDConCleanWindow (one field only) corresponding to the
 * previous loadImage call and allows the user to interactively edit the window.
 * The call argument to #ObitRPCCall is generated by #ObitXMLWindow2XML and the 
 * return value converted into a (single field) #ObitDConCleanWindow  by
 * #ObitXMLXML2Window.  
 * The return from this call may include a Request for a further action.
 *
 *
 * \subsection markPos call arguments and return
 * This sends a celestial position as ("hh mm ss.s dd mm ss.s") to be marked on the display
 * The call argument to #ObitRPCCall is generated by #ObitXMLMarkPosXML.
 *
 * \section ObitDisplayaccess Creators and Destructors
 * An ObitDisplay will usually be created using ObitDisplayCreate which allows 
 * specifying a name for the object as well as other information.
 *
 * A copy of a pointer to an ObitDisplay should always be made using the
 * #ObitDisplayRef function which updates the reference count in the object.
 * Then whenever freeing an ObitDisplay or changing a pointer, the function
 * #ObitDisplayUnref will decrement the reference count and destroy the object
 * when the reference count hits 0.
 * There is no explicit destructor.
 */

/*-------------- enumerations -------------------------------------*/
/**
 * \enum obitDisplayRequest
 * enum Display request coads (MUST be synchronized with server usage
 * which is defined in ObitRPC.h #ObitRPCRequestType)
 * This specifies the request
 */
enum obitDisplayRequest {
  /** Continue program (no request for action)  */
  OBIT_Request_Continue = 0, 
  /** Abort  program current results are assumed of no value */
  OBIT_Request_Abort, 
  /** Quit, graceful shutdown with current results saved */
  OBIT_Request_Quit, 
  /** No more TV display */
  OBIT_Request_NoTV,
  /** View, send another field in the same ObitImageMosaic */
  OBIT_Request_View,
  /** Edit, send window for editing */
  OBIT_Request_Edit
}; /* end enum obitDisplayRequest */
/** typedef for enum for ObitDisplayRequest object status. */
typedef enum obitDisplayRequest ObitDisplayRequest;

/*--------------Class definitions-------------------------------------*/
/** ObitDisplay Class structure. */
typedef struct {
#include "ObitDisplayDef.h"   /* this class definition */
} ObitDisplay;

/*----------------- Macroes ---------------------------*/
/** 
 * Macro to unreference (and possibly destroy) an ObitDisplay
 * returns a ObitDisplay*.
 * in = object to unreference
 */
#define ObitDisplayUnref(in) ObitUnref (in)

/** 
 * Macro to reference (update reference count) an ObitDisplay.
 * returns a ObitDisplay*.
 * in = object to reference
 */
#define ObitDisplayRef(in) ObitRef (in)

/** 
 * Macro to determine if an object is the member of this or a 
 * derived class.
 * Returns TRUE if a member, else FALSE
 * in = object to reference
 */
#define ObitDisplayIsA(in) ObitIsA (in, ObitDisplayGetClass())

/*---------------Public functions---------------------------*/
/** Public: Class initializer. */
void ObitDisplayClassInit (void);

/** Public: Default Constructor. */
ObitDisplay* newObitDisplay (gchar* name);

/** Public: Create/initialize ObitDisplay structures */
ObitDisplay* ObitDisplayCreate (gchar* name, gchar *ServerURL, ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef ObitDisplay* (*ObitDisplayCreateFP) (gchar* name, gchar *ServerURL, 
					     ObitErr *err);

/** Public: ClassInfo pointer */
gconstpointer ObitDisplayGetClass (void);

/** Public: Send Display and Window edit request */
gboolean ObitDisplayShow (ObitDisplay* display, Obit *image, 
			  ObitDConCleanWindow *window, 
			  olong field, ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef gboolean (*ObitDisplayShowFP) (ObitDisplay* display, Obit *image, 
				       ObitDConCleanWindow *window, olong field, 
				       ObitErr *err);
/** Public: Mark position on display */
gboolean ObitDisplayMarkPos (ObitDisplay* display, gchar *pos, ObitErr *err);
/** Typedef for definition of class pointer structure */
typedef gboolean (*ObitDisplayMarkPosFP) (ObitDisplay* display, gchar *pos, ObitErr *err);

/** Public: Turn display on */
void ObitDisplayTurnOn (ObitDisplay* display);
/** Typedef for definition of class pointer structure */
typedef void (*ObitDisplayTurnOnFP) (ObitDisplay* display);

/** Public: Turn display off */
void ObitDisplayTurnOff (ObitDisplay* display);
/** Typedef for definition of class pointer structure */
typedef void (*ObitDisplayTurnOffFP) (ObitDisplay* display);


/*----------- ClassInfo Structure -----------------------------------*/
/**
 * ClassInfo Structure.
 * Contains class name, a pointer to any parent class
 * (NULL if none) and function pointers.
 */
typedef struct  {
#include "ObitDisplayClassDef.h"
} ObitDisplayClassInfo; 

#endif /* OBITDISPLAY_H */ 
