/* $Id$  */
/* help dialog box  for ObitView */
/* adopted from "Power programming Motif" by E. F. Johnson and
   K. Reichard, 1993, MIS Press, New York */
/*-----------------------------------------------------------------------
*  Copyright (C) 1998-2016
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
*-----------------------------------------------------------------------*/

#include <Xm/Xm.h> 
#include <Xm/Form.h> 
#include <Xm/Label.h> 
#include <Xm/List.h> 
#include <Xm/PushB.h> 
#include <Xm/RowColumn.h> 
#include <Xm/Separator.h> 
#include <Xm/Text.h> 
#include <Xm/Text.h> 
#include <Xm/DialogS.h> 
#include <X11/IntrinsicP.h>
#include <glib.h>
#include <stdio.h>
#include "helpbox.h"
#include <ObitVersion.h>

/**
 *  \file helpbox.c
 * displays "help" dialog.
 */

/*--------------- file global data ----------------*/
Widget  help_dialog = (Widget)NULL;
Widget  help_topic  = (Widget)NULL;
Widget  help_text   = (Widget)NULL;

/* Globals for storing help information max 100 */
Boolean HelpBoxActive = False;
int number_topics;      /* how many topics are there data for */
char* topic_title[100]; /* name of item */
char** topic_text[100]; /* help text, each is an array of strings */

/*---------------Private function prototypes----------------*/
void InitHelpText(void);
void HelpBoxSelectCB(Widget parent, XtPointer clientData, XtPointer callData);

/**
 * Callback to unmanage dialog
 * \param w           widget activated
 * \param clientData  client data
 * \param callData    call data
 */
static void unmanage_helpdialogCB(Widget widget, 
				  XtPointer client_data, 
				  XtPointer call_data) {
  Widget dialog = (Widget)client_data;
  
  if (dialog != (Widget)NULL) {
    if (XtIsManaged(dialog)) XtUnmanageChild(dialog);
  }
} /* end unmanage_helpdialogCB */

/**
 * Show Help dialog
 */
void HelpBoxShow(void) {
  if (!XtIsManaged(help_dialog)) {
    /* clear any previous text, add instructions */
    XmTextReplace(help_text,(XmTextPosition)0, 
		  XmTextGetLastPosition(help_text), 
		  "     Click topic on left to display help.");
    XtManageChild(help_dialog);
  }
} /* end HelpBoxShow */

/**
 * Create help box 
 * \param parent           parent widget
 * \param topic_callback   topic callback function
 * \param topic_data       topic callback data
 * \param help_callback    callback for help
 * \param help_text_data   help text
 * \param help_topic_data  help topic
 * \return true if created, else previously existed 
 */
Boolean HelpBoxCreate (Widget parent, 
		       XtCallbackProc topic_callback,
		       XtPointer topic_data,
		       XtCallbackProc help_callback,
		       XtPointer help_text_data,
		       XtPointer help_topic_data) {
  Arg         args[20];
  Cardinal    n;
  Widget      row, sep, dismiss;
  XFontStruct *XFont;
  GC          gc;
  int         chei, cwid;
  Dimension   topicWid, textHei, textWid;
  
  /* If there's a Help box lurking around unseen wake it up instead */
  if (HelpBoxActive) {
    HelpBoxShow();
    return False;
  }
  HelpBoxActive = True; /* soon to be true */
  
  /* see how big characters are to make boxes to size */
  /* create graphics context for box */
  gc = XCreateGC (XtDisplay (parent), 
		  DefaultRootWindow (XtDisplay(parent)),
		  0, NULL);
  
  /* how big are the characters */
  XFont = XQueryFont(XtDisplay (parent), XGContextFromGC(gc));
  chei = XFont->ascent + XFont->descent + 2;
  cwid = XFont->max_bounds.width;
  XFreeFontInfo(NULL, XFont, 0);
  
  /* text ~80 char x 20 lines, topics ~25 char wide */
  topicWid = 35*cwid;
  textWid = 80*cwid;
  textHei = 20*chei;
  
  if (gc) XFreeGC(XtDisplay(parent), gc); gc = NULL;
  
  n = 0;
  XtSetArg(args[n], XmNallowResize, True); n++;
  XtSetArg(args[n], XmNtitle, "FITSview Help"); n++;
  help_dialog = XmCreateFormDialog(parent, "helpbox", args, n);
  
  /* create button area at bottom */
  /* Note: the stuff at the bottom needs to be put in first due to the
     adjustlast policy of the RowColumn widget */
  dismiss = XtVaCreateManagedWidget("dismiss",
				    xmPushButtonWidgetClass, help_dialog,
				    /*XmNtopAttachment,   XmATTACH_WIDGET,
				      XmNtopWidget,       sep,*/
				    XmNleftAttachment,  XmATTACH_FORM,
				    XmNrightAttachment,  XmATTACH_FORM,
				    XmNbottomAttachment, XmATTACH_FORM,
				    NULL);
  XtAddCallback(dismiss, XmNactivateCallback, 
		(XtCallbackProc)unmanage_helpdialogCB, (XtPointer)NULL);
  sep = XtVaCreateManagedWidget("sep",
				xmSeparatorWidgetClass, help_dialog,
				XmNbottomAttachment,   XmATTACH_WIDGET,
				XmNbottomWidget,       dismiss,
				XmNleftAttachment,  XmATTACH_FORM,
				XmNrightAttachment, XmATTACH_FORM,
				NULL);
  
  row = XtVaCreateWidget("row",
			 xmFormWidgetClass, help_dialog,
			 XmNorientation, XmHORIZONTAL,
			 XmNwidth, textWid+topicWid+20,
			 XmNresizeHeight, True,
			 NULL);
  
  /* Create scrolled list of help topics */
  n = 0;
  XtSetArg(args[n], XmNselectionPolicy, XmSINGLE_SELECT); n++;
  XtSetArg(args[n], XmNtopAttachment, XmATTACH_FORM); n++;
  XtSetArg(args[n], XmNbottomAttachment, XmATTACH_FORM); n++;
  XtSetArg(args[n], XmNheight, textHei); n++;
  XtSetArg(args[n], XmNwidth, topicWid); n++;
  XtSetArg(args[n], XmNresizeHeight, True); n++;
  XtSetArg(args[n], XmNeditable, False); n++;
  help_topic = XmCreateScrolledList(row, "help_topic", args, n);
  XtAddCallback(help_topic, XmNsingleSelectionCallback,
		topic_callback, topic_data);
  XtAddCallback(help_topic, XmNhelpCallback,
		help_callback, help_topic_data);
  XtAddCallback(XtParent(help_topic), XmNhelpCallback,
		help_callback, help_topic_data);
  
  /* Create text Widget to display help */
  n = 0;
  XtSetArg(args[n], XmNeditMode, XmMULTI_LINE_EDIT); n++;
  XtSetArg(args[n], XmNbottomAttachment, XmATTACH_FORM); n++;
  XtSetArg(args[n], XmNtopAttachment, XmATTACH_FORM); n++;
  XtSetArg(args[n], XmNrightAttachment, XmATTACH_FORM); n++;
  /*XtSetArg(args[n], XmNleftAttachment, XmATTACH_WIDGET); n++;
    XtSetArg(args[n], XmNleftWidget, help_topic); n++;*/
  XtSetArg(args[n], XmNheight, textHei); n++; 
  XtSetArg(args[n], XmNwidth, textWid); n++;
  XtSetArg(args[n], XmNresizeHeight, True); n++;
  XtSetArg(args[n], XmNeditable, False); n++;
  help_text = XmCreateScrolledText (row, "help_text", args, n);
  XtAddCallback(help_text, XmNhelpCallback,
		help_callback, help_topic_data);
  XtAddCallback(XtParent(help_text), XmNhelpCallback,
		help_callback, help_topic_data);
  
  
  /* set up scrolled area attachments */
  XtVaSetValues (row,
		 XmNtopAttachment,   XmATTACH_FORM,
		 XmNleftAttachment,  XmATTACH_FORM,
		 XmNrightAttachment, XmATTACH_FORM,
		 XmNbottomAttachment, XmATTACH_WIDGET,
		 XmNbottomWidget,       sep,
		 NULL);
  
  XtManageChild(help_topic);
  XtManageChild(help_text);
  XtManageChild(row);
  return True;
} /* end  HelpBoxCreate */

/**
 * Put help text in display 
 * \param array of text strings to display, NULL terminated
 */
void HelpBoxSetText (char** text) {
  int next=0;
  if (help_text == (Widget)NULL) return; /* sanity check */
  if (text == NULL) return;
  
  /* hide text until it's redrawn */
  if (XtIsRealized (help_text)) XtUnmapWidget(help_text);
  
  /* clear any previous text */
  XmTextReplace(help_text,(XmTextPosition)0, 
		XmTextGetLastPosition(help_text), "");
  
  /* copy text */
  while (1) { /* loop till done */
    if (!strncmp (text[next], "*** FINISHED ***", 16)) break; /* Finished */
    XmTextInsert(help_text, XmTextGetLastPosition(help_text), text[next++]);
  }
  /* start at top */
  XmTextShowPosition(help_text, (XmTextPosition)0);
  
  /* show text */
  if (XtIsRealized (help_text)) XtMapWidget(help_text);
  
} /* end  HelpBoxSetText */

/**
 * Adds topic to end of list
 * \param topic topic to add
 */
void HelpBoxAddTopic (char* topic) {
  XmString xmstring = NULL;
  if (help_topic == (Widget)NULL) return; /* sanity check */
  
  xmstring = XmStringCreateSimple(topic);
  XmListAddItemUnselected(help_topic, xmstring, 0); /* add to end */
  if (xmstring) XmStringFree(xmstring); xmstring = NULL;/* release structure */ 
} /* end HelpBoxAddTopic */

/**
 * Delete all topicsKAdds topic to end of list
 */
void HelpBoxDeleteAllTopics(void) {
  if (help_topic == (Widget)NULL) return; /* sanity check */
  
  XmListDeselectAllItems (help_topic);
  XmListDeleteAllItems (help_topic);
} /* end HelpBoxDeleteAllTopics */

/**
 * Callback for help on a given topic
 * \param w           scroll widget activated
 * \param clientData  cpointer to topic string
 * \param callData    call data
 */
void HelpBoxTopicCB (Widget w, XtPointer clientData, XtPointer callData) {
  HelpBoxShowTopic ((char *)clientData);
} /* end HelpBoxTopicCB */

/**
 * Displays help on a given topic 
 * \param topic topic to display
 */
void HelpBoxShowTopic( char* topic) {
  XmString xmstring = NULL;
  char message[100];
  
  /* make sure it exists before proceeding */
  if (!HelpBoxActive || (help_topic == (Widget)NULL)) return; 
  
  /* in case topic is not found */
  /* hide text until it's redrawn */
  if (XtIsRealized (help_text)) XtUnmapWidget(help_text);
  
  g_snprintf (message, 99,
	   "  Requested topic %s not found\n  Click topic on left to display help.\n",
	   topic);
  XmTextReplace(help_text,(XmTextPosition)0, 
		XmTextGetLastPosition(help_text), message);
  
  /* select topic from list */
  xmstring = XmStringCreateSimple(topic);
  XmListSelectItem(help_topic, xmstring, True); /* select to show */
  XmListSetItem(help_topic, xmstring); /* adjust topic list */
  
  /* in case it's not already visible */
  if (!XtIsManaged(help_dialog)) XtManageChild(help_dialog);
  if (XtIsRealized (help_text)) XtMapWidget(help_text);
  
  
  if (xmstring) XmStringFree(xmstring); xmstring = NULL;/* release structure */
} /* end HelpBoxAddTopic */


/*********** generic help up to this point, below FITSview specific *********/

/**
 * Create/Initialize HelpBox
 * \param parent  parent widget
 * \param visible True if to display
 */
/* Create/Initialize HelpBox */
void MakeHelpBox(Widget parent, Boolean visible) {
  int topic;
  
  /* Create dialog - if returns false then it already exists*/
  if (!HelpBoxCreate (parent, 
		      (XtCallbackProc)HelpBoxSelectCB, 
		      (XtPointer)NULL, 
		      (XtCallbackProc)HelpBoxTopicCB, 
		      (XtPointer)"How to use help", 
		      (XtPointer)"How to use help")) return;
  
  /* initialize help table */
  InitHelpText();
  
  /* add topics */
  for (topic=0; topic<number_topics; topic++) {
    
    HelpBoxAddTopic(topic_title[topic]);
  } /* end of loop over topic list */
  
  
  /* show the box if requested */
  if (visible) HelpBoxShow();
} /* end MakeHelpBox */

/**
 * Callback for HelpBox selection
 * \param w           scroll widget activated
 * \param clientData  client data
 * \param callData    call data
 */
void HelpBoxSelectCB(Widget parent, XtPointer clientData, XtPointer callData) {
  XmListCallbackStruct* ptr = (XmListCallbackStruct*)callData;
  char*                 string;
  int                   topic;
  
  /* which item selected */
  XmStringGetLtoR(ptr->item, XmSTRING_DEFAULT_CHARSET, &string);
  
  /* lookup in list */
  for (topic=0; topic<number_topics; topic++) {
    if (!strcmp(string, topic_title[topic])) break;
  } /* end of loop over topic list */
  
  /* did this fall off the end? - if so return */
  if (topic>=number_topics) return;
  
  /* display */
  HelpBoxSetText(topic_text[topic]);
  
  if (string) XtFree(string); string = NULL;/* free string */ 
} /* end MakeHelpBox */

/**
 * Set help topics and titles and text
 */
void InitHelpText(void) {
  /* intro to help display */
  static  char *helponhelp_text[] = {
    "     Click topic on left to display help.\n",
    "*** FINISHED ***"}; /* end of text */
  
  static char *intro_text[] = {
    "                     ObitView\n",
    "  \n",
    "Please relay comments and/or suggestions to Bill Cotton at NRAO \n ",
    "(bcotton@nrao.edu) \n",
    " \n",
    "This program is a viewer for astronomical images in FITS (Flexible \n",
    "Image Transport System) or AIPS format.  FITS Images in normal or\n",
    "gzip compressed form may be viewed.  An image can be displayed in a \n",
    "number of ways including colorizing the display, zoom and scroll.  \n",
    "In addition, celestial positions can be determined for locations in \n",
    "the image by clicking the left mouse button when the cursor is over \n",
    "the desired feature. \n",
    "  \n",
    "*** FINISHED ***"}; /* end of text */

  /* overview of ObitView */
  static  char *over_text[] = {
    "     Overview \n",
    "   \n",
    "  This viewer will display and manipulate astronomical images in FITS \n",
    "or AIPS format.  The file to view  can be specified as a command line \n",
    "argument or using the Open item in the File menu.   Subsequent FITS files\n",
    "can also be selected using the Open function.   Information about the \n",
    "displayed image can be obtained using the Image info item in the File \n",
    "menu. If the FITS file directory contains an appropriate index, then a \n",
    "celestial position can be entered and ObitView will look up the \n",
    "image containing the position (if any) and load this image centered on \n",
    "the requested position. \n",
    "   \n",
    "   This image browser can either load FITS or AIPS files selected from a\n",
    "filebrowser or as specified through an xmlrpc interface. The xmlrpc \n",
    "interface can also be used to interactively edit CLEAN boxes.\n",
    "    \n",
    "   Either interface is internet enabled for FITS files. File names of\n",
    "the form ftp://ftp.here.there.com/myFile.fits.gz or\n",
    "http://www.FaroukU.edu/~someuser/myFavorite.fits will work.\n",
    "   \n",
    "     The header of a FITS file (or the contents of a text file) may be \n",
    "previewed before deciding which image to load.  Once a file is \n",
    "displayed it can be manipulated and examined in a number of ways.  If \n",
    "the image is larger than the display, the scroll bars on the edge of \n",
    "the display will scroll around inside of the image.  Clicking the left \n",
    "mouse button in the display will result in the brightness and \n",
    "celestial position of the pixel under the cursor being displayed at \n",
    "the bottom of the Display control box.  A click on the right mouse \n",
    "button is similar except that a point model is fitted to the image \n",
    "near the selected pixel; the results are given in the Display control \n",
    "box.  These brightness and position displays can be logged to a text \n",
    "file by selecting the 'Log positions' option in the file menu. \n",
    "   \n",
    "     Standard World Coordinate System (WCS) coordinates are\n",
    "supported as well as the astrometric plate parameters of the Digitized\n",
    "Sky Survey (DSS) and IRAF coordinates.  Positions can be displayed and \n",
    "entered in either equinox B1950 or J2000 if the equinox of the image \n",
    "is either of these.\n",
    "   \n",
    "     The brightness and contrast of the image can be adjusted using \n",
    "the horizonal scroll bars  at the top of the Display control box. \n",
    "Moving the slider to the right will increase the contrast or the \n",
    "brightness.  If the range of pixel brightness  of the portion of the \n",
    "image of interest is significantly smaller than the total range, \n",
    "contrast and brightness adjustments may be insufficient.  In this \n",
    "case, a limited range of pixel values can be displayed using the Pixel \n",
    "Range items in the Options control box,  Alternately, the nonlinear \n",
    "option in the Options menu may display the desired range of \n",
    "brightness.  Blanked pixels always appear as black. By default,\n",
    "ObitView will attempt to guess the proper pixel range to display.\n",
    "   \n",
    "     The image can be displayed in color using one of two color \n",
    "schemes, Color Contour and Pseudo Flame in the Colorize menu.  Color \n",
    "Contour is an 8 color scheme  which gives a contouring effect and \n",
    "Pseudo Flame is a continous color pseudo coloring scheme. giving the \n",
    "image a flame like quality.  Option Grayscale is a black and white \n",
    "coloring scheme.  The order of the color table (black becomes white \n",
    "etc.) is reversed using the Reverse item.   Brightness and contract \n",
    "controls also work on colorized images.  The color, contrast and \n",
    "brightness can be reset using the Reset item on the Colorize menu. \n",
    "   \n",
    "     When an image is initially loaded, generally the first plane in \n",
    "the file is displayed.  If the image contains multiple frequency or \n",
    "polarization planes, other planes can be loaded using the Plane number \n",
    "item in the Options control box.  The number and type of planes in the \n",
    "file can be determined using the Image Info item in the File menu.  A \n",
    "cube can be displayed as a movie (using the 'Movie' item in the Movie \n",
    "menu) to show a range of planes in sequence or by selecting planes at \n",
    "random. Locations in dimensions higher than 3 can be selected on either \n",
    "the Option or Movie dialogs using the text fields labeled   \n",
    "'Higher dimensions' \n",
    "   \n",
    "     An image can be zoomed in or out using the Zoom menu and \n",
    "selecting the desired magnification factor.  Zooming in (factor > \n",
    "100%) is done by replicating pixels and zooming out (magnification \n",
    "<100%)  by displaying only a subset of the pixels.  Zooming is \n",
    "centered on the current scroll position controlled by the image scroll \n",
    "bars.  Selecting a zoom factor of 100% undoes the effects of zooming. \n",
    "   \n",
    "     Celestial positions determined from right mouse clicks will be \n",
    "refined by fitting a point model to the position selected.  This will \n",
    "fit an accurate position and flux assuming a point object near the \n",
    "position of the mouse click.  The results will be displayed in the \n",
    "Control Panel.  The Mark Position item in the Position menu will \n",
    "bring up a dialog box in which the celestial coordinates of an object \n",
    "of interest can be entered; alternately a list of positions can be \n",
    "given in a file.  The corresponding location(s) on the image will be \n",
    "marked. If the current FITS directory contains a special index \n",
    "(named 'findex.txt'), then the Lookup Position item in the Position \n",
    "menu can be used to find the FITS image containing that position and \n",
    "load it. \n",
    "   \n",
    "     Two images can be compared using the Blink facility invoked by \n",
    "the Blink menu.  Blinking will alternately display one image and then \n",
    "the other.  The first image is loaded into the display and desired \n",
    "adjustments are made.  It is then copied into the Blink image using \n",
    "the 'Swap Blink and Current' item in the Blink menu.  The second image \n",
    "is then loaded into the display and adjusted as desired.  The Blink \n",
    "images item on the Blink menu will then begin blinking.  The dwell \n",
    "time on each image can be controlled using the scroll bar in the blink \n",
    "dialog box.  The Quit button on the dialog box ends blinking.  If the \n",
    "two images have pixels coincident on the sky, the zoom and scroll used \n",
    "are that for the current display (the one visible before the blink \n",
    "starts).  If the pixels are not aligned, blinking uses the scroll, \n",
    "zoom and display setup for the blink image that were in effect when it \n",
    "was copied to the Blink image and the current setup for the second \n",
    "(current) (normal display) image before the blink began.  The 'Swap \n",
    "Blink and Current' item swaps the current and blink images. \n",
    "*** FINISHED ***"} ; /* end of text */   
  
  /* file menu */
  static  char *file_text[] = {
    "--------------------------- File Menu --------------------------------- \n",
    "   \n",
    "  \n",
    "  This menu includes a number items related to files.  Note: the \n",
    "following may read gzip compressed files. \n",
    "   \n",
    "*** FINISHED ***"} ; /* end of text */   
  
  static char*fileopen_text[] = {
    "  Open FITS\n",
    "   \n",
    "     This item will bring up a file browser dialog box to select the \n",
    "FITS file to load.  When a file is loaded, the previous image is \n",
    "discarded.  The title bar of the main window gives the name of the \n",
    "currently loaded file. \n",
    "     When the default Pixel Range (0, 0), the default, is specified, a  \n",
    "first pass is made through the plane to determine an appropriate range of  \n",
    "pixel values to display.  This decision is based on a histogram to  \n",
    "determine the sky and noise levels. \n",
    "     Select a FITS image from the browser and click the OK button to \n",
    "load the image.  When the file is being loaded to the display, a box \n",
    "appears with a progress message and a cancel button.  If the message \n",
    "'WARNING: BAD PIXEL RANGE' appears then all of the pixels loaded are \n",
    "at one extreme of the range of displayed brightness.  This usually \n",
    "indicates inappropriate values in the Set Pixel Range option in the \n",
    "Options dialog.  This may be the result of a previous image with a very \n",
    "different range of pixel values.  Setting both values to zero will get \n",
    "the default display.  \n",
    "   \n",
    "*** FINISHED ***"} ; /* end of text */   
  
  static char*afileopen_text[] = {
    "  Open AIPS\n",
    "   \n",
    "     This item will bring up an image browser dialog box to select the \n",
    "AIPS image to load.  When an image is loaded, the previous image is \n",
    "discarded.  The title bar of the main window gives the name of the \n",
    "currently loaded file. \n",
    "     To select an image, set the AIPS user number and path to the AIPS \n",
    "directory.  Then either the full specification of the AIPS image can be\n",
    "entered in the textboxes at the bottom or after entering the user number \n",
    "and directory name, activate the directory name box (hit return) if needed\n",
    "to obtain a listing of the images. Clicking on a line in the image list\n",
    "will select it. Loading an image from ObitTalk saves the user & directory\n",
    "The Open button will cause the image to be displayed.\n",
    "The Cancel button will cancel the dialog.\n",
    "The Help button gets this help description.\n",
    "     When the default Pixel Range (0, 0), is specified, a  first pass \n",
    "is made through the plane to determine an appropriate range of  pixel \n",
    "values to display.  This decision is based on a histogram to determine  \n",
    "the sky and noise levels. \n",
    "   \n",
    "*** FINISHED ***"} ; /* end of text */   
  
  static char*filepreview_text[] = {
    "  Preview \n",
    "   \n",
    "     This item will display up to about 1000 lines of the header of a \n",
    "selected FITS file or text file in a scroll box.  This allows deciding \n",
    "which image to load or reading explanatory text.  The first line of \n",
    "the header of a FITS header is:  \n",
    "  SIMPLE  =                     T     /possibly some comment \n",
    " \n",
    "  A FITS image file has a line near the beginning of the form: \n",
    "  NAXIS1   =              nnn        /possibly some comment \n",
    "   \n",
    "  where nnn is an integer larger than 0.  The size of the image is \n",
    "given by the NAXIS1 and NAXIS2 entries.  Information about the image \n",
    "may be contained in HISTORY or COMMENT lines. \n",
    "   \n",
    "*** FINISHED ***"} ; /* end of text */   
  
  static char*filesaveas_text[] = {
    "  Save as \n",
    " \n",
    " NYI \n",
    "     This item will copy the currently displayed FITS image into a \n",
    "FITS file to be specified in a file specification dialog box.  The file \n",
    "is written in uncompressed form irregardless of the compression state of \n",
    "the input file. \n",
    "  \n",
    "*** FINISHED ***"} ; /* end of text */   
  
  static char*fileinfo_text[] = {
    "  Image info \n",
    "   \n",
    "     This item will display information about the current image \n",
    "including positions, frequencies, observation dates, etc.  The Dismiss \n",
    "button clears this display; Refresh updates the display for the \n",
    "currently loaded image. \n",
    "   \n",
    "*** FINISHED ***"} ; /* end of text */   
  
  static char*filelog_text[] = {
    "  Log positions \n",
    "   \n",
    "     This will toggle the logging of brightnesses and positions \n",
    "selected by a mouse button click or fitting a point model to the \n",
    "image.  When this is turned on, a dialog box will allow selection of \n",
    "the text file.  The logging file contains one line per position \n",
    "containing 1) the pixel location on the first three axes, 2) the \n",
    "celestial position and equinox of the first two axes, 3) the corresponding \n",
    "brightness from the image and,  4) the name of the FITS file.  \n",
    "If the menu item is labeled 'Start Position Logging' then logging is \n",
    "currently disabled and selecting this item will turn it on. \n",
    "Conversely, if the item is labeled 'Stop Position Logging' then \n",
    "logging is currently active and selecting this item will turn it off. \n",
    "   \n",
    "*** FINISHED ***"} ; /* end of text */   
  
  static char*filequit_text[] = {
    "  Quit \n",
    "   \n",
    "     This will terminate the program. \n",
    " \n",
    "*** FINISHED ***"} ; /* end of text */   
  
  static char*filesaveimwindow_text[] = {
    "  Save Image Window as \n",
    "   \n",
    "     If an image window has been attached to the currently displayed  \n",
    " image, it is written as python to a text file.   \n",
    " The name of the file is provided via a file browser. \n",
    " \n",
    "*** FINISHED ***"} ; /* end of text */   
  
  static char*filehelp_text[] = {
    "  Help \n",
    "    Displays information about the file menu. \n",
    "*** FINISHED ***"} ; /* end of text */
  
  /* Options menu */
  static char *options_text[] = {
    "--------------------------- Options Menu ------------------------------ \n",
    "   \n",
    "  This menu item will bring up the Options control box  \n"
    "The first item is a text box which allows setting the timeout for editing  \n",
    "requests.  These are requests sent by other software expecting a user  \n",
    "response.  If none is received in the timeout period, the OK button is  \n",
    "activated. The timeout can be entered in the text box in seconds and will\n",
    "become active when the OK button is hit.   \n",
    "A nonpositive value means wait forever; minimum timeout = 5 sec..  \n",
    "  \n",
    "Further items control the range of values loaded into the display and the \n",
    "plane in the image loaded.  These items are discussed below. \n",
    "   \n",
    "Pixel Range \n",
    "     The range of values that are present in the image plane displayed \n",
    "are shown.  The first value entered (labeled 'Minimum pixel value') is \n",
    "the minimum pixel value to be displayed and the second ('Maximum pixel \n",
    "value') is the maximum to be displayed.  Pixel values below the \n",
    "minimum are set to the minimum and above the maximum are set to the \n",
    "maximum.  0, 0 means let ObitView decide the pixel range; the default \n",
    "should be adequate for most purposes but uses assumptions about the\n",
    "image which are true for many astronomical images but may not be\n",
    "valid especially if the features in the image are fainter than the\n",
    "'sky' level (e.g. absorption, Stokes Q, U or V or negative images).\n",
    "   \n",
    "Plane \n",
    "     If the image contains multiple planes (frequency, polarization \n",
    "etc.) the desired plane number can be specified using this item.  The \n",
    "dialog box will initially contain the current plane number and will \n",
    "tell the allowed range of values (1 - n).  Information about the \n",
    "number and type of planes may be obtained from the Image info item in \n",
    "the File menu.  Planes are numbered 1 relative. \n",
    "If many planes are to be viewed, use Movie from the Movie menu.\n",
    "   \n",
    "Higher Dimensions \n",
    "     If the image contains more than three non degenerate (dimension > 1) \n",
    "axes, the (1-rel) location on dimensions 4-6 can be entered in the boxes   \n",
    "labeled 'Higher dimensions'.  The window label also indicates the location  \n",
    "on higher dimensions. \n",
    "   \n",
    "Linear display (radio button) \n",
    "     This option specifies a linear mapping of image pixel values to \n",
    "display colors.  \n",
    "   \n",
    "Nonlinear display (radio button) \n",
    "     This option specifies a nonlinear mapping of pixel values to \n",
    "display colors; the mapping function is the square root of the pixel \n",
    "value.  This option is useful for displaying an image with a large \n",
    "range of pixel values and uses more levels to display low brightness \n",
    "than the linear mapping. \n",
    "   \n",
    "Histogram equalization (radio button)\n",
    "     This option specifies a histogram equalization of the pixels in\n",
    "the relevant pixel range (specified or ObitView default).  This\n",
    "attempts to have equal numbers of pixels (within the range) shown in\n",
    "each of the colors of the display.  This option may be useful for\n",
    "displaying an image with interesting structure over a wide range of\n",
    "brightness. \n",
    "  \n",
    "OK \n",
    "     This reads the current values and saves them for the next image \n",
    "load and dismisses the dialog box. \n",
    " \n",
    "Cancel \n",
    "     Dismisses the dialog box with no changes to the loading \n",
    "parameters.  \n",
    "   \n",
    "Reload \n",
    "     This option causes the currently selected image to be reloaded \n",
    "into the display using the current set of options.  This needs to be \n",
    "done after any of the other options in this dialog box have been \n",
    "changed in order for these changes to take effect. \n",
    "   \n",
    "     When the file is being loaded to the display, a box appears with \n",
    "a progress message and a cancel button.  If the message 'WARNING: BAD \n",
    "PIXEL RANGE' appears then all of the pixels loaded are at one extreme \n",
    "of the range of displayed brightness.  This usually indicates \n",
    "inappropriate values in the Set Pixel Range in the Options menu.  This \n",
    "may be the result of a previous image with a very different range of \n",
    "pixel values.  Setting both values to zero will get the default \n",
    "display.  Hitting the Cancel button cancels loading the image. \n",
    "*** FINISHED ***"} ; /* end of text */   
  
  static char*optionshelp_text[] = {
    "  Help \n",
    "    Displays information about the options menu. \n",
    "   \n",
    "*** FINISHED ***" }; /* end of text */
  
  /* Zoom menu */
  static char *zoom_text[] = {
    "--------------------------- Zoom Menu ------------------------------ \n",
    "   \n",
    "  This menu controls zooming of the image. \n",
    "   \n",
    "  25% \n",
    "   \n",
    "     This will reduce the size of the displayed image to 1/4 of its \n",
    "normal size by discarding 3 out of 4 rows and columns. \n",
    "   \n",
    "  50% \n",
    "   \n",
    "     This will reduce the size of the displayed image to 1/2 of its \n",
    "normal size by discarding alternate rows and columns. \n",
    "   \n",
    "  100% \n",
    "   \n",
    "     Resets the zoom to its initial setting. \n",
    "   \n",
    "  200% \n",
    "   \n",
    "    This magnifies the image by a factor of two by replicating pixels. \n",
    "   \n",
    "   \n",
    "  400% \n",
    "   \n",
    "    This magnifies the image by 400%. \n",
    "   \n",
    "  800% \n",
    "   \n",
    "     This magnifies the image by 800%. \n",
    "   \n",
    "  1600% \n",
    "   \n",
    "    This magnifies the image by 1600%. \n",
    "   \n",
    "*** FINISHED ***"} ; /* end of text */   
  
  static char*zoomhelp_text[] = {
    "  Help \n",
    "    Displays information about the zoom menu. \n",
    "   \n",
    "     \n",
    "*** FINISHED ***" }; /* end of text */
  
  /* Position menu */
  static char *position_text[] = {
    "--------------------------- Position Menu ---------------------------- \n",
    "   \n",
    "  This menu contains functions related to celestial position. \n",
    "  \n",
    "*** FINISHED ***"} ; /* end of text */   
  
  static char* positionsetequinox_text[] = {
    "  Set Equinox \n",
    "    \n",
    "   This option allows specifying the equinox (B1950 or J2000) of the \n",
    "celestial coordinates displayed or entered.  Thus, if the image is in \n",
    "J2000 coordinates and you have a position in B1950, clicking on the \n",
    "B1950 button in the Set Equinox dialog box will cause positions \n",
    "displayed after a mouse click in the image to be equinox B1950. \n",
    "Furthermore, the position specified to Mark Position will be B1950. \n",
    "See the 'Source Info' box in the file menu to determine the equinox of \n",
    "the image. \n",
    "   The dialog box invoked by this menu item has three options: 1) use \n",
    "the equinox of the image (default), 2) equinox J2000 and 3) equinox \n",
    "B1950.  Note if the image does not specify the equinox then this \n",
    "selection has no effect. \n",
    "   \n",
    "*** FINISHED ***"} ; /* end of text */   
  
  static char* positionmark_text[] = {
    "  Mark Position \n",
    "   \n",
    "     This option lets you mark a particular celestial position in the \n",
    "current image. Selecting this option brings up a dialog box into which \n",
    "the desired celestial position is entered; optionally, the name of a \n",
    "text file can be given with a list of positions.  The selected \n",
    "positions are marked in the displayed image by a cross whose inner \n",
    "positions are not shown so as not to obscure the image.  The cross is \n",
    "marked in the display by replacing the previous values with that of \n",
    "the brightest pixel in the image.  These markers will persist until \n",
    "the display is reloaded with the same or another image.  Note that the \n",
    "cross may not be visible on some reduced zoom displays as the marked \n",
    "pixels may not be among those shown.  When the desired values are \n",
    "entered into this dialog box, the 'Mark' button will proceed to mark \n",
    "the image and set the scroll to center the last position marked.  The \n",
    "'Cancel' button dismisses the box without marking any positions. \n",
    "If only 1 value is given for ra and dec, it is assumed to be the pixel\n",
    "numbers to mark.\n",
    "   \n",
    "     The size of the cross can be controlled by the two values in the \n",
    "line labeled 'size'.  The first of these is the 'inner' size of the \n",
    "cross or the distance from the center in pixels in which the cross is \n",
    "not marked.  The second value is the 'outer' size or the distance in \n",
    "pixels from the center over which the cross appears. \n",
    " \n",
    "     The 'Clear' button will blank the position text boxes.  If the 'DEC'  \n",
    "box is blank, then the full RA and Dec strings can be entered in the RA \n",
    "field. \n",
    "   \n",
    "     The 'Swap/Mark' button will swap the current and blink images and\n",
    "mark the indicated position in that image.\n",
    "   \n",
    "     If many positions are to be marked, they can be entered in a \n",
    "text file prepared by a text editor (e.g. emacs).   Each line of this \n",
    "file should have an entry of the form:  \n",
    "   \n",
    "  RA:h RA:m RAs  Dec:d Dec:m Dec:s inner outer \n",
    "   \n",
    "  where RA:h, RA:m, RA:s are the hours, minutes and seconds of the \n",
    "Right Ascension and Dec:d, Dec:m, and Dec:s are the degrees (with \n",
    "-sign if in the south), minutes, and seconds of the Declination. \n",
    "Inner and outer are the inner and outer sizes of the marking cross in \n",
    "pixels.   \n",
    "  An example is: \n",
    "   \n",
    "    12 34 23.7898  -15 23 45.634 3 10 \n",
    "   \n",
    "     To select a file containing positions, hit the 'file' button for \n",
    "a directory browser to specify the text file.  If there are positions \n",
    "in the file that are out of the image then the number of these \n",
    "positions are reported and the remainder marked.  \n",
    "   \n",
    "   Note: if the inner size of the cross is -1 a circle is drawn and the  \n",
    "window editing function is entered.   \n",
    "*** FINISHED ***"} ; /* end of text */   
  
  static char* positionlookup_text[] = {
    "  Lookup Position \n",
    " \n",
    "    If the current FITS file directory contains a special index, this \n",
    "item allows specifying a celestial position and ObitView will \n",
    "determine which image (if any) contains that position, load the image, \n",
    "and center the display on that position.  The equinox of the position \n",
    "can be specified using the 'Set Equinox' option in this menu.  Once \n",
    "the desired position is entered in the dialog, the Lookup button will \n",
    "cause ObitView to attempt to find and load the desired image.  The \n",
    "Cancel button dismisses the dialog with no change to the current image. \n",
    " \n",
    "   The index file in the FITS file directory should have name \n",
    "'findex.txt' and obey the following rules. \n",
    " \n",
    "line 1: \n",
    "   This line must specify the equinox of the image field centers, e.g. \n",
    "equinox 1950. \n",
    " \n",
    "lines 2...: \n",
    "   These lines give the name of the file, the central position, the \n",
    "half width (degrees) in RA and Declination and a priority number.  The \n",
    "file containing the specified position with the highest priority \n",
    "number is the one selected.  Comments may be given by putting an '!' \n",
    "in the first column.  The entry is free format (but must be confined \n",
    "to a single line) and the structure of a file entry is: \n",
    " \n",
    "filename hh mm ss.s -dd mm ss.s hw.ra  hw.dec  priority \n",
    " \n",
    "Examples follow \n",
    "! These images are from the long lost Summarian clay tablets \n",
    "!name              RA             Dec   delt_ra delt_dec priority \n",
    "sky.fit        11 23 56.234  82 17 56.7  1.87     .25      1 \n",
    "south.fit.gz   17 23 54.1   -74  4 53    3.65     1        2 \n",
    " \n",
    "   Note: the half width in RA is half the number of cells in RA times \n",
    "the cell spacing divided by the cosine of the declination to account \n",
    "for the converging lines of RA towards the poles. \n",
    "   \n",
    "*** FINISHED ***"} ; /* end of text */   
  
  static char* positionfit_text[] = {
    "  Fit Position \n",
    " \n",
    "    A point model will be fitted near the position of the last left \n",
    "mouse button click in the image.  Clicking the right mouse button will\n",
    "cause a point model to be fitted at the current position.  In either case,\n",
    "the Display control box will show the fitted position to subcell accuracy \n",
    "as well as the brightness interpolated to that position.  If the fitting is\n",
    "successful, 'fitted' will appear in the Display control else 'fit failed'\n",
    "is shown.  If logging is enabled then the fitted position is logged.\n",
    " \n",
    "*** FINISHED ***"} ; /* end of text */   
  
  static char* positionhelp_text[] = {
    "  Help \n",
    "    Displays information about the position menu. \n",
    "   \n",
    "    \n",
    "*** FINISHED ***" }; /* end of text */
  
  /* Blink menu */
  static char *blink_text[] = {
    "--------------------------- Blink Menu ---------------------------- \n",
    "   \n",
    "  This menu controls blinking. \n",
    "   \n",
    "*** FINISHED ***"} ; /* end of text */   
  
  static char* blinkswap_text[] = {
    "  Swap blink and current \n",
    "   \n",
    "     This will swap the current and blink images.  An image must be \n",
    "copied to the blink image for blinking to be effective.  This item \n",
    "allows changing the color table of the blink image or examining values \n",
    "or positions in that image.  This can be used repeatedly. \n",
    "   \n",
    "*** FINISHED ***"} ; /* end of text */   
  
  static char* blinkblink_text[] = {
    "  Blink images \n",
    "   \n",
    "     This will bring up a dialog control box and start blinking the \n",
    "images.  The dwell time is controlled using a scroll bar and the Quit \n",
    "button terminates blinking.  The title bar of the main window gives \n",
    "the name of the currently displayed file.  When blinking stops, the \n",
    "(previously) current image is displayed.  If the images have aligned \n",
    "pixels on the sky (the only case that makes sense) then the zoom and \n",
    "scroll of the blink image is forced to that of the current image. \n",
    "Otherwise, the zoom and scroll are those set for each of the images \n",
    "and there may be no correspondence between the pixels of the two \n",
    "images. \n",
    " \n",
    "*** FINISHED ***"} ; /* end of text */   
  
  static char* blinkhelp_text[] = {
    "  Help \n",
    "    Displays information about the blink menu. \n",
    "     \n",
    "*** FINISHED ***" }; /* end of text */
  
  /* Movie menu */
  static char *movie_text[] = {
    "--------------------------- Movie Menu ---------------------------- \n",
    "   \n",
    "  Movie \n",
    "   \n",
    "     This will bring up a dialog box which controls displaying planes \n",
    "in a movie-like fashion where the display is periodically updated with \n",
    "the next plane.  Planes can be shown as movies or selected manually \n",
    "using the scroll bar.  The current plane is indicated by the location \n",
    "of the slider in the scroll bar and the text lines under it giving the \n",
    "plane number and the value along this axis in the cube.  The movie \n",
    "function is controlled by the values in the text boxes labeled \n",
    "'planes' (the start plane for the movie), 'to' (the final plane for \n",
    "the movie), and 'Dwell (sec)' (the dwell time on each frame in \n",
    "seconds).  These values may be modified by clicking on the box and \n",
    "typing in new values. The movie is started using the 'Play' button and \n",
    "can be stopped prematurely by the 'Stop' button.  The movie proceeds \n",
    "through the selected range once.  NB: the speed of the movie may be \n",
    "slower than indicated by the 'Dwell' value if it takes longer than \n",
    "this to load the next plane from the disk.  The displayed plane can be \n",
    "controlled manually using the scroll bar.  The selected plane will \n",
    "remain displayed until another plane is selected.  The 'Quit' button \n",
    "exits movie mode and resumes the normal display. \n",
    "   Looping is always over the third dimension but the (1-rel) location\n",
    "on dimensions 4-6 can be entered in the boxes labeled 'Higher dimensions'.\n",
    "The window label also indicates the location on higher dimensions.   \n",
    "   \n",
    "*** FINISHED ***"} ; /* end of text */   
  
  static char* moviehelp_text[] = {
    "  Help \n",
    "    Displays information about the movie menu. \n",
    "   \n",
    "*** FINISHED ***" }; /* end of text */
  
  /* Colorize menu */
  static char *colorize_text[] = {
    "--------------------------- Colorize Menu ---------------------------- \n",
    "   \n",
    "  This menu controls the colorizing of the image which is \n",
    "intrinsically monochromatic.  The brightness and contrast controls \n",
    "will modify the color schemes.  Note: blanked pixels are always \n",
    "displayed as black. \n",
    "   \n",
    "*** FINISHED ***"} ; /* end of text */   
  
  static char* colorizecolor_text[] = {
    "  Color Contour \n",
    "   \n",
    "     This scheme uses a small number of colors to represent the image. \n",
    "This gives a color contour effect. \n",
    "   \n",
    "*** FINISHED ***"} ; /* end of text */   
  
  static char* colorizeflame_text[] = {
    "  Pseudo Flame \n",
    "   \n",
    "     This scheme uses a continous set of colors and intensities to \n",
    "represent the image with a pseudo coloring scheme giving the image a \n",
    "flame like quality.  \n",
    "   \n",
    "*** FINISHED ***"} ; /* end of text */   
  
  static char* colorizegray_text[] = {
    "  Grayscale \n",
    "      \n",
    "     This function uses shades of gray to represent the image. \n",
    "   \n",
    "*** FINISHED ***"} ; /* end of text */   
  
  static char* colorizereverse_text[] = {
    "  Reverse colors \n",
    "      \n",
    "     This function reverses the order of the color table causing \n",
    "(nearly) black to become white etc. Blanked pixels still appear black. \n",
    "   \n",
    "*** FINISHED ***"} ; /* end of text */   
  
  static char* colorizereset_text[] = {
    "  Reset colors \n",
    "   \n",
    "     This resets the colors to shades of gray and resets the \n",
    "brightness and contrast controls.  \n",
    "   \n",
    "*** FINISHED ***"} ; /* end of text */   
  
 static char* cleargraphics_text[] = {
    "Clear Graphics \n",
    " \n",
    "   Clears the graphics overlay (lines over image). \n",
    " \n",
    "*** FINISHED ***"} ; /* end of text */   

 static char* setgraphicscolor_text[] = {
    "Set Graphics Color \n",
    " \n",
    "   Gives 6 choices for the color of graphics overlay: \n",
    "   Magenta, Yellow, Red, Cyan, Blue, Green \n",
    "Any subsequent darwing (mark position or window editing) in this image \n",
    "buffer will be done in this color.  (Note: setting the inner Mark \n",
    "Position size to zero does a test of the window editing). \n",
    " \n",
    "*** FINISHED ***"} ; /* end of text */   

  static char* colorizehelp_text[] = {
    "  Help \n",
    "    Displays information about the Colorize menu. \n",
    "*** FINISHED ***" }; /* end of text */
  
  /* Help menu  */
  static char *help_text[] = {
    "--------------------------- Help Menu ---------------------------- \n",
    "   \n",
    "  This menu controls informative displays about the program. \n",
    " \n",
    "  About ObitView \n",
    "     Gives information about this program. \n",
    "   \n",
    "  Help \n",
    "   \n",
    "     Displays this information. \n",
    "   \n",
    "*** FINISHED ***" }; /* end of text */
  

  /* (Clean) Window editing */
  static char *windowedit_text[] = {
    "                   Interactive Window Editing \n",
    "   ObitView can be used as a data processing image display server to \n",
    "interactively edit a set of windows on an image.  These are frequently\n",
    "used in the CLEAN deconvolution algorithm but is not limited to this use.\n",
    "When a client program displays an image and a window to be possibly  \n",
    "edited, the \"Request\" dialog is displayed and the client program waits \n",
    "for a response.  This box can be used to control the editing of a window \n",
    "and pass instructions back to the client program.  \n",
    "Note: outer windows and unboxes (deselected) are shown in different colors \n",
    "Requests are entered using the following radio buttons at the top of the\n",
    "dialog:\n",
    "  Continue           Continue program\n",
    "  Abort              Abort program - treat as an error\n",
    "  Quit operation     Stop current interation\n",
    "  Turn Off TV        Send no more images to the display\n",
    "  View field         Send the image from the field whose \n",
    "                     number (1-rel) is in the following text field\n",
    "                     and allow editing its window\n",
    "This is followed by a text field in which the, 1-relative, field number\n",
    "can be entered, the maximum possible value is given in the label.\n",
    "\n",
    "   The next row contains two buttons that control the editing of windows.\n",
    "Edit -  Begin interactive edit of window on current image.  \n",
    "        Instructions and other information are displayed in the dialog box.\n",
    "        This editing allows adding new windows (rectangles and circles)\n",
    "        and modification or deletion of existing windows.  \n",
    "        Note: some line segments may not be displayed with zoom factors\n",
    "        25% and 50%.\n",
    "   Unwindows: Normal windows are regions which are selected.\n",
    "        An unwindow is one which is deselected and deselection takes\n",
    "        precedence over selection.  Unwindows are displayed in a different\n",
    "        color than regular inner windows or outer windows.\n",
    "   Outer Window: This is a region of an image in which automatic window.\n",
    "        assignment is allowed to operate.  These are in a different color\n",
    "        from inner windows.\n",
    "        Outer windows cannot be edited.\n",
    "Clear - Delete all inner windows\n",
    "\n",
    "   The interactive session is ended and the client program sent your\n",
    "request using the buttons on the bottom line.\n",
    "OK       Return the request in the radio buttons and the edited window\n",
    "         to the calling client.\n",
    "Cancel   cancel edit - continue process with no change to window. \n",
    "Help     Give this help text\n",
    "\n",
    "Timeout: A timeout period can be set using the Options dialog to specify\n",
    "a timeout period.  If there is no user activity in this period, ObitView \n",
    "will tell the communication program to continue.  \n",
    "A negative timeout value means wait forever.\n",
    "Any user request will disable the timeout for that edit only.\n",
    "*** FINISHED ***" }; /* end of text */

  /* Xmlrpc Interface */
  static char *xmlrpcinterface_text[] = {
    "   ObitView is an xmlrpc server listening to port 8765 (subject to \n",
    "change).  Currently (Jul 07) the following functions exist: \n",
    " \n",
    "ping         sent arbitrary integer, returns ObitView \n",
    "loadFITS     sent filename, returns completion/failure message \n",
    "loadImage    sent file info, returns completion/failure message \n",
    "editWindow   sent Obit Clean window, returns edited version \n",
    "             or failure message \n",
    "copyFile     Copy a file from a remote host to ObitView working\n",
    "             directory \n",
    " \n",
    "clientTest.c contains test client software. \n",
    "*** FINISHED ***" }; /* end of text */
  
  /* Display control */
  static char *displaycontrol_text[] = {
    " \n",
    "   Display control \n",
    "   \n",
    "  The display control contains scroll bars to control the brightness \n",
    "and contrast of the displayed image and information about pixels \n",
    "selected in the image.  \n",
    "   \n",
    "*** FINISHED ***" }; /* end of text */

  /* Graphics overlay */
  static char *graphicsoverlay_text[] = {
    " \n",
    "   Line drawing may be superposed on images.  This is done using the \n",
    "graphics overlay.  These can be erased or color selected using items \n",
    "in the Colorize menu.\n",
    "   \n",
    "*** FINISHED ***" }; /* end of text */
  
  /* Image position */
  static char *imageposition_text[] = {
    " \n",
    "   Image position \n",
    "   \n",
    "     The celestial position and brightness of a given pixel can be \n",
    "determined by clicking the left mouse button when the cursor is on the \n",
    "desired position in the image display.  The results are shown at the \n",
    "bottom of the Display control.  A more accurate position may be \n",
    "obtained for small objects using the right mouse button to select the \n",
    "pixel.  The initial position for the fitting must be within two pixels \n",
    "of the local maximum (or minimum) being fitted. \n",
    "   \n",
    "*** FINISHED ***" }; /* end of text */
  
  /* Image scrolling */
  static char *imagescrolling_text[] = {
    "   \n",
    "     Image scrolling       \n",
    " \n",
    "     If the displayed image is larger than the display area there will \n",
    "be scroll bars on the display area.  These scroll bars can be used to \n",
    "move the visible area around on the image. \n",
    "   \n",
    "*** FINISHED ***" }; /* end of text */
  
  /* Image enhancement */
  static char *imageenhancement_text[] = {
    "   \n",
    "    Image enhancement \n",
    "   \n",
    "     The horizional scroll bars in the Display control box set the \n",
    "brightness and contrast of the image.  Moving the slider towards the \n",
    "right will increase brightness or contrast.  The scroll bars are \n",
    "labeled 'Brightness' and 'Contrast' and the value in parentheses are \n",
    "relative values between 0 and 255.   \n",
    " \n",
    "*** FINISHED ****" }; /* end of text */
  
  /* Browser */
  static char *browser_text[] = {
    "    File Browser (FITS) \n",
    "   \n",
    "File selection uses a standard file browser dialog box.  The selected \n",
    "file should be entered into the item labeled 'Selection'.  This can be \n",
    "done either by clicking on the desired file in the 'Files' box or \n",
    "clicking in the 'Selection' box and entering the name manually.   \n",
    "   The directory displayed can be modified using the 'Directories' box \n",
    "by clicking on the desired subdirectory name (or .. for the parent \n",
    "directory) and hitting the 'Filter' button.  Alternatively, the \n",
    "directory can be changed by modifying the contents of the 'Filter' box \n",
    "to the desired directory and a wildcard match for the file names, e.g. \n",
    "/home/mydir/* \n",
    "selects all files in directory /home/mydir.  Hitting the 'Filter' \n",
    "button will update the contents of the 'Directories' and 'Files' \n",
    "boxes.  The 'filter' box can be used to display only certain files; \n",
    "for example to see only files whose names end in .fits use: \n",
    "/home/myfits/*.fits \n",
    "The 'Filter' button will update the display.  Changing the directory by\n",
    "editing the string in the selection box may NOT have the desired effect.\n",
    "   When the desired file is entered in the 'Selection' box the 'OK' \n",
    "button invokes the requested action.  The 'Cancel' button dismisses \n",
    "the file selection dialog and cancels the requested action.  The \n",
    "'Help' button produces this display. \n",
    " \n", 
    "*** FINISHED ****" }; /* end of text */
  
  /*  Glossary */
  static char *glossary_text[] = {
    "--------------------------- Glossary ---------------------------- \n",
    " \n",
    "  Blanked pixels \n",
    "   \n",
    "     If an image has no measured value associated with a given \n",
    "pixel it is said to be blanked or invalid.  ObitView displays these \n",
    "as black. \n",
    " \n",
    "  Blinking \n",
    "   \n",
    "     Blinking is a technique for comparing images.  If the pixels in \n",
    "two images are aligned and the two are repeatedly displayed one after \n",
    "another, the details of the two can be compared. \n",
    "   \n",
    "  Celestial Position \n",
    "   \n",
    "     Celestial positions on the sky are similar to latitude and \n",
    "longitude used to measure position on the earth.  A celestial position \n",
    "consists of a 'Right Ascension', usually abreviated RA, which \n",
    "corresponds to longitude and a 'Declination', abreviated Dec, \n",
    "corresponds to latitude.  Declination is measured in degrees north \n",
    "and south of the celestial equator, the projection of the earth's \n",
    "equator onto the sky.  Right Ascension is measured in time units, \n",
    "hours, minutes and seconds of sidereal time of transit at Greenwich. \n",
    "     The Earth's rotation axis wobbles with a 25,000 year period due \n",
    "to precession which causes the apparent position of a object to change \n",
    "with time.  Celestial positions are therefore usually expressed in \n",
    "terms of the earth's orientation at a set of standard times called \n",
    "Equinoxes.    The current standard equinoxes are B1950 and J2000 \n",
    "corresponding to the beginnings of the years 1950 and 2000.  The J \n",
    "and B refer to the set of conventions used to 'precess' the \n",
    "coordinates (change them to another time). \n",
    "   \n",
    "   \n",
    "  Color Table \n",
    "   \n",
    "     Image displays show images by representing the value of each \n",
    "pixel by a color or gray shade. The correspondence between the pixel \n",
    "values and the color displayed is called the color table. \n",
    "   \n",
    "   \n",
    "  Image plane \n",
    "   \n",
    "     The simplest images consist of a single two dimensional array of \n",
    "pixels.  An image may contain several (or many) of these 2-D arrays or \n",
    "planes.  Each of the planes can be displayed as an image.  These \n",
    "planes may represent the same region of the sky at different \n",
    "times, frequencies of light, or different polarization states of the \n",
    "light. \n",
    "   \n",
    "  \n",
    "  Pixel \n",
    "   \n",
    "     A pixel is a cell in an image; its value corresponds to the \n",
    "brightness of the image at that position.  In astronomical images a \n",
    "pixel corresponds to a location on the sky.  In images with more than \n",
    "two dimensions, pixels are sometimes called voxels.  \n",
    "   \n",
    " \n",
    "  Precession \n",
    "   \n",
    "     Precession is the wobbling of the earth's rotation axis due to \n",
    "the gravational field of the sun and moon.  This effect is like the \n",
    "wobbling of a top as it slows down.  Earth's precession takes about \n",
    "25,000 years for each cycle. \n",
    "   \n",
    "   \n",
    "  Scrolling \n",
    "   \n",
    "       If an image is larger than the display, only a portion can be \n",
    "seen at once.  Scrolling is the technique of moving the image in the \n",
    "display so that different parts are visible.  \n",
    "     \n",
    "   \n",
    "  Zooming \n",
    "   \n",
    "     Zooming an image on a display gives the visual impression of \n",
    "getting closer or further from the objects.  In this program, zooming \n",
    "in is done by copying the pixels and zooming out by dropping pixels. \n",
    "This technique either blows up a portion of the image for easier \n",
    "examination or increases the region of the image that can be shown at \n",
    "once on the display.   \n",
    "    \n",
    "*** FINISHED ***" }; /* end of text */
  
  /* fill in arrays */
  
  number_topics = 0;
  
  /* help on help display */
  topic_title[number_topics] = "How to use help";
  topic_text[number_topics++] = helponhelp_text;
  
  /* intro to help display */
  topic_title[number_topics] = "Introduction";
  topic_text[number_topics++] = intro_text;
  
  /* overview of ObitView */
  topic_title[number_topics] = "Overview";
  topic_text[number_topics++] =  over_text;
  
  /* Browser  */
  topic_title[number_topics] = "Browser";
  topic_text[number_topics++] = browser_text;
  
  /* file menu */
  topic_title[number_topics] = "File menu";
  topic_text[number_topics++] = file_text;
  topic_title[number_topics] = "File/Open FITS";
  topic_text[number_topics++] = fileopen_text;
  topic_title[number_topics] = "File/Open AIPS";
  topic_text[number_topics++] = afileopen_text;
  topic_title[number_topics] = "File/Preview";
  topic_text[number_topics++] = filepreview_text;
  topic_title[number_topics] = "File/Save as";
  topic_text[number_topics++] = filesaveas_text;
  topic_title[number_topics] = "File/Image info";
  topic_text[number_topics++] = fileinfo_text;
  topic_title[number_topics] = "File/Log position";
  topic_text[number_topics++] = filelog_text;
  topic_title[number_topics] = "File/Quit";
  topic_text[number_topics++] = filequit_text;
  topic_title[number_topics] = "File/Save Window as";
  topic_text[number_topics++] = filesaveimwindow_text;
  topic_title[number_topics] = "File/Help";
  topic_text[number_topics++] = filehelp_text;
  
  /* Options menu */
  topic_title[number_topics] = "Options menu";
  topic_text[number_topics++] = options_text;
  topic_title[number_topics] = "Options/Help";
  topic_text[number_topics++] = optionshelp_text;
  
  /* Zoom menu */
  topic_title[number_topics] = "Zoom menu";
  topic_text[number_topics++] = zoom_text;
  topic_title[number_topics] = "Zoom/Help";
  topic_text[number_topics++] = zoomhelp_text;
  
  /* Position menu */
  topic_title[number_topics] = "Position menu";
  topic_text[number_topics++] = position_text;
  topic_title[number_topics] = "Position/Set equinox";
  topic_text[number_topics++] = positionsetequinox_text;
  topic_title[number_topics] = "Position/Mark position";
  topic_text[number_topics++] = positionmark_text;
  topic_title[number_topics] = "Position/Lookup position";
  topic_text[number_topics++] = positionlookup_text;
  topic_title[number_topics] = "Position/Fit position";
  topic_text[number_topics++] = positionfit_text;
  topic_title[number_topics] = "Position/Help";
  topic_text[number_topics++] = positionhelp_text;
  
  /* Blink menu */
  topic_title[number_topics] = "Blink menu";
  topic_text[number_topics++] = blink_text;
  topic_title[number_topics] = "Blink/Swap";
  topic_text[number_topics++] = blinkswap_text;
  topic_title[number_topics] = "Blink/Blink";
  topic_text[number_topics++] = blinkblink_text;
  topic_title[number_topics] = "Blink/Help";
  topic_text[number_topics++] = blinkhelp_text;
  
  /* Movie menu */
  topic_title[number_topics] = "Movie menu";
  topic_text[number_topics++] = movie_text;
  topic_title[number_topics] = "Movie/Help";
  topic_text[number_topics++] = moviehelp_text;
  
  /* Colorize menu */
  topic_title[number_topics] = "Colorize menu";
  topic_text[number_topics++] = colorize_text;
  topic_title[number_topics] = "Colorize/Color Contour";
  topic_text[number_topics++] = colorizecolor_text;
  topic_title[number_topics] = "Colorize/Pseudo flame";
  topic_text[number_topics++] = colorizeflame_text;
  topic_title[number_topics] = "Colorize/Grayscale";
  topic_text[number_topics++] = colorizegray_text;
  topic_title[number_topics] = "Colorize/Reverse colors";
  topic_text[number_topics++] = colorizereverse_text;
  topic_title[number_topics] = "Colorize/Reset colors";
  topic_text[number_topics++] = colorizereset_text;
  topic_title[number_topics] = "Colorize/Clear Graphics";
  topic_text[number_topics++] = cleargraphics_text;
  topic_title[number_topics] = "Colorize/Set Graphics Color";
  topic_text[number_topics++] = setgraphicscolor_text;
  topic_title[number_topics] = "Colorize/Help";
  topic_text[number_topics++] = colorizehelp_text;
  
  /* Help menu  */
  topic_title[number_topics] = "Help menu";
  topic_text[number_topics++] = help_text;

  /* Window editing */
  topic_title[number_topics] = "Window Editing";
  topic_text[number_topics++] = windowedit_text;
  
  /* xmlrpc interface */
  topic_title[number_topics] = "xmlrpc interface";
  topic_text[number_topics++] = xmlrpcinterface_text;
  
  /* Display control  */
  topic_title[number_topics] = "Display control";
  topic_text[number_topics++] = displaycontrol_text;
    
  /* Graphics overlay  */
  topic_title[number_topics] = "Graphics overlay";
  topic_text[number_topics++] = graphicsoverlay_text;
    
  /* Image position  */
  topic_title[number_topics] = "Image position";
  topic_text[number_topics++] = imageposition_text;
  
  /* Image scrolling  */
  topic_title[number_topics] = "Image scrolling";
  topic_text[number_topics++] = imagescrolling_text;
  
  /* Image enhancement  */
  topic_title[number_topics] = "Image enhancement";
  topic_text[number_topics++] = imageenhancement_text;
  
  /*  Glossary  */
  topic_title[number_topics] = "Glossary";
  topic_text[number_topics++] = glossary_text;
  
  
} /* end InitHelpText*/
