/* $Id: helpbox.h,v 1.1.1.1 2005/06/09 12:45:29 bcotton Exp $ */
/* help dialog box header for XFITSview */
/*-----------------------------------------------------------------------
*  Copyright (C) 1998
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
#ifndef HELPBOX_H
#define HELPBOX_H
/* create HelpBox */
/* note: visible must be true if the dialog is to be displayed 
   immediately, otherwise it will be created bu invisible */
void MakeHelpBox (Widget parent, Boolean visible); 

/* Help callback on a given topic  */
/* MakeHelpBox must have been called first */
/* clientData is pointer to topic string */
void HelpBoxTopicCB (Widget w, XtPointer clientData, XtPointer callData);

/* Displays help on a given topic  */
/* MakeHelpBox must have been called first */
void HelpBoxShowTopic (char* topic);

/* Show Help dialog */
/* MakeHelpBox must have been called first */
void HelpBoxShow(void);

/* Put help text in display */
/* MakeHelpBox must have been called first */
void HelpBoxSetText (char** text);

/* adds topic to end of list */
/* MakeHelpBox must have been called first */
void HelpBoxAddTopic (char* topic);

/* delete all topics */
/* MakeHelpBox must have been called first */
void HelpBoxDeleteAllTopics(void);

#endif /* end HELPBOX_H */


