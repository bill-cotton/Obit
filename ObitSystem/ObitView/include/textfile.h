/* $Id: textfile.h,v 1.2 2005/07/23 00:20:08 bcotton Exp $ */
/*  header file for TextFile (Unix/Xwindow) utilities */ 
/*-----------------------------------------------------------------------
*  Copyright (C) 1996-2008
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
  
#include <stdio.h>
#include "zlib.h" 

#ifndef TEXTFILE_H 
#define TEXTFILE_H 
#define TEXTBUFFER_SIZE 4096  /* I/O buffer size */
 
/* argument is actually TextFilePtr */
typedef void (*TextFileProc)(XPointer);

typedef struct {
Widget w;                  /* Widget to attach Dialog boxes to */
long InBuffer;             /* The number of characters in Buffer */ 
long BufferPos;            /* Next character position in Buffer */
TextFileProc OKfunc;       /* function to be called to process TextFile */
TextFileProc CancelFunc;   /* function to be called if select file canceled
			      or fails*/
char *FileName;            /* file name */
char *directory;           /* file directory */
unsigned char Buffer[TEXTBUFFER_SIZE];/* I/O buffer */
FILE *file;                /* file pointer for write */
gzFile *zFile;             /* zlib file for reading */
int  Good;                 /* file properly selected */
int  HitEOF;               /* if 1 then read an EOF else 0 */
int  State;                /* State=0=> not open, 1=> read, 2=write */
int  FileType;             /* File type 0=unknown, 1=text, 2=FITS, 3=DOStext*/
} TextFileInfo, *TextFilePtr;
 
/* Public file functions */ 

/* create/initialize TextFileInfo structure                               */
/* w is a Widget for dialog boxes to be attached to                       */
/* if non null filename and directory filled in                           */
TextFilePtr TextFileMake (Widget w, char* filename, char* directory);

/* delete TextFileInfo structure                                          */
/* if non null filename and directory filled in                           */
void TextFileKill (TextFilePtr TFilePtr);

/* get file specification, fill in TextFileInfo, call specified function  */
/* inout = 1 => input; 2=> output                                         */
/* TFilePtr->Good indicated if file selected func called                  */
void TextFileFind (int inout, TextFilePtr TFilePtr, TextFileProc OKfunc,
		   TextFileProc CancelFunc);
  
/* open text file                                                         */ 
/* inout = 1 => input; 2=> output                                         */
/* returns 1 if successful                                                */
int TextFileOpen (TextFilePtr TFilePtr, int inout);
  
/* close file specified by TFilePtr                                       */ 
/* returns 1 if successful                                                */ 
int TextFileClose(TextFilePtr TFilePtr); 
  
/* Read a line of text from a file                                        */ 
/* returns 1 if OK, 0 if failed, -1 -> HitEOF                             */ 
/* maxchar is the maximum number of characters allowed in line            */
int TextFileRead(TextFilePtr TFilePtr, char* line, int maxchar); 
  
/* Write a line of text to a file                                         */ 
/* returns 1 if OK, 0 if failed                                           */ 
int TextFileWrite(TextFilePtr TFilePtr, char* line); 
  
#endif /* TEXTFILE_H */ 
  
  
  
  
  
