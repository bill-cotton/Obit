/* $Id: ObitParser.c,v 1.8 2006/06/28 17:20:20 bcotton Exp $ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2008                                          */
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
#include "ObitParser.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitParser.c
 * ObitParser function definitions.
 *
 * This Files contains text file parsing utilities
 */

/*---------------Private function prototypes----------------*/

/** Private: parse next entry */
static ObitIOCode ObitParserEntry(ObitFile *myFile,  gchar* name, ObitInfoType *type, 
				  gint32 *dim, gpointer *data, ObitErr *err);


/*----------------------Public functions---------------------------*/
/**
 * Parse text file
 * \param infile Name of the input text file to parse
 * \param list  ObitInfoList to accept values.
 * \param err ObitErr for reporting errors.
 * \return return code, OBIT_IO_OK => OK
 */
ObitIOCode ObitParserParse(gchar *infile, ObitInfoList *list, ObitErr *err)
{
  ObitFile *myFile=NULL;
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitInfoType type=0;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar name[100];
  gpointer data=NULL;
  olong i;
  gchar *routine = "ObitParserParse";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert(ObitInfoListIsA(list));
  g_assert(infile!=NULL);

  /* Open text file */
  myFile =  newObitFile(infile);
  retCode = ObitFileOpen (myFile, infile, OBIT_IO_ReadOnly, OBIT_IO_Text, 0, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) 
    Obit_traceback_val (err, routine, infile, retCode);

  /* Loop through file parsing */
  while (retCode==OBIT_IO_OK) {
    retCode = ObitParserEntry(myFile, name, &type, dim, &data, err);
    if (retCode!=OBIT_IO_OK) break;

    /* save to InfoList */
    ObitInfoListAlwaysPut (list, name, type, dim, (gconstpointer)data);

    /* May have array of strings */
    if ((type==OBIT_string) && (dim[1]>1)) {
      for (i=0; i<dim[1]; i++) g_free(((gchar**)data)[i]);
      g_free(data); data = NULL; /* release storage */
    } else {
      g_free(data); data = NULL;  /*release storage */
    }

  } /* end loop over file */

  /* Close */
  retCode = ObitFileClose (myFile, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) 
    Obit_traceback_val (err, routine, infile, retCode);

  /* Cleanup */
  myFile = ObitFileUnref(myFile);
  if (data) g_free(data); data = NULL; /* release storage */
 
  return retCode;
} /* end ObitParserParse */

/**
 * Parse text entry into parts
 * \param myFile Open ObitFile to read from
 * \param name    [out] The label (keyword) of the information.
 * \param type    [out] Data type of data element (enum defined in ObitInfoList class.
 * \param dim     [out] Dimensionality of datum.
 *                Note: for strings, the first element is the length in char.
 * \param data Pointer to the data.  The memory should be g_freed after use.
 */
static ObitIOCode ObitParserEntry(ObitFile *myFile, gchar* name, ObitInfoType *type, 
				  gint32 *dim, gpointer *data, ObitErr *err)
{
  gchar line[2049], typeStr[20], temp[60];
  olong size, i, j, i1, i2, iStart, iEnd, nread, idim, iout,  nvalue;
  ObitIOCode retCode = OBIT_IO_SpecErr;
  gboolean done, areMore;
  gchar *cdata=NULL, **carray=NULL;
  ofloat ftemp, *fdata;
  odouble dtemp, *ddata;
  oint *idata;
  olong itemp;
  gboolean *bdata;
  gchar *routine = "ObitParserEntry";

 /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitFileIsA(myFile));
  g_assert(name!=NULL);
  g_assert(type!=NULL);
  g_assert(dim!=NULL);
  g_assert(data!=NULL);

  /* Read header line - look for next as there may be comments */
  done = FALSE;
  while (!done) {
    retCode = ObitFileReadLine (myFile, line, 2048, err);
    if (err->error) Obit_traceback_val (err, routine, "Input parser", retCode);
    if (retCode==OBIT_IO_EOF) return retCode; /* done? */
    done = !strncmp (line, "$Key = ", 7);
  }

  /* Parse header line */
  iStart = 7; 
  iEnd = strlen(line);

  /* name */
  i1 = -1;
  i2 = -1;
  for (i=iStart; i<iEnd; i++) {
    if ((line[i]!=' ') && (i1<0)) i1 = i; /* beginning of name? */
    if ((i1>0) && (line[i]==' ')) { /* end of name */
      i2 = i;
      iStart = i+1;
      break;
    }
  }
  j = 0;
  for (i=i1; i<i2; i++) name[j++] = line[i]; name[j] = 0;

  /* Type code */
  i1 = -1;
  for (i=iStart; i<iEnd; i++) {
    if ((line[i]!=' ') && (i1<0)) i1 = i; /* beginning of type? */
    if ((i1>0) && (line[i]==' ')) { /* end? */
      i2 = i;
      iStart = i+1;
      break;
    }
  }
  j = 0;
  for (i=i1; i<i2; i++) typeStr[j++] = line[i]; typeStr[j] = 0;

  /* Interprete type */
  if (!strncmp (typeStr, "Str", 3)) {
    *type = OBIT_string;
  } else if (!strncmp (typeStr, "Int", 3)) {
    *type = OBIT_oint;
  } else if (!strncmp (typeStr, "Boo", 3)) {
    *type = OBIT_bool;
  } else if (!strncmp (typeStr, "Flt", 3)) {
    *type = OBIT_float;
  } else if (!strncmp (typeStr, "Dbl", 3)) {
    *type = OBIT_double;
  } else { /* unknown */
    Obit_log_error(err, OBIT_Error, "%s: Unknown data type %s", routine, typeStr);
    Obit_log_error(err, OBIT_Error, "%s: Error parsing value in line %s", 
		   routine, line);
    return OBIT_IO_ReadErr;
  }

  /* Get dimensionality */
  for (i=0; i<MAXINFOELEMDIM; i++) dim[i] = 1; /* init */
  i1 = -1;
  idim = 0;
  for (i=iStart; i<iEnd; i++) {
    if ((line[i]=='(') && (i1<0)) i1 = i+1; /* beginning of dimension? */
    if ((i1>0) && ((line[i]==',') || (line[i]==')'))) { /* end of dimension */
      i2 = i;
      j = 0;
      for (i=i1; i<i2; i++) temp[j++] = line[i]; temp[j] = 0;
      nread = sscanf (temp, "%d", &itemp);
      dim[idim++] = itemp;
      if (nread!=1) {
	Obit_log_error(err, OBIT_Error, "%s: Error parsing dimension in line %s", 
		       routine, line);
	return OBIT_IO_ReadErr;
     }
      i1 = i+1;
    } /* end of parse dimension value */
    if (line[i]==')') { /* end of dimensions */
      break;
    }
  } /* End parsing dimension array */

  /* allocate output */
  size = ObitInfoElemSize (*type, dim);
  *data = g_malloc0(size+16);

  /* How many values to parse? */
  nvalue = 1;
  for (i=0; i<MAXINFOELEMDIM; i++) nvalue *= MAX (1, dim[i]);

  /* by data type parse data */
  done = FALSE;
  iout = 0;
  switch (*type) {
  case OBIT_string:

    /* String array? */
    if (dim[1]>1) {
      g_free(*data);
      *data = g_malloc(sizeof(gchar*)*dim[1]);
      carray = (gchar**)*data;
    } else { /* single string */
      cdata = (gchar*)*data;
    }
    /* read whole strings */
    nvalue /= dim[0];
    while (!done) {
      /* read line */
      retCode = ObitFileReadLine (myFile, line, 2048, err);
      if (err->error) Obit_traceback_val (err, routine, "Input parser", retCode);
      if (retCode==OBIT_IO_EOF) return retCode; /* done? */

      /* terminate line at comment delimiter */
      for (i=0; i<strlen(line); i++) if (line[i]=='#') line[i] = 0;

      /* terminate line at line feed */
      for (i=0; i<strlen(line); i++) if (line[i]=='\n') line[i] = 0;

      /* find first non blank */
      i1 = -1;
      for (i=0; i<strlen(line); i++) {
	if (line[i]!=' ') {
	  i1 = i;
	  break;
	}
      }

      /* If only one string use actual size in dimension
	 Oh SHIT - this can blow the allocation
      if ((nvalue==1) && (iout==0)) dim[0] = strlen(line)-i1; */

      /* Copy to end of line */
      j = 0;
      i1 = MAX (0, i1);
      if (dim[1]>1) { /* string array */
	carray[iout] = g_malloc0(sizeof(gchar)*(dim[0]+1));
	for (i=0; i<dim[0]; i++)  carray[iout][i] = ' '; carray[iout][i] = 0;
	for (i=i1; i<MIN(dim[0],strlen(line)); i++)  carray[iout][j++] = line[i];
	iout++;
      } else { /* Single string */
	/* blank fill */
	for (i=0; i<dim[0]; i++) cdata[iout+i] = ' ';
	for (i=i1; i<MIN(dim[0],strlen(line)); i++)  cdata[iout+(j++)] = line[i];
	iout += dim[0]; /* index in output */
      }

      /* Keep track of reads */
      nvalue--;
      done = nvalue<=0;  /* read them all? */
    } /* end loop reading strings */
    break;

  case OBIT_oint:
    idata = (oint*)*data;
    while (!done) {
      /* read line */
      retCode = ObitFileReadLine (myFile, line, 2048, err);
      if (err->error) Obit_traceback_val (err, routine, "Input parser", retCode);
      if (retCode==OBIT_IO_EOF) return retCode; /* done? */

      /* terminate line at comment delimiter */
      for (i=0; i<strlen(line); i++) if (line[i]=='#') line[i] = 0;

      /* Loop over possible entries in line */
      iStart = 0; 
      iEnd = strlen(line);
      areMore = TRUE;
      while (areMore) {
	i1 = -1;
	/* find beginning and end */
	for (i=iStart; i<iEnd; i++) {
	  if (((line[i]!=' ') && (line[i]!='\n')) && (i1<0)) i1 = i; /* beginning? */
	  if ((i1>=0) && ((line[i]==' ') || (line[i]=='\n'))) { /* end?  */
	    i2 = i;
	    iStart = i+1;
	    break;
	  }
	}

	/* Find one? */
	if (i1<0) break;

	/* copy to temp */
	j = 0;
	for (i=i1; i<i2; i++) temp[j++] = line[i]; temp[j] = 0;
	/* get value */
	nread = sscanf (temp, "%d", &itemp);
	idata[iout++] = itemp;
	if (nread!=1) {
	  Obit_log_error(err, OBIT_Error, "%s: Error parsing value in line %s", 
			 routine, line);
	  return OBIT_IO_ReadErr;
	}
	/* Keep track of reads */
	nvalue--;
	if (nvalue<=0) break;
      } /* end loop over line */

      done = nvalue<=0;  /* read them all? */
    } /* end loop reading integers */
    break;

  case OBIT_bool:
    bdata = (gboolean*)*data;
    while (!done) {
      /* read line */
      retCode = ObitFileReadLine (myFile, line, 2048, err);
      if (err->error) Obit_traceback_val (err, routine, "Input parser", retCode);
      if (retCode==OBIT_IO_EOF) return retCode; /* done? */

      /* terminate line at comment delimiter */
      for (i=0; i<strlen(line); i++) if (line[i]=='#') line[i] = 0;

       /* Loop over possible entries in line */
      iStart = 0; 
      iEnd = strlen(line);
      for (i=iStart; i<iEnd; i++) {
	if (line[i]=='T') {
	  bdata[iout++] = TRUE;
	  nvalue--;
	}
	if (line[i]=='F') {
	  bdata[iout++] = FALSE;
	  nvalue--;
	}
	if (nvalue<=0) break;
      }
      
     /* Keep track of reads */
      done = nvalue<=0;  /* read them all? */
    } /* end loop reading Booleans */
    break;

  case OBIT_float:
    fdata = (ofloat*)*data;
    while (!done) {
      /* read line */
      retCode = ObitFileReadLine (myFile, line, 2048, err);
      if (err->error) Obit_traceback_val (err, routine, "Input parser", retCode);
      if (retCode==OBIT_IO_EOF) return retCode; /* done? */
      
      /* terminate line at comment delimiter */
      for (i=0; i<strlen(line); i++) if (line[i]=='#') line[i] = 0;
      
      /* Loop over possible entries in line */
      iStart = 0; 
      iEnd = strlen(line);
      areMore = TRUE;
      while (areMore) {
	i1 = -1;
	/* find beginning and end */
	for (i=iStart; i<iEnd; i++) {
	  if (((line[i]!=' ') && (line[i]!='\n')) && (i1<0)) i1 = i; /* beginning? */
	  if ((i1>=0) && ((line[i]==' ') || (line[i]=='\n'))) { /* end?  */
	    i2 = i;
	    iStart = i+1;
	    break;
	  }
	}
	
	/* Find one? */
	if (i1<0) break;
	
	/* copy to temp */
	j = 0;
	for (i=i1; i<i2; i++) temp[j++] = line[i]; temp[j] = 0;
	
	/* get value */
	nread = sscanf (temp, "%f", &ftemp);
	fdata[iout++] = ftemp;
	if (nread!=1) {
	  Obit_log_error(err, OBIT_Error, "%s: Error parsing value in line %s", 
			 routine, line);
	  return OBIT_IO_ReadErr;
	}
	/* Keep track of reads */
	nvalue--;
	if (nvalue<=0) break;
      } /* end loop over line */
      
      done = nvalue<=0;  /* read them all? */
    } /* end loop reading floatss */
    break;
   case OBIT_double:
    ddata = (odouble*)*data;
    while (!done) {
      /* read line */
      retCode = ObitFileReadLine (myFile, line, 2048, err);
      if (err->error) Obit_traceback_val (err, routine, "Input parser", retCode);
      if (retCode==OBIT_IO_EOF) return retCode; /* done? */
      
      /* terminate line at comment delimiter */
      for (i=0; i<strlen(line); i++) if (line[i]=='#') line[i] = 0;
      
      /* Loop over possible entries in line */
      iStart = 0; 
      iEnd = strlen(line);
      areMore = TRUE;
      while (areMore) {
	i1 = -1;
	/* find beginning and end */
	for (i=iStart; i<iEnd; i++) {
	  if (((line[i]!=' ') && (line[i]!='\n')) && (i1<0)) i1 = i; /* beginning? */
	  if ((i1>=0) && ((line[i]==' ') || (line[i]=='\n'))) { /* end?  */
	    i2 = i;
	    iStart = i+1;
	    break;
	  }
	}
	
	/* Find one? */
	if (i1<0) break;
	
	/* copy to temp */
	j = 0;
	for (i=i1; i<i2; i++) temp[j++] = line[i]; temp[j] = 0;
	
	/* get value */
	nread = sscanf (temp, "%lf", &dtemp);
	ddata[iout++] = dtemp;
	if (nread!=1) {
	  Obit_log_error(err, OBIT_Error, "%s: Error parsing value in line %s", 
			 routine, line);
	  return OBIT_IO_ReadErr;
	}
	/* Keep track of reads */
	nvalue--;
	if (nvalue<=0) break;
      } /* end loop over line */
      
      done = nvalue<=0;  /* read them all? */
    } /* end loop reading doubles */
    break;
 default:
    break;
  }; /* end switch by type */
  
  return retCode;
} /* end ObitParserEntry */
