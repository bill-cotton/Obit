/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2002-2008                                          */
/*;  Associated Universities, Inc. Washington DC, USA.                */
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
#include <glib.h>
#include <string.h>
#include "ObitInfoElem.h"
#include "ObitMem.h"

/**
 * \file ObitInfoElem.c
 * Elements of an ObitInfoList function definition.
 * Stores a label, size and type info and a data array.
 * The limit on the number of dimensions is MAXINFOELEMDIM.
 */

/**
 * Constructor
 * \param label Name of element, use as key.
 * \param type Data type of object
 * \param dim  Dimension array (max no. dimensions MAXINFOELEMDIM)
 * \param data Data array, if NULL only create array and not copy data.
 * \return the new object.
 */
ObitInfoElem* newObitInfoElem (gchar *label, ObitInfoType type, 
			       gint32 *dim, gconstpointer data )
{ 
  ObitInfoElem *me; 
  gchar *name;
 
  /* error checks */
  g_assert (label != NULL);
  g_assert (dim   != NULL);

  name =  g_strconcat ("IElem:", label, NULL);
  me =  ObitMemAlloc0Name(sizeof(ObitInfoElem),name);
  g_free(name);
  me->iname = g_strstrip(g_strdup(label));
  me->itype = type; 
  me->data = NULL;
  me->size = 0;

  /* set data size */
  ObitInfoElemResize (me, type, dim);

  /* copy data */
  if (data!=NULL) ObitInfoElemSave (me, data);

  return me; 
} /* end newObitInfoElem */ 
  
/**
 * Copy constructor
 * \param in   Object to copy
 * \return the new object, a byte-for-byte copy.
 */
ObitInfoElem* ObitInfoElemCopy (ObitInfoElem *in)
{ 
  ObitInfoElem *out; 
  olong i;
  gchar *name;
 
  name =  g_strconcat ("IElem:", in->iname, NULL);
  out =  ObitMemAlloc0Name(sizeof(ObitInfoElem), name);
  g_free(name);
  out->iname = g_strdup(in->iname);
  out->itype = in->itype; 
  out->size  = in->size;
  for (i=0; i<MAXINFOELEMDIM; i++) out->idim[i] = in->idim[i];

  /* allocate storage */
  name =  g_strconcat ("IEData:", in->iname, NULL);
  out->data = ObitMemAlloc0Name(out->size, name); /* allocate new */
  g_free(name);
  
  /* Copy data */
  g_memmove (out->data, in->data, in->size);

  return out; 
} /* end ObitInfoElemCopy */ 
  
  
/**
 * Destructor 
 * \param me Object to delete
 */
  void freeObitInfoElem(ObitInfoElem *me) 
{ 
  g_assert (me != NULL);

  /* deallocate */
  if (me->iname) g_free (me->iname);
  if (me->data)  ObitMemFree (me->data);
  ObitMemFree (me);
}/* end of freeObitInfoElem  */ 
  
  
/**
 * Compare element name with test string. 
 * TRUE=match, else no match.
 * \param me Object to compare
 * \param testname String to compare with
 * \return True if they are the same.
 */
 gboolean ObitInfoElemTest (ObitInfoElem *me, char *testname)
{
  g_assert (me != NULL);
  g_assert (testname != NULL);

  /* are the lengths the same? */
  if (strlen(testname)!=strlen(me->iname)) return FALSE;
  /* do the strings match? */
  return (!strncmp(me->iname, testname, strlen(testname)));
} /* end ObitInfoElemTest */ 
  
/**
 * Compare size and type of object with specified values.
 * \param me Object to compare
 * \param type Type for comparison.
 * \param dim Dimension array to compare (max MAXINFOELEMDIM)
 * \return True if they are the same.
 */
gboolean ObitInfoElemComp (ObitInfoElem *me, ObitInfoType type, 
			   gint32 *dim)
{
  gint32 i;

  /* error checks */
  g_assert (me != NULL);
  g_assert (dim != NULL);

  /* check type */
  if (me->itype != type) return FALSE;

  /* compare dimensions */
  for (i=0; i<MAXINFOELEMDIM; i++)
    if (MAX (1,dim[i]) != MAX (1, me->idim[i])) return FALSE;

  return TRUE; /* must be OK */
} /* end ObitInfoElemComp */ 
  
/**
 * Update contents of an info element.
 * \param me Object to update.
 * \param type data type 
 * \param dim Dimension array (max MAXINFOELEMDIM).
 * \param data Data array
 * \param warn if true, issue detailed warning about mismatches
 * \return TRUE if OK, FALSE if object incompatable with request.
 */
gboolean ObitInfoElemUpdate (ObitInfoElem *me, gint32 type, 
			 gint32 *dim, gconstpointer data, gboolean warn)
{ 
  gboolean OK=FALSE;

 /* error checks */
  g_assert (me != NULL);
  g_assert (dim != NULL);
  g_assert (data != NULL);

  /* Check compatability */
  OK = ObitInfoElemComp(me,type,dim);
  if (!OK && warn) {
    g_warning ("ObitInfoElemUpdate: request incompatable with existing entry");
    g_warning (" %s  Old: type = %d dim = %d %d %d", 
	       me->iname, me->itype, me->idim[0], me->idim[1], me->idim[2]);
    g_warning (" %s  New: type = %d dim = %d %d %d", 
	       me->iname,type, dim[0], dim[1], dim[2]);
  }

  if (!OK) return OK;

  /* copy data */
  ObitInfoElemSave (me, data);

  return TRUE;
} /* end ObitInfoElemUpdate */ 
  
/**
 * Save data on object.
 * \param me Object to store data to.
 * \param data The new data array.
 * \return the new object.
 */
void ObitInfoElemSave (ObitInfoElem *me, gconstpointer data)
{ 
  olong i, j, off;
  const gchar *ddata, **c2data, ***c3data;

  /* error checks */
  g_assert (me != NULL);
  g_assert (data != NULL);
   
  /* Strings are (always) different - an extra byte is allocated for final NULL */
  if (me->itype==OBIT_string) {
    /* blank fill */
    for (i=0; i<me->size; i++) ((gchar*)me->data)[i] = ' ';
    /* Final NULL*/
    ((gchar*)me->data)[me->size-1] = 0;
    ddata = data;
    if (me->idim[1]==1) { /* 1-D strings */
      g_memmove (me->data, data, MIN (me->idim[0], strlen(ddata)));
    } else if (me->idim[1]>1) { /* >=2-D strings */
      off = 0;
      if (me->idim[2]==1) { /* no third axis */
	c2data = (const gchar**)data;
	for (i=0; i<me->idim[1]; i++) {
	  g_memmove (&((gchar*)me->data)[off], c2data[i], MIN (me->idim[0], strlen(c2data[i])));
	  off += me->idim[0];
	}
      } else if (me->idim[2]>1) { /* 3D array */
	c3data = (const gchar***)data;
	for (j=0; j<me->idim[2]; j++) {
	  for (i=0; i<me->idim[1]; i++) {
	    g_memmove (&((gchar*)me->data)[off], c3data[j][i], MIN (me->idim[0], strlen(c3data[j][i])));
	    off += me->idim[0];
	  }
	}
      }
    }
  } else { /* anything else (nonstring) */
    g_memmove (me->data, data, me->size);
  }
}  /*  End of  ObitInfoElemSave  */
  
/**
 * Set the size and type of an ObitInfoElem.
 * Any existing values in data will be lost
 * \param me Object to modify.
 * \param type New data type.
 * \param dim  dimension of the array.
 * \return the new object.
 */
void ObitInfoElemResize  (ObitInfoElem *me, ObitInfoType type, 
			  gint32 *dim)
{
  gint32 i, size; 
  gchar *name;
  
  me->itype = type; /* reset type */

  for (i=0; i<MAXINFOELEMDIM; i++)
    me->idim[i] = dim[i]; /* Copy dim.  */ 

  /* get size */
  size = ObitInfoElemSize (type, dim);
  /* Add an extra byte for strings */
  if (type==OBIT_string) size++;
  if (me->data) ObitMemFree(me->data); /* free old if allocated */
  name =  g_strconcat ("IEData:", me->iname, NULL);
  me->data = ObitMemAlloc0Name(size, name); /* allocate new */
  g_free(name);

  /* remember how big this is */
  me->size = size;
}  /*  End of ObitInfoElemResize  */

/**
 * Determine the size an ObitInfoElem data record.
 * \param type Data type.
 * \param dim  Dimension of the array.
 * \return the size in bytes.
 */
olong ObitInfoElemSize  (ObitInfoType type, gint32 *dim)
{
  gint32 i, size, number; 
  
  number = dim[0];
  for (i=1; i<MAXINFOELEMDIM; i++) {
    number *= MAX (1, dim[i]);
  }

  /* determine element size */
  switch (type) { 
    case OBIT_byte:
      size = sizeof(gint8);
      break;
    case OBIT_short:
      size = sizeof(gint16);
      break;
    case OBIT_int:
      size = sizeof(olong);
      break;
     case OBIT_oint:
      size = sizeof(oint);
      break;
    case OBIT_long:
      size = sizeof(olong);
      break;
    case OBIT_llong:
      size = sizeof(ollong);
      break;
    case OBIT_ubyte:
      size = sizeof(guint8);
      break;
    case OBIT_ushort:
      size = sizeof(guint16);
      break;
    case OBIT_uint:
      size = sizeof(guint);
      break;
    case OBIT_ulong:
      size = sizeof(gulong);
      break;
    case OBIT_float:
      size = sizeof(ofloat);
      break;
    case OBIT_double:
      size = sizeof(odouble);
      break;
    case OBIT_complex:
      size = 2*sizeof(ofloat);
      break;
    case OBIT_dcomplex:
      size = 2*sizeof(odouble);
      break;
    case OBIT_string:
      size = sizeof(gchar);
      break;
    case OBIT_bool:
      size = sizeof(gboolean);
      break;
    case OBIT_bits:
      size = sizeof(olong);
      /* Calculate the number of gints involved */
      number = 1 + ((number-1) / sizeof(olong));
      break;
  default:
    size = 16;
    g_assert_not_reached(); /* unknown, barf */
  }; /* end switch to find size */

  return number*size;
}  /*  End of ObitInfoElemSize  */


/**
 * Print an InfoListElem (GHFunc) 
 * This function is called from g_hash_table_foreach.
 * \param key   pointer to hash key
 * \param in    ObitInfoElem
 * \param file  FILE* to write to
 */
void  ObitInfoElemPrint(ObitInfoElem *me, FILE *file)
{
  /* Should match order of enum obitInfoType in ObitTypes.h */
  gchar *infoType[] = {"byte", "short", "int", "oint", "long",  
		       "ubyte", "ushort", "uint", "ulong", "llong",
		       "float", "double", "complex", "dcomplex",
		       "string", "bool", "bits"};
  olong        *ldata, i, j, more, indx, ltemp, lstr, size;
  ollong       *lldata, lltemp;
  gboolean     *bdata;
  olong        *idata;
  oint         *odata;
  ofloat       *fdata;
  odouble      *ddata;
  gchar        *cdata, line[81], bchar, cstring[65];

  /* Header */
  g_snprintf (line, 80, "item='%s' type=%s, dim=[%d,%d,%d], data= ", 
	   me->iname, infoType[me->itype], me->idim[0], me->idim[1], me->idim[2]);
  indx = strlen (line);

  /* How many? */
  size = MAX (1, me->idim[0]);
  for (i=1; i< MAXINFOELEMDIM; i++)size *= MAX (1, me->idim[i]);

  /* Data by type */
  switch (me->itype) { 
  case OBIT_int:
    idata = (olong*)me->data;
    more = size;
    while (more>0) {
      for (j=0; j<8; j++) {
	g_snprintf (&line[indx], 80-indx, "%d ", *idata++);
	indx = strlen (line);
	more--;                    /* finished? */
	if (more<=0) break;
      }
      fprintf (file, "%s\n", line);
      g_snprintf (line, 80, "    ");
      indx = strlen (line);
    }
    break;
  case OBIT_oint:
    odata = (oint*)me->data;
    more = size;
    while (more>0) {
      for (j=0; j<20; j++) {
	ltemp = (olong)(*odata++);
	g_snprintf (&line[indx], 80-indx, " %d ", ltemp);
	indx = strlen (line);
	more--;                    /* finished? */
	if (more<=0) break;
	if (indx>60) break;        /* Line full? */
      }
      fprintf (file, "%s\n", line);
      g_snprintf (line, 80, "    ");
      indx = strlen (line);
    }
    break;
  case OBIT_long:
    ldata = (olong*)me->data;
    more = size;
    while (more>0) {
      for (j=0; j<20; j++) {
	ltemp = (olong)(*ldata++);
	g_snprintf (&line[indx], 80-indx, " %d ", ltemp);
	indx = strlen (line);
	more--;                    /* finished? */
	if (more<=0) break;
	if (indx>60) break;        /* Line full? */
      }
      fprintf (file, "%s\n", line);
      g_snprintf (line, 80, "    ");
      indx = strlen (line);
    }
    break;
    
  case OBIT_llong:
    lldata = (ollong*)me->data;
    more = size;
    while (more>0) {
      for (j=0; j<20; j++) {
	lltemp = (ollong)(*lldata++);
	g_snprintf (&line[indx], 80-indx, " %ld ", (long)lltemp);
	indx = strlen (line);
	more--;                    /* finished? */
	if (more<=0) break;
	if (indx>60) break;        /* Line full? */
      }
      fprintf (file, "%s\n", line);
      g_snprintf (line, 80, "    ");
      indx = strlen (line);
    }
    break;
    
  case OBIT_float:
    fdata = (ofloat*)me->data;
    more = size;
    while (more>0) {
      for (j=0; j<4; j++) {
	g_snprintf (&line[indx], 80-indx, "%15.5g ", *fdata++);
	indx = strlen (line);
	more--;                    /* finished? */
	if (more<=0) break;
	if (indx>55) break;   /* Line full? */
      }
      fprintf (file, "%s\n", line);
      g_snprintf (line, 80, "    ");
      indx = strlen (line);
    }
    break;
    
  case OBIT_double:
    ddata = (odouble*)me->data;
    more = size;
    while (more>0) {
      for (j=0; j<2; j++) {
	g_snprintf (&line[indx], 80-indx, "%25.12lg ", *ddata++);
	indx = strlen (line);
	more--;                    /* finished? */
	if (more<=0) break;
	if (indx>45) break;   /* Line full? */
      }
      fprintf (file, "%s\n", line);
      g_snprintf (line, 80, "    ");
      indx = strlen (line);
    }
    break;

  case OBIT_string:   /* only 64 char of string */
    cdata = (gchar*)me->data;
    lstr = me->idim[0];  /* length of string */
    if ((80-(indx+lstr))<=0) {  /* Will one fit? */
      fprintf (file, "%s\n", line);
      g_snprintf (line, 80, "    ");
      indx = strlen (line);
    }
    more = (size / lstr);
    while (more>0) {
      for (j=0; j<2; j++) {
	strncpy (cstring, cdata, MIN (lstr, 64));
	cstring[MIN (lstr, 64)] = 0;  /* null terminate */
	cdata += lstr;         /* move down string array */
	g_snprintf (&line[indx], 80-indx, "'%s' ", cstring);
	indx = strlen (line);
	more--;                    /* finished? */
	if (more<=0) break;
	if (indx>40) break;   /* Line full? */
      }
      fprintf (file, "%s\n", line);
      g_snprintf (line, 80, "    ");
      indx = strlen (line);
    }
    break;

  case OBIT_bool:
    bdata = (gboolean*)me->data;
    more = size;
    while (more>0) {
      for (j=0; j<30; j++) {
	if (*bdata++) bchar = 'T';
	else bchar = 'F';
	g_snprintf (&line[indx], 80-indx, "%c ", bchar);
	indx = strlen (line);
	more--;                    /* finished? */
	if (more<=0) break;
	if (indx>60) break;   /* Line full? */
      }
      fprintf (file, "%s\n", line);
      g_snprintf (line, 80, "    ");
      indx = strlen (line);
    }
    break;
    
  default:
    g_assert_not_reached(); /* unknown, barf */
  }; /* end switch to copy by type */

  /* Anything left over? */
  if (indx>5) fprintf (file, "%s\n", line);
} /* end ObitInfoElemPrint */
