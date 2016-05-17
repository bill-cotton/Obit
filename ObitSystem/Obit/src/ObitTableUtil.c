/* $Id$  */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2004-2016                                          */
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

#include <math.h>
#include <string.h>
/*#include "glib/gqsort.h"*/
#include "ObitTableUtil.h"
#include "ObitImage.h"
#include "ObitInfoElem.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitTableUtil.c
 * ObitTableUtil class utility function definitions.
 */

/*-------------------- unions -------------------------------------*/
/** Equivalence for Sort key arrays */
  union ObitSortEquiv { 
    olong   itg;
    ofloat  flt;
    odouble dbl;
    gchar   str[8];
    olong   itg2[2];
    ofloat  flt2[2];
    odouble dbl2[2];
  };
/*-------------------- structure -------------------------------------*/
/** Sort Index plus sort key */
typedef struct {
  olong index;
  union ObitSortEquiv key;
} ObitSortStruct;


/*----------------------Private functions---------------------------*/
/** Private: Form sort structure for a table */
static gpointer 
MakeSortStruct (ObitTable *in, olong which[2], gboolean desc,
		olong *size, olong *number, olong *ncomp,
		ObitInfoType *type, ObitErr *err);

/** Private: Form sort structure for a table for 2 float sort */
static gpointer 
MakeSortStruct2f (ObitTable *in, olong which[4], gboolean desc1,
		  gboolean desc2, olong *size, olong *number, olong *ncomp,
		  ObitInfoType *type, ObitErr *err);

/** Private: Sort comparison function for  integers */
static gint CompareInt (gconstpointer in1, gconstpointer in2, 
			gpointer ncomp);

/** Private: Sort abs ascending comparison function for  integers */
static gint CompareAAInt (gconstpointer in1, gconstpointer in2, 
			  gpointer ncomp);

/** Private: Sort abs descending comparison function for  integers */
static gint CompareADInt (gconstpointer in1, gconstpointer in2, 
			  gpointer ncomp);

/** Private: Sort comparison function for floats */
static gint CompareFloat (gconstpointer in1, gconstpointer in2, 
			  gpointer ncomp);

/** Private: Sort abs ascending comparison function for floats */
static gint CompareAAFloat (gconstpointer in1, gconstpointer in2, 
			    gpointer ncomp);

/** Private: Sort abs descending ccomparison function for floats */
static gint CompareADFloat (gconstpointer in1, gconstpointer in2, 
			    gpointer ncomp);

/** Private: Sort comparison function for doubles */
static gint CompareDouble (gconstpointer in1, gconstpointer in2, 
			   gpointer ncomp);

/** Private: Sort abs ascending comparison function for doubles */
static gint CompareAADouble (gconstpointer in1, gconstpointer in2, 
			     gpointer ncomp);

/** Private: Sort abs descending ccomparison function for doubles */
static gint CompareADDouble (gconstpointer in1, gconstpointer in2, 
			     gpointer ncomp);

/** Private: Sort comparison function for string */
static gint CompareString (gconstpointer in1, gconstpointer in2, 
			   gpointer ncomp);

/** Private: reorder table based on Sort structure */
static ObitIOCode 
ReorderTable(ObitTable *in, gpointer base, olong size, olong number, 
	     olong sortCol1, olong sortCol2, ObitErr *err);
/*----------------------Public functions---------------------------*/


/**
 * Sort the rows of a table on a given column
 * If the table is already marked in this order then it is not resorted.
 * \param in      Table to sort
 * \param colName Column label to sort by (for now first cell)
 * \param desc    If true want descending sort, else ascent
 * \param err     ObitErr error stack.
 * \return I/O Code  OBIT_IO_OK = OK.
 */
ObitIOCode ObitTableUtilSort (ObitTable *in, gchar *colName, gboolean desc,
			      ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  olong i, size, number=0, ncomp, which[2], colNo;
  ObitInfoType type;
  gboolean bugOut;
  gpointer SortStruct = NULL;
  gchar *routine = "ObitTableUtilSort";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitTableIsA(in));

  /* Open table */
  retCode = ObitTableOpen (in, OBIT_IO_ReadWrite, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, in->name, retCode);

  /* Which column is desired? */
  colNo = -1;
  size = strlen(colName);
  for (i=0; i<in->myDesc->nfield; i++) {
    if (!strncmp (colName, in->myDesc->FieldName[i], size)) {
      colNo = i;
      break;
    }
  }

  /* If this is already sorted, button up and go home */
  if (desc) {  /* descending sort */
    bugOut = -(colNo) == (in->myDesc->sort[0]-1);
  } else {  /* ascending sort */
    bugOut = (colNo) == (in->myDesc->sort[0]-1);
  }
  if (bugOut) {
    /* Close table */
    retCode = ObitTableClose (in, err);
    if ((retCode != OBIT_IO_OK) || (err->error)) {
      retCode = OBIT_IO_OK;
      goto cleanup;
    }
    return OBIT_IO_OK;
  }

  /* Found it? */
  if (colNo<0) { /* No column by that name */
     Obit_log_error(err, OBIT_Error, 
		   "%s: NO column %s found in %s", 
		    routine, colName, in->name);
     goto cleanup;
 }

  /* build sort structure from table */
  which[0] = colNo+1;  /* 1-rel */
  which[1] = 1;  /* First cell for now */
  SortStruct = MakeSortStruct (in, which, desc, &size, &number, &ncomp, &type, err);
  if ((SortStruct==NULL) || (err->error)) goto cleanup;

  /* Close table */
  retCode = ObitTableClose (in, err);
  if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;
  
  /* Sort keys by key type */
  switch (type) { 
  case OBIT_int:
  case OBIT_long:
    g_qsort_with_data (SortStruct, number, size, CompareInt, &ncomp);
    break;
  case OBIT_float:
    g_qsort_with_data (SortStruct, number, size, CompareFloat, &ncomp);
    break;
  case OBIT_double:
    g_qsort_with_data (SortStruct, number, size, CompareDouble, &ncomp);
    break;
  case OBIT_string:
    g_qsort_with_data (SortStruct, number, size, CompareString, &ncomp);
    break;
  default:
    g_assert_not_reached(); /* unknown, barf */
  }; /* end switch */

  /* Reorder table */
  retCode = ReorderTable (in, SortStruct, size, number, colNo+1, 0, err);
  if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;
  
  /* Cleanup */
  cleanup: if (SortStruct) g_free(SortStruct);
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);

  return OBIT_IO_OK;
} /* end ObitTableUtilSort */


/**
 * Sort the rows of a table on absolute values of a given column
 * If the table is already marked in this order then it is not resorted.
 * \param in      Table to sort
 * \param colName Column label to sort by (for now first cell)
 * \param desc    If true want descending sort, else ascent
 * \param err     ObitErr error stack.
 * \return I/O Code  OBIT_IO_OK = OK.
 */
ObitIOCode ObitTableUtilAbsSort (ObitTable *in, gchar *colName, gboolean desc,
			     ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  olong i, size, number=0, ncomp, which[2], colNo;
  ObitInfoType type;
  gboolean bugOut;
  gpointer SortStruct = NULL;
  gchar *routine = "ObitTableUtilAbsSort";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitTableIsA(in));

  /* Open table */
  retCode = ObitTableOpen (in, OBIT_IO_ReadWrite, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, in->name, retCode);

  /* Which column is desired? */
  colNo = -1;
  size = strlen(colName);
  for (i=0; i<in->myDesc->nfield; i++) {
    if (!strncmp (colName, in->myDesc->FieldName[i], size)) {
      colNo = i;
      break;
    }
  }

  /* If this is already sorted, button up and go home */
  if (desc) {  /* descending sort */
    bugOut = -(colNo) == (in->myDesc->sort[0]-1);
  } else {  /* ascending sort */
    bugOut = (colNo) == (in->myDesc->sort[0]-1);
  }
  if (bugOut) {
    /* Close table */
    retCode = ObitTableClose (in, err);
    if ((retCode != OBIT_IO_OK) || (err->error)) {
      retCode = OBIT_IO_OK;
      goto cleanup;
    }
    return OBIT_IO_OK;
  }

  /* Found it? */
  if (colNo<0) { /* No column by that name */
     Obit_log_error(err, OBIT_Error, 
		   "%s: NO column %s found in %s", 
		    routine, colName, in->name);
     goto cleanup;
 }

  /* build sort structure from table */
  which[0] = colNo+1;  /* 1-rel */
  which[1] = 1;  /* First cell for now */
  /* ascending/descending order done in comparison */
  SortStruct = MakeSortStruct (in, which, FALSE, &size, &number, &ncomp, &type, err);
  if ((SortStruct==NULL) || (err->error)) goto cleanup;

  /* Close table */
  retCode = ObitTableClose (in, err);
  if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;
  
  /* Sort keys by key type */
  switch (type) { 
  case OBIT_int:
  case OBIT_long:
    if (desc) g_qsort_with_data (SortStruct, number, size, CompareADInt, &ncomp);
    else g_qsort_with_data (SortStruct, number, size, CompareAAInt, &ncomp);
    break;
  case OBIT_float:
    if (desc) g_qsort_with_data (SortStruct, number, size, CompareADFloat, &ncomp);
    else g_qsort_with_data (SortStruct, number, size, CompareAAFloat, &ncomp);
    break;
  case OBIT_double:
    if (desc) g_qsort_with_data (SortStruct, number, size, CompareADDouble, &ncomp);
    else g_qsort_with_data (SortStruct, number, size, CompareAADouble, &ncomp);
    break;
  case OBIT_string:
    g_qsort_with_data (SortStruct, number, size, CompareString, &ncomp);
    break;
  default:
    g_assert_not_reached(); /* unknown, barf */
  }; /* end switch */

  /* Reorder table */
  retCode = ReorderTable (in, SortStruct, size, number, colNo+1, 0, err);
  if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;
  
  /* Cleanup */
  cleanup: if (SortStruct) g_free(SortStruct);
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);

  return OBIT_IO_OK;
} /* end ObitTableUtilAbsSort */


/**
 * Two key sort the rows of a table on two columns as floats.
 * If the table is already marked in this order then it is not resorted.
 * \param in       Table to sort
 * \param colName1 Column label to sort by (for now first cell)
 *                 Most slowly varying key
 * \param cell1    Cell in colName1 of sortkey, 1 rel
 * \param desc1    If true want descending sort, else ascent
 * \param err      ObitErr error stack.
 * \param colName2 Column label for secondary key
 *                 Most rapidly varying key
 * \param cell2    Cell in colName2 of sortkey, 1 rel
 * \param desc2    If true want descending sort key 2, else ascent
 * \return I/O Code  OBIT_IO_OK = OK.
 */
ObitIOCode ObitTableUtilSort2f (ObitTable *in, 
				gchar *colName1, olong cell1, gboolean desc1,
				gchar *colName2, olong cell2, gboolean desc2, 
				ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  olong i, size, number=0, ncomp, which[4], col1No, col2No;
  ObitInfoType type;
  gboolean bugOut;
  gpointer SortStruct = NULL;
  gchar *routine = "ObitTableUtilSort2f";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;
  g_assert (ObitTableIsA(in));

  /* Open table */
  retCode = ObitTableOpen (in, OBIT_IO_ReadWrite, err);
  if ((retCode != OBIT_IO_OK) || (err->error))
    Obit_traceback_val (err, routine, in->name, retCode);

  /* Which column is desired for key1? */
  col1No = -1000;
  size = strlen(colName1);
  for (i=0; i<in->myDesc->nfield; i++) {
    if (!strncmp (colName1, in->myDesc->FieldName[i], size)) {
      col1No = i;
      break;
    }
  }

   /* Which column is desired for key2? */
  col2No = -1000;
  size = strlen(colName2);
  for (i=0; i<in->myDesc->nfield; i++) {
    if (!strncmp (colName2, in->myDesc->FieldName[i], size)) {
      col2No = i;
      break;
    }
  }

 /* If this is already sorted, button up and go home */
  if (desc1) {  /* descending sort key 1*/
    bugOut = -(col1No) == (in->myDesc->sort[0]-1);
  } else {  /* ascending sort */
    bugOut = (col1No) == (in->myDesc->sort[0]-1);
  }
  if (desc2) {  /* descending sort key 2 */
    bugOut = bugOut && (-(col2No) == (in->myDesc->sort[1]-1));
  } else {  /* ascending sort */
    bugOut =  bugOut && ((col2No) == (in->myDesc->sort[1]-1));
  }
  if (bugOut) {
    /* Close table */
    retCode = ObitTableClose (in, err);
    if ((retCode != OBIT_IO_OK) || (err->error)) {
      retCode = OBIT_IO_OK;
      return retCode;
    }
    return OBIT_IO_OK;
  }

  /* Found it? */
  if (col1No<0) { /* No column by that name */
     Obit_log_error(err, OBIT_Error, 
		   "%s: NO column %s found in %s", 
		    routine, colName1, in->name);
     return retCode;
 }

  if (col2No<0) { /* No column by that name */
     Obit_log_error(err, OBIT_Error, 
		   "%s: NO column %s found in %s", 
		    routine, colName2, in->name);
     return retCode;
 }

  /* Check that cell numbers in range */
  Obit_retval_if_fail(((cell2>=1) && (cell1<= in->myDesc->repeat[col1No])),
		      err, retCode, "%s: Cell1 ou2 of range [1, %d]", 
		      routine, in->myDesc->repeat[col1No]);
  Obit_retval_if_fail(((cell2>=1) && (cell2<= in->myDesc->repeat[col1No])),
		      err, retCode, "%s: Cell2 out of range [1, %d]", 
		      routine, in->myDesc->repeat[col1No]);

  /* build sort structure from table */
  which[0] = col1No+1;  /* 1-rel */
  which[1] = cell1;
  which[2] = col2No+1;  /* 1-rel */
  which[3] = cell2; 
  SortStruct = MakeSortStruct2f (in, which, desc1, desc2, &size, &number, &ncomp, 
				 &type, err);
  if ((SortStruct==NULL) || (err->error)) goto cleanup;

  /* Close table */
  retCode = ObitTableClose (in, err);
  if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;
  
  /* Sort only floats */
  g_qsort_with_data (SortStruct, number, size, CompareFloat, &ncomp);

  /* Reorder table */
  retCode = ReorderTable (in, SortStruct, size, number, col1No+1, col2No+1, err);
  if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;
  
  /* Cleanup */
  cleanup: if (SortStruct) g_free(SortStruct);
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);

  return OBIT_IO_OK;
} /* end ObitTableUtilSort2f */


/**
 * Truncate the size of a table to a given number of rows.
 * \param in     Pointer to object
 * \param nrows  Number of rows desired
 * \param err    ObitErr for reporting errors.
 * \return error code, OBIT_IO_OK=> OK
 */
void ObitTableUtilTruncate (ObitTable *in, olong nrows, ObitErr *err)
{
  gchar *routine = "ObitTableUtilTruncate";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;

  in->bufferSize = -1;  /* Don't need to assign buffer here */

  /* Open and close */
  ObitTableOpen(in, OBIT_IO_ReadWrite, err);
  if (err->error)Obit_traceback_msg (err, routine, in->name);

  /* reset count */
  in->myDesc->nrow = MIN (nrows, in->myDesc->nrow);
  /* The one that counts is in the IO */
  ((ObitTableDesc*)(in->myIO->myDesc))->nrow = in->myDesc->nrow;
  /* Mark as changed */
  in->myStatus = OBIT_Modified;
  
  ObitTableClose(in, err);
  if (err->error)Obit_traceback_msg (err, routine, in->name);
  in->bufferSize = 0;  /* May need buffer later */
} /* end ObitTableUtilTruncate */

/*----------------------Private functions---------------------------*/
/**
 * Create/fill sort structure for a table
 * The sort structure has one "entry" per row (the first  dimension of 
 * a 2D array) which contains 1) a olong row number as an index (1-rel), 
 * and 2 one or more entries of the relevant type. 
 * Each valid row in the table has an entry.
 * \param in     Table to sort, assumed already open;
 * \param which  Column and element number to use as key (1-rel)
 * \param desc   If true want descending sort, else ascent
 * \param size   [out] number of bytes in entry
 * \param number [out] number of entries
 * \param type   Data type of key
 * \param ncomp Number of values to compare
 * \return sort structure, use ObitSortStruc to access
 *  should be g_freeed when done.
 */
static gpointer 
MakeSortStruct (ObitTable *in, olong which[2], gboolean desc,
	        olong *size, olong *number, olong *ncomp,
	        ObitInfoType *type, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  gpointer out = NULL;
  ObitTableRow *row;
  ObitSortStruct *entry;
  olong irow, nrow, ttsize, col, cell, byteOffset, i;
  ollong tsize, count;
  ObitInfoType itype;
  gint32 dim[MAXINFOELEMDIM];
  /* Pointers for row data */
  gint8   *ptrint8;
  gint16  *ptrint16;
  olong   *ptrint;
  oint    *ptroint;
  olong   *ptrlong=NULL;
  guint8  *ptruint8;
  guint16 *ptruint16;
  guint   *ptruint;
  gulong  *ptrulong;
  ofloat  *ptrfloat;
  odouble *ptrdouble;
  gchar   *ptrchar, *optr;
  gchar *routine = "MakeSortStruc";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return out;
  g_assert (ObitTableIsA(in));

  /* Get table info */
  nrow = in->myDesc->nrow;

  /* which column? Cell? */
  col  = which[0]-1;
  cell = which[1]-1;

  /* Column type */
  itype = in->myDesc->type[col];
  *type = itype;
  /* Only one integer type for sort */
  if (*type==OBIT_short)       *type = OBIT_long;
  else if (*type==OBIT_oint)   *type = OBIT_long;
  else if (*type==OBIT_long)   *type = OBIT_long;
  else if (*type==OBIT_byte)   *type = OBIT_long;
  else if (*type==OBIT_ubyte)  *type = OBIT_long;
  else if (*type==OBIT_ushort) *type = OBIT_long;
  else if (*type==OBIT_uint)   *type = OBIT_long;
  else if (*type==OBIT_ulong)  *type = OBIT_long;

  /* element size */
  dim[0] = 1; dim[1] = 1;dim[2] = 1;dim[3] = 1;dim[4] = 1;
  /* string size - one element of everything else 
   64 bit OS does something wierd here */
  ttsize = sizeof(ObitSortStruct);
  if (*type == OBIT_string) {
    dim[0] = MIN(8,in->myDesc->dim[col][0]);
    ttsize += dim[0];
  }
  *size = MAX ((sizeof(olong) + ObitInfoElemSize(*type, dim)), ttsize);


  /* Total size of structure in case all rows valid */
  tsize = (*size);
  tsize *= (nrow+10);

  /* If 32 bit check range */
  Obit_retval_if_fail(((sizeof(olong*)>4) || (tsize<1000000000)), err, out, 
		      "%s: Too many records to sort on 32 bit system", routine);
  out = g_malloc(tsize);   /* create output structure */
  
  /* Compare 1 except for strings */
  *ncomp = 1;
  if (*type == OBIT_string) *ncomp = dim[0];

  /* Create table row */
  row = newObitTableRow (in);

  byteOffset = in->myDesc->byteOffset[col];  /* where it starts in row data */

 /* loop over table */
  irow = 0;
  count = 0;
  retCode = OBIT_IO_OK;
  while (retCode==OBIT_IO_OK) {
    irow++;
   retCode = ObitTableReadRow (in, irow, row, err);
    if (retCode == OBIT_IO_EOF) break;
    if ((retCode != OBIT_IO_OK) || (err->error)) 
      Obit_traceback_val (err, routine, in->name, out);

    /* add to structure */
    entry = (ObitSortStruct*)(out + count * (*size));  /* set pointer to entry */
    entry->index = irow;

    /* Add data by type */
    switch (itype) { 
    case OBIT_float:
      ptrfloat = (ofloat*)(row->myRowData+byteOffset);
      if (desc) entry->key.flt = -(ofloat)ptrfloat[cell];
      else entry->key.flt      =  (ofloat)ptrfloat[cell];
      break;
    case OBIT_double:
      ptrdouble = (odouble*)(row->myRowData+byteOffset);
      if (desc) entry->key.dbl = -(odouble)ptrdouble[cell];
      else entry->key.dbl      =  (odouble)ptrdouble[cell];
      break;
    case OBIT_string:
      ptrchar = (gchar*)(row->myRowData+byteOffset);
      ptrchar += cell*dim[0];  /* Entry in column */
      optr = entry->key.str;
      for (i=0; i<dim[0]; i++) { /* copy string */
	optr[i] = ptrchar[i];
      }
      break;
    case OBIT_int:
      ptrint = (olong*)(row->myRowData+byteOffset);
      if (desc) entry->key.itg = -(olong)ptrint[cell];
      else entry->key.itg      =  (olong)ptrint[cell];
      break;
    case OBIT_byte:
      ptrint8 = (gint8*)(row->myRowData+byteOffset);
      if (desc) entry->key.itg = -(olong)ptrint8[cell];
      else entry->key.itg      =  (olong)ptrint8[cell];
      break;
    case OBIT_short:
      ptrint16 = (gint16*)(row->myRowData+byteOffset);
      if (desc) entry->key.itg = -(olong)ptrint16[cell];
      else entry->key.itg      =  (olong)ptrint16[cell];
      break;
     case OBIT_oint:
      ptroint = (oint*)(row->myRowData+byteOffset);
      if (desc) entry->key.itg = -(olong)ptroint[cell];
      else entry->key.itg      =  (olong)ptroint[cell];
      break;
    case OBIT_long:
      ptrlong = (olong*)(row->myRowData+byteOffset);
      if (desc) entry->key.itg = -(olong)ptrlong[cell];
      else entry->key.itg      =  (olong)ptrlong[cell];
      break;
    case OBIT_uint:
      ptruint = (guint*)(row->myRowData+byteOffset);
      if (desc) entry->key.itg = -(olong)ptruint[cell];
      else entry->key.itg      =  (olong)ptruint[cell];
      break;
    case OBIT_ubyte:
      ptruint8 = (guint8*)(row->myRowData+byteOffset);
      if (desc) entry->key.itg = -(olong)ptruint8[cell];
      else entry->key.itg      =  (olong)ptruint8[cell];
      break;
    case OBIT_ushort:
      ptruint16 = (guint16*)(row->myRowData+byteOffset);
      if (desc) entry->key.itg = -(olong)ptruint16[cell];
      else entry->key.itg      =  (olong)ptruint16[cell];
      break;
    case OBIT_ulong:
      ptrulong = (gulong*)(row->myRowData+byteOffset);
      if (desc) entry->key.itg = -(olong)ptrlong[cell];
      else entry->key.itg      =  (olong)ptrlong[cell];
      break;
    case OBIT_complex:
    case OBIT_dcomplex:
    case OBIT_bool:
    case OBIT_bits:
    default:
      g_assert_not_reached(); /* unknown, barf */
    }; /* end switch */

    count++;  /* How many valid */
  } /* end loop over file */
  
  /* check for errors */
  if ((retCode > OBIT_IO_EOF) || (err->error))
    Obit_traceback_val (err, routine, in->name, out);
  
  /* Release table row */
  row = ObitTableRowUnref (row);
  
  /* Actual number */
  *number = count;

  return out;
} /* end MakeSortStruc */ 

/**
 * Create/fill sort structure for a table using two (float) keys.
 * Integers and doubles converted to float for sort.
 * The sort structure has one "entry" per row (the first  dimension of 
 * a 2D array) which contains 1) a olong row number as an index (1-rel), 
 * and 2 one or more entries of the relevant type. 
 * Each valid row in the table has an entry.
 * \param in     Table to sort, assumed already open;
 * \param which  Columns and element numbers to use as key (1-rel)
 *               in order, key 1 col, key1 elem, key2 col, key2 elem
 *               Key 1 the slowly varying, key2 rapidly
 * \param desc1  If true want descending sort key 1, else ascent
 * \param desc2  If true want descending sort key 2, else ascent
 * \param size   [out] number of bytes in entry
 * \param number [out] number of entries
 * \param type   Data type of key
 * \param ncomp Number of values to compare
 * \return sort structure, use ObitSortStruc to access
 *  should be g_freeed when done.
 */
static gpointer 
MakeSortStruct2f (ObitTable *in, olong which[4], gboolean desc1, 
		  gboolean desc2, olong *size, olong *number, olong *ncomp,
		  ObitInfoType *type, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  gpointer out = NULL;
  ObitTableRow *row;
  ObitSortStruct *entry;
  olong irow, nrow, tsize, count, col1, cell1, col2, cell2, 
    byteOffset1,  byteOffset2;
  ObitInfoType itype;
  gint32 dim[MAXINFOELEMDIM];
  /* Pointers for row data */
  gint8   *ptrint8;
  gint16  *ptrint16;
  olong   *ptrint;
  oint    *ptroint;
  olong   *ptrlong=NULL;
  guint8  *ptruint8;
  guint16 *ptruint16;
  guint   *ptruint;
  gulong  *ptrulong;
  ofloat  *ptrfloat;
  odouble *ptrdouble;
  gchar *routine = "MakeSortStruc2f";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return out;
  g_assert (ObitTableIsA(in)); 

  /* Get table info */
  nrow = in->myDesc->nrow;

  /* which column? Cell? */
  col1  = which[0]-1;
  cell1 = which[1]-1;
  col2  = which[2]-1;
  cell2 = which[3]-1;

  /* Key Column type - converted to floats */
  *type = OBIT_float;

  /* Check that both columns are either float, integer or double */
  itype = in->myDesc->type[col1];
  Obit_retval_if_fail(((itype==OBIT_float) || (itype==OBIT_double) ||
		       (itype==OBIT_short) || (itype==OBIT_ushort) ||
		       (itype==OBIT_long)  || (itype==OBIT_ulong)  ||
		       (itype==OBIT_byte)  || (itype==OBIT_ubyte)  ||
		       (itype==OBIT_int)   || (itype==OBIT_uint)   || 
		       (itype==OBIT_oint)),
		      err, out, "%s: Column %s not numeric", 
		      routine, in->myDesc->FieldName[col1]);
  itype = in->myDesc->type[col2];
  Obit_retval_if_fail(((itype==OBIT_float) || (itype==OBIT_double) ||
		       (itype==OBIT_short) || (itype==OBIT_ushort) ||
		       (itype==OBIT_long)  || (itype==OBIT_ulong)  ||
		       (itype==OBIT_byte)  || (itype==OBIT_ubyte)  ||
		       (itype==OBIT_int)   || (itype==OBIT_uint)   || 
		       (itype==OBIT_oint)),
		      err, out, "%s: Column %s not numeric", 
		      routine, in->myDesc->FieldName[col2]);

  /* element size */
  dim[0] = 2; dim[1] = 1;dim[2] = 1;dim[3] = 1;dim[4] = 1;
  *size = MAX ((sizeof(olong) + ObitInfoElemSize(*type, dim)), sizeof(ObitSortStruct));

  /* Total size of structure in case all rows valid */
  tsize = (*size) * (nrow+10);
  out = g_malloc(tsize);   /* create output structure */
  
  /* Comparing 2 floats */
  *ncomp = 2;
  if (*type == OBIT_string) *ncomp = dim[0];

  /* Create table row */
  row = newObitTableRow (in);

  byteOffset1 = in->myDesc->byteOffset[col1];  /* where it starts in row data */
  byteOffset2 = in->myDesc->byteOffset[col2];  /* where it starts in row data */

 /* loop over table */
  irow = 0;
  count = 0;
  retCode = OBIT_IO_OK;
  while (retCode==OBIT_IO_OK) {
    irow++;
   retCode = ObitTableReadRow (in, irow, row, err);
    if (retCode == OBIT_IO_EOF) break;
    if ((retCode != OBIT_IO_OK) || (err->error)) 
      Obit_traceback_val (err, routine, in->name, out);

    /* add to structure */
    entry = (ObitSortStruct*)(out + count * (*size));  /* set pointer to entry */
    entry->index = irow;

    /* Add first key by type */
    itype = in->myDesc->type[col1];
    switch (itype) { 
    case OBIT_float:
      ptrfloat = (ofloat*)(row->myRowData+byteOffset1);
      if (desc1) entry->key.flt2[0] = -(ofloat)ptrfloat[cell1];
      else entry->key.flt2[0]       =  (ofloat)ptrfloat[cell1];
      break;
    case OBIT_double:
      ptrdouble = (odouble*)(row->myRowData+byteOffset1);
      if (desc1) entry->key.flt2[0] = -(ofloat)ptrdouble[cell1];
      else entry->key.flt2[0]       =  (ofloat)ptrdouble[cell1];
      break;
    case OBIT_int:
      ptrint = (olong*)(row->myRowData+byteOffset1);
      if (desc1) entry->key.flt2[0] = -(ofloat)ptrint[cell1];
      else entry->key.flt2[0]       =  (ofloat)ptrint[cell1];
      break;
    case OBIT_byte:
      ptrint8 = (gint8*)(row->myRowData+byteOffset1);
      if (desc1) entry->key.flt2[0] = -(ofloat)ptrint8[cell1];
      else entry->key.flt2[0]       =  (ofloat)ptrint8[cell1];
      break;
    case OBIT_short:
      ptrint16 = (gint16*)(row->myRowData+byteOffset1);
      if (desc1) entry->key.flt2[0] = -(ofloat)ptrint16[cell1];
      else entry->key.flt2[0]       =  (ofloat)ptrint16[cell1];
      break;
     case OBIT_oint:
      ptroint = (oint*)(row->myRowData+byteOffset1);
      if (desc1) entry->key.flt2[0] = -(ofloat)ptroint[cell1];
      else entry->key.flt2[0]       =  (ofloat)ptroint[cell1];
      break;
    case OBIT_long:
      ptrlong = (olong*)(row->myRowData+byteOffset1);
      if (desc1) entry->key.flt2[0] = -(ofloat)ptrlong[cell1];
      else entry->key.flt2[0]       =  (ofloat)ptrlong[cell1];
      break;
    case OBIT_uint:
      ptruint = (guint*)(row->myRowData+byteOffset1);
      if (desc1) entry->key.flt2[0] = -(ofloat)ptruint[cell1];
      else entry->key.flt2[0]       =  (ofloat)ptruint[cell1];
      break;
    case OBIT_ubyte:
      ptruint8 = (guint8*)(row->myRowData+byteOffset1);
      if (desc1) entry->key.flt2[0] = -(ofloat)ptruint8[cell1];
      else entry->key.flt2[0]       =  (ofloat)ptruint8[cell1];
      break;
    case OBIT_ushort:
      ptruint16 = (guint16*)(row->myRowData+byteOffset1);
      if (desc1) entry->key.flt2[0] = -(ofloat)ptruint16[cell1];
      else entry->key.flt2[0]       =  (ofloat)ptruint16[cell1];
      break;
    case OBIT_ulong:
      ptrulong = (gulong*)(row->myRowData+byteOffset1);
      if (desc1) entry->key.flt2[0] = -(ofloat)ptrlong[cell1];
      else entry->key.flt2[0]       =  (ofloat)ptrlong[cell1];
      break;
    case OBIT_complex:
    case OBIT_dcomplex:
    case OBIT_bool:
    case OBIT_bits:
    default:
      g_assert_not_reached(); /* unknown, barf */
    }; /* end switch key 1 */

    /* Add first key by type */
    itype = in->myDesc->type[col2];
    switch (itype) { 
    case OBIT_float:
      ptrfloat = (ofloat*)(row->myRowData+byteOffset2);
      if (desc2) entry->key.flt2[1] = -(ofloat)ptrfloat[cell2];
      else entry->key.flt2[1]       =  (ofloat)ptrfloat[cell2];
      break;
    case OBIT_double:
      ptrdouble = (odouble*)(row->myRowData+byteOffset2);
      if (desc2) entry->key.flt2[1] = -(ofloat)ptrdouble[cell2];
      else entry->key.flt2[1]       =  (ofloat)ptrdouble[cell2];
      break;
    case OBIT_int:
      ptrint = (olong*)(row->myRowData+byteOffset2);
      if (desc2) entry->key.flt2[1] = -(ofloat)ptrint[cell2];
      else entry->key.flt2[1]       =  (ofloat)ptrint[cell2];
      break;
    case OBIT_byte:
      ptrint8 = (gint8*)(row->myRowData+byteOffset2);
      if (desc2) entry->key.flt2[1] = -(ofloat)ptrint8[cell2];
      else entry->key.flt2[1]       =  (ofloat)ptrint8[cell2];
      break;
    case OBIT_short:
      ptrint16 = (gint16*)(row->myRowData+byteOffset2);
      if (desc2) entry->key.flt2[1] = -(ofloat)ptrint16[cell2];
      else entry->key.flt2[1]       =  (ofloat)ptrint16[cell2];
      break;
     case OBIT_oint:
      ptroint = (oint*)(row->myRowData+byteOffset2);
      if (desc2) entry->key.flt2[1] = -(ofloat)ptroint[cell2];
      else entry->key.flt2[1]       =  (ofloat)ptroint[cell2];
      break;
    case OBIT_long:
      ptrlong = (olong*)(row->myRowData+byteOffset2);
      if (desc2) entry->key.flt2[1] = -(ofloat)ptrlong[cell2];
      else entry->key.flt2[1]       =  (ofloat)ptrlong[cell2];
      break;
    case OBIT_uint:
      ptruint = (guint*)(row->myRowData+byteOffset2);
      if (desc2) entry->key.flt2[1] = -(ofloat)ptruint[cell2];
      else entry->key.flt2[1]       =  (ofloat)ptruint[cell2];
      break;
    case OBIT_ubyte:
      ptruint8 = (guint8*)(row->myRowData+byteOffset2);
      if (desc2) entry->key.flt2[1] = -(ofloat)ptruint8[cell2];
      else entry->key.flt2[1]       =  (ofloat)ptruint8[cell2];
      break;
    case OBIT_ushort:
      ptruint16 = (guint16*)(row->myRowData+byteOffset2);
      if (desc2) entry->key.flt2[1] = -(ofloat)ptruint16[cell2];
      else entry->key.flt2[1]       =  (ofloat)ptruint16[cell2];
      break;
    case OBIT_ulong:
      ptrulong = (gulong*)(row->myRowData+byteOffset2);
      if (desc2) entry->key.flt2[1] = -(ofloat)ptrlong[cell2];
      else entry->key.flt2[1]       =  (ofloat)ptrlong[cell2];
      break;
    case OBIT_complex:
    case OBIT_dcomplex:
    case OBIT_bool:
    case OBIT_bits:
    default:
      g_assert_not_reached(); /* unknown, barf */
    }; /* end switch key 2 */

    count++;  /* How many valid */
  } /* end loop over file */
  
  /* check for errors */
  if ((retCode > OBIT_IO_EOF) || (err->error))
    Obit_traceback_val (err, routine, in->name, out);
  
  /* Release table row */
  row = ObitTableRowUnref (row);
  
  /* Actual number */
  *number = count;

  return out;
} /* end MakeSortStruc2f */ 

/**
 * Compare two lists of integers
 * Conformant to function type GCompareDataFunc
 * \param in1   First list, preceeded by olong index
 * \param in2   Second list, preceeded by olong index
 * \param ncomp Number of values to compare, only does 1
 * \return <0 -> in1 < in2; =0 -> in1 == in2; >0 -> in1 > in2; 
 */
static gint CompareInt (gconstpointer in1, gconstpointer in2, 
			gpointer ncomp)
{
  gint out = 0;
  olong nc, i;
  ObitSortStruct *int1, *int2;
  
  /* get correctly typed local values */
  int1 = (ObitSortStruct*)in1;
  int2 = (ObitSortStruct*)in2;
  nc = *(olong*)ncomp;
  
  /* List or single value? */
  if (nc==1) {
    out = int1->key.itg - int2->key.itg;
  } else { /* list */
    for (i=0; i<nc; i++) {
      out = int1->key.itg2[i] - int2->key.itg2[i];
      if (out) break;   /* stop at first not equal */
    }
  }

  return out;
} /* end CompareInt */

/**
 * Compare ascending abs values of two lists of integers
 * Conformant to function type GCompareDataFunc
 * \param in1   First list, preceeded by olong index
 * \param in2   Second list, preceeded by olong index
 * \param ncomp Number of values to compare
 * \return <0 -> in1 < in2; =0 -> in1 == in2; >0 -> in1 > in2; 
 */
static gint CompareAAInt (gconstpointer in1, gconstpointer in2, 
			gpointer ncomp)
{
  gint out = 0;
  olong nc, i;
  ObitSortStruct *int1, *int2;

  /* get correctly typed local values */
  int1 = (ObitSortStruct*)in1;
  int2 = (ObitSortStruct*)in2;
  nc = *(olong*)ncomp;

  /* List or single value? */
  if (nc==1) {
    out = abs(int1->key.itg) - abs(int2->key.itg);
  } else { /* list */
    for (i=0; i<nc; i++) {
      out = abs(int1->key.itg2[i]) - abs(int2->key.itg2[i]);
      if (out) break;   /* stop at first not equal */
    }
  }

  return out;
} /* end CompareAAInt */


/**
 * Compare descending abs values of two lists of integers
 * Conformant to function type GCompareDataFunc
 * \param in1   First list, preceeded by olong index
 * \param in2   Second list, preceeded by olong index
 * \param ncomp Number of values to compare
 * \return <0 -> in1 < in2; =0 -> in1 == in2; >0 -> in1 > in2; 
 */
static gint CompareADInt (gconstpointer in1, gconstpointer in2, 
			gpointer ncomp)
{
  gint out = 0;
  olong nc, i;
  ObitSortStruct *int1, *int2;

  /* get correctly typed local values */
  int1 = (ObitSortStruct*)in1;
  int2 = (ObitSortStruct*)in2;
  nc = *(olong*)ncomp;

  /* List or single value? */
  if (nc==1) {
    out = abs(int2->key.itg) - abs(int1->key.itg);
  } else { /* list */
    for (i=0; i<nc; i++) {
      out = abs(int2->key.itg2[i]) - abs(int2->key.itg2[i]);
      if (out) break;   /* stop at first not equal */
    }
  }

  return out;
} /* end CompareADInt */

/**
 * Compare two lists of floats
 * Conformant to function type GCompareDataFunc
 * \param in1   First list, preceeded by olong index
 * \param in2   Second list, preceeded by olong index
 * \param ncomp Number of values to compare
 * \return <0 -> in1 < in2; =0 -> in1 == in2; >0 -> in1 > in2; 
 */
static gint CompareFloat (gconstpointer in1, gconstpointer in2, 
			  gpointer ncomp)
{
  gint out = 0;
  olong nc, i;
  ObitSortStruct *float1, *float2;

  /* get correctly typed local values */
  float1 = (ObitSortStruct*)in1;
  float2 = (ObitSortStruct*)in2;
  nc = *(olong*)ncomp;

  /* List or single value? */
  if (nc==1) {
    if (float1->key.flt<float2->key.flt)      out = -1;
    else if (float1->key.flt>float2->key.flt) out = 1;
    else                                      out = 0;
  } else { /* list */
    for (i=0; i<nc; i++) {
      if (float1->key.flt2[i]<float2->key.flt2[i])      out = -1;
      else if (float1->key.flt2[i]>float2->key.flt2[i]) out = 1;
      else                                              out = 0;
      if (out) break;   /* stop at first not equal */
    }
  }

  return out;
} /* end CompareFloat */

/*
 * Compare ascending abs values of two lists of floats
 * Conformant to function type GCompareDataFunc
 * \param in1   First list, preceeded by olong index
 * \param in2   Second list, preceeded by olong index
 * \param ncomp Number of values to compare
 * \return <0 -> in1 < in2; =0 -> in1 == in2; >0 -> in1 > in2; 
 */
static gint CompareAAFloat (gconstpointer in1, gconstpointer in2, 
			  gpointer ncomp)
{
  gint out = 0;
  olong nc, i;
  ObitSortStruct *float1, *float2;

  /* get correctly typed local values */
  float1 = (ObitSortStruct*)in1;
  float2 = (ObitSortStruct*)in2;
  nc = *(olong*)ncomp;

  /* List or single value? */
  if (nc==1) {
    if (fabs(float1->key.flt)<fabs(float2->key.flt))      out = -1;
    else if (fabs(float1->key.flt)>fabs(float2->key.flt)) out = 1;
    else                                                  out = 0;
  } else { /* list */
    for (i=0; i<nc; i++) {
      if (fabs(float1->key.flt2[i])<fabs(float2->key.flt2[i]))      out = -1;
      else if (fabs(float1->key.flt2[i])>fabs(float2->key.flt2[i])) out = 1;
      else                                                          out = 0;
      if (out) break;   /* stop at first not equal */
    }
  }

  return out;
} /* end CompareAAFloat */

/*
 * Compare descending abs values of two lists of floats
 * Conformant to function type GCompareDataFunc
 * \param in1   First list, preceeded by olong index
 * \param in2   Second list, preceeded by olong index
 * \param ncomp Number of values to compare
 * \return <0 -> in1 < in2; =0 -> in1 == in2; >0 -> in1 > in2; 
 */
static gint CompareADFloat (gconstpointer in1, gconstpointer in2, 
			  gpointer ncomp)
{
  gint out = 0;
  olong nc, i;
  ObitSortStruct *float1, *float2;

  /* get correctly typed local values */
  float1 = (ObitSortStruct*)in1;
  float2 = (ObitSortStruct*)in2;
  nc = *(olong*)ncomp;

  /* List or single value? */
  if (nc==1) {
    if (fabs(float1->key.flt)<fabs(float2->key.flt))      out =  1;
    else if (fabs(float1->key.flt)>fabs(float2->key.flt)) out = -1;
    else                                                  out = 0;
  } else { /* list */
    for (i=0; i<nc; i++) {
      if (fabs(float1->key.flt2[i])<fabs(float2->key.flt2[i]))      out =  1;
      else if (fabs(float1->key.flt2[i])>fabs(float2->key.flt2[i])) out = -1;
      else                                                          out = 0;
      if (out) break;   /* stop at first not equal */
    }
  }

  return out;
} /* end CompareADFloat */

/**
 * Compare two lists of double
 * Conformant to function type GCompareDataFunc
 * \param in1   First list, preceeded by olong index
 * \param in2   Second list, preceeded by olong index
 * \param ncomp Number of values to compare
 * \return <0 -> in1 < in2; =0 -> in1 == in2; >0 -> in1 > in2; 
 */
static gint CompareDouble (gconstpointer in1, gconstpointer in2, 
			   gpointer ncomp)
{
  gint out = 0;
  olong nc, i;
  ObitSortStruct *double1, *double2;

  /* get correctly typed local values */
  double1 = (ObitSortStruct*)in1;
  double2 = (ObitSortStruct*)in2;
  nc = *(olong*)ncomp;

  /* List or single value? */
  if (nc==1) {
    if (double1->key.dbl<double2->key.dbl)      out = -1;
    else if (double1->key.dbl>double2->key.dbl) out = 1;
    else                                        out = 0;
  } else { /* list */
    for (i=0; i<nc; i++) {
      if (double1->key.dbl2[i]<double2->key.dbl2[i])      out = -1;
      else if (double1->key.dbl2[i]>double2->key.dbl2[i]) out = 1;
      else                                                out = 0;
      if (out) break;   /* stop at first not equal */
    }
  }

  return out;
} /* end CompareDouble */

/**
 * Compare ascending absolute values of two lists of double
 * Conformant to function type GCompareDataFunc
 * \param in1   First list, preceeded by olong index
 * \param in2   Second list, preceeded by olong index
 * \param ncomp Number of values to compare
 * \return <0 -> in1 < in2; =0 -> in1 == in2; >0 -> in1 > in2; 
 */
static gint CompareAADouble (gconstpointer in1, gconstpointer in2, 
			     gpointer ncomp)
{
  gint out = 0;
  olong nc, i;
  ObitSortStruct *double1, *double2;

  /* get correctly typed local values */
  double1 = (ObitSortStruct*)in1;
  double2 = (ObitSortStruct*)in2;
  nc = *(olong*)ncomp;

  /* List or single value? */
  if (nc==1) {
    if (fabs(double1->key.dbl)<fabs(double2->key.dbl))      out = -1;
    else if (fabs(double1->key.dbl)>fabs(double2->key.dbl)) out = 1;
    else                                                    out = 0;
  } else { /* list */
    for (i=0; i<nc; i++) {
      if (fabs(double1->key.dbl2[i])<fabs(double2->key.dbl2[i]))      out = -1;
      else if (fabs(double1->key.dbl2[i])>fabs(double2->key.dbl2[i])) out = 1;
      else                                                            out = 0;
      if (out) break;   /* stop at first not equal */
    }
  }

  return out;
} /* end CompareAADouble */

/**
 * Compare descending absolute values of two lists of double
 * Conformant to function type GCompareDataFunc
 * \param in1   First list, preceeded by olong index
 * \param in2   Second list, preceeded by olong index
 * \param ncomp Number of values to compare
 * \return <0 -> in1 < in2; =0 -> in1 == in2; >0 -> in1 > in2; 
 */
static gint CompareADDouble (gconstpointer in1, gconstpointer in2, 
			     gpointer ncomp)
{
  gint out = 0;
  olong nc, i;
  ObitSortStruct *double1, *double2;

  /* get correctly typed local values */
  double1 = (ObitSortStruct*)in1;
  double2 = (ObitSortStruct*)in2;
  nc = *(olong*)ncomp;

  /* List or single value? */
  if (nc==1) {
    if (fabs(double1->key.dbl)<fabs(double2->key.dbl))      out =  1;
    else if (fabs(double1->key.dbl)>fabs(double2->key.dbl)) out = -1;
    else                                                    out = 0;
  } else { /* list */
    for (i=0; i<nc; i++) {
      if (fabs(double1->key.dbl2[i])<fabs(double2->key.dbl2[i]))      out =  1;
      else if (fabs(double1->key.dbl2[i])>fabs(double2->key.dbl2[i])) out = -1;
      else                                                            out = 0;
      if (out) break;   /* stop at first not equal */
    }
  }

  return out;
} /* end CompareADDouble */

/**
 * Compare two (ascii) character strings
 * Conformant to function type GCompareDataFunc
 * \param in1   First list, preceeded by olong index
 * \param in2   Second list, preceeded by olong index
 * \param ncomp Number of characters to compare
 * \return <0 -> in1 < in2; =0 -> in1 == in2; >0 -> in1 > in2; 
 */
static gint CompareString (gconstpointer in1, gconstpointer in2, 
			   gpointer ncomp)
{
  gint out = 0;
  olong nc, i;
  ObitSortStruct *string1, *string2;

  /* get correctly typed local values */
  string1 = (ObitSortStruct*)in1;
  string2 = (ObitSortStruct*)in2;
  nc = *(olong*)ncomp;

  /* List or single value? */
  if (nc==1) {
    out = string1->key.str[0] - string2->key.str[0];
  } else { /* list */
    for (i=0; i<nc; i++) {
      out = string1->key.str[i] - string2->key.str[i];
      if (out) break;   /* stop at first not equal */
    }
  }

  return out;
} /* end CompareString */

/**
 * Reorders table
 * Makes scratch image with a scratch table for the reordered
 * table and then copies over the input table.
 * \param in       Table to sort
 * \param base     Base address of sort structure
 * \param size     Size in bytes of a sort element
 * \param number   Number of sort elements
 * \param sortCol1 Sort column number key 1 (1-rel)
 * \param sortCol2 Sort column number key 2 (1-rel), 0=not used
 * \param err     ObitErr error stack.
 * \return I/O Code  OBIT_IO_OK = OK.
 */
static ObitIOCode 
ReorderTable(ObitTable *in, gpointer base, olong size, olong number, 
	     olong sortCol1, olong sortCol2, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitData *scrData = NULL;
  ObitTable *scrTable = NULL;
  ObitTableRow *row = NULL;
  ObitSortStruct *entry;
  ObitIOStatus saveStatus;
  olong irow, orow, tabVer;
  ollong count;
  gchar *routine = "ReorderTable";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return retCode;

  /* Save status on host object */
  saveStatus = ((ObitData*)in->myHost)->myStatus;

  /* scratch Data */
  scrData = newObitDataScratch ((ObitData*)in->myHost, err);
  if ((scrData==NULL) || (err->error))  goto cleanup;

  /* Clone input data */
  ObitDataClone ((ObitData*)in->myHost, scrData, err);
  if (err->error) goto cleanup;

  /* Restore status on host object */
  ((ObitData*)in->myHost)->myStatus = saveStatus;

  /* Make scratch table */
  tabVer = 0;
  scrTable = newObitDataTable(scrData, OBIT_IO_WriteOnly, in->myDesc->TableName, 
			       &tabVer, err);
  if ((scrTable==NULL) || (err->error))  goto cleanup;

  /* Clone input table */
  scrTable = ObitTableClone (in, scrTable);
 
 /* need independent descriptors */
  scrTable->myDesc = ObitUnref(scrTable->myDesc); /* release old */
  scrTable->myDesc = ObitTableDescCopy(in->myDesc, scrTable->myDesc, err);
  if (err->error) goto cleanup;

  /* Open input table */
  retCode = ObitTableOpen (in, OBIT_IO_ReadWrite, err);
  if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;

  /* Create table row */
  row = newObitTableRow (in);

  retCode = ObitTableOpen (scrTable, OBIT_IO_WriteOnly, err);
  if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;

  /* Attach row to output buffer */
  ObitTableSetRow (scrTable, row, err);
  if (err->error) goto cleanup;

  /* loop over table */
  irow = 0;
  count = 0;
  while ((retCode==OBIT_IO_OK) && (count<number)) {
 
    /* which row - get from sort structure */
    entry = base + count * (size);  /* set pointer to entry */
    irow = entry->index;
    retCode = ObitTableReadRow (in, irow, row, err);
    if (retCode == OBIT_IO_EOF) break;
    if ((retCode != OBIT_IO_OK) || (err->error))  goto cleanup;

    orow = -1;
    retCode = ObitTableWriteRow (scrTable, orow, row, err);
    if ((retCode != OBIT_IO_OK) || (err->error))  goto cleanup;

    count++;
  } /* end loop reordering table table */

  /* Mark table as sorted */
  scrTable->myDesc->sort[0] = sortCol1;
  scrTable->myDesc->sort[1] = sortCol2;
  ((ObitTableDesc*)scrTable->myIO->myDesc)->sort[0] = sortCol1;
  ((ObitTableDesc*)scrTable->myIO->myDesc)->sort[1] = sortCol2;

  /* Truncate input table in preparation for copy back */
  in->myDesc->nrow = 0;
  ((ObitTableDesc*)in->myIO->myDesc)->nrow = 0;
 
  retCode = ObitTableClose (in, err);
  if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;
  retCode = ObitTableClose (scrTable, err);
  if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;

  /* Copy table back */
  ObitTableCopy (scrTable, in, err);
  if ((retCode != OBIT_IO_OK) || (err->error)) goto cleanup;

  /* Cleanup */
  cleanup: row = ObitTableRowUnref (row); /* Release table row */
  scrTable = ObitTableUnref(scrTable);
  scrData  = ObitDataZap (scrData, err);
  scrData = ObitDataUnref (scrData);
  if (err->error) Obit_traceback_val (err, routine, in->name, retCode);
 
  return OBIT_IO_OK;
} /* end ReorderTable */
