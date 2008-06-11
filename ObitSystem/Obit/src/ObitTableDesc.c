/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2008                                          */
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
#include "Obit.h"
#include "ObitTableDesc.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitTableDesc.c
 * ObitTableDesc Obit Table descriptor class definition.
 * This contains information about the structure of the table.
 */

/*--------------- File Global Variables  ----------------*/
/** name of the class defined in this file */
static gchar *myClassName = "ObitTableDesc";

/**
 * ClassInfo global structure ObitIOClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitTableDescClassInfo myClassInfo = {FALSE};

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitTableDescInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitTableDescClear (gpointer in);

/*---------------Public functions---------------------------*/
/**
 * Construct Object.
 * \return pointer to object created.
 */
ObitTableDesc* newObitTableDesc (gchar *name)
{
  ObitTableDesc *out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitTableDescClassInit();

  /* allocate structure */
  out = g_malloc0(sizeof(ObitTableDesc));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

 /* set classInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitTableDescInit((gpointer)out);

 return out;
} /* end newObitTableDesc */

/**
 * Returns ClassInfo pointer for the class.
 * Initializes class if needed on first call.
 * \return pointer to the class structure.
 */
gconstpointer ObitTableDescGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitTableDescClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitTableDescGetClass */

/**
 * Copy constructor.
 * The output descriptor will have the structure and values of the input
 * \param in Pointer to object to be copied.
 * \param out Pointer to object to be written.  
 *            If NULL then a new structure is created.
 * \param err ObitErr error stack
 * \return Pointer to new object.
 */
ObitTableDesc* ObitTableDescCopy (ObitTableDesc* in, ObitTableDesc* out, 
			    ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  gboolean oldExist;
  gchar *outName;
  olong i, j;

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return out;
  g_assert (ObitIsA(in, &myClassInfo));
  if (out) g_assert (ObitIsA(out, &myClassInfo));

  /* Don't bother it they are the same */
  if (in==out) return out;

  /* Create if it doesn't exist */
  oldExist = out!=NULL;
  if (!oldExist) {
    /* derive object name */
    outName = g_strconcat ("Copy: ",in->name,NULL);
    out = newObitTableDesc(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /* Free any existing array members */
  if ((in->info)&&(out->info)) out->info = ObitInfoListUnref (out->info); 
  if (out->TableName) g_free(out->TableName);   out->TableName= NULL;
  if (out->type)      g_free(out->type);        out->type     = NULL;
  if (out->order)     g_free(out->order);       out->order    = NULL;
  if (out->repeat)    g_free(out->repeat);      out->repeat   = NULL;
  if (out->offset)    g_free(out->offset);      out->offset   = NULL;
  if (out->byteOffset) g_free(out->byteOffset); out->byteOffset = NULL;
  /* Loop over fields */
  if (out->FieldName) {
    for (i=0; i<out->nfield; i++) {
      if (out->FieldName[i]) g_free(out->FieldName[i]);
      out->FieldName[i] = NULL;
    }
    g_free(out->FieldName); out->FieldName = NULL;
  }
  if (out->FieldUnit) {
    for (i=0; i<out->nfield; i++) {
      if (out->FieldUnit[i]) g_free(out->FieldUnit[i]);
      out->FieldUnit[i] = NULL;
    }
    g_free(out->FieldUnit); out->FieldUnit = NULL;
  }
  if (out->dim) {
    for (i=0; i<out->nfield; i++) {
      if (out->dim[i]) g_free(out->dim[i]); out->dim[i] = NULL;
    }
    g_free(out->dim); out->dim = NULL;
  }

  /* initialize/copy */
  out->access  = in->access;
  out->nrow    = in->nrow;
  out->lrow    = in->lrow;
  out->version = in->version;
  out->nfield  = in->nfield;
  out->nkey    = in->nkey;
  out->sort[0] = in->sort[0];
  out->sort[1] = in->sort[1];
  out->startData = in->startData;

  /* Copy array members */
  if (out->TableName) g_free(out->TableName);
  if (in->TableName) out->TableName = g_strdup(in->TableName);
  else out->TableName = g_strdup("Unnamed");
  if (in->nfield>0) {
    out->type      = g_malloc0(in->nfield*sizeof(ObitInfoType));
    out->offset    = g_malloc0(in->nfield*sizeof(olong));
    out->byteOffset= g_malloc0(in->nfield*sizeof(olong));
    out->repeat    = g_malloc0(in->nfield*sizeof(olong));
    out->order     = g_malloc0(in->nfield*sizeof(olong));
    out->dim       = g_malloc0(in->nfield*sizeof(gint32*));
    out->FieldName = g_malloc0(in->nfield*sizeof(gchar*));
    out->FieldUnit = g_malloc0(in->nfield*sizeof(gchar*));
  }
  if (in->info!=NULL) out->info = ObitInfoListCopy (in->info);
  /* Loop over fields */
  for (i=0; i<out->nfield; i++) {
    out->FieldName[i] =  g_strdup(in->FieldName[i]);
    out->FieldUnit[i] =  g_strdup(in->FieldUnit[i]);
    out->type[i]      = in->type[i];
    out->offset[i]    = in->offset[i];
    out->byteOffset[i]= in->byteOffset[i];
    out->repeat[i]    = in->repeat[i];
    out->order[i]     = in->order[i];
    out->dim[i]       = g_malloc0(MAXINFOELEMDIM*sizeof(gint32));
    for (j=0; j<MAXINFOELEMDIM; j++) out->dim[i][j] = in->dim[i][j];
 }

  /* index output */
  ObitTableDescIndex (out);

  return out;
} /* end ObitTableDescCopy */

/**
 * Copy descriptive material (i.e. things that don't define the structure).
 * \param in  Pointer to object to be copied.
 * \param out Pointer to object to be written.  
 * \param err ObitErr error stack
 */
void ObitTableDescCopyDesc (ObitTableDesc* in, ObitTableDesc* out, 
			    ObitErr *err)
{
  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitIsA(out, &myClassInfo));

  /* Copy Sort order */
  out->sort[0] = in->sort[0];
  out->sort[1] = in->sort[1];
 
  /* index output */
  ObitTableDescIndex (out);

  return;
} /* end ObitTableDescCopyDesc */

/**
 * Set column (field) order, and offset..
 * \param in Pointer to object.
 */
void ObitTableDescIndex (ObitTableDesc* in)
{
  olong i, j, size, nbytes, total, order, maxsize, nadd, it, nelem, nf;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  /* OBIT_bits the same as OBIT_oint - DAMN AIPS */
  ObitInfoType otypes[10] = {OBIT_double, OBIT_float, OBIT_string, 
			    OBIT_oint, OBIT_bits, OBIT_long, OBIT_bool,
			    OBIT_int, OBIT_short, OBIT_ubyte};
  gboolean damnit, status_col=FALSE;

  /* error check */
  g_assert (ObitIsA(in, &myClassInfo));

  /* don't bother if there's nothing in the table yet */
  if (in->nfield<=0) return;

  /* create some arrays if not defined */
  if (!in->offset)    in->offset    = g_malloc0(in->nfield*sizeof(olong));
  if (!in->byteOffset)in->byteOffset= g_malloc0(in->nfield*sizeof(olong));
  if (!in->repeat)    in->repeat    = g_malloc0(in->nfield*sizeof(olong));
  if (!in->order)     in->order     = g_malloc0(in->nfield*sizeof(olong));

  /* the number of elements in each column */
  for (i=0; i<in->nfield; i++) {
    in->order[i] = 0;
    nelem = in->dim[i][0]; /* may be zero */
    for (j=1; j<MAXINFOELEMDIM; j++) nelem *= MAX (1, in->dim[i][j]);
    in->repeat[i] = nelem;
  }

  /* determine physical order - by AIPS order (size) of element */
  /* loop over possible types in AIPS order */
  /* _status column is last */
  if (in->FieldName && in->FieldName[in->nfield-1]) {
    status_col = (!strncmp (in->FieldName[in->nfield-1], "_status", 7));
  } else {
    status_col = FALSE;
  }
  nf = in->nfield;
  if (status_col) nf--;
  order = 1;
  maxsize = 1;  /* Maximum size of element in bytes */
  for (i=0; i<10; i++) {
    /* Loop over fields */
    for (j=0; j<nf; j++) {
      /* How big is an element of this one? */
      size = ObitInfoElemSize(in->type[j], dim);
      /* Want it? */
      /* bits and ints the same */
      damnit = (in->type[j]==OBIT_bits) && (otypes[i]==OBIT_oint);
      if (damnit || (in->type[j]==otypes[i])) { /* Yes */
	if (in->order[j]==0) in->order[j] = order++;
	maxsize = MAX (maxsize, size);
      }
    }
  }
  /* AIPS Status column last */
  if (status_col) in->order[in->nfield-1] = order;

  /* determine offsets from beginning of row (in units of that type) */
  total = 0;
  for (i=0; i<in->nfield; i++) {
    /* Count how many bytes preceed this entry */
    nbytes = 0;
    for (j=0; j<in->nfield; j++) {
      if(in->order[i]>in->order[j]) {
	nadd = ObitInfoElemSize(in->type[j], in->dim[j]);
	/* Strings are in 4 char/float */
	if (in->type[j]==OBIT_string) {
	  nadd = in->repeat[j];
	  it = 1 + (nadd-1)/4;
	  nadd = it * sizeof(ofloat);
	  /* bit arrays are 32 /float */
	} else if (in->type[j]==OBIT_bits) {
	  nadd = in->repeat[j];
	  it = 1 + (nadd-1)/32;
	  nadd = it * sizeof(ofloat);
	  /* logical the size of a float */
	} else if (in->type[j]==OBIT_bool) {
	  nadd = sizeof(ofloat) * in->repeat[j];
	} else if (in->type[j]==OBIT_double) {
	  nadd = sizeof(odouble) * in->repeat[j];
	} else if (in->type[j]==OBIT_float) {
	  nadd = sizeof(ofloat) * in->repeat[j];
	} else if (in->type[j]==OBIT_long) {
	  nadd = sizeof(olong) * in->repeat[j];
	} else if (in->type[j]==OBIT_int) {
	  nadd = sizeof(olong) * in->repeat[j];
	} else if (in->type[j]==OBIT_oint) {
	  nadd = sizeof(oint) * in->repeat[j];
	} else if (in->type[j]==OBIT_byte) {
	  nadd = in->repeat[j];
	} else { /* anything else */
	  nadd = ObitInfoElemSize(in->type[j], in->dim[j]);
	}
	nbytes += nadd;
      }
    }
    /* offset in bytes */
    in->byteOffset[i] = nbytes;

    /* offset in units of words of the appropriate type */
    if (in->type[i]==OBIT_string) {
      in->offset[i] = nbytes;
      /* bit arrays are 32 /float */
    } else if (in->type[i]==OBIT_bits) {
      in->offset[i] = nbytes / sizeof (ofloat);
      /* logical the size of a float */
      in->offset[i] = nbytes / sizeof (ofloat);
    } else if (in->type[i]==OBIT_bool) {
      in->offset[i] = nbytes / sizeof (ofloat);
    } else if (in->type[i]==OBIT_double) {
      in->offset[i] = nbytes / sizeof (odouble);
    } else if (in->type[i]==OBIT_float) {
      in->offset[i] = nbytes / sizeof (ofloat);
    } else if (in->type[i]==OBIT_long) {
      in->offset[i] = nbytes / sizeof (olong);
    } else if (in->type[i]==OBIT_int) {
      in->offset[i] = nbytes / sizeof (olong);
    } else if (in->type[i]==OBIT_byte) {
      in->offset[i] = nbytes;
    } else { /* anything else */
      in->offset[i] = nbytes / ObitInfoElemSize(in->type[i], dim);
    }

    /* Accumulate total in bytes */
    total += ObitInfoElemSize(in->type[i], in->dim[i]);
  } /* End loop over fields */

  /* total length */
  /* make row length an integral multiple of the biggest
     element size to avoid trouble with arrays of row data */
  nelem = 1 + ((total - 1)/maxsize);
  in->lrow = nelem * maxsize;

  /* For IO Need the actual size */
  in->lrowIO = total;

  /* Find "_status" column */
  in->statusOff = -1;
  if (in->FieldName) {
    for (i=0; i<in->nfield; i++) {
      if (in->FieldName[i]) {
	if (!strncmp (in->FieldName[i], "_status", 7)) 
	  in->statusOff = in->offset[i];
      }
    }
  }

} /* end ObitTableDescIndex */

/**
 * Free any existing versions of column based arrays
 * and reallocate.
 * \param in     Pointer to object to reallocate.
 * \param nfield number of columns to allocate
 */
void ObitTableDescRealloc (ObitTableDesc* in, olong nfield)
{
  olong i;

  /* error check */
  g_assert (ObitIsA(in, &myClassInfo));
  if (nfield<1) return; /* don't bother */

  in->nfield = nfield;

  /* Free any existing array members  in descriptor */
  if (in->info)       in->info = ObitInfoListUnref (in->info); 
  if (in->type)       g_free(in->type);        in->type     = NULL;
  if (in->order)      g_free(in->order);       in->order    = NULL;
  if (in->repeat)     g_free(in->repeat);      in->repeat   = NULL;
  if (in->offset)     g_free(in->offset);      in->offset   = NULL;
  if (in->byteOffset) g_free(in->byteOffset);  in->byteOffset = NULL;
  /* Loop over fields */
  if (in->FieldName) {
    for (i=0; i<in->nfield; i++) {
      if (in->FieldName[i]) g_free(in->FieldName[i]);
      in->FieldName[i] = NULL;
    }
    if (in->FieldName) g_free(in->FieldName); in->FieldName = NULL;
  }
  if (in->FieldUnit) {
    for (i=0; i<in->nfield; i++) {
      if (in->FieldUnit[i]) g_free(in->FieldUnit[i]);
      in->FieldUnit[i] = NULL;
    }
    if (in->FieldUnit) g_free(in->FieldUnit); in->FieldUnit = NULL;
  }
  if (in->dim) {
    for (i=0; i<in->nfield; i++) {
      if (in->dim[i]) g_free(in->dim[i]); in->dim[i] = NULL;
    }
    if (in->dim) g_free(in->dim); in->dim = NULL;
  }

  /* create new versions as needed */
  in->info       = newObitInfoList();
  in->type       = g_malloc0(in->nfield*sizeof(ObitInfoType));
  in->offset     = g_malloc0(in->nfield*sizeof(olong));
  in->byteOffset = g_malloc0(in->nfield*sizeof(olong));
  in->repeat     = g_malloc0(in->nfield*sizeof(olong));
  in->order      = g_malloc0(in->nfield*sizeof(olong));
  in->FieldName  = g_malloc0(in->nfield*sizeof(gchar*));
  in->FieldUnit  = g_malloc0(in->nfield*sizeof(gchar*));
  in->dim        = g_malloc0(in->nfield*sizeof(gint32*));
  for (i=0; i<in->nfield; i++) 
    in->dim[i]  = g_malloc0(MAXINFOELEMDIM*sizeof(gint32));
}  /* end ObitTableDescRealloc */

/**
 * Determine if two table descriptors indicate the same structure 
 * of table rows. 
 * Checks number, type, label and dimension of each column
 * \param in1     Pointer first object
 * \param in2     Pointer second object
 * \return TRUE if Compatible, else FALSE
 */
gboolean ObitTableDescCompatible (ObitTableDesc* in1, ObitTableDesc* in2)
{
  olong i;

  /* Number of fields */
  if (in1->nfield != in2->nfield) return FALSE;

  /* Overall size */
  if (in1->lrow != in2->lrow) return FALSE;

  /* Loop over fields */
  for (i=0; i<in1->nfield; i++) {
    /* Check type */
    if (in1->type[i] != in2->type[i]) return FALSE;

    /* Check repeat (number of values) */
    if (in1->repeat[i] != in2->repeat[i]) return FALSE;

    /* Check physical position in row */
    if (in1->order[i] != in2->order[i]) return FALSE;

    /* Check field labels */
    if (strncmp (in1->FieldName[i], in2->FieldName[i], 100)) return FALSE;
  } /* end loop over fields */

  /* If it gets here it, they must be the same */
  return TRUE;
}  /* end ObitTableDescCompatible */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitTableDescClassInit (void)
{
  const ObitClassInfo *ParentClass;

  if (myClassInfo.initialized) return;  /* only once */
  myClassInfo.initialized = TRUE;

   /* Initialize (recursively) parent class first */
  ParentClass = ObitGetClass();
  if ((ParentClass!=NULL) && (ParentClass->ObitClassInit!=NULL))
    ParentClass->ObitClassInit();

  /* function pointers etc. for this class */
  myClassInfo.ClassName     = g_strdup(myClassName);
  myClassInfo.ParentClass   = ParentClass;
  myClassInfo.ObitClassInit = 
    (ObitClassInitFP)ObitTableDescClassInit;
  myClassInfo.newObit       = (newObitFP)newObitTableDesc;
  myClassInfo.ObitCopy      = (ObitCopyFP)ObitTableDescCopy;
  myClassInfo.ObitClone     = NULL;
  myClassInfo.ObitRef       = (ObitRefFP)ObitRef;
  myClassInfo.ObitUnref     = (ObitUnrefFP)ObitUnref;
  myClassInfo.ObitIsA       = (ObitIsAFP)ObitIsA;
  myClassInfo.ObitClear     = (ObitClearFP)ObitTableDescClear;
  myClassInfo.ObitInit      = (ObitInitFP)ObitTableDescInit;
} /* end ObitTableSelClassInit */

/*---------------Private functions--------------------------*/
/**
 * Creates empty member objects, initialize reference count.
 * Does (recursive) initialization of base class members before 
 * this class.
 * \param inn Pointer to the object to initialize.
 */
void ObitTableDescInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitTableDesc *in = inn;

  /* error checks */
  g_assert (in != NULL);
  
  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->info      = newObitInfoList();
  in->nrow      = 0;
  in->lrow      = 0;
  in->nkey      = 0;
  in->statusOff = -1;
  in->TableName = NULL;
  in->FieldName = NULL;
  in->FieldUnit = NULL;
  in->type      = NULL;
  in->dim       = NULL;
  in->offset    = NULL;
  in->byteOffset= NULL;
  in->repeat    = NULL;
  in->order     = NULL;
} /* end ObitTableDescInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 */
void ObitTableDescClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitTableDesc *in = inn;
  olong i;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* free this class members */
  if (in->TableName) g_free(in->TableName);   in->TableName=NULL;
  if (in->type)      g_free(in->type);        in->type     = NULL;
  if (in->order)     g_free(in->order);       in->order    = NULL;
  if (in->repeat)    g_free(in->repeat);      in->repeat   = NULL;
  if (in->offset)    g_free(in->offset);      in->offset   = NULL;
  if (in->byteOffset) g_free(in->byteOffset); in->byteOffset = NULL;
  if (in->info) ObitInfoListUnref (in->info); in->info     = NULL;

 /* Loop over fields - field labels */
  if (in->FieldName) {
    for (i=0; i<in->nfield; i++) {
      if (in->FieldName[i]) g_free(in->FieldName[i]);
      in->FieldName[i] = NULL;
    }
    if (in->FieldName) g_free(in->FieldName); in->FieldName = NULL;
  }
  /* Field units string */
  if (in->FieldUnit) {
    for (i=0; i<in->nfield; i++) {
      if (in->FieldUnit[i]) g_free(in->FieldUnit[i]);
      in->FieldUnit[i] = NULL;
    }
    if (in->FieldUnit) g_free(in->FieldUnit); in->FieldUnit = NULL;
  }

  /* now dim array */
  if (in->dim) {
    for (i=0; i<in->nfield; i++) {
      if (in->dim[i]) g_free(in->dim[i]); in->dim[i] = NULL;
    }
    if (in->dim) g_free(in->dim); in->dim = NULL;
  }
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);

} /* end ObitTableDescClear */



