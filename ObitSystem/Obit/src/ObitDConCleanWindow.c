/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2004-2010                                          */
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

#include "ObitDConCleanWindow.h"
#include "ObitMem.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitDConCleanWindow.c
 * ObitDConCleanWindow class function definitions.
 * This class is derived from the Obit base class.
 */

/** name of the class defined in this file */
static gchar *myClassName = "ObitDConCleanWindow";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo structure ObitDConCleanWindowClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitDConCleanWindowClassInfo myClassInfo = {FALSE};

/*-------------- Private Class definitions-------------------------*/
/**  WindowListElem structure */
typedef struct {
  /** Id */
  olong Id;
  /** Window type */
  ObitDConCleanWindowType type;
  /** window definition */
  olong window[4];
} WindowListElem;

/*--------------- File Global Variables  ----------------*/


/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitDConCleanWindowInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitDConCleanWindowClear (gpointer in);

/** Private: Determine window size. */
olong  GetWindowSize(ObitFArray *image, olong *PeakPos, ofloat sigma);

/** Private: WindowListElem constructor */
WindowListElem*  newWindowListElem (olong Id, ObitDConCleanWindowType type, 
				    olong *window);
/** Private: WindowListElem destructor. */
WindowListElem* freeWindowListElem (WindowListElem *in);

/** Private: Find Window of given ID. */
WindowListElem*  ObitDConCleanWindowFind (GSList *glist, olong Id);

/** Private: Set Class function pointers. */
static void ObitDConCleanWindowClassInfoDefFn (gpointer inClass);

/*----------------------Public functions---------------------------*/
/**
 * Constructor.
 * Initializes class if needed on first call.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitDConCleanWindow* newObitDConCleanWindow (gchar* name)
{
  ObitDConCleanWindow* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitDConCleanWindowClassInit();

  /* allocate/init structure */
  out = g_malloc0(sizeof(ObitDConCleanWindow));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

  /* set ClassInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitDConCleanWindowInit((gpointer)out);

 return out;
} /* end newObitDConCleanWindow */

/**
 * Returns ClassInfo pointer for the class.
 * \return pointer to the class structure.
 */
gconstpointer ObitDConCleanWindowGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitDConCleanWindowClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitDConCleanWindowGetClass */

/**
 * Make a deep copy of an ObitDConCleanWindow.
 * \param in  The object to copy
 * \param out An existing object pointer for output or NULL if none exists.
 * \param err Obit error stack object.
 * \return pointer to the new object.
 */
ObitDConCleanWindow* ObitDConCleanWindowCopy  (ObitDConCleanWindow *in, 
					       ObitDConCleanWindow *out, 
					       ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  gboolean oldExist;
  gchar *outName;
  olong i, j, maxId;
  GSList *glist;
  WindowListElem *elem, *newElem;
  gchar *routine = "ObitDConCleanWindowCopy";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return out;
  g_assert (ObitIsA(in, &myClassInfo));
  if (out) g_assert (ObitIsA(out, &myClassInfo));

  /* Create if it doesn't exist */
  oldExist = out!=NULL;
  if (!oldExist) {
    /* derive object name */
    outName = g_strconcat ("Copy: ",in->name,NULL);
    out = newObitDConCleanWindow(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);
  if (err->error) Obit_traceback_val (err, routine, in->name, NULL);

  /* Free old arrays if any */
  for (i=0; i<out->nfield; i++) {
    if ((out->naxis[i]) && (ObitMemValid (out->naxis[i])))
      out->naxis[i] = ObitMemFree (out->naxis[i]);

    /* Free Window list elements */
    ObitDConCleanWindowDel (in, i+1, -1, err);
    if (err->error) Obit_traceback_val (err, routine, in->name, NULL);
    
    g_slist_free (out->Lists[i]);  /* Head of the list */
  } /* end loop over fields */

  if ((out->Lists) && (ObitMemValid (out->Lists)))
    out->Lists = ObitMemFree (out->Lists);
  if ((out->outWindow) && (ObitMemValid (out->outWindow)))
    out->outWindow = ObitMemFree (out->outWindow);
  if ((out->naxis) && (ObitMemValid (out->naxis))) 
    out->naxis = ObitMemFree (out->naxis);
  if ((out->maxId) && (ObitMemValid (out->maxId)))
    out->maxId = ObitMemFree (out->maxId);

   /*  copy this class */
  out->nfield = in->nfield;
  out->ndim   = in->ndim;

  /* define arrays */
  out->naxis = ObitMemAlloc0Name (out->nfield*sizeof(olong*), "Clean Window Naxis");
  out->maxId = ObitMemAlloc0Name (out->nfield*sizeof(olong),  "Clean Window maxId");
  out->Lists = ObitMemAlloc0Name (out->nfield*sizeof(GSList*), "Clean Window Lists");
  out->outWindow = ObitMemAlloc0Name (out->nfield*sizeof(gpointer), "Clean outer Window");

  for (i=0; i<in->nfield; i++) {
    out->naxis = ObitMemAlloc0 (2*sizeof(olong));
    out->naxis[i][0] = in->naxis[i][0];
    out->naxis[i][1] = in->naxis[i][1];
    maxId            = 0;
    /* Copy any outer window */
    elem = (WindowListElem*)in->outWindow[i];
    out->outWindow[i] = freeWindowListElem (out->outWindow[i]);
    if (elem) 
      out->outWindow[i] = (gpointer) newWindowListElem (1, elem->type, elem->window); 
    else 
      out->outWindow[i] = NULL;
    /* Copy window list, possibly changing IDs */
    glist = in->Lists[i];
    out->Lists[i] = NULL;
    j = 1;
    maxId = 0;
    while (glist!=NULL) {
      elem = glist->data;
      newElem = newWindowListElem (j, elem->type, elem->window);
      out->Lists[i] = g_slist_append (out->Lists[i], newElem);
      maxId = MAX (maxId, newElem->Id);
      j++;
      glist = glist->next;
    }
    out->maxId[i]    = maxId;
 }

  return out;
} /* end ObitDConCleanWindowCopy */

/**
 * Make a copy of a object but do not copy the actual data
 * This is useful to create an DConCleanWindow similar to the input one.
 * \param in  The object to copy
 * \param out An existing object pointer for output, must be defined.
 * \param err Obit error stack object.
 */
void ObitDConCleanWindowClone (ObitDConCleanWindow *in, 
			       ObitDConCleanWindow *out, 
			       ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  olong i;
  gchar *routine = "ObitDConCleanWindowClone";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (ObitIsA(out, &myClassInfo));

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /*  copy this class */
  /* Free old arrays if any */
  for (i=0; i<out->nfield; i++) {
    if ((out->naxis[i]) && (ObitMemValid (out->naxis[i])))
      out->naxis[i] = ObitMemFree (out->naxis[i]);

    /* Free Window list elements */
    ObitDConCleanWindowDel (in, i+1, -1, err);
    if (err->error) Obit_traceback_msg (err, routine, in->name);
    
    g_slist_free (out->Lists[i]);  /* Head of the list */
  } /* end loop over fields */

  if ((out->Lists) && (ObitMemValid (out->Lists)))
    out->Lists = ObitMemFree (out->Lists);
  if ((out->outWindow) && (ObitMemValid (out->outWindow)))
    out->outWindow = ObitMemFree (out->outWindow);
  if ((out->naxis) && (ObitMemValid (out->naxis)))
    out->naxis = ObitMemFree (out->naxis);
  if ((out->maxId) && (ObitMemValid (out->maxId)))
    out->maxId = ObitMemFree (out->maxId);

   /*  copy this class */
  out->nfield = in->nfield;
  out->ndim   = in->ndim;

 /* define arrays */
  out->naxis = ObitMemAlloc0Name (out->nfield*sizeof(olong*), "Clean Window Naxis");
  out->maxId = ObitMemAlloc0Name (out->nfield*sizeof(olong),  "Clean Window maxId");
  out->Lists = ObitMemAlloc0Name (out->nfield*sizeof(GSList*), "Clean Window Lists");
  out->outWindow = ObitMemAlloc0Name (out->nfield*sizeof(gpointer), "Clean outer Window");

  for (i=0; i<in->nfield; i++) {
    out->naxis[i]    = ObitMemAlloc0(2*sizeof(olong));
    out->naxis[i][0] = in->naxis[i][0];
    out->naxis[i][1] = in->naxis[i][1];
    out->maxId[i]    = in->maxId[i];
    out->Lists[i]    = in->Lists[i];
    out->outWindow[i]= in->outWindow[i];
 }

} /* end ObitDConCleanWindowClone */

/**
 * Creates an ObitDConCleanWindow 
 * \param name   An optional name for the object.
 * \param mosaic The image mosaic which this object is to describe.
 * \return the new object.
 */
ObitDConCleanWindow* ObitDConCleanWindowCreate (gchar* name,
						ObitImageMosaic *mosaic, 
						ObitErr *err)
{
  ObitDConCleanWindow* out=NULL;
  olong i;
  gchar *routine = "ObitDConCleanWindowCreate";

  /* error checks */
  if (err->error) return out;
  g_assert (ObitImageMosaicIsA(mosaic));
  if (mosaic->numberImages<=0) {
    Obit_log_error(err, OBIT_Error,"%s ImageMosaic %s has no fields",
                   routine, name);
      return out;
  }

  /* Create basic structure */
  out = newObitDConCleanWindow (name);

  /* Add number of fields */
  out->nfield = mosaic->numberImages;
  out->ndim = 2;

  /* define arrays */
  out->naxis = ObitMemAlloc0Name (out->nfield*sizeof(olong*), "Clean Window Naxis");
  out->maxId = ObitMemAlloc0Name (out->nfield*sizeof(olong),  "Clean Window maxId");
  out->Lists = ObitMemAlloc0Name (out->nfield*sizeof(GSList*), "Clean Window Lists");
  out->outWindow = ObitMemAlloc0Name (out->nfield*sizeof(gpointer), "Clean outer Window");

  for (i=0; i<out->nfield; i++) {
    out->naxis[i] = ObitMemAlloc0(2*sizeof(olong));
    out->naxis[i][0] = mosaic->images[i]->myDesc->inaxes[0];
    out->naxis[i][1] = mosaic->images[i]->myDesc->inaxes[1];
    out->Lists[i]     = NULL;
    out->outWindow[i] = NULL;
    out->maxId[i]     = 0;
  } /* end loop over fields */

  return out;
} /* end ObitDConCleanWindowCreate */

/**
 * Creates an ObitDConCleanWindow with one field
 * \param name   An optional name for the object.
 * \param mosaic The image mosaic which this object is to describe.
 * \return the new object.
 */
ObitDConCleanWindow* 
 ObitDConCleanWindowCreate1 (gchar* name, olong naxis[2], ObitErr *err)
{
  ObitDConCleanWindow* out=NULL;

  /* error checks */
  if (err->error) return out;
  g_assert (naxis!=NULL);

  /* Create basic structure */
  out = newObitDConCleanWindow (name);

  /* Add number of fields */
  out->nfield = 1;
  out->ndim   = 2;

  /* define arrays */
  out->naxis = ObitMemAlloc0Name (out->nfield*sizeof(olong*), "Clean Window Naxis");
  out->maxId = ObitMemAlloc0Name (out->nfield*sizeof(olong),  "Clean Window maxId");
  out->Lists = ObitMemAlloc0Name (out->nfield*sizeof(GSList*), "Clean Window Lists");
  out->outWindow = ObitMemAlloc0Name (out->nfield*sizeof(gpointer), "Clean outer Window");

  /* field dependent stuff */
  out->naxis[0] = ObitMemAlloc0(2*sizeof(olong));
  out->naxis[0][0] = naxis[0];
  out->naxis[0][1] = naxis[1];
  out->Lists[0]     = NULL;
  out->outWindow[0] = NULL;
  out->maxId[0]     = 0;

  return out;
} /* end ObitDConCleanWindowCreate1 */

/**
 * Tell the properties of a window
 * \param in     The Window object
 * \param field  Which field (1-rel) is of interest?
 * \param ID     Window Id, -1=outer
 * \param type   [out] Window type
 * \param window [out] parameters, depends on type
 *                     Pointer into list
 * \param err    Obit error stack object.
 * \return TRUE if window found, else False.
 */
gboolean ObitDConCleanWindowInfo (ObitDConCleanWindow *in, 
			      olong field, olong Id,
			      ObitDConCleanWindowType *type,
			      olong **window, ObitErr *err)
{
  gboolean out=FALSE;
  WindowListElem *elem = NULL;
  gchar *routine = "ObitDConCleanWindowInfo";

  /* error checks */
  if (err->error) return out;
  g_assert (ObitIsA(in, &myClassInfo));
  if ((field<=0) || (field>in->nfield)) {
    Obit_log_error(err, OBIT_Error,"%s field %d out of range 0- %d in %s",
                   routine, field, in->nfield, in->name);
      return out;
  }

  /* Look it up */
  if (Id>0)  /* Inner window */
    elem = ObitDConCleanWindowFind (in->Lists[field-1], Id);
  else  /* Outer window */
    elem = (WindowListElem*)in->outWindow[field-1];
  if (!elem) return FALSE;  /* Not there */
  out = TRUE;
  
  /* return values */
  *type = elem->type;
  *window = elem->window;
  
  return out;
} /* end  ObitDConCleanWindowInfo */

/**
 * Search for a window within toler of a given position.
 * Returns first found meeting criteria
 * \param in     The Window object
 * \param field  Which field (1-rel) is of interest?
 * \param pixel  pixel coordinate (1-rel)
 * \param toler  how close in pixels is required
 * \param which  [out] which part selected, 1=blc or center,
 *                2=trc or radius
 * \param err    Obit error stack object.
 * \return iD if window found, else -1.
 */
olong ObitDConCleanWindowSearch (ObitDConCleanWindow *in, 
				 olong field, olong pixel[2], 
				 olong toler, olong *which, 
				 ObitErr *err)
{
  olong out=-1;
  olong Id, nId, dist, dist1, dist2;
  ofloat dx, dy, arg;
  olong lwhich=-1;
  WindowListElem *elem = NULL;
  gchar *routine = "ObitDConCleanWindowInfo";

  /* error checks */
  if (err->error) return out;
  g_assert (ObitIsA(in, &myClassInfo));
  if ((field<=0) || (field>in->nfield)) {
    Obit_log_error(err, OBIT_Error,"%s field %d out of range 0- %d in %s",
                   routine, field, in->nfield, in->name);
      return out;
  }

  /* Loop over windows */
  nId = in->maxId[field-1];
  for (Id=1; Id<=nId; Id++) {  

    /* Look it up */
    elem = ObitDConCleanWindowFind (in->Lists[field-1], Id);
    if (!elem) continue;

    switch (elem->type) {
    case OBIT_DConCleanWindow_rectangle:
    case OBIT_DConCleanWindow_unrectangle:
      /* try blc and trc*/
      dx = (ofloat)(elem->window[0] - pixel[0]);
      dy = (ofloat)(elem->window[1] - pixel[1]);
      arg = dx*dx + dy*dy;
      arg = MAX (arg, 0.0001);
      dist1 = (olong)(sqrt(arg) + 0.5);
      dx = (ofloat)(elem->window[2] - pixel[0]);
      dy = (ofloat)(elem->window[3] - pixel[1]);
      arg = dx*dx + dy*dy;
      arg = MAX (arg, 0.0001);
      dist2 = (olong)(sqrt(arg) + 0.5);
      /* blc or trc closer? */
      if (dist1<dist2) {
	dist = dist1;
	lwhich = 1;
      } else {
	dist = dist2;
	lwhich = 2;
	}
      break;
    case OBIT_DConCleanWindow_round:
    case OBIT_DConCleanWindow_unround:
      /* try center and along circle*/
      dx = (ofloat)(elem->window[1] - pixel[0]);
      dy = (ofloat)(elem->window[2] - pixel[1]);
      arg = dx*dx + dy*dy;
      arg = MAX (arg, 0.0001);
      dist1 = (olong)(sqrt(arg) + 0.5);
      dist2 = abs (dist1 - elem->window[0]);
      /* center or ring closer? */
      if (dist1<dist2) {
	dist = dist1;
	lwhich = 1;
      } else {
	dist = dist2;
	lwhich = 2;
	}
      break;
    default:
      g_error ("Undefined Clean window type");
      return out;
    }; /* end switch by window type */
    
    /* This one OK? */
    if (dist<=toler) {
      out = Id;
      *which = lwhich;  /* Which part matched */
      break;   /* Done */
    }
    
  } /* end loop over Id */
  return out;
} /* end  ObitDConCleanWindowSearch */

/**
 * Add a window to the list for a given field
 * Window types are:
 * \li OBIT_DConCleanWindow_rectangle
 *     a rectangular box defined by the blc (1-rel) 
 *     (window[0],(window[1]) and
 *     trc corners (window[2],(window[3]) inclusive
 * \li OBIT_DConCleanWindow_round
 *     a round box defined by the radius in cells
 *     (window[0]) and the central pixel (1-rel)
 *     (window[1],(window[2])
 * \li OBIT_DConCleanWindow_unrectangle
 *     a rectangular box defined by the blc (1-rel) 
 *     (window[0],(window[1]) and
 *     trc corners (window[2],(window[3]) inclusive
 *     This region is ALWAYS excluded
 * \li OBIT_DConCleanWindow_unround
 *     a round box defined by the radius in cells
 *     (window[0]) and the central pixel (1-rel)
 *     (window[1],(window[2])
 *     This region is ALWAYS excluded
 * \param in     The Window object
 * \param field  Which field (1-rel) is of interest?
 * \param type   Window type
 * \param window parameters, depends on type
 * \param err    Obit error stack object.
 * \return Id of new window, -1 on failure
 */
olong ObitDConCleanWindowAdd (ObitDConCleanWindow *in, 
			     olong field, 
			     ObitDConCleanWindowType type,
			     olong *window, ObitErr *err)
{
  WindowListElem *elem = NULL;
  olong out = -1;
  gboolean trim;
  gchar *routine = "ObitDConCleanWindowAdd";

  /* error checks */
  if (err->error) return out;
  g_assert (ObitIsA(in, &myClassInfo));
  if ((field<=0) || (field>in->nfield)) {
    Obit_log_error(err, OBIT_Error,"%s field %d out of range 1- %d in %s",
                   routine, field, in->nfield, in->name);
      return out;
  }

  /* Trim to fit if necessary */
  trim = FALSE;
  switch (type) {
  case OBIT_DConCleanWindow_rectangle:
  case OBIT_DConCleanWindow_unrectangle:
    if ((window[0]<1) || (window[0]>in->naxis[field-1][0])) {
      trim = TRUE;
      window[0] = MAX (1, MIN (window[0], in->naxis[field-1][0]));
		       }
    if ((window[1]<1) || (window[1]>in->naxis[field-1][1])) {
      trim = TRUE;
      window[1] = MAX (1, MIN (window[1], in->naxis[field-1][1]));
    }
    if ((window[2]<1) || (window[2]>in->naxis[field-1][0])) {
      trim = TRUE;
      window[2] = MAX (1, MIN (window[2], in->naxis[field-1][0]));
    }
    if ((window[3]<1) || (window[3]>in->naxis[field-1][1])) {
      trim = TRUE;
      window[3] = MAX (1, MIN (window[3], in->naxis[field-1][1]));
    }
    break;
  case OBIT_DConCleanWindow_round:
  case OBIT_DConCleanWindow_unround:
    if ((window[1] - window[0])<0) {
      window[0] = window[1];
      trim = TRUE;
   }
    if ((window[1] + window[0])>in->naxis[field-1][0]) {
      window[0] = MAX (0, in->naxis[field-1][0] - window[1]);
      trim = TRUE;
    }
    if ((window[2] - window[0])<0) {
      window[0] = window[2];
      trim = TRUE;
    }
    if ((window[2] + window[0])>in->naxis[field-1][1]) {
      window[0] =  MAX (0, in->naxis[field-1][1] - window[2]);
      trim = TRUE;
   }
    break;
  default:
    g_error ("Undefined Clean window type");
    return out;
   }; /* end switch by window type */

  /* Warn if trim */
  if (trim) 
    Obit_log_error(err, OBIT_InfoWarn, "%s Trimmed CLEAN window to fit", 
		   routine);

  /* Add it to end of list */
  in->maxId[field-1]++;
  elem =  newWindowListElem (in->maxId[field-1], type, window);
  out = elem->Id;  /* Id number to return */
  in->Lists[field-1] = g_slist_append (in->Lists[field-1], elem);

  return out;
} /* end  ObitDConCleanWindowAdd */

/**
 * Delete a window from the list for a given field
 * Deallocates the list and Window structures
 * \param in     The Window object
 * \param field  Which field (1-rel) is of interest?
 * \param Id     Window Id, -1 => all
 * \param err    Obit error stack object, if NULL, no error reporting
 */
void ObitDConCleanWindowDel (ObitDConCleanWindow *in, 
			     olong field, olong Id, 
			     ObitErr *err)
{
  WindowListElem *elem = NULL;
  GSList *tlist, *xlist;
  gchar *routine = "ObitDConCleanWindowDel";

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));
  if (err) {
    if (err->error) return;
    if ((field<=0) || (field>in->nfield)) {
      Obit_log_error(err, OBIT_Error,"%s field %d out of range 0- %d in %s",
		     routine, field, in->nfield, in->name);
      return;
    }
  } else {  /* harsh test */
    g_assert((field>0) && (field<=in->nfield));
  }

  /* Make sure window list defined - if not just return */
  if (in->Lists[field-1]==NULL) return;

  /* Delete them all? */
  if (Id<0) {
    tlist = in->Lists[field-1];
    while (tlist) {                          /* loop over list */
       elem = (WindowListElem*)tlist->data;
       tlist = g_slist_remove (tlist, elem); /* free list entry */
       freeWindowListElem (elem);            /* free data */
    }
    in->Lists[field-1] = tlist;
    return;
  }

  /* Only one, Look it up */
  elem = ObitDConCleanWindowFind (in->Lists[field-1], Id);
  if (!elem) return;  /* Not there */
  xlist = g_slist_remove (in->Lists[field-1], elem); /* free list entry */
  freeWindowListElem (elem);                         /* free data */
  in->Lists[field-1] = xlist;                        /* new head of list */
} /* end  ObitDConCleanWindowDel */

/**
 * Modify an existing window, add if no window with that ID exists
 * \param in     The Window object
 * \param field  Which field (1-rel) is of interest?
 * \param Id     Window Id, -1=outer
 * \param type   Window type
 * \param window Parameters, depends on type
 * \param err    Obit error stack object.
 */
void ObitDConCleanWindowUpdate (ObitDConCleanWindow *in,  
				olong field, olong Id, 
				ObitDConCleanWindowType type,
				olong *window, ObitErr *err)
{
  WindowListElem *elem = NULL;
  gchar *routine = "ObitDConCleanWindowUpdate";

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  if ((field<=0) || (field>in->nfield)) {
    Obit_log_error(err, OBIT_Error,"%s field %d out of range 0- %d in %s",
                   routine, field, in->nfield, in->name);
      return;
  }
  /* Is this the outer window? */
  if (Id==-1) { /* Yes outer window */
    elem = (WindowListElem*)in->outWindow[field-1];
    /* Update */
    elem->type = type;
    switch (elem->type) {
    case OBIT_DConCleanWindow_rectangle:
    case OBIT_DConCleanWindow_unrectangle:
      /* Copy 4 */
      elem->window[0] = window[0];
      elem->window[1] = window[1];
      elem->window[2] = window[2];
      elem->window[3] = window[3];
      break;
    case OBIT_DConCleanWindow_round:
    case OBIT_DConCleanWindow_unround:
      /* Copy 3 */
      elem->window[0] = window[0];
      elem->window[1] = window[1];
      elem->window[2] = window[2];
      break;
    default:
      g_error ("Undefined Clean window type");
      return;
    }; /* end switch by window type */
    
    return;
  }

  /* Inner window */
  elem = ObitDConCleanWindowFind (in->Lists[field-1], Id);
  if (!elem) {  /* Not there - add with Id */
    elem = newWindowListElem(Id, type, window);
    in->Lists[field-1] = g_slist_append (in->Lists[field-1], elem);
    in->maxId[field-1] = MAX (in->maxId[field-1], Id);
    return;
  }

  /* Update */
  elem->type = type;
  switch (elem->type) {
  case OBIT_DConCleanWindow_rectangle:
  case OBIT_DConCleanWindow_unrectangle:
    /* Copy 4 */
    elem->window[0] = window[0];
    elem->window[1] = window[1];
    elem->window[2] = window[2];
    elem->window[3] = window[3];
    break;
  case OBIT_DConCleanWindow_round:
  case OBIT_DConCleanWindow_unround:
    /* Copy 3 */
    elem->window[0] = window[0];
    elem->window[1] = window[1];
    elem->window[2] = window[2];
    break;
  default:
    g_error ("Undefined Clean window type");
    return;
   }; /* end switch by window type */

} /* end ObitDConCleanWindowUpdate */


/**
 * Set the outer window for a given field
 * Outer windows are used to constrain the automatic setting of windows,
 * the center of an automatically generated window will not be outside
 * of the outer window.
 * Window types are:
 * \li OBIT_DConCleanWindow_rectangle
 *     a rectangular box defined by the blc (1-rel) 
 *     (window[0],(window[1]) and
 *     trc corners (window[2],(window[3]) inclusive
 * \li OBIT_DConCleanWindow_round
 *     a round box defined by the radius in cells
 *     (window[0]) and the central pixel (1-rel)
 *     (window[1],(window[2])
 * \param in     The Window object
 * \param field  Which field (1-rel) is of interest?
 * \param type   Window type
 * \param window parameters, depends on type
 * \param err    Obit error stack object.
 */
void ObitDConCleanWindowOuter (ObitDConCleanWindow *in, olong field, 
			       ObitDConCleanWindowType type,
			       olong *window, ObitErr *err)
{
  gboolean trim;
  gchar *routine = "ObitDConCleanWindowOuter";

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  if ((field<=0) || (field>in->nfield)) {
    Obit_log_error(err, OBIT_Error,"%s: field %d out of range 1- %d in %s",
                   routine, field, in->nfield, in->name);
      return;
  }

  /* Trim to fit if necessary */
  trim = FALSE;
  switch (type) {
  case OBIT_DConCleanWindow_rectangle:
  case OBIT_DConCleanWindow_unrectangle:
    if ((window[0]<1) || (window[0]>in->naxis[field-1][0])) {
      trim = TRUE;
      window[0] = MAX (1, MIN (window[0], in->naxis[field-1][0]));
		       }
    if ((window[1]<1) || (window[1]>in->naxis[field-1][1])) {
      trim = TRUE;
      window[1] = MAX (1, MIN (window[1], in->naxis[field-1][1]));
    }
    if ((window[2]<1) || (window[2]>in->naxis[field-1][0])) {
      trim = TRUE;
      window[2] = MAX (1, MIN (window[2], in->naxis[field-1][0]));
    }
    if ((window[3]<1) || (window[3]>in->naxis[field-1][1])) {
      trim = TRUE;
      window[3] = MAX (1, MIN (window[3], in->naxis[field-1][1]));
    }
    break;
  case OBIT_DConCleanWindow_round:
  case OBIT_DConCleanWindow_unround:
    if ((window[1] - window[0])<0) {
      window[0] = window[1];
      trim = TRUE;
   }
    if ((window[1] + window[0])>in->naxis[field-1][0]) {
      window[0] = MAX (0, in->naxis[field-1][0] - window[1]);
      trim = TRUE;
    }
    if ((window[2] - window[0])<0) {
      window[0] = window[2];
      trim = TRUE;
    }
    if ((window[2] + window[0])>in->naxis[field-1][1]) {
      window[0] =  MAX (0, in->naxis[field-1][1] - window[2]);
      trim = TRUE;
   }
    break;
  default:
    g_error ("Undefined Clean window type");
    return;
   }; /* end switch by window type */

  /* Warn if trim */
  if (trim) 
    Obit_log_error(err, OBIT_InfoWarn, "%s Trimmed CLEAN window to fit", 
		   routine);

  /* Add it to object */
  freeWindowListElem ((WindowListElem*)in->outWindow[field-1]);
  in->outWindow[field-1] = (gpointer)newWindowListElem (1, type, window);

} /* end ObitDConCleanWindowOuter */

/**
 * Are there any valid pixels in this field's image?
 * \param in     The Window object
 * \param field  Which field (1-rel) is of interest?
 * \param err    Obit error stack object.
 * \return TRUE if there are valid pixels, else FALSE
 */
gboolean ObitDConCleanWindowImage (ObitDConCleanWindow *in, 
				   olong field, ObitErr *err)
{
  gboolean out=FALSE;
  gchar *routine = "ObitDConCleanWindowImage";

  /* error checks */
  if (err->error) return out;
  g_assert (ObitIsA(in, &myClassInfo));
  if ((field<=0) || (field>in->nfield)) {
    Obit_log_error(err, OBIT_Error,"%s field %d out of range 0- %d in %s",
                   routine, field, in->nfield, in->name);
      return out;
  }

  /* Always true */
  out = TRUE;

  return out;
} /* end ObitDConCleanWindowImage */


/**
 * Are there any pixels in a specified row in a window but not an unwindow?
 * Default behavior depends on the value of in->autoWindow,
 * Returns TRUE if some pixels are selected.
 * If there are no windows and autoWindow is selected then no pixels
 * are selected.  
 * If there are no windows and autoWindow is not selected,
 * then all pixels are selected.
 * \param in     The Window object
 * \param field  Which field (1-rel) is of interest?
 * \param row    Which row (1-rel)
 * \param mask   [in/out] Mask for pixels in inner boxes,  
 *               TRUE indicates in an inner window but not in unwindow.
 *               If NULL it is created.  Should always be 
 *               Should always be allocated/deallocated using ObitMem
 * \param err    Obit error stack object.
 * \return TRUE if there are selected pixels, else FALSE
 */
gboolean ObitDConCleanWindowRow (ObitDConCleanWindow *in, olong field, 
				 olong row, gboolean **mask, ObitErr *err)
{
  gboolean out=FALSE;
  WindowListElem *elem = NULL;
  GSList *tlist;
  ofloat radius2, rad2Max, y2;
  olong xmax, xmin, ymax, ymin, i;
  gchar *routine = "ObitDConCleanWindowRow";

  /* error checks */
  if (err->error) return out;
  g_assert (ObitIsA(in, &myClassInfo));
  if ((field<=0) || (field>in->nfield)) {
    Obit_log_error(err, OBIT_Error,"%s field %d out of range 1- %d in %s",
                   routine, field, in->nfield, in->name);
      return out;
  }

  if ((row<=0) || (row>in->naxis[field-1][1])) {
    Obit_log_error(err, OBIT_Error,"%s row %d out of range 1- %d in %s",
                   routine, row, in->naxis[field-1][1], in->name);
      return out;
  }

  /* Create bit mask if needed */
  if (*mask==NULL) 
    *mask =  ObitMemAlloc0Name(in->naxis[field-1][0]*sizeof(gboolean),
			       "Pixel Window mask");


  /* Default if no window list */
  out = in->Lists[field-1]!=NULL;
  if (!out) {  /* No list? */
    if (in->autoWindow) {  /* Default = nothing selected */
      for (i=0; i<in->naxis[field-1][0]; i++) (*mask)[i] = FALSE;
      return FALSE;
    } else {               /* Default = all selected */
      for (i=0; i<in->naxis[field-1][0]; i++) (*mask)[i] = TRUE;
      return TRUE;
    }
  } /* end of default with no window list */
  
  /* Init mask */
  for (i=0; i<in->naxis[field-1][0]; i++) (*mask)[i] = FALSE;
  out = FALSE;  /* until proven otherwise */
 
  /* Loop through windows filling mask with windows */
    tlist = in->Lists[field-1];
    while (tlist) { 
       elem = (WindowListElem*)tlist->data;

       /* Process by type */
       switch (elem->type) {
       case OBIT_DConCleanWindow_rectangle:
	 /* (0,1) = blc, (2,3) = trc */
	 /* Can be defined in either way */
	 xmin = MIN(elem->window[0], elem->window[2]);
	 xmax = MAX(elem->window[0], elem->window[2]);
	 ymin = MIN(elem->window[1], elem->window[3]);
	 ymax = MAX(elem->window[1], elem->window[3]);
	 xmax = MAX (1, MIN (xmax, in->naxis[field-1][0]));
	 xmin = MAX (1, MIN (xmin, in->naxis[field-1][0]));
	 /* Is this row inside the window range? */
	 if ((row>=ymin) && (row<=ymax)) {
	   out = TRUE;  /* some valid */
	   for (i=xmin; i<=xmax; i++) (*mask)[i-1] = TRUE;
	 }
	 break;
       case OBIT_DConCleanWindow_unrectangle:
	 break;
       case OBIT_DConCleanWindow_round:
	 /* [0] = radius, (1,2) = center */
	 /* Is this row inside the window range? */
	 if ((row>=(elem->window[2]-elem->window[0])) && 
	     (row<=(elem->window[2]+elem->window[0]))) {
	   /* Maximum radius squared */
	   rad2Max = ((ofloat)elem->window[0]) * ((ofloat)elem->window[0]);
	   y2 = ((ofloat)(elem->window[2]-row) * ((ofloat)(elem->window[2]-row)));
	   /* Check if within radius */
	   xmin = elem->window[1]-elem->window[0];
	   xmax = elem->window[1]+elem->window[0];
	   xmax = MAX (1, MIN (xmax, in->naxis[field-1][0]));
	   xmin = MAX (1, MIN (xmin, in->naxis[field-1][0]));
	   for (i=xmin; i<=xmax; i++) {
	     radius2 = (((ofloat)(elem->window[1]-i)) * 
			((ofloat)(elem->window[1]-i))) + y2;
	     if (radius2<rad2Max) {
	       out = TRUE;
	       (*mask)[i-1] = TRUE;
	     }
	   }
	 }
	 break;
       case OBIT_DConCleanWindow_unround:
	 break;
       default:
	 g_error ("Undefined Clean window type");
       }; /* end switch by window type */
      tlist = tlist->next;  /* Next */
    } /* end unwindow list */
  
  /* Loop through windows unsetting mask with unwindows */
    tlist = in->Lists[field-1];
    while (tlist) { 
       elem = (WindowListElem*)tlist->data;

       /* Process by type */
       switch (elem->type) {
       case OBIT_DConCleanWindow_unrectangle:
	 /* (0,1) = blc, (2,3) = trc */
	 /* Can be defined in either way */
	 xmin = MIN(elem->window[0], elem->window[2]);
	 xmax = MAX(elem->window[0], elem->window[2]);
	 ymin = MIN(elem->window[1], elem->window[3]);
	 ymax = MAX(elem->window[1], elem->window[3]);
	 xmax = MAX (1, MIN (xmax, in->naxis[field-1][0]));
	 xmin = MAX (1, MIN (xmin, in->naxis[field-1][0]));
	 /* Is this row inside the window range? */
	 if ((row>=ymin) && (row<=ymax)) {
	   for (i=xmin; i<=xmax; i++) (*mask)[i-1] = FALSE;
	 }
	 break;
       case OBIT_DConCleanWindow_rectangle:
	 break;
       case OBIT_DConCleanWindow_unround:
	 /* [0] = radius, (1,2) = center */
	 /* Is this row inside the window range? */
	 if ((row>=(elem->window[2]-elem->window[0])) && 
	     (row<=(elem->window[2]+elem->window[0]))) {
	   /* Maximum radius squared */
	   rad2Max = ((ofloat)elem->window[0]) * ((ofloat)elem->window[0]);
	   y2 = ((ofloat)(elem->window[2]-row) * ((ofloat)(elem->window[2]-row)));
	   /* Check if within radius */
	   xmin = elem->window[1]-elem->window[0];
	   xmax = elem->window[1]+elem->window[0];
	   xmax = MAX (1, MIN (xmax, in->naxis[field-1][0]));
	   xmin = MAX (1, MIN (xmin, in->naxis[field-1][0]));
	   for (i=xmin; i<=xmax; i++) {
	     radius2 = (((ofloat)(elem->window[1]-i)) * 
			((ofloat)(elem->window[1]-i))) + y2;
	     if (radius2<rad2Max) {
	       (*mask)[i-1] = FALSE;
	     }
	   }
	 }
	 break;
       case OBIT_DConCleanWindow_round:
	 break;
       default:
	 g_error ("Undefined Clean window type");
       }; /* end switch by window type */
      tlist = tlist->next;  /* Next */
    } /* end unwindow list */
  
  return out;
} /* end ObitDConCleanWindowRow */

/**
 * Are there any valid pixels in a specified row in positive Inner boxes?
 * Any unboxes are ignored
 * If there are no windows and autoWindow is selected, then no pixels
 * are selected.  
 * If there are no windows and autoWindow is not selected,
 * then all pixels are selected.
 * \param in     The Window object
 * \param field  Which field (1-rel) is of interest?
 * \param row    Which row (1-rel)
 * \param mask   [in/out] Mask for pixels in inner boxes,  
 *               TRUE indicates pixel is in an inner window.
 *               If NULL it is created.  Should always be 
 *               allocated/deallocated using ObitMem
 * \param err    Obit error stack object.
 * \return TRUE if there are selected pixels, else FALSE
 */
gboolean ObitDConCleanWindowInnerRow (ObitDConCleanWindow *in, olong field, 
				      olong row, gboolean **mask, ObitErr *err)
{
  gboolean out=FALSE;
  WindowListElem *elem = NULL;
  GSList *tlist;
  ofloat radius2, rad2Max, y2;
  olong xmax, xmin, ymax, ymin, i;
  gchar *routine = "ObitDConCleanWindowInnerRow";

  /* error checks */
  if (err->error) return out;
  g_assert (ObitIsA(in, &myClassInfo));
  if ((field<=0) || (field>in->nfield)) {
    Obit_log_error(err, OBIT_Error,"%s field %d out of range 1- %d in %s",
                   routine, field, in->nfield, in->name);
      return out;
  }

  if ((row<=0) || (row>in->naxis[field-1][1])) {
    Obit_log_error(err, OBIT_Error,"%s row %d out of range 1- %d in %s",
                   routine, row, in->naxis[field-1][1], in->name);
      return out;
  }

  /* Create bit mask if needed */
  if (*mask==NULL) 
    *mask =  ObitMemAlloc0Name(in->naxis[field-1][0]*sizeof(gboolean),
			       "Pixel Window mask");

  /* Default if no window list */
  out = in->Lists[field-1]!=NULL;
  if (!out) {  /* No list? */
    if (in->autoWindow) {  /* Default = nothing selected */
      for (i=0; i<in->naxis[field-1][0]; i++) (*mask)[i] = FALSE;
      return FALSE;
    } else {               /* Default = all selected */
      for (i=0; i<in->naxis[field-1][0]; i++) (*mask)[i] = TRUE;
      return TRUE;
    }
  } /* end of default with no window list */
  
  /* Init mask */
  for (i=0; i<in->naxis[field-1][0]; i++) (*mask)[i] = FALSE;
  out = FALSE;  /* until proven otherwise */
 
  /* Loop through windows filling mask with windows */
  tlist = in->Lists[field-1];
  while (tlist) { 
    elem = (WindowListElem*)tlist->data;
    
    /* Process by type */
    switch (elem->type) {
    case OBIT_DConCleanWindow_rectangle:
      /* (0,1) = blc, (2,3) = trc */
      /* Can be defined in either way */
      xmin = MIN(elem->window[0], elem->window[2]);
      xmax = MAX(elem->window[0], elem->window[2]);
      ymin = MIN(elem->window[1], elem->window[3]);
      ymax = MAX(elem->window[1], elem->window[3]);
      xmax = MAX (1, MIN (xmax, in->naxis[field-1][0]));
      xmin = MAX (1, MIN (xmin, in->naxis[field-1][0]));
      /* Is this row inside the window range? */
      if ((row>=ymin) && (row<=ymax)) {
	out = TRUE;  /* some valid */
	for (i=xmin; i<=xmax; i++) (*mask)[i-1] = TRUE;
      }
      break;
    case OBIT_DConCleanWindow_unrectangle:
      break;
    case OBIT_DConCleanWindow_round:
      /* [0] = radius, (1,2) = center */
      /* Is this row inside the window range? */
      if ((row>=(elem->window[2]-elem->window[0])) && 
	  (row<=(elem->window[2]+elem->window[0]))) {
	/* Maximum radius squared */
	rad2Max = ((ofloat)elem->window[0]) * ((ofloat)elem->window[0]);
	y2 = ((ofloat)(elem->window[2]-row) * ((ofloat)(elem->window[2]-row)));
	/* Check if within radius */
	xmin = elem->window[1]-elem->window[0];
	xmax = elem->window[1]+elem->window[0];
	xmax = MAX (1, MIN (xmax, in->naxis[field-1][0]));
	xmin = MAX (1, MIN (xmin, in->naxis[field-1][0]));
	for (i=xmin; i<=xmax; i++) {
	  radius2 = (((ofloat)(elem->window[1]-i)) * 
		     ((ofloat)(elem->window[1]-i))) + y2;
	  if (radius2<rad2Max) {
	    out = TRUE;
	    (*mask)[i-1] = TRUE;
	  }
	}
      }
      break;
    case OBIT_DConCleanWindow_unround:
      break;
    default:
      g_error ("Undefined Clean window type");
    }; /* end switch by window type */
    tlist = tlist->next;  /* Next */
  } /* end window list */
  
  return out;
} /* end ObitDConCleanWindowInnerRow */

/**
 * Are there any pixels in a specified row within unboxes?
 * If there are no unwindows the mask is all FALSE and FALSE is returned.
 * \param in     The Window object
 * \param field  Which field (1-rel) is of interest?
 * \param row    Which row (1-rel)
 * \param mask   [in/out] Mask for pixels in unboxes,  
 *               TRUE indicates pixel is in an unwindow.
 *               If NULL it is created.  Should always be 
 *               allocated/deallocated using ObitMem
 * \param err    Obit error stack object.
 * \return TRUE if there are selected pixels, else FALSE
 */
gboolean ObitDConCleanWindowUnrow (ObitDConCleanWindow *in, olong field, 
				 olong row, gboolean **mask, ObitErr *err)
{
  gboolean out=FALSE;
  WindowListElem *elem = NULL;
  GSList *tlist;
  ofloat radius2, rad2Max, y2;
  olong xmax, xmin, ymax, ymin, i;
  gchar *routine = "ObitDConCleanWindowunrow";

  /* error checks */
  if (err->error) return out;
  g_assert (ObitIsA(in, &myClassInfo));
  if ((field<=0) || (field>in->nfield)) {
    Obit_log_error(err, OBIT_Error,"%s field %d out of range 1- %d in %s",
                   routine, field, in->nfield, in->name);
      return out;
  }

  if ((row<=0) || (row>in->naxis[field-1][1])) {
    Obit_log_error(err, OBIT_Error,"%s row %d out of range 1- %d in %s",
                   routine, row, in->naxis[field-1][1], in->name);
      return out;
  }

  /* Create bit mask if needed */
  if (*mask==NULL) 
    *mask =  ObitMemAlloc0Name(in->naxis[field-1][0]*sizeof(gboolean),
			       "Pixel Window mask");

  /* Default if no window list */
  out = in->Lists[field-1]!=NULL;
  if (!out) {  /* No list? Nothing selected */
    for (i=0; i<in->naxis[field-1][0]; i++) (*mask)[i] = FALSE;
    return FALSE;
  } /* end of default with no window list */
  
  /* Init mask */
  for (i=0; i<in->naxis[field-1][0]; i++) (*mask)[i] = FALSE;
  out = FALSE;  /* until proven otherwise */
 
  /* Loop through unwindows filling mask with windows */
  tlist = in->Lists[field-1];
  while (tlist) { 
    elem = (WindowListElem*)tlist->data;
    
    /* Process by type */
    switch (elem->type) {
    case OBIT_DConCleanWindow_unrectangle:
      /* (0,1) = blc, (2,3) = trc */
      /* Can be defined in either way */
      xmin = MIN(elem->window[0], elem->window[2]);
      xmax = MAX(elem->window[0], elem->window[2]);
      ymin = MIN(elem->window[1], elem->window[3]);
      ymax = MAX(elem->window[1], elem->window[3]);
      xmax = MAX (1, MIN (xmax, in->naxis[field-1][0]));
      xmin = MAX (1, MIN (xmin, in->naxis[field-1][0]));
      /* Is this row inside the window range? */
      if ((row>=ymin) && (row<=ymax)) {
	out = TRUE;  /* some valid */
	for (i=xmin; i<=xmax; i++) (*mask)[i-1] = TRUE;
      }
      break;
    case OBIT_DConCleanWindow_rectangle:
      break;
    case OBIT_DConCleanWindow_unround:
      /* [0] = radius, (1,2) = center */
      /* Is this row inside the window range? */
      if ((row>=(elem->window[2]-elem->window[0])) && 
	  (row<=(elem->window[2]+elem->window[0]))) {
	/* Maximum radius squared */
	rad2Max = ((ofloat)elem->window[0]) * ((ofloat)elem->window[0]);
	y2 = ((ofloat)(elem->window[2]-row) * ((ofloat)(elem->window[2]-row)));
	/* Check if within radius */
	xmin = elem->window[1]-elem->window[0];
	xmax = elem->window[1]+elem->window[0];
	xmax = MAX (1, MIN (xmax, in->naxis[field-1][0]));
	xmin = MAX (1, MIN (xmin, in->naxis[field-1][0]));
	for (i=xmin; i<=xmax; i++) {
	  radius2 = (((ofloat)(elem->window[1]-i)) * 
		     ((ofloat)(elem->window[1]-i))) + y2;
	  if (radius2<rad2Max) {
	    out = TRUE;
	    (*mask)[i-1] = TRUE;
	  }
	}
      }
      break;
    case OBIT_DConCleanWindow_round:
      break;
    default:
      g_error ("Undefined Clean window type");
    }; /* end switch by window type */
    tlist = tlist->next;  /* Next */
  } /* end unwindow list */
  
  return out;
} /* end ObitDConCleanWindowUnrow */

/**
 * Are there any valid pixels from outer window in a specified row?
 * If there are no windows defined, all pixels are allowed except outer 5.
 * Unwindows not supported.
 * \param in     The Window object
 * \param field  Which field (1-rel) is of interest?
 * \param row    Which row (1-rel)
 * \param mask   [in/out] Mask for valid pixels,  If NULL, it is created.
 *               Should always be allocated/deallocated using ObitMem
 * \param err    Obit error stack object.
 * \return TRUE if there are valid pixels, else FALSE
 */
gboolean 
ObitDConCleanWindowOuterRow (ObitDConCleanWindow *in, olong field, 
			     olong row, gboolean **mask, ObitErr *err)
{
  gboolean out=FALSE;
  gboolean fill;
  WindowListElem *elem = NULL;
  ofloat radius2, rad2Max, y2;
  olong xmax, xmin, ymax, ymin, i;
  gchar *routine = "ObitDConCleanWindowOuterRow";

  /* error checks */
  if (err->error) return out;
  g_assert (ObitIsA(in, &myClassInfo));
  if ((field<=0) || (field>in->nfield)) {
    Obit_log_error(err, OBIT_Error,"%s field %d out of range 1- %d in %s",
                   routine, field, in->nfield, in->name);
      return out;
  }

  if ((row<=0) || (row>in->naxis[field-1][1])) {
    Obit_log_error(err, OBIT_Error,"%s row %d out of range 1- %d in %s",
                   routine, row, in->naxis[field-1][1], in->name);
      return out;
  }

  /* Create bit mask if needed */
  if (*mask==NULL) 
    *mask =  ObitMemAlloc0Name(in->naxis[field-1][0]*sizeof(gboolean),
			       "Pixel Window mask");

  /* If list empty use default = all OK except the outer 5 pixels 
     else check further */
  out = in->outWindow[field-1]!=NULL;
  if (!out) {
    if ((row<=5) || (row>(in->naxis[field-1][1]-5))) fill = FALSE;
    else fill = TRUE;
    for (i=0; i<in->naxis[field-1][0]; i++) (*mask)[i] = fill;
    for (i=0; i<5; i++) (*mask)[i] = FALSE;
    for (i=in->naxis[field-1][0]-5; i<in->naxis[field-1][0]; i++) (*mask)[i] = FALSE;
    return fill;
  }

  /* Init mask */
  for (i=0; i<in->naxis[field-1][0]; i++) (*mask)[i] = FALSE;
  out = FALSE;  /* until proven otherwise */
 
  /* Get outer window */
  elem = (WindowListElem*)in->outWindow[field-1];

  /* Process by type */
  switch (elem->type) {
  case OBIT_DConCleanWindow_rectangle:
    /* (0,1) = blc, (2,3) = trc */
    /* Can be defined in either way */
    xmin = MIN(elem->window[0], elem->window[2]);
    xmax = MAX(elem->window[0], elem->window[2]);
    ymin = MIN(elem->window[1], elem->window[3]);
    ymax = MAX(elem->window[1], elem->window[3]);
    xmax = MAX (1, MIN (xmax, in->naxis[field-1][0]));
    xmin = MAX (1, MIN (xmin, in->naxis[field-1][0]));
    /* Is this row inside the window range? */
    if ((row>=ymin) && (row<=ymax)) {
      out = TRUE;  /* some valid */
      for (i=xmin; i<=xmax; i++) (*mask)[i-1] = TRUE;
    }
    break;
  case OBIT_DConCleanWindow_unrectangle:
    break;
  case OBIT_DConCleanWindow_round:
    /* [0] = radius, (1,2) = center */
    /* Is this row inside the window range? */
    if ((row>=(elem->window[2]-elem->window[0])) && 
	(row<=(elem->window[2]+elem->window[0]))) {
      /* Maximum radius squared */
      rad2Max = ((ofloat)elem->window[0]) * ((ofloat)elem->window[0]);
      y2 = ((ofloat)(elem->window[2]-row) * ((ofloat)(elem->window[2]-row)));
      /* Check if within radius */
      xmin = elem->window[1]-elem->window[0];
      xmax = elem->window[1]+elem->window[0];
      xmax = MAX (1, MIN (xmax, in->naxis[field-1][0]));
      xmin = MAX (1, MIN (xmin, in->naxis[field-1][0]));
      for (i=xmin; i<=xmax; i++) {
	radius2 = (((ofloat)(elem->window[1]-i)) * 
		   ((ofloat)(elem->window[1]-i))) + y2;
	if (radius2<rad2Max) {
	  out = TRUE;
	  (*mask)[i-1] = TRUE;
	}
      }
    }
    break;
  case OBIT_DConCleanWindow_unround:
    break;
  default:
    g_error ("Undefined Clean window type");
  }; /* end switch by window type */
  
  return out;
} /* end ObitDConCleanWindowOuterRow */

/**
 *  What is the maximum extent in either x or y covered?
 * \param in     The Window object
 * \param field  Which field (1-rel) is of interest?
 * \param err    Obit error stack object.
 * \return number of pixels in extent
 */
olong ObitDConCleanWindowSize (ObitDConCleanWindow *in, olong field, 
				 ObitErr *err)
{
  olong out=0;
  olong xmax=-1000000, xmin=1000000, ymax=-1000000, ymin=1000000;
    olong i, deltax, deltay, axsize;
  WindowListElem* win;
  gchar *routine = "ObitDConCleanWindowSize";

  /* error checks */
  if (err->error) return out;
  g_assert (ObitIsA(in, &myClassInfo));
  if ((field<=0) || (field>in->nfield)) {
    Obit_log_error(err, OBIT_Error,"%s field %d out of range 1- %d in %s",
                   routine, field, in->nfield, in->name);
      return out;
  }

  /* If list empty return maximum dimension/2, else check further */
  axsize = MIN (in->naxis[field-1][0], in->naxis[field-1][1]);
  out    = axsize;
  if (in->Lists[field-1]==NULL) return out/2;

  /* Loop over entries comparing finding extrema */
  for (i = 1; i<=in->maxId[field-1]; i++) {
    win = ObitDConCleanWindowFind (in->Lists[field-1], i);
    if (win) { /* this one defined? */
      /* extrema by type */
      switch (win->type) {
      case OBIT_DConCleanWindow_rectangle:
	/* Can be defined in either way */
	xmin = MIN (xmin, MIN(win->window[0], win->window[2]));
	xmax = MAX (xmax, MAX(win->window[0], win->window[2]));
	ymin = MIN (ymin, MIN(win->window[1], win->window[3]));
	ymax = MAX (ymax, MAX(win->window[2], win->window[3]));
	break;
      case OBIT_DConCleanWindow_unrectangle:  /* unwindows ignored */
	break;
      case OBIT_DConCleanWindow_round:
	xmin = MIN (xmin, win->window[1] - win->window[0]);
	xmax = MAX (xmax, win->window[1] + win->window[0]);
	ymin = MIN (ymin, win->window[2] - win->window[0]);
	ymax = MAX (ymax, win->window[2] + win->window[0]);
	break;
      case OBIT_DConCleanWindow_unround:  /* unwindows ignored */
	break;
      default:
	g_error ("Undefined Clean window type");
      }; /* end switch by window type */
    }
  }

  /* Find range - beam patch is half width */
  deltax = abs (xmax - xmin);
  deltax = MIN (deltax, axsize);
  deltay = abs (ymax - ymin);
  deltay = MIN (deltay, axsize);
  out = MAX (out, MAX (deltax, deltay));
  out = MIN (out, MAX (deltax, deltay));

  return out/2;
} /* end ObitDConCleanWindowSize */

/**
 *  How many pixels are selected in a plane?
 * \param in     The Window object
 * \param field  Which field (1-rel) is of interest?
 * \param err    Obit error stack object.
 * \return number of valid pixels
 */
olong ObitDConCleanWindowCount (ObitDConCleanWindow *in, olong field, 
				ObitErr *err)
{
  olong out=0;
  olong i, j;
  gboolean *mask=NULL;
  gchar *routine = "ObitDConCleanWindowCount";

  /* error checks */
  if (err->error) return out;
  g_assert (ObitIsA(in, &myClassInfo));
  if ((field<=0) || (field>in->nfield)) {
    Obit_log_error(err, OBIT_Error,"%s field %d out of range 1- %d in %s",
                   routine, field, in->nfield, in->name);
      return out;
  }

  /* If list empty return size of image */
  out = in->naxis[field-1][0] * in->naxis[field-1][1];
  if (in->Lists[field-1]==NULL) return out;
  
  /* Loop over rows counting */
  out = 0;
  for (i=0; i<in->naxis[field-1][1]; i++) {
    if (ObitDConCleanWindowRow(in, field, i+1, &mask, err)) {
      for (j=0; j<in->naxis[field-1][0]; j++) if (mask[j]) out++;
    }
    if (err->error) Obit_traceback_val (err, routine, in->name, out);
  }

  mask = ObitMemFree(mask);  /* Free mask memory */

  return out;
} /* end ObitDConCleanWindowCount */

/**
 * Add autoWindow box to a window if appropriate.
 * If the peak in the image is > n*RMS and occurs inside of the outer window but 
 * outside the previous inner window, a new round box is added at that position.
 * n=4 for small boxes, 3 large.
 * The added window is round and of a size where the structure function 
 * about the center drops to 10% or 3 sigma whichever is less (max=20)
 * \param in         The Window object
 * \param field      Which field (1-rel) is of interest?
 * \param image      pixel array, will be returned blanked outside the outer
 *                   window and inside the inner window
 * \param doAbs      If TRUE look for max. abs., otherwise max.
 * \param PeakIn     [out] Peak value inside of outer window
 * \param PeakInPos  [out] pixel position (1-rel) of PeakIn
 * \param PeakOut    [out] Peak value outside of inner window but within outer 
 *
 * \param RMS        [out] RMS within outer Window
 * \param err        Obit error stack object.
 * \return TRUE if PeakIn occurs outside of the current inner window
 */
gboolean ObitDConCleanWindowAutoWindow (ObitDConCleanWindow *in, 
					olong field, ObitFArray *image,
					gboolean doAbs,
					ofloat *PeakIn, olong *PeakInPos,
					ofloat *PeakOut, ofloat *RMS,
					ObitErr *err)
{
  gboolean addWin = FALSE;
  gboolean noWin = FALSE, *mask=NULL;
  ObitFArray *tmpImage=NULL;
  olong ix, iy, nx, ny, pos[2];
  olong window[4];
  ofloat *data, minFlux, fblank =  ObitMagicF();
  gchar *routine = "ObitDConCleanWindowAutoWindow";

  /* error checks */
  if (err->error) return addWin;
  g_assert (ObitIsA(in, &myClassInfo));
  if ((field<=0) || (field>in->nfield)) {
    Obit_log_error(err, OBIT_Error,"%s field %d out of range 1- %d in %s",
                   routine, field, in->nfield, in->name);
      return addWin;
  }

  /* Image size */
  nx = in->naxis[field-1][0];
  ny = in->naxis[field-1][1];
  /* Check */
  Obit_retval_if_fail (((nx==image->naxis[0]) && (ny==image->naxis[1])),
		       err, addWin, "%s: Window and image different sizes ( %d  %d) ( %d  %d)",
		       routine, nx, ny, image->naxis[0], image->naxis[1]);
  
  /* Copy of image for statistics  */
  tmpImage = ObitFArrayCopy(image, tmpImage, err);
  if (err->error) Obit_traceback_val (err, routine, image->name, addWin);

  /* blank tmpImage outside of outer window or in unboxes */
  for (iy=0; iy<ny; iy++) {
    /* pointer to data */
    pos[0] = 0; pos[1] = iy;
    data = ObitFArrayIndex(tmpImage, pos);
    /* Get and apply window mask */
    if (ObitDConCleanWindowOuterRow(in, field, iy+1, &mask, err)) {
      for (ix=0; ix<nx; ix++) if (!mask[ix]) data[ix] = fblank;
    } else { /* nothing valid - blank all */
      for (ix=0; ix<nx; ix++) data[ix] = fblank;
    }
    /* Blank in any unboxes as well */
    if (ObitDConCleanWindowUnrow(in, field, iy+1, &mask, err)) {
      for (ix=0; ix<nx; ix++) if (mask[ix]) data[ix] = fblank;
    }
    if (err->error) goto clean;
  } /* end loop blanking array */
  
  /* Find RMS, peak and pos in tmpImage = RMS, PeakIn, PeakInPos */
  *RMS = ObitFArrayRMS (tmpImage);
  if (doAbs) 
    *PeakIn = ObitFArrayMaxAbs (tmpImage, PeakInPos);
  else
    *PeakIn = ObitFArrayMax (tmpImage, PeakInPos);

  /* Blank inside the inner window  - if there is one */
  if (in->Lists[field-1]) {
    for (iy=0; iy<ny; iy++) {
      /* pointer to data */
      pos[0] = 0; pos[1] = iy;
      data = ObitFArrayIndex(tmpImage, pos);
      /* Get, invert, and apply window mask */
      if (ObitDConCleanWindowInnerRow(in, field, iy+1, &mask, err)) {
	for (ix=0; ix<nx; ix++) if (mask[ix]) data[ix] = fblank;
      }
    } /* end loop blanking array */
    if (err->error) goto clean;
  } else { /* No previous windows - add one */
    addWin = TRUE;
  }

  /* if PeakInPos not blanked and > 4 RMS  addWin = TRUE; */
  data = ObitFArrayIndex(tmpImage, PeakInPos);
  addWin = addWin || ((*data)!=fblank);
  /* Reduce threshold for more extended regions */
  window[0] = GetWindowSize(image, PeakInPos, *RMS);
  if (window[0]<5)
    minFlux = 4.0*(*RMS);
  else
    minFlux = 3.0*(*RMS);
  /* Window not set because peak too close to noise? */
  noWin  = (fabs(*data) < minFlux); 
  addWin = addWin && (!noWin);

  /* Add new clean box? */
  if (addWin) {
    window[1] = PeakInPos[0]+1;
    window[2] = PeakInPos[1]+1;
    ObitDConCleanWindowAdd (in, field, OBIT_DConCleanWindow_round, 
			    window, err);
    if (err->error) goto clean;

    /* inform user */
    Obit_log_error(err, OBIT_InfoErr,"Added round box radius= %d ( %d, %d) to field  %d",
		   window[0], window[1], window[2], field);


    /* Need stats - blank with new window */
    for (iy=0; iy<ny; iy++) {
      /* pointer to data */
      pos[0] = 0; pos[1] = iy;
      data = ObitFArrayIndex(tmpImage, pos);
      /* Get, invert, and apply window mask */
      if (ObitDConCleanWindowRow(in, field, iy+1, &mask, err)) {
	for (ix=0; ix<nx; ix++) if (mask[ix]) data[ix] = fblank;
      }
    } /* end loop blanking array */
    if (err->error) goto clean;
     /* end add new box */
  } 
  
  /* find peak PeakOut - this is used to set min CLEAN; 
     make small if peak is in the noise */
  if (doAbs)  *PeakOut = ObitFArrayMaxAbs (tmpImage, pos);
  else *PeakOut = ObitFArrayMax (tmpImage, pos);
  if (fabs(*PeakOut)<=minFlux)  *PeakOut = 0.0;
  /* If no window set because peak too close to noise, set to zero */
  if (noWin)  *PeakOut = 0.0;
  /* Trap problem */
  if (fabs(*PeakOut)>1.0e20) *PeakOut = *PeakIn; 

  /* Cleanup */
 clean:
  if (mask) ObitMemFree (mask);
  tmpImage = ObitFArrayUnref(tmpImage);
  if (err->error) Obit_traceback_val (err, routine, image->name, addWin);

  return addWin;
} /* end ObitDConCleanWindowAutoWindow */

/**
 * Replace all windows in a given field with those from another object
 * Outer window not replaced on out if not defined in in.
 * Naxis must correspond between windows
 * \param in         Source Window object
 * \param ifield     Which field (1-rel) is of interest in in?
 * \param out        Source Window object
 * \param ofield     Which field (1-rel) is of interest in out?
 * \param err        Obit error stack object.
 */
void 
ObitDConCleanWindowReplaceField (ObitDConCleanWindow *in,  olong ifield, 
				 ObitDConCleanWindow *out, olong ofield,
				 ObitErr *err)
{
  WindowListElem *elem,*newElem ;
  GSList *glist;
  olong j, maxId;
  gchar *routine = "ObitDConCleanWindowReplaceField";

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  if ((ifield<=0) || (ifield>in->nfield)) {
    Obit_log_error(err, OBIT_Error,"%s ifield %d out of range 1- %d in %s",
                   routine, ifield, in->nfield, in->name);
      return;
  }
  g_assert (ObitIsA(out, &myClassInfo));
  if ((ofield<=0) || (ofield>out->nfield)) {
    Obit_log_error(err, OBIT_Error,"%s ofield %d out of range 1- %d in %s",
                   routine, ofield, out->nfield, in->name);
      return;
  }
  /* Check field size - must be the same */
  if ((in->naxis[ifield-1][0]!=out->naxis[ofield-1][0]) ||
      (in->naxis[ifield-1][1]!=out->naxis[ofield-1][1])) {
    Obit_log_error(err, OBIT_Error,"%s field sizes incompatible [ %d, %d], [ %d, %d]",
                   routine, in->naxis[ifield-1][0], in->naxis[ifield-1][1],
		   out->naxis[ofield-1][0], out->naxis[ofield-1][1]);
      return;
  }

  /* Replace outer window? */
  if (in->outWindow[ifield-1]!=NULL) {
    elem = (WindowListElem*)in->outWindow[ifield-1];
    if (elem!=NULL) {
      out->outWindow[ofield-1] = freeWindowListElem (out->outWindow[ofield-1]);
      out->outWindow[ofield-1] = 
	(gpointer) newWindowListElem (1, elem->type, elem->window); 
    }
  } /* end replace outer window */

  /* Replace windows */
  ObitDConCleanWindowDel (out, ofield, -1, NULL);
 
  /* Copy window list, possibly changing IDs */
  glist = in->Lists[ifield-1];
  maxId = 0;
  out->maxId[ofield-1] = maxId;
  j = 1;
  while (glist!=NULL) {
    elem = glist->data;
    newElem = newWindowListElem (j, elem->type, elem->window);
    out->Lists[ofield-1] = g_slist_append (out->Lists[ofield-1], newElem);
    maxId = MAX (maxId, newElem->Id);
    j++;
    glist = glist->next;
  }
  out->maxId[ofield-1] = maxId;

  return;
} /* end ObitDConCleanWindowReplaceField */

/**
 * Add a field with an empty window to an ObitDConCleanWindow
 * \param in         Source Window object
 * \param inaxes     Dimension of image
 * \param err        Obit error stack object.
 * \return added field number (1-rel)
 */
olong 
ObitDConCleanWindowAddField (ObitDConCleanWindow *in, 
			     olong inaxes[2], ObitErr *err)

{
  olong newField = 0;
  olong i, oldField;
  olong *ltemp, **lltemp;
  GSList **tlist;
  gpointer *tpointer;
  /*gchar *routine = "ObitDConCleanWindowReplaceField";*/

  /* error checks */
  if (err->error) return newField;

  /* field to add */
  oldField = in->nfield;
  newField = oldField+1;
  in->nfield = newField;

 /* Resize arrays - naxis */
  lltemp = ObitMemAlloc0Name (newField*sizeof(olong*), "Clean Window Naxis");
  for (i=0; i<oldField; i++) {
    lltemp[i] = ObitMemAlloc0 (2*sizeof(olong));
    lltemp[i][0] = in->naxis[i][0];
    lltemp[i][1] = in->naxis[i][1];
    in->naxis[i] = ObitMemFree(in->naxis[i]);
  }
  lltemp[newField-1] = ObitMemAlloc0 (2*sizeof(olong));
  lltemp[newField-1][0] = inaxes[0];
  lltemp[newField-1][1] = inaxes[1];
  in->naxis = ObitMemFree(in->naxis);
  in->naxis = lltemp;

  /* maxId */
  ltemp = ObitMemAlloc0Name (newField*sizeof(olong),  "Clean Window maxId");
  for (i=0; i<oldField; i++) ltemp[i] = in->maxId[i]; 
  ltemp[newField-1] = 0;
  in->maxId = ObitMemFree(in->maxId);
  in->maxId = ltemp;
 
  /* Window Lists */
  tlist = ObitMemAlloc0Name (newField*sizeof(GSList*), "Clean Window Lists");
  for (i=0; i<oldField; i++) tlist[i] = in->Lists[i]; 
  tlist[newField-1] = NULL;
  in->Lists = ObitMemFree(in->Lists);
  in->Lists = tlist;
  
  /* outer Window */
  tpointer = ObitMemAlloc0Name (newField*sizeof(gpointer), "Clean outer Window");
  for (i=0; i<oldField; i++) tpointer[i] = in->outWindow[i]; 
  tpointer[newField-1] = NULL;
  in->outWindow = ObitMemFree(in->outWindow);
  in->outWindow = tpointer;

  return newField;
} /* end ObitDConCleanWindowAddField */

/**
 * Determine statistical measures of the first plane of an image with 
 * corresponding window
 * \param in         The Window object
 * \param field      Which field (1-rel) is of interest?
 * \param Image      pixel array, will be returned blanked outside the outer
 *                   window and inside the inner window
 * \param doAbs      If TRUE look for max. abs., otherwise max.
 * \param PeakIn     [out] Peak value inside of outer window
 * \param PeakInPos  [out] pixel position (1-rel) of PeakIn
 * \param PeakOut    [out] Peak value outside of inner window but within outer 
 *
 * \param RMS        [out] RMS within outer Window
 * \param err        Obit error stack object.
 */
void ObitDConCleanWindowStats (ObitDConCleanWindow *in, 
			       olong field, ObitFArray *image,
			       gboolean doAbs,
			       ofloat *PeakIn, olong *PeakInPos,
			       ofloat *PeakOut, ofloat *RMS,
			       ObitErr *err)
{
  gboolean *mask=NULL;
  ObitFArray *tmpImage=NULL;
  olong ix, iy, nx, ny, pos[2];
  ofloat *data, fblank = ObitMagicF();
  gchar *routine = "ObitDConCleanWindowStats";

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  if ((field<=0) || (field>in->nfield)) {
    Obit_log_error(err, OBIT_Error,"%s field %d out of range 1- %d in %s",
                   routine, field, in->nfield, in->name);
      return;
  }

  /* Image size */
  nx = in->naxis[field-1][0];
  ny = in->naxis[field-1][1];
  /* Check */
  Obit_return_if_fail (((nx==image->naxis[0]) && (ny==image->naxis[1])),
		       err, "%s: Window and image different sizes ( %d  %d) ( %d  %d)",
		       routine, nx, ny, image->naxis[0], image->naxis[1]);
  
  /* Copy of image for statistics  */
  tmpImage = ObitFArrayCopy(image, tmpImage, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* blank tmpImage outside of outer window or in unboxes */
  for (iy=0; iy<ny; iy++) {
    /* pointer to data */
    pos[0] = 0; pos[1] = iy;
    data = ObitFArrayIndex(tmpImage, pos);
    /* Get and apply window mask */
    if (ObitDConCleanWindowOuterRow(in, field, iy+1, &mask, err)) {
      for (ix=0; ix<nx; ix++) if (!mask[ix]) data[ix] = fblank;
    } else { /* nothing valid - blank all */
      for (ix=0; ix<nx; ix++) data[ix] = fblank;
    }
    /* Blank in any unboxes as well */
    if (ObitDConCleanWindowUnrow(in, field, iy+1, &mask, err)) {
      for (ix=0; ix<nx; ix++) if (mask[ix]) data[ix] = fblank;
    }
    if (err->error) goto clean;
  } /* end loop blanking array */
  
  /* Find RMS, peak and pos in tmpImage = RMS, PeakIn, PeakInPos */
  *RMS = ObitFArrayRMS (tmpImage);
  if (doAbs) 
    *PeakIn = ObitFArrayMaxAbs (tmpImage, PeakInPos);
  else
    *PeakIn = ObitFArrayMax (tmpImage, PeakInPos);
  *PeakOut = *PeakIn;  /* In case no inner window */

  /* Blank inside the inner window  - if there is one */
  if (in->Lists[field-1]) {
    for (iy=0; iy<ny; iy++) {
      /* pointer to data */
      pos[0] = 0; pos[1] = iy;
      data = ObitFArrayIndex(tmpImage, pos);
      /* Get, invert, and apply window mask */
      if (ObitDConCleanWindowInnerRow(in, field, iy+1, &mask, err)) {
	for (ix=0; ix<nx; ix++) 
	  if (mask[ix]) 
	    data[ix] = fblank;
      }
    } /* end loop blanking array */
    if (err->error) goto clean;

    /* find peak PeakOut in outer window but outside inner window */
    if (doAbs)  *PeakOut = ObitFArrayMaxAbs (tmpImage, pos);
    else *PeakOut = ObitFArrayMax (tmpImage, pos);
    /* Trap problem */
    if (fabs(*PeakOut)>1.0e20) *PeakOut = *PeakIn; 
  } /* end if inner window */


  /* Cleanup */
 clean:
  if (mask) ObitMemFree (mask);
  tmpImage = ObitFArrayUnref(tmpImage);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  return;
} /* end ObitDConCleanWindowStats */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitDConCleanWindowClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitDConCleanWindowClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitDConCleanWindowClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitDConCleanWindowClassInfoDefFn (gpointer inClass)
{
  ObitDConCleanWindowClassInfo *theClass = (ObitDConCleanWindowClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitDConCleanWindowClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitDConCleanWindowClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitDConCleanWindowGetClass;
  theClass->newObit       = (newObitFP)newObitDConCleanWindow;
  theClass->ObitCopy      = (ObitCopyFP)ObitDConCleanWindowCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitDConCleanWindowClear;
  theClass->ObitInit      = (ObitInitFP)ObitDConCleanWindowInit;
  theClass->ObitDConCleanWindowCreate = (ObitDConCleanWindowCreateFP)ObitDConCleanWindowCreate;
  theClass->ObitDConCleanWindowCreate1= (ObitDConCleanWindowCreate1FP)ObitDConCleanWindowCreate1;
  theClass->ObitDConCleanWindowInfo   = (ObitDConCleanWindowInfoFP)ObitDConCleanWindowInfo;
  theClass->ObitDConCleanWindowSearch = (ObitDConCleanWindowSearchFP)ObitDConCleanWindowSearch;
  theClass->ObitDConCleanWindowAdd    = (ObitDConCleanWindowAddFP)ObitDConCleanWindowAdd;
  theClass->ObitDConCleanWindowDel    = (ObitDConCleanWindowDelFP)ObitDConCleanWindowDel;
  theClass->ObitDConCleanWindowUpdate = (ObitDConCleanWindowUpdateFP)ObitDConCleanWindowUpdate;
  theClass->ObitDConCleanWindowOuter  = (ObitDConCleanWindowOuterFP)ObitDConCleanWindowOuter;
  theClass->ObitDConCleanWindowImage  = (ObitDConCleanWindowImageFP)ObitDConCleanWindowImage;
  theClass->ObitDConCleanWindowRow    = (ObitDConCleanWindowRowFP)ObitDConCleanWindowRow;
  theClass->ObitDConCleanWindowOuterRow = 
    (ObitDConCleanWindowOuterRowFP)ObitDConCleanWindowOuterRow;
  theClass->ObitDConCleanWindowInnerRow = 
    (ObitDConCleanWindowInnerRowFP)ObitDConCleanWindowInnerRow;
  theClass->ObitDConCleanWindowUnrow = 
    (ObitDConCleanWindowUnrowFP)ObitDConCleanWindowUnrow;
  theClass->ObitDConCleanWindowSize   = (ObitDConCleanWindowSizeFP)ObitDConCleanWindowSize;
  theClass->ObitDConCleanWindowCount  = (ObitDConCleanWindowCountFP)ObitDConCleanWindowCount;
  theClass->ObitDConCleanWindowAutoWindow  = 
    (ObitDConCleanWindowAutoWindowFP)ObitDConCleanWindowAutoWindow;
  theClass->ObitDConCleanWindowReplaceField  = 
    (ObitDConCleanWindowReplaceFieldFP)ObitDConCleanWindowReplaceField;
  theClass->ObitDConCleanWindowAddField  = 
    (ObitDConCleanWindowAddFieldFP)ObitDConCleanWindowAddField;
  theClass->ObitDConCleanWindowStats  = 
    (ObitDConCleanWindowStatsFP)ObitDConCleanWindowStats;

} /* end ObitDConCleanWindowClassDefFn */

/*---------------Private functions--------------------------*/

/**
 * Creates empty member objects, initialize reference count.
 * Parent classes portions are (recursively) initialized first
 * \param inn Pointer to the object to initialize.
 */
void ObitDConCleanWindowInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitDConCleanWindow *in = inn;

  /* error checks */
  g_assert (in != NULL);

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  /* define arrays */
  in->naxis     = NULL;
  in->maxId     = NULL;
  in->Lists     = NULL;
  in->outWindow = NULL;
  in->autoWindow= FALSE;

} /* end ObitDConCleanWindowInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 *           Actually it should be an ObitDConCleanWindow* cast to an Obit*.
 */
void ObitDConCleanWindowClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitDConCleanWindow *in = inn;
  olong i;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* delete this class members */
  for (i=0; i<in->nfield; i++) {
    if ((in->naxis[i]) && (ObitMemValid (in->naxis[i])))
      in->naxis[i] = ObitMemFree (in->naxis[i]);

    /* free outer windows */
    freeWindowListElem ((WindowListElem*)in->outWindow[i]);
  
    /* Free Window list elements */
    ObitDConCleanWindowDel (in, i+1, -1, NULL);
    
    g_slist_free (in->Lists[i]);  /* Head of the list */
  } /* end loop over fields */

  if ((in->Lists) && (ObitMemValid (in->Lists)))
    in->Lists = ObitMemFree (in->Lists);
  if ((in->outWindow) && (ObitMemValid (in->outWindow)))
    in->outWindow = ObitMemFree (in->outWindow);
  if ((in->naxis) && (ObitMemValid (in->naxis))) 
    in->naxis = ObitMemFree (in->naxis);
  if ((in->maxId) && (ObitMemValid (in->maxId)))
    in->maxId = ObitMemFree (in->maxId);
      
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);
  
} /* end ObitDConCleanWindowClear */

/**
 * Determine radius of CLEAN window
 * Size is determined from a histogram of the average ratio to the 
 * value at PeakPos as a function of distance from  PeakPos.
 * The size is determined by the 10% point or 3 sigma whichever is less.
 * \param image   Image pixel array
 * \param PeakPos position (x, y, 0-rel) of position in image
 * \param sigma   RMS of "noise" level in image
 * \return the radius to use for this window.
 */
olong  GetWindowSize(ObitFArray *image, olong *PeakPos, ofloat sigma)
{
#define NWINSIZHIST  20  /* Number of values in histogram */
  olong size = 3;
  ofloat *PeakPtr, *Offset, Peak, hist[NWINSIZHIST+1], minHist;
  ofloat fblank =  ObitMagicF();
  olong count[NWINSIZHIST+1], i, x, y, x1, x2, y1, y2, icell;

  /* Peak value */
  PeakPtr = ObitFArrayIndex(image, PeakPos);
  if (PeakPtr==NULL) return 0;
  Peak = *PeakPtr;
  if ((Peak==0.0) || (Peak==fblank)) return size;

  /* Get histogram */
  for (i=0; i<=NWINSIZHIST; i++) {count[i] = 0; hist[i] = 0.0;}
  x1 = MAX (0, PeakPos[0]-NWINSIZHIST);
  x2 = MIN (image->naxis[0]-1, PeakPos[0]+NWINSIZHIST);
  y1 = MAX (0, PeakPos[1]-NWINSIZHIST);
  y2 = MIN (image->naxis[1]-1, PeakPos[1]+NWINSIZHIST);
  for (y=y1; y<=y2; y++) {
    for (x=x1; x<=x2; x++) {
      Offset = PeakPtr + (x-PeakPos[0]) + (y-PeakPos[1])*image->naxis[0];
      if (*Offset!=fblank) {
	icell = (olong) (0.5 + sqrt ((ofloat)(x-PeakPos[0])*(ofloat)(x-PeakPos[0]) + 
				     (ofloat)(y-PeakPos[1])*(ofloat)(y-PeakPos[1])));
	icell = MIN (icell, NWINSIZHIST);
	count[icell]++;
	hist[icell] += (*Offset) / Peak;
      }
    } /* end loop in x */
  } /* end loop in y */

  /* Normalize */
  for (i=0; i<=NWINSIZHIST; i++) if (count[i]>0) hist[i] /= count[i];
  
  /* Minimum value */
  minHist = MAX (0.1, 3.0*sigma/Peak);

  /* Determine wise of window to include source */
  for (i=0; i<NWINSIZHIST; i++) if (hist[i]<minHist) break;
  size = MAX (1, i-1);

  return size;
} /* end GetWindowSize*/

/**
 * WindowListElem Constructor.
 * \param Id       Id to use for window
 * \param type   Window type
 * \param window parameters, depends on type
 * \return the new object.
 */
WindowListElem*  newWindowListElem (olong Id, ObitDConCleanWindowType type, 
				    olong *window)
{
  WindowListElem* me;

  /* allocate structure */
  me = ObitMemAlloc0Name(sizeof(WindowListElem),"WindowListElem");

  /* initialize values */
  me->Id = Id;
  me->type = type;
  switch (me->type) {
  case OBIT_DConCleanWindow_rectangle:
  case OBIT_DConCleanWindow_unrectangle:
    /* Copy 4 */
    me->window[0] = window[0];
    me->window[1] = window[1];
    me->window[2] = window[2];
    me->window[3] = window[3];
    break;
  case OBIT_DConCleanWindow_round:
  case OBIT_DConCleanWindow_unround:
    /* Copy 3 */
    me->window[0] = window[0];
    me->window[1] = window[1];
    me->window[2] = window[2];
    break;
  default:
    g_error ("Undefined Clean window type");
    return me;
   }; /* end switch by window type */

  return me;
} /* end newWindowListElem */

/**
 * WindowListElem destructor.
 * \param in Object to delete
 * \return NULL.
 */
WindowListElem* freeWindowListElem (WindowListElem *in)
{
  if (in == NULL) return NULL;

  /* deallocate structure */
  ObitMemFree (in);
  
  return NULL;
} /* end freeWindowListElem */

/**
 * Find the WindowListElem in a GSList with a given Id.
 * \param glist  GSList of WindowListElems
 * \param Id     Window Id
 * \return element pointer if found, else NULL
 */
WindowListElem* ObitDConCleanWindowFind (GSList *glist, olong Id)
{
  GSList *tlist;
  WindowListElem *elem = NULL;

  /* Loop over list */
  tlist = glist;
  while (tlist) {
    elem = (WindowListElem*)tlist->data;
    if (elem->Id==Id) return elem;  /* found */
    tlist = tlist->next;
  }
  return NULL;  /* not found */
} /*  end ObitDConCleanWindowFind */
