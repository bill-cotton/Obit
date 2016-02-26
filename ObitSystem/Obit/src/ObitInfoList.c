/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2002-2016                                          */
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
/*;Correspondence about this software should be addressed as follows: */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
#include <string.h>
#include <stdio.h>
#include "ObitInfoList.h"
#include "ObitMem.h"

/**
 * \file ObitInfoList.c
 * ObitInfoList Linked list of labeled items class function definition.
 * This facility allows storing arrays of values of the same (native) data type
 * and retrieving them by name.
 * Implementation uses the glib GSList class.
 */
/*---------------Private function prototypes----------------*/
/** Private: Find item in a list */
static ObitInfoElem* ObitInfoListFind(ObitInfoList *in, gchar *name);

/*---------------Public functions---------------------------*/

/** name of the class defined in this file */
static gchar *myClassName = "ObitInfoList";

/**
 *Constructor.
 * \return the new object.
 */
ObitInfoList* newObitInfoList (void)
{
  ObitInfoList *me; 

  /* allocate */
  me =  ObitMemAlloc0Name(sizeof(ObitInfoList),"ObitInfoList");

  /* initialize */
  strncpy (me->className, myClassName, 15);
  me->list = NULL;
  me->ReferenceCount = 1;
  me->number = 0;

  return me;
} /* end newObitInfoList */

/**
 * Unconditionally deletes object.
 * \param in Object to delete
 * \return NULL pointer.
 */
ObitInfoList* freeObitInfoList (ObitInfoList *in)
{
  GSList *tmp;

  /* error checks */
  g_assert (ObitInfoListIsA(in));

  /* loop through list deleting elements */
  tmp = in->list;
  while (tmp!=NULL) {
    if (tmp->data) freeObitInfoElem(tmp->data);
    tmp = g_slist_next(tmp);
  }

  /* delete members  */
  g_slist_free(in->list);

  /* delete object */
  ObitMemFree (in);

  return NULL;
} /* end freeObitInfoList */

/**
 * Copy constructor.
 * \param in  The object to copy
 * \return pointer to the new object.
 */
ObitInfoList* ObitInfoListCopy (ObitInfoList* in)
{
  GSList  *tmp;
  ObitInfoElem *elem, *telem;
  ObitInfoList *out=NULL;

  /* error checks */
  g_assert (ObitInfoListIsA(in));

  /* create output */
  out = newObitInfoList();

  /* number of elements */
  out->number = in->number;
  /* loop through list copying elements */
  tmp = in->list;
  while (tmp!=NULL) {
    elem = (ObitInfoElem*)tmp->data;
    /* make copy */
    telem = ObitInfoElemCopy(elem);
    out->list = g_slist_prepend(out->list, telem); /* add to new list */
    tmp = g_slist_next(tmp);
  }

  /* reverse to get into same order */
  out->list = g_slist_reverse(out->list);

  return out;
} /* end ObitInfoListCopy */

/**
 * To reference a pointer, incrementing ReferenceCount and returning 
 * the pointer.
 * This function should always be used to copy pointers as this 
 * will ensure a proper reference count.
 * \param in Pointer to object to link.
 * \return the pointer to me.
 */
ObitInfoList* ObitInfoListRef (ObitInfoList* in)
{
  /* error checks */
  g_assert (ObitInfoListIsA(in));
 
 /* increment reference count */
  in->ReferenceCount++;
  return in;

} /* end ObitInfoListRef */

/**
 * Always use this function to dismiss an object as it will
 * ensure that the object is only deleted when there are no more 
 * pointers to it.
 * \param  in Pointer to object to unreference.
 * \return NULL pointer.
 */
ObitInfoList* ObitInfoListUnref (ObitInfoList* in)
{
  /* error checks */
  if (in==NULL) return NULL;
  g_assert (ObitInfoListIsA(in));

   /* decrement reference count, delete if non positive */
  in->ReferenceCount--;
  if (in->ReferenceCount<=0) freeObitInfoList(in);

 return NULL;
} /* end ObitInfoListUnref */

/**
 * Copy constructor.
 * \param in  The input object to copy
 * \param out The output object; if NULL, new one created
 * \return pointer to the new object.
 */
ObitInfoList* ObitInfoListCopyData (ObitInfoList* in, ObitInfoList* out)
{
  GSList  *tmp;
  ObitInfoElem *elem, *telem, *xelem;

  /* error checks */
  g_assert (ObitInfoListIsA(in));
  if (out!=NULL) g_assert (ObitInfoListIsA(out));

  /* create output if needed */
  if (out==NULL) out = newObitInfoList();

  
  /* loop through list copying elements */
  tmp = in->list;
  while (tmp!=NULL) {
    elem = (ObitInfoElem*)tmp->data;
    /* Check for bad entry */
    if (elem->iname==NULL) goto loop;
    out->number++; /* number of elements */
    /* Check if it's already there and is so just update */
    xelem = ObitInfoListFind(out, elem->iname);
    /* Check for bad entry */
    if (xelem->iname==NULL) goto loop;
    if (xelem==NULL) { /* not there */
      /* make copy to attach to list */
      telem = ObitInfoElemCopy(elem);
      out->list = g_slist_prepend(out->list, telem); /* add to new list */
    } else { /* already there - replace */
      /* Delete Old */
      ObitInfoListRemove (out, xelem->iname);

      /* make copy to attach to list */
      telem = ObitInfoElemCopy(elem);
      out->list = g_slist_prepend(out->list, telem); /* add to new list */
    }
  loop:
    tmp = g_slist_next(tmp);
  }

  /* reverse to get into same order */
  out->list = g_slist_reverse(out->list);

  return out;
} /* end ObitInfoListCopyData */

/**
 * Copy items from one InfoList to another if they are given in a list.
 * \param in   The input object to copy
 * \param out  The output object
 * \param list NULL terminated list of element names
 */
void  ObitInfoListCopyList (ObitInfoList* in, ObitInfoList* out, gchar **list)
{
  GSList  *tmp;
  olong i;
  gboolean wanted;
  ObitInfoElem *elem, *telem, *xelem;

  /* error checks */
  g_assert (ObitInfoListIsA(in));
  g_assert (ObitInfoListIsA(out));

  /* loop through list copying elements */
  tmp = in->list;
  while (tmp!=NULL) {
    elem = (ObitInfoElem*)tmp->data;
    
    /* Check for bad entry */
    if (elem->iname==NULL) goto loop;

    /* Is this one on the list */
    wanted = FALSE;
    i = 0;
    while (list[i]) {
      if (!strcmp(elem->iname, list[i])) { wanted = TRUE; break;};
      i++;
    }
    
    if (wanted) { /* copy */
      out->number++; /* number of elements */
      /* Check if it's already there and is so just update */
      xelem = ObitInfoListFind(out, elem->iname);
      if (xelem==NULL) { /* not there */
	/* make copy to attach to list */
	telem = ObitInfoElemCopy(elem);
	out->list = g_slist_prepend(out->list, telem); /* add to output list */
      } else { /* already there - replace*/
	/* Delete Old */
	ObitInfoListRemove (out, xelem->iname);
	
	/* make copy to attach to list */
	telem = ObitInfoElemCopy(elem);
	out->list = g_slist_prepend(out->list, telem); /* add to output list */
      }
    }
  loop:
    tmp = g_slist_next(tmp);
  }
} /* end ObitInfoListCopyList */

/**
 * opy entries from one list to another controlled by a list with rename.
 * \param in      The input object to copy
 * \param out     The output object
 * \param inList  NULL terminated list of element names
 * \param outList NULL terminated list of element names, entries must correspond to inList
 */
void ObitInfoListCopyListRename(ObitInfoList* in, ObitInfoList* out, 
				gchar **inList, gchar **outList)
{
  GSList  *tmp;
  olong i;
  gboolean wanted;
  ObitInfoElem *elem, *telem, *xelem;

  /* error checks */
  g_assert (ObitInfoListIsA(in));
  g_assert (ObitInfoListIsA(out));

  /* loop through list copying elements */
  tmp = in->list;
  while (tmp!=NULL) {
    elem = (ObitInfoElem*)tmp->data;
    
    /* Check for bad entry */
    if (elem->iname==NULL) goto loop;

    /* Is this one on the list */
    wanted = FALSE;
    i = 0;
    while (inList[i] && outList[i]) {
      if (!strcmp(elem->iname, inList[i])) { wanted = TRUE; break;};
      i++;
    }
    
    if (wanted) { /* copy */
      out->number++; /* number of elements */
      /* Check if it's already there and is so just update */
      xelem = ObitInfoListFind(out, outList[i]);
      if (xelem==NULL) { /* not there */
	/* make copy with rename to attach to list */
	telem = ObitInfoElemCopy(elem);
	/* Change name 
	g_free(telem->iname);*/
	telem->iname = g_strdup(outList[i]);
	out->list = g_slist_prepend(out->list, telem); /* add to output list */
      } else { /* already there - replace*/
	/* Delete Old */
	ObitInfoListRemove (out, xelem->iname);
	
	/* make copy with rename to attach to list */
	telem = ObitInfoElemCopy(elem);
	/* Change name
	g_free(telem->iname); */
	telem->iname = g_strdup(outList[i]);
	out->list = g_slist_prepend(out->list, telem); /* add to output list */
      }
    }
  loop:
    tmp = g_slist_next(tmp);
  }
} /* end ObitInfoListCopyListRename */


/**
 * Copy entries from one list to another adding a prefix
 * \param in     The input object to copy
 * \param out    The output object
 * \param prefix Prefix to add
 */
void ObitInfoListCopyAddPrefix(ObitInfoList* in, ObitInfoList* out, 
			       gchar *prefix)
{
  GSList  *tmp;
  ObitInfoElem *elem, *telem, *xelem;
  gchar *newName = NULL;

  /* error checks */
  g_assert (ObitInfoListIsA(in));
  g_assert (ObitInfoListIsA(out));

  /* loop through list copying elements */
  tmp = in->list;
  while (tmp!=NULL) {
    elem = (ObitInfoElem*)tmp->data;
    /* Check for bad entry */
    if (elem->iname==NULL) goto loop;
    /* New name */
    newName = g_strconcat(prefix, elem->iname, NULL);
    
    out->number++; /* number of elements */
    /* Check if it's already there and is so just update */
    xelem = ObitInfoListFind(out, newName);
    if (xelem==NULL) { /* not there */
      /* make copy to attach to list */
      telem = newObitInfoElem(newName, elem->itype, elem->idim, elem->data);
      out->list = g_slist_prepend(out->list, telem); /* add to output list */
    } else { /* already there - replace*/
      /* Delete Old */
      ObitInfoListRemove (out, xelem->iname);
      
      /* make copy to attach to list */
      telem = newObitInfoElem(newName, elem->itype, elem->idim, elem->data);
      out->list = g_slist_prepend(out->list, telem); /* add to output list */
    }
    g_free(newName);
  loop:
    tmp = g_slist_next(tmp);
  }
} /* end ObitInfoListCopyAddPrefix */


/**
 * Copy entries with a given prefix from one list to another
 * \param in     The input object to copy
 * \param out    The output object
 * \param prefix Prefix to look for
 * \param strip  If TRUE then remove prefix in output list
 */
void ObitInfoListCopyWithPrefix(ObitInfoList* in, ObitInfoList* out, 
				gchar *prefix, gboolean strip)
{
  GSList  *tmp;
  olong nchk;
  ObitInfoElem *elem, *telem, *xelem;
  gboolean wanted;
  gchar *newName = NULL;

  /* error checks */
  g_assert (ObitInfoListIsA(in));
  g_assert (ObitInfoListIsA(out));

  /* How many characters to compare? */
  nchk = strlen(prefix);

  /* loop through list copying elements */
  tmp = in->list;
  while (tmp!=NULL) {
    elem = (ObitInfoElem*)tmp->data;
    /* Check for bad entry */
    if (elem->iname==NULL) goto loop;

    /* Does this begin with the prefix? */
    wanted = !strncmp(prefix, elem->iname, nchk);     

    if(wanted) {
      /* New name */
      if (strip) newName = g_strdup(&elem->iname[nchk]);
      else newName = g_strdup(elem->iname);

      out->number++; /* number of output elements */
      /* Check if it's already there and is so just update */
      xelem = ObitInfoListFind(out, newName);
      if (xelem==NULL) { /* not there */
	/* make copy to attach to list */
	telem = newObitInfoElem(newName, elem->itype, elem->idim, elem->data);
	out->list = g_slist_prepend(out->list, telem); /* add to output list */
      } else { /* already there - replace*/
	/* Delete Old */
	ObitInfoListRemove (out, xelem->iname);
	
	/* make copy to attach to list */
	telem = newObitInfoElem(newName, elem->itype, elem->idim, elem->data);
	out->list = g_slist_prepend(out->list, telem); /* add to output list */
      }
      g_free(newName);
    } /* end if wanted */
  loop:
    tmp = g_slist_next(tmp);
  }
} /* end ObitInfoListCopyWithPrefix */


/**
 * Store information in the infoList.
 * If the requested keyword (name) does not exist then one is created
 * with the type, and dimension specified.
 * If a previous entry exists, then the type and dimensionality provided must 
 * match.
 * \param in   Pointer to InfoList.
 * \param name The label (keyword) of the information.
 * \param type Data type of data element (enum defined in ObitInfoList class.
 * \param dim  Dimensionality of datum.
 *             Note: for strings, the first element is the length in char.
 * \param data Pointer to the data.  If the data is a multidimensional string array
 *             pass this as the pointer the the string array.
 * \param err  Obit Error stack for information.
 */
void ObitInfoListPut(ObitInfoList *in, 
		      gchar* name, ObitInfoType type, gint32 *dim, 
		      gconstpointer data, ObitErr *err) 
{
  ObitInfoElem *elem;
  gchar *routine = "ObitInfoListPut";

  /* error checks */
  g_assert (ObitInfoListIsA(in));
  g_assert (ObitErrIsA(err));
  g_assert (name != NULL);
  g_assert (dim != NULL);
  if (dim[0] <= 0) {
      Obit_log_error(err, OBIT_Error, 
		     "%s:No dimensionality given for %s",routine, name);
      return;
  }

  /* look up */
  elem = ObitInfoListFind(in, name);

  if (elem!=NULL) { /* Found it */
    /* Update */
    if (!ObitInfoElemUpdate(elem, type, dim, data, TRUE)) {
      /* incompatable */
      Obit_log_error(err, OBIT_Error, 
		     "%s: Wrong type or dimension for %s", routine, name);
      return; /* bail out */
    } /* end of mismatch */
    
    return;
  }

  /* add new one */
  elem =  newObitInfoElem(name, type, dim, data);
  in->list = g_slist_append (in->list, elem); /* add to list */
  in->number++; /* keep count */

} /* end ObitInfoListPut */

/**
 * Store information in the infoList, if the type and size do not
 * match then redefine the entry.
 * If the requested keyword (name) does not exist then one is created
 * with the type, and dimension specified.
 * \param in   Pointer to InfoList.
 * \param name The label (keyword) of the information.
 * \param type Data type of data element (enum defined in ObitInfoList class.
 * \param dim  Dimensionality of datum.
 *             Note: for strings, the first element is the length in char.
 * \param data Pointer to the data.
 */
void ObitInfoListAlwaysPut(ObitInfoList *in, 
			   gchar* name, ObitInfoType type, gint32 *dim, 
			   gconstpointer data) 
{
  ObitInfoElem *elem;

  /* error checks */
  g_assert (ObitInfoListIsA(in));
  g_assert (name != NULL);
  g_assert (dim != NULL);
  g_assert (dim[0] > 0);

  /* look up */
  elem = ObitInfoListFind(in, name);

  if (elem!=NULL) { /* Found it */
    /* Update */
    if (!ObitInfoElemUpdate(elem, type, dim, data, FALSE)) {
      /* incompatable - delete old */
      ObitInfoListRemove (in, name);
      elem = NULL;
    } /* end of mismatch */
  }

  /* Have it? */
  if (elem!=NULL) return;

  /* add new one */
  elem =  newObitInfoElem(name, type, dim, data);
  in->list = g_slist_append (in->list, elem); /* add to list */
  in->number++; /* keep count */

} /* end ObitInfoListAlwaysPut */

/**
 * Retrieve information stored in the infoList if it is available.
 * \param in   Pointer to InfoList.
 * \param name The label (keyword) of the information.
 * \param type (output) data type of data element.
 * \param dim  (output) dimensionality of datum.
 *             Note: for strings, the first element is the length in char.
 *             All Strings and string arrays will have a final NULL which
 *             is one more byte than described in dim.
 * \param data (output) pointer to the data; 
 *             note: data will be copied to this location
 *             which should have been allocated sufficiently.
 * \return TRUE if found, else FALSE.
 */
gboolean 
ObitInfoListGetTest(ObitInfoList *in, 
		    gchar* name, ObitInfoType *type, gint32 *dim, 
		    gpointer data)
{
  ObitInfoElem *elem;
  olong size;

  /* error checks */
  g_assert (ObitInfoListIsA(in));
  g_assert (name != NULL);
  g_assert (dim != NULL);

  /* look up */
  elem = ObitInfoListFind(in, name);
  if (elem==NULL) return FALSE; /* not found return */

  /* copy information */
  *type = elem->itype;
  size = MAXINFOELEMDIM * sizeof (gint32);
  g_memmove (dim, elem->idim, size);
  g_memmove (data, elem->data, elem->size);

  return TRUE;
} /* end ObitInfoListGetTest */

/**
 * Retrieve information stored in the infoList by name.
 * NOTE: String arrays are not returned in the form they
 * are written, they are all in one large, rectangular gchar array.
 * \param in   Pointer to InfoList.
 * \param name The label (keyword) of the information.
 * \param type (output) data type of data element.
 * \param dim  (output) dimensionality of datum.
 *             Note: for strings, the first element is the length in char.
 *             All Strings and string arrays will have a final NULL which
 *             is one more byte than described in dim.
 * \param data (output) pointer to the data; 
 *             note: data will be copied to this location
 *             which should have been allocated sufficiently.
 * \param err  (output) Obit Error stack for information.
 * \return TRUE if found, else FALSE.
 */
gboolean ObitInfoListGet(ObitInfoList *in, 
			 gchar* name, ObitInfoType *type, gint32 *dim, 
			 gpointer data, ObitErr *err)
{
  ObitInfoElem *elem;
  olong size;
  gchar *routine = "ObitInfoListGet";

  /* error checks */
  g_assert (ObitInfoListIsA(in));
  g_assert (ObitErrIsA(err));
  g_assert (name != NULL);
  g_assert (dim != NULL);

  /* look up */
  elem = ObitInfoListFind(in, name);
  if (elem==NULL) { /* not found */
      Obit_log_error(err, OBIT_Error, 
		     "%s: Could not find element %s", routine, name);
      return FALSE;
  }

  /* copy information */
  *type = elem->itype;
  size = MAXINFOELEMDIM * sizeof (gint32);
  g_memmove (dim, elem->idim, size);
  g_memmove (data, elem->data, elem->size);

  return TRUE;
} /* end ObitInfoListGet */

/**
 * Retrieve information about an item in the infoList by name.
 * \param in   Pointer to InfoList.
 * \param name The label (keyword) of the information.
 * \param type (output) data type of data element.
 * \param dim  (output) dimensionality of datum.
 *             Note: for strings, the first element is the length in char.
 *             All Strings and string arrays will have a final NULL which
 *             is one more byte than described in dim.
 * \param err  (output) Obit Error stack for information.
 * \return TRUE if found, else FALSE.
 */
gboolean ObitInfoListInfo(ObitInfoList *in, 
			 gchar* name, ObitInfoType *type, gint32 *dim, 
			 ObitErr *err)
{
  ObitInfoElem *elem;
  olong size;

  /* error checks */
  g_assert (ObitInfoListIsA(in));
  g_assert (ObitErrIsA(err));
  g_assert (name != NULL);
  g_assert (dim != NULL);

  /* look up */
  elem = ObitInfoListFind(in, name);
  if (elem==NULL) { /* not found */
      return FALSE;
  }

  /* copy information */
  *type = elem->itype;
  size = MAXINFOELEMDIM * sizeof (gint32);
  g_memmove (dim, elem->idim, size);

  return TRUE;
} /* end ObitInfoListInfo */

/**
 * Return pointer to information stored in the infoList by name.
 * \param in   Pointer to InfoList.
 * \param name The label (keyword) of the information.
 * \param type (output) data type of data element.
 * \param dim  (output) dimensionality of datum.
 *             Note: for strings, the first element is the length in char.
 *             All Strings and string arrays will have a final NULL which
 *             is one more byte than described in dim.
 * \param data (output) pointer to the data;  NULL if not found.
 * \param err  (output) Obit Error stack for information.
 * \return TRUE if found, else FALSE.
 */
gboolean ObitInfoListGetP(ObitInfoList *in, 
			 gchar* name, ObitInfoType *type, gint32 *dim, 
			 gpointer *data)
{
  ObitInfoElem *elem;
  olong size;

  /* error checks */
  g_assert (ObitInfoListIsA(in));
  g_assert (name != NULL);
  g_assert (dim != NULL);

  *data = NULL; /* initialize output */
  /* look up */
  elem = ObitInfoListFind(in, name);
  if (elem==NULL) return FALSE;

  /* copy information */
  *type = elem->itype;
  size = MAXINFOELEMDIM * sizeof (gint32);
  g_memmove (dim, elem->idim, size);
  *data = elem->data; /* just pointer */

  return TRUE;
} /* end ObitInfoListGetP */

/**
 * Retrieve information stored in the infoList by number.
 * \param in     Pointer to InfoList.
 * \param number 1-rel item number
 * \param name   (output) Pointer to the label (keyword) of the information;
 *               the name will be copied into this address.
 * \param type   (output) data type of data element.
 * \param dim    (output) dimensionality of datum.
 *               Note: for strings, the first element is the length in char.
 *               All Strings and string arrays will have a final NULL which
 *               is one more byte than described in dim.
 * \param data   (output) pointer to the data; 
 *               note: data will be copied to this location
 *               which should have been allocated sufficiently.
 * \param err  (output) Obit Error stack for information.
 * \return TRUE if found, else FALSE.
 */
gboolean 
ObitInfoListGetNumber (ObitInfoList *in,  olong number,
		      gchar** name, ObitInfoType *type, gint32 *dim, 
		      gpointer data, ObitErr *err)
{
  ObitInfoElem *elem;
  GSList *tmp;
  olong i, size;
  gchar *routine = "ObitInfoListGetNumber";

  /* error checks */
  g_assert (ObitInfoListIsA(in));
  g_assert (ObitErrIsA(err));
  g_assert (name != NULL);
  g_assert (dim != NULL);
  if ((number<0) || (number>in->number)) { /* number valid */
    Obit_log_error(err, OBIT_Error, 
		   "%s: Invalid item number %d, max %d", routine, number, in->number);
      return FALSE;
  }

  /* look up */
  tmp = in->list;
  for (i=0; i<number-1; i++) {
    tmp = g_slist_next(tmp);
    if (tmp==NULL) break;  /* problem? */
  }
  if (tmp==NULL) { /* not found */
      Obit_log_error(err, OBIT_Error, 
        "%s: I appear to have been corrupted", routine);
      return FALSE;
  }

  elem = (ObitInfoElem*)tmp->data;
  if ((elem==NULL) || (elem->iname==NULL)) { /* not found or bad */
      Obit_log_error(err, OBIT_Error, 
        "%s: I appear to have been corrupted", routine);
      return FALSE;
  }

  /* copy information */
  *type = elem->itype;
  size = strlen (elem->iname);
  *name = elem->iname;  /* set output pointer rather than copy */
  size = MAXINFOELEMDIM * sizeof (gint32);
  g_memmove (dim, elem->idim, size);
  g_memmove (data, elem->data, elem->size);

  return TRUE;
} /* end ObitInfoListGetNumber */

/**
 * Return pointer to information stored in the infoList by number.
 * \param in     Pointer to InfoList.
 * \param number 1-rel item number
 * \param name   (output) Pointer to the label (keyword) of the information;
 *               the name will be copied into this address.
 * \param type   (output) data type of data element.
 * \param dim    (output) dimensionality of datum.
 *               Note: for strings, the first element is the length in char.
 *               All Strings and string arrays will have a final NULL which
 *               is one more byte than described in dim.
 * \param data   (output) pointer to the data; NULL if not found.
 * \param err  (output) Obit Error stack for information.
 * \return TRUE if found, else FALSE.
 */
gboolean 
ObitInfoListGetNumberP (ObitInfoList *in,  olong number,
			gchar** name, ObitInfoType *type, gint32 *dim, 
			gpointer *data)
{
  ObitInfoElem *elem;
  GSList *tmp;
  olong i, size;

  /* error checks */
  g_assert (ObitInfoListIsA(in));
  g_assert (name != NULL);
  g_assert (dim != NULL);
  *data = NULL; /* initialize output */
  if ((number<0) || (number>in->number)) return FALSE;

  /* look up */
  tmp = in->list;
  for (i=0; i<number-1; i++) {
    tmp = g_slist_next(tmp);
    if (tmp==NULL) break;  /* problem? */
  }
  if (tmp==NULL) return FALSE;

  elem = (ObitInfoElem*)tmp->data;
  if (elem==NULL) return FALSE;
  /* Check for bad entry */
  if (elem->iname==NULL) return FALSE;

  /* copy information */
  *type = elem->itype;
  size = strlen (elem->iname);
  *name = elem->iname;  /* set output pointer rather than copy */
  size = MAXINFOELEMDIM * sizeof (gint32);
  g_memmove (dim, elem->idim, size);
  *data = elem->data; /* just pointer */

  return TRUE;
} /* end ObitInfoListGetNumberP */

/**
 * Remove the item with keyword name from the list.
 * Item is destroyed; does nothing if item not found.
 * \param in   Pointer to InfoList.
 * \param name The label (keyword) of the information.
 */
void ObitInfoListRemove (ObitInfoList *in, gchar* name)
{
  ObitInfoElem *elem;

  /* error checks */
  g_assert (ObitInfoListIsA(in));
  g_assert (name != NULL);

  /* look up */
  elem = ObitInfoListFind(in, name);
  if (elem==NULL) return; /* nothing to do */

  /* remove from list */
  in->list = g_slist_remove(in->list, elem);
  in->number--; /* keep count */

  /* delete element */
  freeObitInfoElem(elem);

} /* end ObitInfoListRemove */

/**
 * Change the dimension and/or data type of an entry.
 * Create an entry in list if it doesn't previously exist,
 * any existing data will be lost.
 * \param in   Pointer to InfoList.
 * \param name The label (keyword) of the element to change.
 * \param type Mew data type of data element (enum).
 * \param dim  New dimensionality of datum.
 */
void ObitInfoListResize(ObitInfoList *in, 
			gchar* name, ObitInfoType type, gint32 *dim)
{
  ObitInfoElem *elem;

  /* error checks */
  g_assert (ObitInfoListIsA(in));
  g_assert (name != NULL);
  g_assert (dim != NULL);
  g_assert (dim[0] > 0);

  /* look up */
  elem = ObitInfoListFind(in, name);
  if (elem==NULL) { /* not there create new */
    elem = newObitInfoElem(name, type, dim, NULL);
    in->list = g_slist_append (in->list, elem); /* add to list */
  } else { /* found it, resize */
    ObitInfoElemResize (elem, type, dim); 
  }
} /* end ObitInfoListResize */

/**
 * Print the contents of an InfoList to file
 * \param file  Where to write output
 */
void ObitInfoListPrint (ObitInfoList *in, FILE *file)
{
  GSList *tmp;

  /* error checks */
  g_assert (ObitInfoListIsA(in));
  g_assert (file != NULL);

  fprintf (file, "Listing Of ObitInfoList\n");

  /* loop through list */
  tmp = in->list;
  while (tmp!=NULL) {
    if (tmp->data)
      ObitInfoElemPrint ((ObitInfoElem*)tmp->data, file);
    tmp = g_slist_next(tmp);
  }
  fprintf (file, "\n");
} /* endObitInfoList  */


/**
 * Determines if the input object is a member of this class
 * \param in Pointer to object to test.
 * \return TRUE if member else FALSE.
 */
gboolean ObitInfoListIsA (ObitInfoList* in)
{
  gboolean out;

  /* error checks */
  if (in == NULL) return FALSE;
  if (in->className == NULL) return FALSE;

  /* compare class name member */
  out = !strcmp(in->className, myClassName);

  return out;
} /* end ObitInfoListIsA */

/*---------------Private functions---------------------------*/

/** 
 * Find the element with a given name in an ObitInfoList
 * \param in   Pointer to InfoList.
 * \param name The label (keyword) of the element to change.
 * \return pointer to ObitInfoElem or NULL if not found.
 */
static ObitInfoElem* ObitInfoListFind(ObitInfoList *in, gchar *name)
{
  GSList *tmp;
  ObitInfoElem *elem;
  gchar *tstring;

  /* error checks */
  g_assert (ObitInfoListIsA(in));
  g_assert (name != NULL);

  /* Get copy of name will all leading and trailing whitespace removed */
  tstring = g_strstrip(g_strdup(name));

  /* loop through list testing elements */
  tmp = in->list;
  while (tmp!=NULL) {
    elem = (ObitInfoElem*)tmp->data;
    /* check if this is a match */
    if (ObitInfoElemTest(elem,tstring)) 
      {g_free(tstring); return elem;}
    tmp = g_slist_next(tmp);
  }
  g_free(tstring);
  return NULL; /* didn't find */
} /* end ObitInfoListFind */

