#ifndef _POSIX_SOURCE
#define _POSIX_SOURCE
#endif
#include <stdio.h>
#include <stdlib.h>

/*   Define the maximum number of objects                             */
#define MAXOBJ  10240

/* Linked list element Structure                                      */
/* The memory allocated for this structure will include a data area   */
/* after the end of the structure.                                    */
struct ll_element {
    char          keyword[8];   /* Keyword name 8 char  */
    struct ll_element *next;    /* Pointer to next element, void=> last */
    int           dataoff;      /* offset from struct to beginning of
                                   data array in bytes*/
    int           lendata;      /* Number of bytes in data array */
    int           type;         /* Data type code */
    int           ndim;         /* number of dimensions */
    int           dim[7];       /* dimensionality array */
                                /* data array follows dim[ndim-1]*/
};

/* AIPS object Structure                                              */
struct AIPS_object {
    char  name[32];             /* Object name */
    char  class[8];             /* Object class */
    int   catblk[256];          /* AIPS catalog header */
    struct ll_element *first;          /* First element of linked list */
    struct ll_element *last;           /* Last element of linked list */
    struct ll_element *select;         /* Selected element of linked list */
};


/* Object directory:                                                  */
struct object_dir {
 struct AIPS_object *objpnt[MAXOBJ]; /* list of pointers to AIPS-Objects */
 int       numobj;               /* max number in objpnt assigned */
 char last_obj[32][3];           /* names of last several objects searched */
 int  last_num[3];               /* numbers of last few objects */
} task_dir;                      /* define task_dir for this task */

/* Functional names of the Fortran callable routines:
Function:               Routine name
init_dir                   ZOINTD
create_object              ZOCROB
destroy_object             ZODEOB
find_object                ZOFNOB
copy_object                ZOCPOB
find_llelem                ZOFNLE
info_llelem                ZOINLE
store_data                 ZOSTDT
fetch_data                 ZOFEDT
store_catblk               ZOSTCT
fetch_catblk               ZOFECT
fetch_name                 ZOFENM
*/
/*   Function Prototypes */
#if __STDC__
void zointd_ ();
void zocrob_ (char *name, char *class, int *objnum, int *ierr);
void zodeob_ (int *objnum);
void zofnob_ (char *name, int *objnum);
void zocpob_ (char *namein, char *nameout, int *ierr);
void zofnle_ (int *objnum, char *keyword, int *ierr);
void zoinle_ (int *objnum, char *keyword, int *type, int *ndim,
              int *dim, int *ierr);
void zostdt_ (int *objnum, char *keyword, int *type, int *ndim,
              int *dim, char *data, int *ierr);
void zofedt_ (int *objnum, char *keyword, int *type, int *ndim,
              int *dim, char *data, int *ierr);
void zostct_ (int *objnum, int *catblk, int *ierr);
void zofect_ (int *objnum, int *catblk, int *ierr);
void zofenm_ (int *objnum, char *name, char *class, int *ierr);
struct ll_element* create_llelem_ (char *keyword, int *type,
                                  int *ndim, int *dim, int *ierr);
#else
struct ll_element* create_llelem_ ();
#endif

/*
Functions called from Fortran:
*/

#if __STDC__
   void zointd_ ()
#else
   void zointd_ ()
#endif
/*--------------------------------------------------------------------*/
/*! Initialize AIPS Object directory                                  */
/*# OOPS                                                              */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 1995, 1997, 1999-2001                              */
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
/*;  Correspondence concerning AIPS should be addressed as follows:   */
/*;         Internet email: aipsmail@nrao.edu.                        */
/*;         Postal address: AIPS Project Office                       */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
/*---------------------------------------------------------------------*/
/* Initialize AIPS object directory                                    */
/*   Inputs:                                                           */
/*---------------------------------------------------------------------*/
{
    extern struct object_dir task_dir;   /* task object directory */
    int  i;
/*---------------------------------------------------------------------*/
    task_dir.numobj = 0;
    task_dir.last_num[0] = 0;
    task_dir.last_num[1] = 0;
    task_dir.last_num[2] = 0;
    for (i=0;i<MAXOBJ;i++) task_dir.objpnt[i]=NULL;
    return;
}

#if __STDC__
   void zocrob_ (char *name, char *class, int *objnum,
                        int *ierr)
#else
   void zocrob_ (name, class, objnum, ierr)
   char name[], class[];
   int *objnum, *ierr;
#endif
/*---------------------------------------------------------------------*/
/*! Create Storage for AIPS Object                                     */
/*# OOPS                                                               */
/*---------------------------------------------------------------------*/
/* Create/initialize object structure                                  */
/*   Inputs:                                                           */
/*      name    H*32   Object name                                     */
/*      class   H*8    Object class                                    */
/*   Outputs:                                                          */
/*      objnum  I      Object number                                   */
/*      ierr    I      Error code, 0=>OK                               */
/*                     1=> table full                                  */
/*                     2=> allocation failed                           */
/*---------------------------------------------------------------------*/
{
    extern struct object_dir task_dir;   /* task object directory */
    struct AIPS_object *ao, aot;
    int    i;
/*---------------------------------------------------------------------*/
    *objnum = -1;
    *ierr = 0;
                                       /* find first available object  */
                                       /* number                       */
    for (i=0;task_dir.objpnt[i];i++) ;
                                       /* Check for blown table        */
    if (i>=MAXOBJ) {*ierr=1; return;}
    if ((ao = (struct AIPS_object *) malloc (sizeof(aot)))==NULL)
        {*ierr = 2; return;}
    *objnum = i+1;
    task_dir.objpnt[i] = ao;
    if ((i+1)>task_dir.numobj) task_dir.numobj = i+1;
                                       /* Initialize                  */
    for (i=0; i<32; i++) ao->name[i] = name[i];
    for (i=0; i<8; i++) ao->class[i] = class[i];
    for (i=0; i<256; i++) ao->catblk[i] = 0;
    ao->first = NULL;
    ao->last = NULL;
    ao->select = NULL;
   return;
}

#if __STDC__
   void zodeob_ (int *objnum)
#else
   void zodeob_ (objnum)
   int *objnum;
#endif
/*---------------------------------------------------------------------*/
/*! Delete Storage for AIPS Object                                     */
/*# OOPS                                                               */
/*---------------------------------------------------------------------*/
/* Destroy object structure                                            */
/*   Inputs:                                                           */
/*      objnum   I    Object number (1 rel)                            */
/*---------------------------------------------------------------------*/
{
    extern struct object_dir task_dir;   /* task object directory */
    struct AIPS_object *ao;
    struct ll_element *this, *next;
    int    i, j, k;
/*---------------------------------------------------------------------*/
    ao = task_dir.objpnt[*objnum-1];
    if (ao==NULL) return;
                                       /* Remove from object search    */
                                       /* cache                        */
    for (i=0;i<3;i++)
        {for (j=0; j<32; j++) {
            if (!(task_dir.last_obj[j][i]==ao->name[j])) break;
        }
         if (j>=32) {for (k=0; k<32; k++) task_dir.last_obj[k][i] = 0;
                     break;}
     }
                                       /* Destroy linked list         */
    this = ao->first;
    while (!(this==NULL)) {
        next = this->next;
        free(this);
        this = next;}
                                       /* Destroy AIPS object        */
    free(ao);
    task_dir.objpnt[*objnum-1] = NULL;
    return;
}

#if __STDC__
   void zofnob_ (char *name, int *objnum)
#else
   void zofnob_ (name, objnum)
   char name[];
   int *objnum;
#endif
/*---------------------------------------------------------------------*/
/*! Look up AIPS Object                                                */
/*# OOPS                                                               */
/*---------------------------------------------------------------------*/
/* Look up object                                                      */
/* A cache of the last three object located is kept.                   */
/*   Inputs:                                                           */
/*      name    H*32   Object name                                     */
/*   Outputs:                                                          */
/*      objnum  I      Object number (1 relative)                      */
/*                     0 => not found.                                 */
/*---------------------------------------------------------------------*/
{
    extern struct object_dir task_dir;   /* task object directory */
    int  i,j;
    struct AIPS_object *ao;
/*---------------------------------------------------------------------*/
    *objnum = 0;
                                       /* Check if accessed recently   */
    for (i=0;i<3;i++)
        {for (j=0; j<32; j++) {
            if (!(task_dir.last_obj[j][i]==name[j])) break;
        }
         if (j>=32) {*objnum = task_dir.last_num[i]+1; return;}
     }
    for (i=0; i<task_dir.numobj; i++) {
        ao = task_dir.objpnt[i];
        if (!(ao==NULL)) {
            for (j=0; j<32; j++) {
                if (!(ao->name[j]==name[j])) break;
            }
            if (j>=32) {*objnum = i+1; break;}
        }
    }
    if (*objnum<=0) return;    /* return if not found */
                                       /* Save in recent access list  */
    for (i=0; i<32; i++)
        {task_dir.last_obj[i][2] = task_dir.last_obj[i][2];
         task_dir.last_obj[i][1] = task_dir.last_obj[i][0];
         task_dir.last_obj[i][0] = name[i];}
    task_dir.last_num[2] = task_dir.last_num[1];
    task_dir.last_num[1] = task_dir.last_num[0];
    task_dir.last_num[0] = *objnum - 1;
    return;
}

#if __STDC__
   void zocpob_ (char *namein, char *nameout, int *ierr)
#else
   void zocpob_ (namein, nameout, ierr)
   char namein[], nameout[];
   int *ierr;
#endif
/*---------------------------------------------------------------------*/
/*! Copy AIPS Object                                                   */
/*# OOPS                                                               */
/*---------------------------------------------------------------------*/
/* Copy object                                                         */
/*   Inputs:                                                           */
/*      namein  H*32   Input object name                               */
/*      nameout H*32   Output object name                              */
/*   Outputs:                                                          */
/*      ierr    I      Error code, 0=>OK                               */
/*                     1=> Problem with input object                   */
/*                     2=> Problem creating new object                 */
/*                     3=> Problem creating new linked list element    */
/*---------------------------------------------------------------------*/
{
    extern struct object_dir task_dir;   /* task object directory */
    struct AIPS_object *aoi, *aoo;
    struct ll_element *llin, *llout, *prev;
    int    i, objin, objout, jerr, ncopy;
    char   *pin, *pout;
/*---------------------------------------------------------------------*/
    *ierr = 0;
                                       /* Find indput object           */
    zofnob_ (namein, &objin);
    if (objin<=0) {*ierr = 1; return;}
                                       /* Create new object           */
    aoi = task_dir.objpnt[objin-1];
    if (aoi==NULL) {*ierr = 1; return;}
    zocrob_ (nameout, aoi->class, &objout, &jerr);
    if ((objout<=0) | (jerr!=0)) {*ierr = 2; return;}
    aoo = task_dir.objpnt[objout-1];
                                       /* Copy catblk                 */
    for (i=0; i<256; i++) aoo->catblk[i] = aoi->catblk[i];
                                       /* Copy linked list            */
    llin = aoi->first;
    prev = NULL;
    while (!(llin==NULL))
        {llout = create_llelem_ (llin->keyword, &llin->type, &llin->ndim,
                                 llin->dim, &jerr);
         if ((!(jerr==0)) | (llout==NULL)) {*ierr = 3; return;}
         ncopy = llin->lendata;
         pin = (char *) llin + llin->dataoff;
         pout = (char *) llout + llout->dataoff;
         for (i=0; i<ncopy; i++) *(pout++) = *(pin++); /* copy data */
         if (aoo->first==NULL) aoo->first = llout;  /* first element */
         if (!(prev==NULL)) prev->next = llout; /* link to list */
         prev = llout;
         llin = llin->next;
        }
    aoo->last = llout;
    *ierr = 0;
    return;
}

#if __STDC__
   void zofnle_ (int *objnum, char *keyword, int *ierr)
#else
   void zofnle_ (objnum, keyword, ierr)
   int *objnum, *ierr;
   char keyword[];
#endif
/*---------------------------------------------------------------------*/
/*! Look up AIPS Object keyword                                        */
/*# OOPS                                                               */
/*---------------------------------------------------------------------*/
/* Look up keyword in object linked list; sets selected linked list    */
/* element pointer.                                                    */
/*   Inputs:                                                           */
/*      objnum   I     Object number, 1 rel                            */
/*      keyword  H*8   keyword name                                    */
/*   Outputs:                                                          */
/*      ierr     I     Return code, 0=>found.                          */
/*                     1 = Not found                                   */
/*                     4 = bad objnum                                  */
/*                     5 = Object doesn't exist                        */
/*---------------------------------------------------------------------*/
{
    extern struct object_dir task_dir;   /* task object directory */
    int  i,j;
    struct AIPS_object *ao;
    struct ll_element *this, *next, *select;
/*---------------------------------------------------------------------*/
    if ((*objnum<=0) | (*objnum>task_dir.numobj)) {*ierr = 4; return;}
    *ierr = 1;
    ao = task_dir.objpnt[*objnum-1];
    if (ao==NULL) {*ierr = 5; return;}
                                       /* Check the selected element   */
                                       /* first                        */
    select = ao->select;
    if (!(select==NULL)) {
        for (i=0; i<8; i++) {
            if (!(select->keyword[i]==keyword[i])) break;
        }
        if (i>=8) {*ierr = 0; return;}  /* OK leave it as is */
    }
                                       /* Check list                   */
    this = ao->first;
    while (!(this==NULL)) {
        next = this->next;
        for (i=0; i<8; i++) {
            if (!(this->keyword[i]==keyword[i])) break;
        }
        if (i>=8) {*ierr = 0; ao->select = this; return;}
        this = next;
    }
    return;
}

#if __STDC__
   void zoinle_ (int *objnum, char *keyword, int *type,
                      int *ndim, int *dim, int *ierr)
#else
   void zoinle_ (objnum, keyword, type, ndim, dim, ierr)
   int *objnum, *type, *ndim, dim[], *ierr;
   char keyword[];
#endif
/*---------------------------------------------------------------------*/
/*! Lookup information about a linked list element                     */
/*# OOPS                                                               */
/*---------------------------------------------------------------------*/
/* Lookup information about a linked list element                      */
/*   Inputs:                                                           */
/*      objnum   I     Object number                                   */
/*      keyword  H*8   keyword name                                    */
/*   Outputs:                                                          */
/*      type     I     keyword type                                    */
/*      ndim     I     number of dimensions                            */
/*      dim      I     Dimensionality array                            */
/*      ierr     I     Return code, 0=>OK.                             */
/*                     1 = Not found                                   */
/*                     4 = bad objnum                                  */
/*                     5 = Object doesn't exist                        */
/*---------------------------------------------------------------------*/
{
    extern struct object_dir task_dir;   /* task object directory */
    int  i, jerr;
    struct AIPS_object *ao;
    struct ll_element *select;
/*---------------------------------------------------------------------*/
    *ierr = 0;
    if ((*objnum<=0) | (*objnum>task_dir.numobj)) {*ierr = 4; return;}
    ao = task_dir.objpnt[*objnum-1];
    if (ao==NULL) {*ierr = 5; return;}
                                       /* Find linked list element    */
    zofnle_ (objnum, keyword, &jerr);
    if (jerr==0) {        /* found */
        select = ao->select;
        *type = select->type;
        *ndim = select->ndim;
        for (i=0; i<*ndim; i++) dim[i] = select->dim[i];
        *ierr = 0;
    }
    else {   /* Not found - error */
        *ierr = 1;
    }
    return;
}

#if __STDC__
   void zostdt_ (int *objnum, char *keyword, int *type, int *ndim,
                      int *dim, char *data, int *ierr)
#else
   void zostdt_ (objnum, keyword, type, ndim, dim, data, ierr)
   int *objnum, *type, *ndim, dim[], *ierr;
   char keyword[], data[];
#endif
/*---------------------------------------------------------------------*/
/*! Store data in AIPS Object linked list                              */
/*# OOPS                                                               */
/*---------------------------------------------------------------------*/
/* Store data in AIPS Object linked list                               */
/* Checks compatibility                                                */
/*   Inputs:                                                           */
/*      objnum   I     Object number (1-rel)                           */
/*      keyword  H*8   keyword name                                    */
/*      type     I     keyword type                                    */
/*      ndim     I     number of dimensions                            */
/*      dim      I(*)  Dimensionality array                            */
/*      data     char  data array                                      */
/*   Outputs:                                                          */
/*      ierr     I     Return code, 0=>OK.                             */
/*                     1 = bad type                                    */
/*                     2 = bad ndim                                    */
/*                     3 = bad dim                                     */
/*                     4 = bad objnum                                  */
/*                     5 = Object doesn't exist                        */
/*                     6 = Error creating storage                      */
/*---------------------------------------------------------------------*/
{
    extern struct object_dir task_dir;   /* task object directory */
    int  i, jerr, ok, len;
    char *p;
    struct AIPS_object *ao;
    struct ll_element *select, *last;
/*---------------------------------------------------------------------*/
    if ((*objnum<=0) | (*objnum>task_dir.numobj)) {*ierr = 4; return;}
    ao = task_dir.objpnt[*objnum-1];
    if (ao==NULL) {*ierr = 5; return;}
                                       /* Does it exist?               */
    zofnle_ (objnum, keyword, &jerr);
    if (jerr==0) {        /* previously exists, check compatability */
        select = ao->select;
        if (!(select->type==*type)) {*ierr = 1; return;}
        if (!(select->ndim==*ndim)) {*ierr = 2; return;}
        ok = 1;        len = *ndim;
        for (i=0; i<len; i++) ok = ok & (select->dim[i]==dim[i]);
        if (!ok) {*ierr = 3; return;}
    }
    else {   /* Create, add to end of list */
        select = create_llelem_ (keyword, type, ndim, dim, &jerr);
        if (!(jerr==0)) {*ierr = 6; return;} /* check */
        ao->select = select;
                                       /* Update linked list         */
        if (ao->first==NULL) ao->first = select;
        if (!(ao->last==NULL))  /* update previous last element if any */
            {last = ao->last;
             last->next = select;}
        ao->last = select;
    }
                                       /* store data                 */
    len = select->lendata;
    p = (char *) select + select->dataoff;
    for (i=0; i<len; i++) *(p++) = data[i];
    *ierr = 0;
    return;
}

#if __STDC__
   void zofedt_ (int *objnum, char *keyword, int *type,
                      int *ndim, int *dim, char *data, int *ierr)
#else
   void zofedt_ (objnum, keyword, type, ndim, dim, data, ierr)
   int *objnum, *type, *ndim, dim[], *ierr;
   char keyword[], data[];
#endif
/*---------------------------------------------------------------------*/
/*! Fetch data from AIPS Object linked list                            */
/*# OOPS                                                               */
/*---------------------------------------------------------------------*/
/* Fetch data from AIPS Object linked list.                            */
/*   Inputs:                                                           */
/*      objnum   I     Object number                                   */
/*      keyword  H*8   keyword name                                    */
/*   Outputs:                                                          */
/*      type     I     keyword type                                    */
/*      ndim     I     number of dimensions                            */
/*      dim      I     Dimensionality array                            */
/*      data     char  data array                                      */
/*      ierr     I     Return code, 0=>OK.                             */
/*                     1 = Not found                                   */
/*                     4 = bad objnum                                  */
/*                     5 = Object doesn't exist                        */
/*---------------------------------------------------------------------*/
{
    extern struct object_dir task_dir;   /* task object directory */
    int  i, jerr, len;
    char *p;
    struct AIPS_object *ao;
    struct ll_element *select;
/*---------------------------------------------------------------------*/
    if ((*objnum<=0) | (*objnum>task_dir.numobj)) {*ierr = 4; return;}
    ao = task_dir.objpnt[*objnum-1];
    if (ao==NULL) {*ierr = 5; return;}
                                       /* Find linked list element    */
    zofnle_ (objnum, keyword, &jerr);
    if (jerr==0) {        /* found */
        select = ao->select;
        *type = select->type;
        *ndim = select->ndim;
        for (i=0; i<*ndim; i++) dim[i] = select->dim[i];
                                       /* Copy data                  */
        len = select->lendata;
        p = (char *) select + select->dataoff;
        for (i=0; i<len; i++) data[i] = *(p++);
        *ierr = 0;
    }
    else {   /* Not found - error */
        *ierr = 1;
    }
    return;
}

#if __STDC__
   void zostct_ (int *objnum, int *catblk, int *ierr)
#else
   void zostct_ (objnum, catblk, ierr)
   int *objnum, catblk[], *ierr;
#endif
/*---------------------------------------------------------------------*/
/*! Store CATBLK to AIPS Object                                        */
/*# OOPS                                                               */
/*---------------------------------------------------------------------*/
/* Store CATBLK to AIPS Object.                                        */
/*   Inputs:                                                           */
/*      objnum   I       Object number                                 */
/*      catblk   I(256)   CATBLK array                                 */
/*   Outputs:                                                          */
/*      ierr     I     Return code, 0=>OK.                             */
/*                     4 = bad objnum                                  */
/*                     5 = Object doesn't exist                        */
/*---------------------------------------------------------------------*/
{
    extern struct object_dir task_dir;   /* task object directory */
    int  i;
    struct AIPS_object *ao;
/*---------------------------------------------------------------------*/
    if ((*objnum<=0) | (*objnum>task_dir.numobj)) {*ierr = 4; return;}
    *ierr = 0;
    ao = task_dir.objpnt[*objnum-1];
    if (ao==NULL) {*ierr = 5; return;}
    for (i=0; i<256; i++) ao->catblk[i] = catblk[i];
    return;
}

#if __STDC__
   void zofect_ (int *objnum, int *catblk, int *ierr)
#else
   void zofect_ (objnum, catblk, ierr)
   int *objnum, catblk[], *ierr;
#endif
/*---------------------------------------------------------------------*/
/*! Fetch CATBLK from AIPS Object                                      */
/*# OOPS                                                               */
/*---------------------------------------------------------------------*/
/* Fetch CATBLK from AIPS Object.                                      */
/*   Inputs:                                                           */
/*      objnum   I        Object number                                */
/*   Outputs:                                                          */
/*      catblk   I(256)   CATBLK array                                 */
/*      ierr     I     Return code, 0=>OK.                             */
/*                     4 = bad objnum                                  */
/*                     5 = Object doesn't exist                        */
/*---------------------------------------------------------------------*/
{
    extern struct object_dir task_dir;   /* task object directory */
    int  i;
    struct AIPS_object *ao;
/*---------------------------------------------------------------------*/
    if ((*objnum<=0) | (*objnum>task_dir.numobj)) {*ierr = 4; return;}
    *ierr = 0;
    ao = task_dir.objpnt[*objnum-1];
    if (ao==NULL) {*ierr = 5; return;}
    for (i=0; i<256; i++) catblk[i] = ao->catblk[i];
    return;
}

#if __STDC__
   void zofenm_ (int *objnum, char *name, char *class, int *ierr)
#else
   void zofenm_ (objnum, name, class, ierr)
   int *objnum, *ierr;
   char name[], class[];
#endif
/*---------------------------------------------------------------------*/
/*! Fetch Object name and class from AIPS Object                       */
/*# OOPS                                                               */
/*---------------------------------------------------------------------*/
/* Fetch Object name and class from AIPS Object.                       */
/*   Inputs:                                                           */
/*      objnum   I     Object number                                   */
/*   Outputs:                                                          */
/*      name    H*32   Object name                                     */
/*      class   H*8    Object class                                    */
/*      ierr     I     Return code, 0=>OK.                             */
/*                     4 = bad objnum                                  */
/*                     5 = Object doesn't exist                        */
/*---------------------------------------------------------------------*/
{
    extern struct object_dir task_dir;   /* task object directory */
    struct AIPS_object *ao;
    int    i;
/*---------------------------------------------------------------------*/
    if ((*objnum<=0) | (*objnum>task_dir.numobj)) {*ierr = 4; return;}
    *ierr = 0;
    ao = task_dir.objpnt[*objnum-1];
    if (ao==NULL) {*ierr = 5; return;}
    for (i=0; i<32; i++) name[i] = ao->name[i];
    for (i=0; i<8; i++) class[i] = ao->class[i];
    return;
}
/*
   Functions called from C:
*/

#if __STDC__
   struct ll_element* create_llelem_ (char *keyword, int *type,
                                     int *ndim, int *dim, int *ierr)
#else
   struct ll_element* create_llelem_ (keyword, type, ndim, dim, ierr)
   char keyword[];
   int *type, *ndim, dim[], *ierr;
#endif
/*---------------------------------------------------------------------*/
/*! Create/initialize AIPS Object linked list element                  */
/*# OOPS                                                               */
/*---------------------------------------------------------------------*/
/* Create/initialize AIPS Object linked list element                   */
/* Returns pointer to object; pointer to next set to void.             */
/*   Inputs:                                                           */
/*      keyword  H*8   keyword name                                    */
/*      type     I     keyword type                                    */
/*      ndim     I     number of dimensions                            */
/*      dim      I     Dimensionality array                            */
/*   Outputs:                                                          */
/*      ierr     I     Error code, 0=OK else failed.                   */
/*                     1 = memory allocation failed                    */
/*                     2 = illegal data type code                      */
/*---------------------------------------------------------------------*/
{
    struct ll_element *this;
    int  i, offset, len, el;
    size_t size, sizeh, sized;
/*---------------------------------------------------------------------*/
    *ierr = 0;
                                       /* Determine size of header     */
    sizeh = (4 + *ndim) * sizeof (int ) + 8 * (sizeof (char)) +
        sizeof(this);
                                       /* Determine size of data       */
    len = 1; i = 0;
    while ((i < *ndim) && (dim[i] > 0))
    {
       len = len * dim[i];
       i++;
    }
    el = 0;
    if (*type==1) el = sizeof(double);  /* double */
    if (*type==2) el = sizeof(float);   /* real */
    if (*type==3) el = sizeof(char);    /* Hollerith */
    if (*type==4) el = sizeof(int );    /* int */
    if (*type==5) el = sizeof(int );    /* fortran logical */
    if (el==0) {*ierr = 2; return NULL;} /* illegal type */
    sized = len * el;
    size = sizeh + sized;
                                       /* Allocate                    */
    if ((this = (struct ll_element *) malloc(size))==NULL)
        {*ierr = 1; return NULL;}      /* allocation failed */
                                       /* Initialize                  */
    for (i=0; i<8; i++) this->keyword[i] = keyword[i];
    this->next = NULL;
    this->dataoff = sizeh;
    this->lendata = sized;
    this->type = *type;
    this->ndim = *ndim;
    for (i=0; i<*ndim; i++) this->dim[i] = dim[i];
    return this;
}



