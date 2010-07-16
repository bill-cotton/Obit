/* $Id$       */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2009                                          */
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
#include "Obit.h"
#include "ObitUVSel.h"
#include "ObitTableNX.h"
#include "ObitTableSU.h"
#include "ObitTableSUUtil.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitUVSel.c
 * ObitUVSel Obit uv data selector class definition.
 * This contains information about data selection and calibration.
 */

/*--------------- File Global Variables  ----------------*/
/** name of the class defined in this file */
static gchar *myClassName = "ObitUVSel";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo global structure ObitIOClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitUVSelClassInfo myClassInfo = {FALSE};

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitUVSelInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitUVSelClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitUVSelClassInfoDefFn (gpointer inClass);

/** Private: Suggest subscan length. */
static void SuggestSubScan (ObitUVSel* in, ofloat scanTime);

/*---------------Public functions---------------------------*/
/**
 * Construct Object.
 * \return pointer to object created.
 */
ObitUVSel* newObitUVSel (gchar *name)
{
  ObitUVSel* out;

  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitUVSelClassInit();

  /* allocate structure */
  out = g_malloc0(sizeof(ObitUVSel));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

 /* set classInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitUVSelInit((gpointer)out);

  return out;
} /* end newObitUVSel */

/**
 * Returns ClassInfo pointer for the class.
 * Initializes class if needed on first call.
 * \return pointer to the class structure.
 */
gconstpointer ObitUVSelGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) 
    ObitUVSelClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitUVSelGetClass */

/**
 * Copy constructor.
 * \param in Pointer to object to be copied.
 * \param out Pointer to object to be written.  
 *            If NULL then a new structure is created.
 * \param err ObitErr error stack
 * \return Pointer to new object.
 */
ObitUVSel* ObitUVSelCopy (ObitUVSel* in, ObitUVSel* out, 
			  ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  gboolean oldExist;
  gchar *outName;
  olong i;

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
    out = newObitUVSel(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /* This class members */
  out->FileType    = in->FileType;
  out->nVisPIO     = in->nVisPIO;
  out->lrecUC      = in->lrecUC;
  out->nrparmUC    = in->nrparmUC;
  out->Compress    = in->Compress;
  out->numberVis   = in->numberVis;
  out->numberPoln  = in->numberPoln;
  out->jincs       = in->jincs;
  out->startChann  = in->startChann;
  out->numberChann = in->numberChann;
  out->jincf       = in->jincf;
  out->startIF     = in->startIF;
  out->numberIF    = in->numberIF;
  out->jincif      = in->jincif;
  out->doCalSelect = in->doCalSelect;
  out->transPol    = in->transPol;
  out->bothCorr    = in->bothCorr;
  out->timeRange[0]= in->timeRange[0];
  out->timeRange[1]= in->timeRange[1];
  out->UVRange[0]  = in->UVRange[0];
  out->UVRange[1]  = in->UVRange[1];
  out->FreqID      = in->FreqID;
  out->selectAnts  = in->selectAnts ;
  out->doPolCal    = in->doPolCal;
  out->doBLCal     = in->doBLCal;
  out->BLversion   = in->BLversion;
  out->doBPCal     = in->doBPCal;
  out->doBand      = in->doBand;
  out->BPversion   = in->BPversion;
  out->corrType    = in->corrType;
  out->doCal       = in->doCal;
  out->doCalWt     = in->doCalWt;
  out->calVersion  = in->calVersion;
  out->doFlag      = in->doFlag;
  out->FGversion   = in->FGversion;
  out->passAll     = in->passAll;
  out->alpha       = in->alpha;
  for (i=0; i<5; i++) out->Stokes[i] = in->Stokes[i];
  for (i=0; i<3; i++) out->smooth[i] = in->smooth[i];

  /* (de)Selected antenna list */
  out->selectAnts = in->selectAnts;
  if ((in->ants!=NULL) && (in->numberAntList>0)) {
    if (out->ants) g_free(out->ants);
    out->ants = g_malloc (in->numberAntList*sizeof(olong));
    out->numberAntList = in->numberAntList;
    for (i=0; i<in->numberAntList; i++) out->ants[i] = in->ants[i];
  }

  /* (de)Selected source list */
  out->selectSources = in->selectSources;
  if ((in->sources!=NULL) && (in->numberSourcesList>0)) {
    if (out->sources) g_free(out->sources);
    out->sources = g_malloc (in->numberSourcesList*sizeof(olong));
    out->numberSourcesList = in->numberSourcesList;
    for (i=0; i<in->numberSourcesList; i++) out->sources[i] = in->sources[i];
  }

  return out;
} /* end ObitUVSelCopy */

/**
 * Determines how large a buffer (in floats) is needed
 * for data transfers as described by data members.
 * The buffer is intended for the uncompressed versions
 * of uv data records.
 * \param desc Pointer input descriptor.
 * \param sel UV selector.
 * \return size in floats needed for I/O.
 */
olong ObitUVSelBufferSize (ObitUVDesc* desc, 
			       ObitUVSel* sel)
{
  olong size = 0;

  /* error checks */
  if (desc==NULL) return size; 
  g_assert (ObitIsA(desc, ObitUVDescGetClass()));
  g_assert (ObitIsA(sel, &myClassInfo));

  /* make sure defaults filled in */
  ObitUVSelDefault (desc, sel);

  /* size of uncompressed vis * number of vis */
  size = sel->lrecUC * sel->nVisPIO;

  return size;
} /* end ObitUVSelBufferSize */

/**
 * Enforce defaults (RA in range (0,360)
 * Also indexes structure.
 * \param in Pointer to descriptor.
 * \param sel UV selector, output vis descriptor changed if needed.
 */
void ObitUVSelDefault (ObitUVDesc* in, ObitUVSel* sel)
{

  /* error checks */
  g_assert (ObitIsA(in, ObitUVDescGetClass()));
  g_assert (ObitIsA(sel, &myClassInfo));


  /* Index as well */
  ObitUVDescIndex(in);

  /* Patch AIPS++ Bugs */
  if (in->jlocr>=0) {
    if (in->crval[in->jlocr]<0.0) in->crval[in->jlocr] += 360.0;
  }
  if (in->obsra<0.0) in->obsra += 360.0;
} /* end ObitUVSelDefault */

/**
 * Derive the descriptor for data being written; 
 * also updates defaults on sel.
 * \param in Pointer to input descriptor, this describes the data
 *           as they appear in memory.
 * \param sel UV selector, blc, trc members changed if needed.
 * \param out Pointer to output descriptor, describing form on disk.
 * \param err Obit error stack
 */
void ObitUVSelGetDesc (ObitUVDesc* in, ObitUVSel* sel,
		       ObitUVDesc* out, ObitErr *err)
{

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, ObitUVDescGetClass()));
  g_assert (ObitIsA(sel, &myClassInfo));
  g_assert (ObitIsA(out, ObitUVDescGetClass()));

  /* make sure defaults filled in */
  ObitUVSelDefault (in, sel);

  /* copy most values */
  ObitUVDescCopy (in, out, err);
  if (err->error) /* add traceback, return on error */
      Obit_traceback_msg (err, "ObitUVSelGetDesc", in->name);

  /* Save values  for compression */
  sel->nrparmUC = out->nrparm;
  sel->lrecUC   = out->lrec;

  /* If only one source selected make sure no "SOURCE" 
     random parameter is written */
  if ((sel->numberSourcesList==1) && (out->ilocsu>=0) )
    strncpy (out->ptype[out->ilocsu], "REMOVED ", UVLEN_KEYWORD); 

  /* compress iff sel->Compress */
  if (sel->Compress) {
    out->inaxes[0] = 1; /* not quite true for but gets float count
			  correct */
    /* Make sure there are WEIGHT and SCALE random parameters */
    if (out->ilocws<0) {
      out->ilocws = in->nrparm;
      strncpy (out->ptype[out->nrparm++], "WEIGHT  ", 8);
      strncpy (out->ptype[out->nrparm++], "SCALE   ", 8);
    }
  }

  /* update output descriptor for selection done in ObitUVCal::ObitUVCalSelectInit */
  /* make sure defaults, indices filled in */
  ObitUVSelDefault (in, sel);
  ObitUVSelDefault (out, sel);

} /* end ObitUVSelGetDesc */

/**
 * Apply selection criteria to input descriptor to derive output.
 * Note: many operations associated with data selection are done in
 * ObitUVCalSelectInit.
 * Also sets previously undefined values on sel.
 * \param in Pointer to input descriptor, this describes the data
 *           as they appear on disk (possibly compressed).
 * \param sel UV selector, members changed if needed.
 * \param out Pointer to output descriptor, this describes the data 
 *            after any processing when read, or before any compression
 *            on output.
 * \param err Obit error stack
 */
void ObitUVSelSetDesc (ObitUVDesc* in, ObitUVSel* sel,
			  ObitUVDesc* out, ObitErr *err)
{

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, ObitUVDescGetClass()));
  g_assert (ObitIsA(sel, &myClassInfo));
  g_assert (ObitIsA(out, ObitUVDescGetClass()));

  /* make sure defaults filled in */
  ObitUVSelDefault (in, sel);

  /* copy most values */
  ObitUVDescCopy (in, out, err);
  if (err->error) /* add traceback, return on error */
      Obit_traceback_msg (err, "ObitUVSelSetDesc", 
			  in->name);

  /* if reading compressed data it will be uncompressed. */
  if ((in->access==OBIT_IO_ReadOnly) ||  (in->access==OBIT_IO_ReadWrite) 
      ||(in->access==OBIT_IO_ReadCal)) {
    if ((in->inaxes[0]==1) || (in->inaxes[0]==2)) {
      sel->Compress = TRUE;
      /* If they're last drop WEIGHT and SCALE random parameters */
      if (out->ilocws == out->nrparm-2) {
	out->nrparm -= 2;
	out->ilocws = -1;
	out->inaxes[0] = 3;
      }
      out->lrec =  out->nrparm + (in->lrec - in->nrparm) * 3;
    }
  }

  /* Save this size for decompression */
  sel->nrparmUC = out->nrparm;
  sel->lrecUC   = out->lrec;

  /* make sure defaults, indices filled in */
  ObitUVSelDefault (in, sel);
  ObitUVSelDefault (out, sel);

  /* Save Selector values */
  /* set data increments as float. */
  sel->jincs  = in->incs;
  sel->jincf  = in->incf;
  sel->jincif = in->incif;

  /* Selection */
  if (sel->numberChann<=0) sel->numberChann = in->inaxes[in->jlocf];
  if (sel->startChann<=0)  sel->startChann = 1;
  sel->numberChann = MIN (sel->numberChann, in->inaxes[in->jlocf]);

  if ((sel->numberIF<=0) && (in->jlocif>=0)) 
    sel->numberIF = in->inaxes[in->jlocif];
  if (sel->numberIF<=0) sel->numberIF = 1;
  if (sel->startIF<=0)  sel->startIF = 1;
  if (in->jlocif>=0) 
    sel->numberIF = MIN (sel->numberIF, in->inaxes[in->jlocif]);

} /* end ObitUVSelSetDesc */

/**
 * See if an NX table exists and if so initialize it to use in deciding
 * which visibilities to read.
 * \param  in      Pointer to the object.
 * \param  desc    UV descriptor from IO where the next visibility to
 *                 read and the number will be stored.
 * \param  err     Error stack
 * \return TRUE is finished, else FALSE
 */
void ObitUVSelNextInit (ObitUVSel *in, ObitUVDesc *desc, ObitErr *err)
{
  ObitIOCode retCode;
  gchar *routine="ObitUVSelNextInit";
 
 /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitUVSelIsA(in));
  g_assert (ObitUVDescIsA(desc));


  /* Indexing? */
  if (in->NXTable==NULL) return;

  /* Open Index table  */
  retCode = 
    ObitTableNXOpen ((ObitTableNX*)(in->NXTable), OBIT_IO_ReadOnly, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_msg (err, routine, in->name);
  in->numRow = ((ObitTableNX*)in->NXTable)->myDesc->nrow;
  in->LastRowRead = 0;
  
   /* Initialize */
  in->scanFirstVis = -1;
  in->scanLastVis  = -1;

  /* Create row structure */
  in->NXTableRow = (Obit*)newObitTableNXRow((ObitTableNX*)(in->NXTable));

  /* We want indexing */
  in->doIndex = TRUE; 
  return;
} /* end ObitUVSelNextInit */

/**
 * Uses selector member to decide which visibilities to
 * read next.
 * If doIndex is TRUE, then visibilities are selected from the NX table.
 * \param  in      Pointer to the object.
 * \param  desc    UV descriptor from IO where the next visibility to
 *                 read and the number will be stored.
 *                 0 causes an initialization.&nleft)
 * \param  err     Error stack
 * \return TRUE is finished, else FALSE
 */
gboolean ObitUVSelNext (ObitUVSel *in, ObitUVDesc *desc, ObitErr *err)
{
  ObitIOCode retCode = OBIT_IO_SpecErr;
  ObitTableNXRow *row;
  olong nleft;
  gboolean gotIt, done = FALSE;
  gchar *routine = "ObitUVSelNext";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return done;
  g_assert (ObitUVSelIsA(in));
  g_assert (ObitUVDescIsA(desc));
  
  /* Is this controlled by an index? */
  if (in->doIndex) {
    /* Which visibilities are wanted? */
    if (desc->firstVis < 1) { /* first read? */
      desc->firstVis = 1;
    } else { /* subsequent reads */
      desc->firstVis += in->numVisRead;
    }
    
    /* Need a new scan? */
    if (desc->firstVis>in->scanLastVis) {
      
      /* Read index file until a selected scan found */
      row = (ObitTableNXRow*)in->NXTableRow;
      gotIt = FALSE;
      while ((!gotIt) && (in->LastRowRead<in->numRow)) {
	in->LastRowRead++;
	retCode = ObitTableNXReadRow ((ObitTableNX*)in->NXTable, in->LastRowRead, row, err);
	if (err->error) Obit_traceback_val (err, routine, in->name, retCode);
	
	/* Is this one wanted? */
	if ((in->SubA > 0) && (in->SubA != row->SubA)) continue;
	if ((in->FreqID > 0) && (in->FreqID != row->FreqID)) continue;
	
	/* Any overlap with time range? */
	if (row->Time+row->TimeI < in->timeRange[0]) continue;
	if (row->Time-row->TimeI > in->timeRange[1]) break;
	
	/* Requested source? */
	gotIt = ObitUVSelWantSour(in, row->SourID);
      } /* end loop over file */
      
	/* Save info if scan found */
      if (gotIt) {
	in->scanFirstVis = row->StartVis;
	in->scanLastVis  = row->EndVis;
	desc->firstVis   = in->scanFirstVis; /* beginning of scan */
	/* Suggest subscan length */
	SuggestSubScan(in, row->TimeI);
      } else {
	/* No more scans - Must be finished */
	done = TRUE;
	return done;
      }
    } /* end get new scan */
    
    /* how many to attempt to read? */
    in->numVisRead = in->nVisPIO;
    
    /* but not more than all of scan */
    nleft = in->scanLastVis - desc->firstVis + 1;
    in->numVisRead  = MIN (nleft, in->numVisRead);
    in->numVisRead  = MAX (0, in->numVisRead);
    done = (nleft<=0);
    
  } else {
    /* not indexed */
    /* Which visibilities are wanted? */
    if (desc->firstVis < 1) { /* first read? */
      desc->firstVis = 1;
    } else { /* subsequent reads */
      desc->firstVis += in->nVisPIO;
    }
    
    /* how many? */
    in->numVisRead = in->nVisPIO;
    
    /* but not more than all */
    nleft = desc->nvis - desc->firstVis + 1;
    in->numVisRead = MIN (nleft, in->numVisRead);
    in->numVisRead = MAX (0, in->numVisRead);
    done = (nleft<=0);
  }
  
  return done;
} /* end ObitUVSelNext */

/**
 * Close NX table if open .
 * If doIndex is TRUE, then visibilities are selected from the NX table.
 * \param  in   Pointer to the Selector.
 * \param  err  Error stack
 * \return TRUE is finished, else FALSE
 */
void ObitUVSelShutdown (ObitUVSel *in, ObitErr *err)
{
  ObitIOCode retCode;
  gchar *routine="ObitUVSelShutdown";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitUVSelIsA(in));

  /* Anything to do? */
  if (!in->doIndex) return;

  /* Close table  */
  retCode = ObitTableNXClose ((ObitTableNX*)in->NXTable, err);
  if ((retCode!=OBIT_IO_OK) || (err->error)) /* add traceback,return */
    Obit_traceback_msg (err, routine, in->name);

  /* Release structures  */
  in->NXTable    = ObitTableNXUnref(in->NXTable);
  in->NXTableRow = ObitTableNXRowUnref(in->NXTableRow);
  in->doIndex    = FALSE;  /* Indexing no longer wanted */

} /* end ObitUVSelShutdown */

/**
 * Set selector for source selection
 * \param sel      UV selector.
 * \param inData   Associated UV data (as gpointer to avoid recursive definition)
 * \param Qual     Source qualifier, -1 => any 
 * \param souCode  selection of Source by Calcode,if not specified in Source
 *                  '    ' => any calibrator code selected  
 *                  '*   ' => any non blank code (cal. only)
 *                  '-CAL' => blank codes only (no calibrators)           
 *                  anything else = calibrator code to select.            
 *               NB: The souCode test is applied in addition to the       
 *               other tests, i.e. Sources and Qual, in the               
 *               selection of sources to process
 * \param Sources  Selected source names, [0] blank=> any
 *                 this is passed as a lsou x nsou array of characters
 * \param lsou     length of source name in Sources 
 * \param nsou     maximum number of entries in Sources
 * \param err      Obit error/message stack
 */
void ObitUVSelSetSour (ObitUVSel* sel, gpointer inData, olong Qual, 
		       gchar *souCode, gchar *Sources, olong lsou, olong nsou, 
		       ObitErr *err)
{
  ObitTableSU *SUTable=NULL;
  gint32 idim[MAXINFOELEMDIM];
  olong i, j, count;
  olong iver;
  gchar *name = "Input Data";
  gchar *routine = "ObitUVSelSetSour";

  /* Set Selected sources  */
  if ((nsou>0) && Sources) {
    /* Count actual entries in source list */
    count = 0;  j = 0;
    for (i=0; i<nsou; i++) {
      if ((Sources[j]!=' ') || (Sources[j+1]!=' ')) count++;
      j += lsou;
    }
    sel->numberSourcesList = count;
    if (count>0) {  /* Anything actually specified? */
      /* have to lookup sources - need SU table for this. */
      iver = 1;
      SUTable = newObitTableSUValue (name, (ObitData*)inData, &iver, 
				     OBIT_IO_ReadOnly, 0, err);
      if (SUTable==NULL) {  /* No source table - only one source and it is selected */
	sel->numberSourcesList = 0;
      } else { /* Lookup sources to get numbers */
	/* In case all selected */
	/* Open table */
	ObitTableSUOpen (SUTable, OBIT_IO_ReadOnly, err);
	if (err->error)  Obit_traceback_msg (err, routine, name);
	sel->numberSourcesList = MAX (count, SUTable->myDesc->nrow); 
	ObitTableSUClose (SUTable, err);
	if (err->error)  Obit_traceback_msg (err, routine, name);
	if (sel->sources)
	  sel->sources = g_realloc(sel->sources, sel->numberSourcesList*sizeof(olong));
	else
	  sel->sources = g_malloc0(sel->numberSourcesList*sizeof(olong));
	
	/* Do lookup */
	idim[0] = lsou; idim[1] = count;
	ObitTableSULookup (SUTable, idim, Sources, Qual, souCode, 
			   sel->sources, &sel->selectSources, 
			   &sel->numberSourcesList, err); 
	if (err->error)  Obit_traceback_msg (err, routine, name);
	SUTable = ObitTableSUUnref(SUTable); /* release table */
      }
    if (err->error) Obit_traceback_msg (err, routine, name);
    } else { /* end of sources specified */
      /* Deallocate sources if needed */
      if (sel->sources) {g_free(sel->sources); sel->sources=NULL;}
    }
  } else { /* no "Sources" specified */
    sel->numberSourcesList = 0; /* everything selected */
  }
} /* end ObitUVSelSetSour */

/**
 * Set selector for antenna selection
 * \param sel    UV selector.
 * \param Antennas List of selected Antennas, NULL 
 *                 or all 0 => all, zero entries after
 *                 first non zero are ignored.
 *                 Any negative values means all named 
 *                 are deselected
 * \param nant     Number of entries in Antennas
 */
void ObitUVSelSetAnt (ObitUVSel* sel, olong *Antennas, olong nant)
{
  olong i;

  /* Setup selector on inData Antennas, nant */
  if ((nant>0) && Antennas) {
   sel->numberAntList = nant;
    sel->ants = g_realloc(sel->ants, sel->numberAntList*sizeof(olong));
    /* loop copying, checking for deselection */
    sel->selectAnts = FALSE;
    for (i=0; i<sel->numberAntList; i++) {
      sel->ants[i] = abs (Antennas[i]);
      sel->selectAnts = sel->selectAnts || (sel->ants[i]>0);
    }
    /* If not selecting by antenna free ants */
    if (!sel->selectAnts) {
      sel->numberAntList = 0;
      g_free(sel->ants); sel->ants = NULL;
    }
  } else {  /* No selection */
    sel->selectAnts = FALSE;
    sel->numberAntList = 0;
  }

} /* end ObitUVSelSetAnt */


/**
 * Determine if a given source is selected.
 * \param sel    UV selector.
 * \param SourID Source ID to be tested
 * \return TRUE if source selected.
 */
gboolean ObitUVSelWantSour (ObitUVSel* sel, olong SourID)
{
  olong i;

  /* error checks */
  g_assert (ObitIsA(sel, &myClassInfo));

  /* If array is null - everything selected */
  if (sel->sources == NULL) return TRUE;

  /* ditto no entries */
  if (sel->numberSourcesList == 0) return TRUE;

  /* Check if explicitly selected */
  if (sel->selectSources) {
    for (i=0; i<sel->numberSourcesList; i++) {
      if (SourID==sel->sources[i]) return TRUE;
    }
    /* didn't make the cut */
    return FALSE;

  } else {
    /* check is explicitly deselected */
    for (i=0; i<sel->numberSourcesList; i++) {
      if (SourID==sel->sources[i]) return FALSE;
    }
    /* Survived */
    return TRUE;
  }

  return TRUE; /* shouldn't get here */
} /* end ObitUVSelWantSour */

/**
 * Determine if a given antenna is selected.
 * \param sel    UV selector.
 * \param ant    antenna id to test
 * \return TRUE if antenna selected.
 */
gboolean ObitUVSelWantAnt (ObitUVSel* sel, olong ant)
{
  olong i;

  /* error checks */
  g_assert (ObitIsA(sel, &myClassInfo));

  /* If array is null - everything selected */
  if (sel->ants == NULL) return TRUE;

  /* ditto no entries */
  if (sel->numberAntList == 0) return TRUE;

  /* Check if explicitly selected */
  if (sel->selectAnts) {
    for (i=0; i<sel->numberAntList; i++) {
      if (ant==sel->ants[i]) return TRUE;
    }
    /* didn't make the cut */
    return FALSE;

  } else {
    /* check is explicitly deselected */
    for (i=0; i<sel->numberAntList; i++) {
      if (ant==sel->ants[i]) return FALSE;
    }
    /* Survived */
    return TRUE;
  }

  return TRUE; /* shouldn't get here */
} /* end ObitUVSelWantAnt */

/**
 * Suggest a length for a sub interval of the current scan 
 * such that the scan is evenly divided.  This is based on the
 * target value SubScanTime.
 * This is only useful for ReadSelect operations on an indexed ObitUV.
 * \param sel    UV selector.
 * \return suggested subscan length in days; 
 */
ofloat ObitUVSelSubScan (ObitUVSel* sel)
{
  return sel->SubScanSuggest;
} /* end ObitUVSelSubScan */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitUVSelClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitUVSelClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitUVSelClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitUVSelClassInfoDefFn (gpointer inClass)
{
  ObitUVSelClassInfo *theClass = (ObitUVSelClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitUVSelClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitUVSelClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitUVSelGetClass;
  theClass->newObit       = (newObitFP)newObitUVSel;
  theClass->ObitCopy      = (ObitCopyFP)ObitUVSelCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitUVSelClear;
  theClass->ObitInit      = (ObitInitFP)ObitUVSelInit;

} /* end ObitUVSelClassDefFn */

/*---------------Private functions--------------------------*/
/**
 * Creates empty member objects, initialize reference count.
 * Does (recursive) initialization of base class members before 
 * this class.
 * \param inn Pointer to the object to initialize.
 */
void ObitUVSelInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ofloat fblank = ObitMagicF();
  ObitUVSel *in = inn;

  /* error checks */
  g_assert (in != NULL);
  
  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->nVisPIO       = 1;
  in->numberVis     = 1;
  in->numberPoln    = 1;
  in->startChann    = 1;
  in->numberChann   = 1;
  in->jincf         = 3;
  in->startIF       = 1;
  in->numberIF      = 1;
  in->jincif        = 3;
  in->selectAnts    = TRUE;
  in->ants          = NULL;
  in->selectSources = TRUE;
  in->sources       = NULL;
  in->timeRange[0]  = -1.0e20;
  in->timeRange[1]  = 1.0e20;
  in->UVRange[0]    = 0.0;
  in->UVRange[1]    = 1.0e20;
  in->smooth[0]     = 0.0;
  in->alpha         = fblank;
  in->corrType      = 1;
  in->doBand        = FALSE;
  in->doCal         = FALSE;
  in->doCalWt       = FALSE;
  in->doFlag        = FALSE;
  in->doPolCal      = FALSE;
  in->doBLCal       = FALSE;
  in->doBPCal       = FALSE;
  in->doCalSelect   = FALSE;
  in->transPol      = FALSE;
  in->bothCorr      = FALSE;
  in->doIndex       = FALSE;
  in->passAll       = FALSE;
  in->NXTable       = NULL;
  in->NXTableRow    = NULL;
  in->numRow        = -1;
  in->LastRowRead   = 0;
} /* end ObitUVSelInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 */
void ObitUVSelClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitUVSel *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* free this class members */
  if (in->ants) g_free(in->ants);       in->ants    = NULL;
  if (in->sources) g_free(in->sources); in->sources = NULL;
  in->NXTable    = ObitTableNXUnref(in->NXTable);
  in->NXTableRow = ObitTableNXRowUnref(in->NXTableRow);
  
  /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);

} /* end ObitUVSelClear */

/**
 * Determine a time interval close to in->SubScanTime which
 * evenly divide scanTime into sub scans and save value as
 * in->SubScanSuggest.
 * in->SubScanTime <= 0 means scan average
 * \param  in   Selector object
 */
void SuggestSubScan (ObitUVSel* in, ofloat scanTime)
{

  olong  npiece;

  /* Scan average? */
  if (in->SubScanTime<=0.0) {
    in->SubScanSuggest = scanTime;
    return;
  }

  /* Divide scanTime into equal pieces */
  npiece = 0.5 + scanTime / in->SubScanTime;
  npiece = MAX (1, npiece);
  in->SubScanSuggest = scanTime / npiece;
} /* end SuggestSubScan */

