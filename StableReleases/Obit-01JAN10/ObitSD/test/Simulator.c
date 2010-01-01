#include <stdio.h>
#include <stdlib.h>
#include "ObitAll.h"
#include "ObitOTF.h"
#include "ObitPennArrayUtil.h"
#include "ObitOTFSkyModel.h"
#include "ObitOTFUtil.h"
#include "ObitFArray.h"
#include "ObitFInterpolate.h"
#include "ObitTableOTFTarget.h"
#include "ObitTableOTFIndex.h"

/* Data simulator for GBT Penn Array */

/* Local function prototypes */
/* Define times and positions for a scan */
void DefineScan (ofloat *data, ObitOTFDesc *desc, ObitOTFArrayGeom *geom, olong lenScan, 
		 ofloat *timeStart, ofloat timeStep,
		 ofloat RACenter, ofloat DecCenter, ofloat PA, ofloat skyStep);

/* Define Sky model */
ObitOTFSkyModel* DefineSkyModel (void);

/* Create Atmospheric model */
ObitImage* DefineAtmosModel(ofloat RACenter, ofloat DecCenter, 
			    ofloat skyStep, olong size, ObitErr *err);

/* Get target id */
olong GetTarget (ObitOTF *outData, gboolean isNew, gchar *name, ObitErr *err);
/* Initialize Index table */
void InitScan (ObitOTF *outData, gboolean isNew, ofloat scan, ofloat target, 
	       ObitErr *err);
/* Undate scan info */
void SetScan (ObitOTF *outData, odouble startTime, odouble endTime, 
	      olong startRec, olong endRec, ObitErr *err);

/* Program globals */
#define MAXSCAN 1000       /* Maximum number of scans */
olong  nScan=0;            /* number of scans */
olong  startRec[MAXSCAN];  /* First record number (1-rel) in scan */
olong  endRec[MAXSCAN];    /* End record number (1-rel) in scan */
odouble startTime[MAXSCAN];/* Start time of Scan in days */
odouble endTime[MAXSCAN];  /* End time of Scan in days */
ofloat target;             /* target number */
ofloat scan;               /* scan number */


int main ( int argc, char **argv )
{
  ObitSystem *mySystem;
  ObitOTF *simData;
  ObitOTFDesc *desc;
  ObitOTFSkyModel* sky=NULL;
  ObitTableSkyModel *skyTable=NULL;
  ObitImage *atmos=NULL;
  ObitFInterpolate *atmosInt=NULL;
  ObitErr *err;
  gchar *Filename="GBTSimulate.fits";
  gchar *FITSdir[] = {"FITSdata/"};
  olong i, ver;
  gboolean isNew;
  ofloat *buffer;
  /* define simulation (suitable for 19 Sept 2000) */
  /* Do 5 passes in two orthogonal directions with no overlap among
     parallel passes */
  olong numScan    = 10;                  /* Number of scans */
  olong lenScan    = 283;                 /* Number of samples in scan, 
					     5 samples/beam * 8 beams * width (5) / cos(PA) */
  ofloat timeStart = 0.0;                 /* Initial UTC (days) */
  ofloat timeStep  = 0.000001;            /* Integration time in Days */
  ofloat skyStep   = (32.0/3600.0)/(8*5); /* motion on the sky in one integration 
					     array width = 32/3600 deg, with 8 beams, 5 samples/beam */
  ofloat RACenter  = 18.0 * 15.0;         /* RA of center (deg) */
  ofloat DecCenter = 50.0;                /* Dec of center (deg) */
  ofloat RAScan[]  = {                    /* Ra offset from RACenter of center of scan */
    -0.02, -0.01, 0.0, 0.01, 0.02, -0.02, -0.01, 0.0, 0.01, 0.02};
  ofloat DecScan[] = {                    /* Dec  offset from DecCenter of center of scan */
    0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  ofloat PAScan[]  = {                    /* Position angle (deg) of scan direction */
    45., 45., 45., 45., 45., -45., -45., -45., -45., -45.};
  ofloat raCen, decCen;

  /* Initialize Obit */
  err = newObitErr();
  mySystem = ObitSystemStartup ("Simulator", 1, 0, 0, NULL, 1, FITSdir, (oint)TRUE, (oint)FALSE, err);
  ObitErrLog(err); /* show any error messages on err */

  /* Create ObitOTF for data */
  simData = ObitPennArrayCreate("GBT Simulation");

  /* Define output, I/O size to one scan */
  ObitOTFSetFITS(simData,lenScan,1,Filename,err);

  /* Center of pattern */
  simData->myDesc->obsra  = RACenter;
  simData->myDesc->obsdec = DecCenter;

  /* Beam size (Gaussian FWHM) 7.2 asec. */
  simData->myDesc->beamSize = 7.2 / 3600.0;

  /* Create Sky model */
  sky = DefineSkyModel();
    
  /* Create atmospheric model */
  atmos =  DefineAtmosModel(RACenter, DecCenter, skyStep*6, lenScan, err);
  ObitErrLog(err); /* show any error messages on err */

  /* An interpolator for the model */
  atmosInt = newObitFInterpolateCreate (atmos->name, atmos->image, atmos->myDesc, 2);

  /* Open */
  if ((ObitOTFOpen (simData, OBIT_IO_WriteOnly, err) 
       != OBIT_IO_OK) || (err->error>0))  /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR opening output FITS file %s", Filename);
  /* show any errors */
  ObitErrLog(err);

  /* local pointers */
  buffer = simData->buffer;
  desc   = simData->myDesc;

  /* Get target id number */
  isNew = TRUE;
  target = (ofloat)GetTarget (simData, isNew, "GBT Simulation", err);

  /* Loop over scans */
  nScan = 0;
  for (i=0; i<numScan; i++) {

    /* Scan start info */
    scan = (float)nScan+1.0;
    startRec[nScan]  = desc->nrecord+1;
    startTime[nScan] = timeStart;

    /* Define times/geometry */
    decCen = DecCenter + DecScan[i];
    raCen  = RACenter + RAScan[i]/cos(decCen*1.74533e-2);
    DefineScan (buffer, desc, simData->geom, lenScan, &timeStart, timeStep,
		raCen, decCen , PAScan[i], skyStep);
    simData->myDesc->numRecBuff = lenScan; /* how many */

    /* Add sky model */
    ObitOTFUtilSubSkyModelBuff (simData, sky, -1.0);

    /* Add atmospheric model */
    /* debug  ObitOTFUtilSubImageBuff (simData, atmosInt, -1.0, err);*/
   
    /* Write scan */
    if ((ObitOTFWrite (simData, NULL, err) != OBIT_IO_OK) || (err->error>0))
      Obit_log_error(err, OBIT_Error, "ERROR reading input Table file");
    
    /* Scan end info */
    endRec[nScan]  = desc->nrecord+1;
    endTime[nScan] = timeStart;
    nScan++;

  } /* end scan loop */
  
  /* show any errors */
  ObitErrLog(err);

  /* Loop over scans */
  isNew = TRUE;
  for (i=0; i<nScan; i++) {
    scan = (float)i+1;
    /* Initialize scan in Index table */
    InitScan (simData, isNew, scan, target, err);
    isNew = FALSE;

    /* Update index table */
    SetScan (simData, startTime[i], endTime[i], startRec[i], endRec[i], err);
  }

  /* Close */
  if ((ObitOTFClose (simData, err) != OBIT_IO_OK) || (err->error>0))
    Obit_log_error(err, OBIT_Error, "ERROR closing output file");

  /* show any errors */
  ObitErrLog(err);

  /* Copy sky model to output */
  ver = 1;
  skyTable =  
    newObitTableSkyModelValue ("SkyModelTable", (Obit*)simData, &ver, 
			       OBIT_IO_WriteOnly, err);
  if ((ObitOTFSkyModelWrite (sky, skyTable, err) 
       != OBIT_IO_OK) || (err->error>0))
    Obit_log_error(err, OBIT_Error, "ERROR writing sky model");
  
  /* Clean up*/
  simData = ObitUnref(simData);
  sky = ObitUnref(sky);
  skyTable = ObitUnref(skyTable);

  /* show any errors */
  ObitErrLog(err);

  /* Shutdown Obit */
  mySystem = ObitSystemShutdown (mySystem);
  
  return 0;
} /* end of main */

/**
 * Fill in the times and sky pointing positions for a scan.
 * Assume a flat sky.
 * \param data       Data buffer, must be allocated large enough for lenScan records
 * \param desc       OTF descriptor for the records to be written
 * \param geom       OTF array geometry object
 * \param lenScan    Number of samples in scan
 * \param timeStart  Initial time, end time plus 10 sec will be returned.
 * \param timeStep   Integration time
 * \param RACenter   RA (deg) of center of scan.
 * \param DecCenter  Dec (deg) of center of scan
 * \param PA         Position angle on the sky of scan
 * \param skyStep    Angular seperation (deg) of samples along scan.
 */
void DefineScan (ofloat *data, ObitOTFDesc *desc, ObitOTFArrayGeom *geom, olong lenScan, 
		 ofloat *timeStart, ofloat timeStep,
		 ofloat RACenter, ofloat DecCenter, ofloat PA, ofloat skyStep)
{
  olong i, j, icen;
  ofloat ra, dec, deltaRA, deltaDec;

  /* center integration */
  icen = 1 + lenScan / 2;

  /* Spacing in RA, dec */
  deltaRA  = (skyStep * sin (PA*1.74533e-2)) / cos(DecCenter*1.74533e-2);
  deltaDec = skyStep * cos (PA*1.74533e-2);

  /* Test write a sequence of entries */
  for (i=0; i<lenScan; i++){

    /* position */
    ra  = RACenter  + (i-icen) * deltaRA;
    dec = DecCenter + (i-icen) * deltaDec;
      
    /* Fill record entry in data */
    data[desc->iloct]   = *timeStart;   /* Time */
    data[desc->ilocti]  = timeStep;     /* time interval */
    data[desc->iloctar] = target;       /* target number */
    data[desc->ilocscan]= scan;         /* scan number */
    data[desc->ilocra]  = ra;           /* RA  */
    data[desc->ilocdec] = dec;          /* Dec */
    data[desc->iloccal] = 0.0;          /* No cal yet */

    /* Parallactic angle */
    data[desc->ilocrot] = ObitOTFArrayGeomParAng(geom, *timeStart, ra, dec);

    /* zero data values */
    for (j=0; j<geom->numberDetect; j++) data[desc->ilocdata+j] = 0;

    /* increment time */
    *timeStart += timeStep;
    
    data += desc->lrec;
  } /* end loop filling data */
  
  *timeStart += 20.0/86400.0; /* Add 10 seconds for scan overhead */
  
} /* end DefineScan */

/**
 * Constructor for default GBT Penn Array Sky Model.
 * \param name An optional name for the object.
 * \return the new object.
 */
ObitOTFSkyModel* DefineSkyModel (void)
{
  /* Model must be compatable with data defined in DefineScan */
  ofloat RACenter  = 18.0 * 15.0;         /* RA of center (deg) */
  ofloat DecCenter = 50.0;                /* Dec of center (deg) */
  ObitOTFProj proj = OBIT_OTF_SIN;        /* -SIN projection */
  olong i, numComp = 4;                   /* Number of components */
  ofloat x[] = {                          /* X offsets (deg) */
  0.0, -0.01, 0.01, 0.01};
  ofloat y[] = {                          /* Y offsets (deg) */
  0.0, 0.01, -0.01, 0.01};
  ofloat f[] = {                          /* Flux density (Jy) */
   /* 0.000001, 0.000001, 0.000001, 0.000001};*/
    1.0, 2.0, 3.0, 4.0};

  ObitOTFSkyModel* out=NULL;
  
  /* Create basic structure */
  out = ObitOTFSkyModelCreate(numComp);

  /* Header */
  out->RACenter  = RACenter;
  out->DecCenter = DecCenter;
  out->proj      = proj;

  /* Model */
  for (i=0; i<numComp; i++) {
    out->RAOffset[i]  = x[i];
    out->DecOffset[i] = y[i];
    out->flux[i]      = f[i];
  }

 return out;
} /* end DefineSkyModel */

/**
 * Create atmospheric model as an Image.
 * \param RACenter   RA (deg) of center of scan.
 * \param DecCenter  Dec (deg) of center of scan
 * \param skyStep    Angular seperation (deg) of samples along scan.
 * \param size       Dimension in pixels of axes
 * \param  err       Error stack
 * \return the Image
 */
/* Create Atmospheric model */
ObitImage* DefineAtmosModel(ofloat RACenter, ofloat DecCenter, 
			    ofloat skyStep, olong size, ObitErr *err)
{
  ObitImage *atmos=NULL;
  ObitImageDesc *desc=NULL;
  olong blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  gchar *Outname="!Atmosphere.fits";
  olong naxis[2], i, j, k;
  ofloat *data, x, y;
  gchar *routine = "DefineAtmosModel";

  atmos = newObitImage("Atmospheric model");

  /* Create array */
  naxis[0] = size;
  naxis[1] = size;
  atmos->image = ObitFArrayCreate ("Model atmos", 2, naxis);

  /* Set Image Descriptor */
  desc = atmos->myDesc;
  desc->bitpix = -32;
  desc->naxis = 2;
  desc->inaxes[0] = naxis[0];
  desc->inaxes[1] = naxis[1];
  strncpy (desc->ctype[0], "RA---SIN", IMLEN_KEYWORD);
  strncpy (desc->ctype[1], "DEC--SIN", IMLEN_KEYWORD);
  desc->cdelt[0] = -skyStep;
  desc->cdelt[1] = skyStep;
  desc->crpix[0] = 1 + naxis[0] / 2;
  desc->crpix[1] = 1 + naxis[1] / 2;
  desc->crota[0] = 0.0;
  desc->crota[1] = 0.0;
  desc->crval[0] = RACenter; 
  desc->crval[1] = DecCenter;
  strncpy (desc->bunit, "KELVIN", IMLEN_VALUE);

  /* Create image */
  ObitImageSetFITS(atmos,OBIT_IO_byPlane,1,Outname,blc,trc,err);
  if (err->error) Obit_traceback_val (err, routine, atmos->name, atmos);

  /* Open/Create */
  if ((ObitImageOpen (atmos, OBIT_IO_WriteOnly, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, 
		   "ERROR opening image %s", atmos->name);
    return atmos;
  }

  /* Set image */
  naxis[0] = 0; naxis[1] = 0; 
  data = ObitFArrayIndex (atmos->image, naxis);
  for (k=0; k<atmos->image->arraySize; k++ ) {
    i = k % desc->inaxes[0];
    j = k / desc->inaxes[0];
    x = 8 * (i-desc->crpix[0]/2) / desc->inaxes[0];
    y = 5 * (j-desc->crpix[1]/2) / desc->inaxes[1];
    data[k] = 5.0 + sin (2.0*G_PI * (x * y));
  }

  /* Write image */
  ObitImageWrite (atmos, NULL, err);
  if (err->error) Obit_traceback_val (err, routine, atmos->name, atmos);
  
  /* Close Image */
  ObitImageClose (atmos, err);
  if (err->error) Obit_traceback_val (err, routine, atmos->name, atmos);
  
  return atmos;
} /* end DefineAtmosModel */

olong GetTarget (ObitOTF *outData, gboolean isNew, gchar *name, ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Get target id, look through existing table create new entry if needed.*/
/*  Returns Target id                                                     */
/*   Input:                                                               */
/*      outData  Output OTF object                                        */
/*      isNew    True if output file just created                         */
/*      name     Name of target                                           */
/*   Output:                                                              */
/*       err       Obit return error stack                                */
/*   Return:                                                              */
/*      Target id                                                         */
/*----------------------------------------------------------------------- */
{
  olong targ = -1;
  ObitTableOTFTarget* table;
  ObitTableOTFTargetRow* row;
  olong iRow, ver;
  gboolean doWrite;
  ObitIOAccess access;
  gchar *routine = "GetTarget";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return targ;
  g_assert (ObitOTFIsA(outData));
  g_assert(name!=NULL);

  /* create Target table object */
  ver = 1;
  if (isNew) access = OBIT_IO_WriteOnly;
  else access = OBIT_IO_ReadWrite;
  table = newObitTableOTFTargetValue ("Target table", (Obit*)outData, &ver, access, err);
  if (err->error) Obit_traceback_val (err, routine, outData->name, targ);

  /* Open table */
  if ((ObitTableOTFTargetOpen (table, access, err) 
       != OBIT_IO_OK) || (err->error))  { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR opening output OTFTarget table");
    return targ;
  }

  /* Create Row */
  row = newObitTableOTFTargetRow (table);

  /* attach to table buffer */
  ObitTableOTFTargetSetRow (table, row, err);
  if (err->error) Obit_traceback_val (err, routine, outData->name, targ);

  /* Newly created?  Just write new one */
  doWrite = FALSE;
  if (isNew) {
    targ = 1;
    row->TargID = targ;
    strncpy(row->Target, name, 16);
    doWrite = TRUE;
  } else { /* Existing, see if already exists? */

    /* loop through table */
    for (iRow = 1; iRow<=table->myDesc->nrow; iRow++) {
      if ((ObitTableOTFTargetReadRow (table, iRow, row, err)
	   != OBIT_IO_OK) || (err->error>0)) { 
	Obit_log_error(err, OBIT_Error, "ERROR reading OTFTarget Table file");
	return targ;
      }
      if (!strncmp (row->Target, name, 16)) {
	/* Found match */
	targ = row->TargID;
	break;
      }  
    } /* end loop over table */

    /* Add new entry? */
    if (targ<=0) {
      targ = table->myDesc->nrow + 1;
      row->TargID = targ;
      strncpy(row->Target, name, 16);
      doWrite = TRUE;
    }
  } /* end output table already exists */

  /* need to write new entry? */
  if (doWrite) {
    iRow = table->myDesc->nrow + 1;
    if ((ObitTableOTFTargetWriteRow (table, iRow, row, err)
	 != OBIT_IO_OK) || (err->error>0)) { 
      Obit_log_error(err, OBIT_Error, "ERROR writing OTFTarget Table file");
      return targ;
    }
  }
  
 /* Close  table */
  if ((ObitTableOTFTargetClose (table, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR closing output OTFTarget Table file");
    return targ;
  }

  /* Cleanup */
  row = ObitTableOTFTargetRowUnref(row);
  table = ObitTableOTFTargetUnref(table);

  return targ;
} /* end  GetTarget */

void InitScan (ObitOTF *outData, gboolean isNew, ofloat scan, ofloat target,
	       ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Initializes Index table for this scan creating an entry               */
/*   Input:                                                               */
/*      outData  Output OTF object                                        */
/*      isNew    True if output file just created                         */
/*      scan     Scan number                                              */
/*      target   Target ID number                                         */
/*   Output:                                                              */
/*       err       Obit return error stack                                */
/*----------------------------------------------------------------------- */
{
  ObitTableOTFIndex* table;
  ObitTableOTFIndexRow* row;
  olong iRow, ver;
  olong scanID, targetID;
  ObitIOAccess access;
  gchar *routine = "InitScan";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitOTFIsA(outData));

  /* create Index table object */
  scanID = (olong)(scan+0.5);
  targetID = (olong)(target+0.5);
  ver = 1;
  if (isNew) access = OBIT_IO_WriteOnly;
  else access = OBIT_IO_ReadWrite;
  table = newObitTableOTFIndexValue ("Index table", (Obit*)outData, &ver, access, err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);

  /* Open table */
  if ((ObitTableOTFIndexOpen (table, access, err) 
       != OBIT_IO_OK) || (err->error))  { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR opening output OTFIndex table");
    return;
  }

  /* Create Row */
  row = newObitTableOTFIndexRow (table);

  /* initialize row */
  row->ScanID = scanID;
  row->TargetID = targetID;
  row->Time = 0.0;
  row->TimeI = 0.0;
  row->StartRec = -1;
  row->EndRec = -1;

  /* attach to table buffer */
  ObitTableOTFIndexSetRow (table, row, err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);

  /* Write at end of table */
  iRow = table->myDesc->nrow + 1;
  if (isNew) iRow = 1;
  if ((ObitTableOTFIndexWriteRow (table, iRow, row, err)
       != OBIT_IO_OK) || (err->error>0)) { 
    Obit_log_error(err, OBIT_Error, "ERROR writing OTFIndex Table file");
    return;
  }

 /* Close  table */
  if ((ObitTableOTFIndexClose (table, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR closing output OTFIndex Table file");
    return;
  }

  /* Cleanup */
  row = ObitTableOTFIndexRowUnref(row);
  table = ObitTableOTFIndexUnref(table);

} /* end  InitScan */

void SetScan (ObitOTF *outData, odouble startTime, odouble endTime, 
	      olong startRec, olong endRec, ObitErr *err)
/*----------------------------------------------------------------------- */
/*  Initializes Index table for this scan creating an entry               */
/*   Input:                                                               */
/*      outData   Output OTF object                                       */
/*      isNew     True if output file just created                        */
/*      startTime Start time of scan in days                              */
/*      endTime   End time of scan in days                                */
/*      startRec  First record in scan                                    */
/*      endRec    Last record in scan                                     */
/*   Output:                                                              */
/*       err       Obit return error stack                                */
/*----------------------------------------------------------------------- */
{
  ObitTableOTFIndex* table;
  ObitTableOTFIndexRow* row;
  olong iRow, ver;
  ObitIOAccess access;
  gchar *routine = "SetScan";

  /* error checks */
  g_assert(ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitOTFIsA(outData));

  /* create Index table object */
  ver = 1;
  access = OBIT_IO_ReadWrite;
  table = newObitTableOTFIndexValue ("Index table", (Obit*)outData, &ver, access, err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);

  /* Open table */
  if ((ObitTableOTFIndexOpen (table, access, err) 
       != OBIT_IO_OK) || (err->error))  { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR opening output OTFIndex table");
    return;
  }

  /* Create Row */
  row = newObitTableOTFIndexRow (table);

  /* attach to table buffer */
  ObitTableOTFIndexSetRow (table, row, err);
  if (err->error) Obit_traceback_msg (err, routine, outData->name);

  /* Update last record */
  iRow = table->myDesc->nrow;
  if ((ObitTableOTFIndexReadRow (table, iRow, row, err)
       != OBIT_IO_OK) || (err->error>0)) { 
    Obit_log_error(err, OBIT_Error, "ERROR reading OTFIndex Table file");
    return;
  }

  /* upate row */
  row->Time = 0.5 * (startTime + endTime);
  row->TimeI = (endTime - startTime);
  row->StartRec = startRec;
  row->EndRec   = endRec;

  /* Rewrite at end of table */
  iRow = table->myDesc->nrow;
  if ((ObitTableOTFIndexWriteRow (table, iRow, row, err)
       != OBIT_IO_OK) || (err->error>0)) { 
    Obit_log_error(err, OBIT_Error, "ERROR writing OTFIndex Table file");
    return;
  }

 /* Close  table */
  if ((ObitTableOTFIndexClose (table, err) 
       != OBIT_IO_OK) || (err->error>0)) { /* error test */
    Obit_log_error(err, OBIT_Error, "ERROR closing output OTFIndex Table file");
    return;
  }

  /* Cleanup */
  row = ObitTableOTFIndexRowUnref(row);
  table = ObitTableOTFIndexUnref(table);

} /* end  SetScan */

