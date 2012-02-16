/* $Id$ */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2012                                               */
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
/*; Correspondence about this software should be addressed as follows:*/
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/

#include "ObitTableVLUtil.h"
#include "ObitTableUtil.h"
#include "ObitSkyGeom.h"
#include "ObitPosLabelUtil.h"
#include "ObitPrinter.h"
#include "ObitImageFitData.h"
/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitSurveyUtil.c
 * ObitTableVL class utility function definitions.
 */

/*----------------------Private prototypes---------------------------*/
static void NVSSCorErr (odouble *ra, odouble *dec, 
			ofloat *peak, ofloat *major, ofloat *minor, ofloat *posang,
			ofloat qcent, ofloat ucent, ofloat *pflux, 
			ofloat irms, ofloat prms, ofloat *beam, 
			gboolean fitted, gboolean doraw, ofloat fblank,
			ofloat *flux, ofloat *eflux, ofloat *epflux, 
			gchar chpang[7], gchar chepan[7], ofloat * errra, ofloat *errdec, 
			gchar cmajor[7], gchar cminor[7], gchar cpa[7], 
			gchar emajor[7], gchar eminor[7], gchar epa[7]);
static void VLSSCorErr (odouble *ra, odouble *dec, ofloat *peak, 
			ofloat *major, ofloat *minor, ofloat *posang,
			ofloat irms, ofloat *beam, 
			gboolean fitted, gboolean doraw, ofloat fblank,
			ofloat *flux, ofloat *eflux, 
			ofloat * errra, ofloat *errdec, 
			gchar cmajor[7], gchar cminor[7], gchar cpa[7], 
    		        gchar emajor[7], gchar eminor[7], gchar epa[7]);
static void NVSSpdbias (ofloat *p, ofloat rms);
static void NVSSdeconv (ofloat ptsmaj, ofloat ptsmin, ofloat ptspa, 
			ofloat fitmaj, ofloat ufitmj, ofloat fitmin, ofloat ufitmn, 
			ofloat fitpa, 
			ofloat* srcmaj, ofloat* usrcmj, ofloat* srcmin, ofloat* usrcmn, 
			ofloat* srcpa);
static void NVSSnewPoint (ofloat decrad, ofloat spmjy, 
			  ofloat *psmaj, ofloat *psmin, ofloat *pspa);
static void NVSSbmval (ofloat bmaj, ofloat bmin, ofloat bpa, 
		       ofloat bmaje, ofloat bmine, ofloat bpae, 
		       ofloat cbmaj, ofloat cbmin, ofloat cbpa, 
		       gchar* cmajor, gchar* cminor, gchar* cpa, 
		       gchar* emajor, gchar* eminor, gchar* epa, 
		       gboolean *resolv, gboolean hafres[3], 
		       olong *ier);
gboolean NVSSPrint (ObitPrinter *printer, ObitData *data, olong VLVer, gboolean first, 
		    gboolean last, ObitErr* err);
static olong VLFindRA (odouble ra, olong* VLindex, ObitTableVL *VLTable, 
			 ObitTableVLRow *VLRow, ObitErr* err);
static void GetVLIndex(ObitTableVL* VLTable, olong *numind, olong VLindex[25]);

/*----------------------Public functions---------------------------*/

/**
 * Print raw contents of VL table
 * \param in        VL table to print
 * \param image     Image to which VL table is attached.
 * \param prtFile   Where to write
 * \param err       ObitErr error stack.
 */
void ObitSurveyUtilVLPrint (ObitTableVL *in, ObitImage *image, FILE  *prtFile, 
			    ObitErr *err)
{
  ObitTableVLRow *row = NULL;
  olong i, irow;
  ofloat maj, min, pa;
  gchar rast[19], decst[19], field[9];
  ofloat epeak, errra, errdec, errmaj, errmin, errpa;
  ofloat beam[3], beamas[2], xcell, flux, eflux;
  gchar *routine = "ObitSurveyUtilVLPrint";

  /* error checks */
  if (err->error) return;
  Obit_return_if_fail (ObitTableVLIsA(in), err, 
		       "%s input %s not an VL Table", routine, in->name);
  Obit_return_if_fail (ObitImageIsA(image), err, 
		       "%s image %s not an image", routine, image->name);

  ObitTableVLOpen (in, OBIT_IO_ReadOnly, err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);

  /* Image stuff */
  xcell = fabs (image->myDesc->cdelt[0]);
  beam[0] = image->myDesc->beamMaj / xcell;  /* cells */
  beam[1] = image->myDesc->beamMin / xcell;
  beam[2] = image->myDesc->beamPA * DG2RAD;
  beamas[0] = image->myDesc->beamMaj * 3600.0;  /* asec */
  beamas[1] = image->myDesc->beamMin * 3600.0;

  fprintf (prtFile,"\n Listing of fitted VL table values\n");
  fprintf (prtFile,"Fitted sizes in asec, Peak, Flux, IRMS in mJy, residual values relative to Peak\n");
  fprintf (prtFile,"Error estimates (asec, mJy, deg) given under value\n");
  fprintf (prtFile,
	   "             RA           Dec          Peak    Flux    IRMS  Fit Maj Fit min   PA    res. RMS res Peak  PixX    PixY   Field\n");

  row  = newObitTableVLRow (in);

  /* Loop over table printing */
  for (irow=1; irow<=in->myDesc->nrow; irow++) {
    ObitTableVLReadRow (in, irow, row, err);
   if (err->error) Obit_traceback_msg (err, routine, in->name);
   /*if (VLdesel(row)) continue;  *//* Skip deselected record */

    ObitPosLabelUtilRA2HMS (row->Ra2000, image->myDesc->ctype[0], rast);
    ObitPosLabelUtilDec2DMS (row->Dec2000, image->myDesc->ctype[1], decst);
    maj = row->MajorAxis * 3600.0;
    min = row->MinorAxis * 3600.0;
    pa  = row->PosAngle;
    for (i=0; i<8; i++) field[i] = row->Field[i]; field[i] = 0;

    /* Errors */
    ObitImageFitDataGaussErr (row->PeakInt, 
			      row->MajorAxis/xcell, row->MinorAxis/xcell, row->PosAngle*DG2RAD,
			      row->IRMS, (ofloat*)&beam[0],
			      &epeak, &errra, &errdec, &errmaj, &errmin, &errpa);

    /* Flux */
    flux = row->PeakInt*((maj/beamas[0]) * (min/beamas[1]));
    eflux = epeak * ((maj/beamas[0]) * (min/beamas[1]));

    /* Values */
    fprintf (prtFile," %5d %14s %14s %8.2f %8.2f %7.3f %7.3f %7.3f %6.1f %8.3f %8.3f %7.1f %7.1f %8s\n",
	     irow, rast, decst, row->PeakInt*1000.0, flux*1000.0, row->IRMS*1000.0, maj, min, pa,
	     row->ResRMS/row->PeakInt, row->ResPeak/row->PeakInt,
	     row->CenterX, row->CenterY, field);
    /* errors */
    fprintf (prtFile,"              %6.2f         %6.2f %8.2f %8.2f         %7.3f %7.3f %6.1f \n",
	     errra*xcell*3600.0, errdec*xcell*3600.0, epeak*1000.0, eflux*1000.0,
	     errmaj*xcell*3600.0, errmin*xcell*3600.0, errpa*RAD2DG);
  } /* end loop over table */
  
  /* Close up */
  ObitTableVLClose (in,  err);
  if (err->error) Obit_traceback_msg (err, routine, in->name);
 
  /* release row object */
  row = ObitTableVLRowUnref(row);
} /* end ObitSurveyUtilVLPrint */

/**
 * Print selected NVSS entries from VL table
 * If Search<=0 and Box[0]<=0 and  Box[1]<=0, all entries are printed
 * Routine adapted from the AIPSish NVSSlist.f/PRTVLT  
 * \param printer    Printer object, created and closed externally
 * \param data       Data file to which VL table is attached.
 *                   info member with:
 * \li Object    OBIT_string  Name of object
 * \li equinCode OBIT_long    Epoch code for output, 
 *                            1=>B1950, 2=>J2000, 3=>Galactic [def 2]
 * \li Fitted    OBIT_bool    If TRUE give fitted values [def F]
 * \li doraw     OBIT_bool    If TRUE give raw values from table, else 
 *                            corrected/deconvolved.[def F]
 * \li RA        OBIT_double  RA center of search, degrees in equinCode [def 0]
 * \li Dec       OBIT_double  Dec center of search, degrees in equinCode [def 0]
 * \li Search    OBIT_double  Search radius in arcsec, <= 0 => all selected. [def 15] 
 * \li Box       OBIT_double[2] RA and Dec halfwidth of search 
 *                              rectangle in hr,deg [0,0]
 * \li Silent    OBIT_double  Half width asec of silent search box. [def 720]
 * \li minFlux   OBIT_float   Minimum peak flux density. Jy [def 0]
 * \li maxFlux   OBIT_float   Maximum peak flux density. Jy [def LARGE]
 * \li minPol    OBIT_float   Minimum percent integrated polarization [def 0]
 * \li minGlat   OBIT_float   Minimum abs galactic latitude [def any]
 * \li maxGlat   OBIT_float   Minimum abs galactic latitude [def any]
 * \param VLVer      VL table version number
 * \param first      First call? open printer 
 * \param last       Last call? close printer when done. 
 * \param err        Obit Error/message stack
 * \return   True if user requested Quit.
 */
gboolean ObitSurveyNVSSPrint (ObitPrinter *printer, ObitData *data, olong VLVer, 
			      gboolean first, gboolean last, ObitErr* err)
{
  olong   irow, bc, ec, inc, rahm[2], decdm[2], sort=1, 
    numind, VLindex[25], bcindx, ecindx, irab, irae, ipass, npass, 
    beg[2], end[2], bci, eci, imark,  numcol, itemp, ibadvl, VLnrow, 
     width=132;
  ofloat      ras, decs, beam[3], tcut, flux, pa, eflux, epflux, 
    errra, errdec, pctpol;
  gboolean   wanted, indxed=FALSE, doall, select, norad, nobox,
    found, dosil, dogal, quit=FALSE;
  gchar  line[133], eline[133],  dsig[2];
  gchar cmajor[7], cminor[7], cpa[7], emajor[7], eminor[7], epa[7], 
    chpang[7], chepan[7], cdist[7], cposa[7],  
    chpflx[7], chepfx[7], chbdvl[5],cflux[8], ceflux[8];
  odouble rac, decc, decx, ra0, dec0, rab, rae, radius, radr, dist, 
    rabeg, raend, decbeg, decend, discl, rat, dect, dist2, radr2, 
    boxra, boxdec, odiscl, mind2, mindm, sildeg, glat, glon;
  ofloat   peak, major, minor, posang, qcent, ucent, pflux, l, m;
  ofloat   irms, prms, resrms, respek, resflx, cenx, ceny;
  ofloat   fblank = ObitMagicF();
  gchar Object[201];
  olong equinCode;
  gboolean fitted, doraw;
  odouble ra, dec, search, box[2], silent;
  ofloat minFlux, maxFlux, minPol, minGlat, maxGlat, farr[2];
  ObitTableVL *VLTable=NULL;
  ObitTableVLRow *VLRow=NULL;
  ObitInfoType type;
  union ObitInfoListEquiv InfoReal; 
  gint32       dim[MAXINFOELEMDIM];
  gchar field[9], datevl[21];
  gchar *mark[] = {"   ", " r*", " p*", " s*"};
  ofloat cutt = 0.002;
  /* Things to save between calls */
  static gchar  Title1[133], Title2[133];
  static olong page, pageno, lpage, scount, ecount, ncount;
  gchar *routine = "ObitSurveyNVSSPrint";

  /* Error checks */
  g_assert(ObitPrinterIsA(printer));
  if (err->error) return FALSE;  /* previous error? */

  /* Get parameters - be somewhat tolerant of the types */
  snprintf (Object, 201, "    ");
  if (ObitInfoListGetTest(data->info, "Object",    &type, dim, Object))
    Object[dim[0]] = 0;  /* Make sure null terminated */
  equinCode = 2;
  if (ObitInfoListGetTest(data->info, "equinCode", &type, dim, &InfoReal)) {
    if (type==OBIT_oint) equinCode = InfoReal.itg;
    else equinCode = InfoReal.otg;   
  }
  fitted = FALSE;
  ObitInfoListGetTest(data->info, "Fitted",    &type, dim, &fitted);
  doraw = FALSE;
  ObitInfoListGetTest(data->info, "doraw",     &type, dim, &doraw);
  ra = 0.0;
  if (ObitInfoListGetTest(data->info, "RA",    &type, dim, &InfoReal)) {
    if (type==OBIT_float)       ra = InfoReal.flt;
    else if (type==OBIT_double) ra = InfoReal.dbl; }  
  dec = 0.0;
  if (ObitInfoListGetTest(data->info, "Dec",       &type, dim, &InfoReal)) {
    if (type==OBIT_float)      dec = InfoReal.flt;
    else if (type==OBIT_double) dec = InfoReal.dbl; }  
  search = 15.0;
  if (ObitInfoListGetTest(data->info, "Search",    &type, dim, &InfoReal)) {
    if (type==OBIT_float)       search = InfoReal.flt;
    else if (type==OBIT_double) search = InfoReal.dbl; }  
  box[0] = 0.0; box[1]=0.0;
  if (ObitInfoListGetTest(data->info, "Box",       &type, dim, &box)) {
    if (type==OBIT_float) {     
      ObitInfoListGetTest(data->info, "Box",       &type, dim, &farr);
      box[0] = (odouble)farr[0];
      box[1] = (odouble)farr[1]; } }
  silent = 720.;
  if (ObitInfoListGetTest(data->info, "Silent",    &type, dim, &InfoReal)) {
    if (type==OBIT_float)      silent = InfoReal.flt;
    else if (type==OBIT_double) silent = InfoReal.dbl; }  
  minFlux = 0.0;
  if (ObitInfoListGetTest(data->info, "minFlux",   &type, dim, &InfoReal)) {
     if (type==OBIT_float)      minFlux = InfoReal.flt;
    else if (type==OBIT_double) minFlux = InfoReal.dbl; }  
  maxFlux = 100000.;
  if (ObitInfoListGetTest(data->info, "maxFlux",   &type, dim, &InfoReal)) {
    if (type==OBIT_float)       maxFlux = InfoReal.flt;
    else if (type==OBIT_double) maxFlux = InfoReal.dbl; }  
  minPol = 0.0;
  if (ObitInfoListGetTest(data->info, "minPol",    &type, dim, &InfoReal)) {
    if (type==OBIT_float)      minPol = InfoReal.flt;
    else if (type==OBIT_double) minPol = InfoReal.dbl; }  
  minGlat = 0.0;
  if (ObitInfoListGetTest(data->info, "minGlat",   &type, dim, &InfoReal)) {
    if (type==OBIT_float)      minGlat = InfoReal.flt;
    else if (type==OBIT_double) minGlat = InfoReal.dbl; }  
  maxGlat = 90.;
  if (ObitInfoListGetTest(data->info, "maxGlat",   &type, dim, &InfoReal)) {
    if (type==OBIT_float)       maxGlat = InfoReal.flt;
    else if (type==OBIT_double) maxGlat = InfoReal.dbl; }  
  
  /* initialize */
  snprintf (datevl, 21, "unknown");
  snprintf (eline, 132, "eline not initialized ");
  found = FALSE;
  numcol = width; /* Output line limit */

  if (first) {
    /* Reset page titles */
    snprintf (Title1, 132, "      ");
    snprintf (Title2, 132, "      ");
    scount = 0;
    ecount = 0;
    ncount = 0;
    page   = 1;
    pageno = 1;
    odiscl = 0.0e0;

    /* Print header info */
    snprintf (line, 132, "      ");
    ObitPrinterWrite (printer, line, &quit, err);
    if (quit) goto Quit;
    if (err->error) goto cleanup;
    snprintf (line,  132, "NRAO/VLA Sky Survey (NVSS) catalog search, ver 3.0");
    ObitPrinterWrite (printer, line, &quit, err);
    if (quit) goto Quit;
    if (err->error) goto cleanup;
    
    snprintf (line , 132, "Error estimates appear below the value.");
    ObitPrinterWrite (printer, line, &quit, err);
    if (quit) goto Quit;
    if (err->error) goto cleanup;
    
    snprintf (line, 132, "Using VL table %d ", VLVer);
    ObitPrinterWrite (printer, line, &quit, err);
    if (quit) goto Quit;
    if (err->error) goto cleanup;
  } /* end if first */
    
  /* Open input VL table */
  VLTable = newObitTableVLValue (data->name, data, &VLVer, OBIT_IO_ReadOnly, err);
  if (err->error) goto cleanup;
  /* Should be there */
  Obit_retval_if_fail((VLTable!=NULL), err, quit, "VL table %d does not exist", VLVer);
  
  ObitTableVLOpen (VLTable, OBIT_IO_ReadOnly, err);
  if (err->error) goto cleanup;
  /* Create row structure */
  VLRow = newObitTableVLRow(VLTable);
  VLnrow = VLTable->myDesc->nrow;   /* How many rows */

  /* Get information from header */
  beam[0] = VLTable->BeamMajor;
  beam[1] = VLTable->BeamMinor;
  beam[2] = VLTable->BeamPA;
  numind  = VLTable->numIndexed;
  strncpy (datevl, ((ObitImage*)data)->myDesc->date, 10); datevl[10] = 0;
  GetVLIndex (VLTable, &numind, VLindex);
  
  if (beam[0] <= 0.0) beam[0] = 45.0 / 3600.0;
  if (beam[1] <= 0.0) beam[1] = 45.0 / 3600.0;

  /* Galactic limits?*/
  dogal = (minGlat > 0.0)  ||  (maxGlat < 90.0);

  /* Rows to search */
  bci = 1;
  eci = VLnrow;
  eci = MAX (bci, eci);
  inc = 1;
  if (inc <= 0) inc = 1;

  /* Silent search? */
  dosil = silent > 0.0e0;
  silent = MAX (0.0, silent);
  mind2 = 1.0e20;
  sildeg = silent / 3600.0;
  
  /* Position search box */
  rac  = ra;
  decc = dec;
  ra0  = rac * DG2RAD;
  dec0 = decc * DG2RAD;
  radius = search / 3600.0e0;

  /* Search window? */
  doall = (radius <= 0.0)  &&  (box[0] <= 0.0e0)  &&  (box[1] <= 0.0e0);
  /* No radius specified? */
  norad = radius <= 0.0;
  if (norad) radius = MAX (box[0], box[1]);
  radr  = radius;
  radr2 = radr * radr;

  /* No Box? */
  nobox = ((box[0] <= 0.0)  &&  (box[1] <= 0.0));
  if (box[0] <= 0.0) box[0] = MAX (MAX(box[0], box[1]), radius);
  if (box[1] <= 0.0) box[1] = MAX (MAX(box[0], box[1]), radius);
  if (!nobox) {
    /*???         RADIUS = MAX (BOX(1), BOX(2)) */
    radr2 = 1.0e20;
    norad = TRUE;
  } 

  /* RA box fullwidth in hours */
  boxra  = 2.0e0 * box[0];
  boxdec = 2.0e0 * box[1];

  /* Always give distance in sec and  PA */
  discl = 3600.0e0;
  /* All positions */
  if (doall) {
    decbeg = -100;
    decend = 100;
    bcindx = 1;
    ecindx = VLnrow;
    npass = 1;
    beg[0] = bcindx;
    end[0] = ecindx;
    rab = 0.0;
    rae = 360.0;
  } else {
    /* Select position range */
    decbeg = decc - MAX (box[1], sildeg);
    decend = decc + MAX (box[1], sildeg);
    decx = MIN (89.0e0, decc);
    rab = rac - MAX (MAX(radius, box[0]), sildeg) / cos (decx * DG2RAD);
    if (rab < 0.0) rab = rab + 360.0;
    rae = rac + MAX (MAX(radius, box[0]), sildeg) / cos (decx * DG2RAD);
    if (rae > 360.0) rae = rae - 360.0;
    irab = rab / 15.0;
    irae = rae / 15.0;
    irab = MIN (23, MAX (0, irab));
    irae = MIN (24, MAX (0, irae));
    VLindex[24] = VLnrow;
    /* Table indexed? */
    indxed = VLindex[0] > 0;
    if (indxed) {
      bcindx = VLindex[irab];
      ecindx = VLindex[MIN(23, irae)+2];
    } else {
      bcindx = 1;
      ecindx = VLnrow;
    } 

    /* It takes two passes for wrap in  RA range. */
    if (irab > irae) {
      npass = 2;
      beg[0] = bcindx;
      end[0] = VLnrow;
      beg[1] = 1;
      end[1] = ecindx;
    } else {
      npass = 1;
      beg[0] = bcindx;
      end[0] = ecindx;
    } 
  }
 
  /* Tell selection criteria */
  if (first) {
    /* Minimum flux density */
    if (minFlux > 1.0e-5) {
      snprintf (line ,132,  "Selecting sources brighter than %9.1f mJy", 
	       minFlux*1000.0);
      ObitPrinterWrite (printer, line, &quit, err);
      if (quit) goto Quit;
      if (err->error) goto cleanup;
    }

    /* Maximum flux density */
    if (maxFlux < 1.0e5) {
      snprintf (line ,132,  "Selecting sources fainter than %9.1f mJy", 
	       maxFlux*1000.0);
      ObitPrinterWrite (printer, line, &quit, err);
      if (quit) goto Quit;
      if (err->error) goto cleanup;
    } 

    /* Minimum percentage pol. */
    if (minPol > 1.0e-5) {
      snprintf (line ,132,  "Selecting sources more than %5.1f percent polarized", 
	       minPol);
      ObitPrinterWrite (printer, line, &quit, err);
      if (quit) goto Quit;
      if (err->error) goto cleanup;
    } 

    /* Galactic latitude range */
    snprintf (line ,132,  "abs. galactic latitude range %5.1f to %5.1f",
	     minGlat, maxGlat);
    ObitPrinterWrite (printer, line, &quit, err);
    if (quit) goto Quit;
    if (err->error) goto cleanup;
  
    /* Table data and number of entries */
    snprintf (line ,132,  "catalog made on %s (dd/mm/yy) with %7i entries",
	     datevl, VLnrow);
    ObitPrinterWrite (printer, line, &quit, err);
    if (quit) goto Quit;
    if (err->error) goto cleanup;

    /* Indexing */
    if (!doall) {
      if (indxed) {
	snprintf (line, 132, "Table indexed for faster access");
	ObitPrinterWrite (printer, line, &quit, err);
	if (quit) goto Quit;
	if (err->error) goto cleanup;
      } else {
 	snprintf (line, 132, "NOTE:table not indexed");
	ObitPrinterWrite (printer, line, &quit, err);
	if (quit) goto Quit;
	if (err->error) goto cleanup;
      } 
    }
    
    /* Show fitted/deconvolved sizes? */
    if (fitted) {
      if (doraw) snprintf (line, 132, "Displaying raw component size, peak flux density");
      else snprintf (line, 132, "Displaying fitted component size, peak flux density");
    } else {
      snprintf (line, 132, "Displaying deconvolved component size, integrated flux density");
    } 
    ObitPrinterWrite (printer, line, &quit, err);
    if (quit) goto Quit;
    if (err->error) goto cleanup;

    /* Res codes */
    snprintf (line, 132, "residual (res) code; nonblank indicates complex source structure:");
    ObitPrinterWrite (printer, line, &quit, err);
    if (quit) goto Quit;
    if (err->error) goto cleanup;

    snprintf (line, 132, "   p* => high peak, r* => high rms, s* => high integral");
    ObitPrinterWrite (printer, line, &quit, err);
    if (quit) goto Quit;
    if (err->error) goto cleanup;
  
     snprintf (line, 132, "      ");
     ObitPrinterWrite (printer, line, &quit, err);
     if (quit) goto Quit;
     if (err->error) goto cleanup;
  } /* end if first */

  if (!doall) {
    /* Position at input equinox */
    ObitPosLabelUtilRAHMS  (rac,  &rahm[0],  &rahm[1],  &ras);
    ObitPosLabelUtilDecDMS (decc, &decdm[0], &decdm[1], &decs);
    /* Deal with sign of -0 */
    if (decc > 0.0) dsig[0] = '+';
    if (decc < 0.0) dsig[0] = '-';
    dsig[1] = 0;
    decdm[0] = abs(decdm[0]);

    /* Need to change target equinx to J2000? */
    if (equinCode == 1) {
      ObitSkyGeomBtoJ (&rac, &decc);  /* from B1950 to J2000 */
    } else if (equinCode==3) {
      /* Given galactic which are defined in B1950 */
      ObitSkyGeomGal2Eq (&rac, &decc);  /* from Galactic to B1950 */
      ObitSkyGeomBtoJ   (&rac, &decc);  /* from B1950 to J2000 */
    } 
    
    snprintf (line, 132, "      ");
    ObitPrinterWrite (printer, line, &quit, err);
    if (quit) goto Quit;
    if (err->error) goto cleanup;
    
    /* Object label */
    if (!norad) {
      snprintf (line, 132, "%s: search within %8.1f arcsec of  %2.2d %2.2d %7.3f %s%2.2d %2.2d %6.2f",
		Object, radius*3600.0, rahm[0], rahm[1], ras, dsig, decdm[0], decdm[1], decs);
      ObitPrinterWrite (printer, line, &quit, err);
      if (quit) goto Quit;
      if (err->error) goto cleanup;
    }
    
    if (!nobox) {
      snprintf (line, 132, "selection box: %5.2f hr x %6.2f deg with center %2.2d %2.2d %7.3f %s%2.2d %2.2d %6.2f",
		boxra, boxdec, rahm[0], rahm[1], ras, dsig, decdm[0], decdm[1], decs);
      ObitPrinterWrite (printer, line, &quit, err);
      if (quit) goto Quit;
      if (err->error) goto cleanup;
    } 
    
    snprintf (line, 132, "      ");
    ObitPrinterWrite (printer, line, &quit, err);
    if (quit) goto Quit;
    if (err->error) goto cleanup;
    
    /* Check declination */
    if (dec < -40.0) {
      snprintf (line, 132, "WARNING: the NVSS is complete only to declination -40!");
      ObitPrinterWrite (printer, line, &quit, err);
      if (quit) goto Quit;
      if (err->error) goto cleanup;
      snprintf (line, 132, "      ");
      ObitPrinterWrite (printer, line, &quit, err);
      if (quit) goto Quit;
      if (err->error) goto cleanup;
    } 
  } /* end not doall */
  
  /* Page labels */
  if (!doall) {
    /*if (numcol < 93) numcol = 70;*/
    /* Default with pos select */
    if (fitted) {
      snprintf (Title1, 132, "   RA[2000]  Dec[2000] Dist(\")  Peak  Major Minor    PA   Res P_Flux P_ang   Field    X_pix  Y_pix");
      snprintf (Title2, 132, " h  m    s    d  m   s   ori      mJy   \"     \"      deg       mJy    deg");
    } else{
      snprintf (Title1, 132, "   RA[2000]  Dec[2000] Dist(\")  Flux  Major Minor    PA   Res P_Flux P_ang   Field    X_pix  Y_pix");
      snprintf (Title2, 132, " h  m    s    d  m   s   ori      mJy   \"     \"     deg        mJy   deg");
    }
  } else {
    /*if (numcol < 87) numcol = 64;*/
    /* Default doall listing */
    if (fitted) {
      snprintf (Title1, 132, "   RA[2000]   Dec[2000]   Peak  Major Minor    PA    Res P_Flux  P_ang  Field    X_pix    Y_pix");
      snprintf (Title2, 132, " h  m    s    d  m   s     mJy    \"     \"      deg         mJy    deg");
    } else {
      snprintf (Title1, 132, "   RA[2000]   Dec[2000]   Flux  Major Minor    PA   Res P_Flux P_ang   Field    X_pix    Y_pix");
      snprintf (Title2, 132, " h  m    s    d  m   s     mJy    \"     \"     deg        mJy     deg");
    }
  } /* end doall */
  /* Set epoch */
  if (equinCode == 1) {
    if (fitted) 
      snprintf (Title1, 132, "   RA[1950]   Dec[1950]   Peak  Major Minor   PA Res P_Flux P_ang  Field    X_pix    Y_pix");
    else
      snprintf (Title1, 132, "   RA[1950]   Dec[1950]   Flux  Major Minor   PA Res P_Flux P_ang  Field    X_pix    Y_pix");
  }  else if (equinCode == 3) {
    if (fitted) 
      snprintf (Title1, 132, "   Gal Long   Gal Lat     Peak  Major Minor   PA Res P_Flux P_ang  Field    X_pix    Y_pix");
    else
      snprintf (Title1, 132, "   Gal Long   Gal Lat     Flux  Major Minor   PA Res P_Flux P_ang  Field    X_pix    Y_pix");
  } 
  
  /* Give headers */
  ObitPrinterSetTitle (printer, Title1, Title2, err);  /* Page titles at top */
  ObitPrinterWrite (printer, Title1, &quit, err);
  ObitPrinterWrite (printer, Title2, &quit, err);
  if (err->error) goto cleanup;
  
  odiscl = discl;
  /* Passes */
  for (ipass= 1; ipass<=npass; ipass++) { /* loop 500 */
    /* Set range of rows. */
    bc = MIN (MAX (bci, beg[ipass-1]), VLnrow);
    if (ec <= 0) ec = VLnrow;
    ec = MIN (eci,  end[ipass-1]);
    if (inc <= 0) inc = 1;
    if (npass == 1) {
      /* No RA wrap in search range */
      rabeg = rab;
      raend = rae;
    }  else if (ipass == 1) {
      /* Wrap - high hours */
      rabeg = rab;
      raend = 360.0e0;
    } else {
      /* Wrap - low hours */
      rabeg = 0.0e0;
      raend = rae;
    } 
    
    /* Find starting location  */
    bc = VLFindRA (rabeg, VLindex, VLTable, VLRow, err);
    if (err->error) goto cleanup;
    
    /* Print selected rows. */
    for (irow=bc; irow<=ec; irow++) { /* loop 100 */
      /* Read table row */
      ObitTableVLReadRow (VLTable, irow, VLRow, err);
      if (err->error) goto cleanup;
      
      /* Want this one? */
      if (doall) {
	dist2 = 0.0;
      } else {
	/* Quick test? */
	if ((VLRow->Ra2000 >= rabeg)  &&  (VLRow->Ra2000 <= raend)) {
	  
	  /* Check Declination */
	  if ((VLRow->Dec2000 < decbeg)  ||  (VLRow->Dec2000 > decend)) continue;
	  
	  /* In RA range, use full test. */
	  ObitSkyGeomShiftXY (rac, decc, 0.0, VLRow->Ra2000, VLRow->Dec2000, &l, &m);
	  dist2 = l*l + m*m;
	} else if ((ipass == npass)  &&  (VLRow->Ra2000 > raend)) {
	  /* Past RA Range, quit if sorted */
	  if (sort == 1) break;
	  dist2 = 1.0e10;
	} else {
	  /* Before RA range. */
	  dist2 = 1.0e10;
	} 
      } 
      
      /* Closest source in silent window */
      mind2 = MIN (mind2, dist2);
      wanted = dist2  <=  radr2;
      /* Rectangular box? */
      if (!nobox) wanted = wanted  &&  ((VLRow->Dec2000 >= decbeg)  &&  (VLRow->Dec2000 <= decend));
   
      if (wanted) {
	/* Set output, Need to change equinox? */
	rat  = VLRow->Ra2000;
	dect = VLRow->Dec2000;
	if (equinCode == 1) {
	  ObitSkyGeomJtoB (&rat, &dect);  /* from J2000 to B1950 */
	} else if (equinCode==3) {
	  /* Convert to galactic which are defined in B1950 */
	  ObitSkyGeomJtoB   (&rat, &dect);  /* from J2000 to B1950 */
	  ObitSkyGeomEq2Gal (&rat, &dect);  /* from B1950 to Galactic */
	} 

	/* VL row to local variables */
	peak   = VLRow->PeakInt;
	major  = VLRow->MajorAxis;
	minor  = VLRow->MinorAxis;
	posang = VLRow->PosAngle;
	qcent  = VLRow->QCenter;
	ucent  = VLRow->UCenter;
	pflux  = VLRow->PFlux;
	irms   = VLRow->IRMS;
	prms   = VLRow->PolRMS;
	resrms = VLRow->ResRMS;
	respek = VLRow->ResPeak;
	resflx = VLRow->ResFlux;
	strncpy (field, VLRow->Field,8); field[8] = 0;
	cenx   = VLRow->CenterX;
	ceny   = VLRow->CenterY;
    
	/* Make corrections, get errors */
	NVSSCorErr (&rat, &dect, &peak, &major, &minor, &posang, qcent, ucent, 
		    &pflux, irms, prms, beam, fitted, doraw, fblank, 
		    &flux, &eflux, &epflux, chpang, chepan, &errra, &errdec, 
		    cmajor, cminor, cpa, emajor, eminor, epa);
	ObitPosLabelUtilRAHMS  (rat, &rahm[0], &rahm[1], &ras);
	ObitPosLabelUtilDecDMS (dect, &decdm[0], &decdm[1], &decs);
	/* Deal with sign of -0 */
	if (dect > 0.0) dsig[0] = '+';
	if (dect < 0.0) dsig[0] = '-';
	dsig[1] = 0;
	decdm[0] = abs(decdm[0]);
	
	/* Distance from location with  possible corrections. */
	dist = sqrt (dist2);
	dist = MIN ((discl * dist), 99999.0);
	
	/* As character string */
	if (dist >= 100.0) {
	  itemp = dist + 0.5;
	  snprintf (cdist, 7, "%6d", itemp);
	}  else if (dist >= 10.0) {
	  snprintf (cdist, 7, "%6.1f", dist);
	} else {
	  snprintf (cdist, 7, "%6.2f", dist);
	} 
	
	/* Position angle from center */
	pa = 57.296 * atan2 (l, m+1.0e-20);
	if (pa < 0) {
	  itemp = pa - 0.5;
	} else {
	  itemp = pa + 0.5;
	} 
	/* As character string */
	snprintf (cposa, 7, "%6d", itemp);
	
	/* Clip flux to prevent overflows */
	/*               IF (FLUX.GT.99.99) FLUX = 99.99 */
	/*               IF (EFLUX.GT.99.99) EFLUX = 99.99 */
	if ((pflux != fblank)  &&  (pflux > 0.9999)) pflux = 0.9999;
	if (epflux > 0.9999) eflux = 0.9999;
	
	/* Convert units for output */
	eflux = eflux * 1000.0;
	epflux = epflux * 1000.0;
	
	/* Percent polarization */
	if (pflux != fblank) {
	  pctpol = 100.0 * pflux / flux;
	  snprintf (chpflx, 7, "%6.2f", pflux*1000.0);
	  snprintf (chepfx, 7, "%6.2f", epflux);
	} else {
	  pctpol = 0.0;
	  snprintf (chpflx, 7, " Blank");
	  snprintf (chepfx, 7, "   NA ");
	} 
	
	/* Convert flux density to string */
	if (flux < 99.99) {
	  snprintf (cflux, 8, "%7.1f", flux * 1000.0);
	}  else if (flux < 999.99) {
	  snprintf (cflux, 8, "%7.0f", flux * 1000.0);
	} else {
	  itemp = flux*1000.0 + 0.5;
	  if (itemp > 9999) itemp = 9999;
	  snprintf (cflux, 8, "%7d", itemp);
	} 
	
	/* Convert flux density error to string */
	if (eflux < 99999.9) {
	  snprintf (ceflux, 8, "%7.1f", eflux);
	} else if (flux < 999999.) {
	  snprintf (ceflux, 8, "%7.0f", eflux);
	} else {
	  itemp = eflux + 0.5;
	  if (itemp > 9999999) itemp = 9999;
	  snprintf (ceflux, 8, "%7d", itemp);
	} 
	
	if (pctpol < 0.0) pctpol = 0.0;
	
	/* Goodness of fit code */
	imark = 1;
	tcut   = sqrt (cutt*cutt + (0.01*peak)*(0.01*peak));
	ibadvl = 0;
	if (resrms > tcut) {
	  imark = 2;
	  ibadvl = 10000.0 * resrms + 0.5;
	  if (ibadvl > 9999) ibadvl = 9999;
	} 
	if (fabs(respek) > tcut) {
	  imark = 3;
	  ibadvl = 10000.0 * respek + 0.5;
	  if (ibadvl > 9999) ibadvl = 9999;
	} 
	if (resflx > tcut) {
	  imark = 4;
	  ibadvl = 10000.0 * resflx + 0.5;
	  if (ibadvl > 9999) ibadvl = 9999;
	} 
	if (imark>1) snprintf (chbdvl, 5, "%4d", ibadvl);
	else         snprintf (chbdvl, 5, "    ");
	
	/* Check selection criteria */
	select = (peak  >=  minFlux)  &&  (peak  <=  maxFlux);
	select = select && ((pflux != fblank) || (pctpol >= MAX(0.0,minPol)));
	
	/* Filter out really bad fits  Peak should be at least 3 times  RMS residual */
	select = select && ((peak/resrms) > 3.0);
	
	/* Galactic latitude range */
	if (dogal) {
	glat = VLRow->Ra2000;
	glon = VLRow->Dec2000;
	  ObitSkyGeomJtoB   (&glat, &glon);  /* from J2000 to B1950 */
	  ObitSkyGeomEq2Gal (&glat, &glon);  /* from B1950 to Galactic */
	  glat = fabs (glat);
	  select = select  &&  ((glat >= minGlat)  &&  (glat <= maxGlat));
	} 
	if (select) {
	  found = TRUE;
	  /* Create appropriate line */
	  if (!doall) {
	    /* Default with pos select */
	    snprintf (line, 132, "%2.2d %2.2d%6.2f %s%2.2d %2.2d%5.1f%s %s %s %s %s%s %s%s   %s%8.2f%8.2f", 
		     rahm[0], rahm[1], ras, dsig, decdm[0], decdm[1], decs, cdist, 
		     cflux, cmajor, cminor, cpa, mark[imark-1], 
		     chpflx, chpang, field, cenx, ceny);
	    /* Errors */
	    snprintf (eline, 132, "     %6.2f       %5.1f%s %s %s %s %s%s%s%s", 
		     errra, errdec, cposa, ceflux, emajor, eminor, epa, chbdvl, chepfx, chepan);
	  } else {
	    /* Default doall */
	    snprintf (line, 132, "%2.2d %2.2d%6.2f %s%2.2d %2.2d%5.1f %s %s %s %s %s %s%s   %s%8.2f%8.2f", 
		     rahm[0], rahm[1], ras, dsig, decdm[0], decdm[1], decs, 
		     cflux, cmajor, cminor, cpa, mark[imark-1], 
		     chpflx, chpang, field, cenx, ceny);
	    /* Errors */
	    snprintf (eline, 132, "     %6.2f       %5.1f %s %s %s %s %s%s%s", 
		     errra, errdec, ceflux, emajor, eminor, epa, chbdvl, chepfx, chepan);
	  } 
	  
	  /* Count number printer */
	  scount++;
	  
	  /* Print line and errors */
	  ObitPrinterWrite (printer, line, &quit, err);
	  if (quit) goto Quit;
	  if (err->error) goto cleanup;
	  if (!doraw) ObitPrinterWrite (printer, eline, &quit, err);
	  if (quit) goto Quit;
	  if (err->error) goto cleanup;
	  
	  /* Don't leave widows and orphans */
	  if (page >= lpage) page = 1000;
	} 
      } /* end if wanted */ 
    } /* end loop  L100 over table */;
  } /* end loop  L500: over passes */;

  /* Stuff at end - was it found? */
  if (!found) {
    /* Silent search window? */
    if (dosil) {
      /* Min source distance in arcmin */
      mindm = 60.0e0 * sqrt (mind2) /  DG2RAD;
      if (mindm < 100.0e0) {
	ecount++;
	snprintf (line, 132, "Source not found  closest source is %8.2f arcmin away", mindm);
      } else {
	ncount++;
	snprintf (line, 132, "Source not found  nothing close");
      } 
    } 
    ObitPrinterWrite (printer, line, &quit, err);
    if (quit) goto Quit;
    if (err->error) goto cleanup;
  } /* end if not found */
  
  /* Finished */
  Quit:
  last = last || quit;  /* Done? */
  
  /* Final summary */
  if (last) {
    /* Number found */
    if (scount > 0) {
      snprintf (line, 132, "Found %9d entries", scount);
    } else {
      snprintf (line, 132, "no sources meeting selection criteria found");
    } 
    if (!quit) ObitPrinterWrite (printer, line, &quit, err);
    if (quit) goto cleanup;
    if (err->error) goto cleanup;
    
    /* How many missed? */
    if (ecount > 0) {
      snprintf (line, 132, "%9d fields has no source", ecount);
      if (!quit) ObitPrinterWrite (printer, line, &quit, err);
      if (quit) goto cleanup;
      if (err->error) goto cleanup;
    } 
    
    /* How many not even close? */
    if (ncount > 0) {
      snprintf (line, 132, "%9d fields have nothing close", ncount);
      if (!quit) ObitPrinterWrite (printer, line, &quit, err);
      if (quit) goto cleanup;
      if (err->error) goto cleanup;
    } 
    
  } /* end summary on last */ 
  
    /* Close VL table */
 cleanup:
  ObitTableVLClose (VLTable, err);
  VLTable = ObitTableVLUnref(VLTable);
  VLRow   = ObitTableVLRowUnref(VLRow);
  if (err->error) Obit_traceback_val (err, routine, data->name, quit);
  
  return quit;
} /* end of routine ObitSurveyNVSSPrint */ 
  
/**
 * Print selected VLSS entries from VL table
 * If Search<=0 and Box[0]<=0 and  Box[1]<=0, all entries are printed
 * Routine adapted from the AIPSish VLSSlist.f/PRTVLT  
 * \param printer    Printer object, created and closed externally
 * \param data       Data file to which VL table is attached.
 *                   info member with:
 * \li Object    OBIT_string  Name of object
 * \li equinCode OBIT_long    Epoch code for output, 
 *                            1=>B1950, 2=>J2000, 3=>Galactic [def 2]
 * \li Fitted    OBIT_bool    If TRUE give fitted values [def F]
 * \li doraw     OBIT_bool    If TRUE give raw values from table, else 
 *                            corrected/deconvolved.[def F]
 * \li RA        OBIT_double  RA center of search, degrees in equinCode [def 0]
 * \li Dec       OBIT_double  Dec center of search, degrees in equinCode [def 0]
 * \li Search    OBIT_double  Search radius in arcsec, <= 0 => all selected. [def 15] 
 * \li Box       OBIT_double[2] RA and Dec halfwidth of search 
 *                              rectangle in hr,deg [0,0]
 * \li Silent    OBIT_double  Half width asec of silent search box. [def 720]
 * \li minFlux   OBIT_float   Minimum peak flux density. [def 0]
 * \li maxFlux   OBIT_float   Maximum peak flux density. [def LARGE]
 * \li minGlat   OBIT_float   Minimum abs galactic latitude [def any]
 * \li maxGlat   OBIT_float   Minimum abs galactic latitude [def any]
 * \param VLVer      VL table version number
 * \param first      First call? open printer 
 * \param last       Last call? close printer when done. 
 * \param err        Obit Error/message stack
 * \return   True if user requested Quit.
 */
gboolean ObitSurveyVLSSPrint (ObitPrinter *printer, ObitData *data, olong VLVer, 
			      gboolean first, gboolean last, ObitErr* err)
{
  olong   irow, bc, ec, inc, rahm[2], decdm[2], sort=1, 
    numind, VLindex[25], bcindx, ecindx, irab, irae, ipass, npass, 
    beg[2], end[2], bci, eci, imark,  numcol, itemp, ibadvl, VLnrow, 
     width=132;
  ofloat      ras, decs, beam[3], tcut, flux, pa, eflux, errra, errdec;
  gboolean   wanted, indxed=FALSE, doall, select, norad, nobox,
    found, dosil, dogal, quit=FALSE;
  gchar  line[133], eline[133],  dsig[2];
  gchar cmajor[7], cminor[7], cpa[7], emajor[7], eminor[7], epa[7], 
    cdist[7], cposa[7], chbdvl[5],cflux[8], ceflux[8];
  odouble rac, decc, decx, ra0, dec0, rab, rae, radius, radr, dist, 
    rabeg, raend, decbeg, decend, discl, rat, dect, dist2, radr2, 
    boxra, boxdec, odiscl, mind2, mindm, sildeg, glat, glon;
  ofloat   peak, major, minor, posang, l, m, farr[2];
  ofloat   irms, resrms, respek, resflx, cenx, ceny;
  ofloat   fblank = ObitMagicF();
  gchar Object[201];
  olong equinCode;
  gboolean fitted, doraw;
  odouble ra, dec, search, box[2], silent; 
  ofloat minFlux, maxFlux, minGlat, maxGlat;
  ObitTableVL *VLTable=NULL;
  ObitTableVLRow *VLRow=NULL;
  ObitInfoType type;
  gint32       dim[MAXINFOELEMDIM];
  union ObitInfoListEquiv InfoReal; 
  gchar field[9], datevl[21];
  gchar *mark[] = {"   ", " r*", " p*", " s*"};
  ofloat cutt = 0.002;
  /* Things to save between calls */
  static gchar  Title1[133], Title2[133];
  static olong page, pageno, lpage, scount, ecount, ncount;
  gchar *routine = "ObitSurveyVLSSPrint";

  /* Error checks */
  g_assert(ObitPrinterIsA(printer));
  if (err->error) return FALSE;  /* previous error? */

  /* Get parameters - be somewhat tolerant of the types */
  snprintf (Object, 201, "    ");
  if (ObitInfoListGetTest(data->info, "Object",    &type, dim, Object))
    Object[dim[0]] = 0;  /* Make sure null terminated */
  equinCode = 2;
  if (ObitInfoListGetTest(data->info, "equinCode", &type, dim, &InfoReal)) {
    if (type==OBIT_oint) equinCode = InfoReal.itg;
    else equinCode = InfoReal.otg;   
  }
  fitted = FALSE;
  ObitInfoListGetTest(data->info, "Fitted",    &type, dim, &fitted);
  doraw = FALSE;
  ObitInfoListGetTest(data->info, "doraw",     &type, dim, &doraw);
  ra = 0.0;
  if (ObitInfoListGetTest(data->info, "RA",    &type, dim, &InfoReal)) {
    if (type==OBIT_float)       ra = InfoReal.flt;
    else if (type==OBIT_double) ra = InfoReal.dbl; }  
  dec = 0.0;
  if (ObitInfoListGetTest(data->info, "Dec",       &type, dim, &InfoReal)) {
    if (type==OBIT_float)      dec = InfoReal.flt;
    else if (type==OBIT_double) dec = InfoReal.dbl; }  
  search = 15.0;
  if (ObitInfoListGetTest(data->info, "Search",    &type, dim, &InfoReal)) {
    if (type==OBIT_float)       search = InfoReal.flt;
    else if (type==OBIT_double) search = InfoReal.dbl; }  
  box[0] = 0.0; box[1]=0.0;
  if (ObitInfoListGetTest(data->info, "Box",       &type, dim, &box)) {
    if (type==OBIT_float) {     
      ObitInfoListGetTest(data->info, "Box",       &type, dim, &farr);
      box[0] = (odouble)farr[0];
      box[1] = (odouble)farr[1]; } }
  silent = 720.;
  if (ObitInfoListGetTest(data->info, "Silent",    &type, dim, &InfoReal)) {
    if (type==OBIT_float)      silent = InfoReal.flt;
    else if (type==OBIT_double) silent = InfoReal.dbl; }  
  minFlux = 0.0;
  if (ObitInfoListGetTest(data->info, "minFlux",   &type, dim, &InfoReal)) {
     if (type==OBIT_float)      minFlux = InfoReal.flt;
    else if (type==OBIT_double) minFlux = InfoReal.dbl; }  
  maxFlux = 100000.;
  if (ObitInfoListGetTest(data->info, "maxFlux",   &type, dim, &InfoReal)) {
    if (type==OBIT_float)       maxFlux = InfoReal.flt;
    else if (type==OBIT_double) maxFlux = InfoReal.dbl; }  
  minGlat = 0.0;
  if (ObitInfoListGetTest(data->info, "minGlat",   &type, dim, &InfoReal)) {
    if (type==OBIT_float)      minGlat = InfoReal.flt;
    else if (type==OBIT_double) minGlat = InfoReal.dbl; }  
  maxGlat = 90.;
  if (ObitInfoListGetTest(data->info, "maxGlat",   &type, dim, &InfoReal)) {
    if (type==OBIT_float)       maxGlat = InfoReal.flt;
    else if (type==OBIT_double) maxGlat = InfoReal.dbl; }  

  /* initialize */
  snprintf (datevl, 21, "unknown");
  snprintf (eline, 132, "eline not initialized ");
  found = FALSE;
  numcol = width; /* Output line limit */

  if (first) {
    /* Reset page titles */
    snprintf (Title1, 132, "      ");
    snprintf (Title2, 132, "      ");
    scount = 0;
    ecount = 0;
    ncount = 0;
    page   = 1;
    pageno = 1;
    odiscl = 0.0e0;

    /* Print header info */
    snprintf (line, 132, "      ");
    ObitPrinterWrite (printer, line, &quit, err);
    if (quit) goto Quit;
    if (err->error) goto cleanup;
    snprintf (line,  132, "VLA Low Frequency Sky Survey (VLSS) catalog search, ver 3.0");
    ObitPrinterWrite (printer, line, &quit, err);
    if (quit) goto Quit;
    if (err->error) goto cleanup;
    
    snprintf (line , 132, "Error estimates appear below the value.");
    ObitPrinterWrite (printer, line, &quit, err);
    if (quit) goto Quit;
    if (err->error) goto cleanup;
    
    snprintf (line, 132, "Using VL table %d ", VLVer);
    ObitPrinterWrite (printer, line, &quit, err);
    if (quit) goto Quit;
    if (err->error) goto cleanup;
  } /* end if first */
    
  /* Open input VL table */
  VLTable = newObitTableVLValue (data->name, data, &VLVer, OBIT_IO_ReadOnly, err);
  if (err->error) goto cleanup;
  /* Should be there */
  Obit_retval_if_fail((VLTable!=NULL), err, quit, "VL table %d does not exist", VLVer);
  
  ObitTableVLOpen (VLTable, OBIT_IO_ReadOnly, err);
  if (err->error) goto cleanup;
  /* Create row structure */
  VLRow = newObitTableVLRow(VLTable);
  VLnrow = VLTable->myDesc->nrow;   /* How many rows */

  /* Get information from header */
  beam[0] = VLTable->BeamMajor;
  beam[1] = VLTable->BeamMinor;
  beam[2] = VLTable->BeamPA;
  numind  = VLTable->numIndexed;
  strncpy (datevl, ((ObitImage*)data)->myDesc->date, 10); datevl[10] = 0;
  GetVLIndex (VLTable, &numind, VLindex);
  
  if (beam[0] <= 0.0) beam[0] = 45.0 / 3600.0;
  if (beam[1] <= 0.0) beam[1] = 45.0 / 3600.0;

  /* Galactic limits?*/
  dogal = (minGlat > 0.0)  ||  (maxGlat < 90.0);

  /* Rows to search */
  bci = 1;
  eci = VLnrow;
  eci = MAX (bci, eci);
  inc = 1;
  if (inc <= 0) inc = 1;

  /* Silent search? */
  dosil = silent > 0.0e0;
  silent = MAX (0.0, silent);
  mind2 = 1.0e20;
  sildeg = silent / 3600.0;
  
  /* Position search box */
  rac  = ra;
  decc = dec;
  ra0  = rac * DG2RAD;
  dec0 = decc * DG2RAD;
  radius = search / 3600.0e0;

  /* Search window? */
  doall = (radius <= 0.0)  &&  (box[0] <= 0.0e0)  &&  (box[1] <= 0.0e0);
  /* No radius specified? */
  norad = radius <= 0.0;
  if (norad) radius = MAX (box[0], box[1]);
  radr  = radius;
  radr2 = radr * radr;

  /* No Box? */
  nobox = ((box[0] <= 0.0)  &&  (box[1] <= 0.0));
  if (box[0] <= 0.0) box[0] = MAX (MAX(box[0], box[1]), radius);
  if (box[1] <= 0.0) box[1] = MAX (MAX(box[0], box[1]), radius);
  if (!nobox) {
    /*???         RADIUS = MAX (BOX(1), BOX(2)) */
    radr2 = 1.0e20;
    norad = TRUE;
  } 

  /* RA box fullwidth in hours */
  boxra  = 2.0e0 * box[0];
  boxdec = 2.0e0 * box[1];

  /* Always give distance in sec and  PA */
  discl = 3600.0e0;
  /* All positions */
  if (doall) {
    decbeg = -100;
    decend = 100;
    bcindx = 1;
    ecindx = VLnrow;
    npass = 1;
    beg[0] = bcindx;
    end[0] = ecindx;
    rab = 0.0;
    rae = 360.0;
  } else {
    /* Select position range */
    decbeg = decc - MAX (box[1], sildeg);
    decend = decc + MAX (box[1], sildeg);
    decx = MIN (89.0e0, decc);
    rab = rac - MAX (MAX(radius, box[0]), sildeg) / cos (decx * DG2RAD);
    if (rab < 0.0) rab = rab + 360.0;
    rae = rac + MAX (MAX(radius, box[0]), sildeg) / cos (decx * DG2RAD);
    if (rae > 360.0) rae = rae - 360.0;
    irab = rab / 15.0;
    irae = rae / 15.0;
    irab = MIN (23, MAX (0, irab));
    irae = MIN (24, MAX (0, irae));
    VLindex[24] = VLnrow;
    /* Table indexed? */
    indxed = VLindex[0] > 0;
    if (indxed) {
      bcindx = VLindex[irab];
      ecindx = VLindex[MIN(23, irae)+2];
    } else {
      bcindx = 1;
      ecindx = VLnrow;
    } 

    /* It takes two passes for wrap in  RA range. */
    if (irab > irae) {
      npass = 2;
      beg[0] = bcindx;
      end[0] = VLnrow;
      beg[1] = 1;
      end[1] = ecindx;
    } else {
      npass = 1;
      beg[0] = bcindx;
      end[0] = ecindx;
    } 
  }
 
  /* Tell selection criteria */
  if (first) {
    /* Minimum flux density */
    if (minFlux > 1.0e-5) {
      snprintf (line ,132,  "Selecting sources brighter than %9.1f mJy", 
	       minFlux*1000.0);
      ObitPrinterWrite (printer, line, &quit, err);
      if (quit) goto Quit;
      if (err->error) goto cleanup;
    }

    /* Maximum flux density */
    if (maxFlux < 1.0e5) {
      snprintf (line ,132,  "Selecting sources fainter than %9.1f mJy", 
	       maxFlux*1000.0);
      ObitPrinterWrite (printer, line, &quit, err);
      if (quit) goto Quit;
      if (err->error) goto cleanup;
    } 

    /* Galactic latitude range */
    snprintf (line ,132,  "abs. galactic latitude range %5.1f to %5.1f",
	     minGlat, maxGlat);
    ObitPrinterWrite (printer, line, &quit, err);
    if (quit) goto Quit;
    if (err->error) goto cleanup;
  
    /* Table data and number of entries */
    snprintf (line ,132,  "catalog made on %s (dd/mm/yy) with %7i entries",
	     datevl, VLnrow);
    ObitPrinterWrite (printer, line, &quit, err);
    if (quit) goto Quit;
    if (err->error) goto cleanup;

    /* Indexing */
    if (!doall) {
      if (indxed) {
	snprintf (line, 132, "Table indexed for faster access");
	ObitPrinterWrite (printer, line, &quit, err);
	if (quit) goto Quit;
	if (err->error) goto cleanup;
      } else {
 	snprintf (line, 132, "NOTE:table not indexed");
	ObitPrinterWrite (printer, line, &quit, err);
	if (quit) goto Quit;
	if (err->error) goto cleanup;
      } 
    }
    
    /* Show fitted/deconvolved sizes? */
    if (fitted) {
      if (doraw) snprintf (line, 132, "Displaying raw component size, peak flux density");
      else snprintf (line, 132, "Displaying fitted component size, peak flux density");
    } else {
      snprintf (line, 132, "Displaying deconvolved component size, integrated flux density");
    } 
    ObitPrinterWrite (printer, line, &quit, err);
    if (quit) goto Quit;
    if (err->error) goto cleanup;

    /* Res codes */
    snprintf (line, 132, "residual (res) code; nonblank indicates complex source structure:");
    ObitPrinterWrite (printer, line, &quit, err);
    if (quit) goto Quit;
    if (err->error) goto cleanup;

    snprintf (line, 132, "   p* => high peak, r* => high rms, s* => high integral");
    ObitPrinterWrite (printer, line, &quit, err);
    if (quit) goto Quit;
    if (err->error) goto cleanup;
  
     snprintf (line, 132, "      ");
     ObitPrinterWrite (printer, line, &quit, err);
     if (quit) goto Quit;
     if (err->error) goto cleanup;
  } /* end if first */

  if (!doall) {
    /* Position at input equinox */
    ObitPosLabelUtilRAHMS  (rac,  &rahm[0],  &rahm[1],  &ras);
    ObitPosLabelUtilDecDMS (decc, &decdm[0], &decdm[1], &decs);
    /* Deal with sign of -0 */
    if (decc > 0.0) dsig[0] = '+';
    if (decc < 0.0) dsig[0] = '-';
    dsig[1] = 0;
    decdm[0] = abs(decdm[0]);

    /* Need to change target equinox to J2000? */
    if (equinCode == 1) {
      ObitSkyGeomBtoJ (&rac, &decc);  /* from B1950 to J2000 */
    } else if (equinCode==3) {
      /* Given galactic which are defined in B1950 */
      ObitSkyGeomGal2Eq (&rac, &decc);  /* from Galactic to B1950 */
      ObitSkyGeomBtoJ   (&rac, &decc);  /* from B1950 to J2000 */
    } 
    
    snprintf (line, 132, "      ");
    ObitPrinterWrite (printer, line, &quit, err);
    if (quit) goto Quit;
    if (err->error) goto cleanup;
    
    /* Object label */
    if (!norad) {
      snprintf (line, 132, "%s: search within %8.1f arcsec of  %2.2d %2.2d %7.3f %s%2.2d %2.2d %6.2f",
		Object, radius*3600.0, rahm[0], rahm[1], ras, dsig, decdm[0], decdm[1], decs);
      ObitPrinterWrite (printer, line, &quit, err);
      if (quit) goto Quit;
      if (err->error) goto cleanup;
    }
    
    if (!nobox) {
      snprintf (line, 132, "selection box: %5.2f hr x %6.2f deg with center %2.2d %2.2d %7.3f %s%2.2d %2.2d %6.2f",
		boxra, boxdec, rahm[0], rahm[1], ras, dsig, decdm[0], decdm[1], decs);
      ObitPrinterWrite (printer, line, &quit, err);
      if (quit) goto Quit;
      if (err->error) goto cleanup;
    } 
    
    snprintf (line, 132, "      ");
    ObitPrinterWrite (printer, line, &quit, err);
    if (quit) goto Quit;
    if (err->error) goto cleanup;
    
    /* Check declination */
    if (dec < -30.0) {
      snprintf (line, 132, "WARNING: the VLSS is complete only to declination -30!");
      ObitPrinterWrite (printer, line, &quit, err);
      if (quit) goto Quit;
      if (err->error) goto cleanup;
      snprintf (line, 132, "      ");
      ObitPrinterWrite (printer, line, &quit, err);
      if (quit) goto Quit;
      if (err->error) goto cleanup;
    } 
  } /* end not doall */
  
  /* Page labels */
  if (!doall) {
    /*if (numcol < 93) numcol = 70;*/
    /* Default with pos select */
    if (fitted) {
      snprintf (Title1, 132, "   RA[2000]  Dec[2000] Dist(\")  Peak  Major Minor    PA   Res  Field    X_pix  Y_pix");
      snprintf (Title2, 132, " h  m    s    d  m   s   ori      Jy    \"     \"      deg     ");
    } else{
      snprintf (Title1, 132, "   RA[2000]  Dec[2000] Dist(\")  Flux  Major Minor    PA   Res   Field    X_pix  Y_pix");
      snprintf (Title2, 132, " h  m    s    d  m   s   ori      Jy    \"     \"      deg       ");
    }
  } else {
    /*if (numcol < 87) numcol = 64;*/
    /* Default doall listing */
    if (fitted) {
      snprintf (Title1, 132, "   RA[2000]   Dec[2000]   Peak  Major Minor    PA    Res   Field    X_pix    Y_pix");
      snprintf (Title2, 132, " h  m    s    d  m   s     Jy     \"     \"       deg");
    } else {
      snprintf (Title1, 132, "   RA[2000]   Dec[2000]   Flux  Major Minor    PA   Res    Field    X_pix    Y_pix");
      snprintf (Title2, 132, " h  m    s    d  m   s     Jy     \"     \"      deg ");
    }
  } /* end doall */
  /* Set epoch */
  if (equinCode == 1) {
    if (fitted) 
      snprintf (Title1, 132, "   RA[1950]   Dec[1950]   Peak  Major Minor     PA   Res   Field    X_pix    Y_pix");
    else
      snprintf (Title1, 132, "   RA[1950]   Dec[1950]   Flux  Major Minor     PA   Res   Field    X_pix    Y_pix");
  }  else if (equinCode == 3) {
    if (fitted) 
      snprintf (Title1, 132, "   Gal Long   Gal Lat     Peak  Major Minor     PA   Res   Field    X_pix    Y_pix");
    else
      snprintf (Title1, 132, "   Gal Long   Gal Lat     Flux  Major Minor     PA   Res   Field    X_pix    Y_pix");
  } 
  
  /* give headers */
  ObitPrinterSetTitle (printer, Title1, Title2, err);  /* Page titles at top */
  ObitPrinterWrite (printer, Title1, &quit, err);
  ObitPrinterWrite (printer, Title2, &quit, err);
  if (err->error) goto cleanup;
  
  odiscl = discl;
  /* Passes */
  for (ipass= 1; ipass<=npass; ipass++) { /* loop 500 */
    /* Set range of rows. */
    bc = MIN (MAX (bci, beg[ipass-1]), VLnrow);
    if (ec <= 0) ec = VLnrow;
    ec = MIN (eci,  end[ipass-1]);
    if (inc <= 0) inc = 1;
    if (npass == 1) {
      /* No RA wrap in search range */
      rabeg = rab;
      raend = rae;
    }  else if (ipass == 1) {
      /* Wrap - high hours */
      rabeg = rab;
      raend = 360.0e0;
    } else {
      /* Wrap - low hours */
      rabeg = 0.0e0;
      raend = rae;
    } 
    
    /* Find starting location  */
    bc = VLFindRA (rabeg, VLindex, VLTable, VLRow, err);
    if (err->error) goto cleanup;
    
    /* Print selected rows. */
    for (irow=bc; irow<=ec; irow++) { /* loop 100 */
      /* Read table row */
      ObitTableVLReadRow (VLTable, irow, VLRow, err);
      if (err->error) goto cleanup;
      
      /* Want this one? */
      if (doall) {
	dist2 = 0.0;
      } else {
	/* Quick test? */
	if ((VLRow->Ra2000 >= rabeg)  &&  (VLRow->Ra2000 <= raend)) {
	  
	  /* Check Declination */
	  if ((VLRow->Dec2000 < decbeg)  ||  (VLRow->Dec2000 > decend)) continue;
	  
	  /* In RA range, use full test. */
	  ObitSkyGeomShiftXY (rac, decc, 0.0, VLRow->Ra2000, VLRow->Dec2000, &l, &m);
	  dist2 = l*l + m*m;
	} else if ((ipass == npass)  &&  (VLRow->Ra2000 > raend)) {
	  /* Past RA Range, quit if sorted */
	  if (sort == 1) break;
	  dist2 = 1.0e10;
	} else {
	  /* Before RA range. */
	  dist2 = 1.0e10;
	} 
      } 
      
      /* Closest source in silent window */
      mind2 = MIN (mind2, dist2);
      wanted = dist2  <=  radr2;
      /* Rectangular box? */
      if (!nobox) wanted = wanted  &&  ((VLRow->Dec2000 >= decbeg)  &&  (VLRow->Dec2000 <= decend));
   
      if (wanted) {
	/* Set output, Need to change equinox? */
	rat  = VLRow->Ra2000;
	dect = VLRow->Dec2000;
	if (equinCode == 1) {
	  ObitSkyGeomJtoB (&rat, &dect);  /* from J2000 to B1950 */
	} else if (equinCode==3) {
	  /* Convert to galactic which are defined in B1950 */
	  ObitSkyGeomJtoB   (&rat, &dect);  /* from J2000 to B1950 */
	  ObitSkyGeomEq2Gal (&rat, &dect);  /* from B1950 to Galactic */
	} 

	/* VL row to local variables */
	peak   = VLRow->PeakInt;
	major  = VLRow->MajorAxis;
	minor  = VLRow->MinorAxis;
	posang = VLRow->PosAngle;
	irms   = VLRow->IRMS;
	resrms = VLRow->ResRMS;
	respek = VLRow->ResPeak;
	resflx = VLRow->ResFlux;
	strncpy (field, VLRow->Field,8); field[8] = 0;
	cenx   = VLRow->CenterX;
	ceny   = VLRow->CenterY;
    
	/* Make corrections, get errors */
	VLSSCorErr (&rat, &dect, &peak, &major, &minor, &posang, 
		    irms, beam, fitted, doraw, fblank, 
		    &flux, &eflux, &errra, &errdec, 
		    cmajor, cminor, cpa, emajor, eminor, epa);
	ObitPosLabelUtilRAHMS  (rat, &rahm[0], &rahm[1], &ras);
	ObitPosLabelUtilDecDMS (dect, &decdm[0], &decdm[1], &decs);
	/* Deal with sign of -0 */
	if (dect > 0.0) dsig[0] = '+';
	if (dect < 0.0) dsig[0] = '-';
	dsig[1] = 0;
	decdm[0] = abs(decdm[0]);
	
	/* Distance from location with  possible corrections. */
	dist = sqrt (dist2);
	dist = MIN ((discl * dist), 99999.0);
	
	/* As character string */
	if (dist >= 100.0) {
	  itemp = dist + 0.5;
	  snprintf (cdist, 7, "%6d", itemp);
	}  else if (dist >= 10.0) {
	  snprintf (cdist, 7, "%6.1f", dist);
	} else {
	  snprintf (cdist, 7, "%6.2f", dist);
	} 
	
	/* Position angle from center */
	pa = 57.296 * atan2 (l, m+1.0e-20);
	if (pa < 0) {
	  itemp = pa - 0.5;
	} else {
	  itemp = pa + 0.5;
	} 
	/* As character string */
	snprintf (cposa, 7, "%6d", itemp);
	
	/* Convert flux density in Jy to string */
	if (flux < 99.99) {
	  snprintf (cflux, 8, "%7.3f", flux);
	}  else if (flux < 999.99) {
	  snprintf (cflux, 8, "%7.2f", flux);
	} else {
	  snprintf (cflux, 8, "%7.1f", flux);
	} 
	
	/* Convert flux density error Jy to string */
	if (eflux < 99.9) {
	  snprintf (ceflux, 8, "%7.3f", eflux);
	} else if (flux < 999.) {
	  snprintf (ceflux, 8, "%7.2f", eflux);
	} else {
	  snprintf (ceflux, 8, "%7.1f", eflux);
	} 
	
	/* Goodness of fit code */
	imark = 1;
	tcut   = sqrt (cutt*cutt + (0.01*peak)*(0.01*peak));
	ibadvl = 0;
	if (resrms > tcut) {
	  imark = 2;
	  ibadvl = 10000.0 * resrms + 0.5;
	  if (ibadvl > 9999) ibadvl = 9999;
	} 
	if (fabs(respek) > tcut) {
	  imark = 3;
	  ibadvl = 10000.0 * respek + 0.5;
	  if (ibadvl > 9999) ibadvl = 9999;
	} 
	if (resflx > tcut) {
	  imark = 4;
	  ibadvl = 10000.0 * resflx + 0.5;
	  if (ibadvl > 9999) ibadvl = 9999;
	} 
	if (imark>1) snprintf (chbdvl, 5, "%4d", ibadvl);
	else         snprintf (chbdvl, 5, "    ");
	
	/* Check selection criteria */
	select = (peak  >=  minFlux)  &&  (peak  <=  maxFlux);
	
	/* Filter out really bad fits  Peak should be at least 3 times  RMS residual */
	select = select && ((peak/resrms) > 3.0);
	
	/* Galactic latitude range */
	if (dogal) {
	glat = VLRow->Ra2000;
	glon = VLRow->Dec2000;
	  ObitSkyGeomJtoB   (&glat, &glon);  /* from J2000 to B1950 */
	  ObitSkyGeomEq2Gal (&glat, &glon);  /* from B1950 to Galactic */
	  glat = fabs (glat);
	  select = select  &&  ((glat >= minGlat)  &&  (glat <= maxGlat));
	} 
	if (select) {
	  found = TRUE;
	  /* Create appropriate line */
	  if (!doall) {
	    /* Default with pos select */
	    snprintf (line, 132, "%2.2d %2.2d%6.2f %s%2.2d %2.2d%5.1f%s %s %s %s %s%s   %s%8.2f%8.2f", 
		     rahm[0], rahm[1], ras, dsig, decdm[0], decdm[1], decs, cdist, 
		     cflux, cmajor, cminor, cpa, mark[imark-1], 
		     field, cenx, ceny);
	    /* Errors */
	    snprintf (eline, 132, "     %6.2f       %5.1f%s %s %s %s %s%s", 
		     errra, errdec, cposa, ceflux, emajor, eminor, epa, chbdvl);
	  } else {
	    /* Default doall */
	    snprintf (line, 132, "%2.2d %2.2d%6.2f %s%2.2d %2.2d%5.1f %s %s %s %s %s   %s%8.2f%8.2f", 
		     rahm[0], rahm[1], ras, dsig, decdm[0], decdm[1], decs, 
		     cflux, cmajor, cminor, cpa, mark[imark-1], 
		     field, cenx, ceny);
	    /* Errors */
	    snprintf (eline, 132, "     %6.2f       %5.1f %s %s %s %s %s", 
		     errra, errdec, ceflux, emajor, eminor, epa, chbdvl);
	  } 
	  
	  /* Count number printer */
	  scount++;
	  
	  /* Print line and errors */
	  ObitPrinterWrite (printer, line, &quit, err);
	  if (quit) goto Quit;
	  if (err->error) goto cleanup;
	  if (!doraw) ObitPrinterWrite (printer, eline, &quit, err);
	  if (quit) goto Quit;
	  if (err->error) goto cleanup;
	  
	  /* Don't leave widows and orphans */
	  if (page >= lpage) page = 1000;
	} 
      } /* end if wanted */ 
    } /* end loop  L100 over table */;
  } /* end loop  L500: over passes */;

  /* Stuff at end - was it found? */
  if (!found) {
    /* Silent search window? */
    if (dosil) {
      /* Min source distance in arcmin */
      mindm = 60.0e0 * sqrt (mind2) /  DG2RAD;
      if (mindm < 100.0e0) {
	ecount++;
	snprintf (line, 132, "Source not found  closest source is %8.2f arcmin away", mindm);
      } else {
	ncount++;
	snprintf (line, 132, "Source not found  nothing close");
      } 
    } 
    ObitPrinterWrite (printer, line, &quit, err);
    if (quit) goto Quit;
    if (err->error) goto cleanup;
  } /* end if not found */
  
  /* Finished */
  Quit:
  last = last || quit;  /* Done? */
  
  /* Final summary */
  if (last) {
    /* Number found */
    if (scount > 0) {
      snprintf (line, 132, "Found %9d entries", scount);
    } else {
      snprintf (line, 132, "no sources meeting selection criteria found");
    } 
    if (!quit) ObitPrinterWrite (printer, line, &quit, err);
    if (quit) goto cleanup;
    if (err->error) goto cleanup;
    
    /* How many missed? */
    if (ecount > 0) {
      snprintf (line, 132, "%9d fields has no source", ecount);
      if (!quit) ObitPrinterWrite (printer, line, &quit, err);
      if (quit) goto cleanup;
      if (err->error) goto cleanup;
    } 
    
    /* How many not even close? */
    if (ncount > 0) {
      snprintf (line, 132, "%9d fields have nothing close", ncount);
      if (!quit) ObitPrinterWrite (printer, line, &quit, err);
      if (quit) goto cleanup;
      if (err->error) goto cleanup;
    } 
  } /* end summary on last */ 
  
    /* Close VL table */
 cleanup:
  ObitTableVLClose (VLTable, err);
  VLTable = ObitTableVLUnref(VLTable);
  VLRow   = ObitTableVLRowUnref(VLRow);
  if (err->error) Obit_traceback_val (err, routine, data->name, quit);
  
  return quit;
} /* end of routine ObitSurveyVLSSPrint */ 
  

/*----------------------Private functions ---------------------------*/
/**
 * Apply final calibration and error analysis of NVSS catalog row
 * Routine to apply any corrections to NVSS fitted source parameters  
 * and determine errors.  Unless fitted, the model parameters returned  
 * are deconvolved using the restoring beam.  
 * Now allows for blanking of pflux, qcent and ucent  
 * Now computes effective NVSS point source response rather than using  
 * beam.  
 * Routine translated from the AIPSish corerr.FOR/CORERR  
 * 
 * \param   ra      [in/out] RA (deg)
 * \param   dec     [in/out] Dec (deg)
 * \param   peak    [in/out] Peak Ipol (Jy/beam)
 * \param   major   [in/out] Fitted major axis size (deg)
 * \param   minor   [in/out] Fitted minor axis size (deg)
 * \param   posang  Fitted PA deg
 * \param   qcent   Center Q flux density
 * \param   ucent   Center U flux density
 * \param   pflux   Integrated polarized flux density
 * \param   irms    RMS (sigma) in Ipol.
 * \param   prms    RMS (sigma) in Qpol and Upol.
 * \param   beam    Restoring beam major, minor axes and position angle
 * \param           (all deg) (Now ignored)
 * \param   fitted  If true return fitted values else deconvolved
 * \param   doraw   If true return raw values else bias corrected.
 * \param           Generally, fitted should be true if doraw is.
 * \param   fblank  Magic value blanking value
 * \param   flux    [out] Model peak/integrated Flux density (Jy)
 * \param   eflux   [out] Error in flux
 * \param   pflux   [out] Polarized flux density (Jy)
 * \param           [out] Now derived from qcent, ucent
 * \param   epflux  [out] Error in pflux
 * \param   chpang  [out] Polarization angle or blank if probabliity of a
 * \param           [out] false detection exceeds 2%
 * \param   chepan  [out] Error in chpang or blank
 * \param   errra   [out] Error (sec of time) of Right ascension
 * \param   errdec  [out] Error (asec) of Declination
 * \param   cmajor  [out] Major axis size or limit as string (asec)
 * \param   cminor  [out] Minor axis size or limit as string (asec)
 * \param   cpa     [out] Position angle or blank as string (asec)
 * \param   emajor  [out] Error of major axis size or limit as string (asec)
 * \param   eminor  [out] Error of Minor axis size or limit as string (asec)
 * \param   epa     [out] Error of position angle or blank as string (asec)
 */
static void NVSSCorErr (odouble *ra, odouble *dec, 
			ofloat *peak, ofloat *major, ofloat *minor, ofloat *posang,
			ofloat qcent, ofloat ucent, ofloat *pflux, 
			ofloat irms, ofloat prms, ofloat *beam, 
			gboolean fitted, gboolean doraw, ofloat fblank,
			ofloat *flux, ofloat *eflux, ofloat *epflux, 
			gchar chpang[7], gchar chepan[7], ofloat *errra, ofloat *errdec, 
			gchar cmajor[7], gchar cminor[7], gchar cpa[7], 
			gchar emajor[7], gchar eminor[7], gchar epa[7])
{
  ofloat bemrat, snr, snramp, snrmaj, snrmin, errpek, errpk2, 
    sinc, cosc, tmaj, tmin, errmaj, errmin, errpa, errx2, erry2, 
    tflux, etflux, pamp, errpan, polang, peakp, peakx, perr, 
    decrad, spmjy, psf[3];
  olong   ier;
  gboolean   resolv, hafres[3];

  /* Calibration component of  position error (squared) */
  ofloat  calera = (0.45 / 3600.0) * (0.45 / 3600.0);
  ofloat  calede = (0.56 / 3600.0) * (0.56 / 3600.0);

  /* Mean and uncertainty of CLEAN  bias  */
  ofloat   biasav, biaser;
  /* DnC values */
  ofloat dncbav = 0.00030;
  ofloat dncber = 0.00030;
  /* D values */
  ofloat dbav = 0.00020;
  ofloat dber = 0.00020;

  /* Amplitude calibration uncertainty */
  ofloat calaer = 0.03;

  /* Polarization calibration error (Wild guess) */
  ofloat calper = 0.003;

  /* Position bias */
  odouble biasra = -0.025e0/3600.0;
  odouble biasde =  0.113e0/3600.0;
  
  /* Separate CLEAN bias for D and  DnC  */
  if ((*dec < 77.994)  &&  (*dec > -10.125)) {
    biasav = dbav;
    biaser = dber;
  } else {
    biasav = dncbav;
    biaser = dncber;
  } 

  /* Compute effective PSF */
  decrad = *dec * DG2RAD;
  spmjy = (*peak) * 1000.0;

  NVSSnewPoint (decrad, spmjy, &psf[0], &psf[1], &psf[2]);
  /* Convert sizes to degrees */
  psf[0] /=  3600.0;
  psf[1] /=  3600.0;

  /* Force fitted value to be at least as large as the PSF */
  *major = MAX (*major, psf[0]);
  *minor = MAX (*minor, psf[1]);

  /* Save peak Flux density */
  peakx = (*peak);

  /* Correct position bias */
  if (!doraw) {
    *ra +=  biasra / cos (*dec * DG2RAD);
    *dec +=  biasde;
  } 

  /* Beam ratio */
  bemrat = (*major/psf[0]) * (*minor/psf[1]);

  /* Trap under size beams */
  bemrat = MAX (bemrat, 1.0);

  /* Effective SNRs^2 to account for correlated noise. */
  snr = peakx / irms;

  /* SNR**2 for amplitude errors */
  snramp = ((*major)*(*minor)/(4.0*psf[0]*psf[0])) * 
    pow((1.0+(psf[0]/(*major))*(psf[0]/(*major))), 1.5) *
    pow((1.0+(psf[1]/(*minor))*(psf[1]/(*minor))), 1.5) *
    snr * snr;

  /* SNR**2 for major axis error */
  snrmaj = ((*major)*(*minor)/(4.0*psf[0]*psf[0])) * 
    pow((1.0+(psf[0]/(*major))*(psf[0]/(*major))), 2.5) *
    pow((1.0+(psf[1]/(*minor))*(psf[1]/(*minor))), 0.5) *
    snr * snr;

  /* SNR**2 for minor axis/PA errors */
  snrmin = ((*major)*(*minor)/(4.0*psf[0]*psf[0])) * 
    pow((1.0+(psf[0]/(*major))*(psf[0]/(*major))), 0.5) *
    pow((1.0+(psf[1]/(*minor))*(psf[1]/(*minor))), 2.5) *
    snr * snr;

  /* Bias correct peak flux */
  /*old      IF (.NOT.DORAW) PEAKX = PEAKX + BIASAV - 1.0E-3*(BIASFA / PEAK */
  if (!doraw) peakx += biasav - (irms * irms / peakx);
  
  /* Flux errors, Add bias uncertainty */
  errpek = sqrt ((2.0 * peakx*peakx / snramp) + biaser*biaser);
  
  /* Add Flux calibration  error */
  errpek = sqrt (errpek*errpek + (calaer * peakx)*(calaer * peakx));

  /* Polarization flux, angle and  error */
  if ((qcent != fblank)  &&  (ucent != fblank)) {
    polang = 28.6479 * atan2 (ucent, qcent+1.0e-20);
    pamp = sqrt (qcent*qcent + ucent*ucent);

    /* Bias correction */
    if (!doraw) NVSSpdbias (&pamp, prms);
    /*         PAMP = MAX (1.0E-15, PAMP) */
    errpan = 28.6479 * prms / MAX (1.0e-15, pamp);
    errpan = MIN (errpan, 90.0);
  } else {
    pamp = 0.0;
    polang = 0.0;
    errpan = -1.0;
  } 

  /* Use central polarization rather than integral */
  *pflux = pamp;
  if (!(fitted || doraw)) {
    *pflux *= bemrat;
    prms   *= bemrat;
  } 

  /* Bad hair? (no polarization) */
  if ((qcent == fblank)  ||  (ucent == fblank)) *pflux = fblank;
  
  /* Only quote Polarization angle if  Poln SNR>2.6  Poln uncertainty 
     adding calibration term */
  perr = sqrt (prms*prms * (calper * (*peak))*(calper * (*peak)));
  
  /* Probability of false detection < */
  if (((*pflux) != fblank)  &&  (((*pflux)/perr) > 2.6)) {
    snprintf (chpang, 6, "%6.1f", polang);
    snprintf (chepan, 6, "%6.1f", errpan);
  } else if (*pflux == fblank) {
    snprintf (chpang, 6, "Blank ");
    snprintf (chepan, 6, "      ");
  } else {
    snprintf (chpang, 6, "      ");
    snprintf (chepan, 6, "      ");
  } 
  
  /* Put position angle in range (+/- 90 deg) */
  if (*posang < -90.0) *posang += 180.0;
  if (*posang > 90.0)  *posang -= 180.0;

  /* Error in polarized flux */
  *epflux = prms * 1.41421 * bemrat;
  
  /* Errors */
  sinc = sin (DG2RAD*(*posang));
  cosc = cos (DG2RAD*(*posang));

  /* Trap under size beams */
  tmaj = MAX (*major, psf[1]);
  tmin = MAX (*minor, psf[1]);

  /* Axis sizes include 2% calibration error. */
  errmaj = sqrt (((2.0 * tmaj*tmaj) / snrmaj) +  (0.02*psf[0])*(0.02*psf[0]));
  errmin = sqrt (((2.0 * tmin*tmin) / snrmin) +  (0.02*psf[1])*(0.02*psf[1]));
  errpa = RAD2DG * sqrt (pow(tmaj*tmin/ (tmaj*tmaj - tmin*tmin + 1.0e-20),2) * 
			 4.0 / snrmin);
  errpa = MIN (errpa, 90.0);
  
  /* Position errors */
  errx2 = 2.0 * (2.0 * tmaj*tmaj) / (8.0 * log (2.0) * snrmaj);
  erry2 = 2.0 * (2.0 * tmin*tmin) / (8.0 * log (2.0) * snrmin);

  /* Include calibration */
  *errra  = sqrt (calera + errx2*sinc*sinc + erry2*cosc*cosc);
  *errdec = sqrt (calede + erry2*sinc*sinc + errx2*cosc*cosc);

  /* RA error in seconds of time */
  *errra /= (15.0 * cos (DG2RAD * *dec));
  
  /* Deconvolve */
  NVSSbmval ((*major)*3600.0, (*minor)*3600.0, (*posang),  
	     errmaj*3600.0, errmin*3600.0, errpa,  
	     psf[0]*3600.0, psf[1]*3600.0, psf[2], 
	     cmajor, cminor, cpa, emajor, eminor, epa,  
	     &resolv, hafres, &ier);

  /* Convert to strings */
  if (fitted) {
    /* Fitted values */
    snprintf (cmajor, 7, "%6.1f", (*major)*3600.0);
    snprintf (cminor, 7, "%6.1f", (*minor)*3600.0);
    snprintf (cpa, 7, "%6.1f", (*posang));
    snprintf (emajor, 7, "%6.1f", errmaj*3600.0);
    snprintf (eminor, 7, "%6.1f", errmin*3600.0);
    snprintf (epa, 7, "%6.1f", errpa);
    if (resolv || doraw) {
      *flux  = peakx;
      *eflux = errpek;
    } else {
      /* Correct if unresolved */
      *flux  = peakx * sqrt (bemrat);
      *eflux = sqrt (((calaer * (*flux))*(calaer * (*flux))) + 
		     ((*flux)*(*flux)/snramp) + biaser*biaser);
    } 
  } else {
    /* Total flux density,  Trap resolved on one axis */
    if (hafres[0]) {
      /* Pseudo peak flux */
      if (hafres[1]) {
	/* Major axis only */
	peakp = peakx * sqrt (MAX (1.0, (*minor)/psf[1]));
	tflux = peakp * MAX (1.0, (*major)/psf[0]);
      } else {
	/* Minor axis only */
	peakp = peakx * sqrt (MAX (1.0, (*major)/psf[0]));
	tflux = peakp * MAX (1.0, (*minor)/psf[1]);
      } 
      /* Error in pseudo peak flux */
      errpk2 = sqrt (((calaer * peakp)*(calaer * peakp)) + 
		     (3.0 * peakp*peakp / (2.0 * snramp)) + biaser*biaser);
      /* Integrated flux error */
      etflux = sqrt (tflux*tflux * ((errpk2*errpk2 / (peakp*peakp)) + 
				    (psf[0]/(*major)) * ((errmaj*errmaj/((*major)*(*major))))));
    } else {
      /* Fully resolved values */
            tflux = peakx * bemrat;
            etflux = sqrt (tflux*tflux * ((errpek*errpek/(peakx*peakx)) + 
					  (1.0 / bemrat) * ((errmaj*errmaj/((*major)*(*major))) + 
							    (errmin*errmin/((*minor)*(*minor))))));
    } 
    if (resolv) {
      *flux  = tflux;
      *eflux = etflux;
    } else {
      /* Correct if unresolved */
      *flux  = peakx * sqrt (bemrat);
      *eflux = sqrt (((calaer * (*flux))*(calaer * (*flux))) + 
		     ((*flux)*(*flux)/snramp) + biaser*biaser);
    } 
  } 
  /* Convert units for output */
  *errra  *= 3600.0;
  *errdec *= 3600.0;
} /* end of routine NVSSCorErr */ 

/**
 * Routine to apply any corrections to VLSS Redux fitted source  
 * parameters and determine errors.  
 * Unless fitted, the model parameters returned are deconvolved using  
 * the restoring Beam.  
 * Uses parameters from Wendy 13 Jan. 2012  
 * Routine translated from the AIPSish VLSScorerr.f/CORERR  
 * \param ra      [in/out] RA (deg) 
 * \param dec     [in/out] Dec (deg) 
 * \param peak    [in/out] Peak Ipol (Jy/beam) 
 * \param major   [in/out] Fitted major axis size  (deg)
 * \param minor   [in/out] Fitted minor axis size  (deg)
 * \param posang  Fitted PA 
 * \param irms    RMS (sigma) in Ipol. 
 * \param beam    Restoring beam major, minor axes and position angle 
 *               (all deg) 
 * \param fitted  If true return fitted values else deconvolved 
 * \param doraw   If true return raw values else bias corrected. 
 *                Generally, fitted should be true if doraw is. 
 * \param fblank  Magic value blanking value, arbitrary value for 
 *                internal cionsistency
 * \param flux    [out] Model peak/integrated Flux density (Jy) 
 * \param eflux   [out] Error in flux 
 * \param errra   [out] Error (sec of time) of Right Ascension 
 * \param errdec  [out] Error (asec) of Declination 
 * \param cmajor  [out] Major axis size or limit as string (asec) 
 * \param cminor  [out] Minor axis size or limit as string (asec) 
 * \param cpa     [out] Position angle or blank as string (asec) 
 * \param emajor  [out] Error of major axis size or limit as string (asec) 
 * \param eminor  [out] Error of Minor axis size or limit as string (asec) 
 * \param epa     [out] Error of position angle or blank as string (asec) 
 */
static void VLSSCorErr (odouble *ra, odouble *dec, ofloat *peak, 
			ofloat *major, ofloat *minor, ofloat *posang,
			ofloat irms, ofloat *beam, 
			gboolean fitted, gboolean doraw, ofloat fblank,
			ofloat *flux, ofloat *eflux, 
			ofloat * errra, ofloat *errdec, 
			gchar cmajor[7], gchar cminor[7], gchar cpa[7], 
			gchar emajor[7], gchar eminor[7], gchar epa[7])
{
  ofloat bemrat, snr, snramp, snrmaj, snrmin, errpek, errpk2, sinc, 
    cosc, tmaj, tmin, errmaj, errmin, errpa, errx2, erry2, 
    tflux, etflux, peakx, psf[3], peakp;
  olong   ier;
  gboolean   resolv, hafres[3];

  /* Calibration component of  position error (squared) */
  ofloat     calera = (3.30283069611/3600.0) * (3.30283069611/3600.0);
  ofloat     calede = (3.4931948185/3600.0)  * (3.4931948185/3600.0);

  /* Amplitude calibration uncertainty */
  ofloat calaer = 0.10;

  /* Mean and uncertainty (wild guess) of CLEAN bias Times IRMS */
  ofloat      biasav = 0.66;
  ofloat      biaser = 0.00;

  /* Position bias */
  odouble biasra = -0.609412908554/3600.0;
  odouble biasde = -0.0273719001561/3600.0;

  /* sizes are in deg degrees */
  psf[0] = beam[0];
  psf[1] = beam[1];
  psf[2] = beam[2];
  
   /* Force fitted value to be at least as large as the PSF */
  *major = MAX (*major, psf[0]);
  *minor = MAX (*minor, psf[1]);

  /* Save peak Flux density */
  peakx = (*peak);

  /* Correct position bias */
  if (!doraw) {
    *ra +=  biasra / cos (*dec * DG2RAD);
    *dec +=  biasde;
  } 

  /* Beam ratio */
  bemrat = (*major/psf[0]) * (*minor/psf[1]);

  /* Trap under size beams */
  bemrat = MAX (bemrat, 1.0);

  /* Effective SNRs^2 to account for correlated noise. */
  snr = peakx / irms;

  /* SNR**2 for amplitude errors */
  snramp = ((*major)*(*minor)/(4.0*psf[0]*psf[0])) * 
    pow((1.0+(psf[0]/(*major))*(psf[0]/(*major))), 1.5) *
    pow((1.0+(psf[1]/(*minor))*(psf[1]/(*minor))), 1.5) *
    snr * snr;

  /* SNR**2 for major axis error */
  snrmaj = ((*major)*(*minor)/(4.0*psf[0]*psf[0])) * 
    pow((1.0+(psf[0]/(*major))*(psf[0]/(*major))), 2.5) *
    pow((1.0+(psf[1]/(*minor))*(psf[1]/(*minor))), 0.5) *
    snr * snr;

  /* SNR**2 for minor axis/PA errors */
  snrmin = ((*major)*(*minor)/(4.0*psf[0]*psf[0])) * 
    pow((1.0+(psf[0]/(*major))*(psf[0]/(*major))), 0.5) *
    pow((1.0+(psf[1]/(*minor))*(psf[1]/(*minor))), 2.5) *
    snr * snr;

  /* Bias correct peak flux : DIFFERS from NVSS */
  if (!doraw) peakx += biasav * irms;

  /* Flux errors, Add bias uncertainty : DIFFERS from NVSS */
  errpek = sqrt ((2.0 * peakx*peakx / snramp) + 
		 (biaser*irms)*(biaser*irms));

  /* Add Flux calibration  error */
  errpek = sqrt (errpek*errpek + (calaer * peakx)*(calaer * peakx));

  /* Put position angle in range (+/- 90 deg) */
  if (*posang < -90.0) *posang += 180.0;
  if (*posang > 90.0)  *posang -= 180.0;

  /* Errors */
  sinc = sin (DG2RAD*(*posang));
  cosc = cos (DG2RAD*(*posang));

 /* Trap under size beams */
  tmaj = MAX (*major, psf[1]);
  tmin = MAX (*minor, psf[1]);

  /* Axis sizes include 2% calibration error. */
  errmaj = sqrt (((2.0 * tmaj*tmaj) / snrmaj) +  (0.02*psf[0])*(0.02*psf[0]));
  errmin = sqrt (((2.0 * tmin*tmin) / snrmin) +  (0.02*psf[1])*(0.02*psf[1]));
  errpa = RAD2DG * sqrt (pow(tmaj*tmin/ (tmaj*tmaj - tmin*tmin + 1.0e-20),2) * 
			 4.0 / snrmin);
  errpa = MIN (errpa, 90.0);
  
  /* Position errors */
  errx2 = 2.0 * (2.0 * tmaj*tmaj) / (8.0 * log (2.0) * snrmaj);
  erry2 = 2.0 * (2.0 * tmin*tmin) / (8.0 * log (2.0) * snrmin);

  /* Include calibration, AC fudge factor=1.0, 
     was RA = 1.11863, dec=1.11532 in original processing */
  *errra  = sqrt (calera + errx2*sinc*sinc + erry2*cosc*cosc);
  *errdec = sqrt (calede + erry2*sinc*sinc + errx2*cosc*cosc);

  /* RA error in seconds of time */
  *errra /= (15.0 * cos (DG2RAD * *dec));
  
  /* Deconvolve */
  NVSSbmval ((*major)*3600.0, (*minor)*3600.0, (*posang),  
	     errmaj*3600.0, errmin*3600.0, errpa,  
	     psf[0]*3600.0, psf[1]*3600.0, psf[2], 
	     cmajor, cminor, cpa, emajor, eminor, epa,  
	     &resolv, hafres, &ier);

  /* Convert to strings */
  if (fitted) {
    /* Fitted values */
    snprintf (cmajor, 7, "%6.1f", (*major)*3600.0);
    snprintf (cminor, 7, "%6.1f", (*minor)*3600.0);
    snprintf (cpa, 7, "%6.1f", (*posang));
    snprintf (emajor, 7, "%6.1f", errmaj*3600.0);
    snprintf (eminor, 7, "%6.1f", errmin*3600.0);
    snprintf (epa, 7, "%6.1f", errpa);
    if (resolv || doraw) {
      *flux  = peakx;
      *eflux = errpek;
    } else {
      /* Correct if unresolved : Differs from NVSS?
	 use of biaserr may not be right(?) */
      *flux  = peakx * sqrt (bemrat);
      *eflux = sqrt (((calaer * (*flux))*(calaer * (*flux))) + 
		     ((*flux)*(*flux)/snramp) + 
		     (biaser*irms)*(biaser*irms));
    } 
  } else {
    /* Total flux density,  Trap resolved on one axis */
    if (hafres[0]) {
      /* Pseudo peak flux */
      if (hafres[1]) {
	/* Major axis only */
	peakp = peakx * sqrt (MAX (1.0, (*minor)/psf[1]));
	tflux = peakp * MAX (1.0, (*major)/psf[0]);
      } else {
	/* Minor axis only */
	peakp = peakx * sqrt (MAX (1.0, (*major)/psf[0]));
	tflux = peakp * MAX (1.0, (*minor)/psf[1]);
      } 
      /* Error in pseudo peak flux  : Differs from NVSS? */
      errpk2 = sqrt (((calaer * peakp)*(calaer * peakp)) + 
		     (3.0 * peakp*peakp / (2.0 * snramp)) + 
		     (biaser*irms)*(biaser*irms));
      /* Integrated flux error */
      etflux = sqrt (tflux*tflux * ((errpk2*errpk2 / (peakp*peakp)) + 
				    (psf[0]/(*major)) * ((errmaj*errmaj/((*major)*(*major))))));
    } else {
      /* Fully resolved values */
            tflux = peakx * bemrat;
            etflux = sqrt (tflux*tflux * ((errpek*errpek/(peakx*peakx)) + 
					  (1.0 / bemrat) * ((errmaj*errmaj/((*major)*(*major))) + 
							    (errmin*errmin/((*minor)*(*minor))))));
    } 
    if (resolv) {
      *flux  = tflux;
      *eflux = etflux;
    } else {
      /* Correct if unresolved */
      *flux  = peakx * sqrt (bemrat);
      *eflux = sqrt (((calaer * (*flux))*(calaer * (*flux))) + 
		     ((*flux)*(*flux)/snramp) + (biaser*irms)*(biaser*irms));
    } 
  } 

  /* Convert units for output */
  *errra  *= 3600.0;
  *errdec *= 3600.0;
} /* end of routine VLSSCorErr */ 

/**
 * Looks up the first row in an open VL table with the RA exceeding RA.  
 * Routine translated from the AIPSish NVSSlist.f/FINDRA  
 * \param ra      Desired right ascension (J2000) degrees 
 * \param VLindex Index table giving the first row number with RA 
 *                greater than ra,  for each hour of RA. 
 *                [24] = number of rows 
 * \param err     Error code, 0=> OK 
 * \return the row number in the table
 */
static olong VLFindRA (odouble ra, olong* VLindex, ObitTableVL *VLTable, 
  ObitTableVLRow *VLRow, ObitErr* err) 
{
  olong outrow = 1;
  odouble rahr, frac;
  olong   irow, ihour, inext, tstrow, n;
  gchar *routine = "VLFindRA";

  /* Error checks */
  if (err->error) return outrow;  /* previous error? */
  
  /* First guess  - interpolate */
  rahr   = ra / 15.0e0;
  ihour  = rahr;
  frac   = rahr - ihour;
  inext  = ihour + 1;
  tstrow = VLindex[ihour] + frac * (VLindex[inext] - VLindex[ihour]);
  ObitTableVLReadRow (VLTable, tstrow, VLRow, err);
  if (err->error) Obit_traceback_val (err, routine, VLTable->name, outrow);
  
  /* Coarse Steps of 50 - Go forward or backwards? */
  if (VLRow->Ra2000 < ra) {
    /* Forward */
    n = VLindex[inext] - tstrow + 1;
    for (irow= 1; irow<=n; irow++) {
      if (tstrow > VLindex[24]) {
	tstrow = VLindex[24];
	break;
      } 
      ObitTableVLReadRow (VLTable, tstrow, VLRow, err);
      if (err->error) Obit_traceback_val (err, routine, VLTable->name, outrow);
      /* Past it? */
      if (VLRow->Ra2000 > ra) break;
      tstrow += 50;
    } /* end forward coarse loop */;
  } else {
       
    /* Backwards */
    n = tstrow - VLindex[ihour] + 1;
    for (irow= 1; irow<=n; irow++) {
      if (tstrow < 1) {
	tstrow = 1;
	break;
      } 
      ObitTableVLReadRow (VLTable, tstrow, VLRow, err);
      if (err->error) Obit_traceback_val (err, routine, VLTable->name, outrow);
      /* Past it? */
      if (VLRow->Ra2000 < ra)  break;
      tstrow -= 50;
    } /* end backwards coarse loop */;
  } 

  /* Steps of 1 -  Go forward or backwards? */
  tstrow = MAX (1, MIN (tstrow, VLindex[24]));
  ObitTableVLReadRow (VLTable, tstrow, VLRow, err);
  if (err->error) Obit_traceback_val (err, routine, VLTable->name, outrow);
  if (VLRow->Ra2000 < ra) {
    /* Forward */
    n = VLindex[inext] - tstrow + 1;
    for (irow= 1; irow<=n; irow++) {
      if (tstrow > VLindex[24]) break;
      ObitTableVLReadRow (VLTable, tstrow, VLRow, err);
      if (err->error) Obit_traceback_val (err, routine, VLTable->name, outrow);
      /* Found it? */
      if (VLRow->Ra2000 > ra) return tstrow;
      tstrow++;
    } /* end forward fine loop */;
  } else {
    /* Backwards */
    n = tstrow - VLindex[ihour] + 1;
    for (irow= 1; irow<=n; irow++) {
      if (tstrow < 1) break;
      ObitTableVLReadRow (VLTable, tstrow, VLRow, err);
      if (err->error) Obit_traceback_val (err, routine, VLTable->name, outrow);
      /* Found it? */
      if (VLRow->Ra2000 < ra) return tstrow+1;
      tstrow--;
    } 
  } /* end backwards fine look */;
  
  /* If it gets here punt */
  outrow = VLindex[ihour];
  
  return outrow;
} /* end of routine VLFindRA */ 

/**
 * Extract index information from table
 * \param VLTable    Open VL table
 * \param numind     [out] number of indexed entries
 * \param VLindex    [out] index
 */
static void GetVLIndex(ObitTableVL* VLTable, olong *numind, olong VLindex[25]) 
{
  olong i;

  if (VLTable->numIndexed>0) *numind = VLTable->numIndexed;
  else                       *numind = VLTable->myDesc->nrow;
  VLindex[0]  = MIN (*numind, VLTable->index00);
  VLindex[1]  = MIN (*numind, VLTable->index01);
  VLindex[2]  = MIN (*numind, VLTable->index02);
  VLindex[3]  = MIN (*numind, VLTable->index03);
  VLindex[4]  = MIN (*numind, VLTable->index04);
  VLindex[5]  = MIN (*numind, VLTable->index05);
  VLindex[6]  = MIN (*numind, VLTable->index06);
  VLindex[7]  = MIN (*numind, VLTable->index07);
  VLindex[8]  = MIN (*numind, VLTable->index08);
  VLindex[9]  = MIN (*numind, VLTable->index09);
  VLindex[10] = MIN (*numind, VLTable->index10);
  VLindex[11] = MIN (*numind, VLTable->index11);
  VLindex[12] = MIN (*numind, VLTable->index12);
  VLindex[13] = MIN (*numind, VLTable->index13);
  VLindex[14] = MIN (*numind, VLTable->index14);
  VLindex[15] = MIN (*numind, VLTable->index15);
  VLindex[16] = MIN (*numind, VLTable->index16);
  VLindex[17] = MIN (*numind, VLTable->index17);
  VLindex[18] = MIN (*numind, VLTable->index18);
  VLindex[19] = MIN (*numind, VLTable->index19);
  VLindex[20] = MIN (*numind, VLTable->index20);
  VLindex[21] = MIN (*numind, VLTable->index21);
  VLindex[22] = MIN (*numind, VLTable->index22);
  VLindex[23] = MIN (*numind, VLTable->index23);
  VLindex[24] = *numind;

  /* In case no actual index */
  for (i=0; i<25; i++) VLindex[i] = MAX (1, VLindex[i]);

} /* end GetVLIndex */

/**
 * Estimates the polarization bias in a polarization amplitude, p,  
 * measured in the presence of Q and U RMS noise, rms.  Returns the  
 * corrected value.  
 * The bias correction is such that the average bias is removed;  
 * thus the average in the absence of a signal is zero.  Does table  
 * lookup of values calculated by J. Condon in the range of p/rms of  
 * 1.253 (the no signal limit) and 4 (the high SNR regime).  Does  
 * second order Lagrange interpolation.  At lower values of P/RMS the  
 * bias is a constant 1.253*RMS. Above a signal-to-noise ratio of 4,  
 * use the formula:  
 * normalized bias = 1 / (2 * s) + 1 / (8 * s**3),  
 * where s is the true normalized flux density, iterating once to  
 * estimate s from the normalized map flux density.  "Normalized" means  
 * divided by the rms noise in the q and u maps.  
 * Routine translated from the AIPSish corerr.FOR/PDBIAS  
 *
 * \param p   On output, p is the estimated intrinsic total 
 *             polarized intensity. 
 * \param rms  The standard deviation of the (assumed equal) 
 *             Gaussian distributions of the Stokes Q or U maps. 
 */
static void  NVSSpdbias (ofloat *p, ofloat rms) 
{
  olong   i, index, i1, i2, i3;
  ofloat  pnorm, bias, d1, d2, d3, wt1, wt2, wt3, sum, sumwt;
  /* (map_flux,map_bias) pairs split into table1 and table2 */
  ofloat table1[] = {
      1.253, 1.256, 1.266, 1.281, 1.303, 1.330, 1.364, 1.402, 1.446, 
      1.495, 1.549, 1.606, 1.668, 1.734, 1.803, 1.875, 1.950, 2.027, 
      2.107, 2.189, 2.272, 2.358, 2.444, 2.532, 2.621, 2.711, 2.802, 
      2.894, 2.986, 3.079, 3.173, 3.267, 3.361, 3.456, 3.551, 3.646, 
      3.742, 3.838, 3.934, 4.031};
    ofloat table2[] = {
      1.253,  1.156,  1.066,  0.9814, 0.9030, 0.8304, 0.7636, 0.7023, 
      0.6462, 0.5951, 0.5486, 0.5064, 0.4683, 0.4339, 0.4028, 0.3749, 
      0.3498, 0.3273, 0.3070, 0.2888, 0.2724, 0.2576, 0.2442, 0.2321, 
      0.2212, 0.2112, 0.2021, 0.1938, 0.1861, 0.1791, 0.1726, 0.1666, 
      0.1610, 0.1557, 0.1509, 0.1463, 0.1420, 0.1380, 0.1342, 0.1306};

    /* Check RMS */
    if (rms <= 0.0) return;
    pnorm = (*p) / rms;
    
    /* Which regime? */
    if (pnorm <= table1[0]) {
      
      /* Low (no) SNR case */
      (*p) -= table2[0] * rms;
    } else if (pnorm >= table1[39]) {
      /* High SNR */
      bias = 1.0 / (2.0 * pnorm) + 1.0 / (8.0 * pnorm* pnorm* pnorm);
      pnorm = pnorm - bias;
      bias = 1.0 / (2.0 * pnorm) + 1.0 / (8.0 * pnorm*pnorm*pnorm);
      
      /* Correct for bias */
      *p -= bias * rms;
    } else {
      /* Middle, interpolate in table */
      index = 2;
      for (i= 3; i<=39; i++) {
	if (pnorm < table1[i-1]) break;
	index = i;
      } 
      /* Lagrange interpolation */
      i1 = index - 1;
      i2 = index;
      i3 = index + 1;
      d1 = (table1[i1-1] - table1[i2-1]) * (table1[i1-1] - table1[i3-1]);
      d2 = (table1[i2-1] - table1[i1-1]) * (table1[i2-1] - table1[i3-1]);
      d3 = (table1[i3-1] - table1[i1-1]) * (table1[i3-1] - table1[i2-1]);
      wt1 = (pnorm - table1[i2-1]) * (pnorm - table1[i3-1]) / d1;
      wt2 = (pnorm - table1[i1-1]) * (pnorm - table1[i3-1]) / d2;
      wt3 = (pnorm - table1[i1-1]) * (pnorm - table1[i2-1]) / d3;
      sum = table2[i1-1] * wt1 + table2[i2-1] * wt2 + table2[i3-1] * wt3;
      sumwt = wt1 + wt2 + wt3;
      if (sumwt > 0.0) {
	bias = sum / sumwt;
      } else {
	/* Shouldn't ever get here but do something reasonable. */
	bias = table2[i2-1];
      } 
      /* Correct for bias */
      *p -= bias * rms;
    } 
} /* end of routine NVSSpdbias */ 

/**
 * Subroutine BMVAL deconvolves the fitted beam from the clean beam and  
 * also generates appropriate errors.  
 * Routine translated from the AIPSish corerr.FOR/BMVAL  
 * \param bmaj    Fitted major axis (asec) 
 * \param bmin    Fitted minor axis (asec) 
 * \param bpa     Fitted pos. angle (deg) 
 * \param bmaje   Fitted major axis error (asec) 
 * \param bmine   Fitted minor axis error (asec) 
 * \param bpae    Fitted pos. angle error (deg) 
 * \param cbmaj   Clean beam major axis (asec) 
 * \param cbmin   Clean beam minor axis (asec) 
 * \param cbpa    Clean beam pos. angle (deg) 
 * \param cmajor  string with Major axis in asec 
 * \param cminor  string with Minor axis in asec 
 * \param cpa     string with position angle in deg 
 * \param emajor  string with error of Major axis in asec 
 * \param eminor  string with error of Minor axis in asec 
 * \param epa     string with error of position angle in deg 
 * \param resolv  If true source was resolved. 
 * \param hafres  If true source was resolved on: 
 *                   1) one axis 
 *                   2) Major axis only 
 *                   3) Minor axis only 
 * \param ier     Error return 0-> Can completely deconvolve 
 *                    1-> Cannot deconv some limits 
 *                    2-> Cannot deconv fitted source 
 */
static void NVSSbmval (ofloat bmaj, ofloat bmin, ofloat bpa, 
		       ofloat bmaje, ofloat bmine, ofloat bpae, 
		       ofloat cbmaj, ofloat cbmin, ofloat cbpa, 
		       gchar* cmajor, gchar* cminor, gchar* cpa, 
		       gchar* emajor, gchar* eminor, gchar* epa, 
		       gboolean *resolv, gboolean hafres[3], 
		       olong *ier) 
{
  ofloat  r[3][3], errpa, srcmaj, usrcmj, srcmin, usrcmn, srcpa;
  
  /* Get deconvolved sizes and errors. */
  NVSSdeconv (cbmaj, cbmin, cbpa, bmaj, bmaje, bmin, bmine, bpa, 
	      &srcmaj, &usrcmj, &srcmin, &usrcmn, &srcpa);
  r[0][0] = srcmaj;
  r[1][0] = usrcmj;
  r[2][0] = usrcmj;
  r[0][1] = srcmin;
  r[1][1] = usrcmn;
  r[2][1] = usrcmn;
  r[0][2] = srcpa;
  r[1][2] = bpae;
  r[2][2] = bpae;
  *resolv = FALSE;
  hafres[0] = FALSE;
  hafres[1] = FALSE;
  hafres[2] = FALSE;
  
  /* Convert to strings */
  if ((r[1][0] > 0.0)  &&  (r[1][1] > 0.0)) {
    /* Deconvolved both axes */
    snprintf (cmajor, 7, "%6.1f", r[0][0]);
    snprintf (cminor, 7, "%6.1f", r[0][1]);
    snprintf (cpa, 7, "%6.1f", r[0][2]);
    snprintf (emajor, 7, "%6.1f", (r[1][0] + r[2][0]) * 0.5);
    snprintf (eminor, 7, "%6.1f", (r[1][1] + r[2][1]) * 0.5);
    errpa = MIN (90.0, 0.5 * fabs (r[1][2] + r[2][2]));
    snprintf (epa, 7, "%6.1f", errpa);
    *resolv = TRUE;
  } else if ((r[1][0] <= 0.0)  &&  (r[1][1] <= 0.0)) {
    /* Upper limits only */
    if (r[0][0] <= 99.0) {
      snprintf (cmajor, 7, " <%4.1f", r[0][0]);
    } else {
      snprintf (cmajor, 7, " <%4.0f", r[0][0]);
    } 
    if (r[0][1] <= 99.0) {
      snprintf (cminor, 7, " <%4.1f", r[0][1]);
    } else {
      snprintf (cminor, 7, " <%4.0f", r[0][1]);
    } 
    snprintf (cpa, 7, "      ");
    snprintf (emajor, 7, "      ");
    snprintf (eminor, 7, "      ");
    snprintf (epa, 7, "      ");
 } else {
   /* Some upper limits,  Major axis */
   if (r[1][0] > 0.0) {
     snprintf (cmajor, 7, "%6.1f", r[0][0]);
     snprintf (emajor, 7, "%6.1f", (r[1][0] + r[2][0]) * 0.5);
     *resolv   = TRUE;
     hafres[1] = TRUE;
   } else {
     if (r[0][0] <= 99.0) {
       snprintf (cmajor, 7, " <%4.1f", r[0][0]);
     } else {
       snprintf (cmajor, 7, " <%4.0f", r[0][0]);
     } 
     snprintf (eminor, 7, "      ");
   } 
   /* Minor axis */
   if (r[1][1] > 0.0) {
     snprintf (cminor, 7, "%6.1f", r[0][1]);
     snprintf (eminor, 7, "%6.1f", (r[1][1] + r[2][1]) * 0.5);
     *resolv   = TRUE;
     hafres[2] = TRUE;
   } else {
     if (r[0][1] <= 99.0) {
       snprintf (cminor, 7, " <%4.1f", r[0][1]);
     } else {
       snprintf (cminor, 7, " <%4.0f", r[0][1]);
     } 
     sprintf (eminor, "      ");
   } 
   /* Position angle */
   if (resolv) {
     snprintf (cpa, 7, "%6.1f", r[0][2]);
     errpa = MIN (90.0, 0.5 * abs (r[1][2] + r[2][2]));
     snprintf (epa, 7, "%6.1f", errpa);
   } else {
     snprintf (cpa, 7, "      ");
     snprintf (epa, 7, "      ");
   } 
 } 
  /* One axis resolved? */
  hafres[0] = hafres[1]  ||  hafres[2];
} /* end of routine NVSSbmval */ 

/**
 * To correct the NVSS "Wall problem", NVSSnewPoint  
 * calculates the NVSS point-source response in terms of an  
 * "effective" elliptical Gaussian. 
 * J. J. Condon 2001 Feb 9  
 * Routine translated from the AIPSish corerr.FOR/NEWPS  
 *  
 * \param decrad   J2000 declination (radians) 
 * \param spmjy    raw peak flux density (mJy) 
 * \param psmaj    [out] FWHM major axis (arcsec) 
 * \param psmin    [out] FWHM minor axis (arcsec) 
 * \param pspa     [out] major axis position angle (deg) east of north 
 */
static void NVSSnewPoint (ofloat decrad, ofloat spmjy, 
			  ofloat *psmaj, ofloat *psmin, ofloat *pspa) 
{
  ofloat  zarad, omega, r;
  ofloat vlalat = 0.5948;  /* VLA latitude, radians */


  /* First calculate dirty-beam parameters for the 1 mJy/beam uncleaned "pedestal".
     DnC configuration: */
  if ((decrad <= -0.176719)  ||  (decrad > 1.361255)) {
    /* Minor axis = 51 arcsec, in units of 45 arcsec restoring beam */
    *psmin = 1.133333;

    /* Zenith angle in radians */
    zarad = decrad - vlalat;
    
    /*  Dirty beam solid angle from NVSS raw fits (/aten/ATEN_1/isot.dir/beamsize.sm, 
	BEAMSIZE.OUT) */
    omega = 1.405 + 0.065 / cos (zarad);
    *psmaj = omega / *psmin;
    *pspa  = 0.;
  } else {
    /*  D configuration: minor axis = 54 arcsec, in units of 45 arcsec 
	restoring beam */
    *psmin = 1.2000;
    zarad  = decrad - vlalat;

    /* Beam is elongated north-south near transit */
    *pspa = 0.;

    /* "zenith" band 6, observed off transit: */
    if ((decrad  >  0.437220)  &&  (decrad  <=  0.768818)) {
      zarad = 0.471; /* average zenith angle is 27 deg */
      *pspa = 90.;   /* beam is elongated east-west instead of north-south */
    } 
    omega  = 0.91 + 0.61 / cos (zarad);
    *psmaj = omega / *psmin;
  } 

  /* Next calculate point-source response  for whole source, cleaned down 
     to 1 mJy/beam equation derived in notes 010201 */
  r = spmjy - 1.;
  *psmin = (r + (*psmin)*(*psmin)*(*psmin)) / (r + (*psmin));
  *psmin = sqrt (*psmin) * 45.0;
  *psmaj = (r + (*psmaj)*(*psmaj)*(*psmaj)) / (r + (*psmaj));
  *psmaj = sqrt (*psmaj) * 45.0;
  return;
} /* end of routine NVSSnewPoint */ 

/**
 * Thus subroutine deconvolves a gaussian point-source response from a  
 * fitted elliptical gaussian to yield the gaussian source parameters.  
 * This calculation is from AJ, 109, 2318 and my notes of 010209.  
 * DECONV also determines whether each source axis is significantly  
 * (98% confidence = 2.33 sigma) resolved.  If so, the rms error  
 * in that axis is calculated.  If not, the 98 onfidence upper limit  
 * to that axis is calculated and its rms error is set to -1.  
 * These error calculations are based on the old DECONV subroutine  
 * and the NVSS paper (AJ, 115, 1693).  
 * Routine translated from the AIPSish corerr.FOR/DECONV  
 * 
 * \param ptsmaj   major-axis of pt src response (arcsec) 
 * \param ptsmin   minor-axis of pt src response (arcsec) 
 * \param ptspa    position angle of pt src response (DEG) 
 *                 MUST BE IN RANGE -180 DEG TO +180 DEG 
 * \param fitmaj   major-axis of raw fit (arcsec) 
 * \param ufitmj   rms error in fitted major axis (arcsec)

 * \param fitmin   minor-axis of raw fit (arcsec) 
 * \param ufitmn   rms error in fitted minor axis (arcsec)  

 * \param fitpa    position angle of raw fit (DEG) 
 *                 MUST BE IN RANGE -180 DEG TO +180 DEG 
 * \param srcmaj  [out] major-axis of deconv source (arcsec) 
 *                 or 98% confidence upper limit 
 * \param usrcmj  [out] rms error in srcmaj if resolved,  
 *                      or -1 if srcmaj is an upper limit  
 * \param srcmin  [out] minor-axis of deconv source (arcsec) 
 *                 or 98% confidence upper limit 
 * \param usrcmn  [out] rms error in srcmin if resolved,  
 *                      or -1 if srcmin in an upper limit 
 * \param srcpa   = [out] position angle of deconv source (deg) 
 *                    RETURNED IN RANGE -90 TO +90 deg 
 */
static void NVSSdeconv (ofloat ptsmaj, ofloat ptsmin, ofloat ptspa, 
			ofloat fitmaj, ofloat ufitmj, ofloat fitmin, ofloat ufitmn, 
			ofloat fitpa, 
			ofloat* srcmaj, ofloat* usrcmj, ofloat* srcmin, ofloat* usrcmn, 
			ofloat* srcpa) 
  {
  ofloat   temp, ptprad, fitprd, phirad, ptmjsq, 
    ptmnsq, ftmjsq, ftmnsq, deltar, scmjsq, scmnsq, phird2, cs2phi, 
    sn2phi, denom, xnum, cosphi, sinphi, srcpar,  ptspar, radarg, 
    scmjmx, umsmin, scmnmx, test, umsmaj, upsmaj, upsmin;
  /* Confidence-level parameters for "significant" resolution; 
     2.33 => 98% onfidence */
  ofloat sigmax = 2.33;
  
  /* fix any reversed input major and minor axes */
  if (ptsmaj  <  ptsmin) {
    temp   = ptsmaj;
    ptsmaj = ptsmin;
    ptsmin = temp;
    ptspa += 90.;
    if (ptspa  >  180.) ptspa -= 360.;
  } 
  if (fitmaj  <  fitmin) {
    temp    = fitmaj;
    fitmaj  = fitmin;
    fitmin  = temp;
    fitpa  +=  90.;
    if (fitpa  >  180.) fitpa -=  360.;
  } 
  
  /*  convert PA'S to radians, in range 0 to pi */
  ptprad = ptspa * DG2RAD;
  if (ptprad  <  0.) ptprad += G_PI;
  fitprd = fitpa * DG2RAD;
  if (fitprd  <  0.) fitprd += G_PI;
  
  /* Calculate PA difference phirad (radians) between raw fit 
     and point-source response */
  phirad = fitprd - ptprad;
  
  /*  Make range of phirad from -G_PI/2 to +G_PI/2 */
  if (phirad  >  G_PI/2.)  phirad -= G_PI;
  if (phirad  <  -G_PI/2.) phirad += G_PI;
  cosphi = cos (phirad);
  sinphi = sin (phirad);
  
  /*  Calculate squares of fwhm sizes */
  ptmjsq = ptsmaj * ptsmaj;
  ptmnsq = ptsmin * ptsmin;
  ftmjsq = fitmaj * fitmaj;
  ftmnsq = fitmin * fitmin;
  
  /* Do various special cases for which the general formula is not 
     needed or is not well defined (division by zero) */
  
  /* Ptsmajsq = ptsminsq? (circular pt-src response) */
  if (fabs (ptmjsq - ptmnsq)  <  1.e-3) {
    /* delta is the position angle difference between the deconvolved 
       src major axis and the pt src response major axis */
    deltar = phirad;
    scmjsq = ftmjsq - ptmjsq;
    scmnsq = ftmnsq - ptmnsq;
    goto L90;
  } 
  
  /*  fitmajsq = fitminsq? (circular raw fit) */
  if (fabs (ftmjsq - ftmnsq)  <  1.e-3) {
    deltar = G_PI / 2.;
    scmjsq = ftmjsq - ptmnsq;
    scmnsq = ftmnsq - ptmjsq;
    goto L90;
  }
  
  /*  phirad = 0? (fit parallel to pt-src response) */
  if (fabs (phirad)  <  1.e-5) {
    scmjsq = ftmjsq - ptmjsq;
    scmnsq = ftmnsq - ptmnsq;
    if (scmjsq  >=  scmnsq) {
      deltar = 0.;
    } else {
      temp = scmjsq;
      scmjsq = scmnsq;
      scmnsq = temp;
      deltar = G_PI / 2.;
    } 
    goto L90;
  }
  
  /* phirad = +- pi/2? (fit perpendicular to beam) */
  if ((fabs (phirad - G_PI/2.)  <  1.e-5)  ||  
      (fabs (phirad + G_PI/2.)  <  1.e-5)) {
    scmjsq = ftmjsq - ptmnsq;
    scmnsq = ftmnsq - ptmjsq;
    deltar = phirad;
    goto L90;
  } 
  /*  end of special cases */
  
  /* Calculate deltarad */
  phird2 = 2. * phirad;
  cs2phi = cos (phird2);
  sn2phi = sin (phird2);
  denom = (ftmjsq - ftmnsq) * cs2phi - (ptmjsq - ptmnsq);
  xnum = (ftmjsq - ftmnsq) * sn2phi;
  
  /* Calculate deltarad */
  if (fabs(denom-1)  >=  1.e-5) {
    deltar = 0.5 * atan (xnum / denom);
  } else {
    if (xnum  >  0.)  deltar =  G_PI / 4.;
    if (xnum  <=  0.) deltar = -G_PI / 4.;
  } 
  /* Range is now +- pi/4., resolve ambiguities to make range +- pi/2. */
  if (denom  <  0.) {
    if (xnum  >  0.)  deltar += G_PI / 2.;
    if (xnum  <=  0.) deltar -= G_PI / 2.;
  } 
  /* Calculate srcmajsq */
  scmjsq = (ftmjsq - ftmnsq) * cosphi * sinphi;
  scmjsq = scmjsq / (cos (deltar) * sin (deltar));
  scmjsq = 0.5 * (scmjsq + ftmjsq + ftmnsq - (ptmjsq + ptmnsq));
  
  /* calculate srcminsq */
  scmnsq = ftmjsq + ftmnsq - (ptmjsq + ptmnsq) - scmjsq;
  
  L90:
  /* srcmajsq < srcminsq? */
  if (scmjsq  <  scmnsq) {
    temp   = scmjsq;
    scmjsq = scmnsq;
    scmnsq = temp;
    deltar+= G_PI / 2.;
  } 
  
  /* If deconvolution fails (that is, the square of the source size 
     is negative, set source size negative */
  if (scmjsq  >  0.) {
    *srcmaj = sqrt (scmjsq);
  } else {
    *srcmaj = -sqrt (-scmjsq);
  } 
  if (scmnsq  >  0.) {
    *srcmin = sqrt (scmnsq);
  } else {
    *srcmin = -sqrt (-scmnsq);
  }
  
  /* Calculate source position angle */
  srcpar   = deltar + ptprad;
  *srcpa = srcpar * RAD2DG;
  /* Make sure srcpa in range -90 to +90 deg */
  if (*srcpa  <  -90.) (*srcpa) += 180.;
  if (*srcpa  >=  90.) (*srcpa) -= 180.;
  /* end of fit calculation */
  
  /* Next calculate fit uncertainties 
     test for significant major-axis resolution 
     (see section 5.2.4 in AJ, 115, 1693) */
  ptspar = sqrt (ptsmaj*ptsmaj*cosphi*cosphi + ptsmin*ptsmin*sinphi*sinphi);
  if (fitmaj  <  ptspar) {
    temp = ptspar;
  } else {
    temp = fitmaj;
  } 
  
  /* Is fit major axis > 2.33*sigma + beam projected onto fit major axis? */
  test = temp - ufitmj * sigmax - ptspar;
  if (test  >  0.) {
    /* src major axis is significantly resolved calculate rms error in 
       src major axis NVSS paper, eq. 32a */
    radarg = (temp - ufitmj) * (temp - ufitmj) - ptspar * ptspar;
    umsmaj = (*srcmaj) - sqrt (radarg);
    radarg = (temp + ufitmj) * (temp + ufitmj) - ptspar * ptspar;
    upsmaj = sqrt (radarg) - (*srcmaj);
    *usrcmj = 0.5 * (umsmaj + upsmaj);
  } else {
    /* src major axis is not significantly resolved  calculate source major axis
       98% confidence upper limit */
    scmjmx = temp + sigmax * ufitmj;
    radarg = scmjmx * scmjmx - ptspar * ptspar;
    scmjmx = sqrt (radarg);
    /* Set source size limit to srcmajmax */
    *srcmaj = scmjmx;
    *usrcmj = -1.;
  } 
  /* Test for significant minor-axis resolution */
  ptspar = sqrt (ptsmaj*ptsmaj*sinphi*sinphi + ptsmin*ptsmin*cosphi*cosphi);
  if (fitmin  <  ptspar) {
    temp = ptspar;
  } else {
    temp = fitmin;
  } 
  test = temp - ufitmn * sigmax - ptspar;
  if (test  >  0.) {
    /* src minor axis is significantly resolved calculate rms error in 
       src minor axis */
    radarg = (temp - ufitmn) * (temp - ufitmn) - ptspar * ptspar;
    umsmin = (*srcmin) - sqrt (radarg);
    radarg = (temp + ufitmn) * (temp + ufitmn) - ptspar * ptspar;
    upsmin = sqrt (radarg) - (*srcmin);
    *usrcmn = 0.5 *(upsmin + umsmin);
  } else {
    /* src minor axis is not significantly resolved calculate source 
       minor axis upper limt */
    scmnmx = temp + sigmax * ufitmn;
    radarg = scmnmx * scmnmx - ptspar * ptspar;
    scmnmx = sqrt (radarg);
    
    /*  set source size limit to srcminmax */
    *srcmin = scmnmx;
    *usrcmn = -1.;
  } 
  /*  Make sure larger upper limit is called srcmaj */
  if (((*usrcmn)  ==  -1.)  &&  ((*usrcmj)  ==  -1.)) {
    if (*srcmin  >  *srcmaj) {
      temp = *srcmin;
      *srcmin = *srcmaj;
      *srcmaj = temp;
    } 
  } 
  return;
  } /* end of routine NVSSdeconv */ 
