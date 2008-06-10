/* $Id: ObitAIPSFortran.c,v 1.6 2007/08/31 17:24:03 bcotton Exp $ */
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
/*;Correspondence about this software should be addressed as follows: */
/*;         Internet email: bcotton@nrao.edu.                         */
/*;         Postal address: William Cotton                            */
/*;                         National Radio Astronomy Observatory      */
/*;                         520 Edgemont Road                         */
/*;                         Charlottesville, VA 22903-2475 USA        */
/*--------------------------------------------------------------------*/
#include "ObitAIPSFortran.h"
#include "ObitUV.h"
#include "ObitSystem.h"
#include "ObitImage.h"
#include "ObitTable.h"
#include "ObitImageUtil.h"
#include "ObitUVWeight.h"

/*-------- ObitIO: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitAIPSFortran.c
 * ObitAIPSFortran module function definitions.
 */

/*------------------ Macroes -----------------------------*/
/** 
 * Macro for traceback when an error in a called routine is encountered.
 * Writes traceback info, writes error log, deleted obitErr
 * and returns (no return value).
 * This whould be used when returning from a fortran callable routine.
 * \li err = ObitErr structure to log to.
 * \li me = name of calling routine;
 * \li name = object name.
 */
#define Obit_traceback_end(err,me,name) G_STMT_START{           \
     g_snprintf (err->buffer, OBITERRBUFSIZE,                   \
                "Traceback: routine %s:  object %s", me, name); \
     ObitErrPush (err, OBIT_Traceback, err->buffer);            \
     g_snprintf (err->buffer, OBITERRBUFSIZE,                   \
                 " Occured at file %s: line %d (%s)",           \
                  __FILE__, __LINE__, __PRETTY_FUNCTION__);     \
      ObitErrPush (err, OBIT_Traceback, err->buffer);           \
      ObitErrLog(err); /* log errors */                         \
      err = ObitErrUnref(err); /* free error stack */           \
      return;                                                   \
     }G_STMT_END

/*-----------------File Globals ---------------------------*/
/** ObitSystem info */
static ObitSystem* mySystem = {NULL};

/*---------------Private function prototypes----------------*/
/** Private: Convert AIPS Object to Obit Object */
static Obit* newObitAIPSFortran (AIPSObj AIPSObject, 
				 gchar **inMem, gchar **outMem, ObitErr *err);

/** Private: Copy a list of ObitInfoList items to an AIPS Object */
static void ObitAIPSFortranCopyInfo (ObitInfoList *info, AIPSObj AIPSObject, 
				     gchar **inMem, gchar **outMem, ObitErr *err);

/*---------------AIPS function prototypes----------------*/
/*---------------Public functions---------------------------*/

/**
 * Fortran callable interface routine.
 * Initialize ObitSystem information
 * This should be called ONCE before any Obit routines are called.
 * \li Initialize AIPS disk information if any
 * \li Initialize FITS disk information if any
 * \li Initialize Scratch file list.
 * \param pgmName        Name of program (max 5 char if AIPS)
 * \param lenPgmName     Number of characters in pgmName
 * \param pgmNumber      Version number of program (e.g. POPS number).
 * \param AIPSuser       AIPS user number if using AIPS files
 * \param numberAIPSdisk Number of AIPS disks
 *                       If 0, no AIPS files.
 * \param AIPSdir        List of AIPS directory names as long HOLLERITH string
 * \param lenAIPSdir     Number of characters in each entry in AIPSdir
 * \param numberFITSdisk Number of FITS disks
 *                       If 0, no FITS files.
 * \param FITSdir        List of FITS directory names as long HOLLERITH string
 * \param lenFITSdir     Number of characters in each entry in FITSdir
 * \param F_TRUE         Value of Fortran TRUE (used in Fortran interface)
 * \param F_FALSE        Value of Fortran FALSE
 * \param err            Obit error stack for any error messages.
 */
void obintx_ (const gchar *pgmName, const oint *lenPgmName, const oint *pgmNumber, 
	      const oint *AIPSuser,
	      const oint *numberAIPSdisk, const gchar *AIPSdir, const oint *lenAIPSdir,
	      const oint *numberFITSdisk, const gchar *FITSdir, const oint *lenFITSdir,
	      const oint *F_TRUE, const oint *F_FALSE, oint *ierr)
{
  ObitErr *err = newObitErr();
  gchar *CpgmName, **CAIPSdir, **CFITSdir;
  olong i, j, pgmNo, AUID, nAdisk, nFdisk;

  
  /* convert strings from AIPS HOLLERITH to c */
  CpgmName = g_malloc0(*lenPgmName+1);
  for (i=0; i<*lenPgmName; i++) CpgmName[i] = pgmName[i]; CpgmName[i] = 0;

  pgmNo = *pgmNumber;
  AUID  = *AIPSuser;
  nAdisk = *numberAIPSdisk;
  nFdisk = *numberFITSdisk;

  CAIPSdir = g_malloc0(*numberAIPSdisk*sizeof(gchar*));
  for (i=0; i<*numberAIPSdisk; i++) {
    CAIPSdir[i] = g_malloc0(*lenAIPSdir+1);
    for (j=0; j<*lenAIPSdir; j++) CAIPSdir[i][j] = AIPSdir[i*(*lenAIPSdir)+j];
    CAIPSdir[i*(*lenAIPSdir)+j] = 0;
    /* Convert blanks (should only be trailing) to NULLs */
    for (j=0; j<*lenAIPSdir; j++) if (CAIPSdir[i][j]==' ') CAIPSdir[i][j] = 0;
  }

  CFITSdir = g_malloc0(*numberFITSdisk*sizeof(gchar*));
  for (i=0; i<*numberFITSdisk; i++) {
    CFITSdir[i] = g_malloc0(*lenFITSdir+1);
    for (j=0; j<*lenFITSdir; j++) CFITSdir[i][j] = FITSdir[i*(*lenFITSdir)+j];
    CFITSdir[i*(*lenFITSdir)+j] = 0;
    /* Convert blanks (should only be trailing) to NULLs */
    for (j=0; j<*lenAIPSdir; j++) if (CFITSdir[i][j]==' ') CFITSdir[i][j] = 0;
  }

  /* Call Obit startup routine */
  mySystem = ObitSystemStartup (CpgmName, pgmNo, AUID, nAdisk, CAIPSdir,
				nFdisk, CFITSdir, *F_TRUE, *F_FALSE, err);

  /* error? */
  if (err->error) *ierr = 1; /* something went wrong */
  else *ierr = 0;            /* OK */

  /* Start Obit's AIPS Object system */
  ObitAIPSObjectOBinit (ierr);
  if (*ierr!=0) {
   Obit_log_error(err, OBIT_Error, 
		  "Initialization of Obit AIPS object system failed");
 }

  ObitErrLog(err);           /* Show any error messages */

  /* Cleanup */
  err = ObitErrUnref(err);
  for (i=0; i<*numberAIPSdisk; i++) g_free(CAIPSdir[i]);
  g_free(CAIPSdir);
  for (i=0; i<*numberFITSdisk; i++) g_free(CFITSdir[i]);
  g_free(CFITSdir);
				
} /* end obintx_ */

/**
 * Fortran callable interface routine.
 * Shutdown ObitSystem information
 * This should be called ONCE after all Obit routines are called.
 * Any remaining scratch objects/files are deleted.
 */
void obshtx_ (void)
{
  if (mySystem) ObitSystemShutdown (mySystem);
} /* end obshtx_ */

/**
 * Fortran callable interface routine.
 * Copies uvdata from one object to another optionally applying
 * calibration editing and selection/translation.
 * \param uvin   Name of input uvdata, Members
 * \li   UMAX     R         Maximum acceptable U in wavelengths (default all)
 * \li   VMAX     R         Maximum acceptable V in wavelengths (default all)
 * \li   SOURCS   C(30)*16  Names of up to 30 sources, *=>all
 * \li                      First character of name '-' => all except
 * \li                      those specified.
 * \li   SELQUA   I         Qualifier wanted (-1 => all)
 * \li   SELCOD   C*4       Cal code ('    ')
 * \li   TIMRNG   R(8)      Start day, hour, min, sec, end day, hour,
 * \li                      min, sec. 0's => all
 * \li   UVRNG    R(2)      Minimum and maximum baseline lengths in
 * \li                      1000's wavelengths. 0's => all
 * \li   STOKES   C*4       Stokes types wanted.
 * \li                      'I','Q','U','V','R','L','IQU','IQUV'
 * \li                      '    '=> Leave data in same form as in input.
 * \li   BCHAN    I         First channel number selected, 1 rel. to first
 * \li                      channel in data base. 0 => all
 * \li   ECHAN    I         Last channel selected. 0=>all
 * \li   BIF      I         First IF number selected, 1 rel. to first
 * \li                      IF in data base. 0 => all
 * \li   EIF      I         Last IF selected. 0=>all
 * \li   DOCAL    I         If >0 apply calibration, else not.
 * \li   DOPOL    I         If >0 then correct for feed polarization
 * \li                      based on antenna file info.
 * \li   DOACOR   L         True if autocorrelations wanted (false)
 * \li   DOXCOR   L         True if cross-correlations wanted (true)
 * \li   DOWTCL   L         True if weight calibration wanted.
 * \li   DOFQSL   L         True if FREQSEL random parm present (false)
 * \li   FRQSEL   I         Default FQ table entry to select (-1)
 * \li   SELBAN   R         Bandwidth (Hz) to select (-1.0)
 * \li   SELFRQ   R         Frequency (Hz) to select (-1.0)
 * \li   DOBAND   I         >0 if bandpass calibration. (-1)
 * \li   DOSMTH   L         True if smoothing requested. (false)
 * \li   SMOOTH   R(3)      Smoothing parameters (0.0s)
 * \li   DXTIME   R         Integration time (days). Used when applying
 * \li                      delay corrections to correct for delay error.
 * \li   ANTENS   I(50)     List of antennas selected, 0=>all,
 * \li                      any negative => all except those specified
 * \li   SUBARR   I         Subarray desired, 0=>all
 * \li   FGVER    I         FLAG file version number, if < 0 then
 * \li                      NO flagging is applied. 0 => use highest
 * \li                      numbered table.
 * \li   CLUSE    I         Cal (CL or SN) file version number to apply.
 * \li   BLVER    I         BL Table to apply .le. 0 => none
 * \li   BPVER    I         BP table to apply .le. 0 => none
 * \param uvout  Name of output uvdata
 * \param ierr return code, 0=>OK
*/
void obuvcp_ (AIPSObj uvin, AIPSObj uvout, oint *ierr)
{
  ObitUV  *inUV=NULL, *outUV;
  ObitErr *err = newObitErr();
  oint objnumIn, objnumOut;
  ObitAIPSObjectType AType;
  AIPSKeyDim Adim;
  oint DOACOR, DOXCOR;
  gboolean doACor, doXCor;
  ofloat UMAX, VMAX, UVRANGE[2];
  olong i, CORRTYPE;
  gint32 Odim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar cdum[1], AOnameIn[33], AOnameOut[33];
  gchar *uvinAMem[] = {
    "SOURCES", "SELQUA", "SELCOD", "TIMRNG", "UVRNG", "STOKES", 
    "BCHAN", "ECHAN", "BIF", "EIF", "DOCAL", "DOPOL", "DOACOR", "DOXCOR", 
    "DOWTCL", "DOFQSL", "FRQSEL", "SELBAN", "SELFRQ", "DOBAND", "DOSMTH", "SMOOTH", 
    "DXTIME", "ANTENS", "SUBARR", "FGVER", "CLUSE", "BLVER", "BPVER",
    NULL};
  gchar *uvinOMem[] = {
    "Sources ", "selQua  ", "selCod  ", "timeRange", "UVRange ", "Stokes  ", 
    "BChan   ", "EChan   ", "BIF     ", "EIF     ", "doCalib ", "doPol   ", "doACor  ", "doXCor  ", 
    "doWtCl  ", "doFQSl  ", "freqID  ", "selBan  ", "selFrq  ", "doBand  ", "doSmth  ", "Smooth  ", 
    "DxTime  ", "Antennas", "Subarray", "flagVer ", "gainUse ", "BLVer   ", "BPVer   ",
    NULL};
  gchar *uvoutOMem[] = {NULL};
  gchar *uvoutAMem[] = {NULL};
  gchar *routine = "obuvcp_";
  
  *ierr = 1; /* in case something goes wrong */
  
  /* Input object names to c string */
  for (i=0; i<32; i++) AOnameIn[i] = uvin[i]; AOnameIn[i] = 0;
  for (i=0; i<32; i++) AOnameOut[i] = uvout[i]; AOnameOut[i] = 0;

  /* Get AIPS object numbers */
  ObitAIPSObjectOBname (AOnameIn, &objnumIn, err);
  ObitAIPSObjectOBname (AOnameOut, &objnumOut, err);
  if (err->error) Obit_traceback_end (err, routine, AOnameIn);
  
  /* Get values from AIPS objects */
  ObitAIPSObjectOBget (objnumIn,   "UMAX", &AType, Adim, (gpointer)&UMAX, cdum, err);
  /* Default */
  if (err->error) {
    UMAX = 1.0e20;
    ObitErrClear(err); /* don't really care - we have a default*/
  }
  ObitAIPSObjectOBget (objnumIn,   "VMAX", &AType, Adim, (gpointer)&VMAX, cdum, err);
  /* Default */
  if (err->error) {
    UMAX = 1.0e20;
    ObitErrClear(err);
  }
  ObitAIPSObjectOBget (objnumIn,   "UVRNG", &AType, Adim, (gpointer)UVRANGE, cdum, err);
  /* Default */
  if (err->error) {
    UVRANGE[0] = 0;
    UVRANGE[1] = MAX (UMAX, VMAX);
    ObitErrClear(err);
  }
  ObitAIPSObjectOBget (objnumIn,   "DOACOR", &AType, Adim, (gpointer)&DOACOR, cdum, err);
  /* Default */
  if (err->error) {
    doACor = FALSE;
    ObitErrClear(err);
  } else { /* convert to c */
    doACor = ObitAIPSBooleanF2C (DOACOR);
  }
  ObitAIPSObjectOBget (objnumIn,   "DOXCOR", &AType, Adim, (gpointer)&DOXCOR, cdum, err);
  /* Default */
  if (err->error) {
    doXCor = TRUE;
    ObitErrClear(err);
  } else { /* convert to c */
    doXCor = ObitAIPSBooleanF2C (DOXCOR);
  }

  /* Convert input objects from AIPS to Obit */
  inUV   =  (ObitUV*)newObitAIPSFortran (AOnameIn,  uvinAMem, uvinOMem, err);
  outUV  =  (ObitUV*)newObitAIPSFortran (AOnameOut, uvoutAMem, uvoutOMem, err);
  if (err->error) {
    /* Cleanup */
    inUV   = ObitUVUnref(inUV);
    outUV  = ObitUVUnref(outUV);
    Obit_traceback_end (err, routine, AOnameIn);
  }
  
  /* Convert values to Obit form */
  if (doXCor && doACor) CORRTYPE = 1;
  else if (doXCor && !doACor) CORRTYPE = 0;
  else if (!doXCor && doACor) CORRTYPE = 2;
  else CORRTYPE = 0;
  
  /* Set control values on inUV not done properly in mass copy*/
  Odim[0] = 1;
  ObitInfoListPut (inUV->info, "CorrType",  OBIT_long, Odim, (gpointer)&CORRTYPE, err);
  Odim[0] = 2;
  ObitInfoListPut (inUV->info, "UVRange",  OBIT_float, Odim, (gpointer)UVRANGE, err);
  
  /* Do operation */
  ObitUVCopy (inUV, outUV, err);
  if (err->error) Obit_traceback_end (err, routine, AOnameIn);
  
  /* Convert members back to AIPSish */
  ObitAIPSFortranCopyInfo (inUV->info,  uvin, uvinOMem, uvinAMem, err);
  ObitAIPSFortranCopyInfo (outUV->info, uvout, uvoutOMem, uvoutAMem, err);
  if (err->error) {
    /* Cleanup */
    inUV   = ObitUVUnref(inUV);
    outUV  = ObitUVUnref(outUV);
    Obit_traceback_end (err, routine, AOnameIn);
  }
  
  /* error? */
  if (err->error) *ierr = 1; /* something went wrong */
  else *ierr = 0;            /* OK */
  ObitErrLog(err);           /* Show any error messages */
  
  /* Cleanup */
  err   = ObitErrUnref(err);
  inUV   = ObitUVUnref(inUV);
  outUV  = ObitUVUnref(outUV);
  
} /* end obuvcp_ */

/**
 * Fortran callable interface routine.
 * Determine and apply uniform weighting corrections to uv data
 *  Note: this routine weights by the sum of the weights in a weighting
 *  box rather than the number of counts.
 * \param uv    Name of uvdata, Members
 * \li CHINC     I  Channel increment (1)
 * \li MAXBLINE  R  Maximum baseline length
 * \param image Name of image defining grid size with members
 * \li CHTYPE  C*4  'LINE' or 'SUM ' ('SUM ')
 * \li UVWTFN  C*2  Uniform weighting option ('UN')
 * \li UVTAPER R(2) U, V taper
 * \li UVBOX   I    Uniform weighting box (0)
 * \li UVBXFN  I    box function type (1)
 *                  1=pillbox, 2=linear, 3=exponential, 4=Gaussian
 * \li IMSIZE  I(2) Number of pixels in X,Y direction. no default
 * \li CELLSIZE R(2) Cellspacing in x,y in arcsec no default.
 * \li NCHAV   I    Number of channels to be averaged.
 * \li ROBUST  R    Brigg's weighting factor, -5->Uniform, 5->Natural (0)
 * \param ierr return code, 0=>OK
 */
void obufwt_ (AIPSObj uv, AIPSObj image, oint *ierr)
{
  ObitUV     *UV=NULL;
  ObitErr *err = newObitErr();
  oint objnumUV, objnumImage;
  ObitAIPSObjectType AType;
  AIPSKeyDim Adim;
  oint IDUM, UVBOX, UVBXFN, IMSIZE[2];
  ofloat UVTAPER[2], CELLS[7];
  olong i, nuGrid, nvGrid, WtBox, WtFunc;
  float xCells, yCells, Taper, WtPower, Robust, MaxBaseline, MinBaseline;
  gint32 Odim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar UVWTFN[3], cdum[1], AOname[33];
  gchar *uvOMem[] = {NULL};
  gchar *routine = "obufwt_";

  *ierr = 1; /* in case something goes wrong */

  /* Input object name to c string */
  for (i=0; i<32; i++) AOname[i] = uv[i]; AOname[i] = 0;

  /* Get AIPS object number */
  ObitAIPSObjectOBname (AOname, &objnumUV, err);
  ObitAIPSObjectOBname (image,  &objnumImage, err);
  if (err->error) Obit_traceback_end (err, routine, AOname);

  /* Get values from AIPS objects */
  ObitAIPSObjectOBget (objnumUV,   "MAXBLINE", &AType, Adim, (gpointer)&MaxBaseline, cdum, err);
  /* Default */
  if (err->error) {
    MaxBaseline = 1.0e20;
    ObitErrClear(err);
  }
  ObitAIPSObjectOBget (objnumImage, "UVWTFN",  &AType, Adim, (gpointer)&IDUM, UVWTFN, err);
  /* Default 'UN' */
  if (err->error) {
    UVWTFN[0] = 'U'; UVWTFN[1] = 'N';
    ObitErrClear(err);
  }
  ObitAIPSObjectOBget (objnumImage, "UVTAPER", &AType, Adim, (gpointer)UVTAPER, cdum, err);
  /* Default none*/
  if (err->error) {
    UVTAPER[0] = 0.0;
    UVTAPER[1] = 0.0;
    ObitErrClear(err);
  }
  ObitAIPSObjectOBget (objnumImage, "ROBUST",  &AType, Adim, (gpointer)&Robust, cdum, err);
  /* Default 0 */
  if (err->error) {
    Robust = 0.0;
    ObitErrClear(err);
  }
  ObitAIPSObjectOBget (objnumImage, "UVBOX",   &AType, Adim, (gpointer)&UVBOX, cdum, err);
  /* Default 0 */
  if (err->error) {
    UVBOX = 0;
    ObitErrClear(err);
  }
  ObitAIPSObjectOBget (objnumImage, "UVBXFN",  &AType, Adim, (gpointer)&UVBXFN, cdum, err);
  /* Default  1*/
  if (err->error) {
    UVBXFN = 1;
    ObitErrClear(err);
  }
  ObitAIPSObjectOBget (objnumImage, "CELLSIZE",   &AType, Adim, (gpointer)CELLS, cdum, err);
  /* No Default */
  ObitAIPSObjectOBget (objnumImage, "IMSIZE  ",    &AType, Adim, (gpointer)&IMSIZE, cdum, err);
  /* No Default */
  if (err->error) Obit_traceback_end (err, routine, AOname);

  /* Convert values to Obit form */
  WtFunc = MAX (1, UVBXFN);
  WtPower = 1.0; /* Weighting power */
  if (UVWTFN[1]=='S') WtPower = 0.50;
  if (UVWTFN[1]=='V') WtPower = 0.25;
  if (UVWTFN[1]=='O') WtPower = 0.00;
  if (UVWTFN[0]=='C') WtPower = -WtPower;
  Taper = sqrt (UVTAPER[0]*UVTAPER[0] + UVTAPER[1]*UVTAPER[1]); /* Only one taper */
  xCells = CELLS[0]; /* asec */
  yCells = CELLS[1];
  WtBox  = UVBOX;
  nuGrid = IMSIZE[0];
  nvGrid = IMSIZE[1];

  /* Convert input objects from AIPS to Obit */
  UV  = (ObitUV*)newObitAIPSFortran (AOname, uvOMem, uvOMem, err);
  /* Image only used for member values */
  if (err->error) {
    /* Cleanup */
    UV   = ObitUVUnref(UV);
    Obit_traceback_end (err, routine, AOname);
  }

  /* Open and close UV to test and fully instantiate */
  ObitUVOpen (UV, OBIT_IO_ReadOnly, err);
  ObitUVClose (UV, err);
  if (err->error) {
    /* Cleanup */
    UV   = ObitUVUnref(UV);
    Obit_traceback_end (err, routine, AOname);
  }

  /* Set control values on UV */
  Odim[0] = 1;
  ObitInfoListPut (UV->info, "nuGrid",  OBIT_long,   Odim, (gpointer)&nuGrid, err);
  ObitInfoListPut (UV->info, "nvGrid",  OBIT_long,   Odim, (gpointer)&nvGrid, err);
  ObitInfoListPut (UV->info, "WtBox",   OBIT_long,   Odim, (gpointer)&WtBox, err);
  ObitInfoListPut (UV->info, "WtFunc",  OBIT_long,   Odim, (gpointer)&WtFunc, err);
  ObitInfoListPut (UV->info, "xCells",  OBIT_float, Odim, (gpointer)&xCells, err);
  ObitInfoListPut (UV->info, "yCells",  OBIT_float, Odim, (gpointer)&yCells, err);
  ObitInfoListPut (UV->info, "UVTaper", OBIT_float, Odim, (gpointer)&Taper, err);
  ObitInfoListPut (UV->info, "WtPower", OBIT_float, Odim, (gpointer)&WtPower, err);
  ObitInfoListPut (UV->info, "Robust",  OBIT_float, Odim, (gpointer)&Robust, err);
  ObitInfoListPut (UV->info, "MaxBaseline", OBIT_float, Odim, (gpointer)&MaxBaseline, err);
  ObitInfoListPut (UV->info, "MinBaseline", OBIT_float, Odim, (gpointer)&MinBaseline, err);
  if (err->error) Obit_traceback_end (err, routine, AOname);

  /* Do operation */
  ObitUVWeightData (UV, err);
  if (err->error) {
    /* Cleanup */
    UV   = ObitUVUnref(UV);
    Obit_traceback_end (err, routine, AOname);
  }

  /* Convert members back to AIPSish */
  ObitAIPSFortranCopyInfo (UV->info, AOname, uvOMem, uvOMem, err);
  if (err->error) {
    /* Cleanup */
    UV   = ObitUVUnref(UV);
    Obit_traceback_end (err, routine, AOname);
  }

  /* error? */
  if (err->error) *ierr = 1; /* something went wrong */
  else *ierr = 0;            /* OK */
  ObitErrLog(err);           /* Show any error messages */

  /* Cleanup */
   err   = ObitErrUnref(err);
   UV    = ObitUVUnref(UV);
} /* end obufwt_ */

/**
 * Fortran callable interface routine.
 *   Makes beams or images from a uv data set.  If NFIELD is <= 0 then
 *   a beam is made.  If an image is to be made the normalization factor
 *   is obtained from the beam, if absent, the beam is remade.
 *   The input uvdata is assumed to have been calibrated, selected and
 *   had any uniform weighting corrections applied.
 *   Note: output IMAGE and BEAM objects should exist (full
 *   instantation) prior to call.
 *
 *      Two methods of Fourier transform are available: FFT and DFT.  The
 *   FFT method supports multiple fields and allows images to be a
 *   different size from the beam.  The DFT method does a full 3D DFT but
 *   supports only a single field and the beam must be the same size as
 *   the image.  DFT NYI.
 * \param uvdata  Name of uvdata object
 *                Members:
 * \li STOKES    C*4   Desired Stokes parameter (I)
 * \li UVRANGE   R(2)  UV range in kilo wavelengths (all)
 * \li GUARDBND  R(2)  Fractional guardband around edge of uv grid (0)
 * \param ifield  Field to image: 0 => all nfield (note 1 beam made on
 *                ifield = 1 when doBeam true (1) and DO3DIM false)
 * \param nfield  Number of fields.
 * \param image   Array of image AIPS object names
 *                MUST be previously instantiated.
 *                Members of image(1):
 * \li  FTTYPE    C*4   Fourier transform type 'FFT' or 'DFT'. ('FFT')
 * \li  IMSIZE    I(2,*) Image size per field (no default)
 * \li  CELLSIZE  R(2)  Cellsize in arcseconds in X and Y (no default)
 * \li  CHTYPE    C*4   'LINE',  or 'SUM ' for imaging ('SUM')
 * \li  SHIFT     R(2)  Shift in arcsec (DFT imaging)
 * \li  RASHIFT   R(*)  X position shift in arcseconds per field (0) FFT
 * \li  DECSHIFT  R(*)  Y position shift in arcseconds per field (0) FFT
 * \li  CENTERX   I(*)  Center X pixel position per field (std default)
 * \li  CENTERY   I(*)  Center Y pixel position per field (std default)
 * \li  CTYPX     I     X convolving function type (std default)
 * \li  XPARM     R(10) X convolving function parameters( std default)
 * \li  CTYPY     I     Y convolving function type (std default)
 * \li  YPARM     R(10) Y convolving function parameters (std default)
 * \li  DOZERO    L     IF true do Zero spacing flux (do if value given)
 * \li  ZEROSP    R(5)  Zero spacing parameters (no zero spacing flux)
 * \li  TFLUXG    R     Total flux to be subtracted from ZEROSP (0.0)
 * \li  DOTAPER   L     If true taper (do if non zero taper given)
 * \li  UVTAPER   R(2)  X and Y taper values (no taper)
 * \param beam    Array of beam AIPS object names
  *                MUST be previously instantiated.
*                Members:
 * \li  IMSIZE    I(2)  Size of beam (no default)
 * \li  SUMWTS    R     Sum of weights used for normalization (make beam)
 *                      Set when beam gridded.
 * \param dobeam   if True, make a beam else make an image
 * \param docreate if True, create beams and images underlying AIPS files
 * \param chan     First channel in uv data to image
 * \param nchan    Number of channels to "average" into the image
 * \param imchan   First channel number in output image or beam
 * \param ierr return code, 0=>OK
 */
void obimuv_ (AIPSObj uvdata, oint *ifield, oint *nfield, AIPSObj *image, AIPSObj *beam, 
	      oint *dobeam, oint *docreate, oint *chan, oint* nchan,  oint* imchan, oint *ierr)
{
  ObitUV     *UV=NULL;
  ObitImage  **Image = NULL;
  ObitImage  **Beam  = NULL;
  gboolean  doBeam, doCreate, exist;
  ObitErr *err = newObitErr();
  oint objnumUV, objnumImage;
  olong i, j, channel;
  ObitAIPSObjectType AType;
  AIPSKeyDim Adim;
  olong Adisk, cno, Auser, Aseq;
  gchar Aname[13], Aclass[7], cdum[1], AOUVname[33], AOIMname[33];
  olong nx, ny;
  ofloat yCells, xCells;
  oint IMSIZE[2];
  ofloat CELLS[7];
  olong blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  gint32 Odim[MAXINFOELEMDIM] = {1,1,1,1,1};
  gchar *uvAMem[] = {NULL};
  gchar *uvOMem[] = {NULL};
  gchar *imageAMem[] = {"DISK", "CNO", NULL};
  gchar *imageOMem[] = {"DISK", "CNO", NULL};
  gchar *beamAMem[] = {"DISK", "CNO", "SUMWTS", NULL};
  gchar *beamOMem[] = {"DISK", "CNO", "SUMWTS", NULL};
  gchar *routine = "obimuv_";

  *ierr = 1; /* in case something goes wrong */

   /* Input object names to c string */
  for (i=0; i<32; i++) AOUVname[i] = uvdata[i]; AOUVname[i] = 0;
  for (i=0; i<32; i++) AOIMname[i] = image[0][i]; AOIMname[i] = 0;

 /* Get AIPS object number */
  ObitAIPSObjectOBname (AOUVname, &objnumUV, err);
  ObitAIPSObjectOBname (AOIMname, &objnumImage, err);
  if (err->error) Obit_traceback_end (err, routine, AOUVname);

  /* Translate inputs from the AIPSish */
  doBeam = ObitAIPSBooleanF2C (*dobeam);
  doCreate = ObitAIPSBooleanF2C (*docreate);

  /* Convert input objects from AIPS to Obit */
  UV  = (ObitUV*)newObitAIPSFortran (AOUVname, uvAMem, uvOMem, err);

  /* Open and close UV to test and fully instantiate */
  ObitUVOpen (UV, OBIT_IO_ReadOnly, err);
  ObitUVClose (UV, err);
  if (err->error) {
    /* Cleanup */
    UV   = ObitUVUnref(UV);
    Obit_traceback_end (err, routine, AOUVname);
  }

  /* Images */
  Image = g_malloc0(*nfield*sizeof(ObitImage*));
  Beam  = g_malloc0(*nfield*sizeof(ObitImage*));

  /* Need to create? */
  if (doCreate) {
    /* how big is the image? */
    ObitAIPSObjectOBget (objnumImage, "CELLSIZE", &AType, Adim, (gpointer)CELLS, cdum, err);
    ObitAIPSObjectOBget (objnumImage, "IMSIZE  ", &AType, Adim, (gpointer)&IMSIZE, cdum, err);
    
    /* Convert values to Obit form */
    xCells = CELLS[0]; /* asec  */
    yCells = CELLS[1];
    nx = IMSIZE[0];
    ny = IMSIZE[1];
    
    /* Set control values on UV */
    Odim[0] = 1;
    ObitInfoListPut (UV->info, "nx",      OBIT_long,   Odim, (gpointer)&nx, err);
    ObitInfoListPut (UV->info, "ny",      OBIT_long,   Odim, (gpointer)&ny, err);
    ObitInfoListPut (UV->info, "xCells",  OBIT_float, Odim, (gpointer)&xCells, err);
    ObitInfoListPut (UV->info, "yCells",  OBIT_float, Odim, (gpointer)&yCells, err);
    
    /* Get disk, user info from first disk */
    ObitAIPSObjectOBget (objnumUV, "USERNO", &AType, Adim, (gpointer)&Auser, cdum, err);
    ObitAIPSObjectOBget (objnumImage, "DISK", &AType, Adim, (gpointer)&Adisk, cdum, err);
    ObitAIPSObjectOBget (objnumImage, "IMSEQ", &AType, Adim, (gpointer)&Aseq, cdum, err);
    Aseq = MAX (1, Aseq); /* Mama AIPS can get really bent out of shape */
    ObitAIPSObjectOBget (objnumImage, "NAME", &AType, Adim, (gpointer)Aname, Aname, err);
    Aname[12] = 0;
    if (err->error) {
      /* Cleanup */
      UV   = ObitUVUnref(UV);
      Obit_traceback_end (err, routine, AOUVname);
    }
    
    for (i=0; i<*nfield; i++) { /* loop over fields */
      Image[i] = ObitImageUtilCreateImage (UV, i+1, TRUE, err);
      Beam[i] = ObitImageRef (Image[i]->myBeam);

      /* Get AIPS assignments - image */
      /* Input object name to c string */
      for (j=0; j<32; j++) AOIMname[j] = image[i][j]; AOIMname[j] = 0;
      ObitAIPSObjectOBname (AOIMname, &objnumImage, err);
      /* get AIPS class */
      ObitAIPSObjectOBget (objnumImage, "CLASS", &AType, Adim, (gpointer)Aclass, Aclass, err);
      Aclass[6] = 0;
      
      /* assign CNO */
      cno = ObitAIPSDirAlloc(Adisk, Auser, Aname, Aclass, "MA", Aseq, &exist, err);
      if (cno<0) Obit_log_error(err, OBIT_Error,
				"Failure assigning CNO for %s %s %d",Aname, Aclass, Aseq);
      ObitImageSetAIPS(Image[i],OBIT_IO_byPlane,Adisk,cno,Auser,blc,trc,err);
      if (err->error) {
	/* Cleanup */
	UV   = ObitUVUnref(UV);
	Obit_traceback_end (err, routine, AOIMname);
      }
 
     /* Get AIPS assignments - beam */
      /* Input object name to c string */
      for (j=0; j<32; j++) AOIMname[j] = beam[i][j]; AOIMname[j] = 0;
      ObitAIPSObjectOBname (AOIMname, &objnumImage, err);
      /* get AIPS class */
      ObitAIPSObjectOBget (objnumImage, "CLASS", &AType, Adim, (gpointer)Aclass, Aclass, err);
      Aclass[6] = 0;
      
      /* assign CNO */
      cno = ObitAIPSDirAlloc(Adisk, Auser, Aname, Aclass, "MA", Aseq, &exist, err);
      if (cno<0) Obit_log_error(err, OBIT_Error,
				"Failure assigning CNO for %s %s %d",Aname, Aclass, Aseq);
      ObitImageSetAIPS(Beam[i],OBIT_IO_byPlane,Adisk,cno,Auser,blc,trc,err);
      if (err->error) {
	/* Cleanup */
	UV   = ObitUVUnref(UV);
	Obit_traceback_end (err, routine, AOIMname);
      }
    } /* end loop creating images */

  } else { /* Images/Beams already exist - convert to Obit Objects */
    for (i=0; i<*nfield; i++) { /* loop over fields */
      Image[i] = (ObitImage*)newObitAIPSFortran (image[i], imageAMem, imageOMem, err);
      Beam[i]  = (ObitImage*)newObitAIPSFortran (beam[i], beamAMem, beamOMem, err);
      
      /* Attach beam to image*/
      if (Image[i]) Image[i]->myBeam = ObitImageRef (Beam[i]);
    }
    if (err->error) Obit_traceback_end (err, routine, AOIMname);
  }

  /* Do operation - loop over fields */
  channel = *chan;
  for (i=0; i<*nfield; i++) {
    /* Make image */
    ObitImageUtilMakeImage (UV, Image[i], channel, doBeam, FALSE, err);
  }
  if (err->error) {
    /* Cleanup */
    UV    = ObitUVUnref(UV);
    /* Images */
    for (i=0; i<*nfield; i++) {
      Image[i] = ObitImageUnref(Image[i]);
      Beam[i]  = ObitImageUnref(Beam[i]);
    }
    g_free(Image);
    g_free(Beam);
    Obit_traceback_end (err, routine, AOUVname);
  }

  /* Convert members back to AIPSish */
  ObitAIPSFortranCopyInfo (UV->info, uvdata, uvOMem, uvAMem, err);
  /* Images */
  for (i=0; i<*nfield; i++) {
    ObitAIPSFortranCopyInfo (Image[i]->info, image[0], imageOMem, imageAMem, err);
    ObitAIPSFortranCopyInfo (Beam[i]->info, beam[0], beamOMem, beamAMem, err);
  }
  if (err->error) {
    /* Cleanup */
    UV    = ObitUVUnref(UV);
    /* Images */
    for (i=0; i<*nfield; i++) {
      Image[i] = ObitImageUnref(Image[i]);
      Beam[i]  = ObitImageUnref(Beam[i]);
    }
    g_free(Image);
    g_free(Beam);
    Obit_traceback_end (err, routine, AOUVname);
  }

  /* error? */
  if (err->error) *ierr = 1; /* something went wrong */
  else *ierr = 0;            /* OK */
  ObitErrLog(err);           /* Show any error messages */

  /* Cleanup */
   err   = ObitErrUnref(err);
   UV    = ObitUVUnref(UV);
  /* Images */
  for (i=0; i<*nfield; i++) {
    Image[i] = ObitImageUnref(Image[i]);
    Beam[i]  = ObitImageUnref(Beam[i]);
 }
  g_free(Image);
  g_free(Beam);
  
} /* end obimuv_ */

/**
 * Fortran callable interface routine.
 *   Optionally does Ionospheric calibration for each field using IN
 *   table attached to uvdata.
 *   Makes beams or images from a uv data set.  If nfield is <= 0 then
 *   a beam is made.  If an image is to be made the normalization factor
 *   is obtained from the beam, if absent, the beam is remade.
 *   The input uvdata is assumed to have been calibrated, selected and
 *   had any uniform weighting corrections applied.
 *   Note: output image and beam objects should exist (full
 *   instantation) prior to call.
 *
 *      Two methods of Fourier transform are available: FFT and DFT.  The
 *   FFT method supports multiple fields and allows images to be a
 *   different size from the beam.  The DFT method does a full 3D DFT but
 *   supports only a single field and the beam must be the same size as
 *   the image.  DFT NYI.
 * \param uvdata  Name of uvdata object
 *                Members:
 * \li STOKES    C*4   Desired Stokes parameter (I)
 * \li UVRANGE   R(2)  UV range in kilo wavelengths (all)
 * \li GUARDBND  R(2)  Fractional guardband areound edge of uv grid (0)
 * \li  DOIONS    L     If present and true then do ionospheric
 * \li                  calibration using NI table IONTAB.
 * \li    If DOINONS then the following:
 * \li  IONTAB    C*32  Name of IoNosphere table 
 * \li  UVINSCR   C*32  Name of a scratch uv data file.
 * \li  NSUBARR   I     Number of subarrays in data to image
 * \param ifield  Field to image: 0 => all nfield (note 1 beam made on
 *                ifield = 1 when doBeam true (1) and DO3DIM false)
 * \param nfield  Number of fields.
 * \param image   Array of image AIPS object names
 *                Members of image(1):
 * \li  FTTYPE    C*4   Fourier transform type 'FFT' or 'DFT'. ('FFT')
 * \li  IMSIZE    I(2,*) Image size per field (no default)
 * \li  CELLSIZE  R(2)  Cellsize in arcseconds in X and Y (no default)
 * \li  CHTYPE    C*4   'LINE',  or 'SUM ' for imaging ('SUM')
 * \li  SHIFT     R(2)  Shift in arcsec (DFT imaging)
 * \li  RASHIFT   R(*)  X position shift in arcseconds per field (0) FFT
 * \li  DECSHIFT  R(*)  Y position shift in arcseconds per field (0) FFT
 * \li  CENTERX   I(*)  Center X pixel position per field (std default)
 * \li  CENTERY   I(*)  Center Y pixel position per field (std default)
 * \li  CTYPX     I     X convolving function type (std default)
 * \li  XPARM     R(10) X convolving function parameters( std default)
 * \li  CTYPY     I     Y convolving function type (std default)
 * \li  YPARM     R(10) Y convolving function parameters (std default)
 * \li  DOZERO    L     IF true do Zero spacing flux (do if value given)
 * \li  ZEROSP    R(5)  Zero spacing parameters (no zero spacing flux)
 * \li  TFLUXG    R     Total flux to be subtracted from ZEROSP (0.0)
 * \li  DOTAPER   L     If true taper (do if non zero taper given)
 * \li  UVTAPER   R(2)  X and Y taper values (no taper)
 * \param beam    Array of beam AIPS object names
 *                Members:
 * \li  IMSIZE    I(2)  Size of beam (no default)
 * \li  SUMWTS    R     Sum of weights used for normalization (make beam)
 *                      Set when beam gridded.
 * \param dobeam  if True, make a beam else make an image
 * \param docreate if True, create beams and images underlying AIPS files
 * \param chan    First channel in uv data to image
 * \param nchan   Number of channels to "average" into the image
 * \param imchan  First channel number in output image or beam
 * \param ierr return code, 0=>OK
 */
void obiuvi_ (AIPSObj uvdata, oint *ifield, oint *nfield, AIPSObj *image, AIPSObj *beam, 
	      oint *dobeam,  oint *docreate, oint *chan, oint* nchan,  oint* imchan, oint *ierr)
{
  gchar *routine = "obiuvi_";
  g_error ("%s: Write me", routine);
} /* end obiuvi_ */

/*---------------Private functions---------------------------*/
/**
 * Convert an AIPS Object to an Obit object of the appropriate type.
 * Copies a list of member names from the input to output object, renaming them.
 * Output object name will be AIPSObject.
 * \li AIPS: "UVDATA  " becomes ObitUV
 * \li AIPS: "IMAGE   " becomes ObitImage
 * \li AIPS: "TABLE   " becomes ObitTable
 * \param AIPSObject Name of AIPS object
 * \param inMem      NULL terminated list of member names to copy
 * \param outMem     NULL terminated list output names on ObitInfoList
 * \param err        Error stack
 * \return new Obit Object
 */
static Obit* newObitAIPSFortran (AIPSObj AIPSObject, 
				 gchar **inMem, gchar **outMem, ObitErr *err)
{
  Obit      *out=NULL;
  ObitUV    *outUV=NULL;
  ObitImage *outImage=NULL;
  ObitTable *outTable=NULL;
  ObitInfoList *info=NULL;
  oint i, objnum, classno, disk, cno, ver, nvis, nrow, user;
  olong size, need=0, needStr;
  gchar **strArr=NULL, tab[3], **tstMem, *value, AOname[33];
  gboolean *bvalue, grumble;
  AIPSObjClass cName;
  AIPSKeyDim Adim;
  olong Odim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitAIPSObjectType AType;
  ObitInfoType OType=0;
  olong blc[IM_MAXDIM] = {1,1,1,1,1,1,1};
  olong trc[IM_MAXDIM] = {0,0,0,0,0,0,0};
  gchar *routine = "newObitAIPSFortran";

 /* error tests */
  g_assert(ObitErrIsA(err)); 
  if (err->error) return out;  /* existing error condition */

  /* Input object name to c string */
  for (i=0; i<32; i++) AOname[i] = AIPSObject[i]; AOname[i] = 0;

  /* get disk and cno number. */
  ObitAIPSObjectOBdskc (AIPSObject, &disk, &cno, err);
  if (err->error) Obit_traceback_val (err, routine, AOname, out);
  
  /* Get AIPS object number */
  ObitAIPSObjectOBname (AIPSObject, &objnum,  err);
  if (err->error) Obit_traceback_val (err, routine, AOname, out);
  
  /* Get class */
  ObitAIPSObjectOBclass (objnum, &classno, cName, err);
  if (err->error) Obit_traceback_val (err, routine, AOname, out);

  /* Get AIPS user ID */
  ObitAIPSObjectOBget (objnum, "USERNO", &AType, Adim, (gpointer)&user, value, err);
  if (err->error) Obit_traceback_val (err, routine, AOname, out);
  
  if (!strncmp (cName, "UVDATA  ", 8)) {
    /* UV Data */ 
    outUV = newObitUV (AOname);
    out = (Obit*)outUV;    /* Set return value */
    info = outUV->info;  /* ObitInfoList pointer */

    /* Set file info */
    nvis = 1000; /* number of visibilities per I/O */
    ObitUVSetAIPS(outUV,nvis,disk,cno,user,err);
    if (err->error) Obit_traceback_val (err, routine, AOname, out);
    
  } else if (!strncmp (cName, "IMAGE   ", 8)) {
    /* Image */
    outImage = newObitImage (AOname);
    out = (Obit*)outImage;    /* Set return value */
    info = outImage->info;  /* ObitInfoList pointer */
    
    /* Set file info */
    ObitImageSetAIPS(outImage,1000,disk,cno,user,blc,trc,err);
    if (err->error) Obit_traceback_val (err, routine, AOname, out);
    
  } else if (!strncmp (cName, "TABLE   ", 8)) {
    /* Table */
    outTable = newObitTable (AOname);
    out = (Obit*)outTable;    /* Set return value */
    info = outTable->info;  /* ObitInfoList pointer */
    
    /* Get AIPS table type */
    ObitAIPSObjectOBget (objnum, "TBLTYPE", &AType, Adim, (gpointer)value, tab, err);
    if (err->error) Obit_traceback_val (err, routine, AOname, out);

    /* Get AIPS table version */
    ObitAIPSObjectOBget (objnum, "VER", &AType, Adim, (gpointer)&ver, value, err);
    if (err->error) Obit_traceback_val (err, routine, AOname, out);

    /* Set file info */
    nrow = 1; /* one row at a time */
    ObitTableSetAIPS(outTable,disk,cno,tab,ver,user,nrow,err);
    if (err->error) Obit_traceback_val (err, routine, AOname, out);
    
  } else { /* unknown */
    g_error ("%s: Unknown AIPS class %s", routine, cName);
  }

  /* Read through member list to find largest and create scratch array */
  size = sizeof(odouble);
  needStr = 0;
  tstMem = inMem;
  while (*tstMem) {
    /* Read from AIPS */
    if (ObitAIPSObjectOBinfo (objnum, *tstMem, &AType, Adim, err)) {
      if (err->error) Obit_traceback_val (err, routine, AOname, out);

      /* How big is an element? */
      if (AType == OBIT_AIPSObjectInt) {
	need = sizeof(oint);
      } else if (AType == OBIT_AIPSObjectRe) {
	need = sizeof(ofloat);
      } else if (AType == OBIT_AIPSObjectDP) {
	need = sizeof(odouble);
      } else if (AType == OBIT_AIPSObjectCar) {
	need = 1;
	/* Strings are a fucking pain - only handle 2D arrays here */
	needStr = MAX (needStr, Adim[1]);
      } else if (AType == OBIT_AIPSObjectLog) {
	need = sizeof(gboolean);
      }
      
      /* multiply by dimensions */
      for (i=0; i<5; i++) need *= MAX (1, Adim[i]);
      
      size = MAX (size, need);
    } /* end if exist */
    tstMem++;
  } /* end loop over input members */

  /* allocate scratch array for copying */
  value = g_malloc0(size);
  bvalue = (gboolean*)value;

  /* Array of string pointers */
  strArr = g_malloc0(needStr*sizeof(gchar*));
	     
  /* Copy member values to list */
  while (*inMem) {
    /* Read from AIPS */
    if (ObitAIPSObjectOBinfo (objnum, *inMem, &AType, Adim, err)) {
      ObitAIPSObjectOBget (objnum, *inMem, &AType, Adim, (gpointer)value, value, err);
      if (err->error) {
	/* Cleanup */
	if (strArr) g_free(strArr);
	g_free(value);
	Obit_traceback_val (err, routine, AOname, out); /* bye, bye */
      }

      /* Convert to Obit */
      grumble = FALSE;
      if (AType == OBIT_AIPSObjectInt) {
	OType = OBIT_oint;
      } else if (AType == OBIT_AIPSObjectRe) {
	OType = OBIT_float;
      } else if (AType == OBIT_AIPSObjectDP) {
	OType = OBIT_double;
      } else if (AType == OBIT_AIPSObjectCar) {
	OType = OBIT_string;
	/* grumble, grumble
	   if 2D need to pass an array of string pointers. */
	if (Adim[1]>1) {
	  grumble = TRUE;
	  for (i=0; i<Adim[1]; i++) strArr[i] = &value[i*Adim[0]];
	}
      } else if (AType == OBIT_AIPSObjectLog) {
	OType = OBIT_bool;
	/* Convert single boolean */
	bvalue[0] = ObitAIPSBooleanF2C ((oint)bvalue[0]);
      }
      for (i=0; i<5; i++) Odim[i] = Adim[i];
      
      /* Save in InfoList */
      if (grumble) /* string array */
	ObitInfoListPut (info, *outMem, OType, Odim, (gpointer)strArr, err);
      else /* non-string */
	ObitInfoListPut (info, *outMem, OType, Odim, (gpointer)&value, err);
      if (err->error) {
	/* Cleanup */
	if (strArr) g_free(strArr);
	g_free(value);
	Obit_traceback_val (err, routine, AOname, out);
      }
      
    } /* end copy if exist */
    inMem++;
    outMem++;
  } /* end loop over input members */

  /* Cleanup */
  g_free(value);
  if (strArr) g_free(strArr);

  return out;
} /* end newObitAIPSFortran */

/**
 *  Copies a list of ObitInfoItems to output members, possible renaming.
 * \param list       List to copy from.
 * \param AIPSObject Name of AIPS object
 * \param inMem      NULL terminated list of ObitInfoList names to copy
 * \param outMem     NULL terminated list output member names.
 * \param err        Error stack
 * \return new Obit Object
 */
static void ObitAIPSFortranCopyInfo (ObitInfoList *info, AIPSObj AIPSObject, 
  gchar **inMem, gchar **outMem, ObitErr *err)
{
  oint i, objnum;
  olong size, need;
  gchar **tstMem, *value, AOname[33];
  oint *Lvalue;
  AIPSKeyDim Adim;
  olong Odim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitAIPSObjectType AType=0;
  ObitInfoType OType;
  gchar *routine = "ObitAIPSFortranCopyInfo";

 /* error tests */
  g_assert(ObitErrIsA(err)); 
  if (err->error) return;  /* existing error condition */

  /* Input object name to c string */
  for (i=0; i<32; i++) AOname[i] = AIPSObject[i]; AOname[i] = 0;

  /* Get AIPS object number */
  ObitAIPSObjectOBname (AIPSObject, &objnum, err);
  if (err->error) Obit_traceback_msg (err, routine, AOname);
  
 /* Read through Info list to find largest and create scratch array */
  size = sizeof(odouble);
  tstMem = inMem;
  while (*tstMem) {
    /* Read from Obit */
    if (ObitInfoListInfo (info, *tstMem, &OType, Odim, err)) {
      if (err->error) Obit_traceback_msg (err, routine, AOname);
      
      /* How big? */
      need = ObitInfoElemSize (OType, Odim);
      size = MAX (size, need);

     } /* end if exist */
    tstMem++;
  } /* end loop over input members */
      
  /* allocate scratch array for copying */
  value = g_malloc0(size);
  Lvalue = (oint*)value;    /* Fortran logical equivalent */
	     
  /* Copy member values to list */
  while (*inMem) {
    /* Read from Obit */
    if (ObitInfoListInfo (info, *inMem, &OType, Odim, err)) {
      ObitInfoListGet (info, *inMem, &OType, Odim, (gpointer)value, err);
      if (err->error) {
	/* Cleanup */
	g_free(value);
	Obit_traceback_msg (err, routine, AOname);
      }

      /* Convert to AIPSish */
      if ((OType == OBIT_oint) || (OType == OBIT_long)) {
	AType = OBIT_AIPSObjectInt;
      } else if (OType == OBIT_float) {
	AType = OBIT_AIPSObjectRe;
      } else if (OType == OBIT_double) {
	AType = OBIT_AIPSObjectDP;
      } else if (OType == OBIT_string) {
	AType = OBIT_AIPSObjectCar;
      } else if (OType == OBIT_bool) {
	AType = OBIT_AIPSObjectLog;
	/* Convert single boolean */
	Lvalue[0] = ObitAIPSBooleanC2F((gboolean)Lvalue[0]);
      }
      for (i=0; i<5; i++) Adim[i] = Odim[i];
      
      /* Write to AIPS */
      ObitAIPSObjectOBput (objnum, *outMem, AType, Adim, (gpointer)value, value, err);
      if (err->error) {
	/* Cleanup */
	g_free(value);
	Obit_traceback_msg (err, routine, AOname);
      }
    } /* End copy keyword if exists */

    inMem++;
    outMem++;
  } /* end loop over input members */

  /* Cleanup */
  g_free(value);

} /* end ObitAIPSFortranCopyInfo */
