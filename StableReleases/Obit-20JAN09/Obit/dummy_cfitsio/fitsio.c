/*  DUMMY Version  All routines Stubbed  */
/*  The FITSIO software was written by William Pence at the High Energy    */
/*  Astrophysic Science Archive Research Center (HEASARC) at the NASA      */
/*  Goddard Space Flight Center.                                           */
/*

Copyright (Unpublished--all rights reserved under the copyright laws of
the United States), U.S. Government as represented by the Administrator
of the National Aeronautics and Space Administration.  No copyright is
claimed in the United States under Title 17, U.S. Code.

Permission to freely use, copy, modify, and distribute this software
and its documentation without fee is hereby granted, provided that this
copyright notice and disclaimer of warranty appears in all copies.

DISCLAIMER:

THE SOFTWARE IS PROVIDED 'AS IS' WITHOUT ANY WARRANTY OF ANY KIND,
EITHER EXPRESSED, IMPLIED, OR STATUTORY, INCLUDING, BUT NOT LIMITED TO,
ANY WARRANTY THAT THE SOFTWARE WILL CONFORM TO SPECIFICATIONS, ANY
IMPLIED WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR
PURPOSE, AND FREEDOM FROM INFRINGEMENT, AND ANY WARRANTY THAT THE
DOCUMENTATION WILL CONFORM TO THE SOFTWARE, OR ANY WARRANTY THAT THE
SOFTWARE WILL BE ERROR FREE.  IN NO EVENT SHALL NASA BE LIABLE FOR ANY
DAMAGES, INCLUDING, BUT NOT LIMITED TO, DIRECT, INDIRECT, SPECIAL OR
CONSEQUENTIAL DAMAGES, ARISING OUT OF, RESULTING FROM, OR IN ANY WAY
CONNECTED WITH THIS SOFTWARE, WHETHER OR NOT BASED UPON WARRANTY,
CONTRACT, TORT , OR OTHERWISE, WHETHER OR NOT INJURY WAS SUSTAINED BY
PERSONS OR PROPERTY OR OTHERWISE, AND WHETHER OR NOT LOSS WAS SUSTAINED
FROM, OR AROSE OUT OF THE RESULTS OF, OR USE OF, THE SOFTWARE OR
SERVICES PROVIDED HEREUNDER."

*/

#include "fitsio.h"
#include "glib.h"

/*----------------  FITS file URL parsing routines -------------*/
int fits_get_token(char **ptr, char *delimiter, char *token, int *isanumber)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
char *fits_split_names(char *list)
  {
  g_error ("cfitsio not implemented");
  return NULL;
  }
int ffiurl(char *url,  char *urltype, char *infile,
                    char *outfile, char *extspec, char *rowfilter,
                    char *binspec, char *colspec, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffrtnm(char *url, char *rootname, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffexts(char *extspec, int *extnum,  char *extname, int *extvers,
          int *hdutype, char *colname, char *rowexpress, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffextn(char *url, int *extension_num, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffurlt(fitsfile *fptr, char *urlType, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffbins(char *binspec, int *imagetype, int *haxis, 
                      char colname[4][FLEN_VALUE], double *minin,
                      double *maxin, double *binsizein,
                      char minname[4][FLEN_VALUE], char maxname[4][FLEN_VALUE],
                      char binname[4][FLEN_VALUE], double *weight, char *wtname,
                      int *recip, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffbinr(char **binspec, char *colname, double *minin, 
                        double *maxin, double *binsizein, char *minname,
                        char *maxname, char *binname, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffimport_file( char *filename, char **contents, int *status )
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffrwrg( char *rowlist, long maxrows, int maxranges, int *numranges,
      long *minrow, long *maxrow, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }

/*----------------  FITS file I/O routines -------------*/
int ffomem(fitsfile **fptr, const char *name, int mode, void **buffptr,
           size_t *buffsize, size_t deltasize,
           void *(*mem_realloc)(void *p, size_t newsize),
           int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffopen(fitsfile **fptr, const char *filename, int iomode, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffdopn(fitsfile **fptr, const char *filename, int iomode, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int fftopn(fitsfile **fptr, const char *filename, int iomode, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffiopn(fitsfile **fptr, const char *filename, int iomode, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffreopen(fitsfile *openfptr, fitsfile **newfptr, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  } 
int ffinit(fitsfile **fptr, const char *filename, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffimem(fitsfile **fptr,  void **buffptr,
           size_t *buffsize, size_t deltasize,
           void *(*mem_realloc)(void *p, size_t newsize),
           int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int fftplt(fitsfile **fptr, const char *filename, const char *tempname,
           int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffflus(fitsfile *fptr, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffflsh(fitsfile *fptr, int clearbuf, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffclos(fitsfile *fptr, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffdelt(fitsfile *fptr, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffflnm(fitsfile *fptr, char *filename, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffflmd(fitsfile *fptr, int *filemode, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }

/*---------------- utility routines -------------*/
float ffvers(float *version)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
void ffupch(char *string)
  {
  g_error ("cfitsio not implemented");
  }
void ffgerr(int status, char *errtext)
  {
  g_error ("cfitsio not implemented");
  }
void ffpmsg(const char *err_message)
  {
  g_error ("cfitsio not implemented");
  }
void ffpmrk(void)
  {
  g_error ("cfitsio not implemented");
  }
int  ffgmsg(char *err_message)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
void ffcmsg(void)
  {
  g_error ("cfitsio not implemented");
  }
void ffcmrk(void)
  {
  g_error ("cfitsio not implemented");
  }
void ffrprt(FILE *stream, int status)
  {
  g_error ("cfitsio not implemented");
  }
void ffcmps(char *templt, char *colname, int  casesen, int *match,
           int *exact)
  {
  g_error ("cfitsio not implemented");
  }
int fftkey(char *keyword, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int fftrec(char *card, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffnchk(fitsfile *fptr, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffkeyn(char *keyroot, int value, char *keyname, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffnkey(int value, char *keyroot, char *keyname, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgkcl(char *card)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffdtyp(char *cval, char *dtype, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpsvc(char *card, char *value, char *comm, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgknm(char *card, char *name, int *length, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgthd(char *tmplt, char *card, int *hdtype, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffasfm(char *tform, int *datacode, long *width, int *decim, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffbnfm(char *tform, int *datacode, long *repeat, long *width, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgabc(int tfields, char **tform, int space, long *rowlen, long *tbcol,
           int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int fits_get_section_range(char **ptr,long *secmin,long *secmax,long *incre,
              int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
 
/*----------------- write single keywords --------------*/
int ffpky(fitsfile *fptr, int datatype, char *keyname, void *value,
          char *comm, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffprec(fitsfile *fptr, const char *card, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpcom(fitsfile *fptr, const char *comm, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpunt(fitsfile *fptr, char *keyname, char *unit, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffphis(fitsfile *fptr, const char *history, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpdat(fitsfile *fptr, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgstm(char *timestr, int *timeref, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgsdt(int *day, int *month, int *year, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffdt2s(int year, int month, int day, char *datestr, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int fftm2s(int year, int month, int day, int hour, int minute, double second,
          int decimals, char *datestr, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffs2dt(char *datestr, int *year, int *month, int *day, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffs2tm(char *datestr, int *year, int *month, int *day, int *hour,
          int *minute, double *second, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpkyu(fitsfile *fptr, char *keyname, char *comm, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpkys(fitsfile *fptr, char *keyname, char *value, char *comm,int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpkls(fitsfile *fptr, char *keyname, char *value, char *comm,int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffplsw(fitsfile *fptr, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpkyl(fitsfile *fptr, char *keyname, int  value, char *comm, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpkyj(fitsfile *fptr, char *keyname, long value, char *comm, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpkyf(fitsfile *fptr, char *keyname, float value, int decim, char *comm,
          int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpkye(fitsfile *fptr, char *keyname, float  value, int decim, char *comm,
          int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpkyg(fitsfile *fptr, char *keyname, double value, int decim, char *comm,
          int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpkyd(fitsfile *fptr, char *keyname, double value, int decim, char *comm,
          int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpkyc(fitsfile *fptr, char *keyname, float *value, int decim, char *comm,
          int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpkym(fitsfile *fptr, char *keyname, double *value, int decim, char *comm,
          int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpkfc(fitsfile *fptr, char *keyname, float *value, int decim, char *comm,
          int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpkfm(fitsfile *fptr, char *keyname, double *value, int decim, char *comm,
          int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpkyt(fitsfile *fptr, char *keyname, long intval, double frac, char *comm,
          int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffptdm( fitsfile *fptr, int colnum, int naxis, long naxes[], int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }

/*----------------- write array of keywords --------------*/
int ffpkns(fitsfile *fptr, char *keyroot, int nstart, int nkey, char *value[],
           char *comm[], int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpknl(fitsfile *fptr, char *keyroot, int nstart, int nkey, int *value,
           char *comm[], int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpknj(fitsfile *fptr, char *keyroot, int nstart, int nkey, long *value,
           char *comm[], int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpknf(fitsfile *fptr, char *keyroot, int nstart, int nkey, float *value,
           int decim, char *comm[], int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpkne(fitsfile *fptr, char *keyroot, int nstart, int nkey, float *value,
           int decim, char *comm[], int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpkng(fitsfile *fptr, char *keyroot, int nstart, int nkey, double *value,
           int decim, char *comm[], int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpknd(fitsfile *fptr, char *keyroot, int nstart, int nkey, double *value,
           int decim, char *comm[], int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffcpky(fitsfile *infptr,fitsfile *outfptr,int incol,int outcol,
           char *rootname, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  } 

/*----------------- write required header keywords --------------*/
int ffphps( fitsfile *fptr, int bitpix, int naxis, long naxes[], int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffphpr( fitsfile *fptr, int simple, int bitpix, int naxis, long naxes[],
            long pcount, long gcount, int extend, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffphtb(fitsfile *fptr, long naxis1, long naxis2, int tfields, char **ttype,
          long *tbcol, char **tform, char **tunit, char *extname, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffphbn(fitsfile *fptr, long naxis2, int tfields, char **ttype,
          char **tform, char **tunit, char *extname, long pcount, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }

/*----------------- write template keywords --------------*/
int ffpktp(fitsfile *fptr, const char *filename, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }

/*------------------ get header information --------------*/
int ffghsp(fitsfile *fptr, int *nexist, int *nmore, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffghps(fitsfile *fptr, int *nexist, int *position, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
 
/*------------------ move position in header -------------*/
int ffmaky(fitsfile *fptr, int nrec, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffmrky(fitsfile *fptr, int nrec, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
 
/*------------------ read single keywords -----------------*/
int ffgnxk(fitsfile *fptr, char **inclist, int ninc, char **exclist,
           int nexc, char *card, int  *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgrec(fitsfile *fptr, int nrec,      char *card, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgcrd(fitsfile *fptr, char *keyname, char *card, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgunt(fitsfile *fptr, char *keyname, char *unit, int  *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgkyn(fitsfile *fptr, int nkey, char *keyname, char *keyval, char *comm,
           int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgkey(fitsfile *fptr, char *keyname, char *keyval, char *comm,
           int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
 
int ffgky( fitsfile *fptr, int datatype, char *keyname, void *value,
           char *comm, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgkys(fitsfile *fptr, char *keyname, char *value, char *comm, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgkls(fitsfile *fptr, char *keyname, char **value, char *comm, int *status)

  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgkyl(fitsfile *fptr, char *keyname, int *value, char *comm, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgkyj(fitsfile *fptr, char *keyname, long *value, char *comm, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgkye(fitsfile *fptr, char *keyname, float *value, char *comm,int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgkyd(fitsfile *fptr, char *keyname, double *value,char *comm,int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgkyc(fitsfile *fptr, char *keyname, float *value, char *comm,int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgkym(fitsfile *fptr, char *keyname, double *value,char *comm,int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgkyt(fitsfile *fptr, char *keyname, long *ivalue, double *dvalue,
           char *comm, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgtdm(fitsfile *fptr, int colnum, int maxdim, int *naxis, long naxes[],
           int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffdtdm(fitsfile *fptr, char *tdimstr, int colnum, int maxdim,
           int *naxis, long naxes[], int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }

/*------------------ read array of keywords -----------------*/
int ffgkns(fitsfile *fptr, char *keyname, int nstart, int nmax, char *value[],
           int *nfound,  int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgknl(fitsfile *fptr, char *keyname, int nstart, int nmax, int *value,
           int *nfound, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgknj(fitsfile *fptr, char *keyname, int nstart, int nmax, long *value,
           int *nfound, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgkne(fitsfile *fptr, char *keyname, int nstart, int nmax, float *value,
           int *nfound, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgknd(fitsfile *fptr, char *keyname, int nstart, int nmax, double *value,
           int *nfound, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffh2st(fitsfile *fptr, char **header, int  *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffhdr2str( fitsfile *fptr,  int exclude_comm, char **exclist,
   int nexc, char **header, int *nkeys, int  *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }

/*----------------- read required header keywords --------------*/
int ffghpr(fitsfile *fptr, int maxdim, int *simple, int *bitpix, int *naxis,
          long naxes[], long *pcount, long *gcount, int *extend, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
 
int ffghtb(fitsfile *fptr,int maxfield, long *naxis1, long *naxis2,
           int *tfields, char **ttype, long *tbcol, char **tform, char **tunit,
           char *extname,  int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
 
int ffghbn(fitsfile *fptr, int maxfield, long *naxis2, int *tfields,
           char **ttype, char **tform, char **tunit, char *extname,
           long *pcount, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
 
/*--------------------- update keywords ---------------*/
int ffuky(fitsfile *fptr, int datatype, char *keyname, void *value,
          char *comm, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffucrd(fitsfile *fptr, char *keyname, char *card, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffukyu(fitsfile *fptr, char *keyname, char *comm, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffukys(fitsfile *fptr, char *keyname, char *value, char *comm, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffukls(fitsfile *fptr, char *keyname, char *value, char *comm, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffukyl(fitsfile *fptr, char *keyname, int value, char *comm, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffukyj(fitsfile *fptr, char *keyname, long value, char *comm, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffukyf(fitsfile *fptr, char *keyname, float value, int decim, char *comm,
          int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffukye(fitsfile *fptr, char *keyname, float value, int decim, char *comm,
          int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffukyg(fitsfile *fptr, char *keyname, double value, int decim, char *comm,
          int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffukyd(fitsfile *fptr, char *keyname, double value, int decim, char *comm,
          int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffukyc(fitsfile *fptr, char *keyname, float *value, int decim, char *comm,
          int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffukym(fitsfile *fptr, char *keyname, double *value, int decim, char *comm,
          int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffukfc(fitsfile *fptr, char *keyname, float *value, int decim, char *comm,
          int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffukfm(fitsfile *fptr, char *keyname, double *value, int decim, char *comm,
          int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }

/*--------------------- modify keywords ---------------*/
int ffmrec(fitsfile *fptr, int nkey, char *card, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffmcrd(fitsfile *fptr, char *keyname, char *card, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffmnam(fitsfile *fptr, char *oldname, char *newname, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffmcom(fitsfile *fptr, char *keyname, char *comm, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffmkyu(fitsfile *fptr, char *keyname, char *comm, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffmkys(fitsfile *fptr, char *keyname, char *value, char *comm,int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffmkls(fitsfile *fptr, char *keyname, char *value, char *comm,int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffmkyl(fitsfile *fptr, char *keyname, int value, char *comm, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffmkyj(fitsfile *fptr, char *keyname, long value, char *comm, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffmkyf(fitsfile *fptr, char *keyname, float value, int decim, char *comm,
          int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffmkye(fitsfile *fptr, char *keyname, float value, int decim, char *comm,
          int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffmkyg(fitsfile *fptr, char *keyname, double value, int decim, char *comm,
          int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffmkyd(fitsfile *fptr, char *keyname, double value, int decim, char *comm,
          int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffmkyc(fitsfile *fptr, char *keyname, float *value, int decim, char *comm,
          int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffmkym(fitsfile *fptr, char *keyname, double *value, int decim, char *comm,
          int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffmkfc(fitsfile *fptr, char *keyname, float *value, int decim, char *comm,
          int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffmkfm(fitsfile *fptr, char *keyname, double *value, int decim, char *comm,
          int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
 
/*--------------------- insert keywords ---------------*/
int ffirec(fitsfile *fptr, int nkey, char *card, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffikey(fitsfile *fptr, char *card, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffikyu(fitsfile *fptr, char *keyname, char *comm, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffikys(fitsfile *fptr, char *keyname, char *value, char *comm,int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffikls(fitsfile *fptr, char *keyname, char *value, char *comm,int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffikyl(fitsfile *fptr, char *keyname, int value, char *comm, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffikyj(fitsfile *fptr, char *keyname, long value, char *comm, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffikyf(fitsfile *fptr, char *keyname, float value, int decim, char *comm,
          int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffikye(fitsfile *fptr, char *keyname, float value, int decim, char *comm,
          int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffikyg(fitsfile *fptr, char *keyname, double value, int decim, char *comm,
          int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffikyd(fitsfile *fptr, char *keyname, double value, int decim, char *comm,
          int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffikyc(fitsfile *fptr, char *keyname, float *value, int decim, char *comm,
          int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffikym(fitsfile *fptr, char *keyname, double *value, int decim, char *comm,
          int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffikfc(fitsfile *fptr, char *keyname, float *value, int decim, char *comm,
          int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffikfm(fitsfile *fptr, char *keyname, double *value, int decim, char *comm,
          int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }

/*--------------------- delete keywords ---------------*/
int ffdkey(fitsfile *fptr, char *keyname, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffdrec(fitsfile *fptr, int keypos, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
 
/*--------------------- get HDU information -------------*/
int ffghdn(fitsfile *fptr, int *chdunum)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffghdt(fitsfile *fptr, int *exttype, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffghad(fitsfile *fptr, long *headstart, long *datastart, long *dataend,
           int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffghof(fitsfile *fptr, OFF_T *headstart, OFF_T *datastart, OFF_T *dataend,
           int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgipr(fitsfile *fptr, int maxaxis, int *imgtype, int *naxis,
           long *naxes, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgidt(fitsfile *fptr, int *imgtype, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgiet(fitsfile *fptr, int *imgtype, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgidm(fitsfile *fptr, int *naxis,  int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgisz(fitsfile *fptr, int nlen, long *naxes, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }

/*--------------------- HDU operations -------------*/
int ffmahd(fitsfile *fptr, int hdunum, int *exttype, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffmrhd(fitsfile *fptr, int hdumov, int *exttype, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffmnhd(fitsfile *fptr, int exttype, char *hduname, int hduvers,
           int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffthdu(fitsfile *fptr, int *nhdu, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffcrhd(fitsfile *fptr, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffcrim(fitsfile *fptr, int bitpix, int naxis, long *naxes, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffcrtb(fitsfile *fptr, int tbltype, long naxis2, int tfields, char **ttype,
           char **tform, char **tunit, char *extname, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffiimg(fitsfile *fptr, int bitpix, int naxis, long *naxes, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffitab(fitsfile *fptr, long naxis1, long naxis2, int tfields, char **ttype,
           long *tbcol, char **tform, char **tunit, char *extname, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffibin(fitsfile *fptr,long naxis2, int tfields, char **ttype, char **tform,
           char **tunit, char *extname, long pcount, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffrsim(fitsfile *fptr, int bitpix, int naxis, long *naxes, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffdhdu(fitsfile *fptr, int *hdutype, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffcopy(fitsfile *infptr, fitsfile *outfptr, int morekeys, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffcpfl(fitsfile *infptr, fitsfile *outfptr, int prev, int cur, int follow,
            int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffcphd(fitsfile *infptr, fitsfile *outfptr, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffcpdt(fitsfile *infptr, fitsfile *outfptr, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffchfl(fitsfile *fptr, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffcdfl(fitsfile *fptr, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }

int ffrdef(fitsfile *fptr, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffhdef(fitsfile *fptr, int morekeys, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpthp(fitsfile *fptr, long theap, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
 
int ffcsum(fitsfile *fptr, long nrec, unsigned long *sum, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
void ffesum(unsigned long sum, int complm, char *ascii)
  {
  g_error ("cfitsio not implemented");
  }
unsigned long ffdsum(char *ascii, int complm, unsigned long *sum)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpcks(fitsfile *fptr, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffupck(fitsfile *fptr, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffvcks(fitsfile *fptr, int *datastatus, int *hdustatus, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgcks(fitsfile *fptr, unsigned long *datasum, unsigned long *hdusum,
    int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
 
/*--------------------- define scaling or null values -------------*/
int ffpscl(fitsfile *fptr, double scale, double zero, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpnul(fitsfile *fptr, long nulvalue, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int fftscl(fitsfile *fptr, int colnum, double scale, double zero, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int fftnul(fitsfile *fptr, int colnum, long nulvalue, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffsnul(fitsfile *fptr, int colnum, char *nulstring, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
 
/*--------------------- get column information -------------*/
int ffgcno(fitsfile *fptr, int casesen, char *templt, int  *colnum,
           int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgcnn(fitsfile *fptr, int casesen, char *templt, char *colname,
           int *colnum, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
 
int ffgtcl(fitsfile *fptr, int colnum, int *typecode, long *repeat,
           long *width, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffeqty(fitsfile *fptr, int colnum, int *typecode, long *repeat,
           long *width, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgncl(fitsfile *fptr, int  *ncols, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgnrw(fitsfile *fptr, long *nrows, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgacl(fitsfile *fptr, int colnum, char *ttype, long *tbcol,
           char *tunit, char *tform, double *tscal, double *tzero,
           char *tnull, char *tdisp, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgbcl(fitsfile *fptr, int colnum, char *ttype, char *tunit,
           char *dtype, long *repeat, double *tscal, double *tzero,
           long *tnull, char *tdisp, int  *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgrsz(fitsfile *fptr, long *nrows, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgcdw(fitsfile *fptr, int colnum, int *width, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }

/*--------------------- read primary array or image elements -------------*/
int ffgpxv(fitsfile *fptr, int  datatype, long *firstpix, long nelem,
          void *nulval, void *array, int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgpxf(fitsfile *fptr, int  datatype, long *firstpix, long nelem,
           void *array, char *nullarray, int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgsv(fitsfile *fptr, int datatype, long *blc, long *trc, long *inc,
          void *nulval, void *array, int *anynul, int  *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgpv(fitsfile *fptr, int  datatype, long firstelem, long nelem,
          void *nulval, void *array, int *anynul, int  *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgpf(fitsfile *fptr, int  datatype, long firstelem, long nelem,
          void *array, char *nullarray, int  *anynul, int  *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgpvb(fitsfile *fptr, long group, long firstelem, long nelem, unsigned
           char nulval, unsigned char *array, int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgpvsb(fitsfile *fptr, long group, long firstelem, long nelem, signed
           char nulval, signed char *array, int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgpvui(fitsfile *fptr, long group, long firstelem, long nelem,
           unsigned short nulval, unsigned short *array, int *anynul, 
           int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgpvi(fitsfile *fptr, long group, long firstelem, long nelem,
           short nulval, short *array, int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgpvuj(fitsfile *fptr, long group, long firstelem, long nelem,
           unsigned long nulval, unsigned long *array, int *anynul, 
           int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgpvj(fitsfile *fptr, long group, long firstelem, long nelem,
           long nulval, long *array, int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgpvjj(fitsfile *fptr, long group, long firstelem, long nelem,
           LONGLONG nulval, LONGLONG *array, int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgpvuk(fitsfile *fptr, long group, long firstelem, long nelem,
           unsigned int nulval, unsigned int *array, int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgpvk(fitsfile *fptr, long group, long firstelem, long nelem,
           int nulval, int *array, int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgpve(fitsfile *fptr, long group, long firstelem, long nelem,
           float nulval, float *array, int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgpvd(fitsfile *fptr, long group, long firstelem, long nelem,
           double nulval, double *array, int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
 
int ffgpfb(fitsfile *fptr, long group, long firstelem, long nelem,
           unsigned char *array, char *nularray, int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgpfsb(fitsfile *fptr, long group, long firstelem, long nelem,
           signed char *array, char *nularray, int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgpfui(fitsfile *fptr, long group, long firstelem, long nelem,
           unsigned short *array, char *nularray, int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgpfi(fitsfile *fptr, long group, long firstelem, long nelem,
           short *array, char *nularray, int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgpfuj(fitsfile *fptr, long group, long firstelem, long nelem,
           unsigned long *array, char *nularray, int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgpfj(fitsfile *fptr, long group, long firstelem, long nelem,
           long *array, char *nularray, int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgpfjj(fitsfile *fptr, long group, long firstelem, long nelem,
           LONGLONG *array, char *nularray, int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgpfuk(fitsfile *fptr, long group, long firstelem, long nelem,
           unsigned int *array, char *nularray, int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgpfk(fitsfile *fptr, long group, long firstelem, long nelem,
           int *array, char *nularray, int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgpfe(fitsfile *fptr, long group, long firstelem, long nelem,
           float *array, char *nularray, int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgpfd(fitsfile *fptr, long group, long firstelem, long nelem,
           double *array, char *nularray, int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
 
int ffg2db(fitsfile *fptr, long group, unsigned char nulval, long ncols,
           long naxis1, long naxis2, unsigned char *array,
           int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffg2dsb(fitsfile *fptr, long group, signed char nulval, long ncols,
           long naxis1, long naxis2, signed char *array,
           int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffg2dui(fitsfile *fptr, long group, unsigned short nulval, long ncols,
           long naxis1, long naxis2, unsigned short *array,
           int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffg2di(fitsfile *fptr, long group, short nulval, long ncols,
           long naxis1, long naxis2, short *array,
           int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffg2duj(fitsfile *fptr, long group, unsigned long nulval, long ncols,
           long naxis1, long naxis2, unsigned long *array,
           int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffg2dj(fitsfile *fptr, long group, long nulval, long ncols,
           long naxis1, long naxis2, long *array,
           int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffg2djj(fitsfile *fptr, long group, LONGLONG nulval, long ncols,
           long naxis1, long naxis2, LONGLONG *array,
           int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffg2duk(fitsfile *fptr, long group, unsigned int nulval, long ncols,
           long naxis1, long naxis2, unsigned int *array,
           int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffg2dk(fitsfile *fptr, long group, int nulval, long ncols,
           long naxis1, long naxis2, int *array,
           int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffg2de(fitsfile *fptr, long group, float nulval, long ncols,
           long naxis1, long naxis2, float *array,
           int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffg2dd(fitsfile *fptr, long group, double nulval, long ncols,
           long naxis1, long naxis2, double *array,
           int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
 
int ffg3db(fitsfile *fptr, long group, unsigned char nulval, long ncols,
           long nrows, long naxis1, long naxis2, long naxis3,
           unsigned char *array, int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffg3dsb(fitsfile *fptr, long group, signed char nulval, long ncols,
           long nrows, long naxis1, long naxis2, long naxis3,
           signed char *array, int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffg3dui(fitsfile *fptr, long group, unsigned short nulval, long ncols,
           long nrows, long naxis1, long naxis2, long naxis3,
           unsigned short *array, int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffg3di(fitsfile *fptr, long group, short nulval, long ncols,
           long nrows, long naxis1, long naxis2, long naxis3,
           short *array, int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffg3duj(fitsfile *fptr, long group, unsigned long nulval, long ncols,
           long nrows, long naxis1, long naxis2, long naxis3,
           unsigned long *array, int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffg3dj(fitsfile *fptr, long group, long nulval, long ncols,
           long nrows, long naxis1, long naxis2, long naxis3,
           long *array, int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffg3djj(fitsfile *fptr, long group, LONGLONG nulval, long ncols,
           long nrows, long naxis1, long naxis2, long naxis3,
           LONGLONG *array, int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffg3duk(fitsfile *fptr, long group, unsigned int nulval, long ncols,
           long nrows, long naxis1, long naxis2, long naxis3,
           unsigned int *array, int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffg3dk(fitsfile *fptr, long group, int nulval, long ncols,
           long nrows, long naxis1, long naxis2, long naxis3,
           int *array, int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffg3de(fitsfile *fptr, long group, float nulval, long ncols,
           long nrows, long naxis1, long naxis2, long naxis3,
           float *array, int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffg3dd(fitsfile *fptr, long group, double nulval, long ncols,
           long nrows, long naxis1, long naxis2, long naxis3,
           double *array, int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
 
int ffgsvb(fitsfile *fptr, int colnum, int naxis, long *naxes, long *blc,
  long *trc, long *inc, unsigned char nulval, unsigned char *array,
  int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgsvsb(fitsfile *fptr, int colnum, int naxis, long *naxes, long *blc,
  long *trc, long *inc, signed char nulval, signed char *array,
  int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgsvui(fitsfile *fptr, int colnum, int naxis, long *naxes, long *blc,
  long *trc, long *inc, unsigned short nulval, unsigned short *array, 
  int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgsvi(fitsfile *fptr, int colnum, int naxis, long *naxes, long *blc,
  long *trc, long *inc, short nulval, short *array, int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgsvuj(fitsfile *fptr, int colnum, int naxis, long *naxes, long *blc,
  long *trc, long *inc, unsigned long nulval, unsigned long *array, 
  int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgsvj(fitsfile *fptr, int colnum, int naxis, long *naxes, long *blc,
  long *trc, long *inc, long nulval, long *array, int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgsvjj(fitsfile *fptr, int colnum, int naxis, long *naxes, long *blc,
  long *trc, long *inc, LONGLONG nulval, LONGLONG *array, int *anynul,
  int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgsvuk(fitsfile *fptr, int colnum, int naxis, long *naxes, long *blc,
  long *trc, long *inc, unsigned int nulval, unsigned int *array,
  int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgsvk(fitsfile *fptr, int colnum, int naxis, long *naxes, long *blc,
  long *trc, long *inc, int nulval, int *array, int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgsve(fitsfile *fptr, int colnum, int naxis, long *naxes, long *blc,
  long *trc, long *inc, float nulval, float *array, int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgsvd(fitsfile *fptr, int colnum, int naxis, long *naxes, long *blc,
  long *trc, long *inc, double nulval, double *array, int *anynul,
  int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
 
int ffgsfb(fitsfile *fptr, int colnum, int naxis, long *naxes, long *blc,
  long *trc, long *inc, unsigned char *array, char *flagval,
  int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgsfsb(fitsfile *fptr, int colnum, int naxis, long *naxes, long *blc,
  long *trc, long *inc, signed char *array, char *flagval,
  int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgsfui(fitsfile *fptr, int colnum, int naxis, long *naxes, long *blc,
  long *trc, long *inc, unsigned short *array, char *flagval, int *anynul, 
  int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgsfi(fitsfile *fptr, int colnum, int naxis, long *naxes, long *blc,
  long *trc, long *inc, short *array, char *flagval, int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgsfuj(fitsfile *fptr, int colnum, int naxis, long *naxes, long *blc,
  long  *trc, long *inc, unsigned long *array, char *flagval, int *anynul,
  int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgsfj(fitsfile *fptr, int colnum, int naxis, long *naxes, long *blc,
  long  *trc, long *inc, long *array, char *flagval, int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgsfjj(fitsfile *fptr, int colnum, int naxis, long *naxes, long *blc,
  long  *trc, long *inc, LONGLONG *array, char *flagval, int *anynul,
  int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgsfuk(fitsfile *fptr, int colnum, int naxis, long *naxes, long *blc,
  long  *trc, long *inc, unsigned int *array, char *flagval, int *anynul,
  int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgsfk(fitsfile *fptr, int colnum, int naxis, long *naxes, long *blc,
  long  *trc, long *inc, int *array, char *flagval, int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgsfe(fitsfile *fptr, int colnum, int naxis, long *naxes, long *blc,
  long *trc, long *inc, float *array, char *flagval, int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgsfd(fitsfile *fptr, int colnum, int naxis, long *naxes, long *blc,
  long *trc, long *inc, double *array, char *flagval, int *anynul,
  int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
 
int ffggpb(fitsfile *fptr, long group, long firstelem, long nelem,
           unsigned char *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffggpsb(fitsfile *fptr, long group, long firstelem, long nelem,
           signed char *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffggpui(fitsfile *fptr, long group, long firstelem, long nelem,
           unsigned short *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffggpi(fitsfile *fptr, long group, long firstelem, long nelem,
           short *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffggpuj(fitsfile *fptr, long group, long firstelem, long nelem,
           unsigned long *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffggpj(fitsfile *fptr, long group, long firstelem, long nelem,
           long *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffggpjj(fitsfile *fptr, long group, long firstelem, long nelem,
           LONGLONG *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffggpuk(fitsfile *fptr, long group, long firstelem, long nelem,
           unsigned int *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffggpk(fitsfile *fptr, long group, long firstelem, long nelem,
           int *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffggpe(fitsfile *fptr, long group, long firstelem, long nelem,
           float *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffggpd(fitsfile *fptr, long group, long firstelem, long nelem,
           double *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
 
/*--------------------- read column elements -------------*/
int ffgcv( fitsfile *fptr, int datatype, int colnum, long firstrow,
           long firstelem, long nelem, void *nulval, void *array, int *anynul,
           int  *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgcf( fitsfile *fptr, int datatype, int colnum, long firstrow,
           long firstelem, long nelem, void *array, char *nullarray,
           int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgcvs(fitsfile *fptr, int colnum, long firstrow, long firstelem,
           long nelem, char *nulval, char **array, int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgcl (fitsfile *fptr, int colnum, long firstrow, long firstelem,
           long nelem, char *array, int  *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgcvl (fitsfile *fptr, int colnum, long firstrow, long firstelem,
           long nelem, char nulval, char *array, int *anynul, int  *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgcvb(fitsfile *fptr, int colnum, long firstrow, long firstelem,
           long nelem, unsigned char nulval, unsigned char *array,
           int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgcvsb(fitsfile *fptr, int colnum, long firstrow, long firstelem,
           long nelem, signed char nulval, signed char *array,
           int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgcvui(fitsfile *fptr, int colnum, long firstrow, long firstelem,
           long nelem, unsigned short nulval, unsigned short *array, 
           int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgcvi(fitsfile *fptr, int colnum, long firstrow, long firstelem,
           long nelem, short nulval, short *array, int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgcvuj(fitsfile *fptr, int colnum, long firstrow, long firstelem,
           long nelem, unsigned long nulval, unsigned long *array, int *anynul,
           int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgcvj(fitsfile *fptr, int colnum, long firstrow, long firstelem,
           long nelem, long nulval, long *array, int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgcvjj(fitsfile *fptr, int colnum, long firstrow, long firstelem,
           long nelem, LONGLONG nulval, LONGLONG *array, int *anynul,
           int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgcvuk(fitsfile *fptr, int colnum, long firstrow, long firstelem,
           long nelem, unsigned int nulval, unsigned int *array, int *anynul,
           int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgcvk(fitsfile *fptr, int colnum, long firstrow, long firstelem,
           long nelem, int nulval, int *array, int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgcve(fitsfile *fptr, int colnum, long firstrow, long firstelem,
           long nelem, float nulval, float *array, int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgcvd(fitsfile *fptr, int colnum, long firstrow, long firstelem,
         long nelem, double nulval, double *array, int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgcvc(fitsfile *fptr, int colnum, long firstrow, long firstelem,
           long nelem, float nulval, float *array, int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgcvm(fitsfile *fptr, int colnum, long firstrow, long firstelem,
         long nelem, double nulval, double *array, int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgcx(fitsfile *fptr, int colnum, long firstrow, long firstbit,
            long nbits, char *larray, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgcxui(fitsfile *fptr, int colnum, long firstrow, long nrows,
            long firstbit, int nbits, unsigned short *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgcxuk(fitsfile *fptr, int colnum, long firstrow, long nrows,
            long firstbit, int nbits, unsigned int *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }

int ffgcfs(fitsfile *fptr, int colnum, long firstrow, long firstelem, long
          nelem, char **array, char *nularray, int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgcfl(fitsfile *fptr, int colnum, long firstrow, long firstelem, long
          nelem, char *array, char *nularray, int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgcfb(fitsfile *fptr, int colnum, long firstrow, long firstelem, long
      nelem, unsigned char *array, char *nularray, int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgcfsb(fitsfile *fptr, int colnum, long firstrow, long firstelem, long
      nelem, signed char *array, char *nularray, int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgcfui(fitsfile *fptr, int colnum, long firstrow, long firstelem,
      long nelem, unsigned short *array, char *nularray, int *anynul, 
      int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgcfi(fitsfile *fptr, int colnum, long firstrow, long firstelem,
      long nelem, short *array, char *nularray, int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgcfuj(fitsfile *fptr, int colnum, long firstrow, long firstelem,
      long nelem, unsigned long *array, char *nularray, int *anynul,
      int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgcfj(fitsfile *fptr, int colnum, long firstrow, long firstelem,
      long nelem, long *array, char *nularray, int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgcfjj(fitsfile *fptr, int colnum, long firstrow, long firstelem,
      long nelem, LONGLONG *array, char *nularray, int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgcfuk(fitsfile *fptr, int colnum, long firstrow, long firstelem,
      long nelem, unsigned int *array, char *nularray, int *anynul,
      int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgcfk(fitsfile *fptr, int colnum, long firstrow, long firstelem,
      long nelem, int *array, char *nularray, int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgcfe(fitsfile *fptr, int colnum, long firstrow, long firstelem,
      long nelem, float *array, char *nularray, int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgcfd(fitsfile *fptr, int colnum, long firstrow, long firstelem,
      long nelem, double *array, char *nularray, int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgcfc(fitsfile *fptr, int colnum, long firstrow, long firstelem,
      long nelem, float *array, char *nularray, int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgcfm(fitsfile *fptr, int colnum, long firstrow, long firstelem,
      long nelem, double *array, char *nularray, int *anynul, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
 
int ffgdes(fitsfile *fptr, int colnum, long rownum, long *length,
           long *heapaddr, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
 
int ffgdess(fitsfile *fptr, int colnum, long firstrow, long nrows, long *length,
           long *heapaddr, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
 
int fftheap(fitsfile *fptr, long *heapsize, long *unused, long *overlap,
            int *valid, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffcmph(fitsfile *fptr, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }

int ffgtbb(fitsfile *fptr, long firstrow, long firstchar, long nchars,
           unsigned char *values, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
 
/*------------ write primary array or image elements -------------*/
int ffppx(fitsfile *fptr, int datatype, long  *firstpix, long nelem,
          void *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffppxn(fitsfile *fptr, int datatype, long  *firstpix, long nelem,
          void *array, void *nulval, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffppr(fitsfile *fptr, int datatype, long  firstelem, long nelem,
          void *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpprb(fitsfile *fptr, long group, long firstelem,
           long nelem, unsigned char *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpprsb(fitsfile *fptr, long group, long firstelem,
           long nelem, signed char *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpprui(fitsfile *fptr, long group, long firstelem,
           long nelem, unsigned short *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffppri(fitsfile *fptr, long group, long firstelem,
           long nelem, short *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffppruj(fitsfile *fptr, long group, long firstelem,
           long nelem, unsigned long *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpprj(fitsfile *fptr, long group, long firstelem,
           long nelem, long *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffppruk(fitsfile *fptr, long group, long firstelem,
           long nelem, unsigned int *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpprk(fitsfile *fptr, long group, long firstelem,
           long nelem, int *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffppre(fitsfile *fptr, long group, long firstelem,
           long nelem, float *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpprd(fitsfile *fptr, long group, long firstelem,
           long nelem, double *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpprjj(fitsfile *fptr, long group, long firstelem,
           long nelem, LONGLONG *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }

int ffppru(fitsfile *fptr, long group, long firstelem, long nelem,
           int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpprn(fitsfile *fptr, long firstelem, long nelem, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
 
int ffppn(fitsfile *fptr, int datatype, long  firstelem, long  nelem,
          void  *array, void *nulval, int  *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffppnb(fitsfile *fptr, long group, long firstelem, long nelem,
           unsigned char *array, unsigned char nulval, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffppnsb(fitsfile *fptr, long group, long firstelem, long nelem,
           signed char *array, signed char nulval, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffppnui(fitsfile *fptr, long group, long firstelem,
           long nelem, unsigned short *array, unsigned short nulval,
           int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffppni(fitsfile *fptr, long group, long firstelem,
           long nelem, short *array, short nulval, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffppnj(fitsfile *fptr, long group, long firstelem,
           long nelem, long *array, long nulval, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffppnuj(fitsfile *fptr, long group, long firstelem, long nelem,
           unsigned long *array, unsigned long nulval, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffppnuk(fitsfile *fptr, long group, long firstelem, long nelem,
           unsigned int *array, unsigned int nulval, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffppnk(fitsfile *fptr, long group, long firstelem,
           long nelem, int *array, int nulval, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffppne(fitsfile *fptr, long group, long firstelem,
           long nelem, float *array, float nulval, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffppnd(fitsfile *fptr, long group, long firstelem,
           long nelem, double *array, double nulval, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffppnjj(fitsfile *fptr, long group, long firstelem,
           long nelem, LONGLONG *array, long nulval, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }

int ffp2db(fitsfile *fptr, long group, long ncols, long naxis1,
           long naxis2, unsigned char *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffp2dsb(fitsfile *fptr, long group, long ncols, long naxis1,
           long naxis2, signed char *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffp2dui(fitsfile *fptr, long group, long ncols, long naxis1,
           long naxis2, unsigned short *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffp2di(fitsfile *fptr, long group, long ncols, long naxis1,
           long naxis2, short *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffp2duj(fitsfile *fptr, long group, long ncols, long naxis1,
           long naxis2, unsigned long *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffp2dj(fitsfile *fptr, long group, long ncols, long naxis1,
           long naxis2, long *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffp2duk(fitsfile *fptr, long group, long ncols, long naxis1,
           long naxis2, unsigned int *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffp2dk(fitsfile *fptr, long group, long ncols, long naxis1,
           long naxis2, int *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffp2de(fitsfile *fptr, long group, long ncols, long naxis1,
           long naxis2, float *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffp2dd(fitsfile *fptr, long group, long ncols, long naxis1,
           long naxis2, double *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffp2djj(fitsfile *fptr, long group, long ncols, long naxis1,
           long naxis2, LONGLONG *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }

int ffp3db(fitsfile *fptr, long group, long ncols, long nrows, long naxis1,
           long naxis2, long naxis3, unsigned char *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffp3dsb(fitsfile *fptr, long group, long ncols, long nrows, long naxis1,
           long naxis2, long naxis3, signed char *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffp3dui(fitsfile *fptr, long group, long ncols, long nrows, long naxis1,
           long naxis2, long naxis3, unsigned short *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffp3di(fitsfile *fptr, long group, long ncols, long nrows, long naxis1,
           long naxis2, long naxis3, short *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffp3duj(fitsfile *fptr, long group, long ncols, long nrows, long naxis1,
           long naxis2, long naxis3, unsigned long *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffp3dj(fitsfile *fptr, long group, long ncols, long nrows, long naxis1,
           long naxis2, long naxis3, long *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffp3duk(fitsfile *fptr, long group, long ncols, long nrows, long naxis1,
           long naxis2, long naxis3, unsigned int *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffp3dk(fitsfile *fptr, long group, long ncols, long nrows, long naxis1,
           long naxis2, long naxis3, int *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffp3de(fitsfile *fptr, long group, long ncols, long nrows, long naxis1,
           long naxis2, long naxis3, float *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffp3dd(fitsfile *fptr, long group, long ncols, long nrows, long naxis1,
           long naxis2, long naxis3, double *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffp3djj(fitsfile *fptr, long group, long ncols, long nrows, long naxis1,
           long naxis2, long naxis3, LONGLONG *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }

int ffpss(fitsfile *fptr, int datatype,
           long *fpixel, long *lpixel, void *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpssb(fitsfile *fptr, long group, long naxis, long *naxes,
           long *fpixel, long *lpixel, unsigned char *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpsssb(fitsfile *fptr, long group, long naxis, long *naxes,
           long *fpixel, long *lpixel, signed char *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpssui(fitsfile *fptr, long group, long naxis, long *naxes,
           long *fpixel, long *lpixel, unsigned short *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpssi(fitsfile *fptr, long group, long naxis, long *naxes,
           long *fpixel, long *lpixel, short *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpssuj(fitsfile *fptr, long group, long naxis, long *naxes,
           long *fpixel, long *lpixel, unsigned long *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpssj(fitsfile *fptr, long group, long naxis, long *naxes,
           long *fpixel, long *lpixel, long *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpssuk(fitsfile *fptr, long group, long naxis, long *naxes,
           long *fpixel, long *lpixel, unsigned int *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpssk(fitsfile *fptr, long group, long naxis, long *naxes,
           long *fpixel, long *lpixel, int *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpsse(fitsfile *fptr, long group, long naxis, long *naxes,
           long *fpixel, long *lpixel, float *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpssd(fitsfile *fptr, long group, long naxis, long *naxes,
           long *fpixel, long *lpixel, double *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpssjj(fitsfile *fptr, long group, long naxis, long *naxes,
           long *fpixel, long *lpixel, LONGLONG *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }

int ffpgpb(fitsfile *fptr, long group, long firstelem,
           long nelem, unsigned char *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpgpsb(fitsfile *fptr, long group, long firstelem,
           long nelem, signed char *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpgpui(fitsfile *fptr, long group, long firstelem,
           long nelem, unsigned short *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpgpi(fitsfile *fptr, long group, long firstelem,
           long nelem, short *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpgpuj(fitsfile *fptr, long group, long firstelem,
           long nelem, unsigned long *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpgpj(fitsfile *fptr, long group, long firstelem,
           long nelem, long *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpgpuk(fitsfile *fptr, long group, long firstelem,
           long nelem, unsigned int *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpgpk(fitsfile *fptr, long group, long firstelem,
           long nelem, int *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpgpe(fitsfile *fptr, long group, long firstelem,
           long nelem, float *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpgpd(fitsfile *fptr, long group, long firstelem,
           long nelem, double *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpgpjj(fitsfile *fptr, long group, long firstelem,
           long nelem, LONGLONG *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }

/*--------------------- iterator functions -------------*/
int fits_iter_set_by_name(iteratorCol *col, fitsfile *fptr, char *colname,
          int datatype,  int iotype)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int fits_iter_set_by_num(iteratorCol *col, fitsfile *fptr, int colnum,
          int datatype,  int iotype)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int fits_iter_set_file(iteratorCol *col, fitsfile *fptr)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int fits_iter_set_colname(iteratorCol *col, char *colname)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int fits_iter_set_colnum(iteratorCol *col, int colnum)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int fits_iter_set_datatype(iteratorCol *col, int datatype)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int fits_iter_set_iotype(iteratorCol *col, int iotype)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }

fitsfile * fits_iter_get_file(iteratorCol *col)
  {
  g_error ("cfitsio not implemented");
  return NULL;
  }
char * fits_iter_get_colname(iteratorCol *col)
  {
  g_error ("cfitsio not implemented");
  return NULL;
  }
int fits_iter_get_colnum(iteratorCol *col)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int fits_iter_get_datatype(iteratorCol *col)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int fits_iter_get_iotype(iteratorCol *col)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
void * fits_iter_get_array(iteratorCol *col)
  {
  g_error ("cfitsio not implemented");
  return NULL;
  }
long fits_iter_get_tlmin(iteratorCol *col)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
long fits_iter_get_tlmax(iteratorCol *col)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
long fits_iter_get_repeat(iteratorCol *col)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
char * fits_iter_get_tunit(iteratorCol *col)
  {
  g_error ("cfitsio not implemented");
  return NULL;
  }
char * fits_iter_get_tdisp(iteratorCol *col)
  {
  g_error ("cfitsio not implemented");
  return NULL;
  }

int ffiter(int ncols,  iteratorCol *data, long offset, long nPerLoop,
           int (*workFn)( long totaln, long offset, long firstn,
             long nvalues, int narrays, iteratorCol *data, void *userPointer),
           void *userPointer, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }

/*--------------------- write column elements -------------*/
int ffpcl(fitsfile *fptr, int datatype, int colnum, long firstrow,
          long firstelem, long nelem, void *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpcls(fitsfile *fptr, int colnum, long firstrow, long firstelem,
           long nelem, char **array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpcll(fitsfile *fptr, int colnum, long firstrow, long firstelem,
           long nelem, char *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpclb(fitsfile *fptr, int colnum, long firstrow, long firstelem,
           long nelem, unsigned char *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpclsb(fitsfile *fptr, int colnum, long firstrow, long firstelem,
           long nelem, signed char *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpclui(fitsfile *fptr, int colnum, long firstrow, long firstelem,
           long nelem, unsigned short *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpcli(fitsfile *fptr, int colnum, long firstrow, long firstelem,
           long nelem, short *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpcluj(fitsfile *fptr, int colnum, long firstrow, long firstelem,
           long nelem, unsigned long *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpclj(fitsfile *fptr, int colnum, long firstrow, long firstelem,
           long nelem, long *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpcluk(fitsfile *fptr, int colnum, long firstrow, long firstelem,
           long nelem, unsigned int *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpclk(fitsfile *fptr, int colnum, long firstrow, long firstelem,
           long nelem, int *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpcle(fitsfile *fptr, int colnum, long firstrow, long firstelem,
           long nelem, float *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpcld(fitsfile *fptr, int colnum, long firstrow, long firstelem,
           long nelem, double *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpclc(fitsfile *fptr, int colnum, long firstrow, long firstelem,
           long nelem, float *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpclm(fitsfile *fptr, int colnum, long firstrow, long firstelem,
           long nelem, double *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpclu(fitsfile *fptr, int colnum, long firstrow, long firstelem,
           long nelem, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpclx(fitsfile *fptr, int colnum, long frow, long fbit, long nbit,
            char *larray, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpcljj(fitsfile *fptr, int colnum, long firstrow, long firstelem,
           long nelem, LONGLONG *array, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }

int ffpcn(fitsfile *fptr, int datatype, int colnum, long firstrow,
          long firstelem, long nelem, void *array, void *nulval, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpcns( fitsfile *fptr, int  colnum, long  firstrow, long  firstelem,
            long  nelem, char **array, char  *nulvalue, int  *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpcnl( fitsfile *fptr, int  colnum, long  firstrow, long  firstelem,
            long  nelem, char *array, char  nulvalue,  int  *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpcnb(fitsfile *fptr, int colnum, long firstrow, long firstelem,
           long nelem, unsigned char *array, unsigned char nulvalue,
           int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpcnsb(fitsfile *fptr, int colnum, long firstrow, long firstelem,
           long nelem, signed char *array, signed char nulvalue,
           int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpcnui(fitsfile *fptr, int colnum, long firstrow, long firstelem,
           long nelem, unsigned short *array, unsigned short nulvalue,
           int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpcni(fitsfile *fptr, int colnum, long firstrow, long firstelem,
           long nelem, short *array, short nulvalue, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpcnuj(fitsfile *fptr, int colnum, long firstrow, long firstelem,
           long nelem, unsigned long *array, unsigned long nulvalue,
           int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpcnj(fitsfile *fptr, int colnum, long firstrow, long firstelem,
           long nelem, long *array, long nulvalue, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpcnuk(fitsfile *fptr, int colnum, long firstrow, long firstelem,
           long nelem, unsigned int *array, unsigned int nulvalue,
           int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpcnk(fitsfile *fptr, int colnum, long firstrow, long firstelem,
           long nelem, int *array, int nulvalue, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpcne(fitsfile *fptr, int colnum, long firstrow, long firstelem,
           long nelem, float *array, float nulvalue, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpcnd(fitsfile *fptr, int colnum, long firstrow, long firstelem,
           long nelem, double *array, double nulvalue, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffpcnjj(fitsfile *fptr, int colnum, long firstrow, long firstelem,
           long nelem, LONGLONG *array, LONGLONG nulvalue, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }

int ffpdes(fitsfile *fptr, int colnum, long rownum, long length,
           long heapaddr, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
 
int ffptbb(fitsfile *fptr, long firstrow, long firstchar, long nchars,
           unsigned char *values, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
 
int ffirow(fitsfile *fptr, long firstrow, long nrows, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffdrow(fitsfile *fptr, long firstrow, long nrows, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffdrrg(fitsfile *fptr, char *ranges, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffdrws(fitsfile *fptr, long *rownum,  long nrows, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int fficol(fitsfile *fptr, int numcol, char *ttype, char *tform, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int fficls(fitsfile *fptr, int firstcol, int ncols, char **ttype,
           char **tform, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffmvec(fitsfile *fptr, int colnum, long newveclen, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffdcol(fitsfile *fptr, int numcol, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffcpcl(fitsfile *infptr, fitsfile *outfptr, int incol, int outcol, 
           int create_col, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }

/*--------------------- WCS Utilities ------------------*/
int ffgics(fitsfile *fptr, double *xrval, double *yrval, double *xrpix,
           double *yrpix, double *xinc, double *yinc, double *rot,
           char *type, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgtcs(fitsfile *fptr, int xcol, int ycol, double *xrval,
           double *yrval, double *xrpix, double *yrpix, double *xinc,
           double *yinc, double *rot, char *type, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffwldp(double xpix, double ypix, double xref, double yref,
           double xrefpix, double yrefpix, double xinc, double yinc,
           double rot, char *type, double *xpos, double *ypos, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffxypx(double xpos, double ypos, double xref, double yref, 
           double xrefpix, double yrefpix, double xinc, double yinc,
           double rot, char *type, double *xpix, double *ypix, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }

/*   WCS support routines (provide interface to Doug Mink's WCS library */
int ffgiwcs(fitsfile *fptr,  char **header, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  } 
int ffgtwcs(fitsfile *fptr, int xcol, int ycol, char **header, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }

/*--------------------- lexical parsing routines ------------------*/
int fftexp( fitsfile *fptr, char *expr, int maxdim,
	    int *datatype, long *nelem, int *naxis,
	    long *naxes, int *status )
  {
  g_error ("cfitsio not implemented");
  return 1;
  }

int fffrow( fitsfile *infptr, char *expr,
	    long firstrow, long nrows,
            long *n_good_rows, char *row_status, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }

int ffffrw( fitsfile *fptr, char *expr, long *rownum, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }

int fffrwc( fitsfile *fptr, char *expr, char *timeCol,    
            char *parCol, char *valCol, long ntimes,      
            double *times, char *time_status, int  *status )
  {
  g_error ("cfitsio not implemented");
  return 1;
  }

int ffsrow( fitsfile *infptr, fitsfile *outfptr, char *expr, 
            int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }

int ffcrow( fitsfile *fptr, int datatype, char *expr,
	    long firstrow, long nelements, void *nulval,
	    void *array, int *anynul, int *status )
  {
  g_error ("cfitsio not implemented");
  return 1;
  }

int ffcalc_rng( fitsfile *infptr, char *expr, fitsfile *outfptr,
               char *parName, char *parInfo, int nRngs,
                 long *start, long *end, int *status )
  {
  g_error ("cfitsio not implemented");
  return 1;
  }

int ffcalc( fitsfile *infptr, char *expr, fitsfile *outfptr,
            char *parName, char *parInfo, int *status )
  {
  g_error ("cfitsio not implemented");
  return 1;
  }

  /* ffhist is not really intended as a user-callable routine */
  /* but it may be useful for some specialized applications   */

int ffhist(fitsfile **fptr, char *outfile, int imagetype, int naxis,
           char colname[4][FLEN_VALUE],
           double *minin, double *maxin, double *binsizein,
           char minname[4][FLEN_VALUE], char maxname[4][FLEN_VALUE],
           char binname[4][FLEN_VALUE], 
           double weightin, char wtcol[FLEN_VALUE],
           int recip, char *rowselect, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }

int fits_select_image_section(fitsfile **fptr, char *outfile,
           char *imagesection, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int fits_select_section( fitsfile *infptr, fitsfile *outfptr,
           char *imagesection, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }

/*--------------------- grouping routines ------------------*/

int ffgtcr(fitsfile *fptr, char *grpname, int grouptype, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgtis(fitsfile *fptr, char *grpname, int grouptype, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgtch(fitsfile *gfptr, int grouptype, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgtrm(fitsfile *gfptr, int rmopt, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgtcp(fitsfile *infptr, fitsfile *outfptr, int cpopt, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgtmg(fitsfile *infptr, fitsfile *outfptr, int mgopt, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgtcm(fitsfile *gfptr, int cmopt, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgtvf(fitsfile *gfptr, long *firstfailed, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgtop(fitsfile *mfptr,int group,fitsfile **gfptr,int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgtam(fitsfile *gfptr, fitsfile *mfptr, int hdupos, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgtnm(fitsfile *gfptr, long *nmembers, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgmng(fitsfile *mfptr, long *nmembers, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgmop(fitsfile *gfptr, long member, fitsfile **mfptr, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgmcp(fitsfile *gfptr, fitsfile *mfptr, long member, int cpopt, 
	   int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgmtf(fitsfile *infptr, fitsfile *outfptr,	long member, int tfopt,	       
	   int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int ffgmrm(fitsfile *fptr, long member, int rmopt, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }

/*--------------------- group template parser routines ------------------*/

int	fits_execute_template(fitsfile *ff, char *ngp_template, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }

/*--------------------- image compression routines ------------------*/

int fits_set_compression_type(fitsfile *fptr, int ctype, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int fits_set_tile_dim(fitsfile *fptr, int ndim, long *dims, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int fits_set_noise_bits(fitsfile *fptr, int noisebits, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }

int fits_get_compression_type(fitsfile *fptr, int *ctype, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int fits_get_tile_dim(fitsfile *fptr, int ndim, long *dims, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int fits_get_noise_bits(fitsfile *fptr, int *noisebits, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }

int fits_compress_img(fitsfile *infptr, fitsfile *outfptr, int compress_type,
         long *tilesize, int parm1, int parm2, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int fits_is_compressed_image(fitsfile *fptr, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }
int fits_decompress_img (fitsfile *infptr, fitsfile *outfptr, int *status)
  {
  g_error ("cfitsio not implemented");
  return 1;
  }


