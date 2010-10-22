/* $Id$        */
/*--------------------------------------------------------------------*/
/*;  Copyright (C) 2003-2010                                          */
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
#include "ObitPlot.h"
#include "ObitSkyGeom.h"
#include "ObitPosLabelUtil.h"

/*----------------Obit: Merx mollis mortibus nuper ------------------*/
/**
 * \file ObitPlot.c
 * ObitPlot Obit Graphics class definition.
 * 
 * Implementations in either plplot or pgplot are supported.
 * If HAVE_PLPLOT is defined plplot is used in preference to pgplot
 * pgplot is in Fortran and requires the use of a Fortran linker to 
 * build executables.
 * This class probably is not threadsafe.
 */

/*--------------- File Global Variables  ----------------*/
/** name of the class defined in this file */
static gchar *myClassName = "ObitPlot";

/** Function to obtain parent ClassInfo - Obit */
static ObitGetClassFP ObitParentGetClass = ObitGetClass;

/**
 * ClassInfo global structure ObitIOClassInfo.
 * This structure is used by class objects to access class functions.
 */
static ObitPlotClassInfo myClassInfo = {FALSE};

/*---------------Private function prototypes----------------*/
/** Private: Initialize newly instantiated object. */
void  ObitPlotInit  (gpointer in);

/** Private: Deallocate members. */
void  ObitPlotClear (gpointer in);

/** Private: Set Class function pointers. */
static void ObitPlotClassInfoDefFn (gpointer inClass);

#ifdef HAVE_PLPLOT  /* Only if plplot available */
/** Private: Parse output specification. */
static void parseOutput (gchar *output, gchar **dev, gchar **file);

/** Private: plplot reorderpixel values. */
static PLFLT** reorderPixels (ObitFArray *fa);

/** Private: plplot coordinate transform pixel to world coordinates. */
static void plplotCoord (PLFLT px, PLFLT py, PLFLT* wx, PLFLT *wy, 
			  PLPointer ID);
#endif /* HAVE_PLPLOT */
/*---------------Public functions---------------------------*/
/**
 * Construct Object.
 * \return pointer to object created.
 */
ObitPlot* newObitPlot (gchar *name)
{
  ObitPlot *out;
  gboolean OK=FALSE;  /* Have an implementation? */

   /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitPlotClassInit();

  /* Check that a plotting implementation is available */
#ifdef HAVE_PLPLOT  /* Only if plplot available */
  OK = TRUE;  /* Have a plotting package */
  /* Check plplot data size */
  if (sizeof(PLFLT)!=sizeof(ofloat)) {
    g_error("PLPLOT data incompatible with Obit");
    return NULL;
  }
#endif /* HAVE_PLPLOT */

#ifdef HAVE_PGPLOT  /* Only if pgplot available */
  OK = TRUE;  /* Have a plotting package */
#endif /* not HAVE_PGPLOT */

  /* Complain if plotting not available */
  if (!OK) {
    g_error("No plotting package available");
    return NULL;
  }

  /* allocate structure */
  out = g_malloc0(sizeof(ObitPlot));

  /* initialize values */
  if (name!=NULL) out->name = g_strdup(name);
  else out->name = g_strdup("Noname");

 /* set classInfo */
  out->ClassInfo = (gpointer)&myClassInfo;

  /* initialize other stuff */
  ObitPlotInit((gpointer)out);

 return out;
} /* end newObitPlot */

/**
 * Returns ClassInfo pointer for the class.
 * Initializes class if needed on first call.
 * \return pointer to the class structure.
 */
gconstpointer ObitPlotGetClass (void)
{
  /* Class initialization if needed */
  if (!myClassInfo.initialized) ObitPlotClassInit();

  return (gconstpointer)&myClassInfo;
} /* end ObitPlotGetClass */

/**
 * Initialize plot specifying output device
 * \param in      Pointer to object to be copied.
 * \param output  name and type of output device in form "filename/device"
 *                NULL => ask
 *                This doesn't work for pgPlot
 *                Devices for plplot
 * \li      xwin       X-Window (Xlib)
 * \li      gcw        Gnome Canvas Widget
 * \li      ps         PostScript File (monochrome)
 * \li      psc        PostScript File (color)
 * \li      xfig       Fig file
 * \li      png        PNG file
 * \li      jpeg       JPEG file
 * \li      gif        GIF file
 * \li      null       Null device
 * \param color  background Color index, not available for pgplot.
 *                0=black, 1=red(default), 2=yellow, 3=green, 
 *                4=aquamarine, 5=pink, 6=wheat, 7=gray, 8=brown,
 *                9=blue, 10=BlueViolet, 11=cyan, 12=turquoise
 *                13=magenta, 14=salmon, 15=white
 * \param nx      Number of frames in x on page
 * \param ny      Number of frames in y on page
 * \param err ObitErr error stack
 */
void ObitPlotInitPlot (ObitPlot* in, gchar *output, olong color, 
		       olong nx, olong ny, ObitErr *err)
{
  gboolean OK=FALSE;  /* Have an implementation? */

/****************** plplot implementation *************************/
#ifdef HAVE_PLPLOT  /* Only if plplot available */
  gchar *dev=NULL, *file=NULL;
  PLINT r, g, b, lcolor;
#endif /* end HAVE_PLPLOT */

/****************** pgplot implementation *************************/
#ifdef HAVE_PGPLOT  /* Only if pgplot available */
  olong istat;
  gchar *outstr;
  gchar *question="?";
  gchar *routine = "ObitPlotInitPlot";
#endif /* HAVE_PGPLOT */

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

  /* Save */
  in->nx = MAX (1,nx);
  in->ny = MAX (1,ny);

/****************** plplot implementation *************************/
#ifdef HAVE_PLPLOT  /* Only if plplot available */
  OK = TRUE;  /* Have a plotting package */

  /* Subpages?*/
  plssub(MAX(1,nx), MAX(1,ny));

  /* Output specified? */
  if (output) {
    /* Extract device and file name */
    parseOutput(output, &dev, &file);
    /* Set output device */
    if (dev) plsdev (dev);
    /* Set output file if any */
    if (file) plsetopt ("-o", file);
    /* Question of subpages */
    if (dev) g_free(dev);
    if (file) g_free(file);
  }
  /* Set background color */
  lcolor = (PLINT)color; 
  plgcol0(MAX(0, MIN(15,lcolor)), &r, &g, &b);
  plscolbg(r,g,b);

  plinit();   /* Init PLplot - it will ask if output not specified */

  /* Reset black as color 0
  plscolbg(0,0,0); apparently not */
#endif /* HAVE_PLPLOT */

/****************** pgplot implementation *************************/
#ifdef HAVE_PGPLOT  /* Only if pgplot available */
  OK = TRUE;  /* Have a plotting package */
  /* Create */
  if (output) outstr = output;
  else outstr = question;
  istat = cpgopen ((const char*)outstr);
  if (istat<=0) {
    Obit_log_error(err, OBIT_Error, "%s ERROR opening plot", 
		   routine);
    return;
  }

  /* Prompt for output */
  cpgask(output==NULL);

  /* Subdivide page */
  cpgsubp((int)nx, (int)ny);
  cpgpage();

#endif /* HAVE_PGPLOT */
 
  /* Complain if plotting not available */
  if (!OK) Obit_log_error(err, OBIT_Error, "No plotting package available");
 
} /* end ObitPlotInitPlot */

/**
 * Copy constructor.
 * \param in Pointer to object to be copied.
 * \param out Pointer to object to be written.  
 *            If NULL then a new structure is created.
 * \param err ObitErr error stack
 * \return Pointer to new object.
 */
ObitPlot* ObitPlotCopy (ObitPlot* in, ObitPlot* out, 
			    ObitErr *err)
{
  const ObitClassInfo *ParentClass;
  gboolean oldExist;
  gchar *outName;

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
    out = newObitPlot(outName);
    g_free(outName);
  }

  /* deep copy any base class members */
  ParentClass = myClassInfo.ParentClass;
  g_assert ((ParentClass!=NULL) && (ParentClass->ObitCopy!=NULL));
  ParentClass->ObitCopy (in, out, err);

  /* initialize/copy */
  return out;
} /* end ObitPlotCopy */

/**
 * Plot X vs Y using symbol.
 * Plot should be finalized with ObitPlotFinishPlot after all 
 * drawing on the current frame is finished.
 * This routine draws the frame and adds labels, to only overplot data
 * on the same frame, use ObitPlotXYOver
 * \param in      Pointer to Plot object.
 * \param symbol  Symbol index to use for plotting 
 *                values in the range [1,12] are usable 
 *                if negative, use abs value and connect points
 * \li 0 = line only
 * \li 1 = dot
 * \li 2 = plus
 * \li 3 = *
 * \li 4 = open circle
 * \li 5 = x
 * \li 6 = open square
 * \li 7 = open triangle
 * \li 8 = open star
 * \li 9 = filled triangle
 * \li 10 = filled square
 * \li 11 = filled circle
 * \li 12 = filled star
 *
 * \param n       Number of data points in x, y
 * \param x       Independent variable, if NULL use index
 * \param y       Dependent variable
 * \param err ObitErr error stack
 *
 * Optional parameters on in->info
 * \li XMAX (float) maximum X value (defaults to actual value)
 * \li XMIN (float) minimum X value (defaults to actual value)
 * \li YMAX (float) maximum Y value (defaults to actual value)
 * \li YMIN (float) minimum Y value (defaults to actual value)
 * \li TITLE (string)  Label for the plot (defaults to none), max 120
 * \li XLABEL (string) Label for horizontal axis (defaults to none)
 * \li XOPT   (string) Option for horizontal axis (default "BCNTS")
 *                     See #ObitPlotDrawAxes
 * \li YLABEL (string) Label for vertical axis (defaults to none)
 * \li YOPT   (string) Option for horizontal axis (default "BCNTS")
 *                     See #ObitPlotDrawAxes
 * \li XTICK   (float) world coordinate interval between major tick marks
 *                     on X axis. If xtick=0.0 [def], the interval is chosen.
 * \li NXSUB   (int)   the number of subintervals to divide the major
 *                     coordinate interval into. If xtick=0.0 or nxsub=0,
 *                     the number is chosen. [def 0]
 * \li YTICK  (float)  like xtick for the Y axis.
 * \li NYSUB  (int)    like nxsub for the Y axis
 * \li CSIZE  (int)    Scaling factor for characters(default = 1)
 * \li SSIZE  (int)    Scaling factor for symbols(default = 1)
 * \li LWIDTH (int)    Line width (default = 1)
 * \li JUST   (int)    If !=0 then force X and Y axis scaling to be the same
 */
void ObitPlotXYPlot (ObitPlot* in, olong symbol, 
		     olong n, ofloat *x, ofloat *y, ObitErr *err)
{
  gboolean OK=FALSE;  /* Have an implementation? */
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  gboolean gotRange, gotXMAX, gotYMAX, gotXMIN, gotYMIN, doConnect;
  gboolean doSymbol;
  ofloat xtick, ytick;
  olong nxsub, nysub;
  gchar title[121], xlabel[121], ylabel[121], xopt[31], yopt[31];
  gchar *routine="ObitPlotXYPlot";

/****************** plplot implementation *************************/
#ifdef HAVE_PLPLOT  /* Only if plplot available */
  olong i, just, axis;
  ofloat xmax, xmin, ymax, ymin, cscale, fblank = ObitMagicF();
  olong csize, ssize, lwidth, lsymbol;
  olong syms[12] = {210,225,228,840,227,841,842,844,852,851,850,856};
  gchar *optemp, *xopt_def="BCNTS", *yopt_def="BCNTS";
  PLFLT sscale, *xp, *yp, *xx, *ply=NULL, *plx=NULL, *logx=NULL, *logy=NULL;

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (y!=NULL);

  /* anything to do? */
  if (n<=0) return;
  
  OK = TRUE;  /* Have a plotting package */
  lsymbol = MAX(1, MIN(12, abs(symbol)));
  if (lsymbol>0) lsymbol = syms[lsymbol-1];  /* Find plplot value */
  doConnect = symbol<=0;
  doSymbol  = symbol!=0;
  in->xLog = FALSE;
  in->yLog = FALSE;

  /* using index for x? */
  if (x!=0) xx = x;
  else { /* yes */
    xx = g_malloc0(n*sizeof(PLFLT));
    for (i=0; i<n; i++) xx[i] = (PLFLT)i;
  }

  /* Get inputs from info */
  /* Plot title */
  for (i=0; i<121; i++) title[i] = 0;
  strcpy (title, "    ");
  ObitInfoListGetTest(in->info, "TITLE", &type, dim, title);
 
  /* Horizontal label */
  for (i=0; i<121; i++) xlabel[i] = 0;
  strcpy (xlabel, "    ");
  ObitInfoListGetTest(in->info, "XLABEL", &type, dim, xlabel);
 
  /* Horizontal label options */
  for (i=0; i<31; i++) xopt[i] = 0;
  strcpy (xopt, xopt_def);
  ObitInfoListGetTest(in->info, "XOPT", &type, dim, xopt);
  for (i=0; i<30; i++) if(xopt[i]=='L') in->xLog = TRUE;  /* Log plot? */
  optemp = g_ascii_strdown (xopt,-1);  /* Convert to lower case */
  strcpy (xopt, optemp);
  g_free(optemp);
  xtick = 0.0;
  ObitInfoListGetTest(in->info, "XTICK", &type, dim, &xtick);
  nxsub = 0;
  ObitInfoListGetTest(in->info, "NXSUB", &type, dim, &nxsub);

  /* Vertical label */
  for (i=0; i<121; i++) ylabel[i] = 0;
  strcpy (ylabel, "    ");
  ObitInfoListGetTest(in->info, "YLABEL", &type, dim, ylabel);

  /* Vertical label options */
  for (i=0; i<31; i++) yopt[i] = 0;
  strcpy (yopt, yopt_def);
  ObitInfoListGetTest(in->info, "YOPT", &type, dim, yopt);
  for (i=0; i<30; i++) if(yopt[i]=='L') in->yLog = TRUE;  /* Log plot? */
  optemp = g_ascii_strdown (yopt,-1);  /* Convert to lower case */
  strcpy (yopt, optemp);
  g_free(optemp);
  ytick = 0.0;
  ObitInfoListGetTest(in->info, "YTICK", &type, dim, &ytick);
  nysub = 0;
  ObitInfoListGetTest(in->info, "NYSUB", &type, dim, &nysub);

  /* Sizes */
  csize = 1;
  ObitInfoListGetTest(in->info, "CSIZE", &type, dim, (gpointer*)&csize);
  ssize = 1;
  ObitInfoListGetTest(in->info, "SSIZE", &type, dim, (gpointer*)&ssize);
  lwidth = 1;
  ObitInfoListGetTest(in->info, "LWIDTH", &type, dim, (gpointer*)&lwidth);
  just = 0;
  ObitInfoListGetTest(in->info, "JUST", &type, dim, (gpointer*)&just);
  if (just!=0) just = 1;
 
  /* Plot range */
  xmax = fblank;
  ObitInfoListGetTest(in->info, "XMAX", &type, dim, (gpointer*)&xmax);
  ymax = fblank;
  ObitInfoListGetTest(in->info, "YMAX", &type, dim, (gpointer*)&ymax);
  xmin = fblank;
  ObitInfoListGetTest(in->info, "XMIN", &type, dim, (gpointer*)&xmin);
  ymin = fblank;
  ObitInfoListGetTest(in->info, "YMIN", &type, dim, (gpointer*)&ymin);

  /* Check data for range? */
  gotXMAX = xmax != fblank;
  gotYMAX = ymax != fblank;
  gotXMIN = xmin != fblank;
  gotYMIN = ymin != fblank;
  gotRange =  gotXMAX && gotYMAX && gotXMIN && gotYMIN;
  if (!gotRange) { /* Need to check data? */
    for (i=0; i<n; i++) {
      if ((x!=NULL) && (x[i]!=fblank)) {
	if  (xmax == fblank) xmax = x[i];
	else xmax = MAX (x[i], xmax);
	if  (xmin == fblank) xmin = x[i];
	else xmin = MIN (x[i], xmin);
      } else if (x==NULL) {
	if  (xmax == fblank) xmax = xx[i];
	else xmax = MAX (xx[i], xmax);
	if  (xmin == fblank) xmin = xx[i];
	else xmin = MIN (xx[i], xmin);
      }
      if (y[i]!=fblank) {
	if  (ymax == fblank) ymax = y[i];
	else ymax = MAX (y[i], ymax);
	if  (ymin == fblank) ymin = y[i];
	else ymin = MIN (y[i], ymin);
      }
    }
  } /* end check range */

  /* Are the range values OK? */
  if (xmax<=xmin) { /* Nope - bug out */
    Obit_log_error(err, OBIT_Error, "%s: XMAX(%g) <= XMIN(%g)", 
		   routine, xmax, xmin);
    if (x==NULL) g_free(xx);
    return;
 }

  if (ymax<=ymin) { /* Nope - bug out */
    Obit_log_error(err, OBIT_Error, "%s: YMAX (%g) <= YMIN(%g)", 
		   routine, ymax, ymin);
    if (x==NULL) g_free(xx);
    return;
 }

  /* Set axis options */
  /* Log? */
  if (in->xLog) {
    /* Must be positive */
    if ((xmax<=0.0) || (xmin<=0.0)) {
      Obit_log_error(err, OBIT_Error, 
		     "%s: Values for log of x axis must be positive %f %f", 
		      routine, xmin, xmax);
      if (x==NULL) g_free(xx);
      return;
    }
    xmin = log10(xmin);
    xmax = log10(xmax);
    logx = g_malloc0(n*sizeof(PLFLT));
    for (i=0; i<n; i++) logx[i] = (PLFLT)log10(xx[i]);
    xp = logx;
  } else {
    plx = g_malloc0(n*sizeof(PLFLT));
    for (i=0; i<n; i++) plx[i] = (PLFLT)xx[i];
    xp = plx;
  }

  if (in->yLog) {
    /* Must be positive */
    if ((ymax<=0.0) || (ymin<=0.0)) {
      Obit_log_error(err, OBIT_Error, 
		     "%s: Values for log of y axis must be positive %f %f", 
		     routine, ymin, ymax);
      if (x==NULL) g_free(xx);
      return;
    }
    ymin = log10(ymin);
    ymax = log10(ymax);
    logy = g_malloc0(n*sizeof(PLFLT));
    for (i=0; i<n; i++) {
      if (y[i]!=fblank) logy[i] = (PLFLT)log10(y[i]);
      else logy[i] = (PLFLT)xmax*10.0;  /* Don't plot if blanked */
    }
    yp = logy;
  } else {
    ply = g_malloc0(n*sizeof(PLFLT));
    for (i=0; i<n; i++) ply[i] = (PLFLT)y[i];
    yp = ply;
  }
  
  /* If autosetting the range, expand a bit */
  if (!gotXMAX) xmax += 0.05 * (xmax-xmin);
  if (!gotXMIN) xmin -= 0.05 * (xmax-xmin);
  if (!gotYMAX) ymax += 0.05 * (ymax-ymin);
  if (!gotYMIN) ymin -= 0.05 * (ymax-ymin);
  
  /* Adjust character size */
  cscale = (ofloat)csize;
  cscale = MAX (1.0, cscale);
  ObitPlotSetCharSize (in, cscale, err);
  
  /* set line width */
  lwidth = MAX (1, lwidth);
  ObitPlotSetLineWidth (in, lwidth, err);
  
  /* set plotting area */
  axis = -2;
  if (in->xLog  && !in->yLog) axis = 10;
  if (!in->xLog &&  in->yLog) axis = 20;
  if ( in->xLog &&  in->yLog) axis = 30;
  ObitPlotSetPlot (in, xmin, xmax, ymin, ymax, just, axis, err);
  
  /*  Plot label */
  ObitPlotLabel (in, xlabel, ylabel, title, err);
  
  /* set ticks */
  ObitPlotDrawAxes (in, xopt, xtick, nxsub, yopt, ytick, nysub, err);
  
  /* Connect the dots? */
  if (doConnect) plline((PLINT)n, xp, yp);
  
  /* Adjust symbol size */
  sscale = (PLFLT)ssize;
  sscale = MAX (1.0, sscale);
  c_plssym((PLFLT)0.0, sscale);
  
  /* Plot symbols */
  if (doSymbol) plsym ((PLINT)n, xp, yp, (PLINT)lsymbol);
  
  /* Flush the buffer */
  plflush();

  /* deallocate arrays if allocated */
  if (x==NULL)     g_free(xx);
  if (logx!=NULL)  g_free(logx);
  if (logy!=NULL)  g_free(logy);
  if (plx!=NULL)   g_free(plx);
  if (ply!=NULL)   g_free(ply);
#endif /* HAVE_PLPLOT */

/****************** pgplot implementation *************************/
#ifdef HAVE_PGPLOT  /* Only if pgplot available */
  olong i, just, axis;
  float *ply=NULL, *plx=NULL;
  ofloat sscale, xmax, xmin, ymax, ymin, cscale, *xx, fblank = ObitMagicF();
  olong csize, ssize, lwidth, lsymbol;
  olong syms[12] = {1,2,3,4,5,6,7,12,13,16,17,18};
  gchar *xopt_def="BCNTS", *yopt_def="BCNTS";

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (y!=NULL);


  /* anything to do? */
  if (n<=0) return;
  
  OK = TRUE;  /* Have a plotting package */
  lsymbol = MAX(1, MIN(12, abs(symbol)));
  lsymbol = syms[lsymbol-1];  /* Find plplot value */
  doConnect = symbol<=0;
  doSymbol  = symbol!=0;
  in->xLog = FALSE;
  in->yLog = FALSE;

  /* using index for x? */
  if (x!=0) xx = x;
  else { /* yes */
    xx = g_malloc0(n*sizeof(ofloat));
    for (i=0; i<n; i++) xx[i] = (float)i;
  }

  /* Get inputs from info */
  /* Plot title */
  for (i=0; i<121; i++) title[i] = 0;
  strcpy (title, "    ");
  ObitInfoListGetTest(in->info, "TITLE", &type, dim, title);
 
  /* Horizontal label */
  for (i=0; i<121; i++) xlabel[i] = 0;
  strcpy (xlabel, "    ");
  ObitInfoListGetTest(in->info, "XLABEL", &type, dim, xlabel);
 
  /* Horizontal label options */
  for (i=0; i<31; i++) xopt[i] = 0;
  strcpy (xopt, xopt_def);
  ObitInfoListGetTest(in->info, "XOPT", &type, dim, xopt);
  for (i=0; i<30; i++) if(xopt[i]=='L') in->xLog = TRUE;  /* Log plot? */
  xtick = 0.0;
  ObitInfoListGetTest(in->info, "XTICK", &type, dim, &xtick);
  nxsub = 0;
  ObitInfoListGetTest(in->info, "NXSUB", &type, dim, &nxsub);

  /* Vertical label */
  for (i=0; i<121; i++) ylabel[i] = 0;
  strcpy (ylabel, "    ");
  ObitInfoListGetTest(in->info, "YLABEL", &type, dim, ylabel);

  /* Vertical label options */
  for (i=0; i<31; i++) yopt[i] = 0;
  strcpy (yopt, yopt_def);
  ObitInfoListGetTest(in->info, "YOPT", &type, dim, yopt);
  for (i=0; i<30; i++) if(yopt[i]=='L') in->yLog = TRUE;  /* Log plot? */
  ytick = 0.0;
  ObitInfoListGetTest(in->info, "YTICK", &type, dim, &ytick);
  nysub = 0;
  ObitInfoListGetTest(in->info, "NYSUB", &type, dim, &nysub);

  /* Sizes */
  csize = 1;
  ObitInfoListGetTest(in->info, "CSIZE", &type, dim, (gpointer*)&csize);
  ssize = 1;
  ObitInfoListGetTest(in->info, "SSIZE", &type, dim, (gpointer*)&ssize);
  lwidth = 1;
  ObitInfoListGetTest(in->info, "LWIDTH", &type, dim, (gpointer*)&lwidth);
  just = 0;
  ObitInfoListGetTest(in->info, "JUST", &type, dim, (gpointer*)&just);
  just = (just!=0);
 
  /* Plot range */
  xmax = fblank;
  ObitInfoListGetTest(in->info, "XMAX", &type, dim, (gpointer*)&xmax);
  ymax = fblank;
  ObitInfoListGetTest(in->info, "YMAX", &type, dim, (gpointer*)&ymax);
  xmin = fblank;
  ObitInfoListGetTest(in->info, "XMIN", &type, dim, (gpointer*)&xmin);
  ymin = fblank;
  ObitInfoListGetTest(in->info, "YMIN", &type, dim, (gpointer*)&ymin);

  /* Check data for range? */
  gotXMAX = xmax != fblank;
  gotYMAX = ymax != fblank;
  gotXMIN = xmin != fblank;
  gotYMIN = ymin != fblank;
  gotRange =  gotXMAX && gotYMAX && gotXMIN && gotYMIN;
  if (!gotRange) { /* Need to check data? */
    for (i=0; i<n; i++) {
      if (xx[i]!=fblank) {
	if  (xmax == fblank) xmax = xx[i];
	else xmax = MAX (xx[i], xmax);
	if  (xmin == fblank) xmin = xx[i];
	else xmin = MIN (xx[i], xmin);
      }
      if (y[i]!=fblank) {
	if  (ymax == fblank) ymax = y[i];
	else ymax = MAX (y[i], ymax);
	if  (ymin == fblank) ymin = y[i];
	else ymin = MIN (y[i], ymin);
      }
    }
  } /* end check range */

  /* Are the range values OK? */
  if (xmax<=xmin) { /* Nope - bug out */
    Obit_log_error(err, OBIT_Error, "%s: XMAX(%g) <= XMIN(%g)", 
		   routine, xmax, xmin);
    if (x==NULL) g_free(xx);
    return;
 }

  if (ymax<=ymin) { /* Nope - bug out */
    Obit_log_error(err, OBIT_Error, "%s: YMAX (%g) <= YMIN(%g)", 
		   routine, ymax, ymin);
    if (x==NULL) g_free(xx);
    return;
 }

  /* If autosetting the range, expand a bit */
  if (!gotXMAX) xmax += 0.05 * (xmax-xmin);
  if (!gotXMIN) xmin -= 0.05 * (xmax-xmin);
  if (!gotYMAX) ymax += 0.05 * (ymax-ymin);
  if (!gotYMIN) ymin -= 0.05 * (ymax-ymin);

  /* Adjust character size */
  cscale = (ofloat)csize;
  cscale = MAX (1.0, cscale);
  ObitPlotSetCharSize (in, cscale, err);

  /* set line width */
  lwidth = MAX (1, lwidth);
  ObitPlotSetLineWidth (in, lwidth, err);

  /* set plotting area */
  axis = -2;
  ObitPlotSetPlot (in, xmin, xmax, ymin, ymax, just, axis, err);

  /* Init - use buffering */
  cpgbbuf();

  /*  Plot label */
  ObitPlotLabel (in, xlabel, ylabel, title, err);

  /* set ticks */
  ObitPlotDrawAxes (in, xopt, xtick, nxsub, yopt, ytick, nysub, err);

  /* Adjust symbol size */
  sscale = (ofloat)ssize;
  sscale = MAX (1.0, sscale);
  ObitPlotSetCharSize (in, sscale, err);

   /* Want symbols? */
 if (doSymbol) {
   /* To pgplot data types */
   ply = g_malloc0(n*sizeof(float));
   for (i=0; i<n; i++) ply[i] = (float)y[i];
   plx = g_malloc0(n*sizeof(float));
   for (i=0; i<n; i++) plx[i] = (float)xx[i];
   
   /* Plot */
   cpgpt ((int)n, plx, ply, (int)lsymbol);
 }
  /* Connect the dots? */
  if (doConnect) {
    for (i=1; i<n; i++) 
      ObitPlotDrawLine (in, xx[i-1], y[i-1], xx[i], y[i], err);
  }

  /* Reset symbol size */
  ObitPlotSetCharSize (in, cscale, err);

  /* Flush the buffer */
  cpgebuf();
  cpgupdt();

  /* deallocate x index array if allocated */
  if (x==NULL)   g_free(xx);
  if (plx!=NULL) g_free(plx);
  if (ply!=NULL) g_free(ply);
#endif /* HAVE_PGPLOT */

  /* Complain if plotting not available */
  if (!OK) Obit_log_error(err, OBIT_Error, "No plotting package available");

} /* end  ObitPlotXYPlot */

/**
 * Overplot on plot previously generated by ObitPlotXYPlot
 * Plot should be finalized with ObitPlotFinishPlot after all 
 * drawing on the current frame is finished.
 * \param in      Pointer to Plot object.
 * \param symbol  Symbol index to use for plotting 
 *                values in the range [1,12] are usable 
 *                if negative, use abs value and connect points
 * \li 0 = line only
 * \li 1 = dot
 * \li 2 = plus
 * \li 3 = *
 * \li 4 = open circle
 * \li 5 = x
 * \li 6 = open square
 * \li 7 = open triangle
 * \li 8 = open star
 * \li 9 = filled triangle
 * \li 10 = filled square
 * \li 11 = filled circle
 * \li 12 = filled star
 *
 * \param n       Number of data points in x, y
 * \param x       Independent variable, if NULL use index
 * \param y       Dependent variable
 * \param err ObitErr error stack
 *
 * Optional parameters on in->info
 * \li CSIZE  (int)    Scaling factor for characters(default = 1)
 * \li LWIDTH (int)    Line width (default = 1)
 */
void ObitPlotXYOver (ObitPlot* in, olong symbol, 
		     olong n, ofloat *x, ofloat *y, ObitErr *err)
{
  gboolean OK=FALSE;  /* Have an implementation? */
  gboolean doSymbol;
  /*gchar *routine="ObitPlotXYOver";*/

/****************** plplot implementation *************************/
#ifdef HAVE_PLPLOT  /* Only if plplot available */
  olong i;
  ofloat cscale;
  olong csize, lwidth;
  ofloat fblank = ObitMagicF();
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  gboolean doConnect;
  olong lsymbol, syms[12] = {210,225,228,840,227,841,842,844,852,851,850,856};
  PLFLT *xp, *yp, *xx, *ply=NULL, *plx=NULL, *logx=NULL, *logy=NULL;
  gchar *routine="ObitPlotXYOver";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (y!=NULL);


  /* anything to do? */
  if (n<=0) return;
  
  OK = TRUE;  /* Have a plotting package */
  lsymbol = MAX(1, MIN(12, abs(symbol)));
  lsymbol = syms[lsymbol-1];  /* Find plplot value */
  doConnect = symbol<=0;
  doSymbol  = symbol!=0;

  /* using index for x? */
  if (x!=0) xx = x;
  else { /* yes */
    xx = g_malloc0(n*sizeof(PLFLT));
    for (i=0; i<n; i++) xx[i] = (PLFLT)i;
  }

  /* Sizes */
  csize = 1;
  ObitInfoListGetTest(in->info, "CSIZE", &type, dim, (gpointer*)&csize);
  lwidth = 1;
  ObitInfoListGetTest(in->info, "LWIDTH", &type, dim, (gpointer*)&lwidth);
 
  /* Adjust character size */
  cscale = (ofloat)csize;
  cscale = MAX (1.0, cscale);
  ObitPlotSetCharSize (in, cscale, err);

  /* set line width */
  lwidth = MAX (1, lwidth);
  ObitPlotSetLineWidth (in, lwidth, err);

 /* Log? */
  if (in->xLog) {
    logx = g_malloc0(n*sizeof(PLFLT));
    for (i=0; i<n; i++) {
      /* Must be positive */
      if (xx[i]<=0.0) {
	Obit_log_error(err, OBIT_Error, 
		       "%s: Values for log of x axis must be positive %f", 
		       routine,xx[i]);
	if (logx!=NULL) g_free(logx);
	return;
      }
     logx[i] = (PLFLT)log10(xx[i]);
    }
    xp = logx;
  } else {
    plx = g_malloc0(n*sizeof(PLFLT));
    for (i=0; i<n; i++) plx[i] = (PLFLT)xx[i];
    xp = plx;
  }

  if (in->yLog) {
    logy = g_malloc0(n*sizeof(PLFLT));
    for (i=0; i<n; i++) {
      if (y[i]!=fblank) {
	/* Must be positive */
	if (y[i]<=0.0) {
	  Obit_log_error(err, OBIT_Error, 
			 "%s: Values for log of y axis must be positive %f", 
			 routine, y[i]);
	  if (logx!=NULL) g_free(logx);
	  if (logy!=NULL) g_free(logy);
	  return;
	}
	logy[i] = (PLFLT)log10(y[i]);
      }
      else logy[i] = (PLFLT)1.0e30;  /* Don't plot if blanked */
    }
    yp = logy;
  } else {
    ply = g_malloc0(n*sizeof(PLFLT));
    for (i=0; i<n; i++) ply[i] = (PLFLT)y[i];
    yp = ply;
  }

  /* Connect the dots? */
  if (doConnect) plline(n, xp, yp);

  /* Plot symbols */
  if (doSymbol) plsym (n, xp, yp, lsymbol);

  /* Flush the buffer */
  plflush();

  /* deallocate arrays if allocated */
  if (x==NULL)     g_free(xx);
  if (logx!=NULL)  g_free(logx);
  if (logy!=NULL)  g_free(logy);
  if (plx!=NULL)   g_free(plx);
  if (ply!=NULL)   g_free(ply);
#endif /* HAVE_PLPLOT */

/****************** pgplot implementation *************************/
#ifdef HAVE_PGPLOT  /* Only if pgplot available */
  olong i;
  float *ply=NULL, *plx=NULL;
  ofloat cscale, *xx;
  olong csize, lwidth;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  olong lsymbol, syms[12] = {1,2,3,4,5,6,7,12,13,16,17,18};
  gboolean doConnect;
  /*gchar *routine="ObitPlotXYOver";*/

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (y!=NULL);


  /* anything to do? */
  if (n<=0) return;
  
  OK = TRUE;  /* Have a plotting package */
  lsymbol = MAX(1, MIN(12, abs(symbol)));
  lsymbol = syms[lsymbol-1];  /* Find pgplot value */
  doConnect = symbol<=0;
  doSymbol  = symbol!=0;

  /* using index for x? */
  if (x!=0) xx = x;
  else { /* yes */
    xx = g_malloc0(n*sizeof(ofloat));
    for (i=0; i<n; i++) xx[i] = (float)i;
  }

  /* Sizes */
  csize = 1;
  ObitInfoListGetTest(in->info, "CSIZE", &type, dim, (gpointer*)&csize);
  lwidth = 1;
  ObitInfoListGetTest(in->info, "LWIDTH", &type, dim, (gpointer*)&lwidth);
 
  /* Adjust character size */
  cscale = (ofloat)csize;
  cscale = MAX (1.0, cscale);
  ObitPlotSetCharSize (in, cscale, err);

  /* set line width */
  lwidth = MAX (1, lwidth);
  ObitPlotSetLineWidth (in, lwidth, err);

  /* Want symbols? */
 if (doSymbol) {
   /* To pgplot data types */
   ply = g_malloc0(n*sizeof(float));
   for (i=0; i<n; i++) ply[i] = (float)y[i];
   plx = g_malloc0(n*sizeof(float));
   for (i=0; i<n; i++) plx[i] = (float)xx[i];
   
   /* Plot */
   cpgpt ((int)n, plx, ply, (int)lsymbol);
 }

  /* Connect the dots? */
  if (doConnect) {
    for (i=1; i<n; i++) 
      ObitPlotDrawLine (in, xx[i-1], y[i-1], xx[i], y[i], err);
  }

  /* Flush the buffer */
  cpgebuf();
  cpgupdt();

  /* deallocate x index array if allocated */
  if (x==NULL)   g_free(xx);
  if (plx!=NULL) g_free(plx);
  if (ply!=NULL) g_free(ply);

#endif /* HAVE_PGPLOT */

  /* Complain if plotting not available */
 if (!OK) Obit_log_error(err, OBIT_Error, "No plotting package available");
} /* end  ObitPlotXYOver */

/**
 * Contour plot of an ObitImage
 * Negative contours will be dashed
 * Plot should be finalized with ObitPlotFinishPlot after all 
 * drawing on the current frame is finished.
 * \param in      Pointer to existing ObitPlot object.
 * \param label   Label for plot
 * \param image   Image to plot (first plane in BLC,TRC)
 *                Rotated images aren't done quite right
 * \param lev     basic contour level (def 0.1 peak)
 * \param cntfac  Contour level factor (def sqrt(2)
 * \param err ObitErr error stack
 *
 * Optional parameters on in->info
 * \li XTICK   (float) world coordinate interval between major tick marks
 *                     on X axis. If xtick=0.0 [def], the interval is chosen.
 * \li NXSUB   (int)   the number of subintervals to divide the major
 *                     coordinate interval into. If xtick=0.0 or nxsub=0,
 *                     the number is chosen. [def 0]
 * \li YTICK  (float)  like xtick for the Y axis.
 * \li NYSUB  (int)    like nxsub for the Y axis
 * \li CSIZE  (int)    Scaling factor for characters(default = 1)
 * \li LWIDTH (int)    Line width (default = 1)
 */
void ObitPlotContour (ObitPlot* in, gchar *label, ObitImage *image,
		      ofloat lev, ofloat cntfac, ObitErr *err)
{
  gboolean OK=FALSE;  /* Have an implementation? */
  gchar xlabel[121], ylabel[121];
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ofloat xtick, ytick;
  olong nxsub, nysub;
  ObitInfoType type;
  gchar *routine="ObitPlotContour";

/****************** plplot implementation *************************/
#ifdef HAVE_PLPLOT  /* Only if plplot available */
  ObitImageDesc *id;
  olong i, just, axis, npc, nnc;
  ofloat xmax, xmin, ymax, ymin, cscale;
  PLFLT px, py, wx, wy, **map=NULL, *levs=NULL;
  ofloat maxval, minval;
  olong pos[2];
  olong csize, lwidth, nx, ny;
  gchar rast[16], decst[16];
  gchar units[20], *xopt_def="bcnts", *yopt_def="bcnts", pre;


 /* error checks */
  if (err->error) return;
  g_assert (ObitPlotIsA(in));
  g_assert (ObitImageIsA(image));
  
  OK = TRUE;  /* Have a plotting package */
  in->xLog = FALSE;
  in->yLog = FALSE;

  /* Save info */
  in->myImage = ObitImageRef(image);
  in->myErr   = ObitErrRef(err);

  /* Open, read image */
  ObitImageOpen (image, OBIT_IO_ReadOnly, err);
  ObitImageRead (image, NULL, err);
  if (err->error) Obit_traceback_msg (err, routine, image->name);
  
  /* Give warning if image rotated */
  if (image->myDesc->crota[1]!=0.0) {
    Obit_log_error(err, OBIT_InfoWarn, "%s: Rotated image not handled well", 
		   routine);
  }

  /* Get inputs from info */
  xtick = 0.0;
  ObitInfoListGetTest(in->info, "XTICK", &type, dim, &xtick);
  nxsub = 0;
  ObitInfoListGetTest(in->info, "NXSUB", &type, dim, &nxsub);
  ytick = 0.0;
  ObitInfoListGetTest(in->info, "YTICK", &type, dim, &ytick);
  nysub = 0;
  ObitInfoListGetTest(in->info, "NYSUB", &type, dim, &nysub);
  /* Sizes */
  csize = 1;
  ObitInfoListGetTest(in->info, "CSIZE", &type, dim, &csize);
  lwidth = 1;
  ObitInfoListGetTest(in->info, "LWIDTH", &type, dim, &lwidth);
  /* Adjust character size */
  cscale = (ofloat)csize;
  cscale = MAX (1.0, cscale);
  ObitPlotSetCharSize (in, cscale, err);

  /* set line width */
  lwidth = MAX (1, lwidth);
  ObitPlotSetLineWidth (in, lwidth, err);

  /* What units/scaling for axes? */
  id = image->myDesc;
  if (id->inaxes[0]*fabs(id->cdelt[0])>1.0) {
    /* degrees */
    strcpy (units,"deg");
    in->scalex = 1.0;
    in->scaley = 1.0;
  } else if (id->inaxes[0]*fabs(id->cdelt[0])>0.016667) {
    /* arcmin */
    strcpy (units,"amin");
    in->scalex = 60.0;
    in->scaley = 60.0;
  } else if (id->inaxes[0]*fabs(id->cdelt[0])>0.0027778) {
    /* arcsec */
    strcpy (units,"asec");
    in->scalex = 3600.0;
    in->scaley = 3600.0;
  } else {
    /* milliarcsec  */
    strcpy (units,"mas");
    in->scalex = 3600000.0;
    in->scaley = 3600000.0;
  }

  /* Correct X scaling for the artificial effect of declination */
  if (id->jlocd>=0) 
    in->scalex *= cos(DG2RAD*id->crval[id->jlocd]);


  nx = image->myDesc->inaxes[0];
  ny = image->myDesc->inaxes[1];
  
  /* Define size of plot in world coordinates */
  px = (PLFLT)1.0; py = (PLFLT)1.0;
  plplotCoord (px, py, &wx, &wy, (PLPointer)in);
  xmin = (ofloat)wx;
  ymin = (ofloat)wy;
  px = (PLFLT)nx; py = (PLFLT)ny;
  plplotCoord (px, py, &wx, &wy, (PLPointer)in);
  xmax = (ofloat)wx;
  ymax = (ofloat)wy;

  /* set plotting area */
  axis = -2;
  just = 1;
  ObitPlotSetPlot (in, xmin, xmax, ymin, ymax, just, axis, err);

  /*  Plot labels */
  if (image->myDesc->equinox<1975.0) pre = 'B';
  else pre = 'J';
  if (fabs(image->myDesc->crota[1])==0.0)
    g_snprintf (xlabel,120,"Right Ascension (%s)", units);
  else /* Don't call it RA */
    g_snprintf (xlabel,120,"Rotated Right Ascension (%s)", units);
  if (fabs(image->myDesc->crota[1])==0.0)
    g_snprintf (ylabel,120,"Declination (%s)", units);
  else  /* Don't call it Dec */
    g_snprintf (ylabel,120,"Rotated Declination (%s)", units);
  ObitPlotLabel (in, xlabel, ylabel, label, err);

  /* Get center position string */
  ObitPosLabelUtilRA2HMS  (image->myDesc->crval[0], image->myDesc->ctype[0], rast);
  ObitPosLabelUtilDec2DMS (image->myDesc->crval[1], image->myDesc->ctype[1], decst);

 /* set ticks, labels */
  ObitPlotDrawAxes (in, xopt_def, xtick, nxsub, yopt_def, ytick, nysub, err);

  /* Set contour levels */
  maxval = ObitFArrayMax (image->image, pos);
  minval = ObitFArrayMin (image->image, pos);
  if (lev <= 0.0) lev = 0.1 * MAX (fabs(maxval), fabs(minval));
  if (cntfac <= 0.0) cntfac = sqrt(2);
  if (maxval>0.0) {
    npc = 0.999 + log(maxval/lev) / log(cntfac);
    npc = MAX (0, npc);
  } else  npc = 0;
  if (minval<0.0) {
    nnc = 0.999  + log(-minval/lev) / log(cntfac);
    nnc = MAX (0, nnc);
  } else nnc = 0;
  levs = g_malloc0(MAX (npc, nnc)*sizeof(PLFLT));

  /* Convert Obit image to plplot array */
  map = reorderPixels(image->image);

  /* Plot positive contours */
  if (npc>0) {
    levs[0] = (PLFLT)lev;
    for (i=1; i<npc; i++) levs[i] = levs[i-1]*cntfac;
    pllsty ((PLINT)1);  /* Solid lines */
    plcont (map, (PLINT)nx, (PLINT)ny, (PLINT)1, (PLINT)nx, 
	    (PLINT)1, (PLINT)ny, levs, (PLINT)npc, 
	    plplotCoord, (PLPointer)in);
  }

  /* Plot negative contours */
  if (npc>0) {
    levs[0] = (PLFLT)(-lev);
    for (i=1; i<nnc; i++) levs[i] = levs[i-1]*cntfac;
    pllsty ((PLINT)2); /* dashed lines */
    plcont (map, (PLINT)nx, (PLINT)ny, (PLINT)1, (PLINT)nx, 
	    (PLINT)1, (PLINT)ny, levs, (PLINT)nnc, 
	    plplotCoord, (PLPointer)in);
    pllsty ((PLINT)1); /* back to solid lines */
 }

  /* Set contouring info line */
  if (fabs(image->myDesc->crota[1])>0.001)
    g_snprintf (xlabel,120,"Levs=%g*%6.2f**n, Peak=%g %s, rot=%6.1f", 
		lev, cntfac, maxval, image->myDesc->bunit, 
		image->myDesc->crota[1]);
  else
    g_snprintf (xlabel,120,"Levs=%g*%6.2f**n, Peak=%g %s", 
		lev, cntfac, maxval, image->myDesc->bunit);

  just = 0;
  ObitPlotSetCharSize (in, 0.5, err);
  ObitPlotRelText (in, "B", 8.5, 0.0, just, xlabel, err);
  /* Center position */
  g_snprintf (xlabel,120,"Center = %s, %s %c%6.1f", 
	      rast, decst, pre, image->myDesc->equinox);
  ObitPlotRelText (in, "T", 1.0, 0.0, just, xlabel, err);

  if (levs) g_free(levs); /* deallocate levels array */
  /* Deallocate plplot array */
  if (map) {
    for (i=0; i<nx; i++) g_free(map[i]);
    g_free(map);
  }

  /* Flush the buffer */
  plflush();

  /* Close image */
  ObitImageClose (image, err);
  if (err->error) Obit_traceback_msg (err, routine, image->name);
  /* Free image buffer */
  image->image = ObitFArrayUnref(image->image);
#endif /* HAVE_PLPLOT */

/****************** pgplot implementation *************************/
#ifdef HAVE_PGPLOT  /* Only if pgplot available */
  olong i, just, axis, npc, nnc;
  ofloat xmax, xmin, ymax, ymin, cscale, maxval, minval, fblank = ObitMagicF();
  float tr[6];
  float *levs=NULL, *fmap=NULL;
  olong pos[2];
  ofloat dx, dy, cr, sr;
  olong csize, lwidth, nx, ny;
  gchar *xopt_def="BCNTSZYH", *yopt_def="BCNTSZYD", pre;


 /* error checks */
  if (err->error) return;
  g_assert (ObitPlotIsA(in));
  g_assert (ObitImageIsA(image));
  
  OK = TRUE;  /* Have a plotting package */
  in->xLog = FALSE;
  in->yLog = FALSE;

  /* Save info */
  in->myImage = ObitImageRef(image);
  in->myErr   = ObitErrRef(err);

  /* Open, read image */
  ObitImageOpen (image, OBIT_IO_ReadOnly, err);
  ObitImageRead (image, NULL, err);
  if (err->error) Obit_traceback_msg (err, routine, image->name);
  
  /* Give warning if image rotated */
  if (image->myDesc->crota[1]!=0.0) {
    Obit_log_error(err, OBIT_InfoWarn, "%s: Rotated image not handled well", 
		   routine);
  }

  /* Get inputs from info */
  xtick = 0.0;
  ObitInfoListGetTest(in->info, "XTICK", &type, dim, &xtick);
  nxsub = 0;
  ObitInfoListGetTest(in->info, "NXSUB", &type, dim, &nxsub);
  ytick = 0.0;
  ObitInfoListGetTest(in->info, "YTICK", &type, dim, &ytick);
  nysub = 0;
  ObitInfoListGetTest(in->info, "NYSUB", &type, dim, &nysub);
  /* Sizes */
  csize = 1;
  ObitInfoListGetTest(in->info, "CSIZE", &type, dim, &csize);
  lwidth = 1;
  ObitInfoListGetTest(in->info, "LWIDTH", &type, dim, &lwidth);
  /* Adjust character size */
  cscale = (ofloat)csize;
  cscale = MAX (1.0, cscale);
  ObitPlotSetCharSize (in, cscale, err);

  /* set line width */
  lwidth = MAX (1, lwidth);
  ObitPlotSetLineWidth (in, lwidth, err);

  /* Define size of plot */
  ymin = image->myDesc->crval[1] + 
    (1.0 - image->myDesc->crpix[1]) * image->myDesc->cdelt[1];
  ymax = image->myDesc->crval[1] + 
    ((ofloat)image->myDesc->inaxes[1] - image->myDesc->crpix[1]) *  
    image->myDesc->cdelt[1];
  xmin = image->myDesc->crval[0] + 
    (1.0 - image->myDesc->crpix[0]) * image->myDesc->cdelt[0];
  xmax = image->myDesc->crval[0] + 
    ((ofloat)image->myDesc->inaxes[0] - image->myDesc->crpix[0]) *  
    image->myDesc->cdelt[0];
  /* convert to asec */
  xmin *= 3600.0;
  xmax *= 3600.0;
  ymin *= 3600.0;
  ymax *= 3600.0;

  /* set plotting area */
  axis = -2;
  just = 1;
  ObitPlotSetPlot (in, xmin, xmax, ymin, ymax, just, axis, err);

  /* Adjust for cos dec */
  xmin /=  (15.0);
  xmax /=  (15.0);
  cpgswin ((float)xmin, (float)xmax, (float)ymin, (float)ymax);

  /*  Plot labels */
  if (image->myDesc->equinox<1975.0) pre = 'B';
  else pre = 'J';
  if (fabs(image->myDesc->crota[1])==0.0)
    g_snprintf (xlabel,120,"Right Ascension (%c%6.1f)", 
		pre, image->myDesc->equinox); 
  else /* Don't call it RA */
    g_snprintf (xlabel,120,"Rotated Right Ascension (%c%6.1f)", 
		pre, image->myDesc->equinox); 
  if (fabs(image->myDesc->crota[1])==0.0)
    g_snprintf (ylabel,120,"Declination (%c%6.1f)", 
		pre, image->myDesc->equinox); 
  else  /* Don't call it Dec */
    g_snprintf (ylabel,120,"Rotated Declination (%c%6.1f)", 
		pre, image->myDesc->equinox); 
  ObitPlotLabel (in, xlabel, ylabel, label, err);

  /* set ticks, labels */
  cpgtbox ((const char*)xopt_def, (float)xtick, (int)nxsub, (const char*)yopt_def, (float)ytick, (int)nysub);

  /* Coordinate transformation */
  dx = image->myDesc->cdelt[0] * 3600.0 / (15.0);
  dy = image->myDesc->cdelt[1] * 3600.0;
  cr = cos (image->myDesc->crota[1]* 0.017453293);
  sr = sin (image->myDesc->crota[1]* 0.017453293);
  cr = 1.0;
  sr = 0.0;  /* Fooey - pgpplot blows this */
  tr[0] = xmin;
  tr[0] -= dx;  /* pgplot indexing error */
  tr[1] = dx*cr;
  tr[2] = dy*sr;
  tr[3] = ymin;
  tr[3] -= dy;  /* pgplot indexing error */
  tr[4] = dx*sr;
  tr[5] = dy*cr;

  /* Set contour levels */
  maxval = ObitFArrayMax (image->image, pos);
  minval = ObitFArrayMin (image->image, pos);
  if (lev <= 0.0) lev = 0.1 * MAX (fabs(maxval), fabs(minval));
  if (cntfac <= 0.0) cntfac = sqrt(2);
  if (maxval>0.0) {
    npc = 0.999 + log(maxval/lev) / log(cntfac);
    npc = MAX (0, npc);
  } else  npc = 0;
  if (minval<0.0) {
    nnc = 0.999  + log(-minval/lev) / log(cntfac);
    nnc = MAX (0, nnc);
  } else nnc = 0;
  levs = g_malloc0(MAX (npc, nnc)*sizeof(float));
  nx = image->myDesc->inaxes[0];
  ny = image->myDesc->inaxes[1];

  /* Set info line */
  g_snprintf (xlabel,120,"Levs=%g*%6.2f**n, Peak=%g %s, rot=%6.1f", 
	      lev, cntfac, maxval, image->myDesc->bunit, 
	      image->myDesc->crota[1]);
  cpgmtxt ("T", 0.2, 0.0, 0.0, xlabel);

  /* Init - use buffering */
  cpgbbuf();

  /* To pgplot data type */
  fmap = g_malloc0(image->image->arraySize*sizeof(float));
  for (i=0; i<image->image->arraySize; i++) fmap[i] = (float)image->image->array[i];

  /* Plot positive contours */
  if (npc>0) {
    levs[0] = lev;
    for (i=1; i<npc; i++) levs[i] = levs[i-1]*cntfac;
    cpgsls (1);  /* Solid lines */
    cpgconb (fmap, (int)nx, (int)ny, 1, (int)nx, 1, (int)ny, levs, 
	     (int)npc, tr, (float)fblank);
  }

  /* Plot negative contours */
  if (npc>0) {
    levs[0] = -lev;
    for (i=1; i<nnc; i++) levs[i] = levs[i-1]*cntfac;
    cpgsls (2); /* dashed lines */
    cpgconb (fmap, (int)nx, (int)ny, 1, (int)nx, 1, (int)ny, levs, 
	     (int)nnc, tr, (float)fblank);
    cpgsls (1); /* back to solid lines */
 }

  if (levs) g_free(levs); /* deallocate levels array */
  if (fmap) g_free(fmap); /* deallocate fmap array */

  /* Flush the buffer */
  cpgebuf();
  cpgupdt();

  /* Close image */
  ObitImageClose (image, err);
  if (err->error) Obit_traceback_msg (err, routine, image->name);
  /* Free image buffer */
  image->image = ObitFArrayUnref(image->image);
#endif /* HAVE_PGPLOT */

  /* Complain if plotting not available */
  if (!OK) Obit_log_error(err, OBIT_Error, "No plotting package available");
} /* end  ObitPlotContour */

/**
 * Gray scale plot of an ObitImage
 * Plot should be finalized with ObitPlotFinishPlot after all 
 * drawing on the current frame is finished.
 * \param in      Pointer to existing ObitPlot object.
 * \param label   Label for plot
 * \param image   Image to plot (first plane in BLC,TRC)
 *                Rotated images aren't done quite right
 * \param err ObitErr error stack
 *
 * Optional parameters on in->info
 * \li XTICK   (float) world coordinate interval between major tick marks
 *                     on X axis. If xtick=0.0 [def], the interval is chosen.
 * \li NXSUB   (int)   the number of subintervals to divide the major
 *                     coordinate interval into. If xtick=0.0 or nxsub=0,
 *                     the number is chosen. [def 0]
 * \li YTICK  (float)  like xtick for the Y axis.
 * \li NYSUB  (int)    like nxsub for the Y axis
 * \li CSIZE  (int)    Scaling factor for characters(default = 1)
 * \li SQRT   (bool)   If present and true plot sqrt (pixel_value)
 * \li INVERT (bool)   If present and true ionvert colors
 * \li COLOR  (string) Color scheme "GRAY", CONTOUR", "PHLAME"
 *                     default "GRAY"
 * \li PIX_MAX (float) maximum pixel value [def min in image]
 * \li PIX_MIN (float) minimum pixel value [def max in image]
 */
void ObitPlotGrayScale (ObitPlot* in, gchar *label, ObitImage *image,
			ObitErr *err)
{
  gboolean OK=FALSE;  /* Have an implementation? */
  gchar xlabel[121], ylabel[121];
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ofloat xtick, ytick, pixmax=-1.0e20, pixmin=-1.0e20;
  gboolean doSQRT = FALSE, doINVERT = FALSE;
  olong nxsub, nysub;
  ObitInfoType type;
  gchar color[128];
  /* Color contour (lifted from AIPS) tables */
  olong nlevel_CC = 128;
  unsigned char bCC_tab[128]=         /*blue table  */
    {0,15,15,15,15,15,15,15, 15,15,15,15,15,15,15,15,   
     72,72,72,72,72,72,72, 72,72,72,72,72,72,72,
     127,127,127,127,127,127,127, 127,127,127,127,127,127,127,  
     203,203,203,203,203,203,203, 203,203,203,203,203,203,203,
     0,0,0,0,0,0,0, 0,0,0,0,0,0,0,
     0,0,0,0,0,0,0, 0,0,0,0,0,0,0,
     0,0,0,0,0,0,0, 0,0,0,0,0,0,0, 
     0,0,0,0,0,0,0, 0,0,0,0,0,0,0,
     0,0,0,0,0,0,0, 0,0,0,0,0,0,0};
  unsigned char gCC_tab[128]=         /* green table  */
    {0,15,15,15,15,15,15,15, 15,15,15,15,15,15,15,15, 
     0,0,0,0,0,0,0,  0,0,0,0,0,0,0,
     0,0,0,0,0,0,0,  0,0,0,0,0,0,0,
     76,76,76,76,76,76,76, 76,76,76,76,76,76,76,
     59,59,59,59,59,59,59,  59,59,59,59,59,59,59,  
     229,229,229,229,229,229,229, 229,229,229,229,229,229,229,
     255,255,255,255,255,255,255, 255,255,255,255,255,255,255, 
     89,89,89,89,89,89,89, 89,89,89,89,89,89,89, 
     0,0,0,0,0,0,0,  0,0,0,0,0,0,0};
  unsigned char rCC_tab[128]=          /*red table  */
    {0,15,15,15,15,15,15,15,  15, 15,15,15,15,15,15,15,  
     36,36,36,36,36,36,36, 36,36,36,36,36,36,36,
     0,0,0,0,0,0,0, 0,0,0,0,0,0,0, 
     15,15,15,15,15,15,15, 15,15,15,15,15,15,15,
     0,0,0,0,0,0,0, 0,0,0,0,0,0,0, 
     0,0,0,0,0,0,0, 0,0,0,0,0,0,0,
     255,255,255,255,255,255,255, 255,255,255,255,255,255,255,  
     255,255,255,255,255,255,255, 255,255,255,255,255,255,255,
     255,255,255,255,255,255,255, 255,255,255,255,255,255,255};
  /* PHLAME (lifted from AIPS) tables */
  olong nlevel_PH = 128;
  unsigned char bPH_tab[]=         /* green table  */
    {0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
     0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
     0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
     0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
     0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
     0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
     0,  26,  53,  66,  79,  89,  99, 108, 117, 115, 133, 140, 147, 154, 161, 167,
     174, 180, 186, 191, 197, 202, 208, 212, 219, 224, 229, 233, 238, 242, 248, 252};
  unsigned char gPH_tab[]=         /* green table  */
    {0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
     0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
     0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,   0,
     0,   0,   0,   0,   0,   0,   0,  31,   0,   0,   0,   0,   0,   0,   0,  15,
     31,  38,  46,  51,  57,  62,  68,  72,  77,  81,  86,  89,  93,  97, 101, 104,
     107, 110, 114, 117, 120, 123, 127, 130, 133, 135, 138, 141, 144, 146, 149, 151,
     154, 156, 159, 162, 165, 167, 170, 172, 175, 179, 182, 185, 189, 192, 195, 198,
     202, 205, 209, 212, 215, 218, 222, 225, 228, 231, 234, 237, 240, 243, 246, 251};
  unsigned char rPH_tab[]=          /*red table  */
    {0,   0,   0,   0,  0,   0,    0,   0,   0,   0,   0,   0,   0,   0,   0,  15,
     29,  36,  43,  49, 55,  59,   64,  68,  73,  77,  81,  84,  88,  91,  95,  98, 
     102, 105, 109, 112, 114, 117, 120, 123, 126, 128, 131, 134, 137, 139, 142, 144,
     147, 149, 152, 154, 156, 158, 161, 163, 165, 167, 170, 172, 174, 176, 178, 180,
     183, 185, 187, 189, 191, 193, 195, 196, 198, 200, 202, 204, 207, 209, 211, 212,
     214, 216, 218, 219, 221, 223, 225, 226, 228, 230, 232, 231, 235, 236, 238, 240,
     242, 244, 246, 247, 248, 249, 251, 253, 255, 255, 255, 255, 255, 255, 255, 255,
     255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255, 255};
  gchar *routine="ObitPlotGrayScale";

/****************** plplot implementation *************************/
#ifdef HAVE_PLPLOT  /* Only if plplot available */
  ObitImageDesc *id;
  olong i, just, axis;
  ofloat xmax, xmin, ymax, ymin, cscale;
  PLFLT px, py, wx, wy, **map=NULL, clevel[256], delta;
  ofloat maxval=0.0, minval=0.0;
  olong pos[2];
  olong csize, nx, ny;
  gchar rast[16], decst[16];
  gchar units[20], *xopt_def="bcnts", *yopt_def="bcnts", pre;
  PLINT nlevel, red[255], green[255], blue[255], it;


 /* error checks */
  if (err->error) return;
  g_assert (ObitPlotIsA(in));
  g_assert (ObitImageIsA(image));
  
  OK = TRUE;  /* Have a plotting package */
  in->xLog = FALSE;
  in->yLog = FALSE;

  /* Save info */
  in->myImage = ObitImageRef(image);
  in->myErr   = ObitErrRef(err);

  /* Open, read image */
  ObitImageOpen (image, OBIT_IO_ReadOnly, err);
  ObitImageRead (image, NULL, err);
  if (err->error) Obit_traceback_msg (err, routine, image->name);
  
  /* Give warning if image rotated */
  if (image->myDesc->crota[1]!=0.0) {
    Obit_log_error(err, OBIT_InfoWarn, "%s: Rotated image not handled well", 
		   routine);
  }

  /* Get inputs from info */
  xtick = 0.0;
  ObitInfoListGetTest(in->info, "XTICK", &type, dim, &xtick);
  nxsub = 0;
  ObitInfoListGetTest(in->info, "NXSUB", &type, dim, &nxsub);
  ytick = 0.0;
  ObitInfoListGetTest(in->info, "YTICK", &type, dim, &ytick);
  nysub = 0;
  ObitInfoListGetTest(in->info, "NYSUB", &type, dim, &nysub);
  /* Sizes */
  csize = 1;
  ObitInfoListGetTest(in->info, "CSIZE", &type, dim, &csize);
  ObitInfoListGetTest(in->info, "SQRT",  &type, dim, &doSQRT);
  ObitInfoListGetTest(in->info, "INVERT",  &type, dim, &doINVERT);
  ObitInfoListGetTest(in->info, "PIX_MAX", &type, dim, &pixmax); 
  ObitInfoListGetTest(in->info, "PIX_MIN", &type, dim, &pixmin); 
  sprintf (color, "GRAY");
  ObitInfoListGetTest(in->info, "COLOR", &type, dim, color); 

 /* Adjust character size */
  cscale = (ofloat)csize;
  cscale = MAX (1.0, cscale);
  ObitPlotSetCharSize (in, cscale, err);

  /* What units/scaling for axes? */
  id = image->myDesc;
  if (id->inaxes[0]*fabs(id->cdelt[0])>1.0) {
    /* degrees */
    strcpy (units,"deg");
    in->scalex = 1.0;
    in->scaley = 1.0;
  } else if (id->inaxes[0]*fabs(id->cdelt[0])>0.016667) {
    /* arcmin */
    strcpy (units,"amin");
    in->scalex = 60.0;
    in->scaley = 60.0;
  } else if (id->inaxes[0]*fabs(id->cdelt[0])>0.0027778) {
    /* arcsec */
    strcpy (units,"asec");
    in->scalex = 3600.0;
    in->scaley = 3600.0;
  } else {
    /* milliarcsec  */
    strcpy (units,"mas");
    in->scalex = 3600000.0;
    in->scaley = 3600000.0;
  }

  /* Correct X scaling for the artificial effect of declination */
  if (id->jlocd>=0) 
    in->scalex *= cos(DG2RAD*id->crval[id->jlocd]);

  nx = image->myDesc->inaxes[0];
  ny = image->myDesc->inaxes[1];
  
  /* Define size of plot in world coordinates */
  px = (PLFLT)1.0; py = (PLFLT)1.0;
  plplotCoord (px, py, &wx, &wy, (PLPointer)in);
  xmin = (ofloat)wx;
  ymin = (ofloat)wy;
  px = (PLFLT)nx; py = (PLFLT)ny;
  plplotCoord (px, py, &wx, &wy, (PLPointer)in);
  xmax = (ofloat)wx;
  ymax = (ofloat)wy;

  /* set plotting area */
  axis = -2;
  just = 1;
  ObitPlotSetPlot (in, xmin, xmax, ymin, ymax, just, axis, err);

  /*  Plot labels */
  if (image->myDesc->equinox<1975.0) pre = 'B';
  else pre = 'J';
  if (fabs(image->myDesc->crota[1])==0.0)
    g_snprintf (xlabel,120,"Right Ascension (%s)", units);
  else /* Don't call it RA */
    g_snprintf (xlabel,120,"Rotated Right Ascension (%s)", units);
  if (fabs(image->myDesc->crota[1])==0.0)
    g_snprintf (ylabel,120,"Declination (%s)", units);
  else  /* Don't call it Dec */
    g_snprintf (ylabel,120,"Rotated Declination (%s)", units);
  ObitPlotLabel (in, xlabel, ylabel, label, err);

  /* Get center position string */
  ObitPosLabelUtilRA2HMS  (image->myDesc->crval[0], image->myDesc->ctype[0], rast);
  ObitPosLabelUtilDec2DMS (image->myDesc->crval[1], image->myDesc->ctype[1], decst);

  /* set ticks, labels */
  ObitPlotDrawAxes (in, xopt_def, xtick, nxsub, yopt_def, ytick, nysub, err);

  /* Set pixel range */
  if (pixmax<-1.0e19) maxval = ObitFArrayMax (image->image, pos);
  else maxval = pixmax;
  if (pixmin<-1.0e19) minval = ObitFArrayMin (image->image, pos);
  else minval = pixmin;

  /* Set plotting info line */
  if (fabs(image->myDesc->crota[1])>0.001)
    if (doSQRT) {
      g_snprintf (xlabel,120,"min=%g, max=%g %s, sqrt stretch, rot=%6.1f", 
		  minval, maxval, image->myDesc->bunit, 
		  image->myDesc->crota[1]);
    } else {
      g_snprintf (xlabel,120,"min=%g, max=%g %s, rot=%6.1f", 
		  minval, maxval, image->myDesc->bunit, 
		  image->myDesc->crota[1]);
    }
  else
    if (doSQRT) {
      g_snprintf (xlabel,120,"min=%g, max=%g %s, sqrt stretch",
		  minval, maxval, image->myDesc->bunit);
    } else {
      g_snprintf (xlabel,120,"min=%g, max=%g %s",
		  minval, maxval, image->myDesc->bunit);
    }

  just = 0;
  ObitPlotSetCharSize (in, 0.5, err);
  ObitPlotRelText (in, "B", 8.5, 0.0, just, xlabel, err);
  /* Center position */
  g_snprintf (xlabel,120,"Center = %s, %s %c%6.1f", 
	      rast, decst, pre, image->myDesc->equinox);
  ObitPlotRelText (in, "T", 1.0, 0.0, just, xlabel, err);

  /* Square root? */
  if (doSQRT) {
    /* Modify image buffer */
    ObitFArraySAdd(image->image, -minval+1.0e-20);
    ObitFArraySqrt (image->image);
    maxval = sqrt(maxval - minval);
    minval = 0.0;
  }

  /* Clip range to max/min */
  ObitFArrayInClip (image->image, -1.0e25, minval, minval);
  ObitFArrayInClip (image->image, maxval, 1.0e25, maxval);
  /* blank below minval */
  ObitFArrayDeblank (image->image, minval-1.0);

  /* Convert Obit image to plplot array */
  map = reorderPixels(image->image);

  /* Color map */
  if (strncmp(color, "CONTOUR", 7)==0) {  /* Color contour */
    nlevel =  (PLINT)nlevel_CC;
    for (i=0; i<nlevel; i++) {
      red[i]   = (PLINT)rCC_tab[i];
      green[i] = (PLINT)gCC_tab[i];
      blue[i]  = (PLINT)bCC_tab[i];
    }
  } else if (strncmp(color, "PHLAME", 6)==0) { /* Red-yellow Phlame */
    nlevel =  (PLINT)nlevel_PH;
    for (i=0; i<nlevel; i++) {
      red[i]   = (PLINT)rPH_tab[i];
      green[i] = (PLINT)gPH_tab[i];
      blue[i]  = (PLINT)bPH_tab[i];
    }
  } else { /* default 255 level gray scale */
    nlevel = (PLINT)255;
    for (i=0; i<nlevel; i++) {red[i] = green[i] = blue[i] = (PLINT)i;}
  }

  /* Invert colors? keep first level */
  if (doINVERT) { 
    for (i=1; i<nlevel/2; i++) {
      it = red[i];   red[i]   = red[nlevel-i-1];   red[nlevel-i-1]   = it;
      it = green[i]; green[i] = green[nlevel-i-1]; green[nlevel-i-1] = it;
      it = blue[i];  blue[i]  = blue[nlevel-i-1];  blue[nlevel-i-1]  = it;
    }
  }
  plscmap1(red, green, blue, nlevel);

  /* Levels - reserve one for blanked */
  delta = (maxval-minval) / (nlevel-1);
  for (i=0; i<=nlevel; i++) clevel[i] = minval + (i-1)*delta;;

  /* Plot  */
  plshades (map, (PLINT)nx, (PLINT)ny, NULL, 
	   (PLFLT)xmin, (PLFLT)xmax, (PLFLT)ymin, (PLFLT)ymax, 
	    clevel, nlevel+1, 0, 0, 0,
	    plfill, (PLBOOL)0, 
	    plplotCoord, (PLPointer)in);
	   

  /* Deallocate plplot array */
  if (map) {
    for (i=0; i<nx; i++) g_free(map[i]);
    g_free(map);
  }

  /* Flush the buffer */
  plflush();

  /* Close image */
  ObitImageClose (image, err);
  if (err->error) Obit_traceback_msg (err, routine, image->name);
  /* Free image buffer */
  image->image = ObitFArrayUnref(image->image);
#endif /* HAVE_PLPLOT */

/****************** pgplot implementation *************************/
#ifdef HAVE_PGPLOT  /* Only if pgplot available */
  olong i, just, axis, npc, nnc;
  ofloat xmax, xmin, ymax, ymin, cscale, maxval, minval, fblank = ObitMagicF();
  float tr[6];
  float *levs=NULL, *fmap=NULL;
  olong pos[2];
  ofloat dx, dy, cr, sr;
  olong csize, lwidth, nx, ny;
  gchar *xopt_def="BCNTSZYH", *yopt_def="BCNTSZYD", pre;


 /* error checks */
  if (err->error) return;
  g_assert (ObitPlotIsA(in));
  g_assert (ObitImageIsA(image));
  
  OK = TRUE;  /* Have a plotting package */
  in->xLog = FALSE;
  in->yLog = FALSE;

  /* Save info */
  in->myImage = ObitImageRef(image);
  in->myErr   = ObitErrRef(err);

  /* Open, read image */
  ObitImageOpen (image, OBIT_IO_ReadOnly, err);
  ObitImageRead (image, NULL, err);
  if (err->error) Obit_traceback_msg (err, routine, image->name);
  
  /* Give warning if image rotated */
  if (image->myDesc->crota[1]!=0.0) {
    Obit_log_error(err, OBIT_InfoWarn, "%s: Rotated image not handled well", 
		   routine);
  }

  /* Get inputs from info */
  xtick = 0.0;
  ObitInfoListGetTest(in->info, "XTICK", &type, dim, &xtick);
  nxsub = 0;
  ObitInfoListGetTest(in->info, "NXSUB", &type, dim, &nxsub);
  ytick = 0.0;
  ObitInfoListGetTest(in->info, "YTICK", &type, dim, &ytick);
  nysub = 0;
  ObitInfoListGetTest(in->info, "NYSUB", &type, dim, &nysub);
  /* Sizes */
  csize = 1;
  ObitInfoListGetTest(in->info, "CSIZE", &type, dim, &csize);
  lwidth = 1;
  ObitInfoListGetTest(in->info, "LWIDTH", &type, dim, &lwidth);
  ObitInfoListGetTest(in->info, "SQRT",   &type, dim, &doSQRT);
  ObitInfoListGetTest(in->info, "INVERT",  &type, dim, &doINVERT);
  ObitInfoListGetTest(in->info, "PIX_MAX", &type, dim, &pixmax); 
  ObitInfoListGetTest(in->info, "PIX_MIN", &type, dim, &pixmin); 

  /* Adjust character size */
  cscale = (ofloat)csize;
  cscale = MAX (1.0, cscale);
  ObitPlotSetCharSize (in, cscale, err);

  /* set line width */
  lwidth = MAX (1, lwidth);
  ObitPlotSetLineWidth (in, lwidth, err);

  /* Define size of plot */
  ymin = image->myDesc->crval[1] + 
    (1.0 - image->myDesc->crpix[1]) * image->myDesc->cdelt[1];
  ymax = image->myDesc->crval[1] + 
    ((ofloat)image->myDesc->inaxes[1] - image->myDesc->crpix[1]) *  
    image->myDesc->cdelt[1];
  xmin = image->myDesc->crval[0] + 
    (1.0 - image->myDesc->crpix[0]) * image->myDesc->cdelt[0];
  xmax = image->myDesc->crval[0] + 
    ((ofloat)image->myDesc->inaxes[0] - image->myDesc->crpix[0]) *  
    image->myDesc->cdelt[0];
  /* convert to asec */
  xmin *= 3600.0;
  xmax *= 3600.0;
  ymin *= 3600.0;
  ymax *= 3600.0;

  /* set plotting area */
  axis = -2;
  just = 1;
  ObitPlotSetPlot (in, xmin, xmax, ymin, ymax, just, axis, err);

  /* Adjust for cos dec */
  xmin /=  (15.0);
  xmax /=  (15.0);
  cpgswin ((float)xmin, (float)xmax, (float)ymin, (float)ymax);

  /*  Plot labels */
  if (image->myDesc->equinox<1975.0) pre = 'B';
  else pre = 'J';
  if (fabs(image->myDesc->crota[1])==0.0)
    g_snprintf (xlabel,120,"Right Ascension (%c%6.1f)", 
		pre, image->myDesc->equinox); 
  else /* Don't call it RA */
    g_snprintf (xlabel,120,"Rotated Right Ascension (%c%6.1f)", 
		pre, image->myDesc->equinox); 
  if (fabs(image->myDesc->crota[1])==0.0)
    g_snprintf (ylabel,120,"Declination (%c%6.1f)", 
		pre, image->myDesc->equinox); 
  else  /* Don't call it Dec */
    g_snprintf (ylabel,120,"Rotated Declination (%c%6.1f)", 
		pre, image->myDesc->equinox); 
  ObitPlotLabel (in, xlabel, ylabel, label, err);

  /* set ticks, labels */
  cpgtbox ((const char*)xopt_def, (float)xtick, (int)nxsub, (const char*)yopt_def, (float)ytick, (int)nysub);

  /* Coordinate transformation */
  dx = image->myDesc->cdelt[0] * 3600.0 / (15.0);
  dy = image->myDesc->cdelt[1] * 3600.0;
  cr = cos (image->myDesc->crota[1]* 0.017453293);
  sr = sin (image->myDesc->crota[1]* 0.017453293);
  cr = 1.0;
  sr = 0.0;  /* Fooey - pgpplot blows this */
  tr[0] = xmin;
  tr[0] -= dx;  /* pgplot indexing error */
  tr[1] = dx*cr;
  tr[2] = dy*sr;
  tr[3] = ymin;
  tr[3] -= dy;  /* pgplot indexing error */
  tr[4] = dx*sr;
  tr[5] = dy*cr;

  /* Set pixel range */
  if (maxval<1.0e19) maxval = ObitFArrayMax (image->image, pos);
  else maxval = pixmax;
  if (minval<1.0e19) minval = ObitFArrayMin (image->image, pos);
  else minval = pixmin;


  nx = image->myDesc->inaxes[0];
  ny = image->myDesc->inaxes[1];

  /* Set info line */
  if (fabs(image->myDesc->crota[1])>0.001) {
    if (doSQRT) {
      g_snprintf (xlabel,120,"min=%g, max=%g %s, sqrt stretch, rot=%6.1f", 
		  minval, maxval, image->myDesc->bunit, 
		  image->myDesc->crota[1]);
    } else {
      g_snprintf (xlabel,120,"min=%g, max=%g %s, rot=%6.1f", 
		  minval, maxval, image->myDesc->bunit, 
		  image->myDesc->crota[1]);
    }
  } else {
    if (doSQRT) {
      g_snprintf (xlabel,120,"min=%g, max=%g %s, sqrt stretch",
		  minval, maxval, image->myDesc->bunit);
    } else {
      g_snprintf (xlabel,120,"min=%g, max=%g %s",
		  minval, maxval, image->myDesc->bunit);
    }
  }
  cpgmtxt ("T", 0.2, 0.0, 0.0, xlabel);

  /* Init - use buffering */
  cpgbbuf();

  /* Square root? */
  if (doSQRT) {
    /* Modify image buffer */
    ObitFArraySAdd(image->image, -minval+1.0e-20);
    ObitFArraySqrt (image->image);
    maxval = sqrt(maxval - minval);
    minval = 0.0;
  }

  /* Clip range to max/min */
  ObitFArrayInClip (image->image, -1.0e25, minval, minval);
  ObitFArrayInClip (image->image, maxval, 1.0e25, maxval);
  /* blank below minval */
  ObitFArrayDeblank (image->image, minval-1.0);

  /* To pgplot data type */
  fmap = g_malloc0(image->image->arraySize*sizeof(float));
  for (i=0; i<image->image->arraySize; i++) fmap[i] = (float)image->image->array[i];

  /* Plot */
  cpggray (fmap, (int)nx, (int)ny, 1, (int)nx, 1, (int)ny, 
	   (float)maxval, (float)minval, tr);

  if (fmap) g_free(fmap); /* deallocate fmap array */

  /* Flush the buffer */
  cpgebuf();
  cpgupdt();

  /* Close image */
  ObitImageClose (image, err);
  if (err->error) Obit_traceback_msg (err, routine, image->name);
  /* Free image buffer */
  image->image = ObitFArrayUnref(image->image);
#endif /* HAVE_PGPLOT */

  /* Complain if plotting not available */
  if (!OK) Obit_log_error(err, OBIT_Error, "No plotting package available");
} /* end  ObitPlotGrayScale */

/**
 * X-Y Plot with error bars in Y
 * Plot should be finalized with ObitPlotFinishPlot after all 
 * drawing on the current frame is finished.
 * \param in      Pointer to Plot object.
 * \param in      Pointer to Plot object.
 * \param symbol  Symbol index to use for plotting 
 *                values in the range [1,12] are usable 
 *                if negative, use abs value and connect points
 * \li 1 = dot
 * \li 2 = plus
 * \li 3 = *
 * \li 4 = open circle
 * \li 5 = x
 * \li 6 = open square
 * \li 7 = open triangle
 * \li 8 = open star
 * \li 9 = filled triangle
 * \li 10 = filled square
 * \li 11 = filled circle
 * \li 12 = filled star
 *
 * \param n       Number of data points in x, y
 * \param x       Independent variable, if NULL use index
 * \param y       Dependent variable
 * \param e       if nonNULL, error in y
 * \param err ObitErr error stack
 *
 * Optional parameters on in->info
 * \li XMAX (float) maximum X value (defaults to actual value)
 * \li XMIN (float) minimum X value (defaults to actual value)
 * \li YMAX (float) maximum Y value (defaults to actual value)
 * \li YMIN (float) minimum Y value (defaults to actual value)
 * \li TITLE (string)  Label for the plot (defaults to none), max 120
 * \li XLABEL (string) Label for horizontal axis (defaults to none)
 * \li XOPT   (string) Options for horizontal axis (default "BCNTS")
 *                     See #ObitPlotDrawAxes for details.
 * \li YLABEL (string) Label for vertical axis (defaults to none)
 * \li YOPT   (string) Options for  vertical axis (default "BCNTS")
 *                     See #ObitPlotDrawAxes for details.
 * \li XTICK   (float) world coordinate interval between major tick marks
 *                     on X axis. If xtick=0.0 [def], the interval is chosen.
 * \li NXSUB   (int)   the number of subintervals to divide the major
 *                     coordinate interval into. If xtick=0.0 or nxsub=0,
 *                     the number is chosen. [def 0]
 * \li YTICK  (float)  like xtick for the Y axis.
 * \li NYSUB  (int)    like nxsub for the Y axis
 * \li CSIZE  (int)    Scaling factor for characters(default = 1)
 * \li SSIZE  (int)    Scaling factor for symbols(default = 1)
 * \li LWIDTH (int)    Line width (default = 1)
 */
void ObitPlotXYErr (ObitPlot* in, olong symbol,
		    olong n, ofloat *x, ofloat *y, ofloat *e, 
		    ObitErr *err)
{
  gboolean OK=FALSE;  /* Have an implementation? */
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
  gboolean gotRange, gotXMAX, gotYMAX, gotXMIN, gotYMIN, doConnect;
  ofloat xtick, ytick;
  olong nxsub, nysub;
  gchar title[121], xlabel[121], ylabel[121], xopt[31], yopt[31];
  gchar *routine="ObitPlotXYErr";

/****************** plplot implementation *************************/
#ifdef HAVE_PLPLOT  /* Only if plplot available */
  olong i, just, axis;
  ofloat xmax, xmin, ymax, ymin, cscale, fblank = ObitMagicF();
  olong csize, ssize, lwidth, lsymbol;
  olong syms[12] = {210,225,228,840,227,841,842,844,852,851,850,856};
  gchar *optemp, *xopt_def="BCNTS", *yopt_def="BCNTS";
  PLFLT sscale, *xp, *yp, *xx, *ply=NULL, *plx=NULL, *logx=NULL, *logy=NULL;

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (y!=NULL);

  /* anything to do? */
  if (n<=0) return;
  
  OK = TRUE;  /* Have a plotting package */
  lsymbol = MAX(1, MIN(12, abs(symbol)));
  lsymbol = syms[lsymbol-1];  /* Find plplot value */
  doConnect = symbol<0;
  in->xLog = FALSE;
  in->yLog = FALSE;

  /* using index for x? */
  if (x!=0) xx = x;
  else { /* yes */
    xx = g_malloc0(n*sizeof(PLFLT));
    for (i=0; i<n; i++) xx[i] = (PLFLT)i;
  }

  /* Get inputs from info */
  /* Plot title */
  for (i=0; i<121; i++) title[i] = 0;
  strcpy (title, "    ");
  ObitInfoListGetTest(in->info, "TITLE", &type, dim, title);
 
  /* Horizontal label */
  for (i=0; i<121; i++) xlabel[i] = 0;
  strcpy (xlabel, "    ");
  ObitInfoListGetTest(in->info, "XLABEL", &type, dim, xlabel);
 
  /* Horizontal label options */
  for (i=0; i<31; i++) xopt[i] = 0;
  strcpy (xopt, xopt_def);
  ObitInfoListGetTest(in->info, "XOPT", &type, dim, xopt);
  for (i=0; i<30; i++) if(xopt[i]=='L') in->xLog = TRUE;  /* Log plot? */
  optemp = g_ascii_strdown (xopt,-1);  /* Convert to lower case */
  strcpy (xopt, optemp);
  g_free(optemp);
  xtick = 0.0;
  ObitInfoListGetTest(in->info, "XTICK", &type, dim, &xtick);
  nxsub = 0;
  ObitInfoListGetTest(in->info, "NXSUB", &type, dim, &nxsub);

  /* Vertical label */
  for (i=0; i<121; i++) ylabel[i] = 0;
  strcpy (ylabel, "    ");
  ObitInfoListGetTest(in->info, "YLABEL", &type, dim, ylabel);

  /* Vertical label options */
  for (i=0; i<31; i++) yopt[i] = 0;
  strcpy (yopt, yopt_def);
  ObitInfoListGetTest(in->info, "YOPT", &type, dim, yopt);
  for (i=0; i<30; i++) if(yopt[i]=='L') in->yLog = TRUE;  /* Log plot? */
  optemp = g_ascii_strdown (yopt,-1);  /* Convert to lower case */
  strcpy (yopt, optemp);
  g_free(optemp);
  ytick = 0.0;
  ObitInfoListGetTest(in->info, "YTICK", &type, dim, &ytick);
  nysub = 0;
  ObitInfoListGetTest(in->info, "NYSUB", &type, dim, &nysub);
 
  /* Sizes */
  csize = 1;
  ObitInfoListGetTest(in->info, "CSIZE", &type, dim, (gpointer*)&csize);
  ssize = 1;
  ObitInfoListGetTest(in->info, "SSIZE", &type, dim, (gpointer*)&ssize);
  lwidth = 1;
  ObitInfoListGetTest(in->info, "LWIDTH", &type, dim, (gpointer*)&lwidth);
  just = 0;
  ObitInfoListGetTest(in->info, "JUST", &type, dim, (gpointer*)&just);
  just = (just!=0);
 
  /* Plot range */
  xmax = fblank;
  ObitInfoListGetTest(in->info, "XMAX", &type, dim, (gpointer*)&xmax);
  ymax = fblank;
  ObitInfoListGetTest(in->info, "YMAX", &type, dim, (gpointer*)&ymax);
  xmin = fblank;
  ObitInfoListGetTest(in->info, "XMIN", &type, dim, (gpointer*)&xmin);
  ymin = fblank;
  ObitInfoListGetTest(in->info, "YMIN", &type, dim, (gpointer*)&ymin);

  /* Check data for range? */
  gotXMAX = xmax != fblank;
  gotYMAX = ymax != fblank;
  gotXMIN = xmin != fblank;
  gotYMIN = ymin != fblank;
  gotRange =  gotXMAX && gotYMAX && gotXMIN && gotYMIN;
   if (!gotRange) { /* Need to check data? */
    for (i=0; i<n; i++) {
      if ((x!=NULL) && (x[i]!=fblank)) {
	if  (xmax == fblank) xmax = x[i];
	else xmax = MAX (x[i], xmax);
	if  (xmin == fblank) xmin = x[i];
	else xmin = MIN (x[i], xmin);
      } else if (x==NULL) {
	if  (xmax == fblank) xmax = xx[i];
	else xmax = MAX (xx[i], xmax);
	if  (xmin == fblank) xmin = xx[i];
	else xmin = MIN (xx[i], xmin);
      }
      if (y[i]!=fblank) {
	if  (ymax == fblank) ymax = y[i];
	else ymax = MAX (y[i], ymax);
	if  (ymin == fblank) ymin = y[i];
	else ymin = MIN (y[i], ymin);
      }
    }
  } /* end check range */

  /* Are the range values OK? */
  if (xmax<=xmin) { /* Nope - bug out */
    Obit_log_error(err, OBIT_Error, "%s: XMAX(%g) <= XMIN(%g)", 
		   routine, xmax, xmin);
    if (x==NULL) g_free(xx);
    return;
 }

  if (ymax<=ymin) { /* Nope - bug out */
    Obit_log_error(err, OBIT_Error, "%s: YMAX (%g) <= YMIN(%g)", 
		   routine, ymax, ymin);
    if (x==NULL) g_free(xx);
    return;
 }

  /* Log? */
  if (in->xLog) {
    /* Must be positive */
    if ((xmax<=0.0) || (xmin<=0.0)) {
      Obit_log_error(err, OBIT_Error, 
		     "%s: Values for log of x axis must be positive %f %f", 
		      routine, xmin, xmax);
      if (x==NULL) g_free(xx);
      return;
    }
    xmin = log10(xmin);
    xmax = log10(xmax);
    logx = g_malloc0(n*sizeof(PLFLT));
    for (i=0; i<n; i++) logx[i] = (PLFLT)log10(xx[i]);
    xp = logx;
  } else {
    plx = g_malloc0(n*sizeof(PLFLT));
    for (i=0; i<n; i++) plx[i] = (PLFLT)xx[i];
    xp = plx;
  }

  if (in->yLog) {
    /* Must be positive */
    if ((ymax<=0.0) || (ymin<=0.0)) {
      Obit_log_error(err, OBIT_Error, 
		     "%s: Values for log of y axis must be positive %f %f", 
		      routine, ymin, ymax);
      if (x==NULL) g_free(xx);
      return;
    }
    ymin = log10(ymin);
    ymax = log10(ymax);
    logy = g_malloc0(n*sizeof(ofloat));
    for (i=0; i<n; i++) {
      if (y[i]!=fblank) logy[i] = log10(y[i]);
      else logy[i] = xmax*10.0;  /* Don't plot if blanked */
    }
    yp = logy;
  } else {
    yp = y;
  }

  /* If autosetting the range, expand a bit */
  if (!gotXMAX) xmax += 0.05 * (xmax-xmin);
  if (!gotXMIN) xmin -= 0.05 * (xmax-xmin);
  if (!gotYMAX) ymax += 0.05 * (ymax-ymin);
  if (!gotYMIN) ymin -= 0.05 * (ymax-ymin);

  /* Adjust character size */
  cscale = (ofloat)csize;
  cscale = MAX (1.0, cscale);
  ObitPlotSetCharSize (in, cscale, err);

  /* set line width */
  lwidth = MAX (1, lwidth);
  ObitPlotSetLineWidth (in, lwidth, err);

  /* set plotting area */
  axis = -2;
  if (in->xLog  && !in->yLog) axis = 10;
  if (!in->xLog &&  in->yLog) axis = 20;
  if ( in->xLog &&  in->yLog) axis = 30;
  ObitPlotSetPlot (in, xmin, xmax, ymin, ymax, just, axis, err);

  /*  Plot label */
  ObitPlotLabel (in, xlabel, ylabel, title, err);

  /* set ticks */
  ObitPlotDrawAxes (in, xopt, xtick, nxsub, yopt, ytick, nysub, err);

  /* Connect the dots? */
  if (doConnect) plline((PLINT)n, xp, yp);

  /* Adjust symbol size */
  sscale = (PLFLT)ssize;
  sscale = MAX (1.0, sscale);
  c_plssym((PLFLT)0.0, sscale);

  /* Plot symbols */
  plsym ((PLINT)n, xp, yp, (PLINT)lsymbol);

  /* Error bars */
  if (e!=NULL) {
    for (i=0; i<n; i++) {
      ObitPlotDrawLine (in, xx[i], y[i]-e[i], xx[i], y[i]+e[i], err);
    }
  }

  /* Flush the buffer */
  plflush();

  /* deallocate arrays if allocated */
  if (x==NULL)     g_free(xx);
  if (logx!=NULL)  g_free(logx);
  if (logy!=NULL)  g_free(logy);
  if (plx!=NULL)   g_free(plx);
  if (ply!=NULL)   g_free(ply);

#endif /* HAVE_PLPLOT */

/****************** pgplot implementation *************************/
#ifdef HAVE_PGPLOT  /* Only if pgplot available */
  olong i, just, axis;
  float *ply=NULL, *plx=NULL;
  ofloat xmax, xmin, ymax, ymin, cscale, sscale, *xx, fblank = ObitMagicF();
  olong csize, ssize, lwidth;
  olong lsymbol, syms[12] = {210,225,228,840,227,841,842,844,852,851,850,856};
  gchar *xopt_def="BCNTS", *yopt_def="BCNTS";

  /* error checks */
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (y!=NULL);

  /* anything to do? */
  if (n<=0) return;
  
  OK = TRUE;  /* Have a plotting package */
  lsymbol = MAX(1, MIN(12, abs(symbol)));
  lsymbol = syms[lsymbol-1];  /* Find plplot value */
  doConnect = symbol<0;
  in->xLog = FALSE;
  in->yLog = FALSE;

  /* using index for x? */
  if (x!=0) xx = x;
  else { /* yes */
    xx = g_malloc0(n*sizeof(ofloat));
    for (i=0; i<n; i++) xx[i] = (float)i;
  }

  /* Get inputs from info */
  /* Plot title */
  for (i=0; i<121; i++) title[i] = 0;
  strcpy (title, "    ");
  ObitInfoListGetTest(in->info, "TITLE", &type, dim, title);
 
  /* Horizontal label */
  for (i=0; i<121; i++) xlabel[i] = 0;
  strcpy (xlabel, "    ");
  ObitInfoListGetTest(in->info, "XLABEL", &type, dim, xlabel);
 
  /* Horizontal label options */
  for (i=0; i<31; i++) xopt[i] = 0;
  strcpy (xopt, xopt_def);
  ObitInfoListGetTest(in->info, "XOPT", &type, dim, xopt);
  for (i=0; i<30; i++) if(xopt[i]=='L') in->xLog = TRUE;  /* Log plot? */
  xtick = 0.0;
  ObitInfoListGetTest(in->info, "XTICK", &type, dim, &xtick);
  nxsub = 0;
  ObitInfoListGetTest(in->info, "NXSUB", &type, dim, &nxsub);
  
  /* Vertical label */
  for (i=0; i<121; i++) ylabel[i] = 0;
  strcpy (ylabel, "    ");
  ObitInfoListGetTest(in->info, "YLABEL", &type, dim, ylabel);

  /* Vertical label options */
  for (i=0; i<31; i++) yopt[i] = 0;
  strcpy (yopt, yopt_def);
  ObitInfoListGetTest(in->info, "YOPT", &type, dim, yopt);
  for (i=0; i<30; i++) if(yopt[i]=='L') in->yLog = TRUE;  /* Log plot? */
  ytick = 0.0;
  ObitInfoListGetTest(in->info, "YTICK", &type, dim, &ytick);
  nysub = 0;
  ObitInfoListGetTest(in->info, "NYSUB", &type, dim, &nysub);
 
  /* Sizes */
  csize = 1;
  ObitInfoListGetTest(in->info, "CSIZE", &type, dim, (gpointer*)&csize);
  ssize = 1;
  ObitInfoListGetTest(in->info, "SSIZE", &type, dim, (gpointer*)&ssize);
  lwidth = 1;
  ObitInfoListGetTest(in->info, "LWIDTH", &type, dim, (gpointer*)&lwidth);
  just = 0;
  ObitInfoListGetTest(in->info, "JUST", &type, dim, (gpointer*)&just);
  just = (just!=0);
 
  /* Plot range */
  xmax = fblank;
  ObitInfoListGetTest(in->info, "XMAX", &type, dim, (gpointer*)&xmax);
  ymax = fblank;
  ObitInfoListGetTest(in->info, "YMAX", &type, dim, (gpointer*)&ymax);
  xmin = fblank;
  ObitInfoListGetTest(in->info, "XMIN", &type, dim, (gpointer*)&xmin);
  ymin = fblank;
  ObitInfoListGetTest(in->info, "YMIN", &type, dim, (gpointer*)&ymin);

  /* Check data for range? */
  gotXMAX = xmax != fblank;
  gotYMAX = ymax != fblank;
  gotXMIN = xmin != fblank;
  gotYMIN = ymin != fblank;
  gotRange =  gotXMAX && gotYMAX && gotXMIN && gotYMIN;
  if (!gotRange) { /* Need to check data? */
    for (i=0; i<n; i++) {
      if (xx[i]!=fblank) {
	if  (xmax == fblank) xmax = xx[i];
	else xmax = MAX (xx[i], xmax);
	if  (xmin == fblank) xmin = xx[i];
	else xmin = MIN (xx[i], xmin);
      }
      if (y[i]!=fblank) {
	if  (ymax == fblank) ymax = y[i];
	else ymax = MAX (y[i], ymax);
	if  (ymin == fblank) ymin = y[i];
	else ymin = MIN (y[i], ymin);
      }
    }
  } /* end check range */

  /* Are the range values OK? */
  if (xmax<=xmin) { /* Nope - bug out */
    Obit_log_error(err, OBIT_Error, "%s: XMAX(%g) <= XMIN(%g)", 
		   routine, xmax, xmin);
    if (x==NULL) g_free(xx);
    return;
 }

  if (ymax<=ymin) { /* Nope - bug out */
    Obit_log_error(err, OBIT_Error, "%s: YMAX (%g) <= YMIN(%g)", 
		   routine, ymax, ymin);
    if (x==NULL) g_free(xx);
    return;
 }

  /* If autosetting the range, expand a bit */
  if (!gotXMAX) xmax += 0.05 * (xmax-xmin);
  if (!gotXMIN) xmin -= 0.05 * (xmax-xmin);
  if (!gotYMAX) ymax += 0.05 * (ymax-ymin);
  if (!gotYMIN) ymin -= 0.05 * (ymax-ymin);

  /* Adjust character size */
  cscale = (ofloat)csize;
  cscale = MAX (1.0, cscale);
  ObitPlotSetCharSize (in, cscale, err);

  /* set line width */
  lwidth = MAX (1, lwidth);
  ObitPlotSetLineWidth (in, lwidth, err);

  /* set plotting area */
  axis = -2;
  if (in->xLog  && !in->yLog) axis = 10;
  if (!in->xLog &&  in->yLog) axis = 20;
  if ( in->xLog &&  in->yLog) axis = 30;
  ObitPlotSetPlot (in, xmin, xmax, ymin, ymax, just, axis, err);

  /* Init - use buffering */
  cpgbbuf();

  /*  Plot label */
  ObitPlotLabel (in, xlabel, ylabel, title, err);

  /* set ticks */
  ObitPlotDrawAxes (in, xopt, xtick, nxsub, yopt, ytick, nysub, err);

  /* Adjust symbol size */
  sscale = (ofloat)ssize;
  sscale = MAX (1.0, sscale);
  ObitPlotSetCharSize (in, sscale, err);

  /* To pgplot data types */
  ply = g_malloc0(n*sizeof(float));
  for (i=0; i<n; i++) ply[i] = (float)y[i];
  plx = g_malloc0(n*sizeof(float));
  for (i=0; i<n; i++) plx[i] = (float)xx[i];
  
  /* Plot */
  cpgpt ((int)n, plx, ply, (int)lsymbol);
  
  /* Connect the dots? */
  if (doConnect) ObitPlotDrawCurve (in, n, xx, y, err);

  /* Error bars */
  if (e!=NULL) {
    for (i=0; i<n; i++) {
      ObitPlotDrawLine (in, xx[i], y[i]-e[i], xx[i], y[i]+e[i], err);
    }
  }

  /* Reset symbol size */
  ObitPlotSetCharSize (in, cscale, err);

  /* Flush the buffer */
  cpgebuf();
  cpgupdt();

  /* deallocate x index array if allocated */
  if (x==NULL)   g_free(xx);
  if (plx!=NULL) g_free(plx);
  if (ply!=NULL) g_free(ply);
#endif /* HAVE_PGPLOT */

  /* Complain if plotting not available */
  if (!OK) Obit_log_error(err, OBIT_Error, "No plotting package available");

} /* end  ObitPlotXYErr */

/**
 * Mark positions on a plot drawn by ObitPlotContour
 * Positions marked with a cross the size of size
 * Plot should be finalized with ObitPlotFinishPlot after all 
 * drawing on the current frame is finished.
 * \param in      Pointer to existing ObitPlot object.
 * \param image   Image plotted
 *                Descriptor assumed valid
 * \param n       number of positions to plot
 * \param ra      RAs (deg) to plot
 * \param dec     Declinations to plot
 * \param size    size of symbol in pixels
 * \param err ObitErr error stack
 *
 * Optional parameters on in->info
 * \li CSIZE  (int)    Scaling factor for characters(default = 1)
 * \li LWIDTH (int)    Line width (default = 1)
 */
void ObitPlotMarkCross (ObitPlot* in, ObitImage *image, olong n,
			odouble *ra, odouble *dec, ofloat size, ObitErr *err)
{
  olong i;
  ofloat xcen, ycen, dx, dy, xpixo, ypixo;
  olong csize, cscale, lwidth;
  gint32 dim[MAXINFOELEMDIM] = {1,1,1,1,1};
  ObitInfoType type;
#ifdef HAVE_PGPLOT  /* Only if pgplot available */
  ofloat xmin, ymin;
#endif /* HAVE_PGPLOT */
  /*gchar *routine="ObitPlotMarkCross";*/

 /* error checks */
  if (err->error) return;
  g_assert (ObitPlotIsA(in));
  g_assert (ObitImageIsA(image));
  
  /* Get inputs from info */
  csize = 1;
  ObitInfoListGetTest(in->info, "CSIZE", &type, dim, &csize);
  lwidth = 1;
  ObitInfoListGetTest(in->info, "LWIDTH", &type, dim, &lwidth);

  /* Adjust character size */
  cscale = (ofloat)csize;
  cscale = MAX (1.0, cscale);
  ObitPlotSetCharSize (in, cscale, err);

  /* set line width */
  lwidth = MAX (1, lwidth);
  ObitPlotSetLineWidth (in, lwidth, err);

  /* Plot points */
  for (i=0; i<n; i++) {
    /* calculate pgplot coordinates - Get image pixel */
    ObitSkyGeomXYpix (ra[i], dec[i], 
		      image->myDesc->crval[0], image->myDesc->crval[1], 
		      image->myDesc->crpix[0], image->myDesc->crpix[1],
		      image->myDesc->cdelt[0], image->myDesc->cdelt[1],
		      image->myDesc->crota[1], &image->myDesc->ctype[0][4],
		      &xpixo, &ypixo);
#ifdef HAVE_PLPLOT  /* Only if plplot available */
    plplotCoord (xpixo, ypixo, &xcen, &ycen, (PLPointer)in);
    dx = size * image->myDesc->cdelt[0] * in->scalex;
    dy = size * image->myDesc->cdelt[1] * in->scaley;
#endif /* HAVE_PLPLOT */
#ifdef HAVE_PGPLOT  /* Only if pgplot available */
    /* Get plot coordinates */
    xmin = image->myDesc->crval[0] + 
      (1.0 - image->myDesc->crpix[0]) * image->myDesc->cdelt[0];
    ymin = image->myDesc->crval[1] + 
      (1.0 - image->myDesc->crpix[1]) * image->myDesc->cdelt[1];
    xmin -= image->myDesc->cdelt[0];  /* pgplot indexing error */
    ymin -= image->myDesc->cdelt[1];  /* pgplot indexing error */
    xcen = (xmin + xpixo * image->myDesc->cdelt[0]) * 3600.0 / (15.0);
    ycen = (ymin + ypixo * image->myDesc->cdelt[1]) * 3600.0;

    /* cross size */
    dx = size * image->myDesc->cdelt[0] * 3600.0 / 15.0;
    dy = size * image->myDesc->cdelt[1] * 3600.0;
    
#endif /* HAVE_PGPLOT */
    /* Plot cross as error bars */
    ObitPlotDrawLine (in, xcen-dx, ycen,  xcen+dx, ycen, err);
    ObitPlotDrawLine (in, xcen,   ycen-dy,  xcen, ycen+dy, err);
  } /* end loop over points */
} /* end  ObitPlotMarkCross */

/**
 * Finalize plot to be called after all drawing is complete.
 * \param in  Pointer to Plot object.
 * \param err ObitErr error stack, return if existing error
 */
void ObitPlotFinishPlot (ObitPlot* in, ObitErr *err)
{

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

/****************** plplot implementation *************************/
#ifdef HAVE_PLPLOT  /* Only if plplot available */
  plend();
#endif /* HAVE_PLPLOT */

/****************** pgplot implementation *************************/
#ifdef HAVE_PGPLOT  /* Only if pgplot available */
  cpgend();
#endif /* HAVE_PGPLOT */

} /* end  ObitPlotFinishPlot */

/** 
 * Define plotting area
 * \param in    Pointer to Plot object.
 * \param xmin  the world x-coordinate at the bottom left corner of the viewport.
 * \param xmax  the world x-coordinate at the top right corner of the viewport 
 *                  (note XMAX may be less than XMIN).
 * \param ymin  the world y-coordinate at the bottom left corner
 *                  of the viewport.
 * \param ymax  the world y-coordinate at the top right corner
 *                  of the viewport (note YMAX may be less than YMIN)
 * \param just  if JUST=1, the scales of the x and y axes (in
 *                  world coordinates per inch) will be equal,
 *                  otherwise they will be scaled independently.
 * \param axis  controls the plotting of axes, tick marks, etc:
 * \li axis = -2 : draw no box, axes or labels;
 * \li axis = -1 : draw box only;
 * \li axis =  0 : draw box and label it with coordinates;
 * \li axis =  1 : same as axis=0, but also draw the
 *                coordinate axes (X=0, Y=0);
 * \li axis =  2 : same as axis=1, but also draw grid lines
 *                at major increments of the coordinates;
 * \li axis = 10 : draw box and label X-axis logarithmically;
 * \li axis = 20 : draw box and label Y-axis logarithmically;
 * \li axis = 30 : draw box and label both axes logarithmically.
 * \param err   ObitErr error stack
 */
void  ObitPlotSetPlot (ObitPlot* in, 
		       ofloat xmin, ofloat xmax, ofloat ymin, ofloat ymax, 
		       olong just, olong axis, ObitErr *err)
{
  gboolean OK=FALSE;  /* Have an implementation? */

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

/****************** plplot implementation *************************/
#ifdef HAVE_PLPLOT  /* Only if plplot available */
  OK = TRUE;  /* Have a plotting package */
  /* Call plplot routine */
  plenv ((PLFLT)xmin, (PLFLT)xmax, (PLFLT)ymin, (PLFLT)ymax, (PLINT)just, (PLINT)axis);
#endif /* HAVE_PLPLOT */

/****************** pgplot implementation *************************/
#ifdef HAVE_PGPLOT  /* Only if pgplot available */
  OK = TRUE;  /* Have a plotting package */
  /* Call pgplot routine */
  cpgenv ((float)xmin, (float)xmax, (float)ymin, (float)ymax, (int)just, (int)axis);
#endif /* HAVE_PGPLOT */

  /* Complain if plotting not available */
  if (!OK) Obit_log_error(err, OBIT_Error, "No plotting package available");
} /* end ObitPlotSetPlot */

/** 
 * Front end to pgplot routine cpglab.
 * write labels for x-axis, y-axis, and top of plot
 * Write labels outside the viewport. This routine is a simple
 * interface to PGMTXT, which should be used if PGLAB is inadequate.
 * \param in      Pointer to Plot object.
 * \param xlabel  a label for the x-axis (centered below the
 *                viewport).
 * \param ylabel  a label for the y-axis (centered to the left
 *                  of the viewport, drawn vertically)
 * \param title   a label for the entire plot (centered above the viewport)
 * \param err     ObitErr error stack
 */
void  ObitPlotLabel (ObitPlot* in, gchar *xlabel, gchar *ylabel, gchar *title,
		      ObitErr *err)
{
  gboolean OK=FALSE;  /* Have an implementation? */

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (xlabel!=NULL);
  g_assert (ylabel!=NULL);
  g_assert (title!=NULL);

/****************** plplot implementation *************************/
#ifdef HAVE_PLPLOT  /* Only if plplot available */
  OK = TRUE;  /* Have a plotting package */
  /* Call plplot routine */
  pllab ((char*)xlabel, (char*)ylabel, (char*)title);
#endif /* HAVE_PLPLOT */

/****************** pgplot implementation *************************/
#ifdef HAVE_PGPLOT  /* Only if pgplot available */
  OK = TRUE;  /* Have a plotting package */
  /* Call pgplot routine */
  cpglab ((char*)xlabel, (char*)ylabel, (char*)title);
#endif /* HAVE_PGPLOT */

  /* Complain if plotting not available */
  if (!OK) Obit_log_error(err, OBIT_Error, "No plotting package available");
} /* end ObitPlotLabel */

/** 
 * Draw axes for a plot,  annotate the viewport with frame, axes, numeric labels, etc. 
 * \param in      Pointer to Plot object.
 * \param xopt   string of options for X (horizontal) axis of
 *                  plot. Options are single letters, and may be in
 *                  any order (see below).
 * \param xtick    world coordinate interval between major tick marks
 *                 on X axis. If xtick=0.0, the interval is chosen.
 * \param nxsub   the number of subintervals to divide the major
 *                  coordinate interval into. If xtick=0.0 or nxsub=0,
 *                  the number is chosen.
 * \param yopt   string of options for Y (vertical) axis of plot.
 *               Coding is the same as for xopt.
 * \param ytick  like xtick for the Y axis.
 * \param nysub  like nxsub for the Y axis
 * \param err    ObitErr error stack

 * Axis options:
 * \li A : draw Axis (X axis is horizontal line Y=0, Y axis is vertical line X=0).
 * \li B : draw bottom (X) or left (Y) edge of frame.
 * \li C : draw top (X) or right (Y) edge of frame.
 * \li G : draw Grid of vertical (X) or horizontal (Y) lines
 * \li I : Invert the tick marks; ie draw them outside the viewport instead of inside.
 * \li L : label axis Logarithmically
 * \li N : write Numeric labels in the conventional location below the
 *         viewport (X) or to the left of the viewport (Y).
 * \li M : write numeric labels in the unconventional location above the
 *         viewport (X) or to the right of the viewport (Y).
 * \li P : extend ("Project") major tick marks outside the box (ignored if
 *         option I is specified)
 * \li T : draw major Tick marks at the major coordinate interval.
 * \li S : draw minor tick marks (Subticks).
 */
void  ObitPlotDrawAxes (ObitPlot* in, 
			gchar *xopt, ofloat xtick, olong nxsub, 
			gchar *yopt, ofloat ytick, olong nysub, 
			ObitErr *err)
{
  gboolean OK=FALSE;  /* Have an implementation? */

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (xopt!=NULL);
  g_assert (yopt!=NULL);

/****************** plplot implementation *************************/
#ifdef HAVE_PLPLOT  /* Only if plplot available */
  OK = TRUE;  /* Have a plotting package */
  /* Call plplot routine */
  plbox ((char*)xopt, (PLFLT)xtick, (PLINT)nxsub, (char*)yopt, 
	 (PLFLT)ytick, (PLINT)nysub);
#endif /* HAVE_PLPLOT */

/****************** pgplot implementation *************************/
#ifdef HAVE_PGPLOT  /* Only if pgplot available */
  OK = TRUE;  /* Have a plotting package */
  /* Call pgplot routine */
  cpgbox ((char*)xopt, (float)xtick, (int)nxsub, (char*)yopt, (float)ytick, (int)nysub);
#endif /* HAVE_PGPLOT */

  /* Complain if plotting not available */
  if (!OK) Obit_log_error(err, OBIT_Error, "No plotting package available");
 } /* end ObitPlotDrawAxes */

/** 
 * Set scaling for characters
 * The size affects all text and graph markers drawn later in the program. 
 * \param in      Pointer to Plot object.
 * \param cscale  new character size (dimensionless multiple of the default size).
 * \param err     ObitErr error stack
 */
void  ObitPlotSetCharSize (ObitPlot* in, ofloat cscale, ObitErr *err)
{
  gboolean OK=FALSE;  /* Have an implementation? */

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

#ifdef HAVE_PLPLOT  /* Only if plplot available */
  OK = TRUE;  /* Have a plotting package */
  /* Call plplot routine */
  plschr (0, (PLFLT)cscale);
#endif /* HAVE_PLPLOT */

/****************** pgplot implementation *************************/
#ifdef HAVE_PGPLOT  /* Only if pgplot available */
  OK = TRUE;  /* Have a plotting package */
  /* Call pgplot routine */
  cpgsch ((float)cscale);
#endif /* HAVE_PGPLOT */

  /* Complain if plotting not available */
  if (!OK) Obit_log_error(err, OBIT_Error, "No plotting package available");
 
} /* end ObitPlotSetCharSize */

/** 
 * Set line width
 * \param in      Pointer to Plot object.
 * \param lwidth  Width of line, multiple of default
 * \param err     ObitErr error stack
 */
void  ObitPlotSetLineWidth (ObitPlot* in, olong lwidth, ObitErr *err)
{

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

/****************** plplot implementation *************************/
#ifdef HAVE_PLPLOT  /* Only if plplot available */
  /* Call plplot routine */
  plwid ((PLINT)lwidth);
#endif /* HAVE_PLPLOT */

/****************** pgplot implementation *************************/
#ifdef HAVE_PGPLOT  /* Only if pgplot available */
  /* Call pgplot routine */
  cpgslw ((int)lwidth);
#endif /* HAVE_PGPLOT */

} /* end ObitPlotSetLineWidth */

/** 
 * Set line style, stays in effect until another call
 * \param in      Pointer to Plot object. Actually applies to all
 * \param lstyle  Style of line, 
 *                1 = continious, 2 = dashed, 3=dot dash
 *                4 = dotted, 5 = dash dot dot dot
 * \param err     ObitErr error stack
 */
void  ObitPlotSetLineStyle (ObitPlot* in, olong lstyle, ObitErr *err)
{
#ifdef HAVE_PLPLOT  /* Only if plplot available */
  PLINT n1=0;
  PLINT mark1[1]={0};
  PLINT space1[1]={0};
  PLINT n2=1;
  PLINT mark2[1]={1000};
  PLINT space2[1]={1000};
  PLINT n3=2;
  PLINT mark3[2]={200,1000};
  PLINT space3[2]={500,500};
  PLINT n4=1;
  PLINT mark4[1]={200};
  PLINT space4[1]={500};
  PLINT n5=4;
  PLINT mark5[4]={1000,200,200,200};
  PLINT space5[4]={500,500,500,500};
#endif /* HAVE_PLPLOT */
  gchar *routine = "ObitPlotSetLineStyle";

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  if ((lstyle<1) || (lstyle>5)) {
    Obit_log_error(err, OBIT_Error, "%s ERROR invalud line style %d", 
		   routine, lstyle);
    return;
  }
/****************** plplot implementation *************************/
#ifdef HAVE_PLPLOT  /* Only if plplot available */
  /* Call plplot routine by style */
  if (lstyle==1) { /* continuous */
    plstyl (n1, mark1, space1);
  } else if (lstyle==2) {   /* dash */
    plstyl (n2, mark2, space2);
  } else if (lstyle==3) {   /* dot dash */
    plstyl (n3, mark3, space3);
  } else if (lstyle==4) {   /* dot */
    plstyl (n4, mark4, space4);
  } else if (lstyle==5) {   /* dash dot dot dot */
    plstyl (n5, mark5, space5);
  }
#endif /* HAVE_PLPLOT */

/****************** pgplot implementation *************************/
#ifdef HAVE_PGPLOT  /* Only if pgplot available */
  /* Call pgplot routine */
  cpgsls ((int)lstyle);
#endif /* HAVE_PGPLOT */

} /* end ObitPlotSetLineStyle */

/** 
 * Set foreground color
 * \param in      Pointer to Plot object.
 * \param color   color index (1-15)
 *                0=black, 1=red(default), 2=yellow, 3=green, 
 *                4=aquamarine, 5=pink, 6=wheat, 7=gray, 8=brown,
 *                9=blue, 10=BlueViolet, 11=cyan, 12=turquoise
 *                13=magenta, 14=salmon, 15=white
 * \param err     ObitErr error stack
 */
void  ObitPlotSetColor (ObitPlot* in, olong color, ObitErr *err)
{

  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));

/****************** plplot implementation *************************/
#ifdef HAVE_PLPLOT  /* Only if plplot available */
  /* Call plplot routine */
  plcol0 ((PLINT)color);
#endif /* HAVE_PLPLOT */

/****************** pgplot implementation *************************/
#ifdef HAVE_PGPLOT  /* Only if pgplot available */
  /* Call pgplot routine */
  cpgsci ((int)color);
#endif /* HAVE_PGPLOT */

} /* end ObitPlotSetColor */

/** 
 * Set or advance sub page
 * \param in      Pointer to Plot object.
 * \param sub     if <=0 advance page, if >0 set current subpage to sub
 * \param err     ObitErr error stack
 */
void  ObitPlotSetPage (ObitPlot* in, olong sub, ObitErr *err)
{

/****************** plplot implementation *************************/
  /* error checks */
  g_assert (ObitErrIsA(err));
#ifdef HAVE_PLPLOT  /* Only if plplot available */
  /* Call plplot routine */
  pladv ((PLINT)sub);
#endif /* HAVE_PLPLOT */

/****************** pgplot implementation *************************/
#ifdef HAVE_PGPLOT  /* Only if pgplot available */
  olong ix, iy;
  /* error checks */
  g_assert (ObitErrIsA(err));
  /* Call pgplot routine */
  if (sub<=0) cpgpage ();
  else {
    iy = sub/in->nx;
    ix = sub - iy*in->nx;
    cpgpanl((int)ix, (int)iy);
  }
#endif /* HAVE_PGPLOT */

} /* end ObitPlotSetPage */

/** 
 * Write text at a position specified.
 * This routine is useful for annotating graphs. 
 * The text is written using the current values of attributes color-index, 
 * line-width, character-height, and character-font.
 * \param in      Pointer to Plot object.
 * \param x       Plot x in world coordinates
 * \param y       Plot y in world coordinates
 * \param dx      x component of inclination
 * \param dy      y component of inclination
 * \param just    Controls justification of the string parallel to
 *                the specified edge of the viewport. If
 *                FJUST = 0.0, the left-hand end of the string will
 *                be placed at (x,y); if JUST = 0.5, the center of
 *                the string will be placed at (x,y); if JUST = 1.0,
 *                the right-hand end of the string will be placed at
 *                at (x,y). Other values between 0 and 1 give inter-
 *                mediate placing, but they are not very useful.
 * \param text    The text string to be plotted. Trailing spaces are
 *                ignored when justifying the string, but leading
 *                spaces are significant.
 * \param err     ObitErr error stack
 */
void  ObitPlotText (ObitPlot* in, ofloat x, ofloat y, 
		    ofloat dx, ofloat dy,
		    ofloat just, gchar *text,
		    ObitErr *err)
{
  gboolean OK=FALSE;  /* Have an implementation? */
  
  /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (text!=NULL);

/****************** plplot implementation *************************/
#ifdef HAVE_PLPLOT  /* Only if plplot available */
  OK = TRUE;  /* Have a plotting package */
  /* Call plplot routine */
  plptex ((PLFLT)x, (PLFLT)y, (PLFLT)dx, (PLFLT)dy, (PLFLT)just, (char*)text);
#endif /* HAVE_PLPLOT */

/****************** pgplot implementation *************************/
#ifdef HAVE_PGPLOT  /* Only if pgplot available */
  OK = TRUE;  /* Have a plotting package */
  /* Call pgplot routine */
  cpgptxt ((float)x, (float)y, (float)57.296*atan2(dy,dx+1.0e-20), (float)just, (char*)text);
#endif /* HAVE_PGPLOT */

  /* Complain if plotting not available */
  if (!OK) Obit_log_error(err, OBIT_Error, "No plotting package available");
} /* end ObitPlotText */

/** 
 * Write text at a position specified relative to the viewport (outside
 * or inside).  This routine is useful for annotating graphs. 
 * The text is written using the current values of attributes color-index, 
 * line-width, character-height, and character-font.
 * \param in      Pointer to Plot object.
 * \param side    Must include one of the characters 'B', 'L', 'T',
 *                or 'R' signifying the Bottom, Left, Top, or Right
 *                margin of the viewport. If it includes 'LV' or
 *                'RV', the string is written perpendicular to the
 *                frame rather than parallel to it.
 * \param disp    The displacement of the character string from the
 *                specified edge of the viewport, measured outwards
 *                from the viewport in units of the character
 *                height. Use a negative value to write inside the
 *                viewport, a positive value to write outside.
 * \param coord   The location of the character string along the
 *                specified edge of the viewport, as a fraction of
 *                the length of the edge.
 * \param just    Controls justification of the string parallel to
 *                the specified edge of the viewport. If
 *                just = 0.0, the left-hand end of the string will
 *                be placed at coord; if just = 0.5, the center of
 *                the string will be placed at coord; if just = 1.0,
 *                the right-hand end of the string will be placed at
 *                at coord. Other values between 0 and 1 give inter-
 *                mediate placing, but they are not very useful.
 * \param text    The text string to be plotted. Trailing spaces are
 *                ignored when justifying the string, but leading
 *                spaces are significant.
 * \param err     ObitErr error stack
 */
void  ObitPlotRelText (ObitPlot* in, gchar *side, ofloat disp, 
		       ofloat coord, ofloat fjust, gchar *text,
		       ObitErr *err)
{
   gboolean OK=FALSE;  /* Have an implementation? */
#ifdef HAVE_PLPLOT  /* Only if plplot available */
   char lside[104], *ctemp; 
#endif /* HAVE_PLPLOT */

   /* error checks */
  g_assert (ObitErrIsA(err));
  if (err->error) return;
  g_assert (ObitIsA(in, &myClassInfo));
  g_assert (side!=NULL);
  g_assert (text!=NULL);

/****************** plplot implementation *************************/
#ifdef HAVE_PLPLOT  /* Only if plplot available */
  OK = TRUE;  /* Have a plotting package */
  ctemp = g_ascii_strdown (side,-1);  /* Convert to lower case */
  strncpy (lside, ctemp, 100);
  g_free(ctemp);
  /* Call plplot routine */
  plmtex (lside, (PLFLT)disp, (PLFLT)coord, (PLFLT)fjust, (char*)text);
#endif /* HAVE_PLPLOT */

/****************** pgplot implementation *************************/
#ifdef HAVE_PGPLOT  /* Only if pgplot available */
  OK = TRUE;  /* Have a plotting package */
  /* Call pgplot routine */
  cpgmtxt (side, (float)disp, (float)coord, (float)fjust, (char*)text);
#endif /* HAVE_PGPLOT */

  /* Complain if plotting not available */
  if (!OK) Obit_log_error(err, OBIT_Error, "No plotting package available");
} /* end ObitPlotRelText */

/** 
 * Primitive routine draw a line from (x1,y1) to (x2,y2)
 * \param in      Pointer to Plot object.
 * \param x1      world x-coordinate of the new pen position.
 * \param y1      world y-coordinate of the new pen position.
 * \param x2      world x-coordinate of the new pen position.
 * \param y2      world y-coordinate of the new pen position.
 * \param err     ObitErr error stack
 */
void  ObitPlotDrawLine (ObitPlot* in, ofloat x1, ofloat y1, 
			ofloat x2, ofloat y2, ObitErr *err)
{
  gboolean OK=FALSE;  /* Have an implementation? */
#ifdef HAVE_PLPLOT  /* Only if plplot available */
  PLFLT lx1, ly1, lx2, ly2;
#endif /* HAVE_PLPLOT */

  /* error checks */
  if (err->error) return;

/****************** plplot implementation *************************/
#ifdef HAVE_PLPLOT  /* Only if plplot available */
  OK = TRUE;  /* Have a plotting package */
  /* Plotting logs? */
  if (in->xLog) {lx1 = (PLFLT)log10(x1);lx2 = (PLFLT)log10(x2);}
  else {lx1 = (PLFLT)x1; lx2 = (PLFLT)x2;}
  if (in->yLog) {ly1 = (PLFLT)log10(y1);ly2 = (PLFLT)log10(y2);}
  else {ly1 = (PLFLT)y1; ly2 = (PLFLT)y2;}
  /* Call plplot routine */
  pljoin (lx1, ly1, lx2, ly2);
#endif /* HAVE_PLPLOT */

/****************** pgplot implementation *************************/
#ifdef HAVE_PGPLOT  /* Only if pgplot available */
  OK = TRUE;  /* Have a plotting package */
  /* Call pgplot routine */
  cpgmove ((float)x1, (float)y1);
  cpgdraw ((float)x2, (float)y2);
#endif /* HAVE_PGPLOT */

  /* Complain if plotting not available */
  if (!OK) Obit_log_error(err, OBIT_Error, "No plotting package available");
} /* end ObitPlotDrawLine */

/** 
 * Primitive routine draw a symbol at (x,y)
 * \param in      Pointer to Plot object.
 * \param x       world x-coordinate of the center of the symbol
 * \param y       world y-coordinate of the center of the symbol
 * \param symbol  Symbol index to use for plotting 
 *                values in the range [1,12] are usable 
 * \li 1 = dot
 * \li 2 = plus
 * \li 3 = *
 * \li 4 = open circle
 * \li 5 = x
 * \li 6 = open square
 * \li 7 = open triangle
 * \li 8 = open star
 * \li 9 = filled triangle
 * \li 10 = filled square
 * \li 11 = filled circle
 * \li 12 = filled star
 * \param err     ObitErr error stack
 */
void  ObitPlotDrawSymbol (ObitPlot* in, ofloat x, ofloat y, 
			  olong symbol, ObitErr *err)
{
  gboolean OK=FALSE;  /* Have an implementation? */

/****************** plplot implementation *************************/
#ifdef HAVE_PLPLOT  /* Only if plplot available */
  olong lsymbol;
  olong syms[12] = {210,225,228,840,227,841,842,844,852,851,850,856};
  PLFLT lx[1], ly[1];

  /* error checks */
  if (err->error) return;

  OK = TRUE;  /* Have a plotting package */
  lsymbol = MAX(1, MIN(12, abs(symbol)));
  lsymbol = syms[lsymbol-1];  /* Find plplot value */

  /* Plotting logs? */
  if (in->xLog) lx[0] = (PLFLT)log10(x);
  else lx[0] = (PLFLT)x;
  if (in->yLog) ly[0] = (PLFLT)log10(y);
  else ly[0] = (PLFLT)y;

  /* Call plplot routine */
  plsym ((PLINT)1, lx, ly, (PLINT)lsymbol);
#endif /* HAVE_PLPLOT */

/****************** pgplot implementation *************************/
#ifdef HAVE_PGPLOT  /* Only if pgplot available */
  olong lsymbol;
  olong syms[12] = {1,2,3,4,5,6,7,12,13,16,17,18};

  /* error checks */
  if (err->error) return;

  OK = TRUE;  /* Have a plotting package */
  lsymbol = MAX(1, MIN(12, abs(symbol)));
  lsymbol = syms[lsymbol-1];  /* Find plplot value */

  /* Call pgplot routine */
  cpgpt1 ((float)x, (float)y, (int)lsymbol);
#endif /* HAVE_PGPLOT */

  /* Complain if plotting not available */
  if (!OK) Obit_log_error(err, OBIT_Error, "No plotting package available");
} /* end ObitPlotDrawSymbol */

/** 
 * Primitive routine draw a curve consisting of lines between a sequence of
 * points.
 * \param in      Pointer to Plot object.
 * \param n       Number of points
 * \param x       Array of world x-coordinates of points
 * \param y       Array of world y-coordinates of points
 * \param err     ObitErr error stack
 */
void  ObitPlotDrawCurve (ObitPlot* in, olong n, ofloat *x, ofloat *y, 
			 ObitErr *err)
{
  gboolean OK=FALSE;  /* Have an implementation? */
  olong i;
#ifdef HAVE_PLPLOT  /* Only if plplot available */
  ofloat fblank = ObitMagicF();
  ofloat *logx=NULL, *logy=NULL, *xp, *yp, *ply=NULL, *plx=NULL;
#endif /* HAVE_PLPLOT */
#ifdef HAVE_PGPLOT  /* Only if plplot available */
  float *ply=NULL, *plx=NULL;
#endif /* HAVE_PGPLOT */

  /* error checks */
  if (err->error) return;

/****************** plplot implementation *************************/
#ifdef HAVE_PLPLOT  /* Only if plplot available */
  OK = TRUE;  /* Have a plotting package */
  /* Log? */
  if (in->xLog) {
    logx = g_malloc0(n*sizeof(PLFLT));
    for (i=0; i<n; i++) logx[i] = (PLFLT)log10(x[i]);
    xp = logx;
  } else {
    plx = g_malloc0(n*sizeof(PLFLT));
    for (i=0; i<n; i++) plx[i] = (PLFLT)x[i];
    xp = plx;
  }

  if (in->yLog) {
    logy = g_malloc0(n*sizeof(PLFLT));
    for (i=0; i<n; i++) {
      if (x[i]!=fblank) {
	logy[i] = (PLFLT)log10(y[i]);
      } else logy[i] = (PLFLT)1.0e30;  /* Don't plot if blanked */
    }
    yp = logy;
  } else {
    ply = g_malloc0(n*sizeof(PLFLT));
    for (i=0; i<n; i++) ply[i] = (PLFLT)y[i];
    yp = ply;
  }

  /* Call plplot routine */
  plline ((PLINT)n, xp, yp);

  /* Flush the buffer */
  plflush();

  /* deallocate arrays if allocated */
  if (logx!=NULL)  g_free(logx);
  if (logy!=NULL)  g_free(logy);
  if (plx!=NULL)   g_free(plx);
  if (ply!=NULL)   g_free(ply);
#endif /* HAVE_PLPLOT */

/****************** pgplot implementation *************************/
#ifdef HAVE_PGPLOT  /* Only if pgplot available */
  OK = TRUE;  /* Have a plotting package */
  
  /* To pgplot data types */
  ply = g_malloc0(n*sizeof(float));
  for (i=0; i<n; i++) ply[i] = (float)y[i];
  plx = g_malloc0(n*sizeof(float));
  for (i=0; i<n; i++) plx[i] = (float)x[i];
  
  /* Call pgplot routine */
  for (i=1; i<n; i++) {
    cpgmove (plx[i-1], ply[i-1]);
    cpgdraw (plx[i], ply[i]);
  }

  if (plx!=NULL) g_free(plx);
  if (ply!=NULL) g_free(ply);

  /* Flush the buffer */
  cpgebuf();
  cpgupdt();

#endif /* HAVE_PGPLOT */

  /* Complain if plotting not available */
  if (!OK) Obit_log_error(err, OBIT_Error, "No plotting package available");
} /* end ObitPlotDrawCurve */

/**
 * Initialize global ClassInfo Structure.
 */
void ObitPlotClassInit (void)
{
  if (myClassInfo.initialized) return;  /* only once */
  
  /* Set name and parent for this class */
  myClassInfo.ClassName   = g_strdup(myClassName);
  myClassInfo.ParentClass = ObitParentGetClass();

  /* Set function pointers */
  ObitPlotClassInfoDefFn ((gpointer)&myClassInfo);
 
  myClassInfo.initialized = TRUE; /* Now initialized */
 
} /* end ObitPlotClassInit */

/**
 * Initialize global ClassInfo Function pointers.
 */
static void ObitPlotClassInfoDefFn (gpointer inClass)
{
  ObitPlotClassInfo *theClass = (ObitPlotClassInfo*)inClass;
  ObitClassInfo *ParentClass = (ObitClassInfo*)myClassInfo.ParentClass;

  if (theClass->initialized) return;  /* only once */

  /* Check type of inClass */
  g_assert (ObitInfoIsA(inClass, (ObitClassInfo*)&myClassInfo));

  /* Initialize (recursively) parent class first */
  if ((ParentClass!=NULL) && 
      (ParentClass->ObitClassInfoDefFn!=NULL))
    ParentClass->ObitClassInfoDefFn(theClass);

  /* function pointers defined or overloaded this class */
  theClass->ObitClassInit = (ObitClassInitFP)ObitPlotClassInit;
  theClass->ObitClassInfoDefFn = (ObitClassInfoDefFnFP)ObitPlotClassInfoDefFn;
  theClass->ObitGetClass  = (ObitGetClassFP)ObitPlotGetClass;
  theClass->newObit       = (newObitFP)newObitPlot;
  theClass->ObitCopy      = (ObitCopyFP)ObitPlotCopy;
  theClass->ObitClone     = NULL;
  theClass->ObitClear     = (ObitClearFP)ObitPlotClear;
  theClass->ObitInit      = (ObitInitFP)ObitPlotInit;

} /* end ObitPlotClassDefFn */

/*---------------Private functions--------------------------*/
/**
 * Creates empty member objects, initialize reference count.
 * Does (recursive) initialization of base class members before 
 * this class.
 * \param inn Pointer to the object to initialize.
 */
void ObitPlotInit  (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitPlot *in = inn;
  gboolean OK=FALSE;  /* Have an implementation? */

  /* error checks */
  g_assert (in != NULL);
  
#ifdef HAVE_PLPLOT  /* Only if plplot available */
  OK = TRUE;  /* Have a plotting package */
#endif /* HAVE_PLPLOT */

#ifdef HAVE_PGPLOT  /* Only if pgplot available */
  OK = TRUE;  /* Have a plotting package */
#endif /* not HAVE_PGPLOT */

  if (!OK) g_error ("NO Plotting package available");

  /* recursively initialize parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  if ((ParentClass!=NULL) && ( ParentClass->ObitInit!=NULL)) 
    ParentClass->ObitInit (inn);

  /* set members in this class */
  in->thread  = newObitThread();
  in->info    = newObitInfoList();
  in->myImage = NULL;
  in->myErr   = NULL;

} /* end ObitPlotInit */

/**
 * Deallocates member objects.
 * Does (recursive) deallocation of parent class members.
 * \param  inn Pointer to the object to deallocate.
 */
void ObitPlotClear (gpointer inn)
{
  ObitClassInfo *ParentClass;
  ObitPlot *in = inn;

  /* error checks */
  g_assert (ObitIsA(in, &myClassInfo));

  /* free this class members */
  in->thread    = ObitThreadUnref(in->thread);
  if (in->info) ObitInfoListUnref (in->info); in->info = NULL;
  
  in->myImage = ObitImageUnref(in->myImage);
  in->myErr   = ObitErrUnref(in->myErr);

 /* unlink parent class members */
  ParentClass = (ObitClassInfo*)(myClassInfo.ParentClass);
  /* delete parent class members */
  if ((ParentClass!=NULL) && ( ParentClass->ObitClear!=NULL)) 
    ParentClass->ObitClear (inn);

} /* end ObitPlotClear */

#ifdef HAVE_PLPLOT  /* Only if plplot available */
/**
 * Parse output specifier into device and file
 * \param  output output specifier in form filename/device
 *         if no '/' is found the string is assumed to be the device
 *         with no output file. Max. 100 characters.
 * \param  dev    [out] device name, should be g_freed when done
 * \param  file   [out] output file name, should be g_freed when done
 */
static void parseOutput (gchar *output, gchar **dev, gchar **file)
{
  olong i, off, slen;

  *dev  = NULL;
  *file = NULL;
  if (output==NULL) return;

  /* Find any '/' */
  off = -1;
  slen = MIN (100, strlen(output));
  for (i=0; i<slen; i++) {
    if (output[i]=='/') {off=i; break;}
  }
  
  /* find it? */
  if (off<0) {
    *dev = g_strndup (output, slen);
    return;
  }

  /* Device */
  *dev = g_strndup (&output[off+1], slen);

  /* File? */
  if (off>1) *file = g_strndup (output, off);
} /* end parseOutput */

/**
 * Function for plplot to convert pixel numbers to world coordinates
 * \param  px   [in] X pixel (0-rel)
 * \param  py   [in] Y pixel (0-rel)
 * \param  wx   [out] X world coordinate value
 * \param  wy   [out] Y world coordinate value
 * \param  ID   [in] Image descriptor
 */
static void plplotCoord (PLFLT px, PLFLT py, PLFLT* wx, PLFLT *wy, 
			  PLPointer OP)
{
  ofloat pixel[2];
  odouble pos[2], xo, yo;
  ObitImageDesc *id = ((ObitPlot*)OP)->myImage->myDesc;

  pixel[0] = (ofloat)px;
  pixel[1] = (ofloat)py;
  /* Convert to position */
  ObitImageDescGetPos (id, pixel, pos, ((ObitPlot*)OP)->myErr); 

  /* Offset in this projection */
  if (ObitSkyGeomXYPixLM (pos[0], pos[1], id->crval[0], id->crval[1],
			  id->cdelt[0], id->cdelt[1],id->crota[1],
			  id->ctype[0], &xo, &yo)==0) {
    *wx = (PLFLT)xo*((ObitPlot*)OP)->scalex;
    *wy = (PLFLT)yo*((ObitPlot*)OP)->scaley;
  } else { /* out of bounds */
    
    *wx = (PLFLT)pos[0] - ((ObitPlot*)OP)->myImage->myDesc->crval[0];
    *wy = (PLFLT)pos[1] - ((ObitPlot*)OP)->myImage->myDesc->crval[1];
    *wx *= ((ObitPlot*)OP)->scalex;
    *wy *= ((ObitPlot*)OP)->scaley;
  }
} /* end plplotCoord */

/**
 * Reorder pixel values to plplot order
 * \param  fa   ObitFArray to reorder (2-D)
 * \returns 2D PLFLT array, should be g_freeed when done
 */
static PLFLT** reorderPixels (ObitFArray *fa)
{
  olong i, j, nx=fa->naxis[0], ny=fa->naxis[1], offset;
  ofloat fblank = ObitMagicF();
  PLFLT** out=NULL;

  /* Allocate */
  out = g_malloc(nx*sizeof(PLFLT*));
  for (i=0; i<nx; i++) out[i] = g_malloc(ny*sizeof(PLFLT));

  /* Copy */
  for (i=0; i<nx; i++) {
    for (j=0; j<ny; j++) {
      offset = i + j*nx;
      if (fa->array[offset]!=fblank) {
	out[i][j] = (PLFLT)fa->array[offset];
      } else {
	out[i][j] = (PLFLT)0.0;
      }
    }
  }
  return out;
} /* end reorderPixels */
#endif /* HAVE_PLPLOT */



