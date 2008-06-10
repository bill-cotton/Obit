#include <stdio.h>
#include <stdlib.h>
#include "ObitAll.h"
#include "ObitSource.h"
#include "ObitPrecess.h"

/* program to test Precession */
int main ( int argc, char **argv )
{
  ObitSystem *mySystem;
  ObitErr *err;
  ObitSource *source;
  ObitUVDesc *desc;
  gchar *AIPSdir[] = {"AIPSdata/"};
  gchar *FITSdir[] = {"FITSdata/"};
  olong user;

  /* Initialize Obit */
  err = newObitErr();
  user = 100;
  mySystem = ObitSystemStartup ("test", 1, user, 1, AIPSdir, 1, FITSdir,  
				(oint)TRUE, (oint)FALSE, err);
  ObitErrLog(err); /* show any error messages on err */

  source = newObitSource("Dim bulb");
  /* J2125+0441: App = 3.213810e+02,4.697888e+00 */
  source->RAMean  = 3.213719e+02;
  source->DecMean = 4.693202e+00;

  desc = newObitUVDesc("desc");
  desc->epoch   = 2000.0;
  desc->equinox = 2000.0;
  g_snprintf (desc->obsdat, UVLEN_VALUE-1, "2000-09-19");

  /* precess */
  ObitPrecessUVJPrecessApp(desc, source);

  fprintf (stderr, "Mean %f %f  App %f %f \n",
	   source->RAMean, source->DecMean, source->RAApp, source->DecApp);

  /* zap 'em */
  source = ObitUnref(source);
  desc   = ObitUnref(desc);

  /* show any errors */
  ObitErrLog(err);


  /* Shutdown Obit */
  mySystem = ObitSystemShutdown (mySystem);
  
  return 0;
} /* end of main */

