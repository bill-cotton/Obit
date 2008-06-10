/* $Id$    */  
/*--------------------------------------------------------------------*/
/* Swig module description for AIPS directory utilities               */
/*                                                                    */
/*;  Copyright (C) 2004-2008                                          */
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

%{
#include "ObitAIPS.h"
#include "ObitAIPSDir.h"
%}

// This cleans up the char * string we malloc'd in AIPSDirInfo
%typemap(python,ret) char *AIPSDirInfo {
  if ($source) g_free($source);
}

%inline %{
extern int AIPSDirFindCNO(int disk, int user,  char *Aname, char *Aclass, 
    		          char *Atype, int seq, ObitErr *err)
{
  gchar LAname[13], LAclass[7], LAtype[3];
  olong i, l;

  /* Init AIPS fixed strings - blank fill until end */
  for (i=0; i<12; i++) LAname[i]  = ' ';LAname[i]  = 0;
  for (i=0; i<6; i++)  LAclass[i] = ' ';LAclass[i] = 0;
  for (i=0; i<2; i++)  LAtype[i]  = ' ';LAtype[i]  = 0;

  /* Copy string input into AIPS fixed strings */
  l = MIN (12, strlen(Aname));
  for (i=0; i<l; i++) LAname[i]  = Aname[i];
  l = MIN (6, strlen(Aclass));
  for (i=0; i<l; i++)  LAclass[i] = Aclass[i];
  l = MIN (2, strlen(Atype));
  for (i=0; i<l; i++)  LAtype[i] = Atype[i];

  return ObitAIPSDirFindCNO(disk, user, LAname, LAclass, LAtype, seq, err);
} /* end AIPSDirFindCNO */

extern int AIPSDirHiSeq(int disk, int user,  char *Aname, char *Aclass, 
    		        char *Atype, ObitErr *err)
{
  gchar LAname[13], LAclass[7], LAtype[3];
  olong i, l;

  /* Init AIPS fixed strings - blank fill until end */
  for (i=0; i<12; i++) LAname[i]  = ' ';LAname[i]  = 0;
  for (i=0; i<6; i++)  LAclass[i] = ' ';LAclass[i] = 0;
  for (i=0; i<2; i++)  LAtype[i]  = ' ';LAtype[i]  = 0;

  /* Copy string input into AIPS fixed strings */
  l = MIN (12, strlen(Aname));
  for (i=0; i<l; i++) LAname[i]  = Aname[i];
  l = MIN (6, strlen(Aclass));
  for (i=0; i<l; i++)  LAclass[i] = Aclass[i];
  l = MIN (2, strlen(Atype));
  for (i=0; i<l; i++)  LAtype[i] = Atype[i];

  return ObitAIPSDirHiSeq(disk, user, LAname, LAclass, LAtype, TRUE, err);
} /* end AIPSDirHiSeq */


extern int AIPSDirAlloc(int disk, int user,  char *Aname, char *Aclass, 
		        char *Atype, int seq, ObitErr *err)
{
  gchar LAname[13], LAclass[7], LAtype[3];
  gboolean exist;
  olong i, l;

  /* Init AIPS fixed strings - blank fill until end */
  for (i=0; i<12; i++) LAname[i]  = ' ';LAname[i]  = 0;
  for (i=0; i<6; i++)  LAclass[i] = ' ';LAclass[i] = 0;
  for (i=0; i<2; i++)  LAtype[i]  = ' ';LAtype[i]  = 0;

  /* Copy string input into AIPS fixed strings */
  l = MIN (12, strlen(Aname));
  for (i=0; i<l; i++) LAname[i]  = Aname[i];
  l = MIN (6, strlen(Aclass));
  for (i=0; i<l; i++)  LAclass[i] = Aclass[i];
  l = MIN (2, strlen(Atype));
  for (i=0; i<l; i++)  LAtype[i] = Atype[i];

  return ObitAIPSDirAlloc(disk, user, LAname, LAclass, LAtype, seq, &exist, err);
} /* end AIPSDirAlloc */

extern void AIPSDirRemoveEntry(int disk, int user, int cno, ObitErr *err)
{
  ObitAIPSDirRemoveEntry(disk, user, cno, err);
} /* end  AIPSDirRemoveEntry */

extern int AIPSDirNumber(int disk, int user, ObitErr *err)
{
  return (int)ObitAIPSDirNumber((olong)disk, (olong)user, err);
} /* end AIPSDirNumber */

/* Returns NULL if no entry */
extern char* AIPSDirInfo(int disk, int user, int cno, ObitErr *err)
{
  ObitAIPSDirCatEntry *entry;
  gchar *out=NULL, stat[5];
  olong nout=60;

  entry = ObitAIPSDirGetEntry ((olong)disk, (olong)user, (olong)cno, err);
  if ((err->error) || (entry==NULL)) return out;
  if (entry->user!=user) {
      g_free (entry);
      return out;
  }

  if (entry->status==0) strncpy (stat, "    ", 4);
  else if (entry->status<0) strncpy (stat, "WRIT", 4);
  else if (entry->status>0) strncpy (stat, "READ", 4);

  out = g_malloc(nout);
  g_snprintf (out, nout, "%-12.12s.%-6.6s. %4d %-2.2s          ",
              entry->name, entry->class, entry->seq, entry->type);
  ObitAIPSDirGetAccess (entry, &out[29]);
  g_snprintf (&out[strlen(out)], nout-strlen(out), " %-4.4s", stat);
  g_free (entry);
  return (char*)out;
} /* end AIPSDirInfo */

extern int AIPSDirStatus(int disk, int user, int cno, int code, ObitErr *err)
{
  return ObitAIPSDirStatus (disk, user, cno, code, err);
} /* end AIPSDirStatus */

extern int AIPSSetDirname(int disk, char *dir, ObitErr *err)
{
   return ObitAIPSSetDirname(disk, dir, err);
} /* end AIPSSetDirname */


%}
