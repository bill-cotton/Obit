#include <stdio.h>
#include <stdlib.h>
int main ( int argc, char **argv )
{
  int irecord[256];
  FILE *file;
  char *filename="dump.inp";
  long filepos;

  file = fopen (filename, "rb");
  filepos = 7*1024;
  fseek (file, filepos, SEEK_SET);
  fread ((char*)irecord, 1024, 1, file);

  return 0;
} /* end main */
