/* SWIG interface file for ObitVersion.c */

%{
#include <string.h>
extern char *ObitVersion();
%}

%inline %{
extern char *Version() {
    return ObitVersion();
}
%}
