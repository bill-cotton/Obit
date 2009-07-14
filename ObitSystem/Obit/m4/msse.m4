# Check if compiler option -msse allowed and if so add to $CFLAGS
AC_DEFUN([AC_CHECK_SSE], [
AC_MSG_CHECKING(checking if SSE available and -msse compiler option allowed... )
use_msse=no
ac_msse_saved_CFLAGS="$CFLAGS"
CFLAGS="$CFLAGS -msse"
	AC_COMPILE_IFELSE(
            [AC_LANG_PROGRAM([[#include <xmmintrin.h>\n"]],
                    [[]])],
	   [ac_msse_saved_CFLAGS="$CFLAGS"; use_msse=yes],	
	   [use_msse=no])
CFLAGS="$ac_msse_saved_CFLAGS"
AC_MSG_RESULT($use_msse)
if test $use_msse = yes; then
	AC_DEFINE(HAVE_SSE, 1, [Define to 1 if SSE is available.])
fi
])
