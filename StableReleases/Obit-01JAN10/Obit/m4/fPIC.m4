# Check if compiler option -fPIC allowed and if so add to $CFLAGS
 AC_DEFUN([AC_CHECK_FPIC], [
AC_MSG_CHECKING(checking if -fPIC compiler option allowed... )
use_fPIC=no
ac_fpic_saved_CFLAGS="$CFLAGS"
CFLAGS="$CFLAGS -fPIC"
	#AC_DEFINE([HELLO_WORLD], ["Hello, World\n"])
	AC_COMPILE_IFELSE(
            [AC_LANG_PROGRAM([[const char hw[] = "Hello, World\n";]],
                    [[]])],
	   [ac_fpic_saved_CFLAGS="$CFLAGS"; use_fPIC=yes],	
	   [use_fPIC=no])
CFLAGS="$ac_fpic_saved_CFLAGS"
AC_MSG_RESULT($use_fPIC)
])
