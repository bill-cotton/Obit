# Search for Motif/Lesstif directories
# needs AC_PATH_XTRA run first
AC_DEFUN([AC_PATH_MOTIF], [
	AC_ARG_WITH(motif,
                    AC_HELP_STRING([--with-motif=DIR],
                             [search for MOTIF in DIR/include and DIR/lib]),
                    [for dir in `echo "$withval" | tr : ' '`; do
    if test -d $dir/include; then
      MOTIF_CPPFLAGS="$MOTIF_CPPFLAGS -I$dir/include"
    fi
    if test -d $dir/lib; then
      MOTIF_LDFLAGS="$MOTIF_LDFLAGS -L$dir/lib"
    fi
  done[]])

        AC_ARG_WITH(motif-includes,
                    AC_HELP_STRING([--with-motif-includes=DIR],
	                           [search for MOTIF includes in DIR]),
	            [for dir in `echo "$withval" | tr : ' '`; do
    if test -d $dir; then
      MOTIF_CPPFLAGS="$MOTIF_CPPFLAGS -I$dir"
    fi
  done[]])

ac_motif_saved_CFLAGS="$CPPFLAGS"
ac_motif_saved_LDFLAGS="$LDFLAGS"
ac_motif_saved_LIBS="$LIBS"
CPPFLAGS="$CPPFLAGS $MOTIF_CPPFLAGS $X_CFLAGS"
LDFLAGS="$LDFLAGS $MOTIF_LDFLAGS $X_LIBS"

ac_have_motif=no
ac_have_motifh=no
  	touch /tmp/dummy1_Xm.h
        AC_CHECK_HEADER(/tmp/dummy1_Xm.h, ac_have_motifh=yes[], [ac_have_motifh=no],
		[#include "Xm/Xm.h"])
	rm /tmp/dummy1_Xm.h
 	if test $ac_have_motifh = yes; then
 	       AC_SEARCH_LIBS(XmGetPixmap, Xm, [ac_have_motif=yes], [ac_have_motif=no],
		[-lXt -lX11])
	fi
# List of places to try
testdirs="$HOME/opt/Xm $OBITINSTALL/other"
for dir in $testdirs; do
	if test $ac_have_motif = no; then
		if  test -f $dir/include/Xm/Xm.h; then
			MOTIF_CPPFLAGS="-I$dir/include"
			CPPFLAGS="$ac_Xm_saved_CPPFLAGS $MOTIF_CPPFLAGS"
			MOTIF_LDFLAGS="-L$dir/lib"
			LDFLAGS="$ac_Xm_saved_LDFLAGS $MOTIF_LDFLAGS"
  			touch /tmp/dummy3_Xm.h
	        	AC_CHECK_HEADERS(/tmp/dummy3_Xm.h, [ac_have_motifh=yes], [ac_have_motifh=no],
				[#include "Xm/Xm.h"])
			rm /tmp/dummy3_Xm.h
			if test $ac_have_motifh = yes; then
				# Force check
				ac_cv_search_XmGetPixmap="  "
		 	    	AC_SEARCH_LIBS(XmGetPixmap, Xm, [ac_have_motif=yes], [ac_have_motif=no],
					[-lXt -lX11])
			fi
			if test $ac_have_motif = yes ; then
				if test $ac_have_motifh = yes ; then
					break;
				fi
			fi
		fi
	fi
done[]
if test $ac_have_motifh = no; then
	AC_MSG_ERROR([cannot find Motif(or LessTif) headers])
	ac_have_motif=no
fi
if test $ac_have_motif = no; then
	AC_MSG_ERROR([cannot find Motif(or LessTif) library])
	MOTIF_CPPFLAGS=""
	MOTIF_LDFLAGS=""
	MOTIF_LIBS=""
fi
if test $ac_have_motif = yes; then
	AC_DEFINE(HAVE_MOTIF, 1, [Define to 1 if you have MOTIF.])
fi
MOTIF_LIBS="$LIBS -lXt -lX11"
CPPFLAGS="$ac_motif_saved_CPPFLAGS"
LDFLAGS="$ac_motif_saved_LDFLAGS"
LIBS="$ac_motif_saved_LIBS"
	AC_SUBST(MOTIF_CPPFLAGS)
	AC_SUBST(MOTIF_LDFLAGS)
	AC_SUBST(MOTIF_LIBS)
])
