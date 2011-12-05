# autoconfig script for bojan Nicoklic's libAir (ALMA WVR) library
AC_DEFUN([AC_PATH_WVR], [
	WVR_CFLAGS=""
	WVR_LDFLAGS=""
	AC_ARG_WITH(wvr,
                    AC_HELP_STRING([--with-wvr=DIR],
                                 [search for WVR in DIR/include and DIR/lib]),
                    [for dir in `echo "$withval" | tr : ' '`; do
    if test -d $dir/include; then
      WVR_CFLAGS="$WVR_CFLAGS -I$dir/include/"
    fi
    if test -d $dir/lib; then
      WVR_LDFLAGS="$WVR_LDFLAGS -L$dir/lib"
    fi
  done[]])

        AC_ARG_WITH(wvr-includes,
                    AC_HELP_STRING([--with-wvr-includes=DIR],
	                           [search for WVR includes in DIR]),
	            [for dir in `echo "$withval" | tr : ' '`; do
    if test -d $dir; then
      WVR_CFLAGS="$WVR_CFLAGS -I$dir"
    fi
  done[]])

echo "WVR CFLAGs $WVR_CFLAGS LDFLAGs $WVR_LDFLAGS"
ac_wvr_saved_CFLAGS="$CFLAGS"
ac_wvr_saved_LDFLAGS="$LDFLAGS"
ac_wvr_saved_LIBS="$LIBS"
CFLAGS="$CFLAGS $WVR_CFLAGS"
LDFLAGS="$LDFLAGS $WVR_LDFLAGS"
if ! test WVR_CFLAGS; then
    WVR_CFLAGS="`--cflags`"
fi
# not there WVR_LIBS="`wvr-config --libs`"
ac_have_wvr=no
ac_have_wvrh=no
  	touch /tmp/dummy1_wvr.h
        AC_CHECK_HEADERS([/tmp/dummy1_wvr.h], [ac_have_wvrh=yes], [ac_have_wvrh=no],
			[#include <almawvr/almaabs_c.h>])
	rm /tmp/dummy1_wvr.h
 	if test $ac_have_wvrh = yes; then
	        AC_SEARCH_LIBS(almaabs_ret, [almawvr], [ac_have_wvr=yes], [ac_have_wvr=no], 
			[-lm])
		if test $ac_have_wvr = yes; then
			echo "WVR found1 $LIBS"
			OTHERLIB="-lalmawvr "
		fi
		if test $ac_have_wvr = no; then
			# try other possibility
			ac_cv_search_wvr_almaabs_ret="  "
        		AC_SEARCH_LIBS(wvr_almaabs_ret, [air_cbind], [ac_have_wvr=yes], [ac_have_wvr=no], 
							[-lm -lair -lairapps])
			if test $ac_have_wvr = yes; then
				OTHERLIB="-lair_cbind -lair -lairapps "
				echo "WVR found2 $LIBS"
			fi
		fi
	fi
echo "WVR Test1 $ac_have_wvrh $ac_have_wvr"
# List of places to try
testdirs="$OBITINSTALL/other $HOME/opt/almawvr"
for dir in $testdirs; do
	if test $ac_have_wvr = yes ; then
		if test $ac_have_wvrh = yes ; then
			break;
		fi
	fi
	if test $ac_have_wvr = no; then
		if  test -f $dir/include/almawvr/almaabs_c.h; then
			WVR_CFLAGS="-I$dir/include"
			CPPFLAGS="$ac_wvr_saved_CPPFLAGS $WVR_CFLAGS"
			WVR_LDFLAGS="-L$dir/lib"
			LDFLAGS="$ac_wvr_saved_LDFLAGS $WVR_LDFLAGS"
  			touch /tmp/dummy3_wvr.h
	        	AC_CHECK_HEADERS(/tmp/dummy3_wvr.h, [ac_have_wvrh=yes], [ac_have_wvrh=no],
				[#include "almawvr/almaabs_c.h"])
			rm /tmp/dummy3_wvr.h
			if test $ac_have_wvrh = yes; then
				# Force check
				ac_cv_search_wvr_almaabs_ret="  "
        			AC_SEARCH_LIBS(wvr_almaabs_ret, [almawvr], [ac_have_wvr=yes], [ac_have_wvr=no], 
					[-lm ])
				if test $ac_have_wvr = yes; then
 	                                echo "WVR found3 $LIBS in $dir"
 					OTHERLIB="-lalmawvr"
				fi
				if test $ac_have_wvr = no; then
					# try other possibility
					ac_cv_search_wvr_almaabs_ret="  "
        				AC_SEARCH_LIBS(wvr_almaabs_ret, [air_cbind], [ac_have_wvr=yes], [ac_have_wvr=no], 
						[-lm -lair -lairapps])
					if test $ac_have_wvr = yes; then
	 	                                echo "WVR found4 $LIBS in $dir"
						OTHERLIB=" -lair_cbind -lair -lairapps "
					fi
				fi	
			fi
			if test $ac_have_wvr = yes ; then
				if test $ac_have_wvrh = yes ; then
					break;
				fi
			fi
		fi
	fi
	echo "WVR Testn $ac_have_wvrh $ac_have_wvr"
done[]
if test $ac_have_wvrh = no; then
	AC_MSG_WARN([cannot find WVR headers])
	ac_have_wvr=no
fi
if test $ac_have_wvr = no; then
	AC_MSG_WARN([cannot find WVR library])
	WVR_CFLAGS=""
	WVR_LDFLAGS=""
	WVR_LIBS=""
fi
if test $ac_have_wvr = yes; then
	AC_DEFINE(HAVE_WVR, 1, [Define to 1 if WVR is available.])
fi
WVR_LIBS="$LIBS $OTHERLIB"
CFLAGS="$ac_wvr_saved_CFLAGS"
LDFLAGS="$ac_wvr_saved_LDFLAGS"
LIBS="$ac_wvr_saved_LIBS"
	AC_SUBST(WVR_CFLAGS)
	AC_SUBST(WVR_LDFLAGS)
	AC_SUBST(WVR_LIBS)
])
