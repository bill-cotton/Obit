AC_DEFUN([AC_PATH_XMLRPC], [
	AC_ARG_WITH(xmlrpc,
                    AC_HELP_STRING([--with-xmlrpc=DIR],
                                 [search for XMLRPC in DIR/include and DIR/lib]),
                    [for dir in `echo "$withval" | tr : ' '`; do
    if test -d $dir/include; then
      XMLRPC_CPPFLAGS="$XMLRPC_CPPFLAGS -I$dir/include"
    fi
    if test -d $dir/lib; then
      XMLRPC_LDFLAGS="$XMLRPC_LDFLAGS -L$dir/lib"
    fi
  done[]])

        AC_ARG_WITH(xmlrpc-includes,
                    AC_HELP_STRING([--with-xmlrpc-includes=DIR],
	                           [search for XMLRPC includes in DIR]),
	            [for dir in `echo "$withval" | tr : ' '`; do
    if test -d $dir; then
      XMLRPC_CPPFLAGS="$XMLRPC_CPPFLAGS -I$dir"
    fi
  done[]])

XMLRPC="$dir"

	AC_ARG_WITH(www,
                    AC_HELP_STRING([--with-www=DIR],
                                 [search for WWW in DIR/include and DIR/lib]),
                    [for dir in `echo "$withval" | tr : ' '`; do
    if test -d $dir/include; then
      WWW_CPPFLAGS="$WWW_CPPFLAGS -I$dir/include"
    fi
    if test -d $dir/lib; then
      WWW_LDFLAGS="$WWW_LDFLAGS -L$dir/lib"
      WWW_LIBDIR="$dir/lib"
    fi
  done[]])

        AC_ARG_WITH(www-includes,
                    AC_HELP_STRING([--with-www-includes=DIR],
	                           [search for WWW includes in DIR]),
	            [for dir in `echo "$withval" | tr : ' '`; do
    if test -d $dir; then
    	WWW_CPPFLAGS="$WWW_CPPFLAGS -I$dir"
    fi
  done[]])


ac_xmlrpc_saved_CPPFLAGS="$CPPFLAGS"
ac_xmlrpc_saved_LDFLAGS="$LDFLAGS"
ac_xmlrpc_saved_LIBS="$LIBS"
CPPFLAGS="$CPPFLAGS $XMLRPC_CPPFLAGS"
if ! test XMLRPC_CFLAGS; then
  XMLRPC_CFLAGS="`xmlrpc-c-config --cflags`"
fi
XMLRPC_LDFLAGS="$LDFLAGS "
# Libaries needed 
XMLRPC_LIBS="`xmlrpc-c-config --libs` `xmlrpc-c-config client --libs`  "
# Fall back if this fails
if test "$XMLRPC_LIBS""set" == "set"; then
  XMLRPC_LIBS=" -lxmlrpc -lxmlrpc_util -lxmlrpc_xmlparse -lxmlrpc_xmltok "
fi
#echo " XMLRPC_LIBS .$XMLRPC_LIBS."
XMLRPC_SERVER_LIBS="`xmlrpc-c-config abyss-server --libs`"
# Fall back if this fails
if test "$XMLRPC_SERVER_LIBS""set" == "set"; then
  XMLRPC_SERVER_LIBS="-lxmlrpc_server_abyss -lxmlrpc_server -lxmlrpc_abyss  -lpthread -lxmlrpc -lxmlrpc_util -lxmlrpc_xmlparse -lxmlrpc_xmltok"
fi
#echo " XMLRPC_SERVER_LIBS .$XMLRPC_SERVER_LIBS."

XMLRPC_CLIENT_LIBS="`xmlrpc-c-config client --libs`"
# Fall back if this fails
if test "$XMLRPC_CLIENT_LIBS""set" == "set"; then
  XMLRPC_SERVER_LIBS="-lpthread -lxmlrpc -lxmlrpc_util -lxmlrpc_xmlparse -lxmlrpc_xmltok"
fi
#echo " XMLRPC_CLIENT_LIBS .$XMLRPC_CLIENT_LIBS."

ac_have_xmlrpc=no
ac_have_xmlrpch=no
# In older versions this was xmlrpc.h:
  	touch /tmp/dummy1_xmlrpc.h
        AC_CHECK_HEADERS([/tmp/dummy1_xmlrpc.h], [ac_have_xmlrpch=yes], [ac_have_xmlrpch=no],
		[#include "xmlrpc_client.h"])
	rm /tmp/dummy1_xmlrpc.h
 	if test $ac_have_xmlrpch = yes; then
	        AC_SEARCH_LIBS(xmlrpc_env_init, xmlrpc, [ac_have_xmlrpc=yes], [ac_have_xmlrpc=no],
			[$XMLRPC_LIBS $XMLRPC_LDFLAGS])
	fi
# List of places to try
testdirs="$HOME/opt/xmlrpc $OBITINSTALL/other"
for dir in $testdirs; do
	if test $ac_have_xmlrpc = no; then
		if  test -f $dir/include/xmlrpc_client.h; then
			XMLRPC_CFLAGS="-I$dir/include"
			CPPFLAGS="$ac_xmlrpc_saved_CPPFLAGS $XMLRPC_CFLAGS"
			XMLRPC_LDFLAGS="-L$dir/lib"
			LDFLAGS="$XMLRPC_LDFLAGS $ac_xmlrpc_saved_LDFLAGS "
			#echo "Trying in $dir with $LDFLAGS"
  			touch /tmp/dummy3_xmlrpc.h
		        AC_CHECK_HEADERS([/tmp/dummy3_xmlrpc.h], [ac_have_xmlrpch=yes], [ac_have_xmlrpch=no],
				[#include "xmlrpc_client.h"])
			rm /tmp/dummy3_xmlrpc.h
			if test $ac_have_xmlrpch = yes; then
				# Force check
				ac_cv_search_xmlrpc_env_init="  "
			        AC_SEARCH_LIBS(xmlrpc_env_init, xmlrpc, [ac_have_xmlrpc=yes], 
					[ac_have_xmlrpc=no],
					[$XMLRPC_LIBS $XMLRPC_LDFLAGS])
				fi
			if test $ac_have_xmlrpc = yes ; then
				if test $ac_have_xmlrpch = yes ; then
					break;
				fi
			fi
		fi
	fi
done[]
XMLRPC_SERVER_CPPFLAGS="$CPPFLAGS $XMLRPC_CPPFLAGS"
XMLRPC_SERVER_LDFLAGS="$LDFLAGS $XMLRPC_LDFLAGS"
XMLRPC_SERVER_LIBS="$ac_xmlrpc_saved_LIBS $XMLRPC_SERVER_LIBS "

XMLRPC_CLIENT_CPPFLAGS="$CPPFLAGS $XMLRPC_CPPFLAGS"
XMLRPC_CLIENT_LDFLAGS="$LDFLAGS $XMLRPC_LDFLAGS"
XMLRPC_CLIENT_LIBS="$ac_xmlrpc_saved_LIBS -lxmlrpc_client -lcurl $XMLRPC_LIBS "
XMLRPC_LIBS="$ac_xmlrpc_saved_LIBS $XMLRPC_LIBS "

if test $ac_have_xmlrpch = no; then
	AC_MSG_WARN([cannot find XMLRPC headers])
	ac_have_xmlrpc=no
fi
if test $ac_have_xmlrpc = no; then
	AC_MSG_WARN([cannot find XMLRPC library])
	XMLRPC_SERVER_CPPFLAGS=""
	XMLRPC_SERVER_LDFLAGS=""
	XMLRPC_SERVER_LIBS=""

	XMLRPC_CLIENT_CPPFLAGS=""
	XMLRPC_CLIENT_LDFLAGS=""
	XMLRPC_CLIENT_LIBS=""
fi
AC_MSG_CHECKING(checking if xmlrpc available )
AC_MSG_RESULT($ac_have_xmlrpc)
if test $ac_have_xmlrpc = yes; then
	AC_DEFINE(HAVE_XMLRPC, 1, [Define to 1 if XMLRPC is available.])
fi

AC_SUBST(XMLRPC_SERVER_CPPFLAGS)
AC_SUBST(XMLRPC_SERVER_LDFLAGS)
AC_SUBST(XMLRPC_SERVER_LIBS)
AC_SUBST(XMLRPC_CLIENT_CPPFLAGS)
AC_SUBST(XMLRPC_CLIENT_LDFLAGS)
AC_SUBST(XMLRPC_CLIENT_LIBS)

CPPFLAGS="$ac_xmlrpc_saved_CPPFLAGS"
LDFLAGS="$ac_xmlrpc_saved_LDFLAGS"
LIBS="$ac_xmlrpc_saved_LIBS"
	AC_SUBST(XMLRPC_CPPFLAGS)
	AC_SUBST(XMLRPC_LDFLAGS)
	AC_SUBST(XMLRPC_LIBS)
	AC_SUBST(WWWLIB_WL_RPATH)
	AC_SUBST(WWW_LIBS)
])
