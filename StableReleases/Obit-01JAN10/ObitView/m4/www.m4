AC_DEFUN([AC_PATH_WWW], [
	AC_ARG_WITH(www,
                    AC_HELP_STRING([--with-www=DIR],
                                 [search for WWW in DIR/include and DIR/lib]),
                    [for dir in `echo "$withval" | tr : ' '`; do
    if test -d $dir/include; then
      WWW_CPPFLAGS="$WWW_CPPFLAGS -I$dir/include"
    fi
    if test -d $dir/lib; then
      WWW_LDFLAGS="$WWW_LDFLAGS -L$dir/lib"
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

ac_www_saved_CPPFLAGS="$CPPFLAGS"
ac_www_saved_LDFLAGS="$LDFLAGS"
ac_www_saved_LIBS="$LIBS"
CPPFLAGS="$CPPFLAGS $WWW_CPPFLAGS"
LDFLAGS="$LDFLAGS $WWW_LDFLAGS"
ac_have_www=yes
        AC_SEARCH_LIBS(HTXML_new, [wwwxml], [], [ac_have_www=no
                       AC_MSG_WARN([cannot find WWW library])])
        AC_SEARCH_LIBS(XmlUtf8Encode, [xmltok], [], [ac_have_www=no
                       AC_MSG_WARN([cannot find WWW library])])
        AC_SEARCH_LIBS(xmlparse, [], [], [ac_have_www=no
                       AC_MSG_WARN([cannot find WWW library])])
        AC_SEARCH_LIBS(wwwzip, [], [], [ac_have_www=no
                       AC_MSG_WARN([cannot find WWW library])])
        AC_SEARCH_LIBS(wwwinit, [], [], [ac_have_www=no
                       AC_MSG_WARN([cannot find WWW library])])
        AC_SEARCH_LIBS(wwwapp, [], [], [ac_have_www=no
                       AC_MSG_WARN([cannot find WWW library])])
        AC_SEARCH_LIBS(md5, [], [], [ac_have_www=no
                       AC_MSG_WARN([cannot find WWW library])])
        AC_SEARCH_LIBS(wwwhtml, [], [], [ac_have_www=no
                       AC_MSG_WARN([cannot find WWW library])])
        AC_SEARCH_LIBS(wwwtelnet, [], [], [ac_have_www=no
                       AC_MSG_WARN([cannot find WWW library])])
        AC_SEARCH_LIBS(wwwnews, [], [], [ac_have_www=no
                       AC_MSG_WARN([cannot find WWW library])])
        AC_SEARCH_LIBS(wwwhttp, [], [], [ac_have_www=no
                       AC_MSG_WARN([cannot find WWW library])])
        AC_SEARCH_LIBS(wwwmime, [], [], [ac_have_www=no
                       AC_MSG_WARN([cannot find WWW library])])
        AC_SEARCH_LIBS(wwwgopher, [], [], [ac_have_www=no
                       AC_MSG_WARN([cannot find WWW library])])
        AC_SEARCH_LIBS(wwwftp, [], [], [ac_have_www=no
                       AC_MSG_WARN([cannot find WWW library])])
        AC_SEARCH_LIBS(wwwfile, [], [], [ac_have_www=no
                       AC_MSG_WARN([cannot find WWW library])])
        AC_SEARCH_LIBS(wwwdir, [], [], [ac_have_www=no
                       AC_MSG_WARN([cannot find WWW library])])
        AC_SEARCH_LIBS(wwwcache, [], [], [ac_have_www=no
                       AC_MSG_WARN([cannot find WWW library])])
        AC_SEARCH_LIBS(wwwstream, [], [], [ac_have_www=no
                       AC_MSG_WARN([cannot find WWW library])])
        AC_SEARCH_LIBS(wwwmux, [], [], [ac_have_www=no
                       AC_MSG_WARN([cannot find WWW library])])
        AC_SEARCH_LIBS(wwwtrans, [], [], [ac_have_www=no
                       AC_MSG_WARN([cannot find WWW library])])
        AC_SEARCH_LIBS(wwwcore, [], [], [ac_have_www=no
                       AC_MSG_WARN([cannot find WWW library])])
        AC_SEARCH_LIBS(wwwutils, [], [], [ac_have_www=no
                       AC_MSG_WARN([cannot find WWW library])])
if test $ac_have_www = yes; then
	AC_DEFINE(HAVE_WWW, 1, [Define to 1 if WWW is available.])
fi
WWW_LIBS="$LIBS"
CPPFLAGS="$ac_www_saved_CPPFLAGS"
LDFLAGS="$ac_www_saved_LDFLAGS"
LIBS="$ac_www_saved_LIBS"
	AC_SUBST(WWW_CPPFLAGS)
	AC_SUBST(WWW_LDFLAGS)
	AC_SUBST(WWW_LIBS)
])
