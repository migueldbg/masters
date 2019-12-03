dnl -*- mode: autoconf -*- 
dnl
dnl Autoconf macro to resolve mysql library dependency
dnl Synopsis:
dnl
dnl  MYSQL_DEPS

AC_DEFUN([MYSQL_DEPS],
[
	DEPNAME=mysql
	echo "SEARCHING FOR \"${DEPNAME}\" ..."
	LIBCONFIG="${DEPNAME}_config"
	AC_MSG_CHECKING(${LIBCONFIG})
	if which "${LIBCONFIG}" > /dev/null 2> /dev/null ; then
		AC_MSG_RESULT(yes)
		DEP_CFLAGS="${DEP_CFLAGS} "`"${LIBCONFIG}" --include`
		DEP_LIBS="`${LIBCONFIG} --libs` ${DEP_LIBS}"
	else
		AC_MSG_ERROR([Could not find ]${DEPNAME}[ library!])
	fi
])


#
# EOF
#
