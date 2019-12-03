# Initialization

AC_INIT([ReD PADME-based DAQ], [0], [], [red-daq])

AC_CONFIG_MACRO_DIR([m4])

AM_INIT_AUTOMAKE([-Wall -Werror])


# Libtool init and language selection

# Workaround for newer vs. older versions of Automake:
m4_ifdef([AM_PROG_AR], [AM_PROG_AR])

AC_LIBTOOL_DLOPEN
AC_PROG_LIBTOOL
AC_PROG_CXX

CXX_PROJECT_INIT


# Checks for dependencies
ROOTSYS_DEPS([glibs], [-lSpectrum])
MYSQL_DEPS([])
CXX_PROJECT_DEPS([])


# Output

AC_CONFIG_HEADERS([config.h])

AC_CONFIG_FILES([
	Makefile
	red-daq.pc
	red-daq-config
	src/PadmeRoot/Makefile
	src/Modules/Makefile 
        src/Level1/Makefile        
])

AC_OUTPUT
