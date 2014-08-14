#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.59])
AC_INIT([FBA], [1.0], [andi.dhroso@mail.missouri.edu])
AM_INIT_AUTOMAKE
AC_CONFIG_SRCDIR([src/Biomass.cpp])
AC_CONFIG_HEADERS([config.h])

# Checks for programs.
AC_PROG_CXX
AC_PROG_CC
AC_PROG_RANLIB
AC_PROG_INSTALL
AM_PROG_CC_C_O

# Checks for libraries.
AC_SEARCH_LIBS([glp_create_prob], [glpk], [glpk_exists=yes], [glpk_exists=no])
if test "x${glpk_exists}" = xno; then
	echo "----------------------------------------------------------------------"
    echo "Unable to find glpk on this system."
    echo "Building with internal version.    "
	echo ""

	echo "Note: To use an installed version of glpk" 
	echo " installed in a non-standard location, "
	echo " pass the following to the configure script:"
	echo ""

	echo "./configure CPPFLAGS=-I/opt/local/include LDFLAGS=-L/opt/local/lib"
	echo ""
	echo " Substitute \"/opt/local/include\"  and \"/opt/local/lib\"" 
	echo " your own path header file and lib"
	echo "----------------------------------------------------------------------"
fi

AC_ARG_ENABLE([internal-glpk],
  [AS_HELP_STRING([--enable-internal-glpk], 
    [Build FBA with internal GLPK @<:@default: no@:>@])],
  [internal_glpk=yes],
  [internal_glpk=${enableval}])


# Checks for header files.
AC_CHECK_HEADERS([float.h limits.h stdlib.h string.h unistd.h glpk.h])

# Checks for typedefs, structures, and compiler characteristics.
#AC_CHECK_HEADER_STDBOOL
AC_C_INLINE
AC_TYPE_PID_T
AC_TYPE_SIZE_T

# Checks for library functions.
AC_FUNC_ERROR_AT_LINE
AC_FUNC_FORK
AC_FUNC_MALLOC
AC_FUNC_STRTOD
AC_CHECK_FUNCS([floor memmove memset pow sqrt strchr strerror strrchr utime glp_create_prob])

if test "x${internal_glpk}" = xyes; then
	AM_CONDITIONAL([INTERNAL_GLPK], [test "yes" = "yes"])
	AM_CONDITIONAL([SYSTEM_GLPK], [test "yes" = "no"])
	internal_glpk="yes"
else
	if test "x${glpk_exists}" = xyes; then
		internal_glpk=no
	else
		internal_glpk=yes
	fi
	
	AM_CONDITIONAL([INTERNAL_GLPK], [test "${glpk_exists}" = "no"])
	AM_CONDITIONAL([SYSTEM_GLPK], [test "${glpk_exists}" = "yes"])
fi
 

AC_CONFIG_FILES([Makefile
                 glpk/Makefile
                 src/Makefile])
AC_OUTPUT

AC_MSG_RESULT([  GLPK library support   -- $glpk_exists])
AC_MSG_RESULT([  Use internal GLPK      -- $internal_glpk])