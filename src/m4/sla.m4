# serial 1 sla.m4
AC_DEFUN([AX_CHECK_SLA], [
  AC_PREREQ([2.65])dnl

  AC_ARG_WITH([sla-lib],
              AC_HELP_STRING([--with-sla-lib=DIR],
                            [Location of SLA library]),
              [LIBSLADIR="$withval"],
              [LIBSLADIR="."])

  orig_LDFLAGS="${LDFLAGS}"
  LDFLAGS="${orig_LDFLAGS} -L${LIBSLADIR}"
  AC_CHECK_LIB([sla], [sla_obs_],
              # Found
              AC_SUBST(LIBSLADIR,${LIBSLADIR}),
              # Not found there, error
              AC_MSG_ERROR([SLA library not found]))
  LDFLAGS="${orig_LDFLAGS}"
])
