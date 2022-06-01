# serial 1 blade.m4
AC_DEFUN([AX_CHECK_BLADE], [
  AC_PREREQ([2.65])dnl

  AC_ARG_WITH([blade],
            AC_HELP_STRING([--with-blade=DIR],
                           [Location of BLADE library]),
            [
              BLADEDIR="$withval"
              has_blade=1
            ],
            [
              BLADEDIR=""
              has_blade=0
            ])

  
  if test $has_blade = 0; then
    #AC_MSG_RESULT([no])
    AC_MSG_NOTICE([Library BLADE not provided. BLADE will not be enabled.])
    blade_enabled=0;
  else
    # test blade before enabling
    AC_CHECK_FILE([${BLADEDIR}/include/blade/pipelines/ata/mode_a.hh],
                  # Found
                  AC_SUBST(BLADE_INCDIR,${BLADEDIR}/include),
                  # Not found there: check BLADEDIR
                  AC_CHECK_FILE([${BLADEDIR}/../include/blade/pipelines/ata/mode_a.hh],
                                # Found
                                AC_SUBST(BLADE_INCDIR,${BLADEDIR}/../include),
                                # Not found there: error
                                AC_MSG_ERROR([mode_a.hh header file not found])
                  )
    )

    # orig_LDFLAGS="${LDFLAGS}"
    # LDFLAGS="${orig_LDFLAGS} -L${BLADEDIR}/lib/x86_64-linux-gnu"
    # AC_LANG_PUSH([C++])
    # AX_CXX_CHECK_LIB([blade], [Blade::Pipeline::compute],
    #             # Found
    #             AC_SUBST(BLADE_LIBDIR,${BLADEDIR}/lib/x86_64-linux-gnu),
    #             # Not found there, check BLADEDIR
    #             AS_UNSET(ac_cv_lib_blade_Blade__Pipeline__compute)
    #             LDFLAGS="${orig_LDFLAGS} -L${BLADEDIR}"
    #             AX_CXX_CHECK_LIB([blade], [Blade::Pipeline::compute],
    #                         # Found
    #                         AC_SUBST(BLADE_LIBDIR,${BLADEDIR}),
    #                         # Not found there, error
    #                         AC_MSG_ERROR([BLADE library not found])
    #             )
    # )
    # AC_LANG_POP([C++])
    # LDFLAGS="${LDFLAGS}"

    AC_SUBST(BLADE_LIBDIR,${BLADEDIR}/lib/x86_64-linux-gnu)

    blade_enabled=1;
  fi
  
  AS_IF([test $blade_enabled = 1],
  [
    AM_CONDITIONAL(BLADE_ENABLED, true)
    AC_DEFINE(BLADE_ENABLED,[],[Use BLADE])
    AC_MSG_NOTICE([BLADE will be enabled.])
  ],
  [
    AM_CONDITIONAL(BLADE_ENABLED, false)
  ])
])
