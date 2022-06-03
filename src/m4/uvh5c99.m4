# serial 1 uvh5c99.m4
AC_DEFUN([AX_CHECK_UVH5C99], [
  AC_PREREQ([2.65])dnl

  AC_ARG_WITH([uvh5],
            AC_HELP_STRING([--with-uvh5=DIR],
                           [Location of UVH5C99 install directory]),
            [
              UVH5C99DIR="$withval"
              has_uvh5=1
            ],
            [
              UVH5C99DIR=""
              has_uvh5=0
            ])


  if test $has_uvh5 = 0; then
    AC_MSG_NOTICE([Library UVH5C99 not provided. UVH5 will not be linked.])
    uvh5c99_enabled=0;
  else
    # test uvh5c99 before enabling

    AC_CHECK_FILE([${UVH5C99DIR}/include/uvh5.h],
                  # Found
                  AC_SUBST(UVH5C99_INCDIR,${UVH5C99DIR}/include),
                  # Not found there, check UVH5C99DIR
                  AC_CHECK_FILE([${UVH5C99DIR}/../include/uvh5.h],
                                # Found
                                AC_SUBST(UVH5C99_INCDIR,${UVH5C99DIR}/../include),
                                # Not found there, error
                                AC_MSG_ERROR([uvh5.h header file not found])))

    orig_LDFLAGS="${LDFLAGS}"
    LDFLAGS="${orig_LDFLAGS} -L${UVH5C99DIR}/lib"
    AC_CHECK_LIB([uvh5], [UVH5permute_uvws],
                # Found
                AC_SUBST(UVH5C99_LIBDIR,${UVH5C99DIR}/lib),
                # Not found there, check UVH5C99DIR
                AS_UNSET(ac_cv_lib_uvh5_UVH5permute_uvws)
                LDFLAGS="${orig_LDFLAGS} -L${UVH5C99DIR}"
                AC_CHECK_LIB([uvh5], [UVH5permute_uvws],
                            # Found
                            AC_SUBST(UVH5C99_LIBDIR,${UVH5C99DIR}),
                            # Not found there, error
                            AC_MSG_ERROR([UVH5C99 library not found])))
    LDFLAGS="${orig_LDFLAGS}"

    uvh5c99_enabled=1;
  fi

  AS_IF([test $uvh5c99_enabled = 1],
    [
      AM_CONDITIONAL(UVH5C99_ENABLED, true)
      AC_DEFINE(UVH5C99_ENABLED,[],[Use UVH5C99])
      AC_MSG_NOTICE([UVH5C99 will be enabled.])
    ],
    [
      AM_CONDITIONAL(UVH5C99_ENABLED, false)
    ]
  )
])
