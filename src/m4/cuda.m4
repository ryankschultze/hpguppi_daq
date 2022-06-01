# serial 1 cuda.m4
AC_DEFUN([AX_CHECK_CUDA], [
  AC_PREREQ([2.65])dnl

  AC_ARG_WITH([cuda-include],
              AC_HELP_STRING([--with-cuda-include=DIR],
                            [Location of CUDA include directory]),
              [
                CUDA_INCDIR="$withval"
                has_cuda_include=1
              ],
              [
                CUDA_INCDIR=""
                has_cuda_include=0
              ])

  if test $has_cuda_include = 1; then
    AC_SUBST(CUDA_INCDIR,${CUDA_INCDIR})
    AC_MSG_NOTICE([CUDA include directory set to ${CUDA_INCDIR}])
    cuda_enabled=1;
  else
    cuda_enabled=0;
  fi

  AS_IF([test $cuda_enabled = 1],
    [
      AM_CONDITIONAL(CUDA_ENABLED, true)
      AC_DEFINE(CUDA_ENABLED,[],[Use CUDA])
    ],
    [
      AM_CONDITIONAL(CUDA_ENABLED, false)
    ])

])
