# serial 1 ibverbs.m4
AC_DEFUN([AX_CHECK_IBVERBS], [
  AC_PREREQ([2.65])dnl

  AC_CHECK_HEADER(infiniband/verbs.h,
              # Found
              AC_MSG_NOTICE([IBverbs will be enabled.])
              AM_CONDITIONAL(IBVERBS_ENABLED, true),
              # Not found there, error
              AC_MSG_NOTICE([IBverbs will be disabled.])
              AM_CONDITIONAL(IBVERBS_ENABLED, false)
  )

])
