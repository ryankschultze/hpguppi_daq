#                                               -*- Autoconf -*-
# Process this file with autoconf to produce a configure script.

AC_PREREQ([2.65])
AC_INIT([FULL-PACKAGE-NAME], [VERSION], [BUG-REPORT-ADDRESS])
AM_INIT_AUTOMAKE([foreign])
LT_INIT
AM_SILENT_RULES([yes])
AC_CONFIG_SRCDIR([hpguppi_databuf.h])
AC_CONFIG_HEADERS([config.h])
AC_CONFIG_MACRO_DIR([m4])

# Set CFLAGS to nothing if it is not set by the user.  This prevents AC_PROG_CC
# from setting the (supposedly reserved-for-the-user!) variable CFLAGS in
# Makefile, which prevents AM_CFLAGS in Makefile.am from setting an
# optimization level.  For more details, see
# http://lists.gnu.org/archive/html/autoconf/2006-04/msg00007.html
AS_VAR_SET_IF(CFLAGS,[],[CFLAGS=])

# Checks for programs.
AC_PROG_CC
AC_PROG_CXX

# Check processor flags
AX_CHECK_X86_FEATURES
AM_CONDITIONAL([HAVE_AVX512F],
               [test "x$ax_cv_gcc_x86_cpu_supports_avx512f" = xyes])

# Check for HASHPIPE
AX_CHECK_HASHPIPE

# Check for rawspec
AX_CHECK_RAWSPEC

# Checks for libraries.
#AC_CHECK_LIB([pthread], [pthread_create])
AC_CHECK_LIB([m], [cos])
#AC_CHECK_LIB([rt], [clock_gettime])
#AC_CHECK_LIB([z], [crc32])

AC_ARG_WITH([libsla],
            AC_HELP_STRING([--with-libsla=DIR],
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

AC_ARG_WITH([libblade],
            AC_HELP_STRING([--with-libblade=DIR],
                           [Location of BLADE library]),
            [
              LIBBLADEDIR="$withval"
              has_libblade=1
            ],
            [
              LIBBLADEDIR=""
              has_libblade=0
            ])

AC_ARG_WITH([blade-include],
            AC_HELP_STRING([--with-blade-include=DIR],
                           [Location of BLADE include directory]),
            [BLADEINCDIR="$withval"],
            [BLADEINCDIR="${LIBBLADEDIR}/../include"])

AC_ARG_WITH([cuda-include],
            AC_HELP_STRING([--with-cuda-include=DIR],
                           [Location of CUDA include directory]),
            [
              CUDAINCDIR="$withval"
              has_cuda_include=1
            ],
            [
              CUDAINCDIR=""
              has_cuda_include=0
            ])

AC_ARG_WITH([libxgpu],
            AC_HELP_STRING([--with-libxgpu=DIR],
                           [Location of xGPU library]),
            [
              LIBXGPUDIR="$withval"
              has_libxgpu=1
            ],
            [
              LIBXGPUDIR=""
              has_libxgpu=0
            ])

AC_ARG_WITH([xgpu-include],
            AC_HELP_STRING([--with-xgpu-include=DIR],
                           [Location of xGPU include directory]),
            [XGPUINCDIR="$withval"],
            [XGPUINCDIR="${LIBXGPUDIR}/../include"])

#AC_MSG_CHECKING([for libblade])
#AC_LANG_PUSH(C++)
#orig_LDFLAGS="$LDFLAGS"
#LDFLAGS="${orig_LDFLAGS} -L${LIBBLADEDIR}"
#orig_LIBS="$LIBS"
#LIBS="${orig_LIBS} -lblade"
#AC_TRY_LINK([#include "blade/manager.hh"], 
#        [Blade::Manager object;],
#        has_libblade=1,
#        has_libblade=0)
#AC_LANG_POP(C++)
#LIBS="${orig_LIBS}"
#LDFLAGS="${orig_LDFLAGS}"

if test $has_libblade = 0; then
  #AC_MSG_RESULT([no])
  AC_MSG_NOTICE([Library BLADE not provided. BLADE will not be enabled.])
  blade_enabled=1;
else
  #AC_MSG_RESULT([yes])
  AC_MSG_NOTICE([BLADE will be enabled.])
  AC_SUBST(LIBBLADEDIR,${LIBBLADEDIR})
  
  AC_SUBST(BLADEINCDIR,${BLADEINCDIR})
  AC_MSG_NOTICE([BLADE include directory set to ${BLADEINCDIR}])
  blade_enabled=1;
fi

if test $has_cuda_include = 1; then
  AC_SUBST(CUDAINCDIR,${CUDAINCDIR})
  AC_MSG_NOTICE([CUDA include directory set to ${CUDAINCDIR}])
  cuda_enabled=1;
fi

if test $has_libxgpu = 0; then
  #AC_MSG_RESULT([no])
  AC_MSG_NOTICE([Library xGPU not provided. xGPU will not be enabled.])
  xgpu_enabled=0;
else
  #AC_MSG_RESULT([yes])
  AC_MSG_NOTICE([xGPU will be enabled.])
  AC_SUBST(LIBXGPUDIR,${LIBXGPUDIR})
  
  AC_SUBST(XGPUINCDIR,${XGPUINCDIR})
  AC_MSG_NOTICE([xGPU include directory set to ${XGPUINCDIR}])
  xgpu_enabled=1;
fi

AS_IF([test $blade_enabled = 1],
  [
    AM_CONDITIONAL(BLADE_ENABLED, true)
    AC_DEFINE(BLADE_ENABLED,[],[Use BLADE])
  ],
  [
    AM_CONDITIONAL(BLADE_ENABLED, false)
  ])

AS_IF([test $cuda_enabled = 1],
  [
    AM_CONDITIONAL(CUDA_ENABLED, true)
    AC_DEFINE(CUDA_ENABLED,[],[Use CUDA])
  ],
  [
    AM_CONDITIONAL(CUDA_ENABLED, false)
  ])

AS_IF([test $xgpu_enabled = 1],
  [
    AM_CONDITIONAL(XGPU_ENABLED, true)
    AC_DEFINE(XGPU_ENABLED,[],[Use xGPU])
  ],
  [
    AM_CONDITIONAL(XGPU_ENABLED, false)
  ])

## Checks for header files.
#AC_CHECK_HEADERS([netdb.h stdint.h stdlib.h string.h sys/socket.h sys/time.h unistd.h zlib.h])

# Checks for typedefs, structures, and compiler characteristics.
AC_C_INLINE
AC_TYPE_INT32_T
AC_TYPE_INT64_T
AC_TYPE_OFF_T
AC_TYPE_SIZE_T
AC_TYPE_UINT32_T
AC_TYPE_UINT64_T
AC_TYPE_UINT8_T

## Checks for library functions.
#AC_FUNC_MALLOC
#AC_CHECK_FUNCS([clock_gettime memset socket crc32])

AC_CONFIG_FILES([Makefile])
AC_OUTPUT
