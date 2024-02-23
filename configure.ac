#AC_INIT
AC_CONFIG_SRCDIR([src/grv98.f])
#AM_INIT_AUTOMAKE(pdflib,1.10)
AC_INIT([jam],[1.9095])
AC_PREFIX_DEFAULT($HOME)
AM_INIT_AUTOMAKE(dist-bzip2)
AC_PREFIX_DEFAULT($HOME)
AC_CONFIG_SRCDIR(src/jam.f)
AC_CONFIG_MACRO_DIR([m4])

#AC_PROG_LIBTOOL
AC_PROG_RANLIB

AC_PROG_CXX


prevF77="$F77"
AC_PROG_F77
AC_PROG_FC
if test -n "$prevF77"; then F77="$prevF77"; fi
AC_SUBST(FFLAGS)

FCFLAGS=${FFLAGS}

#$(FFLAGS) += "-mcmodel=medium"

#
# libtool
#
#AC_SUBST(LOADLIB_FLAG)
#AC_SUBST(LOADLIB_LIB)
#AC_SUBST(LOADLIB_LDFLAGS)
#AC_CHECK_HEADER(ltdl.h,
#	LOADLIB_FLAG="-DUSE_LOADLIB"
#	LOADLIB_LIB="-lltdl"
#	LOADLIB_LDFLAGS="-export-dynamic")
 
# CUDA
#AC_MSG_CHECKING(whether to use CUDA)
#AC_ARG_ENABLE(cuda, 
#[  --enable-cuda         compile with cuda],
#[case "${enableval}" in
#  yes)  cuda=true
#        AC_MSG_RESULT(yes)
#        ;;
#  no)   cuda=false
#        AC_MSG_RESULT(no)
#        ;;
#  *) AC_MSG_ERROR(bad value ${enableval} for --enable-cuda) ;;
#esac],
#[
#        cuda=false
#        AC_MSG_RESULT(no)
#])
#AM_CONDITIONAL(CUDA, test x$cuda = xtrue)
#
#if "$cuda" = "true"; then
#  NVCC="nvcc"
#  NVCCFLAGS="$NVCCFLAGS"
#  CUDA_LIBS="-lcudart"
#fi

AC_ARG_WITH([cuda],
   [  --with-cuda=PATH    prefix where cuda is installed [default=/opt/cuda]])

if test -n "$with_cuda"
then
   CUDA_CFLAGS="-I$with_cuda/include"
   CUDA_LIBS="-lcudart"
   CUDA_LDFLAGS="-L$with_cuda/lib$SUFFIX"
   NVCC="$with_cuda/bin/nvcc"
else
   CUDA_CFLAGS="-I/opt/cuda/include"
   CUDA_LIBS="-lcudart"
   CUDA_LDFLAGS="-L/opt/cuda/lib$SUFFIX"
   NVCC="nvcc"
fi
AC_SUBST(CUDA_CFLAGS)
AC_SUBST(CUDA_LIBS)
AC_SUBST(CUDA_LDFLAGS)
AC_SUBST(NVCC)
AM_CONDITIONAL(_CUDA, test x$with_cuda = xtrue)



AC_CONFIG_FILES([Makefile src/Makefile main/Makefile])
AC_OUTPUT
