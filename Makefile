# -----------------------------------------------------------------
# $Revision: 4419 $
# $Date: 2015-03-06 11:49:01 -0800 (Fri, 06 Mar 2015) $
# -----------------------------------------------------------------
# Programmer: Radu Serban @ LLNL
# -----------------------------------------------------------------
# Copyright (c) 2002, The Regents of the University of California.
# Produced at the Lawrence Livermore National Laboratory.
# All rights reserved.
# For details, see the LICENSE file.
# -----------------------------------------------------------------
# Makefile for IDA serial examples.
#
# This file is generated from a template using various variables
# set at configuration time. It can be used as a template for
# other user Makefiles.
#
# Note: if the solver was successfully configured with Blas/Lapack
# support, the Blas/Lapack libraries are specified through the 
# variable LIBRARIES_BL. Otherwise, this variable should contain
# an empty string. We include LIBRARIES_BL in the link line for
# all examples, whether they use the Lapack module or not, to
# address the case in which the SUNDIALS libraries are shared 
# objects. In that case, the solver library references Lapack 
# symbols which must be always resolved by linking against the
# Blas/Lapack libraries. If only static SUNDIALS libraries have 
# been built, it is not required to link the Blas/Lapack libraries
# for examples that do not use that module...
# -----------------------------------------------------------------

SHELL = sh

prefix       = /usr/local
exec_prefix  = /usr/local
includedir   = /usr/local/include
libdir       = /usr/local/lib

CPP      = /usr/bin/cc
CPPFLAGS = -O3 -DNDEBUG
CC       = /usr/bin/cc
CFLAGS   = -O3 -DNDEBUG
LDFLAGS  = 
LIBS     =  -lm /usr/lib/x86_64-linux-gnu/librt.so

LINKFLAGS = -Wl,-rpath,/usr/local/lib

INCLUDES = -I${includedir}
LIBRARIES = -lsundials_ida -lsundials_nvecserial ${LIBS}
LIBRARIES_BL =  /usr/lib/liblapack.so.3gf /usr/lib/libblas.so.3gf
LIBRARIES_SLUMT = 
LIBRARIES_KLU =    

EXAMPLES =  idaHeat1D_bnd

OBJECTS = ${EXAMPLES:=.o}

# -----------------------------------------------------------------------------------------

.SUFFIXES : .o .c

.c.o :
	${CC} ${CPPFLAGS} ${CFLAGS} ${INCLUDES} -c $<

# -----------------------------------------------------------------------------------------

all: ${OBJECTS}
	@for i in ${EXAMPLES} ; do \
	  echo "${CC} -o $${i} $${i}.o ${CFLAGS} ${LDFLAGS} ${INCLUDES} -L${libdir} ${LIBRARIES} ${LIBRARIES_BL} ${LIBRARIES_SLUMT} ${LIBRARIES_KLU}" ${LINKFLAGS} ; \
	  ${CC} -o $${i} $${i}.o ${CFLAGS} ${LDFLAGS} ${INCLUDES} -L${libdir} ${LIBRARIES} ${LIBRARIES_BL} ${LIBRARIES_SLUMT} ${LIBRARIES_KLU} ${LINKFLAGS} ; \
	done

clean:
	rm -f ${OBJECTS}
	rm -f ${EXAMPLES}

# -----------------------------------------------------------------------------------------

