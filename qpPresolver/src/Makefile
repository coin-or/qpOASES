##
##	This file is part of qpPresolver.
##
##	qpPresolver -- An implementation of presolving (= preprocessing) techniques
##  for Quadratic Programming.
##	Copyright (C) 2017 by Dominik Cebulla et al. All rights reserved.
##
##	qpPresolver is free software; you can redistribute it and/or
##	modify it under the terms of the GNU Lesser General Public
##	License as published by the Free Software Foundation; either
##	version 2.1 of the License, or (at your option) any later version.
##
##	qpPresolver is distributed in the hope that it will be useful,
##	but WITHOUT ANY WARRANTY; without even the implied warranty of
##	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. 
##	See the GNU Lesser General Public License for more details.
##
##	You should have received a copy of the GNU Lesser General Public
##	License along with qpPresolver; if not, write to the Free Software
##	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
##

##
##	Filename:  src/Makefile
##	Author:    Dominik Cebulla
##	Version:   1.0 Beta
##	Date:      2017
##

include ../make.mk

QPPRESOLVER_EXTERN_DIR   = ../extern
QPPRESOLVER_LIB_DIR      = ../lib
QPPRESOLVER_INCLUDE_DIR  = ../include

SUITESPARSE_INCLUDE_DIR  = /usr/include/suitesparse

QPPRESOLVER_IFLAGS       = -I../include \
                           -I${QPPRESOLVER_EXTERN_DIR}/mmio/include \
                           -I${SUITESPARSE_INCLUDE_DIR}


QPPRESOLVER_DEPENDS = ${QPPRESOLVER_INCLUDE_DIR}/qpPresolver.h \
                      ${QPPRESOLVER_INCLUDE_DIR}/qpPresolver/cholmod_config.h \
                      ${QPPRESOLVER_INCLUDE_DIR}/qpPresolver/constants.h \
                      ${QPPRESOLVER_INCLUDE_DIR}/qpPresolver/ecrmatrix.h \
                      ${QPPRESOLVER_INCLUDE_DIR}/qpPresolver/errorcodes.h \
                      ${QPPRESOLVER_INCLUDE_DIR}/qpPresolver/minheap.h \
                      ${QPPRESOLVER_INCLUDE_DIR}/qpPresolver/options.h \
                      ${QPPRESOLVER_INCLUDE_DIR}/qpPresolver/presolver.h \
                      ${QPPRESOLVER_INCLUDE_DIR}/qpPresolver/stack.h \
                      ${QPPRESOLVER_INCLUDE_DIR}/qpPresolver/types.h \
                      ${QPPRESOLVER_INCLUDE_DIR}/qpPresolver/utility.h \
                      ${QPPRESOLVER_EXTERN_DIR}/mmio/include/mmio.h


QPPRESOLVER_OBJECTS = ecrmatrix.${OBJEXT} \
	                  minheap.${OBJEXT} \
	                  options.${OBJEXT} \
	                  presolver.${OBJEXT} \
	                  stack.${OBJEXT} \
	                  utility.${OBJEXT}

QPPRESOLVER_EXTERN_OBJECTS = ${QPPRESOLVER_EXTERN_DIR}/mmio/src/mmio.${OBJEXT}


QPPRESOLVER_LIB         = ${QPPRESOLVER_LIB_DIR}/libqpPresolver.${LIBEXT}
QPPRESOLVER_SHARED_LIB  = ${QPPRESOLVER_LIB_DIR}/libqpPresolver.${DLLEXT}

# compiler command for generating objects
C = ${CC} ${CFLAGS} ${DEFINES} ${QPPRESOLVER_IFLAGS}


all: create_lib_dir ${QPPRESOLVER_LIB} ${QPPRESOLVER_SHARED_LIB}
    
create_lib_dir:
	${MKDIR} ${QPPRESOLVER_LIB_DIR}


%.${OBJEXT}: %.c
	@${C} -c $< ${DEF_TARGET}


${QPPRESOLVER_LIB}: ${QPPRESOLVER_OBJECTS} ${QPPRESOLVER_EXTERN_OBJECTS}
	@${ECHO} "Creating static qpPresolver lib."
	@${AR} rs $@ $^


${QPPRESOLVER_SHARED_LIB}: ${QPPRESOLVER_OBJECTS} ${QPPRESOLVER_EXTERN_OBJECTS}
	@${ECHO} "Creating shared qpPresolver lib."
	@${C} ${DEF_TARGET} ${SHARED} $^ ${LINK_LIBRARIES}


clean:
	@${ECHO} "Removing all qpPresolver/src/*.o." 
	@${RM} ${QPPRESOLVER_OBJECTS}

# ==============================================================================
# remove object and library files
# ==============================================================================
purge: clean
	@${ECHO} "Removing (shared and static) qpPresolver libraries."
	@${RM} ${QPPRESOLVER_LIB} ${QPPRESOLVER_SHARED_LIB} 
