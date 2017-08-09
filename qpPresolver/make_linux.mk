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
##	Filename:  make_linux.mk
##	Author:    Dominik Cebulla
##	Version:   1.0 Beta
##	Date:      2017
##


############################ USER CONFIGURATION ############################

# Compilers and commands:
AR      = ar
CD      = cd
CP      = cp
ECHO    = echo
MKDIR   = mkdir -p
RM      = rm -f
CPP     = g++
CC      = gcc
F77     = gfortran

# File extensions:
OBJEXT     = o
LIBEXT     = a
DLLEXT     = so
EXE        =
MEXOCTEXT  = mex
DEF_TARGET = -o $@
SHARED     = -shared

# 32 or 64 depending on target platform
BITS = $(shell getconf LONG_BIT)

# decide on MEX interface extension
ifeq ($(BITS), 32)
	MEXEXT = mexglx
else
	MEXEXT = mexa64
endif

## Compiler flags
CFLAGS += -Wall -Wextra -std=c99 -pedantic -Wno-unused-variable -Wshadow \
          -Wno-unused-function -finline-functions -Wno-conversion -Wno-sign-conversion \
          -fPIC

LINK_LIBRARIES = -lcholmod -lsuitesparseconfig -lumfpack
