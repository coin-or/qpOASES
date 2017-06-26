# Makefile for compiling the qpPresolver library and its dependencies.
# Use only this Makefile or the installation scripts!!!

# ==============================================================================
# integer type
# ==============================================================================
#
# Per default QPP_USE_INT64 is not defined. Call
#
# $ make DEFINES="-DQPP_USE_INT64"
#
# to compile with 64 bit integers.

# ==============================================================================
# single precision
# ==============================================================================
#
# Per default QPP_USE_SINGLE_PRECISION is not defined. Call
#
# $ make DEFINES="-DQPP_USE_SINGLE_PRECISION"
#
# to use single precision floating point numbers.


# ==============================================================================
# C compiler and flags
#
# CFLAGS are additional compiler flags not given to external code
# ==============================================================================

CC = gcc
CFLAGS += -Wall -Wextra -std=c99 -pedantic -Wno-unused-variable -Wshadow \
          -Wno-unused-function -finline-functions -Wno-float-equal \
          -Wno-sign-conversion -Wno-conversion -fPIC

export CC
export CFLAGS
export DEFINES

.PHONY: all mmio src clean purge


# ==============================================================================
# compile libraries
# ==============================================================================
all: src mmio
	
src: mmio
	@cd $@; $(MAKE) -s;

mmio:
	@cd extern/$@; $(MAKE) -s;

# ==============================================================================
# clean removes the object files but keeps the libraries
# ==============================================================================
clean:
	@cd extern/mmio; $(MAKE) -s clean;
	@cd src; $(MAKE) -s clean;


# ==============================================================================
# purge removes object files and library files
# ==============================================================================
purge:
	@cd extern/mmio; $(MAKE) -s purge;
	@cd src; $(MAKE) -s purge;
