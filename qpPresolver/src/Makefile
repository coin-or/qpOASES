# Makefile for creating the qpPresolver library.

EXTDIR = ../extern

LIBDIR = ../lib
QPPDIR = ../include/qpPresolver

MKDIR_P = mkdir -p

# ==============================================================================
# header files for dependencies
# ==============================================================================
HDR = ../include/qpPresolver.h \
      $(QPPDIR)/constants.h \
      $(QPPDIR)/ecrmatrix.h \
      $(QPPDIR)/errorcodes.h \
      $(QPPDIR)/minheap.h \
      $(QPPDIR)/options.h \
      $(QPPDIR)/presolver.h \
      $(QPPDIR)/stack.h \
      $(QPPDIR)/types.h \
      $(QPPDIR)/utility.h \
      $(EXTDIR)/mmio/include/mmio.h
      
INCLUDES = $(EXTDIR)/mmio/include 

# ==============================================================================
# object files
# ==============================================================================
OBJ = $(LIBDIR)/ecrmatrix.o \
	  $(LIBDIR)/minheap.o \
	  $(LIBDIR)/options.o \
	  $(LIBDIR)/presolver.o \
	  $(LIBDIR)/stack.o \
	  $(LIBDIR)/utility.o

EXTERN_OBJ = $(EXTDIR)/mmio/lib/mmio.o

# ==============================================================================
# library files
# ==============================================================================
LIB = $(LIBDIR)/libqpPresolver.a
SLIB = $(LIBDIR)/libqpPresolver.so

# ==============================================================================
# compiler command for generating objects
# ==============================================================================
C = $(CC) $(CFLAGS) $(DEFINES) -I../include -I$(INCLUDES)

all: $(LIB) $(SLIB)
    
# ==============================================================================
# compile all object files
# ==============================================================================
$(LIBDIR)/%.o: %.c
	@$(MKDIR_P) $(LIBDIR)
	@$(C) -c $< -o $@

# ==============================================================================
# build static library
# ==============================================================================
$(LIB): $(OBJ) $(EXTERN_OBJ)
	@$(AR) rs $@ $^

# ==============================================================================
# build shared library
# ==============================================================================
$(SLIB): $(OBJ) $(EXTERN_OBJ)
	@$(C) -o $@ -shared $^

# ==============================================================================
# remove object files
# ==============================================================================
clean:
#	@echo "Cleaning up (lib/*.o)" 
	@$(RM) $(OBJ)

# ==============================================================================
# remove object and library files
# ==============================================================================
purge: clean
#	@echo "Removing (shared and static) libraries"
	@$(RM) $(LIB) $(SLIB)
