##
##	Filename:  cutestDriver/Makefile
##	Author:    Dennis Janka
##

# include directories, relative
IDIR = ../../include
SRCDIR = ../../src
BINDIR = ../../bin
CUTESTLIBPATH = ${CUTEST}/objects/${MYARCH}/double

# don't change that (for now) ...
PROBLEMPATH = ./prob

# file extensions
CPP = g++
CPPFLAGS = -Wall -pedantic -Wshadow -O0 -finline-functions -fPIC -DLINUX -g
FF = gfortran
FFLAGS =

#link against CUTEst library, qpOASES library and compiled test problem
LINKOPTS = -L${CUTESTLIBPATH} -L${PROBLEMPATH} -Wl,-rpath=${PROBLEMPATH} \
		-L${BINDIR} -Wl,-rpath=${BINDIR} -lqpOASES \
		-lcutest -llapack -lm -lgfortran -lprob

IFLAGS      = -I${IDIR} \
	      -I${CUTEST}/include

QPOASES_EXES = \
	qpoasesCutest

##
##	targets
##

all: ${QPOASES_EXES}

qpoasesCutest: qpoasesCutest.o
	@echo "Creating" $@
	${CPP} -o $@ ${CPPFLAGS} $< ${LINKOPTS}

clean:
	rm -f *.o ${QPOASES_EXES}

clobber: clean

%.o: %.cpp
	@echo "Creating" $@
	@${CPP} -c ${IFLAGS} ${CPPFLAGS} $< -o $@

##
##	end of file
##
