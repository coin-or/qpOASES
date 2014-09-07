##
##	This file is part of qpOASES.
##
##	qpOASES -- An Implementation of the Online Active Set Strategy.
##	Copyright (C) 2007-2014 by Hans Joachim Ferreau, Andreas Potschka,
##	Christian Kirches et al. All rights reserved.
##
##	qpOASES is free software; you can redistribute it and/or
##	modify it under the terms of the GNU Lesser General Public
##	License as published by the Free Software Foundation; either
##	version 2.1 of the License, or (at your option) any later version.
##
##	qpOASES is distributed in the hope that it will be useful,
##	but WITHOUT ANY WARRANTY; without even the implied warranty of
##	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##	See the GNU Lesser General Public License for more details.
##
##	You should have received a copy of the GNU Lesser General Public
##	License along with qpOASES; if not, write to the Free Software
##	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
##



##
##	Filename:  testing/Makefile
##	Author:    Hans Joachim Ferreau
##	Version:   3.0
##	Date:      2007-2014
##



include ../make.mk

##
##	flags
##

IDIR        =  ../include

IFLAGS      =  -I./cpp \
               -I${IDIR}

QPOASES_TEST_EXES = \
	${BINDIR}/test_bench${EXE} \
	${BINDIR}/test_matrices${EXE} \
	${BINDIR}/test_matrices2${EXE} \
	${BINDIR}/test_matrices3${EXE} \
	${BINDIR}/test_indexlist${EXE} \
	${BINDIR}/test_example1${EXE} \
	${BINDIR}/test_example1a${EXE} \
	${BINDIR}/test_example1b${EXE} \
	${BINDIR}/test_example2${EXE} \
	${BINDIR}/test_example4${EXE} \
	${BINDIR}/test_example5${EXE} \
	${BINDIR}/test_example6${EXE} \
	${BINDIR}/test_example7${EXE} \
	${BINDIR}/test_exampleLP${EXE} \
	${BINDIR}/test_qrecipe${EXE} \
	${BINDIR}/test_gradientShift${EXE} \
	${BINDIR}/test_runAllOqpExamples${EXE} \
	${BINDIR}/test_sebastien1${EXE} \
	${BINDIR}/test_vanBarelsUnboundedQP${EXE} \
	${BINDIR}/test_janick1${EXE} \
	${BINDIR}/test_janick2${EXE} \
	${BINDIR}/test_constraintProduct1${EXE} \
	${BINDIR}/test_constraintProduct2${EXE} \
	${BINDIR}/test_guessedWS1${EXE}


##
##	targets
##

all: ${QPOASES_TEST_EXES}

runTests: ${QPOASES_TEST_EXES}
	./runUnitTests
	./checkForMemoryLeaks

${BINDIR}/test_bench${EXE}: ./cpp/test_bench.${OBJEXT} ${LINK_DEPENDS}
	@${ECHO} "Creating" $@
	@${CPP} ${DEF_TARGET} ${CPPFLAGS} $< ${QPOASES_LINK} ${LINK_LIBRARIES}

${BINDIR}/test_matrices${EXE}: ./cpp/test_matrices.${OBJEXT} ${LINK_DEPENDS}
	@${ECHO} "Creating" $@
	@${CPP} ${DEF_TARGET} ${CPPFLAGS} $< ${QPOASES_LINK} ${LINK_LIBRARIES}

${BINDIR}/test_matrices2${EXE}: ./cpp/test_matrices2.${OBJEXT} ./cpp/test_qrecipe_data.hpp ${LINK_DEPENDS}
	@${ECHO} "Creating" $@
	@${CPP} ${DEF_TARGET} ${CPPFLAGS} $< ${QPOASES_LINK} ${LINK_LIBRARIES}

${BINDIR}/test_matrices3${EXE}: ./cpp/test_matrices3.${OBJEXT} ./cpp/test_qrecipe_data.hpp ${LINK_DEPENDS}
	@${ECHO} "Creating" $@
	@${CPP} ${DEF_TARGET} ${CPPFLAGS} $< ${QPOASES_LINK} ${LINK_LIBRARIES}

${BINDIR}/test_indexlist${EXE}: ./cpp/test_indexlist.${OBJEXT} ${LINK_DEPENDS}
	@${ECHO} "Creating" $@
	@${CPP} ${DEF_TARGET} ${CPPFLAGS} $< ${QPOASES_LINK} ${LINK_LIBRARIES}

${BINDIR}/test_example1${EXE}: ./cpp/test_example1.${OBJEXT} ${LINK_DEPENDS}
	@${ECHO} "Creating" $@
	@${CPP} ${DEF_TARGET} ${CPPFLAGS} $< ${QPOASES_LINK} ${LINK_LIBRARIES}

${BINDIR}/test_example1a${EXE}: ./cpp/test_example1a.${OBJEXT} ${LINK_DEPENDS}
	@${ECHO} "Creating" $@
	@${CPP} ${DEF_TARGET} ${CPPFLAGS} $< ${QPOASES_LINK} ${LINK_LIBRARIES}

${BINDIR}/test_example1b${EXE}: ./cpp/test_example1b.${OBJEXT} ${LINK_DEPENDS}
	@${ECHO} "Creating" $@
	@${CPP} ${DEF_TARGET} ${CPPFLAGS} $< ${QPOASES_LINK} ${LINK_LIBRARIES}

${BINDIR}/test_example2${EXE}: ./cpp/test_example2.${OBJEXT} ${LINK_DEPENDS}
	@${ECHO} "Creating" $@
	@${CPP} ${DEF_TARGET} ${CPPFLAGS} $< ${QPOASES_LINK} ${LINK_LIBRARIES}

${BINDIR}/test_example4${EXE}: ./cpp/test_example4.${OBJEXT} ${LINK_DEPENDS}
	@${ECHO} "Creating" $@
	@${CPP} ${DEF_TARGET} ${CPPFLAGS} $< ${QPOASES_LINK} ${LINK_LIBRARIES}

${BINDIR}/test_example5${EXE}: ./cpp/test_example5.${OBJEXT} ${LINK_DEPENDS}
	@${ECHO} "Creating" $@
	@${CPP} ${DEF_TARGET} ${CPPFLAGS} $< ${QPOASES_LINK} ${LINK_LIBRARIES}

${BINDIR}/test_example6${EXE}: ./cpp/test_example6.${OBJEXT} ${LINK_DEPENDS}
	@${ECHO} "Creating" $@
	@${CPP} ${DEF_TARGET} ${CPPFLAGS} $< ${QPOASES_LINK} ${LINK_LIBRARIES}

${BINDIR}/test_example7${EXE}: ./cpp/test_example7.${OBJEXT} ${LINK_DEPENDS}
	@${ECHO} "Creating" $@
	@${CPP} ${DEF_TARGET} ${CPPFLAGS} $< ${QPOASES_LINK} ${LINK_LIBRARIES}

${BINDIR}/test_exampleLP${EXE}: ./cpp/test_exampleLP.${OBJEXT} ${LINK_DEPENDS}
	@${ECHO} "Creating" $@
	@${CPP} ${DEF_TARGET} ${CPPFLAGS} $< ${QPOASES_LINK} ${LINK_LIBRARIES}

${BINDIR}/test_qrecipe${EXE}: ./cpp/test_qrecipe.${OBJEXT} ./cpp/test_qrecipe_data.hpp ${LINK_DEPENDS}
	@${ECHO} "Creating" $@
	@${CPP} ${DEF_TARGET} ${CPPFLAGS} $< ${QPOASES_LINK} ${LINK_LIBRARIES}

${BINDIR}/test_gradientShift${EXE}: ./cpp/test_gradientShift.${OBJEXT} ${LINK_DEPENDS}
	@${ECHO} "Creating" $@
	@${CPP} ${DEF_TARGET} ${CPPFLAGS} $< ${QPOASES_LINK} ${LINK_LIBRARIES}

${BINDIR}/test_runAllOqpExamples${EXE}: ./cpp/test_runAllOqpExamples.${OBJEXT} ${LINK_DEPENDS}
	@${ECHO} "Creating" $@
	@${CPP} ${DEF_TARGET} ${CPPFLAGS} $< ${QPOASES_LINK} ${LINK_LIBRARIES}

${BINDIR}/test_sebastien1${EXE}: ./cpp/test_sebastien1.${OBJEXT} ${LINK_DEPENDS}
	@${ECHO} "Creating" $@
	@${CPP} ${DEF_TARGET} ${CPPFLAGS} $< ${QPOASES_LINK} ${LINK_LIBRARIES}

${BINDIR}/test_vanBarelsUnboundedQP${EXE}: ./cpp/test_vanBarelsUnboundedQP.${OBJEXT} ${LINK_DEPENDS}
	@${ECHO} "Creating" $@
	@${CPP} ${DEF_TARGET} ${CPPFLAGS} $< ${QPOASES_LINK} ${LINK_LIBRARIES}

${BINDIR}/test_janick1${EXE}: ./cpp/test_janick1.${OBJEXT} ${LINK_DEPENDS}
	@${ECHO} "Creating" $@
	@${CPP} ${DEF_TARGET} ${CPPFLAGS} $< ${QPOASES_LINK} ${LINK_LIBRARIES}

${BINDIR}/test_janick2${EXE}: ./cpp/test_janick2.${OBJEXT} ${LINK_DEPENDS}
	@${ECHO} "Creating" $@
	@${CPP} ${DEF_TARGET} ${CPPFLAGS} $< ${QPOASES_LINK} ${LINK_LIBRARIES}

${BINDIR}/test_constraintProduct1${EXE}: ./cpp/test_constraintProduct1.${OBJEXT} ${LINK_DEPENDS}
	@${ECHO} "Creating" $@
	@${CPP} ${DEF_TARGET} ${CPPFLAGS} $< ${QPOASES_LINK} ${LINK_LIBRARIES}

${BINDIR}/test_constraintProduct2${EXE}: ./cpp/test_constraintProduct2.${OBJEXT} ${LINK_DEPENDS}
	@${ECHO} "Creating" $@
	@${CPP} ${DEF_TARGET} ${CPPFLAGS} $< ${QPOASES_LINK} ${LINK_LIBRARIES}

${BINDIR}/test_guessedWS1${EXE}: ./cpp/test_guessedWS1.${OBJEXT} ${LINK_DEPENDS}
	@${ECHO} "Creating" $@
	@${CPP} ${DEF_TARGET} ${CPPFLAGS} $< ${QPOASES_LINK} ${LINK_LIBRARIES}


clean:
	@${ECHO} "Cleaning up (testing)"
	@${RM} -f ./cpp/*.${OBJEXT} ${QPOASES_TEST_EXES}

clobber: clean


cpp/test_matrices2.${OBJEXT}: ./cpp/test_matrices2.cpp ./cpp/test_qrecipe_data.hpp
	@${ECHO} "Creating" $@
	@${CPP} ${DEF_TARGET} ${IFLAGS} ${CPPFLAGS} -c $<

cpp/test_matrices3.${OBJEXT}: ./cpp/test_matrices3.cpp ./cpp/test_qrecipe_data.hpp
	@${ECHO} "Creating" $@
	@${CPP} ${DEF_TARGET} ${IFLAGS} ${CPPFLAGS} -c $<

cpp/test_qrecipe.${OBJEXT}: ./cpp/test_qrecipe.cpp ./cpp/test_qrecipe_data.hpp
	@${ECHO} "Creating" $@
	@${CPP} ${DEF_TARGET} ${IFLAGS} ${CPPFLAGS} -c $<

%.${OBJEXT}: %.cpp
	@${ECHO} "Creating" $@
	@${CPP} ${DEF_TARGET} ${IFLAGS} ${CPPFLAGS} -c $<


##
##	end of file
##
