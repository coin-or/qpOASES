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
##	Filename:  interfaces/MATLAB/Makefile
##	Author:    Dominik Cebulla
##	Version:   1.0 Beta
##	Date:      2017
##

include ../../make.mk

# DEFINES will be passed to "main" Makefile (in main directory).
export DEFINES


all:
	@cd ${TOP}; ${MAKE} -s purge; ${MAKE} -s;


clean:
	@cd ${TOP}; ${MAKE} -s clean;


purge:
	@${RM} qpp.mex* qpp.o
	@cd ${TOP}; ${MAKE} -s purge;
