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
##	Filename:  make.mk
##	Author:    Dominik Cebulla
##	Version:   1.0 Beta
##	Date:      2017
##

TOP = $(realpath $(dir $(lastword $(MAKEFILE_LIST))))

include ${TOP}/make_linux.mk
