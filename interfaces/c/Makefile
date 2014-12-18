##
##	This file is part of qpOASES.
##
##	qpOASES -- An Implementation of the Online Active Set Strategy.
##	Copyright (C) 2007-2013 by Hans Joachim Ferreau, Andreas Potschka,
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
##	Filename:  interfaces/c/Makefile
##	Author:    Hans Joachim Ferreau
##	Version:   3.0
##	Date:      2014
##



PATH_LIB = .
PATH_SRC = .
PATH_OBJ = .
PATH_QPOASES_INC = ../../include
PATH_QPOASES_SRC = ../../src
PATH_QPOASES_OBJ = ../../src


##
##	bins
##

CC       = gcc
CPP      = g++
LD       = g++
CP       = cp



##
##	flags
##

CFLAGS   = -I$(PATH_QPOASES_INC) -I$(PATH_QPOASES_SRC) -Wall -pedantic -Wshadow -O3 -finline-functions -fPIC
CPPFLAGS = -I$(PATH_QPOASES_INC) -I$(PATH_QPOASES_SRC) -Wall -pedantic -Wshadow -O3 -finline-functions -fPIC -DLINUX -D__SINGLE_OBJECT__ 
LDFLAGS  = -shared
LIBEXT   = .so


LIB = libqpOASES_wrapper$(LIBEXT)

OBJ = \
	qpOASES_wrapper.o



##
##	targets
##

all: $(PATH_LIB)/$(LIB)

$(PATH_OBJ)/%.o: $(PATH_SRC)/%.c
	$(CC) -c -o $(@) $(CFLAGS) $<

$(PATH_OBJ)/%.o: $(PATH_SRC)/%.cpp
	$(CPP) -c -o $(@) $(CPPFLAGS) $<

$(PATH_LIB)/$(LIB) : $(addprefix $(PATH_OBJ)/, $(OBJ))
	$(LD) $(LDFLAGS) $^ -o $@

clean:
	$(RM) $(PATH_OBJ)/*.o $(PATH_LIB)/$(LIB)



##
##	end of file
##
