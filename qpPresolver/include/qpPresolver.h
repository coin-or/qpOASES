/*
 *	This file is part of qpPresolver.
 *
 *	qpPresolver -- An implementation of presolving (= preprocessing) techniques
 *  for Quadratic Programming.
 *	Copyright (C) 2017 by Dominik Cebulla et al. All rights reserved.
 *
 *	qpPresolver is free software; you can redistribute it and/or
 *	modify it under the terms of the GNU Lesser General Public
 *	License as published by the Free Software Foundation; either
 *	version 2.1 of the License, or (at your option) any later version.
 *
 *	qpPresolver is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *	See the GNU Lesser General Public License for more details.
 *
 *	You should have received a copy of the GNU Lesser General Public
 *	License along with qpPresolver; if not, write to the Free Software
 *	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */


/**
 *	\file include/qpPresolver.h
 *	\author Dominik Cebulla
 *	\version 1.0 Beta
 *	\date 2017
 *
 *  Includes all qpPresolver headers.
 */


#ifndef QPPRESOLVER_H
#define QPPRESOLVER_H


/* qpOASES defines: */

#ifdef __USE_SINGLE_PRECISION__
    #define QPP_USE_SINGLE_PRECISION
#endif

#ifdef __USE_LONG_INTEGERS__
    #define QPP_USE_INT64
#endif

#ifdef __SUPPRESSANYOUTPUT__
    #undef QPP_WRITE_LOGFILE
#endif

#include <qpPresolver/cholmod_config.h>
#include <qpPresolver/constants.h>
#include <qpPresolver/ecrmatrix.h>
#include <qpPresolver/minheap.h>
#include <qpPresolver/options.h>
#include <qpPresolver/presolver.h>
#include <qpPresolver/stack.h>
#include <qpPresolver/types.h>
#include <qpPresolver/utility.h>


#endif /* QPPRESOLVER_H */
