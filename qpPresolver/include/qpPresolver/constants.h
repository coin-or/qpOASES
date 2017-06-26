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
 *	\file include/qpPresolver/constants.h
 *	\author Dominik Cebulla
 *	\version 1.0 Beta
 *	\date 2017
 *
 *  Definition of several values/constants for qpPresolver (primarily as macros).
 */


#ifndef QPPRESOLVER_CONSTANTS_H
#define QPPRESOLVER_CONSTANTS_H

#include <float.h>
#include <inttypes.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>

#include <qpPresolver/types.h>


#ifdef QPP_USE_SINGLE_PRECISION
    #define QPP_EPS FLT_EPSILON     /**< Machine epsilon in single precision. */
    #define QPP_ZERO 1.0e-14        /**< Constant to determine whether a certain value
                                         is sufficiently close to 0. */
#else
    #define QPP_EPS DBL_EPSILON     /**< Machine epsilon in double precision. */
    #define QPP_ZERO 1.0e-25        /**< Constant to determine whether a certain value
                                         is sufficiently close to 0. */
#endif

#ifdef QPP_USE_INT64
    //#define QPP_MIN_INT INT_LEAST64_MIN     /**< Minimum integer value taken by \p int_least64_t. */
    //#define QPP_MAX_INT INT_LEAST64_MAX     /**< Maximum integer value taken by \p int_least64_t. */
    //#define QPP_PRID PRIdLEAST64            /**< Specifier for formatted output for \p int_least64_t. */
    //#define QPP_SCND SCNdLEAST64            /**< Specifier for formatted input for \p int_least64_t. */
    #define QPP_MIN_INT LONG_MIN            /**< Minimum integer value taken by <tt> long int </tt>. */
    #define QPP_MAX_INT LONG_MAX            /**< Maximum integer value taken by <tt> long int </tt>. */
    #define QPP_PRID "ld"                   /**< Specifier for formatted output for <tt> long int </tt>. */
    #define QPP_SCND "ld"                   /**< Specifier for formatted input for <tt> long int </tt>. */
#else
    //#define QPP_MIN_INT INT_LEAST32_MIN     /**< Minimum integer value taken by \p int_least32_t. */
    //#define QPP_MAX_INT INT_LEAST32_MAX     /**< Maximum integer value taken by \p int_least32_t. */
    //#define QPP_PRID PRIdLEAST32            /**< Specifier for formatted output for \p int_least32_t. */
    //#define QPP_SCND SCNdLEAST32            /**< Specifier for formatted input for \p int_least32_t. */
    #define QPP_MIN_INT INT_MIN             /**< Minimum integer value taken by \p int. */
    #define QPP_MAX_INT INT_MAX             /**< Maximum integer value taken by \p int. */
    #define QPP_PRID "d"                    /**< Specifier for formatted output for \p int. */
    #define QPP_SCND "d"                    /**< Specifier for formatted input for \p int. */
#endif

#define QPP_INF INFINITY                    /**< Infinity in floating point arithmetics. */

#define QPP_MAX_STRING_LENGTH 512           /**< Maximum allowed string length, e.g.
                                                 when reading a file with \p fscanf(). */

/* Assure that __FILE__ and __LINE__ are defined! */
#ifndef __FILE__
  #define __FILE__ 0
#endif

#ifndef __LINE__
  #define __LINE__ 0
#endif


#endif /* QPPRESOLVER_CONSTANTS_H */
