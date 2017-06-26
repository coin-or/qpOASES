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
 *	\file include/qpPresolver/types.h
 *	\author Dominik Cebulla
 *	\version 1.0 Beta
 *	\date 2017
 *
 *  Definition of several types for qpPresolver.
 */


#ifndef QPPRESOLVER_TYPES_H
#define QPPRESOLVER_TYPES_H

#include <stdint.h>

#include <qpPresolver/return_values.h>


#ifdef QPP_USE_SINGLE_PRECISION
    typedef float qpp_real_t;       /**< We use single precision as a real type. */
#else
    typedef double qpp_real_t;      /**< We use double precision as a real type. */
#endif


#ifdef QPP_USE_INT64
    //typedef int_least64_t qpp_int_t;    /**< We use integers with at least 64 bit length. */
    typedef long int qpp_int_t;             /**< Integer type. */
    typedef unsigned long int qpp_uint_t;   /**< Unsigned integer type. */
#else
    //typedef int_least32_t qpp_int_t;    /**< We use integers with at least 32 bit length. */
    typedef int qpp_int_t;                  /**< Integer type. */
    typedef unsigned int qpp_uint_t;        /**< Unsigned integer type. */
#endif


/** Boolean type. */
typedef enum bool_type
{
    QPP_BT_FALSE = 0,
    QPP_BT_TRUE = 1
} qpp_bool_type_t;


/** Indicates which type of bounds are passed to the user after presolving. */
typedef enum bound_mode
{
    QPP_BM_TIGHTEST_BOUNDS = 1,    /**< Tightest bounds are passed to the user. */
    QPP_BM_MEDIUM_BOUNDS = 2       /**< Medium bounds (as tight as necessary, but not looser
                                        than original bounds) are passed to the user. */
} qpp_bound_mode_t;


/** Indicates if (and which) new bounds have been derived on a variable. */
typedef enum tighter_bounds_type
{
    QPP_TBT_NONE = 0,       /**< A preprocessing method led to no new bounds on a variable. */
    QPP_TBT_LOWER = 1,      /**< A preprocessing method led to new lower bounds on a variable. */
    QPP_TBT_UPPER = 2,      /**< A preprocessing method led to new upper bounds on a variable. */
    QPP_TBT_BOTH = 3        /**< A preprocessing method led to both new lower and upper
                                 bounds on a variable. */
} qpp_tighter_bounds_type_t;


/** Determines how the presolved matrix (in coordinate format) shall be stored: Either
 *  row-wise or column-wise.
 *
 *  \todo This could (should?) be a part of ecrmatrix.
 */
typedef enum matrix_sort_type
{
    QPP_MST_ROW_WISE = -1,      /**< Row-wise storage. */
    QPP_MST_COLUMN_WISE = 1     /**< Column-wise storage. */
} qpp_matrix_sort_type_t;


/** Indicates if a (linear / bound) constraint is active or inactive. */
typedef enum active_type
{
    QPP_AT_LOWER_BOUND = -1,        /**< Variable / constraint is active at lower bound. */
    QPP_AT_INACTIVE = 0,            /**< Variable / constraint is inactive. */
    QPP_AT_UPPER_BOUND = 1,         /**< Variable / constraint is active at upper bound. */
    QPP_AT_EQ_CONSTRAINT            /**< Bound / Linear constraint is an equality constraint. */
} qpp_active_type_t;


/** Represents the current status / state of the presolver. */
typedef enum status
{
    QPP_STATUS_NEW = 100,           /**< Presolver entity is newly created. */
    QPP_STATUS_MEMORY_ALLOCATED,    /**< All necessary memory has been allocated. */
    QPP_STATUS_DATA_SET,            /**< Data of QP is set. */
    QPP_STATUS_PRESOLVED,           /**< QP has been presolved. */
    QPP_STATUS_PRESOLVED_QP_SET,    /**< The data of the presolved QP has been computed. */
    QPP_STATUS_POSTSOLVED           /**< The QP has been postsolved, i.e. a primal-dual optimal
                                         solution of the original QP has been computed. */
} qpp_status_t;

#endif /* QPPRESOLVER_TYPES_H */
