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
 *	\file include/qpPresolver/return_values.h
 *	\author Dominik Cebulla
 *	\version 1.0 Beta
 *	\date 2017
 *
 *  Definition of several error codes / return values.
 */


#ifndef QPPRESOLVER_ERRORCODES_H
#define QPPRESOLVER_ERRORCODES_H


typedef enum return_value
{
    QPP_FATAL_ERROR = -1,           /**< Presolver is in some invalid state
                                         (e.g. corrupted data). */
    QPP_OK = 0,                     /**< Successful operation. */
    QPP_IO_ERROR,                   /**< Some (not specified) I/O error occurred. */
    QPP_UNABLE_TO_OPEN_FILE,        /**< A file could not be opened. */
    QPP_UNABLE_TO_READ_FILE,        /**< A file could not be read. */
    QPP_UNABLE_TO_WRITE_FILE,       /**< Could not write into a file. */
    QPP_INVALID_ARGUMENT,           /**< Some invalid argument (e.g. input parameter) occurred. */
    QPP_NULL_ARGUMENT,              /**< A NULL pointer was passed where it was not allowed. */
    QPP_OUT_OF_MEMORY,              /**< Not enough memory available. */
    QPP_STATUS_ERROR,               /**< Interface functions were called in wrong order
                                         (e.g. \p qppPostsolve() before \p qppPresolve()). */
    QPP_INVALID_OPTIONS,            /**< Some options passed by the user were invalid. */
    QPP_INFEASIBLE = 30,            /**< QP is infeasible. */
    QPP_PRIMAL_INFEASIBLE,          /**< QP is primal infeasible. */
    QPP_DUAL_INFEASIBLE,            /**< QP is dual infeasible. */
    QPP_BOUNDS_INCOMPATIBLE,        /**< Bounds on the primal / dual variables were inconsistent. */
    QPP_UNBOUNDED,                  /**< QP is unbounded. */
    QPP_BAD_SCALING,                /**< Scaling was not successful. */
    QPP_MATRIX_COLUMN_LIMIT = 50,   /**< Could not add a new column to an \p ecrmatrix. */
    QPP_MATRIX_ROW_LIMIT,           /**< Could not add a new row to an \p ecrmatrix. */
    QPP_HEAP_ELEMENT_NOT_FOUND      /**< The search for a certain element in a \p minheap
                                         was unsuccessful. */
} qpp_return_value_t;


#endif /* QPPRESOLVER_ERRORCODES_H */
