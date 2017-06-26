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
 *	\file include/qpPresolver/options.h
 *	\author Dominik Cebulla
 *	\version 1.0 Beta
 *	\date 2017
 *
 *  Declaration of an option interface for the presolver.
 */


#ifndef QPPRESOLVER_OPTIONS_H
#define QPPRESOLVER_OPTIONS_H

#include <stdio.h>
#include <stdlib.h>

#include <qpPresolver/constants.h>
#include <qpPresolver/return_values.h>
#include <qpPresolver/types.h>


/** \brief Data structure containing several parameters that the user can adjust. */
typedef struct options
{
    qpp_bound_mode_t bound_mode;    /**< Indicates which bounds (tightest or medium)
                                         are passed to the user. */
    qpp_int_t max_iter;             /**< Maximum number of presolver iterations. */
    qpp_real_t eq_tol;              /**< Tolerance when checking for equality of two numbers. */
    qpp_real_t feas_tol;            /**< Tolerance when checking for (in-)feasibility. */
    qpp_real_t stab_tol;            /**< Stability tolerance, used when a value must
                                         not be too small. */
    qpp_int_t log_level;            /**< Determines how much information is written into the
                                         logfile (value between 0 and 5).
                                         \todo Change to an enum. */

    /** \name Flags that determine which preprocessing methods are used. */
    /**@{*/
    qpp_bool_type_t enable_bound_tightening;          /**< Try to compute tighter bounds on
                                                           primal variables. */
    qpp_bool_type_t enable_dual_constraints_method;   /**< Determines if the presolver checks for
                                                           (weakly) forcing and infeasible
                                                           dual constraints (on the basis of the
                                                           implied bounds on a dual constraint). */
    qpp_bool_type_t enable_duplicate_columns_method;  /**< Use duplicate columns method. */
    qpp_bool_type_t enable_empty_columns_method;      /**< Use empty columns method. */
    qpp_bool_type_t enable_primal_constraints_method; /**< Determines if the presolver checks for
                                                           forcing, infeasible and redundant
                                                           primal constraints (on the basis of the
                                                           implied bounds on a primal constraint). */
    qpp_bool_type_t enable_scaling;                   /**< Determines if the QP will be scaled or not. */
    qpp_bool_type_t enable_singleton_columns_method;  /**< Use singleton columns methods. */
    qpp_bool_type_t enable_singleton_rows_method;     /**< Use singleton rows method. */
    qpp_bool_type_t enable_sparsification_method;     /**< Use the sparsification method, i.e.
                                                           try to eliminate further nonzeros of
                                                           the constraint matrix. */
    /**@}*/
} qpp_options_t;


/*  ======================================================================================
    Interface functions.
    ====================================================================================*/

/** \brief Sets options to default values. Can be used to initialize the options structure.
 *
 *  \param opt Valid pointer to \p options. All values are replaced by default values.
 *
 *  \return QPP_OK \n QPP_NULL_ARGUMENT
 */
qpp_return_value_t qppSetDefaultOptions(qpp_options_t *const opt);


/** \brief Checks if some of the given options are invalid.
 *
 *  \param opt Valid pointer to \p options which are checked regarding consistency.
 *
 *  \return QPP_OK \n QPP_INVALID_OPTIONS
 */
qpp_return_value_t qppCheckOptions(qpp_options_t *const opt);


/** \brief Copies options from \p src to \p dest.
 *
 *  \param dest Valid pointer to \p options whose values are replaced by the values of \p src.
 *  \param src Pointer to \p options. If src == NULL, then \p dest is initialized with the
 *      default options, otherwise \p dest will contain the same values as \p src.
 *
 *  \return QPP_OK \n QPP_NULL_ARGUMENT
 */
qpp_return_value_t qppCopyOptions(qpp_options_t *const dest,
                                  const qpp_options_t *const src);


/** \brief Prints the given options in text form on stdout.
 *
 *  \param Valid pointer to \p options.
 */
void qppPrintOptions(const qpp_options_t *const opt);

#endif /* QPPRESOLVER_OPTIONS_H */
