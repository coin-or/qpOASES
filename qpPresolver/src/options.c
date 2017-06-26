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
 *	\file src/options.c
 *	\author Dominik Cebulla
 *	\version 1.0 Beta
 *	\date 2017
 *
 *  Implementation of an option interface for the presolver.
 */


#include <qpPresolver/options.h>

/* Default values */
const qpp_bound_mode_t  DEFAULT_BOUND_MODE = QPP_BM_MEDIUM_BOUNDS;
const qpp_real_t        DEFAULT_EQ_TOL = QPP_EPS * 100.0;

#ifndef QPP_USE_SINGLE_PRECISION
const qpp_real_t        DEFAULT_FEAS_TOL = 1e-8;
#else
const qpp_real_t        DEFAULT_FEAS_TOL = 1e-4;
#endif

const qpp_real_t        DEFAULT_STAB_TOL = 1e-3;
const qpp_int_t         DEFAULT_MAX_ITER = 30;
const qpp_int_t         DEFAULT_LOG_LEVEL = 3;

const qpp_bool_type_t   DEFAULT_ENABLE_SINGLETON_COLUMNS = QPP_BT_TRUE;
const qpp_bool_type_t   DEFAULT_ENABLE_BOUND_TIGHTENING = QPP_BT_TRUE;
const qpp_bool_type_t   DEFAULT_ENABLE_DUAL_METHODS = QPP_BT_TRUE;
const qpp_bool_type_t   DEFAULT_ENABLE_DUPLICATE_COLUMNS = QPP_BT_TRUE;
const qpp_bool_type_t   DEFAULT_ENABLE_EMPTY_COLUMNS = QPP_BT_TRUE;
const qpp_bool_type_t   DEFAULT_ENABLE_PRIMAL_METHODS = QPP_BT_TRUE;
const qpp_bool_type_t   DEFAULT_ENABLE_SCALING = QPP_BT_FALSE;
const qpp_bool_type_t   DEFAULT_ENABLE_SINGLETON_ROWS = QPP_BT_TRUE;
const qpp_bool_type_t   DEFAULT_ENABLE_SPARSIFICATION = QPP_BT_TRUE;


/*  ======================================================================================
    Implementation of interface functions.
    ====================================================================================*/

qpp_return_value_t qppSetDefaultOptions(qpp_options_t *const opt)
{
    if (opt == NULL)
    {
        return QPP_NULL_ARGUMENT;
    }

    /* Set standard options */
    opt->bound_mode     = DEFAULT_BOUND_MODE;
    opt->eq_tol         = DEFAULT_EQ_TOL;
    opt->feas_tol       = DEFAULT_FEAS_TOL;
    opt->stab_tol       = DEFAULT_STAB_TOL;
    opt->max_iter       = DEFAULT_MAX_ITER;
    opt->log_level      = DEFAULT_LOG_LEVEL;

    opt->enable_bound_tightening            = DEFAULT_ENABLE_BOUND_TIGHTENING;
    opt->enable_dual_constraints_method     = DEFAULT_ENABLE_DUAL_METHODS;
    opt->enable_duplicate_columns_method    = DEFAULT_ENABLE_DUPLICATE_COLUMNS;
    opt->enable_empty_columns_method        = DEFAULT_ENABLE_EMPTY_COLUMNS;
    opt->enable_primal_constraints_method   = DEFAULT_ENABLE_PRIMAL_METHODS;
    opt->enable_scaling                     = DEFAULT_ENABLE_SCALING;
    opt->enable_singleton_columns_method    = DEFAULT_ENABLE_SINGLETON_COLUMNS;
    opt->enable_singleton_rows_method       = DEFAULT_ENABLE_SINGLETON_ROWS;
    opt->enable_sparsification_method       = DEFAULT_ENABLE_SPARSIFICATION;

    return QPP_OK;
}



qpp_return_value_t qppCheckOptions(qpp_options_t *const opt)
{
    qpp_int_t invalid_value;
    invalid_value = 0;

    if (! ((opt->bound_mode == QPP_BM_TIGHTEST_BOUNDS) ||
           (opt->bound_mode == QPP_BM_MEDIUM_BOUNDS)) )
    {
        opt->bound_mode = DEFAULT_BOUND_MODE;
        invalid_value = 1;
    }

    if (opt->eq_tol <= 0.0)
    {
        opt->eq_tol = DEFAULT_EQ_TOL;
        invalid_value = 1;
    }

    if (opt->feas_tol <= 0.0)
    {
        opt->feas_tol = DEFAULT_FEAS_TOL;
        invalid_value = 1;
    }

    if (opt->stab_tol <= 0.0)
    {
        opt->stab_tol = DEFAULT_STAB_TOL;
        invalid_value = 1;
    }

    if (opt->max_iter <= 0)
    {
        opt->max_iter = DEFAULT_MAX_ITER;
        invalid_value = 1;
    }

    if ( (opt->log_level < 0) || (opt->log_level > 5) )
    {
        opt->log_level = DEFAULT_LOG_LEVEL;
        invalid_value = 1;
    }

    return (invalid_value == 0 ? QPP_OK : QPP_INVALID_OPTIONS);
}


qpp_return_value_t qppCopyOptions(qpp_options_t *const dest,
                                  const qpp_options_t *const src)
{
    if (dest == NULL)
    {
        return QPP_NULL_ARGUMENT;
    }

    if (src == NULL)
    {
        return qppSetDefaultOptions(dest);
    }
    else
    {
        dest->bound_mode    = src->bound_mode;
        dest->max_iter      = src->max_iter;
        dest->eq_tol        = src->eq_tol;
        dest->feas_tol      = src->feas_tol;
        dest->stab_tol      = src->stab_tol;
        dest->log_level     = src->log_level;

        dest->enable_bound_tightening           = src->enable_bound_tightening;
        dest->enable_dual_constraints_method    = src->enable_dual_constraints_method;
        dest->enable_duplicate_columns_method   = src->enable_duplicate_columns_method;
        dest->enable_empty_columns_method       = src->enable_empty_columns_method;
        dest->enable_primal_constraints_method  = src->enable_primal_constraints_method;
        dest->enable_scaling                    = src->enable_scaling;
        dest->enable_singleton_columns_method   = src->enable_singleton_columns_method;
        dest->enable_singleton_rows_method      = src->enable_singleton_rows_method;
        dest->enable_sparsification_method      = src->enable_sparsification_method;
    }
    return QPP_OK;
}


void qppPrintOptions(const qpp_options_t *const opt)
{
    fprintf(stdout, "\n========== QP Presolver Options ==========\n\n");
    fprintf(stdout, "Print level: %" QPP_PRID "\n", opt->log_level);

    if (opt->bound_mode == QPP_BM_TIGHTEST_BOUNDS)
    {
        fprintf(stdout, "Tightest Bounds will be used\n");
    }
    else if (opt->bound_mode == QPP_BM_MEDIUM_BOUNDS)
    {
        fprintf(stdout, "Medium Bounds will be used\n");
    }
    else
    {
        fprintf(stdout, "Invalid Bound Type! Check your settings!\n");
    }

    fprintf(stdout, "Equality  Tolerance: %.6e\n", opt->eq_tol);
    fprintf(stdout, "Feasible  Tolerance: %.6e\n", opt->feas_tol);
    fprintf(stdout, "Stability Tolerance: %.6e\n", opt->stab_tol);
    fprintf(stdout, "Maximum number of iterations: %" QPP_PRID "\n", opt->max_iter);

    fprintf(stdout, "\nEnabled Preprocessing Methods: enable_\n");
    fprintf(stdout, "\tbound_tightening:          %d\n",
            opt->enable_bound_tightening);
    fprintf(stdout, "\tdual_constraints_method:   %d\n",
            opt->enable_dual_constraints_method);
    fprintf(stdout, "\tduplicate_columns_method:  %d\n",
            opt->enable_duplicate_columns_method);
    fprintf(stdout, "\tempty_columns_method:      %d\n",
            opt->enable_empty_columns_method);
    fprintf(stdout, "\tprimal_constraints_method: %d\n",
            opt->enable_primal_constraints_method);
    fprintf(stdout, "\tscaling:                   %d\n",
            opt->enable_scaling);
    fprintf(stdout, "\tsingleton_columns_method:  %d\n",
            opt->enable_singleton_columns_method);
    fprintf(stdout, "\tsingleton_rows_method:     %d\n",
            opt->enable_singleton_rows_method);
    fprintf(stdout, "\tsparsification_method:     %d\n",
            opt->enable_sparsification_method);
}
