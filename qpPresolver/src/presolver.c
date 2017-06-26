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
 *	\file src/presolver.c
 *	\author Dominik Cebulla
 *	\version 1.0 Beta
 *	\date 2017
 *
 *  Implementation of the presolver's interface functions and all preprocessing
 *  techniques.
 *
 */


#include <qpPresolver/presolver.h>


/*  ======================================================================================
    Forward declaration of static (= private) functions
    ====================================================================================*/

/** \brief Computes the number of variables / constraints and nonzero elements of the
 *      presolved QP (+ number of nonzero diagonal elements of the presolved Hessian).
 *
 *  The computed dimensions are stored within the \p data struct. The number of nonzeros
 *  of the Hessian matrix corresponds to the lower triangular part.
 */
static qpp_return_value_t computeDimensionPresolvedQP(qpp_data_t *const data);

/** \brief Computes the presolved QP and stores it into \p data. */
static qpp_return_value_t computePresolvedQP(qpp_data_t *const data,
                                             const qpp_matrix_sort_type_t sort_type);

/** \brief Stores the \p j-th column of the Hessian into a sparse vector. */
static void getColH(const qpp_ecrmatrix_t *const H,
                    const qpp_int_t j,
                    const qpp_int_t av[],
                    qpp_int_t indices[],
                    qpp_real_t values[],
                    qpp_int_t *const n);

/** \brief Stores the \p i-th row of the constraint matrix into a sparse vector. */
static void getRowA(const qpp_ecrmatrix_t *const A,
                    const qpp_int_t i,
                    const qpp_int_t av[],
                    qpp_int_t indices[],
                    qpp_real_t values[],
                    qpp_int_t *const n);

/** \brief Stores the \p j-th column of the constraint matrix into a sparse vector. */
static void getColA(const qpp_ecrmatrix_t *const A,
                    const qpp_int_t j,
                    const qpp_int_t ac[],
                    qpp_int_t indices[],
                    qpp_real_t values[],
                    qpp_int_t *const n);


/** \name Preprocessing Methods. */
/**@{*/

/** Checks if the bounds on the variables / constraints are inconsistent. */
static qpp_return_value_t checkBounds(qpp_data_t *const data);

/** Checks if some variables have equal upper and lower bounds and hence can be removed. */
static qpp_return_value_t checkForFixVar(qpp_data_t *const data);

/** Transforms a linear inequality constraint to an equality constraint if the corresponding
 *      dual variable is strictly positive / negative. */
static qpp_return_value_t strongLagrangianBounds(qpp_data_t *const data);

/** Computes upper / lower bounds on the dual variables based on redundant constraints. */
static qpp_return_value_t newLagrangianBounds(qpp_data_t *const data);

/** Tries to eliminate variables that are linearly unconstrained and are not coupled with
 *      other variables. */
static qpp_return_value_t emptyColumns(qpp_data_t *const data);

/** Removes rows that do not contain any nonzero element. */
static qpp_return_value_t emptyRows(qpp_data_t *const data);

/** Removes (implied) free variables that occur in a single linear constraint and that
 *      occur only in the linear term of the objective function. */
static qpp_return_value_t singletonColumns(qpp_data_t *const data);

/** Removes rows that contain exactly one nonzero element while transferring the bounds
 *      implied by such constraints on the corresponding variable. */
static qpp_return_value_t singletonRows(qpp_data_t *const data);

/** Tries to identify redundant and forcing linear constraints (and removes them). Also
 *      tries to identify linear constraints that lead to an infeasible QP. */
static qpp_return_value_t primalConstraints(qpp_data_t *const data);

/** Tries to identify and remove dependent variables (i.e. duplicate columns). */
static qpp_return_value_t duplicateColumns(qpp_data_t *const data);

/** Performs sort of Gaussian elimination in order to remove further nonzeros from the constraint matrix. */
static qpp_return_value_t sparsification(qpp_data_t *const data);

/** Tries to identify (weakly) forcing dual constraints. Also, tighter bounds on the multipliers
 *      corresponding to the linear constraints are computed. */
static qpp_return_value_t dualConstraints(qpp_data_t *const data);

/** \brief Analyzes all (remaining) linear constraints in order to deduce tighter bounds
 *      on the primal variables.
 *
 *  The set of medium bounds will not be altered, but if tighter bounds are computed, the
 *  set of tightest bounds is changed. */
static qpp_return_value_t tightenVarBoundsPC(qpp_data_t *const data);

/** \brief Scales the QP by normalizing the infinity norm of the matrices.
 *
 *  The algorithm is described by Schork "A parametric active set method for
 *  general quadratic programming", Master's Thesis, Heidelberg University, 2015.
 *
 *  \param work Temporary storage of length m+n.
 */
static qpp_return_value_t scaleInf(qpp_data_t *const data,
                                   qpp_int_t max_iter,
                                   qpp_real_t *const work);

/**@}*/


/** Fixes a variable at a certain \p value and removes it from the QP. */
static qpp_return_value_t fixVariable(const qpp_int_t j,
                                      const qpp_real_t value,
                                      qpp_data_t *const data,
                                      const qpp_active_type_t type);

/** \brief Computes implied bounds on the primal variables based on the \p i-th
 *      primal constraint.
 *
 *  The set of medium bounds will not be altered, but if tighter bounds are computed, the
 *  set of tightest bounds is changed.
 */
static qpp_return_value_t varBoundsPCi(const qpp_real_t au_i,
                                       const qpp_real_t al_i,
                                       const qpp_real_t xl[],
                                       const qpp_real_t xu[],
                                       qpp_int_t row_indices[],
                                       qpp_real_t row_values[],
                                       const qpp_int_t row_length,
                                       qpp_int_t lb_indices[],
                                       qpp_real_t lb_values[],
                                       qpp_int_t *lb_length,
                                       qpp_int_t ub_indices[],
                                       qpp_real_t ub_values[],
                                       qpp_int_t *ub_length);


/*  ======================================================================================
    Implementation of interface functions.
    ====================================================================================*/

qpp_return_value_t qppNew(qpp_data_t **pdata,
                          const qpp_int_t m,
                          const qpp_int_t n,
                          const qpp_int_t A_nnz,
                          const qpp_int_t H_nnz,
                          const qpp_int_t stack_size)
{
    qpp_int_t memsize;
    qpp_int_t *iptr;
    qpp_real_t *rptr;
    qpp_data_t *data;
    qpp_options_t *opt;
    qpp_minheap_t *heap;
    qpp_stack_t *stack;
    qpp_ecrmatrix_t *A, *H;
    clock_t clock_time;


    clock_time = clock();

    if (pdata == NULL)
    {
        return QPP_NULL_ARGUMENT;
    }
    *pdata = NULL;

    if ((n <= 0) || (m < 0) || (stack_size <= 0) || (A_nnz < 0) || (H_nnz < 0))
    {
        return QPP_INVALID_ARGUMENT;
    }

    data = (qpp_data_t*) malloc(sizeof(qpp_data_t));
    if (data == NULL)
    {
        return QPP_OUT_OF_MEMORY;
    }

    data->status = QPP_STATUS_NEW;

    /* Compute size of memory that must be allocated */
    memsize = 0;
    memsize += 3*n * (qpp_int_t) sizeof(qpp_int_t);     // av, A_nnz_columns, H_nnz_columns
    memsize += 2*m * (qpp_int_t) sizeof(qpp_int_t);     // ac, A_nnz_rows
    memsize += 5*n * (qpp_int_t) sizeof(qpp_real_t);    // xl_tight, xl_medium, xu_tight, xu_medium, g
    memsize += 2*m * (qpp_int_t) sizeof(qpp_real_t);    // yl, yu
    memsize += 2*m * (qpp_int_t) sizeof(qpp_real_t);    // al, au
    memsize += 3*n * (qpp_int_t) sizeof(qpp_real_t);    // zl, zu, H_diag
    memsize += (m+n) * (qpp_int_t) sizeof(qpp_real_t);  // mem_mpn_r
    memsize += (m+n) * (qpp_int_t) sizeof(qpp_int_t);   // mem_mpn_i

    iptr = (qpp_int_t*) malloc((size_t) memsize);
    opt = (qpp_options_t*) malloc(sizeof(qpp_options_t));
    heap = mhAlloc(m);
    stack = qppStackAlloc(stack_size);
    A = ecrAlloc(m, n, A_nnz);
    H = ecrAlloc(n, n, H_nnz);

    if ( (iptr == NULL) || (opt == NULL) || (heap == NULL) || (stack == NULL) ||
         (A == NULL) || (H == NULL) )
    {
        free(iptr);
        free(opt);
        mhFree(heap);
        qppStackFree(stack);
        ecrFree(A);
        ecrFree(H);
        free(data);
        return QPP_OUT_OF_MEMORY;
    }

    qppSetDefaultOptions(opt);

    data->m = m;
    data->n = n;
    data->found_reduction = 1;
    data->num_iter = 0;
    data->options = opt;
    data->eq_constr_queue = heap;
    data->presolve_stack = stack;
    data->presolve_time = 0;
    data->init_time = 0;
    data->postsolve_time = 0;

    data->A = A;
    data->H = H;

    /* Set memory: */
    data->av = iptr; iptr += n;
    data->ac = iptr; iptr += m;
    data->A_nnz_columns = iptr; iptr += n;
    data->A_nnz_rows = iptr; iptr += m;
    data->H_nnz_columns = iptr; iptr += n;
    data->mem_mpn_i = iptr; iptr += m+n;

    rptr = (void*) iptr;
    data->g = rptr; rptr += n;
    data->xl_medium = rptr; rptr += n;
    data->xl_tight = rptr; rptr += n;
    data->xu_medium = rptr; rptr += n;
    data->xu_tight = rptr; rptr += n;
    data->al = rptr; rptr += m;
    data->au = rptr; rptr += m;
    data->yl = rptr; rptr += m;
    data->yu = rptr; rptr += m;
    data->zl = rptr; rptr += n;
    data->zu = rptr; rptr += n;
    data->H_diag = rptr; rptr += n;
    data->mem_mpn_r = rptr; rptr += m+n;

    data->xl = NULL;
    data->xu = NULL;    /* Will be set later (depending on which bounds are used). */

    /*  Data for presolved QP will be set to zero, no memory will be allocated,
        but dimension will be set to the original dimension */
    data->m_ps = m;
    data->n_ps = n;
    data->A_irn_ps = NULL;
    data->A_jcn_ps = NULL;
    data->A_x_ps = NULL;
    data->A_nnz_ps = 0;
    data->H_irn_ps = NULL;
    data->H_jcn_ps = NULL;
    data->H_x_ps = NULL;
    data->H_nnz_ps = 0;
    data->g_ps = NULL;
    data->f_ps = 0;
    data->xl_ps = NULL;
    data->xu_ps = NULL;
    data->al_ps = NULL;
    data->au_ps = NULL;

    /* Working storage will be allocated later (when using qppInit(...)). */
    data->mem_int = NULL;
    data->mem_for_bounds_int = NULL;
    data->mem_real = NULL;
    data->mem_for_bounds_real = NULL;

    data->status = QPP_STATUS_MEMORY_ALLOCATED;
    data->init_time = clock() - clock_time;

    #ifdef QPP_WRITE_LOGFILE
    data->logfile = NULL;
    #endif

    *pdata = data;

    return QPP_OK;
}


qpp_return_value_t qppInit(qpp_data_t *const data,
                           const qpp_int_t A_irn[],
                           const qpp_int_t A_jcn[],
                           const qpp_real_t A_x[],
                           const qpp_int_t A_nnz,
                           const qpp_int_t H_irn[],
                           const qpp_int_t H_jcn[],
                           const qpp_real_t H_x[],
                           const qpp_int_t H_nnz,
                           const qpp_real_t g[],
                           const qpp_real_t f,
                           const qpp_real_t xl[],
                           const qpp_real_t xu[],
                           const qpp_real_t al[],
                           const qpp_real_t au[])
{
    qpp_return_value_t err;
    qpp_int_t k, l, m, n, i, max_nnz_rc, H_diag_nnz;
    qpp_int_t *iptr;
    qpp_real_t *rptr;
    qpp_options_t *opt;
    qpp_ecrmatrix_t *A, *H;
    clock_t clock_time;

    #ifdef QPP_WRITE_LOGFILE
    time_t _time;
    struct tm *atb;
    qpp_int_t log_level;
    #endif

    if (data == NULL)
    {
        return QPP_NULL_ARGUMENT;
    }

    clock_time = clock();
    opt = data->options;

    #ifdef QPP_WRITE_LOGFILE

    /* Write logfile header. */
    log_level = opt->log_level;
    if (log_level > 0)
    {
        time(&_time);
        atb = gmtime(&_time);
        data->logfile = fopen("qpPresolver.log", "w");
        if (data->logfile != NULL)
        {
            fprintf(data->logfile, "*** qpPresolver logfile ***\n");
            fprintf(data->logfile, "Date: %d/%d/%d\t\tUTC-Time: %2d:%02d:%02d\n\n",
                    atb->tm_year+1900, atb->tm_mon+1, atb->tm_mday, atb->tm_hour,
                    atb->tm_min, atb->tm_sec);
            fprintf(data->logfile, "> qppInit:\n");
            fflush(data->logfile);
        }
        else
        {
            log_level = -1;     /* Writing into logfile not possible. */
        }
    }

    #endif  /* QPP_WRITE_LOGFILE */


    /* CHECKING INPUT: */

    if ( (H_nnz < 0) || (A_nnz < 0) )
    {
        #ifdef QPP_WRITE_LOGFILE
        if (log_level > 0)
        {
            fprintf(data->logfile, "\tError: Invalid argument (number of nonzero "
                    "elements: (A_nnz, H_nnz)\n");
            fflush(data->logfile);
        }
        #endif
        return QPP_INVALID_ARGUMENT;
    }

    /* g must not be empty (null pointer) as #variables > 0 */
    if (g == NULL)
    {
        #ifdef QPP_WRITE_LOGFILE
        if (log_level > 0)
        {
            fprintf(data->logfile, "\tError: NULL argument (gradient vector g)\n");
            fflush(data->logfile);
        }
        #endif
        return QPP_NULL_ARGUMENT;
    }

    if (data->status != QPP_STATUS_MEMORY_ALLOCATED)
    {
        #ifdef QPP_WRITE_LOGFILE
        if (log_level > 0)
        {
            fprintf(data->logfile, "\tError: Wrong internal status. Status = %"
                    QPP_PRID "\n", (qpp_int_t) data->status);
            fflush(data->logfile);
        }
        #endif
        return QPP_STATUS_ERROR;
    }

    m = data->m;
    n = data->n;
    A = data->A;
    H = data->H;

    /* ecrCreateFromCoordinate() reallocates memory if necessary! */
    err = ecrCreateFromCoordinate(A, m, n, A_nnz, A_irn, A_jcn, A_x, QPP_MST_ROW_WISE, QPP_BT_FALSE);
    if (err != QPP_OK)
    {
        #ifdef QPP_WRITE_LOGFILE
        if (log_level > 0)
        {
            fprintf(data->logfile, "\tError: Could not create ecrmatrix (A)\n");
            fflush(data->logfile);
        }
        #endif
        return err;
    }

    /* ecrCreateFromCoordinate() reallocates memory if necessary! */
    err = ecrCreateFromCoordinate(H, n, n, H_nnz, H_irn, H_jcn, H_x, QPP_MST_ROW_WISE, QPP_BT_TRUE);
    if (err != QPP_OK)
    {
        #ifdef QPP_WRITE_LOGFILE
        if (log_level > 0)
        {
            fprintf(data->logfile, "\tError: Could not create ecrmatrix (H)\n");
            fflush(data->logfile);
        }
        #endif
        return err;
    }

    /* INITIALIZE DATA: */

    data->f = f;
    for (i = 0; i < n; ++i)
    {
        data->g[i] = g[i];

        if (xl == NULL)
        {
            data->xl_tight[i] = -QPP_INF;
            data->xl_medium[i] = -QPP_INF;
        }
        else
        {
            data->xl_tight[i] = xl[i];
            data->xl_medium[i] = xl[i];
        }

        if (xu == NULL)
        {
            data->xu_tight[i] = QPP_INF;
            data->xu_medium[i] = QPP_INF;
        }
        else
        {
            data->xu_tight[i] = xu[i];
            data->xu_medium[i] = xu[i];
        }

        data->av[i] = 1;
        data->zl[i] = -QPP_INF;
        data->zu[i] = QPP_INF;
        data->H_diag[i] = 0.0;
        data->H_nnz_columns[i] = 0;
        data->A_nnz_columns[i] = 0;
    }

    /* Get diagonal information of Hessian matrix and make adjustments */
    H_diag_nnz = 0; /* Number of nonzero diagonal elements of the Hessian matrix */
    for (i = 0; i < H->nnz; ++i)
    {
        k = H->irn[i];
        l = H->jcn[i];
        if (k == l)
        {
            ++H_diag_nnz;
            data->H_diag[k] = H->x[i];
            --data->H_nnz_columns[k];

            /*  Adjust links, so that we do not iterate through diagonal entries twice.
                This only works when using ecrmatrix! DO NOT TOUCH THIS! */
            H->jcolst[l] = H->linkcl[H->jcolst[l]];
        }
    }

    for (i = 0; i < m; ++i)
    {
        if (al == NULL)
        {
            data->al[i] = -QPP_INF;
        }
        else
        {
            data->al[i] = al[i];
        }

        if (au == NULL)
        {
            data->au[i] = QPP_INF;
        }
        else
        {
            data->au[i] = au[i];
        }

        data->ac[i] = 1;
        data->yl[i] = -QPP_INF;
        data->yu[i] = QPP_INF;
        data->A_nnz_rows[i] = 0;
    }

    /* Getting number of nonzeros of the rows and columns of A. */
    for (i = 0; i < A->nnz; ++i)
    {
        ++data->A_nnz_columns[A->jcn[i]];
        ++data->A_nnz_rows[A->irn[i]];
    }

    /* Getting number of nonzeros of the rows (= columns) of H. */
    for (i = 0; i < H->nnz; ++i)
    {
        ++data->H_nnz_columns[H->irn[i]];
        ++data->H_nnz_columns[H->jcn[i]];
    }

    /*  Reserving memory for tasks within the presolving methods. In max_nnz_rc
        we will store the maximum number of nonzeros of all rows /columns of A / H. */
    max_nnz_rc = 0;
    for (i = 0; i < n; ++i)
    {
        max_nnz_rc = qppMaxi(max_nnz_rc, data->A_nnz_columns[i]);
        max_nnz_rc = qppMaxi(max_nnz_rc, data->H_nnz_columns[i]);
    }
    for (i = 0; i < m; ++i)
    {
        max_nnz_rc = qppMaxi(max_nnz_rc, data->A_nnz_rows[i]);
    }
    ++max_nnz_rc;

    /* max_nnz_rc = qppMaxi(m, n); */

    iptr = (qpp_int_t*) calloc( (size_t)(6 * max_nnz_rc), (sizeof(qpp_int_t) + sizeof(qpp_real_t)));
    if (iptr == NULL)
    {
        #ifdef QPP_WRITE_LOGFILE
        if (log_level > 0)
        {
            fprintf(data->logfile, "\tError: Not enough memory (working storage)\n");
            fflush(data->logfile);
        }
        #endif
        return QPP_OUT_OF_MEMORY;
    }
    data->mem_int = iptr; iptr += 4*max_nnz_rc;
    data->mem_for_bounds_int = iptr; iptr += 2*max_nnz_rc;
    rptr = (void*) iptr;
    data->mem_real = rptr; rptr += 4*max_nnz_rc;
    data->mem_for_bounds_real = rptr; rptr += 2*max_nnz_rc;
    data->mem_length = max_nnz_rc;

    /* We use the heap as a priority-queue for the "sparsification" process */
    for (i = 0; i < m; ++i)
    {
        if (qppIsEqual(al[i], au[i], opt->eq_tol))
        {
            mhInsert(data->eq_constr_queue, data->A_nnz_rows[i], i);
        }
    }

    /* Writing additional information into the logfile */
    #ifdef QPP_WRITE_LOGFILE
    if (log_level > 0)
    {
        if (log_level > 1)
        {
            fprintf(data->logfile, "\tOriginal QP: m = %" QPP_PRID ", n = %" QPP_PRID
                    ", nnz(A) = %" QPP_PRID ", nnz(H) = %" QPP_PRID "\n", m, n,
                    A_nnz, 2*H_nnz - H_diag_nnz);
        }
        fprintf(data->logfile, "\tok\n\n");
        fflush(data->logfile);
    }
    #endif

    data->status = QPP_STATUS_DATA_SET;
    data->init_time += (clock() - clock_time);

    return QPP_OK;
}


qpp_return_value_t qppPresolve(qpp_data_t *const data)
{
    #ifdef QPP_WRITE_LOGFILE
    qpp_int_t log_level;
    #endif

    qpp_return_value_t err;
    qpp_options_t *opt;
    clock_t clock_time;

    clock_time = clock();
    err = QPP_OK;

    if (data == NULL)
    {
        return QPP_NULL_ARGUMENT;
    }

    opt = data->options;

    #ifdef QPP_WRITE_LOGFILE
    log_level = (data->logfile != NULL) ? opt->log_level : -1;
    if (log_level > 0)
    {
        fprintf(data->logfile, "> qppPresolve:\n");
    }
    #endif

    if (data->status != QPP_STATUS_DATA_SET)
    {
        #ifdef QPP_WRITE_LOGFILE
        if (log_level > 0)
        {
            fprintf(data->logfile, "\tError: Wrong internal status. Status = %"
                    QPP_PRID "\n", (qpp_int_t) data->status);
            fflush(data->logfile);
        }
        #endif
        return QPP_STATUS_ERROR;
    }

    /* Scaling could also be applied several times within the presolving loop. */
    if (opt->enable_scaling)
    {
        if (scaleInf(data, 20, data->mem_mpn_r) != QPP_OK)
        {
            #ifdef QPP_WRITE_LOGFILE
            if (log_level > 0)
            {
                fprintf(data->logfile, "\tWarning: Scaling failed\n");
            }
            #endif
        }
    }

    while ( (data->found_reduction != 0) && (data->num_iter < opt->max_iter) )
    {
        ++data->num_iter;
        data->found_reduction = 0;

        #ifdef QPP_WRITE_LOGFILE
        if (log_level > 1)
        {
            fprintf(data->logfile, "%" QPP_PRID ":", data->num_iter);
        }
        #endif

        err = checkBounds(data);
        if (err != QPP_OK)
        {
            break;
        }


        //err = checkForFixVar(data);
        if (err != QPP_OK)
        {
            break;
        }

        /* This method seems to be of little use. */
        /*err = strongLagrangianBounds(data);
        if (err != QPP_OK)
        {
            break;
        }*/

        err = newLagrangianBounds(data);
        if (err != QPP_OK)
        {
            break;
        }

        err = emptyRows(data);
        if (err != QPP_OK)
        {
            break;
        }

        if (opt->enable_singleton_rows_method)
        {
            err = singletonRows(data);
            if (err != QPP_OK)
            {
                break;
            }
        }

        if (opt->enable_empty_columns_method)
        {
            err = emptyColumns(data);
            if (err != QPP_OK)
            {
                break;
            }
        }

        if (opt->enable_singleton_columns_method)
        {
            err = singletonColumns(data);
            if (err != QPP_OK)
            {
                break;
            }
        }

        if (opt->enable_primal_constraints_method)
        {
            err = primalConstraints(data);
            if (err != QPP_OK)
            {
                break;
            }
        }

        err = emptyRows(data);
        if (err != QPP_OK)
        {
            break;
        }

        if ( (data->num_iter <= 5) && (opt->enable_bound_tightening) )
        {
            err = tightenVarBoundsPC(data);
            if (err != QPP_OK)
            {
                break;
            }
        }

        if (opt->enable_sparsification_method)
        {
            err = sparsification(data);
            if (err != QPP_OK)
            {
                break;
            }
        }

        /* In very rare occurrences: slow. */
        if (opt->enable_duplicate_columns_method)
        {
            err = duplicateColumns(data);
            if (err != QPP_OK)
            {
                break;
            }
        }

        /* The dual constraints method is slow if number of variables is large */
        if (opt->enable_dual_constraints_method)
        {
            err = dualConstraints(data);
            if (err != QPP_OK)
            {
                break;
            }
        }
    }

    /* If scaling is enabled then we scale the data once more. */
    if (opt->enable_scaling)
    {
        if (scaleInf(data, 20, data->mem_mpn_r) != QPP_OK)
        {
            #ifdef QPP_WRITE_LOGFILE
            if (log_level > 0)
            {
                fprintf(data->logfile, "\tWarning: Scaling failed\n");
            }
            #endif
        }
    }

    if (err != QPP_OK)
    {
        #ifdef QPP_WRITE_LOGFILE
        if (log_level > 0)
        {
            fprintf(data->logfile, "\n\tError: Presolving failed. Status = %d\n", err);
            fflush(data->logfile);
        }
        #endif
        return err;
    }

    data->status = QPP_STATUS_PRESOLVED;
    data->presolve_time = clock() - clock_time;

    /*  We immediately compute the dimension of the presolved QP, hence the user knows
        how much memory should be provided when retrieving the presolved QP. */
    err = computeDimensionPresolvedQP(data);
    if (err != QPP_OK)
    {
        #ifdef QPP_WRITE_LOGFILE
        if (log_level > 0)
        {
            fprintf(data->logfile, "\tError: computeDimensionPresolvedQP failed.\n");
            fflush(data->logfile);
        }
        #endif
        return err;
    }

    #ifdef QPP_WRITE_LOGFILE
    if (log_level > 0)
    {
        fprintf(data->logfile, "\tok\n\n");
        fflush(data->logfile);
    }
    #endif

    return QPP_OK;
}


qpp_return_value_t qppGetDimensionPresolvedQP(qpp_data_t *const data,
                                              qpp_int_t *m,
                                              qpp_int_t *n,
                                              qpp_int_t *A_nnz,
                                              qpp_int_t *H_nnz)
{
    #ifdef QPP_WRITE_LOGFILE
    qpp_int_t log_level;
    #endif

    if (data == NULL)
    {
        return QPP_NULL_ARGUMENT;
    }

    #ifdef QPP_WRITE_LOGFILE
    log_level = (data->logfile != NULL) ? data->options->log_level : -1;
    if (log_level > 0)
    {
        fprintf(data->logfile, "> qppGetDimensionPresolvedQP:\n");
    }
    #endif

    if (data->status < QPP_STATUS_PRESOLVED)
    {
        #ifdef QPP_WRITE_LOGFILE
        if (log_level > 0)
        {
            fprintf(data->logfile, "\tError: Wrong internal status. Status = %"
                    QPP_PRID "\n", (qpp_int_t) data->status);
            fflush(data->logfile);
        }
        #endif
        return QPP_STATUS_ERROR;
    }

    /* Note that the dimension of the presolved QP is computed within qppPresolve() */
    if (m != NULL)
    {
        *m = data->m_ps;
    }
    if (n != NULL)
    {
        *n = data->n_ps;
    }
    if (A_nnz != NULL)
    {
        *A_nnz = data->A_nnz_ps;
    }
    if (H_nnz != NULL)
    {
        *H_nnz = data->H_nnz_ps;
    }

    #ifdef QPP_WRITE_LOGFILE
    if (log_level > 0)
    {
        fprintf(data->logfile, "\tok\n\n");
        fflush(data->logfile);
    }
    #endif

    return QPP_OK;
}


qpp_return_value_t qppGetPresolvedQP(qpp_data_t *const data,
                                     qpp_int_t *m,
                                     qpp_int_t *n,
                                     qpp_int_t A_irn[],
                                     qpp_int_t A_jcn[],
                                     qpp_real_t A_x[],
                                     qpp_int_t *A_nnz,
                                     qpp_int_t H_irn[],
                                     qpp_int_t H_jcn[],
                                     qpp_real_t H_x[],
                                     qpp_int_t *H_nnz,
                                     qpp_real_t g[],
                                     qpp_real_t *f,
                                     qpp_real_t xl[],
                                     qpp_real_t xu[],
                                     qpp_real_t al[],
                                     qpp_real_t au[],
                                     const qpp_matrix_sort_type_t sort_type)
{
    #ifdef QPP_WRITE_LOGFILE
    qpp_int_t log_level;
    #endif

    qpp_return_value_t err;

    if (data == NULL)
    {
        return QPP_NULL_ARGUMENT;
    }

    #ifdef QPP_WRITE_LOGFILE
    log_level = (data->logfile != NULL) ? data->options->log_level : -1;
    if (log_level > 0)
    {
        fprintf(data->logfile, "> qppGetPresolvedQP:\n");
    }
    #endif

    /*  If the problem has not been presolved yet, this method will fail (also
        if it has already been postsolved). */
    if ( (data->status < QPP_STATUS_PRESOLVED) ||
         (data->status > QPP_STATUS_PRESOLVED_QP_SET) )
    {
        #ifdef QPP_WRITE_LOGFILE
        if (log_level > 0)
        {
            fprintf(data->logfile, "\tError: Wrong internal status. Status = %"
                    QPP_PRID "\n", (qpp_int_t) data->status);
            fflush(data->logfile);
        }
        #endif
        return QPP_STATUS_ERROR;
    }

    /* Note that the dimension of the presolved QP is known. */
    if ( ((data->m_ps > 0) && ((al == NULL) || (au == NULL))) ||
         ((data->n_ps > 0) && ((g == NULL) || (xl == NULL) || (xu == NULL))) ||
         ((data->A_nnz_ps > 0) && ((A_irn == NULL) || (A_jcn == NULL) || (A_x == NULL))) ||
         ((data->H_nnz_ps > 0) && ((H_irn == NULL) || (H_jcn == NULL) || (H_x == NULL))) ||
         (m == NULL) || (n == NULL) || (A_nnz == NULL) || (H_nnz == NULL) )
    {
        #ifdef QPP_WRITE_LOGFILE
        if (log_level > 0)
        {
            fprintf(data->logfile, "\tError: NULL argument\n");
            fflush(data->logfile);
        }
        #endif
        return QPP_NULL_ARGUMENT;
    }

    /* computePresolvedQP() uses g_ps, xl_ps, ... as storage. */
    data->g_ps = g;
    data->xl_ps = xl;
    data->xu_ps = xu;
    data->al_ps = al;
    data->au_ps = au;
    data->A_x_ps = A_x;
    data->H_x_ps = H_x;
    data->A_irn_ps = A_irn;
    data->A_jcn_ps = A_jcn;
    data->H_irn_ps = H_irn;
    data->H_jcn_ps = H_jcn;

    if (data->status == QPP_STATUS_PRESOLVED)
    {
        err = computePresolvedQP(data, sort_type);
        if (err != QPP_OK)
        {
            return err;
        }
    }

    *m = data->m_ps;
    *n = data->n_ps;
    *A_nnz = data->A_nnz_ps;
    *H_nnz = data->H_nnz_ps;

    if (f != NULL)
    {
        *f = data->f_ps;
    }

    /*  We set data->g_ps to NULL since otherwise the presolver tries to deallocate memory
        that has been allocated by the user, see qppFree(). */
    data->g_ps = NULL;

    #ifdef QPP_WRITE_LOGFILE
    if (log_level > 0)
    {
        if (log_level > 1)
        {
            fprintf(data->logfile, "\tPresolved QP: m = %" QPP_PRID ", n = %" QPP_PRID
                    ", nnz(A) = %" QPP_PRID ", nnz(H) = %" QPP_PRID "\n", data->m_ps,
                    data->n_ps, data->A_nnz_ps, 2*data->H_nnz_ps - data->H_diag_nnz_ps);
        }
        fprintf(data->logfile, "\tok\n\n");
        fflush(data->logfile);
    }
    #endif

    return QPP_OK;
}


qpp_return_value_t qppGetPresolvedQPPtr(qpp_data_t *const data,
                                        qpp_int_t *m,
                                        qpp_int_t *n,
                                        qpp_int_t *A_irn[],
                                        qpp_int_t *A_jcn[],
                                        qpp_real_t *A_x[],
                                        qpp_int_t *A_nnz,
                                        qpp_int_t *H_irn[],
                                        qpp_int_t *H_jcn[],
                                        qpp_real_t *H_x[],
                                        qpp_int_t *H_nnz,
                                        qpp_real_t *g[],
                                        qpp_real_t *f,
                                        qpp_real_t *xl[],
                                        qpp_real_t *xu[],
                                        qpp_real_t *al[],
                                        qpp_real_t *au[],
                                        const qpp_matrix_sort_type_t sort_type)
{
    qpp_int_t memsize, A_nnz_ps, H_nnz_ps, m_ps, n_ps;
    qpp_int_t *iptr;
    qpp_real_t *rptr;
    qpp_return_value_t err;

    #ifdef QPP_WRITE_LOGFILE
    qpp_int_t log_level;
    #endif

    if (data == NULL)
    {
        return QPP_NULL_ARGUMENT;
    }

    #ifdef QPP_WRITE_LOGFILE
    log_level = (data->logfile != NULL) ? data->options->log_level : -1;
    if (log_level > 0)
    {
        fprintf(data->logfile, "> qppGetPresolvedQP:\n");
    }
    #endif

    /*  If the problem has not been presolved yet, this method will fail (also
        if it has already been postsolved). */
    if ( (data->status < QPP_STATUS_PRESOLVED) ||
         (data->status > QPP_STATUS_PRESOLVED_QP_SET) )
    {
        #ifdef QPP_WRITE_LOGFILE
        if (log_level > 0)
        {
            fprintf(data->logfile, "\tError: Wrong internal status. Status = %"
                    QPP_PRID "\n", (qpp_int_t) data->status);
            fflush(data->logfile);
        }
        #endif
        return QPP_STATUS_ERROR;
    }

    /* Note that the dimension of the presolved QP is known. */
    if ( ((data->m_ps > 0) && ((al == NULL) || (au == NULL))) ||
         ((data->n_ps > 0) && ((g == NULL) || (xl == NULL) || (xu == NULL))) ||
         ((data->A_nnz_ps > 0) && ((A_irn == NULL) || (A_jcn == NULL) || (A_x == NULL))) ||
         ((data->H_nnz_ps > 0) && ((H_irn == NULL) || (H_jcn == NULL) || (H_x == NULL))) ||
         (m == NULL) || (n == NULL) || (A_nnz == NULL) || (H_nnz == NULL) )
    {
        #ifdef QPP_WRITE_LOGFILE
        if (log_level > 0)
        {
            fprintf(data->logfile, "\tError: NULL argument\n");
            fflush(data->logfile);
        }
        #endif
        return QPP_NULL_ARGUMENT;
    }

    m_ps = data->m_ps;
    n_ps = data->n_ps;
    A_nnz_ps = data->A_nnz_ps;
    H_nnz_ps = data->H_nnz_ps;

    /*  Presolver has to allocate memory to store the presolved QP. The user may not
        deallocate it. */
    memsize = 0;
    memsize += (2 * (A_nnz_ps + H_nnz_ps)) * (qpp_int_t) sizeof(qpp_int_t); // A_irn_ps, A_jcn_ps, H_irn_ps, H_jcn_ps
    memsize += (A_nnz_ps + H_nnz_ps) * (qpp_int_t) sizeof(qpp_real_t);      // A_x_ps, H_x_ps
    memsize += (3*n_ps + 2*m_ps) * (qpp_int_t) sizeof(qpp_real_t);          // g_ps, xl_ps, xu_ps, al_ps, au_ps

    if (memsize > 0)
    {
        rptr = (qpp_real_t*) malloc((size_t) memsize);
        if (rptr == NULL)
        {
            #ifdef QPP_WRITE_LOGFILE
            if (log_level > 1)
            {
                fprintf(data->logfile, "\tError: Not enough memory for presolved QP\n");
                fflush(data->logfile);
            }
            #endif
            return QPP_OUT_OF_MEMORY;
        }

        data->g_ps = rptr; rptr += n_ps;
        data->xl_ps = rptr; rptr += n_ps;
        data->xu_ps = rptr; rptr += n_ps;
        data->al_ps = rptr; rptr += m_ps;
        data->au_ps = rptr; rptr += m_ps;
        data->A_x_ps = rptr; rptr += A_nnz_ps;
        data->H_x_ps = rptr; rptr += H_nnz_ps;
        iptr = (void*) rptr;
        data->A_irn_ps = iptr; iptr += A_nnz_ps;
        data->A_jcn_ps = iptr; iptr += A_nnz_ps;
        data->H_irn_ps = iptr; iptr += H_nnz_ps;
        data->H_jcn_ps = iptr; iptr += H_nnz_ps;
    }

    if (data->status == QPP_STATUS_PRESOLVED)
    {
        err = computePresolvedQP(data, sort_type);
        if (err != QPP_OK)
        {
            return err;
        }
    }

    *m = m_ps;
    *n = n_ps;
    *A_nnz = A_nnz_ps;
    *H_nnz = H_nnz_ps;

    if (A_nnz_ps > 0)
    {
        *A_irn = data->A_irn_ps;
        *A_jcn = data->A_jcn_ps;
        *A_x = data->A_x_ps;
    }

    if (H_nnz_ps > 0)
    {
        *H_irn = data->H_irn_ps;
        *H_jcn = data->H_jcn_ps;
        *H_x = data->H_x_ps;
    }

    if (n_ps > 0)
    {
        *g = data->g_ps;
        *xl = data->xl_ps;
        *xu = data->xu_ps;
    }

    if (m_ps > 0)
    {
        *al = data->al_ps;
        *au = data->au_ps;
    }

    if (f != NULL)
    {
        *f = data->f_ps;
    }

    #ifdef QPP_WRITE_LOGFILE
    if (log_level > 0)
    {
        if (log_level > 1)
        {
            fprintf(data->logfile, "\tPresolved QP: m = %" QPP_PRID ", n = %" QPP_PRID
                    ", nnz(A) = %" QPP_PRID ", nnz(H) = %" QPP_PRID "\n", data->m_ps,
                    data->n_ps, data->A_nnz_ps, 2*data->H_nnz_ps - data->H_diag_nnz_ps);
        }
        fprintf(data->logfile, "\tok\n\n");
        fflush(data->logfile);
    }
    #endif

    return QPP_OK;
}


qpp_return_value_t qppPostsolve(qpp_data_t *const data,
                                qpp_real_t *const x,
                                qpp_real_t *const y,
                                qpp_real_t *const z,
                                qpp_int_t *const wi,
                                qpp_int_t *const wj,
                                const qpp_int_t use_copies)
{
    qpp_int_t m, n, m_ps, n_ps, i, ii, j, jj, k, kk, it, sign, length1,
        length2, vec_length_1, vec_length_2, stack_ptr;
    qpp_int_t *av, *ac, *indices1, *indices2, *vec_indices_1, *vec_indices_2, *iptr;
    qpp_real_t a_ij, value, alpha, new_y, f, eq_tol;
    qpp_real_t *xl, *xu, *R, *C, *values1, *H_diag, *g, *vec_values_1, *vec_values_2, *rptr;

    #ifdef QPP_WRITE_LOGFILE
    qpp_int_t log_level;
    #endif

    qpp_stack_t *stack;
    qpp_stackelement_t *se;
    qpp_se_fixed_variable_t *fv;
    qpp_se_empty_rows_t *er;
    qpp_se_singleton_row_t *sr;
    qpp_se_forcing_pc_t *fpc;
    qpp_se_redundant_pc_t *rpc;
    qpp_se_new_bounds_pc_t *nbpc;
    qpp_se_singleton_column_t *sc;
    qpp_se_duplicate_column_t *dc;
    qpp_se_sparsification_t *spf;
    qpp_se_scaling_t *scal;
    qpp_tighter_bounds_type_t tbt;
    qpp_ecrmatrix_t *A, *H;
    qpp_options_t *opt;
    clock_t time;

    time = clock();

    if (data == NULL)
    {
        return QPP_NULL_ARGUMENT;
    }
    opt = data->options;

    #ifdef QPP_WRITE_LOGFILE
    log_level = (data->logfile != NULL) ? opt->log_level : -1;
    if (log_level > 0)
    {
        fprintf(data->logfile, "> qppPostsolve:\n");
    }
    #endif

    /*  Check status for consistency regarding internal status.
        "<" is important. Otherwise multiple postsolving is not possible */
    if (data->status < QPP_STATUS_PRESOLVED_QP_SET)
    {
        #ifdef QPP_WRITE_LOGFILE
        if (log_level > 0)
        {
            fprintf(data->logfile, "\tError: Wrong internal status. Status = %"
                    QPP_PRID "\n", (qpp_int_t) data->status);
            fflush(data->logfile);
        }
        #endif
        return QPP_STATUS_ERROR;
    }

    // Check if the bound mode is a valid bound mode
    if ((opt->bound_mode != QPP_BM_MEDIUM_BOUNDS) && (opt->bound_mode != QPP_BM_TIGHTEST_BOUNDS))
    {
        #ifdef QPP_WRITE_LOGFILE
        if (log_level > 0)
        {
            fprintf(data->logfile, "\tError: Unknown Bound Mode\n");
            fflush(data->logfile);
        }
        #endif
        return QPP_INVALID_ARGUMENT;
    }

    m = data->m;
    n = data->n;
    m_ps = data->m_ps;
    n_ps = data->n_ps;

    if ( ((x == NULL) || (z == NULL)) && (n_ps > 0) )
    {
        #ifdef QPP_WRITE_LOGFILE
        if (log_level > 0)
        {
            fprintf(data->logfile, "\tError: NULL argument (x, z)\n");
            fflush(data->logfile);
        }
        #endif
        return QPP_NULL_ARGUMENT;
    }

    if ( (y == NULL) && (m_ps > 0) )
    {
        #ifdef QPP_WRITE_LOGFILE
        if (log_level > 0)
        {
            fprintf(data->logfile, "\tError: NULL argument (py)\n");
            fflush(data->logfile);
        }
        #endif
        return QPP_NULL_ARGUMENT;
    }

    if (use_copies) // We need copies of A, H (and H_diag), g, av and ac
    {
        A = ecrCopy(data->A);
        H = ecrCopy(data->H);   // Copy of H only necessary when using scaling
        if ( (A == NULL) || (H == NULL) )
        {
            #ifdef QPP_WRITE_LOGFILE
            if (log_level > 0)
            {
                fprintf(data->logfile, "\tError: ecrCopy failed\n");
                fflush(data->logfile);
            }
            #endif
            return QPP_OUT_OF_MEMORY;
        }
        iptr = (qpp_int_t*) malloc((size_t)((m+n)) * sizeof(qpp_int_t) +
                                   (size_t)(2*n) * sizeof(qpp_real_t));
        if (iptr == NULL)
        {
            #ifdef QPP_WRITE_LOGFILE
            if (log_level > 0)
            {
                fprintf(data->logfile, "\tError: Not enough memory for copies\n");
                fflush(data->logfile);
            }
            #endif
            return QPP_OUT_OF_MEMORY;
        }
        av = iptr; iptr += n;
        ac = iptr; iptr += m;
        rptr = (void*) iptr;
        g = rptr; rptr += n;
        H_diag = rptr;

        for (i = 0; i < n; ++i)
        {
            av[i] = data->av[i];
            g[i] = data->g[i];
            H_diag[i] = data->H_diag[i];
        }
        for (i = 0; i < m; ++i)
        {
            ac[i] = data->ac[i];
        }
    }
    else
    {
        A = data->A;
        H = data->H;
        g = data->g;
        av = data->av;
        ac = data->ac;
        H_diag = data->H_diag;
    }

    stack = data->presolve_stack;
    f = data->f_ps;
    eq_tol = opt->eq_tol;

    /*A_nnz_rows = data->A_nnz_rows;
    A_nnz_columns = data->A_nnz_columns;
    H_nnz_columns = data->H_nnz_columns;
    xl_tight = data->xl_tight;
    xl_medium = data->xl_medium;
    xu_tight = data->xu_tight;
    xu_medium = data->xu_medium;
    al = data->al;
    au = data->au;*/

    xl = data->xl;
    xu = data->xu;

    vec_indices_1 = data->mem_int;
    vec_indices_2 = &vec_indices_1[data->mem_length];
    vec_values_1  = data->mem_real;
    vec_values_2  = &vec_values_1[data->mem_length];

    /*  The elements x[0],...,x[n_ps-1] are the solution of the presolved problem.
        -> Extract solution! */
    it = n_ps - 1;
    for (j = n-1; j >= 0; --j)
    {
        if (av[j])
        {
            x[j] = x[it];
            z[j] = z[it];
            --it;
        }
        else
        {
            x[j] = 0.0;
            z[j] = 0.0;
        }
    }

    if (wj != NULL)
    {
        it = n_ps - 1;
        for (j = n-1; j >= 0; --j)
        {
            if (av[j])
            {
                wj[j] = wj[it];
                --it;
            }
            else
            {
                wj[j] = QPP_AT_INACTIVE;
            }
        }
    }


    it = m_ps - 1;
    for (i = m-1; i >= 0; --i)
    {
        if (ac[i])
        {
            y[i] = y[it];
            --it;
        }
        else
        {
            y[i] = 0.0;
        }
    }

    if (wi != NULL)
    {
        it = m_ps - 1;
        for (i = m-1; i >= 0; --i)
        {
            if (ac[i])
            {
                wi[i] = wi[it];
                --it;
            }
            else
            {
                wi[i] = QPP_AT_INACTIVE;
            }
        }
    }

    /*  Postsolve loop (we do not delete stack here in order to allow multiple
        postsolving). DO NOT TOUCH THIS! */
    stack_ptr = stack->size - 1;
    while (stack_ptr >= 0)
    {
        se = &stack->elem[stack_ptr--];     /* == Stack pop */
        switch (se->type)
        {
        case QPP_SE_EMPTY_ROW:          /* OK */
            er = se->elements.er;
            i = er->i;

            y[i] = 0.0;
            ac[i] = 1;

            /* Determine optimal working set */
            if (wi != NULL)
            {
                wi[i] = QPP_AT_INACTIVE;    /* Otherwise linearly dependent!!! */
            }

            /*free(er);*/
            break;

        case QPP_SE_SINGLETON_ROW:      /* OK */
            sr = se->elements.sr;
            i = sr->i;
            j = sr->j;
            a_ij = sr->a_ij;

            if (opt->bound_mode == QPP_BM_TIGHTEST_BOUNDS)
            {
                tbt = sr->tbt_t;
            }
            else
            {
                tbt = sr->tbt_m;
            }

            /* If the working sets are given, we use and update them. */
            if ( (wj != NULL) && (wi != NULL) )
            {
                if ( ((tbt == QPP_TBT_BOTH) && (wj[j] != QPP_AT_INACTIVE))  ||
                     ((tbt == QPP_TBT_LOWER) && (wj[j] == QPP_AT_LOWER_BOUND)) ||
                     ((tbt == QPP_TBT_UPPER) && (wj[j] == QPP_AT_UPPER_BOUND)) )
                {
                    y[i] = z[j] / a_ij;

                    /*if (y[i] > 0.0)
                    {
                        wi[i] = QPP_AT_LOWER_BOUND;
                    }
                    else if (y[i] <  0.0)
                    {
                        wi[i] = QPP_AT_UPPER_BOUND;
                    }
                    else
                    {*/
                        /* z_j was zero, but x_j was active -> linear constraint is active! */
                        if (a_ij > 0)
                        {
                            wi[i] = wj[j];
                        }
                        else
                        {
                            wi[i] = -wj[j];
                        }
                    //}

                    z[j] = 0.0;
                    wj[j] = QPP_AT_INACTIVE;
                }
                else
                {
                    y[i] = 0.0;
                    wi[i] = QPP_AT_INACTIVE;
                }
            }
            else    /* No working sets. We only consider the sign of the multiplier z_j. */
            {
                if ( ((tbt == QPP_TBT_BOTH) && (fabs(z[j]) > eq_tol)) ||
                     ((tbt == QPP_TBT_LOWER) && (z[j] > eq_tol)) ||
                     ((tbt == QPP_TBT_UPPER) && (z[j] < -eq_tol)) )
                {
                    y[i] = z[j] / a_ij;
                    z[j] = 0.0;
                }
                else
                {
                    y[i] = 0.0;
                }
            }

            /*xl_tight[j] = sr->xl_tight_old;
            xu_tight[j] = sr->xu_tight_old;
            xl_medium[j] = sr->xl_medium_old;
            xu_medium[j] = sr->xu_medium_old;
            ++A_nnz_columns[j];*/

            ac[i] = 1;

            /*free(sr);*/
            break;

        case QPP_SE_FIXED_VARIABLE:     /* OK */
            fv = se->elements.fv;
            j = fv->j;
            value = fv->value;

            getColH(H, j, av, vec_indices_1, vec_values_1, &vec_length_1);
            for (ii = 0; ii < vec_length_1; ++ii)
            {
                i = vec_indices_1[ii];
                /*++H_nnz_columns[i];*/
                g[i] -= vec_values_1[ii] * value;
            }
            av[j] = 1;
            f -= (g[j] * value + 0.5 * value * value * H_diag[j]);

            getColA(A, j, ac, vec_indices_2, vec_values_2, &vec_length_2);
            /*for (ii = 0; ii < vec_length_2; ++ii)
            {
                i = vec_indices_2[ii];
                a_ij = vec_values_2[ii];
                al[i] += a_ij * value;
                au[i] += a_ij * value;
                ++A_nnz_rows[i];
            }*/
            x[j] = value;
            z[j] = g[j] + qppSparseDot(vec_indices_1, vec_values_1, vec_length_1, x) +
                value*H_diag[j] - qppSparseDot(vec_indices_2, vec_values_2, vec_length_2, y);

            /*  Optimal working set. Problems may occur if lower and upper bounds are
                equal -> sign of multiplier is important! */
            if (wj != NULL)
            {
                if (fv->actv_type != QPP_AT_EQ_CONSTRAINT)
                {
                    wj[j] = fv->actv_type;
                }
                else
                {
                    if (z[j] >= 0.0)
                    {
                        wj[j] = QPP_AT_LOWER_BOUND;
                    }
                    else
                    {
                        wj[j] = QPP_AT_UPPER_BOUND;
                    }
                }
            }

            /*free(fv);*/
            break;

        case QPP_SE_DUPLICATE_COLUMN:   /*  In some cases it is impossible to derive
                                            an optimal working set! (QBRANDY)*/
            dc = se->elements.dc;
            j = dc->j;
            k = dc->k;
            alpha = dc->alpha;

            getColA(A, k, ac, vec_indices_1, vec_values_1, &vec_length_1);

            /*for (ii = 0; ii < vec_length_1; ++ii)
            {
                ++A_nnz_rows[vec_indices_1[ii]];
            }*/

            /*xl_tight[k] = dc->xl_tight_old;
            xu_tight[k] = dc->xu_tight_old;
            xl_medium[k] = dc->xl_medium_old;
            xu_medium[k] = dc->xu_medium_old;*/

            if (opt->bound_mode == QPP_BM_MEDIUM_BOUNDS)
            {
                xl[k] = dc->xl_medium_old;
                xu[k] = dc->xu_medium_old;
            }
            else if (opt->bound_mode == QPP_BM_TIGHTEST_BOUNDS)
            {
                xl[k] = dc->xl_tight_old;
                xu[k] = dc->xu_tight_old;
            }

            av[j] = 1;
            value = x[k];

            if (alpha > 0.0)
            {
                x[k] = qppMaxr(xl[k], value - alpha*xu[j]);
                if (isinf(x[k]))
                {
                    // xl_k and xu_j are infinite!
                    x[k] = qppMinr(xu[k], value - alpha*xl[j]);
                    if (isinf(x[k]))
                    {
                        // x_k and x_j are free!
                        x[k] = value;
                    }
                }
            }
            else
            {
                x[k] = qppMaxr(xl[k], value - alpha*xl[j]);
                if (isinf(x[k]))
                {
                    // xl_k and xl_j are infinite!
                    x[k] = qppMinr(xu[k], value - alpha*xu[j]);
                    if (isinf(x[k]))
                    {
                        // x_k and x_j are free!
                        x[k] = value;
                    }
                }
            }
            x[j] = value/alpha - x[k]/alpha;
            z[j] = alpha * z[k];

            /* Optimal working set. Problems with weakly active variables (cf. QBRANDY)! */
            if (wj != NULL)
            {
                //printf("wj[k] = %" QPP_PRID "\n", wj[k]);
                if (alpha > 0.0)
                {
                    wj[j] = wj[k];
                }
                else
                {
                    wj[j] = -wj[k];
                }
            }

            //printf("k = %" QPP_PRID ", j = %" QPP_PRID "\n", k, j);

            /*free(dc);*/
            break;

        case QPP_SE_REDUNDANT_PC:       /* OK */
            rpc = se->elements.rpc;
            i = rpc->i;

            ac[i] = 1;
            y[i] = 0.0;

            /* Optimal working set. */
            if (wi != NULL)
            {
                wi[i] = QPP_AT_INACTIVE;
            }

            /*getRowA(A, i, av, vec_indices_1, vec_values_1, &vec_length_1);
            for (jj = 0; jj < vec_length_1; ++jj)
            {
                ++A_nnz_columns[vec_indices_1[jj]];
            }*/

            /*free(rpc);*/
            break;

        case QPP_SE_NEWBOUNDS_PC:           /* Multipliers OK, Working Set in some rare cases not optimal (?) */
            nbpc = se->elements.nbpc;
            indices1 = nbpc->lb_indices;
            indices2 = nbpc->ub_indices;
            /*values1 = nbpc->lb_values;
            values2 = nbpc->ub_values;*/
            length1 = nbpc->lb_num;
            length2 = nbpc->ub_num;
            i = nbpc->i;

            getRowA(A, i, av, vec_indices_1, vec_values_1, &vec_length_1);

            /* DETERMINING OPTIMAL MULTIPLIERS */

            /* First we consider the tighter lower bounds. */
            for (kk = 0; kk < length1; ++kk)
            {
                k = indices1[kk];

                if (z[k] > opt->feas_tol)
                {
                    it = 0;

                    while (vec_indices_1[it] != k)
                        ++it;

                    alpha = vec_values_1[it];       /* alpha = a_ik */
                    value = z[k];

                    for (jj = 0; jj < vec_length_1; ++jj)
                    {
                        j = vec_indices_1[jj];
                        z[j] -= vec_values_1[jj] / alpha * value;
                    }

                    z[k] = 0.0;
                    y[i] += value / alpha;
                }
            }

            /* Tighter upper bounds: */
            for (kk = 0; kk < length2; ++kk)
            {
                k = indices2[kk];

                if (z[k] < -opt->feas_tol)
                {
                    it = 0;

                    while (vec_indices_1[it] != k)
                        ++it;

                    alpha = vec_values_1[it];       /* alpha = a_ik */
                    value = z[k];

                    for (jj = 0; jj < vec_length_1; ++jj)
                    {
                        j = vec_indices_1[jj];
                        z[j] -= vec_values_1[jj] / alpha * value;
                    }

                    z[k] = 0.0;
                    y[i] += value / alpha;
                }
            }

            /* OPTIMAL WORKING SET */
            if ( (wi != NULL) && (wj != NULL) )
            {
                /* Tighter lower bounds. */
                for (kk = 0; kk < length1; ++kk)
                {
                    int sup = 0;
                    k = indices1[kk];

                    if (wj[k] == QPP_AT_LOWER_BOUND)
                    {
                        it = 0;

                        while (vec_indices_1[it] != k)
                            ++it;

                        alpha = vec_values_1[it];       /* alpha = a_ik */

                        for (jj = 0; jj < vec_length_1; ++jj)
                        {
                            j = vec_indices_1[jj];
                            if (wj[j] != QPP_AT_INACTIVE)
                                ++sup;

                            if (alpha > 0.0)
                            {
                                if (vec_values_1[jj] > 0.0)
                                {
                                    wj[j] = QPP_AT_UPPER_BOUND;
                                }
                                else
                                {
                                    wj[j] = QPP_AT_LOWER_BOUND;
                                }
                            }
                            else
                            {
                                if (vec_values_1[jj] > 0.0)
                                {
                                    wj[j] = QPP_AT_LOWER_BOUND;
                                }
                                else
                                {
                                    wj[j] = QPP_AT_UPPER_BOUND;
                                }
                            }
                        }
                        wj[k] = QPP_AT_INACTIVE;

                        if (wi[i] != QPP_AT_INACTIVE)
                            ++sup;

                        if (alpha > 0.0)
                        {
                            wi[i] = QPP_AT_LOWER_BOUND;
                        }
                        else
                        {
                            wi[i] = QPP_AT_UPPER_BOUND;
                        }

                        /* If constraint i is equality constraint -> no sign restr. */
                        if (y[i] > opt->feas_tol)
                        {
                            wi[i] = QPP_AT_LOWER_BOUND;
                        }

                        if (y[i] < -opt->feas_tol)
                        {
                            wi[i] = QPP_AT_UPPER_BOUND;
                        }
                        printf("1 Number active components / total: %d / %" QPP_PRID "\n", sup, vec_length_1+1);
                    }
                }

                /* Tighter upper bounds. */
                for (kk = 0; kk < length2; ++kk)
                {
                    int sup = 0;
                    int z0 = 0;
                    k = indices2[kk];

                    if (wj[k] == QPP_AT_UPPER_BOUND)
                    {
                        it = 0;

                        while (vec_indices_1[it] != k)
                            ++it;

                        alpha = vec_values_1[it];       /* alpha = a_ik */
                        for (jj = 0; jj < vec_length_1; ++jj)
                        {
                            j = vec_indices_1[jj];

                            if (wj[j] == QPP_AT_INACTIVE)
                                printf("%" QPP_PRID "  ", j);

                            if (wj[j] != QPP_AT_INACTIVE)
                                ++sup;

                            if (fabs(z[j]) <= opt->feas_tol)
                                ++z0;

                            /*if (wj[j] == QPP_AT_INACTIVE)
                            {
                                getColA(A, j, ac, vec_indices_2, vec_values_2, &vec_length_2);
                                int found = 0;
                                for (ii = 0; ii < vec_length_2; ++ii)
                                {
                                    it = A->irowst[vec_indices_2[ii]];
                                    while (it >= 0)
                                    {
                                        if (av[A->jcn[it]] && (fabs(z[A->jcn[it]]) <= eq_tol) && (wj[A->jcn[it]] != QPP_AT_INACTIVE))
                                        {
                                            wj[A->jcn[it]] = QPP_AT_INACTIVE;
                                            found = 1;
                                            break;
                                        }
                                        it = A->linkrw[it];
                                    }
                                    if (found) break;
                                }
                            }*/

                            if (alpha > 0.0)
                            {
                                if (vec_values_1[jj] > 0.0)
                                {
                                    wj[j] = QPP_AT_LOWER_BOUND;
                                }
                                else
                                {
                                    wj[j] = QPP_AT_UPPER_BOUND;
                                }
                            }
                            else
                            {
                                if (vec_values_1[jj] > 0.0)
                                {
                                    wj[j] = QPP_AT_UPPER_BOUND;
                                }
                                else
                                {
                                    wj[j] = QPP_AT_LOWER_BOUND;
                                }
                            }
                        }
                        wj[k] = QPP_AT_INACTIVE;
                        printf("%" QPP_PRID "  ", k);

                        if (wi[i] != QPP_AT_INACTIVE)
                            ++sup;

                        if (fabs(y[i]) <= opt->feas_tol)
                            ++z0;

                        if (alpha > 0.0)
                        {
                            wi[i] = QPP_AT_UPPER_BOUND;
                        }
                        else
                        {
                            wi[i] = QPP_AT_LOWER_BOUND;
                        }

                        if (y[i] > opt->feas_tol)
                        {
                            wi[i] = QPP_AT_LOWER_BOUND;
                        }

                        if (y[i] < -opt->feas_tol)
                        {
                            wi[i] = QPP_AT_UPPER_BOUND;
                        }

                        /*if (vec_length_1-sup > 0)
                        {
                            for (j = 0; j < n; ++j)
                            {
                                if (av[j] && (fabs(z[j]) <= opt->feas_tol) && (wj[j] != QPP_AT_INACTIVE) )
                                {
                                    printf("\nIndex = %" QPP_PRID "\n", j);
                                    wj[j] = QPP_AT_INACTIVE;
                                    break;
                                }
                            }
                        }*/

                        printf("\n2 Number active components / total / Zero Multipliers: %d / %" QPP_PRID " / %d\n", sup, vec_length_1+1, z0);
                    }
                }
                //wj[378] = QPP_AT_INACTIVE;
                //wj[835] = QPP_AT_INACTIVE;
            }


            /*free(indices1);
            free(indices2);
            free(values1);
            free(values2);
            free(nbpc);*/
            break;

        case QPP_SE_FORCING_PC:         /* OK */
            fpc = se->elements.fpc;
            i = fpc->i;
            indices1 = fpc->indices;    /* = {j | a_ij != 0} */
            values1 = fpc->values;      /* coefficients of row i */
            length1 = fpc->num_neg;     /* #(negative coefficients) */
            length2 = fpc->num_pos;     /* #(positive coefficients) */
            sign = fpc->sign;           /* sign == 1 --> u = al[i]; else l = au[i] */

            /* The multiplier y_i is first considered as zero. Note that all variables
                have been introduced by the fixVariable (postsolving) procedure! */
            new_y = 0.0;

            /*  For simplicity we first consider the linear constraint is active. Later,
                one constraint is dropped from the working set. */
            if (wi != NULL)
            {
                if (sign > 0)
                {
                    wi[i] = QPP_AT_LOWER_BOUND;
                }
                else
                {
                    wi[i] = QPP_AT_UPPER_BOUND;
                }
            }

            /* Check which multipliers z_k have the wrong sign. */
            for (jj = 0; jj < length1; ++jj)
            {
                j = indices1[jj];
                a_ij = values1[jj];
                if ((qpp_real_t)(sign) * (z[j] - a_ij*new_y) < 0.0)
                {
                    new_y = z[j] / a_ij;
                }
            }

            for (jj = length1; jj < length1+length2; ++jj)
            {
                j = indices1[jj];
                a_ij = values1[jj];
                if ((qpp_real_t)(sign) * (z[j] - a_ij*new_y) > 0.0)
                {
                    new_y = z[j] / a_ij;
                }
            }

            ac[i] = 1;
            y[i] = new_y;

            /* We have to update the value of all dual variables z. */
            for (jj = 0; jj < length1+length2; ++jj)
            {
                z[indices1[jj]] -= y[i] * values1[jj];
            }

            /*  Optimal working set. At least one multiplier is zero. The corresponding
                constraint will be considered as inactive. */
            if ( (wi != NULL) && (wj != NULL) )
            {
                if (fabs(y[i]) <= eq_tol)
                {
                    wi[i] = QPP_AT_INACTIVE;
                }
                else
                {
                    /*  Linear constraint must be active, so there is a multiplier z_k
                        with z_k = 0. The corresponding variable will be inactive. */
                    it = -1;    /* For testing / debugging */
                    for (jj = 0; jj < length1+length2; ++jj)
                    {
                        if (fabs(z[indices1[jj]]) <= eq_tol)
                        {
                            wj[indices1[jj]] = QPP_AT_INACTIVE;
                            it = indices1[jj];
                            break;
                        }
                    }
                    assert(it >= 0);
                }
            }

            /*free(indices1);
            free(values1);
            free(fpc);*/
            break;

        case QPP_SE_SPARSIFICATION:     /* Multipliers OK, working set not. */
            spf = se->elements.spf;
            alpha = spf->alpha;
            i = spf->i;
            k = spf->k;

            getRowA(A, i, av, vec_indices_1, vec_values_1, &vec_length_1);

            for (jj = 0; jj < vec_length_1; ++jj)
            {
                j = vec_indices_1[jj];
                it = A->irowst[k];

                /*  The following loop always works, but if the entries are sorted
                    (as they should be) we  do not have to reset it = A->irowst[k]
                    in every iteration */
                while (A->jcn[it] != j)
                {
                    it = A->linkrw[it];
                }

                /*if (A->x[it] == 0.0)
                {
                    A_nnz_rows[k]++;
                    A_nnz_columns[j]++;
                }*/
                A->x[it] += alpha * vec_values_1[jj];
            }

            /*al[k] = spf->al;
            au[k] = spf->au;*/

            y[i] -= alpha * y[k];

            /* Remember: i-th linear constraint is an equality constraint! */
            if ( (wi != NULL) && (wj != NULL) )
            {
                if (fabs(y[k]) <= 1e-8)
                {
                    continue;
                }
                /* Problem if constr. k is active and constraint i is not!!! */
                /*if ( wi[k] && !wi[i] && (fabs(y[k]) > 1e-8) )
                {
                    puts("Problem SPF.");
                    printf("z_j[j] = %.6e\n", z[spf->j]);

                    if (y[i] <= 0)
                    {
                        wi[i] = QPP_AT_UPPER_BOUND;
                    }
                    else //if (y[i] > opt->feas_tol)
                    {
                        wi[i] = QPP_AT_LOWER_BOUND;
                    }
                }*/

                if (wi[i] != QPP_AT_INACTIVE)
                {
                    if (y[i] <= -opt->feas_tol)
                    {
                        wi[i] = QPP_AT_UPPER_BOUND;
                    }
                    else if (y[i] > opt->feas_tol)
                    {
                        wi[i] = QPP_AT_LOWER_BOUND;
                    }
                }

            }


            /*free(spf);*/
            break;

        case QPP_SE_SINGLETON_COLUMN:   /* OK */
            sc = se->elements.sc;
            a_ij = sc->a_ij;
            i = sc->i;
            j = sc->j;

            getRowA(A, i, av, vec_indices_1, vec_values_1, &vec_length_1);

            /* Compute "first part" of x_j: */
            for (jj = 0; jj < vec_length_1; ++jj)
            {
                x[j] -= vec_values_1[jj] * x[vec_indices_1[jj]];
                /*++A_nnz_columns[vec_indices_1[jj]];*/
            }
            x[j] /= a_ij;
            value = g[j] / a_ij;    /* = y_i (from original linear constraint). */

            if (sc->split != 0)
            {
                /*al[i] = sc->al;
                au[i] = sc->au;*/

                x[j] += sc->au / a_ij;  /* sc->au == sc->al */
                av[j] = 1;
                z[j] = -a_ij * y[i];
                /*++A_nnz_rows[i];*/

                if (fabs(value) > QPP_ZERO)
                {
                    f -= value * sc->au;
                    for (jj = 0; jj < vec_length_1; ++jj)
                    {
                        g[vec_indices_1[jj]] += value * vec_values_1[jj];
                    }
                    y[i] += value;
                }

                if ( (wi != NULL) && (wj != NULL) )
                {
                    /* Here, wi[i] stores the status of the i-th (inequality!) constraint */
                    if (wi[i] == QPP_AT_LOWER_BOUND)
                    {
                        if (a_ij > 0.0)
                        {
                            wj[j] = QPP_AT_UPPER_BOUND;
                        }
                        else
                        {
                            wj[j] = QPP_AT_LOWER_BOUND;
                        }
                    }
                    else if (wi[i] == QPP_AT_UPPER_BOUND)
                    {
                        if (a_ij > 0.0)
                        {
                            wj[j] = QPP_AT_LOWER_BOUND;
                        }
                        else
                        {
                            wj[j] = QPP_AT_UPPER_BOUND;
                        }
                    }
                    else    /* QPP_AT_INACTIVE */
                    {
                        /* Inequality constraint is inactive and so is x_j. */
                        wj[j] = QPP_AT_INACTIVE;
                    }
                }

                /*  The original linear (eq.) constraint must be considered as active
                    since then the dimension of the null space of the constraints stays
                    the same. */
                if (wi != NULL)
                {
                    if (y[i] >= 0.0)
                    {
                        wi[i] = QPP_AT_LOWER_BOUND;
                    }
                    else
                    {
                        wi[i] = QPP_AT_UPPER_BOUND;
                    }
                }
            }
            else
            {
                y[i] = value;
                av[j] = 1;
                ac[i] = 1;
                z[j] = 0.0;

                /* x_j was (implied) free without splitting an equality constraint */
                if (fabs(y[i]) > QPP_ZERO)
                {
                    for (jj = 0; jj < vec_length_1; ++jj)
                    {
                        g[vec_indices_1[jj]] += y[i] * vec_values_1[jj];
                    }
                    if (y[i] > 0.0)
                    {
                        x[j] += sc->al / a_ij;
                        f -= sc->al * y[i];
                        if (wi != NULL)
                        {
                            wi[i] = QPP_AT_LOWER_BOUND;
                        }
                    }
                    else
                    {
                        x[j] += sc->au / a_ij;
                        f -= sc->au * y[i];
                        if (wi != NULL)
                        {
                            wi[i] = QPP_AT_UPPER_BOUND;
                        }
                    }
                }
                else
                {
                    /* Constraint shall be active (then, reduced Hessian is pos. def.) */
                    if (isfinite(sc->au))
                    {
                        x[j] += sc->au / a_ij;
                        if (wi != NULL)
                        {
                            wi[i] = QPP_AT_UPPER_BOUND;
                        }
                    }
                    else if (isfinite(sc->al))
                    {
                        x[j] += sc->al / a_ij;
                        if (wi != NULL)
                        {
                            wi[i] = QPP_AT_LOWER_BOUND;
                        }
                    }
                    else
                    {
                        /* Actually, this case should not happen (-> redundant row!) */
                        x[j] = 0.0;
                        if (wi != NULL)
                        {
                            wi[i] = QPP_AT_INACTIVE;
                        }
                    }
                }

                if (wj != NULL)
                {
                    wj[j] = QPP_AT_INACTIVE;
                }
            }

            /*free(indices1);
            free(values1);
            free(sc);*/
            break;

        case QPP_SE_SCALING:        /* OK */
            scal = se->elements.scal;
            R = scal->R;
            C = scal->C;

            for (i = 0; i < H->nnz; ++i)
            {
                if ( av[H->irn[i]] && av[H->jcn[i]] )
                {
                    H->x[i] /= C[H->irn[i]] * C[H->jcn[i]];
                }
            }

            for (i = 0; i < A->nnz; ++i)
            {
                if ( ac[A->irn[i]] && av[A->jcn[i]] )
                {
                    A->x[i] /= R[A->irn[i]] * C[A->jcn[i]];
                }
            }

            for (i = 0; i < m; ++i)
            {
                if (ac[i])
                {
                    y[i] *= R[i];

                    /*al[i] /= R[i];
                    au[i] /= R[i];*/
                }
            }

            for (i = 0; i < n; ++i)
            {
                if (av[i])
                {
                    x[i] *= C[i];
                    z[i] /= C[i];

                    /*xl_tight[i] *= C[i];
                    xu_tight[i] *= C[i];
                    xl_medium[i] *= C[i];
                    xu_medium[i] *= C[i];*/
                    g[i] /= C[i];
                    H_diag[i] /= C[i] * C[i];
                }
            }

            /*free(scal->C);
            free(scal->R);
            free(scal);*/
            break;

        default:
            #ifdef QPP_WRITE_LOGFILE
            if (log_level > 0)
            {
                fprintf(data->logfile, "\tError: Unknown stack element type\n");
                fflush(data->logfile);
            }
            #endif
            return QPP_INVALID_ARGUMENT;
        }
    }

    data->status = QPP_STATUS_POSTSOLVED;

    #ifdef QPP_WRITE_LOGFILE
    if (log_level > 0)
    {
        fprintf(data->logfile, "\tok\n\n");
        fflush(data->logfile);
    }
    #endif

    if (use_copies)
    {
        ecrFree(A);
        ecrFree(H);
        free(av);
    }

    data->postsolve_time = clock() - time;

    return QPP_OK;
}


qpp_options_t qppGetOptions(const qpp_data_t *const data)
{
    qpp_options_t opt;

    if (data == NULL)
    {
        qppSetDefaultOptions(&opt);
    }
    else
    {
        opt = *(data->options);
    }
    return opt;
}


qpp_return_value_t qppSetOptions(qpp_data_t *const data,
                                 const qpp_options_t *const opt)
{
    #ifdef QPP_WRITE_LOGFILE
    qpp_int_t log_level;
    #endif

    qpp_return_value_t err;

    if (data == NULL)
    {
        return QPP_NULL_ARGUMENT;
    }

    err = QPP_OK;

    #ifdef QPP_WRITE_LOGFILE
    log_level = (data->logfile != NULL) ? data->options->log_level : -1;
    if (log_level > 0)
    {
        fprintf(data->logfile, "> qppSetOptions\n");
    }
    #endif

    if (data->status >= QPP_STATUS_PRESOLVED_QP_SET)
    {
        #ifdef QPP_WRITE_LOGFILE
        if (log_level > 0)
        {
            fprintf(data->logfile, "\tError: Options can not be changed once the "
                    "problem has been presolved and the presolved problem has been "
                    "set up\n");
            fflush(data->logfile);
        }
        #endif
        return QPP_INVALID_ARGUMENT;
    }

    if (opt == NULL)
    {
        err = qppSetDefaultOptions(data->options);
        return err;
    }

    qppCopyOptions(data->options, opt);

    err = qppCheckOptions(data->options);
    if (err != QPP_OK)
    {
        #ifdef QPP_WRITE_LOGFILE
        if (log_level > 0)
        {
            fprintf(data->logfile, "\tWarning: qppCheckOptions() recognized \
                    invalid options and replaced them by default values.\n");
            fflush(data->logfile);
        }
        #endif
    }

    #ifdef QPP_WRITE_LOGFILE
    if (log_level > 0)
    {
        fprintf(data->logfile, "\tok\n\n");
        fflush(data->logfile);
    }
    #endif

    return err;
}


void qppFree(qpp_data_t **pdata)
{
    qpp_data_t *data;

    if (pdata == NULL)
    {
        return;
    }

    data = *pdata;
    if (data != NULL)
    {
        free(data->av);
        free(data->mem_int);
        free(data->g_ps);
        qppStackFree(data->presolve_stack);
        mhFree(data->eq_constr_queue);
        ecrFree(data->A);
        ecrFree(data->H);

        #ifdef QPP_WRITE_LOGFILE
        if (data->logfile != NULL)
        {
            if (data->options->log_level > 0)
            {
                fprintf(data->logfile, "> qppFree:\n\tok\n\n");
            }
            fclose(data->logfile);
            data->logfile = NULL;
        }
        #endif

        free(data->options);
        free(data);
        *pdata = NULL;
    }
}


/*  ======================================================================================
    PRESOLVING AND FURTHER STATIC FUNCTIONS
    ====================================================================================*/

static qpp_return_value_t computeDimensionPresolvedQP(qpp_data_t *const data)
{
    qpp_int_t m_ps, n_ps, A_nnz_ps, H_nnz_ps, i, m, n, H_diag_nnz_ps;
    qpp_int_t *ac, *av, *A_nnz_rows, *H_nnz_columns;
    qpp_real_t *H_diag;

    #ifdef QPP_WRITE_LOGFILE
    qpp_int_t log_level;
    char indent[] = "\t";
    #endif

    if (data == NULL)
    {
        return QPP_NULL_ARGUMENT;
    }

    #ifdef QPP_WRITE_LOGFILE
    log_level = (data->logfile != NULL) ? data->options->log_level : -1;
    if (log_level > 2)
    {
        fprintf(data->logfile, "%s>> computeDimensionPresolvedQP:\n", indent);
    }
    #endif

    if (data->status < QPP_STATUS_PRESOLVED)
    {
        #ifdef QPP_WRITE_LOGFILE
        if (log_level > 2)
        {
            fprintf(data->logfile, "%s\tError: Wrong internal status. Status = %"
                    QPP_PRID "\n", indent, (qpp_int_t) data->status);
            fflush(data->logfile);
        }
        #endif
        return QPP_STATUS_ERROR;
    }

    m = data->m;
    n = data->n;
    ac = data->ac;
    av = data->av;
    A_nnz_rows = data->A_nnz_rows;
    H_nnz_columns = data->H_nnz_columns;
    H_diag = data->H_diag;

    m_ps = 0;
    n_ps = 0;
    A_nnz_ps = 0;
    H_nnz_ps = 0;
    H_diag_nnz_ps = 0;

    /* Extract "active" data from matrix A. */
    for (i = 0; i < m; ++i)
    {
        if (ac[i])
        {
            A_nnz_ps += A_nnz_rows[i];
            ++m_ps;
        }
    }

    /* Extract the "active" data from matrix H. */
    for (i = 0; i < n; ++i)
    {
        if (av[i])
        {
            H_nnz_ps += H_nnz_columns[i];
            if (H_diag[i] != 0.0)
            {
                ++H_diag_nnz_ps;
            }
            ++n_ps;
        }
    }

    /*  We only save the lower triangular part of the Hessian matrix, hence the
        number of nonzeros has to be adjusted */
    H_nnz_ps = ((H_nnz_ps - H_diag_nnz_ps) / 2) + H_diag_nnz_ps;

    data->m_ps = m_ps;
    data->n_ps = n_ps;
    data->A_nnz_ps = A_nnz_ps;
    data->H_nnz_ps = H_nnz_ps;
    data->H_diag_nnz_ps = H_diag_nnz_ps;

    #ifdef QPP_WRITE_LOGFILE
    if (log_level > 2)
    {
        fprintf(data->logfile, "%s\tok\n\n", indent);
        fflush(data->logfile);
    }
    #endif

    return QPP_OK;
}

/* OK, but maybe pass whole QP data as parameters instead of keeping them in data struct? */
static qpp_return_value_t computePresolvedQP(qpp_data_t *const data,
                                             const qpp_matrix_sort_type_t sort_type)
{
    qpp_int_t i, it, matrix_it, m_ps, n_ps, m, n, row_index, col_index, matrix_length;
    qpp_int_t *ac, *av, *row_index_offset, *col_index_offset, *matrix_link, *matrix_start;
    qpp_real_t eq_tol;
    qpp_real_t *xl, *xu, *al, *au, *g, *H_diag;
    qpp_ecrmatrix_t *A, *H;

    #ifdef QPP_WRITE_LOGFILE
    qpp_int_t log_level;
    char indent[] = "\t";
    #endif

    if (data == NULL)
    {
        return QPP_NULL_ARGUMENT;
    }

    #ifdef QPP_WRITE_LOGFILE
    log_level = (data->logfile != NULL) ? data->options->log_level : -1;
    if (log_level > 1)
    {
        fprintf(data->logfile, "%s>> computePresolvedQP:\n", indent);
    }
    #endif

    if (data->status != QPP_STATUS_PRESOLVED)
    {
        #ifdef QPP_WRITE_LOGFILE
        if (log_level > 1)
        {
            fprintf(data->logfile, "%s\tError: Wrong internal status. Status = %"
                    QPP_PRID "\n", indent, (qpp_int_t) data->status);
            fflush(data->logfile);
        }
        #endif
        return QPP_STATUS_ERROR;
    }

    if (data->options->bound_mode == QPP_BM_TIGHTEST_BOUNDS)
    {
        data->xl = data->xl_tight;
        data->xu = data->xu_tight;
    }
    else if (data->options->bound_mode == QPP_BM_MEDIUM_BOUNDS)
    {
        data->xl = data->xl_medium;
        data->xu = data->xu_medium;
    }
    else
    {
        #ifdef QPP_WRITE_LOGFILE
        if (log_level > 1)
        {
            fprintf(data->logfile, "%s\tError: Unknown bound mode\n", indent);
            fflush(data->logfile);
        }
        #endif
        return QPP_INVALID_ARGUMENT;
    }

    ac = data->ac;
    av = data->av;
    xl = data->xl;
    xu = data->xu;
    g = data->g;
    al = data->al;
    au = data->au;
    H_diag = data->H_diag;
    A = data->A;
    H = data->H;
    m = data->m;
    n = data->n;
    eq_tol = data->options->eq_tol;

    /*  Note that the dimension of the presolved QP is known due to a previous call of
        computeDimensionPresolvedQP(). */
    m_ps = data->m_ps;
    n_ps = data->n_ps;

    /*  If the presolved problem is empty (i.e. m_ps = n_ps = 0), then we
        can leave this function without extracting data. */
    if ( (m_ps == 0) && (n_ps == 0) )
    {
        data->f_ps = data->f;
        data->status = QPP_STATUS_PRESOLVED_QP_SET;
        return QPP_OK;
    }

    /*  Since some rows and columns of A and H might be deleted, we also have to
        compute new row and column subscripts. For this reason, we use temporary
        arrays: row_index_offset and col_index_offset. The subscripts will
        be shifted by these offsets. */
    qppSetArrayi(data->mem_mpn_i, m+n, 0);
    row_index_offset = data->mem_mpn_i;
    col_index_offset = &row_index_offset[m];

    /* Compute row index offsets */
    if ( (m > 0) && (!ac[0]) )
    {
        row_index_offset[0] = 1;
    }
    for (i = 1; i < m; ++i)
    {
        row_index_offset[i] = row_index_offset[i-1];
        if (!ac[i])
        {
            ++row_index_offset[i];
        }
    }

    /* Compute column index offsets */
    if (!av[0])
    {
        col_index_offset[0] = 1;
    }
    for (i = 1; i < n; ++i)
    {
        col_index_offset[i] = col_index_offset[i-1];
        if (!av[i])
        {
            ++col_index_offset[i];
        }
    }

    data->f_ps = data->f;

    /* Extract active components of al and au */
    it = 0;
    for (i = 0; i < m; ++i)
    {
        if (ac[i])
        {
            data->al_ps[it] = al[i];
            data->au_ps[it] = au[i];
            ++it;
        }
    }

    /* Extract active components of g, xl and xu */
    it = 0;
    for (i = 0; i < n; ++i)
    {
        if (av[i])
        {
            data->g_ps[it] = g[i];
            data->xl_ps[it] = xl[i];
            data->xu_ps[it] = xu[i];
            ++it;
        }
    }

    /*  Adjust row and column subscripts of A and H. Depending on parameter sort_type,
        the matrices will be sorted row- or column-wise (coordinate format!) */
    if (sort_type == QPP_MST_ROW_WISE)
    {
        matrix_length = A->nrow;
        matrix_link = A->linkrw;
        matrix_start = A->irowst;
    }
    else
    {
        matrix_length = A->ncol;
        matrix_link = A->linkcl;
        matrix_start = A->jcolst;
    }

    it = 0;
    for (i = 0; i < matrix_length; ++i)
    {
        matrix_it = matrix_start[i];
        while (matrix_it >= 0)
        {
            row_index = A->irn[matrix_it];
            col_index = A->jcn[matrix_it];
            if ( av[col_index] && ac[row_index] && (A->x[matrix_it] != 0.0) )
            {
                data->A_irn_ps[it] = row_index - row_index_offset[row_index];
                data->A_jcn_ps[it] = col_index - col_index_offset[col_index];
                data->A_x_ps[it] = A->x[matrix_it];
                ++it;
            }
            matrix_it = matrix_link[matrix_it];
        }
    }

    /*  Note that we modified the links of the Hessian matrix since we only store the
        lower triangular part of it! */
    if (sort_type == QPP_MST_ROW_WISE)
    {
        matrix_length = H->nrow;
        matrix_link = H->linkrw;
        matrix_start = H->irowst;
    }
    else
    {
        matrix_length = H->ncol;
        matrix_link = H->linkcl;
        matrix_start = H->jcolst;
    }

    it = 0;
    for (i = 0; i < matrix_length; ++i)
    {
        if (sort_type == QPP_MST_COLUMN_WISE)
        {
            /* Diagonal entries must be treated separately. */
            if ( av[i] && (fabs(H_diag[i]) >= eq_tol) )
            {
                data->H_irn_ps[it] = i - col_index_offset[i];
                data->H_jcn_ps[it] = i - col_index_offset[i];
                data->H_x_ps[it] = H_diag[i];
                ++it;
            }
        }

        matrix_it = matrix_start[i];
        while (matrix_it >= 0)
        {
            row_index = H->irn[matrix_it];
            col_index = H->jcn[matrix_it];
            if ( av[col_index] && av[row_index] )
            {
                data->H_irn_ps[it] = row_index - col_index_offset[row_index];
                data->H_jcn_ps[it] = col_index - col_index_offset[col_index];
                data->H_x_ps[it] = H->x[matrix_it];
                ++it;
            }
            matrix_it = matrix_link[matrix_it];
        }
    }

    /*it = 0;
    for (i = 0; i < A->nnz; ++i)
    {
        row_index = A->irn[i];
        col_index = A->jcn[i];
        if ( av[col_index] && ac[row_index] && (A->x[i] != 0.0) )
        {
            data->A_irn_ps[it] = row_index - row_index_offset[row_index];
            data->A_jcn_ps[it] = col_index - col_index_offset[col_index];
            data->A_x_ps[it] = A->x[i];
            ++it;
        }
    }*/

    /*it = 0;
    for (i = 0; i < H->nnz; ++i)
    {
        row_index = H->irn[i];
        col_index = H->jcn[i];
        if ( av[col_index] && av[row_index] )
        {
            data->H_irn_ps[it] = row_index - col_index_offset[row_index];
            data->H_jcn_ps[it] = col_index - col_index_offset[col_index];
            data->H_x_ps[it] = H->x[i];
            ++it;
        }
    }*/

    data->status = QPP_STATUS_PRESOLVED_QP_SET;

    #ifdef QPP_WRITE_LOGFILE
    if (log_level > 1)
    {
        fprintf(data->logfile, "%s\tok\n\n", indent);
        fflush(data->logfile);
    }
    #endif

    return QPP_OK;
}


static void getColH(const qpp_ecrmatrix_t *const H,
                    const qpp_int_t j,
                    const qpp_int_t av[],
                    qpp_int_t indices[],
                    qpp_real_t values[],
                    qpp_int_t *const n)
{
    qpp_int_t it, counter;
    qpp_int_t *H_irn, *H_jcn, *Hlinkrw, *Hlinkcl;
    qpp_real_t *H_x;

    H_irn = H->irn;
    H_jcn = H->jcn;
    Hlinkrw = H->linkrw;
    Hlinkcl = H->linkcl;
    H_x = H->x;

    counter = 0;

    it = H->irowst[j];
    while (it >= 0)
    {
        if (av[H_jcn[it]])
        {
            indices[counter] = H_jcn[it];
            values[counter++] = H_x[it];
        }
        it = Hlinkrw[it];
    }

    /*  Note that the first entry in a column cannot be a diagonal entry since we
        modified the links in qppInit(). */
    it = H->jcolst[j];
    while (it >= 0)
    {
        if (av[H_irn[it]])
        {
            indices[counter] = H_irn[it];
            values[counter++] = H_x[it];
        }
        it = Hlinkcl[it];
    }
    *n = counter;
}


static void getRowA(const qpp_ecrmatrix_t *const A,
                    const qpp_int_t i,
                    const qpp_int_t av[],
                    qpp_int_t indices[],
                    qpp_real_t values[],
                    qpp_int_t *const n)
{
    qpp_int_t it, counter;
    qpp_int_t *A_jcn, *A_linkrw;
    qpp_real_t *A_x;

    A_jcn = A->jcn;
    A_linkrw = A->linkrw;
    A_x = A->x;

    counter = 0;

    it = A->irowst[i];
    while (it >= 0)
    {
        if ( av[A_jcn[it]] && (A_x[it] != 0.0) )
        {
            indices[counter] = A_jcn[it];
            values[counter++] = A_x[it];
        }
        it = A_linkrw[it];
    }
    *n = counter;
}


static void getColA(const qpp_ecrmatrix_t *const A,
                    const qpp_int_t j,
                    const qpp_int_t ac[],
                    qpp_int_t indices[],
                    qpp_real_t values[],
                    qpp_int_t *const n)
{
    qpp_int_t it, counter;
    qpp_int_t *A_irn, *Alinkcl;
    qpp_real_t *A_x;

    A_irn = A->irn;
    Alinkcl = A->linkcl;
    A_x = A->x;

    counter = 0;

    it = A->jcolst[j];
    while (it >= 0)
    {
        if ( ac[A_irn[it]] && (A_x[it] != 0.0) )
        {
            indices[counter] = A_irn[it];
            values[counter++] = A_x[it];
        }
        it = Alinkcl[it];
    }
    *n = counter;
}


static qpp_return_value_t checkBounds(qpp_data_t *const data)
{
    qpp_int_t i;
    qpp_int_t *ac, *av;
    qpp_real_t feas_tol;
    qpp_real_t *al, *au, *yl, *yu, *xl, *xu, *zl, *zu;

    #ifdef QPP_WRITE_LOGFILE
    qpp_int_t log_level;
    char indent[] = "\t";
    #endif

    ac = data->ac;
    av = data->av;
    feas_tol = data->options->feas_tol;
    al = data->al;
    au = data->au;
    yl = data->yl;
    yu = data->yu;
    xl = data->xl_tight;
    xu = data->xu_tight;
    zl = data->zl;
    zu = data->zu;

    #ifdef QPP_WRITE_LOGFILE
    log_level = (data->logfile != NULL) ? data->options->log_level : -1;
    if (log_level > 1)
    {
        fprintf(data->logfile, "%s>> checkBounds:\n", indent);
    }
    #endif

    for (i = 0; i < data->m; ++i)
    {
        if (ac[i])
        {
            if ( qppIsGreater(al[i], au[i], feas_tol) ||
                 qppIsGreater(yl[i], yu[i], feas_tol) )
            {
                #ifdef QPP_WRITE_LOGFILE
                if (log_level > 1)
                {
                    fprintf(data->logfile, "%s\tError: Bounds incompatible "
                            "(al, au, yl, yu)\n", indent);
                    if (log_level > 2)
                    {
                        fprintf(data->logfile, "%s\t\ti = %" QPP_PRID ", al = %.6e, "
                                "au = %.6e, yl = %.6e, yu = %.6e\n", indent, i, al[i],
                                au[i], yl[i], yu[i]);
                    }
                    fflush(data->logfile);
                }
                #endif
                return QPP_BOUNDS_INCOMPATIBLE;
            }

            if ( (al[i] == QPP_INF) || (au[i] == -QPP_INF) ||
                 (yl[i] == QPP_INF) || (yu[i] == -QPP_INF) )
            {
                #ifdef QPP_WRITE_LOGFILE
                if (log_level > 1)
                {
                    fprintf(data->logfile, "%s\tError: Bounds incompatible "
                            "(al, au, yl, yu)\n", indent);
                    if (log_level > 2)
                    {
                        fprintf(data->logfile, "%s\t\ti = %" QPP_PRID ", al = %.6e, "
                                "au = %.6e, yl = %.6e, yu = %.6e\n", indent, i, al[i],
                                au[i], yl[i], yu[i]);
                    }
                    fflush(data->logfile);
                }
                #endif
                return QPP_BOUNDS_INCOMPATIBLE;
            }
        }
    }

    for (i = 0; i < data->n; ++i)
    {
        if (av[i])
        {
            if ( qppIsGreater(xl[i], xu[i], feas_tol) ||
                 qppIsGreater(zl[i], zu[i], feas_tol) )
            {
                #ifdef QPP_WRITE_LOGFILE
                if (log_level > 1)
                {
                    fprintf(data->logfile, "%s\tError: Bounds incompatible "
                            "(xl, xu, zl, zu)\n", indent);
                    if (log_level > 2)
                    {
                        fprintf(data->logfile, "%s\t\tj = %" QPP_PRID ", xl = %.6e, "
                                "xu = %.6e, zl = %.6e, zu = %.6e\n", indent, i, xl[i],
                                xu[i], zl[i], zu[i]);
                    }
                    fflush(data->logfile);
                }
                #endif
                return QPP_BOUNDS_INCOMPATIBLE;
            }

            if ( (xl[i] == QPP_INF) || (xu[i] == -QPP_INF) ||
                 (zl[i] == QPP_INF) || (zu[i] == -QPP_INF))
            {
                #ifdef QPP_WRITE_LOGFILE
                if (log_level > 1)
                {
                    fprintf(data->logfile, "%s\tError: Bounds incompatible "
                            "(xl, xu, zl, zu)\n", indent);
                    if (log_level > 2)
                    {
                        fprintf(data->logfile, "%s\t\tj = %" QPP_PRID ", xl = %.6e, "
                                "xu = %.6e, zl = %.6e, zu = %.6e\n", indent, i, xl[i],
                                xu[i], zl[i], zu[i]);
                    }
                    fflush(data->logfile);
                }
                #endif
                return QPP_BOUNDS_INCOMPATIBLE;
            }
        }
    }

    #ifdef QPP_WRITE_LOGFILE
    if (log_level > 1)
    {
        fprintf(data->logfile, "%s\tok\n\n", indent);
        fflush(data->logfile);
    }
    #endif

    return QPP_OK;
}


static qpp_return_value_t checkForFixVar(qpp_data_t *const data)
{
    qpp_int_t j, err;
    qpp_int_t *av;
    qpp_real_t eq_tol;
    qpp_real_t *xl, *xu;

    #ifdef QPP_WRITE_LOGFILE
    qpp_int_t log_level;
    char indent[] = "\t";
    #endif

    av = data->av;
    xl = data->xl_tight;
    xu = data->xu_tight;
    eq_tol = data->options->eq_tol;

    #ifdef QPP_WRITE_LOGFILE
    log_level = (data->logfile != NULL) ? data->options->log_level : -1;
    if (log_level > 1)
    {
        fprintf(data->logfile, "%s>> checkForFixVar:\n", indent);
    }
    #endif

    for (j = 0; j < data->n; ++j)
    {
        if ( av[j] && qppIsEqual(xl[j], xu[j], eq_tol) )
        {
            //printf("Fix Var j = %" QPP_PRID "\n", j);

            err = fixVariable(j, (xu[j] + xl[j]) / 2, data, QPP_AT_EQ_CONSTRAINT);
            if (err != QPP_OK)
            {
                #ifdef QPP_WRITE_LOGFILE
                if (log_level > 1)
                {
                    fprintf(data->logfile, "%s\tError: fixVariable failed\n", indent);
                    fflush(data->logfile);
                }
                #endif
                return err;
            }
            data->found_reduction = 1;
        }
    }

    #ifdef QPP_WRITE_LOGFILE
    if (log_level > 1)
    {
        fprintf(data->logfile, "%s\tok\n\n", indent);
        fflush(data->logfile);
    }
    #endif

    return QPP_OK;
}


static qpp_return_value_t fixVariable(const qpp_int_t j,
                                      const qpp_real_t value,
                                      qpp_data_t *const data,
                                      const qpp_active_type_t type)
{
    qpp_return_value_t err;
    qpp_int_t i, ii, heap_index, col_length;
    qpp_int_t *H_nnz_columns, *A_nnz_rows, *col_indices;
    qpp_real_t a_ij, eq_tol;
    qpp_real_t *g, *al, *au, *col_values;
    qpp_se_fixed_variable_t *fv;
    qpp_stack_t *stack;
    qpp_minheap_t *heap;

    #ifdef QPP_WRITE_LOGFILE
    qpp_int_t log_level;
    char indent[] = "\t\t";
    #endif

    H_nnz_columns = data->H_nnz_columns;
    A_nnz_rows = data->A_nnz_rows;
    g = data->g;
    al = data->al;
    au = data->au;
    eq_tol = data->options->eq_tol;
    heap = data->eq_constr_queue;
    stack = data->presolve_stack;
    col_indices = data->mem_mpn_i;
    col_values = data->mem_mpn_r;

    #ifdef QPP_WRITE_LOGFILE
    log_level = (data->logfile != NULL) ? data->options->log_level : -1;
    if (log_level > 3)
    {
        fprintf(data->logfile, "%s>>> fixVariable:\n", indent);
    }
    #endif

    err = qppStackReallocIf(stack, 1);
    if (err != QPP_OK)
    {
        return err;
    }

    fv = (qpp_se_fixed_variable_t*) malloc(sizeof(qpp_se_fixed_variable_t));
    if (fv == NULL)
    {
        return QPP_OUT_OF_MEMORY;
    }

    fv->j = j;
    fv->value = value;
    fv->actv_type = type;
    qppStackPush(stack, fv, QPP_SE_FIXED_VARIABLE);

    getColA(data->A, j, data->ac, col_indices, col_values, &col_length);
    assert(col_length == data->A_nnz_columns[j]);

    for (ii = 0; ii < col_length; ++ii)
    {
        a_ij = col_values[ii];
        i = col_indices[ii];
        al[i] -= a_ij * value;
        au[i] -= a_ij * value;
        --A_nnz_rows[i];

        /*  The number of nonzeros in row i decreased. If row i is an equality
            constraint, then we also have to update the heap of equality constraints */
        if (qppIsEqual(al[i], au[i], eq_tol))
        {
            if (mhSearchForValue(heap, i, &heap_index) == QPP_OK)
            {
                mhDecreaseKey(heap, heap_index, A_nnz_rows[i]);
            }
            else
            {
                mhInsert(heap, A_nnz_rows[i], i);
            }
        }
    }

    data->f += g[j] * value + 0.5 * value * value * data->H_diag[j];

    data->av[j] = 0;
    getColH(data->H, j, data->av, col_indices, col_values, &col_length);
    assert(col_length == H_nnz_columns[j] - (data->H_diag[j] != 0.0 ? 1 : 0));

    for (ii = 0; ii < col_length; ++ii)
    {
        i = col_indices[ii];
        --H_nnz_columns[i];
        g[i] += col_values[ii] * value;
    }

    #ifdef QPP_WRITE_LOGFILE
    if (log_level > 3)
    {
        fprintf(data->logfile, "%s\tok\n\n", indent);
        fflush(data->logfile);
    }
    #endif

    return QPP_OK;
}


static qpp_return_value_t strongLagrangianBounds(qpp_data_t *const data)
{
    qpp_int_t i;
    qpp_int_t *ac, *A_nnz_rows;
    qpp_real_t eq_tol;
    qpp_real_t *al, *au, *yl, *yu;
    qpp_minheap_t *heap;

    #ifdef QPP_WRITE_LOGFILE
    qpp_int_t log_level;
    char indent[] = "\t";
    #endif

    eq_tol = data->options->eq_tol;
    ac = data->ac;
    A_nnz_rows = data->A_nnz_rows;
    al = data->al;
    au = data->au;
    yl = data->yl;
    yu = data->yu;
    heap = data->eq_constr_queue;

    #ifdef QPP_WRITE_LOGFILE
    log_level = (data->logfile != NULL) ? data->options->log_level : -1;
    if (log_level > 1)
    {
        fprintf(data->logfile, "%s>> strongLagrangianBounds:\n", indent);
    }
    #endif

    /* Multipliers on linear constraints: */
    for (i = 0; i < data->m; ++i)
    {
        if (ac[i])
        {
            if ( (yu[i] < 0.0) &&  !qppIsEqual(al[i], au[i], eq_tol) )
            {
                /* y(i) < 0 --> Active at upper bound */
                if (isinf(au[i]))
                {
                    #ifdef QPP_WRITE_LOGFILE
                    if (log_level > 1)
                    {
                        fprintf(data->logfile, "%s\tError: QP is dual infeasible\n",
                                indent);
                        fflush(data->logfile);
                    }
                    #endif
                    return QPP_DUAL_INFEASIBLE;
                }
                /* data->found_reduction = 1; */
                mhInsert(heap, A_nnz_rows[i], i);
                al[i] = au[i];
            }
            else if ( (yl[i] > 0.0) && !qppIsEqual(al[i], au[i], eq_tol) )
            {
                /* y(i) > 0 --> Active at lower bound */
                if (isinf(al[i]))
                {
                    #ifdef QPP_WRITE_LOGFILE
                    if (log_level > 1)
                    {
                        fprintf(data->logfile, "%s\tError: QP is dual infeasible\n",
                                indent);
                        fflush(data->logfile);
                    }
                    #endif
                    return QPP_DUAL_INFEASIBLE;
                }
                /* data->found_reduction = 1; */
                mhInsert(heap, A_nnz_rows[i], i);
                au[i] = al[i];
            }
        }
    }

    #ifdef QPP_WRITE_LOGFILE
    if (log_level > 1)
    {
        fprintf(data->logfile, "%s\tok\n\n", indent);
        fflush(data->logfile);
    }
    #endif

    return QPP_OK;
}


static qpp_return_value_t newLagrangianBounds(qpp_data_t *const data)
{
    qpp_int_t i;
    qpp_int_t *ac, *av;
    qpp_real_t stab_tol;
    qpp_real_t *xl_tight, *xu_tight, *xl_medium, *xu_medium, *al, *au, *yl, *yu, *zl, *zu;

    #ifdef QPP_WRITE_LOGFILE
    qpp_int_t log_level;
    char indent[] = "\t";
    #endif

    ac = data->ac;
    av = data->av;
    xl_tight = data->xl_tight;
    xu_tight = data->xu_tight;
    xl_medium = data->xl_medium;
    xu_medium = data->xu_medium;
    al = data->al;
    au = data->au;
    yl = data->yl;
    yu = data->yu;
    zl = data->zl;
    zu = data->zu;
    stab_tol = data->options->stab_tol;

    #ifdef QPP_WRITE_LOGFILE
    log_level = (data->logfile != NULL) ? data->options->log_level : -1;
    if (log_level > 1)
    {
        fprintf(data->logfile, "%s>> newLagrangianBounds:\n", indent);
    }
    #endif

    // Bounds on Lagrange multipliers for the linear constraints:
    for (i = 0; i < data->m; ++i)
    {
        if (ac[i])
        {
            if (isinf(au[i]))
            {
                yl[i] = qppMaxr(yl[i], 0.0);
            }
            if (isinf(al[i]))
            {
                yu[i] = qppMinr(yu[i], 0.0);
            }
        }
    }

    /*  If we can assure that e.g. x_j can not be active at one of its bound, we can
        derive bounds on the dual variable z_j. */
    for (i = 0; i < data->n; ++i)
    {
        if (av[i])
        {
            if ( isinf(xu_medium[i]) || qppIsGreater(xu_medium[i], xu_tight[i], stab_tol) )
            {
                zl[i] = qppMaxr(zl[i], 0.0);
            }
            if ( isinf(xl_medium[i]) || qppIsGreater(xl_tight[i], xl_medium[i], stab_tol) )
            {
                zu[i] = qppMinr(zu[i], 0.0);
            }
        }
    }

    #ifdef QPP_WRITE_LOGFILE
    if (log_level > 1)
    {
        fprintf(data->logfile, "%s\tok\n\n", indent);
        fflush(data->logfile);
    }
    #endif

    return QPP_OK;
}


static qpp_return_value_t emptyRows(qpp_data_t *const data)
{
    qpp_return_value_t err;
    qpp_int_t i;
    qpp_int_t *ac, *A_nnz_rows;
    qpp_real_t feas_tol;
    qpp_real_t *al, *au;
    qpp_se_empty_rows_t *er;
    qpp_stack_t *stack;

    #ifdef QPP_WRITE_LOGFILE
    qpp_int_t log_level;
    char indent[] = "\t";
    #endif

    ac = data->ac;
    A_nnz_rows = data->A_nnz_rows;
    feas_tol = data->options->feas_tol;
    al = data->al;
    au = data->au;
    stack = data->presolve_stack;

    #ifdef QPP_WRITE_LOGFILE
    log_level = (data->logfile != NULL) ? data->options->log_level : -1;
    if (log_level > 1)
    {
        fprintf(data->logfile, "%s>> emptyRows:\n", indent);
    }
    #endif

    for (i = 0; i < data->m; ++i)
    {
        if ( !ac[i] || (A_nnz_rows[i] != 0) )
        {
            continue;
        }

        if ( (al[i] > feas_tol) || (au[i] < -feas_tol) )
        {
            #ifdef QPP_WRITE_LOGFILE
            if (log_level > 1)
            {
                fprintf(data->logfile, "%s\tError: QP is primal infeasible\n", indent);
                if (log_level > 2)
                {
                    fprintf(data->logfile, "%s\t\ti = %" QPP_PRID ", al = %.6e, "
                            "au = %.6e\n", indent, i, al[i], au[i]);
                }
                fflush(data->logfile);
            }
            #endif
            return QPP_PRIMAL_INFEASIBLE;
        }

        err = qppStackReallocIf(stack, 1);
        if (err != QPP_OK)
        {
            return err;
        }

        er = (qpp_se_empty_rows_t*) malloc(sizeof(qpp_se_empty_rows_t));
        if (er == NULL)
        {
            return QPP_OUT_OF_MEMORY;
        }

        er->i = i;
        ac[i] = 0;
        qppStackPush(stack, er, QPP_SE_EMPTY_ROW);
        data->found_reduction = 1;
    }

    #ifdef QPP_WRITE_LOGFILE
    if (log_level > 1)
    {
        fprintf(data->logfile, "%s\tok\n\n", indent);
        fflush(data->logfile);
    }
    #endif

    return QPP_OK;
}


static qpp_return_value_t singletonRows(qpp_data_t *const data)
{
    qpp_return_value_t err;
    qpp_int_t i, j, it;
    qpp_int_t *ac, *A_nnz_rows, *A_nnz_columns, *av, *A_irowst, *A_jcn, *A_linkrw;
    qpp_real_t a_ij, nlbc, nubc, eq_tol, feas_tol, stab_tol;
    qpp_real_t *A_x, *al, *au, *xl_tight, *xu_tight, *xl_medium, *xu_medium, *zl, *zu;
    qpp_se_singleton_row_t *sr;
    qpp_stack_t *stack;

    #ifdef QPP_WRITE_LOGFILE
    qpp_int_t log_level;
    char indent[] = "\t";
    #endif

    ac = data->ac;
    A_nnz_rows = data->A_nnz_rows;
    A_nnz_columns = data->A_nnz_columns;
    av = data->av;
    A_irowst = data->A->irowst;
    A_jcn = data->A->jcn;
    A_linkrw = data->A->linkrw;
    A_x = data->A->x;
    al = data->al;
    au = data->au;
    xl_tight = data->xl_tight;
    xu_tight = data->xu_tight;
    xl_medium = data->xl_medium;
    xu_medium = data->xu_medium;
    zl = data->zl;
    zu = data->zu;
    eq_tol = data->options->eq_tol;
    feas_tol = data->options->feas_tol;
    stab_tol = data->options->stab_tol;
    stack = data->presolve_stack;

    #ifdef QPP_WRITE_LOGFILE
    log_level = (data->logfile != NULL) ? data->options->log_level : -1;
    if (log_level > 1)
    {
        fprintf(data->logfile, "%s>> singletonRows:\n", indent);
    }
    #endif

    for (i = 0; i < data->m; ++i)
    {
        if ( (A_nnz_rows[i] != 1) || !ac[i] )
        {
            continue;
        }

        /* Get column index j (index of the variable) */
        it = A_irowst[i];
        while ( !av[A_jcn[it]] || (A_x[it] == 0.0) )
        {
            it = A_linkrw[it];
        }
        j = A_jcn[it];
        a_ij = A_x[it];

        if (fabs(a_ij) <= stab_tol*stab_tol)
        {
            continue;
        }

        data->found_reduction = 1;
        --A_nnz_columns[j];
        ac[i] = 0;

        /* Compute bounds implied by the singleton row */
        if (a_ij > 0.0)
        {
            nlbc = al[i] / a_ij;
            nubc = au[i] / a_ij;
        }
        else
        {
            nlbc = au[i] / a_ij;
            nubc = al[i] / a_ij;
        }

        /* Check for infeasibility (contradicting bounds): */
        if (qppIsGreater(nlbc, xu_tight[j], feas_tol) || qppIsGreater(xl_tight[j], nubc, feas_tol))
        {
            #ifdef QPP_WRITE_LOGFILE
            if (log_level > 1)
            {
                fprintf(data->logfile, "%s\tError: QP is primal infeasible\n", indent);
                if (log_level > 2)
                {
                    fprintf(data->logfile, "%s\t\ti = %" QPP_PRID ", j = %" QPP_PRID
                            ", xl_tight = %.6e, xu_tight = %.6e, nlbc = %.6e, nubc = %.6e\n",
                            indent, i, j, xl_tight[j], xu_tight[j], nlbc, nubc);
                }
                fflush(data->logfile);
            }
            #endif
            return QPP_PRIMAL_INFEASIBLE;
        }

        err = qppStackReallocIf(stack, 1);
        if (err != QPP_OK)
        {
            return err;
        }

        sr = (qpp_se_singleton_row_t*) malloc(sizeof(qpp_se_singleton_row_t));
        if (sr == NULL)
        {
            return QPP_OUT_OF_MEMORY;
        }

        sr->a_ij = a_ij;
        sr->tbt_t = QPP_TBT_NONE;
        sr->tbt_m = QPP_TBT_NONE;
        sr->i = i;
        sr->j = j;
        /*sr->xl_tight_old = xl_tight[j];
        sr->xu_tight_old = xu_tight[j];
        sr->xl_medium_old = xl_medium[j];
        sr->xu_medium_old = xu_medium[j];*/

        if (nlbc > xl_medium[j])
        {
            sr->tbt_m |= QPP_TBT_LOWER;
            xl_medium[j] = nlbc;

            /*  We have to reset the bounds on the Lagrange multipliers if the new
                bounds are equal to the tightest bounds since then there is a high
                chance that row i implied these bounds. Since row i will be removed,
                there is no implied bound anymore and hence this bound is no longer
                redundant! */
            if (qppIsEqual(nlbc, xl_tight[j], eq_tol))
            {
                zu[j] = QPP_INF;
            }

            if (nlbc > xl_tight[j])
            {
                zu[j] = QPP_INF;
                sr->tbt_t |= QPP_TBT_LOWER;
                xl_tight[j] = nlbc;
            }
        }

        if (nubc < xu_medium[j])
        {
            sr->tbt_m |= QPP_TBT_UPPER;
            xu_medium[j] = nubc;

            if (qppIsEqual(nubc, xu_tight[j], eq_tol))
            {
                zl[j] = -QPP_INF;
            }

            if (nubc < xu_tight[j])
            {
                zl[j] = -QPP_INF;
                sr->tbt_t |= QPP_TBT_UPPER;
                xu_tight[j] = nubc;
            }
        }

        qppStackPush(stack, sr, QPP_SE_SINGLETON_ROW);

        if (qppIsEqual(xl_tight[j], xu_tight[j], eq_tol))
        {
            err = fixVariable(j, (xu_tight[j] + xl_tight[j]) / 2, data, QPP_AT_EQ_CONSTRAINT);
            if (err != QPP_OK)
            {
                #ifdef QPP_WRITE_LOGFILE
                if (log_level > 1)
                {
                    fprintf(data->logfile, "%s\tError: fixVariable failed\n", indent);
                    fflush(data->logfile);
                }
                #endif
                return err;
            }
        }
    }

    #ifdef QPP_WRITE_LOGFILE
    if (log_level > 1)
    {
        fprintf(data->logfile, "%s\tok\n\n", indent);
        fflush(data->logfile);
    }
    #endif

    return QPP_OK;
}


static qpp_return_value_t emptyColumns(qpp_data_t *const data)
{
    qpp_return_value_t err;
    qpp_int_t j;
    qpp_int_t *av, *A_nnz_columns, *H_nnz_columns;
    qpp_real_t value, h_jj, g_j, xl_j, xu_j;
    qpp_real_t *xl, *xu, *g, *H_diag;
    qpp_active_type_t actv_type;

    #ifdef QPP_WRITE_LOGFILE
    qpp_int_t log_level;
    char indent[] = "\t";
    #endif

    xl = data->xl_tight;
    xu = data->xu_tight;
    g = data->g;
    H_diag = data->H_diag;
    av = data->av;
    A_nnz_columns = data->A_nnz_columns;
    H_nnz_columns = data->H_nnz_columns;
    actv_type = QPP_AT_INACTIVE;

    #ifdef QPP_WRITE_LOGFILE
    log_level = (data->logfile != NULL) ? data->options->log_level : -1;
    if (log_level > 1)
    {
        fprintf(data->logfile, "%s>> emptyColumns:\n", indent);
    }
    #endif

    for (j = 0; j < data->n; ++j)
    {
        if ( !av[j] || (A_nnz_columns[j] != 0) || (H_nnz_columns[j] > 1) )
        {
            continue;
        }

        h_jj = H_diag[j];
        g_j = g[j];
        xl_j = xl[j];
        xu_j = xu[j];

        if (H_nnz_columns[j] == 0)
        {
            /* Solve LP: min g_j * x s.t. xl_j <= x <= xu_j */
            if (g_j > QPP_ZERO)
            {
                if (isinf(xl_j))
                {
                    #ifdef QPP_WRITE_LOGFILE
                    if (log_level > 1)
                    {
                        fprintf(data->logfile, "%s\tError: QP is unbounded\n", indent);
                        if (log_level > 2)
                        {
                            fprintf(data->logfile, "%s\t\tj = %" QPP_PRID "g = %.4e, "
                                    "h = %.4e, xl = %.4e, xu = %.4e\n", indent, j, g_j,
                                    h_jj, xl_j, xu_j);
                        }
                        fflush(data->logfile);
                    }
                    #endif
                    return QPP_UNBOUNDED;
                }
                value = xl_j;
                actv_type = QPP_AT_LOWER_BOUND;
            }
            else if (g_j < -QPP_ZERO)
            {
                if (isinf(xu_j))
                {
                    #ifdef QPP_WRITE_LOGFILE
                    if (log_level > 1)
                    {
                        fprintf(data->logfile, "%s\tError: QP is unbounded\n", indent);
                        if (log_level > 2)
                        {
                            fprintf(data->logfile, "%s\t\tj = %" QPP_PRID ", g = %.4e, "
                                    "h = %.4e, xl = %.4e, xu = %.4e\n", indent, j, g_j,
                                    h_jj, xl_j, xu_j);
                        }
                        fflush(data->logfile);
                    }
                    #endif
                    return QPP_UNBOUNDED;
                }
                value = xu_j;
                actv_type = QPP_AT_UPPER_BOUND;
            }
            else
            {
                /*  g_j = 0 --> x_j can be fixed at arbitrary value, as long as
                    xl_j <= x_j <= xu_j */
                if ( isfinite(xl_j) && isfinite(xu_j) )
                {
                    value = (xu_j + xl_j) / 2;
                    actv_type = QPP_AT_INACTIVE;
                }
                else if (isfinite(xl_j))
                {
                    value = xl_j;
                    actv_type = QPP_AT_LOWER_BOUND;
                }
                else if (isfinite(xu_j))
                {
                    value = xu_j;
                    actv_type = QPP_AT_UPPER_BOUND;
                }
                else
                {
                    value = 0.0;
                    actv_type = QPP_AT_INACTIVE;
                }
            }
            data->found_reduction = 1;
            err = fixVariable(j, value, data, actv_type);
            if (err != QPP_OK)
            {
                #ifdef QPP_WRITE_LOGFILE
                if (log_level > 1)
                {
                    fprintf(data->logfile, "%s\tError: fixVariable failed\n", indent);
                    fflush(data->logfile);
                }
                #endif
                return err;
            }
        }
        else if (h_jj != 0.0)
        {
            /* Solve QP: min 0.5 * h_jj * x^2 + g_j * x s.t. xl_j <= x <= xu_j */
            value = -g_j / h_jj;    /* Solution of unconstrained problem */

            if (h_jj > 0.0)
            {
                /* Convex problem */
                if (value > xu_j)
                {
                    value = xu_j;
                    actv_type = QPP_AT_UPPER_BOUND;
                }
                else if (value < xl_j)
                {
                    value = xl_j;
                    actv_type = QPP_AT_LOWER_BOUND;
                }
            }
            else
            {
                /* Concave problem */
                if ( isinf(xl_j) || isinf(xu_j) )
                {
                    #ifdef QPP_WRITE_LOGFILE
                    if (log_level > 1)
                    {
                        fprintf(data->logfile, "%s\tError: QP is unbounded\n", indent);
                        if (log_level > 2)
                        {
                            fprintf(data->logfile, "%s\t\tj = %" QPP_PRID ", g = %.4e, "
                                    "h = %.4e, xl = %.4e, xu = %.4e\n", indent, j,
                                    g_j, h_jj, xl_j, xu_j);
                        }
                        fflush(data->logfile);
                    }
                    #endif
                    return QPP_UNBOUNDED;
                }

                if (fabs(xu_j-value) >= fabs(value-xl_j))
                {
                    value = xu_j;
                    actv_type = QPP_AT_LOWER_BOUND;
                }
                else
                {
                    value = xl_j;
                    actv_type = QPP_AT_UPPER_BOUND;
                }
            }
            data->found_reduction = 1;
            err = fixVariable(j, value, data, actv_type);
            if (err != QPP_OK)
            {
                #ifdef QPP_WRITE_LOGFILE
                if (log_level > 1)
                {
                    fprintf(data->logfile, "%s\tError: fixVariable failed\n", indent);
                    fflush(data->logfile);
                }
                #endif
                return err;
            }
        }
        /*ec->j = j;
        ec->value = value;
        qppStackPush(data->presolve_stack, ec, QPP_ET_EMPTY_COLUMN);
        data->f += 0.5 * H_diag[j] * value * value + g[j] * value;*/
    }

    #ifdef QPP_WRITE_LOGFILE
    if (log_level > 1)
    {
        fprintf(data->logfile, "%s\tok\n\n", indent);
        fflush(data->logfile);
    }
    #endif

    return QPP_OK;
}


static qpp_return_value_t singletonColumns(qpp_data_t *const data)
{
    qpp_return_value_t err;
    qpp_int_t i, j, jj, it, need_to_split, row_length;
    qpp_int_t *ac, *av, *A_nnz_columns, *A_nnz_rows, *H_nnz_columns, *row_indices,
        *A_irn, *Alinkcl;
    qpp_real_t u, l, a_ij, nlbc, nubc, y, eq_tol, stab_tol;
    qpp_real_t *xl_medium, *xu_medium, *al, *au, *yl, *yu, *zl, *zu, *g, *row_values, *A_x;
    qpp_se_singleton_column_t *sc;
    qpp_stack_t *stack;
    qpp_ecrmatrix_t *A;

    #ifdef QPP_WRITE_LOGFILE
    qpp_int_t log_level;
    char indent[] = "\t";
    #endif

    ac = data->ac;
    av = data->av;
    A_nnz_columns = data->A_nnz_columns;
    A_nnz_rows = data->A_nnz_rows;
    H_nnz_columns = data->H_nnz_columns;
    xl_medium = data->xl_medium;
    xu_medium = data->xu_medium;
    al = data->al;
    au = data->au;
    yl = data->yl;
    yu = data->yu;
    zl = data->zl;
    zu = data->zu;
    g = data->g;
    stack = data->presolve_stack;
    A = data->A;
    A_irn = A->irn;
    Alinkcl = A->linkcl;
    A_x = A->x;
    eq_tol = data->options->eq_tol;
    stab_tol = data->options->stab_tol;

    row_indices = data->mem_int;
    row_values = data->mem_real;

    #ifdef QPP_WRITE_LOGFILE
    log_level = (data->logfile != NULL) ? data->options->log_level : -1;
    if (log_level > 1)
    {
        fprintf(data->logfile, "%s>> singletonColumns:\n", indent);
    }
    #endif

    for (j = 0; j < data->n; ++j)
    {
        if ( !av[j] || (H_nnz_columns[j] != 0) || (A_nnz_columns[j] != 1) )
        {
            continue;
        }

        need_to_split = 1;

        /* Search for row index i */
        it = A->jcolst[j];
        while ( !ac[A_irn[it]] || (A_x[it] == 0.0) )
        {
            it = Alinkcl[it];
        }
        i = A_irn[it];
        a_ij = A_x[it];

        if ( (fabs(a_ij) <= stab_tol*stab_tol) || (A_nnz_rows[i] <= 1) )
        {
            continue;
        }

        y = g[j] / a_ij;
        av[j] = 0;
        getRowA(A, i, av, row_indices, row_values, &row_length);
        assert(row_length == data->A_nnz_rows[i]-1);
        av[j] = 1;

        /*  We prepare a stack entry, although it might be deleted if x_j is not
            (implied) free. */
        err = qppStackReallocIf(stack, 1);
        if (err != QPP_OK)
        {
            return err;
        }

        sc = (qpp_se_singleton_column_t*) malloc(sizeof(qpp_se_singleton_column_t));
        if (sc == NULL)
        {
            return QPP_OUT_OF_MEMORY;
        }

        sc->al = al[i];
        sc->au = au[i];
        sc->a_ij = a_ij;
        sc->i = i;
        sc->j = j;

        u = 0.0;
        l = 0.0;
        nlbc = -QPP_INF;
        nubc = QPP_INF;

        if ( (zl[j] >= 0.0) && (zu[j] <= 0.0) )
        //if ( isinf(xl_medium[j]) && isinf(xu_medium[j]) )
        {
            /*  Choosing nlbc and nubc this way lets us determine the variable to be
                (implied) free! */
            nlbc = QPP_INF;
            nubc = -QPP_INF;
        }
        else
        {
            for (jj = 0; jj < row_length; ++jj)
            {
                /*  It is necessary to use the medium bounds here; otherwise it
                    is possible to declare a variable implied free although it is
                    not implied free! */
                if (row_values[jj] > 0.0)
                {
                    u += row_values[jj] * xu_medium[row_indices[jj]];
                    l += row_values[jj] * xl_medium[row_indices[jj]];
                }
                else
                {
                    u += row_values[jj] * xl_medium[row_indices[jj]];
                    l += row_values[jj] * xu_medium[row_indices[jj]];
                }
            }

            if (a_ij > 0.0)
            {
                nlbc = (al[i] - u) / a_ij;
                nubc = (au[i] - l) / a_ij;
            }
            else
            {
                nlbc = (au[i] - l) / a_ij;
                nubc = (al[i] - u) / a_ij;
            }
        }

        if ( (nlbc >= xl_medium[j]) && (nubc <= xu_medium[j]) )
        {
            need_to_split = 0;
            yu[i] = y;
            yl[i] = y;
        }
        else
        {
            /* x_j is not implied free, but maybe we can derive bounds on y_i */
            //if (nlbc > xl_medium[j])
            if (qppIsGreater(nlbc, xl_medium[j], eq_tol))
            {
                zu[j] = qppMinr(0.0, zu[j]);
                if (a_ij < 0.0)
                {
                    yu[i] = qppMinr(y, yu[i]);
                }
                else
                {
                    yl[i] = qppMaxr(y, yl[i]);
                }
            }
            //else if (nubc < xu_medium[j])
            else if (qppIsGreater(xu_medium[j], nubc, eq_tol))
            {
                zl[j] = qppMaxr(0.0, zl[j]);
                if (a_ij < 0.0)
                {
                    yl[i] = qppMaxr(y, yl[i]);
                }
                else
                {
                    yu[i] = qppMinr(y, yu[i]);
                }
            }
        }

        if (need_to_split == 1)
        {
            /* Use split technique only after 2 iterations of the presolving loop. */
            if ((data->num_iter <= 2) || !qppIsEqual(al[i], au[i], eq_tol))
            {
                free(sc);
            }
            else
            {
                sc->split = 1;
                qppStackPush(stack, sc, QPP_SE_SINGLETON_COLUMN);

                /* We can split row i to make x_j free */
                data->found_reduction = 1;
                --A_nnz_rows[i];
                av[j] = 0;          /* x_j is removed from the QP. */

                /*  We do not update the heap (for equality constraints), because
                    within the "sparsification" process it will be checked, whether
                    constraint i is an equality constraint or not! */
                if (fabs(y) > QPP_ZERO)
                {
                    data->f += y * au[i];   /* au[i] == al[i] */
                    for (jj = 0; jj < row_length; ++jj)
                    {
                        g[row_indices[jj]] -= y * row_values[jj];
                    }
                }

                /* (xu_tight and xl_tight will also work) */
                if (a_ij > 0.0)
                {
                    al[i] -= a_ij * xu_medium[j];
                    au[i] -= a_ij * xl_medium[j];
                }
                else
                {
                    al[i] -= a_ij * xl_medium[j];
                    au[i] -= a_ij * xu_medium[j];
                }

                /*  As we changed i from an equality to an inequality constraint,
                    we treat it as a new constraint, which also means that we have
                    to reset the bounds on the corresponding multiplier. Note that
                    the (primal) variable bounds implied by this new constraint are the same
                    as in the equality constraint! Therefore we do not have to
                    update xl_medium and xu_medium! */
                yl[i] = -QPP_INF;
                yu[i] = QPP_INF;
            }
        }
        else
        {
            sc->split = 0;
            qppStackPush(stack, sc, QPP_SE_SINGLETON_COLUMN);

            data->found_reduction = 1;
            av[j] = 0;
            ac[i] = 0;

            for (jj = 0; jj < row_length; ++jj)
            {
                --A_nnz_columns[row_indices[jj]];
            }

            if (y > QPP_ZERO)
            {
                if (isinf(al[i]))
                {
                    #ifdef QPP_WRITE_LOGFILE
                    if (log_level > 1)
                    {
                        fprintf(data->logfile, "%s\tError: QP is dual infeasible\n",
                                indent);
                        if (log_level > 2)
                        {
                            fprintf(data->logfile, "%s\t\ti = %" QPP_PRID ", j = %"
                                    QPP_PRID ", y = %.6e, al = %.6e, au = %.6e\n",
                                    indent, i, j, y, al[i], au[i]);
                        }
                        fflush(data->logfile);
                    }
                    #endif
                    free(sc);
                    return QPP_DUAL_INFEASIBLE;
                }
                data->f += y * al[i];
                for (jj = 0; jj < row_length; ++jj)
                {
                    g[row_indices[jj]] -= y * row_values[jj];
                }
            }
            else if (y < -QPP_ZERO)
            {
                if (isinf(au[i]))
                {
                    #ifdef QPP_WRITE_LOGFILE
                    if (log_level > 1)
                    {
                        fprintf(data->logfile, "%s\tError: QP is dual infeasible\n",
                                indent);
                        if (log_level > 2)
                        {
                            fprintf(data->logfile, "%s\t\ti = %" QPP_PRID ", j = %"
                                    QPP_PRID ", y = %.6e, al = %.6e, au = %.6e\n",
                                    indent, i, j, y, al[i], au[i]);
                        }
                        fflush(data->logfile);
                    }
                    #endif
                    free(sc);
                    return QPP_DUAL_INFEASIBLE;
                }
                data->f += y * au[i];
                for (jj = 0; jj < row_length; ++jj)
                {
                    g[row_indices[jj]] -= y * row_values[jj];
                }
            }
        }
    }

    #ifdef QPP_WRITE_LOGFILE
    if (log_level > 1)
    {
        fprintf(data->logfile, "%s\tok\n\n", indent);
        fflush(data->logfile);
    }
    #endif

    return QPP_OK;
}


static qpp_return_value_t primalConstraints(qpp_data_t *const data)
{
    qpp_return_value_t err;
    qpp_int_t i, j, jj, it, neg_length, pos_length;
    qpp_int_t *ac, *av, *A_nnz_rows, *A_nnz_columns, *all_indices,
        *neg_indices, *pos_indices, *A_linkrw, *A_jcn;
    qpp_real_t l, u, eq_tol, feas_tol, a_ij;
    qpp_real_t *xl, *xu, *al, *au, *yl, *yu, *all_values, *pos_fix_values,
        *neg_fix_values, *neg_values, *pos_values, *A_x;
    qpp_active_type_t neg_actv_type, pos_actv_type;
    qpp_se_forcing_pc_t *fpc;
    qpp_se_redundant_pc_t *rpc;
    qpp_stack_t *stack;
    qpp_ecrmatrix_t *A;

    #ifdef QPP_WRITE_LOGFILE
    qpp_int_t log_level;
    char indent[] = "\t";
    #endif

    eq_tol = data->options->eq_tol;
    feas_tol = data->options->feas_tol;
    ac = data->ac;
    av = data->av;
    A_nnz_rows = data->A_nnz_rows;
    A_nnz_columns = data->A_nnz_columns;
    xl = data->xl_tight;            /* Tightest bounds are necessary and sufficient */
    xu = data->xu_tight;
    al = data->al;
    au = data->au;
    yl = data->yl;
    yu = data->yu;
    stack = data->presolve_stack;
    A = data->A;
    A_linkrw = A->linkrw;
    A_jcn = A->jcn;
    A_x = A->x;

    neg_indices = data->mem_int;
    pos_indices = &neg_indices[data->mem_length];
    neg_values = data->mem_real;
    pos_values = &neg_values[data->mem_length];
    neg_actv_type = QPP_AT_INACTIVE;
    pos_actv_type = QPP_AT_INACTIVE;

    #ifdef QPP_WRITE_LOGFILE
    log_level = (data->logfile != NULL) ? data->options->log_level : -1;
    if (log_level > 1)
    {
        fprintf(data->logfile, "%s>> primalConstraints:\n", indent);
    }
    #endif

    for (i = 0; i < data->m; ++i)
    {
        if ( !ac[i] || (A_nnz_rows[i] <= 1) )
        {
            continue;
        }

        pos_length = 0;
        neg_length = 0;
        u = 0.0;
        l = 0.0;

        /* First we compute lower and upper implied bounds on the constraint */
        it = A->irowst[i];
        while (it >= 0)
        {
            j = A_jcn[it];
            a_ij = A_x[it];

            if ( av[j] && (a_ij > 0.0) )
            {
                pos_values[pos_length] = a_ij;
                pos_indices[pos_length] = j;
                ++pos_length;
                u += a_ij * xu[j];
                l += a_ij * xl[j];
            }
            else if ( av[j] && (a_ij < 0.0) )
            {
                neg_values[neg_length] = a_ij;
                neg_indices[neg_length] = j;
                ++neg_length;
                u += a_ij * xl[j];
                l += a_ij * xu[j];

            }
            it = A_linkrw[it];
        }

        /* Check if row i leads to an infeasible problem: */
        if ( qppIsGreater(al[i], u, feas_tol) || qppIsGreater(l, au[i], feas_tol) )
        {
            #ifdef QPP_WRITE_LOGFILE
            if (log_level > 1)
            {
                fprintf(data->logfile, "%s\tError: QP is primal infeasible\n", indent);
                if (log_level > 2)
                {
                    fprintf(data->logfile, "%s\t\ti = %" QPP_PRID ", al = %.6e, au = %.6e"
                            ", l = %.6e, u = %.6e\n", indent, i, al[i], au[i], l, u);
                }
                fflush(data->logfile);
            }
            #endif
            return QPP_PRIMAL_INFEASIBLE;
        }

        /* Now check for forcing primal constraint: */
        if ( qppIsEqual(u, al[i], eq_tol) || qppIsEqual(l, au[i], eq_tol) )
        {
            data->found_reduction = 1;

            err = qppStackReallocIf(stack, 1);
            if (err != QPP_OK)
            {
                return err;
            }

            fpc = (qpp_se_forcing_pc_t*) malloc(sizeof(qpp_se_forcing_pc_t));
            all_values = (qpp_real_t*) malloc((size_t)(pos_length + neg_length) * sizeof(qpp_real_t));
            all_indices = (qpp_int_t*) malloc((size_t)(pos_length + neg_length) * sizeof(qpp_int_t));
            if ( (fpc == NULL) || (all_values == NULL) || (all_indices == NULL) )
            {
                free(fpc);
                free(all_values);
                free(all_indices);
                return QPP_OUT_OF_MEMORY;
            }

            /*  The first <neg_length> entries will be the negative ones, followed by
                <pos_length> positive entries. */
            memcpy(all_values, neg_values, (size_t)(neg_length) * sizeof(qpp_real_t));
            memcpy(&all_values[neg_length], pos_values, (size_t)(pos_length) * sizeof(qpp_real_t));
            memcpy(all_indices, neg_indices, (size_t)(neg_length) * sizeof(qpp_int_t));
            memcpy(&all_indices[neg_length], pos_indices, (size_t)(pos_length) * sizeof(qpp_int_t));
            fpc->i = i;
            fpc->indices = all_indices;
            fpc->values = all_values;
            fpc->num_neg = neg_length;
            fpc->num_pos = pos_length;

            /* Fix all variables at their resp. bounds: */
            if (qppIsEqual(u, al[i], eq_tol))
            {
                fpc->sign = 1;
                pos_fix_values = xu;
                neg_fix_values = xl;
                pos_actv_type = QPP_AT_UPPER_BOUND;
                neg_actv_type = QPP_AT_LOWER_BOUND;
            }
            else
            {
                fpc->sign = -1;
                pos_fix_values = xl;
                neg_fix_values = xu;
                pos_actv_type = QPP_AT_LOWER_BOUND;
                neg_actv_type = QPP_AT_UPPER_BOUND;
            }

            qppStackPush(stack, fpc, QPP_SE_FORCING_PC);

            for (jj = 0; jj < neg_length; ++jj)
            {
                j = neg_indices[jj];
                err = fixVariable(j, neg_fix_values[j], data, neg_actv_type);
                if (err != QPP_OK)
                {
                    #ifdef QPP_WRITE_LOGFILE
                    if (log_level > 1)
                    {
                        fprintf(data->logfile, "%s\tError: fixVariable failed\n", indent);
                        fflush(data->logfile);
                    }
                    #endif
                    return err;
                }
            }

            for (jj = 0; jj < pos_length; ++jj)
            {
                j = pos_indices[jj];
                err = fixVariable(j, pos_fix_values[j], data, pos_actv_type);
                if (err != QPP_OK)
                {
                    #ifdef QPP_WRITE_LOGFILE
                    if (log_level > 1)
                    {
                        fprintf(data->logfile, "%s\tError: fixVariable failed\n", indent);
                        fflush(data->logfile);
                    }
                    #endif
                    return err;
                }
            }
            ac[i] = 0;
        }
        else
        {
            /*  Redundant primal constraint. Note that "<=" and ">=" is okay (rather
                than "<" and ">"). */
            if ( (al[i] <= l)  && (au[i] >= u) )
            {
                data->found_reduction = 1;

                err = qppStackReallocIf(stack, 1);
                if (err != QPP_OK)
                {
                    return err;
                }

                rpc = (qpp_se_redundant_pc_t*) malloc (sizeof(qpp_se_redundant_pc_t));
                if (rpc == NULL)
                {
                    return QPP_OUT_OF_MEMORY;
                }

                rpc->i = i;
                ac[i] = 0;

                for (jj = 0; jj < neg_length; ++jj)
                {
                    --A_nnz_columns[neg_indices[jj]];
                }
                for (jj = 0; jj < pos_length; ++jj)
                {
                    --A_nnz_columns[pos_indices[jj]];
                }
                qppStackPush(stack, rpc, QPP_SE_REDUNDANT_PC);
            }
            /* Semi-redundant primal constraint: */
            else if (al[i] <= l)
            {
                yu[i] = qppMinr(0.0, yu[i]);
            }
            else if (au[i] >= u)
            {
                yl[i] = qppMaxr(0.0, yl[i]);
            }
        }
    }

    #ifdef QPP_WRITE_LOGFILE
    if (log_level > 1)
    {
        fprintf(data->logfile, "%s\tok\n\n", indent);
        fflush(data->logfile);
    }
    #endif

    return QPP_OK;
}


static qpp_return_value_t duplicateColumns(qpp_data_t *const data)
{
    qpp_return_value_t err;
    qpp_int_t i, ii, j, jj, k, it_k, it_j, heap_index, row_min_count, row_index,
        is_duplicate_column, row_length, col_length_k, col_length_j;
    qpp_int_t *row_indices, *col_indices_k, *col_indices_j, *ac, *av, *H_nnz_columns,
        *A_nnz_columns, *A_nnz_rows, *link_k, *link_j, *ind_k, *ind_j;
    qpp_real_t alpha, coeff_ik, coeff_ij, eq_tol;
    qpp_real_t *row_values, *col_values_k, *col_values_j, *H_diag, *g, *xl_tight, *xu_tight,
        *xl_medium, *xu_medium, *zl, *zu, *al, *au, *A_x, *H_x;
    qpp_stack_t *stack;
    qpp_se_duplicate_column_t *dc;
    qpp_ecrmatrix_t *A, *H;
    qpp_minheap_t *heap;

    #ifdef QPP_WRITE_LOGFILE
    qpp_int_t log_level;
    char indent[] = "\t";
    #endif

    eq_tol = data->options->eq_tol;
    ac = data->ac;
    av = data->av;
    H_nnz_columns = data->H_nnz_columns;
    A_nnz_columns = data->A_nnz_columns;
    A_nnz_rows = data->A_nnz_rows;
    H_diag = data->H_diag;
    g = data->g;
    xl_tight = data->xl_tight;
    xu_tight = data->xu_tight;
    xl_medium = data->xl_medium;
    xu_medium = data->xu_medium;
    zl = data->zl;
    zu = data->zu;
    al = data->al;
    au = data->au;
    stack = data->presolve_stack;
    A = data->A;
    A_x = A->x;
    H = data->H;
    H_x = H->x;
    heap = data->eq_constr_queue;

    k = data->mem_length;
    row_indices = data->mem_int;
    col_indices_j = &row_indices[k];
    col_indices_k = &col_indices_j[k];
    row_values = data->mem_real;
    col_values_j = &row_values[k];
    col_values_k = &col_values_j[k];

    #ifdef QPP_WRITE_LOGFILE
    log_level = (data->logfile != NULL) ? data->options->log_level : -1;
    if (log_level > 1)
    {
        fprintf(data->logfile, "%s>> duplicateColumns:\n", indent);
        fflush(data->logfile);
    }
    #endif

    for (k = 0; k < data->n; ++k)
    {
        if ( !av[k] || ((H_nnz_columns[k] == 1) && (H_diag[k] != 0.0)) ||
            ((A_nnz_columns[k] == 0) && (H_nnz_columns[k] == 0)) )
        {
            continue;
        }

        /* We have to find the row with the least number of entries: */
        row_min_count = QPP_MAX_INT;
        row_index = -1;
        coeff_ik = 0.0;

        if ( ((H_nnz_columns[k] > 0) && (H_nnz_columns[k] < A_nnz_columns[k])) || (A_nnz_columns[k] == 0) )
        {
            getColH(H, k, av, col_indices_k, col_values_k, &col_length_k);
            assert(col_length_k == H_nnz_columns[k]);

            for (ii = 0; ii < col_length_k; ++ii)
            {
                i = col_indices_k[ii];
                if (H_nnz_columns[i] < row_min_count)
                {
                    row_min_count = H_nnz_columns[i];
                    row_index = i;
                    coeff_ik = col_values_k[ii];
                }
            }
            assert(row_index != -1);
            getColH(H, row_index, av, row_indices, row_values, &row_length);
            assert(row_length == H_nnz_columns[row_index]);
        }
        else
        {
            getColA(A, k, ac, col_indices_k, col_values_k, &col_length_k);
            assert(col_length_k == A_nnz_columns[k]);
            for (ii = 0; ii < col_length_k; ++ii)
            {
                i = col_indices_k[ii];
                if (A_nnz_rows[i] <= row_min_count)
                {
                    row_min_count = A_nnz_rows[i];
                    row_index = i;
                    coeff_ik = col_values_k[ii];
                }
            }
            assert(row_index != -1);
            getRowA(A, row_index, av, row_indices, row_values, &row_length);
            assert(row_length == A_nnz_rows[row_index]);
        }

        is_duplicate_column = 1;
        for (jj = 0; jj < row_length; ++jj)
        {
            j = row_indices[jj];
            if ( (A_nnz_columns[k] != A_nnz_columns[j]) || (j <= k) ||
                 (H_nnz_columns[k] != H_nnz_columns[j]) )
            {
                continue;
            }

            coeff_ij = row_values[jj];
            alpha = coeff_ij / coeff_ik;

            /*  At the first sight that column j and k are not duplicate columns
                we will continue with the next candidate. First we check the Hessian
                matrix: */
            link_k = H->linkrw;
            link_j = H->linkrw;
            ind_k = H->jcn;
            ind_j = H->jcn;
            it_k = H->irowst[k];
            it_j = H->irowst[j];

            while (1)
            {
                /*  Iterate through columns k and j simultaneously and compare
                    both entries (check if the values and subscripts are correct) */
                while ( (it_k >= 0) && !av[ind_k[it_k]] )
                {
                    it_k = link_k[it_k];
                }

                if (it_k < 0)
                {
                    if (link_k == H->linkrw)
                    {
                        it_k = H->jcolst[k];     /* This can not be a diagonal entry (see qppInit) */
                        link_k = H->linkcl;
                        ind_k = H->irn;
                        while ( (it_k >= 0) && !av[ind_k[it_k]] )
                        {
                            it_k = link_k[it_k];
                        }
                        if (it_k < 0) break;
                    }
                    else
                    {
                        break;
                    }
                }

                while ( (it_j >= 0) && !av[ind_j[it_j]] )
                {
                    it_j = link_j[it_j];
                }

                if (it_j < 0)
                {
                    if (link_j == H->linkrw)
                    {
                        it_j = H->jcolst[j];
                        link_j = H->linkcl;
                        ind_j = H->irn;
                        while ( (it_j >= 0) && !av[ind_j[it_j]] )
                        {
                            it_j = link_j[it_j];
                        }
                        if (it_j < 0) break;
                    }
                    else
                    {
                        break;
                    }
                }

                if ( (ind_k[it_k] != ind_j[it_j]) ||
                      !qppIsEqual(alpha*H_x[it_k], H_x[it_j], eq_tol) )
                {
                    is_duplicate_column = 0;
                    break;
                }
                it_k = link_k[it_k];
                it_j = link_j[it_j];
            }

            if (!is_duplicate_column)
            {
                continue;
            }

            /*  Check if columns j and k of constraint matrix A are duplicate.
                We extract column k (since we might need it later anyways) and
                check step by step, if the columns are duplicate columns: */
            getColA(A, k, ac, col_indices_k, col_values_k, &col_length_k);
            assert(col_length_k == A_nnz_columns[k]);
            it_j = A->jcolst[j];
            ind_j = A->irn;
            link_j = A->linkcl;

            for (ii = 0; ii < col_length_k; ++ii)
            {
                /*  it_j can not be < 0, as the number of nonzeros in both
                    columns are equal!!! */
                while ( !ac[ind_j[it_j]] || (A_x[it_j] == 0.0) )
                {
                    it_j = link_j[it_j];
                }

                if ( (col_indices_k[ii] != ind_j[it_j]) ||
                      !qppIsEqual(alpha*col_values_k[ii], A_x[it_j], eq_tol) )
                {
                    is_duplicate_column = 0;
                    break;
                }
                it_j = link_j[it_j];
            }

            if (!is_duplicate_column)
            {
                continue;
            }

            /* If we reached this point then we have a duplicate column! */
            if (qppIsEqual(alpha*g[k], g[j], QPP_ZERO))
            {
                /*  Within qpOASES: Do not use this part since it is in some cases
                    not possible to compute optimal working sets. */
                //continue;

                err = qppStackReallocIf(stack, 1);
                if (err != QPP_OK)
                {
                    return err;
                }

                dc = (qpp_se_duplicate_column_t*) malloc(sizeof(qpp_se_duplicate_column_t));
                if (dc == NULL)
                {
                    return QPP_OUT_OF_MEMORY;
                }

                data->found_reduction = 1;
                av[j] = 0;

                for (ii = 0; ii < col_length_k; ++ii)
                {
                    i = col_indices_k[ii];
                    A_nnz_rows[i]--;
                    /* Update heap since the number of nonzeros in row i decreased */
                    if (qppIsEqual(al[i], au[i], eq_tol))
                    {
                        if (mhSearchForValue(heap, i, &heap_index) == QPP_OK)
                        {
                            mhDecreaseKey(heap, heap_index, A_nnz_rows[i]);
                        }
                        else
                        {
                            mhInsert(heap, A_nnz_rows[i], i);
                        }
                    }
                }

                /* Also the number of nonzeros in the Hessian decreases! */
                getColH(H, j, av, col_indices_j, col_values_j, &col_length_j);
                for (ii = 0; ii < col_length_j; ++ii)
                {
                    --H_nnz_columns[col_indices_j[ii]];
                }

                dc->alpha = alpha;
                dc->j = j;
                dc->k = k;
                dc->xl_tight_old = xl_tight[k];
                dc->xu_tight_old = xu_tight[k];
                dc->xl_medium_old = xl_medium[k];
                dc->xu_medium_old = xu_medium[k];

                if (alpha > 0.0)
                {
                    xu_tight[k] += alpha * xu_tight[j];
                    xl_tight[k] += alpha * xl_tight[j];
                    xu_medium[k] += alpha * xu_medium[j];
                    xl_medium[k] += alpha * xl_medium[j];
                }
                else
                {
                    xu_tight[k] += alpha * xl_tight[j];
                    xl_tight[k] += alpha * xu_tight[j];
                    xu_medium[k] += alpha * xl_medium[j];
                    xl_medium[k] += alpha * xu_medium[j];
                }
                qppStackPush(stack, dc, QPP_SE_DUPLICATE_COLUMN);
            }
            else if (qppIsGreater(alpha*g[k], g[j], QPP_ZERO))   /* "else" is important! */
            {
                if ( ((zl[k] >= 0.0) && (alpha < 0.0)) ||
                     ((zu[k] <= 0.0) && (alpha > 0.0)) )
                {
                    data->found_reduction = 1;
                    if (isinf(xu_medium[j]))
                    {
                        #ifdef QPP_WRITE_LOGFILE
                        if (log_level > 1)
                        {
                            fprintf(data->logfile, "%s\tError: QP is dual infeasible\n",
                                    indent);
                            if (log_level > 2)
                            {
                                fprintf(data->logfile, "%s\t\tk = %" QPP_PRID ", j = %"
                                        QPP_PRID ", zl = %.6e, zu = %.6e, alpha = "
                                        "%.6e, >\n", indent, k, j, zl[k], zu[k], alpha);
                            }
                            fflush(data->logfile);
                        }
                        #endif
                        return QPP_DUAL_INFEASIBLE;
                    }
                    else
                    {
                        err = fixVariable(j, xu_medium[j], data, QPP_AT_UPPER_BOUND);
                        if (err != QPP_OK)
                        {
                            #ifdef QPP_WRITE_LOGFILE
                            if (log_level > 1)
                            {
                                fprintf(data->logfile, "%s\tError: fixVariable failed\n",
                                        indent);
                                fflush(data->logfile);
                            }
                            #endif
                            return err;
                        }
                    }
                }
            }
            else if (qppIsGreater(g[j], alpha*g[k], QPP_ZERO))
            {
                if ( ((zl[k] >= 0.0) && (alpha > 0.0)) ||
                     ((zu[k] <= 0.0) && (alpha < 0.0)) )
                {
                    data->found_reduction = 1;
                    if (isinf(xl_medium[j]))
                    {
                        #ifdef QPP_WRITE_LOGFILE
                        if (log_level > 1)
                        {
                            fprintf(data->logfile, "%s\tError: QP is dual infeasible\n",
                                    indent);
                            if (log_level > 2)
                            {
                                fprintf(data->logfile, "%s\t\tk = %" QPP_PRID ", j = %"
                                        QPP_PRID ", zl = %.6e, zu = %.6e, alpha = "
                                        "%.6e, <\n", indent, k, j, zl[k], zu[k], alpha);
                            }
                            fflush(data->logfile);
                        }
                        #endif
                        return QPP_DUAL_INFEASIBLE;
                    }
                    else
                    {
                        err = fixVariable(j, xl_medium[j], data, QPP_AT_LOWER_BOUND);
                        if (err != QPP_OK)
                        {
                            #ifdef QPP_WRITE_LOGFILE
                            if (log_level > 1)
                            {
                                fprintf(data->logfile, "%s\tError: fixVariable failed\n",
                                        indent);
                                fflush(data->logfile);
                            }
                            #endif
                            return err;
                        }
                    }
                }
            }
        }
    }

    #ifdef QPP_WRITE_LOGFILE
    if (log_level > 1)
    {
        fprintf(data->logfile, "%s\tok\n\n", indent);
        fflush(data->logfile);
    }
    #endif

    return QPP_OK;
}


static qpp_return_value_t sparsification(qpp_data_t *const data)
{
    qpp_return_value_t err;
    qpp_int_t i, j, jj, k, kk, it, min_nnz_column_nnz, min_nnz_column_index,
        has_subpattern, pivot_column_nnz, pivot_row_length, col_length, pivot_column_index;
    qpp_int_t *ac, *av, *pivot_row_indices, *col_indices, *A_nnz_rows, *A_nnz_columns,
        *A_irowst, *A_jcn, *A_linkrw;
    qpp_real_t factor, eq_tol, stab_tol;
    qpp_real_t *al, *au, *yl, *yu, *pivot_row_values, *col_values, *A_x;
    qpp_minheap_t *heap;
    qpp_ecrmatrix_t *A;
    qpp_se_sparsification_t *spf;
    qpp_stack_t *stack;

    #ifdef QPP_WRITE_LOGFILE
    qpp_int_t log_level;
    char indent[] = "\t";
    #endif

    stab_tol = data->options->stab_tol;
    eq_tol = data->options->eq_tol;
    ac = data->ac;
    av = data->av;
    A_nnz_rows = data->A_nnz_rows;
    A_nnz_columns = data->A_nnz_columns;
    al = data->al;
    au = data->au;
    yl = data->yl;
    yu = data->yu;
    heap = data->eq_constr_queue;
    A = data->A;
    A_irowst = A->irowst;
    A_jcn = A->jcn;
    A_linkrw = A->linkrw;
    A_x = A->x;
    stack = data->presolve_stack;

    pivot_row_indices = data->mem_int;
    col_indices = &pivot_row_indices[data->mem_length];
    pivot_row_values = data->mem_real;
    col_values = &pivot_row_values[data->mem_length];

    #ifdef QPP_WRITE_LOGFILE
    log_level = (data->logfile != NULL) ? data->options->log_level : -1;
    if (log_level > 1)
    {
        fprintf(data->logfile, "%s>> sparsification:\n", indent);
    }
    #endif

    while (!mhIsEmpty(heap))
    {
        if (mhExtractMin(heap, &i) != QPP_OK)
        {
            #ifdef QPP_WRITE_LOGFILE
            if (log_level > 1)
            {
                fprintf(data->logfile, "%s\tError: mhExtractMin failed\n",
                        indent);
                fflush(data->logfile);
            }
            #endif
            return QPP_FATAL_ERROR;
        }

        if ( !ac[i] || !qppIsEqual(al[i], au[i], eq_tol) || (A_nnz_rows[i] == 0) )
        {
            continue;
        }

        getRowA(A, i, av, pivot_row_indices, pivot_row_values, &pivot_row_length);
        assert(pivot_row_length == A_nnz_rows[i]);

        /*  !!! The pivot element must fulfill (çf. Master's thesis "Preprocessing for
            Quadratic Programming", §3.11.2): |a_ij| >= tol * |a_kj| for some tolerance
            tol > 0 and some j. We prefer a pivot element with as less nonzeros in
            its column as possible, as long as it fulfills the above stated
            property. */

        /*  Search for the column with the least number of nonzero elements. */
        min_nnz_column_index = -1;
        min_nnz_column_nnz = QPP_MAX_INT;
        for (jj = 0; jj < pivot_row_length; ++jj)
        {
            j = pivot_row_indices[jj];
            if (A_nnz_columns[j] <= min_nnz_column_nnz)
            {
                min_nnz_column_index = j;
                min_nnz_column_nnz = A_nnz_columns[j];
            }
        }
        assert(min_nnz_column_index >= 0);

        /* We extract the column with the least number of nonzero elements. */
        ac[i] = 0;
        getColA(A, min_nnz_column_index, ac, col_indices, col_values, &col_length);
        ac[i] = 1;

        /*  We check manually (without getting the whole row k of A at once)
            if row k is a superset of row i */
        for (kk = 0; kk < col_length; ++kk)
        {
            k = col_indices[kk];
            if (A_nnz_rows[i] > A_nnz_rows[k])
            {
                continue;
            }

            has_subpattern = 1;
            factor = 0.0;

            pivot_column_index = -1;
            pivot_column_nnz = QPP_MAX_INT;

            /*  If the constraint matrix is not sorted row-wise, then put the
                following line into the for loop!!! */
            it = A_irowst[k];

            /* Check if the nonzero pattern of row i is a subset of the pattern of row k. */
            for (jj = 0; jj < pivot_row_length; ++jj)
            {
                j = pivot_row_indices[jj];

                while ( (it >= 0) && (A_jcn[it] != j) )
                {
                    it = A_linkrw[it];
                }

                if ( (it < 0) || (A_x[it] == 0.0) || !av[A_jcn[it]] )
                {
                    has_subpattern = 0;
                    break;
                }

                /* Stability criterion satisfied? */
                if (fabs(pivot_row_values[jj] / A_x[it]) >= stab_tol)
                {
                    if (fabs(factor) > 0.0)
                    {
                        /*  We already have a candidate. Only change pivot element if the
                            number of nonzeros in its column is smaller. */
                        if (A_nnz_columns[j] < pivot_column_nnz)
                        {
                            pivot_column_index = j;
                            pivot_column_nnz = A_nnz_columns[j];
                            factor = A_x[it] / pivot_row_values[jj];
                        }
                    }
                    else
                    {
                        /* This is the first element that satisfies the criterion */
                        pivot_column_index = j;
                        pivot_column_nnz = A_nnz_columns[j];
                        factor = A_x[it] / pivot_row_values[jj];
                    }
                }
            }

            if ( !has_subpattern || (pivot_column_index < 0) )
            {
                continue;
            }

            data->found_reduction = 1;

            err = qppStackReallocIf(stack, 1);
            if (err != QPP_OK)
            {
                return err;
            }

            spf = (qpp_se_sparsification_t*) malloc(sizeof(qpp_se_sparsification_t));
            if (spf == NULL)
            {
                return QPP_OUT_OF_MEMORY;
            }

            /* Nonzero pattern of row i is a subset of the pattern of row k! */
            it = A_irowst[k];
            for (jj = 0; jj < pivot_row_length; ++jj)
            {
                j = pivot_row_indices[jj];

                while (A_jcn[it] != j)
                {
                    it = A_linkrw[it];
                }

                /* At least a_kj with j = <pivot_column_index> will vanish. */
                if ( qppIsEqual(A_x[it], factor * pivot_row_values[jj], eq_tol) )
                {
                    A_x[it] = 0.0;
                    --A_nnz_rows[k];
                    --A_nnz_columns[j];
                }
                else
                {
                    A_x[it] -= factor * pivot_row_values[jj];
                }
            }

            /* Update heap of equality constraints */
            if (qppIsEqual(al[k], au[k], eq_tol))
            {
                if (mhSearchForValue(heap, k, &j) == QPP_OK)
                {
                    mhDecreaseKey(heap, j, A_nnz_rows[k]);
                }
                else
                {
                    mhInsert(heap, A_nnz_rows[k], k);
                }
            }

            spf->al = al[k];
            spf->au = au[k];
            spf->i = i;
            spf->k = k;
            spf->j = pivot_column_index;
            spf->alpha = factor;

            al[k] -= factor * al[i];
            au[k] -= factor * au[i];

            if (factor > 0.0)
            {
                yl[i] += factor * yl[k];
                yu[i] += factor * yu[k];
            }
            else
            {
                yl[i] += factor * yu[k];
                yu[i] += factor * yl[k];
            }
            qppStackPush(stack, spf, QPP_SE_SPARSIFICATION);
        }
    }

    #ifdef QPP_WRITE_LOGFILE
    if (log_level > 1)
    {
        fprintf(data->logfile, "%s\tok\n\n", indent);
        fflush(data->logfile);
    }
    #endif

    return QPP_OK;
}


static qpp_return_value_t tightenVarBoundsPC(qpp_data_t *const data)
{
    qpp_return_value_t err;
    qpp_int_t i, j, jj, ub_length, lb_length, row_length, old_ub_length,
        old_lb_length;
    qpp_int_t *ac, *av, *A_nnz_rows, *lb_indices, *ub_indices,
        *row_indices, *old_lb_indices, *old_ub_indices;
    qpp_real_t nubc, nlbc, stab_tol;
    qpp_real_t *xl, *xu, *al, *au, *lb_values, *ub_values, *zl, *zu,
        *row_values, *old_lb_values, *old_ub_values;

    qpp_se_new_bounds_pc_t *nbpc;
    qpp_stack_t *stack;
    qpp_ecrmatrix_t *A;

    #ifdef QPP_WRITE_LOGFILE
    qpp_int_t log_level;
    char indent[] = "\t";
    #endif

    A_nnz_rows = data->A_nnz_rows;
    ac = data->ac;
    av = data->av;
    xl = data->xl_tight;
    xu = data->xu_tight;
    al = data->al;
    au = data->au;
    zl = data->zl;
    zu = data->zu;
    stack = data->presolve_stack;
    A = data->A;
    stab_tol = data->options->stab_tol;

    row_indices = data->mem_int;
    old_lb_indices = &row_indices[data->mem_length];
    old_ub_indices = &old_lb_indices[data->mem_length];
    row_values = data->mem_real;
    old_lb_values = &row_values[data->mem_length];
    old_ub_values = &old_lb_values[data->mem_length];

    lb_indices = data->mem_for_bounds_int;
    ub_indices = &lb_indices[data->mem_length];
    lb_values = data->mem_for_bounds_real;
    ub_values = &lb_values[data->mem_length];

    #ifdef QPP_WRITE_LOGFILE
    log_level = (data->logfile != NULL) ? data->options->log_level : -1;
    if (log_level > 1)
    {
        fprintf(data->logfile, "%s>> tightenVarBoundsPC:\n", indent);
    }
    #endif

    for (i = 0; i < data->m; ++i)
    {
        if ( !ac[i] || (A_nnz_rows[i] <= 1) )
        {
            continue;
        }

        getRowA(A, i, av, row_indices, row_values, &row_length);
        err = varBoundsPCi(au[i], al[i], xl, xu, row_indices, row_values, row_length,
                           lb_indices, lb_values, &lb_length, ub_indices, ub_values,
                           &ub_length);
        if (err != QPP_OK)
        {
            #ifdef QPP_WRITE_LOGFILE
            if (log_level > 1)
            {
                fprintf(data->logfile, "%s\tError: varBoundsPCi failed. Error code = %"
                        QPP_PRID "\n", indent, (qpp_int_t) err);
                fflush(data->logfile);
            }
            #endif
            return err;
        }

        old_lb_length = 0;
        old_ub_length = 0;

        for (jj = 0; jj < lb_length; ++jj)
        {
            j = lb_indices[jj];
            nlbc = lb_values[jj];

            if (qppIsGreater(nlbc, xl[j], stab_tol))
            {
                old_lb_indices[old_lb_length] = j;
                old_lb_values[old_lb_length++] = xl[j];
                xl[j] = nlbc;
                zu[j] = qppMinr(0.0, zu[j]);
            }
        }

        for (jj = 0; jj < ub_length; ++jj)
        {
            j = ub_indices[jj];
            nubc = ub_values[jj];

            if (qppIsGreater(xu[j], nubc, stab_tol))
            {
                old_ub_indices[old_ub_length] = j;
                old_ub_values[old_ub_length++] = xu[j];
                xu[j] = nubc;
                zl[j] = qppMaxr(0.0, zl[j]);
            }
        }

        if (old_lb_length + old_ub_length > 0)
        {
            data->found_reduction = 1;
            err = qppStackReallocIf(stack, 1);
            if (err != QPP_OK)
            {
                return err;
            }

            nbpc = (qpp_se_new_bounds_pc_t*) malloc(sizeof(qpp_se_new_bounds_pc_t));
            if (nbpc == NULL)
            {
                return QPP_OUT_OF_MEMORY;
            }

            nbpc->i = i;
            nbpc->lb_indices = NULL;
            nbpc->ub_indices = NULL;
            nbpc->lb_num = 0;
            nbpc->ub_num = 0;

            /* We do not save the old bounds (their values) since the bounds will not be restored. */
            nbpc->lb_values = NULL;
            nbpc->ub_values = NULL;

            if (old_lb_length > 0)
            {
                nbpc->lb_indices = (qpp_int_t*) malloc((size_t)(old_lb_length) * sizeof(qpp_int_t));
                if (nbpc->lb_indices == NULL)
                {
                    free(nbpc);
                    return QPP_OUT_OF_MEMORY;
                }
                memcpy(nbpc->lb_indices, old_lb_indices, (size_t)(old_lb_length) * sizeof(qpp_int_t));
                nbpc->lb_num = old_lb_length;
            }

            if (old_ub_length > 0)
            {
                nbpc->ub_indices = (qpp_int_t*) malloc((size_t)(old_ub_length) * sizeof(qpp_int_t));
                if (nbpc->ub_indices == NULL)
                {
                    free(nbpc->lb_indices);
                    free(nbpc);
                    return QPP_OUT_OF_MEMORY;
                }
                memcpy(nbpc->ub_indices, old_ub_indices, (size_t)(old_ub_length) * sizeof(qpp_int_t));
                nbpc->ub_num = old_ub_length;
            }
            qppStackPush(stack, nbpc, QPP_SE_NEWBOUNDS_PC);
        }
    }

    #ifdef QPP_WRITE_LOGFILE
    if (log_level > 1)
    {
        fprintf(data->logfile, "%s\tok\n\n", indent);
        fflush(data->logfile);
    }
    #endif

    return QPP_OK;
}


static qpp_return_value_t scaleInf(qpp_data_t *const data,
                                   qpp_int_t max_iter,
                                   qpp_real_t *const work)
{
    qpp_return_value_t err;
    qpp_int_t m, n, num_iter, i, it, scaled;
    qpp_int_t *A_linkrw, *A_irn, *A_jcn, *Hlinkrw, *Hlinkcl, *H_irn, *H_jcn, *ac, *av;
    int E;
    qpp_real_t *H_x, *A_x, *R, *C, *RR, *CC;
    qpp_ecrmatrix_t *A, *H;
    qpp_se_scaling_t *scal;

    if (data == NULL)
    {
        return QPP_NULL_ARGUMENT;
    }

    err = qppStackReallocIf(data->presolve_stack, 1);
    if (err != QPP_OK)
    {
        return err;
    }

    A = data->A;
    A_irn = A->irn;
    A_jcn = A->jcn;
    A_linkrw = A->linkrw;
    A_x = A->x;
    H = data->H;
    H_irn = H->irn;
    H_jcn = H->jcn;
    H_x = H->x;
    Hlinkrw = H->linkrw;
    Hlinkcl = H->linkcl;
    m = A->nrow;
    n = A->ncol;
    ac = data->ac;
    av = data->av;

    RR = work;
    CC = &RR[m];

    /* R stores row scaling factors, C column scaling factors. */
    R = (qpp_real_t*) malloc((size_t)(m) * sizeof(qpp_real_t));
    C = (qpp_real_t*) malloc((size_t)(n) * sizeof(qpp_real_t));
    scal = (qpp_se_scaling_t*) malloc(sizeof(qpp_se_scaling_t));
    if ( (C == NULL) || (R == NULL) || (scal == NULL) )
    {
        free(C);
        free(R);
        free(scal);
        return QPP_OUT_OF_MEMORY;
    }
    qppSetArrayr(R, m, 1.0);
    qppSetArrayr(C, n, 1.0);

    for (num_iter = 0; num_iter < max_iter; ++num_iter)
    {
        scaled = 1;

        for (i = 0; i < n; ++i)
        {
            /* Compute inf-norm of every column of H */
            CC[i] = 0.0;
            if (av[i])
            {
                it = H->irowst[i];
                while (it >= 0)
                {
                    if (av[H_jcn[it]])
                    {
                        CC[i] = qppMaxr(CC[i], fabs(H_x[it]));
                    }
                    it = Hlinkrw[it];
                }
                it = H->jcolst[i];
                while (it >= 0)
                {
                    if (av[H_irn[it]])
                    {
                        CC[i] = qppMaxr(CC[i], fabs(H_x[it]));
                    }
                    it = Hlinkcl[it];
                }
            }
        }

        for (i = 0; i < m; ++i)
        {
            /*  Compute inf-norm of the rows and columns of A. The inf-norms of
                the columns will be compared to the inf-norms of the Hessian matrix */
            RR[i] = 0.0;
            if (ac[i])
            {
                it = A->irowst[i];
                while (it >= 0)
                {
                    if (av[A_jcn[it]])
                    {
                        RR[i] = qppMaxr(RR[i], fabs(A_x[it]));
                        CC[A_jcn[it]] = qppMaxr(CC[A_jcn[it]], fabs(A_x[it]));
                    }
                    it = A_linkrw[it];
                }
            }
        }

        /* Compute scaling factors. */
        for (i = 0; i < n; ++i)
        {
            if (av[i])
            {
                frexp(CC[i], &E);
                if (E < 0)
                {
                    CC[i] = ldexp(1.0, (1-E)/2);
                    scaled = 0;
                }
                else if (E > 1)
                {
                    CC[i] = ldexp(1.0, (-E/2));
                    scaled = 0;
                }
                else
                {
                    CC[i] = 1.0;
                }
                C[i] *= CC[i];
            }
        }

        for (i = 0; i < m; ++i)
        {
            if (ac[i])
            {
                frexp(RR[i], &E);
                if (E < 0)
                {
                    RR[i] = ldexp(1.0, (1-E)/2);
                    scaled = 0;
                }
                else if (E > 1)
                {
                    RR[i] = ldexp(1.0, (-E/2));
                    scaled = 0;
                }
                else
                {
                    RR[i] = 1.0;
                }
                R[i] *= RR[i];
            }
        }

        if (scaled) break;

        for (i = 0; i < H->nnz; ++i)
        {
            if ( av[H_irn[i]] && av[H_jcn[i]] )
            {
                H_x[i] *= CC[H_irn[i]] * CC[H_jcn[i]];
            }
        }

        for (i = 0; i < A->nnz; ++i)
        {
            if ( ac[A_irn[i]] && av[A_jcn[i]] )
            {
                A_x[i] *= RR[A_irn[i]] * CC[A_jcn[i]];
            }
        }
    }

    /*  Undo scaling when it failed. */
    if (num_iter > max_iter)
    {
        for (i = 0; i < H->nnz; ++i)
        {
            if ( av[H_irn[i]] && av[H_jcn[i]] )
            {
                H_x[i] /= C[H_irn[i]] * C[H_jcn[i]];
            }
        }

        for (i = 0; i < A->nnz; ++i)
        {
            if ( ac[A_irn[i]] && av[A_jcn[i]] )
            {
                A_x[i] /= R[A_irn[i]] * C[A_jcn[i]];
            }
        }
        free(R);
        free(C);
        free(scal);
        return QPP_BAD_SCALING;
    }

    /* Scale QP: */
    for (i = 0; i < n; ++i)
    {
        if (av[i])
        {
            data->xl_tight[i] /= C[i];
            data->xu_tight[i] /= C[i];
            data->xl_medium[i] /= C[i];
            data->xu_medium[i] /= C[i];
            data->zl[i] *= C[i];
            data->zu[i] *= C[i];
            data->g[i] *= C[i];
            data->H_diag[i] *= C[i] * C[i];
        }
    }

    for (i = 0; i < m; ++i)
    {
        if (ac[i])
        {
            data->al[i] *= R[i];
            data->au[i] *= R[i];
            data->yl[i] /= R[i];
            data->yu[i] /= R[i];
        }
    }

    scal->C = C;
    scal->R = R;

    qppStackPush(data->presolve_stack, scal, QPP_SE_SCALING);

    return QPP_OK;
}


static qpp_return_value_t varBoundsPCi(const qpp_real_t au_i,
                                       const qpp_real_t al_i,
                                       const qpp_real_t xl[],
                                       const qpp_real_t xu[],
                                       qpp_int_t row_indices[],
                                       qpp_real_t row_values[],
                                       const qpp_int_t row_length,
                                       qpp_int_t lb_indices[],
                                       qpp_real_t lb_values[],
                                       qpp_int_t *lb_length,
                                       qpp_int_t ub_indices[],
                                       qpp_real_t ub_values[],
                                       qpp_int_t *ub_length)
{
    qpp_int_t lb_length__, ub_length__, j, jj, jj_u, jj_l, iubpc, ilbpc, iubnc, ilbnc,
        neg_length, pos_length;
    qpp_real_t a_ij, u, l, nubc, nlbc;

    lb_length__ = 0;
    ub_length__ = 0;
    *lb_length = 0;
    *ub_length = 0;
    /* iubpc = (Number of) Infinite Upper Bound Positive Coefficient, etc. */
    iubpc = 0; ilbpc = 0; iubnc = 0; ilbnc = 0;
    jj_l = -1; jj_u = -1;
    u = 0.0;
    l = 0.0;

    /*  Sort row such that the first entries are negative and the following
        entries are positive */
    qppSortNegPos(row_indices, row_values, row_length, &neg_length, &pos_length);

    for (jj = 0; jj < neg_length; ++jj)
    {
        j = row_indices[jj];
        a_ij = row_values[jj];

        if (isinf(xl[j]))
        {
            ilbnc++;
            jj_u = jj;
        }
        else
        {
            u += a_ij * xl[j];
        }

        if (isinf(xu[j]))
        {
            iubnc++;
            jj_l = jj;
        }
        else
        {
            l += a_ij * xu[j];
        }
    }

    for (jj = neg_length; jj < row_length; ++jj)
    {
        j = row_indices[jj];
        a_ij = row_values[jj];

        if (isinf(xu[j]))
        {
            iubpc++;
            jj_u = jj;
        }
        else
        {
            u += a_ij * xu[j];
        }

        if (isinf(xl[j]))
        {
            ilbpc++;
            jj_l = jj;
        }
        else
        {
            l += a_ij * xl[j];
        }
    }

    if ( qppIsEqual(u, al_i, 1e-8) || qppIsEqual(l, au_i, 1e-8)  || ((au_i >= u) && (al_i <= l)) )
        return QPP_OK;  /* Forcing or redundant constraint! */

    if (iubpc+ilbnc == 0)
    {
        /* We can deduce finite bounds on ALL variables in row <i> */
        for (jj = 0; jj < neg_length; ++jj)
        {
            j = row_indices[jj];
            a_ij = row_values[jj];
            nubc = xl[j] + al_i/a_ij - u/a_ij;
            ub_indices[ub_length__] = j;
            ub_values[ub_length__++] = nubc;
        }

        for (jj = neg_length; jj < row_length; ++jj)
        {
            j = row_indices[jj];
            a_ij = row_values[jj];
            nlbc = xu[j] + al_i/a_ij - u/a_ij;
            lb_indices[lb_length__] = j;
            lb_values[lb_length__++] = nlbc;
        }
    }
    else if (iubpc+ilbnc == 1)
    {
        /* One bound is infinite so we can use Gondzio's method */
        assert(jj_u != -1);
        a_ij = row_values[jj_u];
        j = row_indices[jj_u];

        if (a_ij > 0.0)
        {
            nlbc = al_i/a_ij - u/a_ij;
            lb_indices[lb_length__] = j;
            lb_values[lb_length__++] = nlbc;
        }
        else
        {
            nubc = al_i/a_ij - u/a_ij;
            ub_indices[ub_length__] = j;
            ub_values[ub_length__++] = nubc;
        }
    }

    if (ilbpc+iubnc == 0)
    {
        /* We can deduce finite bounds on ALL variables in row <i> */
        for (jj = 0; jj < neg_length; ++jj)
        {
            j = row_indices[jj];
            a_ij = row_values[jj];
            nlbc = xu[j] + au_i/a_ij - l/a_ij;
            lb_indices[lb_length__] = j;
            lb_values[lb_length__++] = nlbc;
        }

        for (jj = neg_length; jj < row_length; ++jj)
        {
            j = row_indices[jj];
            a_ij = row_values[jj];
            nubc = xl[j] + au_i/a_ij - l/a_ij;
            ub_indices[ub_length__] = j;
            ub_values[ub_length__++] = nubc;
        }
    }
    else if (ilbpc+iubnc == 1)
    {
        /* One bound is infinite -> Gondzio's method */
        assert(jj_l != -1);
        a_ij = row_values[jj_l];
        j = row_indices[jj_l];

        if (a_ij > 0.0)
        {
            nubc = au_i/a_ij - l/a_ij;
            ub_indices[ub_length__] = j;
            ub_values[ub_length__++] = nubc;
        }
        else
        {
            nlbc = au_i/a_ij - l/a_ij;
            lb_indices[lb_length__] = j;
            lb_values[lb_length__++] = nlbc;
        }
    }
    *lb_length = lb_length__;
    *ub_length = ub_length__;

    return QPP_OK;
}


static qpp_return_value_t dualConstraints(qpp_data_t *const data)
{
    qpp_return_value_t err;
    qpp_int_t j, i, ii, A_pos_length, A_neg_length, H_pos_length, H_neg_length,
        H_length_col, A_length_col;
    qpp_int_t *ac, *av, *A_column_indices, *H_column_indices, *A_nnz_columns, *H_nnz_columns;
    qpp_real_t l, u, a, h, eq_tol, nlbc, nubc;
    qpp_real_t *xl, *xu, *yl, *yu, *A_column_values, *H_column_values, *zl, *zu,
        *g, *px_pos, *px_neg, *py_pos, *py_neg;
    qpp_active_type_t neg_actv_type, pos_actv_type;
    qpp_ecrmatrix_t *A, *H;

    #ifdef QPP_WRITE_LOGFILE
    qpp_int_t log_level;
    char indent[] = "\t";
    #endif // QPP_WRITE_LOGFILE

    ac = data->ac;
    av = data->av;
    xl = data->xl_tight;
    xu = data->xu_tight;
    yl = data->yl;
    yu = data->yu;
    zl = data->zl;
    zu = data->zu;
    A = data->A;
    H = data->H;
    g = data->g;
    eq_tol = data->options->eq_tol;
    A_nnz_columns = data->A_nnz_columns;
    H_nnz_columns = data->H_nnz_columns;

    A_column_indices = data->mem_int;;
    H_column_indices = &A_column_indices[data->mem_length];
    A_column_values = data->mem_real;
    H_column_values = &A_column_values[data->mem_length];

    #ifdef QPP_WRITE_LOGFILE
    log_level = (data->logfile != NULL) ? data->options->log_level : -1;
    if (log_level > 1)
    {
        fprintf(data->logfile, "%s>> dualConstraints:\n", indent);
    }
    #endif

    for (j = 0; j < data->n; ++j)
    {
        if (!av[j])
        {
            continue;
        }

        l = 0.0;
        u = 0.0;
        A_pos_length = 0;
        H_pos_length = 0;
        A_neg_length = 0;
        H_neg_length = 0;

        getColA(A, j, ac, A_column_indices, A_column_values, &A_length_col);
        getColH(H, j, av, H_column_indices, H_column_values, &H_length_col);

        qppSortNegPos(A_column_indices, A_column_values, A_length_col, &A_neg_length, &A_pos_length);
        qppSortNegPos(H_column_indices, H_column_values, H_length_col, &H_neg_length, &H_pos_length);

        for (ii = 0; ii < A_neg_length; ++ii)
        {
            i = A_column_indices[ii];
            a = A_column_values[ii];
            u -= a * yu[i];
            l -= a * yl[i];
        }

        for (ii = A_neg_length; ii < A_length_col; ++ii)
        {
            i = A_column_indices[ii];
            a = A_column_values[ii];
            u -= a * yl[i];
            l -= a * yu[i];
        }

        for (ii = 0; ii < H_neg_length; ++ii)
        {
            i = H_column_indices[ii];
            h = H_column_values[ii];
            u += h * xl[i];
            l += h * xu[i];
        }

        for (ii = H_neg_length; ii < H_length_col; ++ii)
        {
            i = H_column_indices[ii];
            h = H_column_values[ii];
            u += h * xu[i];
            l += h * xl[i];
        }

        if ( (qppIsEqual(u + g[j], zl[j], eq_tol) || qppIsEqual(l + g[j], zu[j], eq_tol)) &&
                ((A_nnz_columns[j] > 1) || (H_nnz_columns[j] > 0)) )
        {
            if (H_pos_length + H_neg_length > 0)
            {
                data->found_reduction = 1;
                #ifdef QPP_WRITE_LOGFILE
                if (log_level > 1)
                {
                    fprintf(data->logfile, "%s\t\tHessian Forcing dual constraint\n", indent);
                }
                #endif
            }
            #ifdef QPP_WRITE_LOGFILE
            if (log_level > 1)
            {
                fprintf(data->logfile, "%s\t\tForcing dual constraint\n", indent);
            }
            #endif

            // Dual constraint is forcing! Fix variables x and multipliers y!
            if (qppIsEqual(u + g[j], zl[j], eq_tol))
            {
                px_pos = xu;
                px_neg = xl;
                py_pos = yu;
                py_neg = yl;
                neg_actv_type = QPP_AT_LOWER_BOUND;
                pos_actv_type = QPP_AT_UPPER_BOUND;
            }
            else
            {
                px_pos = xl;
                px_neg = xu;
                py_pos = yl;
                py_neg = yu;
                neg_actv_type = QPP_AT_UPPER_BOUND;
                pos_actv_type = QPP_AT_LOWER_BOUND;
            }

            for (ii = 0; ii < H_neg_length; ++ii)
            {
                i = H_column_indices[ii];
                err = fixVariable(i, px_neg[i], data, neg_actv_type);
                if (err != QPP_OK)
                {
                    #ifdef QPP_WRITE_LOGFILE
                    if (log_level > 1)
                    {
                        fprintf(data->logfile, "%s\tError: fixVariable failed\n", indent);
                        fflush(data->logfile);
                    }
                    #endif
                    return err;
                }
            }

            for (ii = H_neg_length; ii < H_length_col; ++ii)
            {
                i = H_column_indices[ii];
                err = fixVariable(i, px_pos[i], data, pos_actv_type);
                if (err != QPP_OK)
                {
                    #ifdef QPP_WRITE_LOGFILE
                    if (log_level > 1)
                    {
                        fprintf(data->logfile, "%s\tError: fixVariable failed\n", indent);
                        fflush(data->logfile);
                    }
                    #endif
                    return err;
                }
            }

            for (ii = 0; ii < A_neg_length; ++ii)
            {
                i = A_column_indices[ii];
                py_neg[i] = py_pos[i];
            }

            for (ii = A_neg_length; ii < A_length_col; ++ii)
            {
                i = A_column_indices[ii];
                py_pos[i] = py_neg[i];
            }
        }

        if (av[j])
        {
            if (u + g[j] < -eq_tol)
            {
                if (zl[j] >= 0.0)
                {
                    #ifdef QPP_WRITE_LOGFILE
                    if (log_level > 1)
                    {
                        fprintf(data->logfile, "%s\tError: QP is dual infeasible\n",
                                indent);
                        if (log_level > 2)
                        {
                            fprintf(data->logfile, "%s\t\tj = %" QPP_PRID ", zl = %.6e"
                                    ", zu = %.6e <\n", indent, j, zl[j], zu[j]);
                        }
                        fflush(data->logfile);
                    }
                    #endif
                    return QPP_DUAL_INFEASIBLE;
                }
                // z[j] < 0 -> Fix at upper bound (will not be inf since zl[j] < 0)
                data->found_reduction = 1;
                err = fixVariable(j, xu[j], data, QPP_AT_UPPER_BOUND);
                if (err != QPP_OK)
                {
                    #ifdef QPP_WRITE_LOGFILE
                    if (log_level > 1)
                    {
                        fprintf(data->logfile, "%s\tError: fixVariable failed\n", indent);
                        fflush(data->logfile);
                    }
                    #endif
                    return err;
                }
            }
            else if (l + g[j] > eq_tol)
            {
                if (zu[j] <= 0.0)
                {
                    #ifdef QPP_WRITE_LOGFILE
                    if (log_level > 1)
                    {
                        fprintf(data->logfile, "%s\tError: QP is dual infeasible\n",
                                indent);
                        if (log_level > 2)
                        {
                            fprintf(data->logfile, "%s\t\tj = %" QPP_PRID ", zl = %.6e"
                                    ", zu = %.6e >\n", indent, j, zl[j], zu[j]);
                        }
                        fflush(data->logfile);
                    }
                    #endif
                    return QPP_DUAL_INFEASIBLE;
                }
                // z[j] > 0 -> Fix at lower bound (will not be inf since zu[j] > 0)
                data->found_reduction = 1;
                err = fixVariable(j, xl[j], data, QPP_AT_LOWER_BOUND);
                if (err != QPP_OK)
                {
                    #ifdef QPP_WRITE_LOGFILE
                    if (log_level > 1)
                    {
                        fprintf(data->logfile, "%s\tError: fixVariable failed\n", indent);
                        fflush(data->logfile);
                    }
                    #endif
                    return err;
                }
            }

            if ( (zl[j] >= 0.0) && isfinite(u) )
            {
                for (ii = 0; ii < A_neg_length; ++ii)
                {
                    i = A_column_indices[ii];
                    a = A_column_values[ii];
                    nlbc = yu[i] + u/a + g[j]/a;
                    if (nlbc > yl[i])
                    {
                        yl[i] = nlbc;
                    }
                }

                for (ii = A_neg_length; ii < A_length_col; ++ii)
                {
                    i = A_column_indices[ii];
                    a = A_column_values[ii];
                    nubc = yl[i] + u/a + g[j]/a;
                    if (nubc < yu[i])
                    {
                        yu[i] = nubc;
                    }
                }
            }

            if ( (zu[j] <= 0.0) && isfinite(l) )
            {
                for (ii = 0; ii < A_neg_length; ++ii)
                {
                    i = A_column_indices[ii];
                    a = A_column_values[ii];
                    nubc = yl[i] + l/a + g[j]/a;
                    if (nubc < yu[i])
                    {
                        yu[i] = nubc;
                    }
                }

                for (ii = A_neg_length; ii < A_length_col; ++ii)
                {
                    i = A_column_indices[ii];
                    a = A_column_values[ii];
                    nlbc = yu[i] + l/a + g[j]/a;
                    if (nlbc > yl[i])
                    {
                        yl[i] = nlbc;
                    }
                }
            }
        }
    }

    #ifdef QPP_WRITE_LOGFILE
    if (log_level > 1)
    {
        fprintf(data->logfile, "%s\tok\n\n", indent);
        fflush(data->logfile);
    }
    #endif

    return QPP_OK;
}
