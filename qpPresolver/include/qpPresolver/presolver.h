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
 *	\file include/qpPresolver/presolver.h
 *	\author Dominik Cebulla
 *	\version 1.0 Beta
 *	\date 2017
 *
 *  Declaration of the presolver data structure and the corresponding interface
 *  functions.
 *
 */


#ifndef QPPRESOLVER_PRESOLVER_H
#define QPPRESOLVER_PRESOLVER_H

#include <assert.h>
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <string.h>
#include <time.h>

#include <qpPresolver/constants.h>
#include <qpPresolver/ecrmatrix.h>
#include <qpPresolver/minheap.h>
#include <qpPresolver/options.h>
#include <qpPresolver/stack.h>
#include <qpPresolver/types.h>
#include <qpPresolver/utility.h>


/** \brief Struct containing necessary (and helpful) data for presolving a QP.
 *
 *  This struct stores the QP (matrices and vectors), several arrays
 *  containing useful information for better runtime performance, additional
 *  memory, internal status, time, etc.
 */
typedef struct data
{
    /** \name QP data. All computations are performed on this data. */
    /**@{*/
    qpp_ecrmatrix_t *A;         /**< Constraint matrix. */
    qpp_ecrmatrix_t *H;         /**< Hessian matrix. */
    qpp_real_t *g;          /**< Gradient vector. */
    qpp_real_t f;           /**< Constant term of objective function. */
    qpp_real_t *xl;         /**< Lower bounds on the variables. */
    qpp_real_t *xu;         /**< Upper bounds on the variables. */
    qpp_real_t *al;         /**< Lower bounds on linear constraints. */
    qpp_real_t *au;         /**< Upper bounds on linear constraints. */
    qpp_int_t m;            /**< Number of linear constraints (>= 0). */
    qpp_int_t n;            /**< Number of variables (> 0). */
    /**@}*/

    /** \name Additional bounds on the primal and dual variables. */
    /**@{*/
    qpp_real_t *xl_tight;       /**< Tightest lower bounds on the variables. */
    qpp_real_t *xl_medium;      /**< Medium lower bounds on the variables. */
    qpp_real_t *xu_tight;       /**< Tightest upper bounds on the variables. */
    qpp_real_t *xu_medium;      /**< Medium lower bounds on the variables. */
    qpp_real_t *yl;             /**< Lower bounds on the multipliers corresponding to
                                     the linear constraints. */
    qpp_real_t *yu;             /**< Upper bounds on the multipliers corresponding to
                                     the linear constraints. */
    qpp_real_t *zl;             /**< Lower bounds on the multipliers corresponding to
                                     the bound constraints. */
    qpp_real_t *zu;             /**< Upper bounds on the multipliers corresponding to
                                     the bound constraints. */
    /**@}*/

    /** \name Additional information.
     *
     *  E.g. the number of nonzero elements in each row/column of the matrices
     *  and information which variable/constraint is removed.
     */
    /**@{*/
    qpp_int_t *av;                  /**< Array indicating which variables are still
                                         part of the QP (= active variables). */
    qpp_int_t *ac;                  /**< Array indicating which constraints are still
                                         part of the QP (= active constraints). */
    qpp_int_t *A_nnz_rows;          /**< Number of nonzero elements in every row of A. */
    qpp_int_t *A_nnz_columns;       /**< Number of nonzero elements in every column of A. */
    qpp_int_t *H_nnz_columns;       /**< Number of nonzero elements in every column of H. */
    qpp_real_t *H_diag;             /**< Array storing the main diagonal of H. */
    /**@}*/

    /** \name Working storage.
     *
     *  Primarily for storing rows/columns of a matrix in sparse format.
     */
    /**@{*/
    qpp_int_t *mem_mpn_i;               /**< Memory of length m+n for storing integer values. */
    qpp_real_t *mem_mpn_r;              /**< Memory of length m+n for storing real values. */
    qpp_int_t mem_length;               /**< Length of one unit of extra memory. */
    qpp_int_t *mem_int;                 /**< Memory for storing integer values (4 units). */
    qpp_real_t *mem_real;               /**< Memory for storing real values. */
    qpp_int_t *mem_for_bounds_int;      /**< Memory for storing integer values when
                                             computing bounds on variables (2 units). */
    qpp_real_t *mem_for_bounds_real;    /**< Memory for storing real values when
                                             computing bounds on variables (2 units). */
    /**@}*/

    qpp_int_t found_reduction;  /**< Variable indicating whether a reduction during
                                     an iteration of the presolving loop was possible. */
    qpp_int_t num_iter;         /**< Number of iterations of the presolving loop until
                                     termination. */

    qpp_options_t *options;     /**< Options. */


    /** \name Data of the presolved QP.
     *
     *  Depending on which function is used, either the user or the presolver allocates
     *  all necessary memory. In the latter case, the presolver also deallocates the
     *  memory. The user has to allocate memory if \p qppGetPresolvedQP() is used, the
     *  presolver allocates memory if \p qppGetPresolvedQPPtr() is used.
     *  Matrices are stored in coordinate scheme. Only the lower triangular part of
     *  the Hessian matrix is returned.
     */
    /**@{*/
    qpp_int_t m_ps;             /**< Number of constraints of the presolved QP. */
    qpp_int_t n_ps;             /**< Number of variables of the presolved QP. */
    qpp_int_t *A_irn_ps;        /**< Row subscripts of the presolved constraint matrix. */
    qpp_int_t *A_jcn_ps;        /**< Column subscripts of the presolved constraint matrix. */
    qpp_real_t *A_x_ps;         /**< Nonzero elements of the presolved constraint matrix. */
    qpp_int_t A_nnz_ps;         /**< Number of nonzeros in the presolved constraint matrix. */
    qpp_int_t *H_irn_ps;        /**< Row subscripts of the presolved Hessian matrix. */
    qpp_int_t *H_jcn_ps;        /**< Column subscripts of the presolved Hessian matrix. */
    qpp_real_t *H_x_ps;         /**< Nonzero elements of the presolved Hessian matrix. */
    qpp_int_t H_nnz_ps;         /**< Number of nonzeros in the presolved Hessian matrix. */
    qpp_int_t H_diag_nnz_ps;    /**< Number of nonzeros in the main diagonal of the
                                     presolved Hessian matrix. */
    qpp_real_t *g_ps;           /**< Gradient vector of the presolved QP. */
    qpp_real_t f_ps;            /**< Constant term of the objective function of
                                     the presolved QP. */
    qpp_real_t *xl_ps;          /**< Lower bounds on the variables of the presolved QP. */
    qpp_real_t *xu_ps;          /**< Upper bounds on the variables of the presolved QP. */
    qpp_real_t *al_ps;          /**< Lower bounds on the lin. constr. of the presolved QP. */
    qpp_real_t *au_ps;          /**< Upper bounds on the lin. constr. of the presolved QP. */
    /**@}*/

    qpp_minheap_t *eq_constr_queue;     /**< Priority queue containing subscripts of
                                         all equality constraints. */
    qpp_stack_t *presolve_stack;    /**< Presolve stack storing necessary information of
                                         every preprocessing technique. */

    /** \name Quantities for measuring several (CPU) times. */
    /**@{*/
    clock_t init_time;              /**< CPU time for initializing the presolver. */
    clock_t presolve_time;          /**< CPU time for presolving the QP. */
    clock_t postsolve_time;         /**< CPU time for postsolving the QP. */
    /**@}*/

    qpp_status_t status;    /**< Internal status of the presolver. */

    #ifdef QPP_WRITE_LOGFILE
    FILE *logfile;          /**< Handle to logfile. */
    #endif
} qpp_data_t;


/*  ======================================================================================
    Interface functions
    ====================================================================================*/

/** \brief Creates new presolver entity and allocates memory.
 *
 *  All memory is allocated here, except for the presolved QP
 *  and for some temporary storage whose size is dependent on the
 *  maximum number of nonzero elements in every row/column of the matrices.
 *
 *  \param pdata New presolver entity. Will be set to NULL on failure, otherwise will
 *      contain a valid pointer to the new entity.
 *  \param m Number of linear constraints of the QP.
 *  \param n Number of variables of the QP.
 *  \param A_nnz Estimate of the number of nonzero elements in the constraint matrix
 *      of the QP. If necessary, memory will be reallocated by the presolver later.
 *  \param H_nnz Estimate of the number of nonzero elements of the lower triangular part (!)
 *      of the Hessian matrix. If necessary, memory will be reallocated by the presolver later.
 *  \param stack_size Initial size of the presolve stack. Whenever the stack has to be
 *      enlarged, it will be enlarged by \p stack_size elements.
 *
 *  \return QPP_OK \n QPP_OUT_OF_MEMORY \n QPP_INVALID_ARGUMENT \n QPP_NULL_ARGUMENT
 */
qpp_return_value_t qppNew(qpp_data_t **pdata,
                          const qpp_int_t m,
                          const qpp_int_t n,
                          const qpp_int_t A_nnz,
                          const qpp_int_t H_nnz,
                          const qpp_int_t stack_size);


/** \brief Initializes data of presolver entity with given QP.
 *
 *  Matrices must be passed in coordinate format. If a matrix does not contain any nonzero
 *  elements, then the respective parameters can be null pointers. Furthermore, some
 *  additional memory for temporary storage is allocated. If the estimates on the number
 *  of nonzero elements in qppNew() were not sufficient, then memory is reallocated for
 *  the respective matrices if necessary.
 *
 *  \param data Valid presolver entity.
 *  \param A_irn Row subscripts of the constraint matrix (coordinate scheme).
 *  \param A_jcn Column subscripts of the constraint matrix (coordinate scheme).
 *  \param A_x Nonzero elements of the constraint matrix (coordinate scheme).
 *  \param A_nnz Number of nonzero elements of the constraint matrix (length of arrays
 *      \p A_irn, \p A_jcn and \p A_x).
 *  \param H_irn Row subscripts of the lower triangular part (!) of the Hessian matrix
 *      (coordinate scheme).
 *  \param H_jcn Column subscripts of the lower triangular part (!) of the Hessian matrix
 *      (coordinate scheme).
 *  \param H_x Nonzero elements of the lower triangular part (!) of the Hessian matrix
 *      (coordinate scheme).
 *  \param H_nnz Number of nonzero elements of the lower triangular part (!) of the
 *      Hessian matrix (length of arrays \p H_irn, \p H_jcn and \p H_x).
 *  \param g Gradient vector of the QP.
 *  \param f Constant term of the objective function of the QP.
 *  \param xl Lower bounds on the variables of the QP. If this is a null pointer, then
 *      all lower bounds are set to -inf.
 *  \param xu Upper bounds on the variables of the QP. If this is a null pointer, then
 *      all upper bounds are set to +inf.
 *  \param al Lower bounds on the linear constraints of the QP. If this is a null pointer,
 *      then all lower bounds are set to -inf.
 *  \param au Upper bounds on the linear constraints of the QP. If this is a null pointer,
 *      then all upper bounds are set to +inf.
 *
 *  \return QPP_OK \n QPP_OUT_OF_MEMORY \n QPP_STATUS_ERROR \n QPP_INVALID_ARGUMENT
 *      \n QPP_NULL_ARGUMENT
 */
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
                           const qpp_real_t au[]);


/** \brief Presolves the given QP by applying several preprocessing techniques.
 *
 *  Which methods are used can be determined with the presolver options. The presolving
 *  loop is executed as long as a reduction can be found (and the maximum number of
 *  iterations is not reached).
 *
 *  \param data Valid presolver entity.
 *
 *  \return QPP_OK \n QPP_STATUS_ERROR \n QPP_NULL_ARGUMENT \n QPP_PRIMAL_INFEASIBLE
 *      \n QPP_DUAL_INFEASIBLE \n QPP_UNBOUNDED \n QPP_BOUNDS_INCOMPATIBLE
 *      \n QPP_BAD_SCALING
 */
qpp_return_value_t qppPresolve(qpp_data_t *const data);


/** \brief Returns the dimension of the presolved QP to the user.
 *
 *  Can also be directly accessed from \p data (after a call of \a qppPresolve()).
 *
 *  \param m Stores the number of linear constraints of the presolved QP.
 *  \param n Stores the number of variables of the presolved QP.
 *  \param A_nnz Stores the number of nonzero elements of the presolved constraint matrix.
 *  \param H_nnz Stores the number of nonzero elements of the lower triangular part (!)
 *      of the Hessian matrix.
 *
 *  \return QPP_OK \n QPP_NULL_ARGUMENT \n QPP_STATUS_ERROR
 */
qpp_return_value_t qppGetDimensionPresolvedQP(qpp_data_t *const data,
                                              qpp_int_t *m,
                                              qpp_int_t *n,
                                              qpp_int_t *A_nnz,
                                              qpp_int_t *H_nnz);


/** \brief Returns the presolved QP to the user.
 *
 *  Matrices are returned in coordinate scheme. Memory for the presolved QP is allocated
 *  by the user (the dimension can be retrieved by \a qppGetDimensionPresolvedQP()).
 *  Unlike the qppInit() method, the bounds \p xl, \p xu, \p al and \p au may not be
 *  null pointers, at least if the respective dimension is not zero. The given arrays for
 *  storing the matrices are allowed to be null pointers if they do not contain any
 *  nonzero elements. However, \p m, \p n, \p A_nnz and \p H_nnz must not be null pointers.
 *
 *  \param data Valid presolver entity.
 *  \param m Stores the number of constraints of the presolved QP.
 *  \param n Stores the number of variables of the presolved QP.
 *  \param A_irn Stores the row subscripts of the constraint matrix of the presolved QP.
 *  \param A_jcn Stores the column subscripts of the constraint matrix of the presolved QP.
 *  \param A_x Stores the nonzero elements of the constraint matrix of the presolved QP.
 *  \param A_nnz Stores the number of nonzero elements of the constraint matrix of the presolved QP.
 *  \param H_irn Stores the row subscripts of the lower triangular part (!) of the
 *      Hessian matrix of the presolved QP.
 *  \param H_jcn Stores the column subscripts of the lower triangular part (!) of the
 *      Hessian matrix of the presolved QP.
 *  \param H_x Stores the nonzero elements of the lower triangular part (!) of the
 *      Hessian matrix of the presolved QP.
 *  \param H_nnz Stores the number of nonzero elements of the lower triangular part (!)
 *      of the Hessian matrix of the presolved QP.
 *  \param g Stores the gradient vector of the presolved QP.
 *  \param f Stores the constant term of the objective function of the presolved QP.
 *  \param xl Stores the (tightest or medium) lower bounds on the variables of the presolved QP.
 *  \param xu Stores the (tightest or medium) upper bounds on the variables of the presolved QP.
 *  \param al Stores the lower bounds on the linear constraints of the presolved QP.
 *  \param au Stores the upper bounds on the linear constraints of the presolved QP.
 *  \param sort_type Determines if the returned matrices are sorted row-wise
 *      (QPP_MST_ROW_WISE) or column-wise (QPP_MST_COLUMN_WISE).
 *
 *  \return QPP_OK \n QPP_STATUS_ERROR \n QPP_INVALID_ARGUMENT \n QPP_NULL_ARGUMENT
 */
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
                                     const qpp_matrix_sort_type_t sort_type);


/** \brief Returns the presolved QP to the user.
 *
 *  Matrices are returned in coordinate scheme. Data is returned via call-by-reference.
 *  Memory is allocated by the presolver.
 *  The given pointer are allowed to be null pointers if the respective dimension, which
 *  can be retrieved by qppGetDimensionPresolvedQP(), is zero, otherwise they must be valid
 *  pointers. However, \p m, \p n, \p A_nnz and \p H_nnz must not be null pointers.
 *
 *  \param data Valid presolver entity.
 *  \param m Stores the number of constraints of the presolved QP.
 *  \param n Stores the number of variables of the presolved QP.
 *  \param A_irn Stores the row subscripts of the constraint matrix of the presolved QP.
 *  \param A_jcn Stores the column subscripts of the constraint matrix of the presolved QP.
 *  \param A_x Stores the nonzero elements of the constraint matrix of the presolved QP.
 *  \param A_nnz Stores the number of nonzero elements of the constraint matrix of the presolved QP.
 *  \param H_irn Stores the row subscripts of the lower triangular part (!) of the
 *      Hessian matrix of the presolved QP.
 *  \param H_jcn Stores the column subscripts of the lower triangular part (!) of the
 *      Hessian matrix of the presolved QP.
 *  \param H_x Stores the nonzero elements of the lower triangular part (!) of the
 *      Hessian matrix of the presolved QP.
 *  \param H_nnz Stores the number of nonzero elements of the lower triangular part (!)
 *      of the Hessian matrix of the presolved QP.
 *  \param g Stores the gradient vector of the presolved QP.
 *  \param f Stores the constant term of the objective function of the presolved QP.
 *  \param xl Stores the (tightest or medium) lower bounds on the variables of the presolved QP.
 *  \param xu Stores the (tightest or medium) upper bounds on the variables of the presolved QP.
 *  \param al Stores the lower bounds on the linear constraints of the presolved QP.
 *  \param au Stores the upper bounds on the linear constraints of the presolved QP.
 *  \param sort_type Determines if the returned matrices are sorted row-wise
 *      (QPP_MST_ROW_WISE) or column-wise (QPP_MST_COLUMN_WISE).
 *
 *  \return QPP_OK \n QPP_STATUS_ERROR \n QPP_INVALID_ARGUMENT \n QPP_NULL_ARGUMENT
 */
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
                                        const qpp_matrix_sort_type_t sort_type);


/** \brief Based on a primal-dual optimal solution of the presolved QP, a primal-dual
 *      optimal solution of the original QP is computed and returned to the user.
 *
 *  CAVEAT: The QP is not fully restored (only as far as it needs to be), e.g. the bounds
 *  on the linear / bound constraints will not be restored. The computed solution is
 *  stored in parameters x, y and z (and in wi and wj if they are also given).
 *
 *  \param data Valid presolver entity.
 *  \param x Vector of length n (#variables in original QP). The first n_ps
 *      (#variables in presolved QP) entries contain the primal solution of the presolved QP.
 *      Contains the primal solution of the original QP afterwards.
 *  \param y Vector of length m (#constraints in original QP). The first m_ps
 *      (#constraints in presolved QP) entries contain the dual solution of the presolved QP
 *      (w.r.t. linear constraints). Contains the dual solution of the original QP afterwards.
 *  \param z Vector of length n (#variables in original QP). The first n_ps
 *      (#variables in presolved QP) entries contain the dual solution of the presolved QP
 *      (w.r.t. bound constraints). Contains the dual solution of the original QP afterwards.
 *  \param wi Vector of length m (#constraints in original QP). The first
 *      m_ps (#constraints in presolved QP) entries contain the optimal working set w.r.t.
 *      linear constraints of the presolved QP. After postsolving, this parameter contains
 *      the optimal working set of the original QP. Set to NULL if the working set is not
 *      required (or an optimal working set of the presolved QP is not given).
 *  \param wj Vector of length n (#variables in original QP). The first
 *      n_ps (#variables in presolved QP) entries contain the optimal working set w.r.t.
 *      bound constraints of the presolved QP. After postsolving, this parameter contains
 *      the optimal working set of the original QP. Set to NULL if the working set is not
 *      required (or an optimal working set of the presolved QP is not given).
 *  \param use_copies If \p use_copies != 0, then some quantities of the presolved QP
 *      are copied such that multiple postsolving is possible (if desired).
 *
 *  \return QPP_OK \n QPP_STATUS_ERROR \n QPP_OUT_OF_MEMORY \n QPP_INVALID_ARGUMENT
 *      \n QPP_NULL_ARGUMENT
 */
qpp_return_value_t qppPostsolve(qpp_data_t *const data,
                                qpp_real_t *const x,
                                qpp_real_t *const y,
                                qpp_real_t *const z,
                                qpp_int_t *const wi,
                                qpp_int_t *const wj,
                                const qpp_int_t use_copies);


/** \brief Returns current options of the presolver entity.
 *
 *  \param data Presolver entity whose options will be returned. If this is a null
 *      pointer, then the default options will be returned.
 *
 *  \return Struct containing the current options of the given presolver entity
 *      (or default options).
 */
qpp_options_t qppGetOptions(const qpp_data_t *const data);


/** \brief Replaces current options of a presolver entity with new ones.
 *
 *  This function checks if all values are consistent and replaces inconsistent ones
 *  with default values. After the presolved QP has been retrieved by the user, e.g.
 *  via \a qppGetPresolvedQP() or \a qppGetPresolvedQPPtr(), the options can NOT be
 *  changed anymore.
 *
 *  \param data Presolver entity whose options will be changed.
 *  \param opt Struct containing new options. If NULL, then default options are used.
 *
 *  \return QPP_OK \n QPP_INVALID_ARGUMENT \n QPP_NULL_ARGUMENT \n QPP_INVALID_OPTIONS
 */
qpp_return_value_t qppSetOptions(qpp_data_t *const data,
                                 const qpp_options_t *const opt);


/** \brief Deallocates all memory maintained by the presolver entity.
 *
 *  \param pdata Presolver entity. After deallocating memory, *pdata is set to NULL.
 */
void qppFree(qpp_data_t **pdata);

#endif /* QPPRESOLVER_PRESOLVER_H */
