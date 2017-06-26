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
 *	\file include/qpPresolver/ecrmatrix.h
 *	\author Dominik Cebulla
 *	\version 1.0 Beta
 *	\date 2017
 *
 *  Declaration of a sparse matrix class stored in (extended) Curtis-Reid format
 *  (basically based on singly linked lists).
 */


#ifndef QPPRESOLVER_ECRMATRIX_H
#define QPPRESOLVER_ECRMATRIX_H

#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>

#include <qpPresolver/constants.h>
#include <qpPresolver/types.h>


/*  ======================================================================================
    Definition of the data type for the sparse matrix format.
    ====================================================================================*/

/** \brief Struct containing necessary data for the realization of a matrix based on the
 *      Curtis-Reid scheme.
 *
 *  The Curtis-Reid scheme is described in: <tt> Duff, I.S., Erisman, A.M., Reid, J.K.:
 *  Direct Methods for Sparse Matrices. Oxford Univ. Press, Oxford, England, 1986 </tt>.
 *  We extend the scheme by providing direct access to the row and column subscripts
 *  via \p irn and \p jcn (hence the term Extended Curtis-Reid scheme).
 */
typedef struct ecrmatrix
{
    qpp_int_t *irowst;      /**< Subscripts of the first elements in every row. */
    qpp_int_t *jcolst;      /**< Subscripts of the first elements in every column. */
    qpp_int_t *irn;         /**< Row subscripts. */
    qpp_int_t *jcn;         /**< Column subscripts. */
    qpp_int_t *linkrw;      /**< Links to the next elements in a row. */
    qpp_int_t *linkcl;      /**< Links to the next elements in a column. */
    qpp_real_t *x;          /**< Nonzero elements of the matrix. */
    qpp_int_t nnz;          /**< Number of nonzero elements currently stored. */
    qpp_int_t nrow;         /**< Number of rows of the matrix. */
    qpp_int_t ncol;         /**< Number of columns of the matrix. */
    qpp_int_t nnzmax;       /**< Maximum allowed number of nonzero elements, can be reallocated. */
    qpp_int_t nrowmax;      /**< Maximum number of rows of the matrix (unchangeable!). */
    qpp_int_t ncolmax;      /**< Maximum number of columns of the matrix (unchangeable!). */
} qpp_ecrmatrix_t;



/*  ======================================================================================
    Interface functions.
    ====================================================================================*/

/** \brief Allocates memory for a new ecrmatrix of maximum dimension \p nrowmax x \p ncolmax.
 *
 *  Creates a \p nrowmax x \p ncolmax matrix in extended Curtis-Reid format with a
 *  maximum number of \p nrowmax rows, \p ncolmax columns and \p nnzmax nonzero elements.
 *  The actual dimension, i.e. \a nrow and \a ncol is set to 0.
 *
 *  \param nrowmax Maximum number of rows the matrix can store. Can not be changed afterwards.
 *      Matrix is allowed to actually contain less rows.
 *  \param ncolmax Maximum number of columns the matrix can store. Can not be changed afterwards.
 *      Matrix is allowed to actually contain less columns.
 *  \param nnzmax Maximum number of nonzero elements the matrix can store.
 *      Can be changed afterwards by reallocating.
 *
 *  \return Valid pointer to a new ecrmatrix on success, otherwise a null pointer.
 */
qpp_ecrmatrix_t *ecrAlloc(const qpp_int_t nrowmax,
                          const qpp_int_t ncolmax,
                          const qpp_int_t nnzmax);


/** \brief Changes the maximum number of nonzero elements the given matrix can store.
 *
 *  Reallocates memory for matrix \p A. If elements will be deleted, then the links
 *  of \p A will be rearranged such that it still will be a valid matrix. \n
 *
 *  CAVEAT: Reducing the size of \p A is possible, but since the order of the elements
 *  depends only on the links between elements, the resultant matrix might not be
 *  as expected.
 *
 *  \param A Pointer to a valid \p ecrmatrix whose maximum number of nonzero elements will
 *      be changed.
 *  \param nnzmax New maximum number of nonzero elements. May be smaller or larger
 *      than the original maximum number of nonzero elements of \p A.
 *
 *  \return QPP_OK \n QPP_FATAL_ERROR \n QPP_OUT_OF_MEMORY
 *      \n QPP_INVALID_ARGUMENT \n QPP_NULL_ARGUMENT
 */
qpp_return_value_t ecrRealloc(qpp_ecrmatrix_t *A,
                              const qpp_int_t nnzmax);


/** \brief Creates a deep copy of the given matrix.
 *
 *  \param A Valid pointer to an \p ecrmatrix.
 *
 *  \return Pointer to the copy of \p A, or null2 pointer on failure.
 */
qpp_ecrmatrix_t *ecrCopy(const qpp_ecrmatrix_t *A);


/** \brief Deallocates memory allocated by the given matrix.
 *
 *  \param Pointer to an \p ecrmatrix whose memory will be deallocated.
 */
void ecrFree(qpp_ecrmatrix_t *A);


/** \brief Creates a new ecrmatrix from a matrix given in coordinate scheme.
 *
 *  Function takes a matrix in coordinate format (with \p nrow rows,
 *  \p ncol columns, \p nnz nonzero elements, row subscripts \p irn, column subscripts
 *  \p jcn and nonzero elements \p x) and creates an \p ecrmatrix from this data. \n
 *
 *  CAVEAT: \p A has to be allocated first via ecrAlloc().
 *  Matrix subscripts must be zero-based.
 *
 *  \param A Valid pointer to an (allocated) \p ecrmatrix. Data is initialized with the
 *      matrix given in coordinate scheme.
 *  \param nrow Number of rows of the matrix given in coordinate scheme. Must not be
 *      greater than the maximum number of rows of \p A.
 *  \param ncol Number of columns of the matrix given in coordinate scheme. Must not be
 *      greater than the maximum number of columns of \p A.
 *  \param nnz Number of nonzero elements of the matrix given in coordinate scheme,
 *      i.e. length of the arrays \p irn, \p jcn and \p x. If \p nnz is greater than the
 *      maximum allowed number of nonzero elements in \p A, then memory is reallocated.
 *  \param irn Row subscripts of the matrix given in coordinate scheme.
 *  \param jcn Column subscripts of the matrix given in coordinate scheme.
 *  \param x Nonzero elements of the matrix given in coordinate scheme.
 *  \param sort_type Determines how the matrix will be sorted; either row-wise
 *      (QPP_MST_ROW_WISE) or column-wise (QPP_MST_COLUMN_WISE).
 *  \param only_lt If parameter equals QPP_BT_TRUE and the given matrix is a square matrix
 *      (i.e. \p nrow = \p ncol), then only the lower triangular part of the matrix will
 *      be stored.
 *
 *  \return QPP_OK \n QPP_FATAL_ERROR \n QPP_OUT_OF_MEMORY \n
 *      QPP_INVALID_ARGUMENT \n QPP_NULL_ARGUMENT
 */
qpp_return_value_t ecrCreateFromCoordinate(qpp_ecrmatrix_t *A,
                                           const qpp_int_t nrow,
                                           const qpp_int_t ncol,
                                           const qpp_int_t nnz,
                                           const qpp_int_t irn[],
                                           const qpp_int_t jcn[],
                                           const qpp_real_t x[],
                                           const qpp_matrix_sort_type_t sort_type,
                                           const qpp_bool_type_t only_lt);


/** \brief Computes y = a * M * x + y with M = A or M = A^T.
 *
 *  \param transp If \p transp equals 't' or 'T', then y = a * A^T * x + y is computed,
 *      otherwise y = a * A * x + y is computed.
 *  \param A Valid pointer to an \p ecrmatrix which is used in the matrix-vector product.
 *  \param a Scalar multiplier.
 *  \param x Vector which is used in the matrix-vector product.
 *  \param y Vector which is used as an additive term and which stores the computed result.
 */
void ecrAxpy(const char transp,
             const qpp_ecrmatrix_t *A,
             const qpp_real_t a,
             const qpp_real_t x[],
             qpp_real_t y[]);


/** \brief Checks if a given matrix is an upper or lower triangular matrix.
 *
 *  \param A Valid pointer to an ecrmatrix.
 *  \param type If \p type < 0, then it is checked whether \p A is a lower triangular
 *      matrix. Otherwise it is checked whether \p A is an upper triangular matrix.
 *
 *  \return QPP_BT_TRUE if the test is positive, QPP_BT_FALSE otherwise.
 */
qpp_bool_type_t ecrIsTriangular(const qpp_ecrmatrix_t *A,
                                const qpp_int_t type);


/** \brief Scales the given matrix with a scaling algorithm proposed by Curtis and Reid.
 *
 *  The scaling algorithm can be found in: <tt> Curtis, A.R. and Reid, J.K.
 *  "On the automatic scaling of matrices for Gaussian elimination." In: IMA Journal
 *  of Applied Mathematics 10.1 (1972), pp. 118-124 </tt>.
 *
 *  \param A Valid pointer to an \p ecrmatrix. This matrix will be scaled.
 *  \param R Vector in which the row scaling factors are stored. Memory must be allocated
 *      by the user. If the scaling procedure fails, this vector represents an identity matrix.
 *  \param C Vector in which the column scaling factors are stored. Memory must be allocated
 *      by the user. If the scaling procedure fails, this vector represents an identity matrix.
 *  \param max_iter The procedure is an iterative algorithm. The maximum allowed number of
 *      iterations is given by \p max_iter. If no convergence is achieved within this
 *      maximum number of iterations, then matrix \p A is not scaled and the scaling
 *      factors stored in \p R and \p C are set to 1.
 *
 *  \return QPP_OK \n QPP_BAD_SCALING \n QPP_OUT_OF_MEMORY
 *      \n QPP_INVALID_ARGUMENT \n QPP_NULL_ARGUMENT
 */
qpp_return_value_t ecrScaleCR(qpp_ecrmatrix_t *A,
                              qpp_real_t R[],
                              qpp_real_t C[],
                              const qpp_int_t max_iter);

#endif /* QPPRESOLVER_ECRMATRIX_H */
