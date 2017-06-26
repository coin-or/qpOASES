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
 *	\file include/qpPresolver/utility.h
 *	\author Dominik Cebulla
 *	\version 1.0 Beta
 *	\date 2017
 *
 *  Declaration of useful utility functions and some I/O functions for matrices
 *  stored in MatrixMarket format, see http://math.nist.gov/MatrixMarket/mmio-c.html.
 *
 */


#ifndef QPPRESOLVER_UTILITY_H
#define QPPRESOLVER_UTILITY_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "../../extern/mmio/include/mmio.h"

#include <qpPresolver/types.h>
#include <qpPresolver/constants.h>
#include <qpPresolver/return_values.h>


/** \brief Prints error message on stderr and redirects error code.
 *
 *  \param error Error code / return value.
 *  \param filename Name of the file in which the error occurred.
 *  \param line_number Line number in which the error occurred.
 */
qpp_return_value_t throwError(const qpp_return_value_t error,
                              const char* filename,
                              const int line_number);

/** \brief Macro for easy error printing and redirecting. */
#define QPP_THROWERROR(ERROR) throwError((ERROR), __FILE__, __LINE__)


/** \brief Sets all elements of an integer array to a specific value.
 *
 *  \param x Integer array.
 *  \param n Length of array \p x.
 *  \param value Every element in array \p x is set to this \p value.
 */
void qppSetArrayi(qpp_int_t x[],
                  qpp_int_t n,
                  const qpp_int_t value);


/** \brief Sets all elements of an array containing real values to a specific value.
 *
 *  \param x Real array.
 *  \param n Length of array \p x.
 *  \param value Every element in array \p x is set to this \p value.
 */
void qppSetArrayr(qpp_real_t x[],
                  qpp_int_t n,
                  const qpp_real_t value);


/** \brief Swaps two integer elements.
 *
 *  \param a Integer value 1.
 *  \param b Integer value 2.
 */
void qppSwapi(qpp_int_t* a,
              qpp_int_t* b);


/** \brief Swaps two real elements.
 *
 *  \param a Real value 1.
 *  \param b Real value 2.
 */
void qppSwapr(qpp_real_t* a,
              qpp_real_t* b);


/** \brief Prints an integer array on stdout.
 *
 *  \param x Integer array.
 *  \param n Length of array \p x.
 */
void qppPrintArrayi(const qpp_int_t x[],
                    const qpp_int_t n);


/** \brief Prints an array containing real values on stdout.
 *
 *  \param x Real array.
 *  \param n Length of array \p x.
 */
void qppPrintArrayr(const qpp_real_t x[],
                    const qpp_int_t n);


/** \brief Returns the maximum of two integer values.
 *
 *  \param a Integer value 1.
 *  \param b Integer value 2.
 *
 *  \return max(a,b)
 */
qpp_int_t qppMaxi(const qpp_int_t a,
                  const qpp_int_t b);


/** \brief Returns the minimum of two integer values.
 *
 *  \param a Integer value 1.
 *  \param b Integer value 2.
 *
 *  \return min(a,b)
 */
qpp_int_t qppMini(const qpp_int_t a,
                  const qpp_int_t b);


/** \brief Returns the maximum of two real values.
 *
 *  \param a Real value 1.
 *  \param b Real value 2.
 *
 *  \return max(a,b)
 */
qpp_real_t qppMaxr(const qpp_real_t a,
                   const qpp_real_t b);


/** \brief Returns the minimum of two real values.
 *
 *  \param a Real value 1.
 *  \param b Real value 2.
 *
 *  \return min(a,b)
 */
qpp_real_t qppMinr(const qpp_real_t a,
                   const qpp_real_t b);


/** \brief Checks if two real values are equal up to a certain relative error.
 *
 *  \param x Real value 1.
 *  \param y Real value 2.
 *  \param tol Tolerance specifying the maximum allowed relative error.
 *
 *  \return QPP_BT_TRUE if \p a and \p b are equal (up to a certain tolerance),
 *      QPP_BT_FALSE otherwise.
 */
qpp_bool_type_t qppIsEqual(const qpp_real_t x,
                           const qpp_real_t y,
                           const qpp_real_t tol);


/** \brief Checks if a real value is greater than another real value up to a
 *      certain relative error.
 *
 *  \param x Real value 1 (it is checked whether this is the greater value).
 *  \param y Real value 2 (it is checked whether this is the smaller value).
 *  \param tol Tolerance specifying the maximum allowed relative error.
 *
 *  \return QPP_BT_TRUE if \p a is greater than \p b (up to a certain tolerance),
 *      QPP_BT_FALSE otherwise.
 */
qpp_bool_type_t qppIsGreater(const qpp_real_t x,
                             const qpp_real_t y,
                             const qpp_real_t tol);


/** \brief Computes a dot product of a full (dense) vector with a sparse vector.
 *
 *  \param indices Subscripts of the sparse vector.
 *  \param values Nonzero elements of the sparse vector.
 *  \param length Number of nonzero elements in the sparse vector (i.e. length
 *      of \p indices and \p values).
 *  \param x Full (dense) vector.
 *
 *  \return Value of the dot product.
 */
qpp_real_t qppSparseDot(const qpp_int_t indices[],
                        const qpp_real_t values[],
                        const qpp_int_t length,
                        const qpp_real_t x[]);


/** \brief Rearranges a given sparse vector into two parts: The first part contains only
 *      negative elements and the second part contains only positive elements.
 *
 *  After the call of this function, the first \p *n_neg_entries in the sparse vector are
 *  all negative and the following \p *n_pos_entries are all positive.
 *
 *  \param indices Subscripts of the sparse vector.
 *  \param values Nonzero elements of the sparse vector.
 *  \param length Length of arrays \p indices and \p values.
 *  \param n_neg_entries Will contain the number of negative entries in the sparse vector.
 *  \param n_pos_entries Will contain the number of positive entries in the sparse vector.
 *
 *  \return QPP_OK \n QPP_NULL_ARGUMENT
 */
qpp_return_value_t qppSortNegPos(qpp_int_t indices[],
                                 qpp_real_t values[],
                                 const qpp_int_t length,
                                 qpp_int_t *n_neg_entries,
                                 qpp_int_t *n_pos_entries);


/** \name MatrixMarket-I/O-functions. */
/**@{*/

/** \brief Reads a dense vector from a MatrixMarket file or from stdin.
 *
 *  The vector is returned via call-by-reference. The vector is allowed to contain
 *  infinite entries. This function allocates memory for the vector and the user has
 *  to deallocate this memory again.
 *
 *  \param filename Name of the MatrixMarket file which shall be read.
 *  \param x The vector that has been read from the MatrixMarket file is stored in \p x.
 *  \param n Will contain the length of vector \p x.
 *
 *  \return QPP_OK \n QPP_UNABLE_TO_OPEN_FILE \n QPP_UNABLE_TO_READ_FILE \n
 *      QPP_OUT_OF_MEMORY \n QPP_INVALID_ARGUMENT \n QPP_NULL_ARGUMENT
 */
qpp_return_value_t qppMMReadDenseVector(const char *filename,
                                        qpp_real_t *x[],
                                        qpp_int_t *n);


/** \brief Writes a dense vector into a MatrixMarket file or into stdout.
 *
 *  \param filename Name of the file in which the vector shall be written (in MatrixMarket
 *      format). File descriptor stdout is also possible.
 *  \param x Array representing the dense vector.
 *  \param n Length of the vector.
 *
 *  \return QPP_OK \n QPP_UNABLE_TO_OPEN_FILE \n QPP_UNABLE_TO_WRITE_FILE
 */
qpp_return_value_t qppMMWriteDenseVector(const char *filename,
                                         const qpp_real_t x[],
                                         const qpp_int_t n);


/** \brief Reads a sparse matrix in coordinate scheme from a MatrixMarket file
 *      or from stdin.
 *
 *  The sparse matrix is returned via call-by-reference. Every nonzero element must be
 *  finite. This function allocates memory for the matrix and the user has to
 *  deallocate this memory again.
 *
 *  \param filename Name of the MatrixMarket file which shall be read.
 *  \param m Will contain the number of rows of the sparse matrix.
 *  \param n Will contain the number of columns of the sparse matrix.
 *  \param nnz Will contain the number of nonzero elements of the sparse matrix.
 *  \param irn Will contain the row subscripts of the sparse matrix.
 *  \param jcn Will contain the column subscripts of the sparse matrix.
 *  \param x Will contain the nonzero elements of the sparse matrix.
 *  \param store_ut If the matrix is (skew-) symmetric, then the MatrixMarket file
 *      only contains the lower triangular part of the matrix. If \p store_ut != 0,
 *      then we will also store the upper triangular part of the matrix, otherwise
 *      we only store the lower triangular part. If the matrix is not (skew-) symmetric,
 *      then this parameter has no effect.
 *
 *  \return QPP_OK \n QPP_UNABLE_TO_OPEN_FILE \n QPP_UNABLE_TO_READ_FILE \n
 *      QPP_OUT_OF_MEMORY \n QPP_INVALID_ARGUMENT \n QPP_NULL_ARGUMENT
 */
qpp_return_value_t qppMMReadSparseCoordinate(const char *filename,
                                             qpp_int_t *m,
                                             qpp_int_t *n,
                                             qpp_int_t *nnz,
                                             qpp_int_t *irn[],
                                             qpp_int_t *jcn[],
                                             qpp_real_t *x[],
                                             const qpp_int_t store_ut);


/** \brief Writes a sparse matrix in coordinate scheme into a MatrixMarket file
 *      or into stdout.
 *
 *  \param filename Name of the file in which the matrix shall be written (in
 *      MatrixMarket format). File descriptor stdout is also possible.
 *  \param m Number of rows of the sparse matrix.
 *  \param n Number of columns of the sparse matrix.
 *  \param nnz Number of nonzero elements of the sparse matrix (number of entries which
 *      will be written into the file).
 *  \param irn Row subscripts of the sparse matrix.
 *  \param jcn Column subscripts of the sparse matrix.
 *  \param x Nonzero elements of the sparse matrix.
 *  \param is_sym If \p is_sym < 0, then we treat the matrix as skew-symmetric and if
 *      \p is_sym > 0, then we treat the matrix as symmetric. This information is written
 *      into the file and only the lower triangular part of the matrix will be written into
 *      the file. Otherwise, the whole matrix will be written into the file. \n
 *      CAVEAT: It is not checked whether the sparse matrix is really (skew-) symmetric!
 *
 *  \return QPP_OK \n QPP_UNABLE_TO_OPEN_FILE \n QPP_UNABLE_TO_WRITE_FILE \n
 *      QPP_INVALID_ARGUMENT \n QPP_NULL_ARGUMENT
 */
qpp_return_value_t qppMMWriteSparseCoordinate(const char *filename,
                                              const qpp_int_t m,
                                              const qpp_int_t n,
                                              const qpp_int_t nnz,
                                              const qpp_int_t irn[],
                                              const qpp_int_t jcn[],
                                              const qpp_real_t x[],
                                              const qpp_int_t is_sym);
/**@}*/


#endif /* QPPRESOLVER_UTILITY_H */
