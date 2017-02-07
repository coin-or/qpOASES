/*
 *	This file is part of qpOASES.
 *
 *	qpOASES -- An Implementation of the Online Active Set Strategy.
 *	Copyright (C) 2007-2015 by Hans Joachim Ferreau, Andreas Potschka,
 *	Christian Kirches et al. All rights reserved.
 *
 *	qpOASES is free software; you can redistribute it and/or
 *	modify it under the terms of the GNU Lesser General Public
 *	License as published by the Free Software Foundation; either
 *	version 2.1 of the License, or (at your option) any later version.
 *
 *	qpOASES is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *	See the GNU Lesser General Public License for more details.
 *
 *	You should have received a copy of the GNU Lesser General Public
 *	License along with qpOASES; if not, write to the Free Software
 *	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */


/**
 *	\file include/qpOASES/MatrixConversion.hpp
 *	\author Christian Hoffmann
 *	\date 2016--2017
 *
 *	Function interfaces for converting between matrices in different formats.
 */



#ifndef QPOASES_MATRIXCONVERSION_HPP
#define QPOASES_MATRIXCONVERSION_HPP

#include <qpOASES/Matrices.hpp>

// FIXME temporary (chris)
extern "C" {
	#include <spmatrix.h>
}


BEGIN_NAMESPACE_QPOASES


/// Converts a sparse matrix in COO format (A) to CSR format (B).
/// \author Christian Hoffmann
/// \date 2016
/// 
/// Input row and column indices are not assumed to be ordered.
/// Duplicate entries are carried over to the CSR represention.
/// Output arrays Bp, Bj, and Bx must be preallocated.
/// 
/// Complexity: Linear, O(nnz(A) + max(nrow,ncol))
/// 
/// To convert A to CSC format, exchange cols and rows of A, that is, exchange
/// nrow <-> ncol and Ai <-> Aj
void coo_to_csr(
	int_t const nrow, ///< number of rows of A
	int_t const ncol, ///< number of columns of A
	int_t const nnz, ///< number of nonzeros of A
	real_t const * const Ax, ///< nnz-array of nonzeros of A
	int_t const * const Ai, ///< nnz-array of row indices of A
	int_t const * const Aj, ///< nnz-array of column indices of A
	real_t * const Bx, ///< nnz-array of nonzeros of B (output)
	int_t * const Bp, ///< nrow-array of row pointers of B (output)
	int_t * const Bj ///< nnz-array of column indices of B (output)
);


/// Converts a sparse matrix in COO format (A) to dense column-major format (B).
/// \author Christian Hoffmann
/// \date 2016
/// 
/// Input row and column indices are not assumed to be ordered.
/// Duplicate entries overwrite predecessors.
/// The output array Bx must be preallocated.
void coo_to_dense( 
	int_t const nrow, ///< number of rows of A
	int_t const ncol, ///< number of columns of A
	int_t const nnz, ///< number of nonzeros of A
	real_t const * const Ax, ///< nnz-array of nonzeros of A
	int_t const * const Ai, ///< nnz-array of row indices of A
	int_t const * const Aj, ///< nnz-array of column indices of A
	real_t * const Bx ///< nrow*ncol-array of nonzeros of B (output)
);


/// Converts a sparse matrix in CSR format (A) to dense column-major format (B).
/// \author Christian Hoffmann
/// \date 2016
/// 
/// The output array Bx must be preallocated.
void csr_to_dense( 
	int_t const nrow, ///< number of rows of A
	int_t const ncol, ///< number of columns of A
	int_t const nnz, ///< number of nonzeros of A
	real_t const * const Ax, ///< nnz-array of nonzeros of A
	int_t const * const Ap, ///< nrow-array of row pointers of A
	int_t const * const Aj, ///< nnz-array of column indices of A
	real_t * const Bx ///< nrow*ncol-array of nonzeros of B (output)
);


/**	\brief Convert sparse QORE matrix to sparse qpOASES matrix (deep copy).
 * 	\author Christian Hoffmann
 * 	\date 2016
 */
// SparseMatrix fromQoreMatrix(
// 	spmatrix const & /**< sparse matrix in QORE format */
// );


/**	\brief Convert sparse qpOASES matrix to sparse QORE matrix (deep copy).
 * 	\author Christian Hoffmann
 * 	\date 2016
 */
// spmatrix * fromQpoasesMatrix(
// 	Matrix const & /**< sparse matrix in qpOASES format */
// );


END_NAMESPACE_QPOASES

#endif	/* QPOASES_MATRIXCONVERSION_HPP */


/*
 *	end of file
 */
