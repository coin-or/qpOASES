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
 *	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
 *
 */


/**
 *	\file src/MatrixConversion.cpp
 *	\author Christian Hoffmann
 *	\date 2016--2017
 *
 *	Implementations of function for converting between matrices in different formats.
 */

#include <assert.h>
#include <qpOASES/MatrixConversion.hpp>

BEGIN_NAMESPACE_QPOASES


void coo_to_dense( 
	int_t const nrow, int_t const ncol, int_t const nnz,
	real_t const * const Ax, int_t const * const Ai, int_t const * const Aj,
	real_t * const Bx
)
{
	assert( ncol > 0 );
	assert( nnz > 0 ); // TODO allow case nnz==0
	assert( Ai != 0 );
	assert( Aj != 0 );
	assert( Ax != 0 );
	assert( Bx != 0 );
	
	for( int_t kk = 0; kk < nrow*ncol; kk++ ) Bx[kk] = 0;
	
	for( int_t kk = 0; kk < nnz; kk++ ){
		Bx[(Aj[kk]-1)+(Ai[kk]-1)*ncol] = Ax[kk];
	}
}


void csr_to_dense( 
	int_t const nrow, int_t const ncol, int_t const nnz,
	real_t const * const Ax, int_t const * const Ap, int_t const * const Aj,
	real_t * const Bx
)
{
	assert( nrow > 0 );
	assert( ncol > 0 );
	assert( nnz > 0 ); // TODO allow case nnz==0
	assert( Ax != 0 );
	assert( Ap != 0 );
	assert( Aj != 0 );
	assert( Bx != 0 );
	
	for( int_t kk = 0; kk < nrow*ncol; kk++ ) Bx[kk] = 0;
	
	for( int_t ii = 0; ii < nrow; ii++ ){ // loop over rows
		for( int_t kk = Ap[ii]; kk < Ap[ii+1]; kk++ ){ // loop over 
			Bx[ii+1+(Aj[kk]-1)*ncol] = Ax[kk];
		}
	}

}


void coo_to_csr(
	int_t const nrow, int_t const ncol, int_t const nnz, 
	real_t const * const Ax, int_t const * const Ai, int_t const * const Aj,
	real_t * const Bx, int_t * const Bp, int_t * const Bj )
{
	assert( nrow > 0 );
	assert( ncol > 0 );
	assert( nnz > 0 );  // TODO allow case nnz==0
	assert( Ai != 0 );
	assert( Aj != 0 );
	assert( Ax != 0 );
	assert( Bp != 0 );
	assert( Bj != 0 );
	assert( Bx != 0 );
	
	int_t kk;
	
	// compute number of non-zero entries per row of A 
	for( kk = 0; kk < nrow; kk++ ) Bp[kk] = 0;
	for( kk = 0; kk < nnz; kk++ ) Bp[Ai[kk]]++;

	// cumsum the nnz per row to get Bp[]
	int_t cumsum = 0;
	for( kk = 0; kk < nrow; kk++ ){
		int_t temp = Bp[kk];
		Bp[kk] = cumsum;
		cumsum += temp;
	}
	Bp[nrow] = nnz; 

	// write Aj and Ax into Bj and Bx, respectively
	for( kk = 0; kk < nnz; kk++ ){
		int_t row = Ai[kk];
		int_t dest = Bp[row];

		Bj[dest] = Aj[kk];
		Bx[dest] = Ax[kk];

		Bp[row]++;
	}

	int_t last = 0;
	for( kk = 0; kk <= nrow; kk++){
		int_t temp = Bp[kk];
		Bp[kk] = last;
		last = temp;
	}
}


// SparseMatrix fromQoreMatrix( spmatrix const & mat )
// {
// 	return SparseMatrix( mat.nrow, mat.ncol, mat.i, mat.p, mat.x );
// }
// 
// 
// spmatrix * fromQpoasesMatrix(
// 	Matrix const & A /**< sparse matrix in qpOASES format */
// )
// {
// 	long kk;
// 	
// 	int_t const nr = A.getRowNum();
// 	int_t const nc = A.getColNum();
// 	int_t nnz = 0;
// 	returnValue rv;
// 
// 	// step 1: get matrix in triplet / coodinate list format
// 	
// 	int_t * const rows = new int_t[nr]; // array of row indices
// 	int_t * const cols = new int_t[nc]; // array of col indices
// 	for( kk = 1; kk <= nr; kk++ )
// 		rows[kk] = kk;
// 	for( kk = 1; kk <= nc; kk++ )
// 		cols[kk] = kk;
// 	
// 	// get number of nonzeros
// 	rv = A.getSparseSubmatrix( nr, rows, nc, cols, 0, 0, nnz, 0, 0, 0 );
// 	assert( rv == SUCCESSFUL_RETURN );
// 	
// 	// allocate mem for triplet format
// 	int_t * const ridx = new int_t[nnz];
// 	int_t * const cidx = new int_t[nnz];
// 	real_t * const val = new real_t[nnz];
// 	
// 	// retrieve matrix in triplet format
// 	rv = A.getSparseSubmatrix( nr, rows, nc, cols, 0, 0, nnz, ridx, cidx, val );
// 	assert( rv == SUCCESSFUL_RETURN );
// 	
// 	// step 2: convert triplet format to column-compressed format
// 	
// 	
// 	// step 3: create output matrix
// // 	spmatrix * B = spmatrix_alloc( nnz, nc );
// 	spmatrix * B = 0;
// 	
// // 	B->nrow = nr;
// // 	B->ncol = nc;
// // 
// // 	returnValue rv = A.getSparseSubmatrix( nr, rows, nc, cols, 0, 0, nnz, 0, 0, B->x );
// // 	assert( rv == SUCCESSFULL_RETURN );
// // 	
// // 	dupl->ir = new sparse_int_t[nnz];
// // 	dupl->jc = new sparse_int_t[nCols+1];
// // 	dupl->val = new real_t[nnz];
// // 
// // 	for (kk = 0; kk <= nc; i++) B->p[kk] = A.jci[i];
// // 	for (kk = 0; kk < nnz; i++) B->i[kk] = A.ir[kk];
// // 	for (kk = 0; kk < nnz; i++) B->x[kk] = A.val[kk];
// 
// 	delete[] val;
// 	delete[] cidx;
// 	delete[] ridx;
// 	
// 	delete[] cols;
// 	delete[] rows;
// 	
// 	return B;
// 
// }


END_NAMESPACE_QPOASES

/*
 *	end of file
 */
