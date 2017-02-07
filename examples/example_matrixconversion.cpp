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
 *	\file examples/example_matrixconversion.cpp
 *	\author Christian Hoffmann
 *	\date 2016--2017
 *
 *	Example for the matrix conversion functions declared in MatrixConversions.hpp.
 */


#include <iostream>
#include <qpOASES/MatrixConversion.hpp>

USING_NAMESPACE_QPOASES
using namespace std;

void print_dense_mat( int_t nrow, int_t ncol, real_t * A )
{
	for( int_t ii = 0; ii < nrow; ii++ ){
		for( int jj = 0; jj < ncol; jj++ ){
			cout << A[ii*ncol+jj] << "\t ";
		}
		cout << endl;
	}
	cout << endl << endl;
}


void test_coo_to_csr()
{
	int_t kk;
	
	// initialize matrix A in COO format
	int_t const nrow = 4;
	int_t const ncol = 5;
	real_t const Ax[] = {11, 22, 33};
	int_t const Ai[] = {1, 2, 3};
	int_t const Aj[] = {3, 2, 1};
	int_t const nnz = sizeof(Ax)/sizeof(real_t);

	// init dense matrix D
	real_t * Dx = new real_t[nrow*ncol];

	// print matrix A
	for( kk = 0; kk < nnz; kk++ ) cout << Ax[kk] << "\t ";
	cout << endl;
	for( kk = 0; kk < nnz; kk++ ) cout << Ai[kk] << "\t ";
	cout << endl;
	for( kk = 0; kk < nnz; kk++ ) cout << Aj[kk] << "\t ";
	cout << endl << endl;

	// convert A to dense matrix D 
	coo_to_dense( nrow, ncol, nnz, Ax, Ai, Aj, Dx );
	print_dense_mat( nrow, ncol, Dx );
	
	// convert A to CSR matrix B
	real_t * Bx = new real_t[nnz];
	int_t * Bp = new int_t[nrow+1];
	int_t * Bj = new int_t[nnz];
	coo_to_csr( nrow, ncol, nnz, Ax, Ai, Aj, Bx, Bp, Bj );

	// print matrix B
	for( kk = 0; kk < nnz; kk++ ) cout << Bx[kk] << "\t ";
	cout << endl;
	for( kk = 0; kk < nrow+1; kk++ ) cout << Bp[kk] << "\t ";
	cout << endl;
	for( kk = 0; kk < nnz; kk++ ) cout << Bj[kk] << "\t ";
	cout << endl << endl;
	
	// convert CSR matrix B to dense matrix D
	csr_to_dense( nrow, ncol, nnz, Bx, Bp, Bj, Dx );
	print_dense_mat( nrow, ncol, Dx );
	
	delete[] Dx;

	delete[] Bx;
	delete[] Bp;
	delete[] Bj;
	
	return;
}


int main( )
{
	test_coo_to_csr();
	
// 	// load test matrix
// 	spmatrix * A = spmatrix_mmread( "./testmat.mm" );
// 	assert( A != NULL );
// 	
// 	// print test matrix to screen
//     qp_int rv = spmatrix_mmwrite( A, "stdout" );
// 	assert( rv == SPMATRIX_OK );	
// 	cout << "ncmax: " << A->ncmax << endl;
// 	cout << "nzmax: " << A->nzmax << endl;
// 	cout << "nrow: " << A->nrow << endl;
// 	cout << "ncol: " << A->ncol << endl;
//     for (int jj = 0; jj < A->ncol; jj++)
//     {
//         for (int ll = A->p[jj]; ll < A->p[jj+1]; ll++)
//         {
// 			cout << A->i[ll]+1 << " " << jj+1 << " " << A->x[ll] << endl;
//         }
//     }
// 
// 	SparseMatrix B = fromQoreMatrix( *A );
// 	B.print();
// 	
// 	// free matrix memory
// 	spmatrix_free( A );
	
	return 0;
}


/*
 *	end of file
 */
