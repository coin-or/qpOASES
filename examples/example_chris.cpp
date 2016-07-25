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
 *	\file examples/example_chris.cpp
 *	\author Christian Hoffmann
 *	\version 3.2
 *	\date 2016
 *
 *	Example for testing Chris' code during development
 */



#include <iostream>
#include <qpOASES.hpp>
extern "C" { 
	#include <spmatrix.h>
}
#include <qpOASES/ChrisTemp.hpp>

USING_NAMESPACE_QPOASES
using namespace std;

int main( )
{

	// load test matrix
	spmatrix * A = spmatrix_mmread( "./testmat.mm" );
	assert( A != NULL );
	
	// print test matrix to screen
    qp_int rv = spmatrix_mmwrite( A, "stdout" );
	assert( rv == SPMATRIX_OK );
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

	SparseMatrix B = fromQoreMatrix( *A );
	B.print();
	
// 	SparseMatrix qpoasismat();
// 	SparseMatrix qpoasismat( fromQoreMatrix( *A ) );
	
// 	qp_int nr = 6;
// 	qp_int nc = 6;
// 	qp_int nnz = 14;

// 	spmatrix * A = spmatrix_alloc( nnz, nc );
	
// 	qp_int p[] = {1, 2, 4, 7, 10, 12, 15};
// 	qp_int i[] = {1, 1, 2, 1, 2, 3, 1, 3, 4, 3, 4, 3, 5, 6};
// 	qp_int x[] = {11, 12, 22, 13, 23, 33, 14, 34, 44, 35, 45, 36, 56, 66};
// 
// 	for( unsigned int j=0; j<sizeof(p); j++ ) {
// 		A->p[j] = p[j];
// 	}
// 	for( unsigned int j=0; j<sizeof(i); j++ ) {
// 		A->i[j] = i[j];
// 	}
// 	for( unsigned int j=0; j<sizeof(x); j++ ) {
// 		A->x[j] = x[j];
// 	}

// 	A->p[] = {1, 2, 4, 7, 10, 12, 15};
// 	A->i = {1, 1, 2, 1, 2, 3, 1, 3, 4, 3, 4, 3, 5, 6};
// 	A->x = {11, 12, 22, 13, 23, 33, 14, 34, 44, 35, 45, 36, 56, 66};
	
	// free matrix memory
	spmatrix_free( A );
	
// 	real_t A[2*2] = { 1.0, 0.0, 0.0, 0.5 };

	return 0;
}

/** Example for qpOASES main function using the QProblem class. */
// int main( )
// {
// 	USING_NAMESPACE_QPOASES
// 
// 	/* Setup data of first QP. */
// 	real_t H[2*2] = { 1.0, 0.0, 0.0, 0.5 };
// 	real_t A[1*2] = { 1.0, 1.0 };
// 	real_t g[2] = { 1.5, 1.0 };
// 	real_t lb[2] = { 0.5, -2.0 };
// 	real_t ub[2] = { 5.0, 2.0 };
// 	real_t lbA[1] = { -1.0 };
// 	real_t ubA[1] = { 2.0 };
// 
// 	/* Setup data of second QP. */
// 	real_t g_new[2] = { 1.0, 1.5 };
// 	real_t lb_new[2] = { 0.0, -1.0 };
// 	real_t ub_new[2] = { 5.0, -0.5 };
// 	real_t lbA_new[1] = { -2.0 };
// 	real_t ubA_new[1] = { 1.0 };
// 
// 
// 	/* Setting up QProblem object. */
// 	QProblem example( 2,1 );
// 
// 	Options options;
// 	example.setOptions( options );
// 
// 	/* Solve first QP. */
// 	int_t nWSR = 10;
// 	example.init( H,g,A,lb,ub,lbA,ubA, nWSR );
// 
// 	/* Get and print solution of first QP. */
// 	real_t xOpt[2];
// 	real_t yOpt[2+1];
// 	example.getPrimalSolution( xOpt );
// 	example.getDualSolution( yOpt );
// 	printf( "\nxOpt = [ %e, %e ];  yOpt = [ %e, %e, %e ];  objVal = %e\n\n", 
// 			xOpt[0],xOpt[1],yOpt[0],yOpt[1],yOpt[2],example.getObjVal() );
// 	
// 	/* Solve second QP. */
// 	nWSR = 10;
// 	example.hotstart( g_new,lb_new,ub_new,lbA_new,ubA_new, nWSR );
// 
// 	/* Get and print solution of second QP. */
// 	example.getPrimalSolution( xOpt );
// 	example.getDualSolution( yOpt );
// 	printf( "\nxOpt = [ %e, %e ];  yOpt = [ %e, %e, %e ];  objVal = %e\n\n", 
// 			xOpt[0],xOpt[1],yOpt[0],yOpt[1],yOpt[2],example.getObjVal() );
// 
// 	example.printOptions();
// 	/*example.printProperties();*/
// 
// 	/*getGlobalMessageHandler()->listAllMessages();*/
// 
// 	return 0;
// }


/*
 *	end of file
 */
