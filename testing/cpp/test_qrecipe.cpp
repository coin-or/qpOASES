/*
 *	This file is part of qpOASES.
 *
 *	qpOASES -- An Implementation of the Online Active Set Strategy.
 *	Copyright (C) 2007-2014 by Hans Joachim Ferreau, Andreas Potschka,
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
 *	\file testing/cpp/test_qrecipe.cpp
 *	\author Andreas Potschka
 *	\version 3.0
 *	\date 2007-2014
 *
 *	QRECIPE example from the CUTEr test set with sparse matrices.
 */



#include <qpOASES.hpp>
#include <qpOASES/UnitTesting.hpp>

#include "test_qrecipe_data.hpp"



int main( )
{
	USING_NAMESPACE_QPOASES

	long i;
	int nWSR;
	real_t tic, toc;
	real_t errP=0.0, errD=0.0;
	real_t *x1 = new real_t[180];
	real_t *y1 = new real_t[271];
	real_t *x2 = new real_t[180];
	real_t *y2 = new real_t[271];

	/* create sparse matrices */
	SymSparseMat *H = new SymSparseMat(180, 180, H_ir, H_jc, H_val);
	SparseMatrix *A = new SparseMatrix(91, 180, A_ir, A_jc, A_val);

	H->createDiagInfo();

	real_t* H_full = H->full();
	real_t* A_full = A->full();

	SymDenseMat *Hd = new SymDenseMat(180,180,180,H_full);
	DenseMatrix *Ad = new DenseMatrix(91,180,180,A_full);

	/* solve with dense matrices */
	nWSR = 1000;
	QProblem qrecipeD(180, 91);
	tic = getCPUtime();
	qrecipeD.init(Hd, g, Ad, lb, ub, lbA, ubA, nWSR, 0);
	toc = getCPUtime();
	qrecipeD.getPrimalSolution(x1);
	qrecipeD.getDualSolution(y1);

	fprintf(stdFile, "Solved dense problem in %d iterations, %.3f seconds.\n", nWSR, toc-tic);

	/* Compute KKT tolerances */
	real_t stat, feas, cmpl;

	getKKTResidual(	180,91,
					H_full,g,A_full,lb,ub,lbA,ubA,
					x1,y1,
					stat,feas,cmpl
					);
	printf( "stat = %e\nfeas = %e\ncmpl = %e\n\n", stat,feas,cmpl );


	/* solve with sparse matrices */
	nWSR = 1000;
	QProblem qrecipeS(180, 91);
	tic = getCPUtime();
	qrecipeS.init(H, g, A, lb, ub, lbA, ubA, nWSR, 0);
	toc = getCPUtime();
	qrecipeS.getPrimalSolution(x2);
	qrecipeS.getDualSolution(y2);

	fprintf(stdFile, "Solved sparse problem in %d iterations, %.3f seconds.\n", nWSR, toc-tic);
	
	/* Compute KKT tolerances */
	real_t stat2, feas2, cmpl2;
	getKKTResidual(	180,91,
					H_full,g,A_full,lb,ub,lbA,ubA,
					x2,y2,
					stat2,feas2,cmpl2
					);
	printf( "stat = %e\nfeas = %e\ncmpl = %e\n\n", stat2,feas2,cmpl2 );

	/* check distance of solutions */
	for (i = 0; i < 180; i++)
		if (getAbs(x1[i] - x2[i]) > errP)
			errP = getAbs(x1[i] - x2[i]);
	fprintf(stdFile, "Primal error: %9.2e\n", errP);

	for (i = 0; i < 271; i++)
		if (getAbs(y1[i] - y2[i]) > errD)
			errD = getAbs(y1[i] - y2[i]);
	fprintf(stdFile, "Dual error: %9.2e (might not be unique)\n", errD);

	delete H;
	delete A;
	delete[] H_full;
	delete[] A_full;
	delete Hd;
	delete Ad;

	delete[] y2;
	delete[] x2;
	delete[] y1;
	delete[] x1;

	QPOASES_TEST_FOR_TOL( stat,1e-14 );
	QPOASES_TEST_FOR_TOL( feas,1e-14 );
	QPOASES_TEST_FOR_TOL( cmpl,1e-13 );
	
	QPOASES_TEST_FOR_TOL( stat2,1e-14 );
	QPOASES_TEST_FOR_TOL( feas2,1e-14 );
	QPOASES_TEST_FOR_TOL( cmpl2,1e-13 );

	QPOASES_TEST_FOR_TOL( errP,1e-13 );

	return TEST_PASSED;
}


/*
 *	end of file
 */
