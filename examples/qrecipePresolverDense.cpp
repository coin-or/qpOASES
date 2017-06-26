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
 *	\file examples/qrecipePresolverDense.cpp
 *	\author Dominik Cebulla, Andreas Potschka
 *	\version 3.2
 *	\date 2007-2015
 *
 *	QRECIPE example from the CUTEr test set with dense matrices and preprocessing.
 */


#include <qpOASES.hpp>

#include "qrecipe_data.hpp"


int main()
{
    USING_NAMESPACE_QPOASES

    long i;
	int_t nWSR;
	int_t nV = 180;			/* Number of variables. */
	int_t nC = 91;			/* Number of constraints. */
	int_t nVPrs = 0;		/* Number of variables of presolved QP. */
	int_t nCPrs = 0;		/* Number of constraints of presolved QP. */
	real_t err;
	real_t tic;
	real_t toc;
	real_t *x = new real_t[nV];
	real_t *y = new real_t[nC+nV];
	real_t *xPrs = new real_t[nV];
	real_t *yPrs = new real_t[nC+nV];
	DenseMatrix* AdPrs = 0;
	SymDenseMat* HdPrs = 0;

	/* Sparse matrices are only used to create dense matrices out of them. */
	SymSparseMat *H = new SymSparseMat(nV, nV, H_ir, H_jc, H_val);
	SparseMatrix *A = new SparseMatrix(nC, nV, A_ir, A_jc, A_val);
	H->createDiagInfo();


	/* Get dense matrices in row-major format. */
	real_t* H_full = H->full();
	real_t* A_full = A->full();

	SymDenseMat *Hd = new SymDenseMat(nV, nV, nV, H_full);
	DenseMatrix *Ad = new DenseMatrix(nC, nV, nV, A_full);


	/* Solve with dense matrices, but without preprocessing. */
	nWSR = 1000;
	QProblem qrecipe(nV, nC);
	tic = getCPUtime();
	qrecipe.init(Hd, g, Ad, lb, ub, lbA, ubA, nWSR, 0);
	toc = getCPUtime();
	qrecipe.getPrimalSolution(x);
	qrecipe.getDualSolution(y);

	fprintf(stdFile, "Solved dense problem in %d iterations, %.3f seconds.\n",
            (int)nWSR, toc-tic);


    /* Solve with dense matrices and also with preprocessing. */
    nWSR = 1000;
    Presolver prs(nV, nC);  /* Number of nonzeros in A / H can optionally be passed to constructor. */

    /* Optionally, user options can be passed to the presolver: */
    PresolverOptions popt;
    popt.disableDuplicateColumnsMethod();
    popt.setStabilityTol(1e-4);

    prs.setOptions(popt);

    tic = getCPUtime();

    /*  If you want to use preprocessing with dense matrices, then use one of the
        two following methods: */
	/* 1) */
    prs.presolve(&nVPrs, &nCPrs, Hd, &HdPrs, g, Ad, &AdPrs, lb, ub, lbA, ubA);

    /* 2) */
    /*prs.presolve(&nVPrs, &nCPrs, H_full, g, A_full, lb, ub, lbA, ubA);
	HdPrs = new SymDenseMat(nVPrs, nVPrs, nVPrs, H_full);
	AdPrs = new DenseMatrix(nCPrs, nVPrs, nVPrs, A_full);*/

	toc = getCPUtime();

	/* Print information of presolved QP (number of constraints and variables). */
	fprintf(stdFile, "\nNumber of variables in original QP:    %d\n", nV);
	fprintf(stdFile, "Number of constraints in original QP:  %d\n", nC);
	fprintf(stdFile, "Number of variables in presolved QP:   %d\n", nVPrs);
	fprintf(stdFile, "Number of constraints in presolved QP: %d\n\n", nCPrs);


	/* Solve presolved QP. */
	QProblem qrecipePrs(nVPrs, nCPrs);

	tic += getCPUtime();
	qrecipePrs.init(HdPrs, g, AdPrs, lb, ub, lbA, ubA, nWSR, 0);
	toc += getCPUtime();

	qrecipePrs.getPrimalSolution(xPrs);
	qrecipePrs.getDualSolution(yPrs);


	/* Retrieve optimal primal-dual solution of original QP via postprocessing. */
	prs.postsolve(xPrs, yPrs);

	fprintf(stdFile, "Solved presolved dense problem in %d iterations, %.3f seconds "
			"(Time for preprocessing included).\n\n", (int) nWSR, toc-tic);

	/* Check distance of solutions. */
	err = 0.0;
	for (i = 0; i < nV; i++)
		if (getAbs(x[i] - xPrs[i]) > err)
			err = getAbs(x[i] - xPrs[i]);
	fprintf(stdFile, "Primal error: %9.2e\n", err);
	err = 0.0;
	for (i = 0; i < nV+nC; i++)
		if (getAbs(y[i] - yPrs[i]) > err)
			err = getAbs(y[i] - yPrs[i]);
	fprintf(stdFile, "Dual error: %9.2e  (might not be unique)\n", err);


    delete H;
    delete A;
    delete[] H_full;
    delete[] A_full;
    delete Hd;
    delete Ad;
    delete HdPrs;
    delete AdPrs;

    delete[] x;
    delete[] y;
    delete[] xPrs;
    delete[] yPrs;

    return 0;
}

