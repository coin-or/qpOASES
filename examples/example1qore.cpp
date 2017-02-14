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
 *	\file examples/example1qore.cpp
 *	\author Christian Hoffmann
 *	\date 2016-2017
 *
 *	Variant of example1 (without hotstarts) using the QProblemQore class 
 *	instead of QProblem.
 */



#include <qpOASES.hpp>


/** Example for qpOASES main function using the QProblemQore class. */
int main( )
{
	/* Setup data of QP. */
	REFER_NAMESPACE_QPOASES real_t H[2*2] = { 1.0, 0.0, 0.0, 0.5 };
	REFER_NAMESPACE_QPOASES real_t A[1*2] = { 1.0, 1.0 };
	REFER_NAMESPACE_QPOASES real_t g[2] = { 1.5, 1.0 };
	REFER_NAMESPACE_QPOASES real_t lb[2] = { 0.5, -2.0 };
	REFER_NAMESPACE_QPOASES real_t ub[2] = { 5.0, 2.0 };
	REFER_NAMESPACE_QPOASES real_t lbA[1] = { -1.0 };
	REFER_NAMESPACE_QPOASES real_t ubA[1] = { 2.0 };

	/* Setting up QProblem object. */
	REFER_NAMESPACE_QPOASES QProblemQore example( 2, 1 );

	REFER_NAMESPACE_QPOASES Options options;
	example.setOptions( options );

	/* Solve QP. */
	REFER_NAMESPACE_QPOASES int_t nWSR = 10;
	example.init( H, g, A, lb, ub, lbA, ubA, nWSR );

	/* Get and print solution of QP. */
	REFER_NAMESPACE_QPOASES real_t xOpt[2] = {1, 2}; // TODO remove init once getXSolution implemented
	REFER_NAMESPACE_QPOASES real_t yOpt[2+1] = {3, 4, 5}; // TODO remove init once getXSolution implemented
	example.getPrimalSolution( xOpt );
	example.getDualSolution( yOpt );
	printf( "\nxOpt = [ %e, %e ];  yOpt = [ %e, %e, %e ];  objVal = %e\n\n",  xOpt[0], xOpt[1], yOpt[0], yOpt[1], yOpt[2], example.getObjVal() );
	
	example.printOptions();

	return 0;
}


/*
 *	end of file
 */
