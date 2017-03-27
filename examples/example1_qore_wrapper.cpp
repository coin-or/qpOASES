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
 *	\file examples/example1_qore_wrapper.cpp
 *	\example example1_qore_wrapper.cpp
 *	\brief Example of solving a very simple QP with qpOASES using the 
 *	QProblemQore class. The QP is mathematically identical to that in 
 *	'example1'.
 *	\author Christian Hoffmann
 *	\date 2016-2017
 */



#include <qpOASES.hpp>


int main( )
{
	USING_NAMESPACE_QPOASES

	/* Setup data of QP. */
	int_t const nv = 2; // number of variables
	int_t const nc = 1; // number of constraints
	
	real_t H[nv*nv] = { 1.0, 0.0, 0.0, 0.5 };
	real_t A[nc*nv] = { 1.0, 1.0 };
	real_t g[nv] = { 1.5, 1.0 };
	real_t lb[nv] = { 0.5, -2.0 };
	real_t ub[nv] = { 5.0, 2.0 };
	real_t lbA[nc] = { -1.0 };
	real_t ubA[nc] = { 2.0 };

	/* Setup data of second QP. */
	real_t g_new[nv] = { 1.0, 1.5 };
	real_t lb_new[nv] = { 0.0, -1.0 };
	real_t ub_new[nv] = { 5.0, -0.5 };
	real_t lbA_new[nc] = { -2.0 };
	real_t ubA_new[nc] = { 1.0 };

	
	/* Setting up QProblem object. */
	QProblemQore example( nv, nc );

	Options options;
	example.setOptions( options ); // has no effect, only for compatibility with QProblem

	/* Solve first QP. */
	int_t nWSR = 10; // is ignored, only for compatibility with QProblem
	example.init( H, g, A, lb, ub, lbA, ubA, nWSR );

	/* Get and print solution of QP. */
	real_t xOpt[nv];
	real_t yOpt[nv+nc];
	example.getPrimalSolution( xOpt );
	example.getDualSolution( yOpt );
	printf( "\nxOpt = [ %e, %e ]; yOpt = [ %e, %e, %e ]; objVal = %e\n\n", 
			xOpt[0], xOpt[1], yOpt[0], yOpt[1], yOpt[2], example.getObjVal() );
	
	/* Solve second QP. */
	nWSR = 10; // is ignored, only for compatibility with QProblem
	example.hotstart( g_new, lb_new, ub_new, lbA_new, ubA_new, nWSR );

	/* Get and print solution of second QP. */
	example.getPrimalSolution( xOpt );
	example.getDualSolution( yOpt );
	printf( "\nxOpt = [ %e, %e ]; yOpt = [ %e, %e, %e ]; objVal = %e\n\n", 
			xOpt[0], xOpt[1], yOpt[0], yOpt[1], yOpt[2], example.getObjVal() );

	example.printOptions();
	/*example.printProperties();*/

	/*getGlobalMessageHandler()->listAllMessages();*/

	
	return 0;
}


/*
 *	end of file
 */
