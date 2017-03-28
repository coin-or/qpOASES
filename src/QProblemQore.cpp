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
*	\file src/QProblemQore.cpp
*	\author Christian Hoffmann
*	\date 2016-2017
*
*	Implementation of the QProblemQore class, a QProblem-compatible wrapper 
*	around the QORE QP solver.
*/

#include <assert.h> 
#include <stdlib.h>
#include <string.h>
#include <qpOASES/QProblemQore.hpp>

extern "C" {
	#include "matrixconversion.h"
}


BEGIN_NAMESPACE_QPOASES

/// \todo Internally calls QORE's QPNew(), which dynamically allocates memory 
/// and may thus fail unpredictably.  The proper way of dealing with a failed 
/// allocation in the constructor would be to throw an exception, but 
/// exceptions are currently not used in qpOASES. We currently rely on 
/// assertions instead. A failed allocation thus leads to an immediate program 
/// termination!
QProblemQore::QProblemQore ( int_t _nV, int_t _nC ) : pproblem(0), nV(_nV), nC(_nC)
{
	// check preconditions 
	assert( _nV >= 1 );
	assert( _nC >= 0 );
	
	qp_int const rv = QPNew( &pproblem, _nV, _nC, 0, 0 );
	assert( rv == QPSOLVER_OK );
	assert( pproblem != 0 );
	
	/// \todo check postconditions
}


QProblemQore::~QProblemQore( ) 
{
	QPFree( &pproblem );
	assert( pproblem == 0 );
}


returnValue 
QProblemQore::setOptions( const Options & _options )
{
	return SUCCESSFUL_RETURN;
}


returnValue 
QProblemQore::printOptions( ) const
{
	return SUCCESSFUL_RETURN;
}	


returnValue 
QProblemQore::init(	const real_t* const _H, const real_t* const _g,
					const real_t* const _lb, const real_t* const _ub,
					int_t &, real_t* const, const real_t * const,
					const real_t * const, const Bounds * const,
					const real_t * const _R
)
{
	// check preconditions
	assert( _H != 0 );
	assert( _g != 0 );
	assert( _lb != 0 );
	assert( _ub != 0 );
	
	// check invariants
	assert( pproblem != 0 );
	assert( nV > 0 );
	assert( nC >= 0 );

	// convert Hessian matrix to csr format
	qp_int * Hri = 0;
	qp_int * Hcp = 0;
	real_t * Hnz = 0;
	matrix_dense_to_csc( nV, nV, _H, &Hcp, &Hri, &Hnz );

	// init working variables
	qp_int rv = QPSOLVER_OK;
	
	// init QP
	rv = QPSetData( pproblem, nV, nC, 0, 0, 0, Hcp, Hri, Hnz );
	assert( rv == QPSOLVER_OK );

	// free memory for H in csc format -- QPSetData() created deep copies internally
	free( Hri );
	free( Hcp );
	free( Hnz );

	// solve QP 
	rv = QPOptimize( pproblem, _lb, _ub, _g, 0, 0 );
	assert( rv == QPSOLVER_OK );

	/// \todo check postconditions
	
	return SUCCESSFUL_RETURN;
}


returnValue 
QProblemQore::hotstart(	const real_t* const _g_new,
						const real_t* const _lb_new, const real_t* const _ub_new,
						int_t & nWSR, real_t * const, const Bounds * const
) {
	/// \todo check preconditions
	
	qp_int const rv = QPOptimize( pproblem, _lb_new, _ub_new, _g_new, 0, 0 );
	assert( rv == QPSOLVER_OK );
	
	/// \todo check postconditions
	return SUCCESSFUL_RETURN;
}


returnValue 
QProblemQore::init(	const real_t* const _H, const real_t* const _g,
					const real_t* const _A,
					const real_t* const _lb, const real_t* const _ub,
					const real_t* const _lbA, const real_t* const _ubA,
					int_t&, real_t* const, const real_t* const, 
					const real_t* const, const Bounds* const, 
					const Constraints* const, const real_t* const
){
	// check preconditions
	assert( _H != 0 );
	assert( _g != 0 );
	assert( _A != 0 );
	assert( _lb != 0 );
	assert( _ub != 0 );
	assert( _lbA != 0 );
	assert( _ubA != 0 );

	// check invariants
	assert( pproblem != 0 );
	assert( nV > 0 );
	assert( nC >= 0 );

	// convert Hessian matrix to csr format
	qp_int * Hcp = 0;
	qp_int * Hri = 0;
	real_t * Hnz = 0;
	matrix_dense_to_csc( nV, nV, _H, &Hcp, &Hri, &Hnz );

    // convert constraint matrix to csc format and transpose (as required by QORE)
	// TODO with matrix_dense_to_csc_trans() we get a valgrind hit, yet none 
	// with matrix_dense_to_csc()... error in the former?
	qp_int * Atcp = 0;
	qp_int * Atri = 0;
	real_t * Atnz = 0;
	matrix_dense_to_csc_trans( nC, nV, _A, &Atcp, &Atri, &Atnz );

	// create (nV+nC) arrays of bounds to variables AND constraints: [x ; A'*x]
	real_t * const lb = (real_t *) malloc( sizeof(real_t)*(nV+nC) ); // lower bounds
	assert( lb != 0 );
	memcpy( lb, _lb, nV*sizeof(real_t) );
	memcpy( lb+nV, _lbA, nC*sizeof(real_t) );
	real_t * const ub = (real_t *) malloc( sizeof(real_t)*(nV+nC) ); // upper bounds
	assert( ub != 0 );
	memcpy( ub, _ub, nV*sizeof(real_t) );
	memcpy( ub+nV, _ubA, nC*sizeof(real_t) );

	// init working variables
	qp_int rv = QPSOLVER_OK;
	
	// pass problem data to QP solver instance
	rv = QPSetData( pproblem, nV, nC, Atcp, Atri, Atnz, Hcp, Hri, Hnz );
	assert( rv == QPSOLVER_OK );

	// free memory for H and A' in csc format -- QPSetData() created deep copies internally
	free( Hri );
	free( Hcp );
	free( Hnz );
	free( Atri );
	free( Atcp );
	free( Atnz );
	
	// solve QP 
	rv = QPOptimize( pproblem, lb, ub, _g, 0, 0 );
	assert( rv == QPSOLVER_OK );

	free(lb);
	free(ub);
	
	/// \todo check postconditions
	
	return SUCCESSFUL_RETURN;
}	


returnValue 
QProblemQore::hotstart(	const real_t* const _g_new,
			const real_t* const _lb_new, const real_t* const _ub_new,
			const real_t* const _lbA_new, const real_t* const _ubA_new,
			int_t&, real_t* const, const Bounds* const, const Constraints* const 
) 
{
	// check preconditions
	assert( _g_new != 0 );
	assert( _lb_new != 0 );
	assert( _ub_new != 0 );
	assert( _lbA_new != 0 );
	assert( _ubA_new != 0 );

	// check invariants
	assert( pproblem != 0 );
	assert( nV > 0 );
	assert( nC >= 0 );

	// create (nV+nC) arrays of bounds to variables AND constraints: [x ; A'*x]
	real_t * const lb = (real_t *) malloc( sizeof(real_t)*(nV+nC) ); // lower bounds
	assert( lb != 0 );
	memcpy( lb, _lb_new, nV*sizeof(real_t) );
	memcpy( lb+nV, _lbA_new, nC*sizeof(real_t) );
	real_t * const ub = (real_t *) malloc( sizeof(real_t)*(nV+nC) ); // upper bounds
	assert( ub != 0 );
	memcpy( ub, _ub_new, nV*sizeof(real_t) );
	memcpy( ub+nV, _ubA_new, nC*sizeof(real_t) );

	qp_int const rv = QPOptimize( pproblem, lb, ub, _g_new, 0, 0 );
	assert( rv == QPSOLVER_OK );
	
	/// \todo check postconditions
	return SUCCESSFUL_RETURN;
}


returnValue 
QProblemQore::getPrimalSolution( real_t* const xOpt ) 
const
{
	/// \todo check preconditions

	qp_int rv = QPGetDblVector( pproblem, "primalsol", xOpt );
	assert( rv == QPSOLVER_OK );

	/// \todo check postconditions

	return SUCCESSFUL_RETURN;
}	


returnValue 
QProblemQore::getDualSolution( real_t* const yOpt ) 
const
{
	/// \todo check preconditions

	qp_int rv = QPGetDblVector( pproblem, "dualsol", yOpt );
	assert( rv == QPSOLVER_OK );

	/// \todo check postconditions

	return SUCCESSFUL_RETURN;
}	


real_t 
QProblemQore::getObjVal( ) 
const
{
	return NAN;
}	


END_NAMESPACE_QPOASES


/*
*	end of file
*/
