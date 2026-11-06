/*
 *	This file is part of qpOASES.
 *
 *	qpOASES -- An Implementation of the Online Active Set Strategy.
 *	Copyright (C) 2007-2017 by Hans Joachim Ferreau, Andreas Potschka,
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
 *	\file testing/cpp/test_flipper.cpp
 *	\author qpOASES contributors
 *	\version 3.2
 *	\date 2026
 *
 *	Unit test for copying Flipper objects (copy constructor and assignment
 *	must take over the dimensions of the source object, see issue #103).
 */



#include <qpOASES.hpp>
#include <qpOASES/UnitTesting.hpp>


USING_NAMESPACE_QPOASES


/** Returns BT_TRUE if f holds exactly the R, Q and T entries of a 3x2 flipper set up below. */
static BooleanType holdsData( const Flipper& f, const real_t* R, const real_t* Q, const real_t* T )
{
	real_t R2[9], Q2[9], T2[4];
	int_t i;

	for( i=0; i<9; ++i ) { R2[i] = -1.0; Q2[i] = -1.0; }
	for( i=0; i<4; ++i ) T2[i] = -1.0;

	f.get( 0,R2, 0,Q2,T2 );

	for( i=0; i<9; ++i )
		if ( ( isEqual( R2[i],R[i] ) == BT_FALSE ) || ( isEqual( Q2[i],Q[i] ) == BT_FALSE ) )
			return BT_FALSE;
	for( i=0; i<4; ++i )
		if ( isEqual( T2[i],T[i] ) == BT_FALSE )
			return BT_FALSE;

	return BT_TRUE;
}


int main( )
{
	real_t R[9], Q[9], T[4];
	int_t i;

	for( i=0; i<9; ++i ) { R[i] = 1.0 + i; Q[i] = 10.0 + i; }
	for( i=0; i<4; ++i ) T[i] = 100.0 + i;

	/* nV = 3, nC = 2: R and Q are 3x3, T is 2x2. */
	Flipper source( 3,2 );
	source.set( 0,R, 0,Q,T );

	/* Assigning to a flipper set up with smaller dimensions. */
	Flipper assigned( 1,1 );
	assigned = source;
	QPOASES_TEST_FOR_TRUE( holdsData( assigned, R,Q,T ) == BT_TRUE )

	/* Copy constructor. */
	Flipper constructed( source );
	QPOASES_TEST_FOR_TRUE( holdsData( constructed, R,Q,T ) == BT_TRUE )

	return TEST_PASSED;
}


/*
 *	end of file
 */
