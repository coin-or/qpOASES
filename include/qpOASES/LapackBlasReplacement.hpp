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
 *	\file include/qpOASES/LapackBlasReplacement.hpp
 *	\author Andreas Potschka, Hans Joachim Ferreau, Christian Kirches
 *	\version 3.2
 *	\date 2009-2017
 *
 *  Declarations for external LAPACK/BLAS functions.
 */



#ifndef QPOASES_LAPACKBLASREPLACEMENT_HPP
#define QPOASES_LAPACKBLASREPLACEMENT_HPP


extern "C"
{
	/** Performs one of the matrix-matrix operation in double precision. */
	void dgemm_(	const char*, const char*, const la_uint_t*, const la_uint_t*, const la_uint_t*,
					const double*, const double*, const la_uint_t*, const double*, const la_uint_t*,
					const double*, double*, const la_uint_t* );
	/** Performs one of the matrix-matrix operation in single precision. */
	void sgemm_(	const char*, const char*, const la_uint_t*, const la_uint_t*, const la_uint_t*,
					const float*, const float*, const la_uint_t*, const float*, const la_uint_t*,
					const float*, float*, const la_uint_t* );

	/** Performs a symmetric rank 1 operation in double precision. */
	void dsyr_(		const char*, const la_uint_t*, const double*, const double*,
					const la_uint_t*, double*, const la_uint_t* );
	/** Performs a symmetric rank 1 operation in single precision. */
	void ssyr_(		const char*, const la_uint_t*, const float*, const float*,
					const la_uint_t*, float*, const la_uint_t* );

	/** Performs a symmetric rank 2 operation in double precision. */
	void dsyr2_(	const char*, const la_uint_t*, const double*, const double*,
					const la_uint_t*, const double*, const la_uint_t*, double*, const la_uint_t*);
	/** Performs a symmetric rank 2 operation in single precision. */
	void ssyr2_(	const char*, const la_uint_t*, const float*, const float*,
					const la_uint_t*, const float*, const la_uint_t*, float*, const la_uint_t*);

	/** Calculates the Cholesky factorization of a real symmetric positive definite matrix in double precision. */
	void dpotrf_(	const char*, const la_uint_t*, double*, const la_uint_t*, la_int_t* );
	/** Calculates the Cholesky factorization of a real symmetric positive definite matrix in single precision. */
	void spotrf_(	const char*, const la_uint_t*, float*, const la_uint_t*, la_int_t* );
}

#endif	/* QPOASES_LAPACKBLASREPLACEMENT_HPP */
