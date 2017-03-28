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
 *	\file include/qpOASES/LapackBlasDeclarations.hpp
 *	\author Andreas Potschka, Hans Joachim Ferreau, Christian Kirches
 *	\version 3.2
 *	\date 2009-2017
 *
 *  Declarations for external LAPACK/BLAS functions.
 */



#ifndef QPOASES_LAPACKBLASDECLARATIONS_HPP
#define QPOASES_LAPACKBLASDECLARATIONS_HPP


extern "C"
{
	/** Performs one of the matrix-matrix operation in double precision. */
	void dgemm_(	const char*, const char*, const unsigned long*, const unsigned long*, const unsigned long*,
					const double*, const double*, const unsigned long*, const double*, const unsigned long*,
					const double*, double*, const unsigned long* );
	/** Performs one of the matrix-matrix operation in single precision. */
	void sgemm_(	const char*, const char*, const unsigned long*, const unsigned long*, const unsigned long*,
					const float*, const float*, const unsigned long*, const float*, const unsigned long*,
					const float*, float*, const unsigned long* );

	/** Performs a symmetric rank 1 operation in double precision. */
	void dsyr_(		const char *, const unsigned long *, const double *, const double *,
					const unsigned long *, double *, const unsigned long *);
	/** Performs a symmetric rank 1 operation in single precision. */
	void ssyr_(		const char *, const unsigned long *, const float *, const float *,
					const unsigned long *, float *, const unsigned long *);

	/** Performs a symmetric rank 2 operation in double precision. */
	void dsyr2_(	const char *, const unsigned long *, const double *, const double *,
					const unsigned long *, const double *, const unsigned long *, double *, const unsigned long *);
	/** Performs a symmetric rank 2 operation in single precision. */
	void ssyr2_(	const char *, const unsigned long *, const float *, const float *,
					const unsigned long *, const float *, const unsigned long *, float *, const unsigned long *);

	/** Calculates the Cholesky factorization of a real symmetric positive definite matrix in double precision. */
	void dpotrf_(	const char *, const unsigned long *, double *, const unsigned long *, long * );
	/** Calculates the Cholesky factorization of a real symmetric positive definite matrix in single precision. */
	void spotrf_(	const char *, const unsigned long *, float *, const unsigned long *, long * );
}

#endif	/* QPOASES_LAPACKBLASDECLARATIONS_HPP */
