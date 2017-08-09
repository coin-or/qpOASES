/*
 *	This file is part of qpPresolver.
 *
 *	qpPresolver -- An implementation of presolving (= preprocessing) techniques
 *  for Quadratic Programming.
 *	Copyright (C) 2017 by Dominik Cebulla et al. All rights reserved.
 *
 *	qpPresolver is free software; you can redistribute it and/or
 *	modify it under the terms of the GNU Lesser General Public
 *	License as published by the Free Software Foundation; either
 *	version 2.1 of the License, or (at your option) any later version.
 *
 *	qpPresolver is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *	See the GNU Lesser General Public License for more details.
 *
 *	You should have received a copy of the GNU Lesser General Public
 *	License along with qpPresolver; if not, write to the Free Software
 *	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */


 /**
 *	\file include/qpPresolver/cholmod_config.h
 *	\author Dominik Cebulla
 *	\version 1.0 Beta
 *	\date 2017
 *
 *  Definition of several macros to use functions from cholmod, SPQR and umfpack (all part of
 *  SuiteSparse) with 32 or 64 bit integers.
 *
 */

#include <cholmod.h>
#include <umfpack.h>


#ifdef QPP_USE_INT64
    #define CM_START cholmod_l_start
    #define CM_ALLOCATE_SPARSE cholmod_l_allocate_sparse
    #define CM_ANALYZE cholmod_l_analyze
    #define CM_FACTORIZE cholmod_l_factorize
    #define CM_ALLOCATE_DENSE cholmod_l_allocate_dense
    #define CM_SOLVE cholmod_l_solve
    #define CM_FREE_SPARSE cholmod_l_free_sparse
    #define CM_FREE_DENSE cholmod_l_free_dense
    #define CM_FINISH cholmod_l_finish
    #define CM_WRITE_SPARSE cholmod_l_write_sparse
    #define SPQR_BACKSLASH SuiteSparseQR_C_backslash_default
    #define SPQR_BACKSLASH_SPARSE SuiteSparseQR_C_backslash_sparse
#else
    #define CM_START cholmod_start
    #define CM_ALLOCATE_SPARSE cholmod_allocate_sparse
    #define CM_ANALYZE cholmod_analyze
    #define CM_FACTORIZE cholmod_factorize
    #define CM_ALLOCATE_DENSE cholmod_allocate_dense
    #define CM_SOLVE cholmod_solve
    #define CM_FREE_SPARSE cholmod_free_sparse
    #define CM_FREE_DENSE cholmod_free_dense
    #define CM_FINISH cholmod_finish
    #define CM_WRITE_SPARSE cholmod_write_sparse
    #define SPQR_BACKSLASH SuiteSparseQR_C_backslash_default
    #define SPQR_BACKSLASH_SPARSE SuiteSparseQR_C_backslash_sparse
#endif


#ifdef QPP_USE_INT64
	#define UMFPACK_DEFAULTS umfpack_dl_defaults
	#define UMFPACK_SYMBOLIC umfpack_dl_symbolic
	#define UMFPACK_NUMERIC umfpack_dl_numeric
	#define UMFPACK_FREE_SYMBOLIC umfpack_dl_free_symbolic
	#define UMFPACK_FREE_NUMERIC umfpack_dl_free_numeric
	#define UMFPACK_GET_LUNZ umfpack_dl_get_lunz
	#define UMFPACK_GET_NUMERIC umfpack_dl_get_numeric
	#define UMFPACK_WSOLVE umfpack_dl_wsolve
	#define UMFPACK_SOLVE umfpack_dl_solve
#else
	#define UMFPACK_DEFAULTS umfpack_di_defaults
	#define UMFPACK_SYMBOLIC umfpack_di_symbolic
	#define UMFPACK_NUMERIC umfpack_di_numeric
	#define UMFPACK_FREE_SYMBOLIC umfpack_di_free_symbolic
	#define UMFPACK_FREE_NUMERIC umfpack_di_free_numeric
	#define UMFPACK_GET_LUNZ umfpack_di_get_lunz
	#define UMFPACK_GET_NUMERIC umfpack_di_get_numeric
	#define UMFPACK_WSOLVE umfpack_di_wsolve
	#define UMFPACK_SOLVE umfpack_di_solve
#endif
