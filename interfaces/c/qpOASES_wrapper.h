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
 *	\file interfaces/c/qpOASES_wrapper.h
 *	\author Hans Joachim Ferreau
 *	\version 3.1
 *	\date 2014-2015
 *
 *	Interface that enables to call qpOASES from plain C.
 *
 */


#ifndef QPOASES_WRAPPER_H
#define QPOASES_WRAPPER_H


#ifndef QPOASES_TYPES_HPP

	#ifdef __USE_SINGLE_PRECISION__
	typedef float real_t;
	#else
	typedef double real_t;
	#endif /* __USE_SINGLE_PRECISION__ */

	/* dummy definitions, not used when calling from C */
	#define QProblemBClass int
	#define OptionsClass int
	#define returnValue int

	/* HessianType */
	#define HST_ZERO             0
	#define HST_IDENTITY         1
	#define HST_POSDEF           2
	#define HST_POSDEF_NULLSPACE 3
	#define HST_SEMIDEF          4
	#define HST_INDEF            5
	#define HST_UNKNOWN	         6

	/* SubjectToStatus */
	#define ST_LOWER            -1
	#define ST_INACTIVE          0
	#define ST_UPPER             1
	#define ST_INFEASIBLE_LOWER  2
	#define ST_INFEASIBLE_UPPER  3
	#define ST_UNDEFINED         4

	/* PrintLevel */
	#define PL_DEBUG_ITER       -2
	#define PL_TABULAR          -1
	#define PL_NONE              0
	#define PL_LOW               1
	#define PL_MEDIUM            2
	#define PL_HIGH              3

#else 

	#define QProblemBClass QProblemB
	#define OptionsClass REFER_NAMESPACE_QPOASES Options 

	/* only declare when compiling C++ library */
	static QProblem*  globalQProblemObject  = 0;
	static QProblemB* globalQProblemBObject = 0;
	static SQProblem* globalSQProblemObject = 0;
	static Options globalOptionsObject;

#endif /* QPOASES_TYPES_HPP */



/**
 *	\brief Manages all user-specified options for solving QPs.
 *
 *	This struct manages all user-specified options used for solving
 *	quadratic programs.
 *
 *	\author Hans Joachim Ferreau
 *	\version 3.1
 *	\date 2014
 */
typedef struct
{
	int printLevel;						/**< Print level. */

	int enableRamping;					/**< Specifies whether ramping shall be enabled or not. */
	int enableFarBounds;				/**< Specifies whether far bounds shall be used or not. */
	int enableFlippingBounds;			/**< Specifies whether flipping bounds shall be used or not. */
	int enableRegularisation;			/**< Specifies whether Hessian matrix shall be regularised in case semi-definiteness is detected. */
	int enableFullLITests;				/**< Specifies whether condition-hardened LI test shall be used or not. */
	int enableNZCTests;					/**< Specifies whether nonzero curvature tests shall be used. */
	int enableDriftCorrection;			/**< Specifies the frequency of drift corrections (0 = off). */
	int enableCholeskyRefactorisation;	/**< Specifies the frequency of full refactorisation of proj. Hessian (otherwise updates). */
	int enableEqualities;				/**< Specifies whether equalities shall be always treated as active constraints. */

	real_t terminationTolerance;		/**< Termination tolerance. */
	real_t boundTolerance;				/**< Lower/upper (constraints') bound tolerance (an inequality constraint whose lower and upper bounds differ by less is regarded to be an equality constraint). */
	real_t boundRelaxation;				/**< Offset for relaxing (constraints') bounds at beginning of an initial homotopy. It is also as initial value for far bounds. */
	real_t epsNum;						/**< Numerator tolerance for ratio tests. */
	real_t epsDen;						/**< Denominator tolerance for ratio tests. */
	real_t maxPrimalJump;				/**< Maximum allowed jump in primal variables in nonzero curvature tests. */
	real_t maxDualJump;					/**< Maximum allowed jump in dual variables in linear independence tests. */

	real_t initialRamping;				/**< Start value for Ramping Strategy. */
	real_t finalRamping;				/**< Final value for Ramping Strategy. */
	real_t initialFarBounds;			/**< Initial size of Far Bounds. */
	real_t growFarBounds;				/**< Factor to grow Far Bounds. */
	int initialStatusBounds;			/**< Initial status of bounds at first iteration. */
	real_t epsFlipping;					/**< Tolerance of squared Cholesky diagonal factor which triggers flipping bound. */
	int numRegularisationSteps;			/**< Maximum number of successive regularisation steps. */
	real_t epsRegularisation;			/**< Scaling factor of identity matrix used for Hessian regularisation. */
	int numRefinementSteps;				/**< Maximum number of iterative refinement steps. */
	real_t epsIterRef;					/**< Early termination tolerance for iterative refinement. */
	real_t epsLITests;					/**< Tolerance for linear independence tests. */
	real_t epsNZCTests;					/**< Tolerance for nonzero curvature tests. */

	int enableDropInfeasibles;			/**< ... */
	int dropBoundPriority;				/**< ... */
	int dropEqConPriority;				/**< ... */
	int dropIneqConPriority;			/**< ... */

} qpOASES_Options;


int qpOASES_Options_init(	qpOASES_Options* const options,
							int mode
							);

int qpOASES_Options_copy(	const qpOASES_Options* const from,
							OptionsClass* const to
							);


int qpOASES_obtainOutputs(	const QProblemBClass* const globalQpObject,
							returnValue returnvalue,
							real_t* const x,
							real_t* const y,
							real_t* const obj,
							int* const status
							);


int QProblem_setup(	int nV,
					int nC,
					int hessianType
					);

int QProblem_init(	const real_t* const H,
					const real_t* const g,
					const real_t* const A,
					const real_t* const lb,
					const real_t* const ub,
					const real_t* const lbA,
					const real_t* const ubA,
					int* const nWSR,
					real_t* const cputime,
					const qpOASES_Options* const options,
					real_t* const x,
					real_t* const y,
					real_t* const obj,
					int* const status
					);

int QProblem_hotstart(	const real_t* const g,
						const real_t* const lb,
						const real_t* const ub,
						const real_t* const lbA,
						const real_t* const ubA,
						int* const nWSR,
						real_t* const cputime,
						real_t* const x,
						real_t* const y,
						real_t* const obj,
						int* const status
						);

int QProblem_cleanup( );



int QProblemB_setup(	int nV,
						int hessianType
						);

int QProblemB_init(	const real_t* const H,
					const real_t* const g,
					const real_t* const lb,
					const real_t* const ub,
					int* const nWSR,
					real_t* const cputime,
					const qpOASES_Options* const options,
					real_t* const x,
					real_t* const y,
					real_t* const obj,
					int* const status
					);

int QProblemB_hotstart(	const real_t* const g,
						const real_t* const lb,
						const real_t* const ub,
						int* const nWSR,
						real_t* const cputime,
						real_t* const x,
						real_t* const y,
						real_t* const obj,
						int* const status
						);

int QProblemB_cleanup( );



int SQProblem_setup(	int nV,
						int nC,
						int hessianType
						);

int SQProblem_init(	const real_t* const H,
					const real_t* const g,
					const real_t* const A,
					const real_t* const lb,
					const real_t* const ub,
					const real_t* const lbA,
					const real_t* const ubA,
					int* const nWSR,
					real_t* const cputime,
					const qpOASES_Options* const options,
					real_t* const x,
					real_t* const y,
					real_t* const obj,
					int* const status
					);

int SQProblem_hotstart(	const real_t* const H,
						const real_t* const g,
						const real_t* const A,
						const real_t* const lb,
						const real_t* const ub,
						const real_t* const lbA,
						const real_t* const ubA,
						int* const nWSR,
						real_t* const cputime,
						real_t* const x,
						real_t* const y,
						real_t* const obj,
						int* const status
						);

int SQProblem_cleanup( );


#endif /* QPOASES_WRAPPER_H */


/*
 *	end of file
 */
