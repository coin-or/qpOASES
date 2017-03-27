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
 *	\file QProblemQore.hpp
 *	\author Christian Hoffmann
 *	\date 2016-2017
 *
 *	Declaration of the QProblemQore class,  a QProblem-compatible wrapper 
 *	around the QORE QP solver.
 */



#ifndef QPOASES_QPROBLEMQORE_HPP
#define QPOASES_QPROBLEMQORE_HPP


#include <qpOASES/Options.hpp>
// #include <qpOASES/Constraints.hpp>
// #include <qpOASES/ConstraintProduct.hpp>
#include <qpOASES/Matrices.hpp>
#include <qpOASES/QProblemB.hpp> // FIXME overkill, just included for declaration of 'Bounds' type

extern "C" { 
	#include <qpsolver.h>
}

BEGIN_NAMESPACE_QPOASES


/**
 *	\brief Variant of the QProblem class based on the QORE QP solver. 
 *
 *	\author Christian Hoffmann
 *	\date 2016-2017
 */
class QProblemQore // : public QProblem
{
	/*
	 *	PUBLIC MEMBER FUNCTIONS
	 */
	public:
		/** Constructor which takes the QP dimension as input. */
		QProblemQore(	int_t _nV,	/**< Number of variables. */
						int_t _nC	/**< Number of constraints. */
					);

		
		/** Destructor. */
		~QProblemQore( );

		/** Overrides current options with given ones.
 		 *	\return SUCCESSFUL_RETURN */
		returnValue setOptions(	const Options & _options	/**< New options. */
										);

		
		/** Prints a list of all options and their current values.
		 *	\return  SUCCESSFUL_RETURN \n */
		returnValue printOptions( ) const;

		
		/** Initialises a simply bounded QP problem with given QP data and tries to solve it
		 *	using at most nWSR iterations.
		 *
		 *	\return SUCCESSFUL_RETURN \n
					RET_INIT_FAILED \n
					RET_INIT_FAILED_CHOLESKY \n
					RET_INIT_FAILED_HOTSTART \n
					RET_INIT_FAILED_INFEASIBILITY \n
					RET_INIT_FAILED_UNBOUNDEDNESS \n
					RET_MAX_NWSR_REACHED \n
					RET_INVALID_ARGUMENTS */
		returnValue init(	const real_t* const _H, 				/**< Hessian matrix (a shallow copy is made). \n
																		 If Hessian matrix is trivial, a NULL pointer can be passed. */
							const real_t* const _g,					/**< Gradient vector. */
							const real_t* const _lb,				/**< Lower bounds (on variables). \n
																		 If no lower bounds exist, a NULL pointer can be passed. */
							const real_t* const _ub,				/**< Upper bounds (on variables). \n
																		 If no upper bounds exist, a NULL pointer can be passed. */
							int_t& nWSR 							/**< Input: Maximum number of working set recalculations when using initial homotopy. \n
																		 Output: Number of performed working set recalculations. */
							);

		/** Initialises a QP problem with given QP data and tries to solve it
		 *	using at most nWSR iterations. 
		 *	
		 *	\return SUCCESSFUL_RETURN \n
					RET_INIT_FAILED \n
					RET_INIT_FAILED_CHOLESKY \n
					RET_INIT_FAILED_TQ \n
					RET_INIT_FAILED_HOTSTART \n
					RET_INIT_FAILED_INFEASIBILITY \n
					RET_INIT_FAILED_UNBOUNDEDNESS \n
					RET_MAX_NWSR_REACHED \n
					RET_INVALID_ARGUMENTS */
		returnValue init(	const real_t* const _H,							/**< Hessian matrix (a shallow copy is made). \n
																				 If Hessian matrix is trivial, a NULL pointer can be passed. */
							const real_t* const _g,							/**< Gradient vector. */
							const real_t* const _A,							/**< Constraint matrix (a shallow copy is made). */
							const real_t* const _lb,						/**< Lower bound vector (on variables). \n
																				 If no lower bounds exist, a NULL pointer can be passed. */
							const real_t* const _ub,						/**< Upper bound vector (on variables). \n
																				 If no upper bounds exist, a NULL pointer can be passed. */
							const real_t* const _lbA,						/**< Lower constraints' bound vector. \n
																				 If no lower constraints' bounds exist, a NULL pointer can be passed. */
							const real_t* const _ubA,						/**< Upper constraints' bound vector. \n
																				 If no lower constraints' bounds exist, a NULL pointer can be passed. */
							int_t& nWSR									/**< Input: Maximum number of working set recalculations when using initial homotopy.
																				 Output: Number of performed working set recalculations. */
							);

		
		/** Solves an initialised QP sequence using the online active set strategy.
		 *	By default, QP solution is started from previous solution. If a guess
		 *	for the working set is provided, an initialised homotopy is performed.
		 *
		 *  Note: This function internally calls solveQP/solveRegularisedQP
		 *        for solving an initialised QP!
		 *
		 *	\return SUCCESSFUL_RETURN \n
					RET_MAX_NWSR_REACHED \n
					RET_HOTSTART_FAILED_AS_QP_NOT_INITIALISED \n
					RET_HOTSTART_FAILED \n
					RET_SHIFT_DETERMINATION_FAILED \n
					RET_STEPDIRECTION_DETERMINATION_FAILED \n
					RET_STEPLENGTH_DETERMINATION_FAILED \n
					RET_HOMOTOPY_STEP_FAILED \n
					RET_HOTSTART_STOPPED_INFEASIBILITY \n
					RET_HOTSTART_STOPPED_UNBOUNDEDNESS \n
					RET_SETUP_AUXILIARYQP_FAILED */
		returnValue hotstart(	const real_t* const g_new,				/**< Gradient of neighbouring QP to be solved. */
								const real_t* const lb_new,				/**< Lower bounds of neighbouring QP to be solved. \n
													 						 If no lower bounds exist, a NULL pointer can be passed. */
								const real_t* const ub_new,				/**< Upper bounds of neighbouring QP to be solved. \n
													 						 If no upper bounds exist, a NULL pointer can be passed. */
								int_t& nWSR,							/**< Input: Maximum number of working set recalculations; \n
																			 Output: Number of performed working set recalculations. */
								real_t* const cputime = 0,				/**< Input: Maximum CPU time allowed for QP solution. \n
																			 Output: CPU time spent for QP solution (or to perform nWSR iterations). */
								const Bounds* const guessedBounds = 0	/**< Optimal working set of bounds for solution (xOpt,yOpt). \n
																			 (If a null pointer is passed, the previous working set is kept!) */
								);

		
		/** Initialises a QP problem with given QP data and tries to solve it
		 *	using at most nWSR iterations. Depending on the parameter constellation it: \n
		 *	1. 0,    0,    0    : starts with xOpt = 0, yOpt = 0 and gB/gC empty (or all implicit equality bounds), \n
		 *	2. xOpt, 0,    0    : starts with xOpt, yOpt = 0 and obtain gB/gC by "clipping", \n
		 *	3. 0,    yOpt, 0    : starts with xOpt = 0, yOpt and obtain gB/gC from yOpt != 0, \n
		 *	4. 0,    0,    gB/gC: starts with xOpt = 0, yOpt = 0 and gB/gC, \n
		 *	5. xOpt, yOpt, 0    : starts with xOpt, yOpt and obtain gB/gC from yOpt != 0, \n
		 *	6. xOpt, 0,    gB/gC: starts with xOpt, yOpt = 0 and gB/gC, \n
		 *	7. xOpt, yOpt, gB/gC: starts with xOpt, yOpt and gB/gC (assume them to be consistent!)
		 *
		 *  Note: This function internally calls solveInitialQP for initialisation!
		 *
		  *	\return SUCCESSFUL_RETURN \n
					RET_INIT_FAILED \n
					RET_INIT_FAILED_CHOLESKY \n
					RET_INIT_FAILED_TQ \n
					RET_INIT_FAILED_HOTSTART \n
					RET_INIT_FAILED_INFEASIBILITY \n
					RET_INIT_FAILED_UNBOUNDEDNESS \n
					RET_MAX_NWSR_REACHED \n
					RET_INVALID_ARGUMENTS */
		returnValue init(	SymmetricMatrix *_H,							/**< Hessian matrix (a shallow copy is made). */
							const real_t* const _g, 						/**< Gradient vector. */
							Matrix *_A,  									/**< Constraint matrix (a shallow copy is made). */
							const real_t* const _lb,						/**< Lower bound vector (on variables). \n
																				 If no lower bounds exist, a NULL pointer can be passed. */
							const real_t* const _ub,						/**< Upper bound vector (on variables). \n
																				 If no upper bounds exist, a NULL pointer can be passed. */
							const real_t* const _lbA,						/**< Lower constraints' bound vector. \n
																				 If no lower constraints' bounds exist, a NULL pointer can be passed. */
							const real_t* const _ubA,						/**< Upper constraints' bound vector. \n
																				 If no lower constraints' bounds exist, a NULL pointer can be passed. */
							int_t& nWSR,									/**< Input: Maximum number of working set recalculations when using initial homotopy.
																				 Output: Number of performed working set recalculations. */
							real_t* const cputime = 0,						/**< Input: Maximum CPU time allowed for QP initialisation. \n
																				 Output: CPU time spent for QP initialisation (if pointer passed). */
							const real_t* const xOpt = 0,					/**< Optimal primal solution vector. \n
																				 (If a null pointer is passed, the old primal solution is kept!) */
							const real_t* const yOpt = 0,					/**< Optimal dual solution vector. \n
																				 (If a null pointer is passed, the old dual solution is kept!) */
							const Bounds* const guessedBounds = 0,			/**< Optimal working set of bounds for solution (xOpt,yOpt). \n
																				 (If a null pointer is passed, all bounds are assumed inactive!) */
							const Constraints* const guessedConstraints = 0,/**< Optimal working set of constraints for solution (xOpt,yOpt). \n
																				 (If a null pointer is passed, all constraints are assumed inactive!) */
							const real_t* const _R = 0						/**< Pre-computed (upper triangular) Cholesky factor of Hessian matrix.
																			 	 The Cholesky factor must be stored in a real_t array of size nV*nV
																				 in row-major format. Note: Only used if xOpt/yOpt and gB are NULL! \n
																				 (If a null pointer is passed, Cholesky decomposition is computed internally!) */
							);
		
		/** Solves an initialised QP sequence using the online active set strategy.
		 *	By default, QP solution is started from previous solution. If a guess
		 *	for the working set is provided, an initialised homotopy is performed.
		 *
		 *  Note: This function internally calls solveQP/solveRegularisedQP
		 *        for solving an initialised QP!
		 *
		 *	\return SUCCESSFUL_RETURN \n
		 			RET_MAX_NWSR_REACHED \n
		 			RET_HOTSTART_FAILED_AS_QP_NOT_INITIALISED \n
					RET_HOTSTART_FAILED \n
					RET_SHIFT_DETERMINATION_FAILED \n
					RET_STEPDIRECTION_DETERMINATION_FAILED \n
					RET_STEPLENGTH_DETERMINATION_FAILED \n
					RET_HOMOTOPY_STEP_FAILED \n
					RET_HOTSTART_STOPPED_INFEASIBILITY \n
					RET_HOTSTART_STOPPED_UNBOUNDEDNESS */
		returnValue hotstart(	const real_t* const g_new,						/**< Gradient of neighbouring QP to be solved. */
								const real_t* const lb_new,						/**< Lower bounds of neighbouring QP to be solved. \n
													 							 	 If no lower bounds exist, a NULL pointer can be passed. */
								const real_t* const ub_new,						/**< Upper bounds of neighbouring QP to be solved. \n
													 							 	 If no upper bounds exist, a NULL pointer can be passed. */
								const real_t* const lbA_new,					/**< Lower constraints' bounds of neighbouring QP to be solved. \n
													 							 	 If no lower constraints' bounds exist, a NULL pointer can be passed. */
								const real_t* const ubA_new,					/**< Upper constraints' bounds of neighbouring QP to be solved. \n
													 							 	 If no upper constraints' bounds exist, a NULL pointer can be passed. */
								int_t& nWSR,									/**< Input: Maximum number of working set recalculations; \n
																			 		 Output: Number of performed working set recalculations. */
								real_t* const cputime = 0,						/**< Input: Maximum CPU time allowed for QP solution. \n
																				 	 Output: CPU time spent for QP solution (or to perform nWSR iterations). */
								const Bounds* const guessedBounds = 0,			/**< Optimal working set of bounds for solution (xOpt,yOpt). \n
																					 (If a null pointer is passed, the previous working set of bounds is kept!) */
								const Constraints* const guessedConstraints = 0	/**< Optimal working set of constraints for solution (xOpt,yOpt). \n
																					 (If a null pointer is passed, the previous working set of constraints is kept!) */
								);

		/** Returns the primal solution vector.
		 *	\return SUCCESSFUL_RETURN \n
					RET_QP_NOT_SOLVED */
		returnValue getPrimalSolution(	real_t* const xOpt			/**< Output: Primal solution vector (if QP has been solved). */
										) const;

										
		/** Returns the dual solution vector (deep copy).
		 *	\return SUCCESSFUL_RETURN \n
					RET_QP_NOT_SOLVED */
		returnValue getDualSolution(	real_t* const yOpt	/**< Output: Dual solution vector (if QP has been solved). */
												) const;

		/** Returns the optimal objective function value.
		 *	\return finite value: Optimal objective function value (QP was solved) \n
		 			+infinity:	  QP was not yet solved */
		real_t getObjVal( ) const;
		
	/*
	 *	PRIVATE MEMBER VARIABLES
	 */
	private:
		QoreProblem * pproblem; ///< pointer to a QORE QP solver instance
		int_t nV; ///< number of variables
		int_t nC; ///< number of constraints
};


END_NAMESPACE_QPOASES

#endif	/* QPOASES_QPROBLEMQORE_HPP */


/*
 *	end of file
 */
