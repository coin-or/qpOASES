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
 *	\file include/qpOASES/QProblemQore.hpp
 *	\author Christian Hoffmann
 *	\date 2016-2017
 *
 *	Declaration of the QProblemQore class,  a QProblem-compatible wrapper 
 *	around the QORE QP solver.
 */



#ifndef QPOASES_QPROBLEMQORE_HPP
#define QPOASES_QPROBLEMQORE_HPP


#include <qpOASES/Options.hpp>
// #include <qpOASES/QProblemB.hpp>
// #include <qpOASES/Constraints.hpp>
// #include <qpOASES/ConstraintProduct.hpp>
#include <qpOASES/Matrices.hpp>

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
		QoreProblem * pproblem;
};


END_NAMESPACE_QPOASES

#endif	/* QPOASES_QPROBLEMQORE_HPP */


/*
 *	end of file
 */
