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
 *	Declaration of the QProblemQore class, a wrapper around the QORE QP solver
 *	that is compatible to a subset of the interface provided by the QProblem
 *	class.
 */



#ifndef QPOASES_QPROBLEMQORE_HPP
#define QPOASES_QPROBLEMQORE_HPP


#include <qpOASES/Options.hpp>
#include <qpOASES/Matrices.hpp>
#include <qpOASES/QProblemB.hpp> // required only for declaration of 'Bounds'

extern "C" { 
	#include <qpsolver.h>
}

BEGIN_NAMESPACE_QPOASES


/**
 *	\brief A wrapper around the QORE QP solver that is compatible to a subset 
 *	of the interface provided by the QProblem class.
 *	\author Christian Hoffmann
 *	\date 2016-2017
 */
class QProblemQore 
{
	/*
	 *	PUBLIC MEMBER FUNCTIONS
	 */
	public:
		/** \brief Constructor which takes the QP dimension as input. 
		 *  \pre `_nV >= 1`
		 *  \pre `_nC >= 0` 
		 *  \todo Internally calls QORE's QPNew(), which dynamically allocates memory and may
		 *  thus fail unpredictably. The proper way of dealing with a failed allocation in 
		 *  the constructor would be to throw an exception, but exceptions are currently not
		 *  used in qpOASES. We currently rely on assertions instead. A failed allocation thus
		 *  leads to an immediate program termination!
		 */
		QProblemQore(	int_t _nV,	/**< Number of variables */
						int_t _nC	/**< Number of constraints. */
					);


		/** Destructor. */
		~QProblemQore( );


		/** \brief Overrides current options with given ones.
 		 *	\attention Not yet implemented, method has no effect.
 		 *	The underlying QORE implementation currently uses its default options.
 		 *	\return SUCCESSFUL_RETURN */
		returnValue setOptions(	const Options & _options	/**< New options. */
										);


		/** \brief Prints a list of all options and their current values.
 		 *	\attention Not yet implemented, method has no effect.
		 *	\return  SUCCESSFUL_RETURN \n */
		returnValue printOptions( ) const;


		/** \brief Initialises a bound-constrained QP and tries to solve it.
		 *  \pre `H != 0`
		 *  \pre `g != 0`
		 *  \pre `lb != 0`
		 *  \pre `ub != 0`
		 *  \note Arguments 5 to 10 are ignored. They are accepted only to ensure 
		 *  compatibility with the corresponding method of QProblem(B).
		 *  \attention Only minimalistic error handling implemented. On success, the method returns
		 *  SUCCESSFUL_RETURN. If any error occurrs in the method itself or in underlying QORE 
		 *  implementation, execution is terminated due to a violated assertion (in debug mode) 
		 *  or continues, possibly causing undefined behavior (in release mode).
		 *	\return SUCCESSFUL_RETURN
		 */
		returnValue init(	const real_t* const H, 				/**< Hessian matrix (a deep copy is made). \n
																	If Hessian matrix is trivial, a NULL pointer can be passed. */
							const real_t* const g,				/**< Gradient vector. */
							const real_t* const lb,				/**< Lower bounds (on variables). \n
																	If no lower bounds exist, a NULL pointer can be passed. */
							const real_t* const ub,				/**< Upper bounds (on variables). \n
																	If no upper bounds exist, a NULL pointer can be passed. */
							int_t &,	 						/**< ignored, implemented only for compatibility with QProblem and QProblemB */
				 			real_t * const = 0,					/**< ignored, implemented only for compatibility with QProblem and QProblemB */
							const real_t * const = 0,			/**< ignored, implemented only for compatibility with QProblem and QProblemB */
							const real_t * const = 0,			/**< ignored, implemented only for compatibility with QProblem and QProblemB */
							const Bounds * const = 0,			/**< ignored, implemented only for compatibility with QProblem and QProblemB */
							const real_t * const _R = 0			/**< ignored, implemented only for compatibility with QProblem and QProblemB */
							);

		
		/** \brief Solves an initialised bound-constrained QP sequence using the online active set strategy.
		 *	The QP solution is started from previous solution, if available.
		 *  \note Arguments 4 to 6 are ignored. They are accepted only to ensure 
		 *  compatibility with the corresponding method of QProblem(B).
		 *  \attention Only minimalistic error handling implemented. On success, the method returns
		 *  SUCCESSFUL_RETURN. If any error occurrs in the method itself or in underlying QORE 
		 *  implementation, execution is terminated due to a violated assertion (in debug mode) 
		 *  or continues, possibly causing undefined behavior (in release mode).
		 *	\return SUCCESSFUL_RETURN 
		 */
		returnValue hotstart(	const real_t* const g_new,				/**< Gradient of neighbouring QP to be solved. */
								const real_t* const lb_new,				/**< Lower bounds of neighbouring QP to be solved. \n
													 						 If no lower bounds exist, a NULL pointer can be passed. */
								const real_t* const ub_new,				/**< Upper bounds of neighbouring QP to be solved. \n
													 						 If no upper bounds exist, a NULL pointer can be passed. */
								int_t &,								/**< ignored, implemented only for compatibility with QProblem and QProblemB */
								real_t * const = 0,						/**< ignored, implemented only for compatibility with QProblem and QProblemB */
								const Bounds * const = 0				/**< ignored, implemented only for compatibility with QProblem and QProblemB */
								);

		
		/** \brief Initialises a QP and tries to solve it.
		 *  \pre `_H != 0`
		 *  \pre `_g != 0`
		 *  \pre `_A != 0`
		 *  \pre `_lb != 0`
		 *  \pre `_ub != 0`
		 *  \pre `_lbA != 0`
		 *  \pre `_ubA != 0`
		 *  \note Arguments 8 to 14 are ignored. They are accepted only to ensure 
		 *  compatibility with the corresponding method of QProblem(B).
		 *  \attention Only minimalistic error handling implemented. On success, the method returns
		 *  SUCCESSFUL_RETURN. If any error occurrs in the method itself or in underlying QORE 
		 *  implementation, execution is terminated due to a violated assertion (in debug mode) 
		 *  or continues, possibly causing undefined behavior (in release mode).
		 *  \return SUCCESSFUL_RETURN 
		 */
		returnValue init(	const real_t* const _H,							/**< Hessian matrix (a deep copy is made). */
							const real_t* const _g,							/**< Gradient vector. */
							const real_t* const _A,							/**< Constraint matrix (a deep copy is made). */
							const real_t* const _lb,						/**< Lower bound vector (on variables). \n
																				If no lower bounds exist, a NULL pointer can be passed. */
							const real_t* const _ub,						/**< Upper bound vector (on variables). \n
																				If no upper bounds exist, a NULL pointer can be passed. */
							const real_t* const _lbA,						/**< Lower constraints' bound vector. \n
																				If no lower constraints' bounds exist, a NULL pointer can be passed. */
							const real_t* const _ubA,						/**< Upper constraints' bound vector. \n
																				If no lower constraints' bounds exist, a NULL pointer can be passed. */
							int_t& ,										/**< ignored, implemented only for compatibility with QProblem and QProblemB */
							real_t* const  = 0,								/**< ignored, implemented only for compatibility with QProblem and QProblemB */
							const real_t* const = 0,						/**< ignored, implemented only for compatibility with QProblem and QProblemB */
							const real_t* const = 0,						/**< ignored, implemented only for compatibility with QProblem and QProblemB */
							const Bounds* const = 0,						/**< ignored, implemented only for compatibility with QProblem and QProblemB */
							const Constraints* const  = 0,					/**< ignored, implemented only for compatibility with QProblem and QProblemB */
							const real_t* const = 0							/**< ignored, implemented only for compatibility with QProblem and QProblemB */
							);

		
		/** \brief Solves an initialised QP sequence using the online active set strategy.
		 *	The QP solution is started from previous solution, if available.
		 *  \note Arguments 6 to 9 are ignored. They are accepted only to ensure 
		 *  compatibility with the corresponding method of QProblem(B).
		 *  \attention Only minimalistic error handling implemented. On success, the method returns
		 *  SUCCESSFUL_RETURN. If any error occurrs in the method itself or in underlying QORE 
		 *  implementation, execution is terminated due to a violated assertion (in debug mode) 
		 *  or continues, possibly causing undefined behavior (in release mode).
		 *	\return SUCCESSFUL_RETURN 
		 */
		returnValue hotstart(	const real_t* const g_new,						/**< Gradient of neighbouring QP to be solved. */
								const real_t* const lb_new,						/**< Lower bounds of neighbouring QP to be solved. \n
													 							 	 If no lower bounds exist, a NULL pointer can be passed. */
								const real_t* const ub_new,						/**< Upper bounds of neighbouring QP to be solved. \n
													 							 	 If no upper bounds exist, a NULL pointer can be passed. */
								const real_t* const lbA_new,					/**< Lower constraints' bounds of neighbouring QP to be solved. \n
													 							 	 If no lower constraints' bounds exist, a NULL pointer can be passed. */
								const real_t* const ubA_new,					/**< Upper constraints' bounds of neighbouring QP to be solved. \n
													 							 	 If no upper constraints' bounds exist, a NULL pointer can be passed. */
								int_t &,										/**< ignored, implemented only for compatibility with QProblem and QProblemB */
								real_t* const = 0,								/**< ignored, implemented only for compatibility with QProblem and QProblemB */
								const Bounds* const = 0,						/**< ignored, implemented only for compatibility with QProblem and QProblemB */
								const Constraints* const = 0					/**< ignored, implemented only for compatibility with QProblem and QProblemB */
								);

		
		/** \brief Returns the primal solution vector (deep copy).
		 *  \attention Only minimalistic error handling implemented. On success, the method returns
		 *  SUCCESSFUL_RETURN. If any error occurrs in the method itself or in underlying QORE 
		 *  implementation, execution is terminated due to a violated assertion (in debug mode) 
		 *  or continues, possibly causing undefined behavior (in release mode).
		 *	\return SUCCESSFUL_RETURN 
		 */
		returnValue getPrimalSolution(	real_t* const xOpt			/**< Output: Primal solution vector (if QP has been solved). */
										) const;


		/** \brief Returns the dual solution vector (deep copy).
		 *  \attention Only minimalistic error handling implemented. On success, the method returns
		 *  SUCCESSFUL_RETURN. If any error occurrs in the method itself or in underlying QORE 
		 *  implementation, execution is terminated due to a violated assertion (in debug mode) 
		 *  or continues, possibly causing undefined behavior (in release mode).
		 *	\return SUCCESSFUL_RETURN
		 */
		returnValue getDualSolution(	real_t* const yOpt	/**< Output: Dual solution vector (if QP has been solved). */
												) const;

		/** Returns the optimal objective function value.
 		 *	\attention Not yet implemented, method always returns `NaN`.
 		 *	\returns `NAN`
		 */
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
