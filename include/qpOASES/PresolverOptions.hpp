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
 *	\file include/qpOASES/PresolverOptions.hpp
 *	\author Hans Joachim Ferreau, Andreas Potschka, Christian Kirches, Dominik Cebulla
 *	\version 3.2
 *	\date 2007-2015
 *
 *	Declaration of the PresolverOptions class designed to manage user-specified
 *	options for preprocessing a QP.
 */


#ifndef QPOASES_PRESOLVEROPTIONS_HPP
#define QPOASES_PRESOLVEROPTIONS_HPP

extern "C"
{
	#include <qpPresolver/options.h>
	#include <qpPresolver/types.h>
}

#include <qpOASES/MessageHandling.hpp>
#include <qpOASES/Types.hpp>
#include <qpOASES/Utils.hpp>


BEGIN_NAMESPACE_QPOASES


class Presolver;

/* Redefinition of bound types for qpPresolver in order to match the qpOASES
   naming convention. */
typedef qpp_bound_mode_t 	PresolverBoundType;
#define PBT_TIGHTEST 		QPP_BM_TIGHTEST_BOUNDS
#define PBT_MEDIUM 			QPP_BM_MEDIUM_BOUNDS


/** \brief Manages all user-specified options for preprocessing (= presolving) QPs.
 *
 *	This class manages all user-specified options used for preprocessing
 *	quadratic programs. The options can not be changed once the QP has been presolved.
 */
class PresolverOptions : private qpp_options_t
{
	friend class Presolver;

/*
 *	PUBLIC MEMBER FUNCTIONS
 */
public:
	/** Default constructor. Sets options to default values. */
	PresolverOptions();


	/** Copy constructor (deep copy).
	 *
	 *	\param rhs Rhs object.
	 */
	PresolverOptions(const PresolverOptions& rhs);


	/** Assignment operator (deep copy).
	 *
	 *	\param rhs Rhs object.
	 */
	PresolverOptions& operator=(const PresolverOptions& rhs);


	/** Sets all options to default values.
	 *
	 *	\return SUCCESSFUL_RETURN
	 */
	returnValue setToDefault();


	/** Sets options to values such that the QP does not get ill-conditioned by preprocessing.
	 *
	 *	\return SUCCESSFUL_RETURN
	 */
	returnValue setToReliable();


	/** Disables most of the more involved preprocessing techniques.
	 *
	 *	Although only very fast preprocessing techniques are used, the time for
     *  preprocessing may sometimes increase (due to less problem reduction). Also
     *  note that no extra bound tightening on the primal variables is performed.
     *  Hence there is no difference between the medium and the tightest bounds
     *  (we refer to the qpPresolver documentation for further details
     *  regarding the different sets of bounds).
	 *
	 *	\return SUCCESSFUL_RETURN
	 */
	returnValue setToFast();


	/** Enables all implemented preprocessing methods. This leads in most cases to the
	 *		best problem reduction.
	 *
	 *	\return SUCCESSFUL_RETURN
	 */
	returnValue setToMaxReduction();


	/** Prints values of all options.
	 *
	 *	\return SUCCESSFUL_RETURN
	 */
	returnValue print() const;


    /** \name Methods for changing certain option values.
     *
     *  Note that the passed options are not checked for consistency. This will later
     *  be done within the Presolver class (via \a Presolver::setOptions()).
     */
    /**@{*/

    /** Defines which bounds are passed to the user.
     *
     *  Currently supported are the tightest bounds (\a PBT_TIGHTEST) and the medium
     *  bounds (\a PBT_MEDIUM). */
	inline void setBoundMode(const PresolverBoundType& pbt);

	/** Sets the used tolerance when comparing two numbers for equality. */
	inline void setEqualityTol(const real_t eqTol);

	/** Sets the tolerance for ensuring numerical stability */
	inline void setStabilityTol(const real_t stabTol);

	/** Sets the feasibility tolerance, i.e. which constraint violation is tolerated
	 *      until the QP is treated as infeasible. */
	inline void setFeasibilityTol(const real_t feasTol);

	/** Sets the maximum number of iterations of the presolving loop. */
	inline void setMaxIter(const int_t maxIter);

	/** Sets the amount of information written to the logfile.
	 *
	 *  Value must be between 0 (= no information) and 5 (detailed information).
	 *  We refer to the qpPresolver documentation for further details. */
	inline void setLogfileLevel(const int_t logLevel);
	/**@}*/

	/** \name Methods for enabling/disabling certain preprocessing methods.
	 *
	 *  A description of the methods can be found e.g. in: Gould, N.I.M. and Toint, P.L.
	 *  "Preprocessing for quadratic programming." In: Mathematical Programming 100.1
	 *  (2004), pp. 95–132. */
	/**@{*/

	/** Enables bound tightening on the primal variables. */
	inline void enableBoundTightening();

	/** Enables methods based on the dual constraints of the QP; includes bound tightening
	 *      on the dual variables. */
	inline void enableDualConstraintsMethod();

	/** Enables the search for duplicate columns in the Hessian and constraint matrix. */
	inline void enableDuplicateColumnsMethod();

	/** Enables the removal of linearly unconstrained variables that do not occur in
	 *      mixed terms with other variables in the objective function. */
	inline void enableEmptyColumnsMethod();

	/** Enables methods based on the primal (= linear) constraints of the QP. */
	inline void enablePrimalConstraintsMethod();

	/** Applies scaling to the QP before and after the presolving loop. */
	inline void enableScaling();

	/** Enables the removal of (implied) free variables that do not occur in the Hessian
	 *      matrix and that occur in only a single linear constraint. */
	inline void enableSingletonColumnsMethod();

	/** Enables the removal of rows that contain only one nonzero element. */
	inline void enableSingletonRowsMethod();

	/** Enables the sparsification method which removes further nonzero elements from the
	 *      constraint matrix without creating any fill-in (-> Gaussian elimination). */
	inline void enableSparsificationMethod();

	/** Disables bound tightening on the primal variables. */
	inline void disableBoundTightening();

	/** Disables methods based on the dual constraints of the QP; includes bound tightening
	 *      on the dual variables. */
	inline void disableDualConstraintsMethod();

	/** Disables the search for duplicate columns in the Hessian and constraint matrix. */
	inline void disableDuplicateColumnsMethod();

	/** Enables the removal of linearly unconstrained variables that do not occur in
	 *      mixed terms with other variables in the objective function. */
	inline void disableEmptyColumnsMethod();

	/** Disables methods based on the primal (= linear) constraints of the QP. */
	inline void disablePrimalConstraintsMethod();

	/** QP is not scaled by the presolver. */
	inline void disableScaling();

	/** Disables the removal of (implied) free variables that do not occur in the Hessian
	 *      matrix and that occur in only a single linear constraint. */
	inline void disableSingletonColumnsMethod();

	/** Enables the removal of rows that contain only one nonzero element. */
	inline void disableSingletonRowsMethod();

	/** Enables the sparsification method which removes further nonzero elements from the
	 *      constraint matrix without creating any fill-in (-> Gaussian elimination). */
	inline void disableSparsificationMethod();
	/**@}*/

/*
 *	PRIVATE MEMBER FUNCTIONS
 */
private:
	/** Copies all members from given \p rhs object.
	 *
	 *	It is furthermore checked whether the given options are consistent. If this is
	 *  not the case, then the respective option values are set to default values and
	 *  a warning is given.
	 *
	 *	\param rhs Rhs object.
	 *
	 *	\return SUCCESSFUL_RETURN
	 */
	returnValue copy(const PresolverOptions& rhs);
};


END_NAMESPACE_QPOASES

#include "PresolverOptions.ipp"

#endif /* PRESOLVEROPTIONS_HPP */
