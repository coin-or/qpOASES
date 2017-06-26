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
 *	\file include/qpOASES/PresolverOptions.ipp
 *	\author Hans Joachim Ferreau, Andreas Potschka, Christian Kirches, Dominik Cebulla
 *	\version 3.2
 *	\date 2007-2015
 *
 *	Implementation of inlined member functions of the PresolverOptions class which
 *	is designed to manage user-specified options for preprocessing a QP.
 */


BEGIN_NAMESPACE_QPOASES


/*****************************************************************************
 *  P U B L I C                                                              *
 *****************************************************************************/

inline void PresolverOptions::setBoundMode(const PresolverBoundType& pbt)
{
    this->bound_mode = pbt;
}

inline void PresolverOptions::setEqualityTol(const real_t eqTol)
{
    this->eq_tol = eqTol;
}

inline void PresolverOptions::setStabilityTol(const real_t stabTol)
{
    this->stab_tol = stabTol;
}

inline void PresolverOptions::setFeasibilityTol(const real_t feasTol)
{
    this->feas_tol = feasTol;
}

inline void PresolverOptions::setMaxIter(const int_t maxIter)
{
    this->max_iter = maxIter;
}

inline void PresolverOptions::setLogfileLevel(const int_t logLevel)
{
    this->log_level = logLevel;
}

inline void PresolverOptions::enableBoundTightening()
{
    this->enable_bound_tightening = QPP_BT_TRUE;
}

inline void PresolverOptions::enableDualConstraintsMethod()
{
    this->enable_dual_constraints_method = QPP_BT_TRUE;
}

inline void PresolverOptions::enableDuplicateColumnsMethod()
{
    this->enable_duplicate_columns_method = QPP_BT_TRUE;
}

inline void PresolverOptions::enableEmptyColumnsMethod()
{
    this->enable_empty_columns_method = QPP_BT_TRUE;
}

inline void PresolverOptions::enablePrimalConstraintsMethod()
{
    this->enable_primal_constraints_method = QPP_BT_TRUE;
}

inline void PresolverOptions::enableScaling()
{
    this->enable_scaling = QPP_BT_TRUE;
}

inline void PresolverOptions::enableSingletonColumnsMethod()
{
    this->enable_singleton_columns_method = QPP_BT_TRUE;
}

inline void PresolverOptions::enableSingletonRowsMethod()
{
    this->enable_singleton_rows_method = QPP_BT_TRUE;
}

inline void PresolverOptions::enableSparsificationMethod()
{
    this->enable_sparsification_method = QPP_BT_TRUE;
}

inline void PresolverOptions::disableBoundTightening()
{
    this->enable_bound_tightening = QPP_BT_FALSE;
}

inline void PresolverOptions::disableDualConstraintsMethod()
{
    this->enable_dual_constraints_method = QPP_BT_FALSE;
}

inline void PresolverOptions::disableDuplicateColumnsMethod()
{
    this->enable_duplicate_columns_method = QPP_BT_FALSE;
}

inline void PresolverOptions::disableEmptyColumnsMethod()
{
    this->enable_empty_columns_method = QPP_BT_FALSE;
}

inline void PresolverOptions::disablePrimalConstraintsMethod()
{
    this->enable_primal_constraints_method = QPP_BT_FALSE;
}

inline void PresolverOptions::disableScaling()
{
    this->enable_scaling = QPP_BT_FALSE;
}

inline void PresolverOptions::disableSingletonColumnsMethod()
{
    this->enable_singleton_columns_method = QPP_BT_FALSE;
}

inline void PresolverOptions::disableSingletonRowsMethod()
{
    this->enable_singleton_rows_method = QPP_BT_FALSE;
}

inline void PresolverOptions::disableSparsificationMethod()
{
    this->enable_sparsification_method = QPP_BT_FALSE;
}

END_NAMESPACE_QPOASES
