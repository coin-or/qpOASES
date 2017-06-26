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
 *	\file src/PresolverOptions.cpp
 *	\author Hans Joachim Ferreau, Andreas Potschka, Christian Kirches, Dominik Cebulla
 *	\version 3.2
 *	\date 2007-2015
 *
 *	Implementation of the PresolverOptions class which is designed to manage user-specified
 *	options for preprocessing a QP.
 */


#include <qpOASES/PresolverOptions.hpp>


BEGIN_NAMESPACE_QPOASES


/*****************************************************************************
 *  P U B L I C                                                              *
 *****************************************************************************/

PresolverOptions::PresolverOptions() : qpp_options_t()
{
    qppSetDefaultOptions(this);
}


PresolverOptions::PresolverOptions(const PresolverOptions& rhs)
{
    copy(rhs);
}


PresolverOptions& PresolverOptions::operator=(const PresolverOptions& rhs)
{
    if (this != &rhs)
    {
        copy(rhs);
    }
    return *this;
}


returnValue PresolverOptions::setToDefault()
{
    qppSetDefaultOptions(this);

    return SUCCESSFUL_RETURN;
}


returnValue PresolverOptions::setToReliable()
{
    setToDefault();

    stab_tol = 1e-2;
    feas_tol /= 1e2;
    enable_sparsification_method = QPP_BT_FALSE;

    return SUCCESSFUL_RETURN;
}


returnValue PresolverOptions::setToFast()
{
    setToDefault();

    enable_duplicate_columns_method = QPP_BT_FALSE;
    enable_dual_constraints_method = QPP_BT_FALSE;
    enable_bound_tightening = QPP_BT_FALSE;
    max_iter = 10;

    return SUCCESSFUL_RETURN;
}


returnValue PresolverOptions::setToMaxReduction()
{
    setToDefault();

    stab_tol = 1e-4;
    max_iter = 50;
    enable_duplicate_columns_method = QPP_BT_TRUE;

    return SUCCESSFUL_RETURN;
}

returnValue PresolverOptions::print() const
{
    /*  Same print style as in the Options class of qpOASES. Otherwise one can use
        qppPrintOptions() from qpPresolver to print all option values. */

    #ifndef __SUPPRESSANYOUTPUT__

	char myPrintfString[MAX_STRING_LENGTH];
	char boundInfo[32];

	if (bound_mode == PBT_MEDIUM)
    {
        strncpy(boundInfo, "Medium Bounds", 32);
    }
    else
    {
        strncpy(boundInfo, "Tightest Bounds", 32);
    }

	myPrintf("\n###################   qpOASES  --  PRESOLVER OPTIONS   ##################\n");
	myPrintf("\n");

    snprintf(myPrintfString, MAX_STRING_LENGTH, "bound_mode                    =  %s\n",
             boundInfo);
	myPrintf(myPrintfString);

	snprintf(myPrintfString, MAX_STRING_LENGTH, "max_iter                      =  %"
          QPP_PRID "\n", max_iter);
	myPrintf(myPrintfString);

	snprintf(myPrintfString, MAX_STRING_LENGTH, "log_level                     =  %"
          QPP_PRID "\n", log_level);
	myPrintf(myPrintfString);

	snprintf(myPrintfString, MAX_STRING_LENGTH, "eq_tol                        =  %e\n",
          eq_tol);
	myPrintf(myPrintfString);

	snprintf(myPrintfString, MAX_STRING_LENGTH, "feas_tol                      =  %e\n",
          feas_tol);
	myPrintf(myPrintfString);

	snprintf(myPrintfString, MAX_STRING_LENGTH, "stab_tol                      =  %e\n",
          stab_tol);
	myPrintf(myPrintfString);

	myPrintf("\nEnabled Preprocessing Methods (enable_):\n");

	snprintf(myPrintfString, MAX_STRING_LENGTH, "bound_tightening              =  %"
          QPP_PRID "\n", enable_bound_tightening);
	myPrintf(myPrintfString);

	snprintf(myPrintfString, MAX_STRING_LENGTH, "dual_constraints_method       =  %"
          QPP_PRID "\n", enable_dual_constraints_method);
	myPrintf(myPrintfString);

	snprintf(myPrintfString, MAX_STRING_LENGTH, "duplicate_columns_method      =  %"
          QPP_PRID "\n", enable_duplicate_columns_method);
	myPrintf(myPrintfString);

	snprintf(myPrintfString, MAX_STRING_LENGTH, "empty_columns_method          =  %"
          QPP_PRID "\n", enable_empty_columns_method);
	myPrintf(myPrintfString);

	snprintf(myPrintfString, MAX_STRING_LENGTH, "primal_constraints_method     =  %"
          QPP_PRID "\n", enable_primal_constraints_method);
	myPrintf(myPrintfString);

	snprintf(myPrintfString, MAX_STRING_LENGTH, "scaling                       =  %"
          QPP_PRID "\n", enable_scaling);
	myPrintf(myPrintfString);

	snprintf(myPrintfString, MAX_STRING_LENGTH, "singleton_columns_method      =  %"
          QPP_PRID "\n", enable_singleton_columns_method);
	myPrintf(myPrintfString);

	snprintf(myPrintfString, MAX_STRING_LENGTH, "singleton_rows_method         =  %"
          QPP_PRID "\n", enable_singleton_rows_method);
	myPrintf(myPrintfString);

	snprintf(myPrintfString, MAX_STRING_LENGTH, "sparsification_method         =  %"
          QPP_PRID "\n", enable_sparsification_method);
	myPrintf(myPrintfString);

	myPrintf("\n\n");

	#endif  /* __SUPPRESSANYOUTPUT__ */

    return SUCCESSFUL_RETURN;
}


/*****************************************************************************
 *  P R I V A T E                                                            *
 *****************************************************************************/

returnValue PresolverOptions::copy(const PresolverOptions& rhs)
{
    qppCopyOptions(this, &rhs);

    if (qppCheckOptions(this) != QPP_OK)
    {
        THROWWARNING( RET_OPTIONS_ADJUSTED );
    }

    return SUCCESSFUL_RETURN;
}


END_NAMESPACE_QPOASES
