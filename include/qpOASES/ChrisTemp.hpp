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
 *	\file include/qpOASES/ChrisTemp.hpp
 *	\author Christian Hoffmann
 *	\version 3.2
 *	\date 2016
 *
 *	File temporarily containing Chris' declarations during development.
 */



#ifndef QPOASES_CHRISTEMP_HPP
#define QPOASES_CHRISTEMP_HPP


// #include <qpOASES/Flipper.hpp>
// #include <qpOASES/Options.hpp>
#include <qpOASES/Matrices.hpp>

// FIXME temporary (chris)
extern "C" {
	#include <spmatrix.h>
}

BEGIN_NAMESPACE_QPOASES

/**	\brief Convert sparse QORE matrix to sparse qpOASES matrix (deep copy).
 * 	\author Christian Hoffmann
 * 	\date 2016
 */
SparseMatrix fromQoreMatrix(
	spmatrix const & /**< sparse matrix in QORE format */
);

END_NAMESPACE_QPOASES

// #include <qpOASES/QProblemBase.ipp>

#endif	/* QPOASES_CHRISTEMP_HPP */


/*
 *	end of file
 */
