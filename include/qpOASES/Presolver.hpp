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
 *	\file include/qpOASES/Presolver.hpp
 *	\author Hans Joachim Ferreau, Andreas Potschka, Christian Kirches, Dominik Cebulla
 *	\version 3.2
 *	\date 2007-2015
 *
 *	Declaration of the Presolver class designed for preprocessing a QP.
 */


#ifndef PRESOLVER_HPP
#define PRESOLVER_HPP

extern "C"
{
    #include <qpPresolver.h>
}

#include <qpOASES/PresolverOptions.hpp>
#include <qpOASES/Bounds.hpp>
#include <qpOASES/Constraints.hpp>
#include <qpOASES/Constants.hpp>
#include <qpOASES/Matrices.hpp>
#include <qpOASES/Types.hpp>


BEGIN_NAMESPACE_QPOASES


/* Redefinition of matrix sort type for qpPresolver in order to match the qpOASES
   naming convention. */
typedef qpp_matrix_sort_type_t MatrixSortType;
#define MST_ROW_WISE 		QPP_MST_ROW_WISE
#define MST_COLUMN_WISE 	QPP_MST_COLUMN_WISE

/** \brief Provides functionalities to presolve QPs (primarily to reduce their size).
 *
 *	The presolver must be used before solving the QP (via one of its presolve() methods)
 *	which then returns a preprocessed (possibly smaller) QP. This QP is passed to the solver.
 *	In order to obtain a primal-dual optimal solution of the original QP, you must use
 *	the postsolve() method afterwards.
 *
 *	CAVEAT: The impact of preprocessing depends on the QP data, i.e. on the matrices,
 *	gradient vector and bounds. Care must be taken if a sequence of QPs shall be presolved.
 *	We only recommend to presolve the initial QP, solving it, then postsolving it
 *	(i.e. retrieving the optimal working set and primal-dual solution of the
 *	original QP). Finally, the computed working set may be used as a warm-start for
 *	the sequence of QPs.
 */
class Presolver
{
public:

	/** Constructor of Presolver class.
	 *
	 *	Allocates (most of the) memory the presolver requires (e.g. for storing the QP).
	 *	If the number of nonzero elements of the matrices are not known, the presolver
	 *	reallocates memory for the matrices later.
	 *
	 *	\param nV Number of variables of the QP.
	 *	\param nC Number of constraints of the QP.
	 *	\param nzH Number of nonzero elements of the Hessian matrix. Default value: 0.
	 *		Note that the presolver only accesses the lower triangular part of
	 *		the Hessian matrix and hence only the number of nonzero elements of the lower
	 *		triangular part of the Hessian matrix must be passed to the presolver.
	 *	\param nzA Number of nonzero elements of the constraint matrix. Default value: 0.
	 *	\param presolveStackSize Initial size (= number of elements) of the presolve stack.
	 *		Default value: 1000. The presolve stack stores all necessary information of
	 *		the applied preprocessing techniques and will be reallocated if necessary.
	 */
	Presolver(const int_t nV,
			  const int_t nC,
			  const int_t nzH = 0,
			  const int_t nzA = 0,
			  const int_t presolveStackSize = 1000);


	/** Destructor of Presolver class. Deallocates all memory allocated by the presolver. */
	~Presolver();


	/** \brief Presolves the given QP with matrices stored in CSC (Harwell-Boeing) scheme.
	 *
	 *	All input parameters are going to be modified; they will contain the presolved QP.
	 *	A description of all parameters can be found in \a presolveCoordMat(), except
	 *	for those which are specifically introduced here.
	 *
	 *	\param Hcp Column pointers of the Hessian matrix in CSC scheme. Will be overwritten
	 *		with the presolved Hessian matrix. Is allowed to be a null pointer if
	 *		the Hessian matrix is empty. \n
	 *		CAVEAT: Only the lower triangular part of the Hessian matrix will be returned!
	 *	\param Acp Column pointers of the constraint matrix in CSC scheme. Will be
	 *		overwritten with the presolved constraint matrix. Is allowed to be a
	 *		null pointer if the constraint matrix is empty.
	 *
	 *	\return See \a presolveCoordMat()
	 */
	returnValue presolve(int_t* nV,
						 int_t* nC,
						 int_t* const Hirn,
						 int_t* const Hcp,
						 real_t* const Hx,
						 int_t* nzH,
						 real_t* const g,
						 int_t* const Airn,
						 int_t* const Acp,
						 real_t* const Ax,
						 int_t* nzA,
						 real_t* const lb,
						 real_t* const ub,
						 real_t* const lbA,
						 real_t* const ubA);


	/** \brief Presolves the given QP with matrices stored in \p SparseMatrix type.
	 *
	 *	All parameters, except for the given matrices \p H and \p A are going to be modified;
	 *	they contain the presolved QP. A description of all parameters can be found in \a
	 *	presolveCoordMat(), except for those which are specifically introduced here.
	 *
	 *	\param H Hessian matrix. (Internal) Data is not altered. User must provide diagonal
	 *		information via a previous call of H->createDiagInfo(). Is allowed to be a
	 *		null pointer if the Hessian matrix is empty.
	 *	\param Hprs Lower triangular part (!) of the presolved Hessian matrix. Memory is
	 *		allocated in this method and must be deallocated by the user.
	 *		Diagonal information of the matrix is not provided; the user must do this.
	 *		If \p H is a null pointer, then \p Hprs will also be a null pointer.
	 *	\param A Constraint matrix. (Internal) Data is not altered. Is allowed to be a
	 *		null pointer if the constraint matrix is empty.
	 *	\param Aprs Contains the presolved constraint matrix. Memory is allocated by
	 *		this method and must be deallocated by the user. If \p A
	 *		is a null pointer, then \p Aprs will also be a null pointer.
	 *
	 *	\return See \a presolveCoordMat().
	 */
	returnValue presolve(int_t* nV,
						 int_t* nC,
						 const SymSparseMat* const H,
						 SymSparseMat** Hprs,
						 real_t* const g,
						 const SparseMatrix* const A,
						 SparseMatrix** Aprs,
						 real_t* const lb,
						 real_t* const ub,
						 real_t* const lbA,
						 real_t* const ubA);


	/**	\brief Presolves the given QP with matrices given in dense (row-major) format.
	 *
	 *	All input parameters are going to be modified; they contain the presolved QP.
	 *	A description of all parameters can be found in \a presolveCoordMat(), except
	 *	for those which are specifically introduced here.
	 *
	 *	\param H Dense Hessian matrix in row-major format. Will contain the (full!)
	 *		presolved Hessian matrix. Is allowed to be a null pointer if the original
	 *		Hessian matrix is empty.
	 *	\param A Dense constraint matrix in row-major format. Will contain the
	 *		presolved constraint matrix. Is allowed to be a null pointer if the original
	 *		constraint matrix is empty. \n
	 *
	 *	\return See \a presolveCoordMat().
	 */
	returnValue presolve(int_t* nV,
						 int_t* nC,
						 real_t* const H,
						 real_t* const g,
						 real_t* const A,
						 real_t* const lb,
						 real_t* const ub,
						 real_t* const lbA,
						 real_t* const ubA);


	/** \brief Presolves the given QP with matrices stored in \p DenseMatrix type.
	 *
	 *	All parameters, except for the given matrices \p H and \p A are going to be modified;
	 *	they contain the presolved QP. A description of all parameters can be found in \a
	 *	presolveCoordMat(), except for those which are specifically introduced here.
	 *
	 *	\param H Dense Hessian matrix. (Internal) Data is not altered. Is allowed to be
	 *		a null pointer if the Hessian matrix is empty.
	 *	\param Hprs Will point to the presolved Hessian matrix. Both the lower and
	 *		upper triangular part are stored. Memory is allocated in this method and
	 *		must be deallocated by the user. If \p H is a null pointer, then
	 *		\p Hprs will also be a null pointer.
	 *	\param A Dense constraint matrix. (Internal) Data is not altered. Is allowed to
	 *		be a null pointer if the constraint matrix is empty.
	 *	\param Aprs Will point to the presolved constraint matrix. Memory is allocated in
	 *		this method and must be deallocated by the user.
	 *		If \p A is a null pointer, then \p Aprs will also be a null pointer.
	 *
	 *	\return See \a presolveCoordMat().
	 */
	returnValue presolve(int_t* nV,
						 int_t* nC,
						 const SymDenseMat* const H,
						 SymDenseMat** Hprs,
						 real_t* const g,
						 const DenseMatrix* const A,
						 DenseMatrix** Aprs,
						 real_t* const lb,
						 real_t* const ub,
						 real_t* const lbA,
						 real_t* const ubA);


	/**	\brief Postsolves the QP.
	 *
	 *	Computes a primal-dual optimal solution of the original QP on the basis of
	 *	the primal-dual solution of the presolved QP. Tries also to compute optimal
	 *	working sets of the original QP on the basis of the working sets of the presolved
	 *	QP if desired.
	 *
	 *	\param x Primal solution of the presolved QP. Will contain the primal solution of original QP.
	 *	\param y Dual solution of the presolved QP. Will contain the dual solution of original QP.
	 *	\param bounds Working set w.r.t. the bound constraints. If not a null pointer,
	 *		then this parameter will contain the optimal working set of the original QP afterwards.
	 *	\param constraints Working set w.r.t. the linear constraints. If not a null pointer,
	 *		then this parameter will contain the optimal working set of the original QP afterwards.
	 *
	 *	\return SUCCESSFUL_RETURN \n RET_ERROR_UNDEFINED \n RET_INVALID_ARGUMENTS
	 */
	returnValue postsolve(real_t* const x,
						  real_t* const y,
						  Bounds* const bounds = 0,
						  Constraints* const constraints = 0);


	/** Sets the internal options based on the given options.
	 *
	 *	This method identifies invalid option values and sets those values to default values.
	 *
	 *	\param opt Options.
	 *
	 *	\return SUCCESSFUL_RETURN
	 */
	inline returnValue setOptions(const PresolverOptions& opt);


	/** Return the current options of the presolver. */
	inline PresolverOptions getOptions() const;


	/** Print the current options of the presolver. */
	returnValue printOptions() const;

private:
	qpp_data_t* data;	/**< Struct containing all necessary data for qpPresolver. */


	/** \brief Presolves the given QP with matrices stored in coordinate scheme.
	 *
	 *	All input parameters are going to be modified; they contain the presolved QP
	 *	afterwards. Furthermore, all other presolve() methods call this function internally.
	 *
	 *	\param nV Will contain the number of variables of the presolved QP.
	 *	\param nC Will contain the number of constraints of the presolved QP.
	 *	\param Hirn Row subscripts of the Hessian matrix in coordinate format. Will be
	 *		overwritten with presolved Hessian matrix. Is allowed to be a null pointer if
	 *		the Hessian matrix is empty. \n
	 *		CAVEAT: Only the lower triangular part of the Hessian matrix will be returned!
	 *	\param Hjcn Column subscripts of the Hessian matrix in coordinate format. Will be
	 *		overwritten with presolved Hessian matrix. Is allowed to be a null pointer if
	 *		the Hessian matrix is empty. \n
	 *		CAVEAT: Only the lower triangular part of the Hessian matrix will be returned!
	 *	\param Hx Nonzero elements of the Hessian matrix in coordinate format. Will be
	 *		overwritten with presolved Hessian matrix. Is allowed to be a null pointer if
	 *		the Hessian matrix is empty. \n
	 *		CAVEAT: Only the lower triangular part of the Hessian matrix will be returned!
	 *	\param nzH Contains the length of the arrays \p Hirn, \p Hjcn and \p Hx. Will contain
	 *		the length of these arrays after presolving, i.e. the number of nonzero elements
	 *		of the lower triangular part of the Hessian matrix.
	 *	\param g Gradient vector of QP. Will be overwritten with presolved gradient vector.
	 *	\param Airn Row subscripts of the constraint matrix in coordinate format. Will be
	 *		overwritten with presolved constraint matrix. Is allowed to be a null pointer if
	 *		the constraint matrix is empty.
	 *	\param Ajcn Column subscripts of the constraint matrix in coordinate format. Will be
	 *		overwritten with presolved constraint matrix. Is allowed to be a null pointer if
	 *		the constraint matrix is empty.
	 *	\param Ax Nonzero elements of the constraint matrix in coordinate format. Will be
	 *		overwritten with presolved constraint matrix. Is allowed to be a null pointer if
	 *		the constraint matrix is empty.
	 *	\param nzA Contains the length of the arrays \p Airn, \p Ajcn and \p Ax. Will contain
	 *		the length of these arrays after presolving, i.e. the number of nonzero elements
	 *		of the constraint matrix.
	 *	\param lb Lower bounds on the variables. Will be overwritten with presolved bounds.
	 *	\param ub Upper bounds on the variables. Will be overwritten with presolved bounds.
	 *	\param lbA Lower bounds on the constraints. Will be overwritten with presolved bounds.
	 *	\param ubA Upper bounds on the constraints. Will be overwritten with presolved bounds.
	 *	\param sortType MST_ROW_WISE: Matrices of presolved QP are sorted row-wise. \n
	 *					MST_COLUMN_WISE: Matrices of presolved QP are sorted column-wise (default).
	 *
	 *	\return SUCCESSFUL_RETURN \n RET_INVALID_ARGUMENTS \n RET_QP_INFEASIBLE
	 *		\n RET_QP_UNBOUNDED \n RET_OPTIONS_ADJUSTED \n RET_ERROR_UNDEFINED
	 */
	returnValue presolveCoordMat(int_t* nV,
								 int_t* nC,
								 int_t* const Hirn,
								 int_t* const Hjcn,
								 real_t* const Hx,
								 int_t* nzH,
								 real_t* const g,
								 int_t* const Airn,
								 int_t* const Ajcn,
								 real_t* const Ax,
								 int_t* nzA,
								 real_t* const lb,
								 real_t* const ub,
								 real_t* const lbA,
								 real_t* const ubA,
								 const MatrixSortType& sortType = MST_COLUMN_WISE);

	/**	Converts error code from qpPresolver to error code from qpOASES (if possible).
	 *
	 *	\param err Error code from qpPresolver.
	 *
	 *	\return SUCCESSFUL_RETURN \n RET_INVALID_ARGUMENTS \n RET_QP_INFEASIBLE
	 *		\n RET_QP_UNBOUNDED \n RET_ERROR_UNDEFINED \n RET_OPTIONS_ADJUSTED
	 *		\n RET_UNABLE_TO_OPEN_FILE \n RET_UNABLE_TO_READ_FILE
	 *		\n RET_UNABLE_TO_WRITE_FILE
	 */
	returnValue convertErrorCode(const qpp_return_value_t err) const;
};


END_NAMESPACE_QPOASES

#include "Presolver.ipp"

#endif /* PRESOLVER_HPP */
