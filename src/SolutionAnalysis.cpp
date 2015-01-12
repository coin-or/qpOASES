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
 *	\file src/SolutionAnalysis.cpp
 *	\author Hans Joachim Ferreau (thanks to Boris Houska)
 *	\version 3.1
 *	\date 2008-2015
 *
 *	Implementation of the SolutionAnalysis class designed to perform
 *	additional analysis after solving a QP with qpOASES.
 *
 */


#include <qpOASES/extras/SolutionAnalysis.hpp>


BEGIN_NAMESPACE_QPOASES


/*****************************************************************************
 *  P U B L I C                                                              *
 *****************************************************************************/


/*
 *	S o l u t i o n A n a l y s i s
 */
SolutionAnalysis::SolutionAnalysis( )
{

}


/*
 *	S o l u t i o n A n a l y s i s
 */
SolutionAnalysis::SolutionAnalysis( const SolutionAnalysis& rhs )
{

}


/*
 *	~ S o l u t i o n A n a l y s i s
 */
SolutionAnalysis::~SolutionAnalysis( )
{

}


/*
 *	o p e r a t o r =
 */
SolutionAnalysis& SolutionAnalysis::operator=( const SolutionAnalysis& rhs )
{
	if ( this != &rhs )
	{

	}

	return *this;
}



/*
 *	g e t K k t V i o l a t i o n
 */
real_t SolutionAnalysis::getKktViolation(	QProblemB* const qp,
											real_t* const maxStat, real_t* const maxFeas, real_t* const maxCmpl
											) const
{
	int i;
	int nV = qp->getNV();

	if ( qp == 0 )
		return INFTY;

	/* setup Hessian matrix array (or pass NULL pointer) */
	real_t* H_ptr = 0;
	BooleanType hasIdentityHessian = BT_FALSE;

	switch( qp->getHessianType() )
	{
		case HST_ZERO:
			break;

		case HST_IDENTITY:
			hasIdentityHessian = BT_TRUE;
			break;

		default:
			H_ptr = qp->H->full();
			if ( qp->usingRegularisation() == BT_TRUE )
				for( i=0; i<nV; ++i )
					H_ptr[i*nV+i] -= qp->regVal;
	}

	real_t* workingSetB = new real_t[nV];
	qp->getWorkingSetBounds( workingSetB );

	/* determine maximum KKT violation */
	real_t maxKktViolation=0.0, stat=0.0, feas=0.0, cmpl=0.0;

	returnValue returnvalue = REFER_NAMESPACE_QPOASES getKktViolation(	nV,
																		H_ptr,qp->g,
																		qp->lb,qp->ub,
																		qp->x,qp->y,
																		stat,feas,cmpl,
																		workingSetB,hasIdentityHessian
																		);
	if ( workingSetB != 0 )
		delete[] workingSetB;

	if ( H_ptr != 0 )
		delete[] H_ptr;

	if ( returnvalue != SUCCESSFUL_RETURN )
		THROWERROR( returnvalue );

	/* assign return values */
	if ( maxStat != 0 )
		*maxStat = stat;

	if ( maxFeas != 0 )
		*maxFeas = feas;

	if ( maxCmpl != 0 )
		*maxCmpl = cmpl;

	maxKktViolation = getMax( maxKktViolation,stat );
	maxKktViolation = getMax( maxKktViolation,feas );
	maxKktViolation = getMax( maxKktViolation,cmpl );

	return maxKktViolation;
}


/*
 *	g e t K k t V i o l a t i o n
 */
real_t SolutionAnalysis::getKktViolation(	QProblem* const qp,
											real_t* const maxStat, real_t* const maxFeas, real_t* const maxCmpl
											) const
{
	int i;
	int nV = qp->getNV();
	int nC = qp->getNC();

	if ( qp == 0 )
		return INFTY;

	/* setup Hessian matrix array (or pass NULL pointer) */
	real_t* H_ptr = 0;
	BooleanType hasIdentityHessian = BT_FALSE;

	switch( qp->getHessianType() )
	{
		case HST_ZERO:
			break;

		case HST_IDENTITY:
			hasIdentityHessian = BT_TRUE;
			break;

		default:
			H_ptr = qp->H->full();
			if ( qp->usingRegularisation() == BT_TRUE )
				for( i=0; i<nV; ++i )
					H_ptr[i*nV+i] -= qp->regVal;
	}

	/* setup constraint matrix array */
	real_t* A_ptr = qp->A->full();

	real_t* workingSetB = new real_t[nV];
	qp->getWorkingSetBounds( workingSetB );

	real_t* workingSetC = new real_t[nC];
	qp->getWorkingSetConstraints( workingSetC );

	/* determine maximum KKT violation */
	real_t maxKktViolation=0.0, stat=0.0, feas=0.0, cmpl=0.0;

	returnValue returnvalue = REFER_NAMESPACE_QPOASES getKktViolation(	nV,nC,
																		H_ptr,qp->g,A_ptr,
																		qp->lb,qp->ub,qp->lbA,qp->ubA,
																		qp->x,qp->y,
																		stat,feas,cmpl,
																		workingSetB,workingSetC,hasIdentityHessian
																		);

	if ( workingSetC != 0 )
		delete[] workingSetC;

	if ( workingSetB != 0 )
		delete[] workingSetB;

	if ( A_ptr != 0 )
		delete[] A_ptr;
	
	if ( H_ptr != 0 )
		delete[] H_ptr;

	if ( returnvalue != SUCCESSFUL_RETURN )
		THROWERROR( returnvalue );

	/* assign return values */
	if ( maxStat != 0 )
		*maxStat = stat;

	if ( maxFeas != 0 )
		*maxFeas = feas;

	if ( maxCmpl != 0 )
		*maxCmpl = cmpl;

	maxKktViolation = getMax( maxKktViolation,stat );
	maxKktViolation = getMax( maxKktViolation,feas );
	maxKktViolation = getMax( maxKktViolation,cmpl );

	return maxKktViolation;
}


/*
 *	g e t K k t V i o l a t i o n
 */
real_t SolutionAnalysis::getKktViolation(	SQProblem* const qp,
											real_t* const maxStat, real_t* const maxFeas, real_t* const maxCmpl
											) const
{
	return getKktViolation( (QProblem*)qp, maxStat,maxFeas,maxCmpl );
}



/*
 *	g e t V a r i a n c e C o v a r i a n c e
 */
returnValue SolutionAnalysis::getVarianceCovariance(	QProblemB* const qp,
														const real_t* const g_b_bA_VAR, real_t* const Primal_Dual_VAR
														) const
{
	return THROWERROR( RET_NOT_YET_IMPLEMENTED );
}


/*
 *	g e t V a r i a n c e C o v a r i a n c e
 */
returnValue SolutionAnalysis::getVarianceCovariance(	QProblem* qp,
														const real_t* const g_b_bA_VAR, real_t* const Primal_Dual_VAR 
														) const
{

  /* DEFINITION OF THE DIMENSIONS nV AND nC:
   * --------------------------------------- */
  int nV  = qp->getNV( );                      /* dimension of x / the bounds */
  int nC  = qp->getNC( );                      /* dimension of the constraints */
  int dim = 2*nV+nC;                           /* dimension of input and output */
                                               /* variance-covariance matrix */
  int run1, run2, run3;                        /* simple run variables (for loops). */


  /* ALLOCATION OF MEMORY:
   * --------------------- */
  real_t* delta_g_cov    = new real_t[nV];     /* a covariance-vector of g */
  real_t* delta_lb_cov   = new real_t[nV];     /* a covariance-vector of lb */
  real_t* delta_ub_cov   = new real_t[nV];     /* a covariance-vector of ub */
  real_t* delta_lbA_cov  = new real_t[nC];     /* a covariance-vector of lbA */
  real_t* delta_ubA_cov  = new real_t[nC];     /* a covariance-vector of ubA */

  returnValue returnvalue;                     /* the return value */
  BooleanType Delta_bC_isZero = BT_FALSE;      /* (just use FALSE here) */
  BooleanType Delta_bB_isZero = BT_FALSE;      /* (just use FALSE here) */



  /* ASK FOR THE NUMBER OF FREE AND FIXED VARIABLES:
   * (ASSUMES THAT ACTIVE SET IS CONSTANT FOR THE
   *  VARIANCE-COVARIANCE EVALUATION)
   * ----------------------------------------------- */
  int nFR, nFX, nAC;

  nFR = qp->getNFR( );
  nFX = qp->getNFX( );
  nAC = qp->getNAC( );


  /* ASK FOR THE CORRESPONDING INDEX ARRAYS:
   * --------------------------------------- */
  int *FR_idx, *FX_idx, *AC_idx;

  if ( qp->bounds.getFree( )->getNumberArray( &FR_idx ) != SUCCESSFUL_RETURN )
       return THROWERROR( RET_HOTSTART_FAILED );

  if ( qp->bounds.getFixed( )->getNumberArray( &FX_idx ) != SUCCESSFUL_RETURN )
       return THROWERROR( RET_HOTSTART_FAILED );

  if ( qp->constraints.getActive( )->getNumberArray( &AC_idx ) != SUCCESSFUL_RETURN )
       return THROWERROR( RET_HOTSTART_FAILED );



  /* INTRODUCE VARIABLES TO MEASURE THE REACTION OF THE QP-SOLUTION TO
   * THE VARIANCE-COVARIANCE DISTURBANCE:
   * ----------------------------------------------------------------- */
  real_t *delta_xFR = new real_t[nFR];
  real_t *delta_xFX = new real_t[nFX];
  real_t *delta_yAC = new real_t[nAC];
  real_t *delta_yFX = new real_t[nFX];

  real_t* K             = new real_t[dim*dim];  /* matrix to store */
                                                /* an intermediate */
                                                /* result. */

  /* SOME INITIALIZATIONS:
   * --------------------- */
  for( run1 = 0; run1 < dim*dim; run1++ ){
    K              [run1] = 0.0;
    Primal_Dual_VAR[run1] = 0.0;
  }


  /* ================================================================= */

  /* FIRST MATRIX MULTIPLICATION (OBTAINS THE INTERMEDIATE RESULT
   *  K := [ ("ACTIVE" KKT-MATRIX OF THE QP)^(-1) * g_b_bA_VAR ]^T )
   * THE EVALUATION OF THE INVERSE OF THE KKT-MATRIX OF THE QP
   * WITH RESPECT TO THE CURRENT ACTIVE SET
   * USES THE EXISTING CHOLESKY AND TQ-DECOMPOSITIONS. FOR DETAILS
   * cf. THE (protected) FUNCTION determineStepDirection. */

  for( run3 = 0; run3 < dim; run3++ ){


    for( run1 = 0; run1 < nV; run1++ ){
      delta_g_cov  [run1]   = g_b_bA_VAR[run3*dim+run1];
      delta_lb_cov [run1]   = g_b_bA_VAR[run3*dim+nV+run1];         /*  LINE-WISE LOADING OF THE INPUT */
      delta_ub_cov [run1]   = g_b_bA_VAR[run3*dim+nV+run1];         /*  VARIANCE-COVARIANCE            */
    }
    for( run1 = 0; run1 < nC; run1++ ){
      delta_lbA_cov [run1]  = g_b_bA_VAR[run3*dim+2*nV+run1];
      delta_ubA_cov [run1]  = g_b_bA_VAR[run3*dim+2*nV+run1];
    }


    /* EVALUATION OF THE STEP:
     * ------------------------------------------------------------------------------ */

    returnvalue = qp->determineStepDirection( delta_g_cov, delta_lbA_cov, delta_ubA_cov, delta_lb_cov, delta_ub_cov,
                                              Delta_bC_isZero, Delta_bB_isZero, delta_xFX,delta_xFR,
                                              delta_yAC,delta_yFX );

    /* ------------------------------------------------------------------------------ */


    /* STOP THE ALGORITHM IN THE CASE OF NO SUCCESFUL RETURN:
     * ------------------------------------------------------ */
    if ( returnvalue != SUCCESSFUL_RETURN ){

      delete[] delta_g_cov;
      delete[] delta_lb_cov;
      delete[] delta_ub_cov;
      delete[] delta_lbA_cov;
      delete[] delta_ubA_cov;
      delete[] delta_xFR;
      delete[] delta_xFX;
      delete[] delta_yAC;
      delete[] delta_yFX;
      delete[] K;

      THROWERROR( RET_STEPDIRECTION_DETERMINATION_FAILED );
      return returnvalue;
    }



    for( run1=0; run1<nFR; run1++ ){
      run2                  = FR_idx[run1];
      K[run3*dim+run2]      = delta_xFR[run1];
    }                                                               /*  LINE WISE                  */
    for( run1=0; run1<nFX; run1++ ){                                /*  STORAGE OF THE QP-REACTION */
      run2                  = FX_idx[run1];                         /*  (uses the index list)      */
      K[run3*dim+run2]      = delta_xFX[run1];
      K[run3*dim+nV+run2]   = delta_yFX[run1];
    }
    for( run1=0; run1<nAC; run1++ ){
      run2                  = AC_idx[run1];
      K[run3*dim+2*nV+run2] = delta_yAC[run1];
    }

  }


  /* ================================================================= */

  /* SECOND MATRIX MULTIPLICATION (OBTAINS THE FINAL RESULT
   * Primal_Dual_VAR := ("ACTIVE" KKT-MATRIX OF THE QP)^(-1) * K )
   * THE APPLICATION OF THE KKT-INVERSE IS AGAIN REALIZED
   * BY USING THE PROTECTED FUNCTION
   * determineStepDirection */

  for( run3 = 0; run3 < dim; run3++ ){

    for( run1 = 0; run1 < nV; run1++ ){
      delta_g_cov  [run1]   = K[run3+     run1*dim];
      delta_lb_cov [run1]   = K[run3+(nV+run1)*dim];                /*  ROW WISE LOADING OF THE */
      delta_ub_cov [run1]   = K[run3+(nV+run1)*dim];                /*  INTERMEDIATE RESULT K   */
    }
    for( run1 = 0; run1 < nC; run1++ ){
      delta_lbA_cov [run1]  = K[run3+(2*nV+run1)*dim];
      delta_ubA_cov [run1]  = K[run3+(2*nV+run1)*dim];
    }


    /* EVALUATION OF THE STEP:
     * ------------------------------------------------------------------------------ */

    returnvalue = qp->determineStepDirection( delta_g_cov, delta_lbA_cov, delta_ubA_cov, delta_lb_cov, delta_ub_cov,
                                              Delta_bC_isZero, Delta_bB_isZero, delta_xFX,delta_xFR,
                                              delta_yAC,delta_yFX);


    /* ------------------------------------------------------------------------------ */


    /* STOP THE ALGORITHM IN THE CASE OF NO SUCCESFUL RETURN:
     * ------------------------------------------------------ */
    if ( returnvalue != SUCCESSFUL_RETURN ){

      delete[] delta_g_cov;
      delete[] delta_lb_cov;
      delete[] delta_ub_cov;
      delete[] delta_lbA_cov;
      delete[] delta_ubA_cov;
      delete[] delta_xFR;
      delete[] delta_xFX;
      delete[] delta_yAC;
      delete[] delta_yFX;
      delete[] K;

      THROWERROR( RET_STEPDIRECTION_DETERMINATION_FAILED );
      return returnvalue;
    }



    for( run1=0; run1<nFR; run1++ ){
      run2                                = FR_idx[run1];
      Primal_Dual_VAR[run3+run2*dim]      = delta_xFR[run1];
    }
    for( run1=0; run1<nFX; run1++ ){                                 /*  ROW-WISE STORAGE */
      run2                  = FX_idx[run1];                          /*  OF THE RESULT.   */
      Primal_Dual_VAR[run3+run2*dim     ]   = delta_xFX[run1];
      Primal_Dual_VAR[run3+(nV+run2)*dim]   = delta_yFX[run1];
    }
    for( run1=0; run1<nAC; run1++ ){
      run2                                  = AC_idx[run1];
      Primal_Dual_VAR[run3+(2*nV+run2)*dim] = delta_yAC[run1];
    }

  }


  /* DEALOCATE MEMORY:
   * ----------------- */

  delete[] delta_g_cov;
  delete[] delta_lb_cov;
  delete[] delta_ub_cov;
  delete[] delta_lbA_cov;
  delete[] delta_ubA_cov;
  delete[] delta_xFR;
  delete[] delta_xFX;
  delete[] delta_yAC;
  delete[] delta_yFX;
  delete[] K;

  return SUCCESSFUL_RETURN;
}


/*
 *	g e t V a r i a n c e C o v a r i a n c e
 */
returnValue SolutionAnalysis::getVarianceCovariance(	SQProblem* const qp,
														const real_t* const g_b_bA_VAR, real_t* const Primal_Dual_VAR
														) const
{
	/* Call QProblem variant. */
	return getVarianceCovariance( (QProblem*)qp,g_b_bA_VAR,Primal_Dual_VAR );
}


END_NAMESPACE_QPOASES


/*
 *	end of file
 */
