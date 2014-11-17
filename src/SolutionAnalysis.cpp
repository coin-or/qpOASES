/*
 *	This file is part of qpOASES.
 *
 *	qpOASES -- An Implementation of the Online Active Set Strategy.
 *	Copyright (C) 2007-2014 by Hans Joachim Ferreau, Andreas Potschka,
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
 *	\version 3.0
 *	\date 2008-2014
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
 *	g e t M a x K K T v i o l a t i o n
 */
returnValue SolutionAnalysis::getMaxKKTviolation( QProblemB* qp, real_t& maxKKTviolation ) const
{
	int i;
	int nV = qp->getNV( );

	real_t *tmp = new real_t[nV];
	maxKKTviolation = 0.0;


	/* 1) Check for Hx + g - y*A' = 0  (here: A = Id). */
	for( i=0; i<nV; ++i )
		tmp[i] = qp->g[i];

	switch ( qp->getHessianType( ) )
	{
		case HST_ZERO:
			/*tmp += qp->eps * qp->x[i]; */
			break;

		case HST_IDENTITY:
			for( i=0; i<nV; ++i )
				tmp[i] += qp->x[i];
			break;

		default:
			qp->H->times(1, 1.0, qp->x, nV, 1.0, tmp, nV);
			break;
	}

	for( i=0; i<nV; ++i )
	{
		tmp[i] -= qp->y[i];

		if ( getAbs( tmp[i] ) > maxKKTviolation )
			maxKKTviolation = getAbs( tmp[i] );
	}
	delete[] tmp;

	/* 2) Check for lb <= x <= ub. */
	for( i=0; i<nV; ++i )
	{
		if ( qp->lb[i] - qp->x[i] > maxKKTviolation )
			maxKKTviolation = qp->lb[i] - qp->x[i];

		if ( qp->x[i] - qp->ub[i] > maxKKTviolation )
			maxKKTviolation = qp->x[i] - qp->ub[i];
	}

	/* 3) Check for correct sign of y and for complementary slackness. */
	for( i=0; i<nV; ++i )
	{
		switch ( qp->bounds.getStatus( i ) )
		{
			case ST_LOWER:
				if ( -qp->y[i] > maxKKTviolation )
					maxKKTviolation = -qp->y[i];
				if ( getAbs( ( qp->x[i] - qp->lb[i] ) * qp->y[i] ) > maxKKTviolation )
					maxKKTviolation = getAbs( ( qp->x[i] - qp->lb[i] ) * qp->y[i] );
				break;

			case ST_UPPER:
				if ( qp->y[i] > maxKKTviolation )
					maxKKTviolation = qp->y[i];
				if ( getAbs( ( qp->ub[i] - qp->x[i] ) * qp->y[i] ) > maxKKTviolation )
					maxKKTviolation = getAbs( ( qp->ub[i] - qp->x[i] ) * qp->y[i] );
				break;

			default: /* inactive */
			if ( getAbs( qp->y[i] ) > maxKKTviolation )
					maxKKTviolation = getAbs( qp->y[i] );
				break;
		}
	}


	return SUCCESSFUL_RETURN;
}


/*
 *	g e t M a x K K T v i o l a t i o n
 */
returnValue SolutionAnalysis::getMaxKKTviolation( QProblem* qp, real_t& maxKKTviolation ) const
{
	int i;
	int nV = qp->getNV( );
	int nC = qp->getNC( );

	real_t *tmp = new real_t[nV];
	maxKKTviolation = 0.0;


	/* 1) check for Hx + g - [yFX yAC]*[Id A]' = 0. */
	for( i=0; i<nV; ++i )
		tmp[i] = qp->g[i];

	switch ( qp->getHessianType( ) )
	{
		case HST_ZERO:
				/*tmp += qp->eps * qp->x[i]; */
			break;

		case HST_IDENTITY:
			for( i=0; i<nV; ++i )
				tmp[i] += qp->x[i];
			break;

		default:
			qp->H->times(1, 1.0, qp->x, nV, 1.0, tmp, nV);
			break;
	}

	qp->A->transTimes(1, -1.0, qp->y + nV, nC, 1.0, tmp, nV);

	for( i=0; i<nV; ++i )
	{
		tmp[i] -= qp->y[i];

		if ( getAbs( tmp[i] ) > maxKKTviolation )
			maxKKTviolation = getAbs( tmp[i] );
	}

	/* 2) Check for [lb lbA] <= [Id A]*x <= [ub ubA]. */
	/* lbA <= Ax <= ubA */
	real_t* Ax = new real_t[nC];
	qp->A->times(1, 1.0, qp->x, nV, 0.0, Ax, nC);

	for( i=0; i<nC; ++i )
	{
		if ( qp->lbA[i] - Ax[i] > maxKKTviolation )
			maxKKTviolation = qp->lbA[i] - Ax[i];

		if ( Ax[i] - qp->ubA[i] > maxKKTviolation )
			maxKKTviolation = Ax[i] - qp->ubA[i];
	}

	/* lb <= x <= ub */
	for( i=0; i<nV; ++i )
	{
		if ( qp->lb[i] - qp->x[i] > maxKKTviolation )
			maxKKTviolation = qp->lb[i] - qp->x[i];

		if ( qp->x[i] - qp->ub[i] > maxKKTviolation )
			maxKKTviolation = qp->x[i] - qp->ub[i];
	}

	/* 3) Check for correct sign of y and for complementary slackness. */
	/* bounds */
	for( i=0; i<nV; ++i )
	{
		switch ( qp->bounds.getStatus( i ) )
		{
			case ST_LOWER:
				if ( -qp->y[i] > maxKKTviolation )
					maxKKTviolation = -qp->y[i];
				if ( getAbs( ( qp->x[i] - qp->lb[i] ) * qp->y[i] ) > maxKKTviolation )
					maxKKTviolation = getAbs( ( qp->x[i] - qp->lb[i] ) * qp->y[i] );
				break;

			case ST_UPPER:
				if ( qp->y[i] > maxKKTviolation )
					maxKKTviolation = qp->y[i];
				if ( getAbs( ( qp->ub[i] - qp->x[i] ) * qp->y[i] ) > maxKKTviolation )
					maxKKTviolation = getAbs( ( qp->ub[i] - qp->x[i] ) * qp->y[i] );
				break;

			default: /* inactive */
			if ( getAbs( qp->y[i] ) > maxKKTviolation )
					maxKKTviolation = getAbs( qp->y[i] );
				break;
		}
	}

	/* constraints */
	for( i=0; i<nC; ++i )
	{
		switch ( qp->constraints.getStatus( i ) )
		{
			case ST_LOWER:
				if ( -qp->y[nV+i] > maxKKTviolation )
					maxKKTviolation = -qp->y[nV+i];
				if ( getAbs( ( Ax[i] - qp->lbA[i] ) * qp->y[nV+i] ) > maxKKTviolation )
					maxKKTviolation = getAbs( ( Ax[i] - qp->lbA[i] ) * qp->y[nV+i] );
				break;

			case ST_UPPER:
				if ( qp->y[nV+i] > maxKKTviolation )
					maxKKTviolation = qp->y[nV+i];
				if ( getAbs( ( qp->ubA[i] - Ax[i] ) * qp->y[nV+i] ) > maxKKTviolation )
					maxKKTviolation = getAbs( ( qp->ubA[i] - Ax[i] ) * qp->y[nV+i] );
				break;

			default: /* inactive */
			if ( getAbs( qp->y[nV+i] ) > maxKKTviolation )
					maxKKTviolation = getAbs( qp->y[nV+i] );
				break;
		}
	}

	delete[] tmp;
	delete[] Ax;

	return SUCCESSFUL_RETURN;
}


/*
 *	g e t M a x K K T v i o l a t i o n
 */
returnValue SolutionAnalysis::getMaxKKTviolation( SQProblem* qp, real_t& maxKKTviolation ) const
{
	return getMaxKKTviolation( (QProblem*)qp,maxKKTviolation );
}



/*
 *	g e t V a r i a n c e C o v a r i a n c e
 */
returnValue SolutionAnalysis::getVarianceCovariance( QProblemB* qp, real_t* g_b_bA_VAR, real_t* Primal_Dual_VAR ) const
{
	return THROWERROR( RET_NOT_YET_IMPLEMENTED );
}


/*
 *	g e t V a r i a n c e C o v a r i a n c e
 */
returnValue SolutionAnalysis::getVarianceCovariance( QProblem* qp, real_t* g_b_bA_VAR, real_t* Primal_Dual_VAR ) const
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
returnValue SolutionAnalysis::getVarianceCovariance( SQProblem* qp, real_t* g_b_bA_VAR, real_t* Primal_Dual_VAR ) const
{
	/* Call QProblem variant. */
	return getVarianceCovariance( (QProblem*)qp,g_b_bA_VAR,Primal_Dual_VAR );
}


/*
 *	c h e c k C u r v a t u r e O n S e t S
 */
int SolutionAnalysis::checkCurvatureOnStronglyActiveConstraints( SQProblem* qp )
{
  printf("checkCurvatureOnStronglyActiveConstraints( SQProblem* qp ) not yet implemented for standard qpOASES!\n");
  return 1;
}


/*
 *	c h e c k C u r v a t u r e O n S t r o n g l y A c t i v e C o n s t r a i n t s
 */
int SolutionAnalysis::checkCurvatureOnStronglyActiveConstraints( SQProblemSchur* qp )
{
  real_t eps = 1.0e-16;
  returnValue ret;
  Bounds saveBounds;
  QProblemStatus saveStatus;
  int k, fail, neig, nAC, nFX, *FX_idx;

  nFX = qp->getNFX( );
  nAC = qp->getNAC( );

  // If no bounds are active reduced Hessian is positive definite (otherwise qpOASES wouldnt have finished)
  if( nFX == 0 )
    return 0;

  // Get active bounds (deep copy)
  qp->getBounds( saveBounds );
  saveBounds.getFixed( )->getNumberArray( &FX_idx );

  // We have to change the status to modify the active set
  saveStatus = qp->getStatus();
  qp->status = QPS_PERFORMINGHOMOTOPY;

  // If a variable is active now but has not been in the previous major iteration remove it
  for( k=0; k<nFX; k++ )
    if( getAbs(qp->x[FX_idx[k]]) > eps )
      if ( qp->bounds.moveFixedToFree( FX_idx[k] ) != SUCCESSFUL_RETURN )
	return THROWERROR( RET_REMOVEBOUND_FAILED );

  // Do a new factorization and check the inertia
  ret = qp->resetSchurComplement( BT_FALSE );
  neig = qp->sparseSolver->getNegativeEigenvalues( );
  if( ret != SUCCESSFUL_RETURN || neig != nAC )
    fail = 1;
  else
    fail = 0;

  // Add all bounds that have been removed
  for( k=0; k<nFX; k++ )
    if( qp->bounds.getStatus( FX_idx[k] ) == ST_INACTIVE )
      qp->bounds.moveFreeToFixed( FX_idx[k], saveBounds.getStatus( FX_idx[k] ) );

  qp->status = saveStatus;
  return fail;
}


//int SolutionAnalysis::checkCurvatureOnStronglyActiveConstraints( SQProblemSchur* qp )
//{
  //real_t eps = 1.0e-16;
  //real_t oldDet, newDet;
  //int oldNS;
  //returnValue ret;
  //Bounds saveBounds;
  //QProblemStatus saveStatus;
  //int nFX, *FX_idx;
  //int k, fail, neig, rmCnt, nAC;

  //// Get active bounds (deep copy)
  //nFX = qp->getNFX( );
  //nAC = qp->getNAC( );
  //qp->getBounds( saveBounds );
  //saveBounds.getFixed( )->getNumberArray( &FX_idx );

  //// If no bounds are active reduced Hessian is positive definite (otherwise qpOASES wouldnt have finished)
  //if( nFX == 0 )
    //return 0;

  //// We have to modify the status to call removeBound()
  //saveStatus = qp->getStatus();
  //qp->status = QPS_PERFORMINGHOMOTOPY;

  //// If a variable is active but was not active in the previous major iteration
  //// remove it to see if a negative eigenvalue appears
  //rmCnt = 0;
  //fail = 0;
  //for( k=0; k<nFX; k++ )
    //if( getAbs(qp->x[FX_idx[k]]) > eps )
    //{
      //oldDet = qp->detS;
      //oldNS = qp->nS;

      //ret = qp->removeBound( FX_idx[k], BT_TRUE, BT_FALSE, BT_FALSE );
      //if( ret != SUCCESSFUL_RETURN )
      //{
	//fail = 1;
	//break;
      //}

      //newDet = qp->detS;
      //rmCnt++;

      //// Case 1: S has grown by 1 row and column
      //if( qp->nS == oldNS + 1 )
      //{
	//// If the determinant does not change sign, then S has gained a positive eigenvalue.
	//// That means there is a negative eigenvalue in the (extended) reduced Hessian!
	//if ( ( oldDet <= 0.0 && newDet <= 0.0 ) || ( oldDet >= 0.0 && newDet >= 0.0 ) )
	//{
	  //fail = 1;
	  //break;
	//}
      //}
      //// Case 2: S has shrunk by 1 row and column
      //else if( qp->nS == oldNS - 1 )
      //{
	//// If the determinant changes sign, then S has lost a negative eigenvalue.
	//// That means there is a negative eigenvalue in the (extended) reduced Hessian!
	//if ( ( oldDet <= 0.0 && newDet > 0.0 ) || ( oldDet >= 0.0 && newDet < 0.0 ) )
	//{
	  //fail = 1;
	  //break;
	//}
      //}
      //// Case 3: S was reset
      //else if( qp->nS == 0 )
      //{
	//// Check inertia of KKT matrix
	//neig = qp->sparseSolver->getNegativeEigenvalues( );
	//if( neig > nAC )
	//{
	  //fail = 1;
	  //break;
	//}
      //}
      //else
	//printf("ERROR!\n");
    //}

  //// If test is successful, add all bounds that have been removed
  //// If not, don't bother with that because we will discard this QP object anyway
  //if( fail == 0 )
    //for( k=0; k<nFX; k++ )
    //{
      //ret = qp->addBound( FX_idx[k], saveBounds.getStatus( FX_idx[k] ), BT_TRUE, BT_FALSE );
      //if( ret != SUCCESSFUL_RETURN && ret != RET_BOUND_ALREADY_ACTIVE )
	//printf( "addBound() in checkCurvatureOnStronglyActiveConstraints(): %s\n", getGlobalMessageHandler()->getErrorCodeMessage( ret ) );
    //}

  //qp->status = saveStatus;
  //return fail;
//}


END_NAMESPACE_QPOASES


/*
 *	end of file
 */
