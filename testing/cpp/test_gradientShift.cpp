/*
 *	This file is part of qpOASES.
 *
 *	qpOASES -- An Implementation of the Online Active Set Strategy.
 *	Copyright (C) 2007-2009 by Hans Joachim Ferreau et al. All rights reserved.
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
 *	\file testing/cpp/test_gradientShift.cpp
 *	\author Hans Joachim Ferreau
 *	\version 3.1
 *	\date 2007-2015
 *
 *	Simple test case which caused troubles in version 2.0.
 */

 

#include <math.h>
#include <iostream>
#include <qpOASES.hpp>
#include <qpOASES/UnitTesting.hpp>


/** Running simple test case which caused troubles in version 2.0. */
int main( )
{
	USING_NAMESPACE_QPOASES

	
	double H[2*2] 	= { 0.055944055944055944,0, 0, 0 };
	double A[1*2] 	= { 0.70514808036997589, -1 };
	//double g[2]	= { -15.830543073928741, 0};
	double g[2]		= { 0, 0 };
	double lb[2]	= { 137.00242299940646, 154.0 };
	double ub[2]	= { 282.19008595111382, 198.98579740786641};
	double lbA[1]	= {0.0};
	double ubA[1]	= {0.0};
	
	double gStart=-16.4;
 
	/* Setting up QProblem object. */
	QProblem example( 2,1 );
	Options options;
	options.setToMPC();
	options.printLevel = REFER_NAMESPACE_QPOASES PL_NONE;
	example.setOptions( options );
	
	returnValue ret;
	int nWSR = 10;
	double Xopt[2]={0.0,0.0};
	//fprintf(stdFile, "g[0]\t,\tReturn code\n");
	int errorCount=0;
	int i=0;
	double granularity=0.00001;
	g[0]=gStart;
	for( i=0; i< 70000; i++)
	{

		g[0]=g[0]+granularity;
		/* Solve first QP. */
		nWSR = 10;
		ret= example.init( H,g,A,lb,ub,lbA,ubA, nWSR,0 );
		if (ret != SUCCESSFUL_RETURN)
		{
			//fprintf(stdFile, "%f\t,\t%d\n",g[0],ret);
			errorCount++;
		}
		//fprintf(stdFile, "%f\t,\t%d\n",g[0],ret);
	}
	example.printProperties();
	fprintf(stdFile, "#Number of optimizer runs: %d\n",i);
	fprintf(stdFile, "#g[0] test interval: %f < g[0] < %f\n",gStart,g[0]);
	fprintf(stdFile, "#Granularity: %f\n",granularity);
	double errorPercent = double(errorCount)/double(i)*100.0;
	fprintf(stdFile, "#Number of errors (error): %d (%f)\n",errorCount,errorPercent);

	example.getPrimalSolution(Xopt);
	fprintf(stdFile,"#Optimization primary result : LD=%f BD=%f\n",Xopt[0], Xopt[1]);

	double Yopt[3]={0.0,0.0,0.0};
	example.getDualSolution(Yopt);
	fprintf(stdFile,"#Optimization dual result : %f %f %f8\n",Yopt[0], Yopt[1], Yopt[2]);

	int Nc=0;
	Nc=example.getNC();
	fprintf(stdFile,"#Number of constraints : %d\n",Nc);

	int Nec=0;
	Nec=example.getNEC();
	fprintf(stdFile,"#Number of equality constraints : %d\n",Nec);

	int Nac=0;
	Nac=example.getNAC();
	fprintf(stdFile,"#Number of active constraints : %d\n",Nac);

	int Niac=0;
	Niac=example.getNIAC();
	fprintf(stdFile,"#Number of inactive constraints : %d\n",Niac);

	QPOASES_TEST_FOR_TRUE( errorCount == 0 )

	return TEST_PASSED;
}


/*
 *	end of file
 */
