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
 *	\file src/Presolver.cpp
 *	\author Hans Joachim Ferreau, Andreas Potschka, Christian Kirches, Dominik Cebulla
 *	\version 3.2
 *	\date 2007-2015
 *
 *	Implementation of the Presolver class designed for preprocessing a QP.
 */


#include <qpOASES/Presolver.hpp>


BEGIN_NAMESPACE_QPOASES


Presolver::Presolver(const int_t nV,
					 const int_t nC,
					 const int_t nzH,
					 const int_t nzA,
					 const int_t presolveStackSize) : data(0)
{
	qpp_return_value_t err = qppNew(&data, nC, nV, nzA, nzH, presolveStackSize);
	if (err != QPP_OK)
	{
		THROWERROR( convertErrorCode(err) );
	}
}


Presolver::~Presolver()
{
	qppFree(&data);
}


returnValue Presolver::presolveCoordMat(int_t* nV,
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
                                        const MatrixSortType& sortType)
{
	qpp_return_value_t err = qppInit(data, Airn, Ajcn, Ax, *nzA, Hirn, Hjcn, Hx, *nzH, g,
                                 0.0, lb, ub, lbA, ubA);
    if (err != QPP_OK)
    {
        return THROWERROR( convertErrorCode(err) );
    }

    err = qppPresolve(data);
    if (err != QPP_OK)
    {
        return THROWERROR( convertErrorCode(err) );
    }

    real_t pf;
    err = qppGetPresolvedQP(data, nC, nV, Airn, Ajcn, Ax, nzA, Hirn, Hjcn, Hx,
                            nzH, g, &pf, lb, ub, lbA, ubA, sortType);
    if (err != QPP_OK)
    {
        return THROWERROR( convertErrorCode(err) );
    }
    return SUCCESSFUL_RETURN;
}


returnValue Presolver::presolve(int_t* nV,
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
                                real_t* const ubA)
{
    int_t* Hjcn = 0;
    int_t* Ajcn = 0;
    int_t nVar = data->n;
    int_t nCon = data->m;

    /* Converting column pointers from CSC format into row subscripts for coordinate format. */
    if ( (Hcp != 0) && (*nzH > 0) )
    {
        Hjcn = new int_t[*nzH];
        for (int_t j = 0; j < nVar; ++j)
            for (int_t i = Hcp[j]; i < Hcp[j+1]; ++i)
                Hjcn[i] = j;
    }

    if ( (Acp != 0) && (*nzA > 0) )
    {
        Ajcn = new int_t[*nzA];
        for (int_t j = 0; j < nVar; ++j)
            for (int_t i = Acp[j]; i < Acp[j+1]; ++i)
                Ajcn[i] = j;
    }

    /* PRESOLVING: */
    returnValue retval = presolveCoordMat(&nVar, &nCon, Hirn, Hjcn, Hx, nzH, g, Airn,
                                          Ajcn, Ax, nzA, lb, ub, lbA, ubA, MST_COLUMN_WISE);
    if (retval != SUCCESSFUL_RETURN)
    {
        delete[] Hjcn;
        delete[] Ajcn;
        return THROWERROR( retval );
    }

    /* Compute column pointers for CSC format from column subscripts of coordinate format. */
    if (Hcp != 0)
    {
        for (int_t i = 0; i <= nVar; ++i)
            Hcp[i] = 0;

        for (int_t i = 0; i < *nzH; ++i)
            Hcp[Hjcn[i]+1]++;

        for (int_t i = 1; i <= nVar; ++i)
            Hcp[i] += Hcp[i-1];
    }

    if (Acp != 0)
    {
        for (int_t i = 0; i <= nVar; ++i)
            Acp[i] = 0;

        for (int_t i = 0; i < *nzA; ++i)
            Acp[Ajcn[i]+1]++;

        for (int_t i = 1; i <= nVar; ++i)
            Acp[i] += Acp[i-1];
    }

    delete[] Hjcn;
    delete[] Ajcn;

    *nV = nVar;
    *nC = nCon;

    return SUCCESSFUL_RETURN;
}


returnValue Presolver::presolve(int_t* nV,
                                int_t* nC,
                                const SymSparseMat* const H,
                                SymSparseMat** Hprs,
                                real_t* const g,
                                const SparseMatrix* const A,
                                SparseMatrix** Aprs,
                                real_t* const lb,
                                real_t* const ub,
                                real_t* const lbA,
                                real_t* const ubA)
{
    int_t nzA = 0;
    int_t nzH = 0;
    int_t nVar = data->n;
    int_t nCon = data->m;
    returnValue retval = SUCCESSFUL_RETURN;

    int_t* ri = 0; int_t* ci = 0;
    int_t* Airn = 0; int_t* Ajcn = 0; real_t* Ax = 0;
    int_t* Hirn = 0; int_t* Hjcn = 0; real_t* Hx = 0;
    *Aprs = 0;
    *Hprs = 0;

    ri = new int_t[nCon];     /* = 1:nCon */
    ci = new int_t[nVar];     /* = 1:nVar */

    for (int_t i = 0; i < nCon; ++i)
        ri[i] = i;

    for (int_t i = 0; i < nVar; ++i)
        ci[i] = i;

    /* Converting the constraint matrix into coordinate format: */
    if (A != 0)
    {
        /* First get number of nonzero elements: */
        retval = A->getSparseSubmatrix(nCon, ri, nVar, ci, 0, 0, nzA, 0, 0, 0);
        if (retval != SUCCESSFUL_RETURN)
        {
            goto PRS_EXIT;
        }

        /* Memory for storing the constraint matrix in coordinate or CSC format. */
        Airn = new int_t[nzA];
        Ajcn = new int_t[nzA < nVar+1 ? nVar+1 : nzA];
        Ax = new real_t[nzA];

        retval = A->getSparseSubmatrix(nCon, ri, nVar, ci, 0, 0, nzA, Airn, Ajcn, Ax);
        if (retval != SUCCESSFUL_RETURN)
        {
            goto PRS_EXIT;
        }
    }

    /* Converting the Hessian matrix into coordinate format: */
    if (H != 0)
    {
        /* First get the number of nonzero elements of the lower triangular part of the Hessian mat. */
        retval = H->getSparseSubmatrix(nVar, ci, nVar, ci, 0, 0, nzH, 0, 0, 0, BT_TRUE);
        if (retval != SUCCESSFUL_RETURN)
        {
            goto PRS_EXIT;
        }

        /* Memory for storing the lower triangular part of the Hessian matrix in coordinate or CSC format. */
        Hirn = new int_t[nzH];
        Hjcn = new int_t[nzH < nVar+1 ? nVar+1 : nzH];
        Hx = new real_t[nzH];

        retval = H->getSparseSubmatrix(nVar, ci, nVar, ci, 0, 0, nzH, Hirn, Hjcn, Hx, BT_TRUE);
        if (retval != SUCCESSFUL_RETURN)
        {
            goto PRS_EXIT;
        }
    }

    /* PRESOLVING */
    retval = presolveCoordMat(&nVar, &nCon, Hirn, Hjcn, Hx, &nzH, g, Airn, Ajcn, Ax, &nzA,
                              lb, ub, lbA, ubA, MST_COLUMN_WISE);
    if (retval != SUCCESSFUL_RETURN)
    {
        goto PRS_EXIT;
    }

    /*  Convert constraint matrix back into CSC format. A new SparseMatrix will be created if
        matrix A was not a null pointer. */

    /* Compute column pointers for CSC format: */
    if (A != 0)
    {
        for (int_t i = 0; i < nVar; ++i)
            ci[i] = 0;

        for (int_t i = 0; i < nzA; ++i)
            ci[Ajcn[i]]++;

        Ajcn[0] = 0;
        for (int_t i = 1; i < nVar; ++i)
        {
            Ajcn[i] = ci[i-1] + Ajcn[i-1];
        }
        Ajcn[nVar] = nzA;

        *Aprs = new SparseMatrix(nCon, nVar, Airn, Ajcn, Ax);
        (*Aprs)->doFreeMemory();  /* Deep copy! */
    }


    /*  Convert Hessian matrix back into CSC format. A new SymSparseMat will be created
        if matrix H was not a null pointer. Note that only the lower triangular part of
        the Hessian matrix will be stored! */

    /* Compute column pointers for CSC format: */
    if (H != 0)
    {
        for (int_t i = 0; i < nVar; ++i)
            ci[i] = 0;

        for (int_t i = 0; i < nzH; ++i)
            ci[Hjcn[i]]++;

        Hjcn[0] = 0;
        for (int_t i = 1; i < nVar; ++i)
        {
            Hjcn[i] = ci[i-1] + Hjcn[i-1];
        }
        Hjcn[nVar] = nzH;

        *Hprs = new SymSparseMat(nVar, nVar, Hirn, Hjcn, Hx);
        (*Hprs)->doFreeMemory();      /* Deep copy! */
    }

	delete[] ri;
	delete[] ci;

	*nV = nVar;
	*nC = nCon;

	return SUCCESSFUL_RETURN;

PRS_EXIT:
    delete[] ci;
    delete[] ri;
    delete[] Airn;
    delete[] Ajcn;
    delete[] Ax;
    delete[] Hirn;
    delete[] Hjcn;
    delete[] Hx;

    return THROWERROR( retval );
}


returnValue Presolver::presolve(int_t* nV,
                                int_t* nC,
                                real_t* const H,
                                real_t* const g,
                                real_t* const A,
                                real_t* const lb,
                                real_t* const ub,
                                real_t* const lbA,
                                real_t* const ubA)
{
    int_t nVar = data->n;
    int_t nCon = data->m;

    /* Dense matrices (row-wise stored) must be converted into Coordinate format */
    int_t nzH = 0;
    int_t nzA = 0;
    int_t* Airn = 0; int_t* Ajcn = 0; real_t* Ax = 0;
    int_t* Hirn = 0; int_t* Hjcn = 0; real_t* Hx = 0;

    int_t it = 0;

    /* Converting dense constraint matrix into coordinate format: */
    if (A != 0)
    {
        for (int_t i = 0; i < nCon*nVar; ++i)
        {
            if (fabs(A[i]) > ZERO)
                ++nzA;
        }

        Airn = new int_t[nzA];
        Ajcn = new int_t[nzA];
        Ax = new real_t[nzA];

        it = 0;
        for (int_t i = 0; i < nCon; ++i)
        {
            for (int_t j = 0; j < nVar; ++j)
            {
                if (fabs(A[i*nVar+j]) > ZERO)
                {
                    Airn[it] = i;
                    Ajcn[it] = j;
                    Ax[it] = A[i*nVar+j];
                    ++it;
                }
            }
        }
    }

    /* Converting dense Hessian matrix into coordinate format (only lower triangular part): */
    if (H != 0)
    {
        for (int_t i = 0; i < nVar; ++i)
        {
            for (int_t j = 0; j < nVar; ++j)
            {
                if ( (fabs(H[i*nVar+j]) > ZERO) && (i >= j) )
                    ++nzH;
            }
        }

        Hirn = new int_t[nzH];
        Hjcn = new int_t[nzH];
        Hx = new real_t[nzH];

        it = 0;
        for (int_t i = 0; i < nVar; ++i)
        {
            for (int_t j = 0; j < nVar; ++j)
            {
                if ( (fabs(H[i*nVar+j]) > ZERO) && (i >= j) )
                {
                    Hirn[it] = i;
                    Hjcn[it] = j;
                    Hx[it] = H[i*nVar+j];
                    ++it;
                }
            }
        }
    }

    /* PRESOLVING */
    returnValue retval = presolveCoordMat(&nVar, &nCon, Hirn, Hjcn, Hx, &nzH, g, Airn, Ajcn,
                                          Ax, &nzA, lb, ub, lbA, ubA, MST_ROW_WISE);
    if (retval != SUCCESSFUL_RETURN)
    {
        goto PRS_EXIT;
    }

    /* Converting presolved constraint matrix from coordinate to dense format: */
    if (A != 0)
    {
        for (int_t i = 0; i < nCon*nVar; ++i)
            A[i] = 0.0;

        for (int_t i = 0; i < nzA; ++i)
            A[Airn[i]*nVar + Ajcn[i]] = Ax[i];
    }

    /*  Converting presolved Hessian matrix from coordinate to dense format
        (WITH upper triangular part): */
    if (H != 0)
    {
        for (int_t i = 0; i < nVar*nVar; ++i)
            H[i] = 0.0;

        for (int_t i = 0; i < nzH; ++i)
        {
            H[Hirn[i]*nVar + Hjcn[i]] = Hx[i];
            H[Hjcn[i]*nVar + Hirn[i]] = Hx[i];
        }
    }

    delete[] Airn;
    delete[] Ajcn;
    delete[] Ax;
    delete[] Hirn;
    delete[] Hjcn;
    delete[] Hx;

    *nC = nCon;
    *nV = nVar;

    return SUCCESSFUL_RETURN;

PRS_EXIT:
    delete[] Airn;
    delete[] Ajcn;
    delete[] Ax;
    delete[] Hirn;
    delete[] Hjcn;
    delete[] Hx;

    return THROWERROR( retval );
}


returnValue Presolver::presolve(int_t* nV,
                                int_t* nC,
                                const SymDenseMat* const H,
                                SymDenseMat** Hprs,
                                real_t* const g,
                                const DenseMatrix* const A,
                                DenseMatrix** Aprs,
                                real_t* const lb,
                                real_t* const ub,
                                real_t* const lbA,
                                real_t* const ubA)
{
    int_t nzA = 0;
    int_t nzH = 0;
    int_t nCon = data->m;
    int_t nVar = data->n;

    int_t* ri = 0; int_t* ci = 0;
    int_t* Airn = 0; int_t* Ajcn = 0; real_t* Ax = 0;
    int_t* Hirn = 0; int_t* Hjcn = 0; real_t* Hx = 0;
    real_t* Hval = 0;
    real_t* Aval = 0;

    *Hprs = 0;
    *Aprs = 0;

    returnValue retval = SUCCESSFUL_RETURN;

    ri = new int_t[nCon];     /* = 1:nCon */
    ci = new int_t[nVar];     /* = 1:nVar */

    for (int_t i = 0; i < nCon; ++i)
        ri[i] = i;

    for (int_t i = 0; i < nVar; ++i)
        ci[i] = i;

    /* Converting the constraint matrix into coordinate format: */
    if (A != 0)
    {
        /* First get number of nonzero elements of the constraint matrix: */
        retval = A->getSparseSubmatrix(nCon, ri, nVar, ci, 0, 0, nzA, 0, 0, 0);
        if (retval != SUCCESSFUL_RETURN)
        {
            goto PRS_EXIT;
        }

        /* Memory for storing the constraint matrix in coordinate format. */
        Airn = new int_t[nzA];
        Ajcn = new int_t[nzA];
        Ax = new real_t[nzA];

        retval = A->getSparseSubmatrix(nCon, ri, nVar, ci, 0, 0, nzA, Airn, Ajcn, Ax);
        if (retval != SUCCESSFUL_RETURN)
        {
            goto PRS_EXIT;
        }
    }

    /* Converting the lower triangular part of the Hessian matrix into coordinate format: */
    if (H != 0)
    {
        /* First get number of nonzero elements of the lower triangular part of the Hessian matrix: */
        retval = H->getSparseSubmatrix(nVar, ci, nVar, ci, 0, 0, nzH, 0, 0, 0, BT_TRUE);
        if (retval != SUCCESSFUL_RETURN)
        {
            goto PRS_EXIT;
        }

        /* Memory for storing the Hessian matrix in coordinate format. */
        Hirn = new int_t[nzH];
        Hjcn = new int_t[nzH];
        Hx = new real_t[nzH];

        retval = H->getSparseSubmatrix(nVar, ci, nVar, ci, 0, 0, nzH, Hirn, Hjcn, Hx, BT_TRUE);
        if (retval != SUCCESSFUL_RETURN)
        {
            goto PRS_EXIT;
        }
    }

    /* PRESOLVING */
    retval = presolveCoordMat(&nVar, &nCon, Hirn, Hjcn, Hx, &nzH, g, Airn, Ajcn, Ax, &nzA,
                              lb, ub, lbA, ubA, MST_ROW_WISE);
    if (retval != SUCCESSFUL_RETURN)
    {
        goto PRS_EXIT;
    }

    /*  Convert constraint matrix back into dense format. A new DenseMatrix will be created if
        matrix A is not a null pointer. */
    if (A != 0)
    {
        Aval = new real_t[nCon*nVar];

        for (int_t i = 0; i < nCon*nVar; ++i)
            Aval[i] = 0.0;

        for (int_t i = 0; i < nzA; ++i)
        {
            Aval[Airn[i]*nVar + Ajcn[i]] = Ax[i];
        }

        *Aprs = new DenseMatrix(nCon, nVar, nVar, Aval);
        (*Aprs)->doFreeMemory();  /* Deep copy! */
    }

    /*  Convert Hessian matrix back into dense format. A new SymDenseMat will be created if matrix
        H is not a null pointer. We store the full Hessian matrix (i.e. lower AND upper triangular)
        part. */
    if (H != 0)
    {
        Hval = new real_t[nVar*nVar];

        for (int_t i = 0; i < nVar*nVar; ++i)
            Hval[i] = 0.0;

        for (int_t i = 0; i < nzH; ++i)
        {
            Hval[Hirn[i]*nVar + Hjcn[i]] = Hx[i];
            Hval[Hjcn[i]*nVar + Hirn[i]] = Hx[i];
        }

        *Hprs = new SymDenseMat(nVar, nVar, nVar, Hval);
        (*Hprs)->doFreeMemory();      /* Deep copy! */
    }

	delete[] ri;
	delete[] ci;
	delete[] Airn;
    delete[] Ajcn;
    delete[] Ax;
    delete[] Hirn;
    delete[] Hjcn;
    delete[] Hx;

    *nV = nVar;
    *nC = nCon;

	return SUCCESSFUL_RETURN;

PRS_EXIT:
    delete[] ci;
    delete[] ri;
    delete[] Airn;
    delete[] Ajcn;
    delete[] Ax;
    delete[] Hirn;
    delete[] Hjcn;
    delete[] Hx;

    return THROWERROR( retval );
}


returnValue Presolver::postsolve(real_t* const x,
                                 real_t* const y,
                                 Bounds* const bounds,
                                 Constraints* const constraints)
{
    returnValue retval = SUCCESSFUL_RETURN;
    int_t* wi = 0;
    int_t* wj = 0;
    int_t nVar = data->n;
    int_t pnVar = data->n_ps;
    int_t nCon = data->m;
    int_t pnCon = data->m_ps;

    if (pnCon > 0)
    {
        memmove(&y[nVar], &y[pnVar], static_cast<size_t>(pnCon) * sizeof(real_t));
    }

    if (bounds != 0)
    {
        wj = new int_t[nVar];

        /* Extract data from bounds, set entries 1:pnVar */
        for (int_t i = 0; i < pnVar; ++i)
        {
            switch (bounds->getStatus(i))
            {
            case ST_LOWER:
                wj[i] = QPP_AT_LOWER_BOUND;
                break;

            case ST_UPPER:
                wj[i] = QPP_AT_UPPER_BOUND;
                break;

            case ST_INACTIVE:
                wj[i] = QPP_AT_INACTIVE;
                break;

            default:
                /* Other status is currently not allowed. */
                delete[] wj;
                return RET_INVALID_ARGUMENTS;
            }
        }
    }

    if (constraints != 0)
    {
        wi = new int_t[nCon];
        /* Extract data from constraints */
        for (int_t i = 0; i < pnCon; ++i)
        {
            switch (constraints->getStatus(i))
            {
            case ST_LOWER:
                wi[i] = QPP_AT_LOWER_BOUND;
                break;

            case ST_UPPER:
                wi[i] = QPP_AT_UPPER_BOUND;
                break;

            case ST_INACTIVE:
                wi[i] = QPP_AT_INACTIVE;
                break;

            default:
                /* Other status is currently not allowed. */
                delete[] wi;
                delete[] wj;
                return RET_INVALID_ARGUMENTS;
            }
        }
    }

	qpp_return_value_t err = qppPostsolve(data, x, &y[nVar], y, wi, wj, 0);

	if (err != QPP_OK)
	{
		return THROWERROR( convertErrorCode(err) );
	}

	if (bounds != 0)
    {
        /* Initialize bounds with new dimension and working set */
        retval = bounds->init(nVar);
        if (retval != SUCCESSFUL_RETURN)
        {
            return THROWERROR( retval );
        }

        for (int_t i = 0; i < nVar; ++i)
        {
            switch (wj[i])
            {
            case QPP_AT_UPPER_BOUND:
                retval = bounds->setupBound(i, ST_UPPER);
                if (retval != SUCCESSFUL_RETURN)
                {
                    return THROWERROR( retval );
                }
                break;

            case QPP_AT_LOWER_BOUND:
                retval = bounds->setupBound(i, ST_LOWER);
                if (retval != SUCCESSFUL_RETURN)
                {
                    return THROWERROR( retval );
                }
                break;

            case QPP_AT_INACTIVE:
                retval = bounds->setupBound(i, ST_INACTIVE);
                if (retval != SUCCESSFUL_RETURN)
                {
                    return THROWERROR( retval );
                }
                break;

            default:
                /* Other status is not allowed. */
                delete[] wi;
                delete[] wj;
                return THROWERROR( RET_INVALID_ARGUMENTS );
            }
        }
    }

    if (constraints != 0)
    {
        /* Initialize constraints with new dimension and working set */
        retval = constraints->init(nCon);
        if (retval != SUCCESSFUL_RETURN)
        {
            return THROWERROR( retval );
        }

        for (int_t i = 0; i < nCon; ++i)
        {
            switch (wi[i])
            {
            case QPP_AT_UPPER_BOUND:
                retval = constraints->setupConstraint(i, ST_UPPER);
                if (retval != SUCCESSFUL_RETURN)
                {
                    return THROWERROR( retval );
                }
                break;

            case QPP_AT_LOWER_BOUND:
                retval = constraints->setupConstraint(i, ST_LOWER);
                if (retval != SUCCESSFUL_RETURN)
                {
                    return THROWERROR( retval );
                }
                break;

            case QPP_AT_INACTIVE:
                retval = constraints->setupConstraint(i, ST_INACTIVE);
                if (retval != SUCCESSFUL_RETURN)
                {
                    return THROWERROR( retval );
                }
                break;

            default:
                /* Other status is not allowed. */
                delete[] wi;
                delete[] wj;
                return THROWERROR( RET_INVALID_ARGUMENTS );
            }
        }
    }

	delete[] wi;
	delete[] wj;

	return SUCCESSFUL_RETURN;
}


returnValue Presolver::printOptions() const
{
    return ( PresolverOptions(*(static_cast<PresolverOptions*>(data->options))).print() );
}


returnValue Presolver::convertErrorCode(const qpp_return_value_t err) const
{
	switch (err)
    {
    case QPP_OK:
        return SUCCESSFUL_RETURN;

    case QPP_INVALID_ARGUMENT:
        return RET_INVALID_ARGUMENTS;

    case QPP_NULL_ARGUMENT:
        return RET_INVALID_ARGUMENTS;

    case QPP_INFEASIBLE:
        return RET_QP_INFEASIBLE;

    case QPP_PRIMAL_INFEASIBLE:
        return RET_QP_INFEASIBLE;

    case QPP_DUAL_INFEASIBLE:
        return RET_QP_INFEASIBLE;

    case QPP_BOUNDS_INCOMPATIBLE:
        return RET_QP_INFEASIBLE;

    case QPP_UNBOUNDED:
        return RET_QP_UNBOUNDED;

    case QPP_INVALID_OPTIONS:
        return RET_OPTIONS_ADJUSTED;

    case QPP_UNABLE_TO_OPEN_FILE:
        return RET_UNABLE_TO_OPEN_FILE;

    case QPP_UNABLE_TO_READ_FILE:
        return RET_UNABLE_TO_READ_FILE;

    case QPP_UNABLE_TO_WRITE_FILE:
        return RET_UNABLE_TO_WRITE_FILE;

    default:
        return RET_ERROR_UNDEFINED;
    }
}


END_NAMESPACE_QPOASES
