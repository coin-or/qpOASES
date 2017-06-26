/*
 *	This file is part of qpPresolver.
 *
 *	qpPresolver -- An implementation of presolving (= preprocessing) techniques
 *  for Quadratic Programming.
 *	Copyright (C) 2017 by Dominik Cebulla et al. All rights reserved.
 *
 *	qpPresolver is free software; you can redistribute it and/or
 *	modify it under the terms of the GNU Lesser General Public
 *	License as published by the Free Software Foundation; either
 *	version 2.1 of the License, or (at your option) any later version.
 *
 *	qpPresolver is distributed in the hope that it will be useful,
 *	but WITHOUT ANY WARRANTY; without even the implied warranty of
 *	MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
 *	See the GNU Lesser General Public License for more details.
 *
 *	You should have received a copy of the GNU Lesser General Public
 *	License along with qpPresolver; if not, write to the Free Software
 *	Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
 *
 */


/**
 *	\file src/ecrmatrix.c
 *	\author Dominik Cebulla
 *	\version 1.0 Beta
 *	\date 2017
 *
 *  Implementation of a sparse matrix class stored in (extended) Curtis-Reid format
 *  (basically based on singly linked lists).
 */


#include <qpPresolver/ecrmatrix.h>


/*  ======================================================================================
    Forward declaration of static (= private) functions
    ====================================================================================*/

/** \brief Sorts a matrix given in coordinate scheme.
 *
 *  First we sort w.r.t. \p jcn and then w.r.t. \p irn. If \p irn corresponds to the row
 *  subscripts of the matrix (and \p jcn to the column subscripts), then the matrix
 *  is sorted row-wise afterwards (by changing \p irn and \p jcn, the matrix will be
 *  sorted column-wise). We use Counting Sort (stable, linear complexity (fast!),
 *  but unfortunately not an in-place algorithm).
 */
static qpp_return_value_t sortCoordinateCS(const qpp_int_t nrow,
                                           const qpp_int_t ncol,
                                           const qpp_int_t nnz,
                                           qpp_int_t irn[],
                                           qpp_int_t jcn[],
                                           qpp_real_t x[]);


/** \brief Computes the maximum of two integer numbers. */
static inline qpp_int_t max(const qpp_int_t m,
                            const qpp_int_t n);


/*  ======================================================================================
    Implementation of interface functions.
    ====================================================================================*/

qpp_ecrmatrix_t *ecrAlloc(const qpp_int_t nrowmax,
                          const qpp_int_t ncolmax,
                          const qpp_int_t nnzmax)
{
    qpp_ecrmatrix_t *A;
    qpp_int_t memsize, i;
    qpp_int_t *iptr;

    /* Check input: */
    if ( (nnzmax < 0) || (nrowmax < 0) || (ncolmax < 0) ||
         ((nnzmax > 0) && ((ncolmax == 0) || (nrowmax == 0))) )
    {
        return NULL;
    }

    /* A m x n matrix may not contain more than m*n nonzero elements. */
    if (nrowmax > 0)
    {
        if ( ((nnzmax/nrowmax == ncolmax) && (nnzmax%nrowmax > 0)) ||
            (nnzmax/nrowmax > ncolmax) )
        {
            return NULL;
        }
    }

    A = (qpp_ecrmatrix_t*) malloc(sizeof(qpp_ecrmatrix_t));
    if (A == NULL)
    {
        return NULL;
    }

    memsize = (ncolmax + nrowmax + 4*nnzmax) * (qpp_int_t)(sizeof(qpp_int_t)) +
              nnzmax * (qpp_int_t)(sizeof(qpp_real_t));
    if (memsize == 0)
    {
        /* We allow to create an empty matrix! */
        A->irowst = NULL;
        A->jcolst = NULL;
        A->irn    = NULL;
        A->jcn    = NULL;
        A->linkrw = NULL;
        A->linkcl = NULL;
        A->x = NULL;
        A->nnz = 0;
        A->nrow = 0;
        A->ncol = 0;
        A->nnzmax = 0;
        A->nrowmax = 0;
        A->ncolmax = 0;
    }
    else
    {
        iptr = (qpp_int_t*) malloc((size_t)(memsize));

        if (iptr == NULL)
        {
            free(A);
            return NULL;
        }

        A->irowst = iptr; iptr += nrowmax;
        A->jcolst = iptr; iptr += ncolmax;
        A->irn    = iptr; iptr += nnzmax;
        A->jcn    = iptr; iptr += nnzmax;
        A->linkrw = iptr; iptr += nnzmax;
        A->linkcl = iptr; iptr += nnzmax;
        A->x = (void*) iptr;
        A->nnz = 0;
        A->nrow = 0;
        A->ncol = 0;
        A->nnzmax = nnzmax;
        A->nrowmax = nrowmax;
        A->ncolmax = ncolmax;
    }

    for (i = 0; i < nrowmax; ++i)
    {
        A->irowst[i] = -1;
    }

    for (i = 0; i < ncolmax; ++i)
    {
        A->jcolst[i] = -1;
    }

    return A;
}


qpp_return_value_t ecrRealloc(qpp_ecrmatrix_t *A,
                              const qpp_int_t nnzmax)
{
    qpp_int_t i, it1, it2, nnzmax_old, memsize;
    qpp_int_t *A_linkrw, *A_linkcl, *A_irowst, *A_jcolst, *iptr;

    if (A == NULL)
    {
        return QPP_NULL_ARGUMENT;
    }

    if (nnzmax < 0)
    {
        return QPP_INVALID_ARGUMENT;
    }

    if (A->nnzmax == nnzmax)
    {
        return QPP_OK;
    }

    if ( ((nnzmax > 0) && ((A->ncolmax == 0) || (A->nrowmax == 0))) )
    {
        return QPP_INVALID_ARGUMENT;
    }

    /* A m x n matrix may not contain more than m*n nonzero elements. */
    if (A->nrowmax > 0)
    {
        if ( ((nnzmax / A->nrowmax == A->ncolmax) && (nnzmax % A->nrowmax > 0)) ||
            (nnzmax/A->nrowmax > A->ncolmax) )
        {
            return QPP_INVALID_ARGUMENT;
        }
    }

    nnzmax_old = A->nnzmax;

    if (nnzmax == 0)    /* Create an empty <nrowmax> x <ncolmax> matrix */
    {
        A->nnz = 0;
        A->nrow = 0;
        A->ncol = 0;

        for (i = 0; i < A->nrowmax; ++i)
        {
            A->irowst[i] = -1;
        }
        for (i = 0; i < A->ncolmax; ++i)
        {
            A->jcolst[i] = -1;
        }
    }

    if ( (nnzmax < A->nnz) && (nnzmax > 0) )
    {
        /*  We have to adjust the links since the matrix will get "smaller" and only
            the first <nnzmax> nonzero elements will "survive". */
        A_linkrw = A->linkrw;
        A_linkcl = A->linkcl;
        A_irowst = A->irowst;
        A_jcolst = A->jcolst;

        /* Adjust row links */
        for (i = 0; i < A->nrow; ++i)
        {
            it1 = A_irowst[i];
            /* Get the (possibly new) "anchor" element for row i */
            while ( (it1 >= nnzmax) && (it1 >= 0) )
            {
                it1 = A_linkrw[it1];
            }
            A_irowst[i] = it1;   /* Can also be negative (which is ok)! */

            /*  Now we remove row links if one element will be deleted (this happens,
                when the index of one element is >= than nnzmax */
            while (it1 >= 0)
            {
                it2 = A_linkrw[it1];
                while ( (it2 >= nnzmax) && (it2 >= 0) )
                {
                    it2 = A_linkrw[it2];
                }
                A_linkrw[it1] = it2;
                it1 = it2;
            }
        }

        /* Adjust column links */
        for (i = 0; i < A->ncol; ++i)
        {
            it1 = A_jcolst[i];
            while ( (it1 >= nnzmax) && (it1 >= 0) )
            {
                it1 = A_linkcl[it1];
            }
            A_jcolst[i] = it1;

            while (it1 >= 0)
            {
                it2 = A_linkcl[it1];
                while ( (it2 >= nnzmax) && (it2 >= 0) )
                {
                    it2 = A_linkcl[it2];
                }
                A_linkcl[it1] = it2;
                it1 = it2;
            }
        }
        A->nnz = nnzmax;
    }

    /* If we shrink size, then data has to be moved before reallocating */
    if (A->nnzmax > nnzmax)
    {
        A->jcn    = memmove(&A->irn[nnzmax], A->jcn, (size_t)(nnzmax) * sizeof(qpp_int_t));
        A->linkrw = memmove(&A->jcn[nnzmax], A->linkrw, (size_t)(nnzmax) * sizeof(qpp_int_t));
        A->linkcl = memmove(&A->linkrw[nnzmax], A->linkcl, (size_t)(nnzmax) * sizeof(qpp_int_t));
        A->x      = memmove(&A->linkcl[nnzmax], A->x, (size_t)(nnzmax) * sizeof(qpp_real_t));
    }

    /* Reallocate memory */
    memsize = (A->nrowmax + A->ncolmax + 4*nnzmax) * sizeof(qpp_int_t) +
              nnzmax * sizeof(qpp_real_t);
    iptr = (qpp_int_t*) realloc(A->irowst, memsize);

    if ( (iptr == NULL) && (nnzmax < A->nnzmax) )
    {
        /*  Reallocating did not succeed although we shrank the size! The data
            of <A> was changed, so <A> should not be used now */
        return QPP_FATAL_ERROR;
    }

    if ( (iptr == NULL) && (nnzmax > A->nnzmax) )
    {
        return QPP_OUT_OF_MEMORY;
    }

    A->nnzmax = nnzmax;
    A->irowst = iptr; iptr += A->nrowmax;
    A->jcolst = iptr; iptr += A->ncolmax;
    A->irn    = iptr; iptr += nnzmax;
    A->jcn    = iptr; iptr += nnzmax;
    A->linkrw = iptr; iptr += nnzmax;
    A->linkcl = iptr; iptr += nnzmax;
    A->x = (void*) iptr;

    /* If we enlarged the size, then we have to move memory now */
    if (nnzmax_old < nnzmax)
    {
        memmove(A->x, &A->irn[4*nnzmax_old], nnzmax_old * sizeof(qpp_real_t));
        memmove(A->linkcl, &A->irn[3*nnzmax_old], nnzmax_old * sizeof(qpp_int_t));
        memmove(A->linkrw, &A->irn[2*nnzmax_old], nnzmax_old * sizeof(qpp_int_t));
        memmove(A->jcn, &A->irn[nnzmax_old], nnzmax_old * sizeof(qpp_int_t));
    }
    return QPP_OK;
}


qpp_ecrmatrix_t *ecrCopy(const qpp_ecrmatrix_t *A)
{
    qpp_ecrmatrix_t *C;
    qpp_int_t memsize;

    if (A == NULL)
    {
        return NULL;
    }

    memsize = (A->ncolmax + A->nrowmax + 4*A->nnzmax) * sizeof(qpp_int_t) +
              A->nnzmax * sizeof(qpp_real_t);

    C = ecrAlloc(A->nrowmax, A->ncolmax, A->nnzmax);
    if (C == NULL)
    {
        return NULL;
    }

    memcpy(C->irowst, A->irowst, memsize);
    C->ncol = A->ncol;
    C->nrow = A->nrow;
    C->nnz = A->nnz;

    return C;
}


void ecrFree(qpp_ecrmatrix_t *A)
{
    if (A != NULL)
    {
        free(A->irowst);
        free(A);
    }
}


qpp_return_value_t ecrCreateFromCoordinate(qpp_ecrmatrix_t *A,
                                           const qpp_int_t nrow,
                                           const qpp_int_t ncol,
                                           const qpp_int_t nnz,
                                           const qpp_int_t irn[],
                                           const qpp_int_t jcn[],
                                           const qpp_real_t x[],
                                           const qpp_matrix_sort_type_t sort_type,
                                           const qpp_bool_type_t only_lt)
{
    qpp_int_t i, it, is_sorted, row_index, col_index, status, nnz_lt;
    qpp_int_t *A_irn, *A_jcn, *A_irowst, *A_jcolst, *A_linkrw, *A_linkcl,
        *pred_row, *pred_col;
    qpp_real_t *A_x;

    if (A == NULL)
    {
        return QPP_NULL_ARGUMENT;
    }

    if ( (A->nrowmax < nrow) || (A->ncolmax < ncol) || (ncol < 0) ||
         (nrow < 0) || (nnz < 0) || ((only_lt != 0) && (ncol != nrow)) )
    {
        return QPP_INVALID_ARGUMENT;
    }

    if ( (nnz > 0) && ((irn == NULL) || (jcn == NULL) || (x == NULL)) )
    {
        return QPP_NULL_ARGUMENT;
    }

    /*  Check if data is consistent, i.e. no row index given by <irn> may be greater
        than <nrow> (same for columns). If only the lower triangular part of the matrix
        shall be stored, compute the number of nonzeros. */
    nnz_lt = 0;
    for (i = 0; i < nnz; ++i)
    {
        if (irn[i] >= nrow)
        {
            return QPP_INVALID_ARGUMENT;
        }

        if (jcn[i] >= ncol)
        {
            return QPP_INVALID_ARGUMENT;
        }

        if (irn[i] >= jcn[i])
        {
            ++nnz_lt;
        }
    }

    if (only_lt == QPP_BT_TRUE)
    {
        /* Enlarge matrix if necessary */
        if (A->nnzmax < nnz_lt)
        {
            status = ecrRealloc(A, nnz_lt);
            if (status != QPP_OK)
            {
                return status;
            }
        }

        /* If the number of nonzeros can be reduced by more than 10%: Shrink to fit */
        if ( ((qpp_real_t)(A->nnzmax-nnz_lt) / A->nnzmax) >= 0.1)
        {
            status = ecrRealloc(A, nnz_lt);
            if (status != QPP_OK)
            {
                return status;
            }
        }
    }
    else
    {
        if (A->nnzmax < nnz)
        {
            status = ecrRealloc(A, nnz);
            if (status != QPP_OK)
            {
                return status;
            }
        }
    }


    A_x = A->x;
    A_irn = A->irn;
    A_jcn = A->jcn;
    A_irowst = A->irowst;
    A_jcolst = A->jcolst;
    A_linkrw = A->linkrw;
    A_linkcl = A->linkcl;

    if (only_lt == QPP_BT_TRUE)
    {
        it = 0;
        for (i = 0; i < nnz; ++i)
        {
            if (irn[i] >= jcn[i])
            {
                A_x[it] = x[i];
                A_irn[it] = irn[i];
                A_jcn[it] = jcn[i];
                ++it;
            }
        }
    }
    else
    {
        for (i = 0; i < nnz; ++i)
        {
            A_x[i] = x[i];
            A_irn[i] = irn[i];
            A_jcn[i] = jcn[i];
        }
    }

    if (only_lt == QPP_BT_FALSE)
    {
        nnz_lt = nnz;
    }

    /*  Check if the matrix is already sorted (row- or column-wise, depends on parameter
        <sort_type>), because we then do not sort the matrix. */
    is_sorted = 1;
    if (sort_type == QPP_MST_ROW_WISE)
    {
        for (i = 1; i < nnz_lt; ++i)
        {
            if ( ((A_irn[i-1] == A_irn[i]) && (A_jcn[i-1] > A_jcn[i])) || (A_irn[i-1] > A_irn[i]) )
            {
                is_sorted = 0;
                break;
            }
        }
    }
    else
    {
        for (i = 1; i < nnz; ++i)
        {
            if ( ((A_jcn[i-1] == A_jcn[i]) && (A_irn[i-1] > A_irn[i])) || (A_jcn[i-1] > A_jcn[i]) )
            {
                is_sorted = 0;
                break;
            }
        }
    }

    if (!is_sorted)
    {
        if (sort_type == QPP_MST_ROW_WISE)
        {
            status = sortCoordinateCS(nrow, ncol, nnz_lt, A_irn, A_jcn, A_x);
        }
        else
        {
            status = sortCoordinateCS(nrow, ncol, nnz_lt, A_jcn, A_irn, A_x);
        }

        if (status != QPP_OK)
        {
            return status;
        }
    }

    /*  Now with have to create all links between the nonzero elements. */
    pred_row = (qpp_int_t*) malloc((nrow+ncol) * sizeof(qpp_int_t));
    if (pred_row == NULL)
    {
        return QPP_OUT_OF_MEMORY;
    }
    pred_col = &pred_row[nrow];

    for (i = 0; i < nnz_lt; ++i)
    {
        row_index = A_irn[i];
        col_index = A_jcn[i];

        if (A_irowst[row_index] < 0)
        {
            A_irowst[row_index] = i;
        }
        else
        {
            A_linkrw[pred_row[row_index]] = i;
        }
        A_linkrw[i] = -(row_index+1);
        pred_row[row_index] = i;

        if (A_jcolst[col_index] < 0)
        {
            A_jcolst[col_index] = i;
        }
        else
        {
            A_linkcl[pred_col[col_index]] = i;
        }
        A_linkcl[i] = -(col_index+1);
        pred_col[col_index] = i;
    }

    free(pred_row);

    A->ncol = ncol;
    A->nrow = nrow;
    A->nnz = nnz_lt;

    return QPP_OK;
}


void ecrAxpy(const char transp,
             const qpp_ecrmatrix_t *A,
             const qpp_real_t a,
             const qpp_real_t x[],
             qpp_real_t y[])
{
    qpp_int_t i, it, length;
    qpp_int_t *link, *start, *indices;
    qpp_real_t temp;
    qpp_real_t *A_x;

    if ( (A == NULL) || (x == NULL) || (y == NULL) )
    {
        return;
    }
    A_x = A->x;

    if ( (transp == 'T') || (transp == 't') )
    {
        link = A->linkcl;
        start = A->jcolst;
        indices = A->irn;
        length = A->ncol;
    }
    else
    {
        link = A->linkrw;
        start = A->irowst;
        indices = A->jcn;
        length = A->nrow;
    }

    for (i = 0; i < length; ++i)
    {
        temp = 0.0;
        it = start[i];
        while (it >= 0)
        {
            temp += A_x[it] * x[indices[it]];
            it = link[it];
        }
        y[i] += a * temp;
    }
}


qpp_bool_type_t ecrIsTriangular(const qpp_ecrmatrix_t *A,
                                const qpp_int_t type)
{
    qpp_int_t i;
    qpp_int_t *A_irn, *A_jcn;

    if (A == NULL)
    {
        return QPP_BT_FALSE;
    }

    A_irn = A->irn;
    A_jcn = A->jcn;

    if (type >= 0)
    {
        for (i = 0; i < A->nnz; ++i)
        {
            if (A_irn[i] > A_jcn[i])
            {
                return QPP_BT_FALSE;
            }
        }
    }
    else
    {
        for (i = 0; i < A->nnz; ++i)
        {
            if (A_irn[i] < A_jcn[i])
            {
                return QPP_BT_FALSE;
            }
        }
    }
    return QPP_BT_TRUE;
}


qpp_return_value_t ecrScaleCR(qpp_ecrmatrix_t *A,
                              qpp_real_t R[],
                              qpp_real_t C[],
                              const qpp_int_t max_iter)
{
    qpp_int_t m, n, nnz, i, k;
    qpp_int_t *A_irn, *A_jcn;
    qpp_real_t value, s[2], q[2], e[3];
    qpp_real_t *A_x, *rptr, *M, *N, *gamma1, *gamma2, *rho1, *rho2, *r;

    if ( (A == NULL) || (R == NULL) || (C == NULL) )
    {
        return QPP_NULL_ARGUMENT;
    }

    if ( (A->nrow <= 0) || (A->ncol <= 0) || (A->nnz <= 0) )
    {
        return QPP_INVALID_ARGUMENT;
    }

    m = A->nrowmax;
    n = A->ncolmax;
    nnz = A->nnz;

    A_irn = A->irn;
    A_jcn = A->jcn;
    A_x = A->x;
    M = R;
    N = C;

    gamma1 = (qpp_real_t*) malloc(n * sizeof(qpp_real_t));
    gamma2 = (qpp_real_t*) malloc(n * sizeof(qpp_real_t));
    rho1 = (qpp_real_t*) malloc(m * sizeof(qpp_real_t));
    rho2 = (qpp_real_t*) malloc(m * sizeof(qpp_real_t));
    r = (qpp_real_t*) malloc((m+n) * sizeof(qpp_real_t));

    if ( (r == NULL) || (rho2 == NULL) || (rho1 == NULL) ||
         (gamma2 == NULL) || (gamma1 == NULL) )
    {
        free(gamma1);
        free(gamma2);
        free(rho1);
        free(rho2);
        free(r);
        return QPP_OUT_OF_MEMORY;
    }

    for (i = 0; i < n; ++i)
    {
        gamma1[i] = 0.0;
        gamma2[i] = 0.0;
        N[i] = 0.0;
        r[m+i] = 0.0;
    }

    for (i = 0; i < m; ++i)
    {
        M[i] = 0.0;
        r[i] = 0.0;
        rho1[i] = 0.0;
    }

    for (i = 0; i < nnz; ++i)
    {
        M[A_irn[i]] += 1.0;
        N[A_jcn[i]] += 1.0;
        value = log2(fabs(A_x[i]));
        rho1[A_irn[i]] += value;
        r[m+A_jcn[i]] += value;
    }

    for (i = 0; i < m; ++i)
    {
        if (fabs(M[i]) > QPP_ZERO)
        {
            M[i] = 1.0 / M[i];
        }
        rho1[i] *= M[i];
        rho2[i] = rho1[i];
    }

    for (i = 0; i < n; ++i)
    {
        if (fabs(N[i]) > QPP_ZERO)
        {
            N[i] = 1.0 / N[i];
        }
    }

    for (i = 0; i < nnz; ++i)
    {
        r[m+A_jcn[i]] -= rho1[A_irn[i]];
    }

    q[0] = 1.0; q[1] = 0.0;
    s[0] = 0.0; s[1] = 0.0;
    e[0] = 0.0; e[1] = 0.0; e[2] = 0.0;

    for (i = 0; i < n; ++i)
    {
        s[0] += r[m+i] * r[m+i] * N[i];
    }

    k = 0;
    while (k < max_iter)
    {
        value = 0.0;
        if (k%2 == 0)
        {
            for (i = 0; i < m; ++i)
            {
                r[i] *= e[1];
            }
            for (i = 0; i < nnz; ++i)
            {
                r[A_irn[i]] += N[A_jcn[i]] * r[m+A_jcn[i]];
            }
            for (i = 0; i < m; ++i)
            {
                r[i] /= -q[0];
                value += r[i] * r[i] * M[i];
            }
            s[1] = value;
        }
        else
        {
            for (i = 0; i < n; ++i)
            {
                r[m+i] *= e[1];
            }
            for (i = 0; i < nnz; ++i)
            {
                r[m+A_jcn[i]] += M[A_irn[i]] * r[A_irn[i]];
            }
            for (i = 0; i < n; ++i)
            {
                r[m+i] /= -q[0];
                value += r[m+i] * r[m+i] * N[i];
            }
            s[1] = value;
        }

        if (s[1] <= 0.01 * nnz)
        {
            if (k%2 == 0)
            {
                for (i = 0; i < n; ++i)
                {
                    gamma2[i] = gamma2[i] + (N[i]*r[m+i] + e[0]*e[1]*(gamma2[i]-
                                gamma1[i]))/q[0];
                }
            }
            else
            {
                for (i = 0; i < m; ++i)
                {
                    rho2[i] = rho2[i] + (M[i]*r[i] + e[0]*e[1]*(rho2[i]-rho1[i])) / q[0];
                }
            }
            break;
        }

        e[2] = q[0] * s[1] / s[0];
        q[1] = 1.0 - e[2];

        if (k%2 == 0)
        {
            for (i = 0; i < n; ++i)
            {
                gamma1[i] = gamma2[i] + ( (N[i]*r[m+i] + e[0]*e[1]*(gamma2[i]-
                            gamma1[i])) / (q[0]*q[1]) );
            }
            rptr = gamma1;
            gamma1 = gamma2;
            gamma2 = rptr;
        }
        else
        {
            for (i = 0; i < m; ++i)
            {
                rho1[i] = rho2[i] + ( (M[i]*r[i] + e[0]*e[1]*(rho2[i]-rho1[i])) /
                          (q[0]*q[1]) );
            }
            rptr = rho1;
            rho1 = rho2;
            rho2 = rptr;
        }

        s[0] = s[1];
        e[0] = e[1];
        e[1] = e[2];
        q[0] = q[1];
        ++k;
    }

    /* We scale matrix A only if the procedure was successful! */
    if (k < max_iter)
    {
        /* Scaling based on power 2, so there will be no round-off error */
        for (i = 0; i < m; ++i)
        {
            R[i] = exp2(-floor(rho2[i]+0.5));
        }

        for (i = 0; i < n; ++i)
        {
            C[i] = exp2(-floor(gamma2[i]+0.5));
        }

        for (i = 0; i < nnz; ++i)
        {
            A_x[i] *= R[A_irn[i]] * C[A_jcn[i]];
        }
    }
    else
    {
        for (i = 0; i < m; ++i)
        {
            R[i] = 1.0;
        }

        for (i = 0; i < n; ++i)
        {
            C[i] = 1.0;
        }
    }

    free(r);
    free(gamma1);
    free(gamma2);
    free(rho1);
    free(rho2);

    if (k >= max_iter)
    {
        return QPP_BAD_SCALING;
    }

    return QPP_OK;
}


/*  ======================================================================================
    Static procedures
    ====================================================================================*/

static qpp_return_value_t sortCoordinateCS(const qpp_int_t nrow,
                                           const qpp_int_t ncol,
                                           const qpp_int_t nnz,
                                           qpp_int_t irn[],
                                           qpp_int_t jcn[],
                                           qpp_real_t x[])
{
    qpp_int_t i, k, memsize;
    qpp_int_t *rc, *cc, *t, *iptr;
    qpp_real_t *xc;

    if (nrow > ncol)
    {
        k = nrow;
    }
    else
    {
        k = ncol;
    }

    /* Initialize temporary storage (Counting Sort is not an in-place algorithm) */
    memsize = (2*nnz + k + 1) * sizeof(qpp_int_t) + nnz * sizeof(qpp_real_t);
    iptr = (qpp_int_t*) malloc(memsize);

    if (iptr == NULL)
    {
        return QPP_OUT_OF_MEMORY;
    }

    rc = iptr; iptr += nnz;
    cc = iptr; iptr += nnz;
    t = iptr; iptr += k+1;
    xc = (void*) iptr;

    for (i = 0; i <= ncol; ++i)
    {
        t[i] = 0;
    }

    for (i = 0; i < nnz; ++i)
    {
        ++t[jcn[i]];
    }

    for (i = 1; i <= ncol; ++i)
    {
        t[i] += t[i-1];
    }

    for (i = nnz-1; i >= 0; --i)
    {
        cc[t[jcn[i]]-1] = jcn[i];
        rc[t[jcn[i]]-1] = irn[i];
        xc[t[jcn[i]]-1] = x[i];
        --t[jcn[i]];
    }

    for (i = 0; i <= nrow; ++i)
    {
        t[i] = 0;
    }

    for (i = 0; i < nnz; ++i)
    {
        ++t[rc[i]];
    }

    for (i = 1; i <= nrow; ++i)
    {
        t[i] += t[i-1];
    }

    for (i = nnz-1; i >= 0; --i)
    {
        jcn[t[rc[i]]-1] = cc[i];
        irn[t[rc[i]]-1] = rc[i];
        x[t[rc[i]]-1] = xc[i];
        --t[rc[i]];
    }

    free(rc);

    return QPP_OK;
}


static inline qpp_int_t max(const qpp_int_t m,
                            const qpp_int_t n)
{
    return (m >= n ? m : n);
}
