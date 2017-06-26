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
 *	\file src/utility.c
 *	\author Dominik Cebulla
 *	\version 1.0 Beta
 *	\date 2017
 *
 *  Implementation of useful utility functions and some I/O functions for matrices
 *  stored in MatrixMarket format, see http://math.nist.gov/MatrixMarket/mmio-c.html.
 *
 */


 #include <qpPresolver/utility.h>


inline qpp_return_value_t throwError(const qpp_return_value_t error,
                                     const char* filename,
                                     const int line_number)
{
    fprintf(stderr, "\nErrorcode %d occurred in file %s at line %d.\n", error,
            filename, line_number);

    return error;
}


inline void qppSetArrayi(qpp_int_t x[],
                         qpp_int_t n,
                         const qpp_int_t value)
{
    while (n--)
    {
        *(x++) = value;
    }
}


inline void qppSetArrayr(qpp_real_t x[],
                         qpp_int_t n,
                         const qpp_real_t value)
{
    while (n--)
    {
        *(x++) = value;
    }
}


inline void qppSwapi(qpp_int_t* a,
                     qpp_int_t* b)
{
    qpp_int_t t;
    t = *b; *b = *a; *a = t;
}


inline void qppSwapr(qpp_real_t* a,
                     qpp_real_t* b)
{
    qpp_real_t t;
    t = *b; *b = *a; *a = t;
}


inline void qppPrintArrayi(const qpp_int_t x[],
                           const qpp_int_t n)
{
    qpp_int_t i;
    for (i = 0; i < n; ++i)
    {
        fprintf(stdout, "\t%" QPP_PRID "\n", x[i]);
    }
    fprintf(stdout, "\n");
}


inline void qppPrintArrayr(const qpp_real_t x[],
                           const qpp_int_t n)
{
    qpp_int_t i;
    for (i = 0; i < n; ++i)
    {
        fprintf(stdout, "\t%.4e\n", x[i]);
    }
    fprintf(stdout, "\n");
}


inline qpp_real_t qppMaxr(const qpp_real_t a,
                          const qpp_real_t b)
{
    return (a >= b ? a : b);
}


inline qpp_real_t qppMinr(const qpp_real_t a,
                          const qpp_real_t b)
{
    return (a <= b ? a : b);
}


inline qpp_int_t qppMaxi(const qpp_int_t a,
                         const qpp_int_t b)
{
    return (a >= b ? a : b);
}


inline qpp_int_t qppMini(const qpp_int_t a,
                         const qpp_int_t b)
{
    return (a <= b ? a : b);
}


inline qpp_bool_type_t qppIsEqual(const qpp_real_t x,
                                  const qpp_real_t y,
                                  const qpp_real_t tol)
{
    if ( isfinite(x) && isfinite(y) &&
        (fabs(x-y) <= tol * qppMaxr(1.0, qppMaxr(fabs(x), fabs(y)))) )
    {
        return QPP_BT_TRUE;
    }
    return QPP_BT_FALSE;
}


inline qpp_bool_type_t qppIsGreater(const qpp_real_t x,
                                    const qpp_real_t y,
                                    const qpp_real_t tol)
{
    if ( (x == QPP_INF) && (y != QPP_INF) )
    {
        return QPP_BT_TRUE;
    }
    if (isfinite(x) && isfinite(y) && (x-y > tol*qppMaxr(1.0, qppMaxr(fabs(x), fabs(y)))))
    {
        return QPP_BT_TRUE;
    }
    return QPP_BT_FALSE;
}


inline qpp_real_t qppSparseDot(const qpp_int_t indices[],
                               const qpp_real_t values[],
                               const qpp_int_t length,
                               const qpp_real_t x[])
{
    qpp_int_t i;
    qpp_real_t result;

    result = 0.0;

    for (i = 0; i < length; ++i)
    {
        result += values[i] * x[indices[i]];
    }
    return result;
}


qpp_return_value_t qppSortNegPos(qpp_int_t indices[],
                                 qpp_real_t values[],
                                 const qpp_int_t length,
                                 qpp_int_t *n_neg_entries,
                                 qpp_int_t *n_pos_entries)
{
    qpp_int_t i, k;

    if ( (indices == NULL) || (values == NULL) || (n_neg_entries == NULL) ||
         (n_pos_entries == NULL) )
    {
        return QPP_NULL_ARGUMENT;
    }

    i = 0;
    k = length - 1;

    while (1)
    {
        while ( (values[i] < 0.0) && (i < length) )
        {
            ++i;
        }

        while ( (values[k] > 0.0) && (k >= 0) )
        {
            --k;
        }

        if (k <= i)
        {
            break;
        }

        qppSwapi(&indices[i], &indices[k]);
        qppSwapr(&values[i], &values[k]);
        ++i;
        --k;
    }
    *n_neg_entries = i;
    *n_pos_entries = length - i;

    return QPP_OK;
}


qpp_return_value_t qppMMReadDenseVector(const char *filename,
                                        qpp_real_t *x[],
                                        qpp_int_t *n)
{
    FILE *file;
    MM_typecode mmtc;
    qpp_int_t i, num_rows, num_cols, nnz;
    qpp_real_t *mat_x;
    int num_items_read;
    char line[QPP_MAX_STRING_LENGTH];
    double xd;

    if ( (n == NULL) || (x == NULL) )
    {
        return QPP_NULL_ARGUMENT;
    }

    *x = NULL;
    *n = 0;

    if (strcmp(filename, "stdin") == 0)
    {
        file = stdin;
    }
    else
    {
        file = fopen(filename, "r");
        if (file == NULL)
        {
            return QPP_UNABLE_TO_OPEN_FILE;
        }
    }

    if (mm_read_banner(file, &mmtc) != 0)
    {
        if (file != stdin)
        {
            fclose(file);
        }
        return QPP_UNABLE_TO_READ_FILE;
    }

    if ( !mm_is_matrix(mmtc) || !mm_is_array(mmtc) || !mm_is_real(mmtc) )
    {
        if (file != stdin)
        {
            fclose(file);
        }
        return QPP_INVALID_ARGUMENT;
    }

    if ( !(mm_is_general(mmtc) || mm_is_symmetric(mmtc)) )
    {
        if (file != stdin)
        {
            fclose(file);
        }
        return QPP_INVALID_ARGUMENT;
    }

    /*  We ensure that we can read the size of the matrix (i.e. m and n) as 64bit
        integers. Therefore we have to search manually for the desired information. */

    /*  First ignore all comment lines (starting with '%') */
    do
    {
        if (fgets(line, QPP_MAX_STRING_LENGTH, file) == NULL)
        {
            if (file != stdin)
            {
                fclose(file);
            }
            return QPP_UNABLE_TO_READ_FILE;
        }
    } while (line[0] == '%');

    /*  Now line is either blank or contains m and n. */
    num_items_read = sscanf(line, "%" QPP_SCND " %" QPP_SCND "", &num_rows, &num_cols);

    while (num_items_read != 2)
    {
        num_items_read = fscanf(file, "%" QPP_SCND " %" QPP_SCND "", &num_rows, &num_cols);
        if (num_items_read == EOF)
        {
            if (file != stdin)
            {
                fclose(file);
            }
            return QPP_INVALID_ARGUMENT;
        }
    }

    /* Check if the data stored in the file is a vector */
    if ( (num_rows <= 0) || (num_cols <= 0) || ((num_rows > 1) && (num_cols > 1)) )
    {
        if (file != stdin)
        {
            fclose(file);
        }
        return QPP_INVALID_ARGUMENT;
    }

    if (num_rows == 1)
    {
        nnz = num_cols;
    }
    else
    {
        nnz = num_rows;
    }

    mat_x = (qpp_real_t*) malloc(nnz * sizeof(qpp_real_t));

    if (mat_x == NULL)
    {
        if (file != stdin)
        {
            fclose(file);
        }
        return QPP_OUT_OF_MEMORY;
    }

    for (i = 0; i < nnz; ++i)
    {
        if (fgets(line, QPP_MAX_STRING_LENGTH, file) == NULL)
        {
            free(mat_x);
            if (file != stdin)
            {
                fclose(file);
            }
            return QPP_UNABLE_TO_READ_FILE;
        }

        if ( strstr(line, "-Inf") || strstr(line, "-INF") || strstr(line, "-inf") )
        {
            mat_x[i] = -QPP_INF;
        }
        else if ( strstr(line, "Inf") || strstr(line, "INF") || strstr(line, "inf") )
        {
            mat_x[i] = QPP_INF;
        }
        else
        {
            if (sscanf(line, "%lg\n", &xd) != 1)
            {
                free(mat_x);
                if (file != stdin)
                {
                    fclose(file);
                }
                return QPP_INVALID_ARGUMENT;
            }
            mat_x[i] = (qpp_real_t) xd;
        }
    }

    *x = mat_x;
    *n = nnz;

    if (file != stdin)
    {
        fclose(file);
    }

    return QPP_OK;
}


qpp_return_value_t qppMMWriteDenseVector(const char *filename,
                                         const qpp_real_t x[],
                                         const qpp_int_t n)
{
    qpp_int_t i;
    MM_typecode mmtc;
    FILE *file;

    if (strcmp(filename, "stdout") == 0)
    {
        file = stdout;
    }
    else
    {
        file = fopen(filename, "w");
        if (file == NULL)
        {
            return QPP_UNABLE_TO_OPEN_FILE;
        }
    }

    mm_clear_typecode(&mmtc);
    mm_set_matrix(&mmtc);
    mm_set_array(&mmtc);
    mm_set_real(&mmtc);
    mm_set_general(&mmtc);

    if (mm_write_banner(file, mmtc) != 0)
    {
        if (file != stdout)
        {
            fclose(file);
        }
        return QPP_UNABLE_TO_WRITE_FILE;
    }

    /*  Write (dense) matrix size manually */
    if (fprintf(file, "%" QPP_PRID " %d\n", n, 1) <= 0)
    {
        if (file != stdout)
        {
            fclose(file);
        }
        return QPP_UNABLE_TO_WRITE_FILE;
    }

    for (i = 0; i < n; ++i)
    {
        if (fprintf(file, "%lg\n", (double) x[i]) <= 0)
        {
            if (file != stdout)
            {
                fclose(file);
            }
            return QPP_UNABLE_TO_WRITE_FILE;
        }
    }

    if (file != stdout)
    {
        fclose(file);
    }

    return QPP_OK;
}


qpp_return_value_t qppMMReadSparseCoordinate(const char *filename,
                                             qpp_int_t *m,
                                             qpp_int_t *n,
                                             qpp_int_t *nnz,
                                             qpp_int_t *irn[],
                                             qpp_int_t *jcn[],
                                             qpp_real_t *x[],
                                             const qpp_int_t store_ut)
{
    FILE *file;
    MM_typecode mmtc;
    qpp_int_t i, it, num_items_read, sym_type, sym_nnz, mat_m, mat_n, mat_nnz;
    qpp_int_t *mat_irn, *mat_jcn, *sym_irn, *sym_jcn;
    qpp_real_t *mat_x, *sym_x;
    double xd;
    char line[QPP_MAX_STRING_LENGTH];

    if ( (m == NULL) || (n == NULL) || (nnz == NULL) || (irn == NULL) ||
         (jcn == NULL) || (x == NULL) )
    {
        return QPP_NULL_ARGUMENT;
    }

    *n = 0;
    *m = 0;
    *nnz = 0;
    *irn = NULL;
    *jcn = NULL;
    *x = NULL;

    if (strcmp(filename, "stdin") == 0)
    {
        file = stdin;
    }
    else
    {
        file = fopen(filename, "r");
        if (file == NULL)
        {
            return QPP_UNABLE_TO_OPEN_FILE;
        }
    }

    if (mm_read_banner(file, &mmtc))
    {
        if (file != stdin)
        {
            fclose(file);
        }
        return QPP_UNABLE_TO_READ_FILE;
    }

    if ( !mm_is_matrix(mmtc) || !mm_is_coordinate(mmtc) || !mm_is_real(mmtc) )
    {
        if (file != stdin)
        {
            fclose(file);
        }
        return QPP_INVALID_ARGUMENT;
    }

    if (!( mm_is_general(mmtc) || mm_is_symmetric(mmtc) || mm_is_skew(mmtc) ))
    {
        if (file != stdin)
        {
            fclose(file);
        }
        return QPP_INVALID_ARGUMENT;
    }

    sym_type = 0;
    if (mm_is_symmetric(mmtc))
    {
        sym_type = 1;
    }
    else if (mm_is_skew(mmtc))
    {
        sym_type = -1;
    }

    /*  We ensure that we can read the size of the matrix (i.e. m, n, nnz) as 64bit
        integers (if we use 64Bit integers). Therefore we have to search manually
        for the desired information. */

    /*  First ignore all comment lines (starting with '%') */
    do
    {
        if (fgets(line, QPP_MAX_STRING_LENGTH, file) == NULL)
        {
            if (file != stdin)
            {
                fclose(file);
            }
            return QPP_UNABLE_TO_READ_FILE;
        }
    } while (line[0] == '%');

    /*  Now line is either blank or contains m, n and nnz */
    num_items_read = sscanf(line, "%" QPP_SCND " %" QPP_SCND " %" QPP_SCND "", &mat_m,
                            &mat_n, &mat_nnz);

    while (num_items_read != 3)
    {
        num_items_read = fscanf(file, "%" QPP_SCND " %" QPP_SCND " %" QPP_SCND "",
                                &mat_m, &mat_n, &mat_nnz);
        if (num_items_read == EOF)
        {
            if (file != stdin)
            {
                fclose(file);
            }
            return QPP_INVALID_ARGUMENT;
        }
    }

    /*  Now in mat_m, mat_n, and mat_nnz should be the correct size of the sparse matrix */
    mat_irn = (qpp_int_t*)  malloc(mat_nnz * sizeof(qpp_int_t));
    mat_jcn = (qpp_int_t*)  malloc(mat_nnz * sizeof(qpp_int_t));
    mat_x   = (qpp_real_t*) malloc(mat_nnz * sizeof(qpp_real_t));

    if ( (mat_irn == NULL) || (mat_jcn == NULL) || (mat_x == NULL) )
    {
        free(mat_irn);
        free(mat_jcn);
        free(mat_x);

        if (file != stdin)
        {
            fclose(file);
        }
        return QPP_OUT_OF_MEMORY;
    }

    sym_nnz = mat_nnz;   /* Number of nonzeros if we have to write the upper triangular part */

    for (i = 0; i < mat_nnz; ++i)
    {
        if (fscanf(file, "%" QPP_SCND " %" QPP_SCND " %lg\n", &mat_irn[i],
                   &mat_jcn[i], &xd) != 3)
        {
            free(mat_irn);
            free(mat_jcn);
            free(mat_x);

            if (file != stdin)
            {
                fclose(file);
            }
            return QPP_INVALID_ARGUMENT;
        }

        mat_x[i] = (qpp_real_t) xd;

        /* Convert 1-based to 0-based */
        --mat_jcn[i];
        --mat_irn[i];

        if (mat_irn[i] != mat_jcn[i])
        {
            ++sym_nnz;
        }
    }

    /*  Maybe the matrix is symmetric and we must also write the upper
        triangular part of it. */
    if ( (store_ut != 0) && (sym_type != 0) )
    {
        sym_irn = (qpp_int_t*)  realloc(mat_irn, sym_nnz * sizeof(qpp_int_t));
        sym_jcn = (qpp_int_t*)  realloc(mat_jcn, sym_nnz * sizeof(qpp_int_t));
        sym_x   = (qpp_real_t*) realloc(mat_x, sym_nnz * sizeof(qpp_real_t));

        if ( (sym_irn == NULL) || (sym_jcn == NULL) || (sym_x == NULL) )
        {
            if (sym_irn == NULL)
            {
                free(mat_irn);
            }
            else
            {
                free(sym_irn);
            }

            if (sym_jcn == NULL)
            {
                free(mat_jcn);
            }
            else
            {
                free(sym_jcn);
            }

            if (sym_x == NULL)
            {
                free(mat_x);
            }
            else
            {
                free(sym_x);
            }

            if (file != stdin)
            {
                fclose(file);
            }
            return QPP_OUT_OF_MEMORY;
        }

        /* Now we write the upper triangular part: */
        it = 0;
        for (i = 0; i < mat_nnz; ++i)
        {
            if (sym_irn[i] != sym_jcn[i])
            {
                sym_irn[mat_nnz+it] = sym_jcn[i];
                sym_jcn[mat_nnz+it] = sym_irn[i];
                sym_x[mat_nnz+it] = sym_type * sym_x[i];
                ++it;
            }
        }
        mat_nnz = sym_nnz;
        mat_irn = sym_irn;
        mat_jcn = sym_jcn;
        mat_x = sym_x;
    }

    *irn = mat_irn;
    *jcn = mat_jcn;
    *x = mat_x;
    *nnz = mat_nnz;
    *m = mat_m;
    *n = mat_n;

    if (file != stdin)
    {
        fclose(file);
    }

    return QPP_OK;
}


qpp_return_value_t qppMMWriteSparseCoordinate(const char *filename,
                                              const qpp_int_t m,
                                              const qpp_int_t n,
                                              const qpp_int_t nnz,
                                              const qpp_int_t irn[],
                                              const qpp_int_t jcn[],
                                              const qpp_real_t x[],
                                              const qpp_int_t is_sym)
{
    FILE *file;
    MM_typecode mmtc;
    qpp_int_t i, num_upper_elements;

    if ( (m < 0) || (n < 0) || (nnz < 0) )
    {
        return QPP_INVALID_ARGUMENT;
    }

    if ( (irn == NULL) || (jcn == NULL) || (x == NULL) )
    {
        return QPP_NULL_ARGUMENT;
    }

    if ( (is_sym != 0) && (m != n) )
    {
        /* Matrix cannot be treated as (skew-) symmetric if it is not quadratic. */
        return QPP_INVALID_ARGUMENT;
    }

    if (strcmp(filename, "stdout") == 0)
    {
        file = stdout;
    }
    else
    {
        file = fopen(filename, "w");
        if (file == NULL)
        {
            return QPP_UNABLE_TO_OPEN_FILE;
        }
    }

    mm_clear_typecode(&mmtc);
    mm_set_matrix(&mmtc);
    mm_set_coordinate(&mmtc);
    mm_set_real(&mmtc);

    if (is_sym == 0)
    {
        mm_set_general(&mmtc);
    }
    else if (is_sym < 0)
    {
        mm_set_skew(&mmtc);
    }
    else
    {
        mm_set_symmetric(&mmtc);
    }

    if (mm_write_banner(file, mmtc) != 0)
    {
        if (file != stdout)
        {
            fclose(file);
        }
        return QPP_UNABLE_TO_WRITE_FILE;
    }

    num_upper_elements = 0;     /* Number of elements in upper triangular part */
    if (is_sym)
    {
        for (i = 0; i < nnz; ++i)
        {
            if (irn[i] < jcn[i])
            {
                ++num_upper_elements;
            }
        }
    }

    /*  Manually write size information */
    if (fprintf(file, "%" QPP_PRID " %" QPP_PRID " %" QPP_PRID "\n",
                m, n, nnz-num_upper_elements) <= 0)
    {
        if (file != stdout)
        {
            fclose(file);
        }
        return QPP_UNABLE_TO_WRITE_FILE;
    }

    for (i = 0; i < nnz; ++i)
    {
        if ( !is_sym || ( is_sym && (irn[i] >= jcn[i]) ) )
        {
            /*  Convert 0-based to 1-based! */
            if (fprintf(file, "%" QPP_PRID " %" QPP_PRID " %lg\n", irn[i]+1, jcn[i]+1,
                        (double) x[i]) <= 0)
            {
                if (file != stdout)
                {
                    fclose(file);
                }
                return QPP_UNABLE_TO_WRITE_FILE;
            }
        }
    }

    if (file != stdout)
    {
        fclose(file);
    }

    return QPP_OK;
}
