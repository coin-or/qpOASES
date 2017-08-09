#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include <mex.h>
#include <qpPresolver.h>

static qpp_int_t checkID(const mxArray* const ID);
static qpp_return_value_t setupOptions(qpp_options_t *opt,
                                       const mxArray *mat_opt);
static qpp_return_value_t hasOptionsValueD(const mxArray *opt,
                                           const char* const opt_str,
                                           double **opt_val);
static qpp_return_value_t hasOptionsValueC(const mxArray *opt,
                                           const char* const opt_str,
                                           char **opt_val);


#define QPP_MAX_NUM_PROBLEMS 100
static qpp_data_t *qppEntity[QPP_MAX_NUM_PROBLEMS] = {NULL};


/* mxIsScalar is only defined for Matlab R2015a upwards */
#define MX_IS_SCALAR(X) ( (mxGetM((X)) == 1) && (mxGetN((X)) == 1) )
#define MX_IS_DOUBLE(X) ( mxIsDouble((X)) && !mxIsComplex((X)) )


static void exitFcn()
{
    qpp_int_t i;
    for (i = 0; i < QPP_MAX_NUM_PROBLEMS; ++i) 
    {	
        qppFree(&qppEntity[i]);
        qppEntity[i] = NULL;
    }
}


void mexFunction(int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[])
{
    qpp_return_value_t exitflag;
    qpp_int_t task, m, n, id, i, j, Annz, Hnnz, pm, pn;
    qpp_int_t *Airn, *Ajcn, *Acp, *Hirn, *Hjcn, *Hcp, *wi, *wj;
    qpp_data_t *entity;
    double f;
    double *Ax, *Hx, *xl, *xu, *al, *au, *g, *x, *y, *z;
    mxClassID intclass;
    
    mexAtExit(exitFcn);

    if (sizeof(mwIndex) != sizeof(qpp_int_t))
    {
        mexErrMsgTxt("qpPresolver: Integer type mismatch (32bit <-> 64bit)!");
    }

    if (sizeof(qpp_int_t) == 4)
    {
        intclass = mxINT32_CLASS;
    }
    else if (sizeof(qpp_int_t) == 8)
    {
        intclass = mxINT64_CLASS;
    }
    else
    {
        mexErrMsgTxt("qpPresolver: Integer type mismatch (only 32 and 64 Bit supported!");
    }

    if ( (nrhs == 0) || !MX_IS_DOUBLE(prhs[0]) || !MX_IS_SCALAR(prhs[0]) )
    {
        mexErrMsgTxt("qpPresolver: Task missing / Task must be a scalar!");
    }

    /*  To make sure that casting does not lead to a wrong number, we add 0.5 and then
        we truncate. */
    task = (qpp_int_t) (mxGetScalar(prhs[0]) + 0.5);

    switch (task)
    {
        case 0:	 goto QPPRESOLVER_NEW;
        case 1:  goto QPPRESOLVER_INIT;
        case 10: goto QPPRESOLVER_PRESOLVE;
        case 11: goto QPPRESOLVER_POSTSOLVE;
        case 15: goto QPPRESOLVER_GETPRESOLVEDQP;
        case 16: goto QPPRESOLVER_GETSOLUTION;
        case 20: goto QPPRESOLVER_FREE;
        case 30: goto QPPRESOLVER_SETOPTIONS;
        default: mexErrMsgTxt("qpPresolver: Unknown task!");
    }

    return;

QPPRESOLVER_NEW:
    /*	==============================================================================
    prhs[1] = number of constraints m (scalar)
    prhs[2] = number of variables n (scalar)
    prhs[3] = (maximum) number of nonzero elements in A (scalar)
    prhs[4] = (maximum) number of nonzero elements in (lower triangular part of) H (scalar)
    prhs[5] = initial size of the presolve stack (scalar)
    plhs[0] = id of new qpPresolver entity (scalar)
    plhs[1] = exitflag (scalar)
    ==================================================================================*/

    /* First check if we can create a new entity */
    id = 0;
    while ( (qppEntity[id] != NULL) && (id < QPP_MAX_NUM_PROBLEMS) ) 
    {
        ++id;
    }
    if (id >= QPP_MAX_NUM_PROBLEMS)
    {
        mexErrMsgTxt("qpPresolver (NEW): Maximum number of entities reached!");
    }

    /* CHECK INPUT */
    if (nrhs != 6)
    {
        mexErrMsgTxt("qpPresolver (NEW): Invalid number of input arguments!");
    }
    
    if ( (nlhs < 1) || (nlhs > 2) )
    {
        mexErrMsgTxt("qpPresolver (NEW): Invalid number of output arguments!");
    }

    for (i = 1; i <= 5; ++i)
    {
        if (!MX_IS_DOUBLE(prhs[i]))
        {
	        mexErrMsgTxt("qpPresolver (NEW): Invalid (complex or non-double) input!");
        }
    }

    if (!MX_IS_SCALAR(prhs[1]))
    {
        mexErrMsgTxt("qpPresolver (NEW): Number of constraints must be a non-negative scalar!");
    }

    if (!MX_IS_SCALAR(prhs[2]))
    {
        mexErrMsgTxt("qpPresolver (NEW): Number of variables must be a positive scalar!");
    }

    if (!MX_IS_SCALAR(prhs[3]))
    {
        mexErrMsgTxt("qpPresolver (NEW): (Maximum) number of nonzeros in A \
            must be a non-negative scalar!");
    }

    if (!MX_IS_SCALAR(prhs[4]))
    {
        mexErrMsgTxt("qpPresolver (NEW): (Maximum) number of nonzeros in H \
            must be a non-negative scalar!");
    }

    if (!MX_IS_SCALAR(prhs[5]))
    {
        mexErrMsgTxt("qpPresolver (NEW): Initial stack size must be a positive scalar!");
    }

    /* EXTRACT DATA */
    m    = (qpp_int_t) (mxGetScalar(prhs[1]) + 0.5);
    n    = (qpp_int_t) (mxGetScalar(prhs[2]) + 0.5);
    Annz = (qpp_int_t) (mxGetScalar(prhs[3]) + 0.5);
    Hnnz = (qpp_int_t) (mxGetScalar(prhs[4]) + 0.5);

    exitflag = qppNew(&qppEntity[id], m, n, Annz, Hnnz, 
                      (qpp_int_t) (mxGetScalar(prhs[5]) + 0.5));
    
    if (exitflag != QPP_OK)
    {
        mexWarnMsgTxt("qpPresolver (NEW): qppNew() failed (cf. exitflag). Do not use ID!");
        id = -1;
    }

    plhs[0] = mxCreateDoubleScalar((double) id);
    plhs[1] = mxCreateDoubleScalar((double) exitflag);

    return;

QPPRESOLVER_INIT:
    /*	==============================================================================
    prhs[1] = ID of presolver entity
    prhs[2] = A (sparse)
    prhs[3] = H (sparse, only lower triangular part)
    prhs[4] = g (dense column vector)
    prhs[5] = f (scalar)
    prhs[6] = xl (dense column vector)
    prhs[7] = xu (dense column vector)
    prhs[8] = al (dense column vector)
    prhs[9] = au (dense column vector)
    plhs[0] = error code (scalar)
    ==================================================================================*/

    if (nrhs != 10)
    {
        mexErrMsgTxt("qpPresolver (INIT): Invalid number of input arguments!");
    }

    id = checkID(prhs[1]);
    if (id < 0)
    {
        mexErrMsgTxt("qpPresolver (INIT): Invalid ID!");
    }
    
    entity = qppEntity[id];
    m = entity->m;
    n = entity->n;

    for (i = 2; i <= 9; ++i)
    {
        if (!mxIsEmpty(prhs[i]))
        {
            if (!mxIsDouble(prhs[i]) || mxIsComplex(prhs[i]) )
            {
                mexErrMsgTxt("qpPresolver (INIT): Invalid (complex or non-double) input!");
            }
        }
    }

    /* EXTRACT SPARSE MATRICES */
    if ( (!mxIsEmpty(prhs[2]) && !mxIsSparse(prhs[2])) ||
         (!mxIsEmpty(prhs[3]) && !mxIsSparse(prhs[3])) )
    {
        mexErrMsgTxt("qpPresolver (INIT): Input matrices must be sparse (or empty)!");
    }

    Annz = mxGetNzmax(prhs[2]);
    Airn = (qpp_int_t*) mxGetIr(prhs[2]);
    Acp  = (qpp_int_t*) mxGetJc(prhs[2]);
    Ax   = mxGetPr(prhs[2]);

    Hnnz = mxGetNzmax(prhs[3]);
    Hirn = (qpp_int_t*) mxGetIr(prhs[3]);
    Hcp  = (qpp_int_t*) mxGetJc(prhs[3]);
    Hx   = mxGetPr(prhs[3]);
    
    /*  We have to check if a matrix is empty (the maximum number of nonzeros will
     *  always be > 0). */
    if ( (Annz > 0) && (Ax[0] == 0.0) )
    {
        Annz = 0;
    }
    if ( (Hnnz > 0) && (Hx[0] == 0.0) )
    {
        Hnnz = 0;
    }
    
    Ajcn = (qpp_int_t*) malloc(Annz * sizeof(qpp_int_t));
    Hjcn = (qpp_int_t*) malloc(Hnnz * sizeof(qpp_int_t));

    if ( ((Annz > 0) && (Ajcn == NULL)) || ((Hnnz > 0) && (Hjcn == NULL)) )
    {
        free(Ajcn);
        free(Hjcn);
        exitflag = QPP_OUT_OF_MEMORY;
        plhs[0] = mxCreateDoubleScalar((double) exitflag);
        return;
    }

    /* Create column indices for matrices in coordinate format */
    for (i = 0; i < n; ++i)
        for (j = Acp[i]; j < Acp[i+1]; ++j)
            Ajcn[j] = i;

    for (i = 0; i < n; ++i)
        for (j = Hcp[i]; j < Hcp[i+1]; ++j)
            Hjcn[j] = i;

    /* EXTRACT DENSE VECTORS */
    if ( (mxGetN(prhs[4]) != 1) || (mxGetM(prhs[4]) != (size_t)(n)) )
    {
        mexErrMsgTxt("qpPresolver (INIT): g must be a column vector of length n!");
    }
    g = mxGetData(prhs[4]);

	if (mxIsEmpty(prhs[5]))
	{
	    f = 0.0;
	}
	else
	{
        if (!MX_IS_SCALAR(prhs[5]))
        {
            mexErrMsgTxt("qpPresolver (INIT): f must be a scalar!");
        }
        f = mxGetScalar(prhs[5]);
    }

    if (mxIsEmpty(prhs[6]))
    {
        xl = NULL;
    }
    else
    {
        if ( (mxGetN(prhs[6]) != 1) || (mxGetM(prhs[6]) != (size_t)(n)) )
        {
            mexErrMsgTxt("qpPresolver (INIT): xl must be a column vector of length n!");
        }
        xl = mxGetData(prhs[6]);
    }

    if (mxIsEmpty(prhs[7]))
    {
        xu = NULL;
    }
    else
    {
        if ( (mxGetN(prhs[7]) != 1) || (mxGetM(prhs[7]) != (size_t)(n)) )
        {
            mexErrMsgTxt("qpPresolver (INIT): xu must be a column vector of length n!");
        }
        xu = mxGetData(prhs[7]);
    }

    al = NULL;
    au = NULL;
    if (m > 0)
    {
        if (!mxIsEmpty(prhs[8]))
        {
            if ( (mxGetN(prhs[8]) != 1) || (mxGetM(prhs[8]) != (size_t)(m)) )
            {
                mexErrMsgTxt("qpPresolver (INIT): al must be a column vector \
                    of length m!");
            }
            al = mxGetData(prhs[8]);
        }

        if (!mxIsEmpty(prhs[9]))
        {
            if ( (mxGetN(prhs[9]) != 1) || (mxGetM(prhs[9]) != (size_t)(m)) )
            {
                mexErrMsgTxt("qpPresolver (INIT): au must be a column vector \
                    of length m!");
            }
            au = mxGetData(prhs[9]);
        }
    }

    exitflag = qppInit(entity, Airn, Ajcn, Ax, Annz, Hirn, Hjcn, Hx, Hnnz, 
                       g, f, xl, xu, al, au);

    if (Annz > 0)
    {
        free(Ajcn);
    }
    
    if (Hnnz > 0)
    {
        free(Hjcn);
    }

    plhs[0] = mxCreateDoubleScalar((double) exitflag);

    return;

QPPRESOLVER_PRESOLVE:
    /*	==============================================================================
    prhs[1] = ID of presolver entity
    plhs[0] = error code
    plhs[1] = number of iterations of the presolver
    ==================================================================================*/
    if (nrhs != 2)
    {
        mexErrMsgTxt("qpPresolver (PRESOLVE): Invalid number of input arguments!");
    }

    id = checkID(prhs[1]);
    if (id < 0)
    {
        mexErrMsgTxt("qpPresolver (PRESOLVE): Invalid ID!");
    }

    exitflag = qppPresolve(qppEntity[id]);

    plhs[0] = mxCreateDoubleScalar((double) exitflag);
    plhs[1] = mxCreateDoubleScalar((double) qppEntity[id]->num_iter);

    return;

QPPRESOLVER_GETPRESOLVEDQP:
    /*	==============================================================================
    prhs[1]  = ID of presolver entity
    plhs[0]  = error code
    plhs[1]  = m (number of constraints, scalar)
    plhs[2]  = n (number of variables, scalar)
    plhs[3]  = Airn (integer column vector)
    plhs[4]  = Ajcn (integer column vector)
    plhs[5]  = Ax (real column vector)
    plhs[6]  = Hirn (integer column vector)
    plhs[7]  = Hjcn (integer column vector)
    plhs[8]  = Hx (real column vector)
    plhs[9]  = g (real column vector)
    plhs[10] = f (scalar)
    plhs[11] = xl (real column vector)
    plhs[12] = xu (real column vector)
    plhs[13] = al (real column vector)
    plhs[14] = au (real column vector)
    ==================================================================================*/

    if (nrhs != 2)
    {
        mexErrMsgTxt("qpPresolver (GETPRESOLVEDQP): Invalid number \
            of input arguments!");
    }

    if (nlhs != 15)
    {
        mexErrMsgTxt("qpPresolver (GETPRESOLVEDQP): Invalid number \
        of output arguments!");
    }

    id = checkID(prhs[1]);
    if (id < 0)
    {
        mexErrMsgTxt("qpPresolver (GETPRESOLVEDQP): Invalid ID!");
    }

    entity = qppEntity[id];
    exitflag = qppGetDimensionPresolvedQP(entity, &m, &n, &Annz, &Hnnz);
    
    if (exitflag != QPP_OK)
    {
        plhs[0] = mxCreateDoubleScalar((double) exitflag);
        for (i = 1; i < 15; ++i)
        {
            plhs[i] = mxCreateDoubleMatrix(0, 0, mxREAL);
        }
        return;
    }
    
    /* Allocate memory */
    plhs[1] = mxCreateDoubleScalar(m);
    plhs[2] = mxCreateDoubleScalar(n);
    plhs[3] = mxCreateNumericMatrix(Annz, 1, intclass, mxREAL);
    plhs[4] = mxCreateNumericMatrix(Annz, 1, intclass, mxREAL);
    plhs[5] = mxCreateDoubleMatrix(Annz, 1, mxREAL);
    plhs[6] = mxCreateNumericMatrix(Hnnz, 1, intclass, mxREAL);
    plhs[7] = mxCreateNumericMatrix(Hnnz, 1, intclass, mxREAL);
    plhs[8] = mxCreateDoubleMatrix(Hnnz, 1, mxREAL);
    plhs[9] = mxCreateDoubleMatrix(n, 1, mxREAL);
    plhs[10] = mxCreateDoubleScalar(0.0);
    plhs[11] = mxCreateDoubleMatrix(n, 1, mxREAL);
    plhs[12] = mxCreateDoubleMatrix(n, 1, mxREAL);
    plhs[13] = mxCreateDoubleMatrix(m, 1, mxREAL);
    plhs[14] = mxCreateDoubleMatrix(m, 1, mxREAL);
    
    Airn = mxGetData(plhs[3]);
    Ajcn = mxGetData(plhs[4]);
    Ax   = mxGetPr(plhs[5]);
    Hirn = mxGetData(plhs[6]);
    Hjcn = mxGetData(plhs[7]);
    Hx   = mxGetPr(plhs[8]);
    g    = mxGetPr(plhs[9]);
    xl   = mxGetPr(plhs[11]);
    xu   = mxGetPr(plhs[12]);
    al   = mxGetPr(plhs[13]);
    au   = mxGetPr(plhs[14]);
    
    exitflag = qppGetPresolvedQP(entity, &m, &n, Airn, Ajcn, Ax, &Annz, Hirn, Hjcn, 
                                 Hx, &Hnnz, g, mxGetPr(plhs[10]), xl, xu, al, au,
                                 QPP_MST_COLUMN_WISE);
    
    plhs[0] = mxCreateDoubleScalar((double) exitflag);

    if (exitflag != QPP_OK)
    {
        /* Return empty matrices as output if function qppGetPresolvedQP failed */
        for (i = 1; i < 15; ++i)
        {
            mxDestroyArray(plhs[i]);
            plhs[i] = mxCreateDoubleMatrix(0, 0, mxREAL);
        }
    }

    return;

QPPRESOLVER_POSTSOLVE:
    /*	==============================================================================
    prhs[1] = ID of presolver entity
    prhs[2] = Primal solution (x) of presolved problem
    prhs[3] = Dual solution (y) of presolved problem
    prhs[4] = Dual solution (z) of presolved problem
    prhs[5] = Working set for linear constraints (optional)
    prhs[6] = Working set for bound constraints (optional)
    plhs[0] = error code
    plhs[1] = Primal solution (x) of original problem
    plhs[2] = Dual solution (y) of original problem
    plhs[3] = Dual solution (z) of original problem
    plhs[4] = Optimal working set regarding linear constraints of original problem (optional)
    plhs[5] = Optimal working set regarding bound constraints of original problem (optional)
    ==================================================================================*/
    	
    if (nrhs < 7)
    {
        mexErrMsgTxt("qpPresolver (POSTSOLVE): Invalid number of input arguments!");
    }

    id = checkID(prhs[1]);
    if (id < 0)
    {
        mexErrMsgTxt("qpPresolver (POSTSOLVE): Invalid ID!");
    }

    entity = qppEntity[id];
    m = entity->m;
    n = entity->n;
    
    exitflag = qppGetDimensionPresolvedQP(entity, &pm, &pn, NULL, NULL);
    if (exitflag != QPP_OK)
    {
        plhs[0] = mxCreateDoubleScalar((double) exitflag);
        for (i = 1; i <= 3; ++i)
        {
            plhs[i] = mxCreateDoubleMatrix(0, 0, mxREAL);
        }
        return;
    }

    /* Check primal solution */
    if ( !mxIsEmpty(prhs[2]) && ( !MX_IS_DOUBLE(prhs[2]) ||
         (mxGetM(prhs[2]) != (size_t)(pn)) || (mxGetN(prhs[2]) != 1) ) )
    {
        mexErrMsgTxt("qpPresolver (POSTSOLVE): Primal solution (x) of \
            presolved QP must be a (real) column vector!");
    }

    if ( mxIsEmpty(prhs[2]) && (pn > 0) )
    {
        mexErrMsgTxt("qpPresolver (POSTSOLVE): Primal solution (x) of \
            presolved QP must be a (real) column vector!");
    }

    /* Check dual solution y */
    if ( !mxIsEmpty(prhs[3]) && ( !MX_IS_DOUBLE(prhs[3]) ||
         (mxGetM(prhs[3]) != (size_t)(pm)) || (mxGetN(prhs[3]) != 1) ) )
    {
        mexErrMsgTxt("qpPresolver (POSTSOLVE): Dual solution (y) of \
            presolved QP must be a (real) column vector!");
    }

    if ( mxIsEmpty(prhs[3]) && (pm > 0) )
    {
        mexErrMsgTxt("qpPresolver (POSTSOLVE): Dual solution (y) of \
            presolved QP must be a (real) column vector!");
    }

    /* Check dual solution z */
    if ( !mxIsEmpty(prhs[4]) && ( !MX_IS_DOUBLE(prhs[4]) ||
         (mxGetM(prhs[4]) != (size_t)(pn)) || (mxGetN(prhs[4]) != 1) ) )
    {
        mexErrMsgTxt("qpPresolver (POSTSOLVE): Dual solution (z) of \
            presolved QP must be a (real) column vector!");
    }

    if ( mxIsEmpty(prhs[4]) && (pn > 0) )
    {
        mexErrMsgTxt("qpPresolver (POSTSOLVE): Dual solution (z) of \
            presolved QP must be a (real) column vector!");
    }
    wi = NULL;    
    if (!mxIsEmpty(prhs[5]))
    {
        if (!MX_IS_DOUBLE(prhs[5]) || (mxGetM(prhs[5]) != (size_t)(pm)) || (mxGetN(prhs[5]) != 1) )
        {
            mexErrMsgTxt("qpPresolver (POSTSOLVE): Working set for linear constraints of \
                presolved QP must be a (real) column vector!");
        }

        /* Copy working set from input (double) to wi (int). */
        wi = (qpp_int_t*) mxMalloc(m * sizeof(qpp_int_t));
        Ax = mxGetPr(prhs[5]);
        for (i = 0; i < pm; ++i)
        {
            wi[i] = (qpp_int_t)(Ax[i]);
        }
    }
    else
    {
        if ( (pm == 0) && (m > 0) )
        {
            wi = (qpp_int_t*) mxMalloc(m * sizeof(qpp_int_t));
        }
    }

    wj = NULL;
    if (!mxIsEmpty(prhs[6]))
    {
        if ( !MX_IS_DOUBLE(prhs[6]) || (mxGetM(prhs[6]) != (size_t)(pn)) || (mxGetN(prhs[6]) != 1) )
        {
            mexErrMsgTxt("qpPresolver (POSTSOLVE): Working set for bound constraints of \
                presolved QP must be a (real) column vector!");
        }

        /* Copy working set from input (double) to wj (int). */
        wj = (qpp_int_t*) mxMalloc(n * sizeof(qpp_int_t));
        Ax = mxGetPr(prhs[6]);
        for (i = 0; i < pn; ++i)
        {
            wj[i] = (qpp_int_t)(Ax[i]);
        }
    }
    else
    {
        if ( (pn == 0) && (n > 0) )
        {
            wj = (qpp_int_t*) mxMalloc(n * sizeof(qpp_int_t));
        }
    }

    x = mxGetPr(prhs[2]);
    y = mxGetPr(prhs[3]);
    z = mxGetPr(prhs[4]);
  
    if (wi != NULL)
    {
        plhs[4] = mxCreateDoubleMatrix(m, 1, mxREAL);
    }
    else
    {
        plhs[4] = mxCreateDoubleMatrix(0, 0, mxREAL);
    }
    
    if (wj != NULL)
    {
        plhs[5] = mxCreateDoubleMatrix(n, 1, mxREAL);
    }
    else
    {
        plhs[5] = mxCreateDoubleMatrix(0, 0, mxREAL);
    }
    
    /* Allocate memory for x, y and z. */
    plhs[1] = mxCreateDoubleMatrix(n, 1, mxREAL);
    plhs[2] = mxCreateDoubleMatrix(m, 1, mxREAL);
    plhs[3] = mxCreateDoubleMatrix(n, 1, mxREAL);
   
    /* Copy solution of presolved QP. */
    memcpy(mxGetPr(plhs[1]), x, pn * sizeof(double));
    memcpy(mxGetPr(plhs[2]), y, pm * sizeof(double));
    memcpy(mxGetPr(plhs[3]), z, pn * sizeof(double));

    exitflag = qppPostsolve(entity, mxGetPr(plhs[1]), mxGetPr(plhs[2]), mxGetPr(plhs[3]), 
                            wi, wj, 0);

    plhs[0] = mxCreateDoubleScalar((double) exitflag);
    
    /* Convert wi/wj (int) into double again. */
    if (wi != NULL)
    {
        Ax = mxGetPr(plhs[4]);
        for (i = 0; i < m; ++i)
        {
            Ax[i] = (double) wi[i];
        }
    }
    
    if (wj != NULL)
    {
        Ax = mxGetPr(plhs[5]);
        for (i = 0; i < n; ++i)
        {
            Ax[i] = (double) wj[i];
        }
    }

    mxFree(wi);
    mxFree(wj);

    return;

QPPRESOLVER_GETSOLUTION:
    
    return;     /* Deprecated */

QPPRESOLVER_FREE:
    /*	================================================================================
    prhs[1] = ID of presolver entity
    ==================================================================================*/
    if (nrhs != 2)
    {
        mexErrMsgTxt("qpPresolver (FREE): Invalid number of input arguments!");
    }

    id = checkID(prhs[1]);
    if (id < 0)
    {
        mexErrMsgTxt("qpPresolver (FREE): Invalid ID!");
    }

    qppFree(&qppEntity[id]);
    qppEntity[id] = NULL;

    return;

QPPRESOLVER_SETOPTIONS:
    /*	================================================================================
    prhs[1] = ID of presolver entity
    prhs[2] = Options struct
    plhs[0] = error code
    ==================================================================================*/
    if (nrhs != 3)
    {
        mexErrMsgTxt("qpPresolver (SETOPTIONS): Invalid number of input arguments!");
    }
    
    id = checkID(prhs[1]);
    if (id < 0)
    {
        mexErrMsgTxt("qpPresolver (SETOPTIONS): Invalid ID!");
    }
    entity = qppEntity[id];
    
    if (mxIsEmpty(prhs[2]))
    {
        /* Default options */
        exitflag = qppSetOptions(entity, NULL);
    }
    else
    {
        exitflag = QPP_OK;
        if (!mxIsStruct(prhs[2]))
        {
            mexErrMsgTxt("qpPresolver (SETOPTIONS): Options must be passed as a struct \
                          cf. qppSetOptions() for further details)!");
        }
        else
        {
            setupOptions(entity->options, prhs[2]);
            
            if (qppCheckOptions(entity->options) != QPP_OK)
            {
                mexWarnMsgTxt("qpPresolver (SETOPTIONS): Some option parameters were \
                              invalid and have been replaced by default values!");
            }
        }
    }
    
    plhs[0] = mxCreateDoubleScalar((double) exitflag);

    return;
}

/* ========================= STATIC FUNCTIONS ========================== */

static qpp_int_t checkID(const mxArray* const ID)
{
    qpp_int_t id;
    
    if (!MX_IS_DOUBLE(ID) || !MX_IS_SCALAR(ID))
    {
        return -1;
    }
    
    id = (qpp_int_t) (mxGetScalar(ID) + 0.5);
    
    if ( (id < 0) || (id >= QPP_MAX_NUM_PROBLEMS) || (qppEntity[id] == NULL) )
    {
        return -1;
    }
    
    return id;
}


static qpp_return_value_t setupOptions(qpp_options_t *opt, const mxArray *mat_opt)
{
    qpp_int_t i;
    char *str;
    double *value;
    
    if (mxGetNumberOfFields(mat_opt) != 15)
    {
        mexWarnMsgTxt("Options might be set incorrectly as struct has wrong \
                      number of entries!");
    }
    
    if (hasOptionsValueC(mat_opt, "boundMode", &str) == QPP_OK)
    {
        for (i = 0; i < (qpp_int_t)(strlen(str)); ++i)
        {
            str[i] = tolower(str[i]);
        }
        
        if (strstr(str, "medium") != NULL)
        {
            opt->bound_mode = QPP_BM_MEDIUM_BOUNDS;
        }
        else if (strstr(str, "tight") != NULL)
        {
            opt->bound_mode = QPP_BM_TIGHTEST_BOUNDS;
        }
        else
        {
            mexWarnMsgTxt("qpPresolver (SETOPTIONS): Unknown bound type! Using default \
                           value instead!");
        }
        mxFree(str);
    }
    
    if (hasOptionsValueD(mat_opt, "equalityTol", &value) == QPP_OK)
    {
        opt->eq_tol = *value;
    }
    
    if (hasOptionsValueD(mat_opt, "stabilityTol", &value) == QPP_OK)
    {
        opt->stab_tol = *value;
    }
    
    if (hasOptionsValueD(mat_opt, "feasibilityTol", &value) == QPP_OK)
    {
        opt->feas_tol = *value;
    }
    
    if (hasOptionsValueD(mat_opt, "maxIter", &value) == QPP_OK)
    {
        opt->max_iter = (qpp_int_t) *value;
    }
    
    if (hasOptionsValueD(mat_opt, "logfileLevel", &value) == QPP_OK)
    {
        opt->log_level = (qpp_int_t) *value;
    }
    
    if (hasOptionsValueD(mat_opt, "enableBoundTightening", &value) == QPP_OK)
    {
        opt->enable_bound_tightening = (qpp_int_t) *value;
    }
    
    if (hasOptionsValueD(mat_opt, "enableDualConstraintsMethod", &value) == QPP_OK)
    {
        opt->enable_dual_constraints_method = (qpp_int_t) *value;
    }
    
    if (hasOptionsValueD(mat_opt, "enableDuplicateColumnsMethod", &value) == QPP_OK)
    {
        opt->enable_duplicate_columns_method = (qpp_int_t) *value;
    }
    
    if (hasOptionsValueD(mat_opt, "enableEmptyColumnsMethod", &value) == QPP_OK)
    {
        opt->enable_empty_columns_method = (qpp_int_t) *value;
    }
    
    if (hasOptionsValueD(mat_opt, "enablePrimalConstraintsMethod", &value) == QPP_OK)
    {
        opt->enable_primal_constraints_method = (qpp_int_t) *value;
    }
    
    if (hasOptionsValueD(mat_opt, "enableScaling", &value) == QPP_OK)
    {
        opt->enable_scaling = (qpp_int_t) *value;
    }
    
    if (hasOptionsValueD(mat_opt, "enableSingletonColumnsMethod", &value) == QPP_OK)
    {
        opt->enable_singleton_columns_method = (qpp_int_t) *value;
    }
    
    if (hasOptionsValueD(mat_opt, "enableSingletonRowsMethod", &value) == QPP_OK)
    {
        opt->enable_singleton_rows_method = (qpp_int_t) *value;
    }
    
    if (hasOptionsValueD(mat_opt, "enableSparsificationMethod", &value) == QPP_OK)
    {
        opt->enable_sparsification_method = (qpp_int_t) *value;
    }
    return QPP_OK;
}


static qpp_return_value_t hasOptionsValueD(const mxArray *opt,
                                           const char* const opt_str,
                                           double **opt_val)
{
    mxArray *opt_field;
    
    opt_field = mxGetField(opt, 0, opt_str);
    
    if (opt_field == NULL)
    {
        char msg[QPP_MAX_STRING_LENGTH];
        if (strlen(opt_str) <= QPP_MAX_STRING_LENGTH-100)
        {
            sprintf(msg, "qpPresolver (SETOPTIONS): Option struct does not contain \
                    entry %s!", opt_str);
        }
        else
        {
            sprintf(msg, "qpPresolver (SETOPTIONS): Unknown option struct entry!");
        }
        mexWarnMsgTxt(msg);
        return QPP_INVALID_ARGUMENT;
    }
    
    if ( !mxIsEmpty(opt_field) && MX_IS_SCALAR(opt_field) )
    {
        *opt_val = mxGetPr(opt_field);
        return QPP_OK;
    }
    else
    {
        char msg[QPP_MAX_STRING_LENGTH];
        if (strlen(opt_str) <= QPP_MAX_STRING_LENGTH-100)
        {
            sprintf(msg, "qpPresolver (SETOPTIONS): Option %s is not a scalar! \
                    Using default value instead!", opt_str);
        }
        else
        {
            sprintf(msg, "qpPresolver (SETOPTIONS): Some option value is not a scalar! \
                    Using default value instead!");
        }
        mexWarnMsgTxt(msg);
        return QPP_INVALID_ARGUMENT;
    }
}


static qpp_return_value_t hasOptionsValueC(const mxArray *opt,
                                           const char* const opt_str,
                                           char **opt_val)
{
    mxArray *opt_field;
    
    opt_field = mxGetField(opt, 0, opt_str);
    
    if (opt_field == NULL)
    {
        char msg[QPP_MAX_STRING_LENGTH];
        if (strlen(opt_str) <= QPP_MAX_STRING_LENGTH-100)
        {
            sprintf(msg, "qpPresolver (SETOPTIONS): Option struct does not contain \
                    entry %s!", opt_str);
        }
        else
        {
            sprintf(msg, "qpPresolver (SETOPTIONS): Unknown option struct entry!");
        }
        mexWarnMsgTxt(msg);
        return QPP_INVALID_ARGUMENT;
    }
    
    if ( !mxIsEmpty(opt_field) && mxIsChar(opt_field) )
    {
        *opt_val = mxArrayToString(opt_field);
        return QPP_OK;
    }
    else
    {
        char msg[QPP_MAX_STRING_LENGTH];
        if (strlen(opt_str) <= QPP_MAX_STRING_LENGTH-100)
        {
            sprintf(msg, "qpPresolver (SETOPTIONS): Option %s is not a string! \
                    Using default value instead!", opt_str);
        }
        else
        {
            sprintf(msg, "qpPresolver (SETOPTIONS): Some option value is not a string! \
                    Using default value instead!");
        }
        mexWarnMsgTxt(msg);
        return QPP_INVALID_ARGUMENT;
    }
}
