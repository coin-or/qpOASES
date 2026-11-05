#cython: language_level=3
#distutils: language=c++

##
##  This file is part of qpOASES.
##
##  qpOASES -- An Implementation of the Online Active Set Strategy.
##  Copyright (C) 2007-2017 by Hans Joachim Ferreau, Andreas Potschka,
##  Christian Kirches et al. All rights reserved.
##
##  qpOASES is free software; you can redistribute it and/or
##  modify it under the terms of the GNU Lesser General Public
##  License as published by the Free Software Foundation; either
##  version 2.1 of the License, or (at your option) any later version.
##
##  qpOASES is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
##  See the GNU Lesser General Public License for more details.
##
##  You should have received a copy of the GNU Lesser General Public
##  License along with qpOASES; if not, write to the Free Software
##  Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
##

## authors of this file: Sebastian F. Walter (thanks to Felix Lenders)

"""
Python interface to qpOASES
using Cython
:author: Sebastian F. Walter, Manuel Kudruss (thanks to Felix Lenders)
"""

import warnings
warnings.simplefilter("once", DeprecationWarning)
import numpy as np
cimport numpy as np
import scipy as sp

from scipy.sparse import SparseEfficiencyWarning

from cython.operator cimport dereference as deref

from libcpp.memory cimport unique_ptr, make_unique


cimport qpoases


def deprecation_warning_nWSR():
    warnings.warn("\nInteger nWSR will be deprecated in qpOASES 4.0.\nUse nWSR = numpy.array([10]) as input to qp.init() and qp.hotstart()", DeprecationWarning, stacklevel=2)

def deprecation_warning_cputime():
    warnings.warn("\nFloat cputime will be deprecated in qpOASES 4.0.\nUse cputime = numpy.array([2.0]) as input to qp.init() and qp.hotstart()", DeprecationWarning, stacklevel=2)


real_t_type = np.dtype('double')
int_t_type = np.dtype('long')
sparse_int_t_type = np.dtype('long')


cdef str __version__ = "3.2.2"


cdef BooleanType bt_from_bool(bint value):
    return BT_TRUE if value else BT_FALSE


cdef bint bt_to_bool(BooleanType value):
    if value == BT_TRUE:
        return True
    assert value == BT_FALSE
    return False


cdef np.ndarray get_cputime(cputime):
    if isinstance(cputime, float) and cputime != 0:
        deprecation_warning_cputime()

    if isinstance(cputime, float):
        cputime = np.array([cputime], dtype=real_t_type)

    if cputime.shape!= (1,) or cputime.item() < 0:
        raise ValueError("Invalid cputime provided")

    return <np.ndarray> cputime


cdef real_t* get_cputime_view(np.ndarray cputime_arr):
    if cputime_arr.item() == 0.:
        return NULL

    return <real_t*> cputime_arr.data


cdef real_t* get_vec_view(np.ndarray vec):
    assert vec.ndim == 1
    if vec.size == 0:
        return <real_t*> NULL
    return <real_t*> vec.data


cdef sparse_int_t* get_index_view(np.ndarray vec):
    assert vec.ndim == 1
    assert vec.dtype == sparse_int_t_type
    if vec.size == 0:
        return <sparse_int_t*> NULL
    return <sparse_int_t*> vec.data


cdef real_t* get_mat_view(np.ndarray mat):
    assert mat.ndim == 2
    if mat.size == 0:
        return <real_t*> NULL
    return <real_t*> mat.data


cdef np.ndarray get_nWSR(nWSR):
    if isinstance(nWSR, int):
        deprecation_warning_nWSR()
        nWSR = np.array([nWSR], dtype=int_t_type)

    if nWSR.shape != (1,) or nWSR.item() < 0:
        raise ValueError("Invalid nWSR provided")

    return <np.ndarray> nWSR


cdef class PyPrintLevel:
    DEBUG_ITER = PL_DEBUG_ITER
    TABULAR    = PL_TABULAR
    NONE       = PL_NONE
    LOW        = PL_LOW
    MEDIUM     = PL_MEDIUM
    HIGH       = PL_HIGH

cpdef enum PyHessianType:
    ZERO               = HST_ZERO
    IDENTITY           = HST_IDENTITY
    POSDEF             = HST_POSDEF
    POSDEF_NULLSPACE   = HST_POSDEF_NULLSPACE
    SEMIDEF            = HST_SEMIDEF
    INDEF              = HST_INDEF
    UNKNOWN            = HST_UNKNOWN

cpdef enum PySubjectToStatus:
    LOWER              = ST_LOWER
    INACTIVE           = ST_INACTIVE
    UPPER              = ST_UPPER
    INFEASIBLE_LOWER   = ST_INFEASIBLE_LOWER
    INFEASIBLE_UPPER   = ST_INFEASIBLE_UPPER
    UNDEFINED          = ST_UNDEFINED

cdef class PyReturnValue:
    TERMINAL_LIST_ELEMENT                 = -1
    SUCCESSFUL_RETURN                     = 0
    DIV_BY_ZERO                           = RET_DIV_BY_ZERO
    INDEX_OUT_OF_BOUNDS                   = RET_INDEX_OUT_OF_BOUNDS
    INVALID_ARGUMENTS                     = RET_INVALID_ARGUMENTS
    ERROR_UNDEFINED                       = RET_ERROR_UNDEFINED
    WARNING_UNDEFINED                     = RET_WARNING_UNDEFINED
    INFO_UNDEFINED                        = RET_INFO_UNDEFINED
    EWI_UNDEFINED                         = RET_EWI_UNDEFINED
    AVAILABLE_WITH_LINUX_ONLY             = RET_AVAILABLE_WITH_LINUX_ONLY
    UNKNOWN_BUG                           = RET_UNKNOWN_BUG
    PRINTLEVEL_CHANGED                    = RET_PRINTLEVEL_CHANGED
    NOT_YET_IMPLEMENTED                   = RET_NOT_YET_IMPLEMENTED
    INDEXLIST_MUST_BE_REORDERD            = RET_INDEXLIST_MUST_BE_REORDERD
    INDEXLIST_EXCEEDS_MAX_LENGTH          = RET_INDEXLIST_EXCEEDS_MAX_LENGTH
    INDEXLIST_CORRUPTED                   = RET_INDEXLIST_CORRUPTED
    INDEXLIST_OUTOFBOUNDS                 = RET_INDEXLIST_OUTOFBOUNDS
    INDEXLIST_ADD_FAILED                  = RET_INDEXLIST_ADD_FAILED
    INDEXLIST_INTERSECT_FAILED            = RET_INDEXLIST_INTERSECT_FAILED
    INDEX_ALREADY_OF_DESIRED_STATUS       = RET_INDEX_ALREADY_OF_DESIRED_STATUS
    ADDINDEX_FAILED                       = RET_ADDINDEX_FAILED
    REMOVEINDEX_FAILED                    = RET_REMOVEINDEX_FAILED
    SWAPINDEX_FAILED                      = RET_SWAPINDEX_FAILED
    NOTHING_TO_DO                         = RET_NOTHING_TO_DO
    SETUP_BOUND_FAILED                    = RET_SETUP_BOUND_FAILED
    SETUP_CONSTRAINT_FAILED               = RET_SETUP_CONSTRAINT_FAILED
    MOVING_BOUND_FAILED                   = RET_MOVING_BOUND_FAILED
    MOVING_CONSTRAINT_FAILED              = RET_MOVING_CONSTRAINT_FAILED
    SHIFTING_FAILED                       = RET_SHIFTING_FAILED
    ROTATING_FAILED                       = RET_ROTATING_FAILED
    QPOBJECT_NOT_SETUP                    = RET_QPOBJECT_NOT_SETUP
    QP_ALREADY_INITIALISED                = RET_QP_ALREADY_INITIALISED
    NO_INIT_WITH_STANDARD_SOLVER          = RET_NO_INIT_WITH_STANDARD_SOLVER
    RESET_FAILED                          = RET_RESET_FAILED
    INIT_FAILED                           = RET_INIT_FAILED
    INIT_FAILED_TQ                        = RET_INIT_FAILED_TQ
    INIT_FAILED_CHOLESKY                  = RET_INIT_FAILED_CHOLESKY
    INIT_FAILED_HOTSTART                  = RET_INIT_FAILED_HOTSTART
    INIT_FAILED_INFEASIBILITY             = RET_INIT_FAILED_INFEASIBILITY
    INIT_FAILED_UNBOUNDEDNESS             = RET_INIT_FAILED_UNBOUNDEDNESS
    INIT_FAILED_REGULARISATION            = RET_INIT_FAILED_REGULARISATION
    INIT_SUCCESSFUL                       = RET_INIT_SUCCESSFUL
    OBTAINING_WORKINGSET_FAILED           = RET_OBTAINING_WORKINGSET_FAILED
    SETUP_WORKINGSET_FAILED               = RET_SETUP_WORKINGSET_FAILED
    SETUP_AUXILIARYQP_FAILED              = RET_SETUP_AUXILIARYQP_FAILED
    NO_CHOLESKY_WITH_INITIAL_GUESS        = RET_NO_CHOLESKY_WITH_INITIAL_GUESS
    NO_EXTERN_SOLVER                      = RET_NO_EXTERN_SOLVER
    QP_UNBOUNDED                          = RET_QP_UNBOUNDED
    QP_INFEASIBLE                         = RET_QP_INFEASIBLE
    QP_NOT_SOLVED                         = RET_QP_NOT_SOLVED
    QP_SOLVED                             = RET_QP_SOLVED
    UNABLE_TO_SOLVE_QP                    = RET_UNABLE_TO_SOLVE_QP
    INITIALISATION_STARTED                = RET_INITIALISATION_STARTED
    HOTSTART_FAILED                       = RET_HOTSTART_FAILED
    HOTSTART_FAILED_TO_INIT               = RET_HOTSTART_FAILED_TO_INIT
    HOTSTART_FAILED_AS_QP_NOT_INITIALISED = RET_HOTSTART_FAILED_AS_QP_NOT_INITIALISED
    ITERATION_STARTED                     = RET_ITERATION_STARTED
    SHIFT_DETERMINATION_FAILED            = RET_SHIFT_DETERMINATION_FAILED
    STEPDIRECTION_DETERMINATION_FAILED    = RET_STEPDIRECTION_DETERMINATION_FAILED
    STEPLENGTH_DETERMINATION_FAILED       = RET_STEPLENGTH_DETERMINATION_FAILED
    OPTIMAL_SOLUTION_FOUND                = RET_OPTIMAL_SOLUTION_FOUND
    HOMOTOPY_STEP_FAILED                  = RET_HOMOTOPY_STEP_FAILED
    HOTSTART_STOPPED_INFEASIBILITY        = RET_HOTSTART_STOPPED_INFEASIBILITY
    HOTSTART_STOPPED_UNBOUNDEDNESS        = RET_HOTSTART_STOPPED_UNBOUNDEDNESS
    WORKINGSET_UPDATE_FAILED              = RET_WORKINGSET_UPDATE_FAILED
    MAX_NWSR_REACHED                      = RET_MAX_NWSR_REACHED
    CONSTRAINTS_NOT_SPECIFIED             = RET_CONSTRAINTS_NOT_SPECIFIED
    INVALID_FACTORISATION_FLAG            = RET_INVALID_FACTORISATION_FLAG
    UNABLE_TO_SAVE_QPDATA                 = RET_UNABLE_TO_SAVE_QPDATA
    STEPDIRECTION_FAILED_TQ               = RET_STEPDIRECTION_FAILED_TQ
    STEPDIRECTION_FAILED_CHOLESKY         = RET_STEPDIRECTION_FAILED_CHOLESKY
    CYCLING_DETECTED                      = RET_CYCLING_DETECTED
    CYCLING_NOT_RESOLVED                  = RET_CYCLING_NOT_RESOLVED
    CYCLING_RESOLVED                      = RET_CYCLING_RESOLVED
    STEPSIZE                              = RET_STEPSIZE
    STEPSIZE_NONPOSITIVE                  = RET_STEPSIZE_NONPOSITIVE
    SETUPSUBJECTTOTYPE_FAILED             = RET_SETUPSUBJECTTOTYPE_FAILED
    ADDCONSTRAINT_FAILED                  = RET_ADDCONSTRAINT_FAILED
    ADDCONSTRAINT_FAILED_INFEASIBILITY    = RET_ADDCONSTRAINT_FAILED_INFEASIBILITY
    ADDBOUND_FAILED                       = RET_ADDBOUND_FAILED
    ADDBOUND_FAILED_INFEASIBILITY         = RET_ADDBOUND_FAILED_INFEASIBILITY
    REMOVECONSTRAINT_FAILED               = RET_REMOVECONSTRAINT_FAILED
    REMOVEBOUND_FAILED                    = RET_REMOVEBOUND_FAILED
    REMOVE_FROM_ACTIVESET                 = RET_REMOVE_FROM_ACTIVESET
    ADD_TO_ACTIVESET                      = RET_ADD_TO_ACTIVESET
    REMOVE_FROM_ACTIVESET_FAILED          = RET_REMOVE_FROM_ACTIVESET_FAILED
    ADD_TO_ACTIVESET_FAILED               = RET_ADD_TO_ACTIVESET_FAILED
    CONSTRAINT_ALREADY_ACTIVE             = RET_CONSTRAINT_ALREADY_ACTIVE
    ALL_CONSTRAINTS_ACTIVE                = RET_ALL_CONSTRAINTS_ACTIVE
    LINEARLY_DEPENDENT                    = RET_LINEARLY_DEPENDENT
    LINEARLY_INDEPENDENT                  = RET_LINEARLY_INDEPENDENT
    LI_RESOLVED                           = RET_LI_RESOLVED
    ENSURELI_FAILED                       = RET_ENSURELI_FAILED
    ENSURELI_FAILED_TQ                    = RET_ENSURELI_FAILED_TQ
    ENSURELI_FAILED_NOINDEX               = RET_ENSURELI_FAILED_NOINDEX
    ENSURELI_FAILED_CYCLING               = RET_ENSURELI_FAILED_CYCLING
    BOUND_ALREADY_ACTIVE                  = RET_BOUND_ALREADY_ACTIVE
    ALL_BOUNDS_ACTIVE                     = RET_ALL_BOUNDS_ACTIVE
    CONSTRAINT_NOT_ACTIVE                 = RET_CONSTRAINT_NOT_ACTIVE
    BOUND_NOT_ACTIVE                      = RET_BOUND_NOT_ACTIVE
    HESSIAN_NOT_SPD                       = RET_HESSIAN_NOT_SPD
    HESSIAN_INDEFINITE                    = RET_HESSIAN_INDEFINITE
    MATRIX_SHIFT_FAILED                   = RET_MATRIX_SHIFT_FAILED
    MATRIX_FACTORISATION_FAILED           = RET_MATRIX_FACTORISATION_FAILED
    PRINT_ITERATION_FAILED                = RET_PRINT_ITERATION_FAILED
    NO_GLOBAL_MESSAGE_OUTPUTFILE          = RET_NO_GLOBAL_MESSAGE_OUTPUTFILE
    DISABLECONSTRAINTS_FAILED             = RET_DISABLECONSTRAINTS_FAILED
    ENABLECONSTRAINTS_FAILED              = RET_ENABLECONSTRAINTS_FAILED
    ALREADY_ENABLED                       = RET_ALREADY_ENABLED
    ALREADY_DISABLED                      = RET_ALREADY_DISABLED
    NO_HESSIAN_SPECIFIED                  = RET_NO_HESSIAN_SPECIFIED
    USING_REGULARISATION                  = RET_USING_REGULARISATION
    EPS_MUST_BE_POSITVE                   = RET_EPS_MUST_BE_POSITVE
    REGSTEPS_MUST_BE_POSITVE              = RET_REGSTEPS_MUST_BE_POSITVE
    HESSIAN_ALREADY_REGULARISED           = RET_HESSIAN_ALREADY_REGULARISED
    CANNOT_REGULARISE_IDENTITY            = RET_CANNOT_REGULARISE_IDENTITY
    CANNOT_REGULARISE_SPARSE              = RET_CANNOT_REGULARISE_SPARSE
    NO_REGSTEP_NWSR                       = RET_NO_REGSTEP_NWSR
    FEWER_REGSTEPS_NWSR                   = RET_FEWER_REGSTEPS_NWSR
    CHOLESKY_OF_ZERO_HESSIAN              = RET_CHOLESKY_OF_ZERO_HESSIAN
    ZERO_HESSIAN_ASSUMED                  = RET_ZERO_HESSIAN_ASSUMED
    CONSTRAINTS_ARE_NOT_SCALED            = RET_CONSTRAINTS_ARE_NOT_SCALED
    INITIAL_BOUNDS_STATUS_NYI             = RET_INITIAL_BOUNDS_STATUS_NYI
    ERROR_IN_CONSTRAINTPRODUCT            = RET_ERROR_IN_CONSTRAINTPRODUCT
    FIX_BOUNDS_FOR_LP                     = RET_FIX_BOUNDS_FOR_LP
    USE_REGULARISATION_FOR_LP             = RET_USE_REGULARISATION_FOR_LP
    UPDATEMATRICES_FAILED                 = RET_UPDATEMATRICES_FAILED
    UPDATEMATRICES_FAILED_AS_QP_NOT_SOLVED= RET_UPDATEMATRICES_FAILED_AS_QP_NOT_SOLVED
    UNABLE_TO_OPEN_FILE                   = RET_UNABLE_TO_OPEN_FILE
    UNABLE_TO_WRITE_FILE                  = RET_UNABLE_TO_WRITE_FILE
    UNABLE_TO_READ_FILE                   = RET_UNABLE_TO_READ_FILE
    FILEDATA_INCONSISTENT                 = RET_FILEDATA_INCONSISTENT
    UNABLE_TO_ANALYSE_QPROBLEM            = RET_UNABLE_TO_ANALYSE_QPROBLEM
    OPTIONS_ADJUSTED                      = RET_OPTIONS_ADJUSTED
    NWSR_SET_TO_ONE                       = RET_NWSR_SET_TO_ONE
    UNABLE_TO_READ_BENCHMARK              = RET_UNABLE_TO_READ_BENCHMARK
    BENCHMARK_ABORTED                     = RET_BENCHMARK_ABORTED
    INITIAL_QP_SOLVED                     = RET_INITIAL_QP_SOLVED
    QP_SOLUTION_STARTED                   = RET_QP_SOLUTION_STARTED
    BENCHMARK_SUCCESSFUL                  = RET_BENCHMARK_SUCCESSFUL
    NO_DIAGONAL_AVAILABLE                 = RET_NO_DIAGONAL_AVAILABLE
    DIAGONAL_NOT_INITIALISED              = RET_DIAGONAL_NOT_INITIALISED
    ENSURELI_DROPPED                      = RET_ENSURELI_DROPPED
    KKT_MATRIX_SINGULAR                   = RET_KKT_MATRIX_SINGULAR
    QR_FACTORISATION_FAILED               = RET_QR_FACTORISATION_FAILED
    INERTIA_CORRECTION_FAILED             = RET_INERTIA_CORRECTION_FAILED
    NO_SPARSE_SOLVER                      = RET_NO_SPARSE_SOLVER
    SIMPLE_STATUS_P1                      = RET_SIMPLE_STATUS_P1
    SIMPLE_STATUS_P0                      = RET_SIMPLE_STATUS_P0
    SIMPLE_STATUS_M1                      = RET_SIMPLE_STATUS_M1
    SIMPLE_STATUS_M2                      = RET_SIMPLE_STATUS_M2
    SIMPLE_STATUS_M3                      = RET_SIMPLE_STATUS_M3


class QPOASESError(Exception):
    def __init__(self, error_message, error_code):
        self.error_code = PyReturnValue(error_code)
        message = f"QPOASES error <{str(error_code)}>: {error_message}"
        super().__init__(message)


cdef check_return_value(returnValue retval):
    cdef unique_ptr[MessageHandling] message_handling
    if retval != SUCCESSFUL_RETURN:
        message_handling = unique_ptr[MessageHandling](new MessageHandling())
        message = deref(message_handling).getErrorCodeMessage(retval).decode("utf-8")
        raise QPOASESError(message, retval)


def check_shape(arr, name, shape):
    if isinstance(arr, np.ndarray) and not arr.flags.c_contiguous:
        raise ValueError("Array must be c-continuous")

    if arr.shape != shape:
        raise ValueError("Array {0} has invalid shape (expected: {1}, actual: {2})".format(name,
                                                                                           shape,
                                                                                           arr.shape))


cdef class PySymmMat:
    cdef dict __dict__
    cdef unique_ptr[SymmetricMatrix] mat_ptr

    def __cinit__(self, H: np.ndarray | sp.sparse.spmatrix):
        cdef SymSparseMat* sparse_ptr = NULL
        assert H.ndim == 2

        if isinstance(H, np.ndarray):
            H = np.atleast_2d(H)
            if H.ndim != 2:
                raise ValueError("Invalid matrix shape")

            (m, n) = H.shape[:2]

            if m != n:
                raise ValueError("Matrix must be square")

            if not H.flags.c_contiguous:
                raise ValueError("Array must be c-continuous")

            self.H = H

            self.mat_ptr = unique_ptr[SymmetricMatrix](new SymDenseMat(m,
                                                                       n,
                                                                       n,
                                                                       get_mat_view(self.H)))

        else:
            if not sp.sparse.issparse(H):
                raise ValueError("Unknown matrix type")

            if H.format != "csc":
                warnings.warn(
                    "Matrix is required to be CSC",
                    SparseEfficiencyWarning,
                )
                H = H.tocsc()

            (m, n) = H.shape[:2]

            if m != n:
                raise ValueError("Matrix must be square")

            self.H_rows = H.indices.astype(sparse_int_t_type, copy=False)
            self.H_cols = H.indptr.astype(sparse_int_t_type, copy=False)
            self.H_data = H.data.astype(real_t_type, copy=False)

            sparse_ptr = new SymSparseMat(m,
                                          n,
                                          get_index_view(self.H_rows),
                                          get_index_view(self.H_cols),
                                          get_vec_view(self.H_data))

            deref(sparse_ptr).createDiagInfo()
            self.mat_ptr = unique_ptr[SymmetricMatrix](sparse_ptr)

        assert not deref(self.mat_ptr).needToFreeMemory()

    def print(self, name):
        return deref(self.mat_ptr).print(name.encode("utf-8"))


cdef class PyMat:
    cdef dict __dict__
    cdef unique_ptr[Matrix] mat_ptr

    def __cinit__(self, A: np.ndarray | sp.sparse.spmatrix):
        cdef int_t m, n, l
        cdef real_t[:, :] A_view

        if isinstance(A, np.ndarray):
            self.A = np.atleast_2d(A)
            if A.ndim != 2:
                raise ValueError("Invalid matrix shape")

            (m, n) = A.shape[:2]

            if not A.flags.c_contiguous:
                raise ValueError("Array must be c-continuous")

            self.mat_ptr = unique_ptr[Matrix](new DenseMatrix(m,
                                                              n,
                                                              n,
                                                              get_mat_view(self.A)))

        else:
            if not sp.sparse.issparse(A):
                raise ValueError("Unknown matrix type")

            if A.format != "csc":
                warnings.warn(
                    "Matrix is required to be CSC",
                    SparseEfficiencyWarning,
                )
                A = A.tocsc()

            (m, n) = A.shape[:2]

            self.A_rows = A.indices.astype(sparse_int_t_type, copy=False)
            self.A_cols = A.indptr.astype(sparse_int_t_type, copy=False)
            self.A_data = A.data.astype(real_t_type, copy=False)

            self.mat_ptr = unique_ptr[Matrix](new SparseMatrix(m,
                                                               n,
                                                               get_index_view(self.A_rows),
                                                               get_index_view(self.A_cols),
                                                               get_vec_view(self.A_data)))

        assert not deref(self.mat_ptr).needToFreeMemory()

    def print(self, name):
        return deref(self.mat_ptr).print(name.encode("utf-8"))


cdef get_working_set_bounds(QProblemB* problem):
    cdef np.ndarray working_set_bounds = np.empty((problem.getNV(),),
                                                  dtype=real_t_type)

    check_return_value(problem.getWorkingSetBounds(<real_t*> working_set_bounds.data))

    return working_set_bounds


cdef get_working_set_constraints(QProblem* problem):
    cdef np.ndarray working_set_bounds = np.empty((problem.getNC(),),
                                                  dtype=real_t_type)

    check_return_value(problem.getWorkingSetConstraints(<real_t*> working_set_bounds.data))

    return working_set_bounds


cdef get_bounds(QProblemB* problem):
        cdef PyBounds bounds = PyBounds()
        cdef Bounds* bounds_view = bounds.thisptr.get()

        check_return_value(problem.getBounds(deref(bounds_view)))

        return bounds

cdef get_constraints(QProblem* problem):
        cdef PyConstraints constraints = PyConstraints()
        cdef Constraints* constraints_view = constraints.thisptr.get()

        check_return_value(problem.getConstraints(deref(constraints_view)))

        return constraints


cdef class PyOptions:
    cdef Options *thisptr      # hold a C++ instance which we're wrapping
    def __cinit__(self):
        # FIXME: add support for the other constructors
        self.thisptr = new Options()

    def __dealloc__(self):
        del self.thisptr

    def setToDefault(self):
        return self.thisptr.setToDefault()

    def setToReliable(self):
        return self.thisptr.setToReliable()

    def setToMPC(self):
        return self.thisptr.setToMPC()

    def setToFast(self):
        return self.thisptr.setToFast()

    def ensureConsistency(self):
        return self.thisptr.ensureConsistency()

    property printLevel:
        def __get__(self): return self.thisptr.printLevel
        def __set__(self, printLevel): self.thisptr.printLevel = printLevel

    property enableRamping:
        def __get__(self): return bt_to_bool(self.thisptr.enableRamping)
        def __set__(self, bint enableRamping): self.thisptr.enableRamping = bt_from_bool(enableRamping)

    property enableFarBounds:
        def __get__(self): return bt_to_bool(self.thisptr.enableFarBounds)
        def __set__(self, bint enableFarBounds): self.thisptr.enableFarBounds = bt_from_bool(enableFarBounds)

    property enableFlippingBounds:
        def __get__(self): return bt_to_bool(self.thisptr.enableFlippingBounds)
        def __set__(self, bint enableFlippingBounds): self.thisptr.enableFlippingBounds = bt_from_bool(enableFlippingBounds)

    property enableRegularisation:
        def __get__(self): return bt_to_bool(self.thisptr.enableRegularisation)
        def __set__(self, bint enableRegularisation): self.thisptr.enableRegularisation = bt_from_bool(enableRegularisation)

    property enableFullLITests:
        def __get__(self): return bt_to_bool(self.thisptr.enableFullLITests)
        def __set__(self, bint enableFullLITests): self.thisptr.enableFullLITests = bt_from_bool(enableFullLITests)

    property enableNZCTests:
        def __get__(self): return bt_to_bool(self.thisptr.enableNZCTests)
        def __set__(self, bint enableNZCTests): self.thisptr.enableNZCTests = bt_from_bool(enableNZCTests)

    property enableDriftCorrection:
        def __get__(self): return self.thisptr.enableDriftCorrection
        def __set__(self, enableDriftCorrection): self.thisptr.enableDriftCorrection = enableDriftCorrection

    property enableCholeskyRefactorisation:
        def __get__(self): return self.thisptr.enableCholeskyRefactorisation
        def __set__(self, enableCholeskyRefactorisation): self.thisptr.enableCholeskyRefactorisation = enableCholeskyRefactorisation

    property enableEqualities:
        def __get__(self): return bt_to_bool(self.thisptr.enableEqualities)
        def __set__(self, bint enableEqualities): self.thisptr.enableEqualities = bt_from_bool(enableEqualities)

    property terminationTolerance:
        def __get__(self): return self.thisptr.terminationTolerance
        def __set__(self, terminationTolerance): self.thisptr.terminationTolerance = terminationTolerance

    property boundTolerance:
        def __get__(self): return self.thisptr.boundTolerance
        def __set__(self, boundTolerance): self.thisptr.boundTolerance = boundTolerance

    property boundRelaxation:
        def __get__(self): return self.thisptr.boundRelaxation
        def __set__(self, boundRelaxation): self.thisptr.boundRelaxation = boundRelaxation

    property epsNum:
        def __get__(self): return self.thisptr.epsNum
        def __set__(self, epsNum): self.thisptr.epsNum = epsNum

    property epsDen:
        def __get__(self): return self.thisptr.epsDen
        def __set__(self, epsDen): self.thisptr.epsDen = epsDen

    property maxPrimalJump:
        def __get__(self): return self.thisptr.maxPrimalJump
        def __set__(self, maxPrimalJump): self.thisptr.maxPrimalJump = maxPrimalJump

    property maxDualJump:
        def __get__(self): return self.thisptr.maxDualJump
        def __set__(self, maxDualJump): self.thisptr.maxDualJump = maxDualJump

    property initialRamping:
        def __get__(self): return self.thisptr.initialRamping
        def __set__(self, initialRamping): self.thisptr.initialRamping = initialRamping

    property finalRamping:
        def __get__(self): return self.thisptr.finalRamping
        def __set__(self, finalRamping): self.thisptr.finalRamping = finalRamping

    property initialFarBounds:
        def __get__(self): return self.thisptr.initialFarBounds
        def __set__(self, initialFarBounds): self.thisptr.initialFarBounds = initialFarBounds

    property growFarBounds:
        def __get__(self): return self.thisptr.growFarBounds
        def __set__(self, growFarBounds): self.thisptr.growFarBounds = growFarBounds

    property initialStatusBounds:
        def __get__(self): return self.thisptr.initialStatusBounds
        def __set__(self, initialStatusBounds): self.thisptr.initialStatusBounds = initialStatusBounds

    property epsFlipping:
        def __get__(self): return self.thisptr.epsFlipping
        def __set__(self, epsFlipping): self.thisptr.epsFlipping = epsFlipping

    property numRegularisationSteps:
        def __get__(self): return self.thisptr.numRegularisationSteps
        def __set__(self, numRegularisationSteps): self.thisptr.numRegularisationSteps = numRegularisationSteps

    property epsRegularisation:
        def __get__(self): return self.thisptr.epsRegularisation
        def __set__(self, epsRegularisation): self.thisptr.epsRegularisation = epsRegularisation

    property numRefinementSteps:
        def __get__(self): return self.thisptr.numRefinementSteps
        def __set__(self, numRefinementSteps): self.thisptr.numRefinementSteps = numRefinementSteps

    property epsIterRef:
        def __get__(self): return self.thisptr.epsIterRef
        def __set__(self, epsIterRef): self.thisptr.epsIterRef = epsIterRef

    property epsLITests:
        def __get__(self): return self.thisptr.epsLITests
        def __set__(self, epsLITests): self.thisptr.epsLITests = epsLITests

    property epsNZCTests:
        def __get__(self): return self.thisptr.epsNZCTests
        def __set__(self, epsNZCTests): self.thisptr.epsNZCTests = epsNZCTests

    property dropBoundPriority:
        def __get__(self): return self.thisptr.dropBoundPriority
        def __set__(self, dropBoundPriority): self.thisptr.dropBoundPriority = dropBoundPriority

    property dropEqConPriority:
        def __get__(self): return self.thisptr.dropEqConPriority
        def __set__(self, dropEqConPriority): self.thisptr.dropEqConPriority = dropEqConPriority

    property dropIneqConPriority:
        def __get__(self): return self.thisptr.dropIneqConPriority
        def __set__(self, dropIneqConPriority): self.thisptr.dropIneqConPriority = dropIneqConPriority


cdef class PyConstraints:
    cdef unique_ptr[Constraints] thisptr

    def __cinit__(self, n=None):
        if n is not None:
            self.thisptr = make_unique[Constraints](<int_t> n)
        else:
            self.thisptr = make_unique[Constraints]()

    def init(self, n):
        check_return_value(deref(self.thisptr).init(<int_t> n))

    def setupConstraint(self,
                   int_t number,
                   PySubjectToStatus status):
        check_return_value(deref(self.thisptr).setupConstraint(number, <SubjectToStatus> status))

    def setupAllInactive(self):
        check_return_value(deref(self.thisptr).setupAllInactive())

    def setupAllLower(self):
        check_return_value(deref(self.thisptr).setupAllLower())

    def setupAllUpper(self):
        check_return_value(deref(self.thisptr).setupAllUpper())

    def moveActiveToInactive(self, int_t number):
        check_return_value(deref(self.thisptr).moveActiveToInactive(number))

    def moveInactiveToActive(self, int_t number, PySubjectToStatus status):
        check_return_value(deref(self.thisptr).moveInactiveToActive(number,
                                                                    <SubjectToStatus> status))

    def print(self):
        check_return_value(deref(self.thisptr).print())

    def getActive(self):
        cdef PyIndexlist activeList = PyIndexlist()
        activeList.thisptr = deref(self.thisptr).getActive()
        activeList.parent = self
        return activeList

    def getInactive(self):
        cdef PyIndexlist inactiveList = PyIndexlist()
        inactiveList.thisptr = deref(self.thisptr).getInactive()
        inactiveList.parent = self
        return inactiveList

    def getStatus(self, int_t i):
        return PySubjectToStatus(deref(self.thisptr).getStatus(i))

    def setStatus(self, int_t i, PySubjectToStatus value):
        check_return_value(deref(self.thisptr).setStatus(i, <SubjectToStatus> value))

    def copy(self):
        cdef PyConstraints copied = PyConstraints()
        copied.thisptr = make_unique[Constraints](deref(self.thisptr.get()))
        return copied

    def __copy__(self):
        return self.copy()



cdef class PyIndexlist:
    cdef Indexlist* thisptr
    cdef object parent

    def __cinit__(self):
        self.thisptr = NULL

    def getIndex(self, int_t number):
        return self.thisptr.getIndex(number)

    def getNumber(self, int_t index):
        return self.thisptr.getNumber(index)

    @property
    def length(self):
        return self.this.ptr.getLength()

    @property
    def lastNumber(self):
        return self.thisptr.getLastNumber()

    def addNumber(self, int_t number):
        check_return_value(self.thisptr.addNumber(number))

    def removeNumber(self, int_t number):
        check_return_value(self.thisptr.removeNumber(number))

    def swapNumbers(self, int_t first, int_t second):
        check_return_value(self.thisptr.swapNumbers(first, second))

    def isMember(self, int_t number):
        return self.thisptr.isMember(number)

    def __contains__(self, number):
        return self.isMember(number)


cdef class PyBounds:
    cdef unique_ptr[Bounds] thisptr

    def __cinit__(self, n=None):
        if n is not None:
            self.thisptr = make_unique[Bounds](<int_t> n)
        else:
            self.thisptr = make_unique[Bounds]()

    def init(self, n):
        check_return_value(deref(self.thisptr).init(<int_t> n))

    def setupBound(self,
                   int_t number,
                   PySubjectToStatus status):
        check_return_value(deref(self.thisptr).setupBound(number, <SubjectToStatus> status))

    def setupAllFree(self):
        check_return_value(deref(self.thisptr).setupAllFree())

    def setupAllLower(self):
        check_return_value(deref(self.thisptr).setupAllLower())

    def setupAllUpper(self):
        check_return_value(deref(self.thisptr).setupAllUpper())

    def print(self):
        check_return_value(deref(self.thisptr).print())

    def moveFixedToFree(self, int_t number):
        check_return_value(deref(self.thisptr).moveFixedToFree(number))

    def moveFreeToFixed(self, int_t number, PySubjectToStatus status):
        check_return_value(deref(self.thisptr).moveFreeToFixed(number,
                                                               <SubjectToStatus> status))

    def getFree(self):
        cdef PyIndexlist freeList = PyIndexlist()
        freeList.thisptr = deref(self.thisptr).getFree()
        freeList.parent = self
        return freeList

    def getFixed(self):
        cdef PyIndexlist fixedList = PyIndexlist()
        fixedList.thisptr = deref(self.thisptr).getFixed()
        fixedList.parent = self
        return fixedList

    def getStatus(self, int_t i):
        return PySubjectToStatus(deref(self.thisptr).getStatus(i))

    def setStatus(self, int_t i, PySubjectToStatus value):
        check_return_value(deref(self.thisptr).setStatus(i, <SubjectToStatus> value))

    def copy(self):
        cdef PyBounds copied = PyBounds()
        copied.thisptr = make_unique[Bounds](deref(self.thisptr.get()))
        return copied

    def __copy__(self):
        return self.copy()

cdef class PyQProblemB:
    cdef unique_ptr[QProblemB] thisptr      # hold a C++ instance which we're wrapping
    cdef PySymmMat H

    def __cinit__(self,
                  long nV,
                  PyHessianType hessian_type=PyHessianType.UNKNOWN):
        self.thisptr = make_unique[QProblemB](<int_t> nV, <HessianType> hessian_type, BT_TRUE)

    def init(self,
             H: np.ndarray | sp.sparse.spmatrix,
             np.ndarray[np.double_t, ndim=1] g,
             np.ndarray[np.double_t, ndim=1] lb,
             np.ndarray[np.double_t, ndim=1] ub,
             nWSR,
             cputime = 0.0,
             np.ndarray[np.double_t, ndim=1] x_opt=None,
             np.ndarray[np.double_t, ndim=1] y_opt=None,
             PyBounds guessed_bounds=None
        ):
        # FIXME: add asserts
        cdef np.ndarray nWSR_tmp
        cdef np.ndarray cput_tmp
        cdef real_t* x_opt_view = NULL
        cdef real_t* y_opt_view = NULL
        cdef Bounds* guessed_bounds_view = NULL
        cdef real_t* cput_view = NULL

        NV = self.NV

        check_shape(H, "H", (NV, NV))
        check_shape(g, "g", (NV,))
        check_shape(lb, "lb", (NV,))
        check_shape(ub, "ub", (NV,))

        if x_opt is not None:
            check_shape(x_opt, "x_opt", (NV,))
            x_opt_view = get_vec_view(x_opt)

        if y_opt is not None:
            check_shape(y_opt, "y_opt", (NV,))
            y_opt_view = get_vec_view(y_opt)

        if guessed_bounds is not None:
            guessed_bounds.thisptr.get()

        nWSR_tmp = get_nWSR(nWSR)

        self.H = PySymmMat(H)

        cput_tmp = get_cputime(cputime)
        cput_view = get_cputime_view(cput_tmp)

        check_return_value(deref(self.thisptr).init(
            self.H.mat_ptr.get(),
            get_vec_view(g),
            get_vec_view(lb),
            get_vec_view(ub),
            <int_t&>  nWSR_tmp.data[0],
            cput_view,
            x_opt_view,
            y_opt_view,
            guessed_bounds_view))

    def hotstart(self,
             np.ndarray[np.double_t, ndim=1] g,
             np.ndarray[np.double_t, ndim=1] lb,
             np.ndarray[np.double_t, ndim=1] ub,
             nWSR,
             cputime=0.0,
             PyBounds guessed_bounds=None
        ):
        # FIXME: add asserts
        cdef np.ndarray nWSR_tmp
        # cdef np.ndarray cput_tmp
        cdef Bounds* guessed_bounds_view = NULL
        cdef np.ndarray cput_tmp
        cdef real_t* cput_view = NULL

        NV = self.NV
        check_shape(g, "g", (NV,))
        check_shape(lb, "lb", (NV,))
        check_shape(ub, "ub", (NV,))

        if guessed_bounds is not None:
            guessed_bounds_view = guessed_bounds.thisptr.get()

        nWSR_tmp = get_nWSR(nWSR)

        cput_tmp = get_cputime(cputime)
        cput_view = get_cputime_view(cput_tmp)

        check_return_value(deref(self.thisptr).hotstart(
                get_vec_view(g),
                get_vec_view(lb),
                get_vec_view(ub),
                <int_t&>  nWSR_tmp.data[0],
                cput_view,
                guessed_bounds_view))

    cpdef bint isSolved(self):
        return deref(self.thisptr).isSolved()

    cpdef bint isInfeasible(self):
        return deref(self.thisptr).isInfeasible()

    cpdef bint isUnbounded(self):
        return deref(self.thisptr).isUnbounded()

    def getPrimalSolution(self, np.ndarray[np.double_t, ndim=1] xOpt):
        return deref(self.thisptr).getPrimalSolution(<real_t*> xOpt.data)

    def getDualSolution(self, np.ndarray[np.double_t, ndim=1] yOpt):
        NV = self.NV
        check_shape(yOpt, "yOpt", (NV,))
        check_return_value(deref(self.thisptr).getDualSolution(<real_t*> yOpt.data))

    def getObjVal(self):
        return deref(self.thisptr).getObjVal()

    def printOptions(self):
        return deref(self.thisptr).printOptions()

    def getOptions(self):
        # FIXME: memory management? who deallocates o
        cdef Options *o = new Options(deref(self.thisptr).getOptions())
        retval = PyOptions()
        retval.thisptr = o
        return retval

    def setOptions(self, PyOptions options):
        check_return_value(deref(self.thisptr).setOptions(deref(options.thisptr)))

    @property
    def NV(self) -> int:
        return deref(self.thisptr).getNV()

    @property
    def NFR(self) -> int:
        return deref(self.thisptr).getNFR()

    @property
    def NFX(self) -> int:
        return deref(self.thisptr).getNFX()

    @property
    def NFV(self) -> int:
        return deref(self.thisptr).getNFV()

    @property
    def NZ(self) -> int:
        return deref(self.thisptr).getNZ()

    def getBounds(self) -> PyBounds:
        return get_bounds(self.thisptr.get())

    def getWorkingSetBounds(self):
        return get_working_set_bounds(self.thisptr.get())


cdef class PyQProblem:
    cdef unique_ptr[QProblem] thisptr      # hold a C++ instance which we're wrapping
    cdef PySymmMat H
    cdef PyMat A

    def __cinit__(self,
                  long nV,
                  long nC,
                  PyHessianType hessian_type=PyHessianType.UNKNOWN):
        self.thisptr = make_unique[QProblem](nV, nC, <HessianType> hessian_type, BT_TRUE)

    cpdef init(self,
             H: np.ndarray | sp.sparse.spmatrix,
             np.ndarray[np.double_t, ndim=1] g,
             A: np.ndarray | sp.sparse.spmatrix,
             np.ndarray[np.double_t, ndim=1] lb,
             np.ndarray[np.double_t, ndim=1] ub,
             np.ndarray[np.double_t, ndim=1] lbA,
             np.ndarray[np.double_t, ndim=1] ubA,
             nWSR,
             cputime=0.0,
             np.ndarray x_opt=None,
             np.ndarray y_opt=None,
             PyBounds guessed_bounds=None,
             PyConstraints guessed_constraints=None):

        # FIXME: add asserts
        cdef np.ndarray nWSR_tmp
        cdef np.ndarray cput_tmp
        cdef real_t* cput_view = NULL
        cdef real_t* x_opt_view = NULL
        cdef real_t* y_opt_view = NULL
        cdef Bounds* guessed_bounds_view = NULL
        cdef Constraints* guessed_constraints_view = NULL

        NV = self.NV
        NC = self.NC

        check_shape(H, "H", (NV, NV))
        check_shape(A, "A", (NC, NV))
        check_shape(g, "g", (NV,))
        check_shape(lb, "lb", (NV,))
        check_shape(ub, "ub", (NV,))
        
        check_shape(lbA, "lbA", (NC,))
        check_shape(ubA, "ubA", (NC,))

        nWSR_tmp = get_nWSR(nWSR)

        self.H = PySymmMat(H)
        self.A = PyMat(A)

        if x_opt is not None:
            check_shape(x_opt, "x_opt", (NV,))
            x_opt_view = get_vec_view(x_opt)

        if y_opt is not None:
            check_shape(y_opt, "y_opt", (NV + NC,))
            y_opt_view = get_vec_view(y_opt)

        if guessed_bounds is not None:
            guessed_bounds_view = guessed_bounds.thisptr.get()

        if guessed_constraints is not None:
            guessed_constraints_view = guessed_constraints.thisptr.get()

        cput_tmp = get_cputime(cputime)
        cput_view = get_cputime_view(cput_tmp)

        check_return_value(deref(self.thisptr).init(
                    self.H.mat_ptr.get(),
                    get_vec_view(g),
                    self.A.mat_ptr.get(),
                    get_vec_view(lb),
                    get_vec_view(ub),
                    get_vec_view(lbA),
                    get_vec_view(ubA),
                    <int_t&>  nWSR_tmp.data[0],
                    cput_view,
                    x_opt_view,
                    y_opt_view,
                    guessed_bounds_view,
                    guessed_constraints_view))

    cpdef hotstart(self,
             np.ndarray[np.double_t, ndim=1] g,
             np.ndarray[np.double_t, ndim=1] lb,
             np.ndarray[np.double_t, ndim=1] ub,
             np.ndarray[np.double_t, ndim=1] lbA,
             np.ndarray[np.double_t, ndim=1] ubA,
             nWSR,
             cputime=0.0,
             PyBounds guessed_bounds=None,
             PyConstraints guessed_constraints=None
        ):
        # FIXME: add asserts
        cdef np.ndarray nWSR_tmp
        cdef np.ndarray cput_tmp
        cdef real_t* cput_view = NULL

        cdef Bounds* guessed_bounds_view = NULL
        cdef Constraints* guessed_constraints_view = NULL

        NV = self.NV
        NC = self.NC

        check_shape(g, "g", (NV,))
        check_shape(lb, "lb", (NV,))
        check_shape(ub, "ub", (NV,))
        
        check_shape(lbA, "lbA", (NC,))
        check_shape(ubA, "ubA", (NC,))

        if guessed_bounds is not None:
            guessed_bounds_view = guessed_bounds.thisptr.get()

        if guessed_constraints is not None:
            guessed_constraints_view = guessed_constraints.thisptr.get()

        nWSR_tmp = get_nWSR(nWSR)

        cput_tmp = get_cputime(cputime)
        cput_view = get_cputime_view(cput_tmp)

        check_return_value(deref(self.thisptr).hotstart(
                    get_vec_view(g),
                    get_vec_view(lb),
                    get_vec_view(ub),
                    get_vec_view(lbA),
                    get_vec_view(ubA),
                    <int_t&>  nWSR_tmp.data[0],
                    cput_view, #<real_t*> &cput_tmp.data[0],
                    guessed_bounds_view,
                    guessed_constraints_view))

    cpdef bint isSolved(self):
        return deref(self.thisptr).isSolved()

    cpdef bint isInfeasible(self):
        return deref(self.thisptr).isInfeasible()

    cpdef bint isUnbounded(self):
        return deref(self.thisptr).isUnbounded()

    cpdef getPrimalSolution(self, np.ndarray[np.double_t, ndim=1] xOpt):
        check_return_value(deref(self.thisptr).getPrimalSolution(<real_t*> xOpt.data))

    cpdef getDualSolution(self, np.ndarray[np.double_t, ndim=1] yOpt):
        NV = self.NV
        NC = self.NC
        check_shape(yOpt, "yOpt", (NV + NC,))
        check_return_value(deref(self.thisptr).getDualSolution(<real_t*> yOpt.data))

    cpdef getObjVal(self):
        return deref(self.thisptr).getObjVal()

    cpdef printOptions(self):
        check_return_value(deref(self.thisptr).printOptions())

    cpdef setOptions(self, PyOptions options):
        check_return_value(deref(self.thisptr).setOptions(deref(options.thisptr)))

    @property
    def NV(self):
        return deref(self.thisptr).getNV()

    @property
    def NFR(self):
        return deref(self.thisptr).getNFR()

    @property
    def NFX(self):
        return deref(self.thisptr).getNFX()

    @property
    def NFV(self):
        return deref(self.thisptr).getNFV()

    @property
    def NZ(self):
        return deref(self.thisptr).getNZ()

    @property
    def NC(self):
        return deref(self.thisptr).getNC()

    @property
    def NEC(self):
        return deref(self.thisptr).getNEC()

    @property
    def NAC(self):
        return deref(self.thisptr).getNAC()

    @property
    def NIAC(self):
        return deref(self.thisptr).getNIAC()

    def getBounds(self) -> PyBounds:
        return get_bounds(self.thisptr.get())

    def getConstraints(self) -> PyConstraints:
        return get_constraints(self.thisptr.get())

    def getWorkingSetBounds(self):
        return get_working_set_bounds(self.thisptr.get())

    def getWorkingSetConstraints(self):
        return get_working_set_constraints(self.thisptr.get())


cdef class PySQProblem:
    cdef unique_ptr[SQProblem] thisptr      # hold a C++ instance which we're wrapping
    cdef PySymmMat H
    cdef PyMat A

    def __cinit__(self,
                  int_t nV,
                  int_t nC,
                  PyHessianType hessian_type=PyHessianType.UNKNOWN):
        self.thisptr = make_unique[SQProblem](nV, nC, <HessianType> hessian_type, BT_TRUE)

    cpdef init(self,
             H: np.ndarray | sp.sparse.spmatrix,
             np.ndarray[np.double_t, ndim=1] g,
             A: np.ndarray | sp.sparse.spmatrix,
             np.ndarray[np.double_t, ndim=1] lb,
             np.ndarray[np.double_t, ndim=1] ub,
             np.ndarray[np.double_t, ndim=1] lbA,
             np.ndarray[np.double_t, ndim=1] ubA,
             nWSR,
             cputime=0.0,
             np.ndarray x_opt=None,
             np.ndarray y_opt=None,
             PyBounds guessed_bounds=None,
             PyConstraints guessed_constraints=None
    ):
        self._maybe_init(H, g, A, lb, ub, lbA, ubA, nWSR, cputime, x_opt, y_opt,
                         guessed_bounds, guessed_constraints, hotstart=False)

    cpdef _maybe_init(self,
             H: np.ndarray | sp.sparse.spmatrix,
             np.ndarray[np.double_t, ndim=1] g,
             A: np.ndarray | sp.sparse.spmatrix,
             np.ndarray[np.double_t, ndim=1] lb,
             np.ndarray[np.double_t, ndim=1] ub,
             np.ndarray[np.double_t, ndim=1] lbA,
             np.ndarray[np.double_t, ndim=1] ubA,
             nWSR,
             cputime,
             np.ndarray x_opt,
             np.ndarray y_opt,
             PyBounds guessed_bounds,
             PyConstraints guessed_constraints,
             hotstart
    ):
        # FIXME: add asserts
        cdef np.ndarray nWSR_tmp
        # cdef np.ndarray cput_tmp
        cdef real_t* x_opt_view = NULL
        cdef real_t* y_opt_view = NULL
        cdef Bounds* guessed_bounds_view = NULL
        cdef Constraints* guessed_constraints_view = NULL
        cdef np.ndarray cput_tmp
        cdef real_t* cput_view = NULL

        NV = self.NV
        NC = self.NC

        check_shape(H, "H", (NV, NV))
        check_shape(A, "A", (NC, NV))
        check_shape(g, "g", (NV,))
        check_shape(lb, "lb", (NV,))
        check_shape(ub, "ub", (NV,))
        
        check_shape(lbA, "lbA", (NC,))
        check_shape(ubA, "ubA", (NC,))

        nWSR_tmp = get_nWSR(nWSR)

        self.H = PySymmMat(H)
        self.A = PyMat(A)

        if x_opt is not None:
            check_shape(x_opt, "x_opt", (NV,))
            x_opt_view = get_vec_view(x_opt)

        if y_opt is not None:
            check_shape(y_opt, "y_opt", (NV + NC,))
            y_opt_view = get_vec_view(y_opt)

        if guessed_bounds is not None:
            guessed_bounds_view = guessed_bounds.thisptr.get()

        if guessed_constraints is not None:
            guessed_constraints_view = guessed_constraints.thisptr.get()

        cput_tmp = get_cputime(cputime)
        cput_view = get_cputime_view(cput_tmp)

        if hotstart:
            check_return_value(deref(self.thisptr).hotstart(
                        self.H.mat_ptr.get(),
                        get_vec_view(g),
                        self.A.mat_ptr.get(),
                        get_vec_view(lb),
                        get_vec_view(ub),
                        get_vec_view(lbA),
                        get_vec_view(ubA),
                        <int_t&>  nWSR_tmp.data[0],
                        cput_view,
                        guessed_bounds_view,
                        guessed_constraints_view))
        else:
            check_return_value(deref(self.thisptr).init(
                        self.H.mat_ptr.get(),
                        get_vec_view(g),
                        self.A.mat_ptr.get(),
                        get_vec_view(lb),
                        get_vec_view(ub),
                        get_vec_view(lbA),
                        get_vec_view(ubA),
                        <int_t&>  nWSR_tmp.data[0],
                        cput_view,
                        x_opt_view,
                        y_opt_view,
                        guessed_bounds_view,
                        guessed_constraints_view))

    cpdef hotstart(self,
             H: np.ndarray | sp.sparse.spmatrix,
             np.ndarray[np.double_t, ndim=1] g,
             A: np.ndarray | sp.sparse.spmatrix,
             np.ndarray[np.double_t, ndim=1] lb,
             np.ndarray[np.double_t, ndim=1] ub,
             np.ndarray[np.double_t, ndim=1] lbA,
             np.ndarray[np.double_t, ndim=1] ubA,
             nWSR,
             cputime=0.0,
             PyBounds guessed_bounds=None,
             PyConstraints guessed_constraints=None
    ):
        self._maybe_init(H, g, A, lb, ub, lbA, ubA, nWSR, cputime, None, None,
                         guessed_bounds, guessed_constraints, hotstart=True)

    cpdef bint isSolved(self):
        return deref(self.thisptr).isSolved()

    cpdef bint isInfeasible(self):
        return deref(self.thisptr).isInfeasible()

    cpdef bint isUnbounded(self):
        return deref(self.thisptr).isUnbounded()

    cpdef getPrimalSolution(self, np.ndarray[np.double_t, ndim=1] xOpt):
        NV = self.NV
        check_shape(xOpt, "xOpt", (NV,))
        return deref(self.thisptr).getPrimalSolution(<real_t*> xOpt.data)

    cpdef getDualSolution(self, np.ndarray[np.double_t, ndim=1] yOpt):
        NV = self.NV
        NC = self.NC
        check_shape(yOpt, "yOpt", (NV + NC,))
        check_return_value(deref(self.thisptr).getDualSolution(<real_t*> yOpt.data))

    cpdef getObjVal(self):
        return deref(self.thisptr).getObjVal()

    cpdef printOptions(self):
        return deref(self.thisptr).printOptions()

    cpdef setOptions(self, PyOptions options):
        check_return_value(deref(self.thisptr).setOptions(deref(options.thisptr)))

    @property
    def NV(self):
        return deref(self.thisptr).getNV()

    @property
    def NFR(self):
        return deref(self.thisptr).getNFR()

    @property
    def NFX(self):
        return deref(self.thisptr).getNFX()

    @property
    def NFV(self):
        return deref(self.thisptr).getNFV()

    @property
    def NZ(self):
        return deref(self.thisptr).getNZ()

    @property
    def NC(self):
        return deref(self.thisptr).getNC()

    @property
    def NEC(self):
        return deref(self.thisptr).getNEC()

    @property
    def NAC(self):
        return deref(self.thisptr).getNAC()

    @property
    def NIAC(self):
        return deref(self.thisptr).getNIAC()

    def getBounds(self) -> PyBounds:
        return get_bounds(self.thisptr.get())

    def getConstraints(self) -> PyConstraints:
        return get_constraints(self.thisptr.get())

    def getWorkingSetBounds(self):
        return get_working_set_bounds(self.thisptr.get())

    def getWorkingSetConstraints(self):
        return get_working_set_constraints(self.thisptr.get())


cdef class PySolutionAnalysis:
    cdef SolutionAnalysis *thisptr      # hold a C++ instance which we're wrapping
    def __cinit__(self):
        self.thisptr = new SolutionAnalysis()
    def __dealloc__(self):
        del self.thisptr

    cpdef getKktViolation(self, qp,
        np.ndarray[np.double_t, ndim=1] maxStat,
        np.ndarray[np.double_t, ndim=1] maxFeas,
        np.ndarray[np.double_t, ndim=1] maxCmpl
    ):
        """ """
        if isinstance(qp, PySQProblem):
            return self._getKktViolation_SQProblem(qp, maxStat, maxFeas, maxCmpl)

        elif isinstance(qp, PyQProblem):
            return self._getKktViolation_QProblem(qp, maxStat, maxFeas, maxCmpl)

        elif isinstance(qp, PyQProblemB):
            return self._getKktViolation_QProblemB(qp, maxStat, maxFeas, maxCmpl)

        else:
            raise ValueError('argument 1 must be QProblemB, QProblem or SQProblem')

    cpdef _getKktViolation_QProblemB(self,
            PyQProblemB qp,
            np.ndarray[real_t, ndim=1] maxStat,
            np.ndarray[real_t, ndim=1] maxFeas,
            np.ndarray[real_t, ndim=1] maxCmpl
        ):
        return self.thisptr.getKktViolation(
                qp.thisptr.get(),
                <real_t*> maxStat.data,
                <real_t*> maxFeas.data,
                <real_t*> maxCmpl.data
            )

    cpdef _getKktViolation_QProblem(self,
            PyQProblem qp,
            np.ndarray[real_t, ndim=1] maxStat,
            np.ndarray[real_t, ndim=1] maxFeas,
            np.ndarray[real_t, ndim=1] maxCmpl
        ):
        return self.thisptr.getKktViolation(
                qp.thisptr.get(),
                <real_t*> maxStat.data,
                <real_t*> maxFeas.data,
                <real_t*> maxCmpl.data
            )


    cpdef _getKktViolation_SQProblem(self,
            PySQProblem qp,
            np.ndarray[real_t, ndim=1] maxStat,
            np.ndarray[real_t, ndim=1] maxFeas,
            np.ndarray[real_t, ndim=1] maxCmpl
        ):
        return self.thisptr.getKktViolation(
                qp.thisptr.get(),
                <real_t*> maxStat.data,
                <real_t*> maxFeas.data,
                <real_t*> maxCmpl.data
            )

    cpdef getVarianceCovariance(self,
                              qp,
                              np.ndarray[np.double_t, ndim=1] g_b_bA_VAR,
                              np.ndarray[np.double_t, ndim=1] Primal_Dual_VAR ):

        if isinstance(qp, PyQProblemB):
            check_return_value(self._getVarianceCovariance_QProblemB(qp, g_b_bA_VAR, Primal_Dual_VAR))

        elif isinstance(qp, PyQProblem):
            check_return_value(self._getVarianceCovariance_QProblem(qp, g_b_bA_VAR, Primal_Dual_VAR))

        elif isinstance(qp, PySQProblem):
            check_return_value(self._getVarianceCovariance_SQProblem(qp, g_b_bA_VAR, Primal_Dual_VAR))
        else:
            raise ValueError('argument 1 must be QProblemB, QProblem or SQProblem')

    cdef returnValue _getVarianceCovariance_QProblemB(self,
                                                      PyQProblemB qp,
                                                      np.ndarray[np.double_t, ndim=1] g_b_bA_VAR,
                                                      np.ndarray[np.double_t, ndim=1] Primal_Dual_VAR ):
        return self.thisptr.getVarianceCovariance(qp.thisptr.get(),
                                                  <real_t*> g_b_bA_VAR.data,
                                                  <real_t*> Primal_Dual_VAR.data)

    cdef returnValue _getVarianceCovariance_QProblem(self,
                                                     PyQProblem qp,
                                                     np.ndarray[np.double_t, ndim=1] g_b_bA_VAR,
                                                     np.ndarray[np.double_t, ndim=1] Primal_Dual_VAR ):
        return self.thisptr.getVarianceCovariance(qp.thisptr.get(),
                                                  <real_t*> g_b_bA_VAR.data,
                                                  <real_t*> Primal_Dual_VAR.data)

    cdef returnValue _getVarianceCovariance_SQProblem(self,
                                                      PySQProblem qp,
                                                      np.ndarray[np.double_t, ndim=1] g_b_bA_VAR,
                                                      np.ndarray[np.double_t, ndim=1] Primal_Dual_VAR ):
        return self.thisptr.getVarianceCovariance(qp.thisptr.get(),
                                                  <real_t*> g_b_bA_VAR.data,
                                                  <real_t*> Primal_Dual_VAR.data)

# Wrapped some utility functions for unit testing
cpdef py_runOqpBenchmark(path,               # Full path of the benchmark files (without trailing slash!).
                         isSparse,           # Shall convert matrices to sparse format before solution?
                         useHotstarts,       # Shall QP solution be hotstarted?
                         PyOptions options,  # QP solver options to be used while solving benchmark problems.
                         long maxAllowedNWSR, # Maximum number of working set recalculations to be performed.
                         double maxCPUTime,  # Maximum allowed CPU time for qp solving.
                         ):
    """run a QP benchmark example"""
    maxNWSR            = 0.0 # Output: Maximum number of performed working set recalculations.
    avgNWSR            = 0.0 # Output: Average number of performed working set recalculations.
    maxCPUtime         = 0.0 # Output: Maximum CPU time required for solving each QP.
    avgCPUtime         = 0.0 # Output: Average CPU time required for solving each QP.
    maxStationarity    = 0.0 # Output: Maximum residual of stationarity condition.
    maxFeasibility     = 0.0 # Output: Maximum residual of primal feasibility condition.
    maxComplementarity = 0.0 # Output: Maximum residual of complementarity condition.

"""
def py_getKktViolation(long nV,                             # Number of variables.
                       long nC,                             # Number of constraints.
                       np.ndarray[np.double_t, ndim=2] H,   # Hessian matrix.
                       np.ndarray[np.double_t, ndim=1] g,   # Sequence of gradient vectors.
                       np.ndarray[np.double_t, ndim=2] A,   # Constraint matrix.
                       np.ndarray[np.double_t, ndim=1] lb,  # Sequence of lower bound vectors (on variables).
                       np.ndarray[np.double_t, ndim=1] ub,  # Sequence of upper bound vectors (on variables).
                       np.ndarray[np.double_t, ndim=1] lbA, # Sequence of lower constraints' bound vectors.
                       np.ndarray[np.double_t, ndim=1] ubA, # Sequence of upper constraints' bound vectors.
                       np.ndarray[np.double_t, ndim=1] x,   # Sequence of primal trial vectors.
                       np.ndarray[np.double_t, ndim=1] y,   # Sequence of dual trial vectors.
                       ):
    stat = 0.0 # Maximum value of stationarity condition residual.
    feas = 0.0 # Maximum value of primal feasibility violation.
    cmpl = 0.0 # Maximum value of complementarity residual.
    getKktViolation(nV,
                    nC,
                    <real_t*> H.data,
                    <real_t*> g.data,
                    <real_t*> A.data,
                    <real_t*> lb.data,
                    <real_t*> ub.data,
                    <real_t*> lbA.data,
                    <real_t*> ubA.data,
                    <real_t*> x.data,
                    <real_t*> y.data,
                    stat,
                    feas,
                    cmpl
                    )
    return stat, feas, cmpl
"""
