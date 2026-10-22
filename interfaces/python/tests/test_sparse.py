import numpy as np
import scipy as sp

from qpoases import PyQProblemB as QProblemB
from qpoases import PyHessianType as HessianType


def test_sparse_solve():
    # Objective: Track x0 in l2 norm subjecting
    # to variable bounds
    x0 = np.array([3., 2.])

    lb = np.array([0., 0.])
    ub = np.array([1, 1.])

    qp = QProblemB(2, HessianType.POSDEF)

    H = sp.sparse.eye_array(2, format="csc")

    nWSR = np.array([10])

    qp.init(H,
            -2*x0,
            lb,
            ub,
            nWSR)

    x_opt = np.empty_like(x0)
    y_opt = np.empty_like(x0)

    qp.getPrimalSolution(x_opt)
    qp.getDualSolution(y_opt)

    assert np.allclose(x_opt, [1., 1.])
