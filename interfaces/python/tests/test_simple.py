import numpy as np

from qpoases import PyQProblemB as QProblemB
from qpoases import PyQProblem as QProblem
from qpoases import PySubjectToStatus as SubjectToStatus


def test_bound_stationarity():
    # Lagrange convention: L(x, y) = f(x) - y^T c(x)
    num_vars = 2
    H = np.eye(2)
    x0 = np.array([3., 1.])
    g = -x0
    lb = np.array([-2., -2.])
    ub = np.array([2., 2.])

    prob = QProblemB(2)
    nWSR = np.array([10_000])

    prob.init(H, g, lb, ub, nWSR)

    x_opt = np.empty_like(x0)
    y_opt = np.empty_like(x0)

    prob.getPrimalSolution(x_opt)
    prob.getDualSolution(y_opt)

    at_lower = (lb + 1e-6 >= x_opt)
    at_upper = (x_opt + 1e-6 >= ub)

    bounds = prob.getBounds()

    for i in range(num_vars):
        if at_lower[i] and at_upper[i]:
            assert bounds.getStatus(i) in [SubjectToStatus.LOWER,
                                           SubjectToStatus.UPPER]
        elif at_lower[i]:
            assert bounds.getStatus(i) == SubjectToStatus.LOWER
        elif at_upper[i]:
            assert bounds.getStatus(i) == SubjectToStatus.UPPER
        else:
            assert bounds.getStatus(i) == SubjectToStatus.INACTIVE

    active = at_lower | at_upper
    inactive = ~active

    assert np.allclose(y_opt[inactive], 0.)
    assert (y_opt[at_upper] <= 0.).all()
    assert (y_opt[at_lower] >= 0.).all()

    assert np.allclose(H @ x_opt + g - y_opt, 0.)


def test_linear_stationarity():
    num_vars = 2
    num_cons = 1

    H = np.eye(num_vars)
    x0 = np.array([0., 0.])
    g = -x0

    lb = np.full((2,), fill_value=-np.inf)
    ub = np.full((2,), fill_value=np.inf)

    A = np.array([[1., 1.]])
    b = np.array([1.])

    prob = QProblem(num_vars, num_cons)
    nWSR = np.array([10_000])

    prob.init(H, g, A, lb, ub, b, b, nWSR)

    x_opt = np.empty_like(x0)
    z_opt = np.empty((num_vars + num_cons,))

    prob.getPrimalSolution(x_opt)
    prob.getDualSolution(z_opt)

    var_dual = z_opt[:num_vars]
    cons_dual = z_opt[num_vars:]

    assert np.allclose(var_dual, 0.)

    assert np.allclose(H @ x_opt + g - cons_dual.T @ A,
                       0.)

    cons = prob.getConstraints()

    assert cons.getStatus(0) in [SubjectToStatus.LOWER,
                                 SubjectToStatus.UPPER]


def test_linear():
    num_vars = 2
    num_cons = 2

    H = np.eye(num_vars)
    x0 = np.array([1., 3.])
    g = -x0

    lb = np.zeros((2,))
    ub = np.full((2,), fill_value=np.inf)

    A = np.array([[-1., 1.],
                  [1., 1.]])

    ub_b = np.array([0.,
                     2.])

    lb_b = np.full((2,), fill_value=-np.inf)

    prob = QProblem(num_vars, num_cons)
    nWSR = np.array([10_000])

    prob.init(H, g, A, lb, ub, lb_b, ub_b, nWSR)

    x_opt = np.empty_like(x0)
    z_opt = np.empty((num_vars + num_cons,))

    prob.getPrimalSolution(x_opt)
    prob.getDualSolution(z_opt)

    assert np.allclose(x_opt, [1., 1.])

    y_var_opt = z_opt[:num_vars]
    y_cons_opt = z_opt[num_vars:]

    assert np.allclose(y_var_opt, 0.)
    assert np.allclose(y_cons_opt, [-1., -1.])
