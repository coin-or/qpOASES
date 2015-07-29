"""
This file is part of qpOASES.

qpOASES -- An Implementation of the Online Active Set Strategy.
Copyright (C) 2007-2015 by Hans Joachim Ferreau et al. All rights reserved.

qpOASES is free software; you can redistribute it and/or
modify it under the terms of the GNU Lesser General Public
License as published by the Free Software Foundation; either
version 2.1 of the License, or (at your option) any later version.

qpOASES is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
See the GNU Lesser General Public License for more details.

You should have received a copy of the GNU Lesser General Public
License along with qpOASES; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
"""

#TODO add doxygen support
# \author Manuel Kudruss
# \version 3.1
# \date 2013-2015

import os
import numpy as np
from numpy.testing import *
from qpoases import PyQProblem as QProblem
from qpoases import PyBooleanType as BooleanType
from qpoases import PyOptions as Options
from qpoases import PyPrintLevel as PrintLevel

# get qpOASES path
qpoases_path = os.path.dirname(os.path.abspath(__file__))
qpoases_path = os.path.dirname(qpoases_path)
qpoases_path = os.path.dirname(qpoases_path)
qpoases_path = os.path.dirname(qpoases_path)

# set qpOASES testing path
testing_path = os.path.join(qpoases_path, "testing")

class TestIdHessian(TestCase):

    def test_id_hessian(self):
        """Very simple example for testing qpOASES (using QProblem class)."""

        path = os.path.join(testing_path, "dev_idhessian_data")

        #Setup data for QP.
        H   = np.loadtxt(os.path.join(path, "H.txt"))
        g   = np.loadtxt(os.path.join(path, "g.txt"))
        A   = np.loadtxt(os.path.join(path, "A.txt"))
        lb  = np.loadtxt(os.path.join(path, "lb.txt"))
        ub  = np.loadtxt(os.path.join(path, "ub.txt"))
        lbA = np.loadtxt(os.path.join(path, "lbA.txt"))
        ubA = np.loadtxt(os.path.join(path, "ubA.txt"))

        #Setting up QProblem object.
        qp = QProblem(72,144)

        options = Options()
        options.numRefinementSteps   = 1
        options.printLevel = PrintLevel.NONE

        #options.setToMPC()
        #options.setToReliable()
        #options.enableFlippingBounds = BooleanType.FALSE
        options.enableRamping = BooleanType.FALSE
        #options.enableRamping = BooleanType.TRUE
        #options.enableFarBounds = BooleanType.FALSE
        #options.enableRamping = BooleanType.FALSE
        #options.printLevel = PL_LOW
        #options.enableFullLITests = BooleanType.FALSE
        #options.boundRelaxation = 1.0e-1
        qp.setOptions( options )

        #Solve QP.
        nWSR = 1200

        qp.init(H, g, A, lb, ub, lbA, ubA, nWSR)

        # FIXME check against what?
        # Where can I find solution?

if __name__=="__main__":
    try:
        import nose
        nose.runmodule()

    except ImportError:
        sys.stderr.write("Please install nosestests for python unittesting.\n")

