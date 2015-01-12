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

author Manuel Kudruss
version 3.1
date 2013-2015
"""

import os
import numpy as np
from numpy.testing import *
from qpoases import py_runOQPbenchmark as runOQPbenchmark
from qpoases import PyQProblem as QProblem
from qpoases import PyBooleanType as BooleanType
from qpoases import PyReturnValue as ReturnValue
from qpoases import PyOptions as Options
from qpoases import PyPrintLevel as PrintLevel

# get qpOASES path
qpoases_path = os.path.dirname(os.path.abspath(__file__))
qpoases_path = os.path.dirname(qpoases_path)
qpoases_path = os.path.dirname(qpoases_path)
qpoases_path = os.path.dirname(qpoases_path)

# set qpOASES testing path
testing_path = os.path.join(qpoases_path, 'testing')


benchmarks = ('CVXQP1_S',
              'CVXQP2_S',
              'CVXQP3_S',
              'DPKLO1',
              'DUAL1',
              'DUAL2',
              'DUAL3',
              'DUAL4',
              'DUALC1',
              'DUALC2',
              'DUALC5',
              'DUALC8',
              'GENHS28',
              'HS118',
              'HS21',
              'HS268',
              'HS35',
              'HS35MOD',
              'HS51',
              'HS52',
              'HS53',
              'HS76',
              'LOTSCHD',
              'PRIMALC1',
              'PRIMALC2',
              'PRIMALC5',
              'QADLITTL',
              'QAFIRO',
              'QBEACONF',
              'QBRANDY',
              'QE226',
              'QISRAEL',
              'QPCBLEND',
              'QPCBOEI2',
              'QPTEST',
              'QRECIPE',
              'QSC205',
              'QSCAGR7',
              'QSHARE1B',
              'QSHARE2B',
              'S268',
              'TAME',
              'VALUES',
              'ZECEVIC2',
              )

def results2str(results):
    """converts results dictionary to pretty string"""
    hline = '{0:->11}|{0:-<12}|{0:-<12}|{0:-<12}|{0:-<8}|{0:-<8}\n'.format("")

    npass = 0
    nfail = 0

    string = ""
    string += '{:<10} | {: <10} | {: <10} | {: <10} | {: <6} | {:<6}\n'.format('problem', 'stat', 'feas', 'compl', 'nWSR', 'result')
    string += hline
    for key in results:
        line = '{name:<10} | {stat: >10.4e} | {feas: >10.4e} | {comp: >10.4e} | {nwsr: >6d} | {pass!s:<6}\n'.format(**results[key])
        string += line
        if results[key]['pass']:
            npass += 1
        else:
            nfail +=1

    string += hline
    string += '\n'
    string += 'Testbench results:\n'
    string += '==================\n'
    string += 'Pass:  {: >10d}\n'.format(npass)
    string += 'Fail:  {: >10d}\n'.format(nfail)
    string += '------------------\n'
    string += 'Ratio: {: >10.2%}\n'.format(npass/float(len(results)))

    return string

def write_results(name, string):
    """writes results into results dictionary"""
    path = os.path.dirname(os.path.abspath(__file__))

    directory = os.path.join(path, 'results')

    if not os.path.exists(directory):
        os.makedirs(directory)

    with open(os.path.join(directory, name), 'w') as f:
        f.write(string)

def get_nfail(results):
    """get nfail from results dictionary"""
    nfail = 0
    for key in results:
        if not results[key]['pass']:
            nfail +=1

    return nfail

def run_benchmarks(benchmarks, options, isSparse, useHotstarts,
                   nWSR, cpu_time, TOL):
    """run all benchmarks and return results as dictionary"""
    results = {}
    for item in benchmarks:
        # NOTE: c/c++ function epects trailing slash in path!
        path = os.path.join(testing_path, 'problems', item, '')

        # Run QP benchmark
        returnvalue, maxNWSR, avgNWSR, maxCPUtime, avgCPUtime, \
        maxStationarity, maxFeasibility, maxComplementarity \
        = runOQPbenchmark(path, isSparse, useHotstarts,
                          options, nWSR, cpu_time )

        if (returnvalue == ReturnValue.SUCCESSFUL_RETURN
            and maxStationarity    < TOL
            and maxFeasibility     < TOL
            and maxComplementarity < TOL):
            ret_val = True

        else:
            ret_val = False

        tmp_d = {'name': item,
                 'stat': maxStationarity,
                 'feas': maxFeasibility,
                 'comp': maxComplementarity,
                 'nwsr': int(maxNWSR),
                 'pass': bool(ret_val),
                }
        results[item] = tmp_d

    return results


class Testbench(TestCase):

    def setUp(self):
        # Setup global options for every problem
        self.TOL = 1e-5
        self.nWSR = 3.10
        self.cpu_time = 3.1
        self.decimal = 7 # number of decimals for assert

    def test_m44_default_dense(self):
        test_name = 'mm44_default_dense.txt'
        print("Test: ", test_name)

        # QP Options
        options = Options()
        options.setToDefault()
        options.printLevel = PrintLevel.NONE

        isSparse = False
        useHotstarts = False

        # run QP benchmarks
        results = run_benchmarks(benchmarks, options, isSparse, useHotstarts,
                                 self.nWSR, self.cpu_time, self.TOL)

        # print and write results
        string = results2str(results)
        print(string)
        write_results(test_name, string)

        assert get_nfail(results) <= 0, 'One ore more tests failed.'

    def test_m44_default_sparse(self):
        test_name = 'mm44_default_sparse.txt'
        print("Test: ", test_name)

        # QP Options
        options = Options()
        options.setToDefault()
        options.printLevel = PrintLevel.NONE

        isSparse = True
        useHotstarts = False

        # run QP benchmarks
        results = run_benchmarks(benchmarks, options, isSparse, useHotstarts,
                                 self.nWSR, self.cpu_time, self.TOL)

        # print and write results
        string = results2str(results)
        print(string)
        write_results(test_name, string)

        assert get_nfail(results) <= 0, 'One ore more tests failed.'

    def test_m44_mpc_dense(self):
        test_name ='mm44_mpc_dense.txt'
        print("Test: ", test_name)

        # QP Options
        options = Options()
        options.setToMPC()
        options.printLevel = PrintLevel.NONE

        isSparse = False
        useHotstarts = False

        # run QP benchmarks
        results = run_benchmarks(benchmarks, options, isSparse, useHotstarts,
                                 self.nWSR, self.cpu_time, self.TOL)

        # print and write results
        string = results2str(results)
        print(string)
        write_results(test_name, string)

        assert get_nfail(results) <= 2, 'One ore more tests failed.'

    def test_m44_mpc_sparse(self):
        test_name ='mm44_mpc_sparse.txt'
        print("Test: ", test_name)

        # QP Options
        options = Options()
        options.setToMPC()
        options.printLevel = PrintLevel.NONE

        isSparse = True
        useHotstarts = False

        # run QP benchmarks
        results = run_benchmarks(benchmarks, options, isSparse, useHotstarts,
                                 self.nWSR, self.cpu_time, self.TOL)

        # print and write results
        string = results2str(results)
        print(string)
        write_results(test_name, string)

        assert get_nfail(results) <= 19, 'One ore more tests failed.'

    def test_m44_reliable_dense(self):
        test_name = 'mm44_reliable_dense.txt'
        print("Test: ", test_name)

        # QP Options
        options = Options()
        options.setToReliable()
        options.printLevel = PrintLevel.NONE

        isSparse = False
        useHotstarts = False

        # run QP benchmarks
        results = run_benchmarks(benchmarks, options, isSparse, useHotstarts,
                                 self.nWSR, self.cpu_time, self.TOL)

        # print and write results
        string = results2str(results)
        print(string)
        write_results(test_name, string)

        assert get_nfail(results) <= 0, 'One ore more tests failed.'

    def test_m44_reliable_sparse(self):
        test_name = 'mm44_reliable_sparse.txt'
        print(test_name)

        # QP Options
        options = Options()
        options.setToReliable()
        options.printLevel = PrintLevel.NONE

        isSparse = True
        useHotstarts = False

        # run QP benchmarks
        results = run_benchmarks(benchmarks, options, isSparse, useHotstarts,
                                 self.nWSR, self.cpu_time, self.TOL)

        # print and write results
        string = results2str(results)
        print(string)
        write_results(test_name, string)

        assert get_nfail(results) <= 0, 'One ore more tests failed.'


if __name__=='__main__':
    try:
        import nose
        nose.runmodule(argv=['', '-s', '-v'])

    except ImportError:
        sys.stderr.write('Please install nosestests for python unittesting.\n')

