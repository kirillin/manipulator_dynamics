# coding: utf-8
#!/usr/bin/env python3
import os, errno
import time
import sympy
import numpy as np
from datetime import timedelta
from termcolor import colored

from libs.equations import *
from libs.regexps import *
from libs.regressor_stuff import *
from libs.utils import writeFile
from libs.initialization import *

from libs.Watchdog import Watchdog

SIMPIFY_TIMEOUT = 30 * 60 * 60 # [sec]
PATH = 'regressors/' + manipulator + '_xi/'

try:
    os.makedirs(PATH)
    open(PATH + '__init__.py', 'w+')
except OSError as e:
    if e.errno != errno.EEXIST:
        raise


def computeRegressor(path, zeros_in_regressor=[]):
    Xi = Regressor(path, zeros_in_regressor)
    writeFile(PATH + 'Xi.py', str(Xi))


def computeRegressorElements():
    regressorScilab = open('xi', 'w')
    print(colored('Start computing of regressor {:}'.format(time.ctime()), 'magenta'))

    zeros_in_regressor = np.ones((n, n * nL))
    for i in range(n-1, n):
        for j in range(i * nL):
            zeros_in_regressor[i, j] = 0

    for j in range(n):
        for i in range(j, n):
            start_time = time.time()
            regressorElement = RegressorElement(j, i)
            print(colored('-Regressor element {0}{1} computing start at {2:>15}'.format(j, i, time.ctime()), 'green'))
            for k in range(nL):
                st = time.time()
                print('--ji={0}{1}, k={2}'.format(j, i, k))

                # compute operator L of lagrange function for [i, k, j]
                expr_raw = operatorL(L[i][k], j)
                expr_raw = expr_raw[0] if expr_raw.is_Matrix else expr_raw
                len_expr_raw = len(str(expr_raw))
                print('\tL was calculated! ({0})'.format(len_expr_raw), end=' ', flush=True)

                # simplify expression
                expr = expr_raw
                try:
                    with Watchdog(SIMPIFY_TIMEOUT):
                        # expr = expand(expr)
                        # print(colored('expand!({0})'.format(len(str(expr))), 'yellow'), end=' ', flush=True)
                        #
                        # expr = factor(expr)
                        # print(colored('factor!({0})'.format(len(str(expr))), 'yellow'), end=' ', flush=True)

                        expr = trigsimp(expr)
                        print(colored('trigsimp!({0})'.format(len(str(expr))), 'yellow'), end=' ', flush=True)

                        expr = powsimp(expr)
                        print(colored('powsimp!({0})'.format(len(str(expr))), 'yellow'), end=' ', flush=True)

                        expr = combsimp(expr)
                        print(colored('combsimp!({0})'.format(len(str(expr))), 'yellow'), end=' ', flush=True)
                except Watchdog:
                    expr = expr_raw
                    print(colored(' Faild simplify ijk={0}{1}{2}'.format(i, j, k), 'red'), end=' ', flush=True)

                # expr = combsimp(powsimp(trigsimp(expand(expr_raw))))   # alternative method
                len_expr = len(str(expr))

                # make record about zeros elements (for removing zeros columns)
                if expr.is_zero:
                    zeros_in_regressor[j][(i * nL) + k] = 0

                # if simplify was shit ;)
                expr = expr_raw if len_expr_raw < len_expr else expr

                # # подмена функций Theta(q_i) на theta_i
                # for l in range(n):
                #     opL_sym = opL_sym.subs(thi[l], theta[l])

                expr = expr.subs(dict(zip(thi, theta)))

                # generate python code
                py_expr_raw = sympy.printing.lambdarepr.lambdarepr(expr)

                # some replaces, e.g. a_1 to a[0], Derivative(q_1(t), t) to dq[0]
                py_expr = python_gencode(py_expr_raw)  # see in lins/regexps.py

                regressorScilab.write('xi{0}{1}{2} = {3};\n'.format(j+1,i+1,k+1, py_expr))

                regressorElement.addOpL(k, py_expr)

                et = timedelta(seconds=time.time() - st)
                print('\n\tFunction computed for {3:>15}'.format(j, i, k, str(et)))
                print('\tSimplify length of expression from {0} to {1}'.format(len_expr_raw,len_expr))

            writeFile(PATH + 'xi_{0}{1}.py'.format(j, i), str(regressorElement))
            end_time = timedelta(seconds=time.time() - start_time)
            print(colored('-Regressor element {0}{1} computed for {2:>15}'.format(j, i, str(end_time)), 'green'), '\n***')
    print(colored('End computing of regressor {:}'.format(time.ctime()), 'magenta'))

    regressorScilab.close()

    # Файл с нулевыми элементами, в котором можно узреть нулевые столбцы
    np.savetxt(PATH + 'zeros_in_regressor.txt', zeros_in_regressor)


def getZeroCols(zeros_in_regressor):
    zero_cols = []
    for i in range(n * nL):
        if sum(zeros_in_regressor[:, i]) == 0:
            zero_cols.append(str(i))
    return zero_cols


if __name__ == '__main__':
    computeRegressorElements()


    file = open(PATH + 'zeros_in_regressor.txt', 'r')
    zeros_in_regressor = np.loadtxt(file)
    computeRegressor(PATH, getZeroCols(zeros_in_regressor))
