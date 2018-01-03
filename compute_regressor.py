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

    print(colored('Start computing of regressor {:}'.format(time.ctime()), 'magenta'))

    zeros_in_regressor = np.ones((n, n * nL))
    for i in range(n):
        for j in range(i * nL):
            zeros_in_regressor[i, j] = 0

    for j in range(n):
        for i in range(j, n):
            start_time = time.time()
            regressorElement = RegressorElement(j, i)
            for k in range(nL):
                st = time.time()

                # compute operator L of lagrange function for [i, k, j]
                expr_raw = operatorL(L[i][k], j)
                expr_raw = expr_raw[0] if expr_raw.is_Matrix else expr_raw
                len_expr_raw = len(str(expr_raw))

                expr_raw.subs([delta[0], delta[1]], [DELTA[0], DELTA[1]])

                # simplify expression
                expr = (expr_raw)
                # opL_sym = combsimp(powsimp(trigsimp(expand(expr_raw))))   # alternative method
                len_expr = len(str(expr))

                # make record about zeros elements (for removing zeros columns)
                if expr.is_zero:
                    zeros_in_regressor[j][(i * nL) + k] = 0

                # if simplify was shit ;)
                expr = expr_raw if len_expr_raw < len_expr else expr

                # # подмена функций Theta(q_i) на theta_i
                # for l in range(n):
                #     opL_sym = opL_sym.subs(thi[l], theta[l])

                # generate python code
                py_expr_raw = sympy.printing.lambdarepr.lambdarepr(expr)

                # some replaces, e.g. a_1 to a[0], Derivative(q_1(t), t) to dq[0]
                py_expr = python_gencode(py_expr_raw)  # see in lins/regexps.py

                regressorElement.addOpL(k, py_expr)

                et = timedelta(seconds=time.time() - st)
                print('--ji={0}{1}, k={2} function computed for {3:>15}\n\tSimplify length of expression from {4} to {5}'.format(j, i, k, str(et), len_expr_raw,len_expr))
            writeFile(PATH + 'xi_{0}{1}.py'.format(j, i), str(regressorElement))
            end_time = timedelta(seconds=time.time() - start_time)
            print(colored('-Regressor element {0}{1} computed for {2:>15}'.format(j, i, str(end_time)), 'green'), '\n***')
    print(colored('End computing of regressor {:}'.format(time.ctime()), 'magenta'))

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

