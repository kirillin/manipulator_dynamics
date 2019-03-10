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


# Setting element of regressor by user
import argparse

from_row, to_row = 0, n
from_col, to_col = 0, n
from_L, to_L = 0, nL

parser = argparse.ArgumentParser()

parser.add_argument("-r", type=int, default=-1, required=False, help="Set what is the row[1..5] to need compute!")
parser.add_argument("-c", type=int, default=-1, required=False, help="Set what is the column[r..5] to need compute!")
parser.add_argument("-l", type=int, default=-1, required=False, help="Set what is the lagrange element[1..10] to need compute!")

args = parser.parse_args()
row = int(args.r)
col = int(args.c)
l = int(args.l)

if row > 0:
    from_row, to_row = row-1, row
if col > 0:
    from_col, to_col = col-1, col
if l > 0:
    from_L, to_L = l-1, l
##

import mpmath
def nsimplify2(expr, constants=(), tolerance=None, full=False, rational=None, rational_conversion='base10'):
    from sympy.core.compatibility import (iterable, ordered, range, as_int)

    def _real_to_rational(expr, tolerance=None, rational_conversion='base10'):
        from sympy.core.sympify import _sympify
        expr = _sympify(expr)
        inf = Float('inf')
        p = expr
        reps = {}
        reduce_num = None
        if tolerance is not None and tolerance < 1:
            reduce_num = ceiling(1 / tolerance)
        for fl in p.atoms(Float):
            key = fl
            if reduce_num is not None:
                # why
                if fl not in [-1, 0, 1]:
                    r = Float(Rational(fl).limit_denominator(reduce_num), 10)
                else:
                    r = Rational(fl).limit_denominator(reduce_num)
            elif (tolerance is not None and tolerance >= 1 and
                  fl.is_Integer is False):
                r = Rational(tolerance * round(fl / tolerance)).limit_denominator(int(tolerance))

            else:
                if rational_conversion == 'exact':
                    r = Rational(fl)
                    reps[key] = r

                    continue
                elif rational_conversion != 'base10':
                    raise ValueError("rational_conversion must be 'base10' or 'exact'")

                r = nsimplify(fl, rational=False)
                # e.g. log(3).n() -> log(3) instead of a Rational
                if fl and not r:
                    r = Rational(fl)
                elif not r.is_Rational:
                    if fl == inf or fl == -inf:
                        r = S.ComplexInfinity
                    elif fl < 0:
                        fl = -fl
                        d = Pow(10, int((mpmath.log(fl) / mpmath.log(10))))
                        r = -Rational(str(fl / d)) * d
                    elif fl > 0:
                        d = Pow(10, int((mpmath.log(fl) / mpmath.log(10))))
                        r = Rational(str(fl / d)) * d
                    else:
                        r = Integer(0)

            reps[key] = r
            # print(key, r)
        # print(reps)
        return p.subs(reps, simultaneous=True)

    try:
        return sympify(as_int(expr))
    except (TypeError, ValueError):
        pass
    expr = sympify(expr).xreplace({
        Float('inf'): S.Infinity,
        Float('-inf'): S.NegativeInfinity,
        })
    if expr is S.Infinity or expr is S.NegativeInfinity:
        return expr
    if rational or expr.free_symbols:
        return _real_to_rational(expr, tolerance, rational_conversion)

    # SymPy's default tolerance for Rationals is 15; other numbers may have
    # lower tolerances set, so use them to pick the largest tolerance if None
    # was given
    if tolerance is None:
        tolerance = 10**-min([15] +
             [mpmath.libmp.libmpf.prec_to_dps(n._prec)
             for n in expr.atoms(Float)])
    # XXX should prec be set independent of tolerance or should it be computed
    # from tolerance?
    prec = 30
    bprec = int(prec*3.33)

    constants_dict = {}
    for constant in constants:
        constant = sympify(constant)
        v = constant.evalf(prec)
        if not v.is_Float:
            raise ValueError("constants must be real-valued")
        constants_dict[str(constant)] = v._to_mpmath(bprec)

    exprval = expr.evalf(prec, chop=True)
    re, im = exprval.as_real_imag()

    # safety check to make sure that this evaluated to a number
    if not (re.is_Number and im.is_Number):
        return expr

    def nsimplify_real(x):
        orig = mpmath.mp.dps
        xv = x._to_mpmath(bprec)
        try:
            # We'll be happy with low precision if a simple fraction
            if not (tolerance or full):
                mpmath.mp.dps = 15
                rat = mpmath.pslq([xv, 1])
                if rat is not None:
                    return -rat[1]/ rat[0] #Rational(-int(rat[1]), int(rat[0]))
            mpmath.mp.dps = prec
            newexpr = mpmath.identify(xv, constants=constants_dict,
                tol=tolerance, full=full)
            if not newexpr:
                raise ValueError
            if full:
                newexpr = newexpr[0]
            expr = sympify(newexpr)
            if x and not expr:  # don't let x become 0
                raise ValueError
            if expr.is_finite is False and not xv in [mpmath.inf, mpmath.ninf]:
                raise ValueError
            return expr
        finally:
            # even though there are returns above, this is executed
            # before leaving
            mpmath.mp.dps = orig
    try:
        if re:
            re = nsimplify_real(re)
        if im:
            im = nsimplify_real(im)
    except ValueError:
        if rational is None:
            return _real_to_rational(expr, rational_conversion=rational_conversion)
        return expr

    rv = re + im*S.ImaginaryUnit
    # if there was a change or rational is explicitly not wanted
    # return the value, else return the Rational representation
    if rv != expr or rational is False:
        return rv

    return _real_to_rational(expr, rational_conversion=rational_conversion)


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

    for j in range(from_row, to_row):

        if from_col != 0:
            j = from_col

        for i in range(j, to_col):
            start_time = time.time()
            regressorElement = RegressorElement(j, i)
            print(colored('-Regressor element {0}{1} computing start at {2:>15}'.format(j, i, time.ctime()), 'green'))
            for k in range(from_L, to_L):
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
                        expr = nsimplify2(expr, tolerance=1e-14, rational=False)
                        #expr = nfloat(expr, 14)
                        print(colored('nsimp&nfloat!({0})'.format(len(str(expr))), 'yellow'), end=' ', flush=True)

                        expr = expand(expr)
                        print(colored('expand!({0})'.format(len(str(expr))), 'yellow'), end=' ', flush=True)

                        expr = factor(expr)
                        print(colored('factor!({0})'.format(len(str(expr))), 'yellow'), end=' ', flush=True)

                        expr = trigsimp(expr)
                        print(colored('trigsimp!({0})'.format(len(str(expr))), 'yellow'), end=' ', flush=True)

                        expr = powsimp(expr)
                        print(colored('powsimp!({0})'.format(len(str(expr))), 'yellow'), end=' ', flush=True)

                        expr = combsimp(expr)
                        print(colored('combsimp!({0})'.format(len(str(expr))), 'yellow'), end=' ', flush=True)

                        "!!! one more time"
                        expr = nsimplify2(expr, tolerance=1e-14, rational=False)
                        if len(str(expr)) < 10000:
                            expr = simplify(expr)
                            print(colored('simplify!({0})'.format(len(str(expr))), 'yellow'), end=' ', flush=True)
                        expr = nsimplify2(expr, tolerance=1e-14, rational=False)
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

                "Write to file expr for next work with it"
                str_expr_file = open('str_expr/{}{}{}.txt'.format(j, i, k), 'w')
                str_expr_file.write(str(expr))
                str_expr_file.close()

                str_expr_file = open('str_expr/_{}{}{}.txt'.format(j, i, k), 'w')
                str_expr_file.write(str(expr_raw))
                str_expr_file.close()

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
