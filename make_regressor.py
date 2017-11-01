# coding: utf-8
#!/usr/bin/env python3

#
# Module makes regressor XI of
#   tau = XI * CHI
#   tau -- measured from real system
#   CHI -- we want estimate and we do!!1!
#


import datetime as dt
import sympy

from libs.equations import *
from libs.regexps import *

from numpy import pi


def makeXImodules():
    """
        1) Compute 150 components of regressor;
        2) Set it to modules.
        An example of the generating module see in libs/example_module.py

        Function create two directories that contain files with expressions for
        Lagrange-equations without parameters (mass, inertia, friction etc.)
        - <<xi>> -- expressions in sympy form;
        - <<xi_tex>> -- expressions in latex form.
    """
    regressor_zeros = [[1 for j in range(n * nL)] for i in range(5)]
    print('\nRegressor computing start ' + time.ctime())
    for j in range(4, n):
        print('\tRow ' + str(j) + ' start at ' + time.ctime())
        for i in range(j, n):
            print('\t\tElement ' + str(i) + ' start at ' + time.ctime())

            fileName = 'xi/xi_' + str(j) + str(i) + '.py'

            file = open(fileName, 'w')

            # First piece of template of module
            module = '#!/usr/bin/env python3\n'
            module += 'from numpy import cos, sin, pi, sqrt\n'
            module += '\n\n'
            module += 'class XI:\n'
            module += '\t"""XI_' + str(i) + str(j) + '"""\n\n'
            module += '\tdef __init__(self, q=(0,0,0,0,0), dq=(0,0,0,0,0), ' \
                      'ddq=(0,0,0,0,0), a=(0,0,0,0,0), d=(0,0,0,0,0), theta=(0,0,0,0,0)):\n' \
                      '\t\tself.q, self.dq, self.ddq = q, dq, ddq\n' \
                      '\t\tself.a, self.d = a, d\n' \
                      '\t\tself.theta = theta\n\n'
            module += '\tdef setData(self, q, dq, ddq, a, d):\n' \
                      '\t\tself.q, self.dq, self.ddq = q, dq, ddq\n' \
                      '\t\tself.a, self.d = a, d\n\n'
            file.write(module)
            file.close()

            for k in range(0, nL):
                print('\t\t\tSub Element ' + str(k) + ' start at ' + time.ctime())

                file = open(fileName, 'a')

                # compute expression
                opL_sym_raw = operatorL(L[i][k], j)

                # it is work and okay!!1!
                print(type(opL_sym_raw))
                if type(opL_sym_raw) in [ImmutableMatrix]:
                    opL_sym_raw = opL_sym_raw[0]
                len_raw = len(str(opL_sym_raw))
                print(len_raw)

                # simplify expression
                opL_sym = combsimp(powsimp(trigsimp(expand(opL_sym_raw))))
                # opL_sym = simplify(opL_sym_raw)   # alternative method
                # opL_sym = opL_sym_raw             # raw expressions

                # make record about zeros elements (for removing zeros columns)
                if opL_sym == sympy.numbers.Zero():
                    regressor_zeros[j][(i * 10) + k] = 0

                len_ok = len(str(opL_sym))
                print(len_ok)
                time.sleep(0.42)

                # if simplify was shit ;)
                if len_raw < len_ok:
                    opL_sym = opL_sym_raw

                " START CREATING method in MODULE "

                # подмена функций Theta(q_i) на theta_i
                for l in range(5):
                    opL_sym = opL_sym.subs(thi[l], theta[l])

                # generate python code
                opL_py_raw = sympy.printing.lambdarepr.lambdarepr(opL_sym)

                # some replaces, e.g. a_1 to a[0], Derivative(q_1(t), t) to dq[0]
                opL_py_well = python_gencode(opL_py_raw)  # see in lins/regexps.py
                #opL_py_well = python_gencode(opL_py_raw, real_thi)    # see in lins/regexps.py

                opL_sym = opL_py_well
                opNum = str(k)
                # aaaaand template ooof module

                # second piece of template of module
                # methods contains computing operator L
                module = '\tdef opL' + opNum + '(self):\n' \
                            '\t\t"""' + str(i) + str(j) + str(k) + '"""\n' \
                            '\t\tq, dq, ddq = self.q, self.dq, self.ddq\n' \
                            '\t\ta, d = self.a, self.d\n' \
                            '\t\ttheta = self.theta\n' \
                            '\t\topL_' + opNum + ' = ' + opL_sym + '\n' \
                            '\t\treturn opL_' + opNum + '\n\n'

                " END CREATING method in MODULE "

                # writing data to files
                file.write(module)
                file.close()

                print('\t\t\tSub Element ' + str(k) + ' end at ' + time.ctime())
            file = open(fileName, 'a')
            # third piece of template of module
            module = '\tdef getXI(self, q, dq, ddq):\n' \
                      '\t\tself.q = q\n' \
                      '\t\tself.dq = dq\n' \
                      '\t\tself.ddq = ddq\n' \
                      '\t\tXI = [0 for i in range(10)]\n' \
                      '\t\tXI[0] = self.opL0()\n' \
                      '\t\tXI[1] = self.opL1()\n' \
                      '\t\tXI[2] = self.opL2()\n' \
                      '\t\tXI[3] = self.opL3()\n' \
                      '\t\tXI[4] = self.opL4()\n' \
                      '\t\tXI[5] = self.opL5()\n' \
                      '\t\tXI[6] = self.opL6()\n' \
                      '\t\tXI[7] = self.opL7()\n' \
                      '\t\tXI[8] = self.opL8()\n' \
                      '\t\tXI[9] = self.opL9()\n' \
                      '\t\treturn XI\n\n'
            file.write(module)
            file.close()
            print('\t\tElement ' + str(i) + ' end at ' + time.ctime())

        print('\tRow ' + str(j) + ' end at ' + time.ctime())
    print("Все!!")

    # Файл с нулевыми элементами, в котором можно узреть нулевые столбцы
    f = open('xi/regressor_zeros.txt', 'w')
    for i in range(1, n):
        for j in range(1, i * nL ):
            regressor_zeros[i][j] = 0
    for i in range(n):
        f.write(str(regressor_zeros[i])+'\n')
    f.close()


if __name__ == '__main__':
    makeXImodules()
