# coding: utf-8
#!/usr/bin/env python3
"""
    Module contains initialization everything
"""

from sympy import Matrix, Symbol, Function, diff
from numpy import pi

# !!! check manipulator here, please
manipulator = {
    0: 'PLANAR_2DOF',
    1: 'KUKA_YOUBOT'
}[1]

#########################################
##### DON'T TOUCH ANITHING BELOW ########
#########################################

# init constants
nL = 10     # quantity of blablabla functions in Lagrange function
g = 9.82

# manipulator parameters
if manipulator is 'PLANAR_2DOF':
    n = 2  # quantity of links
    G0 = [0, -g, 0]
    # DH-params
    A = (0.5, 0.4)
    D = (0., 0.)
    ALPHA = (0, pi / 2)
    DELTA = (60 * pi / 180, 40 * pi / 180)
elif manipulator is 'KUKA_YOUBOT':
    n = 5  # quantity of links
    G0 = [0, 0, -g]
    # DH-params
    A = (0.033, 0.155, 0.135, 0., 0.)
    D = (0.147, 0, 0, 0, 0.218)
    ALPHA = (pi / 2, 0, 0, pi / 2, 0)
    DELTA = (pi * (169. / 180.),
             pi * (65. / 180.) + (pi / 2.),
             -pi * (146. / 180.),
             pi * (102.5 / 180) + (pi / 2.),
             pi * (167.5 / 180.))

# Объявление различных символьных переменных
t = Symbol('t')
a = [Symbol('a_' + str(i)) for i in range(1, n + 1)]
d = [Symbol('d_' + str(i)) for i in range(1, n + 1)]
m = [Symbol('m_' + str(i)) for i in range(1, n + 1)]
xc = [Symbol('xc' + str(i)) for i in range(0, n + 1)]
yc = [Symbol('yc' + str(i)) for i in range(0, n + 1)]
zc = [Symbol('zc' + str(i)) for i in range(0, n + 1)]
x = [Symbol('x' + str(i)) for i in range(0, n + 1)]
y = [Symbol('y' + str(i)) for i in range(0, n + 1)]
z = [Symbol('z' + str(i)) for i in range(0, n + 1)]
Ixx = [Symbol('I{}xx'.format(i)) for i in range(1, n + 1)]
Iyy = [Symbol('I{}yy'.format(i)) for i in range(1, n + 1)]
Izz = [Symbol('I{}zz'.format(i)) for i in range(1, n + 1)]
Ixy = [Symbol('I{}xy'.format(i)) for i in range(1, n + 1)]
Ixz = [Symbol('I{}xz'.format(i)) for i in range(1, n + 1)]
Iyz = [Symbol('I{}yz'.format(i)) for i in range(1, n + 1)]
theta = [Symbol('theta_' + str(i)) for i in range(1, n+1)]
delta = [Symbol('delta_' + str(i)) for i in range(1, n + 1)]

# Объявление q_i(t), dq_i(t), ddq_i(t) как функций от времени, где i = [0..n)
q, dq, ddq = [list() for i in range(3)]
for i in range(0, n):
    q.append(Function('q_' + str(i + 1))(t))
    dq.append(diff(q[i], t))
    ddq.append(diff(q[i], t, 2))

# DH-parameters 2 dof
ai = [a[i] if A[i] else 0 for i in range(len(A))]
di = [d[i] if D[i] else 0 for i in range(len(D))]
alphai = [ALPHA[i] if ALPHA[i] else 0 for i in range(len(ALPHA))]
# thi = [delta[0] - q[0], delta[1] - q[1]]


""" !!! use this block if don't work simplifying """
class Theta(Function):
    """ Theta(q_i(t)) = delta_i - q_i(t) """
    def fdiff(self, argindex=1):
        if len(self.args) > 0:
            return -1
        else:
            return 0

thi = [Theta(q[i]) for i in range(n)]

# Инициализация простых, но нужных вещей...
O = Matrix([0, 0, 0])           # используется при инициализации Якобианов
g0 = Matrix(G0)  # вектор гравитации 2 dof

# Массивы БОЛЬШИХ матриц для более эффективного расчета рекурсивных функций
Ts = [[-1 for j in range(n + 1)] for i in range(n + 1)]     # transformations
Rs = [[-1 for j in range(n + 1)] for i in range(n + 1)]     # rotations
Zs = [-1 for j in range(n + 1)]                             # z-axis
TENSORs = [-1 for j in range(n + 1)]
# jacobians
Jv = [[O for j in range(n)] for i in range(n)]
Jomega = [[O for j in range(n)] for i in range(n)]
# Lagrange function
L = [[-1 for i in range(0, nL)] for j in range(0, n)]
