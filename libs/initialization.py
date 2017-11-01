# coding: utf-8
#!/usr/bin/env python3

#
# Module contains initialization everything
#

from sympy import Matrix, symbols, Symbol, Function, diff, pi

# Инициализаця констант
n = 5       # quantity of links of Youbot's arm
nL = 10     # quantity of blablabla functions in Lagrange function
g = 9.82


# Объявление различных символьных переменных
t = Symbol('t')
a = [Symbol('a_' + str(i)) for i in range(1, n + 1)]
d = symbols("d_1 d_2 d_3 d_4 d_5")
m = symbols('m_1 m_2 m_3 m_4 m_5')
xc = symbols('xc0 xc1 xc2 xc3 xc4 xc5')
yc = symbols('yc0 yc1 yc2 yc3 yc4 yc5')
zc = symbols('zc0 zc1 zc2 zc3 zc4 zc5')
x = symbols('x0 x1 x2 x3 x4 x5')
y = symbols('y0 y1 y2 y3 y4 y5')
z = symbols('z0 z1 z2 z3 z4 z6')
Ixx = symbols('I1xx I2xx I3xx I4xx I5xx')
Iyy = symbols('I1yy I2yy I3yy I4yy I5yy')
Izz = symbols('I1zz I2zz I3zz I4zz I5zz')
Ixy = symbols('I1xy I2xy I3xy I4xy I5xy')
Ixz = symbols('I1xz I2xz I3xz I4xz I5xz')
Iyz = symbols('I1yz I2yz I3yz I4yz I5yz')
theta = [Symbol('theta_' + str(i)) for i in range(5)]

# Объявление q_i(t), dq_i(t), ddq_i(t) как функций от времени, где i = [0..n)
q, dq, ddq = [list() for i in range(3)]
for i in range(0, n):
    q.append(Function('q_' + str(i + 1))(t))
    dq.append(diff(q[i], t))
    ddq.append(diff(q[i], t, 2))


class Theta(Function):
    def fdiff(self, argindex=1):
        if len(self.args) > 0:
            return -1
        else:
            return 0

# DH-parameters
ai = [a[0], a[1], a[2], 0, 0]
di = [d[0], 0, 0, 0, d[4]]
alphai = [pi / 2, 0, 0, pi / 2, 0]
thi = [q[0], q[1], q[2], q[3], q[4]]
thi = [Theta(q[i]) for i in range(5)]

# Ниже настоящие параметры робота, в силу непреодолимости некоторых обстоятельств,
#  они подменяют предыдущие после этапа упрощения уравений динамики манипулятора
# thi = [pi * (169. / 180.) - q[0],
#         pi * (65. / 180.) + (pi / 2.) - q[1],
#         -pi * (146. / 180.) - q[2],
#         pi * (102.5 / 180) + (pi / 2.) - q[3],
#         pi * (167.5 / 180.) - q[4]]

# Инициализация простых, но нужных вещей...
g0 = Matrix([[0], [0], [-g]])   # вектор графитации
O = Matrix([0, 0, 0])           # используется при инициализации Якобианов


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



if __name__ == '__main__':
    print(O)