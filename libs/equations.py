# coding: utf-8
#!/usr/bin/env python3

#
# Module contains all equation functions using in computations
#
# В функциях, если передать simp=True, то результат будет
#   упрощен в соответствии с содержанием упрощающей функции mySimple()
#
# В этом файле вызувыются функции (внизу!)
#   init_jacobians()
#   compute_lagrange_function()
#


import time
from sympy import *
from libs.initialization import *


def mySimple(expr):
    """ Тут задается какими алгоритмами упрощать """
    expr = simplify(expr)
    return expr


def T(i, j, simp=False):
    """
        Возвращает матрицу трансформации из i в j
        T(0,1) is equal to {}^0 A_1
    """
    if Ts[i][j] != -1:
        return Ts[i][j]
    else:
        H = eye(4)
        for k in range(i, j):
            Rz = Matrix(
                [[cos(thi[k]), -sin(thi[k]), 0, 0],
                 [sin(thi[k]), cos(thi[k]), 0, 0],
                 [0, 0, 1, 0], [0, 0, 0, 1]])
            Rx = Matrix([[1, 0, 0, 0],
                         [0, cos(alphai[k]), -sin(alphai[k]), 0],
                         [0, sin(alphai[k]), cos(alphai[k]), 0],
                         [0, 0, 0, 1]])
            Tz = Matrix([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, di[k]], [0, 0, 0, 1]])
            Tx = Matrix([[1, 0, 0, ai[k]], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])
            dT = Rz * Tz * Tx * Rx
            H = H * dT
        if simp:
            H = mySimple(H)
        Ts[i][j] = H
        return H


def R(i, j, simp=False):
    """ Rotation matrices, i > j """
    if Rs[i][j] != -1:
        return Rs[i][j]
    else:
        r = T(i, j)[:3, :3]
        if simp:
            r = mySimple(r)
        Rs[i][j] = r
        return r


def get_z(i, simp=True):
    """
        get_z(1) = z^0_1
    """
    if Zs[i] != -1:
        return Zs[i]
    else:
        z = R(0, i) * Matrix([[0], [0], [1]])
        if simp:
            z = mySimple(z)
        Zs[i] = z
        return z


def get_r0_0To(i):
    """
        get_r0_0To(1) = r^0_{0,1}
    """
    r = T(0, i)[:3, 3]
    return r  # 3x1


def get_ri_0To(i, simp=False):
    """
        get_ri_0To(1) = r^1_{0, 1}
    """
    r = transpose(R(0, i, simp)) * get_r0_0To(i)
    if simp:
        r = mySimple(r)
    return r  # 3x1


def get_g(i, simp=False):
    """
        get_g(1) = g_1
    """
    g = transpose(R(0, i, simp)) * g0
    if simp:
        g = mySimple(g)
    return g  # 3x1


def getJv(i):
    """
        getJv(0) = J_{v1}
    """
    j = Matrix(Jv[i])
    j = j.reshape(5, 3)
    return transpose(j)


def getJomega(i):
    """
        getJomega(0) = J_{omega1}
    """
    j = Matrix(Jomega[i])
    j = j.reshape(5, 3)
    return transpose(j)


def get_v0(i, simp=False):
    """
        get_v0(1) = v^0_1
    """
    v = - getJv(i-1) * Matrix(dq)
    return v


def get_vi(i, simp=False):
    """
        get_vi(1) = v^1_1
    """
    v = transpose(R(0, i, simp)) * get_v0(i, simp)
    if simp:
        v = mySimple(v)
    return v  # 3x1


def get_omega0(i, simp=False):
    """
        get_omega0(1) = omega^0_1
    """
    omega = - getJomega(i-1) * Matrix(dq)
    return omega


def get_omegai(i, simp=False):
    """
        get_omegai(1) = omega^1_1
    """
    omega = transpose(R(0, i, simp)) * get_omega0(i)
    if simp:
        omega = mySimple(omega)
    return omega  # 3x1


def init_jacobians(simp=False):
    """ MUST run first """
    # for v
    for i in range(0, n):
        for j in range(0, i + 1):
            el = get_z(j, simp).cross(get_r0_0To(i + 1) - get_r0_0To(j))
            Jv[i][j] = el
    # for omega
    for i in range(0, n):
        for j in range(0, i + 1):
            el = get_z(j, simp)
            Jomega[i][j] = el


def compute_lagrange_function(simp=False):
    """ Расчет функции Лагранжа"""
    tm = time.time
    totalT = tm()
    print('Start computig of Lagrange function')
    for i in range(0, n):
        """cols(L) = 0..9, rows(L) = 0..4"""
        startTrow = tm()
        l24 = get_vi(i + 1, simp).cross(get_omegai(i + 1, simp)) + get_g(i + 1, simp)
        omegai = get_omegai(i + 1, simp)
        L[i][0] = Rational(1, 2) * transpose(get_vi(i + 1, simp)) * get_vi(i + 1, simp) + transpose(get_g(i + 1, simp)) * get_ri_0To(
            i + 1, simp)
        # L[i][0] = Rational(1,2) * transpose(get_v0(i)) * get_v0(i) + transpose(g0) * get_r0_0To(i)
        L[i][1] = l24[0]
        L[i][2] = l24[1]
        L[i][3] = l24[2]
        L[i][4] = Rational(1, 2) * omegai[0] ** 2
        L[i][5] = Rational(1, 2) * omegai[1] ** 2
        L[i][6] = Rational(1, 2) * omegai[2] ** 2
        L[i][7] = omegai[0] * omegai[1]
        L[i][8] = omegai[0] * omegai[2]
        L[i][9] = omegai[1] * omegai[2]
        if simp:
            for j in range(0, nL):
                startTel = tm()
                L[i][j] = mySimple(L[i][j])
                print(str(tm() - startTel) + ': j= ' + str(j) + 'was simplified!')
        print(str(tm() - startTrow) + ': For link ' + str(i) + ' Lagrange function was computed!')
    totalT = tm() - totalT
    print('End computing! Total time: ' + str(totalT))


def operatorL(L, j, simp=False):
    """
        L: лагранжиан (L = K - U)
        j: порядковый номер обобщенной координаты по которой
            дифференцируем в уравнении Лагранжа
        :returns уравнение Лагранжа для j звена
        operatorL(L, 0) = ~~~~q_1
    """
    dL_Ddq = diff(L, dq[j])
    dLdq_Dt = diff(dL_Ddq, t)
    dL_Dq = diff(L, q[j])
    opL = dLdq_Dt - dL_Dq
    if simp:
        opL = mySimple(opL)
    return opL

init_jacobians()
compute_lagrange_function()

if __name__ == '__main__':
    print(L[1][3])
