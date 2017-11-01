#!/usr/bin/env python3

import numpy as np
linalg = np.linalg


def independent_columns(A, tol=1e-05):
    Q, R = linalg.qr(A)
    independent = np.where(np.abs(R.diagonal()) > tol)[0]
    return A[:, independent]

if __name__ == '__main__':
    f = open('../data_for_identification/ee/EE1.txt', 'r')
    ee = np.loadtxt(f)
    print(ee.shape)
    M = ee
    U, s, V = np.linalg.svd(M)
    print(s)
    print(np.sum(np.abs(s) > 1e-05))

    import sympy

    M = sympy.Matrix(M)
    reduced_form, inds = M.rref()
    print(np.array(reduced_form).shape)
    print(inds)
