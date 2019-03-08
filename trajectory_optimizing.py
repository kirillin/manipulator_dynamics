#!/usr/bin/env python3
import matplotlib.pyplot as plt

import numpy as np
from numpy import *

from libs.identification import Identification
from libs.initialization import *
from libs.utils import getFileNamesContains

from scipy.optimize import minimize

delta = np.array(DELTA)
ident = Identification(A, D, delta)


nf = 5  # qty harmonics
T = 2
w0 = 2 * pi / T  # base frequency
dt = 0.01


def condBigXi(theta, t_max=2):
    qs, dqs, ddqs = [], [], []
    for t in arange(0, t_max, dt):
        omega = np.matrix([
            [cos(k * w0 * t) if k % 2 == 0 else sin(k * w0 * t) for k in range(2 * nf)],
            [-sin(k * w0 * t) if k % 2 == 0 else cos(k * w0 * t) for k in range(2 * nf)],
            [-cos(k * w0 * t) if k % 2 == 0 else -sin(k * w0 * t) for k in range(2 * nf)],
        ])
        q = omega * theta
        qs.append(q[0].tolist()[0])
        dqs.append(q[1].tolist()[0])
        ddqs.append(q[2].tolist()[0])
    bigXi = ident.getBigXi(qs, dqs, ddqs)
    return np.linalg.cond(bigXi)

theta0 = np.matrix([np.random.rand(n) for i in range(2 * nf)])

func = lambda theta: condBigXi(theta)
theta_opt, fopt, _, _, _ = scipy.optimize.fmin(func=func, x0=theta_0, maxiter=maxiter, full_output=True)

print(theta_opt)
print(fopt)