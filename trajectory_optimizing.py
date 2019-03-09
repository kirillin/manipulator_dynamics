#!/usr/bin/env python3
import matplotlib.pyplot as plt

import numpy as np
from numpy import *
import scipy

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
    theta = reshape(theta, (2 * nf, 5))
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
    c = np.linalg.cond(bigXi)
    print(c)
    return c


theta_0 = np.matrix([np.random.rand(n) for i in range(2 * nf)])

theta_0 = [
    [0.61560672, 1.87854406, 0.10757574, 0.15288912, 0.04464681],
    [0.6857373, 0.02439988, 0.12402139, 0.1266051, 0.13560599],
    [0.09568189, 0.40952968, 0.18764632, 0.2612495, 0.22423524],
    [0.8057715, 0.45198997, 0.11324338, 0.23041487, 0.00880013],
    [1.06401582, 0.30938717, 0.23648248, 0.98165427, 0.18576885],
    [0.01919479, 0.61878327, 2.34644975, 0.12344445, 0.59231682],
    [0.52752144, 0.6752909, 0.01712927, 0.5210106, 0.01393499],
    [0.85676443, 0.01346193, 1.65133943, 0.64816887, 0.0673611],
    [0.398081, 0.35627036, 0.20157632, 0.79555341, 0.43982469],
    [0.11293154, 0.46448769, 0.01545221, 0.34791302, 1.03775898]
]

func = lambda theta: condBigXi(theta)
theta_opt, fopt, _, _, _ = scipy.optimize.fmin(func=func, x0=theta_0, maxiter=1000, full_output=True)

print(theta_opt)
print(fopt)
