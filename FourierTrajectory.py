""" Setup:
    1) spectrum of signal
    2) period of signal
    3) maximum allowed limits
    4)* duration of trajectory

"""
from numpy import *
import matplotlib.pyplot as plt

S0 = [0 for i in range(100)]
S0[1] = 0
S0[10] = 0
S0[20] = 0.001
T0 = 2.5
DURATION = 5


class FourierTrajectory:

    def __init__(self, T, s, allowed_limits=(0, 0)):
        self.s = s  # spectrum
        self.T = T  # period of signal
        self.w0 = 2 * pi / T
        self.allowed_limits = allowed_limits

        self.a, self.b = [-1], [-1]

        self.q0 = 0
        self.scale = 1

        self.__compute_ab()
        self.__compute_adjust()

    def __del__(self):
        pass

    def get_a(self):
        return self.a

    def get_b(self):
        return self.b

    def __compute_ab(self):
        func = fft.ifft(self.s)
        complex_s = fft.fft(func)
        for s in complex_s:
            psi = arctan2(imag(s), real(s))
            self.a.append(2 * abs(s) * cos(psi))
            self.b.append(-2 * abs(s) * sin(psi))

    def __compute_adjust(self):
        if not self.allowed_limits[0] == 0 and not self.allowed_limits[1] == 0:
            __, qs, __, __ = self.get_trajectory(self.T)
            A = amax(qs) - amin(qs)
            A0 = abs(self.allowed_limits[1] - self.allowed_limits[0])
            # TODO: check computation of scale and shift (doesn't work very well)
            self.scale = amin([A, A0]) / amax([A, A0])
            self.q0 = self.allowed_limits[1] - self.scale * amax(qs)

    def get_point(self, t, q0=0, scale=1):
        q, dq, ddq = 0, 0, 0
        for k in range(1, len(self.s) + 1):
            # TODO: or not TODO? Add phase compute. phi = atan2(b[k] / a[k]) -- to make velocity equals to zero at t=0
            q = q + self.a[k] * cos(k * self.w0 * t) + self.b[k] * sin(k * self.w0 * t)
            dq = dq + self.a[k] * k * self.w0 * sin(k * self.w0 * t) \
                    - self.b[k] * k * self.w0 * cos(k * self.w0 * t)
            ddq = ddq - self.a[k] * k**2 * self.w0**2 * cos(k * self.w0 * t) \
                    - self.b[k] * k**2 * self.w0**2 * sin(k * self.w0 * t)
        q = q0 + scale * q
        dq = scale * dq
        ddq = scale * ddq
        return q, dq, ddq

    def get_trajectory(self, duration, dt=0.001, q0=0, scale=1):
        ts, qs, dqs, ddqs = [], [], [], []
        for t in arange(0, duration, dt):
            q, dq, ddq = self.get_point(t, self.q0, self.scale)
            ts.append(t)
            qs.append(q)
            dqs.append(dq)
            ddqs.append(ddq)
        return ts, qs, dqs, ddqs


if __name__ == '__main__':
    ft = FourierTrajectory(T0, S0)

    ts, qs, dqs, ddqs = ft.get_trajectory(DURATION)

    print(amax(qs), amin(qs))
    fig1, (ax1, ax2) = plt.subplots(nrows=2, ncols=1)
    ax1.bar(arange(0, len(S0)), S0, 1)
    ax2.plot(ts, qs, ts, dqs)
    plt.show()
