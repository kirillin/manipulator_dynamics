#!/usr/bin/env python3
"""
    Collect xi matrix (regressor).
"""
import numpy as np

from regressors.KUKA_YOUBOT_xi.Xi import Xi
# from regressors.PLANAR_2DOF_xi.Xi import Xi
from libs.initialization import nL, A, D, DELTA


class XiNum(Xi):

    def __init__(self, a, d, delta):
        super().__init__(a, d, delta)
        self._linerCols = []
        #self._linerCols = [3,4,5,7,8,9,16,18,20,21,22,23,0,1,6]   # good 2dof
        self._linerCols = [2,3,4,7,8,0,1,5,20,28,10,14,15,18,22,32,42,45,60]    # Youbot
        # self._linerCols += [6,0,14]
        # self._linerCols = [3,4,5,7,8,9,16,18,20,21,22,23,6,0,14] #6 10

    def getLinerCols(self):
        return self._linerCols

    def getXiNum(self, q, dq, ddq):
        """ :returns Standart not extended regressor 5x50 in numerical form """
        n = len(q)
        xi = np.empty(0)
        for i in range(n):
            for j in range(n):
                if self._xi[i][j] != 0:
                    xi_ij = self._xi[i][j].getXi(q, dq, ddq)
                    xi = np.concatenate((xi, xi_ij))
                else:
                    xi = np.concatenate((xi, np.zeros(nL)))
        return xi.reshape(n, n * nL)

    def getXiNumEx(self, q, dq, ddq):
        """ :returns Extended regressor 5x5*14 in numerical form """
        n = len(q)
        xi = self.getXiNum(q, dq, ddq)
        basePart = xi.reshape(n * n, nL)
        extendPart = np.empty(0)  # 1 x 100
        for i in range(n):
            for j in range(n):
                if i == j:
                    extendPart = np.concatenate((extendPart, [ddq[j], dq[j], np.sign(dq[j]), 1]))
                else:
                    extendPart = np.concatenate((extendPart, [0, 0, 0, 0]))
        extendPart = extendPart.reshape(n * n, 4)
        extendXi = np.concatenate((basePart, extendPart), axis=1)
        return extendXi.reshape(n, n * (nL + 4))

    def getXiNumExCompressed(self, q, dq, ddq):
        """ :returns Compressed extended regressor 5x5*~? in numerical form """
        n = len(q)
        xi = self.getXiNumEx(q, dq, ddq)
        # xi = self.getXiNum(q, dq, ddq)
        cXi = np.delete(xi, self._linerCols, 1)
        return cXi


if __name__ == '__main__':
    xi_num = XiNum(A, D, DELTA)
    # print(xi_num.getXiNumEx([1.1, -0.1], [0.03, -0.7], [0.4, 0.001]))
    q = dq = ddq = (1,2,3,4,5)
    print(len(xi_num.getXiNumExCompressed(q,dq,ddq)[0]))