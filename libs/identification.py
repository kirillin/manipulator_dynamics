#!/usr/bin/env python3
"""
    Different stuff for identification of vector chi
"""
import time

import numpy as np

from libs.initialization import *
from libs.plot_stuff import plot
from libs.utils import rmZeros, printWastedTime
from regressors.xi_num import XiNum


class Identification:
    """
        m -- quantity rows in input data
        n -- quantity joints
    """

    def __init__(self, a, d, delta):
        self.xi = XiNum(a, d, delta)
        self.wellColNums = self.xi.getWellColNums()
        # DH parameters
        self.a = a
        self.d = d

    def getWellColNums(self):
        return self.wellColNums

    def getBigTau(self, taus):
        """
            TESTED!
            :param taus: [m x n]
            :return: [n*m, 1]
        """
        m = len(taus)
        n = len(taus[0])
        bigTau = np.empty(0)
        for k in range(m):
            bigTau = np.concatenate((bigTau, taus[k]), 0)
        return bigTau.reshape(n * m, 1)

    def getBigTausLite(self, taus):
        """
        TESTED!
        :param taus: [m x n]
        :return: [n x m] !
        """
        return np.array(taus)

    def getCalcTausLite(self, bigXisLite, estChisLite):
        """
        :param bigXisLite: [n x [m/n x ~70-?]]
        :param estChis: [~70 x n]
        :return: [n x [m/n x 1]]
        """
        n = len(bigXisLite)
        calcTaus = []
        for i in range(n):
            # TODO NOT WORK MULTIPLY
            e = bigXisLite[i] * np.array([estChisLite[i]])
            calcTaus.append([e])
        return np.array(calcTaus)

    def getBigXi(self, qs, dqs, ddqs):
        """
        TESTED!
        :param qs: [m x n]
        :param dqs: [m x n]
        :param ddqs: [m x n]
        :return: [n*m x ~70-?]
        """
        print('Getting big xi was started!')
        startTime = time.time()
        m = len(qs)
        bigXi = self.xi.getXiNumExCompressed(qs[0], dqs[0], ddqs[0])
        for k in range(1, m):
            xiK = self.xi.getXiNumExCompressed(qs[k], dqs[k], ddqs[k])
            bigXi = np.concatenate((bigXi, xiK), axis=0)
            msgInfo = 'Collecting of bigXi. Progress: {:.2f} %'.format(k * 100. / m)
            print(msgInfo)   # print(msgInfo, end='\n')
        endTime = time.time()
        printWastedTime(startTime, endTime, startTime - endTime,
                                'Make big Xi:')
        return bigXi

    # def getBigXisLite(self, qs, dqs, ddqs):
    #     print('Getting big xis lite (compute) was started!')
    #     startTime = time.time()
    #     m = len(qs)
    #     n = len(qs[0])
    #     bigXi = np.array(self.xi.getXiNumExCompressed(qs[0], dqs[0], ddqs[0]))
    #     bigXis = [rmZeros(xi_l) for xi_l in bigXi]
    #     for k in range(1, m):
    #         xiK = self.xi.getXiNumExCompressed(qs[k], dqs[k], ddqs[k])
    #         bigXis = [rmZeros(xi_l) for xi_l in bigK]
    #
    #         msgInfo = 'Collecting of bigXi. Progress: {:.2f} %'.format(k * 100. / m)
    #         print(msgInfo, end='\n')
    #
    #     endTime = time.time()
    #     printWastedTime(startTime, endTime, startTime - endTime,
    #                             'Make big Xis lite (compute):')
    #     return bigXis

    def getBigXisLite(self, bigXi):
        """
        TESTED!
        :param bigXi: [m x ~70]
        :return: [n x [m/n x ~70-?]]
        """
        print('Getting big xis lite (split exists) was started!')
        startTime = time.time()
        n = 5
        m = len(bigXi)
        bigXis = []
        for j in range(n):
            bigXisCur = np.array([bigXi[i+j, :] for i in range(m-j) if i % 5 == 0])
            bigXisLite = rmZeros(bigXisCur)

            # zerosCols = []
            # p = len(bigXisCur[0])
            # r = len(bigXisCur)
            # for i in range(p):
            #     flag = True
            #     for j in range(r):
            #         if bigXisCur[j, i] != 0.:
            #             flag = False
            #     if flag == True:
            #         zerosCols.append(i)
            # bigXisCur = np.delete(bigXisCur, zerosCols, 1)
            bigXis.append(bigXisCur)

        endTime = time.time()
        printWastedTime(startTime, endTime, startTime - endTime,
                                'Make big Xis lite (split exists):')
        return bigXis

    def getEstChisLite(self, bigXisLite, bigTausLite):
        """
        :param bigXisLite: [n x [m/n x ~70-?]]
        :param bigTausLite: [n x m]
        :return: [n x ~70]
        """
        print('Start Chis estimate! ')
        n = len(bigXisLite)
        estChis = []
        for i in range(n):
            print(i)
            bigXiT = np.transpose(bigXisLite[i])
            prep = np.matmul(np.linalg.inv(bigXiT.dot(bigXisLite[i])), bigXiT)
            print(prep.shape)
            b = np.array(bigTausLite[0][i])
            print(b.shape)
            e = np.dot(prep, np.transpose(b))
            estChis.append(e)
        return np.array(estChis).T

    def getVariance(self, bigTaus, calsTaus):
        m = len(bigTaus[0])
        n = len(bigTaus)    # 5
        sigma = np.zeros(n)
        for i in range(n):
            diff = bigTaus[i] - calsTaus[i]
            sigma[i] = np.transpose(diff) * diff / m
        return sigma

    def getR(self, sigma, m):
        n = 5
        Isigma = sigma**2 * np.eye(n)
        Z = np.zeros((n, n))
        R = np.empty(0)
        for i in range(m):
            for j in range(m):
                if j == i:
                    Rrow = np.column_stack((Rrow, Isigma))
                else:
                    Rrow = np.column_stack((Rrow, Z))
        if i > 0:
            R = np.row_stack((R, Rrow))
        else:
            R = Rrow
        return R

    def readIdentData(self, fileName, isPlot=False):
        timeStart = time.time()
        file = open(fileName, 'r')
        rawData = np.loadtxt(file)
        rawData = rawData.T
        qs = rawData[0:n].T
        dqs = rawData[n:2*n].T
        ddqs = rawData[2*n:3*n].T
        taus = rawData[3*n:4*n].T
        ts = rawData[4*n].T

        endTime = time.time()
        deltaTime = endTime - timeStart
        printWastedTime(timeStart, endTime, deltaTime, 'Reading Ident. data:')
        if isPlot:
            plot(qs, dqs, ddqs, taus, ts)
        return qs, dqs, ddqs, taus, ts

    def readBigXi(self, fileName):
        file = open(fileName, 'r')
        data = np.loadtxt(file)
        return data

    def readBigTau(self, fileName):
        file = open(fileName, 'r')
        data = np.loadtxt(file)
        return data

    def writeBigXi(self, fileName, bigXi):
        np.savetxt(fileName, bigXi)
        print('Big Xi was wrote!')

    def writeBigTau(self, fileName, bigTau):
        np.savetxt(fileName, bigTau)
        print('Big Tau was wrote!')

    def writeEE2(self, fileName, k=14):
        print('Start generation EE files!')
        bigXi = np.zeros(k*n)

        q = np.random.rand(n) * 4 - 2
        dq = np.random.rand(n) * 4 - 2
        ddq = np.random.rand(n) * 4 - 2
        bigXi = self.xi.getXiNumExCompressed(q, dq, ddq)

        for i in range(1, k):
            print(i*100/k)
            q = np.random.rand(n) * 4 - 2
            dq = np.random.rand(n) * 4 - 2
            ddq = np.random.rand(n) * 4 - 2

            xi = self.xi.getXiNumExCompressed(q, dq, ddq)
            bigXi = np.concatenate((bigXi, xi), axis=0)

        np.savetxt(fileName, bigXi.T)

        print('{:d} EE files were generated success!'.format(k))
        print('Liner cols: {0}'.format(str(self.xi.getLinerCols())))
        print('Well cols: {0}; qty: {1}'.format(str(self.xi.getWellColNums()), len(self.xi.getWellColNums())))

    def writeEE(self, fileName, k=14):
        print('Start generation EE files!')
        status = 0
        for i in range(k):
            print(i*100/k)
            q = np.random.rand(n) * 4 - 2
            dq = np.random.rand(n) * 4 - 2
            ddq = np.random.rand(n) * 4 - 2

            xi = self.xi.getXiNumExCompressed(q, dq, ddq)
            xiT = xi.T

            fname = fileName.format(i)
            np.savetxt(fname, xiT)
        print('{:d} EE files were generated success!'.format(k))
        print('Liner cols: {0}'.format(str(self.xi.getLinerCols())))
        print('Well cols: {0}; qty: {1}'.format(str(self.xi.getWellColNums()), len(self.xi.getWellColNums())))

    def computeDDq(self, dq_prev, dq_next, t_prev, t_next):
        ddq = []
        for i in range(n):
            lim = (dq_next[i] - dq_prev[i]) / (t_next - t_prev)
            ddq.append(lim)
        return ddq

    def readIdentDataRaw(self, fileName):
        """
        data measures from real manipulator;
        :return: q, dq, ddq, tau_e -- matrices of vectors
        """
        file = open(fileName, 'r')
        qs, dqs, ddqs, taus, ts = [], [], [], [], []
        ddqs.append([0. for i in range(n)])    # first ddq
        for i, line in enumerate(file):
            raw = line.split(' ')
            q = list(map(float, raw[0:n]))
            dq = list(map(float, raw[n:2*n]))
            tau = list(map(float, raw[2*n:3*n]))
            t = float(raw[3*n])
            qs.append(q)
            dqs.append(dq)
            taus.append(tau)
            ts.append(t)
            if i > 1:
                dq_prev = dqs[i-2]
                dq_next = dqs[i]
                t_prev = ts[i-2]
                t_next = ts[i]
                ddq = self.computeDDq(dq_prev, dq_next, t_prev, t_next+0.0000000001)
                ddqs.append(ddq)

        print('data from ' + fileName + ' was read!')
        file.close()
        ddqs.append([0, 0, 0, 0, 0])
        return qs, dqs, ddqs, taus, ts
