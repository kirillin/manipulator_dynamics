#!/usr/bin/env python3
"""
    Different stuff for identification of vector chi
"""
import numpy as np
import time

from libs.plot_stuff import plot
from libs.utils import rmZeros, printWastedTime
from xifull import XIFull
from libs.initialization import *



class Identification:
    """
        m -- quantity rows in input data
        n -- 5
    """

    def __init__(self, a, d, theta):
        self.xi = XIFull(a, d, theta)
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
        qs = rawData[0:5].T
        dqs = rawData[5:10].T
        ddqs = rawData[10:15].T
        taus = rawData[15:20].T
        ts = rawData[20].T

        endTime = time.time()
        deltaTime = endTime - timeStart
        printWastedTime(timeStart, endTime, deltaTime,
                             'Reading Ident. data:')
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

    def writeEE(self, fileName, n=14):
        print('Start generation EE files!')
        for i in range(n):
            q = np.random.rand(5) * 4 - 2
            dq = np.random.rand(5) * 4 - 2
            ddq = np.random.rand(5) * 4 - 2

            xi = self.xi.getXiNumExCompressed(q, dq, ddq)
            xiT = xi.T

            fname = fileName.format(i)
            np.savetxt(fname, xiT)
        print('{:d} EE files were generated success!'.format(n))

    def computeDDq(self, dq_prev, dq_next, t_prev, t_next):
        ddq = []
        for i in range(5):
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
        ddqs.append([0., 0., 0., 0., 0.])    # first ddq
        for i, line in enumerate(file):
            raw = line.split(' ')
            q = list(map(float, raw[0:5]))
            dq = list(map(float, raw[5:10]))
            tau = list(map(float, raw[10:15]))
            t = float(raw[15])
            qs.append(q)
            dqs.append(dq)
            taus.append(tau)
            ts.append(t)
            if i > 1:
                dq_prev = dqs[i-2]
                dq_next = dqs[i]
                t_prev = ts[i-2]
                t_next = ts[i]
                ddq = self.computeDDq(dq_prev, dq_next, t_prev, t_next)
                ddqs.append(ddq)

        print('data from ' + fileName + ' was read!')
        file.close()
        ddqs.append([0, 0, 0, 0, 0])
        return qs, dqs, ddqs, taus, ts

if __name__ == '__main__':
    path = 'data_for_identification'

    dataFileName = path + '/data/data_{:d}.txt'
    bigTauFileName = path + '/bigs/big_tau_{:d}.txt'
    bigXiFileName = path+ '/bigs/big_xi_{:d}.txt'
    EEfileName = path + '/ee/EE{:d}.txt'

    A = (0.033, 0.155, 0.135, 0., 0.)
    D = (0.147, 0, 0, 0, 0.218)
    ident = Identification(A, D, thi)
    #
    # q = (0,0,0,0,0)
    # dq = (0, 0, 0, 0, 0)
    # ddq = (0, 0, 0, 0, 0)

    # ar = np.arange(10).reshape(10, 1)
    # taus = np.concatenate((ar,ar,ar,ar,ar), 1)
    #
    # bigXi = ident.getBigXi([q, q], [dq, dq], [ddq, ddq])
    # bigXisLite = ident.getBigXisLite(bigXi)
    # for i in range(2):
    #     print(bigXisLite[i].shape)

    # for i in range(len(bigXi)):
    #     print(bigXisLite[i].shape)

    # q = np.random.rand(5) * 4 - 2
    # dq = np.random.rand(5) * 4 - 2
    # ddq = np.random.rand(5) * 4 - 2
    #
    # taus = np.array([
    #     [1, 2, 3, 4, 5],
    #     [4, 5, 6, 7, 8],
    #     [7, 8, 9, 10, 11],
    #     [10, 11, 12, 12, 13]
    # ])
    #
    # bigXisLite = np.array([
    #     [
    #         [1, 2, 3],
    #         [6, 7, 8],
    #         [11, 12, 13],
    #         [16, 17, 18],
    #         [21, 22, 23]
    #     ],
    #     [
    #         [26, 27, 28],
    #         [31, 32, 33],
    #         [36, 37, 38],
    #         [41, 42, 43],
    #         [46, 47, 48]
    #     ]
    # ])
    #
    # q = (-0.003362, - 0.000146,0.000535, 0.001290,0.004798)
    # dq = (0.001441, -0.000631, - 0.008836,0.051949,0.001979)
    # ddq = (- 0.255993,0.248325, - 0.071348,1.710817,0.779999)
    #
    # bigXi = ident.getBigXi([q,q],[dq,dq],[ddq,ddq])
    #
    # bigXisLite = ident.getBigXisLite(bigXi)
    # print(bigXisLite.shape)
    # print(bigXisLite[0])


    # fileName = 'data_for_identification/data_1_filt.txt'
    # data = ident.readIdentData(fileName, isPlot=False)

    # bigXi = np.ndarray(20).reshape(4,5)
    # fileName = 'data_for_identification/processed_data/big_xi_1.txt'
    # ident.writeBigXi(fileName, bigXi)

    # bigTau = np.ndarray(20).reshape(20,1)
    # fileName = 'data_for_identification/processed_data/big_tau_1.txt'
    # ident.writeBigXi(fileName, bigTau)

    ident.writeEE(EEfileName, 20)