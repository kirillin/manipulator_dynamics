# coding: utf-8
#!/usr/bin/env python3
import numpy as np

from identification import Identification
from libs.plot_stuff import plotTaus2, plotTaus, plotChis
from libs.utils import getFileNamesContains
from libs.initialization import *

A = (0.033, 0.155, 0.135, 0., 0.)
D = (0.147, 0, 0, 0, 0.218)

PATH = 'data_for_identification'

DATA_FILE_NAME = PATH + '/data/data_{:d}{:s}.txt'
BIG_TAUs_FILE_NAME = PATH + '/bigs/{:s}big_tau{:d}.txt'
BIG_XIs_FILE_NAME = PATH + '/bigs/{:s}big_xi{:d}.txt'
EE_FILE_NAME = PATH + '/ee/EE{:d}.txt'

def generateBigXisRaw(ident):
    for i in range(1, 2):
        datafileName = DATA_FILE_NAME.format(i, '')
        bigXiFileName = BIG_XIs_FILE_NAME.format('raw/', i)
        bigTauFileName = BIG_TAUs_FILE_NAME.format('raw/', i)
        print('Read {0} filt file {1}!'.format(i, datafileName))
        # # compute matrix
        Q, dQ, ddQ, Tau, T = ident.readIdentDataRaw(datafileName)
        print('{0} {1} {2} {3}'.format(len(Q), len(dQ), len(ddQ), len(Tau)))

        bigXi = ident.getBigXi(Q, dQ, ddQ)
        bigTau = ident.getBigTau(Tau)

        ident.writeBigXi(bigXiFileName, bigXi)
        ident.writeBigTau(bigTauFileName, bigTau)

def generateBigXisFilt(ident):
    for i in range(1, 2):
        datafileName = DATA_FILE_NAME.format(i, '_filt')
        bigXiFileName = BIG_XIs_FILE_NAME.format('filt/', i)
        bigTauFileName = BIG_TAUs_FILE_NAME.format('filt/', i)
        print('Read {0} filt file {1}!'.format(i, datafileName))
        # # compute matrix
        Q, dQ, ddQ, Tau, T = ident.readIdentData(datafileName)

        # ident.writeEE(EE_FILE_NAME, 50)
        bigXi = ident.getBigXi(Q, dQ, ddQ)
        bigTau = ident.getBigTau(Tau)

        ident.writeBigXi(bigXiFileName, bigXi)
        ident.writeBigTau(bigTauFileName, bigTau)


def generateBigXis(dataPath, output='raw'):
    fileNames = getFileNamesContains('date.*'+output, dataPath)
    for fn in fileNames:
        print('Read {0} file!'.format(fn))
        Q, dQ, ddQ, Tau, T = ident.readIdentData(fn)

        bigXi = ident.getBigXi(Q, dQ, ddQ)
        bigTau = ident.getBigTau(Tau)

        bigXiFileName = outputPath
        ident.writeBigXi(bigXiFileName, bigXi)
        ident.writeBigTau(bigTauFileName, bigTau)


def getSD(chis):
    qtyChis =  len(chis)
    qtyCols = len(chis[0])
    avg = [i for i in range(qtyCols)]
    sd = [i for i in range(qtyCols)]
    for i in range(qtyCols):
        avg[i] = sum([raw[i] for raw in chis]) / qtyChis
        sd[i] = np.sqrt(sum([(row[i] - avg[i])**2 for row in chis]) / qtyChis)

    sdPer = [i for i in range(qtyCols)]
    for i in range(qtyCols):
        sdPer[i] = sd[i]  / abs(avg[i]) * 100

    return sd, sdPer, avg

if __name__ == '__main__':
    ident = Identification(A, D, thi)

    # 1. Выбросить лин. завис. столбцы
    # ident.writeEE(EE_FILE_NAME, 50)

    # 2.
    generateBigXisRaw(ident)

    # Chis = []
    # for i in range(1,11):
    #     datafileName = DATA_FILE_NAME.format(i, '')
    #     bigXiFileName = BIG_XIs_FILE_NAME.format('raw/', i)
    #     bigTauFileName = BIG_TAUs_FILE_NAME.format('raw/', i)
    #     #
    #     Q, dQ, ddQ, Tau, T = ident.readIdentDataRaw(datafileName)
    #     #
    #     bigXi = ident.readBigXi(bigXiFileName)
    #     bigTau = ident.readBigTau(bigTauFileName)
    #
    #     # bigXisLite = ident.getBigXisLite(bigXi)
    #     # bigTausLite = ident.getBigTausLite(Tau)
    #     # estChiLite = ident.getEstChisLite(bigXisLite,bigTausLite)
    #     # caltTaus = ident.getCalcTausLite(bigXisLite, estChiLite)
    #     #
    #     # print(caltTaus.shape)
    #     #
    #     # plotTaus2(Tau, caltTaus, T, T)
    #
    #     # prepareBigXi = np.linalg.inv(bigXiT.dot(bigXi)).dot(bigXiT)
    #     # estChi = prepareBigXi *(bigTau)
    #
    #     bigXiT = np.transpose(bigXi)
    #     prepareBigXi = np.matmul(np.linalg.inv(np.matmul(bigXiT, bigXi)), bigXiT)
    #
    #     estChi = np.matmul(prepareBigXi, bigTau)
    #     print(len(estChi))
    #     Chis.append(estChi)
    #
    #     estTau = bigXi.dot(estChi)
    #     estTau = np.array(np.hsplit(estTau, len(estTau) / 5))
    #
    #     Tau = np.array(Tau)
    #     #plotTaus(Tau, estTau, T, T)

    # fileNameChis = 'chis_filt.txt'
    # # np.savetxt(fileNameChis, np.array(Chis))
    #
    # # FILTRED DATA
    # f = open(fileNameChis, 'r')
    # chisFilt = np.loadtxt(f)
    # sdFilt = getSD(chisFilt)
    # print(sdFilt)
    #
    # # RAW DATA
    # f = open('chis_raw.txt', 'r')
    # chis = np.loadtxt(f)
    # sdRaw = getSD(chis)
    # print(sdRaw)
    #
    # plotChis(chis, ident.getWellColNums(), chisFilt, sd={'raw': sdRaw[0], 'filt': sdFilt[0]},
    #                                                 avg={'raw': sdRaw[2], 'filt': sdFilt[2]},
    #                                                 sdP=(sdRaw[1], sdFilt[1]))

