#!/usr/bin/env python3
import numpy as np
import time
import os
import re

EPS = 1.e-15


def printWastedTime(st, et, dt, title=''):
    prepString = title + '\nStart time: {0}\nEnd time: {1}\nWasted time: {2}'
    start = time.strftime("%b %d %H:%M:%S", time.localtime(st))
    end = time.strftime("%b %d %H:%M:%S", time.localtime(et))
    delta = time.strftime("%b %d %H:%M:%S", time.localtime(dt - 3600 * 24))  # -3 hours for 1970 year
    wasted = prepString.format(start, end, delta)
    print(wasted)


def getZeros(mx):
    """
        axis in [1, 2]
            1 -- rows
            2 -- cols
    :param mx:
    :param axis:
    :return:
    """
    m = len(mx[0])
    zeros = []
    for j in range(m):
        if np.all(mx[:, j] < EPS):
            zeros.append(j)
    return zeros


def rmZeros(mx):
    zeroCols = getZeros(mx)
    mx = np.delete(mx, zeroCols, 1)
    return mx


def getDataFileNames(path):
    for dirname, dirnames, filenames in os.walk(path):
        return filenames


def getFileNamesContains(string, path):
    fileNames = getDataFileNames(path)
    wellFileNames = []
    for fn in fileNames:
        pat = r'.*' + string + '.*'
        if re.match(pat, fn) is not None:
            wellFileNames.append(fn)
    return wellFileNames

if __name__ == '__main__':
    path = '/media/data/evo/robotics_report/ros_packages/youbot_arm_control/calculations/data_for_identification/bigs/raw'
    x = getFileNamesContains('', path)
    print(x)
