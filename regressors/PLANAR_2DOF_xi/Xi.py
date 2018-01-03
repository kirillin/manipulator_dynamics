import numpy as np
from regressors.PLANAR_2DOF_xi.xi_00 import Xi00
from regressors.PLANAR_2DOF_xi.xi_01 import Xi01
from regressors.PLANAR_2DOF_xi.xi_11 import Xi11

class Xi:

	def __init__(self, a, d):
		xi0 = [Xi00(a=a, d=d), Xi01(a=a, d=d)]
		xi1 = [				0, Xi11(a=a, d=d)]
		self._xi = [xi0, xi1]
		self._linerCols = [3, 4, 5, 7, 8, 9, 12, 14, 16, 17, 18, 19]

	def getXi(self, q, dq, ddq):
		n = len(q)
		xi = np.empty(0)
		for i in range(n):
			for j in range(n):
				if self._xi[i][j] != 0:
					xi_ij = self._xi[i][j].getXi(q, dq, ddq)
					xi = np.concatenate((xi, xi_ij))
				else:
					xi = np.concatenate((xi, np.zeros(10)))
		return xi.reshape(n, n * 10)

	def getWellColNums(self):
		all = set(range(20))
		well = all.difference(self.__linerCols)
		return list(well)

