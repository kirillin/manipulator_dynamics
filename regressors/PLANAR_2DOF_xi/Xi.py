import numpy as np
from regressors.PLANAR_2DOF_xi.xi_00 import Xi00
from regressors.PLANAR_2DOF_xi.xi_01 import Xi01
from regressors.PLANAR_2DOF_xi.xi_11 import Xi11

class Xi:

	def __init__(self, a, d, delta):
		xi0 = [Xi00(a=a, d=d, delta=delta), Xi01(a=a, d=d, delta=delta)]
		xi1 = [				0, Xi11(a=a, d=d, delta=delta)]
		self._xi = [xi0, xi1]
		self._linerCols = [3, 4, 5, 7, 8, 9, 14, 17, 18]

	def getWellColNums(self):
		all = set(range(20))
		well = all.difference(self._linerCols)
		return list(well)

