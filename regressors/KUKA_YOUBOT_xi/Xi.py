import numpy as np
from regressors.KUKA_YOUBOT_xi.xi_00 import Xi00
from regressors.KUKA_YOUBOT_xi.xi_01 import Xi01
from regressors.KUKA_YOUBOT_xi.xi_02 import Xi02
from regressors.KUKA_YOUBOT_xi.xi_03 import Xi03
from regressors.KUKA_YOUBOT_xi.xi_04 import Xi04
from regressors.KUKA_YOUBOT_xi.xi_11 import Xi11
from regressors.KUKA_YOUBOT_xi.xi_12 import Xi12
from regressors.KUKA_YOUBOT_xi.xi_13 import Xi13
from regressors.KUKA_YOUBOT_xi.xi_14 import Xi14
from regressors.KUKA_YOUBOT_xi.xi_22 import Xi22
from regressors.KUKA_YOUBOT_xi.xi_23 import Xi23
from regressors.KUKA_YOUBOT_xi.xi_24 import Xi24
from regressors.KUKA_YOUBOT_xi.xi_33 import Xi33
from regressors.KUKA_YOUBOT_xi.xi_34 import Xi34
from regressors.KUKA_YOUBOT_xi.xi_44 import Xi44

class Xi:

	def __init__(self, a, d, delta):
		xi0 = [Xi00(a=a, d=d, delta=delta), Xi01(a=a, d=d, delta=delta), Xi02(a=a, d=d, delta=delta), Xi03(a=a, d=d, delta=delta), Xi04(a=a, d=d, delta=delta)]
		xi1 = [				0, Xi11(a=a, d=d, delta=delta), Xi12(a=a, d=d, delta=delta), Xi13(a=a, d=d, delta=delta), Xi14(a=a, d=d, delta=delta)]
		xi2 = [				0, 				0, Xi22(a=a, d=d, delta=delta), Xi23(a=a, d=d, delta=delta), Xi24(a=a, d=d, delta=delta)]
		xi3 = [				0, 				0, 				0, Xi33(a=a, d=d, delta=delta), Xi34(a=a, d=d, delta=delta)]
		xi4 = [				0, 				0, 				0, 				0, Xi44(a=a, d=d, delta=delta)]
		self._xi = [xi0, xi1, xi2, xi3, xi4]
		self._linerCols = []

	def getWellColNums(self):
		all = set(range(50))
		well = all.difference(self._linerCols)
		return list(well)

