from numpy import cos, sin, sqrt, tan, zeros, array


class Xi00:

	def __init__(self, q=(0, 0, 0, 0, 0), dq=(0, 0, 0, 0, 0), ddq=(0, 0, 0, 0, 0), a=(0, 0, 0, 0, 0), d=(0, 0, 0, 0, 0), delta=(0, 0, 0, 0, 0)):
		self.q = q
		self.dq = dq
		self.ddq = ddq
		self.a = a
		self.d = d
		self.delta = delta

	def opL0(self, q, dq, ddq, a, d, theta):
		opL_0 = 1.0*a[0]**2*ddq[0]
		return opL_0

	def opL1(self, q, dq, ddq, a, d, theta):
		opL_1 = 2.0*a[0]*ddq[0]
		return opL_1

	def opL2(self, q, dq, ddq, a, d, theta):
		opL_2 = 0
		return opL_2

	def opL3(self, q, dq, ddq, a, d, theta):
		opL_3 = 0
		return opL_3

	def opL4(self, q, dq, ddq, a, d, theta):
		opL_4 = 0
		return opL_4

	def opL5(self, q, dq, ddq, a, d, theta):
		opL_5 = 1.0*ddq[0]
		return opL_5

	def opL6(self, q, dq, ddq, a, d, theta):
		opL_6 = 3.74939945665464e-33*ddq[0]
		return opL_6

	def opL7(self, q, dq, ddq, a, d, theta):
		opL_7 = 0
		return opL_7

	def opL8(self, q, dq, ddq, a, d, theta):
		opL_8 = 0
		return opL_8

	def opL9(self, q, dq, ddq, a, d, theta):
		opL_9 = 1.22464679914735e-16*ddq[0]
		return opL_9

	def getXi(self, q, dq, ddq):
		XI = zeros(10)
		theta = array(self.delta) - array(q)
		XI[0] = self.opL0(q, dq, ddq, self.a, self.d, theta)
		XI[1] = self.opL1(q, dq, ddq, self.a, self.d, theta)
		XI[2] = self.opL2(q, dq, ddq, self.a, self.d, theta)
		XI[3] = self.opL3(q, dq, ddq, self.a, self.d, theta)
		XI[4] = self.opL4(q, dq, ddq, self.a, self.d, theta)
		XI[5] = self.opL5(q, dq, ddq, self.a, self.d, theta)
		XI[6] = self.opL6(q, dq, ddq, self.a, self.d, theta)
		XI[7] = self.opL7(q, dq, ddq, self.a, self.d, theta)
		XI[8] = self.opL8(q, dq, ddq, self.a, self.d, theta)
		XI[9] = self.opL9(q, dq, ddq, self.a, self.d, theta)
		return XI

