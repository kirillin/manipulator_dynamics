from numpy import cos, sin, sqrt, tan, zeros


class Xi00:

	def __init__(self, q=(0, 0), dq=(0, 0), ddq=(0, 0), a=(0, 0), d=(0, 0), delta=(0, 0)):
		self.q = q
		self.dq = dq
		self.ddq = ddq
		self.a = a
		self.d = d
		self.delta = delta

	def opL0(self, q, dq, ddq, a, d, delta):
		opL_0 = (-a[0]*sin(delta[0] - q[0])**2 - a[0]*cos(delta[0] - q[0])**2)*(-a[0]*sin(delta[0] - q[0])**2*ddq[0]/2 - a[0]*cos(delta[0] - q[0])**2*ddq[0]/2) + (-a[0]*sin(delta[0] - q[0])**2/2 - a[0]*cos(delta[0] - q[0])**2/2)*(-a[0]*sin(delta[0] - q[0])**2*ddq[0] - a[0]*cos(delta[0] - q[0])**2*ddq[0]) - 9.82*(a[0]*sin(delta[0] - q[0])**2 + a[0]*cos(delta[0] - q[0])**2)*cos(delta[0] - q[0])
		return opL_0

	def opL1(self, q, dq, ddq, a, d, delta):
		opL_1 = a[0]*sin(delta[0] - q[0])**2*ddq[0] + a[0]*cos(delta[0] - q[0])**2*ddq[0] - (-a[0]*sin(delta[0] - q[0])**2 - a[0]*cos(delta[0] - q[0])**2)*ddq[0] - 9.82*cos(delta[0] - q[0])
		return opL_1

	def opL2(self, q, dq, ddq, a, d, delta):
		opL_2 = 9.82*sin(delta[0] - q[0])
		return opL_2

	def opL3(self, q, dq, ddq, a, d, delta):
		opL_3 = 0
		return opL_3

	def opL4(self, q, dq, ddq, a, d, delta):
		opL_4 = 0
		return opL_4

	def opL5(self, q, dq, ddq, a, d, delta):
		opL_5 = 0
		return opL_5

	def opL6(self, q, dq, ddq, a, d, delta):
		opL_6 = ddq[0]
		return opL_6

	def opL7(self, q, dq, ddq, a, d, delta):
		opL_7 = 0
		return opL_7

	def opL8(self, q, dq, ddq, a, d, delta):
		opL_8 = 0
		return opL_8

	def opL9(self, q, dq, ddq, a, d, delta):
		opL_9 = 0
		return opL_9

	def getXi(self, q, dq, ddq):
		XI = zeros(10)
		XI[0] = self.opL0(q, dq, ddq, self.a, self.d, self.delta)
		XI[1] = self.opL1(q, dq, ddq, self.a, self.d, self.delta)
		XI[2] = self.opL2(q, dq, ddq, self.a, self.d, self.delta)
		XI[3] = self.opL3(q, dq, ddq, self.a, self.d, self.delta)
		XI[4] = self.opL4(q, dq, ddq, self.a, self.d, self.delta)
		XI[5] = self.opL5(q, dq, ddq, self.a, self.d, self.delta)
		XI[6] = self.opL6(q, dq, ddq, self.a, self.d, self.delta)
		XI[7] = self.opL7(q, dq, ddq, self.a, self.d, self.delta)
		XI[8] = self.opL8(q, dq, ddq, self.a, self.d, self.delta)
		XI[9] = self.opL9(q, dq, ddq, self.a, self.d, self.delta)
		return XI

