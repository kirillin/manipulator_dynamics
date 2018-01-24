from numpy import cos, sin, sqrt, tan, zeros, array


class Xi11:

	def __init__(self, q=(0, 0, 0, 0, 0), dq=(0, 0, 0, 0, 0), ddq=(0, 0, 0, 0, 0), a=(0, 0, 0, 0, 0), d=(0, 0, 0, 0, 0), delta=(0, 0, 0, 0, 0)):
		self.q = q
		self.dq = dq
		self.ddq = ddq
		self.a = a
		self.d = d
		self.delta = delta

	def opL0(self, q, dq, ddq, a, d, theta):
		opL_0 = -1.0*a[1]*(1.0*a[0]*sin(theta[1])*dq[0]**2 - 6.12323399573677e-17*a[0]*cos(theta[1])*ddq[0] - 3.74939945665464e-33*a[1]*sin(theta[1])**3*cos(theta[0])**2*cos(theta[1])*dq[0]**2 - 3.74939945665464e-33*a[1]*sin(theta[1])*cos(theta[0])**2*cos(theta[1])**3*dq[0]**2 + 1.0*a[1]*sin(theta[1])*cos(theta[1])*dq[0]**2 + 1.23259516440783e-32*a[1]*cos(theta[0])**4*ddq[0] - 1.23259516440783e-32*a[1]*cos(theta[0])**2*ddq[0] + 6.16297582203915e-33*a[1]*cos(theta[1])**4*ddq[0] - 6.16297582203915e-33*a[1]*cos(theta[1])**2*ddq[0] - 6.12323399573676e-17*a[1]*ddq[0] - 1.0*a[1]*ddq[1] + 9.82*cos(theta[1]))
		return opL_0

	def opL1(self, q, dq, ddq, a, d, theta):
		opL_1 = -1.0*a[0]*sin(theta[1])*dq[0]**2 - 2.29584502165847e-49*a[0]*sin(theta[1])*dq[0]*dq[1] + 6.12323399573677e-17*a[0]*cos(theta[1])*ddq[0] - 1.0*a[1]*sin(2*theta[1])*dq[0]**2 + 1.22464679914735e-16*a[1]*ddq[0] + 2.0*a[1]*ddq[1] - 9.82*cos(theta[1])
		return opL_1

	def opL2(self, q, dq, ddq, a, d, theta):
		opL_2 = -6.12323399573677e-17*a[0]*sin(theta[1])*ddq[0] - 1.0*a[0]*cos(theta[1])*dq[0]**2 - 6.16297582203915e-33*a[1]*sin(2*theta[1])*ddq[0] - 1.0*a[1]*cos(2*theta[1])*dq[0]**2 + 9.82*sin(theta[1])
		return opL_2

	def opL3(self, q, dq, ddq, a, d, theta):
		opL_3 = -1.0*a[1]*(1.74402815495387e-113*sin(theta[0])**2*sin(theta[1])**2*cos(theta[1])*dq[0]**2 + 2.84821412372633e-97*sin(theta[1])**3*ddq[0] + 4.65148665837263e-81*sin(theta[1])**2*cos(theta[1])*dq[0]**2 + 7.59645419660784e-65*sin(theta[1])*ddq[0] + 4.65148665837263e-81*cos(theta[1])**3*dq[0]**2)
		return opL_3

	def opL4(self, q, dq, ddq, a, d, theta):
		opL_4 = sin(2*theta[1])*dq[0]**2/2
		return opL_4

	def opL5(self, q, dq, ddq, a, d, theta):
		opL_5 = -0.5*sin(2*theta[1])*dq[0]**2
		return opL_5

	def opL6(self, q, dq, ddq, a, d, theta):
		opL_6 = 1.0*(6.12323399573677e-17*ddq[0] + 1.0*ddq[1])
		return opL_6

	def opL7(self, q, dq, ddq, a, d, theta):
		opL_7 = -7.49879891330929e-33*cos(theta[0])**4*cos(theta[1])**2*dq[1]**2 + 7.49879891330929e-33*cos(theta[0])**2*cos(theta[1])**2*dq[1]**2 + 2.0*cos(theta[1])**2*dq[0]**2 - 1.0*dq[0]**2
		return opL_7

	def opL8(self, q, dq, ddq, a, d, theta):
		opL_8 = 1.23259516440783e-32*sin(theta[0])**2*cos(theta[1])*dq[1]**2 + 1.0*sin(theta[1])*ddq[0] + 6.12323399573677e-17*cos(theta[1])*dq[0]**2
		return opL_8

	def opL9(self, q, dq, ddq, a, d, theta):
		opL_9 = -6.12323399573677e-17*sin(theta[1])*dq[0]**2 + 1.0*cos(theta[1])*ddq[0]
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

