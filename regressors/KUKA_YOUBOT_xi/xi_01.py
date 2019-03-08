from numpy import cos, sin, sqrt, tan, zeros, array


class Xi01:

	def __init__(self, q=(0, 0, 0, 0, 0), dq=(0, 0, 0, 0, 0), ddq=(0, 0, 0, 0, 0), a=(0, 0, 0, 0, 0), d=(0, 0, 0, 0, 0), delta=(0, 0, 0, 0, 0)):
		self.q = q
		self.dq = dq
		self.ddq = ddq
		self.a = a
		self.d = d
		self.delta = delta

	def opL0(self, q, dq, ddq, a, d, theta):
		opL_0 = 1.0*a[0]**2*ddq[0] - 3.50321773037416e-34*a[0]*a[1]*sin(theta[1])*cos(theta[0])**2*dq[0]*dq[1] + 2.5*a[0]*a[1]*sin(theta[1])*dq[0]*dq[1] + 4.9789962505148e-17*a[0]*a[1]*sin(theta[1])*dq[1]**2 - 3.08148791101958e-33*a[0]*a[1]*cos(theta[1])**3*ddq[1] + 2.0*a[0]*a[1]*cos(theta[1])*ddq[0] + 4.9789962505148e-17*a[0]*a[1]*cos(theta[1])*ddq[1] + 3.39907768361723e-33*a[1]**2*sin(theta[1])**3*cos(theta[0])**2*cos(theta[1])*dq[1]**2 - 1.85833372810744e-33*a[1]**2*sin(theta[1])**3*cos(theta[1])*dq[1]**2 - 3.50321773037416e-34*a[1]**2*sin(theta[1])*cos(theta[0])**2*cos(theta[1])*dq[0]*dq[1] + 6.12323399573677e-17*a[1]**2*sin(theta[1])*cos(theta[0])**2*cos(theta[1])*dq[1]**2 + 2.5*a[1]**2*sin(theta[1])*cos(theta[1])*dq[0]*dq[1] - 4.77797361570133e-17*a[1]**2*sin(theta[1])*cos(theta[1])*dq[1]**2 + 1.0*a[1]**2*cos(theta[1])**2*ddq[0] - 5.72118872610983e-18*a[1]**2*cos(theta[1])**2*ddq[1] + 5.55111512312578e-17*a[1]**2*ddq[1]
		return opL_0

	def opL1(self, q, dq, ddq, a, d, theta):
		opL_1 = 2.0*a[0]*sin(theta[1])*dq[0]*dq[1] + 6.12323399573677e-17*a[0]*sin(theta[1])*dq[1]**2 + 2.0*a[0]*cos(theta[1])*ddq[0] + 5.55111512312578e-17*a[0]*cos(theta[1])*ddq[1] + 2.0*a[1]*sin(2*theta[1])*dq[0]*dq[1] - 4.29089154458236e-18*a[1]*sin(2*theta[1])*dq[1]**2 - 7.65404249467096e-18*a[1]*sin(2*theta[0] - 2*theta[1])*dq[1]**2 + 7.65404249467096e-18*a[1]*sin(2*theta[0] + 2*theta[1])*dq[1]**2 + 1.0*a[1]*cos(2*theta[1])*ddq[0] - 2.86059436305491e-18*a[1]*cos(2*theta[1])*ddq[1] + 1.0*a[1]*ddq[0] + 1.13882896825571e-16*a[1]*ddq[1]
		return opL_1

	def opL2(self, q, dq, ddq, a, d, theta):
		opL_2 = -2.0*a[0]*sin(theta[1])*ddq[0] - 5.55111512312578e-17*a[0]*sin(theta[1])*ddq[1] + 2.0*a[0]*cos(theta[1])*dq[0]*dq[1] + 6.12323399573677e-17*a[0]*cos(theta[1])*dq[1]**2 + 6.12323399573677e-17*a[1]*sin(theta[0])**2*sin(theta[1])**2*dq[1]**2 - 4.0*a[1]*sin(theta[1])**2*dq[0]*dq[1] - 2.20343868895191e-17*a[1]*sin(theta[1])**2*dq[1]**2 - 2.0*a[1]*sin(theta[1])*cos(theta[1])*ddq[0] + 5.72118872610982e-18*a[1]*sin(theta[1])*cos(theta[1])*ddq[1] + 2.0*a[1]*dq[0]*dq[1] + 1.23259516440783e-32*a[1]*dq[1]**2
		return opL_2

	def opL3(self, q, dq, ddq, a, d, theta):
		opL_3 = 3.1758977259765e-34*a[0]*dq[1]**2 - 5.55111512312578e-17*a[1]*sin(theta[1])**3*ddq[1] - 1.11022302462516e-16*a[1]*sin(theta[1])*ddq[0] - 1.0*a[1]*sin(theta[1])*ddq[1] + 1.11022302462516e-16*a[1]*cos(theta[1])*dq[0]*dq[1] + 1.0*a[1]*cos(theta[1])*dq[1]**2
		return opL_3

	def opL4(self, q, dq, ddq, a, d, theta):
		opL_4 = -1.0*sin(2*theta[1])*dq[0]*dq[1] - 2.86059436305492e-18*sin(2*theta[1])*dq[1]**2 - 0.5*cos(2*theta[1])*ddq[0] + 0.5*ddq[0]
		return opL_4

	def opL5(self, q, dq, ddq, a, d, theta):
		opL_5 = 1.0*sin(2*theta[1])*dq[0]*dq[1] + 2.86059436305492e-18*sin(2*theta[1])*dq[1]**2 + 0.5*cos(2*theta[1])*ddq[0] + 0.5*ddq[0]
		return opL_5

	def opL6(self, q, dq, ddq, a, d, theta):
		opL_6 = 1.0*(3.74939945665464e-33*ddq[0] + 6.12323399573677e-17*ddq[1])
		return opL_6

	def opL7(self, q, dq, ddq, a, d, theta):
		opL_7 = 1.0*sin(2*theta[1])*ddq[0] - 2.0*cos(2*theta[1])*dq[0]*dq[1] - 5.72118872610983e-18*cos(2*theta[1])*dq[1]**2
		return opL_7

	def opL8(self, q, dq, ddq, a, d, theta):
		opL_8 = 1.22464679914735e-16*sin(theta[1])*ddq[0] + 1.0*sin(theta[1])*ddq[1] - 1.22464679914735e-16*cos(theta[1])*dq[0]*dq[1] - 1.0*cos(theta[1])*dq[1]**2
		return opL_8

	def opL9(self, q, dq, ddq, a, d, theta):
		opL_9 = 1.22464679914735e-16*sin(theta[1])*dq[0]*dq[1] + 1.0*sin(theta[1])*dq[1]**2 + 1.22464679914735e-16*cos(theta[1])*ddq[0] + 1.0*cos(theta[1])*ddq[1]
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

