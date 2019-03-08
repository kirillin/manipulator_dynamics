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
		opL_0 = 4.16333634234434e-17*a[0]**2*sin(2*theta[1])*dq[0]**2 - 1.0*a[0]*a[1]*sin(theta[1])*dq[0]**2 - 696236958976523.0*a[0]*a[1]*sin(theta[1])*dq[0]*dq[1] + 2.08166817117217e-17*a[0]*a[1]*sin(3*theta[1])*dq[0]*dq[1] + 5.55111512312578e-17*a[0]*a[1]*cos(theta[1])*ddq[0] - 7.70371977754894e-34*a[0]*a[1]*cos(3*theta[1])*ddq[0] - 0.5*a[1]**2*sin(2*theta[1])*dq[0]**2 - 348118479488261.0*a[1]**2*sin(2*theta[1])*dq[0]*dq[1] + 0.00199165152019878*a[1]**2*sin(2*theta[1])*dq[1]**2 + 2.6016685033094e-50*a[1]**2*sin(4*theta[1])*dq[1]**2 - 7.70371977754894e-34*a[1]**2*cos(2*theta[1])**2*ddq[1] + 1.04066740132376e-49*a[1]**2*cos(2*theta[1])*ddq[0] + 1.87469972832732e-33*a[1]**2*cos(2*theta[1])*ddq[1] + 5.55111512312578e-17*a[1]**2*ddq[0] + 1.0*a[1]**2*ddq[1] - 9.82*a[1]*cos(theta[1]) - 9.20477566608715e-33*a[1]*cos(3*theta[1])
		return opL_0

	def opL1(self, q, dq, ddq, a, d, theta):
		opL_1 = -1.0*a[0]*sin(theta[1])*dq[0]**2 + 464157972651015.0*a[0]*sin(theta[1])*dq[0]*dq[1] + 5.55111512312578e-17*a[0]*cos(theta[1])*ddq[0] - 1.0*a[1]*sin(2*theta[1])*dq[0]**2 - 0.00132776768013244*a[1]*sin(2*theta[1])*dq[1]**2 + 1.14792251082923e-49*a[1]*cos(2*theta[1])*ddq[0] + 1.87469972832732e-33*a[1]*cos(2*theta[1])*ddq[1] + 1.16743491188625e-16*a[1]*ddq[0] + 2.0*a[1]*ddq[1] - 9.82*cos(theta[1])
		return opL_1

	def opL2(self, q, dq, ddq, a, d, theta):
		opL_2 = -5.55111512312578e-17*a[0]*sin(theta[1])*ddq[0] - 1.0*a[0]*cos(theta[1])*dq[0]**2 + 464157972651014.0*a[0]*cos(theta[1])*dq[0]*dq[1] - 1.69953884180861e-33*a[1]*sin(2*theta[1])*ddq[0] - 5.55111512312578e-17*a[1]*sin(2*theta[1])*ddq[1] - 1.0*a[1]*cos(2*theta[1])*dq[0]**2 - 0.00132776768013249*a[1]*cos(2*theta[1])*dq[1]**2 + 464157972651015.0*a[1]*dq[0]*dq[1] - 0.00132776768013249*a[1]*dq[1]**2 + 9.82*sin(theta[1])
		return opL_2

	def opL3(self, q, dq, ddq, a, d, theta):
		opL_3 = 1.0*a[0]*cos(theta[1])**2*dq[0]**2 + 0.0257659434150245*a[0]*dq[0]*dq[1] - 5.55111512312578e-17*a[1]*sin(theta[1])**3*ddq[0] - 1.0*a[1]*sin(theta[1])*ddq[0] + 5.55111512312578e-17*a[1]*cos(theta[1])**3*dq[0]*dq[1] - 0.0257659434150245*a[1]*cos(theta[1])**3*dq[1]**2 - 5.55111512312578e-17*a[1]*cos(theta[1])*dq[0]**2 + 0.0257659434150245*a[1]*cos(theta[1])*dq[0]*dq[1] + 464157972651015.0*a[1]*cos(theta[1])*dq[1]**2
		return opL_3

	def opL4(self, q, dq, ddq, a, d, theta):
		opL_4 = 1.0*(0.5*dq[0] - 116039493162754.0*dq[1])*sin(2*theta[1])*dq[0]
		return opL_4

	def opL5(self, q, dq, ddq, a, d, theta):
		opL_5 = -1.0*(0.5*dq[0] - 116039493162754.0*dq[1])*sin(2*theta[1])*dq[0]
		return opL_5

	def opL6(self, q, dq, ddq, a, d, theta):
		opL_6 = 1.0*(6.12323399573677e-17*ddq[0] + 1.0*ddq[1])
		return opL_6

	def opL7(self, q, dq, ddq, a, d, theta):
		opL_7 = 1.0*(1.0*dq[0] - 464157972651015.0*dq[1])*cos(2*theta[1])*dq[0]
		return opL_7

	def opL8(self, q, dq, ddq, a, d, theta):
		opL_8 = 1.0*sin(theta[1])*ddq[0] + 6.12323399573677e-17*cos(theta[1])*dq[0]**2 - 0.0284214787752895*cos(theta[1])*dq[0]*dq[1] - 464157972651015.0*cos(theta[1])*dq[1]**2
		return opL_8

	def opL9(self, q, dq, ddq, a, d, theta):
		opL_9 = -6.12323399573677e-17*sin(theta[1])*dq[0]**2 + 0.0284214787752894*sin(theta[1])*dq[0]*dq[1] + 464157972651015.0*sin(theta[1])*dq[1]**2 + 1.0*cos(theta[1])*ddq[0]
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

