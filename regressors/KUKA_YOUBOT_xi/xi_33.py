from numpy import cos, sin, sqrt, tan, zeros, array


class Xi33:

	def __init__(self, q=(0, 0, 0, 0, 0), dq=(0, 0, 0, 0, 0), ddq=(0, 0, 0, 0, 0), a=(0, 0, 0, 0, 0), d=(0, 0, 0, 0, 0), delta=(0, 0, 0, 0, 0)):
		self.q = q
		self.dq = dq
		self.ddq = ddq
		self.a = a
		self.d = d
		self.delta = delta

	def opL0(self, q, dq, ddq, a, d, theta):
		opL_0 = 1.40579962855621e-65*a[1]**2*sin(theta[1])**2*sin(theta[1] + theta[2] + theta[3])*cos(theta[1] + theta[2] + theta[3])*dq[0]**2 - 2.29584502165847e-49*a[1]**2*sin(theta[1])*sin(theta[2] + theta[3])*sin(theta[1] + theta[2] + theta[3])*dq[0]*dq[1] + 2.29584502165847e-49*a[1]**2*sin(theta[1])*cos(theta[2] + theta[3])*cos(theta[1] + theta[2] + theta[3])*dq[0]*dq[1] - 3.74939945665464e-33*a[1]**2*sin(theta[2] + theta[3])*cos(theta[2] + theta[3])*dq[1]**2 - 2.29584502165847e-49*a[1]*a[2]*sin(theta[1])*sin(theta[3])*sin(theta[1] + theta[2] + theta[3])*dq[0]*dq[1] - 2.29584502165847e-49*a[1]*a[2]*sin(theta[1])*sin(theta[3])*sin(theta[1] + theta[2] + theta[3])*dq[0]*dq[2] + 2.29584502165847e-49*a[1]*a[2]*sin(theta[1])*cos(theta[3])*cos(theta[1] + theta[2] + theta[3])*dq[0]*dq[1] + 2.29584502165847e-49*a[1]*a[2]*sin(theta[1])*cos(theta[3])*cos(theta[1] + theta[2] + theta[3])*dq[0]*dq[2] + 4.31804923432697e-66*a[1]*a[2]*sin(theta[1])*cos(theta[1] + theta[2])*cos(theta[1] + theta[2] + theta[3])*dq[0]*dq[1] + 4.31804923432697e-66*a[1]*a[2]*sin(theta[1])*cos(theta[1] + theta[2])*cos(theta[1] + theta[2] + theta[3])*dq[0]*dq[2] - 3.74939945665464e-33*a[1]*a[2]*sin(theta[3])*cos(theta[2] + theta[3])*dq[1]**2 - 3.74939945665464e-33*a[1]*a[2]*sin(theta[3])*cos(theta[2] + theta[3])*dq[1]*dq[2] - 3.74939945665464e-33*a[1]*a[2]*sin(theta[2] + theta[3])*cos(theta[3])*dq[1]**2 - 3.74939945665464e-33*a[1]*a[2]*sin(theta[2] + theta[3])*cos(theta[3])*dq[1]*dq[2] - 7.05190955846758e-50*a[1]*a[2]*sin(theta[2] + theta[3])*cos(theta[1] + theta[2])*dq[1]**2 - 7.05190955846758e-50*a[1]*a[2]*sin(theta[2] + theta[3])*cos(theta[1] + theta[2])*dq[1]*dq[2] - 3.74939945665464e-33*a[2]**2*sin(theta[3])*cos(theta[3])*dq[1]**2 - 7.49879891330929e-33*a[2]**2*sin(theta[3])*cos(theta[3])*dq[1]*dq[2] - 3.74939945665464e-33*a[2]**2*sin(theta[3])*cos(theta[3])*dq[2]**2 - 7.05190955846758e-50*a[2]**2*sin(theta[3])*cos(theta[1] + theta[2])*dq[1]**2 - 1.41038191169352e-49*a[2]**2*sin(theta[3])*cos(theta[1] + theta[2])*dq[1]*dq[2] - 7.05190955846758e-50*a[2]**2*sin(theta[3])*cos(theta[1] + theta[2])*dq[2]**2
		return opL_0

	def opL1(self, q, dq, ddq, a, d, theta):
		opL_1 = -1.0*a[0]*sin(theta[1] + theta[2] + theta[3])*dq[0]**2 - 6.12323399573677e-17*a[0]*sin(theta[1] + theta[2] + theta[3])*dq[0]*dq[1] - 6.12323399573677e-17*a[0]*sin(theta[1] + theta[2] + theta[3])*dq[0]*dq[2] - 6.12323399573677e-17*a[0]*sin(theta[1] + theta[2] + theta[3])*dq[0]*dq[3] + 6.12323399573677e-17*a[0]*cos(theta[1] + theta[2] + theta[3])*ddq[0] + 1.50949372273656e-48*a[1]*sin(theta[0])*sin(theta[1])*sin(theta[1] + theta[2] + theta[3])*cos(theta[0])*cos(theta[1] + theta[2] + theta[3])*dq[0]*dq[1] - 1.50949372273656e-48*a[1]*sin(theta[0])*sin(theta[1])*sin(theta[1] + theta[2] + theta[3])*cos(theta[0])*dq[0]*dq[1] + 1.50949372273656e-48*a[1]*sin(theta[0])*cos(theta[0])*cos(theta[1])*cos(theta[1] + theta[2] + theta[3])*dq[0]*dq[1] - 1.50949372273656e-48*a[1]*sin(theta[0])*cos(theta[0])*cos(theta[1])*dq[0]*dq[1] + 4.59169004331693e-49*a[1]*sin(theta[1])**2*sin(theta[2] + theta[3])*cos(theta[1] + theta[2] + theta[3])*dq[0]*dq[1] + 1.22464679914735e-16*a[1]*sin(theta[1])**2*sin(theta[2] + theta[3])*dq[0]*dq[1] - 2.82984747192602e-81*a[1]*sin(theta[1])*sin(theta[1] + theta[2] + theta[3])**2*dq[0]*dq[1] - 2.82984747192602e-81*a[1]*sin(theta[1])*sin(theta[1] + theta[2] + theta[3])**2*dq[0]*dq[2] - 2.82984747192602e-81*a[1]*sin(theta[1])*sin(theta[1] + theta[2] + theta[3])**2*dq[0]*dq[3] + 6.12323399573677e-17*a[1]*sin(theta[1])*sin(theta[1] + theta[2] + theta[3])*ddq[0] - 1.22464679914735e-16*a[1]*sin(theta[1])*cos(theta[1])*cos(theta[2] + theta[3])*dq[0]*dq[1] + 3.74939945665464e-33*a[1]*sin(theta[1])*cos(theta[1] + theta[2] + theta[3])*dq[0]**2 + 1.22464679914735e-16*a[1]*sin(theta[1])*cos(theta[1] + theta[2] + theta[3])*dq[0]*dq[1] + 4.59169004331693e-49*a[1]*sin(theta[2] + theta[3])*cos(theta[1])**2*cos(theta[1] + theta[2] + theta[3])*dq[0]*dq[1] - 4.59169004331693e-49*a[1]*sin(theta[2] + theta[3])*cos(theta[1])**2*dq[0]*dq[1] - 2.29584502165847e-49*a[1]*sin(theta[2] + theta[3])*cos(theta[1] + theta[2] + theta[3])*dq[0]*dq[1] - 1.22464679914735e-16*a[1]*sin(theta[2] + theta[3])*dq[0]*dq[1] - 1.0*a[1]*sin(theta[2] + theta[3])*dq[1]**2 - 1.0*a[1]*sin(theta[1] + theta[2] + theta[3])*cos(theta[1])*dq[0]**2 - 6.12323399573677e-17*a[1]*sin(theta[1] + theta[2] + theta[3])*cos(theta[1])*dq[0]*dq[1] - 6.12323399573677e-17*a[1]*sin(theta[1] + theta[2] + theta[3])*cos(theta[1])*dq[0]*dq[2] - 6.12323399573677e-17*a[1]*sin(theta[1] + theta[2] + theta[3])*cos(theta[1])*dq[0]*dq[3] - 4.62149163970586e-65*a[1]*sin(theta[1] + theta[2] + theta[3])*cos(theta[2] + theta[3])*dq[1]**2 - 4.62149163970586e-65*a[1]*sin(theta[1] + theta[2] + theta[3])*cos(theta[2] + theta[3])*dq[1]*dq[2] - 4.62149163970586e-65*a[1]*sin(theta[1] + theta[2] + theta[3])*cos(theta[2] + theta[3])*dq[1]*dq[3] + 6.12323399573677e-17*a[1]*cos(theta[1])*cos(theta[1] + theta[2] + theta[3])*ddq[0] + 1.0*a[1]*cos(theta[2] + theta[3])*ddq[1] - 1.50949372273656e-48*a[2]*sin(theta[0])*sin(theta[1])*sin(theta[2])*sin(theta[3])*sin(theta[1] + theta[2])*cos(theta[0])*cos(theta[1] + theta[2] + theta[3])*dq[0]*dq[1] - 1.50949372273656e-48*a[2]*sin(theta[0])*sin(theta[1])*sin(theta[2])*sin(theta[3])*sin(theta[1] + theta[2])*cos(theta[0])*cos(theta[1] + theta[2] + theta[3])*dq[0]*dq[2] + 1.50949372273656e-48*a[2]*sin(theta[0])*sin(theta[1])*sin(theta[2])*sin(theta[3])*sin(theta[1] + theta[2])*cos(theta[0])*dq[0]*dq[1] + 1.50949372273656e-48*a[2]*sin(theta[0])*sin(theta[1])*sin(theta[2])*sin(theta[3])*sin(theta[1] + theta[2])*cos(theta[0])*dq[0]*dq[2] + 3.01898744547311e-48*a[2]*sin(theta[0])*sin(theta[1])*sin(theta[2])*cos(theta[0])*cos(theta[1])*cos(theta[2])*cos(theta[3])*cos(theta[1] + theta[2] + theta[3])*dq[0]*dq[1] + 3.01898744547311e-48*a[2]*sin(theta[0])*sin(theta[1])*sin(theta[2])*cos(theta[0])*cos(theta[1])*cos(theta[2])*cos(theta[3])*cos(theta[1] + theta[2] + theta[3])*dq[0]*dq[2] - 3.01898744547311e-48*a[2]*sin(theta[0])*sin(theta[1])*sin(theta[2])*cos(theta[0])*cos(theta[1])*cos(theta[2])*cos(theta[3])*dq[0]*dq[1] - 3.01898744547311e-48*a[2]*sin(theta[0])*sin(theta[1])*sin(theta[2])*cos(theta[0])*cos(theta[1])*cos(theta[2])*cos(theta[3])*dq[0]*dq[2] + 1.50949372273656e-48*a[2]*sin(theta[0])*sin(theta[1])*sin(theta[1] + theta[3])*cos(theta[0])*cos(theta[2])**2*cos(theta[1] + theta[2] + theta[3])*dq[0]*dq[1] + 1.50949372273656e-48*a[2]*sin(theta[0])*sin(theta[1])*sin(theta[1] + theta[3])*cos(theta[0])*cos(theta[2])**2*cos(theta[1] + theta[2] + theta[3])*dq[0]*dq[2] - 1.50949372273656e-48*a[2]*sin(theta[0])*sin(theta[1])*sin(theta[1] + theta[3])*cos(theta[0])*cos(theta[2])**2*dq[0]*dq[1] - 1.50949372273656e-48*a[2]*sin(theta[0])*sin(theta[1])*sin(theta[1] + theta[3])*cos(theta[0])*cos(theta[2])**2*dq[0]*dq[2] + 1.50949372273656e-48*a[2]*sin(theta[0])*sin(theta[2])*sin(theta[2] + theta[3])*cos(theta[0])*cos(theta[1])**2*cos(theta[1] + theta[2] + theta[3])*dq[0]*dq[1] + 1.50949372273656e-48*a[2]*sin(theta[0])*sin(theta[2])*sin(theta[2] + theta[3])*cos(theta[0])*cos(theta[1])**2*cos(theta[1] + theta[2] + theta[3])*dq[0]*dq[2] - 1.50949372273656e-48*a[2]*sin(theta[0])*sin(theta[2])*sin(theta[2] + theta[3])*cos(theta[0])*cos(theta[1])**2*dq[0]*dq[1] - 1.50949372273656e-48*a[2]*sin(theta[0])*sin(theta[2])*sin(theta[2] + theta[3])*cos(theta[0])*cos(theta[1])**2*dq[0]*dq[2] - 6.12323399573677e-17*a[2]*sin(theta[3])*dq[0]*dq[1] - 6.12323399573677e-17*a[2]*sin(theta[3])*dq[0]*dq[2] - 1.0*a[2]*sin(theta[3])*dq[1]**2 - 2.0*a[2]*sin(theta[3])*dq[1]*dq[2] - 1.0*a[2]*sin(theta[3])*dq[2]**2 + 6.12323399573677e-17*a[2]*sin(theta[1] + theta[2])*cos(theta[1] + theta[2] + theta[3])*dq[0]*dq[1] + 6.12323399573677e-17*a[2]*sin(theta[1] + theta[2])*cos(theta[1] + theta[2] + theta[3])*dq[0]*dq[2] - 7.54746861368278e-49*a[2]*sin(theta[1] + theta[2])*cos(theta[1] + theta[2] + theta[3])*dq[1]**2 - 1.50949372273656e-48*a[2]*sin(theta[1] + theta[2])*cos(theta[1] + theta[2] + theta[3])*dq[1]*dq[2] - 7.54746861368278e-49*a[2]*sin(theta[1] + theta[2])*cos(theta[1] + theta[2] + theta[3])*dq[2]**2 + 7.05190955846758e-50*a[2]*sin(theta[1] + theta[2])*dq[1]**2 + 1.41038191169352e-49*a[2]*sin(theta[1] + theta[2])*dq[1]*dq[2] + 7.05190955846758e-50*a[2]*sin(theta[1] + theta[2])*dq[2]**2 - 4.62149163970586e-65*a[2]*sin(theta[1] + theta[2] + theta[3])*cos(theta[3])*dq[1]**2 - 9.24298327941172e-65*a[2]*sin(theta[1] + theta[2] + theta[3])*cos(theta[3])*dq[1]*dq[2] - 4.62149163970586e-65*a[2]*sin(theta[1] + theta[2] + theta[3])*cos(theta[3])*dq[1]*dq[3] - 4.62149163970586e-65*a[2]*sin(theta[1] + theta[2] + theta[3])*cos(theta[3])*dq[2]**2 - 4.62149163970586e-65*a[2]*sin(theta[1] + theta[2] + theta[3])*cos(theta[3])*dq[2]*dq[3] - 1.0*a[2]*sin(theta[1] + theta[2] + theta[3])*cos(theta[1] + theta[2])*dq[0]**2 - 6.12323399573676e-17*a[2]*sin(theta[1] + theta[2] + theta[3])*cos(theta[1] + theta[2])*dq[0]*dq[1] - 6.12323399573676e-17*a[2]*sin(theta[1] + theta[2] + theta[3])*cos(theta[1] + theta[2])*dq[0]*dq[2] - 6.12323399573677e-17*a[2]*sin(theta[1] + theta[2] + theta[3])*cos(theta[1] + theta[2])*dq[0]*dq[3] + 7.54746861368278e-49*a[2]*sin(theta[1] + theta[2] + theta[3])*cos(theta[1] + theta[2])*dq[1]**2 + 1.50949372273656e-48*a[2]*sin(theta[1] + theta[2] + theta[3])*cos(theta[1] + theta[2])*dq[1]*dq[2] + 7.54746861368278e-49*a[2]*sin(theta[1] + theta[2] + theta[3])*cos(theta[1] + theta[2])*dq[1]*dq[3] + 7.54746861368278e-49*a[2]*sin(theta[1] + theta[2] + theta[3])*cos(theta[1] + theta[2])*dq[2]**2 + 7.54746861368278e-49*a[2]*sin(theta[1] + theta[2] + theta[3])*cos(theta[1] + theta[2])*dq[2]*dq[3] + 1.0*a[2]*cos(theta[3])*ddq[1] + 1.0*a[2]*cos(theta[3])*ddq[2] + 6.12323399573677e-17*a[2]*cos(theta[1] + theta[2])*cos(theta[1] + theta[2] + theta[3])*ddq[0] - 7.54746861368278e-49*a[2]*cos(theta[1] + theta[2])*cos(theta[1] + theta[2] + theta[3])*ddq[1] - 7.54746861368278e-49*a[2]*cos(theta[1] + theta[2])*cos(theta[1] + theta[2] + theta[3])*ddq[2] + 7.05190955846758e-50*a[2]*cos(theta[1] + theta[2])*ddq[1] + 7.05190955846758e-50*a[2]*cos(theta[1] + theta[2])*ddq[2] + 9.82*sin(theta[3])*sin(theta[1] + theta[2]) - 9.82*cos(theta[3])*cos(theta[1] + theta[2])
		return opL_1

	def opL2(self, q, dq, ddq, a, d, theta):
		opL_2 = -3.74939945665464e-33*a[0]*sin(theta[1] + theta[2] + theta[3])*ddq[0] - 6.12323399573677e-17*a[0]*cos(theta[1] + theta[2] + theta[3])*dq[0]**2 - 3.74939945665464e-33*a[0]*cos(theta[1] + theta[2] + theta[3])*dq[0]*dq[1] - 3.74939945665464e-33*a[0]*cos(theta[1] + theta[2] + theta[3])*dq[0]*dq[2] - 3.74939945665464e-33*a[0]*cos(theta[1] + theta[2] + theta[3])*dq[0]*dq[3] - 7.49879891330929e-33*a[1]*sin(theta[1])**2*cos(theta[2] + theta[3])*cos(theta[1] + theta[2] + theta[3])*dq[0]*dq[1] + 7.49879891330929e-33*a[1]*sin(theta[1])**2*cos(theta[2] + theta[3])*dq[0]*dq[1] - 7.49879891330929e-33*a[1]*sin(theta[1])*sin(theta[2] + theta[3])*cos(theta[1])*cos(theta[1] + theta[2] + theta[3])*dq[0]*dq[1] + 7.49879891330929e-33*a[1]*sin(theta[1])*sin(theta[2] + theta[3])*cos(theta[1])*dq[0]*dq[1] - 3.74939945665464e-33*a[1]*sin(theta[1])*sin(theta[1] + theta[2] + theta[3])**2*ddq[0] + 7.49879891330929e-33*a[1]*sin(theta[1])*sin(theta[1] + theta[2] + theta[3])*cos(theta[1] + theta[2] + theta[3])*dq[0]*dq[1] - 2.29584502165847e-49*a[1]*sin(theta[1])*sin(theta[1] + theta[2] + theta[3])*dq[0]**2 - 7.49879891330929e-33*a[1]*sin(theta[1])*sin(theta[1] + theta[2] + theta[3])*dq[0]*dq[1] - 3.74939945665464e-33*a[1]*sin(theta[1])*cos(theta[1] + theta[2] + theta[3])**2*ddq[0] + 3.74939945665464e-33*a[1]*sin(theta[1])*cos(theta[1] + theta[2] + theta[3])*ddq[0] + 7.54746861368278e-49*a[1]*sin(theta[2] + theta[3])*sin(theta[1] + theta[2] + theta[3])*dq[1]**2 - 6.12323399573677e-17*a[1]*sin(theta[2] + theta[3])*sin(theta[1] + theta[2] + theta[3])*dq[1]*dq[2] - 6.12323399573677e-17*a[1]*sin(theta[2] + theta[3])*sin(theta[1] + theta[2] + theta[3])*dq[1]*dq[3] + 6.12323399573677e-17*a[1]*sin(theta[2] + theta[3])*cos(theta[1] + theta[2] + theta[3])*ddq[1] - 6.12323399573677e-17*a[1]*sin(theta[2] + theta[3])*ddq[1] - 3.74939945665464e-33*a[1]*sin(theta[1] + theta[2] + theta[3])*cos(theta[1])*ddq[0] - 6.12323399573677e-17*a[1]*sin(theta[1] + theta[2] + theta[3])*cos(theta[2] + theta[3])*ddq[1] - 6.12323399573677e-17*a[1]*cos(theta[1])*cos(theta[1] + theta[2] + theta[3])*dq[0]**2 - 3.74939945665464e-33*a[1]*cos(theta[1])*cos(theta[1] + theta[2] + theta[3])*dq[0]*dq[1] - 3.74939945665464e-33*a[1]*cos(theta[1])*cos(theta[1] + theta[2] + theta[3])*dq[0]*dq[2] - 3.74939945665464e-33*a[1]*cos(theta[1])*cos(theta[1] + theta[2] + theta[3])*dq[0]*dq[3] + 7.54746861368278e-49*a[1]*cos(theta[2] + theta[3])*cos(theta[1] + theta[2] + theta[3])*dq[1]**2 - 6.12323399573677e-17*a[1]*cos(theta[2] + theta[3])*cos(theta[1] + theta[2] + theta[3])*dq[1]*dq[2] - 6.12323399573677e-17*a[1]*cos(theta[2] + theta[3])*cos(theta[1] + theta[2] + theta[3])*dq[1]*dq[3] - 7.49879891330929e-33*a[1]*cos(theta[2] + theta[3])*dq[0]*dq[1] - 6.12323399573677e-17*a[1]*cos(theta[2] + theta[3])*dq[1]**2 + 7.54746861368278e-49*a[2]*sin(theta[3])*sin(theta[1] + theta[2] + theta[3])*dq[1]**2 + 1.50949372273656e-48*a[2]*sin(theta[3])*sin(theta[1] + theta[2] + theta[3])*dq[1]*dq[2] - 6.12323399573677e-17*a[2]*sin(theta[3])*sin(theta[1] + theta[2] + theta[3])*dq[1]*dq[3] + 7.54746861368278e-49*a[2]*sin(theta[3])*sin(theta[1] + theta[2] + theta[3])*dq[2]**2 - 6.12323399573677e-17*a[2]*sin(theta[3])*sin(theta[1] + theta[2] + theta[3])*dq[2]*dq[3] + 6.12323399573677e-17*a[2]*sin(theta[3])*cos(theta[1] + theta[2] + theta[3])*ddq[1] + 6.12323399573677e-17*a[2]*sin(theta[3])*cos(theta[1] + theta[2] + theta[3])*ddq[2] - 6.12323399573677e-17*a[2]*sin(theta[3])*ddq[1] - 6.12323399573677e-17*a[2]*sin(theta[3])*ddq[2] - 3.74939945665464e-33*a[2]*sin(theta[1] + theta[2])*sin(theta[1] + theta[2] + theta[3])*dq[0]*dq[1] - 3.74939945665464e-33*a[2]*sin(theta[1] + theta[2])*sin(theta[1] + theta[2] + theta[3])*dq[0]*dq[2] + 4.18968671627317e-65*a[2]*sin(theta[1] + theta[2])*sin(theta[1] + theta[2] + theta[3])*dq[1]**2 + 8.37937343254633e-65*a[2]*sin(theta[1] + theta[2])*sin(theta[1] + theta[2] + theta[3])*dq[1]*dq[2] + 4.18968671627317e-65*a[2]*sin(theta[1] + theta[2])*sin(theta[1] + theta[2] + theta[3])*dq[2]**2 - 6.12323399573677e-17*a[2]*sin(theta[1] + theta[2] + theta[3])*cos(theta[3])*ddq[1] - 6.12323399573677e-17*a[2]*sin(theta[1] + theta[2] + theta[3])*cos(theta[3])*ddq[2] - 3.74939945665464e-33*a[2]*sin(theta[1] + theta[2] + theta[3])*cos(theta[1] + theta[2])*ddq[0] + 4.18968671627317e-65*a[2]*sin(theta[1] + theta[2] + theta[3])*cos(theta[1] + theta[2])*ddq[1] + 4.18968671627317e-65*a[2]*sin(theta[1] + theta[2] + theta[3])*cos(theta[1] + theta[2])*ddq[2] + 7.54746861368278e-49*a[2]*cos(theta[3])*cos(theta[1] + theta[2] + theta[3])*dq[1]**2 + 1.50949372273656e-48*a[2]*cos(theta[3])*cos(theta[1] + theta[2] + theta[3])*dq[1]*dq[2] - 6.12323399573677e-17*a[2]*cos(theta[3])*cos(theta[1] + theta[2] + theta[3])*dq[1]*dq[3] + 7.54746861368278e-49*a[2]*cos(theta[3])*cos(theta[1] + theta[2] + theta[3])*dq[2]**2 - 6.12323399573677e-17*a[2]*cos(theta[3])*cos(theta[1] + theta[2] + theta[3])*dq[2]*dq[3] - 3.74939945665464e-33*a[2]*cos(theta[3])*dq[0]*dq[1] - 3.74939945665464e-33*a[2]*cos(theta[3])*dq[0]*dq[2] - 6.12323399573677e-17*a[2]*cos(theta[3])*dq[1]**2 - 1.22464679914735e-16*a[2]*cos(theta[3])*dq[1]*dq[2] - 6.12323399573677e-17*a[2]*cos(theta[3])*dq[2]**2 - 6.12323399573677e-17*a[2]*cos(theta[1] + theta[2])*cos(theta[1] + theta[2] + theta[3])*dq[0]**2 - 3.74939945665464e-33*a[2]*cos(theta[1] + theta[2])*cos(theta[1] + theta[2] + theta[3])*dq[0]*dq[1] - 3.74939945665464e-33*a[2]*cos(theta[1] + theta[2])*cos(theta[1] + theta[2] + theta[3])*dq[0]*dq[2] - 3.74939945665464e-33*a[2]*cos(theta[1] + theta[2])*cos(theta[1] + theta[2] + theta[3])*dq[0]*dq[3] + 4.18968671627317e-65*a[2]*cos(theta[1] + theta[2])*cos(theta[1] + theta[2] + theta[3])*dq[1]**2 + 8.37937343254633e-65*a[2]*cos(theta[1] + theta[2])*cos(theta[1] + theta[2] + theta[3])*dq[1]*dq[2] + 4.18968671627317e-65*a[2]*cos(theta[1] + theta[2])*cos(theta[1] + theta[2] + theta[3])*dq[1]*dq[3] + 4.18968671627317e-65*a[2]*cos(theta[1] + theta[2])*cos(theta[1] + theta[2] + theta[3])*dq[2]**2 + 4.18968671627317e-65*a[2]*cos(theta[1] + theta[2])*cos(theta[1] + theta[2] + theta[3])*dq[2]*dq[3] + 6.0130157838135e-16*sin(theta[3])*cos(theta[1] + theta[2]) + 6.0130157838135e-16*sin(theta[1] + theta[2])*cos(theta[3])
		return opL_2

	def opL3(self, q, dq, ddq, a, d, theta):
		opL_3 = 6.12323399573677e-17*a[0]*sin(theta[1] + theta[2] + theta[3])*ddq[0] + 1.0*a[0]*cos(theta[1] + theta[2] + theta[3])*dq[0]**2 + 6.12323399573677e-17*a[0]*cos(theta[1] + theta[2] + theta[3])*dq[0]*dq[1] + 6.12323399573677e-17*a[0]*cos(theta[1] + theta[2] + theta[3])*dq[0]*dq[2] + 6.12323399573677e-17*a[0]*cos(theta[1] + theta[2] + theta[3])*dq[0]*dq[3] + 1.50949372273656e-48*a[1]*sin(theta[0])*sin(theta[1])*sin(theta[1] + theta[2] + theta[3])**2*cos(theta[0])*dq[0]*dq[1] + 1.50949372273656e-48*a[1]*sin(theta[0])*sin(theta[1] + theta[2] + theta[3])*cos(theta[0])*cos(theta[1])*dq[0]*dq[1] - 4.59169004331693e-49*a[1]*sin(theta[1])**2*cos(theta[2] + theta[3])*cos(theta[1] + theta[2] + theta[3])*dq[0]*dq[1] - 1.22464679914735e-16*a[1]*sin(theta[1])**2*cos(theta[2] + theta[3])*dq[0]*dq[1] - 4.59169004331693e-49*a[1]*sin(theta[1])*sin(theta[2] + theta[3])*cos(theta[1])*cos(theta[1] + theta[2] + theta[3])*dq[0]*dq[1] - 1.22464679914735e-16*a[1]*sin(theta[1])*sin(theta[2] + theta[3])*cos(theta[1])*dq[0]*dq[1] - 2.29584502165847e-49*a[1]*sin(theta[1])*sin(theta[1] + theta[2] + theta[3])**2*ddq[0] + 4.59169004331694e-49*a[1]*sin(theta[1])*sin(theta[1] + theta[2] + theta[3])*cos(theta[1] + theta[2] + theta[3])*dq[0]*dq[1] + 2.82984747192602e-81*a[1]*sin(theta[1])*sin(theta[1] + theta[2] + theta[3])*cos(theta[1] + theta[2] + theta[3])*dq[0]*dq[2] + 2.82984747192602e-81*a[1]*sin(theta[1])*sin(theta[1] + theta[2] + theta[3])*cos(theta[1] + theta[2] + theta[3])*dq[0]*dq[3] + 3.74939945665464e-33*a[1]*sin(theta[1])*sin(theta[1] + theta[2] + theta[3])*dq[0]**2 + 1.22464679914735e-16*a[1]*sin(theta[1])*sin(theta[1] + theta[2] + theta[3])*dq[0]*dq[1] - 2.29584502165847e-49*a[1]*sin(theta[1])*cos(theta[1] + theta[2] + theta[3])**2*ddq[0] - 6.12323399573677e-17*a[1]*sin(theta[1])*cos(theta[1] + theta[2] + theta[3])*ddq[0] + 2.05235457060565e-80*a[1]*sin(theta[2] + theta[3])*sin(theta[1] + theta[2] + theta[3])*dq[1]**2 - 3.74939945665464e-33*a[1]*sin(theta[2] + theta[3])*sin(theta[1] + theta[2] + theta[3])*dq[1]*dq[2] - 3.74939945665464e-33*a[1]*sin(theta[2] + theta[3])*sin(theta[1] + theta[2] + theta[3])*dq[1]*dq[3] + 3.74939945665464e-33*a[1]*sin(theta[2] + theta[3])*cos(theta[1] + theta[2] + theta[3])*ddq[1] + 1.0*a[1]*sin(theta[2] + theta[3])*ddq[1] + 6.12323399573677e-17*a[1]*sin(theta[1] + theta[2] + theta[3])*cos(theta[1])*ddq[0] - 3.74939945665464e-33*a[1]*sin(theta[1] + theta[2] + theta[3])*cos(theta[2] + theta[3])*ddq[1] + 1.0*a[1]*cos(theta[1])*cos(theta[1] + theta[2] + theta[3])*dq[0]**2 + 6.12323399573677e-17*a[1]*cos(theta[1])*cos(theta[1] + theta[2] + theta[3])*dq[0]*dq[1] + 6.12323399573677e-17*a[1]*cos(theta[1])*cos(theta[1] + theta[2] + theta[3])*dq[0]*dq[2] + 6.12323399573677e-17*a[1]*cos(theta[1])*cos(theta[1] + theta[2] + theta[3])*dq[0]*dq[3] + 4.62149163970586e-65*a[1]*cos(theta[2] + theta[3])*cos(theta[1] + theta[2] + theta[3])*dq[1]**2 - 3.74939945665464e-33*a[1]*cos(theta[2] + theta[3])*cos(theta[1] + theta[2] + theta[3])*dq[1]*dq[2] - 3.74939945665464e-33*a[1]*cos(theta[2] + theta[3])*cos(theta[1] + theta[2] + theta[3])*dq[1]*dq[3] + 1.22464679914735e-16*a[1]*cos(theta[2] + theta[3])*dq[0]*dq[1] + 1.0*a[1]*cos(theta[2] + theta[3])*dq[1]**2 + 3.01898744547311e-48*a[2]*sin(theta[0])*sin(theta[1])*sin(theta[2])*sin(theta[1] + theta[2] + theta[3])*cos(theta[0])*cos(theta[1] + theta[2] + theta[3])*dq[0]*dq[1] + 3.01898744547311e-48*a[2]*sin(theta[0])*sin(theta[1])*sin(theta[2])*sin(theta[1] + theta[2] + theta[3])*cos(theta[0])*cos(theta[1] + theta[2] + theta[3])*dq[0]*dq[2] + 1.50949372273656e-48*a[2]*sin(theta[0])*sin(theta[1])*sin(theta[1] + theta[3])*sin(theta[1] + theta[2] + theta[3])*cos(theta[0])*dq[0]*dq[1] + 1.50949372273656e-48*a[2]*sin(theta[0])*sin(theta[1])*sin(theta[1] + theta[3])*sin(theta[1] + theta[2] + theta[3])*cos(theta[0])*dq[0]*dq[2] + 1.50949372273656e-48*a[2]*sin(theta[0])*sin(theta[2])*sin(theta[2] + theta[3])*sin(theta[1] + theta[2] + theta[3])*cos(theta[0])*dq[0]*dq[1] + 1.50949372273656e-48*a[2]*sin(theta[0])*sin(theta[2])*sin(theta[2] + theta[3])*sin(theta[1] + theta[2] + theta[3])*cos(theta[0])*dq[0]*dq[2] + 2.05235457060565e-80*a[2]*sin(theta[3])*sin(theta[1] + theta[2] + theta[3])*dq[1]**2 + 4.10470914121129e-80*a[2]*sin(theta[3])*sin(theta[1] + theta[2] + theta[3])*dq[1]*dq[2] - 3.74939945665464e-33*a[2]*sin(theta[3])*sin(theta[1] + theta[2] + theta[3])*dq[1]*dq[3] + 2.05235457060565e-80*a[2]*sin(theta[3])*sin(theta[1] + theta[2] + theta[3])*dq[2]**2 - 3.74939945665464e-33*a[2]*sin(theta[3])*sin(theta[1] + theta[2] + theta[3])*dq[2]*dq[3] + 3.74939945665464e-33*a[2]*sin(theta[3])*cos(theta[1] + theta[2] + theta[3])*ddq[1] + 3.74939945665464e-33*a[2]*sin(theta[3])*cos(theta[1] + theta[2] + theta[3])*ddq[2] + 1.0*a[2]*sin(theta[3])*ddq[1] + 1.0*a[2]*sin(theta[3])*ddq[2] + 6.12323399573677e-17*a[2]*sin(theta[1] + theta[2])*sin(theta[1] + theta[2] + theta[3])*dq[0]*dq[1] + 6.12323399573677e-17*a[2]*sin(theta[1] + theta[2])*sin(theta[1] + theta[2] + theta[3])*dq[0]*dq[2] - 7.54746861368278e-49*a[2]*sin(theta[1] + theta[2])*sin(theta[1] + theta[2] + theta[3])*dq[1]**2 - 1.50949372273656e-48*a[2]*sin(theta[1] + theta[2])*sin(theta[1] + theta[2] + theta[3])*dq[1]*dq[2] - 7.54746861368278e-49*a[2]*sin(theta[1] + theta[2])*sin(theta[1] + theta[2] + theta[3])*dq[2]**2 - 3.74939945665464e-33*a[2]*sin(theta[1] + theta[2] + theta[3])*cos(theta[3])*ddq[1] - 3.74939945665464e-33*a[2]*sin(theta[1] + theta[2] + theta[3])*cos(theta[3])*ddq[2] + 6.12323399573677e-17*a[2]*sin(theta[1] + theta[2] + theta[3])*cos(theta[1] + theta[2])*ddq[0] - 7.54746861368278e-49*a[2]*sin(theta[1] + theta[2] + theta[3])*cos(theta[1] + theta[2])*ddq[1] - 7.54746861368278e-49*a[2]*sin(theta[1] + theta[2] + theta[3])*cos(theta[1] + theta[2])*ddq[2] + 4.62149163970586e-65*a[2]*cos(theta[3])*cos(theta[1] + theta[2] + theta[3])*dq[1]**2 + 9.24298327941173e-65*a[2]*cos(theta[3])*cos(theta[1] + theta[2] + theta[3])*dq[1]*dq[2] - 3.74939945665464e-33*a[2]*cos(theta[3])*cos(theta[1] + theta[2] + theta[3])*dq[1]*dq[3] + 4.62149163970586e-65*a[2]*cos(theta[3])*cos(theta[1] + theta[2] + theta[3])*dq[2]**2 - 3.74939945665464e-33*a[2]*cos(theta[3])*cos(theta[1] + theta[2] + theta[3])*dq[2]*dq[3] + 6.12323399573677e-17*a[2]*cos(theta[3])*dq[0]*dq[1] + 6.12323399573677e-17*a[2]*cos(theta[3])*dq[0]*dq[2] + 1.0*a[2]*cos(theta[3])*dq[1]**2 + 2.0*a[2]*cos(theta[3])*dq[1]*dq[2] + 1.0*a[2]*cos(theta[3])*dq[2]**2 + 1.0*a[2]*cos(theta[1] + theta[2])*cos(theta[1] + theta[2] + theta[3])*dq[0]**2 + 6.12323399573676e-17*a[2]*cos(theta[1] + theta[2])*cos(theta[1] + theta[2] + theta[3])*dq[0]*dq[1] + 6.12323399573676e-17*a[2]*cos(theta[1] + theta[2])*cos(theta[1] + theta[2] + theta[3])*dq[0]*dq[2] + 6.12323399573677e-17*a[2]*cos(theta[1] + theta[2])*cos(theta[1] + theta[2] + theta[3])*dq[0]*dq[3] - 7.54746861368278e-49*a[2]*cos(theta[1] + theta[2])*cos(theta[1] + theta[2] + theta[3])*dq[1]**2 - 1.50949372273656e-48*a[2]*cos(theta[1] + theta[2])*cos(theta[1] + theta[2] + theta[3])*dq[1]*dq[2] - 7.54746861368278e-49*a[2]*cos(theta[1] + theta[2])*cos(theta[1] + theta[2] + theta[3])*dq[1]*dq[3] - 7.54746861368278e-49*a[2]*cos(theta[1] + theta[2])*cos(theta[1] + theta[2] + theta[3])*dq[2]**2 - 7.54746861368278e-49*a[2]*cos(theta[1] + theta[2])*cos(theta[1] + theta[2] + theta[3])*dq[2]*dq[3] - 9.82*sin(theta[3])*cos(theta[1] + theta[2]) - 9.82*sin(theta[1] + theta[2])*cos(theta[3])
		return opL_3

	def opL4(self, q, dq, ddq, a, d, theta):
		opL_4 = 0.5*sin(2*theta[1] + 2*theta[2] + 2*theta[3])*dq[0]**2 + 3.06161699786838e-17*sin(2*theta[1] + 2*theta[2] + 2*theta[3])*dq[0]*dq[1] + 3.06161699786838e-17*sin(2*theta[1] + 2*theta[2] + 2*theta[3])*dq[0]*dq[2] + 3.06161699786838e-17*sin(2*theta[1] + 2*theta[2] + 2*theta[3])*dq[0]*dq[3] - 4.62149163970586e-65*sin(2*theta[1] + 2*theta[2] + 2*theta[3])*dq[1]**2 - 9.24298327941173e-65*sin(2*theta[1] + 2*theta[2] + 2*theta[3])*dq[1]*dq[2] - 9.24298327941173e-65*sin(2*theta[1] + 2*theta[2] + 2*theta[3])*dq[1]*dq[3] - 4.62149163970586e-65*sin(2*theta[1] + 2*theta[2] + 2*theta[3])*dq[2]**2 - 9.24298327941173e-65*sin(2*theta[1] + 2*theta[2] + 2*theta[3])*dq[2]*dq[3] - 4.62149163970586e-65*sin(2*theta[1] + 2*theta[2] + 2*theta[3])*dq[3]**2 - 3.06161699786838e-17*cos(2*theta[1] + 2*theta[2] + 2*theta[3])*ddq[0] - 1.87469972832732e-33*cos(2*theta[1] + 2*theta[2] + 2*theta[3])*ddq[1] - 1.87469972832732e-33*cos(2*theta[1] + 2*theta[2] + 2*theta[3])*ddq[2] - 1.87469972832732e-33*cos(2*theta[1] + 2*theta[2] + 2*theta[3])*ddq[3] + 3.06161699786838e-17*ddq[0] + 1.87469972832732e-33*ddq[1] + 1.87469972832732e-33*ddq[2] + 1.87469972832732e-33*ddq[3]
		return opL_4

	def opL5(self, q, dq, ddq, a, d, theta):
		opL_5 = -3.74939945665464e-33*sin(theta[1] + theta[2] + theta[3])*cos(theta[1] + theta[2] + theta[3])*dq[0]**2 - 2.29584502165847e-49*sin(theta[1] + theta[2] + theta[3])*cos(theta[1] + theta[2] + theta[3])*dq[0]*dq[1] - 2.29584502165847e-49*sin(theta[1] + theta[2] + theta[3])*cos(theta[1] + theta[2] + theta[3])*dq[0]*dq[2] - 2.29584502165847e-49*sin(theta[1] + theta[2] + theta[3])*cos(theta[1] + theta[2] + theta[3])*dq[0]*dq[3] + 1.5390194223783e-112*sin(theta[1] + theta[2] + theta[3])*cos(theta[1] + theta[2] + theta[3])*dq[1]**2 + 3.0780388447566e-112*sin(theta[1] + theta[2] + theta[3])*cos(theta[1] + theta[2] + theta[3])*dq[1]*dq[2] + 3.0780388447566e-112*sin(theta[1] + theta[2] + theta[3])*cos(theta[1] + theta[2] + theta[3])*dq[1]*dq[3] + 1.5390194223783e-112*sin(theta[1] + theta[2] + theta[3])*cos(theta[1] + theta[2] + theta[3])*dq[2]**2 + 3.0780388447566e-112*sin(theta[1] + theta[2] + theta[3])*cos(theta[1] + theta[2] + theta[3])*dq[2]*dq[3] + 1.5390194223783e-112*sin(theta[1] + theta[2] + theta[3])*cos(theta[1] + theta[2] + theta[3])*dq[3]**2 - 3.74939945665464e-33*sin(theta[1] + theta[2] + theta[3])*dq[0]**2 + 2.51340945560765e-96*sin(theta[1] + theta[2] + theta[3])*dq[0]*dq[1] + 2.51340945560765e-96*sin(theta[1] + theta[2] + theta[3])*dq[0]*dq[2] + 2.51340945560765e-96*sin(theta[1] + theta[2] + theta[3])*dq[0]*dq[3] + 4.10470914121129e-80*sin(theta[1] + theta[2] + theta[3])*dq[1]**2 + 8.20941828242259e-80*sin(theta[1] + theta[2] + theta[3])*dq[1]*dq[2] + 8.20941828242259e-80*sin(theta[1] + theta[2] + theta[3])*dq[1]*dq[3] + 4.10470914121129e-80*sin(theta[1] + theta[2] + theta[3])*dq[2]**2 + 8.20941828242259e-80*sin(theta[1] + theta[2] + theta[3])*dq[2]*dq[3] + 4.10470914121129e-80*sin(theta[1] + theta[2] + theta[3])*dq[3]**2 + 2.29584502165847e-49*cos(theta[1] + theta[2] + theta[3])**2*ddq[0] + 1.40579962855621e-65*cos(theta[1] + theta[2] + theta[3])**2*ddq[1] + 1.40579962855621e-65*cos(theta[1] + theta[2] + theta[3])**2*ddq[2] + 1.40579962855621e-65*cos(theta[1] + theta[2] + theta[3])**2*ddq[3] + 6.12323399573677e-17*cos(theta[1] + theta[2] + theta[3])*ddq[0] + 7.49879891330929e-33*cos(theta[1] + theta[2] + theta[3])*ddq[1] + 7.49879891330929e-33*cos(theta[1] + theta[2] + theta[3])*ddq[2] + 7.49879891330929e-33*cos(theta[1] + theta[2] + theta[3])*ddq[3] + 6.12323399573677e-17*ddq[0] + 1.0*ddq[1] + 1.0*ddq[2] + 1.0*ddq[3]
		return opL_5

	def opL6(self, q, dq, ddq, a, d, theta):
		opL_6 = -1.0*sin(theta[1] + theta[2] + theta[3])*cos(theta[1] + theta[2] + theta[3])*dq[0]**2 - 6.12323399573677e-17*sin(theta[1] + theta[2] + theta[3])*cos(theta[1] + theta[2] + theta[3])*dq[0]*dq[1] - 6.12323399573677e-17*sin(theta[1] + theta[2] + theta[3])*cos(theta[1] + theta[2] + theta[3])*dq[0]*dq[2] - 6.12323399573677e-17*sin(theta[1] + theta[2] + theta[3])*cos(theta[1] + theta[2] + theta[3])*dq[0]*dq[3] + 9.24298327941173e-65*sin(theta[1] + theta[2] + theta[3])*cos(theta[1] + theta[2] + theta[3])*dq[1]**2 + 1.84859665588235e-64*sin(theta[1] + theta[2] + theta[3])*cos(theta[1] + theta[2] + theta[3])*dq[1]*dq[2] + 1.84859665588235e-64*sin(theta[1] + theta[2] + theta[3])*cos(theta[1] + theta[2] + theta[3])*dq[1]*dq[3] + 9.24298327941173e-65*sin(theta[1] + theta[2] + theta[3])*cos(theta[1] + theta[2] + theta[3])*dq[2]**2 + 1.84859665588235e-64*sin(theta[1] + theta[2] + theta[3])*cos(theta[1] + theta[2] + theta[3])*dq[2]*dq[3] + 9.24298327941173e-65*sin(theta[1] + theta[2] + theta[3])*cos(theta[1] + theta[2] + theta[3])*dq[3]**2 + 3.74939945665464e-33*sin(theta[1] + theta[2] + theta[3])*dq[0]**2 - 5.65969494385204e-81*sin(theta[1] + theta[2] + theta[3])*dq[0]*dq[1] - 5.65969494385204e-81*sin(theta[1] + theta[2] + theta[3])*dq[0]*dq[2] - 5.65969494385204e-81*sin(theta[1] + theta[2] + theta[3])*dq[0]*dq[3] - 9.24298327941173e-65*sin(theta[1] + theta[2] + theta[3])*dq[1]**2 - 1.84859665588235e-64*sin(theta[1] + theta[2] + theta[3])*dq[1]*dq[2] - 1.84859665588235e-64*sin(theta[1] + theta[2] + theta[3])*dq[1]*dq[3] - 9.24298327941173e-65*sin(theta[1] + theta[2] + theta[3])*dq[2]**2 - 1.84859665588235e-64*sin(theta[1] + theta[2] + theta[3])*dq[2]*dq[3] - 9.24298327941173e-65*sin(theta[1] + theta[2] + theta[3])*dq[3]**2 + 6.12323399573677e-17*cos(theta[1] + theta[2] + theta[3])**2*ddq[0] + 3.74939945665464e-33*cos(theta[1] + theta[2] + theta[3])**2*ddq[1] + 3.74939945665464e-33*cos(theta[1] + theta[2] + theta[3])**2*ddq[2] + 3.74939945665464e-33*cos(theta[1] + theta[2] + theta[3])**2*ddq[3] - 6.12323399573677e-17*cos(theta[1] + theta[2] + theta[3])*ddq[0] - 7.49879891330929e-33*cos(theta[1] + theta[2] + theta[3])*ddq[1] - 7.49879891330929e-33*cos(theta[1] + theta[2] + theta[3])*ddq[2] - 7.49879891330929e-33*cos(theta[1] + theta[2] + theta[3])*ddq[3] + 2.29584502165847e-49*ddq[0] + 3.74939945665464e-33*ddq[1] + 3.74939945665464e-33*ddq[2] + 3.74939945665464e-33*ddq[3]
		return opL_6

	def opL7(self, q, dq, ddq, a, d, theta):
		opL_7 = -6.12323399573677e-17*sin(theta[1] + theta[2] + theta[3])**2*dq[0]**2 - 3.74939945665464e-33*sin(theta[1] + theta[2] + theta[3])**2*dq[0]*dq[1] - 3.74939945665464e-33*sin(theta[1] + theta[2] + theta[3])**2*dq[0]*dq[2] - 3.74939945665464e-33*sin(theta[1] + theta[2] + theta[3])**2*dq[0]*dq[3] + 1.25670472780382e-96*sin(theta[1] + theta[2] + theta[3])**2*dq[1]**2 + 2.51340945560765e-96*sin(theta[1] + theta[2] + theta[3])**2*dq[1]*dq[2] + 2.51340945560765e-96*sin(theta[1] + theta[2] + theta[3])**2*dq[1]*dq[3] + 1.25670472780382e-96*sin(theta[1] + theta[2] + theta[3])**2*dq[2]**2 + 2.51340945560765e-96*sin(theta[1] + theta[2] + theta[3])**2*dq[2]*dq[3] + 1.25670472780382e-96*sin(theta[1] + theta[2] + theta[3])**2*dq[3]**2 + 7.49879891330929e-33*sin(theta[1] + theta[2] + theta[3])*cos(theta[1] + theta[2] + theta[3])*ddq[0] + 4.59169004331693e-49*sin(theta[1] + theta[2] + theta[3])*cos(theta[1] + theta[2] + theta[3])*ddq[1] + 4.59169004331693e-49*sin(theta[1] + theta[2] + theta[3])*cos(theta[1] + theta[2] + theta[3])*ddq[2] + 4.59169004331693e-49*sin(theta[1] + theta[2] + theta[3])*cos(theta[1] + theta[2] + theta[3])*ddq[3] + 1.0*sin(theta[1] + theta[2] + theta[3])*ddq[0] + 1.22464679914735e-16*sin(theta[1] + theta[2] + theta[3])*ddq[1] + 1.22464679914735e-16*sin(theta[1] + theta[2] + theta[3])*ddq[2] + 1.22464679914735e-16*sin(theta[1] + theta[2] + theta[3])*ddq[3] + 6.12323399573677e-17*cos(theta[1] + theta[2] + theta[3])**2*dq[0]**2 + 3.74939945665464e-33*cos(theta[1] + theta[2] + theta[3])**2*dq[0]*dq[1] + 3.74939945665464e-33*cos(theta[1] + theta[2] + theta[3])**2*dq[0]*dq[2] + 3.74939945665464e-33*cos(theta[1] + theta[2] + theta[3])**2*dq[0]*dq[3] - 2.82984747192602e-81*cos(theta[1] + theta[2] + theta[3])**2*dq[1]**2 - 5.65969494385204e-81*cos(theta[1] + theta[2] + theta[3])**2*dq[1]*dq[2] - 5.65969494385204e-81*cos(theta[1] + theta[2] + theta[3])**2*dq[1]*dq[3] - 2.82984747192602e-81*cos(theta[1] + theta[2] + theta[3])**2*dq[2]**2 - 5.65969494385204e-81*cos(theta[1] + theta[2] + theta[3])**2*dq[2]*dq[3] - 2.82984747192602e-81*cos(theta[1] + theta[2] + theta[3])**2*dq[3]**2 + 6.12323399573677e-17*cos(theta[1] + theta[2] + theta[3])*dq[0]**2 - 4.62149163970586e-65*cos(theta[1] + theta[2] + theta[3])*dq[0]*dq[1] - 4.62149163970586e-65*cos(theta[1] + theta[2] + theta[3])*dq[0]*dq[2] - 4.62149163970586e-65*cos(theta[1] + theta[2] + theta[3])*dq[0]*dq[3] - 7.54746861368278e-49*cos(theta[1] + theta[2] + theta[3])*dq[1]**2 - 1.50949372273656e-48*cos(theta[1] + theta[2] + theta[3])*dq[1]*dq[2] - 1.50949372273656e-48*cos(theta[1] + theta[2] + theta[3])*dq[1]*dq[3] - 7.54746861368278e-49*cos(theta[1] + theta[2] + theta[3])*dq[2]**2 - 1.50949372273656e-48*cos(theta[1] + theta[2] + theta[3])*dq[2]*dq[3] - 7.54746861368278e-49*cos(theta[1] + theta[2] + theta[3])*dq[3]**2
		return opL_7

	def opL8(self, q, dq, ddq, a, d, theta):
		opL_8 = 1.0*sin(theta[1] + theta[2] + theta[3])**2*dq[0]**2 + 6.12323399573677e-17*sin(theta[1] + theta[2] + theta[3])**2*dq[0]*dq[1] + 6.12323399573677e-17*sin(theta[1] + theta[2] + theta[3])**2*dq[0]*dq[2] + 6.12323399573677e-17*sin(theta[1] + theta[2] + theta[3])**2*dq[0]*dq[3] - 4.62149163970586e-65*sin(theta[1] + theta[2] + theta[3])**2*dq[1]**2 - 9.24298327941173e-65*sin(theta[1] + theta[2] + theta[3])**2*dq[1]*dq[2] - 9.24298327941173e-65*sin(theta[1] + theta[2] + theta[3])**2*dq[1]*dq[3] - 4.62149163970586e-65*sin(theta[1] + theta[2] + theta[3])**2*dq[2]**2 - 9.24298327941173e-65*sin(theta[1] + theta[2] + theta[3])**2*dq[2]*dq[3] - 4.62149163970586e-65*sin(theta[1] + theta[2] + theta[3])**2*dq[3]**2 - 1.22464679914735e-16*sin(theta[1] + theta[2] + theta[3])*cos(theta[1] + theta[2] + theta[3])*ddq[0] - 7.49879891330929e-33*sin(theta[1] + theta[2] + theta[3])*cos(theta[1] + theta[2] + theta[3])*ddq[1] - 7.49879891330929e-33*sin(theta[1] + theta[2] + theta[3])*cos(theta[1] + theta[2] + theta[3])*ddq[2] - 7.49879891330929e-33*sin(theta[1] + theta[2] + theta[3])*cos(theta[1] + theta[2] + theta[3])*ddq[3] + 6.12323399573677e-17*sin(theta[1] + theta[2] + theta[3])*ddq[0] + 7.49879891330929e-33*sin(theta[1] + theta[2] + theta[3])*ddq[1] + 7.49879891330929e-33*sin(theta[1] + theta[2] + theta[3])*ddq[2] + 7.49879891330929e-33*sin(theta[1] + theta[2] + theta[3])*ddq[3] - 1.0*cos(theta[1] + theta[2] + theta[3])**2*dq[0]**2 - 6.12323399573677e-17*cos(theta[1] + theta[2] + theta[3])**2*dq[0]*dq[1] - 6.12323399573677e-17*cos(theta[1] + theta[2] + theta[3])**2*dq[0]*dq[2] - 6.12323399573677e-17*cos(theta[1] + theta[2] + theta[3])**2*dq[0]*dq[3] + 4.62149163970586e-65*cos(theta[1] + theta[2] + theta[3])**2*dq[1]**2 + 9.24298327941173e-65*cos(theta[1] + theta[2] + theta[3])**2*dq[1]*dq[2] + 9.24298327941173e-65*cos(theta[1] + theta[2] + theta[3])**2*dq[1]*dq[3] + 4.62149163970586e-65*cos(theta[1] + theta[2] + theta[3])**2*dq[2]**2 + 9.24298327941173e-65*cos(theta[1] + theta[2] + theta[3])**2*dq[2]*dq[3] + 4.62149163970586e-65*cos(theta[1] + theta[2] + theta[3])**2*dq[3]**2 + 3.74939945665464e-33*cos(theta[1] + theta[2] + theta[3])*dq[0]**2 - 2.82984747192602e-81*cos(theta[1] + theta[2] + theta[3])*dq[0]*dq[1] - 2.82984747192602e-81*cos(theta[1] + theta[2] + theta[3])*dq[0]*dq[2] - 2.82984747192602e-81*cos(theta[1] + theta[2] + theta[3])*dq[0]*dq[3] - 4.62149163970586e-65*cos(theta[1] + theta[2] + theta[3])*dq[1]**2 - 9.24298327941173e-65*cos(theta[1] + theta[2] + theta[3])*dq[1]*dq[2] - 9.24298327941173e-65*cos(theta[1] + theta[2] + theta[3])*dq[1]*dq[3] - 4.62149163970586e-65*cos(theta[1] + theta[2] + theta[3])*dq[2]**2 - 9.24298327941173e-65*cos(theta[1] + theta[2] + theta[3])*dq[2]*dq[3] - 4.62149163970586e-65*cos(theta[1] + theta[2] + theta[3])*dq[3]**2
		return opL_8

	def opL9(self, q, dq, ddq, a, d, theta):
		opL_9 = 1.22464679914735e-16*sin(theta[1] + theta[2] + theta[3])*cos(theta[1] + theta[2] + theta[3])*dq[0]**2 + 7.49879891330929e-33*sin(theta[1] + theta[2] + theta[3])*cos(theta[1] + theta[2] + theta[3])*dq[0]*dq[1] + 7.49879891330929e-33*sin(theta[1] + theta[2] + theta[3])*cos(theta[1] + theta[2] + theta[3])*dq[0]*dq[2] + 7.49879891330929e-33*sin(theta[1] + theta[2] + theta[3])*cos(theta[1] + theta[2] + theta[3])*dq[0]*dq[3] - 2.82984747192602e-81*sin(theta[1] + theta[2] + theta[3])*cos(theta[1] + theta[2] + theta[3])*dq[1]**2 - 5.65969494385204e-81*sin(theta[1] + theta[2] + theta[3])*cos(theta[1] + theta[2] + theta[3])*dq[1]*dq[2] - 5.65969494385204e-81*sin(theta[1] + theta[2] + theta[3])*cos(theta[1] + theta[2] + theta[3])*dq[1]*dq[3] - 2.82984747192602e-81*sin(theta[1] + theta[2] + theta[3])*cos(theta[1] + theta[2] + theta[3])*dq[2]**2 - 5.65969494385204e-81*sin(theta[1] + theta[2] + theta[3])*cos(theta[1] + theta[2] + theta[3])*dq[2]*dq[3] - 2.82984747192602e-81*sin(theta[1] + theta[2] + theta[3])*cos(theta[1] + theta[2] + theta[3])*dq[3]**2 + 6.12323399573677e-17*sin(theta[1] + theta[2] + theta[3])*dq[0]**2 - 3.74939945665464e-33*sin(theta[1] + theta[2] + theta[3])*dq[0]*dq[1] - 3.74939945665464e-33*sin(theta[1] + theta[2] + theta[3])*dq[0]*dq[2] - 3.74939945665464e-33*sin(theta[1] + theta[2] + theta[3])*dq[0]*dq[3] - 9.84331363534125e-49*sin(theta[1] + theta[2] + theta[3])*dq[1]**2 - 1.96866272706825e-48*sin(theta[1] + theta[2] + theta[3])*dq[1]*dq[2] - 1.96866272706825e-48*sin(theta[1] + theta[2] + theta[3])*dq[1]*dq[3] - 9.84331363534125e-49*sin(theta[1] + theta[2] + theta[3])*dq[2]**2 - 1.96866272706825e-48*sin(theta[1] + theta[2] + theta[3])*dq[2]*dq[3] - 9.84331363534125e-49*sin(theta[1] + theta[2] + theta[3])*dq[3]**2 - 7.49879891330929e-33*cos(theta[1] + theta[2] + theta[3])**2*ddq[0] - 4.59169004331693e-49*cos(theta[1] + theta[2] + theta[3])**2*ddq[1] - 4.59169004331693e-49*cos(theta[1] + theta[2] + theta[3])**2*ddq[2] - 4.59169004331693e-49*cos(theta[1] + theta[2] + theta[3])**2*ddq[3] - 1.0*cos(theta[1] + theta[2] + theta[3])*ddq[0] - 1.22464679914735e-16*cos(theta[1] + theta[2] + theta[3])*ddq[1] - 1.22464679914735e-16*cos(theta[1] + theta[2] + theta[3])*ddq[2] - 1.22464679914735e-16*cos(theta[1] + theta[2] + theta[3])*ddq[3] + 7.49879891330929e-33*ddq[0] + 1.22464679914735e-16*ddq[1] + 1.22464679914735e-16*ddq[2] + 1.22464679914735e-16*ddq[3]
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

