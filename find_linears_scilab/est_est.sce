std_dev = stdev(Chi, 'c', mean(Chi,'c'));
disp([mean(Chi,'c'), std_dev]);

t = Symbol('t')
q_1 = Function('q_1')(t)
theta_1 = Function('theta_1')(q_1)
#q_1 = Symbol('q_1')
dq_1 = diff(q_1, t)
L = sin(theta_1)**2 * dq_1
dL_Ddq = diff(L, dq_1).doit()
dLdq_Dt = diff(dL_Ddq, t).doit()
dL_Dq = diff(L, q_1).doit()
opL = dLdq_Dt - dL_Dq
diff(theta_1, q_1)
Subs(Derivative(theta_1(_xi_2), _xi_2), (_xi_2), (q_1(t)))

