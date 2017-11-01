"""
    UPDATE 1: MAY BE STRUCTURE HAD BEEN CHANGED!!! I DONT REMEMBER (:
    !!!WARNING!!!
    JUST TEMPELATE;
    Example of module for \xi matrix (regressor).
"""
# zero row
#!/usr/bin/env python3
from numpy import cos, sin


class XI:

    def __init__(self, q=(0,0,0,0,0), dq=(0,0,0,0,0), ddq=(0,0,0,0,0), a=(0,0,0,0,0), d=(0,0,0,0,0)):
        self.q, self.dq, self.ddq = q, dq, ddq
        self.a, self.d = a, d

    def setData(self, q, dq, ddq, a, d):
        self.q, self.dq, self.ddq = q, dq, ddq
        self.a, self.d = a, d

    def opL0(self):
        q, dq, ddq = self.q, self.dq, self.ddq
        a, d = self.a, self.d
        opL000 = 1
        return opL000

    def opL1(self):
        q, dq, ddq = self.q, self.dq, self.ddq
        a, d = self.a, self.d
        opL001 = 1
        return opL001

    """
        and so on... to opL9(self)
    """

    def getXI(self, q, dq, ddq, a, d):
        self.setData(q, dq, ddq, a, d)
        XI = [0 for i in range(10)]
        XI[0] = self.opL0()
        XI[1] = self.opL1()
        # XI[2] = self.opL2()
        # XI[3] = self.opL3()
        # XI[4] = self.opL4()
        # XI[5] = self.opL5()
        # XI[6] = self.opL6()
        # XI[7] = self.opL7()
        # XI[8] = self.opL8()
        # XI[9] = self.opL9()
        return XI

if __name__ == '__main__':

    xi = XI()
    x = xi.getXI((0,0,0,0,0),(0,0,0,0,0),(0,0,0,0,0),(0,0,0,0,0),(0,0,0,0,0))
    print(x)

    i, j, k = 0, 0, 0
    opL_sym = '1'
    opNum = str(k)

    module = '#!/usr/bin/env python3\n'
    module += 'from numpy import cos, sin\n'
    module += '\n\n'
    module += 'class XI:\n'
    module += '\t"""XI_' + str(i) + str(j) + '"""\n\n'
    module += '\tdef __init__(self, q=(0,0,0,0,0), dq=(0,0,0,0,0), ddq=(0,0,0,0,0), a=(0,0,0,0,0), d=(0,0,0,0,0)):\n' \
              '\t\tself.q, self.dq, self.ddq = q, dq, ddq\n' \
              '\t\tself.a, self.d = a, d\n\n'
    module += '\tdef setData(self, q, dq, ddq, a, d):\n' \
              '\t\tself.q, self.dq, self.ddq = q, dq, ddq\n' \
              '\t\tself.a, self.d = a, d\n\n'

    # methods contains computing operator L
    module += '\tdef opL' + opNum + '(self):\n' \
              '\t\t"""' + str(i) + str(j) + str(k) + '"""\n' \
              '\t\tq, dq, ddq = self.q, self.dq, self.ddq\n' \
              '\t\ta, d = self.a, self.d\n' \
              '\t\topL_' + opNum + ' = ' + opL_sym + '\n' \
              '\t\treturn opL_' + opNum + '\n\n'

    module += '\tdef getXI(self, q, dq, ddq, a, d):\n' \
              '\t\tself.setData(q, dq, ddq, a, d)\n' \
              '\t\tXI = [0 for i in range(10)]\n' \
              '\t\tXI[0] = self.opL0()\n' \
              '\t\tXI[1] = self.opL1()\n' \
              '\t\tXI[2] = self.opL2()\n' \
              '\t\tXI[3] = self.opL3()\n' \
              '\t\tXI[4] = self.opL4()\n' \
              '\t\tXI[5] = self.opL5()\n' \
              '\t\tXI[6] = self.opL6()\n' \
              '\t\tXI[7] = self.opL7()\n' \
              '\t\tXI[8] = self.opL8()\n' \
              '\t\tXI[9] = self.opL9()\n' \
              '\t\treturn XI\n\n'

    f = open('rm.py', 'w')
    f.write(module)