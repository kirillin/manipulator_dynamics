# coding: utf-8
#!/usr/bin/env python3

#
# Module contains function for lazy processing expressions
#

import re

def octaveButyfiler(expr):
    """ clear code for scilab """
    ptq = r"[d]*q_[1-5]"
    ptex = r"Derivative\([d]*q_[1-5]+\(t\),[\w ,]+\)"
    ptt = r"\(t\)"
    pt = r"[t]+"

    def repldq(m):
        lala = re.search(ptq, m.group(0))
        t = re.findall(pt, m.group(0))
        if len(t) - 2 == 1:
            return "d" + lala.group(0)
        elif len(t) - 2 == 2:
            return "dd" + lala.group(0)

    def rmt(m):
        return ""

    rerdq = re.sub(ptex, repldq, expr)
    remt = re.sub(ptt, rmt, rerdq)
    return remt


def python_gencode(expr, thetas=None):
    """
        Using for creation XI modules
        python code cleaner

    """

    def replaceDotsQ(m):
        "replace Derivativ and replace q_i(t) on dq_i and next dq[i]"
        lala = re.search(r"[d]*q_[1-5]", m.group(0))
        t = re.findall(r"[t]+", m.group(0))
        num = re.search(r"[1-5]", lala.group(0))
        n = int(num.group(0))-1
        if len(t) - 2 == 1:
            return "dq[" + str(n) + "]"
        elif len(t) - 2 == 2:
            return "ddq[" + str(n) + "]"

    def replaceADi(m):
        paramName = re.search(r"[ad]", m.group(0))
        num = re.search(r"[1-5]", m.group(0))
        n = int(num.group(0))-1
        return str(paramName.group(0)) + "[" + str(n) + "]"

    def deleteT(m):
        return ""

    # def replaceQi(m):
    #     q_i = re.search(r"q_[1-5]", m.group(0))
    #     num = re.search(r"[1-5]", m.group(0))
    #     n = int(num.group(0))-1
    #     return str(thetas[n])

    def replaceQi(m):
        q_i = re.search(r"q_[1-5]", m.group(0))
        num = re.search(r"[1-5]", m.group(0))
        n = int(num.group(0)) - 1
        return "q[" + str(n) + "]"

    # replace all Derivative(.*)
    replace_dotsq = re.sub(r"Derivative\([d]*q_[1-5]+\(t\),[\w ,]+\)", replaceDotsQ, expr)

    # delete all "(t)"
    replace_t = re.sub(r"\(t\)", deleteT, replace_dotsq)

    # replace all q_i to (q[i] + XXX*pi)
    replace_qi = re.sub(r"q_[1-5]", replaceQi, replace_t)

    # replace DH parameters
    replace_adi = re.sub(r"[ad]_[1-5]", replaceADi, replace_qi)
    return replace_adi

if __name__ == '__main__':
    st = 'd_1 * sin(theta_1)'
    print(python_gencode(st))
