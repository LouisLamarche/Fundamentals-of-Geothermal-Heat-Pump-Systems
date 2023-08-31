import numpy as np
from conversion_md import *
from scipy.optimize import newton

nn = 3
p2 = np.zeros(nn+1)

def p_e_1535_2961():   # q est en gpm, h est en ft
    p = np.array([-5.53062205e-06, -8.29116036e-04,  2.72102009e-02,  8.60160052e+01])
    return p
def r_p_e_1535_2961():   # q est en gpm
    p = np.array([ 5.53752811e-06, -6.10863076e-03,  1.02031965e+00,  2.80609741e+01])
    return p
def h_p_e_1535_2442(q):   # q est en gpm, h est en ft
    p = np.array([-7.41889082e-06, -6.72621576e-04,  1.36309324e-02,  5.85747055e+01])
    h = np.polyval(p,q)                            # h est en metres
    return h
def h_p_e_1535_2961(q):   # q est en gpm, h est en ft
    p = np.array([-5.53062205e-06, -8.29116036e-04,  2.72102009e-02,  8.60160052e+01])
    h = np.polyval(p,q)                            # h est en metres
    return h
def rend_p_e_1535_2961(q):   # q est en gpm
    p = np.array([ 5.53752811e-06, -6.10863076e-03,  1.02031965e+00,  2.80609741e+01])
    h = np.polyval(p,q)                            # h est en metres
    return h

def rend_p_e_1535_2442(q):   # q est en l/s, h est en m
    p = np.array([ 1.06138779e-05, -9.12057126e-03,  1.24488925e+00,  2.79501874e+01])
    h = np.polyval(p,q)                            # h est en metres
    return h



P1 = p_e_1535_2961()
r1 = r_p_e_1535_2961()
P3  = np.zeros(nn+1)
r3  = np.zeros(nn+1)


def pertes(q):
    hs = 0.0042*q**2
    return hs
def fct1(q):
    hs = pertes(q)
    hp =  h_p_e_1535_2961(q)
    y = hs - hp
    return y
def fct2(q):
    hs = pertes(q)
    hp =  h_p_e_1535_2442(q)
    y = hs - hp
    return y
#
# Eqs 7.37 7.38
#
yy = 2442/2961
for i in range(0,nn+1):
    ii = i-1
    ij = i - 3
    P3[i] = P1[i]*yy**ii
    r3[i] = r1[i]*yy**ij
def fct3(q):
    hs = pertes(q)
    hp =  np.polyval(P3,q)
    y = hs - hp
    return y
#
# a)
# solve at 2961 rpm
#
qn1 = newton(fct1,100)
h1 = pertes(qn1)
rend1 = rend_p_e_1535_2961(qn1)
print (' head losses at 2961 rpm is : {:.2f} ft'.format(h1))
print (' efficiency at 2961 rpm is : {:.0f} %'.format(rend1))
qn2 = newton(fct2,100)
h2 = pertes(qn2)
qn3 = newton(fct3,100)
h3 = pertes(qn3)
rend2 = rend_p_e_1535_2442(qn2)
rend3 = np.polyval(r3,qn3)
print (' head losses at 2442 rpm from Eqs 7.37 and 7.38 is : {:.2f} ft'.format(h3))
print (' efficiency at 2442 rpm from Eqs 7.37 and 7.38  is : {:.0f} %'.format(rend3))
print (' head losses at 2442 rpm from data sheet is : {:.2f} ft'.format(h2))
print (' efficiency at 2442 rpm from data sheet  is : {:.0f} %'.format(rend2))
