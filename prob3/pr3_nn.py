#coding: utf-8
#
import numpy as np
from conversion_mod import *
from CoolProp.CoolProp import *
from CoolProp.HumidAirProp import *
from tabulate import tabulate

#coding: latin-1
import numpy as np
from properties_mod import *

def shell_tube_eff(NTU,Cr,nshell=1):
#
# shell and tube (nshell shell pass)
# NTU is the total NTU
    NTUn = NTU/nshell
    G = np.sqrt(1+Cr**2)
    y = np.exp(-NTUn*G)
    ep1 = 2/(1+Cr+G*(1+y)/(1-y))
    if nshell > 1:
        if Cr == 1:
            ep = nshell*ep1/(1+ep1*(n-1))
        else:
            z = (1-ep1*Cr)/(1-ep1)
            ep = (z**nshell-1)/(z**nshell-Cr)
    else:
        ep = ep1
    return ep

def shell_tube_NTU(ep,Cr,nshell=1):
#
# shell and tube (nshell shell pass)
# NTU is the total NTU
    G = np.sqrt(1+Cr**2)
    if nshell > 1:
        if Cr ==1:
            ep1 = ep/(nshell - ep*(n-1))
        else:
            F = ((ep*Cr -1)/(ep-1))**(1/nshell)
            ep1 = (F-1)/(F-Cr)
    else:
        ep1 = ep
    E = (2/ep1 - (1+Cr))/G
    if E > 1:
        NTU1 = -np.log((E-1)/(E+1))/G
        NTU = nshell*NTU1
    else:
        print('impossible')
        NTU = -999
    return NTU
def F_coef_shell_tube(Tci,Tco,Thi,Tho,N):
    P = (Tco-Tci)/(Thi-Tci)
    R = (Thi-Tho)/(Tco-Tci)
    A = 2/P - 1 - R
    B = 2/P*np.sqrt((1-P)*(1-P*R))
    num = np.sqrt(R**2 + 1)/(R-1)*np.log((1-P)/(1-P*R))
    if N == 1:
        den = np.log((A+np.sqrt(R**2+1))/(A-np.sqrt(R**2+1)))
        F = num/den
    else:
        den = 2*np.log((A+B+np.sqrt(R**2+1))/(A+B-np.sqrt(R**2+1)))
        F = num/den
    return F


def counter_flow_NTU(ep,Cr):
    if Cr < 1:
        NTU = 1/(Cr-1)*np.log((ep-1)/(ep*Cr-1))
    else:
        NTU = ep/(1.0-ep)
    return NTU

def counter_flow_eff(NTU,Cr):

    if Cr < 1:
        ep = (1-np.exp(-NTU*(1-Cr)))/(1-Cr*np.exp(-NTU*(1-Cr)))
    else:
        ep = NTU/(1+NTU)
    return ep





# data
def KC(x):
    return x - 273.15
def CK(x):
    return x + 273.15
#
# data
#
fluid = 'R134a'
patm = 101.325*1000.0
p1 = p4 = 1.8e5
p2 = p3 = 1e6
etac = 0.8
surch = 0
#
#
T4k = PropsSI('T', 'P', p4, 'Q', 0, fluid)
T4c = KC(T4k)
T1c = T4c + surch
T1k = T4k + surch
h1 = PropsSI('H', 'P', p1, 'Q', 1, fluid)
s1 = PropsSI('S', 'P', p1, 'Q', 1, fluid)
s2s = s1
h2s = PropsSI('H', 'P', p2, 'S', s2s, fluid)
h2 = h1 + (h2s - h1)/etac
T2k = PropsSI('T', 'P', p2, 'H', h2, fluid)
s2 = PropsSI('S', 'P', p2, 'H', h2, fluid)
T2c = KC(T2k)
T3k = PropsSI('T', 'P', p3, 'Q', 0, fluid)
T3c = KC(T3k)
h3 = PropsSI('H', 'P', p3, 'Q', 0, fluid)
s3 = PropsSI('S', 'P', p3, 'Q', 0, fluid)
h4 = h3
x4 = PropsSI('Q', 'P', p4, 'H', h4, fluid)
s4 = PropsSI('S', 'P', p1, 'Q', x4, fluid)
qin = h1 - h4
qout = h2 - h3
win = h2  - h1
COPC =qin/win
COPH =qout/win
print ('COP (cooling) = ', '%.2f' % COPC)
print ('COP (heating) = ', '%.2f' % COPH)
head = ['point','T(C)', 'p(kPa)','h(kJ/kg)','s(kJ/kg-K)']
pt1 = ['1','%.1f' % T1c, '%.1f' % (p1/1000),  '%.0f' % (h1/1000), '%.2f' % (s1/1000)]
pt2 = ['2','%.1f' % T2c, '%.1f' % (p2/1000),  '%.0f' % (h2/1000), '%.2f' % (s2/1000)]
pt3 = ['3','%.1f' % T3c, '%.1f' % (p3/1000),  '%.0f' % (h3/1000), '%.2f' % (s3/1000)]
pt4 = ['4','%.1f' % T4c, '%.1f' % (p4/1000),  '%.0f' % (h4/1000), '%.2f' % (s4/1000)]
print(tabulate([pt1,pt2,pt3,pt4], headers= head))


Thi = 2
Thik = KC(Thi)
Tci = T4c
hfg = h1 - h4
Cp = 4200
D = 0.015
A = pi*D**2/4
u = 1.5
Vp = u*A
mp = Vp*1000
Ch = mp*Cp
Cmin = Ch
qmax = Cmin*(Thi - Tci)
q = 10000
ep = q/qmax
NTU = -log(1-ep)
UA = Cmin*NTU
hc = 2200
hh = 22000
U = (1/hc + 1/hh)**(-1)
A = UA/U
print(U,A)
