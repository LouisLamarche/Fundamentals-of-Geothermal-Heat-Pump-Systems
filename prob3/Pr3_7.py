#coding: latin-1
#
import numpy as np
from scipy.optimize import newton
from conversion_md import *
from CoolProp.CoolProp import *
from CoolProp.HumidAirProp import *
from tabulate import tabulate

def KC(x):
    return x - 273.15
def CK(x):
    return x + 273.15
pi = np.pi
# data
di = 0.016
do = 0.021
D = 0.180
b = 0.035
N = 8
om = 2*pi/b
R = D/2
L = np.sqrt(1 + (om*R)**2)/om*2*pi*N
H = N*b
L2 = np.sqrt(1 + (om*R)**2)*H
Ao = pi*di*L
def calcul_ep(NTU,Cr):
    ep = (1 - np.exp(-NTU*(1-Cr)))/(1 - Cr*np.exp(-NTU*(1-Cr)))
    return ep
fluid = 'R410a'
patm = 101.325*1000.0
p1 = p4 = 1e6 # pascal
p2 = p3 = 2e6 # pascal
ps = p2
etas = 0.6
supp = 3
subc = 2
qevap = W_BTUhr(36000) # W
#
T4k = PropsSI('T', 'P', p4, 'Q', 1, fluid)
T4c = KC(T4k)
T1c = T4c + supp
T1k = T4k + supp
h1 = PropsSI('H', 'T',T1k,'P', p1, fluid)
s1 = PropsSI('S', 'T',T1k,'P', p1, fluid)
h2s = PropsSI('H', 'P', p2, 'S', s1, fluid)
h2 = h1 +  (h2s - h1)/etas
s2 = PropsSI('S', 'P', p2, 'H', h2, fluid)
T2k = PropsSI('T', 'P', p2, 'H', h2, fluid)
T2c = KC(T2k)
T3ks = PropsSI('T', 'P', p3, 'Q', 0, fluid)
T3s = KC(T3ks)
T3k = T3ks - subc
T3c = KC(T3k)
h3 = PropsSI('H', 'T',T3k,'P', p3, fluid)
s3 = PropsSI('S', 'T',T3k,'P', p3, fluid)
h4 = h3
x4 = PropsSI('Q', 'P', p4, 'H', h4, fluid)
s4 = PropsSI('S', 'P', p4, 'Q', x4, fluid)
qin = h1 - h4
qout = h2 - h3
win = h2  - h1
COP = qin/win
mpr = qevap/qin
qr = mpr*(h2 - h3)
Pow = mpr*(h2 - h1)
def calcul_Uo(Re):
    Uo = 110.1*Re**0.22
    return Uo
Tcik = CK(20)
Tcok = CK(25)
Trefk = (Tcik+Tcok)/2
rho = PropsSI('D','T',Trefk,'P',patm,'Water')
Cp = PropsSI('Cpmass','T',Trefk,'P',patm,'Water')
mu = PropsSI('viscosity','T',Trefk,'P',patm,'Water')
nu = mu/rho
Dh = do - di
P = pi*(di+do)
mpw = 0.2
Ac = pi*(do**2 - di**2)/4
um = mpw/(rho*Ac)
Re = Dh*um/nu
Re2 = mpw*Dh/(Ac*mu)
Re3 = 4*mpw/(P*mu)
Uo = calcul_Uo(Re)
UA = Uo*Ao
Cc = mpw*Cp
# initial guess
Tsk = T3k + 3
hs = PropsSI('H', 'T',Tsk,'P', ps, fluid)
hmm = PropsSI('H', 'T',Tcik,'P', p3, fluid)
def calcul_qsurc(hs):
    Tsk = PropsSI('T', 'H',hs,'P', p3, fluid)
    Ch = mpr*(h2 - hs)/(T2k - Tsk)
    Cmin = min(Ch,Cc)
    Cmax = max(Cc,Ch)
    Cr = Cmin/Cmax
    NTU = UA/Cmin
    ep = calcul_ep(NTU,Cr)
    qmax = min(Cc*(T2k - Tcik),mpr*(h2 - hmm))
    q = ep*qmax
    hsn = h2 - q/mpr
    return hsn - hs

hf = newton(calcul_qsurc,hs)
Tsk = PropsSI('T', 'H',hf,'P', ps, fluid)
ss = PropsSI('S', 'H',hf,'P', ps, fluid)
if Tsk <= T3ks:   # saturated state
    print('vapor is saturated, split the HE in two parts')
else:
    print('vapor is superheated, ok')
Ch = mpr*(h2 - hf)/(T2k - Tsk)
Cmin = min(Ch,Cc)
Cmax = max(Cc,Ch)
Cr = Cmin/Cmax
NTU = UA/Cmin
ep = calcul_ep(NTU,Cr)
Tsc = KC(Tsk)
qsurc = mpr*(h2 - hf)
Tcok = Tcik + qsurc/Cc
Tcoc = KC(Tcok)
print ('mdot = ', '%.3f' %  mpr ,' kg/s')
print('qsurc = ', '%.3f' %  (qsurc/1000), ' kW')
print('Twater out = ', '%.2f' % Tcoc, ' C')
print('NTU = ','%.3f' % NTU)
print('ep = ','%.3f' % ep)
head = ['point','T(C)', 'p(kPa)','h(kJ/kg)','s(kJ/kg-K)']
pt1 = ['1','%.1f' % T1c, '%.1f' % (p1/1000),  '%.0f' % (h1/1000), '%.3f' % (s1/1000)]
pt2 = ['2','%.1f' % T2c, '%.1f' % (p2/1000),  '%.0f' % (h2/1000), '%.3f' % (s2/1000)]
pts = ['s','%.1f' % Tsc, '%.1f' % (ps/1000),  '%.0f' % (hs/1000), '%.3f' % (ss/1000)]
pt3 = ['3','%.1f' % T3c, '%.1f' % (p3/1000),  '%.0f' % (h3/1000), '%.3f' % (s3/1000)]
pt4 = ['4','%.1f' % T4c, '%.1f' % (p4/1000),  '%.0f' % (h4/1000), '%.3f' % (s4/1000)]
print(tabulate([pt1,pt2,pts,pt3,pt4], headers= head))

