import numpy as np
from CoolProp.CoolProp import *
from scipy.optimize import newton,bisect
g = 9.81
sig = 5.67e-8
patm = 101.325*1000.0
#
def K_C(Tc):
    return(Tc + 273.15)
def C_K(Tk):
    return(Tk - 273.15)
def calcul_h(T1,T2,L):
    Tf = (T1+T2)/2
    Beta = 1/Tf
    Dt = abs(T1-T2)
    rho = PropsSI('D', 'P', patm, 'T', Tf, 'Air')
    mu = PropsSI('V', 'P', patm, 'T', Tf, 'Air')
    k = PropsSI('L', 'P', patm, 'T', Tf, 'Air')
    Pr = PropsSI('Prandtl', 'P', patm, 'T', Tf, 'Air')
    nu = mu/rho
    al = nu/Pr
    Ra = g*Beta*Dt*L**3/(nu*al)
    Nu = (0.825+(0.387*Ra**(1/6))/((1+(0.492/Pr)**(9/16))**(8/27)))**2
    h = Nu*k/L
    return (h)
cas = 'd'
La = 0.05
Lb = 0.05
ka  = 0.75
kb  = 0.25
H = 2
W = 2
A = H*W
Ra = La/(ka*A)
Rb = Lb/(kb*A)
Ti = 20.0
To = -10.0
Tsky = -20.0
TskyK = K_C(Tsky)
ToK = K_C(To)
TiK = K_C(Ti)
ep = 0.85
if cas == 'a':
    hh = 0
    ep1 = 0.8
    ep2 = 0.8
elif cas == 'b':
    hh = 0
    ep1 = 0.1
    ep2 = 0.1
elif cas == 'c':
    hh = 1.5
    ep1 = 0.8
    ep2 = 0.8
else:
    hh = 1.5
    ep1 = 0.1
    ep2 = 0.1
K12 = 1/(1.0/ep1+1.0/ep2-1.0)
def fct2(T1c,T2c,T3c,T4c):
    Dt = T1c-T2c
    T1K = K_C(T1c)
    T2K = K_C(T2c)
    T3K = K_C(T3c)
    T4K = K_C(T4c)
    # Calcul du coef de convection interne
    hconvi = calcul_h(TiK,T1K,H)
    hconvo = calcul_h(T4K,ToK,H)
    # Calcul du coefficient de radiation interne
    hradi = sig*(T2K+T3K)*(T2K**2+T3K**2)/(1.0/ep1+1.0/ep2-1.0)
    Rradi = 1/(hradi*A)
    # Calcul du coef de convection externe
    hrado = ep*sig*(T4K**2+TskyK**2)*(T4K+TskyK)
    Rrado = 1/(hrado*A)
    Rconvi = 1/(hconvi*A)
    Rconvo = 1/(hconvo*A)
    Rint = 1/((hh+hradi)*A)
    #
    # Bilan radiatif à la surface externe de la vitre
    #
    Ri = Rconvi + Rint + Ra + Rb
    C1 = 1/Ri +  1/Rrado + 1/Rconvo
    T4n = (Ti/Ri + To/Rconvo + Tsky/Rrado)/C1
    qn = (Ti - T4n)/Ri
    T1n = Ti  - qn*Rconvi
    T2n = T1n  - qn*Ra
    T3n = T2n  - qn*Rint
    return T1n,T2n,T3n,T4n

T1i = 18
T2i = 16
T3i = 10
T4i = 0
ok = False
compt = 1
delta = 1e-3
itermax = 200
err = 0
while (ok != True):
    T1n,T2n,T3n,T4n = fct2(T1i,T2i,T3i,T4i)
    l1 = (abs(T1i-T1n)<delta)
    l2 = (abs(T2i-T2n)<delta)
    l3 = (abs(T3i-T3n)<delta)
    l4 = (abs(T4i-T4n)<delta)
    if l1 and l2 and l3 and l4:
        ok = True
    else:
        compt = compt+1
        T1i = T1n
        T2i = T2n
        T3i = T3n
        T4i = T4n
        if compt > itermax:
            err = 1
            ok = True

# verification
T1n,T2n,T3n,T4n = fct2(T1i,T2i,T3i,T4i)
T1K = K_C(T1n)
T2K = K_C(T2n)
T3K = K_C(T3n)
T4K = K_C(T4n)
hconvi = calcul_h(TiK,T1K,H)
hconvo = calcul_h(T4K,ToK,H)
Rconvi = 1/(hconvi*A)
Rconvo = 1/(hconvo*A)
q1 = (Ti - T1n)/Rconvi;print('q = ',q1,' W')
q2 = (T1n - T2n)/Ra;print('q = ',q2,' W')
q3a = K12*sig*A*(T2K**4 - T3K**4)
q3b = hh*A*(T2n - T3n)
q3 = q3a + q3b;print('q = ',q3,' W')
q4 = (T3n - T4n)/Rb;print('q = ',q4,' W')
q5r = ep*sig*A*(T4K**4 - TskyK**4)
q5c = (T4n - To)/Rconvo
q5 = q5r+q5c;print('q = ',q5,' W')
