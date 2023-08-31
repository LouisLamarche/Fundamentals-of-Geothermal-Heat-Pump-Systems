#coding: utf-8
from numpy import *
from scipy.optimize import fsolve
from CoolProp.CoolProp import *
# exemple 2.10
patm = 101.325*1000.0
G = 750
alp = 0.09
tau = 0.83
A = 2.0
Ti = 22.0
To = -10.0
Tenv1 = -10.0
Tenvk = Tenv1 + 273.15
Tok = To + 273.15
Tik = Ti + 273.15
DT = Ti-To
ep = 0.84
sig = 5.67e-8
L1 = 0.005
W = 1
H = 2
u = 9

def calcul_h(Tsok,Tsik):
    Tmok = (Tsok + Tok)/2
    Tmik = (Tsik + Tik)/2
    rhoo = PropsSI('D','T',Tmok,'P',patm,'air')
    muo = PropsSI('viscosity','T',Tmok,'P',patm,'air')
    Pro = PropsSI('Prandtl', 'T', Tmok, 'P', patm,'air')
    kfo = PropsSI('conductivity', 'T', Tmok, 'P', patm, 'air')
    nuo = muo/rhoo
    Re = u*W/nuo
    Nuo  = (0.037*Re**0.8-871)*Pro**(1/3)
    ho  = Nuo*kfo/W
    # internal
    rhoi = PropsSI('D','T',Tmik,'P',patm,'air')
    mui = PropsSI('viscosity','T',Tmik,'P',patm,'air')
    Pri = PropsSI('Prandtl', 'T', Tmik, 'P', patm,'air')
    kfi = PropsSI('conductivity', 'T', Tmik, 'P', patm, 'air')
    nui = mui/rhoi
    alpi = nui/Pri
    Beta = 1/Tmik
    Ra = 9.81*Beta*abs(Tik - Tsik)*H**3/(nui*alpi)
    Nui = (0.825+(0.387*Ra**(1/6))/((1+(0.492/Pri)**(9/16))**(8/27)))**2
    hi  = Nui*kfi/H
    return ho,hi

k1 = 1.0
Rppc = L1/k1
Tsok = 0 + 273.15 # hypothese
Tsik = 0 + 273.15 # hypothese
ho,hi = calcul_h(Tsok,Tsik)
hrado = ep*sig*(Tsok+Tenvk)*(Tsok**2+Tenvk**2)
hradi = ep*sig*(Tsik+Tik)*(Tsik**2+Tik**2)
qsol = 0
hext = ho + hrado
hint = hi + hradi
Rppext = 1.0/hext
Rppint = 1.0/hint
Rpptot = Rppext + Rppint + Rppc
U = 1.0/Rpptot
qppa = U*DT
Tsin = Ti - qppa*Rppint
Tson = Tsin - qppa*Rppc
print ('a) Heat losses =  ',qppa, ' W/m2')
print ('a) Global loss coefficient U  =   ',U ,' W/m2 K')
#
# solution exacte
#
def fct(x):
    y = zeros(2)
    Tsok = x[0]
    Tsik = x[1]
    ho,hi = calcul_h(Tsok,Tsik)
    qout = ho*(Tsok-Tok) + ep*sig*(Tsok**4-Tenvk**4)
    qin1 = (Tsik - Tsok)/Rppc
    qin = hi*(Tik - Tsik) - ep*sig*(Tsik**4-Tik**4)
    y[0] = qout - qin1
    y[1] = qin - qin1
    return y


T = fsolve(fct,[272,273])
Tsokb = T[0]
Tsikb = T[1]
ho,hi = calcul_h(Tsok,Tsik)
Tsonn = Tsokb - 273.15
Tsinn = Tsikb - 273.15
hradob = ep*sig*(Tsokb+Tenvk)*(Tsokb**2+Tenvk**2)
hradib = ep*sig*(Tsikb+Tik)*(Tsikb**2+Tik**2)
hextb = ho + hradob
hintb = hi + hradib
Rppextb = 1.0/hextb
Rppintb = 1.0/hintb
Rpptotb = Rppextb + Rppintb + Rppc
Ub = 1.0/Rpptotb
qppa2 = Ub*DT
Tsin = Ti - qppa2*Rppintb
Tson = Tsin - qppa2*Rppc
print ('a) Heat losses =  ',qppa2, ' W/m2')
print ('a) Global loss coefficient U  =   ',Ub ,' W/m2 K')
#
# b
#
Tso = (alp*G +Ti/(Rppint+Rppc) + To/Rppext)/(1/(Rppint+Rppc) + 1/Rppext)
qppb = (Ti-Tso)/(Rppint+Rppc)
Tsi = Ti - qppb*Rppint
Tsoo = Tsi - qppb*Rppc
qppcon = (Tso - To)/Rppext
print('Approxiamte solution')
print ('b) Tso =  ',Tso,Tsoo)
print ('b) Tsi =  ',Tsi)
print ('b) Heat losses =  ',qppb, ' W/m2')

gain_sol1 = qppa - qppb
gain_sol2 = tau*G
gain_sol = gain_sol1 + gain_sol2
print ('b) Solar gain =  ',gain_sol)

#
# solution exacte
#
def fctb(x):
    y = zeros(2)
    Tsok = x[0]
    Tsik = x[1]
    ho,hi = calcul_h(Tsok,Tsik)
    qout = ho*(Tsok-Tok) + ep*sig*(Tsok**4-Tenvk**4)
    qin1 = (Tsik - Tsok)/Rppc
    qin = hi*(Tik - Tsik) - ep*sig*(Tsik**4-Tik**4)
    y[0] = qout - qin1 - alp*G
    y[1] = qin - qin1
    return y

T = fsolve(fctb,[272,273])
Tsokc = T[0]
Tsikc = T[1]
Tsocn = Tsokc - 273.15
Tsicn = Tsikc - 273.15
hradoc = ep*sig*(Tsokc+Tenvk)*(Tsokc**2+Tenvk**2)
hradic = ep*sig*(Tsikc+Tik)*(Tsikc**2+Tik**2)
hextc = ho + hradoc
hintc = hi + hradic
Rppextc = 1.0/hextc
Rppintc = 1.0/hintc
qppb2 = (Ti-Tsocn)/(Rppintc+Rppc)
Tsi = Ti - qppb2*Rppintc
Tso = Tsi - qppb2*Rppc
print('Exact solution')
print ('b) Tso =  ',Tso)
print ('b) Tsi =  ',Tsi)
print ('b) Heat losses =  ',qppb2, ' W/m2')

gain_sol1 = qppa2 - qppb2
gain_sol2 = tau*G
gain_sol = gain_sol1 + gain_sol2
print ('b) Solar gain =  ',gain_sol)
SHGC = ((qppa - qppb) + tau*G)/G
print ('SHGC  ',SHGC)
SHGCb = tau + alp*U/hext
print ('SHGC  ',SHGCb)
SHGC = ((qppa2 - qppb2) + tau*G)/G
print ('SHGC  ',SHGC)
SHGCb = tau + alp*U/hext
print ('SHGC  ',SHGCb)