#coding: utf-8
import numpy as np
from CoolProp.CoolProp import *
pi = np.pi
from heat_exchanger_md import *

patm = 101.325*1000.0
Tci = 15
Thi = 170
N = 6
Ns = 2
Np = 2
mc = 2.5
D = 0.04
Tco = 45
Thoh = 140
Tcm= (Tci+Tco)/2+273.15
Thm= (Thi+Thoh)/2+273.15
Cp_h = PropsSI('Cpmass','T',Thm,'P',patm,'water')
hhi = PropsSI('H','T',Thi+ 273.15,'P',patm,'water')
Cp_c = PropsSI('Cpmass','T',Tcm,'P',patm,'water')
mu_c = PropsSI('VISCOSITY','T',Tcm,'P',patm,'water')
Pr_c = PropsSI('PRANDTL','T',Tcm,'P',patm,'water')
k_c = PropsSI('CONDUCTIVITY','T',Tcm,'P',patm,'water')
hvs = PropsSI('H','P',patm,'Q',1,'water')
hls = PropsSI('H','P',patm,'Q',0,'water')
hfg = hvs - hls
mc1 = mc/N
Ac = pi*D**2/4
u = mc/(Ac*1000)
Re_c = 4*mc1/(pi*D*mu_c)
print(Re_c)
Nuc = 0.023*Re_c**(4/5)*Pr_c**0.4
hc = k_c*Nuc/D
Cc = mc*Cp_c
q = Cc*(Tco - Tci)
hh2a = PropsSI('H','P',patm,'Q',1,'water')
hh = 375
#Cmin = min(Cc,Ch)
#Cmax = max(Cc,Ch)
#Cr = Cmin/Cmax
#qmax = Cmin*(Thi-Tci)
#eff = q/qmax
#NTU = shell_tube_NTU(eff,Cr,Ns)
#print(eff,NTU)
#Rtot = 1/(NTU*Cmin)
P = N*pi*D
L = 12.06923633217468
print(L/(Ns*N))

#
# augmentation de la surface d'échange
#
Ph = 4*P
eta = 0.75
Ah = eta*Ph*L
Ac = P*L
Tsat = 100
ok = False
compt_max = 100
compt = 1
delta = 0.01
La = 0.75*L
hhb = hh*10
Rpa =1/(eta*Ph*hh) + 1/(P*hc)
Rpb =1/(eta*Ph*hhb) + 1/(P*hc)
qmaxb = Cc*(Tsat - Tci)

while not ok:
    Lb = L - La
    Rb = Rpb/Lb
    NTUb = 1/(Cc*Rb)
    effb = 1 - np.exp(-NTUb)
    qb = qmaxb*effb
    Tcia = Tci + qb/Cc
    mh = qb/hfg
    Ch = mh*Cp_h
    Cmax = max(Cc,Ch)
    Cmin = min(Cc,Ch)
    Cr = Cmin/Cmax
    qmaxa = Cmin*(Thi - Tcia)
    qa = Ch*(Thi-Tsat)
    effa = qa/qmaxa
    NTUa = shell_tube_NTU(effa,Cr,Ns)
    Ra = 1/(NTUa*Cmin)
    Lan = Rpa/Ra
    diff = abs(Lan-La)
    if diff < delta:
        ok = True
    else:
        La =Lan
        compt = compt + 1
        if compt > compt_max:
            err = 1
            print('erreur')
            ok = True
ho = hvs - qb/mh
x = (ho - hls)/hfg
Tcoa = Tcia + qa/Cc
print('La = ', La)
print('Lb = ', Lb)
print('Tcoa = ', Tcoa)
print(' x = ',x)
q1 = Cc*(Thi - Tci)
hi = PropsSI('H','T',Thi+ 273.15,'P',patm,'water')
hom= PropsSI('H','T',Tci+ 273.15,'P',patm,'water')
q2 = mh*(hi - hom)
qmax = min(q1,q2)
eff = (qa+qb)/qmax
print(eff)

