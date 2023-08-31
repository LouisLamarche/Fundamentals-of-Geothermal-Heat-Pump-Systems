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
mh = 3.0
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
Ch = mh*Cp_h
Tho = Thi - q/Ch
hho = PropsSI('H','T',Tho+ 273.15,'P',patm,'water')
print('Tho = ',Tho)
hh = 375
Cmin = min(Cc,Ch)
Cmax = max(Cc,Ch)
Cr = Cmin/Cmax
qmax = Cmin*(Thi-Tci)
eff = q/qmax
NTU = shell_tube_NTU(eff,Cr,Ns)
print(eff,NTU)
Rtot = 1/(NTU*Cmin)
P = N*pi*D
U = (1/hh + 1/hc)**(-1)
Rptot = 1/(P*U)
L = Rptot/Rtot
print('L total = ',L)
print('L 1 tube = ',L/(Ns*N*Np))
