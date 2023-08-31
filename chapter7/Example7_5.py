#coding: utf-8
#-------------------------------------------------------------------------------
# Exemple  7.5
#
import numpy as np
from geothermal_md import *
from hydraulic_md import *
from conversion_md import *
from CoolProp.CoolProp import *


g = 9.81
epsilon = 0
fluide1 = 'water'
fluide2 = 'INCOMP::APG-20%'  # ASHRAE propylene glycl 20 % volume
fluide3 = 'INCOMP::APG-40%'  # ASHRAE propylene glycl 40 % volume
cas = 2
if cas ==1:
    fluide = fluide1
    eta = 0.47
elif cas ==2:
    fluide = fluide2
    eta = 0.47
else:
    fluide = fluide3
    eta = 0.37
patm = 101.325*1000.0
#
Cv = 1.08  # l/s/kPa**1/2
K = 2.5
L = 100.0
D = 0.04
Le = 2.0
T_refk = 278
Ql = 1.5       # debit el l/s
Q = Ql/1000.0   # debit en m3/s
A = pi*D**2/4.0
u = Q/A
ed  = epsilon/D
rhow =PropsSI('D','T',16+273.15,'P',patm,'water')
rho1 =PropsSI('D','T',T_refk,'P',patm,'water')
mu1 = PropsSI('viscosity','T',T_refk,'P',patm,'water')
nu1 = mu1/rho1
Re1 = D*u/nu1
f1 = Colebrook(Re1,ed)
rho =PropsSI('D','T',T_refk,'P',patm,fluide)
Cp =PropsSI('Cpmass','T',T_refk,'P',patm,fluide)
mu = PropsSI('viscosity','T',T_refk,'P',patm,fluide)
Sg = rho/rhow
nu = mu/rho
Re = D*u/nu
f = Colebrook(Re,ed)
corr = f/f1
h_pipe = f*(L/D)*u**2/(2*g) # Dp en metres
h_sing1 = K*u**2/(2*g)*corr
h_sing2 = f*(Le/D)*u**2/(2*g)
DPkPa = Sg*(Ql/Cv)**2     # psi
DPPa = DPkPa*1000.0     # Pa
h_valve = DPPa/(rho*g)*corr # m de fluide
h_sing = h_sing1 + h_sing2
h_tot = h_pipe + h_sing + h_valve
mp = Q*rho
Pc = mp*g*h_pipe
Pv = mp*g*h_valve
Ps = mp*g*h_sing
P = mp*g*h_tot
HP = P/W_hp()
print ('h pipe = ' + str(h_pipe) + ' m')
print ('h sing1 = ' + str(h_sing1) + 'm')
print ('h sing2 = ' + str(h_sing2) + 'm')
print ('h sing = ' + str(h_sing) + 'm')
print ('h valve = ' + str(h_valve) + ' m')
print ('h total = ' + str(h_tot) + ' m')
h_pi = ft_m(h_tot)
print ('P cond = ' + str(Pc) + ' W')
print ('P sing = ' + str(Ps) + ' W')
print ('P valve = ' + str(Pv) + ' W')
print ('P = ' + str(P) + ' W')
print ('P = ' + str(HP) + ' h.p.')
Ps = P/eta
HPs = Ps/W_hp()
print ('BHP = ' + str(Ps) + ' W')
print ('BHP = ' + str(HPs) + ' h.p.')


