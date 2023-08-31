#coding: utf-8
#-------------------------------------------------------------------------------
# Exemple 7.2
#
import  numpy as np
from scipy.optimize import newton
from geothermal_md import *
from hydraulic_md import *
from conversion_md import *
from CoolProp.CoolProp import *

g = 9.8
fluide = 'INCOMP::APG-20%'  # ASHRAE propylene glycl 20 %
patm = 101.325*1000.0
#
L = 30.0
sdr = 11
NomD = 2
epsilon = 0
di,do = sdr_pipe(NomD,sdr)
T_refk1 = 32+273.15
T_refk2 = 0+273.15
Ql = 1.9      # flow rate  l/s
Q = Ql/1000.0   # flow rate  en m3/s
A = pi*di**2/4.0
u = Q/A
pour = 0.2
ed = epsilon/di       # rugosite e/D
rho =PropsSI('D','T',T_refk1,'P',patm,fluide)
Cp =PropsSI('Cpmass','T',T_refk1,'P',patm,fluide)
mu = PropsSI('viscosity','T',T_refk1,'P',patm,fluide)
nu = mu/rho
Re = di*u/nu
f = Colebrook(Re,ed)
h_pipe = f*(L/di)*u**2/(2*g) # Dp en metres
print ('h pipe  (32 C) = ', '%.3f' % h_pipe, ' m')
rho =PropsSI('D','T',T_refk2,'P',patm,fluide)
Cp =PropsSI('Cpmass','T',T_refk2,'P',patm,fluide)
mu = PropsSI('viscosity','T',T_refk2,'P',patm,fluide)
nu = mu/rho
Re = di*u/nu
f = Colebrook(Re,ed)
h_pipe = f*(L/di)*u**2/(2*g) # Dp en metres
print ('h pipe  (0C) = ', '%.3f' % h_pipe, ' m')
