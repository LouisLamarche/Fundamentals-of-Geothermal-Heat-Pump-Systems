import numpy as np
from geothermal_md import *
from conversion_md import *
from hydraulic_md import *
Patm = 101.325*1000.0
#
#
#
# fluide
Trefk =  290.0
#
# 2)  Calcul des pertes de charges entre PAC et collecteur
#
L_seg = 8.0
Ks = 1.5
#
# 3 premiers segments 11/4 po
# 2 autres 1 po
#
fluide  = 'Water'   # propylene - glycol 40 %
# calcul des proprietes
mu = PropsSI('viscosity', 'T', Trefk, 'P', Patm, fluide)
rho = PropsSI('D', 'T', Trefk, 'P', Patm, fluide)
Cp = PropsSI('Cpmass', 'T', Trefk, 'P', Patm, fluide)
Pr = PropsSI('Prandtl', 'T', Trefk, 'P', Patm, fluide)
kf = PropsSI('conductivity', 'T', Trefk, 'P', Patm, fluide)
nu = mu/rho
D = np.array([0.03,0.02])
npac = 3
htotal_int = 0
Qpac_nom = 0.0004
Qt = 0.0012            # m3/s
Qpac = Qt/npac
mp = Qt*rho
g = 9.81
Dp = 0.02
A1 = pi*Dp**2/4.0
Q1 = Qpac  # debit en m3/s
u1 = Q1/A1
Re1 = Dp*u1/nu
f1 = Colebrook(Re1,0.0)
h1 = f1*(L_seg/0.02)*u1**2/(2*g)
Q2 = 2*Qpac  # debit en m3/s
A2 = pi*0.03**2/4.0
u2 = Q2/A2
Re2 = 0.03*u2/nu
f2 = Colebrook(Re2,0.0)
h2 = f2*(L_seg/0.03)*u2**2/(2*g)
Cv2 = 0.291
Dp_valve3 = 2*(Qpac*1000/Cv2)**2  # 2 valves
h_valve = 1000*Dp_valve3/(rho*g)
h_pac =  4.5*(Qpac/Qpac_nom)**2  # en mètres
hs = Ks*u1**2/(2.0*g)
hp = h_pac+hs+h_valve;print (hp)
hi = h_pac + h_valve + hs
ht1 = hi + h1 + h2
ht2 = hi + 2*h2
htotal = hi + 2/3*h1+4/3*h2
print ('htotal = ',htotal)
print ('htotal = ',ht1)
print ('htotal = ',ht2)
Wp = mp*g*htotal
print ('Wp = ',Wp)