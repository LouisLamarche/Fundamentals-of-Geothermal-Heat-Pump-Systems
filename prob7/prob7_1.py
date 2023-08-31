#
# exemple de la methode de la ligne source pour évaluer la conductivité du sol
#
import numpy as np
from geothermal_md import *
from conversion_md import *
from hydraulic_md import *
Patm = 101.325*1000.0
#
#
#
# fluide
Trefk = 280.0
debit_vol = 12  # en l/s
debit_vol = debit_vol/1000.0        # en m3/s
fluide  = 'INCOMP::APG-40%'   # propylene - glycol 40 %
# calcul des proprietes
mu = PropsSI('viscosity', 'T', Trefk, 'P', Patm, fluide)
rho = PropsSI('D', 'T', Trefk, 'P', Patm, fluide)
Cp = PropsSI('Cpmass', 'T', Trefk, 'P', Patm, fluide)
Pr = PropsSI('Prandtl', 'T', Trefk, 'P', Patm, fluide)
kf = PropsSI('conductivity', 'T', Trefk, 'P', Patm, fluide)
nu = mu/rho
mp = debit_vol*rho
g = 9.8
#
# circuit exterieur
#
epp = 0
tube_loop = 1.25
tube_coll = 3
nb = 6
SDR = 11
d1,d2 = sdr_pipe(tube_coll,SDR)
ro = d2/2
ri = d1/2
A1 = pi*d1**2/4.0
Q1 =  debit_vol              # debit en m3/s
u1 = Q1/A1
Re1 = d1*u1/nu
ed = epp/d1
f1 = Colebrook(Re1,ed)
Leq_coude = 2.2
L_seg = 2*40.0
L_1 = L_seg + 6*Leq_coude
h_1 = f1*(L_1/d1)*u1**2/(2*g)  # % Dp en metres
#
# parallel path
#
d12,d22 = sdr_pipe(tube_loop,SDR)     # choix des tuyaux SDR-11 1.25 po nominal

ro = d22/2
ri = d12/2
A2 = pi*d12**2/4.0
Q2 =  debit_vol/nb              # debit en m3/s
u2 = Q2/A2
Re2 = d12*u2/nu
ed2 = epp/d12
f2 = Colebrook(Re2,ed2)
L_seg = 200.0
L_2 = L_seg
h_2 = f2*(L_2/d12)*u2**2/(2*g)  # % Dp en metres
h_t = h_1 + h_2
W = mp*g*h_t
print ('3 in head losses = ',h_1,' m')
print ('1.25 in head losses = ',h_2,' m')
print ('Total losses = ',h_t,' m')
print ('Fluid Power =  ',W,' W')
