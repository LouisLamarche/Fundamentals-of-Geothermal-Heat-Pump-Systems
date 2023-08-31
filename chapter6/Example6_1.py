
#
# Example 6.1 written by Louis Lamarche 22 sept 2017
#
import numpy as np
from geothermal_md  import *
from  matplotlib.pyplot import *

cas = 'ICS'
fluid = 'Water'
patm = 101.325*1000.0
als = 0.1       # m2/day
alhr = als/24.0 # m2/hr
rb = 0.15/2.0
ks = 2.0
To = 10.0
qb = -40.0
H = 100.0
#
# borehole
di,do = sdr_pipe(1.25,11)     #  SDR-11 1.25 po nominal
ro = do/2.0
ri = di/2.0
xc = rb/3
Vp = 0.8                # débit en l/s
Trefk = 30.0+273
rhof = PropsSI('D', 'T', Trefk, 'P', patm, fluid)
Cp = PropsSI('Cpmass', 'T', Trefk, 'P', patm, fluid)
Pr = PropsSI('Prandtl','T',Trefk,'P',patm,fluid)
kf = PropsSI('conductivity','T',Trefk,'P',patm,fluid)
mu = PropsSI('viscosity','T',Trefk,'P',patm,fluid)

mp = Vp*rhof/1000.0
kp = 0.4       # conductivité du plastique
#
# Calcul de la résistance de la conduite
#
kg = 1.7   # coulis    mu = PropsSI('viscosity','T',Trefk,'P',patm,fluid)
# pipe resistance
rcond = np.log(ro/ri)/(2*pi*kp)
Re = 4*mp/(pi*di*mu)
if (Re>2300.0):
        # Gnielienski
        f = (0.79*np.log(Re)- 1.64)**-2
        Nud=((Re-1000.)*f*Pr/8.)/(1.+12.7*np.sqrt(f/8.)*(Pr**(2./3.)-1))
else:
    Nud = 3.6
    disp('Careful laminar')
hf = (Nud*kf)/(di)
rconv = 1/(pi*di*hf)
Rp = rcond + rconv
sigma = (kg-ks)/(kg+ks)


Rg,Ra = Rb_linesource(kg,ks,rb,ro,xc)
Rb = Rg + Rp/2
Ra = Ra + 2*Rp

rbb = rb/H
zob = 0.04
q1 = 1000.0
q2 = 1800.0
q3 = -200.0
q1p = q1/H
q2p = q2/H
q3p = q3/H
# intermediate time in hours



t1 = 1.0
t2 = 3.0
t3 = 5.0
dt = 0.01
t = np.arange(0,t2,dt)
Fo1 = alhr*t1/rb**2
Fo2 = alhr*t2/rb**2
TLS1 = (q1*G_function(Fo2)+(q2-q1)*G_function(Fo2-Fo1))/ks
TLS2 = (q1*G_function_ils(Fo2)+(q2-q1)*G_function_ils(Fo2-Fo1))/ks
CCf = mp*Cp
Tfo = 5.0
Tfi = -q2/CCf + Tfo
Tf = (Tfi + Tfo)/2.0
H1 = (TLS1+q2*Rb)/(To - Tf)  # Tb with ICS
H2 = (TLS2+q2*Rb)/(To - Tf)  # Tb with ILS
if cas == 'ICS':
    H = H1
else:
    H = H2
print('The totla length is = ',H)
q1p = q1/H
q2p = q2/H
q3p = q3/H
nh = len(t)
qp = np.zeros(nh)
DT1 = np.zeros(nh)
DT2 = np.zeros(nh)
#
# Calcululation of Tb
#
for i in range(0,nh):
    Fo = alhr*t[i]/rb**2
    if t[i] < t1:
        DT1[i] = -q1p*G_function(Fo)/ks             # ICS
        DT2[i] = -q1p*G_function_ils(Fo)/ks         # ILS
        qp[i] = q1p
    elif t[i] < t2:
        DT1[i] = -(q1p*G_function(Fo)+(q2p-q1p)*G_function(Fo-Fo1))/ks
        DT2[i] = -(q1p*G_function_ils(Fo)+(q2p-q1p)*G_function_ils(Fo-Fo1))/ks
        qp[i] = q2p
    else:
        DT1[i] = -(q1p*G_function(Fo)+(q2p-q1p)*G_function(Fo-Fo1)+(q3p-q2p)*G_function(Fo-Fo2))/ks
        DT2[i] = -(q1p*G_function_ils(Fo)+(q2p-q1p)*G_function_ils(Fo-Fo1)+(q3p-q2p)*G_function_ils(Fo-Fo2))/ks
        qp[i] = q3p
DT1 = DT1 - qp*Rb
DT2 = DT2 - qp*Rb
Tf1 = To + DT1
Tf2 = To + DT2
q = qp*H
Tfo1 = Tf1 + q/(2*mp*Cp)
Tfi1 = Tf1 - q/(2*mp*Cp)
Tfo2 = Tf2 + q/(2*mp*Cp)
Tfi2 = Tf2 - q/(2*mp*Cp)
plot(t,Tfo1,t,Tfo2)
legend(('ICS','ILS'))
show()
