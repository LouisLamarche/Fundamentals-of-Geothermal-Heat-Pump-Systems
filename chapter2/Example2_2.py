import numpy as np
from CoolProp.CoolProp import *
# example 2.1
pi = np.pi
fluid = 'water'
patm = 101.325*1000.0
mp = 0.1
d = 0.025
Ac = pi*d**2/4
D = 0.5
b = 0.1
Dc = D*(1 + (b/(pi*D))**2)
Tmk = 300
Twk = 290
Cp = PropsSI('Cpmass','T',Tmk,'P',patm,fluid)
rho = PropsSI('D','T',Tmk,'P',patm,fluid)
mu = PropsSI('V','T',Tmk,'P',patm,fluid)
Pr = PropsSI('Prandtl','T',Tmk,'P',patm,fluid)
Prw = PropsSI('Prandtl','T',Twk,'P',patm,fluid)
muw = PropsSI('V','T',Twk,'P',patm,fluid)
k = PropsSI('conductivity','T',Tmk,'P',patm,fluid)
u = mp/(rho*Ac)
la = d/D
Rec1 = 2100*(1 + 12*la**0.5)    # eq 1.65
Rec2 = 2300*(1 + 8.6*la**0.45)  # eq 1.66
Re = rho*u*d/mu
Re2 = mp/(d*mu)
De = Re*np.sqrt(d/D)
# Manlapaz
x1 = (1 + 957/(De**2*Pr))**2
x2 = (1 + 0.477/Pr)
Nu1 = ((3.66 + 4.343/x1)**3 + 1.158*(De/x2)**1.5)**(1/3)*(mu/muw)**0.14
# Gnieleski 2015
m = 0.5 + 0.2903*(d/D)**0.194
Nu2 = 3.66 + 0.08*(1 + 0.8*(d/D)**(0.9))*Re2**m*Pr**(1/3)*(Pr/Prw)**0.14
h1 = Nu1*k/d
h2 = Nu2*k/d
print ('h  Manlapaz = ','%.2f' % h1)
print ('h Gnieleski 2015 = ','%.2f' % h2)
