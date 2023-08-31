import numpy as np
from CoolProp.CoolProp import *
# example 2.3
pi = np.pi
fluid = 'water'
patm = 101.325*1000.0
mp = 0.3
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
Rec = 2100*(1 + 12*la**0.5)
Re = rho*u*d/mu
De = Re*np.sqrt(d/D)
# Manlapaz
Decr = Rec*np.sqrt(d/D)
x1 = (1 + 957/(Decr**2*Pr))**2
x2 = (1 + 0.477/Pr)
Nul = ((3.66 + 4.343/x1)**3 + 1.158*(Decr/x2)**1.5)**(1/3)*(mu/muw)**0.14
# Gnieleski circular pipe
fsp = (0.79*np.log(Re)-1.64)**(-2)
Nut1 = fsp/8*(Re-1000)*Pr/(1 + 12.7*np.sqrt(fsp/8)*(Pr**(2/3) - 1))
# Gnieleski circular pipe
fc = (fsp + 0.03*np.sqrt(la))*(muw/mu)**0.27
Nut2 = fc/8*Re*Pr/(1 + 12.7*np.sqrt(fc/8)*(Pr**(2/3) - 1))*(Pr/Prw)**0.14
gam = (Re - Rec)/(22000 - Rec)
Nu = (1 - gam)*Nul + gam*Nut2
h1 = Nut1*k/d
h2 = Nu*k/d
print(h1,h2)