import numpy as np
from CoolProp.CoolProp import *
# example 2.1
pi = np.pi
fluid = 'water'
patm = 101.325*1000.0
mp = 0.5
do = 0.04
di = 0.02
Ac = pi*(do**2 - di**2)/4
dh = do - di
a = di/do
Tmk = 300
Twk = 291
Cp = PropsSI('Cpmass','T',Tmk,'P',patm,fluid)
rho = PropsSI('D','T',Tmk,'P',patm,fluid)
mu = PropsSI('V','T',Tmk,'P',patm,fluid)
Pr = PropsSI('Prandtl','T',Tmk,'P',patm,fluid)
Pr = 5.83
Prw = PropsSI('Prandtl','T',Twk,'P',patm,fluid)
k = PropsSI('conductivity','T',Tmk,'P',patm,fluid)
u = mp/(rho*Ac)
Re = rho*u*dh/mu

k1 = 1.07 + 900/Re -0.63/(1+10*Pr)
Res = Re*((1+a**2)*np.log(a) + (1 - a**2))/((1-a)**2*np.log(a))
f = (0.781*np.log(Res)-1.5)**(-2)
bc = 1
if bc ==1:
    Fa = 0.75*a**-(0.17)
elif bc == 2:
    Fa = 0.9 - 0.15*a**(0.6)
else:
    Fa = (a*0.75*a**-(0.17) + (0.9 - 0.15*a**(0.6)))/(1+a)
# Gnieleski 2009
Nu1 = f/8*Re*Pr/(k1 + 12.7*np.sqrt(f/8)*(Pr**(2/3) - 1))*Fa
# Gnieleski 2015
Nu2 = f/8*(Re-1000)*Pr/(1 + 12.7*np.sqrt(f/8)*(Pr**(2/3) - 1))*Fa
# Gnieleski circular pipe
f2 = (0.79*np.log(Re)-1.64)**(-2)
Nu3 = f2/8*(Re-1000)*Pr/(1 + 12.7*np.sqrt(f2/8)*(Pr**(2/3) - 1))
# Petukhov correction
if a > 0.2:
    zet = 1.0
else:
    zet = 1 + 7.5*((1/a- 5)/Re)**0.6
if bc ==1:
    Nu4 = Nu3*zet*0.86*a**(-0.16)
elif bc == 2:
    Nu4 = Nu3*(1 - 0.14*a**(0.6))
else:
    Nui = Nu3*zet*0.86*a**(-0.16)
    Nuo = Nu3*(1 - 0.14*a**(0.6))
    Nu4 = (a*Nui + Nuo)/(1+a)
print ('Nu Gnieleski 2009 = ','%.2f' % Nu1)
print ('Nu Gnieleski 2015 = ','%.2f' % Nu2)
print ('Nu Gnieleski pipe = ','%.2f' % Nu3)
print ('NuPetukhov  = ', '%.2f' % Nu4)
h1a = Nu1*k/dh
h2a = Nu2*k/dh
h3a = Nu3*k/dh
h4a = Nu4*k/dh
print ('h Gnieleski 2009 (BC-I) = ','%.2f' % h1a)
print ('h Gnieleski 2015 (BC-I)= ','%.2f' % h2a)
print ('h Gnieleski pipe (BC-I)= ','%.2f' % h3a)
print ('h Petukho= (BC-I) ','%.2f' % h4a)

bc = 2
if bc ==1:
    Fa = 0.75*a**-(0.17)
elif bc == 2:
    Fa = 0.9 - 0.15*a**(0.6)
else:
    Fa = (a*0.75*a**-(0.17) + (0.9 - 0.15*a**(0.6)))/(1+a)
# Gnieleski 2009
Nu1 = f/8*Re*Pr/(k1 + 12.7*np.sqrt(f/8)*(Pr**(2/3) - 1))*Fa
# Gnieleski 2015
Nu2 = f/8*(Re-1000)*Pr/(1 + 12.7*np.sqrt(f/8)*(Pr**(2/3) - 1))*Fa
# Petukhov correction
if a > 0.2:
    zet = 1.0
else:
    zet = 1 + 7.5*((1/a- 5)/Re)**0.6
if bc ==1:
    Nu4 = Nu3*zet*0.86*a**(-0.16)
elif bc == 2:
    Nu4 = Nu3*(1 - 0.14*a**(0.6))
else:
    Nui = Nu3*zet*0.86*a**(-0.16)
    Nuo = Nu3*(1 - 0.14*a**(0.6))
    Nu4 = (a*Nui + Nuo)/(1+a)
print ('Nu Gnieleski 2009 (BC-II)= ','%.2f' % Nu1)
print ('Nu Gnieleski 2015 (BC-II)= ','%.2f' % Nu2)
print ('Nu Gnieleski pipe (BC-II)= ','%.2f' % Nu3)
print ('NuPetukhov (BC-II) = ', '%.2f' % Nu4)
h1b = Nu1*k/dh
h2b = Nu2*k/dh
h3b = Nu3*k/dh
h4b = Nu4*k/dh
print ('h Gnieleski 2009 (BC-2)= ','%.2f' % h1b)
print ('h Gnieleski 2015 (BC-2)= ','%.2f' % h2b)
print ('h Gnieleski pipe (BC-2)= ','%.2f' % h3b)
print ('h Petukhov (BC-2) = ', '%.2f' %  h4b)
bc = 3
if bc ==1:
    Fa = 0.75*a**-(0.17)
elif bc == 2:
    Fa = 0.9 - 0.15*a**(0.6)
else:
    Fa = (a*0.75*a**-(0.17) + (0.9 - 0.15*a**(0.6)))/(1+a)
# Gnieleski 2009
Nu1 = f/8*Re*Pr/(k1 + 12.7*np.sqrt(f/8)*(Pr**(2/3) - 1))*Fa
# Gnieleski 2015
Nu2 = f/8*(Re-1000)*Pr/(1 + 12.7*np.sqrt(f/8)*(Pr**(2/3) - 1))*Fa
# Petukhov correction
if a > 0.2:
    zet = 1.0
else:
    zet = 1 + 7.5*((1/a- 5)/Re)**0.6
if bc ==1:
    Nu4 = Nu3*zet*0.86*a**(-0.16)
elif bc == 2:
    Nu4 = Nu3*(1 - 0.14*a**(0.6))
else:
    Nui = Nu3*zet*0.86*a**(-0.16)
    Nuo = Nu3*(1 - 0.14*a**(0.6))
    Nu4 = (a*Nui + Nuo)/(1+a)
print ('Nu Gnieleski 2009 = ','%.2f' % Nu1)
print ('Nu Gnieleski 2015 = ','%.2f' % Nu2)
print ('Nu Gnieleski pipe = ','%.2f' % Nu3)
print ('NuPetukhov  = ', '%.2f' % Nu4)
h1c = Nu1*k/dh
h2c = Nu2*k/dh
h3c = Nu3*k/dh
h4c = Nu4*k/dh
print ('h Gnieleski 2009(BC-3) = ','%.2f' % h1c)
print ('h Gnieleski 2015 (BC-3)= ','%.2f' % h2c)
print ('h Gnieleski pipe (BC-3)= ','%.2f' % h3c)
print ('h Petukhov  (BC-3)= ', '%.2f' %  h4c)

href = h1c

K = (Pr/Prw)**(0.11)
K2 = (1 + (dh/100)**.667)
