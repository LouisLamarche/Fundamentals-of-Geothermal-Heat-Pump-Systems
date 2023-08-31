from geothermal_md import *
from numpy import *
from  matplotlib.pyplot import *
from conversion_md import *
from hydraulic_md import *
#
#
# 2)  Calcul des pertes de charges entre PAC et collecteur
#
L_seg = 8
#
# 3 premiers segments 11/4 po
# 2 autres 1 po
#
Cp = 4180.
nu = 6.0e-7
rho = 1000
D1 = 0.038
D2 = 0.05
D3 = 0.07
htotal_int = 0
np = 3
Qt = 3.875/1000
Q1 = Qt/np
mp = Qt*rho
g = 9.81
d = 9
Leq1 = 5
Leq2 = 4
Leq3 = 6
L_head  = 150
H = 75
# calcul de h1
A1 = pi*D1**2/4.0
u1 = Q1/A1
print ('u1 = ',u1)
Re1a = D1*u1/nu
Re1 = 4*Q1/(pi*D1*nu)
f1 = Colebrook(Re1,0.0)
print (Re1,f1)
Lt1 = 2*H + d + Leq1
h1 = f1*(Lt1/D1)*u1**2/(2*g) #  en metres
# calcul de h1
A2 = pi*D2**2/4.0
u2 = 2*Q1/A2
Re2 = D2*u2/nu
f2 = Colebrook(Re2,0.0)
Lt2 = d + Leq2
h2 = f2*(Lt2/D2)*u2**2/(2*g) #  en metres
# calcul de h3
A3 = pi*D3**2/4.0
u3 = Qt/A3
Re3 = D3*u3/nu
f3 = Colebrook(Re3,0.0)
Lt3 = L_head + Leq3
h3 = f3*(Lt3/D3)*u3**2/(2*g) #  en metres
htotal  = h1 + h2 + h3
print ('htotal = ',htotal,'m ')
print ('htotal = ',ft_m(htotal),' ft ')
print('gpm  = ',gpm_m3s(Qt),' gpm')
Wp = mp*g*htotal/(0.74*0.83)
print ('Wp = ',Wp)

