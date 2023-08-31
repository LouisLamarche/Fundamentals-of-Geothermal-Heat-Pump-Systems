#coding: latin-1
from geothermal_md import *
import numpy as np
from  matplotlib.pyplot import *
from conversion_md import *
from hydraulic_md import *
from CoolProp.CoolProp import *
import pandas as pd
#
#
NaN = np.NaN
Patm = 101.325*1000.0
cas = 2
if cas == 1:
    fluide = 'Water'
elif cas ==2:
    pour = 25
    fluide  = 'INCOMP::APG-25%'   # propylene - glycol 20 %
#
#
g = 9.81
nb = 12
flow = 360 # l/min
q_tot = flow/60/1000 # m3/s
q_loop = q_tot /nb
dist_building = 45.0
socket_fusion = False
H = 140
d = 6
SDR = 13.5
tube_loop = 1.25
tube_conn = 3
D_loop = sdr_pipe(tube_loop,SDR)[0]
Trefk = 5 + 273.15 # choix arbitraire
mu = PropsSI('viscosity', 'T', Trefk, 'P', Patm, fluide)
rho = PropsSI('D', 'T', Trefk, 'P', Patm, fluide)
Cp = PropsSI('Cpmass', 'T', Trefk, 'P', Patm, fluide)
Pr = PropsSI('Prandtl', 'T', Trefk, 'P', Patm, fluide)
kf = PropsSI('conductivity', 'T', Trefk, 'P', Patm, fluide)
rho_eau = PropsSI('D', 'T', Trefk, 'P', Patm, 'Water')
mu_eau = PropsSI('viscosity', 'T', Trefk, 'P', Patm, 'Water')
Sg = rho/rho_eau
nu = mu/rho     # viscosité cinématique
nu_eau = mu_eau/rho_eau     # viscosité cinématique
mp_total = q_tot*rho
#
#
#Calcul des pertes de charges
#
# boucle geohermique
#
# 4) pertes de charges dans la boucle
#
data =  np.zeros((13,5))

if socket_fusion:
    Leq_ubend = sing_head_loss('Socket U-bend',tube_loop)
else:
    Leq_ubend = sing_head_loss('Butt U-bend',tube_loop)
A_loop = pi*D_loop**2/4
u_loop = q_loop/A_loop
Re_loop = D_loop*u_loop/nu;
f_loop = Colebrook(Re_loop,0.0)
L_loop = 2*H
Lt_loop = L_loop + Leq_ubend
#
#
# boucle seule
#
h_loop = f_loop*(Lt_loop/D_loop)*u_loop**2/(2*g)  # % Dp en metres
data[0,0] = Re_loop
data[0,1] = f_loop
data[0,2] = D_loop
data[0,3] = Lt_loop
data[0,4] = h_loop

#
# rajout des inter connexions
#
tube = [tube_loop,1.5,1.5,2,2,3,3,3,3,3,3]
htotal_loop = 0
hinv = 0
Ll = d*np.ones(nb-1)
Ll[7] = 3*d
hv = np.zeros(nb-1)
for i in range(0,nb-1):
    i1 = i+1.0
    x = i1/nb
    debiti = q_tot*x          #  débit qui passe dans la section i
    D =  sdr_pipe(tube[i],SDR)[0]
    A = pi*D**2/4.0
    u = debiti/A
    Re = D*u/nu
    f = Colebrook(Re,0.0)
    h = f*(Ll[i]/D)*u**2/(2*g) #  en metres
    hv[i] = h
    data[i+1,0] = Re
    data[i+1,1] = f
    data[i+1,2] = D
    data[i+1,3] = Ll[i]
    data[i+1,4] = h
    htotal_loop = htotal_loop + h
ht_loop = h_loop + htotal_loop
#
# rajout de la connextion 3 " au batiment
L_bat = 2*dist_building
D_bat =  sdr_pipe(tube_conn,SDR)[0]
A_bat = pi*D_bat**2/4.0
u_bat = q_tot/A_bat
Re_bat = D_bat*u_bat/nu
f_bat = Colebrook(Re_bat,0.0)
h_bat = f_bat*(L_bat/D_bat)*u_bat**2/(2*g) #  en metres
h_tot = ht_loop + h_bat
data[12,0] = NaN
data[12,1] = NaN
data[12,2] = NaN
data[12,3] = NaN
data[12,4] = ht_loop
#
# rajout du reste supposé le même pour toutes les configurations
#

W_fluid = mp_total*g*h_tot
print ('The head losses are  =',h_tot,' m')
print ('fluid power is  =',W_fluid,' W')
#

indd  = ['loop','1','2','3','4','5','6','7','8','9','10','11','total']
coll = ['Re','f','D','L','h']
def f1(x):
    return '%.0f' % x
def f2(x):
    return '%.3f' % x
def f3(x):
    return '%.2f' % x
df = pd.DataFrame(data,columns = coll,index = indd)
print(df.to_string(index=True,formatters=[f1,f2,f2,f1,f3]).replace('NaN','   '))
