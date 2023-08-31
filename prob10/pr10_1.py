#coding: latin-1
# Exemple 4.4
from geothermal_md import *
from hydraulic_md import *
import  numpy as np
from CoolProp.CoolProp import *
from  matplotlib.pyplot import *

#
#
def Calcul_R_horizontal(X1,zt1,xt2,ks):
    X2 = xt2*X1
    X3 = 2*zt1*X1
    dd = np.sqrt((2*zt1)**2 + xt2**2)
    X4 = dd*X1
    R = ((I_function(X1)+I_function(X2))-(I_function(X3)+I_function(X4)))/(2*pi*ks)     # mK/W
    return R
do = 0.03
di = 0.025
ks = 1.2
alj = 0.04
rho = 1000
g = 9.81
q_bat = -10000
rend = 0.3
To = 12.0
tshift = 5     # Valeur pour Montreal ( en jour)
amp = 14      # amplitude de variation de la temperature dans l'annéee
tcl = 186         #  jour le plus chaud
z = 1.
x = 0.75
tf = 120
Tfo = 30

#q_ls = np.array([0.3,0.4,0.5])
#COPv = np.array([4.5,4.75,5.0])

COPv = np.array([4.,4.4,4.7,5.0,5.2])
q_ls = np.array([0.3,0.35,0.4,0.45,0.5])

Vnom = 0.4  # l/s
DPac = 10.1*(q_ls/Vnom)**2
q_cond = q_ls/1000  #
h_condv = DPac*1000/(rho*g)
fluid = 'water'
Trefk = 300
patm = 101.325*1000.0
mu = PropsSI('viscosity','T',Trefk,'P',patm,fluid)
kf = PropsSI('conductivity','T',Trefk,'P',patm,fluid)
rho = PropsSI('D','T',Trefk,'P',patm,fluid)
Pr = PropsSI('Prandtl','T',Trefk,'P',patm,fluid)
nu = mu/rho
Cp = 4180

ncc = len(q_cond)
COPs = np.zeros(ncc)
L = np.zeros(ncc)
Wcomp = np.zeros(ncc)
WP = np.zeros(ncc)
Rp = np.zeros(ncc)
h_g = np.zeros(ncc)
DPacu = np.zeros(ncc)
u = np.zeros(ncc)
Ac = pi*di**2/4
Rcond = np.log(do/di)/(2*pi*0.4)
dTcl,dTmax = Delta_T_earth(z,tcl,tshift,alj,amp)
Tocl = To + dTcl
rp = do/2
zt = z/rp
xt = x/rp
alhr = alj/24

for i in range(0,ncc):
    Vp = q_cond[i]
    mp = Vp*rho
    Re = 4*mp/(pi*di*mu)
    if (Re>2300.0):
        # Gnielienski
        f = (0.79*np.log(Re)- 1.64)**-2
        Nud=((Re-1000.)*f*Pr/8.)/(1.+12.7*np.sqrt(f/8.)*(Pr**(2./3.)-1))
    else:
        Nud = 3.6
        print('Careful laminar')
    hf = (Nud*kf)/di
    Rconv = 1/(pi*di*hf)
    Rp[i] = Rconv + Rcond
    COP = COPv[i]
    q_sol = q_bat*(COP+1)/(COP)
    CC = mp*Cp
    Tfi = Tfo - q_sol/CC
    Tf_cl = (Tfi+Tfo)/2
    print(Tfi,Tfo)
    Xf = rp/(2*np.sqrt(alhr*tf))
    Rs = Calcul_R_horizontal(Xf,zt,xt,ks)
    u[i] = q_cond[i]/Ac
    L[i] = (q_sol*(Rs+Rp[i]))/(Tocl - Tf_cl)
    h_g[i] = f*(L[i]/di)*u[i]**2/(2*g)
    WP[i] = mp*g*(h_condv[i]+h_g[i])/rend
    Wcomp[i] = -q_bat/COP
    COPs[i] = -q_bat/(Wcomp[i]+WP[i])
print(COPs)
print(WP)
plot(q_ls,COPs)
show()