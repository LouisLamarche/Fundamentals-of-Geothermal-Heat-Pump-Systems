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

def Calcul_R_Claesson(z,x,mu,rp,Rp,ks):

    R1 = 1/(2*pi*ks)*np.log(2*z/rp) + Rp
    R2 = R1
    R12 = 1/(2*pi*ks)*np.log(np.sqrt(x**2+4*z**2)/x)
    Rser = 2*(R1*R2 - R12**2 + mu**2/4)/(R1 + R2 - 2*R12)

    R1d = (R1*R2 - R12**2)/(R2 - R12)
    R2d = (R1*R2 - R12**2)/(R1 - R12)
    R12d = (R1*R2 - R12**2)/(R12)
    k1 = mu/R1d
    k2 = mu/R2d
    k12 = mu/R12d
    km = (k1+k2)/2
    eta = np.sqrt(km**2 + 2*k12*km)
    Req = 2*(R1d*R2d)/(R1d+R2d)*eta/np.tanh(eta)
    return Rser,Req
do = 0.03
di = 0.025
ks = 1.2
alj = 0.04
rho = 1000
g = 9.81
q_bat = -10000
rend = 0.3
To = 12.0
tshift = 5     #
amp = 14      #
tcl = 186         #  jour le plus chaud
z = 1.
x = 0.75
tf = 20000
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
kf = 0.61
rhof = 1000
Prf = 5.86
nuf = 8.57e-7
muf = nuf*rhof
Cpf = 4180


ncc = len(q_cond)
COPs1 = np.zeros(ncc)
COPs2 = np.zeros(ncc)
COPs3 = np.zeros(ncc)
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
delta = 0.5
for i in range(0,ncc):
    Vp = q_cond[i]
    mp = Vp*rhof
    Re = 4*mp/(pi*di*muf)
    if (Re>2300.0):
        # Gnielienski
        f = (0.79*np.log(Re)- 1.64)**-2
        Nud=((Re-1000.)*f*Prf/8.)/(1.+12.7*np.sqrt(f/8.)*(Prf**(2./3.)-1))
    else:
        Nud = 3.6
        print('Careful laminar')
    hf = (Nud*kf)/di
    Rconv = 1/(pi*di*hf)
    Rp[i] = Rconv + Rcond
    COP = COPv[i]
    q_sol = q_bat*(COP+1)/(COP)
    CC = mp*Cpf
    Tfi = Tfo - q_sol/CC
    Tf_cl = (Tfi+Tfo)/2
    print(Tfi,Tfo)
    Xf = rp/(2*np.sqrt(alhr*tf))
    Rs = Calcul_R_horizontal(Xf,zt,xt,ks)
    u[i] = q_cond[i]/Ac
    L[i] = (q_sol*(Rs+Rp[i]))/(Tocl - Tf_cl)
    Ltr = L[i]/2
    ok = False
    while not ok:
        mu = Ltr/CC
        Rs2,Rs3 = Calcul_R_Claesson(z,x,mu,rp,Rp[i],ks)
        Ln = (q_sol*Rs2)/(Tocl - Tf_cl)
        errr = abs(Ln - Ltr)
        if errr < delta:
            ok = True
        else:
            Ltr = Ln
    ok = False
    Ltr2 = Ltr
    while not ok:
        mu = Ltr2/CC
        Rs2,Rs3 = Calcul_R_Claesson(z,x,mu,rp,Rp[i],ks)
        Ln = (q_sol*Rs3)/(Tocl - Tf_cl)
        errr = abs(Ln - Ltr2)
        if errr < delta:
            ok = True
        else:
            Ltr2 = Ln
    h_g[i] = f*(L[i]/di)*u[i]**2/(2*g)
    h_g2 = f*(Ltr/di)*u[i]**2/(2*g)
    h_g3 = f*(Ltr2/di)*u[i]**2/(2*g)
    WP[i] = mp*g*(h_condv[i]+h_g[i])/rend
    WP2 = mp*g*(h_condv[i]+h_g2)/rend
    WP3 = mp*g*(h_condv[i]+h_g3)/rend
    Wcomp[i] = -q_bat/COP
    COPs1[i] = -q_bat/(Wcomp[i]+WP[i])
    COPs2[i] = -q_bat/(Wcomp[i]+WP2)
    COPs3[i] = -q_bat/(Wcomp[i]+WP3)
print(COPs1)
print(COPs2)
print(COPs3)
#print(WP)
plot(q_ls,COPs1,q_ls,COPs2,q_ls,COPs3)
show()