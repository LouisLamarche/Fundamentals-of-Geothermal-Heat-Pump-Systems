#coding: utf-8
#
# exemple 11.6
#
import  numpy as np
from geothermal_mod import *
from conversion_mod import *
from  matplotlib.pyplot import *
from heat_exchanger_md import *

Cpf = 4200
muf= 0.00138
Prf = 10
kf = 0.57
rhof = 1000
g = 9.8


def Calcul_h_ann(di,do,mp):
    Ac = pi*(do**2-di**2)/4
    u = mp/(rhof*Ac)
    dh = do-di
    Re = rhof*u*dh/muf
    a = di/do
    Res = Re*((1+a**2)*np.log(a) + (1 - a**2))/((1-a)**2*np.log(a))
    f = (0.781*np.log(Res)-1.5)**(-2)
    Fa = 0.75*a**-(0.17)
    # Gnieleski 2015
    Nu = f/8*(Re-1000)*Prf/(1 + 12.7*np.sqrt(f/8)*(Prf**(2/3) - 1))*Fa
    hh = Nu*kf/dh
    return hh

def Calcul_h_circ(di,mp):
    Ac = pi*di**2/4
    u = mp/(rhof*Ac)
    Re = rhof*u*di/muf
    f = (0.79*np.log(Re)-1.64)**(-2)
    # Gnieleski 2015
    Nu = f/8*(Re-1000)*Prf/(1 + 12.7*np.sqrt(f/8)*(Prf**(2/3) - 1))
    hh = Nu*kf/di
    return hh
#
Tch = np.array([-1.0,3.0,6.0,7.0,9.0,12.5])
Tcl = np.array([10.0,15.0,20,23,27,30])
COPchv = np.array([3.3,3.9,4.4,4.5,4.7,4.9])
COPclv = np.array([8.,7.,5.5,5,4.5,4.3])

q_ch = 12000.
q_cl = 16000.
Tsi = 9.0   # temperature eau sol
Vprch = 50/60/1000
Vpnom = Vprch
mpch = Vprch*rhof    # debit  kg/s si o a lde l'eau)
Vprcl = 41.5/60/1000
mpcl = Vprcl*rhof    # debit  kg/s si o a lde l'eau)
Cch = mpch*Cpf          # en kW/K
Ccl = mpcl*Cpf          # en kW/K
di = 0.021
do = 0.03
Tcod = 7.0
Thod = 15.0
Tho_min = 4.0
eff = 0.8
COPcl = np.interp(Thod,Tcl,COPclv)
COPch = np.interp(Tcod,Tch,COPchv)
#
# convection coef heat pump side
#
h_load = Calcul_h_circ(di,mpch)
#
# heating
#
Thi = Tsi
Tco = Tcod
Cc = Cch
qech = q_ch*(COPch-1)/COPch
Tci = Tco - qech/Cc
qmax = qech/eff
Cmin = qmax/(Thi-Tci)
if Cmin < Cc:
    Ch = Cmin
else:  # impossible to have efficiency = 0.8
    Cmin = Cc
    qmax = Cmin*(Thi - Tci)
    eff = qech/qmax
    print( ' The new efficiency is = ',eff)
    Ch = Cc
Cmax = Cc
Tho = Thi - qech/Ch
mp_ch = Ch/Cpf
print('Tout (heating) = ',Tho)
if Tho < Tho_min:
    print(' Temperature too low')
print('flow rate (heating mode) = ',mp_ch,' kg/s')
h_source = Calcul_h_ann(di,do,mp_ch)
U_ch = h_load*h_source/(h_load+h_source)
Cr = Cmin/Cmax
NTU_ch = counter_flow_NTU(eff,Cr)
UA = NTU_ch*Cmin
A_ch = UA/U_ch
print('A (heating mode) = ',A_ch, ' m2')



# analyse de la puissance de pompe pour le chauffage uniquement
beta = 1000.0   #
C = 1*3600     #
Ccond = 1.4e7
Q = mp_ch/rhof
st = beta*Q + C*Q**2
ho = 10.0
hp = Ccond*Q**2
ht = ho + st + hp
mgh = mp_ch*9.8*ht
DeltaP = 13*1000*(Vprch/Vpnom)**2  # Pa
hhp = DeltaP/g/rhof
mghp = mpch*g*hhp
rend1 = 0.35
rend2 = 0.55
Wcomp = q_ch/COPch
Wpompe = (mgh)/rend1 + mghp/rend2
print (' h total = ',ht,' m')
print (' W pompe = ',Wpompe,' W')
COPsys = q_ch/(Wcomp+Wpompe)
print (' COP systeme = ',COPsys)
