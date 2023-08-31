#coding: utf-8
#
# exemple 11.6
#
import  numpy as np
from geothermal_mod import *
from conversion_mod import *
from  matplotlib.pyplot import *
#
Tch = np.array([-1.0,3.0,6.0,7.0,9.0,12.5])
Tcl = np.array([10.0,15.0,20,23,27,30])
COPch = np.array([3.3,3.9,4.4,4.5,4.7,4.9])
COPcl = np.array([8.,7.,5.5,5,4.5,4.3])
pch = np.polyfit(Tch,COPch,2)
tch = np.arange(-1,13,0.5)
copc = np.polyval(pch,tch)
pcl = np.polyfit(Tcl,COPcl,2)
tcl = np.arange(10,31,1.0)
copc = np.polyval(pch,tch)
copl = np.polyval(pcl,tcl)

def calcul_eff(NTU,Cr):
    if Cr<1:
        eff = (1-np.exp(-NTU*(1-Cr)))/(1-Cr*np.exp(-NTU*(1-Cr)))
    else:
        eff = NTU/(1+NTU)
    return eff
q_ch = 12000.
Tsi = 9.0   # temperature eau sol
rho = 1000
mpch = 50/60/1000*rho    # debit  kg/s si o a lde l'eau)
Cch = mpch*4200.          # en kW/K
# chauffage
Thi = Tsi
Cc = Cch
Tho_min = 4.0
# hypothese
U = 3000.0
eff = 0.8
# hyp eff = 0.8
Tco = 7.0
COPchv = 4.5
Tho = 15.0
COPclv = 7
Wcomp = q_ch/COPchv
qech = q_ch*(COPchv-1)/COPchv
Tci = Tco - qech/Cc
qmax = qech/eff
Cmin = qmax/(Thi-Tci)
if Cmin < Cc:
    Ch = Cmin
else:
    Cmin = Cc
    qmax = Cmin*(Thi - Tci)
    eff = qech/qmax
    Ch = Cc
Cmax = Cc
Tho = Thi - qech/Ch
mp_ch = Ch/4200.0
print('Tout = ',Tho)
print('flow rate = ',mp_ch,' kg/s')
Cr = Cmin/Cmax
if Cr < 1:
    NTU = 1/(Cr-1)*np.log((eff-1)/(eff*Cr-1))
else:
    NTU = eff/(1.0-eff)
UA = NTU*Cmin
A_ch = UA/U
print(A_ch)

# analyse de la puissance de pompe pour le chauffage uniquement
beta = 1000.0   #
C = 1*3600     #
Ccond = 1.4e7
Q = mp_ch/rho
st = beta*Q + C*Q**2
ho = 10.0
hp = Ccond*Q**2
ht = ho + st + hp
mgh = mp_ch*9.8*ht
rend = 0.35
Wpompe = mgh/rend
print (' h total = ',ht,' m')
print (' W pompe = ',Wpompe,' W')
COPsys = q_ch/(Wcomp+Wpompe)
print (' COP systeme = ',COPsys)

Cc = Cc/2
q_ch = q_ch/2
Cmin = min(Ch,Cc)
Cmax = max(Cc,Ch)
Cr = Cmin/Cmax
NTU = U*A_ch/Cmin
eff = calcul_eff(NTU,Cr)


def fct(Tco):
    COPchn = np.interp(Tco,Tch,COPch)
    qech = q_ch*(COPchn-1)/COPchn
    Tci = Tco - qech/Cc
    qmax = Cmin*(Thi-Tci)
    q = eff*qmax
    return q - qech

q = newton(fct,6)


