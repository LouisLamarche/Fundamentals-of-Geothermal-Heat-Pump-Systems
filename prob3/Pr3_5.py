import  numpy as np
from scipy.optimize import curve_fit,newton,minimize,fsolve

def calcul_qevap(T1,T2):
    y =  11.02 -0.033*T1 - 0.164*T2 + 0.012*T1*T2 - 4.862e-5*T2**2
    return y
def calcul_qcond(T1,T2):
    y =  12.38 - 0.04*T1 - 0.136*T2 + 0.012*T1*T2 - 4.43e-5*T2**2
    return y

Tinf = 0
UA12 = .4
UA1 = 0.6
Qin = 0
Qin1 = 12
ep1 = 0.9
ep2 = 0.9
Ce1 = 0.5*4.2
Ce2 = Ce1
Ca1 = 0.9*1.007
Ca2 = Ca1
Cmin1 = min(Ce1,Ca1)
Cmin2 = min(Ce2,Ca2)
def bilan(T):
    T1 = T[0]
    T2 = T[1]
    Tsc = T[2]
    Tlc = T[3]
    Qcb = calcul_qevap(Tsc,Tlc)
    Qhb = calcul_qcond(Tsc,Tlc)
    y1 = Qin1 + UA12*(T1-T2) - Qcb
    y2 = Qin + Qhb - UA12*(T1-T2) - UA1*(T1 - Tinf)
    y3 = Qhb*(-1/Ce1 + 1/(ep2*Cmin2)) + T1 - Tsc
    y4 = T2 + Qcb*(1/Ce2 - 1/(ep2*Cmin2)) - Tlc
    y = np.array([y1,y2,y3,y4])
    return y
Ti = np.array([19,21,25,10])
y = bilan(Ti)
res = fsolve(bilan, Ti)
print(res)
y = bilan(res)
T1 = res[0]
T2 = res[1]
Tsc = res[2]
Tlc = res[3]
print('T1 = ',T1)
print('T2 = ',T2)
Qcb = calcul_qevap(Tsc,Tlc)
Qhb = calcul_qcond(Tsc,Tlc)
Thi = Tsc + Qhb/Ce1
Qhbb = ep1*Cmin1*(Thi - T1)
Tci = Tlc - Qcb/Ce2
#print('T (water out cond) = ',Thi)
#print('T (water out evap) = ',Tci)
Qcbb = ep2*Cmin2*(T2 - Tci)
Qi =  UA12*(T1 - T2)
print('Q12 = ',Qi)
print('Heat from evaporator',Qcb)
Qout =  UA1*(T1 - Tinf)
print('HEat loss 1 = ',Qout + Qi)
print('Heat gain 1 = ',Qhb)
print('HEat gain 2 = ',Qi + Qin1)
print('HEat loss 2 = ',Qcb)
Wc = (Qhb -Qcb)


