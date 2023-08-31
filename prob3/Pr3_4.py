import numpy as np
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
def bilan(T):
    T1 = T[0]
    T2 = T[1]
    Qcb = calcul_qevap(T1,T2)
    Qhb = calcul_qcond(T1,T2)
    y1 = Qin1 + UA12*(T1-T2) - Qcb
    y2 = Qin + Qhb - UA12*(T1-T2) - UA1*(T1 - Tinf)
    y = np.array([y1,y2])
    return y
Ti = np.array([19,21])
y = bilan(Ti)
res = fsolve(bilan, Ti)
Tsc = res[0]
Tlc = res[1]
T1b = Tsc
T2b = Tlc
print('T1 = ',T1b)
print('T2 = ',T2b)
Qcb = calcul_qevap(Tsc,Tlc)
Qhb = calcul_qcond(Tsc,Tlc)
Qi =  UA12*(T1b - T2b)
print('Q12 = ',Qi)
print('q evap zone 2 ',Qcb, ' kW')
Qout =  UA1*(T1b - Tinf)
print('HEat loss 1 = ',Qout + Qi, ' kW')
print('Heat gain 1 = ',Qhb, ' kW')
print('HEat gain 2 = ',Qi + Qin1, ' kW')
print('HEat loss 2 = ',Qcb, ' kW')
Wc = (Qhb -Qcb)
print('Compressor = ',Wc, ' kW')


