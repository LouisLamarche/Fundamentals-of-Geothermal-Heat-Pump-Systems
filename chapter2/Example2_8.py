import  numpy as np
# example 2.3
patm = 101.325*1000.0
Ti = 22.0
To = -10.0
k1 = 2.1
k2 = 0.04
ho = 25
hi = 8
L1 = 0.24
L2 = 0.16
A = 12
R1 = 1/(A*ho)
R2 = L1/(A*k1)
R3 = L2/(A*k2)
R4 = 1/(A*hi)
Rt = (R1+R2+R3+R4)
UA = 1/Rt
q = UA*(Ti-To)
T1 = Ti - q*R4
T2 = T1 - q*R3
T3 = T2 - q*R2
T4 = T3 - q*R1
print (' q = ' + '%.2f' % q + ' W')
print ('T1 = ' + '%.2f' % T1 + ' C')
print ('T2 = ' + '%.2f' % T2 + ' C')
print ('T3 = ' + '%.2f' % T3 + ' C')
print ('T4 = ' + '%.2f' % T4 + ' C')
