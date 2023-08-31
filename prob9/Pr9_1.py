#coding: latin-1
import numpy as np
#
#
pi = np.pi
gamm = np.euler_gamma
qp = 66
rb = 0.126/2.0
To = 12.0
CC = 1.92e6
m1 = 2.
k1 = qp/(4*pi*m1)
print ('k ( rest ) = ',k1)
p1 = [5.0,25]
p2 = [7.25,30]
m2 = (p1[1]-p2[1])/(p1[0]-p2[0])
k2 = qp/(4*pi*m2)
print ('k ( inj ) = ',k2)
b = p1[1] - m2*p1[0]
als = k1/CC
alm = als*60
Rb = (b-To)/qp - (np.log(4*alm/rb**2)-gamm)/(4*pi*k1)
print ('Rb = ',Rb)

