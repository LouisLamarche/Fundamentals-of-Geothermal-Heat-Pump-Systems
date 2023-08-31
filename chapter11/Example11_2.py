#coding: utf-8
#
# exemple 11-2
import numpy as np
from geothermal_md import *
from  matplotlib.pyplot import *
#
q1  = 100.0     #  m3/hr
q = q1/3600.0
b = 20.0
rw = 0.5/2.0
r1 = 15.0
r2 = 45.0
s1 = 1.7
s2 = 0.6
K = q/(2*pi*b*(s1-s2))*np.log(r2/r1)
Kj =   K*3600*24
print ('K = ',Kj,'m/jr')
T = K*b
Tj = Kj*b
print ('T = ',T,'m2/s')
print ('T = ',Tj,'m2/jr')
sw = s2 + q/(2*pi*T)*np.log(r2/rw)
print ('sw  = ',sw,' m')
ro = np.exp(sw*2*pi*T/q)*rw
print ('ro  = ',ro,' m')
C = 0.8*60*60
beta = np.log(ro/rw)/(2*pi*T)
print ('beta = ',beta,'s/m2')
st = 8.0
C1 = (st - sw)/q1**2
print ('C = ',C1,' hr2/m5')
C = (st - sw)/q**2
print ('C = ',C,' s2/m5')
Cm = C/3600
print ('C = ',Cm,' min2/m5')
st = sw + C*q*q
print ('st  = ',st,' m')
st2 = beta*q + C*q*q
print ('st  = ',st2,' m')

#
#
