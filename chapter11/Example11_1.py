#coding: utf-8
#
#
import numpy as np
from geothermal_md import *
from  matplotlib.pyplot import *
#
q  = 100.0      #  m3/hr
q = q/3600.0    # m3/s
b = 20.0        # thickness
rw = 0.5/2.0
r1 = 15.0
r2 = 45.0
s1 = 1.7
s2 = 0.6
K = q/(2*pi*b*(s1-s2))*np.log(r2/r1)
Kj =   K*3600*24
Kh =   K*3600
print ('K = ',K,'m/s')
print ('K = ',Kj,'m/day')
print ('K = ',Kh,'m/hr')
T = K*b
Tj = Kj*b
Th = Kh*b
print ('T = ',T,'m2/s')
print ('T = ',Tj,'m2/day')
print ('T = ',Th,'m2/hr')
sw = s2 + q/(2*pi*T)*np.log(r2/rw)
print ('sw  = ',sw,' m')
sw = s1 + q/(2*pi*T)*np.log(r1/rw)
print ('sw  = ',sw,' m')
ro = np.exp(sw*2*pi*T/q)*rw
print ('ro  = ',ro,' m')
ro = np.exp(s1*2*pi*T/q)*r1
print ('ro  = ',ro,' m')
ro = np.exp(s2*2*pi*T/q)*r2
print ('ro  = ',ro,' m')
#
#
