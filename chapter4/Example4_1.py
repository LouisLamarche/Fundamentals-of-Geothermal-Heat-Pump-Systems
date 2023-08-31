#coding: utf-8
#
# Example 4.1 written by Louis Lamarche 22 sept 2017
#
import numpy as np
from geothermal_md  import *
#
# data
#
qp = -40.0  # heat flux W/m ( negative is heat rejection)
t = 250.0   # time in days
r = 2.0     # distance (m)
rb = 0.01   # borehole radius (m)
alj = 0.1   # diffusivity in m2/day
ks  = 2.5   # soil conductivity W/mK
#
# fct I , ils 1
X = r/(2*np.sqrt(alj*t))
print ('X = ' + '%.2f' % X)
I = I_function(X)
print ('I = ' + '%.2f' % I)
DT = - qp*I/(2*pi*ks)
print ('DT (ils) = ' + '%.4f' % DT)
# ils2
Fo = alj*t/rb**2
rt = r/rb
print ('Fo = ' + '%.2f' % Fo)
G1 = G_function_ils(Fo,rt)
print ('Gils = ' + '%.4f' % G1)
DT1 = - qp*G1/ks
print ('DT (ils) = ' + '%.4f' % DT1)
# ics
G2 = G_function_ics(Fo,rt)
print ('Gcls = ' + '%.4f' % G2)
DT2 = - qp*G2/ks
print ('DT (ics) = ' + '%.4f' % DT2)
# ils3
Fos  = alj*t/r**2
print ('Fos = ' + '%.2f' % Fos)
G3 = G_function_ils(Fos)  # rt = 1 by default
print ('Gils = ' + '%.4f' % G3)
DT3 = - qp*G3/ks
print ('DT (ils) = ' + '%.4f' % DT3)

