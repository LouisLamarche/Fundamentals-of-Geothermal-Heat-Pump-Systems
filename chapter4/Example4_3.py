import numpy as np
from geothermal_md import *
#
# Example 4.4 written by Louis Lamarche 12 june 2018
#
ks = 2.1  # W/m K
Cs =  2.2e6   # J/m3 K
als = ks/Cs  #  diffusivity in m2/s
alhr = als*3600.0 #  diffusivity in m2/hr
tf = 10*24      # final time (hours)
qp1 = 10.0
qp2 = -8.0
qp3 = 12.0
#
x = np.array([1.0,1.0])
x1 = np.array([0.0,0.0])
x2 = np.array([6.0,0.0])
x3 = np.array([12.0,0.0])
d1 = np.linalg.norm(x-x1)   # distance from well 1
d2 = np.linalg.norm(x-x2)  # distance from well 2
d3 = np.linalg.norm(x-x3)  # distance from well 3
X1  = d1/(2.0*np.sqrt(alhr*tf))
X2  = d2/(2.0*np.sqrt(alhr*tf))
X3  = d3/(2.0*np.sqrt(alhr*tf))
I1  = I_function(X1)
I2  = I_function(X2)
I3  = I_function(X3)
DT1 = -qp1/(2*pi*ks)*I1
DT2 = -qp2/(2*pi*ks)*I2
DT3 = -qp3/(2*pi*ks)*I3
DT = DT1+DT2+DT3
print  (' Delta T1 = ',DT1)
print  (' Delta T2 = ',DT2)
print  (' Delta T3 = ',DT3)
print  (' Delta T (total)  = ',DT)
