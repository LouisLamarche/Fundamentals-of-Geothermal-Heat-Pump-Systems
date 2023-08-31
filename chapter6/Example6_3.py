#coding: utf-8
#
# Exemple 6.3 written by Louis lamarche 29 septemeber 2017
#
from geothermal_md import *
import numpy as np
#
# data
#
als = 0.1 # m2/jr
ks = 2.0
alhr = als/24.0
qp = 8.0
n_rings = 3
d = 6.0
ri = d/2.0
dr = 1.5
n_years = 10
tf = n_years*365
Tps = 0
r1 = ri
for i in range(0,n_rings):
    r2 = r1+dr
    rm = (r1+r2)/2.0
    X  = rm/(2.0*np.sqrt(als*tf))
    Ix  = I_function(X)
    Tm = qp/(2*pi*ks)*Ix
    Tps = Tps+Tm*pi*(r2**2-r1**2)/d**2
    print ('Tps = ',Tps)
    r1 = r2
