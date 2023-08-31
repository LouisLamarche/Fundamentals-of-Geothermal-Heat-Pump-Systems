#
# Exemple 6.4 written by Louis lamarche 29 septemeber 2017
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
    r1 = r2
nx = 3  # nombre de rows
ny = 3  # nombre de columns
Nt = nx*ny
N4 = 1
N3 = 4
N2 = 4
N1 = 0
Tp = Tps*(N4 + N3*0.75+N2*0.5+N1*0.25)/Nt
Tp2 = Tp_ashrae(nx,ny,d,qp,als,ks,n_years,n_rings)  # alpha must be in  m2/day , k in W/mK
print (Tp,Tp2)
nx = 9  # nombre de rows
ny = 1  # nombre de columns
Nt = nx*ny
N4 = 0
N3 = 0
N2 = 7
N1 = 2
Tp = Tps*(N4 + N3*0.75+N2*0.5+N1*0.25)/Nt
Tp2 = Tp_ashrae(nx,ny,d,qp,als,ks,n_years,n_rings)  # alpha must be in  m2/day , k in W/mK
print (Tp,Tp2)

