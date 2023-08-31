#
# Example 4.4 written by Louis Lamarche 22 sept 2017
#
import numpy as np
from geothermal_md  import *

alj = 0.1       # m2/day
alhr = alj/24.0 # m2/hr
rb = 0.075
ks = 2.0
To = 10.0
H = 90.0
zo = 3.6
B = 6.0
Bt = B/H
rr  = rb/H
zot = zo/H
ht = H/rb
qp = 3.0
# intermediate time in hours
t1 = 20.0*365
tc = H**2/(9*alj)
tt = t1/tc
nx = 2
ny = 2
nb = nx*ny
zb = np.zeros([nb,2])
ib = 0
gi1 = g_function_fls(tt,rbb = rr,zob = zot)
gi2 = g_function_fls(tt,rbb = Bt,zob = zot)
gi3 = g_function_fls(tt,rbb = Bt*np.sqrt(2),zob = zot)
gi = (gi1+2*gi2+gi3)
for i in range(0,nx):
    for j in range(0,ny):
        zb[ib] = [i*Bt,j*Bt]
        ib = ib+1
zb2 = np.zeros([1,2])
ib = 0
for i in range(0,1):
    for j in range(0,1):
        zb2[ib] = [i*Bt,j*Bt]
        ib = ib+1
g = compute_g_function(zb,tt,rbb = rr,zob = zot)
print ('g = ',g,gi)

DT = -qp/(2*pi*ks)*g
print ('Delta T = ',DT)
#
# interpolation
#
ib = 0
Bt = 0.05
for i in range(0,nx):
    for j in range(0,ny):
        zb[ib] = [i*Bt,j*Bt]
        ib = ib+1
go = compute_g_function(zb,tt)  # rb = 0.0005, zob = 0.0 default values
ib = 0
Bt = 0.1
for i in range(0,nx):
    for j in range(0,ny):
        zb[ib] = [i*Bt,j*Bt]
        ib = ib+1
g1 = compute_g_function(zb,tt)  # rb = 0.0005, zob = 0.0 default values
goc = go - np.log(rr/0.0005)
g1c = g1 - np.log(rr/0.0005)
gn = goc + (g1c - goc)*(0.067 - 0.05)/0.05
print ('g (approx) = ',gn)
DTn = -qp/(2*pi*ks)*gn
print ('Delta T (approx) = ',DTn)