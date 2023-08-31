
#
# Exemple 6.5 written by Louis lamarche 29 septemeber 2017
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
d = 6.0
n_years = 10
nx = 3  # nombre de rows
ny = 3  # nombre de columns
nt = nx*ny
tf = n_years*365
Tpp = np.zeros(nt)
x = np.zeros(nt)
y = np.zeros(nt)
k = 0
dx = d
dy = d
for i in range(0,nx):
    for j in range(0,ny):
        x[k] = dx*(i-1)
        y[k] = dy*(j-1)
        k = k+1
for i in range (0,nt):
    for j in range(0,nt):
        if(i!=j):
            dist = np.sqrt((x[i]-x[j])**2+(y[i]-y[j])**2)
            xx  = dist/(2*np.sqrt(als*tf))
            Ixx  = I_function(xx)
            Tpp[i] = Tpp[i] + qp/(2*pi*ks)*Ixx
Tp = np.mean(Tpp)
Tp2 = Tp_ils(nx,ny,d,qp,als,ks,n_years)  # alpha must be in  m2/day , k in W/mK
print (Tp,Tp2)
nx = 9  # nombre de rows
ny = 1  # nombre de columns
nt = nx*ny
Tpp = np.zeros(nt)
x = np.zeros(nt)
y = np.zeros(nt)
k = 0
for i in range(0,nx):
    for j in range(0,ny):
        x[k] = dx*(i-1)
        y[k] = dy*(j-1)
        k = k+1
for i in range (0,nt):
    for j in range(0,nt):
        if(i!=j):
            dist = np.sqrt((x[i]-x[j])**2+(y[i]-y[j])**2)
            xx  = dist/(2*np.sqrt(als*tf))
            Ixx  = I_function(xx)
            Tpp[i] = Tpp[i] + qp/(2*pi*ks)*Ixx
Tp = np.mean(Tpp)
Tp2 = Tp_ils(nx,ny,d,qp,als,ks,n_years)  # alpha must be in  m2/day , k in W/mK
print (Tp,Tp2)

