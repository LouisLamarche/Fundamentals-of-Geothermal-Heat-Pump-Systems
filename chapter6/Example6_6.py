#
# Exemple 6.6 written by Louis lamarche 29 septemeber 2017
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
H = 100.0
zo = 0
zot = zo/H
tc = H**2/(9*als)   # characteristic time days
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
            tb = tf/tc
            rr = dist/H
            gxx  = g_function_fls(tb,rr,zot)
            Tpp[i] = Tpp[i] + qp/(2*pi*ks)*gxx
Tp = np.mean(Tpp)
Tp2 = Tp_fls(nx,ny,d,qp,als,ks,n_years,H,zot)  # alpha must be in  m2/day , k in W/mK
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
            tb = tf/tc
            rr = dist/H
            gxx  = g_function_fls(tb,rr,zot)
            Tpp[i] = Tpp[i] + qp/(2*pi*ks)*gxx
Tp = np.mean(Tpp)
Tp2 = Tp_fls(nx,ny,d,qp,als,ks,n_years,H,zot)  # alpha must be in  m2/day , k in W/mK
print (Tp,Tp2)

zo = 4
zot = zo/H
tc = H**2/(9*als)   # characteristic time days
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
            tb = tf/tc
            rr = dist/H
            gxx  = g_function_fls(tb,rr,zot)
            Tpp[i] = Tpp[i] + qp/(2*pi*ks)*gxx
Tp = np.mean(Tpp)
Tp2 = Tp_fls(nx,ny,d,qp,als,ks,n_years,H,zot)  # alpha must be in  m2/day , k in W/mK
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
            tb = tf/tc
            rr = dist/H
            gxx  = g_function_fls(tb,rr,zot)
            Tpp[i] = Tpp[i] + qp/(2*pi*ks)*gxx
Tp = np.mean(Tpp)
Tp2 = Tp_fls(nx,ny,d,qp,als,ks,n_years,H,zot)  # alpha must be in  m2/day , k in W/mK
print (Tp,Tp2)
