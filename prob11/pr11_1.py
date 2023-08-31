#coding: latin-1
#
# example 11.4
#
import numpy as np
from geothermal_md import *
from conversion_md import *
from  matplotlib.pyplot import *
from scipy.optimize import curve_fit
import  pandas as pd
#
# donnees du probleme
#
gam = 0.5772157
#
#
A = array([[0.20,1.76],
    [0.50 ,  2.75],
    [1.00 ,  3.59],
    [2.00 ,  4.26],
    [5.00 ,  5.28],
    [10.00 ,  5.90],
    [20.00 ,  6.47],
    [50.00 ,  6.92],
    [100.00 ,  7.11],
    [200.00 ,  7.20],
    [500.00 ,  7.21],
    [1000.00 ,  7.21]])
#
t = A[:,0]      # time in minutes
t = t*60
sp = M[:,1]     # drawndown in ft
sf = m_ft(sp)
gpm = 1000
qo = m3s_gpm(gpm)      # heat transfer per metre
rw = m_ft(100)
# first approach
#
Ti = 2.98e-3     # m2/min
Si = 1.1e-5
rbi = 0.2
po = [Ti,Si,rbi]
alm = Ti/Si
G_vect = np.vectorize(leaky_function)
def s_theo(x,Tnew,Snew,rb):
    al = Tnew/Snew
    u = rw**2/(4*al*x)
    s = qo/(4*pi*Tnew)*G_vect(u,rb)
    return s
#
nt = len(t)
ni = 3
tn =t[ni:nt]
sn = sf[ni:nt]
params,resn = curve_fit(s_theo,tn,sn,po)
Tn = params[0]
Tip = Tn*gal_m3()/ft_m()*3600*24
Sn = params[1]
rbn = params[2]
print(Tn,Sn,rbn)
v = rbn/2
pm = 4*Tip*v**2/100**2
y0  = s_theo(t,Tn,Sn,rbn)
u = rw**2/(4*alm*t)
un = rw**2/(4*alm*tn)
print(un[0])
x = 1/u
p1 = loglog(x, sf, label='Measured',color = 'black')
p2 = loglog(x, y0, label='Curve fit',color = 'black',  linestyle='none', marker='o')
ll = legend()
show()
