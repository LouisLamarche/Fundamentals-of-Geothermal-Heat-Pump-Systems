#coding: utf-8
#
# example 11.4
#
import numpy as np
from geothermal_md import *
from  matplotlib.pyplot import *
from scipy.optimize import curve_fit
#
# donnees du probleme
#
gam = 0.5772157
M = np.loadtxt("..\\data\\pumping_test2.txt")
t = M[:,0]      # time in minutes
sf = M[:,1]     # drawndown in meters
nt = len(sf)
qo = 17      # heat transfer per metre
rw = 12.2
#
#
#
# first approach
#
Ti = 2.4     # m2/min
Si = 0.004
rbi = 0.01
G_vect = np.vectorize(leaky_function)
Si = 0.003
def s_theo(t,Tnew,Snew,rb):
    al = Tnew/Snew
    u = rw**2/(4*al*t)
    s = qo/(4*pi*Tnew)*G_vect(u,rb)
    return s
#
#
rbi = 0.03
po = [Ti,Si,rbi]
ni = 3
tn =t[ni:nt]
sn = sf[ni:nt]
params,resn = curve_fit(s_theo,tn,sn,po)
Tn = params[0]
Sn = params[1]
rbn = params[2]
print ('T = ',Tn,'m2/min')
print('S = ',Sn)
print('r/b = ',rbn)

s1 = s_theo(t,Ti,Si,rbi)
s2 = s_theo(t,Tn,Sn,rbn)
alh = Tn/Sn
u = rw**2/(4*alh*t)
un = rw**2/(4*alh*tn)
x = 1/u
p1 = loglog(x, sf, label='Measured',color = 'black')
p2 = loglog(x, s2, label='Curve fit',color = 'black',  linestyle='none', marker='o')
ll = legend()

sizeOfFont = 14
fontProperties = {'weight' : 'bold', 'size' : sizeOfFont}
a = gca()
gx = xlabel('1/u')
gy = ylabel('')
setp(gx,'fontsize',15,'fontweight','bold')
setp(gy,'fontsize',15,'fontweight','bold')
setp(a,'yscale','log')
setp(a,'xscale','log')
nx = len(x)
show()