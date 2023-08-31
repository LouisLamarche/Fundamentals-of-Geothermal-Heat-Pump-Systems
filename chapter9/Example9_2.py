#coding: utf-8
#
# Example 9.2 TRT testing analysis
#
import numpy as np
from geothermal_md import *
from  matplotlib.pyplot import *
#
gam = 0.5772157
M = np.loadtxt("..\\data\\TRT_test1.txt")
t = M[:,0]      # time in hours
Tf = M[:,1]     # température Farenheight
Tc = (Tf-32.0)/1.8# température Celsius
qp = 55.4      # heat transfer per metre
nt = len(Tc)    # total number of points

To = 10.4       # Undisturbed temperature
CC = 2.5e6     # Volumetric heat capacity
als = 2.26/CC    # Assumed thermal diffusivity
alhri = als*3600
rb = 0.075      # borehole radius (m)
t_crit = 6*rb**2/alhri   #  tmin
n_crit = np.sum(t<t_crit) #  number of point to take out
#
# first approach
#
x = np.log(t[n_crit:nt])
y = Tc[n_crit:nt]
p = np.polyfit(x,y,1)
m1 = p[0]
b1 = p[1]
ks1 = qp/(4.0*pi*m1)
alhr =   ks1/CC*3600.0 # m2/hr
Rb1 = (b1-To)/qp - (np.log(4*alhr/rb**2)-gam)/(4*pi*ks1)
y2 = np.polyval(p,x)
xx = np.log(t)
print (ks1,Rb1)


#
#
# second approach
n1 = 10        # Minimum number of points to take out
nf = 50         # Minimum number of points to keep for the regressoion
nte = nt - nf - n1 + 1
#
ksol = np.zeros(nte)
Rb = np.zeros(nte)
for i in range(0,nte):
    j1 = n1+i
    x = np.log(t[j1:nt])
    y = Tc[j1:nt]
    p = np.polyfit(x,y,1)
    m = p[0]
    b = p[1]
    ksol[i] = qp/(4.0*pi*m)
    alhr =   ksol[i]/CC*3600.0 # m2/hr
    Rb[i] = (b-To)/qp - (np.log(4*alhr/rb**2)-gam)/(4*pi*ksol[i])
# The mean from the last
nlast = 20
ks = np.mean(ksol[nte-nlast:nte-1])
Rbn = np.mean(Rb[nte-nlast:nte-1])
print (ks,Rbn)
xk = n1+ np.arange(0,nte)
plot(xk,ksol)
ss = 0
show()
