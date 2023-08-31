#coding: utf-8
#
# Example 7.2 TRT testing analysis
#
import numpy as np
from geothermal_md import *
from  matplotlib.pyplot import *
#
gam = np.euler_gamma
M = np.loadtxt("test_trt.txt")
ti = M[:,0]      # time in sec
Ti = M[:,1]     # température Farenheight
qp = 30      # heat transfer per metre
nt = len(Ti)    # total number of points
To = 8       # Undisturbed temperature
rb = 0.07
CC = 2.5e6
als = 2.26/CC    # Assumed thermal diffusivity
alhri = als*3600
rb = 0.126/2.0
t1 = 24
t2 = t1 + 8
n_inj = np.sum(ti < t2)     # nombre de points a enlever
Tc = Ti[n_inj:nt]
t = ti[n_inj:nt]
nt2 = len(Tc)
#
#
#
# second approach
n1 = 20        # Minimum number of points to take out
nf = 10         # Minimum number of points to keep for the regressoion
nte = nt2 - nf - n1 + 1
#
ksol = np.zeros(nte)
Rb = np.zeros(nte)
for i in range(0,nte):
    j1 = n1+i
    x1 = t[j1:nt2]
    x2 = t[j1:nt2] - t1
    x3 = t[j1:nt2] - t2
    x = np.log(x1/x2*x3**2)
    y = Tc[j1:nt2]
    p = np.polyfit(x,y,1)
    m = p[0]
    b = p[1]
    ksol[i] = qp/(4.0*pi*m)
    alhr =   ksol[i]/CC*3600.0 # m2/hr
    Rb[i] = ((b-To)/qp - (np.log(4*alhr/rb**2)-2*gam)/(4*pi*ksol[i]))/2
# The mean from the last
nlast = 10
ks = np.mean(ksol[nte-nlast:nte-1])
Rbn = np.mean(Rb[nte-nlast:nte-1])
print (ks,Rbn)
xk = n1+ np.arange(0,nte)
plot(xk,ksol)
show()
