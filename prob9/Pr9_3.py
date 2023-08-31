#coding: utf-8
#
# example 9.5 (Parameter estimation method)
#
import numpy as np
from geothermal_md import *
from  matplotlib.pyplot import *
from scipy.optimize import curve_fit
import pandas as pd
#
# donnees du probleme
#
gam = np.euler_gamma
M = np.loadtxt('test_trt2.txt')
t = M[:,0]      # time in hours
Tc = M[:,1]     # température Farenheight
Tb = M[:,2]     # température Farenheight
nt = len(Tc)
qp =55     # heat transfer per metre
nt = len(Tc)    # total number of points
To = 11  # Undisturbed temperature
#
#
ks = 2.5        # assumed intial value
rb = 0.065
CC = 2.5e6
#
# first approach
#
# second approach
n1 = 10        # Minimum number of points to take out
nf = 100         # Minimum number of points to keep for the regressoion
nte = nt - nf - n1 + 1
#
ksol = np.zeros(nte)
almin = np.zeros(nte)
Rb = np.zeros(nte)
for i in range(0,nte):
    j1 = n1+i
    x1 = t[j1:nt]
    x = np.log(x1)
    y = Tc[j1:nt]
    y2 = Tb[j1:nt]
    p = np.polyfit(x,y,1)
    p2 = np.polyfit(x,y2,1)
    m = p[0]
    b = p[1]
    b2 = p2[1]
    ksol[i] = qp/(4.0*pi*m)
    almin[i] =   np.exp((b2-To)*4*pi*ksol[i]/qp + gam)*rb**2/4
    Rb[i] = ((b-To)/qp - (np.log(4*almin[i]/rb**2)-gam)/(4*pi*ksol[i]))
nlast = 10
ks = np.mean(ksol[nte-nlast:nte-1])
Rbn = np.mean(Rb[nte-nlast:nte-1])
alsoil = np.mean(almin[nte-nlast:nte-1])
als = alsoil/60
CCs = ks/als
xk = n1+ np.arange(0,nte)
plot(xk,almin)
show()
