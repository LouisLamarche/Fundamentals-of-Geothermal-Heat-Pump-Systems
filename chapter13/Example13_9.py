import numpy as np
from scipy.integrate import quad,dblquad
from numba import float64,jit
from math import erfc
from matplotlib.pyplot import *
from time import *

pi = np.pi



def function_Li(Fo,xoj,phi,rpt,h):
#
    zi = -h
    def g_outer(gam):
        xj = xoj + np.cos(phi)*(1+rpt*np.cos(gam))
        yj = np.sin(phi)*(1+rpt*np.cos(gam))
        zj = -h + rpt*np.sin(gam)
        def g_inner(ome):
            xi = np.cos(ome)
            yi = np.sin(ome)
            dij = np.sqrt((xi  - xj)**2 + (yi  - yj)**2 +  (zi  - zj)**2)
            dij_i = np.sqrt((xi  - xj)**2 + (yi  - yj)**2 +  (-zi  - zj)**2)
            y = erfc(dij/(2*np.sqrt(Fo)))/dij -  erfc(dij_i/(2*np.sqrt(Fo)))/dij_i
            return y
        I = quad(g_inner,0,2*pi)
        return I[0]
    I2 = quad(g_outer,0,2*pi)
    return I2[0]/(8*pi**2)

def function_Xi2(Fo,xoj,phi,gam,rpt,h):
#
    zi = -h
    xj = xoj + np.cos(phi)*(1+rpt*np.cos(gam))
    yj = np.sin(phi)*(1+rpt*np.cos(gam))
    zj = -h + rpt*np.sin(gam)
    def g_inner(ome):
        xi = np.cos(ome)
        yi = np.sin(ome)
        dij = np.sqrt((xi  - xj)**2 + (yi  - yj)**2 +  (zi  - zj)**2)
#        dij_i = np.sqrt(dij**2 + 4*h**2)
        dij_i = np.sqrt((xi  - xj)**2 + (yi  - yj)**2 +  (-zi  - zj)**2)
        y = erfc(dij/(2*np.sqrt(Fo)))/dij -  erfc(dij_i/(2*np.sqrt(Fo)))/dij_i
        return y
    I = quad(g_inner,0,2*pi)
    return I[0]/(4*pi)

def function_Xi1(Fo,xoj,phi,rpt,h):
#
    def  f_inner(ome):
        d1 = np.sqrt((xoj + (1 - rpt)*np.cos(phi) - np.cos(ome))**2 + ((1 - rpt)*np.sin(phi) - np.sin(ome))**2)
        d2 = np.sqrt((xoj + (1 + rpt)*np.cos(phi) - np.cos(ome))**2 + ((1 + rpt)*np.sin(phi) - np.sin(ome))**2)
        dij = (d1 + d2)/2.0
        dij_i = np.sqrt(dij**2 + 4*h*h)
        y = erfc(dij/(2*np.sqrt(Fo)))/dij -  erfc(dij_i/(2*np.sqrt(Fo)))/dij_i
        return y
    I = quad(f_inner,0,2*pi)
    return I[0]/(4*pi)
cas = 'b'
if cas == 'a':
    rpt = 0.05
    zt = 3.5
    tit = '$r_p/R = 0.05,h/R = 3.5$'
    fic = 'fig13_Tmean.png'
else:
    rpt = 0.25
    zt = 0.5
    tit = '$r_p/R = 0.25,h/R = 0.5$'
    fic = 'fig13_Tmeanb.png'
phi = 0
Fov = np.logspace(0,1.5,10)
n = len(Fov)
g1 = np.zeros(n)
g2 = np.zeros(n)
g3 = np.zeros(n)
g4 = np.zeros(n)
xoj = 0
for i in range(0,n):
    Fo = Fov[i]
    q = function_Li(Fo,xoj,phi,rpt,zt)
    h3 = function_Xi2(Fo,xoj,phi,pi/2,rpt,zt)
    h4 = function_Xi2(Fo,xoj,phi,pi,rpt,zt)
    h1 = function_Xi2(Fo,xoj,phi,0,rpt,zt)
    h2 = function_Xi2(Fo,xoj,phi,3*pi/2,rpt,zt)
    r = function_Xi1(Fo,xoj,phi,rpt,zt)
    g1[i] = q
    g2[i] = r
    g3[i] = (h1+h4)/2
    g4[i] = (h1+h2+h3+h4)/4
x = np.log(Fov)
p1 = plot(x,g1,color = 'k',linewidth = 2,label = 'Real mean temperature')
p2 = plot(x,g2,color = 'k',linestyle = 'None',marker = 'o',markersize=8,label = 'Xiong''s approximation')
p3 = plot(x,g3,color = 'k',linestyle = 'None',marker = 'x',markersize=8,label = 'Mean at $\gamma = 0,\gamma = \pi$')
p4 = plot(x,g4,color = 'k',linestyle = 'None',marker = '+',markersize=8,label = 'Mean at $\gamma = 0,\gamma = \pi/2,\gamma = \pi,\gamma = 3\pi/2$')
legend(fontsize = 12)
title(tit,fontsize = 14)
ax = gca()
grid(True,which='both')
fts = 16
ftst = 14

xlabel(' Fo',fontname='Times new Roman',fontsize = fts)
ylabel(' $k_s\Delta T/q\'$',fontsize = fts,fontname='Times new Roman')
xticks(fontsize=ftst)
yticks(fontsize=ftst)
tight_layout()
show()




