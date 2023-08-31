#coding: utf-8
#
# Exemple 10.5
#
from geothermal_md import *
import numpy as np
from  matplotlib.pyplot import *
# soil
ks = 1.5
CC = 2.0e6
al = ks/CC
alhr = al*3600
To = 0.0
rp = 0.02
Ac = pi*rp**2
uf = 0.2
Q = uf*Ac
mp = Q*1000
Cp =  4000
CCf = mp*Cp
power = 4000
L = 200
mu = L/CCf
qp = power/(2*L)
Rp = 0.04
hconv = 1/(2*pi*rp*Rp)
zh = 2
zb = 3
#
# Calcul des Tf
#
#
z1 = zh
z2 = zb
R1 = 1/(2*pi*ks)*np.log(2*z1/rp) + Rp
R2 = 1/(2*pi*ks)*np.log(2*z2/rp) + Rp
R12 = 1/(2*pi*ks)*np.log((z1+z2)/abs(z2-z1))
zeta1sn = 2*(R2 + mu/2 - R12)/(R1+R2-2*R12)
zeta2sn = 2*(R1 - mu/2 - R12)/(R1+R2-2*R12)
q1s = qp*zeta1sn
q2s = qp*zeta2sn
Rser = 2*(R1*R2 - R12**2 + mu**2/4)/(R1 + R2 - 2*R12)
Rins = 2*((R1+mu/2)*(R2+mu/2) - R12*(R12+mu))/(R1 + R2 - 2*R12)
Tfs = -qp*Rser
Tfi1a = Tfs - mu*qp
Tfo2a = Tfs + mu*qp
Tfi1b = -qp*Rins
Tfi1c = -q1s*(R1 + mu/2) - q2s*(R12)
Tfi2a = -q1s*(R12) - q2s*(R2 + mu/2)
Tfo2b = -q1s*(R12) - q2s*(R2 - mu/2)
print ('Tin (series enters higher pipe) = ',Tfi1a,Tfi1b,Tfi1c)
print ('Tout (series enters higher pipe) = ',Tfo2a,Tfo2b)
print ('q1 (series) = ',q1s*L)
print ('q2 (series) = ',q2s*L)
#

x = np.arange(0,2*L+1)
xm1 = np.arange(0,L+1)
xm2 = np.arange(1,L+1)
nx = len(x)
nm = int((nx - 1)/2)
dat2 = np.loadtxt('..\\data\\Exemple10_5.txt')
x2 = dat2[:,0]
Tci = dat2[:,1]
nc = len(Tci)
Tc = np.zeros(nc)

R1d = (R1*R2 - R12**2)/(R2 - R12)
R2d = (R1*R2 - R12**2)/(R1 - R12)
R12d = (R1*R2 - R12**2)/(R12)


k1 = mu/R1d
k2 = mu/R2d
k12 = mu/R12d
km = (k1+k2)/2
eta = np.sqrt(km**2 + 2*k12*km)
Req = 2*(R1d*R2d)/(R1d+R2d)*eta/np.tanh(eta)
den = (eta*np.cosh(eta) + km*np.sinh(eta))
tho = (eta*np.cosh(eta) - km*np.sinh(eta))/den
Req2 = mu*(tho+1)/(1-tho)
Tf = -qp*Req2
Tin = -2*qp*mu/(1-tho)
Tcl = np.zeros(nx)
Tla = np.zeros(nx)
xt = xm1/L
xd = xt[0:nm-1]
Tla[0:nm+1] = Tin*np.exp(xt*(k1-k2)/2)*(eta*np.cosh(eta*(1-xt)) + km*np.sinh(eta*(1-xt)))/den
Td = Tin*np.exp(xt*(k1-k2)/2)*(eta*np.cosh(eta*(1-xt)) - km*np.sinh(eta*(1-xt)))/den
Tcl[0:nm+1] = Tfi1a + (Tfi2a - Tfi1a)*xm1/L
Tcl[nm+1:nx] = Tfi2a + (Tfo2a - Tfi2a)*xm2/L
for i in range(0,nc):
    Tc[i] = Tci[nc-i-1]
for i in range(0,nm+1):
    Tla[nx-i-1] = Td[i]
Tcin = Tc[0]
Tcout = Tc[nc-1]
theoc = (Tcout - To)/(Tcin - To)
Rbce = mu*(1+theoc)/(1 - theoc)
def R_horiz(t):
    X1 = rp/(2.0*np.sqrt(alhr*t))
    X2 = (z2-z1)/(2.0*np.sqrt(alhr*t))
    X3a = z1/(np.sqrt(alhr*t))
    X4 = (z1+z2)/(2.0*np.sqrt(alhr*t))
    Rs1 = ((I_function(X1)+I_function(X2))-(I_function(X3a)+I_function(X4)))/(2*pi*ks)     # mK/W
    X3b = z2/(np.sqrt(alhr*t))
    Rs2 = ((I_function(X1)+I_function(X2))-(I_function(X3b)+I_function(X4)))/(2*pi*ks)     # mK/W
    Rs = (Rs1+Rs2)/2.0
    return Rs
t =  1e7
R =  R_horiz(t)
R = R + Rp
print('Rb image = ',R)
print('Rb comsol = ',Rbce)
print('Rb nouv = ',Req)
print('Rb ans = ',Rser)
Tf_ima = To - qp*R
Tfi_ima  = Tf_ima - qp*mu
Tfo_ima  = Tf_ima + qp*mu
T_ima = Tfi_ima + (Tfo_ima - Tfi_ima)*x/(2*L)
print('Tout comsol = ',Tcout)
print('Tout claesson  = ',Tfo2b)
Tout_n = Tin*tho
print('Tout nouv  = ',Tout_n)
print('Tout imag  = ',Tfo_ima)
print('Tin comsol = ',Tcin)
print('Tin claesson  = ',Tfi1a)
Tout_n = Tin*tho
print('Tin nouv  = ',Tin)
print('Tin imag  = ',Tfi_ima)
ms = 8
p1 = plot(x2[0:len(x2)-1:400],Tc[0:len(x2)-1:400],color = 'k',linestyle = '-',marker = 'x',markersize=ms,label = 'Comsol ')
p2 = plot(x[0:len(x)-1:20],Tcl[0:len(x)-1:20],color = 'k',linestyle = '-',marker = '+',markersize=ms,label = 'Claesson')
p3 = plot(x[0:len(x)-1:20],Tla[0:len(x)-1:20],color = 'k',linestyle = ':',marker = 'o',markersize=ms,label = 'Lamarche')
p4 = plot(x[0:len(x)-1:20],T_ima[0:len(x)-1:20],color = 'k',linestyle = ':',marker = '^',markersize=ms,label = 'Image')
legend(fontsize = 11)
ax = gca()
grid(True,which='both')
fts = 18
ftst = 14
ylabel(' T',fontname='Times new Roman',fontsize = fts)
xlabel('x (m)',fontsize = fts,fontname='Times new Roman')
xticks(fontsize=ftst)
yticks(fontsize=ftst)
tight_layout()
show()



