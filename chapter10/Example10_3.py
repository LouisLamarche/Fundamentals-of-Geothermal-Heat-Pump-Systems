#coding: utf-8
#
# Exemple 10.3
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
Cp = 4000
power = 4000
L = 200
qp = power/(2*L)
#
# Calcul des Tf
#
#
kp = 0.4
rho = 1000
u = 0.2
ro = 0.02
Rp = 0.04
mp = pi*ro**2*u*rho
CCf = mp*Cp
mu = L/CCf

#
# series fluid enters at z1
z1 = 2
z2 = 3
R1 = 1/(2*pi*ks)*np.log(2*z1/ro) + Rp
R2 = 1/(2*pi*ks)*np.log(2*z2/ro) + Rp
R12 = 1/(2*pi*ks)*np.log((z1+z2)/abs(z2-z1))
zeta1 = 2*(R2 + mu/2 - R12)/(R1+R2-2*R12)
zeta2 = 2*(R1 - mu/2 - R12)/(R1+R2-2*R12)
q1 = qp*zeta1
q2 = qp*zeta2
Rser = 2*(R1*R2 - R12**2 + mu**2/4)/(R1 + R2 - 2*R12)
Rins = 2*((R1+mu/2)*(R2+mu/2) - R12*(R12+mu))/(R1 + R2 - 2*R12)
Tfs = -qp*Rser
Tfis = Tfs - mu*qp
Tfos = Tfs + mu*qp
Tfis3 = -qp*Rins
Tfis2 = -q1*(R12) - q2*(R2 + mu/2)
Tfos2 = -q1*(R12) - q2*(R2 - mu/2)
print ('Tin (series enters higher pipe) = ',Tfis,Tfis3)
print ('Tout (series enters higher pipe) = ',Tfos,Tfos2)
print ('q1 (series enters higher pipe) = ',q1*L)
print ('q2 (series enters higher pipe) = ',q2*L)

#
# ebnters at z = 2
R1 = 1/(2*pi*ks)*np.log(2*z2/ro) + Rp
R2 = 1/(2*pi*ks)*np.log(2*z1/ro) + Rp
R12 = 1/(2*pi*ks)*np.log((z1+z2)/abs(z2-z1))
zeta1 = 2*(R2 + mu/2 - R12)/(R1+R2-2*R12)
zeta2 = 2*(R1 - mu/2 - R12)/(R1+R2-2*R12)
q1 = qp*zeta1
q2 = qp*zeta2
Rser = 2*(R1*R2 - R12**2 + mu**2/4)/(R1 + R2 - 2*R12)
Rins = 2*((R1+mu/2)*(R2+mu/2) - R12*(R12+mu))/(R1 + R2 - 2*R12)
Tfs = -qp*Rser
Tfis = Tfs - mu*qp
Tfos = Tfs + mu*qp
Tfis3 = -qp*Rins
Tfis2 = -q1*(R12) - q2*(R2 + mu/2)
Tfos2 = -q1*(R12) - q2*(R2 - mu/2)
print ('Tin (series enters lower pipe) = ',Tfis,Tfis3)
print ('Tout (series enters lower pipe) = ',Tfos,Tfos2)
print ('q1 (series enters lower pipe) = ',q1*L)
print ('q2 (series enters lower pipe) = ',q2*L)

#
#
# parallell
#
z1 = 2
z2 = 3
R1p = 1/(2*pi*ks)*np.log(2*z1/ro) + Rp
R2p = 1/(2*pi*ks)*np.log(2*z2/ro) + Rp
R12p = 1/(2*pi*ks)*np.log((z1+z2)/abs(z2-z1))
Rpar = 2*(R1p*R2p - R12p**2 + mu/2*R12p + mu/4*(R1p+R2p))/(R1p + R2p + mu - 2*R12p)
Rinp = 2*((R1p+mu/2)*(R2p+mu/2) - R12p**2)/(R1p + R2p + mu - 2*R12p)
zeta1p = 2*(R2p + mu/2 - R12p)/(R1p+R2p+ mu -2*R12p)
zeta2p = 2*(R1p + mu/2 - R12p)/(R1p+R2p+ mu -2*R12p)
q1p = qp*zeta1p
q2p = qp*zeta2p
Tfp = -qp*Rpar
Tfip3 = -qp*Rinp
Tfip = Tfp - mu*qp/2
Tfop = Tfp + mu*qp/2
Tfip2 = -q1p*(R12p) - q2p*(R2p + mu/2)
Tfop2a = -q1p*(R12p) - q2p*(R2p - mu/2)
Tfop2b = -q1p*(R1p-mu/2) - q2p*R12p
Tfop2 = (Tfop2a+Tfop2b)/2
print ('Tin (parallell) = ',Tfip,Tfip2,Tfip3)
print ('Tout (parallell) = ',Tfop,Tfop2)
print ('q1 (parallel) = ',q1p*L)
print ('q2 (parallel) = ',q2p*L)
#
#
# parallell with the total flow rate equal
#
mu2 = mu*2
z1 = 2
z2 = 3
R1p = 1/(2*pi*ks)*np.log(2*z1/ro) + Rp
R2p = 1/(2*pi*ks)*np.log(2*z2/ro) + Rp
R12p = 1/(2*pi*ks)*np.log((z1+z2)/abs(z2-z1))
Rpar = 2*(R1p*R2p - R12p**2 + mu2/2*R12p + mu2/4*(R1p+R2p))/(R1p + R2p + mu2 - 2*R12p)
Rinp = 2*((R1p+mu2/2)*(R2p+mu2/2) - R12p**2)/(R1p + R2p + mu2 - 2*R12p)
zeta1p = 2*(R2p + mu2/2 - R12p)/(R1p+R2p+ mu2 -2*R12p)
zeta2p = 2*(R1p + mu2/2 - R12p)/(R1p+R2p+ mu2 -2*R12p)
q1p = qp*zeta1p
q2p = qp*zeta2p
Tfp = -qp*Rpar
Tfip3 = -qp*Rinp
Tfip = Tfp - mu2*qp/2
Tfop = Tfp + mu2*qp/2
Tfip2 = -q1p*(R12p) - q2p*(R2p + mu2/2)
Tfop2a = -q1p*(R12p) - q2p*(R2p - mu2/2)
Tfop2b = -q1p*(R1p-mu2/2) - q2p*R12p
Tfop2 = (Tfop2a+Tfop2b)/2
print ('Tin (parallell b) = ',Tfip,Tfip2,Tfip3)
print ('Tout (parallell b) = ',Tfop,Tfop2)
print ('q1 (parallel b) = ',q1p*L)
print ('q2 (parallel b) = ',q2p*L)

def R_horiz(t,rp):
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
# series
Ra =  R_horiz(t,ro)
print(Ra)

Tp = To - qp*Ra
Tf = Tp - qp*Rp
Tfi = Tf - power/(2*CCf)
Tfo = Tf + power/(2*CCf)
print('Tfi (Bose parker series) = ',Tfi)
print('Tfo (Bose parker series) = ',Tfo)

