#coding: utf-8
#
# Exemple 10.4
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
mp1 = 0.25  # parallell
mp2 = mp1*2 # series
CCf1 = mp1*Cp
CCf2 = mp2*Cp
mu1 = L/CCf1
mu2 = L/CCf2
di1,do1 = sdr_pipe(1.25,11)
ro1 = do1/2
ri1 = di1/2
Ac1 = pi*ri1**2
u1 = mp1/(Ac1*1000)
Rp1,rcond,rconv = Rp_fct(mp1,ro1,ri1,kp,300,'Water')
di2,do2 = sdr_pipe(2.0,11)
ro2 = do2/2
ri2 = di2/2
Ac2 = pi*ri2**2
u2 = mp2/(Ac2*1000)
Rp2,rcond,rconv = Rp_fct(mp2,ro2,ri2,kp,300,'Water')
#
# series
z1 = 2
z2 = 3
R1 = 1/(2*pi*ks)*np.log(2*z1/ro2) + Rp2
R2 = 1/(2*pi*ks)*np.log(2*z2/ro2) + Rp2
R12 = 1/(2*pi*ks)*np.log((z1+z2)/abs(z2-z1))
zeta1 = 2*(R2 + mu2/2 - R12)/(R1+R2-2*R12)
zeta2 = 2*(R1 - mu2/2 - R12)/(R1+R2-2*R12)
q1 = qp*zeta1
q2 = qp*zeta2
Rser = 2*(R1*R2 - R12**2 + mu2**2/4)/(R1 + R2 - 2*R12)
Rins = 2*((R1+mu2/2)*(R2+mu2/2) - R12*(R12+mu2))/(R1 + R2 - 2*R12)
Tfs = -qp*Rser
Tfis = Tfs - mu2*qp
Tfos = Tfs + mu2*qp
Tfis3 = -qp*Rins
Tfis2 = -q1*(R12) - q2*(R2 + mu2/2)
Tfos2 = -q1*(R12) - q2*(R2 - mu2/2)
print ('Tin (series enters higher pipe) = ',Tfis,Tfis3)
print ('Tout (series enters higher pipe) = ',Tfos,Tfos2)
print ('q1 (series) = ',q1*L)
print ('q2 (series) = ',q2*L)
#
#
# parallell
#
z1 = 2
z2 = 3
R1p = 1/(2*pi*ks)*np.log(2*z1/ro1) + Rp1
R2p = 1/(2*pi*ks)*np.log(2*z2/ro1) + Rp1
R12p = 1/(2*pi*ks)*np.log((z1+z2)/abs(z2-z1))
Rpar = 2*(R1p*R2p - R12p**2 + mu1/2*R12p + mu1/4*(R1p+R2p))/(R1p + R2p + mu1 - 2*R12p)
Rinp = 2*((R1p+mu1/2)*(R2p+mu1/2) - R12p**2)/(R1p + R2p + mu1 - 2*R12p)
zeta1p = 2*(R2p + mu1/2 - R12p)/(R1p+R2p+ mu1 -2*R12p)
zeta2p = 2*(R1p + mu1/2 - R12p)/(R1p+R2p+ mu1 -2*R12p)
q1p = qp*zeta1p
q2p = qp*zeta2p
Tfp = -qp*Rpar
Tfip3 = -qp*Rinp
Tfip = Tfp - mu1*qp/2
Tfop = Tfp + mu1*qp/2
Tfip2 = -q1p*(R12p) - q2p*(R2p + mu1/2)
Tfop2a = -q1p*(R12p) - q2p*(R2p - mu1/2)
Tfop2b = -q1p*(R1p-mu1/2) - q2p*R12p
Tfop2 = (Tfop2a+Tfop2b)/2
print ('Tin (parallell) = ',Tfip,Tfip2,Tfip3)
print ('Tout (parallell) = ',Tfop,Tfop2)
print ('q1 (parallel) = ',q1p*L)
print ('q2 (parallel) = ',q2p*L)

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
Ra =  R_horiz(t,ro2)
print(Ra)

Tp = To - qp*Ra
Tf = Tp - qp*Rp2
Tfi = Tf - power/(2*CCf2)
Tfo = Tf + power/(2*CCf2)
print('Tfi (Bose parker series) = ',Tfi)
print('Tfo (Bose parker series) = ',Tfo)

# paralell
Rb =  R_horiz(t,ro1)
print(Rb)
Tp = To - qp*Rb
Tf = Tp - qp*Rp1
Tfi = Tf - power/(4*CCf1)
Tfo = Tf + power/(4*CCf1)
print('Tfi (Bose parker para) = ',Tfi)
print('Tfo (Bose parker para) = ',Tfo)
