#coding: utf-8
#
# Exemple 10.1
#
from geothermal_md import *
import numpy as np
from  matplotlib.pyplot import *
# sol
ksh = 2
CC = 2.0e6
al = ksh/CC
alhr = al*3600
To = 12.0
q = 0.9
mp = q
ri,rp = sdr_pipe(1.25,11)
kp = 0.4
Cp = 4200
Rp,rcond,rconv = Rp_fct(mp,rp,ri,kp,300,'Water')
#
# temps
#
ta = 10.0*8760.0
tm = ta + 730.0
tf = tm + 4.0
qh =  10000
qm =  3000
qa = 1000
#
# debut du design0
# Calcul des Tf
#
Tfo = 4
Tfi = Tfo - qh/(mp*Cp)
Tf = (Tfi+Tfo)/2.0
#
# Calul des résistances du sol
#
z1 = 0.75
z2 = 1.5
def R_horiz(t):
    X1 = rp/(2.0*np.sqrt(alhr*t))
    X2 = (z2-z1)/(2.0*np.sqrt(alhr*t))
    X3a = z1/(np.sqrt(alhr*t))
    X4 = (z1+z2)/(2.0*np.sqrt(alhr*t))
    Rs1 = ((I_function(X1)+I_function(X2))-(I_function(X3a)+I_function(X4)))/(2*pi*ksh)     # mK/W
    X3b = z2/(np.sqrt(alhr*t))
    Rs2 = ((I_function(X1)+I_function(X2))-(I_function(X3b)+I_function(X4)))/(2*pi*ksh)     # mK/W
    Rs = (Rs1+Rs2)/2.0
    return Rs
Rhf =  R_horiz(tf)
Rh1 =  R_horiz(tf-ta)
Rh2 =  R_horiz(tf-tm)
Rsa = Rhf - Rh1
Rsm = Rh1 - Rh2
Rsh = Rh2

# calcul de la longeur Chauffage
# 1ere iteration
Lchh1 = (qa*Rsa)/(To - Tf )
Lchh2 = (qm*Rsm)/(To - Tf )
Lchh3 = (qh*Rsh)/(To - Tf )
Lchh4 = (qh*Rp)/(To - Tf)
L_hori = Lchh1 + Lchh2 + Lchh3 + Lchh4
print ('L horizontal = ',L_hori,' m')
print ('L trench horizontal = ',L_hori/2,' m')
ks = 2.0
kg = 1.5
rb = 0.08
xc = rb/3
sigma = (kg-ks)/(kg+ks)
Rg,Rga = Rb_linesource(kg,ks,rb,rp,xc)
Rpb = Rg + Rp/2
Fof = alhr*tf/rb**2
Fo1 = alhr*(tf-ta)/rb**2
Fo2 = alhr*(tf-tm)/rb**2
Gf = G_function(Fof)
G1 = G_function(Fo1)
G2 = G_function(Fo2)
Rsav = (Gf-G1)/ks
Rsmv = (G1-G2)/ks
Rshv = (G2)/ks

# calcul de la longeur Chauffage
# 1ere iteration
Lchv1 = (qa*Rsav)/(To - Tf )
Lchv2 = (qm*Rsmv)/(To - Tf )
Lchv3 = (qh*Rshv)/(To - Tf )
Lchv4 = (qh*Rpb)/(To - Tf)
H_vert = Lchv1 + Lchv2 + Lchv3 + Lchv4
print ('H vertical = ',H_vert,' m')
print(Lchv1,Lchv2,Lchv3,Lchv4)
print(Lchh1/2,Lchh2/2,Lchh3/2,Lchh4/2)