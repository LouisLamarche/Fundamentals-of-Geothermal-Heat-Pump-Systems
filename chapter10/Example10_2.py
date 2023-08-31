#coding: utf-8
#
# Exemple 10.2
#
from geothermal_md import *
import  numpy as np
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
print(Tfo,Tfi)
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

tshift = 11     # Valeur pour Montreal ( en jour)
amp = 10      # amplitude de variation de la temperature dans l'annéee
tch = 15         #  jour le plus froid
z  = (z1+z2)/2
als = alhr*24
dTch,dTmax = Delta_T_earth(z,tch,tshift,als,amp)  # diffusivite en m2/jr
Ton = To -dTmax
print(Ton)
Tom = To  + dTch
print(Tom)
# calcul de la longeur Chauffage
# 1ere iteration
Lch1 = (qa*Rsa)/(Ton - Tf)
Lch2 = (qm*Rsm)/(Ton - Tf)
Lch3 = (qh*Rsh)/(Ton - Tf)
Lch4 = (qh*Rp)/(Ton - Tf)
L_hori = Lch1 + Lch2 + Lch3 + Lch4
print ('L horizontal = ',L_hori)
Lch1 = (qa*Rsa)/(Tom - Tf)
Lch2 = (qm*Rsm)/(Tom - Tf)
Lch3 = (qh*Rsh)/(Tom - Tf)
Lch4 = (qh*Rp)/(Tom - Tf)
L_hori = Lch1 + Lch2 + Lch3 + Lch4
print ('L horizontal = ',L_hori)
