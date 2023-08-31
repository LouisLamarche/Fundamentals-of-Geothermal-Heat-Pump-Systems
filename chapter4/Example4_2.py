#coding: utf-8
#
# Example 4.2 written by Louis Lamarche 22 sept 2017
#
import numpy as np
from geothermal_md  import *
als = 0.1       # m2/day
alhr = als/24.0 # m2/hr
rb = 0.15/2.0
ks = 2.0
To = 10.0
qb = -40.0
H = 100.0
rbb = rb/H
Ht = H/rb
zob = 0.04
zot = zob*Ht
q1 = 1000.0
q2 = 1800.0
q3 = -200.0
q1p = q1/H
q2p = q2/H
q3p = q3/H
# intermediate time in hours
t1 = 1.0
t2 = 3.0
t3 = 5.0
dt = 0.1
t = np.arange(0,t3+dt,dt)
Fo1 = alhr*t1/rb**2
Fo2 = alhr*t2/rb**2
Fo3 = alhr*t3/rb**2
R1a = G_function(Fo3)
R2a = G_function(Fo3-Fo1)
R3a = G_function(Fo3-Fo2)
R1b = G_function_ils(Fo3)
R2b = G_function_ils(Fo3-Fo1)
R3b = G_function_ils(Fo3-Fo2)
R1c = G_function_fls(Fo3,Ht = Ht,zot = zot)
R2c = G_function_fls(Fo3-Fo1,Ht = Ht,zot = zot)
R3c = G_function_fls(Fo3-Fo2,Ht = Ht,zot = zot)
DT1 = -(q1p*R1a+(q2p-q1p)*R2a+(q3p-q2p)*R3a)/ks
DT2 = -(q1p*R1b+(q2p-q1p)*R2b+(q3p-q2p)*R3b)/ks
DT3 = -(q1p*R1c+(q2p-q1p)*R2c+(q3p-q2p)*R3c)/ks
print ('DTb a) = ',DT1)
print ('DTb b) = ',DT2)
print ('DTb c) = ',DT3)
Tb1 = To + DT1  # Tb with CLS
Tb2 = To + DT2  # Tb with ILS
Tb3 = To + DT3  # Tb with FLS
print ('Tb a) = ',Tb1)
print ('Tb b) = ',Tb2)
print ('Tb c) = ',Tb3)
#
# Eskilson formalism
tc = H**2/(9*alhr)
tb1 = t1/tc
tb2 = t2/tc
tb3 = t3/tc
R2d = g_function_fls(tb3-tb1,rbb = rbb,zob = zob)
R3d = g_function_fls(tb3-tb2,rbb = rbb,zob = zob)
R1d = g_function_fls(tb3,rbb = rbb,zob = zob)
DT4 = (q1p*R1d+(q2p-q1p)*R2d+(q3p-q2p)*R3d)/(2*pi*ks)
Tb4 = To - DT4  # Tb with FLS
print ('Tb c) = ',Tb4)

#
#
#
flag_plot = False
if flag_plot:
    nh = len(t)
    qp = zeros(nh)
    Tb1 = zeros(nh)
    Tb2 = zeros(nh)
    Tb3 = zeros(nh)
    #
    # Calcululation of Tb
    #
    for i in range(0,nh):
        Fo = alhr*t[i]/rb**2
        if t[i] < t1:
            DT1 = q1p*G_function(Fo)/ks             # ICS
            DT2 = q1p*G_function_ils(Fo)/ks         # ILS
            DT3 = q1p*G_function_fls(Fo,Ht = Ht,zot = zot)/ks     # FLS
            qp[i] = q1p
        elif t[i] < t2:
            DT1 = (q1p*G_function(Fo)+(q2p-q1p)*G_function(Fo-Fo1))/ks
            DT2 = (q1p*G_function_ils(Fo)+(q2p-q1p)*G_function_ils(Fo-Fo1))/ks
            DT3 = (q1p*G_function_fls(Fo,Ht = Ht,zot = zot)+(q2p-q1p)*G_function_fls(Fo-Fo1,Ht = Ht,zot = zot))/ks
            qp[i] = q2p
        else:
            DT1 = (q1p*G_function(Fo)+(q2p-q1p)*G_function(Fo-Fo1)+(q3p-q2p)*G_function(Fo-Fo2))/ks
            DT2 = (q1p*G_function_ils(Fo)+(q2p-q1p)*G_function_ils(Fo-Fo1)+(q3p-q2p)*G_function_ils(Fo-Fo2))/ks
            DT3 = (q1p*G_function_fls(Fo,Ht = Ht,zot = zot)+(q2p-q1p)*G_function_fls(Fo-Fo1,Ht = Ht,zot = zot)+(q3p-q2p)*G_function_fls(Fo-Fo2,Ht = Ht,zot = zot))/ks
            qp[i] = q3p
        Tb1[i] = To - DT1  # Tb with ICS
        Tb2[i] = To - DT2  # Tb with ILS
        Tb3[i] = To - DT3  # Tb with FLS
    print ('Tb a) = ',Tb1[nh-1] )
    print ('Tb b) = ',Tb2[nh-1])
    print ('Tb c) = ',Tb3[nh-1] )
    plot(t,Tb1,t,Tb2,t,Tb3,'*')
    legend(('ICS','ILS','FLS'))
    show()
