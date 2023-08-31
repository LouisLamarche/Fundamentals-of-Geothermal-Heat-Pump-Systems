#coding: utf-8
# exemple 1.6
# written by Louis Lamarche
# 21 august 2017
#
import numpy as np
from conversion_mod import *
from CoolProp.CoolProp import *
from CoolProp.HumidAirProp import *
patm = 101.325*1000.0
Toc = 25.0
pa = 1e6
T2p = 160.0
T2pk = T2p + 273.15
T1c = 150.0
T1k = T1c + 273.15
Tok = Toc + 273.15
Cp = 4190.0
T3c = 40.0
T3k = T3c +  273.15
etat = 0.85
etap = 0.80
mp = 12.0     # refrigerant flow rate
mpw  = 80.0   # water flow rate
p1 = PropsSI('P', 'T', T1k, 'Q', 1, 'ipentane')
p4 = p1
p5 = p4
print ('p1 = ', p1/1000,' KPa')
s1 = PropsSI('S', 'P', p1, 'Q', 1, 'ipentane')
print ('s1 = ', s1/1000,' KJ/kg-K')
h1 = PropsSI('H', 'P', p1, 'Q', 1, 'ipentane')
print ('h1 = ', h1/1000,' KJ/kg')
h5 = PropsSI('H', 'P', p1, 'Q', 0, 'ipentane')
print ('h5 = ', h5/1000,' KJ/kg')
s5 = PropsSI('S', 'P', p1, 'Q', 0, 'ipentane')
print ('s5 = ', s5/1000,' KJ/kg K')
p3 = PropsSI('P', 'T', T3k, 'Q', 0, 'ipentane')
p2 = p3
h2s = PropsSI('H', 'P', p2, 'S', s1, 'ipentane')
print ('h2s= ', h2s/1000,' KJ/kg')
T2s = PropsSI('T', 'P', p2, 'S', s1, 'ipentane')
print ('T2s= ', T2s - 273.15,' C')
x2s = PropsSI('Q', 'P', p2, 'S', s1, 'ipentane')
h2 = h1 - etat*(h1 - h2s)
print ('h2 = ', h2/1000,' KJ/kg')
T2k = PropsSI('T', 'P', p2, 'H', h2, 'ipentane')
T2c = T2k - 273.15
print ('T2 = ', T2c,' C')
h3 = PropsSI('H', 'P', p3, 'Q', 0, 'ipentane')
s3 = PropsSI('S', 'P', p3, 'Q', 0, 'ipentane')
rho3 = PropsSI('D', 'P', p3, 'Q', 0, 'ipentane')
h4s = PropsSI('H', 'P', p4, 'S', s3, 'ipentane')
h4s2 = h3 + (p4 - p3)/rho3
print ('h4s = ', h4s/1000,' KJ/kg')
h4 = h3 - (h3 - h4s)/etap
h42 = h3 + (h4s - h3 )/etap
T4k = PropsSI('T', 'P', p4, 'H', h4, 'ipentane')
T4c  = T4k - 273.15
print ('T4 = ', T4c,' C')
print ('h4 = ', h4/1000,' KJ/kg')
Wt = mp*(h1 - h2)
print ('Wturb = ', Wt/1000,' KW')
qc = h2 - h3
wp = h4 - h3
Wp = mp*wp
print ('Wpompe = ', Wp/1000,' KW')
qe = h1 - h4
Qe = mp*qe
Qs = mp*(h5 - h4)
Ql = mp*(h1 - h5)
etath = (Wt - Wp)/(mp*qe)
print ('eta thermique = ', etath)
#
# heat exchanger
#
T5c = T1c
Cmin1 = mpw*Cp
T3p = T2p - Ql/Cmin1
Tp = T3p - T5c
T4p = T3p - Qs/Cmin1
T4pk = T4p + 273.15
Qmax1 = Cmin1*(T2p - T5c)
eff1 = (Ql/Qmax1)
Cr = Qs/(T5c - T4c)
Cmin2 = min(Cmin1,Cr)
Qmax2 = Cmin2*(T3p - T4c)
eff2 = (Qs/Qmax2)

qhm = Cmin1*(T2p - T4c)
hm = PropsSI('H', 'T',T2pk,'P', p4, 'ipentane')
qcm = mp*(hm - h4)
qmax = min(qhm,qcm)

hr = PropsSI('H', 'T',T2pk,'P', pa, 'water')
sr = PropsSI('S', 'T',T2pk, 'Q', 0, 'water')
ho = PropsSI('H', 'T',Tok,'P', patm, 'water')
so = PropsSI('S', 'T',Tok,'P', patm, 'water')
ho = PropsSI('H', 'T',Tok,'Q', 0, 'water')
so = PropsSI('S', 'T',Tok,'Q', 0, 'water')
emax = hr - ho - Tok*(sr - so)
etau = Wt/(mpw*emax)
print ('eta exergy = ', etau)





