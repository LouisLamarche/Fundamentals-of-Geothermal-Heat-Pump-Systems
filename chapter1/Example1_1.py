#coding: utf-8
# exemple 1.1
# written by Louis Lamarche
# 15 august 2017
#
import numpy as np
from conversion_md import *
from CoolProp.CoolProp import *
from CoolProp.HumidAirProp import *
# data
patm = 101.325*1000.0
p1 = 2.0e6 # pascal
p2 = 700.0*1000 # pascal
p3 = 20000.0 # pascal
p4 = 20000.0 # pascal
T1k = PropsSI('T', 'P', p1, 'Q', 1, 'water')
T1c = T1k - 273.15
print ('T1 = ', T1c,' C')
h1 = PropsSI('H', 'P', p1, 'Q', 1, 'water')
print ('h1 = ', h1/1000,' KJ/kg')
h2 = h1
print( 'h2 = ', h2/1000,' KJ/kg')
T2k = PropsSI('T', 'P', p2, 'H', h2, 'water')
T2c = T2k - 273.15
print ('T2 = ', T2c,' C')
s2 = PropsSI('S', 'P', p2, 'H', h2, 'water')
print ('s2 = ', s2/1000,' kJ/kg K')
h2g = PropsSI('H', 'P', p2, 'Q', 1, 'water')
T3k = PropsSI('T', 'P', p3, 'Q', 1, 'water')
T3c = T3k - 273.15
print ('T3 = ', T3c,' C')
h3s = PropsSI('H', 'P', p3, 'S', s2, 'water')
print ('h3s= ', h3s/1000,' KJ/kg')
x3s = PropsSI('Q', 'P', p3, 'S', s2, 'water')
print ('x3s= ', x3s)
eta = 0.85*(1+x3s)/2.0
A = 0.425*(h1-h3s)
eta = 0.80
h3 = h2 - eta*(h2-h3s)
h4 = PropsSI('H', 'P', p3, 'Q', 0, 'water')
print ('h4= ', h4/1000,' KJ/kg')
s4 = PropsSI('S', 'P', p3, 'Q', 0, 'water')
print ('s4= ', s4/1000,' KJ/kg')
hg = PropsSI('H', 'P', p3, 'Q', 1, 'water')
sg = PropsSI('S', 'P', p3, 'Q', 1, 'water')
print ('sg= ', sg/1000,' KJ/kg')
hfg = hg - h4
x3 = (h3-h4)/hfg
print ('h3= ', h3/1000,' KJ/kg')
print ('x3= ', x3)
h33 = (h2 - A*(1 - h4/hfg))/(1 + A/hfg)
eta2 = (h2 - h33)/(h2 - h3s)
mp = np.sqrt(1.0 - (p2/p1)**2)
mp = 25.0
W = (h2-h3)*mp
print ('W= ', W/1e6,' MW')
