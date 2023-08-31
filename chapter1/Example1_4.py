#coding: utf-8
# exemple 1.5
# written by Louis Lamarche
# 15 august 2017
#
import numpy as np
from conversion_mod import *
from CoolProp.CoolProp import *
from CoolProp.HumidAirProp import *
# data
patm = 101.325*1000.0
p1 = 4.0e6 # pascal
p2 = 100.0*1000 # pascal
p3 = 20000.0 # pascal
p4 = 20000.0 # pascal
T1c = 200.0
T1k = T1c + 273.15
print ('T1 = ', T1c,' C')
h1 = PropsSI('H', 'T',T1k,'P', p1, 'water')
print ('h1 = ', h1/1000,' KJ/kg')
s1 = PropsSI('S', 'T',T1k,'P', p1, 'water')
print ('s1 = ', s1/1000,' KJ/kg K')
h2 = h1
print ('h2 = ', h2/1000,' KJ/kg')
h2g = PropsSI('H', 'P', p2, 'Q', 1, 'water')
h2f = PropsSI('H', 'P', p2, 'Q', 0, 'water')
h3 = h2g
x2 = PropsSI('Q', 'P', p2, 'H', h2, 'water')
T2k = PropsSI('T', 'P', p2, 'H', h2, 'water')
T2c = T2k - 273.15
print ('T2 = ', T2c,' C')
s3 = PropsSI('S', 'P', p2, 'Q', 1, 'water')
print ('s3 = ', s3/1000,' kJ/kg K')
h4s = PropsSI('H', 'P', p4, 'S', s3, 'water')
print ('h4s= ', h4s/1000,' KJ/kg')
x4s = PropsSI('Q', 'P', p4, 'S', s3, 'water')
print ('x4s= ', x4s)
etas = 0.85
A = etas*(h3-h4s)/2.0
h5 = PropsSI('H', 'P', p4, 'Q', 0, 'water')
print ('h5= ', h5/1000,' KJ/kg')
s5 = PropsSI('S', 'P', p4, 'Q', 0, 'water')
print ('s5= ', s5/1000,' KJ/kg')
hg = PropsSI('H', 'P', p4, 'Q', 1, 'water')
sg = PropsSI('S', 'P', p4, 'Q', 1, 'water')
print ('sg= ', sg/1000,' KJ/kg')
hfg = hg - h5
h4 = (h3 - A*(1 - h5/hfg))/(1 + A/hfg)
eta3 = (h3 - h4)/(h3- h4s)
x4 = (h4 - h5)/hfg
eta33 = etas*(1 + x4)/2.0
mp =  25.0
W = x2*(h3-h4)*mp
print ('W= ', W/1e6,' MW')
Tok = 25 + 273.15
ho = PropsSI('H', 'T', Tok, 'P', patm, 'water')
so = PropsSI('S', 'T', Tok, 'P', patm, 'water')
Wmax = mp*(h1 - ho - Tok*(s1 - so))
eta1 = W/Wmax

