#coding: utf-8
# exemple 1.4
# written by Louis Lamarche
# 15 august 2017
#
import numpy as np
from conversion_mod import *
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
s1 = PropsSI('S', 'P', p1, 'Q', 1, 'water')
print ('s1 = ', s1/1000,' KJ/kg K')
h2 = h1
print ('h2 = ', h2/1000,' KJ/kg')
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
etas = 0.84
A = etas*(h2-h3s)/2.0
h4 = PropsSI('H', 'P', p3, 'Q', 0, 'water')
print ('h4= ', h4/1000,' KJ/kg')
s4 = PropsSI('S', 'P', p3, 'Q', 0, 'water')
print ('s4= ', s4/1000,' KJ/kg')
hg = PropsSI('H', 'P', p3, 'Q', 1, 'water')
sg = PropsSI('S', 'P', p3, 'Q', 1, 'water')
print ('sg= ', sg/1000,' KJ/kg')
hfg = hg - h4
h3 = (h2 - A*(1 - h4/hfg))/(1 + A/hfg)
eta2 = (h2 - h3)/(h2 - h3s)
x3 = (h3 - h4)/hfg
eta22 = etas*(1 + x3)/2.0
mp =  25.0
W = (h2-h3)*mp
print ('W= ', W/1e6,' MW')
Tok = 25 + 273.15
ho = PropsSI('H', 'T', Tok, 'P', patm, 'water')
so = PropsSI('S', 'T', Tok, 'P', patm, 'water')
Wmax = mp*(h1 - ho - Tok*(s1 - so))
eta1 = W/Wmax
Wmax2 = mp*(h2 - ho - Tok*(s2 - so))
eta2 = W/Wmax2
print ('eta1= ', eta1)
print ('eta2= ', eta2)
