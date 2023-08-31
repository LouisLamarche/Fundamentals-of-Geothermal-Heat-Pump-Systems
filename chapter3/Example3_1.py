#coding: utf-8
# exemple 3.1
# written by Louis Lamarche
# 15 august 2017
#
import numpy as np
from CoolProp.CoolProp import *
from CoolProp.HumidAirProp import *
from scipy.optimize import fsolve,minimize,newton,brentq
from tabulate import tabulate
# data
def KC(x):
    return x - 273.15
def CK(x):
    return x + 273.15
fluid = 'R134a'
patm = 101.325*1000.0
p1 = 1.80e5 # pascal
p2 = 1000.0*1000 # pascal
p4 = p1
p3 = p2
T1k = PropsSI('T', 'P', p1, 'Q', 1, fluid)
T1c = KC(T1k)
h1 = PropsSI('H', 'P', p1, 'Q', 1, fluid)
s1 = PropsSI('S', 'P', p1, 'Q', 1, fluid)
s2 = s1
T2k = PropsSI('T', 'P', p2, 'S', s2, fluid)
T2c = KC(T2k)
h2 = PropsSI('H', 'P', p2, 'S', s2, fluid)
T3k = PropsSI('T', 'P', p2, 'Q', 0, fluid)
T3c = KC(T3k)
h3 = PropsSI('H', 'P', p2, 'Q', 0, fluid)
s3 = PropsSI('S', 'P', p2, 'Q', 0, fluid)
h4 = h3
x4 = PropsSI('Q', 'P', p1, 'H', h4, fluid)
s4 = PropsSI('S', 'P', p1, 'Q', x4, fluid)
T4c = T1c
qin = h1 - h4
qout = h2 - h3
win = h2  - h1
COP =qin/win
print (('COP = ' + str(COP)))
Qin = 5000.0
mp = Qin/qin
print (('mdot = ' + str(mp) + ' kg/s'))
head = ['point','T(C)', 'p(kPa)','h(kJ/kg)','s(kJ/kg-K)']
pt1 = ['1','%.1f' % T1c, '%.1f' % (p1/1000),  '%.0f' % (h1/1000), '%.3f' % (s1/1000)]
pt2 = ['2','%.1f' % T2c, '%.1f' % (p2/1000),  '%.0f' % (h2/1000), '%.3f' % (s2/1000)]
pt3 = ['3','%.1f' % T3c, '%.1f' % (p3/1000),  '%.0f' % (h3/1000), '%.3f' % (s3/1000)]
pt4 = ['4','%.1f' % T4c, '%.1f' % (p4/1000),  '%.0f' % (h4/1000), '%.3f' % (s4/1000)]
print(tabulate([pt1,pt2,pt3,pt4], headers= head))

